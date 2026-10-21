#!/usr/bin/env python3
"""Run selected development work through the workspace's pinned Nix producers."""

from __future__ import annotations

import argparse
import fcntl
import grp
import json
import math
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import uuid
import xml.etree.ElementTree as ET
from contextlib import ExitStack, contextmanager
from pathlib import Path, PurePosixPath

REQUEST_BUILD_EXPRESSION = """let
  flake = builtins.getFlake (builtins.getEnv "GAMMALOOP_DEV_SOURCE");
  system = builtins.getEnv "GAMMALOOP_DEV_SYSTEM";
in flake.lib.${system}.development.run (builtins.fromJSON
  (builtins.readFile (builtins.getEnv "GAMMALOOP_DEV_REQUEST_FILE")))"""
PREPARE_CI_APPLY = """development: development.prepareCi (builtins.fromJSON
  (builtins.readFile (builtins.getEnv "GAMMALOOP_DEV_REQUEST_FILE")))"""
COMMANDS = ("test", "build", "check", "clippy", "prepare-ci")
PROFILES = ("ci_gammaloop", "ci_full", "ci_slow")
# These are execution options only. Cargo selection, config, partitioning and
# profile flags must have a structured owner rather than bypassing the request.
RUNTIME_OPTIONS = {
    "--run-ignored": 1,
    "--test-threads": 1,
    "-j": 1,
    "--retries": 1,
    "--status-level": 1,
    "--final-status-level": 1,
    "--failure-output": 1,
    "--success-output": 1,
    "--color": 1,
    "--fail-fast": 0,
    "--ff": 0,
    "--no-fail-fast": 0,
    "--nff": 0,
    "--no-capture": 0,
    "--nocapture": 0,
}
LIBTEST_OPTIONS = {
    "--ignored": 0,
    "--include-ignored": 0,
    "--exact": 0,
    "--skip": 1,
    "--nocapture": 0,
    "--test-threads": 1,
}
TRACE_ENVIRONMENT = ("GL_DISPLAY_FILTER", "GL_LOGFILE_FILTER", "GL_ALL_LOG_FILTER")
INTEGER_ENVIRONMENT = ("RUST_MIN_STACK", "RAYON_NUM_THREADS")


class DevError(Exception):
    pass


class RuntimeArgument(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        arguments = list(namespace.runtime or [])
        arguments.append(option_string)
        if self.nargs != 0:
            arguments.append(values)
        namespace.runtime = arguments


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("command", choices=COMMANDS)
    if not argv or argv[0] in ("-h", "--help"):
        parser.parse_args(argv)
    command = argv[0]
    if command not in COMMANDS:
        parser.error(f"unknown command: {command}")
    parser = argparse.ArgumentParser(
        prog=f"scripts/dev.py {command}", allow_abbrev=False
    )
    parser.set_defaults(command=command)
    if command != "prepare-ci":
        parser.add_argument("-p", "--package", action="append", default=[])
    parser.add_argument(
        "--dry-run", action="store_true", help="inspect without building"
    )
    parser.add_argument(
        "--private-cache",
        type=Path,
        help="use private storage instead of a provisioned shared registry",
    )
    parser.add_argument(
        "--keep-generations", type=int, help="private history per context (default: 4)"
    )
    parser.add_argument(
        "--max-generations", type=int, help="private total history cap (default: 64)"
    )
    parser.add_argument(
        "--max-age-days", type=float, help="private history age cap (default: 30)"
    )
    parser.add_argument(
        "--system", help="Nix system (defaults to builtins.currentSystem)"
    )
    if command == "build":
        parser.add_argument(
            "--python", action="store_true", help="build the Python module"
        )
    if command == "prepare-ci":
        parser.add_argument("--context", action="append", default=[])
        parser.add_argument(
            "--manifest", type=Path, default=Path("nix/development-seeds.json")
        )
        parser.add_argument(
            "--destination",
            default="https://cache.nix-ci.com?compression=xz&parallel-compression=true",
        )
    if command == "test":
        parser.add_argument("--group")
        parser.add_argument("--profile", choices=PROFILES, default=PROFILES[0])
        parser.add_argument("-E", "--filter-expr", action="append", default=[])
        parser.add_argument("--test", action="append", default=[])
        parser.add_argument("--lib", action="store_true")
        parser.add_argument(
            "--no-tests", choices=("fail", "warn", "pass"), default="fail"
        )
        parser.add_argument("--collect-logs", action="store_true")
        parser.add_argument("--memory-limit-gb", type=float)
        parser.add_argument("--rerun", "--rebuild", action="store_true")
        parser.add_argument("names", nargs="*")
        for option, arity in RUNTIME_OPTIONS.items():
            parser.add_argument(
                option,
                action=RuntimeArgument,
                dest="runtime",
                nargs=None if arity else 0,
            )
    before, separator, after = list(argv[1:]), [], []
    if "--" in before:
        position = before.index("--")
        before, after = before[:position], before[position + 1 :]
        separator = ["--"]
    args = parser.parse_intermixed_args(before)
    if args.keep_generations is not None and args.keep_generations < 1:
        parser.error("--keep-generations must be positive")
    if (args.max_generations is not None and args.max_generations < 1) or (
        args.max_age_days is not None
        and (not math.isfinite(args.max_age_days) or args.max_age_days <= 0)
    ):
        parser.error("--max-generations and --max-age-days must be positive and finite")
    if separator and command != "test":
        parser.error("runtime arguments after -- are supported only for test")
    if command == "test":
        validate_runtime(separator + after)
        args.nextest_args = (args.runtime or []) + args.names + separator + after
        for binary in args.test:
            if not re.fullmatch(r"[A-Za-z0-9_-]+", binary):
                parser.error(f"unsupported test binary name: {binary!r}")
        if any(not expression.strip() for expression in args.filter_expr):
            parser.error("filter expressions must not be empty")
        if args.memory_limit_gb is not None and (
            not math.isfinite(args.memory_limit_gb) or args.memory_limit_gb <= 0
        ):
            parser.error("--memory-limit-gb must be positive and finite")
    return args


def validate_runtime(arguments):
    """Validate tokens without splitting, evaluating or rewriting caller strings."""
    libtest = False
    index = 0
    while index < len(arguments):
        token = arguments[index]
        index += 1
        if token == "--" and not libtest:
            libtest = True
            continue
        if not token.startswith("-") or token == "-":
            continue
        option, equals, value = token.partition("=")
        options = LIBTEST_OPTIONS if libtest else RUNTIME_OPTIONS
        if option not in options:
            raise DevError(
                f"unsupported runtime option {option!r}; use structured package/group/profile/filter selection"
            )
        arity = options[option]
        if equals and (not arity or not value):
            raise DevError(f"invalid runtime option: {token!r}")
        if arity and not equals:
            if index == len(arguments) or arguments[index].startswith("-"):
                raise DevError(f"missing value for runtime option {option}")
            index += 1


def requests_for(args, catalog, environment=None):
    """Route selection through existing contexts; names never narrow packages."""
    packages = catalog["packages"]
    requested = list(dict.fromkeys(args.package))
    if unknown := set(requested) - set(packages):
        raise DevError(f"unknown Cargo packages: {', '.join(sorted(unknown))}")
    if args.command != "test":
        if args.command == "build" and args.python:
            if requested and requested != ["gammaloop-api"]:
                raise DevError("build --python requires only package gammaloop-api")
            return [
                {"command": "build", "packages": ["gammaloop-api"], "variant": "python"}
            ]
        return [{"command": args.command, "packages": requested or packages}]
    environment = os.environ if environment is None else environment
    runtime_env = {}
    for name in TRACE_ENVIRONMENT + INTEGER_ENVIRONMENT + ("INSTA_UPDATE",):
        if name not in environment:
            continue
        value = environment[name]
        if not value or "\0" in value:
            raise DevError(f"invalid runtime environment value for {name}")
        if name in INTEGER_ENVIRONMENT and (
            not re.fullmatch(r"[0-9]+", value) or int(value) <= 0
        ):
            raise DevError(f"{name} must be a positive integer")
        if name == "INSTA_UPDATE" and value != "no":
            raise DevError("development tests require INSTA_UPDATE=no")
        runtime_env[name] = value
    nextest_args = list(args.nextest_args)
    before_separator = (
        nextest_args[: nextest_args.index("--")]
        if "--" in nextest_args
        else nextest_args
    )
    threads_explicit = any(
        arg in ("--test-threads", "-j") or arg.startswith(("--test-threads=", "-j"))
        for arg in before_separator
    )
    if "--" in nextest_args:
        threads_explicit |= any(
            arg == "--test-threads" or arg.startswith("--test-threads=")
            for arg in nextest_args[nextest_args.index("--") + 1 :]
        )
    if not threads_explicit and "NEXTEST_TEST_THREADS" in environment:
        threads = environment["NEXTEST_TEST_THREADS"]
        if not re.fullmatch(r"[0-9]+", threads) or int(threads) <= 0:
            raise DevError("NEXTEST_TEST_THREADS must be a positive integer")
        nextest_args = ["--test-threads", threads] + nextest_args
    rerun = uuid.uuid4().hex if args.rerun else None
    groups = catalog["groups"]
    if args.group:
        groups = [group for group in groups if group["name"] == args.group]
        if not groups:
            raise DevError(f"unknown test group: {args.group}")
    else:
        groups = [group for group in groups if group["default"]]
    coverage = {package for group in groups for package in group["packages"]}
    if unavailable := set(requested) - coverage:
        raise DevError(
            f"packages outside selected test groups: {', '.join(sorted(unavailable))}; choose --group explicitly"
        )
    filters = [f"({expression})" for expression in args.filter_expr]
    expression = " or ".join(filters)
    binaries = [f"(kind(test) and binary(={name}))" for name in args.test]
    if args.lib:
        binaries.append("kind(lib)")
    if binaries:
        binary_filter = "(" + " or ".join(binaries) + ")"
        expression = (
            f"({expression}) and {binary_filter}" if expression else binary_filter
        )
    requests = []
    for group in groups:
        selected = [
            package
            for package in group["packages"]
            if not requested or package in requested
        ]
        if selected:
            requests.append(
                {
                    "command": "test",
                    "packages": selected,
                    "group": group["name"],
                    "profile": args.profile,
                    "filter": expression or None,
                    "nextestArgs": nextest_args,
                    "runtimeEnv": runtime_env,
                    "collectLogs": args.collect_logs
                    or bool(environment.get("GL_TEST_LOG_DIR")),
                    "memoryLimitGb": args.memory_limit_gb,
                    # Empty groups are valid; the frontend enforces the overall policy.
                    "noTests": "pass",
                }
            )
            if rerun:
                requests[-1]["rerun"] = rerun
    if not requests:
        raise DevError("selection contains no packages")
    return requests


class SourceSnapshot:
    """Copy a VCS allowlist from current disk contents, never an old commit."""

    def __init__(self, root):
        self.root = Path(root).resolve()

    def files(self):
        if (self.root / ".jj").is_dir():
            tracked = subprocess.check_output(
                [
                    "jj",
                    "--ignore-working-copy",
                    "--at-op=@",
                    "file",
                    "list",
                    "-T",
                    'path ++ "\\0"',
                ],
                cwd=self.root,
            )
            git_root = (
                subprocess.check_output(
                    [
                        "jj",
                        "--ignore-working-copy",
                        "--at-op=@",
                        "git",
                        "root",
                    ],
                    cwd=self.root,
                )
                .decode()
                .strip()
            )
            git = [
                "git",
                "--no-optional-locks",
                f"--git-dir={git_root}",
                f"--work-tree={self.root}",
            ]
            # The backing index belongs to the colocated workspace, often main.
            # A nonexistent index lets Git apply ignores to every current file.
            with tempfile.TemporaryDirectory(
                prefix="gammaloop-dev-index-"
            ) as temporary:
                others = subprocess.check_output(
                    git
                    + [
                        "ls-files",
                        "--others",
                        "--exclude-standard",
                        "--exclude=.git",
                        "--exclude=.jj",
                        "-z",
                    ],
                    cwd=self.root,
                    env=os.environ | {"GIT_INDEX_FILE": str(Path(temporary) / "index")},
                )
            paths = tracked + others
        else:
            paths = subprocess.check_output(
                [
                    "git",
                    "--no-optional-locks",
                    "ls-files",
                    "--cached",
                    "--others",
                    "--exclude-standard",
                    "--exclude=.git",
                    "--exclude=.jj",
                    "-z",
                ],
                cwd=self.root,
            )
        result = set()
        for raw in paths.split(b"\0"):
            if not raw:
                continue
            name = os.fsdecode(raw)
            path = PurePosixPath(name)
            if path.is_absolute() or ".." in path.parts or not path.parts:
                raise DevError(f"unsafe VCS path: {name!r}")
            if any(part in (".git", ".jj") for part in path.parts):
                continue
            result.add(name)
        return sorted(result)

    @contextmanager
    def parent_fd(self, name):
        with ExitStack() as stack:
            descriptor = os.open(self.root, os.O_RDONLY | os.O_DIRECTORY)
            stack.callback(os.close, descriptor)
            for part in PurePosixPath(name).parts[:-1]:
                descriptor = os.open(
                    part,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                    dir_fd=descriptor,
                )
                stack.callback(os.close, descriptor)
            yield descriptor

    @staticmethod
    def signature(info):
        return (
            info.st_dev,
            info.st_ino,
            info.st_mode,
            info.st_size,
            info.st_mtime_ns,
            info.st_ctime_ns,
        )

    def copy(self, destination):
        names = self.files()
        signatures = {}
        for name in names:
            target = destination / name
            try:
                with self.parent_fd(name) as parent:
                    leaf = PurePosixPath(name).name
                    info = os.stat(leaf, dir_fd=parent, follow_symlinks=False)
                    signatures[name] = self.signature(info)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    if stat.S_ISLNK(info.st_mode):
                        target.symlink_to(os.readlink(leaf, dir_fd=parent))
                    elif stat.S_ISREG(info.st_mode):
                        descriptor = os.open(
                            leaf,
                            os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK,
                            dir_fd=parent,
                        )
                        with (
                            os.fdopen(descriptor, "rb") as source,
                            target.open("xb") as output,
                        ):
                            if (
                                self.signature(os.fstat(source.fileno()))
                                != signatures[name]
                            ):
                                raise DevError(
                                    f"source changed during snapshot: {name}"
                                )
                            shutil.copyfileobj(source, output)
                            if (
                                self.signature(os.fstat(source.fileno()))
                                != signatures[name]
                            ):
                                raise DevError(
                                    f"source changed during snapshot: {name}"
                                )
                        target.chmod(0o755 if info.st_mode & 0o111 else 0o644)
                    elif stat.S_ISDIR(info.st_mode):
                        # A tracked file may have become a directory. Its old
                        # file is deleted; new nonignored children are allowlisted
                        # independently, without recursively copying this path.
                        continue
                    else:
                        raise DevError(
                            f"unsupported source file (directory, submodule or special file): {name}"
                        )
            except (FileNotFoundError, NotADirectoryError):
                # A deletion already present when capture starts belongs in the snapshot.
                if name in signatures:
                    raise DevError(f"source changed during snapshot: {name}")
                signatures[name] = None
        if names != self.files():
            raise DevError("source file list changed during snapshot; retry")
        for name, signature in signatures.items():
            try:
                with self.parent_fd(name) as parent:
                    current = self.signature(
                        os.stat(
                            PurePosixPath(name).name,
                            dir_fd=parent,
                            follow_symlinks=False,
                        )
                    )
            except (FileNotFoundError, NotADirectoryError):
                current = None
            if signature != current:
                raise DevError(f"source changed during snapshot: {name}; retry")

    def store(self):
        with tempfile.TemporaryDirectory(prefix="gammaloop-dev-source-") as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            self.copy(source)
            # The fixed basename makes equivalent worktrees the same NAR input.
            return subprocess.check_output(
                [
                    "nix",
                    "store",
                    "add-path",
                    "--name",
                    "source",
                    str(source),
                ],
                cwd=self.root,
                text=True,
            ).strip()


class NixDevelopment:
    def __init__(self, root, source, system=None):
        self.root = root
        self.source = source
        self.system = (
            system
            or subprocess.check_output(
                [
                    "nix",
                    "eval",
                    "--impure",
                    "--raw",
                    "--expr",
                    "builtins.currentSystem",
                ],
                cwd=root,
                text=True,
            ).strip()
        )
        if not re.fullmatch(r"[A-Za-z0-9_+-]+", self.system):
            raise DevError(f"invalid Nix system: {self.system!r}")
        self.installable = f"{source}#lib.{self.system}.development"

    def catalog(self):
        return json.loads(
            subprocess.check_output(
                [
                    "nix",
                    "eval",
                    "--impure",
                    "--no-update-lock-file",
                    "--json",
                    self.installable + ".catalog",
                ],
                cwd=self.root,
                text=True,
            )
        )

    @staticmethod
    @contextmanager
    def request_environment(request):
        with tempfile.NamedTemporaryFile(
            mode="w", prefix="gammaloop-dev-request-", suffix=".json"
        ) as stream:
            json.dump(request, stream)
            stream.flush()
            yield os.environ | {"GAMMALOOP_DEV_REQUEST_FILE": stream.name}

    def prepare_ci(self, state):
        with self.request_environment(state) as environment:
            return json.loads(
                subprocess.check_output(
                    [
                        "nix",
                        "eval",
                        "--impure",
                        "--no-update-lock-file",
                        "--json",
                        self.installable,
                        "--apply",
                        PREPARE_CI_APPLY,
                    ],
                    cwd=self.root,
                    env=environment,
                    text=True,
                )
            )

    def run(self, request, out_link=None):
        with self.request_environment(request) as environment:
            environment |= {
                "GAMMALOOP_DEV_SOURCE": self.source,
                "GAMMALOOP_DEV_SYSTEM": self.system,
            }
            command = [
                "nix",
                "build",
                "--impure",
                "--no-update-lock-file",
                "--expr",
                REQUEST_BUILD_EXPRESSION,
            ]
            if out_link is None:
                command += ["--dry-run", "--no-link"]
            else:
                command += ["--out-link", str(out_link)]
            subprocess.run(
                command, cwd=self.root, env=environment, check=True, stdout=sys.stderr
            )


def read_state(path):
    state = json.loads(path.read_text())
    if (
        not isinstance(state, dict)
        or state.get("version") != 1
        or not isinstance(state.get("contexts"), dict)
    ):
        raise DevError(f"invalid development state: {path}")
    return state


class CacheRegistry:
    """Own rooted results; leases protect seeds while builds run unlocked."""

    def __init__(self, directory, keep=4, maximum=64, max_age_days=30):
        self.directory = Path(directory).expanduser().resolve()
        self.keep = keep
        self.maximum = maximum
        self.max_age = max_age_days * 86400
        self.directory.mkdir(mode=0o700, parents=True, exist_ok=True)
        mode = self.directory.stat().st_mode
        self.shared = bool(mode & stat.S_ISGID and mode & stat.S_IWGRP)
        self.directory_mode = 0o2770 if self.shared else 0o700
        self.file_mode = 0o660 if self.shared else 0o600
        self.mkdir(self.directory / "runs")

    @classmethod
    def from_environment(cls, args):
        directory = (
            args.private_cache
            or os.environ.get("GAMMALOOP_DEV_CACHE")
            or str(
                Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
                / "gammaloop"
                / "development"
            )
        )
        limits = {
            "keepSuccessfulGenerations": args.keep_generations,
            "maxGenerations": args.max_generations,
            "maxAgeDays": args.max_age_days,
        }
        policy_path = os.environ.get("GAMMALOOP_DEV_CACHE_POLICY")
        if policy_path and args.private_cache is None:
            policy = json.loads(Path(policy_path).read_text())
            fields = {"version", "directory", "group", "preserveActiveLeases", *limits}
            if (
                not isinstance(policy, dict)
                or set(policy) != fields
                or type(policy["version"]) is not int
                or policy["version"] != 1
                or policy["preserveActiveLeases"] is not True
            ):
                raise DevError(
                    "invalid shared cache policy; require version 1, active leases and explicit retention limits"
                )
            for name in ("keepSuccessfulGenerations", "maxGenerations"):
                if type(policy[name]) is not int or policy[name] < 1:
                    raise DevError(f"invalid shared cache policy {name}")
            if (
                type(policy["maxAgeDays"]) not in (int, float)
                or not math.isfinite(policy["maxAgeDays"])
                or policy["maxAgeDays"] <= 0
            ):
                raise DevError("invalid shared cache policy maxAgeDays")
            if (
                not isinstance(policy["directory"], str)
                or not Path(policy["directory"]).is_absolute()
            ):
                raise DevError("shared cache policy directory must be absolute")
            if not isinstance(policy["group"], str) or not policy["group"]:
                raise DevError("shared cache policy requires a Unix group")
            configured = Path(policy["directory"])
            if (
                os.environ.get("GAMMALOOP_DEV_CACHE")
                and Path(directory).expanduser().resolve() != configured.resolve()
            ):
                raise DevError(
                    "GAMMALOOP_DEV_CACHE differs from shared policy; use --private-cache for private storage"
                )
            if (
                configured.is_symlink()
                or not configured.is_dir()
                or configured.resolve() != configured
            ):
                raise DevError(
                    "shared cache policy directory must already be a real provisioned directory"
                )
            try:
                group = grp.getgrnam(policy["group"])
            except KeyError as error:
                raise DevError(
                    f"shared cache policy group does not exist: {policy['group']}"
                ) from error
            info = configured.stat()
            policy_info = Path(policy_path).stat()
            if policy_info.st_uid != info.st_uid or policy_info.st_mode & 0o022:
                raise DevError(
                    "shared cache policy must be owned by the cache provisioner and not writable by group/others"
                )
            if info.st_gid != group.gr_gid or stat.S_IMODE(info.st_mode) != 0o2770:
                raise DevError(
                    "shared cache root must match the policy group and mode 2770"
                )
            for name, requested in limits.items():
                if requested is not None and requested != policy[name]:
                    raise DevError(
                        f"shared cache policy owns {name}; per-command override rejected"
                    )
            return cls(
                configured,
                policy["keepSuccessfulGenerations"],
                policy["maxGenerations"],
                policy["maxAgeDays"],
            )
        path = Path(directory).expanduser()
        if path.is_symlink():
            raise DevError("private cache root must be a real account-owned directory")
        if (
            path.exists()
            and path.stat().st_mode & stat.S_ISGID
            and path.stat().st_mode & stat.S_IWGRP
        ):
            raise DevError(
                "shared cache storage requires a matching GAMMALOOP_DEV_CACHE_POLICY; private storage must be private"
            )
        if path.exists() and (
            path.stat().st_uid != os.geteuid()
            or stat.S_IMODE(path.stat().st_mode) & 0o077
        ):
            raise DevError(
                "private cache root must be account-owned with mode 0700; use a dedicated private directory"
            )
        return cls(
            directory,
            args.keep_generations or 4,
            args.max_generations or 64,
            args.max_age_days or 30,
        )

    def mkdir(self, path):
        try:
            path.mkdir(mode=self.directory_mode)
        except FileExistsError:
            if path.is_symlink() or not path.is_dir():
                raise DevError(f"cache directory is not a real directory: {path}")
            if self.shared and path.stat().st_uid == os.geteuid():
                path.chmod(self.directory_mode)
        else:
            # A restrictive umask must not defeat explicit shared provisioning.
            path.chmod(self.directory_mode)

    @contextmanager
    def managed_file(self, path):
        try:
            descriptor = os.open(
                path, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, self.file_mode
            )
        except FileExistsError:
            descriptor = os.open(path, os.O_RDWR | os.O_NOFOLLOW)
            if self.shared and os.fstat(descriptor).st_uid == os.geteuid():
                os.fchmod(descriptor, self.file_mode)
        else:
            os.fchmod(descriptor, self.file_mode)
        with os.fdopen(descriptor, "r+b") as stream:
            yield stream

    @contextmanager
    def locked(self):
        with self.managed_file(self.directory / ".registry.lock") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            yield

    def read(self):
        index = self.directory / "state-index.json"
        if not index.exists():
            return {"version": 1, "generations": []}
        data = json.loads(index.read_text())
        if (
            not isinstance(data, dict)
            or data.get("version") != 1
            or not isinstance(data.get("generations"), list)
        ):
            raise DevError(f"invalid cache registry: {index}")
        for record in data["generations"]:
            if (
                not isinstance(record, dict)
                or not isinstance(record.get("id"), str)
                or not re.fullmatch(r"[0-9a-f]{32}", record["id"])
            ):
                raise DevError(f"invalid generation in cache registry: {index}")
            state = record.get("state")
            if (
                not isinstance(state, dict)
                or state.get("version") != 1
                or not isinstance(state.get("contexts"), dict)
            ):
                raise DevError(f"invalid state in cache registry: {index}")
            if (
                not isinstance(record.get("source"), str)
                or not isinstance(record.get("successful"), bool)
                or not isinstance(record.get("lastUsed"), (int, float))
                or not math.isfinite(record["lastUsed"])
                or not isinstance(record.get("outputs"), list)
                or not all(isinstance(path, str) for path in record["outputs"])
            ):
                raise DevError(f"invalid result metadata in cache registry: {index}")
        return data

    def write(self, index, path=None):
        destination = (
            self.directory / "state-index.json" if path is None else Path(path)
        )
        with tempfile.NamedTemporaryFile(
            mode="w", dir=destination.parent, prefix=".index-", delete=False
        ) as stream:
            temporary = Path(stream.name)
            try:
                json.dump(index, stream)
                stream.flush()
                os.fsync(stream.fileno())
                os.fchmod(stream.fileno(), self.file_mode if path is None else 0o644)
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)

    @contextmanager
    def previous(self, successful_only=False):
        with ExitStack() as leases:
            state = {"version": 1, "contexts": {}}
            with self.locked():
                for record in self.read()["generations"]:
                    if successful_only and not record.get("successful", False):
                        continue
                    lease = leases.enter_context(
                        (self.directory / "runs" / record["id"] / ".lease").open("rb")
                    )
                    fcntl.flock(lease, fcntl.LOCK_SH)
                    # Context descriptors are validated by Nix, never reconstructed here.
                    state["contexts"].update(record["state"]["contexts"])
            yield state if state["contexts"] else {}

    @contextmanager
    def generation(self):
        path = self.directory / "runs" / uuid.uuid4().hex
        path.mkdir(mode=self.directory_mode)
        path.chmod(self.directory_mode)
        with self.managed_file(path / ".lease") as lease:
            fcntl.flock(lease, fcntl.LOCK_SH)
            try:
                yield path
            finally:
                # Publication transfers ownership to the registry; failures own no roots.
                with self.locked():
                    retained = any(
                        record["id"] == path.name
                        for record in self.read()["generations"]
                    )
                    if not retained:
                        shutil.rmtree(path)

    def publish(self, generation, source, results, successful=True):
        state = {"version": 1, "contexts": {}}
        for result in results:
            state["contexts"].update(
                read_state(result / "development-state.json")["contexts"]
            )
        outputs = [str(result.resolve(strict=True)) for result in results]
        with self.locked(), ExitStack() as pruning:
            index = self.read()
            for record in index["generations"]:
                if (
                    record["source"] == source
                    and record["state"] == state
                    and record.get("outputs") == outputs
                ):
                    if successful and not record.get("successful", False):
                        record["successful"] = True
                    record["lastUsed"] = time.time()
                    index["generations"] = [
                        item
                        for item in index["generations"]
                        if item["id"] != record["id"]
                    ] + [record]
                    break
            else:
                index["generations"].append(
                    {
                        "id": generation.name,
                        "source": source,
                        "state": state,
                        "outputs": outputs,
                        "successful": successful,
                        "lastUsed": time.time(),
                    }
                )
            removable = []
            retained = []
            # Keys are opaque compatibility families owned by Nix. Keep recent
            # successful seeds separately from failed diagnostic generations.
            counts = {}
            wanted = set()
            now = time.time()
            for record in reversed(index["generations"]):
                families = [
                    (record["successful"], key) for key in record["state"]["contexts"]
                ] or [(record["successful"], None)]
                if (
                    len(wanted) < self.maximum
                    and now - record["lastUsed"] <= self.max_age
                    and any(counts.get(family, 0) < self.keep for family in families)
                ):
                    wanted.add(record["id"])
                    for family in families:
                        counts[family] = counts.get(family, 0) + 1
            for record in index["generations"]:
                if record["id"] in wanted:
                    retained.append(record)
                    continue
                path = self.directory / "runs" / record["id"]
                lease = pruning.enter_context((path / ".lease").open("rb"))
                try:
                    fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    retained.append(record)
                else:
                    removable.append(path)
            index["generations"] = retained
            self.write(index)
            for path in removable:
                shutil.rmtree(path)


def junit_count(result):
    reports = list(result.glob("*.xml"))
    if not reports:
        raise DevError(f"missing JUnit reports in {result}")
    count = 0
    for report in reports:
        try:
            root = ET.parse(report).getroot()
        except ET.ParseError as error:
            raise DevError(f"invalid JUnit report: {report}: {error}") from error
        if any(
            node.tag.rsplit("}", 1)[-1] in ("failure", "error") for node in root.iter()
        ):
            raise DevError(f"JUnit report contains test failures: {report}")
        count += sum(node.tag.rsplit("}", 1)[-1] == "testcase" for node in root.iter())
    return count


def test_outcome(results, policy):
    statuses = []
    for result in results:
        path = result / "exit-status"
        if not path.is_file():
            raise DevError(f"missing development test exit status: {path}")
        text = path.read_text().strip()
        if not re.fullmatch(r"[0-9]+", text) or int(text) > 255:
            raise DevError(f"invalid development test exit status: {path}")
        statuses.append(int(text))
    for result, status in zip(results, statuses):
        if status not in (0, 4):
            print(
                f"Tests failed (exit {status}); diagnostics: {result}", file=sys.stderr
            )
            return status
    count = sum(junit_count(result) for result in results)
    print(f"JUnit: {count} test cases", file=sys.stderr)
    if count == 0:
        message = "no tests matched the overall request"
        if policy == "fail":
            print(message, file=sys.stderr)
            return 4
        if policy == "warn":
            print(f"Warning: {message}", file=sys.stderr)
    return 0


def link_logs(directory, generation, requests, results):
    destination = Path(directory).expanduser().resolve() / generation.name
    destination.mkdir(parents=True, exist_ok=False)
    for request, result in zip(requests, results):
        logs = result / "result" / "logs"
        if not logs.is_dir():
            raise DevError(f"missing requested test logs: {logs}")
        link = destination / request["group"]
        link.symlink_to(logs.resolve(), target_is_directory=True)
        print(f"Logs: {link}", file=sys.stderr)


def execute(args, root):
    source = SourceSnapshot(root).store()
    development = NixDevelopment(root, source, args.system)
    requests = requests_for(args, development.catalog())
    registry = CacheRegistry.from_environment(args)
    print(f"Source: {source}", file=sys.stderr)
    if args.dry_run:
        with registry.previous() as previous:
            for request in requests:
                print(json.dumps(request, sort_keys=True), file=sys.stderr)
                development.run(request | {"previousState": previous})
        return 0
    with registry.generation() as generation:
        results = []
        products = []
        outcome = None
        try:
            with registry.previous() as previous:
                for index, request in enumerate(requests):
                    out_link = generation / f"result-{index}"
                    print(
                        f"{request['command']}: {', '.join(request['packages'])}"
                        + (f" [{request['group']}]" if "group" in request else ""),
                        file=sys.stderr,
                    )
                    development.run(request | {"previousState": previous}, out_link)
                    result = out_link / "result"
                    if not result.exists():
                        raise DevError(f"missing development result: {result}")
                    results.append(out_link)
                    if args.command == "build":
                        for package in request["packages"]:
                            product = result / package
                            if not product.exists():
                                raise DevError(
                                    f"missing build product for {package}: {product}"
                                )
                            products.append(product)
                    else:
                        products.append(result)
            outcome = (
                test_outcome([result / "result" for result in results], args.no_tests)
                if args.command == "test"
                else 0
            )
            return outcome
        finally:
            # Completed artifacts survive a later group/setup error or malformed
            # report. Seed leases have ended, so safe pruning can proceed here.
            if results:
                active_error = sys.exc_info()[1]
                try:
                    registry.publish(
                        generation, source, results, successful=outcome == 0
                    )
                    for product in products:
                        print(product.resolve())
                    if args.command == "test" and (
                        directory := os.environ.get("GL_TEST_LOG_DIR")
                    ):
                        link_logs(
                            directory, generation, requests[: len(results)], results
                        )
                except Exception as error:
                    if active_error is None and outcome == 0:
                        raise
                    print(f"Could not retain all diagnostics: {error}", file=sys.stderr)


def prepare_ci(args, root):
    registry = CacheRegistry.from_environment(args)
    manifest = (Path(root) / args.manifest).resolve()
    # Serialize promotions into this workspace across different cache accounts.
    # Directory flock needs no lock artifact in the captured source tree.
    lock = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
    try:
        fcntl.flock(lock, fcntl.LOCK_EX)
        with registry.previous(successful_only=True) as available:
            if not available or not available["contexts"]:
                raise DevError("no successful rooted compiler state to promote")
            names = set(args.context) if args.context else set(available["contexts"])
            if missing := names - available["contexts"].keys():
                raise DevError(
                    f"no successful rooted state for contexts: {', '.join(sorted(missing))}"
                )
            state = {
                "version": 1,
                "contexts": {
                    name: available["contexts"][name] for name in sorted(names)
                },
            }
            source = SourceSnapshot(root).store()
            development = NixDevelopment(root, source, args.system)
            plan = development.prepare_ci(state)
            if (
                not isinstance(plan, dict)
                or not isinstance(plan.get("state"), dict)
                or not isinstance(plan.get("paths"), list)
            ):
                raise DevError("invalid development.prepareCi response")
            promoted = plan["state"]
            if promoted.get("version") != 1 or not isinstance(
                promoted.get("contexts"), dict
            ):
                raise DevError("invalid promoted state from development.prepareCi")
            for name, descriptor in promoted["contexts"].items():
                if (
                    name not in state["contexts"]
                    or descriptor != state["contexts"][name]
                ):
                    raise DevError(
                        "development.prepareCi changed an opaque context descriptor"
                    )
            for omitted in sorted(names - promoted["contexts"].keys()):
                print(
                    f"Omitted {omitted}: no identical selected NixCI consumer",
                    file=sys.stderr,
                )
            if not promoted["contexts"]:
                raise DevError(
                    "selected state has no identical selected NixCI consumers"
                )
            paths = plan["paths"]
            if not paths or not all(
                isinstance(path, str)
                and re.fullmatch(r"/nix/store/[0-9a-z]{32}-[^/\n\0]+", path)
                for path in paths
            ):
                raise DevError("development.prepareCi returned invalid closure roots")
            # State is opaque, but a license field is never a publishable descriptor.
            if '"SYMBOLICA_LICENSE"' in json.dumps(promoted):
                raise DevError(
                    "license fields are forbidden in development seed manifests"
                )
            if args.dry_run:
                print(
                    json.dumps(
                        {
                            "manifest": str(manifest),
                            "contexts": sorted(promoted["contexts"]),
                            "paths": paths,
                        },
                        indent=2,
                    )
                )
                return 0
            netrc = Path(
                os.environ.get("NIXCI_NETRC", str(Path.home() / ".netrc"))
            ).expanduser()
            if not netrc.is_file() or not os.access(netrc, os.R_OK):
                raise DevError(
                    "set up ~/.netrc or NIXCI_NETRC before publishing compiler state"
                )
            with tempfile.TemporaryDirectory(
                prefix="gammaloop-dev-publish-"
            ) as temporary:
                # Avoid remembered unsigned narinfo hiding newly signed cache paths.
                environment = os.environ | {
                    "XDG_CACHE_HOME": temporary,
                    "NIX_CACHE_HOME": str(Path(temporary) / "nix"),
                }
                subprocess.run(
                    [
                        "nix",
                        "copy",
                        "--netrc-file",
                        str(netrc),
                        "--to",
                        args.destination,
                        "--stdin",
                    ],
                    input="\n".join(dict.fromkeys(paths)) + "\n",
                    text=True,
                    cwd=root,
                    env=environment,
                    check=True,
                    stdout=sys.stderr,
                )
            if SourceSnapshot(root).store() != source:
                raise DevError(
                    "workspace changed during publication; seed manifest was not written"
                )
            manifest.parent.mkdir(parents=True, exist_ok=True)
            registry.write(promoted, manifest)
            print(manifest)
            return 0
    finally:
        os.close(lock)


def main(argv=None):
    try:
        args = parse_args(sys.argv[1:] if argv is None else argv)
        root = Path(__file__).resolve().parents[1]
        return (
            prepare_ci(args, root)
            if args.command == "prepare-ci"
            else execute(args, root)
        )
    except subprocess.CalledProcessError as error:
        print(f"Development command failed (exit {error.returncode})", file=sys.stderr)
        return error.returncode if 0 < error.returncode < 256 else 1
    except (DevError, OSError, ValueError) as error:
        print(f"Development command failed: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Development command interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
