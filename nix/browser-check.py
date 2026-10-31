"""Run selected native browser producers and validate their fresh JUnit evidence."""

import argparse
import os
import shutil
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from contextlib import nullcontext
from pathlib import Path


def selected_browser_checks(package, report):
    """Select clients by Nextest's binary ID and successful, executed test name."""
    producers = {
        ("spynso3", "spynso3::tensor_paging", "direct_tensor_paging"): "tensor",
        (
            "feynkit-py",
            "feynkit-py::cff_display",
            "cff_explorer_keeps_native_families_and_graph_identity",
        ): "cff",
        (
            "feynkit-py",
            "feynkit-py::combined_community",
            "browser_display_fixtures",
        ): "feynkit",
        (
            "linnet-render-py",
            "linnet-render-py::mathml_view",
            "mathml_view_browser_producer",
        ): "mathml",
    }
    root = ET.parse(report).getroot()
    if root.tag != "testsuites" or root.get("name") != "nextest-run":
        raise ValueError(f"Not a Nextest JUnit report: {report}")
    selected = set()
    totals = dict.fromkeys(("tests", "skipped", "failures", "errors"), 0)
    for suite in root.findall("testsuite"):
        binary = suite.attrib["name"]
        if binary != package and not binary.startswith(package + "::"):
            raise ValueError(f"Unexpected binary {binary!r} in {report}")
        cases = suite.findall("testcase")
        counts = {
            "tests": len(cases),
            "skipped": sum(case.find("skipped") is not None for case in cases),
            "failures": sum(case.find("failure") is not None for case in cases),
            "errors": sum(case.find("error") is not None for case in cases),
        }
        for key, count in counts.items():
            if int(suite.attrib[key]) != count:
                raise ValueError(f"Inconsistent JUnit {key} count in {report}")
            totals[key] += count
        for case in cases:
            name = case.attrib["name"]
            if not name or case.attrib["classname"] != binary:
                raise ValueError(f"Invalid JUnit test identity in {report}")
            if case.find("skipped") is None:
                check = producers.get((package, binary, name))
                if check is not None:
                    selected.add(check)
    if any(int(root.attrib[key]) != count for key, count in totals.items()):
        raise ValueError(f"Inconsistent JUnit totals in {report}")
    if totals["failures"] or totals["errors"]:
        raise ValueError(f"Native tests failed in {report}")
    return tuple(
        check for check in ("tensor", "cff", "feynkit", "mathml") if check in selected
    )


def run_owned(command, package=None, browser_report=None):
    """Forward cancellation and reap each command, optionally coordinating MathML."""
    output = None
    environment = dict(os.environ)
    if package == "linnet-render-py":
        output = Path(environment["LINNET_MATHML_BROWSER_ROOT"]) / "ci"
        if output.exists():
            raise ValueError(f"Browser producer output must be fresh: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        environment["LINNET_MATHML_BROWSER_DIR"] = str(output)
    native = None
    client = None
    client_ran = False
    cancelled = None
    pending_signals = []
    previous = {
        signum: signal.getsignal(signum) for signum in (signal.SIGTERM, signal.SIGINT)
    }

    def interrupted(signum, _frame):
        nonlocal cancelled
        if cancelled is None:
            cancelled = signum
        if native is None:
            pending_signals.append(signum)
        for child in (native, client):
            if child is not None:
                try:
                    child.send_signal(signum)
                except ProcessLookupError:
                    pass

    log_file = (
        (output.parent / "producer.log").open("w", encoding="utf-8")
        if output
        else nullcontext()
    )
    with log_file as log:
        for signum in previous:
            signal.signal(signum, interrupted)
        try:
            # Isolate terminal signals so each cancellation is forwarded once.
            # Nextest then supervises its separate groups of running tests.
            native = subprocess.Popen(
                command,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT if log else None,
                start_new_session=True,
            )
            for signum in pending_signals:
                native.send_signal(signum)
            if output is not None:
                ready = output / "ready.json"
                deadline = None
                while native.poll() is None and cancelled is None:
                    if ready.exists():
                        command = [
                            sys.executable,
                            "crates/linnet-render-py/tests/mathml_view_browser.py",
                            "--ready",
                            str(ready),
                            "--output",
                            str(output / "client"),
                            "--engines",
                            "chromium",
                        ]
                        client = subprocess.Popen(
                            command, env=environment, start_new_session=True
                        )
                        if cancelled is not None:
                            client.send_signal(cancelled)
                        status = client.wait(timeout=230)
                        if status:
                            raise subprocess.CalledProcessError(status, command)
                        client_ran = True
                        break
                    # Other selected tests must not consume the readiness budget.
                    # serve_browser creates this directory when it starts.
                    if output.exists():
                        if deadline is None:
                            deadline = time.monotonic() + 30
                        if time.monotonic() >= deadline:
                            raise RuntimeError(
                                "native browser producer failed to become ready"
                            )
                    time.sleep(0.05)
            return native.wait(), client_ran
        finally:
            try:
                try:
                    if client is not None:
                        if client.poll() is None:
                            client.send_signal(signal.SIGTERM)
                        client.wait()
                finally:
                    if native is not None:
                        if native.poll() is None and cancelled is None:
                            native.send_signal(signal.SIGTERM)
                        # Nextest controls test grace periods and escalation.
                        # Never kill its supervisor while it may be reaping tests.
                        native.wait()
            finally:
                for signum, handler in previous.items():
                    signal.signal(signum, handler)
                if output is not None:
                    if output.exists():
                        shutil.copytree(output.parent, browser_report)
                    print(log.name)
                    print(Path(log.name).read_text(encoding="utf-8"))
                if cancelled is not None:
                    raise SystemExit(128 + cancelled)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--package",
        required=True,
        choices=("spynso3", "feynkit-py", "linnet-render-py"),
    )
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--browser-report", required=True, type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.command or args.command[0] != "--" or len(args.command) == 1:
        parser.error("expected -- followed by the native Nextest command")
    # Pinned Nextest 0.9.143 binds neither stress option to an environment
    # variable. Tokens after the native separator are libtest filters.
    for argument in args.command[1:]:
        if argument == "--":
            break
        if argument.partition("=")[0] in ("--stress-count", "--stress-duration"):
            parser.error(
                "browser runtime does not support Nextest stress-count or stress-duration modes"
            )

    status, client_ran = run_owned(args.command[1:], args.package, args.browser_report)
    if status != 0:
        return status if status > 0 else 128 - status

    selected = selected_browser_checks(args.package, args.report)
    if args.package == "linnet-render-py" and ("mathml" in selected) != client_ran:
        raise RuntimeError("MathML client and executed producer evidence disagree")
    commands = []
    if "tensor" in selected:
        fixtures = os.environ["SPENSO_MATHML_BROWSER_FIXTURES"]
        chromium = os.environ["LINNET_TEST_CHROMIUM"]
        commands.extend(
            [
                [
                    "crates/spynso3/tests/browser_math_font.py",
                    "--fixtures",
                    fixtures,
                    "--chromium",
                    chromium,
                    "--font",
                    str(
                        Path(os.environ["SYMBOLICA_DATA_DIRECTORY"])
                        / "render/fonts/STIXTwoMath-Regular.woff2"
                    ),
                ],
                [
                    "crates/spynso3/tests/browser_tensor_explorer.py",
                    "--fixtures",
                    fixtures,
                    "--chromium",
                    chromium,
                ],
            ]
        )
    if "cff" in selected:
        commands.append(
            [
                "crates/feynkit-py/tests/cff_display_browser.py",
                os.environ["CFF_DISPLAY_TEST_OUTPUT"],
                os.environ["LINNET_TEST_CHROMIUM"],
                "chromium",
            ]
        )
    if "feynkit" in selected:
        fixtures = Path(os.environ["FEYNKIT_BROWSER_FIXTURES"])
        chromium = os.environ["LINNET_TEST_CHROMIUM"]
        commands.extend(
            [
                [
                    "crates/feynkit-py/tests/ltd_display_browser.py",
                    str(fixtures / "ltd"),
                    chromium,
                    "chromium",
                ],
                [
                    "crates/feynkit-py/tests/energy_scoped_display_browser.py",
                    str(fixtures / "scoped"),
                    chromium,
                    "chromium",
                ],
                [
                    "crates/feynkit-py/tests/generalized_energy_display_browser.py",
                    str(fixtures / "generalized"),
                    str(fixtures / "ltd"),
                    chromium,
                    "chromium",
                ],
            ]
        )
    for command in commands:
        command_status, _ = run_owned([sys.executable, *command])
        if command_status:
            status = command_status if command_status > 0 else 128 - command_status
    return status


if __name__ == "__main__":
    sys.exit(main())
