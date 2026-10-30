"""New migration records; frozen certificate helpers and outputs stay read-only."""

import importlib.util
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[3]
OUT = WORKSPACE / "target/xs-qed-l2-migration"
APPROVED = "c95e7f6e77994779929115a156654cfce3bd1e53"
FROZEN = HERE.parent / "generation-xs-qed-l2"
SPEC = importlib.util.spec_from_file_location(
    "frozen_record_design", FROZEN / "common.py"
)
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
M.HERE, M.OUT, M.APPROVED = HERE, OUT, APPROVED
sha, load, write, require, C = M.sha, M.load, M.write, M.require, M.C
MAIN = M.MAIN
CASES = ("xs_qed_eemm_L2",)
REVISIONS = M.REVISIONS
URI = f"git+file:///shared/lcnbr/dev/gammaloop/main?rev={APPROVED}&shallow=1&allRefs=1"


def run(argv, directory, *, cwd=WORKSPACE, env_extra=None, inputs=()):
    argv = [str(a) for a in argv]
    if Path(argv[0]).name == "nix":
        require(
            len(argv) > 2
            and argv[1] == "develop"
            and argv[2] == URI
            and re.fullmatch(r"[0-9a-f]{40}", APPROVED),
            "nix-installable",
            "explicit immutable Git URI required; no default/private flake",
        )
    return M.run(
        argv,
        directory,
        cwd=cwd,
        env_extra=env_extra,
        inputs=[
            *[
                FROZEN / n
                for n in (
                    "common.py",
                    "boundary.py",
                    "bind_runtime.py",
                    "build.py",
                    "fingerprints.py",
                )
            ],
            HERE.parent / "generation-oracle-family-migration/derive.py",
            *C.source_hashes(),
            *inputs,
        ],
    )
