"""Immutable source construction and finite pre-capture sizing, no physics edits."""

import hashlib
import io
import subprocess
import tarfile

import common as M


def main():
    paths = (
        M.MAIN,
        "crates/feynkit-generator/tests/generation_oracle/cases.rs",
        "crates/feynkit-generator/tests/generation_oracle/classes.rs",
        "crates/feynkit-generator/tests/generation_oracle/census.rs",
        "crates/feynkit-graph/src/symbols.rs",
        "crates/idenso/src/cook.rs",
        "crates/idenso/src/lib.rs",
        "Cargo.lock",
        "crates/feynkit-model/src/lib.rs",
    )
    records = []
    for route, revision in M.REVISIONS.items():
        for path in paths:
            text, record = M.immutable(revision, path)
            records.append(
                {
                    "route": route,
                    "revision": revision,
                    "path": path,
                    "sha256": hashlib.sha256(text.encode()).hexdigest(),
                    "command": record,
                }
            )
        root = M.OUT / "source" / route
        root.mkdir(parents=True, exist_ok=False)
        # Git immutable object reads only; no physical MAIN files/index touched.
        argv = [
            "git",
            "--git-dir=" + str(M.WORKSPACE.parent / "main/.git"),
            "archive",
            revision,
        ]
        env = M.clean_env()
        env["GIT_OPTIONAL_LOCKS"] = "0"
        proc = subprocess.run(argv, env=env, capture_output=True, check=False)
        M.require(proc.returncode == 0, "archive", route)
        with tarfile.open(fileobj=io.BytesIO(proc.stdout)) as archive:
            archive.extractall(root, filter="data")
        M.write(
            M.OUT / "source" / (route + "-archive.json"),
            {
                "argv": argv,
                "exit": proc.returncode,
                "archive_sha256": hashlib.sha256(proc.stdout).hexdigest(),
                "revision": revision,
                "policy": "no license variables; Git object archive; GIT_OPTIONAL_LOCKS0; pristine before observer",
            },
        )
    M.write(M.OUT / "source-bindings.json", records)
    M.write(
        M.OUT / "pre-capture.json",
        {
            "case": M.CASES[0],
            "expected_retained_graphs": 30,
            "topologies": 157,
            "interaction_assignments": 714,
            "grouping": "None",
            "modes": ["None", "Auto"],
            "pool": 1,
            "source_mechanism": "Model::qed; EE->MUMU; forward total loops2; Cuts blobs1 spectators0; plain; self loops false; zero-color filter false",
            "sizing": "Row30 retained graphs, source157 topologies714 assignments; 6 trivalent vertices and11 total edges expected from forward loop/Euler incidence, so20 owners per graph=600. This is a pre-capture estimate checked against actual output, not a proof or allocation budget.",
            "new_obligations": [
                "physical cuts and cut-edge numerator owners versus amplitude domains",
                "all source/group singleton chains",
                "full None/Auto occurrence join",
                "routing/external incidence and cut signs",
                "cross-section projector and live typed interfaces",
                "local/global products with sums whole",
            ],
            "stages": [
                "case/model/resolved request",
                "raw assignment and cuts",
                "historical cooked vs native typed tensor admission",
                "factorized local/global product",
                "Grouping None census",
            ],
            "first_unobserved_boundary": "resolved request to actual raw graph/rule/cut identity; do not assign cause until complete pipeline checks",
        },
    )


if __name__ == "__main__":
    main()
