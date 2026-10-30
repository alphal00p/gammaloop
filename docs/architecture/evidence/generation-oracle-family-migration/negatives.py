"""Actual zero-selection and record-mode executions have zero pass credit."""

import re

import common as M
from runtime import OLD, sandbox


def main():
    rows = []
    for label, elf, selection, flags in (
        (
            "zero-selection",
            M.OUT / "bin/generation_oracle",
            "migration_nonexistent_exact_selection",
            {},
        ),
        ("record-mode-old-elf", OLD, M.CASES[0], {"FEYNKIT_ORACLE_RECORD": "1"}),
    ):
        root = M.OUT / "negatives" / label
        proc, command = sandbox(
            elf,
            [selection, "--exact", "--test-threads", "1", "--nocapture"],
            root,
            {"FEYNKIT_ORACLE_THREADS": "1", **flags},
        )
        found = re.findall(
            r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out",
            proc.stdout,
        )
        M.require(proc.returncode == 0 and len(found) == 1, "negative-setup", label)
        selected = int(found[0][0]) + int(found[0][1])
        if label == "zero-selection":
            M.require(selected == 0, "negative-setup", label)
            reason = "actual libtest exit0 with selected0 violates the ordinary matrix's explicit nonzero selection predicate"
        else:
            M.require(selected == 1 and "Row {" in proc.stdout, "negative-setup", label)
            reason = "actual old ELF exits0 in record mode by skipping the Row assertion; RECORD presence disqualifies pass credit independently of selection"
        rows.append(
            {
                "name": label,
                "command": command,
                "selected": selected,
                "admitted": False,
                "reason": reason,
                "credit": "none; diagnostic counterexample to treating exit0 as test acceptance",
            }
        )
    M.write(M.OUT / "negatives/results.json", rows)


if __name__ == "__main__":
    main()
