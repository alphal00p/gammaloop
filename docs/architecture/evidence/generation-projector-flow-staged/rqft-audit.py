"""Certify the exact source26 fixture patch without evaluating graph numerators."""

import hashlib
import json
from pathlib import Path
import re
import subprocess

path = "tests/resources/graphs/uv_tests/rqft_ghG_3l.dot"
source = "8cf71575a739976295edb42eae2052d6b5e455ee"
base = "d48f9d1aedc766eea5cb185ab6029d1f04462043"
before = subprocess.check_output(["jj", "file", "show", "-r", base, path]).decode()
after = Path(path).read_text()
original = subprocess.check_output(["jj", "file", "show", "-r", source, path]).decode()
assert original == after
assert re.sub(r"^//.*\n", "", re.sub(r"lmb_id=\d+ ", "", after), flags=re.M) == before
assert len(re.findall(r"^digraph ", after, re.M)) == 78
pins = {}
for block in re.split(r"(?=digraph\s)", after):
    match = re.match(r"digraph (\w+)", block)
    if match and "lmb_id=" in block:
        pins[match[1]] = re.findall(
            r'(\w+) -> (\w+) \[lmb_id=(\d+) particle= "([^"]+)"\]', block
        )
assert set(pins) == {"d37", "d42", "d43", "d44", "d45", "d46", "d47"}
assert sum(map(len, pins.values())) == 21
record = {
    "source": source,
    "source_fixture_byte_identical": True,
    "physical_preservation": "removing only the 21 lmb_id attributes and three new comment lines reproduces d48 fixture byte-for-byte",
    "graphs": 78,
    "pins": {
        name: [
            {
                "source": source,
                "target": target,
                "loop_slot": int(slot),
                "particle": particle,
            }
            for source, target, slot, particle in edges
        ]
        for name, edges in pins.items()
    },
    "edge_selection_contract": "exact authored source/target/particle and lmb_id identify the selected edge through the current source map; parsed native edge IDs need not equal DOT statement ordinals",
    "expected_sign": "positive unit loop coefficient, zero external coefficients on each pinned source direction",
    "before_sha256": hashlib.sha256(before.encode()).hexdigest(),
    "after_sha256": hashlib.sha256(after.encode()).hexdigest(),
}
Path(
    "docs/architecture/evidence/generation-projector-flow-staged/rqft-bases.json"
).write_text(json.dumps(record, indent=2) + "\n")
print(
    "78 graphs, seven exact directed bases, 21 pins; all physical input bytes preserved"
)
