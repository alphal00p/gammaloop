from pathlib import Path
import subprocess, re, json, hashlib

root = Path("docs/architecture/evidence/generation-projector-flow-staged")
path = "tests/resources/graphs/qqx_aaa_subtracted.dot"
read = lambda rev, p: subprocess.check_output(
    ["jj", "file", "show", "-r", rev, p]
).decode()
before = read("d48f9d1aedc766eea5cb185ab6029d1f04462043", path)
after = Path(path).read_text()
old = "spenso::g(spenso::cof(3,hedge_1),spenso::dind(spenso::cof(3,hedge_2)))"
new = "spenso::g(spenso::dind(spenso::cof(3,hedge_1)),spenso::cof(3,hedge_2))"
assert before.count(old) == 8
expected = before.replace(old, new).replace("group_master=", "is_group_master=")
inverse = re.sub(r"gammalooprs::hedge\((\d+),1\)", r"hedge_\1", after)
inverse = re.sub(r"gammalooprs::edge\((\d+),(\d+)\)", r"edge_\1_\2", inverse)
assert inverse == expected
approved = read("735d599a665acbce1577b2c40d149c152f01543d", path)
stripped = re.sub(r'(?:source|sink)="\{ufo_order:\d+\}", ', "", approved)
assert stripped == after
result = {
    "staged_base": "d48f9d1aedc766eea5cb185ab6029d1f04462043",
    "old_approved": "735d599a665acbce1577b2c40d149c152f01543d",
    "graphs": len(re.findall("^digraph", after, re.M)),
    "exact_variance_corrections": 8,
    "staged17_delta": "only already-removed source/sink ufo_order metadata; never restored",
    "factor_routing_weight_master_preservation": "inverse encoding exact byte equality against staged base after eight flips and master-key spelling",
    "old_preparation_equivalence": "exact byte equality after removing obsolete source/sink ufo_order metadata",
    "before_sha256": hashlib.sha256(before.encode()).hexdigest(),
    "after_sha256": hashlib.sha256(after.encode()).hexdigest(),
}
assert result["graphs"] == 19
(root / "fixture-reconciliation.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))

# The staged setup already makes topology filters and grouping explicit. Keep
# those two command lines while certifying every remaining bubble-test byte.
bubble_path = "tests/tests/test_runs/test_integrated_uv_cts.rs"
old_bubble = read("735d599a665acbce1577b2c40d149c152f01543d", bubble_path)
staged_bubble = read("d48f9d1aedc766eea5cb185ab6029d1f04462043", bubble_path)
new_bubble = Path(bubble_path).read_text()
anchor = "const INSPECT_POINT:"
old_bubble = old_bubble[old_bubble.index(anchor) :]
new_bubble = new_bubble[new_bubble.index(anchor) :]
expected_bubble = old_bubble
commands = []
for line in old_bubble.splitlines(True):
    if '"generate amp scalar_1 > scalar_1 [{1}]' in line:
        candidates = [
            item for item in staged_bubble.splitlines(True) if line.strip()[:-2] in item
        ]
        assert len(candidates) == 1
        expected_bubble = expected_bubble.replace(line, candidates[0])
        commands.append({"old": line.strip(), "staged": candidates[0].strip()})
assert len(commands) == 2
assert expected_bubble == new_bubble
(root / "bubble-preservation.json").write_text(
    json.dumps(
        {
            "old_approved": "735d599a665acbce1577b2c40d149c152f01543d",
            "path": bubble_path,
            "scope": "entire file from INSPECT_POINT through final test",
            "numeric_formulas_coordinates_precision_tolerances_assertions_byte_identical": True,
            "only_staged_delta": "two existing d48 generation commands with explicit filters/grouping; imports use current staged organization",
            "staged_commands": commands,
            "actual_sha256": hashlib.sha256(new_bubble.encode()).hexdigest(),
        },
        indent=2,
    )
    + "\n"
)
print(
    "bubble formulas, coordinates, precision, tolerances and assertions retained exactly"
)
