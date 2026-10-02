"""Transport routing and grouping into complete native runtime DOT fixtures."""

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from itertools import permutations, product
from pathlib import Path

from symbolica import core
from symbolica.community import feynkit as fk

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--generated",
    type=Path,
    required=True,
    help="Directory of complete native runtime DOT files",
)
parser.add_argument(
    "--model",
    type=Path,
    required=True,
    help="model.json exported by the same CLI generation process",
)
parser.add_argument(
    "--historical",
    type=Path,
    required=True,
    help="Directory of original GL000/GL053/GL148 DOT files",
)
parser.add_argument("--generation-card", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
ROOT, OLD, MODEL, OUT = args.generated, args.historical, args.model, args.output
OUT.mkdir(parents=True, exist_ok=False)
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
model = fk.Model.from_json(MODEL.read_text())


def topology(text):
    # Interpret native IDs by the exact exported model. This projection is
    # solely for checking; production imports the original complete syntax.
    text = re.sub(
        r'particle_id="(\d+)"',
        lambda m: 'particle="' + model.particles[int(m[1])].name + '"',
        text,
    )
    text = re.sub(
        r'vertex_rule_id="(\d+)"',
        lambda m: 'int_id="' + model.vertex_rules[int(m[1])].name + '"',
        text,
    )
    return json.loads(fk.FeynmanDiagram.from_dot(model, text).to_json())


def attribute(text, name):
    return re.search(r"\b" + name + r'\s*=\s*"([^"\\]*(?:\\.[^"\\]*)*)"', text).group(1)


def key(ends, data, mapping, external_labels):
    ext = data["external"]
    external = (
        "e:" + ext["state"] + (":" + str(ext["index"]) if external_labels else "")
        if ext
        else None
    )
    a, b = [
        external if ends[k] is None else "v:" + str(mapping[ends[k]])
        for k in ["source", "target"]
    ]
    if not data["directed"]:
        a, b = sorted([a, b])
    return (a, b, data["particle"], data["directed"])


def maps(a, b, labels):
    if Counter(v["interaction"] for v in a["vertices"]) != Counter(
        v["interaction"] for v in b["vertices"]
    ):
        return
    aa = defaultdict(list)
    bb = defaultdict(list)
    for i, v in enumerate(a["vertices"]):
        aa[v["interaction"]].append(i)
    for i, v in enumerate(b["vertices"]):
        bb[v["interaction"]].append(i)
    target = Counter(
        key(ends, data, range(len(b["vertices"])), labels) for ends, data in b["edges"]
    )
    groups = sorted(aa)
    for options in product(*(permutations(bb[g]) for g in groups)):
        mapping = {
            i: j
            for g, positions in zip(groups, options, strict=True)
            for i, j in zip(aa[g], positions, strict=True)
        }
        if (
            Counter(key(ends, data, mapping, labels) for ends, data in a["edges"])
            == target
        ):
            yield mapping


generated = [
    (path, topology(path.read_text())) for path in sorted(ROOT.glob("GL*.dot"))
]
matches = {}
for name in ("GL000", "GL053", "GL148"):
    old = topology((OLD / (name + ".dot")).read_text())
    found = []
    for path, candidate in generated:
        for mapping in maps(old, candidate, True):
            positions = defaultdict(list)
            for index, (ends, data) in enumerate(candidate["edges"]):
                positions[
                    key(ends, data, range(len(candidate["vertices"])), True)
                ].append(index)
            edges = {
                i: positions[key(ends, data, mapping, True)].pop(0)
                for i, (ends, data) in enumerate(old["edges"])
            }
            found.append(
                {
                    "candidate": path.name,
                    "vertices": mapping,
                    "edges": edges,
                    "fixed_external_indices": True,
                    "transported_lmb": [
                        edges[e] for e in old["loop_momentum_basis"]["loop_edges"]
                    ],
                }
            )
    if len(found) != 1:
        raise ValueError(
            f"Expected one exact fixed-external topology match for {name}, found {len(found)}"
        )
    matches[name] = {"matches": found}
records = {}
for name, entry in matches.items():
    assert len(entry["matches"]) == 1
    match = entry["matches"][0]
    assert match["fixed_external_indices"]
    old_path = OLD / (name + ".dot")
    old_dot = old_path.read_text()
    old = topology(old_dot)
    generated_path = ROOT / match["candidate"]
    generated_dot = generated_path.read_text()
    original = topology(generated_dot)
    vertices = {int(a): b for a, b in match["vertices"].items()}
    edges = {int(a): b for a, b in match["edges"].items()}
    reverse_edges = {b: a for a, b in edges.items()}
    orientations = {}
    for a, b in edges.items():
        old_ends, old_data = old["edges"][a]
        new_ends, new_data = original["edges"][b]
        mapped = {
            key: None if value is None else vertices[value]
            for key, value in old_ends.items()
        }
        if mapped == new_ends:
            orientations[a] = 1
        else:
            assert not old_data["directed"]
            assert mapped == {
                "source": new_ends["target"],
                "target": new_ends["source"],
            }
            orientations[a] = -1
        assert old_data["particle"] == new_data["particle"]
        assert (old_data["external"] is None) == (new_data["external"] is None)
        if old_data["external"] is not None:
            assert {k: v for k, v in old_data["external"].items() if k != "name"} == {
                k: v for k, v in new_data["external"].items() if k != "name"
            }
    loop_signs = [orientations[e] for e in old["loop_momentum_basis"]["loop_edges"]]
    assert loop_signs == [1, 1, 1], (
        "This fixture transport preserves loop coordinates exactly"
    )
    old_routing = {}
    for line in old_dot.splitlines():
        edge = re.search(r"\[id=(\d+)\b", line)
        if edge:
            old_routing[int(edge[1])] = attribute(line, "lmb_rep")
    old_loops = {e: i for i, e in enumerate(old["loop_momentum_basis"]["loop_edges"])}
    lines = []
    for line in generated_dot.splitlines():
        if re.match(r"\s*full_num\s*=", line):
            continue  # Rendered display annotation, never a numerator input.
        if line.startswith("digraph "):
            line = f"digraph {name} {{"
        elif re.match(r"\s*overall_factor\s*=", line):
            line = (
                '    overall_factor = "' + attribute(old_dot, "overall_factor") + '";'
            )
        elif re.match(r"\s*overall_factor_evaluated\s*=", line):
            line = (
                '    overall_factor_evaluated = "'
                + attribute(old_dot, "overall_factor_evaluated")
                + '";'
            )
        edge = re.search(r"\[id=(\d+)\b", line)
        if edge:
            a = reverse_edges[int(edge[1])]
            routing = old_routing[a]
            if orientations[a] < 0:
                routing = "-(" + routing + ")"
            line = re.sub(
                r'lmb_rep="[^"]*"',
                lambda _, routing=routing: 'lmb_rep="' + routing + '"',
                line,
            )
            line = re.sub(r'\blmb_id="\d+"\s*', "", line)
            if a in old_loops:
                line = line.replace(
                    "lmb_rep=", 'lmb_id="' + str(old_loops[a]) + '" lmb_rep=', 1
                )
        lines.append(line)
    output = "\n".join(lines) + "\n"
    new = topology(output)
    for field in ("vertices", "edges", "numerator", "numerator_prefactor", "projector"):
        assert new[field] == original[field], field
    # Stronger than a parsed comparison: every local numerator and endpoint
    # slot annotation remains byte-identical in its original native statement.
    local_num = lambda text: [
        attribute(line, "num")
        for line in text.splitlines()
        if re.search(r"\bnum=", line)
    ]
    assert local_num(output) == local_num(generated_dot)
    for a, b in edges.items():
        before = old["loop_momentum_basis"]["edge_signatures"][str(a)]
        after = new["loop_momentum_basis"]["edge_signatures"][str(b)]
        assert after["loops"] == [orientations[a] * x for x in before["loops"]], (
            name,
            a,
            "loops",
        )
        assert after["external"] == [orientations[a] * x for x in before["external"]], (
            name,
            a,
            "external",
        )
    destination = OUT / (name + ".dot")
    destination.write_text(output)
    records[name] = {
        "historical_sha256": sha(old_path),
        "generated_sha256": sha(generated_path),
        "prepared_sha256": sha(destination),
        "generated": match["candidate"],
        "overall_factor": attribute(old_dot, "overall_factor_evaluated"),
        "old_to_new_vertices": vertices,
        "old_to_new_edges": edges,
        "edge_orientations": orientations,
        "loop_coordinate_signs": loop_signs,
        "transported_lmb": match["transported_lmb"],
        "checked": "Exact labelled directed topology; local numerators byte-identical; historical overall factor transported verbatim; every routing signature verified; native strict importer qualification follows",
    }
(OUT / "manifest.json").write_text(
    json.dumps(
        {
            "core_sha256": sha(core.__file__),
            "source_sha256": sha(__file__),
            "model_sha256": sha(MODEL),
            "generator_card_sha256": sha(args.generation_card),
            "records": records,
        },
        indent=2,
    )
    + "\n"
)
print(
    json.dumps(
        {
            k: {
                "generated": v["generated"],
                "factor": v["overall_factor"],
                "loops": v["loop_coordinate_signs"],
            }
            for k, v in records.items()
        },
        indent=2,
    )
)
