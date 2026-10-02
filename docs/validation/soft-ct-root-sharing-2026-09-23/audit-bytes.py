"""Read-only Symbolica 3.0.0 packed-tree audit; no symbol import or algebra.

Length decoding follows atom/representation.rs::ListIterator::next and
atom/coefficient.rs::skip_rational. All container lengths are checked.
"""

from collections import Counter
from pathlib import Path
import hashlib, json

ROOT = Path("/tmp/soft-ct-root-size-2026-09-23")
SIZES = {0: 0, 1: 1, 2: 2, 3: 4, 4: 8}
NAMES = {1: "Num", 2: "Var", 3: "Fun", 4: "Mul", 5: "Add", 6: "Pow"}


def fraction(data, pos):
    tag = data[pos]
    pos += 1
    ns = SIZES[tag & 15]
    ds = SIZES[(tag & 112) >> 4]
    num = int.from_bytes(data[pos : pos + ns], "little")
    pos += ns
    den = int.from_bytes(data[pos : pos + ds], "little") if ds else 1
    return num, den, pos + ds


def skip_rational(data, pos):
    tag = data[pos] & 127
    num = tag & 15
    den = (tag & 112) >> 4
    if num in [1, 2, 3, 4]:
        return pos + 1 + SIZES[num] + SIZES[den]
    pos += 1
    if num in [10, 15]:
        return skip_rational(data, skip_rational(data, pos))
    if num == 7:
        n, d, pos = fraction(data, pos)
        return pos + n + d
    if num == 8:
        return pos + 4 + int.from_bytes(data[pos : pos + 4], "little")
    if num == 5:
        tag = data[pos]
        return pos + 1 + SIZES[tag & 15] + SIZES[(tag & 112) >> 4]
    if num == 9:
        for _ in range(2):
            pos += 8 + int.from_bytes(data[pos : pos + 8], "little")
        return pos
    assert num in [12, 13], num
    return pos


def end_of(data, pos):
    pending = 1
    while pending:
        kind = data[pos] & 7
        pos += 1
        pending -= 1
        if kind in [1, 2]:
            pos = skip_rational(data, pos)
        elif kind in [3, 4]:
            pos += 4 + int.from_bytes(data[pos : pos + 4], "little")
        elif kind == 5:
            _, size, pos = fraction(data, pos)
            pos += size
        elif kind == 6:
            pending += 2
        else:
            raise ValueError(f"unknown kind {kind}")
    return pos


def children(data):
    kind = data[0] & 7
    if kind == 5:
        count, size, pos = fraction(data, 1)
        assert pos + size == len(data)
    elif kind == 4:
        assert 5 + int.from_bytes(data[1:5], "little") == len(data)
        count, _, pos = fraction(data, 5)
    elif kind == 3:
        assert 5 + int.from_bytes(data[1:5], "little") == len(data)
        _, count, pos = fraction(data, 5)
    elif kind == 6:
        count, pos = 2, 1
    else:
        return []
    result = []
    for _ in range(count):
        end = end_of(data, pos)
        assert end <= len(data)
        result.append(data[pos:end])
        pos = end
    assert pos == len(data), (kind, pos, len(data), count)
    return result


results = {}
for label in ["root", "protected"]:
    raw = (ROOT / (label + ".raw")).read_bytes()
    data = memoryview(raw)
    assert end_of(data, 0) == len(data)
    terms = children(data) if data[0] & 7 == 5 else [data]
    factors = [
        factor
        for term in terms
        for factor in (children(term) if term[0] & 7 == 4 else [term])
    ]
    unique = Counter(factors)
    by_type = {}
    for kind, name in NAMES.items():
        selected = [f for f in factors if f[0] & 7 == kind]
        distinct = [f for f in unique if f[0] & 7 == kind]
        by_type[name] = {
            "occurrences": len(selected),
            "bytes": sum(map(len, selected)),
            "unique": len(distinct),
            "unique_bytes": sum(map(len, distinct)),
        }
    result = {
        "packed_bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "root_type": NAMES[data[0] & 7],
        "top_terms": len(terms),
        "term_min_bytes": min(map(len, terms)),
        "term_max_bytes": max(map(len, terms)),
        "direct_factor_occurrences": len(factors),
        "unique_direct_factors": len(unique),
        "factor_bytes": sum(map(len, factors)),
        "unique_factor_bytes": sum(map(len, unique)),
        "direct_factors_by_type": by_type,
        "largest_direct_factors": [
            {
                "kind": NAMES[f[0] & 7],
                "bytes": len(f),
                "occurrences": count,
                "children": len(children(f)),
            }
            for f, count in sorted(
                unique.items(), key=lambda pair: len(pair[0]), reverse=True
            )[:8]
        ],
    }
    results[label] = result
    print(label, json.dumps(result), flush=True)
(ROOT / "root-structure.json").write_text(json.dumps(results, indent=2) + "\n")

# Intern exact encoded subtrees for a structural sharing estimate. This is not
# an evaluated alias transformation; tensor/guard/formal scopes remain untested.
raw = (ROOT / "root.raw").read_bytes()
data = memoryview(raw)
pending = [data]
seen = set()
nodes = Counter()
edges = 0
local_bytes = 0
large_references = Counter()
while pending:
    value = pending.pop()
    if len(value) >= 1024:
        large_references[value] += 1
    if value in seen:
        continue
    seen.add(value)
    nodes[NAMES[value[0] & 7]] += 1
    parts = children(value)
    edges += len(parts)
    local_bytes += len(value) - sum(map(len, parts))
    pending.extend(parts)
summary = {
    "packed_tree_bytes": len(raw),
    "unique_nodes": len(seen),
    "unique_nodes_by_type": dict(nodes),
    "edges_between_unique_nodes": edges,
    "node_local_bytes": local_bytes,
    "illustrative_dag_bytes_with_16_byte_node_overhead_and_8_byte_references": local_bytes
    + 16 * len(seen)
    + 8 * edges,
    "interpretation": "pure structural interning estimate, not measured application memory or an implemented algebraic optimization; scope-safe evaluation not established",
    "largest_repeated_nodes": [
        {
            "kind": NAMES[value[0] & 7],
            "bytes": len(value),
            "references_from_unique_parents": count,
        }
        for value, count in sorted(
            large_references.items(), key=lambda pair: len(pair[0]), reverse=True
        )
        if count > 1
    ][:12],
}
(ROOT / "structural-sharing.json").write_text(json.dumps(summary, indent=2) + "\n")
print("sharing", json.dumps(summary), flush=True)

# A concrete lossless DAG encoding: immutable node-local bytes plus ordered
# child IDs. Reconstruct its original byte stream only into a hash, never algebra.
import struct

values = sorted(seen, key=lambda value: (len(value), bytes(value)))
indices = {value: index for index, value in enumerate(values)}
encoded = bytearray(b"SCFDAG1\0")
encoded.extend(struct.pack("<II", indices[data], len(values)))
nodes = []
for value in values:
    parts = children(value)
    prefix = bytes(value[: len(value) - sum(map(len, parts))])
    refs = [indices[part] for part in parts]
    nodes.append((prefix, refs))
    encoded.extend(struct.pack("<II", len(prefix), len(refs)))
    encoded.extend(prefix)
    encoded.extend(struct.pack("<" + "I" * len(refs), *refs))
(ROOT / "root-structure.dag").write_bytes(encoded)
hash_ = hashlib.sha256()
length = 0
pending = [indices[data]]
while pending:
    prefix, refs = nodes[pending.pop()]
    hash_.update(prefix)
    length += len(prefix)
    pending.extend(reversed(refs))
assert length == len(raw) and hash_.hexdigest() == hashlib.sha256(raw).hexdigest()
summary["lossless_dag_file_bytes"] = len(encoded)
summary["lossless_dag_roundtrip_sha256"] = hash_.hexdigest()
summary["lossless_dag_roundtrip_matches"] = True
summary["lossless_dag_format"] = (
    "8-byte magic, u32 root ID, u32 node count; each node: u32 prefix length, u32 child count, prefix bytes, ordered u32 child IDs"
)
(ROOT / "structural-sharing.json").write_text(json.dumps(summary, indent=2) + "\n")
print(
    "lossless_dag",
    len(encoded),
    "original_bytes",
    length,
    "roundtrip_sha256",
    hash_.hexdigest(),
    flush=True,
)
