"""Audit a quark --symbol-indices run, using the arguments accepted by ladder.py."""

import json
import runpy
from collections import Counter
from itertools import product
from pathlib import Path
from random import Random

from symbolica import E, S

# Reuse the measured production input and pipeline; the oracle below is independent.
run = runpy.run_path(str(Path(__file__).with_name("ladder.py")))
assert run["args"].case == "quark" and run["args"].symbol_indices
result = run["result"]
aliases = result.aliases
assert aliases, "run the audit with the ordinary extension"
handles = {h.to_expression(): i for i, (h, _) in enumerate(aliases)}
bodies = [b.to_expression() for _, b in aliases]
parents = [[] for _ in aliases]
children = []
for i, body in enumerate([*bodies, result.root.to_expression()]):
    dependencies = [handles[x] for x in body.get_all_indeterminates() if x in handles]
    children.append(dependencies)
    for child in dependencies:
        parents[child].append(i)


def short(expression):
    for handle in expression.get_all_indeterminates():
        if handle in handles:
            expression = expression.replace(
                handle, S(f"walkthrough::A{handles[handle]}")
            )
    return str(expression)


leaf = max(
    (i for i in range(len(aliases)) if not children[i]),
    key=lambda i: len(parents[i]),
)
examples = sorted(parents[leaf], key=lambda i: bodies[i].get_byte_size())[:3]
selected = sorted({leaf, *examples, *(c for p in examples for c in children[p])})
sharing = {
    "definitions": len(aliases),
    "parent_count_histogram": dict(sorted(Counter(map(len, parents)).items())),
    "shared_leaf": leaf,
    "examples": examples,
    "definitions_shown": [
        {
            "id": i,
            "body": short(bodies[i]),
            "parents": parents[i],
            "handle_bytes": aliases[i][0].to_expression().get_byte_size(),
            "body_bytes": bodies[i].get_byte_size(),
        }
        for i in selected
    ],
}

# Multiply in the 16-component Clifford basis using exact integers and (+---).
# No trace recurrence, tensor admission, or Symbolica identity kernel is used.
metric = (1, -1, -1, -1)


def multiply_vector(blades, vector):
    result = [0] * 16
    for mask, coefficient in enumerate(blades):
        for mu, component in enumerate(vector):
            sign = (-1) ** (mask >> (mu + 1)).bit_count()
            if mask & (1 << mu):
                sign *= metric[mu]
            result[mask ^ (1 << mu)] += coefficient * component * sign
    return result


raw = run["raw"]
trace = next(
    x for x in raw.get_all_indeterminates() if x.get_name().split("::")[-1] == "trace"
)
word = []
for gamma in list(trace)[1]:
    argument = list(gamma)[2]
    if argument.get_name().split("::")[-1] == "mink":
        word.append(("index", str(list(argument)[1])))
    else:
        assert argument.get_name().split("::")[-1] == "Q"
        word.append(("vector", int(str(next(iter(argument))))))
dummy_names = list(dict.fromkeys(value for kind, value in word if kind == "index"))
assert len(word) == 16 and len(dummy_names) == 4
assert all(word.count(("index", name)) == 2 for name in dummy_names)
# Strip the common scalar prefactor solely for the component comparison.
# The measured reduction and alias/no-alias comparison retain it.
polynomial = (run["scalar"] / (raw / trace)).expand()
dots = [
    x
    for x in polynomial.get_all_indeterminates()
    if x.get_name().split("::")[-1] == "dot"
]
checks = []
for seed in (17, 29, 43):
    rng = Random(seed)
    momenta = {
        i: tuple(rng.randint(-3, 3) for _ in range(4))
        for i in (2, 3, 4, 5, 7, 8, 9, 10)
    }
    total = 0
    for assignment in product(range(4), repeat=4):
        indices = dict(zip(dummy_names, assignment))
        blades = [1] + [0] * 15
        weight = 1
        for mu in assignment:
            weight *= metric[mu]
        for kind, value in word:
            vector = (
                momenta[value]
                if kind == "vector"
                else tuple(int(mu == indices[value]) for mu in range(4))
            )
            blades = multiply_vector(blades, vector)
        # Only the identity component has nonzero ordinary Dirac trace.
        total += weight * 4 * blades[0]
    actual = polynomial
    for dot in dots:
        left, right = (momenta[int(str(next(iter(q))))] for q in dot)
        value = sum(g * a * b for g, a, b in zip(metric, left, right))
        actual = actual.replace(dot, E(str(value)))
    actual = actual.replace(run["dimension"], E("4"))
    assert actual == E(str(total)), (seed, actual, total)
    checks.append(
        {"seed": seed, "momenta": momenta, "trace": total, "exact_equal": True}
    )
    print(f"Independent component check {seed}: {total}", flush=True)

run["args"].save.with_suffix(".audit.json").write_text(
    json.dumps({"sharing": sharing, "word": word, "component_checks": checks}, indent=2)
    + "\n"
)
