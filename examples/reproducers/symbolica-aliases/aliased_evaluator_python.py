"""Python side of the aliased-evaluator reproducer.

`Expression.evaluator` and `Expression.evaluator_multiple` accept only resolved
expressions; there is no alias input and no aliased expression type. A DAG
result therefore has to be resolved into a tree before Python can evaluate it,
and a chain whose definitions use their predecessor twice doubles per level.

Run with an interpreter that has the Symbolica community build installed:
`python aliased_evaluator_python.py [depth]` (default 18).
"""

import inspect
import sys
import time

import symbolica
from symbolica import Expression, S

depth = int(sys.argv[1]) if len(sys.argv) > 1 else 18
x = S("x")

signature = inspect.signature(Expression.evaluator)
print("symbolica", getattr(symbolica, "__version__", "?"))
print("Expression.evaluator parameters:", ", ".join(signature.parameters))
print("alias input available:", any("alias" in name for name in signature.parameters))
print("aliased expression type available:", any("lias" in name for name in dir(symbolica)))

# The only route: resolve the chain into one tree, then build the evaluator.
start = time.perf_counter()
w = x
for _ in range(depth):
    w = w / 2 + w**2 / 4  # each level uses the previous one twice
resolve = time.perf_counter() - start
print(f"resolved tree for depth {depth}: {resolve:.3f} s, {w.get_byte_size() / 1e6:.1f} MB")

start = time.perf_counter()
evaluator = w.evaluator([x])
build = time.perf_counter() - start
print(f"evaluator from the tree: built in {build:.3f} s")
print("value at x = 0.5:", evaluator.evaluate([[0.5]]))
print()
print("ask: an AliasedExpression (root + (alias, body) pairs) with .evaluator(...) in the")
print("     Python bindings, mirroring AliasedAtom::evaluator_multiple in Rust.")
