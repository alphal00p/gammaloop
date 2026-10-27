= Staged Python graph and exporter corrections

This bounded response to C1, C2, B1, E1, E2, E3 and E5 starts from frozen
`817b5886a50276efea05ab860d22016e1fb0c203`. The
#link("evidence/generation-integration/wp3/graph-fixes.json")[receipt]
records the exact code change, source files, native bindings, commands and results.
A new independent Astra review is required before subsequent source ports.

== Identity and callback owners

C1 uses `element_epoch` for live node, edge and half-edge identity, including map
callbacks, half-edge arguments and drawing selectors. Mutation guards and
snapshot selections retain `revision`. A mapped graph resets both counters;
traversal roots use the existing graph-bound node resolver. The regression matrix
checks every map and selector role after append and orientation changes, fresh
views against integer arguments, invalidating transformations and callback errors.

B1 gives terminal selection its own `GenerationDiagramView` assertion and
requires a nonempty completed-field check count. The native regression executes
four such checks. The renderer's `DiagramView` remains a distinct registered type.

== Retained Python and build surfaces

C2 migrates retained callers to one native node forest and the staged
`Theory`/`Legs`/`Plan` interface. Payload identity, topology, exact canonical forms,
symmetry and the 210, 278 and 845 generation oracles are retained. The tests also
state the native terminal-rooted cut, infeasible-request and DOT escaping contracts.

E1 and C2 use `canonical_stub` for the supported Community package, host and
reference snapshots. Metadata names the actual generation views and imports its
numeric alias dependencies. The native specification classes and typed factory,
codec, drawing and callback signatures remain covered. Three missing staged
display docstrings were completed to pass the existing exporter audit.

E2 keeps one module-scope import of each required exporter-test PyO3 item. E5
restores `IntoPyObject` and removes the observed unused imports in the GammaLoop
Python boundary. E3 declares the supported `native, python_stubgen` feature pair.

== Validation boundary

The receipt includes the negative frozen-binding matrix, rebuilt-binding tests,
isolated callback controls, exact exporter feature and signature checks, and
scoped lint/format results. Existing scope fingerprints were audited without
refreshing unrelated hashes: the frozen candidate's 7 inherited and 20 staged
scope/anchor mismatches remain assigned to the planned documentation rescope.
Final distribution, WASM, browser/manual rendering, full documentation, compiled
Nix source closure and combined physics/integration certification remain with
their respective owners.

== Canonicalization key follow-up

Astra round 2 closed C1, B1, E1, E2, E3 and E5, leaving C2-R2: callback
metadata still required expression labels despite native hashable-key support.
The bounded code change `1a68c9a8a5d16f2a195aa5808173764549d06e01` above
`891d7d569a408fb6e7056254802a45173ef2cf26` makes canonicalization callbacks
return the existing `_Key` contract. Rendering expression-label aliases retain
their own contract. The canonical exporter regenerated all three graph stubs;
product prose now describes hashable Python keys and Symbolica expressions.

Six focused native tests pass with zero skips, including byte, tuple and frozenset
callback keys across canonicalization, isomorphism and edge canonicalization.
They preserve payload identity and symmetry checks. The generated-type probe
connects all three API annotations to `_Key`; two exporter tests, canonical
owner checks, scoped Clippy and formatting pass. The registry title matches and
the scope audit retains 27 mismatches with no new drift. Exact source identities,
commands, logs, host SHA and retained round-2 evidence are in the receipt's
`canonicalization_key_followup_r3` record. Earlier broad controls are retained,
not rerun or newly claimed. No runtime acceptance or physics expectation changed.
A NEW independent Astra review remains required before integration.
