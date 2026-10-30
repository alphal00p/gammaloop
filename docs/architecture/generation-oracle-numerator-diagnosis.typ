= Three one-diagram QED oracle numerator diagnosis

== Scope and outcome

This bounded investigation covers only `xs_qed_eedd_L1`, `xs_qed_eemm_L1`, and
`xs_qed_eeee_1pi_L1` in `generation_oracle`. All use builtin `Model::qed`, their
original requests with filters disabled, `Classes::Auto`, `Naming::prefix("D")`,
and one worker. No expectation was recorded or changed. This is not an all-81
certificate, a pool-8/pool-32 result, an integration, or a broad physical
correctness claim.

The untouched historical source `2b1ef4f0` with locked Symbolica `c78cffc2`
passes all three original assertions. The untouched current source `b9f78892`
with locked Symbolica `1a1829b1` fails all three at only the numerator field,
reproducing `f5a21f1a84870bc7 -> d3ad2f3a0802b6ab`,
`8b7b40bba961420d -> 414e2a7a954beae8`, and
`0db0b7cc52558270 -> 1fa84d9ad4fdf0b0`. Both producers use Rust 1.98.1,
the default native feature set, `dev-optim`, their unchanged lock files and
their own Nix development shells. Thus the old pass is a real old-kernel run,
not historical data replayed through the current kernel.

== Complete pipeline comparison

The compared stages are the builtin model and resolved request, topology and
interaction enumeration, automatic classes, local numerator construction and
typed admission, Named-D build, grouping, `Census::Summary::of` serialization,
and the FNV-1a multiset digest. The model and options captures agree exactly.
Each route has the same one graph, canonical key, ID and name; complete physical
incidence, cuts and signed loop routing; overall and source factors; and group
master, member and ratio. Counts, rational weights, symmetry histogram and all
non-numerator digests therefore agree for a concrete reason, not only by their
final hashes.

The last shared boundary is the complete physical graph plus every local
owner's labelled representation, dimension and variance port multiset. The
first unequal boundary is local numerator representation/export. Historical
payloads print Q and edge constructors with presentation tags; current payloads
use the corresponding notation-owned labels. Current payloads also serialize
an explicit logical interface for selected tensors whose historical form was
inferable after reversible cooking. Consequently the first JSON string supplied
to the numerator FNV is already unequal at the aggregate numerator; hashing and
multiset ordering are downstream and excluded.

The three cases differ slightly within that common boundary. The eemm row has
no newly explicit interfaces and differs only in the audited constructor
spelling. The eedd row has four transported owners (`vertex:0`, `vertex:1`,
`edge:2`, and `edge:3`), with historical-to-current axis permutations
`[4,3,0,1,2]`, `[4,3,0,1,2]`, `[3,2,0,1]`, and `[3,2,0,1]`. The one-particle-
irreducible eeee row transports `vertex:1` and `vertex:2` by `[1,0,2]`.
The complete typed-axis totals are 36, 28 and 28; every port set and all metric
and composite flags agree.

== Factor-preserving certificate

The certificate binds only the five already audited constructor spellings for
Q, K, P, edge and vertex. K, P and vertex do not occur in every row, but keeping
the closed audited map avoids inventing a case-local normalization. Historical
live ports are decoded through that producer's actual reversible
`CookSettings::uncook` owner. Each decoded axis maps injectively to the equal
current labelled axis, retaining representation, dimension, variance, owner
and order. After applying those explicit transports and omitting only the now
redundant serialized interface fields, all three complete current dumps equal
their historical dumps exactly.

This proof keeps every local and aggregate numerator factorized. It performs no
whole-graph expansion, distribution, Dirac trace, numerical sampling or
parameter-at-one substitution. SHA-free mutation controls reject a trailing
numerator factor, a typed dimension change and a group-ratio change. The
#link("evidence/generation-oracle-numerator-diagnosis/receipt.json")[structured receipt],
#link("evidence/generation-oracle-numerator-diagnosis/validate.mjs")[validator],
manifest and complete historical/current payload and live-port captures are
retained together. The pristine/adapted producer patches and pristine ELFs and
logs remain in each task workspace under `target/oracle-numerator-diagnosis/`,
with hashes in the receipt.

== Limits and next boundary

The result identifies a representation/interface serialization mismatch, not a
physical numerator defect, for exactly these three graphs. It does not infer the
same result for the remaining 78 current failures, and it does not establish
thread-pool invariance beyond the original failure at pool 1. A source/kernel
pin-swap counterfactual was deliberately not created, so the two locked
pipelines jointly control their source and kernel histories rather than
assigning every byte difference uniquely to one of them.

The next comparison for another model family or a grouped case is one bounded
full-key and local-owner typed-axis transport followed by exact checks of every
group master, member source factor and ratio. That comparison must precede any
grouped membership/master/ratio or oracle-digest migration. This work stops for
new independent review.
