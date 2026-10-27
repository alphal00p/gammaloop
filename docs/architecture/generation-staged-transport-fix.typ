= Generation staged scalar transport fix

== Review boundary

This bounded response addresses independent Astra scope A finding A1 (P1)
against staged candidate `817b5886a50276efea05ab860d22016e1fb0c203`, whose
approved integration base is `438b61500da4b32841c64211103e5da6fc4c4a72`.
It remains a candidate for a new independent review. Integration, publication,
and source19 continuation require the coordinator's remaining review gates.

== Retained history belongs to transport

Generated finalization must transport every admitted tensor through
`NumeratorInterface::translate`. Its existing metadata transport owns free
ports, contracted-index history and checked proof state. A scalar literal zero
can retain contracted indices; scalar rank and absence of visible index heads
cannot certify that its frame is final. The shortcut which inferred that
certificate from the expression is removed. Imported diagrams retain their
existing frame. The plain overall-factor atom still uses its visible-head
check because it carries no admitted tensor metadata.

The first divergence was this finalization shortcut: admission and checked
zero rewriting already retained the history correctly. The fix changes no
admission, conservation reduction, dynamic-rule or shared-template owner, and
does not expand a graph numerator.

== Regression and evidence

The regression constructs an admitted scalar contraction, then obtains a
literal zero through the existing checked-expression rewrite. It verifies that
the retained source index exists before building a generated diagram. An
internal edge inserted before the external legs moves provisional edge and
half-edge zero to their nonidentity final owners. This also permits a genuine
internal-edge fragment without violating the unit external-edge rule. The test
covers vertex and edge fragments, explicit aggregate, prefactor and projector through their
shared finalization closure, for both owner kinds. It compares the transported
history with explicit transport, rejects a subsequent third index occurrence,
and accepts multiplication by an unrelated vector while retaining its free
port on the zero.

The accompanying #link("evidence/generation-staged-transport-fix/receipt.json")[JSON receipt]
records exact code identities, commands,
fail-before/pass-after evidence, focused controls and retained artifacts.
The canonical CI Cargo profile is `ci-optim`; the test build uses the pinned
Nix environment, configured per-workspace shared Cargo cache, two build jobs
and two test threads. Existing licenses are preserved or loaded silently.

== Limits

No historical frame digest is refreshed. The three disclosed historical
frame assertions and seven inherited plus twenty staged documentation scope
or anchor mismatches remain parent gates. This note attests only its precise
finalization scope, not full documentation or physics certification. The
explicit ignored B2 control remains the original GL028/GL034 identity in
four-dimensional SU3, `N28 = (7/2) N34`; no grouping physics is changed.
Full heavy physics, the four-loop matrix, distribution certification and
other review scopes remain outside this bounded fix.
