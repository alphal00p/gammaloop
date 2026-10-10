= Soft-counterterm unit-test inventory

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<soft-counterterm-unit-test-inventory>
Read-only source inventory, 2026-09-14. No builds or tests were run by this inventory task. Repository root: /common/dev/gammaloop/lcnbr.

== Suggested execution filters
<suggested-execution-filters>
Core unit suite: 118 runnable source-level tests (120 declared, with two existing failing-module diagnostics excluded). Confirm the final count using nextest list because compiled feature selection is authoritative.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --no-fail-fast -p gammalooprs --lib -E 'package(gammalooprs) and test(/^uv::(approx::(direct_3d::kernel::soft_tests|local_4d::tests|integrated::tests)|orchestrator::tests|hedge_poset::tests)::/) and not test(/::failing::/)'
```

Use the environment/build settings selected by the main validation task. No need to rebuild separately per subgroup. Keep serial execution if the Symbolica license permits only one process. The test\_gammaloop profile records runtime without terminating slow tests; runner invocation should preserve external bounds suitable for this task. CI's soft\_ct\_acceptance profile is reserved for the six mandatory slow integration tests and does not select these unit tests by default.

Do not use profile local\_test for an acceptance claim: its configuration uses on-timeout = \"pass\".

Related tensor regression filter:

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --no-fail-fast -p spenso -p idenso --lib -E 'test(/^network::tests::(scalar_sums_accept_closed_lazy_tensors_in_either_order|tensor_powers_preserve_exponents_for_every_leaf_kind)$/) or test(/^network::parsing::test::trace_projectors_preserve_exposed_slots_and_close_with_metrics$/) or test(/^tensor::tests::parsing::parse_problem$/)'
```

Related Vakint regressions (backend dependencies apply):

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --no-fail-fast -p vakint -E 'test(=partly_massless_sunsets_preserve_mass_labels_and_exclude_alphaloop) or test(=partly_massless_sunsets_match_independent_gamma_integrals) or test(=public_parser_initializes_dot_attributes_before_first_input)'
```

== What the groups establish
<what-the-groups-establish>
#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,auto,),
    table.header([Group], [Declared tests], [Scope],),
    table.hline(),
    [Direct 3D soft operator], [16], [d=0/1/2 H=U+S-US; factorized S+U(X-S) equivalence; one signed H atom per selected residue; canonical routes; boundary-only rescaling; fixed cograph tree poles; physical and terminal UV masses; integrated finite-vertex weight and inert disconnected sibling; four signed nested families; OS rejection.],
    [Local 4D operator and reconstruction boundary], [56], [Soft jets and H identities; logarithmic H0=U0; physical-mass retention; top-bubble cancellation of constant and linear soft coefficients; canonical and affine-route invariants; nonlocal/numerator-only carriers; nested U/H matrix and three hard limits; exact toy families; Appendix B formulas; integrated massless/massive/mixed-mass jets; immutable momentum provenance; factorized canonical denominators; color and Dirac closure.],
    [Integrated boundary], [17], [Dirac algebra before epsilon truncation; unresolved tensor rejection; finite/pole projection; positive child epsilon coefficients retained for parent poles; component factorization; squared mass preservation; momentum solve and lost-loop rejection; affine external shifts; nested sunset incidence and factorized numerators.],
    [Scheme orchestration], [7], [IR integration allowed; unsupported OS integration rejected; contextual incompatible route errors; positive-degree MUV over IR child rejected; nested or disjoint IR/PolePart woods rejected; canonical diagnostic comparisons.],
    [HedgePoset forest graph], [24 (22 selected)], [Completed disconnected products; typed 4D sectors and direct-3D replay paths; dependency frontiers; independent UV divergence before admitting a disconnected union; local projected data retained when required; representative forest topology snapshots. Two existing failing-module tests are excluded.],
    [Tensor regressions], [4 selected], [Closed lazy scalar tensors can enter scalar sums in either order; tensor powers preserve original exponent and odd open-vector rank; projected traces retain open slots until metric closure; current GammaLoop Q3/OSE parsing fixture remains accepted.],
    [Vakint backend regressions], [3 selected], [Partly massless sunset mass labels and backend applicability; comparison to independent gamma-function integrals; public-parser dot attributes initialized before first input.],
  )]
  , kind: table
  )

Most unit tests compare exact symbolic expressions or structural invariants. Runtime from them is useful for identifying costly validation cases, but does not by itself measure production evaluation throughput or speedup relative to the pre-rebase baseline.

== Exclusions and expected rejection tests
<exclusions-and-expected-rejection-tests>
- No \#\[ignore\] annotations occur in the five core inventoried source files.
- Existing HedgePoset diagnostics uv::hedge\_poset::tests::failing::lobsided\_double\_dumbell and uv::hedge\_poset::tests::failing::double\_double\_dumbell are excluded by the repository's curated profile and the explicit filter above. They are not soft-CT acceptance coverage.
- The two OS dispatch tests intentionally pass only when the documented deferred-dispatch panic occurs; they are not evidence of OS functionality.
- Scheme-error tests intentionally reject unsupported mixed woods. H/U with a positive-degree ordinary parent is unsupported; logarithmic ordinary parents remain covered.
- integrated.rs has a commented-out test annotation near line 2413; it is not a runnable test.
- HedgePoset bugblatter is a known slow diagnostic under contention; CI grants it up to ten minutes.
- Unit coverage does not replace the six mandatory slow CLI acceptance fixtures in profile soft\_ct\_acceptance or normal CLI fixtures. Those exercise generator, evaluator, physical hard/soft rays, baseline failures and profile JSON. Keep them separately reported.
- Existing documentation's original Phase-1 prohibition on integrated IR is superseded by the implemented shared integration boundary and newer architecture text.

== Exact source inventory
<exact-source-inventory>
=== uv::approx::direct\_3d::kernel::soft\_tests (16)
<uvapproxdirect_3dkernelsoft_tests-16>
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:138")[zero\_contributions\_do\_not\_publish\_projection\_paths]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:184")[orientation\_term\_keeps\_external\_selectors\_and\_adds\_internal\_ones]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:200")[one\_orientation\_cff\_hard\_charts\_match\_direct\_external\_jets]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:353")[soft\_hard\_laurent\_projection\_matches\_external\_jets\_and\_keeps\_mass]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:477")[integrated\_mass\_vertex\_keeps\_its\_soft\_weight\_and\_scalar\_logs]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:603")[disconnected\_finite\_mass\_vertex\_keeps\_its\_sibling\_mass\_inert]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:716")[nested\_cff\_toy\_retains\_four\_signed\_forest\_families]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:874")[soft\_hard\_chart\_coscales\_free\_and\_terminal\_uv\_mass\_from\_inner\_u]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:979")[soft\_laurent\_projection\_expands\_materialized\_internal\_energy]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1052")[soft\_projection\_holds\_cograph\_tree\_denominators]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1101")[nested\_soft\_projection\_scales\_only\_the\_component\_boundary\_after\_routing]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1191")[appendix\_b1\_nested\_routes\_keep\_the\_graph\_canonical\_enclosing\_chart]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1256")[nested\_route\_rejects\_a\_retained\_affine\_graph\_external\_carrier]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1330")[production\_ir\_dispatch\_composes\_h\_for\_d0\_d1\_d2]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1536")[one\_orientation\_ir\_dispatch\_stores\_one\_combined\_h\_atom\_per\_residue]
- #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1720")[local\_3d\_os\_dispatch\_is\_unconditionally\_deferred]

=== uv::approx::local\_4d::tests (56)
<uvapproxlocal_4dtests-56>
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:2997")[paper\_appendix\_b1\_fixture\_has\_the\_1pi\_components\_consistent\_with\_b20]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3054")[appendix\_b1\_compatible\_lmb\_is\_independent\_of\_current\_full\_chart]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3115")[paper\_appendix\_b1\_1pi\_component\_consistent\_with\_b20\_cograph\_pair\_cancels]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3223")[top\_bubble\_boundary\_carrier\_uses\_the\_graph\_canonical\_parent\_route]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3305")[top\_bubble\_actual\_fermion\_child\_h2\_removes\_both\_soft\_jet\_coefficients]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3397")[nested\_soft\_carrier\_remains\_a\_parent\_loop\_coordinate]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3475")[nested\_ordinary\_u\_carrier\_remains\_an\_h\_parent\_loop\_coordinate]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3540")[numerator\_only\_nested\_carrier\_is\_preserved\_from\_the\_completed\_atom]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3613")[degree\_zero\_ir\_reuses\_ordinary\_t\_basis\_metadata]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3661")[deep\_nesting\_propagates\_the\_basis\_selected\_by\_the\_completed\_middle\_atom]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3774")[nested\_expansion\_rejects\_affine\_carrier\_promotion]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3908")[tilde\_taylor\_has\_degree\_d\_minus\_one\_and\_keeps\_physical\_mass]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3931")[soft\_projection\_does\_not\_scale\_padded\_cograph\_externals]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3998")[logarithmic\_tilde\_is\_zero\_so\_hat\_is\_t]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4007")[local\_soft\_provenance\_records\_the\_selected\_route\_and\_branch\_sizes]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4088")[production\_hat\_operator\_satisfies\_complement\_and\_refinement\_identities]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4161")[production\_nested\_scheme\_matrix\_applies\_supported\_local\_operator\_pairs]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4628")[full\_parent\_u\_does\_not\_project\_a\_uv\_finite\_completed\_child]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4696")[full\_parent\_u\_and\_h\_treat\_a\_uv\_finite\_child\_refinement\_consistently]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4803")[full\_parent\_u0\_reexpands\_s1\_child\_and\_u1\_rejects\_it]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5117")[uv\_limit\_negates\_and\_marks\_the\_completed\_hat\_atom\_once]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5165")[ir\_local\_operator\_is\_independent\_of\_integrated\_generation]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5191")[integrated\_massless\_gauge\_soft\_jet\_vanishes\_but\_completed\_hat\_does\_not]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5277")[integrated\_massive\_scalar\_soft\_jet\_retains\_nonzero\_physical\_mass]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5340")[integrated\_mixed\_mass\_soft\_bubble\_matches\_massive\_tadpoles]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5450")[integrated\_serial\_mixed\_mass\_sunset\_matches\_contracted\_topologies]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5560")[integrated\_independent\_mixed\_mass\_tensors\_keep\_laurent\_cross\_terms]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5713")[os\_dispatch\_is\_unconditionally\_deferred]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6037")[paper\_eq\_4\_15\_normalized\_d0\_d1\_d2]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6087")[paper\_eq\_4\_16\_quadratic\_kernel\_has\_the\_physical\_soft\_mass\_split]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6149")[paper\_appendix\_b\_7\_gluon\_t2\_includes\_the\_linear\_numerator\_cross\_term]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6198")[paper\_appendix\_b\_6\_b\_7\_production\_outer\_h2\_matches\_the\_routed\_rhs]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6522")[paper\_appendix\_b\_10\_b\_12\_massless\_child\_and\_local\_b\_16\_nested\_branch]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7145")[production\_dod2\_massive\_denominator\_has\_the\_paper\_mass\_coefficient]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7180")[production\_outer\_u\_grades\_a\_free\_fermion\_mass\_but\_keeps\_muv\_fixed]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7276")[nested\_outer\_u\_grades\_physical\_child\_mass\_and\_preserves\_terminal\_m]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7314")[production\_outer\_u\_grades\_free\_local\_muv\_but\_preserves\_its\_terminal\_denominators]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7363")[paper\_eq\_4\_19\_three\_point\_constant\_and\_linear\_terms\_use\_distinct\_masses]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7399")[production\_nested\_toy\_forest\_terms\_match\_the\_explicit\_derivation]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7708")[nested\_toy\_forest\_terms\_reexpand\_the\_completed\_soft\_child]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7897")[rational\_shell\_groups\_mixed\_multisets\_without\_expanding\_numerators]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7936")[canonical\_uv\_classes\_certify\_signed\_aliases\_and\_preserve\_raw\_sectors]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8019")[canonical\_uv\_signed\_multiloop\_carriers\_preserve\_factorized\_numerators]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8091")[canonical\_uv\_positive\_blocks\_and\_absent\_poles\_keep\_their\_roles]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8151")[canonical\_positive\_denominator\_recovers\_one\_class\_from\_expanded\_coordinates]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8199")[projection\_sector\_grouping\_preserves\_frozen\_domains\_and\_recursion]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8232")[analytic\_uv\_rejects\_gamma5\_before\_simplification\_with\_subgraph\_scope]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8382")[term\_projection\_preserves\_complete\_factorized\_values]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8428")[nested\_uv\_rescaling\_keeps\_child\_momentum\_provenance\_immutable]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8538")[uv\_taylor\_provenance\_erasure\_matches\_plain\_child\_lmb\_expansion]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8756")[gl24\_dod\_two\_q1\_quartic\_taylor\_keeps\_owner\_local\_energy\_families]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9127")[dod\_one\_triangle\_keeps\_separate\_denominator\_topologies]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9296")[factorized\_product\_separates\_active\_and\_completed\_sectors]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9392")[affine\_uv\_rescaling\_preserves\_enclosing\_chart\_and\_owner]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9590")[early\_color\_simplification\_preserves\_open\_and\_nested\_numerator\_boundaries]
- #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9680")[local\_taylor\_retains\_dirac\_traces\_with\_or\_without\_analytic\_addbacks]

=== uv::hedge\_poset::tests (24)
<uvhedge_posettests-24>
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:1714")[three\_d\_compute\_skips\_unused\_four\_d\_atoms\_when\_integration\_is\_disabled]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:1785")[local\_leaf\_operations\_follow\_dependency\_frontiers]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:1897")[union\_terms\_project\_factorized\_typed\_4d\_values]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:1949")[union\_terms\_replay\_component\_paths\_from\_typed\_roots]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2157")[disconnected\_soft\_union\_is\_the\_product\_of\_completed\_component\_atoms]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2378")[heterogeneous\_terminal\_projects\_components\_before\_multiplying]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2474")[triple\_tadpole]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2576")[saclay]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2615")[dumbells]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2666")[bugblatter]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2743")[mercedes]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2779")[sunrise]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2813")[dotted\_sunrise]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2850")[dotted]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2888")[spectacles]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3185")[collective\_two\_component\_region\_requires\_independently\_divergent\_factors]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3381")[collective\_three\_component\_region\_keeps\_only\_its\_divergent\_join]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3662")[spectacles\_typed\_4d\_local\_construction\_matches\_uv\_limit]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3720")[basketball]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3752")[fourloop\_b]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3792")[four\_loop\_a]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3829")[triple\_double\_tadpole]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3864")[lobsided\_double\_dumbell]
- #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3940")[double\_double\_dumbell]

=== uv::approx::integrated::tests (17)
<uvapproxintegratedtests-17>
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:1950")[projected\_dirac\_algebra\_precedes\_single\_and\_double\_pole\_expansion]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2056")[projected\_dirac\_numerator\_matches\_scalar\_master\_before\_backend\_truncation]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2146")[scalar\_dimensions\_are\_substituted\_inside\_functions\_but\_not\_lorentz\_slots]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2162")[analytic\_uv\_spin\_algebra\_rejects\_unsupported\_residual\_contractions]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2207")[residual\_tensor\_contractions\_are\_checked\_per\_analytic\_branch]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2324")[integrated\_counterterm\_projects\_one\_laurent\_expansion]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2357")[finite\_projection\_retains\_child\_epsilon\_terms\_until\_parent\_poles\_multiply]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2395")[factorized\_product\_projects\_each\_component]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2438")[vakint\_dot\_conversion\_keeps\_loop\_momentum\_tagged\_until\_to\_dots]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2472")[underdetermined\_vakint\_momentum\_solve\_tracks\_free\_variables]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2495")[underdetermined\_vakint\_momentum\_solve\_projects\_topology]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2537")[empty\_vakint\_momentum\_solve\_treats\_variables\_as\_free]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2550")[vakint\_conversion\_preserves\_squared\_masses]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2625")[vakint\_conversion\_rejects\_a\_lost\_active\_loop]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2682")[affine\_vakint\_routing\_keeps\_external\_momenta\_fixed]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2796")[nested\_active\_sunset\_keeps\_its\_two\_loop\_incidence]
- #link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2905")[nested\_vacuum\_denominators\_preserve\_incidence\_and\_factorized\_numerators]

=== uv::orchestrator::tests (7)
<uvorchestratortests-7>
- #link("../../../crates/gammalooprs/src/uv/orchestrator.rs:427")[compare\_canonicalizes\_contracted\_uv\_indices]
- #link("../../../crates/gammalooprs/src/uv/orchestrator.rs:462")[integrated\_soft\_scheme\_is\_allowed\_before\_forest\_generation]
- #link("../../../crates/gammalooprs/src/uv/orchestrator.rs:485")[incompatible\_selected\_spinney\_route\_is\_reported\_with\_context]
- #link("../../../crates/gammalooprs/src/uv/orchestrator.rs:540")[integrated\_on\_shell\_scheme\_is\_rejected\_before\_forest\_generation]
- #link("../../../crates/gammalooprs/src/uv/orchestrator.rs:567")[soft\_child\_rejects\_positive\_degree\_muv\_parent]
- #link("../../../crates/gammalooprs/src/uv/orchestrator.rs:655")[soft\_and\_pole\_part\_components\_cannot\_share\_a\_wood]
- #link("../../../crates/gammalooprs/src/uv/orchestrator.rs:698")[disjoint\_soft\_and\_pole\_part\_components\_cannot\_share\_a\_wood]
