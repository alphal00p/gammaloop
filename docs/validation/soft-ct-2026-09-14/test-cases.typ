= Per-case execution results

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<per-case-execution-results>
Times are JUnit testcase wall seconds. Source names and individual assertions are described in #link("unit-case-purposes.typ")[unit purposes] and #link("integration-case-purposes.typ")[integration purposes];. Symbolica startup aborts remain unsuccessful executions.

Status and runtime fields are deliberately unfilled. Populate from the current run's machine-readable report; preserve failures, timeouts and retry details. Source links explain the exact case and its assertions.

=== uv::approx::direct\_3d::kernel::soft\_tests
<uvapproxdirect_3dkernelsoft_tests>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Case], [Status], [Seconds], [Notes / evidence],),
    table.hline(),
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:138")[zero\_contributions\_do\_not\_publish\_projection\_paths];], [Pass], [0.059], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:184")[orientation\_term\_keeps\_external\_selectors\_and\_adds\_internal\_ones];], [Pass], [0.044], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:200")[one\_orientation\_cff\_hard\_charts\_match\_direct\_external\_jets];], [Pass], [0.095], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:353")[soft\_hard\_laurent\_projection\_matches\_external\_jets\_and\_keeps\_mass];], [Pass], [0.064], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:477")[integrated\_mass\_vertex\_keeps\_its\_soft\_weight\_and\_scalar\_logs];], [Pass], [0.064], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:603")[disconnected\_finite\_mass\_vertex\_keeps\_its\_sibling\_mass\_inert];], [Pass], [0.074], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:716")[nested\_cff\_toy\_retains\_four\_signed\_forest\_families];], [Pass], [0.233], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:874")[soft\_hard\_chart\_coscales\_free\_and\_terminal\_uv\_mass\_from\_inner\_u];], [Pass], [0.062], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:979")[soft\_laurent\_projection\_expands\_materialized\_internal\_energy];], [Pass], [0.061], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1052")[soft\_projection\_holds\_cograph\_tree\_denominators];], [Pass], [0.061], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1101")[nested\_soft\_projection\_scales\_only\_the\_component\_boundary\_after\_routing];], [Pass], [0.067], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1191")[appendix\_b1\_nested\_routes\_keep\_the\_graph\_canonical\_enclosing\_chart];], [Pass], [0.089], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1256")[nested\_route\_rejects\_a\_retained\_affine\_graph\_external\_carrier];], [Pass], [0.050], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1330")[production\_ir\_dispatch\_composes\_h\_for\_d0\_d1\_d2];], [Pass], [0.063], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1536")[one\_orientation\_ir\_dispatch\_stores\_one\_combined\_h\_atom\_per\_residue];], [Pass], [0.063], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:1720")[local\_3d\_os\_dispatch\_is\_unconditionally\_deferred];], [Pass], [0.060], [JUnit;],
  )]
  , kind: table
  )

=== uv::approx::local\_4d::tests
<uvapproxlocal_4dtests>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Case], [Status], [Seconds], [Notes / evidence],),
    table.hline(),
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:2997")[paper\_appendix\_b1\_fixture\_has\_the\_1pi\_components\_consistent\_with\_b20];], [Pass], [0.072], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3054")[appendix\_b1\_compatible\_lmb\_is\_independent\_of\_current\_full\_chart];], [Pass], [0.071], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3115")[paper\_appendix\_b1\_1pi\_component\_consistent\_with\_b20\_cograph\_pair\_cancels];], [Pass], [0.082], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3223")[top\_bubble\_boundary\_carrier\_uses\_the\_graph\_canonical\_parent\_route];], [Pass], [0.068], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3305")[top\_bubble\_actual\_fermion\_child\_h2\_removes\_both\_soft\_jet\_coefficients];], [Fail], [0.087], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3397")[nested\_soft\_carrier\_remains\_a\_parent\_loop\_coordinate];], [Pass], [0.077], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3475")[nested\_ordinary\_u\_carrier\_remains\_an\_h\_parent\_loop\_coordinate];], [Pass], [0.076], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3540")[numerator\_only\_nested\_carrier\_is\_preserved\_from\_the\_completed\_atom];], [Pass], [0.050], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3613")[degree\_zero\_ir\_reuses\_ordinary\_t\_basis\_metadata];], [Pass], [0.070], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3661")[deep\_nesting\_propagates\_the\_basis\_selected\_by\_the\_completed\_middle\_atom];], [Pass], [0.054], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3774")[nested\_expansion\_rejects\_affine\_carrier\_promotion];], [Pass], [0.053], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3908")[tilde\_taylor\_has\_degree\_d\_minus\_one\_and\_keeps\_physical\_mass];], [Pass], [0.043], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3931")[soft\_projection\_does\_not\_scale\_padded\_cograph\_externals];], [Pass], [0.062], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3998")[logarithmic\_tilde\_is\_zero\_so\_hat\_is\_t];], [Aborted: test initialization], [1.335], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4007")[local\_soft\_provenance\_records\_the\_selected\_route\_and\_branch\_sizes];], [Pass], [0.062], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4088")[production\_hat\_operator\_satisfies\_complement\_and\_refinement\_identities];], [Pass], [0.069], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4161")[production\_nested\_scheme\_matrix\_applies\_supported\_local\_operator\_pairs];], [Pass], [0.583], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4628")[full\_parent\_u\_does\_not\_project\_a\_uv\_finite\_completed\_child];], [Pass], [0.044], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4696")[full\_parent\_u\_and\_h\_treat\_a\_uv\_finite\_child\_refinement\_consistently];], [Pass], [0.045], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4803")[full\_parent\_u0\_reexpands\_s1\_child\_and\_u1\_rejects\_it];], [Pass], [0.131], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5117")[uv\_limit\_negates\_and\_marks\_the\_completed\_hat\_atom\_once];], [Pass], [0.064], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5165")[ir\_local\_operator\_is\_independent\_of\_integrated\_generation];], [Pass], [0.063], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5191")[integrated\_massless\_gauge\_soft\_jet\_vanishes\_but\_completed\_hat\_does\_not];], [Pass], [0.526], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5277")[integrated\_massive\_scalar\_soft\_jet\_retains\_nonzero\_physical\_mass];], [Pass], [0.095], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5340")[integrated\_mixed\_mass\_soft\_bubble\_matches\_massive\_tadpoles];], [Pass], [0.210], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5450")[integrated\_serial\_mixed\_mass\_sunset\_matches\_contracted\_topologies];], [Pass], [0.560], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5560")[integrated\_independent\_mixed\_mass\_tensors\_keep\_laurent\_cross\_terms];], [Pass], [0.410], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5713")[os\_dispatch\_is\_unconditionally\_deferred];], [Pass], [0.059], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6037")[paper\_eq\_4\_15\_normalized\_d0\_d1\_d2];], [Pass], [0.045], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6087")[paper\_eq\_4\_16\_quadratic\_kernel\_has\_the\_physical\_soft\_mass\_split];], [Pass], [0.046], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6149")[paper\_appendix\_b\_7\_gluon\_t2\_includes\_the\_linear\_numerator\_cross\_term];], [Pass], [0.047], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6198")[paper\_appendix\_b\_6\_b\_7\_production\_outer\_h2\_matches\_the\_routed\_rhs];], [Fail], [0.089], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6522")[paper\_appendix\_b\_10\_b\_12\_massless\_child\_and\_local\_b\_16\_nested\_branch];], [Fail], [0.069], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7145")[production\_dod2\_massive\_denominator\_has\_the\_paper\_mass\_coefficient];], [Pass], [0.061], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7180")[production\_outer\_u\_grades\_a\_free\_fermion\_mass\_but\_keeps\_muv\_fixed];], [Fail], [0.066], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7276")[nested\_outer\_u\_grades\_physical\_child\_mass\_and\_preserves\_terminal\_m];], [Pass], [0.045], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7314")[production\_outer\_u\_grades\_free\_local\_muv\_but\_preserves\_its\_terminal\_denominators];], [Pass], [0.061], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7363")[paper\_eq\_4\_19\_three\_point\_constant\_and\_linear\_terms\_use\_distinct\_masses];], [Pass], [0.045], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7399")[production\_nested\_toy\_forest\_terms\_match\_the\_explicit\_derivation];], [Pass], [3.846], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7708")[nested\_toy\_forest\_terms\_reexpand\_the\_completed\_soft\_child];], [Pass], [5.629], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7897")[rational\_shell\_groups\_mixed\_multisets\_without\_expanding\_numerators];], [Pass], [0.044], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7936")[canonical\_uv\_classes\_certify\_signed\_aliases\_and\_preserve\_raw\_sectors];], [Pass], [0.060], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8019")[canonical\_uv\_signed\_multiloop\_carriers\_preserve\_factorized\_numerators];], [Pass], [0.062], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8091")[canonical\_uv\_positive\_blocks\_and\_absent\_poles\_keep\_their\_roles];], [Pass], [0.061], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8151")[canonical\_positive\_denominator\_recovers\_one\_class\_from\_expanded\_coordinates];], [Pass], [0.062], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8199")[projection\_sector\_grouping\_preserves\_frozen\_domains\_and\_recursion];], [Pass], [0.060], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8232")[analytic\_uv\_rejects\_gamma5\_before\_simplification\_with\_subgraph\_scope];], [Pass], [0.065], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8382")[term\_projection\_preserves\_complete\_factorized\_values];], [Pass], [0.045], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8428")[nested\_uv\_rescaling\_keeps\_child\_momentum\_provenance\_immutable];], [Pass], [0.061], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8538")[uv\_taylor\_provenance\_erasure\_matches\_plain\_child\_lmb\_expansion];], [Pass], [0.052], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:8756")[gl24\_dod\_two\_q1\_quartic\_taylor\_keeps\_owner\_local\_energy\_families];], [Pass], [0.062], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9127")[dod\_one\_triangle\_keeps\_separate\_denominator\_topologies];], [Pass], [0.109], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9296")[factorized\_product\_separates\_active\_and\_completed\_sectors];], [Pass], [0.050], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9392")[affine\_uv\_rescaling\_preserves\_enclosing\_chart\_and\_owner];], [Pass], [0.063], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9590")[early\_color\_simplification\_preserves\_open\_and\_nested\_numerator\_boundaries];], [Pass], [0.064], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:9680")[local\_taylor\_retains\_dirac\_traces\_with\_or\_without\_analytic\_addbacks];], [Pass], [0.063], [JUnit;],
  )]
  , kind: table
  )

=== uv::hedge\_poset::tests
<uvhedge_posettests>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Case], [Status], [Seconds], [Notes / evidence],),
    table.hline(),
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:1714")[three\_d\_compute\_skips\_unused\_four\_d\_atoms\_when\_integration\_is\_disabled];], [Pass], [0.064], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:1785")[local\_leaf\_operations\_follow\_dependency\_frontiers];], [Pass], [0.060], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:1897")[union\_terms\_project\_factorized\_typed\_4d\_values];], [Pass], [0.064], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:1949")[union\_terms\_replay\_component\_paths\_from\_typed\_roots];], [Pass], [0.234], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2157")[disconnected\_soft\_union\_is\_the\_product\_of\_completed\_component\_atoms];], [Pass], [0.714], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2378")[heterogeneous\_terminal\_projects\_components\_before\_multiplying];], [Pass], [0.153], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2474")[triple\_tadpole];], [Pass], [0.099], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2576")[saclay];], [Pass], [0.108], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2615")[dumbells];], [Pass], [0.087], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2666")[bugblatter];], [Pass], [0.398], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2743")[mercedes];], [Pass], [0.088], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2779")[sunrise];], [Pass], [0.088], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2813")[dotted\_sunrise];], [Pass], [0.086], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2850")[dotted];], [Pass], [0.088], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:2888")[spectacles];], [Pass], [0.096], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3185")[collective\_two\_component\_region\_requires\_independently\_divergent\_factors];], [Pass], [0.065], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3381")[collective\_three\_component\_region\_keeps\_only\_its\_divergent\_join];], [Fail], [0.131], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3662")[spectacles\_typed\_4d\_local\_construction\_matches\_uv\_limit];], [Pass], [0.065], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3720")[basketball];], [Pass], [0.100], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3752")[fourloop\_b];], [Pass], [0.106], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3792")[four\_loop\_a];], [Pass], [0.103], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:3829")[triple\_double\_tadpole];], [Pass], [0.274], [JUnit;],
  )]
  , kind: table
  )

=== uv::approx::integrated::tests
<uvapproxintegratedtests>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Case], [Status], [Seconds], [Notes / evidence],),
    table.hline(),
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:1950")[projected\_dirac\_algebra\_precedes\_single\_and\_double\_pole\_expansion];], [Pass], [0.258], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2056")[projected\_dirac\_numerator\_matches\_scalar\_master\_before\_backend\_truncation];], [Pass], [0.218], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2146")[scalar\_dimensions\_are\_substituted\_inside\_functions\_but\_not\_lorentz\_slots];], [Pass], [0.043], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2162")[analytic\_uv\_spin\_algebra\_rejects\_unsupported\_residual\_contractions];], [Pass], [0.048], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2207")[residual\_tensor\_contractions\_are\_checked\_per\_analytic\_branch];], [Pass], [0.050], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2324")[integrated\_counterterm\_projects\_one\_laurent\_expansion];], [Pass], [0.044], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2357")[finite\_projection\_retains\_child\_epsilon\_terms\_until\_parent\_poles\_multiply];], [Pass], [0.044], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2395")[factorized\_product\_projects\_each\_component];], [Pass], [0.045], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2438")[vakint\_dot\_conversion\_keeps\_loop\_momentum\_tagged\_until\_to\_dots];], [Pass], [0.045], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2472")[underdetermined\_vakint\_momentum\_solve\_tracks\_free\_variables];], [Pass], [0.044], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2495")[underdetermined\_vakint\_momentum\_solve\_projects\_topology];], [Pass], [0.044], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2537")[empty\_vakint\_momentum\_solve\_treats\_variables\_as\_free];], [Pass], [0.043], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2550")[vakint\_conversion\_preserves\_squared\_masses];], [Pass], [0.062], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2625")[vakint\_conversion\_rejects\_a\_lost\_active\_loop];], [Pass], [0.066], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2682")[affine\_vakint\_routing\_keeps\_external\_momenta\_fixed];], [Pass], [0.065], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2796")[nested\_active\_sunset\_keeps\_its\_two\_loop\_incidence];], [Pass], [0.063], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/integrated.rs:2905")[nested\_vacuum\_denominators\_preserve\_incidence\_and\_factorized\_numerators];], [Pass], [0.062], [JUnit;],
  )]
  , kind: table
  )

=== uv::orchestrator::tests
<uvorchestratortests>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Case], [Status], [Seconds], [Notes / evidence],),
    table.hline(),
    [#link("../../../crates/gammalooprs/src/uv/orchestrator.rs:427")[compare\_canonicalizes\_contracted\_uv\_indices];], [Pass], [0.051], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/orchestrator.rs:462")[integrated\_soft\_scheme\_is\_allowed\_before\_forest\_generation];], [Pass], [0.053], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/orchestrator.rs:485")[incompatible\_selected\_spinney\_route\_is\_reported\_with\_context];], [Pass], [0.056], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/orchestrator.rs:540")[integrated\_on\_shell\_scheme\_is\_rejected\_before\_forest\_generation];], [Pass], [0.055], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/orchestrator.rs:567")[soft\_child\_rejects\_positive\_degree\_muv\_parent];], [Pass], [0.061], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/orchestrator.rs:655")[soft\_and\_pole\_part\_components\_cannot\_share\_a\_wood];], [Pass], [0.075], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/orchestrator.rs:698")[disjoint\_soft\_and\_pole\_part\_components\_cannot\_share\_a\_wood];], [Pass], [0.067], [JUnit;],
  )]
  , kind: table
  )

=== uv::approx::projected\_4d::tests
<uvapproxprojected_4dtests>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Case], [Status], [Seconds], [Notes / evidence],),
    table.hline(),
    [#link("../../../crates/gammalooprs/src/uv/approx/projected_4d.rs:594")[canonical\_source\_and\_template\_preparation\_share\_retention\_invariant\_budget];], [Pass], [0.082], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/projected_4d.rs:682")[nested\_banana\_quotient\_powered\_component\_has\_the\_analytic\_one\_energy\_sign];], [Pass], [0.074], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/projected_4d.rs:888")[typed\_taylor\_wave\_batches\_genuine\_owner\_relabelled\_terms];], [Pass], [0.055], [JUnit;],
    [#link("../../../crates/gammalooprs/src/uv/approx/projected_4d.rs:1020")[typed\_taylor\_next\_component\_reuses\_one\_topology\_for\_prior\_residue\_states];], [Pass], [0.064], [JUnit;],
  )]
  , kind: table
  )

=== Related tensor and backend regressions
<related-tensor-and-backend-regressions>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Case], [Status], [Seconds], [Notes / evidence],),
    table.hline(),
    [#link("../../../crates/spenso/src/network/tests.rs:276")[network::tests::scalar\_sums\_accept\_closed\_lazy\_tensors\_in\_either\_order];], [Blocked: Symbolica license], [0.774], [JUnit;],
    [#link("../../../crates/spenso/src/network/tests.rs:364")[network::tests::tensor\_powers\_preserve\_exponents\_for\_every\_leaf\_kind];], [Blocked: Symbolica license], [0.802], [JUnit;],
    [#link("../../../crates/spenso/src/network/parsing/test.rs:180")[network::parsing::test::trace\_projectors\_preserve\_exposed\_slots\_and\_close\_with\_metrics];], [Blocked: Symbolica license], [1.025], [JUnit;],
    [#link("../../../crates/idenso/src/tensor/tests/parsing.rs:798")[tensor::tests::parsing::parse\_problem];], [Blocked: Symbolica license], [1.223], [JUnit;],
    [#link("../../../crates/vakint/tests/input_matching_tests.rs:142")[partly\_massless\_sunsets\_preserve\_mass\_labels\_and\_exclude\_alphaloop];], [Blocked: Symbolica license], [1.033], [JUnit;],
    [#link("../../../crates/vakint/tests/integral_evaluation_analytic_tests.rs:261")[partly\_massless\_sunsets\_match\_independent\_gamma\_integrals];], [Blocked: Symbolica license], [1.152], [JUnit;],
    [#link("../../../crates/vakint/tests/tensor_reduction_tests.rs:102")[public\_parser\_initializes\_dot\_attributes\_before\_first\_input];], [Blocked: Symbolica license], [1.046], [JUnit;],
  )]
  , kind: table
  )

=== Mandatory slow CLI acceptance
<mandatory-slow-cli-acceptance>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Case], [Status], [Seconds], [Profiles / evidence],),
    table.hline(),
    [#link("../../../tests/tests/uv.rs:5526")[slow::paper\_figure\_b1\_massless\_bubble\_uses\_a\_complete\_soft\_wood];], [Fail: generation panic], [6326.706], [JUnit;],
    [#link("../../../tests/tests/uv.rs:5555")[slow::paper\_figure\_b1\_double\_triangle\_is\_a\_non\_vacuous\_uv\_fixture];], [Pass], [58.997], [JUnit;],
    [#link("../../../tests/tests/uv.rs:5620")[slow::local\_nested\_soft\_fixture\_generates];], [Pass], [19.888], [JUnit;],
    [#link("../../../tests/tests/uv.rs:5746")[slow::paper\_appendix\_b1\_ir\_specialization\_generates\_and\_profiles];], [Fail], [72.271], [JUnit;],
    [#link("../../../tests/tests/uv.rs:6164")[slow::gamma\_star\_ddbar\_top\_bubble\_child\_only\_has\_two\_power\_soft\_improvement];], [Pass], [57.020], [JUnit;],
    [#link("../../../tests/tests/uv.rs:6299")[slow::gamma\_star\_ddbar\_top\_bubble\_has\_two\_power\_soft\_improvement];], [Pass], [198.364], [JUnit;],
  )]
  , kind: table
  )

=== Regular integration cases
<regular-integration-cases>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Case], [Result], [Seconds], [Evidence],),
    table.hline(),
    [`consecutive_soft_components_have_exact_forest_and_cli_profile`], [Pass], [8.409], [JUnit;],
    [`dgse_local_ir_generation_orientation_filter_profiles_one_orientation`], [Pass], [3.761], [JUnit;],
    [`dgse_local_ir_profiles_are_non_vacuous`], [Pass], [35.451], [JUnit;],
    [`genuine_vacuum_tadpole_matches_all_local_uv_routes_and_vakint_inputs`], [Pass], [4.058], [JUnit;],
    [`local_ir_and_pole_part_wood_is_rejected_by_run_generate`], [Pass], [0.722], [JUnit;],
    [`local_ir_child_under_positive_degree_muv_parent_is_rejected_by_run_generate`], [Pass], [0.672], [JUnit;],
    [`local_ir_disconnected_completed_components_form_nonzero_union`], [Pass], [11.977], [JUnit;],
    [`local_ir_integrated_generation_succeeds_across_local_uv_routes`], [Pass], [35.984], [JUnit;],
    [`local_os_integrated_generation_is_rejected_by_run_generate`], [Pass], [0.730], [JUnit;],
    [`massive_top_self_energy_local_os_dispatch_is_deferred`], [Pass], [0.636], [JUnit;],
    [`massless_gluon_self_energy_local_ir_profile_is_non_vacuous`], [Pass], [4.216], [JUnit;],
    [`massless_quark_self_energy_local_ir_profile_is_non_vacuous`], [Pass], [3.927], [JUnit;],
    [`paper_appendix_b1_os_child_reaches_deferred_dispatch`], [Pass], [0.632], [JUnit;],
    [`scalar_amplitudes_match_across_local_uv_routes`], [Pass], [5.295], [JUnit;],
    [`scalar_spectacles_integrated_matches_across_local_uv_routes`], [Pass], [5.957], [JUnit;],
    [`soft_cff_state_survives_a_fresh_process_reload`], [Pass], [3.565], [JUnit;],
    [`sunrise_pole_part_matches_muv_inspect`], [Pass], [4.703], [JUnit;],
  )]
  , kind: table
  )
