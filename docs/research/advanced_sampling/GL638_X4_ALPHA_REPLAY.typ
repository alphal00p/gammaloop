= GL638 alpha-projection replay
<gl638-alpha-projection-replay>
The alpha implementation at `9d2bdf48367cb6bddfe7039c9331f73d1db39dae` preserves the tested GL638 physical results against frozen source `f2f64fb17dd5dab588eca298dc61f1100fa4962d`. Three raw-point evaluations use the complete #strong[936-orientation sum];, all six cuts, unchanged threshold metadata, 3D local UV and integrated UV. This is a bounded check of the physical projection; it makes no accuracy, variance or bounded-weight improvement claim.

#figure(
  align(center)[#table(
    columns: 6,
    align: (auto,auto,right,right,right,right,),
    table.header([Case], [Old → new precision], [New Re], [New Im], [Relative total change], [Largest relative event change],),
    table.hline(),
    [Hard H/Z, ordinary stack], [Double → Double], [−1.1123062424924e−32], [2.5327509671697e−32], [6.71468e−10], [7.22818e−10],
    [Same point, forced Arb], [Arb → Arb], [−1.1123062428592e−32], [2.5327509657920e−32], [1.75542e−296], [3.81138e−296],
    [Double-soft control, CTs on], [Arb → Arb], [−2.2749964951721e−32], [−4.1401765509497e−31], [9.65202e−289], [9.21304e−304],
  )]
  , kind: table
  )

Values are raw momentum-space integrand components in the unchanged internal convention, without a sampling Jacobian or Monte Carlo weight. Relative change means `|new−old|/max(|new|,|old|)` using the complex norm. The #link("GL638_X4_ALPHA_REPLAY.json")[machine artifact] retains the full native decimal components, every event identity and weight, original binary64 point tokens and precision metadata. All six event identities match in every case; all totals/events are finite and evaluation metadata report no NaN.

The hard point is archived `hz_03`, at `R=0.0002 GeV` on the HplusZminus ray defined in #link("GL638_X2_PHYSICAL_PILOT.typ")[the physical pilot];. The other point is the same archived double-soft `soft22` control, with threshold CTs enabled. Inputs are promoted from the original binary64 values by the existing precise source owner. Events remain enabled throughout.

Double precision\'s discrepancy against its own forced-Arb result increases slightly, from `2.44996e−10` to `5.15401e−10`. Both are below the requested `1e−6` tolerance, and the selected precision stays Double. The agreement therefore supports unchanged physics at these points; it does not show improved numerical accuracy.

== What this does and does not test
<what-this-does-and-does-not-test>
Live inspection confirms that every physical cut group has `max_occurence=1`: group IDs 0--5 contain CutIds `[3]`, `[4]`, `[1]`, `[0]`, `[5]`, `[2]`, respectively. The existing generator requests `max_occurence−1` LU derivatives, so each group\'s IFT backend is `None`. Consequently #strong[this GL638 state does not execute the corrected higher-order IFT];. Independent generated raised-cut tests establish that correction; these replays test the new alpha representation at derivative order zero. Neither the ordinary nor the forced-Arb comparison warrants a claim about untested points.

== Provenance and reproduction
<provenance-and-reproduction>
The diagnostic executable is SHA256 `e92c5ae081ae93334674ce3558158af8a73d124ec6981a7816e331be1d26aac9`. It was linked after the `9d2bdf483` core/native and API/differential gates from the exact Cargo API fingerprint `gammaloop-api-a196f440a26dc6bc` and core fingerprint `gammalooprs-f8d2e69e9c98f363`. Their library hashes, compiler argv, dynamic-library hashes, source and effective card text are embedded in the JSON. The executable was frozen before subsequent development resumed; the optimized X2 checkout and binary were untouched.

The existing precise API driver loaded the original checkpoint once, evaluated the three cases, and retained native totals/events as strings. The call took #strong[536.764 seconds];, dominated by loading the full state; the three evaluation calls took #strong[13.150 seconds];. Exit status was zero. All #strong[35 state-file hashes] still match the recorded X2 baseline, with no additional files. Original generation-card, graph and model hashes are also retained.

The JSON embeds the small driver, fingerprint-selecting linker, run/comparison scripts, input manifest and runtime card. To reproduce, extract those files into a scratch directory, point the manifest at the matching saved checkpoint, and link against the gated `9d2bdf483` build using the recorded API fingerprint. Relocate the manifest\'s card, executable, request and output paths together. `run.py` performs the one-load batch through `evaluate_momentum_configuration_precise`; `compare.py` uses decimal arithmetic on native components. The checkpoint\'s multi-gigabyte state itself is identified by hashes rather than copied into this report.
