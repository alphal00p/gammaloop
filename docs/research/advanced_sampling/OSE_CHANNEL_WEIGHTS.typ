= OSE and inverse-density channel weights
<ose-and-inverse-density-channel-weights>
The runtime API distinguishes `sampling_channel_weight="ose"` from `"map_density"` (also spelled `"inverse_jacobian"`). OSE changes the positive scores that partition a point between channels. It does not change the maps, replace their integration Jacobians, or select physical cuts. Global and per-channel choices use the same canonical catalogue and estimator.

Implementation and acceptance are staged separately from the completed build12 GL638 runs. The examples below are prepared settings, not claims of completed OSE integration or improved variance.

== Ordinary optimized LMBs
<ordinary-optimized-lmbs>
```toml
[sampling]
graphs = "monte_carlo"
orientations = "summed"
sampling_multichanneling = true
sampling_channels = "monte_carlo"
sampling_channel_weight = "ose"
alpha = 3.0
default_channel_selection = ["auto:optimized_lmb"]

[sampling.channel_selection]
GL638 = ["auto:optimized_lmb"]
```

This selects ordinary optimized LMB channels only and restores their historical OSE prefactor. Keep `sampling_multichanneling=true`: disabling it would disable multichanneling, not merely the advanced maps. An existing graph-specific selection overrides the default, hence the explicit GL638 entry above. The optimized GL638 list already includes the two gluons; this preset adds no second copy of its soft LMB. Change the weight to `map_density` for the same maps with inverse-density prefactors. `alpha` defaults to 3 and can be omitted.

The #link("../../../examples/cli/epem_a_ttxh/NNLO/sampling/optimized_lmb_ose.toml")[runtime overlay] can be applied to a loaded process using `set process -p epem_a_tth -i NNLO file examples/cli/epem_a_ttxh/NNLO/sampling/optimized_lmb_ose.toml`. Use fresh integration workspaces when changing channel weights.

== OSE only for the named soft channel
<ose-only-for-the-named-soft-channel>
For a mixed catalogue retain the global exact-density default, then override the named ordinary LMB channel:

```toml
[sampling]
sampling_channel_weight = "map_density"
alpha = 3.0

[sampling.channel_definitions.GL638.soft_6_12]
around = "lmb(6,12,13,14)"
parent_lmb = [3,4,7,10]
channel_weight = "ose"
```

The #link("../../../examples/cli/epem_a_ttxh/NNLO/sampling/mixed_soft_ose.toml")[mixed runtime overlay] selects all six ordinary cut channels, the hosted joint H/Z channel, the Cut-3 right-threshold channel and the soft LMB. First load those channel definitions, for example from the existing #link("gl638_hosted_joint_gate/energy600_preparation/cards/cuts_joint_cut3_soft_p2_muv50.toml")[600 GeV reference card];, then apply the overlay with the same `set process ... file` command. The physics/stability settings inherited from that reference card must be chosen and validated separately; these overlays only steer sampling. For a paired comparison, change just the soft channel\'s `channel_weight` to `map_density`. All other selected maps and their density scores remain the same.

A named `channel_weight` overrides the global strategy. OSE requires an ordinary complete `lmb(...)` map; applying it directly to a surface or composed map is an error. A mixed catalogue with global `ose` must override each advanced channel explicitly with `channel_weight="map_density"` (or its supplied proxy). Those channels can retain exact map densities or explicitly select their supplied `singularity_proxy`. A global proxy choice still requires an explicit proxy for every channel that inherits it. An OSE choice does not silently replace advanced-map densities with guessed energy factors.

== Score and numerical contract
<score-and-numerical-contract>
For L loop momenta, D=3L raw spatial dimensions, and the selected LMB edges B,

```
rho_OSE(k) = E_cm^(-D) product_{e in B} (E_cm/E_e(k))^alpha.
```

E\_e includes that edge\'s mass and its complete momentum routing at the #strong[unrescaled master raw point];. The existing native affine LMB transformation contains external shifts; its inverse supplies the selected edge vectors. No Cutkosky rescaling or threshold-projected point defines this score. The E\_cm factor gives the same units as an inverse D-dimensional map Jacobian. It cancels between all-LMB channels of common L, reproducing the historical `product(E_e)^(-alpha)` partition. It matters when OSE and map-density scores are compared in a mixed catalogue.

Evaluate logarithmic scores and use the existing log-sum-exp partition. `alpha` is finite and nonnegative; zero gives the constant E\_cm^(-D) before inspecting edge energies. For positive alpha an exact zero massless energy is a singular score requiring an explicit diagnostic, not absent support. No finite floor or hidden cutoff is inserted. A positive full-support LMB score and the existing support-aware advanced scores preserve a nonvanishing partition denominator away from those endpoints.

The selected map\'s actual forward/inverse Jacobian consistency is checked independently whenever its score is OSE or a user proxy. Its partition score cannot be mistaken for an inverse density. Foreign restricted maps continue to contribute only on their actual support. MC channel probabilities, outer sampling weights, signed/absolute accumulation and physical cut/CT summation keep their existing owners.

At an elementary massless soft limit, alpha=3 gives an OSE score proportional to r^-3, compared with the usual linear spherical map density\'s r^-2. This alone does not prove bounded mixed-channel weights: an overly singular score can allocate a joint singular region to a map that does not sample its other normal well. OSE is not a normalized density, so alpha need not be less than 3. Establish numerical benefit through the existing Gaussian acceptance target, paired full-integrand runs and local single-soft/double-soft/HZ scans, including all saved maxima. The original map-density results remain the control.

Internally, `SamplingChannelSelection` owns both `weight` and `alpha`. `MultiChannelingSettings` retains only its parameterization; the obsolete `LmbChannelWeight` enum and lossy conversions are removed. Rust, Python and TOML settings therefore retain the same explicit global and local choices.

== Validation status
<validation-status>
Formatting, `cargo check` and Clippy pass. The focused run at `/tmp/gl638-final-physics-screen/tests-build14-attempt2.log` passed all 97 selected tests, including settings serialization, OSE energy-product ratios with affine shifts, zero/invalid alpha, master-mass refresh, mixed Gaussian normalization, selected actual inverse/Jacobian consistency and proxy support. The first run caught two routine fixture mistakes (per-channel sample counting and a model inconsistent with the generated graph); the corrected tests retain their normalization, mass-ownership and partition assertions.

The broader settings checks also exposed an existing serialization defect: `escalate_if_exact_zero` defaults to false, but its skip condition omitted true. The corrected condition preserves an explicit true across a save/load cycle; the existing roundtrip test now checks it. No stability default changes.

The optimized build14 completed with unchanged sources. Its first full-state capture evaluates 108 cases across six settings arms, each with ordinary rescue and forced Arb precision. All 216 calls are valid, and the independent offline audit passes 84,021 checks with no failures. These cover sampled geometry, actual Jacobians, independently reconstructed OSE scores, positive partitions, and physical component reweighting. The 600 GeV state retains all 936 orientations, six cuts and nineteen CT variants. The capture and audit are `/tmp/gl638-final-physics-screen/ose-validation-build14/results-controls-c1-build14/summary.json` and `audit-controls-c1-build14.json` in the parent directory.

For the retained nine-channel catalogue, these are complete real estimators at the same saved Samples, retaining their original outer weights:

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,right,),
    table.header([Saved point], [Inverse density], [Soft-only OSE, alpha=3],),
    table.hline(),
    [MC calibration negative maximum], [-0.7960873], [-5.778631e-7],
    [MC calibration positive maximum], [0.8969759], [0.9028363],
    [c=1 SUM validation negative maximum], [-0.1534277], [-0.1985041],
    [c=1 SUM validation positive maximum], [0.8718419], [0.9303918],
    [c=4 SUM validation negative maximum], [-0.1907428], [-0.1918194],
    [c=4 SUM validation positive maximum], [0.4689843], [1.630806e-6],
  )]
  , kind: table
  )

OSE therefore does not uniformly improve the maxima. These are pointwise counterfactuals, not a new MC distribution or evidence for improved variance; they do not change the failed status of any originating pilot. The preparation narrative incorrectly associated the c=1 soft partition 0.773 with the positive maximum: it belongs to the negative maximum. The positive maximum has partition 0.937. The underlying archived prediction rows and the native audit agree.

For both K6 and mixed K9, both score choices give approximately linear decay of `J*pi*Re(f)` on each tested single-soft ray and cubic decay on the tested double-soft ray, between represented radii 1e-3 and 1e-6 GeV. Both HZ directions retain their previous behavior: one decays and one approaches a finite plateau over the six tested decades. These local rays do not establish a global bound.

The c=4 capture also passes: eight cases, sixteen native calls and 17,703 independent checks, with no failures. Its artifacts are `results-controls-c4-build14/summary.json` and `audit-controls-c4-build14.json` alongside the c=1 results. No originating integration flag is cleared by these maximum-point checks.

The current Python extension builds, installs and imports successfully. Both new Python settings tests pass in the actual nextest subprocess harness (2/2), against that installed extension. All 3,490 tracked files were unchanged over its execution; see `python-tests.log` and `python-execution.json` in the build14 scratch directory.

The selected full-state native precision boundary also passes: six K6 inverse- density cases, six K6 OSE cases and three mixed K9 OSE cases, with 1,990 physical component comparisons. The K9 cases include an ordinary cut map, the soft LMB and the hosted joint H/Z map. At the tight H/Z point, the optional f64 lane correctly rejects uncertain geometry; f128 and required Arb pass. Sources, state files, input requests and the actual test binary remain unchanged. This is the selected 15-case scope, not the original unrun full K9 decks. See `native-boundary-auth-build14-attempt1/native-boundary-audit.json`. Physical OSE pilots are pending.

The initial production Gaussian checks at width 300 GeV and 8192 draws have no stability flags but do #strong[not] pass normalization acceptance. All eight arms miss the predeclared 3% error cap. Both mixed MC arms also miss the unit target by more than five reported standard errors; high variance is a hypothesis to test, not a demonstrated explanation of that discrepancy.

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,right,),
    table.header([Catalogue / score], [MC normalization], [SUM normalization],),
    table.hline(),
    [Optimized K6 / inverse density], [1.157 +/- 0.350], [0.628 +/- 0.087],
    [Optimized K6 / OSE], [1.119 +/- 0.330], [0.721 +/- 0.124],
    [Mixed K9 / inverse density], [0.391 +/- 0.084], [1.373 +/- 0.407],
    [Mixed K9 / soft OSE], [0.175 +/- 0.075], [4.612 +/- 2.440],
  )]
  , kind: table
  )

These are reference integrals with target one, not GL638 cross sections. The second normalized moment has the same qualitative difficulty and is retained in `audit-gaussian-eight-build14.json`. These cold-grid results are preserved; completion or compatibility within a large error bar is not acceptance.

An independent ordinary-map variance calculation reproduces all six native K6 map points and effective Jacobians to about 1e-13. It finds radial scale `b*E_cm=180 GeV`; a master-frame Gaussian induces substantial correlations in the sampled LMBs. The best tested common reference width, 114 GeV, still needs roughly 1.4--1.8 million cold MC draws for both errors to reach 3%. A fresh native confirmation therefore uses eight adaptive iterations, 4,194,304 MC or 2,097,152 SUM draws per score choice, with fresh seeds and unchanged acceptance criteria. K9 confirmation remains a separate requirement. The protocol and independent calculation are `GAUSSIAN_CONFIRMATION_PROTOCOL.md` and `OFFLINE_GAUSSIAN_VARIANCE.md` in the build14 scratch directory.

== Completed K6 native confirmation
<completed-k6-native-confirmation>
The predeclared larger confirmation now passes all four arms and both unit targets, retaining the original 3% SEM, 6% distance and five-SEM criteria. There are no unstable/NaN samples and the saved state and inputs remain unchanged. The Gaussian width is 114 GeV, with the original shifted center; each arm uses eight adaptive iterations on twelve workers. MC has 4,194,304 total draws, SUM has 2,097,152. These are reference observables, not physical cross sections.

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,right,right,right,),
    table.header([Score / mode], [Normalization], [Normalized second moment], [Native wall \[s\]],),
    table.hline(),
    [Inverse density / MC], [1.012154 +/- .007749], [1.015646 +/- .008100], [402.2],
    [Inverse density / SUM], [1.001629 +/- .005612], [1.000732 +/- .005792], [712.5],
    [OSE / MC], [1.024783 +/- .010876], [1.027223 +/- .010954], [242.1],
    [OSE / SUM], [1.000907 +/- .006119], [.998156 +/- .006135], [376.2],
  )]
  , kind: table
  )

The complete sequence took 1,822.6 seconds and retained all 32 completed iteration checkpoints. Evidence is in `results-gaussian-k6-confirm-build14/summary.json` and `audit-gaussian-k6-confirm-build14.json` in the build14 scratch directory. The ordinary-LMB confirmation resolves that catalogue\'s initial statistical shortfall. It does not by itself resolve the mixed K9 discrepancy or establish physical variance improvement.
