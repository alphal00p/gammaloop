# GL134 legacy global-phase diagnosis

## Result

The GL134 DOT imported by the saved states retained `num = "1𝑖"`, an explicit global multiplier of
the full Feynman numerator. The current physical LU convention already supplies
the cut phase. Replacing only that DOT attribute by `num = "1"` and generating
again makes GL134 real. The advanced sampling maps did not introduce a phase.

The experiment cards' decision to train Im accepted the unexpected phase
instead of checking the physical normalization. That steering and its explanatory
comment were mistakes introduced during the sampling study.

## Inputs and pipeline comparison

Both saved sampling states use the same GL134 input, all orientations, Standard
Model parameters, and 600 GeV kinematics. They differ only in sampling maps.
Momentum-space inspect bypasses those maps and gives the same result for both.
The first different input in the fresh control is the explicit global numerator
prefactor; graph topology, couplings, generation settings, and external data are
unchanged. Fresh `num=i` reproduces the saved states, excluding an obsolete
generated-state phase as the explanation. Both fresh generations have 686 CFF
orientations in the full unfiltered catalogue, with an explicit orientation sum;
the older manual report's count of 64 does not describe these controls.

At source revision `88b72f76e`, the production CLI binary SHA-256 is
`45b5db8c8858e1475113af31a43741bc99a29275e86537541effc75ef2f44d8d`.
The original DOT SHA-256 is
`3f77305f0e16d8499b2b8539d3db47f4db4a6081c40cfd86c1b35df7c8a5d2c3`.
The two fresh controls completed successfully, including ordinary and arbitrary
precision inspections at both points. The table gives arbitrary precision
physical momentum-space results, rounded only for JSON reporting:

| Raw spatial loop momenta (GeV) | Original `num=i` | Corrected `num=1` |
|---|---:|---:|
| `(37,126,-109; -91,64,83; 31,-47,59; 19,23,-17)` | `+3.299785432016905e-41 i` | `+3.299785432016905e-41` |
| `(80,-30,45; -120,15,95; 41,73,-28; -11,52,39)` | `+7.6152751462406e-43 i` | `+7.6152751462406e-43` |

The [archived controls](gl134_phase_20260917/README.md) retain commands, cards,
original inspect outputs, console logs, the DOT-only difference, and exact input
and exported DOTs. Large diagnostic states remain under repository-relative
`target/gl134_phase_diagnostics/`. These controls left the existing saved states,
workspaces, and running integrations unchanged.

## Historical and source evidence

- Before this migration, `NNLO_manual/graphs/GL134.dot:2` retained the explicit
  `i` introduced with the manual audit resources, and
  `NNLO_manual/generate_all_71.toml` requested `--global-prefactor-num "1𝑖"`.
- The inherited phase-fix commit `dd134b710c` explicitly removed that old CLI
  override from the ttH GL20 and GL00 acceptance cards, as well as a ddx card.
  The manual NNLO DOT inputs were not migrated along with those acceptances.
- Current `GlobalPrefactor::default` is one (`numerator/mod.rs`). DOT parsing
  honors an explicit `num` (`graph/global.rs`); dropping arbitrary user factors
  automatically would be incorrect.
- `CutGroup::lu_prefactor` in `processes/cross_section.rs` supplies
  `2*pi*i*(-1)^C_R`, hence `-2*pi*i` for a connected right side. This convention
  is inherited from the phase-fix branch.
- GL134's ten interaction vertices and thirteen internal propagators supply an
  imaginary numerator phase before the explicit graph factor, up to real
  signs. Its four-loop CFF conversion has `i^4=1`. The LU factor completes a
  real cut result when the extra user factor is one.

The full convention is documented in `docs/architecture/phase-conventions.md`.

## Consequences for the existing studies

The same explicit `num=i` occurred in the GL297 and GL638 inputs, including the DOT retained
in the new live GL638 state. Source review confirms that the graph global factor
multiplies the bare term, both UV contributions, and threshold counterterms
consistently: the final-integrand owner combines the local and integrated UV
branches and multiplies the common global atom; threshold variants reuse that
owner. No special UV or right-threshold exception removes this factor.

Thus migration changes the complete result by `F_corrected = -i F_current`:
`Re_corrected = Im_current`, `Im_corrected = -Re_current`. This is an input
normalization correction, not a new sampling or LU-phase convention.

The advanced GL134 three-million-point estimate therefore corresponds to the
physical real value `1.9941979499794005e-9 +/- 3.228128801291338e-12 pb`.
Its old imaginary-part training was effectively training the corrected real
part, so those statistics can be rotated without discarding the samples.

GL638 was instead trained on the old real component. Its quoted real results
and real-only optimization targets consequently require reinterpretation and
fresh optimization of the corrected physical real component. A repair should
migrate the curated DOT prefactors and their source generator, regenerate the
states, correct the experiment-card comments and phase selection, and restart
physical-real training. The curated input migration and card corrections are now
complete. Existing saved states and live integrations have deliberately been left
unchanged; their regeneration and restart are separate actions. Old workspace
data remain evidence of the former phase convention.

The diagnostic controls did not change production inputs or stop the live GL638
run. The subsequent curated input migration does not retroactively change a
loaded state or its running integration.

## Final migrated-card check

The corrected `NNLO_experiment/gl134_optimized_lmbs.toml` was also used to
generate a separate diagnostic state with real-phase training. Its first
physical-momentum probe returns `3.299785432016905e-41 - 0 i`. The
[executed card and output](gl134_phase_20260917/current_card/) preserve this
check; all diagnostic commands exited successfully. Fourteen modified TOML
cards parsed successfully, and each of the eleven reusable DOT edits changes
only the global `num` from `i` to one.

Hashes of all 34 recorded files in the currently running GL638 saved state
remain unchanged. Its service retained the same PID and start time. No saved
state, workspace, or historical result was migrated, and the monitor continues
to report the stored components as recorded. The new DOT inputs take effect
only upon regeneration.
