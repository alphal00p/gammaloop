# Same GL638 spike under frozen current proposals

This physical point can still produce a large weight under the current augmented proposal.

The exact canonical source point agrees with the retained forced-Arb trace to all 300 stored digits. The original unit-weight result is `-1.949730998643977 - 21.31455396636133 i`; the retained complete Sample weight is `159.38755153536746`. Both saved states give the same direct physical integrand, `-7.940506581818668e-37 - 8.680600358517779e-36 i`.

Scope: the new absolute Re and Im maxima at iteration 231 have the same graph/channel/unit coordinates. They remain the saved maxima at iterations 232, 233, 260, 261 and 262. The nearby saved extrema reveal no other comparable new point; this is not an inventory of nonrecord-breaking samples. `nearby_maxima.json` preserves this check.

Frozen checkpoints: optimized iteration 261 (21,381,120 samples), augmented iteration 60 (4,915,200 samples), captured on 2026-09-17 at 08:04–08:05 UTC. The saved grids define the next iteration proposals. No live workspace or card was modified.

The augmented joint H/Z channel has **no support** at this point. All 13 supported production replays record its policy as `Compact { dyadic_exponent: 21 }`. The constructor sets maximum radius `600 × 0.3 = 180 GeV` and normal scale one, giving selected radius `180 / 2^21 = 8.58306884765625e-5 GeV` **for this fixed complement**. The radius is chosen from the conditional geometry and is not a global channel constant: it covers the asymptotic intersection neighborhood at this complement, while this finite-distance point is outside. The physical cut-1 normal radius `sqrt(H²+Z²) = 6.55264499974 GeV` lies far outside that certified disk. The ordinary fallback was not selected. `joint_support_evidence.json` retains the policy record and derivation. The raw map density sum rises by **1.29370433×**; including each frozen adaptive channel probability and continuous grid, total sampling density rises by **2.60054571×**. These compare unequal training histories; the latter is not an isolated new-map effect.

There is no unique production weight at a fixed physical point: the production partition uses raw map density, and the chosen channel then contributes its own learned inverse density. All 19 supported conditional replays are finite and agree with the exact-point prediction within `1.20e-13` relative; rounding inverse coordinates to binary64 changes raw components by at most the value recorded per channel.

| Mode | Channel | Discrete probability | Continuous density | Final Re weight | Final Im weight |
|---|---|---:|---:|---:|---:|
| optimized | 0: lmb[1] | 0.0338816 | 4.91234 | -11.7144767 | -128.063228 |
| optimized | 1: lmb[2] | 0.0358674 | 2.02777 | -26.8075269 | -293.061186 |
| optimized | 2: lmb[13] | 0.0273138 | 0.0493788 | -1445.61263 | -15803.5075 |
| optimized | 3: lmb[14] | 0.819415 | 1.42785 | -1.66643651 | -18.2175648 |
| optimized | 4: lmb[21] | 0.0465147 | 3.87807 | -10.808578 | -118.1599 |
| optimized | 5: lmb[22] | 0.0370075 | 0.951729 | -55.3569069 | -605.164394 |
| augmented | 0: lmb[1] | 0.164406 | 0.290813 | -31.5215903 | -344.595556 |
| augmented | 1: lmb[2] | 0.122821 | 0.3342 | -36.7164488 | -401.386001 |
| augmented | 2: lmb[13] | 0.062459 | 0.393531 | -61.3148429 | -670.29684 |
| augmented | 3: lmb[14] | 0.0616257 | 0.173686 | -140.803177 | -1539.26717 |
| augmented | 4: lmb[21] | 0.179692 | 0.363936 | -23.0454269 | -251.933726 |
| augmented | 5: lmb[22] | 0.170057 | 1.83587 | -4.8272806 | -52.7720659 |
| augmented | 6: lu_cut_0 | 0.0192283 | 2.215 | -35.3854443 | -386.835396 |
| augmented | 7: lu_cut_1 | 0.0566282 | 17.1471 | -1.55208443 | -16.9674624 |
| augmented | 8: lu_cut_2 | 0.00722081 | 27.1125 | -7.69810771 | -84.1560874 |
| augmented | 9: lu_cut_3 | 0.041019 | 41.4059 | -0.887344718 | -9.70049555 |
| augmented | 10: lu_cut_4 | 0.0195237 | 27.2663 | -2.83108514 | -30.9495603 |
| augmented | 11: lu_cut_5 | 0.00623845 | 33.5199 | -7.20710484 | -78.7884201 |
| augmented | 12: cut1_joint_HZ | — | No support | — | — |
| augmented | 13: cut3_right_5_10 | 0.0830909 | 12.1922 | -1.48766063 | -16.2631783 |

The augmented checkpoint has recorded maxima `|Re| = 4.7634260343` and `|Im| = 8.2648242984`. The counterfactual original-channel weights below exceed those records by **29.56×** and **186.24×**.

The augmented original channel 3 still gives `-140.803177 - 1539.267171 i`, only 2.207× smaller than the original spike. At its frozen sample count, one such sample has mean/error scale `2.86465e-5` in Re and `3.13165e-4` in Im, so a comparable draw could still cause a substantial jump. This channel represents about 0.12925% of the augmented local sampling density at this raw point; that is a local channel share, not the global probability of hitting the point.

Current optimized channel 3 has already learned the spike and now gives only `-1.66643651 - 18.21756485 i` (186.48× lower than the original real weight). Its channel 2 remains a more severe conditional tail at this point. Therefore these replays do not support a claim that the current augmented proposal eliminates this spike.

Production formula audit: let `q_c = 1/J_c`, `Q = sum_c q_c`, `p_c` be the learned discrete probability, and `h_c` the learned cube density. The bridge uses `alpha_c = q_c/Q`; `Grid::probe` returns `1/(p_c h_c)` (the only graph has probability one). Therefore the returned complete estimator is `W_c = f J_c alpha_c/(p_c h_c) = f/(Q p_c h_c)`. All supported rows satisfy both raw-partition identities within `3.34e-16`. The production owners are `sampling_selection.rs::partition_with_contexts`, `SamplingChannelBridgeEvaluation::selected_factor`, and the existing Havana `Grid::probe`.

For comparison only, replacing the production partition by the full adaptive mixture partition would give `-11.9990997 - 131.174740 i` (optimized) or `-4.61406989 - 50.4412362 i` (augmented). **Those are hypothetical estimators, not the production results above.** This is the conditional mean over channel identity at the fixed existing mixture proposal (a Rao–Blackwell improvement candidate), and would require a separate production change and validation; no such change is made here.

`inverse_replay.json` contains all inverse unit coordinates, map densities, partitions, complete learned weights, raw controls, native metadata, and production replays. `frozen_checkpoint.json` and `build_provenance.json` retain checkpoint and executable provenance.

Reproduction uses the existing generated GL638 payloads and compatible development libraries, which are not duplicated here. `generated_state_hashes.json` identifies those read-only payloads; `build_provenance.json` identifies the exact libraries and executable used. No binary integration checkpoint is needed: the compressed grid exports preserve the complete learned grids, and `frozen_inputs.json` preserves iteration metadata and the original complete maximum Sample.

From the repository root, enter the repository development environment with its compatible Symbolica license, then run the following Python snippet. The saved `build_argv.json` resolves the original host's exact library filenames; a rebuilt library set must use matching dependency identities and be recorded as a new replay. The source and archive are left unchanged; outputs go to `target/gl638_spike_inverse_reproduction`.

```python
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess

archive = Path("examples/cli/epem_a_ttxh/NNLO/RESEARCH/gl638_spike_20260917/inverse_current")
output = Path("target/gl638_spike_inverse_reproduction")
output.mkdir(parents=True, exist_ok=True)
for source in archive.rglob("*"):
    if source.is_file():
        target = output / source.relative_to(archive)
        target.parent.mkdir(parents=True, exist_ok=True)
        if source.suffix == ".gz":
            target.with_suffix("").write_bytes(gzip.decompress(source.read_bytes()))
        else:
            shutil.copyfile(source, target)
argv = json.loads((archive / "build_argv.json").read_text())
argv[-3] = str(archive / "inverse_replay.rs")
argv[-1] = str(output / "inverse_replay")
subprocess.run(argv, check=True)
env = dict(os.environ, GL_DISPLAY_FILTER="off", GL_LOGFILE_FILTER="off")
subprocess.run([argv[-1], str(output)], env=env, check=True)
subprocess.run([".venv/bin/python", str(output / "summarize.py")], check=True)
```

`inverse_replay.json.gz` retains all inverse coordinates, including high-precision values, and native production evaluations. The client rederives the canonical Quad point from the original cube and requires exact equality with the retained 300-digit raw trace before proceeding. All inverse diagnostics use 1000-bit Arb; production source policies remain Quad for optimized and Fixed256 for augmented. The actual production cube is binary64, and its raw round-trip deviation is at most `1.14e-13 GeV`. Both runtime direct-momentum controls agree exactly at the reported binary64 boundary.
