# 3L aa -> aa Performance Smoke

Date: 2026-06-15

This file records the first whole-amplitude three-loop `aa -> aa` smoke test using the lightweight `SingleParametric` runtime-evaluation path. Compilation, summed evaluators, and summed function maps are disabled in the card.

## Card and Artifacts

- Card: `examples/cli/aa_aa/3L/aa_aa.toml`
- DOT export: `examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL000.dot` through `GL338.dot`
- DOT count: 339 graph files
- Runtime evaluator: `SingleParametric`
- Function maps/compilation: disabled (`compile=false`, `summed=false`, `summed_function_map=false`)
- UV/threshold/tropical generation: disabled for this raw runtime smoke test
- Kinematics/helicity: kinematics A, helicity `+-+-`
- Final card integration smoke settings after this pass: `n_start=20`, `n_max=20`, `--batch-size 1`

## Diagram Generation

The card uses feyngen with

```text
generate amp a a > a a | a t t~ g ghG ghG~ QCD==4 QED==4 [{3}]
  --numerator-grouping no_grouping
  --symmetrize-initial-states=true --symmetrize-final-states=true
  -p aa_aa -i 3L --only-diagrams
```

Observed counts from `/tmp/aa_aa_3l_full_sequence_nolimit.log`:

| Stage | Graphs |
|---|---:|
| Symbolica generation | 6192 |
| After vetoed topologies | 4284 |
| After complete-graph filters | 1296 |
| After closed-fermion-chain analysis | 1296 |
| After external-state symmetrization | 339 |
| After canonization | 339 |
| After numerator-aware grouping (`no_grouping`) | 339 |

The 339 graphs contain 189 isomorphically unique graph shapes. The sum of generated graph symmetry factors reported by feyngen is `-444`.

A previous attempt with numerator-aware scalar-rescaling grouping reached the same canonical 339-graph stage but was not kept as the default because it was much heavier. The current 3L card deliberately uses `no_grouping` to keep this first all-graph generation pass tractable.

### Release-Build Scalar-Rescaling Grouping Attempt

A follow-up card was added at `examples/cli/aa_aa/3L/aa_aa_grouped.toml` to isolate the one-time feyngen pass with

```text
--numerator-grouping group_identical_graphs_up_to_scalar_rescaling
```

The card saves DOTs to `examples/cli/aa_aa/3L/graphs_grouped` and would save the grouped state to `examples/cli/aa_aa/3L/gammaloop_state_grouped` if the feyngen grouping pass completed.

Command used with a freshly built release binary:

```bash
target/release/gammaloop --clean-state examples/cli/aa_aa/3L/aa_aa_grouped.toml run generate_diagrams
```

The release binary was built with:

```bash
nix develop --no-write-lock-file -c cargo build -p gammaloop-api --bin gammaloop --features ufo_support --release
```

Observed behavior from `/tmp/aa_aa_3l_grouped_release.log` and process monitoring:

| Quantity | Value |
|---|---:|
| Symbolica generation | 6192 graphs |
| After vetoed topologies | 4284 graphs |
| After complete-graph filters | 1296 graphs |
| After closed-fermion-chain analysis | 1296 graphs |
| After external-state symmetrization | 339 graphs |
| After canonization | 339 graphs |
| Time spent after canonization before manual stop | about 80 min |
| Peak observed RSS during grouping | about 166 GiB |
| CPU usage during grouping | about 9-10 cores |
| Grouped DOT files produced | 0 |
| Usable grouped saved state produced | no |

The run was interrupted after about 80 minutes because it remained inside the scalar-rescaling numerator-aware grouping stage with no progress log after the 339 canonized graph line. A tiny startup state scaffold was removed afterward to avoid confusing it with a usable grouped saved state. This confirms that using the release binary does not make the scalar-rescaling grouping pass tractable for the full 3L `aa -> aa` graph set as currently configured.

## Integrand Generation

Command used for the no-cap run:

```bash
./gammaloop --clean-state -s /tmp/gammaloop_aa_aa_3l_perf_nolimit /tmp/aa_aa_3l_noauto.toml run generate_diagrams generate_integrands integrate_diagrams
```

Generation summary:

| Quantity | Value |
|---|---:|
| Wall time, `generate existing` | about 29m 58s |
| Peak RAM reported by GammaLoop | 9.59 GiB |
| Generation cores | 10 |
| Per-graph timing sum | 291.11 min |
| Mean graph generation time | 51.52 s |
| Median graph generation time | 43.30 s |
| Fastest graph | 119 ms |
| Slowest graph | 3.15 min |
| Graphs above 60 s | 112 |
| Graphs above 120 s | 25 |
| Graphs below 1 s | 38 |

Timing share over the 339 per-graph timings:

| Phase | Sum | Share |
|---|---:|---:|
| Expression build | 277.33 min | 95.27% |
| Spenso | 6.32 min | 2.17% |
| Symbolica eval | 7.46 min | 2.56% |
| Compile | 0 ms | 0.00% |

Expression construction is the dominant generation cost. Spenso and Symbolica runtime evaluation are subdominant in this light `SingleParametric` setup.

## Whole-Amplitude Integration Smoke

Two all-graph integration attempts were made.

1. With a 100 GiB virtual-memory cap, integrand generation completed, but integration aborted before the first iteration with `memory allocation of 16896 bytes failed`. The process had reached the virtual-memory cap during integration setup and still held about 77 GiB RSS when stopped.
2. Without the cap, the same setup reached actual integration over 3 nested discrete grids, 339 graphs, and 9 continuous dimensions using Monte Carlo over graphs, orientations, and LMBs. Integration-time memory rose well beyond the capped limit: about 51 GiB RSS at integration setup, then above 570 GiB RSS before any 1000-sample iteration completed. The run was stopped manually after roughly 7 minutes of sampling to avoid unnecessary pressure on the shared machine.

No completed integration iteration was obtained from the 1000-sample, `--batch-size 100` smoke run, so there is no reliable MC uncertainty or ms/sample/core number yet for the full 339-graph amplitude. The immediate performance conclusion is that all-graph raw 3L `SingleParametric` integration is memory dominated before it is statistics dominated. The checked-in card now uses `n_start=20`, `n_max=20`, and `--batch-size 1` for the next safer whole-amplitude smoke attempt.

## Slowest Generation Graphs

| Rank | Graph | Expr build | Spenso | Symbolica eval | Total |
|---:|---|---:|---:|---:|---:|
| 1 | GL296 | 3m 0.9s | 3.95s | 4.08s | 3.15 min |
| 2 | GL299 | 3m 2.0s | 2.86s | 2.54s | 3.12 min |
| 3 | GL191 | 2m 30.2s | 2.63s | 2.88s | 2.60 min |
| 4 | GL169 | 2m 28.5s | 3.19s | 3.17s | 2.58 min |
| 5 | GL179 | 2m 30.1s | 2.18s | 2.49s | 2.58 min |
| 6 | GL213 | 2m 24.6s | 2.25s | 2.29s | 2.49 min |
| 7 | GL192 | 2m 21.0s | 2.40s | 2.49s | 2.43 min |
| 8 | GL216 | 2m 17.9s | 3.15s | 2.44s | 2.39 min |
| 9 | GL219 | 2m 18.6s | 2.30s | 2.24s | 2.39 min |
| 10 | GL226 | 2m 17.3s | 2.56s | 2.63s | 2.37 min |
| 11 | GL215 | 2m 16.7s | 2.49s | 2.46s | 2.36 min |
| 12 | GL178 | 2m 16.4s | 2.41s | 2.72s | 2.36 min |
| 13 | GL170 | 2m 15.1s | 2.27s | 2.78s | 2.34 min |
| 14 | GL186 | 2m 14.8s | 2.69s | 2.59s | 2.33 min |
| 15 | GL220 | 2m 15.2s | 2.46s | 2.33s | 2.33 min |
| 16 | GL187 | 2m 14.0s | 2.75s | 2.59s | 2.32 min |
| 17 | GL214 | 2m 12.3s | 2.75s | 2.31s | 2.29 min |
| 18 | GL198 | 2m 10.6s | 2.63s | 2.34s | 2.26 min |
| 19 | GL197 | 2m 8.2s | 2.84s | 2.26s | 2.22 min |
| 20 | GL309 | 2m 6.0s | 2.78s | 2.09s | 2.18 min |

## Per-Graph Generation Breakdown

| Graph | Expr build | Spenso | Symbolica eval | Total |
|---|---:|---:|---:|---:|
| GL000 | 22.18s | 1.15s | 1.08s | 24.41 s |
| GL001 | 32.26s | 700ms | 1.23s | 34.19 s |
| GL002 | 32.72s | 659ms | 1.23s | 34.61 s |
| GL003 | 22.35s | 673ms | 991ms | 24.01 s |
| GL004 | 20.86s | 635ms | 1.06s | 22.55 s |
| GL005 | 30.60s | 671ms | 1.21s | 32.48 s |
| GL006 | 29.89s | 667ms | 1.16s | 31.72 s |
| GL007 | 20.91s | 735ms | 973ms | 22.62 s |
| GL008 | 20.71s | 652ms | 939ms | 22.30 s |
| GL009 | 30.11s | 767ms | 1.11s | 31.99 s |
| GL010 | 23.74s | 808ms | 1.07s | 25.62 s |
| GL011 | 33.27s | 912ms | 1.29s | 35.47 s |
| GL012 | 32.51s | 1.05s | 1.39s | 34.95 s |
| GL013 | 32.56s | 1.28s | 1.35s | 35.19 s |
| GL014 | 40.29s | 1.38s | 1.63s | 43.30 s |
| GL015 | 134ms | 0ms | 0ms | 134 ms |
| GL016 | 119ms | 0ms | 0ms | 119 ms |
| GL017 | 119ms | 0ms | 0ms | 119 ms |
| GL018 | 39.82s | 1.26s | 1.55s | 42.63 s |
| GL019 | 39.95s | 1.29s | 1.54s | 42.78 s |
| GL020 | 33.02s | 1.07s | 1.26s | 35.35 s |
| GL021 | 135ms | 0ms | 0ms | 135 ms |
| GL022 | 46.47s | 2.19s | 1.45s | 50.11 s |
| GL023 | 38.50s | 1.05s | 1.28s | 40.83 s |
| GL024 | 37.54s | 1.14s | 1.26s | 39.94 s |
| GL025 | 35.35s | 1.44s | 1.24s | 38.03 s |
| GL026 | 22.69s | 641ms | 775ms | 24.11 s |
| GL027 | 36.76s | 780ms | 1.10s | 38.64 s |
| GL028 | 720ms | 0ms | 0ms | 720 ms |
| GL029 | 33.18s | 742ms | 896ms | 34.82 s |
| GL030 | 31.63s | 733ms | 895ms | 33.26 s |
| GL031 | 31.32s | 1.58s | 1.13s | 34.03 s |
| GL032 | 31.88s | 821ms | 1.11s | 33.81 s |
| GL033 | 31.10s | 768ms | 1.06s | 32.93 s |
| GL034 | 31.29s | 781ms | 1.07s | 33.14 s |
| GL035 | 33.14s | 919ms | 997ms | 35.06 s |
| GL036 | 32.19s | 742ms | 998ms | 33.93 s |
| GL037 | 36.22s | 787ms | 1.12s | 38.13 s |
| GL038 | 37.02s | 812ms | 1.18s | 39.01 s |
| GL039 | 630ms | 0ms | 0ms | 630 ms |
| GL040 | 608ms | 0ms | 0ms | 608 ms |
| GL041 | 613ms | 0ms | 0ms | 613 ms |
| GL042 | 644ms | 0ms | 0ms | 644 ms |
| GL043 | 40.45s | 1.32s | 1.29s | 43.06 s |
| GL044 | 37.91s | 846ms | 1.23s | 39.99 s |
| GL045 | 36.70s | 823ms | 994ms | 38.52 s |
| GL046 | 35.00s | 806ms | 948ms | 36.75 s |
| GL047 | 33.82s | 819ms | 894ms | 35.53 s |
| GL048 | 33.51s | 740ms | 885ms | 35.13 s |
| GL049 | 41.90s | 891ms | 1.23s | 44.02 s |
| GL050 | 667ms | 0ms | 0ms | 667 ms |
| GL051 | 33.17s | 750ms | 1.11s | 35.03 s |
| GL052 | 33.20s | 716ms | 1.14s | 35.06 s |
| GL053 | 27.94s | 712ms | 909ms | 29.56 s |
| GL054 | 27.74s | 641ms | 949ms | 29.33 s |
| GL055 | 688ms | 0ms | 0ms | 688 ms |
| GL056 | 689ms | 0ms | 0ms | 689 ms |
| GL057 | 34.26s | 879ms | 1.13s | 36.27 s |
| GL058 | 34.70s | 819ms | 1.26s | 36.78 s |
| GL059 | 27.29s | 777ms | 993ms | 29.06 s |
| GL060 | 28.11s | 662ms | 927ms | 29.70 s |
| GL061 | 37.22s | 930ms | 1.02s | 39.17 s |
| GL062 | 35.35s | 737ms | 1.15s | 37.24 s |
| GL063 | 632ms | 0ms | 0ms | 632 ms |
| GL064 | 647ms | 0ms | 0ms | 647 ms |
| GL065 | 667ms | 0ms | 0ms | 667 ms |
| GL066 | 670ms | 0ms | 0ms | 670 ms |
| GL067 | 34.38s | 1.57s | 1.32s | 37.27 s |
| GL068 | 35.10s | 888ms | 1.30s | 37.29 s |
| GL069 | 33.34s | 801ms | 1.08s | 35.22 s |
| GL070 | 33.58s | 728ms | 1.12s | 35.43 s |
| GL071 | 35.39s | 905ms | 1.07s | 37.37 s |
| GL072 | 33.88s | 710ms | 1.11s | 35.70 s |
| GL073 | 34.95s | 791ms | 1.00s | 36.74 s |
| GL074 | 34.81s | 763ms | 1.14s | 36.71 s |
| GL075 | 652ms | 0ms | 0ms | 652 ms |
| GL076 | 634ms | 0ms | 0ms | 634 ms |
| GL077 | 654ms | 0ms | 0ms | 654 ms |
| GL078 | 686ms | 0ms | 0ms | 686 ms |
| GL079 | 35.62s | 879ms | 950ms | 37.45 s |
| GL080 | 35.23s | 789ms | 1.09s | 37.11 s |
| GL081 | 28.45s | 783ms | 897ms | 30.13 s |
| GL082 | 28.61s | 657ms | 916ms | 30.18 s |
| GL083 | 711ms | 0ms | 0ms | 711 ms |
| GL084 | 801ms | 0ms | 2ms | 803 ms |
| GL085 | 41.87s | 1.73s | 1.28s | 44.88 s |
| GL086 | 39.24s | 847ms | 1.33s | 41.42 s |
| GL087 | 30.63s | 768ms | 1.12s | 32.52 s |
| GL088 | 31.43s | 683ms | 1.04s | 33.15 s |
| GL089 | 35.76s | 957ms | 1.24s | 37.96 s |
| GL090 | 34.98s | 763ms | 1.29s | 37.03 s |
| GL091 | 36.75s | 863ms | 1.23s | 38.84 s |
| GL092 | 35.59s | 724ms | 1.28s | 37.59 s |
| GL093 | 29.15s | 712ms | 976ms | 30.84 s |
| GL094 | 29.11s | 678ms | 1.10s | 30.89 s |
| GL095 | 40.93s | 979ms | 1.20s | 43.11 s |
| GL096 | 41.30s | 988ms | 1.28s | 43.57 s |
| GL097 | 30.62s | 878ms | 1.05s | 32.55 s |
| GL098 | 30.22s | 773ms | 1.04s | 32.03 s |
| GL099 | 52.10s | 1.07s | 2.24s | 55.41 s |
| GL100 | 52.27s | 1.31s | 2.38s | 55.96 s |
| GL101 | 52.11s | 1.23s | 2.73s | 56.07 s |
| GL102 | 51.91s | 1.35s | 2.28s | 55.54 s |
| GL103 | 53.26s | 1.26s | 2.01s | 56.53 s |
| GL104 | 1m 0.1s | 1.38s | 2.79s | 64.27 s |
| GL105 | 1m 18.1s | 1.56s | 1.96s | 81.62 s |
| GL106 | 1m 19.2s | 1.53s | 2.09s | 82.82 s |
| GL107 | 1.02s | 0ms | 0ms | 1.02 s |
| GL108 | 1m 37.8s | 1.88s | 2.87s | 102.55 s |
| GL109 | 1m 18.0s | 1.64s | 2.18s | 81.82 s |
| GL110 | 1m 21.2s | 1.60s | 2.17s | 84.97 s |
| GL111 | 1m 34.8s | 1.97s | 2.59s | 99.36 s |
| GL112 | 1m 32.5s | 1.75s | 2.60s | 96.85 s |
| GL113 | 1.02s | 0ms | 0ms | 1.02 s |
| GL114 | 1.02s | 0ms | 0ms | 1.02 s |
| GL115 | 1.05s | 0ms | 0ms | 1.05 s |
| GL116 | 1.02s | 0ms | 0ms | 1.02 s |
| GL117 | 1m 32.4s | 2.03s | 2.56s | 96.99 s |
| GL118 | 1m 32.9s | 2.31s | 2.54s | 97.75 s |
| GL119 | 1m 18.6s | 1.58s | 2.29s | 82.47 s |
| GL120 | 1m 17.8s | 1.65s | 2.09s | 81.54 s |
| GL121 | 1.07s | 0ms | 0ms | 1.07 s |
| GL122 | 1m 37.1s | 1.74s | 2.36s | 101.20 s |
| GL123 | 1m 21.2s | 1.38s | 2.11s | 84.69 s |
| GL124 | 1m 21.2s | 1.53s | 2.13s | 84.86 s |
| GL125 | 1m 21.6s | 1.60s | 2.08s | 85.28 s |
| GL126 | 1m 25.2s | 2.51s | 2.45s | 90.16 s |
| GL127 | 1m 27.4s | 1.72s | 2.35s | 91.47 s |
| GL128 | 1m 30.1s | 1.74s | 2.12s | 93.96 s |
| GL129 | 52.61s | 1.20s | 1.49s | 55.30 s |
| GL130 | 53.05s | 1.12s | 2.45s | 56.62 s |
| GL131 | 51.99s | 993ms | 1.48s | 54.46 s |
| GL132 | 51.48s | 1.08s | 2.22s | 54.78 s |
| GL133 | 51.94s | 1.29s | 1.49s | 54.72 s |
| GL134 | 58.59s | 1.33s | 1.56s | 61.48 s |
| GL135 | 17.24s | 465ms | 692ms | 18.40 s |
| GL136 | 17.74s | 499ms | 845ms | 19.08 s |
| GL137 | 17.52s | 514ms | 715ms | 18.75 s |
| GL138 | 17.59s | 482ms | 722ms | 18.79 s |
| GL139 | 17.45s | 638ms | 823ms | 18.91 s |
| GL140 | 20.14s | 551ms | 832ms | 21.52 s |
| GL141 | 41.21s | 1.43s | 1.79s | 44.43 s |
| GL142 | 42.33s | 1.22s | 1.72s | 45.27 s |
| GL143 | 1m 3.9s | 1.47s | 1.45s | 66.82 s |
| GL144 | 1m 51.5s | 2.31s | 2.50s | 116.31 s |
| GL145 | 1m 54.5s | 2.23s | 2.98s | 119.71 s |
| GL146 | 1m 51.2s | 2.75s | 2.95s | 116.90 s |
| GL147 | 1m 4.6s | 1.42s | 1.61s | 67.63 s |
| GL148 | 44.11s | 1.49s | 1.82s | 47.42 s |
| GL149 | 41.75s | 1.44s | 1.86s | 45.05 s |
| GL150 | 41.30s | 1.27s | 1.82s | 44.39 s |
| GL151 | 33.43s | 1.32s | 1.61s | 36.36 s |
| GL152 | 50.51s | 1.10s | 1.69s | 53.30 s |
| GL153 | 34.66s | 943ms | 1.86s | 37.46 s |
| GL154 | 36.14s | 1.42s | 1.89s | 39.45 s |
| GL155 | 37.10s | 1.71s | 1.94s | 40.75 s |
| GL156 | 55.48s | 1.08s | 1.57s | 58.13 s |
| GL157 | 1m 3.2s | 1.44s | 1.35s | 65.99 s |
| GL158 | 1m 4.5s | 1.20s | 1.34s | 67.04 s |
| GL159 | 1m 34.7s | 1.89s | 1.86s | 98.45 s |
| GL160 | 1m 33.5s | 2.23s | 1.86s | 97.59 s |
| GL161 | 35.89s | 1.19s | 1.56s | 38.64 s |
| GL162 | 34.84s | 1.06s | 1.57s | 37.47 s |
| GL163 | 35.73s | 1.07s | 1.66s | 38.46 s |
| GL164 | 35.47s | 1.07s | 1.68s | 38.22 s |
| GL165 | 45.69s | 1.36s | 1.16s | 48.21 s |
| GL166 | 2m 0.9s | 2.18s | 2.73s | 2.10 min |
| GL167 | 998ms | 0ms | 0ms | 998 ms |
| GL168 | 972ms | 0ms | 0ms | 972 ms |
| GL169 | 2m 28.5s | 3.19s | 3.17s | 2.58 min |
| GL170 | 2m 15.1s | 2.27s | 2.78s | 2.34 min |
| GL171 | 55.57s | 1.20s | 1.45s | 58.22 s |
| GL172 | 55.83s | 1.09s | 1.47s | 58.39 s |
| GL173 | 57.53s | 1.67s | 1.82s | 61.02 s |
| GL174 | 56.29s | 1.70s | 1.75s | 59.74 s |
| GL175 | 58.02s | 1.50s | 1.79s | 61.31 s |
| GL176 | 58.40s | 1.15s | 1.50s | 61.05 s |
| GL177 | 56.70s | 1.01s | 1.74s | 59.45 s |
| GL178 | 2m 16.4s | 2.41s | 2.72s | 2.36 min |
| GL179 | 2m 30.1s | 2.18s | 2.49s | 2.58 min |
| GL180 | 1m 1.1s | 1.51s | 1.82s | 64.43 s |
| GL181 | 1m 1.2s | 1.52s | 1.84s | 64.56 s |
| GL182 | 58.36s | 2.18s | 1.84s | 62.38 s |
| GL183 | 58.47s | 2.07s | 1.87s | 62.41 s |
| GL184 | 55.51s | 1.04s | 1.39s | 57.94 s |
| GL185 | 55.03s | 1.17s | 1.56s | 57.76 s |
| GL186 | 2m 14.8s | 2.69s | 2.59s | 2.33 min |
| GL187 | 2m 14.0s | 2.75s | 2.59s | 2.32 min |
| GL188 | 56.20s | 1.58s | 1.91s | 59.69 s |
| GL189 | 48.03s | 1.05s | 1.10s | 50.18 s |
| GL190 | 57.69s | 1.40s | 1.13s | 60.22 s |
| GL191 | 2m 30.2s | 2.63s | 2.88s | 2.60 min |
| GL192 | 2m 21.0s | 2.40s | 2.49s | 2.43 min |
| GL193 | 963ms | 0ms | 0ms | 963 ms |
| GL194 | 1.01s | 0ms | 0ms | 1.01 s |
| GL195 | 954ms | 0ms | 0ms | 954 ms |
| GL196 | 998ms | 0ms | 0ms | 998 ms |
| GL197 | 2m 8.2s | 2.84s | 2.26s | 2.22 min |
| GL198 | 2m 10.6s | 2.63s | 2.34s | 2.26 min |
| GL199 | 1m 1.5s | 1.17s | 1.28s | 63.95 s |
| GL200 | 1m 1.6s | 1.27s | 1.74s | 64.61 s |
| GL201 | 56.55s | 2.18s | 2.19s | 60.92 s |
| GL202 | 57.42s | 1.61s | 2.12s | 61.15 s |
| GL203 | 56.65s | 1.58s | 2.17s | 60.40 s |
| GL204 | 58.07s | 1.69s | 1.90s | 61.66 s |
| GL205 | 59.97s | 1.12s | 1.37s | 62.46 s |
| GL206 | 1m 1.0s | 1.78s | 1.61s | 64.39 s |
| GL207 | 48.18s | 954ms | 1.14s | 50.27 s |
| GL208 | 48.76s | 1.21s | 1.15s | 51.12 s |
| GL209 | 964ms | 0ms | 0ms | 964 ms |
| GL210 | 969ms | 0ms | 0ms | 969 ms |
| GL211 | 1.10s | 0ms | 0ms | 1.10 s |
| GL212 | 938ms | 0ms | 0ms | 938 ms |
| GL213 | 2m 24.6s | 2.25s | 2.29s | 2.49 min |
| GL214 | 2m 12.3s | 2.75s | 2.31s | 2.29 min |
| GL215 | 2m 16.7s | 2.49s | 2.46s | 2.36 min |
| GL216 | 2m 17.9s | 3.15s | 2.44s | 2.39 min |
| GL217 | 1m 3.4s | 1.46s | 1.38s | 66.24 s |
| GL218 | 1m 3.1s | 1.19s | 1.60s | 65.89 s |
| GL219 | 2m 18.6s | 2.30s | 2.24s | 2.39 min |
| GL220 | 2m 15.2s | 2.46s | 2.33s | 2.33 min |
| GL221 | 56.55s | 1.55s | 1.90s | 60.00 s |
| GL222 | 56.04s | 1.67s | 1.91s | 59.62 s |
| GL223 | 57.70s | 2.09s | 1.93s | 61.72 s |
| GL224 | 56.12s | 1.85s | 1.94s | 59.91 s |
| GL225 | 53.25s | 1.15s | 1.20s | 55.60 s |
| GL226 | 2m 17.3s | 2.56s | 2.63s | 2.37 min |
| GL227 | 1m 3.4s | 1.55s | 2.14s | 67.09 s |
| GL228 | 58.39s | 2.08s | 2.23s | 62.70 s |
| GL229 | 1.13s | 0ms | 0ms | 1.13 s |
| GL230 | 1.11s | 0ms | 0ms | 1.11 s |
| GL231 | 53.71s | 1.10s | 1.37s | 56.18 s |
| GL232 | 1m 30.0s | 1.45s | 2.24s | 93.69 s |
| GL233 | 1m 24.3s | 1.54s | 1.80s | 87.64 s |
| GL234 | 2.88s | 0ms | 0ms | 2.88 s |
| GL235 | 1m 38.6s | 2.35s | 2.36s | 103.31 s |
| GL236 | 1m 18.5s | 1.81s | 1.82s | 82.13 s |
| GL237 | 1m 17.9s | 1.63s | 1.93s | 81.46 s |
| GL238 | 1m 33.1s | 1.85s | 2.32s | 97.27 s |
| GL239 | 1m 36.0s | 1.86s | 2.02s | 99.88 s |
| GL240 | 2.73s | 0ms | 0ms | 2.73 s |
| GL241 | 2.65s | 0ms | 0ms | 2.65 s |
| GL242 | 2.64s | 0ms | 0ms | 2.64 s |
| GL243 | 2.62s | 0ms | 0ms | 2.62 s |
| GL244 | 1m 33.5s | 1.81s | 2.29s | 97.60 s |
| GL245 | 1m 35.1s | 1.98s | 2.03s | 99.11 s |
| GL246 | 1m 18.4s | 1.65s | 2.06s | 82.11 s |
| GL247 | 1m 20.2s | 1.42s | 1.78s | 83.40 s |
| GL248 | 2.84s | 0ms | 0ms | 2.84 s |
| GL249 | 1m 42.8s | 1.88s | 2.25s | 106.93 s |
| GL250 | 1m 23.6s | 1.48s | 1.87s | 86.95 s |
| GL251 | 1m 30.5s | 1.65s | 2.02s | 94.17 s |
| GL252 | 1m 25.7s | 1.57s | 2.44s | 89.71 s |
| GL253 | 1m 23.6s | 1.48s | 1.66s | 86.74 s |
| GL254 | 1m 37.9s | 2.58s | 2.17s | 102.65 s |
| GL255 | 1m 32.6s | 1.71s | 2.05s | 96.36 s |
| GL256 | 31.60s | 907ms | 936ms | 33.44 s |
| GL257 | 32.41s | 670ms | 896ms | 33.98 s |
| GL258 | 30.84s | 710ms | 847ms | 32.40 s |
| GL259 | 29.85s | 699ms | 846ms | 31.40 s |
| GL260 | 30.77s | 652ms | 883ms | 32.30 s |
| GL261 | 34.94s | 848ms | 937ms | 36.72 s |
| GL262 | 29.96s | 710ms | 849ms | 31.52 s |
| GL263 | 29.82s | 759ms | 917ms | 31.50 s |
| GL264 | 29.44s | 689ms | 897ms | 31.03 s |
| GL265 | 29.73s | 776ms | 894ms | 31.40 s |
| GL266 | 31.26s | 851ms | 854ms | 32.97 s |
| GL267 | 35.67s | 805ms | 903ms | 37.38 s |
| GL268 | 1m 33.3s | 1.69s | 2.39s | 97.38 s |
| GL269 | 1m 32.0s | 1.65s | 2.13s | 95.78 s |
| GL270 | 34.52s | 1.11s | 1.57s | 37.20 s |
| GL271 | 35.04s | 958ms | 1.61s | 37.61 s |
| GL272 | 34.89s | 1.26s | 1.59s | 37.74 s |
| GL273 | 35.50s | 1.23s | 1.70s | 38.43 s |
| GL274 | 53.17s | 1.28s | 1.49s | 55.94 s |
| GL275 | 52.80s | 991ms | 1.42s | 55.21 s |
| GL276 | 1m 56.1s | 2.15s | 2.56s | 2.01 min |
| GL277 | 42.47s | 1.43s | 1.69s | 45.59 s |
| GL278 | 43.01s | 1.38s | 1.68s | 46.07 s |
| GL279 | 44.59s | 1.16s | 1.65s | 47.40 s |
| GL280 | 1m 4.1s | 1.15s | 1.57s | 66.82 s |
| GL281 | 1m 4.3s | 1.16s | 1.58s | 67.04 s |
| GL282 | 1m 29.9s | 1.76s | 2.17s | 93.83 s |
| GL283 | 1m 37.1s | 1.73s | 2.20s | 101.03 s |
| GL284 | 52.86s | 1.00s | 1.24s | 55.10 s |
| GL285 | 51.91s | 1.11s | 1.41s | 54.43 s |
| GL286 | 19.29s | 544ms | 778ms | 20.61 s |
| GL287 | 19.94s | 571ms | 745ms | 21.26 s |
| GL288 | 21.97s | 576ms | 810ms | 23.36 s |
| GL289 | 22.48s | 897ms | 749ms | 24.13 s |
| GL290 | 22.60s | 584ms | 765ms | 23.95 s |
| GL291 | 19.44s | 543ms | 757ms | 20.74 s |
| GL292 | 21.24s | 881ms | 810ms | 22.93 s |
| GL293 | 21.26s | 566ms | 774ms | 22.60 s |
| GL294 | 21.04s | 542ms | 709ms | 22.29 s |
| GL295 | 52.53s | 1.51s | 1.36s | 55.40 s |
| GL296 | 3m 0.9s | 3.95s | 4.08s | 3.15 min |
| GL297 | 1m 23.5s | 2.34s | 2.15s | 87.99 s |
| GL298 | 1m 20.0s | 2.21s | 2.08s | 84.29 s |
| GL299 | 3m 2.0s | 2.86s | 2.54s | 3.12 min |
| GL300 | 1m 19.9s | 2.08s | 2.05s | 84.03 s |
| GL301 | 1m 57.5s | 2.07s | 2.30s | 2.03 min |
| GL302 | 3.43s | 0ms | 0ms | 3.43 s |
| GL303 | 1m 54.6s | 1.95s | 1.97s | 118.52 s |
| GL304 | 1m 54.0s | 2.12s | 1.97s | 118.09 s |
| GL305 | 3.35s | 0ms | 0ms | 3.35 s |
| GL306 | 3.30s | 0ms | 0ms | 3.30 s |
| GL307 | 1m 57.8s | 2.33s | 2.43s | 2.04 min |
| GL308 | 1m 55.5s | 2.26s | 2.41s | 2.00 min |
| GL309 | 2m 6.0s | 2.78s | 2.09s | 2.18 min |
| GL310 | 3.82s | 0ms | 0ms | 3.82 s |
| GL311 | 51.49s | 1.63s | 1.09s | 54.21 s |
| GL312 | 51.09s | 1.02s | 1.09s | 53.20 s |
| GL313 | 1m 3.7s | 1.50s | 1.34s | 66.54 s |
| GL314 | 327ms | 0ms | 0ms | 327 ms |
| GL315 | 285ms | 0ms | 0ms | 285 ms |
| GL316 | 286ms | 0ms | 0ms | 286 ms |
| GL317 | 1m 3.3s | 1.24s | 1.30s | 65.84 s |
| GL318 | 1m 4.8s | 1.87s | 1.28s | 67.95 s |
| GL319 | 50.41s | 1.02s | 1.08s | 52.51 s |
| GL320 | 289ms | 0ms | 0ms | 289 ms |
| GL321 | 1m 3.6s | 1.55s | 1.28s | 66.43 s |
| GL322 | 56.87s | 1.19s | 1.10s | 59.16 s |
| GL323 | 55.71s | 1.13s | 1.14s | 57.98 s |
| GL324 | 54.53s | 1.20s | 1.18s | 56.91 s |
| GL325 | 49.67s | 1.14s | 1.20s | 52.01 s |
| GL326 | 48.16s | 1.04s | 1.23s | 50.43 s |
| GL327 | 1m 0.1s | 1.38s | 1.37s | 62.85 s |
| GL328 | 274ms | 0ms | 0ms | 274 ms |
| GL329 | 273ms | 0ms | 0ms | 273 ms |
| GL330 | 271ms | 0ms | 0ms | 271 ms |
| GL331 | 1m 2.0s | 1.25s | 1.26s | 64.51 s |
| GL332 | 1m 2.3s | 1.38s | 1.23s | 64.91 s |
| GL333 | 50.36s | 1.04s | 1.17s | 52.57 s |
| GL334 | 295ms | 0ms | 0ms | 295 ms |
| GL335 | 1m 5.2s | 1.24s | 1.40s | 67.84 s |
| GL336 | 55.67s | 1.20s | 1.08s | 57.95 s |
| GL337 | 54.29s | 1.09s | 1.12s | 56.50 s |
| GL338 | 53.39s | 1.16s | 1.11s | 55.66 s |

## Follow-Up

- Re-run the checked-in card with `n_start=20` and `--batch-size 1` to see whether the integration memory blow-up was primarily batch-size driven.
- Use the saved DOT files to build one graph at a time and establish per-graph runtime/memory behavior before returning to the full all-graph Monte Carlo setup.
- If the batch-size-one all-graph run still allocates hundreds of GiB before the first iteration, inspect the integration-state layout for graph/orientation/LMB Monte Carlo to identify what scales with all 339 graph evaluators.

## Grouped Rescaling Follow-Up

This follow-up used the release binary `target/release/gammaloop` built with

```bash
nix develop --no-write-lock-file -c cargo build -p gammaloop-api --bin gammaloop --features ufo_support --release
```

The grouped card is `examples/cli/aa_aa/3L/aa_aa_grouped.toml`. The default symbolic scalar-rescaling comparison remained impractical, so the successful one-time feyngen pass used the same scalar-rescaling grouping mode but with a single fully numerical numerator-comparison sample:

```text
--numerator-grouping group_identical_graphs_up_to_scalar_rescaling
--number-of-samples-for-numerator-comparisons 1
--fully-numerical-substitution-when-comparing-numerators=true
--compare-canonized-numerator=false
```

### Grouped Diagram Generation

The grouped feyngen pass completed and saved reusable DOT/state artifacts:

| Quantity | Value |
|---|---:|
| Symbolica generation | 6192 graphs |
| After vetoed topologies | 4284 graphs |
| After complete-graph filters | 1296 graphs |
| After external-state symmetrization | 339 graphs |
| After canonization | 339 graphs |
| After numerator-aware grouping | 155 graphs |
| Isomorphically unique graphs | 163 |
| Color zeros | 58 |
| Grouped by scalar rescaling | 118 |
| Cancellations | 8 |
| Wall time for grouped feyngen stage | about 7 min |
| Peak observed RSS during numeric grouping | about 37 GiB |
| DOT export | `examples/cli/aa_aa/3L/graphs_grouped` |
| Saved state | `examples/cli/aa_aa/3L/gammaloop_state_grouped` |

The earlier release-build attempt with the default symbolic rescaling comparison remained useful as a negative control: it was stopped after about 80 min and about 166 GiB peak RSS without producing grouped DOTs/state.

### Grouped Integrand Generation

Command:

```bash
target/release/gammaloop --override-state \
  -s examples/cli/aa_aa/3L/gammaloop_state_grouped \
  generate existing -p aa_aa -i 3L
```

Observed grouped integrand generation:

| Quantity | Value |
|---|---:|
| Graph groups | 155 |
| Wall time | 11m36s |
| Peak RAM reported by GammaLoop | 9.15 GiB |
| Generation cores | 10 |
| Saved integrand payload | 153.15 MiB |
| Slowest generated graphs | GL299, GL296, GL191, GL226, GL215 |

### All-Grouped-Amplitude Smoke

The grouped state was then integrated with kinematics A, helicity `+-+-`, `SingleParametric`, Monte Carlo over graphs/orientations/LMBs, and inverse-Jacobian LMB channel weights. The 20-sample all-graph smoke reached the integration loop over 155 graphs and 9 continuous dimensions, but it was stopped before the first iteration completed because RSS kept increasing:

| Quantity | Value |
|---|---:|
| Workspace | `/tmp/aa_aa_3l_grouped_overall_smoke_mc_workspace` |
| Requested samples | 20 |
| Batch size | 1 |
| RSS after setup/sampling start | about 52 GiB |
| RSS when stopped | about 400 GiB |
| Completed iterations | 0 |

Grouping therefore removes the previous 339-graph all-amplitude setup, but the combined grouped evaluator still aggregates too much runtime state to be useful in the current all-at-once `SingleParametric` workflow.

### Per-Diagram Grouped Smoke Scan

A separate single-DOT runner was added at `examples/cli/aa_aa/3L/run_grouped_diagram_smokes.py`. It imports one grouped DOT at a time, generates the `SingleParametric` integrand, and runs a 20-sample integration with kinematics A and helicity `+-+-`. The compact per-graph table is stored in `examples/cli/aa_aa/3L/grouped_diagram_smokes_summary.csv`; phase-probe runs are stored in `examples/cli/aa_aa/3L/grouped_phase_probe_summary.csv`.

All 155 grouped diagrams completed and emitted integration results.

| Quantity | Min | Median | Mean | Max |
|---|---:|---:|---:|---:|
| End-to-end wall time per graph | 15.79 s | 61.26 s | 66.58 s | 169.54 s |
| Peak RSS per graph | 2.88 GiB | 17.79 GiB | 18.11 GiB | 46.55 GiB |
| Avg total runtime per sample | 0.176 ms | 0.402 ms | 0.447 ms | 1.627 ms |
| Avg integrand runtime per sample | 0.102 ms | 0.267 ms | 0.314 ms | 1.417 ms |
| NAN or unstable samples | 0.0% | 0.0% | 0.065% | 5.0% |

Top peak-RSS graphs:

| Graph | Wall s | RSS GiB | Avg ms/sample | Unstable % |
|---|---:|---:|---:|---:|
| GL309 | 149.6 | 46.6 | 0.638 | 0.0 |
| GL301 | 163.7 | 46.5 | 1.627 | 0.0 |
| GL303 | 140.0 | 46.5 | 0.521 | 0.0 |
| GL307 | 135.3 | 46.5 | 0.490 | 0.0 |
| GL244 | 126.7 | 35.9 | 1.605 | 0.0 |

Top runtime-per-sample graphs:

| Graph | Wall s | RSS GiB | Avg ms/sample | Unstable % |
|---|---:|---:|---:|---:|
| GL301 | 163.7 | 46.5 | 1.627 | 0.0 |
| GL244 | 126.7 | 35.9 | 1.605 | 0.0 |
| GL296 | 169.5 | 25.3 | 1.406 | 5.0 |
| GL191 | 137.3 | 24.5 | 1.278 | 0.0 |
| GL000 | 32.9 | 9.5 | 1.204 | 0.0 |

Only `GL296` and `GL299` reported nonzero instability in this 20-sample smoke scan, each with one unstable/NAN-class evaluation out of 20 samples. Both still emitted integration results.

### Real vs Imaginary Training Probe

The user pointed out that the three-loop setup may need real-part rather than imaginary-part training. I ran separate high-stat single-graph probes using the same runner with explicit `integrator.integrated_phase` choices.

The first probes (`GL026`, `GL001`, `GL047`) were not reliable phase discriminators: at higher statistics, rare-weight variance dominated and neither real nor imaginary converged cleanly. A lower-max-weight graph, `GL105`, was pushed to 10M samples with real-part training:

| Probe | Graph | Samples | Re | Re err | Im | Im err |
|---|---|---:|---:|---:|---:|---:|
| real-trained | GL105 | 10,000,000 | -5.31e-7 | 1.19e-6 | 6.84e-7 | 1.37e-6 |

This probe is compatible with zero in both components. I therefore did not restart the full per-diagram sweep with real-part training. The completed smoke scan should be read primarily as a runtime/memory/stability scan; its 20-sample central values are too noisy for a phase conclusion. Both real and imaginary columns are retained in the CSV summaries for follow-up.
