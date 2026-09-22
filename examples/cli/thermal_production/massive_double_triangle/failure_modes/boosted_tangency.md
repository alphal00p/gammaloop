# Boosted radial-tangency precision failure

The stable benchmark also has a reproducible **f64 inspection failure** at a
nearly tangent LL radial surface. This is a coordinate/precision limitation,
not evidence for a physical branch point or an incorrect integrated rate.

Use `MX=1`, `MH=0.6`, `ML=0.25`, `T=0.5`, and incoming momentum
`(4.1231056256176606, 0, 0, 4)`. The saved LMB inspection point is

```
loop0 = (0.4330127018922193, 0, 0.9013878188659973)
loop1 = (0.3, 0.5, -0.2)
```

For the LL two-body cut, the first ray is analytically tangent when
`sin(theta)=sqrt(3)/4`, `cos(theta)=sqrt(13)/4`. Its two radial roots coalesce at
`sqrt(13)/8`. The serialized decimal momenta are slightly off the exact tangent;
arbitrary precision resolves two very close roots.

| Control | Result | LL events |
|---|---|---:|
| Thermal f64, tangent decimals | Root-solver error | No result |
| Vacuum f64, same point | Same root-solver error | No result |
| Thermal `inspect -f` | Finite `0.006677593676996147` | 2 |
| Vacuum `inspect -f` | Finite `0.00569178894347385` | 2 |
| Thermal, sine reduced by one part per million | Finite `4.4299703823229157e-7` | 2 |
| Thermal, sine increased by one part per million | Finite `2.309016625270007e-9` | 0 |

`-f` is `--use_arb_prec`; these controls are labelled arbitrary precision, not
quad precision. All entries were rerun from `boosted_tangency.toml`.
`boosted_results.json` preserves full successful inspect output, exact failing
error text, exit codes, and invocation commands.

The pipeline shares the model, graph, generated CFF, LU selectors, and numerical
input settings. Generation succeeds. The first failure is radial root
certification before residue/event evaluation: `all_roots_in_interval` supplies
a tiny compactified bracket near `[0.31067470987665047, 0.3106747098766789]`,
but the safeguarded solver recomputes the supposed negative endpoint as
`+3.9020534377930704e-17` and rejects it. Vacuum reproduces the failure, excluding
thermal distribution weights as its cause. Higher precision removes the error
for these serialized inputs. The exact mathematical tangent has a vanishing
radial derivative, outside the simple-root solver contract; its angular
Jacobian singularity alone does not imply a divergent physical rate.

Reproduce from the repository root, selecting one block per invocation:

```sh
./gammaloop --read-only-state -s /tmp/gammaloop-boosted-failure-state \
  -n -L off \
  examples/cli/thermal_production/massive_double_triangle/failure_modes/boosted_tangency.toml \
  run generate tangent_thermal
```

Other blocks are `tangent_vacuum`, `tangent_thermal_arb`,
`tangent_vacuum_arb`, `inside`, and `outside`. The two f64 tangent blocks exit 1
as expected; the four controls exit 0.

Additional bounded scans did **not** establish failures: a generic boosted
1000-sample thermal run had zero NaN/unstable samples, and an unstable-vertex
choice `MX=1.6, MH=0.25, ML=0.6` also completed 1000 samples with zero
NaN/unstable samples. Their relative errors were 74% and 125%, respectively,
so they provide neither accurate rates nor evidence of an incorrect value.
Their cards/settings and native results are included in `boosted_results.json`.
No core code or existing tests were modified for this scan.
