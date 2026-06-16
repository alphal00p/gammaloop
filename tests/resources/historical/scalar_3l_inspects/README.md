# Historical scalar three-loop inspection values

These 198 files were retained from `raised_energy_cff_wip` at `91142139e`.
They are historical diagnostic artifacts, not active Insta snapshots or current
regression oracles. The current scalar route tests do not read these values.
Several entries are cancellation-scale or subnormal results, and their original
numerical accuracy has not been independently certified. Do not regenerate or
accept them as correctness baselines without establishing an independent oracle.

Active tests compare complete physical cuts across the direct, explicit-sum and
projected routes, with independent source-zero certificates where applicable.
Moving these files out of the test snapshot directory makes the coverage boundary
explicit while preserving the original data and frontmatter.
