# GL134 phase controls, 17 September 2026

This archive preserves the two-point physical-momentum comparison of both saved
sampling states with fresh `num=i` and `num=1` generations. See
[the diagnosis](../GL134_LEGACY_PHASE_NORMALIZATION_20260917.md) and
[the result matrix](point_matrix.json). Each fresh generation used the full
unfiltered catalogue of **686 CFF orientations**, with an explicit orientation
sum and the sole physical cut. The older manual report's 64-orientation count
does not describe these current generations.

The four mode directories contain the exact executed cards or command blocks,
argv records, original inspect JSON, compressed console logs, and exported DOTs.
Paths in those historical commands remain the paths actually used; they are
provenance, not portable commands to rerun against subsequently migrated inputs.
For replay, use a new diagnostic state/output directory, import the archived
`GL134_num_i.dot` or `GL134_num1.dot`, and substitute the desired runtime-settings
path in the archived fresh-generation card. Invoke its `probe` block with a
compatible checkout development environment and Symbolica license. No integration
is required. Do not use `run all`, which also contains an integration block.

`GL134_num_i.dot` restores only the recorded original `num=i` text in the frozen
`GL134_num1.dot`; its SHA-256 was checked against `dot_control.json`. All three
`num=i` exported graphs are byte-identical. Replacing only their `num` attribute
with one makes them byte-identical to the fresh `num=1` export. The empty
`fresh_export_diff.txt` records that comparison.

The old saved states were unchanged by these controls; their input hashes are in
the saved-mode invocation records. Large binary states and credential-bearing
launch scripts are deliberately excluded. `manifest.json` hashes the copied
evidence. The fresh controls' inherited Im training setting is historical and
does not affect these momentum-space inspections; the migrated production cards
now select Re.
