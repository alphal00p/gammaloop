# Counterterm scalar-growth reproducer

`orientation_58_muv.toml` selects one orientation of the three-loop massless
Figure B.1 graph and applies ordinary UV subtraction. The MUV prescriptions,
relevant evaluator options and compact logging are contained in the card. The
graph is the existing `tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot`.

Run from the repository root in the normal GammaLoop development environment,
using the existing `target/dev-optim/gammaloop` binary:

```bash
repro_state="$(mktemp -d /tmp/gammaloop-scalar-growth.XXXXXX)"
env -u GL_LOGFILE_FILTER -u GL_DISPLAY_FILTER -u GL_ALL_LOG_FILTER \
  RAYON_NUM_THREADS=8 RUST_MIN_STACK=134217728 \
  timeout --signal=TERM --kill-after=5s 60s \
  ./target/dev-optim/gammaloop \
  ./examples/cli/uv_scalar_growth/orientation_58_muv.toml \
  -n -s "$repro_state" -t scalar-growth run generate < /dev/null
printf 'Trace: %s/logs/gammalog-scalar-growth.jsonl\n' "$repro_state"
```

The 60-second timeout deliberately bounds generation; exit 124 means that bound
was reached. Several GiB of RAM are needed: observed process peaks were about
6.8–7.3 GiB in the bounded runs. `-n` disables state autosave; the fresh state
directory receives the requested trace log. Existing `GL_LOGFILE_FILTER`,
`GL_DISPLAY_FILTER` or `GL_ALL_LOG_FILTER` variables can override card logging;
the command above unsets them to use the card's own filters.

The behavior is reproduced once `evaluator_stack_parse_atom_done` reports a
large `result_bytes` value, approximately 394 MB. Verification of this card on
2026-09-14 reached that checkpoint in 17.79 seconds:

- The selected evaluator has 14 scalar summands.
- The largest local Taylor output is 591,473 packed bytes.
- Tensor preprocessing produces a 394,216,192-byte scalar atom in 15.12 seconds.
- Observed process peak memory is 7.27 GiB.

The verification run stopped deliberately at that checkpoint, before evaluator
construction, with a 60-second wall limit and a sampled 12-GiB RSS limit. The
command above uses the wall limit only. Small packed-size differences between
runs are expected; the earlier matched run produced 394,216,189 bytes.

Completing Symbolica evaluator construction is unnecessary for observing this
growth. Numeric summand indices can reorder between runs; compare stage sizes
and stable forest/source-map identities. This card reproduces the current
implementation's ordinary-U behavior; it makes no claim about a pre-change
revision or numerical UV/IR correctness.

The [investigation report](../../../docs/validation/soft-ct-expression-growth-2026-09-14.md)
contains the matched U/H measurements and their limits.
