# GL638 final code review notes

- Threshold subtraction remains enabled and metadata is supplied by `GL638.dot`;
  no `force_cuts` path is used. The process specification selects the relevant
  `e+ e- > t t~ h` graph content at import.
- All six cut channels remain in the sampling catalogue when a joint channel is
  present. A joint map is a proposal and never changes the physical cut sum.
- The soft channel is a full-volume LMB and the only channel with an explicit
  OSE override. OSE is a positive partition score, not a replacement for its
  map Jacobian. Surface and composed maps retain exact density/support checks.
- OSE uses master-frame edge energies, current master masses, and the affine
  inverse. `E_cm^(-3L)` gives the common dimension for this fixed loop order.
- The underflow correction is component-local: a completed weighted value below
  binary64's normal range may be rounded to zero. A remaining factor that would
  promote it into the normal range still raises an error. Overflow, nonfinite
  values and invalid stability operations remain hard failures.
- The accepted c=1 SUM control was `Re=(-2.866±3.029)e-5 pb` and
  `|Re|=(2.632±0.303)e-4 pb`; the c=4 control is provisional because of one
  unresolved sample. These controls do not claim the new OSE result.

Before the 100-worker run, rebuild the native client from current source and
compare per-point timing with the prior 44 ms at 50 workers. If 100-worker
throughput per core collapses, report that scaling anomaly rather than hiding it.
