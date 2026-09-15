# GL638 big-run setup

This card is the reproducible 600 GeV production setup for the next GL638 run.
It keeps all 936 production orientations, six Cutkosky cuts, nineteen threshold
variants, 3D local UV counterterms and integrated UV counterterms. The incoming
momenta are 300 GeV beams; `mu_r=91.188 GeV` and `m_uv=50 GeV`.

Sampling uses six LU-h cut channels, the composed Cut-1/H-Z joint channel, the
Cut-3 right-threshold channel, and the explicit soft LMB `[6,12,13,14]`.
Channels are explicitly summed so real and absolute observables share one
covariance-preserving point. Advanced maps use exact map-density scores; the
soft full-volume LMB uses OSE with `alpha=3`. Maps use spherical coordinates,
linear radial mapping, `b=0.3`, and `power=2`.

The integrator uses five adaptive iterations of 16,384 outer samples (81,920
total), 16 bins, minimum 64 samples per update, learning rates 0.25, variance
training, and seed 62001. The run uses 100 workers, SymJIT O3 compression, real
phase integration, and per-iteration maximum diagnostics. Stability uses
separate real/imaginary components and Double/Quad/Arb escalation, with final
Arb `ecm_relative_tolerance_for_re=1e-100` at integrated dimension `-2`.
This allowance is relative to `E_cm` and does not waive nonfinite or meaningful
real discrepancies.

The K6 Gaussian acceptance passed for inverse-density and OSE scores in MC and
SUM modes. The first physical OSE pilot stopped on a reporting-boundary tail
underflow below binary64's normal range; the code now rounds such suppressed
completed contributions to zero while still rejecting a factor that promotes one
into the normal range. Its focused regression passes. Rebuild and rerun the
pilot before interpreting physical central values.
