# axis_x2 im[-] sharp-drop cause

This note compares the two fixed-bin points bracketing the sharp drop on the
`highstat_iter_0001_imag_minus / axis_x2` approach scan:

```text
local_peak_before_drop:
  signed_side_fraction = -8.685113737513521e-7
  x2 = 0.9441182287966725
  fixed-bin weight = -5.69296409 - 14.56829402 i
  |w| = 15.64113266

right_after_drop:
  signed_side_fraction = -5.689866029018304e-7
  x2 = 0.9441185115837143
  fixed-bin weight = -1.54447404 - 3.95216453 i
  |w| = 4.24323043
```

## What does not change

No threshold counterterm appears or disappears across the drop. The two
suspected dual-cancelling surfaces remain present on both sides:

```text
GL05_isr_p2_ct:T4   edges [9,12,15]  exists = 1
GL05_isr_p1_ct:T34  edges [15,16]    exists = 1
```

The raw symbolic CT evaluators are also smooth across the drop. For example:

```text
GL05_isr_p2_ct:T4 raw_result:
  before = -2.641030269164389e-4 + 4.477115270636740e-3 i
  after  = -2.641084037831862e-4 + 4.477109684226697e-3 i

GL05_isr_p1_ct:T34 raw_result:
  before = -1.095239287840850e-5 + 3.758237432774540e-3 i
  after  = -1.095826875315196e-5 + 3.758233242960828e-3 i
```

The event multiplicative factor is the same on both points:

```text
full_multiplicative_factor = 5.42290327520015e19
```

## What changes

The dominant CTs all inherit the same threshold-helper pole factor. The
selected radial distance to the threshold changes from

```text
before: r_left - rstar_left = 2.394608736722148e-5
after:  r_left - rstar_left = 8.826804378259112e-5
```

The ratio is

```text
(r_left - rstar_left)_before / (r_left - rstar_left)_after
  = 0.2712882980
```

which matches the drop of all large CT contributions:

```text
GL05_isr_p2_ct:T4   |after|/|before| = 0.27128797
GL05_isr_p1_ct:T34  |after|/|before| = 0.27128799
GL05:T31            |after|/|before| = 0.27128825
GL05:T5             |after|/|before| = 0.27128820
```

The other threshold-helper inputs are essentially unchanged:

```text
h_left_th:
  before = 1.0626436737215202e-3
  after  = 1.0626437380909040e-3

deta_left_th:
  before = 2.3502743401656603
  after  = 2.3502744962099950

damp_plus_left:
  before = 0.9999999999999994
  after  = 0.9999999999999922
```

## Mechanism

For this first-order two-loop threshold CT, the local helper contains

```text
J(r) * [ damp_plus_left / (r_left - rstar_left)
       + damp_minus_left / (-r_left - rstar_left) ]
```

with

```text
J(r) = 1 / r_left^5
```

up to the common `(2*pi)^(-6)` normalization and the raw pass-on evaluator
coefficient. In this probe, `damp_plus_left ~= 1`, `damp_minus_left ~= 0.011`,
and `r_left` is identical at the two points, so the change is controlled almost
entirely by `1 / (r_left - rstar_left)`.

Therefore the sharp drop is not caused by an existence-mask transition. It is a
continuous threshold-pole shoulder: the scan point before the drop lies about
3.7 times closer to the selected threshold solution than the point after it.

The reason this becomes visible in the complete fixed-bin weight is that the
large `T4/T34/T31/T5` CTs have the same phase and scale together there; the
nearby CTs that would need to cancel this local pole are not providing an
opposite contribution of comparable size on that side of the shoulder.
