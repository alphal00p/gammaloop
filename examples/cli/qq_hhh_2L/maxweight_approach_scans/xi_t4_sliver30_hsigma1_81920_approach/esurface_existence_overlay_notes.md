# E-surface existence overlay notes

The black/red traces in `highstat_iter_0001_imag_plus_approaches_with_existence.png`
and `highstat_iter_0001_imag_minus_approaches_with_existence.png` use this
event-level convention:

```text
1    if the selected threshold_counterterm key is emitted with nonzero weighted
     absolute value in the additional-weight decomposition;
1e-6 otherwise.
```

This is not a direct query of the internal SOCP existence predicate. It tests
whether the candidate local threshold CT is present in the evaluated event
decomposition at the scan point.

The two selected pairs are those dominating the suspicious high-stat fixed-bin
approach maxima:

```text
im[+], graph 0, LMB channel 14:
  black = GL06_isr_p2_ct threshold_counterterm:36:0, edges [9,14,16]
  red   = GL06_isr_p2_ct threshold_counterterm:13:0, edges [11,15]

im[-], graph 0, LMB channel 9:
  black = GL05_isr_p2_ct threshold_counterterm:4:0, edges [9,12,15]
  red   = GL05_isr_p1_ct threshold_counterterm:34:0, edges [15,16]
```

The overlay scan shows both selected CTs emitted across the approach lines,
except for one far-left center-diagonal im[+] point where both selected
`GL06_isr_p2_ct` CTs are absent. That missing point is away from the central
ridge/peak. The shoulder growth seen in the plots therefore is not explained
by one of these two event-level CTs switching off at the shoulder.

Equations are written as

```text
sum_i OSE(edge_i) + DeltaE_T = 0
```

where `DeltaE_T` is the external-energy shift of the corresponding threshold
surface. The Python API currently exposes the OSE edge set and the graph/LMB
edge signatures, but not the raw `external_shift` list of the `Esurface`.

External signature order is `[Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8]`, with the exact-xi
helper slots `Q2=Q3=Q4=Q5=xi`.

## im[+] pair

For `GL06_isr_p2_ct` in the `[9,12]` loop basis:

```text
T36: OSE(9) + OSE(14) + OSE(16) + DeltaE_GL06_isr_p2_ct,T36 = 0
  edge  9: loop [1, 0],  ext [ 0,  0,  0,  0,  0,  0,  0, 0, 0]
  edge 14: loop [0,-1],  ext [ 0, -1,  0, -1,  1,  1,  0, 0, 0]
  edge 16: loop [1, 1],  ext [ 1,  1,  1,  1, -1, -1, -1, 0, 0]

T13: OSE(11) + OSE(15) + DeltaE_GL06_isr_p2_ct,T13 = 0
  edge 11: loop [1, 1],  ext [ 1,  1,  1,  1, -1, -1,  0, 0, 0]
  edge 15: loop [1, 1],  ext [ 0,  0,  0,  0,  0,  0,  0, 1, 0]
```

## im[-] pair

For `GL05_isr_p2_ct` in the `[11,14]` loop basis:

```text
T4: OSE(9) + OSE(12) + OSE(15) + DeltaE_GL05_isr_p2_ct,T4 = 0
  edge  9: loop [1,-1],  ext [ 1,  0,  1,  0,  0,  0,  0,  0, 0]
  edge 12: loop [0,-1],  ext [ 0, -1,  0, -1,  1,  1,  0,  0, 0]
  edge 15: loop [1, 0],  ext [ 1,  1,  1,  1, -1, -1,  0, -1, 0]
```

For `GL05_isr_p1_ct` in the `[11,14]` loop basis:

```text
T34: OSE(15) + OSE(16) + DeltaE_GL05_isr_p1_ct,T34 = 0
  edge 15: loop [1,0],  ext [1, 1, 1, 1, -1, -1, 0, -1, 0]
  edge 16: loop [1,0],  ext [0, 0, 0, 0,  0,  0, 1,  0, 0]
```
