# `axis_x2` probe of `r - rstar`

This probe fixes graph `0`, LMB channel `9`, all scan coordinates except `x2`,
and dumps the `GL05_isr_p2_ct:T4` threshold parameters.

In `spherical_common_radial`, `x2` is the azimuthal angle of the first loop
three-momentum:

```text
phi0 = 2*pi*x2
```

The dense probe shows that `rstar(x2)` is smooth and essentially linear across
the apparent drop:

```text
d(r - rstar)/dx2 = +2.2745738365e+02
root x2 where r - rstar = 0: 0.9441181235200773
max absolute residual to a linear fit: 3.82e-10
```

The corresponding implicit E-surface derivative is not small:

```text
mean dE/dr = mean eta' = 2.3502743593
implied dE/dx2 = 5.3458725663e+02
```

So the change is not a branch switch or discontinuity. It is a smooth crossing of
the selected threshold root. The apparent abruptness comes from the tiny
denominator gap:

```text
x2 = 0.9441182287966725: r - rstar = +2.394608691e-05
x2 = 0.9441185115837143: r - rstar = +8.826804356e-05
```

The absolute change in `rstar` is only `6.432e-05` on top of
`rstar ~= 1.061860e+03`, but because the local CT contains `1/(r-rstar)`, this
small absolute change gives the observed factor:

```text
2.394608691e-05 / 8.826804356e-05 = 0.271288
```

The probe data are in `rstar_x2_probe.tsv` and `rstar_x2_probe.json`.
