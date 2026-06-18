# axis_x2 CT contribution monitor

Runtime: xi = 4*(p1+p2), sliver_width = 30, gaussian_width = 1, threshold h_sigma = 1.
Fixed bin: target `highstat_iter_0001_imag_minus`, graph `0`, LMB channel `9`.
Reference root from the rstar probe: `x2_root = 9.4411812352007729e-01`.

## Main observation

The sharp feature is not caused by one of the dominant CTs turning on or off. Across the monitored points, `0` term keys are absent from at least one point; the large terms listed below are present on both sides of the root. They all flip sign when `x2 - x2_root` changes sign, and they then decay together as `1 / |x2 - x2_root|`.

Dominant terms:

- `GL05_isr_p2_ct` `threshold_counterterm:4:0` edges `[9, 12, 15]`
- `GL05_isr_p1_ct` `threshold_counterterm:34:0` edges `[15, 16]`
- `GL05` `threshold_counterterm:31:0` edges `[9, 10, 15]`
- `GL05` `threshold_counterterm:5:0` edges `[9, 12, 16]`

## Nearest root probes

- nearest left: `dx2 = -1.00000008e-10`, total `+6.02712742e+03 +1.54236277e+04i |w|=1.65594250e+04`, dominant CT sum `+6.02712748e+03 +1.54236276e+04i |w|=1.65594249e+04`
- nearest right: `dx2 = +1.00000008e-10`, total `-5.95996439e+03 -1.52517548e+04i |w|=1.63748954e+04`, dominant CT sum `-5.95996433e+03 -1.52517548e+04i |w|=1.63748954e+04`

## Saved scan rows around the feature

- row `30` dx2 `-1.99098135e-06` scan `+3.00967451e-01 +7.70384053e-01i`, event total `+3.00967451e-01 +7.70384053e-01i |w|=8.27087054e-01`, dominant CT sum `+3.01032817e-01 +7.70326911e-01i |w|=8.27057620e-01`, other `-6.53653036e-05 +5.71416050e-05i |w|=8.68204235e-05`
- row `31` dx2 `-9.85254806e-07` scan `+6.08246394e-01 +1.55672149e+00i`, event total `+6.08246394e-01 +1.55672149e+00i |w|=1.67133045e+00`, dominant CT sum `+6.08311759e-01 +1.55666435e+00i |w|=1.67130102e+00`, other `-6.53648358e-05 +5.71404971e-05i |w|=8.68193421e-05`
- row `32` dx2 `-3.26374626e-07` scan `+1.83628143e+00 +4.69930561e+00i`, event total `+1.83628143e+00 +4.69930561e+00i |w|=5.04533474e+00`, dominant CT sum `+1.83634679e+00 +4.69924847e+00i |w|=5.04530531e+00`, other `-6.53645294e-05 +5.71397713e-05i |w|=8.68186337e-05`
- row `33` dx2 `+1.05276595e-07` scan `-5.69296409e+00 -1.45682940e+01i`, event total `-5.69296409e+00 -1.45682940e+01i |w|=1.56411327e+01`, dominant CT sum `-5.69289873e+00 -1.45683512e+01i |w|=1.56411621e+01`, other `-6.53643286e-05 +5.71392958e-05i |w|=8.68181696e-05`
- row `34` dx2 `+3.88063637e-07` scan `-1.54447404e+00 -3.95216453e+00i`, event total `-1.54447404e+00 -3.95216453e+00i |w|=4.24323043e+00`, dominant CT sum `-1.54440867e+00 -3.95222167e+00i |w|=4.24325986e+00`, other `-6.53641971e-05 +5.71389843e-05i |w|=8.68178656e-05`
- row `35` dx2 `+5.73325504e-07` scan `-1.04541830e+00 -2.67506366e+00i`, event total `-1.04541830e+00 -2.67506366e+00i |w|=2.87208374e+00`, dominant CT sum `-1.04535294e+00 -2.67512080e+00i |w|=2.87211317e+00`, other `-6.53641109e-05 +5.71387802e-05i |w|=8.68176664e-05`
- row `36` dx2 `+6.94695843e-07` scan `-8.62783268e-01 -2.20769428e+00i`, event total `-8.62783268e-01 -2.20769428e+00i |w|=2.37029724e+00`, dominant CT sum `-8.62717903e-01 -2.20775142e+00i |w|=2.37032667e+00`, other `-6.53640545e-05 +5.71386465e-05i |w|=8.68175359e-05`
- row `37` dx2 `+7.74209018e-07` scan `-7.74179175e-01 -1.98095335e+00i`, event total `-7.74179175e-01 -1.98095335e+00i |w|=2.12685908e+00`, dominant CT sum `-7.74113811e-01 -1.98101048e+00i |w|=2.12688851e+00`, other `-6.53640175e-05 +5.71385589e-05i |w|=8.68174504e-05`
- row `38` dx2 `+8.26300370e-07` scan `-7.25377232e-01 -1.85606749e+00i`, event total `-7.25377232e-01 -1.85606749e+00i |w|=1.99277662e+00`, dominant CT sum `-7.25311868e-01 -1.85612463e+00i |w|=1.99280604e+00`, other `-6.53639933e-05 +5.71385015e-05i |w|=8.68173944e-05`
- row `39` dx2 `+8.60426902e-07` scan `-6.96609347e-01 -1.78244947e+00i`, event total `-6.96609347e-01 -1.78244947e+00i |w|=1.91373737e+00`, dominant CT sum `-6.96543983e-01 -1.78250661e+00i |w|=1.91376679e+00`, other `-6.53639774e-05 +5.71384639e-05i |w|=8.68173577e-05`
- row `40` dx2 `+8.82784168e-07` scan `-6.78968563e-01 -1.73730610e+00i`, event total `-6.78968563e-01 -1.73730610e+00i |w|=1.86526963e+00`, dominant CT sum `-6.78903199e-01 -1.73736324e+00i |w|=1.86529906e+00`, other `-6.53639670e-05 +5.71384393e-05i |w|=8.68173337e-05`
- row `41` dx2 `+8.97431051e-07` scan `-6.67888116e-01 -1.70895085e+00i`, event total `-6.67888116e-01 -1.70895085e+00i |w|=1.83482630e+00`, dominant CT sum `-6.67822752e-01 -1.70900799e+00i |w|=1.83485573e+00`, other `-6.53639602e-05 +5.71384232e-05i |w|=8.68173180e-05`
- row `42` dx2 `+9.07026642e-07` scan `-6.60823015e-01 -1.69087101e+00i`, event total `-6.60823015e-01 -1.69087101e+00i |w|=1.81541506e+00`, dominant CT sum `-6.60757651e-01 -1.69092815e+00i |w|=1.81544449e+00`, other `-6.53639558e-05 +5.71384126e-05i |w|=8.68173076e-05`
- row `43` dx2 `+9.13312988e-07` scan `-6.56274953e-01 -1.67923237e+00i`, event total `-6.56274953e-01 -1.67923237e+00i |w|=1.80291934e+00`, dominant CT sum `-6.56209589e-01 -1.67928950e+00i |w|=1.80294877e+00`, other `-6.53639528e-05 +5.71384057e-05i |w|=8.68173009e-05`
- row `44` dx2 `+9.17431353e-07` scan `-6.53329180e-01 -1.67169403e+00i`, event total `-6.53329180e-01 -1.67169403e+00i |w|=1.79482588e+00`, dominant CT sum `-6.53263816e-01 -1.67175117e+00i |w|=1.79485531e+00`, other `-6.53639509e-05 +5.71384011e-05i |w|=8.68172965e-05`
- row `t0` dx2 `+9.25254727e-07` scan `-6.47805519e-01 -1.65755879e+00i`, event total `-6.47805519e-01 -1.65755879e+00i |w|=1.77964972e+00`, dominant CT sum `-6.47740155e-01 -1.65761593e+00i |w|=1.77967915e+00`, other `-6.53639473e-05 +5.71383925e-05i |w|=8.68172881e-05`

## Pole scaling check

- `GL05_isr_p2_ct` `threshold_counterterm:4:0` left mean `(x2-x2_root)*w = -2.96013093e-07 -7.02671957e-07i |w|=7.62477298e-07` from `17` points
- `GL05_isr_p2_ct` `threshold_counterterm:4:0` right mean `(x2-x2_root)*w = -2.96005024e-07 -7.02660499e-07i |w|=7.62463607e-07` from `28` points
- `GL05_isr_p1_ct` `threshold_counterterm:34:0` left mean `(x2-x2_root)*w = -2.80625310e-07 -5.74015920e-07i |w|=6.38940405e-07` from `17` points
- `GL05_isr_p1_ct` `threshold_counterterm:34:0` right mean `(x2-x2_root)*w = -2.80617618e-07 -5.74007027e-07i |w|=6.38929037e-07` from `28` points
- `GL05` `threshold_counterterm:31:0` left mean `(x2-x2_root)*w = -7.99941586e-09 -1.34382804e-07i |w|=1.34620685e-07` from `17` points
- `GL05` `threshold_counterterm:31:0` right mean `(x2-x2_root)*w = -7.99877344e-09 -1.34380780e-07i |w|=1.34618626e-07` from `28` points
- `GL05` `threshold_counterterm:5:0` left mean `(x2-x2_root)*w = -1.47042022e-08 -1.22658895e-07i |w|=1.23537112e-07` from `17` points
- `GL05` `threshold_counterterm:5:0` right mean `(x2-x2_root)*w = -1.47034102e-08 -1.22657028e-07i |w|=1.23535164e-07` from `28` points

## Outputs

- point summary TSV: `/Users/vjhirsch/Documents/Work/gammaloop_lucien/examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor_points.tsv`
- term TSV: `/Users/vjhirsch/Documents/Work/gammaloop_lucien/examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor_terms.tsv`
- JSON: `/Users/vjhirsch/Documents/Work/gammaloop_lucien/examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor.json`
- component plot: `/Users/vjhirsch/Documents/Work/gammaloop_lucien/examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor_components.png`
- magnitude plot: `/Users/vjhirsch/Documents/Work/gammaloop_lucien/examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor_abs.png`
