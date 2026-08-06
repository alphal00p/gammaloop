# Verified current GL638 smoke document, matching the dressed graph attribute. These deliberately
# illustrative complementary squared-E-surface ratios are not the final physics prescription.
# Each named variant remains one CT in one exact solving subspace and is presented to overlap
# finding as its own projected E-surface instance.
schema_version = 1

[[cuts]]
edges = [2, 6, 12, 13]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, 8, 12, 14)^2 / (eta(effective, 7, 8)^2 + eta(effective, 8, 12, 14)^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, 7, 8)^2 / (eta(effective, 7, 8)^2 + eta(effective, 8, 12, 14)^2)"

[[cuts]]
edges = [2, 6, 10]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, 8, 12, 14)^2 / (eta(effective, 7, 8)^2 + eta(effective, 8, 12, 14)^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, 7, 8)^2 / (eta(effective, 7, 8)^2 + eta(effective, 8, 12, 14)^2)"

[[cuts]]
edges = [2, 4, 12]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, 8, 12, 14)^2 / (eta(effective, 7, 8)^2 + eta(effective, 8, 12, 14)^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, 7, 8)^2 / (eta(effective, 7, 8)^2 + eta(effective, 8, 12, 14)^2)"

[[cuts]]
edges = [2, 4, 10, 13]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, 2, 4, 10, 13)^2 / (eta(effective, 7, 8)^2 + eta(effective, 2, 4, 10, 13)^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(effective, 7, 8)^2 / (eta(effective, 7, 8)^2 + eta(effective, 2, 4, 10, 13)^2)"
