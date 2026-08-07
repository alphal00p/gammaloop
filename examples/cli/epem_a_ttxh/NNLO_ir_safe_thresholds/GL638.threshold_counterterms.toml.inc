# Verified current GL638 cure document, matching the dressed graph attribute. The shared one-loop
# association and complementary squared-E-surface ratios cure the selected-orientation correlated
# q13 soft/threshold limit of cut (2,4,12); their sufficiency for every other cut has not yet been
# demonstrated.
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
      expression = "eta(star, eset(7, 8))^2 / (eta(star, eset(7, 8))^2 + eta(star, eset(8, 12, 14))^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(star, eset(8, 12, 14))^2 / (eta(star, eset(7, 8))^2 + eta(star, eset(8, 12, 14))^2)"

[[cuts]]
edges = [2, 6, 10]

  # The legacy maximal subspace is [3,7] only for this association. Solving it in the shared [7]
  # frame supplies the missing LU partner of the same threshold under cuts (2,6,12,13) and
  # (2,4,12), removing the q13-soft inverse powers.
  [[cuts.thresholds]]
  edges = [8, 12, 14]

    [[cuts.thresholds.counterterms]]
    name = "shared_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(star, eset(7, 8))^2 / (eta(star, eset(7, 8))^2 + eta(star, eset(8, 12, 14))^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(star, eset(8, 12, 14))^2 / (eta(star, eset(7, 8))^2 + eta(star, eset(8, 12, 14))^2)"

[[cuts]]
edges = [2, 4, 12]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(star, eset(7, 8))^2 / (eta(star, eset(7, 8))^2 + eta(star, eset(8, 12, 14))^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(star, eset(8, 12, 14))^2 / (eta(star, eset(7, 8))^2 + eta(star, eset(8, 12, 14))^2)"

[[cuts]]
edges = [2, 4, 10, 13]

  [[cuts.thresholds]]
  edges = [7, 8]

    [[cuts.thresholds.counterterms]]
    name = "intrinsic_1l"
    subspace = [7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(star, eset(7, 8))^2 / (eta(star, eset(7, 8))^2 + eta(star, eset(2, 4, 10, 13))^2)"

    [[cuts.thresholds.counterterms]]
    name = "embedded_2l"
    subspace = [3, 7]
    parent_lmb = [3, 4, 7, 10]

      [cuts.thresholds.counterterms.multiplier]
      expression = "eta(star, eset(2, 4, 10, 13))^2 / (eta(star, eset(7, 8))^2 + eta(star, eset(2, 4, 10, 13))^2)"
