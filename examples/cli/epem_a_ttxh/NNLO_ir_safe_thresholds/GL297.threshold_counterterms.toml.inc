# Verified current GL297 cure document, matching the dressed graph attribute. Only the four
# associations whose legacy CTs are genuinely two-loop are overridden; every other threshold
# remains an implicit legacy/default association.
schema_version = 1

[[cuts]]
edges = [2, 4, 9]

  [[cuts.thresholds]]
  edges = [5, 7]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [5]

  [[cuts.thresholds]]
  edges = [5, 11, 13]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [5]

[[cuts]]
edges = [2, 6, 7]

  [[cuts.thresholds]]
  edges = [3, 9]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [3]

  [[cuts.thresholds]]
  edges = [3, 11, 12]

    [[cuts.thresholds.counterterms]]
    name = "forced_1l"
    subspace = [3]
