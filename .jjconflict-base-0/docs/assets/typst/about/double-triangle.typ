#import "../portal-graphs/figure.typ": render

// Start from the scalar double-triangle topology used by GammaLoop's test
// suite, then attach website-only particle styles and fixed presentation
// coordinates. The topology remains owned by the test resource; Linnest parses,
// lays out, and draws the complete illustration.
#let input = (
  read("../../../../tests/resources/graphs/double_triangle.dot")
  .replace("v1;", "v1 [label=\"\", hidden=true, pos=\"-8,0!\"];")
  .replace("v2;", "v2 [label=\"\", hidden=true, pos=\"8,0!\"];")
  .replace("\nv3 [", "\nv3 [label=\"\", pos=\"-4,0!\", ")
  .replace("\nv4 [", "\nv4 [label=\"\", pos=\"0,2.2!\", ")
  .replace("\nv5 [", "\nv5 [label=\"\", pos=\"4,0!\", ")
  .replace("\nv6 [", "\nv6 [label=\"\", pos=\"0,-2.2!\", ")
  .replace(
    "v1 -> v3 [pdg=1000, name=p1, mom=p1];",
    "v1 -> v3 [particle=\"about-photon\", name=p1, mom=p1, pos=\"-6,0!\"];",
  )
  .replace(
    "v5 -> v2 [pdg=1000, name=p2, mom=p2];",
    "v5 -> v2 [particle=\"about-photon\", name=p2, mom=p2, pos=\"6,0!\"];",
  )
  .replace(
    "v3 -> v4 [pdg=1000, name=q1];",
    "v3 -> v4 [particle=\"about-fermion\", name=q1, pos=\"-2,2.2!\"];",
  )
  .replace(
    "v4 -> v5 [pdg=1000, name=q2];",
    "v4 -> v5 [particle=\"about-fermion\", name=q2, pos=\"2,2.2!\"];",
  )
  .replace(
    "v5 -> v6 [pdg=1000, name=q3, lmb_index=0];",
    "v5 -> v6 [particle=\"about-fermion\", name=q3, lmb_index=0, pos=\"2,-2.2!\"];",
  )
  .replace(
    "v6 -> v3 [pdg=1000, name=q4];",
    "v6 -> v3 [particle=\"about-fermion\", name=q4, pos=\"-2,-2.2!\"];",
  )
  .replace(
    "v4 -> v6 [pdg=1000, name=q5, lmb_index=1];",
    "v4 -> v6 [particle=\"about-gluon\", name=q5, lmb_index=1, pos=\"0,0!\"];",
  )
  .replace(
    "\n}",
    "\ncut_blue [label=\"\", hidden=true, cut_curve=blue, cut_start=\"-0.9,4.3\", cut_end=\"-2.4,-4.3\", cut_control_start=\"-3.1,3.2\", cut_control_end=\"-0.2,-3.2\", pos=\"0,0!\"];\ncut_red [label=\"\", hidden=true, cut_curve=red, cut_start=\"3.1,4.2\", cut_end=\"-1.2,-4.2\", cut_control_start=\"-0.2,3.0\", cut_control_end=\"-0.6,-3.1\", pos=\"0,0!\"];\n}\n",
  )
)

#render(input, momentum-arrows: false, cut-curves: true)
