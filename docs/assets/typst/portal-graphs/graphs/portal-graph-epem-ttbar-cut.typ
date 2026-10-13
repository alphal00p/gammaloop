#import "../figure.typ": render

// Stable drawing-only fixture for the portal's open external-leg illustration.
// Preserve statement order; canonical physics fixture regeneration is independent.
#let input = ```dot
digraph feynkit {
  v0 [id=0];
  v1 [id=1];
  v2 [id=2];
  v3 [id=3];
  v4 [id=4];
  v5 [id=5];
  v6 [id=6];
  v7 [id=7];
  v8 [id=8];
  v9 [id=9];
  v10 [id=10];
  v11 [id=11];
  v12 [id=12];
  v13 [id=13];
  v13 -> v0 [id=0, particle="e-"];
  v1 -> v13 [id=1, particle="e-"];
  v12 -> v2 [id=2, particle="e+"];
  v3 -> v12 [id=3, particle="e+"];
  v4 -> v5 [id=4, particle="t"];
  v4 -> v5 [id=5, particle="H", dir=none];
  v11 -> v4 [id=6, particle="t"];
  v5 -> v10 [id=7, particle="t"];
  v6 -> v7 [id=8, particle="t"];
  v9 -> v6 [id=9, particle="t"];
  v6 -> v13 [id=10, particle="a", dir=none];
  v7 -> v8 [id=11, particle="t"];
  v7 -> v12 [id=12, particle="a", dir=none];
  v8 -> v9 [id=13, particle="g", dir=none];
  v8 -> v11 [id=14, particle="t"];
  v10 -> v9 [id=15, particle="t"];
  v10 -> v11 [id=16, particle="g", dir=none];
}
```.text
#render(
  input,
  cross-section-mode: true,
  momentum-arrows: true,
)
