// One-loop Feynman-diagram drawing helpers (cetz).
// Line types: "F" fermion (solid, arrow), "V" gauge boson (wavy), "S" scalar (dashed, arrow),
//             "G" ghost (dotted, arrow), "K" BRST source (double solid, external only).
// Vertex marks: "dot" ordinary vertex, "delta" insertion of the breaking operator (crossed circle),
//               "ct" counterterm (cross), "src" BRST-source vertex (square).
#import "@preview/cetz:0.4.2": canvas, draw, decorations

#let _stroke(kind) = {
  if kind == "S" { (dash: "dashed", thickness: 0.6pt) }
  else if kind == "G" { (dash: "dotted", thickness: 0.9pt) }
  else { (thickness: 0.6pt) }
}

// straight external leg from point a (outside) to b (vertex); arrow: 0 none, 1 towards b, -1 away from b
#let _leg(a, b, kind, arrow: 0) = {
  import draw: *
  if kind == "V" {
    decorations.wave(line(a, b), amplitude: 0.07, segment-length: 0.16, stroke: 0.6pt)
  } else if kind == "K" {
    line(a, b, stroke: 1.6pt)
    line(a, b, stroke: 0.7pt + white)
  } else {
    line(a, b, stroke: _stroke(kind))
    if arrow != 0 {
      let (p, q) = if arrow == 1 { (a, b) } else { (b, a) }
      line(p, (p, 58%, q), stroke: 0.6pt, mark: (end: ">", fill: black, scale: 0.55))
    }
  }
}

// arc of the loop circle (centre origin, radius r) from angle t0 to angle t1 (degrees, t1 < t0: clockwise)
#let _arc(r, t0, t1, kind, arrow: 0) = {
  import draw: *
  let p0 = (r * calc.cos(t0 * 1deg), r * calc.sin(t0 * 1deg))
  if kind == "V" {
    decorations.wave(arc(p0, start: t0 * 1deg, stop: t1 * 1deg, radius: r), amplitude: 0.07, segment-length: 0.16, stroke: 0.6pt)
  } else {
    arc(p0, start: t0 * 1deg, stop: t1 * 1deg, radius: r, stroke: _stroke(kind))
    if arrow != 0 {
      let tm = (t0 + t1) / 2
      let dt = if arrow == 1 { -6 } else { 6 }   // arrow == 1: along the direction t0 -> t1
      let pa = (r * calc.cos((tm - dt) * 1deg), r * calc.sin((tm - dt) * 1deg))
      let pb = (r * calc.cos((tm + dt) * 1deg), r * calc.sin((tm + dt) * 1deg))
      line(pa, pb, stroke: 0.6pt, mark: (end: ">", fill: black, scale: 0.55))
    }
  }
}

#let _vmark(p, kind) = {
  import draw: *
  if kind == "delta" {
    circle(p, radius: 0.13, fill: white, stroke: 0.6pt)
    line((p.at(0) - 0.09, p.at(1) - 0.09), (p.at(0) + 0.09, p.at(1) + 0.09), stroke: 0.6pt)
    line((p.at(0) - 0.09, p.at(1) + 0.09), (p.at(0) + 0.09, p.at(1) - 0.09), stroke: 0.6pt)
  } else if kind == "ct" {
    line((p.at(0) - 0.1, p.at(1) - 0.1), (p.at(0) + 0.1, p.at(1) + 0.1), stroke: 1pt)
    line((p.at(0) - 0.1, p.at(1) + 0.1), (p.at(0) + 0.1, p.at(1) - 0.1), stroke: 1pt)
  } else if kind == "src" {
    rect((p.at(0) - 0.09, p.at(1) - 0.09), (p.at(0) + 0.09, p.at(1) + 0.09), fill: black)
  } else {
    circle(p, radius: 0.045, fill: black)
  }
}

// One-loop 1PI diagram.
//   verts: array of dictionaries, in clockwise order starting on the left:
//          (legs: ((kind, label, arrow), ...), mark: "dot" | "delta" | "ct" | "src")
//          arrow of a leg: 1 = flows into the vertex, -1 = flows out, 0 = none
//   props: array of (kind, label, arrow), props.at(i) joins vertex i to vertex i+1 (cyclically);
//          arrow 1 = flows from vertex i to vertex i+1, -1 = the opposite, 0 = none.
//          A single vertex with one prop gives a tadpole (closed line attached to that vertex).
#let oneloop(verts, props, r: 0.62, leg: 0.85, caption: none, scale: 1.0) = {
  let n = verts.len()
  let fig = canvas(length: 1cm * scale, {
    import draw: *
    if n == 1 {
      // tadpole: vertex at origin, loop above
      let v = verts.at(0)
      let p = props.at(0)
      _arc(r, 270, -90, p.at(0), arrow: p.at(2))
      let c = (0, -r)
      let k = v.legs.len()
      for (j, l) in v.legs.enumerate() {
        let ang = if k == 1 { 270 } else { 200 + 140 * j / (k - 1) }
        let a = (c.at(0) + leg * calc.cos(ang * 1deg), c.at(1) + leg * calc.sin(ang * 1deg))
        _leg(a, c, l.at(0), arrow: l.at(2))
        content((c.at(0) + (leg + 0.22) * calc.cos(ang * 1deg), c.at(1) + (leg + 0.22) * calc.sin(ang * 1deg)), text(size: 7.5pt, l.at(1)))
      }
      if p.at(1) != none { content((0, r + 0.2), text(size: 7.5pt, p.at(1))) }
      _vmark(c, v.at("mark", default: "dot"))
    } else {
      let ang(i) = 180 - 360 * i / n
      for (i, p) in props.enumerate() {
        _arc(r, ang(i), ang(i) - 360 / n, p.at(0), arrow: p.at(2))
        if p.at(1) != none {
          let tm = ang(i) - 180 / n
          content(((r + 0.24) * calc.cos(tm * 1deg), (r + 0.24) * calc.sin(tm * 1deg)), text(size: 7.5pt, p.at(1)))
        }
      }
      for (i, v) in verts.enumerate() {
        let t = ang(i)
        let c = (r * calc.cos(t * 1deg), r * calc.sin(t * 1deg))
        let k = v.legs.len()
        for (j, l) in v.legs.enumerate() {
          let off = if k == 1 { 0 } else { -28 + 56 * j / (k - 1) }
          let a = (c.at(0) + leg * calc.cos((t + off) * 1deg), c.at(1) + leg * calc.sin((t + off) * 1deg))
          _leg(a, c, l.at(0), arrow: l.at(2))
          content((c.at(0) + (leg + 0.24) * calc.cos((t + off) * 1deg), c.at(1) + (leg + 0.24) * calc.sin((t + off) * 1deg)), text(size: 7.5pt, l.at(1)))
        }
        _vmark(c, v.at("mark", default: "dot"))
      }
    }
  })
  if caption == none { fig } else {
    box(align(center, stack(dir: ttb, spacing: 3pt, fig, text(size: 7.5pt, caption))))
  }
}

// lay several diagrams out in a row
#let diagrow(..ds) = align(center, stack(dir: ltr, spacing: 9mm, ..ds.pos()))
