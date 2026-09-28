#set page(width: 637pt, height: 189pt, margin: 0pt, fill: none)

// The canonical Spenso wordmark is four independently weighted text layers.
// Roboto and the placements reproduce the upstream geometry while keeping the
// wording, overlap, opacity and spacing directly editable.
#place(top + left)[
  #block(width: 637pt, height: 189pt, fill: rgb("#d9d9d9"), radius: 49pt)
]
#place(top + left, dx: 37.5pt, dy: 40pt)[
  #text(font: "Roboto", size: 128pt, fill: black.transparentize(50%))[sp]
]
#place(top + left, dx: 175.5pt, dy: 40pt)[
  #text(font: "Roboto", size: 128pt, fill: black.transparentize(25%))[enso]
]
#place(top + left, dx: 450.5pt, dy: 40pt)[
  #text(font: "Roboto", size: 128pt, fill: black.transparentize(70%))[.rs]
]
#place(top + left, dx: 106.5pt, dy: 40pt)[
  #text(font: "Roboto", size: 128pt, fill: black.transparentize(50%))[d]
]
