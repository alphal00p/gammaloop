#import "../lib.typ": *

#let M = mink(4)
#let B = bis(4)

// Explicit gamma ports are storage ordered: bispinor, bispinor, Minkowski.
#let rejected = gamma(slot(M, 1), slot(B, 1), slot(B, 2))
