#import "@local/tydenso:0.1.0": *

#set page(width: auto, height: auto, margin: 12pt)

#let V = mink(4)
#let p = vector("p")
#let mu = slot(V, 1)

$ #p(mu) = #to-typst(p(mu)) $
