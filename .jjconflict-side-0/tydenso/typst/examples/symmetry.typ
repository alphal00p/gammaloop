#import "../lib.typ": *

#set page(width: auto, height: auto, margin: 12pt)

#let V = mink(4)
#let mu = slot(V, 1)
#let nu = slot(V, 2)
#let F = tensor("F", antisymmetric: true)
#let cancellation = add(F(mu, nu), F(nu, mu))

$ F^(std.sym.mu std.sym.nu) + F^(std.sym.nu std.sym.mu)
    = #to-typst(cancellation) $
