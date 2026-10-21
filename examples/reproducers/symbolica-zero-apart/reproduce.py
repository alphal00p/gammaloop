"""Expose Symbolica939's zero-polynomial multivariate-apart panic."""
from symbolica import E, S

x, y = S("zero_apart_x", "zero_apart_y")
assert E("0").apart(x, y) == E("0")
