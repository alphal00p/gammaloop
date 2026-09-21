"""Covariant tensor identities through the installed shared reducer.

These exercise the external-momentum reduction used by FeynCalc's loop
examples; scalar integration and IBP reduction are separate operations.
"""

from symbolica import E, S
from symbolica.community import feynkit as fk

D, k, p, r = S(
    "tensor_check::D", "tensor_check::k", "tensor_check::p", "tensor_check::r"
)
mink, dot = S("spenso::mink", "spenso::dot")
mu, nu, rho, sigma = S(
    "tensor_check::mu", "tensor_check::nu", "tensor_check::rho", "tensor_check::sigma"
)
kc, pc, rc = k(mink(D)), p(mink(D)), r(mink(D))
reducer = fk.TensorReducer(D).with_integrated_vector(kc).with_external_vector(pc)

parallel = dot(kc, pc) * dot(rc, pc) / dot(pc, pc)
transverse = (dot(kc, kc) - dot(kc, pc) ** 2 / dot(pc, pc)) * (
    dot(rc, rc) - dot(rc, pc) ** 2 / dot(pc, pc)
)
numerator = E("1")
expected = [
    parallel,
    parallel**2 + transverse / (D - 1),
    parallel**3 + 3 * parallel * transverse / (D - 1),
    parallel**4
    + 6 * parallel**2 * transverse / (D - 1)
    + 3 * transverse**2 / ((D - 1) * (D + 1)),
]
for index, result in zip([mu, nu, rho, sigma], expected):
    numerator *= k(mink(D, index)) * r(mink(D, index))
    assert (reducer.reduce(numerator) - result).together() == E("0")

# Scalar loop invariants multiply the projector; they must not themselves be
# replaced by transverse norms during projection.
weighted = dot(kc, kc) * k(mink(D, mu))
assert (
    reducer.reduce(weighted) - dot(kc, kc) * dot(kc, pc) / dot(pc, pc) * p(mink(D, mu))
).together() == E("0")

# An auxiliary null direction makes the two-vector Gram matrix invertible.
# Imposing p²=r²=0 is safe after the covariant projection with p.r nonzero.
null_kin = (
    fk.Kinematics(D).with_scalar_product(p, p, E("0")).with_scalar_product(r, r, E("0"))
)
null_result = reducer.with_external_vector(rc).reduce(k(mink(D, mu)))
null_expected = (dot(kc, rc) * p(mink(D, mu)) + dot(kc, pc) * r(mink(D, mu))) / dot(
    pc, rc
)
assert (null_kin.apply(null_result) - null_expected).together() == E("0")
print("External-basis tensor identities through rank four passed")
