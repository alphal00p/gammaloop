"""Verify discovered loop shifts and scalar numerator/power mappings."""

from symbolica import E, S
from symbolica.community import feynkit as fk

k, q, l, r, p, s = S(
    "map_check::k", "map_check::q", "map_check::l", "map_check::r", "map_check::p", "s"
)
kin = fk.Kinematics(momenta=[k, q, l, r, p]).with_scalar_product(p, p, s)


def square(v):
    return kin.scalar_product(v, v)


source = fk.IntegralFamily(
    [k, q], [p], [square(k) - 1, square(q) - 2, square(k - q + p) - 3], kinematics=kin
)
target = fk.IntegralFamily(
    [l, r],
    [p],
    [square(l + p) - 3, square(l + r + p) - 1, square(r + p) - 2, square(l)],
    kinematics=kin,
)
mapping = source.find_mapping(target)
assert mapping is not None
assert mapping.denominator_map == [1, 2, 0]
assert mapping.map_powers([1, 2, -1]) == [-1, 1, 2, 0]
for i, j in enumerate(mapping.denominator_map):
    assert (
        mapping.apply(source.denominators[i]) - target.denominators[j]
    ).together() == E("0")

explicit = source.mapping_to(target, [l + r + p, r + p])
assert explicit is not None
assert (
    explicit.apply(kin.scalar_product(k, q)) - kin.scalar_product(l + r + p, r + p)
).expand() == E("0")
assert source.mapping_to(target, [2 * l, r]) is None

try:
    source.find_mapping(target, max_candidates=0)
except fk.IntegralFamilyError:
    pass
else:
    raise AssertionError("mapping search silently exceeded its budget")
print("Verified loop shifts, subtopology embeddings and numerator maps passed")
