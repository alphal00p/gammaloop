"""Exercise default-off diagnostics and cooperative budgets in one native host."""

from symbolica import E, S
from symbolica.community import hepkit as fk

k, l = S("resource_policy::k", "resource_policy::l")
kin = fk.Kinematics()
source = fk.IntegralFamily([k], [], [kin.scalar_product(k, k) - 1], kinematics=kin)
target = fk.IntegralFamily([l], [], [kin.scalar_product(l, l) - 1], kinematics=kin)
operations = (
    lambda **kwargs: source.partial_fraction([1], **kwargs),
    lambda **kwargs: source.find_mapping(target, **kwargs),
    lambda **kwargs: fk.IntegralFamily.find_mappings([source, target], **kwargs),
)
for operation in operations:
    for invalid in (float("nan"), float("inf"), -1.0):
        try:
            operation(max_elapsed_seconds=invalid)
        except ValueError as error:
            assert "finite and nonnegative" in str(error)
        else:
            raise AssertionError("invalid elapsed-time setting was accepted")
    try:
        operation(max_elapsed_seconds=0.0)
    except fk.IntegralFamilyError as error:
        assert "elapsed" in str(error).lower()
    else:
        raise AssertionError("explicit zero-time budget was ignored")
    assert operation(max_process_rss_bytes=None, max_elapsed_seconds=None) is not None

assert source.find_mapping(target, max_candidates=None) is not None
assert source.partial_fraction([1], max_states=None) == [(E("1"), [1])]
try:
    source.find_mapping(target, max_candidates=0)
except fk.IntegralFamilyError:
    pass
else:
    raise AssertionError("explicit mapping diagnostic was ignored")

reducer = fk.TensorReducer.feynkit(S("resource_policy::D"))
try:
    reducer.with_resources(max_elapsed_seconds=0.0).reduce(E("1"))
except fk.TensorReductionError as error:
    assert "elapsed" in str(error).lower()
else:
    raise AssertionError("explicit tensor time budget was ignored")
assert reducer.with_resources().reduce(E("1")).to_expression() == E("1")
print("Default-off diagnostics, typed resource errors and healthy reuse passed")
