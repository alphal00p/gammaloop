# pylint: disable=unused-import
from gammaloop.tests.common import scalar_euclidean_integration_test

max_mc_error_dif = 5.0
max_rel_error_dif = 0.01
max_percent_error = 0.01


class TestScalarTopologies:

    def test_box(self):
        assert True
        # TODO
        # gl = get_gamma_loop_interpreter()
        # gl.run(CommandList.from_string("import_model sm"))

    def test_scalar_triangle(self):

        scalar_euclidean_integration_test(
            "massless_triangle", True, 0.00009765455799148221)

    def test_tree_triangle(self):

        scalar_euclidean_integration_test(
            "tree_triangle", True, 0.00009765455799148221 / -49)

    def test_ltd_topology_f(self):

        scalar_euclidean_integration_test(
            "ltd_topology_f", True, -0.00000526647)

    def test_ltd_topology_h(self):
        scalar_euclidean_integration_test(
            "ltd_topology_h", False, -8.36515e-8)
