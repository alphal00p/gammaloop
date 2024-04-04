# pylint: disable=unused-import
from gammaloop.tests.common import get_gamma_loop_interpreter, check_integration_result
import gammaloop.interface.gammaloop_interface as gl_interface
from pathlib import Path


max_mc_error_dif = 5.0
max_rel_error_dif = 0.01
max_percent_error = 0.01


class TestScalarTopologies:

    def test_box(self):
        assert True
        # TODO
        # gl = get_gamma_loop_interpreter()
        # gl.run(CommandList.from_string("import_model sm"))

    def test_scalar_triangle(self, scalar_massless_triangle_export: Path):
        target = 0.00009765455799148221
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(scalar_massless_triangle_export))
        command_list.add_command(
            "set externals.momenta [[1,3,4,5],[-1,-6,-7,-8],]")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical'}")
        command_list.add_command("set n_start 10000")
        command_list.add_command("set n_max 10000")
        command_list.add_command("set seed 1")
        command_list.add_command("integrate massless_triangle -r")

        gl.run(command_list)
        check_integration_result(
            True, target, scalar_massless_triangle_export)

    def test_tree_triangle(self, scalar_tree_triangle_export: Path):
        target = 0.00009765455799148221 / -49

        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(scalar_tree_triangle_export))
        command_list.add_command(
            "set externals.momenta [[0.5,1.5,2,2.5],[0.5,1.5,2,2.5],[-1,-6,-7,-8],]")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical'}")
        command_list.add_command("set n_start 10000")
        command_list.add_command("set n_max 10000")
        command_list.add_command("set seed 1")
        command_list.add_command("integrate tree_triangle -r")

        gl.run(command_list)
        check_integration_result(
            True, target, scalar_tree_triangle_export)

    def test_ltd_topology_f(self, scalar_ltd_topology_f_export: Path):
        target = -0.00000526647
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(scalar_ltd_topology_f_export))
        command_list.add_command(
            "set externals.momenta [[0,0,0,1],]")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical'}")
        command_list.add_command("set n_start 100000")
        command_list.add_command("set n_max 100000")
        command_list.add_command("set seed 1")
        command_list.add_command("integrate ltd_topology_f -r")

        gl.run(command_list)

        check_integration_result(
            True, target, scalar_ltd_topology_f_export)

    def test_ltd_topology_h(self, scalar_ltd_topology_h_export: Path):
        target = -8.36515e-8
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(scalar_ltd_topology_h_export))
        command_list.add_command(
            "set externals.momenta [[0,0,0,1],]")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical'}")
        command_list.add_command("set n_start 100000")
        command_list.add_command("set n_max 100000")
        command_list.add_command("set seed 1")
        command_list.add_command("integrate ltd_topology_h -r")

        gl.run(command_list)
        check_integration_result(
            False, target, scalar_ltd_topology_h_export)
