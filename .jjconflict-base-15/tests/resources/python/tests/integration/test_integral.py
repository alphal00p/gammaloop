# pylint: disable=unused-import
from gammaloop.tests.common import get_gamma_loop_interpreter, check_integration_result
import gammaloop.interface.gammaloop_interface as gl_interface
from pathlib import Path
import pytest


class TestScalarTopologies:

    def test_box(self, scalar_massless_box_export: Path):
        target_re = -7.43707e-8
        target_im = 6.57830e-8

        gl = get_gamma_loop_interpreter()
        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(scalar_massless_box_export))
        command_list.add_command(
            "set externals.data.momenta [[14.0,-6.6,-40.,0.],[43.,-15.2,-33.,0.],[17.9,50.0,-11.8,0.0],]")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical','upcast_on_failure':False,'matrix_stability_test':1.0e-5}")
        command_list.add_command("set continuous_dim_learning_rate 0.0")
        command_list.add_command("set rotation_axis [{'type':'x'}]")
        command_list.add_command("set n_start 1_000_000")
        command_list.add_command("set n_max 1_000_000")
        command_list.add_command("integrate massless_box -r")
        gl.run(command_list)

        check_integration_result(target_re, scalar_massless_box_export)
        check_integration_result(
            target_im, scalar_massless_box_export, imag_phase=True)

    def test_scalar_triangle(self, massless_scalar_triangle_export: Path):
        target = -0.00009765455799148221
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(massless_scalar_triangle_export))
        command_list.add_command(
            "set externals.data.momenta [[1,3,4,5],[-1,-6,-7,-8],]")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical','upcast_on_failure':False,'matrix_stability_test':1.0e-5}")
        command_list.add_command("set n_start 10000")
        command_list.add_command("set n_max 10000")
        command_list.add_command("set seed 1")
        command_list.add_command("integrate massless_triangle -r")

        gl.run(command_list)
        check_integration_result(
            target, massless_scalar_triangle_export, imag_phase=False)

    def test_tree_triangle(self, scalar_tree_triangle_export: Path):
        target = 0.00009765455799148221 / 49

        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(scalar_tree_triangle_export))
        command_list.add_command(
            "set externals.data.momenta [[0.5,1.5,2,2.5],[0.5,1.5,2,2.5],[-1,-6,-7,-8],]")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical','upcast_on_failure':False,'matrix_stability_test':1.0e-5}")
        command_list.add_command("set n_start 10000")
        command_list.add_command("set n_max 10000")
        command_list.add_command("set seed 2")
        command_list.add_command("integrate tree_triangle -r")

        gl.run(command_list)
        check_integration_result(
            target, scalar_tree_triangle_export, imag_phase=False)

    def test_ltd_topology_f(self, scalar_ltd_topology_f_export: Path):
        target = -0.00000526647
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(scalar_ltd_topology_f_export))
        command_list.add_command(
            "set externals.data.momenta [[0,0,0,1],]")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical','upcast_on_failure':False,'matrix_stability_test':1.0e-5}")
        command_list.add_command("set n_start 100000")
        command_list.add_command("set n_max 100000")
        command_list.add_command("set seed 1")
        command_list.add_command("integrate ltd_topology_f -r")

        gl.run(command_list)

        check_integration_result(
            target, scalar_ltd_topology_f_export, imag_phase=False)

    def test_ltd_topology_h(self, scalar_ltd_topology_h_export: Path):
        target = 8.36515e-8
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(scalar_ltd_topology_h_export))
        command_list.add_command(
            "set externals.data.momenta [[0,0,0,1],]")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical','upcast_on_failure':False,'matrix_stability_test':1.0e-5}")
        command_list.add_command("set n_start 100000")
        command_list.add_command("set n_max 100000")
        command_list.add_command("set seed 1")
        command_list.add_command("integrate ltd_topology_h -r")

        gl.run(command_list)
        check_integration_result(
            target, scalar_ltd_topology_h_export, imag_phase=False)

    def test_scalar_3L_6P_topology_A_with_tropical_sampling(self, scalar_3L_6P_topology_A_export: Path):

        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList()

        command_list.add_command(
            "launch {}".format(scalar_3L_6P_topology_A_export))

        command_list.add_command("set_model_param mass_scalar_1 172.0")
        command_list.add_command(
            "set externals.data.momenta [\
[5.,0.,0.,5.],\
[5.,0.,0.,-5.],\
[8.855133305450298e-1,-2.210069028768998e-1,4.008035319168533e-1,-7.580543095693663e-1],\
[3.283294192270986e0,-1.038496118834563e0,-3.019337553895401e0,7.649492138716588e-1],\
[1.523581094674306e0,-1.058809596665922e0,-9.770963832697570e-1,4.954838522679282e-1],\
'dependent',\
]")

        command_list.add_command("set integrated_phase 'real'")
        command_list.add_command("set e_cm 10.")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical','upcast_on_failure':False,'matrix_stability_test':1.0e-5}")
        command_list.add_command("set n_start 40000")
        command_list.add_command("set n_max 40000")
        command_list.add_command("set continuous_dim_learning_rate 0.0")
        command_list.add_command("set seed 1")
        command_list.add_command("integrate scalar_3L_6P_topology_A -r")

        gl.run(command_list)
        check_integration_result(
            0., scalar_3L_6P_topology_A_export, imag_phase=True)
        # PySecDec reference run showed 2.8555(17)e-36
        check_integration_result(
            2.8555e-36, scalar_3L_6P_topology_A_export, imag_phase=False, max_mc_error_diff=5.0, max_rel_error_diff=0.1, max_percent_error=0.1)


class TestPhysicalTopologies:

    @pytest.mark.slow
    def test_physical_1L_6photons_with_tropical_sampling(self, physical_1L_6photons_export: Path):
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(physical_1L_6photons_export))
        command_list.add_command(
            "set externals.data.momenta [\
[500.0,0.,-300.,400.],\
[500.0,0.,300.,-400.],\
[88.551333054502976,-22.100690287689979,40.080353191685333,-75.805430956936632],\
[328.32941922709853,-103.84961188345630,-301.93375538954012,76.494921387165888],\
[152.35810946743061,-105.88095966659220,-97.709638326975707,49.548385226792817],\
]")
        command_list.add_command(
            "set externals.data.helicities [-1,-1,-1,-1,-1,-1]")
        command_list.add_command("set_model_param mz 91.188 -nu")
        command_list.add_command(
            "set_model_param gf 1.19874983504616246e-5 -nu")
        command_list.add_command("set_model_param mt 173.0 -nu")
        command_list.add_command("set_model_param ymt 173.0 -nu")
        command_list.add_command("set_model_param aewm1 128.93 -nu")
        command_list.add_command("set_model_param update_only 0.")

        command_list.add_command("set integrated_phase 'real'")
        command_list.add_command("set e_cm 1.")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical','upcast_on_failure':False}")
        command_list.add_command("set n_start 100000")
        command_list.add_command("set n_max 1000000")
        command_list.add_command("set seed 1")
        command_list.add_command(
            "integrate physical_1L_6photons -r -t (9.27759500687454717e-11,3.68394576249870544e-11)")

        gl.run(command_list)

        check_integration_result(
            9.27759500687454717e-11, physical_1L_6photons_export, imag_phase=False, max_mc_error_diff=5.0, max_rel_error_diff=5.0, max_percent_error=5.0)
        check_integration_result(
            3.68394576249870544e-11, physical_1L_6photons_export, imag_phase=True, max_mc_error_diff=5.0, max_rel_error_diff=10.0, max_percent_error=10.0)

    @pytest.mark.slow
    def test_euclidean_1L_6photons_with_tropical_sampling(self, physical_1L_6photons_export: Path):
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList.from_string(
            "launch {}".format(physical_1L_6photons_export))
        command_list.add_command(
            "set externals.data.momenta [\
[500.0,0.,-300.,400.],\
[500.0,0.,300.,-400.],\
[88.551333054502976,-22.100690287689979,40.080353191685333,-75.805430956936632],\
[328.32941922709853,-103.84961188345630,-301.93375538954012,76.494921387165888],\
[152.35810946743061,-105.88095966659220,-97.709638326975707,49.548385226792817],\
]")
        command_list.add_command(
            "set externals.data.helicities [-1,-1,-1,-1,-1,-1]")
        command_list.add_command("set_model_param mz 91.188 -nu")
        command_list.add_command(
            "set_model_param gf 1.19874983504616246e-5 -nu")
        command_list.add_command("set_model_param mt 1500.0 -nu")
        command_list.add_command("set_model_param ymt 1500.0 -nu")
        command_list.add_command("set_model_param aewm1 128.93 -nu")
        command_list.add_command("set_model_param update_only 0.")

        command_list.add_command("set integrated_phase 'real'")
        command_list.add_command("set e_cm 1.")
        command_list.add_command(
            "set sampling {'type':'discrete_graph_sampling','subtype':'tropical','upcast_on_failure':False}")
        command_list.add_command("set n_start 100000")
        command_list.add_command("set n_max 1000000")
        command_list.add_command("set seed 1")
        command_list.add_command(
            "integrate physical_1L_6photons -r -t (1.22898408452706e-13,3.94362534040412e-13)")

        gl.run(command_list)

        check_integration_result(
            1.22898408452706e-13, physical_1L_6photons_export, imag_phase=False, max_mc_error_diff=5.0, max_rel_error_diff=1.0, max_percent_error=1.0)
        check_integration_result(
            3.94362534040412e-13, physical_1L_6photons_export, imag_phase=True, max_mc_error_diff=5.0, max_rel_error_diff=1.0, max_percent_error=1.0)
