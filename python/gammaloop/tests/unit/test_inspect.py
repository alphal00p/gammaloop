# pylint: disable=unused-import
from gammaloop.tests.common import get_gamma_loop_interpreter, check_inspect_result
import gammaloop.interface.gammaloop_interface as gl_interface
from pathlib import Path


class TestScalarTopologies:

    def test_inspect_scalar_3L_6P_topology_A(self, scalar_3L_6P_topology_A_export: Path):

        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList()
        command_list.add_command(
            "set externals.data.momenta [\
[5.,0.,0.,5.],\
[5.,0.,0.,-5.],\
[8.855133305450298e-1,-2.210069028768998e-1,4.008035319168533e-1,-7.580543095693663e-1],\
[3.283294192270986e0,-1.038496118834563e0,-3.019337553895401e0,7.649492138716588e-1],\
[1.523581094674306e0,-1.058809596665922e0,-9.770963832697570e-1,4.954838522679282e-1],\
]")
        command_list.add_command("set integrated_phase 'imag'")
        command_list.add_command("set e_cm 1.")
        command_list.add_command("set sampling {'type':'default'}")

        command_list.add_command("launch {}".format(
            scalar_3L_6P_topology_A_export))
        command_list.add_command("set_model_param mass_scalar_1 172.0")

        gl.run(command_list)

        inspect_res = gl.do_inspect(
            'scalar_3L_6P_topology_A -p 0.123 0.3242 0.4233 0.14235 0.25122 0.3245 0.12337 0.224237 0.32327')
        check_inspect_result(
            inspect_res['final_result'], complex(2.2265511758628186e-44, 0.), max_relative_diff=1.0e-12)


class TestPhysicalTopologies:

    def test_inspect_euclidean_1L_6photons(self, physical_1L_6photons_export: Path):
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList()
        command_list.add_command(
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
        command_list.add_command("set e_cm 1.")
        command_list.add_command("set sampling {'type':'default'}")

        gl.run(command_list)

        inspect_res = gl.do_inspect(
            'physical_1L_6photons -p 0.123 0.3242 0.4233 0.14235 0.25122 0.3245 0.12337 0.224237 0.32327')
        check_inspect_result(
            inspect_res['final_result'], complex(-6.511498992086646e-16, 5.698855018015016e-16), max_relative_diff=1.0e-12)

    def test_inspect_physical_1L_6photons(self, physical_1L_6photons_export: Path):
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList()
        command_list.add_command(
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
        command_list.add_command("set e_cm 1.")
        command_list.add_command("set sampling {'type':'default'}")
        command_list.add_command(
            "set force_global_center [[115.91508160079663,179.7810253192939,-25.11851705401448]]")

        gl.run(command_list)

        inspect_res = gl.do_inspect(
            'physical_1L_6photons -p 0.123 0.3242 0.4233')
        check_inspect_result(
            inspect_res['final_result'], complex(
                3.2329161713509245e-16, -1.2667260904767698e-16), max_relative_diff=1.0e-12)
