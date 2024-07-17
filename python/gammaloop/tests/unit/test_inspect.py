# pylint: disable=unused-import
from gammaloop.tests.common import get_gamma_loop_interpreter, check_inspect_result
import gammaloop.interface.gammaloop_interface as gl_interface
from pathlib import Path


class TestScalarTopologies:

    def test_inspect_scalar_3L_6P_topology_A(self, scalar_3L_6P_topology_A_export: Path):

        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList()
        command_list.add_command(
            "set externals.momenta [\
[5.,0.,0.,5.],\
[5.,0.,0.,-5.],\
[8.855133305450298e-1,-2.210069028768998e-1,4.008035319168533e-1,-7.580543095693663e-1],\
[3.283294192270986e0,-1.038496118834563e0,-3.019337553895401e0,7.649492138716588e-1],\
[1.523581094674306e0,-1.058809596665922e0,-9.770963832697570e-1,4.954838522679282e-1],\
[4.307611382509676e0,2.318312618377385e0,3.595630405248305e0,-5.023787565702210e-1],\
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
            inspect_res, complex(0, 2.2265511758628186e-44), max_relative_diff=1.0e-12)


class TestPhysicalTopologies:

    def test_inspect_physical_3L_6photons_topology_A(self, physical_3L_6photons_topology_A_export: Path):
        ######################################################################
        # TODO This test will need to be updated once numerators are supported!
        ######################################################################
        gl = get_gamma_loop_interpreter()

        command_list = gl_interface.CommandList()
        command_list.add_command(
            "launch {}".format(physical_3L_6photons_topology_A_export))
        command_list.add_command(
            "set externals.momenta [\
[5.,0.,0.,5.],\
[5.,0.,0.,-5.],\
[8.855133305450298e-1,-2.210069028768998e-1,4.008035319168533e-1,-7.580543095693663e-1],\
[3.283294192270986e0,-1.038496118834563e0,-3.019337553895401e0,7.649492138716588e-1],\
[1.523581094674306e0,-1.058809596665922e0,-9.770963832697570e-1,4.954838522679282e-1],\
[4.307611382509676e0,2.318312618377385e0,3.595630405248305e0,-5.023787565702210e-1],\
]")
        ######################################################################
        # TODO Helicity data choice will need to be supplied here as well.
        ######################################################################
        command_list.add_command("set e_cm 1.")
        command_list.add_command("set sampling {'type':'default'}")

        gl.run(command_list)

        inspect_res = gl.do_inspect(
            'physical_3L_6photons_topology_A -p 0.123 0.3242 0.4233 0.14235 0.25122 0.3245 0.12337 0.224237 0.32327')
        check_inspect_result(
            inspect_res, complex(0, 2.2265511758628186e-44), max_relative_diff=1.0e-12)
