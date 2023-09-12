from tests.common import *
from gammaloop.interface.gammaloop_interface import CommandList
import pytest

class TestLoadModel:

    def test_load_sm_from_ufo(self):
        gl = get_gamma_loop_interpreter()
        gl.run(CommandList.from_string("import_model sm --format ufo"))

        assert(len(gl.model.particles)==43)
        assert(len(gl.model.lorentz_structures)==22)
        assert(len(gl.model.couplings)==108)
        assert(len(gl.model.parameters)==72)
        assert(len(gl.model.vertex_rules)==153)
        assert(len(gl.model.orders)==2)

    # This test uses a session-wide fixture defined in conftest.py
    def test_load_sm_from_yaml(self, sm_model_yaml_file):
        gl = get_gamma_loop_interpreter()
        gl.run(CommandList.from_string(f"import_model {sm_model_yaml_file} --format yaml"))

        assert(len(gl.model.particles)==43)
        assert(len(gl.model.lorentz_structures)==22)
        assert(len(gl.model.couplings)==108)
        assert(len(gl.model.parameters)==72)
        assert(len(gl.model.vertex_rules)==153)
        assert(len(gl.model.orders)==2)

    def test_load_scalars(self):
        gl = get_gamma_loop_interpreter()
        gl.run(CommandList.from_string("import_model scalars"))
        
        assert(len(gl.model.particles)==3)
        assert(len(gl.model.lorentz_structures)==1)
        assert(len(gl.model.couplings)==1)
        assert(len(gl.model.parameters)==12)
        assert(len(gl.model.vertex_rules)==10)
        assert(len(gl.model.orders)==2)

class TestLoadQGraph:

    def test_epem_a_ddx_nlo(self):
        gl = get_gamma_loop_interpreter()
        gl.run(CommandList.from_string("import_model sm; import_qgraph %s --no_compile"%pjoin(RESOURCES_PATH,'qgraf_outputs','epem_a_ddx_NLO.py')))
        assert(True)