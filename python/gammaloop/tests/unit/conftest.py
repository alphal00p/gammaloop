import pytest
from tests.common import *
from gammaloop.interface.gammaloop_interface import CommandList

@pytest.fixture(scope="session")
def sm_model_yaml_file(tmpdir_factory):
    gl = get_gamma_loop_interpreter()
    yaml_model_path = tmpdir_factory.mktemp("sm_model").join("sm_model.yaml")
    gl.run(CommandList.from_string(f"import_model sm; export_model {yaml_model_path} --format yaml"))
    return yaml_model_path