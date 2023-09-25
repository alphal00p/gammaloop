import pytest
from gammaloop.tests.common import get_gamma_loop_interpreter, RESOURCES_PATH, pjoin
from gammaloop.interface.gammaloop_interface import CommandList


@pytest.fixture(scope="session")
def sm_model_yaml_file(tmpdir_factory):
    gloop = get_gamma_loop_interpreter()
    yaml_model_path = tmpdir_factory.mktemp("sm_model").join("sm_model.yaml")
    gloop.run(CommandList.from_string(
        f"import_model sm; export_model {yaml_model_path} --format yaml"))
    return yaml_model_path


@pytest.fixture(scope="session")
def scalar_massless_triangle_export(tmpdir_factory):
    gloop = get_gamma_loop_interpreter()
    output_path = tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_massless_triangle_scalar").join("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','massless_triangle.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def epem_a_ddx_nlo_export(tmpdir_factory):
    gloop = get_gamma_loop_interpreter()
    output_path = tmpdir_factory.mktemp(
        "TEST_CROSS_SECTION_epem_a_ddx_nlo").join("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'epem_a_ddx_NLO.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path
