import json
from subprocess import Popen, PIPE
import pytest
from gammaloop.tests.common import get_gamma_loop_interpreter, RESOURCES_PATH, pjoin
from gammaloop.interface.gammaloop_interface import CommandList
from gammaloop.misc.common import GL_PATH, GammaLoopError, logger


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
def scalar_fishnet_2x2_export(tmpdir_factory):
    gloop = get_gamma_loop_interpreter()
    output_path = tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_fishnet_2x2").join("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','fishnet_2x2.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_fishnet_2x3_export(tmpdir_factory):
    gloop = get_gamma_loop_interpreter()
    output_path = tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_fishnet_2x3").join("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','fishnet_2x3.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path

@pytest.fixture(scope="session")
def scalar_cube_export(tmpdir_factory):
    gloop = get_gamma_loop_interpreter()
    output_path = tmpdir_factory.mktemp("TEST_AMPLITUDE_cube").join("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','cube.py')} -f qgraph --no_compile
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


@pytest.fixture(scope="session")
def compile_rust_tests():
    process = Popen(['cargo', 'build', '--release', '--features=binary', '--no-default-features',
                     '--tests', '--message-format=json'], cwd=GL_PATH, stdout=PIPE, stderr=PIPE)
    logger.critical("Compiling rust tests...")
    output, err = process.communicate()
    if process.returncode != 0:
        raise GammaLoopError(
            "Failed to compile rust tests. Error:\n" + err.decode("utf-8"))

    compiler_artifact = output.decode("utf-8")
    for json_line in reversed(compiler_artifact.split("\n")):
        if json_line == "":
            continue
        try:
            json_obj = json.loads(json_line)
        except json.decoder.JSONDecodeError:
            continue
        if json_obj["reason"] == "compiler-artifact" and json_obj["package_id"].startswith('gammalooprs') and "lib" in json_obj["target"]["kind"] and json_obj["executable"] is not None:
            logger.critical(
                "Rust tests successfully compiled to binary '%s'", json_obj["executable"])
            return json_obj["executable"]
    raise GammaLoopError(
        "Failed to find executable in compiler artifact:\n"+compiler_artifact)
