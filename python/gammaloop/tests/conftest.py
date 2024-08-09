import json
from subprocess import Popen, PIPE
import pytest
import os
from pathlib import Path
from gammaloop.tests.common import get_gamma_loop_interpreter, RESOURCES_PATH, pjoin
from gammaloop.interface.gammaloop_interface import CommandList
from gammaloop.misc.common import GL_PATH, GammaLoopError, logger

# Was intended to run with pytest --mypy but stupidly it won't read any mypy config file so it's unworkable.
# We will use pyright instead.
# def pytest_configure(config: pytest.Config):
#    plugin = config.pluginmanager.getplugin('mypy')
#    plugin.mypy_argv.append(  # type: ignore
#        "--check-untyped-defs")
#    plugin.mypy_argv.extend([  # type: ignore
#        '--config-file', pjoin(GL_PATH, os.path.pardir, os.path.pardir, 'pyproject.toml')])
#    plugin.mypy_argv.extend(  # type: ignore
#        ['--exclude', "'python/gammaloop/data/templates/drawing/combine_pages.py'"])
#    plugin.mypy_argv.extend(  # type: ignore
#        ['--exclude', "'python/gammaloop/data/models/*'"])

# pytest_plugins = ['pytest_profiling']

GIT_REVISION: str = 'N/A'


def pytest_addoption(parser):
    parser.addoption(
        "--runrust", action="store_true", default=False, help="run rust tests"
    )
    parser.addoption(
        "--codecheck", action="store_true", default=False, help="run code checks"
    )


def pytest_collection_modifyitems(config, items):
    run_rust = config.getoption("--runrust")
    run_codecheck = config.getoption("--codecheck")

    skip_rust = pytest.mark.skip(reason="need --runrust option to run")
    skip_codecheck = pytest.mark.skip(reason="need --codecheck option to run")

    for item in items:
        if "rust" in item.keywords and not run_rust:
            item.add_marker(skip_rust)
        if "codecheck" in item.keywords and not run_codecheck:
            item.add_marker(skip_codecheck)


def get_test_directory(tmpdir_factory: pytest.TempPathFactory, test_folder: str, output_for_rust_test: bool = False) -> Path:
    user_specified_test_dir_path = os.environ.get(
        'PYTEST_OUTPUT_PATH_FOR_RUST', None)
    if (not output_for_rust_test) or user_specified_test_dir_path == "TMP":
        test_output_path = Path(tmpdir_factory.mktemp(test_folder))
    elif user_specified_test_dir_path is None:
        test_output_path = Path(os.path.normpath(os.path.join(
            GL_PATH, os.path.pardir, os.path.pardir, 'src', 'test_resources', test_folder)))
    else:
        if not os.path.isdir(user_specified_test_dir_path):
            raise GammaLoopError(
                f"User specified path for pytest output for rust tests '{user_specified_test_dir_path}' does not exist.")
        test_output_path = Path(os.path.join(
            user_specified_test_dir_path, test_folder))

    if not os.path.isdir(test_output_path):
        os.makedirs(test_output_path)

    return test_output_path


@pytest.fixture(scope="session")
def sm_model_yaml_file(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    yaml_model_path = get_test_directory(
        tmpdir_factory, 'gammaloop_models', True).joinpath("sm.yaml")
    gloop.run(CommandList.from_string(
        f"import_model sm; export_model {yaml_model_path} --format yaml"))
    return yaml_model_path


@pytest.fixture(scope="session")
def scalars_model_yaml_file(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    yaml_model_path = get_test_directory(
        tmpdir_factory, 'gammaloop_models', True).joinpath("scalars.yaml")
    gloop.run(CommandList.from_string(
        f"import_model scalars; export_model {yaml_model_path} --format yaml"))
    return yaml_model_path


@pytest.fixture(scope="session")
def massless_scalar_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(
        tmpdir_factory, "TEST_AMPLITUDE_massless_scalar_triangle", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'massless_triangle.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_massless_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_massless_box", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'massless_box.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_fishnet_2x2_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_fishnet_2x2", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'fishnet_2x2.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_fishnet_2x3_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_fishnet_2x3", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'fishnet_2x3.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_cube_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_cube", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'cube.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_bubble_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_bubble", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'bubble.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_sunrise_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_sunrise", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'sunrise.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_double_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_double_triangle", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'double_triangle.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_mercedes_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_mercedes", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'mercedes.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_triangle_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_triangle_box", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'triangle_box.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_isopod_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_isopod", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'isopod.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_tree_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_tree_triangle", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'tree_triangle.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_ltd_topology_f_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_ltd_topology_f", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'ltd_topology_f.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_ltd_topology_h_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_ltd_topology_h", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'ltd_topology_h.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_raised_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_raised_triangle", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'raised_triangle.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def lbl_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_lbl_box", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'lbl_box.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def epem_a_ddx_nlo_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_CROSS_SECTION_epem_a_ddx_nlo", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'epem_a_ddx_NLO.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def massive_epem_a_ddx_nlo_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_CROSS_SECTION_massive_epem_a_ddx_nlo", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'epem_a_ddx_NLO.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_hexagon_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_hexagon", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'hexagon.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_ltd_topology_c_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_ltd_topology_c", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'ltd_topology_c.py')} -f qgraph --no_compile
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_massless_pentabox_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_massless_pentabox", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
generate_graph --name=massless_pentabox -ve=[(1000,1,6),(1000,6,7),(1000,7,2),(1000,2,1),(1000,7,3),(1000,3,4),(1000,4,5),(1000,5,6)] -ee=[("in",1),("in",2),("in",3),("in",4),("out",5)]
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_massless_3l_pentabox_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_massless_3l_pentabox", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
generate_graph --name=massless_3l_pentabox -ve=[(1000,1,6),(1000,6,7),(1000,7,2),(1000,2,1),(1000,7,8),(1000,8,9),(1000,9,6),(1000,8,3),(1000,3,4),(1000,4,5),(1000,5,9)] -ee=[("in",1),("in",2),("in",3),("in",4),("out",5)]
output {output_path} --overwrite_output"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_3L_6P_topology_A_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_3L_6P_topology_A", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'scalar_3L_6P_topology_A.py')} -f qgraph --no_compile
set target_omega 2.0
output {output_path} -exp -ef file"""))
    return output_path


@pytest.fixture(scope="session")
def physical_3L_6photons_topology_A_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_physical_3L_6photons_topology_A", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'physical_3L_6photons_topology_A.py')} -f qgraph --no_compile
output {output_path} --overwrite_output -exp -ef file"""))
    return output_path


@pytest.fixture(scope="session")
def compile_rust_tests() -> Path | None:

    # If you want to bypass the "manual" compilation of the rust tests, then uncomment the line below
    # return None

    cmd_list = ['cargo', 'build', '--release', '--target-dir', os.path.normpath(os.path.join(
        GL_PATH, os.path.pardir, os.path.pardir, 'rust_test_binaries')), '--features=binary,fail-on-warnings', '--no-default-features',
        '--tests', '--message-format=json']
    logger.debug("Compiling rust tests with command: %s", " ".join(cmd_list))
    process = Popen(cmd_list, cwd=GL_PATH, stdout=PIPE, stderr=PIPE)
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
        if json_obj["reason"] == "compiler-artifact" and 'gammalooprs' in json_obj["package_id"] and "lib" in json_obj["target"]["kind"] and json_obj["executable"] is not None:
            logger.critical(
                "Rust tests successfully compiled to binary '%s'", json_obj["executable"])
            return json_obj["executable"]

    return None
