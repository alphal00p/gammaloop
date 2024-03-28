import json
from subprocess import Popen, PIPE
import pytest
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


@pytest.fixture(scope="session")
def sm_model_yaml_file(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    yaml_model_path = Path(tmpdir_factory.mktemp(
        "sm_model")).joinpath("sm_model.yaml")
    gloop.run(CommandList.from_string(
        f"import_model sm; export_model {yaml_model_path} --format yaml"))
    return yaml_model_path


@pytest.fixture(scope="session")
def scalar_massless_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_massless_triangle_scalar")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','massless_triangle.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_massless_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_massless_box_scalar")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','massless_box.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_fishnet_2x2_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_fishnet_2x2")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','fishnet_2x2.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_fishnet_2x3_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_fishnet_2x3")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','fishnet_2x3.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_cube_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_cube")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','cube.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_bubble_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_bubble")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','bubble.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_sunrise_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_bubble")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','sunrise.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_double_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_double_triangle")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','double_triangle.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_mercedes_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_mercedes")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','mercedes.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_triangle_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_triangle_box")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','triangle_box.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_isopod_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_isopod")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','isopod.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_tree_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_tree_triangle")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','tree_triangle.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_ltd_topology_f_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_ltd_topology_f")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','ltd_topology_f.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def scalar_ltd_topology_h_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_ltd_topology_h")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','ltd_topology_h.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def lbl_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_AMPLITUDE_lbl_box")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm;
import_graphs {pjoin(RESOURCES_PATH,'qgraf_outputs','lbl_box.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def epem_a_ddx_nlo_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_CROSS_SECTION_epem_a_ddx_nlo")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'epem_a_ddx_NLO.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def massive_epem_a_ddx_nlo_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = Path(tmpdir_factory.mktemp(
        "TEST_CROSS_SECTION_massive_epem_a_ddx_nlo")).joinpath("OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full
import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'epem_a_ddx_NLO.py')} -f qgraph --no_compile
output {output_path}"""))
    return output_path


@pytest.fixture(scope="session")
def compile_rust_tests() -> Path | None:
    process = Popen(['cargo', 'build', '--release', '--features=binary,fail-on-warnings', '--no-default-features',
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
        if json_obj["reason"] == "compiler-artifact" and 'gammalooprs' in json_obj["package_id"] and "lib" in json_obj["target"]["kind"] and json_obj["executable"] is not None:
            logger.critical(
                "Rust tests successfully compiled to binary '%s'", json_obj["executable"])
            return json_obj["executable"]

    return None
