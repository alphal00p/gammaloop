import functools
import sys
import shutil
import time
from gammaloop.misc.common import GL_PATH, GammaLoopError, logger
from gammaloop.interface.gammaloop_interface import CommandList
from gammaloop.tests.common import get_gamma_loop_interpreter, get_gamma_loop_interpreter_no_compilation, RESOURCES_PATH, pjoin
from pathlib import Path
import json
from subprocess import Popen, PIPE
import pytest
import os

FILELOCK_AVAILABLE = True
try:
    from filelock import FileLock  # type: ignore
except:
    FILELOCK_AVAILABLE = False
    print("Warning: filelock not installed. Please install it with 'pip install filelock'. Runtimes saving will not be thread-safe.")

COLORAMA_AVAILABLE = True
FIRST_COLORAMA_OUTPUT = True
try:
    import colorama  # type: ignore
    # Initialize colorama
    colorama.init()
except:
    COLORAMA_AVAILABLE = False

fixture_setup_times = {}


def get_terminal_width():
    return shutil.get_terminal_size().columns


def write_current_test_name(nodeid):
    # Leave space for the progress marker
    terminal_width = get_terminal_width()-8
    test_name = f"{nodeid}"
    max_test_name_length = min(len(test_name), terminal_width - 1)

    # Truncate test name if necessary
    if len(test_name) > max_test_name_length:
        test_name = '...' + test_name[-(max_test_name_length - 3):]

    # Move cursor to the right position and write the test name
    output = f"\033[s\033[{terminal_width - max_test_name_length}G{test_name}\033[u"  # nopep8
    sys.stderr.write(output)
    sys.stderr.flush()


def clear_current_test_name():
    terminal_width = get_terminal_width()
    # Clear the area where the test name was displayed
    blank_space = ' ' * (terminal_width // 2)
    output = f"\033[s\033[{terminal_width - len(blank_space)}G{blank_space}\033[u"  # nopep8
    sys.stderr.write(output)
    sys.stderr.flush()


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_logstart(nodeid, location):
    global FIRST_COLORAMA_OUTPUT
    if COLORAMA_AVAILABLE:
        if FIRST_COLORAMA_OUTPUT:
            FIRST_COLORAMA_OUTPUT = False
        else:
            write_current_test_name(nodeid)
    yield


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_logfinish(nodeid, location):
    if COLORAMA_AVAILABLE:
        clear_current_test_name()
    yield

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
    parser.addoption(
        "--max-runtime",
        action="store",
        type=float,
        default=None,
        help="Run only tests that last less than the specified time (in seconds).",
    )
    parser.addoption(
        "--update-runtime",
        action="store_true",
        default=False,
        help="Update the stored test runtimes.",
    )


def measure_fixture_setup_time(*fixture_args, **fixture_kwargs):
    """
    Decorator factory to measure the setup time of a fixture.

    Accepts the same arguments as @measure_fixture_setup_time.
    """
    def decorator(fixture_func):
        @pytest.fixture(*fixture_args, **fixture_kwargs)
        @functools.wraps(fixture_func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            result = fixture_func(*args, **kwargs)
            setup_duration = time.time() - start_time
            fixture_name = fixture_func.__name__
            if fixture_name not in fixture_setup_times:
                fixture_setup_times[fixture_name] = setup_duration
            return result
        return wrapper
    return decorator


def get_all_fixtures(item):
    """Recursively collect all fixtures used by a test item."""
    all_fixtures = set()
    stack = list(item._fixtureinfo.name2fixturedefs.keys())
    while stack:
        fixture_name = stack.pop()
        if fixture_name not in all_fixtures:
            all_fixtures.add(fixture_name)
            fixturedefs = item._fixtureinfo.name2fixturedefs.get(
                fixture_name, [])
            for fixturedef in fixturedefs:
                if fixturedef.argnames:
                    stack.extend(fixturedef.argnames)
    return all_fixtures


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    # Execute all other hooks to obtain the report object
    outcome = yield
    report = outcome.get_result()
    if report.when == "call" and report.passed:
        duration = report.duration
        if duration is not None:

            # Collect all fixtures used by the test, including indirect dependencies
            used_fixtures = get_all_fixtures(item)
            for fix in used_fixtures:
                if fix in ['tmpdir_factory', 'request', 'snapshot_check']:
                    continue
                if fix not in fixture_setup_times:
                    print(f"WARNING: setup time for fixture '{fix}' is not recorded. Make sure you decorated it with '@measure_fixture_setup_time(scope=\"session\")'.")  # nopep8
            fixtures_duration = sum(fixture_setup_times.get(fix, 0)
                                    for fix in used_fixtures)
            test_runtime = duration + fixtures_duration

            cache_dir = item.config.cache._cachedir

            # Get the existing runtimes from the cache
            runtimes = item.config.cache.get("test_runtimes", {})
            update_runtime = item.config.getoption("--update-runtime")
            runtime_exists = item.nodeid in runtimes

            # Use a file lock to ensure thread-safe cache access
            if FILELOCK_AVAILABLE:
                lock_file = os.path.join(cache_dir, "test_runtimes.lock")
                with FileLock(lock_file):  # type: ignore
                    # Update the runtime if --update-runtime is set OR there is no existing runtime
                    if update_runtime or not runtime_exists:
                        # Update the runtime for this test
                        runtimes[item.nodeid] = test_runtime
                        # Save the updated runtimes back to the cache
                        item.config.cache.set("test_runtimes", runtimes)
            else:
                update_runtime = item.config.getoption("--update-runtime")
                if update_runtime or not runtime_exists:
                    # Get the existing runtimes from the cache
                    runtimes = item.config.cache.get("test_runtimes", {})
                    # Update the runtime for this test
                    runtimes[item.nodeid] = test_runtime
                    # Save the updated runtimes back to the cache
                    item.config.cache.set("test_runtimes", runtimes)


def pytest_collection_modifyitems(config, items):
    max_runtime = config.getoption("--max-runtime")
    if max_runtime is not None:
        runtimes = config.cache.get("test_runtimes", {})
        selected_items = []
        deselected_items = []

        for item in items:
            duration = runtimes.get(item.nodeid, None)
            if duration is not None and duration <= max_runtime:
                selected_items.append(item)
            else:
                deselected_items.append(item)

        if deselected_items:
            config.hook.pytest_deselected(items=deselected_items)
            items[:] = selected_items

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


@measure_fixture_setup_time(scope="session")
def sm_model_yaml_file(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    yaml_model_path = get_test_directory(
        tmpdir_factory, 'gammaloop_models', True).joinpath("sm.yaml")
    gloop.run(CommandList.from_string(
        f"import_model sm; export_model {yaml_model_path} --format yaml"))
    return yaml_model_path


@measure_fixture_setup_time(scope="session")
def scalars_model_yaml_file(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    yaml_model_path = get_test_directory(
        tmpdir_factory, 'gammaloop_models', True).joinpath("scalars.yaml")
    gloop.run(CommandList.from_string(
        f"import_model scalars; export_model {yaml_model_path} --format yaml"))
    return yaml_model_path


@measure_fixture_setup_time(scope="session")
def massless_scalar_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(
        tmpdir_factory, "TEST_AMPLITUDE_massless_scalar_triangle", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'massless_triangle.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope='session')
def scalar_massless_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_massless_box", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'massless_box.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_fishnet_2x2_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_fishnet_2x2", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'fishnet_2x2.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_fishnet_2x3_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter_no_compilation()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_fishnet_2x3", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'fishnet_2x3.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_cube_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_cube", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'cube.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_bubble_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_bubble", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'bubble.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_sunrise_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_sunrise", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'sunrise.dot')}
output {output_path} --overwrite_output --yaml_only"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_double_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_double_triangle", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'double_triangle.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_mercedes_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_mercedes", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'mercedes.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_triangle_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_triangle_box", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'triangle_box.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_isopod_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_isopod", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'isopod.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_tree_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_tree_triangle", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'tree_triangle.dot')}
output {output_path}"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_ltd_topology_f_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_ltd_topology_f", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'ltd_topology_f.dot')}
output {output_path}"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_ltd_topology_h_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_ltd_topology_h", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'ltd_topology_h.dot')}
output {output_path}"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_raised_triangle_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_raised_triangle", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'raised_triangle.dot')}
output {output_path} --overwrite_output --yaml_only"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def lbl_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_lbl_box", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'lbl_box.dot')}
output {output_path} --overwrite_output --yaml_only"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def epem_a_ddx_nlo_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_CROSS_SECTION_epem_a_ddx_nlo", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'epem_a_ddx_NLO.dot')}
output {output_path} --yaml_only"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def massive_epem_a_ddx_nlo_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_CROSS_SECTION_massive_epem_a_ddx_nlo", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'epem_a_ddx_NLO.dot')}
output {output_path} --yaml_only"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_hexagon_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_hexagon", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'hexagon.dot')}
output {output_path} --overwrite_output --yaml_only"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_ltd_topology_c_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_ltd_topology_c", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'ltd_topology_c.dot')}
output {output_path} --overwrite_output --yaml_only"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_massless_pentabox_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_massless_pentabox", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
define_graph --name=massless_pentabox -ve=[(1000,1,6),(1000,6,7),(1000,7,2),(1000,2,1),(1000,7,3),(1000,3,4),(1000,4,5),(1000,5,6)] -ee=[("in",1),("in",2),("in",3),("in",4),("out",5)]
output {output_path} --overwrite_output --yaml_only"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_massless_3l_pentabox_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_massless_3l_pentabox", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
define_graph --name=massless_3l_pentabox -ve=[(1000,1,6),(1000,6,7),(1000,7,2),(1000,2,1),(1000,7,8),(1000,8,9),(1000,9,6),(1000,8,3),(1000,3,4),(1000,4,5),(1000,5,9)] -ee=[("in",1),("in",2),("in",3),("in",4),("out",5)]
output {output_path} --overwrite_output --yaml_only"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def scalar_3L_6P_topology_A_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_scalar_3L_6P_topology_A", False).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model scalars;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'scalar_3L_6P_topology_A.dot')}
set target_omega 2.0
set panic_on_fail True
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def physical_3L_6photons_topology_A_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_physical_3L_6photons_topology_A", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'physical_3L_6photons_topology_A.dot')}
output {output_path} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def physical_2L_6photons_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_physical_2L_6photons", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'physical_2L_6photons.dot')}
output {output_path} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def physical_1L_6photons_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_physical_1L_6photons", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'physical_1L_6photons.dot')}
output {output_path} --overwrite_output"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def physical_1L_2A_final_4H_top_internal_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_physical_1L_2A_final_4H_top_internal", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'physical_1L_2A_final_4H_top_internal.dot')}
output {output_path} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def top_bubble_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_top_bubble", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'top_bubble.dot')}
output {output_path} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def hairy_glue_box_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     "TEST_AMPLITUDE_hairy_glue_box", True).joinpath("GL_OUTPUT")
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'hairy_glue_box.dot')}
output {output_path} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def ta_ta_tree_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     pjoin('trees', 'ta_ta'), True)
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(output_path, 'tree_amplitude_1_ta_ta.yaml')} --format yaml
output {output_path.joinpath("GL_OUTPUT")} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def th_th_tree_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     pjoin('trees', 'th_th'), True)
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(output_path, 'tree_amplitude_1_th_th.yaml')} --format yaml
output {output_path.joinpath("GL_OUTPUT")} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def t_ta_tree_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     pjoin('trees', 't_ta'), True)
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(output_path, 'tree_amplitude_1_t_ta.yaml')} --format yaml
output {output_path.joinpath("GL_OUTPUT")} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def hh_ttxaa_tree_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     pjoin('trees', 'hh_ttxaa'), True)
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(output_path, 'tree_amplitude_1_hh_ttxaa.yaml')} --format yaml
output {output_path.joinpath("GL_OUTPUT")} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def h_ttxaah_tree_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     pjoin('trees', 'h_ttxaah'), True)
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(output_path, 'tree_amplitude_1_h_ttxaah.yaml')} --format yaml
output {output_path.joinpath("GL_OUTPUT")} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def aa_aahhttx_tree_export(tmpdir_factory: pytest.TempPathFactory) -> Path:
    gloop = get_gamma_loop_interpreter()
    # Specify "True" below for a pytest designed to generate input for a rust test.
    output_path = get_test_directory(tmpdir_factory,
                                     pjoin('trees', 'aa_aahhttx'), True)
    gloop.run(CommandList.from_string(
        f"""import_model sm-full;
import_graphs {pjoin(output_path, 'tree_amplitude_1_aa_aahhttx.yaml')} --format yaml
output {output_path.joinpath("GL_OUTPUT")} --overwrite_output --yaml_only -exp -ef file"""))
    return output_path


@measure_fixture_setup_time(scope="session")
def compile_rust_tests() -> Path | None:

    # If you want to bypass the "manual" compilation of the rust tests, then uncomment the line below
    return None

    cmd_list = ['cargo', 'build', '--release', '--target-dir', os.path.normpath(os.path.join(
        GL_PATH, os.path.pardir, os.path.pardir, 'rust_test_binaries')), '--features=binary', '--no-default-features',
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
