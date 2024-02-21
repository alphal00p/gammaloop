import logging
import os
from pathlib import Path
from subprocess import Popen, PIPE
import gammaloop
from gammaloop.misc.common import GL_PATH, logger
import shutil
import yaml

pjoin = os.path.join

RESOURCES_PATH = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), 'test_data')


def get_gamma_loop_interpreter() -> gammaloop.GammaLoop:
    gloop = gammaloop.GammaLoop()
    gammaloop.GL_DEBUG = True
    gammaloop.GL_CONSOLE_HANDLER.setLevel(logging.CRITICAL)
    return gloop


def run_rust_test(rust_tests_binary: Path, output_path: Path, test_name: str) -> bool:

    new_env: dict[str, str] = os.environ.copy()
    new_env['PYTEST_OUTPUT_PATH_FOR_RUST'] = str(output_path)
    process = Popen([rust_tests_binary, f'pytest_{test_name}', '--test-threads=1', '--ignored',
                    '--nocapture'], cwd=GL_PATH, stdout=PIPE, stderr=PIPE, env=new_env)
    output, error = process.communicate()
    if process.returncode != 0:
        logger.info("RUST TEST STDOUT:\n%s", output.decode("utf-8"))
        logger.info("RUST TEST STDERR:\n%s", error.decode("utf-8"))
        return False
    logger.debug("RUST TEST STDOUT:\n%s", output.decode("utf-8"))
    logger.debug("RUST TEST STDERR:\n%s", error.decode("utf-8"))
    return True


def run_drawing(drawing_path: str) -> bool:

    process = Popen(['make', f'feynman_diagrams.pdf'],
                    cwd=drawing_path, stdout=PIPE, stderr=PIPE)
    output, error = process.communicate()
    if process.returncode != 0 or not os.path.isfile(pjoin(drawing_path, 'feynman_diagrams.pdf')):
        logger.info("DRAWING TEST STDOUT:\n%s", output.decode("utf-8"))
        logger.info("DRAWING TEST STDERR:\n%s", error.decode("utf-8"))
        return False
    logger.debug("DRAWING TEST STDOUT:\n%s", output.decode("utf-8"))
    logger.debug("DRAWING TEST STDERR:\n%s", error.decode("utf-8"))
    return True


def scalar_euclidean_integration_test(graph: str, imag_phase: bool, target: float):
    max_mc_error_dif = 5.0
    max_rel_error_dif = 0.01
    max_percent_error = 0.01

    workspace = os.path.join(RESOURCES_PATH, 'test_{}'.format(graph))
    graph_path = os.path.join(
        RESOURCES_PATH, 'qgraf_outputs', '{}.py'.format(graph))
    test_config_path = os.path.join(
        RESOURCES_PATH, 'test_configs', '{}.yaml'.format(graph))

    # clean up output if it is there
    if os.path.exists(workspace):
        shutil.rmtree(workspace)

    gl = get_gamma_loop_interpreter()

#     # create the output needed for the test
    command_list = gammaloop.CommandList.from_string(
        "import_model scalars-full")
    command_list.add_command(
        f"import_graphs {graph_path} --format=qgraph")
    command_list.add_command(f"output {workspace}")

    gl.run(command_list)

    # copy settings
    shutil.copyfile(test_config_path,
                    os.path.join(workspace, 'cards', 'run_card.yaml'))

    command_list = gammaloop.CommandList.from_string(
        "launch {}".format(workspace))
    command_list.add_command("integrate {}".format(graph))

    gl.run(command_list)

    f = open(os.path.join(workspace, 'runs', 'run.yaml'), "r").read()
    run_yaml = yaml.safe_load(f)

    res_for_check = run_yaml["result"][0]
    error_for_check = run_yaml["error"][0]
    zero_res = run_yaml["result"][1]
    if imag_phase:
        res_for_check = run_yaml["result"][1]
        error_for_check = run_yaml["error"][1]
        zero_res = run_yaml["result"][0]

    assert zero_res == 0.0

    absolute_difference = abs(res_for_check - target)
    relative_difference = absolute_difference / abs(target)

    assert absolute_difference < max_mc_error_dif * error_for_check
    assert relative_difference < max_rel_error_dif
    assert error_for_check < max_percent_error * abs(target)

    shutil.rmtree(workspace)
