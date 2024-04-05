import logging
import os
from pathlib import Path
from subprocess import Popen, PIPE
import gammaloop.misc.common
import gammaloop.interface.gammaloop_interface as gl_interface
from gammaloop.misc.common import GL_PATH, logger
import yaml

pjoin = os.path.join

RESOURCES_PATH = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), 'test_data')


def get_gamma_loop_interpreter() -> gl_interface.GammaLoop:
    gloop = gl_interface.GammaLoop()
    gammaloop.misc.common.GL_DEBUG = True
    gammaloop.misc.common.GL_CONSOLE_HANDLER.setLevel(logging.CRITICAL)
    return gloop


def run_rust_test(rust_tests_binary: Path, output_path: Path, test_name: str) -> bool:

    new_env: dict[str, str] = os.environ.copy()
    new_env['PYTEST_OUTPUT_PATH_FOR_RUST'] = str(output_path)
    if 'SYMBOLICA_LICENSE' not in new_env:
        new_env['SYMBOLICA_LICENSE'] = 'GAMMALOOP_USER'
    process = Popen(['cargo', 'test', f'pytest_{test_name}', '--features=binary,fail-on-warnings', '--no-default-features', '--release', '--target-dir', os.path.join(GL_PATH, os.path.pardir, os.path.pardir, 'rust_test_binaries'), '--', '--test-threads=1', '--ignored',
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

    process = Popen(['make', 'feynman_diagrams.pdf'],
                    cwd=drawing_path, stdout=PIPE, stderr=PIPE)
    output, error = process.communicate()
    if process.returncode != 0 or not os.path.isfile(pjoin(drawing_path, 'feynman_diagrams.pdf')):
        logger.info("DRAWING TEST STDOUT:\n%s", output.decode("utf-8"))
        logger.info("DRAWING TEST STDERR:\n%s", error.decode("utf-8"))
        if os.path.isfile(pjoin(drawing_path, 'compilation.log')):
            with open(pjoin(drawing_path, 'compilation.log'), 'r', encoding="utf-8") as f:
                compilation_log = f.read()
        else:
            compilation_log = "No compilation log found."
        logger.info("DRAWING TEST compilation log:\n%s", compilation_log)
        return False
    logger.debug("DRAWING TEST STDOUT:\n%s", output.decode("utf-8"))
    logger.debug("DRAWING TEST STDERR:\n%s", error.decode("utf-8"))
    return True


def check_integration_result(imag_phase: bool, target: float, process_path: Path):
    max_mc_error_dif = 5.0
    max_rel_error_dif = 0.01
    max_percent_error = 0.01

    # copy settings
    f = open(os.path.join(process_path, 'runs', 'run.yaml'), "r").read()
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
