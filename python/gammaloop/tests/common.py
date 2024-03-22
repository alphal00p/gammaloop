import logging
import os
from pathlib import Path
from subprocess import Popen, PIPE
import gammaloop.misc.common
import gammaloop.interface.gammaloop_interface as gl_interface
from gammaloop.misc.common import GL_PATH, logger
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
