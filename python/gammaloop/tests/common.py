import logging
import os
from subprocess import Popen, PIPE
import gammaloop
from gammaloop.misc.common import GL_PATH, logger
pjoin = os.path.join

RESOURCES_PATH = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), 'test_data')


def get_gamma_loop_interpreter() -> gammaloop.GammaLoop:
    gloop = gammaloop.GammaLoop()
    gammaloop.GL_DEBUG = True
    gammaloop.GL_CONSOLE_HANDLER.setLevel(logging.CRITICAL)
    return gloop


def run_rust_test(rust_tests_binary, output_path, test_name) -> bool:

    new_env = os.environ.copy()
    new_env['PYTEST_OUTPUT_PATH_FOR_RUST'] = output_path
    process = Popen([rust_tests_binary, test_name, '--test-threads=1', '--ignored',
                    '--nocapture'], cwd=GL_PATH, stdout=PIPE, stderr=PIPE, env=new_env)
    output, error = process.communicate()
    if process.returncode != 0:
        logger.info("RUST TEST STDOUT:\n%s", output.decode("utf-8"))
        logger.info("RUST TEST STDERR:\n%s", error.decode("utf-8"))
        return False
    logger.debug("RUST TEST STDOUT:\n%s", output.decode("utf-8"))
    logger.debug("RUST TEST STDERR:\n%s", error.decode("utf-8"))
    return True
