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
    # gloop.run(gl_interface.CommandList.from_string("set compile_cff False"))
    # gloop.run(gl_interface.CommandList.from_string("set load_compiled_cff False"))
    # gloop.run(gl_interface.CommandList.from_string("set export_settings.numerator_settings.eval_settings.compile_options.subtype 'NotCompiled'"))
    # gloop.run(gl_interface.CommandList.from_string("set export_settings.gammaloop_compile_options.inline_asm False"))

    gammaloop.misc.common.GL_DEBUG = True
    gammaloop.misc.common.GL_CONSOLE_HANDLER.setLevel(logging.CRITICAL)
    return gloop


def get_gamma_loop_interpreter_no_compilation() -> gl_interface.GammaLoop:
    gloop = gl_interface.GammaLoop()
    gloop.run(gl_interface.CommandList.from_string("set compile_cff False"))
    gloop.run(gl_interface.CommandList.from_string(
        "set load_compiled_cff False"))
    gloop.run(gl_interface.CommandList.from_string(
        "set export_settings.numerator_settings.eval_settings.compile_options.subtype 'NotCompiled'"))
    # gloop.run(gl_interface.CommandList.from_string("set export_settings.gammaloop_compile_options.inline_asm False"))

    gammaloop.misc.common.GL_DEBUG = True
    gammaloop.misc.common.GL_CONSOLE_HANDLER.setLevel(logging.CRITICAL)
    return gloop


def run_rust_test(rust_tests_binary: Path | None, output_path: Path, test_name: str) -> bool:
    new_env: dict[str, str] = os.environ.copy()
    match new_env.get('PYTEST_OUTPUT_PATH_FOR_RUST', None):
        case "TMP":
            if os.path.isfile(output_path):
                new_env['PYTEST_OUTPUT_PATH_FOR_RUST'] = os.path.normpath(
                    os.path.join(os.path.dirname(output_path), os.path.pardir))
            elif os.path.isdir(output_path):
                new_env['PYTEST_OUTPUT_PATH_FOR_RUST'] = os.path.normpath(
                    os.path.join(output_path, os.path.pardir, os.path.pardir))
            else:
                raise gammaloop.misc.common.GammaLoopError(
                    f"Could not find specified output path '{output_path}' to run the rust test on.")
        case None:
            new_env['PYTEST_OUTPUT_PATH_FOR_RUST'] = os.path.normpath(os.path.join(
                GL_PATH, os.path.pardir, os.path.pardir, 'src', 'test_resources'))
        case _:
            pass

    if 'SYMBOLICA_LICENSE' not in new_env:
        new_env['SYMBOLICA_LICENSE'] = 'GAMMALOOP_USER'

    if rust_tests_binary is None:
        cmd_list = ['cargo', 'test', f'pytest_{test_name}', '--features=binary,fail-on-warnings', '--no-default-features',
                    '--release', '--target-dir', os.path.normpath(os.path.join(
                        GL_PATH, os.path.pardir, os.path.pardir, 'rust_test_binaries')),
                    '--', '--test-threads=1', '--nocapture']
    else:
        cmd_list = [f'{rust_tests_binary}', f'pytest_{test_name}',
                    '--test-threads=1', '--nocapture']
    logger.debug("Running rust test with command: %s", " ".join(cmd_list))
    process = Popen(cmd_list, cwd=GL_PATH, stdout=PIPE,
                    stderr=PIPE, env=new_env)
    output, error = process.communicate()
    if process.returncode != 0:
        logger.info("RUST TEST STDOUT:\n%s", output.decode("utf-8"))
        logger.info("RUST TEST STDERR:\n%s", error.decode("utf-8"))
        return False
    logger.debug("RUST TEST STDOUT:\n%s", output.decode("utf-8"))
    logger.debug("RUST TEST STDERR:\n%s", error.decode("utf-8"))
    return True


def run_drawing(drawing_root_path: str, drawing_modes=['dot','feynmp']) -> bool:

    for drawing_mode in drawing_modes:
        drawing_path = os.path.join(drawing_root_path, drawing_mode)
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


def check_integration_result(target: float, process_path: Path, max_mc_error_diff=5.0, max_rel_error_diff=0.01, max_percent_error=0.01, imag_phase=False):

    with open(os.path.join(process_path, 'runs', 'run.yaml'), "r") as f:
        run_yaml = yaml.safe_load(f.read())

    index_for_check = 1 if imag_phase else 0
    res_for_check = run_yaml["result"][index_for_check]
    error_for_check = run_yaml["error"][index_for_check]

    absolute_difference = abs(res_for_check - target)
    assert absolute_difference <= max_mc_error_diff * error_for_check

    if abs(target) > 0.:
        relative_difference = absolute_difference / abs(target)
        if relative_difference >= max_rel_error_diff:
            logger.debug("Failing test: relative_difference = %.16e (>= max_rel_error_diff = %.16e)",
                         relative_difference, max_rel_error_diff)
        assert relative_difference < max_rel_error_diff
        if error_for_check >= max_percent_error * abs(target):
            logger.debug("Failing test: error_for_check = %.16e (>= max_percent_error * abs(target) = %.16e)",
                         error_for_check, max_percent_error * abs(target))
        assert error_for_check < max_percent_error * abs(target)


def check_inspect_result(inspect_result: complex, target: complex, max_relative_diff=1.0e-12):
    for target, diff in [
        (abs(target.real), abs(inspect_result.real - target.real)),
        (abs(target.imag), abs(inspect_result.imag - target.imag))
    ]:
        if abs(target) > 0.:
            if diff >= max_relative_diff * target:
                logger.debug("Failing test: diff = %.16e (>= max_relative_diff * target = %.16e)",
                             diff, max_relative_diff * target)
            assert diff < max_relative_diff * target
        else:
            if diff >= max_relative_diff:
                logger.debug("Failing test: diff = %.16e (>= max_relative_diff = %.16e)",
                             diff, max_relative_diff)
            assert diff < max_relative_diff
