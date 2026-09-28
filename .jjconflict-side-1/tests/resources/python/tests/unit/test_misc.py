import pytest
import json
import os
import sys
from subprocess import Popen, PIPE
from pprint import pformat
from gammaloop.misc.common import GL_PATH, logger, GammaLoopError
from gammaloop import check_gammaloop_dependencies


class TestCode:

    @pytest.mark.codecheck
    def test_pyright(self):
        check_gammaloop_dependencies()
        custom_env = dict(os.environ)
        custom_env["PYTHONPATH"] = ':'.join(
            sys.path+([custom_env['PYTHONPATH'],] if 'PYTHONPATH' in custom_env else []))
        custom_env["PYRIGHT_PYTHON_FORCE_VERSION"] = "latest"
        process = Popen(['pyright', '--warnings'],
                        cwd=GL_PATH, stdout=PIPE, stderr=PIPE, env=custom_env)
        output, error = process.communicate()
        stdout = output.decode("utf-8")
        stderr = error.decode("utf-8")
        success = (process.returncode ==
                   0 and '0 error' in stdout and '0 warning' in stdout)
        if not success:
            logger.info("PYRIGHT TEST STDOUT:\n%s", stdout)
            logger.info("PYRIGHT TEST STDERR:\n%s", stderr)
        else:
            logger.debug("PYRIGHT TEST STDOUT:\n%s", stdout)
            logger.debug("PYRIGHT TEST STDERR:\n%s", stderr)
        assert process.returncode == 0
        assert '0 error' in stdout
        assert '0 warning' in stdout

    @pytest.mark.codecheck
    def test_clippy(self):
        process = Popen(['cargo', 'clippy', '--message-format=json'],
                        cwd=GL_PATH, stdout=PIPE, stderr=PIPE)
        logger.debug("Running clippy diagnostics...")
        output, err = process.communicate()
        if process.returncode != 0:
            raise GammaLoopError(
                "Failed to run clippy diagnostics. Error:\n" + err.decode("utf-8"))

        compiler_artifact = output.decode("utf-8")
        warning_msg = None
        i_msg = 0
        found_warning = None
        for json_line in reversed(compiler_artifact.split("\n")):
            if json_line == "":
                continue
            try:
                json_obj = json.loads(json_line)
            except json.decoder.JSONDecodeError:
                continue
            if 'message' in json_obj:
                json_msg = json_obj['message']
                if "level" in json_msg and json_msg["level"] in ["warning", "error"]:
                    if warning_msg is None:
                        found_warning = f"Clippy issued at least one warning : {'Unknown' if 'message' not in json_msg else json_msg['message']}"
                        warning_msg = found_warning
                        warning_msg += "\nAll clippy warnings and errors:"
                    else:
                        i_msg += 1
                        for category in ['spans', 'children']:
                            if category in json_msg:
                                del json_msg[category]
                        warning_msg += f"\n\n>>> #{i_msg}:\n\n{pformat(json_msg)}"
        if warning_msg is not None:
            logger.critical(warning_msg)

        assert found_warning is None
