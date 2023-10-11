import pytest
from subprocess import Popen, PIPE
from gammaloop.misc.common import GL_PATH, logger


class TestCode:

    @pytest.mark.typecheck
    def test_pyright(self):
        process = Popen(['pyright'], cwd=GL_PATH, stdout=PIPE, stderr=PIPE)
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
