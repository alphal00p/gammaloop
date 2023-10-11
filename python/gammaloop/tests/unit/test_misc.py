import pytest
from subprocess import Popen, PIPE
from gammaloop.misc.common import GL_PATH, logger


class TestCode:

    @pytest.mark.typecheck
    def test_pyright(self):
        process = Popen(['pyright'], cwd=GL_PATH, stdout=PIPE, stderr=PIPE)
        output, error = process.communicate()
        if process.returncode != 0:
            logger.info("PYRIGHT TEST STDOUT:\n%s", output.decode("utf-8"))
            logger.info("PYRIGHT TEST STDERR:\n%s", error.decode("utf-8"))
        logger.debug("PYRIGHT TEST STDOUT:\n%s", output.decode("utf-8"))
        logger.debug("PYRIGHT TEST STDERR:\n%s", error.decode("utf-8"))
        assert process.returncode == 0
