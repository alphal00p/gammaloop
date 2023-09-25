import logging
import os
import gammaloop

pjoin = os.path.join

RESOURCES_PATH = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), 'test_data')


def get_gamma_loop_interpreter() -> gammaloop.GammaLoop:
    gloop = gammaloop.GammaLoop()
    gammaloop.GL_DEBUG = True
    gammaloop.GL_CONSOLE_HANDLER.setLevel(logging.CRITICAL)
    return gloop
