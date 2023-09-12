from enum import StrEnum

import logging
import logging.handlers
import os
import sys
import symbolica as sb
import gammaloop.misc.common as common
import yaml

class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data):
        return True

def verbose_yaml_dump(data):
    return yaml.dump(data, Dumper=NoAliasDumper, default_flow_style=False, sort_keys=False)

class Colour(StrEnum):
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

def parse_python_expression(expr: str | None) -> sb.Expression | None:
    if expr is None:
        return None
    santized_expr = expr.replace('**', '^')\
        .replace('cmath.sqrt', 'sqrt')\
        .replace('cmath.pi', 'pi')\
        .replace('math.sqrt', 'sqrt')\
        .replace('math.pi', 'pi')
    try:
        return sb.Expression.parse(santized_expr)
    except Exception as e:
        common.logger.critical("Symbolica failed to parse expression:\n{}\nwith exception:\n{}".format(santized_expr, e))
        return None

def expression_to_string(expr: sb.Expression | None) -> str | None:
    if expr is None:
        return None
    try:
        return expr.pretty_str(
            terms_on_new_line = False,
            color_top_level_sum = False,
            color_builtin_functions = False,
            print_finite_field = False,
            explicit_rational_polynomial = False,
            number_thousands_separator = None,
            multiplication_operator = '*',
            square_brackets_for_function = False,
            num_exp_as_superscript = False,
            latex = False)
    except Exception as e:
        common.logger.critical("Symbolica failed to cast expression to string:\n{}\nwith exception:\n{}".format(expr, e))
        return None

def setup_logging() -> logging.StreamHandler:
    console_format = '[{} %(asctime)s {}] @{}%(name)s{} %(levelname)s: %(message)s'.format(
        Colour.GREEN, Colour.END, Colour.BLUE, Colour.END)
    file_format = '[%(asctime)s] %(name)s %(levelname)s: %(message)s'

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(console_format))

    logging.getLogger().handlers = []
    logging.getLogger().addHandler(console_handler)

    if 'GL_ENABLE_FILE_HANDLERS' in os.environ and os.environ['GL_ENABLE_FILE_HANDLERS'].upper() != 'FALSE':
        log_file_name = 'gammaloop_debug.log'
        file_handler = logging.FileHandler(log_file_name, mode='w')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(file_format))
        logging.getLogger().addHandler(file_handler)

        error_file_name = 'gammaloop_error.log'
        error_file_handler = logging.FileHandler(error_file_name, mode='w')
        error_file_handler.setLevel(logging.ERROR)
        error_file_handler.setFormatter(logging.Formatter(file_format))
        logging.getLogger().addHandler(error_file_handler)

    logging.getLogger().setLevel(logging.DEBUG)

    return console_handler
