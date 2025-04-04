from enum import StrEnum
from typing import TextIO, Any, Callable
from collections import deque, defaultdict
import logging
import logging.handlers
import os
import sys

import yaml
import time
import re
from pprint import pformat

import gammaloop.misc.common as common
from gammaloop.misc import LOGGING_PREFIX_FORMAT
import gammaloop._gammaloop as gl_rust

from symbolica import Expression as SBE  # type: ignore # pylint: disable=import-error # nopep8
import symbolica as sb  # pylint: disable=import-error # type: ignore # nopep8


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any):
        return True


def verbose_yaml_dump(data: Any):
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
    GRAY = '\033[21m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def evaluate_graph_overall_factor(overall_factor: str) -> SBE:
    e = parse_python_expression(overall_factor)
    assert e is not None, f"Could not parse overall factor: '{overall_factor}'"  # nopep8
    for header in ["AutG",
                   "CouplingsMultiplicity",
                   "InternalFermionLoopSign",
                   "ExternalFermionOrderingSign",
                   "AntiFermionSpinSumSign",
                   "NumeratorIndependentSymmetryGrouping"]:
        e = e.replace_all(sb.Expression.parse(
            f"{header}(x_)"), sb.Expression.parse("x_"))
    e = e.replace_all(
        sb.Expression.parse("NumeratorDependentGrouping(GraphId_,ratio_,GraphSymmetryFactor_)"), sb.Expression.parse("ratio_*GraphSymmetryFactor_"))
    return e.expand().collect_num()


def recursive_replace_all(expr: sb.Expression, replace: Callable[[SBE], SBE], max_recursion: int = 1000) -> SBE:
    n_iter = 0
    while n_iter < max_recursion:
        replaced_expr = replace(expr)
        if replaced_expr == expr:
            break
        expr = replaced_expr
        n_iter += 1
    if n_iter == max_recursion:
        raise common.GammaLoopError(
            "Symbolica @%s failed to fully replace expression:\n%s\nwith %d recursions.", sb.__file__, expr, max_recursion)
    return expr


def evaluate_symbolica_expression(expr: sb.Expression, evaluation_variables: dict[SBE, complex], evaluation_functions: dict[Callable[[SBE], SBE], Callable[[list[complex]], complex]]) -> complex | None:
    try:
        return evaluate_symbolica_expression_safe(expr, evaluation_variables, evaluation_functions)
    except common.GammaLoopError:
        return None


def evaluate_symbolica_expression_safe(expr: sb.Expression, evaluation_variables: dict[SBE, complex], evaluation_functions: dict[Callable[[SBE], SBE], Callable[[list[complex]], complex]]) -> complex:
    try:
        res: complex = expr.evaluate_complex(  # type: ignore
            evaluation_variables, evaluation_functions)  # type: ignore
        # Small adjustment to avoid havin -0. in either the real or imaginary part
        return complex(0 if abs(res.real) == 0. else res.real,  # type: ignore
                       0 if abs(res.imag) == 0. else res.imag)  # type: ignore
    except BaseException as e:
        raise common.GammaLoopError(
            "Symbolica (@%s) failed to evaluate expression:\n%s\nwith exception:\n%s.\nVariables:\n%s\nFunctions:\n%s",
            sb.__file__,
            expression_to_string(expr), e,
            pformat({expression_to_string(k): v for k,
                    v in evaluation_variables.items()}),
            pformat([expression_to_string(k()).replace('()', '')  # type: ignore
                    for k in evaluation_functions.keys()])  # type: ignore
        )


def parse_python_expression(expr: str | None) -> sb.Expression | None:
    if expr is None:
        return None
    try:
        return parse_python_expression_safe(expr)
    except Exception as exception:  # pylint: disable=broad-except
        common.logger.critical("%s", exception)
        return None


def replace_pseudo_floats(expression: str) -> str:

    def rationalize_float(fl: re.Match[str]) -> sb.Expression:
        fl_eval: float = eval(fl.group())
        # Work around a bug for 0.0 in symbolica
        rationalized_fl = SBE.num(
            fl_eval, 1e-13) if fl_eval != 0. else SBE.num(0)  # type: ignore
        rationalized_fl_eval: float = eval(str(rationalized_fl)+'.')
        if common.GammaLoopWarning.FloatInExpression not in common.GL_WARNINGS_ISSUED:
            common.GL_WARNINGS_ISSUED.add(
                common.GammaLoopWarning.FloatInExpression)
            common.logger.warning(
                f"Expression contains the following floating point values in the expression: {fl.group()} (mapped to {rationalized_fl})")
            common.logger.warning(
                "It is typically best for the UFO model to be adjusted and declare these floats as internal constant parameters instead, or have them written as a fraction.")
        valid_cast = True
        if abs(rationalized_fl_eval-fl_eval) != 0.:
            if abs(fl_eval) != 0.:
                if abs(rationalized_fl_eval-fl_eval)/abs(fl_eval) > 1.0e-12:
                    valid_cast = False
            else:
                if abs(rationalized_fl_eval-fl_eval) > 1.0e-12:
                    valid_cast = False
        if not valid_cast:
            raise common.GammaLoopError(
                "The float value %s in the expression could not be cast to a rational number %s = %.16e", fl.group(), rationalized_fl, rationalized_fl_eval)

        return rationalized_fl

    modified_expression = re.sub(
        r'(\d+\.\d+e[+-]?\d+)', lambda x: f'({rationalize_float(x)})', expression)
    modified_expression = re.sub(
        r'(\d+\.\d+)', lambda x: f'({rationalize_float(x)})', expression)
    modified_expression = re.sub(r'(\d+)\.', r'\1', modified_expression)

    return modified_expression


def replace_to_sqrt(expr: sb.Expression) -> sb.Expression:
    expr = expr.replace_all(sb.Expression.parse(
        'x__^(1/2)'), sb.Expression.parse('sqrt(x__)'))
    expr = expr.replace_all(sb.Expression.parse(
        'x__^(-1/2)'), sb.Expression.parse('(sqrt(x__))^-1'))
    expr = expr.replace_all(sb.Expression.parse(
        'x__^(1/4)'), sb.Expression.parse('sqrt(sqrt(x__))'))
    expr = expr.replace_all(sb.Expression.parse(
        'x__^(-1/4)'), sb.Expression.parse('(sqrt(sqrt(x__)))^-1'))
    str_expr = expression_to_string(expr)
    if str_expr is None or re.match(r'\^\(\d+/\d+\)', str_expr):
        raise common.GammaLoopError(
            "Expoentiation with real arguments not supported in model expressions: %s", str_expr)
    return expr


def replace_from_sqrt(expr: sb.Expression) -> sb.Expression:
    expr = expr.replace_all(sb.Expression.parse(
        'sqrt(x__)'), sb.Expression.parse('x__^(1/2)'))
    str_expr = expression_to_string(expr)
    if str_expr is None or re.match(r'\^\(\d+/\d+\)', str_expr):
        raise common.GammaLoopError(
            "Exponentiation with real arguments not supported in model expressions: %s", str_expr)
    return expr


def parse_python_expression_safe(expr: str) -> sb.Expression:

    sanitized_expr = expr.replace('**', '^')\
        .replace('cmath.sqrt', 'sqrt')\
        .replace('cmath.pi', 'pi')\
        .replace('math.sqrt', 'sqrt')\
        .replace('math.pi', 'pi')
    sanitized_expr = replace_pseudo_floats(sanitized_expr)
    try:
        sb_expr = sb.Expression.parse(sanitized_expr)
        # No longer needed since we automatically include a `complex(x,y)` function in the function map.
        # sb_expr = recursive_replace_all(sb_expr, lambda e: e.replace_all(
        #     SBE.parse('complex(x_,y_)'), SBE.parse('x_+I*y_')), max_recursion=1000)
    except Exception as exception:  # pylint: disable=broad-except
        raise common.GammaLoopError(
            "Symbolica (@%s) failed to parse expression:\n%s\nwith exception:\n%s", sb.__file__, sanitized_expr, exception)

    # sb_expr_processed = replace_to_sqrt(sb_expr)
    sb_expr_processed = replace_from_sqrt(sb_expr)

    return sb_expr_processed


def expression_to_string(expr: sb.Expression | None, canonical=False) -> str | None:
    if expr is None:
        return None
    try:
        return expression_to_string_safe(expr, canonical)
    except Exception as exception:  # pylint: disable=broad-except
        common.logger.critical("%s", exception)
        return None


def expression_to_string_safe(expr: sb.Expression, canonical=False) -> str:
    try:
        expr_str = expr.format(
            terms_on_new_line=False,
            color_top_level_sum=False,
            color_builtin_symbols=False,
            print_finite_field=False,
            explicit_rational_polynomial=False,
            number_thousands_separator=None,
            multiplication_operator='*',
            square_brackets_for_function=False,
            num_exp_as_superscript=False,
            latex=False)
        if canonical:
            return gl_rust.atom_to_canonical_string(expr_str)
        else:
            return expr_str
    except Exception as exception:  # pylint: disable=broad-except
        raise common.GammaLoopError(
            "Symbolica (@%s)failed to cast expression to string:\n%s\nwith exception:\n%s", sb.__file__, expr, exception)


class GammaLoopCustomFormatter(logging.Formatter):
    """Logging colored formatter"""

    def __init__(self, fmt: str, datefmt: str | None = None):
        super().__init__(datefmt=datefmt)
        self.fmt = fmt
        self.datefmt = datefmt

    def format(self, record: logging.LogRecord) -> str:
        formatter = logging.Formatter(self.fmt, self.datefmt)
        if record.levelno != logging.DEBUG:
            if record.name.startswith('_gammaloop'):
                record.name = f"rust.{record.name[11:]}"
            if len(record.name) > 20:
                record.name = f"{record.name[:17]}..."
            record.name = f"{record.name:20}"
        match record.levelno:
            case logging.DEBUG:
                record.levelname = f"{Colour.GRAY}{record.levelname:8}{Colour.END}"  # nopep8
            case logging.INFO:
                record.levelname = f"{record.levelname:8}"  # nopep8
            case logging.WARNING:
                record.levelname = f"{Colour.YELLOW}{record.levelname:8}{Colour.END}"  # nopep8
            case logging.ERROR:
                record.levelname = f"{Colour.RED}{record.levelname:8}{Colour.END}"  # nopep8
            case logging.CRITICAL:
                record.levelname = f"{Colour.RED}{Colour.BOLD}{record.levelname:8}{Colour.END}"  # nopep8
            case _:
                record.levelname = f"{record.levelname:8}"
        record.asctime = self.formatTime(record, self.datefmt)
        return formatter.format(record)

    # def formatTime(self, record: logging.LogRecord, datefmt: str | None = None) -> str:
    #     if datefmt is None:
    #         return time.strftime('%Y-%m-%d %H:%M:%S.%f', time.localtime(record.created))
    #     else:
    #         return time.strftime(datefmt, time.localtime(record.created))


def setup_logging() -> logging.StreamHandler[TextIO]:
    match LOGGING_PREFIX_FORMAT:
        case 'none':
            console_format = f'%(message)s'
            time_format = "%H:%M:%S"
        case 'min':
            console_format = f'%(levelname)s: %(message)s'
            time_format = "%H:%M:%S"
        case 'short':
            console_format = f'[{Colour.GREEN}%(asctime)s{Colour.END}] %(levelname)s: %(message)s'  # nopep8
            time_format = "%H:%M:%S"
        case 'long':
            console_format = f'[{Colour.GREEN}%(asctime)s.%(msecs)03d{Colour.END}] @{Colour.BLUE}%(name)s{Colour.END} %(levelname)s: %(message)s'  # nopep8
            time_format = '%Y-%m-%d %H:%M:%S'
        case _:
            raise common.GammaLoopError(
                "Invalid LOGGING_PREFIX_FORMAT: %s", LOGGING_PREFIX_FORMAT)
    file_format = '[%(asctime)s] %(name)s %(levelname)s: %(message)s'
    console_formatter = GammaLoopCustomFormatter(
        console_format, datefmt=time_format)
    file_formatter = GammaLoopCustomFormatter(file_format)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)

    console_handler.setFormatter(console_formatter)
    logging.getLogger().handlers = []
    logging.getLogger().addHandler(console_handler)

    if 'GL_ENABLE_FILE_HANDLERS' in os.environ and os.environ['GL_ENABLE_FILE_HANDLERS'].upper() != 'FALSE':
        log_file_name = 'gammaloop_debug.log'
        file_handler = logging.FileHandler(log_file_name, mode='w')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logging.getLogger().addHandler(file_handler)

        error_file_name = 'gammaloop_error.log'
        error_file_handler = logging.FileHandler(error_file_name, mode='w')
        error_file_handler.setLevel(logging.ERROR)
        error_file_handler.setFormatter(file_formatter)
        logging.getLogger().addHandler(error_file_handler)

    logging.getLogger().setLevel(logging.DEBUG)

    gl_rust.setup_rust_logging(logging.getLevelName(
        console_handler.level), LOGGING_PREFIX_FORMAT)

    return console_handler


def remove_duplicates(seq: list[Any]):
    seen: set[Any] = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


### useful graph algorithms ###
def generate_spanning_trees(result: list[tuple[int, ...]], edge_map: list[tuple[int, int]], adjacency_map: dict[int, list[tuple[tuple[int, tuple[int, int]], int]]], current_tree: set[int],
                            accumulated_edge_sequence: list[int], excluded_edges: list[int], seen_edge_sequences: set[tuple[int, ...]] | None, target_number_of_st_to_find: int | None = None, follow: bool = False) -> None:
    """Compute all spanning trees of a graph component. Disconnected graphs
    are supported: only the component connected to the vertex in `tree` is considered.
    """
    if target_number_of_st_to_find is not None and len(result) >= target_number_of_st_to_find:
        return
    if seen_edge_sequences is not None:
        sorted_sequence = tuple(sorted(accumulated_edge_sequence))
        if sorted_sequence in seen_edge_sequences:
            return
        seen_edge_sequences.add(sorted_sequence)

    # find all edges that connect the tree to a new node
    # This is the most computationally heavy function, it is written like this for performance reason
    # len(current_tree) > len(edge_map)-len(accumulated_edge_sequence):
    if len(current_tree) > len(edge_map)//2:
        connected_edges = [(i, e) for i, e in enumerate(
            edge_map) if i not in accumulated_edge_sequence and e[0] != e[1] and (e[0] in current_tree) ^ (e[1] in current_tree)]
    else:
        connected_edges = [e for v in current_tree for (
            e, v2) in adjacency_map[v] if v2 not in current_tree]

    edges_filtered = [(i, e)
                      for i, e in connected_edges if i not in excluded_edges]

    if len(edges_filtered) == 0:
        if len(connected_edges) == 0:
            # no more new edges, so we are done
            spanning_tree = tuple(sorted(accumulated_edge_sequence))
            if spanning_tree not in result:
                result.append(spanning_tree)
                if follow:
                    print(
                        f"Current length of spanning trees: {len(result)}", end='\r')

    else:
        # these edges can be added in any order, so only allow an order from lowest to highest
        for e_i, (i, edge_vertices) in enumerate(connected_edges):
            excluded_edges.extend(j for j, _ in connected_edges[:e_i])
            new_vertex = edge_vertices[0] if edge_vertices[1] in current_tree else edge_vertices[1]
            accumulated_edge_sequence.append(i)
            current_tree.add(new_vertex)
            generate_spanning_trees(
                result, edge_map, adjacency_map, current_tree, accumulated_edge_sequence,
                excluded_edges, seen_edge_sequences, target_number_of_st_to_find, follow)
            current_tree.remove(new_vertex)
            accumulated_edge_sequence.pop()
            excluded_edges = excluded_edges[:-e_i]


def find_longest_cycle(edge_map: list[tuple[int, int]]) -> list[int]:
    adjacency_map: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i, (u, v) in enumerate(edge_map):
        adjacency_map[u].append((v, i))
        # Add the reverse connection as it is an undirected graph
        adjacency_map[v].append((u, i))

    longest_cycle: list[int] = []

    def dfs(u: int, visited_edges: set[int], stack: list[int]):
        nonlocal longest_cycle
        stack.append(u)

        for v, edge_id in adjacency_map[u]:
            if edge_id in visited_edges:
                continue

            # A potential cycle is found
            if v in stack:
                cycle_start_index = stack.index(v)
                cycle = stack[cycle_start_index:]

                # Check if the cycle is simple (doesn't contain other cycles)
                if len(cycle) == len(set(cycle)):
                    if len(cycle) > len(longest_cycle):
                        longest_cycle = cycle

            else:
                # Continue DFS
                dfs(v, visited_edges | {edge_id}, stack.copy())

    for vertex in adjacency_map:
        dfs(vertex, set(), [])

    # Convert the cycle of vertices to a cycle of edges
    edge_cycle: list[int] = []
    for i in range(len(longest_cycle)):
        u, v = longest_cycle[i], longest_cycle[(i + 1) % len(longest_cycle)]
        edge_id: int = next(
            edge for edge in adjacency_map[u] if edge[0] == v)[1]
        edge_cycle.append(edge_id)

    return edge_cycle


def edges_cycle_to_vertex_cycle(edges: list[tuple[int, int]]) -> list[int]:
    # Check if the list is empty
    if not edges:
        return []

    # Special case: single unique edge forming a cycle
    unique_edges = set((min(src, dest), max(src, dest))
                       for src, dest in edges)
    if len(unique_edges) == 1:
        # Handle self-loop as a special case
        if edges[0][0] == edges[0][1]:
            return [edges[0][0]]
        return list(edges[0])

    # Initialize the ordered list of nodes with the source node of the first edge
    ordered_nodes = [edges[0][0]]

    # Initialize a set to keep track of visited edges
    visited_edges: set[tuple[int, int]] = set()

    # Start from the source node of the first edge
    current_node = edges[0][0]

    while True:
        next_node: int | None = None  # Define a default value for next_node
        for src, dest in edges:
            # Create edge tuple for checking visited status
            edge = (min(src, dest), max(src, dest))

            if edge in visited_edges:
                continue

            if src == current_node:
                next_node = dest
            elif dest == current_node:
                next_node = src

            if next_node is not None:
                # Mark edge as visited
                visited_edges.add(edge)

                # Break if the cycle is complete
                if next_node == ordered_nodes[0]:
                    return ordered_nodes

                # Add the next node to the ordered list
                ordered_nodes.append(next_node)

                # Move to the next node for the next iteration
                current_node = next_node
                break
        else:
            # Handle the case where next_node is not defined
            raise common.GammaLoopError(f"Invalid cycle {str(edges)}")


def create_adjacency_list(edge_map: list[tuple[int, tuple[int, int]]]) -> dict[int, set[int]]:
    adj_list: dict[int, set[int]] = defaultdict(set)
    for _, (u, v) in edge_map:  # pylint: disable=invalid-name
        adj_list[u].add(v)
        adj_list[v].add(u)
    return adj_list


def find_shortest_path(edge_map_list: list[tuple[int, int]], start: int, targets: list[int]) -> list[int]:
    edge_map = list(enumerate(edge_map_list))
    graph = create_adjacency_list(edge_map)
    visited: set[int] = set()
    # Queue for BFS initialized with the start node and parent
    queue: deque[tuple[int, None | int]] = deque([(start, None)])
    # Dictionary to hold the parent of each node
    parent: dict[int, None | int] = {start: None}

    while queue:
        current_node: int | None = None
        current_node, _ = queue.popleft()

        # If this node is the target, reconstruct and return the path
        if current_node in targets:
            path: list[int] = []
            while current_node is not None:
                path.insert(0, current_node)
                current_node = parent[current_node]
            return path

        visited.add(current_node)

        for neighbor in graph[current_node]:
            if neighbor not in visited:
                parent[neighbor] = current_node
                queue.append((neighbor, current_node))
                visited.add(neighbor)

    raise common.GammaLoopError(
        "Cannot find path between start vertex and target ones.")


def find_all_paths(edge_map: list[tuple[int, int]], start: int, dest: int, excluding: set[int] | None = None) -> list[list[tuple[int, bool]]]:
    if excluding is None:
        excluding = set()
    # find all paths from source to dest
    loop = start == dest
    start_check = 1 if loop else 0
    paths: list[list[tuple[int, bool]]] = [[(start, True)]]  # store direction
    if not loop:
        paths.append([(start, False)])
    res: list[list[tuple[int, bool]]] = []
    while True:
        newpaths: list[list[tuple[int, bool]]] = []
        for p in paths:  # pylint: disable=invalid-name
            if len(p) > 1 and p[-1][0] == dest:
                res.append(p)
                continue
            last_vertex = edge_map[p[-1][0]
                                   ][1] if p[-1][1] else edge_map[p[-1][0]][0]
            for i, x in enumerate(edge_map):  # pylint: disable=invalid-name
                if i not in excluding and all(i != pp[0] for pp in p[start_check:]):
                    if loop and i == start:
                        # if we need a loop, we need to enter from the right direction
                        if x[0] == last_vertex:
                            newpaths.append(p + [(i, True)])
                        continue

                    if x[0] == last_vertex and all(x[1] not in edge_map[pp[0]] for pp in p[start_check:-1]):
                        newpaths.append(p + [(i, True)])
                    if x[1] == last_vertex and all(x[0] not in edge_map[pp[0]] for pp in p[start_check:-1]):
                        newpaths.append(p + [(i, False)])
        paths = newpaths
        if len(paths) == 0:
            break
    return res


def generate_momentum_flow(edge_map: list[tuple[int, int]], loop_momenta: list[int], sink_edge: int, ext: list[int]) -> list[tuple[list[int], list[int]]]:
    """ Specify sing_edge to None for a vacuum graph"""

    flows: list[list[tuple[int, bool]]] = []
    for loop_momentum in loop_momenta:
        paths = find_all_paths(edge_map,
                               loop_momentum, loop_momentum, excluding={lm for lm in loop_momenta if lm != loop_momentum})
        assert len(paths) == 1
        flows.append(paths[0][:-1])

    # now route the external loop_momenta to the sink
    ext_flows: list[tuple[int, list[tuple[int, bool]]]] = []
    for i, e in enumerate(ext):  # pylint: disable=invalid-name
        if e == sink_edge:
            continue
        paths = find_all_paths(edge_map, e, sink_edge,
                               excluding=set(loop_momenta))
        assert len(paths) == 1
        ext_flows.append((i, paths[0]))

    # propagator momenta
    signatures: list[tuple[list[int], list[int]]] = [
        ([0 for _ in range(len(loop_momenta))], [0 for _ in range(len(ext))])
        for _ in range(len(edge_map))]
    for i, _ in enumerate(edge_map):
        if i in ext:
            if i == sink_edge:
                for j, y in ext_flows:  # pylint: disable=invalid-name
                    overall_sign = 1 if y[0][1] else -1
                    for yy in y:  # pylint: disable=invalid-name
                        if yy[0] == i:
                            prop_sign = 1 if yy[1] else -1
                            signatures[i][1][j] += prop_sign * overall_sign
                            break
            else:
                signatures[i][1][ext.index(i)] += 1
        elif i in loop_momenta:
            signatures[i][0][loop_momenta.index(i)] += 1
        else:
            for j, y in enumerate(flows):  # pylint: disable=invalid-name
                for yy in y:  # pylint: disable=invalid-name
                    if yy[0] == i:
                        signatures[i][0][j] += (1 if yy[1] else -1)
                        break
            for j, y in ext_flows:  # pylint: disable=invalid-name
                overall_sign = 1 if y[0][1] else -1
                for yy in y:  # pylint: disable=invalid-name
                    if yy[0] == i:
                        signatures[i][1][j] += overall_sign * \
                            (1 if yy[1] else -1)
                        break

    return signatures


def format_elapsed(elapsed_seconds: float) -> str:
    ms_in_a_day = 86400000
    time_in_ms = round(elapsed_seconds*1000)
    time_remainder = time_in_ms % ms_in_a_day
    n_days = int((time_in_ms-time_remainder)/ms_in_a_day)
    day_prefix = "" if n_days == 0 else f"{n_days} day{'s' if n_days > 1 else ''}, "  # nopep8
    if time_remainder < 60_000:
        hours_suffix = time.strftime('%H:%M:%S.{:03d}'.format(
            time_remainder % 1000), time.gmtime(time_remainder/1000.))
    else:
        hours_suffix = time.strftime(
            '%H:%M:%S', time.gmtime(time_remainder/1000.))
    return day_prefix+hours_suffix
