from __future__ import annotations
from argparse import Namespace
from typing import Literal
from gammaloop.base_objects.model import Particle
from gammaloop.misc.common import GammaLoopError
from gammaloop.base_objects.model import Model
import re


class MalformedProcessError(GammaLoopError):
    pass


class Process(object):
    def __init__(self,
                 initial_particles: list[Particle],
                 final_particles: list[Particle],
                 amplitude_orders: dict[str, int] | None = None,
                 cross_section_orders: dict[str, int] | None = None,
                 perturbative_orders: dict[str, int] | None = None,
                 amplitude_loop_count: tuple[int, int] | None = None,
                 cross_section_loop_count: tuple[int, int] | None = None,
                 particle_vetos: list[Particle] | None = None,
                 ):
        self.initial_states = initial_particles
        self.final_states = final_particles
        self.amplitude_orders = amplitude_orders
        self.cross_section_orders = cross_section_orders
        self.perturbative_orders = perturbative_orders
        self.amplitude_loop_count = amplitude_loop_count
        self.cross_section_loop_count = cross_section_loop_count
        self.particle_vetos = particle_vetos

    def __repr__(self) -> str:
        process_str_pieces: list[str] = []
        for p in self.initial_states:
            process_str_pieces.append(p.name.lower())
        process_str_pieces.append(">")
        for p in self.final_states:
            process_str_pieces.append(p.name.lower())
        if self.particle_vetos is not None:
            process_str_pieces.append("/")
            for p in self.particle_vetos:
                process_str_pieces.append(p.name.lower())
        if self.amplitude_orders is not None:
            for k, v in sorted(self.amplitude_orders.items(), key=lambda x: x[0]):
                process_str_pieces.append(f"{k}={v}")
        if self.perturbative_orders is not None or self.amplitude_loop_count is not None or self.cross_section_loop_count is not None:
            process_str_pieces.append("[")
            if self.amplitude_loop_count is not None:
                if self.amplitude_loop_count[0] == self.amplitude_loop_count[1]:
                    process_str_pieces.append(
                        f"{{{self.amplitude_loop_count[0]}}}")
                else:
                    process_str_pieces.append(
                        f"{{{','.join([str(i) for i in self.amplitude_loop_count])}}}")
            if self.cross_section_loop_count is not None:
                if self.cross_section_loop_count[0] == self.cross_section_loop_count[1]:
                    process_str_pieces.append(
                        f"{{{{{self.cross_section_loop_count[0]}}}}}")
                else:
                    process_str_pieces.append(
                        f"{{{{{','.join([str(i) for i in self.cross_section_loop_count])}}}}}")
            if self.perturbative_orders is not None:
                for k, v in sorted(self.perturbative_orders.items(), key=lambda x: x[0]):
                    if v > 1:
                        process_str_pieces.append(f"{k}={v}")
                    else:
                        process_str_pieces.append(f"{k}")
            process_str_pieces.append("]")

        if self.cross_section_orders is not None:
            for k, v in sorted(self.cross_section_orders.items(), key=lambda x: x[0]):
                process_str_pieces.append(f"{k}^2={v}")

        return ' '.join(process_str_pieces)

    @classmethod
    def process_perturbative_orders(cls,
                                    model_coupling_orders: list[str],
                                    perturbative_order_str: str,
                                    perturbative_orders: dict[str, int]
                                    ) -> tuple[tuple[int, int] | None, tuple[int, int] | None]:
        res_amplitude_loop_count = None
        res_cross_section_loop_count = None
        processed_input = perturbative_order_str.replace(' ', '').replace(
            '}}', 'XSEC_LOOP_COUNT_DEMLIMITER_CLOSE').replace(
            '{{', 'XSEC_LOOP_COUNT_DEMLIMITER_OPEN')
        for subs in [('=', ' = '), ('{', ' { '), ('}', ' } '), ('XSEC_LOOP_COUNT_DEMLIMITER_OPEN', ' {{ '), ('XSEC_LOOP_COUNT_DEMLIMITER_CLOSE', ' }} ')]:
            processed_input = processed_input.replace(subs[0], subs[1])
        loop_count_strs: list[str] = []
        po_parsing_stage = "orders"
        previous_order_before_equal: str | None = None
        previous_was_equal = False
        for t in processed_input.split(' '):
            if t == "":
                continue
            if po_parsing_stage == "loop_count":
                if t.endswith("}"):
                    loop_count_strs[-1] += t
                    po_parsing_stage = "orders"
                else:
                    loop_count_strs[-1] += t
            elif t.startswith("{"):
                if previous_order_before_equal is not None:
                    perturbative_orders[previous_order_before_equal] = 1
                    previous_order_before_equal = None
                loop_count_strs.append(t)
                if not t.endswith("}"):
                    po_parsing_stage = "loop_count"
            elif t == "}":
                loop_count_strs[-1] += t
            else:
                if t == "=":
                    if previous_order_before_equal is None:
                        raise MalformedProcessError(
                            "No perturbative order specified before '='.")
                    previous_was_equal = True
                    continue
                if previous_was_equal:
                    if previous_order_before_equal is None:
                        raise MalformedProcessError(
                            "No perturbative order specified before '='.")
                    try:
                        perturbative_orders[previous_order_before_equal] = int(
                            t)
                    except ValueError:
                        raise MalformedProcessError(
                            f"Invalid perturbative order value {t}.")
                    previous_was_equal = False
                    previous_order_before_equal = None
                    continue

                if t not in model_coupling_orders:
                    raise MalformedProcessError(
                        f"Unknown coupling order {t}.")

                if previous_order_before_equal is not None:
                    perturbative_orders[previous_order_before_equal] = 1

                previous_was_equal = False
                previous_order_before_equal = t

        if previous_order_before_equal is not None:
            perturbative_orders[previous_order_before_equal] = 1
            previous_order_before_equal = None

        for loop_count_specifier in loop_count_strs:
            spec = loop_count_specifier.replace(' ', '')
            loop_count_type_is_amplitude: bool = False
            loop_count_string: str = ""
            if spec.startswith("{{") and spec.endswith("}}"):
                try:
                    loop_count_type_is_amplitude, loop_count_string = False, spec[2:-2]
                except ValueError:
                    raise MalformedProcessError(
                        f"Invalid cross-section loop count specifier {spec}.")
            elif spec.startswith("{") and spec.endswith("}"):
                try:
                    loop_count_type_is_amplitude, loop_count_string = True, spec[1:-1]
                except ValueError:
                    raise MalformedProcessError(
                        f"Invalid amplitude loop count specifier {spec}.")
            else:
                raise MalformedProcessError(
                    f"Invalid loop count specifier {spec}.")
            loop_counts = [lc.strip()
                           for lc in loop_count_string.split(',')]
            if len(loop_counts) > 2:
                raise MalformedProcessError(
                    f"Invalid loop count specifier {spec}.")
            elif len(loop_counts) == 2:
                try:
                    loop_count_range = (
                        int(loop_counts[0]), int(loop_counts[1]))
                except ValueError:
                    raise MalformedProcessError(
                        f"A Invalid loop count specifier {spec}.")
            else:
                try:
                    loop_count_range = (
                        int(loop_counts[0]), int(loop_counts[0]))
                except ValueError:
                    raise MalformedProcessError(
                        f"B Invalid loop count specifier {spec}.")
            if loop_count_type_is_amplitude:
                res_amplitude_loop_count = loop_count_range
            else:
                res_cross_section_loop_count = loop_count_range

        return (res_amplitude_loop_count, res_cross_section_loop_count)

    @ classmethod
    def from_input_args(cls, model: Model, process_args: Namespace) -> Process:
        """ Creates a process object from the user input string representation. """

        initial_particles: list[Particle] = []
        final_particles: list[Particle] = []
        perturbative_orders: dict[str, int] = {}
        amplitude_loop_count: tuple[int, int] | None = None
        cross_section_loop_count: tuple[int, int] | None = None
        particle_vetos: list[Particle] = []
        perturbative_orders: dict[str, int] = {}
        amplitude_orders: dict[str, int] = {}
        cross_section_orders: dict[str, int] = {}

        model_coupling_orders = model.get_coupling_orders()

        def check_process_defined(check_final=True):
            if len(initial_particles) == 0:
                raise MalformedProcessError(
                    "No initial state specified.")
            if check_final and len(final_particles) == 0:
                raise MalformedProcessError(
                    "No final state specified.")

        parsing_stage: str = "initial_states"
        current_coupling_order: str | None = None
        is_coupling_order_for_xsec = False
        current_perturbative_order_str: str = ""
        process_string = ' '.join(process_args.process)
        for l in ['=', '[', ']', '>', '/']:
            process_string = process_string.replace(l, ' '+l+' ')

        tokens = [token.strip()
                  for token in process_string.split(' ') if token != '']
        for i_t, t in enumerate(tokens):
            if i_t+1 < len(tokens):
                next_token = tokens[i_t+1]
            else:
                next_token = None
            match t:
                case ">":
                    check_process_defined(check_final=False)
                    parsing_stage = "final_states"
                case "/":
                    check_process_defined()
                    parsing_stage = "vetos"
                case "[":
                    check_process_defined()
                    if t.endswith("]"):
                        amplitude_loop_count, cross_section_loop_count = Process.process_perturbative_orders(
                            model_coupling_orders, t[1:-1], perturbative_orders)
                    else:
                        parsing_stage = "perturbative_orders"
                case "]":
                    amplitude_loop_count, cross_section_loop_count = Process.process_perturbative_orders(
                        model_coupling_orders, current_perturbative_order_str, perturbative_orders)
                    parsing_stage = "orders"
                case _:
                    match (parsing_stage, next_token == "="):
                        case ("waiting_for_equal", _):
                            if t != "=":
                                raise MalformedProcessError(
                                    f"Expected '=' after coupling order {current_coupling_order}.")
                            parsing_stage = "waiting_on_order_value"
                        case ("waiting_on_order_value", _):
                            try:
                                order_value = int(t)
                            except ValueError:
                                raise MalformedProcessError(
                                    f"Invalid order value {t}.")
                            if current_coupling_order is None:
                                raise MalformedProcessError(
                                    "No coupling order specified before the equal sign.")
                            if is_coupling_order_for_xsec:
                                cross_section_orders[current_coupling_order] = order_value
                            else:
                                amplitude_orders[current_coupling_order] = order_value
                            parsing_stage = "orders"
                        case ("perturbative_orders", _):
                            if t.endswith("]"):
                                current_perturbative_order_str += " "+t[:-1]
                                amplitude_loop_count, cross_section_loop_count = Process.process_perturbative_orders(
                                    model_coupling_orders, current_perturbative_order_str, perturbative_orders)
                                parsing_stage = "orders"
                            else:
                                current_perturbative_order_str += t
                        case ("initial_states", _):
                            try:
                                initial_particles.append(model.get_particle(t))
                            except KeyError:
                                raise MalformedProcessError(
                                    f"Unknown particle {t} in initial state.")
                        case ("final_states", False):
                            try:
                                final_particles.append(model.get_particle(t))
                            except KeyError:
                                raise MalformedProcessError(
                                    f"Unknown particle {t} in final state.")
                        case ("vetos", False):
                            try:
                                veto_particle = model.get_particle(t)
                            except KeyError:
                                raise MalformedProcessError(
                                    f"Unknown particle {t} in vetoed particles specification.")
                            particle_vetos.append(veto_particle)
                        case _:
                            if t.endswith("^2"):
                                current_coupling_order = t[:-2]
                                is_coupling_order_for_xsec = True
                            else:
                                current_coupling_order = t
                                is_coupling_order_for_xsec = False
                            if current_coupling_order not in model_coupling_orders:
                                raise MalformedProcessError(
                                    f"Unknown coupling order {current_coupling_order}.")
                            parsing_stage = "waiting_for_equal"

        check_process_defined()

        return Process(
            initial_particles,
            final_particles,
            None if len(amplitude_orders) == 0 else amplitude_orders,
            None if len(cross_section_orders) == 0 else cross_section_orders,
            None if len(perturbative_orders) == 0 else perturbative_orders,
            amplitude_loop_count,
            cross_section_loop_count,
            None if len(particle_vetos) == 0 else particle_vetos,
        )
