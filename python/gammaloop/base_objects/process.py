from __future__ import annotations
from argparse import Namespace
from typing import Literal
from gammaloop.base_objects.model import Particle
from gammaloop.misc.common import GammaLoopError
from gammaloop.base_objects.model import Model
from gammaloop.base_objects.graph import Graph
import yaml
import gammaloop._gammaloop as gl_rust
from gammaloop.misc.common import logger

COMPARISON_OPS = ["=", "==", "<=", ">="]


class MalformedProcessError(GammaLoopError):
    pass


class Process(object):
    def __init__(self,
                 initial_particles: list[Particle],
                 final_particles_lists: list[list[Particle]],
                 amplitude_orders: dict[str,
                                        tuple[int, int | None]] | None = None,
                 cross_section_orders: dict[str,
                                            tuple[int, int | None]] | None = None,
                 perturbative_orders: dict[str, int] | None = None,
                 amplitude_loop_count: tuple[int, int] | None = None,
                 cross_section_loop_count: tuple[int, int] | None = None,
                 particle_vetos: list[Particle] | None = None,
                 sewed_filter: gl_rust.SewedFilterOptions | None = None,
                 self_energy_filter: gl_rust.SelfEnergyFilterOptions | None = None,
                 tadpole_filter: gl_rust.TadpolesFilterOptions | None = None,
                 zero_snail_filter: gl_rust.SnailFilterOptions | None = None,
                 fermion_loop_count_range: tuple[int, int] | None = None,
                 factorized_loop_topologies_count_range: tuple[int,
                                                               int] | None = None,
                 max_n_bridges: int | None = None,
                 ):
        self.initial_states = initial_particles
        self.final_states_lists = final_particles_lists
        self.amplitude_orders = amplitude_orders
        self.cross_section_orders = cross_section_orders
        self.perturbative_orders: dict[str, int] | None = perturbative_orders
        self.amplitude_loop_count = amplitude_loop_count
        self.cross_section_loop_count = cross_section_loop_count
        self.particle_vetos: list[Particle] | None = particle_vetos
        self.amplitude_filters = gl_rust.FeynGenFilters(
            particle_veto=(None if self.particle_vetos is None else [
                p.get_pdg_code() for p in self.particle_vetos]),
            self_energy_filter=self_energy_filter,
            tadpoles_filter=tadpole_filter,
            zero_snails_filter=zero_snail_filter,
            max_number_of_bridges=max_n_bridges,
            coupling_orders=amplitude_orders,
            loop_count_range=amplitude_loop_count,
            fermion_loop_count_range=fermion_loop_count_range,
            factorized_loop_topologies_count_range=factorized_loop_topologies_count_range
        )

        # Adjust amplitude and cross-section orders given the perturbative orders
        if self.perturbative_orders is not None:
            if self.amplitude_orders is not None:
                amplitude_orders = self.amplitude_orders.copy()
            else:
                amplitude_orders = None
            if self.cross_section_orders is not None:
                cross_section_orders = self.cross_section_orders.copy()
            else:
                cross_section_orders = None

        self.cross_section_filters = gl_rust.FeynGenFilters(
            particle_veto=(None if self.particle_vetos is None else [
                p.get_pdg_code() for p in self.particle_vetos]),
            sewed_filter=sewed_filter,
            self_energy_filter=self_energy_filter,
            tadpoles_filter=tadpole_filter,
            zero_snails_filter=zero_snail_filter,
            max_number_of_bridges=max_n_bridges,
            perturbative_orders=self.perturbative_orders,
            coupling_orders=cross_section_orders,
            loop_count_range=cross_section_loop_count,
            fermion_loop_count_range=fermion_loop_count_range
        )

    def process_shell_name(self, short=True) -> str:
        res = []
        proc_name = (''.join([p.name.lower() for p in self.initial_states]) if len(
            self.initial_states) > 0 else "empty") + '_'
        proc_name += '_or_'.join((''.join([p.name.lower() for p in final_states]
                                          ) if len(self.initial_states) > 0 else "empty") for final_states in self.final_states_lists)
        res.append(proc_name)
        if not short:
            if self.particle_vetos is not None:
                res.append(
                    f"no__{'_'.join([p.name.lower() for p in self.particle_vetos])}")
            if self.amplitude_orders is not None:
                amp_orders = []
                for k, (v_min, v_max) in self.amplitude_orders.items():
                    if v_min == v_max:
                        amp_orders.append(f"{k}_eq_{v_min}")
                    else:
                        if v_min > 0:
                            amp_orders.append(
                                f"{k}_ge_{v_min}")
                        if v_max is not None:
                            amp_orders.append(
                                f"{k}_le_{v_max}")
                res.append('__'.join(amp_orders))
            if self.cross_section_orders is not None:
                xsec_orders = []
                for k, (v_min, v_max) in self.cross_section_orders.items():
                    if v_min == v_max:
                        xsec_orders.append(f"{k}sq_eq_{v_min}")
                    else:
                        if v_min > 0:
                            xsec_orders.append(
                                f"{k}sq_ge_{v_min}")
                        if v_max is not None:
                            xsec_orders.append(
                                f"{k}sq_le_{v_max}")
                res.append('__'.join(xsec_orders))
            if self.perturbative_orders is not None:
                res.append(
                    '__'.join([f"{k}loop_eq_{v}" for k, v in self.perturbative_orders.items()]))
        return '-'.join(s.replace('~', 'x').replace('+', 'p').replace('-', 'm') for s in res)

    def __repr__(self) -> str:
        process_str_pieces: list[str] = []
        if len(self.initial_states) == 0:
            process_str_pieces.append("{}")
        else:
            for p in self.initial_states:
                process_str_pieces.append(p.name.lower())
        process_str_pieces.append(">")
        if len(self.final_states_lists) == 1:
            if len(self.final_states_lists[0]) == 0:
                process_str_pieces.append("{}")
            else:
                for p in self.final_states_lists[0]:
                    process_str_pieces.append(p.name.lower())
        else:
            final_states_str_pieces = [
                ' '.join(p.name.lower() for p in final_states) for final_states in self.final_states_lists
            ]
            process_str_pieces.append(
                f"{{ {', '.join(final_states_str_pieces)} }}")
        if self.particle_vetos is not None:
            process_str_pieces.append("/")
            for p in self.particle_vetos:
                process_str_pieces.append(p.name.lower())
        if self.amplitude_orders is not None:
            amp_orders = []
            for k, (v_min, v_max) in sorted(self.amplitude_orders.items(), key=lambda x: x[0]):
                if v_min == v_max:
                    amp_orders.append(f"{k}=={v_min}")
                else:
                    if v_min > 0:
                        amp_orders.append(
                            f"{k}>={v_min}")
                    if v_max is not None:
                        amp_orders.append(
                            f"{k}<={v_max}")
            process_str_pieces.extend(amp_orders)
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
                for k, v in sorted(self.perturbative_orders.items(), key=lambda x: x[0]):  # type: ignore # nopep8
                    if v > 0:  # type: ignore # nopep8
                        if v > 1:  # type: ignore # nopep8
                            process_str_pieces.append(f"{k}={v}")
                        else:
                            process_str_pieces.append(f"{k}")
            process_str_pieces.append("]")

        if self.cross_section_orders is not None:
            xsec_orders = []
            for k, (v_min, v_max) in sorted(self.cross_section_orders.items(), key=lambda x: x[0]):
                if v_min == v_max:
                    xsec_orders.append(f"{k}^2=={v_min}")
                else:
                    if v_min > 0:
                        xsec_orders.append(
                            f"{k}^2>={v_min}")
                    if v_max is not None:
                        xsec_orders.append(
                            f"{k}^2<={v_max}")
            process_str_pieces.extend(xsec_orders)

        return ' '.join(process_str_pieces)

    def generate_diagrams(self, gl_worker: gl_rust.Worker, model: Model, generation_args: Namespace,
                          global_prefactor_color: str, global_prefactor_colorless: str) -> list[Graph]:

        if generation_args.number_of_samples_for_numerator_comparisons == 1:
            raise GammaLoopError(
                "The specified option 'number_of_samples_for_numerator_comparisons' must be > 1 or set to zero for disabling numberical numerator comparisons.")

        if generation_args.amplitude:
            loop_count_range = self.amplitude_loop_count
        else:
            loop_count_range = self.cross_section_loop_count

        # TODO: Improve automatic detection of ampltidue / cross section coupling orders and loop count
        if loop_count_range is None:
            if self.perturbative_orders is not None:
                loop_count = sum(self.perturbative_orders.values())
            else:
                loop_count = 0
            if not generation_args.amplitude:
                n_final_states = min(len(final_states)
                                     for final_states in self.final_states_lists)
                if n_final_states > 0:
                    loop_count += n_final_states - 1
            loop_count_range = (loop_count, loop_count)

        # Automatically adjust symmetrizations
        if generation_args.symmetrize_left_right_states is None:
            generation_args.symmetrize_left_right_states = False

        if generation_args.amplitude:
            if generation_args.symmetrize_initial_states is None:
                generation_args.symmetrize_initial_states = generation_args.symmetrize_left_right_states
            elif not generation_args.symmetrize_initial_states and generation_args.symmetrize_left_right_states:
                raise GammaLoopError(
                    "Symmetrization of left and right states for amplitudes requires also enabling initial-state symmetrization.")
            if generation_args.symmetrize_final_states is None:
                generation_args.symmetrize_final_states = generation_args.symmetrize_left_right_states
            elif not generation_args.symmetrize_final_states and generation_args.symmetrize_left_right_states:
                raise GammaLoopError(
                    "Symmetrization of left and right states for amplitudes requires also enabling final-state symmetrization.")
        else:
            if generation_args.symmetrize_initial_states is None:
                generation_args.symmetrize_initial_states = False
            if generation_args.symmetrize_final_states is None:
                generation_args.symmetrize_final_states = True
            elif not generation_args.symmetrize_final_states:
                raise GammaLoopError(
                    "Symmetrization of final states is mandatory for cross-section generation.")

        all_graphs: list[str] = gl_worker.generate_diagrams(
            gl_rust.FeynGenOptions(
                "amplitude" if generation_args.amplitude else "cross_section",
                [p.get_pdg_code() for p in self.initial_states],
                [[p.get_pdg_code() for p in final_states]
                 for final_states in self.final_states_lists],
                loop_count_range,
                tuple(generation_args.n_cut_blobs),
                tuple(generation_args.n_cut_spectators),
                generation_args.symmetrize_initial_states,
                generation_args.symmetrize_final_states,
                generation_args.symmetrize_left_right_states,
                generation_args.allow_symmetrization_of_external_fermions_in_amplitudes,
                generation_args.max_multiplicity_for_fast_cut_filter,
                amplitude_filters=self.amplitude_filters,
                cross_section_filters=self.cross_section_filters
            ),
            gl_rust.NumeratorAwareGroupingOption(
                generation_args.numerator_aware_isomorphism_grouping,
                generation_args.compare_canonized_numerator,
                generation_args.number_of_samples_for_numerator_comparisons,
                generation_args.consider_internal_masses_only_in_numerator_isomorphisms,
                generation_args.fully_numerical_substitution_when_comparing_numerators,
                generation_args.numerical_samples_seed
            ),
            generation_args.filter_self_loop,
            generation_args.graph_prefix,
            generation_args.select_graphs,
            generation_args.veto_graphs,
            generation_args.loop_momentum_bases,
            global_prefactor_color,
            global_prefactor_colorless,
            generation_args.num_threads
        )

        return [Graph.from_serializable_dict(model, yaml.safe_load(graph_str)) for graph_str in all_graphs]

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

    @classmethod
    def process_final_state_particles(cls, model: Model, final_states_str: str) -> list[list[Particle]]:
        final_states_str = final_states_str.strip()
        final_state_str_options = []
        if final_states_str.startswith("{") and final_states_str.endswith("}"):
            final_state_str_options = [s.strip()
                                       for s in final_states_str[1:-1].strip().split(",")]
        else:
            final_state_str_options.append(final_states_str.strip())

        if final_state_str_options == [""]:
            return [[]]
        final_state_options = []
        for final_states_str in final_state_str_options:
            final_states = []
            for p_str in final_states_str.split(" "):
                if p_str.strip() == "":
                    continue

                try:
                    final_states.append(model.get_particle(p_str.strip()))
                except KeyError:
                    raise MalformedProcessError(
                        f"Unknown particle '{p_str}' in final state. Particle in the model are:\n{sorted([p.name for p in model.particles])}")
            final_state_options.append(final_states)

        return final_state_options

    @classmethod
    def from_input_args(cls, model: Model, process_args: Namespace) -> Process:
        """ Creates a process object from the user input string representation. """

        initial_particles: list[Particle] = []
        final_particles_lists: list[list[Particle]] = []
        perturbative_orders: dict[str, int] = {}
        amplitude_loop_count: tuple[int, int] | None = None
        cross_section_loop_count: tuple[int, int] | None = None
        particle_vetos: list[Particle] = []
        perturbative_orders: dict[str, int] = {}
        amplitude_orders: dict[str, list[int | None]] = {}
        cross_section_orders: dict[str, list[int | None]] = {}

        model_coupling_orders = model.get_coupling_orders()

        accept_empty_initial_state = False
        current_final_states_tokens: list[str] = []

        def check_process_defined(check_final=True):
            if len(current_final_states_tokens) > 0:
                final_particles_lists.extend(cls.process_final_state_particles(
                    model, ' '.join(current_final_states_tokens)))
                current_final_states_tokens.clear()

            if len(initial_particles) == 0 and not accept_empty_initial_state:
                raise MalformedProcessError(
                    "No initial state specified.")
            if check_final:
                if len(final_particles_lists) == 0:
                    raise MalformedProcessError(
                        "No final state specified.")
                if len(final_particles_lists) > 1:
                    if process_args.amplitude:
                        raise MalformedProcessError(
                            "Multiple final state combinations are not supported for amplitude generation.")
                    elif any(len(final_states) == 0 for final_states in final_particles_lists):
                        raise MalformedProcessError(
                            "When specifying multiple set of final states for cross-section generation, each set must contain at least one particle.")

        parsing_stage: str = "initial_states"
        current_coupling_order: str | None = None
        current_comparison_operator: str | None = None
        is_coupling_order_for_xsec = False
        current_perturbative_order_str: str = ""
        process_string = ' '.join(process_args.process)

        # Handle the comparison operators specially to avoid splitting them with the '=' and '>' below
        special_reps = [('==', 'eq'), ('>=', 'geq'), ('<=', 'leq')]
        for src, trgt in special_reps:
            process_string = process_string.replace(
                src, f"____GLMARKER____{trgt}")
        for l in ['=', '[', ']', '>', '/', '|']:
            process_string = process_string.replace(l, f" {l} ")
        for src, trgt in special_reps:
            process_string = process_string.replace(
                f"____GLMARKER____{trgt}", f" {src} ")

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
                    particle_vetos = []
                    parsing_stage = "vetos"
                case "|":
                    check_process_defined()
                    particle_vetos.extend(
                        [p for p in model.particles if p.pdg_code >= 0])
                    parsing_stage = "vetos_complement"
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
                    match (parsing_stage, next_token in COMPARISON_OPS):
                        case ("waiting_for_equal", _):
                            if t not in COMPARISON_OPS:
                                raise MalformedProcessError(
                                    f"Expected comparison operator ({','.join(COMPARISON_OPS)}) after coupling order {current_coupling_order}, not '{t}'.")
                            current_comparison_operator = t
                            parsing_stage = "waiting_on_order_value"
                        case ("waiting_on_order_value", _):
                            try:
                                order_value = int(t)
                            except ValueError:
                                raise MalformedProcessError(
                                    f"Invalid order value {t}.")
                            if current_coupling_order is None:
                                raise MalformedProcessError(
                                    f"No coupling order specified before the comparison operator ({current_comparison_operator if current_comparison_operator is not None else ','.join(COMPARISON_OPS)})")
                            if current_comparison_operator is None:
                                raise MalformedProcessError(
                                    f"No comparison operator specified before the coupling order value.")
                            if is_coupling_order_for_xsec:
                                xsec_order = cross_section_orders.get(
                                    current_coupling_order, [0, None])
                                match current_comparison_operator:
                                    case "==" | "=":
                                        xsec_order: list[int | None] = [
                                            order_value, order_value]
                                    case "<=":
                                        xsec_order[1] = order_value
                                    case ">=":
                                        xsec_order[0] = order_value
                                    case _:
                                        raise MalformedProcessError(
                                            f"Invalid comparison operator {current_comparison_operator} for cross-section orders.")
                                cross_section_orders[current_coupling_order] = xsec_order
                            else:
                                amp_order = amplitude_orders.get(
                                    current_coupling_order, [0, None])
                                match current_comparison_operator:
                                    case "==" | "=":
                                        amp_order: list[int | None] = [
                                            order_value, order_value]
                                    case "<=":
                                        amp_order[1] = order_value
                                    case ">=":
                                        amp_order[0] = order_value
                                    case _:
                                        raise MalformedProcessError(
                                            f"Invalid comparison operator {current_comparison_operator} for amplitude orders.")
                                amplitude_orders[current_coupling_order] = amp_order
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
                            if t == "{}":
                                accept_empty_initial_state = True
                            else:
                                try:
                                    initial_particles.append(
                                        model.get_particle(t))
                                except KeyError:
                                    raise MalformedProcessError(
                                        f"Unknown particle '{t}' in initial state. Particle in the model are:\n{sorted([p.name for p in model.particles])}")
                        case ("final_states", False):
                            current_final_states_tokens.append(t)
                            # if t == "{}":
                            #     accept_empty_final_state = True
                            # else:
                            #     try:
                            #         final_particles.append(
                            #             model.get_particle(t))
                            #     except KeyError:
                            #         raise MalformedProcessError(
                            #             f"Unknown particle '{t}' in final state. Particle in the model are:\n{sorted([p.name for p in model.particles])}")
                        case ("vetos", False) | ("vetos_complement", False):
                            try:
                                veto_particle = model.get_particle(t)
                                if veto_particle.pdg_code < 0:
                                    logger.warning(
                                        "Particle veto %s in process definition is not a particle but an antiparticle. Automatically converting to particle.", veto_particle.name)
                                    veto_particle = veto_particle.get_anti_particle(
                                        model)
                            except KeyError:
                                raise MalformedProcessError(
                                    f"Unknown particle '{t}' in vetoed particles specification. Particle in the model are:\n{sorted([p.name for p in model.particles])}")
                            if parsing_stage == "vetos":
                                if veto_particle not in particle_vetos:
                                    particle_vetos.append(veto_particle)
                            else:
                                if veto_particle in particle_vetos:
                                    particle_vetos.remove(veto_particle)
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

        # Adjust filter default values
        # Disable these filter for vacuum generation
        is_vaccuum_topology = False
        if len(initial_particles) == 0:
            if not process_args.amplitude:
                is_vaccuum_topology = True
            else:
                if len(final_particles_lists[0]) == 0:
                    is_vaccuum_topology = True

        if is_vaccuum_topology:
            # For vaccuum-like graphs, the user will typically expect the number of bridges to be forced to zero
            if process_args.max_n_bridges is None:
                process_args.max_n_bridges = 0
            if process_args.number_of_factorized_loop_subtopologies is None:
                process_args.number_of_factorized_loop_subtopologies = 1

        if process_args.filter_tadpoles is None:
            process_args.filter_tadpoles = not is_vaccuum_topology
        if process_args.filter_snails is None:
            process_args.filter_snails = not is_vaccuum_topology
        if process_args.filter_selfenergies is None:
            process_args.filter_selfenergies = not is_vaccuum_topology

        if process_args.number_of_factorized_loop_subtopologies is not None and process_args.number_of_factorized_loop_subtopologies < 0:
            process_args.number_of_factorized_loop_subtopologies = None
        if process_args.max_n_bridges is not None and process_args.max_n_bridges < 0:
            process_args.max_n_bridges = None

        if len(amplitude_orders) == 0:
            amplitude_orders_input = None
        else:
            amplitude_orders_input = {k: (0 if v[0] is None else v[0], v[1])
                                      for k, v in amplitude_orders.items()}
        if len(cross_section_orders) == 0:
            cross_section_orders_input = None
        else:
            cross_section_orders_input = {k: (0 if v[0] is None else v[0], v[1])
                                          for k, v in cross_section_orders.items()}
        return Process(
            initial_particles,
            final_particles_lists,
            amplitude_orders_input,
            cross_section_orders_input,
            None if len(perturbative_orders) == 0 else perturbative_orders,
            amplitude_loop_count,
            cross_section_loop_count,
            None if len(particle_vetos) == 0 else particle_vetos,
            sewed_filter=None if not process_args.filter_cross_section_tadpoles else gl_rust.SewedFilterOptions(filter_tadpoles=process_args.filter_cross_section_tadpoles),
            self_energy_filter=None if not process_args.filter_selfenergies else gl_rust.SelfEnergyFilterOptions(
                veto_self_energy_of_massive_lines=process_args.veto_self_energy_of_massive_lines,
                veto_self_energy_of_massless_lines=process_args.veto_self_energy_of_massless_lines,
                veto_only_scaleless_self_energy=process_args.veto_only_scaleless_self_energy
            ),
            tadpole_filter=None if not process_args.filter_tadpoles else gl_rust.TadpolesFilterOptions(
                veto_tadpoles_attached_to_massive_lines=process_args.veto_tadpoles_attached_to_massive_lines,
                veto_tadpoles_attached_to_massless_lines=process_args.veto_tadpoles_attached_to_massless_lines,
                veto_only_scaleless_tadpoles=process_args.veto_only_scaleless_tadpoles,
            ),

            zero_snail_filter=None if not process_args.filter_snails else gl_rust.SnailFilterOptions(
                veto_snails_attached_to_massive_lines=process_args.veto_snails_attached_to_massive_lines,
                veto_snails_attached_to_massless_lines=process_args.veto_snails_attached_to_massless_lines,
                veto_only_scaleless_snails=process_args.veto_only_scaleless_snails
            ),
            max_n_bridges=process_args.max_n_bridges,
            fermion_loop_count_range=(None if process_args.number_of_fermion_loops is None else (
                process_args.number_of_fermion_loops, process_args.number_of_fermion_loops)),
            factorized_loop_topologies_count_range=(None if process_args.number_of_factorized_loop_subtopologies is None else (
                process_args.number_of_factorized_loop_subtopologies, process_args.number_of_factorized_loop_subtopologies))
        )
