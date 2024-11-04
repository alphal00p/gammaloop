from logging import LogRecord, log
from re import A
from typing import Dict, Any
import json
from gammaloop.misc.common import logger, pjoin
from pathlib import Path


def parse_log_impl(file: str) -> Dict[str, Any]:
    log = open(file, 'r')
    res = {}

    # messages that appear twice are overwritten
    # The splitting into multiple files should ensure this never happens.
    for line in log.readlines():
        json_line = json.loads(line)
        msg = json_line['msg']
        if msg in res:
            res[msg].append(json_line['data'])
        else:
            res[msg] = [json_line['data']]
    return res


def build_general_debug_dict(log_file: str) -> Dict[str, Any]:
    general_log_path = pjoin(log_file, "general.jsonl")
    return parse_log_impl(general_log_path)


def build_eval_debug_dict(log_file: str, file_name: str) -> Dict[str, Any]:
    eval_log_path = pjoin(log_file, file_name)
    return parse_log_impl(eval_log_path)


def display_general(general_debug_dict: Dict[str, Any]) -> None:
    display_havana_sample(general_debug_dict)
    display_jacobian(general_debug_dict)
    display_final_result(general_debug_dict)
    logger.info("")


def display_havana_sample(general_debug_dict: Dict[str, Any]) -> None:
    if 'havana_sample' not in general_debug_dict:
        logger.warn("no havana sample in debug info")
        return

    # assume single sample
    sample = general_debug_dict['havana_sample'][0]
    logger.info("havana sample: ")
    logger.info(sample)


def display_jacobian(general_debug_dict: Dict[str, Any]) -> None:
    if 'jacobian' not in general_debug_dict:
        logger.warn("no jacobian in debug info")

    logger.info("jacobian: %s", general_debug_dict['jacobian'][0])


def display_final_result(general_debug_dict: Dict[str, Any]) -> None:
    if 'final_result' not in general_debug_dict:
        logger.warn("no final result in debug info")
        return

    res = general_debug_dict['final_result'][0]
    format_result = format_complex(res['re'], res['im'])

    logger.info("final result: {}".format(format_result))


def display_momenta_samples(eval_debug_dict: Dict[str, Any]) -> None:
    if 'momenta_sample' not in eval_debug_dict:
        logger.warn("no momenta sample in debug info")
        return

    gammaloop_sample = eval_debug_dict['momenta_sample'][0]
    logger.info("momenta sampled: ")
    for index, loop_momentum in enumerate(gammaloop_sample['loop_moms']):
        logger.info("\t loop momentum {}: {}".format(
            index, format_3_momentum(loop_momentum)))

    logger.info("")
    for index, external_momentum in enumerate(gammaloop_sample['external_moms']):
        logger.info("\t external momentum {}: {}".format(
            index, format_4_momentum(external_momentum)))

    logger.info("")
    logger.info("\t jacobian: {}".format(gammaloop_sample['jacobian']))
    logger.info("")


def display_largest_and_smallest_orientation(eval_debug_dict: Dict[str, Any]) -> None:
    if 'orientations' not in eval_debug_dict:
        logger.warn("no orientations in debug info")
        return

    orientation = eval_debug_dict['orientations'][0]
    max_or = max(orientation)
    max_or_ind = orientation.index(max_or)
    min_or = min(orientation)
    min_or_ind = orientation.index(min_or)

    logger.info("largest orienatation: {}, value: {}".format(
        max_or_ind, max_or))
    logger.info("smallest orientation: {}, value: {}".format(
        min_or_ind, min_or))

    logger.info("number of orientations: %s", len(orientation))

    logger.info("")


def display_onshell_energies(eval_debug_dict: Dict[str, Any]) -> None:
    if 'onshell_energies' not in eval_debug_dict:
        logger.warn("onshell energies not logged")
        return

    onshell_energies = eval_debug_dict['onshell_energies'][0]
    logger.info("onshell energies: ")
    logger.info(onshell_energies)
    logger.info('')


def display_rep3d(eval_debug_dict: Dict[str, Any]) -> None:
    if 'rep3d' not in eval_debug_dict:
        logger.warn("no rep3d in debug info")

    rep3d = eval_debug_dict['rep3d'][0]
    logger.info(
        "3-dimensional representation: {}".format(format_complex(rep3d['re'], rep3d['im'])))
    logger.info('')


def display_counterterm(eval_debug_dict: Dict[str, Any]) -> None:
    if 'counter_terms' not in eval_debug_dict:
        logger.warn("no threshold counterterm in debug info")

    ct = eval_debug_dict['counter_terms'][0]
    logger.info("threshold counterterm: {}".format(
        format_complex(ct['re'], ct['im'])))
    logger.info('')


def display_eval_default(eval_debug_dict: Dict[str, Any]) -> None:
    display_momenta_samples(eval_debug_dict)
    display_onshell_energies(eval_debug_dict)
    display_largest_and_smallest_orientation(eval_debug_dict)
    display_rep3d(eval_debug_dict)
    display_counterterm(eval_debug_dict)


def format_complex(re: float, im: float) -> str:
    if im >= 0:
        return "{} + {}i".format(re, im)
    else:
        return "{} - {}i".format(re, abs(im))


def display_subtraction_data(debug_dict: Dict[str, Any]) -> None:
    if 'overlap_structure' not in debug_dict:
        logger.info('no subtraction performed')
        return

    counter = 0
    for overlap_index, overlap_structure in enumerate(debug_dict['overlap_structure']):
        logger.info("overlap_structure: ")
        logger.info("\tesurfaces: {}".format(overlap_structure[0]))
        logger.info("\tcenter:")
        for momentum in overlap_structure[1]:
            logger.info("\t\t{}".format(format_3_momentum(momentum)))

        logger.info("\themispherical_radius: {}".format(
            debug_dict['hemispherical_radius'][overlap_index]))
        logger.info("")

        for index in range(len(overlap_structure[0])):
            esurface_subtraction_data = debug_dict["esurface_subtraction"][counter + index]
            logger.info("subtracting esurface %s",
                        esurface_subtraction_data['esurface_id'])
            logger.info("edges: %s", esurface_subtraction_data['edges'])

            plus_solution, minus_solution = esurface_subtraction_data[
                'plus_solution'], esurface_subtraction_data['minus_solution']

            logger.info("r*+  = %s, r*-  = %s",
                        plus_solution['solution'], minus_solution['solution'])

            logger.info("∇η*+ = %s, ∇η*- = %s",
                        plus_solution['derivative_at_solution'], minus_solution['derivative_at_solution'])

            logger.info("j*+/j = %s, j*-/j = %s",
                        esurface_subtraction_data['jacobian_ratio_plus'], esurface_subtraction_data['jacobian_ratio_minus'])

            logger.info("uv_d+ = %s, uv_d = %s",
                        esurface_subtraction_data['uv_damper_plus'], esurface_subtraction_data['uv_damper_minus'])

            logger.info("ir_d+ = %s, ir_d = %s",
                        esurface_subtraction_data['singularity_dampener_plus'], esurface_subtraction_data['singularity_dampener_minus'])

            logger.info("r*+ res = %s, r*- res = %s",
                        esurface_subtraction_data['r_plus_eval'], esurface_subtraction_data['r_minus_eval'])

            ct_plus_dict = esurface_subtraction_data['ct_plus']
            ct_minus_dict = esurface_subtraction_data['ct_minus']
            ict_plus_dict = esurface_subtraction_data['integrated_ct_plus']
            ict_minus_dict = esurface_subtraction_data['integrated_ct_minus']

            ct_plus = format_complex(ct_plus_dict['re'], ct_plus_dict['im'])
            ct_minus = format_complex(ct_minus_dict['re'], ct_minus_dict['im'])
            ict_plus = format_complex(ict_plus_dict['re'], ict_plus_dict['im'])
            ict_minus = format_complex(
                ict_minus_dict['re'], ict_minus_dict['im'])

            logger.info("ct+ = %s, ct- = %s",
                        ct_plus, ct_minus)

            logger.info("ict+ = %s, ict- = %s",
                        ict_plus, ict_minus)

            logger.info("h_plus = %s, h_minus = %s",
                        esurface_subtraction_data['h_plus'], esurface_subtraction_data['h_minus'])

            r_plus_energy_cache = esurface_subtraction_data['r_plus_energy_cache']
            r_plus_energy_product = esurface_subtraction_data['r_plus_energy_product']
            r_minus_energy_cache = esurface_subtraction_data['r_minus_energy_cache']
            r_minus_energy_product = esurface_subtraction_data['r_minus_energy_product']

            logger.info("Π2E*+  = %s, Π2E*-= %s",
                        r_plus_energy_product, r_minus_energy_product)

            logger.info("r*+ energy cache")
            logger.info(r_plus_energy_cache)

            logger.info("r*- energy cache")
            logger.info(r_minus_energy_cache)

            logger.info("r*+ esurface_cache:")
            logger.info(esurface_subtraction_data['r_plus_esurface_cache'])

            logger.info("r*- esurface_cache:")
            logger.info(esurface_subtraction_data['r_minus_esurface_cache'])

            logger.info("")

        counter += len(overlap_structure[0])


def format_3_momentum(momentum_dict: Dict[str, Any]) -> str:
    return "({}, {}, {})".format(momentum_dict['px'],
                                 momentum_dict['py'], momentum_dict['pz'])


def format_4_momentum(momentum_dict: Dict[str, Any]) -> str:
    energy = momentum_dict['temporal']['value']
    spatial_part = momentum_dict['spatial']
    return "({}, {}, {}, {})".format(energy, spatial_part['px'], spatial_part['py'], spatial_part['pz'])
