from logging import LogRecord, log
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
        res[msg] = json_line['data']
    return res


def build_general_debug_dict(log_file: str) -> Dict[str, Any]:
    general_log_path = pjoin(log_file, "general.jsonl")
    return parse_log_impl(general_log_path)


def build_eval_debug_dict(log_file: str, rotation: str, prec: str) -> Dict[str, Any]:
    file_name = "{}_{}.jsonl".format(rotation, prec)
    eval_log_path = pjoin(log_file, file_name)
    return parse_log_impl(eval_log_path)


def display_general(general_debug_dict: Dict[str, Any]) -> None:
    display_havana_sample(general_debug_dict)
    display_jacobian(general_debug_dict)
    display_final_result(general_debug_dict)


# I assume one sample from the integrator per point
def display_havana_sample(general_debug_dict: Dict[str, Any]) -> None:
    if 'havana_sample' not in general_debug_dict:
        logger.warn("no havana sample in debug info")
        return

    sample = general_debug_dict['havana_sample']
    logger.info("havana sample: ")
    logger.info(sample)


def display_jacobian(general_debug_dict: Dict[str, Any]) -> None:
    if 'jacobian' not in general_debug_dict:
        logger.warn("no jacobian in debug info")

    logger.info("jacobian: '%s%s%s'".format(general_debug_dict['jacobian']))


def display_final_result(general_debug_dict: Dict[str, Any]) -> None:
    if 'final_result' not in general_debug_dict:
        logger.warn("no final result in debug info")
        return

    res = general_debug_dict['final_result']
    format_result = format_complex(res['re'], res['im'])

    logger.info("final result: {}".format(format_result))


def display_momenta_samples(eval_debug_dict: Dict[str, Any]) -> None:
    if 'momenta_sample' not in eval_debug_dict:
        logger.warn("no momenta sample in debug info")
        return

    gammaloop_sample = eval_debug_dict['momenta_sample']
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

    orientation = eval_debug_dict['orientations']
    max_or = max(orientation)
    max_or_ind = orientation.index(max_or)
    min_or = min(orientation)
    min_or_ind = orientation.index(min_or)

    logger.info("largest orienatation: {}, value: {}".format(
        max_or_ind, max_or))
    logger.info("smallest orientation: {}, value: {}".format(
        min_or_ind, min_or))

    logger.info("")


def display_onshell_energies(eval_debug_dict: Dict[str, Any]) -> None:
    if 'onshell_energies' not in eval_debug_dict:
        logger.warn("onshell energies not logged")
        return

    onshell_energies = eval_debug_dict['onshell_energies']
    logger.info("onshell energies: ")
    logger.info(onshell_energies)
    logger.info('')


def display_rep3d(eval_debug_dict: Dict[str, Any]) -> None:
    if 'rep3d' not in eval_debug_dict:
        logger.warn("no rep3d in debug info")

    rep3d = eval_debug_dict['rep3d']
    logger.info(
        "3-dimensional representation: {}".format(format_complex(rep3d['re'], rep3d['im'])))
    logger.info('')


def display_counterterm(eval_debug_dict: Dict[str, Any]) -> None:
    if 'counter_terms' not in eval_debug_dict:
        logger.warn("no threshold counterterm in debug info")

    ct = eval_debug_dict['counter_terms']
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
        logger.info("\tcenter: {}".format(overlap_structure[1]))
        logger.info("\themispherical_radius: {}".format(
            debug_dict['hemispherical_radius'][overlap_index]))
        logger.info("")

        for index in range(len(overlap_structure[0])):
            esurface_subtraction_data = debug_dict["esurface_subtraction"][counter + index]
            logger.info("\t\tsubtraction: ")
            logger.info("\t\t{}".format(esurface_subtraction_data))


def format_3_momentum(momentum_dict: Dict[str, Any]) -> str:
    return "({}, {}, {})".format(momentum_dict['px'],
                                 momentum_dict['py'], momentum_dict['pz'])


def format_4_momentum(momentum_dict: Dict[str, Any]) -> str:
    energy = momentum_dict['temporal']['value']
    spatial_part = momentum_dict['spatial']
    return "({}, {}, {}, {})".format(energy, spatial_part['px'], spatial_part['py'], spatial_part['pz'])
