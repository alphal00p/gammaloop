from logging import LogRecord, log
from typing import Dict, Any
import json
from gammaloop.misc.common import logger


def build_debug_dict(log_file: str) -> Dict[str, Any]:
    log = open(log_file, 'r')
    res = {}

    # certain messaages may appear multiple times due to the stability check or multichanneling,
    # this ensures we keep them all
    for line in log.readlines():
        json_line = json.loads(line)
        msg = json_line['msg']
        if msg in res:
            res[msg].append(json_line['data'])
        else:
            res[msg] = [json_line['data']]

    return res


# I assume one sample from the integrator per point
def display_havana_sample(debug_dict: Dict[str, Any]) -> None:
    if 'havana_sample' not in debug_dict:
        logger.warn("no havana sample in debug info")
        return

    sample = debug_dict['havana_sample'][0]
    logger.info("havana sample: ")
    logger.info(sample)


def display_momenta_samples(debug_dict: Dict[str, Any]) -> None:
    if 'momenta_sample' not in debug_dict:
        logger.warn("no momenta sample in debug info")
        return

    for gammaloop_sample in debug_dict['momenta_sample']:
        logger.info("momenta sampled: ")
        for index, loop_momentum in enumerate(gammaloop_sample['loop_moms']):
            logger.info("\t loop momentum {}: {}".format(index, loop_momentum))

        logger.info("")
        for index, external_momentum in enumerate(gammaloop_sample['external_moms']):
            logger.info("\t external momentum {}: {}".format(
                index, external_momentum))

        logger.info("")
        logger.info("\t jacobian: {}".format(gammaloop_sample['jacobian']))
        logger.info("")


def display_final_result(debug_dict: Dict[str, Any]) -> None:
    if 'final_result' not in debug_dict:
        logger.warn("no final result in debug info")
        return

    res = debug_dict['final_result'][0]
    format_result = format_complex(res['re'], res['im'])

    logger.info("final result: {}".format(format_result))


def display_largest_and_smallest_orientation(debug_dict: Dict[str, Any]) -> None:
    if 'orientations' not in debug_dict:
        logger.warn("no orientations in debug info")
        return

    for (index, orientation) in enumerate(debug_dict['orientations']):
        max_or = max(orientation)
        max_or_ind = orientation.index(max_or)
        min_or = min(orientation)
        min_or_ind = orientation.index(min_or)

        logger.info("largest orienatation: {}, value: {}".format(
            max_or_ind, max_or))
        logger.info("smallest orientation: {}, value: {}".format(
            min_or_ind, min_or))

        logger.info("")


def display_default(debug_dict: Dict[str, Any]) -> None:
    display_havana_sample(debug_dict)
    display_momenta_samples(debug_dict)
    display_onshell_energies(debug_dict)
    display_largest_and_smallest_orientation(debug_dict)
    display_rep3d(debug_dict)
    display_counterterm(debug_dict)
    display_final_result(debug_dict)


def display_onshell_energies(debug_dict: Dict[str, Any]) -> None:
    if 'onshell_energies' not in debug_dict:
        logger.warn("onshell energies not logged")
        return

    for onshell_energies in debug_dict['onshell_energies']:
        logger.info("onshell energies: ")
        logger.info(onshell_energies)
        logger.info('')


def display_rep3d(debug_dict: Dict[str, Any]) -> None:
    if 'rep3d' not in debug_dict:
        logger.warn("no rep3d in debug info")

    for rep3d in debug_dict['rep3d']:
        logger.info(
            "3-dimensional representation: {}".format(format_complex(rep3d['re'], rep3d['im'])))
        logger.info('')


def display_counterterm(debug_dict: Dict[str, Any]) -> None:
    if 'counter_terms' not in debug_dict:
        logger.warn("no threshold counterterm in debug info")

    for ct in debug_dict['counter_terms']:
        logger.info("threshold counterterm: {}".format(
            format_complex(ct['re'], ct['im'])))
        logger.info('')


def format_complex(re: float, im: float) -> str:
    if im >= 0:
        return "{} + {}i".format(re, im)
    else:
        return "{} - {}i".format(re, abs(im))


def display_subtraction_data(debug_dict: Dict[str, Any]) -> None:
    if 'overlap_structure' not in debug_dict:
        logger.info('no subtraction performed')
        return

    logger.warn(
        "this printout has only been checked when the stability test is disabled")

    counter = 0
    for overlap_structure in debug_dict['overlap_structure']:
        logger.info("overlap_structure: ")
        logger.info("\tesurfaces: {}".format(overlap_structure[0]))
        logger.info("\tcenter: {}".format(overlap_structure[1]))
        logger.info("")

        for index in range(len(overlap_structure[0])):
            esurface_subtraction_data = debug_dict["esurface_subtraction"][counter + index]
            logger.info("\t\tsubtraction: ")
            logger.info("\t\t{}".format(esurface_subtraction_data))
