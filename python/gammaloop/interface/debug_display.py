from logging import log
from typing import Dict, Any
import json
from gammaloop.misc.common import logger


def build_debug_dict(log_file: str) -> Dict[str, Any]:
    log = open(log_file, 'r')
    res = {}

    # certain messaages may appear multiple times due to the stability check,
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
    res_real = res['re']
    res_imag = res['im']

    if res_imag >= 0:
        format_result = "{} + {}i".format(res_real, res_imag)
    else:
        format_result = "{} - {}i".format(res_real, abs(res_imag))

    logger.info("final result: {}".format(format_result))


def display_default(debug_dict: Dict[str, Any]) -> None:
    display_havana_sample(debug_dict)
    display_momenta_samples(debug_dict)
    display_final_result(debug_dict)
