from os.path import exists
from typing import Dict, Any
import json

from symbolica import Sample


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
        print("no havana sample in debug info")
        return

    nested_sample = (debug_dict['havana_sample'][0])
    continuous_found = False

    while continuous_found is False:
        if 'Continuous' in nested_sample:
            continuous_found = True
            continuous_sample = nested_sample["Continuous"]
            print(continuous_sample)


def display_default(debug_dict: Dict[str, Any]) -> None:
    display_havana_sample(debug_dict)
