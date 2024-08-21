from os.path import exists
from typing import Dict, Any
import json


def build_debug_dict(log_file: str) -> Dict[str, Any]:
    log = open(log_file, 'r')
    res = {}

    # certain messaages may appear multiple times due to the stability check,
    # this ensures we keep them all
    for line in log.readlines():
        json_line = json.loads(line)
        msg = json_line['msg']
        # not all messages have a value
        if 'value' in json_line:
            if msg in res:
                res[msg].append(json_line['value'])
            else:
                res[msg] = [json_line['value']]

    return res
