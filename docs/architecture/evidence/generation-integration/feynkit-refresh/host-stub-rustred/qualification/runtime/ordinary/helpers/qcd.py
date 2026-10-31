"""Observe the unchanged script's actual failing API boundary, then propagate it."""
import hashlib
import json
from pathlib import Path
import runpy
import sys

import symbolica.core
from symbolica.community import feynkit_native
from symbolica.community import hepkit

script = '/build/source/crates/feynkit-py/tests/installed_feyncalc_qcd_ghost_vertex.py'
receipt = {
    'script': script,
    'sha256': hashlib.sha256(Path(script).read_bytes()).hexdigest(),
    'native_file': symbolica.core.__file__,
    'native_sha256': hashlib.sha256(Path(symbolica.core.__file__).read_bytes()).hexdigest(),
    'hepkit_file': hepkit.__file__,
    'hepkit_exports': sorted(n for n in dir(hepkit) if not n.startswith('_')),
    'native_exports': sorted(n for n in dir(feynkit_native) if not n.startswith('_')),
    'IBPFamily_in_hepkit': hasattr(hepkit, 'IBPFamily'),
    'IBPFamily_in_native': hasattr(feynkit_native, 'IBPFamily'),
    'IntegralFamily_members': sorted(n for n in dir(hepkit.IntegralFamily) if not n.startswith('_')),
    'script_assertions_modified': False,
}
try:
    scope = runpy.run_path(script, run_name='__main__')
except BaseException as error:
    tb = sys.exc_info()[2]
    while tb and tb.tb_frame.f_code.co_filename != script:
        tb = tb.tb_next
    if tb:
        scope = tb.tb_frame.f_globals
        receipt['failure_line'] = tb.tb_lineno
        receipt['targets'] = scope.get('targets')
        receipt['terms_by_diagram'] = [
            {'key': [pdg, kind, list(topology)],
             'terms': [{'powers': powers, 'coefficient': str(coefficient)} for powers, coefficient in terms]}
            for (pdg, kind, topology), terms in scope.get('terms_by_diagram', {}).items()
        ]
        receipt['trees'] = {str(pdg): str(tree) for pdg, tree in scope.get('trees', {}).items()}
    receipt['error_type'] = type(error).__name__
    receipt['error'] = str(error)
    receipt['complete_script_passed'] = False
    Path('/work/qcd-boundary.json').write_text(json.dumps(receipt, indent=2) + '\n')
    raise
else:
    receipt.update({'targets':scope.get('targets'),'residuals':scope['solution'].residuals,'trees_count':len(scope.get('trees',{})),'loop_groups':len(scope.get('terms_by_diagram',{}))})
    receipt['complete_script_passed'] = True
    Path('/work/qcd-boundary.json').write_text(json.dumps(receipt, indent=2) + '\n')
