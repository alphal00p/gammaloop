"""Prepare unchanged qualified capsule orchestration for sixteen exact identities."""
import ast
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];OLD=W/'.local/feynkit-refresh/grouped-identical-migration-20261031'
s=(OLD/'runtime_common.py').read_text();start=s.index('PRIOR12 =');end=s.index("FILTER =",start)
s=s[:start]+"PRIOR13 = json.loads((O / 'source-invariance.json').read_text())['prior13']\nNEW3 = ['xs_qcd_gguu_sym_L1', 'xs_qcd_gguu_L1', 'xs_qcd_gggg_L1']\nCASES = PRIOR13 + NEW3\n"+s[end:]
# Exactly one identity per nextest invocation; existing16-worker config and groups unchanged.
(O/'runtime_common.py').write_text(s)
s=(OLD/'runtime_prepare.py').read_text().replace('len(found) == 13','len(found) == 16').replace('and 13 exact compiled identities.','and 16 exact compiled identities.')
s=s.replace("write(O / 'selected-inventory.json',", "write(O / 'discovery-scope.json', {'authorization':'Bounded inventory discovery only','test_executions':0,'pass_credit':0,'selected_names':CASES,'command_record':str(e/'command.json')})\nwrite(O / 'selected-inventory.json',")
(O/'runtime_prepare.py').write_text(s)
for p in O.glob('*.py'):ast.parse(p.read_bytes())
print('Prepared capsule sources; no execution/license access')
