"""Capture exact pure bundle and controller plan before any private acquisition."""
import os
from pathlib import Path
import sys
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
from qualification import O,W,ASSIGN,PRIVATE,load,save
sys.path.insert(0,str(O))
import supervisor_guard as K
assert not PRIVATE.intersection(os.environ)
assert load(P/'runtime-qualification.json')['status']=='PASS'
b=K._check_bundle('faithfulness-historical--c78')
b['assignment']=K.verify_assignment(ASSIGN,b['supervisor_seal_sha256'])
for p,d in [(b['assignment']['assignment_path'],b['assignment']['assignment_sha256']),(b['assignment']['review_path'],b['assignment']['review_sha256'])]:
    b['protected_inputs'][p]=d
b['assignment_sha256']=b['assignment']['assignment_sha256']
save('validated-context.json',b)
assert not (O/'execution').exists()
print('ALL_TOOLS_PREPARED_AND_CONTEXT_PASS')
