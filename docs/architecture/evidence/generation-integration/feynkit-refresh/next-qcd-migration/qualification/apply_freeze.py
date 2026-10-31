"""Apply the exact reviewed patch, preserving the starting child and accepted base."""
import ast,hashlib,json,os,subprocess
from pathlib import Path
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent
P=W.parent/'generation-oracle-numerator-families/target/raw-next-qcd-migration-preparation-20261031-r1'
FILE='crates/feynkit-generator/tests/generation_oracle/main.rs'
env={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')}
commands=[]
def run(args):
 r=subprocess.run(args,cwd=W,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE);commands.append({'argv':args,'exit':r.returncode,'stdout':r.stdout.decode(),'stderr':r.stderr.decode()});(O/'source-commands.json').write_text(json.dumps(commands,indent=2)+'\n');assert r.returncode==0;return r.stdout.decode().strip()
def jj(*args):return run(['jj','--config','user.name="Lucien Huber"','--config','user.email="im@lcnbr.ch"',*args])
assert json.loads((O/'input-verification.json').read_text())['external_files']==3869
assert jj('log','--no-graph','-r','@','-T','commit_id')=='bd2642ac1a5cfd9d9bae35f5ada68360b62a2535'
jj('bookmark','create','generation-next-qcd-migration-start-20261031','-r','@')
class Context:
 MAIN=FILE
 @staticmethod
 def need(x,*detail):
  if not x:raise AssertionError(detail)
ns={'re':__import__('re'),'C':Context};tree=ast.parse((P/'reconstruct.py').read_bytes());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='apply_strict');exec(compile(ast.Module(body=[fn],type_ignores=[]),str(P/'reconstruct.py'),'exec'),ns)
old=(W/FILE).read_bytes();assert old==(P/'before.rs').read_bytes()
new=ns['apply_strict'](old.decode(),(P/'proposed-not-applied.patch').read_text()).encode();assert new==(P/'after.NOT-APPLIED.rs').read_bytes()
assert hashlib.sha1(b'blob %d\0'%len(new)+new).hexdigest()=='cc4505a6c4f6fc1687c075ed79e154a8e19f05a6'
(W/FILE).write_bytes(new)
jj('metaedit','--author','Lucien Huber <im@lcnbr.ch>','-m','migrate three approved QCD numerator expectations')
code=jj('log','--no-graph','-r','@','-T','commit_id')
assert jj('diff','--from','83276dc0','--to',code,'--name-only')==FILE
assert jj('log','--no-graph','-r',code,'-T','author.name() ++ " <" ++ author.email() ++ ">"')=='Lucien Huber <im@lcnbr.ch>'
jj('bookmark','create','generation-next-qcd-migration-code-20261031','-r',code)
jj('new','-m','qualify postapplication QCD migration evidence')
jj('metaedit','--author','Lucien Huber <im@lcnbr.ch>')
(O/'code.json').write_text(json.dumps({'code':code,'code_bookmark':'generation-next-qcd-migration-code-20261031','before':'83276dc047815d1dcc7f003aeb512c204b26758e','after_blob':'cc4505a6c4f6fc1687c075ed79e154a8e19f05a6','author':'Lucien Huber <im@lcnbr.ch>','patch':'Exact reviewed patch via strict no-offset/no-fuzz parser','start_bookmark':'generation-next-qcd-migration-start-20261031','global_integration_unchanged':True},indent=2)+'\n')
print('Frozen exact3-token code',code)
