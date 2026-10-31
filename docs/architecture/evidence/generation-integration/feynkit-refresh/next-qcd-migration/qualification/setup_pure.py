"""Reuse accepted producer and qualification scripts with task-specific identities."""
import ast,json
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];OLD=W/'.local/feynkit-refresh/grouped-identical-migration-20261031'
code=json.loads((O/'code.json').read_text())['code']
s=(OLD/'pure.py').read_text().replace("REV = '6f41807b064bcf2f579883134daef65c8ac279d3'",f"REV = '{code}'")
s=s.replace("cwd='/tmp'","cwd=str(W)").replace("'cwd':'/tmp'","'cwd':str(W)")
s=s.replace("['nix','--offline','--option','allow-import-from-derivation','false',*args]","['nix',*(['--offline'] if args[0]!='build' else []),'--option','allow-import-from-derivation','false',*args]")
s=s.replace("nix('build',['build','--no-link','--json',drv+'^out'])","assert not (O/'build.command.json').exists(), 'one new producer build only'\n    nix('build',['build','--max-jobs','1','--cores','2','--no-link','--json',drv+'^out'])")
(O/'pure.py').write_text(s)
for name in ['qualify.py','dependency_invariance.py','elf_helpers.py']:(O/name).write_bytes((OLD/name).read_bytes())
s=(O/'dependency_invariance.py').read_text().replace('raw-three-migration-20261031','grouped-identical-migration-20261031').replace('accepted_raw_three','accepted_grouped');(O/'dependency_invariance.py').write_text(s)
# Original452-unit dependency index is reused only after exact unchanged-layer checks.
# Add a direct producer-source comparison with the immediate accepted grouped producer.
s=(O/'qualify.py').read_text()
s=s.replace("write('selected-source-delta.json',", "grouped=json.loads((W/'.local/feynkit-refresh/grouped-identical-migration-20261031/archive-qualification.json').read_text())\ngsrc=Path(grouped['selected_source']);assert set(a)=={str(p.relative_to(gsrc)) for p in gsrc.rglob('*') if p.is_file() or p.is_symlink()}\ngdelta=[p for p in a if (Path(src)/p).read_bytes()!=(gsrc/p).read_bytes()];assert gdelta==delta\nassert (Path(src)/delta[0]).read_bytes()==(W.parent/'generation-oracle-numerator-families/target/raw-next-qcd-migration-preparation-20261031-r1/after.NOT-APPLIED.rs').read_bytes()\nwrite('immediate-producer-source-delta.json',{'accepted_grouped_code':grouped['code'],'old':str(gsrc),'new':src,'paths':len(a),'changed_paths':gdelta,'exact_approved_after_bytes':True})\nwrite('selected-source-delta.json',")
(O/'qualify.py').write_text(s)
for p in O.glob('*.py'):ast.parse(p.read_bytes())
print('Prepared reused pure workflow for',code)
