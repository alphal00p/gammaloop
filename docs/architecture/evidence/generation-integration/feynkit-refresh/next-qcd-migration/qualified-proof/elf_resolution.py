"""Corrected pure ELF breadth-first loaded-library resolution, separate execs."""
from pathlib import Path
import os
import sys
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
from qualification import load,save,elf,PYTHON,PRIVATE,O,W
sys.path.insert(0,str(O))
from lifecycle import bounded
assert not PRIVATE.intersection(os.environ)
old=load(P/'elf-qualification.json')
corrected=[]
for config in old['configurations']:
    mounts=config['readonly_mounts']
    def exposed(path):
        resolved=Path(path).resolve()
        return any(resolved==Path(src).resolve() or (Path(src).is_dir() and resolved.is_relative_to(Path(src).resolve())) for src,dest in mounts if src==dest)
    processes=[]
    objects={}
    for main in [config['program']['path'],PYTHON]:
        loaded={};queue=[];errors=[];edges=[]
        root=elf(main);objects[main]=root
        interp=root['interpreter']
        assert exposed(interp) and Path(interp).is_file()
        loaded[Path(interp).name]=interp
        objects[interp]=elf(interp)
        queue.extend([(main,()),(interp,())])
        processed=set()
        while queue:
            path,inherited=queue.pop(0)
            if path in processed:continue
            processed.add(path)
            info=elf(path);objects[path]=info
            origin=str(Path(path).parent)
            expand=lambda s:s.replace('${ORIGIN}',origin).replace('$ORIGIN',origin)
            own=info['runpath'] or info['rpath']
            dirs=[expand(d) for d in own if d]+list(inherited)
            child_inherited=tuple([expand(d) for d in info['rpath'] if d]+list(inherited))
            for name in info['needed']:
                existing=loaded.get(name)
                found=existing or next((str(Path(d)/name) for d in dirs if (Path(d)/name).is_file() and exposed(Path(d)/name)),None)
                if found:
                    edges.append({'object':path,'needed':name,'resolved':found,'already_loaded':existing is not None})
                    loaded[name]=found
                    queue.append((found,child_inherited))
                else:errors.append({'object':path,'needed':name,'search':dirs})
        processes.append({'main':main,'interpreter':interp,'loaded':loaded,'dependency_edges':edges,'errors':errors})
    corrected.append({**config,'objects':objects,'missing':[e for proc in processes for e in proc['errors']],'separate_exec_processes':processes})
value={'status':'PASS' if not any(c['missing'] for c in corrected) else 'FAIL','configurations':corrected,
       'prior_tool_error':'qualification.py scanned recursive children before registering the complete directly loaded set and shared the two distinct exec namespaces. The first scan is retained; no native launch or missing mount was observed.',
       'mounts_cards_supervisor_unchanged':True}
save('elf-qualification-corrected.json',value)
assert value['status']=='PASS', 'actual unresolved ELF dependency; stop without native/private execution'
with (P/'default-pure.stdout').open('xb') as out,(P/'default-pure.stderr').open('xb') as err:
    result=bounded([PYTHON,'-I','-S','-B',str(O/'supervisor.py'),'--id','faithfulness-historical--c78'],dict(os.environ),W,out,err,seconds=240)
save('default-pure-lifecycle.json',result)
assert result['exit']==0 and result['process_reaped'] and not result['timed_out']
assert not (O/'execution').exists()
print('CORRECTED_STATIC_ELF_AND_DEFAULT_PURE_PASS')
