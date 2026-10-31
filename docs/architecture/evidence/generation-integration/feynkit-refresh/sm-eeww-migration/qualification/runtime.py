"""Fresh networkless capsule; no Nix calls. License enters only the exec child pipe."""
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

EXPORTED=os.environ.pop('SYMBOLICA_LICENSE',None)
os.environ.pop('SYMBOLICA_LICENSE_PATH',None)
os.environ.pop('SYMBOLICA_LICENSE_SIGNED',None)
O=Path(__file__).resolve().parent
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh')
Q=json.loads((O/'archive-qualification.json').read_text())
TOOLS=Q['runtime_tools']; PYTHON=Path(TOOLS['python']); NEXT=Path(TOOLS['nextest']); BWRAP=Path(TOOLS['bwrap'])
ARCH=O/'archive-extracted'; SOURCE=Path(Q['source']); INPUT=O/'capsule-input'
CASES=['xs_qed_eemm_L1','xs_qed_eedd_L1','xs_qed_eeee_1pi_L1','qed_eemm_L1','ctorder_qcd_uugg_L1','sm_eemm_L0','ym_g6_L0','xs_qed_eemm_L2','sm_eeww_L0']
FILTER='package(=feynkit-generator) and binary(=generation_oracle) and ('+' or '.join('test(='+n+')' for n in CASES)+')'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def write(p,data):
    Path(p).write_text(json.dumps(data,indent=2)+'\n')

# This unchanged parser performs no subprocesses, licensing or Nix operations.
spec=importlib.util.spec_from_file_location('elf_helpers',O/'elf_helpers.py');C=importlib.util.module_from_spec(spec);spec.loader.exec_module(C)

def closure():
    queue=[PYTHON,ARCH/'target/nextest/libdirs/host/libstd-64f5f36fb0927694.so',Path(Q['oracle_elf']['path']),NEXT]; seen=set();packages={PYTHON.parents[1]};records={};resolutions=[]
    while queue:
        p=Path(queue.pop()).resolve()
        if p in seen:continue
        seen.add(p);info=C.elf(p);records[str(p)]={'sha256':sha(p),**info}
        if str(p).startswith('/nix/store/'):packages.add(Path('/nix/store')/p.parts[3])
        if info.get('non_elf'):continue
        search=[]
        for d in info.get('runpath',[]):
            if not d:continue
            path=Path(d)
            assert path.is_absolute(),d
            if str(path).startswith('/nix/store/'):
                package=Path('/nix/store')/path.parts[3]
                if package.exists():packages.add(package)
            elif not path.is_relative_to(ARCH):assert not path.exists(),('unexpected host RUNPATH',d)
            search.append(path)
        interp=info.get('interpreter')
        if interp:queue.append(Path(interp));search.append(Path(interp).parent)
        # Exact archive library locations are the only non-store search roots.
        search += [ARCH/'target/nextest/libdirs/host',ARCH/'target/ci-optim/deps']
        for needed in info.get('needed',[]):
            hits=[d/needed for d in search if (d/needed).is_file()]
            if not hits:hits=[d/'lib'/needed for d in sorted(packages) if (d/'lib'/needed).is_file()]
            assert hits,(p,needed)
            hit=hits[0].resolve();assert str(hit).startswith('/nix/store/') or hit.is_relative_to(ARCH)
            resolutions.append({'requester':str(p),'needed':needed,'resolved':str(hit),'sha256':sha(hit)});queue.append(hit)
    return sorted(packages),records,resolutions

def prepare():
    assert (O/'nix-complete.json').is_file(),'all Nix/closure queries must precede this phase'
    assert sha(Q['archive'])==Q['archive_sha256'] and sha(Q['oracle_elf']['path'])==Q['oracle_elf']['sha256']
    INPUT.mkdir(exist_ok=False)
    config=(SOURCE/'.config/nextest.toml').read_bytes();assert b'[store]' not in config
    (INPUT/'nextest.toml').write_bytes(config+b'\n[store]\ndir = "/work/nextest"\n')
    # Only the existing store output location changes; profile, filters, timeouts,
    # retries, worker groups and assertions retain their immutable source bytes.
    (INPUT/'exec.py').write_text('''import os,sys
assert not any(k in os.environ for k in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED'))
assert os.statvfs('/source').f_flag & os.ST_RDONLY
assert os.statvfs('/archive').f_flag & os.ST_RDONLY
assert os.statvfs('/input').f_flag & os.ST_RDONLY
chunks=[]
while data:=os.read(0,8192):chunks.append(data)
null=os.open('/dev/null',os.O_RDONLY);os.dup2(null,0);os.close(null)
value=b''.join(chunks).decode().strip();assert value
env=dict(os.environ);env['SYMBOLICA_LICENSE']=value
os.execve(sys.argv[1],sys.argv[1:],env)
''')
    (INPUT/'preflight.py').write_text('''import os,json
assert not any(k in os.environ for k in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED'))
for p in ('/source','/archive','/input'):assert os.statvfs(p).f_flag & os.ST_RDONLY
for p in ('/run/nix/daemon-socket/socket','/home/lcnbr','/root','/shared'):assert not os.path.exists(p)
print(json.dumps({'netns':os.readlink('/proc/self/ns/net'),'readonly':['/source','/archive','/input'],'host_roots_credentials_nix_socket_absent':True}))
''')
    return closure()

def bwrap(root,packages):
    root.mkdir(exist_ok=False)
    for n in ['tmp','home','cache','state']:(root/n).mkdir()
    argv=[str(BWRAP),'--unshare-all','--die-with-parent','--new-session','--proc','/proc','--dev','/dev','--bind',str(root),'/work','--bind',str(root/'tmp'),'/tmp','--dir','/home','--bind',str(root/'home'),'/home/runtime','--ro-bind',str(INPUT),'/input','--ro-bind',str(ARCH),'/archive','--ro-bind',Q['archive'],'/archive.tar.zst','--ro-bind',str(SOURCE),'/source']
    for p in packages:argv+=['--ro-bind',str(p),str(p)]
    argv+=['--chdir','/source','--']
    return argv

def env(pool=1):
    workers=min(16,len(os.sched_getaffinity(0)))
    cores=int(os.environ.get('NIX_BUILD_CORES','0'))
    if cores>0:workers=min(workers,cores)
    return {'PATH':str(NEXT.parent),'HOME':'/home/runtime','USER':'runtime','LANG':'C.UTF-8','LC_ALL':'C.UTF-8','TMPDIR':'/work/tmp','XDG_CACHE_HOME':'/work/cache','GL_TEST_LOG_DIR':'/work/state','NEXTEST_TEST_THREADS':str(workers),'NEXTEST_NO_INPUT_HANDLER':'1','CARGO_NET_OFFLINE':'true','FEYNKIT_ORACLE_THREADS':str(pool)}

def run(name,command,packages,*,pool=1,licensed=False):
    root=O/name;argv=bwrap(root,packages)+command;e=env(pool)
    assert not any('LICENSE' in k for k in e)
    before={str(p):sha(p) for p in [Path(Q['archive']),Path(Q['oracle_elf']['path']),Path(__file__),*INPUT.iterdir()]}
    r={'argv':argv,'cwd':str(W),'environment_before_exec':e,'start':datetime.datetime.now(datetime.UTC).isoformat(),'license_delivery':'anonymous pipe -> isolated exec.py -> immediate native nextest; no license in argv/log/hash/file/Nix environment' if licensed else 'UNLICENSED','inputs_before':before,'mounts':'empty root; explicit library packages, exact new archive/full immutable source/input read-only; fresh output root only writable; separate network namespace; no Nix socket/host root/home credentials/cache/outlib'}
    write(root/'command.json',r);start=time.monotonic()
    value=None
    if licensed:
        value=EXPORTED
        if not value:
            config=Path(os.environ.get('XDG_CONFIG_HOME',str(Path.home()/'.config')))
            value=(config/'symbolica/license').read_text().strip()
        assert value
    with (root/'stdout').open('wb') as out,(root/'stderr').open('wb') as err:
        p=subprocess.Popen(argv,cwd=W,env=e,stdin=subprocess.PIPE if licensed else subprocess.DEVNULL,stdout=out,stderr=err)
        if licensed:p.stdin.write(value.encode());p.stdin.close();del value
        code=p.wait()
    after={p:sha(p) for p in before};assert after==before
    r.update(exit=code,end=datetime.datetime.now(datetime.UTC).isoformat(),seconds=time.monotonic()-start,inputs_after=after,stdout_sha256=sha(root/'stdout'),stderr_sha256=sha(root/'stderr'));write(root/'command.json',r)
    return root,code

def base():
    return [str(NEXT),'nextest']

def reuse():
    return ['--cargo-metadata','/archive/target/nextest/cargo-metadata.json','--binaries-metadata','/archive/target/nextest/binaries-metadata.json','--workspace-remap','/source','--target-dir-remap','/archive/target','--config-file','/input/nextest.toml','--user-config-file','none']

if __name__=='__main__':
    packages,elfs,resolutions=prepare();write(O/'runtime-closure.json',{'packages':list(map(str,packages)),'ELF_RUNPATH_NEEDED':elfs,'resolutions':resolutions,'tools':{k:{'path':v,'sha256':sha(v)} for k,v in TOOLS.items()},'Nix_queries_complete':True})
    hostnet=os.readlink('/proc/self/ns/net');p,rc=run('preflight', [str(PYTHON),'-I','-S','-B','/input/preflight.py'],packages);assert rc==0
    pre=json.loads((p/'stdout').read_text());assert pre['netns']!=hostnet;write(O/'network-preflight.json',{'host_netns':hostnet,'child':pre,'distinct_network_namespace':True})
    p,rc=run('inventory',base()+['list',*reuse(),'--message-format','json','-E',FILTER],packages);assert rc==0
    inv=json.loads((p/'stdout').read_text());found=[]
    for bid,suite in inv['rust-suites'].items():
        for name,t in suite['testcases'].items():
            if t['filter-match']['status']=='matches':found.append({'binary_id':bid,'name':name,**t})
    assert len(found)==9 and {r['name'] for r in found}==set(CASES),found
    write(O/'selected-inventory.json',{'actual_compiled_selected':found,'count':9,'filter':FILTER,'inventory_sha256':sha(p/'stdout'),'ignored_will_be_explicitly_run':True})
    print('actual compiled inventory: exact nine nonempty identities',flush=True)
    results=[]
    for pool in [1,8,32]:
        cmd=base()+['run',*reuse(),'--profile','ci_gammaloop','--run-ignored','all','--no-tests','fail','--no-fail-fast','--retries','0','--success-output','final','--status-level','all','--final-status-level','all','-E',FILTER]
        p,rc=run('pool'+str(pool),[str(PYTHON),'-I','-S','-B','/input/exec.py',*cmd],packages,pool=pool,licensed=True)
        junit=p/'nextest/ci_gammaloop/junit.xml';assert junit.is_file(),str(junit)
        tree=ET.parse(junit);cases=tree.findall('.//testcase');names=[c.attrib['name'] for c in cases];assert len(names)==9 and set(names)==set(CASES),names
        statuses=[]
        for c in cases:
            status='skipped' if c.find('skipped') is not None else 'failed' if c.find('failure') is not None or c.find('error') is not None else 'passed'
            statuses.append({'case':c.attrib['name'],'status':status,'time':c.attrib.get('time'),'classname':c.attrib.get('classname')})
        result={'pool':pool,'exit':rc,'outcomes':statuses,'actual_executions':sum(x['status']!='skipped' for x in statuses),'junit':str(junit),'junit_sha256':sha(junit),'command':str(p/'command.json')};results.append(result)
        write(O/'runtime-results.json',{'code':Q['code'],'archive_sha256':Q['archive_sha256'],'oracle_elf':Q['oracle_elf'],'results':results,'expected_selected_executions':27,'actual_executions':sum(x['actual_executions'] for x in results),'passes':sum(y['status']=='passed' for x in results for y in x['outcomes']),'failures':sum(y['status']=='failed' for x in results for y in x['outcomes']),'skip_credit':0,'proof_scope':'Prior finite SM proof remains pool1; pools8/32 are runtime evidence only.'})
        print('pool',pool,'exit',rc,[(r['case'],r['status']) for r in statuses],flush=True)
        if rc or any(r['status']!='passed' for r in statuses):raise RuntimeError('actual runtime divergence retained; do not weaken or rerecord')
    EXPORTED=None
    print('complete: 27 actual full-Row executions, zero retries or skip credit',flush=True)
