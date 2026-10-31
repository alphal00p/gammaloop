"""Qualify actual new terminal dependency payloads against accepted raw-three output."""
import sys,tarfile,struct,re,hashlib,json,subprocess,traceback
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import pure
O=pure.O;W=pure.W;env=pure.ENV;sha=pure.sha;write=pure.write
def load(p):return json.loads(Path(p).read_text())
def ar_members(b):
 assert b[:8]==b'!<arch>\n';pos=8;out=[]
 while pos<len(b):
  h=b[pos:pos+60];assert h[-2:]==b'`\n';size=int(h[48:58]);raw=b[pos+60:pos+60+size]
  out.append((h[:16].strip().decode(),raw));pos+=60+size+(size%2)
 assert pos==len(b)
 names=next(data for name,data in out if name=='//') if any(name=='//' for name,data in out) else b''
 return [(names[int(name[1:]):names.index(b'/\n',int(name[1:]))].decode() if re.fullmatch(r'/\d+',name) else name,data) for name,data in out]
def elf_sections(b):
 assert b[:6]==b'\x7fELF\x02\x01'
 off=struct.unpack_from('<Q',b,40)[0];size,n,si=struct.unpack_from('<HHH',b,58)
 rows=[struct.unpack_from('<IIQQQQIIQQ',b,off+i*size) for i in range(n)];names=b[rows[si][4]:rows[si][4]+rows[si][5]]
 return {names[r[0]:names.index(b'\0',r[0])].decode():{'flags':r[2],'type':r[1],'sha256':hashlib.sha256(b[r[4]:r[4]+r[5]] if r[1]!=8 else b'').hexdigest()} for r in rows}
def inventory(path):
 p=subprocess.Popen(['/run/current-system/sw/bin/zstd','-dc',path],env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 out={};libs={}
 try:
  with tarfile.open(fileobj=p.stdout,mode='r|') as tf:
   for m in tf:
    if m.isfile():
     f=tf.extractfile(m)
     if m.name.endswith('.rlib'):
      b=f.read();h=hashlib.sha256(b);members=ar_members(b)
      libs[m.name.removeprefix('./')]=[{'name':name,'sha256':hashlib.sha256(data).hexdigest(),'special':data if name in ['//','lib.rmeta-link/'] else None,'sections':elf_sections(data) if name=='lib.rmeta-link/' else None} for name,data in members]
     else:
      h=hashlib.sha256()
      for b in iter(lambda:f.read(1<<20),b''):h.update(b)
     out[m.name.removeprefix('./')]=h.hexdigest()
 finally:
  p.stdout.close();err=p.stderr.read();rc=p.wait();assert rc==0,err
 return out,libs
def main():
 new=load(O/'archive-qualification.json');old=load(W/'.local/feynkit-refresh/raw-three-migration-20261031/archive-qualification.json')
 layers={'new':new['compiled_layers'][0],'accepted_raw_three':old['compiled_layers'][0]}
 inv={};libs={}
 for name,r in layers.items():assert sha(r['path'])==r['sha256'];inv[name],libs[name]=inventory(r['path'])
 linked=[p for p in inv['new'] if p.endswith(('.rlib','.rmeta','.so','.a','.o'))]
 assert linked
 qualifications=[]
 for label in ['accepted_raw_three']:
  for p in linked:
   if inv['new'][p]==inv[label].get(p):continue
   assert p.endswith('.rlib')
   newmembers,oldmembers=libs['new'][p],libs[label][p];assert len(newmembers)==len(oldmembers)
   mapping={};equal=0;special=[];changed_sections=[]
   for a,b in zip(newmembers,oldmembers):
    if a['name'] in ['//','lib.rmeta-link/']:
     assert a['name']==b['name'];special.append((a,b));continue
    assert a['sha256']==b['sha256'],(label,p,a['name'],'member-content-diff')
    equal+=1
    if a['name']!=b['name']:
     assert re.sub(r'\.[a-z0-9]+\.rcgu\.o$','.SESSION.rcgu.o',a['name'])==re.sub(r'\.[a-z0-9]+\.rcgu\.o$','.SESSION.rcgu.o',b['name'])
     assert len(a['name'])==len(b['name']);mapping[b['name']]=a['name']
   assert len(mapping)==len(set(mapping.values()))
   for a,b in special:
    mapped=b['special']
    for before,after in mapping.items():mapped=mapped.replace(before.encode(),after.encode())
    assert mapped==a['special'],(label,p,a['name'],'not-only-member-name-change')
    if a['sections']:
     assert a['sections'].keys()==b['sections'].keys()
     changed_sections=[k for k in a['sections'] if a['sections'][k]!=b['sections'][k]]
     assert all(not(a['sections'][k]['flags']&6) for k in changed_sections),(label,p,changed_sections)
   qualifications.append({'against':label,'library':p,'raw_archive_byte_equal':False,'identical_non_name_members':equal,'member_name_bijection':mapping,'metadata_name_changes_only':['GNU ar // string table','lib.rmeta-link/'],'changed_link_metadata_sections':changed_sections,'no_allocated_or_executable_section_difference':True})
 diffs={n:[p for p in inv['new'].keys()|v.keys() if inv['new'].get(p)!=v.get(p)] for n,v in inv.items() if n!='new'}
 write('dependency-invariance-probes.json',{'status':'PASS_WITH_EXPLICIT_NEW_TERMINAL_METADATA_QUALIFICATION','layers':layers,'member_sha256':inv,'terminal_dependency_artifacts_checked':linked,'count':len(linked),'whole_rlib_not_byte_identical':True,'new_terminal_library_qualification':qualifications,'all_member_differences':diffs,'meaning':'Exact unchanged whole dependency layers are reused byte-for-byte. Newly rebuilt terminal rlibs are newly qualified: all object payloads, symbol indexes and lib.rmeta are byte-exact; compiler-session rcgu member names differ only in the ar long-name table and nonallocated linker metadata. Explicit bijections verify every such byte. This is not whole-rlib byte identity; no capture payload normalization, donor execution or transplant.'})
 print('PASS',len(linked),'terminal dependency artifacts; rebuilt rlib objects and lib.rmeta byte-identical, explicit metadata-only member-name bijections')
try:main()
except BaseException:write('dependency-invariance-ERROR.json',{'traceback':traceback.format_exc()});raise
