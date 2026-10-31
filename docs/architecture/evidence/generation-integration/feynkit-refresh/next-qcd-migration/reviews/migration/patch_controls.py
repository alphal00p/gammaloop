import sys,re,copy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import common as C
import rows_review as R
import proof_review as J
require=R.require
TOKENS={C.CASES[0]:('0xb9fdc23ed056fa0c','0xe4b209849953f1fa'),C.CASES[1]:('0x5edf8e5f7a11eeb3','0x665e13a4c88da679'),C.CASES[2]:('0x2628c62dc373ac91','0x96e5abcebda25959')}

def exact_after(before,after):
 require(C.blob(before)=='7283198bf0152c82c1cd5f0173703cddfa7f5c0f','immutable-before-blob')
 fields,lines=R.parse_rows(before.decode());newfields,newlines=R.parse_rows(after.decode());expected=before
 reconstructed=C.load(C.O/'all12-Rows-independent.json');changes=[]
 for case,(old,new) in TOKENS.items():
  cr=next(r['Row'] for r in reconstructed if (r['case'],r['mode'],r['era'])==(case,'None','current'))
  require(newfields[case]==cr,'fresh-complete-Row-'+case)
  oldline=lines[case].encode();newliteral=('numerators: '+new).encode();oldliteral=('numerators: '+old).encode()
  require(oldline.count(oldliteral)==1,'unique-old-case-literal')
  new_line=oldline.replace(oldliteral,newliteral,1);expected=expected.replace(oldline,new_line,1)
  start=before.index(oldline)+oldline.index(oldliteral)+len('numerators: ')
  changes.append({'case':case,'old':old,'new':new,'start':start,'end_exclusive':start+18})
 require(after==expected,'exact-three-token-byte-scope')
 require([case for case in lines if lines[case]!=newlines[case]]==[case for case in lines if case in TOKENS],'exact-three-row-domain')
 require(lines['xs_qcd_gggg_sym_L1']==newlines['xs_qcd_gggg_sym_L1'],'same-hash-sibling-untouched')
 return fields,newfields,lines,newlines,changes

def apply_exact(before,patch):
 lines=patch.splitlines(keepends=True);path=C.MAIN.encode();header=[b'diff --git a/'+path+b' b/'+path+b'\n',b'index 7283198bf0152c82c1cd5f0173703cddfa7f5c0f..cc4505a6c4f6fc1687c075ed79e154a8e19f05a6 100644\n',b'--- a/'+path+b'\n',b'+++ b/'+path+b'\n']
 require(lines[:4]==header,'exact-patch-header');original=before.splitlines(keepends=True);out=[];cursor=0;i=4
 while i<len(lines):
  m=re.fullmatch(rb'@@ -(\d+),(\d+) \+(\d+),(\d+) @@\n',lines[i]);require(m is not None,'exact-hunk-header')
  start,n,dest,nn=map(int,m.groups());require(start-1>=cursor,'hunk-order');out.extend(original[cursor:start-1]);cursor=start-1
  require(len(out)==dest-1,'zero-offset-hunk');i+=1;left=right=0
  while i<len(lines) and not lines[i].startswith(b'@@'):
   tag,content=lines[i][:1],lines[i][1:];require(tag in (b' ',b'-',b'+'),'hunk-line')
   if tag in (b' ',b'-'):
    require(cursor<len(original) and original[cursor]==content,'zero-fuzz-context');cursor+=1;left+=1
   if tag in (b' ',b'+'):out.append(content);right+=1
   i+=1
  require((left,right)==(n,nn),'hunk-count')
 out.extend(original[cursor:]);return b''.join(out)

def main():
 before,after,patch=(C.read(C.P/p) for p in ('before.rs','after.NOT-APPLIED.rs','proposed-not-applied.patch'))
 old,new,ol,nl,spans=exact_after(before,after);trial=apply_exact(before,patch)
 require(trial==after and C.blob(trial)=='cc4505a6c4f6fc1687c075ed79e154a8e19f05a6','private-applied-blob')
 C.put('private-patch-application/'+C.MAIN,trial)
 prior=C.load(C.P/'preserved-prior13-full-Rows.json');require(len(prior)==13,'prior13-domain')
 for case,r in prior.items():require(ol[case].rstrip('\n')==r['line']==nl[case].rstrip('\n') and old[case]==r['complete_Row'],'prior-full-Row')
 # Also bind original previous12 from preceding packet, and latest prior grouped full Row.
 prior12=C.load(C.T/'grouped-identical-migration-preparation-20261031-r1/preserved-prior12-Rows.json')
 require(set(prior)==set(prior12)|{'xs_qed_eemm_identical_L2'},'prior13-source-domain')
 for case,r in prior12.items():require(ol[case].rstrip('\n')==r['row'],'original-prior12-source')
 controls=[]
 def positive(name,fn):fn();controls.append({'name':name,'result':'ACCEPT','credit':'pure-review-only'})
 def negative(name,predicate,fn):
  try:fn()
  except R.Rejected as e:
   require(str(e)==predicate,'intended-control-predicate '+name+': '+str(e));controls.append({'name':name,'result':'REJECT','predicate':str(e),'credit':'pure-review-only'})
  else:raise AssertionError('unexpected acceptance '+name)
 positive('exact immutable patch applies without fuzz or offset',lambda:require(apply_exact(before,patch)==after,'valid-patch'))
 positive('independent full Rows and complete byte scope',lambda:exact_after(before,after))
 negative('different source before','immutable-before-blob',lambda:exact_after(before+b'\n',after))
 sibling=nl['xs_qcd_gggg_sym_L1'];changed_sibling=sibling.replace(TOKENS[C.CASES[2]][0],TOKENS[C.CASES[2]][1])
 negative('same-old-hash sibling migration','exact-three-token-byte-scope',lambda:exact_after(before,after.replace(sibling.encode(),changed_sibling.encode(),1)))
 negative('one extra source byte','exact-three-token-byte-scope',lambda:exact_after(before,after+b'\n'))
 for case,(oldliteral,newliteral) in TOKENS.items():
  malformed=after.replace(nl[case].encode(),nl[case].replace(newliteral,'0x0000000000000000').encode(),1)
  negative('wrong reconstructed numerator '+case,'fresh-complete-Row-'+case,lambda malformed=malformed:exact_after(before,malformed))
 case=C.CASES[0];malformed=after.replace(nl[case].encode(),nl[case].replace('diagrams: 5','diagrams: 6').encode(),1)
 negative('selected nonnumerator field','fresh-complete-Row-'+case,lambda:exact_after(before,malformed))
 priorcase=next(iter(prior));bad=after.replace(nl[priorcase].encode(),nl[priorcase].replace('numerators: 0x','numerators: 0x0',1).encode(),1)
 negative('prior migrated full Row changed','exact-three-token-byte-scope',lambda:exact_after(before,bad))
 negative('wrong patch path','exact-patch-header',lambda:apply_exact(before,patch.replace(C.MAIN.encode(),b'wrong.rs',1)))
 negative('fuzzy source context','zero-fuzz-context',lambda:apply_exact(before.replace(b'    fast xs_qcd_gggg_sym_L1:',b'    slow xs_qcd_gggg_sym_L1:'),patch))
 inputs=C.load(C.R1/'prepared/interface-inputs.json')['cells'];native=C.load(C.R3/'execution/records/all-six-cells--typed8eb/results.json')['cells'];cell=inputs[0]
 worlds={e:C.load(C.CAP/f'prepared/captured/{cell["case"]}/{e}-{cell["mode"]}.json') for e in ('historical','current')}
 positive('exact captured model/request/Run cell',lambda:J.cell_input(cell,worlds))
 for key,value,pred in [('case',C.CASES[1],'world-identity'),('mode','Auto','world-identity')]:
  bad=copy.deepcopy(cell);bad[key]=value
  negative('wrong scientific '+key,pred,lambda bad=bad:J.cell_input(bad,worlds))
 bad=copy.deepcopy(cell);bad['model']['name']='review_wrong_model'
 negative('wrong model with otherwise identical selected inputs','proof-model',lambda:J.cell_input(bad,worlds))
 for key,value,pred in [('case',C.CASES[1],'live-case'),('mode','Auto','live-mode')]:
  bad=copy.deepcopy(worlds);bad['current']['live'][key]=value
  negative('wrong current live '+key,pred,lambda bad=bad:J.cell_input(cell,bad))
 g=cell['graphs'][0];rec=g['owners'][0];nr=native[0]['graphs'][0]['owners'][0];gid=g['id'];label=rec['owner']
 raw=lambda era:next(o['tensor'] for d in worlds[era]['raw']['diagrams'] if d['id']==gid for o in d['owners'] if o['owner']==label)
 h,c=raw('historical'),raw('current');public=R.owners(g['current_payload'])[label]
 check=lambda n:J.owner_input(rec,n,h,c,public)
 positive('actual owner with complete operand and output maps',lambda:check(nr))
 bad=copy.deepcopy(nr);bad['owner']='edge:2';negative('wrong actual owner','native-owner',lambda:check(bad))
 bad=copy.deepcopy(nr);bad['transport']=list(range(len(bad['transport'])));negative('replace nontrivial operand map with identity','native-operand-map',lambda:check(bad))
 bad=copy.deepcopy(nr);v=bad['historical']['canonical_to_logical'];v[0],v[1]=v[1],v[0];negative('valid bijection but wrong surviving output map','output-axis-map',lambda:check(bad))
 badrec=copy.deepcopy(rec);badrec['expression']=badrec['expression'].replace('spenso::{linear,spenso::tensor}::gamma','spenso::{spenso::tensor}::gamma')
 negative('same gamma head with changed full definition','owner-expression',lambda:J.owner_input(badrec,nr,h,c,public))
 m=C.load(C.N/'math-member-to-native-row.json')[0]['result'];bad=copy.deepcopy(m);bad['sha256']='0'*64;negative('wrong original result source hash','exact-reference-hash',lambda:J.resolve(bad))
 step=next(s for s in native[0]['graphs'][0]['actual_fold'] if s['physical_plan']['kind']=='Contract');previous=[];by={r['owner']:r for r in g['owners']}
 for s in native[0]['graphs'][0]['actual_fold']:
  if s is step:break
  previous=s['output']['ports']
 negative('corrupt actual integer plan despite valid named endpoints','integer-debug-plan',lambda:J.plan('Contract([])',step['physical_plan'],previous,by[step['owner']]['current']))
 fields=copy.deepcopy(step['physical_plan']);fields['names'][0][1]+='wrong'
 negative('corrupt named fold endpoint despite unchanged output count','named-endpoint-map',lambda:J.plan(step['current_plan'],fields,previous,by[step['owner']]['current']))
 C.save('fresh-controls.json',{'controls':controls,'total':len(controls),'accepted':sum(x['result']=='ACCEPT' for x in controls),'rejected':sum(x['result']=='REJECT' for x in controls),'native_or_scientific_negative_credit':0})
 C.save('patch-scope-independent.json',{'before_blob':C.blob(before),'after_blob':C.blob(after),'patch_sha256':C.sha(patch),'bytes':len(before),'Rows':149,'unchanged_full_Rows':146,'prior13_preserved':list(prior),'only_fields':['numerators'],'token_spans':spans,'actual_changed_byte_positions':[i for i,(x,y) in enumerate(zip(before,after)) if x!=y],'all_other_bytes_equal':True,'same_hash_sibling_preserved':True,'private_patch_fuzz':0,'private_patch_offset':0,'postapplication_runtime':'NOT_RUN','applied_to_tracked_source':False})
 C.flush('controls-consumed.json');print('PASS exact private patch and',len(controls),'fresh controls;',sum(x['result']=='REJECT' for x in controls),'intended rejections.')
if __name__=='__main__':main()
