"""Private case/path adapter retains original strict owner and physical predicates."""

import collections
import copy
import pure_parser as B

class Reject(Exception):
 def __init__(self,stage,detail): self.stage=stage;super().__init__(detail)

def need(ok,stage,detail):
 if not ok: raise Reject(stage,detail)

def parse(text,mapped=False):
 try:return B.parse(text,mapped)
 except B.Reject as e:raise Reject(e.stage,str(e)) from e

def defs(tree,d=None):
 try:return B.definitions(tree,d)
 except B.Reject as e:raise Reject(e.stage,str(e)) from e

def equal(a,b,stage,detail):need(type(a)==type(b) and a==b,stage,detail)

def shape(a,b,path=''):
 need(type(a)==type(b),'schema',path)
 if isinstance(a,dict):
  need(a.keys()==b.keys(),'schema',path+' fields')
  for k in a:shape(a[k],b[k],path+'/'+k)
 elif isinstance(a,list):
  # Fixed occurrence domain; valid owner reorder leaves all homogeneous shapes valid.
  need(len(a)==len(b),'owner-domain',path+' length')
  for i,(x,y) in enumerate(zip(a,b)):shape(x,y,path+'/'+str(i))

def owner_table(w,route):
 p=w['payload'];raw=w['raw'];live=w['live']
 need(set(p)=={'diagrams','groups'},'schema','payload fields')
 ids=[d['id'] for d in p['diagrams']]; rm={d['id']:d for d in raw['diagrams']}
 need(len(ids)==len(set(ids))==len(rm)==len(raw['diagrams']) and set(ids)==set(rm),'owner-domain',route)
 need(live['grouping']=='none' and live['prepared_len']==len(ids),'selected-grouping',route)
 need(live['Run'].find('threads: Some(1)')>=0,'selected-pool',live['Run'])
 tables={}
 for i,g in enumerate(p['diagrams']):
  r=rm[g['id']];equal(r['name'],g['name'],'source-owner','raw/public name')
  po=B.owners(g['payload']);om={o['owner']:o['tensor'] for o in r['owners']}
  need(len(om)==len(r['owners']) and set(om)==set(po),'owner-domain',route+'/'+g['id'])
  for name,(expr,interface) in po.items():
   equal(expr,om[name]['decoded'],'payload-live-link',route+'/'+name)
   if interface is not None:equal(interface,om[name]['ports'],'interface-live-link',route+'/'+name)
   need(set(om[name])=={'raw','decoded','ports','is_metric','is_composite'},'schema','tensor fields')
  group={'master':i,'members':[{'id':g['id'],'name':g['name'],'overall_factor':g['overall_factor'],'ratio':'1','source_index':i}]}
  equal(p['groups'][i],group,'group-source-factor',route+'/'+g['id'])
  equal(r['json']['overall_factor'],g['payload']['overall_factor'],'source-factor','raw/public factor')
  equal(g['overall_factor'],g['payload']['overall_factor'],'source-factor','census/native factor')
  # Exact serialized model selected by original case; physical interaction and particle owner bindings.
  model=live['model'];selected=[]
  for vi,v in enumerate(g['payload']['vertices']):
   need(type(v['interaction'])==int and 0<=v['interaction']<len(model['vertex_rules']),'model-owner','interaction domain')
   rule=model['vertex_rules'][v['interaction']]
   incidence=collections.Counter()
   anti={a['name']:(a.get('antiname') or a['name']) for a in model['particles']}
   base=lambda x:min(x,anti[x])
   for ends,e in g['payload']['edges']:
    need(type(e['particle'])==int and 0<=e['particle']<len(model['particles']),'model-owner','particle domain')
    for side in ('source','target'):
     if ends[side]==vi:incidence[base(model['particles'][e['particle']]['name'])]+=1
   equal(incidence,collections.Counter(base(x) for x in rule['particles']),'model-owner','physical rule incidence')
   allowed={x for row in rule['couplings'] for x in row if x is not None}
   used={name.rsplit('::',1)[-1] for name in defs(parse(v['numerator'])) if name.startswith('UFO::') and name.rsplit('::',1)[-1].startswith('GC_')}
   need(used and used<=allowed,'model-source-factor','vertex coupling ownership')
   selected.append({'vertex':vi,'interaction':v['interaction'],'rule':rule['name'],'couplings':sorted(used)})
  tables[g['id']]=om
  w.setdefault('model_bindings',[]).append({'id':g['id'],'selected':selected})
 return tables

def physical(w,t):
 h,c=w['historical']['payload'],w['current']['payload']
 equal(h['groups'],c['groups'],'cross-source-factor','source/group records')
 need(len(h['diagrams'])==len(c['diagrams']),'cross-owner-domain','counts')
 for a,b in zip(h['diagrams'],c['diagrams']):
  equal({k:v for k,v in a.items() if k!='payload'},{k:v for k,v in b.items() if k!='payload'},'cross-physical','complete occurrence/key/cut/routing')
  x,y=copy.deepcopy(a['payload']),copy.deepcopy(b['payload'])
  for p in (x,y):
   for k in ('numerator','numerator_prefactor','projector','num_interface','numerator_prefactor_interface','projector_interface'):p.pop(k,None)
   for v in p['vertices']:
    v.pop('numerator');v.pop('num_interface',None)
   for e in p['edges']:
    need(type(e)==list and len(e)==2 and set(e[0])=={'source','target'},'edge-tuple','positional endpoint/data tuple')
    e[1].pop('numerator');e[1].pop('num_interface',None)
  equal(x,y,'cross-physical','complete physical owner records/model fingerprint')
  for owner in t['historical'][a['id']]:
   for flag in ('is_metric','is_composite'):equal(t['historical'][a['id']][owner][flag],t['current'][b['id']][owner][flag],'cross-owner-flags',owner)
