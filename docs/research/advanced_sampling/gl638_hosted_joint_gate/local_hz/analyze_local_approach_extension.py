#!/usr/bin/env python3
"""Read retained local-approach JSON only; never loads a state or calls physics."""
from pathlib import Path
from decimal import Decimal as D, getcontext
import json, math, hashlib, sys
getcontext().prec=380
original=Path(sys.argv[1]); extension=Path(sys.argv[2]); out=Path(sys.argv[3]) if len(sys.argv)>3 else extension
folders=[original,extension]
frozen_path=Path(sys.argv[4]) if len(sys.argv)>4 else original.parent/'provenance/local-eight-frozen.json'
frozen=json.loads(frozen_path.read_text())
for filename,digest in frozen['files'].items():
 assert hashlib.sha256((original/Path(filename).name).read_bytes()).hexdigest()==digest,filename
budgets={'Double':D('1e-6'),'Quad':D('1e-10'),'Arb':D('1e-12')}
eps={'Double':D(2)**-52,'Quad':D(2)**-105,'Arb':D(2)**-999}
def num(v):
 d=D(str(v));assert d.is_finite(),v;return d
def complex_(v):return (num(v['re']),num(v['im']))
def norm(v):return (v[0]*v[0]+v[1]*v[1]).sqrt()
def plus(a,b):return(a[0]+b[0],a[1]+b[1])
def minus(a,b):return(a[0]-b[0],a[1]-b[1])
def scale(v,s):return(v[0]*s,v[1]*s)
def summed(vs):
 v=(D(0),D(0))
 for a in vs:v=plus(v,a)
 return v
def relative(a,b):
 n=max(norm(a),norm(b));return norm(minus(a,b))/n if n else D(0)
def strings(v):return {'re':str(v[0]),'im':str(v[1]),'norm':str(norm(v))}
def slope(xs,ys):
 if any(y<=0 for y in ys):return None
 x=[math.log(float(a)) for a in xs];y=[math.log(float(a)) for a in ys]
 xm=sum(x)/len(x);ym=sum(y)/len(y)
 return -sum((a-xm)*(b-ym) for a,b in zip(x,y))/sum((a-xm)**2 for a in x)
checks=[]; points=[]; failures=[]; hashes={}; raw_by_mode={}; raw_events_by_mode={}; source_points_by_mode={}
def check(kind,label,passed,**evidence):checks.append({'kind':kind,'label':label,'passed':bool(passed),**evidence})
comparison_by_mode={}
for path in sorted(p for folder in folders for p in folder.glob('*.approach.json')):
 hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest();row=json.loads(path.read_text());label=f"{row['mode']}:{row['name']}"
 if not row.get('complete'):
  failures.append({'label':label,'error':row.get('error'),'raw_normal_error':row.get('raw_normal_error'),'actual_forward_normal_error':row.get('actual_forward_normal_error')})
  continue
 record={'mode':row['mode'],'name':row['name'],'direction':row['archived_direction'],'archived_R_WH':row['archived_radius_GeV'],'comparison_R_map':row.get('comparison_native_normals',row['raw_normals']).get('R_map'),'raw_R_map':row['raw_normals'].get('R_map'),'actual_R_map':row['actual_forward_normals'].get('R_map'),'certified_rho':row['actual_forward_normals'].get('certified_rho'),'raw_normal_status':row['raw_normals'].get('status'),'source_normal_status':row['actual_forward_normals'].get('status')}
 vals={}
 for source_kind in ['raw_physics','source_physics']:
  pair={entry['forced_arb']:entry['result']['evaluation'] for entry in row[source_kind]}
  check('pair',label+':'+source_kind,set(pair)=={False,True})
  for forced,e in pair.items():
   tag=f'{label}:{source_kind}:{e["precision"]}:{forced}'
   precision=e['precision'];Y=complex_(e['complete_sample_estimator']);f=complex_(e['integrand_result']);J=num(e['parameterization_jacobian']) if e['parameterization_jacobian'] is not None else D(1);W=num(e['integrator_weight'])
   check('valid',tag,e['valid'])
   check('reported_factor_once',tag,relative(Y,scale(f,J*W)) <= eps[precision]*16,relative=str(relative(Y,scale(f,J*W))))
   events=e['events'];ids=sorted(ev['cut_info']['cut_id'] for ev in events)
   check('event_identity',tag,ids==list(range(6)) and all(ev['cut_info']['graph_id']==0 and ev['cut_info']['sampling_channel_id']==(row['channel_id'] if source_kind=='source_physics' else None) for ev in events),cuts=ids)
   total_events=summed(complex_(ev['weight']) for ev in events);rel=relative(total_events,Y)
   check('event_sum_vs_complete_estimator',tag,rel<=budgets[precision],relative=str(rel),budget=str(budgets[precision]))
   for ev in events:
    ct=ev['threshold_counterterms'];assert ct is not None
    terms=[complex_(ct['original']),*[complex_(c['weighted']) for c in ct['components']]]
    total=summed(terms);target=complex_(ev['weight']);absolute_error=norm(minus(total,target));magnitude=sum(norm(t) for t in terms);bound=eps[precision]*32*len(terms)*magnitude
    check('decomposition',tag+f':cut{ev["cut_info"]["cut_id"]}',absolute_error<=bound,relative_to_terms=str(absolute_error/magnitude if magnitude else D(0)),roundoff_bound_relative=str(eps[precision]*32*len(terms)),components=len(ct['components']))
   if forced:
    vals[source_kind]=Y
    if source_kind=='raw_physics':raw_events_by_mode.setdefault(row['name'],{})[row['mode']]={ev['cut_info']['cut_id']:ev for ev in events}
  ordinary=pair[False];arb=pair[True];rel=relative(complex_(ordinary['complete_sample_estimator']),complex_(arb['complete_sample_estimator']));tol=budgets[ordinary['precision']]+budgets['Arb']
  check('ordinary_vs_arb',label+':'+source_kind,rel<=tol,relative=str(rel),budget=str(tol),ordinary_precision=ordinary['precision'])
  record[source_kind+'_ordinary_vs_arb']=str(rel)
 # Canonical q and map factor below belong to the actual source cube. They are
 # diagnostics only; returned physical Y above is the estimator under study.
 jf=num(row['actual_forward']['selected_J_times_w']);q=num(row['actual_forward']['raw_map_density_sum']);W=num(row['grid_probes']['total']);F=scale(vals['source_physics'],D(1)/(jf*W));Fq=scale(F,D(1)/q)
 record.update(source_Y=strings(vals['source_physics']),source_F=strings(F),source_F_over_unweighted_qsum=strings(Fq),source_qsum=str(q),source_Jw=str(jf),source_outer_weight=str(W),raw_F=strings(vals['raw_physics']),comparison_qsum=row['raw_inverse']['raw_map_density_sum'],comparison_density_scope='At comparison_native_point; bare public raw F is at separately reported binary64 point, so raw F divided by this q is not treated as a same-point estimator')
 raw_by_mode.setdefault(row['name'],{})[row['mode']]=vals['raw_physics']
 native_point=[num(x) for x in row['actual_forward']['raw_coordinates']]
 original_point=[D.from_float(float(x)) for x in row.get('raw_binary64_point',row['original_binary64_tokens'])]
 comparison_point=[num(x) for x in row['comparison_native_point']] if 'comparison_native_point' in row else original_point
 comparison_by_mode.setdefault(row['name'],{})[row['mode']]=comparison_point
 record['source_max_coordinate_shift_from_comparison_native_GeV']=str(max(abs(a-b) for a,b in zip(native_point,comparison_point)))
 record['raw_binary64_max_shift_from_comparison_native_GeV']=str(max(abs(a-b) for a,b in zip(original_point,comparison_point)))
 record['source_R_relative_shift_from_comparison']=str(abs(num(record['actual_R_map'])/num(record['comparison_R_map'])-1))
 if 'extension_source' in row:
  ext=row['extension_source'];old=ext['retained_cube'];new=ext['cube'];divisor=int(ext['radial_divisor'])
  check('extension_unchanged_cube',label,all(float(a).hex()==float(b).hex() for i,(a,b) in enumerate(zip(old,new)) if i!=9) and len(old)==len(new)==12)
  radial_error=abs(D.from_float(float(new[9]))*divisor/D.from_float(float(old[9]))-1)
  check('extension_radial_boundary',label,radial_error<=D(2)**-52,relative=str(radial_error),divisor=divisor)
  check('extension_common_anchor',label,comparison_point==[num(x) for x in ext['actual_forward']['raw_coordinates']])
  check('extension_retained_reproduction',label,ext['retained_forward_reproduced_exactly'] is True and hashlib.sha256((original/Path(ext['artifact']).name).read_bytes()).hexdigest()==ext['sha256'])
  check('extension_same_compact_policy',label,ext['actual_forward_normals']['status']=='compact' and ext['actual_forward_normals']['normal_scale']==ext['retained_actual_forward_normals']['normal_scale']=='1' and ext['actual_forward_normals']['certified_rho']==ext['retained_actual_forward_normals']['certified_rho'])
  record['radial_divisor']=divisor
 record['source_max_coordinate_shift_from_original_binary64_GeV']=str(max(abs(a-b) for a,b in zip(native_point,original_point)))
 record['source_R_relative_shift_from_raw']=str(abs(num(record['actual_R_map'])/num(record['raw_R_map'])-1))
 record['source_F_relative_shift_from_raw']=str(relative(F,vals['raw_physics']))
 source_points_by_mode.setdefault(row['name'],{})[row['mode']]=native_point
 for route,key in [('raw_inverse','raw_inverse'),('actual_forward','actual_forward')]:
  mapped=row[key];jq=num(mapped['selected_J_times_w'])*num(mapped['raw_map_density_sum']);error=abs(jq-1)
  check('density_prefactor_consistency',label+':'+route,error<D('1e-13'),absolute_Jw_qsum_minus_one=str(error),scope='same represented map law: selected J*w times unweighted sum of raw map densities; canonical-Arb budget, no outer grid probability included')
 points.append(record)
for name,modes in comparison_by_mode.items():
 reference=next(iter(modes.values()))
 check('common_comparison_native_point',name,all(point==reference for point in modes.values()),modes=sorted(modes))
for name,modes in raw_by_mode.items():
 if 'optimized_lmb' not in modes:continue
 for mode,value in modes.items():
  if mode!='optimized_lmb':check('raw_bare_invariance',mode+':'+name,value==modes['optimized_lmb'],relative=str(relative(value,modes['optimized_lmb'])),comparison='exact native Arb retained values at original binary64 raw input')
for name,modes in raw_events_by_mode.items():
 if 'optimized_lmb' not in modes:continue
 for mode,events in modes.items():
  if mode=='optimized_lmb':continue
  for cut in range(6):
   for key in ['weight','threshold_counterterms']:
    check('raw_event_invariance',f'{mode}:{name}:cut{cut}:{key}',events[cut][key]==modes['optimized_lmb'][cut][key],comparison='exact native Arb raw cut weight/decomposition at the same original binary64 point')
source_differences=[]
for name,modes in source_points_by_mode.items():
 for ma in sorted(modes):
  for mb in sorted(modes):
   if ma>=mb:continue
   source_differences.append({'name':name,'modes':[ma,mb],'max_absolute_coordinate_difference_GeV':str(max(abs(a-b) for a,b in zip(modes[ma],modes[mb]))),'scope':'distinct actual source points due to independent inverse-cube binary64 boundary, not exact common-point physics comparisons'})
curves=[]
for mode in sorted({p['mode'] for p in points}):
 for direction in ['Hplus_Zminus','Hplus_Zplus']:
  ps=sorted([p for p in points if p['mode']==mode and p['direction']==direction],key=lambda p:num(p['actual_R_map']),reverse=True)
  if len(ps)!=6 or any(p['source_normal_status']!='compact' for p in ps):continue
  xs=[num(p['actual_R_map']) for p in ps];last=ps[-1]
  curves.append({'mode':mode,'direction':direction,'radii':[str(x) for x in xs],'last_complete_Y':last['source_Y'],'point_names':[p['name'] for p in ps],'Y_norms':[p['source_Y']['norm'] for p in ps],'last_qsum':last['source_qsum'],'last_F':last['source_F'],'last_to_previous_norm_ratio':str(num(last['source_Y']['norm'])/num(ps[-2]['source_Y']['norm'])),'slopes_last3':{f'{field}_{component}':slope(xs[-3:],[abs(num(p[field][component])) for p in ps[-3:]]) for field in ['source_Y','source_F','source_F_over_unweighted_qsum'] for component in ['re','im','norm']},'qsum_beta_last3':slope(xs[-3:],[num(p['source_qsum']) for p in ps[-3:]])})
report={'scope':'Independent Decimal380 postprocessing of retained native results only; slopes use actual forwarded R_map, not historical WH-scaled radii. Full Y uses actual selected source/grid factors. Local two-direction finite-range evidence, no global bound/MC convergence claim.','precision_epsilon_scope':'QuadFloat wraps pinned Symbolica DoubleFloat with106-bit precision: epsilon2^-105 (the original scan script used a stricter2^-112 rounding oracle; old result files are unchanged). Physical acceptance budgets unchanged.','original_frozen_files_verified':len(frozen['files']),'point_files':len(hashes),'complete_points':len(points),'failed_rows':failures,'all_checks_passed':all(c['passed'] for c in checks),'checks_passed':sum(c['passed'] for c in checks),'checks_failed':[c for c in checks if not c['passed']],'max_ordinary_vs_arb_norm_relative':max((D(c['relative']) for c in checks if c['kind']=='ordinary_vs_arb'),default=D(0)),'max_event_sum_norm_relative':max((D(c['relative']) for c in checks if c['kind']=='event_sum_vs_complete_estimator'),default=D(0)),'max_decomposition_error_relative_to_sum_terms':max((D(c['relative_to_terms']) for c in checks if c['kind']=='decomposition'),default=D(0)),'source_point_differences':source_differences,'points':points,'curves':curves,'checks':checks,'input_sha256':hashes}
out.mkdir(exist_ok=True);(out/'independent_extended_approach_analysis.json').write_text(json.dumps(report,indent=2,default=str)+'\n')
lines=['# Independent local H/Z approach analysis','',f'Completed point files: {len(hashes)}; complete rows: {len(points)}; retained failures: {len(failures)}. Native checks: {report["checks_passed"]}/{len(checks)} passed.','',report['scope'],'','Beta is −d log(abs(value))/d log(R); fits use the last three actual native forward radii. Complete source weights Y are in pb. qsum excludes discrete/grid probabilities, which are included once in Y.', '', '| Mode | Direction | beta F norm | beta qsum | beta Y norm | smallest-R Re Y [pb] | Im Y [pb] | Y norm [pb] |','|---|---|---:|---:|---:|---:|---:|---:|']
for c in curves:
 b=c['slopes_last3'];y=c['last_complete_Y'];lines.append(f"| {c['mode']} | {c['direction']} | {b['source_F_norm']:.6f} | {c['qsum_beta_last3']:.6f} | {b['source_Y_norm']:.6f} | {float(y['re']):.8g} | {float(y['im']):.8g} | {float(y['norm']):.8g} |")
lines+=['', 'Actual native R_map ranges [GeV]:', '']
for c in curves:lines.append(f"- {c['mode']}, {c['direction']}: {float(c['radii'][0]):.12g} → {float(c['radii'][-1]):.12g}.")
lines+=['',f"Largest actual-source radius displacement from the common native anchor: {max((D(p['source_R_relative_shift_from_comparison']) for p in points),default=D(0)):.6E} relative; this is recorded, not treated as identical-point physics.",f"Largest coordinate displacement from its intended binary64 raw target: {max((D(p['source_max_coordinate_shift_from_original_binary64_GeV']) for p in points),default=D(0)):.6E} GeV. Different maps recover slightly different binary64 cubes; actual forwarded points are retained and are not claimed identical.",f"Maximum same-map abs(J*w*qsum−1): {max((D(c['absolute_Jw_qsum_minus_one']) for c in checks if c['kind']=='density_prefactor_consistency'),default=D(0)):.6E}."]
lines+=['','Original24point files remain hash-identical. The four new anchors extend each direction by two decades; fits combine original and extended actual forwarded radii. Raw binary64 momenta, common native comparison anchors and each candidate actual forward point are distinct and retained; no raw F / common-anchor q coincidence is assumed. Quad rounding allowance follows the actual106-bit backend; physical stability budgets remain unchanged.',f"Maximum ordinary/Arb norm difference: {report['max_ordinary_vs_arb_norm_relative']:.6E}.",f"Maximum six-event sum / complete-estimator norm difference: {report['max_event_sum_norm_relative']:.6E}.",f"Maximum decomposition error relative to sum of term norms: {report['max_decomposition_error_relative_to_sum_terms']:.6E}.",'','Detailed check evidence and original-file hashes are in independent_extended_approach_analysis.json. Event/decomposition weights are already normalized by the precise result owner, including outer grid weight; they are not multiplied again.']
(out/'independent_extended_approach_analysis.md').write_text('\n'.join(lines)+'\n')
print(json.dumps({k:v for k,v in report.items() if k in ['point_files','complete_points','all_checks_passed','checks_passed','checks_failed','failed_rows','max_ordinary_vs_arb_norm_relative','max_event_sum_norm_relative','max_decomposition_error_relative_to_sum_terms']},default=str))
for c in curves:print(c['mode'],c['direction'],'betaF/q/Y',c['slopes_last3']['source_F_norm'],c['qsum_beta_last3'],c['slopes_last3']['source_Y_norm'],'lastY',c['last_complete_Y'])
