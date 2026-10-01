from pathlib import Path
from collections import deque
import datetime, hashlib, json, os, resource, shlex, signal, subprocess, time, tomllib
root=Path('/tmp/soft-ct-complexity-2026-09-14');repo=Path('/common/dev/gammaloop/lcnbr');binary=repo/'target/dev-optim/gammaloop'
original=repo/'tests/resources/run_cards/paper_figure_b1_double_triangle_soft_ir.toml'
card_text=original.read_text().replace('./tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot',str(repo/'tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot'))
card_text=card_text.replace('[cli_settings.global.generation.uv]\n','[cli_settings.global.generation.uv]\nfinal_integrand = "ThreeD"\nlocal_uv_cts_from_expanded_4d_integrands = false\nsubtract_uv = true\n').replace('[cli_settings.global.generation.evaluator]\n','[cli_settings.global.generation.evaluator]\ncompile = false\n')
card_text+='\n[cli_settings.global.generation]\nexplicit_orientation_sum_only = false\n'
card=root/'matched-single.toml';card.write_text(card_text);tomllib.loads(card_text)
all_text=card_text.replace('pat = "(0,0,-,-,-,+,-,-,+,-)"','').replace('explicit_orientation_sum_only = false','explicit_orientation_sum_only = true')
all_card=root/'matched-all.toml';all_card.write_text(all_text);tomllib.loads(all_text)
prescription={'log_divergent':'MUV','massive_power_divergent':'MUV','massless_power_divergent':'MUV','overrides':[]}
override=root/'muv.json';override.write_text(json.dumps({'global':{'generation':{'uv':{'renormalization_prescription':prescription}}}},indent=2)+'\n')
log_filter='off,[{#generation,#summary}]=debug,[{#generation,#cff,#profile}]=debug'
(root/'inputs.json').write_text(json.dumps({'original_card':str(original),'original_card_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'graph':str(repo/'tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot'),'graph_sha256':hashlib.sha256((repo/'tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot').read_bytes()).hexdigest(),'binary':str(binary),'binary_mtime_ns':binary.stat().st_mtime_ns,'single_orientation_pattern':'(0,0,-,-,-,+,-,-,+,-)','single_route':'localized direct 3D','all_route':'explicit direct 3D','reason_explicit_single_not_run':'GenerationSettings::validate_explicit_orientation_sum_options rejects explicit_orientation_sum_only=true with a nonempty orientation_pattern','limits':{'wall_seconds':180,'rss_tree_kib':12*1024*1024,'poll_seconds':0.2},'logging':{'GL_DISPLAY_FILTER':'off','GL_LOGFILE_FILTER':log_filter}},indent=2)+'\n')
results=[]
for label,scheme,selected_card in [('single-u','U',card),('single-h','H',card),('all-u-explicit','U',all_card)]:
    if label=='all-u-explicit' and any(x['exit_code']!=0 for x in results):
        (root/'all-u-explicit-skipped.json').write_text(json.dumps({'reason':'Single-orientation comparison did not complete.'},indent=2)+'\n');break
    folder=root/label;folder.mkdir(exist_ok=False)
    commands=(f'set global file {override}; ' if scheme=='U' else '')+f'run generate; save state --path {folder}/state'
    argv=[str(binary),str(selected_card),'-n','-s',str(folder/'state'),'-t',str(folder/'trace.jsonl'),'run','-c',commands]
    env=dict(os.environ,GL_DISPLAY_FILTER='off',GL_LOGFILE_FILTER=log_filter,GL_TEST_LOG_DIR=str(folder/'test-logs'))
    (folder/'command.json').write_text(json.dumps({'argv':argv,'shell_equivalent':shlex.join(argv),'cwd':str(root),'scheme':scheme,'input_card':str(selected_card)},indent=2)+'\n')
    start=time.monotonic();started=datetime.datetime.now(datetime.timezone.utc).isoformat();before=resource.getrusage(resource.RUSAGE_CHILDREN);peak_tree=0;peak_single=0;cap=None;cap_sample=None;last_telemetry=0
    print(f'START {label}',flush=True)
    with (folder/'stdout.log').open('w') as log,(folder/'resources.jsonl').open('w') as telemetry:
        proc=subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while proc.poll() is None:
            pending=deque([proc.pid]);seen=set();rss_tree=0;sample=[]
            while pending:
                pid=pending.popleft()
                if pid in seen:continue
                seen.add(pid)
                try:
                    fields={line.split(':',1)[0]:line.split(':',1)[1].strip() for line in Path(f'/proc/{pid}/status').read_text().splitlines() if ':' in line}
                    rss=int(fields.get('VmRSS','0 kB').split()[0]);hwm=int(fields.get('VmHWM','0 kB').split()[0]);rss_tree+=rss;peak_single=max(peak_single,hwm);sample.append({'pid':pid,'rss_kib':rss,'hwm_kib':hwm})
                    pending.extend(int(x) for x in Path(f'/proc/{pid}/task/{pid}/children').read_text().split())
                except (FileNotFoundError,ProcessLookupError):pass
            elapsed=time.monotonic()-start;peak_tree=max(peak_tree,rss_tree);row={'elapsed_seconds':elapsed,'rss_tree_kib':rss_tree,'processes':sample}
            if elapsed-last_telemetry>=1 or rss_tree>12*1024*1024 or elapsed>=180:telemetry.write(json.dumps(row)+'\n');telemetry.flush();last_telemetry=elapsed
            if rss_tree>12*1024*1024 or elapsed>=180:
                cap='rss_tree_exceeded_12_GiB' if rss_tree>12*1024*1024 else 'wall_exceeded_180_seconds';cap_sample=row
                try:os.killpg(proc.pid,signal.SIGTERM)
                except ProcessLookupError:pass
                try:proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    try:os.killpg(proc.pid,signal.SIGKILL)
                    except ProcessLookupError:pass
                    proc.wait()
                break
            time.sleep(.2)
        exit_code=proc.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result={'label':label,'scheme':scheme,'started_utc':started,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-start,'exit_code':exit_code,'termination_reason':cap,'cap_sample':cap_sample,'peak_sampled_tree_rss_kib':peak_tree,'peak_observed_process_hwm_kib':peak_single,'user_seconds':after.ru_utime-before.ru_utime,'system_seconds':after.ru_stime-before.ru_stime,'argv':argv,'cwd':str(root),'state_saved':(folder/'state/run.toml').exists()}
    (folder/'result.json').write_text(json.dumps(result,indent=2)+'\n');results.append(result);print(json.dumps(result),flush=True)
(root/'comparison-finished.json').write_text(json.dumps({'results':results},indent=2)+'\n')
