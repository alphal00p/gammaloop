from pathlib import Path
import hashlib,json,sys,datetime
O=Path(__file__).resolve().parent
W=O.parent.parent
P=W/'target/grouped-identical-migration-preparation-20261031-r1'
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def load(p):return json.loads(Path(p).read_text())
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
expected=load(P/'complete-input-hashes.json')
assert len(expected)==15899
# Only inert sealed repository/build/evidence files; never open private user settings or license files.
for p in expected:
    assert not p.startswith(('/home/','/root/','/etc/')),(p,'private-input prohibited')
    assert Path(p).name.lower() not in ['.env','.envrc','.netrc'],p
fail=[];observed={}
for s,h in expected.items():
    try:
        observed[s]=digest(s)
        if observed[s]!=h:fail.append({'path':s,'expected':h,'actual':observed[s]})
    except OSError as e:fail.append({'path':s,'error':str(e)})
save('input-hashes-'+sys.argv[1]+'.json',observed)
checks={}
pins={'handoff.json':'6279084309a00acfadede1c160cee1a2715f5b45a2e16e82dfba3470194b099c','seal.json':'49bfcf3c5b9bf9dd67552b5811498f2c3a80d12c1813c9aba973497ad0635b42','manifest.json':'2b767c3244605ab814c5dce3eab09fbe846bd10cc3b74b58b1dae7e45b5fe1f4','SHA256SUMS':'c91611445713af77e9bd8acd0c5de60f0a93318048fc88086761119dcee3c065','proposed-NOT-APPLIED.patch':'2450f56adfe23ea3b157b4d687c4ad6edd39866ff841f6ba8e85fa6f6d864fb7'}
for n,h in pins.items():assert digest(P/n)==h,(n,'pin');checks[n]=h
manifest=load(P/'manifest.json')['files'];assert len(manifest)==71
for n,r in manifest.items():assert digest(P/n)==r['sha256'] and (P/n).stat().st_size==r['bytes'],n
sums={line.split(maxsplit=1)[1].lstrip('*'):line.split()[0] for line in (P/'SHA256SUMS').read_text().splitlines()}
assert sums==({n:r['sha256'] for n,r in manifest.items()} | {'manifest.json':pins['manifest.json']})
assert set(p.relative_to(P).as_posix() for p in P.rglob('*') if p.is_file())==set(manifest)|{'manifest.json','SHA256SUMS','seal.json','seal.sha256'}
assert (P/'seal.sha256').read_text().split()[0]==pins['seal.json']
# The initial review pins and root decisions are independent inputs too.
root=W.parent/'generation-feynkit-refresh/.local/feynkit-refresh'
extra={str(root/n):digest(root/n) for n in ['root-grouped-none-acceptance-20261031.json','root-grouped-auto-acceptance-20261031.json','root-full-generator-pool1-acceptance-20261031.json']}
for folder in ['.local/review-grouped-none-aggregation-astra-20261030-independent-b1','.local/review-grouped-auto-aggregation-astra-20261031-independent-r1-1793413266607828149']:
    s=load(W/folder/'seal.json')
    for n,r in s['files'].items():
        p=W/folder/n
        assert digest(p)==r['sha256'],str(p)
    extra[str(W/folder/'seal.json')]=digest(W/folder/'seal.json')
registry={str(w/'docs/developers.toml'):digest(w/'docs/developers.toml') for w in [W,W.parent/'generation-feynkit-refresh',W.parent/'generation-integration-plan',W.parent/'main']}
save('integrity-'+sys.argv[1]+'.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'expected':len(expected),'verified':len(observed),'mismatches':fail,'candidate_ordinary_files':len(manifest),'candidate_pins':checks,'root_and_review_seals':extra,'registry':registry,'no_VCS_process_or_ref_index_operation':True,'private_config_reads':0,'configuration_note':'Complete frozen repository/config-source bytes are inert hash inputs only; no private license/environment/tool config is inspected or sourced.'})
print(json.dumps({'verified':len(observed),'mismatches':fail,'candidate_files':len(manifest),'extra':extra},indent=2))
