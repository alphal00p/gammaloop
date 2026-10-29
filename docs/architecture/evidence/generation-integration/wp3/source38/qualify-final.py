from pathlib import Path
import json,subprocess,hashlib,os
T=Path('target/wp3-source38');rev='0c65bf7c519764460a15a922340d958a9dfa7429';uri='git+file:///shared/lcnbr/dev/gammaloop/main?rev='+rev+'&shallow=1&allRefs=1'
r=subprocess.run(['nix','eval','--impure','--no-write-lock-file','--raw',uri+'#linnest-wasm.drvPath'],capture_output=True);(T/'runtime-eval.log').write_bytes(r.stderr);assert r.returncode==0,r.stderr.decode();drv=r.stdout.decode();assert drv=='/nix/store/8iq7qvzf38wcw165afivy6ipcw4sdlq9-linnest-wasm-0.1.0.drv',drv
raw=json.loads(subprocess.check_output(['nix','derivation','show',drv],stderr=subprocess.DEVNULL));d=next(iter(raw.get('derivations',raw).values()));src=Path(d['env']['src']);rows=[]
for p in sorted(src.rglob('*')):
 if p.is_file() or p.is_symlink():
  rel=str(p.relative_to(src));data=str(p.readlink()).encode() if p.is_symlink() else p.read_bytes();actual=subprocess.check_output(['git','--git-dir=../main/.git','show',rev+':'+rel]);assert data==actual,rel;rows.append(dict(path=rel,sha256=hashlib.sha256(data).hexdigest()))
assert len(rows)==1455;bundle=Path('/nix/store/88qva12y43c0fig08mslcrmjbb3ydxa9-linnest-wasm-0.1.0');assert bundle.is_dir();assert subprocess.check_output(['nix','hash','path',str(bundle)]).decode().strip()=='sha256-GBMop/R/i6hSu9DGsZlmo7yv7YR7uSdj9QtYhvrizZc=';payload=[]
for p in sorted(bundle.rglob('*')):
 if p.is_file():
  rel=str(p.relative_to(bundle));owner=Path('crates')/p.stem/'typst'/p.name if rel in ['kurvst.wasm','linnest.wasm'] else Path(rel.removeprefix('templates/')) if rel.startswith('templates/crates/') else Path('crates/clinnet')/rel;assert p.read_bytes()==owner.read_bytes(),rel;payload.append(dict(path=rel,owner=str(owner),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
assert len(payload)==29 and not any('tiptoe' in p['path'].lower() for p in payload)
(T/'runtime-source-final.json').write_text(json.dumps(dict(code=rev,input=uri,derivation=drv,source=str(src),source_NAR=subprocess.check_output(['nix','hash','path',str(src)]).decode().strip(),source_count=len(rows),sources=rows,output=str(bundle),output_NAR='sha256-GBMop/R/i6hSu9DGsZlmo7yv7YR7uSdj9QtYhvrizZc=',payload=payload,qualification='Actual combined Git-pinned source selects exactly approved runtime derivation, all1455 source bytes and29 payloads verified. Reuse registered exact producer, not payload-only or ancestry equivalence. No cold dependency determinism/full distribution claim.'),indent=2)+'\n');print('Exact runtime derivation/source/features qualify approved producer reuse:1455 inputs/29payloads')
