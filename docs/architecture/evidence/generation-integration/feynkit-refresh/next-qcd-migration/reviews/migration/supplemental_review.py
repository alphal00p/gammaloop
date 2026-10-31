import sys,re
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import common as C
# Attach complete accepted source and independent-audit hashes, not merely reported counts.
expected={}
def add(p,h):
 p=str(p);assert p not in expected or expected[p]==h,p;expected[p]=h
for rec in C.load(C.S/'supplemental-final-source-rehash.json'):add(rec['path'],rec['sha256'])
for p,r in C.load(C.S/'source-runtime-independent.json')['source_bindings'].items():add(p,r['sha256'])
for p,h in C.load(C.S/'parser-source-independent.json')['source_bindings'].items():add(p,h)
for p,r in C.load(C.S/'associations-independent.json')['unique_files'].items():add(p,r[0])
for p,h in expected.items():assert C.sha(C.read(p))==h,p
# All original historical reviews and root assignment are explicitly attached, with original decisions intact.
history=[]
for d,h in [('review-next-qcd-captures-astra-20261031-IYGt5Jrv','7ccf8f133bff073e5f428d1fabc0ebd6552aa1e04a8b7edd149f526b37330aa7'),('review-next-qcd-proof-preparation-astra-20261031-6xpiqaAK','62dade9062105a077ae6089d7e261742d57d63d9d07737b020c9b3098a71db95'),('review-next-qcd-proof-r2-supervisor-astra-20261031-mtFpBwVu','64e9bf551fc949878e0c6647238fad4daf9eb69d82e0622939c89b138ca8d1d3'),('review-next-qcd-supervisor-r3-astra-20261031-22d61lkb','6ad6242a78d52d0b371fb39d3942e66ed32bbbafec8878a91e6a72d5e3a9208b')]:
 p=C.W/'.local'/d/'review.json';r=C.load(p);assert C.sha(C.read(p))==h;history.append({'path':str(p),'sha256':h,'original_decision':r['decision']})
assignment=C.F/'root-next-qcd-proof-native-execution-assignment-20261031-r3.json';assert C.sha(C.read(assignment))=='650f597bb727c64e53898d676f2884ee2125a5e73be8a8c0f0edf28376ba1382'
# Independently inspect source factor evaluation and serializer feature declarations at immutable target.
for path in ['crates/feynkit-graph/src/expressions.rs','crates/feynkit-generator/Cargo.toml']:
 b=C.git('show',C.REV+':'+path);C.put('source-blobs/target/'+path,b)
C.save('accepted-source-audits-rehashed.json',{'files':len(expected),'bindings':expected,'history':history,'root_native_assignment_sha256':C.sha(C.read(assignment)),'meaning':'Accepted scientific/source premises retained on exact bytes. Historical pending/negative verdicts remain historical; no new native permission or credit.'})
C.flush('supplemental-consumed.json');print('PASS',len(expected),'accepted source/audit dependencies rehashed; four original review decisions preserved.')
