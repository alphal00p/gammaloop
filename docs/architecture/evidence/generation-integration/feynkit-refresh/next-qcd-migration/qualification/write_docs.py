"""Write an additive canonical execution record from actual completed qualification."""
import hashlib,json,os,re,subprocess,tomllib
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2]
def load(n):return json.loads((O/n).read_text())
def save(n,d):(O/n).write_text(json.dumps(d,indent=2)+'\n')
q=load('archive-qualification.json');r=load('runtime-results.json');j=json.loads((O/'joins/capture-joins.json').read_text());assert r['actual_executions']==r['passes']==48 and r['failures']==0 and j['actual_file_count']==15
note=W/'docs/architecture/generation-next-qcd-migration.typ';title='Next three QCD numerator migration qualification'
text=f'''= {title}
<generation-next-qcd-migration>
#quote(block: true)[
#strong[Status:] NEW independent integration review required. Root authorized
only the accepted three-token source/build/runtime/evidence phase. Global
`generation-integration` remains `83276dc047815d1dcc7f003aeb512c204b26758e`.
No self-approval, publication or whole-WP3 completion is claimed.
]

== Source and authorization

Code `{q['code']}` applies precisely the reviewed patch to immutable832.
`xs_qcd_gguu_sym_L1.numerators` changes from `0xb9fdc23ed056fa0c` to
`0xe4b209849953f1fa`; `xs_qcd_gguu_L1` changes from `0x5edf8e5f7a11eeb3` to
`0x665e13a4c88da679`; `xs_qcd_gggg_L1` changes from `0x2628c62dc373ac91` to
`0x96e5abcebda25959`. The sibling `xs_qcd_gggg_sym_L1` keeps its old token.
All149 Rows,146 unchanged full Rows,13 prior migrations and118 ordered source
bindings are preserved; the approved source blob is
`cc4505a6c4f6fc1687c075ed79e154a8e19f05a6`. Production algorithms, dependencies,
features, public APIs, defaults, budgets, stubs and CI false are unchanged.

Root authorization `91b7937e…` joins the independent migration review
`132a73c7…` to scientific/interface review `6d892bd9…` and root acceptance
`225b2c24…`. Historical PENDING strings and rejected preparations remain
immutable. The complete copied reviews, exact source/helpers, root receipts
and lossless reconstruction bindings are retained in the scoped receipt.

== Actual producer and execution

One pure Nix realization built the exact code above with one job/two cores,
normal exact-input substitution and all four license environment keys removed.
The separate six-root core union was untouched. Filtered source, full source,
NARs, derivations, compiler, lock/profile, recursive fingerprints and libraries
were qualified before private runtime handoff. Matching dependency stores were
reused; this is a newly produced archive, not the prior `ba02…` archive.

Archive SHA256: `{q['archive_sha256']}`.
Oracle ELF SHA256: `{q['oracle_elf']['sha256']}`.
Actual profile is `ci-optim`; features are `default,native`; Symbolica,
Numerica and Graphica resolve to kernel `8eb7c3729d6abefc916decc1bdf0c06536f8af80`.
The actual recursive root closure contains{q['closure_units']} units.

All16 approved complete-Row identities passed once at each pool1,8,32:
48 actual selected executions, zero retries, failures or selected-skip credit.
Discovery is retained separately and earns no execution credit. Commands,
streams, JUnit, exact selections, configuration, input before/after hashes,
wall times and available resource observations are retained. The qualified
capsule has no network, read-only source/ELF/libraries and protected records
outside its writable mount. The private license reaches immediate native exec
through an anonymous pipe; no value, value hash or private path is recorded.

The planned three pool1 executions produced15 actual files: each case's None
and Auto public payloads, both source companions and its census. Every full
byte equals the corresponding accepted pristine current capture. All six
cells,65/117/118 owners per case/mode, reconstruct their complete ten-field
Rows and join600 retained paired relations with exact case/mode/request,
128-bit string IDs, full keys and owner chains. OnlyNone asserts a literal
Row; Auto payload success and historical mode report/search counters remain
separate premises. No observer or scientific proof execution was repeated.

== Scientific and interface boundaries

The accepted finite proof retains578 complete-tree identities and22 even-f
cycle relations with144 exact rewrites,96 unexpanded products,1890 selected
pairs,456 actual folds and720 endpoint pairs. Both entering operand maps and
surviving output-axis maps are mandatory. The210 changed owners retain420
independently audited lowering layouts;390 equal imported-input source
arguments keep their deterministic-source premises. Rank8 numerator plus
rank4 projector leaves FOUR coad8 axes. Four of five sym gguu graphs have
factor2. The26 built-in bis no-ops earn no negative credit.

Historical producer `2b1ef4f…` uses Symbolica/Numerica `c78…` and registry
Graphica3.0.1. Pristine current source `ffa76ec…`/archive `8993…`, observer
built `5ea149ff…` plus qualified overlay and proof source `5ba04864…` remain
distinct from this new code/archive. Typed8eb imports both eras; it does not
supply historical c78 consumer execution. Missing historical interface metadata
means infer through its source, never an empty interface. Whole symbol
definitions, raw/cooked/decoded associations, callbacks/userdata, model/module
initialization, physical source factors and profile/compiler/store honesty
remain explicit finite premises. No positional identity or stale cached-plan
interchangeability is inferred.

INFO-1 remains excluded from general validator approval. The independent
exact actual integer-plan, dispatch and layout audits are retained and bound
to the identical actual outputs. No general decoder, locality/UV, full
integrand, denominator-evenness, trace/component execution, evaluator or
physical cross-section theorem is added.

== Durability and remaining gates

#link("evidence/generation-integration/feynkit-refresh/next-qcd-migration/receipt.json")[The receipt]
records exact identities, all restrictions and external dependencies. Committed
manifests are audited against actual immutable Git blobs, including explicitly
tracked ignored logs. Oversized evidence is losslessly compressed/split under
the unchanged1MiB cap. External original captures, source/store artifacts and
retained native records remain required; this packet is not self-contained.
The old425-file grouped manifest is checked at its exact corrected7f revision;
its registry source and the pre-edit832 registry snapshot remain explicit.
The failed7c candidate stays incomplete historical evidence.

Latest actual full-package results remain302passed/69failed/49ignored at the
prior source. These48 bounded runs do not update that result or certify full
Docs/Clippy, core/API/evaluator/host, WP3/M-rebuild, performance/memory policy
or publication. A NEW independent integration review and root disposition are
required before advancing integration832. Original evidence, failed attempts,
all useful artifacts and the active separate core work are retained.
'''
note.write_text(text)
registry=W/'docs/developers.toml';old=registry.read_bytes();assert old==(O/'registry-before-832.toml').read_bytes()
entry=f'''[[section.note]]
id = "generation-next-qcd-migration"
title = "{title}"
summary = "Exactly three accepted QCD numerator tokens; new pure archive, 48 complete-Row executions and six fresh public capture/proof joins."
source = "docs/architecture/generation-next-qcd-migration.typ"
kind = "Bounded migration qualification record"
status = "NEW independent integration review required; integration832 unchanged"
lifecycle = "investigation"
owner = "unassigned"
products = ["feynkit", "gammaloop"]
topics = ["generation", "testing", "numerators", "provenance"]
review_triggers = ["crates/feynkit-generator/tests/generation_oracle/**", "docs/architecture/generation-next-qcd-migration.typ", "docs/architecture/evidence/generation-integration/feynkit-refresh/next-qcd-migration/**"]

'''
reg=old.decode().replace('[[section.note]]',entry+'[[section.note]]',1);registry.write_text(reg)
a=tomllib.loads(old.decode());b=tomllib.loads(reg);removed=[n for sec in b['section'] for n in sec.get('note',[]) if n['id']=='generation-next-qcd-migration'];assert len(removed)==1
for sec in b['section']:sec['note']=[n for n in sec.get('note',[]) if n['id']!='generation-next-qcd-migration']
assert a==b
save('documentation-checks.json',{'title':title,'single_matching_title':re.findall(r'^= (.*)$',text,re.M)==[title],'registry_additive_only':True,'all118_ordered_scopes_preserved':True,'before_registry_revision':'83276dc047815d1dcc7f003aeb512c204b26758e','broader_docs_Clippy_credit':False})
t=load('typst-native-tools.json');env={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};env.update(TYPST_FONT_PATHS=t['fontPaths'],TYPST_PACKAGE_CACHE_PATH=t['packageCache'])
render={'compiler':t,'available':t['existing'],'compiled':False,'standalone_only':True,'license_absent':True}
if t['existing']:
 args=[t['compiler'],'compile','--root',str(W),str(note),str(O/'migration-note.pdf')];render['argv']=args
 with (O/'typst-compile.stdout').open('wb') as out,(O/'typst-compile.stderr').open('wb') as err:rc=subprocess.call(args,cwd=W,env=env,stdout=out,stderr=err)
 render['exit']=rc;assert rc==0,'Standalone Typst failure retained; do not claim render success'
 render.update(compiled=True,PDF_external={'path':str(O/'migration-note.pdf'),'sha256':hashlib.sha256((O/'migration-note.pdf').read_bytes()).hexdigest(),'bytes':(O/'migration-note.pdf').stat().st_size})
save('typst-render.json',render)
print('PASS additive note/registry and standalone Typst',render['compiled'])
