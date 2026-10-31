"""Freeze immutable Git sources and hash the retained evidence without executing it."""
import context as C
import tomllib

def main():
    pins = {
        C.CAP/'handoff.json': '9da40bc2d34381c7c6b67fa81a6d91032ff699222349f6910f9787a18cc7e660',
        C.CAP/'SEAL.json': '166dd40ee056d9bf67da1a84eef200dc5bc9ff0fea55e15a88e9feb8f4e98d4b',
        C.R1/'SEAL.json': '6282d3457b22457e7465b89883db7dd4915967e0bbe6665011d0db3f48bbbe65',
        C.R2/'SEAL.json': 'f97aaf00c938d08940f22a8315358b3656b9b914eb8b968c95543fc356964ac3',
        C.R3/'SEAL.json': 'a7e3a1e8e4e8e6d979d8cf92f8a9de5878e88d6ac1f105bf19a86611e4381cc1',
        C.N/'report.json': '3baf87d1d031f790b753bb944d2a4c180f0620297b7808ad348dfa7ea65c76f0',
        C.N/'COMPLETE-SEAL.json': '22b11020019c4b9f8e058a2951de177f5a6f111fbb4f6976932f98ef0c399626',
        C.W/'.local/review-next-qcd-captures-astra-20261031-IYGt5Jrv/review.json': '7ccf8f133bff073e5f428d1fabc0ebd6552aa1e04a8b7edd149f526b37330aa7',
        C.W/'.local/review-next-qcd-proof-r2-supervisor-astra-20261031-mtFpBwVu/review.json': '64e9bf551fc949878e0c6647238fad4daf9eb69d82e0622939c89b138ca8d1d3',
        C.W/'.local/review-next-qcd-supervisor-r3-astra-20261031-22d61lkb/review.json': '6ad6242a78d52d0b371fb39d3942e66ed32bbbafec8878a91e6a72d5e3a9208b',
        C.REFRESH/'.local/feynkit-refresh/root-next-qcd-proof-native-execution-assignment-20261031-r3.json': '650f597bb727c64e53898d676f2884ee2125a5e73be8a8c0f0edf28376ba1382',
    }
    expected = {}
    manifests = []
    def add(p,h):
        p = str(p)
        C.need(p not in expected or expected[p] == h, 'consistent-original-binding', p)
        expected[p] = h
    for p,h in pins.items():
        C.need(C.sha(C.read(p)) == h, 'supplied-pin', p)
        add(p,h)
    for directory, name in [(C.CAP,'SHA256SUMS'), (C.R1,'SHA256SUMS'), (C.R2,'SHA256SUMS'), (C.R3,'SHA256SUMS'), (C.N,'COMPLETE-SHA256SUMS'), (C.N,'EXECUTION-SHA256SUMS')]:
        p=directory/name
        rows=C.manifest(p)
        manifests.append({'manifest':C.ref(p),'bindings':len(rows)})
        for q,h in rows: add(q,h)
    for p,r in C.load(C.N/'inputs-before.json').items(): add(p,r['sha256'])
    for directory in [C.W/'.local/review-next-qcd-captures-astra-20261031-IYGt5Jrv', C.W/'.local/review-next-qcd-proof-r2-supervisor-astra-20261031-mtFpBwVu', C.W/'.local/review-next-qcd-supervisor-r3-astra-20261031-22d61lkb']:
        for name in ('SHA256SUMS','MANIFEST.sha256'):
            p=directory/name
            if p.exists():
                rows=C.manifest(p); manifests.append({'manifest':C.ref(p),'bindings':len(rows)})
                for q,h in rows: add(q,h)
    originals = {}
    withheld = {}
    for p,h in sorted(expected.items()):
        if C.restricted(p):
            withheld[p]={'recorded_sha256':h,'fresh_read':False,'reason':'Explicit no config/license read boundary; external recorded source-config dependency only.'}
            continue
        r=C.fingerprint(p)
        C.need(r['sha256']==h,'original-input-hash',p)
        originals[p]=r
    for p in [C.W/'AGENTS.md', C.W/'CONTRIBUTING.typ', C.W.parent/'main/AGENTS.md', C.W.parent/'main/CONTRIBUTING.typ', C.W.parent/'generation-integration-plan/docs/architecture/generation-integration-plan.typ', C.W/'.gitignore', C.REFRESH/'.local/feynkit-refresh/root-next-qcd-capture-provenance-acceptance-20261031.json']:
        originals[str(p)]=C.ref(p) | {}
        originals[str(p)].pop('path')
    C.save('original-inputs-before.json', originals)
    C.save('recorded-unread-config-dependencies.json',withheld)
    C.save('original-manifests.json', manifests)
    C.save('supplied-pins.json',{str(p):h for p,h in pins.items()})
    C.save('refs-before.json',{'jj_ignore_working_copy':C.jj_state(),'git_refs':C.git('for-each-ref','--format=%(refname) %(objectname)').decode()})
    tree=C.git('rev-parse',C.TARGET+'^{tree}').decode().strip()
    C.need(tree=='c9de5c3fd962bd5d7c9375ab1f86f75699e8a69b','target-tree')
    selected=C.git('ls-tree','-r',C.TARGET,'crates/feynkit-generator/tests/generation_oracle').decode()
    names=[line.split('\t')[1] for line in selected.splitlines() if '7283198bf0152c82c1cd5f0173703cddfa7f5c0f' in line]
    C.need(names==[C.MAIN], 'immutable-expectation-path')
    freeze=[]
    for revision,label in [(C.TARGET,'target'),(C.CODE,'code'),(C.ACTUAL,'actual-pristine'),(C.HISTORICAL,'historical')]:
        for path in (C.MAIN,'crates/feynkit-generator/tests/generation_oracle/census.rs','crates/feynkit-generator/tests/generation_oracle/cases.rs','crates/feynkit-generator/tests/generation_oracle/classes.rs'):
            b=C.git('show',revision+':'+path)
            C.put('source-blobs/'+label+'/'+path,b)
            freeze.append({'revision':revision,'path':path,'packet_path':'source-blobs/'+label+'/'+path,'git_blob':C.blob(b),'sha256':C.sha(b),'bytes':len(b)})
    C.save('frozen-source-bindings.json',{'target':C.TARGET,'tree':tree,'relevant_code':C.CODE,'path_discovery_ls_tree':selected,'files':freeze})
    registry=C.git('show',C.TARGET+':docs/developers.toml')
    C.put('source-blobs/target/docs/developers.toml',registry)
    parsed=tomllib.loads(registry.decode())
    sources=[]
    for section in parsed['section']:
        for note in section.get('note',[]):
            path=note['source']; b=C.git('show',C.TARGET+':'+path)
            sources.append({'position':len(sources),'id':note['id'],'path':path,'git_blob':C.blob(b),'sha256':C.sha(b),'bytes':len(b),'registry_entry':note})
    C.need(len(sources)==118, 'all118-ordered-source-bindings',len(sources))
    C.save('ordered118-source-bindings.json',{'revision':C.TARGET,'registry_blob':C.blob(registry),'registry_sha256':C.sha(registry),'count':len(sources),'ordered_sources':sources,'candidate_registry_and_sources_changed':False})
    copies=[]
    for p in [C.R1/'pure_parser.py',C.R1/'strict_owners.py']:
        b=C.read(p);C.put('scripts/'+p.name,b);copies.append({'original':C.ref(p),'copy':'scripts/'+p.name,'copy_sha256':C.sha(b),'byte_equal':True,'executed_scope':'Pure syntax/owner/physical predicates only.'})
    C.save('copied-pure-source-bindings.json',copies)
    # Pattern packets are read as patterns only, with no theorem or runtime reuse.
    patterns=[]
    for packet,names in [('raw-three-migration-preparation-20261031-r2',('prepare_candidate.py','handoff.json','patch-binding.json')),('grouped-identical-migration-preparation-20261031-r1',('handoff.json','patch-binding.json','preserved-prior12-Rows.json'))]:
        for name in names: patterns.append(C.ref(C.T/packet/name))
    C.save('pattern-only-inputs.json',patterns)
    C.save('bootstrap-consumed-inputs.json',C.used)
    print('Frozen target; original inputs',len(originals),'unread source-config references',len(withheld),'ordered sources',len(sources))

if __name__=='__main__':
    main()
