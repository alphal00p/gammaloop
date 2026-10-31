import json, os
from pathlib import Path
private = ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')
assert all(k not in os.environ for k in private)
assert os.readlink('/proc/self/ns/net') != os.environ['EXPECTED_HOST_NETNS']
readonly = ['/build/source','/proposal','/work/venv'] + json.loads(Path('/proposal/mounts.json').read_text())
for p in readonly: assert os.statvfs(p).f_flag & os.ST_RDONLY, p
assert not any(Path(p).exists() for p in ('/home/lcnbr','/run/opengl-driver','/nix/var/nix/daemon-socket/socket','/build/source/.git','/dev/dri','/run/dbus/system_bus_socket'))
interfaces = [l.split(':')[0].strip() for l in Path('/proc/net/dev').read_text().splitlines()[2:]]
assert interfaces == ['lo'], interfaces
assert len(Path('/proc/net/route').read_text().splitlines()) == 1
for root in Path('/nix/store').iterdir():
    for name in ('git','nix','cargo','rustc','cc','gcc','g++','c++','clang','clang++','tcc'):
        assert not (root/'bin'/name).exists(), str(root/'bin'/name)
print(json.dumps({'private_variables_absent':True,'distinct_network_namespace':True,'interfaces':interfaces,'no_routes':True,'readonly_inputs_checked':len(readonly),'compiler_vcs_nix_absent':True,'credentials_and_host_driver_tree_absent':True}))
