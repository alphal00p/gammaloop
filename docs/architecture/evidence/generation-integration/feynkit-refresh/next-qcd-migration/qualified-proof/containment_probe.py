"""Trusted stdlib-only unlicensed containment body passed as Python -c text."""
import errno
import hashlib
import json
import os
from pathlib import Path
import socket
import sys

cfg = json.loads(sys.argv[1])
private_names = {'SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'}
assert not private_names.intersection(os.environ)
assert os.getcwd() == '/input'
assert set(os.listdir('/input')) == {'program','exec.py'}
assert not any(Path(p).exists() for p in ['/shared','/root','/home/lcnbr','/run/current-system','/nix/var/nix/daemon-socket/socket','/etc/ld.so.cache','/usr','/bin'])
mount_checks = []
for source, destination in cfg['readonly_mounts']:
    assert os.statvfs(destination).f_flag & os.ST_RDONLY, destination
    mount_checks.append(destination)
byte_checks = []
for path, expected in cfg['files'].items():
    p = Path(path)
    with p.open('rb') as stream:
        actual = hashlib.file_digest(stream,'sha256').hexdigest()
    assert actual == expected, path
    byte_checks.append({'path':path,'sha256':actual})
for path in ['/input/program','/input/exec.py',cfg['input']]:
    try:
        fd = os.open(path,os.O_WRONLY)
    except OSError as exc:
        assert exc.errno == errno.EROFS, (path,exc.errno)
    else:
        os.close(fd)
        raise AssertionError('readonly mount allowed write-open')
for directory in ['/work','/home/runtime','/work/tmp','/work/cache']:
    with (Path(directory)/'owned-output-probe').open('xb') as stream:
        stream.write(b'unlicensed trusted stdlib qualification\n')
    assert not os.statvfs(directory).f_flag & os.ST_RDONLY
netdev = Path('/proc/net/dev').read_text()
interfaces = [s.split(':')[0].strip() for s in netdev.splitlines()[2:] if ':' in s]
assert set(interfaces) <= {'lo'}, interfaces
with socket.socket(socket.AF_INET,socket.SOCK_STREAM) as probe:
    probe.settimeout(1)
    network_status = probe.connect_ex(('192.0.2.1',9))
assert network_status in [errno.ENETUNREACH,errno.EHOSTUNREACH], network_status
out = {'schema':'unlicensed-trusted-runtime-containment/1','status':'PASS',
       'id':cfg['id'],'readonly_mounts':mount_checks,'bytes':byte_checks,
       'all_four_private_names_absent':True,'hostroot_config_credentials_socket_compiler_cache_absent':True,
       'owned_output_roots':['/work','/home/runtime','/work/tmp','/work/cache'],
       'network_interfaces':interfaces,'network_probe_errno':network_status,
       'scientific_helper_executed':False,'private_handoff_qualified':False}
with Path('/work/probe.json').open('x') as stream:
    json.dump(out,stream,indent=2)
print(json.dumps({'id':cfg['id'],'status':'PASS','scientific_helper_executed':False}))
