import json, os, sys
from pathlib import Path

assert all(k not in os.environ for k in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'))
assert os.readlink('/proc/self/ns/net') != os.environ['EXPECTED_HOST_NETNS']
assert all(os.statvfs(p).f_flag & os.ST_RDONLY for p in ('/build/source','/proposal','/work/venv'))
assert not any(Path(p).exists() for p in ('/home/lcnbr','/nix/var/nix/daemon-socket/socket','/build/source/.git','/run/opengl-driver'))
commands = json.loads(Path('/proposal/commands.json').read_text())
assert len(sys.argv) == 2 and sys.argv[1] in commands
command = commands[sys.argv[1]]
policy = json.loads(Path('/proposal/cwd-library-policy.json').read_text())
for directory in ('/work','/build/source','/build/source/crates/feynkit-py'):
    assert all(not Path(directory, name).exists() for name in policy['needed_names'])
parts = []
while block := os.read(0, 8192):
    parts.append(block)
private = b''.join(parts).decode()
assert private
null = os.open('/dev/null', os.O_RDONLY)
os.dup2(null,0)
os.close(null)
environment = dict(os.environ)
environment['SYMBOLICA_LICENSE'] = private
os.execve(command[0], command, environment)
