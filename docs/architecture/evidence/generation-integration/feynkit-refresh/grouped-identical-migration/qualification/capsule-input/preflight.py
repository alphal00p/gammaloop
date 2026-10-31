import os,json
assert not any(k in os.environ for k in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'))
for p in ('/source','/archive','/input'):assert os.statvfs(p).f_flag & os.ST_RDONLY
for p in ('/shared','/root','/home/lcnbr','/run/nix/daemon-socket/socket'):assert not os.path.exists(p)
print(json.dumps({'netns':os.readlink('/proc/self/ns/net'),'readonly':['/source','/archive','/input'],'host_roots_credentials_nix_socket_absent':True,'license_absent':True}))
