import os,sys
assert not any(k in os.environ for k in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'))
for p in ('/source','/archive','/input'):assert os.statvfs(p).f_flag & os.ST_RDONLY
parts=[]
while b:=os.read(0,8192):parts.append(b)
fd=os.open('/dev/null',os.O_RDONLY);os.dup2(fd,0);os.close(fd)
value=b''.join(parts).decode();assert value.strip()
e=dict(os.environ);e['SYMBOLICA_LICENSE']=value
os.execve(sys.argv[1],sys.argv[1:],e)
