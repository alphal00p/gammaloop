from pathlib import Path
from decimal import Decimal,localcontext
import json,hashlib,re,subprocess
out=Path('target/mark-geometry-round3');probe=out/'new-probe';rows=[]
with localcontext() as ctx:
 ctx.prec=100;D=Decimal.from_float
 def cross(a,b):return a[0]*b[1]-a[1]*b[0]
 def sub(a,b):return (a[0]-b[0],a[1]-b[1])
 def add(a,b):return (a[0]+b[0],a[1]+b[1])
 def mul(a,k):return (a[0]*k,a[1]*k)
 def norm(a):return (a[0]*a[0]+a[1]*a[1]).sqrt()
 for scale in [1e-6,1.,1e6]:
  for shape,l,w,i in [('triangle',10.,1e-7,0.),('triangle',10.,1e-9,0.),('stealth',10.,8.,.99999999),('triangle',10.,1e-5,0.),('stealth',10.,8.,.999999)]:
   args=[str(probe),shape,str(l*scale),str(w*scale),str(scale),str(i)];text=subprocess.check_output(args,text=True)
   line=next(l for l in text.splitlines() if l.startswith('path='));pts=[tuple(D(float(a)) for a in xy) for xy in re.findall(r'(?:MoveTo|LineTo)\(\(([^,]+), ([^)]+)\)\)',line)]
   outline=next(l for l in text.splitlines() if l.startswith('outline='));paint=[tuple(D(float(a)) for a in xy) for xy in re.findall(r'(?:MoveTo|LineTo)\(\(([^,]+), ([^)]+)\)\)',outline)]
   miter=D(float(re.search(r'miter=([^ ]+)',text)[1]));outer=[]
   for idx,p in enumerate(pts):
    u=sub(p,pts[idx-1]);v=sub(pts[(idx+1)%len(pts)],p);det=cross(u,v)
    if det<=0:continue
    n1=mul((u[1],-u[0]),D(.5*scale)/norm(u));n2=mul((v[1],-v[0]),D(.5*scale)/norm(v));a=add(p,n1);b=add(p,n2);t=cross(sub(b,a),v)/det;point=add(a,mul(u,t));ratio=norm(sub(point,p))/D(.5*scale);assert ratio<miter
    error=min(norm(sub(q,point)) for q in paint);condition=max(Decimal(1),norm(u)*norm(v)/abs(det));bound=64*D(2.220446049250313e-16)*condition*(norm(p)+norm(sub(point,p))+D(scale));assert error<=bound,(args,idx,error,bound)
    outer.append({'vertex':idx,'expected':[str(x) for x in point],'miter_ratio':str(ratio),'nearest_painted_vertex_error':str(error),'absolute_error_bound':str(bound),'edge_intersection_condition':str(condition),'relative_error':str(error/max(norm(point),D(scale)))})
   rows.append({'args':args,'native_stdout':text,'convex_miters':outer})
(out/'numeric-oracle.json').write_text(json.dumps({'digits':100,'method':'Independent exact-f64 emitted edge offset line intersections, checked against actual stroke path vertices; 64 epsilon times edge intersection condition times vertex + miter-displacement + stroke norms: accounts for rounding the f64 edge differences as well as the local join and cancellation at a nominal zero tip. No production-normal or half-angle algorithm copied.','cases':rows},indent=2)+'\n');print('15 exact-emitted-edge high-precision cases pass')
