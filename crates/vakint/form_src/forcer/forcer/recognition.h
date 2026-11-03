V <[P+p1],[P-p1]>,...,<[P+p11],[P-p11]>;
Autodeclare S M;
Autodeclare V p;
Autodeclare S x,n;
V Q;
CF vx,v(s),acc,EX(s),SYM,ID,TOPO,ex,ex1,replace,topo,vx,vy,vc,vxc,vxcs(s),dv,dx,vz,dy;
Set ppQ:Q,p1,...,p22;

#procedure loadTopologies()
#define PREFIX "FORCER"
#$MAXPROPS = 11;

* Initialize dictionary
#OpenDictionary wildmom
	#do i = 1,`$MAXPROPS'
		#add p`i': "p`i'?$p`i'"
	#enddo
#CloseDictionary

Skip;
L	`PREFIX'notation =
#include- notation.h
+TOPO(Md1)
 *SYM()
 *ex(Q,Q,Q,Q)
		;

#$numtopo = 0;

	$numtopo = $numtopo+1;
	id	TOPO(x?) = TOPO($numtopo,x);

ModuleOption noparallel;
.sort:load-topologies-1;
Skip;
NSkip `PREFIX'notation;
* Load the topologies into dollar variables and set
* TODO: split into multiple dollars according to vertex count?
	id	vx(?a) = acc(vx(?a));
	repeat id ex(?a,Q,?b) = ex(?a,?b); * strip Q
	id ex(?a) = acc(EX(?a));
	repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
	#do i = 1,`$numtopo'
		#$topoID`i' = 0;
		id  TOPO(`i',x1?)*ID(?a$topoID`i') = TOPO(`i',x1);
		id	TOPO(`i',x1?$toponame`i')*acc(x2?$topo`i') = TOPO(`i',x1)*acc(x2);
	#enddo
ModuleOption noparallel;
.sort:load-topologies-2;
Set `PREFIX'tnames:
	#do i = 1,`$numtopo'
		`$toponame`i''
	#enddo
;

* remove signs from reference topologies
#do i = 1,`$numtopo'
#$topofull`i' = $topo`i';
#inside $topo`i'
	repeat id vx(?a,-p?ppQ,?b) = vx(?a,p,?b);
	id	vx(?a) = v(?a);
	$counttopo`i' = count_(v,1);
#endinside
#enddo

Drop `PREFIX'notation;
.sort:load-topologies-3;
#endprocedure

#procedure matchTopologies(WARNNOMATCH)
#define PREFIX "FORCER"

#ifndef `$numtopo'
	#write "Error: call loadTopopologies first."
#endif

* remove signs from topologies to be matched
Multiply acc(1);
repeat id vx(?a) = v(?a)*acc(vx(?a));
id ex(?a) = ex1(?a)*ex(?a);
repeat id ex(?a,Q,?b) = ex(?a,?b); * strip Q
id ex(?a)*ex1(?b) = ex(?b)*acc(EX(?a));
repeat id acc(x1?)*acc(x2?) = acc(x1 * x2);
argument acc;
	repeat id vx(?a,-p?vector_,?b) = vx(?a,p,?b);
	id vx(?a) = v(?a); * make symmetric
endargument;

B+ acc;
.sort:recognition-1;
Keep brackets;

* match propagators
#$catch = -1;
#do i = 1,`$MAXPROPS'
	$p`i' = p`i';
#enddo
Argument acc;
	if ( count(EX,1) == 0 ) Multiply EX;
	$catch = 1;
	#UseDictionary wildmom($)
	#do i = 1,`$numtopo'
		if ( ( count(v,1) == `$counttopo`i'' ) && match(`$topo`i'') );
			$toponum = `i';
			goto caught;
		endif;
	#enddo
	#CloseDictionary
	$catch = 0;
	Label caught;
EndArgument;
if ( $catch );
	id acc(?a) = 1;
	Multiply replace(Q,Q,<$p1,p1>,...,<$p`$MAXPROPS',p`$MAXPROPS'>)*topo($toponum);
	id	topo(n?) = topo(n,`PREFIX'tnames[n]);
elseif ( $catch == 0 );
	#if `WARNNOMATCH'
		Print "WARNING: the following diagram did not match a topology: %t";
	#endif
	id acc(?a) = 1;
endif;
ModuleOption local,$catch,<$p1>,...,<$p{`$MAXPROPS'}>,$toponum;
.sort:recognition-2a;
id	replace(?a) = replace_(?a); * this will cause merges
B+ topo,v,ex;
.sort:recognition-2b;
Keep brackets;

* skip topologies that failed to match
if (count(topo,1) == 0) goto quit;

Multiply replace_(v,vy);

* get a copy of reference topology to fix the signs
#do i=1,`$numtopo'
	id topo(`i',n?) = topo(`i',n)*$topofull`i'*ID($topoID`i');
#enddo
id ID(0) = 1;
Multiply replace_(EX,ex);
id ex(?a) = 1;

* label every vertex and strip signs off a copy
Multiply acc(1);
repeat id vy(p?,?a)*acc(n?) = vy(n,p,?a)*vz(n,p,?a)*acc(n+1); * input
repeat id vx(p?,?a)*acc(n?) = vx(n,p,?a)*vxc(n,p,?a)*acc(n+1); * notation
repeat id vz(?a,-p?vector_,?b) = vz(?a,p,?b);
repeat id vxc(?a,-p?vector_,?b) = vxc(?a,p,?b);
id acc(n?) = 1;

Symmetrize vz;
Symmetrize vxc;

* match vertices in notation file to the input
repeat id vz(n1?,?a)*vxc(n2?,?a)*vy(n1?,?b)*vx(n2?,?c) = vy(n1,?b)*vx(n1,?c);
repeat id vy(n1?,?b,p?,?c)*vx(n1?,?d,-p?,?e) = vy(n1,?b,p,?c)*vx(n1,?d,?e)*replace(p);
id vx(?a) = 1;
id vy(n?,?a) = vx(?a);
id replace(-p?vector_) = replace(p);

label quit;
.sort:set signs-1;
id replace(p?)*replace(p?) = replace_(p,-p);
Multiply replace_(v,vx);
.sort:set signs-2;
#endprocedure

#procedure DotRewrite(DOTS)
* Assumes input expression is called `NAME'
* NOTE: all terms in `NAME' that do not have a valid topology
* label will be dropped
* Apply momentum substitution on all topologies.
* Rewrites all dots to a basis.
#define TOPOS "Mno2,Mno1,Mhaha,Mla4,Mj1,Mj6,Mj2,Mj5,Mj11,Mno3,\
    Mj4,Mno6,Mlala,Mnono,Mcross,Mbebe,Mbubu,Mlastar5,Mnostar0,Mnostar5,\
    Mnostar6,Mbustar5,Mfastar2,Mfastar3,Mt1star05,Mt1star24,Mt1star34,Mt1star45,Mt1star55,\
    <Md1>,...,<Md438>"

#define HASANYTOPO "0"
#do i={`TOPOS'}
	#define HASTOPO`i' "0"
	if (count(`i', 1));
		redefine HASTOPO`i' "1";
		redefine HASANYTOPO "1";
	endif;
#enddo

B+ `TOPOS';
.sort:DotRewrite-1;

#if `HASANYTOPO'

#do i={`TOPOS'}
	#if `HASTOPO`i''
		L DR`i' = `NAME'[`i']*`i';
	#endif
#enddo
ModuleOption noparallel;
.sort:DotRewrite-2;
Drop `NAME';
Hide
#do i={`TOPOS'}
	#if `HASTOPO`i''
		,DR`i'
	#endif
#enddo
;
.sort:DotRewrite-3;

#do i={`TOPOS'}
    #if `HASTOPO`i''
        UnHide DR`i';
        #call momsubs(`i',1,`DOTS')
        .sort:DotRewrite-`i';
        Hide DR`i';
    #endif
#enddo
Drop
#do i={`TOPOS'}
	#if `HASTOPO`i''
		,DR`i'
	#endif
#enddo
;

L `NAME' =
#do i={`TOPOS'}
	#if `HASTOPO`i''
		+DR`i'
	#endif
#enddo
;

ModuleOption noparallel;
.sort:DotRewrite-4;
#endif
#endprocedure

#procedure ToForcer()

* Fix ID
repeat id ID(p1?,p2?,p3?,p4?,?a) = ID(p1,p2)*ID(p3,p4,?a);
repeat id ID(-p1?vector_,p2?) = ID(p1,-p2);
id ID(p1?ppQ[n1?],p2?ppQ[n2?]) = theta_(n1-n2)*ID(p1,p2) + theta_(n2-n1)*ID(p2,p1);
id ID(p1?ppQ[n1?],-p2?ppQ[n2?]) = theta_(n1-n2)*ID(p1,p2) + theta_(n2-n1)*ID(p2,-p1);

* apply ID
id ID(?a) = replace_(?a);
repeat;
	id replace_(p1?vector_,p2?,?a) = replace_(p1,p2)*replace_(?a);
	id replace_(-p1?vector_,p2?,?a) = replace_(p1,-p2)*replace_(?a);
endrepeat;

id topo(n?,n2?) = n2;
id vx(?a) = 1;

#call DotRewrite(1)

#endprocedure
