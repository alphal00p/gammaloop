#ifndef `FORCERH'
#define FORCERH "1"
*--#[ Definitions :
*
*	This is the library file for the forcer system.
*	It contains declarations and procedures for the creation of the
*	reduction equations and the actual reduction program.
*	The actual programs for creating the reduction schemes are separate
*	files: one each for each master topology.
*	The major parts are the reduction schemes in the procedure FORCER.
*	The definitions of the names and notations are in the paper.
*
#define NUMEQ "0"
*#define POS ""
#define POS "pos_"
#define WITHPRESELECTION "1"
#define keep(n) "if ( ( match(Z(?a,1+`~n',?b)) == 0 ) && ( match(Z(?a,2+`~n',?b)) == 0 ) && ( match(Z(?a,3+`~n',?b)) == 0 ) ) Discard;"
#define keepm(n) "if ( ( match(Z(?a,-1+`~n',?b)) == 0 ) && ( match(Z(?a,-2+`~n',?b)) == 0 ) && ( match(Z(?a,-3+`~n',?b)) == 0 ) ) Discard;"
#define notp(n) "( ( match(Z(?a,1+`~n',?b)) == 0 ) && ( match(Z(?a,2+`~n',?b)) == 0 ) && ( match(Z(?a,3+`~n',?b)) == 0 ) )"
#define notm(n) "( ( match(Z(?a,-1+`~n',?b)) == 0 ) && ( match(Z(?a,-2+`~n',?b)) == 0 ) && ( match(Z(?a,-3+`~n',?b)) == 0 ) )"
*#define SKIP4 "if ( ( match(Z(?a,-1+n4,?b)) == 0 ) && ( match(Z(?a,-2+n4,?b)) == 0 ) && ( match(Z(?a,-3+n4,?b)) == 0 ) && ( match(Z(?a,-4+n4,?b)) == 0 ) && ( match(Z(?a,-5+n4,?b)) == 0 ) )"
#define SKIP48 "if ( ( match(Z(?a,-1+n4,?b)) == 0 ) && ( match(Z(?a,-2+n4,?b)) == 0 ) && ( match(Z(?a,-3+n4,?b)) == 0 ) && ( match(Z(?a,-1+n8,?b)) == 0 ) && ( match(Z(?a,-2+n8,?b)) == 0 ) && ( match(Z(?a,-3+n8,?b)) == 0 ) )"
#define SKIP4 "if ( ( match(Z(?a,-1+n4,?b)) == 0 ) && ( match(Z(?a,-2+n4,?b)) == 0 ) && ( match(Z(?a,-3+n4,?b)) == 0 ) && ( match(Z(?a,-4+n4,?b)) == 0 ) && ( match(Z(?a,-5+n4,?b)) == 0 ) \
	&& ( match(Z(?a,-1+n8,?b)) == 0 ) && ( match(Z(?a,-2+n8,?b)) == 0 ) && ( match(Z(?a,-3+n8,?b)) == 0 ) && ( match(Z(?a,-4+n8,?b)) == 0 ) && ( match(Z(?a,-5+n8,?b)) == 0 ) \
  && ( match(Z(?a,1+n10,?b)) == 0 ) && ( match(Z(?a,2+n10,?b)) == 0 ) && ( match(Z(?a,3+n10,?b)) == 0 ) && ( match(Z(?a,4+n10,?b)) == 0 ) && ( match(Z(?a,5+n10,?b)) == 0 ) )"
*#define SKIP4 ""
*#define SKIP48 ""
*#define SKIP57 "if ( ( match(Z(?a,-1+n5,?b)) == 0 ) && ( match(Z(?a,-2+n5,?b)) == 0 ) && ( match(Z(?a,-3+n5,?b)) == 0 ) && ( match(Z(?a,-1+n7,?b)) == 0 ) && ( match(Z(?a,-2+n7,?b)) == 0 ) && ( match(Z(?a,-3+n7,?b)) == 0 ) )"
*#define SKIP67 "if ( ( match(Z(?a,-1+n6,?b)) == 0 ) && ( match(Z(?a,-2+n6,?b)) == 0 ) && ( match(Z(?a,-3+n6,?b)) == 0 ) && ( match(Z(?a,-1+n7,?b)) == 0 ) && ( match(Z(?a,-2+n7,?b)) == 0 ) && ( match(Z(?a,-3+n7,?b)) == 0 ) )"
*#define SKIP67 "if ( ( match(Z(?a,-1+n6,?b)) == 0 ) && ( match(Z(?a,-2+n6,?b)) == 0 ) && ( match(Z(?a,-3+n6,?b)) == 0 ) && ( match(Z(?a,-4+n7,?b)) == 0 ) && ( match(Z(?a,-5+n7,?b)) == 0 ) && ( match(Z(?a,-3+n7,?b)) == 0 ) )"
*#define SKIP11 "if ( ( match(Z(?a,-1+n11,?b)) == 0 ) && ( match(Z(?a,-2+n11,?b)) == 0 ) && ( match(Z(?a,-3+n11,?b)) == 0 ) )"
*#define SKIP11 "if ( ( match(Z(?a,1+n11,?b)) == 0 ) && ( match(Z(?a,2+n11,?b)) == 0 ) && ( match(Z(?a,3+n11,?b)) == 0 ) )"
#define SKIP11 ""
#define SKIP68 "if ( ( match(Z(?a,-1+n6,?b)) == 0 ) && ( match(Z(?a,-2+n6,?b)) == 0 ) && ( match(Z(?a,-3+n6,?b)) == 0 ) && ( match(Z(?a,-1+n8,?b)) == 0 ) && ( match(Z(?a,-2+n8,?b)) == 0 ) && ( match(Z(?a,-3+n8,?b)) == 0 ) )"
#define SKIP6 "if ( ( match(Z(?a,-1+n6,?b)) == 0 ) && ( match(Z(?a,-2+n6,?b)) == 0 ) && ( match(Z(?a,-3+n6,?b)) == 0 ) )"
*#define SKIP6 ""
*#define SKIP5 "if ( ( match(Z(?a,-1+n5,?b)) == 0 ) && ( match(Z(?a,-2+n5,?b)) == 0 ) && ( match(Z(?a,-3+n5,?b)) == 0 ) )"
#define SKIP5 ""
*#define SKIP3 "if ( ( match(Z(?a,-1+n3,?b)) == 0 ) && ( match(Z(?a,-2+n3,?b)) == 0 ) && ( match(Z(?a,-3+n3,?b)) == 0 ) )"
#define SKIP3 ""
#define SKIP8 "if ( ( match(Z(?a,-1+n8,?b)) == 0 ) && ( match(Z(?a,-2+n8,?b)) == 0 ) && ( match(Z(?a,-3+n8,?b)) == 0 ) )"
*#define SKIP8 ""
*#define SKIP2 "if ( ( match(Z(?a,-1+n2,?b)) == 0 ) && ( match(Z(?a,-2+n2,?b)) == 0 ) && ( match(Z(?a,-3+n2,?b)) == 0 ) )"
#define SKIP2 ""
#define SKIP10 "if ( ( match(Z(?a,1+n10,?b)) == 0 ) && ( match(Z(?a,2+n10,?b)) == 0 ) && ( match(Z(?a,3+n10,?b)) == 0 ) )"
#define SKIP ""
#define SKIP1 ""
#define SKIP2 ""
#define SKIP3 ""
#define SKIP4 ""
#define SKIP5 ""
#define SKIP6 ""
#define SKIP7 ""
#define SKIP8 ""
#define SKIP9 ""
*
*--#] Definitions : 
*--#[ Declarations :
*
*	This file is a special copy of forcer2.h, used to redo the derivation
*	of no1(-n11) (in the file redo.frm) as it seems to contain at least one
*	error. It gives a relation for the master integral.
*	Because of its special application we have eliminated much code that
*	is superfluous for the rederivation.
*
#ifndef `N'
	#define N "14"
#endif
Symbols D,ep;
Dimension D;
Vectors p,Q,p0,...,p14,q1,...,q14;
Symbols n,n0,...,n14,k,k0,...,k9,c,x,x1,...,x14,plus,minus,a,a1,...,a20;
Symbols <[G1]>,...,<[G300]>;
Symbols m11,m12,m13,m21,m22,m23,m31,m32,m33,mdet,x21,x22,x23,x33,x31,x32;
CFunction L,C,C1,C2,C3,Master,Z,Y,ZZ,X,der,derp,acc,del,rat,RAT,frat,sign,rep,repp;
CFunction DER,num,NUM,den,DEN,HH,ZH;
Function CC,CC1;
Symbols lala,no1,no;
AutoDeclare Symbols y,h,x,u;
AutoDeclare CFunction rat;
AutoDeclare Index mu,nu;
Index MU,NU;
CFunction V3g,V4g,Dg;
Symbols xi,z1,z2,z3,z4,z5,z6,z7,z8,z9,z53;
*
*	The following symbols indicate topologies
*
*Symbols tr,l1,t1,t2,t3;
*Symbols no,la,be,bu,fa,o1,o2,o3,o4,y1,y2,y3,y4,y5;
*
*	4-loop master topologies
*
Symbols haha,no1,no2,no6,nono,lala,cross,bebe,fastar2;
Cfunction Zhaha,Zno1,Zno2,Zno6,Znono,Zlala,Zcross,Zbebe;
Cfunction Yhaha,Yno1,Yno2,Yno6,Ynono,Ylala,Ycross,Ybebe;
CFunction Yno1min1,Yno1min2,Yno1min3,Yno1min4;
#define NOCOMPLEXITY "0"
#define WITHCOMPLEXITY "1"
#define TOYNOCOMPLEXITY "2"
#define TOYWITHCOMPLEXITY "3"
#do i = 1,`N'
#ifdef `SPURIOUS'
	#define RAT`i' "rat`i'"
#else
	#define RAT`i' "rat"
#endif
#define POS`i' "pos_"
#define NEG`i' "neg0_"
#define extn`i' "?pos_"
#enddo
#ifndef `INSERT'
#define INSERT "0"
#endif
*
*--#] Declarations : 
.global
*--#[ FORCERD :
*
#procedure FORCERD(TOPO,Z,Y,InName,OutName)
*
*	This procedure has a reduction scheme for the various topologies.
*	It takes them to the point where other routines can take over.
*	The numerators are expressed in terms of dotproducts.
*
*	We work with the notation:
*	Z is the integral we are treating
*	Y is an integral that comes from Z but misses a line. Spectator.
*	OutName is the name of the output expression with the spectators.
*
#ifndef `NOSPECTATORS'
CreateSpectator x`TOPO' "x`TOPO'-`PID_'.spec";
#endif
*
#switch `TOPO'
*--#[ bustar5 :
* not a master topology
#case bustar5
#do il = 1,6
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering bustar5 loop `il', pass `$pass' at `time_' sec."
#endif

#if (`il' >= 1)
*
*  Prework to make the reduction of the bu topology faster.
*
* B6
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,-2*n3+2)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n3,1)+`Z'(n1,-1+n2,-1+n3,n4,1+n5,n6,n7,1+n8,n9)*rat(-ep-n5,1)
   +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-n3+1,1)+`Z'(n1,n2,-1+n3,-1+n4,1+n5,n6,n7,1+n8,n9)*rat(ep+n5,1)+`Z'(n1,n2,-1+n3,
   n4,n5,n6,n7,1+n8,n9)*rat(-3*ep-2*n2-n3-n5-n9+5,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9)*rat(-n9,1)+`Z'(n1,n2,n3,n4,n5,-1+
   n6,n7,1+n8,n9)*rat(-n3+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(-1+n3,1));
* B11
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?neg_,n9?neg0_) =
  -rat(1,-2*n7+2)*(`Z'(-1+n1,n2,n3,n4,n5,1+n6,-1+n7,1+n8,n9)*rat(n6,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n7,1)+`Z'(
   n1,-1+n2,n3,n4,1+n5,n6,-1+n7,1+n8,n9)*rat(ep+n5,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,-1+n7,1+n8,n9)*rat(-ep-n5,1)+`Z'(n1,n2,n3,-
   1+n4,n5,1+n6,-1+n7,1+n8,n9)*rat(-n6,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9)*rat(-n7+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+
   n8,n9)*rat(-n7+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(-3*ep-2*n4-n5-n6-n8+3,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,1+
   n9)*rat(n9,1));
* B8
id,ifmatch->1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,-2*n2+2)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n2,1)+`Z'(n1,-1+n2,-1+n3,n4,1+n5,n6,n7,1+n8,n9)*rat(-ep-n5,1)
   +`Z'(n1,-1+n2,n3,n4,1+n5,n6,-1+n7,1+n8,n9)*rat(ep+n5,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,1+n9)*rat(n9,2)+`Z'(n1,-1+n2,
   n3,n4,n5,n6,-1+n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,1+
   n8,n9)*rat(-3*ep-n2-2*n3-n5-n9+5,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+n8,n9)*rat(-n2+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,
   n9)*rat(-n2+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(-1+n2,1));
* B2
id,ifmatch->1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,-2*n1+2)*(`Z'(-1+n1,-1+n2,n3,1+n4,n5,n6,n7,n8,1+n9)*rat(-n4,1)+`Z'(-1+n1,n2,n3,1+n4,-1+n5,n6,n7,n8,1+n9)*rat(n4,1)+
   `Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-2*ep-n1-2*n2-n4-n8+5,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n1+1,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,
   n8,1+n9)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*rat(-n1+1,1));
* B1
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,2*n3-2)*(`Z'(n1,n2,-1+n3,n4,n5,-1+n6,1+n7,n8,1+n9)*rat(-n7,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,1+n7,n8,1+n9)*rat(n7,1)+`Z'(
   n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9)*rat(-8*ep-2*n1-2*n2-2*n3-2*n4-2*n5-2*n6-n7-n8-n9+13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,
   1+n9)*rat(-n3+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*rat(-1+n3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-1+n3,1));
* B4
id,ifmatch->1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,-2*n1+2)*(`Z'(-1+n1,n2,n3,n4,n5,1+n6,-1+n7,1+n8,n9)*rat(n6,1)+`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,1+n8,n9)*rat(-n6,1)+`Z'(-
   1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(8*ep+2*n1+2*n2+2*n3+2*n4+2*n5+n6+2*n7+n8+n9-13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,
   n9)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-n1+1,1));
* B5
id,ifmatch->1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,-2*n2+2)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n2+1,1)+`Z'(n1,-1+n2,-1+n3,n4,1+n5,n6,n7,n8,1+n9)*rat(ep+n5,1)+
   `Z'(n1,-1+n2,n3,n4,1+n5,n6,-1+n7,n8,1+n9)*rat(-ep-n5,1)+`Z'(n1,-1+n2,n3,n4,1+n5,n6,n7,n8,n9)*rat(-2*ep-2*n5,1)+`Z'(n1,-1+
   n2,n3,n4,n5,-1+n6,1+n7,n8,1+n9)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,2+n9)*rat(-n9-1,2)+`Z'(n1,-1+n2,n3,n4,n5,n6,
   -1+n7,n8,2+n9)*rat(1+n9,2)+`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,n8,1+n9)*rat(n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-5*
   ep-2*n1-n2-2*n4-n5-2*n6-n7-n8+9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,2+n9)*rat(-n9-1,2)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9)
   *rat(-1+n2,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,1+n9)*rat(-1+n2,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*rat(-n2+1,1));
* B12
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,-2*n4+2)*(`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,1+n8,-1+n9)*rat(2*ep+2*n5,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,1+n7,1+n8,n9)*rat(
   n7,1)+`Z'(n1,n2,n3,-1+n4,n5,1+n6,-1+n7,1+n8,n9)*rat(-n6,1)+`Z'(n1,n2,n3,-1+n4,n5,1+n6,n7,1+n8,n9)*rat(n6,1)+`Z'(n1,n2,n3,-
   1+n4,n5,n6,1+n7,1+n8,n9)*rat(-n7,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(n9,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9)
   *rat(n6-n7,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,2+n8,n9)*rat(-n8-1,1));
* B10
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,2*n6-2)*(`Z'(n1,-1+n2,n3,1+n4,n5,-1+n6,n7,n8,1+n9)*rat(-n4,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,1+n7,n8,1+n9)*rat(-n7,1)+
   `Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9)*rat(-n6+1,1)+`Z'(n1,n2,n3,1+n4,-1+n5,-1+n6,n7,n8,1+n9)*rat(n4,1)+`Z'(n1,n2,n3,n4,-1+
   n5,-1+n6,1+n7,n8,1+n9)*rat(n7,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9)*rat(-1+n6,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9)*
   rat(-n8,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,1+n9)*rat(4*ep+n4+2*n5+n7+n9-3,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*rat(-1+
   n6,1));
#endif

#if (`il' >= 1)
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
            -rat(1,-(n3-1)) * (
      +`Z'(1+n1,n2,(n3-1),-1+n4,1+n5,n6,n7,n8,n9)*rat(-n1,1)
      +`Z'(1+n1,n2,(n3-1),n4,1+n5,-1+n6,n7,n8,n9)*rat(n1,1)
      +`Z'(n1,n2,1+(n3-1),n4,1+n5,n6,-1+n7,n8,n9)*rat((n3-1),1)
      +`Z'(n1,n2,(n3-1),n4,1+n5,n6,n7,n8,n9)*rat(-6*ep-n1-2*n2-(n3-1)-2*n4-2*n5-n8-n9+6,1))
      ;
#endif

#if (`il' >= 2)
id,ifmatch->1,`Z'(n1?pos_,n2?{>1},1,n4?pos_,n5?,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
            -rat(1,(n2-1)) * (
      +`Z'(1+n1,(n2-1),1,-1+n4,1+n5,n6,n7,n8,n9)*rat(-n1,1)
      +`Z'(1+n1,(n2-1),1,n4,1+n5,-1+n6,n7,n8,n9)*rat(n1,1)
      +`Z'(n1,1+(n2-1),1,-1+n4,1+n5,n6,n7,n8,n9)*rat(-(n2-1),1)
      +`Z'(n1,(n2-1),1,n4,1+n5,n6,n7,n8,n9)*rat(-2*ep-n1-(n2-1)-2*n4-n8+4,1))
      ;
#endif

#if (`il' >= 3)
id,ifmatch->1,`Z'(n1?pos_,1,1,n4?pos_,n5?,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
            -rat(1,-10*ep-4*n1-3*n4-2*n5-2*n6-n7-2*(n8+1)-2*n9+14) * (
      +`Z'(n1,1,1,-1+n4,1+n5,n6,n7,(n8+1),-1+n9)*rat(ep+n5,1)
      +`Z'(n1,1,1,-1+n4,n5,1+n6,n7,(n8+1),n9)*rat(-n6,2)
      +`Z'(n1,1,1,-1+n4,n5,n6,1+n7,(n8+1),n9)*rat(-n7,2)
      +`Z'(n1,1,1,n4,n5,n6,n7,(n8+1),n9)*rat(-3*ep-n4-n5-n7-(n8+1)+3,1)
      +`Z'(-1+n1,1,1,n4,n5,1+n6,-1+n7,(n8+1),n9)*rat(n6,2)
      +`Z'(-1+n1,1,1,n4,n5,n6,n7,(n8+1),n9)*rat(16*ep+4*n1+4*n4+4*n5+2*n6+3*n7+2*(n8+1)+2*n9-19,2)
      +`Z'(n1,1,1,-1+n4,n5,-1+n6,1+n7,(n8+1),n9)*rat(n7,2)
      +`Z'(n1,1,1,-1+n4,n5,n6,n7,(n8+1),n9)*rat(n6-1,2)
      +`Z'(n1,1,1,n4,n5,-1+n6,n7,(n8+1),n9)*rat(-10*ep-4*n1-2*n4-2*n5-2*n6-n7-2*(n8+1)-2*n9+13,2)
      +`Z'(n1,1,1,n4,n5,n6,-1+n7,-1+(n8+1),1+n9)*rat(-n9,2)
      +`Z'(n1,1,1,n4,n5,n6,-1+n7,(n8+1),n9)*rat(13*ep+4*n1+4*n4+3*n5+3*n6+2*n7+3*(n8+1)+2*n9-18,2)
      +`Z'(-1+n1,1,1,n4,n5,n6,n7,1+(n8+1),n9)*rat((n8+1),2)
      +`Z'(n1,1,1,-1+n4,1+n5,n6,-1+n7,(n8+1),n9)*rat(ep+n5,2)
      +`Z'(n1,1,1,n4,n5,-1+n6,n7,1+(n8+1),n9)*rat(-(n8+1),2)
      +`Z'(n1,1,1,n4,n5,n6,-1+n7,(n8+1),1+n9)*rat(n9,2)
      +`Y'(n1,0,1,n4,1+n5,n6,n7,(n8+1),-1+n9)*rat(-ep-n5,1)
      +`Y'(n1,0,1,n4,n5,n6,1+n7,(n8+1),n9)*rat(n7,2)
      +`Y'(n1,1,0,n4,n5,n6,1+n7,(n8+1),n9)*rat(n7,2)
      +`Y'(n1,0,1,n4,n5,-1+n6,1+n7,(n8+1),n9)*rat(-n7,2)
      +`Y'(n1,0,1,n4,n5,n6,n7,(n8+1),n9)*rat(-8*ep-2*n1-2*n4-2*n5-2*n6-n7-(n8+1)-n9+10,2)
      +`Y'(n1,1,0,n4,n5,-1+n6,1+n7,(n8+1),n9)*rat(-n7,2)
      +`Y'(n1,1,0,n4,n5,n6,n7,(n8+1),n9)*rat(-8*ep-2*n1-2*n4-2*n5-2*n6-n7-(n8+1)-n9+10,2)
      +`Y'(n1,0,1,n4,1+n5,n6,-1+n7,(n8+1),n9)*rat(-ep-n5,2)
      +`Y'(n1,1,0,n4,n5,n6,n7,(n8+1),1+n9)*rat(-n9,2))
      ;
#endif

#if (`il' >= 4)
id,ifmatch->1,`Z'(n1?pos_,1,1,n4?pos_,n5?,n6?pos_,n7?pos_,0,n9?neg_) =
            -rat(1,2*ep+2*n5-2) * (
      +`Z'(n1,1,1,n4,-1+n5,n6,1+n7,0,(n9+1))*rat(-n7,1)
      +`Z'(-1+n1,1,1,1+n4,-1+n5,n6,n7,0,(n9+1))*rat(-n4,1)
      +`Z'(n1,1,1,1+n4,-1+n5,-1+n6,n7,0,(n9+1))*rat(n4,1)
      +`Z'(n1,1,1,n4,-1+n5,-1+n6,1+n7,0,(n9+1))*rat(n7,1)
      +`Z'(n1,1,1,n4,-1+n5,n6,n7,0,(n9+1))*rat(3*ep+n4+n5+2*n6+n7-5,1)
      +`Z'(n1,1,1,n4,-1+n5,n6,n7,0,1+(n9+1))*rat((n9+1),2)
      +`Z'(n1,1,1,n4,n5,n6,-1+n7,0,(n9+1))*rat(ep+n5-1,1)
      +`Z'(n1,1,1,n4,-1+n5,-1+n6,n7,0,1+(n9+1))*rat((n9+1),2)
      +`Z'(n1,1,1,n4,-1+n5,n6,-1+n7,0,1+(n9+1))*rat(-(n9+1),2)
      +`Y'(n1,1,0,n4,n5,n6,n7,0,(n9+1))*rat(-ep-n5+1,1))
      ;
#endif

#if (`il' >= 5)
id,ifmatch->1,`Z'(n1?pos_,1,1,n4?pos_,n5?,n6?pos_,n7?{>1},0,0) =
            -rat(1,-3*ep-n5+1) * (
      +`Z'(-1+n1,1,1,n4,1+n5,1+n6,(n7-1),0,0)*rat(-2*ep*n6-n1*n6-2*n4*n6+4*n6,(n7-1))
      +`Z'(-1+n1,2,1,-1+n4,1+n5,1+n6,(n7-1),0,0)*rat(-n6,(n7-1))
      +`Z'(1+n1,1,1,-2+n4,1+n5,1+n6,(n7-1),0,0)*rat(n1*n6,(n7-1))
      +`Z'(1+n1,1,1,-2+n4,1+n5,n6,1+(n7-1),0,0)*rat(n1,1)
      +`Z'(1+n1,1,1,-1+n4,1+n5,-1+n6,1+(n7-1),0,0)*rat(-n1,1)
      +`Z'(1+n1,1,1,-1+n4,1+n5,n6,(n7-1),0,0)*rat(3*ep*n1+2*n1*n4+n1*n5+n1*(n7-1)-4*n1,(n7-1))
      +`Z'(1+n1,1,1,n4,1+n5,-1+n6,(n7-1),0,0)*rat(-3*ep*n1-2*n1*n4-n1*n5-n1*n6-n1*(n7-1)+4*n1,(n7-1))
      +`Z'(n1,1,1,-1+n4,1+n5,1+n6,(n7-1),0,0)*rat(2*ep*n6+2*n4*n6-4*n6,(n7-1))
      +`Z'(n1,1,1,-1+n4,1+n5,n6,1+(n7-1),0,0)*rat(2*ep+n1+2*n4-5,1)
      +`Z'(n1,1,1,n4,1+n5,n6,(n7-1),0,0)*rat(6*ep^2+3*ep*n1+10*ep*n4+2*ep*n5+2*ep*n6+3*ep*(n7-1)-18*ep+2*n1*n4+n1*n5+2*n1*n6+n1
      *(n7-1)-4*n1+4*n4^2+2*n4*n5+2*n4*n6+2*n4*(n7-1)-14*n4+n5*(n7-1)-4*n5-4*n6-3*(n7-1)+12,(n7-1))
      +`Z'(n1,2,1,-2+n4,1+n5,1+n6,(n7-1),0,0)*rat(n6,(n7-1))
      +`Z'(n1,2,1,-2+n4,1+n5,n6,1+(n7-1),0,0)*rat(1,1)
      +`Z'(n1,2,1,-1+n4,1+n5,n6,(n7-1),0,0)*rat(4*ep+2*n4+2*n5+n6+(n7-1)-4,(n7-1))
      +`Y'(n1,2,0,n4,n5,n6,1+(n7-1),0,0)*rat(-1,1)
      +`Y'(n1,1,0,n4,1+n5,n6,1+(n7-1),0,0)*rat(-ep-n5,1))
      ;
#endif

#if (`il' >= 6)
id,ifmatch->1,`Z'(n1?pos_,1,1,n4?pos_,n5?,n6?pos_,1,0,0) =
      	-rat(1,24*ep^2+6*ep*n1+12*ep*n4+14*ep*n5-32*ep+2*n1*n5-2*n1+4*n4*n5-4*n4+2*n5^2-10*n5+8) * (
      +`Z'(1+n1,1,1,-1+n4,n5,n6,1,0,0)*rat(6*ep*n1+2*n1*n5-2*n1,1)
      +`Z'(1+n1,1,1,n4,n5,-1+n6,1,0,0)*rat(-6*ep*n1-2*n1*n5+2*n1,1)

      +`Z'(n1,2,1,-1+n4,n5,n6,1,0,0)*rat(2*ep,1)
      +`Y'(n1,1,2,n4,n5,n6,0,0,0)*rat(-2*ep,1))
      ;
#endif

GoTo 2;
Label 1;
$ired = 1;
Label 2;

id  `Z'(n1?neg0_,n2?,n3?,n4?,n5?,?a) = `Y'(n1,n2,n3,n4,n5,?a);
id  `Z'(n1?,n2?neg0_,n3?,n4?,n5?,?a) = `Y'(n1,n2,n3,n4,n5,?a);
id  `Z'(n1?,n2?,n3?neg0_,n4?,n5?,?a) = `Y'(n1,n2,n3,n4,n5,?a);
id  `Z'(n1?,n2?,n3?,n4?neg0_,n5?,?a) = `Y'(n1,n2,n3,n4,n5,?a);
id  `Z'(n1?,...,n5?,?a,n?neg0_,?b,n8?,n9?) = `Y'(n1,...,n5,?a,n,?b,n8,n9);

ModuleOption,maximum,$ired;
.sort:Main bustar5 loop, pass `$pass';
#if ( `$ired' > 0 )
      #redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
#break
*--#] bustar5 : 
*--#[ fastar2 :
* master topology
#case fastar2
#do il = 1,9
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering fastar2 loop `il', pass `$pass' at `time_' sec."
#endif

#if (`il' >= 1)
*
*  Prework to make the reduction of the fa topology faster.
*
* B5
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?neg_,n9?neg0_) =
  -rat(1,2*n7-2)*(`Z'(-1+n1,1+n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(-ep-n2,1)+`Z'(-1+n1,n2,n3,n4,n5,1+n6,-1+n7,1+n8,n9)*rat(n6,
   1)+`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,1+n9)*rat(-n9,2)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-n7+1,1)+`Z'(n1,-1+n2,n3,
   n4,n5,1+n6,-1+n7,1+n8,n9)*rat(-n6,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-n7+1,1)+`Z'(n1,1+n2,n3,n4,n5,-1+n6,-1+n7,
   1+n8,n9)*rat(ep+n2,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,1+n9)*rat(n9,2)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9)*rat(-1+n7,
   1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(-ep-n2+n6,1));
* B7
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,2*n6-2)*(`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,1+n8,n9)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-n6+1,1)+`Z'(
   n1,1+n2,-1+n3,n4,n5,-1+n6,n7,1+n8,n9)*rat(-ep-n2,1)+`Z'(n1,1+n2,n3,n4,n5,-1+n6,-1+n7,1+n8,n9)*rat(ep+n2,1)+`Z'(n1,n2,-1+
   n3,n4,n5,-1+n6,1+n7,1+n8,n9)*rat(n7,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,
   1+n8,n9)*rat(-n6+1,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,1+n8,1+n9)*rat(n9,2)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,1+n9)*rat(-n9,
   2)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9)*rat(-ep-n2+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(-1+n6,1));
* B2
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,2*n5-2)*(`Z'(-1+n1,n2,n3,n4,-1+n5,1+n6,n7,n8,1+n9)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n5+1,1)+`Z'(
   1+n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9)*rat(-n1,1)+`Z'(1+n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,1+n9)*rat(n1,1)+`Z'(n1,-1+n2,n3,n4,-
   1+n5,1+n6,n7,n8,1+n9)*rat(n6,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)
   *rat(-n5+1,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,1+n9)*rat(n8,2)+`Z'(
   n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9)*rat(-n1+n6,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,1+n9)*rat(-1+n5,1));
* B4
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,2*n6-2)*(`Z'(-1+n1,n2,n3,n4,1+n5,-1+n6,n7,n8,1+n9)*rat(-n5,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n6+1,1)+`Z'(
   1+n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,1+n9)*rat(n1,1)+`Z'(1+n1,n2,n3,n4,n5,-1+n6,n7,n8,1+n9)*rat(-n1,1)+`Z'(n1,n2,-1+n3,n4,n5,-
   1+n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,1+n8,1+n9)*rat(n8,2)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9)*
   rat(-1+n6,1)+`Z'(n1,n2,n3,n4,1+n5,-1+n6,n7,n8,1+n9)*rat(n5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,n2,
   n3,n4,n5,-1+n6,n7,n8,1+n9)*rat(-n1+n5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n6+1,1));
* B10
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,2*n4-2)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,1+n9)*rat(-n3,1)+
   `Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,1+n9)*rat(n7,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,-1+n2,n3,
   n4,n5,n6,n7,n8,1+n9)*rat(-n4+1,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,1+n7,n8,1+n9)*rat(-n7,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+
   n9)*rat(-n4+1,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,-1+n7,n8,1+n9)*rat(n3,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,1+n8,1+n9)*rat(n8,2)+
   `Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,1+n9)*rat(-n3+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*rat(-1+n4,1));
* B12
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?neg0_,n9?neg_) =
  -rat(1,2*n7-2)*(`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,-1+n7,n8,1+n9)*rat(-n4,1)+
   `Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9)*rat(-n7+1,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,-1+n7,n8,1+n9)*rat(n3,1)+`Z'(n1,n2,1+n3,n4,n5,
   n6,-1+n7,n8,1+n9)*rat(-n3,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,1+n9)*rat(-1+n7,1)+`Z'(n1,n2,n3,1+n4,n5,n6,-1+n7,n8,1+n9)*
   rat(n4,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,1+n9)*rat(n8,2)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,n2,
   n3,n4,n5,n6,-1+n7,n8,1+n9)*rat(-n3+n4,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n7+1,1));
* B3
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,-2*n5+2)*(`Z'(1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9)*rat(-2*n1,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,1+n8,n9)*rat(-n7,1)+`Z'(
   n1,1+n2,-1+n3,n4,-1+n5,n6,n7,1+n8,n9)*rat(-ep-n2,1)+`Z'(n1,1+n2,n3,n4,-1+n5,n6,-1+n7,1+n8,n9)*rat(ep+n2,1)+`Z'(n1,n2,-1+
   n3,n4,-1+n5,n6,1+n7,1+n8,n9)*rat(n7,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,
   n7,2+n8,n9)*rat(-n8-1,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n5,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,n7,1+n8,1+n9)*rat(
   n9,2)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9)*rat(-n5+1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,n2,n3,n4,-
   1+n5,n6,n7,1+n8,n9)*rat(-ep-n2+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n5,1));
* B9
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,-2*n4+2)*(`Z'(-1+n1,1+n2,n3,-1+n4,n5,n6,n7,1+n8,n9)*rat(-ep-n2,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,1+n6,n7,1+n8,n9)*rat(n6,
   1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,2+n8,n9)*rat(-n8-1,1)+`Z'(-1+n1,n2,
   n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n4,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,1+n6,n7,1+n8,n9)*rat(-n6,1)+`Z'(n1,1+n2,n3,-1+n4,n5,-1+n6,
   n7,1+n8,n9)*rat(ep+n2,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,n9)*rat(-2*n3,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,n7,1+n8,1+n9)*rat(
   n9,2)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9)*rat(-ep-n2+n6,1)+`Z'(n1,n2,n3,
   n4,-1+n5,n6,n7,1+n8,n9)*rat(-n4+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n4,1));
#endif

#if (`il' >= 1)
id,ifmatch->1,`Z'(n1?pos_,n2?,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
            -rat(1,(n5-1)) * (
      +`Z'(n1,-1+n2,n3,n4,(n5-1),1+n6,n7,n8,n9)*rat(n6,1)
      +`Z'(-1+n1,n2,n3,n4,1+(n5-1),n6,n7,n8,n9)*rat(-(n5-1),1)
      +`Z'(-1+n1,n2,n3,n4,(n5-1),1+n6,n7,n8,n9)*rat(-n6,1)
      +`Z'(n1,n2,n3,n4,(n5-1),n6,n7,n8,n9)*rat(-2*ep-2*n1-(n5-1)-n6-n8+4,1))
      ;
#endif

#if (`il' >= 2)
id,ifmatch->1,`Z'(n1?pos_,n2?,n3?pos_,n4?{>1},1,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
            -rat(1,(n4-1)) * (
      +`Z'(n1,-1+n2,n3,(n4-1),1,n6,1+n7,n8,n9)*rat(n7,1)
      +`Z'(n1,n2,-1+n3,1+(n4-1),1,n6,n7,n8,n9)*rat(-(n4-1),1)
      +`Z'(n1,n2,-1+n3,(n4-1),1,n6,1+n7,n8,n9)*rat(-n7,1)
      +`Z'(n1,n2,n3,(n4-1),1,n6,n7,n8,n9)*rat(-2*ep-2*n3-(n4-1)-n7-n8+4,1))
      ;
#endif

#if (`il' >= 3)
id,ifmatch->1,`Z'(n1?pos_,n2?,n3?pos_,1,1,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
            -rat(1,-8*ep-2*n1-2*n2-2*n3-2*n6-2*n7-2*n8-(n9+1)+11) * (
      +`Z'(1+n1,-1+n2,n3,1,1,n6,n7,n8,(n9+1))*rat(n1,2)
      +`Z'(n1,-2+n2,n3,1,1,1+n6,n7,n8,(n9+1))*rat(-n6,2)
      +`Z'(n1,-2+n2,n3,1,1,n6,1+n7,n8,(n9+1))*rat(-n7,2)
      +`Z'(n1,-1+n2,1+n3,1,1,n6,n7,n8,(n9+1))*rat(n3,2)
      +`Z'(n1,-1+n2,n3,1,1,1+n6,n7,n8,(n9+1))*rat(-n6,2)
      +`Z'(n1,-1+n2,n3,1,1,n6,1+n7,n8,(n9+1))*rat(-n7,2)
      +`Z'(-1+n1,-1+n2,n3,1,1,1+n6,n7,n8,(n9+1))*rat(n6,1)
      +`Z'(-1+n1,n2,n3,1,1,1+n6,n7,n8,(n9+1))*rat(n6,2)
      +`Z'(-1+n1,n2,n3,1,2,n6,n7,n8,(n9+1))*rat(1,2)
      +`Z'(1+n1,n2,n3,1,1,-1+n6,n7,n8,(n9+1))*rat(-n1,2)
      +`Z'(n1,-1+n2,-1+n3,1,1,n6,1+n7,n8,(n9+1))*rat(n7,1)
      +`Z'(n1,-1+n2,n3,1,1,n6,n7,n8,(n9+1))*rat(4*ep+n1+n3+n6+n7+2*n8-4,2)
      +`Z'(n1,n2,-1+n3,1,1,n6,1+n7,n8,(n9+1))*rat(n7,2)
      +`Z'(n1,n2,-1+n3,2,1,n6,n7,n8,(n9+1))*rat(1,2)
      +`Z'(n1,n2,1+n3,1,1,n6,-1+n7,n8,(n9+1))*rat(-n3,2)
      +`Z'(n1,n2,n3,1,1,n6,n7,n8,(n9+1))*rat(n1+n3-n6-n7,2)
      +`Z'(n1,n2,n3,1,2,-1+n6,n7,n8,(n9+1))*rat(-1,2)
      +`Z'(n1,n2,n3,2,1,n6,-1+n7,n8,(n9+1))*rat(-1,2)
      +`Z'(-2+n1,n2,n3,1,1,1+n6,n7,n8,(n9+1))*rat(-n6,2)
      +`Z'(-2+n1,n2,n3,1,2,n6,n7,n8,(n9+1))*rat(-1,2)
      +`Z'(-1+n1,n2,n3,1,1,n6,n7,n8,(n9+1))*rat(-n1+n6+1,2)
      +`Z'(-1+n1,n2,n3,1,2,-1+n6,n7,n8,(n9+1))*rat(1,2)
      +`Z'(n1,-1+n2,n3,1,1,n6,n7,1+n8,(n9+1))*rat(n8,2)
      +`Z'(n1,n2,-2+n3,1,1,n6,1+n7,n8,(n9+1))*rat(-n7,2)
      +`Z'(n1,n2,-2+n3,2,1,n6,n7,n8,(n9+1))*rat(-1,2)
      +`Z'(n1,n2,-1+n3,1,1,n6,n7,n8,(n9+1))*rat(-n3+n7+1,2)
      +`Z'(n1,n2,-1+n3,2,1,n6,-1+n7,n8,(n9+1))*rat(1,2)
      +`Z'(n1,n2,n3,1,1,-1+n6,n7,n8,(n9+1))*rat(n1-1,2)
      +`Z'(n1,n2,n3,1,1,n6,-1+n7,n8,(n9+1))*rat(n3-1,2)
      +`Z'(-1+n1,-1+n2,n3,1,1,n6,n7,1+n8,(n9+1))*rat(-n8,4)
      +`Z'(-1+n1,n2,n3,1,1,n6,n7,1+n8,(n9+1))*rat(n8,4)
      +`Z'(n1,-1+n2,-1+n3,1,1,n6,n7,1+n8,(n9+1))*rat(-n8,4)
      +`Z'(n1,n2,-1+n3,1,1,n6,n7,1+n8,(n9+1))*rat(n8,4)
      +`Z'(n1,n2,n3,1,1,-1+n6,n7,1+n8,(n9+1))*rat(-n8,4)
      +`Z'(n1,n2,n3,1,1,n6,-1+n7,1+n8,(n9+1))*rat(-n8,4)
      +`Z'(-1+n1,n2,-1+n3,1,1,n6,n7,1+n8,(n9+1))*rat(-n8,2)
      +`Z'(-1+n1,n2,n3,1,1,n6,-1+n7,1+n8,(n9+1))*rat(n8,4)
      +`Z'(n1,n2,-1+n3,1,1,-1+n6,n7,1+n8,(n9+1))*rat(n8,4))
      ;
#endif

#if (`il' >= 4)
id,ifmatch->1,`Z'(n1?pos_,n2?,n3?pos_,1,1,n6?pos_,n7?pos_,n8?neg_,0) =
            -rat(1,-4*ep-2*n2-n6-n7+4) * (
      +`Z'(-1+n1,-1+n2,n3,1,1,1+n6,n7,(n8+1),0)*rat(n6,1)
      +`Z'(-1+n1,1+n2,-1+n3,1,1,n6,n7,(n8+1),0)*rat(ep+n2,1)
      +`Z'(-1+n1,1+n2,n3,1,1,n6,-1+n7,(n8+1),0)*rat(-ep-n2,2)
      +`Z'(-1+n1,n2,n3,1,1,1+n6,-1+n7,(n8+1),0)*rat(-n6,2)
      +`Z'(-1+n1,n2,n3,1,1,n6,n7,(n8+1),0)*rat(-n7+1,2)
      +`Z'(n1,-2+n2,n3,1,1,1+n6,n7,(n8+1),0)*rat(-n6,1)
      +`Z'(n1,-2+n2,n3,1,1,n6,1+n7,(n8+1),0)*rat(-n7,1)
      +`Z'(n1,-1+n2,-1+n3,1,1,n6,1+n7,(n8+1),0)*rat(n7,1)
      +`Z'(n1,-1+n2,n3,1,1,-1+n6,1+n7,(n8+1),0)*rat(n7,2)
      +`Z'(n1,-1+n2,n3,1,1,1+n6,-1+n7,(n8+1),0)*rat(n6,2)
      +`Z'(n1,-1+n2,n3,1,1,n6,n7,(n8+1),0)*rat(-2*ep-2*n2+n6+n7+2,2)
      +`Z'(n1,1+n2,-1+n3,1,1,-1+n6,n7,(n8+1),0)*rat(-ep-n2,2)
      +`Z'(n1,n2,-1+n3,1,1,-1+n6,1+n7,(n8+1),0)*rat(-n7,2)
      +`Z'(n1,n2,-1+n3,1,1,n6,n7,(n8+1),0)*rat(-n6+1,2)
      +`Z'(n1,n2,n3,1,1,-1+n6,n7,(n8+1),0)*rat(ep+n2-1,2)
      +`Z'(n1,n2,n3,1,1,n6,-1+n7,(n8+1),0)*rat(ep+n2-1,2))
      ;
#endif

#if (`il' >= 5)
id,ifmatch->1,`Z'(n1?pos_,n2?,n3?pos_,1,1,n6?{>1},n7?pos_,0,0) =
            -rat(1,3*ep+n2+2*(n6-1)-2) * (
      +`Z'(n1,n2,n3,1,1,(n6-1),1+n7,0,0)*rat(-3*ep*n7-n2*n7-2*n7^2+2*n7,(n6-1))
      +`Z'(-1+n1,1+n2,n3,1,1,1+(n6-1),n7,0,0)*rat(-ep-n2,1)
      +`Z'(n1,1+n2,-1+n3,1,1,(n6-1),1+n7,0,0)*rat(ep*n7+n2*n7,(n6-1))
      +`Z'(n1,1+n2,n3,1,1,(n6-1),n7,0,0)*rat(ep*(n6-1)-ep*n7+n2*(n6-1)-n2*n7,(n6-1)))
      ;
#endif

#if (`il' >= 6)
id,ifmatch->1,`Z'(n1?pos_,n2?,n3?pos_,1,1,1,n7?{>1},0,0) =
            -rat(3*ep+n2,-6*ep*(n7-1)-2*n2*(n7-1)-2*(n7-1)^2+2*(n7-1)) * (
      +`Z'(-1+n1,1+n2,n3,1,1,2,(n7-1),0,0)*rat(2*ep,3*ep+n2)
      +`Z'(n1,1+n2,-1+n3,1,1,1,1+(n7-1),0,0)*rat(4*ep*(n7-1)+2*n2*(n7-1),3*ep+n2)
      +`Z'(n1,1+n2,n3,1,1,1,(n7-1),0,0)*rat(-12*ep^2-10*ep*n2-4*ep*(n7-1)+4*ep-2*n2^2-2*n2*(n7-1)+2*n2,3*ep+n2))
      ;
#endif

#if (`il' >= 7)
id,ifmatch->1,`Z'(n1?pos_,n2?,n3?{>1},1,1,1,1,0,0) =
            -rat(1,3*ep*(n3-1)+n2*(n3-1)-(n3-1)) * (
      +`Z'(-2+n1,1+n2,1+(n3-1),1,1,2,1,0,0)*rat(-2*ep*(n3-1)-n2*(n3-1),3*ep+n2)
      +`Z'(-1+n1,1+n2,1+(n3-1),1,1,1,1,0,0)*rat(2*ep*(n3-1)+n2*(n3-1),1)
      +`Z'(-1+n1,1+n2,(n3-1),1,1,1,2,0,0)*rat(-ep*(n3-1),3*ep+n2)
      +`Z'(-1+n1,1+n2,(n3-1),1,1,2,1,0,0)*rat(7*ep^2+3*ep*n2+2*ep*(n3-1)-2*ep+n2*(n3-1)-n2,3*ep+n2)
      +`Z'(-1+n1,n2,(n3-1),1,1,2,1,0,0)*rat(-ep,1)
      +`Z'(n1,1+n2,-1+(n3-1),1,1,1,2,0,0)*rat(5*ep^2+ep*n2+ep*(n3-1)-ep,3*ep+n2)
      +`Z'(n1,1+n2,-1+(n3-1),2,1,1,1,0,0)*rat(-ep-n2,1)
      +`Z'(n1,1+n2,(n3-1),1,1,1,1,0,0)*rat(-10*ep^2-7*ep*n2-3*ep*(n3-1)+3*ep-n2^2-2*n2*(n3-1)+2*n2,1)
      +`Z'(n1,n2,-1+(n3-1),1,1,1,2,0,0)*rat(ep,1)
      +`Z'(n1,n2,-1+(n3-1),2,1,1,1,0,0)*rat(3*ep+n2-1,1)
      +`Z'(n1,n2,(n3-1),1,1,1,1,0,0)*rat(12*ep^2+7*ep*n2+6*ep*(n3-1)-13*ep+n2^2+2*n2*(n3-1)-4*n2-2*(n3-1)+3,1)
      +`Y'(n1,n2,(n3-1),2,0,2,1,0,0)*rat(-1,1)
      +`Y'(-1+n1,1+n2,1+(n3-1),1,1,2,0,0,0)*rat(ep*(n3-1)+n2*(n3-1),3*ep+n2)
      +`Y'(n1,1+n2,1+(n3-1),1,1,1,0,0,0)*rat(3*ep^2*(n3-1)+4*ep*n2*(n3-1)-ep*(n3-1)+n2^2*(n3-1)-n2*(n3-1),3*ep+n2)
      +`Y'(n1,1+n2,(n3-1),2,1,1,0,0,0)*rat(ep+n2,1)
      +`Y'(n1,n2,-1+(n3-1),2,1,2,0,0,0)*rat(1,1)
      +`Y'(n1,n2,(n3-1),1,1,2,0,0,0)*rat(2*ep+2*(n3-1)-3,1))
      ;
#endif

#if (`il' >= 8)
id,ifmatch->1,`Z'(n1?{>1},n2?,1,1,1,1,1,0,0) =
            -rat(1,3*ep*(n1-1)+(n1-1)*n2-(n1-1)) * (
      +`Z'(-1+(n1-1),1+n2,1,1,1,2,1,0,0)*rat(10*ep^2+2*ep*(n1-1)+7*ep*n2-2*ep+(n1-1)*n2+n2^2-n2,3*ep+n2)
      +`Z'(-1+(n1-1),n2,1,1,1,2,1,0,0)*rat(ep,1)
      +`Z'(-1+(n1-1),n2,1,1,2,1,1,0,0)*rat(3*ep+n2-1,1)
      +`Z'((n1-1),1+n2,1,1,1,1,1,0,0)*rat(-10*ep^2-2*ep*(n1-1)-7*ep*n2+2*ep-(n1-1)*n2-n2^2+n2,1)
      +`Z'((n1-1),n2,1,1,1,1,1,0,0)*rat(12*ep^2+6*ep*(n1-1)+7*ep*n2-13*ep+2*(n1-1)*n2-2*(n1-1)+n2^2-4*n2+3,1)
      +`Y'((n1-1),-1+n2,1,0,1,3,1,0,0)*rat(2,1)
      +`Y'((n1-1),-1+n2,1,1,1,3,0,0,0)*rat(-2,1)
      +`Y'(-1+(n1-1),n2,1,0,1,3,1,0,0)*rat(-2,1)
      +`Y'(-1+(n1-1),n2,1,0,2,2,1,0,0)*rat(-1,1)
      +`Y'(-1+(n1-1),n2,1,1,1,3,0,0,0)*rat(2,1)
      +`Y'(-1+(n1-1),n2,1,1,2,2,0,0,0)*rat(1,1)
      +`Y'(1+(n1-1),1+n2,-1,1,1,1,2,0,0)*rat(-ep*(n1-1),3*ep+n2)
      +`Y'(1+(n1-1),1+n2,0,1,1,1,1,0,0)*rat(ep*(n1-1),1)
      +`Y'(1+(n1-1),1+n2,1,1,1,1,0,0,0)*rat(3*ep^2*(n1-1)+4*ep*(n1-1)*n2-ep*(n1-1)+(n1-1)*n2^2-(n1-1)*n2,3*ep+n2)
      +`Y'((n1-1),1+n2,0,1,1,1,2,0,0)*rat(5*ep^2+ep*(n1-1)+ep*n2-ep,3*ep+n2)
      +`Y'((n1-1),1+n2,0,1,1,2,1,0,0)*rat(-3*ep^2-2*ep*(n1-1)-4*ep*n2-(n1-1)*n2-n2^2,3*ep+n2)
      +`Y'((n1-1),1+n2,0,1,2,1,1,0,0)*rat(-ep-n2,1)
      +`Y'((n1-1),1+n2,1,1,1,2,0,0,0)*rat(3*ep^2+ep*(n1-1)+4*ep*n2+(n1-1)*n2+n2^2,3*ep+n2)
      +`Y'((n1-1),1+n2,1,1,2,1,0,0,0)*rat(ep+n2,1)
      +`Y'((n1-1),n2,0,1,1,1,2,0,0)*rat(-ep,1)
      +`Y'((n1-1),n2,1,0,1,2,1,0,0)*rat(-2*ep-2*(n1-1)+1,1)
      +`Y'((n1-1),n2,1,1,1,2,0,0,0)*rat(2*ep+2*(n1-1)-1,1))
      ;

#endif

#if (`il' >= 9)
id,ifmatch->1,`Z'(1,n2?{>1},1,1,1,1,1,0,0) =
		-rat(3*ep+n2-2,-192*ep^3-208*ep^2*n2+232*ep^2-72*ep*n2^2+164*ep*n2-92*ep-8*n2^3+28*n2^2-32*n2+12) * (
      +`Z'(1,-1+n2,1,1,1,1,1,0,0)*rat(128*ep^2+64*ep*n2-112*ep+8*n2^2-28*n2+24,1)
      +`Y'(1,-2+n2,1,0,1,3,1,0,0)*rat(28*ep+8*n2-12,3*ep+n2-2)
      +`Y'(1,-2+n2,1,1,1,3,0,0,0)*rat(-28*ep-8*n2+12,3*ep+n2-2)
      +`Y'(1,-1+n2,1,2,0,2,1,0,0)*rat(-4*ep-2*n2+2,3*ep+n2-2)
      +`Y'(2,-1+n2,1,2,0,1,1,0,0)*rat(2,1)
      +`Y'(-1,n2,2,1,1,2,1,0,0)*rat(-8*ep^2-8*ep*n2+8*ep-2*n2^2+4*n2-2,9*ep^2+6*ep*n2-9*ep+n2^2-3*n2+2)
      +`Y'(0,-1+n2,1,0,1,3,1,0,0)*rat(-28*ep-8*n2+12,3*ep+n2-2)
      +`Y'(0,-1+n2,1,0,2,2,1,0,0)*rat(-14*ep-4*n2+6,3*ep+n2-2)
      +`Y'(0,-1+n2,1,1,1,3,0,0,0)*rat(28*ep+8*n2-12,3*ep+n2-2)
      +`Y'(0,-1+n2,1,1,2,1,1,0,0)*rat(16*ep+4*n2-6,1)
      +`Y'(0,-1+n2,1,1,2,2,0,0,0)*rat(14*ep+4*n2-6,3*ep+n2-2)
      +`Y'(0,n2,1,1,1,1,2,0,0)*rat(-4*ep^2-2*ep*n2+2*ep,9*ep^2+6*ep*n2-9*ep+n2^2-3*n2+2)
      +`Y'(0,n2,1,1,1,2,1,0,0)*rat(168*ep^3+164*ep^2*n2-184*ep^2+48*ep*n2^2-110*ep*n2+62*ep+4*n2^3-14*n2^2+16*n2-6,9*
      ep^2+6*ep*n2-9*ep+n2^2-3*n2+2)
      +`Y'(0,n2,2,1,1,1,1,0,0)*rat(8*ep^2+8*ep*n2-8*ep+2*n2^2-4*n2+2,3*ep+n2-2)
      +`Y'(0,n2,2,1,1,2,0,0,0)*rat(4*ep^2+6*ep*n2-6*ep+2*n2^2-4*n2+2,9*ep^2+6*ep*n2-9*ep+n2^2-3*n2+2)
      +`Y'(1,-1+n2,0,1,1,2,1,0,0)*rat(2*ep,3*ep+n2-2)
      +`Y'(1,-1+n2,0,2,1,1,1,0,0)*rat(16*ep+4*n2-6,1)
      +`Y'(1,-1+n2,0,2,1,2,0,0,0)*rat(4*ep+2*n2-2,3*ep+n2-2)
      +`Y'(1,-1+n2,1,0,1,2,1,0,0)*rat(-28*ep^2-8*ep*n2-2*ep-4*n2+6,3*ep+n2-2)
      +`Y'(1,-1+n2,1,1,1,2,0,0,0)*rat(36*ep^2+12*ep*n2-6*ep+2*n2-4,3*ep+n2-2)
      +`Y'(1,n2,0,1,1,1,2,0,0)*rat(90*ep^3+48*ep^2*n2-58*ep^2+6*ep*n2^2-14*ep*n2+8*ep,9*ep^2+6*ep*n2-9*ep+n2^2-3*n2+2)
      +`Y'(1,n2,0,1,1,2,1,0,0)*rat(-42*ep^3-68*ep^2*n2+46*ep^2-30*ep*n2^2+46*ep*n2-12*ep-4*n2^3+10*n2^2-6*n2,9*ep^2+6*
      ep*n2-9*ep+n2^2-3*n2+2)
      +`Y'(1,n2,0,1,2,1,1,0,0)*rat(-14*ep^2-18*ep*n2+20*ep-4*n2^2+10*n2-6,3*ep+n2-2)
      +`Y'(1,n2,0,2,1,1,1,0,0)*rat(-10*ep^2-14*ep*n2+16*ep-4*n2^2+10*n2-6,3*ep+n2-2)
      +`Y'(1,n2,1,1,1,2,0,0,0)*rat(42*ep^3+68*ep^2*n2-60*ep^2+30*ep*n2^2-50*ep*n2+18*ep+4*n2^3-10*n2^2+6*n2,9*ep^2+6*ep
      *n2-9*ep+n2^2-3*n2+2)
      +`Y'(1,n2,1,1,2,1,0,0,0)*rat(14*ep^2+18*ep*n2-20*ep+4*n2^2-10*n2+6,3*ep+n2-2)
      +`Y'(1,n2,1,2,1,1,0,0,0)*rat(10*ep^2+14*ep*n2-16*ep+4*n2^2-10*n2+6,3*ep+n2-2)
      +`Y'(1,n2,2,1,1,1,0,0,0)*rat(10*ep^2+14*ep*n2-14*ep+4*n2^2-8*n2+4,3*ep+n2-1)
      +`Y'(2,-1+n2,-1,1,1,1,2,0,0)*rat(-2*ep,3*ep+n2-2)
      +`Y'(2,-1+n2,-1,2,1,1,1,0,0)*rat(-2,1)
      +`Y'(2,-1+n2,0,1,1,1,1,0,0)*rat(-8*ep-2*n2+8,1)
      +`Y'(2,n2,-1,1,1,1,2,0,0)*rat(-14*ep^2-4*ep*n2+6*ep,9*ep^2+6*ep*n2-9*ep+n2^2-3*n2+2)
      +`Y'(2,n2,0,1,1,1,1,0,0)*rat(14*ep^2+4*ep*n2-6*ep,3*ep+n2-2)
      +`Y'(2,n2,1,1,1,1,0,0,0)*rat(14*ep^2+18*ep*n2-20*ep+4*n2^2-10*n2+6,3*ep+n2-1))
      ;

id,ifmatch->1,`Z'(1,n2?{<1},1,1,1,1,1,0,0) =
      		-rat(1,128*ep^2+64*ep*n2-48*ep+8*n2^2-12*n2+4) * (      
      +`Z'(1,1+n2,1,1,1,1,1,0,0)*rat(-192*ep^3-208*ep^2*n2+24*ep^2-72*ep*n2^2+20*ep*n2-8*n2^3+4*n2^2,3*ep+n2-1)
      +`Y'(1,-1+n2,1,0,1,3,1,0,0)*rat(28*ep+8*n2-4,3*ep+n2-1)
      +`Y'(1,-1+n2,1,1,1,3,0,0,0)*rat(-28*ep-8*n2+4,3*ep+n2-1)
      +`Y'(1,n2,1,2,0,2,1,0,0)*rat(-4*ep-2*n2,3*ep+n2-1)
      +`Y'(2,n2,1,2,0,1,1,0,0)*rat(2,1)
      +`Y'(-1,1+n2,2,1,1,2,1,0,0)*rat(-8*ep^2-8*ep*n2-2*n2^2,9*ep^2+6*ep*n2-3*ep+n2^2-n2)
      +`Y'(0,1+n2,1,1,1,1,2,0,0)*rat(-4*ep^2-2*ep*n2,9*ep^2+6*ep*n2-3*ep+n2^2-n2)
      +`Y'(0,1+n2,1,1,1,2,1,0,0)*rat(168*ep^3+164*ep^2*n2-20*ep^2+48*ep*n2^2-14*ep*n2+4*n2^3-2*n2^2,9*ep^2+6*ep*n2-3*ep
      +n2^2-n2)
      +`Y'(0,1+n2,2,1,1,1,1,0,0)*rat(8*ep^2+8*ep*n2+2*n2^2,3*ep+n2-1)
      +`Y'(0,1+n2,2,1,1,2,0,0,0)*rat(4*ep^2+6*ep*n2+2*n2^2,9*ep^2+6*ep*n2-3*ep+n2^2-n2)
      +`Y'(0,n2,1,0,1,3,1,0,0)*rat(-28*ep-8*n2+4,3*ep+n2-1)
      +`Y'(0,n2,1,0,2,2,1,0,0)*rat(-14*ep-4*n2+2,3*ep+n2-1)
      +`Y'(0,n2,1,1,1,3,0,0,0)*rat(28*ep+8*n2-4,3*ep+n2-1)
      +`Y'(0,n2,1,1,2,1,1,0,0)*rat(16*ep+4*n2-2,1)
      +`Y'(0,n2,1,1,2,2,0,0,0)*rat(14*ep+4*n2-2,3*ep+n2-1)
      +`Y'(1,1+n2,0,1,1,1,2,0,0)*rat(90*ep^3+48*ep^2*n2-10*ep^2+6*ep*n2^2-2*ep*n2,9*ep^2+6*ep*n2-3*ep+n2^2-n2)
      +`Y'(1,1+n2,0,1,1,2,1,0,0)*rat(-42*ep^3-68*ep^2*n2-22*ep^2-30*ep*n2^2-14*ep*n2+4*ep-4*n2^3-2*n2^2+2*n2,9*ep^2+6*
      ep*n2-3*ep+n2^2-n2)
      +`Y'(1,1+n2,0,1,2,1,1,0,0)*rat(-14*ep^2-18*ep*n2+2*ep-4*n2^2+2*n2,3*ep+n2-1)
      +`Y'(1,1+n2,0,2,1,1,1,0,0)*rat(-10*ep^2-14*ep*n2+2*ep-4*n2^2+2*n2,3*ep+n2-1)
      +`Y'(1,1+n2,1,1,1,2,0,0,0)*rat(42*ep^3+68*ep^2*n2+8*ep^2+30*ep*n2^2+10*ep*n2-2*ep+4*n2^3+2*n2^2-2*n2,9*ep^2+6*ep*
      n2-3*ep+n2^2-n2)
      +`Y'(1,1+n2,1,1,2,1,0,0,0)*rat(14*ep^2+18*ep*n2-2*ep+4*n2^2-2*n2,3*ep+n2-1)
      +`Y'(1,1+n2,1,2,1,1,0,0,0)*rat(10*ep^2+14*ep*n2-2*ep+4*n2^2-2*n2,3*ep+n2-1)
      +`Y'(1,1+n2,2,1,1,1,0,0,0)*rat(10*ep^2+14*ep*n2+4*n2^2,3*ep+n2)
      +`Y'(1,n2,0,1,1,2,1,0,0)*rat(2*ep,3*ep+n2-1)
      +`Y'(1,n2,0,2,1,1,1,0,0)*rat(16*ep+4*n2-2,1)
      +`Y'(1,n2,0,2,1,2,0,0,0)*rat(4*ep+2*n2,3*ep+n2-1)
      +`Y'(1,n2,1,0,1,2,1,0,0)*rat(-28*ep^2-8*ep*n2-10*ep-4*n2+2,3*ep+n2-1)
      +`Y'(1,n2,1,1,1,2,0,0,0)*rat(36*ep^2+12*ep*n2+6*ep+2*n2-2,3*ep+n2-1)
      +`Y'(2,1+n2,-1,1,1,1,2,0,0)*rat(-14*ep^2-4*ep*n2+2*ep,9*ep^2+6*ep*n2-3*ep+n2^2-n2)
      +`Y'(2,1+n2,0,1,1,1,1,0,0)*rat(14*ep^2+4*ep*n2-2*ep,3*ep+n2-1)
      +`Y'(2,1+n2,1,1,1,1,0,0,0)*rat(14*ep^2+18*ep*n2-2*ep+4*n2^2-2*n2,3*ep+n2)
      +`Y'(2,n2,-1,1,1,1,2,0,0)*rat(-2*ep,3*ep+n2-1)
      +`Y'(2,n2,-1,2,1,1,1,0,0)*rat(-2,1)
      +`Y'(2,n2,0,1,1,1,1,0,0)*rat(-8*ep-2*n2+6,1))
      ;
#endif

GoTo 2;
Label 1;
$ired = 1;
Label 2;
id    `Z'(n1?neg0_,?a,n8?,n9?) = `Y'(n1,?a,n8,n9);
id    `Z'(n1?,n2?,?a,x?neg0_,?b,n8?,n9?) = `Y'(n1,n2,?a,x,?b,n8,n9);
ModuleOption,maximum,$ired;
.sort:Main fastar2 loop, pass `$pass';
#if ( `$ired' > 0 )
      #redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
id	`Z'(1,1,1,1,1,1,1,0,0) = Master(fastar2);

#break
*--#] fastar2 : 
*--#[ fastar3 :
#case fastar3
* note: not a master, assumes the n1 is fixed, not the n3!
#do il = 1,9
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering fastar3 loop `il', pass `$pass' at `time_' sec."
#endif

#if (`il' >= 1)
*
*  Prework to make the reduction of the fa topology faster.
*
* B7
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,2*n6-2)*(`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,1+n8,n9)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-n6+1,1)+`Z'(
   n1,1+n2,-1+n3,n4,n5,-1+n6,n7,1+n8,n9)*rat(-n2,1)+`Z'(n1,1+n2,n3,n4,n5,-1+n6,-1+n7,1+n8,n9)*rat(n2,1)+`Z'(n1,n2,-1+n3,n4,
   n5,-1+n6,1+n7,1+n8,n9)*rat(n7,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+n8,
   n9)*rat(-n6+1,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,1+n8,1+n9)*rat(n9,2)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(
   n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9)*rat(-n2+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(-1+n6,1));
* B10
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,2*n4-2)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,1+n9)*rat(-n3,1)+
   `Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,1+n9)*rat(n7,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,-1+n2,n3,
   n4,n5,n6,n7,n8,1+n9)*rat(-n4+1,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,1+n7,n8,1+n9)*rat(-n7,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+
   n9)*rat(-n4+1,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,-1+n7,n8,1+n9)*rat(n3,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,1+n8,1+n9)*rat(n8,2)+
   `Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,1+n9)*rat(-n3+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*rat(-1+n4,1));
* B2
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,2*n5-2)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n5+1,1)+`Z'(1+n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9)*rat(-ep-n1,1)+
   `Z'(1+n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,1+n9)*rat(ep+n1,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,1+n9)*rat(-n7,1)+`Z'(n1,-1+n2,
   n3,n4,-1+n5,n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n5+1,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,1+
   n7,n8,1+n9)*rat(n7,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,1+n9)*rat(
   n8,2)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9)*rat(-3*ep-n1-2*n2-n7-n9+3,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,1+n9)*rat(-1+n5,1)
);
* B5
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?neg_,n9?neg0_) =
  -rat(1,2*n7-2)*(`Z'(-1+n1,1+n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(-n2,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,1+n9)*rat(-n9,2)+
   `Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-n7+1,1)+`Z'(n1,1+n2,n3,n4,n5,-1+n6,-1+n7,1+n8,n9)*rat(n2,1)+`Z'(n1,n2,-1+n3,n4,
   n5,n6,n7,1+n8,n9)*rat(-n7+1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,1+n9)*rat(n9,2)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9)*
   rat(-1+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9)*rat(2*ep+n2+2*n6+n7+
   n9-5,1));
* B12
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?neg0_,n9?neg_) =
  -rat(1,2*n7-2)*(`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n7+1,1)+`Z'(
   n1,n2,1+n3,-1+n4,n5,n6,-1+n7,n8,1+n9)*rat(n3,1)+`Z'(n1,n2,1+n3,n4,n5,n6,-1+n7,n8,1+n9)*rat(-n3,1)+`Z'(n1,n2,n3,-1+n4,n5,
   n6,n7,n8,1+n9)*rat(-1+n7,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,1+n9)*rat(n8,2)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,1+n9)*
   rat(-n8,2)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*rat(2*ep+n3+2*n4+n7+n8-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n7+1,
   1));
* B3
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,-2*n5+2)*(`Z'(1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9)*rat(-2*ep-2*n1,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,1+n8,n9)*rat(-n7,
   1)+`Z'(n1,1+n2,-1+n3,n4,-1+n5,n6,n7,1+n8,n9)*rat(-n2,1)+`Z'(n1,1+n2,n3,n4,-1+n5,n6,-1+n7,1+n8,n9)*rat(n2,1)+`Z'(n1,n2,-1+
   n3,n4,-1+n5,n6,1+n7,1+n8,n9)*rat(n7,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,
   n7,2+n8,n9)*rat(-n8-1,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n5,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,n7,1+n8,1+n9)*rat(
   n9,2)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9)*rat(-n5+1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,n2,n3,n4,-
   1+n5,n6,n7,1+n8,n9)*rat(-n2+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n5,1));
* B4
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,2*n6-2)*(`Z'(-1+n1,n2,n3,n4,1+n5,-1+n6,n7,n8,1+n9)*rat(-n5,1)+`Z'(1+n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,1+n9)*rat(ep+n1,
   1)+`Z'(1+n1,n2,n3,n4,n5,-1+n6,n7,n8,1+n9)*rat(-ep-n1,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,1+n9)*rat(-n7,1)+`Z'(n1,-1+n2,
   n3,n4,n5,n6,n7,n8,1+n9)*rat(-n6+1,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,1+n7,n8,1+n9)*rat(n7,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,
   1+n8,1+n9)*rat(-n8,2)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,1+n8,1+n9)*rat(n8,2)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9)*rat(-1+n6,
   1)+`Z'(n1,n2,n3,n4,1+n5,-1+n6,n7,n8,1+n9)*rat(n5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,1+n9)*rat(-n8,2)+`Z'(n1,n2,n3,n4,n5,-
   1+n6,n7,n8,1+n9)*rat(-3*ep-n1-2*n2+n5-n6-n7-n9+4,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-n6+1,1));
* B8
id,ifmatch->1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
  -rat(1,-2*n2+2)*(`Z'(-1+n1,-1+n2,n3,n4,1+n5,n6,n7,n8,1+n9)*rat(-n5,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,1+n8,1+n9)*rat(-n8,2)
   +`Z'(1+n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9)*rat(ep+n1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*rat(-ep-n1,1)+`Z'(n1,-2+n2,
   n3,n4,n5,n6,1+n7,n8,1+n9)*rat(-n7,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,1+n7,n8,1+n9)*rat(n7,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,
   1+n8,1+n9)*rat(-n8,2)+`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,1+n9)*rat(n3,1)+`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,n8,1+n9)*rat(-n3,
   1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(n8,2)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,1+n8,1+n9)*rat(n8,2)+`Z'(n1,-1+n2,n3,
   n4,1+n5,n6,n7,n8,1+n9)*rat(n5,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,1+n8,1+n9)*rat(-n8,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*
   rat(ep-n1+n3+2*n4+n5+n7+n8-4,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,2+n9)*rat(-n9-1,1));
* B9
id,ifmatch->1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
  -rat(1,-2*n4+2)*(`Z'(-1+n1,1+n2,n3,-1+n4,n5,n6,n7,1+n8,n9)*rat(-n2,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(-n9,2)
   +`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,2+n8,n9)*rat(-n8-1,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n4,1)+`Z'(n1,1+n2,n3,-1+
   n4,n5,-1+n6,n7,1+n8,n9)*rat(n2,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+n8,n9)*rat(-1+n4,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,n9)*
   rat(-2*n3,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,n7,1+n8,1+n9)*rat(n9,2)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9)*rat(-n9,2)+`Z'(n1,
   n2,n3,-1+n4,n5,n6,n7,1+n8,n9)*rat(4*ep+n2+2*n3+n4+2*n6+2*n7+n8+n9-8,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9)*rat(-n4+1,
   1));

#endif

#if (`il' >= 1)
id,ifmatch->1,`Z'(n1?,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
            -rat(1,(n4-1)) * (
      +`Z'(n1,-1+n2,n3,(n4-1),n5,n6,1+n7,n8,n9)*rat(n7,1)
      +`Z'(n1,n2,-1+n3,1+(n4-1),n5,n6,n7,n8,n9)*rat(-(n4-1),1)
      +`Z'(n1,n2,-1+n3,(n4-1),n5,n6,1+n7,n8,n9)*rat(-n7,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9)*rat(-2*ep-2*n3-(n4-1)-n7-n8+4,1))
      ;
#endif

#if (`il' >= 2)
id,ifmatch->1,`Z'(n1?,n2?pos_,n3?pos_,1,n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
            -rat(1,-2*ep-2*n3-n7-n8+3) * (
      +`Z'(n1,-1+n2,1+n3,1,n5,n6,n7,n8,(n9+1))*rat(n3,1)
      +`Z'(n1,-1+n2,n3,2,n5,n6,n7,n8,(n9+1))*rat(1,2)
      +`Z'(n1,n2,-1+n3,2,n5,n6,n7,n8,(n9+1))*rat(1,2)
      +`Z'(n1,n2,1+n3,1,n5,n6,-1+n7,n8,(n9+1))*rat(-n3,2)
      +`Z'(n1,n2,n3,1,n5,n6,n7,n8,(n9+1))*rat(-n7+1,2)
      +`Z'(n1,n2,n3,2,n5,n6,-1+n7,n8,(n9+1))*rat(-1,2)
      +`Z'(-1+n1,n2,n3,1,n5,n6,n7,1+n8,(n9+1))*rat(n8,4)
      +`Z'(n1,-2+n2,n3,1,n5,n6,1+n7,n8,(n9+1))*rat(n7,2)
      +`Z'(n1,-1+n2,-1+n3,2,n5,n6,n7,n8,(n9+1))*rat(-1,2)
      +`Z'(n1,-1+n2,n3,1,n5,n6,n7,n8,(n9+1))*rat(-2*ep-2*n3-n7-n8+3,2)
      +`Z'(n1,n2,-2+n3,1,n5,n6,1+n7,n8,(n9+1))*rat(-n7,2)
      +`Z'(n1,n2,-2+n3,2,n5,n6,n7,n8,(n9+1))*rat(-1,2)
      +`Z'(n1,n2,-1+n3,1,n5,n6,n7,n8,(n9+1))*rat(2*ep+2*n7+n8-2,2)
      +`Z'(n1,n2,-1+n3,2,n5,n6,-1+n7,n8,(n9+1))*rat(1,2)
      +`Z'(n1,n2,n3,1,n5,n6,-1+n7,n8,(n9+1))*rat(n3-1,2)
      +`Z'(-1+n1,-1+n2,n3,1,n5,n6,n7,1+n8,(n9+1))*rat(n8,4)
      +`Z'(-1+n1,n2,-1+n3,1,n5,n6,n7,1+n8,(n9+1))*rat(-n8,2)
      +`Z'(n1,-1+n2,n3,1,n5,n6,n7,1+n8,(n9+1))*rat(n8,2)
      +`Z'(n1,n2,-1+n3,1,n5,n6,n7,1+n8,(n9+1))*rat(-n8,4)
      +`Z'(n1,n2,n3,1,n5,-1+n6,n7,1+n8,(n9+1))*rat(-n8,4)
      +`Z'(n1,-1+n2,-1+n3,1,n5,n6,n7,1+n8,(n9+1))*rat(-n8,4)
      +`Z'(n1,-1+n2,n3,1,-1+n5,n6,n7,1+n8,(n9+1))*rat(-n8,4)
      +`Z'(n1,n2,-1+n3,1,-1+n5,n6,n7,1+n8,(n9+1))*rat(n8,4)
      +`Z'(n1,n2,-1+n3,1,n5,-1+n6,n7,1+n8,(n9+1))*rat(n8,4)
      +`Y'(n1,-1+n2,1+n3,0,n5,n6,n7,n8,(n9+1))*rat(-n3,2)
      +`Y'(n1,-1+n2,n3,0,n5,n6,1+n7,n8,(n9+1))*rat(-n7,2)
      +`Y'(n1,n2,-1+n3,0,n5,n6,1+n7,n8,(n9+1))*rat(n7,2)
      +`Y'(n1,n2,n3,0,n5,n6,n7,n8,(n9+1))*rat(n3-1,2))
      ;
#endif

#if (`il' >= 3)
id,ifmatch->1,`Z'(n1?,n2?pos_,n3?pos_,1,n5?pos_,n6?pos_,n7?pos_,n8?neg_,0) =
            -rat(1,-2*ep-2*n2-n6-n7+4) * (
      +`Z'(-1+n1,-1+n2,n3,1,n5,1+n6,n7,(n8+1),0)*rat(n6,1)
      +`Z'(-1+n1,1+n2,-1+n3,1,n5,n6,n7,(n8+1),0)*rat(n2,1)
      +`Z'(-1+n1,1+n2,n3,1,n5,n6,-1+n7,(n8+1),0)*rat(-n2,2)
      +`Z'(-1+n1,n2,n3,1,n5,1+n6,-1+n7,(n8+1),0)*rat(-n6,2)
      +`Z'(-1+n1,n2,n3,1,n5,n6,n7,(n8+1),0)*rat(-n7+1,2)
      +`Z'(n1,-2+n2,n3,1,n5,1+n6,n7,(n8+1),0)*rat(-n6,1)
      +`Z'(n1,-2+n2,n3,1,n5,n6,1+n7,(n8+1),0)*rat(-n7,1)
      +`Z'(n1,-1+n2,-1+n3,1,n5,n6,1+n7,(n8+1),0)*rat(n7,1)
      +`Z'(n1,-1+n2,n3,1,n5,-1+n6,1+n7,(n8+1),0)*rat(n7,2)
      +`Z'(n1,-1+n2,n3,1,n5,1+n6,-1+n7,(n8+1),0)*rat(n6,2)
      +`Z'(n1,-1+n2,n3,1,n5,n6,n7,(n8+1),0)*rat(-2*n2+n6+n7+2,2)
      +`Z'(n1,1+n2,-1+n3,1,n5,-1+n6,n7,(n8+1),0)*rat(-n2,2)
      +`Z'(n1,n2,-1+n3,1,n5,-1+n6,1+n7,(n8+1),0)*rat(-n7,2)
      +`Z'(n1,n2,-1+n3,1,n5,n6,n7,(n8+1),0)*rat(-n6+1,2)
      +`Z'(n1,n2,n3,1,n5,-1+n6,n7,(n8+1),0)*rat(n2-1,2)
      +`Z'(n1,n2,n3,1,n5,n6,-1+n7,(n8+1),0)*rat(n2-1,2))
      ;
#endif

#if (`il' >= 4)
id,ifmatch->1,`Z'(n1?,n2?pos_,n3?pos_,1,n5?pos_,n6?{>1},n7?pos_,0,0) =
            -rat(1,2*ep+2*(n6-1)-2) * (
      +`Z'(n1,1+n2,-1+n3,1,n5,(n6-1),1+n7,0,0)*rat(2*n2*n7,(n6-1))
      +`Z'(n1,1+n2,n3,1,n5,(n6-1),n7,0,0)*rat(-2*ep*n2-2*n2^2-2*n2*n7+2*n2,(n6-1))
      +`Z'(n1,n2,n3,1,n5,(n6-1),1+n7,0,0)*rat(-2*ep*n7-2*n2*n7-2*n7^2+2*n7,(n6-1)))
      ;
#endif

#if (`il' >= 5)
id,ifmatch->1,`Z'(n1?,n2?pos_,n3?{>1},1,n5?pos_,1,n7?pos_,0,0) =
            -rat(1,-2*ep*(n3-1)-2*(n3-1)^2-2*(n3-1)*n7+2*(n3-1)) * (
      +`Z'(n1,-1+n2,1+(n3-1),1,n5,1,1+n7,0,0)*rat(2*(n3-1)*n7,1)
      +`Z'(n1,n2,(n3-1),1,n5,1,1+n7,0,0)*rat(-2*ep*n7-2*(n3-1)*n7-2*n7^2+2*n7,1)
      +`Z'(n1,-1+n2,(n3-1),1,n5,1,1+n7,0,0)*rat(-2*ep*n7,1)
      +`Z'(n1,n2,-1+(n3-1),1,n5,1,1+n7,0,0)*rat(2*ep*n7,1)
      +`Z'(n1,n2,-1+(n3-1),2,n5,1,n7,0,0)*rat(2*ep,1)
      +`Z'(n1,n2,(n3-1),1,n5,1,n7,0,0)*rat(4*ep^2+4*ep*(n3-1)+2*ep*n7-6*ep,1))
      ;
#endif

#if (`il' >= 6)
id,ifmatch->1,`Z'(n1?,n2?{>1},1,1,n5?pos_,1,n7?pos_,0,0) =
            -rat(ep,2*ep^2*(n2-1)+2*ep*(n2-1)^2+3*ep*(n2-1)*n7-3*ep*(n2-1)+(n2-1)^2*n7-(n2-1)^2+(n2-1)*n7^2-2*(n2-1)*n7+(n2-1)) * (
      +`Z'(n1,-1+(n2-1),2,1,n5,1,1+n7,0,0)*rat(-2*ep^2*n7-3*ep*(n2-1)*n7-3*ep*n7^2+4*ep*n7-(n2-1)^2*n7-2*(n2-1)*n7^2+3*(n2-1)*n7-n7^3+3*
      n7^2-2*n7,ep^2+2*ep*n7-ep+n7^2-n7)
      +`Z'(n1,(n2-1),1,1,n5,1,1+n7,0,0)*rat(4*ep^3*n7+5*ep^2*(n2-1)*n7+8*ep^2*n7^2-9*ep^2*n7+ep*(n2-1)^2*n7+5*ep*(n2-1)*n7^2-6*ep*(n2-1)*
      n7+5*ep*n7^3-11*ep*n7^2+6*ep*n7+(n2-1)*n7^3-2*(n2-1)*n7^2+(n2-1)*n7+n7^4-3*n7^3+3*n7^2-n7,ep^2+ep*n7-ep)
      +`Z'(n1,-1+(n2-1),1,1,-1+n5,1,2+n7,0,0)*rat(ep*n7^2+ep*n7+(n2-1)*n7^2+(n2-1)*n7+n7^3-n7,ep)
      +`Z'(n1,-1+(n2-1),1,1,n5,1,1+n7,0,0)*rat(-2*ep^4*n7+ep^3*(n2-1)*n7-9*ep^3*n7^2+6*ep^3*n7+ep^2*(n2-1)^2*n7-3*ep^2*(n2-1)*n7^2-12*
      ep^2*n7^3+20*ep^2*n7^2-6*ep^2*n7-4*ep*(n2-1)*n7^3+5*ep*(n2-1)*n7^2-ep*(n2-1)*n7-6*ep*n7^4+17*ep*n7^3-13*ep*n7^2+2*ep*n7-(n2-1)*
      n7^4+2*(n2-1)*n7^3-(n2-1)*n7^2-n7^5+4*n7^4-5*n7^3+2*n7^2,ep^3+2*ep^2*n7-ep^2+ep*n7^2-ep*n7)
      +`Z'(n1,1+(n2-1),1,1,-1+n5,1,n7,0,0)*rat(-2*ep^2*(n2-1)-2*ep*(n2-1)^2-3*ep*(n2-1)*n7+3*ep*(n2-1)-(n2-1)^2*n7+(n2-1)^2-(n2-1)*n7^2+2*(n2-1)*n7-(n2-1),ep)
      +`Z'(n1,1+(n2-1),1,1,n5,1,-1+n7,0,0)*rat(8*ep^3*(n2-1)+4*ep^2*(n2-1)^2+12*ep^2*(n2-1)*n7-20*ep^2*(n2-1)+4*ep*(n2-1)^2*n7-6*ep*(n2-1)^2+6*ep*
      (n2-1)*n7^2-20*ep*(n2-1)*n7+16*ep*(n2-1)+(n2-1)^2*n7^2-3*(n2-1)^2*n7+2*(n2-1)^2+(n2-1)*n7^3-5*(n2-1)*n7^2+8*(n2-1)*n7-4*(n2-1),ep^2+ep*n7-ep)
      +`Z'(n1,(n2-1),1,1,-1+n5,1,1+n7,0,0)*rat(-2*ep^2*n7-ep*(n2-1)*n7-3*ep*n7^2+2*ep*n7+(n2-1)^2*n7-(n2-1)*n7-n7^3+n7^2,ep)
      +`Z'(n1,(n2-1),1,1,n5,1,n7,0,0)*rat(4*ep^5-2*ep^4*(n2-1)+20*ep^4*n7-14*ep^4-2*ep^3*(n2-1)^2+ep^3*(n2-1)*n7+5*ep^3*(n2-1)+33*ep^3*n7^2
      -55*ep^3*n7+18*ep^3-3*ep^2*(n2-1)^2*n7+3*ep^2*(n2-1)^2+3*ep^2*(n2-1)*n7^2+ep^2*(n2-1)*n7-4*ep^2*(n2-1)+24*ep^2*n7^3-66*ep^2*n7^2+52*
      ep^2*n7-10*ep^2-3*ep*(n2-1)^2*n7^2+4*ep*(n2-1)^2*n7-ep*(n2-1)^2+ep*(n2-1)*n7^3+ep*(n2-1)*n7^2-3*ep*(n2-1)*n7+ep*(n2-1)+8*ep*n7^4-31*ep*n7^3+
      40*ep*n7^2-19*ep*n7+2*ep-(n2-1)^2*n7^3+2*(n2-1)^2*n7^2-(n2-1)^2*n7+(n2-1)*n7^3-2*(n2-1)*n7^2+(n2-1)*n7+n7^5-5*n7^4+9*n7^3-7*n7^2+2*n7,
      ep^3+2*ep^2*n7-ep^2+ep*n7^2-ep*n7)
      +`Y'(n1,1+(n2-1),0,1,n5,1,1+n7,0,0)*rat(-3*ep*(n2-1)*n7-2*(n2-1)*n7^2+2*(n2-1)*n7,ep)
      +`Y'(n1,1+(n2-1),-1,1,-1+n5,1,2+n7,0,0)*rat((n2-1)*n7^2+(n2-1)*n7,ep)
      +`Y'(n1,1+(n2-1),-1,1,n5,1,1+n7,0,0)*rat(-ep*(n2-1)*n7-(n2-1)*n7^2+(n2-1)*n7,ep)
      +`Y'(n1,1+(n2-1),-1,2,n5,1,n7,0,0)*rat(-(n2-1),1)
      +`Y'(n1,1+(n2-1),0,1,-1+n5,1,1+n7,0,0)*rat(ep*(n2-1)*n7-(n2-1)^2*n7-(n2-1)*n7,ep)
      +`Y'(n1,1+(n2-1),0,1,n5,1,n7,0,0)*rat(-2*ep^3*(n2-1)-ep^2*(n2-1)*n7+3*ep^2*(n2-1)+2*ep*(n2-1)^2*n7-2*ep*(n2-1)^2+2*ep*(n2-1)*n7-2*ep*(n2-1)+(n2-1)^2*
      n7^2-2*(n2-1)^2*n7+(n2-1)^2+(n2-1)*n7^2-2*(n2-1)*n7+(n2-1),ep^2+ep*n7-ep)
      +`Y'(n1,1+(n2-1),0,2,n5,1,-1+n7,0,0)*rat(3*ep*(n2-1)+(n2-1)^2+2*(n2-1)*n7-3*(n2-1),ep+n7-1)
      +`Y'(n1,(n2-1),0,1,-1+n5,1,2+n7,0,0)*rat(-ep*n7^2-ep*n7-2*(n2-1)*n7^2-2*(n2-1)*n7-n7^3+n7,ep)
      +`Y'(n1,(n2-1),0,1,n5,1,1+n7,0,0)*rat(2*ep^4*n7+9*ep^3*n7^2-6*ep^3*n7-ep^2*(n2-1)^2*n7+6*ep^2*(n2-1)*n7^2-2*ep^2*(n2-1)*n7+12*
      ep^2*n7^3-20*ep^2*n7^2+6*ep^2*n7+7*ep*(n2-1)*n7^3-9*ep*(n2-1)*n7^2+2*ep*(n2-1)*n7+6*ep*n7^4-17*ep*n7^3+13*ep*n7^2-2*ep*n7+2*
      (n2-1)*n7^4-4*(n2-1)*n7^3+2*(n2-1)*n7^2+n7^5-4*n7^4+5*n7^3-2*n7^2,ep^3+2*ep^2*n7-ep^2+ep*n7^2-ep*n7)
      +`Y'(n1,(n2-1),0,2,-1+n5,2,n7,0,0)*rat(-1,1)
      +`Y'(n1,(n2-1),0,2,n5,1,n7,0,0)*rat(2*ep^3-ep^2*(n2-1)+8*ep^2*n7-5*ep^2-ep*(n2-1)^2+2*ep*(n2-1)*n7+ep*(n2-1)+9*ep*n7^2-14*ep*n7+3*ep+
      2*(n2-1)*n7^2-2*(n2-1)*n7+3*n7^3-8*n7^2+5*n7,ep^2+2*ep*n7-ep+n7^2-n7)
      +`Y'(n1,(n2-1),0,2,n5,2,-1+n7,0,0)*rat(1,1))
      ;
#endif

#if (`il' >= 7)
id,ifmatch->1,`Z'(n1?,1,1,1,n5?{>1},1,n7?pos_,0,0) =
            -rat(1,-2*ep*(n5-1)-(n5-1)*n7+(n5-1)) * (
      +`Z'(n1,1,1,1,1+(n5-1),1,-1+n7,0,0)*rat(-4*ep^2*(n5-1)-4*ep*(n5-1)*n7+6*ep*(n5-1)-(n5-1)*n7^2+3*(n5-1)*n7-2*(n5-1),ep+n7-1)
      +`Z'(n1,1,1,1,(n5-1),1,n7,0,0)*rat(6*ep^2+2*ep*n1+2*ep*(n5-1)+3*ep*n7-7*ep+n1*n7-n1+(n5-1)*n7-(n5-1)-2*n7+2,1)
      +`Y'(n1,0,2,1,1+(n5-1),1,n7,0,0)*rat(ep*(n5-1),ep+n7-1)
      +`Y'(1+n1,0,1,2,(n5-1),1,n7,0,0)*rat(-ep-n1,1)
      +`Y'(1+n1,1,1,2,(n5-1),0,n7,0,0)*rat(ep+n1,1)
      +`Y'(n1,0,1,1,1+(n5-1),0,1+n7,0,0)*rat(-(n5-1)*n7,1)
      +`Y'(n1,0,1,1,1+(n5-1),1,n7,0,0)*rat(2*ep*(n5-1)*n7-2*ep*(n5-1)+(n5-1)*n7^2-2*(n5-1)*n7+(n5-1),ep+n7-1)
      +`Y'(n1,0,1,1,(n5-1),1,1+n7,0,0)*rat(-3*ep*n7-n1*n7-(n5-1)*n7+2*n7,1)
      +`Y'(n1,1,0,1,1+(n5-1),0,1+n7,0,0)*rat((n5-1)*n7,1)
      +`Y'(n1,1,0,1,1+(n5-1),1,n7,0,0)*rat(-2*ep*(n5-1)*n7+2*ep*(n5-1)-(n5-1)*n7^2+2*(n5-1)*n7-(n5-1),ep+n7-1)
      +`Y'(n1,1,0,1,(n5-1),1,1+n7,0,0)*rat(3*ep*n7+n1*n7+(n5-1)*n7-2*n7,1)
      +`Y'(n1,1,0,2,1+(n5-1),0,n7,0,0)*rat((n5-1),1)
      +`Y'(n1,1,0,2,1+(n5-1),1,-1+n7,0,0)*rat(-2*ep*(n5-1)-(n5-1)*n7+(n5-1),ep+n7-1)
      +`Y'(n1,1,0,2,(n5-1),1,n7,0,0)*rat(3*ep+n1+(n5-1)-2,1)
      +`Y'(n1,1,1,1,1+(n5-1),0,n7,0,0)*rat(2*ep*(n5-1)+(n5-1)*n7-(n5-1),1))
      ;
#endif

#if (`il' >= 8)
id,ifmatch->1,`Z'(n1?,1,1,1,1,1,n7?{>1},0,0) =
            -rat(ep+(n7-1)-1,2*ep^2*(n7-1)+3*ep*(n7-1)^2-ep*(n7-1)+(n7-1)^3-(n7-1)^2) * (
      +`Z'(1+n1,1,1,1,1,1,(n7-1),0,0)*rat(-2*ep^2-2*ep*n1-ep*(n7-1)+ep-n1*(n7-1)+n1,1)
      +`Z'(n1,1,1,1,1,1,(n7-1),0,0)*rat(10*ep^3+2*ep^2*n1+21*ep^2*(n7-1)-23*ep^2+3*ep*n1*(n7-1)-3*ep*n1+12*ep*(n7-1)^2-29*ep*(n7-1)+17*ep+
      n1*(n7-1)^2-2*n1*(n7-1)+n1+2*(n7-1)^3-8*(n7-1)^2+10*(n7-1)-4,ep+(n7-1)-1)
      +`Z'(n1,2,1,1,1,1,-1+(n7-1),0,0)*rat(4*ep^2+4*ep*(n7-1)-6*ep+(n7-1)^2-3*(n7-1)+2,ep+(n7-1)-1)
      +`Y'(n1,0,2,1,1,1,1+(n7-1),0,0)*rat(-ep*(n7-1),ep+(n7-1)-1)
      +`Y'(n1,2,0,1,1,1,1+(n7-1),0,0)*rat(-(n7-1),1)
      +`Y'(1+n1,0,1,1,1,1,1+(n7-1),0,0)*rat(ep*(n7-1)+n1*(n7-1),1)
      +`Y'(1+n1,1,0,1,1,1,1+(n7-1),0,0)*rat(-ep*(n7-1)-n1*(n7-1),1)
      +`Y'(1+n1,1,0,2,1,1,(n7-1),0,0)*rat(-ep-n1,1)
      +`Y'(1+n1,1,1,2,0,1,(n7-1),0,0)*rat(ep+n1,1)
      +`Y'(n1,0,1,1,1,1,1+(n7-1),0,0)*rat(-5*ep^2*(n7-1)-ep*n1*(n7-1)-8*ep*(n7-1)^2+9*ep*(n7-1)-n1*(n7-1)^2+n1*(n7-1)-2*(n7-1)^3+5*(n7-1)^2-3*(n7-1),ep+(n7-1)-1)
      +`Y'(n1,1,0,1,1,1,1+(n7-1),0,0)*rat(5*ep^2*(n7-1)+ep*n1*(n7-1)+8*ep*(n7-1)^2-8*ep*(n7-1)+n1*(n7-1)^2-n1*(n7-1)+2*(n7-1)^3-4*(n7-1)^2+2*(n7-1),ep+(n7-1)-1)
      +`Y'(n1,1,0,2,1,1,(n7-1),0,0)*rat(5*ep^2+ep*n1+8*ep*(n7-1)-9*ep+n1*(n7-1)-n1+2*(n7-1)^2-5*(n7-1)+3,ep+(n7-1)-1)
      +`Y'(n1,2,-1,1,1,1,1+(n7-1),0,0)*rat(-(n7-1),1)
      +`Y'(n1,2,-1,2,1,1,(n7-1),0,0)*rat(-1,1)
      +`Y'(n1,2,0,1,1,1,(n7-1),0,0)*rat(-2*ep^2-ep*(n7-1)+3*ep+2*(n7-1)-2,ep+(n7-1)-1)
      +`Y'(n1,2,0,2,1,1,-1+(n7-1),0,0)*rat(2*ep+(n7-1)-1,ep+(n7-1)-1))
      ;
#endif

#if (`il' >= 9)
id,ifmatch->1,`Z'(n1?,1,1,1,1,1,1,0,0) =
        -rat(1,-6*ep-2*n1+2) * (
      +`Y'(-1+n1,0,2,1,1,1,2,0,0)*rat(1,2*ep)
      +`Y'(-1+n1,0,2,1,2,1,1,0,0)*rat(-1,2*ep)
      +`Y'(n1,0,2,1,2,1,1,0,0)*rat(1,2*ep)
      +`Y'(-1+n1,0,1,1,0,1,3,0,0)*rat(-ep-1,ep^2)
      +`Y'(-1+n1,0,1,1,2,0,2,0,0)*rat(1,2*ep)
      +`Y'(-1+n1,1,0,1,0,1,3,0,0)*rat(ep+2,ep^2)
      +`Y'(-1+n1,1,0,1,2,0,2,0,0)*rat(-1,2*ep)
      +`Y'(-1+n1,1,0,2,0,2,1,0,0)*rat(1,2*ep)
      +`Y'(-1+n1,1,0,2,1,2,0,0,0)*rat(-1,2*ep)
      +`Y'(-1+n1,1,0,2,2,0,1,0,0)*rat(-1,2*ep)
      +`Y'(-1+n1,1,0,2,2,1,0,0,0)*rat(1,ep)
      +`Y'(-1+n1,1,1,1,0,1,2,0,0)*rat(ep+1,ep)
      +`Y'(-1+n1,1,1,1,2,0,1,0,0)*rat(-1,1)
      +`Y'(-1+n1,1,1,1,2,1,0,0,0)*rat(2*ep-1,ep)
      +`Y'(-1+n1,2,-1,1,0,1,3,0,0)*rat(-1,ep^2)
      +`Y'(-1+n1,2,0,1,0,1,2,0,0)*rat(-ep+2,2*ep^2)
      +`Y'(-1+n1,2,0,2,1,1,0,0,0)*rat(-1,2*ep)
      +`Y'(-1+n1,2,1,1,0,1,1,0,0)*rat(ep+1,ep)
      +`Y'(-1+n1,2,1,1,1,1,0,0,0)*rat(-2*ep+1,ep)
      +`Y'(1+n1,0,1,2,1,1,1,0,0)*rat(-ep-n1,2*ep)
      +`Y'(1+n1,1,1,2,1,0,1,0,0)*rat(ep+n1,2*ep)
      +`Y'(n1,0,1,1,1,1,2,0,0)*rat(-4,1)
      +`Y'(n1,0,1,1,2,0,2,0,0)*rat(-1,2*ep)
      +`Y'(n1,0,1,2,1,1,1,0,0)*rat(ep+n1-1,2*ep)
      +`Y'(n1,1,0,1,1,1,2,0,0)*rat(3,1)
      +`Y'(n1,1,0,1,2,0,2,0,0)*rat(1,2*ep)
      +`Y'(n1,1,0,2,1,1,1,0,0)*rat(1,1)
      +`Y'(n1,1,0,2,2,0,1,0,0)*rat(1,2*ep)
      +`Y'(n1,1,0,2,2,1,0,0,0)*rat(-1,ep)
      +`Y'(n1,1,1,1,2,0,1,0,0)*rat(1,1)
      +`Y'(n1,1,1,1,2,1,0,0,0)*rat(-2*ep+1,ep)
      +`Y'(n1,1,1,2,0,1,1,0,0)*rat(ep+n1-1,2*ep)
      +`Y'(n1,1,1,2,1,0,1,0,0)*rat(-ep-n1+1,2*ep))
      ;
#endif

GoTo 2;
Label 1;
$ired = 1;
Label 2;
id    `Z'(n1?,?a,x?neg0_,?b,n8?,n9?) = `Y'(n1,?a,x,?b,n8,n9);

ModuleOption,maximum,$ired;
.sort:Main fastar3 loop, pass `$pass';
#if ( `$ired' > 0 )
      #redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
#break
*--#] fastar3 : 
*--#[ haha :
#case haha
*
#do il = 1,15
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering haha loop `il', pass `$pass' at `time_' sec."
#endif
#if ( `il' >= 1 )
#if(1==1)
*--#[ stat file :
*
*  Prework to make the reduction of the haha topology faster.
*
* B8
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n11)*(`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,n8,
   n9,n10,-1+n11,n12,1+n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n11+1,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,1+
   n13,n14)*rat(n9,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,-1+n11,n12,1+n13,n14)*rat(-2*ep-2*n2-n7-n9-n11-n12+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-
   1+n11,1));
* B20
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n11+1)*(`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n11,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,
   1+n10,-1+n11,1+n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(-n8,1)+`Z'(n1,
   n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,-1+n11,1+n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,
   n14)*rat(-n11+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,
   n9,n10,-1+n11,1+n12,n13,n14)*rat(-2*ep-2*n5-n8-n10-n11-n13+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*
   rat(-1+n11,1));
* B10
id,ifmatch->haha1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,n8,
   n9,n10,-1+n11,n12,1+n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,-1+n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,1+n13,n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,
   n14)*rat(-2*n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,1+n13,n14)*rat(2*ep+n2+n7+n9+2*n11+n12-5,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*
   rat(-n2+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-1+n2,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,
   n11,n12,1+n13,n14)*rat(-n2+1,1));
* B18
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n5)*(`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n5+1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,
   n9,n10,n11,1+n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,
   n3,n4,-1+n5,n6,n7,n8,-1+n9,1+n10,n11,1+n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,-1+n11,1+n12,
   n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(2*n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,
   n8,n9,n10,n11,1+n12,n13,n14)*rat(2*ep+n5+n8+n10+2*n11+n13-5,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)
   *rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(-1+n5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,
   n11,1+n12,n13,n14)*rat(-n5+1,1));
* B3
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(1+n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n1,1)+`Z'(1+n1,n2,n3,n4,n5,-1+n6,-1+n7,
   n8,n9,n10,n11,1+n12,n13,n14)*rat(n1,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,1+
   n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+
   n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(
   n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(2*ep+n1+n6+2*n7+n8+n14-5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,
   n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n6,1)+`Z'(n1,n2,n3,n4,
   n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n6+1,1));
* B5
id,ifmatch->haha1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n1+1)*(`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,1+n6,n7,
   n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(-1+
   n1,n2,n3,n4,n5,1+n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n6,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,-1+n8,n9,n10,n11,n12,1+
   n13,n14)*rat(n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,
   n8,n9,-1+n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(-
   1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(2*ep+n1+n6+n7+2*n8+n14-5,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,
   n10,n11,n12,1+n13,n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,
   n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n1+1,1));
* B7
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n9)*(`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+
   n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,1+n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(n2,1)+`Z'(n1,
   n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n9+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,n11,1+n12,n13,
   1+n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n9+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-
   1+n9,n10,1+n11,n12,n13,1+n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n9,1)+`Z'(n1,n2,
   n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,1+n11,n12,n13,1+n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+
   n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(2*ep+n2+2*n7+n9+n11+n12-5,1)+`Z'(n1,n2,n3,
   n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n9,1));
* B9
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n7)*(`Z'(n1,1+n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(n1,1+n2,n3,n4,n5,n6,-1+n7,n8,-
   1+n9,n10,n11,n12,n13,1+n14)*rat(n2,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,
   n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,
   1+n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-
   1+n9,n10,1+n11,n12,n13,1+n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,-1+n10,1+n11,n12,n13,1+n14)*rat(-n11,1)+`Z'(
   n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+
   n14)*rat(2*ep+n2+n7+2*n9+n11+n12-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-1+n7,1)+`Z'(n1,n2,n3,
   n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n7,1));
* B13
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n4)*(`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,
   1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,1+
   n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,
   -1+n9,1+n10,n11,1+n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,-1+n11,1+n12,n13,n14)*rat(-n10,1)+`Z'(
   n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,
   n14)*rat(2*ep+n3+n4+2*n9+n10+n14-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-1+n4,1)+`Z'(n1,n2,n3,
   n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n4+1,1));
* B15
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n3+1)*(`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,
   n9,-1+n10,n11,n12,1+n13,n14)*rat(n4,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,
   n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+
   n13,1+n14)*rat(n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,-1+n10,n11,n12,1+n13,n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,
   n7,n8,1+n9,n10,-1+n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(2*ep+n3+
   n4+n9+2*n10+n14-5,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,
   n9,n10,n11,n12,1+n13,n14)*rat(-n3+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,n14)*rat(-1+n3,1)+`Z'(n1,n2,
   n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n3+1,1));
* B17
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?{>1},n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n10)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(-n10+1,1)+`Z'(n1,n2,n3,n4,1+n5,-1+n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-n5,1)+`Z'(n1,
   n2,n3,n4,1+n5,n6,n7,-1+n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(n5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,-1+n10,n11,n12,1+
   n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n10+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,
   n8,n9,-1+n10,1+n11,n12,n13,1+n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,-1+n10,1+n11,n12,n13,1+n14)*rat(n11,1)+
   `Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+
   n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(2*ep+n5+2*n8+n10+n11+n13-5,1)+`Z'(n1,
   n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n10,1));
* B19
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n8)*(`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(-n8+1,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5,1)+`Z'(n1,
   n2,n3,-1+n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,1+n5,n6,n7,-1+n8,n9,-1+n10,n11,n12,
   n13,1+n14)*rat(n5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+
   n8,-1+n9,n10,1+n11,n12,n13,1+n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,-1+n10,1+n11,n12,n13,1+n14)*rat(n11,1)+
   `Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,
   1+n14)*rat(2*ep+n5+n8+2*n10+n11+n13-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-1+n8,1)+`Z'(n1,n2,
   n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n8,1));
* B2
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n8)*(`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,-1+
   n8,n9,n10,n11,n12,1+n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n8+1,1)+`Z'(n1,-1+
   n2,n3,n4,n5,n6,1+n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*
   rat(-1+n8,1)+`Z'(n1,n2,n3,n4,n5,1+n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n6,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,
   n11,n12,1+n13,n14)*rat(-2*ep-2*n1-n6-n7-n8-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n8,1))
;
* B14
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n9+1)*(`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n9,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+
   n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,n2,n3,-
   1+n4,n5,n6,n7,n8,-1+n9,1+n10,n11,1+n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*
   rat(-n9+1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,1+n10,n11,1+n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,
   n10,n11,1+n12,n13,n14)*rat(-2*ep-n3-2*n4-n9-n10-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+
   n9,1));
* B6
id,ifmatch->haha1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(1+n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,n14)*rat(-n1,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+
   n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,
   n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(8*ep+n1+2*n2+n3+2*
   n4+2*n5+2*n6+2*n7+2*n8+2*n9+2*n10+2*n11+n12+n13+2*n14-17,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*
   rat(-2*n12-2,1));
* B16
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n5+1)*(`Z'(-1+n1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,
   n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n4,1)+`Z'(n1,
   n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n4,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+
   n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n6,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,
   n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-8*ep-2*n1-2*n2-2*n3-
   n4-2*n5-n6-2*n7-2*n8-2*n9-2*n10-2*n11-n12-n13-2*n14+17,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,2+n13,n14)*rat(
   -2*n13-2,1));
* B1
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n7)*(`Z'(1+n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n1,1)+`Z'(1+n1,n2,n3,n4,n5,n6,-1+n7,n8,
   n9,n10,n11,1+n12,n13,n14)*rat(-n1,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n7+1,1)+`Z'(n1,n2,-1+
   n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+
   n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,
   1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n7,1)+`Z'(n1,n2,
   n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,n14)*
   rat(2*ep+n1+2*n6+n7+n8+n14-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n7+1,1));
* B11
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?{>1},n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n10)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   1+n9,-1+n10,n11,n12,1+n13,n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,n14)*rat(-n4,1)+`Z'(n1,
   n2,-1+n3,n4,n5,n6,n7,n8,1+n9,-1+n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,
   n14)*rat(-n10+1,1)+`Z'(n1,n2,n3,1+n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,n14)*rat(n4,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,
   n9,n10,n11,n12,1+n13,n14)*rat(-1+n10,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,-1+n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,
   n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,
   n14)*rat(-2*ep-2*n3-n4-n9-n10-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n10,1));
* B4
id,ifmatch->haha1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n1+1)*(`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,-
   1+n9,n10,n11,n12,n13,1+n14)*rat(n2,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n6,1)+`Z'(-1+n1,
   n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,
   1+n14)*rat(-n13,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(-n14-1,1)+`Z'(-1+n1,n2,n3,-1+n4,1+n5,
   n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n6,1)+
   `Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,
   n12,1+n13,1+n14)*rat(n13,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(-n14-1,1)+`Z'(-1+n1,n2,n3,n4,
   1+n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(n5,1)+`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-
   n6,1)+`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+
   n12,n13,1+n14)*rat(-n12,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(-1+n1,n2,n3,n4,n5,
   n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(6*ep+n1+n2+n5+2*n6+2*n7+2*n8+2*n9+2*n10+2*n11+n12+n13+n14-12,1)+`Z'(-1+n1,n2,n3,
   n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(1+n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+
   n1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n1+1,1));
* B12
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n3+1)*(`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(-1+
   n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(-n14-1,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+
   n14)*rat(-n3+1,1)+`Z'(n1,1+n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n2,1)+`Z'(n1,n2,-2+n3,1+n4,n5,n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,-
   1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*
   rat(-n5,1)+`Z'(n1,n2,-1+n3,n4,1+n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(n5,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,
   n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,-
   1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(-n14-1,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+
   n14)*rat(-n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,
   n10,n11,n12,n13,1+n14)*rat(6*ep+n2+2*n3+n4+n5+2*n7+2*n8+2*n9+2*n10+2*n11+n12+n13+n14-13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,
   n8,n9,n10,n11,n12,n13,2+n14)*rat(1+n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n3,1)+`Z'(n1,
   n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n3+1,1));
*-----------------------------------------------------------------------------------------
* G234
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n3)*(`Z'(n1,-1+n2,-1+n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,1+n14)*rat(-n9*n14,-1+n4)+`Z'(n1,-1+n2,n3,-1+n4,
   n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n3*n9+n9,-1+n4)+`Z'(n1,n2,-1+n3,-1+n4,n5,-1+n6,n7,n8,1+n9,n10,n11,n12,n13,1+
   n14)*rat(-n9*n14,-1+n4)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,-1+n7,n8,1+n9,n10,n11,n12,n13,1+n14)*rat(n9*n14,-1+n4)+`Z'(n1,n2,-1+
   n3,-1+n4,n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,n13,n14)*rat(-n9*n10,-1+n4)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,
   n12,n13,1+n14)*rat(n9*n14,-1+n4)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(2*ep*n9+n3*n9+2*n9^2+
   n9*n10+n9*n14-3*n9,-1+n4)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n4*n10+n9*n10+n10,-1+n4)+`Z'(
   n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,
   n14)*rat(-2*ep-n3-2*n4-n10-n14+5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n3*n4+n3*n9+n3+n4-n9-1,
   -1+n4));
* G282
id,ifmatch->haha1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(-1+n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8*n14,-1+n1)+`Z'(-1+n1,n2,n3,n4,-1+
   n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8*n14,-1+n1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,
   n14)*rat(-n6*n8+n8,-1+n1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,1+n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(-n7*n8,-1+n1)+`Z'(-1+n1,
   n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n1*n7+n7*n8+n7,-1+n1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,-1+
   n10,n11,n12,n13,1+n14)*rat(n8*n14,-1+n1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(n8*n14,-1+
   n1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(2*ep*n8+n6*n8+n7*n8+2*n8^2+n8*n14-3*n8,-1+n1)+`Z'(-1+
   n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n1*n6+n1+n6*n8+n6-n8-1,-1+n1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,
   n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep-2*n1-n6-n7-n14+5,1))
;
* G82
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(-1+n1,n2,n3,n4,n5,-1+n6,1+n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(-n7*n8,-1+n11)+`Z'(-1+n1,n2,n3,n4,n5,-
   1+n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+
   `Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,1+n8,n9,n10,-1+n11,
   n12,n13,n14)*rat(n7*n8,-1+n11)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,n4,
   n5,-1+n6,n7,1+n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(n8*n9,-1+n11)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,
   n14)*rat(n8,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(-n8*n9,-1+n11)+`Z'(n1,n2,n3,n4,n5,-1+
   n6,n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(2*ep*n8+2*n2*n8+n7*n8+n8*n9+n8*n11+n8*n12-5*n8,-1+n11)+`Z'(n1,n2,n3,n4,n5,-1+
   n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep-2*n1-n6-n7-n8-n14+5,1));
* G94
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n3)*(`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,n13,n14)*rat(-n9*n10,-1+n11)+`Z'(n1,n2,-1+n3,-1+n4,
   n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n10,
   1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(n8*n9,-1+n11)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,
   1+n9,1+n10,-1+n11,n12,n13,n14)*rat(n9*n10,-1+n11)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)
   +`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,1+n9,n10,-1+
   n11,n12,n13,n14)*rat(-n8*n9,-1+n11)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(2*ep*n9+2*n5*n9+n8*
   n9+n9*n10+n9*n11+n9*n13-5*n9,-1+n11)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep-n3-2*n4-n9-n10-
   n14+5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n3+1,1));
* G214
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n3)*(`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,
   n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,-1+n7,1+n8,1+n9,n10,n11,n12,n13,n14)*rat(-n8*n9,-1+
   n5)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(n8*n9,-1+n5)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,
   1+n9,1+n10,-1+n11,n12,n13,n14)*rat(n9*n10,-1+n5)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,-1+n12,1+n13,n14)*rat(-
   n9*n13,-1+n5)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,1+n13,n14)*rat(2*n9*n13,-1+n5)+`Z'(n1,n2,-1+n3,n4,-1+n5,
   n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(2*ep*n9+n5*n9+n8*n9+n9*n10+2*n9*n11+n9*n13-5*n9,-1+n5)+`Z'(n1,n2,-1+n3,n4,-1+n5,
   n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n5*n10-n9*n10-n10,-1+n5)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,n13,
   n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep-n3-2*n4-n9-n10-n14+5,1)+`Z'(n1,n2,n3,-1+
   n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n3+1,1));
* G262
id,ifmatch->haha1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(-1+n1,n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,1+n8,
   n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,-1+n2,n3,
   n4,n5,-1+n6,1+n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(n7*n8,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,
   n13,n14)*rat(n2*n7-n7*n8-n7,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,1+n9,-1+n10,n11,n12,n13,n14)*rat(-n8*n9,-1+n2)+
   `Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(n8*n9,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,n9,
   n10,n11,1+n12,-1+n13,n14)*rat(-n8*n12,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*n8*
   n12,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(2*ep*n8+n2*n8+n7*n8+n8*n9+2*n8*n11+n8*n12-5*
   n8,-1+n2)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,
   n11,n12,n13,n14)*rat(-2*ep-2*n1-n6-n7-n8-n14+5,1));
* G81
id,ifmatch->haha1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,1-n1)*(`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,
   n9,1+n10,-1+n11,n12,n13,n14)*rat(n7*n10,-1+n11)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)+
   `Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(-n7*n8,-1+n11)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,
   n9,1+n10,-1+n11,n12,n13,n14)*rat(-n7*n10,-1+n11)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)
   +`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,1+n7,1+n8,n9,n10,-1+
   n11,n12,n13,n14)*rat(n7*n8,-1+n11)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(-1+n1,n2,n3,
   n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,n12,n13,n14)*rat(-
   2*ep*n7-2*n5*n7-n7*n8-n7*n10-n7*n11-n7*n13+5*n7,-1+n11)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(
   n14,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(2*ep+n1+2*n6+n7+n8+n14-5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,
   n7,n8,n9,n10,n11,n12,n13,n14)*rat(-1+n1,1));
* G91
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n4)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(-n7*n10,-1+n11)+`Z'(-1+n1,n2,n3,-1+n4,
   n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(
   n7*n10,-1+n11)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,n13,n14)*rat(n9*n10,-1+n11)+`Z'(n1,-1+n2,n3,-1+n4,
   n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)+
   `Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,n13,n14)*rat(-n9*n10,-1+n11)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,1+
   n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,-1+
   n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(
   n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(2*ep*n10+2*n2*n10+n7*n10+n9*n10+n10*n11+n10*n12-
   5*n10,-1+n11)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,
   n11,n12,n13,n14)*rat(-2*ep-2*n3-n4-n9-n10-n14+5,1));
* G201
id,ifmatch->haha1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,1-n1)*(`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(-n7*n8,-1+n5)
   +`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,-1+n9,1+n10,n11,n12,n13,n14)*rat(n7*n10,-1+n5)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,
   n9,1+n10,-1+n11,n12,n13,n14)*rat(-n7*n10,-1+n5)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,-1+n12,1+n13,n14)*rat(
   n7*n13,-1+n5)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-2*n7*n13,-1+n5)+`Z'(-1+n1,n2,n3,n4,-1+
   n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep*n7-n5*n7-n7*n8-n7*n10-2*n7*n11-n7*n13+5*n7,-1+n5)+`Z'(-1+n1,n2,n3,n4,-
   1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n5*n8+n7*n8+n8,-1+n5)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,
   n13,n14)*rat(n7,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,
   n8,n9,n10,-1+n11,n12,n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(-1+n1,
   n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(2*ep+n1+2*n6+n7+n8+n14-5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,
   n12,n13,n14)*rat(-1+n1,1));
* G271
id,ifmatch->haha1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n4)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,-
   1+n8,n9,1+n10,n11,n12,n13,n14)*rat(-n7*n10,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(
   n7*n10,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,n13,n14)*rat(n9*n10,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,
   n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n2*n9-n9*n10-n9,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,-
   1+n13,n14)*rat(-n10*n12,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-2*n10*n12,-1+n2)+`Z'(
   n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(2*ep*n10+n2*n10+n7*n10+n9*n10+2*n10*n11+n10*n12-5*n10,-1+
   n2)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,1+n10,
   n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,n3,-1+n4,n5,-
   1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(n10,1)+
   `Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,
   n14)*rat(-2*ep-2*n3-n4-n9-n10-n14+5,1));
* G221
id,ifmatch->haha1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,1-n1)*(`Z'(-1+n1,-1+n2,1+n3,-1+n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n3*n7,-1+n4)+`Z'(-1+n1,-1+n2,n3,-1+n4,
   n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n7*n14,-1+n4)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*
   rat(-n14,1)+`Z'(-1+n1,n2,1+n3,-1+n4,n5,n6,1+n7,n8,-1+n9,n10,n11,n12,n13,n14)*rat(-n3*n7,-1+n4)+`Z'(-1+n1,n2,n3,-1+n4,n5,-
   1+n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n7*n14,-1+n4)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,-1+n9,1+n10,n11,n12,n13,
   n14)*rat(-n7*n10,-1+n4)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(n7*n10,-1+n4)+`Z'(-1+n1,n2,
   n3,-1+n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7*n14,-1+n4)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,n9,n10,n11,n12,
   n13,n14)*rat(-2*ep*n7-n3*n7-n4*n7-2*n7*n9-n7*n10-n7*n14+5*n7,-1+n4)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,
   n13,1+n14)*rat(n4*n14-n7*n14-n14,-1+n4)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(-1+n1,
   n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)
   *rat(n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,n8,-1+n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,
   n11,n12,n13,1+n14)*rat(n14,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(2*ep+n1+2*n6+n7+n8+n14-5,1)+`Z'(
   n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-1+n1,1));
* G194
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n3)*(`Z'(1+n1,-1+n2,-1+n3,n4,n5,-1+n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n1*n9,-1+n6)+`Z'(1+n1,n2,-1+n3,n4,n5,
   -1+n6,-1+n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n1*n9,-1+n6)+`Z'(n1,-1+n2,-1+n3,n4,n5,-1+n6,n7,n8,1+n9,n10,n11,n12,n13,1+
   n14)*rat(-n9*n14,-1+n6)+`Z'(n1,n2,-2+n3,n4,n5,-1+n6,n7,n8,1+n9,n10,n11,n12,n13,1+n14)*rat(-n9*n14,-1+n6)+`Z'(n1,n2,-1+n3,
   -1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(
   -n10,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,-1+n7,1+n8,1+
   n9,n10,n11,n12,n13,n14)*rat(n8*n9,-1+n6)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(-n8*n9,-
   1+n6)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,1+n9,n10,n11,-1+n12,n13,1+n14)*rat(n9*n14,-1+n6)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,
   n8,1+n9,n10,n11,n12,n13,n14)*rat(2*ep*n9+n1*n9+n6*n9+2*n7*n9+n8*n9+n9*n14-5*n9,-1+n6)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(n9*n14,-1+n6)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(
   n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep-n3-2*n4-n9-n10-n14+5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,
   n10,n11,n12,n13,n14)*rat(-n3+1,1));
* G21
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n6)*(`Z'(1+n1,n2,n3,n4,n5,-2+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(n1,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,n9,
   n10,n11,n12,n13,1+n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,1+n2,-1+
   n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,-1+n12,n13,1+n14)*rat(n2,1)+`Z'(n1,1+n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,n11,-1+n12,n13,
   1+n14)*rat(-n2,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,-1+n12,n13,1+n14)*rat(n7,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,
   1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+
   `Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n12-n14-1,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,
   n12,n13,1+n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n5,1)+`Z'(n1,n2,n3,-1+n4,
   n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(
   -n12+n14+1,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(2*n14+2,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(-1+n6,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8,1)+`Z'(n1,
   n2,n3,n4,-1+n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,n3,n4,1+n5,-1+n6,n7,n8,n9,-1+n10,n11,n12,n13,
   1+n14)*rat(-n5,1)+`Z'(n1,n2,n3,n4,n5,-2+n6,1+n7,n8,n9,n10,n11,-1+n12,n13,1+n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,n5,-2+n6,1+n7,
   n8,n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,n5,-2+n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(n8,1)+`Z'(n1,n2,n3,
   n4,n5,-2+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,1+n7,n8,-1+n9,n10,n11,-1+n12,n13,1+n14)*
   rat(-n7,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,1+n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,
   n10,1+n11,-1+n12,n13,1+n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,1+n11,n12,n13,1+n14)*rat(n11,1)+`Z'(n1,
   n2,n3,n4,n5,-1+n6,n7,n8,n9,-1+n10,1+n11,-1+n12,n13,1+n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,-1+n10,1+n11,
   n12,n13,1+n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,-1+n12,n13,1+n14)*rat(-2*ep-n2-n7-2*n9-n11-n12+5,
   1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,
   n13,1+n14)*rat(-2*ep-n5-n8-2*n10-n11+n12-n13+n14+3,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(2*ep+
   n1+2*n6+n7+n8+n14-6,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n6,1));
* G31
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n4)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(n13-n14-1,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(-n14-1,1)+`Z'(
   n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,1+n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,
   n13,n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,1+n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,
   n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,-1+n13,1+n14)*rat(n10,1)+
   `Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,
   n13,1+n14)*rat(-n4+1,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,
   n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)+`Z'(
   n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,-1+n13,1+n14)*rat(n5,1)+`Z'(n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,n9,n10,n11,
   n12,n13,1+n14)*rat(n5,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,-1+n8,n9,n10,n11,n12,-1+n13,1+n14)*rat(-n5,1)+`Z'(n1,n2,n3,-1+n4,
   1+n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,-1+n13,1+n14)*
   rat(n10,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,
   n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n13+n14+1,1)+`Z'(
   n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(1+n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,1+n11,
   n12,-1+n13,1+n14)*rat(n11,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(n11,1)+`Z'(n1,n2,n3,-1+n4,
   n5,n6,n7,-1+n8,n9,1+n10,n11,n12,-1+n13,1+n14)*rat(-n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,-1+n8,n9,1+n10,n11,n12,n13,1+n14)
   *rat(-n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,-1+n13,1+n14)*rat(-n11,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,-1+
   n8,n9,n10,1+n11,n12,n13,1+n14)*rat(-n11,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,-1+n13,1+n14)*rat(-2*ep-n5-2*
   n8-n10-n11-n13+5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(-4*ep-2*n3-n4-n5-2*n8-n9-2*n10-n11-2*n13-2*n14+9,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,
   n10,n11,n12,n13,2+n14)*rat(-n14-1,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep-2*n3-n4-n9-n10-
   n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n4,1));
* G86
id,ifmatch->haha1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,1-n2)*(`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(-n7*n10,-1+n11)+`Z'(n1,-1+n2,n3,-1+n4,
   n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,n13,n14)*rat(-n9*n10,-1+n11)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,
   n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,
   1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(n7*n8,-1+n11)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(
   n7*n10,-1+n11)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,1+n8,
   1+n9,n10,-1+n11,n12,n13,n14)*rat(n8*n9,-1+n11)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(
   n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,n13,n14)*rat(n9*n10,-1+n11)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,
   n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,-1+n2,n3,n4,
   n5,-1+n6,1+n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(-n7*n8,-1+n11)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,
   n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(-n8*n9,-1+n11)+`Z'(n1,-1+n2,n3,
   n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,n12,n13,n14)*rat(
   2*ep*n7+2*n5*n7+n7*n8+n7*n10+n7*n11+n7*n13-5*n7,-1+n11)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,n13,n14)*
   rat(2*ep*n9+2*n5*n9+n8*n9+n9*n10+n9*n11+n9*n13-5*n9,-1+n11)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*
   rat(-2*n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(4*ep+n2+2*n5+n7+n8+n9+n10+2*n11+n12+n13-9,1)+
   `Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(-n2*n10+n10,-1+n11)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,
   n10,-1+n11,n12,n13,n14)*rat(n2*n8-n8,-1+n11)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(n2*n10-
   n10,-1+n11)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-1+n2,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,-
   1+n11,n12,n13,n14)*rat(-n2*n8+n8,-1+n11)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,n14)*rat(2*ep*n2-2*ep+2*n2*
   n5+n2*n8+n2*n10+n2*n11+n2*n13-5*n2-2*n5-n8-n10-n11-n13+5,-1+n11));
* G96
id,ifmatch->haha1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n5)*(`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(n7*n8,-1+n11)+`Z'(-1+n1,n2,n3,n4,-1+
   n5,n6,1+n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(n7*n10,-1+n11)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,
   n14)*rat(n7,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,n8,n9,
   n10,-1+n11,n12,n13,n14)*rat(n5*n7-n7,-1+n11)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(-n7*
   n8,-1+n11)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(-n7*n10,-1+n11)+`Z'(n1,-1+n2,n3,n4,-1+n5,
   n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,1+n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(-n8*
   n9,-1+n11)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,
   1+n10,-1+n11,n12,n13,n14)*rat(-n9*n10,-1+n11)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(
   n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,n12,
   n13,n14)*rat(-n5*n7+n7,-1+n11)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(-n5*n9+n9,-1+n11)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n5+1,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,1+n9,n10,-1+n11,n12,
   n13,n14)*rat(n8*n9,-1+n11)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,n13,n14)*rat(n9*n10,-1+n11)+`Z'(n1,n2,
   -1+n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)
   *rat(n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,n13,n14)*rat(n5*n9-n9,-1+n11)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,
   1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(-2*ep*n8-2*n2*n8-n7*n8-n8*n9-n8*n11-n8*n12+5*n8,-1+n11)+`Z'(n1,n2,n3,n4,-1+n5,n6,
   n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(-2*ep*n10-2*n2*n10-n7*n10-n9*n10-n10*n11-n10*n12+5*n10,-1+n11)+`Z'(n1,n2,n3,n4,-
   1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-2*n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-4*ep-
   2*n2-n5-n7-n8-n9-n10-2*n11-n12-n13+9,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,n14)*rat(-2*ep*n5+2*ep-2*n2*
   n5+2*n2-n5*n7-n5*n9-n5*n11-n5*n12+5*n5+n7+n9+n11+n12-5,-1+n11));
* G206
id,ifmatch->haha1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,1-n2)*(`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,
   1+n10,n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,-1+n7,1+n8,1+n9,n10,n11,n12,n13,n14)*rat(-n8*n9,-1+n5)+
   `Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(n7*n8,-1+n5)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,-1+
   n9,1+n10,n11,n12,n13,n14)*rat(-n7*n10,-1+n5)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(n7*
   n10,-1+n5)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,-1+n12,1+n13,n14)*rat(-n7*n13,-1+n5)+`Z'(n1,-1+n2,n3,n4,-1+n5,
   n6,1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(2*n7*n13,-1+n5)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*
   rat(2*ep*n7+n5*n7+n7*n8+n7*n10+2*n7*n11+n7*n13-5*n7,-1+n5)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,1+n8,1+n9,n10,-1+n11,n12,n13,
   n14)*rat(n8*n9,-1+n5)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(n5*n8-n7*n8-n8,-1+n5)+`Z'(n1,-1+n2,
   n3,n4,-1+n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,n13,n14)*rat(n9*n10,-1+n5)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,-
   1+n12,1+n13,n14)*rat(-n9*n13,-1+n5)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,1+n13,n14)*rat(2*n9*n13,-1+n5)+
   `Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(2*ep*n9+n5*n9+n8*n9+n9*n10+2*n9*n11+n9*n13-5*n9,-1+n5)+
   `Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n5*n10-n9*n10-n10,-1+n5)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+
   n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,n13,
   n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,
   n10,n11,n12,n13,n14)*rat(4*ep+n2+2*n5+n7+n8+n9+n10+2*n11+n12+n13-9,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,n9,n10,n11,
   n12,n13,n14)*rat(-n2*n8+n8,-1+n5)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,-1+n11,n12,n13,n14)*rat(n2*n8-n8,-1+n5)+`Z'(n1,
   n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,1+n10,n11,n12,n13,n14)*rat(-n2*n10+n10,-1+n5)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,-1+
   n11,n12,n13,n14)*rat(n2*n10-n10,-1+n5)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,n14)*rat(-n2*n13+n13,-1+
   n5)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(2*n2*n13-2*n13,-1+n5)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,
   n10,n11,n12,n13,n14)*rat(2*ep*n2-2*ep+n2*n5+n2*n8+n2*n10+2*n2*n11+n2*n13-5*n2-n5-n8-n10-2*n11-n13+5,-1+n5)+`Z'(n1,n2,
   n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,n14)*rat(-1+n2,1));
*--#] stat file : 
#endif
#endif
#if ( `il' >= 2 )
*--#[ n14 :
`SKIP5'
id,ifmatch->haha1,`Z'(n1?`POS1',n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`POS9',n10?`POS10',n11?`POS11',n12?`NEG12',n13?`NEG13',n14?neg_)
			 = -rat(1,-4*ep-3*n1-2*n6-2*n7-2*n8-2*(n14+1)+9)*(
      +`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,-1+(n14+1))*rat(-n6,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,n13,(n14+1))*rat(2*ep+n2+n7+2*n9+n11+n12-5,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,-1+n13,(n14+1))*rat(-2*ep-n5-n8-2*n10-n11-n13+5,1)
      +`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(-n6,1)
      +`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,(n14+1))*rat(-n3,1)
      +`Z'(n1,n2,n3,1+n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,(n14+1))*rat(-n4,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(-2*ep-2*n1+n3+n4-n6-n7-n8-n12-n13-2*(n14+1)+5,1)
      +`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(-n2,1)
      +`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,(n14+1))*rat(n2,1)
      +`Z'(-1+n1,n2,-1+n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(n6,1)
      +`Z'(-1+n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(-n5,1)
      +`Z'(-1+n1,n2,n3,-1+n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(-n6,1)
      +`Z'(-1+n1,n2,n3,n4,1+n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,(n14+1))*rat(n5,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(4*ep+n2+n5+n6+n7+n8+2*n9+2*n10+2*n11+n12+n13-9,1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,-1+n10,n11,n12,n13,(n14+1))*rat(-n9,1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(-4*ep-2*n2-2*n7-3*n9-2*n11-2*n12+9,1)
      +`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,(n14+1))*rat(n4,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,-1+n10,n11,n12,n13,(n14+1))*rat(n9,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(-4*ep+2*n1-n3-3*n4-2*n5-3*n9-2*n10-2*n11-2*n12+(n14+1)+9,1)
      +`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,(n14+1))*rat(n3,1)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,n13,(n14+1))*rat(n10,1)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(-8*ep-2*n1-2*n2-3*n3-n4-2*n7-2*n8-2*n9-3*n10-2*n11-2*
      n13-3*(n14+1)+17,1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,1+n10,n11,n12,n13,(n14+1))*rat(-n10,1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(-4*ep-2*n5-2*n8-3*n10-2*n11-2*n13+9,1)
      +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,(n14+1))*rat(2*ep+2*n1+n6+n7+n8+(n14+1)-5,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,(n14+1))*rat(6*ep+2*n2+n3+2*n4+2*n7+3*n9+n10+2*n11+2*n12+(n14+1)-13,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,(n14+1))*rat(6*ep+2*n3+n4+2*n5+2*n8+n9+3*n10+2*n11+2*n13+(n14+1)-13,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,(n14+1))*rat(-n12,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(n14+1))*rat(n13,1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,(n14+1))*rat(2*n12,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,(n14+1))*rat(2*n12,1)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(n14+1))*rat(-2*n13,1)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+(n14+1))*rat(-3*(n14+1),1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(n14+1))*rat(-2*n13,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,(n14+1))*rat(-2*n12,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,(n14+1))*rat(2*n13,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+(n14+1))*rat((n14+1),1)
      +`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,(n14+1))*rat(-n12,1)
      +`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(n14+1))*rat(-n13,1)
      +`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,(n14+1))*rat(n12,1)
      +`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(n14+1))*rat(n13,1)
      +`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+(n14+1))*rat(-(n14+1),1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+(n14+1))*rat((n14+1),1)
      +`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,(n14+1))*rat(2*n12,1)
      +`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,(n14+1))*rat(-2*n12,1)
      +`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(n14+1))*rat(2*n13,1)
      +`Z'(n1,n2,n3,-1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(n14+1))*rat(-2*n13,1)
      +`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+(n14+1))*rat(3*(n14+1),1)
      +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,-1+n10,n11,n12,n13,1+(n14+1))*rat(-(n14+1),1)
      );
*--#] n14 : 
#endif
#if ( `il' >= 3 )
*--#[ n13 :
`SKIP5'
id,ifmatch->haha1,`Z'(n1?`POS1',n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`POS9',n10?`POS10',n11?`POS11',n12?`NEG12',n13?neg_,0)
			 = -rat(1,-2*ep-2*n5-n8-n10-n11-(1+n13)+5)*(
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,(1+n13),0)*rat(-2*ep-2*n2-n7-n9-n11-n12+5,1)
      +`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n1,1)
      +`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n3,1)
      +`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(n4,1)
      +`Z'(n1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(n6,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(n1+2*n2+n3-n4-2*n5-n6+n7-n8+n9-n10+n12-(1+n13),1)
      +`Z'(-1+n1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n6,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-2*ep-n1-n6-n7-n8+4,1)
      +`Z'(1+n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(n1,1)
      +`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(n3,1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(8*ep+n1+2*n2+n3+2*n4+2*n5+2*n6+2*n7+2*n8+2*n9+2*n10+2*
      n11+n12+(1+n13)-18,1)
      +`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n4,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-2*ep-n3-n4-n9-n10+4,1)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(2*ep+n3+n4+n9+n10-4,1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-8*ep-2*n1-2*n2-2*n3-n4-2*n5-n6-2*n7-2*n8-2*n9-2*n10-2*
      n11-n12-(1+n13)+18,1)
      +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(2*ep+n1+n6+n7+n8-4,1)
      );
*--#] n13 : 
#endif
#if ( `il' >= 4 )
*--#[ n12 :
`SKIP5'
id,ifmatch->haha1,`Z'(n1?`POS1',n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`POS9',n10?`POS10',n11?`POS11',n12?neg_,0,0)
			 = -rat(2*ep+2*n5+n8+n10+n11-5,
      -20*ep^2-4*ep*n1-8*ep*n2-4*ep*n3-2*ep*n4-22*ep*n5-2*ep*
      n6-6*ep*n7-14*ep*n8-6*ep*n9-14*ep*n10-18*ep*n11-4*ep*(n12+1)+96*ep-4*n1*n5-2*n1*n8-2*n1*n10-2*n1*n11+10*n1-6*n2*n5-4*
      n2*n8-4*n2*n10-6*n2*n11+20*n2-4*n3*n5-2*n3*n8-2*n3*n10-2*n3*n11+10*n3-2*n4*n5-n4*n8-n4*n10-n4*n11+5*n4-4*n5^2-2*
      n5*n6-5*n5*n7-6*n5*n8-5*n5*n9-6*n5*n10-7*n5*n11-3*n5*(n12+1)+51*n5-n6*n8-n6*n10-n6*n11+5*n6-3*n7*n8-3*n7*n10-4*n7*n11
      +15*n7-2*n8^2-3*n8*n9-4*n8*n10-5*n8*n11-2*n8*(n12+1)+33*n8-3*n9*n10-4*n9*n11+15*n9-2*n10^2-5*n10*n11-2*n10*(n12+1)+33*n10
      -4*n11^2-3*n11*(n12+1)+43*n11+10*(n12+1)-115)*(
      +`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(n1*n5-n1*n11,2*ep+2*n5+n8+n10+n11-5)
      +`Z'(1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,(n12+1),0,0)*rat(-n1,1)
      +`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(n3*n5-n3*n11,2*ep+2*n5+n8+n10+n11-5)
      +`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,(n12+1),0,0)*rat(-n3,1)
      +`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(2*ep*n4+n4*n5+n4*n8+n4*n10+2*n4*n11-5*n4,2*ep+2*n5+n8+
      n10+n11-5)
      +`Z'(n1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(2*ep*n6+n5*n6+n6*n8+n6*n10+2*n6*n11-5*n6,2*ep+2*n5+n8+
      n10+n11-5)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(8*ep^2+4*ep*n1+8*ep*n2+4*ep*n3-2*ep*n4+8*ep*n5-2*ep*n6+2*ep*
      n7+4*ep*n8+2*ep*n9+4*ep*n10+8*ep*n11+4*ep*(n12+1)-40*ep+3*n1*n5+2*n1*n8+2*n1*n10+3*n1*n11-10*n1+6*n2*n5+4*n2*n8+4*n2*
      n10+6*n2*n11-20*n2+3*n3*n5+2*n3*n8+2*n3*n10+3*n3*n11-10*n3-n4*n5-n4*n8-n4*n10-2*n4*n11+5*n4+2*n5^2-n5*n6+n5*n7+n5
      *n8+n5*n9+n5*n10+2*n5*n11+3*n5*(n12+1)-20*n5-n6*n8-n6*n10-2*n6*n11+5*n6+n7*n8+n7*n10+2*n7*n11-5*n7+n8*n9+n8*n11+2*n8*
      (n12+1)-10*n8+n9*n10+2*n9*n11-5*n9+n10*n11+2*n10*(n12+1)-10*n10+2*n11^2+3*n11*(n12+1)-20*n11-10*(n12+1)+50,2*ep+2*n5+n8+n10+n11-5)
      +`Z'(-1+n1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(-2*ep*n6-n5*n6-n6*n8-n6*n10-2*n6*n11+5*n6,2*ep+2*n5+
      n8+n10+n11-5)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,(n12+1),0,0)*rat(n7,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(2*ep*n5-2*ep*n11+n1*n5-n1*n11+n5*n6+n5*n7+n5*n8-4*n5-n6*
      n11-n7*n11-n8*n11+4*n11,2*ep+2*n5+n8+n10+n11-5)
      +`Z'(1+n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(2*ep*n1+n1*n5+n1*n8+n1*n10+2*n1*n11-5*n1,2*ep+2*n5+
      n8+n10+n11-5)
      +`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(2*ep*n3+n3*n5+n3*n8+n3*n10+2*n3*n11-5*n3,2*ep+2*n5+
      n8+n10+n11-5)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,(n12+1),0,0)*rat(-n7,1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,(n12+1),0,0)*rat(-n9,1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(32*ep^2+4*ep*n1+8*ep*n2+4*ep*n3+6*ep*n4+32*ep*n5+6*ep*n6+
      10*ep*n7+24*ep*n8+10*ep*n9+24*ep*n10+32*ep*n11+4*ep*(n12+1)-152*ep+3*n1*n5+2*n1*n8+2*n1*n10+3*n1*n11-10*n1+6*n2*n5+4*
      n2*n8+4*n2*n10+6*n2*n11-20*n2+3*n3*n5+2*n3*n8+2*n3*n10+3*n3*n11-10*n3+4*n4*n5+3*n4*n8+3*n4*n10+5*n4*n11-15*n4+6*
      n5^2+4*n5*n6+8*n5*n7+10*n5*n8+8*n5*n9+10*n5*n10+12*n5*n11+3*n5*(n12+1)-74*n5+3*n6*n8+3*n6*n10+5*n6*n11-15*n6+5*n7*n8+
      5*n7*n10+7*n7*n11-25*n7+4*n8^2+5*n8*n9+8*n8*n10+10*n8*n11+2*n8*(n12+1)-56*n8+5*n9*n10+7*n9*n11-25*n9+4*n10^2+10*n10*
      n11+2*n10*(n12+1)-56*n10+6*n11^2+3*n11*(n12+1)-74*n11-10*(n12+1)+180,2*ep+2*n5+n8+n10+n11-5)
      +`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(-2*ep*n4-n4*n5-n4*n8-n4*n10-2*n4*n11+5*n4,2*ep+2*n5+
      n8+n10+n11-5)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,(n12+1),0,0)*rat(n9,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(2*ep*n5-2*ep*n11+n3*n5-n3*n11+n4*n5-n4*n11+n5*n9+n5*n10-4
      *n5-n9*n11-n10*n11+4*n11,2*ep+2*n5+n8+n10+n11-5)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(4*ep^2+2*ep*n3+2*ep*n4+2*ep*n5+2*ep*n8+2*ep*n9+4*ep*n10+4
      *ep*n11-18*ep+n3*n5+n3*n8+n3*n10+2*n3*n11-5*n3+n4*n5+n4*n8+n4*n10+2*n4*n11-5*n4+n5*n9+n5*n10-4*n5+n8*n9+n8*n10-4*
      n8+n9*n10+2*n9*n11-5*n9+n10^2+2*n10*n11-9*n10-8*n11+20,2*ep+2*n5+n8+n10+n11-5)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,n9,n10,n11,(n12+1),0,0)*rat(n8,1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,-1+n11,(n12+1),0,0)*rat(-n8,1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,1+n10,n11,(n12+1),0,0)*rat(n10,1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,-1+n11,(n12+1),0,0)*rat(-n10,1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(-16*ep^2-4*ep*n1-4*ep*n2-4*ep*n3-2*ep*n4-12*ep*n5-2*ep*n6
      -4*ep*n7-12*ep*n8-4*ep*n9-12*ep*n10-20*ep*n11-2*ep*(n12+1)+76*ep-2*n1*n5-2*n1*n8-2*n1*n10-4*n1*n11+10*n1-2*n2*n5-2*n2
      *n8-2*n2*n10-4*n2*n11+10*n2-2*n3*n5-2*n3*n8-2*n3*n10-4*n3*n11+10*n3-n4*n5-n4*n8-n4*n10-2*n4*n11+5*n4-2*n5^2-n5*n6
      -2*n5*n7-4*n5*n8-2*n5*n9-4*n5*n10-6*n5*n11-n5*(n12+1)+28*n5-n6*n8-n6*n10-2*n6*n11+5*n6-2*n7*n8-2*n7*n10-4*n7*n11+10*
      n7-2*n8^2-2*n8*n9-4*n8*n10-6*n8*n11-n8*(n12+1)+28*n8-2*n9*n10-4*n9*n11+10*n9-2*n10^2-6*n10*n11-n10*(n12+1)+28*n10-4*n11^2
      -2*n11*(n12+1)+46*n11+5*(n12+1)-90,2*ep+2*n5+n8+n10+n11-5)
      +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,(n12+1),0,0)*rat(4*ep^2+2*ep*n1+2*ep*n5+2*ep*n6+2*ep*n7+4*ep*n8+2*ep*n10+4
      *ep*n11-18*ep+n1*n5+n1*n8+n1*n10+2*n1*n11-5*n1+n5*n6+n5*n7+n5*n8-4*n5+n6*n8+n6*n10+2*n6*n11-5*n6+n7*n8+n7*n10+2*
      n7*n11-5*n7+n8^2+n8*n10+2*n8*n11-9*n8-4*n10-8*n11+20,2*ep+2*n5+n8+n10+n11-5)
      +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,(n12+1),0,0)*rat(-2*ep-n1-n6-n7-n8+4,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,(n12+1),0,0)*rat(-2*ep-n3-n4-n9-n10+4,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,(n12+1),0,0)*rat(-4*ep-2*n2-2*n5-n7-n8-n9-n10-2*n11-(n12+1)+10,1)
      );
*--#] n12 : 
#endif
#if ( `il' >= 5 )
*--#[ n1 :
`SKIP5'
id,ifmatch->haha1,`Z'(n1?{>1},n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`POS9',n10?`POS10',n11?`POS11',0,0,0)
			 = -rat(n6,-2*ep*(n1-1)-2*(n1-1)^2-2*(n1-1)*n7+2*(n1-1))*(
      +`Z'(1+(n1-1),-1+n2,n3,n4,n5,n6,1+n7,n8,n9,n10,n11,0,0,0)*rat(2*(n1-1)*n7,n6)
      +`Z'((n1-1),n2,n3,n4,-1+n5,1+n6,n7,1+n8,n9,n10,n11,0,0,0)*rat(-2*n8,1)
      +`Z'((n1-1),n2,n3,n4,n5,1+n6,n7,n8,n9,n10,n11,0,0,0)*rat(2*ep+2*n6+2*n8-2,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,1+n7,n8,n9,n10,n11,0,0,0)*rat(-2*ep*n7-2*(n1-1)*n7-2*n7^2+2*n7,n6)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,1+n8,n9,n10,n11,0,0,0)*rat(2*ep*n8+2*n6*n8+2*n8^2-2*n8,n6)
      );
*--#] n1 : 
#endif
#if ( `il' >= 6 )
*--#[ n2 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,n2?{>1},n3?`POS3',n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`POS9',n10?`POS10',n11?`POS11',0,0,0)
			 = -rat(n11,-2*ep*(n2-1)-2*(n2-1)^2-(n2-1)*n7-(n2-1)*n9+2*(n2-1))*(
      +`Z'(1,1+(n2-1),-1+n3,n4,n5,n6,n7,n8,1+n9,n10,n11,0,0,0)*rat((n2-1)*n9,n11)
      +`Z'(1,(n2-1),n3,n4,n5,n6,1+n7,-1+n8,n9,n10,1+n11,0,0,0)*rat(-n7,1)
      +`Z'(1,(n2-1),n3,n4,n5,n6,1+n7,n8,n9,n10,n11,0,0,0)*rat(-(n2-1)*n7+n7*n11,n11)
      +`Z'(1,(n2-1),n3,n4,n5,n6,n7,n8,1+n9,-1+n10,1+n11,0,0,0)*rat(-n9,1)
      +`Z'(1,(n2-1),n3,n4,n5,n6,n7,n8,1+n9,n10,n11,0,0,0)*rat(-(n2-1)*n9+n9*n11,n11)
      +`Z'(1,(n2-1),n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,0,0)*rat(2*ep+n7+n9+2*n11-2,1)
      +`Y'(0,1+(n2-1),n3,n4,n5,n6,1+n7,n8,n9,n10,n11,0,0,0)*rat((n2-1)*n7,n11)
      );
*--#] n2 : 
#endif
#if ( `il' >= 7 )
*--#[ n3 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,1,n3?{>1},n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`POS9',n10?`POS10',n11?`POS11',0,0,0)
			 = -rat(1,-2*ep-2*(n3-1)-2*n9+2)*(
      +`Z'(1,1,(n3-1),1+n4,-1+n5,n6,n7,n8,n9,1+n10,n11,0,0,0)*rat(-2*n4*n10,(n3-1))
      +`Z'(1,1,(n3-1),1+n4,n5,n6,n7,n8,n9,n10,n11,0,0,0)*rat(2*ep*n4+2*n4^2+2*n4*n10-2*n4,(n3-1))
      +`Z'(1,1,(n3-1),n4,n5,n6,n7,n8,1+n9,n10,n11,0,0,0)*rat(-2*ep*n9-2*(n3-1)*n9-2*n9^2+2*n9,(n3-1))
      +`Z'(1,1,(n3-1),n4,n5,n6,n7,n8,n9,1+n10,n11,0,0,0)*rat(2*ep*n10+2*n4*n10+2*n10^2-2*n10,(n3-1))
      +`Y'(1,0,1+(n3-1),n4,n5,n6,n7,n8,1+n9,n10,n11,0,0,0)*rat(2*n9,1)
      );
*--#] n3 : 
#endif
#if ( `il' >= 8 )
*--#[ n7 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,1,1,n4?`POS4',n5?`POS5',n6?`POS6',
	n7?{>1},n8?`POS8',n9?`POS9',n10?`POS10',n11?`POS11',0,0,0)
			 = -rat(2*ep+(n7-1)+n9,
       4*ep^2+6*ep*(n7-1)+2*ep*n9+2*ep*n11-2*ep+2*(n7-1)^2+2*(n7-1)*n9+2*(n7-1)*n11-2*(n7-1))*(
      +`Z'(1,1,1,n4,n5,n6,1+(n7-1),-1+n8,n9,n10,1+n11,0,0,0)*rat(-2*ep*n11-2*(n7-1)*n11,2*ep+(n7-1)+n9)
      +`Z'(1,1,1,n4,n5,n6,(n7-1),n8,1+n9,-1+n10,1+n11,0,0,0)*rat(2*ep*n9*n11+2*n9^2*n11,2*ep*(n7-1)+(n7-1)^2+(n7-1)*n9)
      +`Z'(1,1,1,n4,n5,n6,(n7-1),n8,1+n9,n10,n11,0,0,0)*rat(-4*ep^2*n9-2*ep*(n7-1)*n9-6*ep*n9^2-2*ep*n9*n11+2*ep*n9-2*(n7-1)*n9^2-
      2*n9^3-2*n9^2*n11+2*n9^2,2*ep*(n7-1)+(n7-1)^2+(n7-1)*n9)
      +`Z'(1,1,1,n4,n5,n6,(n7-1),n8,n9,n10,1+n11,0,0,0)*rat(4*ep*(n7-1)*n11-4*ep*n9*n11+2*(n7-1)^2*n11+2*(n7-1)*n11^2-2*(n7-1)*n11-2*n9^2*
      n11-2*n9*n11^2+2*n9*n11,2*ep*(n7-1)+(n7-1)^2+(n7-1)*n9)
      +`Y'(0,2,1,n4,n5,n6,1+(n7-1),n8,n9,n10,n11,0,0,0)*rat(-2*ep-2*n9,2*ep+(n7-1)+n9)
      +`Y'(1,2,0,n4,n5,n6,(n7-1),n8,1+n9,n10,n11,0,0,0)*rat(2*ep*n9+2*(n7-1)*n9,2*ep*(n7-1)+(n7-1)^2+(n7-1)*n9)
      );
*--#] n7 : 
#endif
#if ( `il' >= 9 )
*--#[ n8 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,1,1,n4?`POS4',n5?`POS5',n6?`POS6',
	1,n8?{>1},n9?`POS9',n10?`POS10',n11?`POS11',0,0,0)
			 = -rat(1,2*ep+n5+2*(n8-1)+n11-2)*(
      +`Z'(1,1,1,-1+n4,1+n5,n6,1,(n8-1),n9,1+n10,n11,0,0,0)*rat(n5*n10,(n8-1))
      +`Z'(1,1,1,n4,1+n5,-1+n6,1,1+(n8-1),n9,n10,n11,0,0,0)*rat(-n5,1)
      +`Z'(1,1,1,n4,1+n5,n6,1,(n8-1),n9,n10,n11,0,0,0)*rat(n5*(n8-1)-n5*n10,(n8-1))
      +`Z'(1,1,1,n4,n5,n6,1,(n8-1),-1+n9,1+n10,1+n11,0,0,0)*rat(n10*n11,(n8-1))
      +`Z'(1,1,1,n4,n5,n6,1,(n8-1),n9,1+n10,n11,0,0,0)*rat(-2*ep*n10-n5*n10-2*n10^2-n10*n11+2*n10,(n8-1))
      +`Z'(1,1,1,n4,n5,n6,1,(n8-1),n9,n10,1+n11,0,0,0)*rat((n8-1)*n11-n10*n11,(n8-1))
      +`Y'(1,1,1,n4,n5,n6,0,1+(n8-1),n9,n10,1+n11,0,0,0)*rat(-n11,1)
      );
*--#] n8 : 
#endif
#if ( `il' >= 10 )
*--#[ n9 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,1,1,n4?`POS4',n5?`POS5',n6?`POS6',
	1,1,n9?{>1},n10?`POS10',n11?`POS11',0,0,0)
			 = -rat(1,2*ep+n5+n10+n11-2)*(
      +`Z'(1,1,1,n4,1+n5,-1+n6,1,1,1+(n9-1),n10,n11,0,0,0)*rat(-n5,1)
      +`Z'(1,1,1,n4,n5,n6,1,1,(n9-1),1+n10,n11,0,0,0)*rat(-2*ep*n10-(n9-1)*n10-n10*n11+n10,(n9-1))
      +`Y'(0,2,1,n4,n5,n6,1,1,(n9-1),1+n10,n11,0,0,0)*rat(n10,(n9-1))
      +`Y'(1,1,1,n4,1+n5,n6,1,0,1+(n9-1),n10,n11,0,0,0)*rat(n5,1)
      +`Y'(1,1,1,n4,n5,n6,0,1,1+(n9-1),1+n10,n11,0,0,0)*rat(-2*ep*n10-(n9-1)*n10-n10*n11+n10,2*ep+(n9-1))
      +`Y'(1,1,1,n4,n5,n6,0,1,1+(n9-1),n10,1+n11,0,0,0)*rat(-2*ep*n11-(n9-1)*n11+n10*n11,2*ep+(n9-1))
      +`Y'(1,1,1,n4,n5,n6,0,1,(n9-1),1+n10,1+n11,0,0,0)*rat(-4*ep*n10*n11-2*(n9-1)*n10*n11-2*n10*n11^2+2*n10*n11,2*ep*(n9-1)+(n9-1)^2)
      +`Y'(1,1,1,n4,n5,n6,1,0,1+(n9-1),1+n10,n11,0,0,0)*rat(n10,1)
      +`Y'(1,1,1,n4,n5,n6,1,0,1+(n9-1),n10,1+n11,0,0,0)*rat(n11,1)
      +`Y'(1,1,1,n4,n5,n6,1,0,(n9-1),1+n10,1+n11,0,0,0)*rat(n10*n11,(n9-1))
      +`Y'(1,2,0,n4,n5,n6,0,1,1+(n9-1),1+n10,n11,0,0,0)*rat(-n10,2*ep+(n9-1))
      );
*--#] n9 : 
#endif
#if ( `il' >= 11 )
*--#[ n11 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,1,1,n4?`POS4',n5?`POS5',n6?`POS6',
	1,1,1,n10?`POS10',n11?{>1},0,0,0)
			 = -rat(2*ep+n5+(n11-1),
      4*ep^2+2*ep*n5+2*ep*n10+6*ep*(n11-1)-2*ep+2*n5*(n11-1)+2*n10*(n11-1)+2*(n11-1)^2-2*(n11-1))*(
      +`Z'(1,1,1,-1+n4,1+n5,n6,1,1,1,1+n10,(n11-1),0,0,0)*rat(2*ep*n5*n10+2*n5^2*n10,2*ep*(n11-1)+n5*(n11-1)+(n11-1)^2)
      +`Z'(1,1,1,n4,1+n5,-1+n6,1,2,1,n10,(n11-1),0,0,0)*rat(2*ep*n5+2*n5*(n11-1),2*ep*(n11-1)+n5*(n11-1)+(n11-1)^2)
      +`Z'(1,1,1,n4,1+n5,n6,1,1,1,n10,(n11-1),0,0,0)*rat(-4*ep^2*n5-6*ep*n5^2-2*ep*n5*n10-2*ep*n5*(n11-1)+2*ep*n5-2*n5^3-2*
      n5^2*n10-2*n5^2*(n11-1)+2*n5^2,2*ep*(n11-1)+n5*(n11-1)+(n11-1)^2)
      +`Z'(1,1,1,n4,n5,n6,1,1,1,1+n10,(n11-1),0,0,0)*rat(-4*ep*n5*n10+4*ep*n10*(n11-1)-2*n5^2*n10-2*n5*n10^2+2*n5*n10+2*n10^2*
      (n11-1)+2*n10*(n11-1)^2-2*n10*(n11-1),2*ep*(n11-1)+n5*(n11-1)+(n11-1)^2)
      +`Y'(1,1,1,n4,n5,n6,0,2,1,n10,1+(n11-1),0,0,0)*rat(-2*ep-2*n5,2*ep+n5+(n11-1))
      +`Y'(1,1,1,n4,n5,n6,1,1,0,1+n10,1+(n11-1),0,0,0)*rat(-2*ep*n10-2*n10*(n11-1),2*ep+n5+(n11-1))
      );
*--#] n11 : 
#endif
#if ( `il' >= 12 )
*--#[ n6 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,1,1,n4?`POS4',n5?`POS5',n6?{>1},
	1,1,1,n10?`POS10',1,0,0,0)
			 = -rat(n4,-2*ep*(n6-1)-n4*(n6-1)-(n6-1)*n10+(n6-1))*(
      +`Z'(1,1,1,1+n4,n5,(n6-1),1,1,1,n10,1,0,0,0)*rat(2*ep+(n6-1),1)
      +`Y'(1,0,2,n4,n5,1+(n6-1),1,1,1,n10,1,0,0,0)*rat((n6-1),n4)
      +`Y'(1,1,1,1+n4,-1+n5,1+(n6-1),0,2,1,n10,1,0,0,0)*rat(-(n6-1),ep)
      +`Y'(1,1,1,1+n4,-1+n5,1+(n6-1),1,1,0,1+n10,1,0,0,0)*rat((n6-1)*n10,ep)
      +`Y'(1,1,1,1+n4,1+n5,-1+(n6-1),1,2,1,n10,0,0,0,0)*rat(-n5,2*ep+n5)
      +`Y'(1,1,1,1+n4,1+n5,(n6-1),1,1,1,n10,0,0,0,0)*rat(-n5*n10+n5,2*ep+n5)
      +`Y'(1,1,1,1+n4,n5,1+(n6-1),0,1,1,n10,1,0,0,0)*rat(2*ep*(n6-1)+(n6-1)^2,ep)
      +`Y'(1,1,1,1+n4,n5,1+(n6-1),1,1,0,n10,1,0,0,0)*rat(-2*ep*(n6-1)-n4*(n6-1)-(n6-1)*n10+(n6-1),ep)
      +`Y'(1,1,1,1+n4,n5,(n6-1),0,2,1,n10,1,0,0,0)*rat(2*ep+(n6-1),ep)
      +`Y'(1,1,1,1+n4,n5,(n6-1),1,1,1,1+n10,0,0,0,0)*rat(-2*ep*n10-n5*n10-2*n10^2+2*n10,2*ep+n5)
      +`Y'(1,1,1,n4,1+n5,(n6-1),1,1,1,1+n10,0,0,0,0)*rat(n5*n10,2*ep+n5)
      +`Y'(1,1,1,n4,n5,1+(n6-1),1,1,0,1+n10,1,0,0,0)*rat(-2*ep*(n6-1)*n10-n4*(n6-1)*n10-(n6-1)*n10^2+(n6-1)*n10,ep*n4)
      +`Y'(1,1,1,n4,n5,1+(n6-1),1,1,1,1+n10,0,0,0,0)*rat((n6-1)*n10,n4)
      +`Y'(2,0,1,1+n4,n5,(n6-1),1,1,1,n10,1,0,0,0)*rat(-1,1)
      );
*--#] n6 : 
#endif
#if ( `il' >= 13 )
*--#[ n4 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,1,1,n4?{>1},n5?`POS5',1,1,1,1,n10?`POS10',1,0,0,0)
			 = -rat(n5,-2*ep*(n4-1)-(n4-1)*n5-(n4-1)*n10+(n4-1))*(
      +`Z'(1,1,1,(n4-1),1+n5,1,1,1,1,n10,1,0,0,0)*rat(2*ep+(n4-1)+n10-1,1)
      +`Y'(1,0,2,(n4-1),1+n5,1,1,1,1,n10,1,0,0,0)*rat(-1,1)
      +`Y'(1,1,1,1+(n4-1),1+n5,0,1,2,1,n10,0,0,0,0)*rat(-(n4-1),2*ep+n5)
      +`Y'(1,1,1,1+(n4-1),1+n5,1,1,1,0,n10,1,0,0,0)*rat(2*ep*(n4-1)+(n4-1)^2+(n4-1)*n10-(n4-1),ep)
      +`Y'(1,1,1,1+(n4-1),1+n5,1,1,1,1,n10,0,0,0,0)*rat(-2*ep*(n4-1)-(n4-1)*n5-(n4-1)*n10+(n4-1),2*ep+n5)
      +`Y'(1,1,1,1+(n4-1),n5,1,0,2,1,n10,1,0,0,0)*rat((n4-1),n5)
      +`Y'(1,1,1,1+(n4-1),n5,1,1,1,0,1+n10,1,0,0,0)*rat(ep*(n4-1)*n10-(n4-1)*n5*n10,ep*n5)
      +`Y'(1,1,1,1+(n4-1),n5,1,1,1,1,1+n10,0,0,0,0)*rat(-4*ep*(n4-1)*n10-2*(n4-1)*n5*n10-2*(n4-1)*n10^2+2*(n4-1)*n10,2*ep*n5+n5^2)
      +`Y'(1,1,1,(n4-1),1+n5,1,1,1,0,1+n10,1,0,0,0)*rat(2*ep*n10+(n4-1)*n10+n10^2-n10,ep)
      +`Y'(1,1,1,(n4-1),1+n5,1,1,1,1,1+n10,0,0,0,0)*rat(-2*ep*n10+(n4-1)*n10-n5*n10,2*ep+n5)
      );
*--#] n4 : 
#endif
#if ( `il' >= 14 )
*--#[ n10 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,1,1,1,n5?`POS5',1,1,1,1,n10?{>1},1,0,0,0)
			 = -rat(ep+1,6*ep^2+5*ep*(n10-1)+5*ep+(n10-1)^2+2*(n10-1)+1)*(
      +`Z'(1,1,1,1,1+n5,1,1,1,1,(n10-1),1,0,0,0)*rat(4*ep^2*n5+4*ep*n5*(n10-1)+2*ep*n5+n5*(n10-1)^2+n5*(n10-1),ep*(n10-1)+(n10-1))
      +`Y'(0,1,2,1,n5,1,1,1,1,1+(n10-1),1,0,0,0)*rat(-1,1)
      +`Y'(1,0,2,1,1+n5,0,1,1,2,(n10-1),1,0,0,0)*rat(-n5,ep*(n10-1)+(n10-1))
      +`Y'(1,0,2,1,1+n5,1,1,0,2,(n10-1),1,0,0,0)*rat(n5,ep*(n10-1)+(n10-1))
      +`Y'(1,0,2,1,1+n5,1,1,1,1,(n10-1),1,0,0,0)*rat(-2*ep*n5-n5*(n10-1)-n5,ep*(n10-1)+(n10-1))
      +`Y'(1,0,2,1,n5,1,0,1,2,(n10-1),2,0,0,0)*rat(-1,ep*(n10-1)+(n10-1))
      +`Y'(1,0,2,1,n5,1,1,0,2,1+(n10-1),1,0,0,0)*rat(1,ep+1)
      +`Y'(1,0,2,1,n5,1,1,0,2,(n10-1),2,0,0,0)*rat(1,ep*(n10-1)+(n10-1))
      +`Y'(1,0,2,1,n5,1,1,1,1,1+(n10-1),1,0,0,0)*rat(n5-1,ep+1)
      +`Y'(1,0,2,1,n5,1,1,1,2,(n10-1),1,0,0,0)*rat(2*ep+n5+(n10-1)-1,ep*(n10-1)+(n10-1))
      +`Y'(1,1,1,1,1+n5,0,1,1,1,1+(n10-1),1,0,0,0)*rat(-2*ep*n5-n5*(n10-1)-n5,ep+1)
      +`Y'(1,1,1,1,1+n5,1,1,0,1,1+(n10-1),1,0,0,0)*rat(2*ep*n5+n5*(n10-1)+n5,ep+1)
      +`Y'(1,1,1,1,1+n5,1,1,1,0,1+(n10-1),1,0,0,0)*rat(4*ep^2*n5+4*ep*n5*(n10-1)+2*ep*n5+n5*(n10-1)^2+n5*(n10-1),ep^2+ep)
      +`Y'(1,1,1,1,1+n5,1,1,1,1,1+(n10-1),0,0,0,0)*rat(-4*ep^2*n5-2*ep*n5^2-2*ep*n5*(n10-1)-n5^2*(n10-1)-n5^2+n5*(n10-1)+n5,2*ep^2+ep*
      n5+2*ep+n5)
      +`Y'(1,1,1,1,n5,1,0,1,1,1+(n10-1),2,0,0,0)*rat(-2*ep-(n10-1)-1,ep+1)
      +`Y'(1,1,1,1,n5,1,0,1,2,1+(n10-1),1,0,0,0)*rat(-1,1)
      +`Y'(1,1,1,1,n5,1,1,0,1,1+(n10-1),2,0,0,0)*rat(2*ep+(n10-1)+1,ep+1)
      +`Y'(1,1,1,1,n5,1,1,0,1,2+(n10-1),1,0,0,0)*rat(ep*(n10-1)+ep+(n10-1)^2+2*(n10-1)+1,ep+1)
      +`Y'(1,1,1,1,n5,1,1,1,0,2+(n10-1),1,0,0,0)*rat(-2*ep*n5*(n10-1)-2*ep*n5+2*ep*(n10-1)+2*ep-n5*(n10-1)^2-2*n5*(n10-1)-n5+(n10-1)^2+2*(n10-1)+1
      ,ep^2+ep)
      +`Y'(1,1,1,1,n5,1,1,1,1,2+(n10-1),0,0,0,0)*rat(2*ep*n5*(n10-1)+2*ep*n5-2*ep*(n10-1)-2*ep+n5^2*(n10-1)+n5^2-3*n5*(n10-1)-3*n5+2*(n10-1)+2
      ,2*ep^2+ep*n5+ep+n5-1)
      +`Y'(1,1,1,2,-1+n5,1,0,1,1,1+(n10-1),2,0,0,0)*rat(1,ep+1)
      +`Y'(1,1,1,2,-1+n5,1,0,2,1,1+(n10-1),1,0,0,0)*rat(-1,ep+1)
      +`Y'(1,1,1,2,-1+n5,1,1,0,1,1+(n10-1),2,0,0,0)*rat(-1,ep+1)
      +`Y'(1,1,1,2,-1+n5,1,1,0,1,2+(n10-1),1,0,0,0)*rat(-(n10-1)-1,ep+1)
      +`Y'(1,1,1,2,-1+n5,1,1,1,0,2+(n10-1),1,0,0,0)*rat(-ep*(n10-1)-ep+n5*(n10-1)+n5-(n10-1)-1,ep^2+ep)
      +`Y'(1,1,1,2,-1+n5,1,1,1,1,2+(n10-1),0,0,0,0)*rat(4*ep*(n10-1)+4*ep+2*n5*(n10-1)+2*n5+2*(n10-1)^2-2,2*ep^2+ep*n5+ep+n5-1)
      +`Y'(1,1,1,2,1+n5,0,1,1,1,(n10-1),1,0,0,0)*rat(-2*ep*n5-n5*(n10-1)-n5,ep*(n10-1)+(n10-1))
      +`Y'(1,1,1,2,1+n5,0,1,2,1,(n10-1),0,0,0,0)*rat(-2*ep*n5-n5*(n10-1)-n5,2*ep^2*(n10-1)+ep*n5*(n10-1)+2*ep*(n10-1)+n5*(n10-1))
      +`Y'(1,1,1,2,1+n5,1,1,0,1,(n10-1),1,0,0,0)*rat(2*ep*n5+n5*(n10-1)+n5,ep*(n10-1)+(n10-1))
      +`Y'(1,1,1,2,1+n5,1,1,1,0,(n10-1),1,0,0,0)*rat(4*ep^2*n5+4*ep*n5*(n10-1)+2*ep*n5+n5*(n10-1)^2+n5*(n10-1),ep^2*(n10-1)+ep*(n10-1))
      +`Y'(1,1,1,2,1+n5,1,1,1,1,(n10-1),0,0,0,0)*rat(-4*ep^2*n5-2*ep*n5^2-4*ep*n5*(n10-1)-n5^2*(n10-1)-n5^2-n5*(n10-1)^2+n5,2*ep^2*(n10-1)
      +ep*n5*(n10-1)+2*ep*(n10-1)+n5*(n10-1))
      +`Y'(1,1,1,2,n5,0,1,1,1,1+(n10-1),1,0,0,0)*rat(-ep+n5-1,ep+1)
      +`Y'(1,1,1,2,n5,0,1,2,1,1+(n10-1),0,0,0,0)*rat(n5-1,2*ep^2+ep*n5+ep+n5-1)
      +`Y'(1,1,1,2,n5,1,0,1,1,(n10-1),2,0,0,0)*rat(-2*ep-(n10-1)-1,ep*(n10-1)+(n10-1))
      +`Y'(1,1,1,2,n5,1,0,2,1,(n10-1),1,0,0,0)*rat(2*ep+(n10-1)+1,ep*(n10-1)+(n10-1))
      +`Y'(1,1,1,2,n5,1,1,0,1,1+(n10-1),1,0,0,0)*rat(2*ep-n5+(n10-1)+2,ep+1)
      +`Y'(1,1,1,2,n5,1,1,0,1,(n10-1),2,0,0,0)*rat(2*ep+(n10-1)+1,ep*(n10-1)+(n10-1))
      +`Y'(1,1,1,2,n5,1,1,1,0,1+(n10-1),1,0,0,0)*rat(2*ep^2-4*ep*n5+ep*(n10-1)+3*ep-2*n5*(n10-1)-2*n5+(n10-1)+1,ep^2+ep)
      +`Y'(1,1,1,2,n5,1,1,1,1,1+(n10-1),0,0,0,0)*rat(-16*ep^3-12*ep^2*n5-16*ep^2*(n10-1)+4*ep^2-10*ep*n5*(n10-1)-6*ep*n5-4*ep*
      (n10-1)^2+6*ep*(n10-1)+6*ep+n5^3-n5^2*(n10-1)-4*n5^2-2*n5*(n10-1)^2+n5*(n10-1)+5*n5+2*(n10-1)^2-2,4*ep^3+4*ep^2*n5+2*ep^2+ep*n5^2+3*ep*n5
      -2*ep+n5^2-n5)
      );
*--#] n10 : 
#endif
#if ( `il' >= 15 )
*--#[ n5 :
`SKIP5'
id,ifmatch->haha1,`Z'(1,1,1,1,n5?{>1},1,1,1,1,1,1,0,0,0)
			 = -rat(1,-6*ep^2*(n5-1)-5*ep*(n5-1)^2-2*ep*(n5-1)-(n5-1)^3-(n5-1)^2)*(
      +`Z'(1,1,1,1,(n5-1),1,1,1,1,1,1,0,0,0)*rat(20*ep^3+24*ep^2*(n5-1)+4*ep^2+9*ep*(n5-1)^2+4*ep*(n5-1)+(n5-1)^3+(n5-1)^2,1)
      +`Y'(0,1,1,1,1+(n5-1),1,1,1,1,1,1,0,0,0)*rat(ep*(n5-1)^2,1)
      +`Y'(0,1,1,1,1+(n5-1),1,2,1,1,1,0,0,0,0)*rat(ep*(n5-1)+(n5-1)^2,1)
      +`Y'(0,1,1,1,(n5-1),2,1,1,1,1,1,0,0,0)*rat(-2*ep*(n5-1)-(n5-1)^2,2)
      +`Y'(1,0,1,1,1+(n5-1),1,1,1,1,1,1,0,0,0)*rat(6*ep^2*(n5-1)+5*ep*(n5-1)^2+2*ep*(n5-1)+(n5-1)^3+(n5-1)^2,1)
      +`Y'(1,0,1,1,1+(n5-1),1,1,1,2,1,0,0,0,0)*rat(-ep*(n5-1)-(n5-1)^2,1)
      +`Y'(1,0,1,1,1+(n5-1),1,2,1,1,1,0,0,0,0)*rat(-ep*(n5-1)-(n5-1)^2,1)
      +`Y'(1,0,2,0,1+(n5-1),1,1,1,1,1,1,0,0,0)*rat(2*ep*(n5-1)+(n5-1)^2,2)
      +`Y'(1,0,2,1,1+(n5-1),1,1,1,1,1,1,0,0,0)*rat(-(n5-1)^2,2)
      +`Y'(1,0,2,1,(n5-1),2,1,1,1,1,1,0,0,0)*rat(2*ep*(n5-1)+(n5-1)^2,4*ep+2)
      +`Y'(1,1,0,1,1+(n5-1),1,1,1,1,1,1,0,0,0)*rat(ep*(n5-1)^2,1)
      +`Y'(1,1,0,1,1+(n5-1),1,1,1,2,1,0,0,0,0)*rat(ep*(n5-1)+(n5-1)^2,1)
      +`Y'(1,1,0,2,(n5-1),1,1,1,1,1,1,0,0,0)*rat(-2*ep*(n5-1)-(n5-1)^2,2)
      +`Y'(1,1,1,0,1+(n5-1),1,1,1,1,1,1,0,0,0)*rat(2*ep^2*(n5-1)+ep*(n5-1)^2,1)
      +`Y'(1,1,1,1,1+(n5-1),0,1,1,1,1,1,0,0,0)*rat(2*ep^2*(n5-1)+ep*(n5-1)^2,1)
      +`Y'(1,1,1,1,1+(n5-1),1,0,1,1,1,1,0,0,0)*rat(-2*ep^2*(n5-1)-2*ep*(n5-1)^2,1)
      +`Y'(1,1,1,1,1+(n5-1),1,1,1,0,1,1,0,0,0)*rat(-2*ep^2*(n5-1)-2*ep*(n5-1)^2,1)
      +`Y'(1,1,1,1,1+(n5-1),1,1,1,0,2,1,0,0,0)*rat(2*ep*(n5-1)^2+(n5-1)^2,ep)
      +`Y'(1,1,1,1,1+(n5-1),1,1,1,1,1,0,0,0,0)*rat(6*ep^2*(n5-1)+5*ep*(n5-1)^2+2*ep*(n5-1)+2*(n5-1)^2,1)
      +`Y'(1,1,1,1,1+(n5-1),1,1,1,1,2,0,0,0,0)*rat(-8*ep^2*(n5-1)^2-4*ep*(n5-1)^3+2*ep*(n5-1)^2-(n5-1)^3+2*(n5-1)^2,8*ep^2+4*ep*(n5-1)+4*ep+2*(n5-1))
      +`Y'(1,1,1,1,(n5-1),1,0,2,1,1,1,0,0,0)*rat(-10*ep^2-10*ep*(n5-1)-2*ep-(n5-1)^2-2*(n5-1),1)
      +`Y'(1,1,1,1,(n5-1),1,1,1,0,2,1,0,0,0)*rat(-10*ep^2-10*ep*(n5-1)-2*ep-(n5-1)^2-2*(n5-1),1)
      +`Y'(1,1,1,1,(n5-1),1,1,1,1,2,0,0,0,0)*rat(10*ep^2+10*ep*(n5-1)+2*ep+(n5-1)^2+2*(n5-1),1)
      +`Y'(1,1,1,1,(n5-1),1,1,2,1,1,0,0,0,0)*rat(10*ep^2+10*ep*(n5-1)+2*ep+(n5-1)^2+2*(n5-1),1)
      +`Y'(1,1,1,1,(n5-1),2,1,1,0,2,1,0,0,0)*rat(-2*ep*(n5-1)-(n5-1)^2,2*ep)
      +`Y'(1,1,1,1,(n5-1),2,1,1,1,2,0,0,0,0)*rat(2*ep*(n5-1)+(n5-1)^2,4*ep+2)
      +`Y'(1,1,1,2,-1+(n5-1),2,0,2,1,1,1,0,0,0)*rat(-2*ep*(n5-1)-(n5-1)^2,4*ep^2+2*ep)
      +`Y'(1,1,1,2,-1+(n5-1),2,1,1,0,2,1,0,0,0)*rat(2*ep*(n5-1)+(n5-1)^2,4*ep^2+2*ep)
      +`Y'(1,1,1,2,1+(n5-1),0,1,2,1,1,0,0,0,0)*rat(-6*ep*(n5-1)^2-(n5-1)^3-2*(n5-1)^2,8*ep^2+4*ep*(n5-1)+4*ep+2*(n5-1))
      +`Y'(1,1,1,2,1+(n5-1),1,1,1,0,1,1,0,0,0)*rat(2*ep*(n5-1)^2+(n5-1)^2,ep)
      +`Y'(1,1,1,2,1+(n5-1),1,1,1,1,1,0,0,0,0)*rat(-(n5-1)^2,1)
      +`Y'(1,1,1,2,(n5-1),1,0,2,1,1,1,0,0,0)*rat(4*ep*(n5-1)+(n5-1)^2,2*ep)
      +`Y'(1,1,1,2,(n5-1),1,1,1,0,2,1,0,0,0)*rat(ep*(n5-1)-(n5-1)^2,ep)
      +`Y'(1,1,1,2,(n5-1),1,1,1,1,2,0,0,0,0)*rat(-10*ep*(n5-1)-(n5-1)^2-4*(n5-1),4*ep+2)
      +`Y'(1,1,1,2,(n5-1),2,0,1,1,1,1,0,0,0)*rat(2*ep*(n5-1)+(n5-1)^2,2*ep)
      +`Y'(1,1,1,2,(n5-1),2,1,1,0,1,1,0,0,0)*rat(-2*ep*(n5-1)-(n5-1)^2,2*ep)
      +`Y'(1,1,2,1,1+(n5-1),1,1,1,0,1,1,0,0,0)*rat(-ep*(n5-1)-(n5-1)^2,1)
      +`Y'(2,0,1,1,1+(n5-1),0,1,1,1,1,1,0,0,0)*rat(2*ep*(n5-1)+(n5-1)^2,2)
      +`Y'(2,0,1,1,1+(n5-1),1,1,1,1,1,1,0,0,0)*rat((n5-1)^2,2)
      +`Y'(2,0,1,2,(n5-1),1,1,1,1,1,1,0,0,0)*rat(-2*ep*(n5-1)-(n5-1)^2,4*ep+2)
      +`Y'(2,1,1,1,1+(n5-1),1,0,1,1,1,1,0,0,0)*rat(-ep*(n5-1)-(n5-1)^2,1)
      );
*--#] n5 : 
#endif
Goto haha2;
Label haha1;
$ired = 1;
Label haha2;
id	`Z'(?a,x?neg0_,?b,n12?,n13?,n14?) = `Y'(?a,x,?b,n12,n13,n14);
*id	`Y'(?a) = 0;
ModuleOption,maximum,$ired;
*Print +f +s;
.sort:Main haha loop, pass `$pass';
#if ( `$ired' > 0 )
	#redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
id	`Z'(1,1,1,1,1,1,1,1,1,1,1,0,0,0) = Master(haha);
*
#break
* 	#] haha : 
*--#[ no1 :
#case no1
*
* 	#[ Loop 1 :
#do il = 1,4
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering loop 1-`il', pass `$pass' at `time_' sec."
#endif
#if ( `il' >= 1 )
*--#[ stat file :
*
*  Prework to make the reduction of the no1 topology faster.
*
* B3
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n8)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n8+1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,
   n9,n10,n11,n12,1+n13,n14)*rat(-n1,1)+`Z'(1+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(n1,1)+`Z'(n1,-1+
   n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)
   *rat(-n8+1,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,n12,1+n13,n14)*rat(n9,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,
   n10,n11,n12,1+n13,n14)*rat(-2*ep-n1-2*n2-n9-n12-n13+3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-
   1+n8,1));
* B13
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?{>1},n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n10)*(`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,n14)*rat(-n10+1,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,
   n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n10+1,1)+`Z'(n1,n2,n3,n4,n5,n6,1+n7,-1+n8,n9,-1+n10,n11,1+n12,
   n13,n14)*rat(-n7,1)+`Z'(n1,n2,n3,n4,n5,n6,1+n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-
   1+n9,n10,n11,1+n12,n13,n14)*rat(-1+n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,-1+n13,1+n14)*rat(-n14,1)+`Z'(
   n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(2*ep+n3+n7+2*n9+n10+n14-5,1));
* B14
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n7)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n7+1,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,1+
   n9,n10,n11,n12,1+n13,n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,n2,-
   1+n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)
   *rat(-n7+1,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,
   n10,n11,n12,1+n13,n14)*rat(-2*ep-2*n3-n7-n9-n10-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-
   1+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-1+n7,1));
* B15
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n9)*(`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n9+1,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+
   n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,n2,
   n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,
   n14)*rat(-n9+1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,1+n7,
   n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n7,1)+`Z'(n1,n2,n3,n4,n5,n6,1+n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(
   n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,
   n14)*rat(2*ep+n3+n7+n9+2*n10+n14-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(-1+n9,1));
* B12
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n3)*(`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,1+n13,n14)*rat(-n3+1,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n10,1)+`Z'(n1,
   n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+
   n13,n14)*rat(-n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,n12,1+n13,n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,
   -1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(
   n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-2*ep-n3-2*n7-n9-n10-n14+5,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n3+1,1)+`Z'(n1,n2,n3,
   n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-
   1+n3,1));
* B17
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n4)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,1+n13,n14)*rat(-n4+1,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(n1,
   n2,n3,-1+n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n5,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,-1+n11,n12,1+
   n13,n14)*rat(n5,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+
   n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,1+n13,n14)*rat(n12,1)+
   `Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,
   n14)*rat(-2*ep-n4-n5-2*n6-n10-n12+5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,
   n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n4,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*
   rat(-1+n4,1));
* B20
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,
   1+n10,n11,n12,1+n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n5,1)+`Z'(n1,n2,
   n3,-1+n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,
   n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,
   n10,n11,n12,1+n13,n14)*rat(-2*ep-2*n4-n5-n6-n10-n12+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-
   1+n6,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-1+n6,1));
* B8
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,
   n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,
   n3,-1+n4,1+n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)
   *rat(-n6+1,1)+`Z'(n1,n2,n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,1+n7,-1+n8,
   n9,n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,1+n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n7,1)+`Z'(n1,n2,
   n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-6*ep-2*n2-2*n3-2*n4-n5-n6-n7-2*n9-2*n10-n11-2*n12-n13-n14+11,1)+
   `Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n6,1));
* B19
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n5)*(`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,
   n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5+1,1)+`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,
   1+n14)*rat(n4,1)+`Z'(n1,n2,n3,n4,-1+n5,1+n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n6,1)+`Z'(n1,n2,n3,n4,-1+n5,1+n6,n7,
   n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(n6,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(
   n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(2*ep+n4+n5+n6+2*n10+n12-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-
   1+n10,n11,n12,n13,1+n14)*rat(-1+n5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5+1,1));
* B10
id,ifmatch->no11,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,n14)*rat(-n1,1)+`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+
   n2,1+n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,
   n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,
   n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(-2*n12-2,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,
   1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,
   1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+`Z'(n1,-1+n2,n3,n4,
   n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-
   n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(6*ep+n1+2*n2+n3+2*n6+2*n7+2*n8+2*n9+2*n10+n11+n12+
   n13+n14-13,1));
* B6
id,ifmatch->no11,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,2+n13,n14)*rat(1+n13,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n2,1)+`Z'(1+n1,-1+
   n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*
   rat(-n1,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,
   n9,n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,n12,1+n13,n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,1+
   n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,2+n13,n14)*rat(-n13-1,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,
   n7,n8,n9,n10,-1+n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+
   `Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(2*ep+n1+n2+2*n8+n9+n12-5,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,2+n13,n14)*rat(-n13-1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-n2+1,1));
* B11
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n3+1)*(`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,
   1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(n2,1)+`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,
   n13,1+n14)*rat(n4,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,
   n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+
   `Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,n12,n13,1+n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,
   n12,n13,1+n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,
   n6,n7,-1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,1+n14)*rat(
   n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,
   n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(6*ep+n2+2*n3+n4+2*n5+2*n6+2*
   n7+n9+2*n10+2*n11+n12+n13+n14-13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(-2*n14-2,1));
* B16
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?{>1},n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n10)*(`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n10+1,1)+`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,
   n9,-1+n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(n1,n2,n3,1+n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,
   n3,n4,-1+n5,1+n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(n6,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*
   rat(-1+n10,1)+`Z'(n1,n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,-1+n11,n12,n13,1+n14)*rat(-n6,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+
   n10,n11,1+n12,-1+n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(2*ep+n4+2*n5+n6+
   n10+n12-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n10+1,1));
* B2
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n7+1)*(`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,
   1+n9,n10,n11,n12,n13,1+n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,n12,n13,1+n14)*rat(n9,1)+`Z'(n1,n2,
   -1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n7,1)+`Z'(n1,n2,n3,n4,-1+n5,1+n6,-1+n7,n8,n9,n10,n11,n12,n13,1+
   n14)*rat(n6,1)+`Z'(n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-n6,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,
   n9,n10,n11,n12,n13,1+n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-4*ep-2*n1-2*n2-n6-n7-
   n8-n9-2*n11-n12-n13+9,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n7,1));
* B9
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n11)*(`Z'(1+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(n1,1)+`Z'(1+n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,-1+n11,n12,n13,1+n14)*rat(-n1,1)+`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(n1,1+
   n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,-1+n11,n12,n13,1+n14)*rat(n2,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,
   1+n14)*rat(-n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n11+1,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,-1+n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(n6,1)+`Z'(
   n1,n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,-1+n11,n12,n13,1+n14)*rat(-n6,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,
   1+n14)*rat(-1+n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,-1+n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,
   n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(4*ep+n1+n2+n6+
   2*n7+2*n8+2*n9+n11+n12+n13-9,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n11+1,1));
* B1
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?{>1},n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n11)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n11+1,1)+`Z'(1+n1,n2,n3,n4,n5,n6,n7,-1+n8,
   n9,n10,-1+n11,n12,1+n13,n14)*rat(n1,1)+`Z'(1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-n1,1)+`Z'(n1,n2,
   n3,-1+n4,1+n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(n5,1)+`Z'(n1,n2,n3,n4,1+n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,
   n14)*rat(-n5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,-1+n11,n12,1+n13,n14)*rat(8*ep+n1+2*n2+2*n3+2*n4+n5+2*n6+2*n7+2*n8+2*n9+2*n10+2*n11+2*n12+n13+n14-17,1)+`Z'(n1,n2,
   n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n11+1,1));
* B4
id,ifmatch->no11,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n1+1)*(`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,-
   1+n9,n10,n11,n12,n13,1+n14)*rat(n2,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(-1+
   n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9,n10,n11,n12,
   n13,1+n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,
   n7,n8,n9,-1+n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(
   -1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(2*ep+n2+n8+2*n9+n12+n13-4,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n1,1)+`Z'(n1,n2,
   n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,-1+n13,1+n14)*
   rat(-1+n1,1));
* B7
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?pos_,n11?pos_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n9)*(`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,-
   1+n9,n10,n11,n12,n13,1+n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(-1+
   n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,-1+n11,n12,
   n13,1+n14)*rat(n2,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,
   n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n9,1)+
   `Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,
   n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n9+1,1)+`Z'(n1,n2,n3,n4,n5,n6,
   n7,1+n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,-1+n11,1+n12,n13,1+n14)*rat(n12,1)
   +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,-1+n13,1+n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+
   n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,
   n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-2*ep-2*n1-n2-n8-n9-n12+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,-1+n13,
   1+n14)*rat(-n9+1,1));
* B18
id,ifmatch->no11,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n4+1)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,
   n8,n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)+`Z'(n1,
   -1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(-2*n12-2,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+
   n12,n13,n14)*rat(n3,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,n14)*rat(n5,1)+`Z'(n1,n2,n3,-1+n4,1+n5,
   n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,1+n14)*rat(-n14,1)+
   `Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*ep-2*n2-n3-n11-n12-n13+3,1));
* B5
id,ifmatch->no11,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n1+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,1+
   n8,n9,n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+`Z'(-1+
   n1,n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,1+n12,
   n13,n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+
   n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(-1+
   n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*ep-n1-2*n6-n7-n8-n11+5,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,
   n10,n11,1+n12,n13,n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,
   n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,n14)*rat(-
   1+n1,1));
*--#] stat file : 
#endif
#if ( `il' >= 2 )
*--#[ n12 :
id,ifmatch->no11,
	`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_
			,n8?pos_,n9?pos_,n10?pos_,n11?pos_,n12?neg_,n13?neg0_,n14?neg0_) =
         -RAT(4*ep+3*n4+2*n5+2*n6+2*n10+2*n12-7,1)*(
*      +U(-1)*(
         +`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,n14)*rat(-n5,1)
         +`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(n5,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,n14)*rat(-6*ep-2*n2-2*n3-2*n4-n6-n7-2*n9-2*n10-n11-3*n12-
         n13+9,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,-1+n14)*rat(-2*ep-n3-n7-2*n9-n10-n14+5,1)
*         )
*      +U(0)*(
         +`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n1,1)
         +`Z'(n1,n2,n3,n4,1+n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(2*n5,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,1+n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n8,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(4*ep+n1+2*n2+2*n3+n4-n5+n7+n8+2*n9+n10+n11+n12+n13+
         n14-9,1)
*         )
*      +U(1)*(
         +`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(n2,1)
         +`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(n8,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*ep-n1-n2-n6-n7-3*n8+n11+5,1)
         +`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n1,1)
         +`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n3,1)
         +`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,1+n12,n13,n14)*rat(n9,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,1+n12,n13,n14)*rat(-3*n9,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(10*ep+n1+3*n2+2*n3+2*n6+4*n7+2*n8+5*n9+2*n10+2*n11
         +n12+2*n13+2*n14-20,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,2+n12,-1+n13,n14)*rat(n12+1,1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,-1+n11,1+n12,n13,n14)*rat(-n2,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,1+n12,n13,n14)*rat(2*n9,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,1+n14)*rat(-2*n14,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(4*ep+2*n3+2*n7+2*n9+2*n10+2*n14-8,1)
         +`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)
         +`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(-2*n5,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,1+n14)*rat(n14,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(2*n2+n3-n4-n5-n6-n10+n11+n13,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,2+n12,-1+n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(2*ep+n4+n5+n6+n10+n12-3,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n9+1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-6*ep-3*n1-4*n2-n6-n7-n8-n9-3*n11-2*n12-2*n13+11,1
         )
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-2*ep+2*n1+n2-2*n3-2*n7+n8-n9-2*n10+n12-2*n14+4,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,2+n12,-1+n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(-12*ep-3*n2-4*n3-4*n4-2*n5-2*n6-2*n7-4*n9-4*n10-2*
         n11-4*n12-2*n13-2*n14+21,1)
*         )
*      +U(2)*(
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-2*n14,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(2*n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,2+n12,n13,n14)*rat(n12+1,1)
*         )
*      +U(3)*(
         +`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)
         +`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,-1+n11,2+n12,n13,n14)*rat(-n12-1,1)
*         )
      );
*--#] n12 : 
#endif
#if ( `il' >= 3 )
*--#[ n14 :
id,ifmatch->no11,
	`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_
			,n8?pos_,n9?pos_,n10?pos_,n11?pos_,0,n13?neg0_,n14?neg_) =
			-RAT(-2*ep-2*n3-n7-n9-n10-n14+4,1)*(
*      +U(-1)*(
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,-1+n13,1+n14)*rat(-4*ep-2*n1-2*n2-n6-n7-n8-n9-2*n11-n13+9,1)
*         )
*      +U(0)*(
         +`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(n1,1)
         +`Z'(1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,n13,1+n14)*rat(-n1,1)
         +`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(-n4,1)
         +`Z'(n1,n2,n3,n4,1+n5,n6,n7,n8,n9,n10,-1+n11,0,n13,1+n14)*rat(n5,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(4*ep+n1+2*n2+2*n3+n4-n5+n7+n8+2*n9+n10+n13+n14-8,1)
*         )
*      +U(1)*(
         +`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(-n2,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(2*ep+n1+n2+n6+n7+2*n11-5,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,0,n13,1+n14)*rat(-n9,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(2*ep+n1+n2+n8+n9-4,1)
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,n13,1+n14)*rat(n2,1)
         +`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(n4,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,0,n13,1+n14)*rat(n9,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(6*ep+n2+2*n3+n4+2*n5+2*n6+2*n7+n9+2*n10+2*n11+n13+
         n14-13,1)
         +`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,-1+n11,0,n13,1+n14)*rat(-n5,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(-2*ep-n4-n5-n6-n10+4,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,0,n13,1+n14)*rat(2*ep+n4+n5+n6+n10-4,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,0,n13,1+n14)*rat(-4*ep-2*n1-2*n2-n6-n7-n8-n9-2*n11-n13+9,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,n13,1+n14)*rat(-8*ep-n1-3*n2-2*n3-2*n4-n5-2*n6-2*n7-n8-2*n9-2*n10
         -2*n11-n13-n14+17,1)
*         )
*      +U(2)*(
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,1+n14)*rat(n13,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,1+n14)*rat(-n13,1)
*         )
*      +U(3)*(
         +`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,1+n14)*rat(n13,1)
         +`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,1+n14)*rat(-n13,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,0,1+n13,1+n14)*rat(-n13,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,0,1+n13,1+n14)*rat(n13,1)
*         )
      );
*--#] n14 : 
#endif
#if ( `il' >= 4 )
*--#[ n13 :
id,ifnomatch->no12,
	`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_
		,n8?pos_,n9?pos_,n10?pos_,n11?pos_,0,n13?neg_,0) =
			-frat(20*ep^2+4*ep*n1+6*ep*n2+20*ep*n3+2*ep*n4+4*ep*n5+6*ep*n6+
         20*ep*n7+2*ep*n8+14*ep*n9+14*ep*n10+8*ep*n11+4*ep*n13-92*ep+2*n1*n3+4*n1*n7+2*n1*n9+2*n1*n10-10*n1+4*n2*n3+5*
         n2*n7+3*n2*n9+3*n2*n10-15*n2+4*n3^2+2*n3*n4+4*n3*n5+5*n3*n6+7*n3*n7+n3*n8+5*n3*n9+6*n3*n10+6*n3*n11+3*n3*n13-
         44*n3+n4*n7+n4*n9+n4*n10-5*n4+2*n5*n7+2*n5*n9+2*n5*n10-10*n5+4*n6*n7+3*n6*n9+3*n6*n10-15*n6+4*n7^2+2*n7*n8+6*
         n7*n9+5*n7*n10+6*n7*n11+3*n7*n13-44*n7+n8*n9+n8*n10-5*n8+2*n9^2+4*n9*n10+4*n9*n11+2*n9*n13-31*n9+2*n10^2+4*n10
         *n11+2*n10*n13-31*n10-20*n11-10*n13+105,2*ep+2*n3+n7+n9+n10-5)*
			ZZ(n1,...,n11,n13);
id	ZZ(n1?,...,n11?,n13?)*frat(x1?,x2?) = rat(x2,x1)*(
*      +U(0)*(
         +`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(-2*ep*n1-n1*n3-2*n1*n7-n1*n9-n1*n10+5*n1,2*ep+2*n3
         +n7+n9+n10-5)
         +`Z'(1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,1+n13,0)*rat(-n1*n3+n1*n7,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(2*ep*n4+n3*n4+2*n4*n7+n4*n9+n4*n10-5*n4,2*ep+2*n3+
         n7+n9+n10-5)
         +`Z'(n1,n2,n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(-n5,1)
         +`Z'(n1,n2,n3,n4,1+n5,n6,n7,n8,n9,n10,-1+n11,0,1+n13,0)*rat(n3*n5-n5*n7,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9,n10,n11,0,1+n13,0)*rat(n8,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(-12*ep^2-4*ep*n2-12*ep*n3-4*ep*n4-4*ep*n6-12*ep*n7-4*ep
         *n8-10*ep*n9-10*ep*n10-4*ep*n13+52*ep+n1*n3-n1*n7-2*n2*n3-4*n2*n7-2*n2*n9-2*n2*n10+10*n2-2*n3^2-3*n3*n4-n3*n5-
         4*n3*n6-5*n3*n7-3*n3*n8-4*n3*n9-5*n3*n10-3*n3*n13+24*n3-3*n4*n7-2*n4*n9-2*n4*n10+10*n4+n5*n7-2*n6*n7-2*n6*n9-2
         *n6*n10+10*n6-2*n7^2-3*n7*n8-5*n7*n9-4*n7*n10-3*n7*n13+24*n7-2*n8*n9-2*n8*n10+10*n8-2*n9^2-4*n9*n10-2*n9*n13+
         21*n9-2*n10^2-2*n10*n13+21*n10+10*n13-55,2*ep+2*n3+n7+n9+n10-5)
*         )
*      +U(1)*(
         +`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(2*ep*n2+n2*n3+2*n2*n7+n2*n9+n2*n10-5*n2,2*ep+2*
         n3+n7+n9+n10-5)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9,n10,n11,0,1+n13,0)*rat(-n8,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(-20*ep^2-4*ep*n1-6*ep*n2-22*ep*n3-2*ep*n4-4*ep*n5-6*
         ep*n6-18*ep*n7-2*ep*n8-14*ep*n9-14*ep*n10-8*ep*n11-4*ep*n13+92*ep-3*n1*n3-3*n1*n7-2*n1*n9-2*n1*n10+10*n1-5*n2*
         n3-4*n2*n7-3*n2*n9-3*n2*n10+15*n2-4*n3^2-2*n3*n4-4*n3*n5-5*n3*n6-7*n3*n7-2*n3*n8-6*n3*n9-6*n3*n10-6*n3*n11-4*
         n3*n13+47*n3-n4*n7-n4*n9-n4*n10+5*n4-2*n5*n7-2*n5*n9-2*n5*n10+10*n5-4*n6*n7-3*n6*n9-3*n6*n10+15*n6-4*n7^2-n7*
         n8-5*n7*n9-5*n7*n10-6*n7*n11-2*n7*n13+41*n7-n8*n9-n8*n10+5*n8-2*n9^2-4*n9*n10-4*n9*n11-2*n9*n13+31*n9-2*n10^2-
         4*n10*n11-2*n10*n13+31*n10+20*n11+10*n13-105,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(1+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,-1+n11,0,1+n13,0)*rat(n1,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,0,1+n13,0)*rat(2*ep*n9+n3*n9+2*n7*n9+n9^2+n9*n10-5*n9,2*ep+2*
         n3+n7+n9+n10-5)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(-4*ep^2-2*ep*n1-2*ep*n2-2*ep*n3-4*ep*n7-2*ep*n8-4*ep
         *n9-2*ep*n10+18*ep-n1*n3-2*n1*n7-n1*n9-n1*n10+5*n1-n2*n3-2*n2*n7-n2*n9-n2*n10+5*n2-n3*n8-n3*n9+4*n3-2*n7*n8-2*
         n7*n9+8*n7-n8*n9-n8*n10+5*n8-n9^2-n9*n10+9*n9+4*n10-20,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,1+n13,0)*rat(-2*ep*n2-n2*n3-2*n2*n7-n2*n9-n2*n10+5*n2,2*ep+2
         *n3+n7+n9+n10-5)
         +`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(-2*ep*n4-n3*n4-2*n4*n7-n4*n9-n4*n10+5*n4,2*ep+2
         *n3+n7+n9+n10-5)
         +`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,1+n10,n11,0,1+n13,0)*rat(n10,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,0,1+n13,0)*rat(-2*ep*n9-n3*n9-2*n7*n9-n9^2-n9*n10+5*n9,2*ep+2*
         n3+n7+n9+n10-5)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,0,1+n13,0)*rat(-n10,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(-12*ep^2-2*ep*n2-10*ep*n3-2*ep*n4-4*ep*n5-4*ep*n6-16
         *ep*n7-8*ep*n9-10*ep*n10-4*ep*n11-2*ep*n13+56*ep-n2*n3-2*n2*n7-n2*n9-n2*n10+5*n2-2*n3^2-n3*n4-2*n3*n5-2*n3*n6-
         6*n3*n7-3*n3*n9-4*n3*n10-2*n3*n11-n3*n13+23*n3-2*n4*n7-n4*n9-n4*n10+5*n4-4*n5*n7-2*n5*n9-2*n5*n10+10*n5-4*n6*
         n7-2*n6*n9-2*n6*n10+10*n6-4*n7^2-4*n7*n9-6*n7*n10-4*n7*n11-2*n7*n13+36*n7-n9^2-3*n9*n10-2*n9*n11-n9*n13+18*n9-
         2*n10^2-2*n10*n11-n10*n13+23*n10+10*n11+5*n13-65,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,-1+n11,0,1+n13,0)*rat(2*ep*n5+n3*n5+2*n5*n7+n5*n9+n5*n10-5*n5,2*ep+2*
         n3+n7+n9+n10-5)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(4*ep^2+2*ep*n3+2*ep*n4+2*ep*n5+2*ep*n6+4*ep*n7+2*ep*
         n9+4*ep*n10-18*ep+n3*n4+n3*n5+n3*n6+n3*n10-4*n3+2*n4*n7+n4*n9+n4*n10-5*n4+2*n5*n7+n5*n9+n5*n10-5*n5+2*n6*n7+n6
         *n9+n6*n10-5*n6+2*n7*n10-8*n7+n9*n10-4*n9+n10^2-9*n10+20,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,n2,n3,n4,-1+n5,1+n6,-1+n7,n8,n9,n10,n11,0,1+n13,0)*rat(n6,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(2*ep*n3-2*ep*n7+n3*n4+n3*n5+n3*n6+n3*n10-4*n3-n4*n7-
         n5*n7-n6*n7-n7*n10+4*n7,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,0,1+n13,0)*rat(-2*ep-n4-n5-n6-n10+4,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,-1+n11,0,1+n13,0)*rat(-n6,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,0,1+n13,0)*rat(-6*ep-2*n1-2*n2-2*n3-n6-2*n7-n8-2*n9-n10-2*n11-n13+
         13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,0,1+n13,0)*rat(24*ep^2+6*ep*n1+8*ep*n2+24*ep*n3+2*ep*n4+4*ep*n5+6*
         ep*n6+22*ep*n7+4*ep*n8+18*ep*n9+16*ep*n10+8*ep*n11+6*ep*n13-108*ep+4*n1*n3+5*n1*n7+3*n1*n9+3*n1*n10-15*n1+6*n2
         *n3+6*n2*n7+4*n2*n9+4*n2*n10-20*n2+4*n3^2+2*n3*n4+4*n3*n5+5*n3*n6+7*n3*n7+3*n3*n8+7*n3*n9+6*n3*n10+6*n3*n11+5*
         n3*n13-50*n3+n4*n7+n4*n9+n4*n10-5*n4+2*n5*n7+2*n5*n9+2*n5*n10-10*n5+4*n6*n7+3*n6*n9+3*n6*n10-15*n6+4*n7^2+3*n7
         *n8+7*n7*n9+5*n7*n10+6*n7*n11+4*n7*n13-47*n7+2*n8*n9+2*n8*n10-10*n8+3*n9^2+5*n9*n10+4*n9*n11+3*n9*n13-39*n9+2*
         n10^2+4*n10*n11+3*n10*n13-34*n10-20*n11-15*n13+120,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,1+n13,0)*rat(32*ep^2+4*ep*n1+10*ep*n2+32*ep*n3+6*ep*n4+4*ep*n5+10
         *ep*n6+32*ep*n7+6*ep*n8+24*ep*n9+24*ep*n10+8*ep*n11+4*ep*n13-148*ep+3*n1*n3+3*n1*n7+2*n1*n9+2*n1*n10-10*n1+7*
         n2*n3+8*n2*n7+5*n2*n9+5*n2*n10-25*n2+6*n3^2+4*n3*n4+3*n3*n5+8*n3*n6+12*n3*n7+5*n3*n8+10*n3*n9+10*n3*n10+6*n3*
         n11+3*n3*n13-71*n3+5*n4*n7+3*n4*n9+3*n4*n10-15*n4+3*n5*n7+2*n5*n9+2*n5*n10-10*n5+7*n6*n7+5*n6*n9+5*n6*n10-25*
         n6+6*n7^2+4*n7*n8+10*n7*n9+10*n7*n10+6*n7*n11+3*n7*n13-71*n7+3*n8*n9+3*n8*n10-15*n8+4*n9^2+8*n9*n10+4*n9*n11+2
         *n9*n13-54*n9+4*n10^2+4*n10*n11+2*n10*n13-54*n10-20*n11-10*n13+170,2*ep+2*n3+n7+n9+n10-5)
*         )
*      +U(2)*(
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,2+n13,0)*rat(2*n13+2,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,2+n13,0)*rat(-2*ep*n13-2*ep-n3*n13-n3-2*n7*n13-2*n7-n9*n13-n9-n10
         *n13-n10+5*n13+5,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,0,2+n13,0)*rat(2*ep*n13+2*ep+n3*n13+n3+2*n7*n13+2*n7+n9*n13+n9+n10*
         n13+n10-5*n13-5,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,2+n13,0)*rat(-2*n13-2,1)
*         )
*      +U(3)*(
         +`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,2+n13,0)*rat(-2*ep*n13-2*ep-n3*n13-n3-2*n7*n13-2*n7-n9*n13-n9-
         n10*n13-n10+5*n13+5,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,0,2+n13,0)*rat(2*ep*n13+2*ep+n3*n13+n3+2*n7*n13+2*n7+n9*n13+n9+
         n10*n13+n10-5*n13-5,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,0,2+n13,0)*rat(2*ep*n13+2*ep+n3*n13+n3+2*n7*n13+2*n7+n9*n13+n9+
         n10*n13+n10-5*n13-5,2*ep+2*n3+n7+n9+n10-5)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,0,2+n13,0)*rat(-2*ep*n13-2*ep-n3*n13-n3-2*n7*n13-2*n7-n9*n13-n9-
         n10*n13-n10+5*n13+5,2*ep+2*n3+n7+n9+n10-5)
*         )
      );
goto no11;
*--#] n13 : 
#endif
Goto no12;
Label no11;
$ired = 1;
Label no12;
id	`Z'(?a,x?neg0_,?b,n12?,n13?,n14?) = `Y'(?a,x,?b,n12,n13,n14);
ModuleOption,maximum,$ired;
*Print +f +s;
.sort:Main no1 loop, pass `$pass';
#if ( `$ired' > 0 )
	#redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
* 	#] Loop 1 : 
* 	#[ Loop 2 :
#do il = 1,11
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering loop 2-`il', pass `$pass' at `time_' sec."
#endif
*
#if ( `il' >= 1 )
*--#[ n1 :
id,ifmatch->no13
	`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_
			,n8?pos_,n9?pos_,n10?pos_,n11?pos_,0,0,0) =
			-RAT(-2*ep-2*n1-2*n2+4,1)*(
         +`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,0,0)*rat(2*ep*n2+2*n1*n2+2*n2^2-4*n2,-n1+1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,1+n9,n10,n11,0,0,0)*rat(2*n8*n9,-n1+1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9,n10,n11,0,0,0)*rat(-2*ep*n8-2*n8^2-2*n8*n9+2*n8,-n1+1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,0,0,0)*rat(-2*ep*n9-2*n8*n9-2*n9^2+2*n9,-n1+1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,0,0)*rat(2*n2,1)
      );
*--#] n1 : 
#endif
#if ( `il' >= 2 )
*--#[ n11 :
id,ifmatch->no13
	`Z'(1,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_
			,n8?pos_,n9?pos_,n10?pos_,n11?{>1},0,0,0) =
			-RAT(2*ep*n11-2*ep+n2*n11-n2+n6*n11-n6+2*n11^2-6*n11+4,n3)*(
         +`Z'(1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,-1+n11,0,0,0)*rat(-n2*n10,n3)
         +`Z'(1,1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,-1+n11,0,0,0)*rat(n2*n10,n3)
         +`Z'(1,1+n2,n3,n4,n5,n6,1+n7,-1+n8,n9,n10,-1+n11,0,0,0)*rat(n2*n7,n3)
         +`Z'(1,1+n2,n3,n4,n5,n6,1+n7,n8,-1+n9,n10,-1+n11,0,0,0)*rat(-n2*n7,n3)
         +`Z'(1,1+n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,0,0)*rat(-2*ep*n2-2*n2*n3-n2*n7-2*n2*n9-n2*n10+n2*n11+3*n2,n3)
         +`Z'(1,n2,1+n3,-1+n4,n5,1+n6,n7,n8,n9,n10,-1+n11,0,0,0)*rat(n6,1)
         +`Z'(1,n2,1+n3,-1+n4,n5,n6,n7,n8,n9,1+n10,-1+n11,0,0,0)*rat(2*n10,1)
         +`Z'(1,n2,1+n3,n4,n5,1+n6,n7,n8,n9,-1+n10,-1+n11,0,0,0)*rat(-n6,1)
         +`Z'(1,n2,1+n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,0,0,0)*rat(-2*ep-2*n3-2*n10+2,1)
         +`Z'(1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,0,0,0)*rat(-n6*n11+n6,n3)
         +`Z'(1,n2,n3,n4,n5,1+n6,-1+n7,n8,1+n9,n10,-1+n11,0,0,0)*rat(-n6*n9,n3)
         +`Z'(1,n2,n3,n4,n5,1+n6,n7,-1+n8,1+n9,n10,-1+n11,0,0,0)*rat(n6*n9,n3)
         +`Z'(1,n2,n3,n4,n5,1+n6,n7,n8,n9,n10,-1+n11,0,0,0)*rat(-2*ep*n6-n3*n6-2*n6*n7-n6*n9-2*n6*n10+n6*n11+3*n6,n3)
         +`Z'(1,n2,n3,n4,n5,n6,n7,n8,n9,1+n10,-1+n11,0,0,0)*rat(-2*ep*n10-2*n3*n10-2*n10^2+2*n10,n3)
         +`Y'(0,1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,0,0)*rat(-n2*n11+n2,n3)
      );
*--#] n11 : 
#endif
#if ( `il' >= 3 )
*--#[ n5 :
id,ifmatch->no13
	`Z'(1,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_
			,n8?pos_,n9?pos_,n10?pos_,1,0,0,0) =
			-RAT(2*ep*n5-2*ep+2*n5^2+2*n5*n6-6*n5-2*n6+4,n4)*(
         +`Z'(1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,1+n10,1,0,0,0)*rat(2*n10,1)
         +`Z'(1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,1,0,0,0)*rat(-2*ep-2*n4-2*n10+2,1)
         +`Z'(1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,1,0,0,0)*rat(2*ep*n6+2*n5*n6+2*n6^2-4*n6,n4)
         +`Z'(1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,1,0,0,0)*rat(-2*ep*n10-2*n4*n10-2*n10^2+2*n10,n4)
         +`Y'(1,n2,n3,n4,n5,1+n6,n7,n8,n9,n10,0,0,0,0)*rat(-2*n5*n6+2*n6,n4)
      );
*--#] n5 : 
#endif
#if ( `il' >= 4 )
*--#[ n2 :
id,ifmatch->no13
	`Z'(1,n2?{>1},n3?pos_,n4?pos_,1,n6?pos_,n7?pos_
			,n8?pos_,n9?pos_,n10?pos_,1,0,0,0) =
			-RAT(2*ep*n2-2*ep+n2^2+n2*n6-3*n2-n6+2,1)*(
         +`Z'(1,-1+n2,-1+n3,n4,1,n6,-1+n7,1+n8,1+n9,1+n10,1,0,0,0)*rat(n8*n9*n10,ep+n2-1)
         +`Z'(1,-1+n2,-1+n3,n4,1,n6,n7,1+n8,n9,1+n10,1,0,0,0)*rat(-ep*n8*n10-n8^2*n10-n8*n9*n10+n8*n10,ep+n2-1)
         +`Z'(1,-1+n2,-1+n3,n4,1,n6,n7,n8,1+n9,1+n10,1,0,0,0)*rat(-ep*n9*n10-n8*n9*n10-n9^2*n10+n9*n10,ep+n2-1)
         +`Z'(1,-1+n2,n3,-1+n4,1,n6,-1+n7,1+n8,1+n9,1+n10,1,0,0,0)*rat(-n8*n9*n10,ep+n2-1)
         +`Z'(1,-1+n2,n3,-1+n4,1,n6,n7,1+n8,n9,1+n10,1,0,0,0)*rat(ep*n8*n10+n8^2*n10+n8*n9*n10-n8*n10,ep+n2-1)
         +`Z'(1,-1+n2,n3,-1+n4,1,n6,n7,n8,1+n9,1+n10,1,0,0,0)*rat(ep*n9*n10+n8*n9*n10+n9^2*n10-n9*n10,ep+n2-1)
         +`Z'(1,-1+n2,n3,n4,1,n6,-1+n7,1+n8,1+n9,n10,1,0,0,0)*rat(5*ep*n8*n9+2*n2*n8*n9+2*n3*n8*n9+n6*n8*n9+n7*n8*n9+2*
         n8*n9^2+n8*n9*n10-8*n8*n9,ep+n2-1)
         +`Z'(1,-1+n2,n3,n4,1,n6,1+n7,-1+n8,1+n9,n10,1,0,0,0)*rat(2*ep*n7*n9+n2*n7*n9+n7*n8*n9+n7*n9^2-3*n7*n9,ep+n2-1)
         +`Z'(1,-1+n2,n3,n4,1,n6,1+n7,1+n8,-1+n9,n10,1,0,0,0)*rat(-2*ep*n7*n8-n2*n7*n8-n7*n8^2-n7*n8*n9+3*n7*n8,ep+n2-1
         )
         +`Z'(1,-1+n2,n3,n4,1,n6,1+n7,n8,n9,n10,1,0,0,0)*rat(2*ep*n7*n8-2*ep*n7*n9+n2*n7*n8-n2*n7*n9+n7*n8^2-3*n7*n8-n7
         *n9^2+3*n7*n9,ep+n2-1)
         +`Z'(1,-1+n2,n3,n4,1,n6,n7,1+n8,n9,n10,1,0,0,0)*rat(-4*ep^2*n8-ep*n2*n8-2*ep*n3*n8-ep*n6*n8-4*ep*n8^2-7*ep*n8*
         n9-ep*n8*n10+11*ep*n8+n2*n7*n8-n2*n8^2-2*n2*n8*n9+n2*n8-2*n3*n8^2-2*n3*n8*n9+2*n3*n8-n6*n8^2-n6*n8*n9+n6*n8-n7
         *n8^2-n7*n8-2*n8^2*n9-n8^2*n10+7*n8^2-2*n8*n9^2-n8*n9*n10+10*n8*n9+n8*n10-7*n8,ep+n2-1)
         +`Z'(1,-1+n2,n3,n4,1,n6,n7,n8,1+n9,n10,1,0,0,0)*rat(-6*ep^2*n9-3*ep*n2*n9-2*ep*n3*n9-ep*n6*n9-2*ep*n7*n9-5*ep*
         n8*n9-8*ep*n9^2-ep*n9*n10+15*ep*n9-n2*n7*n9-2*n2*n8*n9-3*n2*n9^2+3*n2*n9-2*n3*n8*n9-2*n3*n9^2+2*n3*n9-n6*n8*n9
         -n6*n9^2+n6*n9-2*n7*n8*n9-n7*n9^2+3*n7*n9-2*n8*n9^2-n8*n9*n10+8*n8*n9-2*n9^3-n9^2*n10+11*n9^2+n9*n10-9*n9,ep+
         n2-1)
         +`Y'(1,-1+n2,n3,n4,0,1+n6,-1+n7,1+n8,1+n9,n10,1,0,0,0)*rat(-n6*n8*n9,ep+n2-1)
         +`Y'(1,-1+n2,n3,n4,0,1+n6,n7,1+n8,n9,n10,1,0,0,0)*rat(ep*n6*n8+n6*n8^2+n6*n8*n9-n6*n8,ep+n2-1)
         +`Y'(1,-1+n2,n3,n4,0,1+n6,n7,n8,1+n9,n10,1,0,0,0)*rat(ep*n6*n9+n6*n8*n9+n6*n9^2-n6*n9,ep+n2-1)
         +`Y'(1,n2,n3,n4,0,1+n6,n7,n8,n9,n10,1,0,0,0)*rat(-n2*n6+n6,1)
         +`Y'(2,-1+n2,n3,n4,1,1+n6,n7,n8,n9,n10,0,0,0,0)*rat(-n6,1)
         +`Y'(2,n2,-1+n3,n4,1,n6,n7,n8,n9,1+n10,0,0,0,0)*rat(-n2*n10+n10,ep+n2-1)
         +`Y'(2,n2,n3,-1+n4,1,n6,n7,n8,n9,1+n10,0,0,0,0)*rat(n2*n10-n10,ep+n2-1)
         +`Y'(2,n2,n3,n4,0,1+n6,n7,n8,n9,n10,0,0,0,0)*rat(n2*n6-n6,ep+n2-1)
         +`Y'(2,n2,n3,n4,1,n6,1+n7,-1+n8,n9,n10,0,0,0,0)*rat(n2*n7-n7,ep+n2-1)
         +`Y'(2,n2,n3,n4,1,n6,1+n7,n8,-1+n9,n10,0,0,0,0)*rat(-n2*n7+n7,ep+n2-1)
         +`Y'(2,n2,n3,n4,1,n6,n7,n8,n9,n10,0,0,0,0)*rat(-5*ep*n2+5*ep-2*n2^2-2*n2*n3-n2*n6-n2*n7-2*n2*n9-n2*n10+10*n2+2
         *n3+n6+n7+2*n9+n10-8,ep+n2-1)
      );
*--#] n2 : 
#endif
#if ( `il' >= 5 )
*--#[ n6 :
id,ifmatch->no13
	`Z'(1,1,n3?pos_,n4?pos_,1,n6?{>1},n7?pos_
			,n8?pos_,n9?pos_,n10?pos_,1,0,0,0) =
			-RAT(-2*ep*n6+2*ep-n6^2+2*n6-1,1)*(
         +`Z'(1,1,-1+n3,1+n4,1,-1+n6,-1+n7,n8,1+n9,1+n10,1,0,0,0)*rat(-n4*n9*n10,ep+n6-1)
         +`Z'(1,1,-1+n3,1+n4,1,-1+n6,n7,-1+n8,1+n9,1+n10,1,0,0,0)*rat(n4*n9*n10,ep+n6-1)
         +`Z'(1,1,-1+n3,1+n4,1,-1+n6,n7,n8,n9,1+n10,1,0,0,0)*rat(-5*ep*n4*n10-n3*n4*n10-2*n4*n6*n10-2*n4*n7*n10-n4*n9*
         n10-2*n4*n10^2+7*n4*n10,ep+n6-1)
         +`Z'(1,1,1+n3,-1+n4,1,-1+n6,n7,n8,n9,1+n10,1,0,0,0)*rat(-2*ep*n3*n10-n3*n4*n10-n3*n6*n10-n3*n10^2+3*n3*n10,ep+
         n6-1)
         +`Z'(1,1,1+n3,1+n4,1,-1+n6,n7,n8,n9,-1+n10,1,0,0,0)*rat(2*ep*n3*n4+n3*n4^2+n3*n4*n6+n3*n4*n10-3*n3*n4,ep+n6-1)
         +`Z'(1,1,1+n3,n4,1,-1+n6,n7,n8,n9,n10,1,0,0,0)*rat(-2*ep*n3*n4+2*ep*n3*n10-n3*n4^2-n3*n4*n6+3*n3*n4+n3*n6*n10+
         n3*n10^2-3*n3*n10,ep+n6-1)
         +`Z'(1,1,n3,1+n4,1,-1+n6,-1+n7,n8,1+n9,n10,1,0,0,0)*rat(ep*n4*n9+n4^2*n9+n4*n9*n10-n4*n9,ep+n6-1)
         +`Z'(1,1,n3,1+n4,1,-1+n6,n7,-1+n8,1+n9,n10,1,0,0,0)*rat(-ep*n4*n9-n4^2*n9-n4*n9*n10+n4*n9,ep+n6-1)
         +`Z'(1,1,n3,1+n4,1,-1+n6,n7,n8,n9,n10,1,0,0,0)*rat(4*ep^2*n4+4*ep*n4^2+ep*n4*n6+2*ep*n4*n7+ep*n4*n9+7*ep*n4*
         n10-10*ep*n4+n3*n4^2-n3*n4*n6+n3*n4+n4^2*n6+2*n4^2*n7+n4^2*n9+2*n4^2*n10-6*n4^2+2*n4*n6*n10-n4*n6+2*n4*n7*n10-
         2*n4*n7+n4*n9*n10-n4*n9+2*n4*n10^2-9*n4*n10+6*n4,ep+n6-1)
         +`Z'(1,1,n3,n4,1,-1+n6,-1+n7,n8,1+n9,1+n10,1,0,0,0)*rat(ep*n9*n10+n4*n9*n10+n9*n10^2-n9*n10,ep+n6-1)
         +`Z'(1,1,n3,n4,1,-1+n6,n7,-1+n8,1+n9,1+n10,1,0,0,0)*rat(-ep*n9*n10-n4*n9*n10-n9*n10^2+n9*n10,ep+n6-1)
         +`Z'(1,1,n3,n4,1,-1+n6,n7,n8,n9,1+n10,1,0,0,0)*rat(6*ep^2*n10+2*ep*n3*n10+5*ep*n4*n10+3*ep*n6*n10+2*ep*n7*n10+
         ep*n9*n10+8*ep*n10^2-14*ep*n10+2*n3*n4*n10+n3*n6*n10+n3*n10^2-3*n3*n10+2*n4*n6*n10+2*n4*n7*n10+n4*n9*n10+2*n4*
         n10^2-7*n4*n10+3*n6*n10^2-3*n6*n10+2*n7*n10^2-2*n7*n10+n9*n10^2-n9*n10+2*n10^3-10*n10^2+8*n10,ep+n6-1)
         +`Y'(0,2,n3,n4,2,-1+n6,n7,n8,n9,n10,1,0,0,0)*rat(-1,1)
         +`Y'(1,1,1+n3,-1+n4,2,n6,n7,n8,n9,n10,0,0,0,0)*rat(-n3*n6+n3,ep+n6-1)
         +`Y'(1,1,1+n3,n4,2,n6,n7,n8,n9,-1+n10,0,0,0,0)*rat(n3*n6-n3,ep+n6-1)
         +`Y'(1,1,n3,n4,2,n6,-1+n7,n8,1+n9,n10,0,0,0,0)*rat(n6*n9-n9,ep+n6-1)
         +`Y'(1,1,n3,n4,2,n6,n7,-1+n8,1+n9,n10,0,0,0,0)*rat(-n6*n9+n9,ep+n6-1)
         +`Y'(1,1,n3,n4,2,n6,n7,n8,n9,n10,0,0,0,0)*rat(5*ep*n6-5*ep+n3*n6-n3+2*n6^2+2*n6*n7+n6*n9+2*n6*n10-9*n6-2*n7-n9
         -2*n10+7,ep+n6-1)
         +`Y'(1,2,n3,n4,2,-1+n6,n7,n8,n9,n10,0,0,0,0)*rat(1,1)
      );
*--#] n6 : 
#endif
#if ( `il' >= 6 )
*--#[ n3 :
id,ifmatch->no13
	`Z'(1,1,n3?{>1},n4?pos_,1,1,n7?pos_
			,n8?pos_,n9?pos_,n10?pos_,1,0,0,0) =
			-RAT(2*ep*n3-2*ep+2*n3^2+n3*n9+n3*n10-6*n3-n9-n10+4,n7)*(
         +`Z'(1,1,-1+n3,n4,1,1,1+n7,-1+n8,1+n9,n10,1,0,0,0)*rat(n9,1)
         +`Z'(1,1,-1+n3,n4,1,1,1+n7,n8,n9,n10,1,0,0,0)*rat(-2*ep-2*n7-n9-n10+2,1)
         +`Z'(1,1,-1+n3,n4,1,1,n7,n8,1+n9,n10,1,0,0,0)*rat(n3*n9-n7*n9-n9,n7)
         +`Z'(1,1,-1+n3,n4,1,1,n7,n8,n9,1+n10,1,0,0,0)*rat(n3*n10-n7*n10-n10,n7)
         +`Z'(1,1,n3,-1+n4,1,1,n7,n8,n9,1+n10,1,0,0,0)*rat(-n3*n10+n10,n7)
         +`Y'(1,0,n3,n4,1,1,n7,n8,1+n9,n10,1,0,0,0)*rat(-n3*n9+n9,n7)
         +`Y'(1,1,-1+n3,n4,1,0,1+n7,n8,n9,1+n10,1,0,0,0)*rat(n10,1)
      );
*--#] n3 : 
#endif
#if ( `il' >= 7 )
*--#[ n4 :
id,ifmatch->no13
	`Z'(1,1,1,n4?{>1},1,1,n7?pos_
			,n8?pos_,n9?pos_,n10?pos_,1,0,0,0) =
			-RAT(2*ep*n4-2*ep+n4*n7+n4*n9+n4*n10-2*n4-n7-n9-n10+2,n7)*(
         +`Z'(1,1,1,-1+n4,1,1,1+n7,n8,n9,n10,1,0,0,0)*rat(-2*ep-n4-n10+2,1)
         +`Z'(1,1,1,-1+n4,1,1,n7,n8,n9,1+n10,1,0,0,0)*rat(-n4*n10+n7*n10+n10,n7)
         +`Y'(1,0,1,n4,1,1,n7,n8,1+n9,n10,1,0,0,0)*rat(-n4*n9+n9,n7)
         +`Y'(1,1,0,n4,1,0,1+n7,n8,n9,1+n10,1,0,0,0)*rat(n4*n10-n10,ep)
         +`Y'(1,1,0,n4,1,1,1+n7,n8,n9,n10,1,0,0,0)*rat(n4-1,1)
         +`Y'(1,1,0,n4,1,1,n7,n8,1+n9,n10,1,0,0,0)*rat(n4*n9-n9,n7)
         +`Y'(1,1,0,n4,1,1,n7,n8,n9,1+n10,1,0,0,0)*rat(n4*n10-n10,n7)
         +`Y'(1,1,1,-1+n4,1,0,1+n7,n8,n9,1+n10,1,0,0,0)*rat(-2*ep*n10-n4*n10-n10^2+2*n10,ep)
         +`Y'(1,1,1,-1+n4,2,1,1+n7,n8,n9,n10,0,0,0,0)*rat(1,1)
         +`Y'(1,1,1,n4,1,0,1+n7,n8,n9,n10,1,0,0,0)*rat(-2*ep*n4+2*ep-n4^2-n4*n10+3*n4+n10-2,ep)
      );
*--#] n4 : 
#endif
#if ( `il' >= 8 )
*--#[ n7 :
id,ifmatch->no13
	`Z'(1,1,1,1,1,1,n7?{>1},n8?pos_,n9?pos_,n10?pos_,1,0,0,0) =
			-RAT(2*ep*n7-2*ep+n7*n8+n7*n9-n7-n8-n9+1,n8)*(
         +`Z'(1,1,1,1,1,1,-1+n7,1+n8,n9,n10,1,0,0,0)*rat(-2*ep-n7-n9-n10+3,1)
         +`Y'(1,0,1,1,1,1,-1+n7,1+n8,1+n9,n10,1,0,0,0)*rat(ep*n9-n7*n9+n9,ep)
         +`Y'(1,0,1,1,1,1,n7,1+n8,n9,n10,1,0,0,0)*rat(2*ep*n7-2*ep+n7*n8+n7*n9-n7-n8-n9+1,ep)
         +`Y'(1,0,1,1,1,1,n7,n8,1+n9,n10,1,0,0,0)*rat(2*ep*n7*n9-2*ep*n9+n7*n8*n9+n7*n9^2-n7*n9-n8*n9-n9^2+n9,ep*n8)
         +`Y'(1,1,0,1,1,1,-1+n7,1+n8,1+n9,n10,1,0,0,0)*rat(-n9,1)
         +`Y'(1,1,0,1,1,1,-1+n7,1+n8,n9,1+n10,1,0,0,0)*rat(-n10,1)
         +`Y'(1,1,0,1,1,1,n7,1+n8,n9,n10,1,0,0,0)*rat(-n7+1,1)
         +`Y'(1,1,0,1,1,1,n7,n8,1+n9,n10,1,0,0,0)*rat(-n7*n9+n9,n8)
         +`Y'(1,1,1,0,1,1,-1+n7,1+n8,n9,1+n10,1,0,0,0)*rat(n10,1)
         +`Y'(2,1,1,1,1,1,n7,n8,n9,n10,0,0,0,0)*rat(-n7+1,n8)
      );
*--#] n7 : 
#endif
#if ( `il' >= 9 )
*--#[ n10 :
id,ifmatch->no13
	`Z'(1,1,1,1,1,1,1,n8?pos_,n9?pos_,n10?{>1},1,0,0,0) =
			-RAT(2*ep+2*n10-2,1)*(
         +`Z'(1,1,1,1,1,1,1,1+n8,n9,-1+n10,1,0,0,0)*rat(-4*ep*n8*n9+4*ep*n8*n10-4*ep*n8-2*n8*n9^2+2*n8*n9+2*n8*n10^2-6*
         n8*n10+4*n8,2*ep*n10-2*ep+n8*n10-n8+n9*n10-n9-n10+1)
         +`Z'(1,1,1,1,1,1,1,n8,1+n9,-1+n10,1,0,0,0)*rat(-4*ep^2*n9-6*ep*n9^2-2*n9^3,2*ep*n10-2*ep+n8*n10-n8+n9*n10-n9-
         n10+1)
         +`Y'(1,0,1,1,1,1,1,1+n8,1+n9,-1+n10,1,0,0,0)*rat(2*ep*n8*n9^2-2*ep*n8*n9*n10+2*ep*n8*n9-2*n8*n9^2+2*n8*n9*n10-
         2*n8*n9,2*ep^2*n10-2*ep^2+ep*n8*n10-ep*n8+ep*n9*n10-ep*n9-ep*n10+ep)
         +`Y'(1,0,1,1,1,1,1,n8,2+n9,-1+n10,1,0,0,0)*rat(-2*ep^2*n8*n9^2-2*ep^2*n8*n9+2*ep^2*n9^2+2*ep^2*n9-2*ep*n8*n9^3
         +2*ep*n8*n9+2*ep*n9^3-2*ep*n9+2*n8*n9^3+2*n8*n9^2-2*n9^3-2*n9^2,4*ep^3*n10-4*ep^3+2*ep^2*n8*n10-2*ep^2*n8+4*
         ep^2*n9*n10-4*ep^2*n9+2*ep^2*n10^2-6*ep^2*n10+4*ep^2+ep*n8*n9*n10-ep*n8*n9+ep*n8*n10^2-2*ep*n8*n10+ep*n8+ep*
         n9^2*n10-ep*n9^2+ep*n9*n10^2-3*ep*n9*n10+2*ep*n9-ep*n10^2+2*ep*n10-ep)
         +`Y'(1,0,1,1,1,1,2,-1+n8,2+n9,-1+n10,1,0,0,0)*rat(-2*ep*n9^2-2*ep*n9-2*n9^3-2*n9^2,2*ep^2*n10-2*ep^2+ep*n9*n10
         -ep*n9+ep*n10^2-2*ep*n10+ep)
         +`Y'(1,0,1,1,1,1,2,1+n8,n9,-1+n10,1,0,0,0)*rat(2*n8*n9-2*n8*n10+2*n8,ep*n10-ep)
         +`Y'(1,0,1,1,1,1,2,n8,1+n9,-1+n10,1,0,0,0)*rat(-2*ep*n8*n9+4*ep*n9^2-4*ep*n9*n10+6*ep*n9-2*n8*n9^2+2*n9^3+2*
         n9^2-2*n9*n10^2+4*n9*n10-2*n9,2*ep^2*n10-2*ep^2+ep*n9*n10-ep*n9+ep*n10^2-2*ep*n10+ep)
         +`Y'(1,0,2,1,1,1,1,n8,1+n9,-1+n10,1,0,0,0)*rat(2*ep*n9+2*n9*n10-2*n9,2*ep*n10-2*ep+n9*n10-n9+n10^2-2*n10+1)
         +`Y'(1,1,0,1,1,1,1,1+n8,1+n9,-1+n10,1,0,0,0)*rat(-2*n8*n9^2+2*n8*n9*n10-2*n8*n9,2*ep*n10-2*ep+n8*n10-n8+n9*n10
         -n9-n10+1)
         +`Y'(1,1,0,1,1,1,1,1+n8,n9,n10,1,0,0,0)*rat(-2*n8*n9+2*n8*n10-2*n8,2*ep+n8+n9-1)
         +`Y'(1,1,0,1,1,1,1,n8,1+n9,n10,1,0,0,0)*rat(2*ep*n8*n9-2*ep*n9+2*n8*n9^2-2*n9^2,4*ep^2+2*ep*n8+4*ep*n9+2*ep*
         n10-4*ep+n8*n9+n8*n10-n8+n9^2+n9*n10-2*n9-n10+1)
         +`Y'(1,1,0,1,1,1,1,n8,2+n9,-1+n10,1,0,0,0)*rat(2*ep*n8*n9^2+2*ep*n8*n9-2*ep*n9^2-2*ep*n9+2*n8*n9^3+2*n8*n9^2-2
         *n9^3-2*n9^2,4*ep^2*n10-4*ep^2+2*ep*n8*n10-2*ep*n8+4*ep*n9*n10-4*ep*n9+2*ep*n10^2-6*ep*n10+4*ep+n8*n9*n10-n8*
         n9+n8*n10^2-2*n8*n10+n8+n9^2*n10-n9^2+n9*n10^2-3*n9*n10+2*n9-n10^2+2*n10-1)
         +`Y'(1,1,0,1,1,1,2,-1+n8,2+n9,-1+n10,1,0,0,0)*rat(2*ep*n9^2+2*ep*n9+2*n9^3+2*n9^2,4*ep^2*n10-4*ep^2+2*ep*n8*
         n10-2*ep*n8+4*ep*n9*n10-4*ep*n9+2*ep*n10^2-6*ep*n10+4*ep+n8*n9*n10-n8*n9+n8*n10^2-2*n8*n10+n8+n9^2*n10-n9^2+n9
         *n10^2-3*n9*n10+2*n9-n10^2+2*n10-1)
         +`Y'(1,1,0,1,1,1,2,1+n8,n9,-1+n10,1,0,0,0)*rat(-2*n8*n9+2*n8*n10-2*n8,2*ep*n10-2*ep+n8*n10-n8+n9*n10-n9-n10+1)
         +`Y'(1,1,0,1,1,1,2,n8,1+n9,-1+n10,1,0,0,0)*rat(2*ep*n8*n9-4*ep*n9^2+4*ep*n9*n10-6*ep*n9+2*n8*n9^2-2*n9^3-2*
         n9^2+2*n9*n10^2-4*n9*n10+2*n9,4*ep^2*n10-4*ep^2+2*ep*n8*n10-2*ep*n8+4*ep*n9*n10-4*ep*n9+2*ep*n10^2-6*ep*n10+4*
         ep+n8*n9*n10-n8*n9+n8*n10^2-2*n8*n10+n8+n9^2*n10-n9^2+n9*n10^2-3*n9*n10+2*n9-n10^2+2*n10-1)
         +`Y'(1,1,1,0,1,1,1,1+n8,n9,n10,1,0,0,0)*rat(2*n8*n9-2*n8*n10+2*n8,2*ep+n8+n9-1)
         +`Y'(1,1,1,0,1,1,1,n8,1+n9,n10,1,0,0,0)*rat(-2*ep*n8*n9+2*ep*n9-2*n8*n9^2+2*n9^2,4*ep^2+2*ep*n8+4*ep*n9+2*ep*
         n10-4*ep+n8*n9+n8*n10-n8+n9^2+n9*n10-2*n9-n10+1)
         +`Y'(1,1,1,1,1,0,2,n8,n9,n10,1,0,0,0)*rat(-2*ep-2*n10+2,2*ep+n9+n10-1)
         +`Y'(1,1,2,0,1,1,1,n8,n9,n10,1,0,0,0)*rat(-2*ep-2*n9,2*ep+n9+n10-1)
         +`Y'(2,1,1,1,1,1,2,-1+n8,1+n9,-1+n10,0,0,0,0)*rat(2*ep*n9+2*n9^2,4*ep^2*n10-4*ep^2+2*ep*n8*n10-2*ep*n8+4*ep*n9
         *n10-4*ep*n9+2*ep*n10^2-6*ep*n10+4*ep+n8*n9*n10-n8*n9+n8*n10^2-2*n8*n10+n8+n9^2*n10-n9^2+n9*n10^2-3*n9*n10+2*
         n9-n10^2+2*n10-1)
         +`Y'(2,1,1,1,1,1,2,n8,n9,-1+n10,0,0,0,0)*rat(-2*n9+2*n10-2,2*ep*n10-2*ep+n8*n10-n8+n9*n10-n9-n10+1)
      );
*--#] n10 : 
#endif
#if ( `il' >= 10 )
*--#[ n9 :
id,ifmatch->no13,`Z'(1,1,1,1,1,1,1,1,1,1,1,0,0,0) = Master(no1);
id,ifmatch->no13
	`Z'(1,1,1,1,1,1,1,n8?pos_,n9?pos_,1,1,0,0,0) =
			-RAT(-36*ep^2*n8-16*ep^2*n9+52*ep^2-18*ep*n8^2-30*ep*n8*n9+48*ep*n8-20*ep*
         n9^2+66*ep*n9-46*ep-2*n8^3-7*n8^2*n9+9*n8^2-8*n8*n9^2+23*n8*n9-15*n8-4*n9^3+16*n9^2-20*n9+8,20*ep^2+2*ep*n8+14
         *ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)*(
*      +U(0)*(
         +`Z'(1,1,1,1,1,1,1,1+n8,-1+n9,1,1,0,0,0)*rat(-n8,1)
         +`Y'(1,1,0,2,1,1,1,n8,n9,1,1,0,0,0)*rat(-60*ep^2-10*ep*n8-26*ep*n9+14*ep-n8*n9-n8-2*n9^2+2,20*ep^2+2*ep*n8+14*
         ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,1,1,1,0,1+n8,n9,1,1,0,0,0)*rat(60*ep^2*n8+22*ep*n8^2+38*ep*n8*n9-38*ep*n8+2*n8^3+7*n8^2*n9-7*n8^2+6
         *n8*n9^2-12*n8*n9+6*n8,20*ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,1,2,0,1,n8,n9,1,1,0,0,0)*rat(-40*ep^2-20*ep*n8-24*ep*n9+26*ep-2*n8^2-6*n8*n9+6*n8-4*n9^2+9*n9-5,20*
         ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,1,2,1,1,n8,n9,1,0,0,0,0)*rat(8*ep+2*n8+2*n9-2,2*ep+n9-1)
         +`Y'(1,1,1,2,1,1,1,n8,n9,0,1,0,0,0)*rat(1,1)
         +`Y'(2,0,1,1,1,1,1,n8,n9,1,1,0,0,0)*rat(40*ep^2+8*ep*n8+12*ep*n9-2*ep+2*n8+3*n9-3,20*ep^2+2*ep*n8+14*ep*n9-12*
         ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(2,1,1,1,1,1,1,n8,n9,1,0,0,0,0)*rat(-8*ep-2*n8-2*n9+2,2*ep+n9-1)
*         )
*      +U(1)*(
         +`Z'(1,1,1,1,1,1,1,-1+n8,n9,1,1,0,0,0)*rat(120*ep^3+104*ep^2*n8+136*ep^2*n9-196*ep^2+26*ep*n8^2+74*ep*n8*n9-96
         *ep*n8+50*ep*n9^2-138*ep*n9+88*ep+2*n8^3+9*n8^2*n9-11*n8^2+13*n8*n9^2-33*n8*n9+20*n8+6*n9^3-24*n9^2+30*n9-12,
         20*ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Z'(1,1,1,1,1,1,1,n8,-1+n9,1,1,0,0,0)*rat(-2*ep-n8-n9+2,1)
         +`Y'(0,1,1,1,1,1,0,1+n8,n9,1,1,0,0,0)*rat(-40*ep^2*n8-20*ep*n8^2-24*ep*n8*n9+26*ep*n8-2*n8^3-6*n8^2*n9+6*n8^2-
         4*n8*n9^2+9*n8*n9-5*n8,20*ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(0,1,1,1,1,1,1,n8,n9,1,1,0,0,0)*rat(-16*ep^2-12*ep*n8-12*ep*n9+20*ep-2*n8^2-4*n8*n9+6*n8-2*n9^2+6*n9-4,2*
         ep+n9-1)
         +`Y'(0,2,0,1,1,1,1,n8,n9,1,1,0,0,0)*rat(-40*ep^2-8*ep*n8-12*ep*n9+2*ep-2*n8-3*n9+3,20*ep^2+2*ep*n8+14*ep*n9-12
         *ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,0,1,1,1,1,0,n8,1+n9,1,1,0,0,0)*rat(-40*ep^2*n9-8*ep*n8*n9-12*ep*n9^2+2*ep*n9-2*n8*n9-3*n9^2+3*n9,20*
         ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,0,1,1,1,1,1,n8,n9,1,1,0,0,0)*rat(80*ep^3+56*ep^2*n8+64*ep^2*n9-84*ep^2+8*ep*n8^2+20*ep*n8*n9-14*ep*n8+
         12*ep*n9^2-20*ep*n9-2*ep+2*n8^2+5*n8*n9-7*n8+3*n9^2-9*n9+6,20*ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9
         +1)
         +`Y'(1,1,0,1,1,0,1,n8,n9,2,1,0,0,0)*rat(40*ep^2+20*ep*n8+24*ep*n9-26*ep+2*n8^2+6*n8*n9-6*n8+4*n9^2-9*n9+5,20*
         ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,0,1,1,1,0,n8,1+n9,1,1,0,0,0)*rat(20*ep^2*n9+6*ep*n8*n9-2*ep*n9^2+10*ep*n9-n8*n9^2+3*n8*n9-2*n9^3+6*
         n9^2-4*n9,20*ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,0,1,1,1,0,n8,n9,2,1,0,0,0)*rat(-20*ep^2-18*ep*n8-10*ep*n9+14*ep-2*n8^2-5*n8*n9+5*n8-2*n9^2+6*n9-4,20*
         ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,0,1,1,1,1,-1+n8,1+n9,1,1,0,0,0)*rat(n9,1)
         +`Y'(1,1,0,1,1,1,1,n8,n9,1,1,0,0,0)*rat(120*ep^3+36*ep^2*n8+8*ep^2*n9+60*ep^2+18*ep*n8-14*ep*n9^2+46*ep*n9-24*
         ep-n8*n9^2+3*n8*n9-2*n9^3+6*n9^2-4*n9,20*ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,0,2,0,1,1,n8,n9,1,1,0,0,0)*rat(40*ep^2+8*ep*n8+12*ep*n9-2*ep+2*n8+3*n9-3,20*ep^2+2*ep*n8+14*ep*n9-12*
         ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,0,1,1,0,n8,n9,2,1,0,0,0)*rat(-1,1)
         +`Y'(1,1,1,0,1,1,1,n8,n9,1,1,0,0,0)*rat(-120*ep^3-20*ep^2*n8-52*ep^2*n9+28*ep^2-2*ep*n8*n9-2*ep*n8-4*ep*n9^2+4
         *ep,20*ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,0,2,1,1,n8,n9,1,0,0,0,0)*rat(-40*ep^2-8*ep*n8-12*ep*n9+2*ep-2*n8-3*n9+3,20*ep^2+2*ep*n8+14*ep*n9-12
         *ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,1,0,1,1,n8,n9,1,1,0,0,0)*rat(16*ep^2+4*ep*n8+4*ep*n9-4*ep,2*ep+n9-1)
         +`Y'(1,1,1,1,0,2,0,n8,n9,1,1,0,0,0)*rat(40*ep^2+20*ep*n8+24*ep*n9-26*ep+2*n8^2+6*n8*n9-6*n8+4*n9^2-9*n9+5,20*
         ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,1,1,0,1,n8,n9,1,1,0,0,0)*rat(-80*ep^3-40*ep^2*n8-48*ep^2*n9+52*ep^2-4*ep*n8^2-12*ep*n8*n9+12*ep*n8-
         8*ep*n9^2+18*ep*n9-10*ep,20*ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,1,1,1,0,n8,n9,1,1,0,0,0)*rat(-120*ep^3-128*ep^2*n8-80*ep^2*n9+104*ep^2-30*ep*n8^2-74*ep*n8*n9+78*ep
         *n8-18*ep*n9^2+58*ep*n9-40*ep-2*n8^3-9*n8^2*n9+9*n8^2-11*n8*n9^2+24*n8*n9-13*n8-2*n9^3+11*n9^2-16*n9+7,20*ep^2
         +2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,1,1,1,1,n8,n9,0,1,0,0,0)*rat(2*ep,1)
         +`Y'(1,1,1,1,1,1,1,n8,n9,1,0,0,0,0)*rat(136*ep^2*n8+96*ep^2*n9-232*ep^2+48*ep*n8^2+108*ep*n8*n9-146*ep*n8+56*
         ep*n9^2-180*ep*n9+114*ep+4*n8^3+16*n8^2*n9-18*n8^2+20*n8*n9^2-49*n8*n9+27*n8+8*n9^3-32*n9^2+37*n9-13,20*ep^2+2
         *ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,1,1,2,0,n8,n9,1,0,0,0,0)*rat(-20*ep^2-18*ep*n8-10*ep*n9+14*ep-2*n8^2-5*n8*n9+5*n8-2*n9^2+6*n9-4,20*
         ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,1,1,1,1,2,1,n8,n9,0,0,0,0,0)*rat(-1,1)
         +`Y'(1,2,0,1,1,1,1,n8,n9,1,0,0,0,0)*rat(20*ep^2+6*ep*n8-2*ep*n9+10*ep-n8*n9+3*n8-2*n9^2+6*n9-4,20*ep^2+2*ep*n8
         +14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
         +`Y'(1,2,1,1,1,1,1,n8,-1+n9,1,0,0,0,0)*rat(1,1)
         +`Y'(2,1,1,1,1,1,1,-1+n8,n9,1,0,0,0,0)*rat(40*ep^2+20*ep*n8+24*ep*n9-26*ep+2*n8^2+6*n8*n9-6*n8+4*n9^2-9*n9+5,
         20*ep^2+2*ep*n8+14*ep*n9-12*ep+n8*n9-n8+2*n9^2-3*n9+1)
*         );
      );
*--#] n9 : 
#endif
GoTo no14;
Label no13;
$ired = 1;
Label no14;
id	`Z'(?a,x?neg0_,?b,n12?,n13?,n14?) = `Y'(?a,x,?b,n12,n13,n14);
*id	`Y'(?a) = 0;
ModuleOption,maximum,$ired;
*Print +f +s;
.sort:Main no1 loop, pass `$pass';
#if ( `$ired' > 0 )
	#redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
id	`Z'(1,1,1,1,1,1,1,1,1,1,1,0,0,0) = Master(no1);
*
*--#] Loop 2 : 
#break
*--#] no1 : 
*--#[ no6 :
#case no6
*
*--#[ Loop 1 :
#do il = 1,5
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering no6 loop 1-`il', pass `$pass' at `time_' sec."
#endif
.sort
#if ( `il' >= 1 )
#if(1==1)
*--#[ stat :
*
*  Prework to make the reduction of the no6 topology faster.
*
* B13
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?{>1},n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n10)*(`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,n14)*rat(-n10+1,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,
   n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n10+1,1)+`Z'(n1,n2,n3,n4,n5,n6,1+n7,-1+n8,n9,-1+n10,n11,1+n12,
   n13,n14)*rat(-n7,1)+`Z'(n1,n2,n3,n4,n5,n6,1+n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-
   1+n9,n10,n11,1+n12,n13,n14)*rat(-1+n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,-1+n13,1+n14)*rat(-n14,1)+`Z'(
   n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(2*ep+n3+n7+2*n9+n10+n14-5,1));
* B15
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?pos_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n9)*(`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n9+1,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+
   n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,n2,
   n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,
   n14)*rat(-n9+1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,1+n7,
   n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n7,1)+`Z'(n1,n2,n3,n4,n5,n6,1+n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(
   n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,
   n14)*rat(2*ep+n3+n7+n9+2*n10+n14-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(-1+n9,1));
* B8
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,1+n2,-1+n3,n4,n5,-1+n6,n7,n8,
   n9,n10,n11,1+n12,n13,n14)*rat(n2,1)+`Z'(n1,1+n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n2,1)+`Z'(n1,n2,-
   1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(1+n12,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,2+n12,n13,
   n14)*rat(1+n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,
   n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(n1,
   n2,n3,n4,n5,-1+n6,1+n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n7,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,-1+n8,n9,n10,1+n11,1+n12,
   n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,
   n7,n8,n9,-1+n10,n11,2+n12,n13,n14)*rat(-n12-1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,-1+n13,n14)*rat(-
   n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,n13,-1+n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,
   n11,1+n12,1+n13,-1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*ep-n2-n7-2*n9-n12-
   n13+3,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,
   1+n12,n13,n14)*rat(-1+n6,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(-1+n6,1));
* B19
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n5)*(`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,1+n2,-1+n3,n4,-1+n5,n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(n1,1+n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(n2,1)+`Z'(n1,
   n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,
   n13,1+n14)*rat(-2*n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5+1,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,
   n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5+1,1)+
   `Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,n10,1+n11,
   n12,n13,1+n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,-1+n8,n9,n10,1+n11,n12,n13,1+n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,-1+
   n5,n6,n7,n8,-1+n9,n10,1+n11,n12,n13,1+n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*
   rat(n12,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,
   n10,1+n11,n12,-1+n13,1+n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,
   n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-
   2*ep-n2-2*n3-n4-n11-n14+3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-1+n5,1));
* B10
id,ifmatch->no61,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(n1,-2+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,
   n10,n11,2+n12,n13,n14)*rat(1+n12,1)+`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+
   n2,1+n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,
   n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,
   n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+
   `Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+
   n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)+`Z'(n1,-1+n2,n3,n4,
   n5,1+n6,-1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n6,1)+`Z'(n1,-1+n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*
   rat(n6,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,
   n10,1+n11,1+n12,n13,n14)*rat(-n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,2+n12,n13,n14)*rat(-n12-1,1)+`Z'(n1,-1+
   n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,-1+n14)*rat(n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,
   n14)*rat(-n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,
   n10,n11,1+n12,1+n13,-1+n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+
   n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*
   rat(2*ep+n2+n3+n6+2*n10+n14-5,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n2,1)+`Z'(n1,n2,n3,n4,
   n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n2+1,1));
* B12
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n3+1)*(`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,1+n11,n12,n13,n14)*rat(-1+n3,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,
   n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,1+n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,
   n13,n14)*rat(-n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,1+n9,n10,1+n11,n12,n13,n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,
   -1+n8,n9,n10,1+n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,-1+n13,1+n14)*rat(n14,1)+
   `Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,
   n13,n14)*rat(-2*ep-n3-2*n7-n9-n10-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n3+1,1)+`Z'(n1,
   n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,-1+n13,n14)
   *rat(-1+n3,1));
* B4
id,ifmatch->no61,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n11,1)+`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,
   n8,n9,n10,n11,n12,1+n13,n14)*rat(-n2,1)+`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(n2,1)+`Z'(-1+
   n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,
   1+n13,n14)*rat(-n12,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(-1+n1,n2,n3,n4,n5,
   n6,-1+n7,1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-
   n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,-1+n9,
   n10,n11,n12,1+n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,n12,1+n13,n14)*rat(2*n11,1)+`Z'(-1+n1,n2,
   n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,-1+
   n14)*rat(-n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,n11,n12,1+n13,n14)*rat(2*ep+n2+n8+2*n9+n12+n13-3,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,2+n13,-1+n14)*
   rat(-n13-1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,
   n10,n11,n12,1+n13,n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,
   n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,-1+n14)*rat(-n1+1,1));
* B7
id,ifmatch->no61,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,
   n8,n9,n10,2+n11,n12,n13,n14)*rat(-2*n11-2,1)+`Z'(n1,-2+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(
   n1,-1+n2,1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,-1+n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,1+n11,
   n12,n13,n14)*rat(n4,1)+`Z'(n1,-1+n2,n3,1+n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,
   n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+
   `Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,1+n9,n10,1+n11,n12,n13,n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,1+n11,
   n12,n13,n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,-1+n2,n3,n4,n5,
   n6,n7,-1+n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,
   1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,
   n12,n13,1+n14)*rat(-2*n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(6*ep+2*n2+2*n3+n4+2*n5+2*n6+
   2*n7+n9+2*n10+n11+n12+n13+n14-13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1));
* B6
id,ifmatch->no61,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,
   n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,-1+n14)*rat(-n3,1)+`Z'(
   n1,-1+n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n4,1)+`Z'(n1,-1+n2,n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+
   n13,n14)*rat(-n4,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(n1,-1+n2,n3,n4,n5,1+n6,-
   1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6,1)+`Z'(n1,-1+n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,n11,n12,1+n13,n14)*rat(n6,1)+`Z'(
   n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,
   n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,n4,
   n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,-1+n14)*rat(n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-
   n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+
   n12,1+n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,
   n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-2*n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(4*ep+n2+2*
   n3+n4+2*n5+n6+2*n10+n11+n12+n14-9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,2+n13,-1+n14)*rat(1+n13,1)+`Z'(n1,-1+
   n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,2+n13,n14)*rat(-2*n13-2,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)
   *rat(-1+n2,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-n2+1,1));
* B20
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n6+1)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n6,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,
   n10,1+n11,n12,n13,n14)*rat(-1+n6,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,
   n3,1+n4,-1+n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,n2,n3,1+n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,
   n14)*rat(n4,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,
   n9,n10,1+n11,n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,1+n10,1+n11,n12,n13,-1+n14)*rat(-n10,1)+`Z'(n1,
   n2,n3,n4,n5,-1+n6,n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,-1+
   n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-2*ep-n4-2*n5-n6-n10-n12+5,1));
* B16
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6+1,1)+`Z'(n1,1+n2,-1+n3,n4,n5,-1+n6,n7,n8,
   n9,n10,n11,n12,1+n13,n14)*rat(n2,1)+`Z'(n1,1+n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-n2,1)+`Z'(n1,n2,-
   1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,
   n14)*rat(-n10,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,n11,n12,1+n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,1+n4,-1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n4,1)+`Z'(n1,n2,
   n3,1+n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n4,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,
   n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n6,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,
   n9,n10,1+n11,n12,1+n13,n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,-1+n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,
   n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,-1+n10,n11,1+n12,
   1+n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,-1+n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,
   n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,1+n13,-1+n14)*rat(n11,1)+
   `Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,
   n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(4*ep+n2+2*n3+n4+2*n5+n6+2*n10+n11+n12+
   n14-9,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,2+n13,-1+n14)*rat(1+n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,
   n11,n12,1+n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n6,1)+`Z'(n1,n2,n3,n4,n5,
   n6,n7,n8,n9,-1+n10,n11,n12,1+n13,n14)*rat(-1+n6,1));
* B3
id,ifmatch->no61,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n1+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,1+n11,n12,n13,n14)*rat(-n9,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,2+n11,n12,n13,n14)*rat(-2*n11-2,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,1+n11,n12,n13,n14)*rat(n9,1)+`Z'(-1+
   n1,n2,n3,n4,n5,n6,n7,1+n8,n9,n10,1+n11,n12,-1+n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9,n10,n11,n12,n13,
   n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-2*ep-2*n2-n9-n11-n12-n13+3,1));
* B18
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n4+1)*(`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,n11,2+n12,n13,n14)*rat(-2*n12-2,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,
   n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,n14)*rat(n5,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,n11,n12,n13,
   n14)*rat(-n5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,n14)*rat(-2*ep-2*n2-n3-n11-n12-n13+3,1));
* B17
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n4+1)*(`Z'(-1+n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n5,1)+`Z'(n1,-1+n2,n3,-1+n4,1+n5,n6,n7,
   n8,n9,n10,1+n11,n12,n13,n14)*rat(n5,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(2*n12,1)+`Z'(n1,
   -1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n4,1)+`Z'(n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,n9,n10,1+n11,n12,
   n13,n14)*rat(-n5,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n5,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,
   n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(n10,1)+`Z'(
   n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,-1+n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,n12,
   n13,n14)*rat(-2*ep-n4-n5-2*n6-n10-n12+5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,
   n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,n13,n14)
   *rat(-1+n4,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,-1+n13,n14)*rat(-1+n4,1));
* B5
id,ifmatch->no61,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?pos_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n1+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n2,1)+`Z'(-1+n1,
   1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n2,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,
   n13,n14)*rat(n11,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(1+n12,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,
   n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(1+n12,1)+`Z'(
   -1+n1,n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,1+
   n12,n13,n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,
   -1+n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(n7,1)
   +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,-1+n9,n10,
   n11,1+n12,n13,n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,
   n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(-2*n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,2+n12,n13,n14)*rat(-
   n12-1,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,-1+n14)*rat(n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,
   1+n11,n12,n13,n14)*rat(-n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,-1+n14)*rat(n13,1)+`Z'(-1+n1,n2,n3,
   n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-4*ep-n1-n2-2*n6-n7-2*n8-2*n9-n11-n12-n13+8,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,
   n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n1,1)+`Z'(
   n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,
   n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+
   n9,n10,n11,1+n12,n13,n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,-1+n14)*rat(-1+n1,1));
* B11
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?pos_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n9)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(n8,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+
   n9,n10,n11,n12,1+n13,-1+n14)*rat(-n3,1)+`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(n4,1)+`Z'(n1,
   n2,n3,1+n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-n4,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,1+n8,-1+n9,n10,n11,n12,1+n13,
   n14)*rat(-n8,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,
   n10,n11,n12,1+n13,1+n14)*rat(-2*n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(8*ep+2*n1+2*n2+2*
   n3+n4+2*n5+2*n6+2*n7+n8+2*n9+2*n10+2*n11+2*n12+n13+n14-17,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,-1+n14)*
   rat(-n9+1,1));
* B2
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?pos_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n7)*(`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,n11,n12,1+n13,n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n10,1)+`Z'(n1,n2,
   n3,1+n4,-1+n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n4,1)+`Z'(n1,n2,n3,1+n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,
   n14)*rat(n4,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+
   n8,n9,n10,n11,n12,1+n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,-1+n14)*rat(-n10,1)+`Z'(n1,
   n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,n14)
   *rat(n12,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-8*ep-2*n1-2*n2-2*n3-n4-2*n5-2*n6-2*n7-n8-2*
   n9-2*n10-2*n11-2*n12-n13-n14+17,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n7,1)+`Z'(n1,n2,n3,n4,
   n5,n6,n7,n8,n9,n10,n11,n12,1+n13,-1+n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n7,
   1));
* B14
id,ifmatch->no61,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?pos_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n7+1)*(`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9,n10,1+n11,n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,1+n11,n12,n13,n14)*rat(-1+n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,1+n9,n10,1+n11,n12,n13,n14)*rat(n9,1)+`Z'(n1,-1+n2,
   n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n7,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,1+n11,n12,n13,n14)*
   rat(-n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,
   n10,1+n11,n12,n13,n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,1+n4,-1+n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n4,1)+`Z'(n1,n2,n3,
   1+n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*
   rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9,n10,1+n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,
   1+n11,n12,n13,-1+n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,
   n5,n6,-1+n7,n8,n9,n10,1+n11,1+n12,-1+n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(
   6*ep+2*n1+2*n2+n4+2*n5+2*n6+n7+n8+n9+n10+2*n11+2*n12+n13-11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,-1+n14)
   *rat(-1+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n7+1,1));
*--#] stat : 
#endif
#endif
#if ( `il' >= 2 )
*--#[ n12 :
id,ifmatch->no61
	`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_
    ,n9?pos_,n10?pos_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
			-RAT(-2*ep-n2-n6-2*n7-n9-n11-n13+5,1)*(
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,1+n13,-1+n14)*rat(n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(2*ep+n2+n3+n6+2*n10+n12+n14-4,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,-1+n13,-1+n14)*rat(-n11,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,-1+n14)*rat(n11,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,n14)*rat(-2*ep-n2-n3-n6-2*n10-n12-n14+4,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,-1+n14)*rat(n11-n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,-1+n14)*rat(n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,1+n13,n14)*rat(-n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,-1+n13,n14)*rat(n11,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n11,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,1+n14)*rat(n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n11+n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n11,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)
         +`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)
         +`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(n3,1)
         +`Z'(n1,-1+n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*n6,1)
         +`Z'(n1,-1+n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(2*n6,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,1+n12,n13,n14)*rat(-n9,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,1+n12,n13,n14)*rat(n9,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,-1+n13,n14)*rat(-n11,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n11-n12-1,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,n11,1+n12,n13,n14)*rat(n9,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-n10,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,1+n12,n13,n14)*rat(-n9,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,-1+n13,n14)*rat(n11,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n11,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(4*ep+2*n2+n3+2*n6+2*n7+n9+3*n10+n12+n13+n14-10,1)
         +`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(n3,1)
         +`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(-n3,1)
         +`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-n10,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n10,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,1+n13,n14)*rat(-n13,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,-1+n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(2*ep+2*n7+n9+n13+n14-3,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,2+n12,-1+n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,-1+n11,1+n12,1+n13,n14)*rat(n13,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,1+n14)*rat(n14,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,1+n12,-1+n13,n14)*rat(-n11,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n11-n13,1)
         +`Z'(n1,n2,n3,n4,n5,-1+n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n7,1)
         +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,-1+n13,n14)*rat(n11,1)
         +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n7+n10-n11+n12+1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,1+n12,-1+n13,n14)*rat(-n11,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n11,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*ep-2*n2-n3-2*n9-n12-n13+4,1)
         +`Z'(n1,n2,n3,n4,n5,n6,1+n7,-1+n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(-2*n7,1)
         +`Z'(n1,n2,n3,n4,n5,n6,1+n7,n8,-1+n9,-1+n10,n11,1+n12,n13,n14)*rat(n7,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*ep-n3-2*n10-n12-n14+2,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,-1+n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-4*ep-2*n2-n3-2*n6-3*n7-n9-2*n10-n11-n13-n14+11,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,2+n12,-1+n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,-1+n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,-1+n13,1+n14)*rat(n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(2*ep+2*n2+2*n3+n9+n11+n13-5,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,2+n12,-1+n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,1+n14)*rat(n14,1)
         +`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)
         +`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)
         +`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)
         +`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)
         +`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)
         +`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)
         +`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)
         +`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)
         +`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(-n11,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)
         +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(n11,1)
         +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,n11,2+n12,n13,n14)*rat(n12+1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,-1+n10,n11,1+n12,n13,1+n14)*rat(n14,1)
         );
*--#] n12 : 
#endif
#if ( `il' >= 3 )
*--#[ n11 :
id,ifmatch->no61
	`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_
    ,n9?pos_,n10?pos_,n11?neg_,0,n13?neg0_,n14?neg0_) =
			-RAT(-8*ep-3*n1-2*n2-2*n6-4*n7-2*n8-4*n9-2*n11-2*n13+15,1)*(
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9,n10,1+n11,0,-1+n13,n14)*rat(n8,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+n8,n9,n10,n11,0,n13,n14)*rat(-n8,1)
         +`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,-1+n14)*rat(n3,1)
         +`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,0,n13,-1+n14)*rat(2*n3,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,0,n13,-1+n14)*rat(2*n10,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,-1+n13,n14)*rat(4*ep+2*n1+2*n2+n6+n7+n9+n11+n13-6,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,-1+n14)*rat(2*ep-n3+n7+n9-n10+2*n11+n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,0,1+n13,-1+n14)*rat(2*n13,1)
         +`Z'(n1,-1+n2,n3,1+n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(n4,1)
         +`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(-2*n4,1)
         +`Z'(n1,n2,n3,1+n4,n5,n6,-1+n7,n8,n9,n10,1+n11,0,n13,n14)*rat(-2*n4,1)
         +`Z'(n1,n2,n3,1+n4,n5,n6,n7,n8,-1+n9,n10,1+n11,0,n13,n14)*rat(2*n4,1)
         +`Z'(n1,n2,n3,1+n4,n5,n6,n7,n8,n9,-1+n10,1+n11,0,n13,n14)*rat(2*n4,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,0,n13,n14)*rat(-2*n10,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(-4*ep-2*n2-2*n3-n4-n7-n8-2*n9+n10-n11-n13-n14+7,1)
         +`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,1+n11,0,n13,n14)*rat(-n9,1)
         +`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(-n2,1)
         +`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,0,n13,n14)*rat(n2,1)
         +`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,1+n11,0,n13,n14)*rat(n9,1)
         +`Z'(-1+n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,1+n11,0,n13,n14)*rat(-n6,1)
         +`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,1+n11,0,n13,n14)*rat(n6,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,1+n13,-1+n14)*rat(-n13,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(-n2+n7+n8+n9-n11-2,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,2+n11,0,-1+n13,n14)*rat(n11+1,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,2+n11,0,n13,-1+n14)*rat(-n11-1,1)
         +`Z'(n1,-1+n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(-n4,1)
         +`Z'(n1,-1+n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,1+n11,0,n13,n14)*rat(-n6,1)
         +`Z'(n1,-1+n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,1+n11,0,n13,n14)*rat(n6,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,1+n11,0,n13,n14)*rat(2*n9,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,1+n13,-1+n14)*rat(-n13,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(-4*ep-n2-2*n3-n4-2*n5-n6-2*n10-n11-n14+8,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,2+n11,0,-1+n13,n14)*rat(n11+1,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,2+n11,0,n13,-1+n14)*rat(-n11-1,1)
         +`Z'(n1,1+n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(n2,1)
         +`Z'(n1,1+n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,1+n11,0,n13,n14)*rat(-n2,1)
         +`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(2*n4,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,n10,1+n11,0,n13,n14)*rat(2*n9,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,0,n13,n14)*rat(-2*n10,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,1+n9,n10,1+n11,0,n13,n14)*rat(-4*n9,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(12*ep+3*n2+4*n3+2*n4+4*n5+4*n6+4*n7+2*n9+6*n10+2*
         n11+2*n13+2*n14-27,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(2*ep+n4+2*n5+n6+n10-5,1)
         +`Z'(n1,n2,n3,1+n4,-1+n5,n6,-1+n7,n8,n9,n10,1+n11,0,n13,n14)*rat(2*n4,1)
         +`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,-1+n9,n10,1+n11,0,n13,n14)*rat(-2*n4,1)
         +`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,-1+n10,1+n11,0,n13,n14)*rat(-2*n4,1)
         +`Z'(n1,n2,n3,n4,-1+n5,1+n6,-1+n7,n8,n9,n10,1+n11,0,n13,n14)*rat(n6,1)
         +`Z'(n1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,-1+n10,1+n11,0,n13,n14)*rat(-n6,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,1+n10,1+n11,0,n13,n14)*rat(2*n10,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,0,1+n13,-1+n14)*rat(n13,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,0,n13,n14)*rat(2*ep+n2+2*n3+n4-n10+n11+n14-2,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,2+n11,0,-1+n13,n14)*rat(-n11-1,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,2+n11,0,n13,-1+n14)*rat(n11+1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,0,n13,n14)*rat(4*ep+2*n4+4*n5+2*n6-2*n9+2*n10+2*n11-6,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,0,n13,n14)*rat(6*ep+2*n1+4*n2+n6+n7+n8+4*n9+n11+2*n13-12,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,0,n13,n14)*rat(-12*ep-3*n2-4*n3-2*n4-4*n5-4*n6-4*n7-2*n9-4*n10-4*
         n11-2*n13-2*n14+23,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,0,n13,n14)*rat(-4*ep-2*n4-4*n5-2*n6-2*n10+10,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,1+n14)*rat(2*n14,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,2+n11,0,n13,n14)*rat(-2*n11-2,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,1+n13,n14)*rat(-2*n13,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,n13,1+n14)*rat(-6*n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,0,1+n13,n14)*rat(2*n13,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,0,n13,1+n14)*rat(4*n14,1)
         +`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,0,1+n13,n14)*rat(-2*n13,1)
         +`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,2+n11,0,n13,n14)*rat(-4*n11-4,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,2+n11,0,n13,n14)*rat(-n11-1,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,2+n11,0,n13,n14)*rat(n11+1,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,0,1+n13,n14)*rat(2*n13,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,2+n11,0,n13,n14)*rat(5*n11+5,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,2+n11,0,n13,n14)*rat(-n11-1,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,2+n11,0,n13,n14)*rat(3*n11+3,1)
         +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,2+n11,0,n13,n14)*rat(n11+1,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,0,1+n13,n14)*rat(2*n13,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,n10,2+n11,0,n13,n14)*rat(n11+1,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,-1+n8,n9,n10,2+n11,0,n13,n14)*rat(-n11-1,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,2+n11,0,n13,n14)*rat(-n11-1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,1+n11,0,1+n13,n14)*rat(-2*n13,1)
         );
*--#] n11 : 
#endif
#if ( `il' >= 4 )
*--#[ n13 :
id,ifmatch->no61
	`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_
    ,n9?pos_,n10?pos_,0,0,n13?neg_,n14?neg0_) =
			-RAT(4*ep+n1+n2+n6+2*n7+2*n8+2*n9+n13-8,1)*(
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,1+n13,-1+n14)*rat(n2,1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,1+n13,-1+n14)*rat(-n2,1)
         +`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,1+n13,-1+n14)*rat(-n3,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,0,0,1+n13,-1+n14)*rat(-n6,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,0,0,1+n13,-1+n14)*rat(n6,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,0,0,1+n13,-1+n14)*rat(-n10,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,0,0,1+n13,-1+n14)*rat(-4*ep-n1-n2+n3-n6-2*n7-2*n8-2*n9+n10-2*n13+5,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,0,0,2+n13,-2+n14)*rat(n13+1,1)
         +`Z'(1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,0,0,1+n13,n14)*rat(-n1,1)
         +`Z'(1+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,0,0,1+n13,n14)*rat(n1,1)
         +`Z'(1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,1+n13,n14)*rat(n1,1)
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,1+n13,n14)*rat(-n2,1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,1+n13,n14)*rat(n2,1)
         +`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,0,0,1+n13,n14)*rat(n4,1)
         +`Z'(n1,n2,n3,1+n4,n5,n6,-1+n7,n8,n9,n10,0,0,1+n13,n14)*rat(n4,1)
         +`Z'(n1,n2,n3,1+n4,n5,n6,n7,n8,-1+n9,n10,0,0,1+n13,n14)*rat(-n4,1)
         +`Z'(n1,n2,n3,1+n4,n5,n6,n7,n8,n9,-1+n10,0,0,1+n13,n14)*rat(-n4,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+n8,n9,n10,0,0,1+n13,n14)*rat(-n8,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,0,0,1+n13,n14)*rat(n10,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,1+n8,-1+n9,n10,0,0,1+n13,n14)*rat(n8,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,0,0,1+n13,n14)*rat(2*ep+n2+n8+2*n9-n10+n13-2,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,0,0,2+n13,-1+n14)*rat(-n13-1,1)
         +`Z'(-1+n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,0,0,1+n13,n14)*rat(n6,1)
         +`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,0,0,1+n13,n14)*rat(-n6,1)
         +`Z'(1+n1,n2,n3,n4,n5,n6,-1+n7,-1+n8,n9,n10,0,0,1+n13,n14)*rat(n1,1)
         +`Z'(1+n1,n2,n3,n4,n5,n6,n7,-2+n8,n9,n10,0,0,1+n13,n14)*rat(-n1,1)
         +`Z'(1+n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,0,0,1+n13,n14)*rat(-n1,1)
         +`Z'(n1,1+n2,-2+n3,n4,n5,n6,n7,n8,n9,n10,0,0,1+n13,n14)*rat(-n2,1)
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,0,0,1+n13,n14)*rat(-n2,1)
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,0,0,1+n13,n14)*rat(n2,1)
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,1+n13,n14)*rat(2*n2,1)
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,-1+n10,0,0,1+n13,n14)*rat(n2,1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,0,0,1+n13,n14)*rat(n2,1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,0,0,1+n13,n14)*rat(-n2,1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-2+n9,n10,0,0,1+n13,n14)*rat(-n2,1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,0,0,1+n13,n14)*rat(-n2,1)
         +`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,0,0,1+n13,n14)*rat(-n4,1)
         +`Z'(n1,n2,-1+n3,n4,n5,1+n6,-1+n7,n8,n9,n10,0,0,1+n13,n14)*rat(n6,1)
         +`Z'(n1,n2,-1+n3,n4,n5,1+n6,n7,n8,n9,-1+n10,0,0,1+n13,n14)*rat(-n6,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,1+n13,n14)*rat(-4*ep-n2-2*n3-n4-2*n5-n6-2*n10-n14+10,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,2+n13,-1+n14)*rat(-n13-1,1)
         +`Z'(n1,n2,n3,-1+n4,n5,1+n6,-1+n7,n8,n9,n10,0,0,1+n13,n14)*rat(n6,1)
         +`Z'(n1,n2,n3,-1+n4,n5,1+n6,n7,n8,n9,-1+n10,0,0,1+n13,n14)*rat(-n6,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,1+n10,0,0,1+n13,n14)*rat(n10,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,0,0,1+n13,n14)*rat(-n10+1,1)
         +`Z'(n1,n2,n3,1+n4,-1+n5,n6,-1+n7,n8,n9,n10,0,0,1+n13,n14)*rat(-n4,1)
         +`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,-1+n9,n10,0,0,1+n13,n14)*rat(n4,1)
         +`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,-1+n10,0,0,1+n13,n14)*rat(n4,1)
         +`Z'(n1,n2,n3,n4,-1+n5,1+n6,-1+n7,n8,n9,n10,0,0,1+n13,n14)*rat(-n6,1)
         +`Z'(n1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,-1+n10,0,0,1+n13,n14)*rat(n6,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,1+n10,0,0,1+n13,n14)*rat(-n10,1)
         +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,0,0,1+n13,n14)*rat(n10-1,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,-2+n7,n8,n9,n10,0,0,1+n13,n14)*rat(n6,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,-1+n7,-1+n8,n9,n10,0,0,1+n13,n14)*rat(-n6,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,-1+n7,n8,-1+n9,n10,0,0,1+n13,n14)*rat(-n6,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,-1+n10,0,0,1+n13,n14)*rat(-2*n6,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,n7,-1+n8,n9,-1+n10,0,0,1+n13,n14)*rat(n6,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,n7,n8,-1+n9,-1+n10,0,0,1+n13,n14)*rat(n6,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,n7,n8,n9,-2+n10,0,0,1+n13,n14)*rat(n6,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,0,0,1+n13,n14)*rat(-4*ep-n2-2*n3-n4-2*n5-n6+n8-2*n10-n14+7,1)
         +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,0,0,2+n13,-1+n14)*rat(-n13-1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,0,0,1+n13,n14)*rat(-2*ep-n2-n8-2*n9-n13+4,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,0,0,2+n13,-1+n14)*rat(n13+1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,1+n13,n14)*rat(4*ep+n2+2*n3+n4+2*n5+n6-n8+2*n10+n14-7,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,2+n13,-1+n14)*rat(n13+1,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,0,0,1+n13,n14)*rat(4*ep+n2+2*n3+n4+2*n5+n6+2*n10+n14-10,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,0,0,2+n13,-1+n14)*rat(n13+1,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,1+n13,1+n14)*rat(2*n14,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,2+n13,n14)*rat(2*n13+2,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,1+n13,1+n14)*rat(-2*n14,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,2+n13,n14)*rat(-2*n13-2,1)
         );
*--#] n13 : 
#endif
#if ( `il' >= 5 )
*--#[ n14 :
id,ifmatch->no61
	`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_
    ,n9?pos_,n10?pos_,0,0,0,n14?neg_) =
			-RAT(-4*ep-n2-2*n3-2*n4-n5-n6-2*n10-n14+8,1)*(
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,0,n14)*rat(-n2,1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,0,n14)*rat(n2,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,0,0,0,n14)*rat(n6,1)
         +`Z'(n1,n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,0,0,0,n14)*rat(-n6,1)
         +`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(n2,1)
         +`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,0,1+n14)*rat(-n2,1)
         +`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(n4,1)
         +`Z'(n1,n2,-1+n3,n4,1+n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(n5,1)
         +`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(n5,1)
         +`Z'(n1,n2,n3,1+n4,n5,n6,n7,n8,n9,-1+n10,0,0,0,1+n14)*rat(-n4,1)
         +`Z'(n1,n2,n3,n4,1+n5,n6,n7,n8,n9,-1+n10,0,0,0,1+n14)*rat(-n5,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(2*ep+n2+2*n3+n4+n14-3,1)
         +`Z'(n1,1+n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(-n2,1)
         +`Z'(n1,1+n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,0,0,0,1+n14)*rat(n2,1)
         +`Z'(n1,n2,-1+n3,-1+n4,1+n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(-n5,1)
         +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(-n4+1,1)
         +`Z'(n1,n2,n3,-2+n4,1+n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(-n5,1)
         +`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,n9,-1+n10,0,0,0,1+n14)*rat(n5,1)
         +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,0,0,0,1+n14)*rat(-2*ep-n2-2*n3-n4-n14+4,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,0,0,0,1+n14)*rat(n4-1,1)
         );
*--#] n14 : 
#endif
GoTo no62;
Label no61;
$ired = 1;
Label no62;
id	`Z'(?a,x?neg0_,?b,n11?,...,n14?) = `Y'(?a,x,?b,n11,...,n14);
ModuleOption,maximum,$ired;
.sort:Main no6 loop, pass `$pass';
#if ( `$ired' > 0 )
	#redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
.sort
*
*--#] Loop 1 : 
*--#[ Loop 2 :
#do il = 1,10
*
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering no6 loop 2-`il', pass `$pass' at `time_' sec."
#endif
.sort
#if ( `il' >= 1 )
*--#[ n1 :
id,ifmatch->no63
	`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_
	,n8?pos_,n9?pos_,n10?pos_,0,0,0,0) =
			-RAT(-n1+1,1)*(
         +`Z'(-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,0,0,0,0)*rat(-n2,1)
         +`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,0,0,0,0)*rat(n2,1)
         +`Z'(-1+n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,0,0,0,0)*rat(n6,1)
         +`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,0,0,0,0)*rat(-n6,1)
         +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,0,0,0,0)*rat(4*ep+n1+n2+n6+2*n7+2*n8+2*n9-9,1)
         +`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,0,0,0,0)*rat(n1-1,1)
	);
*--#] n1 : 
#endif
#if ( `il' >= 2 )
*--#[ n5 :
id,ifmatch->no63
	`Z'(1,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_
	,n8?pos_,n9?pos_,n10?pos_,0,0,0,0) =
			-RAT(n5-1,1)*(
         +`Z'(1,1+n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,n10,0,0,0,0)*rat(-n2,1)
         +`Z'(1,1+n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,0,0,0,0)*rat(n2,1)
         +`Z'(1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,0,0,0,0)*rat(-n5+1,1)
         +`Z'(1,n2,n3,n4,-1+n5,1+n6,-1+n7,n8,n9,n10,0,0,0,0)*rat(n6,1)
         +`Z'(1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,-1+n10,0,0,0,0)*rat(-n6,1)
         +`Z'(1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,0,0,0,0)*rat(-4*ep-n2-2*n3-2*n4-n5-n6-2*n10+9,1)
      );
*--#] n5 : 
#endif
#if ( `il' >= 3 )
*--#[ n9 :
id,ifmatch->no63
	`Z'(1,n2?pos_,n3?pos_,n4?pos_,1,n6?pos_,n7?pos_
	,n8?pos_,n9?{>1},n10?pos_,0,0,0,0) =
			-RAT(-2*ep*n9+2*ep-2*n2*n9+2*n2-2*n9^2+6*n9-4,n2)*(
         +`Z'(1,1+n2,-1+n3,n4,1,n6,n7,1+n8,-1+n9,n10,0,0,0,0)*rat(2*n8,1)
         +`Z'(1,1+n2,-1+n3,n4,1,n6,n7,n8,n9,n10,0,0,0,0)*rat(2*n9-2,1)
         +`Z'(1,1+n2,n3,n4,1,n6,n7,1+n8,-2+n9,n10,0,0,0,0)*rat(-2*n8,1)
         +`Z'(1,1+n2,n3,n4,1,n6,n7,n8,-1+n9,n10,0,0,0,0)*rat(-2*ep-2*n2-2*n9+4,1)
         +`Z'(1,n2,n3,n4,1,1+n6,-1+n7,1+n8,-1+n9,n10,0,0,0,0)*rat(-2*n6*n8,n2)
         +`Z'(1,n2,n3,n4,1,1+n6,n7,1+n8,-1+n9,-1+n10,0,0,0,0)*rat(2*n6*n8,n2)
         +`Z'(1,n2,n3,n4,1,n6,n7,1+n8,-1+n9,n10,0,0,0,0)*rat(-6*ep*n8-2*n2*n8-2*n6*n8-4*n7*n8-2*n8^2-4*n8*n9+14*n8,n2)
         +`Z'(1,1+n2,-1+n3,n4,1,n6,n7,n8,-1+n9,n10,0,0,0,0)*rat(-2*ep,1)
         +`Z'(1,1+n2,n3,n4,1,n6,n7,n8,-2+n9,n10,0,0,0,0)*rat(2*ep,1)
         +`Z'(1,n2,n3,n4,1,1+n6,-1+n7,n8,-1+n9,n10,0,0,0,0)*rat(2*ep*n6,n2)
         +`Z'(1,n2,n3,n4,1,1+n6,n7,n8,-1+n9,-1+n10,0,0,0,0)*rat(-2*ep*n6,n2)
         +`Z'(1,n2,n3,n4,1,n6,n7,n8,-1+n9,n10,0,0,0,0)*rat(8*ep^2+2*ep*n2+2*ep*n6+4*ep*n7+4*ep*n8+4*ep*n9-18*ep,n2)
         +`Z'(2,n2,n3,n4,1,n6,n7,-1+n8,-1+n9,n10,0,0,0,0)*rat(2*ep,n2)
      );
*--#] n9 : 
#endif
#if ( `il' >= 4 )
*--#[ n3 :
id,ifmatch->no63
	`Z'(1,n2?pos_,n3?{>1},n4?pos_,1,n6?pos_,n7?pos_
	,n8?pos_,1,n10?pos_,0,0,0,0) =
			-RAT(-2*ep*n3+2*ep-2*n2*n3+2*n2-2*n3^2+6*n3-4,n2)*(
         +`Z'(1,1+n2,-2+n3,1+n4,1,n6,n7,n8,1,n10,0,0,0,0)*rat(-2*n4,1)
         +`Z'(1,1+n2,-1+n3,n4,1,n6,n7,n8,1,n10,0,0,0,0)*rat(-2*ep-2*n2-2*n3+4,1)
         +`Z'(1,n2,-1+n3,1+n4,1,1+n6,-1+n7,n8,1,n10,0,0,0,0)*rat(2*n4*n6,n2)
         +`Z'(1,n2,-1+n3,1+n4,1,1+n6,n7,n8,1,-1+n10,0,0,0,0)*rat(-2*n4*n6,n2)
         +`Z'(1,n2,-1+n3,1+n4,1,n6,n7,n8,1,n10,0,0,0,0)*rat(-6*ep*n4-2*n2*n4-4*n3*n4-2*n4^2-2*n4*n6-4*n4*n10+14*n4,n2)
         +`Y'(1,1+n2,-1+n3,1+n4,1,n6,n7,n8,0,n10,0,0,0,0)*rat(2*n4,1)
         +`Y'(1,1+n2,n3,n4,1,n6,n7,n8,0,n10,0,0,0,0)*rat(2*n3-2,1)
         +`Z'(1,1+n2,-2+n3,n4,1,n6,n7,n8,1,n10,0,0,0,0)*rat(2*ep,1)
         +`Z'(1,n2,-1+n3,-1+n4,2,n6,n7,n8,1,n10,0,0,0,0)*rat(2*ep,n2)
         +`Z'(1,n2,-1+n3,n4,1,1+n6,-1+n7,n8,1,n10,0,0,0,0)*rat(-2*ep*n6,n2)
         +`Z'(1,n2,-1+n3,n4,1,1+n6,n7,n8,1,-1+n10,0,0,0,0)*rat(2*ep*n6,n2)
         +`Z'(1,n2,-1+n3,n4,1,n6,n7,n8,1,n10,0,0,0,0)*rat(8*ep^2+2*ep*n2+4*ep*n3+4*ep*n4+2*ep*n6+4*ep*n10-18*ep,n2)
         +`Y'(1,1+n2,-1+n3,n4,1,n6,n7,n8,0,n10,0,0,0,0)*rat(-2*ep,1)
      );
*--#] n3 : 
#endif
#if ( `il' >= 5 )
*--#[ n10 :
id,ifmatch->no63
	`Z'(1,n2?pos_,1,n4?pos_,1,n6?pos_,n7?pos_,n8?pos_,1,n10?{>1},0,0,0,0) =
			-RAT(-2*ep*n10+2*ep-2*n6*n10+2*n6-2*n10^2+6*n10-4,n6)*(
         +`Z'(1,n2,1,1+n4,1,1+n6,-1+n7,n8,1,-1+n10,0,0,0,0)*rat(2*n4,1)
         +`Z'(1,n2,1,1+n4,1,1+n6,n7,n8,1,-2+n10,0,0,0,0)*rat(-2*n4,1)
         +`Z'(1,n2,1,1+n4,1,n6,n7,n8,1,-1+n10,0,0,0,0)*rat(-6*ep*n4-2*n2*n4-2*n4^2-2*n4*n6-4*n4*n10+10*n4,n6)
         +`Z'(1,n2,1,n4,1,1+n6,-1+n7,n8,1,n10,0,0,0,0)*rat(2*n10-2,1)
         +`Z'(1,n2,1,n4,1,1+n6,n7,n8,1,-1+n10,0,0,0,0)*rat(-2*ep-2*n6-2*n10+4,1)
         +`Y'(1,1+n2,0,1+n4,1,n6,n7,n8,1,-1+n10,0,0,0,0)*rat(-2*n2*n4,n6)
         +`Y'(1,1+n2,1,1+n4,1,n6,n7,n8,0,-1+n10,0,0,0,0)*rat(2*n2*n4,n6)
         +`Z'(1,n2,1,-1+n4,2,n6,n7,n8,1,-1+n10,0,0,0,0)*rat(2*ep,n6)
         +`Z'(1,n2,1,n4,1,1+n6,-1+n7,n8,1,-1+n10,0,0,0,0)*rat(-2*ep,1)
         +`Z'(1,n2,1,n4,1,1+n6,n7,n8,1,-2+n10,0,0,0,0)*rat(2*ep,1)
         +`Z'(1,n2,1,n4,1,n6,n7,n8,1,-1+n10,0,0,0,0)*rat(8*ep^2+2*ep*n2+4*ep*n4+2*ep*n6+4*ep*n10-14*ep,n6)
         +`Y'(1,1+n2,0,n4,1,n6,n7,n8,1,-1+n10,0,0,0,0)*rat(2*ep*n2,n6)
         +`Y'(1,1+n2,1,n4,1,n6,n7,n8,0,-1+n10,0,0,0,0)*rat(-2*ep*n2,n6)
     );
*--#] n10 : 
#endif
#if ( `il' >= 6 )
*--#[ n7 :
id,ifmatch->no63
	`Z'(1,n2?pos_,1,n4?pos_,1,n6?pos_,n7?{>1},n8?pos_,1,1,0,0,0,0) =
			-RAT(-2*ep^2*n7+2*ep^2-ep*n2*n7+ep*n2-3*ep*n4*n7+3*ep*n4-2*ep*n6*n7+2*ep
         *n6+4*ep*n7-4*ep-n2*n4*n7+n2*n4-n2*n6*n7+n2*n6+n2*n7-n2-n4^2*n7+n4^2-n4*n6*n7+n4*n6+3*n4*n7-3*n4+2*n6*n7-2*n6-
         2*n7+2,ep*n6+n6^2)*(
         +`Z'(1,-1+n2,1,n4,1,1+n6,-1+n7,1+n8,1,1,0,0,0,0)*rat(-3*ep*n8-n2*n8-n6*n8-2*n7*n8-n8^2+5*n8,ep+n2-1)
         +`Z'(1,-1+n2,1,n4,1,2+n6,-2+n7,1+n8,1,1,0,0,0,0)*rat(-n6*n8-n8,ep+n2-1)
         +`Z'(1,n2,1,-1+n4,1,2+n6,-2+n7,n8,1,2,0,0,0,0)*rat(ep*n6+ep+n6^2+n6*n7+n7-1,ep^2+2*ep*n6+ep+n6^2+n6)
         +`Z'(1,n2,1,-1+n4,1,2+n6,-1+n7,n8,1,1,0,0,0,0)*rat(-ep*n6-ep-n6^2-n6*n7-n7+1,ep+n6)
         +`Z'(1,n2,1,n4,1,1+n6,-1+n7,n8,1,1,0,0,0,0)*rat(-2*ep^3-ep^2*n2-3*ep^2*n4-4*ep^2*n6+2*ep^2-ep*n2*n4-2*ep*n2*n6
         -ep*n4^2-4*ep*n4*n6-2*ep*n4*n7+4*ep*n4-2*ep*n6^2+4*ep*n6+2*ep*n7-2*ep-n2*n4*n6-n2*n4*n7+n2*n4-n2*n6^2+n2*n7-n2
         -n4^2*n6-n4^2*n7+n4^2-n4*n6^2+2*n4*n6+3*n4*n7-3*n4+2*n6^2-2*n7+2,ep^2+2*ep*n6+ep+n6^2+n6)
         +`Z'(1,n2,1,n4,1,2+n6,-2+n7,n8,1,1,0,0,0,0)*rat(ep*n4*n6+ep*n4-ep*n6-ep+n4*n6^2+n4*n6*n7+n4*n7-n4-n6^2-n6*n7-
         n7+1,ep^2+2*ep*n6+ep+n6^2+n6)
         +`Y'(1,-1+n2,1,n4,1,2+n6,-1+n7,1+n8,1,0,0,0,0,0)*rat(n6*n8+n8,ep+n2-1)
         +`Y'(1,1+n2,0,n4,1,1+n6,-1+n7,n8,1,1,0,0,0,0)*rat(-ep*n2*n4+ep*n2-n2*n4*n6-n2*n4*n7+n2*n4+n2*n6+n2*n7-n2,ep^2+
         2*ep*n6+ep+n6^2+n6)
         +`Y'(1,1+n2,0,n4,1,n6,n7,n8,1,1,0,0,0,0)*rat(-ep*n2*n7+ep*n2-n2*n4*n7+n2*n4-n2*n6*n7+n2*n6+n2*n7-n2,ep*n6+n6^2
         )
         +`Y'(1,1+n2,1,n4,1,1+n6,-1+n7,n8,0,1,0,0,0,0)*rat(ep*n2*n4-ep*n2+n2*n4*n6+n2*n4*n7-n2*n4-n2*n6-n2*n7+n2,ep^2+2
         *ep*n6+ep+n6^2+n6)
         +`Y'(1,1+n2,1,n4,1,n6,n7,n8,0,1,0,0,0,0)*rat(ep*n2*n7-ep*n2+n2*n4*n7-n2*n4+n2*n6*n7-n2*n6-n2*n7+n2,ep*n6+n6^2)
         +`Y'(1,n2,0,n4,1,1+n6,-1+n7,1+n8,1,1,0,0,0,0)*rat(n2*n8-n8,ep+n2-1)
         +`Y'(1,n2,0,n4,1,1+n6,-1+n7,n8,1,2,0,0,0,0)*rat(-1,1)
         +`Y'(1,n2,0,n4,1,1+n6,-1+n7,n8,2,1,0,0,0,0)*rat(-ep,ep+n2-1)
         +`Y'(1,n2,0,n4,1,1+n6,n7,n8,1,1,0,0,0,0)*rat(-n7+1,1)
         +`Y'(1,n2,0,n4,1,n6,n7,n8,1,2,0,0,0,0)*rat(-n7+1,n6)
         +`Y'(1,n2,1,n4,1,1+n6,-1+n7,1+n8,0,1,0,0,0,0)*rat(-n2*n8+n8,ep+n2-1)
         +`Y'(1,n2,1,n4,1,1+n6,n7,n8,1,0,0,0,0,0)*rat(-ep*n7+ep-n4*n7+n4-n6*n7+n6+n7-1,ep+n6)
         +`Y'(1,n2,1,n4,1,2+n6,-1+n7,n8,1,0,0,0,0,0)*rat(-ep*n4*n6-ep*n4+ep*n6+ep-n4*n6^2-n4*n6*n7-n4*n7+n4+n6^2+n6*n7+
         n7-1,ep^2+2*ep*n6+ep+n6^2+n6)
         +`Z'(1,-1+n2,1,n4,1,1+n6,-1+n7,n8,1,1,0,0,0,0)*rat(4*ep^2+ep*n2+ep*n6+2*ep*n7+2*ep*n8-7*ep,ep+n2-1)
         +`Z'(1,-1+n2,1,n4,1,2+n6,-2+n7,n8,1,1,0,0,0,0)*rat(ep*n6+ep,ep+n2-1)
         +`Z'(1,n2,1,-2+n4,2,1+n6,-1+n7,n8,1,1,0,0,0,0)*rat(ep^2+ep*n6+ep*n7-ep,ep^2+2*ep*n6+ep+n6^2+n6)
         +`Z'(1,n2,1,-2+n4,2,n6,n7,n8,1,1,0,0,0,0)*rat(ep*n7-ep,ep*n6+n6^2)
         +`Z'(1,n2,1,-1+n4,1,1+n6,-1+n7,n8,1,1,0,0,0,0)*rat(4*ep^3+ep^2*n2+2*ep^2*n4+5*ep^2*n6+3*ep^2*n7-7*ep^2+ep*n2*
         n6+ep*n2*n7-ep*n2+2*ep*n4*n6+2*ep*n4*n7-2*ep*n4+ep*n6^2-4*ep*n6-5*ep*n7+5*ep,ep^2+2*ep*n6+ep+n6^2+n6)
         +`Z'(1,n2,1,-1+n4,1,2+n6,-2+n7,n8,1,1,0,0,0,0)*rat(-ep^2*n6-ep^2-ep*n6^2-ep*n6*n7-ep*n7+ep,ep^2+2*ep*n6+ep+
         n6^2+n6)
         +`Z'(1,n2,1,-1+n4,1,n6,n7,n8,1,1,0,0,0,0)*rat(4*ep^2*n7-4*ep^2+ep*n2*n7-ep*n2+2*ep*n4*n7-2*ep*n4+ep*n6*n7-ep*
         n6-5*ep*n7+5*ep,ep*n6+n6^2)
         +`Z'(2,-1+n2,1,n4,1,1+n6,-1+n7,-1+n8,1,1,0,0,0,0)*rat(ep,ep+n2-1)
         +`Y'(1,-1+n2,1,n4,1,2+n6,-1+n7,n8,1,0,0,0,0,0)*rat(-ep*n6-ep,ep+n2-1)
         +`Y'(1,1+n2,0,-1+n4,1,1+n6,-1+n7,n8,1,1,0,0,0,0)*rat(ep^2*n2+ep*n2*n6+ep*n2*n7-ep*n2,ep^2+2*ep*n6+ep+n6^2+n6)
         +`Y'(1,1+n2,0,-1+n4,1,n6,n7,n8,1,1,0,0,0,0)*rat(ep*n2*n7-ep*n2,ep*n6+n6^2)
         +`Y'(1,1+n2,1,-1+n4,1,1+n6,-1+n7,n8,0,1,0,0,0,0)*rat(-ep^2*n2-ep*n2*n6-ep*n2*n7+ep*n2,ep^2+2*ep*n6+ep+n6^2+n6)
         +`Y'(1,1+n2,1,-1+n4,1,n6,n7,n8,0,1,0,0,0,0)*rat(-ep*n2*n7+ep*n2,ep*n6+n6^2)
         +`Y'(1,n2,0,n4,1,1+n6,-1+n7,n8,1,1,0,0,0,0)*rat(-ep*n2+ep,ep+n2-1)
         +`Y'(1,n2,1,-1+n4,1,1+n6,n7,n8,1,0,0,0,0,0)*rat(ep*n7-ep,ep+n6)
         +`Y'(1,n2,1,-1+n4,1,2+n6,-1+n7,n8,1,0,0,0,0,0)*rat(ep^2*n6+ep^2+ep*n6^2+ep*n6*n7+ep*n7-ep,ep^2+2*ep*n6+ep+n6^2
         +n6)
         +`Y'(1,n2,1,n4,1,1+n6,-1+n7,n8,0,1,0,0,0,0)*rat(ep*n2-ep,ep+n2-1)
      );
*--#] n7 : 
#endif
#if ( `il' >= 7 )
*--#[ n2 :
id,ifmatch->no63
	`Z'(1,n2?{>1},1,n4?pos_,1,n6?pos_,1,n8?pos_,1,1,0,0,0,0) =
			-RAT(2*ep^2+ep*n2+3*ep*n4+2*ep*n6-4*ep+n2*n4+n2*n6-n2+n4^2+n4*n6-3*n4-2*
         n6+2,ep+n6)*(
         +`Z'(1,-1+n2,1,n4,1,1+n6,1,n8,1,1,0,0,0,0)*rat(-2*ep^2*n6-2*ep*n2*n6-ep*n6^2-3*ep*n6*n8+5*ep*n6-n2*n6^2-n2*n6*
         n8+n2*n6-n6^2*n8+2*n6^2-n6*n8^2+3*n6*n8-2*n6,ep*n2-ep+n2^2-2*n2+1)
         +`Y'(1,-1+n2,1,n4,1,1+n6,0,1+n8,1,1,0,0,0,0)*rat(3*ep*n6*n8+n2*n6*n8+n6^2*n8+n6*n8^2-3*n6*n8,ep*n2-ep+n2^2-2*
         n2+1)
         +`Y'(1,-1+n2,1,n4,1,2+n6,-1,1+n8,1,1,0,0,0,0)*rat(n6^2*n8+n6*n8,ep*n2-ep+n2^2-2*n2+1)
         +`Y'(1,-1+n2,1,n4,1,2+n6,0,1+n8,1,0,0,0,0,0)*rat(-n6^2*n8-n6*n8,ep*n2-ep+n2^2-2*n2+1)
         +`Y'(1,-1+n2,1,n4,1,2+n6,0,n8,1,1,0,0,0,0)*rat(-ep*n6^2-ep*n6-n2*n6^2-n2*n6-n6^2*n8+2*n6^2-n6*n8+2*n6,ep*n2-ep
         +n2^2-2*n2+1)
         +`Y'(1,-1+n2,1,n4,1,2+n6,1,n8,1,0,0,0,0,0)*rat(ep*n6^2+ep*n6+n2*n6^2+n2*n6+n6^2*n8-2*n6^2+n6*n8-2*n6,ep*n2-ep+
         n2^2-2*n2+1)
         +`Y'(1,1+n2,0,n4,1,n6,1,n8,1,1,0,0,0,0)*rat(ep*n2+n2*n4+n2*n6-n2,ep+n6)
         +`Y'(1,1+n2,1,n4,1,n6,1,n8,0,1,0,0,0,0)*rat(-ep*n2-n2*n4-n2*n6+n2,ep+n6)
         +`Y'(1,n2,0,n4,1,1+n6,0,1+n8,1,1,0,0,0,0)*rat(-n6*n8,ep+n2-1)
         +`Y'(1,n2,0,n4,1,1+n6,0,n8,2,1,0,0,0,0)*rat(-n6,ep+n2-1)
         +`Y'(1,n2,0,n4,1,1+n6,1,-1+n8,2,1,0,0,0,0)*rat(n6,ep+n2-1)
         +`Y'(1,n2,0,n4,1,1+n6,1,n8,1,1,0,0,0,0)*rat(ep*n6+n2*n6+n6*n8-2*n6,ep+n2-1)
         +`Y'(1,n2,0,n4,1,n6,1,n8,1,2,0,0,0,0)*rat(1,1)
         +`Y'(1,n2,1,-1+n4,1,1+n6,0,n8,1,2,0,0,0,0)*rat(-n6,ep+n6)
         +`Y'(1,n2,1,n4,1,1+n6,0,1+n8,0,1,0,0,0,0)*rat(n6*n8,ep+n2-1)
         +`Y'(1,n2,1,n4,1,1+n6,0,n8,1,1,0,0,0,0)*rat(-n4*n6+n6,ep+n6)
         +`Y'(1,n2,1,n4,1,1+n6,1,n8,0,1,0,0,0,0)*rat(-ep*n6-n2*n6-n6*n8+2*n6,ep+n2-1)
         +`Y'(1,n2,1,n4,1,1+n6,1,n8,1,0,0,0,0,0)*rat(ep*n6+n4*n6+n6^2-n6,ep+n6)
         +`Z'(1,-1+n2,1,n4,1,1+n6,1,-1+n8,1,1,0,0,0,0)*rat(4*ep^2*n6+ep*n2*n6+ep*n6^2+2*ep*n6*n8-5*ep*n6,ep*n2-ep+n2^2-
         2*n2+1)
         +`Z'(1,n2,1,-2+n4,2,n6,1,n8,1,1,0,0,0,0)*rat(-ep,ep+n6)
         +`Z'(1,n2,1,-1+n4,1,n6,1,n8,1,1,0,0,0,0)*rat(-4*ep^2-ep*n2-2*ep*n4-ep*n6+5*ep,ep+n6)
         +`Z'(2,-1+n2,1,n4,1,1+n6,1,-2+n8,1,1,0,0,0,0)*rat(ep*n6,ep*n2-ep+n2^2-2*n2+1)
         +`Y'(1,-1+n2,1,n4,1,1+n6,0,n8,1,1,0,0,0,0)*rat(-4*ep^2*n6-ep*n2*n6-ep*n6^2-2*ep*n6*n8+5*ep*n6,ep*n2-ep+n2^2-2*
         n2+1)
         +`Y'(1,-1+n2,1,n4,1,2+n6,-1,n8,1,1,0,0,0,0)*rat(-ep*n6^2-ep*n6,ep*n2-ep+n2^2-2*n2+1)
         +`Y'(1,-1+n2,1,n4,1,2+n6,0,-1+n8,1,1,0,0,0,0)*rat(ep*n6^2+ep*n6,ep*n2-ep+n2^2-2*n2+1)
         +`Y'(1,-1+n2,1,n4,1,2+n6,0,n8,1,0,0,0,0,0)*rat(ep*n6^2+ep*n6,ep*n2-ep+n2^2-2*n2+1)
         +`Y'(1,-1+n2,1,n4,1,2+n6,1,-1+n8,1,0,0,0,0,0)*rat(-ep*n6^2-ep*n6,ep*n2-ep+n2^2-2*n2+1)
         +`Y'(1,1+n2,0,-1+n4,1,n6,1,n8,1,1,0,0,0,0)*rat(-ep*n2,ep+n6)
         +`Y'(1,1+n2,1,-1+n4,1,n6,1,n8,0,1,0,0,0,0)*rat(ep*n2,ep+n6)
         +`Y'(1,n2,0,n4,1,1+n6,0,n8,1,1,0,0,0,0)*rat(ep*n6,ep+n2-1)
         +`Y'(1,n2,0,n4,1,1+n6,1,-1+n8,1,1,0,0,0,0)*rat(-ep*n6,ep+n2-1)
         +`Y'(1,n2,1,-1+n4,1,1+n6,0,n8,1,1,0,0,0,0)*rat(ep*n6,ep+n6)
         +`Y'(1,n2,1,-1+n4,1,1+n6,1,n8,1,0,0,0,0,0)*rat(-ep*n6,ep+n6)
         +`Y'(1,n2,1,n4,1,1+n6,0,n8,0,1,0,0,0,0)*rat(-ep*n6,ep+n2-1)
         +`Y'(1,n2,1,n4,1,1+n6,1,-1+n8,0,1,0,0,0,0)*rat(ep*n6,ep+n2-1)
         +`Y'(2,-1+n2,1,n4,1,1+n6,0,-1+n8,1,1,0,0,0,0)*rat(-ep*n6,ep*n2-ep+n2^2-2*n2+1)
      );
*--#] n2 : 
#endif
#if ( `il' >= 8 )
*--#[ n4 :
id,ifmatch->no63
	`Z'(1,1,1,n4?{>1},1,n6?pos_,1,n8?pos_,1,1,0,0,0,0) =
		-RAT(6*ep^2*n4-6*ep^2+2*ep*n4^2+2*ep*n4*n6-4*ep*n4-2*ep*n6+2*ep,ep+n6)*(
         +`Z'(1,1,1,-1+n4,1,2+n6,1,-1+n8,1,1,0,0,0,0)*rat(-2*ep^2*n6^2-2*ep^2*n6-2*ep*n6^3-4*ep*n6^2-2*ep*n6,2*ep^2+3*
         ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,0,1,n4,1,1+n6,1,n8,1,1,0,0,0,0)*rat(-6*ep^3*n4*n6+6*ep^3*n6-11*ep^2*n4^2*n6-8*ep^2*n4*n6^2+9*ep^2*n4*n6
         *n8+13*ep^2*n4*n6+8*ep^2*n6^2-9*ep^2*n6*n8-2*ep^2*n6-6*ep*n4^3*n6-8*ep*n4^2*n6^2+3*ep*n4^2*n6*n8+15*ep*n4^2*n6
         -2*ep*n4*n6^3+6*ep*n4*n6^2*n8+10*ep*n4*n6^2+3*ep*n4*n6*n8^2-12*ep*n4*n6*n8-9*ep*n4*n6+2*ep*n6^3-6*ep*n6^2*n8-2
         *ep*n6^2-3*ep*n6*n8^2+9*ep*n6*n8-n4^4*n6-2*n4^3*n6^2+4*n4^3*n6-n4^2*n6^3+n4^2*n6^2*n8+5*n4^2*n6^2+n4^2*n6*n8^2
         -2*n4^2*n6*n8-5*n4^2*n6+n4*n6^3*n8+n4*n6^3+n4*n6^2*n8^2-4*n4*n6^2*n8-3*n4*n6^2-2*n4*n6*n8^2+4*n4*n6*n8+2*n4*n6
         -n6^3*n8-n6^2*n8^2+3*n6^2*n8+n6*n8^2-2*n6*n8,2*ep^3+3*ep^2*n4+2*ep^2*n6-3*ep^2+ep*n4^2+ep*n4*n6-2*ep*n4-ep*n6+
         ep)
         +`Y'(1,0,1,n4,1,2+n6,0,n8,1,1,0,0,0,0)*rat(2*ep^2*n4*n6^2+2*ep^2*n4*n6-2*ep^2*n6^2-2*ep^2*n6+3*ep*n4^2*n6^2+3*
         ep*n4^2*n6+2*ep*n4*n6^3+3*ep*n4*n6^2*n8-7*ep*n4*n6^2+3*ep*n4*n6*n8-9*ep*n4*n6-2*ep*n6^3-3*ep*n6^2*n8+4*ep*n6^2
         -3*ep*n6*n8+6*ep*n6+n4^3*n6^2+n4^3*n6+n4^2*n6^3+n4^2*n6^2*n8-3*n4^2*n6^2+n4^2*n6*n8-4*n4^2*n6+n4*n6^3*n8-3*n4*
         n6^3-n4*n6^2*n8+2*n4*n6^2-2*n4*n6*n8+5*n4*n6-n6^3*n8+2*n6^3+n6*n8-2*n6,2*ep^3+3*ep^2*n4+2*ep^2*n6-3*ep^2+ep*
         n4^2+ep*n4*n6-2*ep*n4-ep*n6+ep)
         +`Y'(1,0,1,n4,1,2+n6,1,n8,1,0,0,0,0,0)*rat(-2*ep^2*n4*n6^2-2*ep^2*n4*n6+2*ep^2*n6^2+2*ep^2*n6-3*ep*n4^2*n6^2-3
         *ep*n4^2*n6-2*ep*n4*n6^3-3*ep*n4*n6^2*n8+7*ep*n4*n6^2-3*ep*n4*n6*n8+9*ep*n4*n6+2*ep*n6^3+3*ep*n6^2*n8-4*ep*
         n6^2+3*ep*n6*n8-6*ep*n6-n4^3*n6^2-n4^3*n6-n4^2*n6^3-n4^2*n6^2*n8+3*n4^2*n6^2-n4^2*n6*n8+4*n4^2*n6-n4*n6^3*n8+3
         *n4*n6^3+n4*n6^2*n8-2*n4*n6^2+2*n4*n6*n8-5*n4*n6+n6^3*n8-2*n6^3-n6*n8+2*n6,2*ep^3+3*ep^2*n4+2*ep^2*n6-3*ep^2+
         ep*n4^2+ep*n4*n6-2*ep*n4-ep*n6+ep)
         +`Y'(1,1,0,n4,1,1+n6,1,-1+n8,1,2,0,0,0,0)*rat(3*ep*n4*n6-3*ep*n6+n4^2*n6+n4*n6^2-2*n4*n6-n6^2+n6,2*ep^2+3*ep*
         n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,1,0,n4,1,1+n6,1,-1+n8,2,1,0,0,0,0)*rat(3*ep*n4*n6-3*ep*n6+n4^2*n6+n4*n6^2-2*n4*n6-n6^2+n6,2*ep^2+3*ep*
         n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,1,0,n4,1,1+n6,2,-1+n8,1,1,0,0,0,0)*rat(3*ep*n4*n6-3*ep*n6+n4^2*n6+n4*n6^2-2*n4*n6-n6^2+n6,2*ep^2+3*ep*
         n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,1,0,n4,1,n6,2,-1+n8,1,2,0,0,0,0)*rat(3*ep*n4-3*ep+n4^2+n4*n6-2*n4-n6+1,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2
         +n4*n6-2*n4-n6+1)
         +`Y'(1,1,1,-1+n4,1,1+n6,0,n8,1,2,0,0,0,0)*rat(-3*ep*n6-n6^2,ep+n6)
         +`Y'(1,1,1,-1+n4,1,1+n6,1,n8,0,2,0,0,0,0)*rat(-n6,1)
         +`Y'(1,1,1,-1+n4,1,1+n6,2,n8,0,1,0,0,0,0)*rat(-n6,1)
         +`Y'(1,1,1,-1+n4,1,2+n6,0,-1+n8,1,2,0,0,0,0)*rat(2*ep*n6^2+2*ep*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4
         -n6+1)
         +`Y'(1,1,1,-1+n4,1,n6,2,n8,0,2,0,0,0,0)*rat(-1,1)
         +`Y'(1,1,1,n4,1,1+n6,0,n8,1,1,0,0,0,0)*rat(-2*ep*n4*n6+2*ep*n6,ep+n6)
         +`Y'(1,1,1,n4,1,1+n6,1,n8,1,0,0,0,0,0)*rat(2*ep*n4*n6-2*ep*n6,ep+n6)
         +`Y'(1,1,1,n4,1,1+n6,2,-1+n8,1,0,0,0,0,0)*rat(n4*n6-n6,2*ep+n4-1)
         +`Y'(1,1,1,n4,1,2+n6,0,-1+n8,1,1,0,0,0,0)*rat(2*ep*n4*n6^2+2*ep*n4*n6-2*ep*n6^2-2*ep*n6,2*ep^2+3*ep*n4+2*ep*n6
         -3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,1,1,n4,1,2+n6,1,-1+n8,1,0,0,0,0,0)*rat(-2*ep*n4*n6^2-2*ep*n4*n6+2*ep*n6^2+2*ep*n6,2*ep^2+3*ep*n4+2*ep*
         n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,1,2,-1+n4,1,1+n6,1,n8,0,1,0,0,0,0)*rat(-n6,1)
         +`Y'(1,2,0,-1+n4,1,n6,1,n8,1,2,0,0,0,0)*rat(1,1)
         +`Y'(1,2,0,n4,1,1+n6,1,-1+n8,1,1,0,0,0,0)*rat(-2*ep*n4*n6+2*ep*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-
         n6+1)
         +`Y'(1,2,0,n4,1,n6,1,n8,1,1,0,0,0,0)*rat(2*ep*n4-2*ep,ep+n6)
         +`Y'(1,2,0,n4,1,n6,2,-1+n8,1,1,0,0,0,0)*rat(n4-1,2*ep+n4-1)
         +`Y'(1,2,1,-1+n4,1,n6,1,n8,0,2,0,0,0,0)*rat(-1,1)
         +`Y'(1,2,1,n4,1,1+n6,1,-1+n8,0,1,0,0,0,0)*rat(2*ep*n4*n6-2*ep*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-
         n6+1)
         +`Y'(1,2,1,n4,1,n6,1,n8,0,1,0,0,0,0)*rat(-2*ep*n4+2*ep,ep+n6)
         +`Y'(1,2,1,n4,1,n6,2,-1+n8,0,1,0,0,0,0)*rat(-n4+1,2*ep+n4-1)
         +`Z'(1,1,1,-2+n4,2,1+n6,1,-1+n8,1,1,0,0,0,0)*rat(2*ep^2*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Z'(1,1,1,-2+n4,2,n6,1,n8,1,1,0,0,0,0)*rat(-2*ep^2,ep+n6)
         +`Z'(1,1,1,-2+n4,2,n6,2,-1+n8,1,1,0,0,0,0)*rat(2*ep^2,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Z'(1,1,1,-1+n4,1,1+n6,1,-1+n8,1,1,0,0,0,0)*rat(8*ep^3*n6+4*ep^2*n4*n6+2*ep^2*n6^2-8*ep^2*n6,2*ep^2+3*ep*n4+2
         *ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Z'(1,1,1,-1+n4,1,n6,1,n8,1,1,0,0,0,0)*rat(-8*ep^3-4*ep^2*n4-2*ep^2*n6+8*ep^2,ep+n6)
         +`Z'(1,1,1,-1+n4,1,n6,2,-1+n8,1,1,0,0,0,0)*rat(8*ep^3+4*ep^2*n4+2*ep^2*n6-8*ep^2,2*ep^2+3*ep*n4+2*ep*n6-3*ep+
         n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,0,1,-2+n4,2,1+n6,1,n8,1,1,0,0,0,0)*rat(n6,1)
         +`Y'(1,0,1,-1+n4,1,1+n6,1,n8,1,1,0,0,0,0)*rat(4*ep*n6+2*n4*n6+n6^2-4*n6,1)
         +`Y'(1,0,1,-1+n4,1,2+n6,0,n8,1,1,0,0,0,0)*rat(-n6^2-n6,1)
         +`Y'(1,0,1,-1+n4,1,2+n6,1,n8,1,0,0,0,0,0)*rat(n6^2+n6,1)
         +`Y'(1,0,1,n4,1,1+n6,1,-1+n8,1,1,0,0,0,0)*rat(-12*ep^2*n4*n6+12*ep^2*n6-4*ep*n4^2*n6-7*ep*n4*n6^2-6*ep*n4*n6*
         n8+20*ep*n4*n6+7*ep*n6^2+6*ep*n6*n8-16*ep*n6-n4^2*n6^2-2*n4^2*n6*n8+4*n4^2*n6-n4*n6^3-2*n4*n6^2*n8+6*n4*n6^2+4
         *n4*n6*n8-8*n4*n6+n6^3+2*n6^2*n8-5*n6^2-2*n6*n8+4*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,0,1,n4,1,2+n6,0,-1+n8,1,1,0,0,0,0)*rat(-3*ep*n4*n6^2-3*ep*n4*n6+3*ep*n6^2+3*ep*n6-n4^2*n6^2-n4^2*n6-n4*
         n6^3+n4*n6^2+2*n4*n6+n6^3-n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,0,1,n4,1,2+n6,1,-1+n8,1,0,0,0,0,0)*rat(3*ep*n4*n6^2+3*ep*n4*n6-3*ep*n6^2-3*ep*n6+n4^2*n6^2+n4^2*n6+n4*
         n6^3-n4*n6^2-2*n4*n6-n6^3+n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,1,1,-1+n4,1,1+n6,0,n8,1,1,0,0,0,0)*rat(2*ep^2*n6,ep+n6)
         +`Y'(1,1,1,-1+n4,1,1+n6,1,n8,1,0,0,0,0,0)*rat(-2*ep^2*n6,ep+n6)
         +`Y'(1,1,1,-1+n4,1,1+n6,2,-1+n8,1,0,0,0,0,0)*rat(2*ep^2*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,1,1,-1+n4,1,2+n6,0,-1+n8,1,1,0,0,0,0)*rat(-2*ep^2*n6^2-2*ep^2*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6
         -2*n4-n6+1)
         +`Y'(1,1,1,-1+n4,1,2+n6,1,-1+n8,1,0,0,0,0,0)*rat(2*ep^2*n6^2+2*ep^2*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-
         2*n4-n6+1)
         +`Y'(1,2,0,-1+n4,1,1+n6,1,-1+n8,1,1,0,0,0,0)*rat(2*ep^2*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,2,0,-1+n4,1,n6,1,n8,1,1,0,0,0,0)*rat(-2*ep^2,ep+n6)
         +`Y'(1,2,0,-1+n4,1,n6,2,-1+n8,1,1,0,0,0,0)*rat(2*ep^2,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,2,1,-1+n4,1,1+n6,1,-1+n8,0,1,0,0,0,0)*rat(-2*ep^2*n6,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(1,2,1,-1+n4,1,n6,1,n8,0,1,0,0,0,0)*rat(2*ep^2,ep+n6)
         +`Y'(1,2,1,-1+n4,1,n6,2,-1+n8,0,1,0,0,0,0)*rat(-2*ep^2,2*ep^2+3*ep*n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
         +`Y'(2,0,1,n4,1,1+n6,1,-2+n8,1,1,0,0,0,0)*rat(-3*ep*n4*n6+3*ep*n6-n4^2*n6-n4*n6^2+2*n4*n6+n6^2-n6,2*ep^2+3*ep*
         n4+2*ep*n6-3*ep+n4^2+n4*n6-2*n4-n6+1)
      );
*--#] n4 : 
#endif
#if ( `il' >= 9 )
*--#[ n8 :
id,ifmatch->no63
	`Z'(1,1,1,1,1,n6?pos_,1,n8?{>1},1,1,0,0,0,0) =
			-RAT(-4*ep*n8+4*ep-n6*n8+n6,1)*(
         +`Z'(1,1,1,1,1,1+n6,1,-1+n8,1,1,0,0,0,0)*rat(2*ep*n6+n6^2+n6,1)
         +`Y'(0,1,1,1,2,n6,1,n8,1,1,0,0,0,0)*rat(n8-1,1)
         +`Y'(1,0,1,2,1,1+n6,0,n8,1,1,0,0,0,0)*rat(-ep*n6*n8+ep*n6-n6*n8+n6,ep^2+ep*n6)
         +`Y'(1,0,1,2,1,1+n6,1,n8,1,0,0,0,0,0)*rat(ep*n6*n8-ep*n6+n6*n8-n6,ep^2+ep*n6)
         +`Y'(1,0,1,2,1,n6,1,n8,1,1,0,0,0,0)*rat(-3*ep^2*n8+3*ep^2-ep*n6*n8+ep*n6-ep*n8^2+ep-n6*n8+n6-n8^2+3*n8-2,ep^2+
         ep*n6)
         +`Y'(1,1,0,2,1,-1+n6,2,-1+n8,1,2,0,0,0,0)*rat(1,ep+n6)
         +`Y'(1,1,0,2,1,n6,1,-1+n8,1,2,0,0,0,0)*rat(n6-1,ep+n6)
         +`Y'(1,1,0,2,1,n6,1,-1+n8,2,1,0,0,0,0)*rat(-ep-1,ep+n6)
         +`Y'(1,1,0,2,1,n6,2,-1+n8,1,1,0,0,0,0)*rat(n6-1,ep+n6)
         +`Y'(1,1,1,1,1,1+n6,0,-1+n8,1,2,0,0,0,0)*rat(-ep*n6-n6,ep+n6)
         +`Y'(1,1,1,1,1,1+n6,0,n8,1,1,0,0,0,0)*rat(n6*n8-n6,1)
         +`Y'(1,1,1,1,1,1+n6,1,n8,1,0,0,0,0,0)*rat(-n6*n8+n6,1)
         +`Y'(1,1,1,1,1,1+n6,2,-1+n8,1,0,0,0,0,0)*rat(-n6,1)
         +`Y'(1,1,1,2,1,1+n6,0,-1+n8,1,1,0,0,0,0)*rat(-ep*n6-n6,ep+n6)
         +`Y'(1,1,1,2,1,1+n6,1,-1+n8,1,0,0,0,0,0)*rat(ep*n6+n6,ep+n6)
         +`Y'(1,1,1,2,1,n6,2,-1+n8,1,0,0,0,0,0)*rat(-1,1)
         +`Y'(1,2,0,1,1,n6,1,n8,1,1,0,0,0,0)*rat(-n8+1,1)
         +`Y'(1,2,0,2,1,n6,1,-1+n8,1,1,0,0,0,0)*rat(ep+1,ep+n6)
         +`Y'(1,2,1,1,1,n6,1,n8,0,1,0,0,0,0)*rat(n8-1,1)
         +`Y'(1,2,1,2,1,n6,1,-1+n8,0,1,0,0,0,0)*rat(-ep-1,ep+n6)
         +`Z'(1,1,1,1,1,-1+n6,2,-1+n8,1,1,0,0,0,0)*rat(8*ep^2+6*ep*n6-2*ep+n6^2-n6,ep+n6)
         +`Z'(1,1,1,1,1,n6,1,-1+n8,1,1,0,0,0,0)*rat(16*ep^3+32*ep^2*n6+4*ep^2*n8-16*ep^2+15*ep*n6^2+5*ep*n6*n8-18*ep*n6
         +2*ep+2*n6^3+n6^2*n8-4*n6^2+n6,ep+n6)
         +`Y'(1,0,1,0,2,n6,1,-1+n8,2,1,0,0,0,0)*rat(1,1)
         +`Y'(1,0,1,1,1,1+n6,0,-1+n8,2,1,0,0,0,0)*rat(-n6,1)
         +`Y'(1,0,1,1,1,1+n6,1,-1+n8,2,0,0,0,0,0)*rat(n6,1)
         +`Y'(1,0,1,1,1,n6,1,-1+n8,2,1,0,0,0,0)*rat(4*ep+n6-1,1)
         +`Y'(1,0,1,2,1,1+n6,0,-1+n8,1,1,0,0,0,0)*rat(ep*n6+n6,ep+n6)
         +`Y'(1,0,1,2,1,1+n6,1,-1+n8,1,0,0,0,0,0)*rat(-ep*n6-n6,ep+n6)
         +`Y'(1,0,1,2,1,n6,1,-1+n8,1,1,0,0,0,0)*rat(4*ep^2+ep*n6+2*ep*n8-ep+n6+2*n8-5,ep+n6)
         +`Y'(1,1,0,0,2,n6,1,-1+n8,2,1,0,0,0,0)*rat(-1,1)
         +`Y'(1,1,0,1,1,1+n6,0,-1+n8,2,1,0,0,0,0)*rat(n6,1)
         +`Y'(1,1,0,1,1,1+n6,1,-1+n8,2,0,0,0,0,0)*rat(-n6,1)
         +`Y'(1,1,0,1,1,n6,1,-1+n8,2,1,0,0,0,0)*rat(-4*ep-n6+2,1)
         +`Y'(1,1,1,-1,2,n6,1,n8,1,1,0,0,0,0)*rat(-n8+1,1)
         +`Y'(1,1,1,0,1,1+n6,0,n8,1,1,0,0,0,0)*rat(n6*n8-n6,1)
         +`Y'(1,1,1,0,1,1+n6,1,n8,1,0,0,0,0,0)*rat(-n6*n8+n6,1)
         +`Y'(1,1,1,0,1,n6,1,n8,1,1,0,0,0,0)*rat(-4*ep*n8+4*ep-n6*n8+n6+2*n8-2,1)
         +`Y'(1,1,1,0,2,-1+n6,2,-1+n8,1,1,0,0,0,0)*rat(2*ep+n6,ep+n6)
         +`Y'(1,1,1,0,2,n6,1,-1+n8,1,1,0,0,0,0)*rat(4*ep^2+7*ep*n6+ep*n8-4*ep+2*n6^2+n6*n8-3*n6,ep+n6)
         +`Y'(1,1,1,0,2,n6,2,-1+n8,1,0,0,0,0,0)*rat(-1,1)
         +`Y'(1,1,1,1,1,1+n6,0,-1+n8,1,1,0,0,0,0)*rat(-4*ep^2*n6-7*ep*n6^2-ep*n6*n8+4*ep*n6-2*n6^3-n6^2*n8+3*n6^2,ep+n6)
         +`Y'(1,1,1,1,1,1+n6,1,-1+n8,1,0,0,0,0,0)*rat(4*ep^2*n6+7*ep*n6^2+ep*n6*n8-3*ep*n6+2*n6^3+n6^2*n8-2*n6^2,ep+n6)
         +`Y'(1,1,1,1,1,1+n6,2,-1+n8,1,-1,0,0,0,0)*rat(-n6,1)
         +`Y'(1,1,1,1,1,n6,2,-1+n8,1,0,0,0,0,0)*rat(-4*ep^2-3*ep*n6+n6,ep+n6)
         +`Y'(1,2,-1,1,1,n6,1,-1+n8,2,1,0,0,0,0)*rat(-1,1)
         +`Y'(1,2,0,0,1,n6,1,n8,1,1,0,0,0,0)*rat(-n8+1,1)
         +`Y'(1,2,0,1,1,-1+n6,2,-1+n8,1,1,0,0,0,0)*rat(2*ep+n6,ep+n6)
         +`Y'(1,2,0,1,1,n6,1,-1+n8,1,1,0,0,0,0)*rat(4*ep^2+7*ep*n6+ep*n8-3*ep+2*n6^2+n6*n8-2*n6,ep+n6)
         +`Y'(1,2,0,1,1,n6,2,-1+n8,1,0,0,0,0,0)*rat(-1,1)
         +`Y'(1,2,1,0,1,n6,1,n8,0,1,0,0,0,0)*rat(n8-1,1)
         +`Y'(1,2,1,1,1,-1+n6,2,-1+n8,0,1,0,0,0,0)*rat(-2*ep-n6,ep+n6)
         +`Y'(1,2,1,1,1,n6,1,-1+n8,0,1,0,0,0,0)*rat(-4*ep^2-7*ep*n6-ep*n8+4*ep-2*n6^2-n6*n8+3*n6,ep+n6)
         +`Y'(1,2,1,1,1,n6,2,-1+n8,0,0,0,0,0,0)*rat(1,1)
         +`Y'(2,0,1,2,1,n6,1,-2+n8,1,1,0,0,0,0)*rat(ep+1,ep+n6)
      );
*--#] n8 : 
#endif
#if ( `il' >= 10 )
*--#[ n6 :
`SKIP4'
id,ifmatch->no63,`Z'(1,1,1,1,1,n6?{>1},1,1,1,1,0,0,0,0) =
		-RAT(12*ep^4*n6-12*ep^4+10*ep^3*n6^2-10*ep^3*n6+2*ep^2*n6^3-2*ep^2*n6^2,1)*(
         +`Y'(0,1,1,1,2,-1+n6,1,2,1,1,0,0,0,0)*rat(6*ep^3+2*ep^2*n6,1)
         +`Y'(1,0,1,1,1,1+n6,0,2,1,1,0,0,0,0)*rat(-4*ep^2*n6^2+4*ep^2*n6-5*ep*n6^3+10*ep*n6^2-5*ep*n6-n6^4+3*n6^3-3*
         n6^2+n6,1)
         +`Y'(1,0,1,1,1,1+n6,1,2,1,0,0,0,0,0)*rat(4*ep^2*n6^2-4*ep^2*n6+5*ep*n6^3-10*ep*n6^2+5*ep*n6+n6^4-3*n6^3+3*n6^2
         -n6,1)
         +`Y'(1,0,1,1,1,n6,1,2,1,1,0,0,0,0)*rat(-12*ep^3*n6+12*ep^3-19*ep^2*n6^2+34*ep^2*n6-15*ep^2-8*ep*n6^3+19*ep*
         n6^2-14*ep*n6+3*ep-n6^4+3*n6^3-3*n6^2+n6,1)
         +`Y'(1,0,1,2,1,-1+n6,1,2,1,1,0,0,0,0)*rat(90*ep^5+114*ep^4*n6+6*ep^4+55*ep^3*n6^2+32*ep^3*n6-57*ep^3+12*ep^2*
         n6^3+28*ep^2*n6^2-64*ep^2*n6+24*ep^2+ep*n6^4+9*ep*n6^3-24*ep*n6^2+17*ep*n6-3*ep+n6^4-3*n6^3+3*n6^2-n6,ep^2)
         +`Y'(1,0,1,2,1,1+n6,-1,2,1,1,0,0,0,0)*rat(-12*ep^3*n6^2+12*ep^3*n6-10*ep^2*n6^3+10*ep^2*n6^2-2*ep*n6^4-ep*n6^3
         +6*ep*n6^2-3*ep*n6-n6^4+2*n6^3-n6^2,ep^2)
         +`Y'(1,0,1,2,1,1+n6,0,2,1,0,0,0,0,0)*rat(24*ep^3*n6^2-24*ep^3*n6+20*ep^2*n6^3-20*ep^2*n6^2+4*ep*n6^4+2*ep*n6^3
         -12*ep*n6^2+6*ep*n6+2*n6^4-4*n6^3+2*n6^2,ep^2)
         +`Y'(1,0,1,2,1,1+n6,1,2,1,-1,0,0,0,0)*rat(-12*ep^3*n6^2+12*ep^3*n6-10*ep^2*n6^3+10*ep^2*n6^2-2*ep*n6^4-ep*n6^3
         +6*ep*n6^2-3*ep*n6-n6^4+2*n6^3-n6^2,ep^2)
         +`Y'(1,0,1,2,1,n6,0,2,1,1,0,0,0,0)*rat(-6*ep^4*n6+6*ep^4-14*ep^3*n6^2+50*ep^3*n6-36*ep^3-7*ep^2*n6^3+34*ep^2*
         n6^2-33*ep^2*n6+6*ep^2-ep*n6^4+6*ep*n6^3-4*ep*n6^2-4*ep*n6+3*ep+n6^3-2*n6^2+n6,ep^2)
         +`Y'(1,0,1,2,1,n6,1,2,1,0,0,0,0,0)*rat(6*ep^4*n6-6*ep^4+14*ep^3*n6^2-26*ep^3*n6+12*ep^3+7*ep^2*n6^3-14*ep^2*
         n6^2+13*ep^2*n6-6*ep^2+ep*n6^4-2*ep*n6^3+6*ep*n6^2-8*ep*n6+3*ep+n6^3-2*n6^2+n6,ep^2)
         +`Y'(1,1,0,1,1,-1+n6,2,1,1,2,0,0,0,0)*rat(-4*ep^3-5*ep^2*n6+5*ep^2-ep*n6^2+2*ep*n6-ep,1)
         +`Y'(1,1,0,1,1,n6,1,1,1,2,0,0,0,0)*rat(-4*ep^3*n6+4*ep^3-5*ep^2*n6^2+10*ep^2*n6-5*ep^2-ep*n6^3+3*ep*n6^2-3*ep*
         n6+ep,1)
         +`Y'(1,1,0,1,1,n6,1,1,2,1,0,0,0,0)*rat(-4*ep^3*n6+4*ep^3-5*ep^2*n6^2+10*ep^2*n6-5*ep^2-ep*n6^3+3*ep*n6^2-3*ep*
         n6+ep,1)
         +`Y'(1,1,0,1,1,n6,2,1,1,1,0,0,0,0)*rat(-4*ep^3*n6+4*ep^3-5*ep^2*n6^2+10*ep^2*n6-5*ep^2-ep*n6^3+3*ep*n6^2-3*ep*
         n6+ep,1)
         +`Y'(1,1,0,2,1,-2+n6,2,1,1,2,0,0,0,0)*rat(-12*ep^2-7*ep*n6+3*ep-n6^2+n6,1)
         +`Y'(1,1,0,2,1,-1+n6,1,1,1,2,0,0,0,0)*rat(-12*ep^2*n6+24*ep^2-7*ep*n6^2+17*ep*n6-6*ep-n6^3+3*n6^2-2*n6,1)
         +`Y'(1,1,0,2,1,-1+n6,1,1,2,1,0,0,0,0)*rat(-6*ep^3-14*ep^2*n6+24*ep^2-7*ep*n6^2+17*ep*n6-6*ep-n6^3+3*n6^2-2*n6,1)
         +`Y'(1,1,0,2,1,-1+n6,2,1,1,1,0,0,0,0)*rat(-12*ep^2*n6+24*ep^2-7*ep*n6^2+17*ep*n6-6*ep-n6^3+3*n6^2-2*n6,1)
         +`Y'(1,1,1,0,1,1+n6,0,1,1,2,0,0,0,0)*rat(4*ep^2*n6^2-4*ep^2*n6+ep*n6^3-2*ep*n6^2+ep*n6,1)
         +`Y'(1,1,1,0,1,1+n6,1,1,1,1,0,0,0,0)*rat(-4*ep^3*n6^2+4*ep^3*n6-5*ep^2*n6^3+6*ep^2*n6^2-ep^2*n6-ep*n6^4+2*ep*
         n6^3-ep*n6^2,1)
         +`Y'(1,1,1,1,1,1+n6,-1,2,1,1,0,0,0,0)*rat(-12*ep^2*n6^2+12*ep^2*n6-10*ep*n6^3+16*ep*n6^2-6*ep*n6-2*n6^4+4*n6^3
         -2*n6^2,1)
         +`Y'(1,1,1,1,1,1+n6,0,1,1,1,0,0,0,0)*rat(12*ep^3*n6^2-12*ep^3*n6+10*ep^2*n6^3-4*ep^2*n6^2-6*ep^2*n6+2*ep*n6^4+
         6*ep*n6^3-14*ep*n6^2+6*ep*n6+2*n6^4-4*n6^3+2*n6^2,1)
         +`Y'(1,1,1,1,1,1+n6,0,2,1,0,0,0,0,0)*rat(12*ep^2*n6^2-12*ep^2*n6+10*ep*n6^3-16*ep*n6^2+6*ep*n6+2*n6^4-4*n6^3+2
         *n6^2,1)
         +`Y'(1,1,1,1,1,1+n6,1,1,1,0,0,0,0,0)*rat(-12*ep^3*n6^2+12*ep^3*n6-10*ep^2*n6^3+4*ep^2*n6^2+6*ep^2*n6-2*ep*n6^4
         -6*ep*n6^3+14*ep*n6^2-6*ep*n6-2*n6^4+4*n6^3-2*n6^2,1)
         +`Y'(1,1,1,1,1,n6,0,1,1,2,0,0,0,0)*rat(6*ep^3*n6-6*ep^3+2*ep^2*n6^2+10*ep^2*n6-12*ep^2+7*ep*n6^2-10*ep*n6+3*ep
         +n6^3-2*n6^2+n6,1)
         +`Y'(1,1,1,1,1,n6,0,2,1,1,0,0,0,0)*rat(-22*ep^3*n6+22*ep^3-38*ep^2*n6^2+66*ep^2*n6-28*ep^2-16*ep*n6^3+38*ep*
         n6^2-28*ep*n6+6*ep-2*n6^4+6*n6^3-6*n6^2+2*n6,1)
         +`Y'(1,1,1,1,1,n6,1,2,1,0,0,0,0,0)*rat(-14*ep^3*n6+14*ep^3-4*ep^2*n6^2+6*ep^2*n6-2*ep^2,1)
         +`Y'(1,1,1,1,1,n6,2,1,1,0,0,0,0,0)*rat(-6*ep^3*n6+6*ep^3+ep^2*n6^2+6*ep^2*n6-7*ep^2+ep*n6^3+4*ep*n6^2-9*ep*n6+
         4*ep+n6^3-2*n6^2+n6,1)
         +`Y'(1,1,1,2,1,-1+n6,2,1,1,0,0,0,0,0)*rat(-6*ep^3-8*ep^2*n6+12*ep^2-5*ep*n6^2+13*ep*n6-6*ep-n6^3+3*n6^2-2*n6,1)
         +`Y'(1,1,1,2,1,n6,0,1,1,1,0,0,0,0)*rat(-6*ep^3*n6+6*ep^3-8*ep^2*n6^2+20*ep^2*n6-12*ep^2-2*ep*n6^3+6*ep*n6^2-4*
         ep*n6,1)
         +`Y'(1,1,1,2,1,n6,1,1,1,0,0,0,0,0)*rat(6*ep^3*n6-6*ep^3+8*ep^2*n6^2-20*ep^2*n6+12*ep^2+2*ep*n6^3-6*ep*n6^2+4*
         ep*n6,1)
         +`Y'(1,2,0,1,1,-1+n6,1,1,1,2,0,0,0,0)*rat(-12*ep^3-10*ep^2*n6-6*ep^2-2*ep*n6^2-8*ep*n6+6*ep-2*n6^2+2*n6,1)
         +`Y'(1,2,0,1,1,-1+n6,1,1,2,1,0,0,0,0)*rat(-12*ep^3-10*ep^2*n6-2*ep*n6^2-3*ep*n6+3*ep-n6^2+n6,1)
         +`Y'(1,2,0,1,1,-1+n6,1,2,1,1,0,0,0,0)*rat(-14*ep^3-4*ep^2*n6+2*ep^2,1)
         +`Y'(1,2,0,1,1,-1+n6,2,1,1,1,0,0,0,0)*rat(-4*ep^3-5*ep^2*n6+5*ep^2-ep*n6^2+2*ep*n6-ep,1)
         +`Y'(1,2,0,1,1,n6,0,1,2,1,0,0,0,0)*rat(12*ep^2*n6-12*ep^2+10*ep*n6^2-16*ep*n6+6*ep+2*n6^3-4*n6^2+2*n6,1)
         +`Y'(1,2,0,1,1,n6,0,2,1,1,0,0,0,0)*rat(12*ep^2*n6-12*ep^2+10*ep*n6^2-16*ep*n6+6*ep+2*n6^3-4*n6^2+2*n6,1)
         +`Y'(1,2,0,1,1,n6,1,0,2,1,0,0,0,0)*rat(-12*ep^2*n6+12*ep^2-10*ep*n6^2+16*ep*n6-6*ep-2*n6^3+4*n6^2-2*n6,1)
         +`Y'(1,2,0,1,1,n6,1,1,1,1,0,0,0,0)*rat(-12*ep^3*n6+12*ep^3-10*ep^2*n6^2+4*ep^2*n6+6*ep^2-2*ep*n6^3-6*ep*n6^2+
         14*ep*n6-6*ep-2*n6^3+4*n6^2-2*n6,1)
         +`Y'(1,2,0,2,1,-2+n6,2,1,1,1,0,0,0,0)*rat(-6*ep^2-5*ep*n6+3*ep-n6^2+n6,1)
         +`Y'(1,2,0,2,1,-1+n6,1,1,1,1,0,0,0,0)*rat(6*ep^3+8*ep^2*n6-12*ep^2+2*ep*n6^2-4*ep*n6,1)
         +`Y'(1,2,1,0,1,n6,0,1,1,2,0,0,0,0)*rat(12*ep^3*n6-12*ep^3+10*ep^2*n6^2-4*ep^2*n6-6*ep^2+2*ep*n6^3+6*ep*n6^2-14
         *ep*n6+6*ep+2*n6^3-4*n6^2+2*n6,ep+n6-1)
         +`Y'(1,2,1,1,1,-1+n6,1,2,0,1,0,0,0,0)*rat(14*ep^3+4*ep^2*n6-2*ep^2,1)
         +`Y'(1,2,1,1,1,-1+n6,2,1,0,1,0,0,0,0)*rat(4*ep^3+5*ep^2*n6-5*ep^2+ep*n6^2-2*ep*n6+ep,1)
         +`Y'(1,2,1,1,1,n6,0,2,0,1,0,0,0,0)*rat(-12*ep^2*n6+12*ep^2-10*ep*n6^2+16*ep*n6-6*ep-2*n6^3+4*n6^2-2*n6,1)
         +`Y'(1,2,1,1,1,n6,1,1,0,1,0,0,0,0)*rat(12*ep^3*n6-12*ep^3+10*ep^2*n6^2-4*ep^2*n6-6*ep^2+2*ep*n6^3+6*ep*n6^2-14
         *ep*n6+6*ep+2*n6^3-4*n6^2+2*n6,1)
         +`Y'(1,2,1,1,1,n6,1,1,1,0,0,0,0,0)*rat(-12*ep^3*n6+12*ep^3-10*ep^2*n6^2+4*ep^2*n6+6*ep^2-2*ep*n6^3-6*ep*n6^2+
         14*ep*n6-6*ep-2*n6^3+4*n6^2-2*n6,1)
         +`Y'(1,2,1,2,1,-2+n6,2,1,0,1,0,0,0,0)*rat(6*ep^2+5*ep*n6-3*ep+n6^2-n6,1)
         +`Y'(1,2,1,2,1,-1+n6,1,1,0,1,0,0,0,0)*rat(-6*ep^3-8*ep^2*n6+12*ep^2-2*ep*n6^2+4*ep*n6,1)
         +`Y'(1,2,2,1,1,-1+n6,1,1,0,1,0,0,0,0)*rat(-12*ep^3-10*ep^2*n6-2*ep*n6^2-3*ep*n6+3*ep-n6^2+n6,1)
         +`Y'(1,3,0,1,1,-1+n6,1,1,1,1,0,0,0,0)*rat(-24*ep^3-20*ep^2*n6-12*ep^2-4*ep*n6^2-16*ep*n6+12*ep-4*n6^2+4*n6,1)
         +`Y'(1,3,1,1,1,-1+n6,1,1,0,1,0,0,0,0)*rat(24*ep^3+20*ep^2*n6+12*ep^2+4*ep*n6^2+16*ep*n6-12*ep+4*n6^2-4*n6,1)
         +`Z'(1,1,1,1,1,-1+n6,1,1,1,1,0,0,0,0)*rat(64*ep^5+64*ep^4*n6-32*ep^4+20*ep^3*n6^2-24*ep^3*n6+4*ep^3+2*ep^2*
         n6^3-4*ep^2*n6^2+2*ep^2*n6,1)
         +`Y'(1,0,1,0,2,-1+n6,1,1,2,1,0,0,0,0)*rat(6*ep^3+2*ep^2*n6,1)
         +`Y'(1,0,1,0,2,-1+n6,1,2,1,1,0,0,0,0)*rat(-36*ep^4-42*ep^3*n6+12*ep^3-16*ep^2*n6^2+ep^2*n6+9*ep^2-2*ep*n6^3-4*
         ep*n6^2+9*ep*n6-3*ep-n6^3+2*n6^2-n6,ep)
         +`Y'(1,0,1,0,2,n6,0,2,1,1,0,0,0,0)*rat(-12*ep^3*n6+12*ep^3-10*ep^2*n6^2+10*ep^2*n6-2*ep*n6^3-ep*n6^2+6*ep*n6-3
         *ep-n6^3+2*n6^2-n6,ep)
         +`Y'(1,0,1,0,2,n6,1,2,1,0,0,0,0,0)*rat(12*ep^3*n6-12*ep^3+10*ep^2*n6^2-10*ep^2*n6+2*ep*n6^3+ep*n6^2-6*ep*n6+3*
         ep+n6^3-2*n6^2+n6,ep)
         +`Y'(1,0,1,1,1,-1+n6,1,1,2,1,0,0,0,0)*rat(24*ep^4+14*ep^3*n6-12*ep^3+2*ep^2*n6^2-4*ep^2*n6,1)
         +`Y'(1,0,1,1,1,-1+n6,1,2,1,1,0,0,0,0)*rat(-144*ep^5-204*ep^4*n6+120*ep^4-106*ep^3*n6^2+100*ep^3*n6+12*ep^3-24*
         ep^2*n6^3+17*ep^2*n6^2+43*ep^2*n6-30*ep^2-2*ep*n6^4-4*ep*n6^3+25*ep*n6^2-25*ep*n6+6*ep-n6^4+4*n6^3-5*n6^2+2*n6
         ,ep)
         +`Y'(1,0,1,1,1,1+n6,-1,2,1,1,0,0,0,0)*rat(12*ep^3*n6^2-12*ep^3*n6+10*ep^2*n6^3-10*ep^2*n6^2+2*ep*n6^4+ep*n6^3-
         6*ep*n6^2+3*ep*n6+n6^4-2*n6^3+n6^2,ep)
         +`Y'(1,0,1,1,1,1+n6,0,1,1,1,0,0,0,0)*rat(4*ep^3*n6^2-4*ep^3*n6+5*ep^2*n6^3-10*ep^2*n6^2+5*ep^2*n6+ep*n6^4-3*ep
         *n6^3+3*ep*n6^2-ep*n6,1)
         +`Y'(1,0,1,1,1,1+n6,0,2,1,0,0,0,0,0)*rat(-24*ep^3*n6^2+24*ep^3*n6-20*ep^2*n6^3+20*ep^2*n6^2-4*ep*n6^4-2*ep*
         n6^3+12*ep*n6^2-6*ep*n6-2*n6^4+4*n6^3-2*n6^2,ep)
         +`Y'(1,0,1,1,1,1+n6,1,1,1,0,0,0,0,0)*rat(-4*ep^3*n6^2+4*ep^3*n6-5*ep^2*n6^3+10*ep^2*n6^2-5*ep^2*n6-ep*n6^4+3*
         ep*n6^3-3*ep*n6^2+ep*n6,1)
         +`Y'(1,0,1,1,1,1+n6,1,2,1,-1,0,0,0,0)*rat(12*ep^3*n6^2-12*ep^3*n6+10*ep^2*n6^3-10*ep^2*n6^2+2*ep*n6^4+ep*n6^3-
         6*ep*n6^2+3*ep*n6+n6^4-2*n6^3+n6^2,ep)
         +`Y'(1,0,1,1,1,n6,0,1,2,1,0,0,0,0)*rat(-6*ep^3*n6+6*ep^3-2*ep^2*n6^2+2*ep^2*n6,1)
         +`Y'(1,0,1,1,1,n6,0,2,1,1,0,0,0,0)*rat(-12*ep^3*n6+12*ep^3-10*ep^2*n6^2+10*ep^2*n6-2*ep*n6^3-ep*n6^2+6*ep*n6-3
         *ep-n6^3+2*n6^2-n6,1)
         +`Y'(1,0,1,1,1,n6,1,1,1,1,0,0,0,0)*rat(16*ep^4*n6-16*ep^4+24*ep^3*n6^2-48*ep^3*n6+24*ep^3+9*ep^2*n6^3-27*ep^2*
         n6^2+27*ep^2*n6-9*ep^2+ep*n6^4-4*ep*n6^3+6*ep*n6^2-4*ep*n6+ep,1)
         +`Y'(1,0,1,1,1,n6,1,1,2,0,0,0,0,0)*rat(6*ep^3*n6-6*ep^3+2*ep^2*n6^2-2*ep^2*n6,1)
         +`Y'(1,0,1,1,1,n6,1,2,1,0,0,0,0,0)*rat(12*ep^4*n6-12*ep^4+10*ep^3*n6^2-34*ep^3*n6+24*ep^3+2*ep^2*n6^3-19*ep^2*
         n6^2+14*ep^2*n6+3*ep^2-3*ep*n6^3-4*ep*n6^2+13*ep*n6-6*ep-2*n6^3+4*n6^2-2*n6,ep)
         +`Y'(1,0,1,2,1,-1+n6,1,1,1,1,0,0,0,0)*rat(24*ep^4+62*ep^3*n6-108*ep^3+42*ep^2*n6^2-120*ep^2*n6+72*ep^2+11*ep*
         n6^3-43*ep*n6^2+48*ep*n6-12*ep+n6^4-5*n6^3+8*n6^2-4*n6,1)
         +`Y'(1,0,1,2,1,n6,0,1,1,1,0,0,0,0)*rat(6*ep^3*n6-6*ep^3+14*ep^2*n6^2-38*ep^2*n6+24*ep^2+7*ep*n6^3-24*ep*n6^2+
         23*ep*n6-6*ep+n6^4-4*n6^3+5*n6^2-2*n6,1)
         +`Y'(1,0,1,2,1,n6,1,1,1,0,0,0,0,0)*rat(-6*ep^3*n6+6*ep^3-14*ep^2*n6^2+38*ep^2*n6-24*ep^2-7*ep*n6^3+24*ep*n6^2-
         23*ep*n6+6*ep-n6^4+4*n6^3-5*n6^2+2*n6,1)
         +`Y'(1,0,2,1,1,-1+n6,1,1,1,1,0,0,0,0)*rat(48*ep^4+52*ep^3*n6-24*ep^3+18*ep^2*n6^2-8*ep^2*n6-12*ep^2+2*ep*n6^3+
         3*ep*n6^2-13*ep*n6+6*ep+n6^3-3*n6^2+2*n6,1)
         +`Y'(1,0,2,1,1,n6,0,1,1,1,0,0,0,0)*rat(12*ep^3*n6-12*ep^3+10*ep^2*n6^2-10*ep^2*n6+2*ep*n6^3+ep*n6^2-6*ep*n6+3*
         ep+n6^3-2*n6^2+n6,1)
         +`Y'(1,0,2,1,1,n6,1,1,1,0,0,0,0,0)*rat(-12*ep^3*n6+12*ep^3-10*ep^2*n6^2+10*ep^2*n6-2*ep*n6^3-ep*n6^2+6*ep*n6-3
         *ep-n6^3+2*n6^2-n6,1)
         +`Y'(1,1,0,0,2,-1+n6,1,1,2,1,0,0,0,0)*rat(-6*ep^3-2*ep^2*n6,1)
         +`Y'(1,1,0,1,1,-1+n6,1,1,2,1,0,0,0,0)*rat(-24*ep^4-14*ep^3*n6+18*ep^3-2*ep^2*n6^2+6*ep^2*n6,1)
         +`Y'(1,1,0,1,1,n6,0,1,2,1,0,0,0,0)*rat(6*ep^3*n6-6*ep^3+2*ep^2*n6^2-2*ep^2*n6,1)
         +`Y'(1,1,0,1,1,n6,1,1,2,0,0,0,0,0)*rat(-6*ep^3*n6+6*ep^3-2*ep^2*n6^2+2*ep^2*n6,1)
         +`Y'(1,1,1,-1,2,-1+n6,1,2,1,1,0,0,0,0)*rat(-6*ep^3-2*ep^2*n6,1)
         +`Y'(1,1,1,-1,2,-1+n6,2,1,1,1,0,0,0,0)*rat(4*ep^3+ep^2*n6-ep^2,1)
         +`Y'(1,1,1,-1,2,n6,1,1,1,1,0,0,0,0)*rat(4*ep^3*n6-4*ep^3+ep^2*n6^2-2*ep^2*n6+ep^2,1)
         +`Y'(1,1,1,0,1,-1+n6,1,2,1,1,0,0,0,0)*rat(-24*ep^4-14*ep^3*n6+18*ep^3-2*ep^2*n6^2+6*ep^2*n6,1)
         +`Y'(1,1,1,0,1,-1+n6,2,1,1,1,0,0,0,0)*rat(16*ep^4+8*ep^3*n6-16*ep^3+ep^2*n6^2-4*ep^2*n6+3*ep^2,1)
         +`Y'(1,1,1,0,1,1+n6,0,1,1,1,0,0,0,0)*rat(-4*ep^3*n6^2+4*ep^3*n6-ep^2*n6^3+2*ep^2*n6^2-ep^2*n6,1)
         +`Y'(1,1,1,0,1,1+n6,1,1,1,0,0,0,0,0)*rat(4*ep^3*n6^2-4*ep^3*n6+ep^2*n6^3-2*ep^2*n6^2+ep^2*n6,1)
         +`Y'(1,1,1,0,1,n6,0,2,1,1,0,0,0,0)*rat(6*ep^3*n6-6*ep^3+2*ep^2*n6^2-2*ep^2*n6,1)
         +`Y'(1,1,1,0,1,n6,1,1,1,1,0,0,0,0)*rat(16*ep^4*n6-16*ep^4+8*ep^3*n6^2-24*ep^3*n6+16*ep^3+ep^2*n6^3-5*ep^2*n6^2
         +7*ep^2*n6-3*ep^2,1)
         +`Y'(1,1,1,0,1,n6,1,2,1,0,0,0,0,0)*rat(-6*ep^3*n6+6*ep^3-2*ep^2*n6^2+2*ep^2*n6,1)
         +`Y'(1,1,1,0,1,n6,2,1,1,0,0,0,0,0)*rat(4*ep^3*n6-4*ep^3+ep^2*n6^2-2*ep^2*n6+ep^2,1)
         +`Y'(1,1,1,0,2,-1+n6,1,1,1,1,0,0,0,0)*rat(24*ep^4+14*ep^3*n6-6*ep^3+2*ep^2*n6^2-2*ep^2*n6,1)
         +`Y'(1,1,1,0,2,-1+n6,2,1,1,0,0,0,0,0)*rat(-6*ep^3-2*ep^2*n6,1)
         +`Y'(1,1,1,1,1,-1+n6,2,1,1,0,0,0,0,0)*rat(-24*ep^4-14*ep^3*n6+18*ep^3-2*ep^2*n6^2+6*ep^2*n6,1)
         +`Y'(1,1,1,1,1,1+n6,-1,1,1,1,0,0,0,0)*rat(12*ep^3*n6^2-12*ep^3*n6+10*ep^2*n6^3-16*ep^2*n6^2+6*ep^2*n6+2*ep*
         n6^4-4*ep*n6^3+2*ep*n6^2,1)
         +`Y'(1,1,1,1,1,1+n6,0,0,1,1,0,0,0,0)*rat(-12*ep^3*n6^2+12*ep^3*n6-10*ep^2*n6^3+16*ep^2*n6^2-6*ep^2*n6-2*ep*
         n6^4+4*ep*n6^3-2*ep*n6^2,1)
         +`Y'(1,1,1,1,1,1+n6,0,1,1,0,0,0,0,0)*rat(-12*ep^3*n6^2+12*ep^3*n6-10*ep^2*n6^3+16*ep^2*n6^2-6*ep^2*n6-2*ep*
         n6^4+4*ep*n6^3-2*ep*n6^2,1)
         +`Y'(1,1,1,1,1,1+n6,1,0,1,0,0,0,0,0)*rat(12*ep^3*n6^2-12*ep^3*n6+10*ep^2*n6^3-16*ep^2*n6^2+6*ep^2*n6+2*ep*n6^4
         -4*ep*n6^3+2*ep*n6^2,1)
         +`Y'(1,1,1,1,1,n6,0,1,1,1,0,0,0,0)*rat(16*ep^4*n6-16*ep^4+36*ep^3*n6^2-76*ep^3*n6+40*ep^3+16*ep^2*n6^3-48*ep^2
         *n6^2+44*ep^2*n6-12*ep^2+2*ep*n6^4-8*ep*n6^3+10*ep*n6^2-4*ep*n6,1)
         +`Y'(1,1,1,1,1,n6,1,0,1,1,0,0,0,0)*rat(-48*ep^4*n6+48*ep^4-52*ep^3*n6^2+100*ep^3*n6-48*ep^3-18*ep^2*n6^3+52*
         ep^2*n6^2-46*ep^2*n6+12*ep^2-2*ep*n6^4+8*ep*n6^3-10*ep*n6^2+4*ep*n6,1)
         +`Y'(1,1,1,1,1,n6,1,1,1,0,0,0,0,0)*rat(32*ep^4*n6-32*ep^4+16*ep^3*n6^2-18*ep^3*n6+2*ep^3+2*ep^2*n6^3-2*ep^2*
         n6^2,1)
         +`Y'(1,1,1,1,1,n6,2,1,1,-1,0,0,0,0)*rat(-6*ep^3*n6+6*ep^3-2*ep^2*n6^2+2*ep^2*n6,1)
         +`Y'(1,2,-1,1,1,-1+n6,1,1,2,1,0,0,0,0)*rat(-6*ep^3-2*ep^2*n6,1)
         +`Y'(1,2,0,0,1,-1+n6,1,2,1,1,0,0,0,0)*rat(-6*ep^3-2*ep^2*n6,1)
         +`Y'(1,2,0,0,1,-1+n6,2,1,1,1,0,0,0,0)*rat(4*ep^3+ep^2*n6-ep^2,1)
         +`Y'(1,2,0,0,1,n6,1,1,1,1,0,0,0,0)*rat(4*ep^3*n6-4*ep^3+ep^2*n6^2-2*ep^2*n6+ep^2,1)
         +`Y'(1,2,0,1,1,-1+n6,1,1,1,1,0,0,0,0)*rat(32*ep^4+16*ep^3*n6-2*ep^3+2*ep^2*n6^2,1)
         +`Y'(1,2,0,1,1,-1+n6,2,1,1,0,0,0,0,0)*rat(-6*ep^3-2*ep^2*n6,1)
         +`Y'(1,2,0,1,1,n6,0,1,1,1,0,0,0,0)*rat(-12*ep^3*n6+12*ep^3-10*ep^2*n6^2+16*ep^2*n6-6*ep^2-2*ep*n6^3+4*ep*n6^2-
         2*ep*n6,1)
         +`Y'(1,2,0,1,1,n6,1,0,1,1,0,0,0,0)*rat(12*ep^3*n6-12*ep^3+10*ep^2*n6^2-16*ep^2*n6+6*ep^2+2*ep*n6^3-4*ep*n6^2+2
         *ep*n6,1)
         +`Y'(1,2,1,-1,2,-1+n6,1,1,1,1,0,0,0,0)*rat(12*ep^4+10*ep^3*n6+6*ep^3+2*ep^2*n6^2+8*ep^2*n6-6*ep^2+2*ep*n6^2-2*
         ep*n6,ep+n6-1)
         +`Y'(1,2,1,0,1,-1+n6,1,1,1,1,0,0,0,0)*rat(48*ep^5+52*ep^4*n6+18*ep^3*n6^2+18*ep^3*n6-36*ep^3+2*ep^2*n6^3+12*
         ep^2*n6^2-30*ep^2*n6+12*ep^2+2*ep*n6^3-6*ep*n6^2+4*ep*n6,ep+n6-1)
         +`Y'(1,2,1,0,1,-1+n6,1,2,0,1,0,0,0,0)*rat(6*ep^3+2*ep^2*n6,1)
         +`Y'(1,2,1,0,1,-1+n6,2,1,0,1,0,0,0,0)*rat(-4*ep^3-ep^2*n6+ep^2,1)
         +`Y'(1,2,1,0,1,n6,0,1,1,1,0,0,0,0)*rat(-12*ep^4*n6+12*ep^4-10*ep^3*n6^2+4*ep^3*n6+6*ep^3-2*ep^2*n6^3-6*ep^2*
         n6^2+14*ep^2*n6-6*ep^2-2*ep*n6^3+4*ep*n6^2-2*ep*n6,ep+n6-1)
         +`Y'(1,2,1,0,1,n6,1,1,0,1,0,0,0,0)*rat(-4*ep^3*n6+4*ep^3-ep^2*n6^2+2*ep^2*n6-ep^2,1)
         +`Y'(1,2,1,0,1,n6,1,1,1,0,0,0,0,0)*rat(12*ep^4*n6-12*ep^4+10*ep^3*n6^2-4*ep^3*n6-6*ep^3+2*ep^2*n6^3+6*ep^2*
         n6^2-14*ep^2*n6+6*ep^2+2*ep*n6^3-4*ep*n6^2+2*ep*n6,ep+n6-1)
         +`Y'(1,2,1,1,1,-1+n6,1,1,0,1,0,0,0,0)*rat(-32*ep^4-16*ep^3*n6+8*ep^3-2*ep^2*n6^2+2*ep^2*n6,1)
         +`Y'(1,2,1,1,1,-1+n6,2,1,0,0,0,0,0,0)*rat(6*ep^3+2*ep^2*n6,1)
         +`Y'(1,2,1,1,1,n6,0,1,0,1,0,0,0,0)*rat(12*ep^3*n6-12*ep^3+10*ep^2*n6^2-16*ep^2*n6+6*ep^2+2*ep*n6^3-4*ep*n6^2+2
         *ep*n6,1)
         +`Y'(1,2,1,1,1,n6,1,0,0,1,0,0,0,0)*rat(-12*ep^3*n6+12*ep^3-10*ep^2*n6^2+16*ep^2*n6-6*ep^2-2*ep*n6^3+4*ep*n6^2-
         2*ep*n6,1)
         +`Y'(1,3,0,0,1,-1+n6,1,1,1,1,0,0,0,0)*rat(24*ep^4+20*ep^3*n6+12*ep^3+4*ep^2*n6^2+16*ep^2*n6-12*ep^2+4*ep*n6^2-
         4*ep*n6,ep+n6-1)
         +`Y'(1,3,1,0,1,-1+n6,1,1,0,1,0,0,0,0)*rat(-24*ep^4-20*ep^3*n6-12*ep^3-4*ep^2*n6^2-16*ep^2*n6+12*ep^2-4*ep*n6^2
         +4*ep*n6,ep+n6-1)
         +`Y'(2,0,1,1,1,n6,1,0,1,1,0,0,0,0)*rat(4*ep^3*n6-4*ep^3+5*ep^2*n6^2-10*ep^2*n6+5*ep^2+ep*n6^3-3*ep*n6^2+3*ep*
         n6-ep,1)
         +`Y'(2,0,1,2,1,-1+n6,1,0,1,1,0,0,0,0)*rat(6*ep^3+14*ep^2*n6-24*ep^2+7*ep*n6^2-17*ep*n6+6*ep+n6^3-3*n6^2+2*n6,1)
         +`Y'(2,0,2,1,1,-1+n6,1,0,1,1,0,0,0,0)*rat(12*ep^3+10*ep^2*n6+2*ep*n6^2+3*ep*n6-3*ep+n6^2-n6,1)
         +`Y'(2,1,1,1,1,-1+n6,1,0,1,1,0,0,0,0)*rat(-8*ep^4-2*ep^3*n6+2*ep^3,1)
         +`Y'(2,1,1,1,1,n6,0,0,1,1,0,0,0,0)*rat(12*ep^3*n6-12*ep^3+10*ep^2*n6^2-16*ep^2*n6+6*ep^2+2*ep*n6^3-4*ep*n6^2+2
         *ep*n6,1)
         +`Y'(2,1,1,1,1,n6,1,-1,1,1,0,0,0,0)*rat(-12*ep^3*n6+12*ep^3-10*ep^2*n6^2+16*ep^2*n6-6*ep^2-2*ep*n6^3+4*ep*n6^2
         -2*ep*n6,1)
        );
*--#] n6 : 
#endif
GoTo no64;
Label no63;
$ired = 1;
Label no64;
id	`Z'(?a,x?neg0_,?b,n11?,...,n14?) = `Y'(?a,x,?b,n11,...,n14);
*id	`Y'(?a) = 0;
ModuleOption,maximum,$ired;
*Print +f +s;
.sort:Main no6 loop, pass `$pass';
#if ( `$ired' > 0 )
	#redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
id	`Z'(1,1,1,1,1,1,1,1,1,1,0,0,0,0) = Master(no6);
*
*--#] Loop 2 : 
*
#break
*--#] no6 : 
*--#[ lala :
#case lala
*
*--#[ Loop :
#do il = 1,17
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering lala loop `il', pass `$pass' at `time_' sec."
#endif
#if(1==1)
*--#[ stat file :
#if ( `il' >= 1 )
*
*  Prework to make the reduction of the lala topology faster.
*
* B3
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,n8,1+n9,n10,1+n11,n12,n13,n14)*rat(-n9,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,1+n11,n12,n13,n14)*rat(-n6+1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n1,1)+`Z'(1+n1,n2,
   n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,1+n11,n12,n13,n14)*rat(n1,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,1+n9,n10,1+n11,n12,n13,
   n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,1+n11,n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,
   n2,n3,n4,n5,-1+n6,n7,-1+n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,-1+n10,1+n11,n12,1+
   n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,n13,-1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,
   n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n1+n9,1)+`Z'(
   n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,n12,n13,
   n14)*rat(-1+n6,1));
* B15
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n3+1)*(`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,1+n7,
   n8,n9,1+n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,
   n2,-1+n3,1+n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n4,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,1+n7,n8,n9,1+n10,n11,n12,
   n13,n14)*rat(n7,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,
   -1+n8,n9,1+n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,-1+n11,n12,n13,1+n14)*rat(n14,1)+
   `Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,-1+n13,n14)*rat(-n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,
   n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n4+n7,1)+`Z'(n1,n2,-1+n3,n4,n5,
   n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n3+1,1)+`Z'(
   n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n3+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,
   n14)*rat(-1+n3,1));
* B8
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n8)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,-1+n8,1+n9,1+n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,
   1+n9,1+n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n8+1,1)+`Z'(n1,n2,
   n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-1+n8,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,1+n10,n11,n12,n13,n14)*
   rat(-2*ep-2*n2-n8-n9-n11-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-1+n8,1));
* B20
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n8+1)*(`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,1+n7,-1+n8,
   n9,n10,1+n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,1+n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,n2,
   n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n8+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,n13,n14)*
   rat(-2*ep-2*n5-n7-n8-n10-n13+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n8,1));
* B6
id,ifmatch->lala1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(1+n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,1+n11,n12,n13,n14)*rat(-n1,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,-1+
   n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,n9,n10,1+n11,n12,
   n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,
   n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n7,1)+`Z'(
   n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+
   n12,n13,n14)*rat(n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,
   n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(6*ep+n1+2*n2+2*n5+2*n6+n7+2*n8+2*n9+n10+n11+n12+2*n13+n14-13,1)+`Z'(n1,-1+n2,n3,n4,
   n5,n6,n7,n8,n9,n10,2+n11,n12,n13,n14)*rat(-2*n11-2,1));
* B7
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n8)*(`Z'(-1+n1,1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-n2,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,-1+n8,
   1+n9,n10,n11,n12,1+n13,n14)*rat(n9,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,
   n14)*rat(-n8+1,1)+`Z'(n1,1+n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(n2,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,-1+
   n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n8+1,1)+`Z'(n1,
   n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,-1+n12,1+n13,
   1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-n2+n9,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+
   n9,n10,n11,n12,1+n13,n14)*rat(-1+n8,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n8,1));
* B16
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n5+1)*(`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9,1+n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,
   n9,1+n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,
   -1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,1+n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,
   n13,n14)*rat(-n4,1)+`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n4,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,
   n8,n9,1+n10,n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n13,1)+`Z'(
   n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,-1+n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,1+
   n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,
   n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-6*ep-2*n2-2*n3-n4-2*n5-2*n7-2*n8-n9-n10-n11-n12-n13-2*n14+13,1)+`Z'(n1,n2,n3,n4,-
   1+n5,n6,n7,n8,n9,2+n10,n11,n12,n13,n14)*rat(-2*n10-2,1));
* B19
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n8)*(`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(-n8+1,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5,1)+`Z'(n1,
   n2,n3,-1+n4,n5,n6,1+n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(n7,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,-1+n8,n9,1+n10,n11,n12,
   n13,1+n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,
   n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8+1,1)+`Z'(n1,n2,n3,n4,1+n5,n6,-1+n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(n5,1)+`Z'(
   n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n8,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,1+n10,n11,n12,n13,
   1+n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,-1+n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+
   n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n8,1));
* B1
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n9)*(`Z'(1+n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,1+n11,n12,n13,n14)*rat(n1,1)+`Z'(1+n1,n2,n3,n4,n5,n6,n7,n8,-1+
   n9,n10,1+n11,n12,n13,n14)*rat(-n1,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n9+1,1)+`Z'(n1,n2,-1+
   n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,
   n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n9,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,-
   1+n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,n2,
   n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,n12,n13,n14)*rat(2*ep+n1+2*n6+n9+n12+n13-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+
   n11,n12,n13,n14)*rat(-n9+1,1));
* B11
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n7)*(`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,-1+n7,
   n8,n9,1+n10,n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n7+1,1)+`Z'(n1,n2,
   n3,1+n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n4,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*
   rat(-1+n7,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,
   1+n10,-1+n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,
   n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-2*ep-2*n3-n4-n7-n12-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,1+n10,
   n11,n12,n13,n14)*rat(-1+n7,1));
* B9
id,ifmatch->lala1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n2+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n6,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,
   1+n9,n10,n11,n12,n13,1+n14)*rat(-n9,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,n13,1+n14)*rat(-n1,1)+`Z'(n1,-
   2+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,1+n14)*rat(n9,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,1+
   n14)*rat(n6,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,
   n8,n9,n10,1+n11,n12,n13,1+n14)*rat(-n11,1)+`Z'(n1,-1+n2,n3,-1+n4,1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5,1)+`Z'(
   n1,-1+n2,n3,-1+n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n6,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,n9,n10,n11,
   n12,n13,1+n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(n1,-1+n2,n3,-1+n4,
   n5,n6,n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(n11,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-
   2*n12,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,
   n10,n11,n12,n13,2+n14)*rat(-2*n14-2,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7,1)+`Z'(n1,-
   1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,
   n13,1+n14)*rat(n5,1)+`Z'(n1,-1+n2,n3,n4,n5,1+n6,n7,n8,n9,n10,n11,-1+n12,n13,1+n14)*rat(-n6,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-
   1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(
   n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+
   n13,1+n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-2*ep-2*n1-n5-n6+n7-n9-n12-n13+4,1)
);
#endif
#if ( `il' >= 2 )
* B17
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n5+1)*(`Z'(-1+n1,1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n2,1)+`Z'(-1+n1,n2,1+n3,n4,-1+n5,n6,n7,
   n8,n9,n10,n11,n12,1+n13,n14)*rat(-n3,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,1+n13,n14)*rat(n9,1)+`Z'(-1+
   n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n10,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,
   1+n13,n14)*rat(-n11,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-2*n12,1)+`Z'(-1+n1,n2,n3,n4,-1+
   n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,2+n13,n14)*rat(-2*
   n13-2,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,
   n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,1+n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(n2,1)+`Z'(n1,n2,1+
   n3,-1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,n4,-1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,
   n14)*rat(n3,1)+`Z'(n1,n2,1+n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,n14)*rat(-n3,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,1+n7,
   n8,n9,n10,n11,n12,1+n13,n14)*rat(-n7,1)+`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,n14)*rat(-n4,1)+`Z'(n1,
   n2,n3,n4,-2+n5,n6,1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,
   n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-
   1+n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n10,1)+`Z'(n1,
   n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,
   1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-2*ep-n2-n3-2*n4-n7+n9-n12-n14+4,1));
* B10
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n9)*(`Z'(n1,1+n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-n2,1)+`Z'(n1,1+n2,n3,n4,n5,n6,n7,-1+n8,-
   1+n9,n10,n11,n12,1+n13,n14)*rat(n2,1)+`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,n11,n12,1+n13,n14)*rat(-n2,1)+`Z'(n1,1+
   n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-n2,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+
   n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,
   n8,n9,n10,n11,n12,1+n13,n14)*rat(-n9+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n9+1,1)+`Z'(n1,
   n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,
   n14)*rat(-1+n9,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-
   1+n9,n10,1+n11,n12,1+n13,n14)*rat(-2*n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,
   n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(2*ep+n2+2*n8+n9+n11+n14-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+
   n10,n11,n12,1+n13,n14)*rat(-n9+1,1));
* B18
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n7)*(`Z'(n1,-1+n2,n3,n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,
   n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7+1,1)+`Z'(n1,n2,-
   1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,n4,1+n5,n6,-1+n7,-1+n8,n9,n10,n11,n12,n13,1+
   n14)*rat(n5,1)+`Z'(n1,n2,n3,n4,1+n5,n6,-1+n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(n5,1)+`Z'(n1,n2,n3,n4,1+n5,n6,-1+n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(-n5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,
   n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,-1+n11,n12,n13,
   1+n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(2*n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,
   n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(2*ep+n5+n7+2*n8+
   n10+n13-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+
   n11,n12,n13,1+n14)*rat(-1+n7,1));
* B4
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n9)*(`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,n11,n12,n13,1+n14)*rat(-n9+1,1)+`Z'(1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,-1+n12,n13,1+n14)*rat(-n1,1)+`Z'(n1,-1+
   n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n9,1)+`Z'(n1,n2,-1+n3,n4,n5,1+n6,n7,n8,-1+n9,n10,n11,n12,n13,1+
   n14)*rat(n6,1)+`Z'(n1,n2,n3,-1+n4,n5,1+n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-n6,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-
   1+n9,n10,n11,1+n12,n13,1+n14)*rat(-2*n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+
   `Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,1+n6,n7,n8,-1+n9,n10,n11,-1+
   n12,n13,1+n14)*rat(-n6,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,
   n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-2*ep-2*n1-n6-n9-n12-n13+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,n13,
   1+n14)*rat(-n9+1,1));
* B12
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n7)*(`Z'(-1+n1,n2,1+n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n3,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,n8,
   n9,n10,n11,1+n12,1+n13,n14)*rat(-2*n12,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(
   n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,
   n12,1+n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n3,1)+`Z'(n1,n2,1+n3,n4,n5,
   n6,-1+n7,n8,n9,n10,n11,-1+n12,1+n13,n14)*rat(-n3,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n7+1,
   1)+`Z'(n1,n2,n3,1+n4,n5,n6,-1+n7,n8,n9,n10,n11,-1+n12,1+n13,n14)*rat(-n4,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,
   n12,1+n13,n14)*rat(-1+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,
   n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-2*ep-n3-2*n4-n7-n12-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,
   1+n13,n14)*rat(-n7+1,1));
* B5
id,ifmatch->lala1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n1+1)*(`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,
   n8,1+n9,n10,n11,n12,1+n13,n14)*rat(n9,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(-
   1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,2+n13,n14)*rat(-2*n13-2,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,n8,1+n9,n10,n11,
   n12,1+n13,n14)*rat(n9,1)+`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,-1+n10,n11,n12,1+n13,n14)*rat(n6,1)+`Z'(-1+n1,n2,n3,n4,n5,1+
   n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(
   -1+n1,n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9,-1+n10,n11,
   n12,1+n13,n14)*rat(n9,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,
   n8,n9,-1+n10,n11,n12,2+n13,n14)*rat(-n13-1,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n12,1)+`Z'(-
   1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(2*ep+n1+2*n6+n9+n12+n13-4,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,
   n10,n11,n12,1+n13,n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n1+1,1));
* B13
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n3+1)*(`Z'(n1,-1+n2,-1+n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n7,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(-2*n14-2,1)+`Z'(
   n1,n2,-2+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(n1,n2,-2+n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+
   n14)*rat(n7,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,
   n10,n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,-1+
   n3,n4,n5,n6,1+n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,n12,n13,1+
   n14)*rat(-n7,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,-1+n9,
   n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,2+n14)*rat(1+n14,1)+`Z'(n1,n2,-
   1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*
   rat(2*ep+2*n3+n4+n7+n12+n14-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-n3+1,1));
#endif
#if ( `il' >= 3 )
*--------------------------Code from complexity 2-----------------------
* G141
id,ifmatch->lala1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,1-n1)*(`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,1+n7,-1+
   n8,1+n9,n10,n11,n12,n13,n14)*rat(n7*n9,-1+n8)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n12,1)+`Z'(
   -1+n1,n2,n3,n4,-1+n5,n6,1+n7,-1+n8,1+n9,n10,n11,n12,n13,n14)*rat(-n7*n9,-1+n8)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9,
   n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(-1+n1,n2,n3,n4,
   n5,n6,n7,-1+n8,1+n9,n10,n11,n12,n13,n14)*rat(-2*ep*n9-2*n5*n9-n7*n9-n8*n9-n9*n10-n9*n13+5*n9,-1+n8)+`Z'(-1+n1,n2,n3,n4,
   n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n12,
   1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(2*ep+n1+2*n6+n9+n12+n13-5,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,
   n8,n9,n10,n11,n12,n13,n14)*rat(-1+n1,1));
* G151
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n4)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,1+n7,-1+n8,1+n9,n10,n11,n12,n13,n14)*rat(-n7*n9,-1+n8)+`Z'(-1+n1,n2,n3,-1+n4,
   n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,-1+n8,1+n9,n10,n11,n12,n13,n14)*rat(
   n7*n9,-1+n8)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,1+n7,n8,
   n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,n3,-1+
   n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,1+n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(
   2*ep*n7+2*n2*n7+n7*n8+n7*n9+n7*n11+n7*n14-5*n7,-1+n8)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-
   n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,
   n12,n13,n14)*rat(-2*ep-2*n3-n4-n7-n12-n14+5,1));
* G146
id,ifmatch->lala1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,1-n2)*(`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,-1+
   n8,1+n9,n10,n11,n12,n13,n14)*rat(-n7*n9,-1+n8)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+
   `Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,-1+n8,1+n9,n10,
   n11,n12,n13,n14)*rat(n7*n9,-1+n8)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,
   n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-
   n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,n12,n13,n14)*rat(2*ep*n9+2*n5*n9+n7*n9+n8*n9+n9*n10+n9*n13-5*n9,-
   1+n8)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-2*n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,
   n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(4*ep+n2+2*n5+n7+2*n8+n9+n10+n11+
   n13+n14-9,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,1+n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(-n2*n7+n7,-1+n8)+`Z'(n1,n2,n3,n4,-1+n5,n6,1+
   n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(n2*n7-n7,-1+n8)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-1+n2,
   1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(2*ep*n2-2*ep+2*n2*n5+n2*n7+n2*n8+n2*n10+n2*n13-5*n2-2*n5-
   n7-n8-n10-n13+5,-1+n8));
* G156
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n5)*(`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,-1+n8,1+n9,n10,n11,n12,n13,n14)*rat(n7*n9,-1+n8)+`Z'(-1+n1,n2,n3,n4,-1+n5,
   n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(
   -1+n1,n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,n12,n13,n14)*rat(n5*n9-n9,-1+n8)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,-1+n8,1+
   n9,n10,n11,n12,n13,n14)*rat(-n7*n9,-1+n8)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,-
   1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,1+n9,n10,n11,n12,n13,
   n14)*rat(-n5*n9+n9,-1+n8)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n5+1,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,
   n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)+`Z'(
   n1,n2,n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep*n7-2*n2*n7-n7*n8-n7*n9-n7*n11-n7*n14+5*n7,-1+n8)+`Z'(
   n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-2*n10,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,
   n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-4*ep-2*n2-n5-n7-2*n8-n9-n10-n11-n13-n14+9,
   1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep*n5+2*ep-2*n2*n5+2*n2-n5*n8-n5*n9-n5*n11-n5*n14+5*
   n5+n8+n9+n11+n14-5,-1+n8));
* G2132
id,ifmatch->lala1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(n1,-n13+n3*n13)*(`Z'(-1+n1,n2,-1+n3,n4,n5,1+n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-n6*n12,n1)+`Z'(-1+n1,n2,-1+n3,
   n4,n5,n6,n7,n8,1+n9,1+n10,n11,1+n12,n13,n14)*rat(-n9*n12,n1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,1+n9,1+n10,n11,1+n12,
   n13,n14)*rat(n9*n12,n1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,1+n14)*rat(n13*n14,n1)+`Z'(n1,n2,-1+n3,-
   1+n4,n5,n6,1+n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(n7*n12,n1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,1+n7,n8,n9,1+n10,n11,n12,1+
   n13,n14)*rat(n7*n13,n1)+`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n4*n13,n1)+`Z'(n1,n2,-1+n3,1+
   n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n4*n13,n1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,1+n7,n8,n9,1+n10,n11,1+n12,
   n13,n14)*rat(-n7*n12,n1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,1+n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n7*n13,n1)+`Z'(n1,n2,-1+n3,
   n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,1+n14)*rat(n13*n14,n1)+`Z'(n1,n2,-1+n3,n4,n5,1+n6,n7,n8,n9,1+n10,n11,1+n12,
   n13,n14)*rat(n6*n12,n1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,1+n10,n11,n12,1+n13,1+n14)*rat(-n13*n14,n1)+`Z'(n1,n2,-1+n3,
   n4,n5,n6,n7,n8,n9,1+n10,-1+n11,n12,1+n13,1+n14)*rat(-n13*n14,n1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,
   n14)*rat(-2*n1*n12+n3*n12+2*n4*n12-n6*n12+n7*n12-n9*n12+n12*n14-n12,n1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,
   n12,1+n13,1+n14)*rat(n13*n14,n1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n4*n13-n7*n13,n1)+`Z'(n1,
   n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13*n14,n1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,
   n13,n14)*rat(n3*n12-n12,n1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n3*n13-n13,n1)+`Z'(n1,n2,n3,
   n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n3*n13-n13,n1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,
   n14)*rat(-n3*n13+n13,n1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-n3*n12+n12,n1));
*-----------------------------------------------------------------------
#endif
*--#] stat file : 
#endif
#if ( `il' >= 4 )
*--#[ n3 :
id,ifmatch->lala1,
	`Z'(n1?`POS1',n2?`POS2',n3?{>1},n4?`POS4',n5?`POS5',n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?`POS9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(1,(n3-1))*(
      +`Z'(n1,n2,1+(n3-1),-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-(n3-1),1)
      +`Z'(n1,n2,(n3-1),-1+n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)
      +`Z'(n1,n2,(n3-1),n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n7,1)
      +`Z'(n1,n2,(n3-1),n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep-(n3-1)-2*n4-n7-n12-n14+4,1)
      );
*--#] n3 : 
#endif
#if ( `il' >= 5 )
*--#[ n6 :
id,ifmatch->lala1,
	`Z'(n1?`POS1',n2?`POS2',1,n4?`POS4',n5?`POS5',n6?{>1},n7?`POS7',n8?`POS8'
		,n9?`POS9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(1,(n6-1))*(
      +`Z'(-1+n1,n2,1,n4,n5,1+(n6-1),n7,n8,n9,n10,n11,n12,n13,n14)*rat(-(n6-1),1)
      +`Z'(-1+n1,n2,1,n4,n5,(n6-1),n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)
      +`Z'(n1,-1+n2,1,n4,n5,(n6-1),n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)
      +`Z'(n1,n2,1,n4,n5,(n6-1),n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep-2*n1-(n6-1)-n9-n12-n13+4,1)
      );
*--#] n6 : 
#endif
#if ( `il' >= 6 )
*--#[ n11 :
id,ifmatch->lala1,
	`Z'(n1?`POS1',n2?`POS2',1,n4?`POS4',n5?`POS5',1,n7?`POS7',n8?`POS8'
		,n9?`POS9',n10?`NEG10',n11?neg_,n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(1,-2*ep-2*n1-n9-n12-n13+3)*(
      +`Z'(-1+n1,n2,1,n4,n5,2,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(1,1)
      +`Z'(1+n1,-1+n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(2*n1,1)
      +`Z'(1+n1,n2,1,n4,n5,1,n7,n8,-1+n9,n10,(1+n11),n12,n13,n14)*rat(-n1,1)
      +`Z'(n1,-1+n2,1,n4,n5,2,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(1,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,-1+n10,(1+n11),n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,n10,-1+(1+n11),n12,1+n13,n14)*rat(-n13,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),1+n12,n13,-1+n14)*rat(n12,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(-n9+1,1)
      +`Z'(n1,n2,1,n4,n5,2,n7,n8,-1+n9,n10,(1+n11),n12,n13,n14)*rat(-1,1)
      +`Z'(-2+n1,n2,1,n4,n5,1,n7,n8,1+n9,n10,(1+n11),n12,n13,n14)*rat(-n9,1)
      +`Z'(-2+n1,n2,1,n4,n5,2,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(-1,1)
      +`Z'(-1+n1,-1+n2,1,n4,n5,2,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(-1,1)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,-1+n10,(1+n11),n12,1+n13,n14)*rat(-2*n13,1)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,n10,-1+(1+n11),n12,1+n13,n14)*rat(n13,1)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),1+n12,n13,-1+n14)*rat(-n12,1)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(2*ep+2*n9+n12+n13-2,1)
      +`Z'(-1+n1,n2,1,n4,n5,2,n7,n8,-1+n9,n10,(1+n11),n12,n13,n14)*rat(1,1)
      +`Z'(n1,-2+n2,1,n4,n5,1,n7,n8,1+n9,n10,(1+n11),n12,n13,n14)*rat(n9,1)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,-1+n10,(1+n11),n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(-2*ep-2*n1-n9-n12-n13+3,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,-1+n9,n10,(1+n11),n12,n13,n14)*rat(n1-1,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),n12,1+n13,n14)*rat(n13,1)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),1+n12,n13,n14)*rat(n12,1)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),n12,1+n13,n14)*rat(-n13,1)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),1+n12,n13,n14)*rat(-n12,1)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,n10,(1+n11),n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,-1+n8,n9,n10,(1+n11),n12,1+n13,n14)*rat(-n13,1)
      +`Z'(-1+n1,-1+n2,1,n4,n5,1,n7,n8,n9,n10,(1+n11),n12,1+n13,n14)*rat(-n13,1)
      +`Z'(-1+n1,n2,1,-1+n4,n5,1,n7,n8,n9,n10,(1+n11),1+n12,n13,n14)*rat(n12,1)
      +`Z'(-1+n1,n2,1,n4,-1+n5,1,n7,n8,n9,n10,(1+n11),n12,1+n13,n14)*rat(-n13,1)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,-1+n8,n9,n10,(1+n11),n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,-1+n2,1,-1+n4,n5,1,n7,n8,n9,n10,(1+n11),1+n12,n13,n14)*rat(-n12,1)
      +`Y'(-1+n1,n2,1,n4,n5,0,n7,n8,1+n9,n10,(1+n11),n12,n13,n14)*rat(n9,1)
      +`Y'(1+n1,-1+n2,1,n4,n5,0,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(-n1,1)
      +`Y'(n1,-1+n2,1,n4,n5,0,n7,n8,1+n9,n10,(1+n11),n12,n13,n14)*rat(-n9,1)
      +`Y'(n1,n2,1,n4,n5,0,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(n1-1,1)
      +`Y'(-1+n1,n2,0,n4,n5,1,n7,n8,n9,n10,(1+n11),1+n12,n13,n14)*rat(-n12,1)
      +`Y'(n1,-1+n2,0,n4,n5,1,n7,n8,n9,n10,(1+n11),1+n12,n13,n14)*rat(n12,1)
      );
*--#] n11 : 
#endif
#if ( `il' >= 7 )
*--#[ n10 :
id,ifmatch->lala1,
	`Z'(n1?`POS1',n2?`POS2',1,n4?`POS4',n5?`POS5',1,n7?`POS7',n8?`POS8'
		,n9?`POS9',n10?neg_,0,n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(1,-2*ep-2*n4-n7-n12-n14+3)*(
      +`Z'(n1,n2,1,1+n4,-1+n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(-2*n4,1)
      +`Z'(n1,n2,1,1+n4,n5,1,-1+n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(n4,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,-1+(1+n10),0,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,(1+n10),-1,n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,(1+n10),0,1+n12,-1+n13,n14)*rat(-n12,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(n7-1,1)
      +`Z'(n1,n2,2,-1+n4,n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(-1,1)
      +`Z'(n1,n2,2,n4,-1+n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(-1,1)
      +`Z'(n1,n2,2,n4,n5,1,-1+n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(1,1)
      +`Z'(n1,n2,1,-2+n4,n5,1,1+n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(n7,1)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,-1+(1+n10),0,n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,(1+n10),-1,n12,n13,1+n14)*rat(-2*n14,1)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,(1+n10),0,1+n12,-1+n13,n14)*rat(n12,1)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(-2*ep-2*n7-n12-n14+2,1)
      +`Z'(n1,n2,1,n4,-2+n5,1,1+n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(-n7,1)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,(1+n10),-1,n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(2*ep+2*n4+n7+n12+n14-3,1)
      +`Z'(n1,n2,1,n4,n5,1,-1+n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(-n4+1,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,(1+n10),0,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,2,-2+n4,n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(1,1)
      +`Z'(n1,n2,2,-1+n4,-1+n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(1,1)
      +`Z'(n1,n2,2,-1+n4,n5,1,-1+n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(-1,1)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,(1+n10),0,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,(1+n10),0,1+n12,n13,n14)*rat(-n12,1)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,(1+n10),0,n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,(1+n10),0,1+n12,n13,n14)*rat(n12,1)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,(1+n10),0,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,-1+n8,n9,(1+n10),0,n12,n13,1+n14)*rat(n14,1)
      +`Z'(-1+n1,n2,1,-1+n4,n5,1,n7,n8,n9,(1+n10),0,1+n12,n13,n14)*rat(-n12,1)
      +`Z'(-1+n1,n2,1,n4,-1+n5,1,n7,n8,n9,(1+n10),0,1+n12,n13,n14)*rat(n12,1)
      +`Z'(n1,-1+n2,1,-1+n4,n5,1,n7,n8,n9,(1+n10),0,n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,n2,1,-1+n4,-1+n5,1,n7,n8,n9,(1+n10),0,n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,-1+n8,n9,(1+n10),0,n12,n13,1+n14)*rat(-n14,1)
      +`Y'(n1,n2,0,-1+n4,n5,1,1+n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(-n7,1)
      +`Y'(n1,n2,0,1+n4,-1+n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(n4,1)
      +`Y'(n1,n2,0,n4,-1+n5,1,1+n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(n7,1)
      +`Y'(n1,n2,0,n4,n5,1,n7,n8,n9,(1+n10),0,n12,n13,n14)*rat(-n4+1,1)
      +`Y'(n1,n2,1,-1+n4,n5,0,n7,n8,n9,(1+n10),0,1+n12,n13,n14)*rat(n12,1)
      +`Y'(n1,n2,1,n4,-1+n5,0,n7,n8,n9,(1+n10),0,1+n12,n13,n14)*rat(-n12,1)
      );
*--#] n10 : 
#endif
#if ( `il' >= 8 )
*--#[ n14 :
id,ifmatch->lala1,
	`Z'(n1?`POS1',n2?`POS2',1,n4?`POS4',n5?`POS5',1,n7?`POS7',n8?`POS8'
		,n9?`POS9',0,0,n12?`NEG12',n13?`NEG13',n14?neg_)
			 = -rat(1,4*ep+2*n4+2*n7+2*n12+2*(1+n14)-6)*(
      +`Z'(-1+n1,n2,1,n4,n5,2,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-2*ep-2*n2+2*n4-2*n8-n9+n12-n13+(1+n14)+4,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(1+n1,-1+n2,1,n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-6*ep*n1-2*n1^2-4*n1*n2+4*n1*n4-4*n1*n8-3*n1*n9+n1*n12-3*
      n1*n13+2*n1*(1+n14)+11*n1,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(1+n1,n2,1,n4,n5,1,n7,n8,-1+n9,0,0,n12,n13,(1+n14))*rat(2*ep*n1+2*n1*n2-2*n1*n4+2*n1*n8+n1*n9-n1*n12+n1*n13-n1*
      (1+n14)-4*n1,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,-1+n2,1,n4,n5,2,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-2*ep-2*n2+2*n4-2*n8-n9+n12-n13+(1+n14)+4,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,n2,1,1+n4,-1+n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-4*ep*n4-2*n4*n5-2*n4*n7-4*n4*n8-2*n4*n13+10*n4,2*ep+2*n4
      +n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,1+n4,n5,1,-1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(2*ep*n4+n4*n5+n4*n7+2*n4*n8+n4*n13-5*n4,2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,-1,0,n12,1+n13,(1+n14))*rat(2*n1*n13-2*n2*n13+2*n4*n13-2*n8*n13+2*n12*n13+n13*(1+n14)+n13,2
      *ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,-1,0,n12,n13,1+(1+n14))*rat(-2*ep*(1+n14)-n5*(1+n14)-n7*(1+n14)-2*n8*(1+n14)-n13*(1+n14)+5*(1+n14),2*ep+2*n4+n7
      +n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,0,-1,n12,1+n13,(1+n14))*rat(-2*n1*n13+2*n2*n13-2*n4*n13+2*n8*n13-2*n12*n13-n13*(1+n14)-n13,
      2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,0,-1,n12,n13,1+(1+n14))*rat(2*ep*(1+n14)+n5*(1+n14)+n7*(1+n14)+2*n8*(1+n14)+n13*(1+n14)-5*(1+n14),2*ep+2*n4+n7+
      n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,0,0,1+n12,-1+n13,(1+n14))*rat(-2*ep*n12-n5*n12-n7*n12-2*n8*n12-n12*n13+5*n12,2*ep+2*n4+
      n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,0,0,1+n12,n13,-1+(1+n14))*rat(2*n1*n12-2*n2*n12+2*n4*n12-2*n8*n12+2*n12^2+n12*(1+n14)+n12,2
      *ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(16*ep^3+20*ep^2*n1+8*ep^2*n2+12*ep^2*n4+16*ep^2*n7+8*ep^2*n8+
      16*ep^2*n9+16*ep^2*n12+12*ep^2*n13+12*ep^2*(1+n14)-88*ep^2+4*ep*n1^2+8*ep*n1*n2+16*ep*n1*n4+18*ep*n1*n7+8*ep*n1*n8+6*
      ep*n1*n9+12*ep*n1*n12+6*ep*n1*n13+14*ep*n1*(1+n14)-72*ep*n1+8*ep*n2*n4+4*ep*n2*n7+8*ep*n2*n9+8*ep*n2*n12+4*ep*n2*n13+
      4*ep*n2*(1+n14)-28*ep*n2-4*ep*n4^2+2*ep*n4*n7+8*ep*n4*n8+10*ep*n4*n9+4*ep*n4*n12+10*ep*n4*n13+2*ep*n4*(1+n14)-44*ep*n4+2*
      ep*n5*n7-2*ep*n5+4*ep*n7^2+8*ep*n7*n8+12*ep*n7*n9+10*ep*n7*n12+12*ep*n7*n13+4*ep*n7*(1+n14)-60*ep*n7+8*ep*n8*n9+8*ep*
      n8*n12+4*ep*n8*n13+4*ep*n8*(1+n14)-32*ep*n8+4*ep*n9^2+8*ep*n9*n12+6*ep*n9*n13+8*ep*n9*(1+n14)-58*ep*n9+4*ep*n12^2+8*ep*
      n12*n13+8*ep*n12*(1+n14)-58*ep*n12+2*ep*n13^2+8*ep*n13*(1+n14)-46*ep*n13+2*ep*(1+n14)^2-40*ep*(1+n14)+162*ep+4*n1^2*n4+2*n1^2*n7+
      2*n1^2*n12+2*n1^2*(1+n14)-6*n1^2+8*n1*n2*n4+4*n1*n2*n7+4*n1*n2*n12+4*n1*n2*(1+n14)-12*n1*n2-4*n1*n4^2+2*n1*n4*n7+8*n1*n4*
      n8+6*n1*n4*n9+6*n1*n4*n13+2*n1*n4*(1+n14)-32*n1*n4+2*n1*n5*n7-2*n1*n5+4*n1*n7^2+8*n1*n7*n8+3*n1*n7*n9+3*n1*n7*n12+5*
      n1*n7*n13+4*n1*n7*(1+n14)-37*n1*n7+4*n1*n8*n12+4*n1*n8*(1+n14)-16*n1*n8+3*n1*n9*n12+3*n1*n9*(1+n14)-9*n1*n9+n1*n12^2+3*n1*n12
      *n13+3*n1*n12*(1+n14)-22*n1*n12+3*n1*n13*(1+n14)-11*n1*n13+2*n1*(1+n14)^2-25*n1*(1+n14)+67*n1+8*n2*n4*n9+4*n2*n4*n12+4*n2*n4*n13-
      16*n2*n4+4*n2*n7*n9+2*n2*n7*n12+2*n2*n7*n13-8*n2*n7+4*n2*n9*n12+4*n2*n9*(1+n14)-12*n2*n9+2*n2*n12^2+2*n2*n12*n13+2*n2
      *n12*(1+n14)-14*n2*n12+2*n2*n13*(1+n14)-6*n2*n13-8*n2*(1+n14)+24*n2-6*n4^2*n9-2*n4^2*n12-2*n4^2*n13+10*n4^2-n4*n7*n9+n4*n7*
      n12+n4*n7*n13-n4*n7+8*n4*n8*n9+4*n4*n8*n12+4*n4*n8*n13-16*n4*n8+4*n4*n9^2-3*n4*n9*n12+6*n4*n9*n13-3*n4*n9*(1+n14)-23*
      n4*n9-n4*n12^2+n4*n12*n13+n4*n12*(1+n14)-6*n4*n12+2*n4*n13^2+n4*n13*(1+n14)-21*n4*n13+n4*(1+n14)+41*n4+n5*n7*n9+n5*n7*n12+n5*
      n7*n13-3*n5*n7-n5*n9-n5*n12-n5*n13+3*n5+2*n7^2*n9+2*n7^2*n12+2*n7^2*n13-6*n7^2+6*n7*n8*n9+4*n7*n8*n12+4*n7*n8*n13
      -14*n7*n8+2*n7*n9^2+n7*n9*n12+4*n7*n9*n13+n7*n9*(1+n14)-25*n7*n9+n7*n12^2+3*n7*n12*n13+2*n7*n12*(1+n14)-19*n7*n12+2*n7*
      n13^2+2*n7*n13*(1+n14)-24*n7*n13-5*n7*(1+n14)+55*n7+4*n8*n9*n12+4*n8*n9*(1+n14)-14*n8*n9+2*n8*n12^2+2*n8*n12*n13+2*n8*n12*(1+n14)
      -16*n8*n12+2*n8*n13*(1+n14)-8*n8*n13-8*n8*(1+n14)+30*n8+2*n9^2*n12+2*n9^2*(1+n14)-6*n9^2+3*n9*n12*n13-16*n9*n12+3*n9*n13*(1+n14)-
      10*n9*n13-16*n9*(1+n14)+53*n9+n12^2*n13+n12^2*(1+n14)-7*n12^2+n12*n13^2+2*n12*n13*(1+n14)-16*n12*n13+n12*(1+n14)^2-12*n12*(1+n14)+54*
      n12+n13^2*(1+n14)-4*n13^2+n13*(1+n14)^2-15*n13*(1+n14)+44*n13-2*(1+n14)^2+34*(1+n14)-99,4*ep^2+4*ep*n1+4*ep*n4+2*ep*n7+2*ep*n9+4*ep*
      n12+2*ep*n13+2*ep*(1+n14)-12*ep+4*n1*n4+2*n1*n7+2*n1*n12+2*n1*(1+n14)-6*n1+2*n4*n9+2*n4*n12+2*n4*n13-6*n4+n7*n9+n7*n12+n7
      *n13-3*n7+n9*n12+n9*(1+n14)-3*n9+n12^2+n12*n13+n12*(1+n14)-6*n12+n13*(1+n14)-3*n13-3*(1+n14)+9)
      +`Z'(n1,n2,1,n4,n5,2,n7,n8,-1+n9,0,0,n12,n13,(1+n14))*rat(2*ep+2*n2-2*n4+2*n8+n9-n12+n13-(1+n14)-4,2*ep+2*n1+n9+n12+n13-3
      )
      +`Z'(n1,n2,2,-1+n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-2*ep-n5-n7-2*n8-n13+5,2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,2,n4,-1+n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-2*ep-n5-n7-2*n8-n13+5,2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,2,n4,n5,1,-1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(2*ep+n5+n7+2*n8+n13-5,2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(-2+n1,n2,1,n4,n5,1,n7,n8,1+n9,0,0,n12,n13,(1+n14))*rat(2*ep*n9+2*n2*n9-2*n4*n9+2*n8*n9+n9^2-n9*n12+n9*n13-n9*(1+n14)
      -4*n9,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(-2+n1,n2,1,n4,n5,2,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(2*ep+2*n2-2*n4+2*n8+n9-n12+n13-(1+n14)-4,2*ep+2*n1+n9+n12+n13-3
      )
      +`Z'(-1+n1,-1+n2,1,n4,n5,2,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(2*ep+2*n2-2*n4+2*n8+n9-n12+n13-(1+n14)-4,2*ep+2*n1+n9+n12+
      n13-3)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,-1+n8,1+n9,0,0,n12,n13,(1+n14))*rat(n9,1)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,-1,0,n12,1+n13,(1+n14))*rat(6*ep*n13+2*n1*n13+4*n2*n13-4*n4*n13+4*n8*n13+3*n9*n13-
      n12*n13+3*n13^2-2*n13*(1+n14)-11*n13,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,0,-1,n12,1+n13,(1+n14))*rat(-4*ep*n13-2*n1*n13-2*n2*n13+2*n4*n13-2*n8*n13-2*n9*n13-2
      *n13^2+n13*(1+n14)+7*n13,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,0,0,1+n12,n13,-1+(1+n14))*rat(4*ep*n12+2*n1*n12+2*n2*n12-2*n4*n12+2*n8*n12+2*n9*n12+
      2*n12*n13-n12*(1+n14)-7*n12,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-8*ep^2-6*ep*n1-4*ep*n2+4*ep*n4-4*ep*n8-10*ep*n9-4*ep*n12-8
      *ep*n13+2*ep*(1+n14)+24*ep-2*n1^2-3*n1*n9-3*n1*n12-3*n1*n13+9*n1-4*n2*n9-2*n2*n12-2*n2*n13+4*n2+4*n4*n9+2*n4*n12+2*n4
      *n13-4*n4-4*n8*n9-2*n8*n12-2*n8*n13+4*n8-3*n9^2-n9*n12-5*n9*n13+2*n9*(1+n14)+16*n9-2*n12*n13+n12*(1+n14)+8*n12-2*n13^2+
      n13*(1+n14)+12*n13-2*(1+n14)-17,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(-1+n1,n2,1,n4,n5,2,n7,n8,-1+n9,0,0,n12,n13,(1+n14))*rat(-2*ep-2*n2+2*n4-2*n8-n9+n12-n13+(1+n14)+4,2*ep+2*n1+n9+n12+
      n13-3)
      +`Z'(n1,-2+n2,1,n4,n5,1,n7,n8,1+n9,0,0,n12,n13,(1+n14))*rat(-2*ep*n9-2*n2*n9+2*n4*n9-2*n8*n9-n9^2+n9*n12-n9*n13+n9*
      (1+n14)+4*n9,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,-1+n2,1,-1+n4,n5,1,1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(n7,1)
      +`Z'(n1,-1+n2,1,n4,-1+n5,1,1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-n7,1)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,-1+n8,1+n9,0,0,n12,n13,(1+n14))*rat(-n9,1)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,-1,0,n12,1+n13,(1+n14))*rat(-6*ep*n13-4*n1*n13-2*n2*n13+2*n4*n13-2*n8*n13-3*n9*n13-
      n12*n13-3*n13^2+n13*(1+n14)+10*n13,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(8*ep+n1+4*n2-2*n4+n7+4*n8+3*n9+3*n13-14,1)
      +`Z'(n1,n2,1,-2+n4,n5,1,1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(2*ep*n7+n5*n7+n7^2+2*n7*n8+n7*n13-5*n7,2*ep+2*n4+n7+n12+
      (1+n14)-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,1+n7,-1+n8,n9,0,0,n12,n13,(1+n14))*rat(-n7,1)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,-1,0,n12,n13,1+(1+n14))*rat(2*ep*(1+n14)+n5*(1+n14)+n7*(1+n14)+2*n8*(1+n14)+n13*(1+n14)-5*(1+n14),2*ep+2*n4+
      n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,-1,n12,n13,1+(1+n14))*rat(-2*ep*(1+n14)+2*n4*(1+n14)-2*n5*(1+n14)-n7*(1+n14)-4*n8*(1+n14)+n12*(1+n14)-2*
      n13*(1+n14)+(1+n14)^2+7*(1+n14),2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,-1+n13,(1+n14))*rat(2*ep*n12+n5*n12+n7*n12+2*n8*n12+n12*n13-5*n12,2*ep+2*
      n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(6*ep*n4-2*ep*n5-2*ep*n7-4*ep*n8+2*ep*n12-2*ep*n13+2*ep*(1+n14)+
      2*ep+2*n4^2+3*n4*n7+3*n4*n12+3*n4*(1+n14)-9*n4-2*n5*n7-n5*n12-n5*(1+n14)+2*n5-n7^2-4*n7*n8+n7*n12-2*n7*n13+n7*(1+n14)+6*n7-2*
      n8*n12-2*n8*(1+n14)+4*n8+n12^2-n12*n13+2*n12*(1+n14)-n12-n13*(1+n14)+2*n13+(1+n14)^2-(1+n14)-1,2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,-2+n5,1,1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-2*ep*n7-n5*n7-n7^2-2*n7*n8-n7*n13+5*n7,2*ep+2*n4+n7+n12+
      (1+n14)-3)
      +`Z'(n1,n2,1,n4,-1+n5,1,1+n7,-1+n8,n9,0,0,n12,n13,(1+n14))*rat(n7,1)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,-1,n12,n13,1+(1+n14))*rat(2*ep*(1+n14)+n5*(1+n14)+n7*(1+n14)+2*n8*(1+n14)+n13*(1+n14)-5*(1+n14),2*ep+2*n4+
      n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(6*ep+3*n5+3*n7+6*n8+3*n13-15,1)
      +`Z'(n1,n2,1,n4,n5,1,-1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-2*ep*n4+2*ep-n4*n5-n4*n7-2*n4*n8-n4*n13+5*n4+n5+n7+2*n8+
      n13-5,2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,-1+n8,n9,0,0,n12,n13,(1+n14))*rat(-4*ep-2*n2-n7-2*n8-n9-n13-(1+n14)+8,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,-1+n9,-1,0,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,-1+n9,0,0,n12,n13,(1+n14))*rat(-2*ep*n1+2*ep-2*n1*n2+2*n1*n4-2*n1*n8-n1*n9+n1*n12-n1*n13+
      n1*(1+n14)+4*n1+2*n2-2*n4+2*n8+n9-n12+n13-(1+n14)-4,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(2*n1*n13-2*n2*n13+2*n4*n13-2*n8*n13+2*n12*n13+n13*(1+n14)+n13,2*
      ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,0,0,n12,n13,1+(1+n14))*rat(-2*ep*(1+n14)-n5*(1+n14)-n7*(1+n14)-2*n8*(1+n14)-n13*(1+n14)+5*(1+n14),2*ep+2*n4+n7+
      n12+(1+n14)-3)
      +`Z'(n1,n2,2,-2+n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(2*ep+n5+n7+2*n8+n13-5,2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,2,-1+n4,-1+n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(2*ep+n5+n7+2*n8+n13-5,2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,2,-1+n4,n5,1,-1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-2*ep-n5-n7-2*n8-n13+5,2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(-2*ep*n12-2*n2*n12+2*n4*n12-2*n8*n12-n9*n12+n12^2-n12*n13
      +n12*(1+n14)+4*n12,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(4*ep*n13+2*n1*n13+2*n2*n13-2*n4*n13+2*n8*n13+2*n9*n13+2*
      n13^2-n13*(1+n14)-7*n13,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(6*ep*n12+4*n1*n12+2*n2*n12-2*n4*n12+2*n8*n12+3*n9*n12+
      n12^2+3*n12*n13-n12*(1+n14)-10*n12,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(2*n1*n13-2*n2*n13+2*n4*n13-2*n8*n13+2*n12*n13+n13*(1+n14)+n13
      ,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,0,0,n12,n13,1+(1+n14))*rat(2*ep*(1+n14)+4*n4*(1+n14)-n5*(1+n14)+n7*(1+n14)-2*n8*(1+n14)+2*n12*(1+n14)-n13*
      (1+n14)+2*(1+n14)^2-(1+n14),2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(2*n4*n12-n5*n12-2*n8*n12+n12^2-n12*n13+n12*(1+n14)+2*n12,2*ep
      +2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,0,n12,n13,1+(1+n14))*rat(2*ep*(1+n14)+n5*(1+n14)+n7*(1+n14)+2*n8*(1+n14)+n13*(1+n14)-5*(1+n14),2*ep+2*n4+
      n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(2*ep*n12+n5*n12+n7*n12+2*n8*n12+n12*n13-5*n12,2*ep+2*n4+
      n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(2*ep*n13+4*n1*n13-2*n2*n13+2*n4*n13-2*n8*n13+n9*n13+3*n12
      *n13+n13^2+n13*(1+n14)-2*n13,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,0,n12,n13,1+(1+n14))*rat(-2*ep*(1+n14)-n5*(1+n14)-n7*(1+n14)-2*n8*(1+n14)-n13*(1+n14)+5*(1+n14),2*ep+2*n4+
      n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,-1+n8,n9,0,0,n12,1+n13,(1+n14))*rat(-2*n1*n13+2*n2*n13-2*n4*n13+2*n8*n13-2*n12*n13-n13*(1+n14)-
      n13,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,-1+n8,n9,0,0,n12,n13,1+(1+n14))*rat(2*ep*(1+n14)+n5*(1+n14)+n7*(1+n14)+2*n8*(1+n14)+n13*(1+n14)-5*(1+n14),2*ep+2*n4+
      n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,-1+n9,0,0,1+n12,n13,(1+n14))*rat(-n12,1)
      +`Z'(-1+n1,-1+n2,1,n4,n5,1,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(4*ep*n13+2*n1*n13+2*n2*n13-2*n4*n13+2*n8*n13+2*n9*n13+
      2*n13^2-n13*(1+n14)-7*n13,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(-1+n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(-12*ep^2*n12-8*ep*n1*n12-4*ep*n2*n12-4*ep*n4*n12-2*ep*
      n5*n12-6*ep*n7*n12-8*ep*n8*n12-6*ep*n9*n12-6*ep*n12^2-8*ep*n12*n13-2*ep*n12*(1+n14)+42*ep*n12-4*n1*n4*n12-2*n1*n5*n12
      -4*n1*n7*n12-4*n1*n8*n12-2*n1*n12^2-2*n1*n12*n13-2*n1*n12*(1+n14)+16*n1*n12-4*n2*n4*n12-2*n2*n7*n12-2*n2*n12^2-2*n2*
      n12*(1+n14)+6*n2*n12+4*n4^2*n12+2*n4*n7*n12-4*n4*n8*n12-4*n4*n9*n12+2*n4*n12^2-4*n4*n12*n13+4*n4*n12*(1+n14)+8*n4*n12-n5*
      n9*n12-n5*n12^2-n5*n12*n13+3*n5*n12-2*n7*n8*n12-3*n7*n9*n12-n7*n12^2-3*n7*n12*n13+n7*n12*(1+n14)+10*n7*n12-2*n8*n9*
      n12-4*n8*n12^2-2*n8*n12*n13-2*n8*n12*(1+n14)+12*n8*n12-2*n9*n12^2-n9*n12*n13-2*n9*n12*(1+n14)+11*n9*n12-3*n12^2*n13+n12^2
      *(1+n14)+12*n12^2-n12*n13^2-2*n12*n13*(1+n14)+14*n12*n13+n12*(1+n14)^2+4*n12*(1+n14)-36*n12,4*ep^2+4*ep*n1+4*ep*n4+2*ep*n7+2*ep*
      n9+4*ep*n12+2*ep*n13+2*ep*(1+n14)-12*ep+4*n1*n4+2*n1*n7+2*n1*n12+2*n1*(1+n14)-6*n1+2*n4*n9+2*n4*n12+2*n4*n13-6*n4+n7*n9+
      n7*n12+n7*n13-3*n7+n9*n12+n9*(1+n14)-3*n9+n12^2+n12*n13+n12*(1+n14)-6*n12+n13*(1+n14)-3*n13-3*(1+n14)+9)
      +`Z'(-1+n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(2*ep*n12+n5*n12+n7*n12+2*n8*n12+n12*n13-5*n12,2*ep+2*
      n4+n7+n12+(1+n14)-3)
      +`Z'(-1+n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(2*ep*n13+2*n2*n13-2*n4*n13+2*n8*n13+n9*n13-n12*n13+
      n13^2-n13*(1+n14)-4*n13,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,-1+n8,n9,0,0,n12,1+n13,(1+n14))*rat(-4*ep*n13-2*n1*n13-2*n2*n13+2*n4*n13-2*n8*n13-2*n9*n13
      -2*n13^2+n13*(1+n14)+7*n13,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,-1+n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(2*ep*n12+2*n2*n12-2*n4*n12+2*n8*n12+n9*n12-n12^2+n12*
      n13-n12*(1+n14)-4*n12,2*ep+2*n1+n9+n12+n13-3)
      +`Z'(n1,-1+n2,1,-1+n4,n5,1,n7,n8,n9,0,0,n12,n13,1+(1+n14))*rat(-2*ep*(1+n14)-4*n4*(1+n14)+n5*(1+n14)-n7*(1+n14)+2*n8*(1+n14)-2*n12*(1+n14)+
      n13*(1+n14)-2*(1+n14)^2+(1+n14),2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Z'(n1,-1+n2,1,n4,-1+n5,1,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(-2*n13,1)
      +`Z'(n1,n2,1,-1+n4,-1+n5,1,n7,n8,n9,0,0,n12,n13,1+(1+n14))*rat(2*ep*(1+n14)+n5*(1+n14)+n7*(1+n14)+2*n8*(1+n14)+n13*(1+n14)-5*(1+n14),2*ep+2*
      n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,-1+n8,n9,0,0,n12,n13,1+(1+n14))*rat(-2*ep*(1+n14)-n5*(1+n14)-n7*(1+n14)-2*n8*(1+n14)-n13*(1+n14)+5*(1+n14),2*ep+2*
      n4+n7+n12+(1+n14)-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,-1+n9,0,0,1+n12,n13,(1+n14))*rat(n12,1)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,-1+n9,0,0,n12,1+n13,(1+n14))*rat(2*n13,1)
      +`Y'(n1,n2,0,1+n4,n5,1,n7,n8,n9,0,0,n12,n13,-1+(1+n14))*rat(n4,1)
      +`Y'(n1,n2,0,n4,1+n5,1,n7,n8,n9,0,-1,n12,n13,(1+n14))*rat(-n5,1)
      +`Y'(n1,n2,0,1+n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(n4,1)
      +`Y'(n1,n2,0,n4,1+n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(n5,1)
      +`Y'(-1+n1,n2,1,n4,n5,0,n7,n8,1+n9,0,0,n12,n13,(1+n14))*rat(-4*ep*n9-2*n1*n9-2*n2*n9+2*n4*n9-2*n8*n9-2*n9^2-2*n9*n13+
      n9*(1+n14)+7*n9,2*ep+2*n1+n9+n12+n13-3)
      +`Y'(1+n1,-1+n2,1,n4,n5,0,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(4*ep*n1+2*n1^2+2*n1*n2-2*n1*n4+2*n1*n8+2*n1*n9+2*n1*n13-
      n1*(1+n14)-7*n1,2*ep+2*n1+n9+n12+n13-3)
      +`Y'(n1,-1+n2,0,n4,1+n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(n5,1)
      +`Y'(n1,-1+n2,1,n4,n5,0,n7,n8,1+n9,0,0,n12,n13,(1+n14))*rat(4*ep*n9+2*n1*n9+2*n2*n9-2*n4*n9+2*n8*n9+2*n9^2+2*n9*n13-
      n9*(1+n14)-7*n9,2*ep+2*n1+n9+n12+n13-3)
      +`Y'(n1,n2,-1,1+n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-n4,1)
      +`Y'(n1,n2,0,-1+n4,n5,1,1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(2*n4*n7-n5*n7-2*n7*n8+n7*n12-n7*n13+n7*(1+n14)+2*n7,2*ep+2*n4
      +n7+n12+(1+n14)-3)
      +`Y'(n1,n2,0,1+n4,-1+n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(2*ep*n4+n4*n5+n4*n7+2*n4*n8+n4*n13-5*n4,2*ep+2*n4+n7+n12+
      (1+n14)-3)
      +`Y'(n1,n2,0,n4,-1+n5,1,1+n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-2*n4*n7+n5*n7+2*n7*n8-n7*n12+n7*n13-n7*(1+n14)-2*n7,2*ep+2*
      n4+n7+n12+(1+n14)-3)
      +`Y'(n1,n2,0,n4,1+n5,1,n7,-1+n8,n9,0,0,n12,n13,(1+n14))*rat(-n5,1)
      +`Y'(n1,n2,0,n4,n5,1,n7,n8,n9,0,-1,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Y'(n1,n2,0,n4,n5,1,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-8*ep^2-10*ep*n4-2*ep*n5-8*ep*n7-4*ep*n8-6*ep*n12-2*ep*n13-6*
      ep*(1+n14)+28*ep-3*n4*n5-5*n4*n7-6*n4*n8-2*n4*n12-3*n4*n13-2*n4*(1+n14)+19*n4-n5*n7-n5*n12-n5*(1+n14)+4*n5-2*n7^2-2*n7*n8-3*
      n7*n12-n7*n13-3*n7*(1+n14)+14*n7-2*n8*n12-2*n8*(1+n14)+8*n8-n12^2-n12*n13-2*n12*(1+n14)+10*n12-n13*(1+n14)+4*n13-(1+n14)^2+10*(1+n14)-26,
      2*ep+2*n4+n7+n12+(1+n14)-3)
      +`Y'(n1,n2,1,n4,n5,0,n7,n8,n9,-1,0,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(n1,n2,1,n4,n5,0,n7,n8,n9,0,-1,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(n1,n2,1,n4,n5,0,n7,n8,n9,0,0,1+n12,n13,-1+(1+n14))*rat(-n12,1)
      +`Y'(n1,n2,1,n4,n5,0,n7,n8,n9,0,0,n12,n13,(1+n14))*rat(-4*ep*n1+4*ep-2*n1^2-2*n1*n2+2*n1*n4-2*n1*n8-2*n1*n9-2*n1*n13+
      n1*(1+n14)+9*n1+2*n2-2*n4+2*n8+2*n9+2*n13-(1+n14)-7,2*ep+2*n1+n9+n12+n13-3)
      +`Y'(n1,n2,0,n4,n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(-n12,1)
      +`Y'(n1,n2,0,n4,n5,1,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(n1,n2,1,n4,n5,0,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(-1+n1,n2,0,n4,n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(2*ep*n12+2*n2*n12-2*n4*n12+2*n8*n12+n9*n12-n12^2+n12*n13-
      n12*(1+n14)-4*n12,2*ep+2*n1+n9+n12+n13-3)
      +`Y'(n1,-1+n2,0,n4,n5,1,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(-4*ep*n12-2*n1*n12-2*n2*n12+2*n4*n12-2*n8*n12-2*n9*n12-2*
      n12*n13+n12*(1+n14)+7*n12,2*ep+2*n1+n9+n12+n13-3)
      +`Y'(n1,-1+n2,0,n4,n5,1,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(n1,-1+n2,1,n4,n5,0,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(n1,n2,0,n4,n5,0,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(n12,1)
      +`Y'(n1,n2,0,n4,n5,0,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(n1,n2,0,n4,n5,1,n7,n8,-1+n9,0,0,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(n1,n2,1,-1+n4,n5,0,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(-2*n4*n12+n5*n12+2*n8*n12-n12^2+n12*n13-n12*(1+n14)-2*n12,2*
      ep+2*n4+n7+n12+(1+n14)-3)
      +`Y'(n1,n2,1,n4,-1+n5,0,n7,n8,n9,0,0,1+n12,n13,(1+n14))*rat(-2*ep*n12-n5*n12-n7*n12-2*n8*n12-n12*n13+5*n12,2*ep+2*n4+
      n7+n12+(1+n14)-3)
      +`Y'(n1,n2,1,n4,-1+n5,0,n7,n8,n9,0,0,n12,1+n13,(1+n14))*rat(-2*n13,1)
      +`Y'(n1,n2,1,n4,n5,0,n7,-1+n8,n9,0,0,n12,1+n13,(1+n14))*rat(n13,1)
      );
*--#] n14 : 
#endif
#if ( `il' >= 9 )
*--#[ n13 :
id,ifmatch->lala1,
	`Z'(n1?`POS1',n2?`POS2',1,n4?`POS4',n5?`POS5',1,n7?`POS7',n8?`POS8'
		,n9?`POS9',0,0,n12?`NEG12',n13?neg_,0)
			 = -rat(1,2*ep+2*n2+n8+n9-4)*(
      +`Z'(-1+n1,1+n2,1,1+n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(4*n2*n4,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,1+n4,n5,1,-1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*n2*n4,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,n4,n5,1,n7,n8,n9,0,0,1+n12,-1+(1+n13),0)*rat(2*n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-4*ep*n2-4*n2*n4-4*n2*n7-2*n2*n12+8*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,2,-1+n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(2*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,2,n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(2*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,2,n4,n5,1,-1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,1+n4,-1+n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(-2*n2*n4,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,1+n4,n5,1,-1+n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(n2*n4,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,n4,n5,1,n7,n8,-1+n9,0,0,1+n12,-1+(1+n13),0)*rat(-n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,n4,n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(2*ep*n2+2*n2*n4+2*n2*n7+n2*n12-4*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,2,-1+n4,n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(-n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,2,n4,-1+n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(-n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,2,n4,n5,1,-1+n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,1+n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*n4*n9+2*n4,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,1+n4,n5,1,-1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(n4*n9-n4,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,0,0,1+n12,-1+(1+n13),0)*rat(-n9*n12+n12,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(4*ep^2+4*ep*n2+4*ep*n4+2*ep*n7+2*ep*n8+4*ep*n9+2*ep*n12-16*ep+4*
      n2*n4+2*n2*n7+2*n2*n12-6*n2+2*n4*n8+4*n4*n9-10*n4+n7*n8+3*n7*n9-6*n7+n8*n12-3*n8+2*n9*n12-7*n9-5*n12+16,2*ep+2*n4
      +n7+n12-3)
      +`Z'(n1,n2,2,-1+n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-n9+1,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,2,n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-n9+1,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,2,n4,n5,1,-1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(n9-1,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,-2+n4,n5,1,1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*n2*n7,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,-1+(1+n13),0)*rat(-2*n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,-1+n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(4*ep*n2+4*n2*n7+2*n2*n12-4*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,n4,-2+n5,1,1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(2*n2*n7,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-4*n2,1)
      +`Z'(-1+n1,1+n2,1,n4,n5,1,-1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(2*n2*n4-2*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,n4,n5,1,n7,-1+n8,n9,0,0,n12,(1+n13),0)*rat(n2,1)
      +`Z'(-1+n1,1+n2,2,-2+n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,2,-1+n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,2,-1+n4,n5,1,-1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(2*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,-1+n8,1+n9,0,0,n12,(1+n13),0)*rat(n9,1)
      +`Z'(-1+n1,n2,1,n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(2*ep+2*n2+2*n8-4,1)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,-1+n8,1+n9,0,0,n12,(1+n13),0)*rat(-n9,1)
      +`Z'(n1,-1+n2,1,n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-4*ep-2*n2-3*n8-2*n9+9,1)
      +`Z'(n1,1+n2,1,-2+n4,n5,1,1+n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(n2*n7,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,-1+n4,n5,1,n7,n8,-1+n9,0,0,1+n12,-1+(1+n13),0)*rat(n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,-1+n4,n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(-2*ep*n2-2*n2*n7-n2*n12+2*n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,n4,-2+n5,1,1+n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(-n2*n7,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,n4,-1+n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(2*n2,1)
      +`Z'(n1,1+n2,1,n4,n5,1,-1+n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(-n2*n4+n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,2,-2+n4,n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,2,-1+n4,-1+n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,2,-1+n4,n5,1,-1+n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(-n2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,-2+n4,n5,1,1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(n7*n9-n7,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,-1+(1+n13),0)*rat(n9*n12-n12,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*ep*n9+2*ep-2*n7*n9+2*n7-n9*n12+2*n9+n12-2,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,n4,-2+n5,1,1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(-n7*n9+n7,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(2*n9-2,1)
      +`Z'(n1,n2,1,n4,n5,1,-1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(-n4*n9+n4+n9-1,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,n4,n5,1,n7,-1+n8,n9,0,0,n12,(1+n13),0)*rat(-n2+1,1)
      +`Z'(n1,n2,1,n4,n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(2*ep+n2+n8+n9-4,1)
      +`Z'(n1,n2,2,-2+n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(n9-1,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,2,-1+n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(n9-1,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,2,-1+n4,n5,1,-1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(-n9+1,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(2*n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,n4,-1+n5,1,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(-2*n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,-1+n4,n5,1,n7,n8,-1+n9,0,0,1+n12,(1+n13),0)*rat(-n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,1+n2,1,n4,-1+n5,1,n7,n8,-1+n9,0,0,1+n12,(1+n13),0)*rat(n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(-n9*n12+n12,2*ep+2*n4+n7+n12-3)
      +`Z'(n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(n9*n12-n12,2*ep+2*n4+n7+n12-3)
      +`Z'(-2+n1,1+n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(2*n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(-2+n1,1+n2,1,n4,-1+n5,1,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(-2*n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,-1+n4,n5,1,n7,n8,-1+n9,0,0,1+n12,(1+n13),0)*rat(-n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,1+n2,1,n4,-1+n5,1,n7,n8,-1+n9,0,0,1+n12,(1+n13),0)*rat(n2*n12,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,n2,1,-1+n4,n5,1,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(-n9*n12+n12,2*ep+2*n4+n7+n12-3)
      +`Z'(-1+n1,n2,1,n4,-1+n5,1,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(n9*n12-n12,2*ep+2*n4+n7+n12-3)
      +`Y'(-1+n1,1+n2,0,-1+n4,n5,1,1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(2*n2*n7,2*ep+2*n4+n7+n12-3)
      +`Y'(-1+n1,1+n2,0,1+n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*n2*n4,2*ep+2*n4+n7+n12-3)
      +`Y'(-1+n1,1+n2,0,n4,-1+n5,1,1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*n2*n7,2*ep+2*n4+n7+n12-3)
      +`Y'(-1+n1,1+n2,0,n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(2*n2*n4-2*n2,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,1+n2,0,-1+n4,n5,1,1+n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(-n2*n7,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,1+n2,0,1+n4,-1+n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(n2*n4,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,1+n2,0,n4,-1+n5,1,1+n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(n2*n7,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,1+n2,0,n4,n5,1,n7,n8,-1+n9,0,0,n12,(1+n13),0)*rat(-n2*n4+n2,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,n2,0,-1+n4,n5,1,1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(-n7*n9+n7,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,n2,0,1+n4,-1+n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(n4*n9-n4,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,n2,0,n4,-1+n5,1,1+n7,n8,n9,0,0,n12,(1+n13),0)*rat(n7*n9-n7,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,n2,0,n4,n5,1,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-n4*n9+n4+n9-1,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,n2,1,n4,n5,0,n7,n8,n9,0,0,n12,(1+n13),0)*rat(-2*ep-2*n2-n8-n9+4,1)
      +`Y'(-1+n1,1+n2,1,-1+n4,n5,0,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(-2*n2*n12,2*ep+2*n4+n7+n12-3)
      +`Y'(-1+n1,1+n2,1,n4,-1+n5,0,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(2*n2*n12,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,1+n2,1,-1+n4,n5,0,n7,n8,-1+n9,0,0,1+n12,(1+n13),0)*rat(n2*n12,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,1+n2,1,n4,-1+n5,0,n7,n8,-1+n9,0,0,1+n12,(1+n13),0)*rat(-n2*n12,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,n2,1,-1+n4,n5,0,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(n9*n12-n12,2*ep+2*n4+n7+n12-3)
      +`Y'(n1,n2,1,n4,-1+n5,0,n7,n8,n9,0,0,1+n12,(1+n13),0)*rat(-n9*n12+n12,2*ep+2*n4+n7+n12-3)
      );
*--#] n13 : 
#endif
#if ( `il' >= 10 )
*--#[ n8 :
id,ifmatch->lala1,
	`Z'(n1?`POS1',n2?`POS2',1,n4?`POS4',n5?`POS5',1,n7?`POS7',n8?{>1}
		,n9?`POS9',0,0,n12?`NEG12',0,0)
			 = -rat(1,2*ep+2*(n8-1)-2)*(
      +`Z'(n1,n2,1,-1+n4,1+n5,1,1+n7,(n8-1),n9,0,0,n12,0,0)*rat(2*n5*n7,(n8-1))
      +`Z'(n1,n2,1,n4,1+n5,1,n7,(n8-1),n9,0,0,n12,0,0)*rat(-2*ep*n5-2*n5^2-2*n5*n7+2*n5,(n8-1))
      +`Z'(n1,n2,1,n4,n5,1,1+n7,(n8-1),n9,0,0,n12,0,0)*rat(-2*ep*n7-2*n5*n7-2*n7^2+2*n7,(n8-1))
      );
*--#] n8 : 
#endif
#if ( `il' >= 11 )
*--#[ n2 :
id,ifmatch->lala1,
	`Z'(n1?`POS1',n2?{>1},1,n4?`POS4',n5?`POS5',1,n7?`POS7',1
		,n9?`POS9',0,0,n12?`NEG12',0,0)
			 = -rat(n7,2*ep*(n2-1)+2*(n2-1)^2+2*(n2-1)*n9-2*(n2-1))*(
      +`Z'(-1+n1,1+(n2-1),1,n4,n5,1,n7,1,1+n9,0,0,n12,0,0)*rat(-2*(n2-1)*n9,n7)
      +`Z'(n1,(n2-1),1,-1+n4,1+n5,1,1+n7,1,n9,0,0,n12,0,0)*rat(2*n5,1)
      +`Z'(n1,(n2-1),1,n4,1+n5,1,n7,1,n9,0,0,n12,0,0)*rat(-2*ep*n5-2*n5^2-2*n5*n7+2*n5,n7)
      +`Z'(n1,(n2-1),1,n4,n5,1,1+n7,1,n9,0,0,n12,0,0)*rat(-2*ep-2*n5-2*n7+2,1)
      +`Z'(n1,(n2-1),1,n4,n5,1,n7,1,1+n9,0,0,n12,0,0)*rat(2*ep*n9+2*(n2-1)*n9+2*n9^2-2*n9,n7)
      );
*--#] n2 : 
#endif
#if ( `il' >= 12 )
*--#[ n1 :
id,ifmatch->lala1,
	`Z'(n1?{>1},1,1,n4?`POS4',n5?`POS5',1,n7?`POS7',1
		,n9?`POS9',0,0,n12?`NEG12',0,0)
			 = -rat(1,-2*ep*(n1-1)-2*(n1-1)^2-2*(n1-1)*n9+2*(n1-1))*(
      +`Z'((n1-1),1,1,n4,n5,1,n7,1,1+n9,0,0,n12,0,0)*rat(-2*ep*n9-2*(n1-1)*n9-2*n9^2+2*n9,1)
      +`Z'(-1+(n1-1),1,1,n4,n5,1,n7,1,1+n9,0,0,n12,0,0)*rat(2*ep*n9+2*n9*n12,1)
      +`Z'(-1+(n1-1),1,1,n4,n5,2,n7,1,n9,0,0,n12,0,0)*rat(2*ep+2*n12,1)
      +`Z'((n1-1),1,1,n4,n5,1,n7,1,n9,0,0,n12,0,0)*rat(4*ep^2+4*ep*(n1-1)+2*ep*n9+6*ep*n12-6*ep+4*(n1-1)*n12+2*n9*n12+2*n12^2-6*n12,1)
      +`Z'((n1-1),1,1,n4,n5,2,n7,1,n9,0,0,1+n12,0,0)*rat(n12,1)
      +`Z'(-1+(n1-1),1,1,n4,n5,1,n7,1,1+n9,0,0,1+n12,0,0)*rat(n9*n12,1)
      +`Z'(-1+(n1-1),1,1,n4,n5,2,n7,1,n9,0,0,1+n12,0,0)*rat(n12,1)
      +`Z'((n1-1),1,1,-1+n4,n5,2,n7,1,n9,0,0,1+n12,0,0)*rat(2*n12,1)
      +`Z'((n1-1),1,1,n4,n5,1,n7,1,n9,0,0,1+n12,0,0)*rat(2*ep*n12+2*(n1-1)*n12+n9*n12+n12^2-2*n12,1)
      +`Z'((n1-1),1,1,-1+n4,n5,1,n7,1,n9,0,0,2+n12,0,0)*rat(2*n12^2+2*n12,1)
      +`Y'(1+(n1-1),0,1,n4,n5,1,n7,1,1+n9,0,0,n12,0,0)*rat(2*(n1-1)*n9,1)
      +`Y'((n1-1),0,1,n4,n5,1,n7,1,1+n9,0,0,n12,0,0)*rat(-2*ep*n9-2*n9*n12,1)
      +`Y'((n1-1),0,1,n4,n5,1,n7,1,1+n9,0,0,1+n12,0,0)*rat(-n9*n12,1)
      +`Y'((n1-1),1,0,n4,n5,2,n7,1,n9,0,0,1+n12,0,0)*rat(-2*n12,1)
      );
*--#] n1 : 
#endif
#if ( `il' >= 13 )
*--#[ n4 :
id,ifmatch->lala1,
	`Z'(1,1,1,n4?{>1},n5?`POS5',1,n7?`POS7',1
		,n9?`POS9',0,0,n12?`NEG12',0,0)
			 = -rat(1,2*ep*(n4-1)+2*(n4-1)^2+2*(n4-1)*n7-2*(n4-1))*(
      +`Z'(1,1,1,1+(n4-1),-1+n5,1,1+n7,1,n9,0,0,n12,0,0)*rat(-2*(n4-1)*n7,1)
      +`Z'(1,1,1,(n4-1),n5,1,1+n7,1,n9,0,0,n12,0,0)*rat(2*ep*n7+2*(n4-1)*n7+2*n7^2-2*n7,1)
      +`Z'(1,1,1,-1+(n4-1),n5,1,1+n7,1,n9,0,0,n12,0,0)*rat(-2*ep*n7-2*n7*n12,1)
      +`Z'(1,1,1,(n4-1),-1+n5,1,1+n7,1,n9,0,0,n12,0,0)*rat(2*ep*n7+2*n7*n12,1)
      +`Z'(1,1,1,(n4-1),n5,1,n7,1,n9,0,0,n12,0,0)*rat(-4*ep^2-4*ep*(n4-1)-2*ep*n7-6*ep*n12+6*ep-4*(n4-1)*n12-2*n7*n12-2*n12^2+6*n12,1)
      +`Z'(1,1,2,-1+(n4-1),n5,1,n7,1,n9,0,0,n12,0,0)*rat(-2*ep-2*n12,1)
      +`Z'(1,1,2,(n4-1),n5,1,n7,1,n9,0,0,1+n12,0,0)*rat(-n12,1)
      +`Z'(1,1,1,-1+(n4-1),n5,1,1+n7,1,n9,0,0,1+n12,0,0)*rat(-n7*n12,1)
      +`Z'(1,1,1,(n4-1),-1+n5,1,1+n7,1,n9,0,0,1+n12,0,0)*rat(n7*n12,1)
      +`Z'(1,1,1,(n4-1),n5,1,n7,1,n9,0,0,1+n12,0,0)*rat(-2*ep*n12-2*(n4-1)*n12-n7*n12-n12^2+2*n12,1)
      +`Z'(1,1,2,-1+(n4-1),n5,1,n7,1,n9,0,0,1+n12,0,0)*rat(-n12,1)
      +`Y'(0,1,2,(n4-1),n5,1,n7,1,n9,0,0,1+n12,0,0)*rat(-2*n12,1)
      +`Y'(1,1,2,(n4-1),n5,0,n7,1,n9,0,0,1+n12,0,0)*rat(2*n12,1)
      +`Y'(0,1,1,(n4-1),n5,1,n7,1,n9,0,0,2+n12,0,0)*rat(-2*n12^2-2*n12,1)
      );
*--#] n4 : 
#endif
#if ( `il' >= 14 )
*--#[ n5 :
id,ifmatch->lala1,
	`Z'(1,1,1,1,n5?{>1},1,n7?`POS7',1
		,n9?`POS9',0,0,n12?`NEG12',0,0)
		 = -rat(ep+(n5-1)+n7-2,-2*ep^2*(n5-1)-2*ep*(n5-1)^2-3*ep*(n5-1)*n7-ep*(n5-1)*n12+3*ep*(n5-1)-(n5-1)^2*n7-(n5-1)^2*n12+
      (n5-1)^2-(n5-1)*n7^2-(n5-1)*n7*n12+2*(n5-1)*n7+(n5-1)*n12-(n5-1))*(
      +`Z'(1,1,1,1,(n5-1),1,1+n7,1,n9,0,0,n12,0,0)*rat(ep*(n5-1)*n7-ep*n7*n12-ep*n7+(n5-1)^2*n7+(n5-1)*n7^2-(n5-1)*n7*n12-2*(n5-1)*n7-n7^2*n12
      -n7^2+n7*n12+n7,ep+(n5-1)+n7-2)
      +`Z'(1,1,1,1,-1+(n5-1),1,1+n7,1,n9,0,0,n12,0,0)*rat(-6*ep^2*n7-5*ep*(n5-1)*n7-6*ep*n7^2-2*ep*n7*n9+13*ep*n7-(n5-1)^2*n7-3*(n5-1)*
      n7^2-(n5-1)*n7*n9+6*(n5-1)*n7-2*n7^3-n7^2*n9+8*n7^2+2*n7*n9-8*n7,ep+(n5-1)+n7-2)
      +`Z'(1,1,1,1,-1+(n5-1),1,n7,2,n9,0,0,n12,0,0)*rat(ep*n7-ep,ep+(n5-1)+n7-2)
      +`Z'(1,1,1,1,1+(n5-1),1,-1+n7,1,n9,0,0,n12,0,0)*rat(2*ep*(n5-1)+(n5-1)*n7+(n5-1)*n12-2*(n5-1),1)
      +`Z'(1,1,1,1,(n5-1),1,-1+n7,2,n9,0,0,n12,0,0)*rat(-2*ep^2-ep*n7-ep*n12+2*ep,ep+(n5-1)+n7-2)
      +`Z'(1,1,1,1,(n5-1),1,n7,1,n9,0,0,n12,0,0)*rat(12*ep^3+10*ep^2*(n5-1)+18*ep^2*n7+4*ep^2*n9+6*ep^2*n12-32*ep^2+2*ep*(n5-1)^2+
      10*ep*(n5-1)*n7+2*ep*(n5-1)*n9+5*ep*(n5-1)*n12-16*ep*(n5-1)+10*ep*n7^2+4*ep*n7*n9+6*ep*n7*n12-35*ep*n7+2*ep*n9*n12-6*ep*n9-13*ep*
      n12+29*ep+(n5-1)^2*n12+2*(n5-1)*n7^2+(n5-1)*n7*n9+3*(n5-1)*n7*n12-6*(n5-1)*n7+(n5-1)*n9*n12-(n5-1)*n9-6*(n5-1)*n12+4*(n5-1)+2*n7^3+n7^2*n9+2*n7^2*n12
      -10*n7^2+n7*n9*n12-3*n7*n9-8*n7*n12+16*n7-2*n9*n12+2*n9+8*n12-8,ep+(n5-1)+n7-2)
      +`Z'(1,1,1,1,(n5-1),1,n7,2,n9,0,0,1+n12,0,0)*rat(-ep*n12,ep+(n5-1)+n7-2)
      +`Z'(1,1,1,1,(n5-1),1,n7,1,n9,0,0,1+n12,0,0)*rat(-2*ep*n9*n12+2*ep*n12-(n5-1)*n9*n12+(n5-1)*n12-n7*n9*n12+n7*n12+2*n9*n12-2*
      n12,ep+(n5-1)+n7-2)
      +`Z'(1,1,1,1,(n5-1),1,n7,2,-1+n9,0,0,1+n12,0,0)*rat(-ep*n12,ep+(n5-1)+n7-2)
      +`Z'(1,2,1,1,(n5-1),1,n7,1,-1+n9,0,0,1+n12,0,0)*rat(-2*ep*n12-(n5-1)*n12-n7*n12+2*n12,ep+(n5-1)+n7-2)
      +`Y'(1,1,1,0,1+(n5-1),1,1+n7,1,n9,0,0,n12,0,0)*rat(ep*(n5-1)*n7-(n5-1)^2*n7+(n5-1)*n7^2+(n5-1)*n7*n12,ep+(n5-1)+n7-2)
      +`Y'(0,1,1,0,(n5-1),1,1+n7,1,1+n9,0,0,n12,0,0)*rat(-2*ep*n7*n9-(n5-1)*n7*n9-n7^2*n9+2*n7*n9,ep+(n5-1)+n7-2)
      +`Y'(0,1,1,1,-1+(n5-1),1,1+n7,1,1+n9,0,0,n12,0,0)*rat(2*ep*n7*n9+(n5-1)*n7*n9+n7^2*n9-2*n7*n9,ep+(n5-1)+n7-2)
      +`Y'(0,1,1,1,(n5-1),1,n7,1,1+n9,0,0,n12,0,0)*rat(-4*ep^2*n9-2*ep*(n5-1)*n9-4*ep*n7*n9-2*ep*n9*n12+6*ep*n9-(n5-1)*n7*n9-(n5-1)*n9*
      n12+(n5-1)*n9-n7^2*n9-n7*n9*n12+3*n7*n9+2*n9*n12-2*n9,ep+(n5-1)+n7-2)
      +`Y'(0,1,2,0,(n5-1),1,n7,1,1+n9,0,0,n12,0,0)*rat(-2*ep*n9-(n5-1)*n9-n7*n9+2*n9,ep+(n5-1)+n7-2)
      +`Y'(1,0,1,0,(n5-1),1,1+n7,1,1+n9,0,0,n12,0,0)*rat(2*ep*n7*n9+(n5-1)*n7*n9+n7^2*n9-2*n7*n9,ep+(n5-1)+n7-2)
      +`Y'(1,0,1,0,(n5-1),1,1+n7,2,n9,0,0,n12,0,0)*rat(ep*n7,ep+(n5-1)+n7-2)
      +`Y'(1,0,1,1,-1+(n5-1),1,1+n7,1,1+n9,0,0,n12,0,0)*rat(-2*ep*n7*n9-(n5-1)*n7*n9-n7^2*n9+2*n7*n9,ep+(n5-1)+n7-2)
      +`Y'(1,0,1,1,-1+(n5-1),1,1+n7,2,n9,0,0,n12,0,0)*rat(-ep*n7,ep+(n5-1)+n7-2)
      +`Y'(1,0,1,1,(n5-1),1,n7,1,1+n9,0,0,n12,0,0)*rat(4*ep^2*n9+2*ep*(n5-1)*n9+4*ep*n7*n9+2*ep*n9*n12-6*ep*n9+(n5-1)*n7*n9+(n5-1)*n9*
      n12-(n5-1)*n9+n7^2*n9+n7*n9*n12-3*n7*n9-2*n9*n12+2*n9,ep+(n5-1)+n7-2)
      +`Y'(1,0,1,1,(n5-1),1,n7,2,n9,0,0,n12,0,0)*rat(2*ep^2+ep*n7+ep*n12-ep,ep+(n5-1)+n7-2)
      +`Y'(1,0,2,0,(n5-1),1,n7,1,1+n9,0,0,n12,0,0)*rat(2*ep*n9+(n5-1)*n9+n7*n9-2*n9,ep+(n5-1)+n7-2)
      +`Y'(1,0,2,0,(n5-1),1,n7,2,n9,0,0,n12,0,0)*rat(ep,ep+(n5-1)+n7-2)
      +`Y'(1,1,1,-1,1+(n5-1),1,1+n7,1,n9,0,0,n12,0,0)*rat(-(n5-1)*n7,1)
      +`Y'(1,1,1,0,1+(n5-1),1,n7,1,n9,0,0,n12,0,0)*rat(-2*ep*(n5-1)-(n5-1)*n12+2*(n5-1),1)
      +`Y'(1,1,1,0,(n5-1),1,1+n7,1,n9,0,0,n12,0,0)*rat(6*ep^2*n7+6*ep*(n5-1)*n7+6*ep*n7^2+2*ep*n7*n9-13*ep*n7+2*(n5-1)^2*n7+4*(n5-1)*
      n7^2+(n5-1)*n7*n9-8*(n5-1)*n7+2*n7^3+n7^2*n9-8*n7^2-2*n7*n9+8*n7,ep+(n5-1)+n7-2)
      +`Y'(1,1,1,0,(n5-1),1,n7,2,n9,0,0,n12,0,0)*rat(-ep*n7+ep,ep+(n5-1)+n7-2)
      +`Y'(1,1,2,-1,1+(n5-1),1,n7,1,n9,0,0,n12,0,0)*rat(-(n5-1),1)
      +`Y'(1,1,2,0,1+(n5-1),1,-1+n7,1,n9,0,0,n12,0,0)*rat((n5-1),1)
      +`Y'(1,1,2,0,(n5-1),1,-1+n7,2,n9,0,0,n12,0,0)*rat(-ep,ep+(n5-1)+n7-2)
      +`Y'(1,1,2,0,(n5-1),1,n7,1,n9,0,0,n12,0,0)*rat(6*ep^2+5*ep*(n5-1)+6*ep*n7+2*ep*n9-13*ep+(n5-1)^2+3*(n5-1)*n7+(n5-1)*n9-6*(n5-1)+2*n7^2+n7
      *n9-8*n7-2*n9+8,ep+(n5-1)+n7-2)
      +`Y'(0,1,1,1,(n5-1),1,n7,1,1+n9,0,0,1+n12,0,0)*rat(-2*ep*n9*n12-(n5-1)*n9*n12-n7*n9*n12+2*n9*n12,ep+(n5-1)+n7-2)
      +`Y'(0,2,1,1,(n5-1),1,n7,1,n9,0,0,1+n12,0,0)*rat(2*ep*n12+(n5-1)*n12+n7*n12-2*n12,ep+(n5-1)+n7-2)
      +`Y'(1,0,1,1,(n5-1),1,n7,1,1+n9,0,0,1+n12,0,0)*rat(2*ep*n9*n12+(n5-1)*n9*n12+n7*n9*n12-2*n9*n12,ep+(n5-1)+n7-2)
      +`Y'(1,0,1,1,(n5-1),1,n7,2,n9,0,0,1+n12,0,0)*rat(ep*n12,ep+(n5-1)+n7-2)
      +`Y'(1,1,1,1,(n5-1),0,n7,2,n9,0,0,1+n12,0,0)*rat(ep*n12,ep+(n5-1)+n7-2)
      );
*--#] n5 : 
#endif
#if ( `il' >= 15 )
*--#[ n12 :
id,ifmatch->lala1,
	`Z'(1,1,1,1,1,1,n7?`POS7',1,n9?`POS9',0,0,n12?neg_,0,0)
			 = -rat(1,-4*ep-2*n7-2*(1+n12)+4)*(
      +`Z'(1,1,1,1,1,1,n7,1,n9,0,0,1+(1+n12),-1,0)*rat(-2*ep*n9*(1+n12)+2*ep*(1+n12)-n7*n9*(1+n12)+n7*(1+n12)+2*n9*(1+n12)-2*(1+n12),4*ep^2+2*ep*n7
      +2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Z'(1,1,1,1,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(8*ep^3+4*ep^2*n7+16*ep^2*n9-4*ep^2*(1+n12)-20*ep^2+12*ep*n7*n9-4*ep*n7*(1+n12)-
      12*ep*n7+4*ep*n9^2+4*ep*n9*(1+n12)-22*ep*n9-4*ep*(1+n12)^2+18*ep+2*n7^2*n9-2*n7^2+2*n7*n9^2-n7*n9*(1+n12)-10*n7*n9+n7*(1+n12)+8*
      n7+2*n9^2*(1+n12)-2*n9^2-2*n9*(1+n12)^2-4*n9*(1+n12)+8*n9+2*(1+n12)^2+2*(1+n12)-6,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*
      (1+n12)-n9-(1+n12)+1)
      +`Z'(1,1,1,2,1,1,-1+n7,1,n9,0,0,(1+n12),0,0)*rat(2*ep*n9-2*ep+n7*n9-n7-2*n9+2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7
      *n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Z'(1,1,2,1,1,1,-1+n7,1,n9,0,0,(1+n12),0,0)*rat(2*ep*n9-2*ep+n7*n9-n7-2*n9+2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7
      *n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Z'(1,2,1,1,1,1,n7,1,-1+n9,0,0,1+(1+n12),-1,0)*rat(-2*ep*(1+n12)-n7*(1+n12)+2*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9
      -n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Z'(1,2,1,1,1,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(8*ep^2+8*ep*n7+2*ep*n9+4*ep*(1+n12)-12*ep+2*n7^2+n7*n9+n7*(1+n12)-7*n7+n9*
      (1+n12)-n9-3*(1+n12)+5,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Z'(1,2,1,2,1,1,-1+n7,1,-1+n9,0,0,(1+n12),0,0)*rat(2*ep+n7-2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9
      -(1+n12)+1)
      +`Z'(1,2,2,1,1,1,-1+n7,1,-1+n9,0,0,(1+n12),0,0)*rat(2*ep+n7-2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9
      -(1+n12)+1)
      +`Z'(1,1,1,1,1,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(8*ep^2+2*ep*n7+6*ep*n9-14*ep+n7*n9-2*n7+n9^2-5*n9+6,2*ep+n9-1)
      +`Y'(1,1,1,0,1,1,1+n7,1,n9,0,0,-1+(1+n12),0,0)*rat(-n7,1)
      +`Y'(1,1,1,1,0,1,1+n7,1,n9,0,0,-1+(1+n12),0,0)*rat(n7,1)
      +`Y'(1,1,1,2,0,1,n7,1,n9,0,0,-1+(1+n12),0,0)*rat(-1,1)
      +`Y'(0,1,2,1,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-1,1)
      +`Y'(0,2,1,1,1,1,n7,1,n9,0,0,1+(1+n12),-1,0)*rat(4*ep*(1+n12)+2*n7*(1+n12)-4*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,1,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-12*ep^2-14*ep*n7-2*ep*n9-6*ep*(1+n12)+20*ep-4*n7^2-n7*n9-2*n7*(1+n12)+13*n7-n9
      *(1+n12)+n9+5*(1+n12)-9,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,2,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(8*ep+4*n7-8,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)
      +1)
      +`Y'(0,2,1,2,1,1,-1+n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep-2*n7+4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9
      -(1+n12)+1)
      +`Y'(0,2,2,0,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(4*ep+2*n7-4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)
      +1)
      +`Y'(0,2,2,1,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(4*ep+2*n7-4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)
      +1)
      +`Y'(0,2,2,1,1,1,-1+n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep-2*n7+4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9
      -(1+n12)+1)
      +`Y'(1,1,1,-1,1,1,1+n7,1,n9,0,0,-1+(1+n12),0,0)*rat(-n7,1)
      +`Y'(1,1,1,0,1,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(-n7,1)
      +`Y'(1,1,1,0,1,1,n7,1,n9,0,0,-1+(1+n12),0,0)*rat(-2*ep-n7-(1+n12)+4,1)
      +`Y'(1,1,1,1,-1,1,1+n7,1,n9,0,0,-1+(1+n12),0,0)*rat(n7,1)
      +`Y'(1,1,1,1,0,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(n7,1)
      +`Y'(1,1,1,1,0,1,n7,1,n9,0,0,-1+(1+n12),0,0)*rat(-2*ep-n7-(1+n12)+2,1)
      +`Y'(1,1,1,2,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep*n9+4*ep-2*n7*n9+2*n7+4*n9-4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+
      n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,2,-1,1,1,n7,1,n9,0,0,-1+(1+n12),0,0)*rat(-1,1)
      +`Y'(1,1,2,0,0,1,n7,1,n9,0,0,-1+(1+n12),0,0)*rat(-1,1)
      +`Y'(1,1,2,0,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep^2-2*ep*n7-4*ep*n9-2*ep*(1+n12)+6*ep-2*n7*n9+2*n7-n9*(1+n12)+3*n9+(1+n12)-3,4*
      ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,2,1,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-2*ep*n9+2*ep-n7*n9+n7+2*n9-2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*
      n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,2,1,1,0,n7,1,n9,0,0,(1+n12),0,0)*rat(1,1)
      +`Y'(1,2,1,2,0,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(-4*ep-2*n7+4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9
      -(1+n12)+1)
      +`Y'(1,2,2,0,1,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(-2*ep-n7+2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-
      (1+n12)+1)
      +`Y'(1,2,2,1,0,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(-2*ep-n7+2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-
      (1+n12)+1)
      +`Y'(0,1,1,1,0,1,n7,1,1+n9,0,0,(1+n12),0,0)*rat(n9,1)
      +`Y'(0,1,1,1,1,0,n7,1,1+n9,0,0,(1+n12),0,0)*rat(n9,1)
      +`Y'(0,1,1,1,1,1,n7,0,1+n9,0,0,(1+n12),0,0)*rat(2*ep*n9+n7*n9-2*n9,2*ep+n9-1)
      +`Y'(0,1,1,1,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(4*ep^2+2*ep*n7-2*ep*n9-2*ep*(1+n12)-n9^2-n9*(1+n12)+3*n9+(1+n12)-2,2*ep+n9-1)
      +`Y'(0,1,2,1,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-1,1)
      +`Y'(0,2,0,0,1,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(4*ep*n7+2*n7^2-4*n7,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*
      (1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,0,1,0,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep*n7-2*n7^2+4*n7,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9
      *(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,0,2,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep-2*n7+4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-
      (1+n12)+1)
      +`Y'(0,2,1,-1,1,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep*n7-2*n7^2+4*n7,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+
      n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,0,1,1,n7,1,n9,0,0,1+(1+n12),-1,0)*rat(-4*ep*(1+n12)-2*n7*(1+n12)+4*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,0,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(8*ep^2+12*ep*n7+4*ep*(1+n12)-16*ep+4*n7^2+2*n7*(1+n12)-12*n7-4*(1+n12)+8,4*ep^2+2*
      ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,1,-1,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(4*ep*n7+2*n7^2-4*n7,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9
      *(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,1,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-10*ep-4*n7-n9+9,2*ep+n9-1)
      +`Y'(0,2,1,1,1,1,n7,0,n9,0,0,(1+n12),0,0)*rat(2*ep+n7-2,2*ep+n9-1)
      +`Y'(0,2,2,-1,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep-2*n7+4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-
      (1+n12)+1)
      +`Y'(0,2,2,0,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep-2*n7+4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-
      (1+n12)+1)
      +`Y'(0,2,2,0,1,1,-1+n7,1,n9,0,0,(1+n12),0,0)*rat(4*ep+2*n7-4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-
      (1+n12)+1)
      +`Y'(1,0,1,1,0,1,n7,1,1+n9,0,0,(1+n12),0,0)*rat(-n9,1)
      +`Y'(1,0,1,1,1,0,n7,1,1+n9,0,0,(1+n12),0,0)*rat(-n9,1)
      +`Y'(1,0,1,1,1,1,n7,0,1+n9,0,0,(1+n12),0,0)*rat(-2*ep*n9-n7*n9+2*n9,2*ep+n9-1)
      +`Y'(1,0,1,1,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-12*ep^2-4*ep*n7-8*ep*n9+22*ep-2*n7*n9+4*n7-n9^2+7*n9-10,2*ep+n9-1)
      +`Y'(1,1,0,0,1,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(4*ep^2*n7+2*ep*n7^2+2*ep*n7*(1+n12)-2*ep*n7+n7*n9*(1+n12)+n7*n9-n7*(1+n12)-n7,4*
      ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,0,1,0,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep^2*n7-2*ep*n7^2-2*ep*n7*(1+n12)+2*ep*n7-n7*n9*(1+n12)-n7*n9+n7*(1+n12)+n7,4*
      ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,0,1,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(2*ep+n7+(1+n12)-2,1)
      +`Y'(1,1,0,2,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(2*ep*n9-2*ep+n7*n9-n7-2*n9+2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9
      -n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,1,-1,1,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(2*ep*n7*n9-2*ep*n7+n7^2*n9-n7^2-2*n7*n9+2*n7,4*ep^2+2*ep*n7+2*ep*n9+
      2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,1,0,0,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(-n7,1)
      +`Y'(1,1,1,0,1,1,n7,1,n9,0,0,1+(1+n12),-1,0)*rat(2*ep*n9*(1+n12)-2*ep*(1+n12)+n7*n9*(1+n12)-n7*(1+n12)-2*n9*(1+n12)+2*(1+n12),4*ep^2+2*ep*n7+
      2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,1,0,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep^2*n9+4*ep^2-6*ep*n7*n9+6*ep*n7-2*ep*n9*(1+n12)+8*ep*n9+2*ep*(1+n12)-8*ep-
      2*n7^2*n9+2*n7^2-n7*n9*(1+n12)+6*n7*n9+n7*(1+n12)-6*n7+2*n9*(1+n12)-4*n9-2*(1+n12)+4,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,1,1,-1,1,1+n7,1,n9,0,0,(1+n12),0,0)*rat(4*ep^2*n7+2*ep*n7^2+2*ep*n7*(1+n12)-2*ep*n7+n7*n9*(1+n12)+n7*n9-n7*(1+n12)-n7,4*
      ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,1,1,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep^2-2*ep*n7+4*ep*n9-2*ep*(1+n12)-2*ep+n7*n9-n7+n9^2-n9*(1+n12)-5*n9+(1+n12)+4,2
      *ep+n9-1)
      +`Y'(1,1,1,1,1,0,n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep-n7-n9+(1+n12)+2,1)
      +`Y'(1,1,2,-1,1,1,n7,1,n9,0,0,(1+n12),0,0)*rat(2*ep*n9-2*ep+n7*n9-n7-2*n9+2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*
      n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,2,0,0,1,n7,1,n9,0,0,(1+n12),0,0)*rat(-4*ep^2-2*ep*n7-2*ep*(1+n12)+2*ep-n9*(1+n12)-n9+(1+n12)+1,4*ep^2+2*ep*n7+2*ep*n9+2*
      ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,2,0,1,1,-1+n7,1,n9,0,0,(1+n12),0,0)*rat(-2*ep*n9+2*ep-n7*n9+n7+2*n9-2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+
      n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,2,1,0,0,n7,1,n9,0,0,(1+n12),0,0)*rat(1,1)
      +`Y'(1,2,0,0,1,1,1+n7,1,-1+n9,0,0,(1+n12),0,0)*rat(-2*ep*n7-n7^2+2*n7,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+
      n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,0,1,0,1,1+n7,1,-1+n9,0,0,(1+n12),0,0)*rat(2*ep*n7+n7^2-2*n7,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9
      *(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,0,2,0,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(2*ep+n7-2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-
      (1+n12)+1)
      +`Y'(1,2,1,-1,1,1,1+n7,1,-1+n9,0,0,(1+n12),0,0)*rat(2*ep*n7+n7^2-2*n7,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+
      n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,1,0,1,1,n7,1,-1+n9,0,0,1+(1+n12),-1,0)*rat(2*ep*(1+n12)+n7*(1+n12)-2*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,1,0,1,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(-4*ep^2-6*ep*n7-2*ep*(1+n12)+8*ep-2*n7^2-n7*(1+n12)+6*n7+2*(1+n12)-4,4*ep^2+2*ep
      *n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,1,1,-1,1,1+n7,1,-1+n9,0,0,(1+n12),0,0)*rat(-2*ep*n7-n7^2+2*n7,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+
      n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,1,1,0,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(6*ep+2*n7+n9-5,2*ep+n9-1)
      +`Y'(1,2,2,-1,1,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(2*ep+n7-2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-
      (1+n12)+1)
      +`Y'(1,2,2,0,0,1,n7,1,-1+n9,0,0,(1+n12),0,0)*rat(2*ep+n7-2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-
      (1+n12)+1)
      +`Y'(1,2,2,0,1,1,-1+n7,1,-1+n9,0,0,(1+n12),0,0)*rat(-2*ep-n7+2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-
      n9-(1+n12)+1)
      +`Y'(0,1,1,1,1,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(-2*(1+n12),1)
      +`Y'(0,2,1,0,1,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(4*ep*(1+n12)+2*n7*(1+n12)-4*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7
      +n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,1,0,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(-4*ep*(1+n12)-2*n7*(1+n12)+4*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,1,0,1,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(-2*ep*n9*(1+n12)+2*ep*(1+n12)-n7*n9*(1+n12)+n7*(1+n12)+2*n9*(1+n12)-2*(1+n12),4*ep^2+2*ep*n7+
      2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,1,1,0,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(2*ep*n9*(1+n12)-2*ep*(1+n12)+n7*n9*(1+n12)-n7*(1+n12)-2*n9*(1+n12)+2*(1+n12),4*ep^2+2*ep*n7+2
      *ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,1,0,1,1,n7,1,-1+n9,0,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)-n7*(1+n12)+2*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,1,1,0,1,n7,1,-1+n9,0,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)+n7*(1+n12)-2*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(-1,2,1,0,1,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(4*ep*(1+n12)+2*n7*(1+n12)-4*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(-1,2,1,1,0,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(-4*ep*(1+n12)-2*n7*(1+n12)+4*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,1,0,1,1,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(2*(1+n12),1)
      +`Y'(0,1,1,0,1,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(-2*ep*n9*(1+n12)+2*ep*(1+n12)-n7*n9*(1+n12)+n7*(1+n12)+2*n9*(1+n12)-2*(1+n12),4*ep^2+2*ep*n7+
      2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,1,1,1,0,1,n7,1,n9,0,0,1+(1+n12),0,0)*rat(-8*ep^2*(1+n12)-4*ep*n7*(1+n12)-2*ep*n9*(1+n12)-4*ep*(1+n12)^2+6*ep*(1+n12)-n7*n9*(1+n12)+n7*
      (1+n12)-2*n9*(1+n12)^2+2*(1+n12)^2,4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,0,1,0,n7,1,n9,0,0,1+(1+n12),0,0)*rat(-4*ep*(1+n12)-2*n7*(1+n12)+4*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,0,1,1,n7,1,-1+n9,0,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)-n7*(1+n12)+2*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,1,0,0,n7,1,n9,0,0,1+(1+n12),0,0)*rat(4*ep*(1+n12)+2*n7*(1+n12)-4*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7
      +n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(0,2,1,1,0,1,n7,1,-1+n9,0,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)+n7*(1+n12)-2*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,1,0,1,0,n7,1,n9,0,0,1+(1+n12),0,0)*rat(2*ep*n9*(1+n12)-2*ep*(1+n12)+n7*n9*(1+n12)-n7*(1+n12)-2*n9*(1+n12)+2*(1+n12),4*ep^2+2*ep*n7+2
      *ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,1,1,1,0,0,n7,1,n9,0,0,1+(1+n12),0,0)*rat(-2*ep*n9*(1+n12)+2*ep*(1+n12)-n7*n9*(1+n12)+n7*(1+n12)+2*n9*(1+n12)-2*(1+n12),4*ep^2+2*ep*n7+
      2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,1,0,1,0,n7,1,-1+n9,0,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)+n7*(1+n12)-2*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      +`Y'(1,2,1,1,0,0,n7,1,-1+n9,0,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)-n7*(1+n12)+2*(1+n12),4*ep^2+2*ep*n7+2*ep*n9+2*ep*(1+n12)-4*ep+n7*n9-
      n7+n9*(1+n12)-n9-(1+n12)+1)
      );
*--#] n12 : 
#endif
#if ( `il' >= 16 )
*--#[ n9 :
id,ifmatch->lala1,
	`Z'(1,1,1,1,1,1,n7?`POS7',1,n9?{>1},0,0,0,0,0)
			 = -rat(ep+(n9-1)-1,2*ep^2*(n9-1)+3*ep*(n9-1)^2-ep*(n9-1)+(n9-1)^3-(n9-1)^2)*(
      +`Z'(1,1,1,1,1,1,1+n7,1,-1+(n9-1),0,0,0,0,0)*rat(2*ep^3*n7+2*ep^2*n7^2+ep^2*n7*(n9-1)-2*ep^2*n7+ep*n7^2*(n9-1)-2*ep*n7^2,ep^2
      +2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Z'(1,1,1,1,1,1,n7,1,(n9-1),0,0,0,0,0)*rat(8*ep^2+2*ep*n7+8*ep*(n9-1)-10*ep+n7*(n9-1)-n7+2*(n9-1)^2-6*(n9-1)+4,1)
      +`Z'(1,1,1,1,2,1,n7,1,-1+(n9-1),0,0,0,0,0)*rat(2*ep^3+2*ep^2*n7+ep^2*(n9-1)-2*ep^2+ep*n7*(n9-1)-2*ep*n7,ep^2+2*ep*(n9-1)-2*ep+
      (n9-1)^2-2*(n9-1)+1)
      +`Z'(1,2,1,1,1,1,n7,1,-1+(n9-1),0,0,0,0,0)*rat(2*ep+(n9-1)-2,1)
      +`Y'(0,2,1,1,1,1,n7,1,1+(n9-1),0,0,0,0,0)*rat(-ep*(n9-1),ep+(n9-1)-1)
      +`Y'(2,0,1,0,2,1,1+n7,1,(n9-1),0,0,0,0,0)*rat(-n7*(n9-1)+n7,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(2,0,1,1,1,1,1+n7,1,(n9-1),0,0,0,0,0)*rat(ep*n7*(n9-1)-ep*n7+n7^2*(n9-1)-n7^2,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(2,0,1,1,1,1,n7,1,1+(n9-1),0,0,0,0,0)*rat(-(n9-1),1)
      +`Y'(2,0,1,1,2,1,n7,1,(n9-1),0,0,0,0,0)*rat(ep*(n9-1)-ep+n7*(n9-1)-n7,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(-1,2,1,1,1,1,n7,1,1+(n9-1),0,0,0,0,0)*rat(-(n9-1),1)
      +`Y'(-1,2,1,1,1,2,n7,1,(n9-1),0,0,0,0,0)*rat(-1,1)
      +`Y'(0,1,1,0,1,1,1+n7,1,1+(n9-1),0,0,0,0,0)*rat(-n7*(n9-1),1)
      +`Y'(0,1,1,0,1,2,1+n7,1,(n9-1),0,0,0,0,0)*rat(-n7,1)
      +`Y'(0,1,1,0,2,1,1+n7,1,(n9-1),0,0,0,0,0)*rat(-ep*n7*(n9-1)+ep*n7,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(0,1,1,0,2,2,1+n7,1,-1+(n9-1),0,0,0,0,0)*rat(-ep*n7,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(0,1,1,1,0,1,1+n7,1,1+(n9-1),0,0,0,0,0)*rat(n7*(n9-1),1)
      +`Y'(0,1,1,1,0,2,1+n7,1,(n9-1),0,0,0,0,0)*rat(n7,1)
      +`Y'(0,1,1,1,1,1,1+n7,1,(n9-1),0,0,0,0,0)*rat(ep^2*n7*(n9-1)-ep^2*n7+ep*n7^2*(n9-1)-ep*n7^2,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(0,1,1,1,1,1,n7,1,1+(n9-1),0,0,0,0,0)*rat(4*ep*(n9-1)+n7*(n9-1)+2*(n9-1)^2-2*(n9-1),1)
      +`Y'(0,1,1,1,1,2,1+n7,1,-1+(n9-1),0,0,0,0,0)*rat(ep^2*n7+ep*n7^2,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(0,1,1,1,1,2,n7,1,(n9-1),0,0,0,0,0)*rat(4*ep+n7+2*(n9-1)-3,1)
      +`Y'(0,1,1,1,2,1,n7,1,(n9-1),0,0,0,0,0)*rat(ep^2*(n9-1)-ep^2+ep*n7*(n9-1)-ep*n7,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(0,1,1,1,2,2,n7,1,-1+(n9-1),0,0,0,0,0)*rat(ep^2+ep*n7,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(0,2,1,1,1,1,n7,1,(n9-1),0,0,0,0,0)*rat(-2*ep+2,1)
      +`Y'(0,2,1,1,1,2,n7,1,-1+(n9-1),0,0,0,0,0)*rat(1,1)
      +`Y'(1,0,1,0,1,1,1+n7,1,1+(n9-1),0,0,0,0,0)*rat(n7*(n9-1),1)
      +`Y'(1,0,1,0,2,1,1+n7,1,(n9-1),0,0,0,0,0)*rat(ep*n7*(n9-1)-ep*n7,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(1,0,1,1,0,1,1+n7,1,1+(n9-1),0,0,0,0,0)*rat(-n7*(n9-1),1)
      +`Y'(1,0,1,1,1,1,1+n7,1,(n9-1),0,0,0,0,0)*rat(-ep^2*n7*(n9-1)+ep^2*n7-ep*n7^2*(n9-1)+ep*n7^2,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(1,0,1,1,1,1,n7,1,1+(n9-1),0,0,0,0,0)*rat(-4*ep*(n9-1)-n7*(n9-1)-2*(n9-1)^2+3*(n9-1),1)
      +`Y'(1,0,1,1,2,1,n7,1,(n9-1),0,0,0,0,0)*rat(-ep^2*(n9-1)+ep^2-ep*n7*(n9-1)+ep*n7,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(1,1,1,0,1,1,1+n7,1,(n9-1),0,0,0,0,0)*rat(-2*ep*n7-n7*(n9-1)+n7,1)
      +`Y'(1,1,1,0,2,1,1+n7,1,-1+(n9-1),0,0,0,0,0)*rat(-2*ep^2*n7-ep*n7*(n9-1)+2*ep*n7,ep^2+2*ep*(n9-1)-2*ep+(n9-1)^2-2*(n9-1)+1)
      +`Y'(1,1,1,1,0,1,1+n7,1,(n9-1),0,0,0,0,0)*rat(2*ep*n7+n7*(n9-1)-n7,1)
      );
*--#] n9 : 
#endif
#if ( `il' >= 17 )
*--#[ n7 :
id,ifmatch->lala1,
	`Z'(1,1,1,1,1,1,n7?{>1},1,1,0,0,0,0,0)
			 = -rat(1,2*ep*(n7-1)+2*(n7-1)^2)*(
      +`Z'(1,1,1,1,1,1,-1+(n7-1),2,1,0,0,0,0,0)*rat(-4*ep^2-2*ep*(n7-1)+4*ep,2*ep+(n7-1)-1)
      +`Z'(1,1,1,1,1,1,(n7-1),1,1,0,0,0,0,0)*rat(16*ep^3+28*ep^2*(n7-1)-28*ep^2+18*ep*(n7-1)^2-38*ep*(n7-1)+20*ep+4*(n7-1)^3-14*(n7-1)^2+16*(n7-1)-
      6,2*ep+(n7-1)-1)
      +`Z'(1,1,1,1,2,1,-1+(n7-1),1,1,0,0,0,0,0)*rat(4*ep^2+6*ep*(n7-1)-8*ep+2*(n7-1)^2-6*(n7-1)+4,2*ep+(n7-1)-1)
      +`Y'(1,1,1,0,2,1,1+(n7-1),1,1,0,0,0,0,0)*rat(-2*ep*(n7-1),2*ep+(n7-1)-1)
      +`Y'(2,0,1,0,2,1,1+(n7-1),1,1,0,0,0,0,0)*rat((n7-1),ep)
      +`Y'(2,0,1,1,1,1,1+(n7-1),1,1,0,0,0,0,0)*rat(-ep*(n7-1)-(n7-1)^2,ep)
      +`Y'(2,0,1,1,2,1,(n7-1),1,1,0,0,0,0,0)*rat(-ep-(n7-1),ep)
      +`Y'(0,1,1,0,1,1,1+(n7-1),1,2,0,0,0,0,0)*rat(-(n7-1),1)
      +`Y'(0,1,1,0,1,2,1+(n7-1),1,1,0,0,0,0,0)*rat((n7-1),1)
      +`Y'(0,1,1,0,2,2,1+(n7-1),1,0,0,0,0,0,0)*rat(-(n7-1),ep)
      +`Y'(0,1,1,1,0,1,1+(n7-1),1,2,0,0,0,0,0)*rat((n7-1),1)
      +`Y'(0,1,1,1,0,1,(n7-1),2,2,0,0,0,0,0)*rat(-1,1)
      +`Y'(0,1,1,1,0,2,1+(n7-1),1,1,0,0,0,0,0)*rat(-(n7-1),1)
      +`Y'(0,1,1,1,0,2,(n7-1),2,1,0,0,0,0,0)*rat(-1,1)
      +`Y'(0,1,1,1,1,1,(n7-1),1,2,0,0,0,0,0)*rat(-6*ep-3*(n7-1)+3,1)
      +`Y'(0,1,1,1,1,2,1+(n7-1),1,0,0,0,0,0,0)*rat(ep*(n7-1)+(n7-1)^2,ep)
      +`Y'(0,1,1,1,1,2,(n7-1),1,1,0,0,0,0,0)*rat(-2*ep-(n7-1)+1,1)
      +`Y'(0,1,1,1,1,2,(n7-1),2,0,0,0,0,0,0)*rat(1,1)
      +`Y'(0,1,1,1,2,2,(n7-1),1,0,0,0,0,0,0)*rat(ep+(n7-1),ep)
      +`Y'(0,1,2,0,1,1,(n7-1),1,2,0,0,0,0,0)*rat(-2,1)
      +`Y'(1,0,1,0,1,1,1+(n7-1),1,2,0,0,0,0,0)*rat((n7-1),1)
      +`Y'(1,0,1,0,1,1,1+(n7-1),2,1,0,0,0,0,0)*rat(2*ep*(n7-1),2*ep+(n7-1)-1)
      +`Y'(1,0,1,1,0,1,1+(n7-1),1,2,0,0,0,0,0)*rat(-(n7-1),1)
      +`Y'(1,0,1,1,0,1,1+(n7-1),2,1,0,0,0,0,0)*rat(-2*ep*(n7-1),2*ep+(n7-1)-1)
      +`Y'(1,0,1,1,0,1,(n7-1),2,2,0,0,0,0,0)*rat(1,1)
      +`Y'(1,0,1,1,1,1,(n7-1),1,2,0,0,0,0,0)*rat(6*ep+3*(n7-1)-3,1)
      +`Y'(1,0,1,1,1,1,(n7-1),2,1,0,0,0,0,0)*rat(2*ep,1)
      +`Y'(1,0,2,0,1,1,(n7-1),1,2,0,0,0,0,0)*rat(2,1)
      +`Y'(1,0,2,0,1,1,(n7-1),2,1,0,0,0,0,0)*rat(2*ep,2*ep+(n7-1)-1)
      +`Y'(1,1,1,-1,2,1,1+(n7-1),1,1,0,0,0,0,0)*rat(-2*ep*(n7-1)-2*(n7-1)^2+2*(n7-1),2*ep+(n7-1)-1)
      +`Y'(1,1,1,0,1,1,1+(n7-1),1,1,0,0,0,0,0)*rat(16*ep^2*(n7-1)+14*ep*(n7-1)^2-12*ep*(n7-1)+4*(n7-1)^3-6*(n7-1)^2+2*(n7-1),2*ep+(n7-1)-1)
      +`Y'(1,1,1,0,1,1,(n7-1),2,1,0,0,0,0,0)*rat(-2*ep*(n7-1)+2*ep,2*ep+(n7-1)-1)
      +`Y'(1,1,1,0,2,1,1+(n7-1),1,0,0,0,0,0,0)*rat(-2*ep*(n7-1)+(n7-1),ep)
      +`Y'(1,1,1,0,2,1,(n7-1),1,1,0,0,0,0,0)*rat(-4*ep^2-4*ep*(n7-1)+8*ep+4*(n7-1)-4,2*ep+(n7-1)-1)
      +`Y'(1,1,1,1,0,1,1+(n7-1),1,1,0,0,0,0,0)*rat(-16*ep^2*(n7-1)-14*ep*(n7-1)^2+14*ep*(n7-1)-4*(n7-1)^3+8*(n7-1)^2-4*(n7-1),2*ep+(n7-1)-1)
      +`Y'(1,1,1,1,0,1,(n7-1),2,1,0,0,0,0,0)*rat(-4*ep^2,2*ep+(n7-1)-1)
      +`Y'(1,1,1,1,1,1,1+(n7-1),1,0,0,0,0,0,0)*rat(2*ep^2*(n7-1)+2*ep*(n7-1)^2-ep*(n7-1)-(n7-1)^2,ep)
      +`Y'(1,1,1,1,1,1,(n7-1),2,0,0,0,0,0,0)*rat(2*ep-1,1)
      +`Y'(1,1,1,1,2,1,(n7-1),1,0,0,0,0,0,0)*rat(2*ep^2+2*ep*(n7-1)-ep-(n7-1),ep)
      +`Y'(1,1,2,-1,2,1,(n7-1),1,1,0,0,0,0,0)*rat(-2*ep-2*(n7-1)+2,2*ep+(n7-1)-1)
      +`Y'(1,1,2,0,1,1,-1+(n7-1),2,1,0,0,0,0,0)*rat(-2*ep,2*ep+(n7-1)-1)
      +`Y'(1,1,2,0,1,1,(n7-1),1,1,0,0,0,0,0)*rat(12*ep^2+12*ep*(n7-1)-12*ep+4*(n7-1)^2-8*(n7-1)+4,2*ep+(n7-1)-1)
      +`Y'(1,1,2,0,2,1,-1+(n7-1),1,1,0,0,0,0,0)*rat(2*ep+2*(n7-1)-2,2*ep+(n7-1)-1)
      );
*--#] n7 : 
#endif
goto lala2;
Label lala1;
$ired = 1;
Label lala2;
id	`Z'(?a,n?neg0_,?b,n10?,...,n14?) = `Y'(?a,n,?b,n10,...,n14);
ModuleOption,maximum,$ired;
.sort:Main lala loop, pass `$pass';
#if ( `$ired' > 0 )
	#redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
id	`Z'(1,1,1,1,1,1,1,1,1,0,0,0,0,0) = Master(lala);
*
*--#] Loop : 
#break
*--#] lala : 
*--#[ nono :
#case nono
*
#do il = 1,15
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering nono loop `il', pass `$pass' at `time_' sec."
#endif
#if ( `il' >= 1 )
*--#[ stat :
*
*  Prework to make the reduction of the nonoJ topology faster.
*
* B18
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n7)*(`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,1+n13,n14)*rat(-n7+1,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,
   n14)*rat(-n5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,
   n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,1+n5,n6,-1+n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n5,1)+`Z'(
   n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,n11,n12,
   1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-2+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+
   n7,n8,n9,n10,-1+n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,-1+n12,1+n13,n14)*rat(-
   n11,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,
   n12,n13,n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,-1+n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-
   1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(2*ep+n5+
   n7+2*n8+n10+n11-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n7,1));
* B2
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n7)*(`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,
   n9,1+n10,n11,n12,n13,n14)*rat(-n7+1,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,1+n9,1+n10,n11,n12,n13,n14)*rat(n9,1)+`Z'(n1,-1+
   n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(2*n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,
   n14)*rat(-1+n7,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,1+n9,1+n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+
   n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n5,1)+
   `Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,1+n10,n11,
   n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,1+n5,n6,-1+n7,-1+n8,n9,1+n10,n11,n12,n13,n14)*rat(n5,1)+`Z'(n1,n2,n3,n4,n5,-1+
   n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(
   n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-2+n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,1+n10,
   n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,-1+n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,
   n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,-1+n12,n13,1+n14)*
   rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,-1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+
   n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(4*ep+2*n2+n5+n7+2*n8+
   n9+2*n10+n11+n12+n13-7,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-
   1+n8,n9,1+n10,n11,n12,n13,n14)*rat(-1+n7,1));
* B19
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n7)*(`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(n1,1+n2,-1+n3,n4,n5,n6,-1+n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(-n2,1)+`Z'(n1,1+n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(n2,1)+`Z'(n1,
   n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+
   n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,
   -1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,
   1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,n10,1+
   n11,n12,n13,1+n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,
   n5,-1+n6,-1+n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,1+n11,n12,n13,1+n14)*
   rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-2+n7,n8,n9,1+
   n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-2+n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,
   n5,n6,-2+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,1+n11,n12,n13,1+n14)*
   rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+
   n9,n10,n11,n12,n13,2+n14)*rat(2*n14+2,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,-1+n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(
   n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,-1+n12,n13,1+n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,
   n13,1+n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,-1+n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-
   1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n11,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,-1+n12,1+n13,1+n14)*rat(n13,1)+
   `Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+
   n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(2*ep+n2+2*n9+n11+n12+n14-3,1));
* B1
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?pos_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n7)*(`Z'(n1,1+n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n2,1)+`Z'(n1,1+n2,n3,n4,n5,n6,-1+n7,n8,-
   1+n9,n10,1+n11,n12,n13,n14)*rat(-n2,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,
   n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,1+
   n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,1+n8,n9,n10,1+n11,n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+
   n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)
   +`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,n9,n10,1+
   n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,-1+
   n5,n6,-1+n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,1+n14)*
   rat(n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,
   n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,
   n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,2+n11,n12,n13,n14)*rat(2*n11+2,1)+`Z'(n1,n2,n3,n4,n5,n6,-2+n7,n8,n9,1+n10,1+n11,n12,
   n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-2+n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-2+n7,
   n8,n9,n10,1+n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,1+n10,1+n11,n12,n13,n14)*rat(n10,1)+`Z'(
   n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,1+n11,
   n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-
   1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,-1+n12,1+n13,n14)*rat(-n13,1)+
   `Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,1+
   n13,-1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,
   n8,n9,n10,1+n11,n12,n13,n14)*rat(-4*ep-n2-2*n5-2*n7-n8-2*n9-n10-n11-n12-n14+9,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,
   n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1));
* B6
id,ifmatch->nono1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-2+n2,n3,n4,n5,n6,n7,n8,
   n9,1+n10,n11,1+n12,n13,n14)*rat(-n10,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-n10,1)+`Z'(n1,
   -1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,2+n12,
   n13,n14)*rat(1+n12,1)+`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,-1+n2,1+n3,n4,n5,n6,
   n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(
   n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-n10,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+
   n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,2+n12,n13,n14)*rat(-n12-1,1)+`Z'(n1,-1+n2,n3,n4,-
   1+n5,n6,n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-n8,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(
   n10,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,
   n11,1+n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-n10,1)+`Z'(n1,-1+n2,
   n3,n4,n5,n6,n7,n8,-1+n9,1+n10,n11,1+n12,n13,n14)*rat(n10,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,1+n13,
   n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,
   n10,-1+n11,1+n12,1+n13,n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,-1+n14)*rat(n13,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(2*ep+n2+n3+2*n4+n8+n13-5,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,
   n11,2+n12,n13,n14)*rat(-n12-1,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,-1+n3,n4,
   n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n2,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n2+1,
   1));
* B11
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n9)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-
   1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n9+1,1)+`Z'(n1,
   n2,1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(n3,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,
   n14)*rat(-n3,1)+`Z'(n1,n2,n3,-1+n4,1+n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(n5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,n14)*rat(-1+n9,1)+`Z'(n1,n2,n3,n4,1+n5,n6,n7,-1+n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n5,1)+`Z'(n1,
   n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,1+n12,
   n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,-1+n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,
   n8,-1+n9,n10,n11,1+n12,-1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(2*ep+n3+2*
   n4+n5+n9+n14-5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n9+1,1));
* B5
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n6)*(`Z'(1+n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n1,1)+`Z'(1+n1,n2,n3,n4,n5,-1+n6,-1+n7,
   n8,n9,n10,1+n11,n12,n13,n14)*rat(n1,1)+`Z'(1+n1,n2,n3,n4,n5,-1+n6,n7,-1+n8,n9,n10,1+n11,n12,n13,n14)*rat(-n1,1)+`Z'(1+n1,
   n2,n3,n4,n5,-1+n6,n7,n8,n9,-1+n10,1+n11,n12,n13,n14)*rat(-n1,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,1+
   n13,n14)*rat(-2*n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,-1+n7,n8,n9,n10,2+n11,n12,n13,n14)*rat(2*n11+2,1)+`Z'(n1,n2,n3,n4,n5,-1+
   n6,n7,n8,n9,-1+n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(2*ep+
   n1+2*n7+n11+n13+n14-3,1));
* B13
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n4)*(`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,
   n9,1+n10,n11,1+n12,n13,n14)*rat(2*n10,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,
   -1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,
   n14)*rat(n3,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,
   n10,n11,1+n12,n13,n14)*rat(-2*ep-2*n2-n3-n10-n12-n13+3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(
   -1+n4,1));
* B14
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n5+1)*(`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(2*n10,1)+`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,
   n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n4,1)+`Z'(n1,
   n2,n3,n4,-1+n5,n6,n7,1+n8,n9,-1+n10,1+n11,n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,
   n14)*rat(-4*ep-2*n2-2*n3-n4-n5-2*n9-n10-n12-n13-n14+9,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-
   1+n5,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+
   n11,n12,n13,n14)*rat(-n5+1,1));
* B10
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-1+n9)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-
   1+n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(
   n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,1+n11,n12,n13,n14)*rat(-n2,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,1+n11,
   n12,n13,n14)*rat(n3,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,n12,n13,n14)*rat(-n3,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,
   n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n9,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(
   n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n9+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,1+n11,1+n12,
   n13,n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,1+n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,
   n8,-1+n9,n10,1+n11,n12,n13,1+n14)*rat(-2*n14,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n9+1,1)+
   `Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,1+n10,1+n11,n12,n13,n14)*rat(-2*n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,1+
   n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,
   n5,n6,n7,n8,-1+n9,-1+n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,1+n11,n12,1+n13,n14)
   *rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,
   n10,1+n11,n12,n13,n14)*rat(4*ep+n3+2*n4+2*n5+2*n8+n9+n10+n12+n13+n14-9,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+
   n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(n14,1));
* B3
id,ifmatch->nono1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n1+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,n12,1+n13,n14)*rat(-n9,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,1+n11,n12,1+n13,n14)*rat(2*n11,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-2*n14,1)+`Z'(
   -1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,2+n13,n14)*rat(-2*n13-2,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,n11,
   n12,1+n13,n14)*rat(n9,1)+`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,n10,n11,-1+n12,1+n13,n14)*rat(n6,1)+`Z'(-1+n1,n2,n3,n4,n5,1+
   n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,n12,1+n13,n14)*rat(n11,1)+`Z'(
   -1+n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,
   n13,n14)*rat(-n11,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-2*ep-2*n2-n9-n10-n12-n13+3,1)+`Z'(-1+
   n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1));
* B8
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n5+1)*(`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,1+n12,n13,n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,
   n9,1+n10,n11,1+n12,n13,n14)*rat(2*n10,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-2*n14,1)+`Z'(
   n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n5,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,1+n12,
   n13,n14)*rat(n9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n5,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+
   n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(
   n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,-1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,1+n12,
   n13,n14)*rat(-2*ep-2*n2-n9-n10-n12-n13+3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n5+1,1)+`Z'(n1,
   n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(-1+n5,1));
* B7
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{>1},n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n9)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-
   1+n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(
   n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+
   n12,n13,1+n14)*rat(n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,1+n2,n3,n4,n5,
   n6,n7,n8,-1+n9,-1+n10,n11,n12,n13,1+n14)*rat(n2,1)+`Z'(n1,1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,-1+n13,1+n14)*rat(-
   n2,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,1+n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,
   n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,
   n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,n13,2+n14)*rat(2*n14+2,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,1+n10,n11,n12,n13,1+
   n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-
   1+n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,
   n2,n3,n4,n5,n6,n7,n8,-1+n9,-1+n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,-1+
   n13,1+n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,-1+n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,
   n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(2*ep+2*n2+2*n7+n11+n14-3,1));
* B17
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n5+1)*(`Z'(-1+n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,1+
   n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8,1)+`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(n1,
   n2,n3,1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,-1+n7,1+n8,n9,n10,n11,n12,n13,
   1+n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-n8,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,1+
   n8,n9,n10,n11,n12,-1+n13,1+n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,-1+n13,1+n14)*rat(-n10,1)+`Z'(
   n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(6*ep+2*n2+2*n3+n4+n5+2*n7+n8+2*n9+2*n10+n11+n12+n13+2*n14-
   11,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n5+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,
   n12,n13,1+n14)*rat(-1+n5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-1+n5,1));
* B16
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?pos_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n8+1)*(`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n8+1,1)+`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,-1+n8,
   n9,n10,1+n11,n12,n13,n14)*rat(n2,1)+`Z'(n1,1+n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,1+n11,n12,n13,n14)*rat(-n2,1)+`Z'(n1,n2,-
   1+n3,1+n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,n13,n14)*rat(n4,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,1+n10,1+n11,n12,n13,
   n14)*rat(-n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+
   n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,-1+n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(
   n1,n2,n3,-1+n4,n5,n6,n7,-1+n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,-1+n8,n9,n10,1+n11,
   n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n8,1)+`Z'(n1,n2,n3,1+n4,n5,n6,
   n7,-1+n8,n9,n10,1+n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,-1+n8,n9,1+n10,1+n11,n12,n13,n14)*rat(n10,1)+
   `Z'(n1,n2,n3,n4,-1+n5,n6,n7,-1+n8,n9,n10,1+n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,
   n13,n14)*rat(-n8+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,-1+n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,
   n7,-1+n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n8,1)+
   `Z'(n1,n2,n3,n4,n5,n6,-1+n7,-1+n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,-1+n8,n9,n10,1+n11,
   n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n8,1)+`Z'(n1,n2,n3,n4,n5,n6,
   n7,-1+n8,-1+n9,1+n10,1+n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,1+n11,n12,1+n13,n14)*rat(-
   n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,1+n10,1+n11,-1+n12,n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,1+n10,
   1+n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,n4,n5,n6,
   n7,-1+n8,n9,n10,1+n11,-1+n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,
   1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,1+n13,-1+n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,
   n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,n13,n14)*rat(2*ep+n2+2*n3+n4+n10+n13-4,1)+`Z'(
   n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,n12,n13,
   n14)*rat(-n8+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,-1+n12,n13,n14)*rat(-1+n8,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,1+n11,n12,n13,n14)*rat(-n8+1,1));
* B15
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n3+1)*(`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,
   n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(
   n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n3+1,1)+`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,n12,
   n13,n14)*rat(-n2,1)+`Z'(n1,n2,-1+n3,1+n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,
   n7,-1+n8,n9,n10,1+n11,n12,n13,n14)*rat(n4,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+
   `Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,
   n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-2*n10,1)+`Z'(n1,n2,-1+n3,n4,
   n5,n6,n7,-1+n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,1+n11,n12,1+n13,n14)*
   rat(-n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,-
   1+n10,1+n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,-
   1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(2*ep+n4+2*n8+n10+n12+n13-4,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,
   n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n3,1)+`Z'(n1,n2,n3,n4,-1+
   n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n3+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n3,
   1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,n12,
   n13,n14)*rat(-n3+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,-1+n12,n13,n14)*rat(-1+n3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,
   n9,n10,1+n11,n12,n13,n14)*rat(-n3+1,1));
* B4
id,ifmatch->nono1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n1+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,
   n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(
   -1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9,n10,n11,
   n12,n13,1+n14)*rat(n8,1)+`Z'(-1+n1,n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,-1+
   n5,n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(
   n6,1)+`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(n6,1)+`Z'(-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9,n10,
   n11,n12,n13,n14)*rat(-n6,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,1+n11,n12,n13,1+n14)*rat(2*n11,1)+`Z'(-1+n1,n2,n3,
   n4,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*
   rat(n13,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,1+
   n10,n11,-1+n12,n13,1+n14)*rat(n10,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n11,1)+`Z'(-1+n1,n2,
   n3,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,1+n14)*rat(n13,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*
   rat(-2*ep-2*n5-n6-n8-n10-n11+4,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n1+1,1)+`Z'(n1,n2,-1+n3,
   n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+
   n1,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,
   n12,n13,1+n14)*rat(2*n1-2,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,
   n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-1+n1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-n1+1,1)+
   `Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-n1+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,n13,
   1+n14)*rat(-1+n1,1));
* B12
id,ifmatch->nono1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?pos_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n3+1)*(`Z'(-1+n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,1+
   n8,n9,n10,n11,n12,n13,1+n14)*rat(n8,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(-1+
   n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+
   n13,1+n14)*rat(-n13,1)+`Z'(n1,-1+n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,
   n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(n8,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-3*
   n10,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,
   n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(2*n3-2,1)+`Z'(n1,1+n2,-1+n3,
   n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(n2,1)+`Z'(n1,1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,-1+n13,1+n14)*
   rat(-n2,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,-1+n8,
   n9,n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,-
   1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,-1+n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*
   rat(-n4,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,
   1+n10,n11,n12,n13,1+n14)*rat(n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,-
   1+n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,1+n10,n11,n12,n13,
   1+n14)*rat(n10,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-
   1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+n8,n9,n10,n11,n12,-1+n13,1+n14)*rat(-n8,1)+`Z'(
   n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,
   1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,-1+n13,1+n14)*rat(2*n10,1)+`Z'(n1,n2,-1+n3,n4,n5,
   n6,n7,n8,n9,n10,-1+n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(2*n2-
   n8-n10,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-1+n3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,
   n12,-1+n13,1+n14)*rat(-n3+1,1));
*--#] stat : 
#endif
#if ( `il' >= 2 )
*--#[ n1 : 
*`SKIP4'
id,ifmatch->nono1,
	`Z'(n1?{>1},n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',n6?`POS6',n7?`POS7'
			,n8?`POS8',n9?`POS9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(-1,(n1-1))*(
      +`Z'(1+(n1-1),n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat((n1-1),1)
      +`Z'((n1-1),1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n2,1)
      +`Z'((n1-1),1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,n14)*rat(n2,1)
      +`Z'((n1-1),n2,n3,-1+n4,n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)
      +`Z'((n1-1),n2,n3,n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(n8,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(-n10,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,-1+n14)*rat(-n13,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(6*ep+(n1-1)+n2+2*n5+2*n6+2*n7+n8+2*n9+n10+2*n11+n12+n13+2*n14-12,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'((n1-1),n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)
      +`Z'((n1-1),n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)
      +`Z'((n1-1),n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'((n1-1),n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)
      +`Z'((n1-1),n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n12,1)
      +`Z'((n1-1),n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)
      +`Z'((n1-1),n2,n3,n4,-1+n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)
      +`Z'((n1-1),n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'((n1-1),n2,n3,n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)
      +`Z'((n1-1),n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(n13,1)
      );
*--#] n1 : 
#endif
#if ( `il' >= 3 )
*--#[ n6 : 
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',n6?{>1},n7?`POS7'
			,n8?`POS8',n9?`POS9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(1,(n6-1))*(
      +`Z'(1,-1+n2,n3,n4,n5,(n6-1),n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)
      +`Z'(1,n2,-1+n3,n4,n5,(n6-1),n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)
      +`Z'(1,n2,n3,-1+n4,1+n5,(n6-1),n7,n8,n9,n10,n11,n12,n13,n14)*rat(n5,1)
      +`Z'(1,n2,n3,n4,1+n5,(n6-1),n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(-n5,1)
      +`Z'(1,n2,n3,n4,n5,(n6-1),n7,n8,n9,-1+n10,1+n11,n12,n13,n14)*rat(n11,1)
      +`Z'(1,n2,n3,n4,n5,(n6-1),n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(1,n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,n11,n12,-1+n13,1+n14)*rat(n14,1)
      +`Z'(1,n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,n11,n12,n13,n14)*rat(-6*ep-2*n2-n5-(n6-1)-2*n7-2*n8-n9-2*n10-n11-n12-2*n13-n14+10,1)
      +`Z'(1,-1+n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n11,1)
      +`Z'(1,-1+n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-2*n14,1)
      +`Z'(1,n2,n3,n4,n5,(n6-1),-1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n11,1)
      +`Z'(1,n2,n3,n4,n5,(n6-1),n7,-1+n8,n9,n10,1+n11,n12,n13,n14)*rat(n11,1)
      +`Y'(0,n2,n3,n4,n5,1+(n6-1),n7,n8,n9,n10,n11,n12,n13,n14)*rat(-(n6-1),1)
      +`Y'(0,n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n11,1)
	);
*--#] n6 : 
#endif
#if ( `il' >= 4 )
*--#[ n7 : 
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',1,n7?{>1}
			,n8?`POS8',n9?`POS9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(n5,2*ep*(n7-1)+2*(n7-1)^2+(n7-1)*n10+(n7-1)*n11+(n7-1)*n14-2*(n7-1))*(
      +`Z'(1,n2,n3,-1+n4,1+n5,1,(n7-1),1+n8,n9,n10,n11,n12,n13,n14)*rat(2*n8,1)
      +`Z'(1,n2,n3,n4,1+n5,1,(n7-1),n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep-2*n5-2*n8+2,1)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),1+n8,n9,-1+n10,1+n11,n12,n13,n14)*rat(-n8*n11,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),1+n8,n9,n10,1+n11,-1+n12,n13,n14)*rat(n8*n11,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),1+n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep*n8-2*n5*n8-2*n8^2-n8*n11+2*n8,n5)
      +`Z'(1,n2,n3,1+n4,n5,1,(n7-1),n8,n9,n10,1+n11,n12,n13,n14)*rat(-n4*n11,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),1+n8,n9,n10,1+n11,n12,n13,n14)*rat(-n8*n11,n5)
      +`Z'(1,-1+n2,n3,n4,n5,1,1+(n7-1),n8,n9,1+n10,n11,n12,n13,n14)*rat((n7-1)*n10,n5)
      +`Z'(1,-1+n2,n3,n4,n5,1,(n7-1),1+n8,n9,n10,1+n11,n12,n13,n14)*rat(-n8*n11,n5)
      +`Z'(1,-1+n2,n3,n4,n5,1,(n7-1),n8,1+n9,1+n10,n11,n12,n13,n14)*rat(-n9*n10,n5)
      +`Z'(1,1+n2,-1+n3,n4,n5,1,(n7-1),n8,n9,n10,n11,n12,n13,1+n14)*rat(n2*n14,n5)
      +`Z'(1,1+n2,n3,n4,n5,1,(n7-1),n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-n2*n14,n5)
      +`Z'(1,n2,-1+n3,1+n4,n5,1,(n7-1),n8,n9,n10,1+n11,n12,n13,n14)*rat(n4*n11,n5)
      +`Z'(1,n2,-1+n3,n4,n5,1,(n7-1),n8,1+n9,1+n10,n11,n12,n13,n14)*rat(n9*n10,n5)
      +`Z'(1,n2,n3,-1+n4,1+n5,1,(n7-1),n8,n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)
      +`Z'(1,n2,n3,n4,1+n5,1,(n7-1),-1+n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)
      +`Z'(1,n2,n3,n4,n5,1,-1+(n7-1),1+n8,n9,n10,1+n11,n12,n13,n14)*rat(n8*n11,n5)
      +`Z'(1,n2,n3,n4,n5,1,1+(n7-1),-1+n8,n9,1+n10,n11,n12,n13,n14)*rat((n7-1)*n10,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,1+n10,1+n11,-1+n12,n13,n14)*rat(-n10*n11,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,1+n10,n11,n12,-1+n13,1+n14)*rat(n10*n14,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,1+n10,n11,n12,n13,n14)*rat(-2*n2*n10+n5*n10+(n7-1)*n10+2*n8*n10-n9*n10+n10*n11-n10*n12-
      n10*n13,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,n10,-1+n11,n12,1+n13,1+n14)*rat(n13*n14,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,n10,1+n11,n12,n13,n14)*rat(6*ep*n11+2*n2*n11+2*n3*n11+n4*n11+2*n5*n11+2*(n7-1)*n11+n8*
      n11+2*n9*n11+2*n10*n11+n11^2+n11*n12+n11*n13+2*n11*n14-11*n11,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,n10,n11,-1+n12,1+n13,1+n14)*rat(-n13*n14,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,n10,n11,n12,1+n13,n14)*rat(n13*n14,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,n10,n11,n12,n13,1+n14)*rat(-2*ep*n14-n2*n14-2*n9*n14-n10*n14-n12*n14-n14^2+3*n14,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,1+n10,1+n11,n12,n13,n14)*rat(n10*n11,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12*n14,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13*n14,n5)
      +`Z'(1,-1+n2,n3,n4,n5,1,(n7-1),n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n10*n14,n5)
      +`Z'(1,-1+n2,n3,n4,n5,1,(n7-1),n8,n9,2+n10,n11,n12,n13,n14)*rat(2*n10^2+2*n10,n5)
      +`Z'(1,n2,-1+n3,n4,n5,1,(n7-1),n8,n9,1+n10,n11,n12,n13,1+n14)*rat(-n10*n14,n5)
      +`Z'(1,n2,-1+n3,n4,n5,1,(n7-1),n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12*n14,n5)
      +`Z'(1,n2,-1+n3,n4,n5,1,(n7-1),n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13*n14,n5)
      +`Z'(1,n2,n3,-1+n4,n5,1,(n7-1),n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12*n14,n5)
      +`Z'(1,n2,n3,-1+n4,n5,1,(n7-1),n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13*n14,n5)
      +`Z'(1,n2,n3,n4,-1+n5,1,(n7-1),n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13*n14,n5)
      +`Z'(1,n2,n3,n4,n5,1,-1+(n7-1),n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13*n14,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,-1+n9,1+n10,n11,n12,n13,1+n14)*rat(n10*n14,n5)
      +`Z'(1,n2,n3,n4,n5,1,(n7-1),n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(-n13*n14,n5)
      +`Y'(0,n2,n3,n4,n5,1,1+(n7-1),n8,n9,1+n10,n11,n12,n13,n14)*rat(-(n7-1)*n10,n5)
      +`Y'(1,n2,n3,n4,n5,0,(n7-1),1+n8,n9,n10,1+n11,n12,n13,n14)*rat(n8*n11,n5)
      +`Y'(0,n2,n3,n4,n5,1,(n7-1),n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10*n11,n5)
      +`Y'(1,n2,n3,n4,n5,0,(n7-1),n8,n9,1+n10,1+n11,n12,n13,n14)*rat(-n10*n11,n5)
      +`Y'(1,n2,n3,n4,n5,0,(n7-1),n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13*n14,n5)
	);
*--#] n7 : 
#endif
#if ( `il' >= 5 )
*--#[ n14 : 
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',1,1
			,n8?`POS8',n9?`POS9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?neg_)
			 = -rat(1,n13-1)*(
      +`Z'(1,1+n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,n12,-1+n13,(1+n14))*rat(n2,1)
      +`Z'(1,1+n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,-1+n13,(1+n14))*rat(-n2,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,1+n8,n9,-1+n10,n11,n12,n13,(1+n14))*rat(n8,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,1+n8,n9,n10,n11,n12,-1+n13,(1+n14))*rat(n8,1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,1+n8,n9,n10,n11,n12,-1+n13,(1+n14))*rat(-n8,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,-1+n10,n11,n12,n13,(1+n14))*rat(-n8+1,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,1+n10,-1+n11,n12,-1+n13,(1+n14))*rat(n10,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,-2+n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,-1+n11,-1+n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,-1+n11,n12,-1+n13,1+(1+n14))*rat((1+n14),1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,-1+n11,n12,n13,(1+n14))*rat(-2*ep-n5-2*n8-n10-n11+n13+3,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,n11,-1+n12,n13,(1+n14))*rat(n11-n13,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,-1+n13,(1+n14))*rat(-4*ep-n2-2*n5-n8-2*n9-n10-2*n11-n12-(1+n14)+7,1)
      +`Z'(1,n2,n3,1+n4,n5,1,1,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n4,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,1+n10,n11,n12,-1+n13,(1+n14))*rat(n10,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,-1+n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,n11,1+n12,-1+n13,(1+n14))*rat(n12,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n4-n11+n13-1,1)
      +`Z'(1,-1+n2,n3,-1+n4,n5,1,1,1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n8,1)
      +`Z'(1,-1+n2,n3,n4,-1+n5,1,1,1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(n8,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,1+n10,-1+n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,-1+n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,1+n13,-1+(1+n14))*rat(-n13,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,n13,(1+n14))*rat(4*ep+n2+2*n5+n8+2*n9+n10+n11+n12+(1+n14)-7,1)
      +`Z'(1,1+n2,-1+n3,n4,n5,1,1,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(n2,1)
      +`Z'(1,1+n2,n3,n4,n5,1,1,-1+n8,-1+n9,n10,n11,n12,n13,(1+n14))*rat(-n2,1)
      +`Z'(1,n2,-1+n3,1+n4,n5,1,1,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(n4,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,1+n10,n11,n12,-1+n13,(1+n14))*rat(-n10,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,1+n12,-1+n13,(1+n14))*rat(n12,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n2-n4+n13+1,1)
      +`Z'(1,n2,n3,-1+n4,1+n5,1,1,n8,-1+n9,n10,n11,n12,n13,(1+n14))*rat(-n5,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,n8,n9,1+n10,n11,n12,-1+n13,(1+n14))*rat(-n10,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,n8,n9,n10,-1+n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,n8,n9,n10,n11,1+n12,-1+n13,(1+n14))*rat(-n12,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-4*ep-2*n2-2*n3-n4-n5+n8-2*n9-n10-n12-2*n13-(1+n14)+9,1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,n8,n9,1+n10,n11,n12,-1+n13,(1+n14))*rat(n10,1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,n8,n9,n10,-1+n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,n8,n9,n10,n11,n12,-1+n13,1+(1+n14))*rat((1+n14),1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n8+n13,1)
      +`Z'(1,n2,n3,n4,1+n5,1,1,-1+n8,-1+n9,n10,n11,n12,n13,(1+n14))*rat(n5,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,1+n10,-1+n11,n12,n13,(1+n14))*rat(n10,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,n10,-1+n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,n10,-1+n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,n10,n11,n12,1+n13,-1+(1+n14))*rat(n13,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(n2+2*n3+n4-n5-n8-n11+n13-1,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,-1+n9,1+n10,n11,n12,-1+n13,(1+n14))*rat(n10,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,-1+n13,1+(1+n14))*rat(-(1+n14),1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,n13,(1+n14))*rat(n2-n13,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(-n12,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(n12,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,-1+n2,-1+n3,n4,n5,1,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Z'(1,-1+n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(-n12,1)
      +`Z'(1,-1+n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,-1+n2,n3,-1+n4,n5,1,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(3*n10,1)
      +`Z'(1,-1+n2,n3,-1+n4,n5,1,1,n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(n12,1)
      +`Z'(1,-1+n2,n3,-1+n4,n5,1,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,-1+n2,n3,n4,-1+n5,1,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Z'(1,-1+n2,n3,n4,-1+n5,1,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,-1+n2,n3,n4,-1+n5,1,1,n8,n9,n10,n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,-1+n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-2*n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,-1+n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,-1+n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,-1+n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(n12,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,-1+n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,-1+n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(-n12,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,-1+n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,-1+n8,n9,n10,n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,-1+n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,-1+n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,-1+n9,n10,n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Y'(0,1+n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n2,1)
      +`Y'(0,1+n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,n13,(1+n14))*rat(n2,1)
      +`Y'(0,n2,n3,-1+n4,n5,1,1,1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n8,1)
      +`Y'(0,n2,n3,n4,-1+n5,1,1,1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(n8,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,1+n10,-1+n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,-1+n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,-1+n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,1+n13,-1+(1+n14))*rat(-n13,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,n13,(1+n14))*rat(4*ep+n2+2*n5+n8+2*n9+n10+2*n11+n12+(1+n14)-7,1)
      +`Y'(1,n2,n3,-1+n4,1+n5,0,1,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n5,1)
      +`Y'(1,n2,n3,-1+n4,1+n5,1,0,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n5,1)
      +`Y'(1,n2,n3,n4,1+n5,0,1,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n5,1)
      +`Y'(1,n2,n3,n4,1+n5,1,0,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n5,1)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,1+n10,n11,n12,-1+n13,(1+n14))*rat(-n10,1)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,n10,-1+n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,n10,1+n11,n12,-1+n13,(1+n14))*rat(2*n11,1)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,n10,n11,n12,-1+n13,1+(1+n14))*rat(-(1+n14),1)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n11-n13,1)
      +`Y'(1,n2,n3,n4,n5,1,0,n8,n9,1+n10,n11,n12,-1+n13,(1+n14))*rat(-n10,1)
      +`Y'(1,n2,n3,n4,n5,1,0,n8,n9,n10,-1+n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Y'(1,n2,n3,n4,n5,1,0,n8,n9,n10,n11,n12,-1+n13,1+(1+n14))*rat(-(1+n14),1)
      +`Y'(1,n2,n3,n4,n5,1,0,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n13+1,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(-n12,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(0,n2,-1+n3,n4,n5,1,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Y'(0,n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(-n12,1)
      +`Y'(0,n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(0,n2,n3,-1+n4,n5,1,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Y'(0,n2,n3,-1+n4,n5,1,1,n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(n12,1)
      +`Y'(0,n2,n3,-1+n4,n5,1,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(0,n2,n3,n4,-1+n5,1,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Y'(0,n2,n3,n4,-1+n5,1,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(0,n2,n3,n4,-1+n5,1,1,n8,n9,n10,n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Y'(0,n2,n3,n4,n5,0,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Y'(0,n2,n3,n4,n5,0,1,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-2*n11,1)
      +`Y'(0,n2,n3,n4,n5,0,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(0,n2,n3,n4,n5,0,1,n8,n9,n10,n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Y'(0,n2,n3,n4,n5,1,0,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Y'(0,n2,n3,n4,n5,1,0,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(0,n2,n3,n4,n5,1,0,n8,n9,n10,n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,-1+n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Y'(1,-1+n2,n3,n4,n5,0,1,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Y'(1,-1+n2,n3,n4,n5,0,1,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-2*n11,1)
      +`Y'(1,-1+n2,n3,n4,n5,0,1,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(1,-1+n2,n3,n4,n5,0,1,n8,n9,n10,n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Y'(1,-1+n2,n3,n4,n5,1,0,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Y'(1,-1+n2,n3,n4,n5,1,0,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(1,-1+n2,n3,n4,n5,1,0,n8,n9,n10,n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Y'(1,n2,n3,n4,n5,0,1,-1+n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Y'(1,n2,n3,n4,n5,0,1,-1+n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(2*n11,1)
      +`Y'(1,n2,n3,n4,n5,0,1,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(1,n2,n3,n4,n5,0,1,-1+n8,n9,n10,n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Y'(1,n2,n3,n4,n5,1,0,-1+n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Y'(1,n2,n3,n4,n5,1,0,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(1,n2,n3,n4,n5,1,0,-1+n8,n9,n10,n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      );
*--#] n14 : 
#endif
#if ( `il' >= 6 )
*--#[ n13 : 
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',1,1
			,n8?`POS8',n9?`POS9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?neg_,0)
			 = -rat(1,-4*ep-2*(1+n13)+2)*(
      +`Z'(1,n2,n3,-1+n4,n5,1,1,1+n8,n9,-1+n10,n11,n12,(1+n13),0)*rat(-n8,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,-1+n10,n11,n12,(1+n13),0)*rat(-4*ep-2*n5-2*n8-2*n10-2*n11+6,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,n11,-1+n12,(1+n13),0)*rat(2*ep-2*n3-n4+2*n8-n9+2*n10+2*(1+n13),1)
      +`Z'(1,-1+n2,1+n3,n4,n5,1,1,n8,n9,n10,n11,n12,(1+n13),0)*rat(2*n3,1)
      +`Z'(1,n2,1+n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,(1+n13),0)*rat(-n3,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n9+1,1)
      +`Z'(1,-1+n2,1+n3,-1+n4,n5,1,1,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n3,1)
      +`Z'(1,-1+n2,n3,-1+n4,n5,1,1,1+n8,n9,n10,n11,n12,(1+n13),0)*rat(2*n8,1)
      +`Z'(1,-1+n2,n3,n4,-1+n5,1,1,1+n8,n9,n10,n11,n12,(1+n13),0)*rat(-2*n8,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,1+n10,n11,-1+n12,(1+n13),0)*rat(2*n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,-1+n11,n12,1+(1+n13),0)*rat(2*(1+n13),1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,n11,-1+n12,1+(1+n13),0)*rat(-2*(1+n13),1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,1+(1+n13),-1)*rat(2*(1+n13),1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,(1+n13),0)*rat(-12*ep-2*n2-n4-6*n5-6*n8-3*n9-6*n10-4*n11-2*n12+22,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,n12,(1+n13),0)*rat(4*ep+2*n2+2*n3+3*n4+2*n5+2*n9-11,1)
      +`Z'(1,n2,n3,-1+n4,1+n5,1,1,n8,-1+n9,n10,n11,n12,(1+n13),0)*rat(n5,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,n8,n9,n10,n11,n12,(1+n13),0)*rat(2*ep+2*n2+n3+n9+n10+n12+(1+n13)-5,1)
      +`Z'(1,n2,n3,n4,1+n5,1,1,-1+n8,-1+n9,n10,n11,n12,(1+n13),0)*rat(-n5,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,(1+n13),0)*rat(-2*ep-2*n2-n3-n4-n5-n9+6,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,n11,1+n12,(1+n13),0)*rat(2*n12,1)
      +`Z'(1,-2+n2,n3,n4,n5,1,1,n8,n9,1+n10,n11,n12,(1+n13),0)*rat(-2*n10,1)
      +`Z'(1,-1+n2,-1+n3,n4,n5,1,1,n8,n9,1+n10,n11,n12,(1+n13),0)*rat(-2*n10,1)
      +`Z'(1,-1+n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,1+n12,(1+n13),0)*rat(2*n12,1)
      +`Z'(1,-1+n2,-1+n3,n4,n5,1,1,n8,n9,n10,n11,n12,1+(1+n13),0)*rat(2*(1+n13),1)
      +`Z'(1,-1+n2,n3,-1+n4,n5,1,1,n8,n9,1+n10,n11,n12,(1+n13),0)*rat(-4*n10,1)
      +`Z'(1,-1+n2,n3,-1+n4,n5,1,1,n8,n9,n10,n11,1+n12,(1+n13),0)*rat(-2*n12,1)
      +`Z'(1,-1+n2,n3,-1+n4,n5,1,1,n8,n9,n10,n11,n12,1+(1+n13),0)*rat(-2*(1+n13),1)
      +`Z'(1,-1+n2,n3,n4,-1+n5,1,1,n8,n9,1+n10,n11,n12,(1+n13),0)*rat(2*n10,1)
      +`Z'(1,-1+n2,n3,n4,-1+n5,1,1,n8,n9,n10,n11,n12,1+(1+n13),0)*rat(2*(1+n13),1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,-1+n8,n9,1+n10,n11,n12,(1+n13),0)*rat(2*n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,-1+n9,1+n10,n11,n12,(1+n13),0)*rat(2*n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,-1+n9,n10,n11,n12,1+(1+n13),0)*rat(-2*(1+n13),1)
      +`Y'(0,n2,n3,n4,n5,2,1,n8,n9,n10,n11,-1+n12,(1+n13),0)*rat(1,1)
      +`Y'(0,n2,n3,n4,n5,2,1,n8,n9,n10,n11,n12,-1+(1+n13),0)*rat(-1,1)
      +`Y'(0,-1+n2,n3,n4,n5,1,1,n8,1+n9,n10,n11,n12,(1+n13),0)*rat(-n9,1)
      +`Y'(0,n2,-1+n3,n4,n5,1,1,n8,1+n9,n10,n11,n12,(1+n13),0)*rat(n9,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,-1+n10,1+n11,n12,(1+n13),0)*rat(n11,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,1+n11,n12,-1+(1+n13),0)*rat(-n11,1)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,n11,n12,(1+n13),0)*rat(-2*ep-2*n2-n9-n10-n12-(1+n13)+4,1)
      +`Y'(1,-1+n2,n3,n4,n5,1,0,n8,1+n9,n10,n11,n12,(1+n13),0)*rat(-n9,1)
      +`Y'(1,n2,-1+n3,n4,n5,1,0,n8,1+n9,n10,n11,n12,(1+n13),0)*rat(n9,1)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,-1+n10,1+n11,n12,(1+n13),0)*rat(2*n11,1)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,n10,1+n11,n12,-1+(1+n13),0)*rat(-2*n11,1)
      +`Y'(1,n2,n3,n4,n5,1,0,n8,n9,n10,1+n11,-1+n12,(1+n13),0)*rat(n11,1)
      +`Y'(1,n2,n3,n4,n5,1,0,n8,n9,n10,1+n11,n12,-1+(1+n13),0)*rat(-n11,1)
      +`Y'(1,n2,n3,n4,n5,1,0,n8,n9,n10,n11,n12,(1+n13),0)*rat(-2*ep-2*n2-n9-n10-n12-(1+n13)+4,1)
      +`Y'(0,-1+n2,n3,n4,n5,1,1,n8,n9,n10,1+n11,n12,(1+n13),0)*rat(2*n11,1)
      +`Y'(1,-1+n2,n3,n4,n5,0,1,n8,n9,n10,1+n11,n12,(1+n13),0)*rat(4*n11,1)
      +`Y'(1,-1+n2,n3,n4,n5,1,0,n8,n9,n10,n11,n12,1+(1+n13),0)*rat(-2*(1+n13),1)
      );
*--#] n13 : 
#endif
#if ( `il' >= 7 )
*--#[ n11 : 
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',1,1
			,n8?`POS8',n9?`POS9',n10?`NEG10',n11?neg_,n12?`NEG12',0,0)
			 = -rat(1,-2*ep-2*n2-n9-n12+3)*(
      +`Z'(1,1+n2,-1+n3,n4,n5,1,1,n8,n9,-1+n10,(1+n11),n12,0,0)*rat(2*n2,1)
      +`Z'(1,1+n2,n3,n4,n5,1,1,n8,-1+n9,-1+n10,(1+n11),n12,0,0)*rat(-n2,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,1+n8,n9,-1+n10,(1+n11),n12,0,0)*rat(n8,1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,1+n8,n9,-1+n10,(1+n11),n12,0,0)*rat(-n8,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,-1+n10,(1+n11),n12,0,0)*rat(-2*ep-2*n2-n5-2*n9-n12+5,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,(1+n11),-1+n12,0,0)*rat(4*ep+2*n2+2*n3+n4+n5+2*n9+n12-8,1)
      +`Z'(1,-1+n2,1+n3,n4,n5,1,1,n8,n9,n10,(1+n11),n12,0,0)*rat(-2*n3,1)
      +`Z'(1,n2,1+n3,n4,n5,1,1,n8,-1+n9,n10,(1+n11),n12,0,0)*rat(n3,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,-1+n10,(1+n11),1+n12,0,0)*rat(n12,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,(1+n11),n12,0,0)*rat(-2*ep-2*n2-n12+2,1)
      +`Z'(1,-1+n2,1+n3,-1+n4,n5,1,1,n8,n9,n10,(1+n11),n12,0,0)*rat(n3,1)
      +`Z'(1,-1+n2,n3,n4,-1+n5,1,1,n8,1+n9,n10,(1+n11),n12,0,0)*rat(n9,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,-1+n10,(1+n11),1+n12,0,0)*rat(-n12,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,1+n10,-1+(1+n11),n12,0,0)*rat(2*n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,1+n10,(1+n11),-1+n12,0,0)*rat(-2*n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,-1+(1+n11),1+n12,0,0)*rat(-n12,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,(1+n11),n12,0,0)*rat(-2*n2+n4+n5+2*n8-2*n9+2*n10,1)
      +`Z'(1,n2,-1+n3,n4,-1+n5,1,1,n8,1+n9,n10,(1+n11),n12,0,0)*rat(-n9,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,-1+n10,(1+n11),1+n12,0,0)*rat(2*n12,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,n10,-1+(1+n11),1+n12,0,0)*rat(n12,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,n10,(1+n11),n12,0,0)*rat(-6*ep-2*n3-3*n4-3*n5-2*n8-2*n9-2*n10-n12+14,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,n8,n9,-1+n10,(1+n11),1+n12,0,0)*rat(-n12,1)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,n8,n9,n10,(1+n11),n12,0,0)*rat(-n3-n10+2,1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,n8,n9,n10,(1+n11),n12,0,0)*rat(n10-1,1)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,-1+n9,n10,(1+n11),n12,0,0)*rat(2*ep+n3+n4+n5+n9+n10-5,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,1+n10,(1+n11),n12,0,0)*rat(2*n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,(1+n11),1+n12,0,0)*rat(-n12,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,n10,(1+n11),1+n12,0,0)*rat(n12,1)
      +`Z'(1,-2+n2,n3,n4,n5,1,1,n8,n9,1+n10,(1+n11),n12,0,0)*rat(2*n10,1)
      +`Z'(1,-2+n2,n3,n4,n5,1,1,n8,n9,n10,(1+n11),1+n12,0,0)*rat(-n12,1)
      +`Z'(1,-1+n2,-1+n3,n4,n5,1,1,n8,n9,n10,(1+n11),1+n12,0,0)*rat(n12,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,-1+n8,n9,1+n10,(1+n11),n12,0,0)*rat(-2*n10,1)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,-1+n8,n9,n10,(1+n11),1+n12,0,0)*rat(-n12,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,-1+n8,n9,1+n10,(1+n11),n12,0,0)*rat(2*n10,1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,-1+n8,n9,n10,(1+n11),1+n12,0,0)*rat(n12,1)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,n10,(1+n11),n12,0,0)*rat(2*ep+2*n2+n9+n12-3,1)
      +`Y'(1,n2,n3,n4,n5,1,0,n8,n9,n10,(1+n11),n12,0,0)*rat(2*ep+2*n2+n9+n12-3,1)
      +`Y'(1,-1+n2,n3,n4,n5,0,1,n8,n9,1+n10,(1+n11),n12,0,0)*rat(-2*n10,1)
      +`Y'(1,-1+n2,n3,n4,n5,0,1,n8,n9,n10,(1+n11),1+n12,0,0)*rat(n12,1)
      +`Y'(1,-1+n2,n3,n4,n5,1,0,n8,n9,1+n10,(1+n11),n12,0,0)*rat(-2*n10,1)
      +`Y'(1,-1+n2,n3,n4,n5,1,0,n8,n9,n10,(1+n11),1+n12,0,0)*rat(n12,1)
      +`Y'(1,n2,-1+n3,n4,n5,0,1,n8,n9,n10,(1+n11),1+n12,0,0)*rat(-n12,1)
      +`Y'(1,n2,-1+n3,n4,n5,1,0,n8,n9,n10,(1+n11),1+n12,0,0)*rat(-n12,1)
      );
*--#] n11 : 
#endif
#if ( `il' >= 8 )
*--#[ n12 : 
id,ifmatch->nono1,
	`Z'(1,n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',1,1
			,n8?`POS8',n9?`POS9',n10?`NEG10',0,n12?neg_,0,0)
			 = -rat(1,4*ep+n3+2*n4+2*n5+2*n8+n9+(1+n12)-8)*(
      +`Z'(1,1+n2,-1+n3,n4,n5,1,1,n8,n9,-1+n10,0,(1+n12),0,0)*rat(-n2,1)
      +`Z'(1,1+n2,n3,n4,n5,1,1,n8,-1+n9,-1+n10,0,(1+n12),0,0)*rat(4*ep*n2+2*n2^2-n2*n3-2*n2*n4+2*n2*n8+3*n2*n10-4*n2,6*ep+
      2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,1+n8,n9,-1+n10,0,(1+n12),0,0)*rat(-2*n8,1)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,1+n8,n9,-1+n10,0,(1+n12),0,0)*rat(4*ep*n8+2*n2*n8-n3*n8-2*n4*n8+2*n8^2+3*n8*n10-4*n8,6*ep+
      2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,-1+n10,0,(1+n12),0,0)*rat(-24*ep^2-4*ep*n2-12*ep*n3-24*ep*n4-12*ep*n5-4*ep*n8-12*ep*n9-
      8*ep*n10-16*ep*(1+n12)+80*ep-2*n2*n3-4*n2*n4-2*n2*(1+n12)+6*n2-3*n3*n5-2*n3*n8-3*n3*n9-4*n3*n10-2*n3*(1+n12)+16*n3-6*n4*n5-4*
      n4*n8-6*n4*n9-8*n4*n10-4*n4*(1+n12)+32*n4-2*n5^2-4*n5*n9-n5*n10-4*n5*(1+n12)+20*n5-2*n8*(1+n12)+6*n8-2*n9^2-n9*n10-4*n9*(1+n12)+
      20*n9-4*n10*(1+n12)+14*n10-2*(1+n12)^2+24*(1+n12)-62,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,1+n3,n4,n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(2*ep*n3+n3^2+2*n3*n4+n3*n5+n3*n9+n3*(1+n12)-5*n3,6*ep+2*n2+n5+2*
      n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,1+n4,-1+n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(-n4,1)
      +`Z'(1,n2,n3,1+n4,n5,1,1,-1+n8,n9,n10,0,(1+n12),0,0)*rat(14*ep*n4+4*n2*n4+n3*n4+2*n4^2+3*n4*n5+4*n4*n8+3*n4*n9+6*n4*
      n10+3*n4*(1+n12)-23*n4,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,-1+n10,0,1+(1+n12),0,0)*rat(8*ep*(1+n12)+2*n2*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+2*n5*(1+n12)+2*n8*(1+n12)+2*n9*(1+n12)
      +3*n10*(1+n12)+2*(1+n12)^2-14*(1+n12),6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(-8*ep^2-4*ep*n2-6*ep*n3-10*ep*n4-12*ep*n5-4*ep*n8-6*ep*n9-4*ep*
      n10-6*ep*(1+n12)+44*ep-2*n2*n3-4*n2*n4-4*n2*n5-2*n2*n9-2*n2*(1+n12)+12*n2-n3^2-3*n3*n4-2*n3*n5-2*n3*n8-2*n3*n9-2*n3*n10-2
      *n3*(1+n12)+14*n3-2*n4^2-3*n4*n5-4*n4*n8-3*n4*n9-4*n4*n10-3*n4*(1+n12)+23*n4-2*n5^2-4*n5*n8-3*n5*n9-5*n5*n10-3*n5*(1+n12)+24*
      n5-2*n8*n9-2*n8*(1+n12)+12*n8-n9^2-2*n9*n10-2*n9*(1+n12)+15*n9-2*n10*(1+n12)+13*n10-(1+n12)^2+15*(1+n12)-54,6*ep+2*n2+n5+2*n8+n9+3*
      n10+(1+n12)-9)
      +`Z'(2,-1+n2,n3,n4,n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(-2*ep-n3-2*n4-n5-n9-(1+n12)+5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(2,n2,n3,n4,n5,1,1,-1+n8,n9,n10,0,(1+n12),0,0)*rat(-2*ep-n3-2*n4-n5-n9-(1+n12)+5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,1+n3,-1+n4,n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(-2*ep*n3-n3^2-2*n3*n4-n3*n5-n3*n9-n3*(1+n12)+5*n3,6*ep+2*n2+
      n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,n3,n4,-1+n5,1,1,n8,1+n9,n10,0,(1+n12),0,0)*rat(-2*ep*n9-n3*n9-2*n4*n9-n5*n9-n9^2-n9*(1+n12)+5*n9,6*ep+2*n2+
      n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,-1+n10,0,1+(1+n12),0,0)*rat(8*ep*(1+n12)+2*n2*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+2*n5*(1+n12)+2*n8*(1+n12)+2*n9*
      (1+n12)+3*n10*(1+n12)+2*(1+n12)^2-14*(1+n12),6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,-1,1+(1+n12),0,0)*rat(8*ep*(1+n12)+2*n2*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+2*n5*(1+n12)+2*n8*(1+n12)+2*n9*
      (1+n12)+3*n10*(1+n12)+2*(1+n12)^2-14*(1+n12),6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(-32*ep^2-8*ep*n2-12*ep*n3-24*ep*n4-24*ep*n5-24*ep*n8-16*ep*n9-
      14*ep*n10-22*ep*(1+n12)+122*ep-2*n2*n3-4*n2*n4-4*n2*n5-4*n2*n8-2*n2*n9-4*n2*(1+n12)+16*n2-n3^2-4*n3*n4-4*n3*n5-4*n3*n8-3*
      n3*n9-4*n3*n10-3*n3*(1+n12)+23*n3-4*n4^2-8*n4*n5-8*n4*n8-6*n4*n9-8*n4*n10-6*n4*(1+n12)+46*n4-4*n5^2-8*n5*n8-6*n5*n9-7*n5*
      n10-7*n5*(1+n12)+45*n5-4*n8^2-6*n8*n9-6*n8*n10-8*n8*(1+n12)+44*n8-2*n9^2-4*n9*n10-5*n9*(1+n12)+31*n9-7*n10*(1+n12)+29*n10-3*(1+n12)^2
      +40*(1+n12)-117,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,-1+n3,1+n4,n5,1,1,-1+n8,n9,n10,0,(1+n12),0,0)*rat(-8*ep*n4-2*n2*n4-n3*n4-2*n4^2-2*n4*n5-2*n4*n8-2*n4*n9-3*
      n4*n10-2*n4*(1+n12)+14*n4,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,-1+n3,n4,-1+n5,1,1,n8,1+n9,n10,0,(1+n12),0,0)*rat(2*ep*n9+n3*n9+2*n4*n9+n5*n9+n9^2+n9*(1+n12)-5*n9,6*ep+2*n2+n5
      +2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,-1+n10,0,1+(1+n12),0,0)*rat(-(1+n12),1)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,n10,-1,1+(1+n12),0,0)*rat(-2*ep*(1+n12)-n3*(1+n12)-2*n4*(1+n12)-n5*(1+n12)-n9*(1+n12)-(1+n12)^2+5*(1+n12),6*ep+2*
      n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(8*ep^2+6*ep*n3+20*ep*n4+8*ep*n5+4*ep*n8+6*ep*n9+6*ep*n10+6*ep*
      (1+n12)-50*ep+2*n2*n4+2*n2*n10-4*n2+n3^2+5*n3*n4+3*n3*n5+2*n3*n8+2*n3*n9+2*n3*(1+n12)-14*n3+6*n4^2+8*n4*n5+6*n4*n8+6*n4*
      n9+3*n4*n10+6*n4*(1+n12)-42*n4+2*n5^2+2*n5*n8+3*n5*n9+n5*n10+3*n5*(1+n12)-21*n5+2*n8*n9+2*n8*n10+2*n8*(1+n12)-14*n8+n9^2+n9*
      n10+2*n9*(1+n12)-16*n9+3*n10^2+n10*(1+n12)-15*n10+(1+n12)^2-16*(1+n12)+63,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,-1+n4,1+n5,1,1,n8,-1+n9,n10,0,(1+n12),0,0)*rat(8*ep*n5+2*n2*n5+n3*n5+2*n4*n5+2*n5^2+2*n5*n8+2*n5*n9+3*n5
      *n10+2*n5*(1+n12)-14*n5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,n8,n9,-1+n10,0,1+(1+n12),0,0)*rat(4*ep*(1+n12)+2*n2*(1+n12)-n3*(1+n12)-2*n4*(1+n12)+2*n8*(1+n12)+3*n10*(1+n12)-4*
      (1+n12),6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,-1+n4,n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(44*ep^2+28*ep*n2+28*ep*n3+22*ep*n4+24*ep*n5+12*ep*n8+32*ep*n9+
      30*ep*n10+18*ep*(1+n12)-176*ep+4*n2^2+8*n2*n3+8*n2*n4+8*n2*n5+4*n2*n8+10*n2*n9+10*n2*n10+6*n2*(1+n12)-56*n2+3*n3^2+7*n3*
      n4+7*n3*n5+6*n3*n8+8*n3*n9+9*n3*n10+7*n3*(1+n12)-51*n3+2*n4^2+5*n4*n5+4*n4*n8+7*n4*n9+6*n4*n10+5*n4*(1+n12)-41*n4+3*n5^2+
      4*n5*n8+8*n5*n9+8*n5*n10+5*n5*(1+n12)-46*n5+6*n8*n9+4*n8*n10+2*n8*(1+n12)-28*n8+5*n9^2+11*n9*n10+7*n9*(1+n12)-60*n9+6*n10^2+5
      *n10*(1+n12)-60*n10+2*(1+n12)^2-37*(1+n12)+171,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,n4,-1+n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(-12*ep^2-4*ep*n2-6*ep*n3-6*ep*n4-8*ep*n5-4*ep*n8-8*ep*n9-10*ep
      *n10-2*ep*(1+n12)+46*ep-2*n2*n3-2*n2*n4-2*n2*n5-2*n2*n9-2*n2*n10+10*n2-n3*n5-2*n3*n8-n3*n9-2*n3*n10-n3*(1+n12)+8*n3-n4*n5
      -2*n4*n8-n4*n9-n4*n10-n4*(1+n12)+7*n4-n5^2-2*n5*n8-2*n5*n9-3*n5*n10-n5*(1+n12)+13*n5-2*n8*n9-2*n8*n10+10*n8-n9^2-3*n9*n10
      -n9*(1+n12)+13*n9-3*n10^2+19*n10+4*(1+n12)-40,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,n4,1+n5,1,1,-1+n8,-1+n9,n10,0,(1+n12),0,0)*rat(-8*ep*n5-2*n2*n5-n3*n5-2*n4*n5-2*n5^2-2*n5*n8-2*n5*n9-3*
      n5*n10-2*n5*(1+n12)+14*n5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,1+n10,0,-1+(1+n12),0,0)*rat(-2*n10,1)
      +`Z'(1,n2,n3,n4,n5,1,1,-1+n8,n9,n10,0,(1+n12),0,0)*rat(-32*ep^2-24*ep*n2-20*ep*n3-22*ep*n4-10*ep*n5-8*ep*n8-24*ep*n9-
      14*ep*n10-22*ep*(1+n12)+122*ep-4*n2^2-6*n2*n3-8*n2*n4-4*n2*n5-4*n2*n8-8*n2*n9-6*n2*n10-8*n2*(1+n12)+44*n2-2*n3^2-5*n3*n4-
      5*n3*n5-4*n3*n8-6*n3*n9-7*n3*n10-5*n3*(1+n12)+37*n3-2*n4^2-5*n4*n5-4*n4*n8-7*n4*n9-8*n4*n10-5*n4*(1+n12)+41*n4-n5^2-5*n5*
      n9-n5*n10-4*n5*(1+n12)+22*n5-4*n8*n9-4*n8*(1+n12)+16*n8-4*n9^2-7*n9*n10-7*n9*(1+n12)+45*n9-7*n10*(1+n12)+29*n10-3*(1+n12)^2+40*(1+n12)-
      117,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,n3,n4,n5,1,1,n8,-1+n9,n10,0,(1+n12),0,0)*rat(-4*ep*n10+4*ep-2*n2*n10+2*n2+n3*n10-n3+2*n4*n10-2*n4-2*n8*n10+
      2*n8-3*n10^2+7*n10-4,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,n8,n9,n10,0,1+(1+n12),0,0)*rat(8*ep*(1+n12)+2*n2*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+2*n5*(1+n12)+2*n8*(1+n12)+2*n9*(1+n12)
      +3*n10*(1+n12)+2*(1+n12)^2-14*(1+n12),6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,n8,n9,n10,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)-n3*(1+n12)-2*n4*(1+n12)-n5*(1+n12)-n9*(1+n12)-(1+n12)^2+5*(1+n12),6*ep+2*
      n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-2+n2,n3,n4,n5,1,1,n8,n9,n10,0,1+(1+n12),0,0)*rat(8*ep*(1+n12)+2*n2*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+2*n5*(1+n12)+2*n8*(1+n12)+2*n9*(1+n12)
      +3*n10*(1+n12)+2*(1+n12)^2-14*(1+n12),6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,-1+n3,n4,n5,1,1,n8,n9,n10,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)-n3*(1+n12)-2*n4*(1+n12)-n5*(1+n12)-n9*(1+n12)-(1+n12)^2+5*(1+n12),6*ep+
      2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,n3,-1+n4,n5,1,1,n8,n9,1+n10,0,(1+n12),0,0)*rat(-16*ep*n10-4*n2*n10-2*n3*n10-4*n4*n10-4*n5*n10-4*n8*n10-4
      *n9*n10-6*n10^2-4*n10*(1+n12)+28*n10,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,-1+n8,n9,1+n10,0,(1+n12),0,0)*rat(32*ep*n10+8*n2*n10+4*n3*n10+8*n4*n10+8*n5*n10+8*n8*n10+8*
      n9*n10+12*n10^2+8*n10*(1+n12)-56*n10,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,-1+n2,n3,n4,n5,1,1,-1+n8,n9,n10,0,1+(1+n12),0,0)*rat(8*ep*(1+n12)+2*n2*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+2*n5*(1+n12)+2*n8*(1+n12)+2*n9*
      (1+n12)+3*n10*(1+n12)+2*(1+n12)^2-14*(1+n12),6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,-1+n8,n9,1+n10,0,(1+n12),0,0)*rat(-4*ep*n10-2*n3*n10-4*n4*n10-2*n5*n10-2*n9*n10-2*n10*(1+n12)+
      10*n10,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Z'(1,n2,-1+n3,n4,n5,1,1,-1+n8,n9,n10,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)-n3*(1+n12)-2*n4*(1+n12)-n5*(1+n12)-n9*(1+n12)-(1+n12)^2+5*(1+n12),6*ep+
      2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,1+n2,-1+n3,n4,n5,0,1,n8,n9,-1+n10,0,(1+n12),0,0)*rat(2*ep*n2+n2*n3+2*n2*n4+n2*n5+n2*n9+n2*(1+n12)-5*n2,6*ep+2*n2+
      n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,1+n2,n3,n4,n5,0,1,n8,-1+n9,-1+n10,0,(1+n12),0,0)*rat(-2*ep*n2-n2*n3-2*n2*n4-n2*n5-n2*n9-n2*(1+n12)+5*n2,6*ep+2*n2+
      n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,-1+n4,n5,0,1,1+n8,n9,-1+n10,0,(1+n12),0,0)*rat(2*ep*n8+n3*n8+2*n4*n8+n5*n8+n8*n9+n8*(1+n12)-5*n8,6*ep+2*n2+
      n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,-1+n5,0,1,1+n8,n9,-1+n10,0,(1+n12),0,0)*rat(-2*ep*n8-n3*n8-2*n4*n8-n5*n8-n8*n9-n8*(1+n12)+5*n8,6*ep+2*n2+
      n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,-1+n10,0,(1+n12),0,0)*rat(-12*ep^2-2*ep*n2-6*ep*n3-12*ep*n4-10*ep*n5-2*ep*n8-10*ep*n9-2*
      ep*n10-8*ep*(1+n12)+50*ep-n2*n3-2*n2*n4-n2*n5-n2*n9-n2*(1+n12)+5*n2-2*n3*n5-n3*n8-2*n3*n9-n3*n10-n3*(1+n12)+10*n3-4*n4*n5-2*
      n4*n8-4*n4*n9-2*n4*n10-2*n4*(1+n12)+20*n4-2*n5^2-n5*n8-4*n5*n9-n5*n10-3*n5*(1+n12)+20*n5-n8*n9-n8*(1+n12)+5*n8-2*n9^2-n9*n10-
      3*n9*(1+n12)+20*n9-n10*(1+n12)+5*n10-(1+n12)^2+15*(1+n12)-50,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,n10,-1,(1+n12),0,0)*rat(2*ep*n10-2*ep+n3*n10-n3+2*n4*n10-2*n4+n5*n10-n5+n9*n10-n9+n10*
      (1+n12)-5*n10-(1+n12)+5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(2,n2,n3,n4,n5,-1,1,n8,n9,-1+n10,0,(1+n12),0,0)*rat(-2*ep-n3-2*n4-n5-n9-(1+n12)+5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(2,n2,n3,n4,n5,1,0,n8,n9,n10,0,(1+n12),0,0)*rat(2*ep+n3+2*n4+n5+n9+(1+n12)-5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(0,n2,n3,n4,n5,1,1,n8,n9,n10,0,(1+n12),0,0)*rat(-4*ep^2-2*ep*n3-4*ep*n4-2*ep*n5-2*ep*n9-2*ep*(1+n12)+12*ep+n3+2*n4+n5
      +n9+(1+n12)-5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,-1+n2,n3,n4,n5,0,1,n8,1+n9,n10,0,(1+n12),0,0)*rat(2*ep*n9+n3*n9+2*n4*n9+n5*n9+n9^2+n9*(1+n12)-5*n9,6*ep+2*n2+n5+2*
      n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,-1+n2,n3,n4,n5,1,0,n8,1+n9,n10,0,(1+n12),0,0)*rat(2*ep*n9+n3*n9+2*n4*n9+n5*n9+n9^2+n9*(1+n12)-5*n9,6*ep+2*n2+n5+2*
      n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,-1+n3,n4,n5,0,1,n8,1+n9,n10,0,(1+n12),0,0)*rat(-2*ep*n9-n3*n9-2*n4*n9-n5*n9-n9^2-n9*(1+n12)+5*n9,6*ep+2*n2+n5+2
      *n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,-1+n3,n4,n5,1,0,n8,1+n9,n10,0,(1+n12),0,0)*rat(-2*ep*n9-n3*n9-2*n4*n9-n5*n9-n9^2-n9*(1+n12)+5*n9,6*ep+2*n2+n5+2
      *n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,-1+n4,1+n5,0,1,n8,n9,n10,0,(1+n12),0,0)*rat(-2*ep*n5-n3*n5-2*n4*n5-n5^2-n5*n9-n5*(1+n12)+5*n5,6*ep+2*n2+n5+2
      *n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,-1+n4,1+n5,1,0,n8,n9,n10,0,(1+n12),0,0)*rat(-2*ep*n5-n3*n5-2*n4*n5-n5^2-n5*n9-n5*(1+n12)+5*n5,6*ep+2*n2+n5+2
      *n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,1+n5,0,1,-1+n8,n9,n10,0,(1+n12),0,0)*rat(2*ep*n5+n3*n5+2*n4*n5+n5^2+n5*n9+n5*(1+n12)-5*n5,6*ep+2*n2+n5+2*
      n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,1+n5,1,0,-1+n8,n9,n10,0,(1+n12),0,0)*rat(2*ep*n5+n3*n5+2*n4*n5+n5^2+n5*n9+n5*(1+n12)-5*n5,6*ep+2*n2+n5+2*
      n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,-1+n10,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+n5*(1+n12)+n9*(1+n12)+(1+n12)^2-5*(1+n12),6*ep+2*n2
      +n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,n9,n10,0,(1+n12),0,0)*rat(2*ep+n3+2*n4+n5+n9+(1+n12)-5,1)
      +`Y'(1,n2,n3,n4,n5,1,0,n8,n9,n10,0,(1+n12),0,0)*rat(12*ep^2+4*ep*n2+6*ep*n3+12*ep*n4+8*ep*n5+4*ep*n8+8*ep*n9+4*ep*n10
      +8*ep*(1+n12)-48*ep+2*n2*n3+4*n2*n4+2*n2*n5+2*n2*n9+2*n2*(1+n12)-10*n2+n3*n5+2*n3*n8+n3*n9+2*n3*n10+n3*(1+n12)-9*n3+2*n4*n5+4
      *n4*n8+2*n4*n9+4*n4*n10+2*n4*(1+n12)-18*n4+n5^2+2*n5*n8+2*n5*n9+2*n5*n10+2*n5*(1+n12)-14*n5+2*n8*n9+2*n8*(1+n12)-10*n8+n9^2+2
      *n9*n10+2*n9*(1+n12)-14*n9+2*n10*(1+n12)-10*n10+(1+n12)^2-14*(1+n12)+45,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,-1+n3,n4,n5,0,1,n8,n9,-1+n10,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+n5*(1+n12)+n9*(1+n12)+(1+n12)^2-5*(1+n12),6*ep+2
      *n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,-1+n3,n4,n5,0,1,n8,n9,n10,0,(1+n12),0,0)*rat(-2*ep*n10+2*ep-n3*n10+n3-2*n4*n10+2*n4-n5*n10+n5-n9*n10+n9-n10
      *(1+n12)+5*n10+(1+n12)-5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,-1+n4,n5,0,1,n8,n9,-1+n10,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)-n3*(1+n12)-2*n4*(1+n12)-n5*(1+n12)-n9*(1+n12)-(1+n12)^2+5*(1+n12),6*ep+
      2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,-1+n4,n5,0,1,n8,n9,n10,0,(1+n12),0,0)*rat(-2*ep*n10+2*ep-n3*n10+n3-2*n4*n10+2*n4-n5*n10+n5-n9*n10+n9-n10
      *(1+n12)+5*n10+(1+n12)-5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,-1+n5,0,1,n8,n9,n10,0,(1+n12),0,0)*rat(2*ep*n10-2*ep+n3*n10-n3+2*n4*n10-2*n4+n5*n10-n5+n9*n10-n9+n10*
      (1+n12)-5*n10-(1+n12)+5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,n5,-1,1,n8,n9,n10,0,(1+n12),0,0)*rat(-2*ep*n10+2*ep-n3*n10+n3-2*n4*n10+2*n4-n5*n10+n5-n9*n10+n9-n10*
      (1+n12)+5*n10+(1+n12)-5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,n5,0,0,n8,n9,n10,0,(1+n12),0,0)*rat(-2*ep*n10+2*ep-n3*n10+n3-2*n4*n10+2*n4-n5*n10+n5-n9*n10+n9-n10*
      (1+n12)+5*n10+(1+n12)-5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,n3,n4,n5,0,1,n8,-1+n9,n10,0,(1+n12),0,0)*rat(2*ep*n10-2*ep+n3*n10-n3+2*n4*n10-2*n4+n5*n10-n5+n9*n10-n9+n10*
      (1+n12)-5*n10-(1+n12)+5,6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,-1+n2,n3,n4,n5,0,1,n8,n9,n10,0,1+(1+n12),0,0)*rat(-8*ep*(1+n12)-2*n2*(1+n12)-n3*(1+n12)-2*n4*(1+n12)-2*n5*(1+n12)-2*n8*(1+n12)-2*n9*
      (1+n12)-3*n10*(1+n12)-2*(1+n12)^2+14*(1+n12),6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,-1+n2,n3,n4,n5,1,0,n8,n9,n10,0,1+(1+n12),0,0)*rat(-8*ep*(1+n12)-2*n2*(1+n12)-n3*(1+n12)-2*n4*(1+n12)-2*n5*(1+n12)-2*n8*(1+n12)-2*n9*
      (1+n12)-3*n10*(1+n12)-2*(1+n12)^2+14*(1+n12),6*ep+2*n2+n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,-1+n3,n4,n5,0,1,n8,n9,n10,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+n5*(1+n12)+n9*(1+n12)+(1+n12)^2-5*(1+n12),6*ep+2*n2
      +n5+2*n8+n9+3*n10+(1+n12)-9)
      +`Y'(1,n2,-1+n3,n4,n5,1,0,n8,n9,n10,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)+n3*(1+n12)+2*n4*(1+n12)+n5*(1+n12)+n9*(1+n12)+(1+n12)^2-5*(1+n12),6*ep+2*n2
      +n5+2*n8+n9+3*n10+(1+n12)-9)
      );
*--#] n12 : 
#endif
#if ( `il' >= 9 )
*--#[ n3 : 
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,n2?`POS2',n3?{>1},n4?`POS4',n5?`POS5',1,1
			,n8?`POS8',n9?`POS9',n10?`NEG10',0,0,0,0)
			 = -rat(n8,-2*ep*(n3-1)-2*n2*(n3-1)-2*(n3-1)^2+2*(n3-1))*(
      +`Z'(1,1+n2,1+(n3-1),n4,n5,1,1,n8,-1+n9,n10,0,0,0,0)*rat(2*n2*(n3-1),n8)
      +`Z'(1,1+n2,(n3-1),n4,n5,1,1,n8,n9,n10,0,0,0,0)*rat(-2*ep*n2-2*n2^2-2*n2*(n3-1)+2*n2,n8)
      +`Z'(1,n2,(n3-1),1+n4,-1+n5,1,1,1+n8,n9,n10,0,0,0,0)*rat(-2*n4,1)
      +`Z'(1,n2,(n3-1),1+n4,n5,1,1,n8,n9,n10,0,0,0,0)*rat(2*ep*n4+2*n4^2+2*n4*n8-2*n4,n8)
      +`Z'(1,n2,(n3-1),n4,n5,1,1,1+n8,n9,n10,0,0,0,0)*rat(2*ep+2*n4+2*n8-2,1)
      +`Z'(1,-1+n2,1+(n3-1),n4,n5,1,1,n8,n9,1+n10,0,0,0,0)*rat(2*(n3-1)*n10,n8)
      +`Z'(1,n2,1+(n3-1),n4,n5,1,1,n8,-1+n9,1+n10,0,0,0,0)*rat(-2*(n3-1)*n10,n8)
      +`Z'(1,n2,(n3-1),1+n4,-1+n5,1,1,n8,n9,1+n10,0,0,0,0)*rat(2*n4*n10,n8)
      +`Z'(1,n2,(n3-1),1+n4,n5,1,1,-1+n8,n9,1+n10,0,0,0,0)*rat(-2*n4*n10,n8)
      +`Z'(1,n2,(n3-1),n4,n5,1,1,n8,n9,1+n10,0,0,0,0)*rat(4*n2*n10+2*(n3-1)*n10-2*n4*n10-4*n8*n10,n8)
      +`Z'(1,-1+n2,(n3-1),n4,n5,1,1,n8,n9,2+n10,0,0,0,0)*rat(-2*n10^2-2*n10,n8)
      +`Z'(1,n2,(n3-1),n4,n5,1,1,-1+n8,n9,2+n10,0,0,0,0)*rat(2*n10^2+2*n10,n8)
      );
*--#] n3 : 
#endif
#if ( `il' >= 10 )
*--#[ n9 : 
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,n2?`POS3',1,n4?`POS4',n5?`POS5',1,1
			,n8?`POS8',n9?{>1},n10?`NEG10',0,0,0,0)
			 = -rat(1,2*ep+2*n2+2*(n9-1)-2)*(
      +`Z'(1,1+n2,1,n4,n5,1,1,n8,(n9-1),n10,0,0,0,0)*rat(2*ep*n2+2*n2^2+2*n2*(n9-1)-2*n2,(n9-1))
      +`Z'(1,n2,1,-1+n4,1+n5,1,1,1+n8,(n9-1),n10,0,0,0,0)*rat(4*ep*n5*n8,2*ep*(n9-1)+(n9-1)*n10)
      +`Z'(1,n2,1,n4,1+n5,1,1,n8,(n9-1),n10,0,0,0,0)*rat(-4*ep^2*n5-4*ep*n5^2-4*ep*n5*n8+4*ep*n5,2*ep*(n9-1)+(n9-1)*n10)
      +`Z'(1,n2,1,n4,n5,1,1,1+n8,(n9-1),n10,0,0,0,0)*rat(-4*ep^2*n8-4*ep*n5*n8-4*ep*n8^2+4*ep*n8,2*ep*(n9-1)+(n9-1)*n10)
      +`Z'(1,-1+n2,1,n4,n5,1,1,n8,1+(n9-1),1+n10,0,0,0,0)*rat(-2*ep*n10,2*ep+n10)
      +`Z'(1,-1+n2,1,n4,n5,1,2,n8,(n9-1),1+n10,0,0,0,0)*rat(2*ep*n10,2*ep*(n9-1)+(n9-1)*n10)
      +`Z'(1,n2,1,-1+n4,1+n5,1,1,n8,(n9-1),1+n10,0,0,0,0)*rat(-6*ep*n5*n10-2*n5*n10^2,2*ep*(n9-1)+(n9-1)*n10)
      +`Z'(1,n2,1,n4,1+n5,1,1,-1+n8,(n9-1),1+n10,0,0,0,0)*rat(6*ep*n5*n10+2*n5*n10^2,2*ep*(n9-1)+(n9-1)*n10)
      +`Z'(1,n2,1,n4,n5,1,1,n8,(n9-1),1+n10,0,0,0,0)*rat(8*ep^2*n10-4*ep*n2*n10+6*ep*n5*n10+12*ep*n8*n10-2*ep*(n9-1)*n10+8*ep*
      n10^2-10*ep*n10+2*n5*n10^2+4*n8*n10^2+2*n10^3-6*n10^2,2*ep*(n9-1)+(n9-1)*n10)
      +`Z'(1,n2,1,n4,n5,1,2,-1+n8,(n9-1),1+n10,0,0,0,0)*rat(2*ep*n10,2*ep*(n9-1)+(n9-1)*n10)
      +`Z'(1,-1+n2,1,n4,n5,1,1,n8,(n9-1),2+n10,0,0,0,0)*rat(4*ep*n10^2+4*ep*n10,2*ep*(n9-1)+(n9-1)*n10)
      +`Z'(1,n2,1,n4,n5,1,1,-1+n8,(n9-1),2+n10,0,0,0,0)*rat(-2*n10^2-2*n10,(n9-1))
      +`Y'(1,1+n2,0,n4,n5,1,1,n8,1+(n9-1),n10,0,0,0,0)*rat(-2*n2,1)
      +`Y'(0,n2,1,n4,n5,1,2,n8,(n9-1),1+n10,0,0,0,0)*rat(-2*ep*n10,2*ep*(n9-1)+(n9-1)*n10)
      +`Y'(1,n2,0,n4,n5,1,1,n8,1+(n9-1),1+n10,0,0,0,0)*rat(2*ep*n10,2*ep+n10)
      );
*--#] n9 : 
#endif
#if ( `il' >= 11 )
*--#[ n2 : 
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,n2?{>1},1,n4?`POS4',n5?`POS5',1,1,n8?`POS8',1,n10?`NEG10',0,0,0,0)
			 = -rat(ep+(n2-1)-1,-2*ep^2*(n2-1)-3*ep*(n2-1)^2+ep*(n2-1)-(n2-1)^3+(n2-1)^2)*(
      +`Z'(1,-1+(n2-1),1,-1+n4,1+n5,1,1,2+n8,1,n10,0,0,0,0)*rat(2*ep^2*n5*n8^2+2*ep^2*n5*n8+2*ep*n4*n5*n8^2+2*ep*n4*n5*n8+
      2*ep*n5*n8^3-2*ep*n5*n8,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*
      (n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,1+n4,-1+n5,1,1,2+n8,1,n10,0,0,0,0)*rat(2*ep^2*n4*n8^2+2*ep^2*n4*n8+2*ep*n4*n5*n8^2+2*ep*n4*n5*n8+
      2*ep*n4*n8^3-2*ep*n4*n8,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*
      (n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,1+n4,1+n5,1,1,n8,1,n10,0,0,0,0)*rat(-2*ep^3*n4*n5-2*ep^2*n4^2*n5-2*ep^2*n4*n5^2-4*ep^2*n4*n5*n8+4
      *ep^2*n4*n5-2*ep*n4^2*n5^2-2*ep*n4^2*n5*n8+2*ep*n4^2*n5-2*ep*n4*n5^2*n8+2*ep*n4*n5^2-2*ep*n4*n5*n8^2+4*ep*n4*n5*
      n8-2*ep*n4*n5,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,1+n4,n5,1,1,1+n8,1,n10,0,0,0,0)*rat(-2*ep^3*n4*n8-2*ep^2*n4^2*n8-4*ep^2*n4*n8^2+2*ep^2*n4*n8-2*ep
      *n4^2*n5*n8-2*ep*n4^2*n8^2+2*ep*n4^2*n8+2*ep*n4*n5^2*n8-2*ep*n4*n5*n8-2*ep*n4*n8^3+2*ep*n4*n8^2,2*ep^3+4*ep^2*(n2-1)+
      ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,n4,1+n5,1,1,1+n8,1,n10,0,0,0,0)*rat(-2*ep^3*n5*n8-2*ep^2*n5^2*n8-4*ep^2*n5*n8^2+2*ep^2*n5*n8+2*ep
      *n4^2*n5*n8-2*ep*n4*n5^2*n8-2*ep*n4*n5*n8-2*ep*n5^2*n8^2+2*ep*n5^2*n8-2*ep*n5*n8^3+2*ep*n5*n8^2,2*ep^3+4*ep^2*(n2-1)+
      ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,n4,n5,1,1,2+n8,1,n10,0,0,0,0)*rat(-2*ep^3*n8^2-2*ep^3*n8-2*ep^2*n4*n8^2-2*ep^2*n4*n8-2*ep^2*n5*
      n8^2-2*ep^2*n5*n8-4*ep^2*n8^3-2*ep^2*n8^2+2*ep^2*n8-4*ep*n4*n5*n8^2-4*ep*n4*n5*n8-2*ep*n4*n8^3+2*ep*n4*n8-2*ep*n5
      *n8^3+2*ep*n5*n8-2*ep*n8^4+2*ep*n8^2,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep
      +(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,(n2-1),1,-1+n4,1+n5,1,1,1+n8,1,n10,0,0,0,0)*rat(-n5*n8,1)
      +`Z'(1,(n2-1),1,1+n4,-1+n5,1,1,1+n8,1,n10,0,0,0,0)*rat(-n4*n8,1)
      +`Z'(1,(n2-1),1,1+n4,1+n5,1,1,-1+n8,1,n10,0,0,0,0)*rat(n4*n5,1)
      +`Z'(1,(n2-1),1,1+n4,n5,1,1,n8,1,n10,0,0,0,0)*rat(-n4*n5+n4*n8,1)
      +`Z'(1,(n2-1),1,n4,1+n5,1,1,n8,1,n10,0,0,0,0)*rat(-n4*n5+n5*n8,1)
      +`Z'(1,(n2-1),1,n4,n5,1,1,1+n8,1,n10,0,0,0,0)*rat(2*ep*n8+n4*n8+n5*n8+2*n8^2-2*n8,1)
      +`Z'(1,-2+(n2-1),1,1+n4,-1+n5,1,1,1+n8,2,1+n10,0,0,0,0)*rat(ep*n4*n8*n10,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2
      *ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-2+(n2-1),1,1+n4,-1+n5,1,2,1+n8,1,1+n10,0,0,0,0)*rat(-ep*n4*n8*n10,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+
      2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-2+(n2-1),1,1+n4,n5,1,1,n8,2,1+n10,0,0,0,0)*rat(-ep^2*n4*n10-ep*n4^2*n10-ep*n4*n8*n10+ep*n4*n10,2*ep^3+4*ep^2*
      (n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-2+(n2-1),1,1+n4,n5,1,2,n8,1,1+n10,0,0,0,0)*rat(ep^2*n4*n10+ep*n4^2*n10+ep*n4*n8*n10-ep*n4*n10,2*ep^3+4*ep^2*
      (n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-2+(n2-1),1,n4,n5,1,1,1+n8,2,1+n10,0,0,0,0)*rat(-ep^2*n8*n10-ep*n4*n8*n10-ep*n8^2*n10+ep*n8*n10,2*ep^3+4*ep^2*
      (n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-2+(n2-1),1,n4,n5,1,2,1+n8,1,1+n10,0,0,0,0)*rat(ep^2*n8*n10+ep*n4*n8*n10+ep*n8^2*n10-ep*n8*n10,2*ep^3+4*ep^2*
      (n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-2+(n2-1),2,n4,n5,1,1,n8,2,1+n10,0,0,0,0)*rat(-n10,ep+(n2-1)-1)
      +`Z'(1,-1+(n2-1),1,-1+n4,1+n5,1,1,1+n8,1,1+n10,0,0,0,0)*rat(-3*ep^2*n5*n8*n10-3*ep*n4*n5*n8*n10-3*ep*n5*n8^2*n10-ep*
      n5*n8*n10^2+3*ep*n5*n8*n10-n4*n5*n8*n10^2-n5*n8^2*n10^2+n5*n8*n10^2,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*
      ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,1+n4,-1+n5,1,1,1+n8,1,1+n10,0,0,0,0)*rat(-4*ep^2*n4*n8*n10+2*ep*(n2-1)*n4*n8*n10-3*ep*n4*n5*n8*n10-6*
      ep*n4*n8^2*n10-4*ep*n4*n8*n10^2+ep*n4*n8*n10-n4*n5*n8*n10^2-2*n4*n8^2*n10^2-n4*n8*n10^3+2*n4*n8*n10^2,2*ep^3+4*
      ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,1+n4,-1+n5,1,1,n8,2,1+n10,0,0,0,0)*rat(-n4*n10,ep+(n2-1)-1)
      +`Z'(1,-1+(n2-1),1,1+n4,-1+n5,1,2,n8,1,1+n10,0,0,0,0)*rat(-ep*n4*n8*n10,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*
      ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,1+n4,1+n5,1,1,-1+n8,1,1+n10,0,0,0,0)*rat(3*ep^2*n4*n5*n10+3*ep*n4^2*n5*n10+3*ep*n4*n5*n8*n10+ep*n4
      *n5*n10^2-3*ep*n4*n5*n10+n4^2*n5*n10^2+n4*n5*n8*n10^2-n4*n5*n10^2,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep
      *(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,1+n4,n5,1,1,-1+n8,2,1+n10,0,0,0,0)*rat(n4*n10,ep+(n2-1)-1)
      +`Z'(1,-1+(n2-1),1,1+n4,n5,1,1,n8,1,1+n10,0,0,0,0)*rat(4*ep^3*n4*n10-2*ep^2*(n2-1)*n4*n10+4*ep^2*n4^2*n10+3*ep^2*n4*n5*
      n10+10*ep^2*n4*n8*n10+4*ep^2*n4*n10^2-8*ep^2*n4*n10-2*ep*(n2-1)*n4^2*n10-2*ep*(n2-1)*n4*n8*n10+2*ep*(n2-1)*n4*n10+3*ep*n4^2*
      n5*n10+6*ep*n4^2*n8*n10+4*ep*n4^2*n10^2-4*ep*n4^2*n10+ep*n4*n5*n10^2-3*ep*n4*n5*n10+6*ep*n4*n8^2*n10+6*ep*n4*n8*
      n10^2-7*ep*n4*n8*n10+ep*n4*n10^3-7*ep*n4*n10^2+4*ep*n4*n10+n4^2*n5*n10^2+2*n4^2*n8*n10^2+n4^2*n10^3-3*n4^2*n10^2-
      n4*n5*n10^2+2*n4*n8^2*n10^2+n4*n8*n10^3-4*n4*n8*n10^2-n4*n10^3+3*n4*n10^2,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*
      (n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,1+n4,n5,1,2,-1+n8,1,1+n10,0,0,0,0)*rat(ep^2*n4*n10+ep*n4^2*n10+ep*n4*n8*n10-ep*n4*n10,2*ep^3+4*
      ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,n4,1+n5,1,1,n8,1,1+n10,0,0,0,0)*rat(-3*ep^2*n4*n5*n10+3*ep^2*n5*n8*n10-3*ep*n4^2*n5*n10-ep*n4*n5*
      n10^2+3*ep*n4*n5*n10+3*ep*n5*n8^2*n10+ep*n5*n8*n10^2-3*ep*n5*n8*n10-n4^2*n5*n10^2+n4*n5*n10^2+n5*n8^2*n10^2-n5*n8
      *n10^2,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,n4,n5,1,1,1+n8,1,1+n10,0,0,0,0)*rat(4*ep^3*n8*n10-2*ep^2*(n2-1)*n8*n10+4*ep^2*n4*n8*n10+3*ep^2*n5*n8*
      n10+10*ep^2*n8^2*n10+4*ep^2*n8*n10^2-2*ep^2*n8*n10-2*ep*(n2-1)*n4*n8*n10-2*ep*(n2-1)*n8^2*n10+2*ep*(n2-1)*n8*n10+6*ep*n4*n5*
      n8*n10+6*ep*n4*n8^2*n10+4*ep*n4*n8*n10^2-ep*n4*n8*n10+3*ep*n5*n8^2*n10+ep*n5*n8*n10^2-3*ep*n5*n8*n10+6*ep*n8^3*
      n10+6*ep*n8^2*n10^2-4*ep*n8^2*n10+ep*n8*n10^3-5*ep*n8*n10^2-2*ep*n8*n10+2*n4*n5*n8*n10^2+2*n4*n8^2*n10^2+n4*n8*
      n10^3-2*n4*n8*n10^2+n5*n8^2*n10^2-n5*n8*n10^2+2*n8^3*n10^2+n8^2*n10^3-3*n8^2*n10^2-n8*n10^3+n8*n10^2,2*ep^3+4*
      ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),1,n4,n5,1,1,n8,2,1+n10,0,0,0,0)*rat(ep*n10-(n2-1)*n10+n4*n10+2*n8*n10,ep+(n2-1)-1)
      +`Z'(1,-1+(n2-1),1,n4,n5,1,2,n8,1,1+n10,0,0,0,0)*rat(ep^2*n8*n10+ep*n4*n8*n10+ep*n8^2*n10-ep*n8*n10,2*ep^3+4*ep^2*(n2-1)+
      ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-1+(n2-1),2,n4,n5,1,1,n8,1,1+n10,0,0,0,0)*rat(ep*n10+n10,ep+(n2-1)-1)
      +`Z'(1,(n2-1),1,-1+n4,1+n5,1,1,n8,1,1+n10,0,0,0,0)*rat(n5*n10,1)
      +`Z'(1,(n2-1),1,1+n4,-1+n5,1,1,n8,1,1+n10,0,0,0,0)*rat(ep*n4*n10,ep+(n2-1)-1)
      +`Z'(1,(n2-1),1,1+n4,n5,1,1,-1+n8,1,1+n10,0,0,0,0)*rat(-ep*n4*n10,ep+(n2-1)-1)
      +`Z'(1,(n2-1),1,n4,1+n5,1,1,-1+n8,1,1+n10,0,0,0,0)*rat(-n5*n10,1)
      +`Z'(1,(n2-1),1,n4,n5,1,1,n8,1,1+n10,0,0,0,0)*rat(4*ep*(n2-1)*n10-ep*n4*n10-ep*n5*n10-4*ep*n8*n10+2*ep*n10+2*(n2-1)^2*n10-(n2-1)*
      n5*n10-2*(n2-1)*n8*n10-(n2-1)*n10+n5*n10+2*n8*n10-n10,ep+(n2-1)-1)
      +`Z'(1,-2+(n2-1),1,1+n4,-1+n5,1,1,1+n8,1,2+n10,0,0,0,0)*rat(-2*ep*n4*n8*n10^2-2*ep*n4*n8*n10,2*ep^3+4*ep^2*(n2-1)+ep^2*
      n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-2+(n2-1),1,1+n4,n5,1,1,n8,1,2+n10,0,0,0,0)*rat(2*ep^2*n4*n10^2+2*ep^2*n4*n10+2*ep*n4^2*n10^2+2*ep*n4^2*n10+2*
      ep*n4*n8*n10^2+2*ep*n4*n8*n10-2*ep*n4*n10^2-2*ep*n4*n10,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*
      ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-2+(n2-1),1,n4,n5,1,1,1+n8,1,2+n10,0,0,0,0)*rat(2*ep^2*n8*n10^2+2*ep^2*n8*n10+2*ep*n4*n8*n10^2+2*ep*n4*n8*n10+
      2*ep*n8^2*n10^2+2*ep*n8^2*n10-2*ep*n8*n10^2-2*ep*n8*n10,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*
      ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Z'(1,-2+(n2-1),1,n4,n5,1,1,n8,2,2+n10,0,0,0,0)*rat(n10^2+n10,ep+(n2-1)-1)
      +`Z'(1,-1+(n2-1),1,1+n4,-1+n5,1,1,n8,1,2+n10,0,0,0,0)*rat(n4*n8*n10^2+n4*n8*n10,ep^2+2*ep*(n2-1)-2*ep+(n2-1)^2-2*(n2-1)+1)
      +`Z'(1,-1+(n2-1),1,1+n4,n5,1,1,-1+n8,1,2+n10,0,0,0,0)*rat(-ep*n4*n10^2-ep*n4*n10-n4^2*n10^2-n4^2*n10-n4*n8*n10^2-n4*
      n8*n10+n4*n10^2+n4*n10,ep^2+2*ep*(n2-1)-2*ep+(n2-1)^2-2*(n2-1)+1)
      +`Z'(1,-1+(n2-1),1,n4,n5,1,1,-1+n8,2,2+n10,0,0,0,0)*rat(-n10^2-n10,ep+(n2-1)-1)
      +`Z'(1,-1+(n2-1),1,n4,n5,1,1,n8,1,2+n10,0,0,0,0)*rat(-2*ep^2*n10^2-2*ep^2*n10-3*ep*(n2-1)*n10^2-3*ep*(n2-1)*n10-ep*n8*n10^2-
      ep*n8*n10+3*ep*n10^2+3*ep*n10-(n2-1)^2*n10^2-(n2-1)^2*n10+2*(n2-1)*n10^2+2*(n2-1)*n10-n4*n8*n10^2-n4*n8*n10-n8^2*n10^2-n8^2*n10+
      n8*n10^2+n8*n10-n10^2-n10,ep^2+2*ep*(n2-1)-2*ep+(n2-1)^2-2*(n2-1)+1)
      +`Z'(1,(n2-1),1,n4,n5,1,1,-1+n8,1,2+n10,0,0,0,0)*rat(2*ep*n10^2+2*ep*n10+(n2-1)*n10^2+(n2-1)*n10-n10^2-n10,ep+(n2-1)-1)
      +`Y'(1,1+(n2-1),0,n4,n5,1,1,n8,2,n10,0,0,0,0)*rat((n2-1),1)
      +`Y'(1,1+(n2-1),2,n4,n5,1,1,n8,0,n10,0,0,0,0)*rat(ep*(n2-1),ep+(n2-1)-1)
      +`Y'(1,(n2-1),0,1+n4,-1+n5,1,1,1+n8,2,n10,0,0,0,0)*rat((n2-1)*n4*n8-n4*n8,ep^2+2*ep*(n2-1)-2*ep+(n2-1)^2-2*(n2-1)+1)
      +`Y'(1,(n2-1),0,1+n4,n5,1,1,n8,2,n10,0,0,0,0)*rat(-ep*(n2-1)*n4+ep*n4-(n2-1)*n4^2-(n2-1)*n4*n8+(n2-1)*n4+n4^2+n4*n8-n4,ep^2+2*ep*(n2-1)-
      2*ep+(n2-1)^2-2*(n2-1)+1)
      +`Y'(1,(n2-1),0,n4,n5,1,1,1+n8,2,n10,0,0,0,0)*rat(-ep*(n2-1)*n8+ep*n8-(n2-1)*n4*n8-(n2-1)*n8^2+(n2-1)*n8+n4*n8+n8^2-n8,ep^2+2*ep*(n2-1)-
      2*ep+(n2-1)^2-2*(n2-1)+1)
      +`Y'(0,-1+(n2-1),1,1+n4,-1+n5,1,2,1+n8,1,1+n10,0,0,0,0)*rat(ep*n4*n8*n10,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2
      *ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Y'(0,-1+(n2-1),1,1+n4,n5,1,2,n8,1,1+n10,0,0,0,0)*rat(-ep^2*n4*n10-ep*n4^2*n10-ep*n4*n8*n10+ep*n4*n10,2*ep^3+4*ep^2*
      (n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Y'(0,-1+(n2-1),1,n4,n5,1,2,1+n8,1,1+n10,0,0,0,0)*rat(-ep^2*n8*n10-ep*n4*n8*n10-ep*n8^2*n10+ep*n8*n10,2*ep^3+4*ep^2*
      (n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Y'(1,-1+(n2-1),0,1+n4,-1+n5,1,1,1+n8,2,1+n10,0,0,0,0)*rat(-ep*n4*n8*n10,2*ep^3+4*ep^2*(n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+
      2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Y'(1,-1+(n2-1),0,1+n4,n5,1,1,n8,2,1+n10,0,0,0,0)*rat(ep^2*n4*n10+ep*n4^2*n10+ep*n4*n8*n10-ep*n4*n10,2*ep^3+4*ep^2*
      (n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Y'(1,-1+(n2-1),0,n4,n5,1,1,1+n8,2,1+n10,0,0,0,0)*rat(ep^2*n8*n10+ep*n4*n8*n10+ep*n8^2*n10-ep*n8*n10,2*ep^3+4*ep^2*
      (n2-1)+ep^2*n10-4*ep^2+2*ep*(n2-1)^2+2*ep*(n2-1)*n10-4*ep*(n2-1)-2*ep*n10+2*ep+(n2-1)^2*n10-2*(n2-1)*n10+n10)
      +`Y'(1,(n2-1),0,n4,n5,1,1,n8,2,1+n10,0,0,0,0)*rat(-n10,1)
      +`Y'(1,(n2-1),2,n4,n5,1,1,n8,0,1+n10,0,0,0,0)*rat(-ep*n10,ep+(n2-1)-1)
      );
*--#] n2 : 
#endif
#if ( `il' >= 12 )
*--#[ n5 : 
*`SKIP48'
id,ifmatch->nono1,
	`Z'(1,1,1,n4?`POS4',n5?{>1},1,1,n8?`POS8',1,0,0,0,0,0)
			 = -rat(1,-2*ep*(n5-1)-2*(n5-1)^2-2*(n5-1)*n8+2*(n5-1))*(
      +`Z'(1,1,1,-1+n4,1+(n5-1),1,1,1+n8,1,0,0,0,0,0)*rat(2*(n5-1)*n8,1)
      +`Z'(1,1,1,n4,(n5-1),1,1,1+n8,1,0,0,0,0,0)*rat(-2*ep*n8-2*(n5-1)*n8-2*n8^2+2*n8,1)
      +`Z'(1,1,1,-2+n4,1+(n5-1),1,1,1+n8,1,0,0,0,0,0)*rat(-2*(n5-1)*n8,1)
      +`Z'(1,1,1,-1+n4,1+(n5-1),1,1,n8,1,0,0,0,0,0)*rat(10*ep*(n5-1)+4*(n5-1)^2+4*(n5-1)*n8-6*(n5-1),1)
      +`Z'(1,1,1,-1+n4,(n5-1),1,1,1+n8,1,0,0,0,0,0)*rat(10*ep*n8+4*(n5-1)*n8+4*n8^2-6*n8,1)
      +`Z'(1,1,1,n4,-1+(n5-1),1,1,1+n8,1,0,0,0,0,0)*rat(-10*ep*n8-2*(n5-1)*n8-4*n8^2+6*n8,1)
      +`Z'(1,1,1,n4,1+(n5-1),1,1,-1+n8,1,0,0,0,0,0)*rat(-10*ep*(n5-1)-4*(n5-1)^2-2*(n5-1)*n8+6*(n5-1),1)
      +`Z'(1,1,1,n4,(n5-1),1,1,n8,1,0,0,0,0,0)*rat(-48*ep^2-30*ep*(n5-1)-30*ep*n8+68*ep-4*(n5-1)^2-12*(n5-1)*n8+22*(n5-1)-4*n8^2+22*n8-24,1)
      +`Y'(0,1,1,-1+n4,(n5-1),2,1,1+n8,1,0,0,0,0,0)*rat(n8,1)
      +`Y'(0,1,1,n4,-1+(n5-1),2,1,1+n8,1,0,0,0,0,0)*rat(-n8,1)
      +`Y'(0,1,1,n4,(n5-1),2,1,n8,1,0,0,0,0,0)*rat(-4*ep-2*(n5-1)-n8+3,1)
      +`Y'(0,2,0,n4,(n5-1),2,1,n8,1,0,0,0,0,0)*rat(1,1)
      +`Y'(0,2,1,n4,(n5-1),2,1,n8,0,0,0,0,0,0)*rat(-1,1)
      +`Y'(1,0,1,-1+n4,(n5-1),1,1,1+n8,2,0,0,0,0,0)*rat(2*n8,1)
      +`Y'(1,0,1,n4,-1+(n5-1),1,1,1+n8,2,0,0,0,0,0)*rat(-2*n8,1)
      +`Y'(1,0,1,n4,(n5-1),1,1,n8,2,0,0,0,0,0)*rat(-10*ep-4*(n5-1)-2*n8+6,1)
      +`Y'(1,1,0,-1+n4,(n5-1),1,1,1+n8,2,0,0,0,0,0)*rat(-2*n8,1)
      +`Y'(1,1,0,n4,-1+(n5-1),1,1,1+n8,2,0,0,0,0,0)*rat(2*n8,1)
      +`Y'(1,1,0,n4,(n5-1),1,1,n8,2,0,0,0,0,0)*rat(10*ep+4*(n5-1)+2*n8-4,1)
      +`Y'(1,2,-1,n4,(n5-1),1,1,n8,2,0,0,0,0,0)*rat(-2,1)
      +`Y'(1,2,0,-1+n4,1+(n5-1),1,1,n8,1,0,0,0,0,0)*rat(-2*(n5-1),1)
      +`Y'(1,2,0,n4,1+(n5-1),1,1,-1+n8,1,0,0,0,0,0)*rat(2*(n5-1),1)
      +`Y'(1,2,0,n4,(n5-1),1,1,n8,1,0,0,0,0,0)*rat(10*ep+2*(n5-1)+4*n8-4,1)
      +`Y'(1,2,1,-1+n4,1+(n5-1),1,1,n8,0,0,0,0,0,0)*rat(2*(n5-1),1)
      +`Y'(1,2,1,n4,1+(n5-1),1,1,-1+n8,0,0,0,0,0,0)*rat(-2*(n5-1),1)
      +`Y'(1,2,1,n4,(n5-1),1,1,n8,0,0,0,0,0,0)*rat(-10*ep-2*(n5-1)-4*n8+6,1)
      +`Y'(2,0,1,n4,(n5-1),0,1,n8,2,0,0,0,0,0)*rat(-1,1)
      +`Y'(2,1,0,n4,(n5-1),0,1,n8,2,0,0,0,0,0)*rat(1,1)
      +`Y'(2,1,1,-1+n4,1+(n5-1),0,1,n8,1,0,0,0,0,0)*rat((n5-1),1)
      +`Y'(2,1,1,n4,1+(n5-1),0,1,-1+n8,1,0,0,0,0,0)*rat(-(n5-1),1)
      +`Y'(2,1,1,n4,(n5-1),0,1,n8,1,0,0,0,0,0)*rat(-4*ep-(n5-1)-2*n8+3,1)
		);
*`SKIP48'
id,ifmatch->nono1,
	`Z'(1,1,1,n4?`POS4',n5?{>1},1,1,n8?`POS8',1,n10?neg_,0,0,0,0)
			 = -rat(2*ep+n10,-2*ep-2*(n5-1)-2*n8+2)*(
      +`Z'(1,1,1,-1+n4,1+(n5-1),1,1,1+n8,1,n10,0,0,0,0)*rat(2*n8,2*ep+n10)
      +`Z'(1,1,1,n4,(n5-1),1,1,1+n8,1,n10,0,0,0,0)*rat(-2*ep*n8-2*(n5-1)*n8-2*n8^2+2*n8,2*ep*(n5-1)+(n5-1)*n10)
      +`Z'(1,1,1,-1+n4,1+(n5-1),1,1,n8,1,1+n10,0,0,0,0)*rat(2*ep,2*ep+n10)
      +`Z'(1,1,1,n4,1+(n5-1),1,1,-1+n8,1,1+n10,0,0,0,0)*rat(-2*ep,2*ep+n10)
      +`Z'(1,1,1,n4,(n5-1),1,1,n8,1,1+n10,0,0,0,0)*rat(-8*ep^2-2*ep*(n5-1)-4*ep*n8-8*ep*n10+4*ep-2*n10^2,2*ep*(n5-1)+(n5-1)*n10)
      +`Z'(1,1,1,n4,(n5-1),1,2,-1+n8,1,1+n10,0,0,0,0)*rat(-2*ep,2*ep*(n5-1)+(n5-1)*n10)
      +`Y'(0,1,1,n4,(n5-1),1,2,n8,1,1+n10,0,0,0,0)*rat(2*ep,2*ep*(n5-1)+(n5-1)*n10)
      +`Y'(1,0,1,n4,(n5-1),1,1,n8,2,1+n10,0,0,0,0)*rat(-2*ep-2*n10,2*ep*(n5-1)+(n5-1)*n10)
      +`Y'(1,0,1,n4,(n5-1),1,2,n8,1,1+n10,0,0,0,0)*rat(-2*ep,2*ep*(n5-1)+(n5-1)*n10)
      +`Y'(1,1,0,n4,(n5-1),1,1,n8,2,1+n10,0,0,0,0)*rat(2*ep+2*n10,2*ep*(n5-1)+(n5-1)*n10)
      +`Y'(1,0,1,n4,(n5-1),1,1,n8,1,2+n10,0,0,0,0)*rat(2*n10^2+2*n10,2*ep*(n5-1)+(n5-1)*n10)
		);
*--#] n5 : 
#endif
#if ( `il' >= 13 )
*--#[ n10 :
*
*	This works for n8 > 1
*
*#if(1==1)
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,1,1,n4?`POS4',1,1,1,n8?{>1},1,n10?neg_,0,0,0,0)
			 = -rat(1,24*ep^2*n8-24*ep^2+14*ep*n4*n8-14*ep*n4+12*ep*n8^2+6*ep*n8*(1+n10)-42*
      ep*n8-6*ep*(1+n10)+30*ep+2*n4^2*n8-2*n4^2+4*n4*n8^2+2*n4*n8*(1+n10)-14*n4*n8-2*n4*(1+n10)+10*n4+2*n8^3+n8^2*(1+n10)-11*n8^2-3*n8*
      (1+n10)+19*n8+2*(1+n10)-10)*(
      +`Z'(1,1,1,1+n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(12*ep^2*n4-8*ep^2*n8+16*ep^2+7*ep*n4^2+4*ep*n4*n8+6*ep*n4*(1+n10)-11*
      ep*n4-4*ep*n8^2-4*ep*n8*(1+n10)+18*ep*n8+8*ep*(1+n10)-20*ep+n4^3+2*n4^2*n8+2*n4^2*(1+n10)-5*n4^2,1)
      +`Z'(1,1,1,2+n4,1,1,1,-2+n8,1,(1+n10),0,0,0,0)*rat(-16*ep^2*n4-16*ep^2-7*ep*n4^2-6*ep*n4*n8-4*ep*n4*(1+n10)+9*ep*n4-6*ep*
      n8-4*ep*(1+n10)+16*ep-n4^3-2*n4^2*n8-n4^2*(1+n10)+4*n4^2-2*n4*n8-n4*(1+n10)+5*n4,1)
      +`Z'(1,1,1,1+n4,1,1,1,-2+n8,1,1+(1+n10),-1,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
      +`Z'(1,1,1,1+n4,1,1,1,-2+n8,1,1+(1+n10),0,-1,0,0)*rat(16*ep^2*(1+n10)+4*ep*n4*(1+n10)+8*ep*n8*(1+n10)+4*ep*(1+n10)^2-20*ep*(1+n10)+2*n4*
      n8*(1+n10)+n4*(1+n10)^2-5*n4*(1+n10),1)
      +`Z'(1,1,1,1+n4,1,1,1,-2+n8,1,(1+n10),0,0,0,0)*rat(32*ep^3+24*ep^2*n4+16*ep^2*n8+8*ep^2*(1+n10)-40*ep^2+7*ep*n4^2+7*ep*n4
      *n8+6*ep*n4*(1+n10)-17*ep*n4+2*ep*n8^2-10*ep*n8+12*ep+n4^3+n4^2*n8+n4^2*(1+n10)-2*n4^2,1)
      +`Z'(1,1,1,n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(-48*ep^3-40*ep^2*n4-16*ep^2*n8-36*ep^2*(1+n10)+68*ep^2-11*ep*n4^2-9*ep*
      n4*n8-21*ep*n4*(1+n10)+39*ep*n4-2*ep*n8^2-12*ep*n8*(1+n10)+20*ep*n8-6*ep*(1+n10)^2+36*ep*(1+n10)-38*ep-n4^3-n4^2*n8-3*n4^2*(1+n10)+5*
      n4^2-4*n4*n8*(1+n10)+4*n4*n8-2*n4*(1+n10)^2+12*n4*(1+n10)-10*n4-2*n8^2*(1+n10)+2*n8^2-n8*(1+n10)^2+10*n8*(1+n10)-9*n8+2*(1+n10)^2-12*(1+n10)+10,1
      )
      +`Z'(1,1,1,1+n4,1,1,1,-2+n8,1,1+(1+n10),0,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
      +`Z'(1,1,1,n4,1,1,1,-2+n8,1,1+(1+n10),0,0,0,0)*rat(-3*ep*n4*(1+n10)+2*ep*n8*(1+n10)-4*ep*(1+n10)-n4^2*(1+n10),1)
      +`Y'(1,1,1,1+n4,0,1,1,n8,1,-1+(1+n10),0,0,0,0)*rat(-32*ep^2*n8+32*ep^2-18*ep*n4*n8+18*ep*n4-20*ep*n8^2-8*ep*n8*(1+n10)+
      68*ep*n8+8*ep*(1+n10)-48*ep-2*n4^2*n8+2*n4^2-6*n4*n8^2-3*n4*n8*(1+n10)+21*n4*n8+3*n4*(1+n10)-15*n4-4*n8^3-2*n8^2*(1+n10)+22*n8^2+
      6*n8*(1+n10)-38*n8-4*(1+n10)+20,2)
      +`Y'(1,2,0,1+n4,1,1,1,-1+n8,1,-1+(1+n10),0,0,0,0)*rat(48*ep^2+32*ep*n4+32*ep*n8+12*ep*(1+n10)-76*ep+4*n4^2+10*n4*n8+5*n4
      *(1+n10)-25*n4+8*n8^2+4*n8*(1+n10)-36*n8-8*(1+n10)+40,2)
      +`Y'(1,2,1,1+n4,1,1,1,-1+n8,0,-1+(1+n10),0,0,0,0)*rat(-32*ep^2-18*ep*n4-20*ep*n8-8*ep*(1+n10)+48*ep-2*n4^2-6*n4*n8-3*n4*
      (1+n10)+15*n4-4*n8^2-2*n8*(1+n10)+18*n8+4*(1+n10)-20,2)
      +`Y'(1,0,2,1+n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(-16*ep^2-11*ep*n4-14*ep*n8-4*ep*(1+n10)+32*ep-n4^2-4*n4*n8-2*n4*(1+n10)+10
      *n4-4*n8^2-2*n8*(1+n10)+18*n8+4*(1+n10)-20,1)
      +`Y'(1,1,1,2+n4,0,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(16*ep^2*n4+16*ep^2+4*ep*n4^2+8*ep*n4*n8+4*ep*n4*(1+n10)-16*ep*n4+8*ep*
      n8+4*ep*(1+n10)-20*ep+2*n4^2*n8+n4^2*(1+n10)-5*n4^2+2*n4*n8+n4*(1+n10)-5*n4,2)
      +`Y'(1,1,2,1+n4,1,1,1,-1+n8,0,(1+n10),0,0,0,0)*rat(8*ep^2+4*ep*n4+8*ep*n8+2*ep*(1+n10)-18*ep+2*n4*n8+n4*(1+n10)-5*n4+2*n8^2+
      n8*(1+n10)-9*n8-2*(1+n10)+10,1)
      +`Y'(0,1,1,1+n4,0,1,1,n8,1,(1+n10),0,0,0,0)*rat(3*ep*n4*n8-3*ep*n4-2*ep*n8^2+6*ep*n8-4*ep+n4^2*n8-n4^2,1)
      +`Y'(0,1,1,1+n4,1,1,1,-1+n8,1,1+(1+n10),-1,0,0,0)*rat(-3*ep*n4*(1+n10)+2*ep*n8*(1+n10)-4*ep*(1+n10)-n4^2*(1+n10),1)
      +`Y'(0,1,1,1+n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(12*ep^2*n4-8*ep^2*n8+16*ep^2+4*ep*n4^2+3*ep*n4*n8+3*ep*n4*(1+n10)-9*ep
      *n4-2*ep*n8^2-2*ep*n8*(1+n10)+10*ep*n8+4*ep*(1+n10)-12*ep+n4^2*n8+n4^2*(1+n10)-3*n4^2,1)
      +`Y'(0,1,1,n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(-3*ep*n4*n8+3*ep*n4+2*ep*n8^2-6*ep*n8+4*ep-n4^2*n8+n4^2,1)
      +`Y'(0,2,0,1+n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(-3*ep*n4+2*ep*n8-4*ep-n4^2,1)
      +`Y'(0,2,1,1+n4,1,1,1,-1+n8,0,(1+n10),0,0,0,0)*rat(3*ep*n4-2*ep*n8+4*ep+n4^2,1)
      +`Y'(1,0,1,1+n4,0,1,1,-1+n8,2,(1+n10),0,0,0,0)*rat(8*ep^2+7*ep*n4+6*ep*n8+2*ep*(1+n10)-14*ep+n4^2+2*n4*n8+n4*(1+n10)-5*n4+2*
      n8^2+n8*(1+n10)-9*n8-2*(1+n10)+10,1)
      +`Y'(1,0,1,1+n4,0,1,1,n8,1,(1+n10),0,0,0,0)*rat(-3*ep*n4*n8+3*ep*n4+2*ep*n8^2-6*ep*n8+4*ep-n4^2*n8+n4^2,1)
      +`Y'(1,0,1,1+n4,1,0,1,-1+n8,2,(1+n10),0,0,0,0)*rat(-3*ep*n4+2*ep*n8-4*ep-n4^2,1)
      +`Y'(1,0,1,1+n4,1,1,0,-1+n8,2,(1+n10),0,0,0,0)*rat(-3*ep*n4+2*ep*n8-4*ep-n4^2,1)
      +`Y'(1,0,1,1+n4,1,1,1,-1+n8,1,1+(1+n10),-1,0,0,0)*rat(16*ep^2*(1+n10)+11*ep*n4*(1+n10)+14*ep*n8*(1+n10)+4*ep*(1+n10)^2-32*ep*(1+n10)+n4^2
      *(1+n10)+4*n4*n8*(1+n10)+2*n4*(1+n10)^2-10*n4*(1+n10)+4*n8^2*(1+n10)+2*n8*(1+n10)^2-18*n8*(1+n10)-4*(1+n10)^2+20*(1+n10),1)
      +`Y'(1,0,1,1+n4,1,1,1,-1+n8,1,1+(1+n10),0,-1,0,0)*rat(-16*ep^2*(1+n10)-8*ep*n4*(1+n10)-16*ep*n8*(1+n10)-4*ep*(1+n10)^2+36*ep*(1+n10)-4*n4
      *n8*(1+n10)-2*n4*(1+n10)^2+10*n4*(1+n10)-4*n8^2*(1+n10)-2*n8*(1+n10)^2+18*n8*(1+n10)+4*(1+n10)^2-20*(1+n10),1)
      +`Y'(1,0,1,1+n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(32*ep^3+32*ep^2*n4+48*ep^2*n8+24*ep^2*(1+n10)-104*ep^2+14*ep*n4^2+31*
      ep*n4*n8+16*ep*n4*(1+n10)-67*ep*n4+22*ep*n8^2+24*ep*n8*(1+n10)-98*ep*n8+4*ep*(1+n10)^2-52*ep*(1+n10)+108*ep+2*n4^3+5*n4^2*n8+2*
      n4^2*(1+n10)-11*n4^2+8*n4*n8^2+8*n4*n8*(1+n10)-36*n4*n8+2*n4*(1+n10)^2-18*n4*(1+n10)+40*n4+4*n8^3+6*n8^2*(1+n10)-26*n8^2+2*n8*(1+n10)^2-
      26*n8*(1+n10)+56*n8-4*(1+n10)^2+28*(1+n10)-40,1)
      +`Y'(1,0,1,n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(3*ep*n4*n8-3*ep*n4-2*ep*n8^2+6*ep*n8-4*ep+n4^2*n8-n4^2,1)
      +`Y'(1,0,2,n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(8*ep^2+7*ep*n4+6*ep*n8+2*ep*(1+n10)-14*ep+n4^2+2*n4*n8+n4*(1+n10)-5*n4+2*
      n8^2+n8*(1+n10)-9*n8-2*(1+n10)+10,1)
      +`Y'(1,1,0,1+n4,0,1,1,-1+n8,2,(1+n10),0,0,0,0)*rat(-8*ep^2-7*ep*n4-6*ep*n8-2*ep*(1+n10)+14*ep-n4^2-2*n4*n8-n4*(1+n10)+5*n4-2*
      n8^2-n8*(1+n10)+9*n8+2*(1+n10)-10,1)
      +`Y'(1,1,0,1+n4,1,0,1,-1+n8,2,(1+n10),0,0,0,0)*rat(3*ep*n4-2*ep*n8+4*ep+n4^2,1)
      +`Y'(1,1,0,1+n4,1,1,0,-1+n8,2,(1+n10),0,0,0,0)*rat(3*ep*n4-2*ep*n8+4*ep+n4^2,1)
      +`Y'(1,1,0,1+n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(-96*ep^3-136*ep^2*n4-112*ep^2*n8-72*ep^2*(1+n10)+296*ep^2-54*ep*n4^2-
      96*ep*n4*n8-60*ep*n4*(1+n10)+252*ep*n4-48*ep*n8^2-52*ep*n8*(1+n10)+244*ep*n8-12*ep*(1+n10)^2+128*ep*(1+n10)-308*ep-6*n4^3-18*n4^2
      *n8-11*n4^2*(1+n10)+47*n4^2-20*n4*n8^2-20*n4*n8*(1+n10)+100*n4*n8-5*n4*(1+n10)^2+50*n4*(1+n10)-125*n4-8*n8^3-12*n8^2*(1+n10)+60*n8^2-
      4*n8*(1+n10)^2+56*n8*(1+n10)-148*n8+8*(1+n10)^2-64*(1+n10)+120,2)
      +`Y'(1,1,0,2+n4,1,1,1,-2+n8,1,(1+n10),0,0,0,0)*rat(16*ep^2*n4+16*ep^2+10*ep*n4^2+4*ep*n4*n8+4*ep*n4*(1+n10)-2*ep*n4+4*ep*
      n8+4*ep*(1+n10)-12*ep+2*n4^3+2*n4^2*n8+n4^2*(1+n10)-3*n4^2+2*n4*n8+n4*(1+n10)-5*n4,2)
      +`Y'(1,1,1,1+n4,0,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(32*ep^3+24*ep^2*n4+16*ep^2*n8+40*ep^2*(1+n10)-72*ep^2+4*ep*n4^2+6*ep*
      n4*n8+24*ep*n4*(1+n10)-36*ep*n4+4*ep*n8^2+20*ep*n8*(1+n10)-36*ep*n8+8*ep*(1+n10)^2-56*ep*(1+n10)+64*ep+3*n4^2*(1+n10)-3*n4^2+6*n4*n8*
      (1+n10)-6*n4*n8+3*n4*(1+n10)^2-18*n4*(1+n10)+15*n4+4*n8^2*(1+n10)-4*n8^2+2*n8*(1+n10)^2-20*n8*(1+n10)+18*n8-4*(1+n10)^2+24*(1+n10)-20,2)
      +`Y'(1,1,1,1+n4,1,0,1,-1+n8,1,(1+n10),0,0,0,0)*rat(-12*ep^2*n4+8*ep^2*n8-16*ep^2-4*ep*n4^2-6*ep*n4*n8-6*ep*n4*(1+n10)+15*
      ep*n4+4*ep*n8^2+4*ep*n8*(1+n10)-18*ep*n8-8*ep*(1+n10)+20*ep-2*n4^2*n8-2*n4^2*(1+n10)+5*n4^2,1)
      +`Y'(1,1,1,1+n4,1,1,0,-1+n8,1,(1+n10),0,0,0,0)*rat(-12*ep^2*n4+8*ep^2*n8-16*ep^2-4*ep*n4^2-6*ep*n4*n8-6*ep*n4*(1+n10)+18*
      ep*n4+4*ep*n8^2+4*ep*n8*(1+n10)-20*ep*n8-8*ep*(1+n10)+24*ep-2*n4^2*n8-2*n4^2*(1+n10)+6*n4^2,1)
      +`Y'(1,1,1,1+n4,1,1,1,-1+n8,0,(1+n10),0,0,0,0)*rat(32*ep^3+32*ep^2*n4+32*ep^2*n8+40*ep^2*(1+n10)-104*ep^2+8*ep*n4^2+24*ep
      *n4*n8+26*ep*n4*(1+n10)-74*ep*n4+8*ep*n8^2+24*ep*n8*(1+n10)-56*ep*n8+8*ep*(1+n10)^2-64*ep*(1+n10)+88*ep+4*n4^2*n8+4*n4^2*(1+n10)-12*
      n4^2+4*n4*n8^2+8*n4*n8*(1+n10)-24*n4*n8+3*n4*(1+n10)^2-22*n4*(1+n10)+35*n4+4*n8^2*(1+n10)-4*n8^2+2*n8*(1+n10)^2-20*n8*(1+n10)+18*n8-4*
      (1+n10)^2+24*(1+n10)-20,2)
      +`Y'(1,1,1,1+n4,2,0,1,-2+n8,1,(1+n10),0,0,0,0)*rat(-3*ep*n4+2*ep*n8-4*ep-n4^2,1)
      +`Y'(1,1,1,1+n4,2,1,0,-2+n8,1,(1+n10),0,0,0,0)*rat(-3*ep*n4+2*ep*n8-4*ep-n4^2,1)
      +`Y'(1,1,1,1+n4,2,1,1,-2+n8,0,(1+n10),0,0,0,0)*rat(16*ep^2+10*ep*n4+4*ep*n8+4*ep*(1+n10)-12*ep+2*n4^2+2*n4*n8+n4*(1+n10)-5*n4
      ,2)
      +`Y'(1,1,1,n4,2,0,1,-1+n8,1,(1+n10),0,0,0,0)*rat(3*ep*n4-2*ep*n8+4*ep+n4^2,1)
      +`Y'(1,1,1,n4,2,1,0,-1+n8,1,(1+n10),0,0,0,0)*rat(3*ep*n4-2*ep*n8+4*ep+n4^2,1)
      +`Y'(1,1,1,n4,2,1,1,-1+n8,0,(1+n10),0,0,0,0)*rat(-16*ep^2-10*ep*n4-4*ep*n8-4*ep*(1+n10)+12*ep-2*n4^2-2*n4*n8-n4*(1+n10)+5*n4,
      2)
      +`Y'(1,2,0,1+n4,1,1,1,-2+n8,1,(1+n10),0,0,0,0)*rat(3*ep*n4-2*ep*n8+4*ep+n4^2,1)
      +`Y'(1,2,1,1+n4,1,1,1,-2+n8,0,(1+n10),0,0,0,0)*rat(-3*ep*n4+2*ep*n8-4*ep-n4^2,1)
      +`Y'(0,1,1,1+n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-3*ep*n4*(1+n10)+2*ep*n8*(1+n10)-4*ep*(1+n10)-n4^2*(1+n10),1)
      +`Y'(1,0,1,1+n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(16*ep^2*(1+n10)+11*ep*n4*(1+n10)+14*ep*n8*(1+n10)+4*ep*(1+n10)^2-32*ep*(1+n10)+n4^2*
      (1+n10)+4*n4*n8*(1+n10)+2*n4*(1+n10)^2-10*n4*(1+n10)+4*n8^2*(1+n10)+2*n8*(1+n10)^2-18*n8*(1+n10)-4*(1+n10)^2+20*(1+n10),1)
      +`Y'(0,1,0,1+n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
      +`Y'(0,1,1,1+n4,0,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-3*ep*n4*(1+n10)+2*ep*n8*(1+n10)-4*ep*(1+n10)-n4^2*(1+n10),1)
      +`Y'(0,1,1,1+n4,1,0,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
      +`Y'(0,1,1,1+n4,1,1,0,-1+n8,1,1+(1+n10),0,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
      +`Y'(0,1,1,1+n4,1,1,1,-1+n8,0,1+(1+n10),0,0,0,0)*rat(-3*ep*n4*(1+n10)+2*ep*n8*(1+n10)-4*ep*(1+n10)-n4^2*(1+n10),1)
      +`Y'(0,1,1,n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
      +`Y'(1,-1,1,1+n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(16*ep^2*(1+n10)+8*ep*n4*(1+n10)+16*ep*n8*(1+n10)+4*ep*(1+n10)^2-36*ep*(1+n10)+4*n4*
      n8*(1+n10)+2*n4*(1+n10)^2-10*n4*(1+n10)+4*n8^2*(1+n10)+2*n8*(1+n10)^2-18*n8*(1+n10)-4*(1+n10)^2+20*(1+n10),1)
      +`Y'(1,0,0,1+n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-3*ep*n4*(1+n10)+2*ep*n8*(1+n10)-4*ep*(1+n10)-n4^2*(1+n10),1)
      +`Y'(1,0,1,1+n4,0,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
      +`Y'(1,0,1,1+n4,1,0,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-16*ep^2*(1+n10)-11*ep*n4*(1+n10)-14*ep*n8*(1+n10)-4*ep*(1+n10)^2+32*ep*(1+n10)-n4^2
      *(1+n10)-4*n4*n8*(1+n10)-2*n4*(1+n10)^2+10*n4*(1+n10)-4*n8^2*(1+n10)-2*n8*(1+n10)^2+18*n8*(1+n10)+4*(1+n10)^2-20*(1+n10),1)
      +`Y'(1,0,1,1+n4,1,1,0,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-16*ep^2*(1+n10)-11*ep*n4*(1+n10)-14*ep*n8*(1+n10)-4*ep*(1+n10)^2+32*ep*(1+n10)-n4^2
      *(1+n10)-4*n4*n8*(1+n10)-2*n4*(1+n10)^2+10*n4*(1+n10)-4*n8^2*(1+n10)-2*n8*(1+n10)^2+18*n8*(1+n10)+4*(1+n10)^2-20*(1+n10),1)
      +`Y'(1,0,1,1+n4,1,1,1,-2+n8,1,1+(1+n10),0,0,0,0)*rat(-48*ep^2*(1+n10)-28*ep*n4*(1+n10)-24*ep*n8*(1+n10)-12*ep*(1+n10)^2+60*ep*(1+n10)-4*
      n4^2*(1+n10)-8*n4*n8*(1+n10)-4*n4*(1+n10)^2+20*n4*(1+n10)-4*n8^2*(1+n10)-2*n8*(1+n10)^2+18*n8*(1+n10)+4*(1+n10)^2-20*(1+n10),1)
      +`Y'(1,0,1,1+n4,1,1,1,-1+n8,0,1+(1+n10),0,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
      +`Y'(1,0,1,n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(16*ep^2*(1+n10)+7*ep*n4*(1+n10)+6*ep*n8*(1+n10)+4*ep*(1+n10)^2-16*ep*(1+n10)+n4^2*(1+n10)+
      2*n4*n8*(1+n10)+n4*(1+n10)^2-5*n4*(1+n10),1)
      +`Y'(1,1,0,1+n4,1,1,1,-2+n8,1,1+(1+n10),0,0,0,0)*rat(16*ep^2*(1+n10)+11*ep*n4*(1+n10)+14*ep*n8*(1+n10)+4*ep*(1+n10)^2-32*ep*(1+n10)+n4^2*
      (1+n10)+4*n4*n8*(1+n10)+2*n4*(1+n10)^2-10*n4*(1+n10)+4*n8^2*(1+n10)+2*n8*(1+n10)^2-18*n8*(1+n10)-4*(1+n10)^2+20*(1+n10),1)
      +`Y'(1,1,1,1+n4,0,1,1,-2+n8,1,1+(1+n10),0,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
      +`Y'(1,1,1,1+n4,1,0,1,-2+n8,1,1+(1+n10),0,0,0,0)*rat(-3*ep*n4*(1+n10)+2*ep*n8*(1+n10)-4*ep*(1+n10)-n4^2*(1+n10),1)
      +`Y'(1,1,1,1+n4,1,1,0,-2+n8,1,1+(1+n10),0,0,0,0)*rat(-3*ep*n4*(1+n10)+2*ep*n8*(1+n10)-4*ep*(1+n10)-n4^2*(1+n10),1)
      +`Y'(1,1,1,1+n4,1,1,1,-2+n8,0,1+(1+n10),0,0,0,0)*rat(3*ep*n4*(1+n10)-2*ep*n8*(1+n10)+4*ep*(1+n10)+n4^2*(1+n10),1)
		);
*#endif
*#if(1==1)
*`SKIP4'
id,ifmatch->nono1,
	`Z'(1,1,1,n4?`POS4',1,1,1,n8?`POS8',1,n10?neg_,0,0,0,0)
		 = -rat(1,-48*ep^2-32*ep*n8-28*ep*(1+n10)+68*ep-4*n8^2-10*n8*(1+n10)+22*n8-4*(1+n10)^2+20*(1+n10)-24)*(
      +`Z'(1,1,1,-1+n4,1,1,1,1+n8,1,-1+(1+n10),0,0,0,0)*rat(12*ep*n8+4*n8^2+6*n8*(1+n10)-10*n8,1)
      +`Z'(1,1,1,1+n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(-2*ep*n4-2*n4*(1+n10)+2*n4,1)
      +`Z'(1,1,1,n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(8*ep^2+2*ep*n4+4*ep*n8+12*ep*(1+n10)-16*ep+2*n4*(1+n10)-2*n4+4*n8*(1+n10)-4*n8+4*
      (1+n10)^2-12*(1+n10)+8,1)
      +`Z'(2,1,1,n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(-4*ep-2*n8-(1+n10)+3,1)
      +`Z'(1,1,1,-1+n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(-8*ep^2-2*ep*n4+6*ep*n8-20*ep*(1+n10)+16*ep-2*n4*(1+n10)+2*n4+2*n8^2-n8*(1+n10)-
      3*n8-6*(1+n10)^2+15*(1+n10)-7,1)
      +`Z'(1,1,1,n4,1,1,1,-1+n8,1,1+(1+n10),-1,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Z'(1,1,1,n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(-16*ep^2+2*ep*n4-14*ep*n8-8*ep*(1+n10)+26*ep+2*n4*(1+n10)-2*n4-2*n8^2-5*n8*
      (1+n10)+11*n8-(1+n10)^2+8*(1+n10)-11,1)
      +`Z'(1,1,1,n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Z'(1,1,1,-1+n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(1,1,1,n4,0,1,1,1+n8,1,-1+(1+n10),0,0,0,0)*rat(-10*ep*n8-4*n8^2-4*n8*(1+n10)+8*n8,1)
      +`Y'(1,2,0,n4,1,1,1,n8,1,-1+(1+n10),0,0,0,0)*rat(12*ep+4*n8+6*(1+n10)-10,1)
      +`Y'(1,2,1,n4,1,1,1,n8,0,-1+(1+n10),0,0,0,0)*rat(-10*ep-4*n8-4*(1+n10)+8,1)
      +`Y'(1,0,2,n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(-2*ep-2*(1+n10)+2,1)
      +`Y'(1,1,1,-1+n4,1,0,1,1+n8,1,-1+(1+n10),0,0,0,0)*rat(4*ep*n8+2*n8^2+n8*(1+n10)-3*n8,1)
      +`Y'(1,1,1,n4,0,0,1,1+n8,1,-1+(1+n10),0,0,0,0)*rat(-4*ep*n8-2*n8^2-n8*(1+n10)+3*n8,1)
      +`Y'(1,1,1,n4,1,0,1,n8,1,-1+(1+n10),0,0,0,0)*rat(-24*ep^2-16*ep*n8-10*ep*(1+n10)+38*ep-2*n8^2-3*n8*(1+n10)+13*n8-(1+n10)^2+8*(1+n10)-15,1)
      +`Y'(1,1,1,n4,1,0,1,n8,1,(1+n10),-1,0,0,0)*rat(4*ep*(1+n10)-4*ep+2*n8*(1+n10)-2*n8+(1+n10)^2-4*(1+n10)+3,1)
      +`Y'(1,2,0,n4,1,0,1,n8,1,-1+(1+n10),0,0,0,0)*rat(4*ep+2*n8+(1+n10)-3,1)
      +`Y'(1,2,1,n4,1,0,1,n8,0,-1+(1+n10),0,0,0,0)*rat(-4*ep-2*n8-(1+n10)+3,1)
      +`Y'(2,0,1,n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(-4*ep-2*n8-(1+n10)+3,1)
      +`Y'(2,1,1,n4,1,-1,1,n8,1,-1+(1+n10),0,0,0,0)*rat(-4*ep-2*n8-(1+n10)+3,1)
      +`Y'(2,1,1,n4,1,1,0,n8,1,(1+n10),0,0,0,0)*rat(4*ep+2*n8+(1+n10)-3,1)
      +`Y'(0,1,1,-1+n4,1,1,1,1+n8,1,(1+n10),0,0,0,0)*rat(-6*ep*n8-2*n8^2-3*n8*(1+n10)+5*n8,1)
      +`Y'(0,1,1,n4,0,1,1,1+n8,1,(1+n10),0,0,0,0)*rat(6*ep*n8+2*n8^2+3*n8*(1+n10)-5*n8,1)
      +`Y'(0,1,1,n4,1,1,1,n8,1,1+(1+n10),-1,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(0,1,1,n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(16*ep^2+10*ep*n8+16*ep*(1+n10)-22*ep+2*n8^2+5*n8*(1+n10)-7*n8+3*(1+n10)^2-10*(1+n10)+7,1)
      +`Y'(0,2,0,n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(-6*ep-2*n8-3*(1+n10)+5,1)
      +`Y'(0,2,1,n4,1,1,1,n8,0,(1+n10),0,0,0,0)*rat(6*ep+2*n8+3*(1+n10)-5,1)
      +`Y'(1,0,1,-1+n4,1,1,1,1+n8,1,(1+n10),0,0,0,0)*rat(6*ep*n8+2*n8^2+3*n8*(1+n10)-5*n8,1)
      +`Y'(1,0,1,n4,0,1,1,1+n8,1,(1+n10),0,0,0,0)*rat(-6*ep*n8-2*n8^2-3*n8*(1+n10)+5*n8,1)
      +`Y'(1,0,1,n4,0,1,1,n8,2,(1+n10),0,0,0,0)*rat(2*ep+2*(1+n10)-2,1)
      +`Y'(1,0,1,n4,1,0,1,n8,2,(1+n10),0,0,0,0)*rat(-2*ep-2*(1+n10)+2,1)
      +`Y'(1,0,1,n4,1,1,0,n8,2,(1+n10),0,0,0,0)*rat(-2*ep-2*(1+n10)+2,1)
      +`Y'(1,0,1,n4,1,1,1,n8,1,1+(1+n10),-1,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Y'(1,0,1,n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(-16*ep^2+4*ep*n4-10*ep*n8-8*ep*(1+n10)+20*ep+4*n4*(1+n10)-4*n4-2*n8^2-n8*(1+n10)+7*
      n8-(1+n10)^2+2*(1+n10)-5,1)
      +`Y'(1,0,2,-1+n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(2*ep+2*(1+n10)-2,1)
      +`Y'(1,1,0,1+n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(2*ep*n4+2*n4*(1+n10)-2*n4,1)
      +`Y'(1,1,0,n4,0,1,1,n8,2,(1+n10),0,0,0,0)*rat(-2*ep-2*(1+n10)+2,1)
      +`Y'(1,1,0,n4,1,0,1,n8,2,(1+n10),0,0,0,0)*rat(2*ep+2*(1+n10)-2,1)
      +`Y'(1,1,0,n4,1,1,0,n8,2,(1+n10),0,0,0,0)*rat(2*ep+2*(1+n10)-2,1)
      +`Y'(1,1,0,n4,1,1,1,n8,1,(1+n10),0,0,0,0)*rat(-8*ep^2-6*ep*n4-4*ep*n8-20*ep*(1+n10)+30*ep-6*n4*(1+n10)+6*n4-8*n8*(1+n10)+8*n8-6*
      (1+n10)^2+26*(1+n10)-20,1)
      +`Y'(1,1,1,-1+n4,2,0,1,n8,1,(1+n10),0,0,0,0)*rat(2*ep+2*(1+n10)-2,1)
      +`Y'(1,1,1,-1+n4,2,1,0,n8,1,(1+n10),0,0,0,0)*rat(2*ep+2*(1+n10)-2,1)
      +`Y'(1,1,1,-1+n4,2,1,1,n8,0,(1+n10),0,0,0,0)*rat(-2*ep-2*(1+n10)+2,1)
      +`Y'(1,1,1,n4,0,1,1,n8,1,(1+n10),0,0,0,0)*rat(-6*ep*n8+10*ep*(1+n10)-4*ep-2*n8^2+n8*(1+n10)+3*n8+4*(1+n10)^2-9*(1+n10)+3,1)
      +`Y'(1,1,1,n4,1,0,1,n8,1,(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(1,1,1,n4,1,1,0,n8,1,(1+n10),0,0,0,0)*rat(-10*ep*(1+n10)+6*ep-4*n8*(1+n10)+2*n8-4*(1+n10)^2+11*(1+n10)-5,1)
      +`Y'(1,1,1,n4,1,1,1,n8,0,(1+n10),0,0,0,0)*rat(10*ep*(1+n10)-10*ep+4*n8*(1+n10)-4*n8+4*(1+n10)^2-12*(1+n10)+8,1)
      +`Y'(1,1,1,n4,2,0,1,-1+n8,1,(1+n10),0,0,0,0)*rat(-2*ep-2*(1+n10)+2,1)
      +`Y'(1,1,1,n4,2,1,0,-1+n8,1,(1+n10),0,0,0,0)*rat(-2*ep-2*(1+n10)+2,1)
      +`Y'(1,1,1,n4,2,1,1,-1+n8,0,(1+n10),0,0,0,0)*rat(2*ep+2*(1+n10)-2,1)
      +`Y'(1,2,0,n4,1,1,1,-1+n8,1,(1+n10),0,0,0,0)*rat(6*ep+2*n8+3*(1+n10)-5,1)
      +`Y'(1,2,1,n4,1,1,1,-1+n8,0,(1+n10),0,0,0,0)*rat(-6*ep-2*n8-3*(1+n10)+5,1)
      +`Y'(0,1,1,n4,1,1,1,n8,1,1+(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(1,0,1,n4,1,1,1,n8,1,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Y'(1,1,0,n4,1,0,1,n8,1,(1+n10),0,0,0,0)*rat(-4*ep*(1+n10)+4*ep-2*n8*(1+n10)+2*n8-(1+n10)^2+4*(1+n10)-3,1)
      +`Y'(1,1,1,-1+n4,1,0,1,n8,1,(1+n10),0,0,0,0)*rat(-4*ep*(1+n10)+4*ep-2*n8*(1+n10)+2*n8-(1+n10)^2+4*(1+n10)-3,1)
      +`Y'(1,1,1,n4,0,0,1,n8,1,(1+n10),0,0,0,0)*rat(4*ep*(1+n10)-4*ep+2*n8*(1+n10)-2*n8+(1+n10)^2-4*(1+n10)+3,1)
      +`Y'(1,1,1,n4,1,-1,1,n8,1,(1+n10),0,0,0,0)*rat(-4*ep*(1+n10)+4*ep-2*n8*(1+n10)+2*n8-(1+n10)^2+4*(1+n10)-3,1)
      +`Y'(1,1,1,n4,1,0,0,n8,1,(1+n10),0,0,0,0)*rat(-4*ep*(1+n10)+4*ep-2*n8*(1+n10)+2*n8-(1+n10)^2+4*(1+n10)-3,1)
      +`Y'(1,1,1,n4,1,0,1,n8,0,(1+n10),0,0,0,0)*rat(4*ep*(1+n10)-4*ep+2*n8*(1+n10)-2*n8+(1+n10)^2-4*(1+n10)+3,1)
      +`Y'(0,1,0,n4,1,1,1,n8,1,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Y'(0,1,1,-1+n4,1,1,1,n8,1,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Y'(0,1,1,n4,0,1,1,n8,1,1+(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(0,1,1,n4,1,0,1,n8,1,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Y'(0,1,1,n4,1,1,0,n8,1,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Y'(0,1,1,n4,1,1,1,n8,0,1+(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(1,0,0,n4,1,1,1,n8,1,1+(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(1,0,1,-1+n4,1,1,1,n8,1,1+(1+n10),0,0,0,0)*rat(-2*ep*(1+n10)-2*n8*(1+n10)+(1+n10)^2+(1+n10),1)
      +`Y'(1,0,1,n4,0,1,1,n8,1,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Y'(1,0,1,n4,1,0,1,n8,1,1+(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(1,0,1,n4,1,1,0,n8,1,1+(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(1,0,1,n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-8*ep*(1+n10)-8*(1+n10)^2+8*(1+n10),1)
      +`Y'(1,0,1,n4,1,1,1,n8,0,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Y'(1,1,0,n4,1,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-2*ep*(1+n10)-2*n8*(1+n10)+(1+n10)^2+(1+n10),1)
      +`Y'(1,1,1,n4,0,1,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      +`Y'(1,1,1,n4,1,0,1,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(1,1,1,n4,1,1,0,-1+n8,1,1+(1+n10),0,0,0,0)*rat(-6*ep*(1+n10)-2*n8*(1+n10)-3*(1+n10)^2+5*(1+n10),1)
      +`Y'(1,1,1,n4,1,1,1,-1+n8,0,1+(1+n10),0,0,0,0)*rat(6*ep*(1+n10)+2*n8*(1+n10)+3*(1+n10)^2-5*(1+n10),1)
      );
*#endif
*--#] n10 : 
#endif
#if ( `il' >= 14 )
*--#[ n4 :
*
*   G16=
*`SKIP8'
id,ifmatch->nono1,
	`Z'(1,1,1,n4?{>1},1,1,1,n8?`POS8',1,0,0,0,0,0)
		 = -rat(1,-24*ep^2*(n4-1)-14*ep*(n4-1)^2-17*ep*(n4-1)*n8+20*ep*(n4-1)-2*(n4-1)^3-5*(n4-1)^2*n8+6*(n4-1)^2-3*(n4-1)*n8^2+7*(n4-1)*n8-4*(n4-1))*(
      +`Z'(1,1,1,-1+(n4-1),1,1,1,1+n8,1,0,0,0,0,0)*rat(12*ep^2*n8+7*ep*(n4-1)*n8+7*ep*n8^2-10*ep*n8+(n4-1)^2*n8+2*(n4-1)*n8^2-3*(n4-1)*n8+
      n8^3-3*n8^2+2*n8,1)
      +`Z'(1,1,1,(n4-1),1,1,1,1+n8,1,0,0,0,0,0)*rat(-12*ep^2*n8-7*ep*(n4-1)*n8-10*ep*n8^2+10*ep*n8-(n4-1)^2*n8-3*(n4-1)*n8^2+3*(n4-1)*n8-2*
      n8^3+4*n8^2-2*n8,1)
      +`Z'(1,1,1,(n4-1),1,1,1,n8,1,0,0,0,0,0)*rat(48*ep^3+52*ep^2*(n4-1)+40*ep^2*n8-64*ep^2+18*ep*(n4-1)^2+29*ep*(n4-1)*n8-46*ep*(n4-1)+11*ep
      *n8^2-36*ep*n8+28*ep+2*(n4-1)^3+5*(n4-1)^2*n8-8*(n4-1)^2+4*(n4-1)*n8^2-13*(n4-1)*n8+10*(n4-1)+n8^3-5*n8^2+8*n8-4,1)
      +`Y'(0,1,1,-1+(n4-1),1,1,1,2+n8,1,0,0,0,0,0)*rat(3*ep*n8^2+3*ep*n8+(n4-1)*n8^2+(n4-1)*n8+n8^3-n8,1)
      +`Y'(0,1,1,(n4-1),0,1,1,2+n8,1,0,0,0,0,0)*rat(-3*ep*n8^2-3*ep*n8-(n4-1)*n8^2-(n4-1)*n8-n8^3+n8,1)
      +`Y'(0,1,1,(n4-1),1,1,1,1+n8,1,0,0,0,0,0)*rat(-12*ep^2*n8-4*ep*(n4-1)*n8-7*ep*n8^2+7*ep*n8-(n4-1)*n8^2+(n4-1)*n8-n8^3+2*n8^2-n8,1)
      +`Y'(0,2,0,(n4-1),1,1,1,1+n8,1,0,0,0,0,0)*rat(3*ep*n8+(n4-1)*n8+n8^2-n8,1)
      +`Y'(0,2,1,(n4-1),1,1,1,1+n8,0,0,0,0,0,0)*rat(-3*ep*n8-(n4-1)*n8-n8^2+n8,1)
      +`Y'(1,0,1,-1+(n4-1),1,1,1,2+n8,1,0,0,0,0,0)*rat(-3*ep*n8^2-3*ep*n8-(n4-1)*n8^2-(n4-1)*n8-n8^3+n8,1)
      +`Y'(1,0,1,(n4-1),0,1,1,1+n8,2,0,0,0,0,0)*rat(-5*ep*n8-(n4-1)*n8-n8^2+2*n8,1)
      +`Y'(1,0,1,(n4-1),0,1,1,2+n8,1,0,0,0,0,0)*rat(3*ep*n8^2+3*ep*n8+(n4-1)*n8^2+(n4-1)*n8+n8^3-n8,1)
      +`Y'(1,0,1,(n4-1),1,0,1,1+n8,2,0,0,0,0,0)*rat(3*ep*n8+(n4-1)*n8+n8^2-n8,1)
      +`Y'(1,0,1,(n4-1),1,1,0,1+n8,2,0,0,0,0,0)*rat(3*ep*n8+(n4-1)*n8+n8^2-n8,1)
      +`Y'(1,0,1,(n4-1),1,1,1,1+n8,1,0,0,0,0,0)*rat(-8*ep^2*n8-10*ep*(n4-1)*n8-7*ep*n8^2+11*ep*n8-2*(n4-1)^2*n8-3*(n4-1)*n8^2+5*(n4-1)*n8-
      n8^3+4*n8^2-3*n8,1)
      +`Y'(1,0,2,-1+(n4-1),1,1,1,1+n8,1,0,0,0,0,0)*rat(-5*ep*n8-(n4-1)*n8-n8^2+2*n8,1)
      +`Y'(1,0,2,(n4-1),1,1,1,1+n8,1,0,0,0,0,0)*rat(7*ep*n8+(n4-1)*n8+n8^2-3*n8,1)
      +`Y'(1,1,0,1+(n4-1),1,1,1,n8,1,0,0,0,0,0)*rat(48*ep^2*(n4-1)+28*ep*(n4-1)^2+26*ep*(n4-1)*n8-40*ep*(n4-1)+4*(n4-1)^3+8*(n4-1)^2*n8-12*(n4-1)^2+4*(n4-1)*
      n8^2-11*(n4-1)*n8+8*(n4-1),2)
      +`Y'(1,1,0,(n4-1),0,1,1,1+n8,2,0,0,0,0,0)*rat(5*ep*n8+(n4-1)*n8+n8^2-2*n8,1)
      +`Y'(1,1,0,(n4-1),1,0,1,1+n8,2,0,0,0,0,0)*rat(-3*ep*n8-(n4-1)*n8-n8^2+n8,1)
      +`Y'(1,1,0,(n4-1),1,1,0,1+n8,2,0,0,0,0,0)*rat(-3*ep*n8-(n4-1)*n8-n8^2+n8,1)
      +`Y'(1,1,0,(n4-1),1,1,1,1+n8,1,0,0,0,0,0)*rat(48*ep^2*n8+30*ep*(n4-1)*n8+28*ep*n8^2-54*ep*n8+4*(n4-1)^2*n8+8*(n4-1)*n8^2-15*(n4-1)*n8+4
      *n8^3-14*n8^2+14*n8,2)
      +`Y'(1,1,1,-1+(n4-1),2,0,1,1+n8,1,0,0,0,0,0)*rat(-3*ep*n8-(n4-1)*n8-n8^2+n8,1)
      +`Y'(1,1,1,-1+(n4-1),2,1,0,1+n8,1,0,0,0,0,0)*rat(-3*ep*n8-(n4-1)*n8-n8^2+n8,1)
      +`Y'(1,1,1,-1+(n4-1),2,1,1,1+n8,0,0,0,0,0,0)*rat(-2*ep*n8+n8,2)
      +`Y'(1,1,1,1+(n4-1),0,1,1,1+n8,1,0,0,0,0,0)*rat(8*ep*(n4-1)*n8+2*(n4-1)^2*n8+2*(n4-1)*n8^2-3*(n4-1)*n8,2)
      +`Y'(1,1,1,(n4-1),0,1,1,1+n8,1,0,0,0,0,0)*rat(-32*ep^2*n8-16*ep*(n4-1)*n8-14*ep*n8^2+28*ep*n8-2*(n4-1)^2*n8-4*(n4-1)*n8^2+7*(n4-1)*n8-2
      *n8^3+6*n8^2-6*n8,2)
      +`Y'(1,1,1,(n4-1),0,1,1,2+n8,1,-1,0,0,0,0)*rat(2*ep*n8^2+2*ep*n8-n8^2-n8,2)
      +`Y'(1,1,1,(n4-1),1,0,1,1+n8,1,0,0,0,0,0)*rat(12*ep^2*n8+4*ep*(n4-1)*n8+10*ep*n8^2-7*ep*n8+2*(n4-1)*n8^2-(n4-1)*n8+2*n8^3-3*n8^2+n8
      ,1)
      +`Y'(1,1,1,(n4-1),1,1,0,1+n8,1,0,0,0,0,0)*rat(12*ep^2*n8+4*ep*(n4-1)*n8+10*ep*n8^2-10*ep*n8+2*(n4-1)*n8^2-2*(n4-1)*n8+2*n8^3-4*n8^2
      +2*n8,1)
      +`Y'(1,1,1,(n4-1),1,1,1,1+n8,0,0,0,0,0,0)*rat(-8*ep^2*n8-4*ep*(n4-1)*n8+10*ep*n8+2*(n4-1)*n8-3*n8,2)
      +`Y'(1,1,1,(n4-1),2,0,1,n8,1,0,0,0,0,0)*rat(3*ep*n8+(n4-1)*n8+n8^2-n8,1)
      +`Y'(1,1,1,(n4-1),2,1,0,n8,1,0,0,0,0,0)*rat(3*ep*n8+(n4-1)*n8+n8^2-n8,1)
      +`Y'(1,1,1,(n4-1),2,1,1,n8,0,0,0,0,0,0)*rat(2*ep*n8-n8,2)
      +`Y'(1,1,2,(n4-1),1,1,1,1+n8,0,0,0,0,0,0)*rat(-2*ep*n8+n8,1)
      +`Y'(1,2,0,(n4-1),1,1,1,1+n8,1,-1,0,0,0,0)*rat(-12*ep*n8-2*(n4-1)*n8-2*n8^2+5*n8,2)
      +`Y'(1,2,0,(n4-1),1,1,1,n8,1,0,0,0,0,0)*rat(24*ep^2+14*ep*(n4-1)+9*ep*n8-20*ep+2*(n4-1)^2+3*(n4-1)*n8-6*(n4-1)+n8^2-4*n8+4,1)
      +`Y'(1,2,1,(n4-1),1,1,1,1+n8,0,-1,0,0,0,0)*rat(2*ep*n8-n8,2)
      +`Y'(1,2,1,(n4-1),1,1,1,n8,0,0,0,0,0,0)*rat(-24*ep^2-14*ep*(n4-1)-9*ep*n8+20*ep-2*(n4-1)^2-3*(n4-1)*n8+6*(n4-1)-n8^2+4*n8-4,1)
      );
*--#] n4 : 
#endif
#if ( `il' >= 15 )
*--#[ n8 : 
id,ifmatch->nono1,
	`Z'(1,1,1,1,1,1,1,n8?{>1},1,0,0,0,0,0)
		 = -rat(1,-8*ep*(n8-1)-2*(n8-1)^2+2*(n8-1))*(
      +`Z'(1,1,1,1,1,1,1,(n8-1),1,0,0,0,0,0)*rat(-96*ep^3-88*ep^2*(n8-1)+104*ep^2-24*ep*(n8-1)^2+60*ep*(n8-1)-36*ep-2*(n8-1)^3+8*(n8-1)^2-10*(n8-1)+4,ep+(n8-1))
      +`Y'(0,1,1,0,1,1,1,2+(n8-1),1,0,0,0,0,0)*rat(-2*(n8-1)^2-2*(n8-1),1)
      +`Y'(0,1,1,1,0,1,1,2+(n8-1),1,0,0,0,0,0)*rat(2*(n8-1)^2+2*(n8-1),1)
      +`Y'(0,1,1,1,1,1,1,1+(n8-1),1,0,0,0,0,0)*rat(8*ep*(n8-1)+2*(n8-1)^2-2*(n8-1),1)
      +`Y'(0,2,0,1,1,1,1,1+(n8-1),1,0,0,0,0,0)*rat(-2*(n8-1),1)
      +`Y'(0,2,1,1,1,1,1,1+(n8-1),0,0,0,0,0,0)*rat(2*(n8-1),1)
      +`Y'(1,0,1,0,1,1,1,2+(n8-1),1,0,0,0,0,0)*rat(2*(n8-1)^2+2*(n8-1),1)
      +`Y'(1,0,1,1,0,1,1,1+(n8-1),2,0,0,0,0,0)*rat(10*ep*(n8-1)+2*(n8-1)^2-2*(n8-1),3*ep+(n8-1))
      +`Y'(1,0,1,1,0,1,1,2+(n8-1),1,0,0,0,0,0)*rat(-2*(n8-1)^2-2*(n8-1),1)
      +`Y'(1,0,1,1,1,0,1,1+(n8-1),2,0,0,0,0,0)*rat(-2*(n8-1),1)
      +`Y'(1,0,1,1,1,1,0,1+(n8-1),2,0,0,0,0,0)*rat(-2*(n8-1),1)
      +`Y'(1,0,1,1,1,1,1,1+(n8-1),1,0,0,0,0,0)*rat(16*ep^2*(n8-1)+14*ep*(n8-1)^2-2*ep*(n8-1)+2*(n8-1)^3-2*(n8-1)^2,3*ep+(n8-1))
      +`Y'(1,0,2,0,1,1,1,1+(n8-1),1,0,0,0,0,0)*rat(10*ep*(n8-1)+2*(n8-1)^2-2*(n8-1),3*ep+(n8-1))
      +`Y'(1,0,2,1,1,1,1,1+(n8-1),1,0,0,0,0,0)*rat(-14*ep*(n8-1)-2*(n8-1)^2+4*(n8-1),3*ep+(n8-1))
      +`Y'(1,1,0,1,0,1,1,1+(n8-1),2,0,0,0,0,0)*rat(-10*ep*(n8-1)-2*(n8-1)^2+2*(n8-1),3*ep+(n8-1))
      +`Y'(1,1,0,1,1,0,1,1+(n8-1),2,0,0,0,0,0)*rat(2*(n8-1),1)
      +`Y'(1,1,0,1,1,1,0,1+(n8-1),2,0,0,0,0,0)*rat(2*(n8-1),1)
      +`Y'(1,1,0,1,1,1,1,1+(n8-1),1,0,0,0,0,0)*rat(-48*ep^3*(n8-1)-76*ep^2*(n8-1)^2+48*ep^2*(n8-1)-32*ep*(n8-1)^3+47*ep*(n8-1)^2-9*ep*(n8-1)-4*(n8-1)^4+9
      *(n8-1)^3-5*(n8-1)^2,3*ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,1,0,2,1,1,1,(n8-1),1,0,0,0,0,0)*rat(-16*ep^2*(n8-1)-7*ep*(n8-1)^2+5*ep*(n8-1)-(n8-1)^3+(n8-1)^2,3*ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,1,1,0,1,1,1,1+(n8-1),1,0,0,0,0,0)*rat(-8*ep^2*(n8-1)+6*ep*(n8-1)^2+2*ep*(n8-1)+4*(n8-1)^3-2*(n8-1)^2,ep+(n8-1))
      +`Y'(1,1,1,0,1,1,1,2+(n8-1),1,-1,0,0,0,0)*rat(16*ep*(n8-1)^2+16*ep*(n8-1)+6*(n8-1)^3+2*(n8-1)^2-4*(n8-1),ep+(n8-1))
      +`Y'(1,1,1,0,2,0,1,1+(n8-1),1,0,0,0,0,0)*rat(2*(n8-1),1)
      +`Y'(1,1,1,0,2,1,0,1+(n8-1),1,0,0,0,0,0)*rat(2*(n8-1),1)
      +`Y'(1,1,1,0,2,1,1,1+(n8-1),0,0,0,0,0,0)*rat(-22*ep^2*(n8-1)-15*ep*(n8-1)^2+5*ep*(n8-1)-3*(n8-1)^3+(n8-1)^2,3*ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,1,1,0,2,1,1,1+(n8-1),1,0,-1,0,0,0)*rat(-16*ep*(n8-1)-6*(n8-1)^2+4*(n8-1),ep+(n8-1))
      +`Y'(1,1,1,1,0,1,1,1+(n8-1),1,0,0,0,0,0)*rat(-64*ep^3*(n8-1)-70*ep^2*(n8-1)^2+36*ep^2*(n8-1)-30*ep*(n8-1)^3+23*ep*(n8-1)^2-5*ep*(n8-1)-4*(n8-1)^4+5
      *(n8-1)^3-(n8-1)^2,3*ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,1,1,1,0,1,1,2+(n8-1),1,-1,0,0,0,0)*rat(-26*ep^2*(n8-1)^2-26*ep^2*(n8-1)-19*ep*(n8-1)^3-12*ep*(n8-1)^2+7*ep*(n8-1)-3*(n8-1)^4+3*(n8-1)^2,3*
      ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,1,1,1,1,0,1,1+(n8-1),1,0,0,0,0,0)*rat(8*ep*(n8-1)+2*(n8-1)^2-2*(n8-1),1)
      +`Y'(1,1,1,1,1,1,0,1+(n8-1),1,0,0,0,0,0)*rat(8*ep*(n8-1)+2*(n8-1)^2,1)
      +`Y'(1,1,1,1,1,1,1,1+(n8-1),0,0,0,0,0,0)*rat(8*ep^3*(n8-1)+8*ep^2*(n8-1)^2-30*ep^2*(n8-1)-23*ep*(n8-1)^2+7*ep*(n8-1)-3*(n8-1)^3+3*(n8-1)^2,3*ep^2+4
      *ep*(n8-1)+(n8-1)^2)
      +`Y'(1,1,1,1,2,0,1,(n8-1),1,0,0,0,0,0)*rat(16*ep+4*(n8-1)-4,1)
      +`Y'(1,1,1,1,2,1,0,(n8-1),1,0,0,0,0,0)*rat(16*ep+4*(n8-1)-4,1)
      +`Y'(1,1,1,1,2,1,1,(n8-1),0,0,0,0,0,0)*rat(-48*ep^3-60*ep^2*(n8-1)+12*ep^2-25*ep*(n8-1)^2+11*ep*(n8-1)-3*(n8-1)^3+3*(n8-1)^2,3*ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,1,1,2,0,1,1,1+(n8-1),1,0,0,0,0,0)*rat(16*ep^2*(n8-1)+7*ep*(n8-1)^2-5*ep*(n8-1)+(n8-1)^3-(n8-1)^2,3*ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,1,2,1,1,1,1,1+(n8-1),0,0,0,0,0,0)*rat(4*ep*(n8-1)-2*(n8-1),3*ep+(n8-1))
      +`Y'(1,2,0,1,1,1,1,1+(n8-1),1,-1,0,0,0,0)*rat(36*ep^2*(n8-1)+31*ep*(n8-1)^2-9*ep*(n8-1)+5*(n8-1)^3-5*(n8-1)^2,3*ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,2,0,1,1,1,1,(n8-1),1,0,0,0,0,0)*rat(96*ep^3+84*ep^2*(n8-1)-72*ep^2+32*ep*(n8-1)^2-44*ep*(n8-1)+12*ep+4*(n8-1)^3-8*(n8-1)^2+4*(n8-1),3*
      ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,2,1,1,1,1,1,1+(n8-1),0,-1,0,0,0,0)*rat(-26*ep^2*(n8-1)-19*ep*(n8-1)^2+7*ep*(n8-1)-3*(n8-1)^3+3*(n8-1)^2,3*ep^2+4*ep*(n8-1)+(n8-1)^2)
      +`Y'(1,2,1,1,1,1,1,(n8-1),0,0,0,0,0,0)*rat(-96*ep^3-84*ep^2*(n8-1)+72*ep^2-32*ep*(n8-1)^2+44*ep*(n8-1)-12*ep-4*(n8-1)^3+8*(n8-1)^2-4*(n8-1),3*
      ep^2+4*ep*(n8-1)+(n8-1)^2)
      );
*--#] n8 : 
#endif
Goto nono2;
Label nono1;
$ired = 1;
Label nono2;
id	`Z'(?a,x?neg0_,?b,n10?,n11?,n12?,n13?,n14?) = `Y'(?a,x,?b,n10,n11,n12,n13,n14);
*id	`Y'(?a) = 0;
*Print +f +s;
ModuleOption,maximum,$ired;
.sort:Main nono loop, pass `$pass';
#if ( `$ired' > 0 )
	#redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
id	`Z'(1,1,1,1,1,1,1,1,1,0,0,0,0,0) = Master(nono);
*
#break
*--#] nono : 
*--#[ cross :
#case cross
*
*--#[ Loop :
#do il = 1,17
#$pass = 0;
#do ired = 1,1
#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering cross loop `il', pass `$pass' at `time_' sec."
#endif
*#if ( `il' >= 1 )
* causes infinite loop!
#if 0
*--#[ stat file :
*
*	Note: the computer generated file had to be edited by hand in the
*	expressions B5 and B13, because otherwise they would cause loops.
*	This might be helped by making CleanupBasic sensitive to the order
*	of the elimination of the variables. This makes the design of the
*	scheme more work, but on the other hand one might get better schemes.
*
*  Prework to make the reduction of the cross topology faster.
*
* B2
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n5+1)*(`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,1+n11,n12,n13,n14)*rat(-n5+1,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,n2,-1+
   n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n5,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-
   2*ep-2*n1-n5-n7-n9-n13+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-1+n5,1));
* B8
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n6+1)*(`Z'(-1+n1,n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,
   n8,n9,n10,n11,1+n12,n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,
   n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n6,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*
   rat(-2*ep-2*n2-n6-n7-n10-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n6,1));
* B14
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?neg_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n5+1)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-1+n5,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,
   1+n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n5+1,1)+`Z'(n1,n2,n3,-
   1+n4,-1+n5,n6,n7,1+n8,1+n9,n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(
   -2*ep-2*n3-n5-n8-n11-n14+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-1+n5,1));
* B20
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?neg0_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n6+1)*(`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-1+n6,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,
   n9,1+n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,1+n8,n9,1+n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,
   n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*
   rat(-2*ep-2*n4-n6-n8-n12-n13+5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-1+n6,1));
* B6
id,ifmatch->cross1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,n13,n14)*rat(-n1,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+
   n9,1+n10,n11,n12,n13,n14)*rat(-2*n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,-1+n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,-1+n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,
   n14)*rat(8*ep+2*n1+2*n2+2*n3+2*n4+2*n5+2*n6+2*n7+2*n8+n9+n10+n11+n12+2*n13+2*n14-17,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,2+n10,n11,n12,n13,n14)*rat(-2*n10-2,1));
* B16
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n4+1)*(`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+
   n9,n10,n11,1+n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,
   n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,1+n12,n13,n14)*rat(-2*n11,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,
   n14)*rat(8*ep+2*n1+2*n2+2*n3+2*n4+2*n5+2*n6+2*n7+2*n8+n9+n10+n11+n12+2*n13+2*n14-17,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,
   n9,n10,n11,2+n12,n13,n14)*rat(-2*n12-2,1));
* B3
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n5+1)*(`Z'(-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,
   n10,n11,n12,n13,1+n14)*rat(-n5+1,1)+`Z'(1+n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n1,1)+`Z'(1+n1,n2,
   n3,n4,-1+n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n1,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+
   n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(-n5+1,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,
   n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,1+n9,-1+n10,n11,n12,
   n13,1+n14)*rat(-n9,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,
   n7,n8,n9,n10,n11,-1+n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+
   `Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n1+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,
   1+n14)*rat(-1+n5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-1+n5,1));
* B7
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n6+1)*(`Z'(-1+n1,1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n2,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,1+n7,
   n8,n9,n10,n11,n12,1+n13,n14)*rat(n7,1)+`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(-1+
   n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6+1,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,1+
   n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6+1,1)+`Z'(n1,1+n2,n3,n4,n5,-1+n6,-1+
   n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n2,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(
   n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,1+n10,n11,
   n12,1+n13,n14)*rat(-n10,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,-
   1+n6,n7,n8,n9,n10,-1+n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(
   n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n2+n7,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,
   n12,1+n13,n14)*rat(-1+n6,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-1+n6,1));
* B15
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n5+1)*(`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,n7,1+
   n8,n9,n10,n11,n12,1+n13,n14)*rat(-n8,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n5+1,1)+`Z'(n1,n2,
   1+n3,-1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,n4,-1+n5,n6,n7,-1+n8,n9,n10,n11,n12,1+
   n13,n14)*rat(n3,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,n7,1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,-1+n4,-1+n5,n6,
   n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n5+1,1)+`Z'(
   n1,n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,-1+n10,n11,n12,
   1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,-1+n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,-1+n5,
   n6,n7,n8,n9,n10,n11,-1+n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,
   1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n3+n8,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,
   1+n13,n14)*rat(-1+n5,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,n14)*rat(-1+n5,1));
* B19
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n6+1)*(`Z'(-1+n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,-1+n3,1+n4,n5,-1+n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(n8,1)+`Z'(n1,
   n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+
   n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(-n6+1,1)+`Z'(n1,n2,n3,1+n4,n5,-1+n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(n1,
   n2,n3,n4,-1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,-1+n9,n10,n11,n12,1+
   n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,-1+n11,1+n12,n13,1+n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,-1+
   n6,n7,n8,n9,n10,-1+n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,
   1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n4+n8,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,
   n13,1+n14)*rat(-1+n6,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-1+n6,1));
* B1
id,ifmatch->cross1,`Z'(n1?{>1},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n1+1)*(`Z'(-2+n1,n2,n3,n4,n5,n6,1+n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n7,1)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,1+n7,n8,1+
   n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,1+n8,1+n9,n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(-1+n1,n2,
   n3,-1+n4,n5,n6,n7,1+n8,1+n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,n8,1+n9,-1+n10,n11,n12,n13,
   n14)*rat(n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9,
   n10,n11,-1+n12,1+n13,n14)*rat(-n13,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(4*ep+2*n1+2*n3+2*n5+
   n7+n8+n9+n11+n13+n14-9,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,2+n9,n10,n11,n12,n13,n14)*rat(-2*n9-2,1));
* B11
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n3+1)*(`Z'(-1+n1,n2,-1+n3,n4,n5,n6,1+n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,1+n7,
   n8,n9,n10,1+n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,n2,-2+n3,n4,n5,n6,n7,1+n8,n9,n10,1+n11,n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,-
   1+n3,-1+n4,n5,n6,n7,1+n8,n9,n10,1+n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+n8,n9,n10,1+n11,-1+n12,n13,
   n14)*rat(n8,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,-1+
   n10,1+n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(4*ep+2*n1+2*n3+2*n5+
   n7+n8+n9+n11+n13+n14-9,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,2+n11,n12,n13,n14)*rat(-2*n11-2,1));
* B9
id,ifmatch->cross1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-n2+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(1+n1,-1+n2,-1+n3,n4,n5,n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(n1,1)+`Z'(1+n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n1,1)+`Z'(1+n1,
   -1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-n1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,
   1+n14)*rat(-n1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n1,1)+`Z'(n1,-1+n2,-1+n3,1+n4,n5,n6,n7,
   n8,n9,n10,n11,n12,n13,1+n14)*rat(n4,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-2*n13,1)+`Z'(
   n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,2+n14)*rat(-2*n14-2,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,
   n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,1+n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(-n4,1)+`Z'(n1,-1+n2,n3,n4,
   -1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,1+n14)*
   rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,
   n10,-1+n11,n12,n13,1+n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,-1+n11,n12,n13,1+n14)*rat(-n10,1)+`Z'(n1,-1+
   n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,1+n14)*rat(n12,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,
   1+n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(2*ep+n1+2*n3+n4+n11+n14-3,1));
* B17
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-n4+1)*(`Z'(-1+n1,1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n2,1)+`Z'(-1+n1,n2,1+n3,-1+n4,n5,n6,n7,
   n8,n9,n10,n11,n12,1+n13,n14)*rat(n3,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-2*n14,1)+`Z'(-
   1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,2+n13,n14)*rat(-2*n13-2,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,
   n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,1+n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n2,1)+`Z'(n1,n2,-1+n3,-1+
   n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,1+n3,-1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*
   rat(-n3,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,n9,
   n10,-1+n11,n12,1+n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n3,1)+`Z'(n1,n2,n3,-
   1+n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+
   n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,1+n13,n14)*rat(n10,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-
   1+n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(
   n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,
   1+n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,
   n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(2*ep+2*n1+n2+n3+n9+n13-3,1));
* B5
id,ifmatch->cross1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(n1,1)+`Z'(1+n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,-1+n13,n14)*rat(-n1,1)+`Z'(n1,-2+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+
   n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,
   n13,n14)*rat(n3,1)+`Z'(n1,-1+n2,1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-n3,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,
   n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-2*n13,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-2*
   n14,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,
   n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,1+n12,n13,n14)*rat(2*n9,1)+`Z'(n1,-1+n2,n3,
   n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-n9,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,1+n14)*rat(
   n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n10,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,
   n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,
   n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(n11,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)
   +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-6*ep-2*n1-n2-n3-2*n5-2*n6-2*n7-2*n8-n9-n10-n11-n13-2*
   n14+13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,
   n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n2,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,
   n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n2+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n2,1));
* B13
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n4+1)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(1+n1,-1+n2,n3,-1+n4,n5,n6,n7,
   n8,n9,1+n10,n11,n12,n13,n14)*rat(n1,1)+`Z'(1+n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n1,1)+`Z'(n1,-
   1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-2*n13,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,
   n13,1+n14)*rat(-2*n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-1+n4,1)+`Z'(n1,n2,1+n3,-1+n4,n5,
   n6,n7,n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(n3,1)+`Z'(n1,n2,1+n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,-1+n14)*rat(-n3,
   1)+`Z'(n1,n2,n3,-2+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,n8,n9,1+n10,n11,
   n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,-1+n4,n5,
   n6,n7,n8,-1+n9,1+n10,n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(n9,1)+
   `Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,1+n11,n12,n13,n14)*rat(2*n11,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,-1+
   n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n12,1)+`Z'(n1,n2,n3,-1+n4,n5,
   n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-6*ep-n1-
   2*n3-n4-2*n5-2*n6-2*n7-2*n8-n9-n11-n12-2*n13-n14+13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(-
   n11,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,
   n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,
   n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-1+n4,1));
* B10
id,ifmatch->cross1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,-n2+1)*(`Z'(-1+n1,-1+n2,n3,n4,n5,n6,1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(n1,-
   1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,n9,n10,n11,1+n12,
   n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,1+
   n7,n8,n9,n10,n11,1+n12,-1+n13,n14)*rat(n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(n1,
   -1+n2,n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,1+
   n14)*rat(n14,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n10,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,
   n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,1+n12,n13,1+n14)*rat(n14,1)+`Z'(n1,-1+
   n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*ep-n2-2*n6-n7-n10-n14+5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,
   n11,1+n12,n13,n14)*rat(-1+n2,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n2+1,1)+`Z'(n1,n2,n3,n4,n5,
   n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n2,1));
* B18
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,-n4+1)*(`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,1+
   n8,n9,1+n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-n13,1)+`Z'(n1,-
   1+n2,n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-1+n4,1)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,1+n8,n9,1+n10,n11,n12,n13,
   n14)*rat(n8,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,n7,1+n8,n9,1+n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,
   n8,n9,1+n10,n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9,1+n10,-1+n11,n12,n13,n14)*rat(-n8,1)+`Z'(
   n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9,1+n10,n11,n12,n13,-1+n14)*rat(n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9,1+n10,n11,n12,
   n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,-1+
   n9,1+n10,n11,n12,1+n13,n14)*rat(n13,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n12,1)+`Z'(n1,n2,
   n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-2*ep-n4-2*n6-n8-n12-n13+5,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,
   n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n4+1,1)+`Z'(n1,n2,n3,n4,n5,
   n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-1+n4,1));
* B4
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_) =
  -rat(1,-1+n7)*(`Z'(-1+n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n7,1)+`Z'(1+n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(n1,1)+`Z'(1+n1,n2,n3,n4,-1+n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n1,1)+`Z'(1+n1,n2,
   n3,n4,n5,n6,-1+n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-n1,1)+`Z'(1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,-1+n11,n12,n13,1+
   n14)*rat(-n1,1)+`Z'(1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n1,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,1+n8,
   n9,n10,n11,n12,n13,1+n14)*rat(n8,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,-
   1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+n7,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,1+n8,n9,n10,n11,n12,n13,1+
   n14)*rat(-n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,
   n9,n10,n11,n12,n13,1+n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,-1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n13,1)+`Z'(n1,
   n2,n3,n4,n5,n6,-1+n7,n8,1+n9,n10,-1+n11,n12,n13,1+n14)*rat(-n9,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,1+
   n14)*rat(2*ep+n1+2*n3+n8+n11+n14-3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,
   n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-n7+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-1+
   n7,1));
* B12
id,ifmatch->cross1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,-1+n8)*(`Z'(-1+n1,n2,1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n3,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,1+n7,-1+
   n8,n9,n10,n11,n12,1+n13,n14)*rat(n7,1)+`Z'(-1+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(-1+
   n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n8,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,1+n7,-1+n8,n9,n10,n11,n12,1+
   n13,n14)*rat(-n7,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,-1+n3,n4,n5,n6,
   n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-1+n8,1)+`Z'(n1,n2,1+n3,n4,-1+n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-n3,1)+
   `Z'(n1,n2,1+n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,-1+n8,n9,n10,-1+n11,
   n12,1+n13,n14)*rat(-n3,1)+`Z'(n1,n2,1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n3,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,
   n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n8+1,1)+`Z'(n1,n2,n3,n4,n5,n6,-1+n7,-1+n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n14,1)+
   `Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,-1+n9,n10,1+n11,n12,1+n13,n14)*rat(-n11,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,
   1+n13,n14)*rat(2*ep+2*n1+n3+n7+n9+n13-3,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-n8+1,1)+`Z'(n1,
   n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-n8+1,1)+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*
   rat(-1+n8,1));
*--#] stat file : 
#endif
#if ( `il' >= 2 )
*--#[ n2 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(n1?`POS1',n2?{>1},n3?`POS3',n4?`POS4',n5?`POS5',n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(n7,2*ep*(n2-1)+2*(n2-1)^2+2*(n2-1)*n7+(n2-1)*n10+(n2-1)*n14-2*(n2-1))*(
      +`Z'(-1+n1,1+(n2-1),n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*(n2-1),1)
      +`Z'(n1,(n2-1),n3,n4,n5,1+n6,n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(n6*n14,n7)
      +`Z'(n1,(n2-1),n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep*n6-2*n6^2-n6*n10-n6*n14+2*n6,n7)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,1+n7,n8,-1+n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,1+n7,n8,-1+n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,1+n7,n8,n9,n10,-1+n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(2*ep+2*(n2-1)+2*n7+n10+n14-2,1)
      +`Z'(n1,(n2-1),n3,n4,n5,1+n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n6*n10,n7)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)
      +`Z'(-1+n1,(n2-1),n3,n4,n5,n6,1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n10,1)
      +`Z'(-1+n1,(n2-1),n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,-1+(n2-1),n3,n4,n5,1+n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n6*n10,n7)
      +`Z'(n1,-1+(n2-1),n3,n4,n5,n6,1+n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(n10,1)
      +`Z'(n1,(n2-1),-1+n3,n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n6*n14,n7)
      +`Z'(n1,(n2-1),-1+n3,n4,n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,(n2-1),n3,-1+n4,n5,1+n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-n6*n10,n7)
      +`Z'(n1,(n2-1),n3,-1+n4,n5,1+n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n6*n14,n7)
      +`Z'(n1,(n2-1),n3,n4,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,(n2-1),n3,n4,n5,1+n6,n7,-1+n8,n9,n10,n11,n12,n13,1+n14)*rat(n6*n14,n7)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(2*ep*n10+2*(n2-1)*n10+n6*n10+n7*n10+n10^2+n10*n14-3*n10,n7)
      );
*--#] n2 : 
#endif
#if ( `il' >= 3 )
*--#[ n4 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(n1?`POS1',1,n3?`POS3',n4?{>1},n5?`POS5',n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(n8,2*ep*(n4-1)+2*(n4-1)^2+2*(n4-1)*n8+(n4-1)*n12+(n4-1)*n13-2*(n4-1))*(
      +`Z'(n1,1,-1+n3,1+(n4-1),n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(-2*(n4-1),1)
      +`Z'(n1,1,n3,(n4-1),n5,1+n6,n7,n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(n6*n13,n8)
      +`Z'(n1,1,n3,(n4-1),n5,1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-2*ep*n6-2*n6^2-n6*n12-n6*n13+2*n6,n8)
      +`Z'(n1,1,n3,(n4-1),n5,n6,n7,1+n8,-1+n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'(n1,1,n3,(n4-1),n5,n6,n7,1+n8,n9,n10,-1+n11,1+n12,n13,n14)*rat(-n12,1)
      +`Z'(n1,1,n3,(n4-1),n5,n6,n7,1+n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'(n1,1,n3,(n4-1),n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(2*ep+2*(n4-1)+2*n8+n12+n13-2,1)
      +`Z'(n1,1,n3,(n4-1),n5,1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n6*n12,n8)
      +`Z'(n1,1,n3,(n4-1),n5,n6,n7,1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)
      +`Z'(-1+n1,1,n3,(n4-1),n5,1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6*n13,n8)
      +`Z'(-1+n1,1,n3,(n4-1),n5,n6,n7,1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,1,-1+n3,(n4-1),n5,n6,n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)
      +`Z'(n1,1,-1+n3,(n4-1),n5,n6,n7,1+n8,n9,n10,n11,n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,1,n3,-1+(n4-1),n5,1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n6*n12,n8)
      +`Z'(n1,1,n3,-1+(n4-1),n5,n6,n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n12,1)
      +`Z'(n1,1,n3,(n4-1),-1+n5,n6,n7,1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'(n1,1,n3,(n4-1),n5,1+n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n6*n13,n8)
      +`Z'(n1,1,n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(2*ep*n12+2*(n4-1)*n12+n6*n12+n8*n12+n12^2+n12*n13-3*n12,n8)
      +`Y'(n1,0,n3,(n4-1),n5,1+n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n6*n12,n8)
      +`Y'(n1,0,n3,(n4-1),n5,1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n6*n13,n8)
      );
*--#] n4 : 
#endif
#if ( `il' >= 4 )
*--#[ n14 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(n1?`POS1',1,n3?`POS3',1,n5?`POS5',n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?neg_)
			 = -rat(1,2*ep+2*n1+n5+n7+n9+n13-4)*(
      +`Z'(1+n1,1,n3,1,n5,n6,-1+n7,n8,n9,n10,-1+n11,n12,n13,(1+n14))*rat(n1,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,1+n9,-1+n10,-1+n11,n12,n13,(1+n14))*rat(-n9,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,-1+n10,-1+n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,(1+n14))*rat(-2*ep-2*n1-n5-n7-n9-n13+4,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,n10,-1+n11,-1+n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,(1+n14))*rat(n7-1,1)
      +`Z'(1+n1,1,n3,1,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n1,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,1+n9,-1+n10,n11,n12,n13,(1+n14))*rat(n9,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n7+1,1)
      +`Z'(-1+n1,1,n3,1,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n7,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,1+n9,-1+n10,n11,n12,n13,(1+n14))*rat(-n9,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,n13,(1+n14))*rat(n9,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(-2*n13,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(4*ep+2*n1+3*n5+2*n7+2*n9+2*n13-9,1)
      +`Z'(1+n1,1,-1+n3,1,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n1,1)
      +`Z'(n1,1,-1+n3,1,n5,n6,n7,n8,1+n9,-1+n10,n11,n12,n13,(1+n14))*rat(n9,1)
      +`Z'(n1,1,-1+n3,1,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,1,-1+n3,1,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,1,-1+n3,1,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n7+1,1)
      +`Z'(n1,1,n3,1,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n1-1,1)
      +`Z'(n1,1,n3,1,n5,-1+n6,n7,n8,n9,n10,-1+n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,1,n3,1,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-2*ep-n1-n5-n7-n9-n13+4,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,(1+n14))*rat(-2*n9,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,1,-1+n3,1,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,1,n3,1,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(-1+n1,1,-1+n3,1,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(-1+n1,1,n3,1,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,1,-1+n3,1,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(1+n1,0,n3,1,n5,n6,n7,n8,n9,n10,-1+n11,n12,n13,(1+n14))*rat(-2*n1,1)
      +`Y'(1+n1,0,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(2*n1,1)
      +`Y'(1+n1,0,-1+n3,1,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(2*n1,1)
      +`Y'(1+n1,0,n3,1,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n1,1)
      +`Y'(n1,0,n3,1,-1+n5,n6,1+n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n7,1)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,n13,(1+n14))*rat(-n9,1)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n5+n7,1)
      +`Y'(n1,1,n3,0,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,(1+n14))*rat(2*n9,1)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(n1,1,n3,0,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(-1+n1,0,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(-1+n1,1,n3,0,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(2*n13,1)
      +`Y'(n1,0,-1+n3,1,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-2*n13,1)
      +`Y'(n1,0,n3,0,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(n1,0,n3,1,n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(n1,1,-1+n3,0,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      );
*--#] n14 : 
#endif
#if ( `il' >= 5 )
*--#[ n13 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(n1?`POS1',1,n3?`POS3',1,n5?`POS5',n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?neg_,0)
			 = -rat(1,2*ep+n6+n7+n10-2)*(
      +`Z'(-1+n1,2,n3,1,n5,n6,n7,n8,n9,n10,n11,-1+n12,(1+n13),0)*rat(-4*ep-2*n7-2*n10+2,2*ep+2*n7+n10-2)
      +`Z'(n1,1,n3,1,n5,1+n6,-1+n7,n8,n9,n10,n11,-1+n12,(1+n13),0)*rat(2*ep*n6+2*n6^2+n6*n10-2*n6,2*ep+2*n7+n10-2)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,-1+n9,1+n10,n11,-1+n12,(1+n13),0)*rat(-2*ep*n10-n7*n10-n10^2+n10,2*ep+2*n7+n10-2)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,-1+n9,n10,n11,n12,(1+n13),0)*rat(-2*ep-n6-n7-n10+2,1)
      +`Z'(-1+n1,2,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(2,1)
      +`Z'(n1,1,n3,1,n5,1+n6,-1+n7,n8,n9,1+n10,n11,-1+n12,(1+n13),0)*rat(-n6*n10,2*ep+2*n7+n10-2)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,(1+n13),0)*rat(n10,1)
      +`Z'(n1,1,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n7+1,1)
      +`Z'(n1,2,n3,1,n5,n6,-1+n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-1,1)
      +`Z'(-1+n1,1,n3,1,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n7,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,n9,1+n10,n11,-1+n12,(1+n13),0)*rat(-2*ep*n10-n7*n10-n10^2+n10,2*ep+2*n7+n10-2)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n6+n7,1)
      +`Z'(-1+n1,2,n3,1,n5,-1+n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-1,1)
      +`Z'(n1,1,n3,1,n5,n6,-1+n7,n8,n9,1+n10,n11,-1+n12,(1+n13),0)*rat(-2*ep*n10-n6*n10-n7*n10-n10^2+2*n10,2*ep+2*n7+n10-2)
      +`Z'(n1,1,n3,1,n5,n6,-1+n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-2*ep-n6-n7-n10+3,1)
      +`Z'(-1+n1,1,n3,1,n5,n6,n7,n8,n9,1+n10,n11,n12,(1+n13),0)*rat(2*n10,1)
      +`Y'(-1+n1,2,n3,0,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(2,1)
      +`Y'(n1,0,n3,1,n5,-1+n6,1+n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(n7,1)
      +`Y'(n1,0,n3,1,n5,1+n6,-1+n7,n8,n9,1+n10,n11,-1+n12,(1+n13),0)*rat(-n6*n10,2*ep+2*n7+n10-2)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,(1+n13),0)*rat(-n10,1)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,n9,1+n10,n11,-1+n12,(1+n13),0)*rat(2*ep*n10+n7*n10+n10^2-n10,2*ep+2*n7+n10-2)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(4*ep+3*n6+2*n7+2*n10-7,1)
      +`Y'(n1,1,n3,0,n5,1+n6,-1+n7,n8,n9,1+n10,n11,-1+n12,(1+n13),0)*rat(n6*n10,2*ep+2*n7+n10-2)
      +`Y'(n1,1,n3,0,n5,n6,n7,n8,-1+n9,1+n10,n11,n12,(1+n13),0)*rat(n10,1)
      +`Y'(n1,1,n3,0,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-n7+1,1)
      +`Y'(n1,2,n3,0,n5,n6,-1+n7,n8,n9,n10,n11,n12,(1+n13),0)*rat(-1,1)
      +`Y'(n1,0,n3,1,n5,n6,n7,n8,n9,1+n10,n11,n12,(1+n13),0)*rat(-2*n10,1)
      );
*--#] n13 : 
#endif
#if ( `il' >= 6 )
*--#[ n6x :
`SKIP6'
id,ifmatch->cross1,
	`Z'(n1?`POS1',1,n3?`POS3',1,n5?`POS5',n6?{>1},n7?`POS7',n8?`POS8'
		,n9?`NEG9',n10?!{n12?},n11?`NEG11',n12?!{n10?},0,0)
			 = -rat(1,(n6-1)*n10-(n6-1)*n12)*(
      +`Z'(n1,1,n3,1,n5,1+(n6-1),n7,n8,n9,1+n10,n11,n12,0,0)*rat(-(n6-1)*n10,1)
      +`Z'(n1,1,n3,1,n5,1+(n6-1),n7,n8,n9,n10,n11,1+n12,0,0)*rat((n6-1)*n12,1)
      +`Z'(-1+n1,1,n3,1,n5,(n6-1),1+n7,n8,n9,n10,n11,1+n12,0,0)*rat(n7*n12,1)
      +`Z'(n1,1,-1+n3,1,n5,(n6-1),n7,1+n8,n9,1+n10,n11,n12,0,0)*rat(-n8*n10,1)
      +`Z'(n1,1,n3,1,n5,(n6-1),n7,n8,n9,1+n10,n11,n12,0,0)*rat(2*ep*n10+(n6-1)*n10+n8*n10+n10*n12-2*n10,1)
      +`Z'(n1,1,n3,1,n5,(n6-1),n7,n8,n9,n10,n11,1+n12,0,0)*rat(-2*ep*n12-(n6-1)*n12-n7*n12-n10*n12+2*n12,1)
      +`Y'(n1,0,n3,1,n5,1+(n6-1),n7,n8,n9,1+n10,n11,n12,0,0)*rat(-(n6-1)*n10,1)
      +`Y'(n1,0,n3,1,n5,1+(n6-1),n7,n8,n9,n10,n11,1+n12,0,0)*rat(-(n6-1)*n12,1)
      +`Y'(n1,0,n3,1,n5,(n6-1),1+n7,n8,n9,n10,n11,1+n12,0,0)*rat(-n7*n12,1)
      +`Y'(n1,1,n3,0,n5,1+(n6-1),n7,n8,n9,1+n10,n11,n12,0,0)*rat((n6-1)*n10,1)
      +`Y'(n1,1,n3,0,n5,1+(n6-1),n7,n8,n9,n10,n11,1+n12,0,0)*rat((n6-1)*n12,1)
      +`Y'(n1,1,n3,0,n5,(n6-1),n7,1+n8,n9,1+n10,n11,n12,0,0)*rat(n8*n10,1)
      );
*--#] n6x : 
#endif
#if ( `il' >= 7 )
*--#[ n1 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(n1?{>1},1,n3?`POS3',1,n5?`POS5',n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?`NEG12',0,0)
			 = -rat(n7,2*ep*(n1-1)+2*(n1-1)^2+2*(n1-1)*n7-2*(n1-1))*(
      +`Z'((n1-1),1,n3,1,1+n5,n6,n7,n8,n9,n10,n11,n12,0,0)*rat(-2*ep*n5-2*n5^2-n5*n9+2*n5,n7)
      +`Z'((n1-1),1,n3,1,n5,n6,1+n7,n8,n9,n10,n11,n12,0,0)*rat(2*ep+2*(n1-1)+2*n7-2,1)
      +`Z'((n1-1),1,n3,1,1+n5,n6,n7,n8,1+n9,n10,n11,n12,0,0)*rat(n5*n9,n7)
      +`Z'(-1+(n1-1),1,n3,1,1+n5,n6,n7,n8,1+n9,n10,n11,n12,0,0)*rat(n5*n9,n7)
      +`Z'(-1+(n1-1),1,n3,1,n5,n6,1+n7,n8,1+n9,n10,n11,n12,0,0)*rat(2*n9,1)
      +`Z'((n1-1),1,-1+n3,1,1+n5,n6,n7,n8,1+n9,n10,n11,n12,0,0)*rat(-n5*n9,n7)
      +`Z'((n1-1),1,-1+n3,1,n5,n6,n7,1+n8,1+n9,n10,n11,n12,0,0)*rat(n8*n9,n7)
      +`Z'((n1-1),1,n3,1,n5,n6,n7,n8,1+n9,n10,n11,n12,0,0)*rat(6*ep*n9+4*(n1-1)*n9+2*n3*n9+3*n5*n9+2*n7*n9+n8*n9+2*n9^2+n9*n11-
      10*n9,n7)
      +`Z'((n1-1),1,n3,1,n5,n6,n7,n8,2+n9,n10,n11,n12,0,0)*rat(-2*n9^2-2*n9,n7)
      +`Y'(1+(n1-1),0,n3,1,n5,n6,1+n7,n8,n9,n10,n11,n12,0,0)*rat(-2*(n1-1),1)
      +`Y'((n1-1),0,n3,1,n5,n6,1+n7,n8,1+n9,n10,n11,n12,0,0)*rat(-2*n9,1)
      +`Y'((n1-1),1,n3,0,n5,n6,n7,1+n8,1+n9,n10,n11,n12,0,0)*rat(-n8*n9,n7)
      );
*--#] n1 : 
#endif
#if ( `il' >= 8 )
*--#[ n3 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,n3?{>1},1,n5?`POS5',n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?`NEG12',0,0)
			 = -rat(n8,2*ep*(n3-1)+2*(n3-1)^2+2*(n3-1)*n8-2*(n3-1))*(
      +`Z'(1,1,(n3-1),1,1+n5,n6,n7,n8,n9,n10,n11,n12,0,0)*rat(-2*ep*n5-2*n5^2-n5*n11+2*n5,n8)
      +`Z'(1,1,(n3-1),1,n5,n6,n7,1+n8,n9,n10,n11,n12,0,0)*rat(2*ep+2*(n3-1)+2*n8-2,1)
      +`Z'(1,1,(n3-1),1,1+n5,n6,n7,n8,n9,n10,1+n11,n12,0,0)*rat(n5*n11,n8)
      +`Z'(1,1,-1+(n3-1),1,1+n5,n6,n7,n8,n9,n10,1+n11,n12,0,0)*rat(n5*n11,n8)
      +`Z'(1,1,-1+(n3-1),1,n5,n6,n7,1+n8,n9,n10,1+n11,n12,0,0)*rat(2*n11,1)
      +`Z'(1,1,(n3-1),1,n5,n6,n7,n8,n9,n10,1+n11,n12,0,0)*rat(6*ep*n11+4*(n3-1)*n11+3*n5*n11+n7*n11+2*n8*n11+n9*n11+2*n11^2-8*
      n11,n8)
      +`Z'(1,1,(n3-1),1,n5,n6,n7,n8,n9,n10,2+n11,n12,0,0)*rat(-2*n11^2-2*n11,n8)
      +`Y'(1,1,1+(n3-1),0,n5,n6,n7,1+n8,n9,n10,n11,n12,0,0)*rat(-2*(n3-1),1)
      +`Y'(0,1,(n3-1),1,1+n5,n6,n7,n8,n9,n10,1+n11,n12,0,0)*rat(-n5*n11,n8)
      +`Y'(0,1,(n3-1),1,n5,n6,1+n7,n8,n9,n10,1+n11,n12,0,0)*rat(n7*n11,n8)
      +`Y'(1,0,(n3-1),1,n5,n6,1+n7,n8,n9,n10,1+n11,n12,0,0)*rat(-n7*n11,n8)
      +`Y'(1,1,(n3-1),0,n5,n6,n7,1+n8,n9,n10,1+n11,n12,0,0)*rat(-2*n11,1)
      );
*--#] n3 : 
#endif
#if ( `il' >= 9 )
*--#[ n5x :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,1,1,n5?{>1},n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?!{n11?},n10?`NEG10',n11?!{n9?},n12?`NEG12',0,0)
			 = -rat(1,-(n5-1)*n9+(n5-1)*n11)*(
      +`Z'(1,1,1,1,1+(n5-1),n6,n7,n8,1+n9,n10,n11,n12,0,0)*rat((n5-1)*n9,1)
      +`Z'(1,1,1,1,1+(n5-1),n6,n7,n8,n9,n10,1+n11,n12,0,0)*rat(-(n5-1)*n11,1)
      +`Z'(1,1,1,1,(n5-1),n6,n7,n8,1+n9,n10,n11,n12,0,0)*rat(-2*ep*n9-(n5-1)*n9-n8*n9-n9*n11+2*n9,1)
      +`Z'(1,1,1,1,(n5-1),n6,n7,n8,n9,n10,1+n11,n12,0,0)*rat(2*ep*n11+(n5-1)*n11+n7*n11+n9*n11-2*n11,1)
      +`Y'(0,1,1,1,1+(n5-1),n6,n7,n8,1+n9,n10,n11,n12,0,0)*rat((n5-1)*n9,1)
      +`Y'(0,1,1,1,1+(n5-1),n6,n7,n8,n9,n10,1+n11,n12,0,0)*rat((n5-1)*n11,1)
      +`Y'(0,1,1,1,(n5-1),n6,1+n7,n8,n9,n10,1+n11,n12,0,0)*rat(n7*n11,1)
      +`Y'(1,0,1,1,(n5-1),n6,1+n7,n8,n9,n10,1+n11,n12,0,0)*rat(-n7*n11,1)
      +`Y'(1,1,0,1,1+(n5-1),n6,n7,n8,1+n9,n10,n11,n12,0,0)*rat(-(n5-1)*n9,1)
      +`Y'(1,1,0,1,1+(n5-1),n6,n7,n8,n9,n10,1+n11,n12,0,0)*rat(-(n5-1)*n11,1)
      +`Y'(1,1,0,1,(n5-1),n6,n7,1+n8,1+n9,n10,n11,n12,0,0)*rat(-n8*n9,1)
      +`Y'(1,1,1,0,(n5-1),n6,n7,1+n8,1+n9,n10,n11,n12,0,0)*rat(n8*n9,1)
      );
*--#] n5x : 
#endif
#if ( `il' >= 10 )
*--#[ n12 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,1,1,n5?`POS5',n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?neg_,0,0)
			 = -rat(1,2*ep+n6+n8+(1+n12)-3)*(
      +`Z'(1,1,1,1,n5,n6,n7,n8,n9,-1+n10,n11,(1+n12),0,0)*rat(-2*ep-n6-n7-n10+3,1)
      +`Z'(1,1,1,1,n5,n6,n7,n8,n9,n10,n11,(1+n12),0,0)*rat(n7-n8+n10-(1+n12),1)
      +`Y'(0,1,1,1,n5,n6,1+n7,n8,-1+n9,n10,n11,(1+n12),0,0)*rat(n7,1)
      +`Y'(1,0,1,1,n5,n6,1+n7,n8,-1+n9,n10,n11,(1+n12),0,0)*rat(-n7,1)
      +`Y'(1,1,0,1,n5,n6,n7,1+n8,n9,n10,-1+n11,(1+n12),0,0)*rat(-n8,1)
      +`Y'(1,1,1,0,n5,n6,n7,1+n8,n9,n10,-1+n11,(1+n12),0,0)*rat(n8,1)
      +`Y'(1,1,2,0,n5,n6,n7,n8,n9,n10,-1+n11,(1+n12),0,0)*rat(1,1)
      +`Y'(2,0,1,1,n5,n6,n7,n8,-1+n9,n10,n11,(1+n12),0,0)*rat(-1,1)
      +`Y'(0,1,1,1,n5,n6,1+n7,n8,n9,n10,n11,(1+n12),0,0)*rat(-n7,1)
      +`Y'(1,0,1,1,n5,n6,1+n7,n8,n9,n10,n11,(1+n12),0,0)*rat(n7,1)
      +`Y'(1,1,0,1,n5,n6,n7,1+n8,n9,n10,n11,(1+n12),0,0)*rat(n8,1)
      +`Y'(1,1,1,0,n5,n6,n7,1+n8,n9,n10,n11,(1+n12),0,0)*rat(-n8,1)
      +`Y'(-1,1,1,1,n5,n6,1+n7,n8,n9,n10,n11,(1+n12),0,0)*rat(-n7,1)
      +`Y'(0,0,1,1,n5,n6,1+n7,n8,n9,n10,n11,(1+n12),0,0)*rat(n7,1)
      +`Y'(0,1,0,1,n5,n6,1+n7,n8,n9,n10,n11,(1+n12),0,0)*rat(n7,1)
      +`Y'(0,1,0,1,n5,n6,n7,1+n8,n9,n10,n11,(1+n12),0,0)*rat(-n8,1)
      +`Y'(0,1,1,0,n5,n6,n7,1+n8,n9,n10,n11,(1+n12),0,0)*rat(n8,1)
      +`Y'(0,1,1,1,n5,n6,n7,n8,n9,n10,n11,(1+n12),0,0)*rat(-4*ep-2*n5-n7-n8-n9-n11+6,1)
      +`Y'(1,0,0,1,n5,n6,1+n7,n8,n9,n10,n11,(1+n12),0,0)*rat(-n7,1)
      +`Y'(1,0,1,1,n5,n6,n7,n8,n9,n10,n11,(1+n12),0,0)*rat(8*ep+2*n5+2*n6+2*n7+2*n8+n9+n10+n11+(1+n12)-10,1)
      +`Y'(1,1,-1,1,n5,n6,n7,1+n8,n9,n10,n11,(1+n12),0,0)*rat(n8,1)
      +`Y'(1,1,0,0,n5,n6,n7,1+n8,n9,n10,n11,(1+n12),0,0)*rat(-n8,1)
      +`Y'(1,1,0,1,n5,n6,n7,n8,n9,n10,n11,(1+n12),0,0)*rat(4*ep+2*n5+n7+n8+n9+n11-6,1)
      +`Y'(1,1,1,0,n5,n6,n7,n8,n9,n10,n11,(1+n12),0,0)*rat(-8*ep-2*n5-2*n6-2*n7-2*n8-n9-n10-n11-(1+n12)+10,1)
      +`Y'(0,1,1,1,n5,n6,n7,n8,1+n9,n10,n11,(1+n12),0,0)*rat(2*n9,1)
      +`Y'(1,0,1,1,n5,n6,n7,n8,1+n9,n10,n11,(1+n12),0,0)*rat(-2*n9,1)
      +`Y'(1,1,0,1,n5,n6,n7,n8,n9,n10,1+n11,(1+n12),0,0)*rat(-2*n11,1)
      +`Y'(1,1,1,0,n5,n6,n7,n8,n9,n10,1+n11,(1+n12),0,0)*rat(2*n11,1)
      );
*--#] n12 : 
#endif
#if ( `il' >= 11 )
*--#[ n11 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,1,1,n5?`POS5',n6?`POS6',n7?`POS7',n8?`POS8'
		,n9?`NEG9',n10?`NEG10',n11?neg_,0,0,0)
			 = -rat(1,2*ep+n5+n8+(1+n11)-3)*(
      +`Z'(1,1,1,1,n5,n6,n7,n8,-1+n9,n10,(1+n11),0,0,0)*rat(-2*ep-n5-n7-n9+3,1)
      +`Z'(1,1,1,1,n5,n6,n7,n8,n9,n10,(1+n11),0,0,0)*rat(n7-n8+n9-(1+n11),1)
      +`Y'(0,1,1,1,n5,n6,1+n7,n8,-1+n9,n10,(1+n11),0,0,0)*rat(-n7,1)
      +`Y'(1,0,1,1,n5,n6,1+n7,n8,-1+n9,n10,(1+n11),0,0,0)*rat(n7,1)
      +`Y'(1,1,0,1,n5,n6,n7,1+n8,n9,n10,-1+(1+n11),0,0,0)*rat(n8,1)
      +`Y'(1,1,1,0,n5,n6,n7,1+n8,n9,n10,-1+(1+n11),0,0,0)*rat(-n8,1)
      +`Y'(0,1,1,1,n5,n6,1+n7,n8,n9,n10,(1+n11),0,0,0)*rat(n7,1)
      +`Y'(1,0,1,1,n5,n6,1+n7,n8,n9,n10,(1+n11),0,0,0)*rat(-n7,1)
      +`Y'(1,1,0,1,n5,n6,n7,1+n8,n9,n10,(1+n11),0,0,0)*rat(-n8,1)
      +`Y'(1,1,1,0,n5,n6,n7,1+n8,n9,n10,(1+n11),0,0,0)*rat(n8,1)
      +`Y'(-1,1,1,1,n5,n6,1+n7,n8,n9,n10,(1+n11),0,0,0)*rat(n7,1)
      +`Y'(0,0,1,1,n5,n6,1+n7,n8,n9,n10,(1+n11),0,0,0)*rat(-n7,1)
      +`Y'(0,1,0,1,n5,n6,1+n7,n8,n9,n10,(1+n11),0,0,0)*rat(-n7,1)
      +`Y'(0,1,0,1,n5,n6,n7,1+n8,n9,n10,(1+n11),0,0,0)*rat(n8,1)
      +`Y'(0,1,1,0,n5,n6,n7,1+n8,n9,n10,(1+n11),0,0,0)*rat(-n8,1)
      +`Y'(0,1,1,1,n5,n6,n7,n8,n9,n10,(1+n11),0,0,0)*rat(4*ep+2*n5+n7+n8+n9+(1+n11)-6,1)
      +`Y'(1,0,0,1,n5,n6,1+n7,n8,n9,n10,(1+n11),0,0,0)*rat(n7,1)
      +`Y'(1,1,-1,1,n5,n6,n7,1+n8,n9,n10,(1+n11),0,0,0)*rat(-n8,1)
      +`Y'(1,1,0,0,n5,n6,n7,1+n8,n9,n10,(1+n11),0,0,0)*rat(n8,1)
      +`Y'(1,1,0,1,n5,n6,n7,n8,n9,n10,(1+n11),0,0,0)*rat(-4*ep-2*n5-n7-n8-n9-(1+n11)+6,1)
      );
*--#] n11 : 
#endif
#if ( `il' >= 12 )
*--#[ n8 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,1,1,n5?`POS5',n6?`POS6',n7?`POS7',n8?{>1}
		,n9?`NEG9',n10?`NEG10',0,0,0,0)
			 = -rat(ep+(n8-1)-1,2*ep^2*(n8-1)+3*ep*(n8-1)^2-ep*(n8-1)+(n8-1)^3-(n8-1)^2)*(
      +`Z'(1,1,1,1,1+n5,1+n6,-1+n7,(n8-1),n9,n10,0,0,0,0)*rat(2*ep*n5*n6+2*n5^2*n6+n5*n6*n9-2*n5*n6,2*ep+2*n7-2)
      +`Z'(1,1,1,1,1+n5,1+n6,n7,-1+(n8-1),n9,n10,0,0,0,0)*rat(2*ep^2*n5*n6+2*ep*n5^2*n6+ep*n5*n6^2+ep*n5*n6*(n8-1)-4*ep*n5*n6+
      n5^2*n6^2+n5^2*n6*(n8-1)-2*n5^2*n6-n5*n6^2-n5*n6*(n8-1)+2*n5*n6,ep^2+2*ep*(n8-1)-2*ep+(n8-1)^2-2*(n8-1)+1)
      +`Z'(1,1,1,1,n5,1+n6,n7,(n8-1),n9,n10,0,0,0,0)*rat(4*ep*n6+2*n5*n6+n6*n7+n6*(n8-1)+n6*n9-4*n6,1)
      +`Z'(1,1,1,1,1+n5,1+n6,-1+n7,(n8-1),1+n9,n10,0,0,0,0)*rat(-n5*n6*n9,2*ep+2*n7-2)
      +`Z'(1,1,1,1,n5,1+n6,n7,(n8-1),1+n9,n10,0,0,0,0)*rat(-n6*n9,1)
      +`Z'(1,1,1,1,n5,1+n6,-1+n7,(n8-1),1+n9,n10,0,0,0,0)*rat(-6*ep*n6*n9-3*n5*n6*n9-2*n6*n7*n9-n6*(n8-1)*n9-2*n6*n9^2+6*n6*n9,
      2*ep+2*n7-2)
      +`Z'(1,1,1,1,n5,n6,n7,(n8-1),1+n9,n10,0,0,0,0)*rat(2*ep*n9+n6*n9+(n8-1)*n9-2*n9,1)
      +`Z'(1,1,1,1,n5,1+n6,-1+n7,(n8-1),2+n9,n10,0,0,0,0)*rat(n6*n9^2+n6*n9,ep+n7-1)
      +`Y'(1,1,0,2,n5,n6,n7,1+(n8-1),n9,n10,0,0,0,0)*rat(-ep*(n8-1),ep+(n8-1)-1)
      +`Y'(1,1,2,0,n5,1+n6,n7,(n8-1),n9,n10,0,0,0,0)*rat(-ep^2*n6+n6^2*(n8-1)-n6^2-n6*(n8-1)+n6,ep^2+2*ep*(n8-1)-2*ep+(n8-1)^2-2*(n8-1)+1)
      +`Y'(1,1,2,0,n5,n6,n7,1+(n8-1),n9,n10,0,0,0,0)*rat(-(n8-1),1)
      +`Y'(2,0,1,1,n5,1+n6,n7,(n8-1),n9,n10,0,0,0,0)*rat(-ep*n6,ep+n7-1)
      +`Y'(0,1,1,1,1+n5,1+n6,-1+n7,(n8-1),1+n9,n10,0,0,0,0)*rat(-n5*n6*n9,2*ep+2*n7-2)
      +`Y'(0,1,1,1,n5,1+n6,n7,(n8-1),1+n9,n10,0,0,0,0)*rat(-n6*n7*n9+n6*n9,ep+n7-1)
      +`Y'(1,0,1,1,n5,1+n6,n7,(n8-1),1+n9,n10,0,0,0,0)*rat(-ep*n6*n9,ep+n7-1)
      +`Y'(1,1,0,1,1+n5,1+n6,-1+n7,(n8-1),1+n9,n10,0,0,0,0)*rat(n5*n6*n9,2*ep+2*n7-2)
      +`Y'(1,1,0,1,n5,1+n6,-1+n7,1+(n8-1),1+n9,n10,0,0,0,0)*rat(-n6*(n8-1)*n9,2*ep+2*n7-2)
      +`Y'(1,1,0,1,n5,n6,n7,1+(n8-1),1+n9,n10,0,0,0,0)*rat(-(n8-1)*n9,1)
      +`Y'(1,1,1,0,n5,1+n6,-1+n7,1+(n8-1),1+n9,n10,0,0,0,0)*rat(n6*(n8-1)*n9,2*ep+2*n7-2)
      +`Y'(1,1,1,0,n5,1+n6,n7,(n8-1),1+n9,n10,0,0,0,0)*rat(n6*n9,1)
      +`Y'(1,1,1,0,n5,n6,n7,1+(n8-1),1+n9,n10,0,0,0,0)*rat((n8-1)*n9,1)
      );
*--#] n8 : 
#endif
#if ( `il' >= 13 )
*--#[ n7 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,1,1,n5?`POS5',n6?`POS6',n7?{>1},1,n9?`NEG9',n10?`NEG10',0,0,0,0)
		 = -rat(2*ep+2*(n7-1)+n10-2,
      4*ep^2*(n7-1)+6*ep*(n7-1)^2+4*ep*(n7-1)*n10-2*ep*(n7-1)+2*(n7-1)^3+3*(n7-1)^2*n10-2*(n7-1)^2+(n7-1)*n10^2-(n7-1)*n10)*(
      +`Z'(1,1,1,1,1+n5,1+n6,-1+(n7-1),1,n9,n10,0,0,0,0)*rat(4*ep^2*n5*n6+4*ep*n5^2*n6+2*ep*n5*n6^2+2*ep*n5*n6*(n7-1)+2*ep*n5*
      n6*n9+2*ep*n5*n6*n10-8*ep*n5*n6+2*n5^2*n6^2+2*n5^2*n6*(n7-1)+2*n5^2*n6*n10-4*n5^2*n6+n5*n6^2*n9-2*n5*n6^2+n5*n6*(n7-1)*n9
      -2*n5*n6*(n7-1)+n5*n6*n9*n10-2*n5*n6*n9-2*n5*n6*n10+4*n5*n6,2*ep^2+4*ep*(n7-1)+ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Z'(1,1,1,1,1+n5,n6,(n7-1),1,-1+n9,1+n10,0,0,0,0)*rat(-4*ep^2*n5*n10-4*ep*n5^2*n10-2*ep*n5*(n7-1)*n10-2*ep*n5*n9*n10-2*
      ep*n5*n10^2+8*ep*n5*n10-2*n5^2*(n7-1)*n10-2*n5^2*n10^2+2*n5^2*n10-n5*(n7-1)*n9*n10+3*n5*(n7-1)*n10-n5*n9*n10^2+n5*n9*n10+3*n5
      *n10^2-3*n5*n10,4*ep^2+8*ep*(n7-1)+2*ep*n10-4*ep+4*(n7-1)^2+2*(n7-1)*n10-4*(n7-1))
      +`Z'(1,1,1,1,n5,1+n6,(n7-1),1,n9,n10,0,0,0,0)*rat(4*ep*n6+2*n5*n6+n6*(n7-1)+n6*n9-3*n6,1)
      +`Z'(1,1,1,1,1+n5,1+n6,-1+(n7-1),1,1+n9,n10,0,0,0,0)*rat(-2*ep*n5*n6*n9-n5*n6^2*n9-n5*n6*(n7-1)*n9-n5*n6*n9*n10+2*n5*n6*
      n9,2*ep^2+4*ep*(n7-1)+ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Z'(1,1,1,1,1+n5,n6,(n7-1),1,n9,1+n10,0,0,0,0)*rat(2*ep*n5*n9*n10-2*ep*n5*n10+n5*(n7-1)*n9*n10-n5*(n7-1)*n10+n5*n9*n10^2-n5*
      n9*n10-n5*n10^2+n5*n10,4*ep^2+8*ep*(n7-1)+2*ep*n10-4*ep+4*(n7-1)^2+2*(n7-1)*n10-4*(n7-1))
      +`Z'(1,1,1,1,n5,1+n6,(n7-1),1,1+n9,n10,0,0,0,0)*rat(-n6*n9,1)
      +`Z'(1,1,1,1,n5,1+n6,(n7-1),1,n9,1+n10,0,0,0,0)*rat(-n6*(n7-1)*n10+n6*n10,2*ep+2*(n7-1)+n10-2)
      +`Z'(2,1,1,1,n5,1+n6,-1+(n7-1),1,n9,1+n10,0,0,0,0)*rat(-n6*n10,2*ep+2*(n7-1)+n10-2)
      +`Z'(1,1,1,1,n5,1+n6,-1+(n7-1),1,1+n9,n10,0,0,0,0)*rat(-12*ep^2*n6*n9-6*ep*n5*n6*n9-6*ep*n6^2*n9-10*ep*n6*(n7-1)*n9-4*ep*
      n6*n9^2-6*ep*n6*n9*n10+22*ep*n6*n9-3*n5*n6^2*n9-3*n5*n6*(n7-1)*n9-3*n5*n6*n9*n10+6*n5*n6*n9-2*n6^2*(n7-1)*n9-2*n6^2*n9^2+
      5*n6^2*n9-2*n6*(n7-1)^2*n9-2*n6*(n7-1)*n9^2-2*n6*(n7-1)*n9*n10+9*n6*(n7-1)*n9-2*n6*n9^2*n10+4*n6*n9^2+5*n6*n9*n10-10*n6*n9,2*ep^2
      +4*ep*(n7-1)+ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Z'(1,1,1,1,n5,n6,(n7-1),1,1+n9,n10,0,0,0,0)*rat(-2*ep*n9-n6*n9-(n7-1)*n9-n9*n10+2*n9,1)
      +`Z'(1,1,1,1,n5,n6,(n7-1),1,n9,1+n10,0,0,0,0)*rat(-4*ep^2*(n7-1)*n10+12*ep^2*n9*n10-8*ep^2*n10+6*ep*n5*n9*n10-6*ep*n5*n10
      -2*ep*n6*(n7-1)*n10+2*ep*n6*n10-6*ep*(n7-1)^2*n10+10*ep*(n7-1)*n9*n10-2*ep*(n7-1)*n10^2+4*ep*n9^2*n10+6*ep*n9*n10^2-20*ep*n9*n10-
      4*ep*n10^2+12*ep*n10+3*n5*(n7-1)*n9*n10-3*n5*(n7-1)*n10+3*n5*n9*n10^2-3*n5*n9*n10-3*n5*n10^2+3*n5*n10-2*n6*(n7-1)^2*n10+2*n6*
      (n7-1)*n10-2*(n7-1)^3*n10+2*(n7-1)^2*n9*n10-2*(n7-1)^2*n10^2+4*(n7-1)^2*n10+2*(n7-1)*n9^2*n10+2*(n7-1)*n9*n10^2-9*(n7-1)*n9*n10+3*(n7-1)*n10+2*n9^2*
      n10^2-2*n9^2*n10-7*n9*n10^2+7*n9*n10+5*n10^2-5*n10,4*ep^2+8*ep*(n7-1)+2*ep*n10-4*ep+4*(n7-1)^2+2*(n7-1)*n10-4*(n7-1))
      +`Z'(2,1,1,1,n5,n6,-1+(n7-1),1,n9,1+n10,0,0,0,0)*rat(-2*ep*n10-n6*n10-(n7-1)*n10-n10^2+2*n10,2*ep+2*(n7-1)+n10-2)
      +`Z'(1,1,1,1,n5,1+n6,-1+(n7-1),1,2+n9,n10,0,0,0,0)*rat(4*ep*n6*n9^2+4*ep*n6*n9+2*n6^2*n9^2+2*n6^2*n9+2*n6*(n7-1)*n9^2+2*
      n6*(n7-1)*n9+2*n6*n9^2*n10-4*n6*n9^2+2*n6*n9*n10-4*n6*n9,2*ep^2+4*ep*(n7-1)+ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Z'(1,1,1,1,n5,n6,(n7-1),1,1+n9,1+n10,0,0,0,0)*rat(-2*ep*n9^2*n10+2*ep*n9*n10-(n7-1)*n9^2*n10+(n7-1)*n9*n10-n9^2*n10^2+n9^2*
      n10+n9*n10^2-n9*n10,2*ep^2+4*ep*(n7-1)+ep*n10-2*ep+2*(n7-1)^2+(n7-1)*n10-2*(n7-1))
      +`Y'(0,2,1,1,n5,n6,1+(n7-1),1,n9,n10,0,0,0,0)*rat(-2*ep*(n7-1)-(n7-1)*n10,2*ep+2*(n7-1)+n10-2)
      +`Y'(1,1,1,1,1+n5,1+n6,(n7-1),0,n9,n10,0,0,0,0)*rat(ep*n5*n6+n5^2*n6-n5*n6,ep)
      +`Y'(1,1,2,0,n5,1+n6,(n7-1),1,n9,n10,0,0,0,0)*rat(-n6,1)
      +`Y'(2,0,1,1,n5,1+n6,(n7-1),1,n9,n10,0,0,0,0)*rat(-2*ep^2*n6-ep*n6*n10+2*n6^2*(n7-1)-2*n6^2+n6*(n7-1)*n10-2*n6*(n7-1)-n6*n10+2*
      n6,2*ep^2+4*ep*(n7-1)+ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Y'(2,0,1,1,n5,n6,1+(n7-1),1,-1+n9,1+n10,0,0,0,0)*rat(-2*ep*(n7-1)*n10-(n7-1)^2*n10-(n7-1)*n10^2+(n7-1)*n10,2*ep^2+4*ep*(n7-1)+ep*n10-2
      *ep+2*(n7-1)^2+(n7-1)*n10-2*(n7-1))
      +`Y'(2,0,1,1,n5,n6,1+(n7-1),1,n9,n10,0,0,0,0)*rat(-(n7-1),1)
      +`Y'(0,1,1,1,1+n5,1+n6,-1+(n7-1),1,1+n9,n10,0,0,0,0)*rat(-2*ep*n5*n6*n9-n5*n6^2*n9-n5*n6*(n7-1)*n9-n5*n6*n9*n10+2*n5*n6*
      n9,2*ep^2+4*ep*(n7-1)+ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Y'(0,1,1,1,1+n5,n6,(n7-1),1,n9,1+n10,0,0,0,0)*rat(2*ep*n5*n9*n10-2*ep*n5*n10+n5*(n7-1)*n9*n10-n5*(n7-1)*n10+n5*n9*n10^2-n5*
      n9*n10-n5*n10^2+n5*n10,4*ep^2+8*ep*(n7-1)+2*ep*n10-4*ep+4*(n7-1)^2+2*(n7-1)*n10-4*(n7-1))
      +`Y'(0,1,1,1,n5,1+n6,(n7-1),1,1+n9,n10,0,0,0,0)*rat(-4*ep*n6*(n7-1)*n9+4*ep*n6*n9-2*n6^2*(n7-1)*n9+2*n6^2*n9-2*n6*(n7-1)^2*n9-2*
      n6*(n7-1)*n9*n10+6*n6*(n7-1)*n9+2*n6*n9*n10-4*n6*n9,2*ep^2+4*ep*(n7-1)+ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Y'(0,1,1,1,n5,n6,1+(n7-1),1,1+n9,n10,0,0,0,0)*rat((n7-1)*n9,1)
      +`Y'(0,1,1,1,n5,n6,1+(n7-1),1,n9,1+n10,0,0,0,0)*rat(ep*(n7-1)^2*n10+2*ep*(n7-1)*n9*n10-3*ep*(n7-1)*n10+(n7-1)^3*n10+(n7-1)^2*n9*n10-2*
      (n7-1)^2*n10+(n7-1)*n9*n10^2-(n7-1)*n9*n10-(n7-1)*n10^2+(n7-1)*n10,2*ep^2+4*ep*(n7-1)+ep*n10-2*ep+2*(n7-1)^2+(n7-1)*n10-2*(n7-1))
      +`Y'(1,0,1,1,n5,1+n6,(n7-1),1,1+n9,n10,0,0,0,0)*rat(-2*ep^2*n6*n9-ep*n6*n9*n10+2*n6^2*(n7-1)*n9-2*n6^2*n9+n6*(n7-1)*n9*n10-2*
      n6*(n7-1)*n9-n6*n9*n10+2*n6*n9,2*ep^2+4*ep*(n7-1)+ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Y'(1,0,1,1,n5,1+n6,(n7-1),1,n9,1+n10,0,0,0,0)*rat(-n6*(n7-1)*n10+n6*n10,2*ep+2*(n7-1)+n10-2)
      +`Y'(1,0,1,1,n5,n6,1+(n7-1),1,1+n9,n10,0,0,0,0)*rat(-(n7-1)*n9,1)
      +`Y'(1,0,1,1,n5,n6,1+(n7-1),1,n9,1+n10,0,0,0,0)*rat(-ep*(n7-1)^2*n10-2*ep*(n7-1)*n9*n10+3*ep*(n7-1)*n10-(n7-1)^3*n10-(n7-1)^2*n9*n10+2*
      (n7-1)^2*n10-(n7-1)*n9*n10^2+(n7-1)*n9*n10+(n7-1)*n10^2-(n7-1)*n10,2*ep^2+4*ep*(n7-1)+ep*n10-2*ep+2*(n7-1)^2+(n7-1)*n10-2*(n7-1))
      +`Y'(1,1,0,1,1+n5,1+n6,-1+(n7-1),1,1+n9,n10,0,0,0,0)*rat(2*ep*n5*n6*n9+n5*n6^2*n9+n5*n6*(n7-1)*n9+n5*n6*n9*n10-2*n5*n6*n9
      ,2*ep^2+4*ep*(n7-1)+ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Y'(1,1,0,1,1+n5,n6,(n7-1),1,n9,1+n10,0,0,0,0)*rat(-2*ep*n5*n9*n10+2*ep*n5*n10-n5*(n7-1)*n9*n10+n5*(n7-1)*n10-n5*n9*n10^2+n5
      *n9*n10+n5*n10^2-n5*n10,4*ep^2+8*ep*(n7-1)+2*ep*n10-4*ep+4*(n7-1)^2+2*(n7-1)*n10-4*(n7-1))
      +`Y'(1,1,0,1,n5,1+n6,-1+(n7-1),2,1+n9,n10,0,0,0,0)*rat(-2*ep*n6*n9-n6^2*n9-n6*(n7-1)*n9-n6*n9*n10+2*n6*n9,2*ep^2+4*ep*(n7-1)+
      ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Y'(1,1,0,1,n5,n6,(n7-1),2,n9,1+n10,0,0,0,0)*rat(2*ep*n9*n10-2*ep*n10+(n7-1)*n9*n10-(n7-1)*n10+n9*n10^2-n9*n10-n10^2+n10,4*
      ep^2+8*ep*(n7-1)+2*ep*n10-4*ep+4*(n7-1)^2+2*(n7-1)*n10-4*(n7-1))
      +`Y'(1,1,1,0,n5,1+n6,-1+(n7-1),2,1+n9,n10,0,0,0,0)*rat(2*ep*n6*n9+n6^2*n9+n6*(n7-1)*n9+n6*n9*n10-2*n6*n9,2*ep^2+4*ep*(n7-1)+
      ep*n10-4*ep+2*(n7-1)^2+(n7-1)*n10-4*(n7-1)-n10+2)
      +`Y'(1,1,1,0,n5,1+n6,(n7-1),1,1+n9,n10,0,0,0,0)*rat(n6*n9,1)
      +`Y'(1,1,1,0,n5,1+n6,(n7-1),1,n9,1+n10,0,0,0,0)*rat(n6*(n7-1)*n10-n6*n10,2*ep+2*(n7-1)+n10-2)
      +`Y'(1,1,1,0,n5,n6,(n7-1),2,n9,1+n10,0,0,0,0)*rat(-2*ep*n9*n10+2*ep*n10-(n7-1)*n9*n10+(n7-1)*n10-n9*n10^2+n9*n10+n10^2-n10,4*
      ep^2+8*ep*(n7-1)+2*ep*n10-4*ep+4*(n7-1)^2+2*(n7-1)*n10-4*(n7-1))
      +`Y'(2,0,1,1,n5,1+n6,-1+(n7-1),1,n9,1+n10,0,0,0,0)*rat(-n6*n10,2*ep+2*(n7-1)+n10-2)
      +`Y'(2,0,1,1,n5,n6,(n7-1),1,n9,1+n10,0,0,0,0)*rat(-(n7-1)*n10+n10,2*ep+2*(n7-1)+n10-2)
      +`Y'(2,1,1,0,n5,1+n6,-1+(n7-1),1,n9,1+n10,0,0,0,0)*rat(n6*n10,2*ep+2*(n7-1)+n10-2)
      );
*--#] n7 : 
#endif
#if ( `il' >= 14 )
*--#[ n5 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,1,1,n5?{>1},n6?`POS6',1,1,n9?`NEG9',n10?`NEG10',0,0,0,0)
			 = -rat(1,2*ep+n6+n10-1)*(
      +`Z'(1,1,1,1,(n5-1),1+n6,1,1,n9,n10,0,0,0,0)*rat(-2*ep*n6-(n5-1)*n6+n6,(n5-1))
      +`Z'(1,1,1,1,1+(n5-1),n6,1,1,n9,1+n10,0,0,0,0)*rat(-n10,1)
      +`Z'(1,1,1,1,(n5-1),n6,1,1,n9,1+n10,0,0,0,0)*rat(2*ep*n10+(n5-1)*n10-n10,(n5-1))
      +`Y'(0,2,1,1,1+(n5-1),n6,1,1,n9,n10,0,0,0,0)*rat(-1,1)
      +`Y'(1,1,1,1,1+(n5-1),1+n6,0,1,n9,n10,0,0,0,0)*rat(4*ep*n6+2*n6^2+2*n6*n10-2*n6,2*ep+n10)
      +`Y'(1,1,1,1,1+(n5-1),1+n6,1,0,n9,n10,0,0,0,0)*rat(-2*ep*n6-(n5-1)*n6+n6,ep)
      +`Y'(1,1,2,0,(n5-1),1+n6,1,1,n9,n10,0,0,0,0)*rat(n6,(n5-1))
      +`Y'(1,1,1,1,1+(n5-1),1+n6,0,1,n9,1+n10,0,0,0,0)*rat(-n6*n10,2*ep+n10)
      +`Y'(0,1,1,1,1+(n5-1),n6,1,1,n9,1+n10,0,0,0,0)*rat(-n10,1)
      +`Y'(1,0,1,1,1+(n5-1),1+n6,0,1,n9,1+n10,0,0,0,0)*rat(-n6*n10,2*ep+n10)
      +`Y'(1,1,0,1,1+(n5-1),n6,1,1,n9,1+n10,0,0,0,0)*rat(n10,1)
      +`Y'(1,1,0,1,(n5-1),n6,1,2,n9,1+n10,0,0,0,0)*rat(n10,(n5-1))
      +`Y'(1,1,1,0,1+(n5-1),1+n6,0,1,n9,1+n10,0,0,0,0)*rat(n6*n10,2*ep+n10)
      +`Y'(1,1,1,0,(n5-1),n6,1,2,n9,1+n10,0,0,0,0)*rat(-n10,(n5-1))
      +`Y'(1,1,1,1,1+(n5-1),n6,0,1,n9,1+n10,0,0,0,0)*rat(-2*ep*n10-n6*n10-n10^2+n10,2*ep+n10)
      );
*--#] n5 : 
#endif
#if ( `il' >= 15 )
*--#[ n10 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,1,1,1,n6?`POS6',1,1,n9?`NEG9',n10?neg_,0,0,0,0)
		 = -rat(2*ep+n6-2,
       20*ep^2+18*ep*n6+4*ep*n9+4*ep*(1+n10)-36*ep+4*n6^2+2*n6*n9+3*n6*(1+n10)-16*n6-4*n9-5*(1+n10)+16)*(
      +`Z'(1,1,1,1,1,n6,1,1,1+n9,-1+(1+n10),0,0,0,0)*rat(-2*n9,1)
      +`Z'(1,1,1,1,1,n6,1,1,n9,(1+n10),0,0,0,0)*rat(-8*ep^2-8*ep*n6-4*ep*(1+n10)+16*ep-2*n6^2-3*n6*(1+n10)+8*n6+5*(1+n10)-8,2*ep+n6-2)
      +`Z'(1,1,1,1,1,-1+n6,1,1,n9,(1+n10),0,0,0,0)*rat(4*ep+2*n6+(1+n10)-4,1)
      +`Y'(0,1,1,1,1,n6,2,1,-1+n9,(1+n10),0,0,0,0)*rat(-2*ep-2*n6+3,2*ep+n6-2)
      +`Y'(0,1,1,1,1,n6,2,1,n9,-1+(1+n10),0,0,0,0)*rat(n9,2*ep+n9)
      +`Y'(1,0,1,1,1,n6,2,1,-1+n9,(1+n10),0,0,0,0)*rat(2*ep+2*n6-3,2*ep+n6-2)
      +`Y'(1,0,1,1,1,n6,2,1,n9,-1+(1+n10),0,0,0,0)*rat(-n9,2*ep+n9)
      +`Y'(1,1,0,1,1,n6,1,2,1+n9,-1+(1+n10),-1,0,0,0)*rat(-n9,2*ep+n9)
      +`Y'(1,1,0,1,1,n6,1,2,n9,(1+n10),-1,0,0,0)*rat(2*ep+2*n6-3,2*ep+n6-2)
      +`Y'(1,1,0,1,2,n6,0,2,n9,(1+n10),-1,0,0,0)*rat(2*ep-1,4*ep^2)
      +`Y'(1,1,0,2,1,n6,1,1,n9,-1+(1+n10),0,0,0,0)*rat(1,1)
      +`Y'(1,1,1,0,1,n6,1,2,1+n9,-1+(1+n10),-1,0,0,0)*rat(n9,2*ep+n9)
      +`Y'(1,1,1,0,1,n6,1,2,n9,(1+n10),-1,0,0,0)*rat(-2*ep-2*n6+3,2*ep+n6-2)
      +`Y'(1,1,1,0,2,n6,0,2,n9,(1+n10),-1,0,0,0)*rat(-2*ep+1,4*ep^2)
      +`Y'(1,1,1,1,2,n6,0,1,-1+n9,(1+n10),0,0,0,0)*rat(-4*ep^2-2*ep*n9+4*ep+n9-1,4*ep^2)
      +`Y'(1,1,2,0,1,n6,1,1,n9,(1+n10),-1,0,0,0)*rat(-2*ep-2*n6+3,2*ep+n6-2)
      +`Y'(2,0,1,1,1,n6,1,1,-1+n9,(1+n10),0,0,0,0)*rat(2*ep+2*n6-3,2*ep+n6-2)
      +`Y'(2,0,1,1,1,n6,1,1,n9,(1+n10),-1,0,0,0)*rat(4*ep-2,2*ep+n9)
      +`Y'(0,1,1,1,1,n6,2,1,1+n9,-1+(1+n10),0,0,0,0)*rat(-n9,2*ep+n9)
      +`Y'(0,1,1,1,1,n6,2,1,n9,(1+n10),0,0,0,0)*rat(2*ep+2*n6-3,2*ep+n6-2)
      +`Y'(1,0,1,1,1,n6,2,1,1+n9,-1+(1+n10),0,0,0,0)*rat(n9,2*ep+n9)
      +`Y'(1,0,1,1,1,n6,2,1,n9,(1+n10),0,0,0,0)*rat(-2*ep-2*n6+3,2*ep+n6-2)
      +`Y'(1,1,0,1,1,n6,1,2,1+n9,-1+(1+n10),0,0,0,0)*rat(n9,2*ep+n9)
      +`Y'(1,1,0,1,1,n6,1,2,n9,(1+n10),0,0,0,0)*rat(-2*ep-2*n6+3,2*ep+n6-2)
      +`Y'(1,1,0,1,2,n6,0,2,n9,(1+n10),0,0,0,0)*rat(-2*ep+1,4*ep^2)
      +`Y'(1,1,0,2,1,n6,1,1,n9,(1+n10),0,0,0,0)*rat(-1,1)
      +`Y'(1,1,1,0,1,n6,1,2,1+n9,-1+(1+n10),0,0,0,0)*rat(-n9,2*ep+n9)
      +`Y'(1,1,1,0,1,n6,1,2,n9,(1+n10),0,0,0,0)*rat(2*ep+2*n6-3,2*ep+n6-2)
      +`Y'(1,1,1,0,2,n6,0,2,n9,(1+n10),0,0,0,0)*rat(2*ep-1,4*ep^2)
      +`Y'(1,1,1,1,2,n6,0,1,1+n9,(1+n10),-1,0,0,0)*rat(2*ep*n9-n9,4*ep^2+2*ep*n9)
      +`Y'(1,1,1,1,2,n6,0,1,n9,(1+n10),0,0,0,0)*rat(2*ep*n9-2*ep-n9+1,4*ep^2)
      +`Y'(2,0,1,1,1,n6,1,1,n9,(1+n10),0,0,0,0)*rat(-6*ep-n9+2,2*ep+n9)
      +`Y'(2,1,1,1,1,n6,0,1,n9,(1+n10),0,0,0,0)*rat(4*ep+n9-1,2*ep+n9)
      +`Y'(-1,1,1,1,1,n6,2,1,1+n9,-1+(1+n10),0,0,0,0)*rat(-n9,2*ep+n9)
      +`Y'(-1,1,1,1,1,n6,2,1,n9,(1+n10),0,0,0,0)*rat(2*ep+2*n6-3,2*ep+n6-2)
      +`Y'(0,0,1,1,1,n6,2,1,1+n9,-1+(1+n10),0,0,0,0)*rat(n9,2*ep+n9)
      +`Y'(0,0,1,1,1,n6,2,1,n9,(1+n10),0,0,0,0)*rat(-2*ep-2*n6+3,2*ep+n6-2)
      +`Y'(0,1,0,1,1,n6,1,2,1+n9,-1+(1+n10),0,0,0,0)*rat(-n9,2*ep+n9)
      +`Y'(0,1,0,1,1,n6,1,2,n9,(1+n10),0,0,0,0)*rat(2*ep+2*n6-3,2*ep+n6-2)
      +`Y'(0,1,0,1,1,n6,2,1,1+n9,-1+(1+n10),0,0,0,0)*rat(n9,2*ep+n9)
      +`Y'(0,1,0,1,1,n6,2,1,n9,(1+n10),0,0,0,0)*rat(-2*ep-2*n6+3,2*ep+n6-2)
      +`Y'(0,1,0,1,2,n6,0,2,n9,(1+n10),0,0,0,0)*rat(2*ep-1,4*ep^2)
      +`Y'(0,1,1,0,1,n6,1,2,1+n9,-1+(1+n10),0,0,0,0)*rat(n9,2*ep+n9)
      +`Y'(0,1,1,0,1,n6,1,2,n9,(1+n10),0,0,0,0)*rat(-2*ep-2*n6+3,2*ep+n6-2)
      +`Y'(0,1,1,0,2,n6,0,2,n9,(1+n10),0,0,0,0)*rat(-2*ep+1,4*ep^2)
      +`Y'(0,1,1,1,0,n6,2,1,n9,(1+n10),0,0,0,0)*rat(-2*ep+1,2*ep+n9)
      +`Y'(0,1,1,1,1,-1+n6,2,1,n9,(1+n10),0,0,0,0)*rat(-1,1)
      +`Y'(0,1,1,1,1,n6,1,1,1+n9,(1+n10),-1,0,0,0)*rat(-2*ep*n9+n9,2*ep+n9)
      +`Y'(0,1,1,1,1,n6,1,1,n9,(1+n10),0,0,0,0)*rat(-8*ep^3+4*ep^2*n6-4*ep^2*n9+12*ep^2+4*ep*n6*n9+2*ep*n6-12*ep+n6*n9^2-n6
      *n9-2*n6-n9^2+4,4*ep^2+2*ep*n6+2*ep*n9-4*ep+n6*n9-2*n9)
      +`Y'(0,1,1,1,2,n6,0,1,1+n9,(1+n10),-1,0,0,0)*rat(2*ep*n9-n9,4*ep^2+2*ep*n9)
      +`Y'(0,1,1,1,2,n6,0,1,n9,(1+n10),0,0,0,0)*rat(8*ep^2+2*ep*n9-6*ep-n9+1,4*ep^2)
      +`Y'(1,0,0,1,1,n6,2,1,1+n9,-1+(1+n10),0,0,0,0)*rat(-n9,2*ep+n9)
      +`Y'(1,0,0,1,1,n6,2,1,n9,(1+n10),0,0,0,0)*rat(2*ep+2*n6-3,2*ep+n6-2)
      +`Y'(1,0,0,2,1,n6,1,1,n9,(1+n10),0,0,0,0)*rat(-1,1)
      +`Y'(1,0,1,1,0,n6,2,1,n9,(1+n10),0,0,0,0)*rat(2*ep-1,2*ep+n9)
      +`Y'(1,0,1,1,1,-1+n6,2,1,n9,(1+n10),0,0,0,0)*rat(1,1)
      +`Y'(1,0,1,1,1,n6,1,1,1+n9,(1+n10),-1,0,0,0)*rat(2*ep*n9-n9,2*ep+n9)
      +`Y'(1,0,1,1,1,n6,1,1,n9,(1+n10),0,0,0,0)*rat(-28*ep^2-30*ep*n6-4*ep*n9-4*ep*(1+n10)+52*ep-6*n6^2-3*n6*n9-3*n6*(1+n10)+22*n6+
      5*n9+5*(1+n10)-20,2*ep+n6-2)
      +`Y'(1,1,-1,1,1,n6,1,2,1+n9,-1+(1+n10),0,0,0,0)*rat(n9,2*ep+n9)
      +`Y'(1,1,-1,1,1,n6,1,2,n9,(1+n10),0,0,0,0)*rat(-2*ep-2*n6+3,2*ep+n6-2)
      +`Y'(1,1,-1,1,2,n6,0,2,n9,(1+n10),0,0,0,0)*rat(-2*ep+1,4*ep^2)
      +`Y'(1,1,0,0,1,n6,1,2,1+n9,-1+(1+n10),0,0,0,0)*rat(-n9,2*ep+n9)
      +`Y'(1,1,0,0,1,n6,1,2,n9,(1+n10),0,0,0,0)*rat(2*ep+2*n6-3,2*ep+n6-2)
      +`Y'(1,1,0,0,2,n6,0,2,n9,(1+n10),0,0,0,0)*rat(2*ep-1,4*ep^2)
      +`Y'(1,1,0,1,1,n6,0,2,1+n9,(1+n10),-1,0,0,0)*rat(2*ep*n9-n9,4*ep^2+2*ep*n9)
      +`Y'(1,1,0,1,1,n6,1,1,n9,(1+n10),0,0,0,0)*rat(4*ep^2+2*ep*n6-4*ep+2*n6^2-n6*n9-4*n6+n9+2,2*ep+n6-2)
      +`Y'(1,1,0,1,2,n6,0,1,1+n9,(1+n10),-1,0,0,0)*rat(-2*ep*n9+n9,4*ep^2+2*ep*n9)
      +`Y'(1,1,0,1,2,n6,0,1,n9,(1+n10),0,0,0,0)*rat(-8*ep^2-2*ep*n9+6*ep+n9-1,4*ep^2)
      +`Y'(1,1,0,2,1,-1+n6,1,1,n9,(1+n10),0,0,0,0)*rat(1,1)
      +`Y'(1,1,1,0,1,n6,0,2,1+n9,(1+n10),-1,0,0,0)*rat(-2*ep*n9+n9,4*ep^2+2*ep*n9)
      +`Y'(1,1,1,0,1,n6,1,1,n9,(1+n10),0,0,0,0)*rat(12*ep^2+14*ep*n6+2*ep*n9+2*ep*(1+n10)-22*ep+2*n6^2+2*n6*n9+2*n6*(1+n10)-7*n6-3*
      n9-3*(1+n10)+6,2*ep+n6-2)
      +`Y'(1,1,1,1,1,n6,0,1,1+n9,(1+n10),-1,0,0,0)*rat(12*ep^2*n9+4*ep*n9^2-8*ep*n9-2*n9^2+n9,4*ep^2+2*ep*n9)
      +`Y'(1,1,1,1,1,n6,0,1,n9,(1+n10),0,0,0,0)*rat(8*ep^2+6*ep*n9-6*ep+n9^2-2*n9+1,2*ep+n9)
      +`Y'(2,0,0,1,1,n6,1,1,n9,(1+n10),0,0,0,0)*rat(-6*ep-n9+2,2*ep+n9)
      +`Y'(2,0,1,1,0,n6,1,1,n9,(1+n10),0,0,0,0)*rat(2*ep-1,2*ep+n9)
      +`Y'(2,1,0,1,1,n6,0,1,n9,(1+n10),0,0,0,0)*rat(4*ep+n9-1,2*ep+n9)
      +`Y'(0,1,1,1,1,n6,1,1,1+n9,(1+n10),0,0,0,0)*rat(-4*ep*n6*n9-4*ep*n9^2-4*n6*n9^2-2*n6*n9+6*n9^2+4*n9,4*ep^2+2*ep*n6+2*
      ep*n9-4*ep+n6*n9-2*n9)
      +`Y'(1,0,1,1,1,n6,1,1,1+n9,(1+n10),0,0,0,0)*rat(4*ep*n6*n9+4*ep*n9^2+4*n6*n9^2+2*n6*n9-6*n9^2-4*n9,4*ep^2+2*ep*n6+2*
      ep*n9-4*ep+n6*n9-2*n9)
      +`Y'(1,1,1,1,1,n6,0,1,2+n9,(1+n10),-1,0,0,0)*rat(-2*ep*n9^2-2*ep*n9+n9^2+n9,2*ep^2+ep*n9)
      );
*--#] n10 : 
#endif
#if ( `il' >= 16 )
*--#[ n9 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,1,1,1,n6?`POS6',1,1,n9?neg_,0,0,0,0,0)
			 = -rat(1,10*ep+2*n6+2*(1+n9)-6)*(
      +`Z'(1,1,1,1,1,n6,1,1,(1+n9),0,0,0,0,0)*rat(-4*ep-2*(1+n9)+2,1)
      +`Y'(0,1,1,1,1,n6,2,1,-1+(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(0,2,1,1,1,n6,1,1,(1+n9),0,0,-1,0,0)*rat(4*ep+2*n6-4,2*ep+n6-1)
      +`Y'(1,0,1,1,1,n6,2,1,-1+(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(1,1,0,1,1,1+n6,0,2,(1+n9),0,-1,0,0,0)*rat(-2*ep^2*n6-3*ep*n6^2+4*ep*n6-n6^3+3*n6^2-2*n6,4*ep^3+4*ep^2*n6-4*ep^2
      +ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,0,1,1,n6,1,2,(1+n9),0,-1,0,0,0)*rat(-1,1)
      +`Y'(1,1,1,0,1,1+n6,0,2,(1+n9),0,-1,0,0,0)*rat(2*ep^2*n6+3*ep*n6^2-4*ep*n6+n6^3-3*n6^2+2*n6,4*ep^3+4*ep^2*n6-4*ep^2+
      ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,1,0,1,n6,1,2,(1+n9),0,-1,0,0,0)*rat(1,1)
      +`Y'(1,1,1,1,1,1+n6,0,1,(1+n9),-1,0,0,0,0)*rat(-4*ep^3*n6-8*ep^2*n6^2+12*ep^2*n6-5*ep*n6^3+16*ep*n6^2-12*ep*n6-n6^4+
      5*n6^3-8*n6^2+4*n6,4*ep^3+4*ep^2*n6-4*ep^2+ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,2,0,1,1+n6,0,1,(1+n9),0,-1,0,0,0)*rat(2*ep^2*n6+3*ep*n6^2-4*ep*n6+n6^3-3*n6^2+2*n6,4*ep^3+4*ep^2*n6-4*ep^2+
      ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,2,0,1,n6,1,1,-1+(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(1,1,2,0,1,n6,1,1,(1+n9),0,-1,0,0,0)*rat(1,1)
      +`Y'(2,0,1,1,1,1+n6,0,1,-1+(1+n9),0,0,0,0,0)*rat(-2*ep^2*n6-3*ep*n6^2+4*ep*n6-n6^3+3*n6^2-2*n6,4*ep^3+4*ep^2*n6-4*
      ep^2+ep*n6^2-2*ep*n6+ep)
      +`Y'(0,1,1,1,1,n6,2,1,(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(0,2,1,1,1,n6,1,1,(1+n9),0,0,0,0,0)*rat(-6*ep-3*n6+5,2*ep+n6-1)
      +`Y'(1,0,1,1,1,n6,2,1,(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(1,1,0,1,1,1+n6,0,2,(1+n9),0,0,0,0,0)*rat(2*ep^2*n6+3*ep*n6^2-4*ep*n6+n6^3-3*n6^2+2*n6,4*ep^3+4*ep^2*n6-4*ep^2+ep
      *n6^2-2*ep*n6+ep)
      +`Y'(1,1,0,1,1,n6,1,2,(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(1,1,1,0,1,1+n6,0,2,(1+n9),0,0,0,0,0)*rat(-2*ep^2*n6-3*ep*n6^2+4*ep*n6-n6^3+3*n6^2-2*n6,4*ep^3+4*ep^2*n6-4*ep^2+
      ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,1,0,1,n6,1,2,(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(1,1,1,1,1,1+n6,0,1,(1+n9),0,0,0,0,0)*rat(-2*ep^2*n6-3*ep*n6^2+4*ep*n6-n6^3+3*n6^2-2*n6,4*ep^3+4*ep^2*n6-4*ep^2+
      ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,2,0,1,n6,1,1,(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(1,2,1,1,1,n6,0,1,(1+n9),0,0,0,0,0)*rat(4*ep+2*n6-3,2*ep+n6-1)
      +`Y'(-1,1,1,1,1,n6,2,1,(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(0,0,1,1,1,n6,2,1,(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(0,1,0,1,1,1+n6,0,2,(1+n9),0,0,0,0,0)*rat(-2*ep^2*n6-3*ep*n6^2+4*ep*n6-n6^3+3*n6^2-2*n6,4*ep^3+4*ep^2*n6-4*ep^2+
      ep*n6^2-2*ep*n6+ep)
      +`Y'(0,1,0,1,1,n6,1,2,(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(0,1,0,1,1,n6,2,1,(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(0,1,1,0,1,1+n6,0,2,(1+n9),0,0,0,0,0)*rat(2*ep^2*n6+3*ep*n6^2-4*ep*n6+n6^3-3*n6^2+2*n6,4*ep^3+4*ep^2*n6-4*ep^2+ep
      *n6^2-2*ep*n6+ep)
      +`Y'(0,1,1,0,1,n6,1,2,(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(0,1,1,1,0,n6,2,1,(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(0,1,1,1,1,-1+n6,2,1,(1+n9),0,0,0,0,0)*rat(2*ep+n6-2,2*ep+n6-1)
      +`Y'(0,1,1,1,1,1+n6,0,1,(1+n9),0,0,0,0,0)*rat(-8*ep^3*n6-12*ep^2*n6^2-2*ep^2*n6*(1+n9)+22*ep^2*n6-4*ep*n6^3-3*ep*n6^2*(1+n9)+
      21*ep*n6^2+4*ep*n6*(1+n9)-20*ep*n6-n6^3*(1+n9)+3*n6^3+3*n6^2*(1+n9)-9*n6^2-2*n6*(1+n9)+6*n6,4*ep^3+4*ep^2*n6-4*ep^2+ep*n6^2-2*ep*n6+ep)
      +`Y'(0,1,1,1,1,n6,1,1,(1+n9),0,0,0,0,0)*rat(-20*ep^2-10*ep*n6-4*ep*(1+n9)+18*ep-2*n6*(1+n9)+3*n6+2*(1+n9)-3,2*ep+n6-1)
      +`Y'(0,1,2,0,1,n6,1,1,(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(0,2,1,0,1,n6,1,1,(1+n9),0,0,0,0,0)*rat(-6*ep-3*n6+5,2*ep+n6-1)
      +`Y'(0,2,1,1,1,-1+n6,1,1,(1+n9),0,0,0,0,0)*rat(2*ep+n6-2,2*ep+n6-1)
      +`Y'(1,0,0,1,1,n6,2,1,(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(1,0,1,1,0,n6,2,1,(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(1,0,1,1,1,-1+n6,2,1,(1+n9),0,0,0,0,0)*rat(-2*ep-n6+2,2*ep+n6-1)
      +`Y'(1,0,1,1,1,1+n6,0,1,(1+n9),0,0,0,0,0)*rat(16*ep^3*n6+28*ep^2*n6^2+2*ep^2*n6*(1+n9)-40*ep^2*n6+14*ep*n6^3+3*ep*n6^2*(1+n9)
      -44*ep*n6^2-4*ep*n6*(1+n9)+32*ep*n6+2*n6^4+n6^3*(1+n9)-10*n6^3-3*n6^2*(1+n9)+16*n6^2+2*n6*(1+n9)-8*n6,4*ep^3+4*ep^2*n6-4*ep^2+ep*
      n6^2-2*ep*n6+ep)
      +`Y'(1,0,1,1,1,n6,1,1,(1+n9),0,0,0,0,0)*rat(-12*ep^2-14*ep*n6+24*ep-4*n6^2+14*n6-12,2*ep+n6-1)
      +`Y'(1,1,-1,1,1,1+n6,0,2,(1+n9),0,0,0,0,0)*rat(2*ep^2*n6+3*ep*n6^2-4*ep*n6+n6^3-3*n6^2+2*n6,4*ep^3+4*ep^2*n6-4*ep^2+
      ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,-1,1,1,n6,1,2,(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(1,1,0,0,1,1+n6,0,2,(1+n9),0,0,0,0,0)*rat(-2*ep^2*n6-3*ep*n6^2+4*ep*n6-n6^3+3*n6^2-2*n6,4*ep^3+4*ep^2*n6-4*ep^2+
      ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,0,0,1,n6,1,2,(1+n9),0,0,0,0,0)*rat(-1,1)
      +`Y'(1,1,0,1,1,1+n6,0,1,(1+n9),0,0,0,0,0)*rat(8*ep^3*n6+12*ep^2*n6^2+2*ep^2*n6*(1+n9)-22*ep^2*n6+4*ep*n6^3+3*ep*n6^2*(1+n9)-
      21*ep*n6^2-4*ep*n6*(1+n9)+20*ep*n6+n6^3*(1+n9)-3*n6^3-3*n6^2*(1+n9)+9*n6^2+2*n6*(1+n9)-6*n6,4*ep^3+4*ep^2*n6-4*ep^2+ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,0,1,1,n6,1,1,(1+n9),0,0,0,0,0)*rat(2*ep+(1+n9)-1,1)
      +`Y'(1,1,1,0,1,1+n6,0,1,(1+n9),0,0,0,0,0)*rat(-16*ep^3*n6-28*ep^2*n6^2-2*ep^2*n6*(1+n9)+40*ep^2*n6-14*ep*n6^3-3*ep*n6^2*
      (1+n9)+44*ep*n6^2+4*ep*n6*(1+n9)-32*ep*n6-2*n6^4-n6^3*(1+n9)+10*n6^3+3*n6^2*(1+n9)-16*n6^2-2*n6*(1+n9)+8*n6,4*ep^3+4*ep^2*n6-4*ep^2+
      ep*n6^2-2*ep*n6+ep)
      +`Y'(1,1,1,0,1,n6,1,1,(1+n9),0,0,0,0,0)*rat(-2*ep-(1+n9),1)
      +`Y'(1,1,1,1,0,n6,1,1,(1+n9),0,0,0,0,0)*rat(4*ep+(1+n9)-2,1)
      +`Y'(1,1,1,1,1,n6,0,1,(1+n9),0,0,0,0,0)*rat(8*ep^2+8*ep*n6-14*ep+2*n6^2-7*n6+6,2*ep+n6-1)
      +`Y'(1,1,2,0,0,n6,1,1,(1+n9),0,0,0,0,0)*rat(1,1)
      +`Y'(1,2,1,0,1,n6,0,1,(1+n9),0,0,0,0,0)*rat(4*ep+2*n6-3,2*ep+n6-1)
      +`Y'(0,1,1,1,1,1+n6,0,1,1+(1+n9),0,0,0,0,0)*rat(4*ep^2*n6*(1+n9)+6*ep*n6^2*(1+n9)-8*ep*n6*(1+n9)+2*n6^3*(1+n9)-6*n6^2*(1+n9)+4*n6*(1+n9),4*
      ep^3+4*ep^2*n6-4*ep^2+ep*n6^2-2*ep*n6+ep)
      +`Y'(1,0,1,1,1,1+n6,0,1,1+(1+n9),0,0,0,0,0)*rat(-4*ep^2*n6*(1+n9)-6*ep*n6^2*(1+n9)+8*ep*n6*(1+n9)-2*n6^3*(1+n9)+6*n6^2*(1+n9)-4*n6*(1+n9),4*
      ep^3+4*ep^2*n6-4*ep^2+ep*n6^2-2*ep*n6+ep)
      );
*--#] n9 : 
#endif
#if ( `il' >= 17 )
*--#[ n6 :
`SKIP6'
id,ifmatch->cross1,
	`Z'(1,1,1,1,1,n6?{>1},1,1,0,0,0,0,0,0)
			 = -rat(1,3*ep*(n6-1)+(n6-1)^2-(n6-1))*(
      +`Z'(1,1,1,1,1,(n6-1),1,1,0,0,0,0,0,0)*rat(-10*ep^2-7*ep*(n6-1)+9*ep-(n6-1)^2+3*(n6-1)-2,1)
      +`Y'(-1,1,1,1,1,1+(n6-1),2,1,0,0,0,0,0,0)*rat(-8*ep*(n6-1)-2*(n6-1)^2+3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(0,0,1,1,1,1+(n6-1),2,1,0,0,0,0,0,0)*rat(8*ep*(n6-1)+2*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(0,1,0,1,1,1+(n6-1),1,2,0,0,0,0,0,0)*rat(-8*ep*(n6-1)-2*(n6-1)^2+3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(0,1,0,1,1,1+(n6-1),2,1,0,0,0,0,0,0)*rat(8*ep*(n6-1)+2*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(0,1,0,1,2,1+(n6-1),0,2,0,0,0,0,0,0)*rat(2*ep*(n6-1)-(n6-1),8*ep^2)
      +`Y'(0,1,1,0,1,1+(n6-1),1,2,0,0,0,0,0,0)*rat(8*ep*(n6-1)+2*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(0,1,1,0,2,1+(n6-1),0,2,0,0,0,0,0,0)*rat(-2*ep*(n6-1)+(n6-1),8*ep^2)
      +`Y'(0,1,1,1,0,1+(n6-1),2,1,0,0,0,0,0,0)*rat(-2*ep*(n6-1)+(n6-1),4*ep)
      +`Y'(0,1,1,1,1,1+(n6-1),1,1,0,0,0,0,0,0)*rat(-44*ep^3*(n6-1)-14*ep^2*(n6-1)^2+44*ep^2*(n6-1)+9*ep*(n6-1)^2-13*ep*(n6-1)-(n6-1)^2+(n6-1),4*ep^2+2*ep
      *(n6-1)-2*ep)
      +`Y'(0,1,1,1,1,1+(n6-1),2,1,-1,0,0,0,0,0)*rat(8*ep*(n6-1)+2*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(0,1,1,1,1,1+(n6-1),2,1,0,0,0,0,0,0)*rat(-8*ep*(n6-1)-2*(n6-1)^2+3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(0,1,1,1,1,(n6-1),2,1,0,0,0,0,0,0)*rat(10*ep+3*(n6-1)-4,2)
      +`Y'(0,1,1,1,2,1+(n6-1),0,1,0,0,0,0,0,0)*rat(8*ep^2*(n6-1)-6*ep*(n6-1)+(n6-1),8*ep^2)
      +`Y'(1,0,0,1,1,1+(n6-1),2,1,0,0,0,0,0,0)*rat(-8*ep*(n6-1)-2*(n6-1)^2+3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,0,0,2,1,1+(n6-1),1,1,0,0,0,0,0,0)*rat(-(n6-1),2)
      +`Y'(1,0,1,1,0,1+(n6-1),2,1,0,0,0,0,0,0)*rat(2*ep*(n6-1)-(n6-1),4*ep)
      +`Y'(1,0,1,1,1,1+(n6-1),1,1,0,0,0,0,0,0)*rat(16*ep^2*(n6-1)+2*ep*(n6-1)^2-6*ep*(n6-1)-(n6-1)^3+(n6-1)^2,2*ep+(n6-1)-1)
      +`Y'(1,0,1,1,1,1+(n6-1),2,1,-1,0,0,0,0,0)*rat(-8*ep*(n6-1)-2*(n6-1)^2+3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,0,1,1,1,1+(n6-1),2,1,0,0,0,0,0,0)*rat(8*ep*(n6-1)+2*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,0,1,1,1,(n6-1),2,1,0,0,0,0,0,0)*rat(-10*ep-3*(n6-1)+4,2)
      +`Y'(1,1,-1,1,1,1+(n6-1),1,2,0,0,0,0,0,0)*rat(8*ep*(n6-1)+2*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,1,-1,1,2,1+(n6-1),0,2,0,0,0,0,0,0)*rat(-2*ep*(n6-1)+(n6-1),8*ep^2)
      +`Y'(1,1,0,0,1,1+(n6-1),1,2,0,0,0,0,0,0)*rat(-8*ep*(n6-1)-2*(n6-1)^2+3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,1,0,0,2,1+(n6-1),0,2,0,0,0,0,0,0)*rat(2*ep*(n6-1)-(n6-1),8*ep^2)
      +`Y'(1,1,0,1,1,1+(n6-1),1,1,0,0,0,0,0,0)*rat(22*ep^2*(n6-1)+9*ep*(n6-1)^2-19*ep*(n6-1)+(n6-1)^3-4*(n6-1)^2+4*(n6-1),2*ep+(n6-1)-1)
      +`Y'(1,1,0,1,1,1+(n6-1),1,2,0,0,-1,0,0,0)*rat(-8*ep*(n6-1)-2*(n6-1)^2+3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,1,0,1,1,1+(n6-1),1,2,0,0,0,0,0,0)*rat(8*ep*(n6-1)+2*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,1,0,1,2,1+(n6-1),0,1,0,0,0,0,0,0)*rat(-8*ep^2*(n6-1)+6*ep*(n6-1)-(n6-1),8*ep^2)
      +`Y'(1,1,0,1,2,1+(n6-1),0,2,0,0,-1,0,0,0)*rat(2*ep*(n6-1)-(n6-1),8*ep^2)
      +`Y'(1,1,0,1,2,1+(n6-1),0,2,0,0,0,0,0,0)*rat(-2*ep*(n6-1)+(n6-1),8*ep^2)
      +`Y'(1,1,0,2,1,1+(n6-1),1,1,0,-1,0,0,0,0)*rat((n6-1),2)
      +`Y'(1,1,0,2,1,1+(n6-1),1,1,0,0,0,0,0,0)*rat(-(n6-1),2)
      +`Y'(1,1,0,2,1,(n6-1),1,1,0,0,0,0,0,0)*rat((n6-1),2)
      +`Y'(1,1,1,0,1,1+(n6-1),1,1,0,0,0,0,0,0)*rat(-48*ep^2*(n6-1)-20*ep*(n6-1)^2+26*ep*(n6-1)-2*(n6-1)^3+5*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,1,1,0,1,1+(n6-1),1,2,0,0,-1,0,0,0)*rat(8*ep*(n6-1)+2*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,1,1,0,1,1+(n6-1),1,2,0,0,0,0,0,0)*rat(-8*ep*(n6-1)-2*(n6-1)^2+3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(1,1,1,0,2,1+(n6-1),0,2,0,0,-1,0,0,0)*rat(-2*ep*(n6-1)+(n6-1),8*ep^2)
      +`Y'(1,1,1,0,2,1+(n6-1),0,2,0,0,0,0,0,0)*rat(2*ep*(n6-1)-(n6-1),8*ep^2)
      +`Y'(1,1,1,1,1,1+(n6-1),0,1,0,0,0,0,0,0)*rat(8*ep^2*(n6-1)-6*ep*(n6-1)+(n6-1),4*ep)
      +`Y'(1,1,1,1,2,1+(n6-1),0,1,-1,0,0,0,0,0)*rat(-4*ep^2*(n6-1)+4*ep*(n6-1)-(n6-1),8*ep^2)
      +`Y'(1,1,1,1,2,1+(n6-1),0,1,0,0,0,0,0,0)*rat(-2*ep*(n6-1)+(n6-1),8*ep^2)
      +`Y'(1,1,2,0,1,1+(n6-1),1,1,0,0,-1,0,0,0)*rat(8*ep*(n6-1)+2*(n6-1)^2-3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(2,0,0,1,1,1+(n6-1),1,1,0,0,0,0,0,0)*rat(-3*ep*(n6-1)+(n6-1),2*ep)
      +`Y'(2,0,1,1,0,1+(n6-1),1,1,0,0,0,0,0,0)*rat(2*ep*(n6-1)-(n6-1),4*ep)
      +`Y'(2,0,1,1,1,1+(n6-1),1,1,-1,0,0,0,0,0)*rat(-8*ep*(n6-1)-2*(n6-1)^2+3*(n6-1),4*ep+2*(n6-1)-2)
      +`Y'(2,0,1,1,1,1+(n6-1),1,1,0,0,-1,0,0,0)*rat(2*ep*(n6-1)-(n6-1),2*ep)
      +`Y'(2,0,1,1,1,1+(n6-1),1,1,0,0,0,0,0,0)*rat(-3*ep*(n6-1)+(n6-1),2*ep)
      +`Y'(2,1,0,1,1,1+(n6-1),0,1,0,0,0,0,0,0)*rat(4*ep*(n6-1)-(n6-1),4*ep)
      +`Y'(2,1,1,1,1,1+(n6-1),0,1,0,0,0,0,0,0)*rat(4*ep*(n6-1)-(n6-1),4*ep)
      );
*--#] n6 : 
#endif
goto cross2;
Label cross1;
$ired = 1;
Label cross2;
id	`Z'(?a,x?neg0_,?b,n9?,...,n14?) = `Y'(?a,x,?b,n9,...,n14);
*id	`Y'(?a) = 0;
ModuleOption,maximum,$ired;
*Print +f +s;
.sort:Main cross loop, pass `$pass';
#if ( `$ired' > 0 )
	#redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
id	`Z'(1,1,1,1,1,1,1,1,0,0,0,0,0,0) = Master(cross);
*
*--#] Loop : 
#break
*--#] cross : 
*--#[ bebe :
#case bebe
*
*--#[ Loop :
#do il = 1,16
#$pass = 0;

#do ired = 1,1

#if 0

#define SORTCRITA "1"
#define SORTCRITB "1"


#if ( `il' == 4)
* third index should be lowest
	#redefine SORTCRITA "n1"
	#redefine SORTCRITB "n3"
#endif

#if ( `il' == 10)
* third index should be lowest
	#redefine SORTCRITA "n2"
	#redefine SORTCRITB "n5"
#endif

#if ( `il' == 11)
* third index should be lowest
	#redefine SORTCRITA "n8"
	#redefine SORTCRITB "n7"
#endif

#if ( `il' == 12)
* third index should be lowest
	#redefine SORTCRITA "n7"
	#redefine SORTCRITB "n8"
#endif

#if ( `il' == 15)
* third index should be lowest
	#redefine SORTCRITA "n5"
	#redefine SORTCRITB "n2"
#endif

* Apply automorphisms
ifmatch->Lbebeauto `Z'(n1?,...,n8?,n9?{0,},n10?{0,},n11?{0,},n12?{0,},n13?{0,},n14?{0,}) = 
   `Z'(`SORTCRITA',n3,n5,n1,n4,n2,n6,n8,n7,n9,n10,n11,n12,n13,n14,sign_(n11+n13))
 * `Z'(`SORTCRITB',n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,1)
 ;
ifmatch->Lbebeauto `Z'(n1?,...,n8?,n9?{0,},n10?,n11?{0,},n12?{0,},n13?{0,},n14?{0,}) = 
   `Z'(`SORTCRITA',n3,n5,n1,n4,n2,n6,n8,n7,n9,n10,n11,n12,n13,n14,sign_(n11+n13))
 * `Z'(`SORTCRITB',n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,1)
 ;
ifmatch->Lbebeauto `Z'(n1?,...,n8?,n9?{0,},n10?{0,},n11?,n12?{0,},n13?{0,},n14?{0,}) = 
   `Z'(`SORTCRITA',n3,n5,n1,n4,n2,n6,n8,n7,n9,n10,n11,n12,n13,n14,sign_(n11+n13))
 * `Z'(`SORTCRITB',n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,1)
 ;
ifmatch->Lbebeauto `Z'(n1?,...,n8?,n9?{0,},n10?,n11?,n12?{0,},n13?{0,},n14?{0,}) = 
   `Z'(`SORTCRITA',n3,n5,n1,n4,n2,n6,n8,n7,n9,n10,n11,n12,n13,n14,sign_(n11+n13))
 * `Z'(`SORTCRITB',n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,1)
 ;
goto Lbebeautoend;
label Lbebeauto;
repeat id `Z'(?a) * `Z'(?b) = `Z'(?a);
id `Z'(n1?,?a,n?) = n*`Z'(n1,?a);
label Lbebeautoend;
.sort:`Z' automorphism;
#endif

#$ired = 0;
#$pass = $pass+1;
#ifdef `RULEINFO'
#write "Entering bebe loop `il', pass `$pass' at `time_' sec."
#endif
#if ( `il' >= 1 )
*--#[ stat file :
*
*  Prework to make the reduction of the bebe topology faster.
*
*   B20= (1)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?neg0_,n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_)
			 = -rat(-1,(n8-1))*(
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+(n8-1),n9,n10,n11,(1+n12),n13,n14)*rat((n8-1),1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,1+(n8-1),n9,n10,n11,(1+n12),n13,n14)*rat((n8-1),1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,1+n7,(n8-1),n9,n10,n11,(1+n12),n13,n14)*rat(n7,1)
      +`Z'(n1,n2,n3,n4,n5,-1+n6,1+n7,(n8-1),n9,n10,n11,(1+n12),n13,n14)*rat(-n7,1)
      +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+(n8-1),n9,n10,n11,(1+n12),n13,n14)*rat(-(n8-1),1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,n11,(1+n12),n13,n14)*rat(-2*ep-2*n6-n7-(n8-1)-n11-n13+4,1)
      );
 
*   B10= (1)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_)
			 = -rat(-1,(n2-1))*(
      +`Z'(n1,(n2-1),n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,(1+n13),n14)*rat(n5,1)
      +`Z'(n1,(n2-1),n3,n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,(1+n13),n14)*rat(-n5,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,(1+n13),n14)*rat(-n9,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,-1+(1+n13),n14)*rat(-n12,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,(1+n13),n14)*rat(2*ep+n5+2*n6+n11+n12+(1+n13)-4,1)
      +`Z'(-1+n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,(1+n13),n14)*rat(-n12,1)
      +`Z'(n1,-1+(n2-1),n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,(1+n13),n14)*rat(-n12,1)
      +`Z'(n1,(n2-1),n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,(1+n13),n14)*rat(-n12,1)
      +`Z'(n1,(n2-1),n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+(1+n13),n14)*rat(-2*(1+n13),1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,(1+n13),n14)*rat(n12,1)
      );

*   B11= (1)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_)
			 = -rat(1,(n4-1))*(
      +`Z'(n1,n2,n3,1+(n4-1),n5,n6,n7,n8,(1+n9),n10,n11,n12,n13,n14)*rat(-(n4-1),1)
      +`Z'(n1,-1+n2,n3,1+(n4-1),n5,n6,n7,n8,(1+n9),n10,n11,n12,n13,n14)*rat(-(n4-1),1)
      +`Z'(n1,n2,n3,1+(n4-1),-1+n5,n6,n7,n8,(1+n9),n10,n11,n12,n13,n14)*rat((n4-1),1)
      +`Z'(n1,n2,n3,(n4-1),-1+n5,n6,1+n7,n8,(1+n9),n10,n11,n12,n13,n14)*rat(n7,1)
      +`Z'(n1,n2,n3,(n4-1),n5,-1+n6,1+n7,n8,(1+n9),n10,n11,n12,n13,n14)*rat(-n7,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,-1+(1+n9),1+n10,n11,n12,n13,n14)*rat(n10,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,(1+n9),n10,n11,-1+n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,(1+n9),n10,n11,n12,n13,n14)*rat(2*ep+(n4-1)+2*n5+n7+n10+n14-4,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,(1+n9),1+n10,n11,n12,n13,n14)*rat(-2*n10,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,(1+n9),n10,n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(-1+n1,n2,n3,(n4-1),n5,n6,n7,n8,(1+n9),n10,n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,-1+n3,(n4-1),n5,n6,n7,n8,(1+n9),n10,n11,n12,n13,1+n14)*rat(n14,1)
      );
 
*   B9= (1)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_)
			 = -rat(-1,(n2-1))*(
      +`Z'(n1,1+(n2-1),n3,n4,n5,n6,n7,n8,n9,(1+n10),n11,n12,n13,n14)*rat((n2-1),1)
      +`Z'(n1,1+(n2-1),n3,-1+n4,n5,n6,n7,n8,n9,(1+n10),n11,n12,n13,n14)*rat((n2-1),1)
      +`Z'(n1,1+(n2-1),n3,n4,-1+n5,n6,n7,n8,n9,(1+n10),n11,n12,n13,n14)*rat(-(n2-1),1)
      +`Z'(n1,(n2-1),n3,n4,-1+n5,1+n6,n7,n8,n9,(1+n10),n11,n12,n13,n14)*rat(-n6,1)
      +`Z'(n1,(n2-1),n3,n4,n5,1+n6,-1+n7,n8,n9,(1+n10),n11,n12,n13,n14)*rat(n6,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,-1+n9,(1+n10),1+n11,n12,n13,n14)*rat(n11,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,1+n9,-1+(1+n10),n11,n12,n13,n14)*rat(-n9,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,-1+(1+n10),1+n11,n12,n13,n14)*rat(n11,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,-1+(1+n10),n11,n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,(1+n10),n11,1+n12,n13,-1+n14)*rat(-n12,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,(1+n10),n11,n12,n13,n14)*rat(-2*ep-(n2-1)-2*n5-n6-n9-n12+4,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,1+n9,(1+n10),n11,n12,n13,n14)*rat(2*n9,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,(1+n10),1+n11,n12,n13,n14)*rat(-2*n11,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,(1+n10),n11,1+n12,n13,n14)*rat(n12,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,(1+n10),n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'(-1+n1,(n2-1),n3,n4,n5,n6,n7,n8,n9,(1+n10),n11,1+n12,n13,n14)*rat(n12,1)
      +`Z'(n1,-1+(n2-1),n3,n4,n5,n6,n7,n8,n9,(1+n10),n11,n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,(n2-1),-1+n3,n4,n5,n6,n7,n8,n9,(1+n10),n11,1+n12,n13,n14)*rat(-n12,1)
      +`Z'(n1,(n2-1),n3,-1+n4,n5,n6,n7,n8,n9,(1+n10),n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'(n1,(n2-1),n3,n4,-1+n5,n6,n7,n8,n9,(1+n10),n11,n12,1+n13,n14)*rat(2*n13,1)
      +`Z'(n1,(n2-1),n3,n4,n5,-1+n6,n7,n8,n9,(1+n10),n11,n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,(n2-1),n3,n4,n5,n6,-1+n7,n8,n9,(1+n10),n11,n12,1+n13,n14)*rat(-n13,1)
      );
*   B2=  (2)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1}
			,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_)
			 = -rat(1,(n8-1))*(
      +`Z'(n1,n2,n3,1+n4,n5,n6,n7,(n8-1),n9,-1+n10,n11,n12,(1+n13),n14)*rat(-n4,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+(n8-1),n9,n10,n11,n12,(1+n13),n14)*rat(-(n8-1),1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,1+(n8-1),n9,n10,n11,n12,(1+n13),n14)*rat((n8-1),1)
      +`Z'(n1,n2,n3,n4,n5,-1+n6,n7,1+(n8-1),n9,n10,n11,n12,(1+n13),n14)*rat((n8-1),1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,n11,n12,(1+n13),n14)*rat(6*ep+2*n2+2*n4+2*n5+2*n6+2*n7+(n8-1)+n9+n10+n11+n12+2*(1+n13)+n14-12,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,1+n10,n11,n12,(1+n13),n14)*rat(-2*n10,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,n11,n12,(1+n13),1+n14)*rat(-n14,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,n11,n12,(1+n13),1+n14)*rat(-n14,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,(n8-1),n9,n10,n11,n12,(1+n13),1+n14)*rat(n14,1)
      );
*   B19= (2)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_
        ,n8?{>1},n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_)
			 = -rat(1,(n8-1))*(
      +`Z'(n1,n2,n3,n4,n5,n6,n7,1+(n8-1),n9,-1+n10,n11,n12,n13,(1+n14))*rat(-(n8-1),1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+(n8-1),n9,n10,n11,n12,n13,(1+n14))*rat((n8-1),1)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,1+(n8-1),n9,n10,n11,n12,n13,(1+n14))*rat((n8-1),1)
      +`Z'(n1,n2,n3,n4,-1+n5,1+n6,n7,(n8-1),n9,n10,n11,n12,n13,(1+n14))*rat(n6,1)
      +`Z'(n1,n2,n3,n4,n5,1+n6,-1+n7,(n8-1),n9,n10,n11,n12,n13,(1+n14))*rat(-n6,1)
      +`Z'(n1,n2,n3,n4,n5,n6,-1+n7,1+(n8-1),n9,n10,n11,n12,n13,(1+n14))*rat(-(n8-1),1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),-1+n9,n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,-1+n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,-1+n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,n11,n12,n13,(1+n14))*rat(-2*ep-n6-2*n7-(n8-1)-n11-n13+4,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,1+n11,n12,n13,(1+n14))*rat(2*n11,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,n7,(n8-1),n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,n2,n3,n4,-1+n5,n6,n7,(n8-1),n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      );
 
*   B18= (2)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_)
			 = -rat(-1,(n6-1))*(
      +`Z'(n1,n2,n3,n4,1+n5,(n6-1),n7,n8,n9,-1+n10,n11,n12,(1+n13),n14)*rat(n5,1)
      +`Z'(n1,n2,n3,n4,1+n5,(n6-1),n7,n8,n9,n10,n11,n12,(1+n13),n14)*rat(-n5,1)
      +`Z'(n1,-1+n2,n3,n4,1+n5,(n6-1),n7,n8,n9,n10,n11,n12,(1+n13),n14)*rat(n5,1)
      +`Z'(n1,n2,n3,-1+n4,1+n5,(n6-1),n7,n8,n9,n10,n11,n12,(1+n13),n14)*rat(-n5,1)
      +`Z'(n1,n2,n3,n4,n5,(n6-1),n7,n8,-1+n9,n10,1+n11,n12,(1+n13),n14)*rat(-n11,1)
      +`Z'(n1,n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,n11,n12,(1+n13),n14)*rat(2*ep+2*n2+n5+n9+n12+(1+n13)-4,1)
      +`Z'(n1,-1+n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,n11,n12,1+(1+n13),n14)*rat(-2*(1+n13),1)
      );
 
*   B15= (2)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_)
			 = -rat(-1,(n4-1))*(
      +`Z'(n1,n2,n3,1+(n4-1),n5,n6,n7,n8,n9,n10,(1+n11),n12,-1+n13,n14)*rat((n4-1),1)
      +`Z'(n1,n2,n3,1+(n4-1),-1+n5,n6,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat((n4-1),1)
      +`Z'(n1,n2,n3,1+(n4-1),n5,-1+n6,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat((n4-1),1)
      +`Z'(n1,n2,n3,1+(n4-1),n5,n6,-1+n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(-(n4-1),1)
      +`Z'(n1,n2,n3,(n4-1),-1+n5,n6,1+n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(n7,1)
      +`Z'(n1,n2,n3,(n4-1),1+n5,-1+n6,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(n5,1)
      +`Z'(n1,n2,n3,(n4-1),1+n5,n6,-1+n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(-n5,1)
      +`Z'(n1,n2,n3,(n4-1),n5,-1+n6,1+n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(-n7,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,1+n10,-1+(1+n11),n12,n13,n14)*rat(-n10,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,(1+n11),-1+n12,n13,1+n14)*rat(n14,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,(1+n11),n12,-1+n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(n5-n7,1)
      +`Z'(-1+n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,(1+n11),n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,-1+n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,(1+n11),n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,n3,(n4-1),n5,-1+n6,n7,n8,n9,n10,(1+n11),n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,-1+n8,n9,n10,(1+n11),n12,n13,1+n14)*rat(n14,1)
      );
 
*   B17= (3)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_)
			 = -rat(-1,(n7-1))*(
      +`Z'(n1,n2,n3,n4,n5,1+n6,(n7-1),n8,n9,n10,n11,-1+n12,n13,(1+n14))*rat(n6,1)
      +`Z'(n1,n2,n3,n4,n5,n6,1+(n7-1),n8,n9,-1+n10,n11,n12,n13,(1+n14))*rat((n7-1),1)
      +`Z'(-1+n1,n2,n3,n4,n5,1+n6,(n7-1),n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n6,1)
      +`Z'(n1,-1+n2,n3,n4,n5,1+n6,(n7-1),n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n6,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,1+(n7-1),n8,n9,n10,n11,n12,n13,(1+n14))*rat(-(n7-1),1)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,1+(n7-1),n8,n9,n10,n11,n12,n13,(1+n14))*rat(-(n7-1),1)
      +`Z'(n1,n2,n3,n4,n5,1+n6,(n7-1),-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(n6,1)
      +`Z'(n1,n2,n3,n4,n5,n6,1+(n7-1),-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat((n7-1),1)
      +`Z'(n1,n2,n3,n4,n5,n6,(n7-1),n8,-1+n9,n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Z'(n1,n2,n3,n4,n5,n6,(n7-1),n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,n2,n3,n4,n5,n6,(n7-1),n8,n9,n10,n11,n12,n13,(1+n14))*rat(2*ep+n6+(n7-1)+2*n8+n11+n13-4,1)
      +`Z'(n1,n2,n3,n4,n5,n6,(n7-1),n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,(n7-1),n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,(n7-1),n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(2*n13,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,(n7-1),n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      );
 
*   B12= (3)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_)
			 = -rat(-1,(n4-1))*(
      +`Z'(n1,1+n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,-1+n12,n13,(1+n14))*rat(n2,1)
      +`Z'(n1,n2,n3,1+(n4-1),n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,(1+n14))*rat((n4-1),1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-6*ep-2*n2-2*(n4-1)-2*n5-2*n6-2*n7-2*n8-n9-n10-n11-n12-n13-
      (1+n14)+12,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,1+n9,n10,n11,n12,n13,(1+n14))*rat(n9,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,n12,n13,1+(1+n14))*rat((1+n14),1)
      +`Z'(-1+n1,n2,n3,(n4-1),n5,n6,n7,n8,1+n9,n10,n11,n12,n13,(1+n14))*rat(n9,1)
      +`Z'(-1+n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(-n10,1)
      +`Z'(-1+n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(2*n12,1)
      +`Z'(-1+n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(-1+n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Z'(n1,-1+n2,n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,n2,-1+n3,(n4-1),n5,n6,n7,n8,1+n9,n10,n11,n12,n13,(1+n14))*rat(-n9,1)
      +`Z'(n1,n2,-1+n3,(n4-1),n5,n6,n7,n8,n9,1+n10,n11,n12,n13,(1+n14))*rat(n10,1)
      +`Z'(n1,n2,-1+n3,(n4-1),n5,n6,n7,n8,n9,n10,n11,n12,n13,1+(1+n14))*rat(-(1+n14),1)
      +`Z'(n1,n2,n3,(n4-1),n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      );
 
*   B16= (4)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>1},n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_)
			 = -rat(-1,(n6-1))*(
      +`Z'(n1,1+n2,n3,n4,n5,(n6-1),n7,n8,-1+n9,n10,(1+n11),n12,n13,n14)*rat(n2,1)
      +`Z'(n1,n2,n3,n4,1+n5,(n6-1),n7,n8,-1+n9,n10,(1+n11),n12,n13,n14)*rat(n5,1)
      +`Z'(n1,n2,n3,n4,1+n5,(n6-1),n7,n8,n9,-1+n10,(1+n11),n12,n13,n14)*rat(n5,1)
      +`Z'(n1,n2,n3,n4,1+n5,(n6-1),n7,n8,n9,n10,(1+n11),n12,n13,n14)*rat(-2*n5,1)
      +`Z'(n1,n2,n3,n4,n5,(n6-1),n7,n8,-1+n9,n10,(1+n11),n12,1+n13,n14)*rat(-n13,1)
      +`Z'(n1,n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,-1+(1+n11),n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,n2,n3,n4,n5,(n6-1),n7,n8,1+n9,n10,(1+n11),n12,n13,n14)*rat(2*n9,1)
      +`Z'(n1,n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,1+(1+n11),n12,n13,n14)*rat(-2*(1+n11),1)
      +`Z'(n1,n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,(1+n11),1+n12,n13,n14)*rat(n12,1)
      +`Z'(-1+n1,n2,n3,n4,n5,(n6-1),n7,n8,n9,n10,(1+n11),1+n12,n13,n14)*rat(n12,1)
      +`Z'(n1,n2,-1+n3,n4,n5,(n6-1),n7,n8,n9,n10,(1+n11),1+n12,n13,n14)*rat(-n12,1)
      );

*   B1= (4)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?neg0_,n10?neg0_,n11?neg_,n12?neg0_,n13?neg0_,n14?neg0_)
			 = -rat(1,(n8-1))*(
      +`Z'(n1,1+n2,n3,n4,n5,n6,n7,(n8-1),-1+n9,n10,(1+n11),n12,n13,n14)*rat(n2,1)
      +`Z'(n1,n2,n3,1+n4,n5,n6,n7,(n8-1),n9,-1+n10,(1+n11),n12,n13,n14)*rat(-n4,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,1+(n8-1),-1+n9,n10,(1+n11),n12,n13,n14)*rat((n8-1),1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,1+(n8-1),n9,n10,(1+n11),n12,n13,n14)*rat(-(n8-1),1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,1+(n8-1),n9,n10,(1+n11),n12,n13,n14)*rat(-(n8-1),1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,1+(n8-1),n9,n10,(1+n11),n12,n13,n14)*rat((n8-1),1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,-1+(1+n11),n12,1+n13,n14)*rat(n13,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),1+n9,n10,(1+n11),n12,n13,n14)*rat(2*n9,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,1+n10,(1+n11),n12,n13,n14)*rat(-2*n10,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,(1+n11),1+n12,n13,n14)*rat(n12,1)
      +`Z'(n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,(1+n11),n12,n13,1+n14)*rat(-n14,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,(1+n11),1+n12,n13,n14)*rat(n12,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,n7,(n8-1),n9,n10,(1+n11),n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,(n8-1),n9,n10,(1+n11),1+n12,n13,n14)*rat(-n12,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,n7,(n8-1),n9,n10,(1+n11),n12,n13,1+n14)*rat(n14,1)
      );
 
*   B8= (4)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>1},n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_)
			 = -rat(-1,(n7-1))*(
      +`Z'(n1,n2,n3,1+n4,n5,n6,(n7-1),n8,n9,-1+n10,n11,n12,(1+n13),n14)*rat(-n4,1)
      +`Z'(n1,n2,n3,n4,1+n5,n6,(n7-1),n8,n9,-1+n10,n11,n12,(1+n13),n14)*rat(-n5,1)
      +`Z'(n1,n2,n3,n4,n5,n6,1+(n7-1),n8,n9,-1+n10,n11,n12,(1+n13),n14)*rat(-(n7-1),1)
      +`Z'(n1,n2,n3,n4,1+n5,n6,(n7-1),n8,n9,n10,n11,n12,(1+n13),n14)*rat(n5,1)
      +`Z'(n1,n2,n3,n4,n5,n6,1+(n7-1),n8,n9,n10,n11,n12,(1+n13),n14)*rat((n7-1),1)
      +`Z'(n1,-1+n2,n3,n4,1+n5,n6,(n7-1),n8,n9,n10,n11,n12,(1+n13),n14)*rat(-n5,1)
      +`Z'(n1,-1+n2,n3,n4,n5,n6,1+(n7-1),n8,n9,n10,n11,n12,(1+n13),n14)*rat(-(n7-1),1)
      +`Z'(n1,n2,n3,-1+n4,1+n5,n6,(n7-1),n8,n9,n10,n11,n12,(1+n13),n14)*rat(n5,1)
      +`Z'(n1,n2,n3,-1+n4,n5,n6,1+(n7-1),n8,n9,n10,n11,n12,(1+n13),n14)*rat((n7-1),1)
      +`Z'(n1,n2,n3,n4,n5,-1+n6,1+(n7-1),n8,n9,n10,n11,n12,(1+n13),n14)*rat(-(n7-1),1)
      +`Z'(n1,n2,n3,n4,n5,n6,(n7-1),n8,n9,n10,n11,n12,(1+n13),n14)*rat(2*ep+2*n4+n5+(n7-1)+n10+n14-4,1)
      +`Z'(n1,n2,n3,n4,n5,n6,(n7-1),n8,n9,1+n10,n11,n12,(1+n13),n14)*rat(-2*n10,1)
      +`Z'(n1,n2,n3,n4,n5,n6,(n7-1),n8,n9,n10,n11,n12,(1+n13),1+n14)*rat(-n14,1)
      +`Z'(-1+n1,n2,n3,n4,n5,n6,(n7-1),n8,n9,n10,n11,n12,(1+n13),1+n14)*rat(-n14,1)
      +`Z'(n1,n2,-1+n3,n4,n5,n6,(n7-1),n8,n9,n10,n11,n12,(1+n13),1+n14)*rat(n14,1)
      );

*   B7= (5)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>1},n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg_)
			 = -rat(-1,(n5-1))*(
      +`Z'(n1,1+n2,n3,n4,(n5-1),n6,n7,n8,n9,n10,n11,-1+n12,n13,(1+n14))*rat(-n2,1)
      +`Z'(n1,n2,n3,n4,1+(n5-1),n6,n7,n8,n9,-1+n10,n11,n12,n13,(1+n14))*rat((n5-1),1)
      +`Z'(n1,n2,n3,n4,1+(n5-1),n6,n7,n8,n9,n10,n11,-1+n12,n13,(1+n14))*rat(-(n5-1),1)
      +`Z'(n1,n2,n3,n4,(n5-1),1+n6,n7,n8,n9,n10,n11,-1+n12,n13,(1+n14))*rat(-n6,1)
      +`Z'(-1+n1,n2,n3,n4,1+(n5-1),n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat((n5-1),1)
      +`Z'(-1+n1,n2,n3,n4,(n5-1),1+n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n6,1)
      +`Z'(n1,-1+n2,n3,n4,1+(n5-1),n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat((n5-1),1)
      +`Z'(n1,-1+n2,n3,n4,(n5-1),1+n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n6,1)
      +`Z'(n1,n2,-1+n3,n4,1+(n5-1),n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-(n5-1),1)
      +`Z'(n1,n2,n3,-1+n4,1+(n5-1),n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-(n5-1),1)
      +`Z'(n1,n2,n3,n4,(n5-1),1+n6,n7,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n6,1)
      +`Z'(n1,n2,n3,n4,(n5-1),n6,n7,n8,-1+n9,n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(n1,n2,n3,n4,(n5-1),n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(2*n13,1)
      +`Z'(n1,n2,n3,n4,(n5-1),n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(2*ep+2*n2+(n5-1)+n6+n9+n12-4,1)
      +`Z'(n1,n2,n3,n4,(n5-1),n6,n7,n8,1+n9,n10,n11,n12,n13,(1+n14))*rat(-n9,1)
      +`Z'(n1,n2,n3,n4,(n5-1),n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Z'(-1+n1,n2,n3,n4,(n5-1),n6,n7,n8,1+n9,n10,n11,n12,n13,(1+n14))*rat(-n9,1)
      +`Z'(-1+n1,n2,n3,n4,(n5-1),n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Z'(-1+n1,n2,n3,n4,(n5-1),n6,n7,n8,n9,n10,n11,1+n12,n13,(1+n14))*rat(-2*n12,1)
      +`Z'(-1+n1,n2,n3,n4,(n5-1),n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,-1+n2,n3,n4,(n5-1),n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-3*n13,1)
      +`Z'(n1,n2,-1+n3,n4,(n5-1),n6,n7,n8,1+n9,n10,n11,n12,n13,(1+n14))*rat(n9,1)
      +`Z'(n1,n2,-1+n3,n4,(n5-1),n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(n1,n2,n3,n4,(n5-1),-1+n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(n1,n2,n3,n4,(n5-1),n6,n7,-1+n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      );
#if(1==0)
*   B14= (6)
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_)
			 = -rat(-1,(n4-1))*(
      +`Z'(n1,n2,n3,(n4-1),1+n5,n6,n7,n8,-1+n9,(1+n10),n11,n12,n13,n14)*rat(-n5,1)
      +`Z'(n1,n2,n3,(n4-1),1+n5,n6,n7,n8,n9,-1+(1+n10),n11,n12,n13,n14)*rat(-n5,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,1+n7,n8,-1+n9,(1+n10),n11,n12,n13,n14)*rat(-n7,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,1+n7,n8,n9,-1+(1+n10),n11,n12,n13,n14)*rat(-n7,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,1+n7,n8,n9,(1+n10),-1+n11,n12,n13,n14)*rat(-n7,1)
      +`Z'(n1,n2,n3,(n4-1),1+n5,n6,n7,n8,n9,(1+n10),n11,n12,n13,n14)*rat(2*n5,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,1+n7,n8,n9,(1+n10),n11,n12,n13,n14)*rat(2*n7,1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,1+(1+n10),n11,n12,n13,n14)*rat(-2*(1+n10),1)
      +`Z'(n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,(1+n10),n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(-1+n1,n2,n3,(n4-1),n5,n6,n7,n8,n9,(1+n10),n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(n1,n2,-1+n3,(n4-1),n5,n6,n7,n8,n9,(1+n10),n11,n12,n13,1+n14)*rat(n14,1)
      );
#endif
*--#] stat file : 
#endif
#if ( `il' >= 2 )
*--#[ stat2 file :
* G250
id,ifmatch->bebe1,`Z'(n1?pos_,n2?{>1},n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg_,n14?neg0_) =
  -rat(1,1-n2)*(`Z'(-2+n1,-1+n2,-1+n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n5*n14,-1+n3)+`Z'(-2+n1,-1+n2,-1+
   n3,n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n5*n14,-1+n3)+`Z'(-2+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+
   n11,n12,1+n13,1+n14)*rat(n9*n14,-1+n3)+`Z'(-2+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n12*n14,-1+
   n3)+`Z'(-2+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-2*ep*n14-n5*n14-2*n6*n14-n11*n14-n12*n14-n13*
   n14+3*n14,-1+n3)+`Z'(-2+n1,-1+n2,n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n5,1)+`Z'(-2+n1,-1+n2,n3,n4,1+n5,
   n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n5,1)+`Z'(-2+n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,1+n13,n14)*rat(-
   n9,1)+`Z'(-2+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)+`Z'(-2+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,
   n11,n12,1+n13,n14)*rat(2*ep+n5+2*n6+n11+n12+n13-3,1)+`Z'(-1+n1,-1+n2,-2+n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+
   n14)*rat(n5*n14,-1+n3)+`Z'(-1+n1,-1+n2,-2+n3,n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n5*n14,-1+n3)+`Z'(-1+
   n1,-1+n2,-2+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,1+n13,1+n14)*rat(-n9*n14,-1+n3)+`Z'(-1+n1,-1+n2,-2+n3,n4,n5,n6,n7,n8,
   n9,n10,n11,1+n12,n13,1+n14)*rat(-n12*n14,-1+n3)+`Z'(-1+n1,-1+n2,-2+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(2*
   ep*n14+n5*n14+2*n6*n14+n11*n14+n12*n14+n13*n14-3*n14,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,1+n4,1+n5,-1+n6,n7,n8,n9,-1+n10,n11,
   n12,1+n13,n14)*rat(-n4*n5,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,1+n4,1+n5,n6,-1+n7,n8,n9,-1+n10,n11,n12,1+n13,n14)*rat(n4*n5,-1+
   n3)+`Z'(-1+n1,-1+n2,-1+n3,1+n4,n5,n6,n7,n8,1+n9,-1+n10,-1+n11,n12,1+n13,n14)*rat(n4*n9,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,1+n4,
   n5,n6,n7,n8,n9,-1+n10,n11,1+n12,1+n13,n14)*rat(-n4*n12,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,1+n4,n5,n6,n7,n8,n9,-1+n10,n11,1+
   n12,n13,n14)*rat(n4*n12,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,1+n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,n14)*rat(-2*ep*n4-n4*n5-
   2*n4*n6-n4*n11-n4*n12-n4*n13+3*n4,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,1+n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-2*
   n5*n10,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-2*n5*n14,-1+n3)+`Z'(-1+n1,-1+n2,
   -1+n3,n4,1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(8*ep*n5+2*n1*n5+2*n2*n5+n3*n5+2*n4*n5+2*n5^2+2*n5*n6+2*n5*n7+
   2*n5*n8+n5*n9+n5*n10+n5*n11+2*n5*n12+2*n5*n13+2*n5*n14-19*n5,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,1+n5,n6,-1+n7,n8,n9,1+n10,
   n11,n12,1+n13,n14)*rat(2*n5*n10,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(2*n5*
   n14,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-8*ep*n5-2*n1*n5-2*n2*n5-n3*n5-2*n4*
   n5-2*n5^2-2*n5*n6-2*n5*n7-2*n5*n8-n5*n9-n5*n10-n5*n11-2*n5*n12-2*n5*n13-2*n5*n14+19*n5,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,
   n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,1+n13,n14)*rat(2*n9*n10,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,
   n12,1+n13,1+n14)*rat(2*n9*n14,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,1+n13,n14)*rat(-8*ep*n9-
   2*n1*n9-2*n2*n9-n3*n9-2*n4*n9-2*n5*n9-2*n6*n9-2*n7*n9-2*n8*n9-n9^2-n9*n10-n9*n11-2*n9*n12-2*n9*n13-2*n9*n14+19*n9,-1+
   n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n10*n12,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,
   n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-4*ep*n10-2*n5*n10-4*n6*n10-2*n10*n11-2*n10*n12-2*n10*n13+6*n10,-1+n3)+`Z'(-1+n1,
   -1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(8*ep*n12+2*n1*n12+2*n2*n12+n3*n12+2*n4*n12+2*n5*n12+2*n6*
   n12+2*n7*n12+2*n8*n12+n9*n12+n10*n12+n11*n12+2*n12^2+2*n12*n13+2*n12*n14-17*n12,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,
   n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(2*n12*n14,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*
   rat(-8*ep*n12-2*n1*n12-2*n2*n12-n3*n12-2*n4*n12-2*n5*n12-2*n6*n12-2*n7*n12-2*n8*n12-n9*n12-n10*n12-n11*n12-2*n12^2-2*
   n12*n13-2*n12*n14+19*n12,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-4*ep*n14-2*n5*
   n14-4*n6*n14-2*n11*n14-2*n12*n14-2*n13*n14+6*n14,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*
   rat(16*ep^2+4*ep*n1+4*ep*n2+2*ep*n3+4*ep*n4+12*ep*n5+20*ep*n6+4*ep*n7+4*ep*n8+2*ep*n9+2*ep*n10+10*ep*n11+12*ep*n12+
   12*ep*n13+4*ep*n14-62*ep+2*n1*n5+4*n1*n6+2*n1*n11+2*n1*n12+2*n1*n13-6*n1+2*n2*n5+4*n2*n6+2*n2*n11+2*n2*n12+2*n2*n13-
   6*n2+n3*n5+2*n3*n6+n3*n11+n3*n12+n3*n13-3*n3+2*n4*n5+4*n4*n6+2*n4*n11+2*n4*n12+2*n4*n13-6*n4+2*n5^2+6*n5*n6+2*n5*n7+
   2*n5*n8+n5*n9+n5*n10+3*n5*n11+4*n5*n12+4*n5*n13+2*n5*n14-25*n5+4*n6^2+4*n6*n7+4*n6*n8+2*n6*n9+2*n6*n10+4*n6*n11+6*n6*
   n12+6*n6*n13+4*n6*n14-44*n6+2*n7*n11+2*n7*n12+2*n7*n13-6*n7+2*n8*n11+2*n8*n12+2*n8*n13-6*n8+n9*n11+n9*n12+n9*n13-3*
   n9+n10*n11+n10*n12+n10*n13-3*n10+n11^2+3*n11*n12+3*n11*n13+2*n11*n14-22*n11+2*n12^2+4*n12*n13+2*n12*n14-25*n12+2*n13^
   2+2*n13*n14-25*n13-6*n14+57,-1+n3)+`Z'(-1+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(-1+n1,
   n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n2*n14-n14,-1+n3)+`Z'(n1,-2+n2,-1+n3,1+n4,n5,n6,n7,n8,n9,-1+n10,
   n11,1+n12,1+n13,n14)*rat(-n4*n12,-1+n3)+`Z'(n1,-2+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(8*ep*n12+2*
   n1*n12+2*n2*n12+n3*n12+2*n4*n12+2*n5*n12+2*n6*n12+2*n7*n12+2*n8*n12+n9*n12+n10*n12+n11*n12+2*n12^2+2*n12*n13+2*n12*
   n14-17*n12,-1+n3)+`Z'(n1,-2+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,-2+n3,n4,1+n5,-1+
   n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n5*n14,-1+n3)+`Z'(n1,-1+n2,-2+n3,n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+
   n14)*rat(-n5*n14,-1+n3)+`Z'(n1,-1+n2,-2+n3,n4,n5,n6,n7,n8,1+n9,n10,-1+n11,n12,1+n13,1+n14)*rat(-n9*n14,-1+n3)+`Z'(n1,-1+
   n2,-2+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n12*n14,-1+n3)+`Z'(n1,-1+n2,-2+n3,n4,n5,n6,n7,n8,n9,n10,n11,
   n12,1+n13,1+n14)*rat(2*ep*n14+n5*n14+2*n6*n14+n11*n14+n12*n14+n13*n14-3*n14,-1+n3)+`Z'(n1,-1+n2,-1+n3,1+n4,n5,-1+n6,n7,
   n8,n9,-1+n10,n11,1+n12,1+n13,n14)*rat(-n4*n12,-1+n3)+`Z'(n1,-1+n2,-1+n3,1+n4,n5,-1+n6,n7,n8,n9,-1+n10,n11,n12,2+n13,
   n14)*rat(-2*n4*n13-2*n4,-1+n3)+`Z'(n1,-1+n2,-1+n3,1+n4,n5,n6,n7,-1+n8,n9,-1+n10,n11,1+n12,1+n13,n14)*rat(n4*n12,-1+n3)+
   `Z'(n1,-1+n2,-1+n3,n4,1+n5,-1+n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-2*n5*n10,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,1+n5,-1+n6,
   n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-n5*n14,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,1+n5,n6,-1+n7,n8,n9,1+n10,n11,n12,1+n13,n14)
   *rat(2*n5*n10,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,1+n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(n5*n14,-1+n3)+`Z'(n1,-1+n2,-
   1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(8*ep*n12+2*n1*n12+2*n2*n12+n3*n12+2*n4*n12+2*n5*n12+2*n6*n12+
   2*n7*n12+2*n8*n12+n9*n12+n10*n12+n11*n12+2*n12^2+2*n12*n13+2*n12*n14-17*n12,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,n5,-1+n6,n7,
   n8,n9,n10,n11,n12,2+n13,n14)*rat(16*ep*n13+16*ep+4*n1*n13+4*n1+4*n2*n13+4*n2+2*n3*n13+2*n3+4*n4*n13+4*n4+4*n5*n13+4*
   n5+4*n6*n13+4*n6+4*n7*n13+4*n7+4*n8*n13+4*n8+2*n9*n13+2*n9+2*n10*n13+2*n10+2*n11*n13+2*n11+4*n12*n13+4*n12+4*n13^2+4*
   n13*n14-30*n13+4*n14-34,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-8*ep*n12-2*n1*n12-
   2*n2*n12-n3*n12-2*n4*n12-2*n5*n12-2*n6*n12-2*n7*n12-2*n8*n12-n9*n12-n10*n12-n11*n12-2*n12^2-2*n12*n13-2*n12*n14+17*
   n12,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,1+n9,1+n10,-1+n11,n12,1+n13,n14)*rat(2*n9*n10,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,
   n5,n6,n7,n8,1+n9,n10,-1+n11,n12,1+n13,1+n14)*rat(n9*n14,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,
   n13,n14)*rat(2*n10*n12,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,1+n13,n14)*rat(-4*ep*n10-2*n5*n10-4*
   n6*n10-2*n10*n11-2*n10*n12-2*n10*n13+6*n10,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(
   n12*n14,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,1+n14)*rat(-2*ep*n14-n5*n14-2*n6*n14-n11*n14-n12*
   n14-n13*n14+3*n14,-1+n3)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,1+n13,n14)*rat(-n12,1)+`Z'(n1,-1+n2,n3,n4,n5,
   -1+n6,n7,n8,n9,n10,n11,n12,2+n13,n14)*rat(-2*n13-2,1)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,1+n13,n14)*
   rat(n12,1)+`Z'(n1,n2,-2+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n2*n14+n14,-1+n3)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,
   n8,n9,1+n10,n11,n12,n13,n14)*rat(2*n2*n10-2*n10,-1+n3)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n2*
   n14-n14,-1+n3));
* G260
id,ifmatch->bebe1,`Z'(n1?pos_,n2?pos_,n3?{>1},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>1},n9?neg0_,n10?neg0_,n11?neg0_,n12?neg_,n13?neg0_,n14?neg0_) =
  -rat(1,1-n8)*(`Z'(-3+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n8*n14+n14,-1+n3)+`Z'(-3+n1,n2,n3,n4,n5,
   n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n8,1)+`Z'(-2+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-
   n8*n14+n14,-1+n3)+`Z'(-2+n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-1+n8,1)+`Z'(-2+n1,n2,-2+n3,n4,n5,n6,
   n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n8*n14-n14,-1+n3)+`Z'(-2+n1,n2,-1+n3,1+n4,n5,n6,n7,n8,n9,-1+n10,n11,1+n12,n13,
   n14)*rat(-n4*n8+n4,-1+n3)+`Z'(-2+n1,n2,-1+n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n7*n14,-1+n3)+`Z'(-
   2+n1,n2,-1+n3,n4,n5,-1+n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n7*n14,-1+n3)+`Z'(-2+n1,n2,-1+n3,n4,n5,-1+n6,n7,
   n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n8*n14-n14,-1+n3)+`Z'(-2+n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*
   rat(2*ep*n14+2*n6*n14+n7*n14+n8*n14+n11*n14+n13*n14-5*n14,-1+n3)+`Z'(-2+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,
   n13,n14)*rat(-2*n8*n10+2*n10,-1+n3)+`Z'(-2+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-2*n8*n14+2*n14,-
   1+n3)+`Z'(-2+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(8*ep*n8-8*ep+2*n1*n8-2*n1+2*n2*n8-2*n2+n3*n8-n3+
   2*n4*n8-2*n4+2*n5*n8-2*n5+2*n6*n8-2*n6+2*n7*n8-2*n7+2*n8^2+n8*n9+n8*n10+n8*n11+2*n8*n12+2*n8*n13+2*n8*n14-21*n8-n9-
   n10-n11-2*n12-2*n13-2*n14+19,-1+n3)+`Z'(-2+n1,n2,n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(n7,1)+`Z'(-2+n1,
   n2,n3,n4,n5,-1+n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-n7,1)+`Z'(-2+n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,
   n13,n14)*rat(-n8+1,1)+`Z'(-2+n1,n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-2*ep-2*n6-n7-n8-n11-n13+5,1)+`Z'(-
   1+n1,-1+n2,-2+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n8*n14-n14,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,1+n4,n5,n6,n7,
   n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(-n4*n8+n4,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*
   rat(-2*n8*n10+2*n10,-1+n3)+`Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-2*n8*n14+2*n14,-1+n3)+
   `Z'(-1+n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(8*ep*n8-8*ep+2*n1*n8-2*n1+2*n2*n8-2*n2+n3*n8-n3+2*
   n4*n8-2*n4+2*n5*n8-2*n5+2*n6*n8-2*n6+2*n7*n8-2*n7+2*n8^2+n8*n9+n8*n10+n8*n11+2*n8*n12+2*n8*n13+2*n8*n14-21*n8-n9-n10-
   n11-2*n12-2*n13-2*n14+19,-1+n3)+`Z'(-1+n1,n2,-2+n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n7*n14,-1+n3)
   +`Z'(-1+n1,n2,-2+n3,n4,n5,-1+n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n7*n14,-1+n3)+`Z'(-1+n1,n2,-2+n3,n4,n5,-1+n6,
   n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n8*n14+n14,-1+n3)+`Z'(-1+n1,n2,-2+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,1+
   n14)*rat(-2*ep*n14-2*n6*n14-n7*n14-n8*n14-n11*n14-n13*n14+5*n14,-1+n3)+`Z'(-1+n1,n2,-2+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+
   n12,n13,1+n14)*rat(n8*n14-n14,-1+n3)+`Z'(-1+n1,n2,-1+n3,1+n4,-1+n5,n6,1+n7,-1+n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(-n4*
   n7,-1+n3)+`Z'(-1+n1,n2,-1+n3,1+n4,n5,-1+n6,1+n7,-1+n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(n4*n7,-1+n3)+`Z'(-1+n1,n2,-1+n3,1+
   n4,n5,-1+n6,n7,n8,n9,-1+n10,n11,1+n12,n13,n14)*rat(n4*n8-n4,-1+n3)+`Z'(-1+n1,n2,-1+n3,1+n4,n5,n6,n7,-1+n8,n9,-1+n10,
   n11,1+n12,n13,n14)*rat(2*ep*n4+2*n4*n6+n4*n7+n4*n8+n4*n11+n4*n13-5*n4,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,-1+n5,n6,1+n7,-1+n8,
   n9,1+n10,n11,1+n12,n13,n14)*rat(-2*n7*n10,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*
   rat(-2*n7*n14,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(8*ep*n7+2*n1*n7+2*n2*n7+
   n3*n7+2*n4*n7+2*n5*n7+2*n6*n7+2*n7^2+2*n7*n8+n7*n9+n7*n10+n7*n11+2*n7*n12+2*n7*n13+2*n7*n14-19*n7,-1+n3)+`Z'(-1+n1,n2,-
   1+n3,n4,n5,-1+n6,1+n7,-1+n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n7*n10,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,n5,-1+n6,1+n7,-1+n8,
   n9,n10,n11,1+n12,n13,1+n14)*rat(2*n7*n14,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,n5,-1+n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*
   rat(-8*ep*n7-2*n1*n7-2*n2*n7-n3*n7-2*n4*n7-2*n5*n7-2*n6*n7-2*n7^2-2*n7*n8-n7*n9-n7*n10-n7*n11-2*n7*n12-2*n7*n13-2*n7*
   n14+19*n7,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n8*n10-2*n10,-1+n3)+`Z'(-1+n1,n2,
   -1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(2*n8*n14-2*n14,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,
   n10,n11,1+n12,n13,n14)*rat(-8*ep*n8+8*ep-2*n1*n8+2*n1-2*n2*n8+2*n2-n3*n8+n3-2*n4*n8+2*n4-2*n5*n8+2*n5-2*n6*n8+2*n6-2*
   n7*n8+2*n7-2*n8^2-n8*n9-n8*n10-n8*n11-2*n8*n12-2*n8*n13-2*n8*n14+21*n8+n9+n10+n11+2*n12+2*n13+2*n14-19,-1+n3)+`Z'(-1+
   n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,1+n10,n11,1+n12,n13,n14)*rat(4*ep*n10+4*n6*n10+2*n7*n10+2*n8*n10+2*n10*n11+2*n10*
   n13-10*n10,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(4*ep*n14+4*n6*n14+2*n7*n14+2*n8*
   n14+2*n11*n14+2*n13*n14-10*n14,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-16*ep^2-4*ep*
   n1-4*ep*n2-2*ep*n3-4*ep*n4-4*ep*n5-20*ep*n6-12*ep*n7-12*ep*n8-2*ep*n9-2*ep*n10-10*ep*n11-4*ep*n12-12*ep*n13-4*ep*n14+
   78*ep-4*n1*n6-2*n1*n7-2*n1*n8-2*n1*n11-2*n1*n13+10*n1-4*n2*n6-2*n2*n7-2*n2*n8-2*n2*n11-2*n2*n13+10*n2-2*n3*n6-n3*n7-
   n3*n8-n3*n11-n3*n13+5*n3-4*n4*n6-2*n4*n7-2*n4*n8-2*n4*n11-2*n4*n13+10*n4-4*n5*n6-2*n5*n7-2*n5*n8-2*n5*n11-2*n5*n13+
   10*n5-4*n6^2-6*n6*n7-6*n6*n8-2*n6*n9-2*n6*n10-4*n6*n11-4*n6*n12-6*n6*n13-4*n6*n14+48*n6-2*n7^2-4*n7*n8-n7*n9-n7*n10-
   3*n7*n11-2*n7*n12-4*n7*n13-2*n7*n14+29*n7-2*n8^2-n8*n9-n8*n10-3*n8*n11-2*n8*n12-4*n8*n13-2*n8*n14+29*n8-n9*n11-n9*
   n13+5*n9-n10*n11-n10*n13+5*n10-n11^2-2*n11*n12-3*n11*n13-2*n11*n14+24*n11-2*n12*n13+10*n12-2*n13^2-2*n13*n14+29*n13+
   10*n14-95,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-2*n8*n10+2*n10,-1+n3)+`Z'(-1+n1,n2,-
   1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n8*n14+n14,-1+n3)+`Z'(-1+n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,
   n12,n13,1+n14)*rat(n8*n14-n14,-1+n3)+`Z'(n1,-1+n2,-2+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n8*n14-n14,-1+
   n3)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-2*n8*n10+2*n10,-1+n3)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,
   n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n8*n14+n14,-1+n3)+`Z'(n1,n2,-2+n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,1+
   n14)*rat(n7*n14,-1+n3)+`Z'(n1,n2,-2+n3,n4,n5,-1+n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n7*n14,-1+n3)+`Z'(n1,n2,-
   2+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n8*n14+n14,-1+n3)+`Z'(n1,n2,-2+n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,
   1+n12,n13,1+n14)*rat(-2*ep*n14-2*n6*n14-n7*n14-n8*n14-n11*n14-n13*n14+5*n14,-1+n3)+`Z'(n1,n2,-2+n3,n4,n5,n6,n7,n8,n9,
   n10,n11,n12,n13,1+n14)*rat(-n8*n14+n14,-1+n3)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,1+n7,-1+n8,n9,1+n10,n11,1+n12,n13,n14)*rat(-
   2*n7*n10,-1+n3)+`Z'(n1,n2,-1+n3,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(-n7*n14,-1+n3)+`Z'(n1,n2,-1+n3,n4,
   n5,-1+n6,1+n7,-1+n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n7*n10,-1+n3)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,1+n7,-1+n8,n9,n10,n11,
   1+n12,n13,1+n14)*rat(n7*n14,-1+n3)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,1+n10,n11,1+n12,n13,n14)*rat(2*n8*n10-2*n10,-1+
   n3)+`Z'(n1,n2,-1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,1+n12,n13,1+n14)*rat(n8*n14-n14,-1+n3)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,-1+
   n8,n9,1+n10,n11,1+n12,n13,n14)*rat(4*ep*n10+4*n6*n10+2*n7*n10+2*n8*n10+2*n10*n11+2*n10*n13-10*n10,-1+n3)+`Z'(n1,n2,-1+
   n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,1+n12,n13,1+n14)*rat(2*ep*n14+2*n6*n14+n7*n14+n8*n14+n11*n14+n13*n14-5*n14,-1+n3)+`Z'(
   n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(2*n8*n10-2*n10,-1+n3)+`Z'(n1,n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,
   n11,n12,n13,1+n14)*rat(n8*n14-n14,-1+n3));
* G273
id,ifmatch->bebe1,`Z'(n1?pos_,n2?{>1},n3?pos_,n4?{>1},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_,n10?neg0_,n11?neg0_,n12?neg0_,n13?neg0_,n14?neg0_) =
  -rat(1,1-n4)*(`Z'(-1+n1,-1+n2,n3,-1+n4,n5,n6,n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-n8*n12,-1+n2)+`Z'(-1+n1,-1+n2,n3,n4,
   n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n4*n12+n12,-1+n2)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,1+n8,n9,n10,n11,n12,n13,
   n14)*rat(-n8,1)+`Z'(-1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,-2+n2,n3,-1+n4,n5,n6,n7,1+
   n8,n9,n10,n11,1+n12,n13,n14)*rat(-n8*n12,-1+n2)+`Z'(n1,-2+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n4*n13+
   n13,-1+n2)+`Z'(n1,-1+n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n4*n12-n12,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,1+n5,-
   1+n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(n5*n8,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,1+n5,n6,-1+n7,1+n8,n9,n10,n11,n12,n13,n14)*
   rat(-n5*n8,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,n5,-1+n6,n7,1+n8,n9,n10,n11,1+n12,n13,n14)*rat(-n8*n12,-1+n2)+`Z'(n1,-1+n2,n3,-1+
   n4,n5,-1+n6,n7,1+n8,n9,n10,n11,n12,1+n13,n14)*rat(-2*n8*n13,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,1+n8,1+n9,n10,-1+n11,
   n12,n13,n14)*rat(-n8*n9,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,1+n8,n9,n10,n11,1+n12,-1+n13,n14)*rat(-n8*n12,-1+n2)+`Z'(
   n1,-1+n2,n3,-1+n4,n5,n6,n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(2*ep*n8+n2*n8+n5*n8+2*n6*n8+n8*n11+n8*n12+n8*n13-5*n8,-1+
   n2)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n8*n12,-1+n2)+`Z'(n1,-1+n2,n3,-1+n4,n5,n6,n7,n8,n9,
   n10,n11,n12,1+n13,n14)*rat(n4*n13-n13,-1+n2)+`Z'(n1,-1+n2,n3,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(n4*n6-n6,-
   1+n2)+`Z'(n1,-1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-2*n4*n13+2*n13,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,-1+n6,
   n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(-n4*n13+n13,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,1+n6,-1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(
   -n4*n6+n6,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n4*n13-n13,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,
   n6,n7,n8,-1+n9,n10,1+n11,n12,n13,n14)*rat(-n4*n11+n11,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,-1+n10,n11,n12,n13,
   n14)*rat(n4*n9-n9,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-2*n4*n9+2*n9,-1+n2)+`Z'(n1,-1+n2,
   n3,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,n12,n13,n14)*rat(-n4*n11+n11,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,
   1+n13,n14)*rat(-n4*n13+n13,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,n14)*rat(2*n4*n11-2*n11,-1+n2)+`Z'(
   n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,-1+n14)*rat(n4*n12-n12,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,
   n11,1+n12,n13,n14)*rat(-n4*n12+n12,-1+n2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,n14)*rat(n4*n13-n13,-1+
   n2)+`Z'(n1,-1+n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(2*ep*n4-2*ep+n2*n4-n2+2*n4*n5+n4*n6+n4*n9+n4*n12-5*n4-2*
   n5-n6-n9-n12+5,-1+n2)+`Z'(n1,n2,-1+n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)+`Z'(n1,n2,n3,-1+n4,n5,-1+n6,
   n7,1+n8,n9,n10,n11,n12,n13,n14)*rat(n8,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-2*n10,1)+`Z'(n1,
   n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)+`Z'(n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*
   rat(6*ep+2*n2+n4+2*n5+2*n6+2*n7+n8+n9+n10+n11+n12+2*n13+n14-13,1)+`Z'(n1,n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,
   n14)*rat(-1+n4,1));
*--#] stat2 file : 
#endif
#if ( `il' >= 3 )
*--#[ n1 :
*`SKIP2'
id,ifmatch->bebe1,`Z'(n1?{>1},n2?`POS2',n3?`POS3',n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(-1,(n1-1))*(
      +`Z'((n1-1),1+n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,n14)*rat(-n2,1)
      +`Z'(1+(n1-1),n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat((n1-1),1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,n10,-1+n11,n12,1+n13,n14)*rat(-n13,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(8*ep+(n1-1)+2*n2+2*n3+2*n4+2*n5+2*n6+2*n7+2*n8+n9+n10+n11+2*
      n12+2*n13+2*n14-16,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,1+n9,n10,n11,n12,n13,n14)*rat(-2*n9,1)
      +`Z'((n1-1),n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)
      +`Z'(-1+(n1-1),n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(-n12,1)
      +`Z'((n1-1),n2,-1+n3,n4,n5,n6,n7,n8,n9,n10,n11,1+n12,n13,n14)*rat(n12,1)
      );
*--#] n1 : 
#endif
#if ( `il' >= 4 )
*--#[ n3 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',n3?{>1},n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?`NEG14')
			 = -rat(-1,(n3-1))*(
      +`Z'(1,n2,(n3-1),1+n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,n14)*rat(-n4,1)
      +`Z'(1,n2,(n3-1),n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(8*ep+2*n2+(n3-1)+2*n4+2*n5+2*n6+2*n7+2*n8+n9+n10+n11+2*
      n12+2*n13+2*n14-14,1)
      +`Z'(1,n2,(n3-1),n4,n5,n6,n7,n8,n9,1+n10,n11,n12,n13,n14)*rat(-2*n10,1)
      +`Z'(1,n2,(n3-1),n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)
      +`Z'(1,n2,-1+(n3-1),n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(n14,1)
      +`Y'(0,n2,1+(n3-1),n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat((n3-1),1)
      +`Y'(0,n2,(n3-1),n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,1+n14)*rat(-n14,1)
      );
*--#] n3 : 
#endif
#if ( `il' >= 5 )
*--#[ n14 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`NEG9',n10?`NEG10',n11?`NEG11',n12?`NEG12',n13?`NEG13',n14?neg_)
			 = -rat(1,-2*ep-2*n6-n7-n8-n11-n13+4)*(
      +`Z'(1,n2,1,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,-1+n12,n13,(1+n14))*rat(2*n6,1)
      +`Z'(1,n2,1,n4,n5,1+n6,-1+n7,n8,n9,n10,n11,-1+n12,n13,(1+n14))*rat(-n6,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,-1+n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,-1+n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,-1+n10,n11,-1+n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,(1+n14))*rat(2*ep+2*n6+n7+n8+n11+n13-4,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,n13,(1+n14))*rat(-2*ep-2*n6-2*n7-n8-n11-n13+5,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,n10,1+n11,-1+n12,n13,(1+n14))*rat(2*n11,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,-1+n2,1,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-2*n6,1)
      +`Z'(1,-1+n2,1,n4,n5,1+n6,-1+n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n6,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(2*ep+2*n6+2*n7+n8+n11+n13-5,1)
      +`Z'(1,n2,1,-1+n4,n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,n2,1,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-2*ep-2*n6-n7-n8-n11-n13+4,1)
      +`Z'(1,n2,1,n4,-1+n5,1+n6,n7,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(n6,1)
      +`Z'(1,n2,1,n4,-1+n5,n6,1+n7,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(n7,1)
      +`Z'(1,n2,1,n4,-1+n5,n6,n7,n8,-1+n9,n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Z'(1,n2,1,n4,-1+n5,n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(-2*n13,1)
      +`Z'(1,n2,1,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(2*ep+2*n6+2*n8+n11+n13-4,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,1+n7,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n7,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,-1+n9,n10,1+n11,n12,n13,(1+n14))*rat(-2*n11,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,n9,-1+n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,n9,-1+n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,n9,n10,n11,-1+n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-4*ep-2*n6-2*n7-3*n8-2*n11-2*n13+9,1)
      +`Z'(1,n2,1,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(2*ep+n6+n7+n8+n11+n13-4,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,(1+n14))*rat(-n6+1,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-2*n11,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,n2,1,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(3*n11,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,-2+n2,1,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,-1+n2,1,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Z'(1,-1+n2,1,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(3*n13,1)
      +`Z'(1,-1+n2,1,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-3*n13,1)
      +`Z'(1,n2,1,-1+n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Z'(1,n2,1,n4,-1+n5,-1+n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(0,n2,1,n4,-1+n5,1+n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-2*n6,1)
      +`Y'(0,n2,1,n4,n5,1+n6,-1+n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(n6,1)
      +`Y'(0,n2,1,n4,n5,n6,n7,n8,-1+n9,n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Y'(0,n2,1,n4,n5,n6,n7,n8,n9,-1+n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Y'(0,n2,1,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(0,n2,1,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(2*ep+2*n6+2*n7+n8+n11+n13-5,1)
      +`Y'(1,n2,0,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,(1+n14))*rat(-2*ep-2*n6-n7-n8-n11-n13+4,1)
      +`Y'(0,n2,1,n4,n5,n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-2*n11,1)
      +`Y'(0,n2,1,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(0,-1+n2,1,n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(0,n2,1,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(-n13,1)
      +`Y'(0,n2,1,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      +`Y'(0,n2,1,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,1+n13,(1+n14))*rat(n13,1)
      +`Y'(0,n2,1,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Y'(1,n2,0,n4,-1+n5,n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(n11,1)
      +`Y'(1,n2,0,n4,n5,-1+n6,n7,n8,n9,n10,1+n11,n12,n13,(1+n14))*rat(-n11,1)
      );
*--#] n14 : 
#endif
#if ( `il' >= 6 )
*--#[ n11 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,n4?`POS4',n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`NEG9',n10?`NEG10',n11?neg_,n12?`NEG12',n13?`NEG13',0)
			 = -rat(1,-4*ep-2*n2-2*n4-2*n5-n7-2*n9-2*n10-n12+9)*(
      +`Z'(1,n2,1,1+n4,n5,n6,-1+n7,n8,n9,-1+n10,(1+n11),n12,n13,0)*rat(-n4,1)
      +`Z'(1,n2,1,n4,1+n5,n6,-1+n7,n8,-1+n9,n10,(1+n11),n12,n13,0)*rat(-n5,1)
      +`Z'(1,n2,1,n4,1+n5,n6,-1+n7,n8,n9,-1+n10,(1+n11),n12,n13,0)*rat(-n5,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,-1+n9,n10,(1+n11),1+n12,-1+n13,0)*rat(-n12,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,-1+n9,n10,(1+n11),n12,n13,0)*rat(-2*n4+2*n6-n7-2*n10+2*(1+n11)+n12+1,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,-1+n10,(1+n11),n12,n13,0)*rat(-4*ep-2*n4-2*n5-3*n7-2*n10+9,1)
      +`Z'(1,n2,1,n4,1+n5,-1+n6,n7,n8,n9,n10,(1+n11),n12,n13,0)*rat(2*n5,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,-1+n9,1+n10,(1+n11),n12,n13,0)*rat(2*n10,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,1+n9,n10,(1+n11),n12,-1+n13,0)*rat(2*n9,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,1+n10,-1+(1+n11),n12,n13,0)*rat(2*n10,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,1+n10,(1+n11),n12,-1+n13,0)*rat(-2*n10,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,n10,1+(1+n11),n12,-1+n13,0)*rat(-2*(1+n11),1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,n10,(1+n11),1+n12,-1+n13,0)*rat(n12,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,n10,(1+n11),n12,n13,0)*rat(4*ep+4*n4+2*n5+2*n7+4*n10-10,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,-1+n9,n10,(1+n11),1+n12,n13,0)*rat(-n12,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,n9,n10,-1+(1+n11),n12,1+n13,0)*rat(2*n13,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,-1+n9,n10,(1+n11),1+n12,n13,0)*rat(-n12,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,-1+n9,n10,(1+n11),n12,1+n13,0)*rat(-2*n13,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,-1+n8,-1+n9,n10,(1+n11),1+n12,n13,0)*rat(n12,1)
      +`Z'(1,n2,1,n4,n5,n6,n7,n8,n9,1+n10,(1+n11),n12,n13,0)*rat(-2*n10,1)
      +`Z'(1,-1+n2,1,n4,n5,n6,n7,n8,n9,1+n10,(1+n11),n12,n13,0)*rat(-2*n10,1)
      +`Z'(1,n2,1,-1+n4,n5,n6,n7,n8,n9,1+n10,(1+n11),n12,n13,0)*rat(2*n10,1)
      +`Z'(1,n2,1,n4,n5,-1+n6,n7,n8,n9,1+n10,(1+n11),n12,n13,0)*rat(-2*n10,1)
      +`Y'(0,n2,1,n4,n5,n6,n7,n8,-1+n9,n10,(1+n11),1+n12,n13,0)*rat(-n12,1)
      +`Y'(0,n2,1,n4,n5,n6,n7,n8,n9,n10,(1+n11),1+n12,-1+n13,0)*rat(n12,1)
      +`Y'(1,n2,0,n4,n5,n6,n7,n8,n9,n10,(1+n11),1+n12,-1+n13,0)*rat(-n12,1)
      );
*--#] n11 : 
#endif
#if ( `il' >= 7 )
*--#[ n4 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,n4?{>1},n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`NEG9',n10?`NEG10',0,n12?`NEG12',n13?`NEG13',0)
			 = -rat(1,-2*ep-2*(n4-1)-2*n10+2)*(
      +`Z'(1,n2,1,(n4-1),1+n5,-1+n6,1+n7,n8,n9,n10,0,n12,n13,0)*rat(-2*n5*n7,(n4-1))
      +`Z'(1,n2,1,(n4-1),1+n5,n6,n7,n8,n9,n10,0,n12,n13,0)*rat(2*ep*n5+2*n5^2+2*n5*n7-2*n5,(n4-1))
      +`Z'(1,n2,1,(n4-1),n5,n6,1+n7,n8,n9,n10,0,n12,n13,0)*rat(2*ep*n7+2*n5*n7+2*n7^2-2*n7,(n4-1))
      +`Z'(1,n2,1,(n4-1),n5,n6,n7,n8,n9,2+n10,0,n12,n13,0)*rat(-2*n10^2-2*n10,(n4-1))
      );
*--#] n4 : 
#endif
#if ( `il' >= 8 )
*--#[ n13 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,1,n5?`POS5',n6?`POS6',
	n7?`POS7',n8?`POS8',n9?`NEG9',n10?`NEG10',0,n12?`NEG12',n13?neg_,0)
			 = -rat(1,4*ep+2*n6+3*n7+3*n8+2*(1+n13)-10)*(
      +`Z'(1,1+n2,1,1,-1+n5,n6,n7,n8,n9,n10,0,-1+n12,(1+n13),0)*rat(n2,1)
      +`Z'(1,n2,1,1,-1+n5,1+n6,n7,n8,n9,n10,0,-1+n12,(1+n13),0)*rat(10*ep*n6+4*n2*n6+2*n5*n6+6*n6^2+3*n6*n7+3*n6*n8+2*n6*
      n9+4*n6*n12+5*n6*(1+n13)-22*n6,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,1+n5,-1+n6,1+n7,-1+n8,n9,-1+n10,0,n12,(1+n13),0)*rat(n5*n7,ep+n10-1)
      +`Z'(1,n2,1,1,1+n5,-1+n6,n7,n8,n9,-1+n10,0,n12,(1+n13),0)*rat(-n5*n7+n5,ep+n10-1)
      +`Z'(1,n2,1,1,1+n5,n6,-1+n7,n8,n9,-1+n10,0,n12,(1+n13),0)*rat(2*ep*n5+n5^2+n5*n7+n5*n10-3*n5,ep+n10-1)
      +`Z'(1,n2,1,1,1+n5,n6,n7,-1+n8,n9,-1+n10,0,n12,(1+n13),0)*rat(-ep*n5-n5^2-n5*n7+n5,ep+n10-1)
      +`Z'(1,n2,1,1,n5,1+n6,-1+n7,n8,n9,n10,0,-1+n12,(1+n13),0)*rat(-6*ep*n6-2*n2*n6-n5*n6-4*n6^2-2*n6*n7-2*n6*n8-n6*n9-2*
      n6*n12-3*n6*(1+n13)+13*n6,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,n6,1+n7,-1+n8,n9,-1+n10,0,n12,(1+n13),0)*rat(-ep*n7-n5*n7-n7^2+n7,ep+n10-1)
      +`Z'(1,n2,1,1,n5,n6,n7,n8,1+n9,-1+n10,0,-1+n12,(1+n13),0)*rat(n9,1)
      +`Z'(1,n2,1,1,n5,n6,n7,n8,n9,-1+n10,0,-1+n12,1+(1+n13),0)*rat(-4*ep*(1+n13)-2*n2*(1+n13)-n5*(1+n13)-2*n6*(1+n13)-n7*(1+n13)-n8*(1+n13)-n9*
      (1+n13)-2*n12*(1+n13)-2*(1+n13)^2+9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,n6,n7,n8,n9,-1+n10,0,n12,(1+n13),0)*rat(4*ep^2+2*ep*n2+2*ep*n6+2*ep*n7+2*ep*n8+ep*n9+4*ep*n10+2*ep*
      n12+2*ep*(1+n13)-14*ep+2*n2*n10-2*n2+n5*n7-n5+2*n6*n10-2*n6+n7^2+n7*n10-4*n7+2*n8*n10-2*n8+n9*n10-n9+2*n10*n12+2*n10*
      (1+n13)-9*n10-2*n12-2*(1+n13)+11,ep+n10-1)
      +`Z'(1,n2,1,1,n5,n6,n7,n8,n9,n10,0,-1+n12,(1+n13),0)*rat(4*ep*n5-6*ep*n7-2*ep*n8+2*ep*n9-4*ep*(1+n13)+4*ep-2*n2*n7+2*n2+
      4*n5*n6+n5*n7+2*n5*n8+2*n5*(1+n13)-7*n5-4*n6*n7-2*n6*n8+2*n6*n9-4*n6*(1+n13)+2*n6-2*n7^2-3*n7*n8-2*n7*n12-5*n7*(1+n13)+14*n7-
      n8^2+n8*n9-3*n8*(1+n13)+5*n8+n9*(1+n13)-3*n9+2*n12-2*(1+n13)^2+10*(1+n13)-9,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,1+n5,n6,-1+n7,n8,n9,n10,0,n12,(1+n13),0)*rat(-n5,1)
      +`Z'(1,n2,1,1,n5,n6,n7,n8,1+n9,-1+n10,0,n12,(1+n13),0)*rat(-n9,1)
      +`Z'(1,n2,1,1,n5,n6,n7,n8,1+n9,n10,0,-1+n12,(1+n13),0)*rat(-2*n9,1)
      +`Z'(1,n2,1,1,n5,n6,n7,n8,n9,n10,0,-1+n12,1+(1+n13),0)*rat(4*ep*(1+n13)+2*n2*(1+n13)+n5*(1+n13)+2*n6*(1+n13)+n7*(1+n13)+n8*(1+n13)+n9*(1+n13)+2*
      n12*(1+n13)+2*(1+n13)^2-9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,n6,n7,n8,n9,n10,0,n12,(1+n13),0)*rat(-4*ep-2*n2-2*n6-n7-2*n8-n9-2*n12-2*(1+n13)+9,1)
      +`Z'(1,-1+n2,1,1,-1+n5,1+n6,n7,n8,n9,n10,0,n12,(1+n13),0)*rat(-10*ep*n6-4*n2*n6-2*n5*n6-6*n6^2-3*n6*n7-3*n6*n8-2*n6*
      n9-4*n6*n12-5*n6*(1+n13)+22*n6,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,-1+n2,1,1,1+n5,n6,-1+n7,n8,n9,n10,0,n12,(1+n13),0)*rat(n5,1)
      +`Z'(1,-1+n2,1,1,n5,1+n6,-1+n7,n8,n9,n10,0,n12,(1+n13),0)*rat(6*ep*n6+2*n2*n6+n5*n6+4*n6^2+2*n6*n7+2*n6*n8+n6*n9+2*n6
      *n12+3*n6*(1+n13)-13*n6,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,-1+n2,1,1,n5,n6,n7,n8,n9,-1+n10,0,n12,1+(1+n13),0)*rat(4*ep*(1+n13)+2*n2*(1+n13)+n5*(1+n13)+2*n6*(1+n13)+n7*(1+n13)+n8*(1+n13)+n9*(1+n13)+
      2*n12*(1+n13)+2*(1+n13)^2-9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,-1+n2,1,1,n5,n6,n7,n8,n9,n10,0,-1+n12,1+(1+n13),0)*rat(-8*ep*(1+n13)-2*n2*(1+n13)-n5*(1+n13)-6*n6*(1+n13)-3*n7*(1+n13)-3*n8*(1+n13)-n9
      *(1+n13)-2*n12*(1+n13)-4*(1+n13)^2+17*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,-1+n2,1,1,n5,n6,n7,n8,n9,n10,0,n12,(1+n13),0)*rat(16*ep^2+4*ep*n2+24*ep*n6+20*ep*n7+16*ep*n8+2*ep*n9+4*ep*n12+
      20*ep*(1+n13)-72*ep+4*n2*n6+4*n2*n7+2*n2*n8+2*n2*(1+n13)-10*n2+n5*n7-n5+8*n6^2+14*n6*n7+12*n6*n8+2*n6*n9+4*n6*n12+16*n6*
      (1+n13)-54*n6+5*n7^2+9*n7*n8+2*n7*n9+4*n7*n12+12*n7*(1+n13)-44*n7+4*n8^2+n8*n9+2*n8*n12+10*n8*(1+n13)-35*n8+n9*(1+n13)-5*n9+2*n12
      *(1+n13)-10*n12+6*(1+n13)^2-44*(1+n13)+81,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,-1+n5,1+n6,n7,-1+n8,n9,n10,0,n12,(1+n13),0)*rat(4*ep*n6+2*n2*n6+n5*n6+2*n6^2+n6*n7+n6*n8+n6*n9+2*n6*n12
      +2*n6*(1+n13)-9*n6,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,-1+n5,n6,1+n7,-1+n8,n9,n10,0,n12,(1+n13),0)*rat(4*ep*n7+2*n2*n7+n5*n7+2*n6*n7+n7^2+n7*n8+n7*n9+2*n7*n12
      +2*n7*(1+n13)-9*n7,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,-1+n5,n6,n7,n8,n9,n10,0,-1+n12,1+(1+n13),0)*rat(-12*ep*(1+n13)-4*n2*(1+n13)-2*n5*(1+n13)-8*n6*(1+n13)-4*n7*(1+n13)-4*n8*(1+n13)
      -2*n9*(1+n13)-4*n12*(1+n13)-6*(1+n13)^2+26*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,-1+n5,n6,n7,n8,n9,n10,0,n12,(1+n13),0)*rat(4*ep^2+6*ep*n6+8*ep*n8+2*ep*n12+6*ep*(1+n13)-16*ep-2*n2*n7+2*n2*
      n8-n5*n7+n5*n8+2*n6^2+n6*n7+5*n6*n8+2*n6*n12+5*n6*(1+n13)-12*n6+2*n7*n8-n7*n9-n7*n12+n7*(1+n13)+n7+2*n8^2+n8*n9+3*n8*n12+
      5*n8*(1+n13)-17*n8+n12*(1+n13)-4*n12+2*(1+n13)^2-12*(1+n13)+16,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,-1+n6,1+n7,-1+n8,n9,n10,0,n12,(1+n13),0)*rat(-4*ep*n7-2*n2*n7-n5*n7-2*n6*n7-n7^2-n7*n8-n7*n9-2*n7*
      n12-2*n7*(1+n13)+9*n7,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,-1+n6,n7,n8,n9,-1+n10,0,n12,1+(1+n13),0)*rat(-4*ep*(1+n13)-2*n2*(1+n13)-n5*(1+n13)-2*n6*(1+n13)-n7*(1+n13)-n8*(1+n13)-n9*(1+n13)
      -2*n12*(1+n13)-2*(1+n13)^2+9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,-1+n6,n7,n8,n9,n10,0,-1+n12,1+(1+n13),0)*rat(2*n2*(1+n13)+n5*(1+n13)-2*n6*(1+n13)-n7*(1+n13)-n8*(1+n13)+n9*(1+n13)+2*n12*(1+n13)
      -(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,-1+n6,n7,n8,n9,n10,0,n12,(1+n13),0)*rat(-8*ep^2-8*ep*n2-4*ep*n5-4*ep*n6-2*ep*n7-6*ep*n8-4*ep*n9-8*ep
      *n12-8*ep*(1+n13)+36*ep-4*n2*n6-4*n2*n7-6*n2*n8-4*n2*(1+n13)+18*n2-2*n5*n6-2*n5*n7-3*n5*n8-2*n5*(1+n13)+9*n5+2*n6*n7-2*n6*n9-
      4*n6*n12-2*n6*(1+n13)+8*n6+n7^2+n7*n8-2*n7*n9-4*n7*n12-n7*(1+n13)+5*n7-3*n8*n9-6*n8*n12-3*n8*(1+n13)+14*n8-2*n9*(1+n13)+9*n9-4*
      n12*(1+n13)+18*n12-2*(1+n13)^2+18*(1+n13)-41,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,n6,-1+n7,n8,n9,n10,0,-1+n12,1+(1+n13),0)*rat((1+n13),1)
      +`Z'(1,n2,1,1,n5,n6,-1+n7,n8,n9,n10,0,n12,(1+n13),0)*rat(4*ep*n2-2*ep*n6-2*ep*n7+2*ep*n9-2*ep*n10+4*ep*n12+2*ep*(1+n13)-4
      *ep+2*n2*n6+2*n2*n7+2*n2*n8+2*n2*(1+n13)-8*n2-n5*n6-2*n6*n7+n6*n9-2*n6*n10+2*n6*n12+n6*(1+n13)+n6-n7^2-n7*n8+n7*n9-n7*n10
      +2*n7*n12+2*n7+n8*n9-n8*n10+2*n8*n12+n8*(1+n13)-2*n8+n9*(1+n13)-4*n9-n10*(1+n13)+4*n10+2*n12*(1+n13)-8*n12+(1+n13)^2-6*(1+n13)+8,2*ep+2*
      n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,n6,n7,-1+n8,n9,-1+n10,0,n12,1+(1+n13),0)*rat((1+n13),1)
      +`Z'(1,n2,1,1,n5,n6,n7,-1+n8,n9,n10,0,n12,(1+n13),0)*rat(8*ep^2+4*ep*n2+4*ep*n5+6*ep*n6+6*ep*n7+4*ep*n8+2*ep*n9+2*ep*
      n10+2*ep*n12+6*ep*(1+n13)-26*ep+2*n2*n6+2*n2*n7+2*n2*n8+2*n2*(1+n13)-6*n2+3*n5*n6+2*n5*n7+2*n5*n8+2*n5*(1+n13)-7*n5+2*n6*n7+
      n6*n9+2*n6*n10+n6*(1+n13)-7*n6+n7^2+n7*n8+n7*n9+n7*n10+n7*n12+2*n7*(1+n13)-10*n7+n8*n9+n8*n10+n8*n12+n8*(1+n13)-6*n8+n9*(1+n13)-3
      *n9+n10*(1+n13)-4*n10+n12*(1+n13)-2*n12+(1+n13)^2-9*(1+n13)+19,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,n6,n7,n8,1+n9,n10,0,n12,(1+n13),0)*rat(n9,1)
      +`Z'(1,-1+n2,1,1,n5,n6,n7,n8,1+n9,n10,0,n12,(1+n13),0)*rat(-n9,1)
      +`Z'(1,-1+n2,1,1,n5,n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(-4*ep*(1+n13)-2*n2*(1+n13)-n5*(1+n13)-2*n6*(1+n13)-n7*(1+n13)-n8*(1+n13)-n9*(1+n13)-2*
      n12*(1+n13)-2*(1+n13)^2+9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,-1+n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(4*ep*(1+n13)+2*n2*(1+n13)+n5*(1+n13)+2*n6*(1+n13)+n7*(1+n13)+n8*(1+n13)+n9*(1+n13)+2*
      n12*(1+n13)+2*(1+n13)^2-9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,n2,1,1,n5,n6,-1+n7,n8,n9,1+n10,0,n12,(1+n13),0)*rat(2*ep*n10+n10^2-n10,ep+n10-1)
      +`Z'(1,n2,1,1,n5,n6,n7,-1+n8,n9,1+n10,0,n12,(1+n13),0)*rat(-2*ep*n10-n10^2+n10,ep+n10-1)
      +`Z'(1,n2,1,1,n5,n6,n7,-1+n8,n9,n10,0,n12,1+(1+n13),0)*rat(-(1+n13),1)
      +`Z'(1,-2+n2,1,1,n5,n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(8*ep*(1+n13)+2*n2*(1+n13)+n5*(1+n13)+6*n6*(1+n13)+3*n7*(1+n13)+3*n8*(1+n13)+n9*(1+n13)
      +2*n12*(1+n13)+4*(1+n13)^2-17*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,-1+n2,1,1,-1+n5,n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(18*ep*(1+n13)+6*n2*(1+n13)+3*n5*(1+n13)+12*n6*(1+n13)+6*n7*(1+n13)+6*n8*(1+n13)
      +3*n9*(1+n13)+6*n12*(1+n13)+9*(1+n13)^2-39*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,-1+n2,1,1,n5,-1+n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(-8*ep*(1+n13)-6*n2*(1+n13)-3*n5*(1+n13)-2*n6*(1+n13)-n7*(1+n13)-n8*(1+n13)-3*n9
      *(1+n13)-6*n12*(1+n13)-4*(1+n13)^2+19*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Z'(1,-1+n2,1,1,n5,n6,-1+n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(-2*(1+n13),1)
      +`Z'(1,-1+n2,1,1,n5,n6,n7,-1+n8,n9,n10,0,n12,1+(1+n13),0)*rat(-(1+n13),1)
      +`Z'(1,n2,1,1,-1+n5,-1+n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(-2*ep*(1+n13)-2*n2*(1+n13)-n5*(1+n13)-n9*(1+n13)-2*n12*(1+n13)-(1+n13)^2+5*(1+n13),
      2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(0,n2,1,1,-1+n5,1+n6,n7,n8,n9,n10,0,n12,(1+n13),0)*rat(-6*ep*n6-4*n2*n6-2*n5*n6-2*n6^2-n6*n7-n6*n8-2*n6*n9-4*n6*
      n12-3*n6*(1+n13)+14*n6,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(0,n2,1,1,n5,1+n6,-1+n7,n8,n9,n10,0,n12,(1+n13),0)*rat(2*ep*n6+2*n2*n6+n5*n6+n6*n9+2*n6*n12+n6*(1+n13)-5*n6,2*ep+2*n6
      +n7+n8+(1+n13)-4)
      +`Y'(0,n2,1,1,n5,n6,n7,n8,1+n9,-1+n10,0,n12,(1+n13),0)*rat(-n9,1)
      +`Y'(0,n2,1,1,n5,n6,n7,n8,n9,-1+n10,0,1+n12,(1+n13),0)*rat(-2*n12,1)
      +`Y'(0,n2,1,1,n5,n6,n7,n8,n9,-1+n10,0,n12,1+(1+n13),0)*rat(2*n2*(1+n13)+n5*(1+n13)-2*n6*(1+n13)-n7*(1+n13)-n8*(1+n13)+n9*(1+n13)+2*n12*(1+n13)-
      (1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(0,n2,1,1,n5,n6,n7,n8,n9,n10,0,n12,(1+n13),0)*rat(2*ep*n7-2*ep+2*n2*n7-2*n2+n5*n7-n5+n7*n9+2*n7*n12+n7*(1+n13)-5*n7-
      n9-2*n12-(1+n13)+5,2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(1,n2,0,1,n5,n6,n7,n8,1+n9,-1+n10,0,n12,(1+n13),0)*rat(n9,1)
      +`Y'(1,n2,1,0,1+n5,n6,-1+n7,n8,n9,n10,0,n12,(1+n13),0)*rat(-n5,1)
      +`Y'(1,n2,1,0,n5,n6,n7,n8,n9,n10,0,-1+n12,1+(1+n13),0)*rat(4*ep*(1+n13)+2*n2*(1+n13)+n5*(1+n13)+2*n6*(1+n13)+n7*(1+n13)+n8*(1+n13)+n9*(1+n13)+2*
      n12*(1+n13)+2*(1+n13)^2-9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(1,n2,1,0,n5,n6,n7,n8,n9,n10,0,n12,(1+n13),0)*rat(-4*ep-2*n2-2*n6-n7-2*n8-n9-2*n12-2*(1+n13)+9,1)
      +`Y'(0,n2,1,1,n5,n6,n7,n8,1+n9,n10,0,n12,(1+n13),0)*rat(n9,1)
      +`Y'(0,n2,1,1,n5,n6,n7,n8,n9,n10,0,1+n12,(1+n13),0)*rat(2*n12,1)
      +`Y'(0,n2,1,1,n5,n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(-2*n2*(1+n13)-n5*(1+n13)+2*n6*(1+n13)+n7*(1+n13)+n8*(1+n13)-n9*(1+n13)-2*n12*(1+n13)+(1+n13),
      2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(1,n2,0,1,n5,n6,n7,n8,1+n9,n10,0,n12,(1+n13),0)*rat(-n9,1)
      +`Y'(1,n2,1,0,n5,n6,n7,n8,1+n9,n10,0,n12,(1+n13),0)*rat(n9,1)
      +`Y'(0,-1+n2,1,1,n5,n6,n7,n8,1+n9,n10,0,n12,(1+n13),0)*rat(-n9,1)
      +`Y'(0,-1+n2,1,1,n5,n6,n7,n8,n9,n10,0,1+n12,(1+n13),0)*rat(-2*n12,1)
      +`Y'(0,-1+n2,1,1,n5,n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(4*ep*(1+n13)+2*n2*(1+n13)+n5*(1+n13)+2*n6*(1+n13)+n7*(1+n13)+n8*(1+n13)+n9*(1+n13)+2*
      n12*(1+n13)+2*(1+n13)^2-9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(0,n2,1,0,n5,n6,n7,n8,1+n9,n10,0,n12,(1+n13),0)*rat(n9,1)
      +`Y'(0,n2,1,0,n5,n6,n7,n8,n9,n10,0,1+n12,(1+n13),0)*rat(2*n12,1)
      +`Y'(0,n2,1,0,n5,n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(-2*n2*(1+n13)-n5*(1+n13)+2*n6*(1+n13)+n7*(1+n13)+n8*(1+n13)-n9*(1+n13)-2*n12*(1+n13)+(1+n13),
      2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(0,n2,1,1,-1+n5,n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(2*ep*(1+n13)+2*n2*(1+n13)+n5*(1+n13)+n9*(1+n13)+2*n12*(1+n13)+(1+n13)^2-5*(1+n13),2*ep
      +2*n6+n7+n8+(1+n13)-4)
      +`Y'(1,-1+n2,0,1,n5,n6,n7,n8,1+n9,n10,0,n12,(1+n13),0)*rat(n9,1)
      +`Y'(1,-1+n2,1,0,n5,n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(-4*ep*(1+n13)-2*n2*(1+n13)-n5*(1+n13)-2*n6*(1+n13)-n7*(1+n13)-n8*(1+n13)-n9*(1+n13)-2*
      n12*(1+n13)-2*(1+n13)^2+9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(1,n2,0,0,n5,n6,n7,n8,1+n9,n10,0,n12,(1+n13),0)*rat(-n9,1)
      +`Y'(1,n2,1,0,n5,-1+n6,n7,n8,n9,n10,0,n12,1+(1+n13),0)*rat(4*ep*(1+n13)+2*n2*(1+n13)+n5*(1+n13)+2*n6*(1+n13)+n7*(1+n13)+n8*(1+n13)+n9*(1+n13)+2*
      n12*(1+n13)+2*(1+n13)^2-9*(1+n13),2*ep+2*n6+n7+n8+(1+n13)-4)
      +`Y'(1,n2,1,0,n5,n6,n7,-1+n8,n9,n10,0,n12,1+(1+n13),0)*rat(-(1+n13),1)
      );
*--#] n13 : 
#endif
#if ( `il' >= 9 )
*--#[ n6 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,1,n5?`POS5',n6?{>1},
	n7?`POS7',n8?`POS8',n9?`NEG9',n10?`NEG10',0,n12?`NEG12',0,0)
			 = -rat(n7,-2*ep*(n6-1)-2*(n6-1)^2-2*(n6-1)*n7+2*(n6-1))*(
      +`Z'(1,n2,1,1,-1+n5,1+(n6-1),1+n7,n8,n9,n10,0,n12,0,0)*rat(2*(n6-1),1)
      +`Z'(1,n2,1,1,n5,(n6-1),1+n7,n8,n9,n10,0,n12,0,0)*rat(-2*ep-2*(n6-1)-2*n7+2,1)
      +`Z'(1,n2,1,1,n5,(n6-1),n7,1+n8,n9,n10,0,n12,0,0)*rat(2*ep*n8+2*n8^2-2*n8,n7)
      );
*--#] n6 : 
#endif
#if ( `il' >= 10 )
*--#[ n5 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,1,n5?{>1},1,
	n7?`POS7',n8?`POS8',n9?`NEG9',n10?`NEG10',0,n12?`NEG12',0,0)
			 = -rat(1,-2*ep-2*(n5-1)-2*n7+2)*(
      +`Z'(1,1+n2,1,1,(n5-1),1,n7,n8,n9,n10,0,n12,0,0)*rat(2*ep*n2+2*n2^2+2*n2*n9+2*n2*n12-2*n2,(n5-1))
      +`Z'(1,n2,1,1,-1+(n5-1),2,1+n7,n8,n9,n10,0,n12,0,0)*rat(-2*ep^2*n7-2*ep*(n5-1)*n7-2*ep*n7^2+2*ep*n7,ep^2*(n5-1)+2*ep*(n5-1)*n7-
      ep*(n5-1)+(n5-1)*n7^2-(n5-1)*n7)
      +`Z'(1,n2,1,1,1+(n5-1),1,-1+n7,1+n8,n9,n10,0,n12,0,0)*rat(2*ep*n8+2*n8^2-2*n8,ep+n7-1)
      +`Z'(1,n2,1,1,(n5-1),1,1+n7,n8,n9,n10,0,n12,0,0)*rat(2*ep^2*n7+2*ep*(n5-1)*n7+2*ep*n7^2-2*ep*n7,ep*(n5-1)+(n5-1)*n7-(n5-1))
      +`Z'(1,n2,1,1,(n5-1),1,n7,1+n8,n9,n10,0,n12,0,0)*rat(-2*ep^3*n8-2*ep^2*(n5-1)*n8-2*ep^2*n7*n8-2*ep^2*n8^2+4*ep^2*n8-2*ep
      *(n5-1)*n8^2+2*ep*(n5-1)*n8-2*ep*n7*n8^2+2*ep*n7*n8+2*ep*n8^2-2*ep*n8,ep^2*(n5-1)+2*ep*(n5-1)*n7-ep*(n5-1)+(n5-1)*n7^2-(n5-1)*n7)
      +`Z'(1,n2,1,1,(n5-1),1,n7,n8,1+n9,n10,0,1+n12,0,0)*rat(2*n9*n12,(n5-1))
      +`Z'(1,n2,1,1,(n5-1),1,n7,n8,2+n9,n10,0,n12,0,0)*rat(2*n9^2+2*n9,(n5-1))
      +`Y'(0,n2,1,1,(n5-1),1,n7,n8,1+n9,n10,0,1+n12,0,0)*rat(2*n9*n12,(n5-1))
      +`Y'(0,n2,1,1,(n5-1),1,n7,n8,n9,n10,0,2+n12,0,0)*rat(2*n12^2+2*n12,(n5-1))
      +`Y'(1,n2,0,1,(n5-1),1,n7,n8,1+n9,n10,0,1+n12,0,0)*rat(-2*n9*n12,(n5-1))
      );
*--#] n5 : 
#endif
#if ( `il' >= 11 )
*--#[ n7 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,1,1,1,
	n7?{>1},n8?`POS8',n9?`NEG9',n10?`NEG10',0,n12?`NEG12',0,0)
			 = -rat(ep*n8+(n7-1)*n8-n8,-2*ep^2*(n7-1)-3*ep*(n7-1)^2+ep*(n7-1)-(n7-1)^3+(n7-1)^2)*(
      +`Z'(1,1+n2,1,1,1,1,(n7-1),n8,n9,n10,0,n12,0,0)*rat(2*ep*n2+n2*(n7-1)+n2*n8-2*n2,n8)
      +`Z'(1,n2,1,1,1,1,(n7-1),1+n8,n9,n10,0,n12,0,0)*rat(4*ep^2+ep*n2+3*ep*(n7-1)+2*ep*n8+ep*n9+ep*n12-6*ep+n2*(n7-1)-n2+(n7-1)*n8+(n7-1)
      *n9+(n7-1)*n12-2*(n7-1)-n8-n9-n12+2,ep+(n7-1)-1)
      +`Z'(1,n2,1,1,2,1,-1+(n7-1),1+n8,n9,n10,0,n12,0,0)*rat(-2*ep-(n7-1)-n8+2,ep+(n7-1)-1)
      +`Z'(1,n2,1,1,1,1,(n7-1),1+n8,1+n9,n10,0,n12,0,0)*rat(-n9,1)
      +`Y'(0,1+n2,1,1,1,1,(n7-1),1+n8,n9,n10,0,n12,0,0)*rat(-n2,1)
      +`Y'(1,1+n2,1,1,0,1,1+(n7-1),n8,n9,n10,0,n12,0,0)*rat(-n2*(n7-1),n8)
      +`Y'(1,1+n2,1,1,1,0,1+(n7-1),n8,n9,n10,0,n12,0,0)*rat(n2*(n7-1),n8)
      +`Y'(1,1+n2,1,1,1,0,(n7-1),1+n8,n9,n10,0,n12,0,0)*rat(n2,1)
      +`Y'(1,n2,1,1,0,2,1+(n7-1),n8,n9,n10,0,n12,0,0)*rat(ep*(n7-1),ep*n8+(n7-1)*n8-n8)
      +`Y'(1,n2,1,1,2,0,1+(n7-1),n8,n9,n10,0,n12,0,0)*rat((n7-1),n8)
      +`Y'(1,n2,1,1,2,0,(n7-1),1+n8,n9,n10,0,n12,0,0)*rat(1,1)
      +`Y'(0,n2,1,1,1,1,(n7-1),1+n8,1+n9,n10,0,n12,0,0)*rat(-n9,1)
      +`Y'(0,n2,1,1,1,1,(n7-1),1+n8,n9,n10,0,1+n12,0,0)*rat(-2*n12,1)
      +`Y'(1,n2,0,1,1,1,(n7-1),1+n8,1+n9,n10,0,n12,0,0)*rat(n9,1)
      );
*--#] n7 : 
#endif
#if ( `il' >= 12 )
*--#[ n8 :
*`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,1,1,1,
	1,n8?{>1},n9?`NEG9',n10?`NEG10',0,n12?`NEG12',0,0)
		 = -rat(ep+n10,
      12*ep^2*(n8-1)+3*ep*n2*(n8-1)+7*ep*(n8-1)^2+3*ep*(n8-1)*n9+7*ep*(n8-1)*n10+3*ep*(n8-1)*
      n12-10*ep*(n8-1)+n2*(n8-1)^2+n2*(n8-1)*n10-n2*(n8-1)+(n8-1)^3+(n8-1)^2*n9+2*(n8-1)^2*n10+(n8-1)^2*n12-3*(n8-1)^2+(n8-1)*n9*n10-(n8-1)*n9+(n8-1)*n10^2+(n8-1)*n10*n12-
      3*(n8-1)*n10-(n8-1)*n12+2*(n8-1))*(
      +`Z'(1,1+n2,1,1,1,1,1,(n8-1),n9,n10,0,n12,0,0)*rat(8*ep^2*n2+2*ep*n2^2+6*ep*n2*(n8-1)+2*ep*n2*n9+2*ep*n2*n10+2*ep*n2*n12
      -8*ep*n2+n2^2*(n8-1)-n2^2+n2*(n8-1)^2+n2*(n8-1)*n9+n2*(n8-1)*n10+n2*(n8-1)*n12-3*n2*(n8-1)-n2*n9-n2*n10-n2*n12+2*n2,ep+n10)
      +`Z'(1,n2,1,1,1,1,1,1+(n8-1),1+n9,n10,0,n12,0,0)*rat(-3*ep*(n8-1)*n9-(n8-1)^2*n9-(n8-1)*n9*n10+(n8-1)*n9,ep+n10)
      +`Z'(1,n2,1,1,1,1,1,1+(n8-1),n9,1+n10,0,n12,0,0)*rat(-(n8-1)*n10,1)
      +`Z'(1,n2,1,1,1,1,1,(n8-1),1+n9,n10,0,1+n12,0,0)*rat(2*ep*n9*n12+(n8-1)*n9*n12-n9*n12,ep+n10)
      +`Z'(1,n2,1,1,1,1,1,(n8-1),2+n9,n10,0,n12,0,0)*rat(2*ep*n9^2+2*ep*n9+(n8-1)*n9^2+(n8-1)*n9-n9^2-n9,ep+n10)
      +`Z'(1,n2,1,1,1,1,1,(n8-1),n9,2+n10,0,n12,0,0)*rat(-2*ep*n10^2-2*ep*n10-(n8-1)*n10^2-(n8-1)*n10+n10^2+n10,ep+n10)
      +`Y'(0,1+n2,1,1,1,1,1,1+(n8-1),n9,n10,0,n12,0,0)*rat(-3*ep*n2*(n8-1)-n2*(n8-1)^2-n2*(n8-1)*n10+n2*(n8-1),ep+n10)
      +`Y'(1,1+n2,1,1,0,1,2,(n8-1),n9,n10,0,n12,0,0)*rat(-3*ep*n2-n2*(n8-1)-n2*n10+n2,ep+n10)
      +`Y'(1,1+n2,1,1,1,0,1,1+(n8-1),n9,n10,0,n12,0,0)*rat(3*ep*n2*(n8-1)+n2*(n8-1)^2+n2*(n8-1)*n10-n2*(n8-1),ep+n10)
      +`Y'(1,1+n2,1,1,1,0,2,(n8-1),n9,n10,0,n12,0,0)*rat(3*ep*n2+n2*(n8-1)+n2*n10-n2,ep+n10)
      +`Y'(1,n2,0,2,1,1,1,1+(n8-1),n9,n10,0,n12,0,0)*rat(-(n8-1),1)
      +`Y'(1,n2,1,1,0,2,2,(n8-1),n9,n10,0,n12,0,0)*rat(-ep,ep+n10)
      +`Y'(1,n2,1,1,2,0,1,1+(n8-1),n9,n10,0,n12,0,0)*rat(2*ep*(n8-1)+(n8-1)^2-(n8-1),ep+n10)
      +`Y'(0,n2,1,1,1,1,1,1+(n8-1),1+n9,n10,0,n12,0,0)*rat(-3*ep*(n8-1)*n9-(n8-1)^2*n9-(n8-1)*n9*n10+(n8-1)*n9,ep+n10)
      +`Y'(0,n2,1,1,1,1,1,1+(n8-1),n9,1+n10,0,n12,0,0)*rat((n8-1)*n10,1)
      +`Y'(0,n2,1,1,1,1,1,1+(n8-1),n9,n10,0,1+n12,0,0)*rat(-6*ep*(n8-1)*n12-2*(n8-1)^2*n12-2*(n8-1)*n10*n12+2*(n8-1)*n12,ep+n10)
      +`Y'(1,n2,0,1,1,1,1,1+(n8-1),1+n9,n10,0,n12,0,0)*rat(3*ep*(n8-1)*n9+(n8-1)^2*n9+(n8-1)*n9*n10-(n8-1)*n9,ep+n10)
      +`Y'(1,n2,0,1,1,1,1,1+(n8-1),n9,1+n10,0,n12,0,0)*rat(-(n8-1)*n10,1)
      +`Y'(1,n2,1,1,0,2,1,(n8-1),n9,2+n10,0,n12,0,0)*rat(n10^2+n10,ep+n10)
      +`Y'(1,n2,1,1,1,1,0,1+(n8-1),n9,2+n10,0,n12,0,0)*rat(-(n8-1)*n10^2-(n8-1)*n10,ep+n10)
      +`Y'(1,n2,1,1,1,2,0,(n8-1),n9,2+n10,0,n12,0,0)*rat(-n10^2-n10,ep+n10)
      +`Y'(0,n2,1,1,1,1,1,(n8-1),1+n9,n10,0,1+n12,0,0)*rat(2*ep*n9*n12+(n8-1)*n9*n12-n9*n12,ep+n10)
      +`Y'(0,n2,1,1,1,1,1,(n8-1),n9,n10,0,2+n12,0,0)*rat(2*ep*n12^2+2*ep*n12+(n8-1)*n12^2+(n8-1)*n12-n12^2-n12,ep+n10)
      +`Y'(1,n2,0,1,1,1,1,(n8-1),1+n9,n10,0,1+n12,0,0)*rat(-2*ep*n9*n12-(n8-1)*n9*n12+n9*n12,ep+n10)
      );
*--#] n8 : 
#endif
#if ( `il' >= 13 )
*--#[ n10 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,1,1,1,
	1,1,n9?`NEG9',n10?neg_,0,n12?`NEG12',0,0)
			 = -rat(2,8*ep+2*n2+n9+2*(1+n10)+2*n12-6)*(
      +`Z'(1,n2,1,1,1,1,1,1,-1+n9,1+(1+n10),0,-1+n12,0,0)*rat(-n9*(1+n10)+(1+n10),4*ep+2*(1+n10))
      +`Z'(1,n2,1,1,1,1,1,1,-1+n9,(1+n10),0,n12,0,0)*rat(16*ep^2+4*ep*n2+4*ep*n9+4*ep*(1+n10)+4*ep*n12-12*ep+n9*(1+n10)-(1+n10),4*ep+2*(1+n10))
      +`Z'(1,n2,1,1,1,1,1,1,1+n9,-1+(1+n10),0,-1+n12,0,0)*rat(n9,2)
      +`Z'(1,n2,1,1,1,1,1,1,n9,(1+n10),0,-1+n12,0,0)*rat(-32*ep^2-8*ep*n2-6*ep*n9-24*ep*(1+n10)-8*ep*n12+24*ep-4*n2*(1+n10)-4*n9*
      (1+n10)-4*(1+n10)^2-4*(1+n10)*n12+13*(1+n10),4*ep+2*(1+n10))
      +`Z'(1,n2,1,1,1,1,1,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(4*ep*(1+n10)+2*n2*(1+n10)+n9*(1+n10)+2*(1+n10)*n12-5*(1+n10),4*ep+2*(1+n10))
      +`Z'(1,n2,1,1,1,1,1,1,1+n9,-1+(1+n10),0,n12,0,0)*rat(-n9,2)
      +`Z'(1,n2,1,1,1,1,1,1,1+n9,(1+n10),0,-1+n12,0,0)*rat(n9,1)
      +`Z'(1,n2,1,1,1,1,1,1,n9,1+(1+n10),-1,n12,0,0)*rat(2*ep*(1+n10)-(1+n10),2*ep+(1+n10))
      +`Z'(1,n2,1,1,1,1,1,1,n9,1+(1+n10),0,-1+n12,0,0)*rat(n9*(1+n10)-(1+n10),2*ep+(1+n10))
      +`Z'(1,n2,1,1,1,1,1,1,n9,1+(1+n10),0,n12,-1,0)*rat(-4*ep*(1+n10)+2*(1+n10),2*ep+(1+n10))
      +`Z'(1,n2,1,1,1,1,1,1,n9,1+(1+n10),0,n12,0,-1)*rat(-2*(1+n10),1)
      +`Z'(1,n2,1,1,1,1,1,1,n9,(1+n10),0,n12,0,0)*rat(-8*ep^2-4*ep*n2-6*ep*n9+12*ep-2*n2*(1+n10)-2*n9*(1+n10)+2*(1+n10)^2+5*(1+n10),4*ep+2*(1+n10))
      +`Z'(1,-1+n2,1,1,1,1,1,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(-8*ep*(1+n10)-2*n2*(1+n10)-n9*(1+n10)-2*(1+n10)*n12+7*(1+n10),4*ep+2*(1+n10))
      +`Z'(1,-1+n2,1,1,1,1,1,1,n9,1+(1+n10),-1,n12,0,0)*rat(-2*ep*(1+n10)+(1+n10),2*ep+(1+n10))
      +`Z'(1,-1+n2,1,1,1,1,1,1,n9,(1+n10),0,n12,0,0)*rat(8*ep+2*n2+n9+4*(1+n10)+2*n12-6,2)
      +`Z'(1,n2,1,1,1,1,1,1,1+n9,(1+n10),0,n12,0,0)*rat(n9,2)
      +`Z'(1,n2,1,1,1,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(-n9*(1+n10)+(1+n10),4*ep+2*(1+n10))
      +`Z'(1,-1+n2,1,1,1,1,1,1,1+n9,(1+n10),0,n12,0,0)*rat(-n9,2)
      +`Z'(1,-1+n2,1,1,1,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(n9*(1+n10)-(1+n10),4*ep+2*(1+n10))
      +`Z'(1,-1+n2,1,1,1,1,1,1,n9,(1+n10),0,1+n12,0,0)*rat(-n12,1)
      +`Y'(1,1+n2,1,0,1,1,1,1,n9,(1+n10),0,-1+n12,0,0)*rat(-2*n2,1)
      +`Y'(1,1+n2,1,1,0,1,1,1,-1+n9,1+(1+n10),0,-1+n12,0,0)*rat(-n2*(1+n10),4*ep+2*(1+n10))
      +`Y'(1,1+n2,1,1,0,1,1,1,n9,(1+n10),0,-1+n12,0,0)*rat(n2,2)
      +`Y'(1,1+n2,1,1,0,1,2,0,-1+n9,(1+n10),0,n12,0,0)*rat(-n2*(1+n10),4*ep^2+6*ep*(1+n10)+2*(1+n10)^2)
      +`Y'(1,1+n2,1,1,0,1,2,0,n9,-1+(1+n10),0,n12,0,0)*rat(-n2,2*ep+2*(1+n10)-2)
      +`Y'(1,1+n2,1,1,1,0,2,0,-1+n9,(1+n10),0,n12,0,0)*rat(n2*(1+n10),4*ep^2+6*ep*(1+n10)+2*(1+n10)^2)
      +`Y'(1,1+n2,1,1,1,0,2,0,n9,-1+(1+n10),0,n12,0,0)*rat(n2,2*ep+2*(1+n10)-2)
      +`Y'(1,1+n2,1,1,1,1,1,0,-1+n9,(1+n10),0,n12,0,0)*rat(4*ep^2*n2+9*ep*n2*(1+n10)+n2^2*(1+n10)+n2*n9*(1+n10)+2*n2*(1+n10)^2+n2*(1+n10)*n12-
      3*n2*(1+n10),4*ep^2+6*ep*(1+n10)+2*(1+n10)^2)
      +`Y'(1,1+n2,1,1,1,1,1,0,n9,-1+(1+n10),0,n12,0,0)*rat(3*ep*n2+n2^2+n2*n9+n2*n12-2*n2,2*ep+2*(1+n10)-2)
      +`Y'(1,n2,1,1,0,2,1,1,-1+n9,1+(1+n10),0,-1+n12,0,0)*rat(-5*ep*(1+n10)-2*n2*(1+n10)-n9*(1+n10)-2*(1+n10)*n12+5*(1+n10),4*ep^2+2*ep*(1+n10))
      +`Y'(1,n2,1,1,0,2,1,1,n9,(1+n10),0,-1+n12,0,0)*rat(13*ep+2*n2+n9+4*(1+n10)+2*n12-8,2*ep)
      +`Y'(1,n2,1,1,1,1,0,2,-1+n9,1+(1+n10),0,-1+n12,0,0)*rat(6*ep*(1+n10)+2*n2*(1+n10)+n9*(1+n10)+2*(1+n10)*n12-5*(1+n10),8*ep^2+4*ep*(1+n10))
      +`Y'(1,n2,1,1,1,1,0,2,n9,(1+n10),0,-1+n12,0,0)*rat(-14*ep-2*n2-n9-4*(1+n10)-2*n12+8,4*ep)
      +`Y'(1,n2,1,1,2,0,1,1,-1+n9,(1+n10),0,n12,0,0)*rat(-4*ep+2,2*ep+(1+n10))
      +`Y'(1,n2,1,1,2,1,0,1,-1+n9,(1+n10),0,n12,0,0)*rat(4*ep^2+2*ep*(1+n10)-2*ep-(1+n10)^2-2*(1+n10),4*ep^2+6*ep*(1+n10)+2*(1+n10)^2)
      +`Y'(1,n2,1,1,2,1,0,1,n9,-1+(1+n10),0,n12,0,0)*rat(2*ep+(1+n10)-1,2*ep+2*(1+n10)-2)
      +`Y'(1,n2,1,1,2,0,1,1,n9,(1+n10),0,n12,0,0)*rat(4*ep-2,2*ep+(1+n10))
      +`Y'(1,n2,1,1,2,1,0,1,-1+n9,1+(1+n10),0,n12,0,0)*rat((1+n10),4*ep+2*(1+n10))
      +`Y'(1,n2,1,1,2,1,0,1,n9,(1+n10),0,n12,0,0)*rat(-6*ep-(1+n10)+2,4*ep+2*(1+n10))
      +`Y'(0,n2,1,1,0,2,1,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(3*ep*(1+n10)+2*n2*(1+n10)+n9*(1+n10)+2*(1+n10)*n12-5*(1+n10),4*ep^2+2*ep*(1+n10))
      +`Y'(0,n2,1,1,0,2,1,1,n9,(1+n10),0,n12,0,0)*rat(-11*ep-2*n2-n9-4*(1+n10)-2*n12+8,2*ep)
      +`Y'(0,n2,1,1,1,1,1,1,-1+n9,(1+n10),0,1+n12,0,0)*rat(-2*ep*n12,2*ep+(1+n10))
      +`Y'(0,n2,1,1,1,1,1,1,1+n9,-1+(1+n10),0,n12,0,0)*rat(-n9,2)
      +`Y'(0,n2,1,1,1,1,1,1,n9,-1+(1+n10),0,1+n12,0,0)*rat(-n12,1)
      +`Y'(0,n2,1,1,1,1,1,1,n9,(1+n10),0,n12,0,0)*rat(8*ep*n2+8*ep*(1+n10)+4*ep*n12-12*ep+4*n2*(1+n10)+n9*(1+n10)+4*(1+n10)^2+2*(1+n10)*n12-7*
      (1+n10),4*ep+2*(1+n10))
      +`Y'(0,n2,1,1,1,2,0,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(-2*ep*(1+n10)-2*n2*(1+n10)-n9*(1+n10)-2*(1+n10)*n12+5*(1+n10),8*ep^2+4*ep*(1+n10))
      +`Y'(0,n2,1,1,1,2,0,1,n9,(1+n10),0,n12,0,0)*rat(10*ep+2*n2+n9+4*(1+n10)+2*n12-8,4*ep)
      +`Y'(1,-1+n2,1,1,0,2,1,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(5*ep*(1+n10)+2*n2*(1+n10)+n9*(1+n10)+2*(1+n10)*n12-5*(1+n10),4*ep^2+2*ep*(1+n10))
      +`Y'(1,-1+n2,1,1,0,2,1,1,n9,(1+n10),0,n12,0,0)*rat(-13*ep-2*n2-n9-4*(1+n10)-2*n12+8,2*ep)
      +`Y'(1,-1+n2,1,1,1,2,0,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(-6*ep*(1+n10)-2*n2*(1+n10)-n9*(1+n10)-2*(1+n10)*n12+5*(1+n10),8*ep^2+4*ep*(1+n10))
      +`Y'(1,-1+n2,1,1,1,2,0,1,n9,(1+n10),0,n12,0,0)*rat(14*ep+2*n2+n9+4*(1+n10)+2*n12-8,4*ep)
      +`Y'(1,-1+n2,1,1,2,0,1,1,n9,(1+n10),0,n12,0,0)*rat(4*ep-2,2*ep+(1+n10))
      +`Y'(1,-1+n2,1,1,2,1,0,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(-(1+n10),4*ep+2*(1+n10))
      +`Y'(1,-1+n2,1,1,2,1,0,1,n9,(1+n10),0,n12,0,0)*rat(-2*ep+(1+n10)+2,4*ep+2*(1+n10))
      +`Y'(1,n2,0,1,1,1,1,1,1+n9,-1+(1+n10),0,n12,0,0)*rat(n9,2)
      +`Y'(1,n2,0,1,1,1,1,1,n9,(1+n10),0,n12,0,0)*rat(-32*ep^2-8*ep*n2-24*ep*(1+n10)-12*ep*n12+20*ep-4*n2*(1+n10)-n9*(1+n10)-4*(1+n10)^2-6*
      (1+n10)*n12+11*(1+n10),4*ep+2*(1+n10))
      +`Y'(1,n2,1,0,0,1,2,1,n9,(1+n10),0,n12,0,0)*rat(2*ep-1,2*ep+(1+n10))
      +`Y'(1,n2,1,0,1,0,2,1,n9,(1+n10),0,n12,0,0)*rat(-2*ep+1,2*ep+(1+n10))
      +`Y'(1,n2,1,0,1,1,1,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(4*ep*(1+n10)+2*n2*(1+n10)+n9*(1+n10)+2*(1+n10)*n12-5*(1+n10),4*ep+2*(1+n10))
      +`Y'(1,n2,1,0,1,1,1,1,n9,(1+n10),0,n12,0,0)*rat(12*ep+6*n2+3*n9+2*n12-8,2)
      +`Y'(1,n2,1,0,2,0,1,1,n9,(1+n10),0,n12,0,0)*rat(-2*ep+1,2*ep+(1+n10))
      +`Y'(1,n2,1,0,2,1,0,1,-1+n9,1+(1+n10),0,n12,0,0)*rat((1+n10),4*ep+2*(1+n10))
      +`Y'(1,n2,1,0,2,1,0,1,n9,(1+n10),0,n12,0,0)*rat(-1,2)
      +`Y'(1,n2,1,1,0,1,1,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(2*ep*(1+n10)-(1+n10)*n12-(1+n10),4*ep+2*(1+n10))
      +`Y'(1,n2,1,1,0,1,1,1,n9,1+(1+n10),-1,n12,0,0)*rat(2*ep*(1+n10)-(1+n10),2*ep+(1+n10))
      +`Y'(1,n2,1,1,0,1,1,1,n9,(1+n10),0,n12,0,0)*rat(36*ep^2+26*ep*(1+n10)+2*ep*n12-26*ep+4*(1+n10)^2+(1+n10)*n12-9*(1+n10)+4,4*ep+2*(1+n10))
      +`Y'(1,n2,1,1,0,1,2,0,-1+n9,1+(1+n10),0,n12,0,0)*rat(-4*ep*(1+n10)-2*n2*(1+n10)-n9*(1+n10)-2*(1+n10)*n12+5*(1+n10),8*ep^2+4*ep*(1+n10))
      +`Y'(1,n2,1,1,0,1,2,0,n9,(1+n10),0,n12,0,0)*rat(12*ep+2*n2+n9+4*(1+n10)+2*n12-8,4*ep)
      +`Y'(1,n2,1,1,0,2,1,0,-1+n9,1+(1+n10),0,n12,0,0)*rat(-4*ep*(1+n10)-2*n2*(1+n10)-n9*(1+n10)-2*(1+n10)*n12+5*(1+n10),8*ep^2+4*ep*(1+n10))
      +`Y'(1,n2,1,1,0,2,1,0,n9,(1+n10),0,n12,0,0)*rat(12*ep+2*n2+n9+4*(1+n10)+2*n12-8,4*ep)
      +`Y'(1,n2,1,1,1,0,1,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(8*ep^2*(1+n10)+4*ep*n2*(1+n10)+2*ep*n9*(1+n10)+4*ep*(1+n10)*n12-14*ep*(1+n10)-2*n2*
      (1+n10)-n9*(1+n10)-2*(1+n10)*n12+5*(1+n10),4*ep^2+2*ep*(1+n10))
      +`Y'(1,n2,1,1,1,0,1,1,n9,(1+n10),0,n12,0,0)*rat(-28*ep^2-4*ep*n2-2*ep*n9-8*ep*(1+n10)-4*ep*n12+30*ep+2*n2+n9+4*(1+n10)+2*n12-
      8,2*ep)
      +`Y'(1,n2,1,1,1,0,2,0,-1+n9,1+(1+n10),0,n12,0,0)*rat(4*ep*(1+n10)+2*n2*(1+n10)+n9*(1+n10)+2*(1+n10)*n12-5*(1+n10),8*ep^2+4*ep*(1+n10))
      +`Y'(1,n2,1,1,1,0,2,0,n9,(1+n10),0,n12,0,0)*rat(-12*ep-2*n2-n9-4*(1+n10)-2*n12+8,4*ep)
      +`Y'(1,n2,1,1,1,1,0,1,-1+n9,1+(1+n10),0,n12,0,0)*rat(-8*ep^2*(1+n10)-4*ep*n2*(1+n10)-2*ep*n9*(1+n10)+2*ep*(1+n10)^2-4*ep*(1+n10)*n12+16*
      ep*(1+n10)+2*n2*(1+n10)+n9*(1+n10)+2*(1+n10)*n12-5*(1+n10),8*ep^2+4*ep*(1+n10))
      +`Y'(1,n2,1,1,1,1,0,1,n9,(1+n10),0,n12,0,0)*rat(16*ep^3+8*ep^2*n2+4*ep^2*n9+20*ep^2*(1+n10)+8*ep^2*n12-32*ep^2+4*ep*n2*
      (1+n10)-4*ep*n2+2*ep*n9*(1+n10)-2*ep*n9+6*ep*(1+n10)^2+4*ep*(1+n10)*n12-28*ep*(1+n10)-4*ep*n12+12*ep-2*n2*(1+n10)-n9*(1+n10)-4*(1+n10)^2-2*(1+n10)*
      n12+8*(1+n10),8*ep^2+4*ep*(1+n10))
      +`Y'(1,n2,1,1,1,1,1,0,-1+n9,1+(1+n10),0,n12,0,0)*rat(-4*ep*(1+n10)-2*n2*(1+n10)-n9*(1+n10)-(1+n10)^2-(1+n10)*n12+3*(1+n10),4*ep+2*(1+n10))
      +`Y'(1,n2,1,1,1,1,1,0,n9,(1+n10),0,n12,0,0)*rat(-12*ep-6*n2-3*n9-3*(1+n10)-3*n12+9,2)
      +`Y'(0,n2,1,1,1,1,1,1,-1+n9,1+(1+n10),0,1+n12,0,0)*rat(-(1+n10)*n12,2*ep+(1+n10))
      +`Y'(0,n2,1,1,1,1,1,1,1+n9,(1+n10),0,n12,0,0)*rat(n9,2)
      +`Y'(0,n2,1,1,1,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(-n9*(1+n10)+(1+n10),4*ep+2*(1+n10))
      +`Y'(0,n2,1,1,1,1,1,1,n9,(1+n10),0,1+n12,0,0)*rat(n12,1)
      +`Y'(1,n2,0,1,1,1,1,1,1+n9,(1+n10),0,n12,0,0)*rat(-n9,2)
      +`Y'(1,n2,0,1,1,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(n9*(1+n10)-(1+n10),4*ep+2*(1+n10))
      +`Y'(1,n2,1,0,1,1,1,1,1+n9,(1+n10),0,n12,0,0)*rat(-3*n9,2)
      +`Y'(1,n2,1,0,1,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(-8*ep*(1+n10)-n9*(1+n10)-4*(1+n10)^2+(1+n10),4*ep+2*(1+n10))
      +`Y'(1,n2,1,1,0,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(-4*ep*(1+n10)+2*(1+n10),2*ep+(1+n10))
      +`Y'(1,n2,1,1,1,0,1,1,n9,1+(1+n10),0,n12,0,0)*rat(-4*ep*(1+n10)+2*(1+n10),2*ep+(1+n10))
      +`Y'(1,n2,1,1,1,1,0,1,-1+n9,2+(1+n10),0,n12,0,0)*rat(-(1+n10)^2-(1+n10),2*ep+2*(1+n10))
      +`Y'(1,n2,1,1,1,1,0,1,n9,1+(1+n10),0,n12,0,0)*rat(12*ep^2*(1+n10)+12*ep*(1+n10)^2-14*ep*(1+n10)+(1+n10)^3-5*(1+n10)^2+4*(1+n10),4*ep^2+6*ep*
      (1+n10)-4*ep+2*(1+n10)^2-2*(1+n10))
      +`Y'(1,n2,1,1,1,1,1,0,-1+n9,2+(1+n10),0,n12,0,0)*rat((1+n10)^2+(1+n10),2*ep+2*(1+n10))
      +`Y'(1,n2,1,1,1,1,1,0,1+n9,-1+(1+n10),0,1+n12,0,0)*rat(n9*n12,2*ep+2*(1+n10)-2)
      +`Y'(1,n2,1,1,1,1,1,0,1+n9,(1+n10),0,n12,0,0)*rat(8*ep^2*n9+12*ep*n9*(1+n10)+n9^2*(1+n10)+4*n9*(1+n10)^2-n9*(1+n10),4*ep^2+6*ep*(1+n10)+2
      *(1+n10)^2)
      +`Y'(1,n2,1,1,1,1,1,0,2+n9,-1+(1+n10),0,n12,0,0)*rat(n9^2+n9,2*ep+2*(1+n10)-2)
      +`Y'(1,n2,1,1,1,1,1,0,n9,1+(1+n10),0,n12,0,0)*rat(2*ep*(1+n10)+(1+n10)^2-(1+n10),2*ep+2*(1+n10)-2)
      +`Y'(1,n2,1,1,1,1,1,0,n9,(1+n10),0,1+n12,0,0)*rat(4*ep^2*n12+6*ep*(1+n10)*n12+n9*(1+n10)*n12+2*(1+n10)^2*n12-(1+n10)*n12,4*ep^2+6*ep*
      (1+n10)+2*(1+n10)^2)
      +`Y'(0,-1+n2,1,1,1,1,1,1,-1+n9,1+(1+n10),0,1+n12,0,0)*rat((1+n10)*n12,2*ep+(1+n10))
      +`Y'(0,-1+n2,1,1,1,1,1,1,1+n9,(1+n10),0,n12,0,0)*rat(3*n9,2)
      +`Y'(0,-1+n2,1,1,1,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(n9*(1+n10)-(1+n10),4*ep+2*(1+n10))
      +`Y'(0,-1+n2,1,1,1,1,1,1,n9,(1+n10),0,1+n12,0,0)*rat(2*n12,1)
      +`Y'(0,n2,1,0,1,1,1,1,-1+n9,1+(1+n10),0,1+n12,0,0)*rat(-(1+n10)*n12,2*ep+(1+n10))
      +`Y'(0,n2,1,0,1,1,1,1,1+n9,(1+n10),0,n12,0,0)*rat(-3*n9,2)
      +`Y'(0,n2,1,0,1,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(-n9*(1+n10)+(1+n10),4*ep+2*(1+n10))
      +`Y'(0,n2,1,0,1,1,1,1,n9,(1+n10),0,1+n12,0,0)*rat(-3*n12,1)
      +`Y'(0,n2,1,1,1,1,1,0,-1+n9,(1+n10),0,2+n12,0,0)*rat((1+n10)*n12^2+(1+n10)*n12,4*ep^2+6*ep*(1+n10)+2*(1+n10)^2)
      +`Y'(0,n2,1,1,1,1,1,0,1+n9,-1+(1+n10),0,1+n12,0,0)*rat(n9*n12,2*ep+2*(1+n10)-2)
      +`Y'(0,n2,1,1,1,1,1,0,n9,-1+(1+n10),0,2+n12,0,0)*rat(n12^2+n12,2*ep+2*(1+n10)-2)
      +`Y'(0,n2,1,1,1,1,1,0,n9,(1+n10),0,1+n12,0,0)*rat(4*ep^2*n12+6*ep*(1+n10)*n12+n9*(1+n10)*n12+2*(1+n10)^2*n12-(1+n10)*n12,4*ep^2+6*ep*
      (1+n10)+2*(1+n10)^2)
      +`Y'(1,-1+n2,0,1,1,1,1,1,1+n9,(1+n10),0,n12,0,0)*rat(-3*n9,2)
      +`Y'(1,-1+n2,0,1,1,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(-n9*(1+n10)+(1+n10),4*ep+2*(1+n10))
      +`Y'(1,-1+n2,0,1,1,1,1,1,n9,(1+n10),0,1+n12,0,0)*rat(n12,1)
      +`Y'(1,n2,0,0,1,1,1,1,1+n9,(1+n10),0,n12,0,0)*rat(3*n9,2)
      +`Y'(1,n2,0,0,1,1,1,1,n9,1+(1+n10),0,n12,0,0)*rat(n9*(1+n10)-(1+n10),4*ep+2*(1+n10))
      +`Y'(1,n2,0,1,1,1,1,0,1+n9,-1+(1+n10),0,1+n12,0,0)*rat(-n9*n12,2*ep+2*(1+n10)-2)
      +`Y'(1,n2,0,1,1,1,1,0,n9,(1+n10),0,1+n12,0,0)*rat(-4*ep^2*n12-6*ep*(1+n10)*n12-n9*(1+n10)*n12-2*(1+n10)^2*n12+(1+n10)*n12,4*ep^2+6*ep
      *(1+n10)+2*(1+n10)^2)
      );
*--#] n10 : 
#endif
#if ( `il' >= 14 )
*--#[ n9 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?`POS2',1,1,1,1,1,1,n9?neg_,0,0,n12?`NEG12',0,0)
			 = -rat(1,-8*ep-2*n2-2*(1+n9)-2*n12+6)*(
      +`Z'(1,n2,1,1,1,1,1,1,(1+n9),0,0,-1+n12,0,0)*rat(8*ep+2*n2+2*(1+n9)+2*n12-6,1)
      +`Z'(1,n2,1,1,1,1,1,1,1+(1+n9),0,0,-1+n12,0,0)*rat(-2*(1+n9),1)
      +`Z'(1,n2,1,1,1,1,1,1,(1+n9),0,0,n12,0,0)*rat(2*ep+2*n2+2*(1+n9)-4,1)
      +`Z'(1,-1+n2,1,1,1,1,1,1,1+(1+n9),0,0,n12,0,0)*rat(2*(1+n9),1)
      +`Z'(1,-1+n2,1,1,1,1,1,1,(1+n9),0,0,1+n12,0,0)*rat(2*n12,1)
      +`Y'(1,1+n2,1,0,1,1,1,1,(1+n9),0,0,-1+n12,0,0)*rat(2*n2,1)
      +`Y'(1,1+n2,1,1,0,1,2,0,(1+n9),-1,0,n12,0,0)*rat(n2,ep-1)
      +`Y'(1,1+n2,1,1,1,0,1,1,-1+(1+n9),0,0,n12,0,0)*rat(n2,1)
      +`Y'(1,1+n2,1,1,1,0,2,0,(1+n9),-1,0,n12,0,0)*rat(-n2,ep-1)
      +`Y'(1,1+n2,1,1,1,1,1,0,-1+(1+n9),0,0,n12,0,0)*rat(-n2,1)
      +`Y'(1,1+n2,1,1,1,1,1,0,(1+n9),-1,0,n12,0,0)*rat(-3*ep*n2-n2^2-n2*(1+n9)-n2*n12+2*n2,ep-1)
      +`Y'(1,n2,1,1,0,2,1,1,(1+n9),0,0,-1+n12,0,0)*rat(-4*ep+2,ep)
      +`Y'(1,n2,1,1,1,1,0,2,(1+n9),0,0,-1+n12,0,0)*rat(2*ep-1,ep)
      +`Y'(1,n2,1,1,2,0,1,1,-1+(1+n9),0,0,n12,0,0)*rat(5*ep-2,ep)
      +`Y'(1,n2,1,1,2,0,1,1,(1+n9),-1,0,n12,0,0)*rat(1,1)
      +`Y'(1,n2,1,1,2,1,0,1,-1+(1+n9),0,0,n12,0,0)*rat(-2*ep+1,ep)
      +`Y'(1,n2,1,1,2,0,1,1,(1+n9),0,0,n12,0,0)*rat(-6*ep+2,ep)
      +`Y'(1,n2,1,1,2,1,0,1,(1+n9),0,0,n12,0,0)*rat(2*ep-1,ep)
      +`Y'(0,n2,1,1,0,2,1,1,(1+n9),0,0,n12,0,0)*rat(4*ep-2,ep)
      +`Y'(0,n2,1,1,1,1,1,1,-1+(1+n9),0,0,1+n12,0,0)*rat(2*n12,1)
      +`Y'(0,n2,1,1,1,1,1,1,(1+n9),0,0,n12,0,0)*rat(-2*ep-2*n2-2*n12+4,1)
      +`Y'(0,n2,1,1,1,2,0,1,(1+n9),0,0,n12,0,0)*rat(-2*ep+1,ep)
      +`Y'(1,-1+n2,1,1,0,2,1,1,(1+n9),0,0,n12,0,0)*rat(4*ep-2,ep)
      +`Y'(1,-1+n2,1,1,1,2,0,1,(1+n9),0,0,n12,0,0)*rat(-2*ep+1,ep)
      +`Y'(1,-1+n2,1,1,2,0,1,1,(1+n9),0,0,n12,0,0)*rat(-4*ep+2,ep)
      +`Y'(1,-1+n2,1,1,2,1,0,1,(1+n9),0,0,n12,0,0)*rat(2*ep-1,ep)
      +`Y'(1,n2,0,1,1,1,1,1,(1+n9),0,0,n12,0,0)*rat(10*ep+2*n2+4*n12-6,1)
      +`Y'(1,n2,1,0,0,1,2,1,(1+n9),0,0,n12,0,0)*rat(-2*ep+1,ep)
      +`Y'(1,n2,1,0,1,0,2,1,(1+n9),0,0,n12,0,0)*rat(2*ep-1,ep)
      +`Y'(1,n2,1,0,1,1,1,1,(1+n9),0,0,n12,0,0)*rat(-8*ep-4*n2-2*(1+n9)-2*n12+6,1)
      +`Y'(1,n2,1,0,2,0,1,1,(1+n9),0,0,n12,0,0)*rat(2*ep-1,ep)
      +`Y'(1,n2,1,1,0,1,1,1,(1+n9),0,0,n12,0,0)*rat(-12*ep^2+10*ep-2,ep)
      +`Y'(1,n2,1,1,0,1,2,0,(1+n9),0,0,n12,0,0)*rat(-2*ep+1,ep)
      +`Y'(1,n2,1,1,0,2,1,0,(1+n9),0,0,n12,0,0)*rat(-2*ep+1,ep)
      +`Y'(1,n2,1,1,1,0,1,1,(1+n9),0,0,n12,0,0)*rat(12*ep^2-10*ep+2,ep)
      +`Y'(1,n2,1,1,1,0,2,0,(1+n9),0,0,n12,0,0)*rat(2*ep-1,ep)
      +`Y'(1,n2,1,1,1,1,1,0,(1+n9),0,0,n12,0,0)*rat(8*ep+4*n2+2*(1+n9)+2*n12-6,1)
      +`Y'(1,n2,1,0,1,1,1,1,1+(1+n9),0,0,n12,0,0)*rat(2*(1+n9),1)
      +`Y'(1,n2,1,1,1,0,1,1,1+(1+n9),0,0,n12,0,0)*rat(2*(1+n9),1)
      +`Y'(1,n2,1,1,1,0,1,1,(1+n9),0,0,1+n12,0,0)*rat(n12,1)
      +`Y'(1,n2,1,1,1,1,1,0,1+(1+n9),-1,0,1+n12,0,0)*rat(-(1+n9)*n12,ep-1)
      +`Y'(1,n2,1,1,1,1,1,0,1+(1+n9),0,0,n12,0,0)*rat(-2*(1+n9),1)
      +`Y'(1,n2,1,1,1,1,1,0,2+(1+n9),-1,0,n12,0,0)*rat(-(1+n9)^2-(1+n9),ep-1)
      +`Y'(1,n2,1,1,1,1,1,0,(1+n9),0,0,1+n12,0,0)*rat(-n12,1)
      +`Y'(0,-1+n2,1,1,1,1,1,1,1+(1+n9),0,0,n12,0,0)*rat(-2*(1+n9),1)
      +`Y'(0,-1+n2,1,1,1,1,1,1,(1+n9),0,0,1+n12,0,0)*rat(-2*n12,1)
      +`Y'(0,n2,1,0,1,1,1,1,1+(1+n9),0,0,n12,0,0)*rat(2*(1+n9),1)
      +`Y'(0,n2,1,0,1,1,1,1,(1+n9),0,0,1+n12,0,0)*rat(4*n12,1)
      +`Y'(0,n2,1,1,1,0,1,1,(1+n9),0,0,1+n12,0,0)*rat(n12,1)
      +`Y'(0,n2,1,1,1,1,1,0,1+(1+n9),-1,0,1+n12,0,0)*rat(-(1+n9)*n12,ep-1)
      +`Y'(0,n2,1,1,1,1,1,0,(1+n9),-1,0,2+n12,0,0)*rat(-n12^2-n12,ep-1)
      +`Y'(0,n2,1,1,1,1,1,0,(1+n9),0,0,1+n12,0,0)*rat(-n12,1)
      +`Y'(1,-1+n2,0,1,1,1,1,1,1+(1+n9),0,0,n12,0,0)*rat(2*(1+n9),1)
      +`Y'(1,-1+n2,0,1,1,1,1,1,(1+n9),0,0,1+n12,0,0)*rat(-2*n12,1)
      +`Y'(1,n2,0,0,1,1,1,1,1+(1+n9),0,0,n12,0,0)*rat(-2*(1+n9),1)
      +`Y'(1,n2,0,1,1,0,1,1,(1+n9),0,0,1+n12,0,0)*rat(-n12,1)
      +`Y'(1,n2,0,1,1,1,1,0,1+(1+n9),-1,0,1+n12,0,0)*rat((1+n9)*n12,ep-1)
      +`Y'(1,n2,0,1,1,1,1,0,(1+n9),0,0,1+n12,0,0)*rat(n12,1)
      );
*--#] n9 : 
#endif
#if ( `il' >= 15 )
*--#[ n12 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,1,1,1,1,1,1,1,0,0,0,-1,0,0)
		 = -rat(8*ep-3,80*ep^2-72*ep+16)*(
      +`Z'(1,1,1,1,1,1,1,1,0,0,0,0,0,0)*rat(-36*ep^2+24*ep-4,8*ep-3)
      +`Y'(0,1,1,1,0,2,1,1,0,0,0,0,0,0)*rat(1152*ep^4-1840*ep^3+1060*ep^2-245*ep+16,128*ep^4-144*ep^3+52*ep^2-6*ep)
      +`Y'(0,1,1,1,1,1,1,1,0,0,0,0,0,0)*rat(-1152*ep^4+1600*ep^3-820*ep^2+181*ep-14,64*ep^3-72*ep^2+26*ep-3)
      +`Y'(0,1,1,1,1,1,1,2,0,-1,0,0,0,0)*rat(-1,2*ep-1)
      +`Y'(0,1,1,1,1,1,1,2,0,0,0,-1,0,0)*rat(-1,4*ep-2)
      +`Y'(0,1,1,1,1,1,1,2,0,0,0,0,0,0)*rat(96*ep^2-16*ep-15,32*ep^2-28*ep+6)
      +`Y'(0,1,1,1,1,2,0,1,0,0,0,0,0,0)*rat(-2560*ep^5+5824*ep^4-5120*ep^3+2192*ep^2-455*ep+36,512*ep^5-960*ep^4+640*ep^3
      -180*ep^2+18*ep)
      +`Y'(0,1,1,2,0,2,1,1,0,-1,0,0,0,0)*rat(-1,ep)
      +`Y'(0,1,1,2,1,1,1,1,0,-1,0,0,0,0)*rat(1,1)
      +`Y'(0,1,1,2,1,2,0,1,0,-1,0,0,0,0)*rat(1,2*ep)
      +`Y'(0,1,2,1,1,1,1,1,0,-1,0,0,0,0)*rat(1,1)
      +`Y'(0,1,2,1,1,1,1,1,0,0,0,0,0,-1)*rat(-1,1)
      +`Y'(0,2,1,0,1,1,1,1,0,0,0,0,0,0)*rat(-12*ep+5,32*ep^2-20*ep+3)
      +`Y'(0,2,1,1,-1,3,1,1,0,0,0,-1,0,0)*rat(24*ep^2-28*ep+8,32*ep^4+12*ep^3-17*ep^2+3*ep)
      +`Y'(0,2,1,1,0,1,1,1,0,0,0,0,0,0)*rat(-16*ep^2+8*ep,32*ep^2-20*ep+3)
      +`Y'(0,2,1,1,0,2,1,1,-1,0,0,0,0,0)*rat(-1,8*ep-3)
      +`Y'(0,2,1,1,0,2,1,1,0,-1,0,0,0,0)*rat(-8*ep+4,32*ep^3-20*ep^2+3*ep)
      +`Y'(0,2,1,1,0,2,1,1,0,0,0,-1,0,0)*rat(-24*ep^3+16*ep^2+6*ep-4,32*ep^3-20*ep^2+3*ep)
      +`Y'(0,2,1,1,0,2,1,1,0,0,0,0,-1,0)*rat(-16*ep+8,32*ep^2-20*ep+3)
      +`Y'(0,2,1,1,0,2,1,1,0,0,0,0,0,0)*rat(48*ep^3-164*ep^2+116*ep-29,256*ep^3-160*ep^2+24*ep)
      +`Y'(0,2,1,1,1,1,1,1,-1,0,0,0,0,0)*rat(10*ep-5,8*ep-3)
      +`Y'(0,2,1,1,1,1,1,1,0,-1,0,0,0,0)*rat(4*ep-3,32*ep^2-20*ep+3)
      +`Y'(0,2,1,1,1,1,1,1,0,0,0,0,0,0)*rat(4*ep^2-2*ep+3,32*ep^2-20*ep+3)
      +`Y'(0,2,1,1,1,1,1,2,0,-1,0,-1,0,0)*rat(1,4*ep-2)
      +`Y'(0,2,1,1,1,1,1,2,0,0,0,-1,0,0)*rat(-6,8*ep-3)
      +`Y'(0,2,1,1,1,2,0,1,-1,0,0,0,0,0)*rat(1,8*ep-3)
      +`Y'(0,2,1,1,1,2,0,1,0,-1,0,0,0,0)*rat(16*ep-8,128*ep^3-176*ep^2+72*ep-9)
      +`Y'(0,2,1,1,1,2,0,1,0,0,0,-1,0,0)*rat(64*ep^5-80*ep^4-40*ep^3+68*ep^2-6*ep-6,256*ep^5-224*ep^4-32*ep^3+54*ep^2-9*
      ep)
      +`Y'(0,2,1,1,1,2,0,1,0,0,0,0,-1,0)*rat(8*ep-4,32*ep^2-20*ep+3)
      +`Y'(0,2,1,1,1,2,0,1,0,0,0,0,0,0)*rat(-64*ep^4+352*ep^3-468*ep^2+144*ep+3,1024*ep^4-1408*ep^3+576*ep^2-72*ep)
      +`Y'(0,2,1,1,1,3,-1,1,0,0,0,-1,0,0)*rat(32*ep^2-48*ep+16,256*ep^4-224*ep^3-32*ep^2+54*ep-9)
      +`Y'(1,0,0,2,1,1,1,2,0,0,0,0,0,0)*rat(96*ep^2-64*ep+9,384*ep^3-432*ep^2+156*ep-18)
      +`Y'(1,0,1,1,1,1,1,1,0,0,0,0,0,0)*rat(-8*ep+4,1)
      +`Y'(1,0,1,1,1,2,0,1,0,0,0,0,0,0)*rat(-48*ep^3+76*ep^2-33*ep+4,16*ep^3-12*ep^2+2*ep)
      +`Y'(1,0,1,1,2,1,0,1,0,0,0,0,0,0)*rat(128*ep^3-112*ep^2+36*ep-5,128*ep^3-144*ep^2+52*ep-6)
      +`Y'(1,0,1,2,1,1,1,1,0,-1,0,0,0,0)*rat(1,1)
      +`Y'(1,0,1,2,1,2,0,1,0,-1,0,0,0,0)*rat(1,2*ep)
      +`Y'(1,1,0,1,1,1,1,1,0,0,0,0,0,0)*rat(160*ep^3-224*ep^2+97*ep-13,16*ep^2-14*ep+3)
      +`Y'(1,1,0,1,1,1,1,2,0,0,0,-1,0,0)*rat(1,4*ep-2)
      +`Y'(1,1,0,2,1,1,1,1,0,-1,0,0,0,0)*rat(-1,1)
      +`Y'(1,1,0,2,1,1,1,2,0,-1,0,-1,0,0)*rat(-1,4*ep-2)
      +`Y'(1,1,0,2,1,1,1,2,0,0,0,-1,0,0)*rat(-2,8*ep-3)
      +`Y'(1,1,1,0,1,1,1,1,0,0,0,0,0,0)*rat(320*ep^3-336*ep^2+118*ep-14,16*ep^2-14*ep+3)
      +`Y'(1,1,1,1,-1,3,1,1,0,0,0,-1,0,0)*rat(40*ep-20,32*ep^3+12*ep^2-17*ep+3)
      +`Y'(1,1,1,1,0,1,1,1,0,0,0,0,0,0)*rat(-640*ep^3+688*ep^2-276*ep+45,64*ep^2-40*ep+6)
      +`Y'(1,1,1,1,0,1,2,0,0,0,0,0,0,0)*rat(-128*ep^3+152*ep^2-56*ep+6,32*ep^3-20*ep^2+3*ep)
      +`Y'(1,1,1,1,0,1,2,1,0,0,0,0,0,0)*rat(96*ep^2-64*ep+9,32*ep^2-28*ep+6)
      +`Y'(1,1,1,1,0,2,1,0,0,0,0,0,0,0)*rat(-128*ep^3+152*ep^2-56*ep+6,32*ep^3-20*ep^2+3*ep)
      +`Y'(1,1,1,1,0,2,1,1,-1,0,0,0,0,0)*rat(1,8*ep-3)
      +`Y'(1,1,1,1,0,2,1,1,0,0,0,-1,0,0)*rat(-1056*ep^4+1632*ep^3-884*ep^2+199*ep-16,128*ep^4-144*ep^3+52*ep^2-6*ep)
      +`Y'(1,1,1,1,0,2,1,1,0,0,0,0,-1,0)*rat(-16*ep+8,32*ep^2-20*ep+3)
      +`Y'(1,1,1,1,0,2,1,1,0,0,0,0,0,0)*rat(336*ep^3-348*ep^2+68*ep+5,256*ep^3-160*ep^2+24*ep)
      +`Y'(1,1,1,1,0,2,2,1,0,-1,0,-1,0,0)*rat(-ep,4*ep^2-6*ep+2)
      +`Y'(1,1,1,1,0,2,2,1,0,0,0,-1,0,0)*rat(-2,8*ep-3)
      +`Y'(1,1,1,1,1,0,1,1,0,0,0,0,0,0)*rat(512*ep^4-880*ep^3+568*ep^2-152*ep+12,32*ep^3-20*ep^2+3*ep)
      +`Y'(1,1,1,1,1,0,1,2,0,0,0,0,0,0)*rat(-96*ep^2+64*ep-9,32*ep^2-28*ep+6)
      +`Y'(1,1,1,1,1,0,2,1,0,0,0,0,0,0)*rat(-96*ep^2+64*ep-9,32*ep^2-28*ep+6)
      +`Y'(1,1,1,1,1,1,0,1,0,0,0,0,0,0)*rat(-384*ep^4+544*ep^3-300*ep^2+73*ep-6,32*ep^3-20*ep^2+3*ep)
      +`Y'(1,1,1,1,1,1,0,2,0,0,0,-1,0,0)*rat(-160*ep^3+144*ep^2-37*ep+2,32*ep^3-28*ep^2+6*ep)
      +`Y'(1,1,1,1,1,1,1,0,0,0,0,0,0,0)*rat(-16*ep^2+4*ep+1,4*ep-2)
      +`Y'(1,1,1,1,1,1,2,0,0,-1,0,0,0,0)*rat(4*ep^2+2*ep-2,32*ep^3-52*ep^2+23*ep-3)
      +`Y'(1,1,1,1,1,2,0,1,-1,0,0,0,0,0)*rat(-1,8*ep-3)
      +`Y'(1,1,1,1,1,2,0,1,0,0,0,-1,0,0)*rat(2240*ep^5-2448*ep^4+112*ep^3+572*ep^2-182*ep+15,256*ep^5-224*ep^4-32*ep^3+54
      *ep^2-9*ep)
      +`Y'(1,1,1,1,1,2,0,1,0,0,0,0,-1,0)*rat(8*ep-4,32*ep^2-20*ep+3)
      +`Y'(1,1,1,1,1,2,0,1,0,0,0,0,0,0)*rat(-304*ep^3+260*ep^2-44*ep-3,256*ep^3-160*ep^2+24*ep)
      +`Y'(1,1,1,1,1,3,-1,1,0,0,0,-1,0,0)*rat(96*ep^2-48*ep,256*ep^4-224*ep^3-32*ep^2+54*ep-9)
      +`Y'(1,1,1,1,2,0,1,1,-1,0,0,0,0,0)*rat(32*ep^2-16*ep+1,32*ep^2-28*ep+6)
      +`Y'(1,1,1,1,2,0,1,1,0,-1,0,0,0,0)*rat(32*ep^2-16*ep+1,32*ep^2-28*ep+6)
      +`Y'(1,1,1,1,2,0,1,1,0,0,0,0,0,0)*rat(-128*ep^3+120*ep^2-44*ep+7,64*ep^3-72*ep^2+26*ep-3)
      +`Y'(1,1,1,1,2,0,1,2,0,-1,0,-1,0,0)*rat(-1,2*ep-1)
      +`Y'(1,1,1,1,2,0,1,2,0,0,0,-1,0,0)*rat(4,8*ep-3)
      +`Y'(1,1,1,1,2,1,0,1,-1,0,0,0,0,0)*rat(-1,8*ep-3)
      +`Y'(1,1,1,1,2,1,0,1,0,-1,0,0,0,0)*rat(128*ep^3-112*ep^2+36*ep-5,64*ep^3-104*ep^2+46*ep-6)
      +`Y'(1,1,1,1,2,1,0,1,0,0,0,-1,0,0)*rat(24*ep^2-16*ep+2,32*ep^2-20*ep+3)
      +`Y'(1,1,1,1,2,1,0,1,0,0,0,0,0,0)*rat(-608*ep^4+504*ep^3-124*ep^2+10*ep-1,512*ep^4-576*ep^3+208*ep^2-24*ep)
      +`Y'(1,1,1,1,2,1,0,2,0,-1,0,-1,0,0)*rat(1,ep-1)
      +`Y'(1,1,1,1,2,1,1,0,0,-1,0,0,0,0)*rat(4*ep^2+2*ep-2,32*ep^3-52*ep^2+23*ep-3)
      +`Y'(1,1,1,1,2,2,-1,1,0,0,0,-1,0,0)*rat(8*ep-4,128*ep^3-176*ep^2+72*ep-9)
      +`Y'(1,1,1,2,0,1,1,1,0,-1,0,0,0,0)*rat(1,1)
      +`Y'(1,1,1,2,0,1,2,0,0,-1,0,0,0,0)*rat(1,2*ep)
      +`Y'(1,1,1,2,0,2,1,0,0,-1,0,0,0,0)*rat(1,2*ep)
      +`Y'(1,1,1,2,0,2,1,1,0,-1,0,-1,0,0)*rat(1,ep)
      +`Y'(1,1,1,2,1,0,1,1,0,-1,0,0,0,0)*rat(-2*ep+1,ep)
      +`Y'(1,1,1,2,1,1,0,1,0,-1,0,0,0,0)*rat(2*ep-1,2*ep)
      +`Y'(1,1,1,2,1,2,0,1,0,-1,0,-1,0,0)*rat(-1,2*ep)
      +`Y'(1,2,0,1,1,1,1,1,-1,0,0,0,0,0)*rat(-2*ep+3,8*ep-3)
      +`Y'(1,2,0,1,1,1,1,1,0,-1,0,0,0,0)*rat(1,8*ep-3)
      +`Y'(1,2,0,1,1,1,1,1,0,0,0,0,0,0)*rat(-36*ep^2+2*ep+1,32*ep^2-20*ep+3)
      +`Y'(1,2,1,0,1,0,1,1,0,0,0,0,0,0)*rat(16*ep^2-24*ep+8,32*ep^2-20*ep+3)
      +`Y'(1,2,1,0,1,0,2,1,0,0,0,0,0,0)*rat(-6*ep+3,32*ep^2-20*ep+3)
      +`Y'(1,2,1,0,1,1,1,0,0,0,0,0,0,0)*rat(-16*ep^2+16*ep-4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,0,1,1,1,1,-1,0,0,0,0,0)*rat(-12*ep+2,8*ep-3)
      +`Y'(1,2,1,0,1,1,1,1,0,0,0,-1,0,0)*rat(-48*ep^2+32*ep-5,16*ep^2-14*ep+3)
      +`Y'(1,2,1,0,1,1,1,1,0,0,0,0,0,0)*rat(36*ep^2-2*ep+2,32*ep^2-20*ep+3)
      +`Y'(1,2,1,0,2,0,1,1,0,0,0,-1,0,0)*rat(16*ep^2-16*ep+4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,0,2,0,1,1,0,0,0,0,0,0)*rat(-16*ep^3-4*ep^2+4*ep+1,128*ep^3-80*ep^2+12*ep)
      +`Y'(1,2,1,1,-1,2,1,1,0,0,0,-1,0,0)*rat(-8*ep^2+8*ep-2,32*ep^3-20*ep^2+3*ep)
      +`Y'(1,2,1,1,-1,2,2,0,0,0,0,-1,0,0)*rat(-2*ep+1,8*ep^3+5*ep^2-3*ep)
      +`Y'(1,2,1,1,-1,3,1,0,0,0,0,-1,0,0)*rat(-4*ep+2,8*ep^3+5*ep^2-3*ep)
      +`Y'(1,2,1,1,-1,3,1,1,0,0,0,-2,0,0)*rat(-40*ep+20,32*ep^3+12*ep^2-17*ep+3)
      +`Y'(1,2,1,1,0,1,1,0,0,0,0,0,0,0)*rat(16*ep^2-16*ep+4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,0,1,1,1,-1,0,0,0,0,0)*rat(2*ep,8*ep-3)
      +`Y'(1,2,1,1,0,1,1,1,0,0,0,-1,0,0)*rat(-192*ep^5+352*ep^4+192*ep^3-452*ep^2+223*ep-37,128*ep^4-16*ep^3-92*ep^2+46*
      ep-6)
      +`Y'(1,2,1,1,0,1,1,1,0,0,0,0,0,0)*rat(96*ep^4+200*ep^3-196*ep^2+34*ep+1,256*ep^3-160*ep^2+24*ep)
      +`Y'(1,2,1,1,0,1,2,0,0,0,0,-1,0,0)*rat(16*ep^4+8*ep^3-4*ep^2-8*ep+3,32*ep^4+12*ep^3-17*ep^2+3*ep)
      +`Y'(1,2,1,1,0,1,2,0,0,0,0,0,-1,0)*rat(8*ep-4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,0,1,2,0,0,0,0,0,0,0)*rat(-16*ep^3+44*ep^2-12*ep-1,128*ep^3-80*ep^2+12*ep)
      +`Y'(1,2,1,1,0,1,2,1,0,-1,0,-1,0,0)*rat(1,4*ep-2)
      +`Y'(1,2,1,1,0,1,2,1,0,0,0,-1,0,0)*rat(-6,8*ep-3)
      +`Y'(1,2,1,1,0,2,1,0,0,0,0,-1,0,0)*rat(16*ep^4-8*ep^3-8*ep^2+2*ep+1,32*ep^4+12*ep^3-17*ep^2+3*ep)
      +`Y'(1,2,1,1,0,2,1,0,0,0,0,0,-1,0)*rat(8*ep-4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,0,2,1,0,0,0,0,0,0,0)*rat(-16*ep^3+44*ep^2-12*ep-1,128*ep^3-80*ep^2+12*ep)
      +`Y'(1,2,1,1,0,2,1,1,-1,0,0,-1,0,0)*rat(-1,8*ep-3)
      +`Y'(1,2,1,1,0,2,1,1,0,-1,0,-1,0,0)*rat(-16*ep+8,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,0,2,1,1,0,0,0,-2,0,0)*rat(40*ep^2-16*ep-2,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,0,2,1,1,0,0,0,-1,0,0)*rat(-80*ep^3+316*ep^2-108*ep-5,256*ep^3-160*ep^2+24*ep)
      +`Y'(1,2,1,1,0,2,2,0,0,-1,0,-1,0,0)*rat(4*ep-2,32*ep^3-52*ep^2+23*ep-3)
      +`Y'(1,2,1,1,1,0,1,1,-1,0,0,0,0,0)*rat(32*ep^2-16*ep+1,32*ep^2-28*ep+6)
      +`Y'(1,2,1,1,1,0,1,1,0,-1,0,0,0,0)*rat(-16*ep^2+24*ep-8,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,1,0,1,1,0,0,0,-1,0,0)*rat(-32*ep^4+32*ep^3-12*ep+4,32*ep^3-20*ep^2+3*ep)
      +`Y'(1,2,1,1,1,0,1,1,0,0,0,0,-1,0)*rat(-32*ep^2+48*ep-16,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,1,0,1,1,0,0,0,0,0,0)*rat(16*ep^4-76*ep^3+40*ep^2+3*ep-1,64*ep^3-40*ep^2+6*ep)
      +`Y'(1,2,1,1,1,0,1,2,0,-1,0,-1,0,0)*rat(-1,4*ep-2)
      +`Y'(1,2,1,1,1,0,1,2,0,0,0,-1,0,0)*rat(6,8*ep-3)
      +`Y'(1,2,1,1,1,0,2,1,0,-1,0,-1,0,0)*rat(-1,4*ep-2)
      +`Y'(1,2,1,1,1,0,2,1,0,0,0,-1,0,0)*rat(6,8*ep-3)
      +`Y'(1,2,1,1,1,1,0,1,-1,0,0,0,0,0)*rat(4*ep-2,8*ep-3)
      +`Y'(1,2,1,1,1,1,0,1,0,0,0,-1,0,0)*rat(256*ep^6-320*ep^5-96*ep^4+272*ep^3-112*ep^2+6,256*ep^5-224*ep^4-32*ep^3+54*
      ep^2-9*ep)
      +`Y'(1,2,1,1,1,1,0,1,0,0,0,0,-1,0)*rat(16*ep^2-24*ep+8,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,1,1,0,1,0,0,0,0,0,0)*rat(-16*ep^4+12*ep^3+4*ep^2-5*ep+1,64*ep^3-40*ep^2+6*ep)
      +`Y'(1,2,1,1,1,1,0,2,-1,0,0,-1,0,0)*rat(10*ep-3,16*ep^2-6*ep)
      +`Y'(1,2,1,1,1,1,0,2,0,-1,0,-1,0,0)*rat(32*ep^3-64*ep^2+24*ep,128*ep^4-304*ep^3+248*ep^2-81*ep+9)
      +`Y'(1,2,1,1,1,1,0,2,0,0,0,-2,0,0)*rat(-192*ep^5+112*ep^4+152*ep^3-76*ep^2-14*ep+6,256*ep^5-224*ep^4-32*ep^3+54*
      ep^2-9*ep)
      +`Y'(1,2,1,1,1,1,0,2,0,0,0,-1,0,0)*rat(48*ep^3-100*ep^2+20*ep+3,256*ep^3-160*ep^2+24*ep)
      +`Y'(1,2,1,1,1,1,1,0,-1,0,0,0,0,0)*rat(48*ep^2-24*ep+1,32*ep^2-28*ep+6)
      +`Y'(1,2,1,1,1,1,1,0,0,-1,0,0,0,0)*rat(16*ep^2-16*ep+4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,1,1,1,0,0,0,0,-1,0,0)*rat(32*ep^3-8*ep,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,1,1,1,0,0,0,0,0,0,0)*rat(-64*ep^4-352*ep^3+156*ep^2-48*ep-1,256*ep^3-160*ep^2+24*ep)
      +`Y'(1,2,1,1,1,2,-1,1,0,0,0,-1,0,0)*rat(-32*ep^2+48*ep-16,256*ep^4-224*ep^3-32*ep^2+54*ep-9)
      +`Y'(1,2,1,1,1,2,0,1,-1,0,0,-1,0,0)*rat(-1,2*ep)
      +`Y'(1,2,1,1,1,2,0,1,0,0,0,-1,0,0)*rat(-8*ep+4,32*ep^2-36*ep+9)
      +`Y'(1,2,1,1,1,3,-1,1,0,0,0,-2,0,0)*rat(-96*ep^2+48*ep,256*ep^4-224*ep^3-32*ep^2+54*ep-9)
      +`Y'(1,2,1,1,2,0,1,1,-2,0,0,0,0,0)*rat(-1,8*ep-3)
      +`Y'(1,2,1,1,2,0,1,1,-1,-1,0,0,0,0)*rat(-1,8*ep-3)
      +`Y'(1,2,1,1,2,0,1,1,-1,0,0,0,0,0)*rat(-6*ep+4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,2,0,1,1,0,-1,0,-1,0,0)*rat(-16*ep^2+16*ep-4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,2,0,1,1,0,-1,0,0,0,0)*rat(16*ep^3-28*ep^2+8*ep-1,128*ep^3-80*ep^2+12*ep)
      +`Y'(1,2,1,1,2,0,1,1,0,0,0,-1,0,0)*rat(16*ep^2-16*ep+4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,2,0,1,1,0,0,0,0,0,0)*rat(-16*ep^3+84*ep^2-32*ep+1,128*ep^3-80*ep^2+12*ep)
      +`Y'(1,2,1,1,2,1,-1,2,0,-1,0,-1,0,0)*rat(16*ep^2-8*ep,128*ep^4-304*ep^3+248*ep^2-81*ep+9)
      +`Y'(1,2,1,1,2,1,0,1,-1,-1,0,0,0,0)*rat(-2*ep+1,8*ep^2-11*ep+3)
      +`Y'(1,2,1,1,2,1,0,1,-1,0,0,0,0,0)*rat(10*ep-4,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,2,1,0,1,0,-1,0,-1,0,0)*rat(192*ep^4-432*ep^3+376*ep^2-148*ep+22,128*ep^4-304*ep^3+248*ep^2-81*ep+9)
      +`Y'(1,2,1,1,2,1,0,1,0,-1,0,0,0,0)*rat(-96*ep^4+184*ep^3-116*ep^2+30*ep-3,256*ep^4-416*ep^3+184*ep^2-24*ep)
      +`Y'(1,2,1,1,2,1,0,1,0,0,0,-1,0,0)*rat(-40*ep^2+32*ep-6,32*ep^2-20*ep+3)
      +`Y'(1,2,1,1,2,1,0,1,0,0,0,0,0,0)*rat(80*ep^3-140*ep^2+52*ep-3,256*ep^3-160*ep^2+24*ep)
      +`Y'(1,2,1,1,2,2,-1,1,0,0,0,-1,0,0)*rat(-8*ep+4,128*ep^3-176*ep^2+72*ep-9)
      +`Y'(1,3,1,0,1,1,1,1,-1,0,0,-1,0,0)*rat(4,8*ep-3)
      +`Y'(1,3,1,0,1,1,1,1,0,0,0,-1,0,0)*rat(-16*ep,32*ep^2-20*ep+3)
      +`Y'(1,3,1,1,-1,2,1,1,0,0,0,-2,0,0)*rat(-8*ep+4,32*ep^3-20*ep^2+3*ep)
      +`Y'(1,3,1,1,0,1,1,1,-1,0,0,-1,0,0)*rat(-2,8*ep-3)
      +`Y'(1,3,1,1,0,1,1,1,0,0,0,-2,0,0)*rat(16*ep^2-4,32*ep^2-20*ep+3)
      +`Y'(1,3,1,1,0,1,1,1,0,0,0,-1,0,0)*rat(-16*ep^3+12*ep^2+4*ep-1,128*ep^3-80*ep^2+12*ep)
      +`Y'(1,3,1,1,0,1,2,0,0,-1,0,-1,0,0)*rat(16*ep^2-4,32*ep^3-52*ep^2+23*ep-3)
      +`Y'(1,3,1,1,0,1,2,0,0,-1,0,0,0,0)*rat(-16*ep^3-20*ep^2+4*ep-1,128*ep^4-208*ep^3+92*ep^2-12*ep)
      +`Y'(1,3,1,1,1,0,1,1,-2,0,0,0,0,0)*rat(-2,8*ep-3)
      +`Y'(1,3,1,1,1,0,1,1,-1,0,0,0,0,0)*rat(-4*ep,32*ep^2-20*ep+3)
      +`Y'(1,3,1,1,1,1,1,0,-2,0,0,0,0,0)*rat(-2,8*ep-3)
      +`Y'(1,3,1,1,1,1,1,0,-1,0,0,0,0,0)*rat(8*ep,32*ep^2-20*ep+3)
      +`Y'(1,3,1,1,1,1,1,0,0,-1,0,-1,0,0)*rat(-48*ep^3+16*ep^2+12*ep-4,32*ep^3-52*ep^2+23*ep-3)
      +`Y'(1,3,1,1,1,1,1,0,0,-1,0,0,0,0)*rat(48*ep^3+60*ep^2-12*ep+3,128*ep^3-208*ep^2+92*ep-12)
      +`Y'(2,1,0,1,1,1,1,1,0,0,0,0,0,-1)*rat(-1,1)
      );
`SKIP2'
id,ifmatch->bebe1,`Z'(1,1,1,1,1,1,1,1,0,0,0,n12?{<-1},0,0)
		 = -rat(16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3,
       80*ep^2*(1+n12)+36*ep*(1+n12)^2-72*ep*(1+n12)+4*(1+n12)^3-16*(1+n12)^2+16*(1+n12))*(
      +`Z'(1,1,1,1,1,1,1,1,-1,0,0,1+(1+n12),0,0)*rat(-12*ep^2*(1+n12)-10*ep*(1+n12)^2+3*ep*(1+n12)-2*(1+n12)^3+(1+n12)^2,16*ep^2+4*ep*(1+n12)-14*ep
      -2*(1+n12)+3)
      +`Z'(1,1,1,1,1,1,1,1,0,-1,0,1+(1+n12),0,0)*rat((1+n12),8*ep+2*(1+n12)-3)
      +`Z'(1,1,1,1,1,1,1,1,0,0,0,1+(1+n12),-1,0)*rat((1+n12),1)
      +`Z'(1,1,1,1,1,1,1,1,0,0,0,(1+n12),0,0)*rat(-28*ep^2*(1+n12)-16*ep*(1+n12)^2+22*ep*(1+n12)-2*(1+n12)^3+7*(1+n12)^2-4*(1+n12),16*ep^2+4*ep*(1+n12)
      -14*ep-2*(1+n12)+3)
      +`Z'(1,1,1,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(-(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,1,0,2,1,1,-1,0,0,(1+n12),0,0)*rat((1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,1,0,2,1,1,0,0,0,-1+(1+n12),0,0)*rat(40*ep^2+ep*(1+n12)-31*ep-2*(1+n12)^2+6,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,1,1,0,2,-1,0,0,(1+n12),0,0)*rat(-(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,1,1,1,0,2,0,0,0,-1+(1+n12),0,0)*rat(-24*ep^2-ep*(1+n12)+17*ep+(1+n12)^2-3,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,2,0,1,1,-2,0,0,1+(1+n12),0,0)*rat(-4*ep*(1+n12)-2*(1+n12)^2+(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(1,1,1,1,2,0,1,1,-1,-1,0,1+(1+n12),0,0)*rat(-4*ep*(1+n12)-2*(1+n12)^2+(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(1,1,1,1,2,0,1,1,-1,0,0,(1+n12),0,0)*rat(32*ep^2+23*ep*(1+n12)-28*ep+4*(1+n12)^2-10*(1+n12)+6,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,2,0,1,1,0,-1,0,(1+n12),0,0)*rat(-(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,1,2,1,0,1,-1,-1,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)-(1+n12),8*ep^2+2*ep*(1+n12)-11*ep-2*(1+n12)+3)
      +`Y'(1,1,1,1,2,1,0,1,-1,0,0,(1+n12),0,0)*rat(-3*ep-(1+n12)+1,ep)
      +`Y'(1,1,1,1,2,1,0,1,0,-1,0,(1+n12),0,0)*rat(-2*ep*(1+n12)+(1+n12),8*ep^2+2*ep*(1+n12)-11*ep-2*(1+n12)+3)
      +`Y'(1,2,1,0,1,1,1,1,-1,0,0,(1+n12),0,0)*rat(4*ep*(1+n12)+2*(1+n12)^2-(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,2,1,0,1,1,1,1,0,0,0,-1+(1+n12),0,0)*rat(20*ep*(1+n12)+4*(1+n12)^2-8*(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,2,1,1,0,1,1,1,-1,0,0,(1+n12),0,0)*rat((1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,2,1,1,0,1,1,1,0,0,0,-1+(1+n12),0,0)*rat(8*ep+(1+n12)-3,8*ep+2*(1+n12)-3)
      +`Y'(1,2,1,1,0,1,2,0,-1,-1,0,1+(1+n12),0,0)*rat((1+n12),4*ep^2-6*ep+2)
      +`Y'(1,2,1,1,0,1,2,0,0,-1,0,(1+n12),0,0)*rat(2*ep+(1+n12)-1,2*ep^2-3*ep+1)
      +`Y'(1,2,1,1,1,0,1,1,-2,0,0,1+(1+n12),0,0)*rat(-4*ep*(1+n12)-2*(1+n12)^2+(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(1,2,1,1,1,0,1,1,-1,0,0,(1+n12),0,0)*rat(-(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,2,1,1,1,0,2,0,-1,-1,0,1+(1+n12),0,0)*rat(-(1+n12),4*ep^2-6*ep+2)
      +`Y'(1,2,1,1,1,0,2,0,0,-1,0,(1+n12),0,0)*rat(-2*ep-(1+n12)+1,2*ep^2-3*ep+1)
      +`Y'(1,2,1,1,1,1,1,0,-2,0,0,1+(1+n12),0,0)*rat(-4*ep*(1+n12)-2*(1+n12)^2+(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(1,2,1,1,1,1,1,0,-1,-1,0,1+(1+n12),0,0)*rat(-3*ep*(1+n12)-(1+n12)^2+(1+n12),4*ep^2-6*ep+2)
      +`Y'(1,2,1,1,1,1,1,0,-1,0,0,(1+n12),0,0)*rat(-10*ep*(1+n12)-2*(1+n12)^2+4*(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,2,1,1,1,1,1,0,0,-1,0,(1+n12),0,0)*rat(-6*ep^2-5*ep*(1+n12)+5*ep-(1+n12)^2+2*(1+n12)-1,2*ep^2-3*ep+1)
      +`Y'(1,1,1,1,2,0,1,1,-1,0,0,1+(1+n12),0,0)*rat(4*ep*(1+n12)+2*(1+n12)^2-(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,1,1,1,2,0,1,1,0,0,0,(1+n12),0,0)*rat(-16*ep^2-18*ep*(1+n12)+22*ep-4*(1+n12)^2+10*(1+n12)-6,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,2,1,0,1,-1,0,0,1+(1+n12),0,0)*rat(-(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,1,2,1,0,1,0,0,0,(1+n12),0,0)*rat(8*ep^2+11*ep*(1+n12)-11*ep+2*(1+n12)^2-5*(1+n12)+3,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(0,1,1,1,0,2,1,1,-1,0,0,1+(1+n12),0,0)*rat((1+n12),8*ep+2*(1+n12)-3)
      +`Y'(0,1,1,1,0,2,1,1,0,0,0,(1+n12),0,0)*rat(-24*ep^2+ep*(1+n12)+25*ep+2*(1+n12)^2-6,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(0,1,1,1,1,1,1,1,-1,-1,0,2+(1+n12),0,0)*rat(-2*(1+n12)^2-2*(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(0,1,1,1,1,1,1,1,-1,0,0,1+(1+n12),0,0)*rat(-12*ep^2*(1+n12)+6*ep*(1+n12)^2+13*ep*(1+n12)+2*(1+n12)^3-(1+n12)^2-3*(1+n12),16*ep^2+4*ep*(1+n12)
      -14*ep-2*(1+n12)+3)
      +`Y'(0,1,1,1,1,1,1,1,0,-1,0,1+(1+n12),0,0)*rat(-16*ep*(1+n12)-2*(1+n12)^2+7*(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(0,1,1,1,1,1,1,1,0,0,0,1+(1+n12),-1,0)*rat((1+n12),1)
      +`Y'(0,1,1,1,1,1,1,1,0,0,0,(1+n12),0,0)*rat(4*ep^2*(1+n12)-8*ep*(1+n12)^2+6*ep*(1+n12)-2*(1+n12)^3+5*(1+n12)^2-3*(1+n12),16*ep^2+4*ep*(1+n12)-14*
      ep-2*(1+n12)+3)
      +`Y'(0,1,1,1,1,2,0,1,-1,0,0,1+(1+n12),0,0)*rat(-(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(0,1,1,1,1,2,0,1,0,0,0,(1+n12),0,0)*rat(8*ep^2-ep*(1+n12)-11*ep-(1+n12)^2+3,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,0,1,1,0,2,1,1,-1,0,0,1+(1+n12),0,0)*rat(-(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,0,1,1,0,2,1,1,0,0,0,(1+n12),0,0)*rat(-40*ep^2-ep*(1+n12)+31*ep+2*(1+n12)^2-6,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,0,1,1,1,1,1,1,0,0,0,(1+n12),0,0)*rat(64*ep^2+24*ep*(1+n12)-56*ep+2*(1+n12)^2-10*(1+n12)+12,8*ep+2*(1+n12)-3)
      +`Y'(1,0,1,1,1,2,0,1,-1,0,0,1+(1+n12),0,0)*rat((1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,0,1,1,1,2,0,1,0,0,0,(1+n12),0,0)*rat(24*ep^2+ep*(1+n12)-17*ep-(1+n12)^2+3,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,0,1,1,2,0,1,1,0,0,0,(1+n12),0,0)*rat(-4*ep-2*(1+n12)+2,ep)
      +`Y'(1,0,1,1,2,1,0,1,-1,0,0,1+(1+n12),0,0)*rat((1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,0,1,1,2,1,0,1,0,0,0,(1+n12),0,0)*rat(24*ep^2+13*ep*(1+n12)-17*ep+2*(1+n12)^2-5*(1+n12)+3,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,0,1,1,1,1,1,-1,0,0,1+(1+n12),0,0)*rat(12*ep^2*(1+n12)+10*ep*(1+n12)^2-3*ep*(1+n12)+2*(1+n12)^3-(1+n12)^2,16*ep^2+4*ep*(1+n12)-14*ep-
      2*(1+n12)+3)
      +`Y'(1,1,0,1,1,1,1,1,0,-1,0,1+(1+n12),0,0)*rat(-(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,0,1,1,1,1,1,0,0,0,1+(1+n12),-1,0)*rat(-(1+n12),1)
      +`Y'(1,1,0,1,1,1,1,1,0,0,0,(1+n12),0,0)*rat(76*ep^2*(1+n12)+44*ep*(1+n12)^2-58*ep*(1+n12)+6*(1+n12)^3-17*(1+n12)^2+11*(1+n12),16*ep^2+4*ep*
      (1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,1,1,0,0,1,2,1,0,0,0,(1+n12),0,0)*rat(-2*ep-(1+n12)+1,ep)
      +`Y'(1,1,1,0,1,0,2,1,0,0,0,(1+n12),0,0)*rat(2*ep+(1+n12)-1,ep)
      +`Y'(1,1,1,0,1,1,1,1,-1,0,0,1+(1+n12),0,0)*rat(-24*ep^2*(1+n12)-16*ep*(1+n12)^2+14*ep*(1+n12)-2*(1+n12)^3+5*(1+n12)^2-2*(1+n12),16*ep^2+4*ep*
      (1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,1,1,0,1,1,1,1,0,0,0,(1+n12),0,0)*rat(-64*ep^3-120*ep^2*(1+n12)+88*ep^2-40*ep*(1+n12)^2+88*ep*(1+n12)-40*ep-4*(1+n12)^3+14*
      (1+n12)^2-16*(1+n12)+6,16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,1,1,0,2,0,1,1,0,0,0,(1+n12),0,0)*rat(2*ep+(1+n12)-1,ep)
      +`Y'(1,1,1,0,2,1,0,1,-1,0,0,1+(1+n12),0,0)*rat(-(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,0,2,1,0,1,0,0,0,(1+n12),0,0)*rat(-8*ep-(1+n12)+3,8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,1,0,1,1,1,-1,0,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)-(1+n12)^2+(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,1,0,1,1,1,0,0,0,(1+n12),0,0)*rat(-48*ep^3-58*ep^2*(1+n12)+74*ep^2-11*ep*(1+n12)^2+48*ep*(1+n12)-37*ep+4*(1+n12)^2-10*(1+n12)+6
      ,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,0,1,2,0,0,0,0,(1+n12),0,0)*rat(16*ep^2-14*ep-(1+n12)^2+3,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,0,2,1,0,0,0,0,(1+n12),0,0)*rat(16*ep^2-14*ep-(1+n12)^2+3,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,1,0,1,1,0,0,0,(1+n12),0,0)*rat(-32*ep^3+24*ep^2*(1+n12)+60*ep^2+8*ep*(1+n12)^2-10*ep*(1+n12)-34*ep-2*(1+n12)^2+6,8*ep^2+
      2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,1,0,2,0,0,0,0,(1+n12),0,0)*rat(-16*ep^2+14*ep+(1+n12)^2-3,8*ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,1,1,0,1,-1,0,0,1+(1+n12),0,0)*rat(-4*ep*(1+n12)+2*(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,1,1,1,0,1,0,0,0,(1+n12),0,0)*rat(32*ep^3+20*ep^2*(1+n12)-60*ep^2+2*ep*(1+n12)^2-20*ep*(1+n12)+34*ep-(1+n12)^2+5*(1+n12)-6,8*
      ep^2+2*ep*(1+n12)-3*ep)
      +`Y'(1,1,1,1,1,1,1,0,-1,0,0,1+(1+n12),0,0)*rat(40*ep^2*(1+n12)+18*ep*(1+n12)^2-28*ep*(1+n12)+2*(1+n12)^3-6*(1+n12)^2+5*(1+n12),16*ep^2+4*ep*
      (1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,1,1,1,1,1,1,0,0,0,0,(1+n12),0,0)*rat(64*ep^3+104*ep^2*(1+n12)-72*ep^2+38*ep*(1+n12)^2-72*ep*(1+n12)+26*ep+4*(1+n12)^3-13*(1+n12)^2
      +12*(1+n12)-3,16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(0,1,1,1,1,1,1,1,-1,0,0,2+(1+n12),0,0)*rat(2*(1+n12)^2+2*(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(0,1,1,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(16*ep*(1+n12)+2*(1+n12)^2-7*(1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,0,1,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(8*ep*(1+n12)^2+6*ep*(1+n12)+2*(1+n12)^3-(1+n12)^2-2*(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,1,0,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat((1+n12),8*ep+2*(1+n12)-3)
      +`Y'(1,1,1,0,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(-6*ep*(1+n12)-2*(1+n12)^2+2*(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,1,1,1,1,0,1,1,-1,0,0,2+(1+n12),0,0)*rat(-4*ep*(1+n12)^2-4*ep*(1+n12)-2*(1+n12)^3-(1+n12)^2+(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(1,1,1,1,1,0,1,1,0,0,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)^2+4*ep*(1+n12)+3*(1+n12)^2-(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,1,1,1,1,1,1,0,-1,0,0,2+(1+n12),0,0)*rat(-4*ep*(1+n12)^2-4*ep*(1+n12)-2*(1+n12)^3-(1+n12)^2+(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(1,1,1,1,1,1,1,0,0,-1,0,2+(1+n12),0,0)*rat((1+n12)^2+(1+n12),4*ep^2-6*ep+2)
      +`Y'(1,1,1,1,1,1,1,0,0,0,0,1+(1+n12),0,0)*rat(-10*ep*(1+n12)^2+4*ep*(1+n12)-2*(1+n12)^3+6*(1+n12)^2-(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+
      3)
      +`Y'(0,0,1,1,1,1,1,1,-1,0,0,2+(1+n12),0,0)*rat(-12*ep*(1+n12)^2-12*ep*(1+n12)-4*(1+n12)^3+4*(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(0,0,1,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(-32*ep^2*(1+n12)-36*ep*(1+n12)^2+34*ep*(1+n12)-6*(1+n12)^3+17*(1+n12)^2-8*(1+n12),16*ep^2+4*ep*
      (1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(0,1,1,0,1,1,1,1,-1,0,0,2+(1+n12),0,0)*rat(12*ep*(1+n12)^2+12*ep*(1+n12)+4*(1+n12)^3-4*(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(0,1,1,0,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(32*ep^2*(1+n12)+44*ep*(1+n12)^2-34*ep*(1+n12)+8*(1+n12)^3-20*(1+n12)^2+8*(1+n12),16*ep^2+4*ep*
      (1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(0,1,1,1,1,0,1,1,-1,0,0,2+(1+n12),0,0)*rat(-4*ep*(1+n12)^2-4*ep*(1+n12)-2*(1+n12)^3-(1+n12)^2+(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(0,1,1,1,1,0,1,1,0,0,0,1+(1+n12),0,0)*rat(-(1+n12)^2,8*ep+2*(1+n12)-3)
      +`Y'(0,1,1,1,1,1,1,0,-1,-1,0,3+(1+n12),0,0)*rat(-(1+n12)^3-3*(1+n12)^2-2*(1+n12),4*ep^2-6*ep+2)
      +`Y'(0,1,1,1,1,1,1,0,-1,0,0,2+(1+n12),0,0)*rat(-4*ep*(1+n12)^2-4*ep*(1+n12)-2*(1+n12)^3-(1+n12)^2+(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(0,1,1,1,1,1,1,0,0,-1,0,2+(1+n12),0,0)*rat(-4*ep*(1+n12)^2-4*ep*(1+n12)-2*(1+n12)^3+(1+n12)^2+3*(1+n12),4*ep^2-6*ep+2)
      +`Y'(0,1,1,1,1,1,1,0,0,0,0,1+(1+n12),0,0)*rat(-10*ep*(1+n12)^2-2*(1+n12)^3+4*(1+n12)^2,16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,0,0,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(-8*ep*(1+n12)^2-6*ep*(1+n12)-2*(1+n12)^3+(1+n12)^2+2*(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3
      )
      +`Y'(1,1,0,0,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(6*ep*(1+n12)+2*(1+n12)^2-2*(1+n12),16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      +`Y'(1,1,0,1,1,0,1,1,-1,0,0,2+(1+n12),0,0)*rat(4*ep*(1+n12)^2+4*ep*(1+n12)+2*(1+n12)^3+(1+n12)^2-(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(1,1,0,1,1,0,1,1,0,0,0,1+(1+n12),0,0)*rat((1+n12)^2,8*ep+2*(1+n12)-3)
      +`Y'(1,1,0,1,1,1,1,0,-1,0,0,2+(1+n12),0,0)*rat(4*ep*(1+n12)^2+4*ep*(1+n12)+2*(1+n12)^3+(1+n12)^2-(1+n12),32*ep^2+8*ep*(1+n12)-28*ep-4*(1+n12)+6)
      +`Y'(1,1,0,1,1,1,1,0,0,-1,0,2+(1+n12),0,0)*rat(-(1+n12)^2-(1+n12),4*ep^2-6*ep+2)
      +`Y'(1,1,0,1,1,1,1,0,0,0,0,1+(1+n12),0,0)*rat(10*ep*(1+n12)^2+2*(1+n12)^3-4*(1+n12)^2,16*ep^2+4*ep*(1+n12)-14*ep-2*(1+n12)+3)
      );

`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?{>1},1,1,1,1,1,1,0,0,0,n12?neg_,0,0)
			 = -rat(2*ep-1,4*ep*n2-4*ep+n2^2+n2*(1+n12)-3*n2-(1+n12)+2)*(
      +`Z'(1,n2,1,1,1,1,1,1,0,0,0,(1+n12),0,0)*rat(-12*ep^2*n2+12*ep^2-7*ep*n2^2-7*ep*n2*(1+n12)+21*ep*n2+7*ep*(1+n12)-14*ep-n2^3-2
      *n2^2*(1+n12)+5*n2^2-n2*(1+n12)^2+6*n2*(1+n12)-8*n2+(1+n12)^2-4*(1+n12)+4,8*ep^2+2*ep*n2+2*ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Z'(1,-1+n2,1,1,1,1,1,1,0,-1,0,1+(1+n12),0,0)*rat(2*ep*(1+n12),2*ep-1)
      +`Z'(1,-1+n2,1,1,1,1,1,1,0,0,0,1+(1+n12),-1,0)*rat(-(1+n12),1)
      +`Z'(1,-1+n2,1,1,1,1,1,1,0,0,0,(1+n12),0,0)*rat(4*ep+n2+(1+n12)-2,1)
      +`Z'(1,-1+n2,1,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12),2*ep-1)
      +`Z'(1,-2+n2,1,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat((1+n12),2*ep-1)
      +`Y'(1,1+n2,1,0,1,1,1,1,0,0,0,-1+(1+n12),0,0)*rat(4*ep*n2^2-4*ep*n2+n2^3+n2^2*(1+n12)-3*n2^2-n2*(1+n12)+2*n2,8*ep^2+2*ep*n2+
      2*ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Y'(1,1+n2,1,1,-1,2,1,1,0,0,0,-1+(1+n12),0,0)*rat(n2*(1+n12),4*ep^2-2*ep)
      +`Y'(1,1+n2,1,1,0,1,1,1,0,0,0,-1+(1+n12),0,0)*rat(8*ep^2*n2+2*ep*n2^2+2*ep*n2*(1+n12)-6*ep*n2-n2^2+n2,16*ep^2+4*ep*n2+4*
      ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,1+n2,1,1,0,1,2,0,0,-1,0,(1+n12),0,0)*rat(2*ep*n2+n2^2-n2,4*ep^2-6*ep+2)
      +`Y'(1,1+n2,1,1,0,2,0,1,0,0,0,-1+(1+n12),0,0)*rat(-n2*(1+n12),8*ep^2-10*ep+3)
      +`Y'(1,1+n2,1,1,1,0,1,1,-1,0,0,(1+n12),0,0)*rat(-4*ep*n2^2+4*ep*n2-n2^3-n2^2*(1+n12)+3*n2^2+n2*(1+n12)-2*n2,16*ep^2+4*ep*n2+
      4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,1+n2,1,1,1,0,2,0,0,-1,0,(1+n12),0,0)*rat(-2*ep*n2-n2^2+n2,4*ep^2-6*ep+2)
      +`Y'(1,1+n2,1,1,1,1,1,0,-1,0,0,(1+n12),0,0)*rat(-4*ep*n2^2+4*ep*n2-n2^3-n2^2*(1+n12)+3*n2^2+n2*(1+n12)-2*n2,16*ep^2+4*ep*n2+
      4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,1+n2,1,1,1,1,1,0,0,-1,0,(1+n12),0,0)*rat(-6*ep^2*n2-5*ep*n2^2-2*ep*n2*(1+n12)+7*ep*n2-n2^3-n2^2*(1+n12)+3*n2^2+n2*(1+n12)
      -2*n2,4*ep^2-6*ep+2)
      +`Y'(1,n2,1,1,-1,3,1,1,0,0,0,-1+(1+n12),0,0)*rat(5*ep*(1+n12)+2*n2*(1+n12)+2*(1+n12)^2-2*(1+n12),2*ep^3+ep^2-ep)
      +`Y'(1,n2,1,1,0,2,1,1,0,-1,0,(1+n12),0,0)*rat(2*ep*(1+n12)+n2*(1+n12)+(1+n12)^2-(1+n12),2*ep^2-ep)
      +`Y'(1,n2,1,1,0,2,1,1,0,0,0,-1+(1+n12),0,0)*rat(40*ep^3+10*ep^2*n2+26*ep^2*(1+n12)-46*ep^2+4*ep*n2*(1+n12)-5*ep*n2+4*ep*
      (1+n12)^2-24*ep*(1+n12)+13*ep-2*n2*(1+n12)-2*(1+n12)^2+6*(1+n12),16*ep^3+4*ep^2*n2+4*ep^2*(1+n12)-20*ep^2-2*ep*n2-2*ep*(1+n12)+6*ep)
      +`Y'(1,n2,1,1,0,2,2,0,0,-1,0,(1+n12),0,0)*rat(-(1+n12),4*ep^2-6*ep+2)
      +`Y'(1,n2,1,1,0,3,0,1,0,0,0,-1+(1+n12),0,0)*rat(-88*ep^3*(1+n12)-32*ep^2*n2*(1+n12)-32*ep^2*(1+n12)^2-12*ep*n2*(1+n12)-12*ep*(1+n12)^2+
      36*ep*(1+n12)+6*n2*(1+n12)+6*(1+n12)^2-3*(1+n12),32*ep^5+8*ep^4-32*ep^3-2*ep^2+6*ep)
      +`Y'(1,n2,1,1,1,1,0,2,0,-1,0,(1+n12),0,0)*rat(-4*ep^2*(1+n12)-2*ep*n2*(1+n12)-2*ep*(1+n12)^2+8*ep*(1+n12)+2*n2*(1+n12)+2*(1+n12)^2-2*(1+n12),8*
      ep^3-18*ep^2+13*ep-3)
      +`Y'(1,n2,1,1,1,1,0,2,0,0,0,-1+(1+n12),0,0)*rat(-192*ep^4-48*ep^3*n2-112*ep^3*(1+n12)+256*ep^3-16*ep^2*n2*(1+n12)+36*ep^2*n2
      -16*ep^2*(1+n12)^2+44*ep^2*(1+n12)-36*ep^2-24*ep*n2*(1+n12)+12*ep*n2-24*ep*(1+n12)^2+78*ep*(1+n12)-64*ep-4*n2^2*(1+n12)-8*n2*(1+n12)^2+22*n2*
      (1+n12)-9*n2-4*(1+n12)^3+22*(1+n12)^2-36*(1+n12)+21,128*ep^4+32*ep^3*n2+32*ep^3*(1+n12)-192*ep^3-24*ep^2*n2-24*ep^2*(1+n12)+40*ep^2-8*ep*
      n2-8*ep*(1+n12)+48*ep+6*n2+6*(1+n12)-18)
      +`Y'(1,n2,1,1,1,3,-1,1,0,0,0,-1+(1+n12),0,0)*rat(12*ep*(1+n12)+4*n2*(1+n12)+4*(1+n12)^2-4*(1+n12),16*ep^3-12*ep^2-4*ep+3)
      +`Y'(1,n2,1,1,2,0,1,1,-1,0,0,(1+n12),0,0)*rat(-4*ep*n2+4*ep-n2^2-n2*(1+n12)+3*n2+(1+n12)-2,16*ep^2+4*ep*n2+4*ep*(1+n12)-20*ep-2*
      n2-2*(1+n12)+6)
      +`Y'(1,n2,1,1,2,0,1,1,0,-1,0,(1+n12),0,0)*rat(-16*ep^2-8*ep*n2-4*ep*(1+n12)+24*ep-n2^2-n2*(1+n12)+5*n2+3*(1+n12)-8,16*ep^2+4*ep*
      n2+4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,n2,1,1,2,1,-1,2,0,-1,0,(1+n12),0,0)*rat(-2*ep*(1+n12),8*ep^3-18*ep^2+13*ep-3)
      +`Y'(1,n2,1,1,2,1,0,1,0,-1,0,(1+n12),0,0)*rat(96*ep^3+24*ep^2*n2+24*ep^2*(1+n12)-208*ep^2-38*ep*n2-42*ep*(1+n12)+154*ep-2*n2
      *(1+n12)+15*n2-2*(1+n12)^2+18*(1+n12)-39,32*ep^3+8*ep^2*n2+8*ep^2*(1+n12)-80*ep^2-14*ep*n2-14*ep*(1+n12)+66*ep+6*n2+6*(1+n12)-18)
      +`Y'(1,n2,1,1,0,2,1,1,0,0,0,(1+n12),0,0)*rat(-2*ep*(1+n12)-n2*(1+n12)-(1+n12)^2+(1+n12),2*ep^2-ep)
      +`Y'(1,n2,1,1,1,2,0,1,0,0,0,(1+n12),0,0)*rat(4*ep*(1+n12)+2*n2*(1+n12)+2*(1+n12)^2-3*(1+n12),8*ep^2-10*ep+3)
      +`Y'(1,n2,1,1,2,0,1,1,0,0,0,(1+n12),0,0)*rat(8*ep^2+6*ep*n2+2*ep*(1+n12)-14*ep+n2^2+n2*(1+n12)-4*n2-2*(1+n12)+5,8*ep^2+2*ep*n2+2*
      ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Y'(1,n2,1,1,2,1,0,1,0,0,0,(1+n12),0,0)*rat(-40*ep^2-10*ep*n2-10*ep*(1+n12)+46*ep+5*n2+4*(1+n12)-13,16*ep^2+4*ep*n2+4*ep*(1+n12)
      -20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,n2,1,1,2,2,-1,1,0,0,0,(1+n12),0,0)*rat((1+n12),8*ep^2-10*ep+3)
      +`Y'(0,n2,1,1,-1,3,1,1,0,0,0,(1+n12),0,0)*rat(-3*ep*(1+n12)-2*n2*(1+n12)-2*(1+n12)^2+4*(1+n12),2*ep^3+ep^2-ep)
      +`Y'(0,n2,1,1,0,2,1,1,0,-1,0,1+(1+n12),0,0)*rat(-(1+n12)^2,2*ep^2-ep)
      +`Y'(0,n2,1,1,0,2,1,1,0,0,0,1+(1+n12),-1,0)*rat(2*(1+n12),2*ep-1)
      +`Y'(0,n2,1,1,0,2,1,1,0,0,0,(1+n12),0,0)*rat(-24*ep^2-6*ep*n2-22*ep*(1+n12)+34*ep-4*n2*(1+n12)+3*n2-4*(1+n12)^2+16*(1+n12)-11,16*ep^2
      +4*ep*n2+4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(0,n2,1,1,0,3,0,1,0,0,0,(1+n12),0,0)*rat(40*ep^3*(1+n12)+32*ep^2*n2*(1+n12)+32*ep^2*(1+n12)^2-48*ep^2*(1+n12)+12*ep*n2*(1+n12)+12*ep*
      (1+n12)^2-24*ep*(1+n12)-6*n2*(1+n12)-6*(1+n12)^2+15*(1+n12),32*ep^5+8*ep^4-32*ep^3-2*ep^2+6*ep)
      +`Y'(0,n2,1,1,1,1,1,1,0,-1,0,1+(1+n12),0,0)*rat(n2*(1+n12)-(1+n12),8*ep^2+2*ep*n2+2*ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Y'(0,n2,1,1,1,1,1,1,0,0,0,(1+n12),0,0)*rat(4*ep^2*n2-4*ep^2-3*ep*n2^2+ep*n2*(1+n12)+5*ep*n2-ep*(1+n12)-2*ep-n2^3-n2^2*(1+n12)+4
      *n2^2+2*n2*(1+n12)-5*n2-(1+n12)+2,8*ep^2+2*ep*n2+2*ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Y'(0,n2,1,1,1,2,0,1,0,-1,0,1+(1+n12),0,0)*rat(2*(1+n12)^2,8*ep^2-10*ep+3)
      +`Y'(0,n2,1,1,1,2,0,1,0,0,0,1+(1+n12),-1,0)*rat(-(1+n12),2*ep-1)
      +`Y'(0,n2,1,1,1,2,0,1,0,0,0,(1+n12),0,0)*rat(64*ep^4+16*ep^3*n2+80*ep^3*(1+n12)-128*ep^3+16*ep^2*n2*(1+n12)-12*ep^2*n2+16*
      ep^2*(1+n12)^2-68*ep^2*(1+n12)+44*ep^2+16*ep*n2*(1+n12)-4*ep*n2+16*ep*(1+n12)^2-58*ep*(1+n12)+32*ep+4*n2^2*(1+n12)+8*n2*(1+n12)^2-26*n2*(1+n12)+3
      *n2+4*(1+n12)^3-26*(1+n12)^2+48*(1+n12)-15,128*ep^4+32*ep^3*n2+32*ep^3*(1+n12)-192*ep^3-24*ep^2*n2-24*ep^2*(1+n12)+40*ep^2-8*ep*n2-8*
      ep*(1+n12)+48*ep+6*n2+6*(1+n12)-18)
      +`Y'(0,n2,1,1,1,3,-1,1,0,0,0,(1+n12),0,0)*rat(-4*ep*(1+n12)-4*n2*(1+n12)-4*(1+n12)^2+8*(1+n12),16*ep^3-12*ep^2-4*ep+3)
      +`Y'(1,-1+n2,1,1,-1,3,1,1,0,0,0,(1+n12),0,0)*rat(-5*ep*(1+n12)-2*n2*(1+n12)-2*(1+n12)^2+2*(1+n12),2*ep^3+ep^2-ep)
      +`Y'(1,-1+n2,1,1,0,2,1,1,0,0,0,1+(1+n12),-1,0)*rat(2*(1+n12),2*ep-1)
      +`Y'(1,-1+n2,1,1,0,2,1,1,0,0,0,(1+n12),0,0)*rat(-40*ep^3-10*ep^2*n2+6*ep^2*(1+n12)+46*ep^2+12*ep*n2*(1+n12)+5*ep*n2+12*ep*
      (1+n12)^2-16*ep*(1+n12)-13*ep+2*n2^2*(1+n12)+4*n2*(1+n12)^2-8*n2*(1+n12)+2*(1+n12)^3-8*(1+n12)^2+6*(1+n12),16*ep^3+4*ep^2*n2+4*ep^2*(1+n12)-20*ep^2-2
      *ep*n2-2*ep*(1+n12)+6*ep)
      +`Y'(1,-1+n2,1,1,0,3,0,1,0,0,0,(1+n12),0,0)*rat(88*ep^3*(1+n12)+32*ep^2*n2*(1+n12)+32*ep^2*(1+n12)^2+12*ep*n2*(1+n12)+12*ep*(1+n12)^2-36*
      ep*(1+n12)-6*n2*(1+n12)-6*(1+n12)^2+3*(1+n12),32*ep^5+8*ep^4-32*ep^3-2*ep^2+6*ep)
      +`Y'(1,-1+n2,1,1,1,2,0,1,0,0,0,1+(1+n12),-1,0)*rat(-(1+n12),2*ep-1)
      +`Y'(1,-1+n2,1,1,1,2,0,1,0,0,0,(1+n12),0,0)*rat(192*ep^4+48*ep^3*n2-16*ep^3*(1+n12)-256*ep^3-48*ep^2*n2*(1+n12)-36*ep^2*n2-48
      *ep^2*(1+n12)^2+68*ep^2*(1+n12)+36*ep^2-8*ep*n2^2*(1+n12)-16*ep*n2*(1+n12)^2+36*ep*n2*(1+n12)-12*ep*n2-8*ep*(1+n12)^3+36*ep*(1+n12)^2-50*ep*
      (1+n12)+64*ep+9*n2+6*(1+n12)-21,128*ep^4+32*ep^3*n2+32*ep^3*(1+n12)-192*ep^3-24*ep^2*n2-24*ep^2*(1+n12)+40*ep^2-8*ep*n2-8*ep*(1+n12)+
      48*ep+6*n2+6*(1+n12)-18)
      +`Y'(1,-1+n2,1,1,1,3,-1,1,0,0,0,(1+n12),0,0)*rat(-12*ep*(1+n12)-4*n2*(1+n12)-4*(1+n12)^2+4*(1+n12),16*ep^3-12*ep^2-4*ep+3)
      +`Y'(1,-1+n2,1,1,2,1,0,1,0,0,0,(1+n12),0,0)*rat(24*ep^2+6*ep*n2+6*ep*(1+n12)-26*ep-3*n2-2*(1+n12)+7,16*ep^2+4*ep*n2+4*ep*(1+n12)-
      20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,-1+n2,1,1,2,2,-1,1,0,0,0,(1+n12),0,0)*rat(-(1+n12),8*ep^2-10*ep+3)
      +`Y'(1,n2,0,1,1,1,1,1,0,0,0,(1+n12),0,0)*rat(12*ep^2*n2-12*ep^2+7*ep*n2^2+7*ep*n2*(1+n12)-21*ep*n2-7*ep*(1+n12)+14*ep+n2^3+2*
      n2^2*(1+n12)-5*n2^2+n2*(1+n12)^2-6*n2*(1+n12)+8*n2-(1+n12)^2+4*(1+n12)-4,8*ep^2+2*ep*n2+2*ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Y'(1,n2,1,0,0,2,1,1,0,0,0,(1+n12),0,0)*rat(-2*ep*(1+n12)-n2*(1+n12)-(1+n12)^2+(1+n12),2*ep^2-ep)
      +`Y'(1,n2,1,0,1,1,1,1,0,0,0,(1+n12),0,0)*rat(-24*ep^2*n2+24*ep^2-14*ep*n2^2-10*ep*n2*(1+n12)+42*ep*n2+10*ep*(1+n12)-28*ep-2*
      n2^3-3*n2^2*(1+n12)+10*n2^2-n2*(1+n12)^2+9*n2*(1+n12)-16*n2+(1+n12)^2-6*(1+n12)+8,8*ep^2+2*ep*n2+2*ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Y'(1,n2,1,0,1,2,0,1,0,0,0,(1+n12),0,0)*rat(4*ep*(1+n12)+2*n2*(1+n12)+2*(1+n12)^2-3*(1+n12),8*ep^2-10*ep+3)
      +`Y'(1,n2,1,0,2,0,1,1,0,0,0,(1+n12),0,0)*rat(1,1)
      +`Y'(1,n2,1,0,2,1,0,1,0,0,0,(1+n12),0,0)*rat(-40*ep^2-10*ep*n2-10*ep*(1+n12)+46*ep+5*n2+4*(1+n12)-13,16*ep^2+4*ep*n2+4*ep*(1+n12)
      -20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,n2,1,0,2,2,-1,1,0,0,0,(1+n12),0,0)*rat((1+n12),8*ep^2-10*ep+3)
      +`Y'(1,n2,1,1,-1,2,1,1,0,0,0,(1+n12),0,0)*rat(2*ep*(1+n12)+(1+n12)^2,4*ep^2-2*ep)
      +`Y'(1,n2,1,1,-1,2,2,0,0,0,0,(1+n12),0,0)*rat(4*ep*(1+n12)+2*n2*(1+n12)+2*(1+n12)^2-3*(1+n12),8*ep^3+4*ep^2-4*ep)
      +`Y'(1,n2,1,1,-1,3,1,0,0,0,0,(1+n12),0,0)*rat(4*ep*(1+n12)+2*n2*(1+n12)+2*(1+n12)^2-3*(1+n12),4*ep^3+2*ep^2-2*ep)
      +`Y'(1,n2,1,1,0,1,1,1,0,0,0,(1+n12),0,0)*rat(-48*ep^4-60*ep^3*n2-36*ep^3*(1+n12)+92*ep^3-12*ep^2*n2^2-18*ep^2*n2*(1+n12)+20*
      ep^2*n2-6*ep^2*(1+n12)^2-8*ep^2*(1+n12)+28*ep^2-6*ep*n2^2-33*ep*n2*(1+n12)+55*ep*n2-26*ep*(1+n12)^2+87*ep*(1+n12)-85*ep-4*n2^2*(1+n12)+6*
      n2^2-8*n2*(1+n12)^2+27*n2*(1+n12)-25*n2-4*(1+n12)^3+22*(1+n12)^2-39*(1+n12)+27,16*ep^3+4*ep^2*n2+4*ep^2*(1+n12)-4*ep^2+2*ep*n2+2*ep*(1+n12)-
      14*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,n2,1,1,0,1,2,0,0,0,0,1+(1+n12),-1,0)*rat(-(1+n12),2*ep-1)
      +`Y'(1,n2,1,1,0,1,2,0,0,0,0,(1+n12),0,0)*rat(8*ep^3+4*ep^2*(1+n12)+4*ep^2-4*ep-2*n2*(1+n12)-2*(1+n12)^2+3*(1+n12),8*ep^3+4*ep^2-4*ep)
      +`Y'(1,n2,1,1,0,2,0,1,0,0,0,(1+n12),0,0)*rat(-8*ep^4*(1+n12)+16*ep^3*n2*(1+n12)+12*ep^3*(1+n12)^2-56*ep^3*(1+n12)-8*ep^2*n2*(1+n12)-14*
      ep^2*(1+n12)^2+26*ep^2*(1+n12)-10*ep*n2*(1+n12)-12*ep*(1+n12)^2+19*ep*(1+n12)-6*(1+n12),32*ep^5+8*ep^4-32*ep^3-2*ep^2+6*ep)
      +`Y'(1,n2,1,1,0,2,1,0,0,0,0,1+(1+n12),-1,0)*rat(-(1+n12),2*ep-1)
      +`Y'(1,n2,1,1,0,2,1,0,0,0,0,(1+n12),0,0)*rat(8*ep^3+12*ep^2*(1+n12)+4*ep^2+4*ep*n2*(1+n12)+2*ep*(1+n12)^2-4*ep+2*n2*(1+n12)-5*(1+n12),8*
      ep^3+4*ep^2-4*ep)
      +`Y'(1,n2,1,1,0,3,0,0,0,0,0,(1+n12),0,0)*rat(-8*ep*(1+n12)-4*n2*(1+n12)-4*(1+n12)^2+6*(1+n12),16*ep^3-12*ep^2-4*ep+3)
      +`Y'(1,n2,1,1,1,0,1,1,0,-1,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)-2*(1+n12),2*ep-1)
      +`Y'(1,n2,1,1,1,0,1,1,0,0,0,1+(1+n12),-1,0)*rat(4*ep*(1+n12)-4*(1+n12),2*ep-1)
      +`Y'(1,n2,1,1,1,0,1,1,0,0,0,(1+n12),0,0)*rat(-4*ep^2+2*ep*n2-4*ep*(1+n12)+2*ep-n2+3*(1+n12),2*ep-1)
      +`Y'(1,n2,1,1,1,0,2,0,0,0,0,1+(1+n12),-1,0)*rat((1+n12),2*ep-1)
      +`Y'(1,n2,1,1,1,0,2,0,0,0,0,(1+n12),0,0)*rat(-2*ep-(1+n12)+1,2*ep-1)
      +`Y'(1,n2,1,1,1,1,0,1,0,0,0,1+(1+n12),-1,0)*rat(-2*ep*(1+n12)+2*(1+n12),2*ep-1)
      +`Y'(1,n2,1,1,1,1,0,1,0,0,0,(1+n12),0,0)*rat(128*ep^5+160*ep^4*n2+160*ep^4*(1+n12)-384*ep^4+32*ep^3*n2^2+64*ep^3*n2*(1+n12)-
      248*ep^3*n2+32*ep^3*(1+n12)^2-232*ep^3*(1+n12)+312*ep^3-24*ep^2*n2^2-8*ep^2*n2*(1+n12)+56*ep^2*n2+16*ep^2*(1+n12)^2-28*ep^2*(1+n12)+8
      *ep*n2^2*(1+n12)-8*ep*n2^2+16*ep*n2*(1+n12)^2-74*ep*n2*(1+n12)+62*ep*n2+8*ep*(1+n12)^3-66*ep*(1+n12)^2+156*ep*(1+n12)-86*ep-4*n2^2*(1+n12)+6*
      n2^2-8*n2*(1+n12)^2+33*n2*(1+n12)-24*n2-4*(1+n12)^3+27*(1+n12)^2-54*(1+n12)+24,64*ep^4+16*ep^3*n2+16*ep^3*(1+n12)-96*ep^3-12*ep^2*n2-12*
      ep^2*(1+n12)+20*ep^2-4*ep*n2-4*ep*(1+n12)+24*ep+3*n2+3*(1+n12)-9)
      +`Y'(1,n2,1,1,1,1,1,0,0,-1,0,1+(1+n12),0,0)*rat(-(1+n12),1)
      +`Y'(1,n2,1,1,1,1,1,0,0,0,0,(1+n12),0,0)*rat(32*ep^3+56*ep^2*n2+16*ep^2*(1+n12)-80*ep^2+28*ep*n2^2+22*ep*n2*(1+n12)-86*ep*n2+
      2*ep*(1+n12)^2-28*ep*(1+n12)+62*ep+4*n2^3+6*n2^2*(1+n12)-20*n2^2+2*n2*(1+n12)^2-17*n2*(1+n12)+31*n2-2*(1+n12)^2+11*(1+n12)-15,16*ep^2+4*ep*n2
      +4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,n2,1,1,1,2,-1,1,0,0,0,(1+n12),0,0)*rat(-4*ep*n2*(1+n12)-4*ep*(1+n12)^2+8*ep*(1+n12)+2*n2*(1+n12)+2*(1+n12)^2-6*(1+n12),16*ep^3-12*ep^2
      -4*ep+3)
      +`Y'(1,n2,1,1,1,2,0,0,0,0,0,(1+n12),0,0)*rat(-8*ep^2*(1+n12)-4*ep*n2*(1+n12)-2*ep*(1+n12)^2+6*ep*(1+n12)+(1+n12)^2,16*ep^3-12*ep^2-4*ep+3
      )
      +`Y'(0,n2,1,1,0,2,1,1,0,0,0,1+(1+n12),0,0)*rat((1+n12)^2,2*ep^2-ep)
      +`Y'(0,n2,1,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(-n2*(1+n12)+(1+n12),8*ep^2+2*ep*n2+2*ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Y'(0,n2,1,1,1,2,0,1,0,0,0,1+(1+n12),0,0)*rat(-2*(1+n12)^2,8*ep^2-10*ep+3)
      +`Y'(1,n2,1,1,1,0,1,1,0,0,0,1+(1+n12),0,0)*rat(-16*ep^2*(1+n12)-8*ep*n2*(1+n12)-4*ep*(1+n12)^2+32*ep*(1+n12)-n2^2*(1+n12)-n2*(1+n12)^2+7*n2*
      (1+n12)+5*(1+n12)^2-14*(1+n12),16*ep^2+4*ep*n2+4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,n2,1,1,1,1,1,0,0,0,0,1+(1+n12),0,0)*rat(16*ep^2*(1+n12)+4*ep*(1+n12)^2-16*ep*(1+n12)-n2^2*(1+n12)-n2*(1+n12)^2+n2*(1+n12)-(1+n12)^2+4*(1+n12),
      16*ep^2+4*ep*n2+4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(0,-1+n2,1,1,0,2,1,1,0,0,0,1+(1+n12),0,0)*rat(-(1+n12)^2,2*ep^2-ep)
      +`Y'(0,-1+n2,1,1,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(-8*ep^2*(1+n12)-10*ep*n2*(1+n12)-2*ep*(1+n12)^2+18*ep*(1+n12)-2*n2^2*(1+n12)-2*n2*(1+n12)^2
      +8*n2*(1+n12)+3*(1+n12)^2-8*(1+n12),8*ep^2+2*ep*n2+2*ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Y'(0,-1+n2,1,1,1,2,0,1,0,0,0,1+(1+n12),0,0)*rat(2*(1+n12)^2,8*ep^2-10*ep+3)
      +`Y'(0,n2,1,0,0,2,1,1,0,0,0,1+(1+n12),0,0)*rat((1+n12)^2,2*ep^2-ep)
      +`Y'(0,n2,1,0,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(8*ep*n2*(1+n12)-8*ep*(1+n12)+2*n2^2*(1+n12)+2*n2*(1+n12)^2-7*n2*(1+n12)-2*(1+n12)^2+5*(1+n12),8*
      ep^2+2*ep*n2+2*ep*(1+n12)-10*ep-n2-(1+n12)+3)
      +`Y'(0,n2,1,0,1,2,0,1,0,0,0,1+(1+n12),0,0)*rat(-2*(1+n12)^2,8*ep^2-10*ep+3)
      +`Y'(0,n2,1,1,0,1,1,1,0,0,0,1+(1+n12),0,0)*rat(2*ep*(1+n12),2*ep-1)
      +`Y'(0,n2,1,1,1,0,1,1,0,0,0,1+(1+n12),0,0)*rat(-4*ep*n2*(1+n12)+4*ep*(1+n12)-n2^2*(1+n12)-n2*(1+n12)^2+3*n2*(1+n12)+(1+n12)^2-2*(1+n12),16*ep^2+4
      *ep*n2+4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(0,n2,1,1,1,1,1,0,0,-1,0,2+(1+n12),0,0)*rat(-2*ep*(1+n12)^2-2*ep*(1+n12)-n2*(1+n12)^2-n2*(1+n12)+(1+n12)^2+(1+n12),4*ep^2-6*ep+2)
      +`Y'(0,n2,1,1,1,1,1,0,0,0,0,1+(1+n12),0,0)*rat(-4*ep*n2*(1+n12)+4*ep*(1+n12)-n2^2*(1+n12)-n2*(1+n12)^2+3*n2*(1+n12)+(1+n12)^2-2*(1+n12),16*ep^2+4
      *ep*n2+4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,-1+n2,1,0,1,1,1,1,0,0,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12),2*ep-1)
      +`Y'(1,-1+n2,1,1,0,1,1,1,0,0,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)+2*(1+n12),2*ep-1)
      +`Y'(1,-1+n2,1,1,1,0,1,1,0,0,0,1+(1+n12),0,0)*rat(-3*(1+n12),2*ep-1)
      +`Y'(1,n2,0,1,1,0,1,1,0,0,0,1+(1+n12),0,0)*rat(4*ep*n2*(1+n12)-4*ep*(1+n12)+n2^2*(1+n12)+n2*(1+n12)^2-3*n2*(1+n12)-(1+n12)^2+2*(1+n12),16*ep^2+4*
      ep*n2+4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,n2,0,1,1,1,1,0,0,0,0,1+(1+n12),0,0)*rat(4*ep*n2*(1+n12)-4*ep*(1+n12)+n2^2*(1+n12)+n2*(1+n12)^2-3*n2*(1+n12)-(1+n12)^2+2*(1+n12),16*ep^2+4*
      ep*n2+4*ep*(1+n12)-20*ep-2*n2-2*(1+n12)+6)
      +`Y'(1,n2,1,0,1,0,1,1,0,0,0,1+(1+n12),0,0)*rat(-2*ep*(1+n12)+2*(1+n12),2*ep-1)
      +`Y'(1,n2,1,0,1,1,1,0,0,0,0,1+(1+n12),0,0)*rat((1+n12),1)
      +`Y'(1,n2,1,1,0,0,1,1,0,0,0,1+(1+n12),0,0)*rat(2*ep*(1+n12)-2*(1+n12),2*ep-1)
      +`Y'(1,n2,1,1,0,1,1,0,0,0,0,1+(1+n12),0,0)*rat(-(1+n12),1)
      );
*--#] n12 : 
#endif
#if ( `il' >= 16 )
*--#[ n2 :
`SKIP2'
id,ifmatch->bebe1,`Z'(1,n2?{>1},1,1,1,1,1,1,0,0,0,0,0,0)
			 = -rat(12*ep+3*(n2-1)-6,6*ep*(n2-1)+2*(n2-1)^2-2*(n2-1))*(
      +`Z'(1,(n2-1),1,1,1,1,1,1,0,0,0,0,0,0)*rat(-10*ep-2*(n2-1)+4,3)
      +`Y'(0,1+(n2-1),1,1,0,2,1,1,0,0,0,0,0,0)*rat(24*ep^2+6*ep*(n2-1)-28*ep-3*(n2-1)+8,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(0,1+(n2-1),1,1,1,1,1,1,0,0,0,0,0,0)*rat(-2*ep*(n2-1)+2*(n2-1)^2,12*ep+3*(n2-1)-6)
      +`Y'(0,1+(n2-1),1,1,1,1,1,2,0,0,0,-1,0,0)*rat(-(n2-1),4*ep+(n2-1)-2)
      +`Y'(0,1+(n2-1),1,1,1,2,0,1,0,0,0,0,0,0)*rat(-8*ep^2-2*ep*(n2-1)+12*ep+(n2-1)-4,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(0,(n2-1),1,1,1,1,1,2,0,0,0,0,0,0)*rat(4*ep+2*(n2-1)-1,4*ep+(n2-1)-2)
      +`Y'(1,-1+(n2-1),0,2,1,1,1,2,0,0,0,0,0,0)*rat(1,12*ep+3*(n2-1)-6)
      +`Y'(1,-1+(n2-1),1,1,0,2,2,1,0,0,0,0,0,0)*rat(1,12*ep+3*(n2-1)-6)
      +`Y'(1,-1+(n2-1),1,1,2,0,1,2,0,0,0,0,0,0)*rat(-2,12*ep+3*(n2-1)-6)
      +`Y'(1,1+(n2-1),0,1,1,1,1,1,0,0,0,0,0,0)*rat(-6*ep*(n2-1)-2*(n2-1)^2+2*(n2-1),12*ep+3*(n2-1)-6)
      +`Y'(1,1+(n2-1),1,0,1,1,1,1,0,0,0,0,0,0)*rat(12*ep*(n2-1)+4*(n2-1)^2-4*(n2-1),12*ep+3*(n2-1)-6)
      +`Y'(1,1+(n2-1),1,0,2,0,1,1,0,0,0,0,0,0)*rat(-4*ep+2,12*ep+3*(n2-1)-3)
      +`Y'(1,1+(n2-1),1,1,0,1,1,1,0,0,0,0,0,0)*rat(48*ep^3+60*ep^2*(n2-1)-80*ep^2+12*ep*(n2-1)^2-56*ep*(n2-1)+44*ep-6*(n2-1)^2+13*(n2-1)-8,48*
      ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,1+(n2-1),1,1,0,1,2,0,0,0,0,0,0,0)*rat(-4*ep+2,12*ep+3*(n2-1)-3)
      +`Y'(1,1+(n2-1),1,1,0,1,2,1,0,0,0,-1,0,0)*rat(-(n2-1),4*ep+(n2-1)-2)
      +`Y'(1,1+(n2-1),1,1,0,2,1,0,0,0,0,0,0,0)*rat(-4*ep+2,12*ep+3*(n2-1)-3)
      +`Y'(1,1+(n2-1),1,1,0,2,1,1,0,0,0,-1,0,0)*rat(-40*ep^2-10*ep*(n2-1)+36*ep+5*(n2-1)-8,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,1+(n2-1),1,1,1,0,1,1,0,0,0,0,0,0)*rat(8*ep^2-4*ep*(n2-1)-8*ep+2*(n2-1)+2,12*ep+3*(n2-1)-3)
      +`Y'(1,1+(n2-1),1,1,1,0,1,2,0,0,0,-1,0,0)*rat((n2-1),4*ep+(n2-1)-2)
      +`Y'(1,1+(n2-1),1,1,1,0,2,1,0,0,0,-1,0,0)*rat((n2-1),4*ep+(n2-1)-2)
      +`Y'(1,1+(n2-1),1,1,1,1,0,1,0,0,0,0,0,0)*rat(-32*ep^3-40*ep^2*(n2-1)+48*ep^2-8*ep*(n2-1)^2+36*ep*(n2-1)-24*ep+4*(n2-1)^2-8*(n2-1)+4,48*ep^2
      +24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,1+(n2-1),1,1,1,1,0,2,0,0,0,-1,0,0)*rat(24*ep^2+6*ep*(n2-1)-20*ep-3*(n2-1)+4,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,1+(n2-1),1,1,1,1,1,0,0,0,0,0,0,0)*rat(-32*ep^3-56*ep^2*(n2-1)+24*ep^2-28*ep*(n2-1)^2+30*ep*(n2-1)-4*ep-4*(n2-1)^3+8*(n2-1)^2-3*(n2-1),48
      *ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,1+(n2-1),1,1,2,0,1,1,-1,0,0,0,0,0)*rat((n2-1),12*ep+3*(n2-1)-6)
      +`Y'(1,1+(n2-1),1,1,2,0,1,1,0,-1,0,0,0,0)*rat(16*ep^2+8*ep*(n2-1)-16*ep+(n2-1)^2-3*(n2-1)+4,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,1+(n2-1),1,1,2,0,1,1,0,0,0,0,0,0)*rat(-16*ep^2-12*ep*(n2-1)+16*ep-2*(n2-1)^2+4*(n2-1)-4,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6
      )
      +`Y'(1,1+(n2-1),1,1,2,1,0,1,0,-1,0,0,0,0)*rat(-48*ep^3-12*ep^2*(n2-1)+80*ep^2+16*ep*(n2-1)-44*ep-5*(n2-1)+8,48*ep^3+24*ep^2*(n2-1)-84*
      ep^2+3*ep*(n2-1)^2-33*ep*(n2-1)+42*ep-3*(n2-1)^2+9*(n2-1)-6)
      +`Y'(1,1+(n2-1),1,1,2,1,0,1,0,0,0,0,0,0)*rat(40*ep^2+10*ep*(n2-1)-36*ep-5*(n2-1)+8,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,2+(n2-1),1,0,1,1,1,1,0,0,0,-1,0,0)*rat(-2*(n2-1)^2-2*(n2-1),12*ep+3*(n2-1)-6)
      +`Y'(1,2+(n2-1),1,1,0,1,1,1,0,0,0,-1,0,0)*rat(-8*ep^2*(n2-1)-8*ep^2-2*ep*(n2-1)^2+2*ep*(n2-1)+4*ep+(n2-1)^2+(n2-1),48*ep^2+24*ep*(n2-1)-36*ep+3
      *(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,2+(n2-1),1,1,0,1,2,0,0,-1,0,0,0,0)*rat(-2*ep*(n2-1)-2*ep-(n2-1)^2-(n2-1),12*ep^2+3*ep*(n2-1)-15*ep-3*(n2-1)+3)
      +`Y'(1,2+(n2-1),1,1,1,0,1,1,-1,0,0,0,0,0)*rat((n2-1)^2+(n2-1),12*ep+3*(n2-1)-6)
      +`Y'(1,2+(n2-1),1,1,1,1,1,0,-1,0,0,0,0,0)*rat((n2-1)^2+(n2-1),12*ep+3*(n2-1)-6)
      +`Y'(1,2+(n2-1),1,1,1,1,1,0,0,-1,0,0,0,0)*rat(6*ep^2*(n2-1)+6*ep^2+5*ep*(n2-1)^2+3*ep*(n2-1)-2*ep+(n2-1)^3-(n2-1),12*ep^2+3*ep*(n2-1)-15*ep-3*
      (n2-1)+3)
      +`Y'(1,(n2-1),0,2,1,1,1,2,0,0,0,-1,0,0)*rat(-1,12*ep+3*(n2-1)-6)
      +`Y'(1,(n2-1),1,1,0,1,2,1,0,0,0,0,0,0)*rat(4*ep+2*(n2-1)-3,4*ep+(n2-1)-2)
      +`Y'(1,(n2-1),1,1,0,2,1,1,0,0,0,0,0,0)*rat(40*ep^2+10*ep*(n2-1)-36*ep-5*(n2-1)+8,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,(n2-1),1,1,0,2,2,1,0,0,0,-1,0,0)*rat(-1,12*ep+3*(n2-1)-6)
      +`Y'(1,(n2-1),1,1,1,0,1,2,0,0,0,0,0,0)*rat(-4*ep-2*(n2-1)+3,4*ep+(n2-1)-2)
      +`Y'(1,(n2-1),1,1,1,0,2,1,0,0,0,0,0,0)*rat(-4*ep-2*(n2-1)+3,4*ep+(n2-1)-2)
      +`Y'(1,(n2-1),1,1,1,2,0,1,0,0,0,0,0,0)*rat(-24*ep^2-6*ep*(n2-1)+20*ep+3*(n2-1)-4,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      +`Y'(1,(n2-1),1,1,2,0,1,2,0,0,0,-1,0,0)*rat(2,12*ep+3*(n2-1)-6)
      +`Y'(1,(n2-1),1,1,2,1,0,1,0,0,0,0,0,0)*rat(-24*ep^2-6*ep*(n2-1)+20*ep+3*(n2-1)-4,48*ep^2+24*ep*(n2-1)-36*ep+3*(n2-1)^2-9*(n2-1)+6)
      );
*--#] n2 : 
#endif
goto bebe2;
Label bebe1;
$ired = 1;
Label bebe2;
id,`Z'(?a,x?neg0_,?b,n9?,...,n14?) = `Y'(?a,x,?b,n9,...,n14);
*id	`Y'(?a) = 0;
if ( count(`Y',1) );
id	`Y'(n1?neg0_,n2?neg0_,?a) = 0;
id	`Y'(n1?,n2?,n3?neg0_,n4?neg0_,?a) = 0;
id	`Y'(n1?neg0_,n2?,n3?neg0_,?a) = 0;
id	`Y'(n1?neg0_,n2?,n3?,n4?,n5?,n6?,n7?,n8?neg0_,?a) = 0;
id	`Y'(n1?,n2?,n3?neg0_,n4?,n5?,n6?,n7?,n8?neg0_,?a) = 0;
id	`Y'(n1?,n2?,n3?,n4?,n5?,n6?neg0_,n7?,n8?neg0_,?a) = 0;
id	`Y'(n1?,n2?,n3?,n4?,n5?,n6?,n7?neg0_,n8?neg0_,?a) = 0;
id	`Y'(n1?,n2?,n3?,n4?,n5?,n6?neg0_,n7?neg0_,n8?,?a) = 0;
id	`Y'(n1?,n2?neg0_,n3?,n4?,n5?neg0_,n6?,n7?,n8?,?a) = 0;
id	`Y'(n1?,n2?neg0_,n3?,n4?,n5?,n6?neg0_,n7?,n8?,?a) = 0;
id	`Y'(n1?,n2?,n3?,n4?neg0_,n5?neg0_,n6?,n7?,n8?,?a) = 0;
id	`Y'(n1?,n2?,n3?,n4?neg0_,n5?,n6?,n7?neg0_,n8?,?a) = 0;
id	`Y'(n1?,n2?,n3?,n4?,n5?neg0_,n6?neg0_,n7?,n8?,?a) = 0;
id	`Y'(n1?,n2?,n3?,n4?,n5?neg0_,n6?,n7?neg0_,n8?,?a) = 0;
endif;
ModuleOption,maximum,$ired;
.sort:Main bebe loop, pass `$pass';
#if ( `$ired' > 0 )
	#redefine ired "0"
#endif
#ifndef `NOSPECTATORS'
if ( count(`Z',1) == 0 ) ToSpectator x`TOPO';
#endif
#enddo
#enddo
id,`Z'(1,1,1,1,1,1,1,1,0,0,0,0,0,0) = Master(bebe);
*
*--#] Loop : 
#break
*--#] bebe : 
#endswitch
.sort:FORCERD-`TOPO';
#ifndef `NOSPECTATORS'
#call FromSpectator(x`TOPO',`InName',`OutName',`Y')
#endif
#endprocedure
*
*--#] FORCERD : 
*--#[ FromSpectator :
*
#procedure FromSpectator(SPECTATOR,INNAME,OUTNAME,Y)
*
#ifndef `NOSPECTATORS'
PolyRatFun;
CopySpectator `OUTNAME' = `SPECTATOR';
AB+	rat;
.sort:FromSpectator-`SPECTATOR';

#ifndef `POLYRATFUNEXPAND'
  #define POLYRATFUNEXPAND "-1"
#endif

#if `POLYRATFUNEXPAND' < 0
  PolyRatFun rat,RAT;
#elseif `POLYRATFUNEXPAND' == 0
  PolyRatFun rat,RAT(divergence,ep);
#else
  PolyRatFun rat,RAT(expand,ep,`POLYRATFUNEXPAND');
#endif

RemoveSpectator `SPECTATOR';
#else
Local `OUTNAME' = `INNAME';
if ( expression(`INNAME') && count(`Y',1) ) Discard;
if ( expression(`OUTNAME') && ( count(`Y',1) == 0 ) ) Discard;
#endif
*
#endprocedure
*
*--#] FromSpectator : 
#endif
