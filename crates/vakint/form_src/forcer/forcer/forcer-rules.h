#ifndef `FORCERRULESHFILE'
#define FORCERRULESHFILE

#include- forcer-manual.h

* Override the default reductions.
#define HASSCHEMEZno3

*--#[ decl :

S n1,...,n19;
S x1,...,x19;
*
* BasicT1Integral and BasicNOIntegral are defined as in MincerExact:
*   Zt1star5(1,1,1,1,1) = BasicT1Integral * ep^2 / Q^2
*   Zno(1,1,1,1,1,1,1,1,0) = BasicNOIntegral / (Q^2)^2
* On the other hand, others are defined in a natural way, e.g.,
*   Zt1star55(1,1,1,1,1) = BasicT1Integral55 / Q^2
*
* TODO: A more unified convention for these integrals?
*
S BasicT1Integral,BasicNOIntegral;
S BasicT1Integral24,BasicT1Integral45,BasicT1Integral55;
S BasicNOIntegral5,BasicNOIntegral6;
CF Z,Y;
S Q2;

*--#] decl : 
*--#[ proc :
*--#[ Reduce :

**
* The main entry point of the reduction routine.
*
* We assume that `Q2' is correctly defined from the value of `SETQONE'.
*
#procedure Reduce(?TOPO)
  #call Timing()

  #define TOPO "`?TOPO'"
  #switch `TOPO'
* for these topologies the FORCERD works.
    #case no6
    #case bebe
    #case cross
    #case haha
    #case no1
    #case bustar5
    #case fastar3
    #case lala
    #case nono
      #call ReduceByFORCERD(`TOPO')
      #break
*--#[ no2 :
    #case no2
      #call ReduceImpl2(`TOPO',Z`TOPO',14)
      id Z`TOPO'(1,1,1,1,1,1,1,1,1,1,1,0,0,0) = Master(no2) / `Q2'^3;
      .sort:`TOPO'-end;
      #break
*--#] no2 : 
*--#[ bubu :
    #case bubu
      #call ReduceImpl2(`TOPO',Z`TOPO',14)
      .sort:`TOPO'-end;
      #break
*--#] bubu : 
*--#[ no4 :
    #case no4
*[(0, 1), (1, 8), (1, 6), (1, 5), (4, 8), (4, 6), (4, 5), (7, 8), (6, 7), (5, 7), (4, 9)]
*[(0, 1), (0, 2), (0, 4), (0, 5), (2, 4)]
repeat;
id  Zno4(n1?{>0},n2?{>0},n3?{>0},n4?{>0},n5?{>0},n6?{>0},n7?{>0},n8?{>0},n9?{>0},
    n10?,n11?,n12?,n13?,n14?) = 
        -rat(1,-4*ep-n1-n2-n3-2*n7-2*n8-2*n9-n10-n11-n14+8) * (
      +Zno4(1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(n1,1)
      +Zno4(1+n1,n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n1,1)
      +Zno4(n1,1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(n2,1)
      +Zno4(n1,1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(-n2,1)
      +Zno4(n1,n2,1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(n3,1)
      +Zno4(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,n14)*rat(-n3,1))
   ;
endrepeat;
.sort
#break
*--#] no4 : 
*--#[ no5 :
    #case no5
repeat;
id      Zno5(n1?{>0},n2?{>0},n3?{>0},n4?{>0},n5?{>0},n6?{>0},n7?{>0},n8?{>0},n9?{>0},
    n10?{>0},n11?,n12?,n13?,n14?) = 
        -rat(1,-4*ep-n2-n3-n4-2*n8-2*n9-2*n10-n11-n12-n14+8) * (
      +Zno5(n1,1+n2,n3,n4,-1+n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(n2,1)
      +Zno5(n1,1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(-n2,1)
      +Zno5(n1,n2,1+n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(n3,1)
      +Zno5(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,n14)*rat(-n3,1)
      +Zno5(n1,n2,n3,1+n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n4,1)
      +Zno5(n1,n2,n3,1+n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,n14)*rat(-n4,1))
      ;
endrepeat;
.sort
#break
*--#] no5 : 
*--#[ no7 :
      #case no7
repeat;
id  Zno7(n1?{>0},n2?{>0},n3?{>0},n4?{>0},n5?{>0},n6?{>0},n7?{>0},n8?{>0},n9?{>0},
    n10?{>0},n11?,n12?,n13?,n14?) = 
        -rat(1,-4*ep-n1-n2-n3-2*n6-2*n8-2*n10-n11-n12-n14+8) * (
      +Zno7(1+n1,n2,n3,-1+n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(n1,1)
      +Zno7(1+n1,n2,n3,n4,n5,-1+n6,n7,n8,n9,n10,n11,n12,n13,n14)*rat(-n1,1)
      +Zno7(n1,1+n2,n3,n4,n5,n6,-1+n7,n8,n9,n10,n11,n12,n13,n14)*rat(n2,1)
      +Zno7(n1,1+n2,n3,n4,n5,n6,n7,-1+n8,n9,n10,n11,n12,n13,n14)*rat(-n2,1)
      +Zno7(n1,n2,1+n3,n4,n5,n6,n7,n8,-1+n9,n10,n11,n12,n13,n14)*rat(n3,1)
      +Zno7(n1,n2,1+n3,n4,n5,n6,n7,n8,n9,-1+n10,n11,n12,n13,n14)*rat(-n3,1))
      ;
endrepeat;
.sort
#break
*--#] no7 : 
*--#[ t1star5 :
    #case t1star5
      id Z`TOPO'(n1?,...,n5?) = Z(0,0,0,0,1,n1,...,n5) * `Q2'^(-n1-...-n5);
      #call ReduceImpl(t1,5)
      id Z(0,0,0,0,1,1,1,1,1,1) = BasicT1Integral * rat(ep^2,1) / `Q2' * `Q2'^5;
      id Z(0,0,0,0,1,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      id Y(0,0,0,0,1,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      .sort:`TOPO'-end;
      #break
*--#] t1star5 : 
*--#[ no/nostar0 :
    #case no
    #case nostar0
*     NOSCHEME: 0, 1, 2 (default: 1)
      #ifndef `NOSCHEME'
        #define NOSCHEME "1"
      #endif
#if `NOSCHEME' == 0
      id Z`TOPO'(n1?,...,n9?) = Z(0,0,0,0,0,0,0,0,0,n1,...,n9) * `Q2'^(-n1-...-n9);
      #call ReduceImpl(no,none)
      id Z(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0) = BasicNOIntegral / `Q2'^2 * `Q2'^8;
      id Z(0,0,0,0,0,0,0,0,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
      id Y(0,0,0,0,0,0,0,0,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
#else
      #call ReduceImpl2(`TOPO',Z`TOPO',9)
      id Z`TOPO'(1,1,1,1,1,1,1,1,0) = BasicNOIntegral / `Q2'^2;
#endif
      .sort:`TOPO'-end;
      #break
*--#] no/nostar0 : 
*--#[ lastar5 :
    #case lastar5
*     LASTAR5SCHEME: 0, 1, 2, 3 (default: 3)
      #ifndef `LASTAR5SCHEME'
        #define LASTAR5SCHEME "3"
      #endif
#if (`LASTAR5SCHEME' == 0) || (`LASTAR5SCHEME' == 1)
      #if `LASTAR5SCHEME' == 0
        #define ALTLA135 "1"
      #else
        #define ALTLA135 ""
      #endif
      id Z`TOPO'(n1?,...,n9?) = Z(0,0,0,0,1,0,0,0,0,n1,...,n9) * `Q2'^(-n1-...-n9);
      #call ReduceImpl(la13,5)
      id Z(0,0,0,0,1,0,0,0,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
      id Y(0,0,0,0,1,0,0,0,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
#else
      #call ReduceImpl2(`TOPO',Z`TOPO',9)
#endif
      .sort:`TOPO'-end;
      #break
*--#] lastar5 : 
*--#[ fastar2 :
    #case fastar2
*     FASTAR2SCHEME: 0, 1 (default: 1)
      #ifndef `FASTAR2SCHEME'
        #define FASTAR2SCHEME "1"
      #endif
#if `FASTAR2SCHEME' == 0
      #call ReduceByFORCERD(`TOPO')
#else
      #call ReduceImpl2(`TOPO',Z`TOPO',9)
      id Z`TOPO'(1,1,1,1,1,1,1,0,0) = Master(fastar2) / `Q2';
#endif
      .sort:`TOPO'-end;
      #break
*--#] fastar2 : 
*--#[ nostar5 :
    #case nostar5
*     NOSTAR5SCHEME: 0, 1 (default: 1)
      #ifndef `NOSTAR5SCHEME'
        #define NOSTAR5SCHEME "1"
      #endif
#if `NOSTAR5SCHEME' == 0
      id Z`TOPO'(n1?,...,n9?) = Z(0,0,0,0,1,0,0,0,0,n1,...,n9) * `Q2'^(-n1-...-n9);
      #call ReduceImpl(no,5)
      id Z(0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,0) = BasicNOIntegral5 / `Q2'^2 * `Q2'^8;  * XXX: normalization?
      id Z(0,0,0,0,1,0,0,0,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
      id Y(0,0,0,0,1,0,0,0,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
#else
      #call ReduceImpl2(`TOPO',Z`TOPO',9)
      id Z`TOPO'(1,1,1,1,1,1,1,1,0) = BasicNOIntegral5 / `Q2'^2;
*     id Z`TOPO'(1,1,1,1,1,1,1,1,0) = Master(nostar5) / `Q2'^2;
#endif
      .sort:`TOPO'-end;
      #break
*--#] nostar5 : 
*--#[ nostar6 :
    #case nostar6
*     NOSTAR6SCHEME: 0, 1 (default: 1)
      #ifndef `NOSTAR6SCHEME'
        #define NOSTAR6SCHEME "1"
      #endif
#if `NOSTAR6SCHEME' == 0
      id Z`TOPO'(n1?,...,n9?) = Z(0,0,0,0,0,1,0,0,0,n1,...,n9) * `Q2'^(-n1-...-n9);
      #call ReduceImpl(no,6)
      id Z(0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,0) = BasicNOIntegral6 / `Q2'^2 * `Q2'^8;  * XXX: normalization?
      id Z(0,0,0,0,0,1,0,0,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
      id Y(0,0,0,0,0,1,0,0,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
#else
      #call ReduceImpl2(`TOPO',Z`TOPO',9)
      id Z`TOPO'(1,1,1,1,1,1,1,1,0) = BasicNOIntegral6 / `Q2'^2;
*     id Z`TOPO'(1,1,1,1,1,1,1,1,0) = Master(nostar6) / `Q2'^2;
#endif
      .sort:`TOPO'-end;
      #break
*--#] nostar6 : 
*--#[ nostar8 :
    #case nostar8
      id Z`TOPO'(n1?,...,n9?) = Z(0,0,0,0,0,0,0,1,0,n1,...,n9) * `Q2'^(-n1-...-n9);
      #call ReduceImpl(no,8)
      id Z(0,0,0,0,0,0,0,1,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
      id Y(0,0,0,0,0,0,0,1,0,n1?,...,n9?) = Z`TOPO'(n1,...,n9) * `Q2'^(n1+...+n9);
*     Mapped to nostar5 when m9=0.
      id Znostar8(n1?,n2?,n3?,n4?,n5?,n6?,n7?,n8?,0) = Znostar5(n1,n7,n4,n3,n8,n6,n2,n5,0);
      .sort:`TOPO'-end;
      #break
*--#] nostar8 : 
*--#[ t1star05 :
    #case t1star05
      id Z`TOPO'(n1?,...,n5?) = Z(0,0,0,0,1,n1,...,n5) * `Q2'^(-n1-...-n5);
      #call ReduceImpl(t1,5)
      id Z(0,0,0,0,1,1,1,1,1,1) = BasicT1Integral * rat(ep^2,1) / `Q2' * `Q2'^5;
      id Z(0,0,0,0,1,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      id Y(0,0,0,0,1,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      .sort:`TOPO'-end;
      #break
*--#] t1star05 : 
*--#[ t1star24 :
    #case t1star24
      id Z`TOPO'(n1?,...,n5?) = Z(0,1,0,1,0,n1,...,n5) * `Q2'^(-n1-...-n5);
      #call ReduceImpl(t1,24)
      id Z(0,1,0,1,0,1,1,1,1,1) = BasicT1Integral24 / `Q2' * `Q2'^5;  * XXX: normalization?
      id Z(0,1,0,1,0,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      id Y(0,1,0,1,0,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      .sort:`TOPO'-end;
      #break
*--#] t1star24 : 
*--#[ t1star34 :
    #case t1star34
      id Z`TOPO'(n1?,...,n5?) = Z(0,0,1,1,0,n1,...,n5) * `Q2'^(-n1-...-n5);
      #call ReduceImpl(t1,34)
      id Z(0,0,1,1,0,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      id Y(0,0,1,1,0,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      .sort:`TOPO'-end;
      #break
*--#] t1star34 : 
*--#[ t1star45 :
    #case t1star45
      id Z`TOPO'(n1?,...,n5?) = Z(0,0,0,1,1,n1,...,n5) * `Q2'^(-n1-...-n5);
      #call ReduceImpl(t1,45)
      id Z(0,0,0,1,1,1,1,1,1,1) = BasicT1Integral45 / `Q2' * `Q2'^5;  * XXX: normalization?
      id Z(0,0,0,1,1,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      id Y(0,0,0,1,1,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      .sort:`TOPO'-end;
      #break
*--#] t1star45 : 
*--#[ t1star55 :
    #case t1star55
      id Z`TOPO'(n1?,...,n5?) = Z(0,0,0,0,2,n1,...,n5) * `Q2'^(-n1-...-n5);
      #call ReduceImpl(t1,5)
      id Z(0,0,0,0,2,1,1,1,1,1) = BasicT1Integral55 / `Q2' * `Q2'^5;  * XXX: normalization?
      id Z(0,0,0,0,2,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      id Y(0,0,0,0,2,n1?,...,n5?) = Z`TOPO'(n1,...,n5) * `Q2'^(n1+...+n5);
      .sort:`TOPO'-end;
      #break
*--#] t1star55 : 
*--#[ default :
    #default
*     Assume that we have SchemeZ`TOPO'().
      #call ReduceImpl2(`TOPO',Z`TOPO',)
      .sort:`TOPO'-end;
      #break
*--#] default : 
  #endswitch

  #call Timing(Reduce(`TOPO'))
#endprocedure

*--#] Reduce : 
*--#[ ReduceByFORCERD :

**
* Applies FORCERD rules.
*
#procedure ReduceByFORCERD(TOPO)
* NOTE: Q^2 = 1 assumed.
  #if `SETQONE' == 0
    #message Warning: FORCERD does not support SETQONE=0
  #endif

  #call FORCERD(`TOPO',Z`TOPO',Y`TOPO',F`TOPO'res,F`TOPO'res)

  #ifndef `NOSPECTATORS'
*   if we know the name of the input, merge all expressions
    .sort:Merge spectator `TOPO';
    #ifdef `NAME'
      Drop `NAME'`TOPO',F`TOPO'res;
      L `NAME'tmp = `NAME'`TOPO' + F`TOPO'res;
      #exchange `NAME'`TOPO',`NAME'tmp
      ModuleOption noparallel;
      .sort:spectator merge `TOPO';
    #endif
  #endif

  Multiply replace_(Y`TOPO', Z`TOPO');
#endprocedure

*--#] ReduceByFORCERD : 
*--#[ ReduceImpl :

**
* Reduction schemes for topologies with insertions. Non-integer parts of indices
* are expressed like Z(x1,x2,x3,x4,n1,n2,n3,n4).
*
#procedure ReduceImpl(?TOPO)
  #define TOPO "`?TOPO'"
  #define IDSPEC "ifmatch->1"
  #define Z "Z"
  #define Y "Y"
  #switch `TOPO'
*--#[ la13,5
    #case la13,5
      #define nloop "1"
      #do loop=1,1
        #$doloop = 0;
        if (count(`Z',1) == 0) goto 2;
        #ifndef `ALTLA135'
*--#[ la13,5-rules-a
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* DEN = den(-1+n4)
* MAX = Z(0,0,0,0,0,0,0,1,0)
* MIN = Z(0,-1,-1,-1,0,0,0,0,0)
* nterms = 4 (Z: 4, Y: 0)
* NonIntegerIndices: 5
* Signs: o,+,+,+,o,o,o,o,o
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?,n2?pos_,n3?pos_,n4?{>=2},n5?,n6?,n7?,n8?,n9?) =
*      #include la.rules # la,noints=5,signs=o,+,+,+,o,o,o,o,o,level=1,stop=0
*
*--#[ la,noints=5,signs=o,+,+,+,o,o,o,o,o,level=1,stop=0 :
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,-1+n4,n5,n6,n7,1+n8,n9)*(rat(-n8,-1+
         n4))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,-1+n4,n5,n6,n7,1+n8,n9)*(rat(n8,-1+n4
         ))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,n4,n5,n6,n7,n8,n9)*(rat(1,1))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,-1+n4,n5,n6,n7,n8,n9)*(rat(-5+n9+n8+n4+2
         *n3+2*ep,-1+n4));
*--#] la,noints=5,signs=o,+,+,+,o,o,o,o,o,level=1,stop=0 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* DEN = den(-1+n6)
* MAX = Z(0,0,0,0,0,0,1,0,0)
* MIN = Z(-1,-1,0,0,0,-1,0,0,0)
* nterms = 4 (Z: 4, Y: 0)
* NonIntegerIndices: 5
* Signs: +,+,o,o,o,+,o,o,o
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?pos_,n2?pos_,n3?,n4?,n5?,n6?{>=2},n7?,n8?,n9?) =
*      #include la.rules # la,noints=5,signs=+,+,o,o,o,+,o,o,o,level=1,stop=0
*
*--#[ la,noints=5,signs=+,+,o,o,o,+,o,o,o,level=1,stop=0 :
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,-1+n6,1+n7,n8,n9)*(rat(n7,-1+n6
         ))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,n7,n8,n9)*(rat(1,1))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,n9)*(rat(-n7,-1+
         n6))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,n5,-1+n6,n7,n8,n9)*(rat(-5+n9+n7+n6+2
         *n1+2*ep,-1+n6));
*--#] la,noints=5,signs=+,+,o,o,o,+,o,o,o,level=1,stop=0 : 
* LHS = Z(n1,n2,n3,1,n5,1,n7,n8,n9)
* DEN = den(-4+n9+n7+2*n1+2*ep)*den(-4+n9+n8+2*n3+2*ep)*den(-6+n8+n7+n3+2*n2+n1+2*ep)*den(-5+n9+n8+n7+n3+n2+n1+2*ep)
* MAX = Z(0,0,0,0,0,0,1,1,1)
* MIN = Z(-1,-1,-1,-1,0,-1,-1,-1,0)
* nterms = 18 (Z: 12, Y: 6)
* NonIntegerIndices: 5
* Signs: +,+,+,1,o,1,+,+,-?neg0_
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?pos_,n2?pos_,n3?pos_,1,n5?,1,n7?pos_,n8?pos_,n9?neg0_) =
*      #include la.rules # la,noints=5,signs=+,+,+,1,o,1,+,+,-?neg0_,level=1,stop=1
*
*--#[ la,noints=5,signs=+,+,+,1,o,1,+,+,-?neg0_,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,-1+n2,n3,1,n5,1,n7,n8,1+n9)*(rat(4*n9-n8*n9
         -n7*n9-n5*n9-n2*n9-2*ep*n9-ep*n9*x5,60-12*n9-22*n8+2*n8*n9+2*n8^2-22*
         n7+2*n7*n9+4*n7*n8+2*n7^2-22*n3+2*n3*n9+4*n3*n8+4*n3*n7+2*n3^2-32*n2+4
         *n2*n9+6*n2*n8+6*n2*n7+6*n2*n3+4*n2^2-22*n1+2*n1*n9+4*n1*n8+4*n1*n7+4*
         n1*n3+6*n1*n2+2*n1^2-44*ep+4*ep*n9+8*ep*n8+8*ep*n7+8*ep*n3+12*ep*n2+8*
         ep*n1+8*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,-1+n3,1,n5,1,n7,n8,1+n9)*(rat(9*n9-n9^2-
         2*n8*n9-2*n7*n9-n5*n9-n3*n9-2*n2*n9-n1*n9-4*ep*n9-ep*n9*x5,30-6*n9-11*
         n8+n8*n9+n8^2-11*n7+n7*n9+2*n7*n8+n7^2-11*n3+n3*n9+2*n3*n8+2*n3*n7+
         n3^2-16*n2+2*n2*n9+3*n2*n8+3*n2*n7+3*n2*n3+2*n2^2-11*n1+n1*n9+2*n1*n8+
         2*n1*n7+2*n1*n3+3*n1*n2+n1^2-22*ep+2*ep*n9+4*ep*n8+4*ep*n7+4*ep*n3+6*
         ep*n2+4*ep*n1+4*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,1,n5,1,n7,-1+n8,1+n9)*(rat(-4*n9+n8*
         n9+n7*n9+n5*n9+n2*n9+2*ep*n9+ep*n9*x5,60-12*n9-22*n8+2*n8*n9+2*n8^2-22
         *n7+2*n7*n9+4*n7*n8+2*n7^2-22*n3+2*n3*n9+4*n3*n8+4*n3*n7+2*n3^2-32*n2+
         4*n2*n9+6*n2*n8+6*n2*n7+6*n2*n3+4*n2^2-22*n1+2*n1*n9+4*n1*n8+4*n1*n7+4
         *n1*n3+6*n1*n2+2*n1^2-44*ep+4*ep*n9+8*ep*n8+8*ep*n7+8*ep*n3+12*ep*n2+8
         *ep*n1+8*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,1,n5,1,n7,n8,1+n9)*(rat(-n9,-12+2*n8+
         2*n7+2*n3+4*n2+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,1,n5,1,n7,n8,n9)*(rat(27-12*n9+n9^2-6
         *n8+2*n8*n9-15*n7+3*n7*n9+2*n7*n8+2*n7^2-3*n5+n5*n9+n5*n7-3*n3+n3*n9+
         n3*n7-6*n2+2*n2*n9+2*n2*n7-12*n1+2*n1*n9+2*n1*n8+3*n1*n7+n1*n5+n1*n3+2
         *n1*n2+n1^2-30*ep-3*ep*x5+6*ep*n9+ep*n9*x5+4*ep*n8+8*ep*n7+ep*n7*x5+2*
         ep*n5+2*ep*n3+4*ep*n2+6*ep*n1+ep*n1*x5+8*ep^2+2*ep^2*x5,30-6*n9-11*n8+
         n8*n9+n8^2-11*n7+n7*n9+2*n7*n8+n7^2-11*n3+n3*n9+2*n3*n8+2*n3*n7+n3^2-
         16*n2+2*n2*n9+3*n2*n8+3*n2*n7+3*n2*n3+2*n2^2-11*n1+n1*n9+2*n1*n8+2*n1*
         n7+2*n1*n3+3*n1*n2+n1^2-22*ep+2*ep*n9+4*ep*n8+4*ep*n7+4*ep*n3+6*ep*n2+
         4*ep*n1+4*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,-1+n3,1,n5,1,n7,n8,1+n9)*(rat(4*n9-n8*n9
         -n7*n9-n5*n9-n2*n9-2*ep*n9-ep*n9*x5,60-12*n9-22*n8+2*n8*n9+2*n8^2-22*
         n7+2*n7*n9+4*n7*n8+2*n7^2-22*n3+2*n3*n9+4*n3*n8+4*n3*n7+2*n3^2-32*n2+4
         *n2*n9+6*n2*n8+6*n2*n7+6*n2*n3+4*n2^2-22*n1+2*n1*n9+4*n1*n8+4*n1*n7+4*
         n1*n3+6*n1*n2+2*n1^2-44*ep+4*ep*n9+8*ep*n8+8*ep*n7+8*ep*n3+12*ep*n2+8*
         ep*n1+8*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,1,n5,1,n7,n8,n9)*(rat(48-8*n9-24*n8+2
         *n8*n9+3*n8^2-24*n7+2*n7*n9+6*n7*n8+3*n7^2-20*n5+2*n5*n9+5*n5*n8+5*n5*
         n7+2*n5^2-4*n3+n3*n8+n3*n7+n3*n5-20*n2+2*n2*n9+5*n2*n8+5*n2*n7+4*n2*n5
         +n2*n3+2*n2^2-4*n1+n1*n8+n1*n7+n1*n5+n1*n2-48*ep-20*ep*x5+4*ep*n9+2*ep
         *n9*x5+12*ep*n8+5*ep*n8*x5+12*ep*n7+5*ep*n7*x5+10*ep*n5+4*ep*n5*x5+2*
         ep*n3+ep*n3*x5+10*ep*n2+4*ep*n2*x5+2*ep*n1+ep*n1*x5+12*ep^2+10*ep^2*x5
         +2*ep^2*x5^2,30-6*n9-11*n8+n8*n9+n8^2-11*n7+n7*n9+2*n7*n8+n7^2-11*n3+
         n3*n9+2*n3*n8+2*n3*n7+n3^2-16*n2+2*n2*n9+3*n2*n8+3*n2*n7+3*n2*n3+2*
         n2^2-11*n1+n1*n9+2*n1*n8+2*n1*n7+2*n1*n3+3*n1*n2+n1^2-22*ep+2*ep*n9+4*
         ep*n8+4*ep*n7+4*ep*n3+6*ep*n2+4*ep*n1+4*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,1,n5,1,-1+n7,n8,1+n9)*(rat(-4*n9+n8*
         n9+n7*n9+n5*n9+n2*n9+2*ep*n9+ep*n9*x5,60-12*n9-22*n8+2*n8*n9+2*n8^2-22
         *n7+2*n7*n9+4*n7*n8+2*n7^2-22*n3+2*n3*n9+4*n3*n8+4*n3*n7+2*n3^2-32*n2+
         4*n2*n9+6*n2*n8+6*n2*n7+6*n2*n3+4*n2^2-22*n1+2*n1*n9+4*n1*n8+4*n1*n7+4
         *n1*n3+6*n1*n2+2*n1^2-44*ep+4*ep*n9+8*ep*n8+8*ep*n7+8*ep*n3+12*ep*n2+8
         *ep*n1+8*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,1,n5,1,n7,n8,1+n9)*(rat(-n9,-12+2*n8+
         2*n7+2*n3+4*n2+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,1,n5,1,n7,n8,n9)*(rat(27-12*n9+n9^2-
         15*n8+3*n8*n9+2*n8^2-6*n7+2*n7*n9+2*n7*n8-3*n5+n5*n9+n5*n8-12*n3+2*n3*
         n9+3*n3*n8+2*n3*n7+n3*n5+n3^2-6*n2+2*n2*n9+2*n2*n8+2*n2*n3-3*n1+n1*n9+
         n1*n8+n1*n3-30*ep-3*ep*x5+6*ep*n9+ep*n9*x5+8*ep*n8+ep*n8*x5+4*ep*n7+2*
         ep*n5+6*ep*n3+ep*n3*x5+4*ep*n2+2*ep*n1+8*ep^2+2*ep^2*x5,30-6*n9-11*n8+
         n8*n9+n8^2-11*n7+n7*n9+2*n7*n8+n7^2-11*n3+n3*n9+2*n3*n8+2*n3*n7+n3^2-
         16*n2+2*n2*n9+3*n2*n8+3*n2*n7+3*n2*n3+2*n2^2-11*n1+n1*n9+2*n1*n8+2*n1*
         n7+2*n1*n3+3*n1*n2+n1^2-22*ep+2*ep*n9+4*ep*n8+4*ep*n7+4*ep*n3+6*ep*n2+
         4*ep*n1+4*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,1,n5,1,-1+n7,n8,n9)*(rat(-12+4*n9+3*n8-
         n8*n9+7*n7-n7*n9-n7*n8-n7^2+3*n5-n5*n9-n5*n7+3*n2-n2*n9-n2*n7+4*n1-n1*
         n8-n1*n7-n1*n5-n1*n2+14*ep+3*ep*x5-2*ep*n9-ep*n9*x5-2*ep*n8-4*ep*n7-ep
         *n7*x5-2*ep*n5-2*ep*n2-2*ep*n1-ep*n1*x5-4*ep^2-2*ep^2*x5,30-6*n9-11*n8
         +n8*n9+n8^2-11*n7+n7*n9+2*n7*n8+n7^2-11*n3+n3*n9+2*n3*n8+2*n3*n7+n3^2-
         16*n2+2*n2*n9+3*n2*n8+3*n2*n7+3*n2*n3+2*n2^2-11*n1+n1*n9+2*n1*n8+2*n1*
         n7+2*n1*n3+3*n1*n2+n1^2-22*ep+2*ep*n9+4*ep*n8+4*ep*n7+4*ep*n3+6*ep*n2+
         4*ep*n1+4*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,1,n5,1,n7,-1+n8,n9)*(rat(-12+4*n9+7*n8-
         n8*n9-n8^2+3*n7-n7*n9-n7*n8+3*n5-n5*n9-n5*n8+4*n3-n3*n8-n3*n7-n3*n5+3*
         n2-n2*n9-n2*n8-n2*n3+14*ep+3*ep*x5-2*ep*n9-ep*n9*x5-4*ep*n8-ep*n8*x5-2
         *ep*n7-2*ep*n5-2*ep*n3-ep*n3*x5-2*ep*n2-4*ep^2-2*ep^2*x5,30-6*n9-11*n8
         +n8*n9+n8^2-11*n7+n7*n9+2*n7*n8+n7^2-11*n3+n3*n9+2*n3*n8+2*n3*n7+n3^2-
         16*n2+2*n2*n9+3*n2*n8+3*n2*n7+3*n2*n3+2*n2^2-11*n1+n1*n9+2*n1*n8+2*n1*
         n7+2*n1*n3+3*n1*n2+n1^2-22*ep+2*ep*n9+4*ep*n8+4*ep*n7+4*ep*n3+6*ep*n2+
         4*ep*n1+4*ep^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,0,n5,1,n7,n8,1+n9)*(rat(n9,-12+2*n8+2
         *n7+2*n3+4*n2+2*n1+4*ep))
      +`Y'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,1,n5,0,1+n7,n8,n9)*(rat(-3*n7+n7*n9+
         n7^2+n1*n7+2*ep*n7,24-6*n9-4*n8+n8*n9-10*n7+n7*n9+n7*n8+n7^2-4*n3+n3*
         n9+n3*n7-8*n2+2*n2*n9+2*n2*n7-16*n1+n1*n9+2*n1*n8+3*n1*n7+2*n1*n3+4*n1
         *n2+2*n1^2-20*ep+2*ep*n9+2*ep*n8+4*ep*n7+2*ep*n3+4*ep*n2+6*ep*n1+4*
         ep^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,0,n5,1,n7,1+n8,n9)*(rat(3*n8-n8*n9-
         n8^2-n3*n8-2*ep*n8,24-6*n9-10*n8+n8*n9+n8^2-4*n7+n7*n9+n7*n8-16*n3+n3*
         n9+3*n3*n8+2*n3*n7+2*n3^2-8*n2+2*n2*n9+2*n2*n8+4*n2*n3-4*n1+n1*n9+n1*
         n8+2*n1*n3-20*ep+2*ep*n9+4*ep*n8+2*ep*n7+6*ep*n3+4*ep*n2+2*ep*n1+4*
         ep^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,1,n5,0,1+n7,n8,n9)*(rat(3*n7-n7*n9-
         n7^2-n1*n7-2*ep*n7,24-6*n9-4*n8+n8*n9-10*n7+n7*n9+n7*n8+n7^2-4*n3+n3*
         n9+n3*n7-8*n2+2*n2*n9+2*n2*n7-16*n1+n1*n9+2*n1*n8+3*n1*n7+2*n1*n3+4*n1
         *n2+2*n1^2-20*ep+2*ep*n9+2*ep*n8+4*ep*n7+2*ep*n3+4*ep*n2+6*ep*n1+4*
         ep^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,0,n5,1,n7,1+n8,n9)*(rat(-3*n8+n8*n9+
         n8^2+n3*n8+2*ep*n8,24-6*n9-10*n8+n8*n9+n8^2-4*n7+n7*n9+n7*n8-16*n3+n3*
         n9+3*n3*n8+2*n3*n7+2*n3^2-8*n2+2*n2*n9+2*n2*n8+4*n2*n3-4*n1+n1*n9+n1*
         n8+2*n1*n3-20*ep+2*ep*n9+4*ep*n8+2*ep*n7+6*ep*n3+4*ep*n2+2*ep*n1+4*
         ep^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,1,n5,0,n7,n8,1+n9)*(rat(n9,-12+2*n8+2
         *n7+2*n3+4*n2+2*n1+4*ep));
*--#] la,noints=5,signs=+,+,+,1,o,1,+,+,-?neg0_,level=1,stop=1 : 
*--#] la13,5-rules-a
        #else
*--#[ la13,5-rules-b
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* DEN = den(-4+n9+n7+n6+2*n1+2*ep)*den(-4+n9+n8+n4+n3+n1+2*ep)
* MAX = Z(0,0,0,0,0,1,1,0,1)
* MIN = Z(-1,-1,-1,-1,0,-1,-1,-1,0)
* nterms = 14 (Z: 14, Y: 0)
* NonIntegerIndices: 5
* Signs: +,+,+,o,o,x,x,o,-
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 2
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?pos_,n2?pos_,n3?pos_,n4?,n5?,n6?,n7?,n8?,n9?neg_) =
*      #include la.rules # la,noints=5,signs=+,+,+,o,o,x,x,o,-,level=2,stop=1
*
*--#[ la,noints=5,signs=+,+,+,o,o,x,x,o,-,level=2,stop=1 :
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,-1+n2,n3,n4,n5,n6,1+n7,n8,1+n9)*(rat(n7,-8+
         2*n9+2*n8+2*n4+2*n3+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,-1+n3,n4,n5,1+n6,n7,n8,1+n9)*(rat(n6,-8+
         2*n9+2*n8+2*n4+2*n3+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,-1+n3,n4,n5,n6,1+n7,n8,1+n9)*(rat(n7,-8+
         2*n9+2*n8+2*n4+2*n3+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,-1+n4,n5,1+n6,n7,n8,1+n9)*(rat(-n6,-8
         +2*n9+2*n8+2*n4+2*n3+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,-1+n6,1+n7,n8,1+n9)*(rat(-4*n7+
         n7*n9+n7*n8+n7^2+n3*n7+n2*n7+n1*n7+2*ep*n7,32-16*n9+2*n9^2-8*n8+2*n8*
         n9-8*n7+2*n7*n9+2*n7*n8-8*n6+2*n6*n9+2*n6*n8-8*n4+2*n4*n9+2*n4*n7+2*n4
         *n6-8*n3+2*n3*n9+2*n3*n7+2*n3*n6-24*n1+6*n1*n9+4*n1*n8+2*n1*n7+2*n1*n6
         +4*n1*n4+4*n1*n3+4*n1^2-32*ep+8*ep*n9+4*ep*n8+4*ep*n7+4*ep*n6+4*ep*n4+
         4*ep*n3+12*ep*n1+8*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,1+n6,n7,n8,1+n9)*(rat(n6,-8+2*
         n9+2*n8+2*n4+2*n3+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,1+n6,n7,n8,n9)*(rat(-n6,-4+n9+
         n8+n4+n3+n1+2*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,1+n7,-1+n8,1+n9)*(rat(-n7,-8
         +2*n9+2*n8+2*n4+2*n3+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,1+n7,n8,n9)*(rat(-n7,-4+n9+
         n8+n4+n3+n1+2*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(-16+4*n9-n8
         +12*n7-2*n7*n9-2*n7^2+9*n6-n6*n9+n6*n8-3*n6*n7-2*n6^2+4*n5-n5*n9-n5*n7
         -n5*n6-5*n3+n3*n9+n3*n7+2*n3*n6-n2+n2*n6+13*n1-n1*n9-5*n1*n7-4*n1*n6-2
         *n1*n5+2*n1*n3-2*n1^2+16*ep+4*ep*x5-2*ep*n9-ep*n9*x5-6*ep*n7-ep*n7*x5-
         4*ep*n6-ep*n6*x5-2*ep*n5+2*ep*n3-6*ep*n1-2*ep*n1*x5-4*ep^2-2*ep^2*x5,
         32-16*n9+2*n9^2-8*n8+2*n8*n9-8*n7+2*n7*n9+2*n7*n8-8*n6+2*n6*n9+2*n6*n8
         -8*n4+2*n4*n9+2*n4*n7+2*n4*n6-8*n3+2*n3*n9+2*n3*n7+2*n3*n6-24*n1+6*n1*
         n9+4*n1*n8+2*n1*n7+2*n1*n6+4*n1*n4+4*n1*n3+4*n1^2-32*ep+8*ep*n9+4*ep*
         n8+4*ep*n7+4*ep*n6+4*ep*n4+4*ep*n3+12*ep*n1+8*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,n4,n5,-1+n6,1+n7,n8,1+n9)*(rat(4*n7-
         n7*n9-n7*n8-n7^2-n3*n7-n2*n7-n1*n7-2*ep*n7,32-16*n9+2*n9^2-8*n8+2*n8*
         n9-8*n7+2*n7*n9+2*n7*n8-8*n6+2*n6*n9+2*n6*n8-8*n4+2*n4*n9+2*n4*n7+2*n4
         *n6-8*n3+2*n3*n9+2*n3*n7+2*n3*n6-24*n1+6*n1*n9+4*n1*n8+2*n1*n7+2*n1*n6
         +4*n1*n4+4*n1*n3+4*n1^2-32*ep+8*ep*n9+4*ep*n8+4*ep*n7+4*ep*n6+4*ep*n4+
         4*ep*n3+12*ep*n1+8*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(4-n8-n7-n5-
         n2-2*ep-ep*x5,-8+2*n9+2*n8+2*n4+2*n3+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*(rat(-4+n8+n7+n5
         +n2+2*ep+ep*x5,-8+2*n9+2*n8+2*n4+2*n3+2*n1+4*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(12-7*n9+n9^2-3
         *n8+n8*n9-7*n7+2*n7*n9+n7*n8+n7^2-3*n3+n3*n9+n3*n7-3*n2+n2*n9+n2*n7-11
         *n1+3*n1*n9+2*n1*n8+3*n1*n7+2*n1*n3+2*n1*n2+2*n1^2-14*ep+4*ep*n9+2*ep*
         n8+4*ep*n7+2*ep*n3+2*ep*n2+6*ep*n1+4*ep^2,32-16*n9+2*n9^2-8*n8+2*n8*n9
         -8*n7+2*n7*n9+2*n7*n8-8*n6+2*n6*n9+2*n6*n8-8*n4+2*n4*n9+2*n4*n7+2*n4*
         n6-8*n3+2*n3*n9+2*n3*n7+2*n3*n6-24*n1+6*n1*n9+4*n1*n8+2*n1*n7+2*n1*n6+
         4*n1*n4+4*n1*n3+4*n1^2-32*ep+8*ep*n9+4*ep*n8+4*ep*n7+4*ep*n6+4*ep*n4+4
         *ep*n3+12*ep*n1+8*ep^2));
*--#] la,noints=5,signs=+,+,+,o,o,x,x,o,-,level=2,stop=1 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(2+n8-n7-n2-n1-ep)*den(-5+n8+n7+n3+n2+n1+2*ep)
* MAX = Z(0,0,1,0,1,0,1,0,0)
* MIN = Z(-1,-1,0,0,0,0,-1,-1,0)
* nterms = 10 (Z: 10, Y: 0)
* NonIntegerIndices: 5
* Signs: +,+,x,o,x,o,+,+,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?pos_,n2?pos_,n3?,n4?,n5?,n6?,n7?pos_,n8?pos_,0) =
*      #include la.rules # la,noints=5,signs=+,+,x,o,x,o,+,+,0,level=1,stop=0
*
*--#[ la,noints=5,signs=+,+,x,o,x,o,+,+,0,level=1,stop=0 :
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,1+n7,-1+n8,0)*(rat(4*n7-n7*
         n8-n4*n7-n3*n7-2*ep*n7,-10-3*n8+n8^2+7*n7-n7^2+2*n3+n3*n8-n3*n7+7*n2-2
         *n2*n7-n2*n3-n2^2+7*n1-2*n1*n7-n1*n3-2*n1*n2-n1^2+9*ep+ep*n8-3*ep*n7-
         ep*n3-3*ep*n2-3*ep*n1-2*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,n7,n8,0)*(rat(36-4*n8-17*n7+
         n7*n8+2*n7^2-13*n6+n6*n8+3*n6*n7+n6^2-4*n5+n5*n7+n5*n6-4*n4+n4*n7+n4*
         n6-4*n3+n3*n7+n3*n6-4*n2+n2*n7+n2*n6-13*n1+n1*n8+3*n1*n7+2*n1*n6+n1*n5
         +n1*n4+n1*n3+n1*n2+n1^2-34*ep-4*ep*x5+2*ep*n8+8*ep*n7+ep*n7*x5+6*ep*n6
         +ep*n6*x5+2*ep*n5+2*ep*n4+2*ep*n3+2*ep*n2+6*ep*n1+ep*n1*x5+8*ep^2+2*
         ep^2*x5,10+3*n8-n8^2-7*n7+n7^2-2*n3-n3*n8+n3*n7-7*n2+2*n2*n7+n2*n3+
         n2^2-7*n1+2*n1*n7+n1*n3+2*n1*n2+n1^2-9*ep-ep*n8+3*ep*n7+ep*n3+3*ep*n2+
         3*ep*n1+2*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,1+n3,n4,n5,n6,n7,n8,0)*(rat(n3,2+n8-n7-
         n2-n1-ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,n4,1+n5,n6,n7,-1+n8,0)*(rat(-4*n5+n5*
         n8+n4*n5+n3*n5-4*ep*x5+ep*n8*x5+2*ep*n5+ep*n4*x5+ep*n3*x5+2*ep^2*x5,-
         10-3*n8+n8^2+7*n7-n7^2+2*n3+n3*n8-n3*n7+7*n2-2*n2*n7-n2*n3-n2^2+7*n1-2
         *n1*n7-n1*n3-2*n1*n2-n1^2+9*ep+ep*n8-3*ep*n7-ep*n3-3*ep*n2-3*ep*n1-2*
         ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,n4,n5,n6,1+n7,-1+n8,0)*(rat(-4*n7+n7*
         n8+n4*n7+n3*n7+2*ep*n7,-10-3*n8+n8^2+7*n7-n7^2+2*n3+n3*n8-n3*n7+7*n2-2
         *n2*n7-n2*n3-n2^2+7*n1-2*n1*n7-n1*n3-2*n1*n2-n1^2+9*ep+ep*n8-3*ep*n7-
         ep*n3-3*ep*n2-3*ep*n1-2*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,n4,n5,n6,n7,n8,0)*(rat(52-12*n8-21*n7
         +3*n7*n8+2*n7^2-6*n6+n6*n8+n6*n7-15*n5+2*n5*n8+3*n5*n7+n5*n6+n5^2-15*
         n4+2*n4*n8+3*n4*n7+n4*n6+2*n4*n5+n4^2-15*n3+2*n3*n8+3*n3*n7+n3*n6+2*n3
         *n5+2*n3*n4+n3^2-15*n2+2*n2*n8+3*n2*n7+n2*n6+2*n2*n5+2*n2*n4+2*n2*n3+
         n2^2-6*n1+n1*n8+n1*n7+n1*n5+n1*n4+n1*n3+n1*n2-50*ep-15*ep*x5+6*ep*n8+2
         *ep*n8*x5+10*ep*n7+3*ep*n7*x5+3*ep*n6+ep*n6*x5+7*ep*n5+2*ep*n5*x5+7*ep
         *n4+2*ep*n4*x5+7*ep*n3+2*ep*n3*x5+7*ep*n2+2*ep*n2*x5+3*ep*n1+ep*n1*x5+
         12*ep^2+7*ep^2*x5+ep^2*x5^2,10+3*n8-n8^2-7*n7+n7^2-2*n3-n3*n8+n3*n7-7*
         n2+2*n2*n7+n2*n3+n2^2-7*n1+2*n1*n7+n1*n3+2*n1*n2+n1^2-9*ep-ep*n8+3*ep*
         n7+ep*n3+3*ep*n2+3*ep*n1+2*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,1+n3,n4,n5,n6,n7,-1+n8,0)*(rat(-n3,2+n8-n7-
         n2-n1-ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,1+n5,n6,n7,-1+n8,0)*(rat(4*n5-n5*n8-
         n4*n5-n3*n5+4*ep*x5-ep*n8*x5-2*ep*n5-ep*n4*x5-ep*n3*x5-2*ep^2*x5,-10-3
         *n8+n8^2+7*n7-n7^2+2*n3+n3*n8-n3*n7+7*n2-2*n2*n7-n2*n3-n2^2+7*n1-2*n1*
         n7-n1*n3-2*n1*n2-n1^2+9*ep+ep*n8-3*ep*n7-ep*n3-3*ep*n2-3*ep*n1-2*ep^2)
         )
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,n5,n6,-1+n7,n8,0)*(rat(-24+4*n8+10*n7
         -n7*n8-n7^2+6*n6-n6*n8-n6*n7+4*n5-n5*n7-n5*n6+4*n4-n4*n7-n4*n6+4*n3-n3
         *n7-n3*n6+4*n2-n2*n7-n2*n6+6*n1-n1*n8-n1*n7-n1*n5-n1*n4-n1*n3-n1*n2+24
         *ep+4*ep*x5-2*ep*n8-5*ep*n7-ep*n7*x5-3*ep*n6-ep*n6*x5-2*ep*n5-2*ep*n4-
         2*ep*n3-2*ep*n2-3*ep*n1-ep*n1*x5-6*ep^2-2*ep^2*x5,10+3*n8-n8^2-7*n7+
         n7^2-2*n3-n3*n8+n3*n7-7*n2+2*n2*n7+n2*n3+n2^2-7*n1+2*n1*n7+n1*n3+2*n1*
         n2+n1^2-9*ep-ep*n8+3*ep*n7+ep*n3+3*ep*n2+3*ep*n1+2*ep^2))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,n5,n6,n7,-1+n8,0)*(rat(28-11*n8+n8^2-
         4*n7+n7*n8-4*n5+n5*n8-11*n4+2*n4*n8+n4*n7+n4*n5+n4^2-11*n3+2*n3*n8+n3*
         n7+n3*n5+2*n3*n4+n3^2-8*n2+2*n2*n8+2*n2*n4+2*n2*n3-26*ep-4*ep*x5+5*ep*
         n8+ep*n8*x5+2*ep*n7+2*ep*n5+5*ep*n4+ep*n4*x5+5*ep*n3+ep*n3*x5+4*ep*n2+
         6*ep^2+2*ep^2*x5,-10-3*n8+n8^2+7*n7-n7^2+2*n3+n3*n8-n3*n7+7*n2-2*n2*n7
         -n2*n3-n2^2+7*n1-2*n1*n7-n1*n3-2*n1*n2-n1^2+9*ep+ep*n8-3*ep*n7-ep*n3-3
         *ep*n2-3*ep*n1-2*ep^2));
*--#] la,noints=5,signs=+,+,x,o,x,o,+,+,0,level=1,stop=0 : 
*--#] la13,5-rules-b
        #endif
        goto 2;
label 1;
        #call ReduceZ2Y(Z,Y,9,ifmatch->2,1,2,3,4,6,7,8)
        $doloop = 1;
label 2;
        ModuleOption maximum, $doloop;
        .sort:`TOPO'-`nloop++';
        #if `$doloop'
          #redefine loop "0"
        #endif
      #enddo
      #break
*--#] la13,5
*--#[ no,none
    #case no,none
      #define nloop "1"
      #do loop=1,1
        #$doloop = 0;
        if (count(`Z',1) == 0) goto 2;
*--#[ no,none-rules
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* DEN = den(-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
* MAX = Z(1,0,1,0,0,0,0,0,2)
* MIN = Z(0,-1,0,0,0,0,-1,-1,1)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: none
* Signs: x,x,x,o,o,o,x,x,-*?neg_
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 2
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,n1?,n2?,n3?,n4?,n5?,n6?,n7?,n8?,n9?neg_) =
*      #include no.rules # no,noints=none,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0
*
*--#[ no,noints=none,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0 :
      +`Z'(0,0,0,0,0,0,0,0,0,1+n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(n1,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep))
      +`Z'(0,0,0,0,0,0,0,0,0,1+n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*(rat(-n1,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,-1+n2,1+n3,n4,n5,n6,n7,n8,1+n9)*(rat(n3,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,-1+n2,n3,n4,n5,n6,n7,n8,2+n9)*(rat(1+n9,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,1+n3,n4,n5,n6,n7,-1+n8,1+n9)*(rat(-n3,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(-3+n9+n3+2*n2+
         n1+2*ep,-14+2*n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep));
*--#] no,noints=none,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(-1+n6)
* MAX = Z(0,0,0,0,1,0,0,1,0)
* MIN = Z(-1,-1,-1,-1,0,-1,-1,0,0)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: none
* Signs: +,x,+,+,x,+,x,x,0
* Order: 6,1,2,3,4,5,6,7,8
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,n1?pos_,n2?,n3?pos_,n4?pos_,n5?,n6?{>=2},n7?,n8?,0) =
*      #include no.rules # no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=6,1,2,3,4,5,6,7,8,level=1,stop=0
*
*--#[ no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=6,1,2,3,4,5,6,7,8,level=1,stop=0 :
      +`Z'(0,0,0,0,0,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,n7,n8,0)*(rat(1,1))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,0)*(rat(n8,-1+n6))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,0)*(rat(-n8,-1+n6)
         )
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,0)*(rat(-n5,-1+n6)
         )
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,n3,n4,1+n5,-1+n6,-1+n7,n8,0)*(rat(n5,-1+n6))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,n3,n4,n5,-1+n6,n7,n8,0)*(rat(-9+n8+2*n7+n6+
         n5+2*n2+2*n1+4*ep,-1+n6));
*--#] no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=6,1,2,3,4,5,6,7,8,level=1,stop=0 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(-1+n4)
* MAX = Z(0,0,0,0,1,0,1,0,0)
* MIN = Z(-1,-1,-1,-1,0,-1,0,-1,0)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: none
* Signs: +,x,+,+,x,+,x,x,0
* Order: 4,1,2,3,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,n1?pos_,n2?,n3?pos_,n4?{>=2},n5?,n6?pos_,n7?,n8?,0) =
*      #include no.rules # no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=4,1,2,3,5,6,7,8,9,level=1,stop=1
*
*--#[ no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=4,1,2,3,5,6,7,8,9,level=1,stop=1 :
      +`Z'(0,0,0,0,0,0,0,0,0,-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,0)*(rat(-n7,-1+n4)
         )
      +`Z'(0,0,0,0,0,0,0,0,0,n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,0)*(rat(n7,-1+n4))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,-1+n3,n4,n5,n6,n7,n8,0)*(rat(1,1))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,0)*(rat(-n5,-1+n4)
         )
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,n3,-1+n4,1+n5,n6,n7,-1+n8,0)*(rat(n5,-1+n4))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,n3,-1+n4,n5,n6,n7,n8,0)*(rat(-9+2*n8+n7+n5+
         n4+2*n3+2*n2+4*ep,-1+n4));
*--#] no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=4,1,2,3,5,6,7,8,9,level=1,stop=1 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(-1+n3)
* MAX = Z(0,1,0,0,0,0,0,1,0)
* MIN = Z(-1,0,-1,-1,-1,-1,-1,0,0)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: none
* Signs: +,x,+,+,x,+,x,x,0
* Order: 3,1,2,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,n1?pos_,n2?,n3?{>=2},n4?pos_,n5?,n6?pos_,n7?,n8?,0) =
*      #include no.rules # no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=3,1,2,4,5,6,7,8,9,level=1,stop=1
*
*--#[ no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=3,1,2,4,5,6,7,8,9,level=1,stop=1 :
      +`Z'(0,0,0,0,0,0,0,0,0,-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,0)*(rat(-n2,-1+n3)
         )
      +`Z'(0,0,0,0,0,0,0,0,0,n1,1+n2,-1+n3,n4,n5,n6,-1+n7,n8,0)*(rat(n2,-1+n3))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,0)*(rat(n8,-1+n3))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,0)*(rat(-n8,-1+n3)
         )
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,-1+n3,n4,n5,n6,n7,n8,0)*(rat(-9+n8+2*n7+2*n5
         +2*n4+n3+n2+4*ep,-1+n3))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,n3,-1+n4,n5,n6,n7,n8,0)*(rat(1,1));
*--#] no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=3,1,2,4,5,6,7,8,9,level=1,stop=1 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(-1+n1)
* MAX = Z(0,1,0,0,0,0,1,0,0)
* MIN = Z(-1,0,-1,-1,-1,-1,0,-1,0)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: none
* Signs: +,x,+,+,x,+,x,x,0
* Order: 1,2,4,5,6,3,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,n1?{>=2},n2?,n3?pos_,n4?pos_,n5?,n6?pos_,n7?,n8?,0) =
*      #include no.rules # no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=1,2,4,5,6,3,7,8,9,level=1,stop=1
*
*--#[ no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=1,2,4,5,6,3,7,8,9,level=1,stop=1 :
      +`Z'(0,0,0,0,0,0,0,0,0,-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,0)*(rat(-n2,-1+n1)
         )
      +`Z'(0,0,0,0,0,0,0,0,0,-1+n1,1+n2,n3,n4,n5,n6,n7,-1+n8,0)*(rat(n2,-1+n1))
      +`Z'(0,0,0,0,0,0,0,0,0,-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,0)*(rat(-n7,-1+n1)
         )
      +`Z'(0,0,0,0,0,0,0,0,0,-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,0)*(rat(n7,-1+n1))
      +`Z'(0,0,0,0,0,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,n7,n8,0)*(rat(-9+2*n8+n7+2*n6
         +2*n5+n2+n1+4*ep,-1+n1))
      +`Z'(0,0,0,0,0,0,0,0,0,n1,n2,n3,n4,n5,-1+n6,n7,n8,0)*(rat(1,1));
*--#] no,noints=none,signs=+,x,+,+,x,+,x,x,0,order=1,2,4,5,6,3,7,8,9,level=1,stop=1 : 
* LHS = Z(1,n2,1,1,n5,1,n7,n8,0)
* DEN = den(-1+n5)*den(-1+n8)
* MAX = Z(0,1,0,0,1,0,1,1,0)
* MIN = Z(-1,-1,-1,-1,-1,-1,-1,-1,0)
* nterms = 16 (Z: 8, Y: 8)
* NonIntegerIndices: none
* Signs: 1,+,1,1,+,1,+,+,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,1,n2?pos_,1,1,n5?{>=2},1,n7?pos_,n8?{>=2},0) =
*      #include no.rules # no,noints=none,signs=1,+,1,1,+,1,+,+,0,level=1,stop=0
*
*--#[ no,noints=none,signs=1,+,1,1,+,1,+,+,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,-1+n5,1,1+n7,n8,0)*(rat(-n7,-1+n5))
      +`Z'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,-1+n5,1,n7,1+n8,0)*(rat(-n8,-1+n5))
      +`Z'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,n5,1,n7,n8,0)*(rat(-1,1))
      +`Z'(0,0,0,0,0,0,0,0,0,1,1+n2,1,1,n5,1,-1+n7,-1+n8,0)*(rat(-n2,-1+n8))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,-1+n5,1,n7,n8,0)*(rat(5-n8-n7-n5-2*n2-2*
         ep,-1+n5))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1+n5,1,-1+n7,-1+n8,0)*(rat(-n5,-1+n8))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,n5,1,-1+n7,n8,0)*(rat(-1,1))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,n5,1,n7,-1+n8,0)*(rat(5-n8-2*n7-n5-n2-2*
         ep,-1+n8))
      +`Y'(0,0,0,0,0,0,0,0,0,0,1+n2,1,1,n5,1,n7,-1+n8,0)*(rat(n2,-1+n8))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,1,-1+n5,1,1+n7,n8,0)*(rat(n7,-1+n5))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,1,n5,1,n7,n8,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,-1+n5,1,n7,1+n8,0)*(rat(n8,-1+n5))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,n5,1,n7,n8,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,1+n5,1,n7,-1+n8,0)*(rat(n5,-1+n8))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,n5,1,n7,n8,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,n5,0,n7,n8,0)*(rat(1,1));
*--#] no,noints=none,signs=1,+,1,1,+,1,+,+,0,level=1,stop=0 : 
* LHS = Z(1,n2,1,1,1,1,n7,n8,0)
* DEN = den(-1+n8)*den(-1+n7+ep)
* MAX = Z(1,1,0,0,1,0,1,1,0)
* MIN = Z(-1,-1,-1,-1,-1,-1,-1,-1,0)
* nterms = 24 (Z: 5, Y: 19)
* NonIntegerIndices: none
* Signs: 1,+,1,1,1,1,+,+,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,1,n2?pos_,1,1,1,1,n7?pos_,n8?{>=2},0) =
*      #include no.rules # no,noints=none,signs=1,+,1,1,1,1,+,+,0,level=1,stop=0
*
*--#[ no,noints=none,signs=1,+,1,1,1,1,+,+,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,1,1,n7,n8,0)*(rat(-1+n8+ep,-1+n7+ep))
      +`Z'(0,0,0,0,0,0,0,0,0,1,1+n2,1,1,1,1,-1+n7,-1+n8,0)*(rat(6*n2-n2*n8-3*n2
         *n7-2*n2^2-4*ep*n2,2-2*n8-2*n7+2*n7*n8-2*ep+2*ep*n8))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1,1,-1+n7,n8,0)*(rat(3-2*n7-n2-4*ep,-2+2*
         n7+2*ep))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1,1,n7,-1+n8,0)*(rat(-17+2*n8+21*n7-2*n7*
         n8-6*n7^2+11*n2-7*n2*n7-2*n2^2+24*ep-2*ep*n8-14*ep*n7-8*ep*n2-8*ep^2,2
         -2*n8-2*n7+2*n7*n8-2*ep+2*ep*n8))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,2,1,-1+n7,-1+n8,0)*(rat(4+n8-3*n7-2*n2-4*
         ep,2-2*n8-2*n7+2*n7*n8-2*ep+2*ep*n8))
      +`Y'(0,0,0,0,0,0,0,0,0,0,1+n2,0,1,1,1,n7,n8,0)*(rat(-n2,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,0,1+n2,1,1,1,1,n7,-1+n8,0)*(rat(-6*n2+n2*n8+3*n2*
         n7+2*n2^2+4*ep*n2,2-2*n8-2*n7+2*n7*n8-2*ep+2*ep*n8))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,0,1,1,1+n7,n8,0)*(rat(n7,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,1,1,1,n7,n8,0)*(rat(-2+n7+n2+2*ep,-2+2*n7+2
         *ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,0,1,1,1+n7,n8,0)*(rat(-n7,-1+n7+ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,0,1,1+n7,n8,0)*(rat(n7,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,0,1,n7,1+n8,0)*(rat(-n8,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,1,0,n7,1+n8,0)*(rat(n8,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,1+n2,0,1,1,1,-1+n7,n8,0)*(rat(n2,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,0,1,n7,1+n8,0)*(rat(n8,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,1,0,n7,1+n8,0)*(rat(-n8,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,1,1,n7,n8,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,1,1,n7,n8,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,2,1,n7,-1+n8,0)*(rat(-5+3*n7+2*n2+4*ep,2-
         2*n8-2*n7+2*n7*n8-2*ep+2*ep*n8))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,0,1,n7,n8,0)*(rat(6-n8-3*n7-2*n2-4*ep,-2+
         2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1,0,n7,n8,0)*(rat(-5+n8+2*n7+2*n2+4*ep,-2
         +2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,2,0,-1+n7,n8,0)*(rat(-1,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,2,-1+n2,1,1,1,0,n7,n8,0)*(rat(1,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,2,n2,1,1,1,0,-1+n7,n8,0)*(rat(-1,-2+2*n7+2*ep));
*--#] no,noints=none,signs=1,+,1,1,1,1,+,+,0,level=1,stop=0 : 
* LHS = Z(1,n2,1,1,n5,1,n7,1,0)
* DEN = den(-1+n5)*den(-3+n7+2*n2+2*ep)
* MAX = Z(1,1,0,0,1,0,1,1,0)
* MIN = Z(-1,-1,-1,-1,-1,-1,-1,-1,0)
* nterms = 23 (Z: 5, Y: 18)
* NonIntegerIndices: none
* Signs: 1,+,1,1,+,1,+,1,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,1,n2?pos_,1,1,n5?{>=2},1,n7?pos_,1,0) =
*      #include no.rules # no,noints=none,signs=1,+,1,1,+,1,+,1,0,level=1,stop=0
*
*--#[ no,noints=none,signs=1,+,1,1,+,1,+,1,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,-1+n5,1,1+n7,1,0)*(rat(7*n7-3*n7^2-2*
         n5*n7-2*n2*n7-4*ep*n7,3-n7-3*n5+n5*n7-2*n2+2*n2*n5-2*ep+2*ep*n5))
      +`Z'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,-1+n5,1,n7,2,0)*(rat(5-3*n7-2*n2-4*ep,
         3-n7-3*n5+n5*n7-2*n2+2*n2*n5-2*ep+2*ep*n5))
      +`Z'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,n5,1,n7,1,0)*(rat(7-3*n7-2*n5-2*n2-6*
         ep,-3+n7+2*n2+2*ep))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,-1+n5,1,n7,1,0)*(rat(6-3*n7-n5-2*n2-4*ep,
         -1+n5))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,n5,1,-1+n7,1,0)*(rat(1-n7+n5-n2,-3+n7+2*
         n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,0,1+n2,0,1,n5,1,n7,1,0)*(rat(n2,-3+n7+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,0,n5,1,1+n7,1,0)*(rat(-n7,-3+n7+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,1,-1+n5,1,1+n7,1,0)*(rat(-6*n7+3*n7^2+n5*n7
         +2*n2*n7+4*ep*n7,3-n7-3*n5+n5*n7-2*n2+2*n2*n5-2*ep+2*ep*n5))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,1,n5,1,n7,1,0)*(rat(-3+2*n7+n2+2*ep,-3+n7+2
         *n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,0,n5,1,1+n7,1,0)*(rat(2*n7,-3+n7+2*n2+2*
         ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,n5,0,n7,2,0)*(rat(-1,-3+n7+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,1+n2,0,1,n5,1,-1+n7,1,0)*(rat(-n2,-3+n7+2*n2+2*
         ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,-1+n5,1,n7,2,0)*(rat(-5+3*n7+2*n2+4*ep,3-
         n7-3*n5+n5*n7-2*n2+2*n2*n5-2*ep+2*ep*n5))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,n5,0,n7,2,0)*(rat(1,-3+n7+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,n5,1,n7,1,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,1+n5,1,n7,0,0)*(rat(n5,-3+n7+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,n5,1,n7,1,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1+n5,0,-1+n7,1,0)*(rat(n5,-3+n7+2*n2+2*ep
         ))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1+n5,1,-1+n7,0,0)*(rat(-2*n5,-3+n7+2*n2+2
         *ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,n5,0,n7,1,0)*(rat(-1+n7,-3+n7+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,n5,1,n7,0,0)*(rat(7-3*n7-n5-3*n2-4*ep,-3+
         n7+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,2,-1+n2,1,1,n5,0,n7,1,0)*(rat(-1,-3+n7+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,2,n2,1,1,n5,0,-1+n7,1,0)*(rat(1,-3+n7+2*n2+2*ep));
*--#] no,noints=none,signs=1,+,1,1,+,1,+,1,0,level=1,stop=0 : 
* LHS = Z(1,n2,1,1,1,1,n7,1,0)
* DEN = den(-1+n7)*den(-1+n7+ep)
* MAX = Z(1,1,0,0,1,0,1,1,0)
* MIN = Z(-1,-1,-1,-1,-1,-1,-1,-1,0)
* nterms = 23 (Z: 2, Y: 21)
* NonIntegerIndices: none
* Signs: 1,+,1,1,1,1,+,1,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,1,n2?pos_,1,1,1,1,n7?{>=2},1,0) =
*      #include no.rules # no,noints=none,signs=1,+,1,1,1,1,+,1,0,level=1,stop=0
*
*--#[ no,noints=none,signs=1,+,1,1,1,1,+,1,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,1,1,n7,1,0)*(rat(ep,-1+n7+ep))
      +`Z'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1,1,-1+n7,1,0)*(rat(-4+4*n7-n7^2+4*n2-2*
         n2*n7-n2^2+9*ep-5*ep*n7-4*ep*n2-4*ep^2,1-2*n7+n7^2-ep+ep*n7))
      +`Y'(0,0,0,0,0,0,0,0,0,0,1+n2,0,1,1,1,n7,1,0)*(rat(-n2,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,0,1,1,1+n7,1,0)*(rat(n7,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,1,1,1,n7,1,0)*(rat(-2+n7+n2+2*ep,-2+2*n7+2*
         ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,0,1,1,1+n7,1,0)*(rat(-n7,-1+n7+ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,0,1,1+n7,1,0)*(rat(n7,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,0,1,n7,2,0)*(rat(-1,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,1,0,n7,2,0)*(rat(1,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,1+n2,0,1,1,1,-1+n7,1,0)*(rat(-3*n2+2*n2*n7+n2^2+
         2*ep*n2,1-2*n7+n7^2-ep+ep*n7))
      +`Y'(0,0,0,0,0,0,0,0,0,1,1+n2,1,1,1,1,-1+n7,0,0)*(rat(5*n2-3*n2*n7-2*n2^2
         -4*ep*n2,2-4*n7+2*n7^2-2*ep+2*ep*n7))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,0,1,n7,2,0)*(rat(1,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,1,0,n7,2,0)*(rat(-1,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,1,1,n7,1,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,1,1,n7,1,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,2,1,n7,0,0)*(rat(-1,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,0,1,n7,1,0)*(rat(5-3*n7-2*n2-4*ep,-2+2*n7
         +2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1,0,n7,1,0)*(rat(-2+n7+n2+2*ep,-1+n7+ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1,1,n7,0,0)*(rat(-1+n2,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,2,0,-1+n7,1,0)*(rat(-2+n7+n2+2*ep,1-2*n7+
         n7^2-ep+ep*n7))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,2,1,-1+n7,0,0)*(rat(3-n7-2*n2-4*ep,2-4*n7
         +2*n7^2-2*ep+2*ep*n7))
      +`Y'(0,0,0,0,0,0,0,0,0,2,-1+n2,1,1,1,0,n7,1,0)*(rat(1,-2+2*n7+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,2,n2,1,1,1,0,-1+n7,1,0)*(rat(-1,-2+2*n7+2*ep));
*--#] no,noints=none,signs=1,+,1,1,1,1,+,1,0,level=1,stop=0 : 
* LHS = Z(1,n2,1,1,1,1,1,1,0)
* DEN = den(-1+n2)*den(-1+n2+ep)
* MAX = Z(1,1,0,0,1,0,1,1,0)
* MIN = Z(-1,-1,-1,-1,-1,-1,-1,-1,0)
* nterms = 20 (Z: 1, Y: 19)
* NonIntegerIndices: none
* Signs: 1,+,1,1,1,1,1,1,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,0,0,1,n2?{>=2},1,1,1,1,1,1,0) =
*      #include no.rules # no,noints=none,signs=1,+,1,1,1,1,1,1,0,level=1,stop=0
*
*--#[ no,noints=none,signs=1,+,1,1,1,1,1,1,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,1,1,1,1,0)*(rat(1-n2-4*ep,-1+n2))
      +`Y'(0,0,0,0,0,0,0,0,0,0,1+n2,0,1,1,1,1,1,0)*(rat(n2,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,0,1,1,2,1,0)*(rat(-1,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,0,n2,1,1,1,1,1,1,0)*(rat(-1+n2+2*ep,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,0,1,1,2,1,0)*(rat(2,-1+n2))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,0,1,1,2,0)*(rat(1-n2-4*ep,2-4*n2+2*
         n2^2-2*ep+2*ep*n2))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,0,1,2,1,0)*(rat(3-3*n2-4*ep,2-4*n2+2*
         n2^2-2*ep+2*ep*n2))
      +`Y'(0,0,0,0,0,0,0,0,0,1,-1+n2,1,1,1,0,1,2,0)*(rat(-1+n2+4*ep,2-4*n2+2*
         n2^2-2*ep+2*ep*n2))
      +`Y'(0,0,0,0,0,0,0,0,0,1,1+n2,0,1,1,1,0,1,0)*(rat(-n2,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,0,1,1,2,0)*(rat(-1,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,1,0,1,2,0)*(rat(1,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,0,1,1,1,1,1,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,1,1,1,1,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,0,2,1,1,0,0)*(rat(1,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1,1,0,1,0)*(rat(1-n2,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,1,1,1,0,0)*(rat(3-3*n2-4*ep,-2+2*n2+2*ep)
         )
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,2,0,0,1,0)*(rat(1,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,1,n2,1,1,2,1,0,0,0)*(rat(-1,-1+n2+ep))
      +`Y'(0,0,0,0,0,0,0,0,0,2,-1+n2,1,1,1,0,1,1,0)*(rat(-1,-2+2*n2+2*ep))
      +`Y'(0,0,0,0,0,0,0,0,0,2,n2,1,1,1,0,0,1,0)*(rat(1,-2+2*n2+2*ep));
*--#] no,noints=none,signs=1,+,1,1,1,1,1,1,0,level=1,stop=0 : 
*--#] no,none-rules
        goto 2;
label 1;
        #call ReduceZ2Y(Z,Y,9,ifmatch->2,1,2,3,4,5,6,7,8)
        $doloop = 1;
label 2;
        ModuleOption maximum, $doloop;
        .sort:`TOPO'-`nloop++';
        #if `$doloop'
          #redefine loop "0"
        #endif
      #enddo
      #break
*--#] no,none
*--#[ no,5
    #case no,5
      #define nloop "1"
      #do loop=1,1
        #$doloop = 0;
        if (count(`Z',1) == 0) goto 2;
*--#[ no,5-rules
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* DEN = den(-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
* MAX = Z(1,0,1,0,0,0,0,0,2)
* MIN = Z(0,-1,0,0,0,0,-1,-1,1)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: 5
* Signs: x,x,x,o,o,o,x,x,-*?neg_
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 2
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?,n2?,n3?,n4?,n5?,n6?,n7?,n8?,n9?neg_) =
*      #include no.rules # no,noints=5,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0
*
*--#[ no,noints=5,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0 :
      +`Z'(0,0,0,0,x5,0,0,0,0,1+n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(n1,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x5))
      +`Z'(0,0,0,0,x5,0,0,0,0,1+n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*(rat(-n1,-14+2
         *n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x5))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,1+n3,n4,n5,n6,n7,n8,1+n9)*(rat(n3,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x5))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,n4,n5,n6,n7,n8,2+n9)*(rat(1+n9,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x5))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,1+n3,n4,n5,n6,n7,-1+n8,1+n9)*(rat(-n3,-14+2
         *n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x5))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(-3+n9+n3+2*n2+
         n1+2*ep,-14+2*n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x5)
         );
*--#] no,noints=5,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(-1+n6)
* MAX = Z(0,0,0,0,1,0,0,1,0)
* MIN = Z(-1,-1,-1,-1,0,-1,-1,0,0)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: 5
* Signs: +,x,+,+,x,+,x,x,0
* Order: 6,1,2,3,4,5,6,7,8
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?pos_,n2?,n3?pos_,n4?pos_,n5?,n6?{>=2},n7?,n8?,0) =
*      #include no.rules # no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=6,1,2,3,4,5,6,7,8,level=1,stop=0
*
*--#[ no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=6,1,2,3,4,5,6,7,8,level=1,stop=0 :
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,n7,n8,0)*(rat(1,1))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,n4,n5,-1+n6,n7,1+n8,0)*(rat(n8,-1+n6)
         )
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,0)*(rat(-n8,-1+n6
         ))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,0)*(rat(-n5-ep*x5
         ,-1+n6))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,1+n5,-1+n6,-1+n7,n8,0)*(rat(n5+ep*x5,
         -1+n6))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,n5,-1+n6,n7,n8,0)*(rat(-9+n8+2*n7+n6+
         n5+2*n2+2*n1+4*ep+ep*x5,-1+n6));
*--#] no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=6,1,2,3,4,5,6,7,8,level=1,stop=0 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(-1+n4)
* MAX = Z(0,0,0,0,1,0,1,0,0)
* MIN = Z(-1,-1,-1,-1,0,-1,0,-1,0)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: 5
* Signs: +,x,+,+,x,+,x,x,0
* Order: 4,1,2,3,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?pos_,n2?,n3?pos_,n4?{>=2},n5?,n6?pos_,n7?,n8?,0) =
*      #include no.rules # no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=4,1,2,3,5,6,7,8,9,level=1,stop=1
*
*--#[ no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=4,1,2,3,5,6,7,8,9,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,0)*(rat(-n7,-1+n4
         ))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,0)*(rat(n7,-1+n4)
         )
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,n4,n5,n6,n7,n8,0)*(rat(1,1))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,0)*(rat(-n5-ep*x5
         ,-1+n4))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,-1+n4,1+n5,n6,n7,-1+n8,0)*(rat(n5+ep*x5,
         -1+n4))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,-1+n4,n5,n6,n7,n8,0)*(rat(-9+2*n8+n7+n5+
         n4+2*n3+2*n2+4*ep+ep*x5,-1+n4));
*--#] no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=4,1,2,3,5,6,7,8,9,level=1,stop=1 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(-1+n3)
* MAX = Z(0,1,0,0,0,0,0,1,0)
* MIN = Z(-1,0,-1,-1,-1,-1,-1,0,0)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: 5
* Signs: +,x,+,+,x,+,x,x,0
* Order: 3,1,2,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?pos_,n2?,n3?{>=2},n4?pos_,n5?,n6?pos_,n7?,n8?,0) =
*      #include no.rules # no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=3,1,2,4,5,6,7,8,9,level=1,stop=1
*
*--#[ no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=3,1,2,4,5,6,7,8,9,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,0)*(rat(-n2,-1+n3
         ))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,1+n2,-1+n3,n4,n5,n6,-1+n7,n8,0)*(rat(n2,-1+n3)
         )
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,n4,-1+n5,n6,n7,1+n8,0)*(rat(n8,-1+n3)
         )
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,n4,n5,-1+n6,n7,1+n8,0)*(rat(-n8,-1+n3
         ))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,-1+n3,n4,n5,n6,n7,n8,0)*(rat(-9+n8+2*n7+2*
         n5+2*n4+n3+n2+4*ep+2*ep*x5,-1+n3))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,-1+n4,n5,n6,n7,n8,0)*(rat(1,1));
*--#] no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=3,1,2,4,5,6,7,8,9,level=1,stop=1 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(-1+n1)
* MAX = Z(0,1,0,0,0,0,1,0,0)
* MIN = Z(-1,0,-1,-1,-1,-1,0,-1,0)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: 5
* Signs: +,x,+,+,x,+,x,x,0
* Order: 1,2,4,5,6,3,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,n1?{>=2},n2?,n3?pos_,n4?pos_,n5?,n6?pos_,n7?,n8?,0) =
*      #include no.rules # no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=1,2,4,5,6,3,7,8,9,level=1,stop=1
*
*--#[ no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=1,2,4,5,6,3,7,8,9,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,1+n2,-1+n3,n4,n5,n6,n7,n8,0)*(rat(-n2,-1+n1
         ))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,1+n2,n3,n4,n5,n6,n7,-1+n8,0)*(rat(n2,-1+n1)
         )
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,0)*(rat(-n7,-1+n1
         ))
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,-1+n5,n6,1+n7,n8,0)*(rat(n7,-1+n1)
         )
      +`Z'(0,0,0,0,x5,0,0,0,0,-1+n1,n2,n3,n4,n5,n6,n7,n8,0)*(rat(-9+2*n8+n7+2*
         n6+2*n5+n2+n1+4*ep+2*ep*x5,-1+n1))
      +`Z'(0,0,0,0,x5,0,0,0,0,n1,n2,n3,n4,n5,-1+n6,n7,n8,0)*(rat(1,1));
*--#] no,noints=5,signs=+,x,+,+,x,+,x,x,0,order=1,2,4,5,6,3,7,8,9,level=1,stop=1 : 
* LHS = Z(1,n2,1,1,n5,1,n7,n8,0)
* DEN = den(-1+n7)*den(-2+n7+n2+ep)*den(-3+n8+n7+n2+2*ep)*den(-6+n8+2*n7+n5+2*n2+4*ep)
* MAX = Z(0,1,0,0,1,0,0,1,0)
* MIN = Z(-1,-1,-1,-1,0,-1,-1,-1,0)
* nterms = 12 (Z: 3, Y: 9)
* NonIntegerIndices: 5
* Signs: 1,+,1,1,o,1,+,+,0
* Order: 7,1,2,3,4,5,6,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,1,n2?pos_,1,1,n5?,1,n7?{>=2},n8?pos_,0) =
*      #include no.rules # no,noints=5,signs=1,+,1,1,o,1,+,+,0,order=7,1,2,3,4,5,6,8,9,level=1,stop=1
*
*--#[ no,noints=5,signs=1,+,1,1,o,1,+,+,0,order=7,1,2,3,4,5,6,8,9,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,0,0,0,0,1,-1+n2,1,1,n5,1,n7,n8,0)*(rat(4-n8-n7-n5-n2-2*ep
         -ep*x5,-3+n8+n7+n2+2*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,1,n2,1,1,1+n5,1,-1+n7,-1+n8,0)*(rat(-24*n5+10*n5*
         n8-n5*n8^2+14*n5*n7-3*n5*n7*n8-2*n5*n7^2+10*n5^2-2*n5^2*n8-3*n5^2*n7-
         n5^3+14*n2*n5-3*n2*n5*n8-4*n2*n5*n7-3*n2*n5^2-2*n2^2*n5-24*ep*x5+10*ep
         *n8*x5-ep*n8^2*x5+14*ep*n7*x5-3*ep*n7*n8*x5-2*ep*n7^2*x5+28*ep*n5+20*
         ep*n5*x5-6*ep*n5*n8-4*ep*n5*n8*x5-8*ep*n5*n7-6*ep*n5*n7*x5-6*ep*n5^2-3
         *ep*n5^2*x5+14*ep*n2*x5-3*ep*n2*n8*x5-4*ep*n2*n7*x5-8*ep*n2*n5-6*ep*n2
         *n5*x5-2*ep*n2^2*x5+28*ep^2*x5+10*ep^2*x5^2-6*ep^2*n8*x5-2*ep^2*n8*
         x5^2-8*ep^2*n7*x5-3*ep^2*n7*x5^2-8*ep^2*n5-12*ep^2*n5*x5-3*ep^2*n5*
         x5^2-8*ep^2*n2*x5-3*ep^2*n2*x5^2-8*ep^3*x5-6*ep^3*x5^2-ep^3*x5^3,-6+2*
         n8+11*n7-3*n7*n8-6*n7^2+n7^2*n8+n7^3+5*n2-n2*n8-7*n2*n7+n2*n7*n8+2*n2*
         n7^2-n2^2+n2^2*n7+7*ep-ep*n8-10*ep*n7+ep*n7*n8+3*ep*n7^2-3*ep*n2+3*ep*
         n2*n7-2*ep^2+2*ep^2*n7))
      +`Z'(0,0,0,0,x5,0,0,0,0,1,n2,1,1,n5,1,-1+n7,n8,0)*(rat(64-52*n8+13*n8^2-
         n8^3-52*n7+30*n7*n8-4*n7*n8^2+13*n7^2-4*n7^2*n8-n7^3-28*n5+16*n5*n8-2*
         n5*n8^2+16*n5*n7-5*n5*n7*n8-2*n5*n7^2+3*n5^2-n5^2*n8-n5^2*n7-60*n2+32*
         n2*n8-4*n2*n8^2+32*n2*n7-9*n2*n7*n8-4*n2*n7^2+18*n2*n5-5*n2*n5*n8-5*n2
         *n5*n7-n2*n5^2+19*n2^2-5*n2^2*n8-5*n2^2*n7-3*n2^2*n5-2*n2^3-100*ep-28*
         ep*x5+55*ep*n8+16*ep*n8*x5-7*ep*n8^2-2*ep*n8^2*x5+55*ep*n7+16*ep*n7*x5
         -16*ep*n7*n8-5*ep*n7*n8*x5-7*ep*n7^2-2*ep*n7^2*x5+27*ep*n5+6*ep*n5*x5-
         8*ep*n5*n8-2*ep*n5*n8*x5-8*ep*n5*n7-2*ep*n5*n7*x5-ep*n5^2+63*ep*n2+18*
         ep*n2*x5-17*ep*n2*n8-5*ep*n2*n8*x5-17*ep*n2*n7-5*ep*n2*n7*x5-9*ep*n2*
         n5-2*ep*n2*n5*x5-10*ep*n2^2-3*ep*n2^2*x5+50*ep^2+27*ep^2*x5+3*ep^2*
         x5^2-14*ep^2*n8-8*ep^2*n8*x5-ep^2*n8*x5^2-14*ep^2*n7-8*ep^2*n7*x5-ep^2
         *n7*x5^2-6*ep^2*n5-2*ep^2*n5*x5-16*ep^2*n2-9*ep^2*n2*x5-ep^2*n2*x5^2-8
         *ep^3-6*ep^3*x5-ep^3*x5^2,-6+2*n8+11*n7-3*n7*n8-6*n7^2+n7^2*n8+n7^3+5*
         n2-n2*n8-7*n2*n7+n2*n7*n8+2*n2*n7^2-n2^2+n2^2*n7+7*ep-ep*n8-10*ep*n7+
         ep*n7*n8+3*ep*n7^2-3*ep*n2+3*ep*n2*n7-2*ep^2+2*ep^2*n7))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,1+n2,0,1,n5,1,n7,n8,0)*(rat(14*n2-9*n2*n8+n2*
         n8^2-6*n2*n7+3*n2*n7*n8-9*n2*n5+2*n2*n5*n8+3*n2*n5*n7+n2*n5^2-4*n2^2+2
         *n2^2*n8+2*n2^2*n5-22*ep*n2-9*ep*n2*x5+6*ep*n2*n8+2*ep*n2*n8*x5+6*ep*
         n2*n7+3*ep*n2*n7*x5+6*ep*n2*n5+2*ep*n2*n5*x5+4*ep*n2^2+2*ep*n2^2*x5+8*
         ep^2*n2+6*ep^2*n2*x5+ep^2*n2*x5^2,-6+2*n8+11*n7-3*n7*n8-6*n7^2+n7^2*n8
         +n7^3+5*n2-n2*n8-7*n2*n7+n2*n7*n8+2*n2*n7^2-n2^2+n2^2*n7+7*ep-ep*n8-10
         *ep*n7+ep*n7*n8+3*ep*n7^2-3*ep*n2+3*ep*n2*n7-2*ep^2+2*ep^2*n7))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,1+n2,1,1,n5,1,n7,-1+n8,0)*(rat(2*n2-n2*n8-n2*n5
         -2*ep*n2-ep*n2*x5,6-2*n8-5*n7+n7*n8+n7^2-5*n2+n2*n8+2*n2*n7+n2^2-7*ep+
         ep*n8+3*ep*n7+3*ep*n2+2*ep^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,n2,1,1,n5,1,n7,n8,0)*(rat(-9+3*n8+2*n7+2*n5+2*
         n2+6*ep+2*ep*x5,-3+n8+n7+n2+2*ep))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,-1+n2,1,1,n5,0,n7,1+n8,0)*(rat(-26*n8+10*n8^2-
         n8^3+17*n7*n8-3*n7*n8^2-3*n7^2*n8+10*n5*n8-2*n5*n8^2-3*n5*n7*n8-n5^2*
         n8+15*n2*n8-3*n2*n8^2-5*n2*n7*n8-3*n2*n5*n8-2*n2^2*n8+35*ep*n8+10*ep*
         n8*x5-7*ep*n8^2-2*ep*n8^2*x5-11*ep*n7*n8-3*ep*n7*n8*x5-7*ep*n5*n8-2*ep
         *n5*n8*x5-10*ep*n2*n8-3*ep*n2*n8*x5-12*ep^2*n8-7*ep^2*n8*x5-ep^2*n8*
         x5^2,-12+2*n8+22*n7-3*n7*n8-12*n7^2+n7^2*n8+2*n7^3+2*n5-3*n5*n7+n5*
         n7^2+10*n2-n2*n8-14*n2*n7+n2*n7*n8+4*n2*n7^2-n2*n5+n2*n5*n7-2*n2^2+2*
         n2^2*n7+14*ep+2*ep*x5-ep*n8-20*ep*n7-3*ep*n7*x5+ep*n7*n8+6*ep*n7^2+ep*
         n7^2*x5-ep*n5+ep*n5*n7-6*ep*n2-ep*n2*x5+6*ep*n2*n7+ep*n2*n7*x5-4*ep^2-
         ep^2*x5+4*ep^2*n7+ep^2*n7*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1+n2,0,1,n5,1,-1+n7,n8,0)*(rat(-12*n2+8*n2*n8-
         n2*n8^2+4*n2*n7-2*n2*n7*n8+8*n2*n5-2*n2*n5*n8-2*n2*n5*n7-n2*n5^2+4*
         n2^2-2*n2^2*n8-2*n2^2*n5+20*ep*n2+8*ep*n2*x5-6*ep*n2*n8-2*ep*n2*n8*x5-
         4*ep*n2*n7-2*ep*n2*n7*x5-6*ep*n2*n5-2*ep*n2*n5*x5-4*ep*n2^2-2*ep*n2^2*
         x5-8*ep^2*n2-6*ep^2*n2*x5-ep^2*n2*x5^2,-6+2*n8+11*n7-3*n7*n8-6*n7^2+
         n7^2*n8+n7^3+5*n2-n2*n8-7*n2*n7+n2*n7*n8+2*n2*n7^2-n2^2+n2^2*n7+7*ep-
         ep*n8-10*ep*n7+ep*n7*n8+3*ep*n7^2-3*ep*n2+3*ep*n2*n7-2*ep^2+2*ep^2*n7)
         )
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,0,1,n5,0,n7,1+n8,0)*(rat(26*n8-10*n8^2+n8^3-
         17*n7*n8+3*n7*n8^2+3*n7^2*n8-10*n5*n8+2*n5*n8^2+3*n5*n7*n8+n5^2*n8-15*
         n2*n8+3*n2*n8^2+5*n2*n7*n8+3*n2*n5*n8+2*n2^2*n8-35*ep*n8-10*ep*n8*x5+7
         *ep*n8^2+2*ep*n8^2*x5+11*ep*n7*n8+3*ep*n7*n8*x5+7*ep*n5*n8+2*ep*n5*n8*
         x5+10*ep*n2*n8+3*ep*n2*n8*x5+12*ep^2*n8+7*ep^2*n8*x5+ep^2*n8*x5^2,-12+
         2*n8+22*n7-3*n7*n8-12*n7^2+n7^2*n8+2*n7^3+2*n5-3*n5*n7+n5*n7^2+10*n2-
         n2*n8-14*n2*n7+n2*n7*n8+4*n2*n7^2-n2*n5+n2*n5*n7-2*n2^2+2*n2^2*n7+14*
         ep+2*ep*x5-ep*n8-20*ep*n7-3*ep*n7*x5+ep*n7*n8+6*ep*n7^2+ep*n7^2*x5-ep*
         n5+ep*n5*n7-6*ep*n2-ep*n2*x5+6*ep*n2*n7+ep*n2*n7*x5-4*ep^2-ep^2*x5+4*
         ep^2*n7+ep^2*n7*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,0,1,n5,1,n7,n8,0)*(rat(60-52*n8+13*n8^2-n8^3
         -44*n7+30*n7*n8-4*n7*n8^2+8*n7^2-4*n7^2*n8-52*n5+26*n5*n8-3*n5*n8^2+30
         *n5*n7-8*n5*n7*n8-4*n5*n7^2+13*n5^2-3*n5^2*n8-4*n5^2*n7-n5^3-32*n2+22*
         n2*n8-3*n2*n8^2+12*n2*n7-6*n2*n7*n8+22*n2*n5-6*n2*n5*n8-6*n2*n5*n7-3*
         n2*n5^2+4*n2^2-2*n2^2*n8-2*n2^2*n5-136*ep-52*ep*x5+74*ep*n8+26*ep*n8*
         x5-9*ep*n8^2-3*ep*n8^2*x5+72*ep*n7+30*ep*n7*x5-22*ep*n7*n8-8*ep*n7*n8*
         x5-8*ep*n7^2-4*ep*n7^2*x5+74*ep*n5+26*ep*n5*x5-18*ep*n5*n8-6*ep*n5*n8*
         x5-22*ep*n5*n7-8*ep*n5*n7*x5-9*ep*n5^2-3*ep*n5^2*x5+52*ep*n2+22*ep*n2*
         x5-16*ep*n2*n8-6*ep*n2*n8*x5-12*ep*n2*n7-6*ep*n2*n7*x5-16*ep*n2*n5-6*
         ep*n2*n5*x5-4*ep*n2^2-2*ep*n2^2*x5+100*ep^2+74*ep^2*x5+13*ep^2*x5^2-26
         *ep^2*n8-18*ep^2*n8*x5-3*ep^2*n8*x5^2-28*ep^2*n7-22*ep^2*n7*x5-4*ep^2*
         n7*x5^2-26*ep^2*n5-18*ep^2*n5*x5-3*ep^2*n5*x5^2-20*ep^2*n2-16*ep^2*n2*
         x5-3*ep^2*n2*x5^2-24*ep^3-26*ep^3*x5-9*ep^3*x5^2-ep^3*x5^3,-6+2*n8+11*
         n7-3*n7*n8-6*n7^2+n7^2*n8+n7^3+5*n2-n2*n8-7*n2*n7+n2*n7*n8+2*n2*n7^2-
         n2^2+n2^2*n7+7*ep-ep*n8-10*ep*n7+ep*n7*n8+3*ep*n7^2-3*ep*n2+3*ep*n2*n7
         -2*ep^2+2*ep^2*n7))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,0,1+n5,0,n7,n8,0)*(rat(5*n5-n5*n8-n5*n7-
         n5^2-2*n2*n5+5*ep*x5-ep*n8*x5-ep*n7*x5-4*ep*n5-2*ep*n5*x5-2*ep*n2*x5-4
         *ep^2*x5-ep^2*x5^2,6-n8-8*n7+n7*n8+2*n7^2-n5+n5*n7-2*n2+2*n2*n7-4*ep-
         ep*x5+4*ep*n7+ep*n7*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,1,1+n5,0,-1+n7,n8,0)*(rat(-174*n5+100*n5*
         n8-18*n5*n8^2+n5*n8^3+163*n5*n7-64*n5*n7*n8+6*n5*n7*n8^2-50*n5*n7^2+10
         *n5*n7^2*n8+5*n5*n7^3+90*n5^2-34*n5^2*n8+3*n5^2*n8^2-57*n5^2*n7+11*
         n5^2*n7*n8+9*n5^2*n7^2-16*n5^3+3*n5^3*n8+5*n5^3*n7+n5^4+169*n2*n5-66*
         n2*n5*n8+6*n2*n5*n8^2-105*n2*n5*n7+21*n2*n5*n7*n8+16*n2*n5*n7^2-57*n2*
         n5^2+11*n2*n5^2*n8+18*n2*n5^2*n7+5*n2*n5^3-55*n2^2*n5+11*n2^2*n5*n8+17
         *n2^2*n5*n7+9*n2^2*n5^2+6*n2^3*n5-174*ep*x5+100*ep*n8*x5-18*ep*n8^2*x5
         +ep*n8^3*x5+163*ep*n7*x5-64*ep*n7*n8*x5+6*ep*n7*n8^2*x5-50*ep*n7^2*x5+
         10*ep*n7^2*n8*x5+5*ep*n7^3*x5+323*ep*n5+180*ep*n5*x5-124*ep*n5*n8-68*
         ep*n5*n8*x5+11*ep*n5*n8^2+6*ep*n5*n8^2*x5-202*ep*n5*n7-114*ep*n5*n7*x5
         +40*ep*n5*n7*n8+22*ep*n5*n7*n8*x5+31*ep*n5*n7^2+18*ep*n5*n7^2*x5-111*
         ep*n5^2-48*ep*n5^2*x5+21*ep*n5^2*n8+9*ep*n5^2*n8*x5+35*ep*n5^2*n7+15*
         ep*n5^2*n7*x5+10*ep*n5^3+4*ep*n5^3*x5+169*ep*n2*x5-66*ep*n2*n8*x5+6*ep
         *n2*n8^2*x5-105*ep*n2*n7*x5+21*ep*n2*n7*n8*x5+16*ep*n2*n7^2*x5-209*ep*
         n2*n5-114*ep*n2*n5*x5+41*ep*n2*n5*n8+22*ep*n2*n5*n8*x5+65*ep*n2*n5*n7+
         36*ep*n2*n5*n7*x5+35*ep*n2*n5^2+15*ep*n2*n5^2*x5-55*ep*n2^2*x5+11*ep*
         n2^2*n8*x5+17*ep*n2^2*n7*x5+34*ep*n2^2*n5+18*ep*n2^2*n5*x5+6*ep*n2^3*
         x5+323*ep^2*x5+90*ep^2*x5^2-124*ep^2*n8*x5-34*ep^2*n8*x5^2+11*ep^2*
         n8^2*x5+3*ep^2*n8^2*x5^2-202*ep^2*n7*x5-57*ep^2*n7*x5^2+40*ep^2*n7*n8*
         x5+11*ep^2*n7*n8*x5^2+31*ep^2*n7^2*x5+9*ep^2*n7^2*x5^2-198*ep^2*n5-222
         *ep^2*n5*x5-48*ep^2*n5*x5^2+38*ep^2*n5*n8+42*ep^2*n5*n8*x5+9*ep^2*n5*
         n8*x5^2+62*ep^2*n5*n7+70*ep^2*n5*n7*x5+15*ep^2*n5*n7*x5^2+34*ep^2*n5^2
         +30*ep^2*n5^2*x5+6*ep^2*n5^2*x5^2-209*ep^2*n2*x5-57*ep^2*n2*x5^2+41*
         ep^2*n2*n8*x5+11*ep^2*n2*n8*x5^2+65*ep^2*n2*n7*x5+18*ep^2*n2*n7*x5^2+
         64*ep^2*n2*n5+70*ep^2*n2*n5*x5+15*ep^2*n2*n5*x5^2+34*ep^2*n2^2*x5+9*
         ep^2*n2^2*x5^2-198*ep^3*x5-111*ep^3*x5^2-16*ep^3*x5^3+38*ep^3*n8*x5+21
         *ep^3*n8*x5^2+3*ep^3*n8*x5^3+62*ep^3*n7*x5+35*ep^3*n7*x5^2+5*ep^3*n7*
         x5^3+40*ep^3*n5+68*ep^3*n5*x5+30*ep^3*n5*x5^2+4*ep^3*n5*x5^3+64*ep^3*
         n2*x5+35*ep^3*n2*x5^2+5*ep^3*n2*x5^3+40*ep^4*x5+34*ep^4*x5^2+10*ep^4*
         x5^3+ep^4*x5^4,36-18*n8+2*n8^2-78*n7+33*n7*n8-3*n7*n8^2+58*n7^2-18*
         n7^2*n8+n7^2*n8^2-18*n7^3+3*n7^3*n8+2*n7^4-6*n5+2*n5*n8+11*n5*n7-3*n5*
         n7*n8-6*n5*n7^2+n5*n7^2*n8+n5*n7^3-42*n2+15*n2*n8-n2*n8^2+74*n2*n7-21*
         n2*n7*n8+n2*n7*n8^2-38*n2*n7^2+6*n2*n7^2*n8+6*n2*n7^3+5*n2*n5-n2*n5*n8
         -7*n2*n5*n7+n2*n5*n7*n8+2*n2*n5*n7^2+16*n2^2-3*n2^2*n8-22*n2^2*n7+3*
         n2^2*n7*n8+6*n2^2*n7^2-n2^2*n5+n2^2*n5*n7-2*n2^3+2*n2^3*n7-66*ep-6*ep*
         x5+21*ep*n8+2*ep*n8*x5-ep*n8^2+118*ep*n7+11*ep*n7*x5-30*ep*n7*n8-3*ep*
         n7*n8*x5+ep*n7*n8^2-62*ep*n7^2-6*ep*n7^2*x5+9*ep*n7^2*n8+ep*n7^2*n8*x5
         +10*ep*n7^3+ep*n7^3*x5+7*ep*n5-ep*n5*n8-10*ep*n5*n7+ep*n5*n7*n8+3*ep*
         n5*n7^2+52*ep*n2+5*ep*n2*x5-9*ep*n2*n8-ep*n2*n8*x5-72*ep*n2*n7-7*ep*n2
         *n7*x5+9*ep*n2*n7*n8+ep*n2*n7*n8*x5+20*ep*n2*n7^2+2*ep*n2*n7^2*x5-3*ep
         *n2*n5+3*ep*n2*n5*n7-10*ep*n2^2-ep*n2^2*x5+10*ep*n2^2*n7+ep*n2^2*n7*x5
         +40*ep^2+7*ep^2*x5-6*ep^2*n8-ep^2*n8*x5-56*ep^2*n7-10*ep^2*n7*x5+6*
         ep^2*n7*n8+ep^2*n7*n8*x5+16*ep^2*n7^2+3*ep^2*n7^2*x5-2*ep^2*n5+2*ep^2*
         n5*n7-16*ep^2*n2-3*ep^2*n2*x5+16*ep^2*n2*n7+3*ep^2*n2*n7*x5-8*ep^3-2*
         ep^3*x5+8*ep^3*n7+2*ep^3*n7*x5));
*--#] no,noints=5,signs=1,+,1,1,o,1,+,+,0,order=7,1,2,3,4,5,6,8,9,level=1,stop=1 : 
* LHS = Z(1,n2,1,1,n5,1,1,n8,0)
* DEN = den(-1+n8)*den(-2+n8+n2+2*ep)*den(-3+n8+n5+n2+3*ep)*den(-4+2*n8+n5+n2+2*ep)
* MAX = Z(0,1,0,0,1,0,0,0,0)
* MIN = Z(-1,-1,-1,-1,0,-1,-1,-1,0)
* nterms = 12 (Z: 2, Y: 10)
* NonIntegerIndices: 5
* Signs: 1,+,1,1,o,1,1,+,0
* Order: 8,1,2,3,4,5,6,7,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,1,n2?pos_,1,1,n5?,1,1,n8?{>=2},0) =
*      #include no.rules # no,noints=5,signs=1,+,1,1,o,1,1,+,0,order=8,1,2,3,4,5,6,7,9,level=1,stop=1
*
*--#[ no,noints=5,signs=1,+,1,1,o,1,1,+,0,order=8,1,2,3,4,5,6,7,9,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,0,0,0,0,1,-1+n2,1,1,n5,1,1,n8,0)*(rat(3-n8-n5-n2-2*ep-ep*
         x5,-2+n8+n2+2*ep))
      +`Z'(0,0,0,0,x5,0,0,0,0,1,n2,1,1,n5,1,1,-1+n8,0)*(rat(-12+7*n8-n8^2+7*n5-
         2*n5*n8-n5^2+10*n2-3*n2*n8-3*n2*n5-2*n2^2+20*ep+7*ep*x5-6*ep*n8-2*ep*
         n8*x5-6*ep*n5-2*ep*n5*x5-8*ep*n2-3*ep*n2*x5-8*ep^2-6*ep^2*x5-ep^2*x5^2
         ,2-3*n8+n8^2-n2+n2*n8-2*ep+2*ep*n8))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,1+n2,0,1,n5,1,1,n8,0)*(rat(8*n2-3*n2*n8-3*n2*n5
         -2*n2^2-8*ep*n2-3*ep*n2*x5,6-5*n8+n8^2-2*n5+n5*n8-5*n2+2*n2*n8+n2*n5+
         n2^2-12*ep-2*ep*x5+5*ep*n8+ep*n8*x5+2*ep*n5+5*ep*n2+ep*n2*x5+6*ep^2+2*
         ep^2*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,1+n2,1,1,n5,1,1,-1+n8,0)*(rat(-6*n2+3*n2*n8+n2*
         n5+2*n2^2+4*ep*n2+ep*n2*x5,2-3*n8+n8^2-n2+n2*n8-2*ep+2*ep*n8))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,n2,1,1,n5,1,1,n8,0)*(rat(-7+3*n8+2*n5+2*n2+6*ep
         +2*ep*x5,-2+n8+n2+2*ep))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1+n2,0,1,n5,1,0,n8,0)*(rat(11*n2-12*n2*n8+3*n2*
         n8^2-7*n2*n5+4*n2*n5*n8+n2*n5^2-3*n2^2+2*n2^2*n8+n2^2*n5-10*ep*n2-7*ep
         *n2*x5+7*ep*n2*n8+4*ep*n2*n8*x5+3*ep*n2*n5+2*ep*n2*n5*x5+ep*n2^2*x5+3*
         ep^2*n2*x5+ep^2*n2*x5^2,-24+32*n8-14*n8^2+2*n8^3+14*n5-13*n5*n8+3*n5*
         n8^2-2*n5^2+n5^2*n8+26*n2-23*n2*n8+5*n2*n8^2-11*n2*n5+5*n2*n5*n8+n2*
         n5^2-9*n2^2+4*n2^2*n8+2*n2^2*n5+n2^3+60*ep+14*ep*x5-54*ep*n8-13*ep*n8*
         x5+12*ep*n8^2+3*ep*n8^2*x5-24*ep*n5-4*ep*n5*x5+11*ep*n5*n8+2*ep*n5*n8*
         x5+2*ep*n5^2-42*ep*n2-11*ep*n2*x5+19*ep*n2*n8+5*ep*n2*n8*x5+9*ep*n2*n5
         +2*ep*n2*n5*x5+7*ep*n2^2+2*ep*n2^2*x5-48*ep^2-24*ep^2*x5-2*ep^2*x5^2+
         22*ep^2*n8+11*ep^2*n8*x5+ep^2*n8*x5^2+10*ep^2*n5+4*ep^2*n5*x5+16*ep^2*
         n2+9*ep^2*n2*x5+ep^2*n2*x5^2+12*ep^3+10*ep^3*x5+2*ep^3*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1+n2,1,1,n5,1,0,-1+n8,0)*(rat(51*n2-59*n2*n8+23
         *n2*n8^2-3*n2*n8^3-41*n2*n5+31*n2*n5*n8-6*n2*n5*n8^2+11*n2*n5^2-4*n2*
         n5^2*n8-n2*n5^3-53*n2^2+41*n2^2*n8-8*n2^2*n8^2+29*n2^2*n5-11*n2^2*n5*
         n8-4*n2^2*n5^2+18*n2^3-7*n2^3*n8-5*n2^3*n5-2*n2^4-118*ep*n2-41*ep*n2*
         x5+89*ep*n2*n8+31*ep*n2*n8*x5-17*ep*n2*n8^2-6*ep*n2*n8^2*x5+65*ep*n2*
         n5+22*ep*n2*n5*x5-24*ep*n2*n5*n8-8*ep*n2*n5*n8*x5-9*ep*n2*n5^2-3*ep*n2
         *n5^2*x5+82*ep*n2^2+29*ep*n2^2*x5-31*ep*n2^2*n8-11*ep*n2^2*n8*x5-23*ep
         *n2^2*n5-8*ep*n2^2*n5*x5-14*ep*n2^3-5*ep*n2^3*x5+92*ep^2*n2+65*ep^2*n2
         *x5+11*ep^2*n2*x5^2-34*ep^2*n2*n8-24*ep^2*n2*n8*x5-4*ep^2*n2*n8*x5^2-
         26*ep^2*n2*n5-18*ep^2*n2*n5*x5-3*ep^2*n2*n5*x5^2-32*ep^2*n2^2-23*ep^2*
         n2^2*x5-4*ep^2*n2^2*x5^2-24*ep^3*n2-26*ep^3*n2*x5-9*ep^3*n2*x5^2-ep^3*
         n2*x5^3,24-56*n8+46*n8^2-16*n8^3+2*n8^4-14*n5+27*n5*n8-16*n5*n8^2+3*n5
         *n8^3+2*n5^2-3*n5^2*n8+n5^2*n8^2-26*n2+49*n2*n8-28*n2*n8^2+5*n2*n8^3+
         11*n2*n5-16*n2*n5*n8+5*n2*n5*n8^2-n2*n5^2+n2*n5^2*n8+9*n2^2-13*n2^2*n8
         +4*n2^2*n8^2-2*n2^2*n5+2*n2^2*n5*n8-n2^3+n2^3*n8-60*ep-14*ep*x5+114*ep
         *n8+27*ep*n8*x5-66*ep*n8^2-16*ep*n8^2*x5+12*ep*n8^3+3*ep*n8^3*x5+24*ep
         *n5+4*ep*n5*x5-35*ep*n5*n8-6*ep*n5*n8*x5+11*ep*n5*n8^2+2*ep*n5*n8^2*x5
         -2*ep*n5^2+2*ep*n5^2*n8+42*ep*n2+11*ep*n2*x5-61*ep*n2*n8-16*ep*n2*n8*
         x5+19*ep*n2*n8^2+5*ep*n2*n8^2*x5-9*ep*n2*n5-2*ep*n2*n5*x5+9*ep*n2*n5*
         n8+2*ep*n2*n5*n8*x5-7*ep*n2^2-2*ep*n2^2*x5+7*ep*n2^2*n8+2*ep*n2^2*n8*
         x5+48*ep^2+24*ep^2*x5+2*ep^2*x5^2-70*ep^2*n8-35*ep^2*n8*x5-3*ep^2*n8*
         x5^2+22*ep^2*n8^2+11*ep^2*n8^2*x5+ep^2*n8^2*x5^2-10*ep^2*n5-4*ep^2*n5*
         x5+10*ep^2*n5*n8+4*ep^2*n5*n8*x5-16*ep^2*n2-9*ep^2*n2*x5-ep^2*n2*x5^2+
         16*ep^2*n2*n8+9*ep^2*n2*n8*x5+ep^2*n2*n8*x5^2-12*ep^3-10*ep^3*x5-2*
         ep^3*x5^2+12*ep^3*n8+10*ep^3*n8*x5+2*ep^3*n8*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,0,1,n5,1,1,n8,0)*(rat(-6+2*n8+2*n5+2*n2+6*ep
         +2*ep*x5,-2+n8+n2+2*ep))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,0,1+n5,0,1,n8,0)*(rat(n5+ep*x5,-3+n8+n5+n2
         +3*ep+ep*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,0,1+n5,1,1,-1+n8,0)*(rat(-4*n5+n5*n8+n5^2+
         2*n2*n5-4*ep*x5+ep*n8*x5+4*ep*n5+2*ep*n5*x5+2*ep*n2*x5+4*ep^2*x5+ep^2*
         x5^2,2-3*n8+n8^2-n2+n2*n8-2*ep+2*ep*n8))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,1,1+n5,0,0,n8,0)*(rat(-17*n5+14*n5*n8-3*n5
         *n8^2+8*n5^2-3*n5^2*n8-n5^3+12*n2*n5-5*n2*n5*n8-3*n2*n5^2-2*n2^2*n5-17
         *ep*x5+14*ep*n8*x5-3*ep*n8^2*x5+28*ep*n5+16*ep*n5*x5-11*ep*n5*n8-6*ep*
         n5*n8*x5-7*ep*n5^2-3*ep*n5^2*x5+12*ep*n2*x5-5*ep*n2*n8*x5-10*ep*n2*n5-
         6*ep*n2*n5*x5-2*ep*n2^2*x5+28*ep^2*x5+8*ep^2*x5^2-11*ep^2*n8*x5-3*ep^2
         *n8*x5^2-12*ep^2*n5-14*ep^2*n5*x5-3*ep^2*n5*x5^2-10*ep^2*n2*x5-3*ep^2*
         n2*x5^2-12*ep^3*x5-7*ep^3*x5^2-ep^3*x5^3,-24+32*n8-14*n8^2+2*n8^3+14*
         n5-13*n5*n8+3*n5*n8^2-2*n5^2+n5^2*n8+26*n2-23*n2*n8+5*n2*n8^2-11*n2*n5
         +5*n2*n5*n8+n2*n5^2-9*n2^2+4*n2^2*n8+2*n2^2*n5+n2^3+60*ep+14*ep*x5-54*
         ep*n8-13*ep*n8*x5+12*ep*n8^2+3*ep*n8^2*x5-24*ep*n5-4*ep*n5*x5+11*ep*n5
         *n8+2*ep*n5*n8*x5+2*ep*n5^2-42*ep*n2-11*ep*n2*x5+19*ep*n2*n8+5*ep*n2*
         n8*x5+9*ep*n2*n5+2*ep*n2*n5*x5+7*ep*n2^2+2*ep*n2^2*x5-48*ep^2-24*ep^2*
         x5-2*ep^2*x5^2+22*ep^2*n8+11*ep^2*n8*x5+ep^2*n8*x5^2+10*ep^2*n5+4*ep^2
         *n5*x5+16*ep^2*n2+9*ep^2*n2*x5+ep^2*n2*x5^2+12*ep^3+10*ep^3*x5+2*ep^3*
         x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,1,1+n5,1,0,-1+n8,0)*(rat(39*n5-37*n5*n8+11
         *n5*n8^2-n5*n8^3-34*n5^2+21*n5^2*n8-3*n5^2*n8^2+10*n5^3-3*n5^3*n8-n5^4
         -46*n2*n5+31*n2*n5*n8-5*n2*n5*n8^2+27*n2*n5^2-9*n2*n5^2*n8-4*n2*n5^3+
         17*n2^2*n5-6*n2^2*n5*n8-5*n2^2*n5^2-2*n2^3*n5+39*ep*x5-37*ep*n8*x5+11*
         ep*n8^2*x5-ep*n8^3*x5-104*ep*n5-68*ep*n5*x5+69*ep*n5*n8+42*ep*n5*n8*x5
         -11*ep*n5*n8^2-6*ep*n5*n8^2*x5+61*ep*n5^2+30*ep*n5^2*x5-20*ep*n5^2*n8-
         9*ep*n5^2*n8*x5-9*ep*n5^3-4*ep*n5^3*x5-46*ep*n2*x5+31*ep*n2*n8*x5-5*ep
         *n2*n8^2*x5+78*ep*n2*n5+54*ep*n2*n5*x5-27*ep*n2*n5*n8-18*ep*n2*n5*n8*
         x5-23*ep*n2*n5^2-12*ep*n2*n5^2*x5+17*ep*n2^2*x5-6*ep*n2^2*n8*x5-14*ep*
         n2^2*n5-10*ep*n2^2*n5*x5-2*ep*n2^3*x5-104*ep^2*x5-34*ep^2*x5^2+69*ep^2
         *n8*x5+21*ep^2*n8*x5^2-11*ep^2*n8^2*x5-3*ep^2*n8^2*x5^2+88*ep^2*n5+122
         *ep^2*n5*x5+30*ep^2*n5*x5^2-30*ep^2*n5*n8-40*ep^2*n5*n8*x5-9*ep^2*n5*
         n8*x5^2-26*ep^2*n5^2-27*ep^2*n5^2*x5-6*ep^2*n5^2*x5^2+78*ep^2*n2*x5+27
         *ep^2*n2*x5^2-27*ep^2*n2*n8*x5-9*ep^2*n2*n8*x5^2-32*ep^2*n2*n5-46*ep^2
         *n2*n5*x5-12*ep^2*n2*n5*x5^2-14*ep^2*n2^2*x5-5*ep^2*n2^2*x5^2+88*ep^3*
         x5+61*ep^3*x5^2+10*ep^3*x5^3-30*ep^3*n8*x5-20*ep^3*n8*x5^2-3*ep^3*n8*
         x5^3-24*ep^3*n5-52*ep^3*n5*x5-27*ep^3*n5*x5^2-4*ep^3*n5*x5^3-32*ep^3*
         n2*x5-23*ep^3*n2*x5^2-4*ep^3*n2*x5^3-24*ep^4*x5-26*ep^4*x5^2-9*ep^4*
         x5^3-ep^4*x5^4,24-56*n8+46*n8^2-16*n8^3+2*n8^4-14*n5+27*n5*n8-16*n5*
         n8^2+3*n5*n8^3+2*n5^2-3*n5^2*n8+n5^2*n8^2-26*n2+49*n2*n8-28*n2*n8^2+5*
         n2*n8^3+11*n2*n5-16*n2*n5*n8+5*n2*n5*n8^2-n2*n5^2+n2*n5^2*n8+9*n2^2-13
         *n2^2*n8+4*n2^2*n8^2-2*n2^2*n5+2*n2^2*n5*n8-n2^3+n2^3*n8-60*ep-14*ep*
         x5+114*ep*n8+27*ep*n8*x5-66*ep*n8^2-16*ep*n8^2*x5+12*ep*n8^3+3*ep*n8^3
         *x5+24*ep*n5+4*ep*n5*x5-35*ep*n5*n8-6*ep*n5*n8*x5+11*ep*n5*n8^2+2*ep*
         n5*n8^2*x5-2*ep*n5^2+2*ep*n5^2*n8+42*ep*n2+11*ep*n2*x5-61*ep*n2*n8-16*
         ep*n2*n8*x5+19*ep*n2*n8^2+5*ep*n2*n8^2*x5-9*ep*n2*n5-2*ep*n2*n5*x5+9*
         ep*n2*n5*n8+2*ep*n2*n5*n8*x5-7*ep*n2^2-2*ep*n2^2*x5+7*ep*n2^2*n8+2*ep*
         n2^2*n8*x5+48*ep^2+24*ep^2*x5+2*ep^2*x5^2-70*ep^2*n8-35*ep^2*n8*x5-3*
         ep^2*n8*x5^2+22*ep^2*n8^2+11*ep^2*n8^2*x5+ep^2*n8^2*x5^2-10*ep^2*n5-4*
         ep^2*n5*x5+10*ep^2*n5*n8+4*ep^2*n5*n8*x5-16*ep^2*n2-9*ep^2*n2*x5-ep^2*
         n2*x5^2+16*ep^2*n2*n8+9*ep^2*n2*n8*x5+ep^2*n2*n8*x5^2-12*ep^3-10*ep^3*
         x5-2*ep^3*x5^2+12*ep^3*n8+10*ep^3*n8*x5+2*ep^3*n8*x5^2));
*--#] no,noints=5,signs=1,+,1,1,o,1,1,+,0,order=8,1,2,3,4,5,6,7,9,level=1,stop=1 : 
* LHS = Z(1,n2,1,1,n5,1,1,1,0)
* DEN = den(-1+n2)*den(-1+n2+ep)*den(-1+n2+2*ep)*den(-1+n5+2*ep)*den(-2+n2+3*ep)*den(-2+n5+n2+2*ep)*den(-3+n5+n2+3*ep)*den(-4+n5+2*n2+4*ep)*den(-4+2*n5+n2+4*ep)
* MAX = Z(1,1,0,0,1,0,0,0,0)
* MIN = Z(-1,-1,-1,-1,-1,-1,-1,-1,0)
* nterms = 23 (Z: 2, Y: 21)
* NonIntegerIndices: 5
* Signs: 1,+,1,1,x,1,1,1,0
* Order: 9,8,7,6,5,4,3,2,1
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 2
* Stop: 2
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,1,n2?{>=2},1,1,n5?,1,1,1,0) =
*      #include no.rules # no,noints=5,signs=1,+,1,1,x,1,1,1,0,order=9,8,7,6,5,4,3,2,1,level=2,stop=2
*
*--#[ no,noints=5,signs=1,+,1,1,x,1,1,1,0,order=9,8,7,6,5,4,3,2,1,level=2,stop=2 :
      +`Z'(0,0,0,0,x5,0,0,0,0,1,-1+n2,1,1,-1+n5,1,1,1,0)*(rat(81-108*n5+54*n5^2
         -12*n5^3+n5^4-135*n2+135*n2*n5-45*n2*n5^2+5*n2*n5^3+81*n2^2-54*n2^2*n5
         +9*n2^2*n5^2-21*n2^3+7*n2^3*n5+2*n2^4-351*ep-108*ep*x5+351*ep*n5+108*
         ep*n5*x5-117*ep*n5^2-36*ep*n5^2*x5+13*ep*n5^3+4*ep*n5^3*x5+432*ep*n2+
         135*ep*n2*x5-288*ep*n2*n5-90*ep*n2*n5*x5+48*ep*n2*n5^2+15*ep*n2*n5^2*
         x5-171*ep*n2^2-54*ep*n2^2*x5+57*ep*n2^2*n5+18*ep*n2^2*n5*x5+22*ep*n2^3
         +7*ep*n2^3*x5+558*ep^2+351*ep^2*x5+54*ep^2*x5^2-372*ep^2*n5-234*ep^2*
         n5*x5-36*ep^2*n5*x5^2+62*ep^2*n5^2+39*ep^2*n5^2*x5+6*ep^2*n5^2*x5^2-
         450*ep^2*n2-288*ep^2*n2*x5-45*ep^2*n2*x5^2+150*ep^2*n2*n5+96*ep^2*n2*
         n5*x5+15*ep^2*n2*n5*x5^2+88*ep^2*n2^2+57*ep^2*n2^2*x5+9*ep^2*n2^2*x5^2
         -384*ep^3-372*ep^3*x5-117*ep^3*x5^2-12*ep^3*x5^3+128*ep^3*n5+124*ep^3*
         n5*x5+39*ep^3*n5*x5^2+4*ep^3*n5*x5^3+152*ep^3*n2+150*ep^3*n2*x5+48*
         ep^3*n2*x5^2+5*ep^3*n2*x5^3+96*ep^4+128*ep^4*x5+62*ep^4*x5^2+13*ep^4*
         x5^3+ep^4*x5^4,4-2*n5-13*n2+6*n2*n5+15*n2^2-6*n2^2*n5-7*n2^3+2*n2^3*n5
         +n2^4-16*ep-2*ep*x5+6*ep*n5+39*ep*n2+6*ep*n2*x5-12*ep*n2*n5-30*ep*n2^2
         -6*ep*n2^2*x5+6*ep*n2^2*n5+7*ep*n2^3+2*ep*n2^3*x5+20*ep^2+6*ep^2*x5-4*
         ep^2*n5-34*ep^2*n2-12*ep^2*n2*x5+4*ep^2*n2*n5+14*ep^2*n2^2+6*ep^2*n2^2
         *x5-8*ep^3-4*ep^3*x5+8*ep^3*n2+4*ep^3*n2*x5))
      +`Z'(0,0,0,0,x5,0,0,0,0,1,-1+n2,1,1,n5,1,1,1,0)*(rat(-29+54*n5-36*n5^2+10
         *n5^3-n5^4+54*n2-78*n2*n5+36*n2*n5^2-5*n2*n5^3-36*n2^2+36*n2^2*n5-9*
         n2^2*n5^2+10*n2^3-5*n2^3*n5-n2^4+127*ep+54*ep*x5-177*ep*n5-72*ep*n5*x5
         +79*ep*n5^2+30*ep*n5^2*x5-11*ep*n5^3-4*ep*n5^3*x5-177*ep*n2-78*ep*n2*
         x5+169*ep*n2*n5+72*ep*n2*n5*x5-39*ep*n2*n5^2-15*ep*n2*n5^2*x5+79*ep*
         n2^2+36*ep*n2^2*x5-39*ep*n2^2*n5-18*ep*n2^2*n5*x5-11*ep*n2^3-5*ep*n2^3
         *x5-202*ep^2-177*ep^2*x5-36*ep^2*x5^2+188*ep^2*n5+158*ep^2*n5*x5+30*
         ep^2*n5*x5^2-42*ep^2*n5^2-33*ep^2*n5^2*x5-6*ep^2*n5^2*x5^2+188*ep^2*n2
         +169*ep^2*n2*x5+36*ep^2*n2*x5^2-90*ep^2*n2*n5-78*ep^2*n2*n5*x5-15*ep^2
         *n2*n5*x5^2-42*ep^2*n2^2-39*ep^2*n2^2*x5-9*ep^2*n2^2*x5^2+136*ep^3+188
         *ep^3*x5+79*ep^3*x5^2+10*ep^3*x5^3-64*ep^3*n5-84*ep^3*n5*x5-33*ep^3*n5
         *x5^2-4*ep^3*n5*x5^3-64*ep^3*n2-90*ep^3*n2*x5-39*ep^3*n2*x5^2-5*ep^3*
         n2*x5^3-32*ep^4-64*ep^4*x5-42*ep^4*x5^2-11*ep^4*x5^3-ep^4*x5^4,4-2*n5-
         13*n2+6*n2*n5+15*n2^2-6*n2^2*n5-7*n2^3+2*n2^3*n5+n2^4-16*ep-2*ep*x5+6*
         ep*n5+39*ep*n2+6*ep*n2*x5-12*ep*n2*n5-30*ep*n2^2-6*ep*n2^2*x5+6*ep*
         n2^2*n5+7*ep*n2^3+2*ep*n2^3*x5+20*ep^2+6*ep^2*x5-4*ep^2*n5-34*ep^2*n2-
         12*ep^2*n2*x5+4*ep^2*n2*n5+14*ep^2*n2^2+6*ep^2*n2^2*x5-8*ep^3-4*ep^3*
         x5+8*ep^3*n2+4*ep^3*n2*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,1+n2,1,0,n5,1,1,0,0)*(rat(28*n2-30*n2*n5+10*n2*
         n5^2-n2*n5^3-30*n2^2+23*n2^2*n5-4*n2^2*n5^2+8*n2^3-4*n2^3*n5-98*ep*n2-
         30*ep*n2*x5+71*ep*n2*n5+20*ep*n2*n5*x5-12*ep*n2*n5^2-3*ep*n2*n5^2*x5+
         68*ep*n2^2+23*ep*n2^2*x5-27*ep*n2^2*n5-8*ep*n2^2*n5*x5-8*ep*n2^3-4*ep*
         n2^3*x5+114*ep^2*n2+71*ep^2*n2*x5+10*ep^2*n2*x5^2-42*ep^2*n2*n5-24*
         ep^2*n2*n5*x5-3*ep^2*n2*n5*x5^2-38*ep^2*n2^2-27*ep^2*n2^2*x5-4*ep^2*
         n2^2*x5^2-44*ep^3*n2-42*ep^3*n2*x5-12*ep^3*n2*x5^2-ep^3*n2*x5^3,4-2*n5
         -12*n2+5*n2*n5+13*n2^2-4*n2^2*n5-6*n2^3+n2^3*n5+n2^4-22*ep-2*ep*x5+9*
         ep*n5+49*ep*n2+5*ep*n2*x5-15*ep*n2*n5-35*ep*n2^2-4*ep*n2^2*x5+6*ep*
         n2^2*n5+8*ep*n2^3+ep*n2^3*x5+44*ep^2+9*ep^2*x5-13*ep^2*n5-65*ep^2*n2-
         15*ep^2*n2*x5+11*ep^2*n2*n5+23*ep^2*n2^2+6*ep^2*n2^2*x5-38*ep^3-13*
         ep^3*x5+6*ep^3*n5+28*ep^3*n2+11*ep^3*n2*x5+12*ep^4+6*ep^4*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,n2,0,1,n5,1,1,1,0)*(rat(-225+456*n5-364*n5^2+
         142*n5^3-27*n5^4+2*n5^5+468*n2-768*n2*n5+468*n2*n5^2-124*n2*n5^3+12*n2
         *n5^4-374*n2^2+464*n2^2*n5-192*n2^2*n5^2+26*n2^2*n5^3+144*n2^3-119*
         n2^3*n5+25*n2^3*n5^2-27*n2^4+11*n2^4*n5+2*n2^5+1392*ep+456*ep*x5-2255*
         ep*n5-728*ep*n5*x5+1347*ep*n5^2+426*ep*n5^2*x5-349*ep*n5^3-108*ep*n5^3
         *x5+33*ep*n5^4+10*ep*n5^4*x5-2298*ep*n2-768*ep*n2*x5+2832*ep*n2*n5+936
         *ep*n2*n5*x5-1150*ep*n2*n5^2-372*ep*n2*n5^2*x5+152*ep*n2*n5^3+48*ep*n2
         *n5^3*x5+1364*ep*n2^2+464*ep*n2^2*x5-1133*ep*n2^2*n5-384*ep*n2^2*n5*x5
         +235*ep*n2^2*n5^2+78*ep*n2^2*n5^2*x5-346*ep*n2^3-119*ep*n2^3*x5+144*ep
         *n2^3*n5+50*ep*n2^3*n5*x5+32*ep*n2^4+11*ep*n2^4*x5-3421*ep^2-2255*ep^2
         *x5-364*ep^2*x5^2+4155*ep^2*n5+2694*ep^2*n5*x5+426*ep^2*n5*x5^2-1651*
         ep^2*n5^2-1047*ep^2*n5^2*x5-162*ep^2*n5^2*x5^2+213*ep^2*n5^3+132*ep^2*
         n5^3*x5+20*ep^2*n5^3*x5^2+4200*ep^2*n2+2832*ep^2*n2*x5+468*ep^2*n2*
         x5^2-3458*ep^2*n2*n5-2300*ep^2*n2*n5*x5-372*ep^2*n2*n5*x5^2+702*ep^2*
         n2*n5^2+456*ep^2*n2*n5^2*x5+72*ep^2*n2*n5^2*x5^2-1645*ep^2*n2^2-1133*
         ep^2*n2^2*x5-192*ep^2*n2^2*x5^2+687*ep^2*n2^2*n5+470*ep^2*n2^2*n5*x5+
         78*ep^2*n2^2*n5*x5^2+206*ep^2*n2^3+144*ep^2*n2^3*x5+25*ep^2*n2^3*x5^2+
         4174*ep^3+4155*ep^3*x5+1347*ep^3*x5^2+142*ep^3*x5^3-3380*ep^3*n5-3302*
         ep^3*n5*x5-1047*ep^3*n5*x5^2-108*ep^3*n5*x5^3+670*ep^3*n5^2+639*ep^3*
         n5^2*x5+198*ep^3*n5^2*x5^2+20*ep^3*n5^2*x5^3-3386*ep^3*n2-3458*ep^3*n2
         *x5-1150*ep^3*n2*x5^2-124*ep^3*n2*x5^3+1398*ep^3*n2*n5+1404*ep^3*n2*n5
         *x5+456*ep^3*n2*n5*x5^2+48*ep^3*n2*n5*x5^3+656*ep^3*n2^2+687*ep^3*n2^2
         *x5+235*ep^3*n2^2*x5^2+26*ep^3*n2^2*x5^3-2528*ep^4-3380*ep^4*x5-1651*
         ep^4*x5^2-349*ep^4*x5^3-27*ep^4*x5^4+1024*ep^4*n5+1340*ep^4*n5*x5+639*
         ep^4*n5*x5^2+132*ep^4*n5*x5^3+10*ep^4*n5*x5^4+1016*ep^4*n2+1398*ep^4*
         n2*x5+702*ep^4*n2*x5^2+152*ep^4*n2*x5^3+12*ep^4*n2*x5^4+608*ep^5+1024*
         ep^5*x5+670*ep^5*x5^2+213*ep^5*x5^3+33*ep^5*x5^4+2*ep^5*x5^5,-12+10*n5
         -2*n5^2+43*n2-33*n2*n5+6*n2*n5^2-58*n2^2+39*n2^2*n5-6*n2^2*n5^2+36*
         n2^3-19*n2^3*n5+2*n2^3*n5^2-10*n2^4+3*n2^4*n5+n2^5+72*ep+10*ep*x5-50*
         ep*n5-4*ep*n5*x5+8*ep*n5^2-203*ep*n2-33*ep*n2*x5+122*ep*n2*n5+12*ep*n2
         *n5*x5-16*ep*n2*n5^2+201*ep*n2^2+39*ep*n2^2*x5-94*ep*n2^2*n5-12*ep*
         n2^2*n5*x5+8*ep*n2^2*n5^2-81*ep*n2^3-19*ep*n2^3*x5+22*ep*n2^3*n5+4*ep*
         n2^3*n5*x5+11*ep*n2^4+3*ep*n2^4*x5-168*ep^2-50*ep^2*x5-2*ep^2*x5^2+90*
         ep^2*n5+16*ep^2*n5*x5-10*ep^2*n5^2+351*ep^2*n2+122*ep^2*n2*x5+6*ep^2*
         n2*x5^2-145*ep^2*n2*n5-32*ep^2*n2*n5*x5+10*ep^2*n2*n5^2-228*ep^2*n2^2-
         94*ep^2*n2^2*x5-6*ep^2*n2^2*x5^2+55*ep^2*n2^2*n5+16*ep^2*n2^2*n5*x5+45
         *ep^2*n2^3+22*ep^2*n2^3*x5+2*ep^2*n2^3*x5^2+192*ep^3+90*ep^3*x5+8*ep^3
         *x5^2-70*ep^3*n5-20*ep^3*n5*x5+4*ep^3*n5^2-265*ep^3*n2-145*ep^3*n2*x5-
         16*ep^3*n2*x5^2+56*ep^3*n2*n5+20*ep^3*n2*n5*x5+85*ep^3*n2^2+55*ep^3*
         n2^2*x5+8*ep^3*n2^2*x5^2-108*ep^4-70*ep^4*x5-10*ep^4*x5^2+20*ep^4*n5+8
         *ep^4*n5*x5+74*ep^4*n2+56*ep^4*n2*x5+10*ep^4*n2*x5^2+24*ep^5+20*ep^5*
         x5+4*ep^5*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,n2,1,0,1+n5,0,1,1,0)*(rat(36*n5-42*n5^2+16*n5^3
         -2*n5^4-96*n2*n5+86*n2*n5^2-15*n2*n5^3-4*n2*n5^4+n2*n5^5+103*n2^2*n5-
         77*n2^2*n5^2+11*n2^2*n5^3+n2^2*n5^4-52*n2^3*n5+30*n2^3*n5^2-3*n2^3*
         n5^3+12*n2^4*n5-4*n2^4*n5^2-n2^5*n5+36*ep*x5-258*ep*n5-84*ep*n5*x5+317
         *ep*n5^2+48*ep*n5^2*x5-149*ep*n5^3-8*ep*n5^3*x5+33*ep*n5^4-3*ep*n5^5-
         96*ep*n2*x5+504*ep*n2*n5+172*ep*n2*n5*x5-434*ep*n2*n5^2-45*ep*n2*n5^2*
         x5+109*ep*n2*n5^3-16*ep*n2*n5^3*x5-7*ep*n2*n5^4+5*ep*n2*n5^4*x5+103*ep
         *n2^2*x5-376*ep*n2^2*n5-154*ep*n2^2*n5*x5+224*ep*n2^2*n5^2+33*ep*n2^2*
         n5^2*x5-27*ep*n2^2*n5^3+4*ep*n2^2*n5^3*x5-52*ep*n2^3*x5+119*ep*n2^3*n5
         +60*ep*n2^3*n5*x5-39*ep*n2^3*n5^2-9*ep*n2^3*n5^2*x5+12*ep*n2^4*x5-13*
         ep*n2^4*n5-8*ep*n2^4*n5*x5-ep*n2^5*x5-258*ep^2*x5-42*ep^2*x5^2+664*
         ep^2*n5+634*ep^2*n5*x5+48*ep^2*n5*x5^2-707*ep^2*n5^2-447*ep^2*n5^2*x5-
         12*ep^2*n5^2*x5^2+260*ep^2*n5^3+132*ep^2*n5^3*x5-33*ep^2*n5^4-15*ep^2*
         n5^4*x5+504*ep^2*n2*x5+86*ep^2*n2*x5^2-926*ep^2*n2*n5-868*ep^2*n2*n5*
         x5-45*ep^2*n2*n5*x5^2+616*ep^2*n2*n5^2+327*ep^2*n2*n5^2*x5-24*ep^2*n2*
         n5^2*x5^2-98*ep^2*n2*n5^3-28*ep^2*n2*n5^3*x5+10*ep^2*n2*n5^3*x5^2-376*
         ep^2*n2^2*x5-77*ep^2*n2^2*x5^2+438*ep^2*n2^2*n5+448*ep^2*n2^2*n5*x5+33
         *ep^2*n2^2*n5*x5^2-148*ep^2*n2^2*n5^2-81*ep^2*n2^2*n5^2*x5+6*ep^2*n2^2
         *n5^2*x5^2+119*ep^2*n2^3*x5+30*ep^2*n2^3*x5^2-66*ep^2*n2^3*n5-78*ep^2*
         n2^3*n5*x5-9*ep^2*n2^3*n5*x5^2-13*ep^2*n2^4*x5-4*ep^2*n2^4*x5^2+664*
         ep^3*x5+317*ep^3*x5^2+16*ep^3*x5^3-794*ep^3*n5-1414*ep^3*n5*x5-447*
         ep^3*n5*x5^2-8*ep^3*n5*x5^3+624*ep^3*n5^2+780*ep^3*n5^2*x5+198*ep^3*
         n5^2*x5^2-126*ep^3*n5^3-132*ep^3*n5^3*x5-30*ep^3*n5^3*x5^2-926*ep^3*n2
         *x5-434*ep^3*n2*x5^2-15*ep^3*n2*x5^3+718*ep^3*n2*n5+1232*ep^3*n2*n5*x5
         +327*ep^3*n2*n5*x5^2-16*ep^3*n2*n5*x5^3-266*ep^3*n2*n5^2-294*ep^3*n2*
         n5^2*x5-42*ep^3*n2*n5^2*x5^2+10*ep^3*n2*n5^2*x5^3+438*ep^3*n2^2*x5+224
         *ep^3*n2^2*x5^2+11*ep^3*n2^2*x5^3-164*ep^3*n2^2*n5-296*ep^3*n2^2*n5*x5
         -81*ep^3*n2^2*n5*x5^2+4*ep^3*n2^2*n5*x5^3-66*ep^3*n2^3*x5-39*ep^3*n2^3
         *x5^2-3*ep^3*n2^3*x5^3-794*ep^4*x5-707*ep^4*x5^2-149*ep^4*x5^3-2*ep^4*
         x5^4+448*ep^4*n5+1248*ep^4*n5*x5+780*ep^4*n5*x5^2+132*ep^4*n5*x5^3-192
         *ep^4*n5^2-378*ep^4*n5^2*x5-198*ep^4*n5^2*x5^2-30*ep^4*n5^2*x5^3+718*
         ep^4*n2*x5+616*ep^4*n2*x5^2+109*ep^4*n2*x5^3-4*ep^4*n2*x5^4-200*ep^4*
         n2*n5-532*ep^4*n2*n5*x5-294*ep^4*n2*n5*x5^2-28*ep^4*n2*n5*x5^3+5*ep^4*
         n2*n5*x5^4-164*ep^4*n2^2*x5-148*ep^4*n2^2*x5^2-27*ep^4*n2^2*x5^3+ep^4*
         n2^2*x5^4+448*ep^5*x5+624*ep^5*x5^2+260*ep^5*x5^3+33*ep^5*x5^4-96*ep^5
         *n5-384*ep^5*n5*x5-378*ep^5*n5*x5^2-132*ep^5*n5*x5^3-15*ep^5*n5*x5^4-
         200*ep^5*n2*x5-266*ep^5*n2*x5^2-98*ep^5*n2*x5^3-7*ep^5*n2*x5^4+ep^5*n2
         *x5^5-96*ep^6*x5-192*ep^6*x5^2-126*ep^6*x5^3-33*ep^6*x5^4-3*ep^6*x5^5,
         48-64*n5+28*n5^2-4*n5^3-172*n2+206*n2*n5-80*n2*n5^2+10*n2*n5^3+244*
         n2^2-251*n2^2*n5+81*n2^2*n5^2-8*n2^2*n5^3-175*n2^3+144*n2^3*n5-34*n2^3
         *n5^2+2*n2^3*n5^3+67*n2^4-39*n2^4*n5+5*n2^4*n5^2-13*n2^5+4*n2^5*n5+
         n2^6-312*ep-64*ep*x5+352*ep*n5+56*ep*n5*x5-126*ep*n5^2-12*ep*n5^2*x5+
         14*ep*n5^3+934*ep*n2+206*ep*n2*x5-911*ep*n2*n5-160*ep*n2*n5*x5+273*ep*
         n2*n5^2+30*ep*n2*n5^2*x5-24*ep*n2*n5^3-1061*ep*n2^2-251*ep*n2^2*x5+835
         *ep*n2^2*n5+162*ep*n2^2*n5*x5-186*ep*n2^2*n5^2-24*ep*n2^2*n5^2*x5+10*
         ep*n2^2*n5^3+570*ep*n2^3+144*ep*n2^3*x5-319*ep*n2^3*n5-68*ep*n2^3*n5*
         x5+39*ep*n2^3*n5^2+6*ep*n2^3*n5^2*x5-145*ep*n2^4-39*ep*n2^4*x5+43*ep*
         n2^4*n5+10*ep*n2^4*n5*x5+14*ep*n2^5+4*ep*n2^5*x5+792*ep^2+352*ep^2*x5+
         28*ep^2*x5^2-704*ep^2*n5-252*ep^2*n5*x5-12*ep^2*n5*x5^2+182*ep^2*n5^2+
         42*ep^2*n5^2*x5-12*ep^2*n5^3-1940*ep^2*n2-911*ep^2*n2*x5-80*ep^2*n2*
         x5^2+1419*ep^2*n2*n5+546*ep^2*n2*n5*x5+30*ep^2*n2*n5*x5^2-282*ep^2*n2*
         n5^2-72*ep^2*n2*n5^2*x5+12*ep^2*n2*n5^3+1675*ep^2*n2^2+835*ep^2*n2^2*
         x5+81*ep^2*n2^2*x5^2-886*ep^2*n2^2*n5-372*ep^2*n2^2*n5*x5-24*ep^2*n2^2
         *n5*x5^2+100*ep^2*n2^2*n5^2+30*ep^2*n2^2*n5^2*x5-604*ep^2*n2^3-319*
         ep^2*n2^3*x5-34*ep^2*n2^3*x5^2+171*ep^2*n2^3*n5+78*ep^2*n2^3*n5*x5+6*
         ep^2*n2^3*n5*x5^2+77*ep^2*n2^4+43*ep^2*n2^4*x5+5*ep^2*n2^4*x5^2-984*
         ep^3-704*ep^3*x5-126*ep^3*x5^2-4*ep^3*x5^3+608*ep^3*n5+364*ep^3*n5*x5+
         42*ep^3*n5*x5^2-84*ep^3*n5^2-36*ep^3*n5^2*x5+1910*ep^3*n2+1419*ep^3*n2
         *x5+273*ep^3*n2*x5^2+10*ep^3*n2*x5^3-906*ep^3*n2*n5-564*ep^3*n2*n5*x5-
         72*ep^3*n2*n5*x5^2+84*ep^3*n2*n5^2+36*ep^3*n2*n5^2*x5-1134*ep^3*n2^2-
         886*ep^3*n2^2*x5-186*ep^3*n2^2*x5^2-8*ep^3*n2^2*x5^3+298*ep^3*n2^2*n5+
         200*ep^3*n2^2*n5*x5+30*ep^3*n2^2*n5*x5^2+208*ep^3*n2^3+171*ep^3*n2^3*
         x5+39*ep^3*n2^3*x5^2+2*ep^3*n2^3*x5^3+600*ep^4+608*ep^4*x5+182*ep^4*
         x5^2+14*ep^4*x5^3-192*ep^4*n5-168*ep^4*n5*x5-36*ep^4*n5*x5^2-876*ep^4*
         n2-906*ep^4*n2*x5-282*ep^4*n2*x5^2-24*ep^4*n2*x5^3+192*ep^4*n2*n5+168*
         ep^4*n2*n5*x5+36*ep^4*n2*n5*x5^2+276*ep^4*n2^2+298*ep^4*n2^2*x5+100*
         ep^4*n2^2*x5^2+10*ep^4*n2^2*x5^3-144*ep^5-192*ep^5*x5-84*ep^5*x5^2-12*
         ep^5*x5^3+144*ep^5*n2+192*ep^5*n2*x5+84*ep^5*n2*x5^2+12*ep^5*n2*x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,n2,1,0,1+n5,1,1,0,0)*(rat(-504*n5+984*n5^2-722*
         n5^3+246*n5^4-38*n5^5+2*n5^6+1776*n2*n5-3082*n2*n5^2+1934*n2*n5^3-509*
         n2*n5^4+35*n2*n5^5+7*n2*n5^6-n2*n5^7-2594*n2^2*n5+3919*n2^2*n5^2-2073*
         n2^2*n5^3+437*n2^2*n5^4-23*n2^2*n5^5-2*n2^2*n5^6+2030*n2^3*n5-2566*
         n2^3*n5^2+1081*n2^3*n5^3-167*n2^3*n5^4+6*n2^3*n5^5-919*n2^4*n5+909*
         n2^4*n5^2-271*n2^4*n5^3+23*n2^4*n5^4+241*n2^5*n5-165*n2^5*n5^2+26*n2^5
         *n5^3-34*n2^6*n5+12*n2^6*n5^2+2*n2^7*n5-504*ep*x5+3828*ep*n5+1968*ep*
         n5*x5-6544*ep*n5^2-2166*ep*n5^2*x5+4315*ep*n5^3+984*ep*n5^3*x5-1462*ep
         *n5^4-190*ep*n5^4*x5+302*ep*n5^5+12*ep*n5^5*x5-42*ep*n5^6+3*ep*n5^7+
         1776*ep*n2*x5-11640*ep*n2*n5-6164*ep*n2*n5*x5+17025*ep*n2*n5^2+5802*ep
         *n2*n5^2*x5-8970*ep*n2*n5^3-2036*ep*n2*n5^3*x5+2079*ep*n2*n5^4+175*ep*
         n2*n5^4*x5-204*ep*n2*n5^5+42*ep*n2*n5^5*x5+6*ep*n2*n5^6-7*ep*n2*n5^6*
         x5-2594*ep*n2^2*x5+14255*ep*n2^2*n5+7838*ep*n2^2*n5*x5-17314*ep*n2^2*
         n5^2-6219*ep*n2^2*n5^2*x5+7076*ep*n2^2*n5^3+1748*ep*n2^2*n5^3*x5-1101*
         ep*n2^2*n5^4-115*ep*n2^2*n5^4*x5+48*ep*n2^2*n5^5-12*ep*n2^2*n5^5*x5+
         2030*ep*n2^3*x5-8972*ep*n2^3*n5-5132*ep*n2^3*n5*x5+8511*ep*n2^3*n5^2+
         3243*ep*n2^3*n5^2*x5-2433*ep*n2^3*n5^3-668*ep*n2^3*n5^3*x5+200*ep*n2^3
         *n5^4+30*ep*n2^3*n5^4*x5-919*ep*n2^4*x5+3060*ep*n2^4*n5+1818*ep*n2^4*
         n5*x5-2013*ep*n2^4*n5^2-813*ep*n2^4*n5^2*x5+303*ep*n2^4*n5^3+92*ep*
         n2^4*n5^3*x5+241*ep*n2^5*x5-537*ep*n2^5*n5-330*ep*n2^5*n5*x5+183*ep*
         n2^5*n5^2+78*ep*n2^5*n5^2*x5-34*ep*n2^6*x5+38*ep*n2^6*n5+24*ep*n2^6*n5
         *x5+2*ep*n2^7*x5+3828*ep^2*x5+984*ep^2*x5^2-12212*ep^2*n5-13088*ep^2*
         n5*x5-2166*ep^2*n5*x5^2+17780*ep^2*n5^2+12945*ep^2*n5^2*x5+1476*ep^2*
         n5^2*x5^2-10087*ep^2*n5^3-5848*ep^2*n5^3*x5-380*ep^2*n5^3*x5^2+3023*
         ep^2*n5^4+1510*ep^2*n5^4*x5+30*ep^2*n5^4*x5^2-533*ep^2*n5^5-252*ep^2*
         n5^5*x5+45*ep^2*n5^6+21*ep^2*n5^6*x5-11640*ep^2*n2*x5-3082*ep^2*n2*
         x5^2+31296*ep^2*n2*n5+34050*ep^2*n2*n5*x5+5802*ep^2*n2*n5*x5^2-37063*
         ep^2*n2*n5^2-26910*ep^2*n2*n5^2*x5-3054*ep^2*n2*n5^2*x5^2+15369*ep^2*
         n2*n5^3+8316*ep^2*n2*n5^3*x5+350*ep^2*n2*n5^3*x5^2-2715*ep^2*n2*n5^4-
         1020*ep^2*n2*n5^4*x5+105*ep^2*n2*n5^4*x5^2+185*ep^2*n2*n5^5+36*ep^2*n2
         *n5^5*x5-21*ep^2*n2*n5^5*x5^2+14255*ep^2*n2^2*x5+3919*ep^2*n2^2*x5^2-
         30939*ep^2*n2^2*n5-34628*ep^2*n2^2*n5*x5-6219*ep^2*n2^2*n5*x5^2+28353*
         ep^2*n2^2*n5^2+21228*ep^2*n2^2*n5^2*x5+2622*ep^2*n2^2*n5^2*x5^2-7966*
         ep^2*n2^2*n5^3-4404*ep^2*n2^2*n5^3*x5-230*ep^2*n2^2*n5^3*x5^2+680*ep^2
         *n2^2*n5^4+240*ep^2*n2^2*n5^4*x5-30*ep^2*n2^2*n5^4*x5^2-8972*ep^2*n2^3
         *x5-2566*ep^2*n2^3*x5^2+14706*ep^2*n2^3*n5+17022*ep^2*n2^3*n5*x5+3243*
         ep^2*n2^3*n5*x5^2-9320*ep^2*n2^3*n5^2-7299*ep^2*n2^3*n5^2*x5-1002*ep^2
         *n2^3*n5^2*x5^2+1358*ep^2*n2^3*n5^3+800*ep^2*n2^3*n5^3*x5+60*ep^2*n2^3
         *n5^3*x5^2+3060*ep^2*n2^4*x5+909*ep^2*n2^4*x5^2-3361*ep^2*n2^4*n5-4026
         *ep^2*n2^4*n5*x5-813*ep^2*n2^4*n5*x5^2+1105*ep^2*n2^4*n5^2+909*ep^2*
         n2^4*n5^2*x5+138*ep^2*n2^4*n5^2*x5^2-537*ep^2*n2^5*x5-165*ep^2*n2^5*
         x5^2+296*ep^2*n2^5*n5+366*ep^2*n2^5*n5*x5+78*ep^2*n2^5*n5*x5^2+38*ep^2
         *n2^6*x5+12*ep^2*n2^6*x5^2-12212*ep^3*x5-6544*ep^3*x5^2-722*ep^3*x5^3+
         21076*ep^3*n5+35560*ep^3*n5*x5+12945*ep^3*n5*x5^2+984*ep^3*n5*x5^3-
         25112*ep^3*n5^2-30261*ep^3*n5^2*x5-8772*ep^3*n5^2*x5^2-380*ep^3*n5^2*
         x5^3+11498*ep^3*n5^3+12092*ep^3*n5^3*x5+3020*ep^3*n5^3*x5^2+40*ep^3*
         n5^3*x5^3-2636*ep^3*n5^4-2665*ep^3*n5^4*x5-630*ep^3*n5^4*x5^2+270*ep^3
         *n5^5+270*ep^3*n5^5*x5+63*ep^3*n5^5*x5^2+31296*ep^3*n2*x5+17025*ep^3*
         n2*x5^2+1934*ep^3*n2*x5^3-44000*ep^3*n2*n5-74126*ep^3*n2*n5*x5-26910*
         ep^3*n2*n5*x5^2-2036*ep^3*n2*n5*x5^3+39580*ep^3*n2*n5^2+46107*ep^3*n2*
         n5^2*x5+12474*ep^3*n2*n5^2*x5^2+350*ep^3*n2*n5^2*x5^3-11494*ep^3*n2*
         n5^3-10860*ep^3*n2*n5^3*x5-2040*ep^3*n2*n5^3*x5^2+140*ep^3*n2*n5^3*
         x5^3+1142*ep^3*n2*n5^4+925*ep^3*n2*n5^4*x5+90*ep^3*n2*n5^4*x5^2-35*
         ep^3*n2*n5^4*x5^3-30939*ep^3*n2^2*x5-17314*ep^3*n2^2*x5^2-2073*ep^3*
         n2^2*x5^3+33058*ep^3*n2^2*n5+56706*ep^3*n2^2*n5*x5+21228*ep^3*n2^2*n5*
         x5^2+1748*ep^3*n2^2*n5*x5^3-20336*ep^3*n2^2*n5^2-23898*ep^3*n2^2*n5^2*
         x5-6606*ep^3*n2^2*n5^2*x5^2-230*ep^3*n2^2*n5^2*x5^3+2950*ep^3*n2^2*
         n5^3+2720*ep^3*n2^2*n5^3*x5+480*ep^3*n2^2*n5^3*x5^2-40*ep^3*n2^2*n5^3*
         x5^3+14706*ep^3*n2^3*x5+8511*ep^3*n2^3*x5^2+1081*ep^3*n2^3*x5^3-10574*
         ep^3*n2^3*n5-18640*ep^3*n2^3*n5*x5-7299*ep^3*n2^3*n5*x5^2-668*ep^3*
         n2^3*n5*x5^3+3362*ep^3*n2^3*n5^2+4074*ep^3*n2^3*n5^2*x5+1200*ep^3*n2^3
         *n5^2*x5^2+60*ep^3*n2^3*n5^2*x5^3-3361*ep^3*n2^4*x5-2013*ep^3*n2^4*
         x5^2-271*ep^3*n2^4*x5^3+1216*ep^3*n2^4*n5+2210*ep^3*n2^4*n5*x5+909*
         ep^3*n2^4*n5*x5^2+92*ep^3*n2^4*n5*x5^3+296*ep^3*n2^5*x5+183*ep^3*n2^5*
         x5^2+26*ep^3*n2^5*x5^3+21076*ep^4*x5+17780*ep^4*x5^2+4315*ep^4*x5^3+
         246*ep^4*x5^4-21060*ep^4*n5-50224*ep^4*n5*x5-30261*ep^4*n5*x5^2-5848*
         ep^4*n5*x5^3-190*ep^4*n5*x5^4+19292*ep^4*n5^2+34494*ep^4*n5^2*x5+18138
         *ep^4*n5^2*x5^2+3020*ep^4*n5^2*x5^3+30*ep^4*n5^2*x5^4-6372*ep^4*n5^3-
         10544*ep^4*n5^3*x5-5330*ep^4*n5^3*x5^2-840*ep^4*n5^3*x5^3+828*ep^4*
         n5^4+1350*ep^4*n5^4*x5+675*ep^4*n5^4*x5^2+105*ep^4*n5^4*x5^3-44000*
         ep^4*n2*x5-37063*ep^4*n2*x5^2-8970*ep^4*n2*x5^3-509*ep^4*n2*x5^4+33936
         *ep^4*n2*n5+79160*ep^4*n2*n5*x5+46107*ep^4*n2*n5*x5^2+8316*ep^4*n2*n5*
         x5^3+175*ep^4*n2*n5*x5^4-20628*ep^4*n2*n5^2-34482*ep^4*n2*n5^2*x5-
         16290*ep^4*n2*n5^2*x5^2-2040*ep^4*n2*n5^2*x5^3+105*ep^4*n2*n5^2*x5^4+
         3156*ep^4*n2*n5^3+4568*ep^4*n2*n5^3*x5+1850*ep^4*n2*n5^3*x5^2+120*ep^4
         *n2*n5^3*x5^3-35*ep^4*n2*n5^3*x5^4+33058*ep^4*n2^2*x5+28353*ep^4*n2^2*
         x5^2+7076*ep^4*n2^2*x5^3+437*ep^4*n2^2*x5^4-17332*ep^4*n2^2*n5-40672*
         ep^4*n2^2*n5*x5-23898*ep^4*n2^2*n5*x5^2-4404*ep^4*n2^2*n5*x5^3-115*
         ep^4*n2^2*n5*x5^4+5372*ep^4*n2^2*n5^2+8850*ep^4*n2^2*n5^2*x5+4080*ep^4
         *n2^2*n5^2*x5^2+480*ep^4*n2^2*n5^2*x5^3-30*ep^4*n2^2*n5^2*x5^4-10574*
         ep^4*n2^3*x5-9320*ep^4*n2^3*x5^2-2433*ep^4*n2^3*x5^3-167*ep^4*n2^3*
         x5^4+2808*ep^4*n2^3*n5+6724*ep^4*n2^3*n5*x5+4074*ep^4*n2^3*n5*x5^2+800
         *ep^4*n2^3*n5*x5^3+30*ep^4*n2^3*n5*x5^4+1216*ep^4*n2^4*x5+1105*ep^4*
         n2^4*x5^2+303*ep^4*n2^4*x5^3+23*ep^4*n2^4*x5^4-21060*ep^5*x5-25112*
         ep^5*x5^2-10087*ep^5*x5^3-1462*ep^5*x5^4-38*ep^5*x5^5+12008*ep^5*n5+
         38584*ep^5*n5*x5+34494*ep^5*n5*x5^2+12092*ep^5*n5*x5^3+1510*ep^5*n5*
         x5^4+12*ep^5*n5*x5^5-7552*ep^5*n5^2-19116*ep^5*n5^2*x5-15816*ep^5*n5^2
         *x5^2-5330*ep^5*n5^2*x5^3-630*ep^5*n5^2*x5^4+1368*ep^5*n5^3+3312*ep^5*
         n5^3*x5+2700*ep^5*n5^3*x5^2+900*ep^5*n5^3*x5^3+105*ep^5*n5^3*x5^4+
         33936*ep^5*n2*x5+39580*ep^5*n2*x5^2+15369*ep^5*n2*x5^3+2079*ep^5*n2*
         x5^4+35*ep^5*n2*x5^5-13512*ep^5*n2*n5-41256*ep^5*n2*n5*x5-34482*ep^5*
         n2*n5*x5^2-10860*ep^5*n2*n5*x5^3-1020*ep^5*n2*n5*x5^4+42*ep^5*n2*n5*
         x5^5+4168*ep^5*n2*n5^2+9468*ep^5*n2*n5^2*x5+6852*ep^5*n2*n5^2*x5^2+
         1850*ep^5*n2*n5^2*x5^3+90*ep^5*n2*n5^2*x5^4-21*ep^5*n2*n5^2*x5^5-17332
         *ep^5*n2^2*x5-20336*ep^5*n2^2*x5^2-7966*ep^5*n2^2*x5^3-1101*ep^5*n2^2*
         x5^4-23*ep^5*n2^2*x5^5+3552*ep^5*n2^2*n5+10744*ep^5*n2^2*n5*x5+8850*
         ep^5*n2^2*n5*x5^2+2720*ep^5*n2^2*n5*x5^3+240*ep^5*n2^2*n5*x5^4-12*ep^5
         *n2^2*n5*x5^5+2808*ep^5*n2^3*x5+3362*ep^5*n2^3*x5^2+1358*ep^5*n2^3*
         x5^3+200*ep^5*n2^3*x5^4+6*ep^5*n2^3*x5^5+12008*ep^6*x5+19292*ep^6*x5^2
         +11498*ep^6*x5^3+3023*ep^6*x5^4+302*ep^6*x5^5+2*ep^6*x5^6-3520*ep^6*n5
         -15104*ep^6*n5*x5-19116*ep^6*n5*x5^2-10544*ep^6*n5*x5^3-2665*ep^6*n5*
         x5^4-252*ep^6*n5*x5^5+1152*ep^6*n5^2+4104*ep^6*n5^2*x5+4968*ep^6*n5^2*
         x5^2+2700*ep^6*n5^2*x5^3+675*ep^6*n5^2*x5^4+63*ep^6*n5^2*x5^5-13512*
         ep^6*n2*x5-20628*ep^6*n2*x5^2-11494*ep^6*n2*x5^3-2715*ep^6*n2*x5^4-204
         *ep^6*n2*x5^5+7*ep^6*n2*x5^6+2144*ep^6*n2*n5+8336*ep^6*n2*n5*x5+9468*
         ep^6*n2*n5*x5^2+4568*ep^6*n2*n5*x5^3+925*ep^6*n2*n5*x5^4+36*ep^6*n2*n5
         *x5^5-7*ep^6*n2*n5*x5^6+3552*ep^6*n2^2*x5+5372*ep^6*n2^2*x5^2+2950*
         ep^6*n2^2*x5^3+680*ep^6*n2^2*x5^4+48*ep^6*n2^2*x5^5-2*ep^6*n2^2*x5^6-
         3520*ep^7*x5-7552*ep^7*x5^2-6372*ep^7*x5^3-2636*ep^7*x5^4-533*ep^7*
         x5^5-42*ep^7*x5^6+384*ep^7*n5+2304*ep^7*n5*x5+4104*ep^7*n5*x5^2+3312*
         ep^7*n5*x5^3+1350*ep^7*n5*x5^4+270*ep^7*n5*x5^5+21*ep^7*n5*x5^6+2144*
         ep^7*n2*x5+4168*ep^7*n2*x5^2+3156*ep^7*n2*x5^3+1142*ep^7*n2*x5^4+185*
         ep^7*n2*x5^5+6*ep^7*n2*x5^6-ep^7*n2*x5^7+384*ep^8*x5+1152*ep^8*x5^2+
         1368*ep^8*x5^3+828*ep^8*x5^4+270*ep^8*x5^5+45*ep^8*x5^6+3*ep^8*x5^7,48
         -112*n5+92*n5^2-32*n5^3+4*n5^4-220*n2+490*n2*n5-378*n2*n5^2+122*n2*
         n5^3-14*n2*n5^4+416*n2^2-873*n2^2*n5+618*n2^2*n5^2-179*n2^2*n5^3+18*
         n2^2*n5^4-419*n2^3+814*n2^3*n5-510*n2^3*n5^2+125*n2^3*n5^3-10*n2^3*
         n5^4+242*n2^4-425*n2^4*n5+222*n2^4*n5^2-41*n2^4*n5^3+2*n2^4*n5^4-80*
         n2^5+123*n2^5*n5-48*n2^5*n5^2+5*n2^5*n5^3+14*n2^6-18*n2^6*n5+4*n2^6*
         n5^2-n2^7+n2^7*n5-456*ep-112*ep*x5+904*ep*n5+184*ep*n5*x5-626*ep*n5^2-
         96*ep*n5^2*x5+180*ep*n5^3+16*ep*n5^3*x5-18*ep*n5^4+1858*ep*n2+490*ep*
         n2*x5-3427*ep*n2*n5-756*ep*n2*n5*x5+2164*ep*n2*n5^2+366*ep*n2*n5^2*x5-
         555*ep*n2*n5^3-56*ep*n2*n5^3*x5+48*ep*n2*n5^4-3071*ep*n2^2-873*ep*n2^2
         *x5+5150*ep*n2^2*n5+1236*ep*n2^2*n5*x5-2859*ep*n2^2*n5^2-537*ep*n2^2*
         n5^2*x5+618*ep*n2^2*n5^3+72*ep*n2^2*n5^3*x5-42*ep*n2^2*n5^4+2644*ep*
         n2^3+814*ep*n2^3*x5-3894*ep*n2^3*n5-1020*ep*n2^3*n5*x5+1787*ep*n2^3*
         n5^2+375*ep*n2^3*n5^2*x5-291*ep*n2^3*n5^3-40*ep*n2^3*n5^3*x5+12*ep*
         n2^3*n5^4-1266*ep*n2^4-425*ep*n2^4*x5+1549*ep*n2^4*n5+444*ep*n2^4*n5*
         x5-523*ep*n2^4*n5^2-123*ep*n2^4*n5^2*x5+48*ep*n2^4*n5^3+8*ep*n2^4*n5^3
         *x5+332*ep*n2^5+123*ep*n2^5*x5-305*ep*n2^5*n5-96*ep*n2^5*n5*x5+57*ep*
         n2^5*n5^2+15*ep*n2^5*n5^2*x5-43*ep*n2^6-18*ep*n2^6*x5+23*ep*n2^6*n5+8*
         ep*n2^6*n5*x5+2*ep*n2^7+ep*n2^7*x5+1824*ep^2+904*ep^2*x5+92*ep^2*x5^2-
         2992*ep^2*n5-1252*ep^2*n5*x5-96*ep^2*n5*x5^2+1672*ep^2*n5^2+540*ep^2*
         n5^2*x5+24*ep^2*n5^2*x5^2-370*ep^2*n5^3-72*ep^2*n5^3*x5+26*ep^2*n5^4-
         6502*ep^2*n2-3427*ep^2*n2*x5-378*ep^2*n2*x5^2+9638*ep^2*n2*n5+4328*
         ep^2*n2*n5*x5+366*ep^2*n2*n5*x5^2-4729*ep^2*n2*n5^2-1665*ep^2*n2*n5^2*
         x5-84*ep^2*n2*n5^2*x5^2+881*ep^2*n2*n5^3+192*ep^2*n2*n5^3*x5-48*ep^2*
         n2*n5^4+9154*ep^2*n2^2+5150*ep^2*n2^2*x5+618*ep^2*n2^2*x5^2-11810*ep^2
         *n2^2*n5-5718*ep^2*n2^2*n5*x5-537*ep^2*n2^2*n5*x5^2+4788*ep^2*n2^2*
         n5^2+1854*ep^2*n2^2*n5^2*x5+108*ep^2*n2^2*n5^2*x5^2-674*ep^2*n2^2*n5^3
         -168*ep^2*n2^2*n5^3*x5+22*ep^2*n2^2*n5^4-6461*ep^2*n2^3-3894*ep^2*n2^3
         *x5-510*ep^2*n2^3*x5^2+6821*ep^2*n2^3*n5+3574*ep^2*n2^3*n5*x5+375*ep^2
         *n2^3*n5*x5^2-2033*ep^2*n2^3*n5^2-873*ep^2*n2^3*n5^2*x5-60*ep^2*n2^3*
         n5^2*x5^2+163*ep^2*n2^3*n5^3+48*ep^2*n2^3*n5^3*x5+2390*ep^2*n2^4+1549*
         ep^2*n2^4*x5+222*ep^2*n2^4*x5^2-1842*ep^2*n2^4*n5-1046*ep^2*n2^4*n5*x5
         -123*ep^2*n2^4*n5*x5^2+302*ep^2*n2^4*n5^2+144*ep^2*n2^4*n5^2*x5+12*
         ep^2*n2^4*n5^2*x5^2-435*ep^2*n2^5-305*ep^2*n2^5*x5-48*ep^2*n2^5*x5^2+
         185*ep^2*n2^5*n5+114*ep^2*n2^5*n5*x5+15*ep^2*n2^5*n5*x5^2+30*ep^2*n2^6
         +23*ep^2*n2^6*x5+4*ep^2*n2^6*x5^2-3984*ep^3-2992*ep^3*x5-626*ep^3*x5^2
         -32*ep^3*x5^3+5200*ep^3*n5+3344*ep^3*n5*x5+540*ep^3*n5*x5^2+16*ep^3*n5
         *x5^3-2194*ep^3*n5^2-1110*ep^3*n5^2*x5-108*ep^3*n5^2*x5^2+330*ep^3*
         n5^3+104*ep^3*n5^3*x5-12*ep^3*n5^4+12166*ep^3*n2+9638*ep^3*n2*x5+2164*
         ep^3*n2*x5^2+122*ep^3*n2*x5^3-13835*ep^3*n2*n5-9458*ep^3*n2*n5*x5-1665
         *ep^3*n2*n5*x5^2-56*ep^3*n2*n5*x5^3+4857*ep^3*n2*n5^2+2643*ep^3*n2*
         n5^2*x5+288*ep^3*n2*n5^2*x5^2-558*ep^3*n2*n5^3-192*ep^3*n2*n5^3*x5+12*
         ep^3*n2*n5^4-14071*ep^3*n2^2-11810*ep^3*n2^2*x5-2859*ep^3*n2^2*x5^2-
         179*ep^3*n2^2*x5^3+13089*ep^3*n2^2*n5+9576*ep^3*n2^2*n5*x5+1854*ep^3*
         n2^2*n5*x5^2+72*ep^3*n2^2*n5*x5^3-3410*ep^3*n2^2*n5^2-2022*ep^3*n2^2*
         n5^2*x5-252*ep^3*n2^2*n5^2*x5^2+228*ep^3*n2^2*n5^3+88*ep^3*n2^2*n5^3*
         x5+7644*ep^3*n2^3+6821*ep^3*n2^3*x5+1787*ep^3*n2^3*x5^2+125*ep^3*n2^3*
         x5^3-5167*ep^3*n2^3*n5-4066*ep^3*n2^3*n5*x5-873*ep^3*n2^3*n5*x5^2-40*
         ep^3*n2^3*n5*x5^3+747*ep^3*n2^3*n5^2+489*ep^3*n2^3*n5^2*x5+72*ep^3*
         n2^3*n5^2*x5^2-1937*ep^3*n2^4-1842*ep^3*n2^4*x5-523*ep^3*n2^4*x5^2-41*
         ep^3*n2^4*x5^3+713*ep^3*n2^4*n5+604*ep^3*n2^4*n5*x5+144*ep^3*n2^4*n5*
         x5^2+8*ep^3*n2^4*n5*x5^3+182*ep^3*n2^5+185*ep^3*n2^5*x5+57*ep^3*n2^5*
         x5^2+5*ep^3*n2^5*x5^3+5136*ep^4+5200*ep^4*x5+1672*ep^4*x5^2+180*ep^4*
         x5^3+4*ep^4*x5^4-5008*ep^4*n5-4388*ep^4*n5*x5-1110*ep^4*n5*x5^2-72*
         ep^4*n5*x5^3+1416*ep^4*n5^2+990*ep^4*n5^2*x5+156*ep^4*n5^2*x5^2-108*
         ep^4*n5^3-48*ep^4*n5^3*x5-13054*ep^4*n2-13835*ep^4*n2*x5-4729*ep^4*n2*
         x5^2-555*ep^4*n2*x5^3-14*ep^4*n2*x5^4+10542*ep^4*n2*n5+9714*ep^4*n2*n5
         *x5+2643*ep^4*n2*n5*x5^2+192*ep^4*n2*n5*x5^3-2274*ep^4*n2*n5^2-1674*
         ep^4*n2*n5^2*x5-288*ep^4*n2*n5^2*x5^2+108*ep^4*n2*n5^3+48*ep^4*n2*n5^3
         *x5+11724*ep^4*n2^2+13089*ep^4*n2^2*x5+4788*ep^4*n2^2*x5^2+618*ep^4*
         n2^2*x5^3+18*ep^4*n2^2*x5^4-6956*ep^4*n2^2*n5-6820*ep^4*n2^2*n5*x5-
         2022*ep^4*n2^2*n5*x5^2-168*ep^4*n2^2*n5*x5^3+858*ep^4*n2^2*n5^2+684*
         ep^4*n2^2*n5^2*x5+132*ep^4*n2^2*n5^2*x5^2-4376*ep^4*n2^3-5167*ep^4*
         n2^3*x5-2033*ep^4*n2^3*x5^2-291*ep^4*n2^3*x5^3-10*ep^4*n2^3*x5^4+1422*
         ep^4*n2^3*n5+1494*ep^4*n2^3*n5*x5+489*ep^4*n2^3*n5*x5^2+48*ep^4*n2^3*
         n5*x5^3+570*ep^4*n2^4+713*ep^4*n2^4*x5+302*ep^4*n2^4*x5^2+48*ep^4*n2^4
         *x5^3+2*ep^4*n2^4*x5^4-3912*ep^5-5008*ep^5*x5-2194*ep^5*x5^2-370*ep^5*
         x5^3-18*ep^5*x5^4+2536*ep^5*n5+2832*ep^5*n5*x5+990*ep^5*n5*x5^2+104*
         ep^5*n5*x5^3-360*ep^5*n5^2-324*ep^5*n5^2*x5-72*ep^5*n5^2*x5^2+7936*
         ep^5*n2+10542*ep^5*n2*x5+4857*ep^5*n2*x5^2+881*ep^5*n2*x5^3+48*ep^5*n2
         *x5^4-3936*ep^5*n2*n5-4548*ep^5*n2*n5*x5-1674*ep^5*n2*n5*x5^2-192*ep^5
         *n2*n5*x5^3+360*ep^5*n2*n5^2+324*ep^5*n2*n5^2*x5+72*ep^5*n2*n5^2*x5^2-
         4992*ep^5*n2^2-6956*ep^5*n2^2*x5-3410*ep^5*n2^2*x5^2-674*ep^5*n2^2*
         x5^3-42*ep^5*n2^2*x5^4+1400*ep^5*n2^2*n5+1716*ep^5*n2^2*n5*x5+684*ep^5
         *n2^2*n5*x5^2+88*ep^5*n2^2*n5*x5^3+968*ep^5*n2^3+1422*ep^5*n2^3*x5+747
         *ep^5*n2^3*x5^2+163*ep^5*n2^3*x5^3+12*ep^5*n2^3*x5^4+1632*ep^6+2536*
         ep^6*x5+1416*ep^6*x5^2+330*ep^6*x5^3+26*ep^6*x5^4-528*ep^6*n5-720*ep^6
         *n5*x5-324*ep^6*n5*x5^2-48*ep^6*n5*x5^3-2472*ep^6*n2-3936*ep^6*n2*x5-
         2274*ep^6*n2*x5^2-558*ep^6*n2*x5^3-48*ep^6*n2*x5^4+528*ep^6*n2*n5+720*
         ep^6*n2*n5*x5+324*ep^6*n2*n5*x5^2+48*ep^6*n2*n5*x5^3+840*ep^6*n2^2+
         1400*ep^6*n2^2*x5+858*ep^6*n2^2*x5^2+228*ep^6*n2^2*x5^3+22*ep^6*n2^2*
         x5^4-288*ep^7-528*ep^7*x5-360*ep^7*x5^2-108*ep^7*x5^3-12*ep^7*x5^4+288
         *ep^7*n2+528*ep^7*n2*x5+360*ep^7*n2*x5^2+108*ep^7*n2*x5^3+12*ep^7*n2*
         x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,n2,1,0,n5,1,1,1,0)*(rat(-162+297*n5-216*n5^2+78
         *n5^3-14*n5^4+n5^5+375*n2-569*n2*n5+325*n2*n5^2-83*n2*n5^3+8*n2*n5^4-
         323*n2^2+373*n2^2*n5-145*n2^2*n5^2+19*n2^2*n5^3+132*n2^3-102*n2^3*n5+
         20*n2^3*n5^2-26*n2^4+10*n2^4*n5+2*n2^5+816*ep+297*ep*x5-1169*ep*n5-432
         *ep*n5*x5+623*ep*n5^2+234*ep*n5^2*x5-147*ep*n5^3-56*ep*n5^3*x5+13*ep*
         n5^4+5*ep*n5^4*x5-1573*ep*n2-569*ep*n2*x5+1776*ep*n2*n5+650*ep*n2*n5*
         x5-673*ep*n2*n5^2-249*ep*n2*n5^2*x5+86*ep*n2*n5^3+32*ep*n2*n5^3*x5+
         1035*ep*n2^2+373*ep*n2^2*x5-796*ep*n2^2*n5-290*ep*n2^2*n5*x5+155*ep*
         n2^2*n5^2+57*ep*n2^2*n5^2*x5-284*ep*n2^3-102*ep*n2^3*x5+110*ep*n2^3*n5
         +40*ep*n2^3*n5*x5+28*ep*n2^4+10*ep*n2^4*x5-1562*ep^2-1169*ep^2*x5-216*
         ep^2*x5^2+1624*ep^2*n5+1246*ep^2*n5*x5+234*ep^2*n5*x5^2-558*ep^2*n5^2-
         441*ep^2*n5^2*x5-84*ep^2*n5^2*x5^2+64*ep^2*n5^3+52*ep^2*n5^3*x5+10*
         ep^2*n5^3*x5^2+2402*ep^2*n2+1776*ep^2*n2*x5+325*ep^2*n2*x5^2-1792*ep^2
         *n2*n5-1346*ep^2*n2*n5*x5-249*ep^2*n2*n5*x5^2+338*ep^2*n2*n5^2+258*
         ep^2*n2*n5^2*x5+48*ep^2*n2*n5^2*x5^2-1082*ep^2*n2^2-796*ep^2*n2^2*x5-
         145*ep^2*n2^2*x5^2+416*ep^2*n2^2*n5+310*ep^2*n2^2*n5*x5+57*ep^2*n2^2*
         n5*x5^2+150*ep^2*n2^3+110*ep^2*n2^3*x5+20*ep^2*n2^3*x5^2+1388*ep^3+
         1624*ep^3*x5+623*ep^3*x5^2+78*ep^3*x5^3-912*ep^3*n5-1116*ep^3*n5*x5-
         441*ep^3*n5*x5^2-56*ep^3*n5*x5^3+148*ep^3*n5^2+192*ep^3*n5^2*x5+78*
         ep^3*n5^2*x5^2+10*ep^3*n5^2*x5^3-1572*ep^3*n2-1792*ep^3*n2*x5-673*ep^3
         *n2*x5^2-83*ep^3*n2*x5^3+580*ep^3*n2*n5+676*ep^3*n2*n5*x5+258*ep^3*n2*
         n5*x5^2+32*ep^3*n2*n5*x5^3+368*ep^3*n2^2+416*ep^3*n2^2*x5+155*ep^3*
         n2^2*x5^2+19*ep^3*n2^2*x5^3-544*ep^4-912*ep^4*x5-558*ep^4*x5^2-147*
         ep^4*x5^3-14*ep^4*x5^4+160*ep^4*n5+296*ep^4*n5*x5+192*ep^4*n5*x5^2+52*
         ep^4*n5*x5^3+5*ep^4*n5*x5^4+368*ep^4*n2+580*ep^4*n2*x5+338*ep^4*n2*
         x5^2+86*ep^4*n2*x5^3+8*ep^4*n2*x5^4+64*ep^5+160*ep^5*x5+148*ep^5*x5^2+
         64*ep^5*x5^3+13*ep^5*x5^4+ep^5*x5^5,-12+22*n5-12*n5^2+2*n5^3+31*n2-54*
         n2*n5+27*n2*n5^2-4*n2*n5^3-27*n2^2+43*n2^2*n5-18*n2^2*n5^2+2*n2^2*n5^3
         +9*n2^3-12*n2^3*n5+3*n2^3*n5^2-n2^4+n2^4*n5+72*ep+22*ep*x5-98*ep*n5-24
         *ep*n5*x5+38*ep*n5^2+6*ep*n5^2*x5-4*ep*n5^3-155*ep*n2-54*ep*n2*x5+185*
         ep*n2*n5+54*ep*n2*n5*x5-58*ep*n2*n5^2-12*ep*n2*n5^2*x5+4*ep*n2*n5^3+
         108*ep*n2^2+43*ep*n2^2*x5-102*ep*n2^2*n5-36*ep*n2^2*n5*x5+20*ep*n2^2*
         n5^2+6*ep*n2^2*n5^2*x5-27*ep*n2^3-12*ep*n2^3*x5+15*ep*n2^3*n5+6*ep*
         n2^3*n5*x5+2*ep*n2^4+ep*n2^4*x5-156*ep^2-98*ep^2*x5-12*ep^2*x5^2+140*
         ep^2*n5+76*ep^2*n5*x5+6*ep^2*n5*x5^2-28*ep^2*n5^2-12*ep^2*n5^2*x5+272*
         ep^2*n2+185*ep^2*n2*x5+27*ep^2*n2*x5^2-198*ep^2*n2*n5-116*ep^2*n2*n5*
         x5-12*ep^2*n2*n5*x5^2+28*ep^2*n2*n5^2+12*ep^2*n2*n5^2*x5-134*ep^2*n2^2
         -102*ep^2*n2^2*x5-18*ep^2*n2^2*x5^2+58*ep^2*n2^2*n5+40*ep^2*n2^2*n5*x5
         +6*ep^2*n2^2*n5*x5^2+18*ep^2*n2^3+15*ep^2*n2^3*x5+3*ep^2*n2^3*x5^2+144
         *ep^3+140*ep^3*x5+38*ep^3*x5^2+2*ep^3*x5^3-64*ep^3*n5-56*ep^3*n5*x5-12
         *ep^3*n5*x5^2-196*ep^3*n2-198*ep^3*n2*x5-58*ep^3*n2*x5^2-4*ep^3*n2*
         x5^3+64*ep^3*n2*n5+56*ep^3*n2*n5*x5+12*ep^3*n2*n5*x5^2+52*ep^3*n2^2+58
         *ep^3*n2^2*x5+20*ep^3*n2^2*x5^2+2*ep^3*n2^2*x5^3-48*ep^4-64*ep^4*x5-28
         *ep^4*x5^2-4*ep^4*x5^3+48*ep^4*n2+64*ep^4*n2*x5+28*ep^4*n2*x5^2+4*ep^4
         *n2*x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,n2,1,1,1+n5,1,0,0,0)*(rat(27*n5-36*n5^2+15*n5^3
         -2*n5^4-36*n2*n5+33*n2*n5^2-7*n2*n5^3+15*n2^2*n5-7*n2^2*n5^2-2*n2^3*n5
         +27*ep*x5-99*ep*n5-72*ep*n5*x5+87*ep*n5^2+45*ep*n5^2*x5-18*ep*n5^3-8*
         ep*n5^3*x5-36*ep*n2*x5+87*ep*n2*n5+66*ep*n2*n5*x5-39*ep*n2*n5^2-21*ep*
         n2*n5^2*x5+15*ep*n2^2*x5-18*ep*n2^2*n5-14*ep*n2^2*n5*x5-2*ep*n2^3*x5-
         99*ep^2*x5-36*ep^2*x5^2+120*ep^2*n5+174*ep^2*n5*x5+45*ep^2*n5*x5^2-52*
         ep^2*n5^2-54*ep^2*n5^2*x5-12*ep^2*n5^2*x5^2+87*ep^2*n2*x5+33*ep^2*n2*
         x5^2-52*ep^2*n2*n5-78*ep^2*n2*n5*x5-21*ep^2*n2*n5*x5^2-18*ep^2*n2^2*x5
         -7*ep^2*n2^2*x5^2+120*ep^3*x5+87*ep^3*x5^2+15*ep^3*x5^3-48*ep^3*n5-104
         *ep^3*n5*x5-54*ep^3*n5*x5^2-8*ep^3*n5*x5^3-52*ep^3*n2*x5-39*ep^3*n2*
         x5^2-7*ep^3*n2*x5^3-48*ep^4*x5-52*ep^4*x5^2-18*ep^4*x5^3-2*ep^4*x5^4,2
         -3*n5+n5^2-5*n2+7*n2*n5-2*n2*n5^2+4*n2^2-5*n2^2*n5+n2^2*n5^2-n2^3+n2^3
         *n5-12*ep-3*ep*x5+13*ep*n5+2*ep*n5*x5-3*ep*n5^2+23*ep*n2+7*ep*n2*x5-20
         *ep*n2*n5-4*ep*n2*n5*x5+3*ep*n2*n5^2-13*ep*n2^2-5*ep*n2^2*x5+7*ep*n2^2
         *n5+2*ep*n2^2*n5*x5+2*ep*n2^3+ep*n2^3*x5+26*ep^2+13*ep^2*x5+ep^2*x5^2-
         18*ep^2*n5-6*ep^2*n5*x5+2*ep^2*n5^2-34*ep^2*n2-20*ep^2*n2*x5-2*ep^2*n2
         *x5^2+14*ep^2*n2*n5+6*ep^2*n2*n5*x5+10*ep^2*n2^2+7*ep^2*n2^2*x5+ep^2*
         n2^2*x5^2-24*ep^3-18*ep^3*x5-3*ep^3*x5^2+8*ep^3*n5+4*ep^3*n5*x5+16*
         ep^3*n2+14*ep^3*n2*x5+3*ep^3*n2*x5^2+8*ep^4+8*ep^4*x5+2*ep^4*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,n2,1,1,n5,1,1,0,0)*(rat(78-128*n5+76*n5^2-20*
         n5^3+2*n5^4-197*n2+262*n2*n5-116*n2*n5^2+20*n2*n5^3-n2*n5^4+195*n2^2-
         199*n2^2*n5+59*n2^2*n5^2-5*n2^2*n5^3-94*n2^3+66*n2^3*n5-10*n2^3*n5^2+
         22*n2^4-8*n2^4*n5-2*n2^5-462*ep-128*ep*x5+596*ep*n5+152*ep*n5*x5-260*
         ep*n5^2-60*ep*n5^2*x5+44*ep*n5^3+8*ep*n5^3*x5-2*ep*n5^4+932*ep*n2+262*
         ep*n2*x5-912*ep*n2*n5-232*ep*n2*n5*x5+264*ep*n2*n5^2+60*ep*n2*n5^2*x5-
         22*ep*n2*n5^3-4*ep*n2*n5^3*x5-691*ep*n2^2-199*ep*n2^2*x5+460*ep*n2^2*
         n5+118*ep*n2^2*n5*x5-67*ep*n2^2*n5^2-15*ep*n2^2*n5^2*x5+222*ep*n2^3+66
         *ep*n2^3*x5-76*ep*n2^3*n5-20*ep*n2^3*n5*x5-26*ep*n2^4-8*ep*n2^4*x5+
         1086*ep^2+596*ep^2*x5+76*ep^2*x5^2-1032*ep^2*n5-520*ep^2*n5*x5-60*ep^2
         *n5*x5^2+294*ep^2*n5^2+132*ep^2*n5^2*x5+12*ep^2*n5^2*x5^2-24*ep^2*n5^3
         -8*ep^2*n5^3*x5-1641*ep^2*n2-912*ep^2*n2*x5-116*ep^2*n2*x5^2+1050*ep^2
         *n2*n5+528*ep^2*n2*n5*x5+60*ep^2*n2*n5*x5^2-149*ep^2*n2*n5^2-66*ep^2*
         n2*n5^2*x5-6*ep^2*n2*n5^2*x5^2+810*ep^2*n2^2+460*ep^2*n2^2*x5+59*ep^2*
         n2^2*x5^2-264*ep^2*n2^2*n5-134*ep^2*n2^2*n5*x5-15*ep^2*n2^2*n5*x5^2-
         130*ep^2*n2^3-76*ep^2*n2^3*x5-10*ep^2*n2^3*x5^2-1266*ep^3-1032*ep^3*x5
         -260*ep^3*x5^2-20*ep^3*x5^3+788*ep^3*n5+588*ep^3*n5*x5+132*ep^3*n5*
         x5^2+8*ep^3*n5*x5^3-110*ep^3*n5^2-72*ep^3*n5^2*x5-12*ep^3*n5^2*x5^2+
         1274*ep^3*n2+1050*ep^3*n2*x5+264*ep^3*n2*x5^2+20*ep^3*n2*x5^3-400*ep^3
         *n2*n5-298*ep^3*n2*n5*x5-66*ep^3*n2*n5*x5^2-4*ep^3*n2*n5*x5^3-314*ep^3
         *n2^2-264*ep^3*n2^2*x5-67*ep^3*n2^2*x5^2-5*ep^3*n2^2*x5^3+732*ep^4+788
         *ep^4*x5+294*ep^4*x5^2+44*ep^4*x5^3+2*ep^4*x5^4-224*ep^4*n5-220*ep^4*
         n5*x5-72*ep^4*n5*x5^2-8*ep^4*n5*x5^3-368*ep^4*n2-400*ep^4*n2*x5-149*
         ep^4*n2*x5^2-22*ep^4*n2*x5^3-ep^4*n2*x5^4-168*ep^5-224*ep^5*x5-110*
         ep^5*x5^2-24*ep^5*x5^3-2*ep^5*x5^4,-2+3*n5-n5^2+7*n2-10*n2*n5+3*n2*
         n5^2-9*n2^2+12*n2^2*n5-3*n2^2*n5^2+5*n2^3-6*n2^3*n5+n2^3*n5^2-n2^4+
         n2^4*n5+14*ep+3*ep*x5-16*ep*n5-2*ep*n5*x5+4*ep*n5^2-40*ep*n2-10*ep*n2*
         x5+40*ep*n2*n5+6*ep*n2*n5*x5-8*ep*n2*n5^2+40*ep*n2^2+12*ep*n2^2*x5-32*
         ep*n2^2*n5-6*ep*n2^2*n5*x5+4*ep*n2^2*n5^2-16*ep*n2^3-6*ep*n2^3*x5+8*ep
         *n2^3*n5+2*ep*n2^3*n5*x5+2*ep*n2^4+ep*n2^4*x5-38*ep^2-16*ep^2*x5-ep^2*
         x5^2+31*ep^2*n5+8*ep^2*n5*x5-5*ep^2*n5^2+83*ep^2*n2+40*ep^2*n2*x5+3*
         ep^2*n2*x5^2-52*ep^2*n2*n5-16*ep^2*n2*n5*x5+5*ep^2*n2*n5^2-57*ep^2*
         n2^2-32*ep^2*n2^2*x5-3*ep^2*n2^2*x5^2+21*ep^2*n2^2*n5+8*ep^2*n2^2*n5*
         x5+12*ep^2*n2^3+8*ep^2*n2^3*x5+ep^2*n2^3*x5^2+50*ep^3+31*ep^3*x5+4*
         ep^3*x5^2-26*ep^3*n5-10*ep^3*n5*x5+2*ep^3*n5^2-74*ep^3*n2-52*ep^3*n2*
         x5-8*ep^3*n2*x5^2+22*ep^3*n2*n5+10*ep^3*n2*n5*x5+26*ep^3*n2^2+21*ep^3*
         n2^2*x5+4*ep^3*n2^2*x5^2-32*ep^4-26*ep^4*x5-5*ep^4*x5^2+8*ep^4*n5+4*
         ep^4*n5*x5+24*ep^4*n2+22*ep^4*n2*x5+5*ep^4*n2*x5^2+8*ep^5+8*ep^5*x5+2*
         ep^5*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,-1+n2,1,0,1+n5,0,1,1,0)*(rat(18*n5-12*n5^2+2*
         n5^3-33*n2*n5+14*n2*n5^2-n2*n5^3+24*n2^2*n5-6*n2^2*n5^2-8*n2^3*n5+n2^3
         *n5^2+n2^4*n5+18*ep*x5-126*ep*n5-24*ep*n5*x5+66*ep*n5^2+6*ep*n5^2*x5-8
         *ep*n5^3-33*ep*n2*x5+172*ep*n2*n5+28*ep*n2*n5*x5-55*ep*n2*n5^2-3*ep*n2
         *n5^2*x5+3*ep*n2*n5^3+24*ep*n2^2*x5-81*ep*n2^2*n5-12*ep*n2^2*n5*x5+12*
         ep*n2^2*n5^2-8*ep*n2^3*x5+13*ep*n2^3*n5+2*ep*n2^3*n5*x5+ep*n2^4*x5-126
         *ep^2*x5-12*ep^2*x5^2+292*ep^2*n5+132*ep^2*n5*x5+6*ep^2*n5*x5^2-102*
         ep^2*n5^2-24*ep^2*n5^2*x5+6*ep^2*n5^3+172*ep^2*n2*x5+14*ep^2*n2*x5^2-
         266*ep^2*n2*n5-110*ep^2*n2*n5*x5-3*ep^2*n2*n5*x5^2+44*ep^2*n2*n5^2+9*
         ep^2*n2*n5^2*x5-81*ep^2*n2^2*x5-6*ep^2*n2^2*x5^2+62*ep^2*n2^2*n5+24*
         ep^2*n2^2*n5*x5+13*ep^2*n2^3*x5+ep^2*n2^3*x5^2+292*ep^3*x5+66*ep^3*
         x5^2+2*ep^3*x5^3-280*ep^3*n5-204*ep^3*n5*x5-24*ep^3*n5*x5^2+48*ep^3*
         n5^2+18*ep^3*n5^2*x5-266*ep^3*n2*x5-55*ep^3*n2*x5^2-ep^3*n2*x5^3+128*
         ep^3*n2*n5+88*ep^3*n2*n5*x5+9*ep^3*n2*n5*x5^2+62*ep^3*n2^2*x5+12*ep^3*
         n2^2*x5^2-280*ep^4*x5-102*ep^4*x5^2-8*ep^4*x5^3+96*ep^4*n5+96*ep^4*n5*
         x5+18*ep^4*n5*x5^2+128*ep^4*n2*x5+44*ep^4*n2*x5^2+3*ep^4*n2*x5^3+96*
         ep^5*x5+48*ep^5*x5^2+6*ep^5*x5^3,-6+2*n5+23*n2-7*n2*n5-34*n2^2+9*n2^2*
         n5+24*n2^3-5*n2^3*n5-8*n2^4+n2^4*n5+n2^5+33*ep+2*ep*x5-9*ep*n5-102*ep*
         n2-7*ep*n2*x5+24*ep*n2*n5+114*ep*n2^2+9*ep*n2^2*x5-21*ep*n2^2*n5-54*ep
         *n2^3-5*ep*n2^3*x5+6*ep*n2^3*n5+9*ep*n2^4+ep*n2^4*x5-66*ep^2-9*ep^2*x5
         +13*ep^2*n5+157*ep^2*n2+24*ep^2*n2*x5-24*ep^2*n2*n5-120*ep^2*n2^2-21*
         ep^2*n2^2*x5+11*ep^2*n2^2*n5+29*ep^2*n2^3+6*ep^2*n2^3*x5+57*ep^3+13*
         ep^3*x5-6*ep^3*n5-96*ep^3*n2-24*ep^3*n2*x5+6*ep^3*n2*n5+39*ep^3*n2^2+
         11*ep^3*n2^2*x5-18*ep^4-6*ep^4*x5+18*ep^4*n2+6*ep^4*n2*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,-1+n2,1,0,1+n5,1,1,0,0)*(rat(4*n5-10*n5^2+8*
         n5^3-2*n5^4-12*n2*n5+32*n2*n5^2-29*n2*n5^3+10*n2*n5^4-n2*n5^5+11*n2^2*
         n5-30*n2^2*n5^2+28*n2^2*n5^3-10*n2^2*n5^4+n2^2*n5^5-3*n2^3*n5+8*n2^3*
         n5^2-7*n2^3*n5^3+2*n2^3*n5^4+4*ep*x5+146*ep*n5-20*ep*n5*x5-209*ep*n5^2
         +24*ep*n5^2*x5+113*ep*n5^3-8*ep*n5^3*x5-29*ep*n5^4+3*ep*n5^5-12*ep*n2*
         x5-364*ep*n2*n5+64*ep*n2*n5*x5+406*ep*n2*n5^2-87*ep*n2*n5^2*x5-153*ep*
         n2*n5^3+40*ep*n2*n5^3*x5+25*ep*n2*n5^4-5*ep*n2*n5^4*x5-2*ep*n2*n5^5+11
         *ep*n2^2*x5+365*ep*n2^2*n5-60*ep*n2^2*n5*x5-303*ep*n2^2*n5^2+84*ep*
         n2^2*n5^2*x5+66*ep*n2^2*n5^3-40*ep*n2^2*n5^3*x5-2*ep*n2^2*n5^4+5*ep*
         n2^2*n5^4*x5-3*ep*n2^3*x5-183*ep*n2^3*n5+16*ep*n2^3*n5*x5+108*ep*n2^3*
         n5^2-21*ep*n2^3*n5^2*x5-13*ep*n2^3*n5^3+8*ep*n2^3*n5^3*x5+44*ep*n2^4*
         n5-14*ep*n2^4*n5^2-4*ep*n2^5*n5+146*ep^2*x5-10*ep^2*x5^2-894*ep^2*n5-
         418*ep^2*n5*x5+24*ep^2*n5*x5^2+1056*ep^2*n5^2+339*ep^2*n5^2*x5-12*ep^2
         *n5^2*x5^2-433*ep^2*n5^3-116*ep^2*n5^3*x5+70*ep^2*n5^4+15*ep^2*n5^4*x5
         -3*ep^2*n5^5-364*ep^2*n2*x5+32*ep^2*n2*x5^2+1786*ep^2*n2*n5+812*ep^2*
         n2*n5*x5-87*ep^2*n2*n5*x5^2-1596*ep^2*n2*n5^2-459*ep^2*n2*n5^2*x5+60*
         ep^2*n2*n5^2*x5^2+443*ep^2*n2*n5^3+100*ep^2*n2*n5^3*x5-10*ep^2*n2*n5^3
         *x5^2-39*ep^2*n2*n5^4-10*ep^2*n2*n5^4*x5+365*ep^2*n2^2*x5-30*ep^2*n2^2
         *x5^2-1318*ep^2*n2^2*n5-606*ep^2*n2^2*n5*x5+84*ep^2*n2^2*n5*x5^2+791*
         ep^2*n2^2*n5^2+198*ep^2*n2^2*n5^2*x5-60*ep^2*n2^2*n5^2*x5^2-107*ep^2*
         n2^2*n5^3-8*ep^2*n2^2*n5^3*x5+10*ep^2*n2^2*n5^3*x5^2-183*ep^2*n2^3*x5+
         8*ep^2*n2^3*x5^2+426*ep^2*n2^3*n5+216*ep^2*n2^3*n5*x5-21*ep^2*n2^3*n5*
         x5^2-132*ep^2*n2^3*n5^2-39*ep^2*n2^3*n5^2*x5+12*ep^2*n2^3*n5^2*x5^2+44
         *ep^2*n2^4*x5-50*ep^2*n2^4*n5-28*ep^2*n2^4*n5*x5-4*ep^2*n2^5*x5-894*
         ep^3*x5-209*ep^3*x5^2+8*ep^3*x5^3+2054*ep^3*n5+2112*ep^3*n5*x5+339*
         ep^3*n5*x5^2-8*ep^3*n5*x5^3-1829*ep^3*n5^2-1299*ep^3*n5^2*x5-174*ep^3*
         n5^2*x5^2+498*ep^3*n5^3+280*ep^3*n5^3*x5+30*ep^3*n5^3*x5^2-39*ep^3*
         n5^4-15*ep^3*n5^4*x5+1786*ep^3*n2*x5+406*ep^3*n2*x5^2-29*ep^3*n2*x5^3-
         3080*ep^3*n2*n5-3192*ep^3*n2*n5*x5-459*ep^3*n2*n5*x5^2+40*ep^3*n2*n5*
         x5^3+1860*ep^3*n2*n5^2+1329*ep^3*n2*n5^2*x5+150*ep^3*n2*n5^2*x5^2-10*
         ep^3*n2*n5^2*x5^3-262*ep^3*n2*n5^3-156*ep^3*n2*n5^3*x5-20*ep^3*n2*n5^3
         *x5^2-1318*ep^3*n2^2*x5-303*ep^3*n2^2*x5^2+28*ep^3*n2^2*x5^3+1510*ep^3
         *n2^2*n5+1582*ep^3*n2^2*n5*x5+198*ep^3*n2^2*n5*x5^2-40*ep^3*n2^2*n5*
         x5^3-462*ep^3*n2^2*n5^2-321*ep^3*n2^2*n5^2*x5-12*ep^3*n2^2*n5^2*x5^2+
         10*ep^3*n2^2*n5^2*x5^3+426*ep^3*n2^3*x5+108*ep^3*n2^3*x5^2-7*ep^3*n2^3
         *x5^3-242*ep^3*n2^3*n5-264*ep^3*n2^3*n5*x5-39*ep^3*n2^3*n5*x5^2+8*ep^3
         *n2^3*n5*x5^3-50*ep^3*n2^4*x5-14*ep^3*n2^4*x5^2+2054*ep^4*x5+1056*ep^4
         *x5^2+113*ep^4*x5^3-2*ep^4*x5^4-2318*ep^4*n5-3658*ep^4*n5*x5-1299*ep^4
         *n5*x5^2-116*ep^4*n5*x5^3+1376*ep^4*n5^2+1494*ep^4*n5^2*x5+420*ep^4*
         n5^2*x5^2+30*ep^4*n5^2*x5^3-186*ep^4*n5^3-156*ep^4*n5^3*x5-30*ep^4*
         n5^3*x5^2-3080*ep^4*n2*x5-1596*ep^4*n2*x5^2-153*ep^4*n2*x5^3+10*ep^4*
         n2*x5^4+2318*ep^4*n2*n5+3720*ep^4*n2*n5*x5+1329*ep^4*n2*n5*x5^2+100*
         ep^4*n2*n5*x5^3-5*ep^4*n2*n5*x5^4-702*ep^4*n2*n5^2-786*ep^4*n2*n5^2*x5
         -234*ep^4*n2*n5^2*x5^2-20*ep^4*n2*n5^2*x5^3+1510*ep^4*n2^2*x5+791*ep^4
         *n2^2*x5^2+66*ep^4*n2^2*x5^3-10*ep^4*n2^2*x5^4-568*ep^4*n2^2*n5-924*
         ep^4*n2^2*n5*x5-321*ep^4*n2^2*n5*x5^2-8*ep^4*n2^2*n5*x5^3+5*ep^4*n2^2*
         n5*x5^4-242*ep^4*n2^3*x5-132*ep^4*n2^3*x5^2-13*ep^4*n2^3*x5^3+2*ep^4*
         n2^3*x5^4-2318*ep^5*x5-1829*ep^5*x5^2-433*ep^5*x5^3-29*ep^5*x5^4+1296*
         ep^5*n5+2752*ep^5*n5*x5+1494*ep^5*n5*x5^2+280*ep^5*n5*x5^3+15*ep^5*n5*
         x5^4-384*ep^5*n5^2-558*ep^5*n5^2*x5-234*ep^5*n5^2*x5^2-30*ep^5*n5^2*
         x5^3+2318*ep^5*n2*x5+1860*ep^5*n2*x5^2+443*ep^5*n2*x5^3+25*ep^5*n2*
         x5^4-ep^5*n2*x5^5-648*ep^5*n2*n5-1404*ep^5*n2*n5*x5-786*ep^5*n2*n5*
         x5^2-156*ep^5*n2*n5*x5^3-10*ep^5*n2*n5*x5^4-568*ep^5*n2^2*x5-462*ep^5*
         n2^2*x5^2-107*ep^5*n2^2*x5^3-2*ep^5*n2^2*x5^4+ep^5*n2^2*x5^5+1296*ep^6
         *x5+1376*ep^6*x5^2+498*ep^6*x5^3+70*ep^6*x5^4+3*ep^6*x5^5-288*ep^6*n5-
         768*ep^6*n5*x5-558*ep^6*n5*x5^2-156*ep^6*n5*x5^3-15*ep^6*n5*x5^4-648*
         ep^6*n2*x5-702*ep^6*n2*x5^2-262*ep^6*n2*x5^3-39*ep^6*n2*x5^4-2*ep^6*n2
         *x5^5-288*ep^7*x5-384*ep^7*x5^2-186*ep^7*x5^3-39*ep^7*x5^4-3*ep^7*x5^5
         ,-16+16*n5-4*n5^2+84*n2-78*n2*n5+18*n2*n5^2-184*n2^2+155*n2^2*n5-32*
         n2^2*n5^2+217*n2^3-160*n2^3*n5+28*n2^3*n5^2-148*n2^4+90*n2^4*n5-12*
         n2^4*n5^2+58*n2^5-26*n2^5*n5+2*n2^5*n5^2-12*n2^6+3*n2^6*n5+n2^7+120*ep
         +16*ep*x5-104*ep*n5-8*ep*n5*x5+22*ep*n5^2-542*ep*n2-78*ep*n2*x5+425*ep
         *n2*n5+36*ep*n2*n5*x5-80*ep*n2*n5^2+993*ep*n2^2+155*ep*n2^2*x5-680*ep*
         n2^2*n5-64*ep*n2^2*n5*x5+108*ep*n2^2*n5^2-940*ep*n2^3-160*ep*n2^3*x5+
         530*ep*n2^3*n5+56*ep*n2^3*n5*x5-64*ep*n2^3*n5^2+482*ep*n2^4+90*ep*n2^4
         *x5-200*ep*n2^4*n5-24*ep*n2^4*n5*x5+14*ep*n2^4*n5^2-126*ep*n2^5-26*ep*
         n2^5*x5+29*ep*n2^5*n5+4*ep*n2^5*n5*x5+13*ep*n2^6+3*ep*n2^6*x5-368*ep^2
         -104*ep^2*x5-4*ep^2*x5^2+264*ep^2*n5+44*ep^2*n5*x5-44*ep^2*n5^2+1398*
         ep^2*n2+425*ep^2*n2*x5+18*ep^2*n2*x5^2-874*ep^2*n2*n5-160*ep^2*n2*n5*
         x5+122*ep^2*n2*n5^2-2068*ep^2*n2^2-680*ep^2*n2^2*x5-32*ep^2*n2^2*x5^2+
         1063*ep^2*n2^2*n5+216*ep^2*n2^2*n5*x5-112*ep^2*n2^2*n5^2+1481*ep^2*
         n2^3+530*ep^2*n2^3*x5+28*ep^2*n2^3*x5^2-560*ep^2*n2^3*n5-128*ep^2*n2^3
         *n5*x5+34*ep^2*n2^3*n5^2-510*ep^2*n2^4-200*ep^2*n2^4*x5-12*ep^2*n2^4*
         x5^2+107*ep^2*n2^4*n5+28*ep^2*n2^4*n5*x5+67*ep^2*n2^5+29*ep^2*n2^5*x5+
         2*ep^2*n2^5*x5^2+592*ep^3+264*ep^3*x5+22*ep^3*x5^2-328*ep^3*n5-88*ep^3
         *n5*x5+38*ep^3*n5^2-1830*ep^3*n2-874*ep^3*n2*x5-80*ep^3*n2*x5^2+833*
         ep^3*n2*n5+244*ep^3*n2*n5*x5-72*ep^3*n2*n5^2+2065*ep^3*n2^2+1063*ep^3*
         n2^2*x5+108*ep^3*n2^2*x5^2-692*ep^3*n2^2*n5-224*ep^3*n2^2*n5*x5+34*
         ep^3*n2^2*n5^2-1002*ep^3*n2^3-560*ep^3*n2^3*x5-64*ep^3*n2^3*x5^2+187*
         ep^3*n2^3*n5+68*ep^3*n2^3*n5*x5+175*ep^3*n2^4+107*ep^3*n2^4*x5+14*ep^3
         *n2^4*x5^2-528*ep^4-328*ep^4*x5-44*ep^4*x5^2+200*ep^4*n5+76*ep^4*n5*x5
         -12*ep^4*n5^2+1262*ep^4*n2+833*ep^4*n2*x5+122*ep^4*n2*x5^2-354*ep^4*n2
         *n5-144*ep^4*n2*n5*x5+12*ep^4*n2*n5^2-978*ep^4*n2^2-692*ep^4*n2^2*x5-
         112*ep^4*n2^2*x5^2+154*ep^4*n2^2*n5+68*ep^4*n2^2*n5*x5+244*ep^4*n2^3+
         187*ep^4*n2^3*x5+34*ep^4*n2^3*x5^2+248*ep^5+200*ep^5*x5+38*ep^5*x5^2-
         48*ep^5*n5-24*ep^5*n5*x5-420*ep^5*n2-354*ep^5*n2*x5-72*ep^5*n2*x5^2+48
         *ep^5*n2*n5+24*ep^5*n2*n5*x5+172*ep^5*n2^2+154*ep^5*n2^2*x5+34*ep^5*
         n2^2*x5^2-48*ep^6-48*ep^6*x5-12*ep^6*x5^2+48*ep^6*n2+48*ep^6*n2*x5+12*
         ep^6*n2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,-1+n2,1,1,1+n5,0,0,1,0)*(rat(-18*n5+39*n5^2-29*
         n5^3+9*n5^4-n5^5+18*n2*n5-33*n2*n5^2+18*n2*n5^3-3*n2*n5^4-4*n2^2*n5+6*
         n2^2*n5^2-2*n2^2*n5^3-18*ep*x5+132*ep*n5+78*ep*n5*x5-191*ep*n5^2-87*ep
         *n5^2*x5+88*ep*n5^3+36*ep*n5^3*x5-13*ep*n5^4-5*ep*n5^4*x5+18*ep*n2*x5-
         128*ep*n2*n5-66*ep*n2*n5*x5+126*ep*n2*n5^2+54*ep*n2*n5^2*x5-30*ep*n2*
         n5^3-12*ep*n2*n5^3*x5-4*ep*n2^2*x5+38*ep*n2^2*n5+12*ep*n2^2*n5*x5-18*
         ep*n2^2*n5^2-6*ep*n2^2*n5^2*x5-4*ep*n2^3*n5+132*ep^2*x5+39*ep^2*x5^2-
         306*ep^2*n5-382*ep^2*n5*x5-87*ep^2*n5*x5^2+280*ep^2*n5^2+264*ep^2*n5^2
         *x5+54*ep^2*n5^2*x5^2-62*ep^2*n5^3-52*ep^2*n5^3*x5-10*ep^2*n5^3*x5^2-
         128*ep^2*n2*x5-33*ep^2*n2*x5^2+214*ep^2*n2*n5+252*ep^2*n2*n5*x5+54*
         ep^2*n2*n5*x5^2-98*ep^2*n2*n5^2-90*ep^2*n2*n5^2*x5-18*ep^2*n2*n5^2*
         x5^2+38*ep^2*n2^2*x5+6*ep^2*n2^2*x5^2-36*ep^2*n2^2*n5-36*ep^2*n2^2*n5*
         x5-6*ep^2*n2^2*n5*x5^2-4*ep^2*n2^3*x5-306*ep^3*x5-191*ep^3*x5^2-29*
         ep^3*x5^3+288*ep^3*n5+560*ep^3*n5*x5+264*ep^3*n5*x5^2+36*ep^3*n5*x5^3-
         128*ep^3*n5^2-186*ep^3*n5^2*x5-78*ep^3*n5^2*x5^2-10*ep^3*n5^2*x5^3+214
         *ep^3*n2*x5+126*ep^3*n2*x5^2+18*ep^3*n2*x5^3-104*ep^3*n2*n5-196*ep^3*
         n2*n5*x5-90*ep^3*n2*n5*x5^2-12*ep^3*n2*n5*x5^3-36*ep^3*n2^2*x5-18*ep^3
         *n2^2*x5^2-2*ep^3*n2^2*x5^3+288*ep^4*x5+280*ep^4*x5^2+88*ep^4*x5^3+9*
         ep^4*x5^4-96*ep^4*n5-256*ep^4*n5*x5-186*ep^4*n5*x5^2-52*ep^4*n5*x5^3-5
         *ep^4*n5*x5^4-104*ep^4*n2*x5-98*ep^4*n2*x5^2-30*ep^4*n2*x5^3-3*ep^4*n2
         *x5^4-96*ep^5*x5-128*ep^5*x5^2-62*ep^5*x5^3-13*ep^5*x5^4-ep^5*x5^5,-8+
         8*n5-2*n5^2+30*n2-27*n2*n5+6*n2*n5^2-43*n2^2+33*n2^2*n5-6*n2^2*n5^2+29
         *n2^3-17*n2^3*n5+2*n2^3*n5^2-9*n2^4+3*n2^4*n5+n2^5+40*ep+8*ep*x5-32*ep
         *n5-4*ep*n5*x5+6*ep*n5^2-120*ep*n2-27*ep*n2*x5+81*ep*n2*n5+12*ep*n2*n5
         *x5-12*ep*n2*n5^2+129*ep*n2^2+33*ep*n2^2*x5-66*ep*n2^2*n5-12*ep*n2^2*
         n5*x5+6*ep*n2^2*n5^2-58*ep*n2^3-17*ep*n2^3*x5+17*ep*n2^3*n5+4*ep*n2^3*
         n5*x5+9*ep*n2^4+3*ep*n2^4*x5-72*ep^2-32*ep^2*x5-2*ep^2*x5^2+40*ep^2*n5
         +12*ep^2*n5*x5-4*ep^2*n5^2+166*ep^2*n2+81*ep^2*n2*x5+6*ep^2*n2*x5^2-70
         *ep^2*n2*n5-24*ep^2*n2*n5*x5+4*ep^2*n2*n5^2-122*ep^2*n2^2-66*ep^2*n2^2
         *x5-6*ep^2*n2^2*x5^2+30*ep^2*n2^2*n5+12*ep^2*n2^2*n5*x5+28*ep^2*n2^3+
         17*ep^2*n2^3*x5+2*ep^2*n2^3*x5^2+56*ep^3+40*ep^3*x5+6*ep^3*x5^2-16*
         ep^3*n5-8*ep^3*n5*x5-92*ep^3*n2-70*ep^3*n2*x5-12*ep^3*n2*x5^2+16*ep^3*
         n2*n5+8*ep^3*n2*n5*x5+36*ep^3*n2^2+30*ep^3*n2^2*x5+6*ep^3*n2^2*x5^2-16
         *ep^4-16*ep^4*x5-4*ep^4*x5^2+16*ep^4*n2+16*ep^4*n2*x5+4*ep^4*n2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,-1+n2,1,1,1+n5,1,0,0,0)*(rat(48*n5-4*n5^2-70*
         n5^3+58*n5^4-18*n5^5+2*n5^6-250*n2*n5+143*n2*n5^2+101*n2*n5^3-97*n2*
         n5^4+25*n2*n5^5-2*n2*n5^6+495*n2^2*n5-357*n2^2*n5^2-11*n2^2*n5^3+49*
         n2^2*n5^4-8*n2^2*n5^5-499*n2^3*n5+353*n2^3*n5^2-36*n2^3*n5^3-11*n2^3*
         n5^4+n2^3*n5^5+282*n2^4*n5-173*n2^4*n5^2+19*n2^4*n5^3+n2^4*n5^4-90*
         n2^5*n5+42*n2^5*n5^2-3*n2^5*n5^3+15*n2^6*n5-4*n2^6*n5^2-n2^7*n5+48*ep*
         x5+156*ep*n5-8*ep*n5*x5-950*ep*n5^2-210*ep*n5^2*x5+1118*ep*n5^3+232*ep
         *n5^3*x5-530*ep*n5^4-90*ep*n5^4*x5+110*ep*n5^5+12*ep*n5^5*x5-8*ep*n5^6
         -250*ep*n2*x5+168*ep*n2*n5+286*ep*n2*n5*x5+1549*ep*n2*n5^2+303*ep*n2*
         n5^2*x5-1917*ep*n2*n5^3-388*ep*n2*n5^3*x5+777*ep*n2*n5^4+125*ep*n2*
         n5^4*x5-123*ep*n2*n5^5-12*ep*n2*n5^5*x5+6*ep*n2*n5^6+495*ep*n2^2*x5-
         1032*ep*n2^2*n5-714*ep*n2^2*n5*x5-558*ep*n2^2*n5^2-33*ep*n2^2*n5^2*x5+
         1065*ep*n2^2*n5^3+196*ep*n2^2*n5^3*x5-332*ep*n2^2*n5^4-40*ep*n2^2*n5^4
         *x5+27*ep*n2^2*n5^5-499*ep*n2^3*x5+1186*ep*n2^3*n5+706*ep*n2^3*n5*x5-
         228*ep*n2^3*n5^2-108*ep*n2^3*n5^2*x5-223*ep*n2^3*n5^3-44*ep*n2^3*n5^3*
         x5+46*ep*n2^3*n5^4+5*ep*n2^3*n5^4*x5+282*ep*n2^4*x5-594*ep*n2^4*n5-346
         *ep*n2^4*n5*x5+179*ep*n2^4*n5^2+57*ep*n2^4*n5^2*x5+11*ep*n2^4*n5^3+4*
         ep*n2^4*n5^3*x5-90*ep*n2^5*x5+138*ep*n2^5*n5+84*ep*n2^5*n5*x5-29*ep*
         n2^5*n5^2-9*ep*n2^5*n5^2*x5+15*ep*n2^6*x5-12*ep*n2^6*n5-8*ep*n2^6*n5*
         x5-ep*n2^7*x5+156*ep^2*x5-4*ep^2*x5^2-2368*ep^2*n5-1900*ep^2*n5*x5-210
         *ep^2*n5*x5^2+5458*ep^2*n5^2+3354*ep^2*n5^2*x5+348*ep^2*n5^2*x5^2-4150
         *ep^2*n5^3-2120*ep^2*n5^3*x5-180*ep^2*n5^3*x5^2+1342*ep^2*n5^4+550*
         ep^2*n5^4*x5+30*ep^2*n5^4*x5^2-176*ep^2*n5^5-48*ep^2*n5^5*x5+6*ep^2*
         n5^6+168*ep^2*n2*x5+143*ep^2*n2*x5^2+3900*ep^2*n2*n5+3098*ep^2*n2*n5*
         x5+303*ep^2*n2*n5*x5^2-8879*ep^2*n2*n5^2-5751*ep^2*n2*n5^2*x5-582*ep^2
         *n2*n5^2*x5^2+5729*ep^2*n2*n5^3+3108*ep^2*n2*n5^3*x5+250*ep^2*n2*n5^3*
         x5^2-1380*ep^2*n2*n5^4-615*ep^2*n2*n5^4*x5-30*ep^2*n2*n5^4*x5^2+106*
         ep^2*n2*n5^5+36*ep^2*n2*n5^5*x5-1032*ep^2*n2^2*x5-357*ep^2*n2^2*x5^2-
         1713*ep^2*n2^2*n5-1116*ep^2*n2^2*n5*x5-33*ep^2*n2^2*n5*x5^2+4680*ep^2*
         n2^2*n5^2+3195*ep^2*n2^2*n5^2*x5+294*ep^2*n2^2*n5^2*x5^2-2376*ep^2*
         n2^2*n5^3-1328*ep^2*n2^2*n5^3*x5-80*ep^2*n2^2*n5^3*x5^2+313*ep^2*n2^2*
         n5^4+135*ep^2*n2^2*n5^4*x5+1186*ep^2*n2^3*x5+353*ep^2*n2^3*x5^2-175*
         ep^2*n2^3*n5-456*ep^2*n2^3*n5*x5-108*ep^2*n2^3*n5*x5^2-857*ep^2*n2^3*
         n5^2-669*ep^2*n2^3*n5^2*x5-66*ep^2*n2^3*n5^2*x5^2+301*ep^2*n2^3*n5^3+
         184*ep^2*n2^3*n5^3*x5+10*ep^2*n2^3*n5^3*x5^2-594*ep^2*n2^4*x5-173*ep^2
         *n2^4*x5^2+259*ep^2*n2^4*n5+358*ep^2*n2^4*n5*x5+57*ep^2*n2^4*n5*x5^2+
         21*ep^2*n2^4*n5^2+33*ep^2*n2^4*n5^2*x5+6*ep^2*n2^4*n5^2*x5^2+138*ep^2*
         n2^5*x5+42*ep^2*n2^5*x5^2-41*ep^2*n2^5*n5-58*ep^2*n2^5*n5*x5-9*ep^2*
         n2^5*n5*x5^2-12*ep^2*n2^6*x5-4*ep^2*n2^6*x5^2-2368*ep^3*x5-950*ep^3*
         x5^2-70*ep^3*x5^3+7912*ep^3*n5+10916*ep^3*n5*x5+3354*ep^3*n5*x5^2+232*
         ep^3*n5*x5^3-12362*ep^3*n5^2-12450*ep^3*n5^2*x5-3180*ep^3*n5^2*x5^2-
         180*ep^3*n5^2*x5^3+6486*ep^3*n5^3+5368*ep^3*n5^3*x5+1100*ep^3*n5^3*
         x5^2+40*ep^3*n5^3*x5^3-1332*ep^3*n5^4-880*ep^3*n5^4*x5-120*ep^3*n5^4*
         x5^2+84*ep^3*n5^5+36*ep^3*n5^5*x5+3900*ep^3*n2*x5+1549*ep^3*n2*x5^2+
         101*ep^3*n2*x5^3-12448*ep^3*n2*n5-17758*ep^3*n2*n5*x5-5751*ep^3*n2*n5*
         x5^2-388*ep^3*n2*n5*x5^3+16071*ep^3*n2*n5^2+17187*ep^3*n2*n5^2*x5+4662
         *ep^3*n2*n5^2*x5^2+250*ep^3*n2*n5^2*x5^3-6157*ep^3*n2*n5^3-5520*ep^3*
         n2*n5^3*x5-1230*ep^3*n2*n5^3*x5^2-40*ep^3*n2*n5^3*x5^3+704*ep^3*n2*
         n5^4+530*ep^3*n2*n5^4*x5+90*ep^3*n2*n5^4*x5^2-1713*ep^3*n2^2*x5-558*
         ep^3*n2^2*x5^2-11*ep^3*n2^2*x5^3+6490*ep^3*n2^2*n5+9360*ep^3*n2^2*n5*
         x5+3195*ep^3*n2^2*n5*x5^2+196*ep^3*n2^2*n5*x5^3-6329*ep^3*n2^2*n5^2-
         7128*ep^3*n2^2*n5^2*x5-1992*ep^3*n2^2*n5^2*x5^2-80*ep^3*n2^2*n5^2*x5^3
         +1333*ep^3*n2^2*n5^3+1252*ep^3*n2^2*n5^3*x5+270*ep^3*n2^2*n5^3*x5^2-
         175*ep^3*n2^3*x5-228*ep^3*n2^3*x5^2-36*ep^3*n2^3*x5^3-1242*ep^3*n2^3*
         n5-1714*ep^3*n2^3*n5*x5-669*ep^3*n2^3*n5*x5^2-44*ep^3*n2^3*n5*x5^3+742
         *ep^3*n2^3*n5^2+903*ep^3*n2^3*n5^2*x5+276*ep^3*n2^3*n5^2*x5^2+10*ep^3*
         n2^3*n5^2*x5^3+259*ep^3*n2^4*x5+179*ep^3*n2^4*x5^2+19*ep^3*n2^4*x5^3+
         56*ep^3*n2^4*n5+42*ep^3*n2^4*n5*x5+33*ep^3*n2^4*n5*x5^2+4*ep^3*n2^4*n5
         *x5^3-41*ep^3*n2^5*x5-29*ep^3*n2^5*x5^2-3*ep^3*n2^5*x5^3+7912*ep^4*x5+
         5458*ep^4*x5^2+1118*ep^4*x5^3+58*ep^4*x5^4-12768*ep^4*n5-24724*ep^4*n5
         *x5-12450*ep^4*n5*x5^2-2120*ep^4*n5*x5^3-90*ep^4*n5*x5^4+13898*ep^4*
         n5^2+19458*ep^4*n5^2*x5+8052*ep^4*n5^2*x5^2+1100*ep^4*n5^2*x5^3+30*
         ep^4*n5^2*x5^4-4644*ep^4*n5^3-5328*ep^4*n5^3*x5-1760*ep^4*n5^3*x5^2-
         160*ep^4*n5^3*x5^3+462*ep^4*n5^4+420*ep^4*n5^4*x5+90*ep^4*n5^4*x5^2-
         12448*ep^4*n2*x5-8879*ep^4*n2*x5^2-1917*ep^4*n2*x5^3-97*ep^4*n2*x5^4+
         16014*ep^4*n2*n5+32142*ep^4*n2*n5*x5+17187*ep^4*n2*n5*x5^2+3108*ep^4*
         n2*n5*x5^3+125*ep^4*n2*n5*x5^4-12416*ep^4*n2*n5^2-18471*ep^4*n2*n5^2*
         x5-8280*ep^4*n2*n5^2*x5^2-1230*ep^4*n2*n5^2*x5^3-30*ep^4*n2*n5^2*x5^4+
         2244*ep^4*n2*n5^3+2816*ep^4*n2*n5^3*x5+1060*ep^4*n2*n5^3*x5^2+120*ep^4
         *n2*n5^3*x5^3+6490*ep^4*n2^2*x5+4680*ep^4*n2^2*x5^2+1065*ep^4*n2^2*
         x5^3+49*ep^4*n2^2*x5^4-6188*ep^4*n2^2*n5-12658*ep^4*n2^2*n5*x5-7128*
         ep^4*n2^2*n5*x5^2-1328*ep^4*n2^2*n5*x5^3-40*ep^4*n2^2*n5*x5^4+2564*
         ep^4*n2^2*n5^2+3999*ep^4*n2^2*n5^2*x5+1878*ep^4*n2^2*n5^2*x5^2+270*
         ep^4*n2^2*n5^2*x5^3-1242*ep^4*n2^3*x5-857*ep^4*n2^3*x5^2-223*ep^4*n2^3
         *x5^3-11*ep^4*n2^3*x5^4+730*ep^4*n2^3*n5+1484*ep^4*n2^3*n5*x5+903*ep^4
         *n2^3*n5*x5^2+184*ep^4*n2^3*n5*x5^3+5*ep^4*n2^3*n5*x5^4+56*ep^4*n2^4*
         x5+21*ep^4*n2^4*x5^2+11*ep^4*n2^4*x5^3+ep^4*n2^4*x5^4-12768*ep^5*x5-
         12362*ep^5*x5^2-4150*ep^5*x5^3-530*ep^5*x5^4-18*ep^5*x5^5+11164*ep^5*
         n5+27796*ep^5*n5*x5+19458*ep^5*n5*x5^2+5368*ep^5*n5*x5^3+550*ep^5*n5*
         x5^4+12*ep^5*n5*x5^5-7768*ep^5*n5^2-13932*ep^5*n5^2*x5-7992*ep^5*n5^2*
         x5^2-1760*ep^5*n5^2*x5^3-120*ep^5*n5^2*x5^4+1260*ep^5*n5^3+1848*ep^5*
         n5^3*x5+840*ep^5*n5^3*x5^2+120*ep^5*n5^3*x5^3+16014*ep^5*n2*x5+16071*
         ep^5*n2*x5^2+5729*ep^5*n2*x5^3+777*ep^5*n2*x5^4+25*ep^5*n2*x5^5-9624*
         ep^5*n2*n5-24832*ep^5*n2*n5*x5-18471*ep^5*n2*n5*x5^2-5520*ep^5*n2*n5*
         x5^3-615*ep^5*n2*n5*x5^4-12*ep^5*n2*n5*x5^5+3532*ep^5*n2*n5^2+6732*
         ep^5*n2*n5^2*x5+4224*ep^5*n2*n5^2*x5^2+1060*ep^5*n2*n5^2*x5^3+90*ep^5*
         n2*n5^2*x5^4-6188*ep^5*n2^2*x5-6329*ep^5*n2^2*x5^2-2376*ep^5*n2^2*x5^3
         -332*ep^5*n2^2*x5^4-8*ep^5*n2^2*x5^5+1948*ep^5*n2^2*n5+5128*ep^5*n2^2*
         n5*x5+3999*ep^5*n2^2*n5*x5^2+1252*ep^5*n2^2*n5*x5^3+135*ep^5*n2^2*n5*
         x5^4+730*ep^5*n2^3*x5+742*ep^5*n2^3*x5^2+301*ep^5*n2^3*x5^3+46*ep^5*
         n2^3*x5^4+ep^5*n2^3*x5^5+11164*ep^6*x5+13898*ep^6*x5^2+6486*ep^6*x5^3+
         1342*ep^6*x5^4+110*ep^6*x5^5+2*ep^6*x5^6-5104*ep^6*n5-15536*ep^6*n5*x5
         -13932*ep^6*n5*x5^2-5328*ep^6*n5*x5^3-880*ep^6*n5*x5^4-48*ep^6*n5*x5^5
         +1728*ep^6*n5^2+3780*ep^6*n5^2*x5+2772*ep^6*n5^2*x5^2+840*ep^6*n5^2*
         x5^3+90*ep^6*n5^2*x5^4-9624*ep^6*n2*x5-12416*ep^6*n2*x5^2-6157*ep^6*n2
         *x5^3-1380*ep^6*n2*x5^4-123*ep^6*n2*x5^5-2*ep^6*n2*x5^6+2240*ep^6*n2*
         n5+7064*ep^6*n2*n5*x5+6732*ep^6*n2*n5*x5^2+2816*ep^6*n2*n5*x5^3+530*
         ep^6*n2*n5*x5^4+36*ep^6*n2*n5*x5^5+1948*ep^6*n2^2*x5+2564*ep^6*n2^2*
         x5^2+1333*ep^6*n2^2*x5^3+313*ep^6*n2^2*x5^4+27*ep^6*n2^2*x5^5-5104*
         ep^7*x5-7768*ep^7*x5^2-4644*ep^7*x5^3-1332*ep^7*x5^4-176*ep^7*x5^5-8*
         ep^7*x5^6+960*ep^7*n5+3456*ep^7*n5*x5+3780*ep^7*n5*x5^2+1848*ep^7*n5*
         x5^3+420*ep^7*n5*x5^4+36*ep^7*n5*x5^5+2240*ep^7*n2*x5+3532*ep^7*n2*
         x5^2+2244*ep^7*n2*x5^3+704*ep^7*n2*x5^4+106*ep^7*n2*x5^5+6*ep^7*n2*
         x5^6+960*ep^8*x5+1728*ep^8*x5^2+1260*ep^8*x5^3+462*ep^8*x5^4+84*ep^8*
         x5^5+6*ep^8*x5^6,48-64*n5+28*n5^2-4*n5^3-268*n2+334*n2*n5-136*n2*n5^2+
         18*n2*n5^3+636*n2^2-727*n2^2*n5+269*n2^2*n5^2-32*n2^2*n5^3-835*n2^3+
         852*n2^3*n5-276*n2^3*n5^2+28*n2^3*n5^3+661*n2^4-578*n2^4*n5+154*n2^4*
         n5^2-12*n2^4*n5^3-322*n2^5+226*n2^5*n5-44*n2^5*n5^2+2*n2^5*n5^3+94*
         n2^6-47*n2^6*n5+5*n2^6*n5^2-15*n2^7+4*n2^7*n5+n2^8-408*ep-64*ep*x5+480
         *ep*n5+56*ep*n5*x5-182*ep*n5^2-12*ep*n5^2*x5+22*ep*n5^3+1998*ep*n2+334
         *ep*n2*x5-2155*ep*n2*n5-272*ep*n2*n5*x5+741*ep*n2*n5^2+54*ep*n2*n5^2*
         x5-80*ep*n2*n5^3-4073*ep*n2^2-727*ep*n2^2*x5+3923*ep*n2^2*n5+538*ep*
         n2^2*n5*x5-1180*ep*n2^2*n5^2-96*ep*n2^2*n5^2*x5+108*ep*n2^2*n5^3+4464*
         ep*n2^3+852*ep*n2^3*x5-3690*ep*n2^3*n5-552*ep*n2^3*n5*x5+914*ep*n2^3*
         n5^2+84*ep*n2^3*n5^2*x5-64*ep*n2^3*n5^3-2830*ep*n2^4-578*ep*n2^4*x5+
         1882*ep*n2^4*n5+308*ep*n2^4*n5*x5-342*ep*n2^4*n5^2-36*ep*n2^4*n5^2*x5+
         14*ep*n2^4*n5^3+1034*ep*n2^5+226*ep*n2^5*x5-491*ep*n2^5*n5-88*ep*n2^5*
         n5*x5+49*ep*n2^5*n5^2+6*ep*n2^5*n5^2*x5-201*ep*n2^6-47*ep*n2^6*x5+51*
         ep*n2^6*n5+10*ep*n2^6*n5*x5+16*ep*n2^7+4*ep*n2^7*x5+1464*ep^2+480*ep^2
         *x5+28*ep^2*x5^2-1472*ep^2*n5-364*ep^2*n5*x5-12*ep^2*n5*x5^2+462*ep^2*
         n5^2+66*ep^2*n5^2*x5-44*ep^2*n5^3-6188*ep^2*n2-2155*ep^2*n2*x5-136*
         ep^2*n2*x5^2+5559*ep^2*n2*n5+1482*ep^2*n2*n5*x5+54*ep^2*n2*n5*x5^2-
         1524*ep^2*n2*n5^2-240*ep^2*n2*n5^2*x5+122*ep^2*n2*n5^3+10581*ep^2*n2^2
         +3923*ep^2*n2^2*x5+269*ep^2*n2^2*x5^2-8171*ep^2*n2^2*n5-2360*ep^2*n2^2
         *n5*x5-96*ep^2*n2^2*n5*x5^2+1845*ep^2*n2^2*n5^2+324*ep^2*n2^2*n5^2*x5-
         112*ep^2*n2^2*n5^3-9331*ep^2*n2^3-3690*ep^2*n2^3*x5-276*ep^2*n2^3*x5^2
         +5814*ep^2*n2^3*n5+1828*ep^2*n2^3*n5*x5+84*ep^2*n2^3*n5*x5^2-966*ep^2*
         n2^3*n5^2-192*ep^2*n2^3*n5^2*x5+34*ep^2*n2^3*n5^3+4457*ep^2*n2^4+1882*
         ep^2*n2^4*x5+154*ep^2*n2^4*x5^2-1991*ep^2*n2^4*n5-684*ep^2*n2^4*n5*x5-
         36*ep^2*n2^4*n5*x5^2+183*ep^2*n2^4*n5^2+42*ep^2*n2^4*n5^2*x5-1089*ep^2
         *n2^5-491*ep^2*n2^5*x5-44*ep^2*n2^5*x5^2+261*ep^2*n2^5*n5+98*ep^2*n2^5
         *n5*x5+6*ep^2*n2^5*n5*x5^2+106*ep^2*n2^6+51*ep^2*n2^6*x5+5*ep^2*n2^6*
         x5^2-2880*ep^3-1472*ep^3*x5-182*ep^3*x5^2-4*ep^3*x5^3+2368*ep^3*n5+924
         *ep^3*n5*x5+66*ep^3*n5*x5^2-574*ep^3*n5^2-132*ep^3*n5^2*x5+38*ep^3*
         n5^3+10276*ep^3*n2+5559*ep^3*n2*x5+741*ep^3*n2*x5^2+18*ep^3*n2*x5^3-
         7279*ep^3*n2*n5-3048*ep^3*n2*n5*x5-240*ep^3*n2*n5*x5^2+1453*ep^3*n2*
         n5^2+366*ep^3*n2*n5^2*x5-72*ep^3*n2*n5^3-14229*ep^3*n2^2-8171*ep^3*
         n2^2*x5-1180*ep^3*n2^2*x5^2-32*ep^3*n2^2*x5^3+8163*ep^3*n2^2*n5+3690*
         ep^3*n2^2*n5*x5+324*ep^3*n2^2*n5*x5^2-1202*ep^3*n2^2*n5^2-336*ep^3*
         n2^2*n5^2*x5+34*ep^3*n2^2*n5^3+9514*ep^3*n2^3+5814*ep^3*n2^3*x5+914*
         ep^3*n2^3*x5^2+28*ep^3*n2^3*x5^3-3935*ep^3*n2^3*n5-1932*ep^3*n2^3*n5*
         x5-192*ep^3*n2^3*n5*x5^2+323*ep^3*n2^3*n5^2+102*ep^3*n2^3*n5^2*x5-3057
         *ep^3*n2^4-1991*ep^3*n2^4*x5-342*ep^3*n2^4*x5^2-12*ep^3*n2^4*x5^3+683*
         ep^3*n2^4*n5+366*ep^3*n2^4*n5*x5+42*ep^3*n2^4*n5*x5^2+376*ep^3*n2^5+
         261*ep^3*n2^5*x5+49*ep^3*n2^5*x5^2+2*ep^3*n2^5*x5^3+3360*ep^4+2368*
         ep^4*x5+462*ep^4*x5^2+22*ep^4*x5^3-2112*ep^4*n5-1148*ep^4*n5*x5-132*
         ep^4*n5*x5^2+350*ep^4*n5^2+114*ep^4*n5^2*x5-12*ep^4*n5^3-9804*ep^4*n2-
         7279*ep^4*n2*x5-1524*ep^4*n2*x5^2-80*ep^4*n2*x5^3+5023*ep^4*n2*n5+2906
         *ep^4*n2*n5*x5+366*ep^4*n2*n5*x5^2-618*ep^4*n2*n5^2-216*ep^4*n2*n5^2*
         x5+12*ep^4*n2*n5^3+10391*ep^4*n2^2+8163*ep^4*n2^2*x5+1845*ep^4*n2^2*
         x5^2+108*ep^4*n2^2*x5^3-3870*ep^4*n2^2*n5-2404*ep^4*n2^2*n5*x5-336*
         ep^4*n2^2*n5*x5^2+268*ep^4*n2^2*n5^2+102*ep^4*n2^2*n5^2*x5-4716*ep^4*
         n2^3-3935*ep^4*n2^3*x5-966*ep^4*n2^3*x5^2-64*ep^4*n2^3*x5^3+959*ep^4*
         n2^3*n5+646*ep^4*n2^3*n5*x5+102*ep^4*n2^3*n5*x5^2+769*ep^4*n2^4+683*
         ep^4*n2^4*x5+183*ep^4*n2^4*x5^2+14*ep^4*n2^4*x5^3-2328*ep^5-2112*ep^5*
         x5-574*ep^5*x5^2-44*ep^5*x5^3+992*ep^5*n5+700*ep^5*n5*x5+114*ep^5*n5*
         x5^2-84*ep^5*n5^2-36*ep^5*n5^2*x5+5294*ep^5*n2+5023*ep^5*n2*x5+1453*
         ep^5*n2*x5^2+122*ep^5*n2*x5^3-1674*ep^5*n2*n5-1236*ep^5*n2*n5*x5-216*
         ep^5*n2*n5*x5^2+84*ep^5*n2*n5^2+36*ep^5*n2*n5^2*x5-3870*ep^5*n2^2-3870
         *ep^5*n2^2*x5-1202*ep^5*n2^2*x5^2-112*ep^5*n2^2*x5^3+682*ep^5*n2^2*n5+
         536*ep^5*n2^2*n5*x5+102*ep^5*n2^2*n5*x5^2+904*ep^5*n2^3+959*ep^5*n2^3*
         x5+323*ep^5*n2^3*x5^2+34*ep^5*n2^3*x5^3+888*ep^6+992*ep^6*x5+350*ep^6*
         x5^2+38*ep^6*x5^3-192*ep^6*n5-168*ep^6*n5*x5-36*ep^6*n5*x5^2-1452*ep^6
         *n2-1674*ep^6*n2*x5-618*ep^6*n2*x5^2-72*ep^6*n2*x5^3+192*ep^6*n2*n5+
         168*ep^6*n2*n5*x5+36*ep^6*n2*n5*x5^2+564*ep^6*n2^2+682*ep^6*n2^2*x5+
         268*ep^6*n2^2*x5^2+34*ep^6*n2^2*x5^3-144*ep^7-192*ep^7*x5-84*ep^7*x5^2
         -12*ep^7*x5^3+144*ep^7*n2+192*ep^7*n2*x5+84*ep^7*n2*x5^2+12*ep^7*n2*
         x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,0,0,1+n5,0,1,1,0)*(rat(-828*n5+2118*n5^2-
         2252*n5^3+1260*n5^4-388*n5^5+62*n5^6-4*n5^7+2550*n2*n5-5617*n2*n5^2+
         5036*n2*n5^3-2292*n2*n5^4+540*n2*n5^5-59*n2*n5^6+2*n2*n5^7-3314*n2^2*
         n5+6101*n2^2*n5^2-4433*n2^2*n5^3+1545*n2^2*n5^4-249*n2^2*n5^5+14*n2^2*
         n5^6+2355*n2^3*n5-3465*n2^3*n5^2+1911*n2^3*n5^3-455*n2^3*n5^4+38*n2^3*
         n5^5-990*n2^4*n5+1085*n2^4*n5^2-402*n2^4*n5^3+49*n2^4*n5^4+247*n2^5*n5
         -178*n2^5*n5^2+33*n2^5*n5^3-34*n2^6*n5+12*n2^6*n5^2+2*n2^7*n5-828*ep*
         x5+7182*ep*n5+4236*ep*n5*x5-15789*ep*n5^2-6756*ep*n5^2*x5+13942*ep*
         n5^3+5040*ep*n5^3*x5-6180*ep*n5^4-1940*ep*n5^4*x5+1400*ep*n5^5+372*ep*
         n5^5*x5-143*ep*n5^6-28*ep*n5^6*x5+4*ep*n5^7+2550*ep*n2*x5-18940*ep*n2*
         n5-11234*ep*n2*n5*x5+35036*ep*n2*n5^2+15108*ep*n2*n5^2*x5-25182*ep*n2*
         n5^3-9168*ep*n2*n5^3*x5+8576*ep*n2*n5^4+2700*ep*n2*n5^4*x5-1338*ep*n2*
         n5^5-354*ep*n2*n5^5*x5+72*ep*n2*n5^6+14*ep*n2*n5^6*x5-3314*ep*n2^2*x5+
         20420*ep*n2^2*n5+12202*ep*n2^2*n5*x5-30447*ep*n2^2*n5^2-13299*ep*n2^2*
         n5^2*x5+16700*ep*n2^2*n5^3+6180*ep*n2^2*n5^3*x5-3891*ep*n2^2*n5^4-1245
         *ep*n2^2*n5^4*x5+314*ep*n2^2*n5^5+84*ep*n2^2*n5^5*x5+2355*ep*n2^3*x5-
         11520*ep*n2^3*n5-6930*ep*n2^3*n5*x5+12927*ep*n2^3*n5^2+5733*ep*n2^3*
         n5^2*x5-4802*ep*n2^3*n5^3-1820*ep*n2^3*n5^3*x5+575*ep*n2^3*n5^4+190*ep
         *n2^3*n5^4*x5-990*ep*n2^4*x5+3594*ep*n2^4*n5+2170*ep*n2^4*n5*x5-2682*
         ep*n2^4*n5^2-1206*ep*n2^4*n5^2*x5+504*ep*n2^4*n5^3+196*ep*n2^4*n5^3*x5
         +247*ep*n2^5*x5-590*ep*n2^5*n5-356*ep*n2^5*n5*x5+218*ep*n2^5*n5^2+99*
         ep*n2^5*n5^2*x5-34*ep*n2^6*x5+40*ep*n2^6*n5+24*ep*n2^6*n5*x5+2*ep*n2^7
         *x5+7182*ep^2*x5+2118*ep^2*x5^2-26542*ep^2*n5-31578*ep^2*n5*x5-6756*
         ep^2*n5*x5^2+48672*ep^2*n5^2+41826*ep^2*n5^2*x5+7560*ep^2*n5^2*x5^2-
         34229*ep^2*n5^3-24720*ep^2*n5^3*x5-3880*ep^2*n5^3*x5^2+11259*ep^2*n5^4
         +7000*ep^2*n5^4*x5+930*ep^2*n5^4*x5^2-1665*ep^2*n5^5-858*ep^2*n5^5*x5-
         84*ep^2*n5^5*x5^2+81*ep^2*n5^6+28*ep^2*n5^6*x5-18940*ep^2*n2*x5-5617*
         ep^2*n2*x5^2+58338*ep^2*n2*n5+70072*ep^2*n2*n5*x5+15108*ep^2*n2*n5*
         x5^2-86826*ep^2*n2*n5^2-75546*ep^2*n2*n5^2*x5-13752*ep^2*n2*n5^2*x5^2+
         46847*ep^2*n2*n5^3+34304*ep^2*n2*n5^3*x5+5400*ep^2*n2*n5^3*x5^2-10604*
         ep^2*n2*n5^4-6690*ep^2*n2*n5^4*x5-885*ep^2*n2*n5^4*x5^2+821*ep^2*n2*
         n5^5+432*ep^2*n2*n5^5*x5+42*ep^2*n2*n5^5*x5^2+20420*ep^2*n2^2*x5+6101*
         ep^2*n2^2*x5^2-50147*ep^2*n2^2*n5-60894*ep^2*n2^2*n5*x5-13299*ep^2*
         n2^2*n5*x5^2+56640*ep^2*n2^2*n5^2+50100*ep^2*n2^2*n5^2*x5+9270*ep^2*
         n2^2*n5^2*x5^2-20817*ep^2*n2^2*n5^3-15564*ep^2*n2^2*n5^3*x5-2490*ep^2*
         n2^2*n5^3*x5^2+2430*ep^2*n2^2*n5^4+1570*ep^2*n2^2*n5^4*x5+210*ep^2*
         n2^2*n5^4*x5^2-11520*ep^2*n2^3*x5-3465*ep^2*n2^3*x5^2+21078*ep^2*n2^3*
         n5+25854*ep^2*n2^3*n5*x5+5733*ep^2*n2^3*n5*x5^2-15992*ep^2*n2^3*n5^2-
         14406*ep^2*n2^3*n5^2*x5-2730*ep^2*n2^3*n5^2*x5^2+2996*ep^2*n2^3*n5^3+
         2300*ep^2*n2^3*n5^3*x5+380*ep^2*n2^3*n5^3*x5^2+3594*ep^2*n2^4*x5+1085*
         ep^2*n2^4*x5^2-4342*ep^2*n2^4*n5-5364*ep^2*n2^4*n5*x5-1206*ep^2*n2^4*
         n5*x5^2+1650*ep^2*n2^4*n5^2+1512*ep^2*n2^4*n5^2*x5+294*ep^2*n2^4*n5^2*
         x5^2-590*ep^2*n2^5*x5-178*ep^2*n2^5*x5^2+352*ep^2*n2^5*n5+436*ep^2*
         n2^5*n5*x5+99*ep^2*n2^5*n5*x5^2+40*ep^2*n2^6*x5+12*ep^2*n2^6*x5^2-
         26542*ep^3*x5-15789*ep^3*x5^2-2252*ep^3*x5^3+54158*ep^3*n5+97344*ep^3*
         n5*x5+41826*ep^3*n5*x5^2+5040*ep^3*n5*x5^3-79425*ep^3*n5^2-102687*ep^3
         *n5^2*x5-37080*ep^3*n5^2*x5^2-3880*ep^3*n5^2*x5^3+41673*ep^3*n5^3+
         45036*ep^3*n5^3*x5+14000*ep^3*n5^3*x5^2+1240*ep^3*n5^3*x5^3-9035*ep^3*
         n5^4-8325*ep^3*n5^4*x5-2145*ep^3*n5^4*x5^2-140*ep^3*n5^4*x5^3+653*ep^3
         *n5^5+486*ep^3*n5^5*x5+84*ep^3*n5^5*x5^2+58338*ep^3*n2*x5+35036*ep^3*
         n2*x5^2+5036*ep^3*n2*x5^3-95332*ep^3*n2*n5-173652*ep^3*n2*n5*x5-75546*
         ep^3*n2*n5*x5^2-9168*ep^3*n2*n5*x5^3+106855*ep^3*n2*n5^2+140541*ep^3*
         n2*n5^2*x5+51456*ep^3*n2*n5^2*x5^2+5400*ep^3*n2*n5^2*x5^3-38438*ep^3*
         n2*n5^3-42416*ep^3*n2*n5^3*x5-13380*ep^3*n2*n5^3*x5^2-1180*ep^3*n2*
         n5^3*x5^3+4335*ep^3*n2*n5^4+4105*ep^3*n2*n5^4*x5+1080*ep^3*n2*n5^4*
         x5^2+70*ep^3*n2*n5^4*x5^3-50147*ep^3*n2^2*x5-30447*ep^3*n2^2*x5^2-4433
         *ep^3*n2^2*x5^3+61309*ep^3*n2^2*n5+113280*ep^3*n2^2*n5*x5+50100*ep^3*
         n2^2*n5*x5^2+6180*ep^3*n2^2*n5*x5^3-46540*ep^3*n2^2*n5^2-62451*ep^3*
         n2^2*n5^2*x5-23346*ep^3*n2^2*n5^2*x5^2-2490*ep^3*n2^2*n5^2*x5^3+8587*
         ep^3*n2^2*n5^3+9720*ep^3*n2^2*n5^3*x5+3140*ep^3*n2^2*n5^3*x5^2+280*
         ep^3*n2^2*n5^3*x5^3+21078*ep^3*n2^3*x5+12927*ep^3*n2^3*x5^2+1911*ep^3*
         n2^3*x5^3-17082*ep^3*n2^3*n5-31984*ep^3*n2^3*n5*x5-14406*ep^3*n2^3*n5*
         x5^2-1820*ep^3*n2^3*n5*x5^3+6558*ep^3*n2^3*n5^2+8988*ep^3*n2^3*n5^2*x5
         +3450*ep^3*n2^3*n5^2*x5^2+380*ep^3*n2^3*n5^2*x5^3-4342*ep^3*n2^4*x5-
         2682*ep^3*n2^4*x5^2-402*ep^3*n2^4*x5^3+1744*ep^3*n2^4*n5+3300*ep^3*
         n2^4*n5*x5+1512*ep^3*n2^4*n5*x5^2+196*ep^3*n2^4*n5*x5^3+352*ep^3*n2^5*
         x5+218*ep^3*n2^5*x5^2+33*ep^3*n2^5*x5^3+54158*ep^4*x5+48672*ep^4*x5^2+
         13942*ep^4*x5^3+1260*ep^4*x5^4-65886*ep^4*n5-158850*ep^4*n5*x5-102687*
         ep^4*n5*x5^2-24720*ep^4*n5*x5^3-1940*ep^4*n5*x5^4+72376*ep^4*n5^2+
         125019*ep^4*n5^2*x5+67554*ep^4*n5^2*x5^2+14000*ep^4*n5^2*x5^3+930*ep^4
         *n5^2*x5^4-25170*ep^4*n5^3-36140*ep^4*n5^3*x5-16650*ep^4*n5^3*x5^2-
         2860*ep^4*n5^3*x5^3-140*ep^4*n5^3*x5^4+2696*ep^4*n5^4+3265*ep^4*n5^4*
         x5+1215*ep^4*n5^4*x5^2+140*ep^4*n5^4*x5^3-95332*ep^4*n2*x5-86826*ep^4*
         n2*x5^2-25182*ep^4*n2*x5^3-2292*ep^4*n2*x5^4+87140*ep^4*n2*n5+213710*
         ep^4*n2*n5*x5+140541*ep^4*n2*n5*x5^2+34304*ep^4*n2*n5*x5^3+2700*ep^4*
         n2*n5*x5^4-65308*ep^4*n2*n5^2-115314*ep^4*n2*n5^2*x5-63624*ep^4*n2*
         n5^2*x5^2-13380*ep^4*n2*n5^2*x5^3-885*ep^4*n2*n5^2*x5^4+11740*ep^4*n2*
         n5^3+17340*ep^4*n2*n5^3*x5+8210*ep^4*n2*n5^3*x5^2+1440*ep^4*n2*n5^3*
         x5^3+70*ep^4*n2*n5^3*x5^4+61309*ep^4*n2^2*x5+56640*ep^4*n2^2*x5^2+
         16700*ep^4*n2^2*x5^3+1545*ep^4*n2^2*x5^4-37296*ep^4*n2^2*n5-93080*ep^4
         *n2^2*n5*x5-62451*ep^4*n2^2*n5*x5^2-15564*ep^4*n2^2*n5*x5^3-1245*ep^4*
         n2^2*n5*x5^4+14250*ep^4*n2^2*n5^2+25761*ep^4*n2^2*n5^2*x5+14580*ep^4*
         n2^2*n5^2*x5^2+3140*ep^4*n2^2*n5^2*x5^3+210*ep^4*n2^2*n5^2*x5^4-17082*
         ep^4*n2^3*x5-15992*ep^4*n2^3*x5^2-4802*ep^4*n2^3*x5^3-455*ep^4*n2^3*
         x5^4+5170*ep^4*n2^3*n5+13116*ep^4*n2^3*n5*x5+8988*ep^4*n2^3*n5*x5^2+
         2300*ep^4*n2^3*n5*x5^3+190*ep^4*n2^3*n5*x5^4+1744*ep^4*n2^4*x5+1650*
         ep^4*n2^4*x5^2+504*ep^4*n2^4*x5^3+49*ep^4*n2^4*x5^4-65886*ep^5*x5-
         79425*ep^5*x5^2-34229*ep^5*x5^3-6180*ep^5*x5^4-388*ep^5*x5^5+47788*
         ep^5*n5+144752*ep^5*n5*x5+125019*ep^5*n5*x5^2+45036*ep^5*n5*x5^3+7000*
         ep^5*n5*x5^4+372*ep^5*n5*x5^5-34928*ep^5*n5^2-75510*ep^5*n5^2*x5-54210
         *ep^5*n5^2*x5^2-16650*ep^5*n5^2*x5^3-2145*ep^5*n5^2*x5^4-84*ep^5*n5^2*
         x5^5+6036*ep^5*n5^3+10784*ep^5*n5^3*x5+6530*ep^5*n5^3*x5^2+1620*ep^5*
         n5^3*x5^3+140*ep^5*n5^3*x5^4+87140*ep^5*n2*x5+106855*ep^5*n2*x5^2+
         46847*ep^5*n2*x5^3+8576*ep^5*n2*x5^4+540*ep^5*n2*x5^5-42236*ep^5*n2*n5
         -130616*ep^5*n2*n5*x5-115314*ep^5*n2*n5*x5^2-42416*ep^5*n2*n5*x5^3-
         6690*ep^5*n2*n5*x5^4-354*ep^5*n2*n5*x5^5+15860*ep^5*n2*n5^2+35220*ep^5
         *n2*n5^2*x5+26010*ep^5*n2*n5^2*x5^2+8210*ep^5*n2*n5^2*x5^3+1080*ep^5*
         n2*n5^2*x5^4+42*ep^5*n2*n5^2*x5^5-37296*ep^5*n2^2*x5-46540*ep^5*n2^2*
         x5^2-20817*ep^5*n2^2*x5^3-3891*ep^5*n2^2*x5^4-249*ep^5*n2^2*x5^5+9028*
         ep^5*n2^2*n5+28500*ep^5*n2^2*n5*x5+25761*ep^5*n2^2*n5*x5^2+9720*ep^5*
         n2^2*n5*x5^3+1570*ep^5*n2^2*n5*x5^4+84*ep^5*n2^2*n5*x5^5+5170*ep^5*
         n2^3*x5+6558*ep^5*n2^3*x5^2+2996*ep^5*n2^3*x5^3+575*ep^5*n2^3*x5^4+38*
         ep^5*n2^3*x5^5+47788*ep^6*x5+72376*ep^6*x5^2+41673*ep^6*x5^3+11259*
         ep^6*x5^4+1400*ep^6*x5^5+62*ep^6*x5^6-19136*ep^6*n5-69856*ep^6*n5*x5-
         75510*ep^6*n5*x5^2-36140*ep^6*n5*x5^3-8325*ep^6*n5*x5^4-858*ep^6*n5*
         x5^5-28*ep^6*n5*x5^6+6976*ep^6*n5^2+18108*ep^6*n5^2*x5+16176*ep^6*n5^2
         *x5^2+6530*ep^6*n5^2*x5^3+1215*ep^6*n5^2*x5^4+84*ep^6*n5^2*x5^5-42236*
         ep^6*n2*x5-65308*ep^6*n2*x5^2-38438*ep^6*n2*x5^3-10604*ep^6*n2*x5^4-
         1338*ep^6*n2*x5^5-59*ep^6*n2*x5^6+8480*ep^6*n2*n5+31720*ep^6*n2*n5*x5+
         35220*ep^6*n2*n5*x5^2+17340*ep^6*n2*n5*x5^3+4105*ep^6*n2*n5*x5^4+432*
         ep^6*n2*n5*x5^5+14*ep^6*n2*n5*x5^6+9028*ep^6*n2^2*x5+14250*ep^6*n2^2*
         x5^2+8587*ep^6*n2^2*x5^3+2430*ep^6*n2^2*x5^4+314*ep^6*n2^2*x5^5+14*
         ep^6*n2^2*x5^6-19136*ep^7*x5-34928*ep^7*x5^2-25170*ep^7*x5^3-9035*ep^7
         *x5^4-1665*ep^7*x5^5-143*ep^7*x5^6-4*ep^7*x5^7+3264*ep^7*n5+13952*ep^7
         *n5*x5+18108*ep^7*n5*x5^2+10784*ep^7*n5*x5^3+3265*ep^7*n5*x5^4+486*
         ep^7*n5*x5^5+28*ep^7*n5*x5^6+8480*ep^7*n2*x5+15860*ep^7*n2*x5^2+11740*
         ep^7*n2*x5^3+4335*ep^7*n2*x5^4+821*ep^7*n2*x5^5+72*ep^7*n2*x5^6+2*ep^7
         *n2*x5^7+3264*ep^8*x5+6976*ep^8*x5^2+6036*ep^8*x5^3+2696*ep^8*x5^4+653
         *ep^8*x5^5+81*ep^8*x5^6+4*ep^8*x5^7,48-160*n5+204*n5^2-124*n5^3+36*
         n5^4-4*n5^5-172*n2+550*n2*n5-664*n2*n5^2+376*n2*n5^3-100*n2*n5^4+10*n2
         *n5^5+244*n2^2-739*n2^2*n5+827*n2^2*n5^2-421*n2^2*n5^3+97*n2^2*n5^4-8*
         n2^2*n5^5-175*n2^3+494*n2^3*n5-497*n2^3*n5^2+214*n2^3*n5^3-38*n2^3*
         n5^4+2*n2^3*n5^5+67*n2^4-173*n2^4*n5+150*n2^4*n5^2-49*n2^4*n5^3+5*n2^4
         *n5^4-13*n2^5+30*n2^5*n5-21*n2^5*n5^2+4*n2^5*n5^3+n2^6-2*n2^6*n5+n2^6*
         n5^2-552*ep-160*ep*x5+1584*ep*n5+408*ep*n5*x5-1714*ep*n5^2-372*ep*n5^2
         *x5+870*ep*n5^3+144*ep*n5^3*x5-206*ep*n5^4-20*ep*n5^4*x5+18*ep*n5^5+
         1746*ep*n2+550*ep*n2*x5-4681*ep*n2*n5-1328*ep*n2*n5*x5+4633*ep*n2*n5^2
         +1128*ep*n2*n5^2*x5-2093*ep*n2*n5^3-400*ep*n2*n5^3*x5+425*ep*n2*n5^4+
         50*ep*n2*n5^4*x5-30*ep*n2*n5^5-2157*ep*n2^2-739*ep*n2^2*x5+5286*ep*
         n2^2*n5+1654*ep*n2^2*n5*x5-4612*ep*n2^2*n5^2-1263*ep*n2^2*n5^2*x5+1742
         *ep*n2^2*n5^3+388*ep*n2^2*n5^3*x5-271*ep*n2^2*n5^4-40*ep*n2^2*n5^4*x5+
         12*ep*n2^2*n5^5+1325*ep*n2^3+494*ep*n2^3*x5-2880*ep*n2^3*n5-994*ep*
         n2^3*n5*x5+2089*ep*n2^3*n5^2+642*ep*n2^3*n5^2*x5-586*ep*n2^3*n5^3-152*
         ep*n2^3*n5^3*x5+52*ep*n2^3*n5^4+10*ep*n2^3*n5^4*x5-425*ep*n2^4-173*ep*
         n2^4*x5+785*ep*n2^4*n5+300*ep*n2^4*n5*x5-427*ep*n2^4*n5^2-147*ep*n2^4*
         n5^2*x5+67*ep*n2^4*n5^3+20*ep*n2^4*n5^3*x5+67*ep*n2^5+30*ep*n2^5*x5-98
         *ep*n2^5*n5-42*ep*n2^5*n5*x5+31*ep*n2^5*n5^2+12*ep*n2^5*n5^2*x5-4*ep*
         n2^6-2*ep*n2^6*x5+4*ep*n2^6*n5+2*ep*n2^6*n5*x5+2736*ep^2+1584*ep^2*x5+
         204*ep^2*x5^2-6624*ep^2*n5-3428*ep^2*n5*x5-372*ep^2*n5*x5^2+5916*ep^2*
         n5^2+2610*ep^2*n5^2*x5+216*ep^2*n5^2*x5^2-2402*ep^2*n5^3-824*ep^2*n5^3
         *x5-40*ep^2*n5^3*x5^2+432*ep^2*n5^4+90*ep^2*n5^4*x5-26*ep^2*n5^5-7482*
         ep^2*n2-4681*ep^2*n2*x5-664*ep^2*n2*x5^2+16370*ep^2*n2*n5+9266*ep^2*n2
         *n5*x5+1128*ep^2*n2*n5*x5^2-12779*ep^2*n2*n5^2-6279*ep^2*n2*n5^2*x5-
         600*ep^2*n2*n5^2*x5^2+4318*ep^2*n2*n5^3+1700*ep^2*n2*n5^3*x5+100*ep^2*
         n2*n5^3*x5^2-593*ep^2*n2*n5^4-150*ep^2*n2*n5^4*x5+22*ep^2*n2*n5^5+7814
         *ep^2*n2^2+5286*ep^2*n2^2*x5+827*ep^2*n2^2*x5^2-14894*ep^2*n2^2*n5-
         9224*ep^2*n2^2*n5*x5-1263*ep^2*n2^2*n5*x5^2+9537*ep^2*n2^2*n5^2+5226*
         ep^2*n2^2*n5^2*x5+582*ep^2*n2^2*n5^2*x5^2-2380*ep^2*n2^2*n5^3-1084*
         ep^2*n2^2*n5^3*x5-80*ep^2*n2^2*n5^3*x5^2+187*ep^2*n2^2*n5^4+60*ep^2*
         n2^2*n5^4*x5-3935*ep^2*n2^3-2880*ep^2*n2^3*x5-497*ep^2*n2^3*x5^2+6176*
         ep^2*n2^3*n5+4178*ep^2*n2^3*n5*x5+642*ep^2*n2^3*n5*x5^2-2891*ep^2*n2^3
         *n5^2-1758*ep^2*n2^3*n5^2*x5-228*ep^2*n2^3*n5^2*x5^2+398*ep^2*n2^3*
         n5^3+208*ep^2*n2^3*n5^3*x5+20*ep^2*n2^3*n5^3*x5^2+987*ep^2*n2^4+785*
         ep^2*n2^4*x5+150*ep^2*n2^4*x5^2-1154*ep^2*n2^4*n5-854*ep^2*n2^4*n5*x5-
         147*ep^2*n2^4*n5*x5^2+299*ep^2*n2^4*n5^2+201*ep^2*n2^4*n5^2*x5+30*ep^2
         *n2^4*n5^2*x5^2-112*ep^2*n2^5-98*ep^2*n2^5*x5-21*ep^2*n2^5*x5^2+76*
         ep^2*n2^5*n5+62*ep^2*n2^5*n5*x5+12*ep^2*n2^5*n5*x5^2+4*ep^2*n2^6+4*
         ep^2*n2^6*x5+ep^2*n2^6*x5^2-7632*ep^3-6624*ep^3*x5-1714*ep^3*x5^2-124*
         ep^3*x5^3+15168*ep^3*n5+11832*ep^3*n5*x5+2610*ep^3*n5*x5^2+144*ep^3*n5
         *x5^3-10738*ep^3*n5^2-7206*ep^3*n5^2*x5-1236*ep^3*n5^2*x5^2-40*ep^3*
         n5^2*x5^3+3264*ep^3*n5^3+1728*ep^3*n5^3*x5+180*ep^3*n5^3*x5^2-394*ep^3
         *n5^4-130*ep^3*n5^4*x5+12*ep^3*n5^5+17538*ep^3*n2+16370*ep^3*n2*x5+
         4633*ep^3*n2*x5^2+376*ep^3*n2*x5^3-30109*ep^3*n2*n5-25558*ep^3*n2*n5*
         x5-6279*ep^3*n2*n5*x5^2-400*ep^3*n2*n5*x5^3+17412*ep^3*n2*n5^2+12954*
         ep^3*n2*n5^2*x5+2550*ep^3*n2*n5^2*x5^2+100*ep^3*n2*n5^2*x5^3-3913*ep^3
         *n2*n5^3-2372*ep^3*n2*n5^3*x5-300*ep^3*n2*n5^3*x5^2+272*ep^3*n2*n5^4+
         110*ep^3*n2*n5^4*x5-14841*ep^3*n2^2-14894*ep^3*n2^2*x5-4612*ep^3*n2^2*
         x5^2-421*ep^3*n2^2*x5^3+20671*ep^3*n2^2*n5+19074*ep^3*n2^2*n5*x5+5226*
         ep^3*n2^2*n5*x5^2+388*ep^3*n2^2*n5*x5^3-8663*ep^3*n2^2*n5^2-7140*ep^3*
         n2^2*n5^2*x5-1626*ep^3*n2^2*n5^2*x5^2-80*ep^3*n2^2*n5^2*x5^3+1073*ep^3
         *n2^2*n5^3+748*ep^3*n2^2*n5^3*x5+120*ep^3*n2^2*n5^3*x5^2+5725*ep^3*
         n2^3+6176*ep^3*n2^3*x5+2089*ep^3*n2^3*x5^2+214*ep^3*n2^3*x5^3-5782*
         ep^3*n2^3*n5-5782*ep^3*n2^3*n5*x5-1758*ep^3*n2^3*n5*x5^2-152*ep^3*n2^3
         *n5*x5^3+1317*ep^3*n2^3*n5^2+1194*ep^3*n2^3*n5^2*x5+312*ep^3*n2^3*n5^2
         *x5^2+20*ep^3*n2^3*n5^2*x5^3-992*ep^3*n2^4-1154*ep^3*n2^4*x5-427*ep^3*
         n2^4*x5^2-49*ep^3*n2^4*x5^3+552*ep^3*n2^4*n5+598*ep^3*n2^4*n5*x5+201*
         ep^3*n2^4*n5*x5^2+20*ep^3*n2^4*n5*x5^3+60*ep^3*n2^5+76*ep^3*n2^5*x5+31
         *ep^3*n2^5*x5^2+4*ep^3*n2^5*x5^3+13104*ep^4+15168*ep^4*x5+5916*ep^4*
         x5^2+870*ep^4*x5^3+36*ep^4*x5^4-20544*ep^4*n5-21476*ep^4*n5*x5-7206*
         ep^4*n5*x5^2-824*ep^4*n5*x5^3-20*ep^4*n5*x5^4+10812*ep^4*n5^2+9792*
         ep^4*n5^2*x5+2592*ep^4*n5^2*x5^2+180*ep^4*n5^2*x5^3-2184*ep^4*n5^3-
         1576*ep^4*n5^3*x5-260*ep^4*n5^3*x5^2+132*ep^4*n5^4+60*ep^4*n5^4*x5-
         24282*ep^4*n2-30109*ep^4*n2*x5-12779*ep^4*n2*x5^2-2093*ep^4*n2*x5^3-
         100*ep^4*n2*x5^4+30722*ep^4*n2*n5+34824*ep^4*n2*n5*x5+12954*ep^4*n2*n5
         *x5^2+1700*ep^4*n2*n5*x5^3+50*ep^4*n2*n5*x5^4-11718*ep^4*n2*n5^2-11739
         *ep^4*n2*n5^2*x5-3558*ep^4*n2*n5^2*x5^2-300*ep^4*n2*n5^2*x5^3+1314*
         ep^4*n2*n5^3+1088*ep^4*n2*n5^3*x5+220*ep^4*n2*n5^3*x5^2+15584*ep^4*
         n2^2+20671*ep^4*n2^2*x5+9537*ep^4*n2^2*x5^2+1742*ep^4*n2^2*x5^3+97*
         ep^4*n2^2*x5^4-14136*ep^4*n2^2*n5-17326*ep^4*n2^2*n5*x5-7140*ep^4*n2^2
         *n5*x5^2-1084*ep^4*n2^2*n5*x5^3-40*ep^4*n2^2*n5*x5^4+2916*ep^4*n2^2*
         n5^2+3219*ep^4*n2^2*n5^2*x5+1122*ep^4*n2^2*n5^2*x5^2+120*ep^4*n2^2*
         n5^2*x5^3-4080*ep^4*n2^3-5782*ep^4*n2^3*x5-2891*ep^4*n2^3*x5^2-586*
         ep^4*n2^3*x5^3-38*ep^4*n2^3*x5^4+1996*ep^4*n2^3*n5+2634*ep^4*n2^3*n5*
         x5+1194*ep^4*n2^3*n5*x5^2+208*ep^4*n2^3*n5*x5^3+10*ep^4*n2^3*n5*x5^4+
         364*ep^4*n2^4+552*ep^4*n2^4*x5+299*ep^4*n2^4*x5^2+67*ep^4*n2^4*x5^3+5*
         ep^4*n2^4*x5^4-14184*ep^5-20544*ep^5*x5-10738*ep^5*x5^2-2402*ep^5*x5^3
         -206*ep^5*x5^4-4*ep^5*x5^5+16464*ep^5*n5+21624*ep^5*n5*x5+9792*ep^5*n5
         *x5^2+1728*ep^5*n5*x5^3+90*ep^5*n5*x5^4-5728*ep^5*n5^2-6552*ep^5*n5^2*
         x5-2364*ep^5*n5^2*x5^2-260*ep^5*n5^2*x5^3+576*ep^5*n5^3+528*ep^5*n5^3*
         x5+120*ep^5*n5^3*x5^2+19860*ep^5*n2+30722*ep^5*n2*x5+17412*ep^5*n2*
         x5^2+4318*ep^5*n2*x5^3+425*ep^5*n2*x5^4+10*ep^5*n2*x5^5-16492*ep^5*n2*
         n5-23436*ep^5*n2*n5*x5-11739*ep^5*n2*n5*x5^2-2372*ep^5*n2*n5*x5^3-150*
         ep^5*n2*n5*x5^4+3116*ep^5*n2*n5^2+3942*ep^5*n2*n5^2*x5+1632*ep^5*n2*
         n5^2*x5^2+220*ep^5*n2*n5^2*x5^3-8580*ep^5*n2^2-14136*ep^5*n2^2*x5-8663
         *ep^5*n2^2*x5^2-2380*ep^5*n2^2*x5^3-271*ep^5*n2^2*x5^4-8*ep^5*n2^2*
         x5^5+3812*ep^5*n2^2*n5+5832*ep^5*n2^2*n5*x5+3219*ep^5*n2^2*n5*x5^2+748
         *ep^5*n2^2*n5*x5^3+60*ep^5*n2^2*n5*x5^4+1140*ep^5*n2^3+1996*ep^5*n2^3*
         x5+1317*ep^5*n2^3*x5^2+398*ep^5*n2^3*x5^3+52*ep^5*n2^3*x5^4+2*ep^5*
         n2^3*x5^5+9456*ep^6+16464*ep^6*x5+10812*ep^6*x5^2+3264*ep^6*x5^3+432*
         ep^6*x5^4+18*ep^6*x5^5-7232*ep^6*n5-11456*ep^6*n5*x5-6552*ep^6*n5*x5^2
         -1576*ep^6*n5*x5^3-130*ep^6*n5*x5^4+1248*ep^6*n5^2+1728*ep^6*n5^2*x5+
         792*ep^6*n5^2*x5^2+120*ep^6*n5^2*x5^3-8888*ep^6*n2-16492*ep^6*n2*x5-
         11718*ep^6*n2*x5^2-3913*ep^6*n2*x5^3-593*ep^6*n2*x5^4-30*ep^6*n2*x5^5+
         3640*ep^6*n2*n5+6232*ep^6*n2*n5*x5+3942*ep^6*n2*n5*x5^2+1088*ep^6*n2*
         n5*x5^3+110*ep^6*n2*n5*x5^4+1936*ep^6*n2^2+3812*ep^6*n2^2*x5+2916*ep^6
         *n2^2*x5^2+1073*ep^6*n2^2*x5^3+187*ep^6*n2^2*x5^4+12*ep^6*n2^2*x5^5-
         3552*ep^7-7232*ep^7*x5-5728*ep^7*x5^2-2184*ep^7*x5^3-394*ep^7*x5^4-26*
         ep^7*x5^5+1344*ep^7*n5+2496*ep^7*n5*x5+1728*ep^7*n5*x5^2+528*ep^7*n5*
         x5^3+60*ep^7*n5*x5^4+1680*ep^7*n2+3640*ep^7*n2*x5+3116*ep^7*n2*x5^2+
         1314*ep^7*n2*x5^3+272*ep^7*n2*x5^4+22*ep^7*n2*x5^5+576*ep^8+1344*ep^8*
         x5+1248*ep^8*x5^2+576*ep^8*x5^3+132*ep^8*x5^4+12*ep^8*x5^5))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,0,1,1+n5,0,0,1,0)*(rat(60*n5-20*n5^2-93*n5^3
         +112*n5^4-54*n5^5+12*n5^6-n5^7-377*n2*n5+398*n2*n5^2-44*n2*n5^3-123*n2
         *n5^4+72*n2*n5^5-15*n2*n5^6+n2*n5^7+815*n2^2*n5-955*n2^2*n5^2+396*n2^2
         *n5^3-40*n2^2*n5^4-15*n2^2*n5^5+3*n2^2*n5^6-859*n2^3*n5+908*n2^3*n5^2-
         368*n2^3*n5^3+64*n2^3*n5^4-3*n2^3*n5^5+497*n2^4*n5-415*n2^4*n5^2+123*
         n2^4*n5^3-13*n2^4*n5^4-162*n2^5*n5+92*n2^5*n5^2-14*n2^5*n5^3+28*n2^6*
         n5-8*n2^6*n5^2-2*n2^7*n5+60*ep*x5-330*ep*n5-40*ep*n5*x5-87*ep*n5^2-279
         *ep*n5^2*x5+668*ep*n5^3+448*ep*n5^3*x5-547*ep*n5^4-270*ep*n5^4*x5+185*
         ep*n5^5+72*ep*n5^5*x5-26*ep*n5^6-7*ep*n5^6*x5+ep*n5^7-377*ep*n2*x5+
         2049*ep*n2*n5+796*ep*n2*n5*x5-1535*ep*n2*n5^2-132*ep*n2*n5^2*x5-273*ep
         *n2*n5^3-492*ep*n2*n5^3*x5+590*ep*n2*n5^4+360*ep*n2*n5^4*x5-194*ep*n2*
         n5^5-90*ep*n2*n5^5*x5+19*ep*n2*n5^6+7*ep*n2*n5^6*x5+815*ep*n2^2*x5-
         3935*ep*n2^2*n5-1910*ep*n2^2*n5*x5+3597*ep*n2^2*n5^2+1188*ep*n2^2*n5^2
         *x5-974*ep*n2^2*n5^3-160*ep*n2^2*n5^3*x5-15*ep*n2^2*n5^4-75*ep*n2^2*
         n5^4*x5+29*ep*n2^2*n5^5+18*ep*n2^2*n5^5*x5-859*ep*n2^3*x5+3421*ep*n2^3
         *n5+1816*ep*n2^3*n5*x5-2752*ep*n2^3*n5^2-1104*ep*n2^3*n5^2*x5+735*ep*
         n2^3*n5^3+256*ep*n2^3*n5^3*x5-58*ep*n2^3*n5^4-15*ep*n2^3*n5^4*x5+497*
         ep*n2^4*x5-1501*ep*n2^4*n5-830*ep*n2^4*n5*x5+862*ep*n2^4*n5^2+369*ep*
         n2^4*n5^2*x5-131*ep*n2^4*n5^3-52*ep*n2^4*n5^3*x5-162*ep*n2^5*x5+326*ep
         *n2^5*n5+184*ep*n2^5*n5*x5-96*ep*n2^5*n5^2-42*ep*n2^5*n5^2*x5+28*ep*
         n2^6*x5-28*ep*n2^6*n5-16*ep*n2^6*n5*x5-2*ep*n2^7*x5-330*ep^2*x5-20*
         ep^2*x5^2+666*ep^2*n5-174*ep^2*n5*x5-279*ep^2*n5*x5^2+806*ep^2*n5^2+
         2004*ep^2*n5^2*x5+672*ep^2*n5^2*x5^2-1770*ep^2*n5^3-2188*ep^2*n5^3*x5-
         540*ep^2*n5^3*x5^2+972*ep^2*n5^4+925*ep^2*n5^4*x5+180*ep^2*n5^4*x5^2-
         208*ep^2*n5^5-156*ep^2*n5^5*x5-21*ep^2*n5^5*x5^2+14*ep^2*n5^6+7*ep^2*
         n5^6*x5+2049*ep^2*n2*x5+398*ep^2*n2*x5^2-4389*ep^2*n2*n5-3070*ep^2*n2*
         n5*x5-132*ep^2*n2*n5*x5^2+1831*ep^2*n2*n5^2-819*ep^2*n2*n5^2*x5-738*
         ep^2*n2*n5^2*x5^2+1285*ep^2*n2*n5^3+2360*ep^2*n2*n5^3*x5+720*ep^2*n2*
         n5^3*x5^2-867*ep^2*n2*n5^4-970*ep^2*n2*n5^4*x5-225*ep^2*n2*n5^4*x5^2+
         128*ep^2*n2*n5^5+114*ep^2*n2*n5^5*x5+21*ep^2*n2*n5^5*x5^2-3935*ep^2*
         n2^2*x5-955*ep^2*n2^2*x5^2+7411*ep^2*n2^2*n5+7194*ep^2*n2^2*n5*x5+1188
         *ep^2*n2^2*n5*x5^2-4803*ep^2*n2^2*n5^2-2922*ep^2*n2^2*n5^2*x5-240*ep^2
         *n2^2*n5^2*x5^2+635*ep^2*n2^2*n5^3-60*ep^2*n2^2*n5^3*x5-150*ep^2*n2^2*
         n5^3*x5^2+75*ep^2*n2^2*n5^4+145*ep^2*n2^2*n5^4*x5+45*ep^2*n2^2*n5^4*
         x5^2+3421*ep^2*n2^3*x5+908*ep^2*n2^3*x5^2-5054*ep^2*n2^3*n5-5504*ep^2*
         n2^3*n5*x5-1104*ep^2*n2^3*n5*x5^2+2720*ep^2*n2^3*n5^2+2205*ep^2*n2^3*
         n5^2*x5+384*ep^2*n2^3*n5^2*x5^2-350*ep^2*n2^3*n5^3-232*ep^2*n2^3*n5^3*
         x5-30*ep^2*n2^3*n5^3*x5^2-1501*ep^2*n2^4*x5-415*ep^2*n2^4*x5^2+1506*
         ep^2*n2^4*n5+1724*ep^2*n2^4*n5*x5+369*ep^2*n2^4*n5*x5^2-444*ep^2*n2^4*
         n5^2-393*ep^2*n2^4*n5^2*x5-78*ep^2*n2^4*n5^2*x5^2+326*ep^2*n2^5*x5+92*
         ep^2*n2^5*x5^2-164*ep^2*n2^5*n5-192*ep^2*n2^5*n5*x5-42*ep^2*n2^5*n5*
         x5^2-28*ep^2*n2^6*x5-8*ep^2*n2^6*x5^2+666*ep^3*x5-87*ep^3*x5^2-93*ep^3
         *x5^3-490*ep^3*n5+1612*ep^3*n5*x5+2004*ep^3*n5*x5^2+448*ep^3*n5*x5^3-
         1927*ep^3*n5^2-5310*ep^3*n5^2*x5-3282*ep^3*n5^2*x5^2-540*ep^3*n5^2*
         x5^3+2217*ep^3*n5^3+3888*ep^3*n5^3*x5+1850*ep^3*n5^3*x5^2+240*ep^3*
         n5^3*x5^3-749*ep^3*n5^4-1040*ep^3*n5^4*x5-390*ep^3*n5^4*x5^2-35*ep^3*
         n5^4*x5^3+77*ep^3*n5^5+84*ep^3*n5^5*x5+21*ep^3*n5^5*x5^2-4389*ep^3*n2*
         x5-1535*ep^3*n2*x5^2-44*ep^3*n2*x5^3+4559*ep^3*n2*n5+3662*ep^3*n2*n5*
         x5-819*ep^3*n2*n5*x5^2-492*ep^3*n2*n5*x5^3-214*ep^3*n2*n5^2+3855*ep^3*
         n2*n5^2*x5+3540*ep^3*n2*n5^2*x5^2+720*ep^3*n2*n5^2*x5^3-1597*ep^3*n2*
         n5^3-3468*ep^3*n2*n5^3*x5-1940*ep^3*n2*n5^3*x5^2-300*ep^3*n2*n5^3*x5^3
         +404*ep^3*n2*n5^4+640*ep^3*n2*n5^4*x5+285*ep^3*n2*n5^4*x5^2+35*ep^3*n2
         *n5^4*x5^3+7411*ep^3*n2^2*x5+3597*ep^3*n2^2*x5^2+396*ep^3*n2^2*x5^3-
         6717*ep^3*n2^2*n5-9606*ep^3*n2^2*n5*x5-2922*ep^3*n2^2*n5*x5^2-160*ep^3
         *n2^2*n5*x5^3+2590*ep^3*n2^2*n5^2+1905*ep^3*n2^2*n5^2*x5-90*ep^3*n2^2*
         n5^2*x5^2-150*ep^3*n2^2*n5^2*x5^3-47*ep^3*n2^2*n5^3+300*ep^3*n2^2*n5^3
         *x5+290*ep^3*n2^2*n5^3*x5^2+60*ep^3*n2^2*n5^3*x5^3-5054*ep^3*n2^3*x5-
         2752*ep^3*n2^3*x5^2-368*ep^3*n2^3*x5^3+3262*ep^3*n2^3*n5+5440*ep^3*
         n2^3*n5*x5+2205*ep^3*n2^3*n5*x5^2+256*ep^3*n2^3*n5*x5^3-868*ep^3*n2^3*
         n5^2-1050*ep^3*n2^3*n5^2*x5-348*ep^3*n2^3*n5^2*x5^2-30*ep^3*n2^3*n5^2*
         x5^3+1506*ep^3*n2^4*x5+862*ep^3*n2^4*x5^2+123*ep^3*n2^4*x5^3-500*ep^3*
         n2^4*n5-888*ep^3*n2^4*n5*x5-393*ep^3*n2^4*n5*x5^2-52*ep^3*n2^4*n5*x5^3
         -164*ep^3*n2^5*x5-96*ep^3*n2^5*x5^2-14*ep^3*n2^5*x5^3-490*ep^4*x5+806*
         ep^4*x5^2+668*ep^4*x5^3+112*ep^4*x5^4-190*ep^4*n5-3854*ep^4*n5*x5-5310
         *ep^4*n5*x5^2-2188*ep^4*n5*x5^3-270*ep^4*n5*x5^4+2108*ep^4*n5^2+6651*
         ep^4*n5^2*x5+5832*ep^4*n5^2*x5^2+1850*ep^4*n5^2*x5^3+180*ep^4*n5^2*
         x5^4-1330*ep^4*n5^3-2996*ep^4*n5^3*x5-2080*ep^4*n5^3*x5^2-520*ep^4*
         n5^3*x5^3-35*ep^4*n5^3*x5^4+212*ep^4*n5^4+385*ep^4*n5^4*x5+210*ep^4*
         n5^4*x5^2+35*ep^4*n5^4*x5^3+4559*ep^4*n2*x5+1831*ep^4*n2*x5^2-273*ep^4
         *n2*x5^3-123*ep^4*n2*x5^4-2182*ep^4*n2*n5-428*ep^4*n2*n5*x5+3855*ep^4*
         n2*n5*x5^2+2360*ep^4*n2*n5*x5^3+360*ep^4*n2*n5*x5^4-916*ep^4*n2*n5^2-
         4791*ep^4*n2*n5^2*x5-5202*ep^4*n2*n5^2*x5^2-1940*ep^4*n2*n5^2*x5^3-225
         *ep^4*n2*n5^2*x5^4+630*ep^4*n2*n5^3+1616*ep^4*n2*n5^3*x5+1280*ep^4*n2*
         n5^3*x5^2+380*ep^4*n2*n5^3*x5^3+35*ep^4*n2*n5^3*x5^4-6717*ep^4*n2^2*x5
         -4803*ep^4*n2^2*x5^2-974*ep^4*n2^2*x5^3-40*ep^4*n2^2*x5^4+2866*ep^4*
         n2^2*n5+5180*ep^4*n2^2*n5*x5+1905*ep^4*n2^2*n5*x5^2-60*ep^4*n2^2*n5*
         x5^3-75*ep^4*n2^2*n5*x5^4-428*ep^4*n2^2*n5^2-141*ep^4*n2^2*n5^2*x5+450
         *ep^4*n2^2*n5^2*x5^2+290*ep^4*n2^2*n5^2*x5^3+45*ep^4*n2^2*n5^2*x5^4+
         3262*ep^4*n2^3*x5+2720*ep^4*n2^3*x5^2+735*ep^4*n2^3*x5^3+64*ep^4*n2^3*
         x5^4-770*ep^4*n2^3*n5-1736*ep^4*n2^3*n5*x5-1050*ep^4*n2^3*n5*x5^2-232*
         ep^4*n2^3*n5*x5^3-15*ep^4*n2^3*n5*x5^4-500*ep^4*n2^4*x5-444*ep^4*n2^4*
         x5^2-131*ep^4*n2^4*x5^3-13*ep^4*n2^4*x5^4-190*ep^5*x5-1927*ep^5*x5^2-
         1770*ep^5*x5^3-547*ep^5*x5^4-54*ep^5*x5^5+540*ep^5*n5+4216*ep^5*n5*x5+
         6651*ep^5*n5*x5^2+3888*ep^5*n5*x5^3+925*ep^5*n5*x5^4+72*ep^5*n5*x5^5-
         1104*ep^5*n5^2-3990*ep^5*n5^2*x5-4494*ep^5*n5^2*x5^2-2080*ep^5*n5^2*
         x5^3-390*ep^5*n5^2*x5^4-21*ep^5*n5^2*x5^5+308*ep^5*n5^3+848*ep^5*n5^3*
         x5+770*ep^5*n5^3*x5^2+280*ep^5*n5^3*x5^3+35*ep^5*n5^3*x5^4-2182*ep^5*
         n2*x5-214*ep^5*n2*x5^2+1285*ep^5*n2*x5^3+590*ep^5*n2*x5^4+72*ep^5*n2*
         x5^5+260*ep^5*n2*n5-1832*ep^5*n2*n5*x5-4791*ep^5*n2*n5*x5^2-3468*ep^5*
         n2*n5*x5^3-970*ep^5*n2*n5*x5^4-90*ep^5*n2*n5*x5^5+436*ep^5*n2*n5^2+
         1890*ep^5*n2*n5^2*x5+2424*ep^5*n2*n5^2*x5^2+1280*ep^5*n2*n5^2*x5^3+285
         *ep^5*n2*n5^2*x5^4+21*ep^5*n2*n5^2*x5^5+2866*ep^5*n2^2*x5+2590*ep^5*
         n2^2*x5^2+635*ep^5*n2^2*x5^3-15*ep^5*n2^2*x5^4-15*ep^5*n2^2*x5^5-440*
         ep^5*n2^2*n5-856*ep^5*n2^2*n5*x5-141*ep^5*n2^2*n5*x5^2+300*ep^5*n2^2*
         n5*x5^3+145*ep^5*n2^2*n5*x5^4+18*ep^5*n2^2*n5*x5^5-770*ep^5*n2^3*x5-
         868*ep^5*n2^3*x5^2-350*ep^5*n2^3*x5^3-58*ep^5*n2^3*x5^4-3*ep^5*n2^3*
         x5^5+540*ep^6*x5+2108*ep^6*x5^2+2217*ep^6*x5^3+972*ep^6*x5^4+185*ep^6*
         x5^5+12*ep^6*x5^6-320*ep^6*n5-2208*ep^6*n5*x5-3990*ep^6*n5*x5^2-2996*
         ep^6*n5*x5^3-1040*ep^6*n5*x5^4-156*ep^6*n5*x5^5-7*ep^6*n5*x5^6+224*
         ep^6*n5^2+924*ep^6*n5^2*x5+1272*ep^6*n5^2*x5^2+770*ep^6*n5^2*x5^3+210*
         ep^6*n5^2*x5^4+21*ep^6*n5^2*x5^5+260*ep^6*n2*x5-916*ep^6*n2*x5^2-1597*
         ep^6*n2*x5^3-867*ep^6*n2*x5^4-194*ep^6*n2*x5^5-15*ep^6*n2*x5^6+80*ep^6
         *n2*n5+872*ep^6*n2*n5*x5+1890*ep^6*n2*n5*x5^2+1616*ep^6*n2*n5*x5^3+640
         *ep^6*n2*n5*x5^4+114*ep^6*n2*n5*x5^5+7*ep^6*n2*n5*x5^6-440*ep^6*n2^2*
         x5-428*ep^6*n2^2*x5^2-47*ep^6*n2^2*x5^3+75*ep^6*n2^2*x5^4+29*ep^6*n2^2
         *x5^5+3*ep^6*n2^2*x5^6-320*ep^7*x5-1104*ep^7*x5^2-1330*ep^7*x5^3-749*
         ep^7*x5^4-208*ep^7*x5^5-26*ep^7*x5^6-ep^7*x5^7+64*ep^7*n5+448*ep^7*n5*
         x5+924*ep^7*n5*x5^2+848*ep^7*n5*x5^3+385*ep^7*n5*x5^4+84*ep^7*n5*x5^5+
         7*ep^7*n5*x5^6+80*ep^7*n2*x5+436*ep^7*n2*x5^2+630*ep^7*n2*x5^3+404*
         ep^7*n2*x5^4+128*ep^7*n2*x5^5+19*ep^7*n2*x5^6+ep^7*n2*x5^7+64*ep^8*x5+
         224*ep^8*x5^2+308*ep^8*x5^3+212*ep^8*x5^4+77*ep^8*x5^5+14*ep^8*x5^6+
         ep^8*x5^7,12-34*n5+34*n5^2-14*n5^3+2*n5^4-55*n2+153*n2*n5-149*n2*n5^2+
         59*n2*n5^3-8*n2*n5^4+101*n2^2-274*n2^2*n5+257*n2^2*n5^2-96*n2^2*n5^3+
         12*n2^2*n5^4-94*n2^3+246*n2^3*n5-218*n2^3*n5^2+74*n2^3*n5^3-8*n2^3*
         n5^4+46*n2^4-114*n2^4*n5+92*n2^4*n5^2-26*n2^4*n5^3+2*n2^4*n5^4-11*n2^5
         +25*n2^5*n5-17*n2^5*n5^2+3*n2^5*n5^3+n2^6-2*n2^6*n5+n2^6*n5^2-120*ep-
         34*ep*x5+282*ep*n5+68*ep*n5*x5-228*ep*n5^2-42*ep*n5^2*x5+74*ep*n5^3+8*
         ep*n5^3*x5-8*ep*n5^4+495*ep*n2+153*ep*n2*x5-1114*ep*n2*n5-298*ep*n2*n5
         *x5+847*ep*n2*n5^2+177*ep*n2*n5^2*x5-252*ep*n2*n5^3-32*ep*n2*n5^3*x5+
         24*ep*n2*n5^4-808*ep*n2^2-274*ep*n2^2*x5+1716*ep*n2^2*n5+514*ep*n2^2*
         n5*x5-1196*ep*n2^2*n5^2-288*ep*n2^2*n5^2*x5+312*ep*n2^2*n5^3+48*ep*
         n2^2*n5^3*x5-24*ep*n2^2*n5^4+658*ep*n2^3+246*ep*n2^3*x5-1288*ep*n2^3*
         n5-436*ep*n2^3*n5*x5+786*ep*n2^3*n5^2+222*ep*n2^3*n5^2*x5-164*ep*n2^3*
         n5^3-32*ep*n2^3*n5^3*x5+8*ep*n2^3*n5^4-276*ep*n2^4-114*ep*n2^4*x5+478*
         ep*n2^4*n5+184*ep*n2^4*n5*x5-232*ep*n2^4*n5^2-78*ep*n2^4*n5^2*x5+30*ep
         *n2^4*n5^3+8*ep*n2^4*n5^3*x5+55*ep*n2^5+25*ep*n2^5*x5-78*ep*n2^5*n5-34
         *ep*n2^5*n5*x5+23*ep*n2^5*n5^2+9*ep*n2^5*n5^2*x5-4*ep*n2^6-2*ep*n2^6*
         x5+4*ep*n2^6*n5+2*ep*n2^6*n5*x5+504*ep^2+282*ep^2*x5+34*ep^2*x5^2-954*
         ep^2*n5-456*ep^2*n5*x5-42*ep^2*n5*x5^2+598*ep^2*n5^2+222*ep^2*n5^2*x5+
         12*ep^2*n5^2*x5^2-142*ep^2*n5^3-32*ep^2*n5^3*x5+10*ep^2*n5^4-1839*ep^2
         *n2-1114*ep^2*n2*x5-149*ep^2*n2*x5^2+3233*ep^2*n2*n5+1694*ep^2*n2*n5*
         x5+177*ep^2*n2*n5*x5^2-1825*ep^2*n2*n5^2-756*ep^2*n2*n5^2*x5-48*ep^2*
         n2*n5^2*x5^2+371*ep^2*n2*n5^3+96*ep^2*n2*n5^3*x5-20*ep^2*n2*n5^4+2599*
         ep^2*n2^2+1716*ep^2*n2^2*x5+257*ep^2*n2^2*x5^2-4126*ep^2*n2^2*n5-2392*
         ep^2*n2^2*n5*x5-288*ep^2*n2^2*n5*x5^2+1997*ep^2*n2^2*n5^2+936*ep^2*
         n2^2*n5^2*x5+72*ep^2*n2^2*n5^2*x5^2-316*ep^2*n2^2*n5^3-96*ep^2*n2^2*
         n5^3*x5+10*ep^2*n2^2*n5^4-1777*ep^2*n2^3-1288*ep^2*n2^3*x5-218*ep^2*
         n2^3*x5^2+2425*ep^2*n2^3*n5+1572*ep^2*n2^3*n5*x5+222*ep^2*n2^3*n5*x5^2
         -911*ep^2*n2^3*n5^2-492*ep^2*n2^3*n5^2*x5-48*ep^2*n2^3*n5^2*x5^2+87*
         ep^2*n2^3*n5^3+32*ep^2*n2^3*n5^3*x5+597*ep^2*n2^4+478*ep^2*n2^4*x5+92*
         ep^2*n2^4*x5^2-634*ep^2*n2^4*n5-464*ep^2*n2^4*n5*x5-78*ep^2*n2^4*n5*
         x5^2+141*ep^2*n2^4*n5^2+90*ep^2*n2^4*n5^2*x5+12*ep^2*n2^4*n5^2*x5^2-88
         *ep^2*n2^5-78*ep^2*n2^5*x5-17*ep^2*n2^5*x5^2+56*ep^2*n2^5*n5+46*ep^2*
         n2^5*n5*x5+9*ep^2*n2^5*n5*x5^2+4*ep^2*n2^6+4*ep^2*n2^6*x5+ep^2*n2^6*
         x5^2-1152*ep^3-954*ep^3*x5-228*ep^3*x5^2-14*ep^3*x5^3+1686*ep^3*n5+
         1196*ep^3*n5*x5+222*ep^3*n5*x5^2+8*ep^3*n5*x5^3-768*ep^3*n5^2-426*ep^3
         *n5^2*x5-48*ep^3*n5^2*x5^2+118*ep^3*n5^3+40*ep^3*n5^3*x5-4*ep^3*n5^4+
         3633*ep^3*n2+3233*ep^3*n2*x5+847*ep^3*n2*x5^2+59*ep^3*n2*x5^3-4744*
         ep^3*n2*n5-3650*ep^3*n2*n5*x5-756*ep^3*n2*n5*x5^2-32*ep^3*n2*n5*x5^3+
         1829*ep^3*n2*n5^2+1113*ep^3*n2*n5^2*x5+144*ep^3*n2*n5^2*x5^2-214*ep^3*
         n2*n5^3-80*ep^3*n2*n5^3*x5+4*ep^3*n2*n5^4-4282*ep^3*n2^2-4126*ep^3*
         n2^2*x5-1196*ep^3*n2^2*x5^2-96*ep^3*n2^2*x5^3+4736*ep^3*n2^2*n5+3994*
         ep^3*n2^2*n5*x5+936*ep^3*n2^2*n5*x5^2+48*ep^3*n2^2*n5*x5^3-1398*ep^3*
         n2^2*n5^2-948*ep^3*n2^2*n5^2*x5-144*ep^3*n2^2*n5^2*x5^2+96*ep^3*n2^2*
         n5^3+40*ep^3*n2^2*n5^3*x5+2305*ep^3*n2^3+2425*ep^3*n2^3*x5+786*ep^3*
         n2^3*x5^2+74*ep^3*n2^3*x5^3-1946*ep^3*n2^3*n5-1822*ep^3*n2^3*n5*x5-492
         *ep^3*n2^3*n5*x5^2-32*ep^3*n2^3*n5*x5^3+337*ep^3*n2^3*n5^2+261*ep^3*
         n2^3*n5^2*x5+48*ep^3*n2^3*n5^2*x5^2-548*ep^3*n2^4-634*ep^3*n2^4*x5-232
         *ep^3*n2^4*x5^2-26*ep^3*n2^4*x5^3+268*ep^3*n2^4*n5+282*ep^3*n2^4*n5*x5
         +90*ep^3*n2^4*n5*x5^2+8*ep^3*n2^4*n5*x5^3+44*ep^3*n2^5+56*ep^3*n2^5*x5
         +23*ep^3*n2^5*x5^2+3*ep^3*n2^5*x5^3+1548*ep^4+1686*ep^4*x5+598*ep^4*
         x5^2+74*ep^4*x5^3+2*ep^4*x5^4-1644*ep^4*n5-1536*ep^4*n5*x5-426*ep^4*n5
         *x5^2-32*ep^4*n5*x5^3+484*ep^4*n5^2+354*ep^4*n5^2*x5+60*ep^4*n5^2*x5^2
         -36*ep^4*n5^3-16*ep^4*n5^3*x5-4086*ep^4*n2-4744*ep^4*n2*x5-1825*ep^4*
         n2*x5^2-252*ep^4*n2*x5^3-8*ep^4*n2*x5^4+3656*ep^4*n2*n5+3658*ep^4*n2*
         n5*x5+1113*ep^4*n2*n5*x5^2+96*ep^4*n2*n5*x5^3-822*ep^4*n2*n5^2-642*
         ep^4*n2*n5^2*x5-120*ep^4*n2*n5^2*x5^2+36*ep^4*n2*n5^3+16*ep^4*n2*n5^3*
         x5+3790*ep^4*n2^2+4736*ep^4*n2^2*x5+1997*ep^4*n2^2*x5^2+312*ep^4*n2^2*
         x5^3+12*ep^4*n2^2*x5^4-2572*ep^4*n2^2*n5-2796*ep^4*n2^2*n5*x5-948*ep^4
         *n2^2*n5*x5^2-96*ep^4*n2^2*n5*x5^3+338*ep^4*n2^2*n5^2+288*ep^4*n2^2*
         n5^2*x5+60*ep^4*n2^2*n5^2*x5^2-1432*ep^4*n2^3-1946*ep^4*n2^3*x5-911*
         ep^4*n2^3*x5^2-164*ep^4*n2^3*x5^3-8*ep^4*n2^3*x5^4+560*ep^4*n2^3*n5+
         674*ep^4*n2^3*n5*x5+261*ep^4*n2^3*n5*x5^2+32*ep^4*n2^3*n5*x5^3+180*
         ep^4*n2^4+268*ep^4*n2^4*x5+141*ep^4*n2^4*x5^2+30*ep^4*n2^4*x5^3+2*ep^4
         *n2^4*x5^4-1224*ep^5-1644*ep^5*x5-768*ep^5*x5^2-142*ep^5*x5^3-8*ep^5*
         x5^4+840*ep^5*n5+968*ep^5*n5*x5+354*ep^5*n5*x5^2+40*ep^5*n5*x5^3-120*
         ep^5*n5^2-108*ep^5*n5^2*x5-24*ep^5*n5^2*x5^2+2580*ep^5*n2+3656*ep^5*n2
         *x5+1829*ep^5*n2*x5^2+371*ep^5*n2*x5^3+24*ep^5*n2*x5^4-1360*ep^5*n2*n5
         -1644*ep^5*n2*n5*x5-642*ep^5*n2*n5*x5^2-80*ep^5*n2*n5*x5^3+120*ep^5*n2
         *n5^2+108*ep^5*n2*n5^2*x5+24*ep^5*n2*n5^2*x5^2-1696*ep^5*n2^2-2572*
         ep^5*n2^2*x5-1398*ep^5*n2^2*x5^2-316*ep^5*n2^2*x5^3-24*ep^5*n2^2*x5^4+
         520*ep^5*n2^2*n5+676*ep^5*n2^2*n5*x5+288*ep^5*n2^2*n5*x5^2+40*ep^5*
         n2^2*n5*x5^3+340*ep^5*n2^3+560*ep^5*n2^3*x5+337*ep^5*n2^3*x5^2+87*ep^5
         *n2^3*x5^3+8*ep^5*n2^3*x5^4+528*ep^6+840*ep^6*x5+484*ep^6*x5^2+118*
         ep^6*x5^3+10*ep^6*x5^4-176*ep^6*n5-240*ep^6*n5*x5-108*ep^6*n5*x5^2-16*
         ep^6*n5*x5^3-824*ep^6*n2-1360*ep^6*n2*x5-822*ep^6*n2*x5^2-214*ep^6*n2*
         x5^3-20*ep^6*n2*x5^4+176*ep^6*n2*n5+240*ep^6*n2*n5*x5+108*ep^6*n2*n5*
         x5^2+16*ep^6*n2*n5*x5^3+296*ep^6*n2^2+520*ep^6*n2^2*x5+338*ep^6*n2^2*
         x5^2+96*ep^6*n2^2*x5^3+10*ep^6*n2^2*x5^4-96*ep^7-176*ep^7*x5-120*ep^7*
         x5^2-36*ep^7*x5^3-4*ep^7*x5^4+96*ep^7*n2+176*ep^7*n2*x5+120*ep^7*n2*
         x5^2+36*ep^7*n2*x5^3+4*ep^7*n2*x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,0,1,1+n5,1,0,0,0)*(rat(8*n5-5*n5^2+n5^3-14*
         n2*n5+4*n2*n5^2+9*n2^2*n5-n2^2*n5^2-2*n2^3*n5+8*ep*x5-26*ep*n5-10*ep*
         n5*x5+11*ep*n5^2+3*ep*n5^2*x5-ep*n5^3-14*ep*n2*x5+27*ep*n2*n5+8*ep*n2*
         n5*x5-3*ep*n2*n5^2+9*ep*n2^2*x5-8*ep*n2^2*n5-2*ep*n2^2*n5*x5-2*ep*n2^3
         *x5-26*ep^2*x5-5*ep^2*x5^2+28*ep^2*n5+22*ep^2*n5*x5+3*ep^2*n5*x5^2-6*
         ep^2*n5^2-3*ep^2*n5^2*x5+27*ep^2*n2*x5+4*ep^2*n2*x5^2-12*ep^2*n2*n5-6*
         ep^2*n2*n5*x5-8*ep^2*n2^2*x5-ep^2*n2^2*x5^2+28*ep^3*x5+11*ep^3*x5^2+
         ep^3*x5^3-10*ep^3*n5-12*ep^3*n5*x5-3*ep^3*n5*x5^2-12*ep^3*n2*x5-3*ep^3
         *n2*x5^2-10*ep^4*x5-6*ep^4*x5^2-ep^4*x5^3,1-n5-3*n2+3*n2*n5+3*n2^2-3*
         n2^2*n5-n2^3+n2^3*n5-6*ep-ep*x5+4*ep*n5+14*ep*n2+3*ep*n2*x5-8*ep*n2*n5
         -10*ep*n2^2-3*ep*n2^2*x5+4*ep*n2^2*n5+2*ep*n2^3+ep*n2^3*x5+13*ep^2+4*
         ep^2*x5-5*ep^2*n5-21*ep^2*n2-8*ep^2*n2*x5+5*ep^2*n2*n5+8*ep^2*n2^2+4*
         ep^2*n2^2*x5-12*ep^3-5*ep^3*x5+2*ep^3*n5+10*ep^3*n2+5*ep^3*n2*x5+4*
         ep^4+2*ep^4*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,0,1,n5,1,0,1,0)*(rat(78-128*n5+76*n5^2-20*
         n5^3+2*n5^4-197*n2+262*n2*n5-116*n2*n5^2+20*n2*n5^3-n2*n5^4+195*n2^2-
         199*n2^2*n5+59*n2^2*n5^2-5*n2^2*n5^3-94*n2^3+66*n2^3*n5-10*n2^3*n5^2+
         22*n2^4-8*n2^4*n5-2*n2^5-462*ep-128*ep*x5+596*ep*n5+152*ep*n5*x5-260*
         ep*n5^2-60*ep*n5^2*x5+44*ep*n5^3+8*ep*n5^3*x5-2*ep*n5^4+932*ep*n2+262*
         ep*n2*x5-912*ep*n2*n5-232*ep*n2*n5*x5+264*ep*n2*n5^2+60*ep*n2*n5^2*x5-
         22*ep*n2*n5^3-4*ep*n2*n5^3*x5-691*ep*n2^2-199*ep*n2^2*x5+460*ep*n2^2*
         n5+118*ep*n2^2*n5*x5-67*ep*n2^2*n5^2-15*ep*n2^2*n5^2*x5+222*ep*n2^3+66
         *ep*n2^3*x5-76*ep*n2^3*n5-20*ep*n2^3*n5*x5-26*ep*n2^4-8*ep*n2^4*x5+
         1086*ep^2+596*ep^2*x5+76*ep^2*x5^2-1032*ep^2*n5-520*ep^2*n5*x5-60*ep^2
         *n5*x5^2+294*ep^2*n5^2+132*ep^2*n5^2*x5+12*ep^2*n5^2*x5^2-24*ep^2*n5^3
         -8*ep^2*n5^3*x5-1641*ep^2*n2-912*ep^2*n2*x5-116*ep^2*n2*x5^2+1050*ep^2
         *n2*n5+528*ep^2*n2*n5*x5+60*ep^2*n2*n5*x5^2-149*ep^2*n2*n5^2-66*ep^2*
         n2*n5^2*x5-6*ep^2*n2*n5^2*x5^2+810*ep^2*n2^2+460*ep^2*n2^2*x5+59*ep^2*
         n2^2*x5^2-264*ep^2*n2^2*n5-134*ep^2*n2^2*n5*x5-15*ep^2*n2^2*n5*x5^2-
         130*ep^2*n2^3-76*ep^2*n2^3*x5-10*ep^2*n2^3*x5^2-1266*ep^3-1032*ep^3*x5
         -260*ep^3*x5^2-20*ep^3*x5^3+788*ep^3*n5+588*ep^3*n5*x5+132*ep^3*n5*
         x5^2+8*ep^3*n5*x5^3-110*ep^3*n5^2-72*ep^3*n5^2*x5-12*ep^3*n5^2*x5^2+
         1274*ep^3*n2+1050*ep^3*n2*x5+264*ep^3*n2*x5^2+20*ep^3*n2*x5^3-400*ep^3
         *n2*n5-298*ep^3*n2*n5*x5-66*ep^3*n2*n5*x5^2-4*ep^3*n2*n5*x5^3-314*ep^3
         *n2^2-264*ep^3*n2^2*x5-67*ep^3*n2^2*x5^2-5*ep^3*n2^2*x5^3+732*ep^4+788
         *ep^4*x5+294*ep^4*x5^2+44*ep^4*x5^3+2*ep^4*x5^4-224*ep^4*n5-220*ep^4*
         n5*x5-72*ep^4*n5*x5^2-8*ep^4*n5*x5^3-368*ep^4*n2-400*ep^4*n2*x5-149*
         ep^4*n2*x5^2-22*ep^4*n2*x5^3-ep^4*n2*x5^4-168*ep^5-224*ep^5*x5-110*
         ep^5*x5^2-24*ep^5*x5^3-2*ep^5*x5^4,-2+3*n5-n5^2+7*n2-10*n2*n5+3*n2*
         n5^2-9*n2^2+12*n2^2*n5-3*n2^2*n5^2+5*n2^3-6*n2^3*n5+n2^3*n5^2-n2^4+
         n2^4*n5+14*ep+3*ep*x5-16*ep*n5-2*ep*n5*x5+4*ep*n5^2-40*ep*n2-10*ep*n2*
         x5+40*ep*n2*n5+6*ep*n2*n5*x5-8*ep*n2*n5^2+40*ep*n2^2+12*ep*n2^2*x5-32*
         ep*n2^2*n5-6*ep*n2^2*n5*x5+4*ep*n2^2*n5^2-16*ep*n2^3-6*ep*n2^3*x5+8*ep
         *n2^3*n5+2*ep*n2^3*n5*x5+2*ep*n2^4+ep*n2^4*x5-38*ep^2-16*ep^2*x5-ep^2*
         x5^2+31*ep^2*n5+8*ep^2*n5*x5-5*ep^2*n5^2+83*ep^2*n2+40*ep^2*n2*x5+3*
         ep^2*n2*x5^2-52*ep^2*n2*n5-16*ep^2*n2*n5*x5+5*ep^2*n2*n5^2-57*ep^2*
         n2^2-32*ep^2*n2^2*x5-3*ep^2*n2^2*x5^2+21*ep^2*n2^2*n5+8*ep^2*n2^2*n5*
         x5+12*ep^2*n2^3+8*ep^2*n2^3*x5+ep^2*n2^3*x5^2+50*ep^3+31*ep^3*x5+4*
         ep^3*x5^2-26*ep^3*n5-10*ep^3*n5*x5+2*ep^3*n5^2-74*ep^3*n2-52*ep^3*n2*
         x5-8*ep^3*n2*x5^2+22*ep^3*n2*n5+10*ep^3*n2*n5*x5+26*ep^3*n2^2+21*ep^3*
         n2^2*x5+4*ep^3*n2^2*x5^2-32*ep^4-26*ep^4*x5-5*ep^4*x5^2+8*ep^4*n5+4*
         ep^4*n5*x5+24*ep^4*n2+22*ep^4*n2*x5+5*ep^4*n2*x5^2+8*ep^5+8*ep^5*x5+2*
         ep^5*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,0,1+n5,0,0,1,0)*(rat(84*n5-76*n5^2+22*n5^3
         -2*n5^4-160*n2*n5+108*n2*n5^2-21*n2*n5^3+n2*n5^4+113*n2^2*n5-51*n2^2*
         n5^2+5*n2^2*n5^3-35*n2^3*n5+8*n2^3*n5^2+4*n2^4*n5+84*ep*x5-378*ep*n5-
         152*ep*n5*x5+254*ep*n5^2+66*ep*n5^2*x5-48*ep*n5^3-8*ep*n5^3*x5+2*ep*
         n5^4-160*ep*n2*x5+539*ep*n2*n5+216*ep*n2*n5*x5-241*ep*n2*n5^2-63*ep*n2
         *n5^2*x5+23*ep*n2*n5^3+4*ep*n2*n5^3*x5+113*ep*n2^2*x5-253*ep*n2^2*n5-
         102*ep*n2^2*n5*x5+57*ep*n2^2*n5^2+15*ep*n2^2*n5^2*x5-35*ep*n2^3*x5+39*
         ep*n2^3*n5+16*ep*n2^3*n5*x5+4*ep*n2^4*x5-378*ep^2*x5-76*ep^2*x5^2+636*
         ep^2*n5+508*ep^2*n5*x5+66*ep^2*n5*x5^2-282*ep^2*n5^2-144*ep^2*n5^2*x5-
         12*ep^2*n5^2*x5^2+26*ep^2*n5^3+8*ep^2*n5^3*x5+539*ep^2*n2*x5+108*ep^2*
         n2*x5^2-603*ep^2*n2*n5-482*ep^2*n2*n5*x5-63*ep^2*n2*n5*x5^2+134*ep^2*
         n2*n5^2+69*ep^2*n2*n5^2*x5+6*ep^2*n2*n5^2*x5^2-253*ep^2*n2^2*x5-51*
         ep^2*n2^2*x5^2+141*ep^2*n2^2*n5+114*ep^2*n2^2*n5*x5+15*ep^2*n2^2*n5*
         x5^2+39*ep^2*n2^3*x5+8*ep^2*n2^3*x5^2+636*ep^3*x5+254*ep^3*x5^2+22*
         ep^3*x5^3-474*ep^3*n5-564*ep^3*n5*x5-144*ep^3*n5*x5^2-8*ep^3*n5*x5^3+
         104*ep^3*n5^2+78*ep^3*n5^2*x5+12*ep^3*n5^2*x5^2-603*ep^3*n2*x5-241*
         ep^3*n2*x5^2-21*ep^3*n2*x5^3+224*ep^3*n2*n5+268*ep^3*n2*n5*x5+69*ep^3*
         n2*n5*x5^2+4*ep^3*n2*n5*x5^3+141*ep^3*n2^2*x5+57*ep^3*n2^2*x5^2+5*ep^3
         *n2^2*x5^3-474*ep^4*x5-282*ep^4*x5^2-48*ep^4*x5^3-2*ep^4*x5^4+132*ep^4
         *n5+208*ep^4*n5*x5+78*ep^4*n5*x5^2+8*ep^4*n5*x5^3+224*ep^4*n2*x5+134*
         ep^4*n2*x5^2+23*ep^4*n2*x5^3+ep^4*n2*x5^4+132*ep^5*x5+104*ep^5*x5^2+26
         *ep^5*x5^3+2*ep^5*x5^4,-16+12*n5-2*n5^2+56*n2-36*n2*n5+5*n2*n5^2-76*
         n2^2+39*n2^2*n5-4*n2^2*n5^2+50*n2^3-18*n2^3*n5+n2^3*n5^2-16*n2^4+3*
         n2^4*n5+2*n2^5+104*ep+12*ep*x5-66*ep*n5-4*ep*n5*x5+9*ep*n5^2-288*ep*n2
         -36*ep*n2*x5+147*ep*n2*n5+10*ep*n2*n5*x5-15*ep*n2*n5^2+290*ep*n2^2+39*
         ep*n2^2*x5-105*ep*n2^2*n5-8*ep*n2^2*n5*x5+6*ep*n2^2*n5^2-126*ep*n2^3-
         18*ep*n2^3*x5+24*ep*n2^3*n5+2*ep*n2^3*n5*x5+20*ep*n2^4+3*ep*n2^4*x5-
         264*ep^2-66*ep^2*x5-2*ep^2*x5^2+132*ep^2*n5+18*ep^2*n5*x5-13*ep^2*n5^2
         +544*ep^2*n2+147*ep^2*n2*x5+5*ep^2*n2*x5^2-195*ep^2*n2*n5-30*ep^2*n2*
         n5*x5+11*ep^2*n2*n5^2-362*ep^2*n2^2-105*ep^2*n2^2*x5-4*ep^2*n2^2*x5^2+
         69*ep^2*n2^2*n5+12*ep^2*n2^2*n5*x5+78*ep^2*n2^3+24*ep^2*n2^3*x5+ep^2*
         n2^3*x5^2+328*ep^3+132*ep^3*x5+9*ep^3*x5^2-114*ep^3*n5-26*ep^3*n5*x5+6
         *ep^3*n5^2-448*ep^3*n2-195*ep^3*n2*x5-15*ep^3*n2*x5^2+84*ep^3*n2*n5+22
         *ep^3*n2*n5*x5+148*ep^3*n2^2+69*ep^3*n2^2*x5+6*ep^3*n2^2*x5^2-200*ep^4
         -114*ep^4*x5-13*ep^4*x5^2+36*ep^4*n5+12*ep^4*n5*x5+136*ep^4*n2+84*ep^4
         *n2*x5+11*ep^4*n2*x5^2+48*ep^5+36*ep^5*x5+6*ep^5*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,0,1+n5,1,0,0,0)*(rat(-84*n5+76*n5^2-22*
         n5^3+2*n5^4+160*n2*n5-108*n2*n5^2+21*n2*n5^3-n2*n5^4-113*n2^2*n5+51*
         n2^2*n5^2-5*n2^2*n5^3+35*n2^3*n5-8*n2^3*n5^2-4*n2^4*n5-84*ep*x5+378*ep
         *n5+152*ep*n5*x5-254*ep*n5^2-66*ep*n5^2*x5+48*ep*n5^3+8*ep*n5^3*x5-2*
         ep*n5^4+160*ep*n2*x5-539*ep*n2*n5-216*ep*n2*n5*x5+241*ep*n2*n5^2+63*ep
         *n2*n5^2*x5-23*ep*n2*n5^3-4*ep*n2*n5^3*x5-113*ep*n2^2*x5+253*ep*n2^2*
         n5+102*ep*n2^2*n5*x5-57*ep*n2^2*n5^2-15*ep*n2^2*n5^2*x5+35*ep*n2^3*x5-
         39*ep*n2^3*n5-16*ep*n2^3*n5*x5-4*ep*n2^4*x5+378*ep^2*x5+76*ep^2*x5^2-
         636*ep^2*n5-508*ep^2*n5*x5-66*ep^2*n5*x5^2+282*ep^2*n5^2+144*ep^2*n5^2
         *x5+12*ep^2*n5^2*x5^2-26*ep^2*n5^3-8*ep^2*n5^3*x5-539*ep^2*n2*x5-108*
         ep^2*n2*x5^2+603*ep^2*n2*n5+482*ep^2*n2*n5*x5+63*ep^2*n2*n5*x5^2-134*
         ep^2*n2*n5^2-69*ep^2*n2*n5^2*x5-6*ep^2*n2*n5^2*x5^2+253*ep^2*n2^2*x5+
         51*ep^2*n2^2*x5^2-141*ep^2*n2^2*n5-114*ep^2*n2^2*n5*x5-15*ep^2*n2^2*n5
         *x5^2-39*ep^2*n2^3*x5-8*ep^2*n2^3*x5^2-636*ep^3*x5-254*ep^3*x5^2-22*
         ep^3*x5^3+474*ep^3*n5+564*ep^3*n5*x5+144*ep^3*n5*x5^2+8*ep^3*n5*x5^3-
         104*ep^3*n5^2-78*ep^3*n5^2*x5-12*ep^3*n5^2*x5^2+603*ep^3*n2*x5+241*
         ep^3*n2*x5^2+21*ep^3*n2*x5^3-224*ep^3*n2*n5-268*ep^3*n2*n5*x5-69*ep^3*
         n2*n5*x5^2-4*ep^3*n2*n5*x5^3-141*ep^3*n2^2*x5-57*ep^3*n2^2*x5^2-5*ep^3
         *n2^2*x5^3+474*ep^4*x5+282*ep^4*x5^2+48*ep^4*x5^3+2*ep^4*x5^4-132*ep^4
         *n5-208*ep^4*n5*x5-78*ep^4*n5*x5^2-8*ep^4*n5*x5^3-224*ep^4*n2*x5-134*
         ep^4*n2*x5^2-23*ep^4*n2*x5^3-ep^4*n2*x5^4-132*ep^5*x5-104*ep^5*x5^2-26
         *ep^5*x5^3-2*ep^5*x5^4,-16+12*n5-2*n5^2+56*n2-36*n2*n5+5*n2*n5^2-76*
         n2^2+39*n2^2*n5-4*n2^2*n5^2+50*n2^3-18*n2^3*n5+n2^3*n5^2-16*n2^4+3*
         n2^4*n5+2*n2^5+104*ep+12*ep*x5-66*ep*n5-4*ep*n5*x5+9*ep*n5^2-288*ep*n2
         -36*ep*n2*x5+147*ep*n2*n5+10*ep*n2*n5*x5-15*ep*n2*n5^2+290*ep*n2^2+39*
         ep*n2^2*x5-105*ep*n2^2*n5-8*ep*n2^2*n5*x5+6*ep*n2^2*n5^2-126*ep*n2^3-
         18*ep*n2^3*x5+24*ep*n2^3*n5+2*ep*n2^3*n5*x5+20*ep*n2^4+3*ep*n2^4*x5-
         264*ep^2-66*ep^2*x5-2*ep^2*x5^2+132*ep^2*n5+18*ep^2*n5*x5-13*ep^2*n5^2
         +544*ep^2*n2+147*ep^2*n2*x5+5*ep^2*n2*x5^2-195*ep^2*n2*n5-30*ep^2*n2*
         n5*x5+11*ep^2*n2*n5^2-362*ep^2*n2^2-105*ep^2*n2^2*x5-4*ep^2*n2^2*x5^2+
         69*ep^2*n2^2*n5+12*ep^2*n2^2*n5*x5+78*ep^2*n2^3+24*ep^2*n2^3*x5+ep^2*
         n2^3*x5^2+328*ep^3+132*ep^3*x5+9*ep^3*x5^2-114*ep^3*n5-26*ep^3*n5*x5+6
         *ep^3*n5^2-448*ep^3*n2-195*ep^3*n2*x5-15*ep^3*n2*x5^2+84*ep^3*n2*n5+22
         *ep^3*n2*n5*x5+148*ep^3*n2^2+69*ep^3*n2^2*x5+6*ep^3*n2^2*x5^2-200*ep^4
         -114*ep^4*x5-13*ep^4*x5^2+36*ep^4*n5+12*ep^4*n5*x5+136*ep^4*n2+84*ep^4
         *n2*x5+11*ep^4*n2*x5^2+48*ep^5+36*ep^5*x5+6*ep^5*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,0,n5,1,1,0,0)*(rat(20-24*n5+9*n5^2-n5^3-26
         *n2+29*n2*n5-10*n2*n5^2+n2*n5^3+8*n2^2-8*n2^2*n5+2*n2^2*n5^2-84*ep-24*
         ep*x5+76*ep*n5+18*ep*n5*x5-19*ep*n5^2-3*ep*n5^2*x5+ep*n5^3+74*ep*n2+29
         *ep*n2*x5-56*ep*n2*n5-20*ep*n2*n5*x5+10*ep*n2*n5^2+3*ep*n2*n5^2*x5-12*
         ep*n2^2-8*ep*n2^2*x5+6*ep*n2^2*n5+4*ep*n2^2*n5*x5+132*ep^2+76*ep^2*x5+
         9*ep^2*x5^2-80*ep^2*n5-38*ep^2*n5*x5-3*ep^2*n5*x5^2+10*ep^2*n5^2+3*
         ep^2*n5^2*x5-68*ep^2*n2-56*ep^2*n2*x5-10*ep^2*n2*x5^2+26*ep^2*n2*n5+20
         *ep^2*n2*n5*x5+3*ep^2*n2*n5*x5^2+4*ep^2*n2^2+6*ep^2*n2^2*x5+2*ep^2*
         n2^2*x5^2-92*ep^3-80*ep^3*x5-19*ep^3*x5^2-ep^3*x5^3+28*ep^3*n5+20*ep^3
         *n5*x5+3*ep^3*n5*x5^2+20*ep^3*n2+26*ep^3*n2*x5+10*ep^3*n2*x5^2+ep^3*n2
         *x5^3+24*ep^4+28*ep^4*x5+10*ep^4*x5^2+ep^4*x5^3,4-2*n5-12*n2+5*n2*n5+
         13*n2^2-4*n2^2*n5-6*n2^3+n2^3*n5+n2^4-22*ep-2*ep*x5+9*ep*n5+49*ep*n2+5
         *ep*n2*x5-15*ep*n2*n5-35*ep*n2^2-4*ep*n2^2*x5+6*ep*n2^2*n5+8*ep*n2^3+
         ep*n2^3*x5+44*ep^2+9*ep^2*x5-13*ep^2*n5-65*ep^2*n2-15*ep^2*n2*x5+11*
         ep^2*n2*n5+23*ep^2*n2^2+6*ep^2*n2^2*x5-38*ep^3-13*ep^3*x5+6*ep^3*n5+28
         *ep^3*n2+11*ep^3*n2*x5+12*ep^4+6*ep^4*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,n2,1,1,n5,1,0,0,0)*(rat(1956-4138*n5+3652*n5^2-
         1724*n5^3+460*n5^4-66*n5^5+4*n5^6-6256*n2+11319*n2*n5-8288*n2*n5^2+
         3120*n2*n5^3-626*n2*n5^4+61*n2*n5^5-2*n2*n5^6+8489*n2^2-12815*n2^2*n5+
         7483*n2^2*n5^2-2105*n2^2*n5^3+282*n2^2*n5^4-14*n2^2*n5^5-6323*n2^3+
         7683*n2^3*n5-3362*n2^3*n5^2+628*n2^3*n5^3-42*n2^3*n5^4+2785*n2^4-2570*
         n2^4*n5+752*n2^4*n5^2-70*n2^4*n5^3-723*n2^5+454*n2^5*n5-67*n2^5*n5^2+
         102*n2^6-33*n2^6*n5-6*n2^7-14202*ep-4138*ep*x5+26003*ep*n5+7304*ep*n5*
         x5-19384*ep*n5^2-5172*ep*n5^2*x5+7480*ep*n5^3+1840*ep*n5^3*x5-1556*ep*
         n5^4-330*ep*n5^4*x5+161*ep*n5^5+24*ep*n5^5*x5-6*ep*n5^6+38713*ep*n2+
         11319*ep*n2*x5-59009*ep*n2*n5-16576*ep*n2*n5*x5+35062*ep*n2*n5^2+9360*
         ep*n2*n5^2*x5-10110*ep*n2*n5^3-2504*ep*n2*n5^3*x5+1401*ep*n2*n5^4+305*
         ep*n2*n5^4*x5-73*ep*n2*n5^5-12*ep*n2*n5^5*x5-43488*ep*n2^2-12815*ep*
         n2^2*x5+53158*ep*n2^2*n5+14966*ep*n2^2*n5*x5-23637*ep*n2^2*n5^2-6315*
         ep*n2^2*n5^2*x5+4529*ep*n2^2*n5^3+1128*ep*n2^2*n5^3*x5-314*ep*n2^2*
         n5^4-70*ep*n2^2*n5^4*x5+25723*ep*n2^3+7683*ep*n2^3*x5-23751*ep*n2^3*n5
         -6724*ep*n2^3*n5*x5+7040*ep*n2^3*n5^2+1884*ep*n2^3*n5^2*x5-672*ep*n2^3
         *n5^3-168*ep*n2^3*n5^3*x5-8429*ep*n2^4-2570*ep*n2^4*x5+5259*ep*n2^4*n5
         +1504*ep*n2^4*n5*x5-782*ep*n2^4*n5^2-210*ep*n2^4*n5^2*x5+1446*ep*n2^5+
         454*ep*n2^5*x5-461*ep*n2^5*n5-134*ep*n2^5*n5*x5-101*ep*n2^6-33*ep*n2^6
         *x5+43250*ep^2+26003*ep^2*x5+3652*ep^2*x5^2-66552*ep^2*n5-38768*ep^2*
         n5*x5-5172*ep^2*n5*x5^2+40159*ep^2*n5^2+22440*ep^2*n5^2*x5+2760*ep^2*
         n5^2*x5^2-11839*ep^2*n5^3-6224*ep^2*n5^3*x5-660*ep^2*n5^3*x5^2+1695*
         ep^2*n5^4+805*ep^2*n5^4*x5+60*ep^2*n5^4*x5^2-93*ep^2*n5^5-36*ep^2*n5^5
         *x5-97665*ep^2*n2-59009*ep^2*n2*x5-8288*ep^2*n2*x5^2+120192*ep^2*n2*n5
         +70124*ep^2*n2*n5*x5+9360*ep^2*n2*n5*x5^2-54198*ep^2*n2*n5^2-30330*
         ep^2*n2*n5^2*x5-3756*ep^2*n2*n5^2*x5^2+10598*ep^2*n2*n5^3+5604*ep^2*n2
         *n5^3*x5+610*ep^2*n2*n5^3*x5^2-755*ep^2*n2*n5^4-365*ep^2*n2*n5^4*x5-30
         *ep^2*n2*n5^4*x5^2+87205*ep^2*n2^2+53158*ep^2*n2^2*x5+7483*ep^2*n2^2*
         x5^2-80746*ep^2*n2^2*n5-47274*ep^2*n2^2*n5*x5-6315*ep^2*n2^2*n5*x5^2+
         24228*ep^2*n2^2*n5^2+13587*ep^2*n2^2*n5^2*x5+1692*ep^2*n2^2*n5^2*x5^2-
         2360*ep^2*n2^2*n5^3-1256*ep^2*n2^2*n5^3*x5-140*ep^2*n2^2*n5^3*x5^2-
         38421*ep^2*n2^3-23751*ep^2*n2^3*x5-3362*ep^2*n2^3*x5^2+23906*ep^2*n2^3
         *n5+14080*ep^2*n2^3*n5*x5+1884*ep^2*n2^3*n5*x5^2-3587*ep^2*n2^3*n5^2-
         2016*ep^2*n2^3*n5^2*x5-252*ep^2*n2^3*n5^2*x5^2+8332*ep^2*n2^4+5259*
         ep^2*n2^4*x5+752*ep^2*n2^4*x5^2-2630*ep^2*n2^4*n5-1564*ep^2*n2^4*n5*x5
         -210*ep^2*n2^4*n5*x5^2-709*ep^2*n2^5-461*ep^2*n2^5*x5-67*ep^2*n2^5*
         x5^2-71458*ep^3-66552*ep^3*x5-19384*ep^3*x5^2-1724*ep^3*x5^3+88627*
         ep^3*n5+80318*ep^3*n5*x5+22440*ep^3*n5*x5^2+1840*ep^3*n5*x5^3-40539*
         ep^3*n5^2-35517*ep^3*n5^2*x5-9336*ep^3*n5^2*x5^2-660*ep^3*n5^2*x5^3+
         8103*ep^3*n5^3+6780*ep^3*n5^3*x5+1610*ep^3*n5^3*x5^2+80*ep^3*n5^3*x5^3
         -597*ep^3*n5^4-465*ep^3*n5^4*x5-90*ep^3*n5^4*x5^2+128271*ep^3*n2+
         120192*ep^3*n2*x5+35062*ep^3*n2*x5^2+3120*ep^3*n2*x5^3-119295*ep^3*n2*
         n5-108396*ep^3*n2*n5*x5-30330*ep^3*n2*n5*x5^2-2504*ep^3*n2*n5*x5^3+
         36216*ep^3*n2*n5^2+31794*ep^3*n2*n5^2*x5+8406*ep^3*n2*n5^2*x5^2+610*
         ep^3*n2*n5^2*x5^3-3592*ep^3*n2*n5^3-3020*ep^3*n2*n5^3*x5-730*ep^3*n2*
         n5^3*x5^2-40*ep^3*n2*n5^3*x5^3-85350*ep^3*n2^2-80746*ep^3*n2^2*x5-
         23637*ep^3*n2^2*x5^2-2105*ep^3*n2^2*x5^3+53104*ep^3*n2^2*n5+48456*ep^3
         *n2^2*n5*x5+13587*ep^3*n2^2*n5*x5^2+1128*ep^3*n2^2*n5*x5^3-8042*ep^3*
         n2^2*n5^2-7080*ep^3*n2^2*n5^2*x5-1884*ep^3*n2^2*n5^2*x5^2-140*ep^3*
         n2^2*n5^2*x5^3+24908*ep^3*n2^3+23906*ep^3*n2^3*x5+7040*ep^3*n2^3*x5^2+
         628*ep^3*n2^3*x5^3-7815*ep^3*n2^3*n5-7174*ep^3*n2^3*n5*x5-2016*ep^3*
         n2^3*n5*x5^2-168*ep^3*n2^3*n5*x5^3-2683*ep^3*n2^4-2630*ep^3*n2^4*x5-
         782*ep^3*n2^4*x5^2-70*ep^3*n2^4*x5^3+68946*ep^4+88627*ep^4*x5+40159*
         ep^4*x5^2+7480*ep^4*x5^3+460*ep^4*x5^4-64568*ep^4*n5-81078*ep^4*n5*x5-
         35517*ep^4*n5*x5^2-6224*ep^4*n5*x5^3-330*ep^4*n5*x5^4+19892*ep^4*n5^2+
         24309*ep^4*n5^2*x5+10170*ep^4*n5^2*x5^2+1610*ep^4*n5^2*x5^3+60*ep^4*
         n5^2*x5^4-2022*ep^4*n5^3-2388*ep^4*n5^3*x5-930*ep^4*n5^3*x5^2-120*ep^4
         *n5^3*x5^3-92155*ep^4*n2-119295*ep^4*n2*x5-54198*ep^4*n2*x5^2-10110*
         ep^4*n2*x5^3-626*ep^4*n2*x5^4+57487*ep^4*n2*n5+72432*ep^4*n2*n5*x5+
         31794*ep^4*n2*n5*x5^2+5604*ep^4*n2*n5*x5^3+305*ep^4*n2*n5*x5^4-8798*
         ep^4*n2*n5^2-10776*ep^4*n2*n5^2*x5-4530*ep^4*n2*n5^2*x5^2-730*ep^4*n2*
         n5^2*x5^3-30*ep^4*n2*n5^2*x5^4+40608*ep^4*n2^2+53104*ep^4*n2^2*x5+
         24228*ep^4*n2^2*x5^2+4529*ep^4*n2^2*x5^3+282*ep^4*n2^2*x5^4-12707*ep^4
         *n2^2*n5-16084*ep^4*n2^2*n5*x5-7080*ep^4*n2^2*n5*x5^2-1256*ep^4*n2^2*
         n5*x5^3-70*ep^4*n2^2*n5*x5^4-5889*ep^4*n2^3-7815*ep^4*n2^3*x5-3587*
         ep^4*n2^3*x5^2-672*ep^4*n2^3*x5^3-42*ep^4*n2^3*x5^4-38636*ep^5-64568*
         ep^5*x5-40539*ep^5*x5^2-11839*ep^5*x5^3-1556*ep^5*x5^4-66*ep^5*x5^5+
         24276*ep^5*n5+39784*ep^5*n5*x5+24309*ep^5*n5*x5^2+6780*ep^5*n5*x5^3+
         805*ep^5*n5*x5^4+24*ep^5*n5*x5^5-3780*ep^5*n5^2-6066*ep^5*n5^2*x5-3582
         *ep^5*n5^2*x5^2-930*ep^5*n5^2*x5^3-90*ep^5*n5^2*x5^4+34128*ep^5*n2+
         57487*ep^5*n2*x5+36216*ep^5*n2*x5^2+10598*ep^5*n2*x5^3+1401*ep^5*n2*
         x5^4+61*ep^5*n2*x5^5-10694*ep^5*n2*n5-17596*ep^5*n2*n5*x5-10776*ep^5*
         n2*n5*x5^2-3020*ep^5*n2*n5*x5^3-365*ep^5*n2*n5*x5^4-12*ep^5*n2*n5*x5^5
         -7464*ep^5*n2^2-12707*ep^5*n2^2*x5-8042*ep^5*n2^2*x5^2-2360*ep^5*n2^2*
         x5^3-314*ep^5*n2^2*x5^4-14*ep^5*n2^2*x5^5+11536*ep^6+24276*ep^6*x5+
         19892*ep^6*x5^2+8103*ep^6*x5^3+1695*ep^6*x5^4+161*ep^6*x5^5+4*ep^6*
         x5^6-3648*ep^6*n5-7560*ep^6*n5*x5-6066*ep^6*n5*x5^2-2388*ep^6*n5*x5^3-
         465*ep^6*n5*x5^4-36*ep^6*n5*x5^5-5036*ep^6*n2-10694*ep^6*n2*x5-8798*
         ep^6*n2*x5^2-3592*ep^6*n2*x5^3-755*ep^6*n2*x5^4-73*ep^6*n2*x5^5-2*ep^6
         *n2*x5^6-1392*ep^7-3648*ep^7*x5-3780*ep^7*x5^2-2022*ep^7*x5^3-597*ep^7
         *x5^4-93*ep^7*x5^5-6*ep^7*x5^6,-48+64*n5-28*n5^2+4*n5^3+220*n2-270*n2*
         n5+108*n2*n5^2-14*n2*n5^3-416*n2^2+457*n2^2*n5-161*n2^2*n5^2+18*n2^2*
         n5^3+419*n2^3-395*n2^3*n5+115*n2^3*n5^2-10*n2^3*n5^3-242*n2^4+183*n2^4
         *n5-39*n2^4*n5^2+2*n2^4*n5^3+80*n2^5-43*n2^5*n5+5*n2^5*n5^2-14*n2^6+4*
         n2^6*n5+n2^7+408*ep+64*ep*x5-480*ep*n5-56*ep*n5*x5+182*ep*n5^2+12*ep*
         n5^2*x5-22*ep*n5^3-1590*ep*n2-270*ep*n2*x5+1675*ep*n2*n5+216*ep*n2*n5*
         x5-559*ep*n2*n5^2-42*ep*n2*n5^2*x5+58*ep*n2*n5^3+2483*ep*n2^2+457*ep*
         n2^2*x5-2248*ep*n2^2*n5-322*ep*n2^2*n5*x5+621*ep*n2^2*n5^2+54*ep*n2^2*
         n5^2*x5-50*ep*n2^2*n5^3-1981*ep*n2^3-395*ep*n2^3*x5+1442*ep*n2^3*n5+
         230*ep*n2^3*n5*x5-293*ep*n2^3*n5^2-30*ep*n2^3*n5^2*x5+14*ep*n2^3*n5^3+
         849*ep*n2^4+183*ep*n2^4*x5-440*ep*n2^4*n5-78*ep*n2^4*n5*x5+49*ep*n2^4*
         n5^2+6*ep*n2^4*n5^2*x5-185*ep*n2^5-43*ep*n2^5*x5+51*ep*n2^5*n5+10*ep*
         n2^5*n5*x5+16*ep*n2^6+4*ep*n2^6*x5-1464*ep^2-480*ep^2*x5-28*ep^2*x5^2+
         1472*ep^2*n5+364*ep^2*n5*x5+12*ep^2*n5*x5^2-462*ep^2*n5^2-66*ep^2*n5^2
         *x5+44*ep^2*n5^3+4724*ep^2*n2+1675*ep^2*n2*x5+108*ep^2*n2*x5^2-4087*
         ep^2*n2*n5-1118*ep^2*n2*n5*x5-42*ep^2*n2*n5*x5^2+1062*ep^2*n2*n5^2+174
         *ep^2*n2*n5^2*x5-78*ep^2*n2*n5^3-5857*ep^2*n2^2-2248*ep^2*n2^2*x5-161*
         ep^2*n2^2*x5^2+4084*ep^2*n2^2*n5+1242*ep^2*n2^2*n5*x5+54*ep^2*n2^2*n5*
         x5^2-783*ep^2*n2^2*n5^2-150*ep^2*n2^2*n5^2*x5+34*ep^2*n2^2*n5^3+3474*
         ep^2*n2^3+1442*ep^2*n2^3*x5+115*ep^2*n2^3*x5^2-1730*ep^2*n2^3*n5-586*
         ep^2*n2^3*n5*x5-30*ep^2*n2^3*n5*x5^2+183*ep^2*n2^3*n5^2+42*ep^2*n2^3*
         n5^2*x5-983*ep^2*n2^4-440*ep^2*n2^4*x5-39*ep^2*n2^4*x5^2+261*ep^2*n2^4
         *n5+98*ep^2*n2^4*n5*x5+6*ep^2*n2^4*n5*x5^2+106*ep^2*n2^5+51*ep^2*n2^5*
         x5+5*ep^2*n2^5*x5^2+2880*ep^3+1472*ep^3*x5+182*ep^3*x5^2+4*ep^3*x5^3-
         2368*ep^3*n5-924*ep^3*n5*x5-66*ep^3*n5*x5^2+574*ep^3*n5^2+132*ep^3*
         n5^2*x5-38*ep^3*n5^3-7396*ep^3*n2-4087*ep^3*n2*x5-559*ep^3*n2*x5^2-14*
         ep^3*n2*x5^3+4911*ep^3*n2*n5+2124*ep^3*n2*n5*x5+174*ep^3*n2*n5*x5^2-
         879*ep^3*n2*n5^2-234*ep^3*n2*n5^2*x5+34*ep^3*n2*n5^3+6833*ep^3*n2^2+
         4084*ep^3*n2^2*x5+621*ep^3*n2^2*x5^2+18*ep^3*n2^2*x5^3-3252*ep^3*n2^2*
         n5-1566*ep^3*n2^2*n5*x5-150*ep^3*n2^2*n5*x5^2+323*ep^3*n2^2*n5^2+102*
         ep^3*n2^2*n5^2*x5-2681*ep^3*n2^3-1730*ep^3*n2^3*x5-293*ep^3*n2^3*x5^2-
         10*ep^3*n2^3*x5^3+683*ep^3*n2^3*n5+366*ep^3*n2^3*n5*x5+42*ep^3*n2^3*n5
         *x5^2+376*ep^3*n2^4+261*ep^3*n2^4*x5+49*ep^3*n2^4*x5^2+2*ep^3*n2^4*
         x5^3-3360*ep^4-2368*ep^4*x5-462*ep^4*x5^2-22*ep^4*x5^3+2112*ep^4*n5+
         1148*ep^4*n5*x5+132*ep^4*n5*x5^2-350*ep^4*n5^2-114*ep^4*n5^2*x5+12*
         ep^4*n5^3+6444*ep^4*n2+4911*ep^4*n2*x5+1062*ep^4*n2*x5^2+58*ep^4*n2*
         x5^3-2911*ep^4*n2*n5-1758*ep^4*n2*n5*x5-234*ep^4*n2*n5*x5^2+268*ep^4*
         n2*n5^2+102*ep^4*n2*n5^2*x5-3947*ep^4*n2^2-3252*ep^4*n2^2*x5-783*ep^4*
         n2^2*x5^2-50*ep^4*n2^2*x5^3+959*ep^4*n2^2*n5+646*ep^4*n2^2*n5*x5+102*
         ep^4*n2^2*n5*x5^2+769*ep^4*n2^3+683*ep^4*n2^3*x5+183*ep^4*n2^3*x5^2+14
         *ep^4*n2^3*x5^3+2328*ep^5+2112*ep^5*x5+574*ep^5*x5^2+44*ep^5*x5^3-992*
         ep^5*n5-700*ep^5*n5*x5-114*ep^5*n5*x5^2+84*ep^5*n5^2+36*ep^5*n5^2*x5-
         2966*ep^5*n2-2911*ep^5*n2*x5-879*ep^5*n2*x5^2-78*ep^5*n2*x5^3+682*ep^5
         *n2*n5+536*ep^5*n2*n5*x5+102*ep^5*n2*n5*x5^2+904*ep^5*n2^2+959*ep^5*
         n2^2*x5+323*ep^5*n2^2*x5^2+34*ep^5*n2^2*x5^3-888*ep^6-992*ep^6*x5-350*
         ep^6*x5^2-38*ep^6*x5^3+192*ep^6*n5+168*ep^6*n5*x5+36*ep^6*n5*x5^2+564*
         ep^6*n2+682*ep^6*n2*x5+268*ep^6*n2*x5^2+34*ep^6*n2*x5^3+144*ep^7+192*
         ep^7*x5+84*ep^7*x5^2+12*ep^7*x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,2,-1+n2,0,1,n5,0,1,1,0)*(rat(-162+297*n5-216*n5^2
         +78*n5^3-14*n5^4+n5^5+375*n2-569*n2*n5+325*n2*n5^2-83*n2*n5^3+8*n2*
         n5^4-323*n2^2+373*n2^2*n5-145*n2^2*n5^2+19*n2^2*n5^3+132*n2^3-102*n2^3
         *n5+20*n2^3*n5^2-26*n2^4+10*n2^4*n5+2*n2^5+816*ep+297*ep*x5-1169*ep*n5
         -432*ep*n5*x5+623*ep*n5^2+234*ep*n5^2*x5-147*ep*n5^3-56*ep*n5^3*x5+13*
         ep*n5^4+5*ep*n5^4*x5-1573*ep*n2-569*ep*n2*x5+1776*ep*n2*n5+650*ep*n2*
         n5*x5-673*ep*n2*n5^2-249*ep*n2*n5^2*x5+86*ep*n2*n5^3+32*ep*n2*n5^3*x5+
         1035*ep*n2^2+373*ep*n2^2*x5-796*ep*n2^2*n5-290*ep*n2^2*n5*x5+155*ep*
         n2^2*n5^2+57*ep*n2^2*n5^2*x5-284*ep*n2^3-102*ep*n2^3*x5+110*ep*n2^3*n5
         +40*ep*n2^3*n5*x5+28*ep*n2^4+10*ep*n2^4*x5-1562*ep^2-1169*ep^2*x5-216*
         ep^2*x5^2+1624*ep^2*n5+1246*ep^2*n5*x5+234*ep^2*n5*x5^2-558*ep^2*n5^2-
         441*ep^2*n5^2*x5-84*ep^2*n5^2*x5^2+64*ep^2*n5^3+52*ep^2*n5^3*x5+10*
         ep^2*n5^3*x5^2+2402*ep^2*n2+1776*ep^2*n2*x5+325*ep^2*n2*x5^2-1792*ep^2
         *n2*n5-1346*ep^2*n2*n5*x5-249*ep^2*n2*n5*x5^2+338*ep^2*n2*n5^2+258*
         ep^2*n2*n5^2*x5+48*ep^2*n2*n5^2*x5^2-1082*ep^2*n2^2-796*ep^2*n2^2*x5-
         145*ep^2*n2^2*x5^2+416*ep^2*n2^2*n5+310*ep^2*n2^2*n5*x5+57*ep^2*n2^2*
         n5*x5^2+150*ep^2*n2^3+110*ep^2*n2^3*x5+20*ep^2*n2^3*x5^2+1388*ep^3+
         1624*ep^3*x5+623*ep^3*x5^2+78*ep^3*x5^3-912*ep^3*n5-1116*ep^3*n5*x5-
         441*ep^3*n5*x5^2-56*ep^3*n5*x5^3+148*ep^3*n5^2+192*ep^3*n5^2*x5+78*
         ep^3*n5^2*x5^2+10*ep^3*n5^2*x5^3-1572*ep^3*n2-1792*ep^3*n2*x5-673*ep^3
         *n2*x5^2-83*ep^3*n2*x5^3+580*ep^3*n2*n5+676*ep^3*n2*n5*x5+258*ep^3*n2*
         n5*x5^2+32*ep^3*n2*n5*x5^3+368*ep^3*n2^2+416*ep^3*n2^2*x5+155*ep^3*
         n2^2*x5^2+19*ep^3*n2^2*x5^3-544*ep^4-912*ep^4*x5-558*ep^4*x5^2-147*
         ep^4*x5^3-14*ep^4*x5^4+160*ep^4*n5+296*ep^4*n5*x5+192*ep^4*n5*x5^2+52*
         ep^4*n5*x5^3+5*ep^4*n5*x5^4+368*ep^4*n2+580*ep^4*n2*x5+338*ep^4*n2*
         x5^2+86*ep^4*n2*x5^3+8*ep^4*n2*x5^4+64*ep^5+160*ep^5*x5+148*ep^5*x5^2+
         64*ep^5*x5^3+13*ep^5*x5^4+ep^5*x5^5,12-34*n5+34*n5^2-14*n5^3+2*n5^4-31
         *n2+85*n2*n5-81*n2*n5^2+31*n2*n5^3-4*n2*n5^4+27*n2^2-70*n2^2*n5+61*
         n2^2*n5^2-20*n2^2*n5^3+2*n2^2*n5^4-9*n2^3+21*n2^3*n5-15*n2^3*n5^2+3*
         n2^3*n5^3+n2^4-2*n2^4*n5+n2^4*n5^2-96*ep-34*ep*x5+214*ep*n5+68*ep*n5*
         x5-160*ep*n5^2-42*ep*n5^2*x5+46*ep*n5^3+8*ep*n5^3*x5-4*ep*n5^4+217*ep*
         n2+85*ep*n2*x5-448*ep*n2*n5-162*ep*n2*n5*x5+297*ep*n2*n5^2+93*ep*n2*
         n5^2*x5-70*ep*n2*n5^3-16*ep*n2*n5^3*x5+4*ep*n2*n5^4-162*ep*n2^2-70*ep*
         n2^2*x5+296*ep*n2^2*n5+122*ep*n2^2*n5*x5-158*ep*n2^2*n5^2-60*ep*n2^2*
         n5^2*x5+24*ep*n2^2*n5^3+8*ep*n2^2*n5^3*x5+45*ep*n2^3+21*ep*n2^3*x5-66*
         ep*n2^3*n5-30*ep*n2^3*n5*x5+21*ep*n2^3*n5^2+9*ep*n2^3*n5^2*x5-4*ep*
         n2^4-2*ep*n2^4*x5+4*ep*n2^4*n5+2*ep*n2^4*n5*x5+300*ep^2+214*ep^2*x5+34
         *ep^2*x5^2-492*ep^2*n5-320*ep^2*n5*x5-42*ep^2*n5*x5^2+244*ep^2*n5^2+
         138*ep^2*n5^2*x5+12*ep^2*n5^2*x5^2-36*ep^2*n5^3-16*ep^2*n5^3*x5-582*
         ep^2*n2-448*ep^2*n2*x5-81*ep^2*n2*x5^2+840*ep^2*n2*n5+594*ep^2*n2*n5*
         x5+93*ep^2*n2*n5*x5^2-342*ep^2*n2*n5^2-210*ep^2*n2*n5^2*x5-24*ep^2*n2*
         n5^2*x5^2+36*ep^2*n2*n5^3+16*ep^2*n2*n5^3*x5+350*ep^2*n2^2+296*ep^2*
         n2^2*x5+61*ep^2*n2^2*x5^2-396*ep^2*n2^2*n5-316*ep^2*n2^2*n5*x5-60*ep^2
         *n2^2*n5*x5^2+98*ep^2*n2^2*n5^2+72*ep^2*n2^2*n5^2*x5+12*ep^2*n2^2*n5^2
         *x5^2-72*ep^2*n2^3-66*ep^2*n2^3*x5-15*ep^2*n2^3*x5^2+48*ep^2*n2^3*n5+
         42*ep^2*n2^3*n5*x5+9*ep^2*n2^3*n5*x5^2+4*ep^2*n2^4+4*ep^2*n2^4*x5+ep^2
         *n2^4*x5^2-456*ep^3-492*ep^3*x5-160*ep^3*x5^2-14*ep^3*x5^3+488*ep^3*n5
         +488*ep^3*n5*x5+138*ep^3*n5*x5^2+8*ep^3*n5*x5^3-120*ep^3*n5^2-108*ep^3
         *n5^2*x5-24*ep^3*n5^2*x5^2+740*ep^3*n2+840*ep^3*n2*x5+297*ep^3*n2*x5^2
         +31*ep^3*n2*x5^3-656*ep^3*n2*n5-684*ep^3*n2*n5*x5-210*ep^3*n2*n5*x5^2-
         16*ep^3*n2*n5*x5^3+120*ep^3*n2*n5^2+108*ep^3*n2*n5^2*x5+24*ep^3*n2*
         n5^2*x5^2-320*ep^3*n2^2-396*ep^3*n2^2*x5-158*ep^3*n2^2*x5^2-20*ep^3*
         n2^2*x5^3+168*ep^3*n2^2*n5+196*ep^3*n2^2*n5*x5+72*ep^3*n2^2*n5*x5^2+8*
         ep^3*n2^2*n5*x5^3+36*ep^3*n2^3+48*ep^3*n2^3*x5+21*ep^3*n2^3*x5^2+3*
         ep^3*n2^3*x5^3+336*ep^4+488*ep^4*x5+244*ep^4*x5^2+46*ep^4*x5^3+2*ep^4*
         x5^4-176*ep^4*n5-240*ep^4*n5*x5-108*ep^4*n5*x5^2-16*ep^4*n5*x5^3-440*
         ep^4*n2-656*ep^4*n2*x5-342*ep^4*n2*x5^2-70*ep^4*n2*x5^3-4*ep^4*n2*x5^4
         +176*ep^4*n2*n5+240*ep^4*n2*n5*x5+108*ep^4*n2*n5*x5^2+16*ep^4*n2*n5*
         x5^3+104*ep^4*n2^2+168*ep^4*n2^2*x5+98*ep^4*n2^2*x5^2+24*ep^4*n2^2*
         x5^3+2*ep^4*n2^2*x5^4-96*ep^5-176*ep^5*x5-120*ep^5*x5^2-36*ep^5*x5^3-4
         *ep^5*x5^4+96*ep^5*n2+176*ep^5*n2*x5+120*ep^5*n2*x5^2+36*ep^5*n2*x5^3+
         4*ep^5*n2*x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,2,n2,0,1,n5,0,0,1,0)*(rat(162-297*n5+216*n5^2-78*
         n5^3+14*n5^4-n5^5-375*n2+569*n2*n5-325*n2*n5^2+83*n2*n5^3-8*n2*n5^4+
         323*n2^2-373*n2^2*n5+145*n2^2*n5^2-19*n2^2*n5^3-132*n2^3+102*n2^3*n5-
         20*n2^3*n5^2+26*n2^4-10*n2^4*n5-2*n2^5-816*ep-297*ep*x5+1169*ep*n5+432
         *ep*n5*x5-623*ep*n5^2-234*ep*n5^2*x5+147*ep*n5^3+56*ep*n5^3*x5-13*ep*
         n5^4-5*ep*n5^4*x5+1573*ep*n2+569*ep*n2*x5-1776*ep*n2*n5-650*ep*n2*n5*
         x5+673*ep*n2*n5^2+249*ep*n2*n5^2*x5-86*ep*n2*n5^3-32*ep*n2*n5^3*x5-
         1035*ep*n2^2-373*ep*n2^2*x5+796*ep*n2^2*n5+290*ep*n2^2*n5*x5-155*ep*
         n2^2*n5^2-57*ep*n2^2*n5^2*x5+284*ep*n2^3+102*ep*n2^3*x5-110*ep*n2^3*n5
         -40*ep*n2^3*n5*x5-28*ep*n2^4-10*ep*n2^4*x5+1562*ep^2+1169*ep^2*x5+216*
         ep^2*x5^2-1624*ep^2*n5-1246*ep^2*n5*x5-234*ep^2*n5*x5^2+558*ep^2*n5^2+
         441*ep^2*n5^2*x5+84*ep^2*n5^2*x5^2-64*ep^2*n5^3-52*ep^2*n5^3*x5-10*
         ep^2*n5^3*x5^2-2402*ep^2*n2-1776*ep^2*n2*x5-325*ep^2*n2*x5^2+1792*ep^2
         *n2*n5+1346*ep^2*n2*n5*x5+249*ep^2*n2*n5*x5^2-338*ep^2*n2*n5^2-258*
         ep^2*n2*n5^2*x5-48*ep^2*n2*n5^2*x5^2+1082*ep^2*n2^2+796*ep^2*n2^2*x5+
         145*ep^2*n2^2*x5^2-416*ep^2*n2^2*n5-310*ep^2*n2^2*n5*x5-57*ep^2*n2^2*
         n5*x5^2-150*ep^2*n2^3-110*ep^2*n2^3*x5-20*ep^2*n2^3*x5^2-1388*ep^3-
         1624*ep^3*x5-623*ep^3*x5^2-78*ep^3*x5^3+912*ep^3*n5+1116*ep^3*n5*x5+
         441*ep^3*n5*x5^2+56*ep^3*n5*x5^3-148*ep^3*n5^2-192*ep^3*n5^2*x5-78*
         ep^3*n5^2*x5^2-10*ep^3*n5^2*x5^3+1572*ep^3*n2+1792*ep^3*n2*x5+673*ep^3
         *n2*x5^2+83*ep^3*n2*x5^3-580*ep^3*n2*n5-676*ep^3*n2*n5*x5-258*ep^3*n2*
         n5*x5^2-32*ep^3*n2*n5*x5^3-368*ep^3*n2^2-416*ep^3*n2^2*x5-155*ep^3*
         n2^2*x5^2-19*ep^3*n2^2*x5^3+544*ep^4+912*ep^4*x5+558*ep^4*x5^2+147*
         ep^4*x5^3+14*ep^4*x5^4-160*ep^4*n5-296*ep^4*n5*x5-192*ep^4*n5*x5^2-52*
         ep^4*n5*x5^3-5*ep^4*n5*x5^4-368*ep^4*n2-580*ep^4*n2*x5-338*ep^4*n2*
         x5^2-86*ep^4*n2*x5^3-8*ep^4*n2*x5^4-64*ep^5-160*ep^5*x5-148*ep^5*x5^2-
         64*ep^5*x5^3-13*ep^5*x5^4-ep^5*x5^5,12-34*n5+34*n5^2-14*n5^3+2*n5^4-31
         *n2+85*n2*n5-81*n2*n5^2+31*n2*n5^3-4*n2*n5^4+27*n2^2-70*n2^2*n5+61*
         n2^2*n5^2-20*n2^2*n5^3+2*n2^2*n5^4-9*n2^3+21*n2^3*n5-15*n2^3*n5^2+3*
         n2^3*n5^3+n2^4-2*n2^4*n5+n2^4*n5^2-96*ep-34*ep*x5+214*ep*n5+68*ep*n5*
         x5-160*ep*n5^2-42*ep*n5^2*x5+46*ep*n5^3+8*ep*n5^3*x5-4*ep*n5^4+217*ep*
         n2+85*ep*n2*x5-448*ep*n2*n5-162*ep*n2*n5*x5+297*ep*n2*n5^2+93*ep*n2*
         n5^2*x5-70*ep*n2*n5^3-16*ep*n2*n5^3*x5+4*ep*n2*n5^4-162*ep*n2^2-70*ep*
         n2^2*x5+296*ep*n2^2*n5+122*ep*n2^2*n5*x5-158*ep*n2^2*n5^2-60*ep*n2^2*
         n5^2*x5+24*ep*n2^2*n5^3+8*ep*n2^2*n5^3*x5+45*ep*n2^3+21*ep*n2^3*x5-66*
         ep*n2^3*n5-30*ep*n2^3*n5*x5+21*ep*n2^3*n5^2+9*ep*n2^3*n5^2*x5-4*ep*
         n2^4-2*ep*n2^4*x5+4*ep*n2^4*n5+2*ep*n2^4*n5*x5+300*ep^2+214*ep^2*x5+34
         *ep^2*x5^2-492*ep^2*n5-320*ep^2*n5*x5-42*ep^2*n5*x5^2+244*ep^2*n5^2+
         138*ep^2*n5^2*x5+12*ep^2*n5^2*x5^2-36*ep^2*n5^3-16*ep^2*n5^3*x5-582*
         ep^2*n2-448*ep^2*n2*x5-81*ep^2*n2*x5^2+840*ep^2*n2*n5+594*ep^2*n2*n5*
         x5+93*ep^2*n2*n5*x5^2-342*ep^2*n2*n5^2-210*ep^2*n2*n5^2*x5-24*ep^2*n2*
         n5^2*x5^2+36*ep^2*n2*n5^3+16*ep^2*n2*n5^3*x5+350*ep^2*n2^2+296*ep^2*
         n2^2*x5+61*ep^2*n2^2*x5^2-396*ep^2*n2^2*n5-316*ep^2*n2^2*n5*x5-60*ep^2
         *n2^2*n5*x5^2+98*ep^2*n2^2*n5^2+72*ep^2*n2^2*n5^2*x5+12*ep^2*n2^2*n5^2
         *x5^2-72*ep^2*n2^3-66*ep^2*n2^3*x5-15*ep^2*n2^3*x5^2+48*ep^2*n2^3*n5+
         42*ep^2*n2^3*n5*x5+9*ep^2*n2^3*n5*x5^2+4*ep^2*n2^4+4*ep^2*n2^4*x5+ep^2
         *n2^4*x5^2-456*ep^3-492*ep^3*x5-160*ep^3*x5^2-14*ep^3*x5^3+488*ep^3*n5
         +488*ep^3*n5*x5+138*ep^3*n5*x5^2+8*ep^3*n5*x5^3-120*ep^3*n5^2-108*ep^3
         *n5^2*x5-24*ep^3*n5^2*x5^2+740*ep^3*n2+840*ep^3*n2*x5+297*ep^3*n2*x5^2
         +31*ep^3*n2*x5^3-656*ep^3*n2*n5-684*ep^3*n2*n5*x5-210*ep^3*n2*n5*x5^2-
         16*ep^3*n2*n5*x5^3+120*ep^3*n2*n5^2+108*ep^3*n2*n5^2*x5+24*ep^3*n2*
         n5^2*x5^2-320*ep^3*n2^2-396*ep^3*n2^2*x5-158*ep^3*n2^2*x5^2-20*ep^3*
         n2^2*x5^3+168*ep^3*n2^2*n5+196*ep^3*n2^2*n5*x5+72*ep^3*n2^2*n5*x5^2+8*
         ep^3*n2^2*n5*x5^3+36*ep^3*n2^3+48*ep^3*n2^3*x5+21*ep^3*n2^3*x5^2+3*
         ep^3*n2^3*x5^3+336*ep^4+488*ep^4*x5+244*ep^4*x5^2+46*ep^4*x5^3+2*ep^4*
         x5^4-176*ep^4*n5-240*ep^4*n5*x5-108*ep^4*n5*x5^2-16*ep^4*n5*x5^3-440*
         ep^4*n2-656*ep^4*n2*x5-342*ep^4*n2*x5^2-70*ep^4*n2*x5^3-4*ep^4*n2*x5^4
         +176*ep^4*n2*n5+240*ep^4*n2*n5*x5+108*ep^4*n2*n5*x5^2+16*ep^4*n2*n5*
         x5^3+104*ep^4*n2^2+168*ep^4*n2^2*x5+98*ep^4*n2^2*x5^2+24*ep^4*n2^2*
         x5^3+2*ep^4*n2^2*x5^4-96*ep^5-176*ep^5*x5-120*ep^5*x5^2-36*ep^5*x5^3-4
         *ep^5*x5^4+96*ep^5*n2+176*ep^5*n2*x5+120*ep^5*n2*x5^2+36*ep^5*n2*x5^3+
         4*ep^5*n2*x5^4));
*--#] no,noints=5,signs=1,+,1,1,x,1,1,1,0,order=9,8,7,6,5,4,3,2,1,level=2,stop=2 : 
* LHS = Z(1,1,1,1,n5,1,1,1,0)
* DEN = den(-1+n5)*den(-1+n5+ep)*den(-1+n5+2*ep)*den(-1+n5+3*ep)*den(-1+n5+4*ep)
* MAX = Z(0,1,0,0,1,0,1,1,0)
* MIN = Z(-1,-1,-1,-1,-1,-1,-1,-1,0)
* nterms = 14 (Z: 1, Y: 13)
* NonIntegerIndices: 5
* Signs: 1,1,1,1,+,1,1,1,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,1,1,1,1,n5?{>=2},1,1,1,0) =
*      #include no.rules # no,noints=5,signs=1,1,1,1,+,1,1,1,0,level=1,stop=1
*
*--#[ no,noints=5,signs=1,1,1,1,+,1,1,1,0,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,0,0,0,0,1,1,1,1,-1+n5,1,1,1,0)*(rat(1-n5-4*ep-ep*x5,-1+n5
         +ep*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,1,1,0,n5,1,2,1,0)*(rat(-2+6*n5-6*n5^2+2*n5^3+10
         *ep+6*ep*x5-20*ep*n5-12*ep*n5*x5+10*ep*n5^2+6*ep*n5^2*x5-14*ep^2-20*
         ep^2*x5-6*ep^2*x5^2+14*ep^2*n5+20*ep^2*n5*x5+6*ep^2*n5*x5^2+8*ep^3+14*
         ep^3*x5+10*ep^3*x5^2+2*ep^3*x5^3,1-4*n5+6*n5^2-4*n5^3+n5^4-7*ep-4*ep*
         x5+21*ep*n5+12*ep*n5*x5-21*ep*n5^2-12*ep*n5^2*x5+7*ep*n5^3+4*ep*n5^3*
         x5+14*ep^2+21*ep^2*x5+6*ep^2*x5^2-28*ep^2*n5-42*ep^2*n5*x5-12*ep^2*n5*
         x5^2+14*ep^2*n5^2+21*ep^2*n5^2*x5+6*ep^2*n5^2*x5^2-8*ep^3-28*ep^3*x5-
         21*ep^3*x5^2-4*ep^3*x5^3+8*ep^3*n5+28*ep^3*n5*x5+21*ep^3*n5*x5^2+4*
         ep^3*n5*x5^3+8*ep^4*x5+14*ep^4*x5^2+7*ep^4*x5^3+ep^4*x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,2,0,1,n5,1,1,1,0)*(rat(1-n5-4*ep-ep*x5,1-2*n5+
         n5^2-3*ep-2*ep*x5+3*ep*n5+2*ep*n5*x5+3*ep^2*x5+ep^2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,2,1,1,n5,1,1,0,0)*(rat(4*ep,1-2*n5+n5^2-2*ep-2*
         ep*x5+2*ep*n5+2*ep*n5*x5+2*ep^2*x5+ep^2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,0,1,0,n5,1,2,1,0)*(rat(2-6*n5+6*n5^2-2*n5^3-10*
         ep-6*ep*x5+20*ep*n5+12*ep*n5*x5-10*ep*n5^2-6*ep*n5^2*x5+14*ep^2+20*
         ep^2*x5+6*ep^2*x5^2-14*ep^2*n5-20*ep^2*n5*x5-6*ep^2*n5*x5^2-8*ep^3-14*
         ep^3*x5-10*ep^3*x5^2-2*ep^3*x5^3,1-4*n5+6*n5^2-4*n5^3+n5^4-7*ep-4*ep*
         x5+21*ep*n5+12*ep*n5*x5-21*ep*n5^2-12*ep*n5^2*x5+7*ep*n5^3+4*ep*n5^3*
         x5+14*ep^2+21*ep^2*x5+6*ep^2*x5^2-28*ep^2*n5-42*ep^2*n5*x5-12*ep^2*n5*
         x5^2+14*ep^2*n5^2+21*ep^2*n5^2*x5+6*ep^2*n5^2*x5^2-8*ep^3-28*ep^3*x5-
         21*ep^3*x5^2-4*ep^3*x5^3+8*ep^3*n5+28*ep^3*n5*x5+21*ep^3*n5*x5^2+4*
         ep^3*n5*x5^3+8*ep^4*x5+14*ep^4*x5^2+7*ep^4*x5^3+ep^4*x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,0,1,1,-1+n5,1,1,2,0)*(rat(2-2*n5-4*ep-2*ep*x5,1
         -2*n5+n5^2-ep-2*ep*x5+ep*n5+2*ep*n5*x5+ep^2*x5+ep^2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,0,1,1,n5,0,1,2,0)*(rat(10*ep-20*ep*n5+10*ep*
         n5^2-42*ep^2-20*ep^2*x5+42*ep^2*n5+20*ep^2*n5*x5+40*ep^3+42*ep^3*x5+10
         *ep^3*x5^2,1-4*n5+6*n5^2-4*n5^3+n5^4-7*ep-4*ep*x5+21*ep*n5+12*ep*n5*x5
         -21*ep*n5^2-12*ep*n5^2*x5+7*ep*n5^3+4*ep*n5^3*x5+14*ep^2+21*ep^2*x5+6*
         ep^2*x5^2-28*ep^2*n5-42*ep^2*n5*x5-12*ep^2*n5*x5^2+14*ep^2*n5^2+21*
         ep^2*n5^2*x5+6*ep^2*n5^2*x5^2-8*ep^3-28*ep^3*x5-21*ep^3*x5^2-4*ep^3*
         x5^3+8*ep^3*n5+28*ep^3*n5*x5+21*ep^3*n5*x5^2+4*ep^3*n5*x5^3+8*ep^4*x5+
         14*ep^4*x5^2+7*ep^4*x5^3+ep^4*x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1,0,1,n5,0,1,2,0)*(rat(-2+6*n5-6*n5^2+2*n5^3+6*
         ep+6*ep*x5-12*ep*n5-12*ep*n5*x5+6*ep*n5^2+6*ep*n5^2*x5+2*ep^2-12*ep^2*
         x5-6*ep^2*x5^2-2*ep^2*n5+12*ep^2*n5*x5+6*ep^2*n5*x5^2-8*ep^3-2*ep^3*x5
         +6*ep^3*x5^2+2*ep^3*x5^3,1-4*n5+6*n5^2-4*n5^3+n5^4-7*ep-4*ep*x5+21*ep*
         n5+12*ep*n5*x5-21*ep*n5^2-12*ep*n5^2*x5+7*ep*n5^3+4*ep*n5^3*x5+14*ep^2
         +21*ep^2*x5+6*ep^2*x5^2-28*ep^2*n5-42*ep^2*n5*x5-12*ep^2*n5*x5^2+14*
         ep^2*n5^2+21*ep^2*n5^2*x5+6*ep^2*n5^2*x5^2-8*ep^3-28*ep^3*x5-21*ep^3*
         x5^2-4*ep^3*x5^3+8*ep^3*n5+28*ep^3*n5*x5+21*ep^3*n5*x5^2+4*ep^3*n5*
         x5^3+8*ep^4*x5+14*ep^4*x5^2+7*ep^4*x5^3+ep^4*x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1,1,0,1+n5,0,1,1,0)*(rat(-2*n5+6*n5^2-6*n5^3+2*
         n5^4-2*ep*x5+12*ep*n5+12*ep*n5*x5-24*ep*n5^2-18*ep*n5^2*x5+12*ep*n5^3+
         8*ep*n5^3*x5+12*ep^2*x5+6*ep^2*x5^2-28*ep^2*n5-48*ep^2*n5*x5-18*ep^2*
         n5*x5^2+28*ep^2*n5^2+36*ep^2*n5^2*x5+12*ep^2*n5^2*x5^2-28*ep^3*x5-24*
         ep^3*x5^2-6*ep^3*x5^3+32*ep^3*n5+56*ep^3*n5*x5+36*ep^3*n5*x5^2+8*ep^3*
         n5*x5^3+32*ep^4*x5+28*ep^4*x5^2+12*ep^4*x5^3+2*ep^4*x5^4,1-4*n5+6*n5^2
         -4*n5^3+n5^4-9*ep-4*ep*x5+27*ep*n5+12*ep*n5*x5-27*ep*n5^2-12*ep*n5^2*
         x5+9*ep*n5^3+4*ep*n5^3*x5+26*ep^2+27*ep^2*x5+6*ep^2*x5^2-52*ep^2*n5-54
         *ep^2*n5*x5-12*ep^2*n5*x5^2+26*ep^2*n5^2+27*ep^2*n5^2*x5+6*ep^2*n5^2*
         x5^2-24*ep^3-52*ep^3*x5-27*ep^3*x5^2-4*ep^3*x5^3+24*ep^3*n5+52*ep^3*n5
         *x5+27*ep^3*n5*x5^2+4*ep^3*n5*x5^3+24*ep^4*x5+26*ep^4*x5^2+9*ep^4*x5^3
         +ep^4*x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1,1,0,1+n5,1,1,0,0)*(rat(2*n5-6*n5^2+6*n5^3-2*
         n5^4+2*ep*x5-8*ep*n5-12*ep*n5*x5+16*ep*n5^2+18*ep*n5^2*x5-8*ep*n5^3-8*
         ep*n5^3*x5-8*ep^2*x5-6*ep^2*x5^2+2*ep^2*n5+32*ep^2*n5*x5+18*ep^2*n5*
         x5^2-2*ep^2*n5^2-24*ep^2*n5^2*x5-12*ep^2*n5^2*x5^2+2*ep^3*x5+16*ep^3*
         x5^2+6*ep^3*x5^3+8*ep^3*n5-4*ep^3*n5*x5-24*ep^3*n5*x5^2-8*ep^3*n5*x5^3
         +8*ep^4*x5-2*ep^4*x5^2-8*ep^4*x5^3-2*ep^4*x5^4,1-4*n5+6*n5^2-4*n5^3+
         n5^4-7*ep-4*ep*x5+21*ep*n5+12*ep*n5*x5-21*ep*n5^2-12*ep*n5^2*x5+7*ep*
         n5^3+4*ep*n5^3*x5+14*ep^2+21*ep^2*x5+6*ep^2*x5^2-28*ep^2*n5-42*ep^2*n5
         *x5-12*ep^2*n5*x5^2+14*ep^2*n5^2+21*ep^2*n5^2*x5+6*ep^2*n5^2*x5^2-8*
         ep^3-28*ep^3*x5-21*ep^3*x5^2-4*ep^3*x5^3+8*ep^3*n5+28*ep^3*n5*x5+21*
         ep^3*n5*x5^2+4*ep^3*n5*x5^3+8*ep^4*x5+14*ep^4*x5^2+7*ep^4*x5^3+ep^4*
         x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1,1,1,1+n5,0,0,1,0)*(rat(6*n5-18*n5^2+18*n5^3-6
         *n5^4+6*ep*x5-36*ep*n5-36*ep*n5*x5+72*ep*n5^2+54*ep*n5^2*x5-36*ep*n5^3
         -24*ep*n5^3*x5-36*ep^2*x5-18*ep^2*x5^2+66*ep^2*n5+144*ep^2*n5*x5+54*
         ep^2*n5*x5^2-66*ep^2*n5^2-108*ep^2*n5^2*x5-36*ep^2*n5^2*x5^2+66*ep^3*
         x5+72*ep^3*x5^2+18*ep^3*x5^3-40*ep^3*n5-132*ep^3*n5*x5-108*ep^3*n5*
         x5^2-24*ep^3*n5*x5^3-40*ep^4*x5-66*ep^4*x5^2-36*ep^4*x5^3-6*ep^4*x5^4,
         1-4*n5+6*n5^2-4*n5^3+n5^4-7*ep-4*ep*x5+21*ep*n5+12*ep*n5*x5-21*ep*n5^2
         -12*ep*n5^2*x5+7*ep*n5^3+4*ep*n5^3*x5+14*ep^2+21*ep^2*x5+6*ep^2*x5^2-
         28*ep^2*n5-42*ep^2*n5*x5-12*ep^2*n5*x5^2+14*ep^2*n5^2+21*ep^2*n5^2*x5+
         6*ep^2*n5^2*x5^2-8*ep^3-28*ep^3*x5-21*ep^3*x5^2-4*ep^3*x5^3+8*ep^3*n5+
         28*ep^3*n5*x5+21*ep^3*n5*x5^2+4*ep^3*n5*x5^3+8*ep^4*x5+14*ep^4*x5^2+7*
         ep^4*x5^3+ep^4*x5^4))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1,1,1,1+n5,1,0,0,0)*(rat(6*n5-12*n5^2+6*n5^3+6*
         ep*x5-24*ep*n5-24*ep*n5*x5+24*ep*n5^2+18*ep*n5^2*x5-24*ep^2*x5-12*ep^2
         *x5^2+16*ep^2*n5+48*ep^2*n5*x5+18*ep^2*n5*x5^2+16*ep^3*x5+24*ep^3*x5^2
         +6*ep^3*x5^3,-1+3*n5-3*n5^2+n5^3+5*ep+3*ep*x5-10*ep*n5-6*ep*n5*x5+5*ep
         *n5^2+3*ep*n5^2*x5-6*ep^2-10*ep^2*x5-3*ep^2*x5^2+6*ep^2*n5+10*ep^2*n5*
         x5+3*ep^2*n5*x5^2+6*ep^3*x5+5*ep^3*x5^2+ep^3*x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,2,0,1,n5,1,0,1,0)*(rat(4*ep,1-2*n5+n5^2-2*ep-2*
         ep*x5+2*ep*n5+2*ep*n5*x5+2*ep^2*x5+ep^2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,2,1,1,n5,1,0,0,0)*(rat(1-2*n5+n5^2+2*ep-2*ep*x5
         -2*ep*n5+2*ep*n5*x5-16*ep^2-2*ep^2*x5+ep^2*x5^2,-1+3*n5-3*n5^2+n5^3+5*
         ep+3*ep*x5-10*ep*n5-6*ep*n5*x5+5*ep*n5^2+3*ep*n5^2*x5-6*ep^2-10*ep^2*
         x5-3*ep^2*x5^2+6*ep^2*n5+10*ep^2*n5*x5+3*ep^2*n5*x5^2+6*ep^3*x5+5*ep^3
         *x5^2+ep^3*x5^3));
*--#] no,noints=5,signs=1,1,1,1,+,1,1,1,0,level=1,stop=1 : 
* LHS = Z(1,1,1,1,n5,1,1,1,0)
* DEN = den(ep)*den(n5+ep)*den(n5+2*ep)*den(n5+4*ep)
* MAX = Z(1,1,1,0,1,0,1,1,0)
* MIN = Z(-1,-1,-1,-1,1,-1,-1,-1,0)
* nterms = 15 (Z: 1, Y: 14)
* NonIntegerIndices: 5
* Signs: 1,1,1,1,-,1,1,1,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 2
* Stop: 2
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},0,0,0,0,1,1,1,1,n5?neg0_,1,1,1,0) =
*      #include no.rules # no,noints=5,signs=1,1,1,1,-,1,1,1,0,level=2,stop=2
*
*--#[ no,noints=5,signs=1,1,1,1,-,1,1,1,0,level=2,stop=2 :
      +`Z'(0,0,0,0,x5,0,0,0,0,1,1,1,1,1+n5,1,1,1,0)*(rat(-n5-ep*x5,n5+4*ep+ep*
         x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,1,1,0,1+n5,1,2,1,0)*(rat(-n5^3-3*ep*n5^2-3*ep*
         n5^2*x5+3*ep^2*n5-6*ep^2*n5*x5-3*ep^2*n5*x5^2+4*ep^3+3*ep^3*x5-3*ep^3*
         x5^2-ep^3*x5^3,ep*n5^3+7*ep^2*n5^2+3*ep^2*n5^2*x5+14*ep^3*n5+14*ep^3*
         n5*x5+3*ep^3*n5*x5^2+8*ep^4+14*ep^4*x5+7*ep^4*x5^2+ep^4*x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,1,1,1,1+n5,1,1,1,0)*(rat(-6*n5^2-24*ep*n5-12*ep
         *n5*x5-16*ep^2-24*ep^2*x5-6*ep^2*x5^2,n5^2+5*ep*n5+2*ep*n5*x5+4*ep^2+5
         *ep^2*x5+ep^2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,2,0,1,1+n5,1,1,1,0)*(rat(5*n5+4*ep+5*ep*x5,n5^2
         +5*ep*n5+2*ep*n5*x5+4*ep^2+5*ep^2*x5+ep^2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,0,2,1,1,1+n5,1,1,0,0)*(rat(-6*n5^2-20*ep*n5-12*ep
         *n5*x5-12*ep^2-20*ep^2*x5-6*ep^2*x5^2,n5^3+7*ep*n5^2+3*ep*n5^2*x5+14*
         ep^2*n5+14*ep^2*n5*x5+3*ep^2*n5*x5^2+8*ep^3+14*ep^3*x5+7*ep^3*x5^2+
         ep^3*x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,0,1,0,1+n5,1,2,1,0)*(rat(n5^3+3*ep*n5^2+3*ep*
         n5^2*x5-3*ep^2*n5+6*ep^2*n5*x5+3*ep^2*n5*x5^2-4*ep^3-3*ep^3*x5+3*ep^3*
         x5^2+ep^3*x5^3,ep*n5^3+7*ep^2*n5^2+3*ep^2*n5^2*x5+14*ep^3*n5+14*ep^3*
         n5*x5+3*ep^3*n5*x5^2+8*ep^4+14*ep^4*x5+7*ep^4*x5^2+ep^4*x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,0,1,1,1+n5,0,1,2,0)*(rat(-2*n5^2-7*ep*n5-4*ep*
         n5*x5-4*ep^2-7*ep^2*x5-2*ep^2*x5^2,ep*n5^2+6*ep^2*n5+2*ep^2*n5*x5+8*
         ep^3+6*ep^3*x5+ep^3*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,0,1,1,1+n5,1,1,1,0)*(rat(2*n5+4*ep+2*ep*x5,n5+4
         *ep+ep*x5))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,0,2,0,1+n5,1,1,1,0)*(rat(n5^3+7*ep*n5^2+3*ep*
         n5^2*x5+15*ep^2*n5+14*ep^2*n5*x5+3*ep^2*n5*x5^2+8*ep^3+15*ep^3*x5+7*
         ep^3*x5^2+ep^3*x5^3,ep*n5^3+7*ep^2*n5^2+3*ep^2*n5^2*x5+14*ep^3*n5+14*
         ep^3*n5*x5+3*ep^3*n5*x5^2+8*ep^4+14*ep^4*x5+7*ep^4*x5^2+ep^4*x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1,0,1,1+n5,0,1,2,0)*(rat(2*n5^2+7*ep*n5+4*ep*n5
         *x5+4*ep^2+7*ep^2*x5+2*ep^2*x5^2,ep*n5^2+6*ep^2*n5+2*ep^2*n5*x5+8*ep^3
         +6*ep^3*x5+ep^3*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,1,2,0,1+n5,1,1,0,0)*(rat(-n5^3-7*ep*n5^2-3*ep*
         n5^2*x5-15*ep^2*n5-14*ep^2*n5*x5-3*ep^2*n5*x5^2-8*ep^3-15*ep^3*x5-7*
         ep^3*x5^2-ep^3*x5^3,ep*n5^3+7*ep^2*n5^2+3*ep^2*n5^2*x5+14*ep^3*n5+14*
         ep^3*n5*x5+3*ep^3*n5*x5^2+8*ep^4+14*ep^4*x5+7*ep^4*x5^2+ep^4*x5^3))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,2,0,1,1+n5,1,0,1,0)*(rat(4*ep,n5^2+6*ep*n5+2*ep
         *n5*x5+8*ep^2+6*ep^2*x5+ep^2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,1,2,1,1,1+n5,1,0,0,0)*(rat(n5+ep*x5,n5^2+5*ep*n5+
         2*ep*n5*x5+4*ep^2+5*ep^2*x5+ep^2*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,2,0,1,1,1+n5,0,1,1,0)*(rat(-2*n5^2-9*ep*n5-4*ep*
         n5*x5-8*ep^2-9*ep^2*x5-2*ep^2*x5^2,ep*n5^2+6*ep^2*n5+2*ep^2*n5*x5+8*
         ep^3+6*ep^3*x5+ep^3*x5^2))
      +`Y'(0,0,0,0,x5,0,0,0,0,2,1,1,1,1+n5,0,0,1,0)*(rat(2*n5^2+9*ep*n5+4*ep*n5
         *x5+8*ep^2+9*ep^2*x5+2*ep^2*x5^2,ep*n5^2+6*ep^2*n5+2*ep^2*n5*x5+8*ep^3
         +6*ep^3*x5+ep^3*x5^2));
*--#] no,noints=5,signs=1,1,1,1,-,1,1,1,0,level=2,stop=2 : 
*--#] no,5-rules
        goto 2;
label 1;
        #call ReduceZ2Y(Z,Y,9,ifmatch->2,1,2,3,4,6,7,8)
        $doloop = 1;
label 2;
        ModuleOption maximum, $doloop;
        .sort:`TOPO'-`nloop++';
        #if `$doloop'
          #redefine loop "0"
        #endif
      #enddo
      #break
*--#] no,5
*--#[ no,6
    #case no,6
      #define nloop "1"
      #do loop=1,1
        #$doloop = 0;
        if (count(`Z',1) == 0) goto 2;
*--#[ no,6-rules
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* DEN = den(-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
* MAX = Z(1,0,1,0,0,0,0,0,2)
* MIN = Z(0,-1,0,0,0,0,-1,-1,1)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: 6
* Signs: x,x,x,o,o,o,x,x,-*?neg_
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 2
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,n1?,n2?,n3?,n4?,n5?,n6?,n7?,n8?,n9?neg_) =
*      #include no.rules # no,noints=6,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0
*
*--#[ no,noints=6,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,1+n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(n1,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x6))
      +`Z'(0,0,0,0,0,x6,0,0,0,1+n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*(rat(-n1,-14+2
         *n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x6))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,-1+n2,1+n3,n4,n5,n6,n7,n8,1+n9)*(rat(n3,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x6))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,-1+n2,n3,n4,n5,n6,n7,n8,2+n9)*(rat(1+n9,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x6))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,1+n3,n4,n5,n6,n7,-1+n8,1+n9)*(rat(-n3,-14+2
         *n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x6))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(-3+n9+n3+2*n2+
         n1+2*ep,-14+2*n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x6)
         );
*--#] no,noints=6,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0 : 
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* DEN = den(-1+n4)
* MAX = Z(0,0,0,0,1,0,1,0,0)
* MIN = Z(-1,-1,-1,-1,0,-1,0,-1,0)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: 6
* Signs: +,x,+,+,x,x,x,x,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,n1?pos_,n2?,n3?pos_,n4?{>=2},n5?,n6?,n7?,n8?,0) =
*      #include no.rules # no,noints=6,signs=+,x,+,+,x,x,x,x,0,level=1,stop=0
*
*--#[ no,noints=6,signs=+,x,+,+,x,x,x,x,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,-1+n1,n2,n3,-1+n4,n5,n6,1+n7,n8,0)*(rat(-n7,-1+n4
         ))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,-1+n2,n3,-1+n4,n5,n6,1+n7,n8,0)*(rat(n7,-1+n4)
         )
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,-1+n3,n4,n5,n6,n7,n8,0)*(rat(1,1))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,n3,-1+n4,1+n5,-1+n6,n7,n8,0)*(rat(-n5,-1+n4
         ))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,n3,-1+n4,1+n5,n6,n7,-1+n8,0)*(rat(n5,-1+n4)
         )
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,n3,-1+n4,n5,n6,n7,n8,0)*(rat(-9+2*n8+n7+n5+
         n4+2*n3+2*n2+4*ep,-1+n4));
*--#] no,noints=6,signs=+,x,+,+,x,x,x,x,0,level=1,stop=0 : 
* LHS = Z(n1,n2,n3,1,n5,n6,n7,n8,0)
* DEN = den(-1+n3)
* MAX = Z(0,1,0,0,0,0,0,1,0)
* MIN = Z(-1,0,-1,-1,-1,-1,-1,0,0)
* nterms = 6 (Z: 5, Y: 1)
* NonIntegerIndices: 6
* Signs: +,x,+,1,x,x,x,x,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,n1?pos_,n2?,n3?{>=2},1,n5?,n6?,n7?,n8?,0) =
*      #include no.rules # no,noints=6,signs=+,x,+,1,x,x,x,x,0,level=1,stop=0
*
*--#[ no,noints=6,signs=+,x,+,1,x,x,x,x,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,-1+n1,1+n2,-1+n3,1,n5,n6,n7,n8,0)*(rat(-n2,-1+n3)
         )
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,1+n2,-1+n3,1,n5,n6,-1+n7,n8,0)*(rat(n2,-1+n3))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,-1+n3,1,-1+n5,n6,n7,1+n8,0)*(rat(n8,-1+n3))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,-1+n3,1,n5,-1+n6,n7,1+n8,0)*(rat(-n8,-1+n3)
         )
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,-1+n3,1,n5,n6,n7,n8,0)*(rat(-7+n8+2*n7+2*n5
         +n3+n2+4*ep,-1+n3))
      +`Y'(0,0,0,0,0,x6,0,0,0,n1,n2,n3,0,n5,n6,n7,n8,0)*(rat(1,1));
*--#] no,noints=6,signs=+,x,+,1,x,x,x,x,0,level=1,stop=0 : 
* LHS = Z(n1,n2,1,1,n5,n6,n7,n8,0)
* DEN = den(-1+n1)*den(-4+n7+n5+n2+n1+2*ep)
* MAX = Z(0,0,0,0,0,1,1,0,0)
* MIN = Z(-1,-1,-1,-1,-1,0,-1,-1,0)
* nterms = 11 (Z: 7, Y: 4)
* NonIntegerIndices: 6
* Signs: +,+,1,1,+,o,+,+,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,n1?{>=2},n2?pos_,1,1,n5?pos_,n6?,n7?pos_,n8?pos_,0) =
*      #include no.rules # no,noints=6,signs=+,+,1,1,+,o,+,+,0,level=1,stop=0
*
*--#[ no,noints=6,signs=+,+,1,1,+,o,+,+,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,-1+n1,n2,1,1,-1+n5,1+n6,n7,n8,0)*(rat(-8*n6+n6*n8
         +2*n6*n7+n6^2+2*n5*n6+n2*n6+n1*n6-8*ep*x6+ep*n8*x6+2*ep*n7*x6+4*ep*n6+
         2*ep*n6*x6+2*ep*n5*x6+ep*n2*x6+ep*n1*x6+4*ep^2*x6+ep^2*x6^2,4-n7-n5-n2
         -5*n1+n1*n7+n1*n5+n1*n2+n1^2-2*ep+2*ep*n1))
      +`Z'(0,0,0,0,0,x6,0,0,0,-1+n1,n2,1,1,-1+n5,n6,1+n7,n8,0)*(rat(-8*n7+n7*n8
         +2*n7^2+n6*n7+2*n5*n7+n2*n7+n1*n7+4*ep*n7+ep*n7*x6,4-n7-n5-n2-5*n1+n1*
         n7+n1*n5+n1*n2+n1^2-2*ep+2*ep*n1))
      +`Z'(0,0,0,0,0,x6,0,0,0,-1+n1,n2,1,1,n5,1+n6,n7,-1+n8,0)*(rat(8*n6-n6*n8-
         2*n6*n7-n6^2-2*n5*n6-n2*n6-n1*n6+8*ep*x6-ep*n8*x6-2*ep*n7*x6-4*ep*n6-2
         *ep*n6*x6-2*ep*n5*x6-ep*n2*x6-ep*n1*x6-4*ep^2*x6-ep^2*x6^2,4-n7-n5-n2-
         5*n1+n1*n7+n1*n5+n1*n2+n1^2-2*ep+2*ep*n1))
      +`Z'(0,0,0,0,0,x6,0,0,0,-1+n1,n2,1,1,n5,n6,n7,n8,0)*(rat(40-5*n8-18*n7+n7
         *n8+2*n7^2-13*n6+n6*n8+3*n6*n7+n6^2-26*n5+2*n5*n8+6*n5*n7+4*n5*n6+4*
         n5^2-5*n2+n2*n7+n2*n6+2*n2*n5-13*n1+n1*n8+3*n1*n7+2*n1*n6+4*n1*n5+n1*
         n2+n1^2-36*ep-13*ep*x6+2*ep*n8+ep*n8*x6+8*ep*n7+3*ep*n7*x6+6*ep*n6+2*
         ep*n6*x6+12*ep*n5+4*ep*n5*x6+2*ep*n2+ep*n2*x6+6*ep*n1+2*ep*n1*x6+8*
         ep^2+6*ep^2*x6+ep^2*x6^2,4-n7-n5-n2-5*n1+n1*n7+n1*n5+n1*n2+n1^2-2*ep+2
         *ep*n1))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,-1+n2,1,1,n5,n6,n7,n8,0)*(rat(8-n8-2*n7-n6-2*
         n5-n2-n1-4*ep-ep*x6,-4+n7+n5+n2+n1+2*ep))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,1,1,-1+n5,n6,n7,n8,0)*(rat(-4+n8+n7+n5+n2+2
         *ep,-4+n7+n5+n2+n1+2*ep))
      +`Z'(0,0,0,0,0,x6,0,0,0,n1,n2,1,1,n5,n6,-1+n7,n8,0)*(rat(-4+n7+n6+n5+n1+2
         *ep+ep*x6,-4+n7+n5+n2+n1+2*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,-1+n1,n2,1,0,n5,n6,1+n7,n8,0)*(rat(10*n7-n7*n8-2*
         n7^2-n6*n7-2*n5*n7-n2*n7-3*n1*n7-4*ep*n7-ep*n7*x6,4-n7-n5-n2-5*n1+n1*
         n7+n1*n5+n1*n2+n1^2-2*ep+2*ep*n1))
      +`Y'(0,0,0,0,0,x6,0,0,0,n1,-1+n2,1,0,n5,n6,1+n7,n8,0)*(rat(2*n7,-4+n7+n5+
         n2+n1+2*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,n1,n2,0,1,n5,n6,n7,n8,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,x6,0,0,0,n1,n2,1,0,n5,n6,n7,n8,0)*(rat(-2+n8+2*n2-n1+2*ep,
         -4+n7+n5+n2+n1+2*ep));
*--#] no,noints=6,signs=+,+,1,1,+,o,+,+,0,level=1,stop=0 : 
* LHS = Z(1,n2,1,1,n5,n6,n7,n8,0)
* DEN = den(-1+n5)*den(-3+n7+n5+n2+2*ep)
* MAX = Z(0,0,0,0,0,1,1,1,0)
* MIN = Z(-1,-1,-1,-1,-1,0,-1,-1,0)
* nterms = 14 (Z: 5, Y: 9)
* NonIntegerIndices: 6
* Signs: 1,+,1,1,+,o,+,+,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,1,n2?pos_,1,1,n5?{>=2},n6?,n7?pos_,n8?pos_,0) =
*      #include no.rules # no,noints=6,signs=1,+,1,1,+,o,+,+,0,level=1,stop=0
*
*--#[ no,noints=6,signs=1,+,1,1,+,o,+,+,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,1,-1+n5,n6,1+n7,n8,0)*(rat(7*n7-n7*n8-2
         *n7^2-n6*n7-2*n5*n7-n2*n7-4*ep*n7-ep*n7*x6,3-n7-4*n5+n5*n7+n5^2-n2+n2*
         n5-2*ep+2*ep*n5))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,1,-1+n5,n6,n7,1+n8,0)*(rat(7*n8-n8^2-2*
         n7*n8-n6*n8-2*n5*n8-n2*n8-4*ep*n8-ep*n8*x6,3-n7-4*n5+n5*n7+n5^2-n2+n2*
         n5-2*ep+2*ep*n5))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,1,n5,n6,n7,n8,0)*(rat(7-n8-2*n7-n6-2*n5
         -n2-4*ep-ep*x6,-3+n7+n5+n2+2*ep))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,n2,1,1,-1+n5,n6,n7,n8,0)*(rat(-31+11*n8-n8^2+16
         *n7-3*n7*n8-2*n7^2+5*n6-n6*n8-n6*n7+12*n5-2*n5*n8-3*n5*n7-n5*n6-n5^2+
         18*n2-3*n2*n8-5*n2*n7-2*n2*n6-4*n2*n5-2*n2^2+32*ep+5*ep*x6-6*ep*n8-ep*
         n8*x6-8*ep*n7-ep*n7*x6-2*ep*n6-6*ep*n5-ep*n5*x6-10*ep*n2-2*ep*n2*x6-8*
         ep^2-2*ep^2*x6,3-n7-4*n5+n5*n7+n5^2-n2+n2*n5-2*ep+2*ep*n5))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,n2,1,1,n5,n6,-1+n7,n8,0)*(rat(4-n8-n7-n5-n2-2*
         ep,-3+n7+n5+n2+2*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,0,n5,n6,1+n7,n8,0)*(rat(-3*n7,-6+2*n7+2*n5
         +2*n2+4*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,-1+n5,1+n6,n7,n8,0)*(rat(-n6-ep*x6,-6+2*
         n7+2*n5+2*n2+4*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,-1+n5,n6,1+n7,n8,0)*(rat(-13*n7+2*n7*n8+
         4*n7^2+2*n6*n7+3*n5*n7+2*n2*n7+8*ep*n7+2*ep*n7*x6,6-2*n7-8*n5+2*n5*n7+
         2*n5^2-2*n2+2*n2*n5-4*ep+4*ep*n5))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,n5,1+n6,n7,-1+n8,0)*(rat(n6+ep*x6,-6+2*
         n7+2*n5+2*n2+4*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,n5,n6,n7,n8,0)*(rat(-10+2*n8+3*n7+n6+2*
         n5+2*n2+6*ep+ep*x6,-6+2*n7+2*n5+2*n2+4*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,0,n5,n6,1+n7,n8,0)*(rat(2*n7,-3+n7+n5+
         n2+2*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,0,1,-1+n5,n6,n7,1+n8,0)*(rat(-7*n8+n8^2+2*n7
         *n8+n6*n8+2*n5*n8+n2*n8+4*ep*n8+ep*n8*x6,3-n7-4*n5+n5*n7+n5^2-n2+n2*n5
         -2*ep+2*ep*n5))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,0,1,n5,n6,n7,n8,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,1,0,n5,n6,n7,n8,0)*(rat(-10+2*n8+2*n7+n6+2*
         n5+3*n2+6*ep+ep*x6,-3+n7+n5+n2+2*ep));
*--#] no,noints=6,signs=1,+,1,1,+,o,+,+,0,level=1,stop=0 : 
* LHS = Z(1,n2,1,1,1,n6,n7,n8,0)
* DEN = den(-1+n8)*den(-1+n7+ep)*den(-2+n7+n2+2*ep)
* MAX = Z(0,1,0,0,1,1,1,1,0)
* MIN = Z(-1,-1,-1,-1,-1,0,-1,-1,0)
* nterms = 21 (Z: 5, Y: 16)
* NonIntegerIndices: 6
* Signs: 1,+,1,1,1,o,+,+,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,1,n2?pos_,1,1,1,n6?,n7?pos_,n8?{>=2},0) =
*      #include no.rules # no,noints=6,signs=1,+,1,1,1,o,+,+,0,level=1,stop=0
*
*--#[ no,noints=6,signs=1,+,1,1,1,o,+,+,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,1,1,n6,n7,n8,0)*(rat(10-7*n8+n8^2-4*n7+
         2*n7*n8-7*n6+2*n6*n8+2*n6*n7+n6^2-2*n2+n2*n8+n2*n6-13*ep-7*ep*x6+5*ep*
         n8+2*ep*n8*x6+2*ep*n7+2*ep*n7*x6+5*ep*n6+2*ep*n6*x6+ep*n2+ep*n2*x6+4*
         ep^2+5*ep^2*x6+ep^2*x6^2,2-3*n7+n7^2-n2+n2*n7-4*ep+3*ep*n7+ep*n2+2*
         ep^2))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,1+n2,1,1,1,n6,-1+n7,-1+n8,0)*(rat(-10*n2+7*n2*
         n8-n2*n8^2+4*n2*n7-2*n2*n7*n8+2*n2*n6-n2*n6*n8+7*n2^2-2*n2^2*n8-2*n2^2
         *n7-n2^2*n6-n2^3+13*ep*n2+2*ep*n2*x6-5*ep*n2*n8-ep*n2*n8*x6-2*ep*n2*n7
         -ep*n2*n6-5*ep*n2^2-ep*n2^2*x6-4*ep^2*n2-ep^2*n2*x6,-2+2*n8+3*n7-3*n7*
         n8-n7^2+n7^2*n8+n2-n2*n8-n2*n7+n2*n7*n8+4*ep-4*ep*n8-3*ep*n7+3*ep*n7*
         n8-ep*n2+ep*n2*n8-2*ep^2+2*ep^2*n8))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,n2,1,1,1,n6,-1+n7,n8,0)*(rat(-18+9*n8-n8^2+10*
         n7-3*n7*n8-n7^2+8*n6-2*n6*n8-2*n6*n7-n6^2+9*n2-2*n2*n8-3*n2*n7-2*n2*n6
         -n2^2+27*ep+8*ep*x6-7*ep*n8-2*ep*n8*x6-7*ep*n7-2*ep*n7*x6-6*ep*n6-2*ep
         *n6*x6-7*ep*n2-2*ep*n2*x6-10*ep^2-6*ep^2*x6-ep^2*x6^2,2-3*n7+n7^2-n2+
         n2*n7-4*ep+3*ep*n7+ep*n2+2*ep^2))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,n2,1,1,1,n6,n7,-1+n8,0)*(rat(15-3*n8-16*n7+2*n7
         *n8+4*n7^2-3*n6+2*n6*n7-23*n2+4*n2*n8+20*n2*n7-2*n2*n7*n8-4*n2*n7^2+4*
         n2*n6-2*n2*n6*n7+9*n2^2-n2^2*n8-4*n2^2*n7-n2^2*n6-n2^3-37*ep-3*ep*x6+5
         *ep*n8+28*ep*n7+2*ep*n7*x6-2*ep*n7*n8-4*ep*n7^2+5*ep*n6-2*ep*n6*n7+36*
         ep*n2+4*ep*n2*x6-3*ep*n2*n8-16*ep*n2*n7-2*ep*n2*n7*x6-3*ep*n2*n6-7*ep*
         n2^2-ep*n2^2*x6+30*ep^2+5*ep^2*x6-2*ep^2*n8-12*ep^2*n7-2*ep^2*n7*x6-2*
         ep^2*n6-14*ep^2*n2-3*ep^2*n2*x6-8*ep^3-2*ep^3*x6,-2+2*n8+3*n7-3*n7*n8-
         n7^2+n7^2*n8+n2-n2*n8-n2*n7+n2*n7*n8+4*ep-4*ep*n8-3*ep*n7+3*ep*n7*n8-
         ep*n2+ep*n2*n8-2*ep^2+2*ep^2*n8))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,n2,1,1,2,n6,-1+n7,-1+n8,0)*(rat(-5+n8+2*n7+n6+6
         *n2-n2*n8-2*n2*n7-n2*n6-n2^2+9*ep+ep*x6-ep*n8-2*ep*n7-ep*n6-5*ep*n2-ep
         *n2*x6-4*ep^2-ep^2*x6,-2+2*n8+3*n7-3*n7*n8-n7^2+n7^2*n8+n2-n2*n8-n2*n7
         +n2*n7*n8+4*ep-4*ep*n8-3*ep*n7+3*ep*n7*n8-ep*n2+ep*n2*n8-2*ep^2+2*ep^2
         *n8))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,1+n2,0,1,1,n6,n7,n8,0)*(rat(5*n2-n2*n8-2*n2*n7-
         n2*n6-n2^2-4*ep*n2-ep*n2*x6,4-6*n7+2*n7^2-2*n2+2*n2*n7-8*ep+6*ep*n7+2*
         ep*n2+4*ep^2))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,1+n2,1,1,1,n6,n7,-1+n8,0)*(rat(15*n2-8*n2*n8+n2
         *n8^2-6*n2*n7+2*n2*n7*n8-3*n2*n6+n2*n6*n8-13*n2^2+3*n2^2*n8+4*n2^2*n7+
         2*n2^2*n6+2*n2^3-22*ep*n2-3*ep*n2*x6+6*ep*n2*n8+ep*n2*n8*x6+4*ep*n2*n7
         +2*ep*n2*n6+10*ep*n2^2+2*ep*n2^2*x6+8*ep^2*n2+2*ep^2*n2*x6,-4+4*n8+6*
         n7-6*n7*n8-2*n7^2+2*n7^2*n8+2*n2-2*n2*n8-2*n2*n7+2*n2*n7*n8+8*ep-8*ep*
         n8-6*ep*n7+6*ep*n7*n8-2*ep*n2+2*ep*n2*n8-4*ep^2+4*ep^2*n8))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,0,1,n6,1+n7,n8,0)*(rat(-2*n7+n7*n8-n7^2+n6
         *n7+n2*n7+ep*n7+ep*n7*x6,4-6*n7+2*n7^2-2*n2+2*n2*n7-8*ep+6*ep*n7+2*ep*
         n2+4*ep^2))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,0,1+n6,n7,n8,0)*(rat(-n6-ep*x6,-4+2*n7+2
         *n2+4*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,0,n6,1+n7,n8,0)*(rat(-n7,-4+2*n7+2*n2+4*
         ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,1,1+n6,n7,-1+n8,0)*(rat(n6+ep*x6,-4+2*n7
         +2*n2+4*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,1,n6,n7,n8,0)*(rat(8-2*n8-6*n7+n7*n8+
         n7^2-n6-7*n2+n2*n8+3*n2*n7+n2*n6+n2^2-14*ep-ep*x6+2*ep*n8+5*ep*n7+ep*
         n6+6*ep*n2+ep*n2*x6+6*ep^2+ep^2*x6,4-6*n7+2*n7^2-2*n2+2*n2*n7-8*ep+6*
         ep*n7+2*ep*n2+4*ep^2))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,0,1,n6,1+n7,n8,0)*(rat(3*n7-n7*n8-n6*n7
         -n2*n7-2*ep*n7-ep*n7*x6,2-3*n7+n7^2-n2+n2*n7-4*ep+3*ep*n7+ep*n2+2*ep^2
         ))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,1,0,n6,1+n7,n8,0)*(rat(-5*n7+n7*n8+2*
         n7^2+n6*n7+n2*n7+4*ep*n7+ep*n7*x6,4-6*n7+2*n7^2-2*n2+2*n2*n7-8*ep+6*ep
         *n7+2*ep*n2+4*ep^2))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,1,0,n6,n7,1+n8,0)*(rat(5*n8-n8^2-2*n7*
         n8-n6*n8-n2*n8-4*ep*n8-ep*n8*x6,4-6*n7+2*n7^2-2*n2+2*n2*n7-8*ep+6*ep*
         n7+2*ep*n2+4*ep^2))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1+n2,0,1,1,n6,-1+n7,n8,0)*(rat(-5*n2+n2*n8+2*n2
         *n7+n2*n6+n2^2+4*ep*n2+ep*n2*x6,2-3*n7+n7^2-n2+n2*n7-4*ep+3*ep*n7+ep*
         n2+2*ep^2))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,0,1,0,n6,n7,1+n8,0)*(rat(-5*n8+n8^2+2*n7*n8+
         n6*n8+n2*n8+4*ep*n8+ep*n8*x6,4-6*n7+2*n7^2-2*n2+2*n2*n7-8*ep+6*ep*n7+2
         *ep*n2+4*ep^2))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,0,1,1,n6,n7,n8,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,1,0,1,n6,n7,n8,0)*(rat(-8+2*n8+2*n7+n6+3*n2+
         6*ep+ep*x6,-2+n7+n2+2*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,1,0,2,n6,n7,-1+n8,0)*(rat(5-n8-2*n7-n6-6*n2+
         n2*n8+2*n2*n7+n2*n6+n2^2-9*ep-ep*x6+ep*n8+2*ep*n7+ep*n6+5*ep*n2+ep*n2*
         x6+4*ep^2+ep^2*x6,-2+2*n8+3*n7-3*n7*n8-n7^2+n7^2*n8+n2-n2*n8-n2*n7+n2*
         n7*n8+4*ep-4*ep*n8-3*ep*n7+3*ep*n7*n8-ep*n2+ep*n2*n8-2*ep^2+2*ep^2*n8)
         )
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,1,1,0,n6,n7,n8,0)*(rat(-24+9*n8-n8^2+19*n7-3
         *n7*n8-4*n7^2+6*n6-n6*n8-3*n6*n7+14*n2-3*n2*n8-5*n2*n7-2*n2*n6-2*n2^2+
         34*ep+6*ep*x6-6*ep*n8-ep*n8*x6-14*ep*n7-3*ep*n7*x6-4*ep*n6-10*ep*n2-2*
         ep*n2*x6-12*ep^2-4*ep^2*x6,4-6*n7+2*n7^2-2*n2+2*n2*n7-8*ep+6*ep*n7+2*
         ep*n2+4*ep^2));
*--#] no,noints=6,signs=1,+,1,1,1,o,+,+,0,level=1,stop=0 : 
* LHS = Z(1,n2,1,1,1,n6,n7,1,0)
* DEN = den(-1+n2)*den(-2+n7+n2+2*ep)
* MAX = Z(0,0,1,0,0,1,1,0,0)
* MIN = Z(-1,-1,-1,-1,-1,0,-1,-1,0)
* nterms = 14 (Z: 3, Y: 11)
* NonIntegerIndices: 6
* Signs: 1,+!,1,1,1,o,+,1,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,1,n2?{>=2},1,1,1,n6?,n7?pos_,1,0) =
*      #include no.rules # no,noints=6,signs=1,+!,1,1,1,o,+,1,0,level=1,stop=0
*
*--#[ no,noints=6,signs=1,+!,1,1,1,o,+,1,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,1,1,n6,n7,1,0)*(rat(4-6*n7+2*n7^2-n6+n6
         *n7-n2+n2*n7-12*ep-ep*x6+8*ep*n7+ep*n7*x6+2*ep*n6+2*ep*n2+8*ep^2+2*
         ep^2*x6,2-n7-3*n2+n2*n7+n2^2-2*ep+2*ep*n2))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,-1+n2,2,1,1,n6,n7,1,0)*(rat(4-2*n7-n6-n2-4*ep-
         ep*x6,2-n7-3*n2+n2*n7+n2^2-2*ep+2*ep*n2))
      +`Z'(0,0,0,0,0,x6,0,0,0,1,n2,1,1,1,n6,-1+n7,1,0)*(rat(-2+n7+n6+2*ep+ep*x6
         ,-2+n7+n2+2*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,0,1,n6,1+n7,1,0)*(rat(-3*n7,-4+2*n7+2*n2+4
         *ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,0,1+n6,n7,1,0)*(rat(-n6-ep*x6,-4+2*n7+2*
         n2+4*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,0,n6,1+n7,1,0)*(rat(-n7,-4+2*n7+2*n2+4*
         ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,1,1+n6,n7,0,0)*(rat(n6+ep*x6,-4+2*n7+2*
         n2+4*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,n2,1,1,1,n6,n7,1,0)*(rat(2-n7-n6-2*ep-ep*x6,-4+
         2*n7+2*n2+4*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,0,1,n6,1+n7,1,0)*(rat(-6*n7+2*n7^2+n6*
         n7+3*n2*n7+4*ep*n7+ep*n7*x6,2-n7-3*n2+n2*n7+n2^2-2*ep+2*ep*n2))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,-1+n2,1,1,0,n6,1+n7,1,0)*(rat(4*n7-2*n7^2-n6*n7
         -n2*n7-4*ep*n7-ep*n7*x6,2-n7-3*n2+n2*n7+n2^2-2*ep+2*ep*n2))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,-1+n2,2,0,1,n6,n7,1,0)*(rat(-4+2*n7+n6+n2+4*ep+
         ep*x6,2-n7-3*n2+n2*n7+n2^2-2*ep+2*ep*n2))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,0,1,1,n6,n7,1,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,1,0,1,n6,n7,1,0)*(rat(-6+2*n7+n6+3*n2+6*ep+
         ep*x6,-2+n7+n2+2*ep))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,n2,1,1,0,n6,n7,1,0)*(rat(2-n7-n6-2*ep-ep*x6,-2+
         n7+n2+2*ep));
*--#] no,noints=6,signs=1,+!,1,1,1,o,+,1,0,level=1,stop=0 : 
* LHS = Z(1,1,1,1,1,n6,n7,1,0)
* DEN = den(-1+n7)*den(-1+n7+2*ep)*den(-2+n7+n6+3*ep)
* MAX = Z(0,1,0,0,1,1,1,1,0)
* MIN = Z(-1,-1,-1,-1,-1,0,-1,-1,0)
* nterms = 21 (Z: 1, Y: 20)
* NonIntegerIndices: 6
* Signs: 1,1,1,1,1,o,+,1,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,1,1,1,1,1,n6?,n7?{>=2},1,0) =
*      #include no.rules # no,noints=6,signs=1,1,1,1,1,o,+,1,0,level=1,stop=0
*
*--#[ no,noints=6,signs=1,1,1,1,1,o,+,1,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,1,1,1,1,1,n6,-1+n7,1,0)*(rat(-4+4*n7-n7^2+4*n6-2*
         n6*n7-n6^2+12*ep+4*ep*x6-6*ep*n7-2*ep*n7*x6-6*ep*n6-2*ep*n6*x6-8*ep^2-
         6*ep^2*x6-ep^2*x6^2,1-2*n7+n7^2-2*ep+2*ep*n7))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,1,1,0,1,n6,1+n7,1,0)*(rat(2*n7^2-n7^3-2*n6*n7+
         n6^2*n7+ep*n7-2*ep*n7*x6-5*ep*n7^2+4*ep*n6*n7+2*ep*n6*n7*x6+4*ep^2*n7*
         x6+ep^2*n7*x6^2,-4+10*n7-8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+2*
         ep*x6-24*ep*n7-4*ep*n7*x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7-
         12*ep^2-4*ep^2*x6+12*ep^2*n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,1,1,1,0,1+n6,n7,1,0)*(rat(-2*n6+3*n6*n7-n6*n7^2
         +n6^2-n6^2*n7-2*ep*x6+3*ep*n7*x6-ep*n7^2*x6+9*ep*n6+2*ep*n6*x6-7*ep*n6
         *n7-2*ep*n6*n7*x6-2*ep*n6^2+9*ep^2*x6+ep^2*x6^2-7*ep^2*n7*x6-ep^2*n7*
         x6^2-8*ep^2*n6-4*ep^2*n6*x6-8*ep^3*x6-2*ep^3*x6^2,-4+10*n7-8*n7^2+2*
         n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+2*ep*x6-24*ep*n7-4*ep*n7*x6+10*ep*
         n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7-12*ep^2-4*ep^2*x6+12*ep^2*n7+4*
         ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,1,1,1,0,n6,1+n7,1,0)*(rat(-2*n7+3*n7^2-n7^3+n6*
         n7-n6*n7^2+9*ep*n7+ep*n7*x6-7*ep*n7^2-ep*n7^2*x6-2*ep*n6*n7-8*ep^2*n7-
         2*ep^2*n7*x6,-4+10*n7-8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+2*ep*
         x6-24*ep*n7-4*ep*n7*x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7-12*
         ep^2-4*ep^2*x6+12*ep^2*n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,1,1,1,1,1+n6,n7,0,0)*(rat(2*n6-3*n6*n7+n6*n7^2-
         n6^2+n6^2*n7+2*ep*x6-3*ep*n7*x6+ep*n7^2*x6-9*ep*n6-2*ep*n6*x6+7*ep*n6*
         n7+2*ep*n6*n7*x6+2*ep*n6^2-9*ep^2*x6-ep^2*x6^2+7*ep^2*n7*x6+ep^2*n7*
         x6^2+8*ep^2*n6+4*ep^2*n6*x6+8*ep^3*x6+2*ep^3*x6^2,-4+10*n7-8*n7^2+2*
         n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+2*ep*x6-24*ep*n7-4*ep*n7*x6+10*ep*
         n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7-12*ep^2-4*ep^2*x6+12*ep^2*n7+4*
         ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,1,1,1,1,n6,n7,1,0)*(rat(-2+5*n7-4*n7^2+n7^3+n6-
         2*n6*n7+n6*n7^2+16*ep+ep*x6-21*ep*n7-2*ep*n7*x6+7*ep*n7^2+ep*n7^2*x6-
         11*ep*n6+7*ep*n6*n7+2*ep*n6^2-30*ep^2-11*ep^2*x6+18*ep^2*n7+7*ep^2*n7*
         x6+12*ep^2*n6+4*ep^2*n6*x6+16*ep^3+12*ep^3*x6+2*ep^3*x6^2,-4+10*n7-8*
         n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+2*ep*x6-24*ep*n7-4*ep*n7*x6+
         10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7-12*ep^2-4*ep^2*x6+12*ep^2*
         n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,2,0,1,1,n6,n7,1,0)*(rat(-6+7*n7-2*n7^2+5*n6-3*
         n6*n7-n6^2+14*ep+5*ep*x6-8*ep*n7-3*ep*n7*x6-6*ep*n6-2*ep*n6*x6-8*ep^2-
         6*ep^2*x6-ep^2*x6^2,-4+10*n7-8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*
         ep+2*ep*x6-24*ep*n7-4*ep*n7*x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6
         *n7-12*ep^2-4*ep^2*x6+12*ep^2*n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,2,1,1,1,n6,n7,0,0)*(rat(2-3*n7+n7^2-n6+n6*n7-9*
         ep-ep*x6+7*ep*n7+ep*n7*x6+2*ep*n6+8*ep^2+2*ep^2*x6,-2+5*n7-4*n7^2+n7^3
         +n6-2*n6*n7+n6*n7^2+7*ep+ep*x6-12*ep*n7-2*ep*n7*x6+5*ep*n7^2+ep*n7^2*
         x6-2*ep*n6+2*ep*n6*n7-6*ep^2-2*ep^2*x6+6*ep^2*n7+2*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,0,1,0,1,n6,1+n7,1,0)*(rat(-2*n7+n7^2+3*n6*n7-n6
         *n7^2-n6^2*n7+8*ep*n7+3*ep*n7*x6-2*ep*n7^2-ep*n7^2*x6-6*ep*n6*n7-2*ep*
         n6*n7*x6-8*ep^2*n7-6*ep^2*n7*x6-ep^2*n7*x6^2,-2+5*n7-4*n7^2+n7^3+n6-2*
         n6*n7+n6*n7^2+7*ep+ep*x6-12*ep*n7-2*ep*n7*x6+5*ep*n7^2+ep*n7^2*x6-2*ep
         *n6+2*ep*n6*n7-6*ep^2-2*ep^2*x6+6*ep^2*n7+2*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,0,1,1,0,n6,1+n7,1,0)*(rat(6*n7-7*n7^2+2*n7^3-5*
         n6*n7+3*n6*n7^2+n6^2*n7-26*ep*n7-5*ep*n7*x6+16*ep*n7^2+3*ep*n7^2*x6+10
         *ep*n6*n7+2*ep*n6*n7*x6+24*ep^2*n7+10*ep^2*n7*x6+ep^2*n7*x6^2,-4+10*n7
         -8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+2*ep*x6-24*ep*n7-4*ep*n7*
         x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7-12*ep^2-4*ep^2*x6+12*
         ep^2*n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,0,1,1,0,n6,n7,2,0)*(rat(-6+7*n7-2*n7^2+5*n6-3*
         n6*n7-n6^2+14*ep+5*ep*x6-8*ep*n7-3*ep*n7*x6-6*ep*n6-2*ep*n6*x6-8*ep^2-
         6*ep^2*x6-ep^2*x6^2,-4+10*n7-8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*
         ep+2*ep*x6-24*ep*n7-4*ep*n7*x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6
         *n7-12*ep^2-4*ep^2*x6+12*ep^2*n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,0,1,1,1,n6,n7,1,0)*(rat(3-2*n7-4*n6+2*n6*n7+
         n6^2-10*ep-4*ep*x6+4*ep*n7+2*ep*n7*x6+6*ep*n6+2*ep*n6*x6+8*ep^2+6*ep^2
         *x6+ep^2*x6^2,1-2*n7+n7^2-2*ep+2*ep*n7))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,0,1,0,n6,n7,2,0)*(rat(6-7*n7+2*n7^2-5*n6+3*n6
         *n7+n6^2-14*ep-5*ep*x6+8*ep*n7+3*ep*n7*x6+6*ep*n6+2*ep*n6*x6+8*ep^2+6*
         ep^2*x6+ep^2*x6^2,-4+10*n7-8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+
         2*ep*x6-24*ep*n7-4*ep*n7*x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7
         -12*ep^2-4*ep^2*x6+12*ep^2*n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,0,1,1,n6,n7,1,0)*(rat(-2+5*n7-4*n7^2+n7^3+n6-
         2*n6*n7+n6*n7^2-2*ep+ep*x6-2*ep*n7*x6+ep*n7^2+ep*n7^2*x6+4*ep*n6-2*ep*
         n6*n7-ep*n6^2+18*ep^2+4*ep^2*x6-10*ep^2*n7-2*ep^2*n7*x6-8*ep^2*n6-2*
         ep^2*n6*x6-16*ep^3-8*ep^3*x6-ep^3*x6^2,-2+5*n7-4*n7^2+n7^3+n6-2*n6*n7+
         n6*n7^2+7*ep+ep*x6-12*ep*n7-2*ep*n7*x6+5*ep*n7^2+ep*n7^2*x6-2*ep*n6+2*
         ep*n6*n7-6*ep^2-2*ep^2*x6+6*ep^2*n7+2*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,1,0,1,n6,n7,1,0)*(rat(-6+13*n7-9*n7^2+2*n7^3+
         5*n6-8*n6*n7+3*n6*n7^2-n6^2+n6^2*n7+27*ep+5*ep*x6-40*ep*n7-8*ep*n7*x6+
         14*ep*n7^2+3*ep*n7^2*x6-14*ep*n6-2*ep*n6*x6+12*ep*n6*n7+2*ep*n6*n7*x6+
         ep*n6^2-38*ep^2-14*ep^2*x6-ep^2*x6^2+30*ep^2*n7+12*ep^2*n7*x6+ep^2*n7*
         x6^2+8*ep^2*n6+2*ep^2*n6*x6+16*ep^3+8*ep^3*x6+ep^3*x6^2,-2+5*n7-4*n7^2
         +n7^3+n6-2*n6*n7+n6*n7^2+7*ep+ep*x6-12*ep*n7-2*ep*n7*x6+5*ep*n7^2+ep*
         n7^2*x6-2*ep*n6+2*ep*n6*n7-6*ep^2-2*ep^2*x6+6*ep^2*n7+2*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,1,0,2,n6,n7,0,0)*(rat(-2+n7+3*n6-n6*n7-n6^2+8
         *ep+3*ep*x6-2*ep*n7-ep*n7*x6-6*ep*n6-2*ep*n6*x6-8*ep^2-6*ep^2*x6-ep^2*
         x6^2,-4+10*n7-8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+2*ep*x6-24*ep
         *n7-4*ep*n7*x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7-12*ep^2-4*
         ep^2*x6+12*ep^2*n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,1,1,0,n6,n7,1,0)*(rat(14-29*n7+19*n7^2-4*n7^3
         -13*n6+20*n6*n7-7*n6*n7^2+3*n6^2-3*n6^2*n7-58*ep-13*ep*x6+80*ep*n7+20*
         ep*n7*x6-26*ep*n7^2-7*ep*n7^2*x6+36*ep*n6+6*ep*n6*x6-28*ep*n6*n7-6*ep*
         n6*n7*x6-4*ep*n6^2+76*ep^2+36*ep^2*x6+3*ep^2*x6^2-52*ep^2*n7-28*ep^2*
         n7*x6-3*ep^2*n7*x6^2-24*ep^2*n6-8*ep^2*n6*x6-32*ep^3-24*ep^3*x6-4*ep^3
         *x6^2,-4+10*n7-8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+2*ep*x6-24*
         ep*n7-4*ep*n7*x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7-12*ep^2-4*
         ep^2*x6+12*ep^2*n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,1,1,1,n6,n7,0,0)*(rat(-2+3*n7-n7^2+3*n6-4*n6*
         n7+n6*n7^2-n6^2+n6^2*n7+13*ep+3*ep*x6-13*ep*n7-4*ep*n7*x6+2*ep*n7^2+ep
         *n7^2*x6-13*ep*n6-2*ep*n6*x6+9*ep*n6*n7+2*ep*n6*n7*x6+2*ep*n6^2-26*
         ep^2-13*ep^2*x6-ep^2*x6^2+14*ep^2*n7+9*ep^2*n7*x6+ep^2*n7*x6^2+12*ep^2
         *n6+4*ep^2*n6*x6+16*ep^3+12*ep^3*x6+2*ep^3*x6^2,-2+5*n7-4*n7^2+n7^3+n6
         -2*n6*n7+n6*n7^2+7*ep+ep*x6-12*ep*n7-2*ep*n7*x6+5*ep*n7^2+ep*n7^2*x6-2
         *ep*n6+2*ep*n6*n7-6*ep^2-2*ep^2*x6+6*ep^2*n7+2*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,1,1,2,n6,-1+n7,0,0)*(rat(2-n7-3*n6+n6*n7+n6^2
         -8*ep-3*ep*x6+2*ep*n7+ep*n7*x6+6*ep*n6+2*ep*n6*x6+8*ep^2+6*ep^2*x6+
         ep^2*x6^2,-4+10*n7-8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*n7^2+14*ep+2*ep*x6-
         24*ep*n7-4*ep*n7*x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6+4*ep*n6*n7-12*ep^2
         -4*ep^2*x6+12*ep^2*n7+4*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,2,0,1,1,n6,-1+n7,1,0)*(rat(6-7*n7+2*n7^2-5*n6+3
         *n6*n7+n6^2-14*ep-5*ep*x6+8*ep*n7+3*ep*n7*x6+6*ep*n6+2*ep*n6*x6+8*ep^2
         +6*ep^2*x6+ep^2*x6^2,-2+5*n7-4*n7^2+n7^3+n6-2*n6*n7+n6*n7^2+7*ep+ep*x6
         -12*ep*n7-2*ep*n7*x6+5*ep*n7^2+ep*n7^2*x6-2*ep*n6+2*ep*n6*n7-6*ep^2-2*
         ep^2*x6+6*ep^2*n7+2*ep^2*n7*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,2,1,1,1,n6,-1+n7,0,0)*(rat(-10+13*n7-4*n7^2+7*
         n6-5*n6*n7-n6^2+32*ep+7*ep*x6-22*ep*n7-5*ep*n7*x6-10*ep*n6-2*ep*n6*x6-
         24*ep^2-10*ep^2*x6-ep^2*x6^2,-4+10*n7-8*n7^2+2*n7^3+2*n6-4*n6*n7+2*n6*
         n7^2+14*ep+2*ep*x6-24*ep*n7-4*ep*n7*x6+10*ep*n7^2+2*ep*n7^2*x6-4*ep*n6
         +4*ep*n6*n7-12*ep^2-4*ep^2*x6+12*ep^2*n7+4*ep^2*n7*x6));
*--#] no,noints=6,signs=1,1,1,1,1,o,+,1,0,level=1,stop=0 : 
* LHS = Z(1,1,1,1,1,n6,1,1,0)
* DEN = den(-1+n6)
* MAX = Z(0,0,0,0,1,0,0,1,0)
* MIN = Z(-1,-1,-1,-1,0,-1,-1,0,0)
* nterms = 6 (Z: 1, Y: 5)
* NonIntegerIndices: 6
* Signs: 1,1,1,1,1,+,1,1,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,1,1,1,1,1,n6?{>=2},1,1,0) =
*      #include no.rules # no,noints=6,signs=1,1,1,1,1,+,1,1,0,level=1,stop=0
*
*--#[ no,noints=6,signs=1,1,1,1,1,+,1,1,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,1,1,1,1,1,-1+n6,1,1,0)*(rat(-1+n6+4*ep+ep*x6,-1+
         n6+ep*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,1,1,1,1,n6,1,1,0)*(rat(1,1))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,0,1,1,1,-1+n6,1,2,0)*(rat(1,-1+n6+ep*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,0,1,1,-1+n6,1,2,0)*(rat(-1,-1+n6+ep*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,1,0,2,-1+n6,1,1,0)*(rat(-1,-1+n6+ep*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,1,1,2,-1+n6,0,1,0)*(rat(1,-1+n6+ep*x6));
*--#] no,noints=6,signs=1,1,1,1,1,+,1,1,0,level=1,stop=0 : 
* LHS = Z(1,1,1,1,1,n6,1,1,0)
* DEN = den(n6+4*ep)
* MAX = Z(0,0,0,0,1,1,0,1,0)
* MIN = Z(-1,-1,-1,-1,0,0,-1,0,0)
* nterms = 6 (Z: 1, Y: 5)
* NonIntegerIndices: 6
* Signs: 1,1,1,1,1,-,1,1,0
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,x6?!{0,},0,0,0,1,1,1,1,1,n6?neg0_,1,1,0) =
*      #include no.rules # no,noints=6,signs=1,1,1,1,1,-,1,1,0,level=1,stop=0
*
*--#[ no,noints=6,signs=1,1,1,1,1,-,1,1,0,level=1,stop=0 :
      +`Z'(0,0,0,0,0,x6,0,0,0,1,1,1,1,1,1+n6,1,1,0)*(rat(n6+ep*x6,n6+4*ep+ep*x6
         ))
      +`Y'(0,0,0,0,0,x6,0,0,0,0,1,1,1,1,1+n6,1,1,0)*(rat(-n6-ep*x6,n6+4*ep+ep*
         x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,0,1,1,1,n6,1,2,0)*(rat(-1,n6+4*ep+ep*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,0,1,1,n6,1,2,0)*(rat(1,n6+4*ep+ep*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,1,0,2,n6,1,1,0)*(rat(1,n6+4*ep+ep*x6))
      +`Y'(0,0,0,0,0,x6,0,0,0,1,1,1,1,2,n6,0,1,0)*(rat(-1,n6+4*ep+ep*x6));
*--#] no,noints=6,signs=1,1,1,1,1,-,1,1,0,level=1,stop=0 : 
*--#] no,6-rules
        goto 2;
label 1;
        #call ReduceZ2Y(Z,Y,9,ifmatch->2,1,2,3,4,5,7,8)
        $doloop = 1;
label 2;
        ModuleOption maximum, $doloop;
        .sort:`TOPO'-`nloop++';
        #if `$doloop'
          #redefine loop "0"
        #endif
      #enddo
      #break
*--#] no,6
*--#[ no,8
    #case no,8
      #define nloop "1"
      #do loop=1,1
        #$doloop = 0;
        if (count(`Z',1) == 0) goto 2;
*--#[ no,8-rules
* LHS = Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* DEN = den(-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
* MAX = Z(1,0,1,0,0,0,0,0,2)
* MIN = Z(0,-1,0,0,0,0,-1,-1,1)
* nterms = 6 (Z: 6, Y: 0)
* NonIntegerIndices: 8
* Signs: x,x,x,o,o,o,x,x,-*?neg_
* Order: 1,2,3,4,5,6,7,8,9
* Seed: 0,0,0,0,0,0,0,0,0
* Level: 2
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,0,0,0,x8?!{0,},0,n1?,n2?,n3?,n4?,n5?,n6?,n7?,n8?,n9?neg_) =
*      #include no.rules # no,noints=8,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0
*
*--#[ no,noints=8,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0 :
      +`Z'(0,0,0,0,0,0,0,x8,0,1+n1,-1+n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(n1,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x8))
      +`Z'(0,0,0,0,0,0,0,x8,0,1+n1,n2,n3,n4,n5,n6,-1+n7,n8,1+n9)*(rat(-n1,-14+2
         *n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x8))
      +`Z'(0,0,0,0,0,0,0,x8,0,n1,-1+n2,1+n3,n4,n5,n6,n7,n8,1+n9)*(rat(n3,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x8))
      +`Z'(0,0,0,0,0,0,0,x8,0,n1,-1+n2,n3,n4,n5,n6,n7,n8,2+n9)*(rat(1+n9,-14+2*
         n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x8))
      +`Z'(0,0,0,0,0,0,0,x8,0,n1,n2,1+n3,n4,n5,n6,n7,-1+n8,1+n9)*(rat(-n3,-14+2
         *n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x8))
      +`Z'(0,0,0,0,0,0,0,x8,0,n1,n2,n3,n4,n5,n6,n7,n8,1+n9)*(rat(-3+n9+n3+2*n2+
         n1+2*ep,-14+2*n9+2*n8+2*n7+2*n6+2*n5+2*n4+2*n3+2*n2+2*n1+8*ep+2*ep*x8)
         );
*--#] no,noints=8,signs=x,x,x,o,o,o,x,x,-*?neg_,level=2,stop=0 : 
*--#] no,8-rules
        goto 2;
label 1;
        #call ReduceZ2Y(Z,Y,9,ifmatch->2,1,2,3,4,5,6,7)
        $doloop = 1;
label 2;
        ModuleOption maximum, $doloop;
        .sort:`TOPO'-`nloop++';
        #if `$doloop'
          #redefine loop "0"
        #endif
      #enddo
      #break
*--#] no,8
*--#[ t1,5
    #case t1,5
      #define nloop "1"
      #do loop=1,1
        #$doloop = 0;
        if (count(`Z',1) == 0) goto 2;
*--#[ t1,5-rules
* LHS = Z(n1,n2,n3,n4,n5)
* DEN = den(-1+n4)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,0,-1,0)
* nterms = 4 (Z: 4, Y: 0)
* NonIntegerIndices: 5
* Signs: +,+,+,+,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?) =
*      #include t1.rules # t1,noints=5,signs=+,+,+,+,o,level=1,stop=0
*
*--#[ t1,noints=5,signs=+,+,+,+,o,level=1,stop=0 :
      +`Z'(0,0,0,0,x5,-1+n1,n2,n3,-1+n4,1+n5)*(rat(n5+ep*x5,-1+n4))
      +`Z'(0,0,0,0,x5,-1+n1,n2,n3,n4,n5)*(rat(1,1))
      +`Z'(0,0,0,0,x5,n1,-1+n2,n3,-1+n4,1+n5)*(rat(-n5-ep*x5,-1+n4))
      +`Z'(0,0,0,0,x5,n1,n2,n3,-1+n4,n5)*(rat(-5+n5+n4+2*n1+2*ep+ep*x5,-1+n4));
*--#] t1,noints=5,signs=+,+,+,+,o,level=1,stop=0 : 
* LHS = Z(n1,n2,n3,1,n5)
* DEN = den(-1+n1)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,0,-1,-1,0)
* nterms = 4 (Z: 2, Y: 2)
* NonIntegerIndices: 5
* Signs: +,+,+,1,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},n1?{>=2},n2?pos_,n3?pos_,1,n5?) =
*      #include t1.rules # t1,noints=5,signs=+,+,+,1,o,level=1,stop=0
*
*--#[ t1,noints=5,signs=+,+,+,1,o,level=1,stop=0 :
      +`Z'(0,0,0,0,x5,-1+n1,n2,-1+n3,1,1+n5)*(rat(-n5-ep*x5,-1+n1))
      +`Z'(0,0,0,0,x5,-1+n1,n2,n3,1,n5)*(rat(-3+n5+n1+2*ep+ep*x5,-1+n1))
      +`Y'(0,0,0,0,x5,-1+n1,n2,n3,0,1+n5)*(rat(n5+ep*x5,-1+n1))
      +`Y'(0,0,0,0,x5,n1,n2,n3,0,n5)*(rat(1,1));
*--#] t1,noints=5,signs=+,+,+,1,o,level=1,stop=0 : 
* LHS = Z(1,n2,n3,1,n5)
* DEN = den(-1+n2)*den(-1+n3)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,-1,-1,1)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 5
* Signs: 1,+,+,1,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},1,n2?{>=2},n3?{>=2},1,n5?) =
*      #include t1.rules # t1,noints=5,signs=1,+,+,1,o,level=1,stop=1
*
*--#[ t1,noints=5,signs=1,+,+,1,o,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,1,-1+n2,-1+n3,1,1+n5)*(rat(-12+7*n5-n5^2+7*n3-2*n3*n5-
         n3^2+7*n2-2*n2*n5-2*n2*n3-n2^2+10*ep+7*ep*x5-3*ep*n5-2*ep*n5*x5-3*ep*
         n3-2*ep*n3*x5-3*ep*n2-2*ep*n2*x5-2*ep^2-3*ep^2*x5-ep^2*x5^2,1-n3-n2+n2
         *n3))
      +`Y'(0,0,0,0,x5,0,n2,-1+n3,1,1+n5)*(rat(-3+n3+n2+ep,-1+n3))
      +`Y'(0,0,0,0,x5,1,-1+n2,n3,0,1+n5)*(rat(-3+n3+n2+ep,-1+n2));
*--#] t1,noints=5,signs=1,+,+,1,o,level=1,stop=1 : 
* LHS = Z(1,1,n3,1,n5)
* DEN = den(-1+n3)*den(-2+n5+n3+ep)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,-1,-1,0)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 5
* Signs: 1,1,+,1,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},1,1,n3?{>=2},1,n5?) =
*      #include t1.rules # t1,noints=5,signs=1,1,+,1,o,level=1,stop=1
*
*--#[ t1,noints=5,signs=1,1,+,1,o,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,1,1,-1+n3,1,n5)*(rat(-3+n5+n3+2*ep+ep*x5,-1+n3))
      +`Y'(0,0,0,0,x5,0,1,-1+n3,1,1+n5)*(rat(-n5-ep*x5,-1+n3))
      +`Y'(0,0,0,0,x5,1,0,n3,0,1+n5)*(rat(n5+ep*x5,-2+n5+n3+ep+ep*x5));
*--#] t1,noints=5,signs=1,1,+,1,o,level=1,stop=1 : 
* LHS = Z(1,n2,1,1,n5)
* DEN = den(-1+n2)*den(-2+n5+n2+ep)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,-1,-1,0)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 5
* Signs: 1,+,1,1,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},1,n2?{>=2},1,1,n5?) =
*      #include t1.rules # t1,noints=5,signs=1,+,1,1,o,level=1,stop=1
*
*--#[ t1,noints=5,signs=1,+,1,1,o,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,1,-1+n2,1,1,n5)*(rat(-3+n5+n2+2*ep+ep*x5,-1+n2))
      +`Y'(0,0,0,0,x5,0,n2,0,1,1+n5)*(rat(n5+ep*x5,-2+n5+n2+ep+ep*x5))
      +`Y'(0,0,0,0,x5,1,-1+n2,1,0,1+n5)*(rat(-n5-ep*x5,-1+n2));
*--#] t1,noints=5,signs=1,+,1,1,o,level=1,stop=1 : 
* LHS = Z(1,1,1,1,n5)
* DEN = den(-1+n5+ep)*den(-2+n5+2*ep)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,-1,-1,-1)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 5
* Signs: 1,1,1,1,+
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},1,1,1,1,n5?{>=2}) =
*      #include t1.rules # t1,noints=5,signs=1,1,1,1,+,level=1,stop=1
*
*--#[ t1,noints=5,signs=1,1,1,1,+,level=1,stop=1 :
      +`Z'(0,0,0,0,x5,1,1,1,1,-1+n5)*(rat(2-n5-2*ep-ep*x5,-1+n5+ep+ep*x5))
      +`Y'(0,0,0,0,x5,0,1,0,1,1+n5)*(rat(-3*n5+2*n5^2-3*ep*x5+3*ep*n5+4*ep*n5*
         x5+3*ep^2*x5+2*ep^2*x5^2,2-3*n5+n5^2-4*ep-3*ep*x5+3*ep*n5+2*ep*n5*x5+2
         *ep^2+3*ep^2*x5+ep^2*x5^2))
      +`Y'(0,0,0,0,x5,1,0,1,0,1+n5)*(rat(-3*n5+2*n5^2-3*ep*x5+3*ep*n5+4*ep*n5*
         x5+3*ep^2*x5+2*ep^2*x5^2,2-3*n5+n5^2-4*ep-3*ep*x5+3*ep*n5+2*ep*n5*x5+2
         *ep^2+3*ep^2*x5+ep^2*x5^2));
*--#] t1,noints=5,signs=1,1,1,1,+,level=1,stop=1 : 
* LHS = Z(1,1,1,1,n5)
* DEN = den(-1+n5+2*ep)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,0,0,1)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 5
* Signs: 1,1,1,1,-
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 2
* Stop: 2
*
    id `IDSPEC' `Z'(0,0,0,0,x5?!{0,},1,1,1,1,n5?neg0_) =
*      #include t1.rules # t1,noints=5,signs=1,1,1,1,-,level=2,stop=2
*
*--#[ t1,noints=5,signs=1,1,1,1,-,level=2,stop=2 :
      +`Z'(0,0,0,0,x5,1,1,1,1,1+n5)*(rat(-n5-ep-ep*x5,-1+n5+2*ep+ep*x5))
      +`Y'(0,0,0,0,x5,0,1,1,1,1+n5)*(rat(-1+2*n5+3*ep+2*ep*x5,-1+n5+2*ep+ep*x5)
         )
      +`Y'(0,0,0,0,x5,1,0,1,1,1+n5)*(rat(-1+2*n5+3*ep+2*ep*x5,-1+n5+2*ep+ep*x5)
         );
*--#] t1,noints=5,signs=1,1,1,1,-,level=2,stop=2 : 
*--#] t1,5-rules
        goto 2;
label 1;
        #call ReduceZ2Y(Z,Y,5,ifmatch->2,1,2,3,4)
        $doloop = 1;
label 2;
        ModuleOption maximum, $doloop;
        .sort:`TOPO'-`nloop++';
        #if `$doloop'
          #redefine loop "0"
        #endif
      #enddo
      #break
*--#] t1,5
*--#[ t1,24
    #case t1,24
      #define nloop "1"
      #do loop=1,1
        #$doloop = 0;
        if (count(`Z',1) == 0) goto 2;
*--#[ t1,24-rules
* LHS = Z(n1,n2,n3,n4,n5)
* DEN = den(-1+n1)*den(-3+n5+n2+n1+ep)*den(-4+n5+n2+n1+2*ep)
* MAX = Z(0,0,0,1,0)
* MIN = Z(-1,0,-1,0,-1)
* nterms = 4 (Z: 4, Y: 0)
* NonIntegerIndices: 2,4
* Signs: +,o,+,o,+
* Order: 1,2,4,5,3
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 2
*
    id `IDSPEC' `Z'(0,x2?!{0,},0,x4?!{0,},0,n1?{>=2},n2?,n3?pos_,n4?,n5?pos_) =
*      #include t1.rules # t1,noints=2,4,signs=+,o,+,o,+,order=1,2,4,5,3,level=1,stop=2
*
*--#[ t1,noints=2,4,signs=+,o,+,o,+,order=1,2,4,5,3,level=1,stop=2 :
      +`Z'(0,x2,0,x4,0,-1+n1,n2,-1+n3,1+n4,n5)*(rat(-28*n4+15*n4*n5-2*n4*n5^2+4
         *n4^2-n4^2*n5+11*n3*n4-3*n3*n4*n5-n3*n4^2-n3^2*n4+11*n2*n4-3*n2*n4*n5-
         n2*n4^2-2*n2*n3*n4-n2^2*n4+4*n1*n4-n1*n4*n5-n1*n3*n4-n1*n2*n4-28*ep*x4
         +15*ep*n5*x4-2*ep*n5^2*x4+26*ep*n4+8*ep*n4*x4+11*ep*n4*x2-7*ep*n4*n5-2
         *ep*n4*n5*x4-3*ep*n4*n5*x2-2*ep*n4^2-ep*n4^2*x2+11*ep*n3*x4-3*ep*n3*n5
         *x4-5*ep*n3*n4-2*ep*n3*n4*x4-2*ep*n3*n4*x2-ep*n3^2*x4+11*ep*n2*x4-3*ep
         *n2*n5*x4-5*ep*n2*n4-2*ep*n2*n4*x4-2*ep*n2*n4*x2-2*ep*n2*n3*x4-ep*n2^2
         *x4+4*ep*n1*x4-ep*n1*n5*x4-2*ep*n1*n4-ep*n1*n4*x2-ep*n1*n3*x4-ep*n1*n2
         *x4+26*ep^2*x4+4*ep^2*x4^2+11*ep^2*x2*x4-7*ep^2*n5*x4-ep^2*n5*x4^2-3*
         ep^2*n5*x2*x4-6*ep^2*n4-4*ep^2*n4*x4-5*ep^2*n4*x2-2*ep^2*n4*x2*x4-ep^2
         *n4*x2^2-5*ep^2*n3*x4-ep^2*n3*x4^2-2*ep^2*n3*x2*x4-5*ep^2*n2*x4-ep^2*
         n2*x4^2-2*ep^2*n2*x2*x4-2*ep^2*n1*x4-ep^2*n1*x2*x4-6*ep^3*x4-2*ep^3*
         x4^2-5*ep^3*x2*x4-ep^3*x2*x4^2-ep^3*x2^2*x4,-12+7*n5-n5^2+7*n2-2*n2*n5
         -n2^2+19*n1-9*n1*n5+n1*n5^2-9*n1*n2+2*n1*n2*n5+n1*n2^2-8*n1^2+2*n1^2*
         n5+2*n1^2*n2+n1^3+10*ep+7*ep*x2-3*ep*n5-2*ep*n5*x2-3*ep*n2-2*ep*n2*x2-
         13*ep*n1-9*ep*n1*x2+3*ep*n1*n5+2*ep*n1*n5*x2+3*ep*n1*n2+2*ep*n1*n2*x2+
         3*ep*n1^2+2*ep*n1^2*x2-2*ep^2-3*ep^2*x2-ep^2*x2^2+2*ep^2*n1+3*ep^2*n1*
         x2+ep^2*n1*x2^2))
      +`Z'(0,x2,0,x4,0,-1+n1,n2,n3,1+n4,-1+n5)*(rat(28*n4-15*n4*n5+2*n4*n5^2-4*
         n4^2+n4^2*n5-11*n3*n4+3*n3*n4*n5+n3*n4^2+n3^2*n4-11*n2*n4+3*n2*n4*n5+
         n2*n4^2+2*n2*n3*n4+n2^2*n4-4*n1*n4+n1*n4*n5+n1*n3*n4+n1*n2*n4+28*ep*x4
         -15*ep*n5*x4+2*ep*n5^2*x4-26*ep*n4-8*ep*n4*x4-11*ep*n4*x2+7*ep*n4*n5+2
         *ep*n4*n5*x4+3*ep*n4*n5*x2+2*ep*n4^2+ep*n4^2*x2-11*ep*n3*x4+3*ep*n3*n5
         *x4+5*ep*n3*n4+2*ep*n3*n4*x4+2*ep*n3*n4*x2+ep*n3^2*x4-11*ep*n2*x4+3*ep
         *n2*n5*x4+5*ep*n2*n4+2*ep*n2*n4*x4+2*ep*n2*n4*x2+2*ep*n2*n3*x4+ep*n2^2
         *x4-4*ep*n1*x4+ep*n1*n5*x4+2*ep*n1*n4+ep*n1*n4*x2+ep*n1*n3*x4+ep*n1*n2
         *x4-26*ep^2*x4-4*ep^2*x4^2-11*ep^2*x2*x4+7*ep^2*n5*x4+ep^2*n5*x4^2+3*
         ep^2*n5*x2*x4+6*ep^2*n4+4*ep^2*n4*x4+5*ep^2*n4*x2+2*ep^2*n4*x2*x4+ep^2
         *n4*x2^2+5*ep^2*n3*x4+ep^2*n3*x4^2+2*ep^2*n3*x2*x4+5*ep^2*n2*x4+ep^2*
         n2*x4^2+2*ep^2*n2*x2*x4+2*ep^2*n1*x4+ep^2*n1*x2*x4+6*ep^3*x4+2*ep^3*
         x4^2+5*ep^3*x2*x4+ep^3*x2*x4^2+ep^3*x2^2*x4,-12+7*n5-n5^2+7*n2-2*n2*n5
         -n2^2+19*n1-9*n1*n5+n1*n5^2-9*n1*n2+2*n1*n2*n5+n1*n2^2-8*n1^2+2*n1^2*
         n5+2*n1^2*n2+n1^3+10*ep+7*ep*x2-3*ep*n5-2*ep*n5*x2-3*ep*n2-2*ep*n2*x2-
         13*ep*n1-9*ep*n1*x2+3*ep*n1*n5+2*ep*n1*n5*x2+3*ep*n1*n2+2*ep*n1*n2*x2+
         3*ep*n1^2+2*ep*n1^2*x2-2*ep^2-3*ep^2*x2-ep^2*x2^2+2*ep^2*n1+3*ep^2*n1*
         x2+ep^2*n1*x2^2))
      +`Z'(0,x2,0,x4,0,-1+n1,n2,n3,n4,n5)*(rat(-168+146*n5-42*n5^2+4*n5^3+59*n4
         -31*n4*n5+4*n4*n5^2-5*n4^2+n4^2*n5+59*n3-38*n3*n5+6*n3*n5^2-17*n3*n4+5
         *n3*n4*n5+n3*n4^2-5*n3^2+2*n3^2*n5+n3^2*n4+59*n2-38*n2*n5+6*n2*n5^2-17
         *n2*n4+5*n2*n4*n5+n2*n4^2-10*n2*n3+4*n2*n3*n5+2*n2*n3*n4-5*n2^2+2*n2^2
         *n5+n2^2*n4+87*n1-46*n1*n5+6*n1*n5^2-21*n1*n4+5*n1*n4*n5+n1*n4^2-21*n1
         *n3+6*n1*n3*n5+3*n1*n3*n4+n1*n3^2-21*n1*n2+6*n1*n2*n5+3*n1*n2*n4+2*n1*
         n2*n3+n1*n2^2-16*n1^2+4*n1^2*n5+2*n1^2*n4+2*n1^2*n3+2*n1^2*n2+n1^3+212
         *ep+59*ep*x4+59*ep*x2-124*ep*n5-31*ep*n5*x4-38*ep*n5*x2+18*ep*n5^2+4*
         ep*n5^2*x4+6*ep*n5^2*x2-49*ep*n4-10*ep*n4*x4-17*ep*n4*x2+13*ep*n4*n5+2
         *ep*n4*n5*x4+5*ep*n4*n5*x2+2*ep*n4^2+ep*n4^2*x2-49*ep*n3-17*ep*n3*x4-
         10*ep*n3*x2+16*ep*n3*n5+5*ep*n3*n5*x4+4*ep*n3*n5*x2+7*ep*n3*n4+2*ep*n3
         *n4*x4+2*ep*n3*n4*x2+2*ep*n3^2+ep*n3^2*x4-49*ep*n2-17*ep*n2*x4-10*ep*
         n2*x2+16*ep*n2*n5+5*ep*n2*n5*x4+4*ep*n2*n5*x2+7*ep*n2*n4+2*ep*n2*n4*x4
         +2*ep*n2*n4*x2+4*ep*n2*n3+2*ep*n2*n3*x4+2*ep*n2^2+ep*n2^2*x4-75*ep*n1-
         21*ep*n1*x4-21*ep*n1*x2+20*ep*n1*n5+5*ep*n1*n5*x4+6*ep*n1*n5*x2+9*ep*
         n1*n4+2*ep*n1*n4*x4+3*ep*n1*n4*x2+9*ep*n1*n3+3*ep*n1*n3*x4+2*ep*n1*n3*
         x2+9*ep*n1*n2+3*ep*n1*n2*x4+2*ep*n1*n2*x2+7*ep*n1^2+2*ep*n1^2*x4+2*ep*
         n1^2*x2-88*ep^2-49*ep^2*x4-5*ep^2*x4^2-49*ep^2*x2-17*ep^2*x2*x4-5*ep^2
         *x2^2+26*ep^2*n5+13*ep^2*n5*x4+ep^2*n5*x4^2+16*ep^2*n5*x2+5*ep^2*n5*x2
         *x4+2*ep^2*n5*x2^2+10*ep^2*n4+4*ep^2*n4*x4+7*ep^2*n4*x2+2*ep^2*n4*x2*
         x4+ep^2*n4*x2^2+10*ep^2*n3+7*ep^2*n3*x4+ep^2*n3*x4^2+4*ep^2*n3*x2+2*
         ep^2*n3*x2*x4+10*ep^2*n2+7*ep^2*n2*x4+ep^2*n2*x4^2+4*ep^2*n2*x2+2*ep^2
         *n2*x2*x4+16*ep^2*n1+9*ep^2*n1*x4+ep^2*n1*x4^2+9*ep^2*n1*x2+3*ep^2*n1*
         x2*x4+ep^2*n1*x2^2+12*ep^3+10*ep^3*x4+2*ep^3*x4^2+10*ep^3*x2+7*ep^3*x2
         *x4+ep^3*x2*x4^2+2*ep^3*x2^2+ep^3*x2^2*x4,-12+7*n5-n5^2+7*n2-2*n2*n5-
         n2^2+19*n1-9*n1*n5+n1*n5^2-9*n1*n2+2*n1*n2*n5+n1*n2^2-8*n1^2+2*n1^2*n5
         +2*n1^2*n2+n1^3+10*ep+7*ep*x2-3*ep*n5-2*ep*n5*x2-3*ep*n2-2*ep*n2*x2-13
         *ep*n1-9*ep*n1*x2+3*ep*n1*n5+2*ep*n1*n5*x2+3*ep*n1*n2+2*ep*n1*n2*x2+3*
         ep*n1^2+2*ep*n1^2*x2-2*ep^2-3*ep^2*x2-ep^2*x2^2+2*ep^2*n1+3*ep^2*n1*x2
         +ep^2*n1*x2^2))
      +`Z'(0,x2,0,x4,0,n1,n2,n3,n4,-1+n5)*(rat(12-7*n5+n5^2-7*n3+2*n3*n5+n3^2-7
         *n2+2*n2*n5+2*n2*n3+n2^2-10*ep-7*ep*x2+3*ep*n5+2*ep*n5*x2+3*ep*n3+2*ep
         *n3*x2+3*ep*n2+2*ep*n2*x2+2*ep^2+3*ep^2*x2+ep^2*x2^2,12-7*n5+n5^2-7*n2
         +2*n2*n5+n2^2-7*n1+2*n1*n5+2*n1*n2+n1^2-10*ep-7*ep*x2+3*ep*n5+2*ep*n5*
         x2+3*ep*n2+2*ep*n2*x2+3*ep*n1+2*ep*n1*x2+2*ep^2+3*ep^2*x2+ep^2*x2^2));
*--#] t1,noints=2,4,signs=+,o,+,o,+,order=1,2,4,5,3,level=1,stop=2 : 
* LHS = Z(1,n2,n3,n4,n5)
* DEN = den(-1+n3)*den(-2+n5+n2+ep)*den(-3+n5+n4+n3+ep)
* MAX = Z(0,0,0,0,0)
* MIN = Z(-1,0,-1,0,-1)
* nterms = 3 (Z: 2, Y: 1)
* NonIntegerIndices: 2,4
* Signs: 1,o,+,o,+
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,x2?!{0,},0,x4?!{0,},0,1,n2?,n3?{>=2},n4?,n5?pos_) =
*      #include t1.rules # t1,noints=2,4,signs=1,o,+,o,+,level=1,stop=1
*
*--#[ t1,noints=2,4,signs=1,o,+,o,+,level=1,stop=1 :
      +`Z'(0,x2,0,x4,0,1,n2,-1+n3,n4,n5)*(rat(24-14*n5+2*n5^2-4*n4+n4*n5-10*n3+
         3*n3*n5+n3*n4+n3^2-10*n2+3*n2*n5+n2*n4+2*n2*n3+n2^2-24*ep-4*ep*x4-10*
         ep*x2+7*ep*n5+ep*n5*x4+3*ep*n5*x2+2*ep*n4+ep*n4*x2+5*ep*n3+ep*n3*x4+2*
         ep*n3*x2+5*ep*n2+ep*n2*x4+2*ep*n2*x2+6*ep^2+2*ep^2*x4+5*ep^2*x2+ep^2*
         x2*x4+ep^2*x2^2,3-n5-n4-4*n3+n3*n5+n3*n4+n3^2-ep-ep*x4+ep*n3+ep*n3*x4)
         )
      +`Z'(0,x2,0,x4,0,1,n2,n3,n4,-1+n5)*(rat(12-7*n5+n5^2-4*n4+n4*n5-3*n3+n3*
         n5+n3*n4-3*n2+n2*n5+n2*n4-14*ep-4*ep*x4-3*ep*x2+4*ep*n5+ep*n5*x4+ep*n5
         *x2+2*ep*n4+ep*n4*x2+2*ep*n3+ep*n3*x4+2*ep*n2+ep*n2*x4+4*ep^2+2*ep^2*
         x4+2*ep^2*x2+ep^2*x2*x4,6-5*n5+n5^2-2*n4+n4*n5-2*n3+n3*n5-3*n2+n2*n5+
         n2*n4+n2*n3-5*ep-2*ep*x4-3*ep*x2+2*ep*n5+ep*n5*x4+ep*n5*x2+ep*n4+ep*n4
         *x2+ep*n3+ep*n3*x2+ep*n2+ep*n2*x4+ep^2+ep^2*x4+ep^2*x2+ep^2*x2*x4))
      +`Y'(0,x2,0,x4,0,0,n2,n3,n4,n5)*(rat(-18+12*n5-2*n5^2+9*n4-3*n4*n5-n4^2+3
         *n3-n3*n5-n3*n4+3*n2-n2*n5-n2*n4+21*ep+9*ep*x4+3*ep*x2-7*ep*n5-3*ep*n5
         *x4-ep*n5*x2-5*ep*n4-2*ep*n4*x4-ep*n4*x2-2*ep*n3-ep*n3*x4-2*ep*n2-ep*
         n2*x4-6*ep^2-5*ep^2*x4-ep^2*x4^2-2*ep^2*x2-ep^2*x2*x4,2-n5-2*n3+n3*n5-
         n2+n2*n3-ep-ep*x2+ep*n3+ep*n3*x2));
*--#] t1,noints=2,4,signs=1,o,+,o,+,level=1,stop=1 : 
* LHS = Z(1,n2,1,n4,n5)
* DEN = den(-1+n5)*den(-2+n5+n2+ep)*den(-3+n5+n2+2*ep)*den(-4+2*n5+n4+2*ep)
* MAX = Z(0,0,0,1,0)
* MIN = Z(-1,0,-1,0,-1)
* nterms = 4 (Z: 2, Y: 2)
* NonIntegerIndices: 2,4
* Signs: 1,o,1,x,+
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,x2?!{0,},0,x4?!{0,},0,1,n2?,1,n4?,n5?{>=2}) =
*      #include t1.rules # t1,noints=2,4,signs=1,o,1,x,+,level=1,stop=1
*
*--#[ t1,noints=2,4,signs=1,o,1,x,+,level=1,stop=1 :
      +`Z'(0,x2,0,x4,0,1,n2,1,1+n4,-1+n5)*(rat(5*n4-2*n4*n5-n4^2-n2*n4+5*ep*x4-
         2*ep*n5*x4-3*ep*n4-2*ep*n4*x4-ep*n4*x2-ep*n2*x4-3*ep^2*x4-ep^2*x4^2-
         ep^2*x2*x4,2-3*n5+n5^2-n2+n2*n5-ep-ep*x2+ep*n5+ep*n5*x2))
      +`Z'(0,x2,0,x4,0,1,n2,1,n4,-1+n5)*(rat(12-7*n5+n5^2-7*n4+2*n4*n5+n4^2-3*
         n2+n2*n5+n2*n4-17*ep-7*ep*x4-3*ep*x2+5*ep*n5+2*ep*n5*x4+ep*n5*x2+5*ep*
         n4+2*ep*n4*x4+ep*n4*x2+2*ep*n2+ep*n2*x4+6*ep^2+5*ep^2*x4+ep^2*x4^2+2*
         ep^2*x2+ep^2*x2*x4,2-3*n5+n5^2-n2+n2*n5-ep-ep*x2+ep*n5+ep*n5*x2))
      +`Y'(0,x2,0,x4,0,0,n2,0,1+n4,n5)*(rat(30*n4-22*n4*n5+4*n4*n5^2-11*n4^2+4*
         n4^2*n5+n4^3-11*n2*n4+4*n2*n4*n5+2*n2*n4^2+n2^2*n4+30*ep*x4-22*ep*n5*
         x4+4*ep*n5^2*x4-38*ep*n4-22*ep*n4*x4-11*ep*n4*x2+14*ep*n4*n5+8*ep*n4*
         n5*x4+4*ep*n4*n5*x2+7*ep*n4^2+3*ep*n4^2*x4+2*ep*n4^2*x2-11*ep*n2*x4+4*
         ep*n2*n5*x4+7*ep*n2*n4+4*ep*n2*n4*x4+2*ep*n2*n4*x2+ep*n2^2*x4-38*ep^2*
         x4-11*ep^2*x4^2-11*ep^2*x2*x4+14*ep^2*n5*x4+4*ep^2*n5*x4^2+4*ep^2*n5*
         x2*x4+12*ep^2*n4+14*ep^2*n4*x4+3*ep^2*n4*x4^2+7*ep^2*n4*x2+4*ep^2*n4*
         x2*x4+ep^2*n4*x2^2+7*ep^2*n2*x4+2*ep^2*n2*x4^2+2*ep^2*n2*x2*x4+12*ep^3
         *x4+7*ep^3*x4^2+ep^3*x4^3+7*ep^3*x2*x4+2*ep^3*x2*x4^2+ep^3*x2^2*x4,-24
         +32*n5-14*n5^2+2*n5^3+6*n4-5*n4*n5+n4*n5^2+20*n2-18*n2*n5+4*n2*n5^2-5*
         n2*n4+2*n2*n4*n5-4*n2^2+2*n2^2*n5+n2^2*n4+40*ep+6*ep*x4+20*ep*x2-36*ep
         *n5-5*ep*n5*x4-18*ep*n5*x2+8*ep*n5^2+ep*n5^2*x4+4*ep*n5^2*x2-7*ep*n4-5
         *ep*n4*x2+3*ep*n4*n5+2*ep*n4*n5*x2-22*ep*n2-5*ep*n2*x4-8*ep*n2*x2+10*
         ep*n2*n5+2*ep*n2*n5*x4+4*ep*n2*n5*x2+3*ep*n2*n4+2*ep*n2*n4*x2+2*ep*
         n2^2+ep*n2^2*x4-22*ep^2-7*ep^2*x4-22*ep^2*x2-5*ep^2*x2*x4-4*ep^2*x2^2+
         10*ep^2*n5+3*ep^2*n5*x4+10*ep^2*n5*x2+2*ep^2*n5*x2*x4+2*ep^2*n5*x2^2+2
         *ep^2*n4+3*ep^2*n4*x2+ep^2*n4*x2^2+6*ep^2*n2+3*ep^2*n2*x4+4*ep^2*n2*x2
         +2*ep^2*n2*x2*x4+4*ep^3+2*ep^3*x4+6*ep^3*x2+3*ep^3*x2*x4+2*ep^3*x2^2+
         ep^3*x2^2*x4))
      +`Y'(0,x2,0,x4,0,0,n2,1,1+n4,-1+n5)*(rat(-30*n4+22*n4*n5-4*n4*n5^2+16*
         n4^2-6*n4^2*n5-2*n4^3+21*n2*n4-13*n2*n4*n5+2*n2*n4*n5^2-10*n2*n4^2+3*
         n2*n4^2*n5+n2*n4^3-3*n2^2*n4+n2^2*n4*n5+n2^2*n4^2-30*ep*x4+22*ep*n5*x4
         -4*ep*n5^2*x4+68*ep*n4+32*ep*n4*x4+21*ep*n4*x2-36*ep*n4*n5-12*ep*n4*n5
         *x4-13*ep*n4*n5*x2+4*ep*n4*n5^2+2*ep*n4*n5^2*x2-26*ep*n4^2-6*ep*n4^2*
         x4-10*ep*n4^2*x2+6*ep*n4^2*n5+3*ep*n4^2*n5*x2+2*ep*n4^3+ep*n4^3*x2+21*
         ep*n2*x4-13*ep*n2*n5*x4+2*ep*n2*n5^2*x4-29*ep*n2*n4-20*ep*n2*n4*x4-6*
         ep*n2*n4*x2+9*ep*n2*n4*n5+6*ep*n2*n4*n5*x4+2*ep*n2*n4*n5*x2+7*ep*n2*
         n4^2+3*ep*n2*n4^2*x4+2*ep*n2*n4^2*x2-3*ep*n2^2*x4+ep*n2^2*n5*x4+2*ep*
         n2^2*n4+2*ep*n2^2*n4*x4+68*ep^2*x4+16*ep^2*x4^2+21*ep^2*x2*x4-36*ep^2*
         n5*x4-6*ep^2*n5*x4^2-13*ep^2*n5*x2*x4+4*ep^2*n5^2*x4+2*ep^2*n5^2*x2*x4
         -50*ep^2*n4-52*ep^2*n4*x4-6*ep^2*n4*x4^2-29*ep^2*n4*x2-20*ep^2*n4*x2*
         x4-3*ep^2*n4*x2^2+14*ep^2*n4*n5+12*ep^2*n4*n5*x4+9*ep^2*n4*n5*x2+6*
         ep^2*n4*n5*x2*x4+ep^2*n4*n5*x2^2+10*ep^2*n4^2+6*ep^2*n4^2*x4+7*ep^2*
         n4^2*x2+3*ep^2*n4^2*x2*x4+ep^2*n4^2*x2^2-29*ep^2*n2*x4-10*ep^2*n2*x4^2
         -6*ep^2*n2*x2*x4+9*ep^2*n2*n5*x4+3*ep^2*n2*n5*x4^2+2*ep^2*n2*n5*x2*x4+
         10*ep^2*n2*n4+14*ep^2*n2*n4*x4+3*ep^2*n2*n4*x4^2+4*ep^2*n2*n4*x2+4*
         ep^2*n2*n4*x2*x4+2*ep^2*n2^2*x4+ep^2*n2^2*x4^2-50*ep^3*x4-26*ep^3*x4^2
         -2*ep^3*x4^3-29*ep^3*x2*x4-10*ep^3*x2*x4^2-3*ep^3*x2^2*x4+14*ep^3*n5*
         x4+6*ep^3*n5*x4^2+9*ep^3*n5*x2*x4+3*ep^3*n5*x2*x4^2+ep^3*n5*x2^2*x4+12
         *ep^3*n4+20*ep^3*n4*x4+6*ep^3*n4*x4^2+10*ep^3*n4*x2+14*ep^3*n4*x2*x4+3
         *ep^3*n4*x2*x4^2+2*ep^3*n4*x2^2+2*ep^3*n4*x2^2*x4+10*ep^3*n2*x4+7*ep^3
         *n2*x4^2+ep^3*n2*x4^3+4*ep^3*n2*x2*x4+2*ep^3*n2*x2*x4^2+12*ep^4*x4+10*
         ep^4*x4^2+2*ep^4*x4^3+10*ep^4*x2*x4+7*ep^4*x2*x4^2+ep^4*x2*x4^3+2*ep^4
         *x2^2*x4+ep^4*x2^2*x4^2,24-56*n5+46*n5^2-16*n5^3+2*n5^4-6*n4+11*n4*n5-
         6*n4*n5^2+n4*n5^3-20*n2+38*n2*n5-22*n2*n5^2+4*n2*n5^3+5*n2*n4-7*n2*n4*
         n5+2*n2*n4*n5^2+4*n2^2-6*n2^2*n5+2*n2^2*n5^2-n2^2*n4+n2^2*n4*n5-40*ep-
         6*ep*x4-20*ep*x2+76*ep*n5+11*ep*n5*x4+38*ep*n5*x2-44*ep*n5^2-6*ep*n5^2
         *x4-22*ep*n5^2*x2+8*ep*n5^3+ep*n5^3*x4+4*ep*n5^3*x2+7*ep*n4+5*ep*n4*x2
         -10*ep*n4*n5-7*ep*n4*n5*x2+3*ep*n4*n5^2+2*ep*n4*n5^2*x2+22*ep*n2+5*ep*
         n2*x4+8*ep*n2*x2-32*ep*n2*n5-7*ep*n2*n5*x4-12*ep*n2*n5*x2+10*ep*n2*
         n5^2+2*ep*n2*n5^2*x4+4*ep*n2*n5^2*x2-3*ep*n2*n4-2*ep*n2*n4*x2+3*ep*n2*
         n4*n5+2*ep*n2*n4*n5*x2-2*ep*n2^2-ep*n2^2*x4+2*ep*n2^2*n5+ep*n2^2*n5*x4
         +22*ep^2+7*ep^2*x4+22*ep^2*x2+5*ep^2*x2*x4+4*ep^2*x2^2-32*ep^2*n5-10*
         ep^2*n5*x4-32*ep^2*n5*x2-7*ep^2*n5*x2*x4-6*ep^2*n5*x2^2+10*ep^2*n5^2+3
         *ep^2*n5^2*x4+10*ep^2*n5^2*x2+2*ep^2*n5^2*x2*x4+2*ep^2*n5^2*x2^2-2*
         ep^2*n4-3*ep^2*n4*x2-ep^2*n4*x2^2+2*ep^2*n4*n5+3*ep^2*n4*n5*x2+ep^2*n4
         *n5*x2^2-6*ep^2*n2-3*ep^2*n2*x4-4*ep^2*n2*x2-2*ep^2*n2*x2*x4+6*ep^2*n2
         *n5+3*ep^2*n2*n5*x4+4*ep^2*n2*n5*x2+2*ep^2*n2*n5*x2*x4-4*ep^3-2*ep^3*
         x4-6*ep^3*x2-3*ep^3*x2*x4-2*ep^3*x2^2-ep^3*x2^2*x4+4*ep^3*n5+2*ep^3*n5
         *x4+6*ep^3*n5*x2+3*ep^3*n5*x2*x4+2*ep^3*n5*x2^2+ep^3*n5*x2^2*x4));
*--#] t1,noints=2,4,signs=1,o,1,x,+,level=1,stop=1 : 
* LHS = Z(1,n2,1,n4,1)
* DEN = den(-1+n2+ep)*den(-2+n2+2*ep)
* MAX = Z(0,0,0,1,0)
* MIN = Z(-1,-1,-1,0,-1)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 2,4
* Signs: 1,+*,1,o,1
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,x2?!{0,},0,x4?!{0,},0,1,n2?{>=2},1,n4?,1) =
*      #include t1.rules # t1,noints=2,4,signs=1,+*,1,o,1,level=1,stop=1
*
*--#[ t1,noints=2,4,signs=1,+*,1,o,1,level=1,stop=1 :
      +`Z'(0,x2,0,x4,0,1,-1+n2,1,n4,1)*(rat(-3+n4+n2+3*ep+ep*x4+ep*x2,-1+n2+ep+
         ep*x2))
      +`Y'(0,x2,0,x4,0,0,n2,0,1+n4,1)*(rat(-3*n4+n4^2+n2*n4-3*ep*x4+3*ep*n4+2*
         ep*n4*x4+ep*n4*x2+ep*n2*x4+3*ep^2*x4+ep^2*x4^2+ep^2*x2*x4,2-3*n2+n2^2-
         4*ep-3*ep*x2+3*ep*n2+2*ep*n2*x2+2*ep^2+3*ep^2*x2+ep^2*x2^2))
      +`Y'(0,x2,0,x4,0,1,n2,1,n4,0)*(rat(2-n4-2*ep-ep*x4,-1+n2+ep+ep*x2));
*--#] t1,noints=2,4,signs=1,+*,1,o,1,level=1,stop=1 : 
* LHS = Z(1,n2,1,n4,1)
* DEN = den(-1+n2+2*ep)*den(-2+n4+n2+3*ep)
* MAX = Z(0,1,0,1,0)
* MIN = Z(-1,1,-1,0,-1)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 2,4
* Signs: 1,-,1,o,1
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 2
* Stop: 3
*
    id `IDSPEC' `Z'(0,x2?!{0,},0,x4?!{0,},0,1,n2?neg0_,1,n4?,1) =
*      #include t1.rules # t1,noints=2,4,signs=1,-,1,o,1,level=2,stop=3
*
*--#[ t1,noints=2,4,signs=1,-,1,o,1,level=2,stop=3 :
      +`Z'(0,x2,0,x4,0,1,1+n2,1,n4,1)*(rat(n2+ep+ep*x2,-2+n4+n2+3*ep+ep*x4+ep*
         x2))
      +`Y'(0,x2,0,x4,0,0,1+n2,0,1+n4,1)*(rat(-n4-ep*x4,-1+n2+2*ep+ep*x2))
      +`Y'(0,x2,0,x4,0,1,1+n2,1,n4,0)*(rat(-2+n4+2*ep+ep*x4,-2+n4+n2+3*ep+ep*x4
         +ep*x2));
*--#] t1,noints=2,4,signs=1,-,1,o,1,level=2,stop=3 : 
* LHS = Z(1,n2,1,n4,1)
* DEN = den(-1+n2+ep)*den(-1+n4+ep)*den(-2+n4+2*ep)
* MAX = Z(0,0,0,1,0)
* MIN = Z(-1,0,-1,-1,-1)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 2,4
* Signs: 1,o,1,+,1
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,x2?!{0,},0,x4?!{0,},0,1,n2?,1,n4?{>=2},1) =
*      #include t1.rules # t1,noints=2,4,signs=1,o,1,+,1,level=1,stop=1
*
*--#[ t1,noints=2,4,signs=1,o,1,+,1,level=1,stop=1 :
      +`Z'(0,x2,0,x4,0,1,n2,1,-1+n4,1)*(rat(-3+n4+n2+3*ep+ep*x4+ep*x2,-1+n4+ep+
         ep*x4))
      +`Y'(0,x2,0,x4,0,0,n2,0,1+n4,1)*(rat(-3*n4+n4^2+n2*n4-3*ep*x4+3*ep*n4+2*
         ep*n4*x4+ep*n4*x2+ep*n2*x4+3*ep^2*x4+ep^2*x4^2+ep^2*x2*x4,2-n4-2*n2+n2
         *n4-4*ep-ep*x4-2*ep*x2+ep*n4+ep*n4*x2+2*ep*n2+ep*n2*x4+2*ep^2+ep^2*x4+
         2*ep^2*x2+ep^2*x2*x4))
      +`Y'(0,x2,0,x4,0,1,n2,1,n4,0)*(rat(2-n2-2*ep-ep*x2,-1+n4+ep+ep*x4));
*--#] t1,noints=2,4,signs=1,o,1,+,1,level=1,stop=1 : 
* LHS = Z(1,n2,1,n4,1)
* DEN = den(-1+n4+2*ep)*den(-2+n4+n2+3*ep)
* MAX = Z(0,1,0,1,0)
* MIN = Z(-1,0,-1,1,-1)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 2,4
* Signs: 1,o,1,-,1
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 2
* Stop: 3
*
    id `IDSPEC' `Z'(0,x2?!{0,},0,x4?!{0,},0,1,n2?,1,n4?neg0_,1) =
*      #include t1.rules # t1,noints=2,4,signs=1,o,1,-,1,level=2,stop=3
*
*--#[ t1,noints=2,4,signs=1,o,1,-,1,level=2,stop=3 :
      +`Z'(0,x2,0,x4,0,1,n2,1,1+n4,1)*(rat(n4+ep+ep*x4,-2+n4+n2+3*ep+ep*x4+ep*
         x2))
      +`Y'(0,x2,0,x4,0,0,1+n2,0,1+n4,1)*(rat(-n2-ep*x2,-1+n4+2*ep+ep*x4))
      +`Y'(0,x2,0,x4,0,1,n2,1,1+n4,0)*(rat(-2+n2+2*ep+ep*x2,-2+n4+n2+3*ep+ep*x4
         +ep*x2));
*--#] t1,noints=2,4,signs=1,o,1,-,1,level=2,stop=3 : 
*--#] t1,24-rules
        goto 2;
label 1;
        #call ReduceZ2Y(Z,Y,5,ifmatch->2,1,3,5)
        $doloop = 1;
label 2;
        ModuleOption maximum, $doloop;
        .sort:`TOPO'-`nloop++';
        #if `$doloop'
          #redefine loop "0"
        #endif
      #enddo
      #break
*--#] t1,24
*--#[ t1,34
    #case t1,34
      #define nloop "1"
      #do loop=1,1
        #$doloop = 0;
        if (count(`Z',1) == 0) goto 2;
*--#[ t1,34-rules
* LHS = Z(n1,n2,n3,n4,n5)
* DEN = den(-3+n5+n2+n1+ep)*den(-4+n5+n2+n1+2*ep)
* MAX = Z(0,0,0,0,0)
* MIN = Z(-1,-1,0,0,-1)
* nterms = 3 (Z: 3, Y: 0)
* NonIntegerIndices: 3,4
* Signs: +,+,o,o,+
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,x3?!{0,},x4?!{0,},0,n1?pos_,n2?pos_,n3?,n4?,n5?pos_) =
*      #include t1.rules # t1,noints=3,4,signs=+,+,o,o,+,level=1,stop=0
*
*--#[ t1,noints=3,4,signs=+,+,o,o,+,level=1,stop=0 :
      +`Z'(0,0,x3,x4,0,-1+n1,n2,n3,n4,n5)*(rat(28-15*n5+2*n5^2-11*n4+3*n4*n5+
         n4^2-4*n3+n3*n5+n3*n4-4*n2+n2*n5+n2*n4-11*n1+3*n1*n5+2*n1*n4+n1*n3+n1*
         n2+n1^2-26*ep-11*ep*x4-4*ep*x3+7*ep*n5+3*ep*n5*x4+ep*n5*x3+5*ep*n4+2*
         ep*n4*x4+ep*n4*x3+2*ep*n3+ep*n3*x4+2*ep*n2+ep*n2*x4+5*ep*n1+2*ep*n1*x4
         +ep*n1*x3+6*ep^2+5*ep^2*x4+ep^2*x4^2+2*ep^2*x3+ep^2*x3*x4,12-7*n5+n5^2
         -7*n2+2*n2*n5+n2^2-7*n1+2*n1*n5+2*n1*n2+n1^2-10*ep+3*ep*n5+3*ep*n2+3*
         ep*n1+2*ep^2))
      +`Z'(0,0,x3,x4,0,n1,-1+n2,n3,n4,n5)*(rat(28-15*n5+2*n5^2-4*n4+n4*n5-11*n3
         +3*n3*n5+n3*n4+n3^2-11*n2+3*n2*n5+n2*n4+2*n2*n3+n2^2-4*n1+n1*n5+n1*n3+
         n1*n2-26*ep-4*ep*x4-11*ep*x3+7*ep*n5+ep*n5*x4+3*ep*n5*x3+2*ep*n4+ep*n4
         *x3+5*ep*n3+ep*n3*x4+2*ep*n3*x3+5*ep*n2+ep*n2*x4+2*ep*n2*x3+2*ep*n1+ep
         *n1*x3+6*ep^2+2*ep^2*x4+5*ep^2*x3+ep^2*x3*x4+ep^2*x3^2,12-7*n5+n5^2-7*
         n2+2*n2*n5+n2^2-7*n1+2*n1*n5+2*n1*n2+n1^2-10*ep+3*ep*n5+3*ep*n2+3*ep*
         n1+2*ep^2))
      +`Z'(0,0,x3,x4,0,n1,n2,n3,n4,-1+n5)*(rat(-16+8*n5-n5^2+4*n4-n4*n5+4*n3-n3
         *n5-n3*n4+4*n2-n2*n5-n2*n4+4*n1-n1*n5-n1*n3-n1*n2+16*ep+4*ep*x4+4*ep*
         x3-4*ep*n5-ep*n5*x4-ep*n5*x3-2*ep*n4-ep*n4*x3-2*ep*n3-ep*n3*x4-2*ep*n2
         -ep*n2*x4-2*ep*n1-ep*n1*x3-4*ep^2-2*ep^2*x4-2*ep^2*x3-ep^2*x3*x4,12-7*
         n5+n5^2-7*n2+2*n2*n5+n2^2-7*n1+2*n1*n5+2*n1*n2+n1^2-10*ep+3*ep*n5+3*ep
         *n2+3*ep*n1+2*ep^2));
*--#] t1,noints=3,4,signs=+,+,o,o,+,level=1,stop=0 : 
*--#] t1,34-rules
        goto 2;
label 1;
        #call ReduceZ2Y(Z,Y,5,ifmatch->2,1,2,5)
        $doloop = 1;
label 2;
        ModuleOption maximum, $doloop;
        .sort:`TOPO'-`nloop++';
        #if `$doloop'
          #redefine loop "0"
        #endif
      #enddo
      #break
*--#] t1,34
*--#[ t1,45
    #case t1,45
      #define nloop "1"
      #do loop=1,1
        #$doloop = 0;
        if (count(`Z',1) == 0) goto 2;
*--#[ t1,45-rules
* LHS = Z(n1,n2,n3,n4,n5)
* DEN = den(-1+n3)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,-1,0,0)
* nterms = 4 (Z: 4, Y: 0)
* NonIntegerIndices: 4,5
* Signs: +,+,+,o,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 0
*
    id `IDSPEC' `Z'(0,0,0,x4?!{0,},x5?!{0,},n1?pos_,n2?pos_,n3?{>=2},n4?,n5?) =
*      #include t1.rules # t1,noints=4,5,signs=+,+,+,o,o,level=1,stop=0
*
*--#[ t1,noints=4,5,signs=+,+,+,o,o,level=1,stop=0 :
      +`Z'(0,0,0,x4,x5,-1+n1,n2,-1+n3,n4,1+n5)*(rat(-n5-ep*x5,-1+n3))
      +`Z'(0,0,0,x4,x5,n1,-1+n2,-1+n3,n4,1+n5)*(rat(n5+ep*x5,-1+n3))
      +`Z'(0,0,0,x4,x5,n1,-1+n2,n3,n4,n5)*(rat(1,1))
      +`Z'(0,0,0,x4,x5,n1,n2,-1+n3,n4,n5)*(rat(-5+n5+n3+2*n2+2*ep+ep*x5,-1+n3))
         ;
*--#] t1,noints=4,5,signs=+,+,+,o,o,level=1,stop=0 : 
* LHS = Z(n1,n2,1,n4,n5)
* DEN = den(-1+n1)*den(-4+n5+2*n2+2*ep)*den(-3+n5+n2+n1+ep)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,-1,0,0)
* nterms = 4 (Z: 2, Y: 2)
* NonIntegerIndices: 4,5
* Signs: +,+,1,o,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,x4?!{0,},x5?!{0,},n1?{>=2},n2?pos_,1,n4?,n5?) =
*      #include t1.rules # t1,noints=4,5,signs=+,+,1,o,o,level=1,stop=1
*
*--#[ t1,noints=4,5,signs=+,+,1,o,o,level=1,stop=1 :
      +`Z'(0,0,0,x4,x5,-1+n1,n2,1,n4,n5)*(rat(12-7*n5+n5^2-7*n4+2*n4*n5+n4^2-7*
         n1+2*n1*n5+2*n1*n4+n1^2-10*ep-7*ep*x5-7*ep*x4+3*ep*n5+2*ep*n5*x5+2*ep*
         n5*x4+3*ep*n4+2*ep*n4*x5+2*ep*n4*x4+3*ep*n1+2*ep*n1*x5+2*ep*n1*x4+2*
         ep^2+3*ep^2*x5+ep^2*x5^2+3*ep^2*x4+2*ep^2*x4*x5+ep^2*x4^2,3-n5-n2-4*n1
         +n1*n5+n1*n2+n1^2-ep-ep*x5+ep*n1+ep*n1*x5))
      +`Z'(0,0,0,x4,x5,n1,-1+n2,1,n4,n5)*(rat(-3+n5+n2+2*ep+ep*x5,-3+n5+n2+n1+
         ep+ep*x5))
      +`Y'(0,0,0,x4,x5,-1+n1,n2,0,n4,1+n5)*(rat(3*n5-n5^2-n2*n5+3*ep*x5-2*ep*n5
         -2*ep*n5*x5-ep*n2*x5-2*ep^2*x5-ep^2*x5^2,4-n5-2*n2-4*n1+n1*n5+2*n1*n2-
         2*ep-ep*x5+2*ep*n1+ep*n1*x5))
      +`Y'(0,0,0,x4,x5,n1,-1+n2,0,n4,1+n5)*(rat(-3*n5+n5^2+n2*n5-3*ep*x5+2*ep*
         n5+2*ep*n5*x5+ep*n2*x5+2*ep^2*x5+ep^2*x5^2,4-n5-2*n2-4*n1+n1*n5+2*n1*
         n2-2*ep-ep*x5+2*ep*n1+ep*n1*x5));
*--#] t1,noints=4,5,signs=+,+,1,o,o,level=1,stop=1 : 
* LHS = Z(1,n2,1,n4,n5)
* DEN = den(-1+n2)*den(-2+n5+n2+ep)*den(-3+n5+n2+ep)*den(-3+n5+n2+2*ep)
* MAX = Z(0,0,0,0,0)
* MIN = Z(-1,-2,-1,0,0)
* nterms = 3 (Z: 2, Y: 1)
* NonIntegerIndices: 4,5
* Signs: 1,+*,1,o,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 2
* Stop: 2
*
    id `IDSPEC' `Z'(0,0,0,x4?!{0,},x5?!{0,},1,n2?{>=2},1,n4?,n5?) =
*      #include t1.rules # t1,noints=4,5,signs=1,+*,1,o,o,level=2,stop=2
*
*--#[ t1,noints=4,5,signs=1,+*,1,o,o,level=2,stop=2 :
      +`Z'(0,0,0,x4,x5,1,-2+n2,1,n4,n5)*(rat(-20+9*n5-n5^2+4*n4-n4*n5+9*n2-2*n2
         *n5-n2*n4-n2^2+22*ep+9*ep*x5+4*ep*x4-5*ep*n5-2*ep*n5*x5-ep*n5*x4-2*ep*
         n4-ep*n4*x5-5*ep*n2-2*ep*n2*x5-ep*n2*x4-6*ep^2-5*ep^2*x5-ep^2*x5^2-2*
         ep^2*x4-ep^2*x4*x5,2-n5-3*n2+n2*n5+n2^2-ep-ep*x5+ep*n2+ep*n2*x5))
      +`Z'(0,0,0,x4,x5,1,-1+n2,1,n4,n5)*(rat(16-8*n5+n5^2-2*n4+n4*n5-12*n2+3*n2
         *n5+n2*n4+2*n2^2-13*ep-8*ep*x5-2*ep*x4+3*ep*n5+2*ep*n5*x5+ep*n5*x4+ep*
         n4*x5+6*ep*n2+3*ep*n2*x5+ep*n2*x4+2*ep^2+3*ep^2*x5+ep^2*x5^2+ep^2*x4*
         x5,2-n5-3*n2+n2*n5+n2^2-ep-ep*x5+ep*n2+ep*n2*x5))
      +`Y'(0,0,0,x4,x5,0,n2,0,n4,n5)*(rat(-90+63*n5-14*n5^2+n5^3+33*n4-17*n4*n5
         +2*n4*n5^2-3*n4^2+n4^2*n5+78*n2-36*n2*n5+4*n2*n5^2-20*n2*n4+5*n2*n4*n5
         +n2*n4^2-22*n2^2+5*n2^2*n5+3*n2^2*n4+2*n2^3+139*ep+63*ep*x5+33*ep*x4-
         64*ep*n5-28*ep*n5*x5-17*ep*n5*x4+7*ep*n5^2+3*ep*n5^2*x5+2*ep*n5^2*x4-
         31*ep*n4-17*ep*n4*x5-6*ep*n4*x4+8*ep*n4*n5+4*ep*n4*n5*x5+2*ep*n4*n5*x4
         +ep*n4^2+ep*n4^2*x5-83*ep*n2-36*ep*n2*x5-20*ep*n2*x4+19*ep*n2*n5+8*ep*
         n2*n5*x5+5*ep*n2*n5*x4+10*ep*n2*n4+5*ep*n2*n4*x5+2*ep*n2*n4*x4+12*ep*
         n2^2+5*ep*n2^2*x5+3*ep*n2^2*x4-71*ep^2-64*ep^2*x5-14*ep^2*x5^2-31*ep^2
         *x4-17*ep^2*x4*x5-3*ep^2*x4^2+16*ep^2*n5+14*ep^2*n5*x5+3*ep^2*n5*x5^2+
         8*ep^2*n5*x4+4*ep^2*n5*x4*x5+ep^2*n5*x4^2+7*ep^2*n4+8*ep^2*n4*x5+2*
         ep^2*n4*x5^2+2*ep^2*n4*x4+2*ep^2*n4*x4*x5+22*ep^2*n2+19*ep^2*n2*x5+4*
         ep^2*n2*x5^2+10*ep^2*n2*x4+5*ep^2*n2*x4*x5+ep^2*n2*x4^2+12*ep^3+16*
         ep^3*x5+7*ep^3*x5^2+ep^3*x5^3+7*ep^3*x4+8*ep^3*x4*x5+2*ep^3*x4*x5^2+
         ep^3*x4^2+ep^3*x4^2*x5,-18+21*n5-8*n5^2+n5^3+21*n2-16*n2*n5+3*n2*n5^2-
         8*n2^2+3*n2^2*n5+n2^3+27*ep+21*ep*x5-21*ep*n5-16*ep*n5*x5+4*ep*n5^2+3*
         ep*n5^2*x5-21*ep*n2-16*ep*n2*x5+8*ep*n2*n5+6*ep*n2*n5*x5+4*ep*n2^2+3*
         ep*n2^2*x5-13*ep^2-21*ep^2*x5-8*ep^2*x5^2+5*ep^2*n5+8*ep^2*n5*x5+3*
         ep^2*n5*x5^2+5*ep^2*n2+8*ep^2*n2*x5+3*ep^2*n2*x5^2+2*ep^3+5*ep^3*x5+4*
         ep^3*x5^2+ep^3*x5^3));
*--#] t1,noints=4,5,signs=1,+*,1,o,o,level=2,stop=2 : 
* LHS = Z(1,1,1,n4,n5)
* DEN = den(-1+n4)*den(-2+n5+n4+ep)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,-1,-1,0)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 4,5
* Signs: 1,1,1,+,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,x4?!{0,},x5?!{0,},1,1,1,n4?{>=2},n5?) =
*      #include t1.rules # t1,noints=4,5,signs=1,1,1,+,o,level=1,stop=1
*
*--#[ t1,noints=4,5,signs=1,1,1,+,o,level=1,stop=1 :
      +`Z'(0,0,0,x4,x5,1,1,1,-1+n4,n5)*(rat(-3+n5+n4+2*ep+ep*x5+ep*x4,-1+n4+ep*
         x4))
      +`Y'(0,0,0,x4,x5,0,1,0,n4,1+n5)*(rat(n5+ep*x5,-2+n5+n4+ep+ep*x5+ep*x4))
      +`Y'(0,0,0,x4,x5,1,0,1,-1+n4,1+n5)*(rat(-n5-ep*x5,-1+n4+ep*x4));
*--#] t1,noints=4,5,signs=1,1,1,+,o,level=1,stop=1 : 
* LHS = Z(1,1,1,n4,n5)
* DEN = den(-1+n5+n4+ep)*den(-2+n5+n4+2*ep)
* MAX = Z(0,0,0,1,1)
* MIN = Z(-1,-1,-1,1,0)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 4,5
* Signs: 1,1,1,-,o
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 2
* Stop: 2
*
    id `IDSPEC' `Z'(0,0,0,x4?!{0,},x5?!{0,},1,1,1,n4?neg0_,n5?) =
*      #include t1.rules # t1,noints=4,5,signs=1,1,1,-,o,level=2,stop=2
*
*--#[ t1,noints=4,5,signs=1,1,1,-,o,level=2,stop=2 :
      +`Z'(0,0,0,x4,x5,1,1,1,1+n4,n5)*(rat(n4+ep*x4,-2+n5+n4+2*ep+ep*x5+ep*x4))
      +`Y'(0,0,0,x4,x5,0,1,0,1+n4,1+n5)*(rat(-n4*n5-ep*n5*x4-ep*n4*x5-ep^2*x4*
         x5,2-3*n5+n5^2-3*n4+2*n4*n5+n4^2-4*ep-3*ep*x5-3*ep*x4+3*ep*n5+2*ep*n5*
         x5+2*ep*n5*x4+3*ep*n4+2*ep*n4*x5+2*ep*n4*x4+2*ep^2+3*ep^2*x5+ep^2*x5^2
         +3*ep^2*x4+2*ep^2*x4*x5+ep^2*x4^2))
      +`Y'(0,0,0,x4,x5,1,0,1,1+n4,n5)*(rat(-2+n5+2*ep+ep*x5,-2+n5+n4+2*ep+ep*x5
         +ep*x4));
*--#] t1,noints=4,5,signs=1,1,1,-,o,level=2,stop=2 : 
* LHS = Z(1,1,1,n4,n5)
* DEN = den(-1+n5+ep)*den(-2+n5+2*ep)*den(-2+n5+n4+ep)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,-1,0,-1)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 4,5
* Signs: 1,1,1,o,+
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 1
* Stop: 1
*
    id `IDSPEC' `Z'(0,0,0,x4?!{0,},x5?!{0,},1,1,1,n4?,n5?{>=2}) =
*      #include t1.rules # t1,noints=4,5,signs=1,1,1,o,+,level=1,stop=1
*
*--#[ t1,noints=4,5,signs=1,1,1,o,+,level=1,stop=1 :
      +`Z'(0,0,0,x4,x5,1,1,1,n4,-1+n5)*(rat(3-n5-n4-2*ep-ep*x5-ep*x4,-1+n5+ep+
         ep*x5))
      +`Y'(0,0,0,x4,x5,0,1,0,n4,1+n5)*(rat(-4*n5+2*n5^2+n4*n5-4*ep*x5+3*ep*n5+4
         *ep*n5*x5+ep*n5*x4+ep*n4*x5+3*ep^2*x5+2*ep^2*x5^2+ep^2*x4*x5,4-4*n5+
         n5^2-2*n4+n4*n5-6*ep-4*ep*x5-2*ep*x4+3*ep*n5+2*ep*n5*x5+ep*n5*x4+2*ep*
         n4+ep*n4*x5+2*ep^2+3*ep^2*x5+ep^2*x5^2+2*ep^2*x4+ep^2*x4*x5))
      +`Y'(0,0,0,x4,x5,1,0,1,n4,n5)*(rat(-4+2*n5+n4+3*ep+2*ep*x5+ep*x4,-1+n5+ep
         +ep*x5));
*--#] t1,noints=4,5,signs=1,1,1,o,+,level=1,stop=1 : 
* LHS = Z(1,1,1,n4,n5)
* DEN = den(-1+n5+2*ep)*den(-2+n5+n4+2*ep)
* MAX = Z(0,0,0,0,1)
* MIN = Z(-1,-1,0,0,1)
* nterms = 3 (Z: 1, Y: 2)
* NonIntegerIndices: 4,5
* Signs: 1,1,1,o,-
* Order: 1,2,3,4,5
* Seed: 0,0,0,0,0
* Level: 2
* Stop: 2
*
    id `IDSPEC' `Z'(0,0,0,x4?!{0,},x5?!{0,},1,1,1,n4?,n5?neg0_) =
*      #include t1.rules # t1,noints=4,5,signs=1,1,1,o,-,level=2,stop=2
*
*--#[ t1,noints=4,5,signs=1,1,1,o,-,level=2,stop=2 :
      +`Z'(0,0,0,x4,x5,1,1,1,n4,1+n5)*(rat(-n5-ep-ep*x5,-2+n5+n4+2*ep+ep*x5+ep*
         x4))
      +`Y'(0,0,0,x4,x5,0,1,1,n4,1+n5)*(rat(-2+2*n5+n4+3*ep+2*ep*x5+ep*x4,-1+n5+
         2*ep+ep*x5))
      +`Y'(0,0,0,x4,x5,1,0,1,n4,1+n5)*(rat(-2+2*n5+n4+3*ep+2*ep*x5+ep*x4,-2+n5+
         n4+2*ep+ep*x5+ep*x4));
*--#] t1,noints=4,5,signs=1,1,1,o,-,level=2,stop=2 : 
*--#] t1,45-rules
        goto 2;
label 1;
        #call ReduceZ2Y(Z,Y,5,ifmatch->2,1,2,3)
        $doloop = 1;
label 2;
        ModuleOption maximum, $doloop;
        .sort:`TOPO'-`nloop++';
        #if `$doloop'
          #redefine loop "0"
        #endif
      #enddo
      #break
*--#] t1,45
  #endswitch
#endprocedure

*--#] ReduceImpl : 
*--#[ ReduceZ2Y :

**
* Changes Z to Y if specified index is non-positive. The replacement is
* performed like
*   id Z(x1?,0,x3?,n1?,n2?neg0_,n3?) = Y(x1,0,x3,n1,n2,n3);
*
#procedure ReduceZ2Y(Z,Y,N,IDSPEC,?ii)
  #do i={`?ii',}
    #ifdef `i'
      #if (`i' >= 1) && (`i' <= `N')
        #if `N' == 1
          id `IDSPEC' `Z'(0,n1?neg0_)
            = `Y'(0,n1);
        #elseif `i' == 1
          id `IDSPEC' `Z'(0,x2?,...,x`N'?,n1?neg0_,n2?,...,n`N'?)
            = `Y'(0,x2,...,x`N',n1,...,n`N');
        #elseif `i' == `N'
          id `IDSPEC' `Z'(x1?,...,x{`N'-1}?,0,n1?,...,n{`N'-1}?,n`N'?neg0_)
            = `Y'(x1,...,x{`N'-1},0,n1,...,n`N');
        #else
          id `IDSPEC' `Z'(x1?,...,x{`i'-1}?,0,x{`i'+1}?,...,x`N'?,n1?,...,n{`i'-1}?,n`i'?neg0_,n{`i'+1}?,...,n`N'?)
            = `Y'(x1,...,x{`i'-1},0,x{`i'+1},...,x`N',n1,...,n`N');
        #endif
      #else
        #message Error: illegal arguments : RudeceZ2Y(`Z',`Y',`N',`IDSPEC',`?ii')
        #terminate
      #endif
    #endif
  #enddo
#endprocedure

*--#] ReduceZ2Y : 
*--#[ ReduceImpl2 :

**
* New reduction schemes. Integrals are represented as Z(n1,n2,n3,n4,n5) and
* the non-integer parts of indices are understood by `TOPO'.
*
* NOTE: Because of the current usage of spectators, the active expression must
*       have a name `NAME'`TOPO'.
*
#procedure ReduceImpl2(TOPO,Z,N)
  #define Y "Y"
  #ifdef `USEZZ'
    #define ZZ "ZZ"
  #endif
  #define SPECTATOR ""
#ifndef `NOSPECTATORS'
  CreateSpectator x`TOPO' "x`TOPO'-`PID_'.spec";
  #redefine SPECTATOR "x`TOPO'"
#endif
  #call ExtractQ2(`Z',`N',`Q2')
  #switch `TOPO'
    #case nostar0
      #call SchemeZno()
      #break
    #default
      #call SchemeZ`TOPO'()
      #break
  #endswitch
#ifndef `NOSPECTATORS'
  .sort:`TOPO'-pre_from_spectator;
  #call FromSpectator(x`TOPO',`NAME'`TOPO'res,`NAME'`TOPO'res,`Y')
  .sort:`TOPO'-from_spectator;
  Drop `NAME'`TOPO',`NAME'`TOPO'res;
  L `NAME'tmp = `NAME'`TOPO' + `NAME'`TOPO'res;
  #exchange `NAME'`TOPO',`NAME'tmp
  ModuleOption noparallel;
  .sort:`TOPO'-merge_spectator;
#endif
  multiply replace_(`Y',`Z');
  #call UnextractQ2(`Z',`N',`Q2')
#endprocedure

*--#] ReduceImpl2 : 
*--#[ StdReduceProc :

**
* The standard reduction procedure. It calls ReduceLoop() with procedures named
* as "reduce*", "turbo" with "z2y" and optionally "sym".
*
* Parameters:
*   arg1,... - reduce     : Use reduce and z2y.
*              turbo      : Use turbo and z2y.
*              sym        : Use sym.
*              <otherwise>: Use reduce`arg1' and z2y.
*
* Variables:
*   TOPO - The name of the topology (used in commentary output).
*
#procedure StdReduceProc(?a)
  #procedure DoReduceLoop(i,sym)
    #ifdef `i'
      #ifdef `sym'
        #define sympostfix "-s"
      #else
        #define sympostfix ""
      #endif
      #if "`i'" == "reduce"
        #call ReduceLoop(reduce,z2y,`sym',`TOPO'-reduce`sympostfix')
      #elseif "`i'" == "turbo"
        #define flag "1"
        #ifdef `TURBOTHRESHOLD'
          #if `TURBOTHRESHOLD' == 0
            #redefine flag "0"
          #endif
        #endif
        #if `flag'
          #call ReduceLoop(turbo,z2y,`sym',`TOPO'-turbo`sympostfix')
        #endif
      #else
        #call ReduceLoop(reduce`i',z2y,`sym',`TOPO'-`i'`sympostfix')
      #endif
    #endif
  #endprocedure

  #define prev ""
  #do a={`?a',}
    #ifdef `a'
      #if "`a'" == "sym"
        #ifdef `prev'
          #call DoReduceLoop(`prev',sym)
          #redefine prev ""
        #else
          #call sym()
          .sort:`TOPO'-sym;
        #endif
      #else
        #call DoReduceLoop(`prev',)
        #redefine prev "`a'"
      #endif
    #endif
  #enddo
  #call DoReduceLoop(`prev',)
#endprocedure

*--#] StdReduceProc : 
*--#[ ReduceLoop :

**
* Makes a reduction loop.
*
* Parameters:
*   proc1   - The name of procedure to be called for Z-integrals. Typically it
*             contains rules like:
*
*               `id' `Z'(n1?pos_,...) = `Z'(n1-1,...) + ...;
*
*   proc2   - The name of procedure to be called for Z-integrals for which rules
*             are applied in proc1. Typically it contains rules like:
*
*               `id' `Z'(n1?neg0_,...) = `Y'(n1,...);
*
*   proc3   - The name of procedure to be called for integrals that will not
*             match to any rules. Typically it contains symmetry rules like:
*
*               `id' `Z'(n1?,n2?,...) = `Z'(n1,n2,...)*`Z'(n2,n1,...);
*               repeat `id `Z'(?a) * `Z'(?b) = `Z'(?a);
*
*   text... - Text printed at .sort.
*
* Variables:
*   `Z'         - The function representing Z-integrals.
*   `Y'         - The function representing Y-integrals.
*   `ZZ'        - The function temporarily used for Z-integrals that don't match
*                 to any rules. (optional)
*   `SPECTATOR' - The name of spectator that Y-integrals will be sent into.
*                 (optional)
*
#procedure ReduceLoop(proc1,proc2,proc3,?text)
  #define id
  #define pass "1"
  #do loop=1,1
    if (count(`Z',1) == 0) goto 3;
    #ifdef `proc1'
      #redefine id "id ifmatch->1"
      #call `proc1'
    #endif
    #ifdef `proc3'
      #redefine id "id"
      #call `proc3'
    #endif
    #ifdef `ZZ'
      multiply replace_(`Z',`ZZ');
    #endif
    goto 2;
label 1;
    #ifdef `proc2'
      #redefine id "id ifmatch->2"
      #call `proc2'
    #endif
    redefine loop "0";
    goto 3;
label 2;
    #ifdef `proc3'
      #redefine id "id"
      #define Z "`Y'"
      #call `proc3'
      #undefine Z
    #endif
label 3;
    .sort:`?text',pass=`pass++';
    #ifdef `ZZ'
      #if `loop'
        multiply replace_(`ZZ',`Z');
      #endif
    #endif
#ifdef `SPECTATOR'
    if (count(`Z',1) == 0) ToSpectator `SPECTATOR';
#endif
  #enddo
#endprocedure

*--#] ReduceLoop : 
*--#[ Z2Y :

**
* Changes Z to Y if one of indices for propagators is non-positive.
*
* Parameters:
*   N         - The number of indices both for propagators and dot products.
*   M         - The number of indices only for propagators.
*   star1,... - The indices associated with non-integer powers. (optional)
*
* Variables:
*   `Z'  - The function representing Z-integrals.
*   `Y'  - The function representing Y-integrals.
*   `id' - The "identify" statement.
*
* Prerequisites:
*   S n1,...,n`N';
*
* Examples:
*   * t1
*   #call Z2Y(5,5)
*
*   * no
*   #call Z2Y(9,8)
*
*   * nostar5
*   #call Z2Y(9,8,5)
*
#procedure Z2Y(N,M,?STARS)
  #define flag
  #do i=1,`M'
    #redefine flag "1"
    #do j={`?STARS',}
      #ifdef `j'
        #if `i' == `j'
          #redefine flag "0"
          #breakdo
        #endif
      #endif
    #enddo
    #if `flag'
      #if `N' != `M'
        #if (`i' > 1) && (`i' < `M')
          `id' `Z'(<n1?>,...,<n{`i'-1}?>,n`i'?neg0_,<n{`i'+1}?>,...,<n`M'?>,<n{`M'+1}?neg0_>,...,<n`N'?neg0_>)
        #elseif (`i' == 1) && (`M' == 1)
          `id' `Z'(n`i'?neg0_,<n{`M'+1}?neg0_>,...,<n`N'?neg0_>)
        #elseif `i' == 1
          `id' `Z'(n`i'?neg0_,<n{`i'+1}?>,...,<n`M'?>,<n{`M'+1}?neg0_>,...,<n`N'?neg0_>)
        #else
          `id' `Z'(<n1?>,...,<n{`i'-1}?>,n`i'?neg0_,<n{`M'+1}?neg0_>,...,<n`N'?neg0_>)
        #endif
      #else
        #if (`i' > 1) && (`i' < `M')
          `id' `Z'(<n1?>,...,<n{`i'-1}?>,n`i'?neg0_,<n{`i'+1}?>,...,<n`M'?>)
        #elseif (`i' == 1) && (`M' == 1)
          `id' `Z'(n`i'?neg0_)
        #elseif `i' == 1
          `id' `Z'(n`i'?neg0_,<n{`i'+1}?>,...,<n`M'?>)
        #else
          `id' `Z'(<n1?>,...,<n{`i'-1}?>,n`i'?neg0_)
        #endif
      #endif
      = `Y'(n1,...,n`N');
    #endif
  #enddo
#endprocedure

*--#] Z2Y : 
*--#[ ExtractQ2 :

**
* Moves to the dimensionless notation. The mass dimension is expressed in terms
* of `Q2'. Does nothing when SETQONE=1.
*
#procedure ExtractQ2(Z,N,Q2)
  #define flag "1"
  #ifdef `SETQONE'
    #if `SETQONE'
      #redefine flag "0"
    #endif
  #endif
  #if `flag'
    #ifdef `N'
      id `Z'(n1?,...,n`N'?) = `Z'(n1,...,n`N') * `Q2'^(-n1-...-n`N');
    #else
*     N is not specified. Consider all cases for 4-loop p-integrals.
      #define id "id"
      #do N={14,9,5,2}
        `id' `Z'(n1?,...,n`N'?) = `Z'(n1,...,n`N') * `Q2'^(-n1-...-n`N');
        #redefine id "al"
      #enddo
    #endif
  #endif
#endprocedure

*--#] ExtractQ2 : 
*--#[ UnextractQ2 :

**
* Moves back to the dimensionful notation.
*
#procedure UnextractQ2(Z,N,Q2)
  #define flag "1"
  #ifdef `SETQONE'
    #if `SETQONE'
      #redefine flag "0"
    #endif
  #endif
  #if `flag'
    #ifdef `N'
      id `Z'(n1?,...,n`N'?) = `Z'(n1,...,n`N') * `Q2'^(n1+...+n`N');
    #else
*     N is not specified. Consider all cases for 4-loop p-integrals.
      #define id "id"
      #do N={14,9,5,2}
        `id' `Z'(n1?,...,n`N'?) = `Z'(n1,...,n`N') * `Q2'^(-n1-...-n`N');
        #redefine id "al"
      #enddo
    #endif
  #endif
#endprocedure

*--#] UnextractQ2 : 
*--#[ Timing :

#define FORCERoldtime "0.0"

**
* Prints the timing when FORCERTIMING is defined. If this procedure is called
* without any arguments, it does not print the timing but resets the timer.
*
#procedure Timing(?message)
  #ifdef `FORCERTIMING'
    #ifdef `?message'
      #define t "{(`timer_'+5)/10}"
      #define sec "{`t'/100}"
      #define csec "{`t'%100}"
      #if `csec' < 10
        #redefine csec "0`csec'"
      #endif
      #message Timing: `?message': `sec'.`csec' sec (`FORCERoldtime' -> `time_' sec)
    #endif
    #reset timer
    #redefine FORCERoldtime "`time_'"
  #endif
#endprocedure

*--#] Timing : 
*--#[ CheckTopoHash :

**
* Checks the topology hash. It detects possible version conflicts in future,
* e.g., changes of unnamed topologies.
*
#procedure CheckTopoHash(?hash)
  #if `?hash' != `topohash'
    #message Error: hash value mismatch
    #message Expected:`?hash'
    #message Actual: `topohash'
    #terminate
  #endif
#endprocedure

*--#] CheckTopoHash : 
*--#] proc : 
#endif
* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4
