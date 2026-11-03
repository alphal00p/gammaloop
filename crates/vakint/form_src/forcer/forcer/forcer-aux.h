#ifndef `FORCERAUXHFILE'
#define FORCERAUXHFILE
*--#[ decl :

#ifndef `SETQONE'
  #define SETQONE "1"
#endif

S ep;
S Q2;  * (= Q.Q)
CF rat,RAT,GschemeConstants;

**
* POLYRATFUNEXPAND - Select options for PolyRatFun.
*                    Positive integer : Expand rat-functions up to ep^n.
*                    Zero             : Keep only the most divergent term in ep.
*                    Negative integer : Exact rational arithmetic. (default)
*
* See EnablePolyratfun().
*
*#define POLYRATFUNEXPAND "-1"

**
* ALLTIMING - Record timings into the file `NAME'.time if defined.
* MINTIME   - Do not record the duration for each topology if < MINTIME
*             in milliseconds. (default: 0)
*
*#define ALLTIMING ""
*#define MINTIME "0"

**
* HASHBASE - An integer used for hashing isomorphic graphs.
*
* See hash().
*
#define HASHBASE "20"

**
* Private variables.
*
#define PREFIX "FORCER"
#define NMAX "29"

S `PREFIX'n1,...,`PREFIX'n`NMAX';
S `PREFIX'x1,...,`PREFIX'x`NMAX';
V `PREFIX'q1,...,`PREFIX'q`NMAX';
S `PREFIX'isum1,...,`PREFIX'isum4;
S `PREFIX'xrhs;
CF `PREFIX'G,`PREFIX'fnum,`PREFIX'del,`PREFIX'distrib;
CF `PREFIX'poch,`PREFIX'pochinv,`PREFIX'num,`PREFIX'den;
CT `PREFIX'ftensor;

NTable,sparse,OneFill,`PREFIX'dotmaptabl(2);

#define `PREFIX'LABELNO "0"
#define `PREFIX'MAXINDEX "14"

#undefine PREFIX
#undefine NMAX

*--#] decl : 
*--#[ proc :
*--#[ EnablePolyratfun :

**
* Enables PolyRatFun with options specified by `POLYRATFUNEXPAND'.
*
* NOTE: forcer-manual.h#FromSpectator contains the same code.
*
#procedure EnablePolyratfun()
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
#endprocedure

*--#] EnablePolyratfun : 
*--#[ ReduceEdge :

**
* Reduces the given topologies by calling Reduce`toponame' for each topology.
*
* Parameters:
*   nedge         - the number of edges (used in commentary output)
*   toponame1,... - the names of topologies to be processed
*
* Remarks:
*   BeforeReduce(TOPO) and AfterReduce(TOPO) are user customization points.
*   They are always called even if the corresponding module is not generated,
*   i.e., `HASTOPO`TOPO'' is 0. They need to check `HASTOPO`TOPO'' when they
*   contain any FORM statements.
*
#procedure ReduceEdge(nedge,?topolist)
  #define prefix "FORCER"
  #define n      "`prefix'n"
  #define oldn   "`prefix'oldn"
  #define MAX    "``prefix'MAXINDEX'"

* Check if we have any topologies to be processed at this level.
  #define HASTOPO "0"
  #do toponame={`?topolist'}
    #if `HASTOPO`toponame''
      #redefine HASTOPO "1"
      #breakdo
    #endif
  #enddo

  #if `HASTOPO'
    #call EnablePolyratfun()
    IntoHide `NAME';

    B+
      #do toponame={`?topolist'}
        M`toponame'
      #enddo
    ;

#ifdef `FORCERDEBUG'
    Print +f +s;
#endif

    .sort:starting with `nedge' edges;

#ifdef `ALLTIMING'
    #$`n' = termsin_(`NAME');
    #write <`NAME'.time> "\n starting with `nedge' edges after `time_' sec. (`$`n'' terms)"
    #reset timer

    #ifndef `MINTIME'
      #define MINTIME "0"
    #endif
#endif

    #define processedtopolist ""

    #do toponame={`?topolist'}
      #if `HASTOPO`toponame''
        L `NAME'`toponame' = `NAME'[M`toponame'] * M`toponame';
#ifdef `FORCERSTATS'
*       NOTE: Computing the complexity of integrals based on the difference of
*       indices to Z(1,...,1,0,...,0) requires information of each topology (the
*       numbers of propagators and irreducible numerators). Instead, we just
*       compute the sum of the absolute values of indices.
        #$`prefix'minindex = 9999;
        #$`prefix'maxindex = -9999;
        #$`prefix'minindexabssum = 9999;
        #$`prefix'maxindexabssum = -9999;
        if (match(M`toponame'/<i1^`n'1?$`n'1>/.../<i`MAX'^`n'`MAX'?$`n'`MAX'>));
          $`prefix'minindex = min_($`prefix'minindex,<$`n'1>,...,<$`n'`MAX'>);
          $`prefix'maxindex = max_($`prefix'maxindex,<$`n'1>,...,<$`n'`MAX'>);
          $`n' = <abs_($`n'1)>+...+<abs_($`n'`MAX')>;
          $`prefix'minindexabssum = min_($`prefix'minindexabssum,$`n');
          $`prefix'maxindexabssum = max_($`prefix'maxindexabssum,$`n');
        endif;
        ModuleOption local $`n',<$`n'1>,...,<$`n'`MAX'>;
        ModuleOption minimum $`prefix'minindex, $`prefix'minindexabssum;
        ModuleOption maximum $`prefix'maxindex, $`prefix'maxindexabssum;
#endif
        ModuleOption noparallel;
        .sort:begin `toponame';
#ifdef `ALLTIMING'
        #$`oldn' = termsin_(`NAME'`toponame');
#endif
#ifdef `FORCERSTATS'
        #if `$`prefix'minindex' != 9999
          #message Stats: Complexity(`toponame'): index=(`$`prefix'minindex',`$`prefix'maxindex'), indexabssum=(`$`prefix'minindexabssum',`$`prefix'maxindexabssum')
        #endif
#endif
      #endif
      #call BeforeReduce(`toponame')
      #if `HASTOPO`toponame''
        #call Reduce`toponame'()
      #endif
      #call AfterReduce(`toponame')
      #if `HASTOPO`toponame''
        Hide `NAME'`toponame';
        #redefine processedtopolist "`processedtopolist',`toponame'"
#ifdef `ALLTIMING'
      #if `timer_' >= `MINTIME'
        #$`n' = termsin_(`NAME'`toponame');
        #write <`NAME'.time> "\t {`timer_'/1000} sec. spent reducing `toponame' (`$`oldn'' -> `$`n'' terms)"
      #endif
      #reset timer
#endif
      #endif
    #enddo

    .sort:ending with `nedge' edges;

    PolyRatFun;
    Drop
      `NAME',
      #do toponame={`processedtopolist'}
        #ifdef `toponame'
          `NAME'`toponame',
        #endif
      #enddo
    ;
    L `NAME'tmp = `NAME'[1]
      #do toponame={`processedtopolist'}
        #ifdef `toponame'
          + `NAME'`toponame'
        #endif
      #enddo
    ;
    #exchange `NAME',`NAME'tmp
    AB+ rat;
    ModuleOption noparallel;
    .sort:term merge with `nedge' edges;
  #else
*   HASTOPO = 0. Still call user customization points.
    #do toponame={`?topolist'}
      #call BeforeReduce(`toponame')
      #call AfterReduce(`toponame')
    #enddo
  #endif
#endprocedure

*--#] ReduceEdge : 
*--#[ replace :

#procedure replace(lhs,rhs,j,?vars,M,msg)
* Replace dot products and filter tadpoles
  #define PREFIX "FORCER"
  #$maxlhs = 0;
  if (count(`M',1) && ($maxlhs < count(`lhs',1)))
    $maxlhs = count_(`lhs',1);
  AB+ `lhs',`M',`?vars';
  ModuleOption maximum,$maxlhs;
  .sort:`msg'-prepare;

  #$maxlhs = {`$maxlhs'/4}*4;
  #do i=`$maxlhs'-3,0,-4
    if (count(`M',1));
      id `lhs'^n?{>`i'} = `lhs'^`i' * `PREFIX'xrhs^n/`PREFIX'xrhs^`i';
      id	`PREFIX'xrhs = `rhs';
    endif;
    #call tadpole(`M',`j')
    AB+ `lhs',`M',`?vars';
    .sort:`msg'-`i';
  #enddo
  if (count(`M',1))
    id `lhs' = `rhs';
  #call tadpole(`M',`j')
#endprocedure

*--#] replace : 
*--#[ hash :

#procedure hash(iso,?args)
* hashes an isomorphism. Arguments should be given from lowest to highest priority.
  #define HPOW "0"
  $iso`iso' =
  #do i={`?args'}
    + abs_(`i')*`HASHBASE'^`HPOW++'
  #enddo
  ;
#endprocedure

*--#] hash : 
*--#[ CheckRemovedPropagators :

**
* Checks if propagators are removed, i.e., indices are zero or negative.
* When the i-th line is missing, M`TOPO' is replaced with M`TOPO'M`i', or
* dropped if it gives a tadpole (specified by the optional arguments).
*
* Parameters:
*   TOPO                    - the name of the topology to be checked
*   numv                    - the number of propagators (i1,...,i`numv')
*   nstars                  - the number of indices with non-integer powers
*   star1,...,nstar`nstars' - the indices associated with non-integer powers.
*                             It skips checks for those indices
*   tad1,...                - the indices for propagators that give a tadpole
*                             when they are removed (optional)
*
* Results:
*   $DOACTION       - set to 1 when there are M`TOPO' without missing lines,
*                     which should be reduced in the next action. Otherwise 0
*   $M1,...#M`numv' - set to 1 when there are M`TOPO' with the i-th line
*                     missing, i.e., M`TOPO'M`i' in the resultant expressions.
*                     Otherwise 0
*
* Examples:
*   * t1: ni <= 0 (i = 1,...,5) gives a simpler topology.
*   #call CheckRemovedPropagators(t1,5,0)
*
*   * t1star45: n1 <= 0 or n2 <= 0 or n3 <= 0 gives a simpler topology.
*   #call CheckRemovedPropagators(t1star45,5,2,4,5)
*
*   * l1 (one-loop): n1 <= 0 or n2 <= 0 gives a tadpole.
*   #call CheckRemovedPropagators(l1,2,0,1,2)
*
*   * l1star1: n2 <= 0 gives a tadpole.
*   #call CheckRemovedPropagators(l1star1,2,1,1,2)
*
*   * l1star12: missing lines gives no simplifications.
*   #call CheckRemovedPropagators(l1star12,2,2,1,2)
*
#procedure CheckRemovedPropagators(TOPO,numv,numskip,?specials)
  #$DOACTION = 0;

  #do i=1,`numv'
    #define LINEACTION`i' "0"
    #$M`i' = 0;
  #enddo

  #define SPECCOUNT "0"
  #do i={`?specials',}
    #ifdef `i'
      #if `SPECCOUNT++' < `numskip'
        #redefine LINEACTION`i' "-1"
      #else
        #redefine LINEACTION`i' "1"
      #endif
    #endif
  #enddo

  #do p=1,`numv'
    if (count(M`TOPO'M`p',1)) $M`p'=1;
  #enddo

  if (count(M`TOPO',1));
    #define ELSE ""
    #do p=1,`numv'
      #if `LINEACTION`p'' >= 0
        `ELSE'if (count(i`p',1) >= 0);
          #if `LINEACTION`p''
            Discard;
          #else
            Multiply replace_(M`TOPO',M`TOPO'M`p');
            $M`p'=1;
          #endif
        #redefine ELSE "else"
      #endif
    #enddo
    #if `ELSE' == else
      else;
        $DOACTION = 1;
      endif;
    #else
      $DOACTION = 1;
    #endif
  endif;
  ModuleOption maximum,<$M1>,...,<$M`numv'>,$DOACTION;
  .sort:relabel-`TOPO';
#endprocedure

*--#] CheckRemovedPropagators : 
*--#[ MtoZ :

**
* Converts M/i1^n1/i2^n2/... into Z(n1,n2,...).
*
* Parameters:
*   M - the symbol for the topology
*   Z - the function for the topology
*   i - the prefix for the symbols for the indices
*   N - the number of the indices
*
#procedure MtoZ(M,Z,i,N)
  #define prefix "FORCER"
  #define n      "`prefix'n"

#ifdef `FORCERSTATS'
  #$`prefix'nterms = 0;
  #$`prefix'nspecs = 0;
  if (count(`M',1));
    $`prefix'nterms = $`prefix'nterms + 1;
  else;
    $`prefix'nspecs = $`prefix'nspecs + 1;
  endif;
#endif
  id `M'*<`i'1^`n'1?>*...*<`i'`N'^`n'`N'?> = `Z'(<-`n'1>,...,<-`n'`N'>);
#ifdef `FORCERSTATS'
  ModuleOption sum, $`prefix'nterms, $`prefix'nspecs;
#endif
  .sort:MtoZ-`M'-`Z';
#ifdef `FORCERSTATS'
  #message Stats: MtoZ(`M'->`Z'): nterms=`$`prefix'nterms', nspecs=`$`prefix'nspecs'
#endif
#endprocedure

*--#] MtoZ : 
*--#[ ZtoM :

**
* Converts Z(n1,n2,...) into M/i1^n1/i2^n2/...
*
* Parameters:
*   Z - the function for the topology
*   M - the symbol for the topology
*   i - the prefix for the symbols for the indices
*   N - the number of the indices
*
#procedure ZtoM(Z,M,i,N)
  #define prefix "FORCER"
  #define n      "`prefix'n"

  id `Z'(<`n'1?>,...,<`n'`N'?>) = `M'*<1/`i'1^`n'1>*...*<1/`i'`N'^`n'`N'>;
#ifdef `FORCERSTATS'
  #$`prefix'nterms = 0;
  #$`prefix'nspecs = 0;
  if (count(`M',1));
    $`prefix'nterms = $`prefix'nterms + 1;
  else;
    $`prefix'nspecs = $`prefix'nspecs + 1;
  endif;
  ModuleOption sum, $`prefix'nterms, $`prefix'nspecs;
#endif
  .sort:ZtoM-`Z'-`M';
#ifdef `FORCERSTATS'
  #message Stats: ZtoM(`Z'->`M'): nterms=`$`prefix'nterms', nspecs=`$`prefix'nspecs'
#endif
#endprocedure

*--#] ZtoM : 
*--#[ triangle :

#define TRISUM

**
* Applies the triangle rule recursively. The assignment of indices is as
* follows:
*
*         a4         a3         a5
*       ----------------------------
*               \          /
*                \        /
*        a1+e1*ep \      / a2+e2*ep
*                  \    /
*                   \  /
*                    \/
*                    /\
*
* with possible dot products in the numerator (p3.Q1)^(-aK1), ...
* So the integrand looks like
*
*   1/p1.p1^(a1+e1*ep)
*    /p2.p2^(a2+e2*ep)
*    /p3.p3^a3
*    /p4.p4^a4
*    /p5.p5^a5
*    /p3.Q1^aK1 ...
*
* Note that the non-integer parts e3,...,e5 are assumed to be zero, because
* otherwise the triangle rule can not be applied. (Actually e4 and e5 are not
* needed due to no their appearance in the right-hand side.) eK's are also
* assumed to be zero. It may be possible to introduce eK's as well, but
* in practice we do not need non-integer parts for dot products if they are in
* the numerator. The caller must be sure that the integral passed to this
* procedure does not have non-integer parts e3,...,e5,eK1,... = 0.
*
* Parameters:
*   Z           - the name of the integral function
*   i1,...,i5   - the positions of the indices for the edges depicted above
*   iK1,iK2,... - the positions of the index for the numerator (optional)
*   e1,e2       - represent the non-integer parts of the two lines depicted
*                 above
*
* i1,...,i5,iK1,... are specified by one-based numbering. e1 and e2 are expected
* be non-negative integers. No extra signs ("+") or no extra parentheses ("("
* and ")") are allowed.
*
* This procedure will change the indices specified by i1,...,i5. a1 and a2 may
* be increased while a3, a4 and a5 are decreased until one of them becomes zero.
*
* The current version uses the explicit summed formula by default. The recursive
* version can be used by
*   #undefine TRISUM
*
#procedure triangle(Z,i1,i2,i3,i4,i5,?iKj,e1,e2)
  #call Timing()

  #define PREFIX "FORCER"

  #define n "`PREFIX'n"

* Determine the maximum index.
  #define N "`i1'"
  #do i={`i2',`i3',`i4',`i5',`?iKj'}
    #ifdef `i'
      #if `i' > `N'
        #redefine N "`i'"
      #endif
    #endif
  #enddo

* TODO: dynamically switch from the recursive version to the summed version?

#ifdef `TRISUM'

  #define k1 "`PREFIX'isum1"
  #define k2 "`PREFIX'isum2"
  #define k3 "`PREFIX'isum3"
  #define k4 "`PREFIX'isum4"

  #define n1 "`n'`i1'"
  #define n2 "`n'`i2'"
  #define n3 "`n'`i3'"
  #define n4 "`n'`i4'"
  #define n5 "`n'`i5'"
  #define a1 "(`n1'+`e1'*ep)"
  #define a2 "(`n2'+`e2'*ep)"
  #define a3 "`n3'"
  #define a4 "`n4'"
  #define a5 "`n5'"

  #define poch "`PREFIX'poch"
  #define num "`PREFIX'num"
  #define den "`PREFIX'den"

  if (count(`Z',1));

    #do i=1,`N'
      #define nL`i'  "`n'`i'?"
      #define nR`i'a "`n'`i'"
      #define nR`i'b "`n'`i'"
      #define nR`i'c "`n'`i'"
    #enddo

    #define M "0"
    #do i={`?iKj',}
      #ifdef `i'
        #redefine M "`M'-`n'`i'"
      #endif
    #enddo
    #redefine M "(`M')"

    #redefine nL`i3' "`n3'?pos_"
    #redefine nL`i4' "`n4'?pos_"
    #redefine nL`i5' "`n5'?pos_"
    #redefine nR`i1'a "`n1'+`k1'+`k2'"
    #redefine nR`i2'a "`n2'+`n3'-`k1'+`k4'"
    #redefine nR`i3'a "0"
    #redefine nR`i4'a "`n4'-`k2'"
    #redefine nR`i5'a "`n5'-`k4'"
    #redefine nR`i1'b "`n1'+`n4'+`k1'"
    #redefine nR`i2'b "`n2'+`k3'+`k4'"
    #redefine nR`i3'b "`n3'-`k1'-`k3'"
    #redefine nR`i4'b "0"
    #redefine nR`i5'b "`n5'-`k4'"
    #redefine nR`i1'c "`n1'+`k1'+`k2'"
    #redefine nR`i2'c "`n2'+`n5'+`k3'"
    #redefine nR`i3'c "`n3'-`k1'-`k3'"
    #redefine nR`i4'c "`n4'-`k2'"
    #redefine nR`i5'c "0"

    id `Z'(
        #do i=1,`N'
          `nL`i'',
        #enddo
      ?nn) = sum_(`k1',0,`a3',sum_(`k2',0,`a4'-1,sum_(`k4',0,`a5'-1,`Z'(
        #do i=1,`N'
          `nR`i'a',
        #enddo
      ?nn) * sign_(`k2'+`k4')
           * `a3'
           * fac_(`a3'+`k2'+`k4'-1)
           * invfac_(`k1')
           * invfac_(`k2')
           * invfac_(`k4')
           * invfac_(`a3'-`k1')
           * `poch'(-`a3'-`k2'-`k4',4-2*ep+`M'-`a1'-`a2'-`a3')
           * `poch'(`k1'+`k2',`a1')
           * `poch'(`a3'-`k1'+`k4',`a2')
      ))) + sum_(`k1',0,`a3'-1,sum_(`k3',0,`a3'-1-`k1',sum_(`k4',0,`a5'-1,`Z'(
        #do i=1,`N'
          `nR`i'b',
        #enddo
      ?nn) * sign_(`a4'+`k4')
           * fac_(`a4'+`k1'+`k3'+`k4'-1)
           * invfac_(`k1')
           * invfac_(`k3')
           * invfac_(`k4')
           * invfac_(`a4'-1)
           * `poch'(-`a4'-`k1'-`k3'-`k4',4-2*ep+`M'-`a1'-`a2'-2*`a3'+`k1'+`k3'+1)
           * `poch'(`a4'+`k1',`a1')
           * `poch'(`k3'+`k4',`a2')
      ))) + sum_(`k1',0,`a3'-1,sum_(`k2',0,`a4'-1,sum_(`k3',0,`a3'-1-`k1',`Z'(
        #do i=1,`N'
          `nR`i'c',
        #enddo
      ?nn) * sign_(`a5'+`k2')
           * fac_(`a5'+`k1'+`k2'+`k3'-1)
           * invfac_(`k1')
           * invfac_(`k2')
           * invfac_(`k3')
           * invfac_(`a5'-1)
           * `poch'(-`a5'-`k1'-`k2'-`k3',4-2*ep+`M'-`a1'-`a2'-2*`a3'+`k1'+`k3'+1)
           * `poch'(`k1'+`k2',`a1')
           * `poch'(`a5'+`k3',`a2')
      )));

    #define n "`PREFIX'n1"
    #define x "`PREFIX'x1"

    id `poch'(0,`x'?) = 1;
    id `poch'(`n'?pos_,`x'?pos_) = fac_(`n'-1+`x') * invfac_(`x'-1);
  endif;

  B `poch';
  .sort:triangle-`Z'-1;
  Keep Brackets;

  repeat id `poch'(`n'?pos_,`x'?) = `poch'(`n'-1,`x')*`num'(`n'-1+`x');
  repeat id `poch'(`n'?neg_,`x'?) = `poch'(`n'+1,`x')*`den'(`n'+`x');

  id `poch'(0,`x'?) = 1;
* id `num'(`x'?)*`den'(`x'?) = 1;
* id `num'(`x'?number_) = `x';
* id `den'(`x'?number_) = 1/`x';
  id `num'(`x'?) = rat(`x',1);
  id `den'(`x'?) = rat(1,`x');

  .sort:triangle-`Z'-2;
#else

  #do i=1,`N'
    #define nL`i'  "`n'`i'?"
    #define nR`i'a "`n'`i'"
    #define nR`i'b "`n'`i'"
    #define nR`i'c "`n'`i'"
    #define nR`i'd "`n'`i'"
  #enddo

  #define Rshift ""
  #do i={`?iKj',}
    #ifdef `i'
      #define Rshift "`Rshift'-`n'`i'"
    #endif
  #enddo

  #redefine nL`i3' "`n'`i3'?pos_"
  #redefine nL`i4' "`n'`i4'?pos_"
  #redefine nL`i5' "`n'`i5'?pos_"
  #redefine nR`i1'a "`n'`i1'+1"
  #redefine nR`i3'a "`n'`i3'-1"
  #redefine nR`i1'b "`n'`i1'+1"
  #redefine nR`i4'b "`n'`i4'-1"
  #redefine nR`i2'c "`n'`i2'+1"
  #redefine nR`i3'c "`n'`i3'-1"
  #redefine nR`i2'd "`n'`i2'+1"
  #redefine nR`i5'd "`n'`i5'-1"

  #define den "4-`n'`i1'-`n'`i2'-2*`n'`i3'`Rshift'-{2+`e1'+`e2'}*ep"

  if (count(`Z',1));
    repeat id `Z'(
        #do i=1,`N'
          `nL`i'',
        #enddo
      ?nn) = `Z'(
        #do i=1,`N'
          `nR`i'a',
        #enddo
      ?nn) * rat( `n'`i1'+`e1'*ep,`den') + `Z'(
        #do i=1,`N'
          `nR`i'b',
        #enddo
      ?nn) * rat(-`n'`i1'-`e1'*ep,`den') + `Z'(
        #do i=1,`N'
          `nR`i'c',
        #enddo
      ?nn) * rat( `n'`i2'+`e2'*ep,`den') + `Z'(
        #do i=1,`N'
          `nR`i'd',
        #enddo
      ?nn) * rat(-`n'`i2'-`e2'*ep,`den');
  endif;

  .sort:triangle-`Z';
#endif

  #call Timing(Triangle(`Z'))
#endprocedure

*--#] triangle : 
*--#[ fastt1 :

**
* Applies the fast triangle rule of the T1 topology recursively. The assignment
* of the indices is as follows.
*
*                     a1+e1*ep         a2
*                     --<-------------<--
*                    /         |         \
*                   /          |          \
*                  /           |           \
*                 /            |            \
*          Q --<--             ^ a5          --<-- Q
*                 \            |            /
*                  \           |           /
*                   \          |          /
*                    \         |         /
*                     -->------------->--
*                     a4+e4*ep         a3
*
* Parameters:
*   Z           - the name of the integral function
*   i1,...,i5   - the positions of the indices for the edges depicted above
*   e1,e4       - represent the non-integer parts of the two lines depicted
*                 above
*
* i1,...,i5 are specified by one-based numbering. e1 and e4 are expected
* be non-negative integers. No extra signs ("+") or no extra parentheses ("("
* and ")") are allowed.
*
* This procedure will change the indices specified by i2, i3 and i5; a3, a4 and
* a5 will be decreased until one of them becomes zero.
*
#procedure fastt1(Z,i1,i2,i3,i4,i5,e1,e4)
  #call Timing()

  #define PREFIX "FORCER"

  #define n "`PREFIX'n"

  #define x "`PREFIX'x1"

  #define num "`PREFIX'num"
  #define den "`PREFIX'den"

* Determine the maximum index.
  #define N "`i1'"
  #do i={`i2',`i3',`i4',`i5'}
    #ifdef `i'
      #if `i' > `N'
        #redefine N "`i'"
      #endif
    #endif
  #enddo

  #do i=1,`N'
    #define nL`i'  "`n'`i'?"
    #define nR`i'a "`n'`i'"
    #define nR`i'b "`n'`i'"
    #define nR`i'c "`n'`i'"
  #enddo

  #redefine nL`i2' "`n'`i2'?pos_"
  #redefine nL`i3' "`n'`i3'?pos_"
  #redefine nL`i5' "`n'`i5'?pos_"
  #redefine nR`i2'a "`n'`i2'-1"
  #redefine nR`i3'b "`n'`i3'-1"
  #redefine nR`i5'c "`n'`i5'-1"

  #define n1 "`n'`i1'"
  #define n2 "`n'`i2'"
  #define n3 "`n'`i3'"
  #define n4 "`n'`i4'"
  #define n5 "`n'`i5'"
  #define a1 "(`n1'+`e1'*ep)"
  #define a2 "`n2'"
  #define a3 "`n3'"
  #define a4 "(`n4'+`e4'*ep)"
  #define a5 "`n5'"

  #define A  "(14-6*ep-2*`a1'-2*`a2'-2*`a3'-2*`a4'-4*`a5')"
  #define B  "(6-2*ep-2*`a1'-2*`a4'-2*`a5')"
  #define C  "(4-2*ep-`a1'-`a4'-`a5')"
  #define D2 "(6-2*ep-2*`a1'-2*`a2'-2*`a5')"
  #define D3 "(6-2*ep-2*`a3'-2*`a4'-2*`a5')"

  #if `SETQONE'
    #define Q2 "1"
  #else
    #define Q2 "Q2"
  #endif

  #define nloop "1"
  B `Z';
  .sort:fastt1-`Z'-`nloop++';
  #do loop=1,1
    Keep Brackets;
    #$doloop = 0;
    id ifnomatch->1, `Z'(
        #do i=1,`N'
          `nL`i'',
        #enddo
      ?nn) = `Z'(
        #do i=1,`N'
          `nR`i'a',
        #enddo
      ?nn) * `num'(`A') * `den'(`D2') / `Q2' + `Z'(
        #do i=1,`N'
          `nR`i'b',
        #enddo
      ?nn) * `num'(`A') * `den'(`D3') / `Q2' + `Z'(
        #do i=1,`N'
          `nR`i'c',
        #enddo
      ?nn) * (-2) * `num'(`B') * `num'(`C') * `den'(`D2') * `den'(`D3') / `Q2';
    $doloop = 1;
    id `num'(`x'?)*`den'(`x'?) = 1;
*   id `num'(`x'?number_) = `x';
*   id `den'(`x'?number_) = 1/`x';
    id `num'(`x'?) = rat(`x',1);
    id `den'(`x'?) = rat(1,`x');
label 1;
    B `Z';
    ModuleOption maximum, $doloop;
    .sort:fastt1-`Z'-`nloop++';
    #if `$doloop'
      #redefine loop "0"
    #endif
  #enddo

  #call Timing(FastT1(`Z'))
#endprocedure

*--#] fastt1 : 
*--#[ insertion :

**
* Performs massless one-loop insertion integrals:
*
*                     p1
*               ------<------
*              /   a1+e1*ep  \
*             /               \
*            /                 \
*     Q --<--                   --<-- Q
*            \                 /
*             \               /
*              \   a2+e2*ep  /
*               ------<------
*                     p2
*
* with possible dot products in the numerator (p1.Qj)^(-aj), ...
*
* So the integrand looks like
*
*   1/p1.p1^(a1+e1*ep)
*    /p2.p2^(a2+e2*ep)
*    /p1.Q3^(a3) ...
*
* Dot products containing p2 are NOT allowed.
*
* The trivial spherical factor per one-loop
*
*   1/(4 pi)^2
*
* is omitted as a convention.
*
* Note that we need more slots for the output than those for the input.
* For example, with 3 dot products, the indices are:
*
*           input               output
*           =====               ======
*   i1      p1.p1               Q.Q
*   i2      p2.p2               (removed)
*   id+0    p1.Q3               Q.Q3
*   id+1    p1.Q4               Q.Q4
*   id+2    p1.Q5 (last index)  Q.Q5
*   id+3                        Q3.Q4
*   id+4                        Q3.Q5
*   id+5                        Q4.Q5 (last index)
*   i3                          Q3.Q3
*   i4                          Q4.Q4
*   i5                          Q5.Q5
*
* In general, indices for M dot products affected by the integration are put
* in the last indices of Z function. Then, this procedure adds M(M-1)/2 indices
* for additionally generated dot products, removes the index specified by i2,
* and changes Z to Zb.
*
* Parameters:
*   Z         - the name of the integral function
*   i1,i2     - the positions of the indices for the edges depicted above.
*               i1 is reused for Q.Q in the output
*   e1,e2     - represent the non-integer parts of the two lines depicted
*               above
*   id        - the position of the index for the first dot product affected
*               by the integration (optional)
*   i3,i4,... - the positions of the indices for p3^2, p4^2, ... for the
*               output (optional)
*
* All the indices are specified by one-based numbering. e1 and e2 are expected
* be non-negative integers. No extra signs ("+") or no extra parentheses
* ("(" and ")") are allowed.
*
#procedure insertion(Z,i1,i2,e1,e2,?iv)
  #call Timing()

  #define PREFIX "FORCER"
  #define LABEL "`PREFIX'LABEL"

  #define n "`PREFIX'n"
  #define x "`PREFIX'x"
  #define q "`PREFIX'q"

  #define n1 "`PREFIX'n1"
  #define n2 "`PREFIX'n2"
  #define n3 "`PREFIX'n3"
  #define isum "`PREFIX'isum1"
  #define G "`PREFIX'G"
  #define fnum "`PREFIX'fnum"
  #define del "`PREFIX'del"
  #define ftensor "`PREFIX'ftensor"
  #define distrib "`PREFIX'distrib"

  #define Q "`q'1"

* Parse the variadic arguments.
  #define M "0"
  #do i={`?iv',}
    #ifdef `i'
      #if `M' == 0
        #define id "`i'"
      #else
        #define i{`M'+2} "`i'"
      #endif
      #redefine M "{`M'+1}"
    #endif
  #enddo
  #if `M' == 0
    #define id "-9999"
  #endif
  #if `M' == 1
    #message Error: invalid number of variadic arguments of insertion
    #terminate
  #endif
  #if `M' >= 2
    #redefine M "{`M'-1}"
  #endif

* Determine the maximum index.
  #if `M' == 0
    #define N "`i1'"
    #if `i2' > `N'
      #redefine N "`i2'"
    #endif
  #else
    #define N "{`id'+`M'-1}"
  #endif

* TODO: Don't use the general formula when the number of dot products is small.
*       The special cases can be faster.

* TODO: Is it possible to optimize for the index notation?

  if (count(`Z',1));

    #do i=1,`N'
      #define nL`i' "`n'`i'?"
      #define nR`i' "`n'`i'"
    #enddo
    #redefine nR`i1' "0"
    #redefine nR`i2' "0"
*   the sum of the numbers of (p1.qj)^n
    #redefine nsum   "0"
    #do j=1,`M'
      #define i "{`id'+`j'-1}"
        #redefine nL`i' "`n'`i'?neg0_"
        #redefine nR`i' "0"
        #redefine nsum "`nsum'-`n'`i'"
      #undefine i
    #enddo

    id ifnomatch->`LABEL'``LABEL'NO', `Z'(
        #do i=1,`N'
          `nL`i'',
        #enddo
      ?nn) = `Z'(
        #do i=1,`N'
          `nR`i'',
        #enddo
      ?nn) * `fnum'(
        #do j=1,`M'
          #define i "{`id'+`j'-1}"
            -`n'`i',
          #undefine i
        #enddo
      `n'`i1',`n'`i2',`nsum');
    #do i=3,{`M'+2}
      repeat id `fnum'(`n1'?pos_,?nn) = `fnum'(`n1'-1,?nn,`q'`i');
      id `fnum'(0,?nn) = `fnum'(?nn);
    #enddo
*   Now we have fnum(n1,n2,nsum,q1,...).
*   [dal(p)/4]^n -> fac_(n) / 2^n * distrib(1,2*n,dd,ftensor,...)
    id `fnum'(`n1'?,`n2'?,`n3'?,?pp) =
      sum_(`isum',0,integer_(`n3'/2),
           `G'(`n1',`e1',`n2',`e2',`n3',`isum')
           * `distrib'(1,2*`isum',`del',`ftensor',?pp)
           * `Q'.`Q'^`isum'
           / 2^`isum')
      * `Q'.`Q'^2
      / `Q'.`Q'^`n1'
      / `Q'.`Q'^`n2';
    id `distrib'(?a) = distrib_(?a);
    tovector `ftensor',`Q';
    id `del'(?a) = dd_(?a);

    #do i=1,`N'
      #define nL`i' "`n'`i'?"
      #define nR`i' "`n'`i'"
    #enddo
    #redefine nL`i1' "0"
    #redefine nR`i1' "-`x'`i1'"
    #redefine nL`i2' "0"
    #redefine nR`i2' "nan"
    #do j=1,`M'
      #define i "{`id'+`j'-1}"
        #redefine nL`i' "0"
        #redefine nR`i' "-`x'`i'"
      #undefine i
    #enddo

    id `Z'(
        #do i=1,`N'
          `nL`i'',
        #enddo
      ?nn) * `Q'.`Q'^`x'`i1'?
      #do j=1,`M'
        #define i "{`id'+`j'-1}"
          * `Q'.`q'{`j'+2}^`x'`i'?
        #undefine i
      #enddo
      #define i "{`id'+`M'}"
        #do j1=3,{`M'+2}
          #do j2={`j1'+1},{`M'+2}
            * `q'`j1'.`q'`j2'^`x'`i++'?
          #enddo
        #enddo
      #undefine i
      #do j=3,{`M'+2}
        #if `i`j'' >= 1
          * `q'`j'.`q'`j'^`x'`i`j''?
          #redefine nR`i`j'' "'nR`i`j'''-`x'`i`j''"
        #endif
      #enddo
      = `Z'b(
        #do i=1,`N'
          #if `i' != `i2'
            `nR`i'',
          #endif
        #enddo
        #define i "{`id'+`M'}"
          #do j1=3,{`M'+2}
            #do j2={`j1'+1},{`M'+2}
              -`x'`i++',
            #enddo
          #enddo
        #undefine i
      ?nn);  * Note that ?nn must be empty if M >= 1.
    #do j=3,{`M'+2}
      #if `i`j'' == 0
        #if `SETQONE'
          id `q'`j'.`q'`j'^`x'1? = 1;
        #else
          id `q'`j'.`q'`j'^`x'1? = Q2^`x'1;
        #endif
        #breakdo
      #endif
    #enddo
label `LABEL'``LABEL'NO++';
  endif;
  B `G',GschemeConstants;
  .sort:insertion-`Z'-1;
  Keep Brackets;
  #call reduceG()
  .sort:insertion-`Z'-2;

  #call Timing(Insertion(`Z'))
#endprocedure

*--#] insertion : 
*--#[ carpet :

**
* Performs massless one-loop outer integrals with blobs:
*
*                   -----
*                  /*****\
*               ---*******---
*     a1+e1*ep /   \*****/   \
*             v p1  -----     ^
*            /                 \
*     Q --<--                   --<-- Q
*            \                 /
*             \               /
*              \   a2+e2*ep  /
*               ------<------
*                     p2
*
* with possible dot products of inner loop momentum pj and Q, (pj.Q)^(-bj), ...
*
* So the integrand looks like
*
*   1/p1.p1^(a1+e1*ep)
*    /p2.p2^(a2+e2*ep)
*    /p3.p3^(a3+e3*ep) ...
*    /p3.Q^(b3) ...
*
* Dot products containing p1 or p2 are NOT allowed.
*
* The trivial spherical factor per one-loop
*
*   1/(4 pi)^2
*
* is omitted as a convention.
*
* Note that we need more slots for the output than those for the input.
* For example, with 3 dot products of an inner loop momentum and Q, the indices
* are:
*
*           input               output
*           =====               ======
*   i1      p1.p1               Q.Q (unless i1 = 0)
*   i2      p2.p2               (removed)
*   id+0    p3.Q                p3.Q
*   id+1    p4.Q                p4.Q
*   id+2    p5.Q (last index)   p5.Q
*   id+3                        p3.p4
*   id+4                        p3.p5
*   id+5                        p4.p5 (last index)
*   i3                          p3.p3
*   i4                          p4.p4
*   i5                          p5.p5
*
* In general, indices for M dot products of an inner loop momenta and Q are put
* in the last indices of Z function. All the indices except i1 and i2 before id
* are regarded as indices for squared inner loop momenta and dot products of
* two inner loop momenta (the formula needs all the indices). Then, this
* procedure adds M(M-1)/2 indices for additionally generated dot products,
* removes the index specified by i2, and changes Z to Zb.
*
* Parameters:
*   Z         - the name of the integral function
*   i1,i2     - the positions of the indices for the edges depicted above
*               (i1 is for the line with the blob). i1 is reused for Q.Q in
*               the output. If i1 is specified as 0, meaning no line
*               corresponding to this, and Q.Q^n is directly generated in
*               the output
*   e1,e2     - represent the non-integer parts of the two lines depicted
*               above
*   eb        - represent the sum of non-integer parts of lines in the blob
*   nl        - the number of loops in the blob
*   id        - the position of the index for the first dot product of an
*               internal loop momentum and Q (or the last index + 1 if no dot
*               products)
*   i3,i4,... - the positions of the indices for p3^2, p4^2, ... for the
*               output (optional)
*
* All the indices are specified by one-based numbering. e1, e2, eb and nl are
* expected be non-negative integers. No extra signs ("+") or no extra
* parentheses ("(" and ")") are allowed.
*
#procedure carpet(Z,i1,i2,e1,e2,eb,nl,id,?iv)
  #call Timing()

  #define PREFIX "FORCER"
  #define LABEL "`PREFIX'LABEL"

  #define n "`PREFIX'n"
  #define x "`PREFIX'x"
  #define q "`PREFIX'q"

  #define n1 "`PREFIX'n1"
  #define n2 "`PREFIX'n2"
  #define n3 "`PREFIX'n3"
  #define n4 "`PREFIX'n4"
  #define isum1 "`PREFIX'isum1"
  #define isum2 "`PREFIX'isum2"
  #define G "`PREFIX'G"
  #define fnum "`PREFIX'fnum"
  #define del "`PREFIX'del"
  #define ftensor "`PREFIX'ftensor"
  #define distrib "`PREFIX'distrib"
  #define poch "`PREFIX'poch"
  #define pochinv "`PREFIX'pochinv"

  #define Q "`q'1"

* Parse the variadic arguments.
  #define M "0"
  #do i={`?iv',}
    #ifdef `i'
      #define i{`M'+3} "`i'"
      #redefine M "{`M'+1}"
    #endif
  #enddo

* Determine the maximum index.
  #if `M' == 0
    #define N "`i1'"
    #if `i2' > `N'
      #redefine N "`i2'"
    #endif
  #else
    #define N "{`id'+`M'-1}"
  #endif

* TODO: Don't use the general formula when the number of dot products is small.
*       The special cases can be faster.

* TODO: Is it possible to optimize for the index notation?

  if (count(`Z',1));

    #do i=1,`N'
      #define nL`i' "`n'`i'?"
      #define nR`i' "`n'`i'"
    #enddo
    #redefine nR`i1' "0"
    #redefine nR`i2' "0"
*   the sum of the numbers of (pj.Q)^n
    #redefine nsum   "0"
    #do i=`id',{`id'+`M'-1}
      #redefine nL`i' "`n'`i'?neg0_"
      #redefine nR`i' "0"
      #redefine nsum "`nsum'-`n'`i'"
    #enddo
*   the sum of the numbers of 1/pj.pj^n
    #define dsum "0"
    #do i=1,{`id'-1}
      #if (`i' != `i1') && (`i' != `i2')
        #redefine dsum "`dsum'+`n'`i'"
      #endif
    #enddo

    id ifnomatch->`LABEL'``LABEL'NO', `Z'(
        #do i=1,`N'
          `nL`i'',
        #enddo
      ?nn) = `Z'(
        #do i=1,`N'
          `nR`i'',
        #enddo
      ?nn) * `fnum'(
        #do i=`id',{`id'+`M'-1}
          -`n'`i',
        #enddo
        #if `i1' == 0
          0,
        #else
          `n'`i1',
        #endif
      `n'`i2',`nsum',`dsum');
    #do i=3,{`M'+2}
      repeat id `fnum'(`n1'?pos_,?nn) = `fnum'(`n1'-1,?nn,`q'`i');
      id `fnum'(0,?nn) = `fnum'(?nn);
    #enddo

*   Now we have fnum(n1,n2,nsum,dsum,p1,...).
*   [dal(Q)/4]^n -> fac_(n) / 2^n * distrib(1,2*n,dd,ftensor,...)
    id `fnum'(`n1'?,`n2'?,`n3'?,`n4'?,?pp) =
      1 / Q2^`n1' / Q2^`n2' * Q2^2
      * sum_(`isum1',0,integer_(`n3'/2),
            `G'(`n1'+`n4'-`isum1'-{2*`nl'},{`e1'+`eb'+`nl'},`n2',`e2',`n3'-2*`isum1',0)
            * `poch'(-`isum1',2-ep+`n3'-`isum1')
            * `Q'.`Q'^`isum1'
            * invfac_(`isum1')
            / 2^`isum1'
            * sum_(`isum2',0,integer_(`n3'/2)-`isum1',
                   `distrib'(1,2*(`isum1'+`isum2'),`del',`ftensor',?pp) * fac_(`isum1'+`isum2')
                   * `poch'(-`isum2',1-ep+`n3'-2*`isum1')
                   * sign_(`isum2')
                   * `Q'.`Q'^`isum2'
                   * invfac_(`isum2')
                   / 2^`isum2'
            )
       );
    id `distrib'(?a) = distrib_(?a);
    tovector `ftensor',`Q';
    id `del'(?a) = dd_(?a);

    #do i=1,`N'
      #define nL`i' "`n'`i'?"
      #define nR`i' "`n'`i'"
    #enddo
    #redefine nL`i1' "0"
    #redefine nR`i1' "-`x'`i1'"
    #redefine nL`i2' "0"
    #redefine nR`i2' "nan"
    #do i=`id',{`id'+`M'-1}
      #redefine nL`i' "0"
      #redefine nR`i' "-`x'`i'"
    #enddo

    id `Z'(
        #do i=1,`N'
          `nL`i'',
        #enddo
      ?nn)
      #if `i1' != 0
        * `Q'.`Q'^`x'`i1'?
      #endif
      #do j=1,`M'
        #define i "{`id'+`j'-1}"
          * `Q'.`q'{`j'+2}^`x'`i'?
        #undefine i
      #enddo
      #define i "{`id'+`M'}"
        #do j1=3,{`M'+2}
          #do j2={`j1'+1},{`M'+2}
            * `q'`j1'.`q'`j2'^`x'`i++'?
          #enddo
        #enddo
      #undefine i
      #do j=3,{`M'+2}
        #if `i`j'' >= 1
          * `q'`j'.`q'`j'^`x'`i`j''?
          #redefine nR`i`j'' "'nR`i`j'''-`x'`i`j''"
        #endif
      #enddo
      = `Z'b(
        #do i=1,`N'
          #if `i' != `i2'
            `nR`i'',
          #endif
        #enddo
        #define i "{`id'+`M'}"
          #do j1=3,{`M'+2}
            #do j2={`j1'+1},{`M'+2}
              -`x'`i++',
            #enddo
          #enddo
        #undefine i
      ?nn);  * Note that ?nn must be empty if M >= 1.
    #if `i1' == 0
      #if `SETQONE'
        id `Q'.`Q'^`x'1? = 1;
      #else
        id `Q'.`Q'^`x'1? = Q2^`x'1;
      #endif
    #endif

*   #define n "`PREFIX'n1"
    #define x "`PREFIX'x1"
    id `poch'(0,`x'?) = 1;
*   id `poch'(`n'?pos_,`x'?pos_) = fac_(`n'-1+`x') * invfac_(`x'-1);
label `LABEL'``LABEL'NO++';
  endif;
  #if `SETQONE'
	id	Q2^n? = 1;
	id	Q.Q^n? = 1;
  #else
	#message  Warning: SETQONE = `SETQONE'
  #endif
  B `poch',`G';
  .sort:carpet-`Z'-1;
  Keep Brackets;
  #call reduceG()
  .sort:carpet-`Z'-2;

  #call Timing(Carpet(`Z'))
#endprocedure

*--#] carpet : 
*--#[ reduceG :

**
* Handles the G function:
*
*                          Gamma(a+b-s-D/2) Gamma(D/2-a+n-s) Gamma(D/2-b+s)
*   G(a,b,n,s) = (4 pi)^ep ------------------------------------------------
*                                 Gamma(a) Gamma(b) Gamma(D-a-b+n)
*
* The input convention:
*
*   `PREFIX'G(n1,x1,n2,x2,n3,n4) = G(n1+x1*ep,n2+x2*ep,n3,n4)
*
* where n1,...,n4,x1,x2 are integers.
*
* The output convention:
*
*   GschemeConstants(x1,x2) = (1+x1+x2) G(1+x1*ep,1+x2*ep,0,0)
*
* Note that GschemeConstants(0,0) is set to 1/ep in the G-scheme.
*
#procedure reduceG
  #define PREFIX "FORCER"
  #define LABEL "`PREFIX'LABEL"

  #define n "`PREFIX'n1"
  #define x "`PREFIX'x1"
  #define n1 "`PREFIX'n1"
  #define n2 "`PREFIX'n2"
  #define n3 "`PREFIX'n3"
  #define n4 "`PREFIX'n4"
  #define x1 "`PREFIX'x1"
  #define x2 "`PREFIX'x2"
  #define G "`PREFIX'G"
  #define poch "`PREFIX'poch"
  #define pochinv "`PREFIX'pochinv"
  #define num "`PREFIX'num"
  #define den "`PREFIX'den"

* NOTE: poch(n,x) here is poch(n,x) := Gamma(n+x)/Gamma(x).

  id ifnomatch->`LABEL'``LABEL'NO', `G'(`n1'?,`x1'?,`n2'?,`x2'?,`n3'?,`n4'?)
    = GschemeConstants(`x1',`x2') / (1+`x1'+`x2')
    * `poch'(`n1'+`n2'-`n4'-2,ep+`x1'*ep+`x2'*ep)
    * `poch'(1-`n1'+`n3'-`n4',1-ep-`x1'*ep)
    * `poch'(1-`n2'+`n4',1-ep-`x2'*ep)
    * `pochinv'(`n1'-1,1+`x1'*ep)
    * `pochinv'(`n2'-1,1+`x2'*ep)
    * `pochinv'(2-`n1'-`n2'+`n3',2-2*ep-`x1'*ep-`x2'*ep);

  repeat id `poch'(`n'?pos_,`x'?) = `poch'(`n'-1,`x')*`num'(`n'-1+`x');
  repeat id `poch'(`n'?neg_,`x'?) = `poch'(`n'+1,`x')*`den'(`n'+`x');
  repeat id `pochinv'(`n'?pos_,`x'?) = `pochinv'(`n'-1,`x')*`den'(`n'-1+`x');
  repeat id `pochinv'(`n'?neg_,`x'?) = `pochinv'(`n'+1,`x')*`num'(`n'+`x');

  id GschemeConstants(`n1'?pos0_,`n2'?pos_) =
    + GschemeConstants(`n1',`n2') * theta_ (`n1'-`n2')
    + GschemeConstants(`n2',`n1') * thetap_(`n2'-`n1');

  id GschemeConstants(1,0)*GschemeConstants(2,1)
    = GschemeConstants(1,1)*GschemeConstants(3,0)*rat(1-4*ep,1-3*ep);
  id GschemeConstants(0,0)*GschemeConstants(2,1)
    = GschemeConstants(2,0)*GschemeConstants(3,0)*rat(1-4*ep,1-2*ep);
  id GschemeConstants(1,1)*GschemeConstants(0,0)
    = GschemeConstants(1,0)*GschemeConstants(2,0)*rat(1-3*ep,1-2*ep);

  id `poch'(0,`x'?) = 1;
  id `pochinv'(0,`x'?) = 1;
* id `num'(`x'?)*`den'(`x'?) = 1;
* id `num'(`x'?number_) = `x';
* id `den'(`x'?number_) = 1/`x';
  id `num'(`x'?) = rat(`x',1);
  id `den'(`x'?) = rat(1,`x');
label `LABEL'``LABEL'NO++';
#endprocedure

*--#] reduceG : 
*--#[ triangleS :

**
* Applies the triangle rule until one of the edges vanishes. The assignment of
* indices is as follows:
*
*         n4         n3         n5
*       ----------------------------
*               \          /
*                \        /
*        n1+e1*ep \      / n2+e2*ep
*                  \    /
*                   \  /
*                    \/
*                    /\
*
* with possible dot products in the numerator (p3.Q1)^(-nK1), ...
*
* The formula is not sensitive to the direction of p3.
*
* Dot products containing p1 or p2 are NOT allowed.
*
* As the input, each factor is represented by a symbol like i1 := p1.p1
* multiplied by a label for the integration Z, so it looks like
*
*   Z / i1^(n1+e1*ep)
*     / i2^(n2+e2*ep)
*     / i3^n3
*     / i4^n4
*     / i5^n5
*     / iK1^nK1 ...
*
* Note that the non-integer parts e3,...,e5 are assumed to be zero, because
* otherwise the triangle rule can not be applied. (Actually e4 and e5 are not
* needed due to no their appearance in the right-hand side.) eK's are also
* assumed to be zero. It may be possible to introduce eK's as well, but
* in practice we do not need non-integer parts for dot products if they are in
* the numerator. The caller must be sure that integrals passed to this
* procedure do not have non-integer parts e3,...,e5,eK1,... = 0.
*
* This procedure will change the indices (their powers with a sign flip)
* specified by i1,...,i5. n1 and n2 may be increased while n3, n4 and n5 are
* decreased until one of them becomes zero. In the output, Z is replaced to Zb.
*
* Parameters:
*   Z           - the symbol representing the integral
*   i1,...,i5   - the symbols representing the propagators depicted above
*   iK1,iK2,... - the symbols for the irreducible numerators (optional)
*   e1,e2       - non-negative integers representing the non-integer parts of
*                 the two lines depicted above
*
* e1 and e2 are expected be non-negative integers. No extra signs ("+") or no
* extra parentheses ("(" and ")") are allowed.
*
* Note that 5 + n symbols are needed to be specified for the triangle rule with
* n dot products.
*
* For example, with 3 dot products, the symbols are:
*
*        input/output
*        =====
*   i1   p1.p1
*   i2   p2.p2
*   i3   p3.p3
*   i4   p4.p4
*   i5   p5.p5
*   i6   p3.Q1
*   i7   p3.Q2
*   i8   p3.Q3
*
#ifndef `TRIANGLESMETHOD'
  #define TRIANGLESMETHOD "3"
#endif

#procedure triangleS(Z,?a)
  #call Timing()

  #call triangleS`TRIANGLESMETHOD'(`Z',`?a')

  #call Timing(Triangle(`Z'))
#endprocedure

*--#[ triangleS1 :

* An implementation by recursion.
#procedure triangleS1(Z,i1,i2,i3,i4,i5,?iKj,e1,e2)
  #define prefix "FORCER"
  #define a1     "`prefix'a1"
  #define a2     "`prefix'a2"
  #define E      "`prefix'E"
  #define doloop "`prefix'doloop"

  #define BRA "`Z',`i1',`i2',`i3',`i4',`i5'"
  #do i={`?iKj',}
    #ifdef `i'
      #redefine BRA "`BRA',`i'"
    #endif
  #enddo

  #define BRA "B `BRA';"
  #define KB  "Keep Brackets;"

  #define nloop "1"
  `BRA'
  .sort:triangleS-`Z'-`nloop++';
  #do loop=1,1
    `KB'
    #$`doloop' = 0;
    if (count(`Z',1));
      if (count(`i3',-1) >= 1);
        if (count(`i4',-1) >= 1);
          if (count(`i5',-1) >= 1);
            $`a1' = count_(`i1',-1)+`e1'*ep;
            $`a2' = count_(`i2',-1)+`e2'*ep;
            $`E' = 4-2*ep-$`a1'-$`a2'
              #do i={`?iKj',}
                #ifdef `i'
                  + count_(`i',1)
                #endif
              #enddo
              + count_(`i3',2)
            ;
            multiply rat($`a1',$`E') / `i1' * (`i3' - `i4')
                   + rat($`a2',$`E') / `i2' * (`i3' - `i5');
            if (count(`i3',-1) >= 1)
              if (count(`i4',-1) >= 1)
                if (count(`i5',-1) >= 1)
                  $`doloop' = 1;
          endif;
        endif;
      endif;
    endif;
    `BRA'
    ModuleOption maximum, $`doloop';
    ModuleOption local, $`a1', $`a2', $`E';
    .sort:triangleS-`Z'-`nloop++';
    #if `$doloop'
      #redefine loop "0"
    #endif
  #enddo
#endprocedure

*--#] triangleS1 : 
*--#[ triangleS2 :

* An implementation by the summation formula.
#procedure triangleS2(Z,i1,i2,i3,i4,i5,?iKj,e1,e2)
  #define prefix "FORCER"

  #define n1 "`prefix'n1"
  #define n2 "`prefix'n2"
  #define n3 "`prefix'n3"
  #define n4 "`prefix'n4"
  #define n5 "`prefix'n5"

  #define k1 "`prefix'isum1"
  #define k2 "`prefix'isum2"
  #define k3 "`prefix'isum3"
  #define k4 "`prefix'isum4"

  #define poch "`prefix'poch"
  #define num  "`prefix'num"
  #define den  "`prefix'den"

  #define a1 "`prefix'a1"
  #define a2 "`prefix'a2"
  #define E  "`prefix'E"

  if (count(`Z',1));
    if (count(`i3',-1) >= 1);
      if (count(`i4',-1) >= 1);
        if (count(`i5',-1) >= 1);
          $`a1' = count_(`i1',-1)+`e1'*ep;
          $`a2' = count_(`i2',-1)+`e2'*ep;
*         D + N - a1 - a2
          $`E' = 4-2*ep-$`a1'-$`a2'
            #do i={`?iKj',}
              #ifdef `i'
                + count_(`i',1)
              #endif
            #enddo
          ;
          id 1 / `i1'^`n1'?
               / `i2'^`n2'?
               / `i3'^`n3'?
               / `i4'^`n4'?
               / `i5'^`n5'? =
            + sum_(`k1',0,`n3',sum_(`k2',0,`n4'-1,sum_(`k4',0,`n5'-1,
                   1 / `i1'^(`n1'+`k1'+`k2')
                     / `i2'^(`n2'+`n3'-`k1'+`k4')
                     / `i4'^(`n4'-`k2')
                     / `i5'^(`n5'-`k4')
                   * sign_(`k2'+`k4')
                   * `n3'
                   * fac_(`n3'+`k2'+`k4'-1)
                   * invfac_(`k1')
                   * invfac_(`k2')
                   * invfac_(`k4')
                   * invfac_(`n3'-`k1')
                   * `poch'(-`n3'-`k2'-`k4',$`E'-`n3')
                   * `poch'(`k1'+`k2',$`a1')
                   * `poch'(`n3'-`k1'+`k4',$`a2')
            )))
            + sum_(`k1',0,`n3'-1,sum_(`k3',0,`n3'-1-`k1',sum_(`k4',0,`n5'-1,
                   1 / `i1'^(`n1'+`n4'+`k1')
                     / `i2'^(`n2'+`k3'+`k4')
                     / `i3'^(`n3'-`k1'-`k3')
                     / `i5'^(`n5'-`k4')
                   * sign_(`n4'+`k4')
                   * fac_(`n4'+`k1'+`k3'+`k4'-1)
                   * invfac_(`k1')
                   * invfac_(`k3')
                   * invfac_(`k4')
                   * invfac_(`n4'-1)
                   * `poch'(-`n4'-`k1'-`k3'-`k4',$`E'-2*`n3'+`k1'+`k3'+1)
                   * `poch'(`n4'+`k1',$`a1')
                   * `poch'(`k3'+`k4',$`a2')
            )))
            + sum_(`k1',0,`n3'-1,sum_(`k2',0,`n4'-1,sum_(`k3',0,`n3'-1-`k1',
                   1 / `i1'^(`n1'+`k1'+`k2')
                     / `i2'^(`n2'+`n5'+`k3')
                     / `i3'^(`n3'-`k1'-`k3')
                     / `i4'^(`n4'-`k2')
                   * sign_(`n5'+`k2')
                   * fac_(`n5'+`k1'+`k2'+`k3'-1)
                   * invfac_(`k1')
                   * invfac_(`k2')
                   * invfac_(`k3')
                   * invfac_(`n5'-1)
                   * `poch'(-`n5'-`k1'-`k2'-`k3',$`E'-2*`n3'+`k1'+`k3'+1)
                   * `poch'(`k1'+`k2',$`a1')
                   * `poch'(`n5'+`k3',$`a2')
            )))
            ;

          #define n "`prefix'n1"
          #define x "`prefix'x1"

          id `poch'(0,`x'?) = 1;
          id `poch'(`n'?pos_,`x'?pos_) = fac_(`n'-1+`x') * invfac_(`x'-1);
        endif;
      endif;
    endif;
  endif;
  B `poch';
  ModuleOption local $`a1',$`a2',$`E';
  .sort:triangleS-`Z'-1;
  Keep Brackets;
  repeat id `poch'(`n'?pos_,`x'?) = `poch'(`n'-1,`x')*`num'(`n'-1+`x');
  repeat id `poch'(`n'?neg_,`x'?) = `poch'(`n'+1,`x')*`den'(`n'+`x');
  id `poch'(0,`x'?) = 1;
* id `num'(`x'?)*`den'(`x'?) = 1;
* id `num'(`x'?number_) = `x';
* id `den'(`x'?number_) = 1/`x';
  id `num'(`x'?) = rat(`x',1);
  id `den'(`x'?) = rat(1,`x');
  .sort:triangleS-`Z'-2;
#endprocedure

*--#] triangleS2 : 
*--#[ triangleS3 :

* An implementation by the recursion/summation formula.
#procedure triangleS3(Z,i1,i2,i3,i4,i5,?iKj,e1,e2)
  #define prefix "FORCER"

  #define n1 "`prefix'n1"
  #define n2 "`prefix'n2"
  #define n3 "`prefix'n3"
  #define n4 "`prefix'n4"
  #define n5 "`prefix'n5"

  #define k1 "`prefix'isum1"
  #define k2 "`prefix'isum2"
  #define k3 "`prefix'isum3"
  #define k4 "`prefix'isum4"

  #define poch "`prefix'poch"
  #define num  "`prefix'num"
  #define den  "`prefix'den"

  #define a1 "`prefix'a1"
  #define a2 "`prefix'a2"
  #define N  "`prefix'E"
  #define E1 "`prefix'E"
  #define haspoch "`prefix'haspoch"

#ifdef `FORCERSTATS'
  #$`prefix'nrec = 0;
  #$`prefix'nsum = 0;
#endif

  #$`haspoch' = 0;

  if (count(`Z',1));
    if (count(`i3',-1) >= 1);
      if (count(`i4',-1) >= 1);
        if (count(`i5',-1) >= 1);
          if (count(`i3',-1,`i4',-1,`i5',-1) <= 4);
#ifdef `FORCERSTATS'
            $`prefix'nrec = $`prefix'nrec + 1;
#endif
            $`N' =
              #do i={`?iKj',}
                #ifdef `i'
                  + count_(`i',1)
                #endif
              #enddo
            ;
            repeat id 1 / `i1'^`n1'?
                        / `i2'^`n2'?
                        / `i3'^`n3'?pos_
                        / `i4'^`n4'?pos_
                        / `i5'^`n5'?pos_ =
              1 / `i1'^`n1'
                / `i2'^`n2'
                / `i3'^`n3'
                / `i4'^`n4'
                / `i5'^`n5' * (
                  + rat(`n1'+`e1'*ep,1) * (`i3' - `i4') / `i1'
                  + rat(`n2'+`e2'*ep,1) * (`i3' - `i5') / `i2'
              ) * rat(1,4+$`N'-`n1'-`n2'-2*`n3'-{2+`e1'+`e2'}*ep);
          else;
#ifdef `FORCERSTATS'
            $`prefix'nsum = $`prefix'nsum + 1;
#endif
            $`a1' = count_(`i1',-1)+`e1'*ep;
            $`a2' = count_(`i2',-1)+`e2'*ep;
*           E1 = D + N - a1 - a2 (not including -2*b)
            $`E1' = 4-2*ep-$`a1'-$`a2'
              #do i={`?iKj',}
                #ifdef `i'
                  + count_(`i',1)
                #endif
              #enddo
            ;
            id 1 / `i1'^`n1'?
                 / `i2'^`n2'?
                 / `i3'^`n3'?
                 / `i4'^`n4'?
                 / `i5'^`n5'? =
              + sum_(`k1',0,`n3',sum_(`k2',0,`n4'-1,sum_(`k4',0,`n5'-1,
                     1 / `i1'^(`n1'+`k1'+`k2')
                       / `i2'^(`n2'+`n3'-`k1'+`k4')
                       / `i4'^(`n4'-`k2')
                       / `i5'^(`n5'-`k4')
                     * sign_(`k2'+`k4')
                     * `n3'
                     * fac_(`n3'+`k2'+`k4'-1)
                     * invfac_(`k1')
                     * invfac_(`k2')
                     * invfac_(`k4')
                     * invfac_(`n3'-`k1')
                     * `poch'(-`n3'-`k2'-`k4',$`E1'-`n3')
                     * `poch'(`k1'+`k2',$`a1')
                     * `poch'(`n3'-`k1'+`k4',$`a2')
              )))
              + sum_(`k1',0,`n3'-1,sum_(`k3',0,`n3'-1-`k1',sum_(`k4',0,`n5'-1,
                     1 / `i1'^(`n1'+`n4'+`k1')
                       / `i2'^(`n2'+`k3'+`k4')
                       / `i3'^(`n3'-`k1'-`k3')
                       / `i5'^(`n5'-`k4')
                     * sign_(`n4'+`k4')
                     * fac_(`n4'+`k1'+`k3'+`k4'-1)
                     * invfac_(`k1')
                     * invfac_(`k3')
                     * invfac_(`k4')
                     * invfac_(`n4'-1)
                     * `poch'(-`n4'-`k1'-`k3'-`k4',$`E1'-2*`n3'+`k1'+`k3'+1)
                     * `poch'(`n4'+`k1',$`a1')
                     * `poch'(`k3'+`k4',$`a2')
              )))
              + sum_(`k1',0,`n3'-1,sum_(`k2',0,`n4'-1,sum_(`k3',0,`n3'-1-`k1',
                     1 / `i1'^(`n1'+`k1'+`k2')
                       / `i2'^(`n2'+`n5'+`k3')
                       / `i3'^(`n3'-`k1'-`k3')
                       / `i4'^(`n4'-`k2')
                     * sign_(`n5'+`k2')
                     * fac_(`n5'+`k1'+`k2'+`k3'-1)
                     * invfac_(`k1')
                     * invfac_(`k2')
                     * invfac_(`k3')
                     * invfac_(`n5'-1)
                     * `poch'(-`n5'-`k1'-`k2'-`k3',$`E1'-2*`n3'+`k1'+`k3'+1)
                     * `poch'(`k1'+`k2',$`a1')
                     * `poch'(`n5'+`k3',$`a2')
              )))
            ;

            #define n "`prefix'n1"
            #define x "`prefix'x1"

            id `poch'(0,`x'?) = 1;
            id `poch'(`n'?pos_,`x'?pos_) = fac_(`n'-1+`x') * invfac_(`x'-1);
            $`haspoch' = 1;
          endif;
        endif;
      endif;
    endif;
  endif;
  B `poch';
  ModuleOption local $`a1',$`a2',$`E1';
  ModuleOption maximum $`haspoch';
#ifdef `FORCERSTATS'
  ModuleOption sum $`prefix'nrec, $`prefix'nsum;
#endif
  .sort:triangleS-`Z'-1;
#ifdef `FORCERSTATS'
  #message Stats: Triangle(`Z'): nrec=`$`prefix'nrec', nsum=`$`prefix'nsum'
#endif
  #if `$`haspoch''
    Keep Brackets;
    repeat id `poch'(`n'?pos_,`x'?) = `poch'(`n'-1,`x')*`num'(`n'-1+`x');
    repeat id `poch'(`n'?neg_,`x'?) = `poch'(`n'+1,`x')*`den'(`n'+`x');
    id `poch'(0,`x'?) = 1;
*   id `num'(`x'?)*`den'(`x'?) = 1;
*   id `num'(`x'?number_) = `x';
*   id `den'(`x'?number_) = 1/`x';
    id `num'(`x'?) = rat(`x',1);
    id `den'(`x'?) = rat(1,`x');
    .sort:triangleS-`Z'-2;
  #endif
#endprocedure

*--#] triangleS3 : 
*--#] triangleS : 
*--#[ diamondS :

**
* Applies the diamond rule until one of edges vanishes.
*
* Parameters:
*   Z           - the symbol representing the integral
*   n           - the number of A- (upper), B- (lower) and C- (external) lines.
*                 Must be an integer >= 3
*   k           - the number of S- (spectator) lines. Must be an integer >= 0
*   a1,...,an   - the symbols for A- (upper) lines
*   b1,...,bn   - the symbols for B- (lower) lines
*   c1,...,cn   - the symbols for C- (external) lines
*   s1,...,sk   - the symbols for S- (spectator) lines
*   d1,...,dd   - the symbols for dot products containing one of B- or S-lines.
*                 The number "d" of them is determined by the number of
*                 arguments passed to this procedure. If there are dot products
*                 consists of two vectors from B- and S-lines, they need to be
*                 specified twice.
*   e1,...,en   - non-negative integers representing the non-integer parts of
*                 A-lines.
*
* Remarks:
*   In general, two of external lines can be connected, i.e., the rule allows
*   duplications in c1,...,cn, twice at most. But the 2nd method uses the
*   explicit summation formula that is valid only if no external lines are
*   connected.
*
#ifndef `DIAMONDSMETHOD'
  #define DIAMONDSMETHOD "3"
#endif

#procedure diamondS(Z,n,k,?va)
  #call Timing()

  #if (`n' >= 3) == 0
    #message Error: diamondS is called with n=`n'.
    #terminate
  #endif
  #if (`k' >= 0) == 0
    #message Error: diamondS is called with k=`k'.
    #terminate
  #endif
* Parse the variadic arguments.
  #define i "0"
  #define j "0"
  #define p0 "a`~j'"
  #define p1 "b`~j'"
  #define p2 "c`~j'"
  #define p3 "s`~j'"
  #define p4 "d`~j'"
  #define n0 "`n'"
  #define n1 "`n'"
  #define n2 "`n'"
  #define n3 "`k'"
  #define n4 "9999"
  #do a={`?va',}
    #ifdef `a'
      #redefine j "{`j'+1}"
      #define `p`i'' "`a'"
      #if `j' == `n`i''
        #if `n{`i'+1}' >= 1
          #redefine i "{`i'+1}"
        #else
*         For the case that k = 0.
          #redefine i "{`i'+2}"
        #endif
        #redefine j "0"
      #endif
    #endif
  #enddo
  #if (`i' < 4) || (`j' < `n')
    #message Error: not enough argumentis: diamondS(`Z',`n',`k',`?va')
    #terminate
  #endif
  #define d "{`j'-`n'}"
  #do i=1,`n'
    #define e`i' "`d{`i'+`d'}'"
    #undefine d{`i'+`d'}
  #enddo

* Now we have n, k, d, a1,...,an, b1,...,bn, c1,...,cn, k1,...,kn, d1,...,dd
* and e1,...,en.

  #define BRA "`Z'"
  #do i=1,`n'
    #redefine BRA "`BRA',`a`i''"
  #enddo
  #do i=1,`n'
    #redefine BRA "`BRA',`b`i''"
  #enddo
  #do i=1,`n'
    #redefine BRA "`BRA',`c`i''"
  #enddo
  #do i=1,`k'
    #redefine BRA "`BRA',`s`i''"
  #enddo
  #do i=1,`d'
    #redefine BRA "`BRA',`d`i''"
  #enddo

  #define BRA "B `BRA';"
  #define KB  "Keep Brackets;"

  #call diamondS`DIAMONDSMETHOD'()

  #call Timing(Diamond(`Z'))
#endprocedure

*--#[ diamondS1 :

* An implementation by recursion.
#procedure diamondS1()
  #define prefix "FORCER"
  #define a      "`prefix'a"
  #define E      "`prefix'E"
  #define doloop "`prefix'doloop"

  #define nloop "1"
  `BRA'
  .sort:diamondS-`Z'-`nloop++';
  #do loop=1,1
    `KB'
    #$`doloop' = 0;
    if (count(`Z',1));
      #do i=1,`n'
        if (count(`b`i'',-1) >= 1);
          if (count(`c`i'',-1) >= 1);
      #enddo
            #do i=1,`n'
              $`a'`i' = count_(`a`i'',-1)+`e`i''*ep;
            #enddo
            $`E' = {`n'+`k'-1} * (4-2*ep)
              #do i=1,`n'
                   - $`a'`i'
              #enddo
              #do i=1,`n'
                   - 2 * count_(`b`i'',-1)
              #enddo
              #do i=1,`k'
                   - 2 * count_(`s`i'',-1)
              #enddo
              #do i=1,`d'
                   + count_(`d`i'',1)
              #enddo
            ;
            multiply
              #do i=1,`n'
                + rat($`a'`i',$`E') / `a`i'' * (`b`i'' - `c`i'')
              #enddo
            ;
            #do i=1,`n'
              if (count(`b`i'',-1) >= 1)
                if (count(`c`i'',-1) >= 1)
            #enddo
                  $`doloop' = 1;
      #do i=1,`n'
          endif;
        endif;
      #enddo
    endif;
    `BRA'
    ModuleOption maximum, $`doloop';
    ModuleOption local, <$`a'1>,...,<$`a'`n'>,$`E';
    .sort:diamondS-`Z'-`nloop++';
    #if `$`doloop''
      #redefine loop "0"
    #endif
  #enddo
#endprocedure

*--#] diamondS1 : 
*--#[ diamondS2 :

* An implementation by the summation formula.
#procedure diamondS2()
  #define prefix "FORCER"
  #define a      "`prefix'a"
  #define b      "`prefix'b"
  #define c      "`prefix'c"
  #define E      "`prefix'E"

  #define poch   "`prefix'poch"

  #define kb     "`prefix'kb"
  #define kc     "`prefix'kc"

  #define kbsum "(`kb'1+...+`kb'`n')"
  #define kcsum "(`kc'1+...+`kc'`n')"
  #define ksum  "(`kbsum'+`kcsum')"

  `BRA'
  .sort:diamondS-`Z'-1;
  S `kb'1,...,`kb'`n';
  S `kc'1,...,`kc'`n';
  `KB'

  if (count(`Z',1));
    #do i=1,`n'
      if (count(`b`i'',-1) >= 1);
        if (count(`c`i'',-1) >= 1);
    #enddo
          #do i=1,`n'
            $`a'`i' = count_(`a`i'',-1)+`e`i''*ep;
            $`b'`i' = count_(`b`i'',-1);
            $`c'`i' = count_(`c`i'',-1);
          #enddo
          $`E' = {`n'+`k'-1} * (4-2*ep)
            #do i=1,`n'
                 - $`a'`i' - 2 * $`b'`i'
            #enddo
            #do i=1,`k'
                 - 2 * count_(`s`i'',-1)
            #enddo
            #do i=1,`d'
                 + count_(`d`i'',1)
            #enddo
          ;
          multiply
            #do r=1,`n'
              +
              #do i=1,`n'
                #if `i' != `r'
                  sum_(`kb'`i',0,$`b'`i'-1,
                #else
                  sum_(`kb'`i',$`b'`i',$`b'`i',
                #endif
                sum_(`kc'`i',0,$`c'`i'-1,
              #enddo
                  sign_(`kcsum')
                  * `kb'`r'
                  * fac_(`ksum'-1)
                  #do i=1,`n'
                    * invfac_(`kb'`i')
                    * invfac_(`kc'`i')
                  #enddo
                  * `poch'(-`ksum',$`E'+`kbsum')
                  #do i=1,`n'
                    * `poch'(`kb'`i'+`kc'`i',$`a'`i')
                  #enddo
                  #do i=1,`n'
                    / `a`i''^`kb'`i'
                    / `a`i''^`kc'`i'
                    * `b`i''^`kb'`i'
                    * `c`i''^`kc'`i'
                  #enddo
              #do i=1,{2*`n'}
                )
              #enddo
              +
              #do i=1,`n'
                sum_(`kb'`i',0,$`b'`i'-1,
                #if `i' != `r'
                  sum_(`kc'`i',0,$`c'`i'-1,
                #else
                  sum_(`kc'`i',$`c'`i',$`c'`i',
                #endif
              #enddo
                  sign_(`kcsum')
                  * `kc'`r'
                  * fac_(`ksum'-1)
                  #do i=1,`n'
                    * invfac_(`kb'`i')
                    * invfac_(`kc'`i')
                  #enddo
                  * `poch'(-`ksum',$`E'+`kbsum'+1)
                  #do i=1,`n'
                    * `poch'(`kb'`i'+`kc'`i',$`a'`i')
                  #enddo
                  #do i=1,`n'
                    / `a`i''^`kb'`i'
                    / `a`i''^`kc'`i'
                    * `b`i''^`kb'`i'
                    * `c`i''^`kc'`i'
                  #enddo
              #do i=1,{2*`n'}
                )
              #enddo
            #enddo
          ;
          #define n "`prefix'n1"
          #define x "`prefix'x1"
          id `poch'(0,`x'?) = 1;
          id `poch'(`n'?pos_,`x'?pos_) = fac_(`n'-1+`x') * invfac_(`x'-1);
          #undefine n
          #undefine x
    #do i=1,`n'
        endif;
      endif;
    #enddo
  endif;
  B `poch';
  ModuleOption local <$`a'1>,...,<$`a'`n'>,
                     <$`b'1>,...,<$`b'`n'>,
                     <$`c'1>,...,<$`c'`n'>,$`E';
  .sort:diamondS-`Z'-2;
  Keep Brackets;
  #define num "`prefix'num"
  #define den "`prefix'den"
  #define n "`prefix'n1"
  #define x "`prefix'x1"
  repeat id `poch'(`n'?pos_,`x'?) = `poch'(`n'-1,`x')*`num'(`n'-1+`x');
  repeat id `poch'(`n'?neg_,`x'?) = `poch'(`n'+1,`x')*`den'(`n'+`x');
  id `poch'(0,`x'?) = 1;
* id `num'(`x'?)*`den'(`x'?) = 1;
* id `num'(`x'?number_) = `x';
* id `den'(`x'?number_) = 1/`x';
  id `num'(`x'?) = rat(`x',1);
  id `den'(`x'?) = rat(1,`x');
  .sort:diamondS-`Z'-3;
#endprocedure

*--#] diamondS2 : 
*--#[ diamondS3 :

* An implementation by the recursion/summation formula.
#procedure diamondS3()
  #define prefix "FORCER"
  #define a      "`prefix'a"
  #define b      "`prefix'b"
  #define c      "`prefix'c"
  #define E      "`prefix'E"
  #define E1     "`prefix'E"
  #define haspoch "`prefix'haspoch"

  #define poch   "`prefix'poch"

  #define ka     "`prefix'ka"
  #define kb     "`prefix'kb"
  #define kc     "`prefix'kc"

  #define kbsum "(`kb'1+...+`kb'`n')"
  #define kcsum "(`kc'1+...+`kc'`n')"
  #define ksum  "(`kbsum'+`kcsum')"
  #define esum  ""
  #do i=1,`n'
    #redefine esum "`esum'+`e`i''"
  #enddo

  `BRA'
  .sort:diamondS-`Z'-1;
  S `ka'1,...,`ka'`n';
  S `kb'1,...,`kb'`n';
  S `kc'1,...,`kc'`n';
  `KB'

#ifdef `FORCERSTATS'
  #$`prefix'nrec = 0;
  #$`prefix'nsum = 0;
#endif

  #$`haspoch' = 0;

  if (count(`Z',1));
    #do i=1,`n'
      if (count(`b`i'',-1) >= 1);
        if (count(`c`i'',-1) >= 1);
    #enddo
          if(count(
            #do i=1,`n'
              #if `i' >= 2
                ,
              #endif
              `b`i'',-1,
              `c`i'',-1
            #enddo
          ) <= {2*`n'+1});
#ifdef `FORCERSTATS'
            $`prefix'nrec = $`prefix'nrec + 1;
#endif
*           E1 = (L+S) D + Ni + Ri - 2 si
            $`E1' = {(`n'+`k'-1)*4} - {(`n'+`k'-1)*2+`esum'}*ep
              #do i=1,`k'
                - 2 * count_(`s`i'',-1)
              #enddo
              #do i=1,`d'
                + count_(`d`i'',1)
              #enddo
            ;
            repeat id 1
              #do i=1,`n'
                / `a`i''^`ka'`i'?
                / `b`i''^`kb'`i'?pos_
                / `c`i''^`kc'`i'?pos_
              #enddo
              = 1
              #do i=1,`n'
                / `a`i''^`ka'`i'
                / `b`i''^`kb'`i'
                / `c`i''^`kc'`i'
              #enddo
              * (
                #do i=1,`n'
                  + rat(`ka'`i' + `e`i''*ep,1) * (`b`i'' - `c`i'') / `a`i''
                #enddo
              ) * rat(1,$`E1'
                #do i=1,`n'
                  - `ka'`i' - 2*`kb'`i'
                #enddo
              );
          else;
#ifdef `FORCERSTATS'
            $`prefix'nsum = $`prefix'nsum + 1;
#endif
            #do i=1,`n'
              $`a'`i' = count_(`a`i'',-1)+`e`i''*ep;
              $`b'`i' = count_(`b`i'',-1);
              $`c'`i' = count_(`c`i'',-1);
            #enddo
            $`E' = {`n'+`k'-1} * (4-2*ep)
              #do i=1,`n'
                   - $`a'`i' - 2 * $`b'`i'
              #enddo
              #do i=1,`k'
                   - 2 * count_(`s`i'',-1)
              #enddo
              #do i=1,`d'
                   + count_(`d`i'',1)
              #enddo
            ;
            multiply
              #do r=1,`n'
                +
                #do i=1,`n'
                  #if `i' != `r'
                    sum_(`kb'`i',0,$`b'`i'-1,
                  #else
                    sum_(`kb'`i',$`b'`i',$`b'`i',
                  #endif
                  sum_(`kc'`i',0,$`c'`i'-1,
                #enddo
                    sign_(`kcsum')
                    * `kb'`r'
                    * fac_(`ksum'-1)
                    #do i=1,`n'
                      * invfac_(`kb'`i')
                      * invfac_(`kc'`i')
                    #enddo
                    * `poch'(-`ksum',$`E'+`kbsum')
                    #do i=1,`n'
                      * `poch'(`kb'`i'+`kc'`i',$`a'`i')
                    #enddo
                    #do i=1,`n'
                      / `a`i''^`kb'`i'
                      / `a`i''^`kc'`i'
                      * `b`i''^`kb'`i'
                      * `c`i''^`kc'`i'
                    #enddo
                #do i=1,{2*`n'}
                  )
                #enddo
                +
                #do i=1,`n'
                  sum_(`kb'`i',0,$`b'`i'-1,
                  #if `i' != `r'
                    sum_(`kc'`i',0,$`c'`i'-1,
                  #else
                    sum_(`kc'`i',$`c'`i',$`c'`i',
                  #endif
                #enddo
                    sign_(`kcsum')
                    * `kc'`r'
                    * fac_(`ksum'-1)
                    #do i=1,`n'
                      * invfac_(`kb'`i')
                      * invfac_(`kc'`i')
                    #enddo
                    * `poch'(-`ksum',$`E'+`kbsum'+1)
                    #do i=1,`n'
                      * `poch'(`kb'`i'+`kc'`i',$`a'`i')
                    #enddo
                    #do i=1,`n'
                      / `a`i''^`kb'`i'
                      / `a`i''^`kc'`i'
                      * `b`i''^`kb'`i'
                      * `c`i''^`kc'`i'
                    #enddo
                #do i=1,{2*`n'}
                  )
                #enddo
              #enddo
            ;
            #define n "`prefix'n1"
            #define x "`prefix'x1"
            id `poch'(0,`x'?) = 1;
            id `poch'(`n'?pos_,`x'?pos_) = fac_(`n'-1+`x') * invfac_(`x'-1);
            $`haspoch' = 1;
            #undefine n
            #undefine x
          endif;
    #do i=1,`n'
        endif;
      endif;
    #enddo
  endif;
  B `poch';
  ModuleOption local <$`a'1>,...,<$`a'`n'>,
                     <$`b'1>,...,<$`b'`n'>,
                     <$`c'1>,...,<$`c'`n'>,$`E';
  ModuleOption maximum $`haspoch';
#ifdef `FORCERSTATS'
  ModuleOption sum $`prefix'nrec, $`prefix'nsum;
#endif
  .sort:diamondS-`Z'-2;
#ifdef `FORCERSTATS'
  #message Stats: Diamond(`Z'): nrec=`$`prefix'nrec', nsum=`$`prefix'nsum'
#endif
  #if `$`haspoch''
    Keep Brackets;
    #define num "`prefix'num"
    #define den "`prefix'den"
    #define n "`prefix'n1"
    #define x "`prefix'x1"
    repeat id `poch'(`n'?pos_,`x'?) = `poch'(`n'-1,`x')*`num'(`n'-1+`x');
    repeat id `poch'(`n'?neg_,`x'?) = `poch'(`n'+1,`x')*`den'(`n'+`x');
    id `poch'(0,`x'?) = 1;
*   id `num'(`x'?)*`den'(`x'?) = 1;
*   id `num'(`x'?number_) = `x';
*   id `den'(`x'?number_) = 1/`x';
    id `num'(`x'?) = rat(`x',1);
    id `den'(`x'?) = rat(1,`x');
    .sort:diamondS-`Z'-3;
  #endif
#endprocedure

*--#] diamondS3 : 
*--#] diamondS : 
*--#[ insertionS :

**
* Performs massless one-loop insertion integrals:
*
*                     p1
*               ------<------
*              /   n1+e1*ep  \
*             /               \
*            /                 \
*     Q --<--                   --<-- Q
*            \                 /
*             \               /
*              \   n2+e2*ep  /
*               ------<------
*                     p2
*
* with possible dot products in the numerator (p1.Qj)^(-aj), ...
*
* Dot products containing p2 are NOT allowed.
*
* The trivial spherical factor per one-loop
*
*   1/(4 pi)^2
*
* is omitted as a convention.
*
* As the input, each factor is represented by a symbol like i1 := p1.p1
* multiplied by a label for the integration Z, so it looks like
*
*   Z / i1^(n1+e1*ep)
*     / i2^(n2+e2*ep)
*     / i3^(n3) ...
*
* In the output, Z is replaced to Zb. Many combination of dot products may
* appear, for which symbols are specified. In general, we do not assume any
* reuse of symbols for the input and output slots.
* 2 + n + 1 + n + n(n-1)/2 theta(n>0) = (n+2) + (n+1)(n+2)/2
* = (n+2)(n+3)/2 symbols are needed to be specified for the insertion integral
* with n dot products. Of course, the user may reuse the input symbols for some
* of the output symbols by specifying so.
*
* For example, with 3 dot products, 15 symbols are used as:
*
*        input  output
*        =====  ======
*   i1   p1.p1
*   i2   p2.p2
*
*   i3   p1.Q3
*   i4   p1.Q4
*   i5   p1.Q5
*
*   i6          Q.Q
*
*   i7          Q3.Q3
*   i8          Q4.Q4
*   i9          Q5.Q5
*
*   i10         Q.Q3
*   i11         Q.Q4
*   i12         Q.Q5
*
*   i13         Q3.Q4
*   i14         Q3.Q5
*   i15         Q4.Q5
*
* Here p1 is the loop momentum and p2 = Q - p1, Q3,Q4,Q5 are momenta outside
* the loop. (Some of) i10,...,i15 are expected to be rewritten in terms of
* propagators and numerators in the new topology.
*
* Parameters:
*   Z         - the symbol representing the integral
*   i1,i2,... - the symbols for the input and the output described above
*   e1,e2     - non-negative integers representing the non-integer parts of
*               the two lines depicted above
*
* e1 and e2 are expected be non-negative integers. No extra signs ("+") or no
* extra parentheses ("(" and ")") are allowed.
*
* Remarks:
*   The convention of passing arguments (in the ordering) has been changed from
*   the previous "function" version. The index "0" for the output must be
*   translated to "Q2" (or "1" without keeping the mass dimension).
*   There are no shift between the input and output due to removing "i2" in
*   the output.
*
#procedure insertionS(Z,i1,i2,?iv,e1,e2)
  #call Timing()

  #define prefix  "FORCER"
  #define LABEL   "`prefix'LABEL"

  #define n       "`prefix'n"
  #define q       "`prefix'q"

  #define n1      "`prefix'n1"
  #define n2      "`prefix'n2"
  #define n3      "`prefix'n3"
  #define isum    "`prefix'isum1"
  #define G       "`prefix'G"
  #define fnum    "`prefix'fnum"
  #define del     "`prefix'del"
  #define ftensor "`prefix'ftensor"
  #define distrib "`prefix'distrib"

  #define Q       "`q'1"

* Parse the variadic arguments.
  #define j "2"
  #do i={`?iv',}
    #ifdef `i'
      #redefine j "{`j'+1}"
      #define i`j' "`i'"
    #endif
  #enddo
  #define x "0"
  #do m=0,99
    #redefine x "{(`m'+2)*(`m'+3)/2}"
    #if `x' == `j'
      #define M "`m'"
      #breakdo
    #endif
    #if `x' > `j'
      #message `x' > `j'
      #breakdo
    #endif
  #enddo
  #if `x' != `j'
    #message Error: `j' symbols for insertionS.
    #terminate
  #endif
  #undefine j
  #undefine x

  if (count(`Z',1));
*   E.g., 1 / i1^1 / i2^2 * i3^3 * i4^4 * i5^5 (M = 3) is translated to
*         fnum(3,4,5,1,2,3+4+5).
    id ifnomatch->`LABEL'1, 1 / `i1'^`n1'? / `i2'^`n2'?
      #do j=3,{`M'+2}
        * `i`j''^`n'`j'?pos0_
      #enddo
      = `fnum'(
        #do j=3,{`M'+2}
          `n'`j',
        #enddo
          `n1',`n2',0
        #do j=3,{`M'+2}
          +`n'`j'
        #enddo
      );
*   E.g., fnum(3,4,5,1,2,3+4+5) becomes
*         fnum(1,2,3+4+5,q3,q3,q3,q4,q4,q4,q4,q5,q5,q5,q5,q5).
    #do j=3,{`M'+2}
      repeat id `fnum'(`n1'?pos_,?nn) = `fnum'(`n1'-1,?nn,`q'`j');
      id `fnum'(0,?nn) = `fnum'(?nn);
    #enddo
*   Now we have fnum(n1,n2,nsum,q3,...).
*   [dal(p)/4]^n -> fac_(n) / 2^n * distrib(1,2*n,dd,ftensor,...)
    id `fnum'(`n1'?,`n2'?,`n3'?,?pp) =
      sum_(`isum',0,integer_(`n3'/2),
           `G'(`n1',`e1',`n2',`e2',`n3',`isum')
           * `distrib'(1,2*`isum',`del',`ftensor',?pp)
           * `Q'.`Q'^`isum'
           / 2^`isum')
      * `Q'.`Q'^2
      / `Q'.`Q'^`n1'
      / `Q'.`Q'^`n2';
    id `distrib'(?a) = distrib_(?a);
    tovector `ftensor',`Q';
    id `del'(?a) = dd_(?a);
*   Replace Z by Zb.
    id `Z' = `Z'b;
*   Convert Q.Q^n to the first output symbol.
    id `Q'.`Q'^`n1'? = `i{`M'+3}'^`n1';
*   Convert all dot products qj1.qj2 in the numerator to corresponding output
*   symbols.
    #do j=3,{`M'+2}
      id `q'`j'.`q'`j'^`n1'? = `i{`j'+`M'+1}'^`n1';
    #enddo
    #do j=3,{`M'+2}
      id `Q'.`q'`j'^`n1'? = `i{`j'+2*`M'+1}'^`n1';
    #enddo
    #define j "{3*`M'+4}"
    #do j1=3,{`M'+2}
      #do j2={`j1'+1},{`M'+2}
        id `q'`j1'.`q'`j2' = `i`j'';
        #redefine j "{`j'+1}"
      #enddo
    #enddo
    #undefine j
    #if `SETQONE'
      id Q2^`n1'? = 1;
    #endif
  endif;
label `LABEL'1;
  B `G',GschemeConstants;
  .sort:insertionS-`Z'-1;
  Keep Brackets;
  #call reduceG()
  .sort:insertionS-`Z'-2;

  #call Timing(Insertion(`Z'))
#endprocedure

*--#] insertionS : 
*--#[ UseTurbo :

* This section controls the variable USETURBO, which is related to turbo rules
* for non-manual reduction rules (e.g., triangle, insertion, etc.).
*
* Our strategy for the turbo rules is now:
* - Start with enabling the turbo rules.
* - After reductions for topologies with 10 edges, they are disabled except
*   few topologies.
*
#define USETURBO "1"

#procedure SetUseTurbo(TOPO,VALUE)
* After one of them, the turbo rules are disabled except for them.
* d386 is the last topology with 10 edges.
* NOTE: these numberings may be changed in future.
  #do n={386,384,345,339,324,308,305}
    #if `TOPO' = d`n'
      #redefine USETURBO "`VALUE'"
      #breakdo
    #endif
  #enddo
#endprocedure

#procedure BeforeReduce(TOPO)
* Called just before Reduce`TOPO'.
  #call SetUseTurbo(`TOPO',1)
#endprocedure

#procedure AfterReduce(TOPO)
* Called just after Reduce`TOPO'.
  #call SetUseTurbo(`TOPO',0)
#endprocedure

*--#] UseTurbo : 
*--#[ counting :

#procedure count7(i,kA,kB,kC,kD,kE,kF,kG)
#do k = {A,B,C,D,E,F,G}
	if ( count(`k`k'',1) > $MaxCount`i'`k' ) $MaxCount`i'`k' = count_(`k`k'',1);
#enddo
#endprocedure

#procedure count6(i,kA,kB,kC,kD,kE,kF)
#do k = {A,B,C,D,E,F}
	if ( count(`k`k'',1) > $MaxCount`i'`k' ) $MaxCount`i'`k' = count_(`k`k'',1);
#enddo
#endprocedure

#procedure count5(i,kA,kB,kC,kD,kE)
#do k = {A,B,C,D,E}
	if ( count(`k`k'',1) > $MaxCount`i'`k' ) $MaxCount`i'`k' = count_(`k`k'',1);
#enddo
#endprocedure

#procedure count4(i,kA,kB,kC,kD)
#do k = {A,B,C,D}
	if ( count(`k`k'',1) > $MaxCount`i'`k' ) $MaxCount`i'`k' = count_(`k`k'',1);
#enddo
#endprocedure

*--#] counting : 
*--#[ filling :

#procedure filling(n,X,lhs,rhs)
 #if ( `$MaxCount`n'`X'' >= 7 )
  #do i = `$MaxCount`n'`X''-3,4,-4
	#$MaxSum`n' = $MaxSum`n'+1;
	Fill `PREFIX'dotmaptabl(`n',`$MaxSum`n'') = id_(`lhs'^`i',`lhs'^`i',`lhs',`rhs');
  #enddo
  #$MaxSum`n' = $MaxSum`n'+1;
  Fill `PREFIX'dotmaptabl(`n',`$MaxSum`n'') = id_(`lhs'^6,`lhs'^3*(`rhs')^3,`lhs'^5,`lhs'^3*(`rhs')^2,`lhs'^4,`lhs'^3*(`rhs'));
  #$MaxSum`n' = $MaxSum`n'+1;
  Fill `PREFIX'dotmaptabl(`n',`$MaxSum`n'') = id_(`lhs',`rhs');
 #elseif ( `$MaxCount`n'`X'' > 0 )
  #$MaxSum`n' = $MaxSum`n'+1;
  Fill `PREFIX'dotmaptabl(`n',`$MaxSum`n'') = id_(`lhs',`rhs');
 #endif
#endprocedure

*--#] filling : 
*--#] proc : 
#endif
* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4
