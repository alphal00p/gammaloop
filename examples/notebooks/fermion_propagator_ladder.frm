#-
* Massless fermion-ring propagator ladder, with scalar external projection.
* LOOPS=3: ring v1->v2->v3->v4->v5->v6->v1; external q enters v1,
* leaves v4. Vector rungs v6->v2 and v5->v3 carry k2 and k3.
* Clockwise fermion momenta p1,...,p6 are
* (k1,k1+k2,k1+k2+k3,k1+k2+k3-q,k1+k2-q,k1-q).
* Thus E=8,V=6,L=E-V+1=3. All vertices conserve momentum.
* LOOPS=4: extend to v1,...,v8 with q leaving v5 and rungs
* v8->v2, v7->v3, v6->v4 carrying k2,k3,k4. Then E=11,V=8,L=4.
* Repeated vertex indices sew the Feynman-gauge vector numerators and
* external g(mu,nu). Includes the closed-fermion-loop minus; omits couplings,
* color, propagator denominators and overall i factors. No on-shell limits.
*
* form -q -d MODE=trace4 -d OUTPUT=result4.txt fermion_propagator_ladder.frm
* form -q -d MODE=tracen -d OUTPUT=resultD.txt fermion_propagator_ladder.frm
* form -q -d MODE=tracen -d DIAGNOSTIC=0 -d REPEATS=32 -d BATCHES=3
*      -d ORDER=trace-first -d LOOPS=3 fermion_propagator_ladder.frm
*
* MODE=trace4 uses dimension 4; MODE=tracen uses symbolic D. tr(1)=4 in both.
* ORDER=trace-first (default) traces before substituting routed momenta;
* ORDER=routing-first substitutes routed momenta first.
* ROUTEGROUP=1 (default) sorts after each momentum substitution. Set 2 to
* sort after every pair of substitutions, or 0 to defer routing collection.
* COLLECT=1 (default) adds a sort between the trace and routing stages;
* COLLECT=0 omits that boundary sort during timing. Routing-group sorts still
* apply. Diagnostic stage counts add that sort even with COLLECT=0, so their
* execution cost is not a timing substitute.
* DIAGNOSTIC=1 (default) enables FORM statistics, reports stages and optionally
* writes OUTPUT. Statistics remain off in DIAGNOSTIC=0 timing runs.
* DIAGNOSTIC=0 times fresh independent expressions in each batch, without
* intermediate reporting/export.
* Input declaration/initial sort are outside the body clock; trace, routing
* expansion, contraction and collection are inside. Output counts, equality
* checks and disposal are outside. No batch is silently discarded.
* BATCH_CPU_MS / REPEATS is amortized body CPU milliseconds, not wall time.
* Whole-process wall time additionally includes setup, checking and disposal.

#ifndef `LOOPS'
  #define LOOPS "3"
#endif
#ifndef `MODE'
  #define MODE "trace4"
#endif
#ifndef `ORDER'
  #define ORDER "trace-first"
#endif
#ifndef `ROUTEGROUP'
  #define ROUTEGROUP "1"
#endif
#ifndef `COLLECT'
  #define COLLECT "1"
#endif
#ifndef `DIAGNOSTIC'
  #define DIAGNOSTIC "1"
#endif
#ifndef `REPEATS'
  #if `DIAGNOSTIC' == 1
    #define REPEATS "1"
  #else
    #define REPEATS "32"
  #endif
#endif
#ifndef `BATCHES'
  #if `DIAGNOSTIC' == 1
    #define BATCHES "1"
  #else
    #define BATCHES "3"
  #endif
#endif
#if `REPEATS' < 1
  #message REPEATS and BATCHES must be positive
  #terminate
#endif
#if `BATCHES' < 1
  #message REPEATS and BATCHES must be positive
  #terminate
#endif
#if `DIAGNOSTIC' < 0
  #message DIAGNOSTIC must be 0 or 1
  #terminate
#endif
#if `DIAGNOSTIC' > 1
  #message DIAGNOSTIC must be 0 or 1
  #terminate
#endif
#if `COLLECT' < 0
  #message COLLECT must be 0 or 1
  #terminate
#endif
#if `COLLECT' > 1
  #message COLLECT must be 0 or 1
  #terminate
#endif

#if `ROUTEGROUP' < 0
  #message ROUTEGROUP must be 0, 1 or 2
  #terminate
#endif
#if `ROUTEGROUP' > 2
  #message ROUTEGROUP must be 0, 1 or 2
  #terminate
#endif

#if `DIAGNOSTIC' == 1
  On Statistics;
#else
  Off Statistics;
#endif
Format nospaces;
Symbols D;
#if "`MODE'" == "trace4"
  Dimension 4;
#elseif "`MODE'" == "tracen"
  Dimension D;
#else
  #message MODE must be trace4 or tracen
  #terminate
#endif
Indices mu,a,c,e;
Vectors p1,...,p8,k1,...,k4,q;
UnitTrace 4;

#procedure Route
  id p1 = k1;
  #if `ROUTEGROUP' == 1
    .sort:route-p1;
  #endif
  id p2 = k1+k2;
  #if `ROUTEGROUP' == 1
    .sort:route-p2;
  #elseif `ROUTEGROUP' == 2
    .sort:route-p2;
  #endif
  id p3 = k1+k2+k3;
  #if `ROUTEGROUP' == 1
    .sort:route-p3;
  #endif
  #if `LOOPS' == 3
    id p4 = k1+k2+k3-q;
    #if `ROUTEGROUP' == 1
      .sort:route-p4;
    #elseif `ROUTEGROUP' == 2
      .sort:route-p4;
    #endif
    id p5 = k1+k2-q;
    #if `ROUTEGROUP' == 1
      .sort:route-p5;
    #endif
    id p6 = k1-q;
    #if `ROUTEGROUP' == 1
      .sort:route-p6;
    #elseif `ROUTEGROUP' == 2
      .sort:route-p6;
    #endif
  #elseif `LOOPS' == 4
    id p4 = k1+k2+k3+k4;
    #if `ROUTEGROUP' == 1
      .sort:route-p4;
    #elseif `ROUTEGROUP' == 2
      .sort:route-p4;
    #endif
    id p5 = k1+k2+k3+k4-q;
    #if `ROUTEGROUP' == 1
      .sort:route-p5;
    #endif
    id p6 = k1+k2+k3-q;
    #if `ROUTEGROUP' == 1
      .sort:route-p6;
    #elseif `ROUTEGROUP' == 2
      .sort:route-p6;
    #endif
    id p7 = k1+k2-q;
    #if `ROUTEGROUP' == 1
      .sort:route-p7;
    #endif
    id p8 = k1-q;
    #if `ROUTEGROUP' == 1
      .sort:route-p8;
    #elseif `ROUTEGROUP' == 2
      .sort:route-p8;
    #endif
  #else
    #message LOOPS must be 3 or 4
    #terminate
  #endif
#endprocedure

#do batch=1,`BATCHES'
  #do j=1,`REPEATS'
    #if `LOOPS' == 3
      Local F`j' = -g_(1,mu,p1,a,p2,c,p3,mu,p4,c,p5,a,p6);
    #else
      Local F`j' = -g_(1,mu,p1,a,p2,c,p3,e,p4,mu,p5,e,p6,c,p7,a,p8);
    #endif
  #enddo
  .sort
  #if `DIAGNOSTIC' == 0
    #reset timer
  #else
    #$terms = termsin_(F1);
    #write "STAGE=word TERMS=%$",$terms
  #endif
  #if "`ORDER'" == "routing-first"
    #call Route
    #if `COLLECT' == 1
      .sort:routed;
    #elseif `DIAGNOSTIC' == 1
      .sort:routed;
    #endif
    #if `DIAGNOSTIC' == 1
      #$terms = termsin_(F1);
      #write "STAGE=routed TERMS=%$",$terms
    #endif
  #elseif "`ORDER'" != "trace-first"
    #message ORDER must be trace-first or routing-first
    #terminate
  #endif
  `MODE',1;
  #if "`ORDER'" == "trace-first"
    #if `COLLECT' == 1
      .sort:traced;
    #elseif `DIAGNOSTIC' == 1
      .sort:traced;
    #endif
    #if `DIAGNOSTIC' == 1
      #$terms = termsin_(F1);
      #write "STAGE=traced TERMS=%$",$terms
    #endif
    #call Route
  #endif
  .sort:scalar;
  #if `DIAGNOSTIC' == 0
    #write "BATCH=`batch' REPEATS=`REPEATS' BATCH_CPU_MS=`timer_'"
  #endif
  #$terms = termsin_(F1);
  #write "ROUTEGROUP=`ROUTEGROUP'"
  #write "MODE=`MODE' LOOPS=`LOOPS' ORDER=`ORDER' COLLECT=`COLLECT' STAGE=scalar TERMS=%$",$terms
  #if `DIAGNOSTIC' == 1
    #ifdef `OUTPUT'
      #write <`OUTPUT'> "%E",F1
    #endif
  #endif
  #do j=2,`REPEATS'
    Local Check`j' = F`j'-F1;
  #enddo
  .sort
  #do j=2,`REPEATS'
    #$difference = termsin_(Check`j');
    #if `$difference' != 0
      #message Independent batch expressions disagree
      #terminate
    #endif
  #enddo
  #write "BATCH=`batch' EQUAL_EXPRESSIONS=`REPEATS'"
  Drop;
  .sort
#enddo
.end
