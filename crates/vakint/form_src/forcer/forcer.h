#ifndef `FORCERHFILE'
#define FORCERHFILE

**
* The input can consist of the following symbols in a proper format.
* Any use of other symbols is at your own risk.
*
CF vx,ex;
V p1,...,p11,Q;

**
* The output will be given in terms of the following symbols.
*
S ep;
CF rat,GschemeConstants,Master;
S haha,no1,no2,no3,no6,lala,nono,bebe,cross,fastar2,nostar5,nostar6,
  no,t105,t124,t145,t155;
S z2,z3,z4,z5,z6,z7,z8,z5z3,z9,z10,z7z3,z11,z5z3z3,z12,z9z3,z6z4z1z1;

**
* The internals.
*
#appendpath forcer
#define FORCERH
#include pathforcer.h
#include recognition.h
#include forcer-manual.h # FORCERD
#include forcer-manual.h # FromSpectator
#do i=1,14
  #define POS`i' "pos_"
  #define NEG`i' "neg0_"
#enddo
#do i=1,9
  #define SKIP`i' ""
#enddo
S G10,G20,G30,G40,G50;
CF frat;
AutoDeclare S ALARM;

S order,xxx1,xxx2;

*--#[ Forcer : Forcer([options...])

**
* The main user interface of Forcer.
*
* Options
* -------
* gschemeexpand=off
* msbarexpand=off
*   Do not expand the master integrals in ep (default).
*
* gschemeexpand=on
* gschemeexpand[=d]: d is int
*   Expand the result in the G-scheme up to depth d (default: d=4, i.e.,
*   up to ep^0 for 4-loops).
*
* msbarexpand=on
* msbarexpand[=d]: d is int
*   Expand the result in the MSbar scheme up to depth d (default: d=4).
*
* keepzeta2[=on]
* keepzeta2=off
*   Keep z2 for the expansions in the MSbar scheme (default: off).
*
* polyratfunexpand=off
*   Do not expand rational coefficients at intermediate steps (default).
*
* polyratfunexpand=div
*   Keep only the lowest power of ep in coefficients at intermediate steps.
*
* polyratfunexpand=n: n is int
*   Expand coefficients with respect to ep up to ep^n at intermediate steps.
*   Note that the result is not guaranteed to be correct up to ep^n due to
*   (spurious) poles. The user should experimentally determine the depth that
*   gives the correct result.
*
#procedure Forcer(?options)
* NOTE: this procedure overwrites MSBAR, KEEPZETA2 and POLYRATFUNEXPAND.

* Parse the arguments.
  #$ForcerExpandScheme = 0;
  #$ForcerExpandDepth = 0;
  #$ForcerPolyRatFunExpand = -1;

  #define KEEPZETA2 ""

  #define gschemeexpand(n) "
    #case gschemeexpand=`~n'
      #$ForcerExpandScheme = 1;
      #$ForcerExpandDepth = `~n';
      #break
  "

  #define msbarexpand(n) "
    #case msbarexpand=`~n'
      #$ForcerExpandScheme = 2;
      #$ForcerExpandDepth = `~n';
      #break
  "

  #define polyratfunexpand(n) "
    #case polyratfunexpand=`~n'
      #$ForcerPolyRatFunExpand = `~n';
      #break
  "

  #do a={`?options',}
    #ifdef `a'
      #switch `a'
        #case gschemeexpand=off
        #case msbarexpand=off
          #$ForcerExpandScheme = 0;
          #$ForcerExpandDepth = 0;
          #break
        #case gschemeexpand
        #case gschemeexpand=on
          #$ForcerExpandScheme = 1;
          #$ForcerExpandDepth = 4;
          #break
        `gschemeexpand(0)'
        `gschemeexpand(1)'
        `gschemeexpand(2)'
        `gschemeexpand(3)'
        `gschemeexpand(4)'
        `gschemeexpand(5)'
        `gschemeexpand(6)'
        `gschemeexpand(7)'
        `gschemeexpand(8)'
        `gschemeexpand(9)'
        #case msbarexpand
        #case msbarexpand=on
          #$ForcerExpandScheme = 2;
          #$ForcerExpandDepth = 4;
          #break
        `msbarexpand(0)'
        `msbarexpand(1)'
        `msbarexpand(2)'
        `msbarexpand(3)'
        `msbarexpand(4)'
        `msbarexpand(5)'
        `msbarexpand(6)'
        `msbarexpand(7)'
        `msbarexpand(8)'
        `msbarexpand(9)'
        #case keepzeta2
        #case keepzeta2=on
          #redefine KEEPZETA2
          #break
        #case keepzeta2=off
          #redefine KEEPZETA2 ""
          #break
        #case polyratfunexpand=off
          #$ForcerPolyRatFunExpand = -1;
          #break
        #case polyratfunexpand=div
          #$ForcerPolyRatFunExpand = 0;
          #break
        `polyratfunexpand(1)'
        `polyratfunexpand(2)'
        `polyratfunexpand(3)'
        `polyratfunexpand(4)'
        `polyratfunexpand(5)'
        `polyratfunexpand(6)'
        `polyratfunexpand(7)'
        `polyratfunexpand(8)'
        `polyratfunexpand(9)'
        `polyratfunexpand(10)'
        `polyratfunexpand(11)'
        `polyratfunexpand(12)'
        `polyratfunexpand(13)'
        `polyratfunexpand(14)'
        `polyratfunexpand(15)'
        `polyratfunexpand(16)'
        `polyratfunexpand(17)'
        `polyratfunexpand(18)'
        `polyratfunexpand(19)'
        `polyratfunexpand(20)'
        `polyratfunexpand(21)'
        `polyratfunexpand(22)'
        `polyratfunexpand(23)'
        `polyratfunexpand(24)'
        `polyratfunexpand(25)'
        `polyratfunexpand(26)'
        `polyratfunexpand(27)'
        `polyratfunexpand(28)'
        `polyratfunexpand(29)'
        `polyratfunexpand(30)'
        `polyratfunexpand(31)'
        `polyratfunexpand(32)'
        `polyratfunexpand(33)'
        `polyratfunexpand(34)'
        `polyratfunexpand(35)'
        `polyratfunexpand(36)'
        `polyratfunexpand(37)'
        `polyratfunexpand(38)'
        `polyratfunexpand(39)'
        `polyratfunexpand(40)'
        `polyratfunexpand(41)'
        `polyratfunexpand(42)'
        `polyratfunexpand(43)'
        `polyratfunexpand(44)'
        `polyratfunexpand(45)'
        `polyratfunexpand(46)'
        `polyratfunexpand(47)'
        `polyratfunexpand(48)'
        `polyratfunexpand(49)'
        `polyratfunexpand(50)'
        #default
          #message Error: Forcer(`?options'): unknown option `a'
          #terminate
          #break
      #endswitch
    #endif
  #enddo

  #ifdef `KEEPZETA2'
    #if `$ForcerExpandScheme' != 2
      #message Warning: Forcer(`?options'): keepzeta2 affects only in MSbar
    #endif
  #endif

  #define POLYRATFUNEXPAND "`$ForcerPolyRatFunExpand'"

  .sort:Forcer-start;

* Parform the reduction for each expression.
  #if `numactiveexprs_' == 0
  #elseif `numactiveexprs_' == 1
    #define NAME "`activeexprnames_'"
    #call ForcerImpl()
  #else
    #define ForcerExprNames "`activeexprnames_'"
    Hide `ForcerExprNames';
    .sort:Forcer-hide;
    #do NAME={`ForcerExprNames'}
      #ifdef `NAME'
        Unhide `NAME';
        #call ForcerImpl()
        .sort:Forcer-end-`NAME';
        Hide `NAME';
      #endif
    #enddo
    .sort:Forcer-unhide;
    Unhide `ForcerExprNames';
  #endif
#endprocedure

#define ForcerNotationLoaded "0"

#procedure ForcerImpl()
  #if `$ForcerExpandScheme' >= 1
*   Count the number of loops.
    id order = 1;
    id vx(?a) = vx(?a)*xxx1*xxx2^nargs_(?a);
    id xxx1^n1?*xxx2^n2? = order^(n2 / 2 - n1);
    id ex(?a) = ex(?a)*order^nargs_(?a);
  #endif

* Match the topology.
  #if `ForcerNotationLoaded' == 0
    .sort:Forcer-notation;
    #call loadTopologies()
    #redefine ForcerNotationLoaded "1"
  #endif
  #call matchTopologies(1)

* Convert dot products into the Forcer notation.
  #call ToForcer()

* Compute the integrals.
  #call DoForcer()

* Convert the output into the standard form.
* #call propagators(1)
  id BasicNOIntegral = Master(no);
  id BasicT1Integral = Master(t105)*rat(1,ep^2);
  id BasicNOIntegral5 = Master(nostar5);
  id BasicNOIntegral6 = Master(nostar6);
  id BasicT1Integral55 = Master(t155);
  id BasicT1Integral45 = Master(t145);
  id BasicT1Integral24 = Master(t124);

* GC(0,0) GC(n,1) = GC(n,0) GC(n+1,0) (1-(n+2)ep) / (1-2ep)
  id GschemeConstants(0,0)*GschemeConstants(2,1)
    = GschemeConstants(2,0)*GschemeConstants(3,0)
      *rat(1-4*ep,1-2*ep);
  id GschemeConstants(0,0)*GschemeConstants(1,1)
    = GschemeConstants(1,0)*GschemeConstants(2,0)
      *rat(1-3*ep,1-2*ep);

  #if `$ForcerExpandScheme' >= 1
*   Expand the master inegrals.
    #if `$ForcerExpandScheme' == 1
      #define SchemeName "GScheme"
    #elseif `$ForcerExpandScheme' == 2
      #define SchemeName "MSBAR"
    #endif
    #do i=1,5
      id GschemeConstants(`i',0)
        = G`i'0*GschemeConstants(0,0)*rat(1-2*ep,1-{`i'+2}*ep);
    #enddo
    id GschemeConstants(0,0) = rat(1,ep);
    .sort:To`SchemeName'-1;
    PolyRatFun rat(expand,ep,12);
    .sort:To`SchemeName'-2;
    PolyRatFun;
    id rat(x?) = x;
    .sort:To`SchemeName'-3;
    #if `$ForcerExpandScheme' == 2
      #define MSBAR
    #else
      #define MSBAR ""
    #endif
    #call Masters(order,`$ForcerExpandDepth')
    id order = 1;
    .sort:To`SchemeName'-4;
  #endif
#endprocedure

*--#] Forcer :
#endif
* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4
