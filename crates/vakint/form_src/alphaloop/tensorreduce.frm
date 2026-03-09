*--#[ Tensor reduction :
#-

S dim,D,[D-4],ep,type,i;
Autodeclare S n,x;
Autodeclare I mu;
Autodeclare V Q,p,q;

CF vx,vxx,vxxs(s),vxxx,K,f,gg,gg1,TASK,dotacc,map,ex,ex1,exx,acc,acc1,acc2,loops,flip,DI,replace,
    cutv,fnum,fnum1,origct,f1,...,f7,vec,vec1,vec2,tmp,tmp1,g(s),accs(s),map,
    map1,map2,deg,REM,coeff,ext;

CF tt,tt1;
Set ns:n1,...,n32;
Set xs:x1,...,x22;
Set pp:p1,...,p32;
Set ppp:pp1,...,pp22;
Set pu:p20,...,p100;
Set Qu:Q20,...,Q100;
Set Qs:Q,Q1,...,Q10;
Set indices:mu1,...,mu100;

Auto V pzero;
Set pzeros : pzero1,...,pzero300;
Set ppzeros : ppzero1,...,ppzero300;
Set pzerospm : pzero1,...,pzero300,<-pzero1>,...,<-pzero300>;
Set Qpzeros : Q,pzero1,...,pzero300;

CF rat;


* The maximum rank of the preloaded tensor decomposition tables.
#ifndef `MAXRANK'
  #define MAXRANK "8"
#endif

* The maximum rank actually loaded.
#$MAXRANK = 0;

*--#] decl :
*--#[ pvgtab2 :

#$npvg2 = 1;

CTable pvgtab2(1:1,n1?,n2?);
Fill pvgtab2(1) = g(n1,n2);

*--#] pvgtab2 :
*--#[ pvctab2 :

CTable pvctab2(1:1);
Fill pvctab2(1) = rat(1,4-2*ep);

*--#] pvctab2 :
*--#[ pvgtab4 :

#$npvg4 = 2;

CTable pvgtab4(1:2,n1?,...,n4?);
Fill pvgtab4(1) = g(n1,n2)*g(n3,n4);
Fill pvgtab4(2) = g(n1,n3)*g(n2,n4)+g(n1,n4)*g(n2,n3);

*--#] pvgtab4 :
*--#[ pvctab4 :

CTable pvctab4(1:2);
Fill pvctab4(1) = rat(5-2*ep,72-108*ep+52*ep^2-8*ep^3);
Fill pvctab4(2) = rat(-1,72-108*ep+52*ep^2-8*ep^3);

*--#] pvctab4 :
*--#[ pvgtab6 :

#$npvg6 = 3;

CTable pvgtab6(1:3,n1?,...,n6?);
Fill pvgtab6(1) = g(n1,n2)*g(n3,n4)*g(n5,n6);
Fill pvgtab6(2) = g(n1,n2)*g(n3,n5)*g(n4,n6)+g(n1,n2)*g(n3,n6)*g(n4,n5)+
   g(n1,n3)*g(n2,n4)*g(n5,n6)+g(n1,n4)*g(n2,n3)*g(n5,n6)+g(n1,n5)*g(n2,
   n6)*g(n3,n4)+g(n1,n6)*g(n2,n5)*g(n3,n4);
Fill pvgtab6(3) = g(n1,n3)*g(n2,n5)*g(n4,n6)+g(n1,n3)*g(n2,n6)*g(n4,n5)+
   g(n1,n4)*g(n2,n5)*g(n3,n6)+g(n1,n4)*g(n2,n6)*g(n3,n5)+g(n1,n5)*g(n2,
   n3)*g(n4,n6)+g(n1,n5)*g(n2,n4)*g(n3,n6)+g(n1,n6)*g(n2,n3)*g(n4,n5)+g(
   n1,n6)*g(n2,n4)*g(n3,n5);

*--#] pvgtab6 :
*--#[ pvctab6 :

CTable pvctab6(1:3);
Fill pvctab6(1) = rat(26-22*ep+4*ep^2,1152-3168*ep+3280*ep^2-1600*ep^3+
   368*ep^4-32*ep^5);
Fill pvctab6(2) = rat(-1,192-464*ep+392*ep^2-136*ep^3+16*ep^4);
Fill pvctab6(3) = rat(2,1152-3168*ep+3280*ep^2-1600*ep^3+368*ep^4-32*ep^
   5);

*--#] pvctab6 :
*--#[ pvgtab8 :

#$npvg8 = 5;

CTable pvgtab8(1:5,n1?,...,n8?);
Fill pvgtab8(1) = g(n1,n2)*g(n3,n4)*g(n5,n6)*g(n7,n8);
Fill pvgtab8(2) = g(n1,n2)*g(n3,n4)*g(n5,n7)*g(n6,n8)+g(n1,n2)*g(n3,n4)*
   g(n5,n8)*g(n6,n7)+g(n1,n2)*g(n3,n5)*g(n4,n6)*g(n7,n8)+g(n1,n2)*g(n3,
   n6)*g(n4,n5)*g(n7,n8)+g(n1,n2)*g(n3,n7)*g(n4,n8)*g(n5,n6)+g(n1,n2)*g(
   n3,n8)*g(n4,n7)*g(n5,n6)+g(n1,n3)*g(n2,n4)*g(n5,n6)*g(n7,n8)+g(n1,n4)
   *g(n2,n3)*g(n5,n6)*g(n7,n8)+g(n1,n5)*g(n2,n6)*g(n3,n4)*g(n7,n8)+g(n1,
   n6)*g(n2,n5)*g(n3,n4)*g(n7,n8)+g(n1,n7)*g(n2,n8)*g(n3,n4)*g(n5,n6)+g(
   n1,n8)*g(n2,n7)*g(n3,n4)*g(n5,n6);
Fill pvgtab8(3) = g(n1,n3)*g(n2,n4)*g(n5,n7)*g(n6,n8)+g(n1,n3)*g(n2,n4)*
   g(n5,n8)*g(n6,n7)+g(n1,n4)*g(n2,n3)*g(n5,n7)*g(n6,n8)+g(n1,n4)*g(n2,
   n3)*g(n5,n8)*g(n6,n7)+g(n1,n5)*g(n2,n6)*g(n3,n7)*g(n4,n8)+g(n1,n5)*g(
   n2,n6)*g(n3,n8)*g(n4,n7)+g(n1,n6)*g(n2,n5)*g(n3,n7)*g(n4,n8)+g(n1,n6)
   *g(n2,n5)*g(n3,n8)*g(n4,n7)+g(n1,n7)*g(n2,n8)*g(n3,n5)*g(n4,n6)+g(n1,
   n7)*g(n2,n8)*g(n3,n6)*g(n4,n5)+g(n1,n8)*g(n2,n7)*g(n3,n5)*g(n4,n6)+g(
   n1,n8)*g(n2,n7)*g(n3,n6)*g(n4,n5);
Fill pvgtab8(4) = g(n1,n2)*g(n3,n5)*g(n4,n7)*g(n6,n8)+g(n1,n2)*g(n3,n5)*
   g(n4,n8)*g(n6,n7)+g(n1,n2)*g(n3,n6)*g(n4,n7)*g(n5,n8)+g(n1,n2)*g(n3,
   n6)*g(n4,n8)*g(n5,n7)+g(n1,n2)*g(n3,n7)*g(n4,n5)*g(n6,n8)+g(n1,n2)*g(
   n3,n7)*g(n4,n6)*g(n5,n8)+g(n1,n2)*g(n3,n8)*g(n4,n5)*g(n6,n7)+g(n1,n2)
   *g(n3,n8)*g(n4,n6)*g(n5,n7)+g(n1,n3)*g(n2,n5)*g(n4,n6)*g(n7,n8)+g(n1,
   n3)*g(n2,n6)*g(n4,n5)*g(n7,n8)+g(n1,n3)*g(n2,n7)*g(n4,n8)*g(n5,n6)+g(
   n1,n3)*g(n2,n8)*g(n4,n7)*g(n5,n6)+g(n1,n4)*g(n2,n5)*g(n3,n6)*g(n7,n8)
   +g(n1,n4)*g(n2,n6)*g(n3,n5)*g(n7,n8)+g(n1,n4)*g(n2,n7)*g(n3,n8)*g(n5,
   n6)+g(n1,n4)*g(n2,n8)*g(n3,n7)*g(n5,n6)+g(n1,n5)*g(n2,n3)*g(n4,n6)*g(
   n7,n8)+g(n1,n5)*g(n2,n4)*g(n3,n6)*g(n7,n8)+g(n1,n5)*g(n2,n7)*g(n3,n4)
   *g(n6,n8)+g(n1,n5)*g(n2,n8)*g(n3,n4)*g(n6,n7)+g(n1,n6)*g(n2,n3)*g(n4,
   n5)*g(n7,n8)+g(n1,n6)*g(n2,n4)*g(n3,n5)*g(n7,n8)+g(n1,n6)*g(n2,n7)*g(
   n3,n4)*g(n5,n8)+g(n1,n6)*g(n2,n8)*g(n3,n4)*g(n5,n7)+g(n1,n7)*g(n2,n3)
   *g(n4,n8)*g(n5,n6)+g(n1,n7)*g(n2,n4)*g(n3,n8)*g(n5,n6)+g(n1,n7)*g(n2,
   n5)*g(n3,n4)*g(n6,n8)+g(n1,n7)*g(n2,n6)*g(n3,n4)*g(n5,n8)+g(n1,n8)*g(
   n2,n3)*g(n4,n7)*g(n5,n6)+g(n1,n8)*g(n2,n4)*g(n3,n7)*g(n5,n6)+g(n1,n8)
   *g(n2,n5)*g(n3,n4)*g(n6,n7)+g(n1,n8)*g(n2,n6)*g(n3,n4)*g(n5,n7);
Fill pvgtab8(5) = g(n1,n3)*g(n2,n5)*g(n4,n7)*g(n6,n8)+g(n1,n3)*g(n2,n5)*
   g(n4,n8)*g(n6,n7)+g(n1,n3)*g(n2,n6)*g(n4,n7)*g(n5,n8)+g(n1,n3)*g(n2,
   n6)*g(n4,n8)*g(n5,n7)+g(n1,n3)*g(n2,n7)*g(n4,n5)*g(n6,n8)+g(n1,n3)*g(
   n2,n7)*g(n4,n6)*g(n5,n8)+g(n1,n3)*g(n2,n8)*g(n4,n5)*g(n6,n7)+g(n1,n3)
   *g(n2,n8)*g(n4,n6)*g(n5,n7)+g(n1,n4)*g(n2,n5)*g(n3,n7)*g(n6,n8)+g(n1,
   n4)*g(n2,n5)*g(n3,n8)*g(n6,n7)+g(n1,n4)*g(n2,n6)*g(n3,n7)*g(n5,n8)+g(
   n1,n4)*g(n2,n6)*g(n3,n8)*g(n5,n7)+g(n1,n4)*g(n2,n7)*g(n3,n5)*g(n6,n8)
   +g(n1,n4)*g(n2,n7)*g(n3,n6)*g(n5,n8)+g(n1,n4)*g(n2,n8)*g(n3,n5)*g(n6,
   n7)+g(n1,n4)*g(n2,n8)*g(n3,n6)*g(n5,n7)+g(n1,n5)*g(n2,n3)*g(n4,n7)*g(
   n6,n8)+g(n1,n5)*g(n2,n3)*g(n4,n8)*g(n6,n7)+g(n1,n5)*g(n2,n4)*g(n3,n7)
   *g(n6,n8)+g(n1,n5)*g(n2,n4)*g(n3,n8)*g(n6,n7)+g(n1,n5)*g(n2,n7)*g(n3,
   n6)*g(n4,n8)+g(n1,n5)*g(n2,n7)*g(n3,n8)*g(n4,n6)+g(n1,n5)*g(n2,n8)*g(
   n3,n6)*g(n4,n7)+g(n1,n5)*g(n2,n8)*g(n3,n7)*g(n4,n6)+g(n1,n6)*g(n2,n3)
   *g(n4,n7)*g(n5,n8)+g(n1,n6)*g(n2,n3)*g(n4,n8)*g(n5,n7)+g(n1,n6)*g(n2,
   n4)*g(n3,n7)*g(n5,n8)+g(n1,n6)*g(n2,n4)*g(n3,n8)*g(n5,n7)+g(n1,n6)*g(
   n2,n7)*g(n3,n5)*g(n4,n8)+g(n1,n6)*g(n2,n7)*g(n3,n8)*g(n4,n5)+g(n1,n6)
   *g(n2,n8)*g(n3,n5)*g(n4,n7)+g(n1,n6)*g(n2,n8)*g(n3,n7)*g(n4,n5)+g(n1,
   n7)*g(n2,n3)*g(n4,n5)*g(n6,n8)+g(n1,n7)*g(n2,n3)*g(n4,n6)*g(n5,n8)+g(
   n1,n7)*g(n2,n4)*g(n3,n5)*g(n6,n8)+g(n1,n7)*g(n2,n4)*g(n3,n6)*g(n5,n8)
   +g(n1,n7)*g(n2,n5)*g(n3,n6)*g(n4,n8)+g(n1,n7)*g(n2,n5)*g(n3,n8)*g(n4,
   n6)+g(n1,n7)*g(n2,n6)*g(n3,n5)*g(n4,n8)+g(n1,n7)*g(n2,n6)*g(n3,n8)*g(
   n4,n5)+g(n1,n8)*g(n2,n3)*g(n4,n5)*g(n6,n7)+g(n1,n8)*g(n2,n3)*g(n4,n6)
   *g(n5,n7)+g(n1,n8)*g(n2,n4)*g(n3,n5)*g(n6,n7)+g(n1,n8)*g(n2,n4)*g(n3,
   n6)*g(n5,n7)+g(n1,n8)*g(n2,n5)*g(n3,n6)*g(n4,n7)+g(n1,n8)*g(n2,n5)*g(
   n3,n7)*g(n4,n6)+g(n1,n8)*g(n2,n6)*g(n3,n5)*g(n4,n7)+g(n1,n8)*g(n2,n6)
   *g(n3,n7)*g(n4,n5);

*--#] pvgtab8 :
*--#[ pvctab8 :

CTable pvctab8(1:5);
Fill pvctab8(1) = rat(287-278*ep+84*ep^2-8*ep^3,28800-125280*ep+199504*
   ep^2-159680*ep^3+71152*ep^4-17888*ep^5+2368*ep^6-128*ep^7);
Fill pvctab8(2) = rat(-166+198*ep-72*ep^2+8*ep^3,57600-308160*ep+649568*
   ep^2-718368*ep^3+461664*ep^4-178080*ep^5+40512*ep^6-4992*ep^7+256*ep^
   8);
Fill pvctab8(3) = rat(54-26*ep+4*ep^2,57600-308160*ep+649568*ep^2-
   718368*ep^3+461664*ep^4-178080*ep^5+40512*ep^6-4992*ep^7+256*ep^8);
Fill pvctab8(4) = rat(2,1800-8280*ep+13864*ep^2-11016*ep^3+4432*ep^4-
   864*ep^5+64*ep^6);
Fill pvctab8(5) = rat(-26+10*ep,57600-308160*ep+649568*ep^2-718368*ep^3+
   461664*ep^4-178080*ep^5+40512*ep^6-4992*ep^7+256*ep^8);

*--#] pvctab8 :
*--#[ LoadPVTable :

#procedure LoadPVTable(rank)
  #do r=2,`rank',2
    #if `r' > `$MAXRANK'
      #if `r' >= 10
        #include pvtab`r'.h
      #endif
      CTable dd`r'tab(n1?,...,n`r'?);
      Fill dd`r'tab =
        #do i=1,`$npvg`r''
          + pvgtab`r'(`i',n1,...,n`r')
        #enddo
      ;
      #$MAXRANK = `r';
    #endif
  #enddo
#endprocedure

*--#] LoadPVTable :
*--#[ init :

#call LoadPVTable(`MAXRANK')

#ifdef `TENSORSYMTABLE'
  #include pvtab8sym.h
#endif

* The maximum rank actually we have encountered.
#$maxrank = 0;

S xpv;
CF g1;

* Expected input: all vectors that belong to the subgraph that needs to be projected written
* as vec(p1,1)*vec(p2,2)... All external vectors can be written as vec1(p3,1)*vec1(p4,2),...
#procedure TensorReduce()
    B+ vec,g;
    .sort:tensor-start;
    Keep brackets;

    repeat id g(n1?,n2?)*g(n2?,n3?) = g(n1,n3);

*   Contract dot products internal to the subgraph
    id vec(p1?,n?)*vec(p2?,n?) = p1.p2;

* note that the following operation is not safe, since we always consider the metric to be on the outside
*    repeat id vec(p1?,n?)*g(n?,n1?) = vec(p1,n1);

    id vec(p?,n?) = vec(p,n)*tt(vec(p,n));
    chainin tt;

*   Find the rank.
    id tt(?a) = tt(?a) * xpv^nargs_(?a);
    id xpv^n?odd_ = 0;

    if (count(xpv,1) > $maxrank) $maxrank = count_(xpv,1);

    #define r "2"
    #do loop=1,1
      if (count(xpv,1) == `r');
        id xpv = 1;

          Multiply acc(1);
          repeat id acc(x1?)*f?{tt,vec,vec1}(?a) = acc(x1*f(?a));
          id acc(x?$expr) = 1;

          inside $expr;
          id tt(?a) = tt(?a)*acc1(?a);
          argument acc1;
            id vec(p?,n?) = p;
          endargument;
          transform tt mulargs(1,last);

*         Generate all contracted structures that have a different
*         coefficient under the internal symmetry
          id acc1(?a) = partitions_(0,gg,2,?a);
          DropCoefficient;

*         generate representative for projection operator
          id gg(p1?,p2?)*vec(p1?,n1?)*vec(p2?,n2?) = g1(n1,n2)*gg1(p1,p2)*vec(p1,n1)*vec(p2,n2);
          id gg1(?a) = gg(?a);
          chainin g1;

          #if (isdefined(TENSORSYMTABLE)) && (`r' >= 8)
*         construct pattern for internal symmetry
            id g1(?a) = g1(?a)*acc1(?a);
            repeat id acc1(?a,n?,?b)*vec(p?,n?) = acc1(?a,p,?b)*vec(p,n);
            id acc1(?a) = acc1(1,?a)*acc2(?a);
            transform acc1 dedup(1,last);
            repeat id acc1(n?,p?,?a) = map(p,n)*acc1(n+1,?a);
            repeat id acc2(?a,p?,?b)*map(p?,n?) = acc2(?a,n,?b)*map(p,n);
            id map(?a) = 1;

            id g1(?a)*tt(x?)*acc1(n?)*acc2(?b) = acc(g1(?a)*tt(x)*acc1(n)*acc2(?b));
            id acc(x?$proj) = 1;

            inside $proj;
               id g1(?a)*tt(x?)*acc1(n?)*acc2(?b) = x*pvgtab`r'(1,?a) * pvctab`r'(1) + sum_(i,2,`$npvg`r'', x*pvgstab`r'(i,?b,?a) * pvctab`r'(i));
            endinside;
          #else
            id g1(?a)*tt(x?) = acc(g1(?a)*tt(x));
            id acc(x?$proj) = 1;
            inside $proj;
              id g1(?a)*tt(x?) = sum_(i,1,`$npvg`r'', x*pvgtab`r'(i,?a) * pvctab`r'(i));
            endinside;
          #endif

            inside $proj;
*             if not found in symmetric table, use old table
              #if (isdefined(TENSORSYMTABLE)) && (`r' >= 8)
                id pvgstab`r'(n?,n1?,...,n`r'?,?a) = pvgtab`r'(n,?a);
              #endif
              repeat id vec(p?,n1?)*g(n1?,n2?) = vec(p,n2);
              id vec(p1?,n1?)*vec(p2?,n1?) = p1.p2;
            endinside;

*         now create the full tensor structure, while considering the outside symmetry
            repeat id vec(p?,n?,?a)*vec(p?,n1?,?b) = vec(p,n,?a,n1,?b);

* now merge in the externals, which could be in metrics or vec1(Q,..)
            repeat id vec(p?,?a,n?,?b)*vec1(Q?,n?) = vec(p,?a,Q,?b)*vec1; * mark that we have an external Q

            Multiply tmp(1);
            repeat id g(n1?,n2?)*tmp(n?) = vxx(ppp[n],n1)*vxx(ppp[n],n2)*tmp(n+1);
            repeat id vxx(-p?vector_,n?)*tmp(n1?) = vxx(ppp[n1],n)*map(ppp[n1],-p)*tmp(n1+1); * rename -vertex to positive one
            repeat id vec(p?,?a,n?,?b)*vxx(p1?,n?) = vec(p,?a,p1,?b);

* any remaining indices are Taylor indices or in a Gtrace
            repeat id vec(?a,n?,?b)*tmp(n1?) = vec(?a,ppp[n1],?b)*map(ppp[n1],n)*tmp(n1+1);
            id tmp(n?) = 1;

            #do i=2,15
               symmetrize vec 2,...,`i'; * symmetrize input
            #enddo

* generate without duplicates the complete g(mu,nu) structure for a given contraction
* while taking into account symmetries on the outside of the projected structure
* this could be considered the inverse of dd_
            id gg(p1?,p1?) = gg1(1,p1);
            repeat id gg1(n1?,p?)*gg1(n2?,p?) = gg1(n1+n2,p); * treat p1.p1^n
            repeat;
              id ifnomatch->endsubs1`r' once gg1(n1?,p?)*vec(p?,?a) = distrib_(1,n1*2,f1,f2,?a)*gg(p);
              id f2(?a)*gg(p?) = vec(p,?a);
              id vec(p?) = 1;
              id f1(?a) = dd_(?a);
            endrepeat;
            label endsubs1`r';

            id gg(p1?,p2?) = gg1(1,p1,p2); * treat p1.p2^n
            repeat id gg1(n1?,p1?,p2?)*gg1(n2?,p1?,p2?) = gg1(n1+n2,p1,p2);

            repeat;
              id ifnomatch->endsubs2`r' once gg1(n1?,p1?,p2?)*vec(p1?,?a)*vec(p2?,?b) =
                  distrib_(1,n1,f1,f2,?a)*distrib_(1,n1,f3,f4,?b)*gg1(p1,p2);
              id f2(?a)*f4(?b)*gg1(p1?,p2?) = vec(p1,?a)*vec(p2,?b);
* go through all permutations
               repeat;
                id once f3(p?,?b) = distrib_(1,1,f4,f3,p,?b);
                id f1(p1?,?a)*f4(p2?) = p1.p2*f1(?a);
               endrepeat;
              id f?{f1,f3} = 1;
            endrepeat;
            label endsubs2`r';
            id vec(?a) = 1;

            repeat id p1?.p2?ppp*p3?.p2?ppp = p1.p3;
            id p?ppp.p?ppp = rat(4-2*ep,1);
            id p1?.p?*map(p?,n?) = vxx(p1,n);
            id p1?.p2? = vxx(p1,p2);
            id vxx(p?,n?)*map(p?,n1?) = g(n,n1);
            repeat id vxx(?a,p?,?b)*map(p?,p1?) = vxx(?a,p1,?b);
            if (count(vec1,1)); * TODO: why this?
              id vxx(p1?,p2?) = p1.p2; * NEW
              id vec1 = 1;
            endif;

            repeat id g(n1?,n2?)*g(n2?,n3?) = g(n1,n3);
            id g(n1?,n1?) = rat(4-2*ep,1);
            id vxx(p?,n?)*g(n?,n1?) = vxx(p,n1);
            id vxx(p?,n?)*vxx(p1?,n?) = vxx(p,p1);

* TODO: is this correct?
            id vxx(p1?,p2?) = p1.p2;
            id vxx(p1?,n?) = vec1(p1,n);

            if (count(map,1)); Print "Mapping failed %t"; exit ""; endif;
            Multiply $proj;
          endinside;
          Multiply $expr;
      endif;

    B+ vec,g,tt,vec1,xpv,vxx;
    ModuleOption local,$expr,$proj;
    ModuleOption maximum,$maxrank;
    .sort:tensor-`r';
    #if `$maxrank' > `$MAXRANK'
      #message Load PV tables up to r=`$maxrank'
      #call LoadPVTable(`$maxrank')

      #if (isdefined(TENSORSYMTABLE)) && (`$maxrank' == 10)
        #message Load rank 10 PV symmetry table
        #include pvtab10sym.h
      #endif
    #endif

    #if `r' < `$MAXRANK'
      #define loop "0"
      #define r "{`r'+2}"
      Keep brackets;
    #endif
  #enddo
  #undefine r

  if (count(xpv,1)); Print "Not all tensors reduced: %t"; exit ""; endif;
#endprocedure

* Example input:
*Polyratfun rat;
*L F = vec(p1,1)*vec(p2,2)*vec(p3,3)*vec(p4,4)*vec(p8,5)*vec(p8,6)*vec1(p5,1)*vec1(p5,2)*vec1(Q,3)*vec1(Q,4)*vec1(p9,5)*vec1(p10,6);
*L F = vec(p1,1)*vec(p2,2)*vec1(p3,1);
*#call TensorReduce()
