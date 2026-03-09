* 
*   MATAD-ng - exact version of M.Steinhauser MATAD package
*   rat() function used at intermediate steps and intermediate
*   results do not expanded in ep.
* 
*   For details see: https://github.com/apik/matad-ng
* 
*   Andrey Pikelner, pikelner[at]theor.jinr.ru
* 
**************************************************************************************
* 
*   Options and flags
*   -------------
* 
*   REDBNTAB   - use table for integrals 
*                BN(1,1,1,1), BN(1,1,1,2), BN(1,1,2,2), 
*                BN(1,2,2,2), BN(2,2,2,2) and dalambertian 
*                application upto dala^11
* 
*   Integration routines
*   --------------------
* 
*   IntOne     - one-loop massless self-energy
*   TadpoleM0  - one-loop massive tadpole
*   TadpoleMM0 - two-loop tadpole with two massive and one massless line
* 
* 
*   Conventions used
*   ----------------
* 
*   All integrals are expressed in terms of master integrals:
* 
*       (1) - miT111 at two-loop and
*       (9) - miD6,miD5,miD4,miD3,miDN,miDM,miE3,miBN,miBN1 at three-loop level
* 
*   Final answer contain symbol *d* instead of *n* in original MATAD package
* 
*   When expansion used all master integrals divided by Exp(-ep*EulerGamma) 
*   for each loop. Gamma functions present in final result expanded according to:
* 
*       Gam(n,x)=Gamma(n+ep*x)*Exp(ep*x*EulerGamma)
*      iGam(n,x)=Exp(-ep*x*EulerGamma)/Gamma(n+ep*x)
* 
* 
*   Routines from exact version of MINCER package
*   ---------------------------------------------
* 
*   IntOne    - One-loop massless self-energy (G-function)
*   DoG       - Reduction of 1-loop inegral and Pochhammer symbols
*   expansion - New name exp4d with corrected expansion for denominators with ep^0
* 
*
*   Space-time dimension shifts
*   ---------------------------
*
*   Using procedure shift4plus(n=2*k) it is possible express general d-dimensional
*   result of reduction as (d+2*k)-dimensional in terms of d=4 dimensional master integrals
* 
**************************************************************************************

* By default we use reduction with tables for BN topology 
#define REDBNTAB

**************************************************************************************
S d;
dimension d;

S  M;
V  P,p1,...,p6;
S  s1m,s2m,s3m,s4m,s5m,s6m;
I  i1,i2;

S x1,...,x6;
S x,y,z,[sqrt(x)],[x],s;

CF G,poch,po,poinv;
S k1,...,k6;
S n,n1,...,n6;
T del;

*
* used in MATAD
*
S S2,OepS2,Oep2S2,T111ep,T111ep2,B4,D3,D3ep,D4,D4ep,D5,D5ep,D6,D6ep,DM,DMep;
S DN,DNep,E3,E3ep;

* Higher weight
S D6ep2,D5ep2,D4ep2,D3ep2,DNep2,DMep2,E3ep2,E3ep3,BNep3,BNep4;

S intm1,...,intm5,intn1,intt1;
S intbm,intbmbm,intbm1,intbm2,intbn,intbnbn,intbn1,intbn2,intbn3;
S intd4,intd5,intd6,intdm,intdn,inte3,inte4;

S  diff,[p1^2],[M^2+p1^2];

S  dala,test5,test6;

S  agam,bgam,cgam,dgam,egam,fgam;
CF gm2,gm2norm,gm3,gm3norm;
CF nom,deno,Gam,iGam;
CF BN1;

* Exact version
S d0;
S  ep,epp,epp1,...,epp6;
CT ftensor,dd;
S  isum1,isum2,isum3,xpower,n0,n8;

CF rat,acc,den,dena,num,ftriangle;
CF Pochhammer,PochhammerINV,GschemeConstants;

* two-loop factorized topologies (1-loop)x(1-loop)
S intMxM;
* two-loop integral topologies
S intM00,intMM0,intMMM;
* one-loop topology
S intM0;
* Final topolgy 
S int0;

* eps1m=1/(p1.p1+M^2)^ep
S eps1m,eps2m,eps3m,eps4m,eps5m,eps6m;
* For lower topologies
S epx1,epx2,epx3,epx4,epx5,epx6;

* Zeta values
Symbols z2,z3,z4,z5,z6;

* Master integrals
CF mi;
CF miT111,miD6,miD5,miD4,miD3,miDN,miDM,miE3;
CF miBN,miBN1;
* And its truncation flags
CF Oep;
S T111tr,D6tr,D5tr,D4tr,DNtr,DMtr,E3tr;
S BNtr,D3tr,BN1tr;
S iGamtr,Gamtr;

set trunc:T111tr,D6tr,D5tr,D4tr,DNtr,DMtr,E3tr,
          BNtr,D3tr,BN1tr,iGamtr,Gamtr;

* Mass distributions
* Two-loop 
Symbols   [000],[M00],[0M0],[00M],[MM0],[M0M],[0MM],[MMM];
* Three-loop 
Symbols   [000000], [00000M], [0000M0], [0000MM], 
[000M00], [000M0M], [000MM0], [000MMM], [00M000], 
[00M00M], [00M0M0], [00M0MM], [00MM00], [00MM0M], 
[00MMM0], [00MMMM], [0M0000], [0M000M], [0M00M0], 
[0M00MM], [0M0M00], [0M0M0M], [0M0MM0], [0M0MMM], 
[0MM000], [0MM00M], [0MM0M0], [0MM0MM], [0MMM00], 
[0MMM0M], [0MMMM0], [0MMMMM], [M00000], [M0000M], 
[M000M0], [M000MM], [M00M00], [M00M0M], [M00MM0], 
[M00MMM], [M0M000], [M0M00M], [M0M0M0], [M0M0MM], 
[M0MM00], [M0MM0M], [M0MMM0], [M0MMMM], [MM0000], 
[MM000M], [MM00M0], [MM00MM], [MM0M00], [MM0M0M], 
[MM0MM0], [MM0MMM], [MMM000], [MMM00M], [MMM0M0], 
[MMM0MM], [MMMM00], [MMMM0M], [MMMMM0], [MMMMMM];

CF tad3l,tad2l,tad1l;

PolyRatFun rat;

.global



#procedure FromAuxTopo

*
*
*       AUX topo is T2
*
*       /----1->--\
*      |           |
*       ---<-3-----   
*      |           |
*       \--<-2----/        
*
*         
        id tad2l([000],n1?,n2?,n3?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3;
        id tad2l([M00],n1?,n2?,n3?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3;
        id tad2l([0M0],n1?,n2?,n3?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3;
        id tad2l([00M],n1?,n2?,n3?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3;
        id tad2l([MM0],n1?,n2?,n3?) = 1*s1m^n1*s2m^n2/p3.p3^n3;
        id tad2l([M0M],n1?,n2?,n3?) = 1*s1m^n1/p2.p2^n2*s3m^n3;
        id tad2l([0MM],n1?,n2?,n3?) = 1/p1.p1^n1*s2m^n2*s3m^n3;
        id tad2l([MMM],n1?,n2?,n3?) = 1*s1m^n1*s2m^n2*s3m^n3;                
        

* 
*          Our AUX topo is D6
* 
*           
*           _______6_______
*          |\      \      /|
*          |  \       _3/  |
*          |   _\|    /|   |
*         /|1   4 \_/     5|\
*          |      / \      |
*          |    /     \    |
*          |  /         \  |
*          |/______/______\|
*                  2
* 
* 
*         

        id tad3l([000000],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([00000M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

        id tad3l([0000M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

        id tad3l([0000MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4*s5m^n5*s6m^n6;

        id tad3l([000M00],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([000M0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3*s4m^n4/p5.p5^n5*s6m^n6;

        id tad3l([000MM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3*s4m^n4*s5m^n5/p6.p6^n6;

        id tad3l([000MMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3*s4m^n4*s5m^n5*s6m^n6;

        id tad3l([00M000],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([00M00M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

        id tad3l([00M0M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

        id tad3l([00M0MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3/p4.p4^n4*s5m^n5*s6m^n6;

        id tad3l([00MM00],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([00MM0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3*s4m^n4/p5.p5^n5*s6m^n6;

        id tad3l([00MMM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3*s4m^n4*s5m^n5/p6.p6^n6;

        id tad3l([00MMMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3*s4m^n4*s5m^n5*s6m^n6;

        id tad3l([0M0000],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([0M000M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

        id tad3l([0M00M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

        id tad3l([0M00MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3/p4.p4^n4*s5m^n5*s6m^n6;

        id tad3l([0M0M00],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([0M0M0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3*s4m^n4/p5.p5^n5*s6m^n6;

        id tad3l([0M0MM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3*s4m^n4*s5m^n5/p6.p6^n6;

        id tad3l([0M0MMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3*s4m^n4*s5m^n5*s6m^n6;

        id tad3l([0MM000],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([0MM00M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

        id tad3l([0MM0M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

        id tad3l([0MM0MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3/p4.p4^n4*s5m^n5*s6m^n6;

        id tad3l([0MMM00],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([0MMM0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3*s4m^n4/p5.p5^n5*s6m^n6;

        id tad3l([0MMMM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3*s4m^n4*s5m^n5/p6.p6^n6;

        id tad3l([0MMMMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3*s4m^n4*s5m^n5*s6m^n6;

        id tad3l([M00000],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([M0000M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

        id tad3l([M000M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

        id tad3l([M000MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4*s5m^n5*s6m^n6;

        id tad3l([M00M00],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([M00M0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3*s4m^n4/p5.p5^n5*s6m^n6;

        id tad3l([M00MM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3*s4m^n4*s5m^n5/p6.p6^n6;

        id tad3l([M00MMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3*s4m^n4*s5m^n5*s6m^n6;

        id tad3l([M0M000],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([M0M00M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

        id tad3l([M0M0M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

        id tad3l([M0M0MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3/p4.p4^n4*s5m^n5*s6m^n6;

        id tad3l([M0MM00],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([M0MM0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3*s4m^n4/p5.p5^n5*s6m^n6;

        id tad3l([M0MMM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3*s4m^n4*s5m^n5/p6.p6^n6;

        id tad3l([M0MMMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3*s4m^n4*s5m^n5*s6m^n6;

        id tad3l([MM0000],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([MM000M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

        id tad3l([MM00M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

        id tad3l([MM00MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3/p4.p4^n4*s5m^n5*s6m^n6;

        id tad3l([MM0M00],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([MM0M0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3*s4m^n4/p5.p5^n5*s6m^n6;

        id tad3l([MM0MM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3*s4m^n4*s5m^n5/p6.p6^n6;

        id tad3l([MM0MMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3*s4m^n4*s5m^n5*s6m^n6;

        id tad3l([MMM000],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([MMM00M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

        id tad3l([MMM0M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

        id tad3l([MMM0MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3/p4.p4^n4*s5m^n5*s6m^n6;

        id tad3l([MMMM00],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

        id tad3l([MMMM0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3*s4m^n4/p5.p5^n5*s6m^n6;

        id tad3l([MMMMM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3*s4m^n4*s5m^n5/p6.p6^n6;

        id tad3l([MMMMMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3*s4m^n4*s5m^n5*s6m^n6;

#endprocedure



#procedure Conv2exact()

        id po(1,?a) = 1;
        id poinv(1,?a) = 1;
        id po(x1?pos_,0) = fac_(x1-1);
        id poinv(x1?pos_,0) = 1/(fac_(x1-1));
        
        id po(n?,x?)                = Pochhammer(n,x*(2-d/2))*den(x*(2-d/2));
        id poinv(n?,x?)             = PochhammerINV(n,x*(2-d/2))*num(x*(2-d/2));
        
        id nom(x1?,x2?)             = num(x1+x2*(2-d/2));
        id deno(x1?,x2?)            = den(x1+x2*(2-d/2));

        id nom(x?,y?,z?)            = x + y*num((2-d/2))*num(1+(2-d/2)*z/y);

        id	num(x?)*den(x?) = 1;
        id	den(x?number_) = 1/x;
        id	num(x?number_) = x;
        id num(x?)                  = rat(x,1);
        id den(x?)                  = rat(1,x);        
                
#endprocedure



#procedure averts(P,in)

        if( count(int`in',1));        
        totensor,nosquare,`P',ftensor;

        id ftensor(?a) = [sqrt(x)]^nargs_(?a)*ftensor(?a);
        id [sqrt(x)]^n?odd_ = 0;
        id,many, [sqrt(x)]*[sqrt(x)] = [x];

        id ftensor(?a) = dd_(?a);
        endif;        
        .sort

        if( count(int`in',1));        
        if ( count([x],1) != 0 );
        id [x]^s? =  num(1-(2-d/2))*poinv(s+2,-1)*`P'.`P'^(s)/2^s;
        endif;
        endif;        
        #call Conv2exact
        .sort
#endprocedure

#procedure ACCU(TEXT)
        #call Conv2exact        
        .sort:`TEXT';
#endprocedure



#procedure tad1l
*
* {tad1l;1;1;0;1; ;(p1:1,1);1;0}
*
        Multiply intM0;
        #call averts(p1,M0)
        #call ACCU{}
        
        #call TadpoleM0(s1m,p1,M0,0)        
        .sort
        
        #call Conv2exact
        #call DoG
        #call subSimple
        #call GammaArgToOne

#endprocedure        



#procedure partfrac(p1,xxx)
*
* this procedure does the partial fractioning from the term
* 1/p1.p1^a xxx^b where xxx = 1/M^2+p1.p1
*
        id 1/`p1'.`p1' = 1/[p1^2];
        id `xxx' = 1/[M^2+p1^2];  
        ratio [p1^2],[M^2+p1^2],diff;
        id 1/[p1^2] = 1/`p1'.`p1';   
        id 1/[M^2+p1^2] = `xxx';     
        id 1/diff = 1/M^2;
#endprocedure





* Exact version
#procedure TadpoleMM0(x,y,z,in,out)
        
* 1/(p1.p1 + M^2) is "x" and 1/(p1.p1 + M^2)^ep  is "epx";
* 1/(p2.p2 + M^2) is "y" and 1/(p2.p2 + M^2)^ep  is "epy";
* 1/(p3.p3) is "z" and 1/p3.p3^ep  is "epz";
        if(count (int`in',1));
        if ( (count(`x',1) <= 0) && (count(ep`x',1) == 0) ) discard;
        if ( (count(`y',1) <= 0) && (count(ep`y',1) == 0) ) discard;
        
        if ( count(ep`x',1)  < count(ep`y',1) ) multiply,replace_(`x',`y',ep`x',ep`y');
        
        id `x'^k1?*ep`x'^k2?*`y'^k3?*ep`y'^k4?/`z'.`z'^k5?*ep`z'^k6? = gm3(k1,k2,k3,k4,k5,k6);

        id gm3(k1?,k2?,k3?,k4?,k5?,k6?) =
        po(2 - k5    ,      -k6 - 1)
        *po(k1 + k5 - 2, k2 + k6 + 1)
        *po(k3 + k5 - 2, k4 + k6 + 1)
        *po(k1 + k3 + k5 - 4, 2 + k2 +  k4 + k6)
        *poinv(k1                ,k2)
        *poinv(k3                ,k4)
        *poinv(k1 + k3 + 2*k5 - 4,2 + k2 +k4 + 2*k6)
        *(1 + k2 + k6)
        *(1 + k4 + k6)
        *(-1)*nom(1 ,-(2 + k2 + k4 + k6))
        *(2 + k2 + k4 + k6)
        /(2 + k2 + k4 + 2*k6)
        *nom(0,1)*nom(0,1)
        *M^(8 - 2*k1 - 2*k3 - 2*k5)
        *gm3norm(k2,k4,k6)
        ;
        Multiply int`out'/int`in';
        endif;
        
* 
        .sort:TadpoleMM0-`in'-1;        
*
        #call Conv2exact
        #call DoG
*         
        .sort:TadpoleMM0-`in'-2;        
*         
#endprocedure



#procedure TadpoleM0(x,y,in,out)

* 1/(Q.Q + M^2) is "x" and Q is "y";
* (1/(Q.Q + M^2))^ep is "epx" and (1/Q.Q)^ep is "epy";

        if( count(int`in',1));
        if ( (count(`x',1) <= 0) && (count(ep`x',1) == 0) ) discard;

        id `x'^k1?*ep`x'^k2?/`y'.`y'^k3?*ep`y'^k4?  = gm2(k1,k2,k3,k4);

        id  gm2(k1?,k2?,k3?,k4?) =
        po(2 - k3    , -k4 - 1)*
        po(k1 + k3 - 2, k2 + k4 + 1)*
        poinv(k1     ,k2)*
        M^(4 - 2*k1 - 2*k3)*
        gm2norm(k2,k4)
        ;   
        
        Multiply int`out'/int`in';        
        endif;        
* 
        .sort:TadpoleM0-`in'-1;        
*
        #call Conv2exact
        #call DoG
*         
        .sort:TadpoleM0-`in'-2;        
*         
#endprocedure

*--#[ DoG :
*
#procedure DoG
*
*	The only objects left are the G(1,x1,1,x2,0,0)
*	which have been written as GschemeConstants(x1,x2)
*
*#$vc = 0;
*Print +f "<1> %t";
id	G(n1?,x1?,n2?,x2?,n3?,n4?) = GschemeConstants(x1,x2)/(1+x1+x2)*
	Pochhammer(n1+n2-n4-2,(1+x1+x2)*(2-d/2))*
	Pochhammer(1-n1+n3-n4,1-(1+x1)*(2-d/2))*
	Pochhammer(1-n2+n4,1-(1+x2)*(2-d/2))*
	PochhammerINV(n1-1,1+x1*(2-d/2))*
	PochhammerINV(n2-1,1+x2*(2-d/2))*
	PochhammerINV(2-n1-n2+n3,2-(2+x1+x2)*(2-d/2));
*Print +f "<2> %t";
repeat id Pochhammer(n?pos_,x?) = Pochhammer(n-1,x)*num(n-1+x);
repeat id Pochhammer(n?neg_,x?) = Pochhammer(n+1,x)*den(n+x);
repeat id PochhammerINV(n?pos_,x?) = PochhammerINV(n-1,x)*den(n-1+x);
repeat id PochhammerINV(n?neg_,x?) = PochhammerINV(n+1,x)*num(n+x);
id	GschemeConstants(0,n?pos_) = GschemeConstants(n,0);
*	GschemeConstants(0,1)*GschemeConstants(0,2)
*			= GschemeConstants(1,1)*(2*D-6)/(3*D-10)
*			= GschemeConstants(1,1)*(1-2*ep)/(1-3*ep)
id	GschemeConstants(1,1)*GschemeConstants(0,0) =
			GschemeConstants(1,0)*GschemeConstants(2,0)*rat(1-3*(2-d/2),1-2*(2-d/2));
id	Pochhammer(0,x?) = 1;
id	PochhammerINV(0,x?) = 1;
id	num(x?)*den(x?) = 1;
id	den(x?number_) = 1/x;
id	num(x?number_) = x;
*Print +f "<3> %t";
id	num(x?) = rat(x,1);
*Print +f "<4> %t";
id	den(x?) = rat(1,x);
*Print +f "<5> %t";
*$vc = $vc+1;
*Print +f "<$vc = %$>",$vc;
*
#endprocedure
*
*--#] DoG : 

*--#[ IntOne :
*
#procedure IntOne(p3,p4,Q,in,out)
*
if ( count(int`in',1) );
  if ( ( count(ep`p3',1) == 0 ) && ( count(`p3'.`p3',1) >= 0 ) ) Discard;
  if ( ( count(ep`p4',1) == 0 ) && ( count(`p4'.`p4',1) >= 0 ) ) Discard;
  ToTensor,nosquare,ftensor,`p3';
  if ( count(ftensor,1) == 0 );
	id	int`in'*ep`p3'^x3?*ep`p4'^x4?/`p3'.`p3'^n3?/`p4'.`p4'^n4? =
			int`out'*G(n3,x3,n4,x4,0,0)*`Q'.`Q'^2/`Q'.`Q'^n3/`Q'.`Q'^n4*ep`Q'^x3*ep`Q'^x4*ep`Q';
  elseif ( match(ftensor(i1?)) );
	id	int`in'*ep`p3'^x3?*ep`p4'^x4?/`p3'.`p3'^n3?/`p4'.`p4'^n4?*ftensor(i1?) = int`out'*`Q'(i1)
			*G(n3,x3,n4,x4,1,0)*`Q'.`Q'^2/`Q'.`Q'^n3/`Q'.`Q'^n4*ep`Q'^x3*ep`Q'^x4*ep`Q';
  elseif ( match(ftensor(i1?,i2?)) );
	id	int`in'*ep`p3'^x3?*ep`p4'^x4?/`p3'.`p3'^n3?/`p4'.`p4'^n4?*ftensor(i1?,i2?) =
				int`out'*`Q'.`Q'^2*ep`Q'*ep`Q'^x3*ep`Q'^x4/`Q'.`Q'^n3/`Q'.`Q'^n4*(
			+G(n3,x3,n4,x4,2,0)*`Q'(i1)*`Q'(i2)
			+G(n3,x3,n4,x4,2,1)*d_(i1,i2)*`Q'.`Q'/2);
  else;
	id	int`in'*ep`p3'^x3?*ep`p4'^x4?/`p3'.`p3'^n3?/`p4'.`p4'^n4?*ftensor(?a) = int`out'*ftensor(?a)
			*sum_(isum2,0,integer_(nargs_(?a)/2),G(n3,x3,n4,x4,nargs_(?a),isum2)
				*y^isum2*`Q'.`Q'^isum2/2^isum2)*`Q'.`Q'^2
				*ep`Q'*ep`Q'^x3*ep`Q'^x4/`Q'.`Q'^n3/`Q'.`Q'^n4;
    id  y^isum2?*ftensor(?a) = distrib_(1,2*isum2,del,ftensor,?a);
    tovector,ftensor,`Q';
    id  del(?a) = dd_(?a);
  endif;
  id  P.P = 0;
endif;
*
.sort:IntOne-`in'-1;
*
#call DoG
*
.sort:IntOne-`in'-2;
*
#endprocedure
*
*--#] IntOne : 




#procedure tad2l
*
* {tad2l;3;2;0;1; ;(p1:1,2)(p2:2,1)(p3:2,1);111;110;101;011;100;010;001;000}
*

*
*          /------\
*         /    |   \
*      p1v   p3^    ^p2
*         \    |   /
*          \------/
*             

        #call partfrac{p1|s1m}
        #call partfrac{p2|s2m}
        #call partfrac{p3|s3m}
        
        if (count(s1m,1,s2m,1,s3m,1) == 0) discard;
        
* map 010 -> 100
        if ( (count(s1m,1) == 0) && (count(s2m,1) != 0) && (count(s3m,1) == 0) );
        id p1 = -p1;
        id p2 = -p2;
        multiply,replace_(p1,p2,p2,p1,s1m,s2m,s2m,s1m);
        endif;

* map 001 -> 100
        if ( (count(s1m,1) == 0) && (count(s2m,1) == 0) && (count(s3m,1) != 0) );
        id p1 = -p1;
        id p3 = -p3;
        multiply,replace_(p1,p3,p3,p1,s1m,s3m,s3m,s1m);
        endif;

* map 101 -> 110
        if ( (count(s1m,1) != 0) && (count(s2m,1) == 0) && (count(s3m,1) != 0) );
        multiply,replace_(p2,p3,p3,p2,s2m,s3m,s3m,s2m);
        endif;

* map 011 -> 110
        if ( (count(s1m,1) == 0) && (count(s2m,1) != 0) && (count(s3m,1) != 0) );
        id p1 = -p1;
        id p3 = -p3;
        multiply,replace_(p1,p3,p3,p1,s1m,s3m,s3m,s1m);
        endif;

* the rest is taken from topT111 with line 5 replaced by line 3

*
* Now the integration is done:
* (only the momenta p1, p2 and p3 appear)
*

*
* Warning: change direction of line 1 if lines 1 and 2 are massive
*
        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) ); 
        id p1=-p1;
        endif;
        .sort

        #message decompose numerator (M|M|0)

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) ); 
        id p3 = -p1-p2;
        endif;

        #call ACCU{nomgm3 1}

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) ); 
        id p1.p1 = 1/test5 - M^2;
        id p2.p2 = 1/test6 - M^2;
        endif;

        #call ACCU{nomgm3 2}

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) ); 
        id p1.p2 = (p3.p3 - 1/test5 - 1/test6 + 2*M^2)/2;
        endif;

        if (count(s3m,1)>0) id, p3.p3 = 1/s3m - M^2;
        id 1/s3m = p3.p3 + M^2;

        #call ACCU{nomgm3 3}

        #message rec. rel. for the case with 3 massive lines

        #do i=1,1
                #message loop
                if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) ); 
                multiply replace_(test5,s1m,test6,s2m);
                if ( (count(s1m,1)>1) && (count(s2m,1)>=1) && (count(s3m,1)>=1) );
                id s1m^n1?*s2m^n2?*s3m^n5? = 
                s1m^n1*s2m^n2*s3m^n5 * (-1) * 1/3/(n1-1)/M^2 * (
                nom(4+3-3*n1, -2)/s1m
                + 2*n2*s2m/s1m*(1/s3m-1/s1m)
                - (n1-1)*(1/s3m-1/s2m)
                );
                redefine i "0";
                endif;
                if ( count(s1m,1) < count(s2m,1) );
                multiply replace_(s1m,s2m,s2m,s1m);
                redefine i "0";
                endif;
                if ( count(s1m,1) < count(s3m,1) );
                multiply replace_(s1m,s3m,s3m,s1m);
                redefine i "0";
                endif;
                if ( count(s2m,1) < count(s3m,1) );
                multiply replace_(s2m,s3m,s3m,s2m);
                redefine i "0";
                endif;
                endif;

                id	num(x?)*den(x?) = 1;
                id	den(x?number_) = 1/x;
                id	num(x?number_) = x;
                id	num(x?) = rat(x,1);
                id	den(x?) = rat(1,x);
                
                .sort                
        #enddo

        #message done
        
        .sort

        #message perform integration

        id 1/s3m = p3.p3 + M^2;

        #call ACCU{T111}

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) ); 
        multiply replace_(test5,s1m,test6,s2m);
        Multiply intMM0;
        
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) );
        Multiply intM00;
 
        elseif ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s3m,1)==1) ); 
        id s1m*s2m*s3m = M^2*miT111;
        Multiply int0;
        
        else;
        multiply 1/(1-1);
        endif;
        .sort 
       
* MM0 case
        #call TadpoleMM0(s1m,s2m,p3,MM0,0)
* M00 case        
        #call IntOne(p2,p3,p1,M00,M0)
        #call averts(p1,M0)

* M0 case        
        #call TadpoleM0(s1m,p1,M0,0)

        #call Conv2exact
        #call DoG
        #call subSimple
        #call GammaArgToOne
        
#endprocedure




#procedure bnm2m(TYPE)
*
* Identify M1, M2, ... in the output of BN/BM
*
        if (count(int`TYPE',1) );

        if (count(intm1,1,intm2,1,intm3,1,intm4,1,intt1,1,intn1,1,
        intbm1,1,intbm2,1)==0);
        
        if (match(1/p1.p1/p2.p2*x5*x6)>0); 
        multiply intm1/int`TYPE';
        elseif (match(1/p1.p1/p2.p2/p5.p5*x6)>0); 
        multiply intm2/int`TYPE';
        elseif (match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6)>0); 
        multiply intm3/int`TYPE';
        elseif (match(1/p3.p3/p4.p4*x5*x6)>0); 
        multiply intm4/int`TYPE';
        elseif( (count(p1.p1,1)==0) && (count(p2.p2,1)==0)
        && (count(x3,1)==0) && (count(x4,1)==1)
        && (count(x5,1)==1) && (count(x6,1)==1) );
        multiply intt1/int`TYPE';
        elseif( (count(p1.p1,1)==0) && (count(p2.p2,1)==0)
        && (count(x3,1)==1) && (count(x4,1)==1)
        && (count(x5,1)==1) && (count(x6,1)>0) );
        multiply intn1/int`TYPE';
        endif;
        endif;
        
* BN BM
        endif;

#endprocedure        




#procedure dorec3l
*
* dorec3l -> master procedure for three loop tadpole recurrence relations
* 
************************************************************ 

* Note:
* - Treat first d6, d5, ...; then bn, ... and at last bm, ...
* - MINCER calculates the diagrams originally in the G-scheme
*   and changes at the end to the MSbar-scheme. In addition 
*   a multiplication of exp(z2*\ep^2/2) per loop is done.
* - The topologies BN1 and BN2 are reduced to topologies BM, BM1 
*   and simpler functions BN1(...) (of course topology BN1 only, 
*   this should be changed in future), so terms proportional to
*   intbm, intbm1 and BN1(...) are present inside intbn1 and intbn2
*   which we have to be set equal to zero, in order not to count anything
*   twice. (BN1(...) are added to 'diarest'; see below.) 

************************************************************ 
        
        #do type = {d6|d5|d4|dm|dn|e4|e3}
                
                
                #message Recursion of type `type'
                
                #call top`type'
                .sort        
        #enddo
        
        #do type = {bn|bn1|bn2|bn3|bm|bm1|bm2}
        
                #message Recursion of type `type'
                multiply replace_(s1m,x1,s2m,x2,s3m,x3,s4m,x4,s5m,x5,s6m,x6);
                
                #call top`type'
                #call bnm2m(`type')
                .sort        
                
        #enddo
        
        #do type = {m1|m2|m3|m4|m5|t1|n1}
                #message Recursion of type `type'
                multiply replace_(s1m,x1,s2m,x2,s3m,x3,s4m,x4,s5m,x5,s6m,x6);
                
                #call top`type'
        #enddo

#endprocedure






************************************************************

#procedure mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)
        
* discard massless tadpoles
        if ( count(intd6,1) );
        if ( (count(`s1m',1)<=0)&&(count(`s2m',1)<=0)&&(count(`s4m',1)<=0) ) discard;
        if ( (count(`s1m',1)<=0)&&(count(`s3m',1)<=0)&&(count(`s6m',1)<=0) ) discard;
        if ( (count(`s2m',1)<=0)&&(count(`s3m',1)<=0)&&(count(`s5m',1)<=0) ) discard;
        if ( (count(`s4m',1)<=0)&&(count(`s5m',1)<=0)&&(count(`s6m',1)<=0) ) discard;
        endif;        
        .sort

#endprocedure

#procedure redD6n5(s1m,s2m,s3m,s4m,s5m,s6m)

* reduce n5 from >1 to =1

        if (match(`s1m'*`s2m'*`s3m'*`s4m'*`s5m'^2*`s6m') > 0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?*`s6m'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5*`s6m'^n6 * 
        (-1)/4/(n5-1)/M^2/`s5m' * (
        num(d)-4*(n5-1)
        -3*( n3*`s3m'*(1/`s5m'-1/`s6m') + n2*`s2m'*(1/`s5m'-1/`s4m') )
        +n2*`s2m'*(1/`s3m'-1/`s1m') 
        +n3*`s3m'*(1/`s2m'-1/`s1m') 
        +(n5-1)*`s5m'*(1/`s3m'-1/`s6m'+1/`s2m'-1/`s4m')
        );

* sort: n5 >= n1,n2,n3,n4,n6

        if (count(`s1m',1) > count(`s5m',1)) 
        multiply replace_(`s1m',`s3m',`s3m',`s1m',`s4m',`s5m',`s5m',`s4m');
        if (count(`s5m',1) < count(`s6m',1)) 
        multiply replace_(`s5m',`s6m',`s6m',`s5m',`s1m',`s2m',`s2m',`s1m');
        if (count(`s5m',1) < count(`s4m',1)) 
        multiply replace_(`s5m',`s4m',`s4m',`s5m',`s1m',`s3m',`s3m',`s1m');
        if (count(`s5m',1) < count(`s3m',1)) 
        multiply replace_(`s5m',`s3m',`s3m',`s5m',`s1m',`s4m',`s4m',`s1m');
        if (count(`s5m',1) < count(`s2m',1)) 
        multiply replace_(`s5m',`s2m',`s2m',`s5m',`s1m',`s6m',`s6m',`s1m');

        redefine i "0";

        endif;

#endprocedure



************************************************************
#procedure topd6

* treat the scalar products


*
* this is topd6
*

        #message this is topd6

        #message numerator

        if ( count(intd6,1) );
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p5.p5 = 1/s5m - M^2;
        id p6.p6 = 1/s6m - M^2;
        endif;
        #call ACCU(D6 0)
        #call mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)

        if ( count(intd6,1) );        
        id p2 = p1+p3;
        id p4 = p1+p6;
        id p5 = p6-p3;
        endif;        
        #call ACCU(D6 1)

        if ( count(intd6,1) );        
        id p1.p3 = 1/2 * (   1/s2m - 1/s1m - 1/s3m + M^2 );
        id p1.p6 = 1/2 * (   1/s4m - 1/s1m - 1/s6m + M^2 );
        id p3.p6 = 1/2 * ( - 1/s5m + 1/s3m + 1/s6m - M^2 );
        endif;        
        #call ACCU(D6 2)
        
        #call mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)

        if ( count(intd6,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p5.p5 = 1/s5m - M^2;
        id p6.p6 = 1/s6m - M^2;
        endif;        
        #call ACCU(D6 3)
        
        #call mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)
        
        #call ACCU(D6)
        
************************************************************
        
* do recursion
        
        #message do recursion
        
* sort: n5 >= n1,n2,n3,n4,n6
        if ( count(intd6,1) );
        if (count(s1m,1) > count(s5m,1)) 
        multiply replace_(s1m,s3m,s3m,s1m,s4m,s5m,s5m,s4m);
        if (count(s5m,1) < count(s6m,1)) 
        multiply replace_(s5m,s6m,s6m,s5m,s1m,s2m,s2m,s1m);
        if (count(s5m,1) < count(s4m,1)) 
        multiply replace_(s5m,s4m,s4m,s5m,s1m,s3m,s3m,s1m);
        if (count(s5m,1) < count(s3m,1)) 
        multiply replace_(s5m,s3m,s3m,s5m,s1m,s4m,s4m,s1m);
        if (count(s5m,1) < count(s2m,1)) 
        multiply replace_(s5m,s2m,s2m,s5m,s1m,s6m,s6m,s1m);
        endif;        

        #do i=1,1
                #call redD6n5(s1m,s2m,s3m,s4m,s5m,s6m)
                .sort
        #enddo
        
        #call mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)
        
        #call ACCU(D6)
        
************************************************************
        
* identify simpler integrals
        if ( count(intd6,1) );
        if ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2; 
        id p1=p4-p6;
        id p2=-p2;
        id p4=-p4;
        multiply replace_(p2,p1,p3,p4,p4,p3,p5,p2,p6,p5,
	s2m,s1m,s3m,s4m,s4m,s3m,s5m,s2m,s6m,s5m);
        multiply, intd5/intd6;

        elseif ( (count(s2m,1)<=0) );
        id 1/s2m=p2.p2+M^2; 
        id p2=p1+p3;
        id p1=-p1;
        id p3=-p3;
        multiply replace_(p1,p3,p3,p5,p4,p1,p5,p4,p6,p2,
	s1m,s3m,s3m,s5m,s4m,s1m,s5m,s4m,s6m,s2m);
        multiply, intd5/intd6;

        elseif ( (count(s3m,1)<=0) );
        id 1/s3m=p3.p3+M^2; 
        id p3=p6-p5;
        id p1=-p1;
        id p2=-p2;
        id p4=-p4;
        id p6=-p6;
        multiply replace_(p2,p4,p4,p2,p6,p3,
        s2m,s4m,s4m,s2m,s6m,s3m);
        multiply, intd5/intd6;

        elseif ( (count(s4m,1)<=0) );
        id 1/s4m=p4.p4+M^2; 
        id p4=p2+p5;
        id p2=-p2;
        id p3=-p3;
        multiply replace_(p1,p3,p2,p1,p3,p2,p5,p4,p6,p5,
        s1m,s3m,s2m,s1m,s3m,s2m,s5m,s4m,s6m,s5m);
        multiply, intd5/intd6;

        elseif ( (count(s5m,1)<=0) );
        id 1/s5m=p5.p5+M^2; 
        id p5=p6-p3;
        id p1=-p1;
        id p2=-p2;
        multiply replace_(p1,p2,p2,p3,p3,p1,p4,p5,p6,p4,
        s1m,s2m,s2m,s3m,s3m,s1m,s4m,s5m,s6m,s4m);
        multiply, intd5/intd6;

        elseif ( (count(s6m,1)<=0) );
        id 1/s6m=p6.p6+M^2; 
        id p6=p4-p1;
        multiply, intd5/intd6;

        elseif ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s3m,1)==1) && 
        (count(s4m,1)==1) && (count(s5m,1)==1) && (count(s6m,1)==1) );
        id s1m*s2m*s3m*s4m*s5m*s6m = miD6;
        Multiply int0/intd6;

        else;
        exit "D6: Unknown simpler topology";
        endif;
        endif;
        
        #message - done

#endprocedure        



************************************************************

#procedure mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

* discard massless tadpoles
        if ( count(intd5,1) );
        if ( (count(`s1m',1)<=0) && (count(`s3m',1)<=0) ) discard;
        if ( (count(`s4m',1)<=0) && (count(`s5m',1)<=0) ) discard;
        endif;        
        .sort
        
#endprocedure

#procedure redD5n6p(s1m,s2m,s3m,s4m,s5m,p6)

* reduce n6 from >0 to =0
        if ( count(intd5,1) );        
        if (match(`s1m'*`s2m'*`s3m'*`s4m'*`s5m'/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 
        
        (-1)*deno(4-2*n6-n3-n1,-2)*(
        + n1*`s1m'*`s4m'^-1
        + n3*`s3m'*`s5m'^-1
        - `p6'.`p6'*n1*`s1m'
        - `p6'.`p6'*n3*`s3m'
        );
        
        endif;
* topd5                
        endif;        
#endprocedure

#procedure redD5n6m(s1m,s2m,s3m,s4m,s5m,p6)

* reduce n6 from <0 to =0
        if ( count(intd5,1) );        
        if (match(`s1m'*`s2m'*`s3m'*`s4m'*`s5m'*`p6'.`p6')>0);
        
        if (count(`s3m',1)==1);
        
        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 
        
        (-1)*deno(4-2*n3-n6-n1,-2)*(
        + n1*M^2*`s1m'                 
        + n1*`s1m'*`s2m'^-1            
        - n1*`s1m'*`s3m'^-1            
        + 2*n3*M^2*`s3m'               
        - `p6'.`p6'^-1*n6*`s3m'^-1     
        + `p6'.`p6'^-1*n6*`s5m'^-1     
        );
        
* sort: n3 >= n4,n5,n1
        
        repeat;
                if (count(`s3m',1) < count(`s4m',1)) 
                multiply replace_(`s3m',`s4m',`s4m',`s3m',`s1m',`s5m',`s5m',`s1m');
                if (count(`s3m',1) < count(`s5m',1)) 
                multiply replace_(`s3m',`s5m',`s5m',`s3m',`s1m',`s4m',`s4m',`s1m');
                if (count(`s3m',1) < count(`s1m',1)) 
                multiply replace_(`s3m',`s1m',`s1m',`s3m',`s5m',`s4m',`s4m',`s5m');
        endrepeat;
        
        redefine i "0";
        
        endif;

        endif;


        if (match(`s1m'*`s2m'*`s3m'*`s4m'*`s5m'*`p6'.`p6')>0);
        
* 29Oct04: added
        if (count(`s3m',1)!=1);
        
        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 
        
        (-1)/(n3-1)*(
        + `p6'.`p6'^-1*n2*M^2*`s2m'*`s3m'^-1     
        + `p6'.`p6'^-1*n2*`s2m'*`s3m'^-1*`s4m'^-1 
        - `p6'.`p6'^-1*n2*`s2m'*`s3m'^-1*`s5m'^-1
        - `p6'.`p6'^-1*n2*`s3m'^-1
        + 2*`p6'.`p6'^-1*n5*M^2*`s3m'^-1*`s5m'
        - 2*`p6'.`p6'^-1*n5*`s3m'^-1
        + 2*`p6'.`p6'^-1*M^2*(n3-1)
        - `p6'.`p6'^-1*(n3-1)*`s3m'^-1
        - `p6'.`p6'^-1*(n3-1)*`s5m'^-1
        + `p6'.`p6'^-1*`s3m'^-1*num(d)
        );
        
* 29Oct04: added
        endif;
        
* 29Oct04: added
        redefine i "0";
        
        endif;
* topd5        
        endif;        
#endprocedure


#procedure redD5n2(s1m,s2m,s3m,s4m,s5m,p6)

* reduce n2 from >1 to =1
* applied for n6=0

        if ( count(intd5,1) );        
        if ( (match(`s1m'*`s2m'^2*`s3m'*`s4m'*`s5m')>0) && (count(`p6',1)==0) );
        
        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 
        
        (-1)/(n2-1)*(
        - n1*M^-2*`s2m'^-1
        + n1*`s1m'*`s2m'^-1
        - n3*M^-2*`s2m'^-1
        + n3*`s2m'^-1*`s3m'
        - n4*M^-2*`s2m'^-1
        + n4*`s2m'^-1*`s4m'
        - n5*M^-2*`s2m'^-1
        + n5*`s2m'^-1*`s5m'
        - n6*M^-2*`s2m'^-1
        - M^-2*(n2-1)*`s2m'^-1
        + 3/2*M^-2*`s2m'^-1*num(d)
        );
        
        redefine i "0";
        
        endif;
* topd5        
        endif;
#endprocedure


#procedure redD5n1345(s1m,s2m,s3m,s4m,s5m,p6)

* reduce n1,n3,n4 and n5 from >1 to =1
* requirements: n2=1 and n6=0 

        if ( count(intd5,1) );        
        if ( (match(`s1m'^2*`s2m'*`s3m'*`s4m'*`s5m')>0) && 
        (count(`p6',1)==0) && (count(`s2m',1)==1) );

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 

        (-1)/(n1-1)*(
        - 2/3*n3*M^-2*`s1m'^-2*`s3m'
        + 2/3*n3*M^-2*`s1m'^-1*`s2m'^-1*`s3m'
        - 1/3*n6*M^-2*`s1m'^-1
        - M^-2*(n1-1)*`s1m'^-1
        - 1/3*M^-2*(n1-1)*`s2m'^-1
        + 1/3*M^-2*(n1-1)*`s3m'^-1
        + 1/3*M^-2*`s1m'^-1*num(d)
        - 2/3*`p6'.`p6'^-1*n6*M^-2*`s1m'^-2
        + 1/3*`p6'.`p6'^-1*n6*M^-2*`s1m'^-1*`s3m'^-1
        + 2/3*`p6'.`p6'^-1*n6*M^-2*`s1m'^-1*`s4m'^-1
        - 1/3*`p6'.`p6'^-1*n6*M^-2*`s1m'^-1*`s5m'^-1
        );

* sort: n1 >= n3,n4,n5

        repeat;
                if (count(`s1m',1) < count(`s3m',1)) 
                multiply replace_(`s1m',`s3m',`s3m',`s1m',`s4m',`s5m',`s5m',`s4m');
                if (count(`s1m',1) < count(`s4m',1)) 
                multiply replace_(`s1m',`s4m',`s4m',`s1m',`s3m',`s5m',`s5m',`s3m');
                if (count(`s1m',1) < count(`s5m',1)) 
                multiply replace_(`s1m',`s5m',`s5m',`s1m',`s3m',`s4m',`s4m',`s3m');
        endrepeat;

        if ( (count(`s2m',1)==0) && (count(`p6',1)==0) );
        if (count(`s1m',1)<count(`s3m',1)) multiply replace_(`s1m',`s3m',`s3m',`s1m');
        if (count(`s1m',1)<count(`s4m',1)) multiply replace_(`s1m',`s4m',`s4m',`s1m');
        if (count(`s1m',1)<count(`s5m',1)) multiply replace_(`s1m',`s5m',`s5m',`s1m');
        if (count(`s3m',1)<count(`s4m',1)) multiply replace_(`s3m',`s4m',`s4m',`s3m');
        if (count(`s3m',1)<count(`s5m',1)) multiply replace_(`s3m',`s5m',`s5m',`s3m');
        if (count(`s4m',1)<count(`s5m',1)) multiply replace_(`s4m',`s5m',`s5m',`s4m');
        endif;

        if ( (count(`s3m',1)==0) && (count(`p6',1)==0) );
        if (count(`s2m',1)<count(`s4m',1)) multiply replace_(`s2m',`s4m',`s4m',`s2m');
        if (count(`s2m',1)<count(`s5m',1)) multiply replace_(`s2m',`s5m',`s5m',`s2m');
        if (count(`s4m',1)<count(`s5m',1)) multiply replace_(`s4m',`s5m',`s5m',`s4m');
        endif;

        redefine i "0";

        endif;
* topd5  
        endif;
#endprocedure

#procedure topd5
*
* this is topd5
*
        #message this is topd5

************************************************************

* treat the scalar products

        #message numerator

        if ( count(intd5,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p5.p5 = 1/s5m - M^2;
        endif;        
        #call ACCU(D5 0)

        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)
        
        if ( count(intd5,1) );        
        id p2 = p1+p3;
        id p4 = p1+p6;
        id p5 = p6-p3;
        endif;        
        #call ACCU(D5 1)

        if ( count(intd5,1) );                
        id  p1.p3 = 1/2 * (   1/s2m - 1/s1m - 1/s3m + M^2 );
        id  p1.p6 = 1/2 * (   1/s4m - 1/s1m - p6.p6 );
        id  p3.p6 = 1/2 * ( - 1/s5m + 1/s3m + p6.p6 );
        endif;        
        #call ACCU(D5 2) 

        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

        if ( count(intd5,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p5.p5 = 1/s5m - M^2;
        endif;        
        #call ACCU(D5 3)

        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

        #call ACCU(D5)


************************************************************

* do recursion

        #message do recursion

* sort: n3 >= n4,n5,n1
* n3 must be the largest. This is needed for the reduction of
* n6 from <0 to =0. 
        if ( count(intd5,1) );        
        repeat;
                if (count(s3m,1) < count(s4m,1)) 
                multiply replace_(s3m,s4m,s4m,s3m,s1m,s5m,s5m,s1m);
                if (count(s3m,1) < count(s5m,1)) 
                multiply replace_(s3m,s5m,s5m,s3m,s1m,s4m,s4m,s1m);
                if (count(s3m,1) < count(s1m,1)) 
                multiply replace_(s3m,s1m,s1m,s3m,s5m,s4m,s4m,s5m);
        endrepeat;
        endif;        
        .sort

        if ( count(intd5,1) );                
        repeat;
                #call redD5n6p(s1m,s2m,s3m,s4m,s5m,p6)
        endrepeat;
        endif;        
        .sort

        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

        #call ACCU(D5)

        #do i=1,1

                #call redD5n6m(s1m,s2m,s3m,s4m,s5m,p6)

* sort: n3 >= n4,n5,n1
* n3 must be the largest. This is needed for the reduction of
* n6 from <0 to =0. 
                if ( count(intd5,1) );        
                repeat;
                        if (count(s3m,1) < count(s4m,1)) 
                        multiply replace_(s3m,s4m,s4m,s3m,s1m,s5m,s5m,s1m);
                        if (count(s3m,1) < count(s5m,1)) 
                        multiply replace_(s3m,s5m,s5m,s3m,s1m,s4m,s4m,s1m);
                        if (count(s3m,1) < count(s1m,1)) 
                        multiply replace_(s3m,s1m,s1m,s3m,s5m,s4m,s4m,s5m);
                endrepeat;
                endif;
                #call ACCU(D5)
        #enddo

        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)
        
        #call ACCU(D5)
        
        #do i=1,1
                
                #call redD5n2(s1m,s2m,s3m,s4m,s5m,p6)
                #call ACCU(D5)
                
        #enddo
        
* sort: n1 >= n3,n4,n5
        if ( count(intd5,1) );        
        repeat;
                if (count(s1m,1) < count(s3m,1)) 
                multiply replace_(s1m,s3m,s3m,s1m,s4m,s5m,s5m,s4m);
                if (count(s1m,1) < count(s4m,1)) 
                multiply replace_(s1m,s4m,s4m,s1m,s3m,s5m,s5m,s3m);
                if (count(s1m,1) < count(s5m,1)) 
                multiply replace_(s1m,s5m,s5m,s1m,s3m,s4m,s4m,s3m);
        endrepeat;
        endif;
        .sort
        
        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)
        
        #call Conv2exact        
        
        #do i=1,1
                
                #call redD5n1345(s1m,s2m,s3m,s4m,s5m,p6)

                #call ACCU(D5)
                
        #enddo
        
        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)
        
        #call ACCU(D5)
        
        if ( count(intd5,1) );        
        if ( (count(s2m,1)==0) && (count(p6,1)==0) );
        if (count(s1m,1) < count(s3m,1)) multiply replace_(s1m,s3m,s3m,s1m);
        if (count(s1m,1) < count(s4m,1)) multiply replace_(s1m,s4m,s4m,s1m);
        if (count(s1m,1) < count(s5m,1)) multiply replace_(s1m,s5m,s5m,s1m);
        if (count(s3m,1) < count(s4m,1)) multiply replace_(s3m,s4m,s4m,s3m);
        if (count(s3m,1) < count(s5m,1)) multiply replace_(s3m,s5m,s5m,s3m);
        if (count(s4m,1) < count(s5m,1)) multiply replace_(s4m,s5m,s5m,s4m);
        endif;
        if ( (count(s3m,1)==0) && (count(p6,1)==0) );
        if (count(s2m,1) < count(s4m,1)) multiply replace_(s2m,s4m,s4m,s2m);
        if (count(s2m,1) < count(s5m,1)) multiply replace_(s2m,s5m,s5m,s2m);
        if (count(s4m,1) < count(s5m,1)) multiply replace_(s4m,s5m,s5m,s4m);
        endif;
* topd5
        endif;
        
        #call ACCU(D5)
        
* transform D5(1,1,1,1,1,0) to D5(1,1,1,1,1,1)
*
* or: compute D5(1,1,1,1,1,0) once and use it as master integral ?
        
        if ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s3m,1)==1) &&
        (count(s4m,1)==1) && (count(s5m,1)==1) && (count(p6,1)==0) );
        id s1m*s2m*s3m*s4m*s5m = 
        
        M^2*nom(0,-2)*deno(2-8/3,4/3) * (
        
        s1m*s2m*s3m*s4m*s5m/p6.p6 - (
        
        + s1m^2*s2m*s3m*s5m * (
        - 2*deno(0,-2)*p6.p6^-1
        )
        
        + s1m^2*s2m*s4m*s5m * (
        + 2/3*deno(0,-2)*M^-2
        )
        
        + s1m^2*s3m*s4m*s5m * (
        - 2/3*deno(0,-2)*M^-2
        )
        )
        
        );
        endif;
        .sort
        
* sort: n5 <= n1,n3,n4
        if ( count(intd5,1) );        
        repeat;
                if (count(s5m,1) > count(s1m,1)) 
                multiply replace_(s5m,s1m,s1m,s5m,s4m,s3m,s3m,s4m);
                if (count(s5m,1) > count(s3m,1)) 
                multiply replace_(s5m,s3m,s3m,s5m,s1m,s4m,s4m,s1m);
                if (count(s5m,1) > count(s4m,1)) 
                multiply replace_(s5m,s4m,s4m,s5m,s3m,s1m,s1m,s3m);
        endrepeat;
        endif;
        .sort
        
        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)
        
        #call ACCU(D5)
        
************************************************************
        
* identify simple integrals
        if ( count(intd5,1) );        
        if ( (count(s2m,1)<=0) );
        id 1/s2m=p2.p2+M^2; 
        id p2=p1+p3;
        id p1=-p1;
        id p4=-p4;
        id p5=-p5;
        multiply replace_(p1,p4,p3,p6,p4,p5,p5,p3,p6,p1,
        s1m,s4m,s3m,s6m,s4m,s5m,s5m,s3m);
        
*** needed because of topBN
        id 1/s3m=p3.p3+M^2;
        multiply, intbn/intd5;

        elseif ( (count(s5m,1)<=0) );
        id 1/s5m=p5.p5+M^2; 
        id p5=-p5;
        multiply replace_(p2,p4,p3,p6,p4,p2,p6,p3,
        s2m,s4m,s3m,s6m,s4m,s2m);
        multiply, intd4/intd5;

        elseif ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s3m,1)==1) && 
        (count(s4m,1)==1) && (count(s5m,1)==1) && (count(p6.p6,1)==-1) );
        id s1m*s2m*s3m*s4m*s5m/p6.p6 = miD5;
        Multiply int0/intd5;

        else;
        exit "D5: Unknown simpler topology";
        endif;
        endif;
        
        #message - done
#endprocedure











#procedure mltadD4(s1m,s2m,p3,s4m,p5,s6m)
        
* discard massless tadpoles
        if ( count(intd4,1) );
        if ( (count(`s1m',1)<=0) && (count(`s6m',1)<=0) ) discard;
        if ( (count(`s4m',1)<=0) && (count(`s6m',1)<=0) ) discard;
        if ( (count(`s2m',1)<=0) && (count(`p3'.`p3',1)>=0) ) discard;
        if ( (count(`s2m',1)<=0) && (count(`p5'.`p5',1)>=0) ) discard;

        if ( (count(`p3'.`p3',1)>=0) && (count(`p5'.`p5',1)>=0) &&
        (count(`s1m',1)<=0)     && (count(`s6m',1)<=0) ) discard;
        endif;
        .sort

#endprocedure

#procedure redD4n35m(s1m,s2m,p3,s4m,p5,s6m)

* reduce n3,n5 from <0 to =0
* (rec. rel. implemented for n3)

        if ( count(intd4,1) );
        if ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m'^2)>0) &&
        (count(p5.p5,1)!=0) );


* for the special case (n3=-2, n5=-1) reverse order of n3 and n5
        if ( (count(`p3'.`p3',1)==2) && (count(`p5'.`p5',1)==1) ) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');


        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n6-1)*(
        - `p3'.`p3'^-1*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        - `p3'.`p3'^-1*`p5'.`p5'*(n6-1)
        + `p3'.`p3'^-1*n4*`s2m'^-1*`s4m'*`s6m'^-1
        - `p3'.`p3'^-1*n4*`s6m'^-1
        - 2*`p3'.`p3'^-1*n5*`s6m'^-1
        + `p3'.`p3'^-1*M^2*(n6-1)
        - `p3'.`p3'^-1*(n6-1)*`s6m'^-1
        + `p3'.`p3'^-1*`s6m'^-1*num(d)
        );

        redefine i "0";

        elseif ( (match(`s1m'^2*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s6m',1)==1) && (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n1-1)*(
        + `p3'.`p3'^-1*`p5'.`p5'*n4*`s1m'^-1*`s4m'
        + 2*`p3'.`p3'^-1*n3*`s1m'^-1
        - `p3'.`p3'^-1*n4*`s1m'^-1*`s2m'^-1*`s4m'
        + `p3'.`p3'^-1*n4*`s1m'^-1
        + 2*`p3'.`p3'^-1*n5*`s1m'^-1
        - 2*`p3'.`p3'^-1*n6*M^2*`s1m'^-1*`s6m'
        + 2*`p3'.`p3'^-1*n6*`s1m'^-1
        + `p3'.`p3'^-1*(n1-1)*`s1m'^-1
        - `p3'.`p3'^-1*(n1-1)*`s2m'^-1
        - 2*`p3'.`p3'^-1*`s1m'^-1*num(d)
        + 2*`p3'.`p3'^-1*`s1m'^-1
        );

        redefine i "0";

        elseif ( (match(`s1m'*`s2m'^2*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s6m',1)==1) && (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n2-1)*(
        + `p3'.`p3'^-1*n1*M^2*`s1m'*`s2m'^-1
        - `p3'.`p3'^-1*n1*`s2m'^-1
        + `p3'.`p3'^-1*n3*`s2m'^-1
        - `p3'.`p3'^-1*n4*`s1m'^-1*`s2m'^-1*`s4m'
        + `p3'.`p3'^-1*n4*`s2m'^-1*`s4m'*`s6m'^-1
        + `p3'.`p3'^-1*n5*`s2m'^-1
        - `p3'.`p3'^-1*n6*M^2*`s2m'^-1*`s6m'
        + `p3'.`p3'^-1*n6*`s2m'^-1
        + `p3'.`p3'^-1*M^2*(n2-1)
        - `p3'.`p3'^-1*(n2-1)*`s1m'^-1
        - 1/2*`p3'.`p3'^-1*`s2m'^-1*num(d)
        + `p3'.`p3'^-1*`s2m'^-1
        );

        redefine i "0";

        elseif ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) && (count(`s6m',1)==1) && 
        (count(`p5'.`p5',1)!=-1) && (count(p5.p5,1)!=0));

*** for the special case (n3=-2, n5=-1) reverse order of n3 and n5
        if ( (count(`p3'.`p3',1)==2) && (count(`p5'.`p5',1)==1) ) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');


        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n5-1)*(
        + `p3'.`p3'^-1*`p5'.`p5'
        - `p3'.`p3'^-1*`p5'.`p5'*n1*M^2*`s1m'
        + `p3'.`p3'^-1*`p5'.`p5'*n1
        - `p3'.`p3'^-1*`p5'.`p5'*n2*M^2*`s2m'
        + `p3'.`p3'^-1*`p5'.`p5'*n2
        + `p3'.`p3'^-1*`p5'.`p5'*n3
        + `p3'.`p3'^-1*`p5'.`p5'*n4*`s1m'^-1*`s4m'
        - `p3'.`p3'^-1*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        + `p3'.`p3'^-1*`p5'.`p5'*n6*M^2*`s6m'
        - `p3'.`p3'^-1*`p5'.`p5'*n6
        - 1/2*`p3'.`p3'^-1*`p5'.`p5'*num(d)
        + `p3'.`p3'^-1*M^2*(n5-1)
        - `p3'.`p3'^-1*(n5-1)*`s6m'^-1
        );

        redefine i "0";

        elseif ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) && (count(`s6m',1)==1) && 
        (count(`p5'.`p5',1)==-1) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)*deno(4-2*n5-n2-n3,-2)*(
        + n2*`s2m'*`s4m'^-1
        - `p3'.`p3'^-1*`p5'.`p5'*n3
        - `p3'.`p3'^-1*n3*M^2
        + `p3'.`p3'^-1*n3*`s6m'^-1
        - `p5'.`p5'*n2*`s2m'
        );

        redefine i "0";

        endif;

        if ( ( (count(`p3',1) < count(`p5',1)) ) && (count(`p3',1)<0) )
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');

* special case (n5=-1):
        if ( (count(`p5'.`p5',1)==1) )
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure

****************************************

* 24.07.98: new realization of "redD4n35m"; not (yet) used 

#procedure redD4n35ma(s1m,s2m,p3,s4m,p5,s6m)

* reduce n3,n5 from <0 to =0
* (rec. rel. implemented for n3)
        if ( count(intd4,1) );
        if ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m'^2)>0) &&
        (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n6-1)*(
        - `p3'.`p3'^-1*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        - `p3'.`p3'^-1*`p5'.`p5'*(n6-1)
        + `p3'.`p3'^-1*n4*`s2m'^-1*`s4m'*`s6m'^-1
        - `p3'.`p3'^-1*n4*`s6m'^-1
        - 2*`p3'.`p3'^-1*n5*`s6m'^-1
        + `p3'.`p3'^-1*M^2*(n6-1)
        - `p3'.`p3'^-1*(n6-1)*`s6m'^-1
        + `p3'.`p3'^-1*`s6m'^-1*num(d)
        );

        elseif ( (match(`s1m'^2*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s6m',1)==1) && (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

*** for this special case (n3=-2, n5=-1) reverse order of n3 and n5
        if ( (count(`p3'.`p3',1)==2) && (count(`p5'.`p5',1)==1) ) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');

        (-1)/(n1-1)*(
        + `p3'.`p3'^-1*`p5'.`p5'*n4*`s1m'^-1*`s4m'
        + 2*`p3'.`p3'^-1*n3*`s1m'^-1
        - `p3'.`p3'^-1*n4*`s1m'^-1*`s2m'^-1*`s4m'
        + `p3'.`p3'^-1*n4*`s1m'^-1
        + 2*`p3'.`p3'^-1*n5*`s1m'^-1
        - 2*`p3'.`p3'^-1*n6*M^2*`s1m'^-1*`s6m'
        + 2*`p3'.`p3'^-1*n6*`s1m'^-1
        + `p3'.`p3'^-1*(n1-1)*`s1m'^-1
        - `p3'.`p3'^-1*(n1-1)*`s2m'^-1
        - 2*`p3'.`p3'^-1*`s1m'^-1*num(d)
        + 2*`p3'.`p3'^-1*`s1m'^-1
        );

        elseif ( (match(`s1m'*`s2m'^2*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s6m',1)==1) && (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n2-1)*(
        + `p3'.`p3'^-1*n1*M^2*`s1m'*`s2m'^-1
        - `p3'.`p3'^-1*n1*`s2m'^-1
        + `p3'.`p3'^-1*n3*`s2m'^-1
        - `p3'.`p3'^-1*n4*`s1m'^-1*`s2m'^-1*`s4m'
        + `p3'.`p3'^-1*n4*`s2m'^-1*`s4m'*`s6m'^-1
        + `p3'.`p3'^-1*n5*`s2m'^-1
        - `p3'.`p3'^-1*n6*M^2*`s2m'^-1*`s6m'
        + `p3'.`p3'^-1*n6*`s2m'^-1
        + `p3'.`p3'^-1*M^2*(n2-1)
        - `p3'.`p3'^-1*(n2-1)*`s1m'^-1
        - 1/2*`p3'.`p3'^-1*`s2m'^-1*num(d)
        + `p3'.`p3'^-1*`s2m'^-1
        );

        endif;

        if (count(`p3',1) < count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure

#procedure redD4n35md(s1m,s2m,p3,s4m,p5,s6m)
        if ( count(intd4,1) );
        if ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) && (count(`s6m',1)==1) && 
        (count(`p5'.`p5',1)!=-1) && (count(`p5'.`p5',1)!=0));

* for this special case (n3=-2, n5=-1) reverse order of n3 and n5
        if ( (count(`p3'.`p3',1)==2) && (count(`p5'.`p5',1)==1) ) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n5-1)*(
        + `p3'.`p3'^-1*`p5'.`p5'
        - `p3'.`p3'^-1*`p5'.`p5'*n1*M^2*`s1m'
        + `p3'.`p3'^-1*`p5'.`p5'*n1
        - `p3'.`p3'^-1*`p5'.`p5'*n2*M^2*`s2m'
        + `p3'.`p3'^-1*`p5'.`p5'*n2
        + `p3'.`p3'^-1*`p5'.`p5'*n3
        + `p3'.`p3'^-1*`p5'.`p5'*n4*`s1m'^-1*`s4m'
        - `p3'.`p3'^-1*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        + `p3'.`p3'^-1*`p5'.`p5'*n6*M^2*`s6m'
        - `p3'.`p3'^-1*`p5'.`p5'*n6
        - 1/2*`p3'.`p3'^-1*`p5'.`p5'*num(d)
        + `p3'.`p3'^-1*M^2*(n5-1)
        - `p3'.`p3'^-1*(n5-1)*`s6m'^-1
        );

        endif;

        if (count(`p3',1) < count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure


#procedure redD4n35me(s1m,s2m,p3,s4m,p5,s6m)
        if ( count(intd4,1) );
        if ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) && (count(`s6m',1)==1) && 
        (count(`p5'.`p5',1)==-1) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)*deno(4-2*n5-n2-n3,-2)*(
        + n2*`s2m'*`s4m'^-1
        - `p3'.`p3'^-1*`p5'.`p5'*n3
        - `p3'.`p3'^-1*n3*M^2
        + `p3'.`p3'^-1*n3*`s6m'^-1
        - `p5'.`p5'*n2*`s2m'
        );

        endif;

        if (count(`p3',1) < count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure

****************************************

#procedure redD4n35(s1m,s2m,p3,s4m,p5,s6m)

* reduce n3,n5 from >1 to =1
* (rec. rel. implemented for n3)

        if ( count(intd4,1) );
        if (match(`s1m'*`s2m'/`p3'.`p3'^2*`s4m'/`p5'.`p5'*`s6m')>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/M^2/(n3-1)*(
        - (n3-1)*`s6m'^-1
        + `p3'.`p3'*`p5'.`p5'*n2*`s2m'
        - `p3'.`p3'*n2*`s2m'*`s4m'^-1
        + `p3'.`p3'*n2
        + 2*`p3'.`p3'*n5
        + `p3'.`p3'*(n3-1)
        - `p3'.`p3'*num(d)
        + `p5'.`p5'*(n3-1)
        );

        redefine i "0";

        endif;

        if (count(`p3',1) > count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure


#procedure redD4n2(s1m,s2m,p3,s4m,p5,s6m)

* reduce n2 from >1 to =1
        if ( count(intd4,1) );
        if (match(`s1m'*`s2m'^2/`p3'.`p3'*`s4m'/`p5'.`p5'*`s6m')>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/M^2/(n2-1)*(
        - 1/2*n3*`s2m'^-1
        - 1/2*n5*`s2m'^-1
        - (n2-1)*`s2m'^-1
        + 1/2*`s2m'^-1*num(d)
        + 1/2*`p3'.`p3'^-1*n3*`s1m'^-1*`s2m'^-1
        - 1/2*`p3'.`p3'^-1*n3*`s2m'^-2
        - 1/2*`p5'.`p5'^-1*n5*`s2m'^-2
        + 1/2*`p5'.`p5'^-1*n5*`s2m'^-1*`s4m'^-1
        );

        redefine i "0";

        endif;

        if (count(`p3',1) > count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure


#procedure redD4n14(s1m,s2m,p3,s4m,p5,s6m)

* reduce n1,n4 from >1 to =1
* (rec. rel. implemented for n1)
        if ( count(intd4,1) );
        if (match(`s1m'^2*`s2m'/`p3'.`p3'*`s4m'/`p5'.`p5'*`s6m')>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/M^2/(n1-1)*(
        + n2*`s1m'^-1*`s2m'*`s4m'^-1
        - n2*`s1m'^-1
        - n4*`s1m'^-1*`s2m'^-1*`s4m'
        + n4*`s1m'^-1
        - (n1-1)*`s2m'^-1
        + (n1-1)*`s4m'^-1
        - (n1-1)*`s6m'^-1
        + `p3'.`p3'*(n1-1)
        - `p5'.`p5'*n2*`s1m'^-1*`s2m'
        + `p5'.`p5'*n4*`s1m'^-1*`s4m'
        );

        redefine i "0";

        endif;

***if (count(`s4m',1) > count(`s1m',1)) 
***  multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
*** 11Mar04: Changed to (see email from M. Faisst)
        if (count(`s4m',1) > count(`s1m',1)) ;
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        redefine i "0" ;
        endif;
        endif;
#endprocedure


#procedure redD4n6(s1m,s2m,p3,s4m,p5,s6m)

* reduce n6 from >1 to =1
        if ( count(intd4,1) );
        if (match(`s1m'*`s2m'/`p3'.`p3'*`s4m'/`p5'.`p5'*`s6m'^2)>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/M^2/(n6-1)*(
        + 1/2*n1*`s1m'*`s2m'^-1*`s6m'^-1
        - 1/2*n1*`s6m'^-1
        - n3*`s6m'^-1
        + 1/2*n4*`s2m'^-1*`s4m'*`s6m'^-1
        - 1/2*n4*`s6m'^-1
        - n5*`s6m'^-1
        - (n6-1)*`s6m'^-1
        + `s6m'^-1*num(d)
        - 1/2*`p3'.`p3'*n1*`s1m'*`s6m'^-1
        - 1/2*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        );

        redefine i "0";

        endif;
        endif;
#endprocedure



#procedure topd4
*
* this is topd4
*
* 20.Jun.02: splitting of diad4 commented
*
        #message this is topd4

************************************************************
*
* Jul. 24th 1998:
* - split expression before the repeat-endrepeat-loop of "redD4n35m"
* - split furthermore the procedure "redD4n35m" into smaller parts 
*   (not extensively tested)
*
************************************************************
        
************************************************************

* treat the scalar products

        #message numerator
        if ( count(intd4,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p6.p6 = 1/s6m - M^2;
        endif;        
        #call ACCU(D4 0)

        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

        if ( count(intd4,1) );                
        id p2 = p1+p3;
        id p4 = p1+p6;
        id p5 = p6-p3;
        endif;        
        #call ACCU(D4 1)

        if ( count(intd4,1) );                
        id  p1.p3 = 1/2 * (   1/s2m - 1/s1m - p3.p3 );
        id  p1.p6 = 1/2 * (   1/s4m - 1/s1m - 1/s6m + M^2 );
        id  p3.p6 = 1/2 * ( - p5.p5 + p3.p3 + 1/s6m - M^2 );
        endif;        
        #call ACCU(D4 2)

        if ( count(intd4,1) );                
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p6.p6 = 1/s6m - M^2;
        endif;        
        #call ACCU(D4 3)

        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

        #call ACCU(D4)

************************************************************
*
* do recursion
*

        #message do recursion

*
* massless indices <0
*
        #do i=1,1
                
                #call redD4n35m(s1m,s2m,p3,s4m,p5,s6m)
                #call ACCU(D4n35m)
                
        #enddo
        
        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)
        
        #do i=1,1
                
                #call redD4n35m(s1m,s2m,p3,s4m,p5,s6m)
                #call ACCU(D4n35m)
                
        #enddo
        
        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)
        
        
*
* massless indices >0
*
        
        #do i=1,1
                
                #call redD4n35(s1m,s2m,p3,s4m,p5,s6m)
                #call redD4n2(s1m,s2m,p3,s4m,p5,s6m)
                #call ACCU(D4n235)
                
        #enddo
        
        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

        #call Conv2exact

        #do i=1,1
                
                #call redD4n14(s1m,s2m,p3,s4m,p5,s6m)
                .sort
                
        #enddo
        
        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)
        
        #call ACCU(D4)
        
        #do i=1,1
                
                #call redD4n35(s1m,s2m,p3,s4m,p5,s6m)
                #call redD4n2(s1m,s2m,p3,s4m,p5,s6m)
                #call ACCU(D4n235)
                
        #enddo
        
        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)
        
        #call Conv2exact
        
        #do i=1,1
                
                #call redD4n6(s1m,s2m,p3,s4m,p5,s6m)
                .sort
                
        #enddo

        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)
        
        #call ACCU(D4)
        
        if ( count(intd4,1) );        
        if (count(s4m,1) < count(s1m,1)) 
        multiply replace_(s1m,s4m,s4m,s1m,p3,p5,p5,p3);
        
        if ( (count(s1m,1)>0) && (count(s4m,1)>0) );
        if (count(p3,1) < count(p5,1)) 
        multiply replace_(s1m,s4m,s4m,s1m,p3,p5,p5,p3);
        endif;
        
        if (count(p5,1)==0) multiply replace_(s1m,s4m,s4m,s1m,p3,p5,p5,p3);
        endif;
        .sort
        
        #call ACCU(D4)
        
************************************************************
        
* identify simple integrals
        
        if ( count(intd4,1) );        
        if ( (count(s6m,1)<=0) && (count(s2m,1)>0) );
        id 1/s6m=p6.p6+M^2; 
        id p6=p4-p1;
        id p5=-p5;
        multiply replace_(p1,p5,p2,p4,p3,p2,p4,p6,p5,p1,
        s1m,s5m,s2m,s4m,s4m,s6m);
        multiply, intbm/intd4;

        elseif ( (count(s2m,1)<=0) );
        id 1/s2m=p2.p2+M^2;
        id p2=p1+p3; 
        id p3=-p3;
        id p4=-p4;
        id p6=-p6;
        multiply replace_(p3,p4,p4,p3,p6,p2,
        s4m,s3m,s6m,s2m);
        multiply, intdm/intd4;

        elseif ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2;
        id p1=p4-p6; 
        id p5=-p5;
        multiply replace_(p2,p4,p3,p5,p4,p6,p5,p2,p6,p3,
        s2m,s4m,s4m,s6m,s6m,s3m);
        
* needed because of topBN1
        id 1/s3m=p3.p3+M^2;
        id 1/s4m=p4.p4+M^2;
        id 1/s6m=p6.p6+M^2;
        
        multiply, intbn1/intd4;

        elseif ( (count(p3,1)==0) );
        id p5=-p5;
        multiply replace_(p1,p3,p2,p4,p4,p2,p6,p1,
        s1m,s3m,s2m,s4m,s4m,s2m,s6m,s1m);
        multiply, inte4/intd4;

        elseif ( (count(s1m,1)==1) && (count(s2m,1)==1)    && (count(p3.p3,1)==-1) && 
        (count(s4m,1)==1) && (count(p5.p5,1)==-1) && (count(s6m,1)==1) );
        id s1m*s2m/p3.p3*s4m/p5.p5*s6m = miD4;
        Multiply int0/intd4;

        else;
        exit "D4: Unknown simpler topology";
        endif;
        endif;
        #call ACCU(D4)
        
        #message - done
#endprocedure






************************************************************

#procedure mltadDM(s1m,s2m,s3m,p4,p5,p6)

* discard massless tadpoles
*
* two massless indices >=0 => D_M = 0
        if ( count(intdm,1) );                
        if ( (count(`p6'.`p6',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;
        if ( (count(`p6'.`p6',1)>=0) && (count(`p5'.`p5',1)>=0) ) discard;
        if ( (count(`p5'.`p5',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;

* two massive indices >=0 => D_M = 0

        if ( (count(`s1m',1)<=0) && (count(`s2m',1)<=0) ) discard;
        if ( (count(`s1m',1)<=0) && (count(`s3m',1)<=0) ) discard;
        if ( (count(`s2m',1)<=0) && (count(`s3m',1)<=0) ) discard;
        endif;        
        .sort

#endprocedure

#procedure redDMn4m(s1m,s2m,s3m,p4,p5,p6)

* reduce massless indices from <0 to =0

* sort: n5,n6 >= n4
        if ( count(intdm,1) );                
        if (count(`p5'.`p5',1) > count(`p4'.`p4',1)) 
        multiply replace_(`p4',`p5',`p5',`p4',`s1m',`s3m',`s3m',`s1m');
        if (count(`p4'.`p4',1) < count(`p6'.`p6',1)) 
        multiply replace_(`p6',`p4',`p4',`p6',`s3m',`s2m',`s2m',`s3m');


*** n1=n2=1:
        if ( (match(`s1m'*`s2m'*`s3m'*`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) );

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)*deno(4-n4-3/2*n1-3/2*n2,-2)*(
        + 3/2*n1*M^2*`s1m'
        - 1/2*n1*`s1m'*`s2m'^-1
        + 1/2*n1*`s1m'*`s3m'^-1
        + 3/2*n2*M^2*`s2m'
        - 1/2*n2*`s1m'^-1*`s2m'
        + 1/2*n2*`s2m'*`s3m'^-1
        + 1/2*`p4'.`p4'^-1*`p5'.`p5'*n4
        + 1/2*`p4'.`p4'^-1*`p6'.`p6'*n4
        + `p4'.`p4'^-1*n4*M^2
        - 1/2*`p4'.`p4'^-1*n4*`s1m'^-1
        - 1/2*`p4'.`p4'^-1*n4*`s2m'^-1
        );

        redefine i "0";

        endif;

        if (match(`s1m'*`s2m'^2*`s3m'*`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)/(n2-1)*(
        - `p4'.`p4'^-1*`p5'.`p5'*n3*`s2m'^-1*`s3m'
        - `p4'.`p4'^-1*`p5'.`p5'*(n2-1)
        + `p4'.`p4'^-1*`p6'.`p6'*n3*`s2m'^-1*`s3m'
        - `p4'.`p4'^-1*n1*M^2*`s1m'*`s2m'^-1
        + `p4'.`p4'^-1*n1*`s2m'^-1
        + `p4'.`p4'^-1*n4*`s2m'^-1
        - `p4'.`p4'^-1*n5*`s2m'^-1
        + `p4'.`p4'^-1*n6*`s2m'^-1
        - 1/2*`p4'.`p4'^-1*`s2m'^-1*num(d)
        + `p4'.`p4'^-1*`s2m'^-1
        );

        redefine i "0";

        endif;

        if ( (count(`p5'.`p5',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;
        if ( (count(`p6'.`p6',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;


        if ( (match(`s1m'^2*`s2m'*`s3m'*`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0) &&
        (count(`s2m',1)==1) );

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)/(n1-1)*(
        + `p4'.`p4'^-1*`p5'.`p5'*n3*`s1m'^-1*`s3m'
        - `p4'.`p4'^-1*`p6'.`p6'*n3*`s1m'^-1*`s3m'
        - `p4'.`p4'^-1*`p6'.`p6'*(n1-1)
        - `p4'.`p4'^-1*n2*M^2*`s1m'^-1*`s2m'
        + `p4'.`p4'^-1*n2*`s1m'^-1
        + `p4'.`p4'^-1*n4*`s1m'^-1
        + `p4'.`p4'^-1*n5*`s1m'^-1
        - `p4'.`p4'^-1*n6*`s1m'^-1
        - 1/2*`p4'.`p4'^-1*`s1m'^-1*num(d)
        + `p4'.`p4'^-1*`s1m'^-1
        );

        redefine i "0";

        endif;

        if ( (count(`p5'.`p5',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;
        if ( (count(`p6'.`p6',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;
        endif;
#endprocedure

#procedure redDMn456(s1m,s2m,s3m,p4,p5,p6)

* reduce massless indices from >1 to =1

* sort: n5 >= n4 >= n6
        if ( count(intdm,1) );                
        if (count(`p5'.`p5',1) > count(`p4'.`p4',1)) 
        multiply replace_(`p4',`p5',`p5',`p4',`s1m',`s3m',`s3m',`s1m');
        if (count(`p5'.`p5',1) > count(`p6'.`p6',1)) 
        multiply replace_(`p6',`p5',`p5',`p6',`s1m',`s2m',`s2m',`s1m');
        if (count(`p4'.`p4',1) > count(`p6'.`p6',1)) 
        multiply replace_(`p6',`p4',`p4',`p6',`s3m',`s2m',`s2m',`s3m');

        if (match(`s1m'*`s2m'*`s3m'/`p4'.`p4'/`p5'.`p5'^2/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 
        (-1)/M^2/(n5-1)*(
        - 1/2*(n5-1)*`s2m'^-1
        - 1/2*(n5-1)*`s3m'^-1
        - 3/2*`p4'.`p4'*`p5'.`p5'*n2*`s2m'
        + 1/2*`p4'.`p4'*(n5-1)
        + 3/2*`p5'.`p5'^2*n2*`s2m'
        + 3/2*`p5'.`p5'^2*n3*`s3m'
        - 3/2*`p5'.`p5'*`p6'.`p6'*n3*`s3m'
        + 1/2*`p5'.`p5'*n2*`s1m'^-1*`s2m'
        - 1/2*`p5'.`p5'*n2*`s2m'*`s3m'^-1
        + 1/2*`p5'.`p5'*n3*`s1m'^-1*`s3m'
        - 1/2*`p5'.`p5'*n3*`s2m'^-1*`s3m'
        + 2*`p5'.`p5'*(n5-1)
        - 1/2*`p5'.`p5'*num(d)
        + 1/2*`p6'.`p6'*(n5-1)
        );

        redefine i "0";

        endif;
        endif;
#endprocedure

#procedure redDMn123(s1m,s2m,s3m,p4,p5,p6)

* reduce massive indices from >1 to =1

* sort: n1 >= n2 >= n3
        if ( count(intdm,1) );                
        if (count(`s1m',1) < count(`s2m',1)) 
        multiply replace_(`p6',`p5',`p5',`p6',`s1m',`s2m',`s2m',`s1m');
        if (count(`s1m',1) < count(`s3m',1)) 
        multiply replace_(`p4',`p5',`p5',`p4',`s1m',`s3m',`s3m',`s1m');
        if (count(`s2m',1) < count(`s3m',1)) 
        multiply replace_(`p6',`p4',`p4',`p6',`s3m',`s2m',`s2m',`s3m');

        if (match(`s1m'^2*`s2m'*`s3m'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 
        (-1)/M^2/(n1-1)*(
        - n4*`s1m'^-1
        + n5*`s1m'^-1
        - n6*`s1m'^-1
        - (n1-1)*`s1m'^-1
        + 1/2*`s1m'^-1*num(d)
        - `p4'.`p4'*n2*`s1m'^-1*`s2m'
        + `p5'.`p5'*n2*`s1m'^-1*`s2m'
        + `p5'.`p5'*n3*`s1m'^-1*`s3m'
        - `p6'.`p6'*n3*`s1m'^-1*`s3m'
        );

        redefine i "0";

        endif;
        endif;
#endprocedure




#procedure topdm

*
* this is topdm
*
        #message this is topdm


************************************************************

* treat the scalar products

        #message numerator

        if ( count(intdm,1) );        
        id p2=p1+p3;
        id p6=p3+p5;
        endif;        
        #call ACCU(DM 1)

        if ( count(intdm,1) );        
        id  p1.p3 = 1/2 * (   p2.p2 - p1.p1 - p3.p3 );
        id  p1.p4 = 1/2 * ( - p6.p6 + p1.p1 + p4.p4 );
        endif;        
        #call ACCU(DM 2)

        if ( count(intdm,1) );                
        id  p1.p5 = 1/2 * (   p3.p3 + p4.p4 - p2.p2 - p6.p6 );
        id  p3.p4 = 1/2 * (   p2.p2 + p6.p6 - p1.p1 - p5.p5 );
        endif;        
        #call ACCU(DM 3)

        if ( count(intdm,1) );                
        id  p3.p5 = 1/2 * (   p6.p6 - p3.p3 - p5.p5 );
        id  p4.p5 = 1/2 * ( - p2.p2 + p4.p4 + p5.p5 );
        endif;        
        #call ACCU(DM 4)

        if ( count(intdm,1) );                
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        endif;        
        #call ACCU(DM 5)

        #call mltadDM(s1m,s2m,s3m,p4,p5,p6)

        #call ACCU(DM)

************************************************************

* do recursion

        #message do recursion

* sort: n5 >= n4 >= n6

        if ( count(intdm,1) );                
        if (count(p5.p5,1) > count(p4.p4,1)) 
        multiply replace_(p4,p5,p5,p4,s1m,s3m,s3m,s1m);
        if (count(p5.p5,1) > count(p6.p6,1)) 
        multiply replace_(p6,p5,p5,p6,s1m,s2m,s2m,s1m);
        if (count(p4.p4,1) > count(p6.p6,1)) 
        multiply replace_(p6,p4,p4,p6,s3m,s2m,s2m,s3m);
        endif;        
        .sort

************************************************************

* massless index <0

        #do i=1,1
                #call redDMn4m(s1m,s2m,s3m,p4,p5,p6)
                .sort
        #enddo

        #call mltadDM(s1m,s2m,s3m,p4,p5,p6)

        #call ACCU(DM)

* massless index >1

        #do i=1,1
                #call redDMn456(s1m,s2m,s3m,p4,p5,p6)
                .sort
        #enddo

        #call mltadDM(s1m,s2m,s3m,p4,p5,p6)

        #call ACCU(DM)

        #do i=1,1
                #call redDMn123(s1m,s2m,s3m,p4,p5,p6)
                .sort
        #enddo

        #call mltadDM(s1m,s2m,s3m,p4,p5,p6)

        #call ACCU(DM)

************************************************************

* identify simple integrals

* sort: n5 >= n4 >= n6
* (this order is needed for the identification of the simple integrals)

        if ( count(intdm,1) );        
        if (count(p5.p5,1) > count(p4.p4,1)) 
        multiply replace_(p4,p5,p5,p4,s1m,s3m,s3m,s1m);
        if (count(p5.p5,1) > count(p6.p6,1)) 
        multiply replace_(p6,p5,p5,p6,s1m,s2m,s2m,s1m);
        if (count(p4.p4,1) > count(p6.p6,1)) 
        multiply replace_(p6,p4,p4,p6,s3m,s2m,s2m,s3m);
        endif;
        .sort

* if one of the massive indices is absent
* sort: n1 <= n2 <= n3

        if ( count(intdm,1) );        
        if (match(s1m*s2m*s3m)==0);
        if (count(s1m,1) > count(s2m,1)) 
        multiply replace_(p6,p5,p5,p6,s1m,s2m,s2m,s1m);
        if (count(s1m,1) > count(s3m,1)) 
        multiply replace_(p4,p5,p5,p4,s1m,s3m,s3m,s1m);
        if (count(s2m,1) > count(s3m,1)) 
        multiply replace_(p6,p4,p4,p6,s3m,s2m,s2m,s3m);
        endif;
        endif;        
        .sort

        #call ACCU(DM)

        if ( count(intdm,1) );        
        if ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2; 
        id p1=p2-p3;
        id p5=-p5;
        multiply replace_(p2,p6,p3,p4,p4,p3,p5,p1,p6,p5,s2m,s6m,s3m,s4m);
        multiply, intbn2/intdm;

        elseif ( (count(s1m,1) >0)   && (count(s2m,1) >0)   && (count(s3m,1) >0) && 
        (count(p4.p4,1) <0) && (count(p5.p5,1) <0) && (count(p6.p6,1)==0) ); 

* use symmetry for E_3:

        if (count(s1m,1) > count(s3m,1))
        multiply replace_(p1,p3,p3,p1,s1m,s3m,s3m,s1m);
        multiply, inte3/intdm;

        elseif ( (count(s1m,1)==1)   && (count(s2m,1)==1)   && (count(s3m,1)==1) && 
        (count(p4.p4,1)==-1) && (count(p5.p5,1)==-1) && (count(p6.p6,1)==-1));
        id s1m*s2m*s3m/p4.p4/p5.p5/p6.p6 = miDM*int0/intdm;

        else;
        exit "DM: Unknown simpler topology";
        endif;
        endif;
        #call ACCU(DM)

        #message - done
        
#endprocedure









************************************************************

#procedure mltadDN(s1m,s2m,p3,p4,p5,p6)

* discard massless tadpoles
        if ( count(intdn,1) );
        if ( (count(`s1m',1)<=0) && (count(`s2m',1)<=0) ) discard;
        if ( (count(`s1m',1)<=0) && 
        ( (count(`p3'.`p3',1)>=0) || (count(`p4'.`p4',1)>=0) ||
        (count(`p5'.`p5',1)>=0) || (count(`p6'.`p6',1)>=0) ) ) discard;
        if ( (count(`s2m',1)<=0) && 
        ( (count(`p3'.`p3',1)>=0) || (count(`p4'.`p4',1)>=0) ||
        (count(`p5'.`p5',1)>=0) || (count(`p6'.`p6',1)>=0) ) ) discard;
        endif;        
        .sort

#endprocedure


#procedure redDNn6(s1m,s2m,p3,p4,p5,p6)

* reduce massless indices from >1 to =1
        if ( count(intdn,1) );
        if (match(`s1m'*`s2m'/`p3'.`p3'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

* sort: n6 >=  n3,n4,n5

        if (count(`p6'.`p6',1) > count(`p5'.`p5',1)) 
        multiply replace_(`p6',`p5',`p5',`p6',`s1m',`s2m',`s2m',`s1m');
        if (count(`p6'.`p6',1) > count(`p4'.`p4',1)) 
        multiply replace_(`p6',`p4',`p4',`p6',`p3',`p5',`p5',`p3');
        if (count(`p6'.`p6',1) > count(`p3'.`p3',1)) 
        multiply replace_(`p6',`p3',`p3',`p6',`p4',`p5',`p5',`p4');
        endif;

        if (match(`s1m'*`s2m'/`p3'.`p3'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6'^2)>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?
        /`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2/`p3'.`p3'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)/M^2/(n6-1)*(
        - 1/2*(n6-1)*`s1m'^-1
        - 1/2*(n6-1)*`s2m'^-1
        + 1/2*`p3'.`p3'*`p6'.`p6'*n2*`s2m'
        + 1/2*`p3'.`p3'*(n6-1)
        + 1/2*`p4'.`p4'*`p6'.`p6'*n1*`s1m'
        + 1/2*`p4'.`p4'*(n6-1)
        - 1/2*`p5'.`p5'*`p6'.`p6'*n1*`s1m'
        - 1/2*`p5'.`p5'*`p6'.`p6'*n2*`s2m'
        + 1/2*`p6'.`p6'*n3
        + 1/2*`p6'.`p6'*n4
        - 1/2*`p6'.`p6'*n5
        + 1/2*`p6'.`p6'*(n6-1)
        - 1/4*`p6'.`p6'*num(d)
        );

        redefine i "0";
        redefine ii "0";
        endif;
        endif;
#endprocedure

#procedure redDNn1(s1m,s2m,p3,p4,p5,p6)

* reduce massive indices from >1 to =1
        if ( count(intdn,1) );
        if (match(`s1m'*`s2m'/`p3'.`p3'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

* sort: n1 >= n2

        if (count(`s1m',1) < count(`s2m',1)) 
        multiply replace_(`p6',`p5',`p5',`p6',`s1m',`s2m',`s2m',`s1m');
        endif;

        if (match(`s1m'^2*`s2m'/`p3'.`p3'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?
        /`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2/`p3'.`p3'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)/M^2/(n1-1)*(
        + 1/2*n3*`s1m'^-1
        - 3/2*n4*`s1m'^-1
        - 1/2*n5*`s1m'^-1
        - 1/2*n6*`s1m'^-1
        - (n1-1)*`s1m'^-1
        + 3/4*`s1m'^-1*num(d)
        + 1/2*`p3'.`p3'*`p6'.`p6'^-1*n6*`s1m'^-1
        + 1/2*`p3'.`p3'*n2*`s1m'^-1*`s2m'
        - 1/2*`p4'.`p4'*`p6'.`p6'^-1*n6*`s1m'^-1
        - 1/2*`p4'.`p4'*(n1-1)
        - 1/2*`p5'.`p5'*n2*`s1m'^-1*`s2m'
        + 1/2*`p5'.`p5'*(n1-1)
        - 1/2*`p6'.`p6'^-1*n6*`s1m'^-2
        + 1/2*`p6'.`p6'^-1*n6*`s1m'^-1*`s2m'^-1
        );

        redefine i "0";
        redefine ii "0" ;
        endif;
        endif;
#endprocedure



#procedure topdn
*
* this is topdn
*
        #message this is topdn


************************************************************

* treat the scalar products     CHECK !!!!!!!!!!

        #message numerator
        if ( count(intdn,1) );
        id p3=p5-p2;
        id p4=p2+p6;
        endif;        
        #call ACCU(DN 1)

        if ( count(intdn,1) );        
        id  p1.p2 = 1/2 * ( p4.p4 + p3.p3 - p5.p5 - p6.p6 );
        id  p1.p5 = 1/2 * ( p4.p4 - p5.p5 - p1.p1 );
        endif;        
        #call ACCU(DN 2)

        if ( count(intdn,1) );        
        id  p1.p6 = 1/2 * (-p3.p3 + p6.p6 + p1.p1 );
        id  p2.p5 = 1/2 * (-p3.p3 + p5.p5 + p2.p2 );
        endif;        
        #call ACCU(DN 3)

        if ( count(intdn,1) );        
        id  p2.p6 = 1/2 * ( p4.p4 - p6.p6 - p2.p2 );
        id  p5.p6 = 1/2 * ( p3.p3 + p4.p4 - p2.p2 - p1.p1);
        endif;        
        #call ACCU(DN 4)

        if ( count(intdn,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        endif;        
        .sort

        #call mltadDN(s1m,s2m,p3,p4,p5,p6)

        #call ACCU(DN)


************************************************************

* do recursion

        #message do recursion

        #define ii "0"
        #do i=1,1
                #call redDNn6(s1m,s2m,p3,p4,p5,p6)
                .sort
        #enddo
        #undefine ii

        #call mltadDN(s1m,s2m,p3,p4,p5,p6)

        #call ACCU(DN)

        #do i=1,1

                #do ii=1,1
                        #call redDNn1(s1m,s2m,p3,p4,p5,p6)
                        .sort
                #enddo 

                #call mltadDN(s1m,s2m,p3,p4,p5,p6)

                #do ii=1,1
                        #call redDNn6(s1m,s2m,p3,p4,p5,p6)
                        .sort
                #enddo

*                 id n=num(4-2*ep);
*                 id acc(x1?)*acc(x2?)=acc(x1*x2);
                #call ACCU(DN)
                
        #enddo
        
        #call mltadDN(s1m,s2m,p3,p4,p5,p6)
        
        #call ACCU(DN)
        
************************************************************
        
* identify simple integrals
        
        if ( count(intdn,1) );
        if (count(s2m,1)<=0) multiply replace_(p6,p5,p5,p6,s1m,s2m,s2m,s1m);
        if (count(p4.p4,1)>=0) multiply replace_(p4,p3,p3,p4,s1m,s2m,s2m,s1m);
        if (count(p5.p5,1)>=0) multiply replace_(p5,p3,p3,p5,p4,p6,p6,p4);
        if (count(p6.p6,1)>=0) multiply replace_(p6,p3,p3,p6,p4,p5,p5,p4);
        endif;
        .sort
        
        if ( count(intdn,1) );
        if ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2; 
        id p1=p4-p5;
        id p3=-p3;
        id p5=-p5;
        multiply replace_(p2,p6,p4,p2,p5,p4,p6,p1,
        s2m,s6m);
        multiply, intm3/intdn;

        elseif ( (count(p3.p3,1)>=0) );
        id p3=p5-p2;
        id p6=-p6;
        multiply replace_(p1,p3,p4,p5,p2,p4,p5,p2,p6,p1,
        s1m,s3m,s2m,s4m);
        multiply, intbn1/intdn;

        elseif ( (count(s1m,1)==1)   && (count(s2m,1)==1)   && (count(p3.p3,1)==-1) && 
        (count(p4.p4,1)==-1) && (count(p5.p5,1)==-1) && (count(p6.p6,1)==-1));
        id s1m*s2m/p3.p3/p4.p4/p5.p5/p6.p6 = miDN*int0/intdn;

        else;
        exit "DN: Unknown simpler topology";
        endif;
        endif;
        #call ACCU(DN)
        
        #message - done
        
#endprocedure        









#procedure tope4

*
* this is tope4
*
        #message this is tope4

* As this topology results form other ones by shrinking one line
* in principle no scalar products may appear. However, one has to
* take care, that, e.g., in the topology D5 the index of line 2 is 
* really reduced to 0 and not to negative values. Otherwise unwanted 
* scalar products appear.
* Remark: It can be avoided to reduce D5 to E4 ...
        if ( count(inte4,1) );
        if ( count(p1.p5,1,p1.p4,1,p3.p4,1,p3.p5,1)>0 ) multiply 1/(1-1);
        endif;        
        .sort
        
************************************************************
*
* treat the scalar products
*
        if ( count(inte4,1) );
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        endif;        
        .sort

        if ( count(inte4,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        endif;        
        #call ACCU(E4_0)

        if ( count(inte4,1) );        
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        endif;        
        .sort

        #message numerator
        if ( count(inte4,1) );
        id  p1.p2 = 1/2 * (-1/s3m + 1/s1m + 1/s2m-M^2);
        id  p1.p3 = 1/2 * ( 1/s2m - 1/s1m - 1/s3m+M^2);
        id  p2.p3 = 1/2 * (-1/s1m + 1/s2m + 1/s3m-M^2);
        endif;        
        #call ACCU(E4_1)
        if ( count(inte4,1) );
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        endif;        
        .sort

        if ( count(inte4,1) );        
        id  p2.p4 = 1/2 * (-p5.p5 + 1/s2m + 1/s4m-2*M^2);
        id  p2.p5 = 1/2 * ( 1/s4m - 1/s2m - p5.p5);
        id  p4.p5 = 1/2 * (-1/s2m + 1/s4m + p5.p5);
        endif;        
        #call ACCU(E4_2)

        if ( count(inte4,1) );        
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        endif;        
        .sort

        if ( count(inte4,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        endif;        
        .sort


************************************************************

* do recursion

        #message do recursion

        if ( count(inte4,1) );
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        if ( count(s1m,1) < count(s3m,1) ) multiply replace_(s1m,s3m,s3m,s1m);
        endif;        
        .sort

        #do i=1,1
                if ( count(inte4,1) );
                if ( match(s1m*s2m*s3m*s4m/p5.p5) > 0 );
                id s1m^n1? * s2m^n2? * s3m^n3? * s4m^n4? / p5.p5^n5? =
                s1m^n1 * s2m^n2 * s3m^n3 * s4m^n4 * (1/p5.p5)^n5 *
                deno(4-2*n5-n4,-2) * n4 * s4m * ( p5.p5 - 1/s2m )
                ;

                redefine i "0";
                
                endif;
                
                if (count(s4m,1)<=0) discard;
                if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
                endif;        
                #call ACCU(E3)
                
        #enddo

        .sort
************************************************************

* identify simple integrals

        if ( count(inte4,1) );
        if ( (count(s1m,1) >0) && (count(s2m,1)<=0) && (count(s3m,1) >0) && 
        (count(s4m,1) >0) );
        id 1/s2m=p2.p2+M^2; 
        id p1=-p1;
        multiply replace_(p1,p3,p2,p1,p3,p6,s1m,s3m,s3m,s6m);
        multiply, intbn1/inte4;

        elseif ( (count(s1m,1) >0) && (count(s2m,1) >0) && (count(s3m,1) >0) && 
        (count(s4m,1) >0) && (count(p5.p5,1)>=0) );
        id p5=p4-p2;
        id p2=-p2;
        id p3=-p3;
        id p4=-p4;
        multiply replace_(p1,p5,p3,p1,s1m,s5m,s3m,s1m);
        multiply, intm5/inte4;

        elseif ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2;
        id p1=p2-p3;
        id p5=-p5;
        multiply replace_(p2,p4,p3,p5,p4,p6,p5,p1,s2m,s4m,s3m,s5m,s4m,s6m);
        multiply, intbm/inte4;

        elseif ( (count(s3m,1)<=0) );
        id 1/s3m=p3.p3+M^2;
        id p3=p2-p1;
        id p1=-p1;
        id p5=-p5;
        multiply replace_(p1,p5,p2,p4,p4,p6,p5,p1,s1m,s5m,s2m,s4m,s4m,s6m);
        multiply, intbm/inte4;

        else;
        exit "E4: Unknown simpler topology";
        endif;
        endif;
        #call ACCU(topE4)

        #message - done
        
#endprocedure




************************************************************

#procedure mltadE3(s1m,s2m,s3m,p4,p5)
        
* discard massless tadpoles
        if ( count(inte3,1) );
        if ( (count(`s1m',1)<=0) && (count(`s2m',1)<=0) ) discard;
        if ( (count(`s3m',1)<=0) && (count(`s2m',1)<=0) ) discard;
        if ( (count(`p4'.`p4',1)>=0) ) discard;
        if ( (count(`p5'.`p5',1)>=0) ) discard;
        endif;        
        .sort

#endprocedure

#procedure redE3n2(s1m,s2m,s3m,p4,p5)

* reduce n2>0 to n2=1

        if ( count(inte3,1) );
        if (match(`s1m'*`s2m'^2*`s3m'/`p4'.`p4'/`p5'.`p5')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5 * 
        (-1)/M^2/(n2-1)*(
        + n1*M^2*`s1m'*`s2m'^-1
        - n1*`s2m'^-1
        + n3*M^2*`s2m'^-1*`s3m'
        - n3*`s2m'^-1
        - n4*`s2m'^-1
        - n5*`s2m'^-1
        - (n2-1)*`s2m'^-1
        + 3/2*`s2m'^-1*num(d)
        );

        redefine i "0";

        endif;
        endif;
#endprocedure

#procedure redE3n45(s1m,s2m,s3m,p4,p5)

* reduce massless lines 4 and 5 to n4=1, n5=1;

* sort: n4 > n5
        if ( count(inte3,1) );
        if ( count(`p4'.`p4',1) > count(`p5'.`p5',1) )
        multiply replace_(`p4',`p5',`p5',`p4');
        endif;
        .sort

        if ( count(inte3,1) );        
        if (match(`s1m'*`s2m'*`s3m'/`p4'.`p4'^2/`p5'.`p5')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5 * 
        (-1)/M^2/(n4-1)*(
        - (n4-1)*`s2m'^-1
        + 2*`p4'.`p4'*n5
        + `p4'.`p4'*(n4-1)
        - `p4'.`p4'*num(d)
        + `p5'.`p5'*(n4-1)
        );

        redefine i "0";

        endif;
        endif;
#endprocedure

#procedure redE3n13(s1m,s2m,s3m,p4,p5)

* reduce lines 1 and 3 to n1=1, n3=1;

* sort: n3 > n1
        if ( count(inte3,1) );
        if ( count(`s1m',1) > count(`s3m',1) )
        multiply replace_(`s1m',`s3m',`s3m',`s1m');
        endif;
        .sort

        if ( count(inte3,1) );        
        if (match(`s1m'*`s2m'*`s3m'^2/`p4'.`p4'/`p5'.`p5')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5 * 
        (-1)/M^2/(n3-1)*(
        + 2/3*n1*`s1m'*`s2m'^-1*`s3m'^-1
        - 2/3*n1*`s1m'*`s3m'^-2
        + 1/3*(n3-1)*`s1m'^-1
        - 1/3*(n3-1)*`s2m'^-1
        - (n3-1)*`s3m'^-1
        + 1/3*`s3m'^-1*num(d)
        );

        redefine i "0";

        endif;
        endif;
#endprocedure



#procedure tope3
*
* this is topE3
*
        #message this is topE3

* As this topology results form other ones by shrinking one line
* in principle no scalar products may appear. However, one has to
* take care, that, e.g., in the topology DM the index of line 6 is 
* really reduced to 0 and not to negative values. Otherwise unwanted 
* scalar products appear.


************************************************************

* treat the scalar products

        #message numerator

        if ( count(inte3,1) );
        if ( count(p1.p5,1,p1.p4,1,p3.p4,1,p3.p5,1)>0 ) multiply 1/(1-1);
        endif;
        .sort

        if ( count(inte3,1) );
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        endif;        
        #call ACCU(E3_0)

        #call mltadE3(s1m,s2m,s3m,p5,p5)

        if ( count(inte3,1) );        
        id  p1.p2 = 1/2 * (-1/s3m + 1/s1m + p2.p2);
        id  p1.p3 = 1/2 * ( 1/s2m - 1/s1m - p3.p3);
        id  p2.p3 = 1/2 * (-1/s1m + 1/s2m + p3.p3);
        endif;        
        #call ACCU(E3_1)

        if ( count(inte3,1) );        
        id  p2.p4 = 1/2 * (-p5.p5 + p2.p2 + p4.p4);
        id  p2.p5 = 1/2 * ( p4.p4 - p2.p2 - p5.p5);
        id  p4.p5 = 1/2 * (-p2.p2 + p4.p4 + p5.p5);
        endif;        
        #call ACCU(E3_2)

        if ( count(inte3,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        endif;        
        #call ACCU(E3_3)

        #call mltadE3(s1m,s2m,s3m,p5,p5)


************************************************************

* do recursion

        #message do recursion

* use symmetry: n1 < n3
*               n5 < n4

        if ( count(inte3,1) );        
        if ( count(s1m,1) > count(s3m,1) ) multiply replace_(s1m,s3m,s3m,s1m);
        if ( count(p4.p4,1) > count(p5.p5,1) ) multiply replace_(p4,p5,p5,p4);
        endif;        
        .sort

        #do i=1,1
                #call redE3n2(s1m,s2m,s3m,p4,p5)
                #call ACCU(E3)
                .sort
        #enddo

        #call mltadE3(s1m,s2m,s3m,p5,p5)

        #do i=1,1
                #call redE3n45(s1m,s2m,s3m,p4,p5)
                #call ACCU(E3)
                .sort
        #enddo

        #call mltadE3(s1m,s2m,s3m,p5,p5)

        #call ACCU(E3)

        #do i=1,1
                #call redE3n13(s1m,s2m,s3m,p4,p5)
                #call ACCU(E3)
                .sort
        #enddo

        #call mltadE3(s1m,s2m,s3m,p5,p5)

        #call ACCU(E3)
        .sort

************************************************************

* identify simple integrals

        if ( count(inte3,1) );
        if ( (count(s2m,1)<=0) );
        id 1/s2m=p2.p2+M^2; 
        id p2=p1+p3;
        id p1=-p1;
        multiply replace_(p1,p4,p3,p6,p4,p3,s1m,s4m,s3m,s6m);
        multiply, intbn2/inte3;

        elseif ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2;
        id p1=p2-p3;
        id p2=-p2;
        id p3=-p3;
        id p4=-p4;
        id p5=-p5;
        multiply replace_(p2,p6,p3,p5,p5,p1,s2m,s6m,s3m,s5m);
        multiply, intbm1/inte3;

        elseif ( (count(s3m,1)<=0) );
        id 1/s3m=p3.p3+M^2;
        id p3=p2-p1;
        id p2=-p2;
        id p4=-p4;
        id p5=-p5;
        multiply replace_(p1,p5,p2,p6,p5,p1,s1m,s5m,s2m,s6m);
        multiply, intbm1/inte3;

        elseif ( (count(s1m,1)==1)   && (count(s2m,1)==1)   && (count(s3m,1)==1) && 
        (count(p4.p4,1)==-1) && (count(p5.p5,1)==-1) );
        id s1m*s2m*s3m/p4.p4/p5.p5 = M^2 * miE3;
        Multiply int0/inte3;

        else;
        multiply 1/(1-1);
        endif;
        endif;

        #call ACCU(E3)

        #message - done
        
#endprocedure        



*************************************************
*
*      Recursion for simpler integrals
* 
*************************************************

#procedure symBNnom(p1,p2,p3,p4,p5,p6,x3,x4,x5,x6)
*
* sort: n6>=n3,n4,n5
*
        if ( count(intbn,1) );                
        if ( (count(`x3',1) > count(`x6',1)) && (count(`x3',1) > count(`x4',1))
        && (count(`x3',1) > count(`x5',1)) );
        id `p2'=-`p2';
        multiply replace_(`x3',`x6',`x6',`x4',`x4',`x5',`x5',`x3',
        `p3',`p6',`p6',`p4',`p4',`p5',`p5',`p3',
        `p2',`p1',`p1',`p2');
        endif;
        if ( (count(`x4',1) > count(`x6',1)) && (count(`x4',1) > count(`x5',1))
        && (count(`x4',1) > count(`x3',1)) );
        id `p1'=-`p1';
        multiply replace_(`x4',`x6',`x6',`x3',`x3',`x5',`x5',`x4',
        `p4',`p6',`p6',`p3',`p3',`p5',`p5',`p4',
        `p2',`p1',`p1',`p2');
        endif;
        if ( (count(`x5',1) > count(`x6',1)) && (count(`x5',1) > count(`x4',1))
        && (count(`x5',1) > count(`x3',1)) );
        id `p1'=-`p1';
        id `p2'=-`p2';
        multiply replace_(`x5',`x6',`x6',`x5',`x3',`x4',`x4',`x3',
        `p5',`p6',`p6',`p5',`p3',`p4',`p4',`p3');
        endif;
*
* sort: n4>n3 or, if n3=n4: n1>=n2
*
        if (count(`x3',1) > count(`x4',1));
        id `p1'=-`p1';
        id `p2'=-`p2';
        multiply replace_(`x3',`x4',`x4',`x3',
        `p3',`p4',`p4',`p3',`p2',`p1',`p1',`p2');
        endif;
        if ( (count(`x3',1) == count(`x4',1)) 
        && (count(`p1'.`p1',1) > count(`p2'.`p2',1)) );
        id `p1'=-`p1';
        id `p2'=-`p2';
        multiply replace_(`x3',`x4',`x4',`x3',
        `p3',`p4',`p4',`p3',`p2',`p1',`p1',`p2');
        endif;
        endif;
#endprocedure


#procedure symBN (p1,p2,x3,x4,x5,x6)
*
* sort: n6>=n3,n4,n5
*       n4>=n3
*
        if ( count(intbn,1) );        
        if (count(`x3',1) > count(`x6',1)) 
        multiply replace_(`x3',`x6',`x6',`x3',`x4',`x5',`x5',`x4');
        if (count(`x4',1) > count(`x6',1)) 
        multiply replace_(`x4',`x6',`x6',`x4',`x3',`x5',`x5',`x3');
        if (count(`x3',1) > count(`x4',1)) 
        multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
        if (count(`x5',1) > count(`x6',1)) 
        multiply replace_(`x5',`x6',`x6',`x5',`p1',`p2',`p2',`p1');
        endif;        
#endprocedure

#procedure BNtoBM (p1,p2,p3,p4,p5,p6,x3,x4,x5,x6)
        
* Changes the notation of the BM's which result when reducing BN.
*
* n6>=1 (after "symBN.prc" is used)
* There are several possibilities for the other n's:
*
* n3=0: 
*      n5!=0 (otherwise the result is ==0)
*      n4=0 or n4!=0 possible
* n4=0:
*      n5!=0 (otherwise the result is ==0)
*      n3=0  (because of "symBN.prc")
*      (-> see case n3=0)
* n5=0:
*      n3!=0 and n4!=0 (otherwise the result is ==0)
*
* So there are two different cases needed: n3=0 or n5=0.
        if ( count(intbn,1) );
        id `x3'^n3?neg_=(`p3'.`p3'+M^2)^(-n3);
        id `x4'^n4?neg_=(`p4'.`p4'+M^2)^(-n4);
        id `x5'^n5?neg_=(`p5'.`p5'+M^2)^(-n5);
        id `x6'^n6?neg_=(`p6'.`p6'+M^2)^(-n6);
        endif;        
        .sort
        
        if ( count(intbn,1) );        
        if((count(`x3',1) = 0) &&
        (count(`x4',1) >= 0)&&(count(`x5',1) >= 0)&&(count(`x6',1) >= 0) 
        ); 
        id `p3'=`p6'-`p1';
        multiply replace_(`p2',`p3',`x6',`x5',`x5',`x4',`x4',`x6',
        `p6',`p5',`p5',`p4',`p4',`p6');
        id `p1'=-`p1';
        elseif((count(`x5',1) = 0)&&
        (count(`x3',1) >= 0)&&(count(`x4',1) >= 0)&&(count(`x6',1) >= 0) 
        ); 
        id `p5'=`p2'+`p3';
        multiply replace_(`p1',`p3',`p2',`p1',`x3',`x5',
        `p3',`p5');
        endif;
        endif;
#endprocedure;



#procedure nomBN
*
* Decomposition of the numerator for type BN
*
        
        #call symBNnom(p1,p2,p3,p4,p5,p6,x3,x4,x5,x6)
        .sort
        if ( count(intbn,1) );
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;        
        .sort
        
        #do i=1,10
                
                #message pi.pi
                if ( count(intbn,1) );                
                id,once  p3.p3 = 1/x3 - M^2;
                id,once  p4.p4 = 1/x4 - M^2;
                id,once  p5.p5 = 1/x5 - M^2;
                id,once  p6.p6 = 1/x6 - M^2;
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;                
                
                #call ACCU(pi.pi)
                
        #enddo

        if ( count(intbn,1) );        
        id  p3.p3 = 1/x3 - M^2;
        id  p4.p4 = 1/x4 - M^2;
        id  p5.p5 = 1/x5 - M^2;
        id  p6.p6 = 1/x6 - M^2;
        
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        
        #call ACCU(pi.pi)

        #do i=1,10

                #message p1

                if ( count(intbn,1) );        
                id,once  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
                id,once  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
                id,once  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
                id,once  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;        
                
                #call ACCU(p1)
                
        #enddo

        if ( count(intbn,1) );
        id  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
        id  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
        id  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
        id  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );

        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;

        #call ACCU(p1)

        #do i=1,10

                #message p1,p2

                if ( count(intbn,1) );        
                id,once  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
                id,once  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
                id,once  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
                id,once  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;        
                
                #call ACCU(p1 p2)
                
        #enddo
        
        if ( count(intbn,1) );
        id  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
        id  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
        id  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
        id  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
        
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        
        #call ACCU(p1 p2)
        
        #do i=1,10
                
                #message p2,p3
                
                if ( count(intbn,1) );        
                id,once  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
                id,once  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
                id,once  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
                id,once  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;
                #call ACCU(p2 p3)
                
        #enddo
        
        if ( count(intbn,1) );
        id  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
        id  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
        id  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
        id  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
        
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        
        #call ACCU(p2 p3)
        
        #do i=1,10
                
                #message p4,p5,p6
                
                if ( count(intbn,1) );
                id,once  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
                id,once  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
                id,once  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;
                
                #call ACCU(p4 p5 p6)
                
        #enddo
        
        if ( count(intbn,1) );
        id  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
        id  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
        id  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
        
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        
        #call ACCU(p4 p5 p6)
        
        #call symBN(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );
        if ( (count(x3,1)>0)  && (count(x4,1)>0) 
        && (count(x5,1)<=0) && (count(x6,1)<=0)
        && ( (count(p1.p1,1)>=0)    
        || (count(p2.p2,1)>=0) )
        ) discard;
        endif;
        .sort
        
        #message numerator decomposition done (BN)

#endprocedure






#procedure redBNn135 (p1,p2,x3,x4,x5,x6)
        if ( count(intbn,1) );
        repeat;
                if ( (count(`p1'.`p1',1) < 0)
                && (count(`x3',1) >= 1)  &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *deno(4-2*n1-n4-n6,-2)
                *(
                n4*`x4'*( `p1'.`p1' - 1/`x5' )
                +n6*`x6'*( `p1'.`p1' - 1/`x3' )
                )
                ;
                endif;
        endrepeat;
        endif;
        #call Conv2exact
#endprocedure

#procedure redBNn146 (p1,p2,x3,x4,x5,x6)
        
* sort: n6 < n4 < n3,n5
***if (count(`p1'.`p1',1) == 0)
***                 multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
        if ( count(intbn,1) );
        repeat;
                if ( (count(`p1'.`p1',1) < 0)
                && (count(`x3',1) > 0)   &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) > 0) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *deno(4-2*n1-n5-n3,-2)
                *(
                n5*`x5'*( `p1'.`p1' - 1/`x4' )
                +n3*`x3'*( `p1'.`p1' - 1/`x6' )
                )
                ;
                endif;
        endrepeat;
        endif;
#endprocedure

#procedure redBNn1p (p1,p2,x3,x4,x5,x6)

        #do i=1,1
                if ( count(intbn,1) );
                if ( (count(`p1'.`p1',1) > 0)
                && (count(`x3',1) >= 1)  &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *2*M^2*deno(3*4-2*(n1+n2+n3+n4+n5+n6),-6)
                *(
                -2*nom(4,-2)+4*(n1+1)+n3+n4+n5+n6
                -(n3*`x3'/`x6' + n4*`x4'/`x5' + n5*`x5'/`x4' + n6*`x6'/`x3')
                )/`p1'.`p1'
                ;
                
                redefine i "0";
                
                endif;
                endif;        
                .sort
                #call Conv2exact
        #enddo
        
#endprocedure

#procedure redBNn236 (p1,p2,x3,x4,x5,x6)
        if ( count(intbn,1) );
        repeat;
                
*
* sort: n2 < n1
*
                
                if ( (count(`p2'.`p2',1) < 0)
                && (count(`x3',1) > 0)   &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) > 0) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *deno(4-2*n2-n4-n5,-2)
                *(
                n4*`x4'*( `p2'.`p2' - 1/`x6' )
                +n5*`x5'*( `p2'.`p2' - 1/`x3' )
                )
                ;
                endif;
                
        endrepeat;
        endif;
        #call Conv2exact
#endprocedure

#procedure redBNn245 (p1,p2,x3,x4,x5,x6)
        if ( count(intbn,1) );
        repeat;
                if ( (count(`p2'.`p2',1) < 0)
                && (count(`x3',1) >= 1)  &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *deno(4-2*n2-n3-n6,-2)
                *(
                n6*`x6'*( `p2'.`p2' - 1/`x4' )
                +n3*`x3'*( `p2'.`p2' - 1/`x5' )
                )
                ;
                endif;
                
        endrepeat;
        endif;
        #call Conv2exact
#endprocedure

#procedure redBNn2p (p1,p2,x3,x4,x5,x6)
        
        #do i=1,1
                if ( count(intbn,1) );
                if ( (count(`p2'.`p2',1) > 0)
                && (count(`x3',1) >= 1)  &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *2*M^2*deno(3*4-2*(n1+n2+n3+n4+n5+n6),-6)
                *(
                -2*nom(4,-2)+4*(n2+1)+n3+n4+n5+n6
                -(n4*`x4'/`x6' + n3*`x3'/`x5' + n5*`x5'/`x3' + n6*`x6'/`x4')
                )/`p2'.`p2'
                ;

                redefine i "0";

                endif;
                endif;        
                .sort
                #call Conv2exact
        #enddo

#endprocedure

#procedure redBNn342 (p1,p2,x3,x4,x5,x6)
        
        if ( count(intbn,1) );        
        repeat;
                
                if((count(`x3',1) >= 1)      &&  (count(`x4',1) >= 1) &&
                (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1)  );
                
* sort: n4 >= n3
                
                if (count(`x3',1) > count(`x4',1)) 
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                
                
* Now do the reduction via eq. (N1) resp. (5)
                
                if ( (count(`x4',1) > 1) && (count(`p2'.`p2',1) < 0) ); 
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^(n4-1) *(-1) * `x5'^n5  * `x6'^n6
                *1/2/(n4-1)/M^2
                *(
                +nom(4-2*(n4-1)-n1-n6,-2)
                +n1/`p1'.`p1'*(1/`x5'-1/`x4')
                +n6*`x6'     *(`p2'.`p2'-1/`x4'+2*M^2)
                )
                ;
                
                endif;
                endif;
                
        endrepeat;
        endif;
        #call Conv2exact
#endprocedure

#procedure redBNn3456 (p1,p2,x3,x4,x5,x6)

* do-enddo loop or not?
* example: reduction of BN(1,2,1,3,2,3)
*          with do-enddo loop: ~ 80s
*          without: ~ 375s

        #do i = 1,1
                
                if ( count(intbn,1) );
                if (   (count(`p1'.`p1',1) == 0) && (count(`p2'.`p2',1) == 0)
                && (count(`x3',1) >= 3)      && (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)      && (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *(
                ( -2*nom(4,-2) + 4*(n3-1) )/4/M^2/(n3-1) /x3
                -1/4/M^2/(n3-2)/(n3-1) * dala/x3^2
                )
                ;
                
                redefine i "0";
                
                endif;
                endif;        
                .sort
                #call Conv2exact
        #enddo
        
#endprocedure

#procedure redBNn5 (p1,p2,x3,x4,x5,x6)

        #do i=1,1
                
* sort: n6>=n3,n4,n5>=1
*       n4>=n3
                if ( count(intbn,1) );     
                if((count(`x3',1) > 0)  &&  (count(`x4',1) > 0) &&
                (count(`x5',1) > 0)  &&  (count(`x6',1) > 0) );
                
                if (count(`x3',1) > count(`x6',1)) 
                multiply replace_(`x3',`x6',`x6',`x3',`x4',`x5',`x5',`x4');
                if (count(`x4',1) > count(`x6',1)) 
                multiply replace_(`x4',`x6',`x6',`x4',`x3',`x5',`x5',`x3');
                if (count(`x5',1) > count(`x6',1)) 
                multiply replace_(`x5',`x6',`x6',`x5',`p1',`p2',`p2',`p1');
                if (count(`x3',1) > count(`x4',1)) 
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                
* Now do the reduction via eq. (N15) resp. (6)
                
                if (count(`x5',1)>1);
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^(n5-1) *(-1) * `x6'^n6
                *(-1)/(n5-1)
                *(
                ( n3*`x3'+n4*`x4'+n6*`x6' )*(-1)
                -1/M^2*nom(6-n1-n2-n3-n4-(n5-1)-n6,-3)
                )
                ;
                
                redefine i "0";
                endif;
                
                
                endif;
                endif;
                
                #call ACCU(BNn5)
                #call Conv2exact
        #enddo

#endprocedure

#procedure redBNn34 (p1,p2,x3,x4,x5,x6)
        
* at this stage: n5=1
*                n1,n2>=0
*                n3,n4,n6>0
        
        #do i=1,1
                if ( count(intbn,1) );
                if((count(`x3',1) >= 1)      &&  (count(`x4',1) >= 1) &&
                (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1)  );
                
* sort: n4 >= n3
                
                if (count(`x3',1) > count(`x4',1)) 
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                
* Now do the reduction via eq. (N1) resp. (5)
                
                if((count(`x4',1) > 1)); 
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^(n4-1) *(-1) * `x5'^n5  * `x6'^n6
                *1/2/(n4-1)/M^2
                *(
                +nom(4-2*(n4-1)-n1-n6,-2)
                +n1/`p1'.`p1'*(1/`x5'-1/`x4')
                +n6*`x6'     *(`p2'.`p2'-1/`x4'+2*M^2)
                )
                ;
               
                redefine i "0"; 
                endif;
               
                
                endif;
                endif;
                #call ACCU(BNn34)
                #call Conv2exact
        #enddo
        
#endprocedure

#procedure redBNn12 (p1,p2,x3,x4,x5,x6)
        
* at this stage: n3=n4=n5=1
*                n1,n2 <=>0
*                n6>0
        
        #do i=1,1
                if ( count(intbn,1) );                
                if ( (count(`p1'.`p1',1) > 0)
                && (count(`x3',1) = 1)       &&  (count(`x4',1) = 1)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *M^2*deno( 
                + 2*n1 + 2*n2 + 3*n6 - 7,6
                )        
                *(-1)
                *(
                + p1.p1^-1*p2.p2^-1*x3^-1*x6^-1 * (  - n2*M^-2 )
                
                + p1.p1^-1*p2.p2^-1*x5^-1*x6^-1 * ( n2*M^-2 )
                
                - p1.p1^-1*x3^-1*x6 * ( 2*n6 )
                
                + p1.p1^-1*x3^-1 * (  - n6*M^-2 + M^-2 )
                
                - p1.p1^-1*x4^-1*x5 * ( 2 )
                
                - p1.p1^-1*x4*x5^-1 * ( 2 )
                
                + p1.p1^-1*x6^-1 * 1/M^2*nom( - n2- n6 + 3, -2)
                
                - p1.p1^-1 * nom( 4 - 8*n1 - 4*n6, -8 )
                
                
                ) ;
* N19b
                
                redefine i "0";
                
                endif;
                endif;                

                #call Conv2exact
                .sort

        #enddo
        
        
        #do i=1,1
                
                if ( count(intbn,1) );        
                if ( (count(`p1'.`p1',1)<0 )
                && (count(`x3',1) = 1)       &&  (count(`x4',1) = 1)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *deno(
                + 12  - 8*n1 - 4*n6, -8
                )
                *(-1)
                *(
                + x3^-1*x6 * ( 2*n6 )
                
                - x3^-1 * (  - n6*M^-2 + M^-2 )
                
                + x4^-1*x5 * ( 2 )
                
                + x4*x5^-1 * ( 2 )
                
                - x6^-1 * 1/M^2*nom( - n2- n6 + 3 , -2 )
                
                - p1.p1 * 1/M^2*nom( + 2*n1 + 2*n2 + 3*n6 - 9,6 )
                
                - p2.p2^-1*x3^-1*x6^-1 * (  - n2*M^-2 )
                
                - p2.p2^-1*x5^-1*x6^-1 * ( n2*M^-2 )
                
                
                );
* N18b
                
                redefine i "0";
                
                endif;
                endif;        

                #call Conv2exact
                .sort
        #enddo
        
        #do i=1,1
                
                if ( count(intbn,1) );        
                if ( (count(`p2'.`p2',1)>0 )
                && (count(`x3',1) = 1)       &&  (count(`x4',1) = 1)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1) );  
                
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *M^2*deno(
                + 2*n1 + 2*n2 + 3*n6 - 7, 6
                )
                *(-1)
                *(
                + p1.p1^-1*p2.p2^-1*x4^-1*x6^-1 * (  - n1*M^-2 )
                
                + p1.p1^-1*p2.p2^-1*x5^-1*x6^-1 * ( n1*M^-2 )
                
                - p2.p2^-1*x3^-1*x5 * ( 2 )
                
                - p2.p2^-1*x3*x5^-1 * ( 2 )
                
                - p2.p2^-1*x4^-1*x6 * ( 2*n6 )
                
                + p2.p2^-1*x4^-1 * (  - n6*M^-2 + M^-2 )
                
                + p2.p2^-1*x6^-1 * 1/M^2*nom(  - n1 - n6 + 3, -2 )
                
                - p2.p2^-1 * nom( 4 - 8*n2 - 4*n6 ,-8 )
                
                
                ) ;
* N19a
                
                redefine i "0";
                
                endif;
                endif;        

                #call Conv2exact
                .sort
        #enddo
        
        
        #do i=1,1
                if ( count(intbn,1) );
                if ( (count(`p2'.`p2',1)<0 )
                && (count(`x3',1) = 1)       &&  (count(`x4',1) = 1)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *deno(
                + 12 - 8*n2 - 4*n6, -8
                )
                *(-1)
                *(
                + x3^-1*x5 * ( 2 )
                
                + x3*x5^-1 * ( 2 )
                
                + x4^-1*x6 * ( 2*n6 )
                
                - x4^-1 * (  - n6*M^-2 + M^-2 )
                
                - x6^-1 * 1/M^2*nom(   - n1 - n6 + 3 , -2)
                
                - p1.p1^-1*x4^-1*x6^-1 * (  - n1*M^-2 )
                
                - p1.p1^-1*x5^-1*x6^-1 * ( n1*M^-2 )
                
                - p2.p2 * 1/M^2*nom(  + 2*n1 + 2*n2 + 3*n6 - 9 ,6 )
                
                
                );
* N18a
                
                redefine i "0";
                
                endif;
                endif;        
                #call Conv2exact
                .sort
        #enddo
        
#endprocedure

#procedure redBNn6 (p1,p2,x3,x4,x5,x6)
        
* at this stage: n1=n2=0
*                n3=n4=n5=1
* 		 n6>0    (?)             
        
        
        #do i=1,1

***
*** We do not use tables for BN(1,1,1,n)
***
                
*                 if ( count(intbn,1) );        
*                 #ifdef `TABINT'
*                         if((count(`x3',1) = 1)      &&  (count(`x4',1) = 1) &&
*                         (count(`x5',1) = 1)      &&  (count(`x6',1) > 10) &&
*                         (count(`p1'.`p1',1) = 0) &&  (count(`p2'.`p2',1) = 0)  );
*                 #endif
                
*                 #ifndef `TABINT'
                        
* old version without TabBN:
                        
                        if((count(`x3',1) = 1)      &&  (count(`x4',1) = 1) &&
                        (count(`x5',1) = 1)      &&  (count(`x6',1) > 1) &&
                        (count(`p1'.`p1',1) = 0) &&  (count(`p2'.`p2',1) = 0)  );
*                 #endif
                
*
* do the reduction via eq. (N20)
*
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *1/8/M^2/(n6-1)*deno(4-n6,-2)*deno(n6-3,2)*deno(n6-3,2)
                *(-1)
                *(
                3*M**2*`x3'*`x4'*nom(3 - n6,-2)*nom(-4 + n6,2)/(`x5'*`x6') 
                +
                3*M**2*`x4'*`x5'*nom(3 - n6,-2)*nom(-4 + n6,2)/(`x3'*`x6') 
                +
                M**2*(-1 + n6)*`x3'*nom(-4 + n6,2)*nom(-3 + n6,2)/`x5' 
                +
                M**2*(-1 + n6)*`x5'*nom(-4 + n6,2)*nom(-3 + n6,2)/`x3' 
                -
                (-2 + n6)*`x3'*nom(-4 + n6,2)*nom(-3 + n6,2)/(`x5'*`x6') 
                -
                (-2 + n6)*`x5'*nom(-4 + n6,2)*nom(-3 + n6,2)/(`x3'*`x6') 
                +
                M**2*(-1 + n6)*n6*`x6'*nom(-4 + n6,2)*nom(-3 + n6,2)/`x4' 
                +
                (-2 + n6)*nom(-5 + n6,3)*nom(-4 + n6,2)*nom(-3 + n6,2)/
                (M**2*`x4'*`x6') 
                +       
                nom(-5 + n6,3)*nom(-4 + n6,2)**2*nom(-3 + n6,2)/
                (M^2*`x6'^2) 
                -    
                (-1 + n6)*nom(-4 + n6,2)*nom(-3 + n6,2)*
                nom(-6 + 2*n6,3)/`x4' 
                -
                nom(-4 + n6,2)*nom(-3 + n6,2)*
                nom(-58 + 41*n6 - 7*n6^2,54 - 20*n6,-12)/`x6'
                )
                
                ;
* N20a
                
                redefine i "0";
                
                endif;
                endif;
                .sort
                
                if ( count(intbn,1) );
                if ( (count(`x3',1)<=0)&&(count(`x5',1)<=0) ) discard;
                if ( (count(`x3',1)<=0)&&(count(`x6',1)<=0) ) discard;
                if ( (count(`x4',1)<=0)&&(count(`x5',1)<=0) ) discard;
                if ( (count(`x4',1)<=0)&&(count(`x6',1)<=0) ) discard;
                endif;

                #call ACCU(n6)
                
        #enddo
        
#endprocedure


#procedure symmetryBN
*
* use symmetry for BN
*
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );            
*
* sort: n6 >= n5 >= n4 >= n3
*
        if (  ( count(x3,1) > count(x6,1) ) && ( count(x3,1) >= count(x5,1) )
        && ( count(x3,1) >= count(x4,1) ) ) 
        multiply replace_(x3,x6,x6,x3);
        if (  ( count(x4,1) > count(x6,1) ) && ( count(x4,1) >= count(x5,1) )
        && ( count(x4,1) >= count(x3,1) ) ) 
        multiply replace_(x4,x6,x6,x4);
        if (  ( count(x5,1) > count(x6,1) ) && ( count(x5,1) >= count(x4,1) )
        && ( count(x5,1) >= count(x3,1) ) ) 
        multiply replace_(x5,x6,x6,x5);
        if ( ( count(x3,1) > count(x5,1) ) && ( count(x3,1) >= count(x4,1) ) )
        multiply replace_(x3,x5,x5,x3);
        if ( ( count(x4,1) > count(x5,1) ) && ( count(x4,1) >= count(x3,1) ) )
        multiply replace_(x4,x5,x5,x4);
        if ( count(x3,1) > count(x4,1) ) 
        multiply replace_(x3,x4,x4,x3);
        endif;
        endif;        
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) );
        if ( count(x3,1) > count(x5,1) )  multiply replace_(x3,x5,x5,x3);
        if ( count(x4,1) > count(x6,1) )  multiply replace_(x4,x6,x6,x4);
        if ( count(x5,1) > count(x6,1) )  multiply replace_(x3,x4,x4,x3,x5,x6,x6,x5);
        endif;
        endif;        
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p2.p2,1)==0) );
        if ( count(x3,1) > count(x6,1) )  multiply replace_(x3,x6,x6,x3);
        if ( count(x4,1) > count(x5,1) )  multiply replace_(x4,x5,x5,x4);
        if ( count(x5,1) > count(x6,1) )  multiply replace_(x3,x4,x4,x3,x5,x6,x6,x5);
        endif;
        endif;        
        .sort
        
#endprocedure        


****************************************************************************
* 
* New tables for BN(1,1,1,1) BN(1,1,1,2) BN(1,1,2,2) BN(1,2,2,2) BN(2,2,2,2)
* and it's dalambertian upto dala^11
* 
#procedure BNd0Exact
id,only intbn*dala^0*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^4*miBN*rat(1,1) );


id,only intbn*dala^0*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^2*miBN*rat( - 3*d + 8
         ,8) );


id,only intbn*dala^0*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( miBN*rat(9*d^3 - 81*d^2
          + 242*d - 240,64*d - 256) + Gam(1,1)^3*rat(2,d^4 - 16*d^3 + 96*d^2
          - 256*d + 256) );


id,only intbn*dala^0*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-2*miBN*rat( - 27*d^4
          + 297*d^3 - 1212*d^2 + 2172*d - 1440,512*d - 2048) + Gam(1,1)^3*M^-2
         *rat( - 11*d + 38,4*d^4 - 64*d^3 + 384*d^2 - 1024*d + 1024) );


id,only intbn*dala^0*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-4*miBN*rat(81*d^6 - 
         1512*d^5 + 11115*d^4 - 40224*d^3 + 71700*d^2 - 50400*d,4096*d^2 - 
         40960*d + 98304) + Gam(1,1)^3*M^-4*rat(81*d^3 - 1039*d^2 + 4402*d - 
         6144,32*d^5 - 704*d^4 + 6144*d^3 - 26624*d^2 + 57344*d - 49152) );
#endprocedure

#procedure BNdExact
id,only intbn*dala^1*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^2*miBN*rat(3*d^3 - 11
         *d^2 - 10*d + 48,16*d - 64) + Gam(1,1)^3*M^2*rat(24,d^4 - 16*d^3 + 96
         *d^2 - 256*d + 256) );


id,only intbn*dala^1*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( miBN*rat( - 9*d^4 + 63*
         d^3 - 80*d^2 - 244*d + 480,128*d - 512) + Gam(1,1)^3*rat( - 9*d + 30,
         d^4 - 16*d^3 + 96*d^2 - 256*d + 256) );


id,only intbn*dala^1*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-2*miBN*rat(27*d^6 - 
         432*d^5 + 2481*d^4 - 5424*d^3 - 1284*d^2 + 21792*d - 23040,1024*d^2
          - 10240*d + 24576) + Gam(1,1)^3*M^-2*rat(27*d^3 - 317*d^2 + 1190*d
          - 1440,8*d^5 - 176*d^4 + 1536*d^3 - 6656*d^2 + 14336*d - 12288) );


id,only intbn*dala^1*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-4*miBN*rat( - 81*d^7
          + 1512*d^6 - 11115*d^5 + 40224*d^4 - 71700*d^3 + 50400*d^2,8192*d^2
          - 81920*d + 196608) + Gam(1,1)^3*M^-4*rat( - 81*d^4 + 1039*d^3 - 
         4402*d^2 + 6144*d,64*d^5 - 1408*d^4 + 12288*d^3 - 53248*d^2 + 114688*
         d - 98304) );


id,only intbn*dala^1*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-6*miBN*rat(243*d^9
          - 7047*d^8 + 85239*d^7 - 568701*d^6 + 2425926*d^5 - 7943412*d^4 + 
         23639160*d^3 - 58066848*d^2 + 88646400*d - 58060800,65536*d^3 - 
         1179648*d^2 + 6815744*d - 12582912) + Gam(1,1)^3*M^-6*rat(243*d^6 - 
         4860*d^5 + 28933*d^4 + 19800*d^3 - 842332*d^2 + 3104016*d - 3628800,
         512*d^6 - 15360*d^5 + 188416*d^4 - 1212416*d^3 + 4325376*d^2 - 
         8126464*d + 6291456) );


id,only intbn*dala^2*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( miBN*rat(9*d^6 - 72*d^5
          + 35*d^4 + 272*d^3 + 4148*d^2 - 19872*d + 23040,256*d^2 - 2560*d + 
         6144) + Gam(1,1)^3*rat(9*d^3 + 9*d^2 - 558*d + 1440,2*d^5 - 44*d^4 + 
         384*d^3 - 1664*d^2 + 3584*d - 3072) );


id,only intbn*dala^2*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-2*miBN*rat( - 27*d^6
          + 216*d^5 - 105*d^4 - 816*d^3 - 12444*d^2 + 59616*d - 69120,2048*d
          - 12288) + Gam(1,1)^3*M^-2*rat( - 27*d^3 - 27*d^2 + 1674*d - 4320,16
         *d^4 - 288*d^3 + 1920*d^2 - 5632*d + 6144) );


id,only intbn*dala^2*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-4*miBN*rat(81*d^9 - 
         1917*d^8 + 16893*d^7 - 65775*d^6 + 119874*d^5 - 549180*d^4 + 4551720*
         d^3 - 17205216*d^2 + 29548800*d - 19353600,16384*d^3 - 294912*d^2 + 
         1703936*d - 3145728) + Gam(1,1)^3*M^-4*rat(81*d^6 - 1188*d^5 + 647*
         d^4 + 74408*d^3 - 501364*d^2 + 1296816*d - 1209600,128*d^6 - 3840*d^5
          + 47104*d^4 - 303104*d^3 + 1081344*d^2 - 2031616*d + 1572864) );


id,only intbn*dala^2*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-6*miBN*rat( - 243*
         d^10 + 6561*d^9 - 71145*d^8 + 398223*d^7 - 1288524*d^6 + 3091560*d^5
          - 7752336*d^4 + 10788528*d^3 + 27487296*d^2 - 119232000*d + 
         116121600,131072*d^3 - 2359296*d^2 + 13631488*d - 25165824) + Gam(1,1
         )^3*M^-6*rat( - 243*d^7 + 4374*d^6 - 19213*d^5 - 77666*d^4 + 802732*
         d^3 - 1419352*d^2 - 2579232*d + 7257600,1024*d^6 - 30720*d^5 + 376832
         *d^4 - 2424832*d^3 + 8650752*d^2 - 16252928*d + 12582912) );


id,only intbn*dala^2*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-8*miBN*rat(729*d^12
          - 29160*d^11 + 495558*d^10 - 4747896*d^9 + 29179449*d^8 - 123024960*
         d^7 + 280266840*d^6 + 993951648*d^5 - 14076897840*d^4 + 69070707072*
         d^3 - 181925602560*d^2 + 253232179200*d - 146313216000,1048576*d^4 - 
         29360128*d^3 + 297795584*d^2 - 1291845632*d + 2013265920) + Gam(1,1)^
         3*M^-8*rat(729*d^9 - 22599*d^8 + 279045*d^7 - 2164605*d^6 + 19758954*
         d^5 - 199774188*d^4 + 1359604968*d^3 - 5309789376*d^2 + 10895283072*d
          - 9144576000,8192*d^7 - 327680*d^6 + 5472256*d^5 - 49545216*d^4 + 
         263192576*d^3 - 822083584*d^2 + 1400897536*d - 1006632960) );


id,only intbn*dala^3*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-2*miBN*rat(27*d^9 - 
         351*d^8 + 591*d^7 + 3963*d^6 + 17286*d^5 + 55596*d^4 - 3119880*d^3 + 
         17590368*d^2 - 38534400*d + 30412800,4096*d^3 - 73728*d^2 + 425984*d
          - 786432) + Gam(1,1)^3*M^-2*rat(27*d^6 - 108*d^5 - 867*d^4 - 24312*
         d^3 + 350052*d^2 - 1431792*d + 1900800,32*d^6 - 960*d^5 + 11776*d^4
          - 75776*d^3 + 270336*d^2 - 507904*d + 393216) );


id,only intbn*dala^3*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-4*miBN*rat( - 81*
         d^10 + 1431*d^9 - 6687*d^8 - 3615*d^7 + 3624*d^6 + 75216*d^5 + 
         10137984*d^4 - 96449424*d^3 + 361868352*d^2 - 630720000*d + 425779200
         ,32768*d^3 - 589824*d^2 + 3407872*d - 6291456) + Gam(1,1)^3*M^-4*rat(
          - 81*d^7 + 702*d^6 + 1089*d^5 + 60798*d^4 - 1390524*d^3 + 9196104*
         d^2 - 25747488*d + 26611200,256*d^6 - 7680*d^5 + 94208*d^4 - 606208*
         d^3 + 2162688*d^2 - 4063232*d + 3145728) );


id,only intbn*dala^3*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-6*miBN*rat(243*d^12
          - 7776*d^11 + 93258*d^10 - 488592*d^9 + 849099*d^8 + 1157712*d^7 - 
         34392120*d^6 + 640660704*d^5 - 5398794384*d^4 + 23666752896*d^3 - 
         57489027840*d^2 + 73943193600*d - 39481344000,262144*d^4 - 7340032*
         d^3 + 74448896*d^2 - 322961408*d + 503316480) + Gam(1,1)^3*M^-6*rat(
         243*d^9 - 5589*d^8 + 38583*d^7 - 217911*d^6 + 5670606*d^5 - 79226532*
         d^4 + 528775032*d^3 - 1862844096*d^2 + 3367361664*d - 2467584000,2048
         *d^7 - 81920*d^6 + 1368064*d^5 - 12386304*d^4 + 65798144*d^3 - 
         205520896*d^2 + 350224384*d - 251658240) );


id,only intbn*dala^3*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-8*miBN*rat( - 729*
         d^13 + 26244*d^12 - 378918*d^11 + 2765664*d^10 - 10187865*d^9 + 
         6307164*d^8 + 211833000*d^7 - 2115019008*d^6 + 10101091248*d^5 - 
         12763115712*d^4 - 94357225728*d^3 + 474470231040*d^2 - 866615500800*d
          + 585252864000,2097152*d^4 - 58720256*d^3 + 595591168*d^2 - 
         2583691264*d + 4026531840) + Gam(1,1)^3*M^-8*rat( - 729*d^10 + 19683*
         d^9 - 188649*d^8 + 1048425*d^7 - 11100534*d^6 + 120738372*d^5 - 
         560508216*d^4 - 128630496*d^3 + 10343874432*d^2 - 34436556288*d + 
         36578304000,16384*d^7 - 655360*d^6 + 10944512*d^5 - 99090432*d^4 + 
         526385152*d^3 - 1644167168*d^2 + 2801795072*d - 2013265920) );


id,only intbn*dala^3*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-10*miBN*rat(2187*
         d^15 - 112995*d^14 + 2488320*d^13 - 30690090*d^12 + 233341911*d^11 - 
         1056919419*d^10 + 914368158*d^9 + 25853022648*d^8 - 99950514912*d^7
          - 1704085546032*d^6 + 24215500975968*d^5 - 149873783099904*d^4 + 
         537525560885760*d^3 - 1152067285094400*d^2 + 1373435191296000*d - 
         702303436800000,16777216*d^5 - 671088640*d^4 + 10401873920*d^3 - 
         77846282240*d^2 + 280246616064*d - 386547056640) + Gam(1,1)^3*M^-10*
         rat(2187*d^12 - 93312*d^11 + 1609146*d^10 - 14353200*d^9 + 44877171*
         d^8 + 1142803824*d^7 - 29378517864*d^6 + 368205556416*d^5 - 
         2794923757392*d^4 + 13267456636416*d^3 - 38540796509952*d^2 + 
         62747225026560*d - 43893964800000,131072*d^8 - 6815744*d^7 + 
         150470656*d^6 - 1843396608*d^5 + 13723762688*d^4 - 63686311936*d^3 + 
         180254408704*d^2 - 285078454272*d + 193273528320) );


id,only intbn*dala^4*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-4*miBN*rat(81*d^12
          - 1512*d^11 + 4230*d^10 + 48120*d^9 - 170751*d^8 + 1875792*d^7 - 
         29960616*d^6 - 111171552*d^5 + 3910322256*d^4 - 26660373888*d^3 + 
         85368119040*d^2 - 135415756800*d + 85929984000,65536*d^4 - 1835008*
         d^3 + 18612224*d^2 - 80740352*d + 125829120) + Gam(1,1)^3*M^-4*rat(81
         *d^9 - 783*d^8 - 4275*d^7 + 30219*d^6 - 1371798*d^5 + 43059348*d^4 - 
         464207640*d^3 + 2334165696*d^2 - 5662942848*d + 5370624000,512*d^7 - 
         20480*d^6 + 342016*d^5 - 3096576*d^4 + 16449536*d^3 - 51380224*d^2 + 
         87556096*d - 62914560) );


id,only intbn*dala^4*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-6*miBN*rat( - 243*
         d^13 + 5832*d^12 - 36882*d^11 - 76680*d^10 + 1282173*d^9 - 8359392*
         d^8 + 119894520*d^7 - 145855200*d^6 - 13509711600*d^5 + 142546277760*
         d^4 - 682670339328*d^3 + 1772137175040*d^2 - 2424442060800*d + 
         1374879744000,524288*d^4 - 14680064*d^3 + 148897792*d^2 - 645922816*d
          + 1006632960) + Gam(1,1)^3*M^-6*rat( - 243*d^10 + 3645*d^9 + 297*d^8
          - 159057*d^7 + 4598898*d^6 - 151126812*d^5 + 2081572488*d^4 - 
         14429819328*d^3 + 54335479680*d^2 - 106718957568*d + 85929984000,4096
         *d^7 - 163840*d^6 + 2736128*d^5 - 24772608*d^4 + 131596288*d^3 - 
         411041792*d^2 + 700448768*d - 503316480) );


id,only intbn*dala^4*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-8*miBN*rat(729*d^15
          - 29889*d^14 + 456192*d^13 - 2829006*d^12 - 221283*d^11 + 110369079*
         d^10 - 1066533750*d^9 + 7165439208*d^8 + 16357988448*d^7 - 
         946495921680*d^6 + 9500913239328*d^5 - 50585129436672*d^4 + 
         162036446270976*d^3 - 314046340116480*d^2 + 340642249113600*d - 
         159188779008000,4194304*d^5 - 167772160*d^4 + 2600468480*d^3 - 
         19461570560*d^2 + 70061654016*d - 96636764160) + Gam(1,1)^3*M^-8*rat(
         729*d^12 - 23328*d^11 + 233118*d^10 - 252720*d^9 - 20371215*d^8 + 
         633538704*d^7 - 12501583608*d^6 + 144168451392*d^5 - 1002014245488*
         d^4 + 4295686181376*d^3 - 11155592975616*d^2 + 16117693894656*d - 
         9949298688000,32768*d^8 - 1703936*d^7 + 37617664*d^6 - 460849152*d^5
          + 3430940672*d^4 - 15921577984*d^3 + 45063602176*d^2 - 71269613568*d
          + 48318382080) );


id,only intbn*dala^4*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-10*miBN*rat( - 2187*
         d^16 + 99873*d^15 - 1810350*d^14 + 15760170*d^13 - 49201371*d^12 - 
         343132047*d^11 + 5427148356*d^10 - 31339231596*d^9 - 55167620976*d^8
          + 2303788635504*d^7 - 13990987699776*d^6 + 4580777244096*d^5 + 
         361717137713664*d^4 - 2073086080220160*d^3 + 5538968519270400*d^2 - 
         7538307710976000*d + 4213820620800000,33554432*d^5 - 1342177280*d^4
          + 20803747840*d^3 - 155692564480*d^2 + 560493232128*d - 773094113280
         ) + Gam(1,1)^3*M^-10*rat( - 2187*d^13 + 80190*d^12 - 1049274*d^11 + 
         4698324*d^10 + 41242029*d^9 - 1412066850*d^8 + 22521694920*d^7 - 
         191934449232*d^6 + 585690418896*d^5 + 3502085907936*d^4 - 
         41063943308544*d^3 + 168497554033152*d^2 - 332589385359360*d + 
         263363788800000,262144*d^8 - 13631488*d^7 + 300941312*d^6 - 
         3686793216*d^5 + 27447525376*d^4 - 127372623872*d^3 + 360508817408*
         d^2 - 570156908544*d + 386547056640) );


id,only intbn*dala^4*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-12*miBN*rat(6561*
         d^17 - 367416*d^16 + 8474625*d^15 - 102789000*d^14 + 651752487*d^13
          - 557195112*d^12 - 31704255045*d^11 + 259910425800*d^10 + 
         541449337068*d^9 - 14863395882528*d^8 - 121739988976560*d^7 + 
         3815602724822400*d^6 - 36323239444583616*d^5 + 193194912017005056*d^4
          - 633239229292139520*d^3 + 1274374072606924800*d^2 - 
         1448842040819712000*d + 713821213163520000,268435456*d^5 - 
         12348030976*d^4 + 217969590272*d^3 - 1831803551744*d^2 + 
         7267084664832*d - 10823317585920) + Gam(1,1)^3*M^-12*rat(6561*d^14 - 
         308367*d^13 + 5581224*d^12 - 46482498*d^11 + 115207029*d^10 - 
         433073127*d^9 + 102149117058*d^8 - 3390454992216*d^7 + 57135461050272
         *d^6 - 585766082751600*d^5 + 3859367113195680*d^4 - 16463017297873152
         *d^3 + 44029209591922176*d^2 - 67170496982999040*d + 
         44613825822720000,2097152*d^8 - 121634816*d^7 + 2961178624*d^6 - 
         39510343680*d^5 + 316418293760*d^4 - 1561757483008*d^3 + 
         4655744548864*d^2 - 7692286427136*d + 5411658792960) );


id,only intbn*dala^5*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-6*miBN*rat(243*d^14
          - 4617*d^13 - 3942*d^12 + 401058*d^11 - 1129461*d^10 + 11036079*d^9
          - 41980776*d^8 - 3763887096*d^7 + 21152057808*d^6 + 383256272688*d^5
          - 5404645775232*d^4 + 29550824671488*d^3 - 84100568279040*d^2 + 
         123647982796800*d - 74392141824000,1048576*d^4 - 35651584*d^3 + 
         436207616*d^2 - 2248146944*d + 4026531840) + Gam(1,1)^3*M^-6*rat(243*
         d^11 - 2430*d^10 - 30186*d^9 + 192564*d^8 + 1209171*d^7 - 107062158*
         d^6 + 4567191012*d^5 - 73486290216*d^4 + 579872609472*d^3 - 
         2443288959360*d^2 + 5295651084288*d - 4649508864000,8192*d^7 - 376832
         *d^6 + 7143424*d^5 - 72351744*d^4 + 423624704*d^3 - 1438646272*d^2 + 
         2634022912*d - 2013265920) );


id,only intbn*dala^5*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-8*miBN*rat( - 729*
         d^15 + 18225*d^14 - 71280*d^13 - 1274130*d^12 + 10607427*d^11 - 
         53438535*d^10 + 324591750*d^9 + 10536007320*d^8 - 131206141152*d^7 - 
         769031777520*d^6 + 23112550234080*d^5 - 185936097968640*d^4 + 
         784216548923904*d^3 - 1884754177413120*d^2 + 2448840115814400*d - 
         1339058552832000,8388608*d^4 - 285212672*d^3 + 3489660928*d^2 - 
         17985175552*d + 32212254720) + Gam(1,1)^3*M^-8*rat( - 729*d^12 + 
         11664*d^11 + 46818*d^10 - 1121040*d^9 - 161361*d^8 + 342951552*d^7 - 
         15628691880*d^6 + 302668308864*d^5 - 3062371052304*d^4 + 
         17767573848576*d^3 - 59866154521344*d^2 + 109270246109184*d - 
         83691159552000,65536*d^7 - 3014656*d^6 + 57147392*d^5 - 578813952*d^4
          + 3388997632*d^3 - 11509170176*d^2 + 21072183296*d - 16106127360) );


id,only intbn*dala^5*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-10*miBN*rat(2187*
         d^18 - 110808*d^17 + 2064771*d^16 - 14081040*d^15 - 56712771*d^14 + 
         1674178344*d^13 - 13691393343*d^12 + 34768195200*d^11 + 918312757956*
         d^10 - 11512751811264*d^9 - 41959411694352*d^8 + 2213094593502720*d^7
          - 24955415148267072*d^6 + 157292459107273728*d^5 - 
         631104495033240576*d^4 + 1652600320416645120*d^3 - 
         2746826647481548800*d^2 + 2638184132247552000*d - 1116943385886720000
         ,67108864*d^6 - 3623878656*d^5 + 79188459520*d^4 - 893890068480*d^3
          + 5480378269696*d^2 - 17239998726144*d + 21646635171840) + Gam(1,1)^
         3*M^-10*rat(2187*d^15 - 91125*d^14 + 1205280*d^13 - 1418310*d^12 - 
         96153129*d^11 - 124449693*d^10 + 61730637438*d^9 - 1966424387160*d^8
          + 34849542704160*d^7 - 391253559826512*d^6 + 2911132991157600*d^5 - 
         14578076120884224*d^4 + 48665879186932224*d^3 - 103913989348478976*
         d^2 + 128397731444490240*d - 69808961617920000,524288*d^9 - 34603008*
         d^8 + 983564288*d^7 - 15799943168*d^6 + 158125260800*d^5 - 
         1023275958272*d^4 + 4287451103232*d^3 - 11234560704512*d^2 + 
         16737487552512*d - 10823317585920) );


id,only intbn*dala^5*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-12*miBN*rat( - 6561*
         d^18 + 314928*d^17 - 5535297*d^16 + 34992000*d^15 + 170559513*d^14 - 
         4656824784*d^13 + 36161815941*d^12 - 6276385440*d^11 - 2620732743468*
         d^10 + 10531801185984*d^9 + 240647156036784*d^8 - 2841682813009920*
         d^7 + 5798417646004416*d^6 + 97391003539663872*d^5 - 
         912320066843900928*d^4 + 3791539761730191360*d^3 - 
         8746150540035686400*d^2 + 10876915113394176000*d - 
         5710569705308160000,536870912*d^5 - 24696061952*d^4 + 435939180544*
         d^3 - 3663607103488*d^2 + 14534169329664*d - 21646635171840) + Gam(1,
         1)^3*M^-12*rat( - 6561*d^15 + 255879*d^14 - 3114288*d^13 + 1832706*
         d^12 + 256652955*d^11 - 488583105*d^10 - 98684532042*d^9 + 
         2573262055752*d^8 - 30011821112544*d^7 + 128682394349424*d^6 + 
         826761548817120*d^5 - 14411919607692288*d^4 + 87674928791063040*d^3
          - 285063179752378368*d^2 + 492750150041272320*d - 356910606581760000
         ,4194304*d^8 - 243269632*d^7 + 5922357248*d^6 - 79020687360*d^5 + 
         632836587520*d^4 - 3123514966016*d^3 + 9311489097728*d^2 - 
         15384572854272*d + 10823317585920) );


id,only intbn*dala^5*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-14*miBN*rat(19683*
         d^21 - 1515591*d^20 + 49240305*d^19 - 858487167*d^18 + 7928588691*
         d^17 - 14881071897*d^16 - 586846673325*d^15 + 7490595070971*d^14 - 
         23682221677746*d^13 - 185132725619148*d^12 - 1483335193481640*d^11 + 
         52574629671692784*d^10 + 231635411393831712*d^9 - 
         17640357491038895424*d^8 + 258775405412924338560*d^7 - 
         2144043203767542425088*d^6 + 11628418305466626140160*d^5 - 
         42969004169148102205440*d^4 + 107761276964693694873600*d^3 - 
         176009889416894152704000*d^2 + 169135359688059125760000*d - 
         72576149527461888000000,4294967296*d^7 - 300647710720*d^6 + 
         8778913153024*d^5 - 138297946931200*d^4 + 1266087639384064*d^3 - 
         6715267266641920*d^2 + 19039143346569216*d - 22166154415964160) + 
         Gam(1,1)^3*M^-14*rat(19683*d^18 - 1338444*d^17 + 36840015*d^16 - 
         501260400*d^15 + 2667834549*d^14 + 20114862516*d^13 - 667149898707*
         d^12 + 20713940389944*d^11 - 746779038794820*d^10 + 19455927510604752
         *d^9 - 343391872235997456*d^8 + 4216875606272830464*d^7 - 
         36849278535731948736*d^6 + 231145682015079744768*d^5 - 
         1034621296481450769408*d^4 + 3227324876453857996800*d^3 - 
         6666788615102266245120*d^2 + 8195267231191439769600*d - 
         4536009345466368000000,33554432*d^10 - 2751463424*d^9 + 98381594624*
         d^8 - 2018366193664*d^7 + 26299158495232*d^6 - 227409928388608*d^5 + 
         1322231451877376*d^4 - 5109361814798336*d^3 + 12575389364781056*d^2
          - 17831879579271168*d + 11083077207982080) );


id,only intbn*dala^6*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-8*miBN*rat(729*d^18
          - 23328*d^17 + 118665*d^16 + 3024864*d^15 - 36989937*d^14 + 
         115062624*d^13 + 1048206195*d^12 - 45453455712*d^11 + 266112165132*
         d^10 + 7520947004736*d^9 - 81470677926960*d^8 - 821202070560768*d^7
          + 22307185784437056*d^6 - 206065496349637632*d^5 + 
         1078845337746201600*d^4 - 3496440014440464384*d^3 + 
         6972043537232363520*d^2 - 7864254361003622400*d + 3847948622364672000
         ,16777216*d^6 - 905969664*d^5 + 19797114880*d^4 - 223472517120*d^3 + 
         1370094567424*d^2 - 4309999681536*d + 5411658792960) + Gam(1,1)^3*
         M^-8*rat(729*d^15 - 16767*d^14 - 45360*d^13 + 2976750*d^12 - 9954243*
         d^11 - 33555447*d^10 - 10947036006*d^9 + 693209196792*d^8 - 
         20098955584800*d^7 + 329179413474576*d^6 - 3330359927329248*d^5 + 
         21695563148608512*d^4 - 91525637479684608*d^3 + 242140818806083584*
         d^2 - 365608792195006464*d + 240496788897792000,131072*d^9 - 8650752*
         d^8 + 245891072*d^7 - 3949985792*d^6 + 39531315200*d^5 - 255818989568
         *d^4 + 1071862775808*d^3 - 2808640176128*d^2 + 4184371888128*d - 
         2705829396480) );


id,only intbn*dala^6*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-10*miBN*rat( - 2187*
         d^19 + 84564*d^18 - 822555*d^17 - 6701292*d^16 + 171467091*d^15 - 
         1084986612*d^14 - 843366105*d^13 + 157324491036*d^12 - 1707405609636*
         d^11 - 17240597711568*d^10 + 394830973875600*d^9 + 834192653143104*
         d^8 - 83345598764526528*d^7 + 1064340204737654016*d^6 - 
         7357845940231357440*d^5 + 32066226798245425152*d^4 - 
         90844930900506378240*d^3 + 163033633827658137600*d^2 - 
         168828933087166464000*d + 76958972447293440000,134217728*d^6 - 
         7247757312*d^5 + 158376919040*d^4 - 1787780136960*d^3 + 
         10960756539392*d^2 - 34479997452288*d + 43293270343680) + Gam(1,1)^3*
         M^-10*rat( - 2187*d^16 + 64881*d^15 - 199260*d^14 - 9837450*d^13 + 
         89397729*d^12 - 98418519*d^11 + 32169999078*d^10 - 2298568310496*d^9
          + 74161050690240*d^8 - 1389517352119728*d^7 + 16574668051479264*d^6
          - 131693887992410496*d^5 + 708488175411224064*d^4 - 
         2556935206011942912*d^3 + 5939642752706691072*d^2 - 
         8033666210593505280*d + 4809935777955840000,1048576*d^9 - 69206016*
         d^8 + 1967128576*d^7 - 31599886336*d^6 + 316250521600*d^5 - 
         2046551916544*d^4 + 8574902206464*d^3 - 22469121409024*d^2 + 
         33474975105024*d - 21646635171840) );


id,only intbn*dala^6*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-12*miBN*rat(6561*
         d^21 - 400221*d^20 + 8855163*d^19 - 63228357*d^18 - 687535695*d^17 + 
         17083906749*d^16 - 127275124239*d^15 - 219250631079*d^14 + 
         15629662753434*d^13 - 96248715700356*d^12 - 1656461219580216*d^11 + 
         23085994374682704*d^10 + 193518107895332448*d^9 - 7556945421774845376
         *d^8 + 92746571191450169472*d^7 - 674560960832797843968*d^6 + 
         3265687118050569418752*d^5 - 10858680215674760503296*d^4 + 
         24622799773686455992320*d^3 - 36489742314143337676800*d^2 + 
         31909672310283436032000*d - 12496803086016184320000,1073741824*d^7 - 
         75161927680*d^6 + 2194728288256*d^5 - 34574486732800*d^4 + 
         316521909846016*d^3 - 1678816816660480*d^2 + 4759785836642304*d - 
         5541538603991040) + Gam(1,1)^3*M^-12*rat(6561*d^18 - 341172*d^17 + 
         5666517*d^16 - 5563728*d^15 - 859971897*d^14 + 9680327820*d^13 - 
         119764617057*d^12 + 7770358145160*d^11 - 326992781262060*d^10 + 
         8155736859033072*d^9 - 133317455718986928*d^8 + 1503275266378538496*
         d^7 - 11998567240402752576*d^6 + 68418377029704939264*d^5 - 
         277129970605456856064*d^4 + 778868441228281909248*d^3 - 
         1443526852246727786496*d^2 + 1585549763518558371840*d - 
         781050192876011520000,8388608*d^10 - 687865856*d^9 + 24595398656*d^8
          - 504591548416*d^7 + 6574789623808*d^6 - 56852482097152*d^5 + 
         330557862969344*d^4 - 1277340453699584*d^3 + 3143847341195264*d^2 - 
         4457969894817792*d + 2770769301995520) );


id,only intbn*dala^6*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-14*miBN*rat( - 19683
         *d^22 + 1318761*d^21 - 34084395*d^20 + 366084117*d^19 + 656282979*
         d^18 - 64404815013*d^17 + 735657392295*d^16 - 1622128337721*d^15 - 
         51223729031964*d^14 + 421954942396608*d^13 + 3334662449673120*d^12 - 
         37741277736876384*d^11 - 757381708110759552*d^10 + 
         15324003377100578304*d^9 - 82371830502535384320*d^8 - 
         443710850361700960512*d^7 + 9812013732208798110720*d^6 - 
         73315178885518159196160*d^5 + 321928764726787327180800*d^4 - 
         901602880230042796032000*d^3 + 1590963534480882401280000*d^2 - 
         1618777447353129369600000*d + 725761495274618880000000,8589934592*d^7
          - 601295421440*d^6 + 17557826306048*d^5 - 276595893862400*d^4 + 
         2532175278768128*d^3 - 13430534533283840*d^2 + 38078286693138432*d - 
         44332308831928320) + Gam(1,1)^3*M^-14*rat( - 19683*d^19 + 1141614*
         d^18 - 23455575*d^17 + 132860250*d^16 + 2344769451*d^15 - 46793208006
         *d^14 + 466001273547*d^13 - 14042441402874*d^12 + 539639634895380*
         d^11 - 11988137122656552*d^10 + 148832597129949936*d^9 - 
         782956883912855904*d^8 - 5319477526996355904*d^7 + 
         137347103342239742592*d^6 - 1276835523669346678272*d^5 + 
         7118888088360649697280*d^4 - 25606460149436313722880*d^3 + 
         58472618919831222681600*d^2 - 77416662966448029696000*d + 
         45360093454663680000000,67108864*d^10 - 5502926848*d^9 + 196763189248
         *d^8 - 4036732387328*d^7 + 52598316990464*d^6 - 454819856777216*d^5
          + 2644462903754752*d^4 - 10218723629596672*d^3 + 25150778729562112*
         d^2 - 35663759158542336*d + 22166154415964160) );


id,only intbn*dala^6*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-16*miBN*rat(59049*
         d^24 - 5353776*d^23 + 203207292*d^22 - 4027719168*d^21 + 37698320646*
         d^20 + 90720235872*d^19 - 7138582730676*d^18 + 82532293364160*d^17 - 
         191337103781127*d^16 - 3828076488223920*d^15 + 8771160561929136*d^14
          + 82207087637266368*d^13 + 12274439542118834400*d^12 - 
         199852526305882927872*d^11 - 1522572833056963283712*d^10 + 
         85107163722451798496256*d^9 - 1332853181249024375394048*d^8 + 
         12459054672037995795357696*d^7 - 79164737600394055794647040*d^6 + 
         356163902498462783461392384*d^5 - 1141160131189932655777873920*d^4 + 
         2554346154770715394179072000*d^3 - 3801423216058473798696960000*d^2
          + 3380750183334554432962560000*d - 1359150310127150353612800000,
         68719476736*d^8 - 6047313952768*d^7 + 227049151135744*d^6 - 
         4741094138970112*d^5 + 60087210946330624*d^4 - 472077516408881152*d^3
          + 2238623266337980416*d^2 - 5837931754467360768*d + 
         6383852471797678080) + Gam(1,1)^3*M^-16*rat(59049*d^21 - 4822335*d^20
          + 158743395*d^19 - 2507502663*d^18 + 11949984513*d^17 + 252299026071
         *d^16 - 5190078599559*d^15 + 9724415782131*d^14 + 2756966734812690*
         d^13 - 143559324370962324*d^12 + 4734888147824330376*d^11 - 
         109894659267031815600*d^10 + 1848053366389595032416*d^9 - 
         22955818180294982731968*d^8 + 213143383992533046148224*d^7 - 
         1485231500101659170761728*d^6 + 7731200928826175442112512*d^5 - 
         29609953771178275723689984*d^4 + 80958307183032759520886784*d^3 - 
         149455951843095594801561600*d^2 + 166812799903220196546969600*d - 
         84946894382946897100800000,536870912*d^11 - 53687091200*d^10 + 
         2366526980096*d^9 - 60627758350336*d^8 + 1002075999698944*d^7 - 
         11212716500844544*d^6 + 86649762605957120*d^5 - 462552447177457664*
         d^4 + 1672702432498417664*d^3 - 3907022210325282816*d^2 + 
         5312910554157809664*d - 3191926235898839040) );


id,only intbn*dala^7*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-10*miBN*rat(2187*
         d^21 - 86751*d^20 + 552825*d^19 + 18633321*d^18 - 298195317*d^17 + 
         558520767*d^16 + 25955109675*d^15 - 480387320973*d^14 + 793287423054*
         d^13 + 96048006807060*d^12 - 329547587722920*d^11 - 22895838989653392
         *d^10 + 201796284845855520*d^9 + 3473995061902857408*d^8 - 
         89505224671765714560*d^7 + 930536035791618783744*d^6 - 
         5819181975847951570944*d^5 + 23837507875276989014016*d^4 - 
         64742319306552938987520*d^3 + 112726977624043801804800*d^2 - 
         114153708475024146432000*d + 51166327001233489920000,268435456*d^7 - 
         18790481920*d^6 + 548682072064*d^5 - 8643621683200*d^4 + 
         79130477461504*d^3 - 419704204165120*d^2 + 1189946459160576*d - 
         1385384650997760) + Gam(1,1)^3*M^-10*rat(2187*d^18 - 67068*d^17 - 
         90153*d^16 + 19204128*d^15 - 126791379*d^14 - 958652892*d^13 + 
         20736029061*d^12 - 1624849997976*d^11 + 109083562772796*d^10 - 
         4129089071598000*d^9 + 95981978073975408*d^8 - 1459382904741447936*
         d^7 + 15104998875867746112*d^6 - 108609344608482429696*d^5 + 
         543692649481401664512*d^4 - 1861181913866998763520*d^3 + 
         4157100622720272531456*d^2 - 5460121339639210967040*d + 
         3197895437577093120000,2097152*d^10 - 171966464*d^9 + 6148849664*d^8
          - 126147887104*d^7 + 1643697405952*d^6 - 14213120524288*d^5 + 
         82639465742336*d^4 - 319335113424896*d^3 + 785961835298816*d^2 - 
         1114492473704448*d + 692692325498880) );


id,only intbn*dala^7*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-12*miBN*rat( - 6561*
         d^22 + 308367*d^21 - 3566997*d^20 - 43737813*d^19 + 1304519013*d^18
          - 8235859275*d^17 - 65577872151*d^16 + 2012174375769*d^15 - 
         12948383330568*d^14 - 270691697113992*d^13 + 3101698912924080*d^12 + 
         61437470039055936*d^11 - 1109097312309941184*d^10 - 
         5982466919099750784*d^9 + 344943565377160006656*d^8 - 
         4760723050153702071552*d^7 + 37929338714959467955200*d^6 - 
         199534527094485901602816*d^5 + 718652131175752575270912*d^4 - 
         1762511957616296063139840*d^3 + 2822454633154036079001600*d^2 - 
         2664880567454231691264000*d + 1125659194027136778240000,2147483648*
         d^7 - 150323855360*d^6 + 4389456576512*d^5 - 69148973465600*d^4 + 
         633043819692032*d^3 - 3357633633320960*d^2 + 9519571673284608*d - 
         11083077207982080) + Gam(1,1)^3*M^-12*rat( - 6561*d^19 + 249318*d^18
          - 1205037*d^17 - 59595750*d^16 + 802864953*d^15 + 86548338*d^14 - 
         83298450807*d^13 + 5330742633270*d^12 - 362997388273860*d^11 + 
         14787105595795512*d^10 - 378785893797082224*d^9 + 6489752231851802784
         *d^8 - 77421420531915092928*d^7 + 658138009094537703552*d^6 - 
         4020483529830818446848*d^5 + 17544784030191832909824*d^4 - 
         53417303973234790391808*d^3 + 107836577718763628593152*d^2 - 
         129716355784793920634880*d + 70353699626696048640000,16777216*d^10 - 
         1375731712*d^9 + 49190797312*d^8 - 1009183096832*d^7 + 13149579247616
         *d^6 - 113704964194304*d^5 + 661115725938688*d^4 - 2554680907399168*
         d^3 + 6287694682390528*d^2 - 8915939789635584*d + 5541538603991040) )
         ;


id,only intbn*dala^7*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-14*miBN*rat(19683*
         d^24 - 1417176*d^23 + 36505404*d^22 - 263227320*d^21 - 5719820022*
         d^20 + 140993726328*d^19 - 956793282372*d^18 - 7861524712488*d^17 + 
         213681578188203*d^16 - 864884342285568*d^15 - 22163945025500208*d^14
          + 106932190797117504*d^13 + 5916430281551472288*d^12 - 
         65160786203483929344*d^11 - 1048052594650080432384*d^10 + 
         35055527211570387649536*d^9 - 463678726258441596333312*d^8 + 
         3820773788623236177653760*d^7 - 21722859246912856794685440*d^6 + 
         88078063785720129556512768*d^5 - 255388725017482877800611840*d^4 + 
         518812164522634280435712000*d^3 - 702360145416656364503040000*d^2 + 
         569459957555740457041920000*d - 209194240963444840857600000,
         17179869184*d^8 - 1511828488192*d^7 + 56762287783936*d^6 - 
         1185273534742528*d^5 + 15021802736582656*d^4 - 118019379102220288*d^3
          + 559655816584495104*d^2 - 1459482938616840192*d + 
         1595963117949419520) + Gam(1,1)^3*M^-14*rat(19683*d^21 - 1240029*d^20
          + 24990849*d^19 - 14414517*d^18 - 6377803029*d^17 + 84971689605*d^16
          - 68710448205*d^15 - 21179736402183*d^14 + 1337690912369670*d^13 - 
         62644641619733532*d^12 + 1983393526222546776*d^11 - 
         43437775641982272144*d^10 + 680640736598568166176*d^9 - 
         7815715992589115184192*d^8 + 66696637619385097937280*d^7 - 
         425094276869066032081920*d^6 + 2015164329339759202455552*d^5 - 
         7000050936580769086783488*d^4 + 17290842897630135606018048*d^3 - 
         28726739615261544225177600*d^2 + 28745546133193138844467200*d - 
         13074640060215302553600000,134217728*d^11 - 13421772800*d^10 + 
         591631745024*d^9 - 15156939587584*d^8 + 250518999924736*d^7 - 
         2803179125211136*d^6 + 21662440651489280*d^5 - 115638111794364416*d^4
          + 418175608124604416*d^3 - 976755552581320704*d^2 + 
         1328227638539452416*d - 797981558974709760) );


id,only intbn*dala^7*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-16*miBN*rat( - 59049
         *d^25 + 4645188*d^24 - 138961980*d^23 + 1589231664*d^22 + 10634309370
         *d^21 - 543100083624*d^20 + 6049939900212*d^19 + 3130699403952*d^18
          - 799050416588793*d^17 + 6124121733597444*d^16 + 37165757296757904*
         d^15 - 187461014380416000*d^14 - 13260924593766030816*d^13 + 
         52559251800456915072*d^12 + 3920803148727558418176*d^11 - 
         66836289725768239091712*d^10 + 311567216579602793438976*d^9 + 
         3535183502950296709370880*d^8 - 70343918464061893749645312*d^7 + 
         593812948706265886074372096*d^6 - 3132806698791620745758834688*d^5 + 
         11139575419508476475155415040*d^4 - 26850730641190110931451904000*d^3
          + 42236328409367131151400960000*d^2 - 39209851889887502841937920000*
         d + 16309803721525804243353600000,137438953472*d^8 - 12094627905536*
         d^7 + 454098302271488*d^6 - 9482188277940224*d^5 + 120174421892661248
         *d^4 - 944155032817762304*d^3 + 4477246532675960832*d^2 - 
         11675863508934721536*d + 12767704943595356160) + Gam(1,1)^3*M^-16*
         rat( - 59049*d^22 + 4113747*d^21 - 100875375*d^20 + 602581923*d^19 + 
         18140047443*d^18 - 395698840227*d^17 + 2162490286707*d^16 + 
         52556527412577*d^15 - 2873659724198262*d^14 + 110475723553210044*d^13
          - 3012176255372782488*d^12 + 53076001493139851088*d^11 - 
         529317455185213245216*d^10 + 779177783619842342976*d^9 + 
         62326434171006746635392*d^8 - 1072489107808737383016960*d^7 + 
         10091577072393734607028224*d^6 - 63164457374735829581660160*d^5 + 
         274361138071106549163393024*d^4 - 822043734353297519449079808*d^3 + 
         1626658622213926941071769600*d^2 - 1916806704455695461462835200*d + 
         1019362732595362765209600000,1073741824*d^11 - 107374182400*d^10 + 
         4733053960192*d^9 - 121255516700672*d^8 + 2004151999397888*d^7 - 
         22425433001689088*d^6 + 173299525211914240*d^5 - 925104894354915328*
         d^4 + 3345404864996835328*d^3 - 7814044420650565632*d^2 + 
         10625821108315619328*d - 6383852471797678080) );


id,only intbn*dala^7*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-18*miBN*rat(177147*
         d^27 - 18600435*d^26 + 810900234*d^25 - 17939794788*d^24 + 
         161306292186*d^23 + 1682833711230*d^22 - 66586642336152*d^21 + 
         729494669607228*d^20 + 36273515353875*d^19 - 78032367182890179*d^18
          + 569765341109379942*d^17 - 2856145565032373376*d^16 + 
         121250177978574518784*d^15 - 14607146198728627488*d^14 - 
         67326073443568154453184*d^13 + 914015658569612157400320*d^12 + 
         10996818542821241718121728*d^11 - 527861802276298196496421632*d^10 + 
         8822999527802188753569338880*d^9 - 91893886220254586959832592384*d^8
          + 670339903437583223508910817280*d^7 - 
         3566344118364856449332485226496*d^6 + 
         13989498105481528582246054625280*d^5 - 
         40140135664381531718456573952000*d^4 + 
         82026150652381035546383745024000*d^3 - 
         113089388781707217031554662400000*d^2 + 
         94259693281675671355627929600000*d - 35851858021068383490932736000000
         ,1099511627776*d^9 - 118747255799808*d^8 + 5567926883057664*d^7 - 
         148513234586959872*d^6 + 2478545499611725824*d^5 - 
         26781147765367898112*d^4 + 186882777512249655296*d^3 - 
         809766353299631505408*d^2 + 1970279800978318295040*d - 
         2042832790975256985600) + Gam(1,1)^3*M^-18*rat(177147*d^24 - 17006112
         *d^23 + 654656580*d^22 - 11727603792*d^21 + 42800617890*d^20 + 
         2317410600048*d^19 - 47024941815468*d^18 + 263501455224528*d^17 + 
         635647735738395*d^16 + 353229166149646320*d^15 - 28083606424460723952
         *d^14 + 1199621255701316268096*d^13 - 35422060569945761376864*d^12 + 
         771820415543449829170944*d^11 - 12746948952399167838719232*d^10 + 
         161873987969216090509028352*d^9 - 1592687159887034872721509632*d^8 + 
         12165943999083690456686223360*d^7 - 71879084158189251453629718528*d^6
          + 325027527964676560465777852416*d^5 - 
         1102940347982367322533999476736*d^4 + 2715784853761057850458861731840
         *d^3 - 4577036994825904570118543769600*d^2 + 
         4717856923336025394017992704000*d - 2240741126316773968183296000000,
         8589934592*d^12 - 1030792151040*d^11 + 55044300865536*d^10 - 
         1727332767236096*d^9 + 35434098667290624*d^8 - 500067783917174784*d^7
          + 4974465481965568000*d^6 - 35128763188745601024*d^5 + 
         174780022016761135104*d^4 - 597777133764698177536*d^3 + 
         1335253676170615455744*d^2 - 1751202197104880517120*d + 
         1021416395487628492800) );


id,only intbn*dala^8*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-12*miBN*rat(6561*
         d^24 - 314928*d^23 + 2458188*d^22 + 103471344*d^21 - 2062513530*d^20
          + 1899595152*d^19 + 328337031996*d^18 - 4657373125488*d^17 - 
         7701528149631*d^16 + 1037782363385952*d^15 + 311024384408304*d^14 - 
         272280267635342400*d^13 - 419034507962967072*d^12 + 
         93565660850417578752*d^11 - 641636800780932373248*d^10 - 
         18801313774505936243712*d^9 + 475519228321890628253952*d^8 - 
         5485567915258430484860928*d^7 + 39894621127373926032629760*d^6 - 
         198173700703198255004319744*d^5 + 686075209411917681763614720*d^4 - 
         1635129607765701424840704000*d^3 + 2563192756153807853322240000*d^2
          - 2381195803928950902620160000*d + 993366516508742438092800000,
         4294967296*d^8 - 377957122048*d^7 + 14190571945984*d^6 - 
         296318383685632*d^5 + 3755450684145664*d^4 - 29504844775555072*d^3 + 
         139913954146123776*d^2 - 364870734654210048*d + 398990779487354880)
          + Gam(1,1)^3*M^-12*rat(6561*d^21 - 255879*d^20 + 37179*d^19 + 
         108936657*d^18 - 1096294743*d^17 - 10052173953*d^16 + 264219753537*
         d^15 - 2109785552901*d^14 - 204942488096430*d^13 + 19350604272416940*
         d^12 - 915134464281570552*d^11 + 27746163197359561296*d^10 - 
         576605341462218471072*d^9 + 8505597508004532900672*d^8 - 
         90983714983625433174912*d^7 + 713382485412574192241664*d^6 - 
         4101373349485949851029504*d^5 + 17092504905488634206306304*d^4 - 
         50241513283849095516094464*d^3 + 98723746362042399109939200*d^2 - 
         116314142646903744076185600*d + 62085407281796402380800000,33554432*
         d^11 - 3355443200*d^10 + 147907936256*d^9 - 3789234896896*d^8 + 
         62629749981184*d^7 - 700794781302784*d^6 + 5415610162872320*d^5 - 
         28909527948591104*d^4 + 104543902031151104*d^3 - 244188888145330176*
         d^2 + 332056909634863104*d - 199495389743677440) );


id,only intbn*dala^8*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-14*miBN*rat( - 19683
         *d^24 + 944784*d^23 - 7374564*d^22 - 310414032*d^21 + 6187540590*d^20
          - 5698785456*d^19 - 985011095988*d^18 + 13972119376464*d^17 + 
         23104584448893*d^16 - 3113347090157856*d^15 - 933073153224912*d^14 + 
         816840802906027200*d^13 + 1257103523888901216*d^12 - 
         280696982551252736256*d^11 + 1924910402342797119744*d^10 + 
         56403941323517808731136*d^9 - 1426557684965671884761856*d^8 + 
         16456703745775291454582784*d^7 - 119683863382121778097889280*d^6 + 
         594521102109594765012959232*d^5 - 2058225628235753045290844160*d^4 + 
         4905388823297104274522112000*d^3 - 7689578268461423559966720000*d^2
          + 7143587411786852707860480000*d - 2980099549526227314278400000,
         34359738368*d^7 - 2748779069440*d^6 + 91534343012352*d^5 - 
         1638272325386240*d^4 + 16937426870075392*d^3 - 100539343243837440*d^2
          + 314996887218290688*d - 398990779487354880) + Gam(1,1)^3*M^-14*rat(
          - 19683*d^21 + 767637*d^20 - 111537*d^19 - 326809971*d^18 + 
         3288884229*d^17 + 30156521859*d^16 - 792659260611*d^15 + 
         6329356658703*d^14 + 614827464289290*d^13 - 58051812817250820*d^12 + 
         2745403392844711656*d^11 - 83238489592078683888*d^10 + 
         1729816024386655413216*d^9 - 25516792524013598702016*d^8 + 
         272951144950876299524736*d^7 - 2140147456237722576724992*d^6 + 
         12304120048457849553088512*d^5 - 51277514716465902618918912*d^4 + 
         150724539851547286548283392*d^3 - 296171239086127197329817600*d^2 + 
         348942427940711232228556800*d - 186256221845389207142400000,268435456
         *d^10 - 24696061952*d^9 + 985694994432*d^8 - 22428319219712*d^7 + 
         321611446091776*d^6 - 3033466681688064*d^5 + 19057147849474048*d^4 - 
         78819040792936448*d^3 + 205798889905717248*d^2 - 307119985916903424*d
          + 199495389743677440) );


id,only intbn*dala^8*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-16*miBN*rat(59049*
         d^27 - 4940433*d^26 + 146008494*d^25 - 1033462476*d^24 - 39425289570*
         d^23 + 1012151480490*d^22 - 6072106996872*d^21 - 117623197015020*d^20
          + 2526610608590769*d^19 - 7628615094716289*d^18 - 281654327691503838
         *d^17 + 868040595659949696*d^16 + 60273284079590419968*d^15 + 
         198047309798454954912*d^14 - 29221283012939305487424*d^13 + 
         243453366451877756643072*d^12 + 6764323038547227931752192*d^11 - 
         211117518328792040437969152*d^10 + 2998517940278420446590497280*d^9
          - 27622280184926081990620188672*d^8 + 180765219695460994266868334592
         *d^7 - 868587971766057462476241174528*d^6 + 
         3088857568059903844421599887360*d^5 - 8054344122124838243548240281600
         *d^4 + 14984696804486065720145215488000*d^3 - 
         18839152632254829658625802240000*d^2 + 
         14342418474814005415457587200000*d - 4991769752505118019813376000000,
         274877906944*d^9 - 29686813949952*d^8 + 1391981720764416*d^7 - 
         37128308646739968*d^6 + 619636374902931456*d^5 - 6695286941341974528*
         d^4 + 46720694378062413824*d^3 - 202441588324907876352*d^2 + 
         492569950244579573760*d - 510708197743814246400) + Gam(1,1)^3*M^-16*
         rat(59049*d^24 - 4408992*d^23 + 105264684*d^22 - 1994544*d^21 - 
         41628368394*d^20 + 642584141712*d^19 + 506348228700*d^18 - 
         127130812871760*d^17 + 13339182772809*d^16 + 201471914528811216*d^15
          - 12944112392946647376*d^14 + 511269820731094401216*d^13 - 
         14224836757461291143712*d^12 + 291211338195507495419136*d^11 - 
         4491447557729303567515392*d^10 + 52960814581768352307790848*d^9 - 
         481423340416129014838374144*d^8 + 3382432336151769197191569408*d^7 - 
         18306446744514059354469113856*d^6 + 75539303229487959818265231360*d^5
          - 233049097417579966120038432768*d^4 + 
         519820907322381900445148774400*d^3 - 790746109768971192004352409600*
         d^2 + 733035042248877328790716416000*d - 
         311985609531569876238336000000,2147483648*d^12 - 257698037760*d^11 + 
         13761075216384*d^10 - 431833191809024*d^9 + 8858524666822656*d^8 - 
         125016945979293696*d^7 + 1243616370491392000*d^6 - 
         8782190797186400256*d^5 + 43695005504190283776*d^4 - 
         149444283441174544384*d^3 + 333813419042653863936*d^2 - 
         437800549276220129280*d + 255354098871907123200) );


id,only intbn*dala^8*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-18*miBN*rat( - 
         177147*d^28 + 16120377*d^27 - 550494144*d^26 + 6587191512*d^25 + 
         89850834846*d^24 - 3941121801834*d^23 + 43026970378932*d^22 + 
         202718323098900*d^21 - 10249198889855067*d^20 + 77524537967935929*
         d^19 + 522687799451082564*d^18 - 5120569210498945812*d^17 - 
         81264140068121291520*d^16 - 1682895345501314635488*d^15 + 
         67530573490350355238016*d^14 + 28549369640342004944256*d^13 - 
         23793037762795811921726208*d^12 + 373906342676800812442717440*d^11 - 
         1432934295934014002619436032*d^10 - 31628107168976055590138151936*d^9
          + 616174503645980993928745476096*d^8 - 
         5818414529761308679792266215424*d^7 + 
         35939319551626461708408738545664*d^6 - 
         155712837812359868432988190801920*d^5 + 
         479935748648960408512008290304000*d^4 - 
         1035276720351627280617817767936000*d^3 + 
         1488991749662225367086137344000000*d^2 - 
         1283783847922391015487858278400000*d + 
         501926012294957368873058304000000,2199023255552*d^9 - 237494511599616
         *d^8 + 11135853766115328*d^7 - 297026469173919744*d^6 + 
         4957090999223451648*d^5 - 53562295530735796224*d^4 + 
         373765555024499310592*d^3 - 1619532706599263010816*d^2 + 
         3940559601956636590080*d - 4085665581950513971200) + Gam(1,1)^3*M^-18
         *rat( - 177147*d^25 + 14526054*d^24 - 416571012*d^23 + 2562411672*
         d^22 + 121385835198*d^21 - 2916619250508*d^20 + 14581193414796*d^19
          + 394847730192024*d^18 - 4324668108881787*d^17 - 362128234449983850*
         d^16 + 23138398098365675472*d^15 - 806450765758866132768*d^14 + 
         18627362990127333623520*d^13 - 275911567564209169894848*d^12 + 
         1941463134790870230326016*d^11 + 16583297364372259233040896*d^10 - 
         673548671681990394404887296*d^9 + 10131676239334797761414911488*d^8
          - 98444131828982414939977408512*d^7 + 681279650249972959885038206976
         *d^6 - 3447445043523104523986890457088*d^5 + 
         12725380017992084665017130942464*d^4 - 
         33443950957828905336305520476160*d^3 + 
         59360661004226638587641620070400*d^2 - 
         63809255800387581548068601856000*d + 31370375768434835554566144000000
         ,17179869184*d^12 - 2061584302080*d^11 + 110088601731072*d^10 - 
         3454665534472192*d^9 + 70868197334581248*d^8 - 1000135567834349568*
         d^7 + 9948930963931136000*d^6 - 70257526377491202048*d^5 + 
         349560044033522270208*d^4 - 1195554267529396355072*d^3 + 
         2670507352341230911488*d^2 - 3502404394209761034240*d + 
         2042832790975256985600) );


id,only intbn*dala^8*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-20*miBN*rat(531441*
         d^30 - 63772920*d^29 + 3151740375*d^28 - 76782595680*d^27 + 
         623145081150*d^26 + 15744728047248*d^25 - 525249952944786*d^24 + 
         5405328074644704*d^23 + 24990013613201373*d^22 - 1205348767912691640*
         d^21 + 10147224459133004427*d^20 - 18601288303595073984*d^19 + 
         96523663959795509316*d^18 + 24064270515169035108192*d^17 - 
         798047237302635321424800*d^16 - 2583336554265661469517312*d^15 + 
         427265485262225221379220864*d^14 - 5056552013089652528443462656*d^13
          - 93719477058052940807564282112*d^12 + 
         4107862737748102813062008463360*d^11 - 
         73044850339828153494856357899264*d^10 + 
         838839001784101744872798319435776*d^9 - 
         6915808226442932375636591955247104*d^8 + 
         42580903771318461915038435979558912*d^7 - 198493974343608416079486894\
         707834880*d^6 + 699443899936925414663975214514176000*d^5 - 1835644936\
         762352665915387004583936000*d^4 + 34786563800225517374556910229913600\
         00*d^3 - 4496693213486652330883691210342400000*d^2 + 3545720093514987\
         510467194257408000000*d - 1285428342978114275324681256960000000,
         17592186044416*d^10 - 2286984185774080*d^9 + 130885864170455040*d^8
          - 4336122016227655680*d^7 + 91933386568397488128*d^6 - 
         1300946380109213859840*d^5 + 12417088453605494620160*d^4 - 
         78738999337105982750720*d^3 + 316562233177123382624256*d^2 - 
         726223814599972151623680*d + 719077142423290458931200) + Gam(1,1)^3*
         M^-20*rat(531441*d^27 - 58989951*d^26 + 2611264878*d^25 - 52176877380
         *d^24 + 102392422542*d^23 + 17763632607222*d^22 - 369643773407400*
         d^21 + 1745521355233548*d^20 + 48501334791708345*d^19 - 
         1215689194199380239*d^18 + 73350260396268958530*d^17 - 
         6012682516067786887824*d^16 + 317758889974603930011648*d^15 - 
         11658561036318548204580000*d^14 + 318047169595829704854436800*d^13 - 
         6680377116743427419369723136*d^12 + 110053300059134127278970621696*
         d^11 - 1436673642953319760975684022016*d^10 + 
         14938486555050223378960066859520*d^9 - 
         123849936744401922001861475893248*d^8 + 
         816193813149755223339425965129728*d^7 - 
         4242589294029498272388437454422016*d^6 + 
         17156931377663566064118033979932672*d^5 - 
         52791489378332387576036267088936960*d^4 + 119264416799306991483371435\
         694489600*d^3 - 186292362797979760129276555296768000*d^2 + 1795377572\
         50944548707055519662080000*d - 80339271436132142207792578560000000,
         137438953472*d^13 - 19516331393024*d^12 + 1243547651014656*d^11 - 
         47012918180446208*d^10 + 1174966712743755776*d^9 - 
         20473887273561096192*d^8 + 255615307650294611968*d^7 - 
         2313072060671809552384*d^6 + 15161804994706629722112*d^5 - 
         71087001890135090397184*d^4 + 231781609903903605784576*d^3 - 
         498028529165734728695808*d^2 + 632765835708719997911040*d - 
         359538571211645229465600) );


id,only intbn*dala^9*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-14*miBN*rat(19683*
         d^27 - 1121931*d^26 + 10563210*d^25 + 531572220*d^24 - 12839292342*
         d^23 + 1995219054*d^22 + 3218392843272*d^21 - 42225998339556*d^20 - 
         202275388251909*d^19 + 11179318823525925*d^18 + 8712847184982534*d^17
          - 2551717020557187072*d^16 - 23635148933049474816*d^15 + 
         1027739405609578393056*d^14 + 8518369824438497384256*d^13 - 
         487609935875422935245568*d^12 + 2441556003671823335853312*d^11 + 
         127460690539379930811831552*d^10 - 3231608121250210102324234752*d^9
          + 40930524665819858258189438976*d^8 - 339467503040876649124788092928
         *d^7 + 1986681394016356832724997177344*d^6 - 
         8402279104654430081817219563520*d^5 + 
         25638857352663505523233652736000*d^4 - 
         55156205324558308696717787136000*d^3 + 
         79427186452354514208738508800000*d^2 - 
         68711948046943064944961126400000*d + 26984492399813926525599744000000
         ,68719476736*d^9 - 7421703487488*d^8 + 347995430191104*d^7 - 
         9282077161684992*d^6 + 154909093725732864*d^5 - 1673821735335493632*
         d^4 + 11680173594515603456*d^3 - 50610397081226969088*d^2 + 
         123142487561144893440*d - 127677049435953561600) + Gam(1,1)^3*M^-14*
         rat(19683*d^24 - 944784*d^23 + 1705860*d^22 + 565505712*d^21 - 
         7835243886*d^20 - 79133369904*d^19 + 2684002726260*d^18 - 
         16849368139248*d^17 - 408442114384893*d^16 - 16967174090244480*d^15
          + 3455852247444570480*d^14 - 214667070995876577408*d^13 + 
         8230586502715411233312*d^12 - 220503110018218924875264*d^11 + 
         4318934882671448224904448*d^10 - 63290173962541150900574208*d^9 + 
         703094470919637530992780032*d^8 - 5957672761457640314506911744*d^7 + 
         38484942322060884453132902400*d^6 - 188002542038642752639625920512*
         d^5 + 682349758177690426775767547904*d^4 - 
         1781964005573365338179669852160*d^3 + 3162405108451541648482920038400
         *d^2 - 3411347553247719354947076096000*d + 
         1686530774988370407849984000000,536870912*d^12 - 64424509440*d^11 + 
         3440268804096*d^10 - 107958297952256*d^9 + 2214631166705664*d^8 - 
         31254236494823424*d^7 + 310904092622848000*d^6 - 2195547699296600064*
         d^5 + 10923751376047570944*d^4 - 37361070860293636096*d^3 + 
         83453354760663465984*d^2 - 109450137319055032320*d + 
         63838524717976780800) );


id,only intbn*dala^9*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-16*miBN*rat( - 59049
         *d^28 + 3877551*d^27 - 60859836*d^26 - 1320073200*d^25 + 52338754746*
         d^24 - 339807258054*d^23 - 9603302834412*d^22 + 210356208943740*d^21
          - 491049792072729*d^20 - 38797116565127409*d^19 + 264523747856726448
         *d^18 + 7881685088481107100*d^17 + 4560804264661560576*d^16 - 
         3697732089088021524384*d^15 + 1166115072533546066688*d^14 + 
         1684307423061669737727360*d^13 - 20002526343776466323944704*d^12 - 
         318901615522672385703308544*d^11 + 13008802317774508508080324608*d^10
          - 206813385149965037434998420480*d^9 + 
         2082596150433946262087289692160*d^8 - 
         14786199261111863375419481948160*d^7 + 
         76860553558388567896301585301504*d^6 - 
         295375828779005698696948666859520*d^5 + 
         832078907142926069694228332544000*d^4 - 
         1672342897795579568740877991936000*d^3 + 
         2271242691902046564262084608000000*d^2 - 
         1867464126419961468145788518400000*d + 
         701596802395162089665593344000000,549755813888*d^9 - 59373627899904*
         d^8 + 2783963441528832*d^7 - 74256617293479936*d^6 + 
         1239272749805862912*d^5 - 13390573882683949056*d^4 + 
         93441388756124827648*d^3 - 404883176649815752704*d^2 + 
         985139900489159147520*d - 1021416395487628492800) + Gam(1,1)^3*M^-16*
         rat( - 59049*d^25 + 3346110*d^24 - 29681964*d^23 - 1652164776*d^22 + 
         38208880170*d^21 + 33683768676*d^20 - 10109475796284*d^19 + 
         120332175300504*d^18 + 787242771534231*d^17 + 40282027296726222*d^16
          - 10808703268680067920*d^15 + 733853371421188564704*d^14 - 
         30273103354039024712544*d^13 + 875504579125257466691904*d^12 - 
         18689885508488036721470208*d^11 + 302162828837081106549238272*d^10 - 
         3754827935784982516393269504*d^9 + 36153474528283496749333016064*d^8
          - 270354318764081301536578412544*d^7 + 
         1564616126489511253700333223936*d^6 - 6935115367537782848957576577024
         *d^5 + 23086985729340047110708965801984*d^4 - 
         55818279470262123738120176271360*d^3 + 
         92456575479483240925397149286400*d^2 - 
         93754628709405814452173930496000*d + 43849800149697630604099584000000
         ,4294967296*d^12 - 515396075520*d^11 + 27522150432768*d^10 - 
         863666383618048*d^9 + 17717049333645312*d^8 - 250033891958587392*d^7
          + 2487232740982784000*d^6 - 17564381594372800512*d^5 + 
         87390011008380567552*d^4 - 298888566882349088768*d^3 + 
         667626838085307727872*d^2 - 875601098552440258560*d + 
         510708197743814246400) );


id,only intbn*dala^9*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-18*miBN*rat(177147*
         d^30 - 17006112*d^29 + 570157461*d^28 - 3870780048*d^27 - 
         241038477270*d^26 + 6569859764448*d^25 - 33055153605558*d^24 - 
         1311783575581152*d^23 + 26183005143436287*d^22 - 48766874684823072*
         d^21 - 3889757772029602863*d^20 + 22188019425712497648*d^19 + 
         431047373815411463052*d^18 + 7268102323547812643808*d^17 - 
         268533723431481597812448*d^16 - 3370414690948323143089152*d^15 + 
         177392786352278519906750592*d^14 - 1099410364219350778082491392*d^13
          - 52776301849017055129151417088*d^12 + 
         1601100551285134769386624094208*d^11 - 
         24345798182807118557505199211520*d^10 + 
         248125171921977680418383270830080*d^9 - 
         1840287322174961418769471984533504*d^8 + 
         10258968049486481241080685706346496*d^7 - 
         43451589283780730068730767704588288*d^6 + 139413130313140744551085774\
         460682240*d^5 - 333628929012360738209294331150336000*d^4 + 5771902163\
         36631299508955263270912000*d^3 - 681899281657091732101374620467200000
         *d^2 + 492022601840227737125188848844800000*d - 163459179834300601009\
         918967808000000,4398046511104*d^10 - 571746046443520*d^9 + 
         32721466042613760*d^8 - 1084030504056913920*d^7 + 
         22983346642099372032*d^6 - 325236595027303464960*d^5 + 
         3104272113401373655040*d^4 - 19684749834276495687680*d^3 + 
         79140558294280845656064*d^2 - 181555953649993037905920*d + 
         179769285605822614732800) + Gam(1,1)^3*M^-18*rat(177147*d^27 - 
         15411789*d^26 + 428262714*d^25 + 275168340*d^24 - 247316566950*d^23
          + 4360904187138*d^22 + 10875991860936*d^21 - 1308083154478812*d^20
          + 14389636195639827*d^19 - 92958095245241853*d^18 + 
         30857739664066348374*d^17 - 2750752767506296646736*d^16 + 
         137491522565684507172864*d^15 - 4759049728189920059925984*d^14 + 
         122472848478938153848840512*d^13 - 2419068795156692467022211840*d^12
          + 37311526751469354353009362176*d^11 - 
         453970065189016931317798452480*d^10 + 4380701589957948385133617251840
         *d^9 - 33571128029869622824033960894464*d^8 + 
         203735391171549674987592999714816*d^7 - 
         971742094966625859401160273887232*d^6 + 
         3593317022475846134080522022682624*d^5 - 
         10075507653623081427478378067263488*d^4 + 
         20671743629933895817451098205061120*d^3 - 
         29223603116552286244997096590540800*d^2 + 
         25401743002601103791662030651392000*d - 
         10216198739643787563119935488000000,34359738368*d^13 - 4879082848256*
         d^12 + 310886912753664*d^11 - 11753229545111552*d^10 + 
         293741678185938944*d^9 - 5118471818390274048*d^8 + 
         63903826912573652992*d^7 - 578268015167952388096*d^6 + 
         3790451248676657430528*d^5 - 17771750472533772599296*d^4 + 
         57945402475975901446144*d^3 - 124507132291433682173952*d^2 + 
         158191458927179999477760*d - 89884642802911307366400) );


id,only intbn*dala^9*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-20*miBN*rat( - 
         531441*d^31 + 55269864*d^30 - 2131373655*d^29 + 26354749680*d^28 + 
         605376449730*d^27 - 25715049345648*d^26 + 273334304188818*d^25 + 
         2998671172471872*d^24 - 111475262807516637*d^23 + 805508550101469672*
         d^22 + 9138355827470061813*d^21 - 143754303042532996848*d^20 + 
         201096948897725674428*d^19 - 25608649138525763257248*d^18 + 
         413018909059930759693728*d^17 + 15352092351107826612314112*d^16 - 
         385932100393974637866943872*d^15 - 1779695751105951013624071168*d^14
          + 174624309267487381262659684608*d^13 - 
         2608351104819255760140979949568*d^12 + 
         7319046535858508485864222485504*d^11 + 
         329878603653148711044903406952448*d^10 - 
         6505615802102695542328181155725312*d^9 + 
         68072027851768456095147035304394752*d^8 - 482800485997486974561128080\
         965107712*d^7 + 2476459689560809242607815100811182080*d^6 - 935545746\
         2228453968708216427642880000*d^5 + 2589166260817509091719050105035161\
         6000*d^4 - 51161808866874175468407365157519360000*d^3 + 6840137132227\
         1449783671865108070400000*d^2 - 55446093153261685892150426861568000000
         *d + 20566853487649828405194900111360000000,35184372088832*d^10 - 
         4573968371548160*d^9 + 261771728340910080*d^8 - 8672244032455311360*
         d^7 + 183866773136794976256*d^6 - 2601892760218427719680*d^5 + 
         24834176907210989240320*d^4 - 157477998674211965501440*d^3 + 
         633124466354246765248512*d^2 - 1452447629199944303247360*d + 
         1438154284846580917862400) + Gam(1,1)^3*M^-20*rat( - 531441*d^28 + 
         50486895*d^27 - 1667425662*d^26 + 10396639332*d^25 + 732437615538*
         d^24 - 19401911367894*d^23 + 85425651691848*d^22 + 4168779019284852*
         d^21 - 76429676475445113*d^20 + 439667837532046719*d^19 - 
         53899233289078874706*d^18 + 4839078349727483551344*d^17 - 
         221555969717519339806464*d^16 + 6574418796724885324393632*d^15 - 
         131510193014732933581156800*d^14 + 1591622403210152141698734336*d^13
          - 3167266191239288569055051520*d^12 - 324179157992826275487845925120
         *d^11 + 8048291732202892796650877492736*d^10 - 
         115165848136401652061499593859072*d^9 + 
         1165405174760675528690357649162240*d^8 - 
         8816511716366585301042377987653632*d^7 + 
         50724497326808406294096965290819584*d^6 - 221719412664284669449852276\
         589985792*d^5 + 725399413254011209733208837728501760*d^4 - 1721938305\
         990932103604666415815065600*d^3 + 28011400475167316133613693650862080\
         00*d^2 - 2792264844578980637105095736033280000*d + 128542834297811427\
         5324681256960000000,274877906944*d^13 - 39032662786048*d^12 + 
         2487095302029312*d^11 - 94025836360892416*d^10 + 2349933425487511552*
         d^9 - 40947774547122192384*d^8 + 511230615300589223936*d^7 - 
         4626144121343619104768*d^6 + 30323609989413259444224*d^5 - 
         142174003780270180794368*d^4 + 463563219807807211569152*d^3 - 
         996057058331469457391616*d^2 + 1265531671417439995822080*d - 
         719077142423290458931200) );


id,only intbn*dala^9*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-22*miBN*rat(1594323*
         d^33 - 216296487*d^32 + 11992320459*d^31 - 318469266945*d^30 + 
         2141793913932*d^29 + 116300920724262*d^28 - 3678669422062434*d^27 + 
         34449324253709310*d^26 + 439646641607371995*d^25 - 
         14832189774550164555*d^24 + 120603738142952598735*d^23 + 
         492696863101210209075*d^22 - 18449922693041230104810*d^21 + 
         424868402282752079093340*d^20 - 5419530369703998307574280*d^19 - 
         159799549823119947062808000*d^18 + 4812033880194816962123478720*d^17
          + 42766885240756286874047061120*d^16 - 
         3296506008689608885342863690240*d^15 + 
         33518906216649162194989321440000*d^14 + 
         944690564166498891743847025946112*d^13 - 
         39230147116691074138122571081724928*d^12 + 73970040594341906249431811\
         7104695296*d^11 - 9285271713714186165251604421891399680*d^10 + 
         85446343340344329242672420169373974528*d^9 - 598807443530847393005708\
         451110494666752*d^8 + 3244873503628576749247738245385658302464*d^7 - 
         13627195375927998083884927094513215733760*d^6 + 440007000454324192844\
         60982688531297075200*d^5 - 107201731090699492832437318082705227776000
         *d^4 + 190566478044216319524599251349966684160000*d^3 - 2330542617840\
         04792821827657208850022400000*d^2 + 175107834329880291816996646114295\
         808000000*d - 60860354971376857373237623746723840000000,
         281474976710656*d^11 - 43347146413441024*d^10 + 2972375754064527360*
         d^9 - 119638124101097226240*d^8 + 3136005039325779591168*d^7 - 
         56117562524012057198592*d^6 + 698236825219626036101120*d^5 - 
         6027985955578205658152960*d^4 + 35300771476282671498264576*d^3 - 
         133179478573614933353693184*d^2 + 290375179085161953566392320*d - 
         276125622690543536229580800) + Gam(1,1)^3*M^-22*rat(1594323*d^30 - 
         201947580*d^29 + 10146094425*d^28 - 223391814840*d^27 - 65834320590*
         d^26 + 120362638902408*d^25 - 2604942425799198*d^24 + 
         8735356587916800*d^23 + 572964811658606007*d^22 - 9951921433323855612
         *d^21 - 43444567950340432491*d^20 + 14554366602834074473944*d^19 - 
         1395567425071440545211180*d^18 + 88367541143158549820695824*d^17 - 
         3942303036783843040738416864*d^16 + 131527809704549845501533255168*
         d^15 - 3409717565271129539916731859072*d^14 + 
         70247151559868727353934935646720*d^13 - 
         1165367291133726260783791175077632*d^12 + 
         15684953261070738482195751537887232*d^11 - 17189986585225425346766242\
         4654380032*d^10 + 1534812774689594750226439356157464576*d^9 - 
         11135533118106326834005340700906455040*d^8 + 652619926099715802790210\
         76439379869696*d^7 - 305855945624184669019046465727731859456*d^6 + 
         1128432070181443841443775587733197553664*d^5 - 3200298812206958289114\
         078754251158323200*d^4 + 6723526932897599623223891011716513792000*d^3
          - 9843129295299006104537253967728476160000*d^2 + 8952392911345809753\
         113276444741468160000*d - 3803772185711053585827351484170240000000,
         2199023255552*d^14 - 365037860421632*d^13 + 27391033671155712*d^12 - 
         1229728988876767232*d^11 + 36852427985191436288*d^10 - 
         778769414070579757056*d^9 + 11951817635452174729216*d^8 - 
         135165431108462083833856*d^7 + 1130808551213280943669248*d^6 - 
         6959525148209507259645952*d^5 + 31005914484274332405071872*d^4 - 
         96972594669750740280410112*d^3 + 201367208570981655785766912*d^2 - 
         248734698051534802869288960*d + 138062811345271768114790400) );


id,only intbn*dala^10*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-16*miBN*rat(59049*
         d^30 - 3936600*d^29 + 44188335*d^28 + 2570337360*d^27 - 73881457074*
         d^26 - 33014660400*d^25 + 26975140722270*d^24 - 358075613243520*d^23
          - 2878447325248395*d^22 + 121160300157472680*d^21 - 
         37364788377945405*d^20 - 24272430617592331920*d^19 - 
         255515633029927999260*d^18 + 7628030404823419454880*d^17 + 
         219143297124730645665120*d^16 - 4763578160234637517086720*d^15 - 
         83510635732959490240302720*d^14 + 3148500426174013906319201280*d^13
          - 9762583220744105287682945280*d^12 - 
         1063119869762866384123351019520*d^11 + 
         27351205910980466575567935845376*d^10 - 
         377170100677721408098665503416320*d^9 + 
         3514311730681145630504451648552960*d^8 - 
         23720625889122547783013509292359680*d^7 + 118973577949007423793920509\
         089153024*d^6 - 445298565367078329501455391412715520*d^5 + 1229511074\
         584391659777942735552512000*d^4 - 24331352763430865078783218506792960\
         00*d^3 + 3264729737131237963470378054451200000*d^2 - 2658850759729759\
         497737380193894400000*d + 991394907155680530343976239104000000,
         1099511627776*d^10 - 142936511610880*d^9 + 8180366510653440*d^8 - 
         271007626014228480*d^7 + 5745836660524843008*d^6 - 
         81309148756825866240*d^5 + 776068028350343413760*d^4 - 
         4921187458569123921920*d^3 + 19785139573570211414016*d^2 - 
         45388988412498259476480*d + 44942321401455653683200) + Gam(1,1)^3*
         M^-16*rat(59049*d^27 - 3405159*d^26 + 12479022*d^25 + 2748665340*d^24
          - 49578148386*d^23 - 530004702330*d^22 + 23287339051704*d^21 - 
         140674524999444*d^20 - 4614713804295567*d^19 + 83393453015437833*d^18
          - 3190461848470734750*d^17 + 643611599582282832528*d^16 - 
         52146357465274269228288*d^15 + 2511069881172464524224096*d^14 - 
         84174578436890283358330944*d^13 + 2092810774647861209220468480*d^12
          - 39765370869484591609255958784*d^11 + 
         586798693696647721532229306624*d^10 - 6784671452309799447332941106688
         *d^9 + 61696057829872080710041441247232*d^8 - 
         440838657060540007680625857675264*d^7 + 
         2460290240050967501413755502657536*d^6 - 
         10593218942954848778774635770937344*d^5 + 
         34455415519244098799554066712100864*d^4 - 
         81771200911486605929881933719797760*d^3 + 133449535473167708675271595\
         288166400*d^2 - 133731692467179524272320336101376000*d + 
         61962181697230033146498514944000000,8589934592*d^13 - 1219770712064*
         d^12 + 77721728188416*d^11 - 2938307386277888*d^10 + 
         73435419546484736*d^9 - 1279617954597568512*d^8 + 
         15975956728143413248*d^7 - 144567003791988097024*d^6 + 
         947612812169164357632*d^5 - 4442937618133443149824*d^4 + 
         14486350618993975361536*d^3 - 31126783072858420543488*d^2 + 
         39547864731794999869440*d - 22471160700727826841600) );


id,only intbn*dala^10*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-18*miBN*rat( - 
         177147*d^31 + 13463172*d^30 - 242789805*d^29 - 6473738700*d^28 + 
         293613817302*d^27 - 1969636816872*d^26 - 81849832658010*d^25 + 
         1829530779954120*d^24 - 1390775195073375*d^23 - 444077425579373100*
         d^22 + 3504582769543071255*d^21 + 71771077778194524420*d^20 + 
         86918841797198704020*d^19 - 30038528939308242343920*d^18 - 
         443845040039136192258720*d^17 + 20426746800196370629883520*d^16 + 
         117151718712308620242480000*d^15 - 11783799079044907445686080000*d^14
          + 117445761595104705239986471680*d^13 + 
         2916007279107764204314930590720*d^12 - 
         111820974086301658482157636082688*d^11 + 
         1897344067540617288411898713919488*d^10 - 
         21103698011019636318275989041315840*d^9 + 169562606126439721003165174\
         036561920*d^8 - 1021098258742453609306139787453530112*d^7 + 466715587\
         8673442854734140428734431232*d^6 - 1615689305403136820537457916621357\
         0560*d^5 + 41725715917392225997417362147508224000*d^4 - 7792197694900\
         0136111004145982373888000*d^3 + 99388984918863941470382726106316800000
         *d^2 - 77422005993900307527678574146355200000*d + 2775905740035905484\
         9631334694912000000,8796093022208*d^10 - 1143492092887040*d^9 + 
         65442932085227520*d^8 - 2168061008113827840*d^7 + 
         45966693284198744064*d^6 - 650473190054606929920*d^5 + 
         6208544226802747310080*d^4 - 39369499668552991375360*d^3 + 
         158281116588561691312128*d^2 - 363111907299986075811840*d + 
         359538571211645229465600) + Gam(1,1)^3*M^-18*rat( - 177147*d^28 + 
         11868849*d^27 - 132781518*d^26 - 7896583404*d^25 + 225697074678*d^24
          + 201825952182*d^23 - 84702148820352*d^22 + 1074069068446044*d^21 + 
         9905254712902269*d^20 - 379392345566589375*d^19 + 
         11906402229844463574*d^18 - 2020167730504029070584*d^17 + 
         174460197184126726995648*d^16 - 8993307652545073111064352*d^15 + 
         322833691983499856753267520*d^14 - 8635320520176511561694671872*d^13
          + 177894814298593888685940993792*d^12 - 
         2873826465435511729655854765824*d^11 + 
         36784377780435534544901243905536*d^10 - 
         375058974154290626655446674728960*d^9 + 
         3050005590418038282923037927948288*d^8 - 
         19724353117848022719298790522880000*d^7 + 100667783550291636375909061\
         387223040*d^6 - 399976376960468062204352001722548224*d^5 + 1210065237\
         273294584177159669098217472*d^4 - 26899422319411280920625089300188364\
         80*d^3 + 4137782070650234415724565676372787200*d^2 - 3930373934172716\
         779064464955670528000*d + 1734941087522440928101958418432000000,
         68719476736*d^13 - 9758165696512*d^12 + 621773825507328*d^11 - 
         23506459090223104*d^10 + 587483356371877888*d^9 - 
         10236943636780548096*d^8 + 127807653825147305984*d^7 - 
         1156536030335904776192*d^6 + 7580902497353314861056*d^5 - 
         35543500945067545198592*d^4 + 115890804951951802892288*d^3 - 
         249014264582867364347904*d^2 + 316382917854359998955520*d - 
         179769285605822614732800) );


id,only intbn*dala^10*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-20*miBN*rat(531441*
         d^33 - 57927069*d^32 + 2183454873*d^31 - 13947078555*d^30 - 
         1352941159356*d^29 + 39490641377154*d^28 - 153047563642998*d^27 - 
         12263438275013430*d^26 + 241518270617384265*d^25 - 52240032934279785*
         d^24 - 52623150151865639355*d^23 + 436068271032809011425*d^22 + 
         3532021579701309441330*d^21 + 44257461509756890396980*d^20 - 
         994910765579434650697560*d^19 - 80669889631294624430740800*d^18 + 
         1458954265833752630583435840*d^17 + 34372490195467094999681888640*
         d^16 - 1298373327119188036319156743680*d^15 + 
         5377315211042248582606110009600*d^14 + 
         496212442115205041526477338786304*d^13 - 
         14941234986937575292154130127911936*d^12 + 24245422900000516730119505\
         0573555712*d^11 - 2710485176843974998886378462694891520*d^10 + 
         22503266257388559104596849842009931776*d^9 - 143161703725217666446045\
         298949616041984*d^8 + 706593265116179006227249040622585643008*d^7 - 
         2707985391664322673527088262248105246720*d^6 + 7988962903099926797055\
         322354476279398400*d^5 - 17798931391217794744652988951142858752000*
         d^4 + 28954565104678728880854453679627960320000*d^3 - 324293137976608\
         91455268547350809804800000*d^2 + 223353290655352895438840781770588160\
         00000*d - 7123998758363062278421138510970880000000,70368744177664*
         d^11 - 10836786603360256*d^10 + 743093938516131840*d^9 - 
         29909531025274306560*d^8 + 784001259831444897792*d^7 - 
         14029390631003014299648*d^6 + 174559206304906509025280*d^5 - 
         1506996488894551414538240*d^4 + 8825192869070667874566144*d^3 - 
         33294869643403733338423296*d^2 + 72593794771290488391598080*d - 
         69031405672635884057395200) + Gam(1,1)^3*M^-20*rat(531441*d^30 - 
         53144100*d^29 + 1695592035*d^28 + 2312358840*d^27 - 1366340912730*
         d^26 + 27243059058936*d^25 + 117830235155094*d^24 - 11783455124595840
         *d^23 + 134493247121388669*d^22 + 1388836133912889756*d^21 - 
         75120991261919595657*d^20 + 6570155838007304716968*d^19 - 
         628726740321556889584740*d^18 + 38461016384012288474597808*d^17 - 
         1631214667459823666372772768*d^16 + 51557169745146753336953857536*
         d^15 - 1263182309372076314576150553984*d^14 + 
         24518815956443787517173636564480*d^13 - 
         381838036527159864265783066662144*d^12 + 
         4806221979716093536509005720745984*d^11 - 
         49077959291837522377724257453200384*d^10 + 40682064409561651199104384\
         1338355712*d^9 - 2730878081752167664945487925167554560*d^8 + 14758776\
         776134354665830651159939776512*d^7 - 63575785914523251792797370080719\
         798272*d^6 + 214899614935232126440467503220714897408*d^5 - 5565922914\
         26660732712648157515743232000*d^4 + 106443806053122381841134064804823\
         0400000*d^3 - 1413853072166854318638339042626764800000*d^2 + 
         1162803380771972595948506043666923520000*d - 445249922397691392401321\
         156935680000000,549755813888*d^14 - 91259465105408*d^13 + 
         6847758417788928*d^12 - 307432247219191808*d^11 + 9213106996297859072
         *d^10 - 194692353517644939264*d^9 + 2987954408863043682304*d^8 - 
         33791357777115520958464*d^7 + 282702137803320235917312*d^6 - 
         1739881287052376814911488*d^5 + 7751478621068583101267968*d^4 - 
         24243148667437685070102528*d^3 + 50341802142745413946441728*d^2 - 
         62183674512883700717322240*d + 34515702836317942028697600) );


id,only intbn*dala^10*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-22*miBN*rat( - 
         1594323*d^34 + 187598673*d^33 - 8098983693*d^32 + 102607498683*d^31
          + 3590652891078*d^30 - 154853211175038*d^29 + 1585252849025718*d^28
          + 31766725343414502*d^27 - 1059734478174139575*d^26 + 
         6918550225617468645*d^25 + 146375677798950363255*d^24 - 
         2663564149674356986305*d^23 + 9581379157219446341460*d^22 - 
         92769793808009937206760*d^21 - 2228100871385539116105840*d^20 + 
         257351096477791916599145040*d^19 - 1935641983378657914992934720*d^18
          - 129383495084262992192269678080*d^17 + 
         2526702074355995721610016590080*d^16 + 
         25818201939763797741182224984320*d^15 - 
         1548030876066183811253654811866112*d^14 + 
         22225716961694094086733324614694912*d^13 - 
         33557757842979728008111837633646592*d^12 - 40293355932673569596461216\
         85993115648*d^11 + 81688547506511021731856459424671219712*d^10 - 
         939226736595350533362395111938236874752*d^9 + 75336604799266763248550\
         13874603245699072*d^8 - 44780527689386383402574361322428633710592*d^7
          + 201288816721271546225467705012706586132480*d^6 - 68481086972708405\
         4287860370310858119577600*d^5 + 1739064681588374551459272474138727415\
         808000*d^4 - 3197142343011888958620958867090550292480000*d^3 + 
         4019868877782205978975901183645004595200000*d^2 - 3091080662966468395\
         332702006310600704000000*d + 1095486389484783432718277227441029120000\
         000,562949953421312*d^11 - 86694292826882048*d^10 + 
         5944751508129054720*d^9 - 239276248202194452480*d^8 + 
         6272010078651559182336*d^7 - 112235125048024114397184*d^6 + 
         1396473650439252072202240*d^5 - 12055971911156411316305920*d^4 + 
         70601542952565342996529152*d^3 - 266358957147229866707386368*d^2 + 
         580750358170323907132784640*d - 552251245381087072459161600) + Gam(1,
         1)^3*M^-22*rat( - 1594323*d^31 + 173249766*d^30 - 6511037985*d^29 + 
         40762115190*d^28 + 4086886987710*d^27 - 119177621131788*d^26 + 
         438414925555854*d^25 + 38153607076468764*d^24 - 730201230241108407*
         d^23 - 361445176531052514*d^22 + 222579153750169833507*d^21 - 
         13772364379727946689106*d^20 + 1133588826220427204680188*d^19 - 
         63247327491872620006894584*d^18 + 2351687296206989143965892032*d^17
          - 60566355042440670768241751616*d^16 + 
         1042216990589232320889133266048*d^15 - 
         8872235384988395635433762183424*d^14 - 
         99081436943910831587037666563328*d^13 + 
         5291657979336334211912489613510144*d^12 - 110429292847019039211861103\
         027590144*d^11 + 1559384810650981812191484287621376000*d^10 - 
         16491096826306378670070567709927907328*d^9 + 135177603515942302733075\
         056176936321024*d^8 - 868859921355303776003332910181105795072*d^7 + 
         4376974951053880200899060795365975916544*d^6 - 1711147845105903085687\
         3881824946397642752*d^5 + 50881851686827649580829526564804336025600*
         d^4 - 111180355496857787113492784243168772096000*d^3 + 16822393440403\
         6300128557294974371102720000*d^2 - 1573393002185135219702116245211761\
         86880000*d + 68467899342798964544892326715064320000000,4398046511104*
         d^14 - 730075720843264*d^13 + 54782067342311424*d^12 - 
         2459457977753534464*d^11 + 73704855970382872576*d^10 - 
         1557538828141159514112*d^9 + 23903635270904349458432*d^8 - 
         270330862216924167667712*d^7 + 2261617102426561887338496*d^6 - 
         13919050296419014519291904*d^5 + 62011828968548664810143744*d^4 - 
         193945189339501480560820224*d^3 + 402734417141963311571533824*d^2 - 
         497469396103069605738577920*d + 276125622690543536229580800) );


id,only intbn*dala^10*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-24*miBN*rat(4782969
         *d^36 - 727011288*d^35 + 44837677170*d^34 - 1287945137736*d^33 + 
         6210962717751*d^32 + 751284212749488*d^31 - 23550796575886140*d^30 + 
         190493107327161936*d^29 + 5222011774804439175*d^28 - 
         153886344195227216760*d^27 + 1092708197096891681250*d^26 + 
         15947812823492931534360*d^25 - 417610806559803984842055*d^24 + 
         5121150773414311962214560*d^23 - 2564719884622589483532600*d^22 - 
         2586884247508748217874066080*d^21 + 31517346242836145021877146160*
         d^20 + 1407293234010384263370011479680*d^19 - 
         32557569852309910533394809408000*d^18 - 
         579890418489353094637026647946240*d^17 + 
         30600198516316355042032093028564736*d^16 - 25962896450839075297674597\
         5392899072*d^15 - 11201066937131357390048358159790725120*d^14 + 
         451265013107736507177523789248305750016*d^13 - 8996167370541668853727\
         024319812767952896*d^12 + 122547091136862651648086780392373547466752*
         d^11 - 1246006814063429422192798621964590890024960*d^10 + 98082520795\
         78395098791815478068483070623744*d^9 - 607417841900783355618017271580\
         43553001635840*d^8 + 297330904309497346301279140066138735132016640*
         d^7 - 1146187394382909449776616716541352541264281600*d^6 + 3439372213\
         670062474271745884568961501102080000*d^5 - 78647233707107728424345104\
         09208741725470720000*d^4 + 132288894420654093932704582216458134618112\
         00000*d^3 - 15413513495349106255801730184230184419328000000*d^2 + 
         11098402950037034090449460775416192040960000000*d - 37153689574482942\
         35744224407574491955200000000,4503599627370496*d^12 - 
         810647932926689280*d^11 + 65590424973023903744*d^10 - 
         3150718299308399001600*d^9 + 99945540255268919574528*d^8 - 
         2202459096743717225103360*d^7 + 34516695213503032372232192*d^6 - 
         386914294580615721548513280*d^5 + 3072454501141056297763864576*d^4 - 
         16815992591311430276937154560*d^3 + 60048665951986403532198641664*d^2
          - 125214084462476069263292497920*d + 114868259039266111071505612800)
          + Gam(1,1)^3*M^-24*rat(4782969*d^33 - 683964567*d^32 + 38595902625*
         d^31 - 927888014385*d^30 - 2884422009060*d^29 + 744496035014574*d^28
          - 16844282969888646*d^27 + 24834458975220702*d^26 + 
         5797521571416316641*d^25 - 102941612234349106059*d^24 + 
         56916126411401590605*d^23 + 7445149776298435375971*d^22 + 
         2576823280776431071307514*d^21 - 336602902567551890089859388*d^20 + 
         25639252895841652464610010472*d^19 - 1371060330920522607318674842080*
         d^18 + 55006060537807804646318071461312*d^17 - 
         1723764430509197527464094879512192*d^16 + 
         43284615605606508438452975583588864*d^15 - 88491318401928192390214688\
         5168179456*d^14 + 14874030689017733168920850338826973696*d^13 - 
         206720943700116408170611244265644534784*d^12 + 2381970690676408361458\
         760074419167471616*d^11 - 22760175466553532792382406991874004140032*
         d^10 + 179955900183100952191061281937640903180288*d^9 - 1171835670874\
         423856608160239663596886032384*d^8 + 62363016190932960049779804719370\
         83219902464*d^7 - 26811299623419606732673442108140808275230720*d^6 + 
         91557080876424454797783980545779468746096640*d^5 - 242239949694627769\
         946803939256112979850035200*d^4 + 47810699355535663738350413071437352\
         0883712000*d^3 - 661662100754500898482835849113134862172160000*d^2 + 
         572053076388838522184400267547392973209600000*d - 2322105598405183897\
         34014025473405747200000000,35184372088832*d^15 - 6755399441055744*
         d^14 + 590112288673890304*d^13 - 31070333829229051904*d^12 + 
         1101206107135798149120*d^11 - 27790920666968913608704*d^10 + 
         515197158420595974602752*d^9 - 7134603034083498028695552*d^8 + 
         74321756160532721973592064*d^7 - 581768759676076988720742400*d^6 + 
         3391257093403544338493865984*d^5 - 14450021940174134124996460544*d^4
          + 43562474719752014449222877184*d^3 - 87748513934352925652787658752*
         d^2 + 105682639370962826283460853760*d - 
         57434129519633055535752806400) );


id,only intbn*dala^11*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-18*miBN*rat(177147*
         d^33 - 13640319*d^32 + 180788355*d^31 + 11838714327*d^30 - 
         399689637732*d^29 - 353312120106*d^28 + 201931233586830*d^27 - 
         2842731613548882*d^26 - 32409682889246685*d^25 + 1268838781069378365*
         d^24 - 1857245214661054425*d^23 - 263400318706881119445*d^22 - 
         1205384956084555419690*d^21 + 45876785977207449001980*d^20 + 
         2572861378109654208871800*d^19 - 22193730873241968234997440*d^18 - 
         1840077725596082788633139520*d^17 + 25155329228488814933972211840*
         d^16 + 824861328706488453010731302400*d^15 - 
         24716304067458222353130676542720*d^14 + 
         24075999980684996523957706157568*d^13 + 
         10732410077710482795521526362532864*d^12 - 28218589239957324558401023\
         9592417280*d^11 + 4207510722322785873451977950380351488*d^10 - 
         43524348946705121228012328935545110528*d^9 + 333325485333094355478254\
         878062530789376*d^8 - 1939614742992217220790679987851354439680*d^7 + 
         8641510969508368141605364053756107292672*d^6 - 2933490408070261463357\
         2554898961659330560*d^5 + 74610609812412563676215058086385156096000*
         d^4 - 137668074734603082225957702853761957888000*d^3 + 17392595119789\
         9087011259550933817753600000*d^2 - 1344560207546755791929243616078200\
         83200000*d + 47914610996207014069416395865587712000000,17592186044416
         *d^11 - 2709196650840064*d^10 + 185773484629032960*d^9 - 
         7477382756318576640*d^8 + 196000314957861224448*d^7 - 
         3507347657750753574912*d^6 + 43639801576226627256320*d^5 - 
         376749122223637853634560*d^4 + 2206298217267666968641536*d^3 - 
         8323717410850933334605824*d^2 + 18148448692822622097899520*d - 
         17257851418158971014348800) + Gam(1,1)^3*M^-18*rat(177147*d^30 - 
         12045996*d^29 + 69185745*d^28 + 12692385720*d^27 - 287480122110*d^26
          - 3172421486808*d^25 + 179472931452882*d^24 - 1182219198748704*d^23
          - 46619228648243457*d^22 + 880290637028208468*d^21 + 
         6911894104111440045*d^20 - 922132684578486711672*d^19 + 
         137189384703227533184532*d^18 - 13245969090163287215610864*d^17 + 
         786497250813449806364007456*d^16 - 32417123653840292566879610880*d^15
          + 995459355435585230476854892416*d^14 - 
         23644788904897095194754363588096*d^13 + 
         443439002920206501700102482446592*d^12 - 
         6643052534821999906228337927129088*d^11 + 
         79997430137468858844297251571975168*d^10 - 77628174690878338165847292\
         6964150272*d^9 + 6063930499712268375259153378526527488*d^8 - 37953130\
         000315669574781725849687556096*d^7 + 18861295947644307210178237684300\
         9327104*d^6 - 733341040176056125593599420130558738432*d^5 + 217987111\
         5173936874544545230218479009792*d^4 - 4777008113717235383840987505034\
         033889280*d^3 + 7263615045057725839761396094572114739200*d^2 - 
         6835345759285825807417566090208739328000*d + 299466318726293837933852\
         4741599232000000,137438953472*d^14 - 22814866276352*d^13 + 
         1711939604447232*d^12 - 76858061804797952*d^11 + 2303276749074464768*
         d^10 - 48673088379411234816*d^9 + 746988602215760920576*d^8 - 
         8447839444278880239616*d^7 + 70675534450830058979328*d^6 - 
         434970321763094203727872*d^5 + 1937869655267145775316992*d^4 - 
         6060787166859421267525632*d^3 + 12585450535686353486610432*d^2 - 
         15545918628220925179330560*d + 8628925709079485507174400) );


id,only intbn*dala^11*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-20*miBN*rat( - 
         531441*d^33 + 40920957*d^32 - 542365065*d^31 - 35516142981*d^30 + 
         1199068913196*d^29 + 1059936360318*d^28 - 605793700760490*d^27 + 
         8528194840646646*d^26 + 97229048667740055*d^25 - 3806516343208135095*
         d^24 + 5571735643983163275*d^23 + 790200956120643358335*d^22 + 
         3616154868253666259070*d^21 - 137630357931622347005940*d^20 - 
         7718584134328962626615400*d^19 + 66581192619725904704992320*d^18 + 
         5520233176788248365899418560*d^17 - 75465987685466444801916635520*
         d^16 - 2474583986119465359032193907200*d^15 + 
         74148912202374667059392029628160*d^14 - 
         72227999942054989571873118472704*d^13 - 
         32197230233131448386564579087598592*d^12 + 84655767719871973675203071\
         8777251840*d^11 - 12622532166968357620355933851141054464*d^10 + 
         130573046840115363684036986806635331584*d^9 - 99997645599928306643476\
         4634187592368128*d^8 + 5818844228976651662372039963554063319040*d^7
          - 25924532908525104424816092161268321878016*d^6 + 880047122421078439\
         00717664696884977991680*d^5 - 223831829437237691028645174259155468288\
         000*d^4 + 413004224203809246677873108561285873664000*d^3 - 5217778535\
         93697261033778652801453260800000*d^2 + 403368062264026737578773084823\
         460249600000*d - 143743832988621042208249187596763136000000,
         140737488355328*d^10 - 20266198323167232*d^9 + 1283525893800591360*
         d^8 - 46983803112542699520*d^7 + 1098164488537462800384*d^6 - 
         17077136376631400595456*d^5 + 178347048843499012096000*d^4 - 
         1230522489354112708116480*d^3 + 5345160844600208667967488*d^2 - 
         13138130840805379997171712*d + 13806281134527176811479040) + Gam(1,1)
         ^3*M^-20*rat( - 531441*d^30 + 36137988*d^29 - 207557235*d^28 - 
         38077157160*d^27 + 862440366330*d^26 + 9517264460424*d^25 - 
         538418794358646*d^24 + 3546657596246112*d^23 + 139857685944730371*
         d^22 - 2640871911084625404*d^21 - 20735682312334320135*d^20 + 
         2766398053735460135016*d^19 - 411568154109682599553596*d^18 + 
         39737907270489861646832592*d^17 - 2359491752440349419092022368*d^16
          + 97251370961520877700638832640*d^15 - 
         2986378066306755691430564677248*d^14 + 
         70934366714691285584263090764288*d^13 - 
         1330317008760619505100307447339776*d^12 + 
         19929157604465999718685013781387264*d^11 - 23999229041240657653289175\
         4715925504*d^10 + 2328845240726350144975418780892450816*d^9 - 
         18191791499136805125777460135579582464*d^8 + 113859390000947008724345\
         177549062668288*d^7 - 565838878429329216305347130529027981312*d^6 + 
         2200023120528168376780798260391676215296*d^5 - 6539613345521810623633\
         635690655437029376*d^4 + 14331024341151706151522962515102101667840*
         d^3 - 21790845135173177519284188283716344217600*d^2 + 205060372778574\
         77422252698270626217984000*d - 89839895617888151380155742247976960000\
         00,1099511627776*d^13 - 171523813933056*d^12 + 11980278696247296*d^11
          - 495061707475910656*d^10 + 13475596917836611584*d^9 - 
         254628737856923762688*d^8 + 3429621439156849737728*d^7 - 
         33286501162662544539648*d^6 + 232539263980015026438144*d^5 - 
         1154369934304603365441536*d^4 + 3959257899091132548120576*d^3 - 
         8893718343964044658999296*d^2 + 11746420845850381302890496*d - 
         6903140567263588405739520) );


id,only intbn*dala^11*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-22*miBN*rat(1594323
         *d^36 - 195570288*d^35 + 8227060974*d^34 - 48669602976*d^33 - 
         7116323951979*d^32 + 223357571773776*d^31 - 569433392565444*d^30 - 
         101104430185196064*d^29 + 2017866154384935261*d^28 + 
         4017519156609815472*d^27 - 646932236273292302370*d^26 + 
         6298798697799355942560*d^25 + 50563463138512620360075*d^24 - 
         605415122291976033499920*d^23 + 9173739160733459406640920*d^22 - 
         867689553019849170661939680*d^21 + 1257518319686780954500740240*d^20
          + 722149679096298021689411195520*d^19 - 
         8533523714288378843326304125440*d^18 - 
         366089879261824814758207043773440*d^17 + 
         11369763697247886350617258396015872*d^16 - 
         21443400458795669451533140394852352*d^15 - 55662748913999554099430843\
         11811180544*d^14 + 168356845402091603504227825892829536256*d^13 - 
         2906121996092959686320241831983719059456*d^12 + 353797930396313230040\
         04361848483410509824*d^11 - 325483821516047119178140090983459431448576
         *d^10 + 2332113061338228520623804746996237159890944*d^9 - 13188027525\
         113366130368970538859202358542336*d^8 + 59053213694903135002782614439\
         767947207507968*d^7 - 208458408396218918411804964961641452556779520*
         d^6 + 573162417251572399177568459857836703639142400*d^5 - 12014649510\
         17175213445683536025121657454592000*d^4 + 185335976666979110993194167\
         6527291329413120000*d^3 - 1981346134880789919038530166498360347852800\
         000*d^2 + 1309869196965031153283010062987575689216000000*d - 
         402965366102369913894935370396472442880000000,1125899906842624*d^12
          - 202661983231672320*d^11 + 16397606243255975936*d^10 - 
         787679574827099750400*d^9 + 24986385063817229893632*d^8 - 
         550614774185929306275840*d^7 + 8629173803375758093058048*d^6 - 
         96728573645153930387128320*d^5 + 768113625285264074440966144*d^4 - 
         4203998147827857569234288640*d^3 + 15012166487996600883049660416*d^2
          - 31303521115619017315823124480*d + 28717064759816527767876403200)
          + Gam(1,1)^3*M^-22*rat(1594323*d^33 - 181221381*d^32 + 6567370731*
         d^31 + 13826264301*d^30 - 7122914351820*d^29 + 159371233216938*d^28
          + 998048372529726*d^27 - 95474817207741798*d^26 + 
         1147699166119619787*d^25 + 16203849672609655263*d^24 - 
         527325704281840607793*d^23 - 4322934105873186041559*d^22 + 
         1432684149990959653004286*d^21 - 155956629658589008708796820*d^20 + 
         11266222983172926823175001336*d^19 - 574239566096047402187816049696*
         d^18 + 21905520105088058278439102190912*d^17 - 
         651351782006051408806691911540608*d^16 + 
         15483321597984550268188494522816000*d^15 - 29883123070969040664545333\
         9462093568*d^14 + 4727222153036155925200469690580395520*d^13 - 
         61631950824285733347703896987835130880*d^12 + 66402753510171413342305\
         7084460453378048*d^11 - 5913695499689708361040693661107546669056*d^10
          + 43442880511703011717809739357918484987904*d^9 - 262029880316893711\
         888863738344943091187712*d^8 + 12877254411132065026010791411501690933\
         73952*d^7 - 5096990277238378039203829802582724128538624*d^6 + 
         15976138410653519678047846864199776456409088*d^5 - 386794672504808741\
         46626205820813935693004800*d^4 + 696408779836636086855065328027529044\
         29568000*d^3 - 87639960139612217951772711751303424901120000*d^2 + 
         68678527448904425012201848630672065822720000*d - 25185335381398119618\
         433460649779527680000000,8796093022208*d^15 - 1688849860263936*d^14
          + 147528072168472576*d^13 - 7767583457307262976*d^12 + 
         275301526783949537280*d^11 - 6947730166742228402176*d^10 + 
         128799289605148993650688*d^9 - 1783650758520874507173888*d^8 + 
         18580439040133180493398016*d^7 - 145442189919019247180185600*d^6 + 
         847814273350886084623466496*d^5 - 3612505485043533531249115136*d^4 + 
         10890618679938003612305719296*d^3 - 21937128483588231413196914688*d^2
          + 26420659842740706570865213440*d - 14358532379908263883938201600) )
         ;


id,only intbn*dala^11*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-24*miBN*rat( - 
         4782969*d^37 + 631351908*d^36 - 30297451410*d^35 + 391191594336*d^34
          + 19547940036969*d^33 - 875503467104508*d^32 + 8525112320896380*d^31
          + 280522824190560864*d^30 - 9031873921347677895*d^29 + 
         49446108699138433260*d^28 + 1985018686807652653950*d^27 - 
         37801976765430765159360*d^26 + 98654550089945354154855*d^25 + 
         3231065357781767734626540*d^24 - 99858295583663649760758600*d^23 + 
         2638178645201200007544718080*d^22 + 20220338707338819335604175440*
         d^21 - 2037640158867107163807554402880*d^20 + 
         4411705172102225265994579814400*d^19 + 
         1231041815535551305304922836106240*d^18 - 
         19002390146529293149291560069639936*d^17 - 35237500581793634786389588\
         5178395648*d^16 + 16393646227299172449583277667648706560*d^15 - 
         227243674365109359376556626052491247616*d^14 - 2913289161306128982345\
         1465153347047424*d^13 + 57376256273970725426453706003881811591168*
         d^12 - 1204935008673823610768936985882880059310080*d^11 + 15111884201\
         690193345064156961223334729875456*d^10 - 1354232574014895664140345824\
         03326108410839040*d^9 + 917504779492069364934755403094732324900700160
         *d^8 - 4800430691807037476248966084781422161376051200*d^7 + 
         19484375673988126521260588446258089324183552000*d^6 - 609227209026904\
         76643000407282170488296570880000*d^5 + 144065577972150047455419749962\
         529021047603200000*d^4 - 24916427534595908160960743424868608481689600\
         0000*d^3 + 297171866956945091025585142909187496345600000000*d^2 - 
         218252690043292387573244991100749348864000000000*d + 7430737914896588\
         4714884488151489839104000000000,9007199254740992*d^12 - 
         1621295865853378560*d^11 + 131180849946047807488*d^10 - 
         6301436598616798003200*d^9 + 199891080510537839149056*d^8 - 
         4404918193487434450206720*d^7 + 69033390427006064744464384*d^6 - 
         773828589161231443097026560*d^5 + 6144909002282112595527729152*d^4 - 
         33631985182622860553874309120*d^3 + 120097331903972807064397283328*
         d^2 - 250428168924952138526584995840*d + 
         229736518078532222143011225600) + Gam(1,1)^3*M^-24*rat( - 4782969*
         d^34 + 588305187*d^33 - 24916611285*d^32 + 155969961885*d^31 + 
         21442182296760*d^30 - 686807594833374*d^29 + 1954362269597166*d^28 + 
         312051200422552218*d^27 - 6294210750920730681*d^26 - 
         13008819193977226761*d^25 + 2001916118275580530575*d^24 - 
         8583472304526467188071*d^23 - 2725726276302399778826934*d^22 + 
         285066436952023268663709108*d^21 - 18907194844490614662812822712*d^20
          + 858275273003689558026474632640*d^19 - 
         27584853919397352499944574619712*d^18 + 
         623643219753041434537733450285952*d^17 - 
         8809326995422557889171077993345024*d^16 + 
         19220871907151755133087373496402176*d^15 + 28242329913679053091220873\
         64536615424*d^14 - 90759670080238255207805762510894939136*d^13 + 
         1752448183325919801953464810893723224064*d^12 - 248792383469746344367\
         92794496509345292288*d^11 + 275247609147969703656586857899839179620352
         *d^10 - 2427282332787595187213065399089221177573376*d^9 + 17200411798\
         395181127185224321334854500745216*d^8 - 97914732758446313366886167330\
         600856122818560*d^7 + 444668911591967679855684861617036696758517760*
         d^6 - 1588901667833861326008875671659476395071897600*d^5 + 4366692000\
         337198761552574654407886076116992000*d^4 - 89004777703526318491872467\
         65174335555502080000*d^3 + 126611889387011794474723167147153042702336\
         00000*d^2 - 11208850967936252053953991325474453716992000000*d + 
         4644211196810367794680280509468114944000000000,70368744177664*d^15 - 
         13510798882111488*d^14 + 1180224577347780608*d^13 - 
         62140667658458103808*d^12 + 2202412214271596298240*d^11 - 
         55581841333937827217408*d^10 + 1030394316841191949205504*d^9 - 
         14269206068166996057391104*d^8 + 148643512321065443947184128*d^7 - 
         1163537519352153977441484800*d^6 + 6782514186807088676987731968*d^5
          - 28900043880348268249992921088*d^4 + 87124949439504028898445754368*
         d^3 - 175497027868705851305575317504*d^2 + 
         211365278741925652566921707520*d - 114868259039266111071505612800) );
#endprocedure

#procedure reduceBNBN
*
* reduce BN to BM and simpler integrals
*

* new method: 
* 1. Reduce indices of massless lines to zero:
*    n1 or n2 >0: "ordinary" triangle rule
*    n1 or n2 <0: "inverse" triangle rule
* 2. Use special recurrence relation for the resulting
*    BN(0,0,n3,n4,n5,n6) integrals.
*
*
* before each rec.rel.: use symmetry and expand (or use 'accun.prc')
*
************************************************************


************************************************************

* n1,n2 > 0

        #call symmetryBN
        if ( count(intbn,1) );
        if ( ( (count(x3,1)==count(x4,1)) || (count(x5,1)==count(x6,1)) )
        && (count(p1.p1,1) > count(p2.p2,1)) )
        multiply replace_(p1,p2,p2,p1);
        endif;

        #call Conv2exact        
        .sort
        

* n1>0:
        
        if ( count(intbn,1) );        
        if ( count(x3,1,x5,1) < count(x4,1,x6,1) );
        #call redBNn135(p1,p2,x3,x4,x5,x6)
        else;
        #call redBNn146(p1,p2,x3,x4,x5,x6)
        endif;
        endif;        
        .sort
        
        #call symmetryBN
        
        #call Conv2exact
        .sort

* n2>0:
        if ( count(intbn,1) );
        if ( count(x4,1,x5,1) < count(x3,1,x6,1) );
        #call redBNn245(p1,p2,x3,x4,x5,x6)
        else;
        #call redBNn236(p1,p2,x3,x4,x5,x6)
        endif;
        endif;        
        .sort

* n1,n2 < 0
        
        #call symmetryBN
        if ( count(intbn,1) );        
        if ( ( (count(x3,1)==count(x4,1)) || (count(x5,1)==count(x6,1)) )
        && (count(p1.p1,1) > count(p2.p2,1)) )
        multiply replace_(p1,p2,p2,p1);
        endif;        

        #call Conv2exact
        .sort

* n1<0:

        #call redBNn1p(p1,p2,x3,x4,x5,x6)
        .sort
        
        #call symmetryBN
        
        
************************************************************
        
* n2<0:

        #call redBNn2p(p1,p2,x3,x4,x5,x6)
        .sort

* now: BN(0,0,n3,n4,n5,n6) or BM-type integrals

        #call symmetryBN
        
        if ( count(intbn,1) );        
        if ( ( (count(x3,1)==count(x4,1)) || (count(x5,1)==count(x6,1)) )
        && (count(p1.p1,1) > count(p2.p2,1)) )
        multiply replace_(p1,p2,p2,p1);
        endif;        

        #call Conv2exact
        .sort
        

        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );
        multiply replace_(x3,x4,x4,x3);
        endif;
        endif;

        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );
        multiply replace_(x3,x5,x5,x3);
        endif;
        endif;

        #call Conv2exact
        .sort        


        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );
        multiply replace_(x3,x6,x6,x3);
        endif;
        endif;

        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort

        #call symmetryBN
        
        if ( count(intbn,1) );        
        if ( ( (count(x3,1)==count(x4,1)) || (count(x5,1)==count(x6,1)) )
        && (count(p1.p1,1) > count(p2.p2,1)) )
        multiply replace_(p1,p2,p2,p1);
        endif;
        #call Conv2exact        
        .sort
        
************************************************************

        #message Use table for BN


* Check if table is too small:
        
        if ( count(intbn,1) );   
        
        if ( count(dala,1) > 11) exit "ERROR: Table is too small for dala > 11 in BN reduction"; 

        #call BNdExact

        endif;        
        .sort

        if ( count(intbn,1) );                
        #call BNd0Exact
        
        endif;        
        .sort
        
#endprocedure




*
*
* Old routine to generate tables and for cases
* which are not tabulated 
*
* 
#procedure reduceBNnotab

        #message This is BN reduction without tables        
        
        #call symmetryBN
        
        #call ACCU{BN}
        
        #message 1
        
        #CALL redBNn5{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        
        #call ACCU{BN}
        

        #message 2
        
        #CALL redBNn34{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN
        
        #call ACCU{BN}
        
        #message 3
        .sort
        
        #CALL redBNn12{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN


                
        
        #call ACCU{BN}
        
        if( count(intbn,1));
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;        
        .sort
        #message 5
        
        
        #CALL redBNn6{p1|p2|x3|x4|x5|x6}

        
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN

        
                
        #call ACCU{BN}

        if( count(intbn,1));
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        .sort
        
                
        
        #call ACCU{BN}

        #message 6
        
        
        .sort
        
        #CALL redBNn12{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN
        
        
        #call ACCU{BN}
        
        #call ACCU(123)
        
        if( count(intbn,1));
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        .sort
        #message 7
        
        
        #CALL redBNn6{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN
        
        #message n6 done
        
        
        #call ACCU{BN}
        
        if( count(intbn,1));
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        .sort
        #message 8
        
        #call symBN{p1,p2,x3,x4,x5,x6}
        
        #call BNtoBM{p1,p2,p3,p4,p5,p6,x3,x4,x5,x6}
        .sort
        
        #call ACCU(topBN)
        
        if ( count(intbn,1) );                
        if ( (count(x3,1)==0) ) multiply intbm/intbn;
        endif;

        if (count(intbm,1)==1);
        
        id  p4.p4 = 1/x4 - M^2;
        id  p5.p5 = 1/x5 - M^2;
        id  p6.p6 = 1/x6 - M^2;
        
        id  p1.p5 = 1/2 * ( p3.p3 + 1/x4  - p2.p2 - 1/x6);
        
        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
*
* sort n6>=n5>=n4
*  
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) ) 
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        
        if ( (count(x4,1)<=0) || (count(x5,1)<=0) || (count(x6,1)<=0) );
*
* sort n6>n5>n4
*
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) )
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;

        endif;
        .sort

        #message - done
#endprocedure









#procedure topbn
*
* this is topbnbn
*

        #message this is topbnbn

        #message numerator
        #call nomBN

        #message do recursion

* Modified: applied only to intbn        
        #ifdef `REDBNTAB'
                #call reduceBNBN
                #else                
                #call reduceBNnotab
        #endif
        
        #call symBN{p1,p2,x3,x4,x5,x6}
        
        #call BNtoBM{p1,p2,p3,p4,p5,p6,x3,x4,x5,x6}
        .sort


        #call ACCU(topBN)

        if ( count(intbn,1) );                
        if ( (count(x3,1)==0));

        id  p4.p4 = 1/x4 - M^2;
        id  p5.p5 = 1/x5 - M^2;
        id  p6.p6 = 1/x6 - M^2;

        id  p1.p5 = 1/2 * ( p3.p3 + 1/x4  - p2.p2 - 1/x6);

        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
*
* sort n6>=n5>=n4
*  
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) ) 
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);

        if ( (count(x4,1)<=0) || (count(x5,1)<=0) || (count(x6,1)<=0) );
*
* sort n6>n5>n4
*
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) )
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;

        Multiply intbm/intbn;
        endif;
        endif;
        .sort

        #message - done


#endprocedure        


#procedure symBN1 (p1,p2,x3,x4,p5,x6)
        if( count(intbn1,1));
        if (count(`x3',1) > count(`x4',1))
        multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
        if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) );
        if (count(`x3',1) < count(`x4',1))
        multiply replace_(`x3',`x4',`x4',`x3');
        if (count(`x3',1) < count(`x6',1))
        multiply replace_(`x3',`x6',`x6',`x3');
        if (count(`x4',1) < count(`x6',1))
        multiply replace_(`x4',`x6',`x6',`x4');
        endif;
        endif;
#endprocedure





#procedure redBN1n12to0 (p1,p2,x3,x4,p5,x6)

* reduces n1 and n2 to zero, if:
* a. n1<=0, n2>0
* b. n1<=0, n2<0

* a.
*
* reduces n2 to zero if n1<=0 (n3,n4>=1, n5,n6=1)
* result: BN1(n1,n2,n3,n4,n5,n6) with n1<=0, n2=0
*         or BM's.

* sort: n2>=n1

        if( count(intbn1,1) );        
        if ( ( count(`p2'.`p2',1) > count(`p1'.`p1',1) ) &&
        (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)==1) &&
        (count(`p5'.`p5',1)<=-1) )
        multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
        endif;
        
        #call redBN1n6 (`p1',`p2',`x3',`x4',`p5',`x6')
        
        #do i=1,10
                
                if( count(intbn1,1) );        
                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)<=-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)==1) &&
                (count(`p5'.`p5',1)<=-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/( - 2*n5 - n1 + 2*n2 + n4 )
                
                *(
                + `x4'*`x6'^-1 * (  - n4 )
                
                + `p1'.`p1'^-1*`p5'.`p5' * (  - 2*n1 )
                
                + `p1'.`p1'^-1*`x3'^-1 * ( n1 )
                
                + `p1'.`p1'^-1*`x4'^-1 * ( 2*n1 )
                
                + `p1'.`p1'^-1*`x6'^-1 * (  - n1 )
                
                - `p1'.`p1'^-1 * ( 2*n1*M^2 )
                
                + `p2'.`p2'*`x3' * ( 2*n3 )
                
                + `p2'.`p2'*`x4' * ( n4 )
                
                + `p5'.`p5'*`x3' * (  - 2*n3 )
                
                );
                
                redefine i "0";
                
                endif;
                endif;
                #call ACCU(BN1n12to0)
                
        #enddo
        

* b.
*
* reduces n2 to zero if n1<=0,n2<0 (n3,n4,n6>=1, n5=1)
* result: BN1(n1,n2,n3,n4,n5,n6) with n1<=0, n2=0
*         or BM's.

        #do i=1,1

* sort: n2>=n1
                if( count(intbn1,1) );
                if ( ( count(`p2'.`p2',1) < count(`p1'.`p1',1) ) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<=-1) )
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');

                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)>0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<=-1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)*deno( + 3/2*4 - n1 - n2 - n3 - n4 - n5 - n6,-3)
                *(
                - `p1'.`p1'^-1*`p2'.`p2'^-1*`x3'^-1 * ( n1*M^2 )

                - `p1'.`p1'^-1*`p2'.`p2'^-1*`x6'^-1 * (  - n1*M^2 )

                - `p2'.`p2'^-1*`p5'.`p5'*`x3' * (  - n3*M^2 )

                + `p2'.`p2'^-1*`x3' * ( 3*n3*M^4 )

                - `p2'.`p2'^-1*`x4'^-1*`x6' * (  - n6*M^2 )

                - `p2'.`p2'^-1*`x4'*`x6'^-1 * (  - n4*M^2 )

                - `p2'.`p2'^-1 *(- 3*num(d)*M^2 + n1*M^2 + 4*n2*M^2 + 3*n3*M^2 + n4*M^2 + 2*
                n5*M^2 + n6*M^2 + 4*M^2 )

                );

                redefine i "0";

                endif;
                endif;
                
                .sort
        #enddo
        
        #call ACCU(BN1n12to0)
        
#endprocedure

#procedure redBN1n12to1 (p1,p2,x3,x4,p5,x6)

        #do i=1,1
                
* sort: n2>=n1
                if( count(intbn1,1) );
                if ( ( count(`p2'.`p2',1) > count(`p1'.`p1',1) ) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<=-1) )
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                
                if ( (count(`p1'.`p1',1)<=-1) && (count(`p2'.`p2',1)<-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<=-1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)
                *(
                - `x3'^-1 * ( M^-2 )

                - `p1'.`p1'*`p2'.`p2'*`x4' * ( 1/( - 1 + n2)*n4*M^-2 )

                - `p2'.`p2'*`p5'.`p5'*`x4' * (  - 1/( - 1 + n2)*n4*M^-2 )

                + `p2'.`p2'*`x4' * (  - 1/( - 1 + n2)*n4 )

                - `p2'.`p2' * (  - M^-2 + 1/( - 1 + n2)*nom(4,-2)*M^-2 
                - 1/( - 1 + n2)*n4*M^-2 - 2/( - 1 + n2)*n5*M^-2 )

                - `p5'.`p5' * (  - M^-2 )

                );

                redefine i "0";

                endif;
                endif;
                #call ACCU(BN1)
                
        #enddo
        
#endprocedure

#procedure redBN1n32 (p1,p2,x3,x4,p5,x6)

* reduce n5 to 1 (if n1=n2=0)
* (take the procedure from redBN1n5.prc)

        #do i=1,1
                
                if( count(intbn1,1) );        
                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2?
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)*M^2*deno(4*n5 - 4 - n1*n5 + n1 + 2*n5 - 2*n5^2,-2*n5+2)
                *(
                + `p1'.`p1'^-1*`p5'.`p5'*`x3'*`x6'^-1 * (  - n1*n3*M^-2 )

                + `p1'.`p1'^-1*`p5'.`p5' * ( n1*n3*M^-2 - n1*n5*M^-2 + n1*M^-2 )

                + `p1'.`p1'^-1*`x4'^-1 * ( n1*n5*M^-2 - n1*M^-2 )

                - `p1'.`p1'^-1 * ( n1*n5 - n1 )

                + `p5'.`p5'*`x3' *(-nom(4,-2)*n3*M^-2+n1*n3*M^-2+2*n3*M^-2+2*n3^2*M^-2 )

                - `p5'.`p5'*`x3'^2 * ( 2*n3 + 2*n3^2 )
                );

                redefine i "0";

                endif;
                endif;
                #call ACCU(BN1n32_0)
                
        #enddo
        
* reduce n3, n4 and n6 to 2, take care if n6=n3-1!!!
* n1=n2=0!
        
* The condition n6=n3-1 is asked for with the following trick; 
* it only works if n5=1!
        
        #do i=1,1
                
                if( count(intbn1,1) );        
                if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) );
                if (count(`x3',1) < count(`x4',1))
                multiply replace_(`x3',`x4',`x4',`x3');
                if (count(`x3',1) < count(`x6',1))
                multiply replace_(`x3',`x6',`x6',`x3');
                if (count(`x4',1) < count(`x6',1))
                multiply replace_(`x4',`x6',`x6',`x4');
                endif;
                
                if ( ( count(`x3',1) > count(`x6',1,`p5'.`p5',-1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(- 6*n3*M^2 + 2*n3^2*M^2 + 4*M^2)
                *(
                - `x3'^-2*`x6' * (  - nom(4,-2)*n6 + n3*n6 + 2*n5*n6 - 2*n6 )
                
                - `x3'^-1*`x4'^-1*`x6' * (  - n3*n6 + 2*n6 )
                
                + `x3'^-1*`x6' * ( n3*n6*M^2 - 2*n6*M^2 )
                
                - `x3'^-1*(4-2*nom(4,-2)*n3+4*nom(4,-2)+2*n3*n5+n3*n6-6*n3+2*n3^2 - 4*n5
                - 2*n6 )
                
                - `p5'.`p5'*`x3'^-1*`x6' * ( n3*n6 - 2*n6 )
                );

                redefine i "0";
                
                endif;
                
                if ( ( count(`x3',1) == count(`x6',1,`p5'.`p5',-1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/( (- 6*n3*M^2 + 2*n3^2*M^2 + 4*M^2) + M^2*n6*(n3-2) )
                *(
                - `x3'^-2*`x6' * (  - nom(4,-2)*n6 + n3*n6 + 2*n5*n6 - 2*n6 )
                
                - `x3'^-1*`x4'^-1*`x6' * (  - n3*n6 + 2*n6 )
                
                - `x3'^-1 * ( 4 - 2*nom(4,-2)*n3 + 4*nom(4,-2) 
                + 2*n3*n5 + n3*n6 - 6*n3 + 2*n3^2 - 4*n5
                - 2*n6 )
                
                - `p5'.`p5'*`x3'^-1*`x6' * ( n3*n6 - 2*n6 )
                );
                
                redefine i "0";

                endif;
                
                if ( ( count(`x3',1) == count(`x6',1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(M^2 * ( 3 - 9/2*n3 + 3/2*n3^2 ) )
                *(
                - 1/`x3'/`x3'*`x6' * (  - nom(4,-2)*n3 + n3^2 )
                
                - 1/`x3'/`x4'*`x6' * ( 2*n3 - n3^2 )
                
                - 1/`x3'*`p5'.`p5'*`x6' * (  - 2*n3 + n3^2 )

                - 1/`x3' * ( 1 - 3/2*nom(4,-2)*n3 + 3*nom(4,-2) - 11/2*n3 + 5/2*n3^2 )
                
                - 1/`x4' * ( 1 - 3/2*n3 + 1/2*n3^2 )
                
                - `p5'.`p5' * (  - 1 + 3/2*n3 - 1/2*n3^2 )
                
                - 1/`x6' * ( 1 + nom(4,-2)*n3 - 2*nom(4,-2) + 5/2*n3 - 3/2*n3^2 )
                
                );
                                
                redefine i "0";
                
                endif;
                endif;
                
                #call ACCU(BN1n32_0)
                
        #enddo

#endprocedure

#procedure redBN1n32s (p1,p2,x3,x4,p5,x6)

* reduce n5 to 1 (if n1=n2=0)
* (take the procedure from redBN1n5.prc)
        
        #do ii=1,1
                
                if( count(intbn1,1) );        
                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2?
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)*M^2*deno(4*n5 - 4 - n1*n5 + n1 + 2*n5 - 2*n5^2,-2*n5+2)
                *(
                + `p1'.`p1'^-1*`p5'.`p5'*`x3'*`x6'^-1 * (  - n1*n3*M^-2 )
                
                + `p1'.`p1'^-1*`p5'.`p5' * ( n1*n3*M^-2 - n1*n5*M^-2 + n1*M^-2 )
                
                + `p1'.`p1'^-1*`x4'^-1 * ( n1*n5*M^-2 - n1*M^-2 )
                
                - `p1'.`p1'^-1 * ( n1*n5 - n1 )
                
                + `p5'.`p5'*`x3' *(-nom(4,-2)*n3*M^-2+n1*n3*M^-2+2*n3*M^-2+2*n3^2*M^-2 )
                
                - `p5'.`p5'*`x3'^2 * ( 2*n3 + 2*n3^2 )
                );
                
                redefine j "0";
                redefine ii "0";
                
                endif;
                endif;        

                #call Conv2exact
                .sort
                
        #enddo
        
        #call ACCU(BN1n32_1)
        
* reduce n3, n4 and n6 to 2, take care if n6=n3-1!!!
* n1=n2=0!
        
* The condition n6=n3-1 is asked for with the following trick; 
* it only works if n5=1!
                
        #do i=1,1

                if( count(intbn1,1) );
                if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) );
                if (count(`x3',1) < count(`x4',1))
                multiply replace_(`x3',`x4',`x4',`x3');
                if (count(`x3',1) < count(`x6',1))
                multiply replace_(`x3',`x6',`x6',`x3');
                if (count(`x4',1) < count(`x6',1))
                multiply replace_(`x4',`x6',`x6',`x4');
                endif;
                endif;
                #call ACCU(BN1n32_2)
                
                if( count(intbn1,1) );        
                if ( ( count(`x3',1) > count(`x6',1,`p5'.`p5',-1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(- 6*n3*M^2 + 2*n3^2*M^2 + 4*M^2)
                *(
                - `x3'^-2*`x6' * (  - nom(4,-2)*n6 + n3*n6 + 2*n5*n6 - 2*n6 )
                
                - `x3'^-1*`x4'^-1*`x6' * (  - n3*n6 + 2*n6 )
                
                + `x3'^-1*`x6' * ( n3*n6*M^2 - 2*n6*M^2 )
                
                - `x3'^-1 * ( 4 - 2*nom(4,-2)*n3 + 4*nom(4,-2) + 2*n3*n5 + n3*n6 
                 - 6*n3 + 2*n3^2 - 4*n5
                - 2*n6 )
                
                - `p5'.`p5'*`x3'^-1*`x6' * ( n3*n6 - 2*n6 )
                );
                
                redefine i "0";
                redefine j "0";
                
                endif;
                endif;
                #call ACCU(BN1n32_3)
                
                if( count(intbn1,1) );        
                if ( ( count(`x3',1) == count(`x6',1,`p5'.`p5',-1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/( (- 6*n3*M^2 + 2*n3^2*M^2 + 4*M^2) + M^2*n6*(n3-2) )
                *(
                - `x3'^-2*`x6' * (  - nom(4,-2)*n6 + n3*n6 + 2*n5*n6 - 2*n6 )
                
                - `x3'^-1*`x4'^-1*`x6' * (  - n3*n6 + 2*n6 )
                
                - `x3'^-1 * ( 4 - 2*nom(4,-2)*n3 + 4*nom(4,-2) 
                + 2*n3*n5 + n3*n6 - 6*n3 + 2*n3^2 - 4*n5
                - 2*n6 )
                
                - `p5'.`p5'*`x3'^-1*`x6' * ( n3*n6 - 2*n6 )
                );
                
                redefine i "0";
                redefine j "0";
                
                endif;
                endif;
                
                #call ACCU(BN1n32_4)
                
                if( count(intbn1,1) );        
                if ( ( count(`x3',1) == count(`x6',1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(M^2 * ( 3 - 9/2*n3 + 3/2*n3^2 ) )
                *(
                - 1/`x3'/`x3'*`x6' * (  - nom(4,-2)*n3 + n3^2 )
                
                - 1/`x3'/`x4'*`x6' * ( 2*n3 - n3^2 )
                
                - 1/`x3'*`p5'.`p5'*`x6' * (  - 2*n3 + n3^2 )
                
                - 1/`x3' * ( 1 - 3/2*nom(4,-2)*n3 + 3*nom(4,-2) - 11/2*n3 + 5/2*n3^2 )
                
                - 1/`x4' * ( 1 - 3/2*n3 + 1/2*n3^2 )
                
                - `p5'.`p5' * (  - 1 + 3/2*n3 - 1/2*n3^2 )
                
                - 1/`x6' * ( 1 + nom(4,-2)*n3 - 2*nom(4,-2) + 5/2*n3 - 3/2*n3^2 )
                
                );
                
                redefine i "0";
                redefine j "0";
                
                endif;
                endif;
                
                #call ACCU(BN1n32_5)
                
* added Jul. '98
                
                if( count(intbn1,1) );        
                if ( (count(`p1',1)==0) && (count(`p2',1)==0) );
                if ((count(`x3',1)<=0) || (count(`x4',1)<=0) || (count(`x6',1)<=0)) discard;
                if (count(`x3',1) > count(`x4',1)) multiply replace_(`x3',`x4',`x4',`x3');
                if (count(`x3',1) > count(`x6',1)) multiply replace_(`x3',`x6',`x6',`x3');
                if (count(`x4',1) > count(`x6',1)) multiply replace_(`x4',`x6',`x6',`x4');
                endif;
                endif;
                
                #call ACCU(BN1n32_6)

#enddo

#endprocedure



#procedure redBN1n34 (p1,p2,x3,x4,p5,x6)
        
* Use this procedure to reduce n3 and n4 to 1.
* n1=n2=1, n5=n6=1;
        
        #do i=1,1
                
                if( count(intbn1,1) );        
                if ( ( count(`x3',1) < count(`x4',1) ) &&
                (count(`x3',1)>=1) && (count(`x4',1)>1) && (count(`x6',1)==1) &&
                (count(`p5'.`p5',1)==-1) ) 
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');

                if ( (count(`p1'.`p1',1)==-1) && (count(`p2'.`p2',1)==-1) &&
                (count(`x3',1)>1) && (count(`x4',1)>=1) && (count(`x6',1)==1) &&
                (count(`p5'.`p5',1)==-1) );


                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(n3-1)
                *(
                + `x3'^-1*`x4'^-1 * (  - 3/2*nom(4,-2)*M^-4 
                + n1*M^-4 + n2*M^-4 + n3*M^-4 + n4*M^-4
                + n5*M^-4 + n6*M^-4 - 2*M^-4 )

                - `x3'^-1 * (  - nom(4,-2)*M^-2 + 2*n2*M^-2 
                + n3*M^-2 + n4*M^-2 + n6*M^-2 - 2*M^-2 )

                - `x4'^-1 * ( n3*M^-2 - M^-2 )
                
                - `p2'.`p2'*`x3'^-1*`x4' * (  - n4*M^-2 )
                
                + `p2'.`p2'*`x3'^-1 * ( 3/2*nom(4,-2)*M^-4 
                - n1*M^-4 - n2*M^-4 - n3*M^-4 - n4*M^-4 -
                n5*M^-4 - n6*M^-4 + 2*M^-4 )
                
                - `p5'.`p5' * (  - n3*M^-2 + M^-2 )
                );

                redefine i "0";

                endif;
                endif;
        
                #call ACCU(BN1)
                
        #enddo

#endprocedure


#procedure redBN1n5 (p1,p2,x3,x4,p5,x6)
        
        #do i=1,1
                if( count(intbn1,1) );
                if ( (count(`p1'.`p1',1)<=-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(1 - n5)
                *(
                - `x4'^-1 * (  - n5*M^-2 + M^-2 )

                - `p1'.`p1' * ( n5*M^-2 - M^-2 )

                - `p2'.`p2'^-1*`p5'.`p5'*`x4'^-1 * (  - n2*M^-2 )

                - `p2'.`p2'^-1*`p5'.`p5'*`x6'^-1 * ( n2*M^-2 )

                + `p5'.`p5'*`x4' * (  - 2*n4 )

                - `p5'.`p5' * ( nom(4,-2)*M^-2 - n2*M^-2 - 2*n4*M^-2 - n5*M^-2 + M^-2 )
                );

                redefine i "0";

                endif;

                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)<=-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(1 - n5)
                *(
                - `x3'^-1 * (  - n5*M^-2 + M^-2 )

                - `p1'.`p1'^-1*`p5'.`p5'*`x3'^-1 * (  - n1*M^-2 )

                - `p1'.`p1'^-1*`p5'.`p5'*`x6'^-1 * ( n1*M^-2 )

                - `p2'.`p2' * ( n5*M^-2 - M^-2 )

                + `p5'.`p5'*`x3' * (  - 2*n3 )

                - `p5'.`p5' * ( nom(4,-2)*M^-2 - n1*M^-2 - 2*n3*M^-2 - n5*M^-2 + M^-2 )
                );

                redefine i "0";

                endif;
* topBN1        
                endif;

                #call ACCU(BN1)
                
        #enddo
        
        
        #do i=1,1
                
                if( count(intbn1,1) );        
                if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)>=0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                if (count(`p1'.`p1',1) < count(`p2'.`p2',1))
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                endif;
                
                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)*M^2*deno(4*n5 - 4 - n1*n5 + n1 + 2*n5 - 2*n5^2,-2*n5+2)
                *(
                + `p1'.`p1'^-1*`p5'.`p5'*`x3'*`x6'^-1 * (  - n1*n3*M^-2 )
                
                + `p1'.`p1'^-1*`p5'.`p5' * ( n1*n3*M^-2 - n1*n5*M^-2 + n1*M^-2 )
                
                + `p1'.`p1'^-1*`x4'^-1 * ( n1*n5*M^-2 - n1*M^-2 )
                
                - `p1'.`p1'^-1 * ( n1*n5 - n1 )
                
                + `p5'.`p5'*`x3' *(-nom(4,-2)*n3*M^-2 
                + n1*n3*M^-2 + 2*n3*M^-2 + 2*n3^2*M^-2 )
                
                - `p5'.`p5'*`x3'^2 * ( 2*n3 + 2*n3^2 )
                );
                
                redefine i "0";
                
                endif;
                endif;
                
                #call Conv2exact
                .sort
        #enddo
        
        if( count(intbn1,1) );
        repeat;
                
                if ( (count(`p1'.`p1',1)<=-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(1 - n5)
                *(
                - `x4'^-1 * (  - n5*M^-2 + M^-2 )
                
                - `p1'.`p1' * ( n5*M^-2 - M^-2 )
                
                - `p2'.`p2'^-1*`p5'.`p5'*`x4'^-1 * (  - n2*M^-2 )
                
                - `p2'.`p2'^-1*`p5'.`p5'*`x6'^-1 * ( n2*M^-2 )
                
                + `p5'.`p5'*`x4' * (  - 2*n4 )
                
                - `p5'.`p5' * ( nom(4,-2)*M^-2 - n2*M^-2 - 2*n4*M^-2 - n5*M^-2 + M^-2 )
                );
                endif;
                
                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)<=-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(1 - n5)
                *(
                - `x3'^-1 * (  - n5*M^-2 + M^-2 )
                
                - `p1'.`p1'^-1*`p5'.`p5'*`x3'^-1 * (  - n1*M^-2 )
                
                - `p1'.`p1'^-1*`p5'.`p5'*`x6'^-1 * ( n1*M^-2 )
                
                - `p2'.`p2' * ( n5*M^-2 - M^-2 )
                
                + `p5'.`p5'*`x3' * (  - 2*n3 )
                
                - `p5'.`p5' * ( nom(4,-2)*M^-2 - n1*M^-2 - 2*n3*M^-2 - n5*M^-2 + M^-2 )
                );
                endif;
                
        endrepeat;
        endif;
        

        #call Conv2exact
        
        if( count(intbn1,1) );
        repeat;
                
                if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)>=0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                if (count(`p1'.`p1',1) < count(`p2'.`p2',1))
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                endif;
                
                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)*M^2*deno(4*n5 - 4 - n1*n5 + n1 + 2*n5 - 2*n5^2,-2*n5+2)
                *(
                + `p1'.`p1'^-1*`p5'.`p5'*`x3'*`x6'^-1 * (  - n1*n3*M^-2 )
                
                + `p1'.`p1'^-1*`p5'.`p5' * ( n1*n3*M^-2 - n1*n5*M^-2 + n1*M^-2 )
                
                + `p1'.`p1'^-1*`x4'^-1 * ( n1*n5*M^-2 - n1*M^-2 )
                
                - `p1'.`p1'^-1 * ( n1*n5 - n1 )
                
                + `p5'.`p5'*`x3' *(-nom(4,-2)*n3*M^-2 
                + n1*n3*M^-2 + 2*n3*M^-2 + 2*n3^2*M^-2 )
                
                - `p5'.`p5'*`x3'^2 * ( 2*n3 + 2*n3^2 )
                );
                endif;
                
        endrepeat;
        endif;
        
#endprocedure

#procedure redBN1n6 (p1,p2,x3,x4,p5,x6)
        
        if( count(intbn1,1) );        
        repeat;
                if ( (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>1) &&
                (count(`p5'.`p5',1)<=-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^(n6-1)
                *1/(n6-1)
                *(
                - n3*`x3' - n4*`x4'
                - 1/M^2*nom(6-n1-n2-n3-n4-n5-(n6-1),-3)
                );
                endif;
        endrepeat;
        endif;
        
#endprocedure





#procedure reduceBN1
        
        #call redBN1n5(p1,p2,x3,x4,p5,x6)
        .sort

        
        #call redBN1n12to1(p1,p2,x3,x4,p5,x6)
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
        #call redBN1n6 (p1,p2,x3,x4,p5,x6)
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
        #call redBN1n34 (p1,p2,x3,x4,p5,x6)
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
        #call redBN1n12to0 (p1,p2,x3,x4,p5,x6)
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
* set massless tadpoles to zero:
        
        if( count(intbn1,1) );        
        if ( (count(x3,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x4,1)<=0) &&
        ( (count(p1.p1,1)>=0) || (count(p2.p2,1)>=0) || (count(p5.p5,1)>=0) )
        ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) && (count(p5.p5,1)<0) &&
        ( (count(x3,1)<=0)    || (count(x4,1)<=0) ||
        (count(x6,1)<=0) )
        ) discard;
        endif;        
        .sort
        

        #do j=1,1
                #call redBN1n32s (p1,p2,x3,x4,p5,x6)
        #enddo
        
        #call redBN1n32 (p1,p2,x3,x4,p5,x6)
        .sort
        
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
* now treat the integrals BN1(0,0,2,2,1,1), BN1(0,0,2,2,1,2)
* and BN1(0,0,2,1,1,1) separate:
        
* BN1(0,0,2,2,1,1) = -1/3/M^2*BN1(0,0,2,1,0,2);
        
        if( count(intbn1,1) );
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) && (count(x3,1)==2) && 
        (count(x4,1)==2) && (count(p5.p5,1)==-1) && (count(x6,1)==1) )
        id x3^2*x4^2/p5.p5*x6 = -1/3/M^2 * x3^2*x4*x6^2;
        endif;
        .sort
        
*   BN(0,0,2,2,1,2) =
*       + 7/6*BN1(0,0,2,1,0,2)*n*M^-4 - 38/9*BN1(0,0,2,1,0,2)*M^-4 
*       + 19/3*BN1(0,0,2,1,1,1)*n*M^-4 - BN1(0,0,2,1,1,1)*n^2*M^-4 
*       - 10*BN1(0,0,2,1,1,1)*M^-4 + 4/3*BN1(0,0,3,1,0,2)*M^-2;
        
        if( count(intbn1,1) );
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) && (count(x3,1)==2) && 
        (count(x4,1)==2) && (count(p5.p5,1)==-1) && (count(x6,1)==2) )
        id x3^2*x4^2/p5.p5*x6^2 =
        + 7/6* x3^2*x4*x6^2 *num(d)*M^-4      - 38/9* x3^2*x4*x6^2 *M^-4 
        + 19/3* x3^2*x4/p5.p5*x6 *num(d)*M^-4 - x3^2*x4/p5.p5*x6 *num(d^2)*M^-4 
        - 10* x3^2*x4/p5.p5*x6 *M^-4     + 4/3* x3^3*x4*x6^2 *M^-2
        ;
        endif;
        .sort
        
*   BN(0,0,2,1,1,1) =
*      -1/3/M^2*(3/2*n-4)*BN1(0,0,1,1,1,1);
        
        if( count(intbn1,1) );
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) && (count(x3,1)==2) && 
        (count(x4,1)==1) && (count(p5.p5,1)==-1) && (count(x6,1)==1) )
        id x3^2*x4/p5.p5*x6 =
        -1/3/M^2*(3/2*nom(4,-2)-4)* x3*x4/p5.p5*x6
        ;
        endif;
        
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
#endprocedure        





#procedure topbn1
*
* this is topbn1
*
        #-
        #message this is topbn1
        
        #message numerator
        
        if( count(intbn1,1) );        
        id  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
        id  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
        id  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
        endif;        
        #call ACCU(BN1 1)
        
        if( count(intbn1,1) );                
        id  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );
        id  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
        endif;        
        #call ACCU(BN1 2)
        
        if( count(intbn1,1) );                
        id  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
        id  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
        endif;        
        #call ACCU(BN1 3)
        
        if( count(intbn1,1) );                
        id  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
        id  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
        endif;        
        #call ACCU(BN1 4)
        
        if( count(intbn1,1) );        
        id  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
        id  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
        endif;        
        #call ACCU(BN1 5)
        
        if( count(intbn1,1) );        
        id  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
        id  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN1 6)
        
        if( count(intbn1,1) );        
        id  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
        id  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN1 7)
        
*
* Warning!
*
        if( count(intbn1,1) );
        id  1/x5 = M^2 + p5.p5;
        
        id  p3.p3 = 1/x3 - M^2;
        endif;
        
        #call ACCU(BN1 6)
        if( count(intbn1,1) );        
        id  p4.p4 = 1/x4 - M^2;
        id  p6.p6 = 1/x6 - M^2;
        endif;
        
        #call ACCU(BN1 7)
        
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
        if( count(intbn1,1) );        
        id x3^n3?neg_=(p3.p3+M^2)^-n3;
        id x4^n4?neg_=(p4.p4+M^2)^-n4;
        id x6^n6?neg_=(p6.p6+M^2)^-n6;
        
        
        if ( (count(x3,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x4,1)<=0) && 
        ( (count(p1.p1,1)>=0) || (count(p2.p2,1)>=0) || (count(p5.p5,1)>=0) )
        ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) && (count(p5.p5,1)<0) &&
        ( (count(x3,1)<=0)    || (count(x4,1)<=0) || 
        (count(x6,1)<=0) )
        ) discard;
        endif;        
        .sort
        
        #message do recursion
        
        #call reduceBN1
        
        if( count(intbn1,1) );        
        if ( (count(x3,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x4,1)<=0) && 
        ( (count(p1.p1,1)>=0) || (count(p2.p2,1)>=0) || (count(p5.p5,1)>=0) )
        ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) && (count(p5.p5,1)<0) &&
        ( (count(x3,1)<=0)    || (count(x4,1)<=0) || 
        (count(x6,1)<=0) )
        ) discard;
        
        if (count(x6,1)==0);
        id p6=p1+p3;
        id p1=-p1;
        multiply replace_(p5,p4,p4,p6,p3,p5,x3,x5,x4,x6);
        multiply intbm1/intbn1;
        elseif (count(x3,1)==0);
        id p3=p5-p2;
        id p1=-p1;
        multiply replace_(p5,p4,p4,p6,p6,p5,p2,p3,x4,x6,x6,x5);
        multiply intbm1/intbn1;
        elseif (count(x4,1)==0);
        id p4=p6-p2;
        multiply replace_(p5,p4,p3,p5,p1,p3,x3,x5);
        multiply intbm1/intbn1;
        elseif (count(p5.p5,1)>=0);
        id p5=p4-p1;
        multiply replace_(p3,p5,p1,p3,p2,p1,x3,x5);
        multiply intbm/intbn1;
        elseif ( (count(x3,1)!=0) && (count(x4,1)!=0) && (count(x6,1)!=0) );
        id 1/p1.p1^n1? /p2.p2^n2? * x3^n3? * x4^n4? /p5.p5^n5? * x6^n6? =
        BN1(n1,n2,n3,n4,n5,n6);        
        endif;
* topBN1        
        endif;        
        .sort

* insert the expansion for BN1(0,0,1,1,1,1) and BN1(1,1,1,1,1,1).

        if(count(intbn1,1));        
        id BN1(0,0,1,1,1,1) = +  M^4*miBN1*int0/intbn1;
        
        id BN1(1,1,1,1,1,1) = + miD3*int0/intbn1;
        endif;
        .sort
        
        #call ACCU();
        
        #message - done
        
#endprocedure        


#procedure topbn2
*
* this is topbn2
*
        #-
        #message this is topbn2

        #message numerator

        if( count(intbn2,1));        
        id  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
        id  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
        endif;        
        #call ACCU(BN2 1)

        if( count(intbn2,1));        
        id  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
        id  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );
        endif;        
        #call ACCU(BN2 2)

        if( count(intbn2,1));        
        id  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
        id  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
        endif;        
        #call ACCU(BN2 3)

        if( count(intbn2,1));        
        id  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
        id  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
        endif;        
        #call ACCU(BN2 4)

        if( count(intbn2,1));        
        id  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
        id  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN2 5)

        if( count(intbn2,1));        
        id  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
        id  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN2 6)

        if( count(intbn2,1));        
        id  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
        id  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
        id  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN2 7)

*
* Warning! 
*
        if( count(intbn2,1));        
        id  1/x3 = M^2 + p3.p3;
        id  1/x5 = M^2 + p5.p5;
        endif;        
        #call ACCU(BN2 8)

        if( count(intbn2,1));        
        id  p6.p6 = 1/x6 - M^2;
        id  p4.p4 = 1/x4 - M^2;
        endif;        
        #call ACCU(BN2 7)

        #message do recursion

        if( count(intbn2,1));        
        id p5=-p5;
        multiply replace_(p6,p5,p4,p6,p1,p4,p3,p2,p5,p1,p2,p3,x6,x5,x4,x6);
        multiply intbm1/intbn2;
        endif;        
        .sort

        #message - done
        
#endprocedure        



************************************************************

#procedure symBN3 (p1,p2,p3,p4,p5,x6)
        if( count(intbn3,1));
        if (count(`p1'.`p1',1) >= count(`p2'.`p2',1))
        multiply replace_(`p1',`p2',`p2',`p1',`p4',`p3',`p3',`p4');
        endif;
#endprocedure




*--#[ ntriangle :
*
#procedure ntriangle(p,pa,pb,p1,p2)
*
*   Routine solves the triangle recursion
*                     p
*      n1 -------------------------- n2
*          p1    \    n0    /    p2
*                 \        /
*               n3 \      / n4
*                pa \    / pb
*                    \  /
*                     \/
*
*   n3,n4,n0,n1,n2 are the powers of the denominators
*   eppa and eppb are 1/pa.pa^ep and 1/pb.pb^ep
*   p is the momentum in n0. We need it here to determine the extra
*   momenta in the numerator.
*
id	`p' = xpower*`p';
if ( count(`p1'.`p1',1,`p2'.`p2',1,`p'.`p',1) >= -4 );
*
*	Here we can just use the recursion
*
	repeat;
		id	xpower^n8?/`p1'.`p1'^n1?pos_/`p2'.`p2'^n2?pos_/`p'.`p'^n?pos_/`pa'.`pa'^n3?/`pb'.`pb'^n4?
			*ep`pa'^x1?*ep`pb'^x2? = xpower^n8/`p1'.`p1'^n1/`p2'.`p2'^n2/`p'.`p'^n/`pa'.`pa'^n3/`pb'.`pb'^n4
				*ep`pa'^x1*ep`pb'^x2*(
					+num(n3+x1*(2-d/2))*(`p'.`p'/`pa'.`pa'-`p1'.`p1'/`pa'.`pa')
					+num(n4+x2*(2-d/2))*(`p'.`p'/`pb'.`pb'-`p2'.`p2'/`pb'.`pb')
				)*den(d+n8-2*n-n3-n4-x1*(2-d/2)-x2*(2-d/2));
	endrepeat;
	id	num(x?number_) = x;
	id	den(x?number_) = 1/x;
	id	num(x?)*den(x?) = 1;
	id	num(x?) = rat(x,1);
	id	den(x?) = rat(1,x);
else;
*
*	Here we have to use the general formula
*
	id	xpower^n8?/`p1'.`p1'^n1?pos_/`p2'.`p2'^n2?pos_/`p'.`p'^n?pos_/`pa'.`pa'^n3?/`pb'.`pb'^n4?
			*ep`pa'^x1?*ep`pb'^x2? = ftriangle(n8,n,n3,x1,n4,x2,n1,n2);
	id	ftriangle(n?,n0?,n3?,x1?,n4?,x2?,n1?,n2?) =
		+sum_(isum1,0,n2-1,sum_(isum2,0,n0-1,sum_(isum3,0,n0-1-isum2,
			ftriangle(n,n0-isum2-isum3,n3+n1+isum2,x1,n4+isum1+isum3,x2,0,n2-isum1)*sign_(n1+isum1)
				*fac_(n1+isum1+isum2+isum3-1)
				*Pochhammer(n1+isum2,n3+x1*(2-d/2))
				*Pochhammer(isum3+isum1,n4+x2*(2-d/2))
				*Pochhammer(-n1-isum1-isum2-isum3,isum2+isum3+5-2*(2-d/2)+n-2*n0-n3-n4-x1*(2-d/2)-x2*(2-d/2))
				*invfac_(isum1)*invfac_(isum2)*invfac_(isum3)*invfac_(n1-1)
		)))
		+sum_(isum1,0,n1-1,sum_(isum2,0,n2-1,sum_(isum3,0,n0,
			ftriangle(n,0,n3+isum1+isum3,x1,n4+n0+isum2-isum3,x2,n1-isum1,n2-isum2)*sign_(isum1+isum2)
				*fac_(n0+isum1+isum2-1)
				*Pochhammer(isum3+isum1,n3+x1*(2-d/2))
				*Pochhammer(n0-isum3+isum2,n4+x2*(2-d/2))
				*Pochhammer(-isum1-isum2-n0,4-2*(2-d/2)+n-n0-n3-n4-x1*(2-d/2)-x2*(2-d/2))
				*n0*invfac_(n0-isum3)
				*invfac_(isum1)*invfac_(isum2)*invfac_(isum3)
		)))
		+sum_(isum1,0,n1-1,sum_(isum2,0,n0-1,sum_(isum3,0,n0-1-isum2,
			ftriangle(n,n0-isum2-isum3,n3+isum1+isum3,x1,n4+n2+isum2,x2,n1-isum1,0)*sign_(n2+isum1)
				*fac_(n2+isum1+isum2+isum3-1)
				*Pochhammer(n2+isum2,n4+x2*(2-d/2))
*                                ++ -> +        
				*Pochhammer(isum3+isum1,n3+x1*(2-d/2))
				*Pochhammer(-n2-isum1-isum2-isum3,isum2+isum3+5-2*(2-d/2)+n-2*n0-n3-n4-x1*(2-d/2)-x2*(2-d/2))
				*invfac_(isum1)*invfac_(isum2)*invfac_(isum3)*invfac_(n2-1)
		)));
	repeat id Pochhammer(n?pos_,x?) = Pochhammer(n-1,x)*num(n-1+x);
	repeat id Pochhammer(n?neg_,x?) = Pochhammer(n+1,x)*den(n+x);
	id	Pochhammer(0,x?) = 1;
	id	num(x?number_) = x;
	id	den(x?number_) = 1/x;
	id	num(x?)*den(x?) = 1;
	id	num(x?) = rat(x,1);
	id	den(x?) = rat(1,x);
	id	ftriangle(n?,n0?,n3?,x1?,n4?,x2?,n1?,n2?) =
		ep`pa'^x1*ep`pb'^x2/`p'.`p'^n0/`p1'.`p1'^n1/`p2'.`p2'^n2/`pa'.`pa'^n3/`pb'.`pb'^n4;
endif;
id	xpower = 1;
*
#endprocedure
*
*--#] ntriangle : 




************************************************************
#procedure topbn3
*
* this is topbn3
*
        #-
        
        #message this is topbn3

************************************************************

        #message numerator
        
        if( count(intbn3,1));
        id  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
        id  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
        endif;
        #call ACCU(BN3 1)

        if( count(intbn3,1));
        id  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
        id  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );
        endif;
        #call ACCU(BN3 2)

        if( count(intbn3,1));
        id  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
        id  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
        endif;
        #call ACCU(BN3 3)

        if( count(intbn3,1));
        id  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
        id  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
        endif;
        #call ACCU(BN3 4)

        if( count(intbn3,1));
        id  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
        id  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
        endif;
        #call ACCU(BN3 5)

        if( count(intbn3,1));
        id  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
        id  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
        id  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
        endif;
        #call ACCU(BN3 6)

        if( count(intbn3,1));
        id  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
        id  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
        endif;
        #call ACCU(BN3 7)

*
* Warning! 
*
        if( count(intbn3,1));
        id  1/x3 = M^2 + p3.p3;
        id  1/x4 = M^2 + p4.p4;
        id  1/x5 = M^2 + p5.p5;

        id  p6.p6 = 1/x6 - M^2;
        endif;
        #call ACCU(BN3 5)

        #message do recursion

* no massive denominator:

        if( count(intbn3,1));
        if ( (count(x6,1)<=0) ) discard;
        endif;

        #call symBN3(p1,p2,p3,p4,p5,x6)

        if( count(intbn3,1));
        if ( (count(p5.p5,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        endif;
        .sort

* use the procedure triangle from MINCER to reduce the lines
* 2, 4 and 5.

* in order to avoid 1/fac_(-x)-terms:
        if( count(intbn3,1));
        if (count(p3.p3,1)>=0) multiply replace_(p1,p2,p2,p1,p4,p3,p3,p4);
        endif;
        .sort
        if( count(intbn3,1));
        if (match(1/p2.p2/p4.p4/p5.p5)>0);
        #call ntriangle(p5,p3,p1,p2,p4)
        endif;




        if ( (count(p5.p5,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        endif;
        .sort

        if( count(intbn3,1));
        if ( ( match(1/p1.p1/p2.p2/p5.p5*x6) > 0 ) && (count(p4.p4,1)>=0) );
        id p4=p1+p5;
        if ( match(1/p1.p1/p2.p2/p5.p5*x6) <= 0 ) discard;
        multiply replace_(p1,p5,p5,p2,p2,p1);

        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6) > 0 ) && (count(p5.p5,1)>=0) );
        id p5=p4-p1;
        if ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6) <= 0 ) discard;
        id p1=-p1;
        id p4=-p4;
        multiply replace_(p3,p2,p2,p3);

        elseif ( ( match(1/p3.p3/p4.p4/p5.p5*x6) > 0 ) && (count(p2.p2,1)>=0) );
        id p2=p5-p3;
        if ( match(1/p3.p3/p4.p4/p5.p5*x6) <= 0 ) discard;
        multiply replace_(p3,p5,p1,p3,p4,p2,p5,p1);

        else;

* If more massless lines than in the three cases above are >= zero,
* a massless tadpole appears.

        exit "BN3: Unknown simpler topology";
        endif;

        endif;
        .sort

        #call Conv2exact

        #message - done
#endprocedure



************************************************************

#procedure symBM (p1,p2,p3,x4,x5,x6)

        if( count(intbm,1));
        if (count(`x4',1) > count(`x6',1)) 
        multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');
        if (count(`x5',1) > count(`x6',1)) 
        multiply replace_(`x5',`x6',`x6',`x5',`p1',`p2',`p2',`p1');
        if (count(`x4',1) > count(`x5',1)) 
        multiply replace_(`x4',`x5',`x5',`x4',`p1',`p3',`p3',`p1');
        endif;        
#endprocedure




#procedure symBMnom (p1,p2,p3,p4,p5,p6,x4,x5,x6)

* sort: n6>=n5>=n4
        if( count(intbm,1));        
        if ( (count(`x4',1) > count(`x6',1)) && (count(`x4',1) > count(`x5',1)) );
        id `p1'=-`p1';
        multiply replace_(`x4',`x6',`x6',`x4',
        `p4',`p6',`p6',`p4',
        `p2',`p3',`p3',`p2');
        endif;
        if ( (count(`x5',1) > count(`x6',1)) && (count(`x5',1) > count(`x4',1)) );
        id `p3'=-`p3';
        multiply replace_(`x5',`x6',`x6',`x5',
        `p5',`p6',`p6',`p5',
        `p1',`p2',`p2',`p1');
        endif;

        if (count(`x4',1) > count(`x5',1));
        id `p4'=-`p4';
        id `p5'=-`p5';
        id `p6'=-`p6';
        multiply replace_(`x5',`x4',`x4',`x5',
        `p5',`p4',`p4',`p5',
        `p3',`p1',`p1',`p3');
        endif;
        endif;
#endprocedure


#procedure nomBM
*
* Decomposition of the numerator for type BM
*

************************************************************

        #call symBMnom(p1,p2,p3,p4,p5,p6,x4,x5,x6)
        .sort

************************************************************

        if( count(intbm,1));        
        id p4.p4 = 1/x4 - M^2;
        id p5.p5 = 1/x5 - M^2;
        id p6.p6 = 1/x6 - M^2;

        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        endif;

        #call ACCU(pi.pi)

************************************************************

        if( count(intbm,1));        
        id  p1.p2 = 1/2 * (-p3.p3 + p1.p1 + p2.p2);
        id  p1.p3 = 1/2 * ( p2.p2 - p1.p1 - p3.p3);
        id  p1.p4 = 1/2 * (-1/x6  + 1/x4  + p1.p1);
        id  p1.p5 = 1/2 * ( p3.p3 + 1/x4  - p2.p2 - 1/x6);

        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        endif;

        #call ACCU(p1)

************************************************************

        if( count(intbm,1));        
        id  p1.p6 = 1/2 * ( 1/x4  - p1.p1 - 1/x6);
        id  p2.p3 = 1/2 * (-p1.p1 + p2.p2 + p3.p3);
        id  p2.p4 = 1/2 * (-1/x5  + p2.p2 + 1/x4);
        id  p2.p5 = 1/2 * ( 1/x4  - p2.p2 - 1/x5);

        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        endif;

        #call ACCU(p1 p2)

************************************************************

        if( count(intbm,1));        
        id  p2.p6 = 1/2 * ( p3.p3 + 1/x4  - p1.p1 - 1/x5);
        id  p3.p4 = 1/2 * ( p2.p2 + 1/x6  - p1.p1 - 1/x5);
        id  p3.p5 = 1/2 * ( 1/x6  - p3.p3 - 1/x5);
        id  p3.p6 = 1/2 * (-1/x5  + p3.p3 + 1/x6);

        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        endif;

        #call ACCU(p2 p3)

************************************************************

        if( count(intbm,1));        
        id  p4.p5 = 1/2 * (-p2.p2 + 1/x4  + 1/x5 - 2*M^2);
        id  p4.p6 = 1/2 * (-p1.p1 + 1/x4  + 1/x6 - 2*M^2);
        id  p5.p6 = 1/2 * (-p3.p3 + 1/x5  + 1/x6 - 2*M^2);

        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        endif;

        #call ACCU(p4 p5 p6)

************************************************************

        #call symBM(p1,p2,p3,x4,x5,x6)


************************************************************
        
#endprocedure        






#procedure redBMn456 (p1,p2,p3,x4,x5,x6)

* Warning! In this procedure x3 is used for historical reasons!
*          Therefore not with `x3'

        if( count(intbm,1));        
        repeat;

* sort: n6>=n4,n5

                if((count(`x4',1) > 0)  &&  (count(x3,1) = 0)  &&
                (count(`x5',1) > 0)  &&  (count(`x6',1) > 0) );

                if (count(`x4',1) > count(`x6',1)) 
                multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');
                if (count(`x5',1) > count(`x6',1)) 
                multiply replace_(`x5',`x6',`x6',`x5',`p1',`p2',`p2',`p1');

* Now do the reduction via eq. (M7) resp. (4)

                if (count(`x6',1)>1);

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3?
                * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3
                * `x4'^n4  * `x5'^n5  * `x6'^(n6-1) * (-1)
                *1/(n6-1)/M^2/2
                *(
                nom(4-2*(n6-1)-n1-n3,-2)
                +n1/`p1'.`p1'*(1/`x4'-1/`x6')
                +n3/`p3'.`p3'*(1/`x5'-1/`x6')
                )
                ;
                endif;

                endif;

endrepeat;
endif;

#endprocedure

#procedure redBMn123 (p1,p2,p3,x4,x5,x6)

* Warning! In this procedure x3 is used for historical reasons!
*          Therefore not with `x3'

* at this stage: n4=n5=n6=1,0
*                n1,n2,n3 <>=0

        if( count(intbm,1));                
        repeat;

                if (count(`p3'.`p3',1) > count(`p1'.`p1',1)) 
                multiply replace_(`x4',`x5',`x5',`x4',`p1',`p3',`p3',`p1');
                if (count(`p3'.`p3',1) > count(`p2'.`p2',1)) 
                multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');


                if ( (count(`p3'.`p3',1) < 0) 
                && (count(`x4',1) = 1)       &&  (count(x3,1) = 0)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) = 1) );  

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
                * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
                * `x4'^n4  * `x5'^n5  * `x6'^n6
                *1/2*deno(n1+n2-3,2)*deno(3-2*n3,-2)*deno(n1+n2+n3-3,2)
                *(-1)
                *(
                -(`x5'*nom(-3 + n1 + n2,2)**2/`x4') - `x6'*nom(-3 + n1 + n2,2)**2/`x4' 
                +(-1)*
                n3*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x5'*`p3'.`p3') 
                -(-1)*
                n3*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x6'*`p3'.`p3') 
                -(-1)*
                n3*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x5'*`p3'.`p3') 
                +(-1)*
                n3*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x6'*`p3'.`p3') 
                +(-1)*
                n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(3 - n1 - n3,-2)/
                (2*M**2*`x5'*`p2'.`p2') 
                +(-1)*
                n2*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(3 - n1 - n3,-2)/
                (2*M**2*`x4'*`p2'.`p2') 
                +(-1)*
                n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(3 - n2 - n3,-2)/
                (2*M**2*`x6'*`p1'.`p1') 
                +(-1)*
                n1*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(3 - n2 - n3,-2)/
                (2*M**2*`x4'*`p1'.`p1') 
                +
                `x4'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/`x6' 
                +
                `x5'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/`x6' 
                +(-1)*
                n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/
                (2*M**2*`x4'*`p2'.`p2') 
                +(-1)*
                n2*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/
                (2*M**2*`x5'*`p2'.`p2') 
                +
                `x4'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/`x5' 
                +
                `x6'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/`x5' 
                +(-1)*
                n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/
                (2*M**2*`x4'*`p1'.`p1') 
                +(-1)*
                n1*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/
                (2*M**2*`x6'*`p1'.`p1') 
                +(-1)*
                `p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)*
                nom(-3 + n2 + n3,2)/M**2 
                +(-1)*
                nom(-3 + n1 + n2,2)*nom(-6 + 9*n2 - n1*n2 - 2*n2**2 + n3 + n1*n3 -
                2*n2*n3,4 - 4*n2)/(2*M**2*`x5') 
                +(-1)*
                nom(-3 + n1 + n2,2)*nom(-6 + 9*n1 - 2*n1**2 - n1*n2 + n3 - 2*n1*n3 +
                n2*n3,4 - 4*n1)/(2*M**2*`x6') 
                +(-1)*
                nom(-3 + n1 + n2,2)*nom(12 - 9*n1 + 2*n1**2 - 9*n2 + 2*n1*n2 +
                2*n2**2 - 2*n3 + n1*n3 + n2*n3,-8 + 4*n1 + 4*n2)/(2*M**2*`x4')
                );
                endif;

                id nom(x?,0)=x;
                if ( (count(`x4',1)==0) && (count(x3,1)==0) && 
                ((count(`p2'.`p2',1)>=0) || (count(`p1'.`p1',1)>=0)) ) discard;

endrepeat;
endif;
.sort

if( count(intbm,1));        
repeat;

        if (count(`p3'.`p3',1) < count(`p1'.`p1',1)) 
        multiply replace_(`x4',`x5',`x5',`x4',`p1',`p3',`p3',`p1');
        if (count(`p3'.`p3',1) < count(`p2'.`p2',1)) 
        multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');

        if ( (count(`p3'.`p3',1) > 0)
        && (count(`x4',1) = 1)       &&  (count(x3,1) = 0)  
        && (count(`x5',1) = 1)       &&  (count(`x6',1) = 1) );  

        id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
        * `x4'^n4? * `x5'^n5? * `x6'^n6?
        =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
        * `x4'^n4  * `x5'^n5  * `x6'^n6
        *deno(n1+n2-3,2)*deno(n1+n3-2,2)*deno(n2+n3-2,2)
        *(-1)
        *(
        (1 + n3)*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*`x5'*`p3'.`p3'**2) 
        -
        (1 + n3)*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*`x6'*`p3'.`p3'**2) 
        -
        (1 + n3)*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*`x5'*`p3'.`p3'**2) 
        +
        (1 + n3)*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*`x6'*`p3'.`p3'**2) 
        +
        M**2*`x5'*nom(-3 + n1 + n2,2)**2/(`x4'*`p3'.`p3') 
        +
        M**2*`x6'*nom(-3 + n1 + n2,2)**2/(`x4'*`p3'.`p3') 
        +
        n2*nom(-3 + n1 + n2,2)*nom(2 - n1 - n3,-2)/(2*`x4'*`p2'.`p2') 
        +
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(2 - n1 - n3,-2)/
        (2*`x5'*`p2'.`p2'*`p3'.`p3') 
        +
        n1*nom(-3 + n1 + n2,2)*nom(2 - n2 - n3,-2)/(2*`x4'*`p1'.`p1') 
        +
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(2 - n2 - n3,-2)/
        (2*`x6'*`p1'.`p1'*`p3'.`p3') 
        +
        n2*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(2*`x5'*`p2'.`p2') 
        -
        M**2*`x4'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(`x6'*`p3'.`p3') 
        -
        M**2*`x5'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(`x6'*`p3'.`p3') 
        +
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/
        (2*`x4'*`p2'.`p2'*`p3'.`p3') 
        +
        n1*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(2*`x6'*`p1'.`p1') 
        -
        M**2*`x4'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(`x5'*`p3'.`p3') 
        -
        M**2*`x6'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(`x5'*`p3'.`p3') 
        +
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/
        (2*`x4'*`p1'.`p1'*`p3'.`p3') 
        -
        2*M**2*nom(-3 + n1 + n2,2)*nom(1 - 2*n3,-2)*nom(-2 + n1 + n2 + n3,2)/
        `p3'.`p3' 
        +       nom(-3 + n1 + n2,2)*
        nom(-5 + n1 + 7*n2 - n1*n2 - 2*n2**2 + n3 + n1*n3 - 2*n2*n3,
        4 - 4*n2)/(2*`x5'*`p3'.`p3') 
        +
        nom(-3 + n1 + n2,2)*nom(-5 + 7*n1 - 2*n1**2 + n2 - n1*n2 + n3 -
        2*n1*n3 + n2*n3,4 - 4*n1)/(2*`x6'*`p3'.`p3') 
        +
        nom(-3 + n1 + n2,2)*nom(10 - 8*n1 + 2*n1**2 - 8*n2 + 2*n1*n2 +
        2*n2**2 - 2*n3 + n1*n3 + n2*n3,-8 + 4*n1 + 4*n2)/
        (2*`x4'*`p3'.`p3')
        );

        endif;

        id nom(x?,0)=x;
        if ( (count(`x4',1)==0) && (count(x3,1)==0) && 
        ((count(`p2'.`p2',1)>=0) || (count(`p1'.`p1',1)>=0)) ) discard;

endrepeat;
endif;

#endprocedure

#procedure redBMn25 (p1,p2,p3,x4,x5,x6)

        if( count(intbm,1));                
        repeat;
                if ( (count(`p2'.`p2',1) < 0) && (count(`p1'.`p1',1) == 0)
                &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
                * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
                * `x4'^n4  * `x5'^n5  * `x6'^n6
                *deno(4-2*n2-n4,-2)
                *(
                n4*`x4'*( `p2'.`p2' - 1/`x5' )
                )
                ;
                endif;

endrepeat;
endif;

#endprocedure

#procedure redBMn3 (p1,p2,p3,x4,x5,x6)

* Warning! In this procedure x3 is used for historical reasons!
*          Therefore not with `x3'
        
* at this stage: n4=n5=n6=1,0
*                n1,n2,n3 <>=0

        if( count(intbm,1));                
        repeat;
                
                
                if ( (count(`p3'.`p3',1) < 0)
                && (count(`p1'.`p1',1) != 0) &&  (count(`p2'.`p2',1) != 0)
                && (count(`x4',1) = 1)       &&  (count(x3,1) = 0)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) = 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
                * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
                * `x4'^n4  * `x5'^n5  * `x6'^n6
                *1/2*deno(n1+n2-3,2)*deno(3-2*n3,-2)*deno(n1+n2+n3-3,2)
                *(-1)
                *(
                -(`x5'*nom(-3 + n1 + n2,2)**2/`x4') - `x6'*nom(-3 + n1 + n2,2)**2/`x4' 
                +(-1)*
                n3*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x5'*`p3'.`p3') 
                -(-1)*
                n3*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x6'*`p3'.`p3') 
                -(-1)*
                n3*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x5'*`p3'.`p3') 
                +(-1)*
                n3*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x6'*`p3'.`p3') 
                +(-1)*
                n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(3 - n1 - n3,-2)/
                (2*M**2*`x5'*`p2'.`p2') 
                +(-1)*
                n2*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(3 - n1 - n3,-2)/
                (2*M**2*`x4'*`p2'.`p2') 
                +(-1)*
                n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(3 - n2 - n3,-2)/
                (2*M**2*`x6'*`p1'.`p1') 
                +(-1)*
                n1*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(3 - n2 - n3,-2)/
                (2*M**2*`x4'*`p1'.`p1') 
                +
                `x4'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/`x6' 
                +
                `x5'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/`x6' 
                +(-1)*
                n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/
                (2*M**2*`x4'*`p2'.`p2') 
                +(-1)*
                n2*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/
                (2*M**2*`x5'*`p2'.`p2') 
                +
                `x4'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/`x5' 
                +
                `x6'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/`x5' 
                +(-1)*
                n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/
                (2*M**2*`x4'*`p1'.`p1') 
                +(-1)*
                n1*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/
                (2*M**2*`x6'*`p1'.`p1') 
                +(-1)*
                `p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)*
                nom(-3 + n2 + n3,2)/M**2 
                +(-1)*
                nom(-3 + n1 + n2,2)*nom(-6 + 9*n2 - n1*n2 - 2*n2**2 + n3 + n1*n3 -
                2*n2*n3,4 - 4*n2)/(2*M**2*`x5') 
                +(-1)*
                nom(-3 + n1 + n2,2)*nom(-6 + 9*n1 - 2*n1**2 - n1*n2 + n3 - 2*n1*n3 +
                n2*n3,4 - 4*n1)/(2*M**2*`x6') 
                +(-1)*
                nom(-3 + n1 + n2,2)*nom(12 - 9*n1 + 2*n1**2 - 9*n2 + 2*n1*n2 +
                2*n2**2 - 2*n3 + n1*n3 + n2*n3,-8 + 4*n1 + 4*n2)/(2*M**2*`x4')
                );
                endif;
                
                id nom(x?,0)=x;
                if ( (count(`x4',1)==0) && (count(x3,1)==0) && 
                ((count(`p2'.`p2',1)>=0) || (count(`p1'.`p1',1)>=0)) ) discard;
                
        endrepeat;
        endif;
        
#endprocedure

#procedure redBMn35 (p1,p2,p3,x4,x5,x6)

        if( count(intbm,1));                
        repeat;
                if ( (count(`p3'.`p3',1) < 0) && (count(`p1'.`p1',1) == 0)
                &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
                * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
                * `x4'^n4  * `x5'^n5  * `x6'^n6
                *deno(4-2*n3-n6,-2)
                *(
                n6*`x6'*( `p3'.`p3' - 1/`x5' )
                )
                ;
                endif;
                
        endrepeat;
        endif;

        #call Conv2exact        
#endprocedure

#procedure redBMn3p (p1,p2,p3,x4,x5,x6)

        if( count(intbm,1));                
        repeat;

                if (count(`p3'.`p3',1) < count(`p1'.`p1',1)) 
                multiply replace_(`x4',`x5',`x5',`x4',`p1',`p3',`p3',`p1');
                if (count(`p3'.`p3',1) < count(`p2'.`p2',1)) 
                multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');

                if ( (count(`p3'.`p3',1) > 0)
                && (count(`x4',1) = 1)       &&  (count(x3,1) = 0)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) = 1) );  

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
                * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
                * `x4'^n4  * `x5'^n5  * `x6'^n6
                *deno(n1+n2-3,2)*deno(n1+n3-2,2)*deno(n2+n3-2,2)
                *(-1)
                *(
                (1 + n3)*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*`x5'*`p3'.`p3'**2) 
                -
                (1 + n3)*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*`x6'*`p3'.`p3'**2) 
                -
                (1 + n3)*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*`x5'*`p3'.`p3'**2) 
                +
                (1 + n3)*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*`x6'*`p3'.`p3'**2) 
                +
                M**2*`x5'*nom(-3 + n1 + n2,2)**2/(`x4'*`p3'.`p3') 
                +
                M**2*`x6'*nom(-3 + n1 + n2,2)**2/(`x4'*`p3'.`p3') 
                +
                n2*nom(-3 + n1 + n2,2)*nom(2 - n1 - n3,-2)/(2*`x4'*`p2'.`p2') 
                +
                n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(2 - n1 - n3,-2)/
                (2*`x5'*`p2'.`p2'*`p3'.`p3') 
                +
                n1*nom(-3 + n1 + n2,2)*nom(2 - n2 - n3,-2)/(2*`x4'*`p1'.`p1') 
                +
                n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(2 - n2 - n3,-2)/
                (2*`x6'*`p1'.`p1'*`p3'.`p3') 
                +
                n2*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(2*`x5'*`p2'.`p2') 
                -
                M**2*`x4'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(`x6'*`p3'.`p3') 
                -
                M**2*`x5'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(`x6'*`p3'.`p3') 
                +
                n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/
                (2*`x4'*`p2'.`p2'*`p3'.`p3') 
                +
                n1*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(2*`x6'*`p1'.`p1') 
                -
                M**2*`x4'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(`x5'*`p3'.`p3') 
                -
                M**2*`x6'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(`x5'*`p3'.`p3') 
                +
                n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/
                (2*`x4'*`p1'.`p1'*`p3'.`p3') 
                -
                2*M**2*nom(-3 + n1 + n2,2)*nom(1 - 2*n3,-2)*nom(-2 + n1 + n2 + n3,2)/
                `p3'.`p3' 
                +       nom(-3 + n1 + n2,2)*
                nom(-5 + n1 + 7*n2 - n1*n2 - 2*n2**2 + n3 + n1*n3 - 2*n2*n3,
                4 - 4*n2)/(2*`x5'*`p3'.`p3') 
                +
                nom(-3 + n1 + n2,2)*nom(-5 + 7*n1 - 2*n1**2 + n2 - n1*n2 + n3 -
                2*n1*n3 + n2*n3,4 - 4*n1)/(2*`x6'*`p3'.`p3') 
                +
                nom(-3 + n1 + n2,2)*nom(10 - 8*n1 + 2*n1**2 - 8*n2 + 2*n1*n2 +
                2*n2**2 - 2*n3 + n1*n3 + n2*n3,-8 + 4*n1 + 4*n2)/
                (2*`x4'*`p3'.`p3')
                );

                endif;

                id nom(x?,0)=x;
                if ( (count(`x4',1)==0) && (count(x3,1)==0) && 
                ((count(`p2'.`p2',1)>=0) || (count(`p1'.`p1',1)>=0)) ) discard;

endrepeat;
endif;

#endprocedure

#procedure redBMn6 (p1,p2,p3,x4,x5,x6)

        #do i = 1,1

* Now do the reduction via eq. (M7) resp. (4)

                if( count(intbm,1));        
                if ( (count(`x4',1)>=1) && (count(`x5',1)>=1) && (count(`x6',1)>1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3?
                * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3
                * `x4'^n4  * `x5'^n5  * `x6'^(n6-1) * (-1)
                *1/(n6-1)/M^2/2
                *(
                nom(4-2*(n6-1)-n1-n3, -2)
                +n1/`p1'.`p1'*(1/`x4'-1/`x6')
                +n3/`p3'.`p3'*(1/`x5'-1/`x6')
                )
                ;

                redefine i "0";

                endif;
                endif;        
                .sort
                
        #enddo
        
#endprocedure

#procedure redBMn6exp (p1,p2,p3,x4,x5,x6)
        
        #do i = 1,1

* Now do the reduction via eq. (M7) resp. (4)

                if( count(intbm,1));                
                if ( (count(`x4',1)>=1) && (count(`x5',1)>=1) && (count(`x6',1)>1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3?
                * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3
                * `x4'^n4  * `x5'^n5  * `x6'^(n6-1) * (-1)
                *1/(n6-1)/M^2/2
                *(
                num(4-2*(n6-1)-n1-n3 -2*(2-d/2))
                +n1/`p1'.`p1'*(1/`x4'-1/`x6')
                +n3/`p3'.`p3'*(1/`x5'-1/`x6')
                )
                ;
                redefine i "0";

                endif;
                endif;
               
                #call ACCU{BMn6}
                
        #enddo
        
#endprocedure


#procedure symmetryBM
*
* symmetryBM
*
        if( count(intbm,1));                        
        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        endif;        
        .sort

*
* sort n6>=n5>=n4
*
        if( count(intbm,1));                        
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) ) multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;        
        .sort

        if( count(intbm,1));                        
        if ( (count(x4,1)<=0) || (count(x5,1)<=0) || (count(x6,1)<=0) );
        endif;        
*
* sort n6>n5>n4
*
        if( count(intbm,1));                        
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) )
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;
        endif;        
        .sort
        
#endprocedure        


#procedure reduceBMBM
*
* reduceBMBM
*

* This program reduces the integrals of type BM to easy integrals of
* the type M1 and T.

************************************************************


************************************************************

* massive lines

        #call symmetryBM

* reduce n6 (expand in ep)
        #call redBMn6exp(p1,p2,p3,x4,x5,x6)
        .sort
        
        #call symmetryBM

* now: n4=1 (or one of the massive lines = 0)

        #call redBMn6exp(p1,p2,p3,x4,x5,x6)
        .sort

        
        #call symmetryBM
        
* now: n4=n5=1

        #call redBMn6exp(p1,p2,p3,x4,x5,x6)
        .sort

        
        #call symmetryBM

************************************************************

* massless lines

* first treat all (!) massless indices which are <0

        #call redBMn3p(p1,p2,p3,x4,x5,x6)
        .sort
        
        #call symmetryBM

        if( count(intbm,1));                        
        if ( (count(x4,1)>0) && (count(x5,1)>0) && (count(x6,1)>0) );
        
* sort: |n3|<=|n1|,|n2| ?
        
        
        if (count(p3.p3,1) < count(p1.p1,1)) 
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        if (count(p3.p3,1) < count(p2.p2,1)) 
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        endif;
        endif;
* the following procedure is only used, if all 3 indices are <0
        
        #call redBMn3(p1,p2,p3,x4,x5,x6)
        .sort

        #call symmetryBM

        if( count(intbm,1));                                
        if ( (count(x4,1)>0) && (count(x5,1)>0) && (count(x6,1)>0) );
        if (count(p3.p3,1)==0) multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;
        endif;        
        .sort

* now one massless index is =0; here n1=0


* Now treat the two other massless lines with simpler
* rec. relations.

                
        #call redBMn35(p1,p2,p3,x4,x5,x6)
        .sort

        
        
        
* Don't use the file 'symmetryBM' here because the next procedure 
* treats n2>0.

        if( count(intbm,1));                                
        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        endif;        
        .sort

        if( count(intbm,1));                                
        if ( (count(x4,1)<=0) || (count(x5,1)<=0) || (count(x6,1)<=0) );
*
* sort n6>n5>n4
*
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) )
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;
        endif;        
        .sort


        #call redBMn25(p1,p2,p3,x4,x5,x6)
        .sort

        
        #call symmetryBM

* if n1=n2=n3=0 the following rec. rel. is very simple

        #CALL redBMn456(p1,p2,p3,x4,x5,x6)

        #call symmetryBM

        #call Conv2exact        
#endprocedure        





************************************************************
#procedure topbm
*
* this is topbmbm
*

************************************************************

        #message this is topbmbm
        
        #message numerator
        
        #call nomBM
        
        #message do recursion
        
        #call reduceBMBM
        
        #call Conv2exact
        
        
        if( count(intbm,1));        
        id x4^n4?neg_=(p4.p4+M^2)^-n4;
        id p4 = p2+p5;
        id p5.p5 = 1/x5 - M^2;
        
        if ( (count(x4,1)=0) && (count(x3,1)=0) 
        && ((count(p2.p2,1)>=0) || (count(p1.p1,1)>=0)) ) discard;
* this is necessary, because it is possible, that n1 and n2 shrink to
* a point!
        if (count(x5,1)<=0) discard;
        if (count(x6,1)<=0) discard;
        endif;        
        .sort

        #call Conv2exact()        


        .sort

        #message - done
        
#endprocedure        



#procedure symBM1 (p1,p2,p3,p4,x5,x6)
        if( count(intbm1,1));
        if (count(`x6',1) > count(`x5',1))
        multiply replace_(`p1',`p2',`p2',`p1',`x5',`x6',`x6',`x5');
        endif;
#endprocedure


#procedure redBM1n12 (p1,p2,p3,p4,x5,x6)
        
        #do i=1,1
                if( count(intbm1,1));
                if ( count(`p1'.`p1',1) > count(`p2'.`p2',1) ) 
                multiply replace_(`p1',`p2',`p2',`p1',`x5',`x6',`x6',`x5');
                endif;
                #call ACCU{BM1}
                
                if( count(intbm1,1));        
                if (match(1/`p1'.`p1'^2/`p2'.`p2'/`p4'.`p4'*`x5'*`x6')>0);
                
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4  * `x5'^n5  * `x6'^n6
                *(-1)/(n1-1)/M^2*(
                - `p1'.`p1'*n3/`p3'.`p3'*(1/`x6'-1/`x5')
                - (n1-1)                *(1/`x6'-`p4'.`p4')
                + `p1'.`p1'*num(d-2*n6-n3-(n1-1))
                +2*n6*M^2*`x6'*`p1'.`p1'
                );
                
                redefine i "0";
                
                endif;
                
                if (match(1/`p1'.`p1'/`p2'.`p2'^2/`p4'.`p4'*`x5'*`x6')>0);
                
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4  * `x5'^n5  * `x6'^n6
                *(-1)/(n2-1)/M^2*(
                - `p2'.`p2'*n3/`p3'.`p3'*(1/`x5'-1/`x6')
                - (n2-1)                *(1/`x5'-`p4'.`p4')
                + `p2'.`p2'*num(d-2*n5-n3-(n2-1))
                +2*n5*M^2*`x5'*`p2'.`p2'
                );
                
                redefine i "0";

                endif;

* topBM1        
                endif;
                #call ACCU{BM1}
                
        #enddo
        
#endprocedure



#procedure redBM1n124 (p1,p2,p3,p4,x5,x6)

        #do i=1,1
                if( count(intbm1,1));
                if (match(1/`p1'.`p1'/`p2'.`p2'/`p4'.`p4'*`x5'*`x6')>0);
                
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4  * `x5'^n5  * `x6'^n6
                *(-1)*deno(2+n4-n1-n2-n3,-1)
                *(
                n5*`x5'*(`p4'.`p4'-`p2'.`p2')
                +n6*`x6'*(`p4'.`p4'-`p1'.`p1')
                )
                ;
                
                redefine i "0";
                
                endif;
                endif;
                #call ACCU{BM1}
                
        #enddo
        
#endprocedure



#procedure redBM1n3 (p1,p2,p3,p4,x5,x6)

        if( count(intbm1,1));
        if ( (match(1/`p2'.`p2'/`p3'.`p3'/`p4'.`p4'*`x5'*`x6')>0) 
        && (count(`p1'.`p1',1)>=0) );
        if ( (count(`p3'.`p3',1) < -1 ) && (count(`p1'.`p1',1) > 0 ) );
        id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
        /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
        =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
        /`p3'.`p3'^n3  /`p4'.`p4'^n4 * `x5'^n5 * `x6'^n6
        *(-1)/(n3-1)
        *(
        `p3'.`p3'/`p1'.`p1'*nom(4-2*n2-(n3-1)-n5,-2)
        -(n3-1)*`p2'.`p2'/`p1'.`p1'
        +`p3'.`p3'/`p1'.`p1'*n5*`x5'*(`p4'.`p4'-`p2'.`p2'+M^2)
        )
        ;
        
        redefine i "0";
        
        endif;
        endif;
        
        if ( (match(1/`p1'.`p1'/`p3'.`p3'/`p4'.`p4'*`x5'*`x6')>0) 
        && (count(`p2'.`p2',1)>=0) );
        if ( (count(`p3'.`p3',1) < -1 ) && (count(`p2'.`p2',1) > 0 ) );
        id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
        /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
        =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
        /`p3'.`p3'^n3  /`p4'.`p4'^n4 * `x5'^n5 * `x6'^n6
        *(-1)/(n3-1)
        *(
        `p3'.`p3'/`p2'.`p2'*nom(4-2*n1-(n3-1)-n6,-2)
        -(n3-1)*`p1'.`p1'/`p2'.`p2'
        +`p3'.`p3'/`p2'.`p2'*n6*`x6'*(`p4'.`p4'-`p1'.`p1'+M^2)
        )
        ;
        endif;
        endif;
* topBM1        
        endif;        
#endprocedure

#procedure redBM1n35 (p1,p2,p3,p4,x5,x6)
        
        #do i=1,10
                
                if( count(intbm1,1));        
                if ( (match(1/`p2'.`p2'/`p3'.`p3'/`p4'.`p4'*`x5'*`x6')>0) 
                && (count(`p1'.`p1',1)>=0) );
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4 * `x5'^n5 * `x6'^n6
                *(-1)*deno(4-2*n3-n1-n6,-2)
                *(
                n1/`p1'.`p1'*(`p2'.`p2'-`p3'.`p3')
                +n6*`x6'*(1/`x5'-`p3'.`p3')
                )
                ;
                
                redefine i "0";
                
                endif;
                endif;
                #call ACCU{BM1}
                
        #enddo
        
#endprocedure


#procedure redBM1n36 (p1,p2,p3,p4,x5,x6)
        
* with this procedure n3 or n6 is reduced to zero.
        
        #do i=1,1
                
                if( count(intbm1,1));        
                if ( (match(1/`p1'.`p1'/`p3'.`p3'/`p4'.`p4'*`x5'*`x6')>0) 
                && (count(p2.p2,1)>=0) );
                
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4 * `x5'^n5 * `x6'^n6
                *(-1)*deno(4-2*n3-n2-n5,-2)
                *(
                n2/`p2'.`p2'*(`p1'.`p1'-`p3'.`p3')
                +n5*`x5'*(1/`x6'-`p3'.`p3')
                )
                ;
                
                redefine i "0";
                
                endif;
                endif;
                
                #call ACCU{BM1}
                
        #enddo
        
#endprocedure





************************************************************
#procedure topbm1
*
* this is topbm1
*
        #-
        #message this is topbm1
        
************************************************************
        
        
************************************************************
        
        #message numerator
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        ID  p6.p6 = 1/x6 - M^2;
        ID  p5.p5 = 1/x5 - M^2;
        endif;        
        #call ACCU(BM1 0)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        id p1=p2-p3;
        id p6=p3+p5;
        endif;        
        #call ACCU(BM1 1)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        id  p2.p3 = 1/2 * (-p1.p1 + p2.p2 + p3.p3);
        id  p2.p4 = 1/2 * (-1/x5  + p2.p2 + M^2+p4.p4);
        endif;        
        #call ACCU(BM1 2)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        id  p3.p4 = 1/2 * ( p2.p2 + 1/x6  - p1.p1 - 1/x5);
        endif;        
        #call ACCU(BM1 3)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        id  p2.p5 = 1/2 * ( M^2+p4.p4  - p2.p2 - 1/x5);
        id  p3.p5 = 1/2 * ( 1/x6  - p3.p3 - 1/x5);
        endif;        
        #call ACCU(BM1 5)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        
        #call ACCU(BM1 6)
        
        if( count(intbm1,1));
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;
        .sort
        if( count(intbm1,1));
        id  p4.p5 = 1/2 * (-p2.p2 + p4.p4 + 1/x5 - M^2);
        endif;        
        #call ACCU(BM1 7)

        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort

*
* Warning!
*
        if( count(intbm1,1));        
        ID  p6.p6 = 1/x6 - M^2;
        ID  p5.p5 = 1/x5 - M^2;
        endif;        
        #call ACCU(BM1 8)
        
        if( count(intbm1,1));       
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort

        #call ACCU(BM1)

        #message do recursion


        #call symBM1(p1,p2,p3,p4,x5,x6)
        #call redBM1n12(p1,p2,p3,p4,x5,x6)

        #call symBM1(p1,p2,p3,p4,x5,x6)
        #call redBM1n124(p1,p2,p3,p4,x5,x6)
        .sort

* from now on n1, n2 or n4 is <= zero. 
        if( count(intbm1,1));
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        

        #call Conv2exact
        .sort

        #do i=1,1
                #call symBM1(p1,p2,p3,p4,x5,x6)
                #call redBM1n3(p1,p2,p3,p4,x5,x6)
                if( count(intbm1,1));        
                if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
                if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
                if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
                if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
                if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
                if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
                if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
                if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
                if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
                if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
                endif;        
                #call ACCU(BM1_n3)
        #enddo
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        #call ACCU(BM1_n3)

        #call redBM1n35(p1,p2,p3,p4,x5,x6)

        
        #call redBM1n36(p1,p2,p3,p4,x5,x6)

        if( count(intbm1,1));        
        if ( (match(1/p1.p1/p2.p2*x5*x6)>0) && (count(p4.p4,1)>=0) );
        id p4=p2+p5;
        if (match(1/p1.p1/p2.p2*x5*x6)<=0) discard;
        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x5) > 0 )  && (count(x6,1)<=0) );
        id 1/x6 = M^2 + p6.p6;
        id p6=p3+p5;
        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6) > 0 )  && (count(x5,1)<=0) );
        id 1/x5 = M^2 + p5.p5;
        id p5=p2+p4;
***else;
***  multiply 1/(1-1);
        endif;
        endif;        
        .sort

* treat negative powers of massive denominators:

        if( count(intbm1,1));        
        id x5^n5?neg_=(p4.p4-2*p4.p2+p2.p2+M^2)^-n5;
        id x6^n6?neg_=(p4.p4-2*p4.p1+p1.p1+M^2)^-n6;
        endif;        
        .sort

        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort        
        
* now, we have reduced diagrams with the following lines != 0:

* 1,2,(3),5,6 : M1
* 2,4,5,6   : M4
* 2,3,4,6   : M2
* 1,4,5,6   : M4
* 1,3,4,5   : M2
* 1,2,3,4,5   : M2
* 1,2,3,4,6   : M2

* the notation has to be adjusted

        if( count(intbm1,1));        
        if ( (match(1/p2.p2/p4.p4*x5*x6) > 0) && (count(p1.p1,1)>=0) 
        && (count(p3.p3,1)>=0) );
        id p1=p4-p6;
        id p3=p6-p5;
        id p2=-p2;
        if ( match(1/p2.p2/p4.p4*x5*x6) <= 0 ) discard;
        multiply replace_(p6,p5,p5,p6,p2,p3,x6,x5,x5,x6);
        Multiply intm4/intbm1;

        elseif ( (match(1/p2.p2/p3.p3/p4.p4*x6) > 0)  && (count(p1.p1,1)>=0) 
        && (count(x5,1)<=0) );
        id p1=p4-p6;
        id p5=p6-p3;
        if ( match(1/p2.p2/p3.p3/p4.p4*x6) <= 0 ) discard;
        multiply replace_(p2,p1,p3,p2,p4,p5);
        Multiply intm2/intbm1;

        elseif ( (match(1/p1.p1/p4.p4*x5*x6) > 0)  && (count(p3.p3,1)>=0) 
        && (count(p2.p2,1)>=0) );
        id p2=p4-p5;
        id p3=p6-p5;
        id p6=-p6;
        id p4=-p4;
        if ( match(1/p1.p1/p4.p4*x5*x6) <= 0 ) discard;
        multiply replace_(p1,p3);
        Multiply intm4/intbm1;
        
        elseif ( (match(1/p1.p1/p3.p3/p4.p4*x5) > 0 ) && (count(p2.p2,1)>=0) 
        && (count(x6,1)<=0) );;
        id p2=p1+p3;
        id p6=p3+p5;
        id p3=-p3;
        if ( match(1/p1.p1/p3.p3/p4.p4*x5) <= 0 ) discard;
        multiply replace_(p5,p6,p3,p2,p4,p5,x5,x6);
        Multiply intm2/intbm1;

        elseif ( ( match(1/p1.p1/p2.p2*x5*x6) > 0 )  && (count(p4.p4,1)>=0) );

* This is type M1 and the notation is already o.k.
        Multiply intm1/intbm1;

        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x5) > 0 )  && (count(x6,1)<=0) );
        
* This is type M2

        id p1=-p1;
        id p2=-p2;
        id p3=-p3;
        id p4=-p4;
        id p5=-p5;
        multiply replace_(p4,p5,p1,p2,p3,p1,p2,p3,p5,p6,x5,x6);
        Multiply intm2/intbm1;

        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6) > 0 )  && (count(x5,1)<=0) );

* This is type M2

        id p2=-p2;
        id p4=-p4;
        id p6=-p6;
        multiply replace_(p4,p5,p1,p3,p3,p1);
        Multiply intm2/intbm1;

        else;
        exit "BM1: Unknown simpler topology";
        endif;
* topBM1        
        endif;
        #call ACCU(BM1)

        #message - done
        
#endprocedure        


#procedure symBM2 (p1,p2,p3,p4,p5,x6)

        if( count(intbm2,1));        
        if (count(`p3'.`p3',1) >= count(`p1'.`p1',1))
        multiply replace_(`p1',`p3',`p3',`p1',`p4',`p5',`p5',`p4');
        endif;
#endprocedure




************************************************************
#procedure topbm2
*
* this is topbm2 
*
        #-
        
        #message this is topbm2
        
************************************************************
        
        
************************************************************

        #message numerator
        
        if( count(intbm2,1));        
        id  p6.p6 = 1/x6 - M^2;
        endif;        
        #call ACCU(BM2 0)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;
        .sort
        
        if( count(intbm2,1));                
        id  p1.p2 = 1/2 * (-p3.p3 + p1.p1 + p2.p2);
        id  p1.p3 = 1/2 * ( p2.p2 - p1.p1 - p3.p3);
        id  p1.p4 = 1/2 * (-1/x6  + 1/x4  + p1.p1);
        endif;        
        #call ACCU(BM2 1)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm2,1));                
        id  p1.p5 = 1/2 * ( p3.p3 + 1/x4  - p2.p2 - 1/x6);
        id  p1.p6 = 1/2 * ( 1/x4  - p1.p1 - 1/x6);
        id  p2.p3 = 1/2 * (-p1.p1 + p2.p2 + p3.p3);
        endif;        
        #call ACCU(BM2 2)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm2,1));                
        id  p2.p4 = 1/2 * (-1/x5  + p2.p2 + 1/x4);
        id  p2.p5 = 1/2 * ( 1/x4  - p2.p2 - 1/x5);
        id  p2.p6 = 1/2 * ( p3.p3 + 1/x4  - p1.p1 - 1/x5);
        endif;        
        #call ACCU(BM2 3)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort        
        
        if( count(intbm2,1));                
        id  p3.p4 = 1/2 * ( p2.p2 + 1/x6  - p1.p1 - 1/x5);
        id  p3.p5 = 1/2 * ( 1/x6  - p3.p3 - 1/x5);
        id  p3.p6 = 1/2 * (-1/x5  + p3.p3 + 1/x6);
        endif;        
        #call ACCU(BM2 4)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm2,1));                
        id  p4.p5 = 1/2 * (-p2.p2 + 1/x4  + 1/x5 - 2*M^2);
        id  p4.p6 = 1/2 * (-p1.p1 + 1/x4  + 1/x6 - 2*M^2);
        id  p5.p6 = 1/2 * (-p3.p3 + 1/x5  + 1/x6 - 2*M^2);
        endif;        
        #call ACCU(BM2 5)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort

*
* Warning! 
*
        if( count(intbm2,1));                
        id  1/x4 = M^2 + p4.p4;
        id  1/x5 = M^2 + p5.p5;
        endif;        
        #call ACCU(BM2 6)
        
        if( count(intbm2,1));                
        id  p6.p6 = 1/x6 - M^2;
        endif;        
        #call ACCU(BM2 7)
        
        #message do recursion

* no massive denominator:
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        endif;
        
        #call symBM2(p1,p2,p3,p4,p5,x6)
        
        if( count(intbm2,1));                
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
* Use the procedure triangle from MINCER to reduce the lines
* 1, 2 and 4.
        
* in order to avoid 1/fac_(-x)-terms:

        if( count(intbm2,1));                
        if (count(p5.p5,1)>=0) multiply replace_(p5,p4,p4,p5,p1,p3,p3,p1);
        
        if (match(1/p1.p1/p2.p2/p4.p4)>0);
        #call ntriangle(p2,p3,p5,p1,p4)        
        endif;
        
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm2,1));        
        if ( ( match(1/p2.p2/p3.p3/p4.p4*x6) > 0 ) && (count(p1.p1,1)>=0) );
        id p1=p2-p3;
        if ( match(1/p2.p2/p3.p3/p4.p4*x6) <= 0 ) discard;
        multiply replace_(p5,p3,p3,p5,p4,p2,p2,p1);
        elseif ( ( match(1/p1.p1/p3.p3/p4.p4/p5.p5*x6) > 0 ) && (count(p2.p2,1)>=0) );
        id p2=p1+p3;
        id p5=-p5;
        if ( match(1/p1.p1/p3.p3/p4.p4/p5.p5*x6) <= 0 ) discard;
        multiply replace_(p4,p2,p5,p4);
        elseif ( ( match(1/p1.p1/p2.p2/p5.p5*x6) > 0 ) && (count(p4.p4,1)>=0) );
        id p4=p2+p5;
        if ( match(1/p1.p1/p2.p2/p5.p5*x6) <= 0 ) discard;
        else;
        
* If more massless lines than in the three cases above are >= zero,
* a massless tadpole appears.
        
***  discard;
        multiply 1/(1-1);
        endif;
        endif;        
        .sort
        
        #message - done
        
#endprocedure        



#procedure nomgm3(v1,v2,v3,x,y,z,in)
* Apply only to topo "in"
        if(count(int`in',1)) id  `v3' = -`v1'-`v2';
        
        #call ACCU{nomgm3 1}
        
        if(count(int`in',1)) id  `v1'.`v1' = 1/`x' - M^2;
        
        #call ACCU{nomgm3 2}
        
        if(count(int`in',1)) id  `v2'.`v2' = 1/`y' - M^2;
        
        #call ACCU{nomgm3 3}
        
        if(count(int`in',1)) id  `v1'.`v2' = (1/(`z') - 1/`x' - 1/`y' +2*M^2)/2;
        
        #call ACCU{nomgm3 4}
        
#endprocedure



************************************************************
#procedure topm1
*
* topm1
* (simple topology)
*

        #message this is topm1
        
        if( count(intm1,1));                
        id p1=-p1;
        endif;                
        #call IntOne(p2,p1,p3,m1,MM0)
        .sort
        
        #call ACCU(one)
        
        if( count(intMM0,1));
        id p6=-p6;
        endif;                
        #call nomgm3(p5,p6,p3,x5,x6,1/p3.p3,MM0)
        
        .sort
        
        #call TadpoleMM0(x5,x6,p3,MM0,0)
        .sort

        #call ACCU(gm3)

#endprocedure        


************************************************************
#procedure topm2
*
* topm2
* (simple topology)
*
        
        #message this is topm2
*         
* Procedure m2 modified        
* IntOne can deal only with numerator containing momentum at first position
*         
        if( count(intm2,1));        
***         id p1=-p1;
        id p1=p2-p3;
        endif;        
        #call IntOne(p2,p1,p3,m2,M00)        
        .sort

        if( count(intM00,1));        
        id p5=p6-p3;
        endif;        
        
        #call IntOne(p3,p5,p6,M00,M0)        
        .sort

        #call averts(p6,M0)
        .sort

        #call TadpoleM0(x6,p6,M0,0)        
        .sort

#endprocedure        



************************************************************
#procedure topm3
*
* topm3
* (simple topology)
*

        #message this is topm3

        if( count(intm3,1)) id p1=p2-p6;

        #call IntOne(p2,p1,p6,m3,M00)                
        .sort

        if( count(intM00,1)) id p4=p3-p6;
        
        #call IntOne(p3,p4,p6,M00,M0)        
        .sort
        #call averts(p6,M0)
        .sort

        #call TadpoleM0(x6,p6,M0,0)        
        .sort

#endprocedure        


************************************************************
#procedure topm4
*
* topm4
* (simple topology)
*
        
        #message this is topm4

* Numerator p4->p3,p6
        if(count(intm4,1)) id p4=p6-p3;        
        #call IntOne(p3,p4,p6,m4,MxM)
        .sort
        
        #call averts(p5,MxM)
        .sort
        #call TadpoleM0(x5,p5,MxM,M0)        
        .sort
        
        #call averts(p6,M0)
        .sort
        
        #call TadpoleM0(x6,p6,M0,0)        
        .sort
        
#endprocedure        


************************************************************
#procedure topm5
*
* this is topm5
*

        #message this is topm5

************************************************************

* treat 1-loop integral

        if( count(intm5,1)) multiply replace_(x1,s1m,x2,s2m,x3,s3m,x4,s4m,x5,s5m,x6,s6m);

        #call averts(p4,m5)
        #call TadpoleM0(s4m,p4,m5,MMM)
        .sort

* treat numerator of 2-loop integral

        
        if( count(intMMM,1)) id p1.p1 = 1/s1m - M^2;
        .sort
        if( count(intMMM,1)) id p1 = p2+p5;
        .sort

        if( count(intMMM,1));        
        id p2.p2 = 1/s2m - M^2;
        id p5.p5 = 1/s5m - M^2;
        endif;        
        .sort

        if( count(intMMM,1)) id p2.p5 = 1/2 * (p1.p1-p2.p2-p5.p5);
        .sort

        if( count(intMMM,1));        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p5.p5 = 1/s5m - M^2;
        endif;        
        .sort

* recursion for 2-loop integral
        
* the following part is identical to topT111
        
        #do i=1,5
                
                if( count(intMMM,1));        
                if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s5m,1)!=0) ); 
                multiply replace_(test5,s1m,test6,s2m);
                if ( (count(s1m,1)>1) && (count(s2m,1)>=1) && (count(s5m,1)>=1) );
                id s1m^n1?*s2m^n2?*s5m^n5? = 
                s1m^n1*s2m^n2*s5m^n5 * (-1) * 1/3/(n1-1)/M^2 * (
                nom(4+3-3*n1,-2)/s1m
                +2*n2*s2m/s1m*(1/s5m-1/s1m)
                -(n1-1)*(1/s5m-1/s2m)
                );
                endif;
                if ( count(s1m,1) < count(s2m,1) ) multiply replace_(s1m,s2m,s2m,s1m);
                if ( count(s1m,1) < count(s5m,1) ) multiply replace_(s1m,s5m,s5m,s1m);
                if ( count(s2m,1) < count(s5m,1) ) multiply replace_(s2m,s5m,s5m,s2m);
                endif;
                endif;        

                #call Conv2exact
                .sort
                
        #enddo

        if( count(intMMM,1));
        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s5m,1)!=0) ); 
        multiply replace_(test5,s1m,test6,s2m);
        repeat;
                if ( (count(s1m,1)>1) && (count(s2m,1)>=1) && (count(s5m,1)>=1) );
                id s1m^n1?*s2m^n2?*s5m^n5? = 
                s1m^n1*s2m^n2*s5m^n5 * (-1) * 1/3/(n1-1)/M^2 * (
                nom(4+3-3*n1,-2)/s1m
                +2*n2*s2m/s1m*(1/s5m-1/s1m)
                -(n1-1)*(1/s5m-1/s2m)
                );
                endif;
                if ( count(s1m,1) < count(s2m,1) ) multiply replace_(s1m,s2m,s2m,s1m);
                if ( count(s1m,1) < count(s5m,1) ) multiply replace_(s1m,s5m,s5m,s1m);
                if ( count(s2m,1) < count(s5m,1) ) multiply replace_(s2m,s5m,s5m,s2m);
        endrepeat;
        endif;
        endif;

        .sort

        #message perform integration

        if( count(intMMM,1));  
        if ( count(s1m,1) < count(s2m,1) ) multiply replace_(s1m,s2m,s2m,s1m);
        if ( count(s1m,1) < count(s5m,1) ) multiply replace_(s1m,s5m,s5m,s1m);
        if ( count(s2m,1) < count(s5m,1) ) multiply replace_(s2m,s5m,s5m,s2m);

        id 1/s5m = p5.p5 + M^2;
        endif;  
        .sort

        if( count(intMMM,1));  
        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s5m,1)==0) ); 
***  #call nomgm3T(p1,p2,p5,test5,test6,1/p5.p5)
        multiply replace_(test5,s1m,test6,s2m)*intMM0/intMMM;

        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s5m,1)==0) );
        multiply intM00/intm5;
        
        elseif ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s5m,1)==1) ); 
***  id s1m*s2m*s5m = -3*deno(1,-2)*M^2*(1/2/ep^2 + 1/2/ep) +  M^2*g3fl1;
*** constant form Davydychev, Tausk, NPB397(93)123
        id s1m*s2m*s5m = M^2*miT111*int0/intMMM;

        else;
        multiply 1/(1-1);
        endif;
        endif;
*       MM0  

        #call TadpoleMM0(s1m,s2m,p5,MM0,0)
*       M00
        #call IntOne(p2,p5,p1,M00,M0)  
        #call averts(p1,M0)
        #call TadpoleM0(s1m,p1,M0,0)
        .sort

************************************************************
#endprocedure



************************************************************
#procedure topt1
        #message this is topt1
        
        if( count(intt1,1));
        if( (count(x3,1)=0) && (count(x4,1)=1)
        && (count(x5,1)=1) && (count(x6,1)=1) )
        id x4*x5*x6 =  + (M^2*Gam(-1,1))^3;
        
        Multiply int0/intt1;
        endif;
        
#endprocedure        




************************************************************        
#procedure topn1
        #message this is topn1
        
* integrals of type N1
        
        if(count(intn1,1));        
        if( (count(x3,1)=1) && (count(x4,1)=1)
        && (count(x5,1)=1) && (count(x6,1)=1) )
        id x3*x4*x5*x6 =M^4*miBN*int0/intn1;
        
        endif;
        
        #call ACCU()
        
#endprocedure        


#procedure tad3l
*
* {top3l;6;3;0;1; ;(p1:3,1)(p2:4,3)(p3:3,2)(p4:1,4)(p5:4,2)(p6:2,1)}
* with all possible mass distributions, i.e.
* 
* {top3l;6;3;0;1; ;(p1:3,1)(p2:4,3)(p3:3,2)(p4:1,4)(p5:4,2)(p6:2,1);
*                  111111;011111;101111;110111;111011;111101;111110;
*                  001111;010111;011011;011101;011110;100111;101011;
*                  101101;101110;110011;110101;110110;111001;111010;
*                  111100;000111;001011;001101;001110;010011;010101;
*                  010110;011001;011010;011100;100011;100101;100110;
*                  101001;101010;101100;110001;110010;110100;111000;
*                  000011;000101;000110;001001;001010;001100;010001;
*                  010010;010100;011000;100001;100010;100100;101000;
*                  110000;000001;000010;000100;001000;010000;100000;
*                  000000\}

*
* do partial fractioning for all lines
*

        #do i = 1, 6
                #call partfrac(p`i',s`i'm)
        #enddo
        .sort

*
* discard massless tadpoles
*

        if (count(s1m,1,s2m,1,s3m,1,s4m,1,s5m,1,s6m,1)==0) discard;
        .sort

*
* map onto MATAD master topologies; mapping computed with EXP
* there are 2^6 - 1 = 63 cases
* (we dont care about symmetry ... )
*

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,intd6;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p6, p2,p1, p3,p4, p4,-p3, p5,-p2, p6,p5,
        s2m,s1m, s3m,s4m, s4m,s3m, s5m,s2m, s6m,s5m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p6, p3,-p4, p4,-p3, p5,p5, p6,-p2,
        s1m,s1m, s3m,s4m, s4m,s3m, s5m,s5m, s6m,s2m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p6, p4,p2, p5,-p5, p6,p3, 
        s1m,s1m, s2m,s4m, s4m,s2m, s5m,s5m, s6m,s3m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p2, p4,-p6, p5,-p5, p6,-p4, 
        s1m,s1m, s2m,s3m, s3m,s2m, s5m,s5m, s6m,s4m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p2, p2,p1, p3,-p3, p4,p4, p5,p6, p6,p5, 
        s1m,s2m, s2m,s1m, s3m,s3m, s4m,s4m, s6m,s5m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p2, p3,p3, p4,p4, p5,p5, p6,p6, 
        s1m,s1m, s2m,s2m, s3m,s3m, s4m,s4m, s5m,s5m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p3, p2,-p5, p3,-p6, p4,p2, p5,p4, p6,p1, 
        s3m,s6m, s4m,s2m, s5m,s4m, s6m,s1m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,-p6, p3,-p5, p4,p1, p5,p4, p6,p2, 
        s2m,s6m, s4m,s1m, s5m,s4m, s6m,s2m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p3, p2,p2, p3,p1, p4,-p5, p5,-p4, p6,-p6, 
        s2m,s2m, s3m,s1m, s5m,s4m, s6m,s6m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s2m,s3m, s3m,s6m, s4m,s5m, s6m,s4m);
        multiply,intbn;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p2, p4,-p6, p5,-p4, p6,-p5, 
        s2m,s1m, s3m,s2m, s4m,s6m, s5m,s4m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p6, p2,-p3, p3,p5, p4,p1, p5,p2, p6,p4, 
        s1m,s6m, s4m,s1m, s5m,s2m, s6m,s4m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p2, p2,p3, p3,-p1, p4,-p5, p5,-p6, p6,-p4, 
        s1m,s2m, s3m,s1m, s5m,s6m, s6m,s4m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p2, p4,-p6, p5,-p5, p6,-p4, 
        s1m,s1m, s3m,s2m, s4m,s6m, s6m,s4m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p6, p4,-p5, p5,-p4, p6,-p2, 
        s1m,s3m, s3m,s6m, s4m,s5m, s5m,s4m);
        multiply,intbn;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,-p5, p3,-p2, p4,p1, p5,p4, p6,p6, 
        s1m,s3m, s2m,s5m, s5m,s4m, s6m,s6m);
        multiply,intbn;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p2, p3,p3, p4,p4, p5,p5, p6,p6, 
        s1m,s1m, s2m,s2m, s4m,s4m, s6m,s6m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p2, p2,p1, p3,-p3, p4,p4, p5,p6, p6,p5, 
        s1m,s2m, s2m,s1m, s4m,s4m, s5m,s6m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p6, p3,-p4, p4,-p3, p5,p5, p6,-p2, 
        s1m,s1m, s2m,s6m, s3m,s4m, s6m,s2m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p6, p2,p1, p3,p4, p4,-p3, p5,-p2, p6,p5, 
        s1m,s6m, s2m,s1m, s3m,s4m, s5m,s2m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p6, p4,p2, p5,-p5, p6,p3, 
        s1m,s1m, s2m,s4m, s3m,s6m, s4m,s2m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p2, p3,p3, p4,p4, p5,p5, p6,p6, 
        s4m,s4m, s5m,s5m, s6m,s6m);
        multiply,intbm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p4, p2,p5, p3,-p2, p4,p6, p5,p3, p6,-p1, 
        s3m,s2m, s5m,s3m, s6m,s1m);
        multiply,intdm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p5, p3,-p4, p4,-p3, p5,p2, p6,-p6, 
        s3m,s4m, s4m,s3m, s6m,s6m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p5, p2,p1, p3,p4, p4,-p3, p5,-p6, p6,p2, 
        s3m,s4m, s4m,s3m, s5m,s6m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p5, p2,-p3, p3,p2, p4,p1, p5,p6, p6,p4, 
        s2m,s3m, s5m,s6m, s6m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p5, p4,p6, p5,-p2, p6,p3, 
        s2m,s4m, s4m,s6m, s6m,s3m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p1, p3,-p6, p4,p2, p5,p3, p6,-p5, 
        s2m,s1m, s4m,s2m, s5m,s3m);
        multiply,intdm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s2m,s3m, s3m,s6m, s6m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p6, p4,p2, p5,-p5, p6,p3, 
        s2m,s4m, s3m,s6m, s5m,s5m);
        multiply,intbm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p6, p3,p3, p4,p4, p5,p2, p6,p5, 
        s2m,s6m, s3m,s3m, s4m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,-p5, p3,-p2, p4,p1, p5,p4, p6,p6, 
        s1m,s3m, s5m,s4m, s6m,s6m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p6, p4,p2, p5,-p5, p6,p3,  
        s1m,s1m, s4m,s2m, s6m,s3m);
        multiply,intdm;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p1, p3,-p5, p4,p6, p5,p3, p6,-p2, 
        s1m,s4m, s4m,s6m, s5m,s3m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p4, p2,p1, p3,-p6, p4,p2, p5,p3, p6,-p5, 
        s1m,s4m, s3m,s6m, s6m,s5m);
        multiply,intbm;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p6, p4,-p5, p5,-p4, p6,-p2, 
        s1m,s3m, s3m,s6m, s5m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p6, p2,p1, p3,-p3, p4,p4, p5,p5, p6,p2, 
        s1m,s6m, s3m,s3m, s4m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p6, p2,p4, p3,p2, p4,p1, p5,-p5, p6,-p3, 
        s1m,s6m, s2m,s4m, s6m,s3m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p6, p3,-p2, p4,p1, p5,-p3, p6,-p5, 
        s1m,s4m, s2m,s6m, s5m,s3m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p5, p3,-p2, p4,p6, p5,p3, p6,-p1, 
        s1m,s4m, s2m,s5m, s4m,s6m);
        multiply,intbm;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p2, p3,p3, p4,p4, p5,p5, p6,p6, 
        s1m,s1m, s2m,s2m, s3m,s3m);
        multiply,intdm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,-p5, p3,-p2, p4,p1, p5,p4, p6,p6, 
        s5m,s4m, s6m,s6m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p2, p2,-p3, p3,-p5, p4,-p6, p5,-p1, p6,-p4, 
        s4m,s6m, s6m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p2, p3,p5, p4,-p6, p5,-p4, p6,-p1, 
        s4m,s6m, s5m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s3m,s6m, s6m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p6, p4,-p5, p5,-p4, p6,-p2, 
        s3m,s6m, s5m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,-p5, p3,-p2, p4,p1, p5,p4, p6,p6, 
        s3m,s2m, s4m,s1m);
        multiply,intdn;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p6, p4,-p5, p5,-p4, p6,-p2, 
        s2m,s1m, s6m,s2m);
        multiply,intdn;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,-p6, p3,-p1, p4,p2, p5,p4, p6,p5, 
        s2m,s6m, s5m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p5, p4,p6, p5,-p2, p6,p3, 
        s2m,s4m, s4m,s6m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p2, p2,-p6, p3,-p4, p4,-p3, p5,p1, p6,-p5, 
        s2m,s6m, s3m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p6, p2,-p3, p3,p1, p4,p2, p5,p5, p6,p4, 
        s1m,s6m, s6m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s1m,s1m, s5m,s2m);
        multiply,intdn;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p1, p3,-p5, p4,p6, p5,p3, p6,-p2, 
        s1m,s4m, s4m,s6m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p6, p2,p2, p3,p4, p4,-p3, p5,-p5, p6,p1, 
        s1m,s6m, s3m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p6, p3,-p2, p4,p1, p5,-p3, p6,-p5, 
        s1m,s4m, s2m,s6m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p5, p3,-p4, p4,-p3, p5,p2, p6,-p6, 
        s6m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p5, p2,p1, p3,p4, p4,-p3, p5,-p6, p6,p2, 
        s5m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p5, p4,p6, p5,-p2, p6,p3, 
        s4m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s3m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p6, p3,p3, p4,p4, p5,p2, p6,p5, 
        s2m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p6, p2,p1, p3,-p3, p4,p4, p5,p5, p6,p2, 
        s1m,s6m);
        multiply,intbn3;
        endif;
        .sort

*
* apply recurrence relations
* 

        #call dorec3l
        

        #call Conv2exact

        #call DoG
        #call subSimple
        #call GammaArgToOne
        
#endprocedure



#procedure matad(LOOPS)

        #message "                                                     "
        #message "                                                     "
        #message "      _____ _____ _____ _____ ____                   "
        #message "     |     |  _  |_   _|  _  |    \ ___ ___ ___      "
        #message "     | | | |     | | | |     |  |  |___|   | . |     "
        #message "     |_|_|_|__|__| |_| |__|__|____/    |_|_|_  |     "
        #message "                                           |___|     "
        #message "                                                     "
        #message "                                                     "
        #message "     Originally written by M.Steinhauser             "
        #message "     About exact version contact to:                 "
        #message "     pikelner[at]theor.jinr.ru, Andrey Pikelner      "
        #message "                                                     "
        #message MATAD called for `LOOPS'-loop integral

        .sort
*       Switch on rational arithmetics        
        PolyRatFun rat;        


        #call tad`LOOPS'l
        

        if(count(int0,1));
        Multiply 1/int0;  
        
        else;
        exit "Not all integrals reduced";

        endif;        
#endprocedure


#procedure subSimple

        id gm3norm(x1?,x2?,x3?) = 
        Gam(-1,2+x1+x2+x3)*
        Gam(0,1+x1+x3)*
        Gam(0,1+x2+x3)*
        Gam(1,-1-x3)*
        iGam(1,x1 )*
        iGam(1,x2 )*
        iGam(0,2+x1+x2+2*x3)*
        iGam(2,-1);        
        
        id gm2norm(k2?,k4?) =
        Gam(1,-1-k4)*
        Gam(1,1+k2+k4)*
        iGam(2,-1)*
        iGam(1,k2);


        id GschemeConstants(x1?,x2?)= 
        Gam(1,1+x1+x2)*
        Gam(1,-1-x1)*
        Gam(1,-1-x2)*
        iGam(1,x1)*
        iGam(1,x2)*
        iGam(2,-2-x1-x2)*
        rat(1,2-d/2);

#endprocedure

#procedure GammaArgToOne
* 
* Reduce all Gam[a,b] to Gam[1,x] 
*         
        id Gam(x?,y?) * iGam(x?,y?) = 1;
        
        id Gam(x?pos_,0)  = fac_(x-1);
        id iGam(x?pos_,0) = 1/fac_(x-1);
        id Gam(0,y?)      = rat(1,(2-d/2)*y) * Gam(1,y);
        id iGam(0,y?)     = rat((2-d/2)*y,1) * iGam(1,y);
        
        .sort
        
        repeat;
                id Gam(x?pos_,y?)         = rat(x-1+(2-d/2)*y,1) * Gam(x-1,y);
                id Gam(1,y?) * iGam(1,y?) = 1;
        endrepeat;
        
        .sort
        
        repeat; 
                id iGam(x?pos_,y?)        = rat(1,x-1+(2-d/2)*y) * iGam(x-1,y);
                id Gam(1,y?) * iGam(1,y?) = 1;
        endrepeat;
        
        .sort
        
        repeat;
                id Gam(x?neg0_,y?)        = rat(1,x+(2-d/2)*y) * Gam(x+1,y);
                id Gam(1,y?) * iGam(1,y?) = 1;
        endrepeat;
        
        .sort

        
        repeat;
                id iGam(x?neg0_,y?)       = rat(x+(2-d/2)*y,1) * iGam(x+1,y);
                id Gam(1,y?) * iGam(1,y?) = 1;
        endrepeat;
        
        .sort:All Gam and iGam reduced;
#endprocedure


#procedure cutep(x)
multiply, 1/ep^('x');
id ep=0;
multiply, ep^('x');
#endprocedure



#procedure subvalues
        
* Expansions from original MATAD        
        repeat  id,once    iGam(1,y?) = 
        1 - ep^2*y^2*z2/2 + ep^3*y^3*z3/3 + ep^4*y^4*(z2^2 - 2*z4)/8 +
        ep^5*y^5*(-5*z2*z3 + 6*z5)/30 +
        ep^6*y^6*(-3*z2^3 + 8*z3^2 + 18*z2*z4 - 24*z6)/144 +
        ep^7*Oep(7,iGamtr);
        
        
        repeat  id,once     Gam(1,y?) = 
        1 + ep^2*y^2*z2/2 - ep^3*y^3*z3/3 + ep^4*y^4*(z2^2 + 2*z4)/8 
        - ep^5*y^5*(5*z2*z3 + 6*z5)/30 +
        ep^6*y^6*(3*z2^3 + 8*z3^2 + 18*z2*z4 + 24*z6)/144 + 
        ep^7*Oep(7,Gamtr);
        


* Two-loop fully massive tadpole        
* MMM(1,1,1)/M^2        
        id miT111 = 
        -21/2 - 3/(2*ep^2) - 9/(2*ep) + (27*S2)/2 - (3*z2)/2 
        + ep * T111ep + ep^2 * T111ep2 + ep^3*Oep(3,T111tr);        
        
* MMMMMM(1,1,1,1,1,1)        
        id miD6 = 2*z3/ep + D6 + ep*D6ep + ep^2*D6ep2 + ep^3*Oep(3,D6tr);
        
* MMMMM0(1,1,1,1,1,1)        
        id miD5 = 2*z3/ep + D5 + ep*D5ep + ep^2*D5ep2 + ep^3*Oep(3,D5tr);
        
* 00MMMM(1,1,1,1,1,1)        
        id miD4 = 2*z3/ep + D4 + ep*D4ep + ep^2*D4ep2 + ep^3*Oep(3,D4tr);
        
* M000MM(1,1,1,1,1,1)
        id miD3 = 2*z3/ep + D3 + ep*D3ep + ep^2*D3ep2 + ep^3*Oep(3,D3tr);

* MM0000(1,1,1,1,1,1)        
        id miDN = 2*z3/ep + DN + ep*DNep + ep^2*DNep2 + ep^3*Oep(3,DNtr);
        
* MMM000(1,1,1,1,1,1)        
        id miDM = 2*z3/ep + DM + ep*DMep + ep^2*DMep2 + ep^3*Oep(3,DMtr);

* MMM000(1,1,1,1,1,0)/M^2        
        id miE3 = -2/3/ep^3-11/3/ep^2 + (-14 + (27*S2)/2 - 2*z2)/ep + E3 + ep*E3ep + ep^2*E3ep2 + ep^3*E3ep3 + ep^4*Oep(4,E3tr);

* MM00MM(1,1,0,0,1,1)/M^4                
        id miBN = agam/ep^3+bgam/ep^2+cgam/ep+dgam+egam*ep
        +fgam*ep^2 + BNep3*ep^3 + BNep4*ep^4 + ep^5*Oep(5,BNtr);

        id agam=2;
        id bgam=23/3;
        id cgam=35/2+3*z2;
        id dgam=275/192*16+23/2*z2-2*z3;
        id egam=-16*(189/128 - 89/48*z3 - 3/32*z4 - 105/64*z2 - 9/64*z2^2);
        id fgam = -384*(
        14917/18432 + 1/128*z3*z2 - 175/256*z3 + 649/1536*z4 + 1/320*z5
        - 275/3072*z2 - 23/1024*z2^2
        )
        +16*B4;
        
* M000MM(1,1,0,0,1,1)
        id miBN1 =
        + ep^-3 + ep^-2 * ( 15/4 ) + ep^-1 * ( 65/8 + 12/8*z2 )
        + 81/4*S2 - z3 + 135/16 + 90/16 *z2
        + ep*OepS2
        + ep^2*Oep2S2
        + ep^3*Oep(3,BN1tr);

        

#endprocedure


*--#[ exp4d :
*
#procedure exp4d(maxeppow)
* Expansion near d=4-2*ep      
        Multiply replace_(d,4-2*ep);
*
*	Expands the PolyRatFun to sufficient powers in ep.
*
        .sort:expansion-1;
        S DUMMYSYMBOL;        
        PolyRatFun;
* 
* We introduce DUMMYSYMBOL for proper terms ordering in SplitArg
* First term has smallest power of ep        
*         
        id	den(x?) = den(DUMMYSYMBOL*x);
        id	rat(x1?,x2?) = num(x1)*den(DUMMYSYMBOL*x2);
        
        SplitArg,den;
        Multiply replace_(DUMMYSYMBOL,1);
        
        repeat id den(x1?,x2?,x3?,?a) = den(x1,x2+x3,?a);
        id	den(x1?,x2?) = den(1,x2/x1)/x1;
        id	den(x1?) = 1/x1;
        
        .sort:expansion-2;
        id	num(x1?) = x1;
*
*       In masters we have minimal power ep^-3
*       And need to expand rat(ep) up to maxpow+3
*         
        if ( count(ep,1) > {`maxeppow' + 3} ) discard;
        repeat;
	        id den(1,x?) = 1-x*den(1,x);
	        if ( count(ep,1) > {`maxeppow'+ 3} ) discard;
        endrepeat;
        .sort:expansion-3;
        Symbol ep(:{`maxeppow'+3});
*         
        #call subvalues

        if ( count(ep,1) > `maxeppow' ) discard;
        .sort
* Switch back to ep without truncation        
        S ep;        
#endprocedure
*
*--#] exp4d : 


* Make step from d to d-2
#procedure dimstepdown

id miD6(d0?{,>4}) = Gam(4 + ep - 1/2*d0)^3*rat(-96,d0^9 - 18*d0^8*ep - 41*d0^8
       + 144*d0^7*ep^2 + 656*d0^7*ep + 739*d0^7 - 672*d0^6*ep^3 - 4592*d0^6*
      ep^2 - 10346*d0^6*ep - 7679*d0^6 + 2016*d0^5*ep^4 + 18368*d0^5*ep^3 + 
      62076*d0^5*ep^2 + 92148*d0^5*ep + 50644*d0^5 - 4032*d0^4*ep^5 - 45920*
      d0^4*ep^4 - 206920*d0^4*ep^3 - 460740*d0^4*ep^2 - 506440*d0^4*ep - 
      219584*d0^4 + 5376*d0^3*ep^6 + 73472*d0^3*ep^5 + 413840*d0^3*ep^4 + 
      1228640*d0^3*ep^3 + 2025760*d0^3*ep^2 + 1756672*d0^3*ep + 625056*d0^3 - 
      4608*d0^2*ep^7 - 73472*d0^2*ep^6 - 496608*d0^2*ep^5 - 1842960*d0^2*ep^4
       - 4051520*d0^2*ep^3 - 5270016*d0^2*ep^2 - 3750336*d0^2*ep - 1124496*
      d0^2 + 2304*d0*ep^8 + 41984*d0*ep^7 + 331072*d0*ep^6 + 1474368*d0*ep^5
       + 4051520*d0*ep^4 + 7026688*d0*ep^3 + 7500672*d0*ep^2 + 4497984*d0*ep
       + 1157760*d0 - 512*ep^9 - 10496*ep^8 - 94592*ep^7 - 491456*ep^6 - 
      1620608*ep^5 - 3513344*ep^4 - 5000448*ep^3 - 4497984*ep^2 - 2315520*ep
       - 518400) + Gam(4 + ep - 1/2*d0)*miT111( - 2 + d0)*rat(-24,d0^5 - 10*d0^4
      *ep - 20*d0^4 + 40*d0^3*ep^2 + 160*d0^3*ep + 155*d0^3 - 80*d0^2*ep^3 - 
      480*d0^2*ep^2 - 930*d0^2*ep - 580*d0^2 + 80*d0*ep^4 + 640*d0*ep^3 + 1860
      *d0*ep^2 + 2320*d0*ep + 1044*d0 - 32*ep^5 - 320*ep^4 - 1240*ep^3 - 2320*
      ep^2 - 2088*ep - 720) + miD6( - 2 + d0)*rat(-4,d0^3 - 6*d0^2*ep - 9*d0^2
       + 12*d0*ep^2 + 36*d0*ep + 26*d0 - 8*ep^3 - 36*ep^2 - 52*ep - 24) + 
      miD5( - 2 + d0)*rat( - 9*d0 + 18*ep + 54,d0^4 - 8*d0^3*ep - 14*d0^3 + 24
      *d0^2*ep^2 + 84*d0^2*ep + 71*d0^2 - 32*d0*ep^3 - 168*d0*ep^2 - 284*d0*ep
       - 154*d0 + 16*ep^4 + 112*ep^3 + 284*ep^2 + 308*ep + 120) + miBN( - 2 + 
      d0)*rat(1,4*d0^3 - 24*d0^2*ep - 48*d0^2 + 48*d0*ep^2 + 192*d0*ep + 188*
      d0 - 32*ep^3 - 192*ep^2 - 376*ep - 240);

id miD5(d0?{,>4}) = Gam(4 + ep - 1/2*d0)^3*rat(-48,d0^9 - 18*d0^8*ep - 41*d0^8
       + 144*d0^7*ep^2 + 656*d0^7*ep + 739*d0^7 - 672*d0^6*ep^3 - 4592*d0^6*
      ep^2 - 10346*d0^6*ep - 7679*d0^6 + 2016*d0^5*ep^4 + 18368*d0^5*ep^3 + 
      62076*d0^5*ep^2 + 92148*d0^5*ep + 50644*d0^5 - 4032*d0^4*ep^5 - 45920*
      d0^4*ep^4 - 206920*d0^4*ep^3 - 460740*d0^4*ep^2 - 506440*d0^4*ep - 
      219584*d0^4 + 5376*d0^3*ep^6 + 73472*d0^3*ep^5 + 413840*d0^3*ep^4 + 
      1228640*d0^3*ep^3 + 2025760*d0^3*ep^2 + 1756672*d0^3*ep + 625056*d0^3 - 
      4608*d0^2*ep^7 - 73472*d0^2*ep^6 - 496608*d0^2*ep^5 - 1842960*d0^2*ep^4
       - 4051520*d0^2*ep^3 - 5270016*d0^2*ep^2 - 3750336*d0^2*ep - 1124496*
      d0^2 + 2304*d0*ep^8 + 41984*d0*ep^7 + 331072*d0*ep^6 + 1474368*d0*ep^5
       + 4051520*d0*ep^4 + 7026688*d0*ep^3 + 7500672*d0*ep^2 + 4497984*d0*ep
       + 1157760*d0 - 512*ep^9 - 10496*ep^8 - 94592*ep^7 - 491456*ep^6 - 
      1620608*ep^5 - 3513344*ep^4 - 5000448*ep^3 - 4497984*ep^2 - 2315520*ep
       - 518400) + Gam(4 + ep - 1/2*d0)*miT111( - 2 + d0)*rat(-12,d0^5 - 10*d0^4
      *ep - 20*d0^4 + 40*d0^3*ep^2 + 160*d0^3*ep + 155*d0^3 - 80*d0^2*ep^3 - 
      480*d0^2*ep^2 - 930*d0^2*ep - 580*d0^2 + 80*d0*ep^4 + 640*d0*ep^3 + 1860
      *d0*ep^2 + 2320*d0*ep + 1044*d0 - 32*ep^5 - 320*ep^4 - 1240*ep^3 - 2320*
      ep^2 - 2088*ep - 720) + miD5( - 2 + d0)*rat( - 9*d0 + 18*ep + 54,2*d0^4
       - 16*d0^3*ep - 28*d0^3 + 48*d0^2*ep^2 + 168*d0^2*ep + 142*d0^2 - 64*d0*
      ep^3 - 336*d0*ep^2 - 568*d0*ep - 308*d0 + 32*ep^4 + 224*ep^3 + 568*ep^2
       + 616*ep + 240) + miBN( - 2 + d0)*rat( - 13*d0 + 26*ep + 74,24*d0^4 - 
      192*d0^3*ep - 336*d0^3 + 576*d0^2*ep^2 + 2016*d0^2*ep + 1704*d0^2 - 768*
      d0*ep^3 - 4032*d0*ep^2 - 6816*d0*ep - 3696*d0 + 384*ep^4 + 2688*ep^3 + 
      6816*ep^2 + 7392*ep + 2880);

id miD4(d0?{,>4}) = Gam( - 2 - ep + 1/2*d0)*Gam(4 + ep - 1/2*d0)*Gam(7 + 2*ep
       - d0)^2*Gam(10 + 3*ep - 3/2*d0)*iGam(13 + 4*ep - 2*d0)*rat(1408*d0^3 - 
      8448*d0^2*ep - 20544*d0^2 + 16896*d0*ep^2 + 82176*d0*ep + 99584*d0 - 
      11264*ep^3 - 82176*ep^2 - 199168*ep - 160512,81*d0^12 - 1944*d0^11*ep - 
      4455*d0^11 + 21384*d0^10*ep^2 + 98010*d0^10*ep + 111609*d0^10 - 142560*
      d0^9*ep^3 - 980100*d0^9*ep^2 - 2232180*d0^9*ep - 1683585*d0^9 + 641520*
      d0^8*ep^4 + 5880600*d0^8*ep^3 + 20089620*d0^8*ep^2 + 30304530*d0^8*ep + 
      17024958*d0^8 - 2052864*d0^7*ep^5 - 23522400*d0^7*ep^4 - 107144640*d0^7*
      ep^3 - 242436240*d0^7*ep^2 - 272399328*d0^7*ep - 121536720*d0^7 + 
      4790016*d0^6*ep^6 + 65862720*d0^6*ep^5 + 375006240*d0^6*ep^4 + 
      1131369120*d0^6*ep^3 + 1906795296*d0^6*ep^2 + 1701514080*d0^6*ep + 
      627746112*d0^6 - 8211456*d0^5*ep^7 - 131725440*d0^5*ep^6 - 900014976*
      d0^5*ep^5 - 3394107360*d0^5*ep^4 - 7627181184*d0^5*ep^3 - 10209084480*
      d0^5*ep^2 - 7532953344*d0^5*ep - 2362424400*d0^5 + 10264320*d0^4*ep^8 + 
      188179200*d0^4*ep^7 + 1500024960*d0^4*ep^6 + 6788214720*d0^4*ep^5 + 
      19067952960*d0^4*ep^4 + 34030281600*d0^4*ep^3 + 37664766720*d0^4*ep^2 + 
      23624244000*d0^4*ep + 6424976736*d0^4 - 9123840*d0^3*ep^9 - 188179200*
      d0^3*ep^8 - 1714314240*d0^3*ep^7 - 9050952960*d0^3*ep^6 - 30508724736*
      d0^3*ep^5 - 68060563200*d0^3*ep^4 - 100439377920*d0^3*ep^3 - 94496976000
      *d0^3*ep^2 - 51399813888*d0^3*ep - 12305692800*d0^3 + 5474304*d0^2*ep^10
       + 125452800*d0^2*ep^9 + 1285735680*d0^2*ep^8 + 7757959680*d0^2*ep^7 + 
      30508724736*d0^2*ep^6 + 81672675840*d0^2*ep^5 + 150659066880*d0^2*ep^4
       + 188993952000*d0^2*ep^3 + 154199441664*d0^2*ep^2 + 73834156800*d0^2*ep
       + 15740987904*d0^2 - 1990656*d0*ep^11 - 50181120*d0*ep^10 - 571438080*
      d0*ep^9 - 3878979840*d0*ep^8 - 17433556992*d0*ep^7 - 54448450560*d0*ep^6
       - 120527253504*d0*ep^5 - 188993952000*d0*ep^4 - 205599255552*d0*ep^3 - 
      147668313600*d0*ep^2 - 62963951616*d0*ep - 12060887040*d0 + 331776*ep^12
       + 9123840*ep^11 + 114287616*ep^10 + 861995520*ep^9 + 4358389248*ep^8 + 
      15556700160*ep^7 + 40175751168*ep^6 + 75597580800*ep^5 + 102799627776*
      ep^4 + 98445542400*ep^3 + 62963951616*ep^2 + 24121774080*ep + 4180377600
      ) + Gam(4 + ep - 1/2*d0)^3*rat( - 32*d0^2 + 128*d0*ep + 128*d0 - 128*
      ep^2 - 256*ep + 64,3*d0^11 - 66*d0^10*ep - 147*d0^10 + 660*d0^9*ep^2 + 
      2940*d0^9*ep + 3249*d0^9 - 3960*d0^8*ep^3 - 26460*d0^8*ep^2 - 58482*d0^8
      *ep - 42741*d0^8 + 15840*d0^7*ep^4 + 141120*d0^7*ep^3 + 467856*d0^7*ep^2
       + 683856*d0^7*ep + 371700*d0^7 - 44352*d0^6*ep^5 - 493920*d0^6*ep^4 - 
      2183328*d0^6*ep^3 - 4786992*d0^6*ep^2 - 5203800*d0^6*ep - 2242800*d0^6
       + 88704*d0^5*ep^6 + 1185408*d0^5*ep^5 + 6549984*d0^5*ep^4 + 19147968*
      d0^5*ep^3 + 31222800*d0^5*ep^2 + 26913600*d0^5*ep + 9576096*d0^5 - 
      126720*d0^4*ep^7 - 1975680*d0^4*ep^6 - 13099968*d0^4*ep^5 - 47869920*
      d0^4*ep^4 - 104076000*d0^4*ep^3 - 134568000*d0^4*ep^2 - 95760960*d0^4*ep
       - 28914864*d0^4 + 126720*d0^3*ep^8 + 2257920*d0^3*ep^7 + 17466624*d0^3*
      ep^6 + 76591872*d0^3*ep^5 + 208152000*d0^3*ep^4 + 358848000*d0^3*ep^3 + 
      383043840*d0^3*ep^2 + 231318912*d0^3*ep + 60463872*d0^3 - 84480*d0^2*
      ep^9 - 1693440*d0^2*ep^8 - 14971392*d0^2*ep^7 - 76591872*d0^2*ep^6 - 
      249782400*d0^2*ep^5 - 538272000*d0^2*ep^4 - 766087680*d0^2*ep^3 - 
      693956736*d0^2*ep^2 - 362783232*d0^2*ep - 83317248*d0^2 + 33792*d0*ep^10
       + 752640*d0*ep^9 + 7485696*d0*ep^8 + 43766784*d0*ep^7 + 166521600*d0*
      ep^6 + 430617600*d0*ep^5 + 766087680*d0*ep^4 + 925275648*d0*ep^3 + 
      725566464*d0*ep^2 + 333268992*d0*ep + 68014080*d0 - 6144*ep^11 - 150528*
      ep^10 - 1663488*ep^9 - 10941696*ep^8 - 47577600*ep^7 - 143539200*ep^6 - 
      306435072*ep^5 - 462637824*ep^4 - 483710976*ep^3 - 333268992*ep^2 - 
      136028160*ep - 24883200) + Gam(4 + ep - 1/2*d0)*miT111( - 2 + d0)*rat(-8,
      d0^5 - 10*d0^4*ep - 20*d0^4 + 40*d0^3*ep^2 + 160*d0^3*ep + 155*d0^3 - 80
      *d0^2*ep^3 - 480*d0^2*ep^2 - 930*d0^2*ep - 580*d0^2 + 80*d0*ep^4 + 640*
      d0*ep^3 + 1860*d0*ep^2 + 2320*d0*ep + 1044*d0 - 32*ep^5 - 320*ep^4 - 
      1240*ep^3 - 2320*ep^2 - 2088*ep - 720) + miD4( - 2 + d0)*rat(2,d0^3 - 6*
      d0^2*ep - 9*d0^2 + 12*d0*ep^2 + 36*d0*ep + 26*d0 - 8*ep^3 - 36*ep^2 - 52
      *ep - 24) + miE3( - 2 + d0)*rat(-1,d0^3 - 6*d0^2*ep - 9*d0^2 + 12*d0*
      ep^2 + 36*d0*ep + 26*d0 - 8*ep^3 - 36*ep^2 - 52*ep - 24) + miBN1( - 2
       + d0)*rat(2,3*d0^4 - 24*d0^3*ep - 42*d0^3 + 72*d0^2*ep^2 + 252*d0^2*ep
       + 213*d0^2 - 96*d0*ep^3 - 504*d0*ep^2 - 852*d0*ep - 462*d0 + 48*ep^4 + 
      336*ep^3 + 852*ep^2 + 924*ep + 360);

id miDN(d0?{,>4}) = Gam( - 2 - ep + 1/2*d0)^3*Gam(4 + ep - 1/2*d0)^2*Gam( - 8
       - 3*ep + 3/2*d0)*Gam(10 + 3*ep - 3/2*d0)*iGam( - 5 - 2*ep + d0)^2*rat(
      -64,d0^9 - 18*d0^8*ep - 41*d0^8 + 144*d0^7*ep^2 + 656*d0^7*ep + 739*d0^7
       - 672*d0^6*ep^3 - 4592*d0^6*ep^2 - 10346*d0^6*ep - 7679*d0^6 + 2016*
      d0^5*ep^4 + 18368*d0^5*ep^3 + 62076*d0^5*ep^2 + 92148*d0^5*ep + 50644*
      d0^5 - 4032*d0^4*ep^5 - 45920*d0^4*ep^4 - 206920*d0^4*ep^3 - 460740*d0^4
      *ep^2 - 506440*d0^4*ep - 219584*d0^4 + 5376*d0^3*ep^6 + 73472*d0^3*ep^5
       + 413840*d0^3*ep^4 + 1228640*d0^3*ep^3 + 2025760*d0^3*ep^2 + 1756672*
      d0^3*ep + 625056*d0^3 - 4608*d0^2*ep^7 - 73472*d0^2*ep^6 - 496608*d0^2*
      ep^5 - 1842960*d0^2*ep^4 - 4051520*d0^2*ep^3 - 5270016*d0^2*ep^2 - 
      3750336*d0^2*ep - 1124496*d0^2 + 2304*d0*ep^8 + 41984*d0*ep^7 + 331072*
      d0*ep^6 + 1474368*d0*ep^5 + 4051520*d0*ep^4 + 7026688*d0*ep^3 + 7500672*
      d0*ep^2 + 4497984*d0*ep + 1157760*d0 - 512*ep^9 - 10496*ep^8 - 94592*
      ep^7 - 491456*ep^6 - 1620608*ep^5 - 3513344*ep^4 - 5000448*ep^3 - 
      4497984*ep^2 - 2315520*ep - 518400) + Gam( - 2 - ep + 1/2*d0)*Gam(4 + ep
       - 1/2*d0)^2*Gam(7 + 2*ep - d0)*rat(128,d0^9 - 18*d0^8*ep - 40*d0^8 + 
      144*d0^7*ep^2 + 640*d0^7*ep + 703*d0^7 - 672*d0^6*ep^3 - 4480*d0^6*ep^2
       - 9842*d0^6*ep - 7120*d0^6 + 2016*d0^5*ep^4 + 17920*d0^5*ep^3 + 59052*
      d0^5*ep^2 + 85440*d0^5*ep + 45760*d0^5 - 4032*d0^4*ep^5 - 44800*d0^4*
      ep^4 - 196840*d0^4*ep^3 - 427200*d0^4*ep^2 - 457600*d0^4*ep - 193360*
      d0^4 + 5376*d0^3*ep^6 + 71680*d0^3*ep^5 + 393680*d0^3*ep^4 + 1139200*
      d0^3*ep^3 + 1830400*d0^3*ep^2 + 1546880*d0^3*ep + 536592*d0^3 - 4608*
      d0^2*ep^7 - 71680*d0^2*ep^6 - 472416*d0^2*ep^5 - 1708800*d0^2*ep^4 - 
      3660800*d0^2*ep^3 - 4640640*d0^2*ep^2 - 3219552*d0^2*ep - 941760*d0^2 + 
      2304*d0*ep^8 + 40960*d0*ep^7 + 314944*d0*ep^6 + 1367040*d0*ep^5 + 
      3660800*d0*ep^4 + 6187520*d0*ep^3 + 6439104*d0*ep^2 + 3767040*d0*ep + 
      946944*d0 - 512*ep^9 - 10240*ep^8 - 89984*ep^7 - 455680*ep^6 - 1464320*
      ep^5 - 3093760*ep^4 - 4292736*ep^3 - 3767040*ep^2 - 1893888*ep - 414720)
       + Gam( - 2 - ep + 1/2*d0)*Gam(4 + ep - 1/2*d0)*Gam(7 + 2*ep - d0)^2*
      Gam(10 + 3*ep - 3/2*d0)*iGam(13 + 4*ep - 2*d0)*rat( - 10240*d0^2 + 40960
      *d0*ep + 105472*d0 - 40960*ep^2 - 210944*ep - 270336,81*d0^11 - 1782*
      d0^10*ep - 4050*d0^10 + 17820*d0^9*ep^2 + 81000*d0^9*ep + 91359*d0^9 - 
      106920*d0^8*ep^3 - 729000*d0^8*ep^2 - 1644462*d0^8*ep - 1226790*d0^8 + 
      427680*d0^7*ep^4 + 3888000*d0^7*ep^3 + 13155696*d0^7*ep^2 + 19628640*
      d0^7*ep + 10891008*d0^7 - 1197504*d0^6*ep^5 - 13608000*d0^6*ep^4 - 
      61393248*d0^6*ep^3 - 137400480*d0^6*ep^2 - 152474112*d0^6*ep - 67081680*
      d0^6 + 2395008*d0^5*ep^6 + 32659200*d0^5*ep^5 + 184179744*d0^5*ep^4 + 
      549601920*d0^5*ep^3 + 914844672*d0^5*ep^2 + 804980160*d0^5*ep + 
      292337712*d0^5 - 3421440*d0^4*ep^7 - 54432000*d0^4*ep^6 - 368359488*d0^4
      *ep^5 - 1374004800*d0^4*ep^4 - 3049482240*d0^4*ep^3 - 4024900800*d0^4*
      ep^2 - 2923377120*d0^4*ep - 900735840*d0^4 + 3421440*d0^3*ep^8 + 
      62208000*d0^3*ep^7 + 491145984*d0^3*ep^6 + 2198407680*d0^3*ep^5 + 
      6098964480*d0^3*ep^4 + 10733068800*d0^3*ep^3 + 11693508480*d0^3*ep^2 + 
      7205886720*d0^3*ep + 1921297536*d0^3 - 2280960*d0^2*ep^9 - 46656000*d0^2
      *ep^8 - 420982272*d0^2*ep^7 - 2198407680*d0^2*ep^6 - 7318757376*d0^2*
      ep^5 - 16099603200*d0^2*ep^4 - 23387016960*d0^2*ep^3 - 21617660160*d0^2*
      ep^2 - 11527785216*d0^2*ep - 2699205120*d0^2 + 912384*d0*ep^10 + 
      20736000*d0*ep^9 + 210491136*d0*ep^8 + 1256232960*d0*ep^7 + 4879171584*
      d0*ep^6 + 12879682560*d0*ep^5 + 23387016960*d0*ep^4 + 28823546880*d0*
      ep^3 + 23055570432*d0*ep^2 + 10796820480*d0*ep + 2244962304*d0 - 165888*
      ep^11 - 4147200*ep^10 - 46775808*ep^9 - 314058240*ep^8 - 1394049024*ep^7
       - 4293227520*ep^6 - 9354806784*ep^5 - 14411773440*ep^4 - 15370380288*
      ep^3 - 10796820480*ep^2 - 4489924608*ep - 836075520) + miDN( - 2 + d0)*
      rat(4,d0^3 - 6*d0^2*ep - 9*d0^2 + 12*d0*ep^2 + 36*d0*ep + 26*d0 - 8*ep^3
       - 36*ep^2 - 52*ep - 24);

id miDM(d0?{,>4}) = Gam( - 2 - ep + 1/2*d0)*Gam(4 + ep - 1/2*d0)^2*Gam(7 + 2*
      ep - d0)*rat(-96,d0^9 - 18*d0^8*ep - 40*d0^8 + 144*d0^7*ep^2 + 640*d0^7*
      ep + 703*d0^7 - 672*d0^6*ep^3 - 4480*d0^6*ep^2 - 9842*d0^6*ep - 7120*
      d0^6 + 2016*d0^5*ep^4 + 17920*d0^5*ep^3 + 59052*d0^5*ep^2 + 85440*d0^5*
      ep + 45760*d0^5 - 4032*d0^4*ep^5 - 44800*d0^4*ep^4 - 196840*d0^4*ep^3 - 
      427200*d0^4*ep^2 - 457600*d0^4*ep - 193360*d0^4 + 5376*d0^3*ep^6 + 71680
      *d0^3*ep^5 + 393680*d0^3*ep^4 + 1139200*d0^3*ep^3 + 1830400*d0^3*ep^2 + 
      1546880*d0^3*ep + 536592*d0^3 - 4608*d0^2*ep^7 - 71680*d0^2*ep^6 - 
      472416*d0^2*ep^5 - 1708800*d0^2*ep^4 - 3660800*d0^2*ep^3 - 4640640*d0^2*
      ep^2 - 3219552*d0^2*ep - 941760*d0^2 + 2304*d0*ep^8 + 40960*d0*ep^7 + 
      314944*d0*ep^6 + 1367040*d0*ep^5 + 3660800*d0*ep^4 + 6187520*d0*ep^3 + 
      6439104*d0*ep^2 + 3767040*d0*ep + 946944*d0 - 512*ep^9 - 10240*ep^8 - 
      89984*ep^7 - 455680*ep^6 - 1464320*ep^5 - 3093760*ep^4 - 4292736*ep^3 - 
      3767040*ep^2 - 1893888*ep - 414720) + Gam( - 2 - ep + 1/2*d0)*Gam(4 + ep
       - 1/2*d0)*Gam(7 + 2*ep - d0)^2*Gam(10 + 3*ep - 3/2*d0)*iGam(13 + 4*ep
       - 2*d0)*rat(3328*d0^2 - 13312*d0*ep - 33664*d0 + 13312*ep^2 + 67328*ep
       + 84480,27*d0^11 - 594*d0^10*ep - 1350*d0^10 + 5940*d0^9*ep^2 + 27000*
      d0^9*ep + 30453*d0^9 - 35640*d0^8*ep^3 - 243000*d0^8*ep^2 - 548154*d0^8*
      ep - 408930*d0^8 + 142560*d0^7*ep^4 + 1296000*d0^7*ep^3 + 4385232*d0^7*
      ep^2 + 6542880*d0^7*ep + 3630336*d0^7 - 399168*d0^6*ep^5 - 4536000*d0^6*
      ep^4 - 20464416*d0^6*ep^3 - 45800160*d0^6*ep^2 - 50824704*d0^6*ep - 
      22360560*d0^6 + 798336*d0^5*ep^6 + 10886400*d0^5*ep^5 + 61393248*d0^5*
      ep^4 + 183200640*d0^5*ep^3 + 304948224*d0^5*ep^2 + 268326720*d0^5*ep + 
      97445904*d0^5 - 1140480*d0^4*ep^7 - 18144000*d0^4*ep^6 - 122786496*d0^4*
      ep^5 - 458001600*d0^4*ep^4 - 1016494080*d0^4*ep^3 - 1341633600*d0^4*ep^2
       - 974459040*d0^4*ep - 300245280*d0^4 + 1140480*d0^3*ep^8 + 20736000*
      d0^3*ep^7 + 163715328*d0^3*ep^6 + 732802560*d0^3*ep^5 + 2032988160*d0^3*
      ep^4 + 3577689600*d0^3*ep^3 + 3897836160*d0^3*ep^2 + 2401962240*d0^3*ep
       + 640432512*d0^3 - 760320*d0^2*ep^9 - 15552000*d0^2*ep^8 - 140327424*
      d0^2*ep^7 - 732802560*d0^2*ep^6 - 2439585792*d0^2*ep^5 - 5366534400*d0^2
      *ep^4 - 7795672320*d0^2*ep^3 - 7205886720*d0^2*ep^2 - 3842595072*d0^2*ep
       - 899735040*d0^2 + 304128*d0*ep^10 + 6912000*d0*ep^9 + 70163712*d0*ep^8
       + 418744320*d0*ep^7 + 1626390528*d0*ep^6 + 4293227520*d0*ep^5 + 
      7795672320*d0*ep^4 + 9607848960*d0*ep^3 + 7685190144*d0*ep^2 + 
      3598940160*d0*ep + 748320768*d0 - 55296*ep^11 - 1382400*ep^10 - 15591936
      *ep^9 - 104686080*ep^8 - 464683008*ep^7 - 1431075840*ep^6 - 3118268928*
      ep^5 - 4803924480*ep^4 - 5123460096*ep^3 - 3598940160*ep^2 - 1496641536*
      ep - 278691840) + miDM( - 2 + d0)*rat(2,d0^3 - 6*d0^2*ep - 9*d0^2 + 12*
      d0*ep^2 + 36*d0*ep + 26*d0 - 8*ep^3 - 36*ep^2 - 52*ep - 24) + miE3( - 2
       + d0)*rat(3,d0^3 - 6*d0^2*ep - 9*d0^2 + 12*d0*ep^2 + 36*d0*ep + 26*d0
       - 8*ep^3 - 36*ep^2 - 52*ep - 24);

id miE3(d0?{,>4}) = Gam( - 2 - ep + 1/2*d0)*Gam(4 + ep - 1/2*d0)^2*Gam(7 + 2*
      ep - d0)*rat(16,d0^9 - 18*d0^8*ep - 39*d0^8 + 144*d0^7*ep^2 + 624*d0^7*
      ep + 667*d0^7 - 672*d0^6*ep^3 - 4368*d0^6*ep^2 - 9338*d0^6*ep - 6561*
      d0^6 + 2016*d0^5*ep^4 + 17472*d0^5*ep^3 + 56028*d0^5*ep^2 + 78732*d0^5*
      ep + 40876*d0^5 - 4032*d0^4*ep^5 - 43680*d0^4*ep^4 - 186760*d0^4*ep^3 - 
      393660*d0^4*ep^2 - 408760*d0^4*ep - 167136*d0^4 + 5376*d0^3*ep^6 + 69888
      *d0^3*ep^5 + 373520*d0^3*ep^4 + 1049760*d0^3*ep^3 + 1635040*d0^3*ep^2 + 
      1337088*d0^3*ep + 448128*d0^3 - 4608*d0^2*ep^7 - 69888*d0^2*ep^6 - 
      448224*d0^2*ep^5 - 1574640*d0^2*ep^4 - 3270080*d0^2*ep^3 - 4011264*d0^2*
      ep^2 - 2688768*d0^2*ep - 759024*d0^2 + 2304*d0*ep^8 + 39936*d0*ep^7 + 
      298816*d0*ep^6 + 1259712*d0*ep^5 + 3270080*d0*ep^4 + 5348352*d0*ep^3 + 
      5377536*d0*ep^2 + 3036096*d0*ep + 736128*d0 - 512*ep^9 - 9984*ep^8 - 
      85376*ep^7 - 419904*ep^6 - 1308032*ep^5 - 2674176*ep^4 - 3585024*ep^3 - 
      3036096*ep^2 - 1472256*ep - 311040) + Gam( - 2 - ep + 1/2*d0)*Gam(4 + ep
       - 1/2*d0)*Gam(7 + 2*ep - d0)^2*Gam(10 + 3*ep - 3/2*d0)*iGam(13 + 4*ep
       - 2*d0)*rat(2176*d0^3 - 13056*d0^2*ep - 27840*d0^2 + 26112*d0*ep^2 + 
      111360*d0*ep + 114944*d0 - 17408*ep^3 - 111360*ep^2 - 229888*ep - 152064
      ,81*d0^12 - 1944*d0^11*ep - 4293*d0^11 + 21384*d0^10*ep^2 + 94446*d0^10*
      ep + 103509*d0^10 - 142560*d0^9*ep^3 - 944460*d0^9*ep^2 - 2070180*d0^9*
      ep - 1500867*d0^9 + 641520*d0^8*ep^4 + 5666760*d0^8*ep^3 + 18631620*d0^8
      *ep^2 + 27015606*d0^8*ep + 14571378*d0^8 - 2052864*d0^7*ep^5 - 22667040*
      d0^7*ep^4 - 99368640*d0^7*ep^3 - 216124848*d0^7*ep^2 - 233142048*d0^7*ep
       - 99754704*d0^7 + 4790016*d0^6*ep^6 + 63467712*d0^6*ep^5 + 347790240*
      d0^6*ep^4 + 1008582624*d0^6*ep^3 + 1631994336*d0^6*ep^2 + 1396565856*
      d0^6*ep + 493582752*d0^6 - 8211456*d0^5*ep^7 - 126935424*d0^5*ep^6 - 
      834696576*d0^5*ep^5 - 3025747872*d0^5*ep^4 - 6527977344*d0^5*ep^3 - 
      8379395136*d0^5*ep^2 - 5922993024*d0^5*ep - 1777748976*d0^5 + 10264320*
      d0^4*ep^8 + 181336320*d0^4*ep^7 + 1391160960*d0^4*ep^6 + 6051495744*d0^4
      *ep^5 + 16319943360*d0^4*ep^4 + 27931317120*d0^4*ep^3 + 29614965120*d0^4
      *ep^2 + 17777489760*d0^4*ep + 4623505056*d0^4 - 9123840*d0^3*ep^9 - 
      181336320*d0^3*ep^8 - 1589898240*d0^3*ep^7 - 8068660992*d0^3*ep^6 - 
      26111909376*d0^3*ep^5 - 55862634240*d0^3*ep^4 - 78973240320*d0^3*ep^3 - 
      71109959040*d0^3*ep^2 - 36988040448*d0^3*ep - 8463097728*d0^3 + 5474304*
      d0^2*ep^10 + 120890880*d0^2*ep^9 + 1192423680*d0^2*ep^8 + 6915995136*
      d0^2*ep^7 + 26111909376*d0^2*ep^6 + 67035161088*d0^2*ep^5 + 118459860480
      *d0^2*ep^4 + 142219918080*d0^2*ep^3 + 110964121344*d0^2*ep^2 + 
      50778586368*d0^2*ep + 10342577664*d0^2 - 1990656*d0*ep^11 - 48356352*d0*
      ep^10 - 529966080*d0*ep^9 - 3457997568*d0*ep^8 - 14921091072*d0*ep^7 - 
      44690107392*d0*ep^6 - 94767888384*d0*ep^5 - 142219918080*d0*ep^4 - 
      147952161792*d0*ep^3 - 101557172736*d0*ep^2 - 41370310656*d0*ep - 
      7570962432*d0 + 331776*ep^12 + 8792064*ep^11 + 105993216*ep^10 + 
      768443904*ep^9 + 3730272768*ep^8 + 12768602112*ep^7 + 31589296128*ep^6
       + 56887967232*ep^5 + 73976080896*ep^4 + 67704781824*ep^3 + 41370310656*
      ep^2 + 15141924864*ep + 2508226560) + miE3( - 2 + d0)*rat(3,2*d0^3 - 12*
      d0^2*ep - 16*d0^2 + 24*d0*ep^2 + 64*d0*ep + 42*d0 - 16*ep^3 - 64*ep^2 - 
      84*ep - 36);

id miBN(d0?{,>4}) = Gam(4 + ep - 1/2*d0)^3*rat(11264*d0 - 22528*ep - 38912,27*
      d0^10 - 540*d0^9*ep - 1107*d0^9 + 4860*d0^8*ep^2 + 19926*d0^8*ep + 20166
      *d0^8 - 25920*d0^7*ep^3 - 159408*d0^7*ep^2 - 322656*d0^7*ep - 214896*
      d0^7 + 90720*d0^6*ep^4 + 743904*d0^6*ep^3 + 2258592*d0^6*ep^2 + 3008544*
      d0^6*ep + 1483200*d0^6 - 217728*d0^5*ep^5 - 2231712*d0^5*ep^4 - 9034368*
      d0^5*ep^3 - 18051264*d0^5*ep^2 - 17798400*d0^5*ep - 6926640*d0^5 + 
      362880*d0^4*ep^6 + 4463424*d0^4*ep^5 + 22585920*d0^4*ep^4 + 60170880*
      d0^4*ep^3 + 88992000*d0^4*ep^2 + 69266400*d0^4*ep + 22161504*d0^4 - 
      414720*d0^3*ep^7 - 5951232*d0^3*ep^6 - 36137472*d0^3*ep^5 - 120341760*
      d0^3*ep^4 - 237312000*d0^3*ep^3 - 277065600*d0^3*ep^2 - 177292032*d0^3*
      ep - 47954304*d0^3 + 311040*d0^2*ep^8 + 5101056*d0^2*ep^7 + 36137472*
      d0^2*ep^6 + 144410112*d0^2*ep^5 + 355968000*d0^2*ep^4 + 554131200*d0^2*
      ep^3 + 531876096*d0^2*ep^2 + 287725824*d0^2*ep + 67143168*d0^2 - 138240*
      d0*ep^9 - 2550528*d0*ep^8 - 20649984*d0*ep^7 - 96273408*d0*ep^6 - 
      284774400*d0*ep^5 - 554131200*d0*ep^4 - 709168128*d0*ep^3 - 575451648*d0
      *ep^2 - 268572672*d0*ep - 54908928*d0 + 27648*ep^10 + 566784*ep^9 + 
      5162496*ep^8 + 27506688*ep^7 + 94924800*ep^6 + 221652480*ep^5 + 
      354584064*ep^4 + 383634432*ep^3 + 268572672*ep^2 + 109817856*ep + 
      19906560) + miBN( - 2 + d0)*rat( - 128*d0 + 256*ep + 512,27*d0^4 - 216*
      d0^3*ep - 297*d0^3 + 648*d0^2*ep^2 + 1782*d0^2*ep + 1212*d0^2 - 864*d0*
      ep^3 - 3564*d0*ep^2 - 4848*d0*ep - 2172*d0 + 432*ep^4 + 2376*ep^3 + 4848
      *ep^2 + 4344*ep + 1440);

id miD3(d0?{,>4}) = Gam( - 2 - ep + 1/2*d0)^2*Gam(7 + 2*ep - d0)*Gam(10 + 
      3*ep - 3/2*d0)*rat(32,27*d0^9 - 486*d0^8*ep - 1107*d0^8 + 3888*d0^7*ep^2
       + 17712*d0^7*ep + 19977*d0^7 - 18144*d0^6*ep^3 - 123984*d0^6*ep^2 - 
      279678*d0^6*ep - 208077*d0^6 + 54432*d0^5*ep^4 + 495936*d0^5*ep^3 + 
      1678068*d0^5*ep^2 + 2496924*d0^5*ep + 1377108*d0^5 - 108864*d0^4*ep^5 - 
      1239840*d0^4*ep^4 - 5593560*d0^4*ep^3 - 12484620*d0^4*ep^2 - 13771080*
      d0^4*ep - 5998008*d0^4 + 145152*d0^3*ep^6 + 1983744*d0^3*ep^5 + 11187120
      *d0^3*ep^4 + 33292320*d0^3*ep^3 + 55084320*d0^3*ep^2 + 47984064*d0^3*ep
       + 17166288*d0^3 - 124416*d0^2*ep^7 - 1983744*d0^2*ep^6 - 13424544*d0^2*
      ep^5 - 49938480*d0^2*ep^4 - 110168640*d0^2*ep^3 - 143952192*d0^2*ep^2 - 
      102997728*d0^2*ep - 31071888*d0^2 + 62208*d0*ep^8 + 1133568*d0*ep^7 + 
      8949696*d0*ep^6 + 39950784*d0*ep^5 + 110168640*d0*ep^4 + 191936256*d0*
      ep^3 + 205995456*d0*ep^2 + 124287552*d0*ep + 32201280*d0 - 13824*ep^9 - 
      283392*ep^8 - 2557056*ep^7 - 13316928*ep^6 - 44067456*ep^5 - 95968128*
      ep^4 - 137330304*ep^3 - 124287552*ep^2 - 64402560*ep - 14515200) + Gam(
       - 2 - ep + 1/2*d0)*Gam(4 + ep - 1/2*d0)^2*Gam(7 + 2*ep - d0)*rat(-32,
      d0^10 - 20*d0^9*ep - 45*d0^9 + 180*d0^8*ep^2 + 810*d0^8*ep + 903*d0^8 - 
      960*d0^7*ep^3 - 6480*d0^7*ep^2 - 14448*d0^7*ep - 10635*d0^7 + 3360*d0^6*
      ep^4 + 30240*d0^6*ep^3 + 101136*d0^6*ep^2 + 148890*d0^6*ep + 81360*d0^6
       - 8064*d0^5*ep^5 - 90720*d0^5*ep^4 - 404544*d0^5*ep^3 - 893340*d0^5*
      ep^2 - 976320*d0^5*ep - 422160*d0^5 + 13440*d0^4*ep^6 + 181440*d0^4*ep^5
       + 1011360*d0^4*ep^4 + 2977800*d0^4*ep^3 + 4881600*d0^4*ep^2 + 4221600*
      d0^4*ep + 1503392*d0^4 - 15360*d0^3*ep^7 - 241920*d0^3*ep^6 - 1618176*
      d0^3*ep^5 - 5955600*d0^3*ep^4 - 13017600*d0^3*ep^3 - 16886400*d0^3*ep^2
       - 12027136*d0^3*ep - 3624720*d0^3 + 11520*d0^2*ep^8 + 207360*d0^2*ep^7
       + 1618176*d0^2*ep^6 + 7146720*d0^2*ep^5 + 19526400*d0^2*ep^4 + 33772800
      *d0^2*ep^3 + 36081408*d0^2*ep^2 + 21748320*d0^2*ep + 5655744*d0^2 - 5120
      *d0*ep^9 - 103680*d0*ep^8 - 924672*d0*ep^7 - 4764480*d0*ep^6 - 15621120*
      d0*ep^5 - 33772800*d0*ep^4 - 48108544*d0*ep^3 - 43496640*d0*ep^2 - 
      22622976*d0*ep - 5149440*d0 + 1024*ep^10 + 23040*ep^9 + 231168*ep^8 + 
      1361280*ep^7 + 5207040*ep^6 + 13509120*ep^5 + 24054272*ep^4 + 28997760*
      ep^3 + 22622976*ep^2 + 10298880*ep + 2073600) + Gam( - 2 - ep + 1/2*d0)*
      Gam(4 + ep - 1/2*d0)*Gam(7 + 2*ep - d0)^2*Gam(10 + 3*ep - 3/2*d0)*iGam(
      13 + 4*ep - 2*d0)*rat( - 4736*d0^3 + 28416*d0^2*ep + 70848*d0^2 - 56832*
      d0*ep^2 - 283392*d0*ep - 352384*d0 + 37888*ep^3 + 283392*ep^2 + 704768*
      ep + 582912,81*d0^12 - 1944*d0^11*ep - 4455*d0^11 + 21384*d0^10*ep^2 + 
      98010*d0^10*ep + 111609*d0^10 - 142560*d0^9*ep^3 - 980100*d0^9*ep^2 - 
      2232180*d0^9*ep - 1683585*d0^9 + 641520*d0^8*ep^4 + 5880600*d0^8*ep^3 + 
      20089620*d0^8*ep^2 + 30304530*d0^8*ep + 17024958*d0^8 - 2052864*d0^7*
      ep^5 - 23522400*d0^7*ep^4 - 107144640*d0^7*ep^3 - 242436240*d0^7*ep^2 - 
      272399328*d0^7*ep - 121536720*d0^7 + 4790016*d0^6*ep^6 + 65862720*d0^6*
      ep^5 + 375006240*d0^6*ep^4 + 1131369120*d0^6*ep^3 + 1906795296*d0^6*ep^2
       + 1701514080*d0^6*ep + 627746112*d0^6 - 8211456*d0^5*ep^7 - 131725440*
      d0^5*ep^6 - 900014976*d0^5*ep^5 - 3394107360*d0^5*ep^4 - 7627181184*d0^5
      *ep^3 - 10209084480*d0^5*ep^2 - 7532953344*d0^5*ep - 2362424400*d0^5 + 
      10264320*d0^4*ep^8 + 188179200*d0^4*ep^7 + 1500024960*d0^4*ep^6 + 
      6788214720*d0^4*ep^5 + 19067952960*d0^4*ep^4 + 34030281600*d0^4*ep^3 + 
      37664766720*d0^4*ep^2 + 23624244000*d0^4*ep + 6424976736*d0^4 - 9123840*
      d0^3*ep^9 - 188179200*d0^3*ep^8 - 1714314240*d0^3*ep^7 - 9050952960*d0^3
      *ep^6 - 30508724736*d0^3*ep^5 - 68060563200*d0^3*ep^4 - 100439377920*
      d0^3*ep^3 - 94496976000*d0^3*ep^2 - 51399813888*d0^3*ep - 12305692800*
      d0^3 + 5474304*d0^2*ep^10 + 125452800*d0^2*ep^9 + 1285735680*d0^2*ep^8
       + 7757959680*d0^2*ep^7 + 30508724736*d0^2*ep^6 + 81672675840*d0^2*ep^5
       + 150659066880*d0^2*ep^4 + 188993952000*d0^2*ep^3 + 154199441664*d0^2*
      ep^2 + 73834156800*d0^2*ep + 15740987904*d0^2 - 1990656*d0*ep^11 - 
      50181120*d0*ep^10 - 571438080*d0*ep^9 - 3878979840*d0*ep^8 - 17433556992
      *d0*ep^7 - 54448450560*d0*ep^6 - 120527253504*d0*ep^5 - 188993952000*d0*
      ep^4 - 205599255552*d0*ep^3 - 147668313600*d0*ep^2 - 62963951616*d0*ep
       - 12060887040*d0 + 331776*ep^12 + 9123840*ep^11 + 114287616*ep^10 + 
      861995520*ep^9 + 4358389248*ep^8 + 15556700160*ep^7 + 40175751168*ep^6
       + 75597580800*ep^5 + 102799627776*ep^4 + 98445542400*ep^3 + 62963951616
      *ep^2 + 24121774080*ep + 4180377600) + Gam(4 + ep - 1/2*d0)^3*rat(32*
      d0^2 - 128*d0*ep - 128*d0 + 128*ep^2 + 256*ep - 64,3*d0^11 - 66*d0^10*ep
       - 147*d0^10 + 660*d0^9*ep^2 + 2940*d0^9*ep + 3249*d0^9 - 3960*d0^8*ep^3
       - 26460*d0^8*ep^2 - 58482*d0^8*ep - 42741*d0^8 + 15840*d0^7*ep^4 + 
      141120*d0^7*ep^3 + 467856*d0^7*ep^2 + 683856*d0^7*ep + 371700*d0^7 - 
      44352*d0^6*ep^5 - 493920*d0^6*ep^4 - 2183328*d0^6*ep^3 - 4786992*d0^6*
      ep^2 - 5203800*d0^6*ep - 2242800*d0^6 + 88704*d0^5*ep^6 + 1185408*d0^5*
      ep^5 + 6549984*d0^5*ep^4 + 19147968*d0^5*ep^3 + 31222800*d0^5*ep^2 + 
      26913600*d0^5*ep + 9576096*d0^5 - 126720*d0^4*ep^7 - 1975680*d0^4*ep^6
       - 13099968*d0^4*ep^5 - 47869920*d0^4*ep^4 - 104076000*d0^4*ep^3 - 
      134568000*d0^4*ep^2 - 95760960*d0^4*ep - 28914864*d0^4 + 126720*d0^3*
      ep^8 + 2257920*d0^3*ep^7 + 17466624*d0^3*ep^6 + 76591872*d0^3*ep^5 + 
      208152000*d0^3*ep^4 + 358848000*d0^3*ep^3 + 383043840*d0^3*ep^2 + 
      231318912*d0^3*ep + 60463872*d0^3 - 84480*d0^2*ep^9 - 1693440*d0^2*ep^8
       - 14971392*d0^2*ep^7 - 76591872*d0^2*ep^6 - 249782400*d0^2*ep^5 - 
      538272000*d0^2*ep^4 - 766087680*d0^2*ep^3 - 693956736*d0^2*ep^2 - 
      362783232*d0^2*ep - 83317248*d0^2 + 33792*d0*ep^10 + 752640*d0*ep^9 + 
      7485696*d0*ep^8 + 43766784*d0*ep^7 + 166521600*d0*ep^6 + 430617600*d0*
      ep^5 + 766087680*d0*ep^4 + 925275648*d0*ep^3 + 725566464*d0*ep^2 + 
      333268992*d0*ep + 68014080*d0 - 6144*ep^11 - 150528*ep^10 - 1663488*ep^9
       - 10941696*ep^8 - 47577600*ep^7 - 143539200*ep^6 - 306435072*ep^5 - 
      462637824*ep^4 - 483710976*ep^3 - 333268992*ep^2 - 136028160*ep - 
      24883200) + miD3( - 2 + d0)*rat(2,d0^3 - 6*d0^2*ep - 9*d0^2 + 12*d0*
      ep^2 + 36*d0*ep + 26*d0 - 8*ep^3 - 36*ep^2 - 52*ep - 24) + miBN1( - 2
       + d0)*rat(-2,d0^3 - 6*d0^2*ep - 9*d0^2 + 12*d0*ep^2 + 36*d0*ep + 26*d0
       - 8*ep^3 - 36*ep^2 - 52*ep - 24);

id miBN1(d0?{,>4}) = Gam(4 + ep - 1/2*d0)^3*rat(1920*d0 - 3840*ep - 6656,9*
      d0^10 - 180*d0^9*ep - 369*d0^9 + 1620*d0^8*ep^2 + 6642*d0^8*ep + 6722*
      d0^8 - 8640*d0^7*ep^3 - 53136*d0^7*ep^2 - 107552*d0^7*ep - 71632*d0^7 + 
      30240*d0^6*ep^4 + 247968*d0^6*ep^3 + 752864*d0^6*ep^2 + 1002848*d0^6*ep
       + 494400*d0^6 - 72576*d0^5*ep^5 - 743904*d0^5*ep^4 - 3011456*d0^5*ep^3
       - 6017088*d0^5*ep^2 - 5932800*d0^5*ep - 2308880*d0^5 + 120960*d0^4*ep^6
       + 1487808*d0^4*ep^5 + 7528640*d0^4*ep^4 + 20056960*d0^4*ep^3 + 29664000
      *d0^4*ep^2 + 23088800*d0^4*ep + 7387168*d0^4 - 138240*d0^3*ep^7 - 
      1983744*d0^3*ep^6 - 12045824*d0^3*ep^5 - 40113920*d0^3*ep^4 - 79104000*
      d0^3*ep^3 - 92355200*d0^3*ep^2 - 59097344*d0^3*ep - 15984768*d0^3 + 
      103680*d0^2*ep^8 + 1700352*d0^2*ep^7 + 12045824*d0^2*ep^6 + 48136704*
      d0^2*ep^5 + 118656000*d0^2*ep^4 + 184710400*d0^2*ep^3 + 177292032*d0^2*
      ep^2 + 95908608*d0^2*ep + 22381056*d0^2 - 46080*d0*ep^9 - 850176*d0*ep^8
       - 6883328*d0*ep^7 - 32091136*d0*ep^6 - 94924800*d0*ep^5 - 184710400*d0*
      ep^4 - 236389376*d0*ep^3 - 191817216*d0*ep^2 - 89524224*d0*ep - 18302976
      *d0 + 9216*ep^10 + 188928*ep^9 + 1720832*ep^8 + 9168896*ep^7 + 31641600*
      ep^6 + 73884160*ep^5 + 118194688*ep^4 + 127878144*ep^3 + 89524224*ep^2
       + 36605952*ep + 6635520) + miBN1( - 2 + d0)*rat(18*d0 - 36*ep - 72,9
      *d0^4 - 72*d0^3*ep - 99*d0^3 + 216*d0^2*ep^2 + 594*d0^2*ep + 404*d0^2 - 
      288*d0*ep^3 - 1188*d0*ep^2 - 1616*d0*ep - 724*d0 + 144*ep^4 + 792*ep^3
       + 1616*ep^2 + 1448*ep + 480);

id miT111(d0?{,>4}) = Gam(4 + ep - 1/2*d0)^2*rat(-48,d0^6 - 12*d0^5*ep - 25*d0^5
       + 60*d0^4*ep^2 + 250*d0^4*ep + 254*d0^4 - 160*d0^3*ep^3 - 1000*d0^3*
      ep^2 - 2032*d0^3*ep - 1340*d0^3 + 240*d0^2*ep^4 + 2000*d0^2*ep^3 + 6096*
      d0^2*ep^2 + 8040*d0^2*ep + 3864*d0^2 - 192*d0*ep^5 - 2000*d0*ep^4 - 8128
      *d0*ep^3 - 16080*d0*ep^2 - 15456*d0*ep - 5760*d0 + 64*ep^6 + 800*ep^5 + 
      4064*ep^4 + 10720*ep^3 + 15456*ep^2 + 11520*ep + 3456) + miT111( - 2 + d0)
      *rat(3,d0^2 - 4*d0*ep - 5*d0 + 4*ep^2 + 10*ep + 6);
#endprocedure

#procedure dimstepup

id miD6(d0?{,<4}) = Gam(2 + ep - 1/2*d0)^3*rat( - 11*d0 + 22*ep + 16,3*d0^4 - 
      24*d0^3*ep - 24*d0^3 + 72*d0^2*ep^2 + 144*d0^2*ep + 72*d0^2 - 96*d0*ep^3
       - 288*d0*ep^2 - 288*d0*ep - 96*d0 + 48*ep^4 + 192*ep^3 + 288*ep^2 + 192
      *ep + 48) + miD6(2 + d0)*rat( - d0^3 + 6*d0^2*ep + 3*d0^2 - 12*d0*ep^2
       - 12*d0*ep - 2*d0 + 8*ep^3 + 12*ep^2 + 4*ep,4) + miD5(2 + d0)*rat(d0^3
       - 6*d0^2*ep - 3*d0^2 + 12*d0*ep^2 + 12*d0*ep + 2*d0 - 8*ep^3 - 12*ep^2
       - 4*ep,2) + miBN(2 + d0)*rat( - 9*d0^4 + 72*d0^3*ep + 27*d0^3 - 216*
      d0^2*ep^2 - 162*d0^2*ep - 26*d0^2 + 288*d0*ep^3 + 324*d0*ep^2 + 104*d0*
      ep + 8*d0 - 144*ep^4 - 216*ep^3 - 104*ep^2 - 16*ep,128*d0 - 256*ep - 256
      );

id miD5(d0?{,<4}) = Gam(2 + ep - 1/2*d0)^3*rat(863*d0^3 - 5178*d0^2*ep - 6061*
      d0^2 + 10356*d0*ep^2 + 24244*d0*ep + 13920*d0 - 6904*ep^3 - 24244*ep^2
       - 27840*ep - 10368,108*d0^6 - 1296*d0^5*ep - 1620*d0^5 + 6480*d0^4*ep^2
       + 16200*d0^4*ep + 9936*d0^4 - 17280*d0^3*ep^3 - 64800*d0^3*ep^2 - 79488
      *d0^3*ep - 31968*d0^3 + 25920*d0^2*ep^4 + 129600*d0^2*ep^3 + 238464*d0^2
      *ep^2 + 191808*d0^2*ep + 57024*d0^2 - 20736*d0*ep^5 - 129600*d0*ep^4 - 
      317952*d0*ep^3 - 383616*d0*ep^2 - 228096*d0*ep - 53568*d0 + 6912*ep^6 + 
      51840*ep^5 + 158976*ep^4 + 255744*ep^3 + 228096*ep^2 + 107136*ep + 20736
      ) + Gam(2 + ep - 1/2*d0)*miT111(2 + d0)*rat(4*d0^2 - 16*d0*ep - 4*d0 + 16*
      ep^2 + 8*ep,9*d0 - 18*ep - 36) + miD5(2 + d0)*rat( - 2*d0^4 + 16*d0^3*ep
       + 12*d0^3 - 48*d0^2*ep^2 - 72*d0^2*ep - 22*d0^2 + 64*d0*ep^3 + 144*d0*
      ep^2 + 88*d0*ep + 12*d0 - 32*ep^4 - 96*ep^3 - 88*ep^2 - 24*ep,9*d0 - 18*
      ep - 36) + miBN(2 + d0)*rat(117*d0^5 - 1170*d0^4*ep - 783*d0^4 + 4680*
      d0^3*ep^2 + 6264*d0^3*ep + 1634*d0^3 - 9360*d0^2*ep^3 - 18792*d0^2*ep^2
       - 9804*d0^2*ep - 1352*d0^2 + 9360*d0*ep^4 + 25056*d0*ep^3 + 19608*d0*
      ep^2 + 5408*d0*ep + 384*d0 - 3744*ep^5 - 12528*ep^4 - 13072*ep^3 - 5408*
      ep^2 - 768*ep,4608*d0^2 - 18432*d0*ep - 27648*d0 + 18432*ep^2 + 55296*ep
       + 36864);

id miD4(d0?{,<4}) = Gam(2 + ep - 1/2*d0)^3*rat( - 234*d0^3 + 1404*d0^2*ep + 
      1428*d0^2 - 2808*d0*ep^2 - 5712*d0*ep - 2812*d0 + 1872*ep^3 + 5712*ep^2
       + 5624*ep + 1848,27*d0^6 - 324*d0^5*ep - 378*d0^5 + 1620*d0^4*ep^2 + 
      3780*d0^4*ep + 2187*d0^4 - 4320*d0^3*ep^3 - 15120*d0^3*ep^2 - 17496*d0^3
      *ep - 6696*d0^3 + 6480*d0^2*ep^4 + 30240*d0^2*ep^3 + 52488*d0^2*ep^2 + 
      40176*d0^2*ep + 11448*d0^2 - 5184*d0*ep^5 - 30240*d0*ep^4 - 69984*d0*
      ep^3 - 80352*d0*ep^2 - 45792*d0*ep - 10368*d0 + 1728*ep^6 + 12096*ep^5
       + 34992*ep^4 + 53568*ep^3 + 45792*ep^2 + 20736*ep + 3888) + Gam(2 + ep
       - 1/2*d0)*miT111(2 + d0)*rat( - 2*d0^2 + 8*d0*ep + 2*d0 - 8*ep^2 - 4*ep,3
      *d0 - 6*ep - 9) + Gam( - ep + 1/2*d0)*Gam(2 + ep - 1/2*d0)^2*Gam(3 + 2*
      ep - d0)*rat(-8,3*d0^3 - 18*d0^2*ep - 18*d0^2 + 36*d0*ep^2 + 72*d0*ep + 
      36*d0 - 24*ep^3 - 72*ep^2 - 72*ep - 24) + Gam( - ep + 1/2*d0)*Gam(2 + ep
       - 1/2*d0)*Gam(3 + 2*ep - d0)^2*Gam(4 + 3*ep - 3/2*d0)*iGam(5 + 4*ep - 2
      *d0)*rat( - 169*d0^3 + 1014*d0^2*ep + 1066*d0^2 - 2028*d0*ep^2 - 4264*d0
      *ep - 2124*d0 + 1352*ep^3 + 4264*ep^2 + 4248*ep + 1296,54*d0^6 - 648*
      d0^5*ep - 693*d0^5 + 3240*d0^4*ep^2 + 6930*d0^4*ep + 3663*d0^4 - 8640*
      d0^3*ep^3 - 27720*d0^3*ep^2 - 29304*d0^3*ep - 10206*d0^3 + 12960*d0^2*
      ep^4 + 55440*d0^2*ep^3 + 87912*d0^2*ep^2 + 61236*d0^2*ep + 15804*d0^2 - 
      10368*d0*ep^5 - 55440*d0*ep^4 - 117216*d0*ep^3 - 122472*d0*ep^2 - 63216*
      d0*ep - 12888*d0 + 3456*ep^6 + 22176*ep^5 + 58608*ep^4 + 81648*ep^3 + 
      63216*ep^2 + 25776*ep + 4320) + miD4(2 + d0)*rat(d0^3 - 6*d0^2*ep - 3*
      d0^2 + 12*d0*ep^2 + 12*d0*ep + 2*d0 - 8*ep^3 - 12*ep^2 - 4*ep,2) + miE3(
      2 + d0)*rat(d0^3 - 6*d0^2*ep - 2*d0^2 + 12*d0*ep^2 + 8*d0*ep + d0 - 8*
      ep^3 - 8*ep^2 - 2*ep,3) + miBN1(2 + d0)*rat( - 9*d0^4 + 72*d0^3*ep + 
      27*d0^3 - 216*d0^2*ep^2 - 162*d0^2*ep - 26*d0^2 + 288*d0*ep^3 + 324*d0*
      ep^2 + 104*d0*ep + 8*d0 - 144*ep^4 - 216*ep^3 - 104*ep^2 - 16*ep,54*d0^2
       - 216*d0*ep - 270*d0 + 216*ep^2 + 540*ep + 324);

id miDN(d0?{,<4}) = Gam( - ep + 1/2*d0)*Gam(2 + ep - 1/2*d0)^2*Gam(3 + 2*ep - 
      d0)*rat(-16,d0^3 - 6*d0^2*ep - 6*d0^2 + 12*d0*ep^2 + 24*d0*ep + 12*d0 - 
      8*ep^3 - 24*ep^2 - 24*ep - 8) + Gam( - ep + 1/2*d0)*Gam(2 + ep - 1/2*d0)
      *Gam(3 + 2*ep - d0)^2*Gam(4 + 3*ep - 3/2*d0)*iGam(5 + 4*ep - 2*d0)*rat(
      40*d0 - 80*ep - 112,6*d0^4 - 48*d0^3*ep - 51*d0^3 + 144*d0^2*ep^2 + 306*
      d0^2*ep + 162*d0^2 - 192*d0*ep^3 - 612*d0*ep^2 - 648*d0*ep - 228*d0 + 96
      *ep^4 + 408*ep^3 + 648*ep^2 + 456*ep + 120) + Gam( - ep + 1/2*d0)^3*Gam(
       - 2 - 3*ep + 3/2*d0)*Gam(2 + ep - 1/2*d0)^2*Gam(4 + 3*ep - 3/2*d0)*
      iGam( - 1 - 2*ep + d0)^2*rat(-16,d0^3 - 6*d0^2*ep - 6*d0^2 + 12*d0*ep^2
       + 24*d0*ep + 12*d0 - 8*ep^3 - 24*ep^2 - 24*ep - 8) + miDN(2 + d0)*rat(
      d0^3 - 6*d0^2*ep - 3*d0^2 + 12*d0*ep^2 + 12*d0*ep + 2*d0 - 8*ep^3 - 12*
      ep^2 - 4*ep,4);

id miDM(d0?{,<4}) = Gam( - ep + 1/2*d0)*Gam(2 + ep - 1/2*d0)^2*Gam(3 + 2*ep - 
      d0)*rat(32,d0^3 - 6*d0^2*ep - 6*d0^2 + 12*d0*ep^2 + 24*d0*ep + 12*d0 - 8
      *ep^3 - 24*ep^2 - 24*ep - 8) + Gam( - ep + 1/2*d0)*Gam(2 + ep - 1/2*d0)*
      Gam(3 + 2*ep - d0)^2*Gam(4 + 3*ep - 3/2*d0)*iGam(5 + 4*ep - 2*d0)*rat(
       - 164*d0^2 + 656*d0*ep + 704*d0 - 656*ep^2 - 1408*ep - 672,18*d0^5 - 
      180*d0^4*ep - 177*d0^4 + 720*d0^3*ep^2 + 1416*d0^3*ep + 690*d0^3 - 1440*
      d0^2*ep^3 - 4248*d0^2*ep^2 - 4140*d0^2*ep - 1332*d0^2 + 1440*d0*ep^4 + 
      5664*d0*ep^3 + 8280*d0*ep^2 + 5328*d0*ep + 1272*d0 - 576*ep^5 - 2832*
      ep^4 - 5520*ep^3 - 5328*ep^2 - 2544*ep - 480) + miDM(2 + d0)*rat(d0^3 - 
      6*d0^2*ep - 3*d0^2 + 12*d0*ep^2 + 12*d0*ep + 2*d0 - 8*ep^3 - 12*ep^2 - 4
      *ep,2) + miE3(2 + d0)*rat( - d0^3 + 6*d0^2*ep + 2*d0^2 - 12*d0*ep^2 - 8*
      d0*ep - d0 + 8*ep^3 + 8*ep^2 + 2*ep,1);

id miE3(d0?{,<4}) = Gam( - ep + 1/2*d0)*Gam(2 + ep - 1/2*d0)^2*Gam(3 + 2*ep - 
      d0)*rat(-16,3*d0^3 - 18*d0^2*ep - 18*d0^2 + 36*d0*ep^2 + 72*d0*ep + 36*
      d0 - 24*ep^3 - 72*ep^2 - 72*ep - 24) + Gam( - ep + 1/2*d0)*Gam(2 + ep - 
      1/2*d0)*Gam(3 + 2*ep - d0)^2*Gam(4 + 3*ep - 3/2*d0)*iGam(5 + 4*ep - 2*d0
      )*rat( - 140*d0^2 + 560*d0*ep + 440*d0 - 560*ep^2 - 880*ep - 288,54*d0^5
       - 540*d0^4*ep - 531*d0^4 + 2160*d0^3*ep^2 + 4248*d0^3*ep + 2070*d0^3 - 
      4320*d0^2*ep^3 - 12744*d0^2*ep^2 - 12420*d0^2*ep - 3996*d0^2 + 4320*d0*
      ep^4 + 16992*d0*ep^3 + 24840*d0*ep^2 + 15984*d0*ep + 3816*d0 - 1728*ep^5
       - 8496*ep^4 - 16560*ep^3 - 15984*ep^2 - 7632*ep - 1440) + miE3(2 + d0)*
      rat(2*d0^3 - 12*d0^2*ep - 4*d0^2 + 24*d0*ep^2 + 16*d0*ep + 2*d0 - 16*
      ep^3 - 16*ep^2 - 4*ep,3);

id miBN(d0?{,<4}) = Gam(2 + ep - 1/2*d0)^3*rat( - 11*d0 + 22*ep + 16,d0^4 - 8*
      d0^3*ep - 8*d0^3 + 24*d0^2*ep^2 + 48*d0^2*ep + 24*d0^2 - 32*d0*ep^3 - 96
      *d0*ep^2 - 96*d0*ep - 32*d0 + 16*ep^4 + 64*ep^3 + 96*ep^2 + 64*ep + 16)
       + miBN(2 + d0)*rat( - 27*d0^4 + 216*d0^3*ep + 81*d0^3 - 648*d0^2*ep^2
       - 486*d0^2*ep - 78*d0^2 + 864*d0*ep^3 + 972*d0*ep^2 + 312*d0*ep + 24*d0
       - 432*ep^4 - 648*ep^3 - 312*ep^2 - 48*ep,128*d0 - 256*ep - 256);

id miD3(d0?{,<4}) = Gam(2 + ep - 1/2*d0)^3*rat(126*d0^3 - 756*d0^2*ep - 
      908*d0^2 + 1512*d0*ep^2 + 3632*d0*ep + 2100*d0 - 1008*ep^3 - 3632*ep^2
       - 4200*ep - 1512,9*d0^6 - 108*d0^5*ep - 126*d0^5 + 540*d0^4*ep^2 + 1260
      *d0^4*ep + 729*d0^4 - 1440*d0^3*ep^3 - 5040*d0^3*ep^2 - 5832*d0^3*ep - 
      2232*d0^3 + 2160*d0^2*ep^4 + 10080*d0^2*ep^3 + 17496*d0^2*ep^2 + 13392*
      d0^2*ep + 3816*d0^2 - 1728*d0*ep^5 - 10080*d0*ep^4 - 23328*d0*ep^3 - 
      26784*d0*ep^2 - 15264*d0*ep - 3456*d0 + 576*ep^6 + 4032*ep^5 + 11664*
      ep^4 + 17856*ep^3 + 15264*ep^2 + 6912*ep + 1296) + Gam( - ep + 1/2*d0)*
      Gam(2 + ep - 1/2*d0)^2*Gam(3 + 2*ep - d0)*rat(8,d0^4 - 8*d0^3*ep - 9*
      d0^3 + 24*d0^2*ep^2 + 54*d0^2*ep + 30*d0^2 - 32*d0*ep^3 - 108*d0*ep^2 - 
      120*d0*ep - 44*d0 + 16*ep^4 + 72*ep^3 + 120*ep^2 + 88*ep + 24) + Gam( - 
      ep + 1/2*d0)*Gam(2 + ep - 1/2*d0)*Gam(3 + 2*ep - d0)^2*Gam(4 + 3*ep - 3/
      2*d0)*iGam(5 + 4*ep - 2*d0)*rat(37*d0^2 - 148*d0*ep - 202*d0 + 148*ep^2
       + 404*ep + 276,6*d0^5 - 60*d0^4*ep - 69*d0^4 + 240*d0^3*ep^2 + 552*d0^3
      *ep + 315*d0^3 - 480*d0^2*ep^3 - 1656*d0^2*ep^2 - 1890*d0^2*ep - 714*
      d0^2 + 480*d0*ep^4 + 2208*d0*ep^3 + 3780*d0*ep^2 + 2856*d0*ep + 804*d0
       - 192*ep^5 - 1104*ep^4 - 2520*ep^3 - 2856*ep^2 - 1608*ep - 360) + Gam(
       - ep + 1/2*d0)^2*Gam(3 + 2*ep - d0)*Gam(4 + 3*ep - 3/2*d0)*rat(56*d0 - 
      112*ep - 144,3*d0^4 - 24*d0^3*ep - 27*d0^3 + 72*d0^2*ep^2 + 162*d0^2*ep
       + 90*d0^2 - 96*d0*ep^3 - 324*d0*ep^2 - 360*d0*ep - 132*d0 + 48*ep^4 + 
      216*ep^3 + 360*ep^2 + 264*ep + 72) + miD3(2 + d0)*rat(d0^3 - 6*d0^2*
      ep - 3*d0^2 + 12*d0*ep^2 + 12*d0*ep + 2*d0 - 8*ep^3 - 12*ep^2 - 4*ep,2)
       + miBN1(2 + d0)*rat(9*d0^4 - 72*d0^3*ep - 27*d0^3 + 216*d0^2*ep^2 + 
      162*d0^2*ep + 26*d0^2 - 288*d0*ep^3 - 324*d0*ep^2 - 104*d0*ep - 8*d0 + 
      144*ep^4 + 216*ep^3 + 104*ep^2 + 16*ep,18*d0 - 36*ep - 36);

id miBN1(d0?{,<4}) = Gam(2 + ep - 1/2*d0)^3*rat(120*d0 - 240*ep - 176,9*
      d0^4 - 72*d0^3*ep - 72*d0^3 + 216*d0^2*ep^2 + 432*d0^2*ep + 216*d0^2 - 
      288*d0*ep^3 - 864*d0*ep^2 - 864*d0*ep - 288*d0 + 144*ep^4 + 576*ep^3 + 
      864*ep^2 + 576*ep + 144) + miBN1(2 + d0)*rat(9*d0^4 - 72*d0^3*ep - 27
      *d0^3 + 216*d0^2*ep^2 + 162*d0^2*ep + 26*d0^2 - 288*d0*ep^3 - 324*d0*
      ep^2 - 104*d0*ep - 8*d0 + 144*ep^4 + 216*ep^3 + 104*ep^2 + 16*ep,18*d0
       - 36*ep - 36);

id miT111(d0?{,<4}) = Gam(2 + ep - 1/2*d0)^2*rat(4,d0^2 - 4*d0*ep - 4*d0 + 4*
      ep^2 + 8*ep + 4) + miT111(2 + d0)*rat(d0^2 - 4*d0*ep - d0 + 4*ep^2 + 2*ep,
      3);

#endprocedure



*--#[ shift4plus :
*
        #procedure shift4plus(DIM)

                #if({`DIM' % 2} == 0)
                        
* Expansion near d=4 + DIM -2*ep      
                        Multiply replace_(d,{4+`DIM'}-2*ep);
                        
* First argument dimension d=4-2ep for n+x*ep        
                        
                        id  Gam(n?,x?) =  Gam(n+x*(2-({4+`DIM'}-2*ep)/2));
                        id iGam(n?,x?) = iGam(n+x*(2-({4+`DIM'}-2*ep)/2));
                        
                        id mi?{miT111,miD6,miD5,miD4,miD3,miDN,miDM,miE3,miBN,miBN1} = mi({4+`DIM'});
                        
                        #if(`DIM' > 0)
                                #message Lowering DRR applied
                                
                                #do i=0,1                        
                                        #call dimstepdown  
                                        if(match(mi?{miT111,miD6,miD5,miD4,miD3,miDN,miDM,miE3,miBN,miBN1}(d0?{,>4}))) redefine i "0";
                                        .sort
                                #enddo
                                
                                #else
                                #message Raising DRR applied
                                
                                #do i=0,1                        
                                        #call dimstepup 
                                        if(match(mi?{miT111,miD6,miD5,miD4,miD3,miDN,miDM,miE3,miBN,miBN1}(d0?{,<4}))) redefine i "0";
                                        .sort
                                #enddo
                                
                        #endif        

* First remove numerical values
                        id Gam(x?pos_)  = fac_(x-1);
                        id iGam(x?pos_) = 1/fac_(x-1);
                        
                        SplitArg (ep),Gam,iGam;
                        id Gam(x?symbol_) = Gam(0,x);                        
                        id Gam(n?,x?)     = Gam(n,x/ep);

*                       When Gamma args are free from ep convert to d=4-2*ep
                        Multiply replace_(ep,2-d/2);                        
                        #call GammaArgToOne
*        Change notation to apply exp4d
                        id mi?{miT111,miD6,miD5,miD4,miD3,miDN,miDM,miE3,miBN,miBN1}(4) = mi;

                        #else
                        #message Expansion near Odd d not implemented                        
                #endif                        
        #endprocedure
*
*--#] shift4plus : 


*--#[ fillbntab : 
* 
* Example program to fill BNd and BNd0 table
* Can be used if tables for dala^n, n>11
* Inverse of eq.4 used with old routine for
* reduction.        
*         
* #-
* #include matad-ng.hh
* #undefine REDBNTAB
* S n3m2;

* * Fill table upto weight 11
* #do dp=0,11
* #do i6=1,2
* #do i5=1,`i6'
* #do i4=1,`i5'
* #do i3=1,`i4'
*         #message dala^`dp'-- `i3' `i4' `i5' `i6'
*         L dalaBN`dp' =dala^`dp'*x3^`i3'*x4^`i4'*x5^`i5'*x6^`i6';
        
*         #do i=1,1
*                 id,once dala^n1?pos_*x3^n3m2?*x4^n4?*x5^n5?*x6^n6? = 
*                 ((4*((n3m2+2)-1)-2*rat(d,1))/4/M^2/((n3m2+2)-1)*x3^((n3m2+2)-1)*x4^n4*x5^n5*x6^n6 -
*                 x3^(n3m2+2)*x4^n4*x5^n5*x6^n6)*4*M^2*(n3m2)*((n3m2+2)-1)*dala^(n1-1);
                
*                 if(count(dala,1)) redefine i "0";
*                 .sort
*         #enddo
        
*         Multiply replace_(x3,s1m,x4,s2m,x5,s5m,x6,s6m);
*         #call matad(3)
*         Multiply int0;
*         b int0;
*         .sort
        
*         #write<BNtbl> "id,only intbn*dala^`dp'*x3^`i3'*x4^`i4'*x5^`i5'*x6^`i6' = %e\n",dalaBN`dp'
        
*         drop;
*         .sort
* #enddo
* #enddo
* #enddo
* #enddo
* #enddo
* .end
*--#] fillbntab : 
