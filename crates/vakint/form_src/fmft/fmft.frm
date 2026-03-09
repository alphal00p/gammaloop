* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* -                                                                           - *
* -                                                                           - *
* -                    ███████╗███╗   ███╗███████╗████████╗                   - *
* -                    ██╔════╝████╗ ████║██╔════╝╚══██╔══╝                   - *
* -                    █████╗  ██╔████╔██║█████╗     ██║                      - *
* -                    ██╔══╝  ██║╚██╔╝██║██╔══╝     ██║                      - *
* -                    ██║     ██║ ╚═╝ ██║██║        ██║                      - *
* -                    ╚═╝     ╚═╝     ╚═╝╚═╝        ╚═╝                      - *
* -                                                                           - *
* -                                                                           - *
* -                                                                           - *
* -                     Fully Massive Four-loop Tadpoles                      - *
* -                                                                           - *
* -                                                                           - *
* -                                                                           - *
* -                                                                           - *
* -                   Contact:                                                - *
* -                           Andrey Pikelner                                 - *
* -                           pikelner[at]theor.jinr.ru                       - *
* -                                                                           - *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *




* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* -                                                                           - *
* ================================== Options ================================== *
* -                                                                           - *
* -           We do not load any other program files and not                  - * 
* -           need to set any other options in this version                   - * 
* -                                                                           - *
* -                                                                           - *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* -                                                                           - *
* ======================= Topologies and momenta mapping ====================== *
* -                                                                           - *
* -                                                                           - *
* -                   Topo H                         Topo X                   - *
* -                                                   p4                      - *
* -              _______\_______               _______/________               - *
* -             /   |   p1   |   \            /      \  /      \              - *
* -            /    ^p8    p9|    \          /        /\        \             - *
* -           /     |        V     \        /    p7|/   |\p6     \            - *
* -          /|p2   |___/____|     |        |     /____\___\   p3|\           - *
* -           |     |   p5   |   p3|/      \|p2   |   p5   |     |            - *
* -           \     |      p7^     /        \    /|p9      |     /            - *
* -            \    V p6     |    /          \    |        V p8 /             - *
* -             \ __|____/___|___/            \ __|____/___|___/              - *
* -                     p4                            p10                     - *
* -                                                                           - *
* -                                                                           - *
* -                                                                           - *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* -                                                                           - *
* =========================== Input and examples ============================== *
* -                                                                           - *
* -  d1,...,d10 - massive denominators d_i = p_i^2-1                          - *
* -  p1,...,p10 - momentum to form scalar products in numerator               - *
* -                                                                           - *
* -                                                                           - *
* -  Sample input:                                                            - *
* -  L ex = p2.p3/d1^-2/d2^2/d3^0/d4^1/d5^1/d6^1/d7^1/d8^1/d9^3/d10^0;        - *
* -                                                                           - *
* -  Start reduction:                                                         - *
* -  #call fmft                                                               - *
* -                                                                           - *
* -  Expand in ep near d=4 up to ep^1:                                        - *
* -  #call exp4d(1)                                                           - *
* -                                                                           - *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


* Switching off statistics from each thread
Off threadstats;

* Massive denominators
S d1,...,d10;
* Scalar products in numerator
V p1,...,p10;
S n1,...,n10;

CF Fab;

CTensor ftensor,dd;
* Loop momenta in numerator
V k1,k2,k3,k4;

V tk4,mtk4,mtp;

* TOPO symbols
S intFMFT,intX,intH,intBMW,intFG,int1d;

* TOPO symbols for red2l reduction
S intF,intV,intJ,intT2,intGxG,intGxT,intTxG,intTxT;

* averts
S [sqrt(x)],[x],n,s;
V [k-q],[tk1-tk2];

CF num,poinv;

CF Oep;
* Master integrals
S PR0,PR1,PR2,PR3,PR4,PR4d,PR5,PR6,PR7,PR8,PR9,PR9d,PR9x,PR10,PR11,PR11d,PR12,PR13,PR14,PR15;

* Expansion in 4d
S PR0ep0,PR12ep0,PR15ep0,z3,ep,D6,z4,PR11ep0,z5,PR11dep0,PR14ep0,S2,z2,T1ep,PR13ep0,PR10ep0,PR9ep0,PR9ep1,PR9dep0,PR9dep1,PR8ep0,PR4dep0,PR4ep0,PR7ep1,T1ep2,PR6ep0,PR5ep0,PR4ep1,PR4dep1,PR3ep1,PR2ep1,PR1ep1;

S T2x111,Gx11,T1x1;

S x1,x2,x3;
S j3;
*
****
* Tarcer declarations
**** 
*
* momentum rules:
*   k1->tk1, k2->tk2, p->tp;
*   tarC1,...,tarC5 with common mass m^2=tmm
*   C_i=p_i^2 - tmm
*
*   Eucledian: tmm = -m^2
*   Minkowski: tmm =  m^2
*

S v,w,x,y,z,n1,...,n5,n;
S tarC1,...,tarC5,pp;
S i,j,k,r;
S d,dp;
V tk1,tk2,tp;
CF TFI;
* Functions with dim shift first argument like F(dp,...)
CF T2,G,T1;

CF Pochhammer,PochhammerINV,num,den,rat;
CF sdim;
* New denominators
S [pp-tmm],[pp-3*tmm],[pp-4*tmm],[pp-9*tmm],mm,mden,tmm;

PolyRatFun rat;


* One dim rec rel 

CF m0denJ0001X01110, m0denJ0001X11110, m0denJ0001X00110,
m0denJ0001X01110, m0denJ0001X11110, m0denJ0011X01100,
m0denJ0011X01110, m0denJ0011X10110, m0denJ0011X11100,
m0denJ0011X11110, m0denJ0011X00110, m0denJ0011X01100,
m0denJ0011X01110, m0denJ0011X10110, m0denJ0011X11100,
m0denJ0011X11110, m0denJ0021X00110, m0denJ0021X10110,
m0denJ0031X10110, m0denJ1000X01100, m0denJ1000X01100,
m0denJ1001X01110, m0denJ1001X11110, m0denJ1001X01110,
m0denJ1001X11110, m0denJ1010X01100, m0denJ1010X11100,
m0denJ1010X01100, m0denJ1010X11100, m0denJ1011X01110,
m0denJ1011X11110, m0denJ1011X01110, m0denJ1011X11110,
m0denJ2010X01100, m0denJ2010X11100, m0denJ2011X11110;

CF m1denJ0001X01110, m1denJ0001X11110, m1denJ0001X00110,
m1denJ0001X01110, m1denJ0001X11110, m1denJ0011X01100,
m1denJ0011X01110, m1denJ0011X10110, m1denJ0011X11100,
m1denJ0011X11110, m1denJ0011X00110, m1denJ0011X01100,
m1denJ0011X01110, m1denJ0011X10110, m1denJ0011X11100,
m1denJ0011X11110, m1denJ0021X00110, m1denJ0021X10110,
m1denJ0031X10110, m1denJ1000X01100, m1denJ1000X01100,
m1denJ1001X01110, m1denJ1001X11110, m1denJ1001X01110,
m1denJ1001X11110, m1denJ1010X01100, m1denJ1010X11100,
m1denJ1010X01100, m1denJ1010X11100, m1denJ1011X01110,
m1denJ1011X11110, m1denJ1011X01110, m1denJ1011X11110,
m1denJ2010X01100, m1denJ2010X11100, m1denJ2011X11110;

CF m3denJ0001X01110, m3denJ0001X11110, m3denJ0001X00110,
m3denJ0001X01110, m3denJ0001X11110, m3denJ0011X01100,
m3denJ0011X01110, m3denJ0011X10110, m3denJ0011X11100,
m3denJ0011X11110, m3denJ0011X00110, m3denJ0011X01100,
m3denJ0011X01110, m3denJ0011X10110, m3denJ0011X11100,
m3denJ0011X11110, m3denJ0021X00110, m3denJ0021X10110,
m3denJ0031X10110, m3denJ1000X01100, m3denJ1000X01100,
m3denJ1001X01110, m3denJ1001X11110, m3denJ1001X01110,
m3denJ1001X11110, m3denJ1010X01100, m3denJ1010X11100,
m3denJ1010X01100, m3denJ1010X11100, m3denJ1011X01110,
m3denJ1011X11110, m3denJ1011X01110, m3denJ1011X11110,
m3denJ2010X01100, m3denJ2010X11100, m3denJ2011X11110;

CF m4denJ0001X01110, m4denJ0001X11110, m4denJ0001X00110,
m4denJ0001X01110, m4denJ0001X11110, m4denJ0011X01100,
m4denJ0011X01110, m4denJ0011X10110, m4denJ0011X11100,
m4denJ0011X11110, m4denJ0011X00110, m4denJ0011X01100,
m4denJ0011X01110, m4denJ0011X10110, m4denJ0011X11100,
m4denJ0011X11110, m4denJ0021X00110, m4denJ0021X10110,
m4denJ0031X10110, m4denJ1000X01100, m4denJ1000X01100,
m4denJ1001X01110, m4denJ1001X11110, m4denJ1001X01110,
m4denJ1001X11110, m4denJ1010X01100, m4denJ1010X11100,
m4denJ1010X01100, m4denJ1010X11100, m4denJ1011X01110,
m4denJ1011X11110, m4denJ1011X01110, m4denJ1011X11110,
m4denJ2010X01100, m4denJ2010X11100, m4denJ2011X11110;

CF m9denJ0001X01110, m9denJ0001X11110, m9denJ0001X00110,
m9denJ0001X01110, m9denJ0001X11110, m9denJ0011X01100,
m9denJ0011X01110, m9denJ0011X10110, m9denJ0011X11100,
m9denJ0011X11110, m9denJ0011X00110, m9denJ0011X01100,
m9denJ0011X01110, m9denJ0011X10110, m9denJ0011X11100,
m9denJ0011X11110, m9denJ0021X00110, m9denJ0021X10110,
m9denJ0031X10110, m9denJ1000X01100, m9denJ1000X01100,
m9denJ1001X01110, m9denJ1001X11110, m9denJ1001X01110,
m9denJ1001X11110, m9denJ1010X01100, m9denJ1010X11100,
m9denJ1010X01100, m9denJ1010X11100, m9denJ1011X01110,
m9denJ1011X11110, m9denJ1011X01110, m9denJ1011X11110,
m9denJ2010X01100, m9denJ2010X11100, m9denJ2011X11110;

CF numJ0001X01110, numJ0001X11110, numJ0001X00110, numJ0001X01110, 
numJ0001X11110, numJ0011X01100, numJ0011X01110, numJ0011X10110, 
numJ0011X11100, numJ0011X11110, numJ0011X00110, numJ0011X01100, 
numJ0011X01110, numJ0011X10110, numJ0011X11100, numJ0011X11110, 
numJ0021X00110, numJ0021X10110, numJ0031X10110, numJ1000X01100, 
numJ1000X01100, numJ1001X01110, numJ1001X11110, numJ1001X01110, 
numJ1001X11110, numJ1010X01100, numJ1010X11100, numJ1010X01100, 
numJ1010X11100, numJ1011X01110, numJ1011X11110, numJ1011X01110, 
numJ1011X11110, numJ2010X01100, numJ2010X11100, numJ2011X11110;

.global




*--#[ dzero :
* If facet shrinked to point integral is zero
* 3 edges
#procedure dzero3(N1,N2,N3)
        if((count(tarC`N1',-1) <= 0) && (count(tarC`N2',-1) <= 0) && (count(tarC`N3',-1) <= 0)) Discard;
#endprocedure        
* 4 edges
#procedure dzero4(N1,N2,N3,N4)
        if((count(tarC`N1',-1) <= 0) && (count(tarC`N2',-1) <= 0) && (count(tarC`N3',-1) <= 0) && (count(tarC`N4',-1) <= 0)) Discard;
#endprocedure        
*--#] dzero :

*--#[ zeroTFI :
#procedure zeroTFI(topo)
        if(count(int`topo',1));
        #call dzero3(2,4,5)        
        #call dzero3(1,3,5)                
        #call dzero4(1,2,3,4)                
        #call dzero4(2,3,4,5)                        
        #call dzero4(1,2,3,5)                                        
        #call dzero4(1,2,4,5)                                                
        #call dzero4(1,3,4,5)
        endif;        
#endprocedure        
*--#] zeroTFI :



#procedure GGdecouple(topo)
*       <-9->        
*       B()*B()
*       Gamma[k + a_] -> Po[a, k]*Gamma[a]
        if(count(int`topo',1));
        if(count(tarC5,-1) == 0);
        id tk1.tk1^v?*tk2.tk2^w?*tp.tk1^x?*tp.tk2^y?*tk1.tk2^z?pos_/tarC1^n1?/tarC2^n2?/tarC3^n3?/tarC4^n4? =
        sum_(k, 0, integer_(z/2), sum_(i, 0, k, sum_(j, 0, k, sign_(i + j)*tk1.tk1^(i + v)*tk2.tk2^(j + w)*pp^(i + j - z)*tp.tk1^(-2*i + x + z)*tp.tk2^(-2*j + y + z)*binom_(k, i)*binom_(k, j)*binom_(z, 2*k)*Pochhammer(k,1/2)*PochhammerINV(k,-1/2 + d/2)/(tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4))));
        endif;
        endif;        
#endprocedure


* Here we cancel each scalar product with its corresponding denominator
*--#[ spcontract1 :
#procedure spcontract1(topo,nc)
* 
* 
*       Red SP 1
* 
*         

        #$repcount = 1;                
        #do irep=1,1       
                #$irep = 1;
                if(count(int`topo',1));
                
*               (tk1^1)--
                if(count(tk1.tk1,1) > count(tarC1,-1));
                id,ifmatch->sortsp tk1.tk1^v?pos_*tk2.tk2^w?*tp.tk1^x?*tp.tk2^y?*tk1.tk2^z?/tarC1^n1?pos_/tarC2^n2?/tarC3^n3?/tarC4^n4?/tarC5^n5? =
                sum_(i, 0, n1, (tarC1^(i - n1)*tk1.tk1^(-n1 + v)*tk1.tk2^z*tk2.tk2^w*(tmm)^(-i + n1)*tp.tk1^x*tp.tk2^y*binom_(n1, i))/(tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5));
                endif;
                
                if(count(tk1.tk1,1) <= count(tarC1,-1));
                id,ifmatch->sortsp tk1.tk1^v?pos_*tk2.tk2^w?*tp.tk1^x?*tp.tk2^y?*tk1.tk2^z?/tarC1^n1?pos_/tarC2^n2?/tarC3^n3?/tarC4^n4?/tarC5^n5? =
                sum_(i, 0, v, (tarC1^(i - n1)*tk1.tk2^z*tk2.tk2^w*(tmm)^(-i + v)*tp.tk1^x*tp.tk2^y*binom_(v, i))/(tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5));
                endif;
                
                
                
*               (tk2^2)--
                if(count(tk2.tk2,1) > count(tarC2,-1));
                id,ifmatch->sortsp tk1.tk1^v?*tk2.tk2^w?pos_*tp.tk1^x?*tp.tk2^y?*tk1.tk2^z?/tarC1^n1?/tarC2^n2?pos_/tarC3^n3?/tarC4^n4?/tarC5^n5? =
                sum_(i, 0, n2, (tarC2^(i - n2)*tk1.tk1^v*tk1.tk2^z*tk2.tk2^(-n2 + w)*(tmm)^(-i + n2)*tp.tk1^x*tp.tk2^y*binom_(n2, i))/(tarC1^n1*tarC3^n3*tarC4^n4*tarC5^n5));
                endif;
                
                
                if(count(tk2.tk2,1) <= count(tarC2,-1));
                id,ifmatch->sortsp tk1.tk1^v?*tk2.tk2^w?pos_*tp.tk1^x?*tp.tk2^y?*tk1.tk2^z?/tarC1^n1?/tarC2^n2?pos_/tarC3^n3?/tarC4^n4?/tarC5^n5? =
                sum_(i, 0, w, (tarC2^(i - n2)*tk1.tk1^v*tk1.tk2^z*(tmm)^(-i + w)*tp.tk1^x*tp.tk2^y*binom_(w, i))/(tarC1^n1*tarC3^n3*tarC4^n4*tarC5^n5));



*               (tk1.tk2)--
                id,ifmatch->sortsp tk1.tk1^v?*tk2.tk2^w?*tp.tk1^x?*tp.tk2^y?*tk1.tk2^z?pos_/tarC1^n1?/tarC2^n2?/tarC3^n3?/tarC4^n4?/tarC5^n5?pos_ =
                - (tarC5^(1 - n5)*tk1.tk1^v*tk1.tk2^(-1 + z)*tk2.tk2^w*tp.tk1^x*tp.tk2^y)/(2*tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4)
                + (tk1.tk1^(1 + v)*tk1.tk2^(-1 + z)*tk2.tk2^w*tp.tk1^x*tp.tk2^y)/(2*tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5)
                + (tk1.tk1^v*tk1.tk2^(-1 + z)*tk2.tk2^(1 + w)*tp.tk1^x*tp.tk2^y)/(2*tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5)
                - (tk1.tk1^v*tk1.tk2^(-1 + z)*tk2.tk2^w*tmm*tp.tk1^x*tp.tk2^y)/(2*tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5);
                endif;
                
                endif;
                goto endrec;                
                la sortsp;
                $irep = 0;
                la endrec;
                
                ModuleOption,minimum,$irep;
                .sort:red-SP1-`nc'-`$repcount++';
                #redefine irep "`$irep'"
        #enddo        
#endprocedure
*--#] spcontract1 :



* We reduce all k1.p and k2.p to C1C2C5
#procedure redTad125(topo)
        #$repcount = 1;                
        #do irep=1,1       
                #$irep = 1;
                if(count(int`topo',1));
*       zero if odd(x+y)=true
                if((count(tarC3,-1)) == 0 && (count(tarC4,-1) == 0));
                id tk1.tk1^v?*tk2.tk2^w?*tp.tk1^x?pos0_*tp.tk2^y?pos0_*tk1.tk2^z?/tarC1^n1?pos_/tarC2^n2?pos_/tarC5^n5? =
                mod_(x+y+1,2)*tk1.tk1^v*tk2.tk2^w*tp.tk1^x*tp.tk2^y*tk1.tk2^z/tarC1^n1/tarC2^n2/tarC5^n5;
                endif;
                
*       even(x+y)
                if((count(tarC3,-1)) == 0 && (count(tarC4,-1) == 0));
                id,ifmatch->sortsp tk1.tk1^v?*tk2.tk2^w?*tp.tk1^x?{>1}*tp.tk2^y?pos_*tk1.tk2^z?/tarC1^n1?/tarC2^n2?/tarC5^n5?pos_ =
                mod_(x+y+1,2)*(pp*((tk1.tk1^(1 + v)*tk1.tk2^z*tk2.tk2^w*tp.tk1^(-2 + x)*tp.tk2^y*(-1 + x))/(tarC1^n1*tarC2^n2*tarC5^n5) + (tk1.tk1^v*tk1.tk2^(1 + z)*tk2.tk2^w*tp.tk1^(-1 + x)*tp.tk2^(-1 + y)*y)/(tarC1^n1*tarC2^n2*tarC5^n5)))*rat(1,-2 + d + x + y);
                endif;
                
*       even(x+y)
                if((count(tarC3,-1)) == 0 && (count(tarC4,-1) == 0));
                id,ifmatch->sortsp tk1.tk1^v?*tk2.tk2^w?*tp.tk1^x?pos_*tp.tk2^y?{>1}*tk1.tk2^z?/tarC1^n1?/tarC2^n2?/tarC5^n5?pos_ =
                mod_(x+y+1,2)*(pp*((tk1.tk1^v*tk1.tk2^(1 + z)*tk2.tk2^w*tp.tk1^(-1 + x)*tp.tk2^(-1 + y)*x)/(tarC1^n1*tarC2^n2*tarC5^n5) + (tk1.tk1^v*tk1.tk2^z*tk2.tk2^(1 + w)*tp.tk1^x*tp.tk2^(-2 + y)*(-1 + y))/(tarC1^n1*tarC2^n2*tarC5^n5)))*rat(1,-2 + d + x + y);
                endif;

                endif;
                goto endrec;                
                la sortsp;
                $irep = 0;
                la endrec;
                
                ModuleOption,minimum,$irep;
                .sort:red-Tad125-`$repcount++';
                #redefine irep "`$irep'"
        #enddo        
#endprocedure

#procedure convto124(topo)
        if(count(int`topo',1));
*       First convert from 135 to 245
*       <-21->
        if((count(tk2.tk2,1) > 0) && (count(tk1.tk2,1) == 0) && (count(tarC2,-1) == 0) && (count(tarC4,-1) == 0))
        Multiply replace_(tk1,tk2,tk2,tk1, tarC1,tarC2, tarC3,tarC4);

*       Second convert from 245 to 124
*       <-17->
        if((count(tp.tk1,1) == 0) && (count(tk1.tk2,1) == 0) && (count(tarC1,-1) == 0)   && (count(tarC3,-1) == 0));
        id tk1.tk1^v?pos_*tk2.tk2^w?*tp.tk2^y?/tarC2^n2?/tarC4^n4?/tarC5^n5? =
        sum_(i, 0, v, sum_(j, 0, i, (2^(-1 + j)*(1 + (-1)^j)*tk1.tk1^(-i + v)*tk1.tk2^j*tk2.tk2^(i - j + w)*tp.tk2^y*binom_(i, j)*binom_(v, i))/(tarC1^n5*tarC2^n2*tarC4^n4)));
        endif;

*       <-20->
        if((count(tk1.tk2,1) == 0) && (count(tarC1,-1) == 0)   && (count(tarC3,-1) == 0));
        id tk1.tk1^v?pos_*tk2.tk2^w?*tp.tk1^x?*tp.tk2^y?/tarC2^n2?/tarC4^n4?/tarC5^n5?pos_ =
        sum_(j, 0, x, sum_(i, 0, v, sum_(r, 0, i, (2^r*tk1.tk1^(-i + v)*tk1.tk2^r*tk2.tk2^(i - r + w)*tp.tk1^j*tp.tk2^(-j + x + y)*binom_(i, r)*binom_(v, i)*binom_(x, j))/(tarC1^n5*tarC2^n2*tarC4^n4))));
        endif;        
        endif;
        .sort conv124;        
#endprocedure





*--#[ redF :
#procedure redF

* eq.35 (1-)
        id,ifmatch->dopartfrac sdim(dp?)/tarC1^n1?{>1}/tarC2^n2?pos_/tarC3^n3?pos_/tarC4^n4?pos_/tarC5^n5?pos_=
        ((tarC5^(1 - n5)*rat(-1, 2))/(tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*[pp-3*tmm]) +
        (tarC3^(1 - n3)*((tmm*rat(-3, 2))/(pp*[pp-3*tmm]) + rat(1, 1)/[pp-3*tmm]))/
        (tarC1^n1*tarC2^n2*tarC4^n4*tarC5^n5) + (tarC2^(1 - n2)*rat(1, 2))/
        (tarC1^n1*tarC3^n3*tarC4^n4*tarC5^n5*[pp-3*tmm]) + (tarC1^(1 - n1)*tarC3^(-1 - n3)*tarC5^(1 - n5)*
        rat(-n3, 2*(-1 + n1)))/(tarC2^n2*tarC4^n4*[pp-3*tmm]) +
        (tarC1^(1 - n1)*tarC3^(-1 - n3)*tarC4^(1 - n4)*rat(n3, 2*(-1 + n1)))/(tarC2^n2*tarC5^n5*[pp-3*tmm]) +
        (tarC1^(2 - n1)*tarC3^(-1 - n3)*tmm*rat(3*n3, 2*(-1 + n1)))/
        (tarC2^n2*tarC4^n4*tarC5^n5*pp*[pp-3*tmm]) +
        (tarC1^(1 - n1)*((tmm*rat(-3 + 3*n1 - 3*n3, 2*(-1 + n1)))/(pp*[pp-3*tmm]) +
        rat(-1 - d - dp + n1 + 3*n3, 2*(-1 + n1))/[pp-3*tmm]))/(tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5) +
        (tarC1^(1 - n1)*tarC2^(1 - n2)*tarC5^(-1 - n5)*tmm*rat(-3*n5, 2*(-1 + n1)))/
        (tarC3^n3*tarC4^n4*pp*[pp-3*tmm]) + (tarC1^(1 - n1)*tarC3^(1 - n3)*tarC5^(-1 - n5)*
        ((tmm*rat(-3*n5, 2*(-1 + n1)))/(pp*[pp-3*tmm]) + rat(n5, -1 + n1)/[pp-3*tmm]))/
        (tarC2^n2*tarC4^n4) + (tarC1^(2 - n1)*tarC5^(-1 - n5)*tmm*rat(3*n5, 2*(-1 + n1)))/
        (tarC2^n2*tarC3^n3*tarC4^n4*pp*[pp-3*tmm]) + (tarC1^(1 - n1)*tarC4^(1 - n4)*tarC5^(-1 - n5)*
        (rat(-n5, -1 + n1)/[pp-3*tmm] + (tmm*rat(3*n5, 2*(-1 + n1)))/(pp*[pp-3*tmm])))/
        (tarC2^n2*tarC3^n3))*sdim(dp);
        
* eq.35 (2-)
        id,ifmatch->dopartfrac sdim(dp?)/tarC1^n1?pos_/tarC2^n2?{>1}/tarC3^n3?pos_/tarC4^n4?pos_/tarC5^n5?pos_=
        ((tarC5^(1 - n5)*rat(-1, 2))/(tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*[pp-3*tmm]) +
        (tarC4^(1 - n4)*((tmm*rat(-3, 2))/(pp*[pp-3*tmm]) + rat(1, 1)/[pp-3*tmm]))/
        (tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5) + (tarC1^(1 - n1)*rat(1, 2))/
        (tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5*[pp-3*tmm]) + (tarC2^(1 - n2)*tarC4^(-1 - n4)*tarC5^(1 - n5)*
        rat(-n4, 2*(-1 + n2)))/(tarC1^n1*tarC3^n3*[pp-3*tmm]) +
        (tarC2^(1 - n2)*tarC3^(1 - n3)*tarC4^(-1 - n4)*rat(n4, 2*(-1 + n2)))/(tarC1^n1*tarC5^n5*[pp-3*tmm]) +
        (tarC2^(2 - n2)*tarC4^(-1 - n4)*tmm*rat(3*n4, 2*(-1 + n2)))/
        (tarC1^n1*tarC3^n3*tarC5^n5*pp*[pp-3*tmm]) +
        (tarC2^(1 - n2)*((tmm*rat(-3 + 3*n2 - 3*n4, 2*(-1 + n2)))/(pp*[pp-3*tmm]) +
        rat(-1 - d - dp + n2 + 3*n4, 2*(-1 + n2))/[pp-3*tmm]))/(tarC1^n1*tarC3^n3*tarC4^n4*tarC5^n5) +
        (tarC1^(1 - n1)*tarC2^(1 - n2)*tarC5^(-1 - n5)*tmm*rat(-3*n5, 2*(-1 + n2)))/
        (tarC3^n3*tarC4^n4*pp*[pp-3*tmm]) + (tarC2^(1 - n2)*tarC4^(1 - n4)*tarC5^(-1 - n5)*
        ((tmm*rat(-3*n5, 2*(-1 + n2)))/(pp*[pp-3*tmm]) + rat(n5, -1 + n2)/[pp-3*tmm]))/
        (tarC1^n1*tarC3^n3) + (tarC2^(2 - n2)*tarC5^(-1 - n5)*tmm*rat(3*n5, 2*(-1 + n2)))/
        (tarC1^n1*tarC3^n3*tarC4^n4*pp*[pp-3*tmm]) + (tarC2^(1 - n2)*tarC3^(1 - n3)*tarC5^(-1 - n5)*
        (rat(-n5, -1 + n2)/[pp-3*tmm] + (tmm*rat(3*n5, 2*(-1 + n2)))/(pp*[pp-3*tmm])))/
        (tarC1^n1*tarC4^n4))*sdim(dp);
        
* eq.35 (3-)
        id,ifmatch->dopartfrac sdim(dp?)/tarC1^n1?pos_/tarC2^n2?pos_/tarC3^n3?{>1}/tarC4^n4?pos_/tarC5^n5?pos_=
        ((tarC5^(1 - n5)*rat(-1, 2))/(tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*[pp-3*tmm]) +
        (tarC1^(1 - n1)*((tmm*rat(-3, 2))/(pp*[pp-3*tmm]) + rat(1, 1)/[pp-3*tmm]))/
        (tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5) + (tarC4^(1 - n4)*rat(1, 2))/
        (tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5*[pp-3*tmm]) + (tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC5^(1 - n5)*
        rat(-n1, 2*(-1 + n3)))/(tarC2^n2*tarC4^n4*[pp-3*tmm]) +
        (tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3)*rat(n1, 2*(-1 + n3)))/(tarC4^n4*tarC5^n5*[pp-3*tmm]) +
        (tarC1^(-1 - n1)*tarC3^(2 - n3)*tmm*rat(3*n1, 2*(-1 + n3)))/
        (tarC2^n2*tarC4^n4*tarC5^n5*pp*[pp-3*tmm]) +
        (tarC3^(1 - n3)*(rat(-1 - d - dp + 3*n1 + n3, 2*(-1 + n3))/[pp-3*tmm] +
        (tmm*rat(-3 - 3*n1 + 3*n3, 2*(-1 + n3)))/(pp*[pp-3*tmm])))/
        (tarC1^n1*tarC2^n2*tarC4^n4*tarC5^n5) + (tarC3^(1 - n3)*tarC4^(1 - n4)*tarC5^(-1 - n5)*tmm*
        rat(-3*n5, 2*(-1 + n3)))/(tarC1^n1*tarC2^n2*pp*[pp-3*tmm]) +
        (tarC1^(1 - n1)*tarC3^(1 - n3)*tarC5^(-1 - n5)*((tmm*rat(-3*n5, 2*(-1 + n3)))/(pp*[pp-3*tmm]) +
        rat(n5, -1 + n3)/[pp-3*tmm]))/(tarC2^n2*tarC4^n4) +
        (tarC3^(2 - n3)*tarC5^(-1 - n5)*tmm*rat(3*n5, 2*(-1 + n3)))/
        (tarC1^n1*tarC2^n2*tarC4^n4*pp*[pp-3*tmm]) + (tarC2^(1 - n2)*tarC3^(1 - n3)*tarC5^(-1 - n5)*
        (rat(-n5, -1 + n3)/[pp-3*tmm] + (tmm*rat(3*n5, 2*(-1 + n3)))/(pp*[pp-3*tmm])))/
        (tarC1^n1*tarC4^n4))*sdim(dp);
        
* eq.35 (4-)
        id,ifmatch->dopartfrac sdim(dp?)/tarC1^n1?pos_/tarC2^n2?pos_/tarC3^n3?pos_/tarC4^n4?{>1}/tarC5^n5?pos_=
        ((tarC5^(1 - n5)*rat(-1, 2))/(tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*[pp-3*tmm]) +
        (tarC2^(1 - n2)*((tmm*rat(-3, 2))/(pp*[pp-3*tmm]) + rat(1, 1)/[pp-3*tmm]))/
        (tarC1^n1*tarC3^n3*tarC4^n4*tarC5^n5) + (tarC3^(1 - n3)*rat(1, 2))/
        (tarC1^n1*tarC2^n2*tarC4^n4*tarC5^n5*[pp-3*tmm]) + (tarC2^(-1 - n2)*tarC4^(1 - n4)*tarC5^(1 - n5)*
        rat(-n2, 2*(-1 + n4)))/(tarC1^n1*tarC3^n3*[pp-3*tmm]) +
        (tarC1^(1 - n1)*tarC2^(-1 - n2)*tarC4^(1 - n4)*rat(n2, 2*(-1 + n4)))/(tarC3^n3*tarC5^n5*[pp-3*tmm]) +
        (tarC2^(-1 - n2)*tarC4^(2 - n4)*tmm*rat(3*n2, 2*(-1 + n4)))/
        (tarC1^n1*tarC3^n3*tarC5^n5*pp*[pp-3*tmm]) +
        (tarC4^(1 - n4)*(rat(-1 - d - dp + 3*n2 + n4, 2*(-1 + n4))/[pp-3*tmm] +
        (tmm*rat(-3 - 3*n2 + 3*n4, 2*(-1 + n4)))/(pp*[pp-3*tmm])))/
        (tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5) + (tarC3^(1 - n3)*tarC4^(1 - n4)*tarC5^(-1 - n5)*tmm*
        rat(-3*n5, 2*(-1 + n4)))/(tarC1^n1*tarC2^n2*pp*[pp-3*tmm]) +
        (tarC2^(1 - n2)*tarC4^(1 - n4)*tarC5^(-1 - n5)*((tmm*rat(-3*n5, 2*(-1 + n4)))/(pp*[pp-3*tmm]) +
        rat(n5, -1 + n4)/[pp-3*tmm]))/(tarC1^n1*tarC3^n3) +
        (tarC4^(2 - n4)*tarC5^(-1 - n5)*tmm*rat(3*n5, 2*(-1 + n4)))/
        (tarC1^n1*tarC2^n2*tarC3^n3*pp*[pp-3*tmm]) + (tarC1^(1 - n1)*tarC4^(1 - n4)*tarC5^(-1 - n5)*
        (rat(-n5, -1 + n4)/[pp-3*tmm] + (tmm*rat(3*n5, 2*(-1 + n4)))/(pp*[pp-3*tmm])))/
        (tarC2^n2*tarC3^n3))*sdim(dp);
        

* eq.41 (5-)
        id,ifmatch->dopartfrac sdim(dp?)/tarC1^n1?pos_/tarC2^n2?pos_/tarC3^n3?pos_/tarC4^n4?pos_/tarC5^n5?{>1}=
        ((tarC3^(1 - n3)*rat(-1, 2))/(tarC1^n1*tarC2^n2*tarC4^n4*tarC5^n5*[pp-3*tmm]) +
        (tarC1^(1 - n1)*rat(-1, 2))/(tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5*[pp-3*tmm]) +
        (tarC4^(1 - n4)*rat(1, 2))/(tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5*[pp-3*tmm]) +
        (tarC2^(1 - n2)*rat(1, 2))/(tarC1^n1*tarC3^n3*tarC4^n4*tarC5^n5*[pp-3*tmm]) +
        (tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC5^(1 - n5)*rat(-n1, 2*(-1 + n5)))/(tarC2^n2*tarC4^n4*[pp-3*tmm]) +
        (tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC5^(1 - n5)*(rat(-2*n1, -1 + n5)/[pp-3*tmm] +
        (pp*rat(n1, 2*(-1 + n5)))/(tmm*[pp-3*tmm])))/(tarC3^n3*tarC4^n4) +
        (tarC1^(-1 - n1)*tarC5^(2 - n5)*((pp*rat(-n1, 2*(-1 + n5)))/(tmm*[pp-3*tmm]) +
        rat(2*n1, -1 + n5)/[pp-3*tmm]))/(tarC2^n2*tarC3^n3*tarC4^n4) +
        (tarC1^(1 - n1)*tarC3^(-1 - n3)*tarC5^(1 - n5)*rat(-n3, 2*(-1 + n5)))/(tarC2^n2*tarC4^n4*[pp-3*tmm]) +
        (tarC3^(-1 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5)*(rat(-2*n3, -1 + n5)/[pp-3*tmm] +
        (pp*rat(n3, 2*(-1 + n5)))/(tmm*[pp-3*tmm])))/(tarC1^n1*tarC2^n2) +
        (tarC3^(-1 - n3)*tarC5^(2 - n5)*((pp*rat(-n3, 2*(-1 + n5)))/(tmm*[pp-3*tmm]) +
        rat(2*n3, -1 + n5)/[pp-3*tmm]))/(tarC1^n1*tarC2^n2*tarC4^n4) +
        (tarC5^(1 - n5)*((pp*rat(2 + d + dp - n1 - n3 - 2*n5, 2*(-1 + n5)))/(tmm*[pp-3*tmm]) +
        rat(-6 - 2*(d + dp) + n1 + n3 + 6*n5, 2*(-1 + n5))/[pp-3*tmm]))/
        (tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4))*sdim(dp);
        
#endprocedure        
*--#] redF :

*--#[ redT2 :
#procedure redT2
        Symmetrize T2 4,3,2;
* eq.91-4
        id T2(dp?,n1?{>1},n2?pos_,n3?pos_)=
        (rat(-2*n2, 3*(-1 + n1))*T2(dp, -2 + n1, 1 + n2, n3))/tmm + (rat(3 + d - 3*n1, 3*(-1 + n1))*T2(dp, -1 + n1, n2, n3))/tmm + 
        (rat(2*n2, 3*(-1 + n1))*T2(dp, -1 + n1, 1 + n2, -1 + n3))/tmm + (rat(1, 3)*T2(dp, n1, -1 + n2, n3))/tmm + 
        (rat(-1, 3)*T2(dp, n1, n2, -1 + n3))/tmm;

*       Master integral        
        id T2(0,1,1,1) = T2x111;
#endprocedure
*--#] redT2 :


*--#[ drrG :
#procedure drrG(topo)

*       Now we need symmetrize only once, because structure (dp,n1 < n2) is kept
        Symmetrize G 2,3;

        #$repcount = 1;        
        #do irep=1,1
                #$irep = 1;                
                if(count(int`topo',1));
               
* eq.95
                id,ifmatch->sortme G(dp?pos_,n1?{>1},n2?pos_) = 
                G(-2 + dp, -1 + n1, -1 + n2)*rat(-1, 2*(-1 + n1))/pp + 
                G(-2 + dp, -1 + n1, n2)*rat(-1, 2*(-1 + n1)) + 
                G(-2 + dp, -2 + n1, n2)*rat(1, 2*(-1 + n1))/pp;
                
                endif;
                
                goto endrec;                
                la sortme;
                $irep = 0;
                la endrec;
                
                ModuleOption,minimum,$irep;
                .sort:drrG-1-`$repcount++';
                #redefine irep "`$irep'"
        #enddo

        #$repcount = 1;        
        #do irep=1,1
                #$irep = 1;                
                if(count(int`topo',1));
                
* eq.96
                id,ifmatch->sortme G(dp?pos_,1,n2?pos_) = 
                T1(-2 + dp, n2)*rat(-1, 2*(-2 + d + dp - n2)) + 
                G(-2 + dp, 1, -1 + n2)*rat(-1, 2*(-2 + d + dp - n2)) + 
                G(-2 + dp, 1, n2)*(tmm*rat(-2, -2 + d + dp - n2) + pp*rat(1, 2*(-2 + d + dp - n2)));                
                endif;
                
                goto endrec;                
                la sortme;
                $irep = 0;
                la endrec;
                
                ModuleOption,minimum,$irep;
                .sort:drrG-2-`$repcount++';
                #redefine irep "`$irep'"
        #enddo
#endprocedure
*--#] drrG :

*--#[ redG :
#procedure redG(topo)
        #$repcount = 1;        
        #do irep=1,1
                #$irep = 1;                
                if(count(int`topo',1));
                
                Symmetrize G 3,2;
* eq.94
                id,ifmatch->sortme G(dp?,n1?{>1},n2?pos_)=
                G(dp, n1, -1 + n2)*((tmm*rat(-2, 1))/(pp*[pp-4*tmm]) + rat(1, 1)/[pp-4*tmm]) + 
                (tmm*G(dp, -2 + n1, 1 + n2)*rat(2*n2, -1 + n1))/(pp*[pp-4*tmm]) + 
                G(dp, -1 + n1, n2)*((tmm*rat(2*(-1 + n1 - n2), -1 + n1))/(pp*[pp-4*tmm]) + 
                rat(-1 - d + dp + n1 + 2*n2, -1 + n1)/[pp-4*tmm]);
                endif;
               
                goto endrec;                
                la sortme;
                $irep = 0;
                la endrec;

                b pp,[pp-4*tmm];
                ModuleOption,minimum,$irep;                
                .sort:redG-`$repcount++';
                #redefine irep "`$irep'"
                Keep Brackets;
                ratio [pp-4*tmm],pp,mden;        
                Multiply replace_(mden,4*tmm);
                .sort                
        #enddo
#endprocedure
*--#] redG :

*--#[ drrT1 :
#procedure drrT1
* eq.98        
        repeat id T1(dp?pos_,n1?pos_) =        
        tmm*rat(-2, d + dp - 2*n1)*T1(-2 + dp, n1);
#endprocedure
*--#] drrT1 :

*--#[ redT1 :
#procedure redT1
* eq.97        
        repeat id T1(dp?,n1?{>1}) =        
        (rat(2 + d + dp - 2*n1, 2*(-1 + n1))*T1(dp, -1 + n1))/tmm;
        
        id T1(0,1) = T1x1;
#endprocedure
*--#] redT1 :

*--#[ partfrac :
#procedure partfrac
* Do partial fractioning for propagators connecting 
* two-loop diagram with one-loop
        b pp,[pp-tmm],[pp-3*tmm],[pp-4*tmm],[pp-9*tmm];
        .sort
        Keep Brackets;        
*       [pp]        
        #do ms={tmm,3*tmm,4*tmm,9*tmm}
                ratio [pp-`ms'],pp,mden;        
                Multiply replace_(mden,`ms');
        #enddo
*       [pp-ms]        
        #do ms={3*tmm,4*tmm,9*tmm}
                ratio [pp-`ms'],[pp-tmm],mden;        
                Multiply replace_(mden,`ms'-tmm);
        #enddo
*       [pp-3ms]                
        #do ms={4*tmm,9*tmm}
                ratio [pp-`ms'],[pp-3*tmm],mden;        
                Multiply replace_(mden,`ms'-3*tmm);
        #enddo
*       [pp-4ms]                
        ratio [pp-9*tmm],[pp-4*tmm],mden;        
        Multiply replace_(mden,5*tmm);
        .sort        
#endprocedure        
*--#] partfrac :


*--#[ subpoch :
#procedure subpoch
*         Pochhammer(n,x)    = Pochhammer[x,n]
*         PochhammerINV(n,x) = 1/Pochhammer[x,n]        
        
        repeat id Pochhammer(n?pos_,x?) = Pochhammer(n-1,x)*num(n-1+x);
        repeat id Pochhammer(n?neg_,x?) = Pochhammer(n+1,x)*den(n+x);
        repeat id PochhammerINV(n?pos_,x?) = PochhammerINV(n-1,x)*den(n-1+x);
        repeat id PochhammerINV(n?neg_,x?) = PochhammerINV(n+1,x)*num(n+x);
        id	  Pochhammer(0,x?) = 1;
        id	  PochhammerINV(0,x?) = 1;
        id	num(x?)*den(x?) = 1;
        id	den(x?number_) = 1/x;
        id	num(x?number_) = x;
*Print +f "<3> %t";
        id	num(x?) = rat(x,1);
*Print +f "<4> %t";
        id	den(x?) = rat(1,x);
#endprocedure
*--#] subpoch :


*--#[ tensG :
#procedure tensG(D1,D2,K,Q)
* Tensor reduction of G with tensor structure in numerator
* Loop momentum is K
*         
*    D1 = K^2 - tm^2        
*    D2 = (K-Q)^2 - tm^2
*         
*    /
*    |       (K.PX)^n
*    | ------------------ dP1
*    |   (D1)^n1 * (D2)^n2
*    /     

        id `D1'^n1?pos_ =  (`K'.`K' - tmm)^n1;        
        id `D2'^n2?pos_ =  (`K'.`K' - 2*`K'.`Q' + `Q'.`Q' - tmm)^n2;        

*       <1> reduce all squares of loop momenta in numerator

        repeat id `K'.`Q'^n?pos_/`D2'^n2?pos_ = `K'.`Q'^(n-1)*(-`D2' + `K'.`K' + `Q'.`Q' - tmm)/2/`D2'^n2;

        repeat id `K'.`K'^n?pos_/`D1'^n1?pos_/`D2'^n2?pos0_ = `K'.`K'^(n-1)*(`D1' + tmm)/`D1'^n1/`D2'^n2;
        repeat id `K'.`K'^n?pos_/`D2'^n2?pos_ = `K'.`K'^(n-1)*(`D2' + 2*`K'.`Q' - `Q'.`Q' + tmm)/`D2'^n2;


        if((count(`D1',-1) > 0) && (count(`D2',-1) > 0));
*         
*       Self-energy        
*         
        totensor,nosquare,`K',ftensor;
*       Davydychev, eq.11, momenta flow differen G[n1,n2]->G[n2,n1]        
        id ftensor(?a)/`D1'^n2?pos_/`D2'^n1?pos_=
        sum_(j,0,integer_(nargs_(?a)/2),sign_(j)/2^j*distrib_(1,2*j,dd,ftensor,?a)*fac_(n1 + nargs_(?a) - 2*j- 1)*invfac_(n1 - 1)*G(2*(nargs_(?a)-j),n1+nargs_(?a)-2*j,n2));
        
        id 1/`D1'^n1?pos_/`D2'^n2?pos_ = G(0, n1,n2);
*       Change external momenta sign according to Davydychev        
        id ftensor(?a) = sign_(nargs_(?a))*ftensor(?a);
        
        tovector, ftensor, `Q';
        id dd(?a)=dd_(?a);

        elseif ((count(`D1',-1) > 0) && (count(`D2',1) >= 0));
*         
*       Tadpole D1        
*         
        totensor,nosquare,`K',ftensor;
        id ftensor(?a) = [sqrt(x)]^nargs_(?a)*ftensor(?a);
        id [sqrt(x)]^n?odd_ = 0;
        id,many, [sqrt(x)]*[sqrt(x)] = [x];
        
        id ftensor(?a) = dd_(?a);
*         .sort
        if ( count([x],1) != 0 );
        id [x]^s? =  rat(1-(2-d/2),1)*PochhammerINV(s+2,d/2-2)*rat(d-4,2)*`K'.`K'^(s)/2^s;
        endif;        
        id `K'.`K'^n?pos_ = (`D1' + tmm)^n;
        if(count(`D1',1) >= 0) Discard;        
        id 1/`D1'^n?pos_ = T1(0,n);


        elseif((count(`D1',1) >= 0) && (count(`D2',-1) > 0));
*         
*       Tadpole D2        
*         
        Multiply replace_(`K', [k-q] + `Q');
        totensor,nosquare,[k-q],ftensor;
        id ftensor(?a) = [sqrt(x)]^nargs_(?a)*ftensor(?a);
        id [sqrt(x)]^n?odd_ = 0;
        id,many, [sqrt(x)]*[sqrt(x)] = [x];
        
        id ftensor(?a) = dd_(?a);
*         .sort
        if ( count([x],1) != 0 );
        id [x]^s? =  rat(1-(2-d/2),1)*PochhammerINV(s+2,d/2-2)*rat(d-4,2)*[k-q].[k-q]^(s)/2^s;
        endif;        
        id [k-q].[k-q] = `D2' + tmm;
        
        Multiply replace_([k-q],`K'-`Q');
        if(count(`D2',1) >= 0) Discard;        
        id 1/`D2'^n?pos_ = T1(0,n);

        else;
        Discard;        
        endif;
#endprocedure
*--#] tensG :


*--#[ tensT1 :
#procedure tensT1(D1,K,MUL)
* Tensor reduction of T1 with tensor structure in numerator
* Loop momentum is K
*         
*    D1 = K^2 - tm^2        
*         
*    /
*    |   (K.PX)^n
*    | ----------- dK
*    |   (D1)^n1 
*    /     

        id `D1'^n1?pos_ =  (`K'.`K' - tmm)^n1;        

*       <1> reduce all squares of loop momenta in numerator

        repeat id `K'.`K'^n?pos_/`D1'^n1?pos_ = `K'.`K'^(n-1)*(`D1' + tmm)/`D1'^n1;

        if (count(`D1',-1) > 0);
*         
*       Tadpole D1        
*         
        totensor,nosquare,`K',ftensor;
        id ftensor(?a) = [sqrt(x)]^nargs_(?a)*ftensor(?a);
        id [sqrt(x)]^n?odd_ = 0;
        id,many, [sqrt(x)]*[sqrt(x)] = [x];
        
        id ftensor(?a) = dd_(?a);
*         .sort
        if ( count([x],1) != 0 );
        id [x]^s? =  rat(1-(2-d/2),1)*PochhammerINV(s+2,d/2-2)*rat(d-4,2)*`K'.`K'^(s)/2^s;
        endif;        
        id `K'.`K'^n?pos_ = (`D1' + tmm)^n;
        if(count(`D1',1) >= 0) Discard;        
        id 1/`D1'^n?pos_ = T1(0,n);
        else;
        Multiply `MUL';        
        endif;
#endprocedure
*--#] tensT1 :

*
****
*   FMFT part
****

*
* scalar products to denominators
*

*--#[ sp2den :
#procedure sp2den(TOPO)

if(count(int`TOPO',1));        
id p1.p1 = 1 + d1;
id p2.p2 = 1 + d2;
id p3.p3 = 1 + d3;
id p4.p4 = 1 + d4;
id p5.p5 = 1 + d5;
id p6.p6 = 1 + d6;
id p7.p7 = 1 + d7;
id p8.p8 = 1 + d8;
id p9.p9 = 1 + d9;
id p10.p10 = 1 + d10;
endif;
.sort:sp2den-pp;
        
if(count(int`TOPO',1));
id p1.p2 = 1/2 + d1/2 + d2/2 - d8/2;
id p1.p3 = 1/2 + d1/2 + d3/2 - d9/2;
id p1.p4 = 1/2 + d1/2 + d4/2 - d5/2;
id p1.p5 = 1/2 + d1/2 - d4/2 + d5/2;
id p1.p6 = d2/2 - d4/2 + d5/2 - d8/2;
id p1.p7 = d3/2 - d4/2 + d5/2 - d9/2;
id p1.p8 = 1/2 + d1/2 - d2/2 + d8/2;
id p1.p9 = 1/2 + d1/2 - d3/2 + d9/2;
id p1.p10 = -d2/2 - d3/2 + d8/2 + d9/2;
endif;
.sort:sp2den-p1;

if(count(int`TOPO',1));
id p2.p3 = d1/2 - d8/2 - d9/2 + d10/2;
id p2.p4 = 1/2 + d2/2 + d4/2 - d6/2;
id p2.p5 = d1/2 - d4/2 + d6/2 - d8/2;
id p2.p6 = 1/2 + d2/2 - d4/2 + d6/2;
id p2.p7 = -1/2 + d1/2 - d2/2 - d4/2 + d6/2 - d8/2 - d9/2 + d10/2;
id p2.p8 = -1/2 + d1/2 - d2/2 - d8/2;
id p2.p9 = 1/2 + d2/2 + d9/2 - d10/2;
id p2.p10 = -1/2 - d2/2 + d9/2 - d10/2;
endif;
.sort:sp2den-p2;

if(count(int`TOPO',1));
id p3.p4 = 1/2 + d3/2 + d4/2 - d7/2;
id p3.p5 = d1/2 - d4/2 + d7/2 - d9/2;
id p3.p6 = -1/2 + d1/2 - d3/2 - d4/2 + d7/2 - d8/2 - d9/2 + d10/2;
id p3.p7 = 1/2 + d3/2 - d4/2 + d7/2;
id p3.p8 = 1/2 + d3/2 + d8/2 - d10/2;
id p3.p9 = -1/2 + d1/2 - d3/2 - d9/2;
id p3.p10 = -1/2 - d3/2 + d8/2 - d10/2;
endif;
.sort:sp2den-p3;

if(count(int`TOPO',1));
id p4.p5 = -1/2 + d1/2 - d4/2 - d5/2;
id p4.p6 = -1/2 + d2/2 - d4/2 - d6/2;
id p4.p7 = -1/2 + d3/2 - d4/2 - d7/2;
id p4.p8 = d1/2 - d2/2 - d5/2 + d6/2;
id p4.p9 = d1/2 - d3/2 - d5/2 + d7/2;
id p4.p10 = -1/2 + d1/2 - d2/2 - d3/2 - d4/2 - d5/2 + d6/2 + d7/2;
endif;
.sort:sp2den-p4;

if(count(int`TOPO',1));
id p5.p6 = 1/2 + d5/2 + d6/2 - d8/2;
id p5.p7 = 1/2 + d5/2 + d7/2 - d9/2;
id p5.p8 = 1/2 + d5/2 - d6/2 + d8/2;
id p5.p9 = 1/2 + d5/2 - d7/2 + d9/2;
id p5.p10 = 1/2 - d1/2 + d4/2 + d5/2 - d6/2 - d7/2 + d8/2 + d9/2;
endif;
.sort:sp2den-p5;

if(count(int`TOPO',1));
id p6.p7 = d1/2 - d2/2 - d3/2 + d6/2 + d7/2 - d8/2 - d9/2 + d10/2;
id p6.p8 = -1/2 + d5/2 - d6/2 - d8/2;
id p6.p9 = 1/2 - d1/2 + d2/2 + d3/2 + d5/2 - d7/2 + d9/2 - d10/2;
id p6.p10 = -d1/2 + d3/2 + d4/2 + d5/2 - d6/2 - d7/2 + d9/2 - d10/2;
endif;
.sort:sp2den-p6;

if(count(int`TOPO',1));
id p7.p8 = 1/2 - d1/2 + d2/2 + d3/2 + d5/2 - d6/2 + d8/2 - d10/2;
id p7.p9 = -1/2 + d5/2 - d7/2 - d9/2;
id p7.p10 = -d1/2 + d2/2 + d4/2 + d5/2 - d6/2 - d7/2 + d8/2 - d10/2;
endif;
.sort:sp2den-p7;

if(count(int`TOPO',1));
id p8.p9 = d1/2 - d2/2 - d3/2 + d10/2;
id p8.p10 = 1/2 - d3/2 + d8/2 + d10/2;
endif;
.sort:sp2den-p8;

if(count(int`TOPO',1));
id p9.p10 = 1/2 - d2/2 + d9/2 + d10/2;
endif;
.sort:sp2den-p9;

#endprocedure
*--#] sp2den :


*
* 
* Mapping topologies
* 
*

*--#[ mapBMW :
#procedure mapBMW(TOPO)        
        if(count(int`TOPO',1));

*       BMW ::       
        if(match(1/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_))
        Multiply intBMW/int`TOPO';
 
*       [0, 1, 0, 1, 1, 1, 1, 1, 1, 1] -> BMW
        if(match(1/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_));
        if ((count(d1,1) == 0) && (count(d3,1) == 0));
        Multiply replace_(d4,d10, d6,d3, d8,d4, d9,d5, d10,d6, d5,d7, d2,d8, d7,d9)*intBMW/int`TOPO';         
        elseif (count(d1,1) == 0);
        Multiply replace_(d4,d10, d3,d2, d6,d3, d8,d4, d9,d5, d10,d6, d5,d7, d2,d8, d7,d9)*intBMW/int`TOPO';         
        elseif (count(d3,1) == 0);
        Multiply replace_(d4,d10, d6,d3, d8,d4, d9,d5, d10,d6, d5,d7, d2,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d7,d9)*intBMW/int`TOPO';         
        else;
        Multiply replace_(d4,d10, d3,d2, d6,d3, d8,d4, d9,d5, d10,d6, d5,d7, d2,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d7,d9)*intBMW/int`TOPO';         
        endif;
        endif;

*       [0, 1, 1, 0, 1, 1, 1, 1, 1, 1] -> BMW
        if(match(1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_));
        if ((count(d1,1) == 0) && (count(d4,1) == 0));
        Multiply replace_(d2,d10, d10,d3, d3,d4, d7,d6, d8,d7, d9,d8, d6,d9)*intBMW/int`TOPO';         
        elseif (count(d1,1) == 0);
        Multiply replace_(d2,d10, d4,d2, d10,d3, d3,d4, d7,d6, d8,d7, d9,d8, d6,d9)*intBMW/int`TOPO';         
        elseif (count(d4,1) == 0);
        Multiply replace_(d2,d10, d10,d3, d3,d4, d7,d6, d8,d7, d9,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d6,d9)*intBMW/int`TOPO';         
        else;
        Multiply replace_(d2,d10, d4,d2, d10,d3, d3,d4, d7,d6, d8,d7, d9,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d6,d9)*intBMW/int`TOPO';         
        endif;        
        endif;

*       [0, 1, 1, 1, 0, 1, 1, 1, 1, 1] -> BMW
        if(match(1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_));
        if ((count(d1,1) == 0) && (count(d5,1) == 0));
        Multiply replace_(d9,d10, d10,d3, d8,d4, d4,d5, d3,d7, d2,d8, d7,d9)*intBMW/int`TOPO';         
        elseif (count(d1,1) == 0);
        Multiply replace_(d9,d10, d5,d2, d10,d3, d8,d4, d4,d5, d3,d7, d2,d8, d7,d9)*intBMW/int`TOPO';         
        elseif (count(d5,1) == 0);
        Multiply replace_(d9,d10, d10,d3, d8,d4, d4,d5, d3,d7, d2,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d7,d9)*intBMW/int`TOPO';         
        else;
        Multiply replace_(d9,d10, d5,d2, d10,d3, d8,d4, d4,d5, d3,d7, d2,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d7,d9)*intBMW/int`TOPO';         
        endif;        
        endif;

*       [0, 1, 1, 1, 1, 0, 1, 1, 1, 1] -> BMW
        if(match(1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_));
        if ((count(d1,1) == 0) && (count(d6,1) == 0));
        Multiply replace_(d2,d10, d9,d3, d5,d4, d3,d5, d8,d6, d10,d8, d4,d9)*intBMW/int`TOPO';         
        elseif (count(d1,1) == 0);
        Multiply replace_(d2,d10, d6,d2, d9,d3, d5,d4, d3,d5, d8,d6, d10,d8, d4,d9)*intBMW/int`TOPO';         
        elseif (count(d6,1) == 0);
        Multiply replace_(d2,d10, d9,d3, d5,d4, d3,d5, d8,d6, d10,d8, d4,d9, d1,(1 - d1 + d3 + d4 + d5 - d7 + d9))*intBMW/int`TOPO';         
        else;
        Multiply replace_(d2,d10, d6,d2, d9,d3, d5,d4, d3,d5, d8,d6, d10,d8, d4,d9, d1,(1 - d1 + d3 + d4 + d5 - d7 + d9))*intBMW/int`TOPO';         
        endif;        
        endif;

*       [0, 1, 1, 1, 1, 1, 0, 1, 1, 1] -> BMW
        if(match(1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_));
        if ((count(d1,1) == 0) && (count(d7,1) == 0));
        Multiply replace_(d3,d10, d8,d3, d5,d4, d2,d5, d9,d6, d6,d7, d10,d8, d4,d9)*intBMW/int`TOPO';         
        elseif (count(d1,1) == 0);
        Multiply replace_(d3,d10, d7,d2, d8,d3, d5,d4, d2,d5, d9,d6, d6,d7, d10,d8, d4,d9)*intBMW/int`TOPO';         
        elseif (count(d7,1) == 0);
        Multiply replace_(d3,d10, d8,d3, d5,d4, d2,d5, d9,d6, d6,d7, d10,d8, d4,d9, d1,(1 - d1 + d3 + d4 + d5 - d7 + d9))*intBMW/int`TOPO';         
        else;
        Multiply replace_(d3,d10, d7,d2, d8,d3, d5,d4, d2,d5, d9,d6, d6,d7, d10,d8, d4,d9, d1,(1 - d1 + d3 + d4 + d5 - d7 + d9))*intBMW/int`TOPO';         
        endif;        
        endif;

*       [0, 1, 1, 1, 1, 1, 1, 0, 1, 1] -> BMW
        if(match(1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_));
        if ((count(d1,1) == 0) && (count(d8,1) == 0));
        Multiply replace_(d3,d10, d7,d3, d5,d4, d2,d5, d9,d7, d4,d8, d10,d9)*intBMW/int`TOPO';         
        elseif (count(d1,1) == 0);
        Multiply replace_(d3,d10, d8,d2, d7,d3, d5,d4, d2,d5, d9,d7, d4,d8, d10,d9)*intBMW/int`TOPO';         
        elseif (count(d8,1) == 0);
        Multiply replace_(d3,d10, d7,d3, d5,d4, d2,d5, d9,d7, d4,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d10,d9)*intBMW/int`TOPO';         
        else;
        Multiply replace_(d3,d10, d8,d2, d7,d3, d5,d4, d2,d5, d9,d7, d4,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d10,d9)*intBMW/int`TOPO';         
        endif;        
        endif;

*       [0, 1, 1, 1, 1, 1, 1, 1, 0, 1] -> BMW
        if(match(1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_));
        if ((count(d1,1) == 0) && (count(d9,1) == 0));
        Multiply replace_(d2,d10, d6,d3, d5,d4, d3,d5, d7,d6, d8,d7, d4,d8, d10,d9)*intBMW/int`TOPO';         
        elseif (count(d1,1) == 0);
        Multiply replace_(d2,d10, d9,d2, d6,d3, d5,d4, d3,d5, d7,d6, d8,d7, d4,d8, d10,d9)*intBMW/int`TOPO';         
        elseif (count(d9,1) == 0);
        Multiply replace_(d2,d10, d6,d3, d5,d4, d3,d5, d7,d6, d8,d7, d4,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d10,d9)*intBMW/int`TOPO';         
        else;
        Multiply replace_(d2,d10, d9,d2, d6,d3, d5,d4, d3,d5, d7,d6, d8,d7, d4,d8, d1,(1 - d1 + d2 + d4 + d5 - d6 + d8), d10,d9)*intBMW/int`TOPO';         
        endif;        
        endif;

*       [0, 1, 1, 1, 1, 1, 1, 1, 1, 0] -> BMW
        if(match(1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_));
        if ((count(d1,1) == 0) && (count(d10,1) == 0));
        Multiply replace_(d2,d10, d4,d3, d3,d4, d8,d6, d6,d8)*intBMW/int`TOPO';
        elseif (count(d1,1) == 0);
        Multiply replace_(d2,d10, d10,d2, d4,d3, d3,d4, d8,d6, d6,d8)*intBMW/int`TOPO';
        elseif (count(d10,1) == 0);
        Multiply replace_(d2,d10, d4,d3, d3,d4, d8,d6, d6,d8, d1,(1 - d1 + d3 + d4 + d5 - d7 + d9))*intBMW/int`TOPO';
        else;
        Multiply replace_(d2,d10, d10,d2, d4,d3, d3,d4, d8,d6, d6,d8, d1,(1 - d1 + d3 + d4 + d5 - d7 + d9))*intBMW/int`TOPO';         
        endif;        
        endif;

* PL
*       [1, 1, 1, 0, 1, 1, 1, 1, 1, 0] -> BMW
        if(match(1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_));
        if ((count(d4,1) == 0) && (count(d10,1) == 0));
        Multiply replace_(d2,d10, d8,d3, d6,d4, d9,d5, d3,d6, d5,d7, d1,d8, d7,d9)*intBMW/int`TOPO';         
        elseif (count(d4,1) == 0);
        Multiply replace_(d2,d10, d8,d3, d6,d4, d9,d5, d3,d6, d5,d7, d1,d8, d10,(d1 + d10 - d4 + d6 + d7 - d8 - d9), d7,d9)*intBMW/int`TOPO';         
        elseif (count(d10,1) == 0);
        Multiply replace_(d2,d10, d8,d3, d6,d4, d9,d5, d3,d6, d5,d7, d4,(d1 + d10 - d2 - d3 - d5 + d6 + d7), d1,d8, d7,d9)*intBMW/int`TOPO';         
        else;
        Multiply replace_(d2,d10, d8,d3, d6,d4, d9,d5, d3,d6, d5,d7, d4,(d1 + d10 - d2 - d3 - d5 + d6 + d7), d1,d8, d10,(d1 + d10 - d4 + d6 + d7 - d8 - d9), d7,d9)*intBMW/int`TOPO';         
        endif;        
        endif;

*       [1, 1, 1, 1, 0, 1, 1, 1, 1, 0] -> BMW
        if(match(1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_));
        if ((count(d5,1) == 0) && (count(d10,1) == 0));
        Multiply replace_(d9,d10, d7,d4, d2,d5, d8,d6, d4,d7, d1,d8, d6,d9)*intBMW/int`TOPO';        
        elseif (count(d5,1) == 0);
        Multiply replace_(d9,d10, d7,d4, d2,d5, d8,d6, d4,d7, d1,d8, d10,(d1 + d10 - d4 + d6 + d7 - d8 - d9), d6,d9)*intBMW/int`TOPO';        
        elseif (count(d10,1) == 0);
        Multiply replace_(d9,d10, d7,d4, d2,d5, d8,d6, d4,d7, d5,(d1 + d10 - d2 - d3 - d5 + d6 + d7), d1,d8, d6,d9)*intBMW/int`TOPO';        
        else;
        Multiply replace_(d9,d10, d7,d4, d2,d5, d8,d6, d4,d7, d5,(d1 + d10 - d2 - d3 - d5 + d6 + d7), d1,d8, d10,(d1 + d10 - d4 + d6 + d7 - d8 - d9), d6,d9)*intBMW/int`TOPO';        
        endif;        
        endif;
        
        endif;
#endprocedure
*--#] mapBMW :


*--#[ mapFG :
#procedure mapFG(TOPO)
id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d1^n4*d10^n6*d4^n7*d5^n3*d7^n8*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1))/(d4^n6*d5^n3*d7^n8*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 + 4*tk1.tk4 + 4*tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d4^n6*d7^n7*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 + 4*tk1.tk4 + 4*tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4))/(d4^n6*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1))/(d3^n6*d4^n5*d5^n3*d7^n8*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n2*d3^n9*d4^n5*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d1^n3*d3^n8*d4^n5*d7^n6*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + pp + 4*tk1.tk1 + 4*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d1^n3*d3^n8*d4^n5*d7^n6*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n7)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n5))/(d4^n4*d5^n3*d7^n8*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n8)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d3^n3*d4^n4*d7^n7*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d3^n2*d4^n4*d7^n6*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d1^n10*d3^n2*d4^n4*d7^n6*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n5)*
num(-1 + pp + 4*tk1.tk1 + 4*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n9*d3^n2*d4^n4*d7^n6*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + 4*tk2.tk2 + 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d3^n2*d4^n4*d7^n6*d8^n7*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n1*d3^n2*d4^n4*d5^n5*d6^n9*d7^n6*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 4*tk1.tk2 - 2*tk1.tk4 + 4*tk2.tk2 + 4*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d10^n3*d3^n1*d4^n4*d7^n5*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 2*tk1.tk4 + 4*tk2.tk2 - 4*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d3^n1*d4^n4*d7^n5*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n10)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n3*d3^n1*d4^n4*d5^n7*d7^n5*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d3^n1*d4^n4*d7^n5*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n2*d3^n1*d4^n4*d5^n6*d7^n5*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n1*d4^n4*d6^n3*d7^n5*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n3*d3^n1*d4^n4*d5^n7*d7^n5*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 + 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n1*d4^n4*d7^n5*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 4*tk1.tk2 - 2*tk1.tk4 + 4*tk2.tk2 + 4*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d3^n4*d4^n3*d7^n7*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + 4*tk2.tk2 + 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d3^n4*d4^n3*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n7))/(d4^n3*d7^n6*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d4^n3*d7^n6*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n8)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + 4*tk1.tk1 + 4*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d1^n5*d4^n3*d7^n6*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk1.tk1 - 4*tk1.tk2 + 2*tk1.tk4 + 4*tk2.tk2 - 4*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d4^n3*d7^n6*d8^n7*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n9))/(d4^n3*d5^n5*d7^n6*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d4^n3*d7^n5*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d1^n4*d4^n3*d5^n7*d7^n5*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d4^n3*d7^n5*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d4^n3*d5^n6*d7^n5*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d1^n4*d4^n3*d5^n7*d7^n5*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d1^n1*d4^n3*d5^n9*d7^n5*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + 4*tk2.tk2 + 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d4^n3*d7^n4*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 + 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d4^n3*d7^n4*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n8)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2))/(d4^n3*d7^n4*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 - 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n3*d7^n4*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n9))/(d4^n3*d7^n4*d8^n6*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n8)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1))/(d3^n7*d4^n3*d7^n4*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n9))/(d4^n3*d7^n4*d8^n5*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n9))/(d4^n3*d5^n2*d7^n4*d8^n5*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d4^n2*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n6))/(d4^n2*d5^n5*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 4*tk1.tk2 - 2*tk1.tk4 + 4*tk2.tk2 + 4*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d3^n4*d4^n2*d7^n6*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + 4*tk2.tk2 + 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3))/(d3^n4*d4^n2*d7^n6*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n8)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + 4*tk1.tk1 + 4*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d3^n4*d4^n2*d7^n6*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n8))/(d3^n4*d4^n2*d7^n6*d8^n7*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk1.tk1 - 4*tk1.tk2 + 2*tk1.tk4 + 4*tk2.tk2 - 4*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d3^n4*d4^n2*d7^n6*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d4^n2*d7^n5*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d4^n2*d7^n5*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d1^n4*d4^n2*d5^n6*d7^n5*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d4^n2*d7^n5*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d4^n2*d5^n7*d7^n5*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n10)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d4^n2*d7^n5*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + 4*tk2.tk2 + 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d4^n2*d7^n4*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 + 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3))/(d4^n2*d7^n4*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n8)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3))/(d4^n2*d7^n4*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n2*d7^n4*d8^n7*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 - 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 + 4*tk4.tp, -n9))/(d4^n2*d7^n4*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n2*d7^n4*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n2*d7^n4*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n2*d5^n3*d7^n4*d8^n5*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d1^n1*d4^n2*d5^n8*d7^n3*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d4^n2*d5^n1*d7^n3*d8^n7*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n9)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 - 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 + 4*tk4.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d4^n2*d7^n3*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d1^n1*d4^n2*d5^n8*d6^n5*d7^n3*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d4^n2*d5^n1*d7^n3*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n9)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d4^n2*d7^n3*d8^n6*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n1*d4^n2*d5^n8*d7^n3*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d4^n2*d5^n1*d7^n3*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n9)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 + 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d4^n2*d7^n3*d8^n5*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n6)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n2*d5^n4*d7^n3*d8^n5*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n4)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 + 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + pp + 4*tk2.tk2 + 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n2*d7^n3*d8^n5*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n1*d2^n5*d4^n2*d5^n8*d7^n3*d8^n4*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d4^n2*d5^n1*d6^n5*d7^n3*d8^n4*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n9)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d4^n2*d7^n3*d8^n4*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + tk1.tk1 - 4*tk1.tk2 - 2*tk1.tk4 + 4*tk2.tk2 + 4*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + 4*tk2.tk2 + 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d4^n2*d7^n3*d8^n4*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n4*d4^n1*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d4^n1*d7^n6*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n8))/(d4^n1*d7^n6*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n9))/(d1^n4*d4^n1*d5^n5*d7^n6*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d3^n4*d4^n1*d7^n5*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n10)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d3^n4*d4^n1*d7^n5*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n10)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d3^n4*d4^n1*d7^n5*d8^n6*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d1^n3*d3^n4*d4^n1*d5^n9*d7^n5*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n1*d7^n4*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n9))/(d4^n1*d7^n4*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n1*d7^n4*d8^n6*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n8)*
num(-1 + 9*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 6*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 6*tk2.tp + tk4.tk4 + 6*tk4.tp, -n10))/(d4^n1*d5^n3*d7^n4*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 - 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 + 4*tk4.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n5))/(d4^n1*d5^n10*d7^n3*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 + 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d4^n1*d5^n10*d7^n3*d8^n6*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n8)*
num(-1 + 9*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 6*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 6*tk2.tp + tk4.tk4 + 6*tk4.tp, -n10))/(d4^n1*d5^n4*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d4^n1*d5^n10*d7^n3*d8^n5*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n1*d7^n3*d8^n5*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d4^n1*d5^n10*d7^n3*d8^n4*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n1*d7^n3*d8^n4*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 + 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d4^n1*d5^n10*d7^n2*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + 4*tk2.tk2 - 4*tk2.tk4 - 8*tk2.tp + tk4.tk4 + 4*tk4.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n5))/(d4^n1*d5^n10*d7^n2*d8^n6*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n5)*
num(-1 + 9*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 6*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 6*tk2.tp + tk4.tk4 + 6*tk4.tp, -n10))/(d4^n1*d7^n2*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d4^n1*d5^n10*d7^n2*d8^n5*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n1*d7^n2*d8^n5*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d4^n1*d5^n10*d7^n2*d8^n4*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n1*d7^n2*d8^n4*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + 4*tk2.tk2 + 4*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n1*d7^n2*d8^n3*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n1*d5^n4*d6^n7*d7^n2*d8^n3*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk2 - 4*tk1.tp + 4*tk2.tk2 - 8*tk2.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n1*d7^n2*d8^n3*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + pp + 4*tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n1*d7^n2*d8^n3*d9^n4);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n10*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n6)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n8*d5^n5*d7^n7*d8^n10*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n7)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n9*d5^n5*d7^n6*d8^n10*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d1^n8*d4^n6*d5^n5*d7^n7*d8^n10*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d1^n9*d4^n7*d5^n5*d7^n6*d8^n10*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1))/(d4^n9*d6^n4*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1))/(d4^n8*d6^n4*d7^n6*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk4 - 4*tk1.tp + 4*tk4.tk4 - 8*tk4.tp, -n1)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d4^n7*d6^n4*d7^n6*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + 4*pp + tk1.tk1 + 4*tk1.tk4 - 4*tk1.tp + 4*tk4.tk4 - 8*tk4.tp, -n1)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d4^n6*d6^n4*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3))/(d2^n6*d3^n4*d4^n5*d6^n8*d7^n1*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n1*d10^n8*d3^n3*d4^n4*d5^n5*d7^n7*d8^n10*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n2*d2^n1*d3^n9*d4^n5*d6^n4*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d1^n3*d2^n1*d3^n8*d4^n5*d6^n4*d7^n6*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d1^n1*d10^n9*d3^n2*d4^n4*d5^n5*d7^n6*d8^n10*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n7))/(d1^n1*d3^n2*d4^n4*d5^n5*d7^n6*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + 4*tk1.tk1 - 4*tk1.tk2 - 4*tk1.tp + tk2.tk2 + 2*tk2.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d10^n3*d2^n8*d3^n4*d4^n5*d6^n6*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n8*d7^n7*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n1*d3^n2*d4^n4*d5^n5*d6^n9*d7^n6*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d1^n4*d10^n6*d4^n7*d5^n3*d7^n8*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1))/(d4^n6*d5^n3*d7^n8*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d1^n8*d3^n3*d4^n6*d5^n5*d6^n7*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1))/(d3^n5*d4^n6*d5^n3*d7^n8*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2))/(d4^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1))/(d3^n6*d4^n5*d5^n3*d7^n8*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d4^n3*d5^n5*d7^n7*d8^n10*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4))/(d4^n5*d5^n3*d7^n8*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n6))/(d4^n3*d5^n5*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4))/(d2^n1*d4^n9*d6^n3*d7^n5*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n7))/(d4^n3*d5^n5*d7^n6*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n2*d2^n4*d3^n9*d4^n7*d6^n3*d7^n5*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n2))/(d4^n3*d5^n5*d6^n8*d7^n7*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n9))/(d4^n3*d5^n5*d7^n6*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n7)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n5))/(d4^n4*d5^n3*d7^n8*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d1^n1*d4^n9*d5^n3*d6^n8*d7^n4*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1))/(d1^n10*d4^n8*d5^n3*d7^n4*d8^n9*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n8)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n5*d2^n7*d4^n3*d6^n4*d7^n6*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n7)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d4^n4*d5^n3*d7^n8*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n5)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1))/(d3^n2*d4^n6*d5^n3*d7^n4*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d10^n10*d3^n2*d4^n6*d5^n3*d7^n4*d8^n9*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n2*d4^n6*d5^n3*d6^n10*d7^n4*d8^n8*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n6*d2^n8*d4^n10*d6^n3*d7^n4*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n7*d4^n4*d5^n3*d7^n8*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n8*d6^n3*d7^n4*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d3^n1*d4^n5*d5^n3*d6^n8*d7^n4*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n10))/(d3^n1*d4^n5*d5^n3*d7^n4*d8^n8*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n1*d4^n4*d6^n3*d7^n5*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n2)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d4^n6*d6^n3*d7^n4*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n9))/(d4^n3*d5^n5*d7^n6*d8^n4*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n8)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d3^n1*d4^n5*d5^n3*d7^n4*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1))/(d4^n7*d5^n2*d7^n9*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d1^n4*d10^n7*d4^n6*d5^n2*d7^n9*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1))/(d1^n4*d4^n6*d5^n2*d7^n9*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n9*d3^n2*d4^n6*d6^n7*d7^n4*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3))/(d4^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1))/(d3^n7*d4^n5*d5^n2*d7^n9*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4))/(d2^n1*d4^n8*d6^n2*d7^n5*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n6))/(d4^n2*d5^n5*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4))/(d4^n5*d5^n2*d7^n9*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n7))/(d4^n2*d5^n5*d7^n6*d8^n10*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n7))/(d4^n2*d5^n5*d7^n6*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n3*d2^n4*d3^n8*d4^n6*d6^n2*d7^n5*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n8))/(d4^n2*d5^n5*d7^n7*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n3))/(d4^n2*d5^n5*d6^n9*d7^n6*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n6)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n5))/(d4^n4*d5^n2*d7^n9*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n6)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5))/(d4^n4*d5^n2*d7^n9*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n5*d2^n6*d4^n2*d6^n4*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5)/(d1^n10*d10^n1*d2^n3*d3^n2*d4^n4*d6^n7*d7^n6*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3))/(d1^n1*d4^n8*d5^n2*d6^n9*d7^n4*d8^n10*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n5)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n10))/(d1^n1*d4^n8*d5^n2*d7^n4*d8^n9*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n5)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1))/(d3^n3*d4^n7*d5^n2*d7^n4*d8^n10*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n3*d4^n7*d5^n2*d6^n10*d7^n4*d8^n9*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d10^n10*d3^n3*d4^n7*d5^n2*d7^n4*d8^n8*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d4^n4*d5^n2*d7^n9*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n7*d2^n9*d4^n10*d6^n2*d7^n4*d8^n5*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n9*d6^n2*d7^n4*d8^n5*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n1*d4^n4*d6^n2*d7^n5*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n2*d5^n5*d7^n7*d8^n4*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n3)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d4^n7*d6^n2*d7^n4*d8^n5*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3))/(d3^n1*d4^n5*d5^n2*d6^n9*d7^n4*d8^n10*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n10))/(d3^n1*d4^n5*d5^n2*d7^n4*d8^n9*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d3^n1*d4^n5*d5^n2*d7^n4*d8^n7*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d3^n1*d4^n3*d5^n2*d7^n9*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d3^n1*d4^n2*d5^n3*d7^n8*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d3^n4*d4^n7*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n7))/(d3^n1*d4^n3*d5^n2*d7^n9*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d3^n1*d4^n2*d5^n3*d7^n8*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d2^n4*d4^n6*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + 9*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 6*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 6*tk2.tp + tk4.tk4 - 6*tk4.tp, -n1)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n5))/(d4^n10*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n1*d4^n9*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d2^n1*d4^n8*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d3^n1*d4^n3*d5^n2*d7^n9*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d3^n1*d4^n2*d5^n3*d7^n8*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n4*d4^n5*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d1^n6*d4^n3*d6^n2*d7^n5*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d10^n6*d4^n2*d5^n5*d7^n7*d8^n3*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n1*d3^n2*d4^n3*d5^n9*d6^n5*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n7*d4^n2*d6^n3*d7^n5*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d1^n10*d10^n4*d3^n2*d4^n3*d5^n8*d6^n5*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n7))/(d4^n2*d5^n5*d7^n6*d8^n3*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n1*d4^n5*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d10^n5*d3^n1*d4^n3*d5^n2*d6^n6*d7^n9*d8^n4*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d10^n5*d3^n1*d4^n2*d5^n3*d6^n7*d7^n8*d8^n4*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n4*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n8))/(d4^n2*d5^n3*d7^n4*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n2*d5^n3*d7^n4*d8^n9*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n2*d5^n3*d6^n10*d7^n4*d8^n8*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n9))/(d4^n3*d5^n2*d7^n4*d8^n10*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n3*d5^n2*d6^n10*d7^n4*d8^n9*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n9))/(d4^n3*d5^n2*d7^n4*d8^n8*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n9*d2^n7*d3^n2*d4^n3*d5^n1*d6^n4*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d10^n8*d3^n6*d4^n2*d6^n3*d7^n4*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d10^n9*d3^n7*d4^n3*d6^n2*d7^n4*d8^n5*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n8))/(d4^n2*d5^n3*d7^n4*d8^n5*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n9))/(d4^n3*d5^n2*d7^n4*d8^n5*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n4*d4^n1*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n4*d4^n1*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n4*d4^n1*d5^n5*d6^n7*d7^n6*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n8))/(d1^n4*d4^n1*d5^n5*d7^n7*d8^n6*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n9))/(d1^n4*d4^n1*d5^n5*d7^n6*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d10^n10*d2^n2*d4^n8*d6^n1*d7^n4*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d10^n10*d2^n3*d4^n9*d6^n1*d7^n4*d8^n6*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n2*d4^n6*d6^n1*d7^n4*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n3*d4^n7*d6^n1*d7^n4*d8^n6*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n6)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3))/(d1^n2*d4^n8*d5^n1*d7^n4*d8^n9*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2))/(d3^n3*d4^n7*d5^n1*d7^n4*d8^n8*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3))/(d3^n2*d4^n6*d5^n1*d7^n4*d8^n9*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk1.tk1 - 4*tk1.tk2 + 2*tk1.tk4 + 4*tk2.tk2 - 4*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3))/(d3^n2*d4^n6*d5^n1*d7^n4*d8^n7*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d3^n4*d4^n7*d5^n1*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n6*d5^n1*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 2*tk1.tk4 + 4*tk2.tk2 - 4*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4))/(d4^n6*d5^n1*d7^n3*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n2*d3^n1*d4^n3*d6^n6*d7^n9*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d1^n4*d4^n5*d5^n1*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d2^n2*d3^n10*d4^n8*d5^n4*d6^n1*d7^n3*d8^n5*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n7)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n4*d4^n5*d5^n1*d7^n3*d8^n6*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d1^n4*d2^n7*d4^n1*d5^n5*d7^n6*d8^n3*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n2*d4^n6*d5^n4*d6^n1*d7^n3*d8^n5*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d4^n4*d5^n1*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n1*d5^n3*d7^n4*d8^n8*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n4*d5^n1*d7^n3*d8^n6*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d4^n3*d6^n1*d7^n4*d8^n6*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n8)*
num(-1 + 9*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 6*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 6*tk2.tp + tk4.tk4 + 6*tk4.tp, -n10))/(d4^n1*d5^n3*d7^n4*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n2))/(d4^n3*d5^n1*d7^n4*d8^n8*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2))/(d4^n3*d5^n1*d7^n4*d8^n6*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n7*d5^n1*d7^n2*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d3^n4*d4^n6*d5^n1*d7^n2*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5))/(d1^n3*d2^n4*d3^n1*d4^n2*d6^n6*d7^n8*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + tk1.tk1 - 4*tk1.tk2 + 2*tk1.tk4 + 4*tk2.tk2 - 4*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n4*d4^n6*d5^n1*d7^n2*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n6)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d1^n4*d4^n5*d5^n1*d7^n2*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n8))/(d1^n4*d2^n6*d4^n1*d5^n5*d7^n7*d8^n2*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n6)*
num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d1^n4*d4^n5*d5^n1*d7^n2*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d2^n3*d3^n10*d4^n9*d5^n4*d6^n1*d7^n2*d8^n5*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d1^n3*d4^n7*d5^n4*d6^n1*d7^n2*d8^n5*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n7))/(d4^n4*d5^n1*d7^n2*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3))/(d4^n2*d6^n1*d7^n4*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n4*d5^n1*d7^n2*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n1*d5^n2*d7^n4*d8^n9*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n9)*
num(-1 + 9*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 6*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 6*tk2.tp + tk4.tk4 + 6*tk4.tp, -n10))/(d4^n1*d5^n2*d7^n4*d8^n7*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n3))/(d4^n2*d5^n1*d7^n4*d8^n9*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3))/(d4^n2*d5^n1*d7^n4*d8^n7*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d4^n2*d5^n1*d7^n3*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d4^n3*d5^n1*d7^n2*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d4^n2*d5^n1*d7^n3*d8^n6*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n7))/(d4^n3*d5^n1*d7^n2*d8^n6*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + 9*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 6*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 6*tk2.tp + tk4.tk4 - 6*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8))/(d4^n1*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n6))/(d4^n2*d5^n1*d7^n3*d8^n5*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 + 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n7))/(d4^n3*d5^n1*d7^n2*d8^n5*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n2*d5^n4*d6^n1*d7^n3*d8^n5*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n7)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n3*d5^n4*d6^n1*d7^n2*d8^n5*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d4^n2*d5^n1*d6^n5*d7^n3*d8^n4*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d4^n3*d5^n1*d6^n5*d7^n2*d8^n4*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n1*d5^n3*d7^n4*d8^n2*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 - 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n1*d5^n2*d7^n4*d8^n3*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d4^n2*d5^n1*d6^n9*d7^n4*d8^n3*d9^n5);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n10*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n6*d3^n5*d4^n9*d6^n4*d7^n7*d8^n10*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n2)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d10^n6*d3^n5*d4^n9*d6^n4*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n3)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d10^n7*d3^n5*d4^n8*d6^n4*d7^n6*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n8*d2^n3*d3^n5*d4^n7*d5^n2*d6^n4*d7^n9*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n9*d2^n2*d3^n5*d4^n6*d5^n3*d6^n4*d7^n8*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n2)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5))/(d10^n4*d3^n3*d4^n10*d6^n6*d7^n8*d8^n7*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n4)*
num(-1 + 4*pp + tk2.tk2 - 4*tk2.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2))/(d3^n5*d4^n7*d5^n3*d6^n8*d7^n9*d8^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4))/(d2^n1*d4^n9*d5^n8*d6^n3*d7^n5*d8^n10*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d3^n7*d4^n9*d6^n3*d7^n5*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n2*d2^n4*d3^n9*d4^n7*d5^n8*d6^n3*d7^n5*d8^n10*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2))/(d4^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n1*d4^n9*d5^n3*d6^n4*d7^n8*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d1^n1*d3^n3*d4^n8*d5^n2*d6^n4*d7^n10*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n8*d3^n3*d4^n7*d6^n6*d7^n4*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n9*d2^n2*d3^n5*d4^n6*d5^n3*d6^n4*d7^n8*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n2)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7))/(d10^n6*d3^n8*d4^n10*d5^n1*d6^n4*d7^n3*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d2^n8*d4^n10*d5^n7*d6^n3*d7^n5*d8^n4*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2))/(d2^n1*d4^n5*d5^n3*d6^n4*d7^n8*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n8*d5^n7*d6^n3*d7^n5*d8^n4*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d1^n3*d2^n1*d3^n8*d4^n5*d6^n4*d7^n6*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n1*d4^n4*d5^n8*d6^n3*d7^n5*d8^n10*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n2*d2^n10*d3^n6*d4^n8*d5^n7*d6^n3*d7^n5*d8^n4*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n2*d2^n9*d3^n4*d4^n7*d5^n8*d6^n5*d7^n3*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d3^n8*d4^n6*d5^n7*d6^n3*d7^n5*d8^n4*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n3)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d3^n2*d4^n10*d6^n6*d7^n9*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n4))/(d2^n1*d4^n8*d5^n9*d6^n2*d7^n5*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n4)*
num(-1 + 4*pp + tk2.tk2 - 4*tk2.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3))/(d3^n5*d4^n6*d5^n2*d6^n9*d7^n8*d8^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n3*d2^n4*d3^n8*d4^n6*d5^n9*d6^n2*d7^n5*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d3^n6*d4^n8*d6^n2*d7^n5*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3))/(d4^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n1*d3^n2*d4^n9*d5^n3*d6^n4*d7^n10*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n1*d4^n8*d5^n2*d6^n4*d7^n9*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n8*d2^n3*d3^n5*d4^n7*d5^n2*d6^n4*d7^n9*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n9*d3^n2*d4^n6*d6^n7*d7^n4*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5)/(d1^n10*d10^n1*d2^n3*d3^n2*d4^n4*d6^n7*d7^n6*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n3)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6))/(d10^n7*d3^n9*d4^n10*d5^n1*d6^n4*d7^n2*d8^n5*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n1*d4^n4*d5^n9*d6^n2*d7^n5*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n2*d2^n1*d3^n9*d4^n5*d6^n4*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n3*d2^n10*d3^n7*d4^n9*d5^n6*d6^n2*d7^n5*d8^n4*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3))/(d2^n1*d4^n5*d5^n2*d6^n4*d7^n9*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d2^n9*d4^n10*d5^n6*d6^n2*d7^n5*d8^n4*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n9*d5^n6*d6^n2*d7^n5*d8^n4*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n3*d2^n8*d3^n4*d4^n6*d5^n9*d6^n5*d7^n2*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n3*d3^n2*d4^n4*d5^n5*d6^n7*d7^n6*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d3^n9*d4^n7*d5^n6*d6^n2*d7^n5*d8^n4*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d4^n7*d5^n10*d6^n2*d7^n3*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d4^n6*d5^n10*d6^n2*d7^n3*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n5*d2^n8*d3^n9*d4^n10*d6^n3*d7^n2*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n5*d2^n9*d3^n8*d4^n10*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d1^n2*d3^n10*d4^n3*d6^n6*d7^n8*d8^n7*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n6)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d4^n5*d5^n10*d6^n2*d7^n3*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n1*d4^n3*d5^n9*d6^n2*d7^n5*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n4*d2^n1*d3^n3*d4^n8*d6^n2*d7^n10*d8^n5*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n6*d3^n5*d4^n7*d6^n2*d7^n9*d8^n3*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d1^n4*d2^n1*d3^n2*d4^n9*d6^n3*d7^n10*d8^n5*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n4)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d1^n1*d4^n2*d5^n8*d6^n3*d7^n5*d8^n10*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n7*d2^n4*d3^n5*d4^n6*d6^n2*d7^n8*d8^n3*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n3*d2^n4*d3^n8*d4^n6*d6^n2*d7^n5*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d1^n8*d2^n10*d3^n5*d4^n9*d6^n2*d7^n7*d8^n3*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n9*d2^n1*d3^n5*d4^n8*d6^n2*d7^n6*d8^n3*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n5*d2^n6*d3^n7*d4^n4*d5^n10*d6^n2*d7^n3*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n6)*
num(-1 + 4*pp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d10^n8*d3^n3*d4^n4*d5^n2*d6^n9*d7^n7*d8^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d3^n1*d4^n2*d5^n3*d6^n4*d7^n8*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d3^n4*d4^n7*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d3^n1*d4^n3*d5^n2*d6^n4*d7^n9*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n7)*
num(-1 + 4*pp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d10^n9*d3^n2*d4^n4*d5^n3*d6^n8*d7^n6*d8^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d3^n4*d4^n6*d5^n10*d6^n3*d7^n2*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + 9*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 6*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 6*tk2.tp + tk4.tk4 - 6*tk4.tp, -n1)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n5))/(d4^n10*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n1*d4^n9*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n9))/(d2^n1*d4^n8*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)/(d10^n6*d2^n8*d3^n2*d4^n10*d5^n7*d6^n3*d7^n9*d8^n4*d9^n5);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)/(d10^n7*d2^n9*d3^n3*d4^n10*d5^n6*d6^n2*d7^n8*d8^n4*d9^n5);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7))/(d1^n2*d2^n10*d4^n8*d5^n1*d6^n3*d7^n4*d8^n9*d9^n5);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d1^n8*d2^n6*d3^n3*d4^n4*d6^n2*d7^n7*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d10^n8*d4^n2*d5^n7*d6^n3*d7^n5*d8^n4*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n10*d3^n3*d4^n7*d6^n2*d7^n4*d8^n5*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n9*d2^n7*d3^n2*d4^n4*d6^n3*d7^n6*d8^n5*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n10*d3^n2*d4^n6*d6^n3*d7^n4*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d10^n9*d4^n3*d5^n6*d6^n2*d7^n5*d8^n4*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n1*d4^n5*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n4*d4^n1*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)/(d1^n6*d10^n10*d2^n2*d4^n8*d5^n5*d6^n1*d7^n7*d8^n4*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)/(d1^n7*d10^n10*d2^n3*d4^n9*d5^n5*d6^n1*d7^n6*d8^n4*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n6*d5^n5*d6^n1*d7^n7*d8^n4*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n7*d5^n5*d6^n1*d7^n6*d8^n4*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n2*d3^n1*d4^n3*d6^n6*d7^n9*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d2^n2*d3^n6*d4^n8*d5^n9*d6^n1*d7^n5*d8^n3*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk2.tk2 - 4*tk2.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n4)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + 4*pp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n10))/(d3^n5*d4^n6*d5^n1*d6^n3*d7^n8*d8^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n8)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n6*d5^n9*d6^n1*d7^n5*d8^n3*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d2^n2*d3^n6*d4^n8*d5^n9*d6^n1*d7^n5*d8^n3*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n10*d2^n2*d4^n8*d5^n3*d6^n1*d7^n4*d8^n9*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d1^n10*d3^n3*d4^n9*d5^n2*d6^n4*d7^n1*d8^n6*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n2*d4^n6*d5^n3*d6^n1*d7^n4*d8^n9*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n3*d4^n7*d6^n1*d7^n4*d8^n6*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7))/(d1^n2*d2^n10*d4^n8*d5^n1*d6^n3*d7^n4*d8^n9*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d2^n2*d3^n10*d4^n8*d5^n4*d6^n1*d7^n3*d8^n5*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7))/(d10^n10*d3^n2*d4^n6*d5^n1*d6^n3*d7^n4*d8^n9*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n2))/(d4^n3*d5^n5*d6^n1*d7^n6*d8^n4*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n2*d4^n6*d5^n4*d6^n1*d7^n3*d8^n5*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5))/(d1^n3*d2^n4*d3^n1*d4^n2*d6^n6*d7^n8*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n4)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d3^n1*d4^n2*d5^n5*d6^n7*d7^n8*d8^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d2^n3*d3^n7*d4^n9*d5^n8*d6^n1*d7^n5*d8^n2*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d2^n3*d3^n7*d4^n9*d5^n8*d6^n1*d7^n5*d8^n2*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 + 4*tk4.tp, -n10))/(d4^n7*d5^n8*d6^n1*d7^n5*d8^n2*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n10*d3^n2*d4^n8*d5^n3*d6^n4*d7^n1*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d1^n10*d2^n3*d4^n9*d5^n2*d6^n1*d7^n4*d8^n8*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 - 4*tk1.tk2 + 4*tk2.tk2, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n2*d4^n6*d6^n1*d7^n4*d8^n7*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n3*d4^n7*d5^n2*d6^n1*d7^n4*d8^n8*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n7)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6))/(d1^n3*d2^n10*d4^n9*d5^n1*d6^n2*d7^n4*d8^n8*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n3))/(d4^n2*d5^n5*d6^n1*d7^n7*d8^n4*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6))/(d10^n10*d3^n3*d4^n7*d5^n1*d6^n2*d7^n4*d8^n8*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d2^n3*d3^n10*d4^n9*d5^n4*d6^n1*d7^n2*d8^n5*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d1^n3*d4^n7*d5^n4*d6^n1*d7^n2*d8^n5*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d3^n4*d4^n7*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d2^n4*d4^n6*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n1*d4^n9*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n1*d4^n8*d6^n3*d7^n2*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n4*d4^n5*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d4^n2*d5^n9*d6^n1*d7^n5*d8^n3*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d2^n10*d3^n2*d4^n8*d6^n3*d7^n1*d8^n5*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d2^n10*d3^n3*d4^n9*d6^n2*d7^n1*d8^n5*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n10)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n4)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d4^n3*d5^n8*d6^n1*d7^n5*d8^n2*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n1*d4^n5*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n6)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n4*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n2*d5^n3*d6^n1*d7^n4*d8^n9*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n6))/(d10^n10*d3^n3*d4^n4*d5^n1*d6^n2*d7^n7*d8^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n7))/(d10^n10*d3^n2*d4^n4*d5^n1*d6^n3*d7^n6*d8^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n3*d5^n2*d6^n1*d7^n4*d8^n8*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n5)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n9)*
num(-1 + 9*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 6*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 6*tk2.tp + tk4.tk4 - 6*tk4.tp, -n10)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8))/(d4^n1*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n2*d5^n1*d6^n3*d7^n4*d8^n9*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 - 4*tk1.tp + tk2.tk2 - 4*tk2.tp, -n7)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d4^n3*d5^n1*d6^n2*d7^n4*d8^n8*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n2*d5^n4*d6^n1*d7^n3*d8^n5*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + 4*pp + tk1.tk1 - 4*tk1.tp, -n7)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d4^n3*d5^n4*d6^n1*d7^n2*d8^n5*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n7*d3^n5*d4^n8*d6^n4*d7^n6*d8^n10*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d3^n3*d4^n10*d5^n5*d6^n7*d7^n8*d8^n9*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n8*d3^n3*d4^n7*d6^n6*d7^n4*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n1*d3^n5*d4^n9*d5^n3*d6^n8*d7^n7*d8^n10*d9^n4);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n2)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d1^n3*d2^n1*d3^n8*d4^n5*d6^n4*d7^n6*d8^n9*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n9*d2^n2*d3^n7*d4^n4*d6^n6*d7^n3*d8^n10*d9^n5);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n9*d2^n2*d3^n5*d4^n6*d5^n3*d6^n4*d7^n8*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d3^n2*d4^n10*d5^n5*d6^n6*d7^n9*d8^n8*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n9*d3^n2*d4^n6*d6^n7*d7^n4*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n3)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n2*d2^n1*d3^n9*d4^n5*d6^n4*d7^n7*d8^n8*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d1^n1*d3^n5*d4^n8*d5^n2*d6^n9*d7^n6*d8^n10*d9^n4);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n8*d2^n3*d3^n5*d4^n7*d5^n2*d6^n4*d7^n9*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d1^n8*d2^n3*d3^n6*d4^n4*d6^n7*d7^n2*d8^n10*d9^n5);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n3*d3^n2*d4^n4*d5^n5*d6^n7*d7^n6*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n1)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4))/(d1^n2*d3^n10*d4^n3*d6^n6*d7^n8*d8^n7*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n1)/(d10^n4*d2^n6*d3^n3*d4^n8*d5^n9*d6^n5*d7^n10*d8^n7*d9^n2);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n1)/(d10^n4*d2^n7*d3^n2*d4^n9*d5^n8*d6^n5*d7^n10*d8^n6*d9^n3);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n2*d2^n4*d3^n9*d4^n7*d6^n3*d7^n5*d8^n6*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n3*d2^n4*d3^n8*d4^n6*d6^n2*d7^n5*d8^n7*d9^n10);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n6*d3^n5*d4^n7*d6^n2*d7^n9*d8^n3*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n1*d3^n2*d4^n9*d5^n3*d6^n4*d7^n10*d8^n7*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n5)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d1^n1*d3^n3*d4^n8*d5^n2*d6^n4*d7^n10*d8^n6*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n8*d2^n5*d3^n3*d4^n7*d5^n2*d6^n9*d7^n4*d8^n10*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1)/(d10^n9*d2^n5*d3^n2*d4^n6*d5^n3*d6^n8*d7^n4*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n1)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n2*d2^n10*d3^n4*d4^n3*d6^n8*d7^n7*d8^n9*d9^n6);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d1^n2*d3^n10*d4^n8*d5^n1*d6^n4*d7^n3*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)/(d10^n6*d2^n8*d3^n2*d4^n10*d5^n7*d6^n3*d7^n9*d8^n4*d9^n5);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n10*d3^n3*d4^n7*d6^n2*d7^n4*d8^n5*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)/(d1^n5*d2^n6*d3^n7*d4^n4*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?pos_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n10*d3^n2*d4^n6*d6^n3*d7^n4*d8^n5*d9^n9);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?pos_ = intFG*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n1)/(d10^n7*d2^n9*d3^n3*d4^n10*d5^n6*d6^n2*d7^n8*d8^n4*d9^n5);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)/(d1^n5*d2^n7*d3^n6*d4^n4*d5^n1*d6^n3*d7^n2*d8^n9*d9^n8);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?pos_ = intFG*
(num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tk4 + 2*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n1))/(d1^n3*d2^n8*d3^n4*d4^n6*d6^n5*d7^n2*d8^n10*d9^n7);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n1)/(d10^n8*d2^n10*d3^n5*d4^n9*d5^n4*d6^n2*d7^n7*d8^n6*d9^n3);

id int`TOPO'/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n1)/(d10^n9*d2^n10*d3^n5*d4^n8*d5^n4*d6^n3*d7^n6*d8^n7*d9^n2);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d10^n2*d3^n1*d4^n3*d5^n5*d6^n6*d7^n9*d8^n8*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n2)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7)/(d1^n2*d2^n10*d3^n6*d4^n8*d5^n1*d6^n3*d7^n5*d8^n9*d9^n4);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n2)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n8*d7^n7*d8^n6*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n2)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n6*d4^n8*d5^n4*d6^n3*d7^n5*d8^n7*d9^n1);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n4)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n10))/(d10^n3*d3^n1*d4^n2*d5^n5*d6^n7*d7^n8*d8^n9*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n3*d3^n2*d4^n4*d5^n5*d6^n7*d7^n6*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6)/(d1^n3*d2^n10*d3^n7*d4^n9*d5^n1*d6^n2*d7^n5*d8^n8*d9^n4);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk1.tk1 + 2*tk1.tk4 + tk4.tk4, -n3)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk4.tk4 + 2*tk4.tp, -n7))/(d1^n1*d3^n2*d4^n4*d5^n5*d7^n6*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n3)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk4 - 4*tk1.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n7*d4^n9*d5^n4*d6^n2*d7^n5*d8^n6*d9^n1);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n3)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n1*d3^n2*d4^n4*d5^n5*d6^n9*d7^n6*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n4)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n5))/(d1^n2*d3^n1*d4^n3*d6^n6*d7^n9*d8^n7*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n4)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10))/(d10^n6*d3^n5*d4^n7*d5^n1*d6^n2*d7^n9*d8^n8*d9^n3);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n4)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10))/(d10^n7*d3^n5*d4^n6*d5^n1*d6^n3*d7^n8*d8^n9*d9^n2);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n7*d3^n9*d4^n3*d5^n4*d6^n2*d7^n1*d8^n6*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?neg0_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n6*d3^n8*d4^n2*d5^n4*d6^n3*d7^n1*d8^n7*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n6))/(d3^n4*d4^n7*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n5)*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n7))/(d3^n4*d4^n6*d5^n1*d6^n3*d7^n2*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d3^n2*d4^n6*d5^n3*d6^n1*d7^n4*d8^n9*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?neg0_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n5)*
num(-1 + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n9))/(d3^n3*d4^n7*d5^n2*d6^n1*d7^n4*d8^n8*d9^n6);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n6)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7))/(d1^n2*d2^n10*d4^n8*d5^n1*d6^n3*d7^n4*d8^n9*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n6))/(d1^n1*d2^n8*d3^n3*d4^n4*d5^n5*d7^n7*d8^n2*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n9)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n6))/(d10^n10*d3^n3*d4^n7*d5^n1*d6^n2*d7^n4*d8^n8*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 - 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 + 2*tk2.tp, -n8)*
num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n7))/(d10^n10*d3^n2*d4^n6*d5^n1*d6^n3*d7^n4*d8^n9*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n10)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 + 2*tk4.tp, -n7))/(d1^n1*d2^n9*d3^n2*d4^n4*d5^n5*d7^n6*d8^n3*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + pp + tk1.tk1 + 2*tk1.tk2 - 2*tk1.tp + tk2.tk2 - 2*tk2.tp, -n9)*
num(-1 + 4*pp + tk1.tk1 + 2*tk1.tk2 + 2*tk1.tk4 - 4*tk1.tp + tk2.tk2 + 2*tk2.tk4 - 4*tk2.tp + tk4.tk4 - 4*tk4.tp, -n10)*
num(-1 + pp + tk1.tk1 + 2*tk1.tk4 - 2*tk1.tp + tk4.tk4 - 2*tk4.tp, -n8))/(d1^n1*d4^n5*d5^n4*d6^n2*d7^n3*d8^n6*d9^n7);

id int`TOPO'/d1^n1?pos_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?neg0_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
1/(d1^n1*d10^n10*d2^n3*d3^n2*d4^n4*d5^n5*d6^n7*d7^n6*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)/(d1^n5*d2^n6*d3^n7*d4^n4*d5^n1*d6^n2*d7^n3*d8^n8*d9^n9);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ = intFG*
num(-1 + pp + tk2.tk2 + 2*tk2.tk4 - 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10)/(d1^n5*d2^n7*d3^n6*d4^n4*d5^n1*d6^n3*d7^n2*d8^n9*d9^n8);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?neg0_/d9^n9?pos_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n8)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n7*d3^n9*d4^n3*d5^n4*d6^n2*d7^n1*d8^n6*d9^n5);

id int`TOPO'/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?neg0_/d10^n10?neg0_ = intFG*
(num(-1 + tk2.tk2 - 2*tk2.tk4 + tk4.tk4, -n9)*
num(-1 + pp + tk1.tk1 - 2*tk1.tk2 + 2*tk1.tk4 - 2*tk1.tp + tk2.tk2 - 2*tk2.tk4 + 2*tk2.tp + tk4.tk4 - 2*tk4.tp, -n10))/(d1^n6*d3^n8*d4^n2*d5^n4*d6^n3*d7^n1*d8^n7*d9^n5);
#endprocedure
*--#] mapFG :



* Zero rules
#procedure zerofacet3(N1,N2,N3)
        if( (count(d`N1',-1) <= 0 ) && (count(d`N2',-1) <= 0 ) && (count(d`N3',-1) <= 0 )) Discard;
#endprocedure

#procedure zerofacet4(N1,N2,N3,N4)
        if( (count(d`N1',-1) <= 0 ) && (count(d`N2',-1) <= 0 ) && (count(d`N3',-1) <= 0 ) && (count(d`N4',-1) <= 0 )) Discard;
#endprocedure

#procedure zerofacet5(N1,N2,N3,N4,N5)
        if( (count(d`N1',-1) <= 0 ) && (count(d`N2',-1) <= 0 ) && (count(d`N3',-1) <= 0 ) && (count(d`N4',-1) <= 0 ) && (count(d`N5',-1) <= 0 )) Discard;
#endprocedure

#procedure zerofacet6(N1,N2,N3,N4,N5,N6)
        if( (count(d`N1',-1) <= 0 ) && (count(d`N2',-1) <= 0 ) && (count(d`N3',-1) <= 0 ) && (count(d`N4',-1) <= 0 ) && (count(d`N5',-1) <= 0 ) && (count(d`N6',-1) <= 0 )) Discard;
#endprocedure

*--#[ zeroX :
#procedure zeroX
        if((count(intX,1) > 0) && (count(d1,1) >= 0));
        #call zerofacet4(2,3,4,10)
        #call zerofacet4(2,4,7,9)
        #call zerofacet4(2,5,6,9)
        #call zerofacet4(2,6,8,10)
        #call zerofacet4(3,4,6,8)
        #call zerofacet4(3,5,7,8)
        #call zerofacet4(3,7,9,10)
        #call zerofacet4(4,5,6,7)
        #call zerofacet4(5,8,9,10)
        #call zerofacet6(2,3,4,5,8,9)
        #call zerofacet6(2,3,5,6,7,10)
        #call zerofacet6(2,3,6,7,8,9)
        #call zerofacet6(2,4,5,7,8,10)
        #call zerofacet6(3,4,5,6,9,10)
        #call zerofacet6(4,6,7,8,9,10)
        endif;
#endprocedure
*--#] zeroX :

*--#[ zeroH :
#procedure zeroH
        if((count(intH,1) > 0) && (count(d10,1) >= 0));
        #call zerofacet3(2,6,8)
        #call zerofacet3(3,7,9)
        #call zerofacet4(1,2,3,4)
        #call zerofacet4(1,5,8,9)
        #call zerofacet4(4,5,6,7)
        #call zerofacet5(1,2,4,7,9)
        #call zerofacet5(1,2,5,6,9)
        #call zerofacet5(1,3,4,6,8)
        #call zerofacet5(1,3,5,7,8)
        #call zerofacet5(2,4,5,7,8)
        #call zerofacet5(3,4,5,6,9)
        #call zerofacet6(1,2,3,5,6,7)
        #call zerofacet6(1,4,6,7,8,9)
        #call zerofacet6(2,3,4,5,8,9)
        endif;
#endprocedure
*--#] zeroH :

*--#[ zeroBMW :
#procedure zeroBMW
        if((count(intBMW,1) > 0) && (count(d1,1) >= 0) && (count(d2,1) >= 0));
        #call zerofacet3(3,4,10)
        #call zerofacet3(4,7,9)
        #call zerofacet3(5,6,9)
        #call zerofacet3(6,8,10)
        #call zerofacet4(3,4,6,8)
        #call zerofacet4(3,5,7,8)
        #call zerofacet4(3,7,9,10)
        #call zerofacet4(4,5,6,7)
        #call zerofacet4(5,8,9,10)
        #call zerofacet5(3,4,5,8,9)
        #call zerofacet5(3,5,6,7,10)
        #call zerofacet5(3,6,7,8,9)
        #call zerofacet5(4,5,7,8,10)
        endif;
#endprocedure
*--#] zeroBMW :

* 
* Non-planar 9 line
*         
*--#[ redX :
#procedure redX
#$repcount = 1;        
#do irep=1,1
        #$irep = 1;                
        if(count(intX,1));

* den: minus
* n6 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?{>1}/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        (d2^(-1 - n2)*d4^(1 - n4))/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d6^(1 - n6)*d9^(1 - n9)*n1)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d6^(1 - n6)*d8^(1 - n8)*n1)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d5^(1 - n5)*d6^(1 - n6)*n1)/(d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*d6^(1 - n6)*n1)/(d10^n10*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d3^(1 - n3)*d6^(1 - n6)*n1)/(d10^n10*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d6^(1 - n6)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d10^(-1 - n10)*d2^(-1 - n2)*d6^(1 - n6)*d8^(1 - n8)*n10)/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (d10^(-1 - n10)*d2^(-1 - n2)*d3^(1 - n3)*d6^(1 - n6)*n10)/(d1^n1*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-2 - n2)*d4^(1 - n4)*d6^(1 - n6)*(-1 - n2))/(d1^n1*d10^n10*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-2 - n2)*d6^(2 - n6)*(1 + n2))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d2^(-1 - n2)*d3^(-1 - n3)*d6^(1 - n6)*d8^(1 - n8)*n3)/(d1^n1*d10^n10*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(-1 - n3)*d6^(1 - n6)*n3)/(d1^n1*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d4^(-1 - n4)*d6^(1 - n6)*n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d2^(-1 - n2)*d4^(-1 - n4)*d6^(2 - n6)*n4)/(d1^n1*d10^n10*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d2^(-1 - n2)*d5^(-1 - n5)*d6^(1 - n6)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n6)) + 
        (d2^(-1 - n2)*d5^(-1 - n5)*d6^(1 - n6)*d7^(1 - n7)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-1 - n2)*d6^(1 - n6)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d8^n8*(-1 + n6)) - 
        (d2^(-1 - n2)*d5^(1 - n5)*d6^(1 - n6)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d1^n1*d10^n10*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (d2^(-1 - n2)*d6^(1 - n6)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d8^n8*(-1 + n6)) + 
        (d2^(-1 - n2)*d5^(1 - n5)*d6^(1 - n6)*d9^(-1 - n9)*n9)/(d1^n1*d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n6));

* n10 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?{>1} =
        
        (d2^(-1 - n2)*d9^(1 - n9))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) - 
        (d2^(-1 - n2)*d8^(1 - n8))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d10^(1 - n10)*d2^(-2 - n2)*d9^(1 - n9)*(-1 - n2))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-2 - n2)*d6^(1 - n6)*(-1 - n2))/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-2 - n2)*d4^(1 - n4)*(1 + n2))/(d1^n1*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(2 - n10)*d2^(-2 - n2)*(1 + n2))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(-1 - n3)*d8^(1 - n8)*n3)/(d1^n1*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(-1 - n3)*d7^(1 - n7)*n3)/(d1^n1*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d1^n1*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(2 - n10)*d2^(-1 - n2)*d3^(-1 - n3)*n3)/(d1^n1*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d4^(-1 - n4)*d7^(1 - n7)*n4)/(d1^n1*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(d1^n1*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(1 - n3)*d4^(-1 - n4)*n4)/(d1^n1*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d4^(-1 - n4)*n4)/(d1^n1*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d5^(-1 - n5)*d8^(1 - n8)*n5)/(d1^n1*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d5^(-1 - n5)*d7^(1 - n7)*n5)/(d1^n1*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d5^(-1 - n5)*d6^(1 - n6)*n5)/(d1^n1*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(d1^n1*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d4^(1 - n4)*d6^(-1 - n6)*n6)/(d1^n1*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d6^(-1 - n6)*n6)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d1^n1*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d4^(1 - n4)*d7^(-1 - n7)*n7)/(d1^n1*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(1 - n3)*d7^(-1 - n7)*n7)/(d1^n1*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(d1^n1*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(1 - n3)*d8^(-1 - n8)*n8)/(d1^n1*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) + 
        (d10^(2 - n10)*d2^(-1 - n2)*d8^(-1 - n8)*n8)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d5^(1 - n5)*d9^(-1 - n9)*n9)/(d1^n1*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n10)) - 
        (d10^(2 - n10)*d2^(-1 - n2)*d9^(-1 - n9)*n9)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d10^(1 - n10)*d9^(-1 - n9)*n9)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10));

* * n6 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?{>1}/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        -d8^(1 - n8)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) + 
        d5^(1 - n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (2*d4^(1 - n4))/(3*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (2*d2^(1 - n2))/(3*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d6^(1 - n6)*d9^(1 - n9)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) + 
        (4*d1^(-1 - n1)*d6^(1 - n6)*d8^(1 - n8)*n1)/(3*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (2*d1^(-1 - n1)*d5^(1 - n5)*d6^(1 - n6)*n1)/(3*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (2*d1^(-1 - n1)*d4^(1 - n4)*d6^(1 - n6)*n1)/(3*d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d6^(1 - n6)*n1)/(d10^n10*d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (4*d1^(-1 - n1)*d2^(1 - n2)*d6^(1 - n6)*n1)/(3*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d6^(1 - n6)*n1)/(3*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (2*d10^(-1 - n10)*d6^(1 - n6)*d8^(1 - n8)*n10)/(3*d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) + 
        (2*d10^(-1 - n10)*d3^(1 - n3)*d6^(1 - n6)*n10)/(3*d1^n1*d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-1 - n2)*d6^(1 - n6)*d9^(1 - n9)*n2)/(3*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) + 
        (2*d2^(-1 - n2)*d4^(1 - n4)*d6^(1 - n6)*n2)/(3*d1^n1*d10^n10*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d6^(1 - n6)*n2)/(3*d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d2^(-1 - n2)*d6^(1 - n6)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (2*d2^(-1 - n2)*d6^(2 - n6)*n2)/(3*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (2*d3^(-1 - n3)*d6^(1 - n6)*d8^(1 - n8)*n3)/(3*d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (2*d10^(1 - n10)*d3^(-1 - n3)*d6^(1 - n6)*n3)/(3*d1^n1*d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d4^(-1 - n4)*d6^(1 - n6)*d7^(1 - n7)*n4)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d8^n8*d9^n9*(-1 + n6)) - 
        (d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d1^n1*d10^n10*d2^n2*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (2*d2^(1 - n2)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d1^n1*d10^n10*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d4^(-1 - n4)*d6^(1 - n6)*n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (2*d4^(-1 - n4)*d6^(2 - n6)*n4)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (2*d5^(-1 - n5)*d6^(1 - n6)*d9^(1 - n9)*n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n6)) + 
        (2*d5^(-1 - n5)*d6^(1 - n6)*d8^(1 - n8)*n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d9^n9*(-1 + n6)) - 
        (2*d5^(-1 - n5)*d6^(1 - n6)*d7^(1 - n7)*n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n6)) - 
        (2*d5^(-1 - n5)*d6^(2 - n6)*n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (2*d6^(1 - n6)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d8^n8*(-1 + n6)) + 
        (2*d5^(1 - n5)*d6^(1 - n6)*d7^(-1 - n7)*n7)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n6)) - 
        (2*d4^(1 - n4)*d6^(1 - n6)*d7^(-1 - n7)*n7)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d8^n8*d9^n9*(-1 + n6)) + 
        (2*d3^(1 - n3)*d6^(1 - n6)*d7^(-1 - n7)*n7)/(3*d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d8^n8*d9^n9*(-1 + n6)) - 
        (2*d3^(1 - n3)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(3*d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) + 
        (2*d10^(1 - n10)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(3*d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) + 
        (2*d6^(1 - n6)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d8^n8*(-1 + n6)) - 
        (2*d5^(1 - n5)*d6^(1 - n6)*d9^(-1 - n9)*n9)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n6)) - 
        (2*d2^(1 - n2)*d6^(1 - n6)*d9^(-1 - n9)*n9)/(3*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) + 
        (2*d10^(1 - n10)*d6^(1 - n6)*d9^(-1 - n9)*n9)/(3*d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) + 
        (d6^(1 - n6)*rat(3 + 2*d - n1 - 3*n2 - 3*n4 - 3*n6, 3*(-1 + n6)))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9);
        
* * n1 != 0 && n10 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?{>1} =
        
        -((d1^(-1 - n1)*d7^(1 - n7))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9)) - 
        (d1^(-1 - n1)*d6^(1 - n6))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d5^(1 - n5))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d4^(1 - n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d3^(1 - n3))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d2^(1 - n2))/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d10^(1 - n10)*d4^(1 - n4)*(-1 - n1))/(d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d1^(-2 - n1)*d10^(1 - n10)*(-1 - n1))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d1^(-2 - n1)*d10^(1 - n10)*d5^(1 - n5)*(1 + n1))/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d4^(1 - n4)*d5^(-1 - n5)*n5)/(d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d5^(-1 - n5)*n5)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d5^(-1 - n5)*n5)/(d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*(-1 - n1 + n5))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d2^(1 - n2)*d8^(-1 - n8)*n8)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d8^(-1 - n8)*n8)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d5^(1 - n5)*d9^(-1 - n9)*n9)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d3^(1 - n3)*d9^(-1 - n9)*n9)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10)) - 
        (d10^(1 - n10)*d9^(-1 - n9)*n9)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10));

* * n7 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?{>1}/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =

        -d9^(1 - n9)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) + 
        d5^(1 - n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (2*d4^(1 - n4))/(3*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (2*d3^(1 - n3))/(3*d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d7^(1 - n7)*d9^(1 - n9)*n1)/(3*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) - 
        (d1^(-1 - n1)*d7^(1 - n7)*d8^(1 - n8)*n1)/(3*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d7^(1 - n7)*n1)/(3*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d1^(-1 - n1)*d2^(1 - n2)*d7^(1 - n7)*n1)/(3*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d2^(-1 - n2)*d7^(1 - n7)*d9^(1 - n9)*n2)/(3*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d7^(1 - n7)*n2)/(3*d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d2^(-1 - n2)*d7^(1 - n7)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d3^(-1 - n3)*d7^(1 - n7)*d8^(1 - n8)*n3)/(3*d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) - 
        (d10^(1 - n10)*d3^(-1 - n3)*d7^(1 - n7)*n3)/(3*d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d3^(-1 - n3)*d7^(1 - n7)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d4^(-1 - n4)*d6^(1 - n6)*d7^(1 - n7)*n4)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d8^n8*d9^n9*(-1 + n7)) + 
        (d3^(1 - n3)*d4^(-1 - n4)*d7^(1 - n7)*n4)/(3*d1^n1*d10^n10*d2^n2*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d2^(1 - n2)*d4^(-1 - n4)*d7^(1 - n7)*n4)/(3*d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d4^(-1 - n4)*d7^(2 - n7)*n4)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (2*d5^(-1 - n5)*d7^(1 - n7)*d9^(1 - n9)*n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*(-1 + n7)) - 
        (2*d5^(-1 - n5)*d7^(1 - n7)*d8^(1 - n8)*n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d9^n9*(-1 + n7)) + 
        (2*d5^(-1 - n5)*d6^(1 - n6)*d7^(1 - n7)*n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n7)) - 
        (2*d5^(-1 - n5)*d7^(2 - n7)*n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d6^(-1 - n6)*d7^(1 - n7)*d8^(1 - n8)*n6)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d9^n9*(-1 + n7)) - 
        (d5^(1 - n5)*d6^(-1 - n6)*d7^(1 - n7)*n6)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n7)) - 
        (2*d4^(1 - n4)*d6^(-1 - n6)*d7^(1 - n7)*n6)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d8^n8*d9^n9*(-1 + n7)) + 
        (2*d2^(1 - n2)*d6^(-1 - n6)*d7^(1 - n7)*n6)/(3*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d8^n8*d9^n9*(-1 + n7)) + 
        (d6^(-1 - n6)*d7^(1 - n7)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d8^n8*d9^n9*(-1 + n7)) + 
        (d7^(1 - n7)*(1 + n2 - n3 + n6 - n7))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (2*d3^(1 - n3)*d7^(1 - n7)*d8^(-1 - n8)*n8)/(3*d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) + 
        (2*d10^(1 - n10)*d7^(1 - n7)*d8^(-1 - n8)*n8)/(3*d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) + 
        (2*d2^(1 - n2)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(3*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) - 
        (2*d10^(1 - n10)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(3*d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7));

* * n8 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?{>1}/d9^n9?pos_/d10^n10?pos_ =

        d3^(1 - n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        d10^(1 - n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d1^(-1 - n1)*d8^(1 - n8)*d9^(1 - n9)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n8)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d8^(1 - n8)*n1)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d1^(-1 - n1)*d2^(1 - n2)*d8^(1 - n8)*n1)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d1^(-1 - n1)*d8^(2 - n8)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d10^(1 - n10)*d3^(-1 - n3)*d8^(1 - n8)*n3)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d3^(-1 - n3)*d8^(1 - n8)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d3^(-1 - n3)*d8^(2 - n8)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d4^(-1 - n4)*d6^(1 - n6)*d8^(1 - n8)*n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d9^n9*(-1 + n8)) + 
        (d2^(1 - n2)*d4^(-1 - n4)*d8^(1 - n8)*n4)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d4^(-1 - n4)*d8^(1 - n8)*n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d5^(-1 - n5)*d8^(1 - n8)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*(-1 + n8)) + 
        (d5^(-1 - n5)*d7^(1 - n7)*d8^(1 - n8)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d9^n9*(-1 + n8)) - 
        (d5^(-1 - n5)*d8^(1 - n8)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d4^(1 - n4)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d9^n9*(-1 + n8)) - 
        (d2^(1 - n2)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n8)) - 
        (d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n8)) + 
        (d7^(-1 - n7)*d8^(1 - n8)*d9^(1 - n9)*n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*(-1 + n8)) - 
        (d5^(1 - n5)*d7^(-1 - n7)*d8^(1 - n8)*n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d9^n9*(-1 + n8)) + 
        (d7^(-1 - n7)*d8^(1 - n8)*n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n8)) + 
        (d8^(1 - n8)*(1 + n3 + n4 - n5 - n6 + n7 - n8))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8));

* * n1 != 0 && n6 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?{>1}/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =

        (d1^(-1 - n1)*d9^(1 - n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) + 
        (d1^(-1 - n1)*d8^(1 - n8))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) - 
        (d1^(-1 - n1)*d7^(1 - n7))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d4^(1 - n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d3^(1 - n3))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d1^(-1 - n1)*d10^(1 - n10))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d6^(1 - n6)*d8^(1 - n8)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d6^(1 - n6)*n10)/(d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d6^(1 - n6)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d6^(1 - n6)*d9^(1 - n9)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d6^(1 - n6)*d8^(1 - n8)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d2^(-1 - n2)*d6^(1 - n6)*n2)/(d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d2^(-1 - n2)*d6^(1 - n6)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d6^(1 - n6)*(-n10 + n8))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6));

* * n10 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?{>1} =
        
        -(d1^(-1 - n1)*d10^(1 - n10)*d9^(1 - n9)*n1)/(2*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d3^(1 - n3)*n1)/(2*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*n1)/(2*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*n2)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(2*d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d10^(1 - n10)*d5^(-1 - n5)*d7^(1 - n7)*n5)/(2*d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d5^(-1 - n5)*n5)/(2*d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d4^(1 - n4)*d6^(-1 - n6)*n6)/(2*d1^n1*d2^n2*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(1 - n2)*d6^(-1 - n6)*n6)/(2*d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d6^(-1 - n6)*n6)/(2*d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d3^(1 - n3)*d8^(-1 - n8)*n8)/(2*d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d8^(-1 - n8)*n8)/(2*d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d10^(2 - n10)*d8^(-1 - n8)*n8)/(2*d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d9^(-1 - n9)*n9)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d10^(1 - n10)*rat(2 + 2*d - n1 - 2*n10 - 2*n2 - n5 - n6 - n8 - 2*n9, 2*(-1 + n10)))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);
                
* * n9 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?{>1}/d10^n10?pos_ =
        
        -(d2^(1 - n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9)) + 
        d10^(1 - n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (2*d1^(-1 - n1)*d8^(1 - n8)*d9^(1 - n9)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) - 
        (2*d1^(-1 - n1)*d2^(1 - n2)*d9^(1 - n9)*n1)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (2*d2^(-1 - n2)*d9^(1 - n9)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d3^(-1 - n3)*d8^(1 - n8)*d9^(1 - n9)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) - 
        (d10^(1 - n10)*d3^(-1 - n3)*d9^(1 - n9)*n3)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d3^(-1 - n3)*d9^(1 - n9)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d4^(-1 - n4)*d6^(1 - n6)*d9^(1 - n9)*n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d8^n8*(-1 + n9)) - 
        (d2^(1 - n2)*d4^(-1 - n4)*d9^(1 - n9)*n4)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d4^(-1 - n4)*d9^(1 - n9)*n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d5^(-1 - n5)*d8^(1 - n8)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*(-1 + n9)) - 
        (d5^(-1 - n5)*d6^(1 - n6)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n9)) + 
        (d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d5^(1 - n5)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*(-1 + n9)) - 
        (d4^(1 - n4)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d8^n8*(-1 + n9)) + 
        (d3^(1 - n3)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n9)) - 
        (d7^(-1 - n7)*d9^(2 - n9)*n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n9)) + 
        (2*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) + 
        (d9^(1 - n9)*(1 - 2*n2 + n3 - n4 + n5 + 2*n8 - n9))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9));
        
* * n1 != 0 && n7 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?{>1}/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        (d1^(-1 - n1)*d9^(1 - n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) + 
        (d1^(-1 - n1)*d8^(1 - n8))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) - 
        (d1^(-1 - n1)*d6^(1 - n6))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d4^(1 - n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d2^(1 - n2))/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d1^(-1 - n1)*d10^(1 - n10))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d7^(1 - n7)*d9^(1 - n9)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d7^(1 - n7)*n10)/(d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d7^(1 - n7)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d1^(-1 - n1)*d3^(-1 - n3)*d7^(1 - n7)*d9^(1 - n9)*n3)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) + 
        (d1^(-1 - n1)*d3^(-1 - n3)*d7^(1 - n7)*d8^(1 - n8)*n3)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d3^(-1 - n3)*d7^(1 - n7)*n3)/(d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d3^(-1 - n3)*d7^(1 - n7)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d1^(-1 - n1)*d2^(1 - n2)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) + 
        (d1^(-1 - n1)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) + 
        (d1^(-1 - n1)*d7^(1 - n7)*(-n10 + n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7));
 
* * n1 != 0 && n2 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?{>1}/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =

        (d1^(-1 - n1)*d8^(1 - n8))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) + 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d8^(1 - n8)*n10)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n2)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*n10)/(d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n2)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*n10)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n2)) + 
        (d1^(-1 - n1)*d2^(1 - n2)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n2)) - 
        (d1^(-1 - n1)*d2^(1 - n2)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n2)) + 
        (d1^(-1 - n1)*d2^(1 - n2)*d6^(-1 - n6)*n6)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n2)) + 
        (2*d1^(-1 - n1)*d2^(1 - n2)*d8^(-1 - n8)*n8)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n2)) + 
        (d1^(-1 - n1)*d2^(1 - n2)*rat(-1 - d + n10 + n2 + n6 + 2*n8, -1 + n2))/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);
        
* * n1 != 0 && n3 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?pos_/d3^n3?{>1}/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =

        (d1^(-1 - n1)*d9^(1 - n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) - 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d3^(1 - n3)*d8^(1 - n8)*(-1 - n1))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n3)) + 
        (d1^(-2 - n1)*d3^(1 - n3)*(-1 - n1))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) + 
        (d1^(-2 - n1)*d2^(1 - n2)*d3^(1 - n3)*(1 + n1))/(d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d9^(1 - n9)*n10)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n3)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*n10)/(d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*n10)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(-1 - n5)*d8^(1 - n8)*n5)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n3)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(-1 - n5)*d6^(1 - n6)*n5)/(d10^n10*d2^n2*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n3)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(-1 - n5)*n5)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n3)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n3)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d4^(1 - n4)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n3)) - 
        (d1^(-1 - n1)*d3^(2 - n3)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n3)) - 
        (2*d1^(-1 - n1)*d3^(1 - n3)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n3)) + 
        (d1^(-1 - n1)*d2^(1 - n2)*d3^(1 - n3)*d9^(-1 - n9)*n9)/(d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n3)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d3^(1 - n3)*d9^(-1 - n9)*n9)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n3)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d9^(-1 - n9)*n9)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n3)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*rat(d - n1 - n10 - n3 - n5 - 2*n8 + n9, -1 + n3))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

* * n1 != 0 && n4 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?pos_/d3^n3?pos_/d4^n4?{>1}/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =

        (d1^(-1 - n1)*d5^(1 - n5))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) - 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d4^(1 - n4)*(-2 - 2*n1))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d4^(1 - n4)*d9^(1 - n9)*n10)/(d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n4)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d4^(1 - n4)*d8^(1 - n8)*n10)/(d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n4)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d4^(1 - n4)*n10)/(d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d4^(1 - n4)*n10)/(d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*d8^(1 - n8)*n2)/(d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n4)) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*n2)/(d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4)*n2)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d1^(-1 - n1)*d3^(-1 - n3)*d4^(1 - n4)*d9^(1 - n9)*n3)/(d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n4)) - 
        (d1^(-1 - n1)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) - 
        (d3^(-1 - n3)*d4^(1 - n4)*n3)/(d1^n1*d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*rat(-1 + d - 2*n1 - n2 - n3 - n4, -1 + n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

* * n1 != 0 && n5 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?{>1}/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =

        (d1^(-1 - n1)*d4^(1 - n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d5^(1 - n5)*(-2 - 2*n1))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d5^(1 - n5)*d9^(1 - n9)*n10)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n5)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d5^(1 - n5)*d8^(1 - n8)*n10)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d5^(1 - n5)*n10)/(d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d5^(1 - n5)*n10)/(d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d5^(1 - n5)*d8^(1 - n8)*n2)/(d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d5^(1 - n5)*n2)/(d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) - 
        (d2^(-1 - n2)*d5^(1 - n5)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d3^(-1 - n3)*d5^(1 - n5)*d9^(1 - n9)*n3)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n5)) - 
        (d1^(-1 - n1)*d3^(-1 - n3)*d5^(1 - n5)*n3)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) - 
        (d3^(-1 - n3)*d5^(1 - n5)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d5^(1 - n5)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(d10^n10*d2^n2*d3^n3*d7^n7*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d2^(1 - n2)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d5^(2 - n5)*d6^(-1 - n6)*n6)/(d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d5^(1 - n5)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*(-1 + n5)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d3^n3*d6^n6*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d5^(2 - n5)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d5^(1 - n5)*rat(-1 + d - 2*n1 - n2 - n3 - n5, -1 + n5))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9);

* * n1 != 0 && n8 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?{>1}/d9^n9?pos_/d10^n10?pos_ =
        
        (d1^(-1 - n1)*d2^(1 - n2))/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d8^(1 - n8)*d9^(1 - n9)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n8)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d8^(1 - n8)*n10)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d8^(1 - n8)*n10)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d8^(2 - n8)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d8^(1 - n8)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d2^(-1 - n2)*d8^(1 - n8)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d8^(2 - n8)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d1^(-1 - n1)*d5^(1 - n5)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d9^n9*(-1 + n8)) - 
        (d1^(-1 - n1)*d4^(1 - n4)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d9^n9*(-1 + n8)) + 
        (d1^(-1 - n1)*d2^(1 - n2)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n8)) - 
        (d1^(-1 - n1)*d6^(-1 - n6)*d8^(2 - n8)*n6)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n8)) + 
        (d1^(-1 - n1)*d8^(1 - n8)*(1 + n2 - n8))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8));

* n1 != 0 && n9 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?{>1}/d10^n10?pos_ =

        (d1^(-1 - n1)*d3^(1 - n3))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d9^(1 - n9)*(-2 - 2*n1))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d8^(1 - n8)*d9^(1 - n9)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d9^(1 - n9)*n10)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d9^(1 - n9)*n10)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d9^(2 - n9)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d10^n10*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d1^(-1 - n1)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d1^(-1 - n1)*d2^(1 - n2)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) - 
        (d1^(-1 - n1)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) - 
        (d8^(-1 - n8)*d9^(1 - n9)*n8)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) + 
        (d1^(-1 - n1)*d9^(1 - n9)*rat(-1 + d - 2*n1 - n5 - n8 - n9, -1 + n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8);

* n1 == 0 && n10 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n5 != 1 && n5 != 2
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?{>2}/d6/d7/d8/d9/d10 =

        -1/(3*d10*d2^n2*d3^n3*d4^n4*d5^n5*d6*d7*d9) + 
        1/(3*d10*d2^n2*d3^n3*d4^n4*d5^n5*d7*d8*d9) + 
        d5^(1 - n5)/(3*d10*d2^n2*d3^n3*d4^n4*d6*d7^2*d8*(-1 + n5)) + 
        d5^(1 - n5)/(3*d10*d2^n2*d3^n3*d4^n4*d6*d8*d9^2*(-1 + n5)) + 
        (d2^(1 - n2)*d5^(1 - n5))/(3*d10*d3^n3*d4^n4*d6*d7*d8*d9^2*(-1 + n5)) - 
        d5^(1 - n5)/(3*d2^n2*d3^n3*d4^n4*d6*d7*d8*d9^2*(-1 + n5)) - 
        d5^(2 - n5)/(3*d10*d2^n2*d3^n3*d4^n4*d6*d7*d8*d9^2*(-1 + n5)) + 
        (2*d5^(1 - n5))/(3*d10*d2^n2*d3^n3*d4^n4*d6^2*d7*d9*(-1 + n5)) - 
        d5^(1 - n5)/(3*d10^2*d2^n2*d3^n3*d4^n4*d6*d7*d9*(-1 + n5)) + 
        (d4^(1 - n4)*d5^(1 - n5))/(3*d10*d2^n2*d3^n3*d6*d7^2*d8*d9*(-1 + n5)) - 
        (d3^(1 - n3)*d5^(1 - n5))/(3*d10*d2^n2*d4^n4*d6*d7^2*d8*d9*(-1 + n5)) - 
        d5^(2 - n5)/(3*d10*d2^n2*d3^n3*d4^n4*d6*d7^2*d8*d9*(-1 + n5)) + 
        (d3^(1 - n3)*d5^(1 - n5))/(3*d10^2*d2^n2*d4^n4*d6*d7*d8*d9*(-1 + n5)) - 
        (d2^(-1 - n2)*d5^(1 - n5)*n2)/(6*d10*d3^n3*d4^n4*d6*d7*d8*(-1 + n5)) - 
        (d2^(-1 - n2)*d5^(1 - n5)*n2)/(3*d10*d3^n3*d4^n4*d7*d8*d9*(-1 + n5)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*d5^(1 - n5)*n2)/(3*d10*d3^n3*d6*d7*d8*d9*(-1 + n5)) + 
        (d2^(-1 - n2)*d5^(1 - n5)*n2)/(6*d3^n3*d4^n4*d6*d7*d8*d9*(-1 + n5)) + 
        (d2^(-1 - n2)*d5^(1 - n5)*n2)/(2*d10*d3^n3*d4^n4*d6*d7*d8*d9*(-1 + n5)) + 
        (d3^(-1 - n3)*d5^(1 - n5)*n3)/(6*d10*d2^n2*d4^n4*d6*d7*d9*(-1 + n5)) - 
        (d3^(-1 - n3)*d5^(1 - n5)*n3)/(6*d2^n2*d4^n4*d6*d7*d8*d9*(-1 + n5)) + 
        (d3^(-1 - n3)*d5^(1 - n5)*n3)/(2*d10*d2^n2*d4^n4*d6*d7*d8*d9*(-1 + n5)) + 
        (d4^(-1 - n4)*d5^(1 - n5)*n4)/(2*d10*d2^n2*d3^n3*d6*d7*d9*(-1 + n5)) - 
        (d4^(-1 - n4)*d5^(1 - n5)*n4)/(6*d10*d2^n2*d3^n3*d6*d8*d9*(-1 + n5)) - 
        (d4^(-1 - n4)*d5^(1 - n5)*n4)/(3*d10*d2^n2*d3^n3*d7*d8*d9*(-1 + n5)) + 
        (d3^(1 - n3)*d4^(-1 - n4)*d5^(1 - n5)*n4)/(6*d10*d2^n2*d6*d7*d8*d9*(-1 + n5)) - 
        (d2^(1 - n2)*d4^(-1 - n4)*d5^(1 - n5)*n4)/(6*d10*d3^n3*d6*d7*d8*d9*(-1 + n5)) + 
        (d4^(-1 - n4)*d5^(1 - n5)*n4)/(2*d10*d2^n2*d3^n3*d6*d7*d8*d9*(-1 + n5)) - 
        (d2^(1 - n2)*d4^(-1 - n4)*d5^(2 - n5)*n4)/(2*d10*d3^n3*d6*d7*d8*d9^2*(-2 + n5)*(-1 + n5)) + 
        (d4^(-1 - n4)*d5^(2 - n5)*n4)/(2*d2^n2*d3^n3*d6*d7*d8*d9^2*(-2 + n5)*(-1 + n5)) - 
        (d4^(-1 - n4)*d5^(2 - n5)*n4)/(2*d10*d2^n2*d3^n3*d6*d7*d8*d9^2*(-2 + n5)*(-1 + n5)) - 
        (d2^(1 - n2)*d4^(-1 - n4)*d5^(2 - n5)*n4)/(2*d10*d3^n3*d6^2*d7*d8*d9*(-2 + n5)*(-1 + n5)) - 
        (d4^(-1 - n4)*d5^(2 - n5)*n4)/(2*d10*d2^n2*d3^n3*d6^2*d7*d8*d9*(-2 + n5)*(-1 + n5)) - 
        (d2^(-1 - n2)*d4^(-1 - n4)*d5^(2 - n5)*n2*n4)/(d10*d3^n3*d6*d7*d8*d9*(-2 + n5)*(-1 + n5)) + 
        (d5^(2 - n5)*(8 + 3*n4 - 4*n5))/(6*d10*d2^n2*d3^n3*d4^n4*d6^2*d7*d8*d9*(-2 + n5)*(-1 + n5)) + 
        (d5^(1 - n5)*rat(6 - d + 3*n2 + 3*n3 + 6*n4 - 6*n5, 6*(-1 + n5)))/(d10*d2^n2*d3^n3*d4^n4*d6*d7*d8*d9) + 
        (d4^(-1 - n4)*d5^(2 - n5)*rat(d*n4 - 2*n2*n4 - n4*n5, 2*(-2 + n5)*(-1 + n5)))/(d10*d2^n2*d3^n3*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n5 == 2 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n4 != 1
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?pos_/d3^n3?pos_/d4^n4?{>1}/d5^2/d6/d7/d8/d9/d10 =
        
        d2^(-1 - n2)/(d10*d3^n3*d4^n4*d5^2*d6*d8*d9) + 
        11/(6*d10*d2^n2*d3^n3*d4^n4*d5^2*d6*d8*d9) - 
        d2^(-1 - n2)/(2*d10*d3^n3*d4^n4*d5^2*d7*d8*d9) - 
        1/(3*d10*d2^n2*d3^n3*d4^n4*d5^2*d7*d8*d9) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(d10*d4^n4*d5^2*d6*d7*d8*d9) - 
        (11*d3^(1 - n3))/(6*d10*d2^n2*d4^n4*d5^2*d6*d7*d8*d9) + 
        d2^(1 - n2)/(3*d10*d3^n3*d4^n4*d5^2*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5^3*d6*d7*d8*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10^2*d3^n3*d5^2*d6*d7*d8*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5^2*d6^2*d7*d9*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5*d6*d7*d8^2*d9*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(2 - n4))/(d10*d3^n3*d5^2*d6*d7^2*d8*d9*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10*d3^n3*d5^2*d6*d7^2*d8*(-1 + n4)) + 
        d4^(1 - n4)/(3*d10*d2^n2*d3^n3*d5^2*d6*d7^2*d8*(-1 + n4)) + 
        (2*d4^(1 - n4))/(3*d10*d2^n2*d3^n3*d5^3*d6*d7*d8*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10*d3^n3*d5^2*d6*d8*d9^2*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d3^n3*d5^2*d6*d7*d8*d9^2*(-1 + n4)) + 
        (d2^(1 - n2)*d4^(1 - n4))/(3*d10*d3^n3*d5^2*d6*d7*d8*d9^2*(-1 + n4)) - 
        d4^(1 - n4)/(3*d2^n2*d3^n3*d5^2*d6*d7*d8*d9^2*(-1 + n4)) + 
        d4^(1 - n4)/(6*d10*d2^n2*d3^n3*d5^2*d6^2*d7*d9*(-1 + n4)) + 
        (2*d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5^3*d6*d7*d9*(-1 + n4)) - 
        (2*d4^(1 - n4))/(3*d10*d2^n2*d3^n3*d5^3*d6*d7*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10^2*d3^n3*d5^2*d6*d7*d9*(-1 + n4)) + 
        (2*d4^(1 - n4))/(d10^2*d2^n2*d3^n3*d5^2*d6*d7*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5^2*d7*d8^2*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*d4^(1 - n4))/(2*d10*d5^2*d6*d7*d8^2*d9*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(3*d10*d2^n2*d5^2*d6*d7*d8^2*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d3^n3*d5^2*d6*d7*d8^2*d9*(-1 + n4)) + 
        d4^(1 - n4)/(3*d2^n2*d3^n3*d5^2*d6*d7*d8^2*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5^3*d6*d8*d9*(-1 + n4)) - 
        (2*d4^(1 - n4))/(3*d10*d2^n2*d3^n3*d5^3*d6*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*d4^(1 - n4))/(d10*d5^2*d6*d7^2*d8*d9*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(3*d10*d2^n2*d5^2*d6*d7^2*d8*d9*(-1 + n4)) + 
        d4^(2 - n4)/(3*d10*d2^n2*d3^n3*d5^2*d6*d7^2*d8*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10*d3^n3*d5*d6*d7^2*d8*d9*(-1 + n4)) - 
        d4^(1 - n4)/(3*d10*d2^n2*d3^n3*d5*d6*d7^2*d8*d9*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5^3*d7*d8*d9*(-1 + n4)) + 
        (2*d4^(1 - n4))/(3*d10*d2^n2*d3^n3*d5^3*d7*d8*d9*(-1 + n4)) + 
        (d2^(1 - n2)*d4^(1 - n4))/(3*d10*d3^n3*d5^2*d6^2*d7*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(2 - n4))/(2*d10*d3^n3*d5^2*d6^2*d7*d8*d9*(-1 + n4)) - 
        d4^(2 - n4)/(3*d10*d2^n2*d3^n3*d5^2*d6^2*d7*d8*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*d4^(1 - n4))/(2*d10^2*d5^2*d6*d7*d8*d9*(-1 + n4)) - 
        (2*d3^(1 - n3)*d4^(1 - n4))/(d10^2*d2^n2*d5^2*d6*d7*d8*d9*(-1 + n4)) + 
        (d2^(-2 - n2)*d4^(1 - n4)*(-1 - n2))/(d3^n3*d5^2*d6*d7*d8*d9*(-1 + n4)) + 
        (d2^(-2 - n2)*d4^(2 - n4)*(-1 - n2))/(2*d10*d3^n3*d5^2*d6*d7*d8*d9*(-1 + n4)) + 
        (11*d2^(-1 - n2)*d4^(1 - n4)*n2)/(6*d10*d3^n3*d5^2*d6*d7*d8*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(1 - n4)*n2)/(d3^n3*d5*d6*d7*d8*d9^2*(-1 + n4)) + 
        (2*d4^(1 - n4)*n2)/(d10*d2^n2*d3^n3*d5*d6*d7*d8*d9^2*(-1 + n4)) + 
        (2*d2^(-1 - n2)*d4^(1 - n4)*n2)/(d10*d3^n3*d5^2*d7*d8*d9*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d5*d6^2*d7*d8*d9*(-1 + n4)) - 
        (11*d2^(-1 - n2)*d4^(1 - n4)*n2)/(6*d3^n3*d5^2*d6*d7*d8*d9*(-1 + n4)) + 
        (3*d2^(-1 - n2)*d4^(1 - n4)*n2)/(2*d10*d3^n3*d5^2*d6*d7*d8*d9*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d5^2*d6*d7*d8*d9*(-1 + n4)) + 
        (d2^(-2 - n2)*d4^(1 - n4)*(1 + n2))/(d10*d3^n3*d5^2*d6*d7*d8*(-1 + n4)) + 
        (d2^(-2 - n2)*d4^(1 - n4)*(1 + n2))/(2*d10*d3^n3*d5^2*d7*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*(1 + 2*n2))/(d10*d3^n3*d5*d6^2*d7*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*(1 + 4*n2))/(2*d10*d3^n3*d5*d6*d7*d8*d9^2*(-1 + n4)) + 
        (d4^(1 - n4)*(-1 + 12*n2))/(6*d10*d2^n2*d3^n3*d5*d6^2*d7*d8*d9*(-1 + n4)) + 
        (d2^(-2 - n2)*d4^(1 - n4)*(4*n2 + 4*n2^2))/(d10*d3^n3*d5*d6*d7*d8*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(2*d10*d5^2*d6*d7*d9*(-1 + n4)) - 
        (d3^(-1 - n3)*d4^(1 - n4)*n3)/(3*d10*d2^n2*d5^2*d6*d7*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d10*d5^2*d6*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(2*d5^2*d6*d7*d8*d9*(-1 + n4)) + 
        (d3^(-1 - n3)*d4^(1 - n4)*n3)/(3*d2^n2*d5^2*d6*d7*d8*d9*(-1 + n4)) - 
        (4*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d10*d2^n2*d5^2*d6*d7*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(2 - n4)*n3)/(d10*d5^2*d6*d7*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*rat(8*n2 - 2*d*n2 + 4*n2^2, -1 + n4))/(d10*d3^n3*d5*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(4 + 2*d + 3*n2 - 8*n3 - 3*n4, 2*(-1 + n4)))/(d10*d2^n2*d3^n3*d5^2*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n4 != 1 && n4 != 2
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?pos_/d3^n3?pos_/d4^n4?{>2}/d5/d6/d7/d8/d9/d10 =

        d2^(-1 - n2)/(d10*d3^n3*d4^n4*d5*d6*d8*d9) + 
        11/(6*d10*d2^n2*d3^n3*d4^n4*d5*d6*d8*d9) - 
        d2^(-1 - n2)/(2*d10*d3^n3*d4^n4*d5*d7*d8*d9) - 
        1/(3*d10*d2^n2*d3^n3*d4^n4*d5*d7*d8*d9) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(d10*d4^n4*d5*d6*d7*d8*d9) - 
        (11*d3^(1 - n3))/(6*d10*d2^n2*d4^n4*d5*d6*d7*d8*d9) + 
        d2^(1 - n2)/(3*d10*d3^n3*d4^n4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10^2*d3^n3*d5*d6*d7*d8*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5*d6^2*d7*d9*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d6*d7*d8^2*d9*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(2 - n4))/(d10*d3^n3*d5*d6*d7^2*d8*d9*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5^2*d7*d8*d9*(1 - n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10*d3^n3*d5*d6*d7^2*d8*(-1 + n4)) + 
        d4^(1 - n4)/(3*d10*d2^n2*d3^n3*d5*d6*d7^2*d8*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10*d3^n3*d5^2*d6*d7*d8*(-1 + n4)) + 
        d4^(1 - n4)/(3*d10*d2^n2*d3^n3*d5^2*d6*d7*d8*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10*d3^n3*d5*d6*d8*d9^2*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10*d3^n3*d6*d7*d8*d9^2*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d3^n3*d5*d6*d7*d8*d9^2*(-1 + n4)) + 
        (d2^(1 - n2)*d4^(1 - n4))/(3*d10*d3^n3*d5*d6*d7*d8*d9^2*(-1 + n4)) - 
        d4^(1 - n4)/(3*d2^n2*d3^n3*d5*d6*d7*d8*d9^2*(-1 + n4)) + 
        d4^(1 - n4)/(6*d10*d2^n2*d3^n3*d5*d6^2*d7*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5^2*d6*d7*d9*(-1 + n4)) - 
        d4^(1 - n4)/(3*d10*d2^n2*d3^n3*d5^2*d6*d7*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10^2*d3^n3*d5*d6*d7*d9*(-1 + n4)) + 
        (2*d4^(1 - n4))/(d10^2*d2^n2*d3^n3*d5*d6*d7*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d5*d7*d8^2*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*d4^(1 - n4))/(2*d10*d5*d6*d7*d8^2*d9*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(3*d10*d2^n2*d5*d6*d7*d8^2*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d3^n3*d5*d6*d7*d8^2*d9*(-1 + n4)) + 
        d4^(1 - n4)/(3*d2^n2*d3^n3*d5*d6*d7*d8^2*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10*d3^n3*d5^2*d6*d8*d9*(-1 + n4)) - 
        d4^(1 - n4)/(3*d10*d2^n2*d3^n3*d5^2*d6*d8*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4))/(2*d10*d3^n3*d6*d7^2*d8*d9*(-1 + n4)) - 
        d4^(1 - n4)/(3*d10*d2^n2*d3^n3*d6*d7^2*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*d4^(1 - n4))/(d10*d5*d6*d7^2*d8*d9*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(3*d10*d2^n2*d5*d6*d7^2*d8*d9*(-1 + n4)) + 
        d4^(2 - n4)/(3*d10*d2^n2*d3^n3*d5*d6*d7^2*d8*d9*(-1 + n4)) + 
        d4^(1 - n4)/(3*d10*d2^n2*d3^n3*d5^2*d7*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d10*d3^n3*d6^2*d7*d8*d9*(-1 + n4)) - 
        d4^(1 - n4)/(6*d10*d2^n2*d3^n3*d6^2*d7*d8*d9*(-1 + n4)) + 
        (d2^(1 - n2)*d4^(1 - n4))/(3*d10*d3^n3*d5*d6^2*d7*d8*d9*(-1 + n4)) - 
        d4^(2 - n4)/(3*d10*d2^n2*d3^n3*d5*d6^2*d7*d8*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*d4^(1 - n4))/(2*d10^2*d5*d6*d7*d8*d9*(-1 + n4)) - 
        (2*d3^(1 - n3)*d4^(1 - n4))/(d10^2*d2^n2*d5*d6*d7*d8*d9*(-1 + n4)) + 
        (d2^(-2 - n2)*d4^(1 - n4)*(-1 - n2))/(d3^n3*d5*d6*d7*d8*d9*(-1 + n4)) + 
        (11*d2^(-1 - n2)*d4^(1 - n4)*n2)/(6*d10*d3^n3*d5*d6*d7*d8*(-1 + n4)) + 
        (2*d2^(-1 - n2)*d4^(1 - n4)*n2)/(d10*d3^n3*d5*d6*d7*d9*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(1 - n4)*n2)/(d10*d3^n3*d6*d7*d8*d9*(-1 + n4)) - 
        (11*d2^(-1 - n2)*d4^(1 - n4)*n2)/(6*d3^n3*d5*d6*d7*d8*d9*(-1 + n4)) + 
        (7*d2^(-1 - n2)*d4^(1 - n4)*n2)/(2*d10*d3^n3*d5*d6*d7*d8*d9*(-1 + n4)) + 
        (d2^(-2 - n2)*d4^(1 - n4)*(1 + n2))/(d10*d3^n3*d5*d6*d7*d8*(-1 + n4)) + 
        (d2^(-2 - n2)*d4^(1 - n4)*(1 + n2))/(2*d10*d3^n3*d5*d7*d8*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(2*d10*d5*d6*d7*d9*(-1 + n4)) - 
        (d3^(-1 - n3)*d4^(1 - n4)*n3)/(3*d10*d2^n2*d5*d6*d7*d9*(-1 + n4)) - 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d10*d5*d6*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(2*d5*d6*d7*d8*d9*(-1 + n4)) + 
        (d3^(-1 - n3)*d4^(1 - n4)*n3)/(3*d2^n2*d5*d6*d7*d8*d9*(-1 + n4)) - 
        (4*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d10*d2^n2*d5*d6*d7*d8*d9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(2 - n4)*n3)/(d10*d5*d6*d7*d8*d9*(-1 + n4)) + 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d5*d6*d7^2*d8*(-2 + n4)*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d3^n3*d5*d6*d7*d8*d9^2*(-2 + n4)*(-1 + n4)) + 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d5*d6*d7*d8*d9^2*(-2 + n4)*(-1 + n4)) + 
        (2*d4^(2 - n4)*n2)/(d10*d2^n2*d3^n3*d5*d6*d7*d8*d9^2*(-2 + n4)*(-1 + n4)) + 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d5*d6^2*d7*d9*(-2 + n4)*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d5^2*d6*d7*d9*(-2 + n4)*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d6*d7^2*d8*d9*(-2 + n4)*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d3^(1 - n3)*d4^(2 - n4)*n2)/(d10*d5*d6*d7^2*d8*d9*(-2 + n4)*(-1 + n4)) + 
        (2*d2^(-1 - n2)*d4^(3 - n4)*n2)/(d10*d3^n3*d5*d6*d7^2*d8*d9*(-2 + n4)*(-1 + n4)) + 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d5^2*d7*d8*d9*(-2 + n4)*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d6^2*d7*d8*d9*(-2 + n4)*(-1 + n4)) - 
        (2*d2^(-1 - n2)*d4^(2 - n4)*n2)/(d10*d3^n3*d5^2*d6*d7*d8*d9*(-2 + n4)*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(2 - n4)*(-2 + 4*n2 + n4))/(2*d10*d3^n3*d5*d6^2*d7*d8*d9*(-2 + n4)*(-1 + n4)) + 
        (d2^(-2 - n2)*d4^(2 - n4)*(2 + 10*n2 + 8*n2^2 - n4 - n2*n4))/(2*d10*d3^n3*d5*d6*d7*d8*d9*(-2 + n4)*(-1 + n4)) + 
        (d4^(1 - n4)*rat(4 + 2*d + 3*n2 - 8*n3 - 3*n4, 2*(-1 + n4)))/(d10*d2^n2*d3^n3*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*d4^(2 - n4)*rat(2*n2 - 2*d*n2 + 4*n2^2 + 2*n2*n4, (-2 + n4)*(-1 + n4)))/(d10*d3^n3*d5*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n4 == 1 && n5 == 2 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n3 != 1
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?pos_/d3^n3?{>1}/d4/d5^2/d6/d7/d8/d9/d10 =

        -d2^(-1 - n2)/(4*d10*d3^n3*d4*d5^2*d6*d7*d9) - 
        1/(6*d10*d2^n2*d3^n3*d4*d5^2*d6*d7*d9) - 
        d2^(-1 - n2)/(4*d10*d3^n3*d4*d5^2*d6*d8*d9) + 
        d2^(-1 - n2)/(4*d10*d3^n3*d5^2*d6*d7*d8*d9) + 
        d2^(-1 - n2)/(4*d3^n3*d4*d5^2*d6*d7*d8*d9) + 
        1/(6*d2^n2*d3^n3*d4*d5^2*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5^2*d6*d7^2*d8*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d2^n2*d4*d5^2*d6*d7^2*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4*d5^3*d6*d7*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10^2*d4*d5^2*d6*d7*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5^2*d6*d8*d9^2*(-1 + n3)) - 
        d3^(1 - n3)/(12*d10*d2^n2*d4*d5^2*d6*d8*d9^2*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d4*d5^2*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d2^(1 - n2)*d3^(1 - n3))/(6*d10*d4*d5^2*d6*d7*d8*d9^2*(-1 + n3)) - 
        d3^(1 - n3)/(6*d2^n2*d4*d5^2*d6*d7*d8*d9^2*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5^2*d6^2*d7*d9*(-1 + n3)) + 
        d3^(1 - n3)/(12*d10*d2^n2*d4*d5^2*d6^2*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4*d5^3*d6*d7*d9*(-1 + n3)) - 
        d3^(1 - n3)/(3*d10*d2^n2*d4*d5^3*d6*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10^2*d4*d5^2*d6*d7*d9*(-1 + n3)) + 
        (7*d3^(1 - n3))/(12*d10^2*d2^n2*d4*d5^2*d6*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5^2*d7*d8^2*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d4*d5^2*d6*d7*d8^2*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3))/(4*d10*d4*d5^2*d6*d7*d8^2*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5*d6*d7*d8^2*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4*d5^3*d6*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5^2*d6*d8*d9*(-1 + n3)) + 
        (5*d3^(1 - n3))/(12*d10*d2^n2*d4^2*d5^2*d6*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d5^2*d6*d7^2*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d2^n2*d5^2*d6*d7^2*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3))/(4*d10*d4*d5^2*d6*d7^2*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(6*d10*d2^n2*d4*d5^2*d6*d7^2*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5*d6*d7^2*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d2^n2*d4*d5*d6*d7^2*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4*d5^3*d7*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(3*d10*d2^n2*d4*d5^3*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5^2*d7*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d2^n2*d4^2*d5^2*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d5^2*d6^2*d7*d8*d9*(-1 + n3)) + 
        (d2^(1 - n2)*d3^(1 - n3))/(6*d10*d4^2*d5^2*d6*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3))/(4*d10*d4^2*d5^2*d6*d7*d8*d9*(-1 + n3)) - 
        (5*d3^(2 - n3))/(12*d10*d2^n2*d4^2*d5^2*d6*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3))/(4*d10^2*d4*d5^2*d6*d7*d8*d9*(-1 + n3)) - 
        (7*d3^(2 - n3))/(12*d10^2*d2^n2*d4*d5^2*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(-1 - n2))/(4*d10*d5^2*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(-1 - n2))/(4*d4*d5^2*d6*d7*d8*d9*(-1 + n3)) + 
        (5*d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d10*d4*d5^2*d6*d7*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        (7*d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d10*d4*d5^2*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d5*d6^2*d7*d8*d9*(-1 + n3)) - 
        (7*d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d10*d5^2*d6*d7*d8*d9*(-1 + n3)) - 
        (5*d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d4*d5^2*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5^2*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(1 + n2))/(4*d10*d4*d5^2*d6*d7*d8*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(1 + n2))/(4*d10*d4*d5^2*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*(1 + 2*n2))/(4*d10*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*(1 + 2*n2))/(4*d10*d4*d5*d6^2*d7*d8*d9*(-1 + n3)) + 
        (d3^(1 - n3)*(-1 + 6*n2))/(12*d10*d2^n2*d4*d5*d6^2*d7*d8*d9*(-1 + n3)) + 
        (d3^(1 - n3)*(1 + 6*n2))/(12*d10*d2^n2*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(n2 + n2^2))/(d10*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*rat(4*n2 - d*n2 + 2*n2^2, 2*(-1 + n3)))/(d10*d4*d5*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(6 + d + 3*n2 - 6*n3, 6*(-1 + n3)))/(d10*d2^n2*d4*d5^2*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n4 == 2 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n3 != 1
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?pos_/d3^n3?{>1}/d4^2/d5/d6/d7/d8/d9/d10 =

        -d2^(-1 - n2)/(4*d10*d3^n3*d4^2*d5*d6*d7*d9) - 
        1/(6*d10*d2^n2*d3^n3*d4^2*d5*d6*d7*d9) - 
        d2^(-1 - n2)/(4*d10*d3^n3*d4^2*d5*d6*d8*d9) + 
        d2^(-1 - n2)/(4*d3^n3*d4^2*d5*d6*d7*d8*d9) + 
        1/(6*d2^n2*d3^n3*d4^2*d5*d6*d7*d8*d9) + 
        d2^(-1 - n2)/(4*d10*d3^n3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5*d6*d7^2*d8*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d2^n2*d4^2*d5*d6*d7^2*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5^2*d6*d7*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10^2*d4^2*d5*d6*d7*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5*d6*d8*d9^2*(-1 + n3)) - 
        d3^(1 - n3)/(12*d10*d2^n2*d4^2*d5*d6*d8*d9^2*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d6*d7*d8*d9^2*(-1 + n3)) + 
        d3^(1 - n3)/(12*d10*d2^n2*d4^2*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d4^2*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d2^(1 - n2)*d3^(1 - n3))/(6*d10*d4^2*d5*d6*d7*d8*d9^2*(-1 + n3)) - 
        d3^(1 - n3)/(6*d2^n2*d4^2*d5*d6*d7*d8*d9^2*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5*d6^2*d7*d9*(-1 + n3)) + 
        d3^(1 - n3)/(12*d10*d2^n2*d4^2*d5*d6^2*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5^2*d6*d7*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d2^n2*d4^2*d5^2*d6*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10^2*d4^2*d5*d6*d7*d9*(-1 + n3)) + 
        (7*d3^(1 - n3))/(12*d10^2*d2^n2*d4^2*d5*d6*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5*d7*d8^2*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d6*d7*d8^2*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d4^2*d5*d6*d7*d8^2*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3))/(4*d10*d4^2*d5*d6*d7*d8^2*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5^2*d6*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4^3*d5*d6*d8*d9*(-1 + n3)) + 
        (5*d3^(1 - n3))/(6*d10*d2^n2*d4^3*d5*d6*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d6*d7^2*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d2^n2*d4^2*d6*d7^2*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3))/(4*d10*d4^2*d5*d6*d7^2*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(6*d10*d2^n2*d4^2*d5*d6*d7^2*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5*d6*d7^2*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d2^n2*d4*d5*d6*d7^2*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5^2*d7*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d2^n2*d4^2*d5^2*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4^3*d5*d7*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(3*d10*d2^n2*d4^3*d5*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d6^2*d7*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(12*d10*d2^n2*d4^2*d6^2*d7*d8*d9*(-1 + n3)) + 
        (d2^(1 - n2)*d3^(1 - n3))/(3*d10*d4^3*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3))/(2*d10*d4^3*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (5*d3^(2 - n3))/(6*d10*d2^n2*d4^3*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3))/(4*d10^2*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (7*d3^(2 - n3))/(12*d10^2*d2^n2*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(-1 - n2))/(4*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5*d6*d7^2*d8*(-1 + n3)) + 
        (5*d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d10*d4^2*d5*d6*d7*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d3^(1 - n3)*n2)/(2*d10*d2^n2*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5*d6^2*d7*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5^2*d6*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4^2*d5*d6*d7*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d6*d7^2*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d5*d6*d7^2*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d10*d4*d5*d6*d7^2*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5^2*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d10*d4^2*d5*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d6^2*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4^2*d6*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5^2*d6*d7*d8*d9*(-1 + n3)) - 
        (5*d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(d10*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(1 + n2))/(4*d10*d4^2*d5*d6*d7*d8*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(1 + n2))/(4*d10*d4^2*d5*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*(1 + 2*n2))/(4*d10*d4*d5*d6^2*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(-1 + 3*n2 + 4*n2^2))/(4*d10*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*rat(23*n2 - 6*d*n2 + 12*n2^2, 12*(-1 + n3)))/(d10*d4*d5*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(9 + 2*d + 6*n2 - 12*n3, 12*(-1 + n3)))/(d10*d2^n2*d4^2*d5*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n4 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n3 != 1 && n3 != 2
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?pos_/d3^n3?{>2}/d4/d5/d6/d7/d8/d9/d10 =

        -d2^(-1 - n2)/(4*d10*d3^n3*d4*d5*d6*d7*d9) - 
        1/(6*d10*d2^n2*d3^n3*d4*d5*d6*d7*d9) - 
        d2^(-1 - n2)/(4*d10*d3^n3*d4*d5*d6*d8*d9) + 
        d2^(-1 - n2)/(4*d10*d3^n3*d5*d6*d7*d8*d9) + 
        d2^(-1 - n2)/(4*d3^n3*d4*d5*d6*d7*d8*d9) + 
        1/(6*d2^n2*d3^n3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5*d6*d7^2*d8*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d2^n2*d4*d5*d6*d7^2*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5^2*d6*d7*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10^2*d4*d5*d6*d7*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5*d6*d8*d9^2*(-1 + n3)) - 
        d3^(1 - n3)/(12*d10*d2^n2*d4*d5*d6*d8*d9^2*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d6*d7*d8*d9^2*(-1 + n3)) + 
        d3^(1 - n3)/(12*d10*d2^n2*d4*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        (d2^(1 - n2)*d3^(1 - n3))/(6*d10*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) - 
        d3^(1 - n3)/(6*d2^n2*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5*d6^2*d7*d9*(-1 + n3)) + 
        d3^(1 - n3)/(12*d10*d2^n2*d4*d5*d6^2*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5^2*d6*d7*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d2^n2*d4*d5^2*d6*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10^2*d4*d5*d6*d7*d9*(-1 + n3)) + 
        (7*d3^(1 - n3))/(12*d10^2*d2^n2*d4*d5*d6*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5*d7*d8^2*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d6*d7*d8^2*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d4*d5*d6*d7*d8^2*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5^2*d6*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5*d6*d8*d9*(-1 + n3)) + 
        (5*d3^(1 - n3))/(12*d10*d2^n2*d4^2*d5*d6*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d6*d7^2*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d2^n2*d4*d6*d7^2*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d5*d6*d7^2*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d2^n2*d5*d6*d7^2*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3))/(4*d10*d4*d5*d6*d7^2*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(6*d10*d2^n2*d4*d5*d6*d7^2*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d5^2*d7*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d2^n2*d4*d5^2*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4^2*d5*d7*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d2^n2*d4^2*d5*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d4*d6^2*d7*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(12*d10*d2^n2*d4*d6^2*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(4*d10*d5*d6^2*d7*d8*d9*(-1 + n3)) + 
        (d2^(1 - n2)*d3^(1 - n3))/(6*d10*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3))/(4*d10*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (5*d3^(2 - n3))/(12*d10*d2^n2*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(-1 - n2))/(4*d10*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(-1 - n2))/(4*d4*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d10*d4*d5*d6*d7*d8*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5*d6*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d10*d4*d5*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d10*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (5*d2^(-1 - n2)*d3^(1 - n3)*n2)/(12*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(d10*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(1 + n2))/(4*d10*d4*d5*d6*d7*d8*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(1 + n2))/(4*d10*d4*d5*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d10*d4*d5*d6*d7^2*d8*(-2 + n3)*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d10^2*d4*d5*d6*d7*d8*(-2 + n3)*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d4*d5*d6*d7*d8*d9^2*(-2 + n3)*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d10*d4*d5*d6*d7*d8*d9^2*(-2 + n3)*(-1 + n3)) - 
        (d3^(2 - n3)*n2)/(2*d10*d2^n2*d4*d5*d6*d7*d8*d9^2*(-2 + n3)*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d10*d4*d5^2*d6*d7*d9*(-2 + n3)*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d10*d4*d6*d7^2*d8*d9*(-2 + n3)*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d10*d5*d6*d7^2*d8*d9*(-2 + n3)*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(3 - n3)*n2)/(2*d10*d4*d5*d6*d7^2*d8*d9*(-2 + n3)*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d10*d4*d5^2*d7*d8*d9*(-2 + n3)*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3)*n2)/(2*d10*d4*d5^2*d6*d7*d8*d9*(-2 + n3)*(-1 + n3)) + 
        (d3^(2 - n3)*(14 + 6*n2 - 7*n3))/(12*d10^2*d2^n2*d4*d5*d6*d7*d8*d9*(-2 + n3)*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3)*(2 + 2*n2 - n3))/(4*d10^2*d4*d5*d6*d7*d8*d9*(-2 + n3)*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3)*(-2 + 4*n2 + n3))/(4*d10*d4*d5*d6*d7*d8^2*d9*(-2 + n3)*(-1 + n3)) + 
        (d3^(1 - n3)*rat(6 + d + 3*n2 - 6*n3, 6*(-1 + n3)))/(d10*d2^n2*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*d3^(2 - n3)*rat(n2 - d*n2 + n2*n3, 2*(-2 + n3)*(-1 + n3)))/(d10*d4*d5*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n3 == 1 && n4 == 1 && n5 == 2 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?pos_/d3/d4/d5^2/d6/d7/d8/d9/d10 =

        d2^(1 - n2)/(3*d10*d3*d4*d5^2*d6*d7^2*d8) - 
        1/(2*d10*d2^n2*d3*d4*d5^2*d6*d7^2*d8) + 
        2/(3*d10*d2^n2*d3*d4*d5*d6*d7^2*d8) - 
        (4*d2^(1 - n2))/(3*d10*d3*d4*d5^3*d6*d7*d8) - 
        1/(d10*d2^n2*d3*d4*d5^3*d6*d7*d8) - 
        1/(2*d10^2*d2^n2*d3*d4*d5^2*d6*d7*d8) - 
        d2^(1 - n2)/(2*d10*d3*d4*d5^2*d6*d8*d9^2) - 
        3/(2*d10*d2^n2*d3*d4*d5^2*d6*d8*d9^2) - 
        2/(3*d10*d2^n2*d3*d4*d5*d6*d8*d9^2) + 
        2/(3*d10*d2^n2*d3*d4*d6*d7*d8*d9^2) - 
        d2^(1 - n2)/(3*d3*d4*d5^2*d6*d7*d8*d9^2) + 
        d2^(2 - n2)/(3*d10*d3*d4*d5^2*d6*d7*d8*d9^2) - 
        1/(2*d2^n2*d3*d4*d5^2*d6*d7*d8*d9^2) + 
        (5*d2^(1 - n2))/(3*d10*d3*d4*d5*d6*d7*d8*d9^2) - 
        7/(6*d2^n2*d3*d4*d5*d6*d7*d8*d9^2) + 
        d2^(1 - n2)/(6*d10*d3*d4*d5^2*d6^2*d7*d9) + 
        3/(2*d10*d2^n2*d3*d4*d5^2*d6^2*d7*d9) + 
        4/(3*d10*d2^n2*d3*d4*d5*d6^2*d7*d9) - 
        (2*d2^(1 - n2))/(3*d10*d3*d4*d5^3*d6*d7*d9) - 
        1/(d10*d2^n2*d3*d4*d5^3*d6*d7*d9) - 
        d2^(1 - n2)/(3*d10*d3^2*d4*d5^2*d6*d7*d9) + 
        1/(2*d10*d2^n2*d3^2*d4*d5^2*d6*d7*d9) + 
        d2^(1 - n2)/(2*d10^2*d3*d4*d5^2*d6*d7*d9) + 
        3/(2*d10^2*d2^n2*d3*d4*d5^2*d6*d7*d9) - 
        7/(6*d10*d2^n2*d3*d4*d5^2*d6*d7*d9) - 
        2/(3*d10*d2^n2*d3^2*d4*d5*d6*d7*d9) + 
        2/(3*d10^2*d2^n2*d3*d4*d5*d6*d7*d9) + 
        1/(2*d10*d2^n2*d3*d4*d5^2*d7*d8^2*d9) + 
        (2*d2^(1 - n2))/(3*d10*d4*d5^2*d6*d7*d8^2*d9) + 
        1/(2*d10*d2^n2*d4*d5^2*d6*d7*d8^2*d9) - 
        (2*d2^(1 - n2))/(3*d3*d4*d5^2*d6*d7*d8^2*d9) - 
        1/(2*d2^n2*d3*d4*d5^2*d6*d7*d8^2*d9) + 
        7/(6*d10*d2^n2*d4*d5*d6*d7*d8^2*d9) - 
        7/(6*d2^n2*d3*d4*d5*d6*d7*d8^2*d9) + 
        (4*d2^(1 - n2))/(3*d10*d3*d4*d5^3*d6*d8*d9) + 
        1/(d10*d2^n2*d3*d4*d5^3*d6*d8*d9) - 
        d2^(1 - n2)/(6*d10*d3*d4^2*d5^2*d6*d8*d9) - 
        3/(2*d10*d2^n2*d3*d4^2*d5^2*d6*d8*d9) - 
        1/(2*d10*d2^n2*d3^2*d4*d5^2*d6*d8*d9) + 
        7/(6*d10*d2^n2*d3*d4*d5^2*d6*d8*d9) - 
        4/(3*d10*d2^n2*d3*d4^2*d5*d6*d8*d9) - 
        2/(3*d10*d2^n2*d3*d4*d6*d7^2*d8*d9) + 
        d2^(1 - n2)/(3*d10*d3*d5^2*d6*d7^2*d8*d9) + 
        1/(2*d10*d2^n2*d3*d5^2*d6*d7^2*d8*d9) - 
        d2^(1 - n2)/(3*d10*d4*d5^2*d6*d7^2*d8*d9) - 
        1/(2*d10*d2^n2*d4*d5^2*d6*d7^2*d8*d9) + 
        7/(6*d10*d2^n2*d3*d5*d6*d7^2*d8*d9) - 
        7/(6*d10*d2^n2*d4*d5*d6*d7^2*d8*d9) - 
        d2^(1 - n2)/(3*d10*d3*d4*d5*d6*d7^2*d8*d9) + 
        (2*d2^(1 - n2))/(3*d10*d3*d4*d5^3*d7*d8*d9) + 
        1/(d10*d2^n2*d3*d4*d5^3*d7*d8*d9) - 
        d2^(1 - n2)/(3*d10*d3*d4^2*d5^2*d7*d8*d9) + 
        1/(2*d10*d2^n2*d3*d4^2*d5^2*d7*d8*d9) - 
        2/(3*d10*d2^n2*d3*d4^2*d5*d7*d8*d9) - 
        4/(3*d10*d2^n2*d3*d4*d6^2*d7*d8*d9) + 
        (2*d2^(1 - n2))/(3*d10*d3*d5^2*d6^2*d7*d8*d9) + 
        1/(2*d10*d2^n2*d3*d5^2*d6^2*d7*d8*d9) - 
        (2*d2^(2 - n2))/(3*d10*d3*d4*d5^2*d6^2*d7*d8*d9) + 
        7/(6*d10*d2^n2*d3*d5*d6^2*d7*d8*d9) - 
        (4*d2^(1 - n2))/(3*d10*d3*d4*d5*d6^2*d7*d8*d9) + 
        1/(2*d10*d2^n2*d3^2*d5^2*d6*d7*d8*d9) + 
        d2^(1 - n2)/(6*d10*d4^2*d5^2*d6*d7*d8*d9) + 
        3/(2*d10*d2^n2*d4^2*d5^2*d6*d7*d8*d9) + 
        d2^(2 - n2)/(3*d10*d3*d4^2*d5^2*d6*d7*d8*d9) - 
        d2^(1 - n2)/(2*d10^2*d4*d5^2*d6*d7*d8*d9) - 
        3/(2*d10^2*d2^n2*d4*d5^2*d6*d7*d8*d9) + 
        d2^(1 - n2)/(3*d3^2*d4*d5^2*d6*d7*d8*d9) - 
        1/(2*d2^n2*d3^2*d4*d5^2*d6*d7*d8*d9) + 
        4/(3*d10*d2^n2*d4^2*d5*d6*d7*d8*d9) + 
        (2*d2^(1 - n2))/(3*d10*d3*d4^2*d5*d6*d7*d8*d9) + 
        2/(d10*d2^n2*d3*d4^2*d5*d6*d7*d8*d9) - 
        2/(3*d10^2*d2^n2*d4*d5*d6*d7*d8*d9) + 
        2/(3*d2^n2*d3^2*d4*d5*d6*d7*d8*d9) + 
        1/(d10*d2^n2*d3^2*d4*d5*d6*d7*d8*d9) + 
        (-6 - n2)/(6*d10*d2^n2*d3*d4*d5^2*d6*d7*d8) + 
        (1 - n2)/(2*d10*d2^n2*d3*d5^2*d6*d7*d8*d9) + 
        (-1 + n2)/(6*d2^n2*d3*d4*d5^2*d6*d7*d8*d9) - 
        (3*d2^(-1 - n2)*n2)/(2*d10*d3*d4*d5^2*d6*d7*d8) - 
        (4*d2^(-1 - n2)*n2)/(3*d10*d3*d4*d5*d6*d7*d8) + 
        (3*d2^(-1 - n2)*n2)/(2*d10*d3*d4*d5^2*d7*d8*d9) + 
        (2*d2^(-1 - n2)*n2)/(3*d10*d3*d4*d5*d7*d8*d9) - 
        (3*d2^(-1 - n2)*n2)/(2*d10*d3*d5^2*d6*d7*d8*d9) + 
        (3*d2^(-1 - n2)*n2)/(2*d3*d4*d5^2*d6*d7*d8*d9) - 
        (2*d2^(-1 - n2)*n2)/(3*d10*d3*d5*d6*d7*d8*d9) + 
        (4*d2^(-1 - n2)*n2)/(3*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*n2)/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (4 + 3*n2)/(6*d10*d2^n2*d3*d4*d5^2*d7*d8*d9) + 
        rat(-6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d6*d7*d8^2) + 
        rat(-6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d6^2*d7*d8) + 
        rat(-6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d6*d8*d9^2) + 
        rat(-6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d7*d8*d9^2) + 
        rat(-6, -4 + d)/(d10*d2^n2*d3*d4*d5^2*d6^2*d7^2*d9) + 
        rat(-6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d6^2*d7*d9) + 
        rat(-6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d7^2*d8*d9) + 
        rat(-3, -4 + d)/(d10*d2^n2*d3*d4*d5^2*d6^2*d7*d8^2) + 
        rat(-3, -4 + d)/(d10*d2^n2*d3*d4*d5^2*d7^2*d8*d9^2) + 
        rat(-3, -4 + d)/(d10*d2^n2*d3*d4^2*d5^2*d7*d8*d9^2) + 
        rat(-3, -4 + d)/(d10*d2^n2*d3*d4^2*d5^2*d6*d7^2*d9) + 
        rat(-3, -4 + d)/(d10*d2^n2*d3^2*d4*d5^2*d6*d7^2*d9) + 
        rat(-3, -4 + d)/(d10*d2^n2*d3^2*d4*d5^2*d6^2*d7*d9) + 
        rat(-3, -4 + d)/(d10*d2^n2*d3^2*d4^2*d5^2*d6*d7*d9) + 
        rat(-3, -4 + d)/(d2^n2*d3*d4^2*d5^2*d6*d7*d8^2*d9) + 
        rat(-3, -4 + d)/(d10*d2^n2*d3^2*d5^2*d6*d7^2*d8*d9) + 
        rat(-3, -4 + d)/(d10*d2^n2*d3^2*d5^2*d6^2*d7*d8*d9) + 
        rat(-3, -4 + d)/(d2^n2*d3*d4^2*d5^2*d6^2*d7*d8*d9) + 
        rat(3, -4 + d)/(d10*d2^n2*d3*d4^2*d5^2*d6^2*d7*d8) + 
        rat(3, -4 + d)/(d10*d2^n2*d3*d4*d5^2*d6*d7^2*d9^2) + 
        rat(3, -4 + d)/(d2^n2*d3*d4^2*d5^2*d6*d7*d8*d9^2) + 
        rat(3, -4 + d)/(d10*d2^n2*d3*d4*d5^2*d6^2*d8^2*d9) + 
        rat(3, -4 + d)/(d10*d2^n2*d3*d4^2*d5^2*d6*d8^2*d9) + 
        rat(3, -4 + d)/(d10*d2^n2*d3^2*d4*d5^2*d6^2*d8*d9) + 
        rat(3, -4 + d)/(d10*d2^n2*d3^2*d4^2*d5^2*d6*d8*d9) + 
        rat(3, -4 + d)/(d2^n2*d3*d4^2*d5^2*d6*d7^2*d8*d9) + 
        rat(3, -4 + d)/(d2^n2*d3^2*d4*d5^2*d6*d7^2*d8*d9) + 
        rat(3, -4 + d)/(d2^n2*d3^2*d4*d5^2*d6^2*d7*d8*d9) + 
        rat(3, -4 + d)/(d2^n2*d3^2*d4^2*d5^2*d6*d7*d8*d9) + 
        rat(6, -4 + d)/(d10*d2^n2*d3*d4*d5^2*d6^2*d7^2*d8) + 
        rat(6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d6*d7^2*d8) + 
        rat(6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d6*d7*d9^2) + 
        rat(6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d6*d7^2*d9) + 
        rat(6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d6*d8^2*d9) + 
        rat(6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d7*d8^2*d9) + 
        rat(6, -4 + d)/(d10*d2^n2*d3*d4*d5^3*d6^2*d8*d9) + 
        (d2^(-1 - n2)*rat(-3*n2, -4 + d))/(d10*d3*d4*d5^2*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(-3*n2, -4 + d))/(d3*d4*d5^2*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(-3*n2, -4 + d))/(d10*d3*d4^2*d5^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-3*n2, -4 + d))/(d3*d4*d5^2*d6^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-3*n2, -4 + d))/(d3*d4^2*d5^2*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(3*n2, -4 + d))/(d10*d3*d4*d5^2*d6*d7^2*d8) + 
        (d2^(-1 - n2)*rat(3*n2, -4 + d))/(d10*d3*d4*d5^2*d6^2*d7*d8) + 
        (d2^(-1 - n2)*rat(3*n2, -4 + d))/(d10*d3*d4^2*d5^2*d6*d7*d8) + 
        (d2^(-1 - n2)*rat(3*n2, -4 + d))/(d10*d3*d5^2*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(3*n2, -4 + d))/(d10*d3*d5^2*d6^2*d7*d8*d9) + 
        rat(6 - 2*d + 3*n2, 3)/(d10*d2^n2*d3*d4*d5*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n3 == 1 && n4 == 2 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n2 != 1
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?{>1}/d3/d4^2/d5/d6/d7/d8/d9/d10 =

        (2*d2^(1 - n2))/(3*d3*d4^2*d5*d6*d7*d8*d9^2) - 
        (2*d2^(2 - n2))/(3*d10*d3*d4^2*d5*d6*d7*d8*d9^2) + 
        (2*d2^(1 - n2))/(3*d10*d3*d4^2*d5^2*d6*d7*d9) + 
        (2*d2^(1 - n2))/(3*d10*d4^2*d5*d6*d7^2*d8*d9) - 
        (2*d2^(1 - n2))/(3*d10*d3*d4*d5*d6*d7^2*d8*d9) - 
        (2*d2^(1 - n2))/(3*d10*d3*d4^2*d5^2*d7*d8*d9) - 
        d2^(1 - n2)/(d10*d3*d4^2*d5*d6*d7*d8*d9) + 
        d2^(1 - n2)/(d10*d3*d4^2*d5*d6^2*d7*d9^2*(-1 + n2)) + 
        (d2^(1 - n2)*rat(-4, -5 + d))/(d10^2*d3*d4^2*d5*d6*d7*d9^2) + 
        (d2^(1 - n2)*rat(-2, -5 + d))/(d10*d3*d4^2*d5*d6*d7^2*d9^2) + 
        (d2^(1 - n2)*rat(-2, -5 + d))/(d10*d3*d4^2*d5^2*d6*d7*d9^2) + 
        (d2^(1 - n2)*rat(-2, -5 + d))/(d10*d3*d4^2*d5^2*d6*d7^2*d9) + 
        (d2^(1 - n2)*rat(2, -5 + d))/(d10*d3*d4^2*d5*d7^2*d8*d9^2) + 
        (d2^(1 - n2)*rat(2, -5 + d))/(d10*d3*d4^2*d5^2*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(2, -5 + d))/(d10*d3*d4^2*d5^2*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(4, -5 + d))/(d10^2*d4^2*d5*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(30 - 2*d - 20*n2, 3*(-5 + d)*(-1 + n2)))/(d3*d4^3*d5*d6*d7*d8*d9^2) + 
        rat(6 + 2*d - 16*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4^2*d5*d6*d7*d9^2) + 
        rat(-6 + 4*d - 14*n2, 3*(-5 + d))/(d10^2*d2^n2*d3*d4^2*d5*d6*d7*d9) + 
        rat(-12 + 4*d - 8*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4^2*d5*d6*d7^2*d8) + 
        (d2^(1 - n2)*rat(-12 + 4*d - 8*n2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d5*d6*d8*d9^3) + 
        (d2^(1 - n2)*rat(9 - d - 4*n2, 3*(-5 + d)*(-1 + n2)))/(d3^2*d4^2*d5*d6*d7*d8*d9^2) + 
        rat(-6 + 2*d - 4*n2, 3*(-5 + d))/(d2^n2*d3*d4*d5*d6*d7*d8*d9^2) + 
        rat(-6 + 2*d - 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4*d5^2*d6*d7*d9) + 
        rat(-6 + 2*d - 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4*d6*d7^2*d8*d9) + 
        rat(-6 + 2*d - 4*n2, 3*(-5 + d))/(d10*d2^n2*d4*d5*d6*d7^2*d8*d9) + 
        rat(-6 + 2*d - 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4*d6^2*d7*d8*d9) + 
        rat(-6 + 2*d - 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4^2*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-6 + 2*d - 4*n2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^3*d5*d6^2*d7*d8) + 
        (d2^(1 - n2)*rat(-6 + 2*d - 4*n2, 3*(-5 + d)*(-1 + n2)))/(d10^2*d3*d4^2*d5^2*d6*d7*d9) + 
        (d2^(1 - n2)*rat(2 - 2*n2, -5 + d))/(d10*d3*d4^2*d6*d7*d8*d9^2) + 
        rat(2 - 2*n2, -5 + d)/(d2^n2*d3*d4^2*d5*d6*d7*d8*d9^2) + 
        rat(2 - 2*n2, -5 + d)/(d10*d2^n2*d3^2*d4*d5*d6*d7*d8*d9) + 
        rat(12 - 2*d - 2*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4^2*d5*d6*d8*d9^2) + 
        (d2^(1 - n2)*rat(12 - 2*d - 2*n2, (-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5*d6*d7^2*d8) + 
        (d2^(1 - n2)*rat(12 - 2*d - 2*n2, (-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^3*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(12 - 2*d - 2*n2, (-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-3 + d - 2*n2, (-5 + d)*(-4 + d - n2)))/(d10*d3*d4^2*d5*d6^2*d7^2*d9) + 
        (d2^(1 - n2)*rat(-3 + d - 2*n2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d5*d6^2*d7^2*d8) + 
        rat(-18 + 4*d - 2*n2, 3*(-5 + d))/(d2^n2*d3^2*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(6 - d - n2, (-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^2*d5*d6*d7^2*d8) + 
        (d2^(1 - n2)*rat(-6 + d + n2, (-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^2*d6*d7^2*d8*d9) + 
        rat(-2 + 2*n2, -5 + d)/(d10^2*d2^n2*d3*d4^2*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-2 + 2*n2, -5 + d))/(d10*d3*d4^2*d5*d6*d8*d9^2) + 
        rat(3 - d + 2*n2, -5 + d)/(d2^n2*d3*d4^2*d5*d6*d7*d8^2*d9) + 
        rat(-12 + 2*d + 2*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4^2*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-12 + 2*d + 2*n2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d5^2*d6*d7^2*d8) + 
        (d2^(1 - n2)*rat(-12 + 2*d + 2*n2, (-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-12 + 2*d + 2*n2, (-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-12 + 2*d + 2*n2, (-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^3*d6*d7*d8*d9) + 
        rat(6 - 2*d + 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4*d5*d6*d7^2*d8) + 
        rat(6 - 2*d + 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4^2*d5*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(6 - 2*d + 4*n2, 3*(-5 + d)))/(d10*d3*d4*d5*d6*d7*d8*d9^2) + 
        rat(6 - 2*d + 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4*d5*d6^2*d7*d9) + 
        rat(6 - 2*d + 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4^2*d5*d6*d7*d9) + 
        rat(6 - 2*d + 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d5*d6*d7^2*d8*d9) + 
        rat(6 - 2*d + 4*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4*d5^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(6 - 2*d + 4*n2, 3*(-5 + d)*(-1 + n2)))/(d3*d4^3*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(6 - 2*d + 4*n2, 3*(-5 + d)*(-1 + n2)))/(d10^2*d4^2*d5^2*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-9 + d + 4*n2, 3*(-5 + d)*(-1 + n2)))/(d10*d3^2*d4^2*d5*d6*d7*d9^2) + 
        (d2^(1 - n2)*rat(-24 + 4*d + 4*n2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^3*d5^2*d6*d7*d8) + 
        rat(12 - 4*d + 8*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4^3*d5*d7*d8*d9) + 
        (d2^(1 - n2)*rat(12 - 4*d + 8*n2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d6*d7*d8*d9^3) + 
        (d2^(1 - n2)*rat(-18 + 2*d + 8*n2, 3*(-5 + d)*(-1 + n2)))/(d3*d4^3*d5*d6*d7*d8^2*d9) + 
        rat(-2*d + 10*n2, 3*(-5 + d))/(d10*d2^n2*d3*d4^2*d5^2*d6*d7*d8) + 
        rat(6 - 4*d + 14*n2, 3*(-5 + d))/(d10^2*d2^n2*d4^2*d5*d6*d7*d8*d9) + 
        rat(-6 - 2*d + 16*n2, 3*(-5 + d))/(d10*d2^n2*d4^2*d5*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-30 + 2*d + 20*n2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^3*d5*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-42 + 2*d + 32*n2, 3*(-5 + d)*(-1 + n2)))/(d10*d4^3*d5*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-60 + 28*d - 4*d^2 - 20*n2 + 12*d*n2 - 20*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10^2*d4^2*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-30 + 22*d - 4*d^2 - 50*n2 + 18*d*n2 - 20*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10^2*d3^2*d4^2*d5*d6*d7*d9) + 
        (d2^(2 - n2)*rat(10 - 6*d + 16*n2 + 4*d*n2 - 16*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^3*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-36 + 14*d - 2*d^2 + 2*n2 + 6*d*n2 - 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5*d6*d7^2*d9) + 
        (d2^(1 - n2)*rat(-36 + 14*d - 2*d^2 + 2*n2 + 6*d*n2 - 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^3*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(-66 + 30*d - 4*d^2 - 18*n2 + 10*d*n2 - 16*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d4^3*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-6 + 8*d - 2*d^2 - 28*n2 + 12*d*n2 - 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d5*d6^3*d7*d9) + 
        (d2^(1 - n2)*rat(174 - 58*d + 4*d^2 - 58*n2 + 18*d*n2 - 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5*d6^2*d7*d9) + 
        rat(-84 + 34*d - 4*d^2 - 2*n2 + 6*d*n2 - 14*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d4^2*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-54 + 28*d - 4*d^2 - 32*n2 + 12*d*n2 - 14*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d5^2*d6*d8*d9^2) + 
        rat(-9 + 19*d - 4*d^2 - 77*n2 + 21*d*n2 - 14*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3^2*d4^2*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(34 - 10*d - 12*n2 + 8*d*n2 - 12*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^3*d5*d6*d8*d9) + 
        rat(-15 + 11*d - 2*d^2 - 25*n2 + 9*d*n2 - 10*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d2^n2*d3*d4^2*d5*d6^2*d7*d8*d9) + 
        rat(-75 + 38*d - 5*d^2 - 40*n2 + 12*d*n2 - 10*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d6*d7^2*d8*d9) + 
        rat(90 - 25*d + d^2 - 55*n2 + 15*d*n2 - 10*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d5*d6^2*d7*d9) + 
        (d2^(1 - n2)*rat(42 - 20*d + 2*d^2 + 16*n2 - 8*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d4^4*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(42 - 20*d + 2*d^2 + 16*n2 - 8*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10^2*d4^3*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-5 - d + 18*n2 - 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d4^2*d5*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-5 - d + 18*n2 - 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d5^2*d6*d8*d9) + 
        (d2^(1 - n2)*rat(-5 - d + 18*n2 - 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(5 - 3*d + 8*n2 + 2*d*n2 - 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d5*d6*d7^2*d8) + 
        (d2^(1 - n2)*rat(5 - 3*d + 8*n2 + 2*d*n2 - 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d3^2*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-18 + 7*d - d^2 + n2 + 3*d*n2 - 8*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^2*d5*d6*d7^2*d9) + 
        (d2^(1 - n2)*rat(-78 + 34*d - 4*d^2 - 14*n2 + 6*d*n2 - 8*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5*d6^2*d8*d9) + 
        (d2^(1 - n2)*rat(-3 + 4*d - d^2 - 14*n2 + 6*d*n2 - 8*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10^2*d3*d4^2*d5*d6^2*d7*d9) + 
        (d2^(1 - n2)*rat(-18 + 12*d - 2*d^2 - 24*n2 + 8*d*n2 - 8*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5*d6*d8^2*d9) + 
        rat(132 - 38*d + 2*d^2 - 74*n2 + 18*d*n2 - 8*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^3*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(-2*d + 16*n2 - 6*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10^2*d3*d4^2*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(11 - 3*d - 10*n2 + 4*d*n2 - 6*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d6^2*d7*d8*d9) + 
        rat(17 - 5*d - 6*n2 + 4*d*n2 - 6*n2^2, 3*(-5 + d))/(d10*d2^n2*d3*d4^2*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-54 + 20*d - 2*d^2 + 8*n2 - 4*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5^2*d6*d8*d9) + 
        rat(-24 + 9*d - d^2 + 3*n2 + d*n2 - 4*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3^2*d4^2*d5*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-24 + 9*d - d^2 + 3*n2 + d*n2 - 4*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^2*d6*d7*d8*d9^2) + 
        rat(6 - 2*d - 2*n2 + 2*d*n2 - 4*n2^2, (-5 + d)*(-4 + d - n2))/(d2^n2*d3*d4^3*d5*d6*d7*d8*d9) + 
        rat(-44 + 18*d - 2*d^2 - 2*n2 + 2*d*n2 - 4*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^3*d5*d6*d8*d9) + 
        (d2^(1 - n2)*rat(-39 + 17*d - 2*d^2 - 7*n2 + 3*d*n2 - 4*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d5*d6^2*d8*d9^2) + 
        rat(-14 + 7*d - d^2 - 7*n2 + 3*d*n2 - 4*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d5*d6^2*d8*d9) + 
        rat(-9 + 6*d - d^2 - 12*n2 + 4*d*n2 - 4*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d5*d6*d7*d8^2) + 
        (d2^(1 - n2)*rat(-9 + 6*d - d^2 - 12*n2 + 4*d*n2 - 4*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^2*d5*d6^2*d7*d9) + 
        (d2^(1 - n2)*rat(-9 + 6*d - d^2 - 12*n2 + 4*d*n2 - 4*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d6^2*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(156 - 52*d + 4*d^2 - 52*n2 + 12*d*n2 - 4*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(2*n2 - 2*n2^2, -5 + d))/(d10*d3*d4^2*d5*d6*d7*d8) + 
        rat(-2 + 4*n2 - 2*n2^2, -5 + d)/(d10*d2^n2*d3*d4^2*d5*d7*d8*d9) + 
        rat(-32 + 11*d - d^2 + 9*n2 - d*n2 - 2*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d5*d7*d8^2*d9) + 
        rat(27 - 10*d + d^2 - 4*n2 + 2*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d2^n2*d3*d4^2*d5*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(-2*n2 + 2*n2^2, -5 + d))/(d3*d4^2*d5*d6*d7*d8*d9) + 
        rat(-108 + 32*d - 2*d^2 + 56*n2 - 12*d*n2 + 2*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d5^2*d6*d7*d9) + 
        rat(7 - 6*d + d^2 + 16*n2 - 4*d*n2 + 2*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d5^2*d6*d8*d9) + 
        (d2^(-1 - n2)*rat(8*n2 - 2*d*n2 + 2*n2^2, -5 + d))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        rat(32 - 11*d + d^2 - 9*n2 + d*n2 + 2*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d6*d7*d8^2*d9) + 
        rat(32 - 11*d + d^2 - 9*n2 + d*n2 + 2*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3^2*d4^2*d5*d6*d8*d9) + 
        (d2^(-1 - n2)*rat(-12*n2 + 2*d*n2 + 2*n2^2, 3*(-5 + d)))/(d10*d3*d4^2*d5*d7*d8*d9) + 
        rat(132 - 46*d + 4*d^2 - 34*n2 + 6*d*n2 + 2*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d5^2*d7*d8*d9) + 
        rat(-21 + 10*d - d^2 - 8*n2 + 4*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d2^n2*d3*d4^2*d5^2*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-21 + 10*d - d^2 - 8*n2 + 4*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d6^2*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-21 + 10*d - d^2 - 8*n2 + 4*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d6^2*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-156 + 52*d - 4*d^2 + 52*n2 - 12*d*n2 + 4*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5^2*d6*d7*d9) + 
        rat(-21 + 5*d + 17*n2 - 5*d*n2 + 4*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3^2*d4^2*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(9 - 6*d + d^2 + 12*n2 - 4*d*n2 + 4*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d5^2*d6^2*d7*d9) + 
        (d2^(1 - n2)*rat(9 - 6*d + d^2 + 12*n2 - 4*d*n2 + 4*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^2*d6^2*d7*d8*d9) + 
        rat(14 - 7*d + d^2 + 7*n2 - 3*d*n2 + 4*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d4^2*d5*d6^2*d7*d8*d9) + 
        rat(-6 + 2*d + 2*n2 - 2*d*n2 + 4*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^3*d5*d6*d7*d8) + 
        rat(44 - 18*d + 2*d^2 + 2*n2 - 2*d*n2 + 4*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d4^3*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(24 - 9*d + d^2 - 3*n2 - d*n2 + 4*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^2*d5*d6*d8*d9^2) + 
        (d2^(1 - n2)*rat(2*d - 16*n2 + 6*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10^2*d4^2*d5*d6*d7*d8*d9) + 
        rat(21 - 13*d + 2*d^2 + 23*n2 - 7*d*n2 + 6*n2^2, (-5 + d)*(-4 + d - n2))/(d10*d2^n2*d4^2*d5*d6*d7*d8^2*d9) + 
        rat(-17 + 5*d + 6*n2 - 4*d*n2 + 6*n2^2, 3*(-5 + d))/(d2^n2*d3*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-11 + 3*d + 10*n2 - 4*d*n2 + 6*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d5*d6^2*d7*d9) + 
        (d2^(1 - n2)*rat(5 + d - 18*n2 + 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d5^2*d6*d7*d8) + 
        (d2^(1 - n2)*rat(5 + d - 18*n2 + 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d3*d4^2*d5*d6*d7*d8^2*d9) + 
        (d2^(2 - n2)*rat(5 + d - 18*n2 + 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-42 + 20*d - 2*d^2 - 16*n2 + 8*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^4*d5*d6*d8*d9) + 
        (d2^(1 - n2)*rat(-42 + 20*d - 2*d^2 - 16*n2 + 8*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10^2*d3*d4^3*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(-42 + 20*d - 2*d^2 - 16*n2 + 8*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d4^3*d5^2*d6*d7*d8*d9) + 
        rat(-132 + 38*d - 2*d^2 + 74*n2 - 18*d*n2 + 8*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^3*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(3 - 4*d + d^2 + 14*n2 - 6*d*n2 + 8*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10^2*d4^2*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(78 - 34*d + 4*d^2 + 14*n2 - 6*d*n2 + 8*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d4^3*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(18 - 7*d + d^2 - n2 - 3*d*n2 + 8*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d3^2*d4^2*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-5 + 3*d - 8*n2 - 2*d*n2 + 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3^2*d4^2*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(-5 + 3*d - 8*n2 - 2*d*n2 + 8*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^2*d6*d7^2*d8*d9) + 
        rat(-90 + 25*d - d^2 + 55*n2 - 15*d*n2 + 10*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d6^2*d7*d8*d9) + 
        rat(15 - 11*d + 2*d^2 + 25*n2 - 9*d*n2 + 10*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d5*d6^2*d7*d8) + 
        (d2^(1 - n2)*rat(-34 + 10*d + 12*n2 - 8*d*n2 + 12*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d4^3*d5*d6*d7*d8*d9) + 
        rat(18 - 10*d + 2*d^2 + 8*n2 - 8*d*n2 + 14*n2^2, 3*(-5 + d))/(d10*d2^n2*d3*d4*d5*d6*d7*d8*d9) + 
        rat(84 - 34*d + 4*d^2 + 2*n2 - 6*d*n2 + 14*n2^2, 3*(-5 + d)*(-4 + d - n2))/(d10*d2^n2*d3*d4^2*d5*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-174 + 58*d - 4*d^2 + 58*n2 - 18*d*n2 + 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-174 + 58*d - 4*d^2 + 58*n2 - 18*d*n2 + 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(6 - 8*d + 2*d^2 + 28*n2 - 12*d*n2 + 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^3*d5*d6*d8*d9^2) + 
        (d2^(1 - n2)*rat(6 - 8*d + 2*d^2 + 28*n2 - 12*d*n2 + 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d6*d7^2*d8*d9^2) + 
        (d2^(1 - n2)*rat(6 - 8*d + 2*d^2 + 28*n2 - 12*d*n2 + 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d6^3*d7*d8*d9) + 
        (d2^(1 - n2)*rat(66 - 30*d + 4*d^2 + 18*n2 - 10*d*n2 + 16*n2^2, (-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d3^2*d4^3*d5*d6*d8*d9) + 
        (d2^(1 - n2)*rat(36 - 14*d + 2*d^2 - 2*n2 - 6*d*n2 + 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d3*d4^3*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(36 - 14*d + 2*d^2 - 2*n2 - 6*d*n2 + 16*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d3^2*d4^3*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-10 + 6*d - 16*n2 - 4*d*n2 + 16*n2^2, 3*(-5 + d)*(-1 + n2)))/(d10*d3*d4^3*d5*d7*d8*d9) + 
        (d2^(1 - n2)*rat(30 - 22*d + 4*d^2 + 50*n2 - 18*d*n2 + 20*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10^2*d4^2*d5*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(60 - 28*d + 4*d^2 + 20*n2 - 12*d*n2 + 20*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10^2*d3*d4^2*d5*d6*d7^2*d9) + 
        (d2^(1 - n2)*rat(-18 - 10*d + 4*d^2 + 86*n2 - 30*d*n2 + 32*n2^2, 3*(-5 + d)*(-4 + d - n2)*(-1 + n2)))/(d10*d4^3*d5*d6*d7*d8^2*d9);

* n1 == 0 && n10 == 1 && n3 == 2 && n4 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?pos_/d3^2/d4/d5/d6/d7/d8/d9/d10 =

        (5*d2^(-1 - n2))/(d10*d3^2*d4*d5*d6*d7*d8^2) + 
        (2*d2^(-1 - n2))/(d10*d3*d4*d5*d6*d7^2*d8) - 
        d2^(-1 - n2)/(d10*d3*d4*d5^2*d6*d7*d8) + 
        (2*d2^(-1 - n2))/(d10*d3^3*d4*d5*d6*d7*d8) - 
        d2^(-1 - n2)/(d10^2*d3^2*d4*d5*d6*d7*d8) + 
        (3*d2^(-1 - n2))/(d10*d3^2*d4*d5*d6*d7*d9^2) + 
        (2*d2^(-1 - n2))/(d10*d3^2*d4*d5*d6*d8*d9^2) - 
        (4*d2^(-1 - n2))/(d10*d3^2*d4*d5*d7*d8*d9^2) - 
        (2*d2^(-1 - n2))/(d10*d3^2*d4*d6*d7*d8*d9^2) + 
        (4*d2^(-1 - n2))/(d10*d3^2*d5*d6*d7*d8*d9^2) + 
        d2^(-1 - n2)/(d3^2*d4*d5*d6*d7*d8*d9^2) + 
        d2^(-1 - n2)/(d3*d4*d5*d6*d7*d8*d9^2) - 
        (2*d2^(-1 - n2))/(d10*d3^2*d4*d5*d6^2*d7*d9) - 
        (2*d2^(-1 - n2))/(d10*d3*d4*d5*d6^2*d7*d9) + 
        d2^(-1 - n2)/(d10*d3*d4*d5^2*d6*d7*d9) + 
        d2^(-1 - n2)/(d10*d3^2*d4^2*d5*d6*d7*d9) - 
        (2*d2^(-1 - n2))/(d10*d3^3*d4*d5*d6*d7*d9) - 
        (2*d2^(-1 - n2))/(d10*d3^2*d4*d5*d7*d8^2*d9) + 
        (2*d2^(-1 - n2))/(d10*d3^2*d4*d6*d7*d8^2*d9) + 
        d2^(-1 - n2)/(d10*d4*d5*d6*d7*d8^2*d9) - 
        (8*d2^(-1 - n2))/(d3^2*d4*d5*d6*d7*d8^2*d9) - 
        d2^(-1 - n2)/(d3*d4*d5*d6*d7*d8^2*d9) - 
        (3*d2^(-1 - n2))/(d10*d3^2*d4*d5*d6^2*d8*d9) + 
        d2^(-1 - n2)/(d10*d3*d4*d5^2*d6*d8*d9) + 
        (2*d2^(-1 - n2))/(d10*d3*d4^2*d5*d6*d8*d9) + 
        (2*d2^(-1 - n2))/(d10*d3^3*d4*d5*d6*d8*d9) - 
        d2^(-1 - n2)/(d10^2*d3^2*d4*d5*d6*d8*d9) - 
        (2*d2^(-1 - n2))/(d10*d3^2*d4*d6*d7^2*d8*d9) - 
        (2*d2^(-1 - n2))/(d10*d3*d4*d6*d7^2*d8*d9) - 
        d2^(-1 - n2)/(d10*d3*d5*d6*d7^2*d8*d9) + 
        d2^(-1 - n2)/(d10*d4*d5*d6*d7^2*d8*d9) - 
        d2^(-1 - n2)/(d10*d3*d4*d5^2*d7*d8*d9) - 
        (2*d2^(-1 - n2))/(d10*d3*d4^2*d5*d7*d8*d9) - 
        (2*d2^(-1 - n2))/(d10*d3^3*d4*d5*d7*d8*d9) - 
        d2^(-1 - n2)/(d10^2*d3^2*d4*d5*d7*d8*d9) + 
        (2*d2^(-1 - n2))/(d10*d3^2*d4*d6^2*d7*d8*d9) + 
        (2*d2^(-1 - n2))/(d10*d3*d4*d6^2*d7*d8*d9) + 
        d2^(-1 - n2)/(d10*d3*d5*d6^2*d7*d8*d9) - 
        d2^(-1 - n2)/(d10*d3^2*d4^2*d6*d7*d8*d9) + 
        d2^(-1 - n2)/(d10^2*d3^2*d4*d6*d7*d8*d9) + 
        d2^(-1 - n2)/(d10^2*d3^2*d5*d6*d7*d8*d9) - 
        (2*d2^(-1 - n2))/(d10*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(-2 - n2)*(-3 - 3*n2))/(d10*d3^2*d4*d5*d7*d8*d9) + 
        (d2^(-2 - n2)*(-3 - 3*n2))/(d3^2*d4*d5*d6*d7*d8*d9) + 
        (d2^(-2 - n2)*(-2 - 2*n2))/(d3*d4*d5*d6*d7*d8*d9) + 
        7/(2*d10*d2^n2*d3^2*d4*d5*d6^2*d7*d8^2*n2) + 
        7/(2*d10*d2^n2*d3^2*d4*d5^2*d6*d7*d8^2*n2) + 
        2/(d10*d2^n2*d3^2*d4^2*d5*d6*d7^2*d8*n2) + 
        4/(d10*d2^n2*d3^3*d4*d5*d6*d7^2*d8*n2) + 
        6/(d10^2*d2^n2*d3^2*d4*d5*d6*d7^2*d8*n2) + 
        7/(2*d10*d2^n2*d3^2*d4*d5^2*d6^2*d7*d8*n2) - 
        7/(d10^2*d2^n2*d3^2*d4*d5*d6^2*d7*d8*n2) - 
        5/(d10^2*d2^n2*d3^2*d4*d5^2*d6*d7*d8*n2) + 
        4/(d10*d2^n2*d3^3*d4^2*d5*d6*d7*d8*n2) + 
        8/(d10^2*d2^n2*d3^2*d4^2*d5*d6*d7*d8*n2) + 
        6/(d10*d2^n2*d3^2*d4*d5*d6*d8*d9^3*n2) - 
        6/(d10*d2^n2*d3^2*d4*d6*d7*d8*d9^3*n2) + 
        3/(d10*d2^n2*d3^2*d4*d5*d6^2*d7*d9^2*n2) - 
        4/(d10*d2^n2*d3^3*d4*d5*d6*d7*d9^2*n2) - 
        2/(d10*d2^n2*d3^2*d4*d5*d6*d8^2*d9^2*n2) - 
        1/(d10*d2^n2*d3^2*d4*d5*d7*d8^2*d9^2*n2) + 
        3/(d10*d2^n2*d3^2*d4*d6*d7*d8^2*d9^2*n2) - 
        3/(d10*d2^n2*d3^2*d4*d5*d6^2*d8*d9^2*n2) + 
        3/(d10*d2^n2*d3^2*d4^2*d5*d6*d8*d9^2*n2) - 
        2/(d10*d2^n2*d3^3*d4*d5*d6*d8*d9^2*n2) - 
        1/(d10^2*d2^n2*d3^2*d4*d5*d6*d8*d9^2*n2) + 
        3/(d10*d2^n2*d3^2*d4*d6*d7^2*d8*d9^2*n2) - 
        1/(d10^2*d2^n2*d3^2*d4*d5*d7*d8*d9^2*n2) - 
        3/(d10*d2^n2*d3^2*d4^2*d6*d7*d8*d9^2*n2) + 
        2/(d10*d2^n2*d3^3*d4*d6*d7*d8*d9^2*n2) + 
        1/(d10^2*d2^n2*d3^2*d4*d6*d7*d8*d9^2*n2) + 
        1/(d10^2*d2^n2*d3^2*d5*d6*d7*d8*d9^2*n2) + 
        4/(d2^n2*d3^3*d4*d5*d6*d7*d8*d9^2*n2) + 
        9/(2*d10*d2^n2*d3^2*d4^2*d5*d6*d7^2*d9*n2) + 
        9/(d10*d2^n2*d3^3*d4*d5*d6*d7^2*d9*n2) + 
        4/(d10^2*d2^n2*d3^2*d4*d5*d6*d7^2*d9*n2) - 
        6/(d10*d2^n2*d3^2*d4*d5*d6^3*d7*d9*n2) - 
        3/(2*d10*d2^n2*d3^2*d4*d5^2*d6^2*d7*d9*n2) + 
        1/(d10*d2^n2*d3^2*d4^2*d5*d6^2*d7*d9*n2) - 
        5/(d10*d2^n2*d3^3*d4*d5*d6^2*d7*d9*n2) - 
        3/(d10^2*d2^n2*d3^2*d4*d5*d6^2*d7*d9*n2) + 
        5/(d10^2*d2^n2*d3^2*d4*d5^2*d6*d7*d9*n2) + 
        9/(d10*d2^n2*d3^3*d4^2*d5*d6*d7*d9*n2) - 
        2/(d10^2*d2^n2*d3^2*d4^2*d5*d6*d7*d9*n2) + 
        8/(d10^2*d2^n2*d3^3*d4*d5*d6*d7*d9*n2) - 
        7/(2*d10*d2^n2*d3^2*d4*d5*d6^2*d8^2*d9*n2) - 
        7/(2*d10*d2^n2*d3^2*d4*d5^2*d6*d8^2*d9*n2) + 
        3/(2*d10*d2^n2*d3^2*d4^2*d5*d6*d8^2*d9*n2) - 
        4/(d10*d2^n2*d3^3*d4*d5*d6*d8^2*d9*n2) - 
        2/(d10^2*d2^n2*d3^2*d4*d5*d6*d8^2*d9*n2) + 
        9/(2*d10*d2^n2*d3^2*d4*d5^2*d7*d8^2*d9*n2) - 
        6/(d10*d2^n2*d3^2*d4^2*d5*d7*d8^2*d9*n2) - 
        4/(d10*d2^n2*d3^3*d4*d5*d7*d8^2*d9*n2) - 
        2/(d10^2*d2^n2*d3^2*d4*d5*d7*d8^2*d9*n2) - 
        3/(d10*d2^n2*d3^2*d4*d6^2*d7*d8^2*d9*n2) + 
        6/(d10*d2^n2*d3^2*d4^2*d6*d7*d8^2*d9*n2) + 
        4/(d10*d2^n2*d3^3*d4*d6*d7*d8^2*d9*n2) + 
        2/(d10^2*d2^n2*d3^2*d4*d6*d7*d8^2*d9*n2) + 
        4/(d10*d2^n2*d3^3*d5*d6*d7*d8^2*d9*n2) + 
        2/(d10^2*d2^n2*d3^2*d5*d6*d7*d8^2*d9*n2) - 
        3/(2*d2^n2*d3^2*d4^2*d5*d6*d7*d8^2*d9*n2) - 
        8/(d2^n2*d3^3*d4*d5*d6*d7*d8^2*d9*n2) - 
        7/(2*d10*d2^n2*d3^2*d4*d5^2*d6^2*d8*d9*n2) - 
        21/(d10*d2^n2*d3^3*d4*d5*d6^2*d8*d9*n2) + 
        5/(d10^2*d2^n2*d3^2*d4*d5^2*d6*d8*d9*n2) + 
        19/(d10*d2^n2*d3^3*d4^2*d5*d6*d8*d9*n2) + 
        2/(d10^2*d2^n2*d3^2*d4^2*d5*d6*d8*d9*n2) - 
        4/(d10^2*d2^n2*d3^3*d4*d5*d6*d8*d9*n2) - 
        3/(d10*d2^n2*d3^2*d4*d6^2*d7^2*d8*d9*n2) - 
        2/(d10*d2^n2*d3^2*d4^2*d6*d7^2*d8*d9*n2) - 
        4/(d10*d2^n2*d3^3*d4*d6*d7^2*d8*d9*n2) - 
        6/(d10^2*d2^n2*d3^2*d4*d6*d7^2*d8*d9*n2) - 
        19/(d10*d2^n2*d3^3*d5*d6*d7^2*d8*d9*n2) - 
        4/(d10^2*d2^n2*d3^2*d5*d6*d7^2*d8*d9*n2) - 
        9/(2*d2^n2*d3^2*d4^2*d5*d6*d7^2*d8*d9*n2) - 
        9/(d2^n2*d3^3*d4*d5*d6*d7^2*d8*d9*n2) - 
        5/(d10^2*d2^n2*d3^2*d4*d5^2*d7*d8*d9*n2) - 
        8/(d10^2*d2^n2*d3^2*d4^2*d5*d7*d8*d9*n2) - 
        4/(d10^2*d2^n2*d3^3*d4*d5*d7*d8*d9*n2) + 
        6/(d10*d2^n2*d3^2*d4*d6^3*d7*d8*d9*n2) - 
        1/(d10*d2^n2*d3^2*d4^2*d6^2*d7*d8*d9*n2) - 
        2/(d10*d2^n2*d3^3*d4*d6^2*d7*d8*d9*n2) + 
        3/(d10^2*d2^n2*d3^2*d4*d6^2*d7*d8*d9*n2) + 
        21/(d10*d2^n2*d3^3*d5*d6^2*d7*d8*d9*n2) + 
        7/(d10^2*d2^n2*d3^2*d5*d6^2*d7*d8*d9*n2) + 
        7/(d2^n2*d3^3*d4*d5*d6^2*d7*d8*d9*n2) - 
        4/(d10*d2^n2*d3^3*d4^2*d6*d7*d8*d9*n2) + 
        4/(d10^2*d2^n2*d3^3*d4*d6*d7*d8*d9*n2) + 
        4/(d10^2*d2^n2*d3^3*d5*d6*d7*d8*d9*n2) - 
        9/(d2^n2*d3^3*d4^2*d5*d6*d7*d8*d9*n2) - 
        (2*d2^(-1 - n2)*n2)/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-2 - n2)*(2 + 2*n2))/(d10*d3*d4*d5*d6*d7*d8) + 
        (d2^(-2 - n2)*(3 + 3*n2))/(d10*d3^2*d4*d5*d6*d7*d8) + 
        (d2^(-2 - n2)*(3 + 3*n2))/(d10*d3^2*d5*d6*d7*d8*d9) + 
        rat(-2, -4 + d - n2)/(d10^2*d2^n2*d3^2*d4*d6*d7*d8*d9) + 
        rat(-2, -4 + d - n2)/(d10^2*d2^n2*d3^2*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(1, -4 + d - n2))/(d10*d3^2*d4*d5*d6^2*d8*d9) + 
        rat(1, -4 + d - n2)/(d10*d2^n2*d3^2*d4*d5*d6^2*d8*d9) + 
        (d1*rat(1, -4 + d - n2))/(d10*d2^n2*d3^2*d4*d5^2*d6*d8*d9) + 
        (d1*rat(1, -4 + d - n2))/(d10*d2^n2*d3^2*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(1, -4 + d - n2))/(d3^2*d4*d5*d6^2*d7*d8*d9) + 
        rat(1, -4 + d - n2)/(d2^n2*d3^2*d4*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(1, 4 - d + n2))/(d10*d3^2*d4*d5*d6^2*d7*d8) + 
        rat(1, 4 - d + n2)/(d10*d2^n2*d3^2*d4*d5*d6^2*d7*d8) + 
        (d1*rat(1, 4 - d + n2))/(d10*d2^n2*d3^2*d4*d5^2*d6*d7*d8) + 
        (d1*rat(1, 4 - d + n2))/(d2^n2*d3^2*d4*d5*d6*d7*d8^2*d9) + 
        rat(1, 4 - d + n2)/(d10*d2^n2*d3^2*d5^2*d6*d7*d8*d9) + 
        rat(2, -4 + d - n2)/(d10^2*d2^n2*d3^2*d4*d5*d6*d8*d9) + 
        rat(2, -4 + d - n2)/(d10^2*d2^n2*d3^2*d4*d5*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(8 - 2*d, -4 + d - n2))/(d10*d3^2*d4*d5*d6*d7*d9) + 
        rat(-12 + 3*d - 10*n2, 12*n2)/(d2^n2*d3*d4*d5*d6*d7*d8*d9^2) + 
        rat(-12 + 3*d - 10*n2, 12*n2)/(d10*d2^n2*d3*d4*d5^2*d6*d7*d9) + 
        rat(-12 + 3*d - 10*n2, 12*n2)/(d10*d2^n2*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(-36 + 9*d - 10*n2, 12))/(d10*d3*d5*d6*d7*d8*d9) + 
        rat(-36 + 9*d - 10*n2, 12*n2)/(d10*d2^n2*d3*d4*d5*d6*d8*d9^2) + 
        rat(-36 + 9*d - 10*n2, 12*n2)/(d10^2*d2^n2*d4*d5*d6*d7*d8*d9) + 
        rat(-60 + 15*d - 10*n2, 12*n2)/(d10*d2^n2*d3*d4*d5^2*d6*d7*d8) + 
        rat(-60 + 15*d - 10*n2, 12*n2)/(d2^n2*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-60 + 15*d - 10*n2, 12*n2))/(d10*d3*d4*d5*d6^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(46 - 13*d - 6*n2, 2*(-4 + d)))/(d10*d3^2*d4^2*d5*d7*d8*d9) + 
        rat(6 - 3*d - 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5^2*d6*d8*d9^2) + 
        (d2^(-1 - n2)*rat(-26 + 5*d - 6*n2, 2*(-4 + d)))/(d3^2*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(-42 + 9*d - 6*n2, 2*(-4 + d)))/(d3^2*d4*d5*d6^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-66 + 15*d - 6*n2, 2*(-4 + d)))/(d3^2*d4^2*d5*d6*d7*d8*d9) + 
        rat(-66 + 15*d - 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4^2*d5*d7*d8*d9^2) + 
        rat(-66 + 15*d - 6*n2, 2*(-4 + d)*n2)/(d2^n2*d3^2*d4^2*d5*d6^2*d7*d8*d9) + 
        rat(-82 + 19*d - 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5*d7^2*d8*d9^2) + 
        rat(-82 + 19*d - 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5^2*d7*d8*d9^2) + 
        rat(-82 + 19*d - 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5^2*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(-114 + 27*d - 6*n2, 2*(-4 + d)))/(d10*d3^2*d4*d5*d7^2*d8*d9) + 
        rat(-15 + 3*d - 3*n2, (-4 + d)*n2)/(d10*d2^n2*d3^2*d4^2*d5^2*d7*d8*d9) + 
        rat(-31 + 7*d - 3*n2, (-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5*d6^2*d7^2*d9) + 
        rat(-39 + 9*d - 3*n2, (-4 + d)*n2)/(d10*d2^n2*d3^2*d4^2*d5^2*d6*d8*d9) + 
        (d2^(-1 - n2)*rat(12 - 3*d - 2*n2, 12))/(d10*d3*d4*d5*d6*d7*d8) + 
        rat(12 - 3*d - 2*n2, 12*n2)/(d10*d2^n2*d3*d4^2*d5*d6*d8*d9) + 
        rat(12 - 3*d - 2*n2, 12*n2)/(d10*d2^n2*d3*d4*d6^2*d7*d8*d9) + 
        rat(-7 + 2*d - 2*n2, -4 + d - n2)/(d10^2*d2^n2*d3*d4*d5*d6*d7*d8) + 
        rat(-12 + 3*d - 2*n2, 2*n2)/(d2^n2*d3^3*d4*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-12 + 3*d - 2*n2, 3*n2))/(d10*d3^2*d4*d5^2*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-12 + 3*d - 2*n2, 3*n2))/(d10*d3^3*d4*d5*d6*d7*d9) + 
        (d2^(2 - n2)*rat(-12 + 3*d - 2*n2, 3*n2))/(d10*d3^2*d4*d5*d6^2*d7*d8*d9) + 
        rat(-12 + 3*d - 2*n2, 4*n2)/(d10*d2^n2*d3^2*d4*d5*d6*d7^2*d8) + 
        rat(-12 + 3*d - 2*n2, 4*n2)/(d10*d2^n2*d3^2*d4*d5^2*d6*d7*d8) + 
        rat(-12 + 3*d - 2*n2, 4*n2)/(d10^2*d2^n2*d3^2*d4*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-12 + 3*d - 2*n2, 4*n2))/(d10*d3^2*d4*d5*d6*d8*d9^2) + 
        rat(-12 + 3*d - 2*n2, 4*n2)/(d2^n2*d3^2*d4*d5*d6*d7*d8*d9^2) + 
        rat(-12 + 3*d - 2*n2, 4*n2)/(d10*d2^n2*d3^2*d4*d5^2*d6*d7*d9) + 
        (d2^(1 - n2)*rat(-12 + 3*d - 2*n2, 6*n2))/(d3^2*d4*d5*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-12 + 3*d - 2*n2, 6*n2))/(d10*d3^2*d4*d5^2*d6*d7*d9) + 
        (d2^(1 - n2)*rat(-12 + 3*d - 2*n2, 6*n2))/(d10*d3^2*d4*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-12 + 3*d - 2*n2, 6*n2))/(d10*d3*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-12 + 3*d - 2*n2, 6*n2))/(d10*d3^2*d4^2*d5*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-12 + 3*d - 2*n2, 12*n2))/(d10*d3^2*d4^2*d5*d6*d8*d9) + 
        (d2^(1 - n2)*rat(3 - d + n2, -4 + d - n2))/(d10*d3^2*d4*d5*d6*d7*d8*d9) + 
        rat(12 - 3*d + 2*n2, 2*n2)/(d10*d2^n2*d3^3*d4*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(12 - 3*d + 2*n2, 3*n2))/(d10*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(12 - 3*d + 2*n2, 3*n2))/(d10*d3^2*d4*d5^2*d6*d8*d9) + 
        (d2^(1 - n2)*rat(12 - 3*d + 2*n2, 3*n2))/(d3^3*d4*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(12 - 3*d + 2*n2, 4*n2))/(d10*d3^2*d4*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(12 - 3*d + 2*n2, 4*n2))/(d10^2*d3^2*d4*d5*d6*d7*d9) + 
        rat(12 - 3*d + 2*n2, 4*n2)/(d10*d2^n2*d3^2*d4*d5^2*d6*d8*d9) + 
        rat(12 - 3*d + 2*n2, 4*n2)/(d10*d2^n2*d3^2*d4*d6*d7^2*d8*d9) + 
        rat(12 - 3*d + 2*n2, 4*n2)/(d10*d2^n2*d3^2*d4*d5^2*d7*d8*d9) + 
        rat(12 - 3*d + 2*n2, 4*n2)/(d10*d2^n2*d3^2*d4^2*d5*d7*d8*d9) + 
        (d2^(1 - n2)*rat(12 - 3*d + 2*n2, 6*n2))/(d10*d3^2*d4*d5*d6*d7^2*d8) + 
        (d2^(2 - n2)*rat(12 - 3*d + 2*n2, 6*n2))/(d10*d3^2*d4*d5*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(12 - 3*d + 2*n2, 6*n2))/(d10*d3^2*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(12 - 3*d + 2*n2, 6*n2))/(d10*d3^2*d4*d5^2*d7*d8*d9) + 
        (d2^(2 - n2)*rat(12 - 3*d + 2*n2, 6*n2))/(d10*d3^2*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(12 - 3*d + 2*n2, 12*n2))/(d10*d3^2*d4*d5*d6^2*d7*d9) + 
        (d2^(-1 - n2)*rat(-12 + 3*d + 2*n2, 12))/(d3*d4*d5*d6*d7*d8*d9) + 
        rat(-12 + 3*d + 2*n2, 12*n2)/(d10*d2^n2*d3*d4*d5*d6^2*d7*d9) + 
        rat(-12 + 3*d + 2*n2, 12*n2)/(d10*d2^n2*d4^2*d5*d6*d7*d8*d9) + 
        rat(39 - 9*d + 3*n2, (-4 + d)*n2)/(d10*d2^n2*d3^2*d4^2*d5^2*d6*d7*d8) + 
        rat(19 - 4*d + 3*n2, (-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5*d6^2*d7^2*d8) + 
        rat(15 - 3*d + 3*n2, (-4 + d)*n2)/(d10*d2^n2*d3^2*d4^2*d5^2*d6*d7*d9) + 
        (d2^(-1 - n2)*rat(90 - 21*d + 6*n2, 2*(-4 + d)))/(d10*d3^2*d5*d6*d7^2*d8*d9) + 
        rat(82 - 19*d + 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5*d6*d7^2*d9^2) + 
        rat(82 - 19*d + 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5^2*d6*d7*d9^2) + 
        rat(82 - 19*d + 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5^2*d6*d7^2*d9) + 
        (d2^(-1 - n2)*rat(66 - 15*d + 6*n2, 2*(-4 + d)))/(d10*d3^2*d4^2*d5*d6*d7*d8) + 
        rat(66 - 15*d + 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4^2*d5*d6^2*d7*d8) + 
        rat(66 - 15*d + 6*n2, 2*(-4 + d)*n2)/(d2^n2*d3^2*d4^2*d5*d6*d7*d8*d9^2) + 
        (d2^(-1 - n2)*rat(42 - 9*d + 6*n2, 2*(-4 + d)))/(d10*d3^2*d4*d5*d6^2*d7*d8) + 
        rat(12 - 3*d + 6*n2, 4*n2)/(d10*d2^n2*d3^2*d4*d6^2*d7*d8*d9) + 
        rat(18 - 3*d + 6*n2, 2*(-4 + d)*n2)/(d10*d2^n2*d3^2*d4*d5^2*d6*d7^2*d8) + 
        (d2^(-1 - n2)*rat(10 - d + 6*n2, 2*(-4 + d)))/(d10*d3^2*d4*d5*d6*d7^2*d8) + 
        (d2^(-1 - n2)*rat(-46 + 13*d + 6*n2, 2*(-4 + d)))/(d10*d3^2*d5*d6^2*d7*d8*d9) + 
        rat(60 - 15*d + 10*n2, 12*n2)/(d10*d2^n2*d4*d5*d6*d7*d8^2*d9) + 
        rat(60 - 15*d + 10*n2, 12*n2)/(d10*d2^n2*d3*d4*d5^2*d6*d8*d9) + 
        (d2^(-1 - n2)*rat(36 - 9*d + 10*n2, 12))/(d10*d3*d4*d5*d7*d8*d9) + 
        rat(36 - 9*d + 10*n2, 12*n2)/(d10*d2^n2*d3*d4*d6*d7*d8*d9^2) + 
        rat(36 - 9*d + 10*n2, 12*n2)/(d10^2*d2^n2*d3*d4*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(12 - 3*d + 10*n2, 12*n2))/(d10*d3*d4*d5*d6*d7*d8*d9^2) + 
        rat(12 - 3*d + 10*n2, 12*n2)/(d10*d2^n2*d3*d5*d6*d7^2*d8*d9) + 
        rat(12 - 3*d + 10*n2, 12*n2)/(d10*d2^n2*d3*d4*d5^2*d7*d8*d9) + 
        rat(-240 + 120*d - 15*d^2 - 112*n2 + 25*d*n2 - 10*n2^2, 12*(-4 + d - n2)*n2)/(d10*d2^n2*d3*d5*d6^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-48 + 24*d - 3*d^2 - 44*n2 + 9*d*n2 - 6*n2^2, 4*(-4 + d - n2)))/(d10*d3^2*d4*d5*d6*d7*d8) + 
        (d2^(-1 - n2)*rat(-48 + 24*d - 3*d^2 - 40*n2 + 9*d*n2 - 6*n2^2, 4*(-4 + d - n2)))/(d10*d3^2*d5*d6*d7*d8*d9) + 
        rat(-48 + 24*d - 3*d^2 - 25*n2 + 7*d*n2 - 4*n2^2, 3*(-4 + d - n2)*n2)/(d10*d2^n2*d3*d4*d5*d6*d7^2*d8) + 
        rat(48 - 24*d + 3*d^2 + 7*n2 - d*n2 - 2*n2^2, 3*(-4 + d - n2)*n2)/(d10*d2^n2*d3^2*d4*d5*d6*d7*d9) + 
        rat(48 - 24*d + 3*d^2 + 7*n2 - d*n2 - 2*n2^2, 3*(-4 + d - n2)*n2)/(d10*d2^n2*d3*d4^2*d5*d7*d8*d9) + 
        rat(12 - 3*d - 10*n2 + 3*d*n2 - 2*n2^2, 4*n2)/(d10*d2^n2*d3^2*d5*d6*d7*d8*d9) + 
        rat(-48 + 24*d - 3*d^2 - 24*n2 + 5*d*n2 - 2*n2^2, 2*(-4 + d - n2)*n2)/(d10*d2^n2*d3^3*d5*d6*d7*d8*d9) + 
        rat(-48 + 24*d - 3*d^2 - 24*n2 + 5*d*n2 - 2*n2^2, 4*(-4 + d - n2)*n2)/(d10*d2^n2*d3^2*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-48 + 24*d - 3*d^2 - 23*n2 + 5*d*n2 - 2*n2^2, 3*(-4 + d - n2)*n2))/(d10*d3^2*d5*d6^2*d7*d8*d9) + 
        rat(-48 + 24*d - 3*d^2 - 16*n2 + 5*d*n2 - 2*n2^2, 4*(-4 + d - n2)*n2)/(d10*d2^n2*d3^2*d4*d5*d7*d8^2*d9) + 
        rat(-48 + 24*d - 3*d^2 - 16*n2 + 5*d*n2 - 2*n2^2, 4*(-4 + d - n2)*n2)/(d10*d2^n2*d3^2*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(48 - 24*d + 3*d^2 + 8*n2 - 5*d*n2 + 2*n2^2, 12*(-4 + d - n2)*n2))/(d10*d3^2*d4*d6^2*d7*d8*d9) + 
        rat(48 - 24*d + 3*d^2 + 16*n2 - 5*d*n2 + 2*n2^2, 4*(-4 + d - n2)*n2)/(d10*d2^n2*d3^2*d4*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(48 - 24*d + 3*d^2 + 23*n2 - 5*d*n2 + 2*n2^2, 3*(-4 + d - n2)*n2))/(d3^2*d4*d5*d6*d7*d8^2*d9) + 
        rat(48 - 24*d + 3*d^2 + 24*n2 - 5*d*n2 + 2*n2^2, 2*(-4 + d - n2)*n2)/(d10*d2^n2*d3^3*d4*d5*d6*d8*d9) + 
        rat(48 - 24*d + 3*d^2 + 24*n2 - 5*d*n2 + 2*n2^2, 4*(-4 + d - n2)*n2)/(d2^n2*d3^2*d4*d5*d6*d7*d8^2*d9) + 
        rat(-12 + 3*d + 10*n2 - 3*d*n2 + 2*n2^2, 4*n2)/(d10*d2^n2*d3^2*d4*d5*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-80 + 40*d - 5*d^2 - 16*n2 + 3*d*n2 + 2*n2^2, 4*(-4 + d - n2)*n2))/(d10*d3*d4^2*d5*d6*d7*d8*d9) + 
        rat(48 - 24*d + 3*d^2 + 25*n2 - 7*d*n2 + 4*n2^2, 3*(-4 + d - n2)*n2)/(d10*d2^n2*d3*d4*d6*d7^2*d8*d9) + 
        rat(48 - 24*d + 3*d^2 + 32*n2 - 9*d*n2 + 6*n2^2, 4*(-4 + d - n2)*n2)/(d10*d2^n2*d3^2*d4*d6*d7*d8*d9^2) + 
        rat(48 - 24*d + 3*d^2 + 32*n2 - 9*d*n2 + 6*n2^2, 4*(-4 + d - n2)*n2)/(d10*d2^n2*d3^2*d4*d5*d6^2*d7*d9) + 
        rat(48 - 24*d + 3*d^2 + 32*n2 - 9*d*n2 + 6*n2^2, 4*(-4 + d - n2)*n2)/(d10^2*d2^n2*d3^2*d4*d5*d6*d7*d9) + 
        rat(48 - 24*d + 3*d^2 + 32*n2 - 9*d*n2 + 6*n2^2, 4*n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4*d5*d6*d8*d9^2) + 
        (d2^(-1 - n2)*rat(16 - 16*d + 3*d^2 + 36*n2 - 9*d*n2 + 6*n2^2, 4*(-4 + d - n2)))/(d3^2*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(48 - 24*d + 3*d^2 + 40*n2 - 9*d*n2 + 6*n2^2, 4*(-4 + d - n2)))/(d10*d3^2*d4*d5*d7*d8*d9) + 
        rat(48 - 24*d + 3*d^2 + 40*n2 - 9*d*n2 + 6*n2^2, 4*n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4^2*d5*d6*d8*d9) + 
        rat(15 - 11*d + 2*d^2 + 24*n2 - 8*d*n2 + 6*n2^2, 3*(-4 + d - n2))/(d10*d2^n2*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(48 - 24*d + 3*d^2 + 48*n2 - 13*d*n2 + 10*n2^2, 4*(-4 + d - n2)*n2))/(d10^2*d3*d4*d5*d6*d7*d8*d9) + 
        rat(-144 + 72*d - 9*d^2 - 68*n2 + 23*d*n2 - 3*d^2*n2 + 2*n2^2 + 5*d*n2^2 - 2*n2^3, 12*(-4 + d - n2)*n2)/(d2^n2*d3^2*d4*d5*d6*d7*d8*d9) + 
        rat(-48 + 24*d - 3*d^2 - 20*n2 + 5*d*n2 + 3*d^2*n2 - 18*n2^2 - 5*d*n2^2 + 2*n2^3, 12*(-4 + d - n2)*n2)/(d10*d2^n2*d3^2*d4*d5*d6*d7*d8);
        
* n1 == 0 && n10 == 1 && n3 == 1 && n4 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n2 != 1
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?{>1}/d3/d4/d5/d6/d7/d8/d9/d10 =

        (d2^(1 - n2)*rat(-16, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d7*d8^3*d9) + 
        (d2^(1 - n2)*rat(-16, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d5*d6*d7^2*d8*d9) + 
        rat(-10, 9*(3 - d + n2))/(d10*d2^n2*d3^2*d4*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-10, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d5*d6*d7^2*d9) + 
        (d2^(1 - n2)*rat(-10, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6*d7^2*d9) + 
        (d2^(1 - n2)*rat(-10, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d5^2*d6*d7*d9) + 
        (d2^(1 - n2)*rat(-10, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4^2*d5*d6*d7*d9) + 
        rat(-8, 9*(3 - d + n2))/(d10*d2^n2*d3^2*d4*d5*d7*d8*d9) + 
        rat(-8, 9*(3 - d + n2))/(d10*d2^n2*d3*d5^2*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d6*d7^2*d8^2) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6*d7*d8^2) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5^2*d6*d7^2*d8) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d3*d4*d5*d6*d7*d8^2*d9^2) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d5*d6*d7^2*d8*d9^2) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d5^2*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d5*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d5^2*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^3*d5*d6*d7*d8*d9) + 
        rat(-4, 9*(3 - d + n2))/(d10*d2^n2*d3*d4*d5*d6*d7*d8^2) + 
        rat(-4, 9*(3 - d + n2))/(d10*d2^n2*d3*d4^2*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-4, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d5^2*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-4, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d6^2*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-4, 3*(-1 + n2)*(3 - d + n2)))/(d3*d4^2*d5*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-4, 3*(-1 + n2)*(3 - d + n2)))/(d10^2*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d6^2*d7*d8^2) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5^2*d6*d7*d8^2) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5^2*d6^2*d7*d8) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5^2*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10^2*d3^2*d4*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d6^2*d7^2*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5^2*d6^2*d7*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6^2*d7*d9) + 
        (d2^(2 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6^2*d7*d9) + 
        (d2^(2 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10^2*d3^2*d4*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d3*d4*d5*d6^2*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d3*d4*d5^2*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d6^2*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d5*d6^2*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d3*d4*d5^2*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4^3*d5*d6*d7*d8*d9) + 
        rat(-2, 9*(3 - d + n2))/(d2^n2*d3*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-2, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d5*d6*d7^2*d8) + 
        (d2^(1 - n2)*rat(-2, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6*d7^2*d8) + 
        (d2^(1 - n2)*rat(-2, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4^2*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-2, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d5*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-2, 9*(-1 + n2)*(3 - d + n2)))/(d10^2*d4^2*d5*d6*d7*d8*d9) + 
        rat(2, 3*(3 - d + n2))/(d2^n2*d3*d4*d5*d6*d7^2*d8*d9) + 
        rat(2, 9*(3 - d + n2))/(d10*d2^n2*d3*d4^2*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(2, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(2, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(2, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(2, 3*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4^2*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d5*d6*d8*d9^2) + 
        (d2^(1 - n2)*rat(2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d5*d6^2*d7*d9) + 
        (d2^(1 - n2)*rat(2, 9*(-1 + n2)*(3 - d + n2)))/(d10^2*d3*d4^2*d5*d6*d7*d9) + 
        rat(4, 3*(3 - d + n2))/(d10*d2^n2*d3*d4*d5*d7^2*d8*d9) + 
        rat(4, 9*(3 - d + n2))/(d10*d2^n2*d3*d4^2*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(4, 3*(-1 + n2)*(3 - d + n2)))/(d10^2*d3*d4*d5*d6*d7^2*d9) + 
        (d2^(1 - n2)*rat(4, 3*(-1 + n2)*(3 - d + n2)))/(d10*d4^2*d5*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d6^2*d8^2*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5^2*d6*d8^2*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4*d5*d6^2*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4*d5^2*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5^2*d6^2*d8*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6^2*d8*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5^2*d6*d8*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^3*d5*d6*d8*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4*d5*d6^2*d7^2*d8*d9) + 
        (d2^(2 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4*d5^2*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(4, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4^2*d5^2*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d6^2*d7^2*d8) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4*d5*d6*d7*d8^2*d9^2) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6*d8*d9^2) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d6*d7^2*d8*d9^2) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4*d5*d6*d7^2*d8*d9^2) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4*d5^2*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6*d8^2*d9) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d6*d7^2*d8^2*d9) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d5^2*d6*d8*d9) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^3*d4*d5*d6*d8*d9) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4*d5^2*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(8, 9*(-1 + n2)*(3 - d + n2)))/(d3^2*d4*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(10, 9*(-1 + n2)*(3 - d + n2)))/(d3*d4^2*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(10, 9*(-1 + n2)*(3 - d + n2)))/(d3^2*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(10, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4^2*d5^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(10, 9*(-1 + n2)*(3 - d + n2)))/(d3^2*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(16, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d6*d7*d8^3*d9) + 
        (d2^(1 - n2)*rat(16, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5^2*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(16, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4^2*d5*d6*d8*d9) + 
        (d2^(1 - n2)*rat(16 - 4*d, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d5*d7*d8*d9) + 
        rat(8 - 2*d, 3*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d5*d6^2*d7*d9) + 
        rat(-8 + 2*d, 3*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d5^2*d6*d7*d9) + 
        rat(-8 + 2*d, 3*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d4*d5*d6*d7^2*d8*d9) + 
        rat(-8 + 2*d, 3*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d6^2*d7*d8*d9) + 
        (d2^(2 - n2)*rat(-16 + 4*d, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-16 + 4*d, 9*(-3 + d - n2)*(-1 + n2)))/(d3^2*d4*d5*d6*d7*d8*d9) + 
        rat(32 - 6*d - 8*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d5^2*d7*d8*d9) + 
        rat(35 - 7*d - 7*n2, 9*(-3 + d - n2)*(-1 + n2))/(d2^n2*d3^2*d4*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-6 + 4*d - 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-6 + 4*d - 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10^2*d3*d4*d5*d6*d7*d9) + 
        rat(-6 + 4*d - 6*n2, 9*(3 - d + n2))/(d10*d2^n2*d3*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-30 + 10*d - 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d4*d5*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(-30 + 10*d - 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4*d5^2*d6*d8*d9) + 
        (d2^(1 - n2)*rat(-30 + 10*d - 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d5*d6^2*d7*d8*d9) + 
        (d2^(2 - n2)*rat(12 - 4*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10^2*d3*d4*d5*d6*d7*d8*d9) + 
        rat(28 - 6*d - 4*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d6*d7*d8^2*d9) + 
        rat(28 - 6*d - 4*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3^2*d4*d5*d6*d8*d9) + 
        (d2^(1 - n2)*rat(-8 + 4*d - 4*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6*d7*d9) + 
        (d2^(1 - n2)*rat(-8 + 4*d - 4*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d6*d7^2*d8*d9) + 
        (d2^(2 - n2)*rat(-8 + 4*d - 4*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6*d7*d8*d9) + 
        rat(-20 + 6*d - 4*n2, 9*(-3 + d - n2)*(-1 + n2))/(d2^n2*d3*d4*d5*d6*d7*d8*d9^2) + 
        rat(-20 + 6*d - 4*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4^2*d5*d6*d8*d9) + 
        rat(-20 + 6*d - 4*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3^2*d5*d6*d7*d8*d9) + 
        rat(31 - 7*d - 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d6*d7*d8*d9^2) + 
        rat(31 - 7*d - 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10^2*d2^n2*d3*d4*d5*d6*d7*d9) + 
        rat(31 - 7*d - 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d4*d5*d6*d7*d8^2*d9) + 
        rat(31 - 7*d - 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d5^2*d6*d8*d9) + 
        rat(31 - 7*d - 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-21 + 5*d - 3*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4^2*d5*d6*d8*d9) + 
        (d2^(1 - n2)*rat(-21 + 5*d - 3*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4*d6^2*d7*d8*d9) + 
        rat(-21 + 5*d - 3*n2, 9*(3 - d + n2))/(d2^n2*d3*d4*d5*d6*d7*d8*d9) + 
        rat(-25 + 7*d - 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3^2*d4*d5*d6*d7*d9) + 
        (d2^(2 - n2)*rat(-26 + 6*d - 2*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(-26 + 6*d - 2*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(-26 + 6*d - 2*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5^2*d7*d8*d9) + 
        rat(29 - 7*d - n2, 9*(-1 + n2)*(3 - d + n2))/(d10*d2^n2*d3*d4*d6*d7^2*d8*d9) + 
        rat(27 - 7*d + n2, 9*(-1 + n2)*(3 - d + n2))/(d2^n2*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d2^(1 - n2)*rat(26 - 6*d + 2*n2, 9*(-1 + n2)*(3 - d + n2)))/(d3*d4*d5*d6*d7*d8*d9^2) + 
        (d2^(1 - n2)*rat(26 - 6*d + 2*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5^2*d6*d7*d9) + 
        (d2^(1 - n2)*rat(26 - 6*d + 2*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(1 - n2)*rat(21 - 5*d + 3*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4*d5*d6^2*d7*d9) + 
        (d2^(1 - n2)*rat(21 - 5*d + 3*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d4^2*d5*d6*d7*d8*d9) + 
        rat(21 - 5*d + 3*n2, 9*(3 - d + n2))/(d10*d2^n2*d3*d4*d5*d6*d7*d8) + 
        rat(-31 + 7*d + 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d5^2*d6*d7*d8) + 
        rat(-31 + 7*d + 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d5*d6*d8*d9^2) + 
        rat(-31 + 7*d + 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4^2*d5*d7*d8*d9) + 
        rat(-31 + 7*d + 3*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10^2*d2^n2*d4*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(-8 + 4*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3^2*d4*d5*d6*d7*d8) + 
        (d2^(1 - n2)*rat(-8 + 4*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10^2*d3*d4*d5*d6*d7*d8) + 
        rat(20 - 6*d + 4*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10^2*d2^n2*d3*d4*d5*d6*d7*d8) + 
        rat(20 - 6*d + 4*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(1 - n2)*rat(8 - 4*d + 4*n2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d6*d7^2*d8) + 
        rat(4 - 2*d + 4*n2, 3*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d5*d6*d7^2*d8*d9) + 
        rat(-28 + 6*d + 4*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d5*d7*d8^2*d9) + 
        rat(23 - 7*d + 5*n2, 9*(-3 + d - n2)*(-1 + n2))/(d10*d2^n2*d3*d4*d5*d6*d7^2*d8) + 
        (d2^(1 - n2)*rat(30 - 10*d + 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4*d5^2*d6*d7*d8) + 
        (d2^(1 - n2)*rat(30 - 10*d + 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d3*d4*d5*d6*d7*d8^2*d9) + 
        (d2^(2 - n2)*rat(30 - 10*d + 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4*d5*d6^2*d7*d8*d9) + 
        (d2^(1 - n2)*rat(6 - 4*d + 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10*d3*d4*d5*d6*d8*d9^2) + 
        (d2^(1 - n2)*rat(6 - 4*d + 6*n2, 9*(-3 + d - n2)*(-1 + n2)))/(d10^2*d4*d5*d6*d7*d8*d9) + 
        rat(6 - 4*d + 6*n2, 9*(3 - d + n2))/(d10*d2^n2*d3*d4*d5*d7*d8*d9) + 
        (d2^(1 - n2)*rat(19 - 3*d - 2*d^2 - 14*n2 + 11*d*n2 - 9*n2^2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-20*n2 + 6*d*n2 - 4*n2^2, 9*(-1 + n2)*(3 - d + n2)))/(d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(31*n2 - 7*d*n2 - 3*n2^2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-31*n2 + 7*d*n2 + 3*n2^2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(20*n2 - 6*d*n2 + 4*n2^2, 9*(-1 + n2)*(3 - d + n2)))/(d10*d3*d4*d5*d6*d7*d8);

* n10 == 1 && n2 == 1 && n3 == 1 && n4 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n1 != 0
        id,only,ifmatch->sortme 1/d1^n1?neg_/d2/d3/d4/d5/d6/d7/d8/d9/d10 =
        
        rat(-3, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d4*d5^2*d6*d7*d8) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d10*d2*d3*d4*d5*d6^2*d7*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d10*d2^2*d3*d4*d5*d6*d7*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d10^2*d2*d3*d4*d5*d6*d7*d9) + 
        rat(-3, -9 + 2*d - n1)/(d1^n1*d2*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d10*d2*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d10*d2*d3*d4*d5*d7^2*d8*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d2*d3*d4*d5*d6*d7^2*d8*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d10*d2*d3*d4*d5^2*d7*d8*d9) + 
        rat(-3, -9 + 2*d - n1)/(d1^n1*d10*d3*d4*d5*d6^2*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d10*d2*d3*d4*d5*d6^2*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d2*d3^2*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d10*d2^2*d3*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-3, -9 + 2*d - n1))/(d10^2*d2*d3*d4*d5*d6*d7*d8*d9) + 
        rat(-2, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d4*d5*d6*d8*d9^2) + 
        rat(-2, -9 + 2*d - n1)/(d1^n1*d10*d3*d4*d5*d6*d7*d8*d9^2) + 
        rat(-2, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d5*d6*d7^2*d8*d9) + 
        rat(-2, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d4*d5^2*d7*d8*d9) + 
        rat(-2, -9 + 2*d - n1)/(d1^n1*d10*d2^2*d3*d5*d6*d7*d8*d9) + 
        rat(-2, -9 + 2*d - n1)/(d1^n1*d10^2*d2*d4*d5*d6*d7*d8*d9) + 
        rat(1, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d4*d5*d6*d7^2*d8) + 
        rat(1, -9 + 2*d - n1)/(d1^n1*d10*d2^2*d3*d4*d5*d6*d7*d8) + 
        rat(1, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d4^2*d5*d6*d8*d9) + 
        rat(1, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d4*d6^2*d7*d8*d9) + 
        rat(1, -9 + 2*d - n1)/(d1^n1*d10*d3*d4^2*d5*d6*d7*d8*d9) + 
        rat(1, -9 + 2*d - n1)/(d1^n1*d2*d3^2*d4*d5*d6*d7*d8*d9) + 
        rat(1, 9 - 2*d + n1)/(d1^n1*d10*d2*d3*d4*d5*d6^2*d7*d9) + 
        rat(1, 9 - 2*d + n1)/(d1^n1*d10*d2*d3^2*d4*d5*d6*d7*d9) + 
        rat(1, 9 - 2*d + n1)/(d1^n1*d10*d2*d3*d4*d6*d7^2*d8*d9) + 
        rat(1, 9 - 2*d + n1)/(d1^n1*d10*d2*d3*d4^2*d5*d7*d8*d9) + 
        rat(1, 9 - 2*d + n1)/(d1^n1*d10*d2*d4^2*d5*d6*d7*d8*d9) + 
        rat(1, 9 - 2*d + n1)/(d1^n1*d2^2*d3*d4*d5*d6*d7*d8*d9) + 
        rat(2, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d4*d6*d7*d8*d9^2) + 
        rat(2, -9 + 2*d - n1)/(d1^n1*d2*d3*d4*d5*d6*d7*d8*d9^2) + 
        rat(2, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d4*d5^2*d6*d7*d9) + 
        rat(2, -9 + 2*d - n1)/(d1^n1*d10^2*d2*d3*d4*d5*d6*d7*d9) + 
        rat(2, -9 + 2*d - n1)/(d1^n1*d10*d2*d4*d5*d6*d7^2*d8*d9) + 
        rat(2, -9 + 2*d - n1)/(d1^n1*d10*d2^2*d3*d4*d5*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d2*d3*d4*d5*d6*d7^2*d8) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d2*d3^2*d4*d5*d6*d7*d8) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d2*d3*d4*d5*d6*d7^2*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d2*d3*d4*d5^2*d6*d7*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d2*d3^2*d4*d5*d6*d7*d9) + 
        rat(3, -9 + 2*d - n1)/(d1^n1*d10*d2*d4*d5*d6*d7*d8^2*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d3*d4*d5*d6*d7*d8^2*d9) + 
        rat(3, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d4*d5^2*d6*d8*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d2*d3*d5*d6*d7^2*d8*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d3*d4*d5*d6*d7^2*d8*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d2*d3*d4*d5*d6*d7^2*d8*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d2*d3*d4*d6^2*d7*d8*d9) + 
        rat(3, -9 + 2*d - n1)/(d1^n1*d10*d2*d3*d5*d6^2*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10*d2*d3*d5^2*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(3, -9 + 2*d - n1))/(d10^2*d2*d4*d5*d6*d7*d8*d9) + 
        (d1^(-2 - n1)*rat(-3 - 3*n1, 9 - 2*d + n1))/(d10*d2*d3*d4*d5*d6*d7*d9) + 
        (d1^(-1 - n1)*rat(-15 + 3*d - 2*n1, -9 + 2*d - n1))/(d10*d2*d3*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-2*n1, 9 - 2*d + n1))/(d10*d2*d3*d4*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-2*n1, 9 - 2*d + n1))/(d10*d2*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-n1, 9 - 2*d + n1))/(d10*d2*d3*d4*d5*d6*d7*d9) + 
        (d1^(-1 - n1)*rat(n1, 9 - 2*d + n1))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(2*n1, 9 - 2*d + n1))/(d10*d2*d3*d4*d5*d6*d7*d8) + 
        (d1^(-1 - n1)*rat(2*n1, 9 - 2*d + n1))/(d10*d2*d3*d5*d6*d7*d8*d9) + 
        (d1^(-2 - n1)*rat(3 + 3*n1, 9 - 2*d + n1))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d1^(-2 - n1)*rat(3 + 3*n1, 9 - 2*d + n1))/(d10*d2*d3*d4*d5*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n2 == 1 && n3 == 1 && n4 == 2 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1
        if(count(d1,1)==0) id,only 1/d2/d3/d4^2/d5/d6/d7/d8/d9/d10 = 1/(d10*d2*d3^2*d4*d5*d6*d7*d8*d9);
        
        if(count(d1,1)==0) id,only 1/d2/d3/d4/d5/d6/d7/d8/d9/d10 = PR0/intX;                

        #call zeroX

        endif;        
        
        goto endrec;         
        la sortme;
        $irep = 0;
        la endrec;
        
        ModuleOption,minimum,$irep;
        .sort:redX-`$repcount++';
        #redefine irep "`$irep'"
#enddo
#endprocedure
*--#] redX :

* 
* Planar 9 line
* 
*--#[ redH :
#procedure redH
#$repcount = 1;        
#do irep=1,1
        #$irep = 1;                
        if(count(intH,1));
* den: minus
* [1111111110]                  
* n9 (-)                        
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?{>1}/d10^n10?neg0_ =
        -((d1^(-1 - n1)*d7^(1 - n7))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9)) +
        (d1^(-1 - n1)*d5^(1 - n5))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d1^(-1 - n1)*d3^(1 - n3))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d1^(-2 - n1)*d4^(1 - n4)*d9^(1 - n9)*(-2 - 2*n1))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*
        d7^n7*d8^n8*(-1 + n9)) + (d1^(-2 - n1)*d5^(1 - n5)*d9^(1 - n9)*(2 + 2*n1))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d8^(1 - n8)*d9^(1 - n9)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d9^(1 - n9)*n10)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d9^(1 - n9)*n10)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d9^(2 - n9)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d8^n8*(-1 + n9)) - (d1^(-1 - n1)*d2^(-1 - n2)*d8^(1 - n8)*d9^(1 - n9)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) +
        (d1^(-1 - n1)*d2^(-1 - n2)*d6^(1 - n6)*d9^(1 - n9)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n9)) -
        (d1^(-1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*d9^(1 - n9)*n2)/(d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) +
        (d2^(-1 - n2)*d9^(1 - n9)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        (-1 + n9)) + (d1^(-1 - n1)*d3^(-1 - n3)*d7^(1 - n7)*d9^(1 - n9)*n3)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n9)) -
        (d1^(-1 - n1)*d3^(-1 - n3)*d4^(1 - n4)*d9^(1 - n9)*n3)/(d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) +
        (d3^(-1 - n3)*d9^(1 - n9)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        (-1 + n9)) - (d1^(-1 - n1)*d3^(-1 - n3)*d9^(2 - n9)*n3)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) -
        (2*d1^(-1 - n1)*d4^(-1 - n4)*d5^(1 - n5)*d9^(1 - n9)*n4)/(d10^n10*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*(-1 + n9)) +
        (2*d4^(-1 - n4)*d9^(1 - n9)*n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*
        (-1 + n9)) + (2*d1^(-1 - n1)*d4^(1 - n4)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d10^n10*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*(-1 + n9)) -
        (2*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*
        (-1 + n9)) + (d1^(-1 - n1)*d6^(-1 - n6)*d8^(1 - n8)*d9^(1 - n9)*n6)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*(-1 + n9)) -
        (d1^(-1 - n1)*d5^(1 - n5)*d6^(-1 - n6)*d9^(1 - n9)*n6)/(d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n9)) +
        (d1^(-1 - n1)*d4^(1 - n4)*d6^(-1 - n6)*d9^(1 - n9)*n6)/(d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d8^n8*(-1 + n9)) -
        (d1^(-1 - n1)*d2^(1 - n2)*d6^(-1 - n6)*d9^(1 - n9)*n6)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n9)) -
        (d1^(-1 - n1)*d5^(1 - n5)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*(-1 + n9)) +
        (d1^(-1 - n1)*d4^(1 - n4)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d8^n8*(-1 + n9)) -
        (d1^(-1 - n1)*d3^(1 - n3)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n9)) +
        (d1^(-1 - n1)*d7^(-1 - n7)*d9^(2 - n9)*n7)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*
        d8^n8*(-1 + n9)) - (d1^(-1 - n1)*d6^(1 - n6)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*(-1 + n9)) +
        (d1^(-1 - n1)*d5^(1 - n5)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*(-1 + n9)) +
        (d1^(-1 - n1)*d2^(1 - n2)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) -
        (d8^(-1 - n8)*d9^(1 - n9)*n8)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        (-1 + n9));

* n1 (-)
        id,ifmatch->sortme 1/d1^n1?{>1}/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ =
        d5^(1 - n5)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) -
        d4^(1 - n4)/(3*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d1^(1 - n1)*d10^(-1 - n10)*d9^(1 - n9)*n10)/(12*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d8^n8*(-1 + n1)) + (d1^(1 - n1)*d10^(-1 - n10)*d8^(1 - n8)*n10)/(12*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n1)) -
        (d1^(1 - n1)*d10^(-1 - n10)*d7^(1 - n7)*n10)/(2*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*
        d9^n9*(-1 + n1)) - (d1^(1 - n1)*d10^(-1 - n10)*d6^(1 - n6)*n10)/(2*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n1)) +
        (d1^(1 - n1)*d10^(-1 - n10)*d5^(1 - n5)*n10)/(2*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n1)) + (d1^(1 - n1)*d10^(-1 - n10)*d4^(1 - n4)*n10)/(2*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n1)) +
        (5*d1^(1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*n10)/(12*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n1)) + (5*d1^(1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*n10)/(12*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n1)) +
        (2*d1^(1 - n1)*d10^(-1 - n10)*n10)/(3*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n1)) - (d1^(2 - n1)*d10^(-1 - n10)*n10)/(2*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d8^n8*d9^n9*(-1 + n1)) + (d1^(1 - n1)*d2^(-1 - n2)*d8^(1 - n8)*n2)/(4*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n1)) -
        (d1^(1 - n1)*d2^(-1 - n2)*d6^(1 - n6)*n2)/(12*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*
        d9^n9*(-1 + n1)) + (d1^(1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*n2)/(12*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n1)) -
        (d1^(2 - n1)*d2^(-1 - n2)*n2)/(4*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n1)) + (d1^(1 - n1)*d3^(-1 - n3)*d9^(1 - n9)*n3)/(4*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n1)) -
        (d1^(1 - n1)*d3^(-1 - n3)*d7^(1 - n7)*n3)/(12*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*
        d9^n9*(-1 + n1)) + (d1^(1 - n1)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(12*d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n1)) -
        (d1^(2 - n1)*d3^(-1 - n3)*n3)/(4*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n1)) + (2*d1^(1 - n1)*d4^(1 - n4)*d5^(-1 - n5)*n5)/(3*d10^n10*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n1)) -
        (2*d1^(2 - n1)*d5^(-1 - n5)*n5)/(3*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n1)) + (d1^(1 - n1)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(4*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n1)) -
        (d1^(1 - n1)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(4*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*
        d9^n9*(-1 + n1)) + (d1^(1 - n1)*d4^(1 - n4)*d6^(-1 - n6)*n6)/(12*d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n1)) -
        (d1^(1 - n1)*d2^(1 - n2)*d6^(-1 - n6)*n6)/(12*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*
        d9^n9*(-1 + n1)) + (d1^(1 - n1)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(4*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n1)) -
        (d1^(1 - n1)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(4*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*
        d9^n9*(-1 + n1)) + (d1^(1 - n1)*d4^(1 - n4)*d7^(-1 - n7)*n7)/(12*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n1)) -
        (d1^(1 - n1)*d3^(1 - n3)*d7^(-1 - n7)*n7)/(12*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*
        d9^n9*(-1 + n1)) - (5*d1^(1 - n1)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(12*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n1)) +
        (5*d1^(1 - n1)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(12*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*
        d9^n9*(-1 + n1)) + (7*d1^(1 - n1)*d2^(1 - n2)*d8^(-1 - n8)*n8)/(12*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n1)) -
        (7*d1^(2 - n1)*d8^(-1 - n8)*n8)/(12*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*
        (-1 + n1)) - (5*d1^(1 - n1)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(12*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n1)) +
        (5*d1^(1 - n1)*d5^(1 - n5)*d9^(-1 - n9)*n9)/(12*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*
        d8^n8*(-1 + n1)) + (7*d1^(1 - n1)*d3^(1 - n3)*d9^(-1 - n9)*n9)/(12*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n1)) -
        (7*d1^(2 - n1)*d9^(-1 - n9)*n9)/(12*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        (-1 + n1)) + (d1^(1 - n1)*rat(6 + d - 6*n1 + n10, 6*(-1 + n1)))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);
        
*  3: n2      (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?{>1}/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ =

        d8^(1 - n8)/(4*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) +
        d6^(1 - n6)/(4*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9) -
        d4^(1 - n4)/(4*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        d1^(1 - n1)/(4*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d10^(-1 - n10)*d2^(1 - n2)*d9^(1 - n9)*n10)/(2*d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d8^n8*(-1 + n2)) + (d10^(-1 - n10)*d2^(1 - n2)*d8^(1 - n8)*n10)/(4*d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n2)) +
        (d10^(-1 - n10)*d2^(1 - n2)*d7^(1 - n7)*n10)/(4*d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*
        d9^n9*(-1 + n2)) + (d10^(-1 - n10)*d2^(1 - n2)*d6^(1 - n6)*n10)/(4*d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n2)) -
        (d10^(-1 - n10)*d2^(1 - n2)*d5^(1 - n5)*n10)/(4*d1^n1*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n2)) - (d10^(-1 - n10)*d2^(1 - n2)*d4^(1 - n4)*n10)/(4*d1^n1*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n2)) -
        (d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*n10)/(2*d1^n1*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n2)) + (d1^(1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*n10)/(4*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n2)) -
        (d10^(-1 - n10)*d2^(1 - n2)*n10)/(2*d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n2)) - (3*d10^(-1 - n10)*d2^(2 - n2)*n10)/(4*d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*
        d7^n7*d8^n8*d9^n9*(-1 + n2)) + (d2^(1 - n2)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(4*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n2)) -
        (d2^(1 - n2)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(4*d1^n1*d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*
        d9^n9*(-1 + n2)) + (3*d2^(1 - n2)*d4^(1 - n4)*d6^(-1 - n6)*n6)/(4*d1^n1*d10^n10*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n2)) -
        (3*d2^(2 - n2)*d6^(-1 - n6)*n6)/(4*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*
        (-1 + n2)) + (d2^(1 - n2)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(4*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n2)) -
        (d2^(1 - n2)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(4*d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*
        d9^n9*(-1 + n2)) + (3*d1^(1 - n1)*d2^(1 - n2)*d8^(-1 - n8)*n8)/(4*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n2)) -
        (3*d2^(2 - n2)*d8^(-1 - n8)*n8)/(4*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*
        (-1 + n2)) + (d2^(1 - n2)*rat(4 + d - n10 - 4*n2, 4*(-1 + n2)))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);
        
*  4: n3      (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?{>1}/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ =
        
        d9^(1 - n9)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) -
        d1^(1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(1 - n5)*n1)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n3)) - (d1^(-1 - n1)*d3^(1 - n3)*d4^(1 - n4)*n1)/(d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) -
        (3*d1^(-1 - n1)*d3^(1 - n3)*n1)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n3)) + (d10^(-1 - n10)*d3^(1 - n3)*d9^(1 - n9)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n3)) +
        (d10^(-1 - n10)*d3^(1 - n3)*d8^(1 - n8)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*
        d9^n9*(-1 + n3)) - (d10^(-1 - n10)*d3^(1 - n3)*d7^(1 - n7)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n3)) -
        (d10^(-1 - n10)*d3^(1 - n3)*d6^(1 - n6)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*
        d9^n9*(-1 + n3)) + (d10^(-1 - n10)*d3^(1 - n3)*d5^(1 - n5)*n10)/(d1^n1*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) +
        (d10^(-1 - n10)*d3^(1 - n3)*d4^(1 - n4)*n10)/(d1^n1*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n3)) - (d1^(1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*n10)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) +
        (d10^(-1 - n10)*d3^(1 - n3)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n3)) + (d2^(-1 - n2)*d3^(1 - n3)*d8^(1 - n8)*n2)/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n3)) -
        (d1^(1 - n1)*d2^(-1 - n2)*d3^(1 - n3)*n2)/(d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n3)) - (d2^(-1 - n2)*d3^(1 - n3)*n2)/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*
        d8^n8*d9^n9*(-1 + n3)) + (2*d3^(1 - n3)*d4^(1 - n4)*d5^(-1 - n5)*n5)/(d1^n1*d10^n10*d2^n2*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) -
        (2*d1^(1 - n1)*d3^(1 - n3)*d5^(-1 - n5)*n5)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n3)) + (d3^(1 - n3)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n3)) -
        (d3^(1 - n3)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(d1^n1*d10^n10*d2^n2*d4^n4*d7^n7*d8^n8*d9^n9*
        (-1 + n3)) + (d3^(1 - n3)*d4^(1 - n4)*d6^(-1 - n6)*n6)/(d1^n1*d10^n10*d2^n2*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n3)) -
        (d2^(1 - n2)*d3^(1 - n3)*d6^(-1 - n6)*n6)/(d1^n1*d10^n10*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*
        (-1 + n3)) + (d3^(1 - n3)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n3)) -
        (d3^(1 - n3)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d2^n2*d4^n4*d6^n6*d8^n8*d9^n9*
        (-1 + n3)) + (d3^(1 - n3)*d4^(1 - n4)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d2^n2*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n3)) -
        (d3^(2 - n3)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*
        (-1 + n3)) - (d3^(1 - n3)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n3)) +
        (d3^(1 - n3)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(d1^n1*d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d9^n9*
        (-1 + n3)) + (d2^(1 - n2)*d3^(1 - n3)*d8^(-1 - n8)*n8)/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n3)) -
        (d1^(1 - n1)*d3^(1 - n3)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*
        (-1 + n3)) - (d3^(1 - n3)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n3)) +
        (d3^(1 - n3)*d5^(1 - n5)*d9^(-1 - n9)*n9)/(d1^n1*d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*
        (-1 + n3)) - (d1^(1 - n1)*d3^(1 - n3)*d9^(-1 - n9)*n9)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n3)) +
        (d3^(2 - n3)*d9^(-1 - n9)*n9)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        (-1 + n3)) + (d3^(1 - n3)*rat(1 + d - 3*n1 - n2 - n3, -1 + n3))/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);
        
*  5: n10, n2 (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?{>1}/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg_ =
        
        (d10^(-1 - n10)*d9^(1 - n9))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) +
        (d10^(-1 - n10)*d8^(1 - n8))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) -
        (d1^(1 - n1)*d10^(-1 - n10))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d9^(1 - n9)*n1)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n2)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*n1)/(d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n2)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*n1)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n2)) + (d10^(-1 - n10)*d2^(1 - n2)*d4^(-1 - n4)*d7^(1 - n7)*n4)/(d1^n1*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n2)) -
        (d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*d4^(-1 - n4)*n4)/(d1^n1*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n2)) -
        (d10^(-1 - n10)*d2^(1 - n2)*d4^(-1 - n4)*n4)/(d1^n1*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n2)) - (d10^(-1 - n10)*d2^(1 - n2)*d4^(1 - n4)*d7^(-1 - n7)*n7)/(d1^n1*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n2)) +
        (d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*d7^(-1 - n7)*n7)/(d1^n1*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n2)) +
        (d10^(-1 - n10)*d2^(1 - n2)*d7^(-1 - n7)*n7)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*
        d9^n9*(-1 + n2)) + (d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*d9^(-1 - n9)*n9)/(d1^n1*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n2)) -
        (d1^(1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d9^(-1 - n9)*n9)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n2)) +
        (d10^(-1 - n10)*d2^(1 - n2)*d9^(-1 - n9)*n9)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d8^n8*(-1 + n2)) + (d10^(-1 - n10)*d2^(1 - n2)*(-n1 - n4 + n7 + n9))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n2));
        
*  6: n7      (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?{>1}/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ =
        
        d4^(1 - n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        d3^(1 - n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (2*d1^(-1 - n1)*d7^(1 - n7)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*
        (-1 + n7)) + (d10^(-1 - n10)*d7^(1 - n7)*d8^(1 - n8)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) -
        (d10^(-1 - n10)*d3^(1 - n3)*d7^(1 - n7)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*
        d9^n9*(-1 + n7)) - (d10^(-1 - n10)*d7^(1 - n7)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) -
        (d2^(-1 - n2)*d7^(1 - n7)*d8^(1 - n8)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*
        (-1 + n7)) + (d1^(1 - n1)*d2^(-1 - n2)*d7^(1 - n7)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) +
        (d2^(-1 - n2)*d7^(1 - n7)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*
        (-1 + n7)) - (2*d3^(-1 - n3)*d7^(1 - n7)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*
        d8^n8*d9^n9*(-1 + n7)) - (d4^(1 - n4)*d5^(-1 - n5)*d7^(1 - n7)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d6^n6*d8^n8*d9^n9*(-1 + n7)) +
        (d1^(1 - n1)*d5^(-1 - n5)*d7^(1 - n7)*n5)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*
        (-1 + n7)) + (d5^(-1 - n5)*d7^(1 - n7)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*
        d8^n8*d9^n9*(-1 + n7)) - (d6^(-1 - n6)*d7^(1 - n7)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d9^n9*(-1 + n7)) +
        (d5^(1 - n5)*d6^(-1 - n6)*d7^(1 - n7)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d8^n8*d9^n9*
        (-1 + n7)) - (d4^(1 - n4)*d6^(-1 - n6)*d7^(1 - n7)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d8^n8*d9^n9*(-1 + n7)) +
        (d2^(1 - n2)*d6^(-1 - n6)*d7^(1 - n7)*n6)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d8^n8*d9^n9*
        (-1 + n7)) + (d7^(1 - n7)*(1 + 2*n1 - n10 + n2 - 2*n3 + n5 - n7))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) -
        (2*d3^(1 - n3)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*
        d8^n8*(-1 + n7)) + (2*d1^(1 - n1)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7));
        
*  7: n10, n6 (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?{>1}/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg_ =
        
        (d10^(-1 - n10)*d9^(1 - n9))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) +
        (d10^(-1 - n10)*d8^(1 - n8))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) -
        (d10^(-1 - n10)*d7^(1 - n7))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9) +
        (d10^(-1 - n10)*d4^(1 - n4))/(d1^n1*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d10^(-1 - n10)*d3^(1 - n3))/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        (d1^(1 - n1)*d10^(-1 - n10))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        d10^(-1 - n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d6^(1 - n6)*d9^(1 - n9)*n1)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d6^(1 - n6)*n1)/(d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d6^(1 - n6)*n1)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*
        d9^n9*(-1 + n6)) + (d10^(-1 - n10)*d2^(-1 - n2)*d6^(1 - n6)*d9^(1 - n9)*n2)/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) +
        (d10^(-1 - n10)*d2^(-1 - n2)*d6^(1 - n6)*d8^(1 - n8)*n2)/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) -
        (d1^(1 - n1)*d10^(-1 - n10)*d2^(-1 - n2)*d6^(1 - n6)*n2)/(d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) -
        (d2^(-1 - n2)*d6^(1 - n6)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*
        (-1 + n6)) + (d10^(-1 - n10)*d5^(-1 - n5)*d6^(1 - n6)*d9^(1 - n9)*n5)/(d1^n1*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n6)) -
        (d10^(-1 - n10)*d5^(-1 - n5)*d6^(1 - n6)*d7^(1 - n7)*n5)/(d1^n1*d2^n2*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n6)) +
        (d10^(-1 - n10)*d4^(1 - n4)*d5^(-1 - n5)*d6^(1 - n6)*n5)/(d1^n1*d2^n2*d3^n3*d7^n7*d8^n8*d9^n9*(-1 + n6)) -
        (d1^(1 - n1)*d10^(-1 - n10)*d5^(-1 - n5)*d6^(1 - n6)*n5)/(d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n6)) +
        (d10^(-1 - n10)*d3^(1 - n3)*d6^(1 - n6)*d9^(-1 - n9)*n9)/(d1^n1*d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) -
        (d1^(1 - n1)*d10^(-1 - n10)*d6^(1 - n6)*d9^(-1 - n9)*n9)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) +
        (d10^(-1 - n10)*d6^(1 - n6)*d9^(-1 - n9)*n9)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*
        d8^n8*(-1 + n6)) + (d10^(-1 - n10)*d6^(1 - n6)*(-n1 + n9))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6));
        
*  8: n10, n8 (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?{>1}/d9^n9?pos_/d10^n10?neg_ =
        
        (d10^(-1 - n10)*d3^(1 - n3))/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        d10^(-1 - n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d8^(1 - n8)*d9^(1 - n9)*n1)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n8)) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d8^(1 - n8)*n1)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d8^(1 - n8)*n1)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d9^n9*(-1 + n8)) + (d10^(-2 - n10)*d3^(1 - n3)*d8^(1 - n8)*(-1 - n10))/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) +
        (d10^(-2 - n10)*d8^(1 - n8)*(-1 - n10))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d9^n9*(-1 + n8)) + (d10^(-2 - n10)*d8^(2 - n8)*(1 + n10))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) -
        (d10^(-1 - n10)*d5^(-1 - n5)*d8^(1 - n8)*d9^(1 - n9)*n5)/(d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*(-1 + n8)) +
        (d10^(-1 - n10)*d5^(-1 - n5)*d7^(1 - n7)*d8^(1 - n8)*n5)/(d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d9^n9*(-1 + n8)) -
        (d10^(-1 - n10)*d4^(1 - n4)*d5^(-1 - n5)*d8^(1 - n8)*n5)/(d1^n1*d2^n2*d3^n3*d6^n6*d7^n7*d9^n9*(-1 + n8)) +
        (d1^(1 - n1)*d10^(-1 - n10)*d5^(-1 - n5)*d8^(1 - n8)*n5)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n8)) +
        (d10^(-1 - n10)*d8^(1 - n8)*(-2 + n1 - n10 + n8 - n9))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) -
        (d10^(-1 - n10)*d3^(1 - n3)*d8^(1 - n8)*d9^(-1 - n9)*n9)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n8)) +
        (d1^(1 - n1)*d10^(-1 - n10)*d8^(1 - n8)*d9^(-1 - n9)*n9)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n8)) -
        (d10^(-1 - n10)*d8^(1 - n8)*d9^(-1 - n9)*n9)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*
        d7^n7*(-1 + n8));
*  9: n10, n3 (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?{>1}/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg_ =
        
        (d10^(-1 - n10)*d9^(1 - n9))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) +
        (d10^(-1 - n10)*d8^(1 - n8))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) -
        (d1^(1 - n1)*d10^(-1 - n10))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d8^(1 - n8)*n1)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n3)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*n1)/(d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*n1)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n3)) + (d10^(-1 - n10)*d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(d1^n1*d2^n2*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n3)) -
        (d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*d4^(-1 - n4)*n4)/(d1^n1*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) -
        (d10^(-1 - n10)*d3^(1 - n3)*d4^(-1 - n4)*n4)/(d1^n1*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n3)) - (d10^(-1 - n10)*d3^(1 - n3)*d4^(1 - n4)*d6^(-1 - n6)*n6)/(d1^n1*d2^n2*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n3)) +
        (d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*d6^(-1 - n6)*n6)/(d1^n1*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n3)) +
        (d10^(-1 - n10)*d3^(1 - n3)*d6^(-1 - n6)*n6)/(d1^n1*d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*
        d9^n9*(-1 + n3)) + (d10^(-1 - n10)*d2^(1 - n2)*d3^(1 - n3)*d8^(-1 - n8)*n8)/(d1^n1*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n3)) -
        (d1^(1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d8^(-1 - n8)*n8)/(d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n3)) +
        (d10^(-1 - n10)*d3^(1 - n3)*d8^(-1 - n8)*n8)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*
        d9^n9*(-1 + n3)) + (d10^(-1 - n10)*d3^(1 - n3)*(-n1 - n4 + n6 + n8))/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3));
        
* 10: n10, n7 (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?{>1}/d8^n8?pos_/d9^n9?pos_/d10^n10?neg_ =
        
        (d10^(-1 - n10)*d9^(1 - n9))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) +
        (d10^(-1 - n10)*d8^(1 - n8))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) -
        (d10^(-1 - n10)*d6^(1 - n6))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9) +
        (d10^(-1 - n10)*d4^(1 - n4))/(d1^n1*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d10^(-1 - n10)*d2^(1 - n2))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        (d1^(1 - n1)*d10^(-1 - n10))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        d10^(-1 - n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d7^(1 - n7)*d8^(1 - n8)*n1)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d7^(1 - n7)*n1)/(d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d7^(1 - n7)*n1)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*
        d9^n9*(-1 + n7)) + (d10^(-1 - n10)*d3^(-1 - n3)*d7^(1 - n7)*d9^(1 - n9)*n3)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) +
        (d10^(-1 - n10)*d3^(-1 - n3)*d7^(1 - n7)*d8^(1 - n8)*n3)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) -
        (d1^(1 - n1)*d10^(-1 - n10)*d3^(-1 - n3)*d7^(1 - n7)*n3)/(d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) -
        (d3^(-1 - n3)*d7^(1 - n7)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*
        (-1 + n7)) + (d10^(-1 - n10)*d5^(-1 - n5)*d7^(1 - n7)*d8^(1 - n8)*n5)/(d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d9^n9*(-1 + n7)) -
        (d10^(-1 - n10)*d5^(-1 - n5)*d6^(1 - n6)*d7^(1 - n7)*n5)/(d1^n1*d2^n2*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n7)) +
        (d10^(-1 - n10)*d4^(1 - n4)*d5^(-1 - n5)*d7^(1 - n7)*n5)/(d1^n1*d2^n2*d3^n3*d6^n6*d8^n8*d9^n9*(-1 + n7)) -
        (d1^(1 - n1)*d10^(-1 - n10)*d5^(-1 - n5)*d7^(1 - n7)*n5)/(d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n7)) +
        (d10^(-1 - n10)*d2^(1 - n2)*d7^(1 - n7)*d8^(-1 - n8)*n8)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) -
        (d1^(1 - n1)*d10^(-1 - n10)*d7^(1 - n7)*d8^(-1 - n8)*n8)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) +
        (d10^(-1 - n10)*d7^(1 - n7)*d8^(-1 - n8)*n8)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*
        d9^n9*(-1 + n7)) + (d10^(-1 - n10)*d7^(1 - n7)*(-n1 + n8))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7));
        
* 11: n6      (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?{>1}/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ =
        
        d8^(1 - n8)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) -
        d5^(1 - n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) +
        (2*d4^(1 - n4))/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        (2*d2^(1 - n2))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        (2*d10^(-1 - n10)*d6^(1 - n6)*d9^(1 - n9)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*
        d8^n8*(-1 + n6)) + (d10^(-1 - n10)*d6^(1 - n6)*d8^(1 - n8)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) -
        (d10^(-1 - n10)*d3^(1 - n3)*d6^(1 - n6)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*
        d9^n9*(-1 + n6)) - (2*d10^(-1 - n10)*d2^(1 - n2)*d6^(1 - n6)*n10)/(d1^n1*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) -
        (d10^(-1 - n10)*d6^(1 - n6)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*
        (-1 + n6)) + (d2^(-1 - n2)*d6^(1 - n6)*d8^(1 - n8)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) -
        (d1^(1 - n1)*d2^(-1 - n2)*d6^(1 - n6)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*
        (-1 + n6)) - (3*d2^(-1 - n2)*d6^(1 - n6)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*
        d8^n8*d9^n9*(-1 + n6)) - (2*d2^(1 - n2)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) +
        (2*d1^(1 - n1)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*
        d9^n9*(-1 + n6)) + (d6^(1 - n6)*rat(1 + d - n10 - 3*n2 - n6, -1 + n6))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9);
        
* 12: n10, n9 (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?{>1}/d10^n10?neg_ =
        
        (d10^(-1 - n10)*d2^(1 - n2))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) +
        d10^(-1 - n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        (d1^(-1 - n1)*d10^(-1 - n10)*d8^(1 - n8)*d9^(1 - n9)*n1)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d2^(1 - n2)*d9^(1 - n9)*n1)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) +
        (d1^(-1 - n1)*d10^(-1 - n10)*d9^(1 - n9)*n1)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d8^n8*(-1 + n9)) + (d10^(-2 - n10)*d2^(1 - n2)*d9^(1 - n9)*(-1 - n10))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) +
        (d10^(-2 - n10)*d9^(1 - n9)*(-1 - n10))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d8^n8*(-1 + n9)) + (d10^(-2 - n10)*d9^(2 - n9)*(1 + n10))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) -
        (d10^(-1 - n10)*d5^(-1 - n5)*d8^(1 - n8)*d9^(1 - n9)*n5)/(d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*(-1 + n9)) +
        (d10^(-1 - n10)*d5^(-1 - n5)*d6^(1 - n6)*d9^(1 - n9)*n5)/(d1^n1*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n9)) -
        (d10^(-1 - n10)*d4^(1 - n4)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*(-1 + n9)) +
        (d1^(1 - n1)*d10^(-1 - n10)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) -
        (d10^(-1 - n10)*d2^(1 - n2)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) +
        (d1^(1 - n1)*d10^(-1 - n10)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) -
        (d10^(-1 - n10)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*
        d7^n7*(-1 + n9)) + (d10^(-1 - n10)*d9^(1 - n9)*(-2 + n1 - n10 - n8 + n9))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9));
        
* 13: n4      (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?{>1}/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ =
        
        d5^(1 - n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) -
        d1^(1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        (2*d1^(-1 - n1)*d4^(1 - n4)*n1)/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n4)) + (d10^(-1 - n10)*d4^(1 - n4)*d9^(1 - n9)*n10)/(d1^n1*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n4)) +
        (d10^(-1 - n10)*d4^(1 - n4)*d8^(1 - n8)*n10)/(d1^n1*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*
        d9^n9*(-1 + n4)) - (d10^(-1 - n10)*d3^(1 - n3)*d4^(1 - n4)*n10)/(d1^n1*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) -
        (d10^(-1 - n10)*d2^(1 - n2)*d4^(1 - n4)*n10)/(d1^n1*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n4)) + (d2^(-1 - n2)*d4^(1 - n4)*d8^(1 - n8)*n2)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n4)) -
        (d1^(1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*n2)/(d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n4)) - (d2^(-1 - n2)*d4^(1 - n4)*n2)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*
        d8^n8*d9^n9*(-1 + n4)) + (d3^(-1 - n3)*d4^(1 - n4)*d9^(1 - n9)*n3)/(d1^n1*d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n4)) -
        (d1^(1 - n1)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n4)) - (d3^(-1 - n3)*d4^(1 - n4)*n3)/(d1^n1*d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*
        d8^n8*d9^n9*(-1 + n4)) + (d4^(1 - n4)*rat(1 + d - 2*n1 - n2 - n3 - n4, -1 + n4))/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);
        
* 14: n5      (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?{>1}/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?neg0_ =
        
        d4^(1 - n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        d1^(1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        (2*d1^(-1 - n1)*d5^(1 - n5)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n5)) + (d10^(-1 - n10)*d5^(1 - n5)*d9^(1 - n9)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n5)) +
        (d10^(-1 - n10)*d5^(1 - n5)*d8^(1 - n8)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*
        d9^n9*(-1 + n5)) - (d10^(-1 - n10)*d3^(1 - n3)*d5^(1 - n5)*n10)/(d1^n1*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) -
        (d10^(-1 - n10)*d2^(1 - n2)*d5^(1 - n5)*n10)/(d1^n1*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*
        d9^n9*(-1 + n5)) + (d2^(-1 - n2)*d5^(1 - n5)*d8^(1 - n8)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n5)) -
        (d1^(1 - n1)*d2^(-1 - n2)*d5^(1 - n5)*n2)/(d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n5)) - (d2^(-1 - n2)*d5^(1 - n5)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*
        d8^n8*d9^n9*(-1 + n5)) + (d3^(-1 - n3)*d5^(1 - n5)*d9^(1 - n9)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n5)) -
        (d1^(1 - n1)*d3^(-1 - n3)*d5^(1 - n5)*n3)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*
        (-1 + n5)) - (d3^(-1 - n3)*d5^(1 - n5)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*
        d8^n8*d9^n9*(-1 + n5)) + (d5^(1 - n5)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d9^n9*(-1 + n5)) +
        (d4^(1 - n4)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d7^n7*d8^n8*d9^n9*
        (-1 + n5)) - (d2^(1 - n2)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(d1^n1*d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n5)) -
        (d5^(2 - n5)*d6^(-1 - n6)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*
        (-1 + n5)) + (d5^(1 - n5)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*(-1 + n5)) +
        (d4^(1 - n4)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d6^n6*d8^n8*d9^n9*
        (-1 + n5)) - (d3^(1 - n3)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d2^n2*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n5)) -
        (d5^(2 - n5)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*
        (-1 + n5)) + (d5^(1 - n5)*rat(1 + d - 2*n1 - n2 - n3 - n5, -1 + n5))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9);
        
* 15: n8      (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?{>1}/d9^n9?pos_/d10^n10?neg0_ =
        
        d2^(1 - n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        d1^(1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        (d10^(-1 - n10)*d8^(1 - n8)*d9^(1 - n9)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*
        d7^n7*(-1 + n8)) + (d10^(-1 - n10)*d3^(1 - n3)*d8^(1 - n8)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) +
        (d10^(-1 - n10)*d2^(1 - n2)*d8^(1 - n8)*n10)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        d9^n9*(-1 + n8)) - (d10^(-1 - n10)*d8^(2 - n8)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) +
        (d1^(1 - n1)*d2^(-1 - n2)*d8^(1 - n8)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*
        (-1 + n8)) + (d2^(-1 - n2)*d8^(1 - n8)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*
        d7^n7*d9^n9*(-1 + n8)) - (d2^(-1 - n2)*d8^(2 - n8)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) +
        (d5^(1 - n5)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d9^n9*
        (-1 + n8)) - (d4^(1 - n4)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d9^n9*(-1 + n8)) +
        (d2^(1 - n2)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*
        (-1 + n8)) - (d6^(-1 - n6)*d8^(2 - n8)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*
        d7^n7*d9^n9*(-1 + n8)) + (d8^(1 - n8)*(1 + n2 - n8))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8));
        
* 16: n9      (-)
        id,ifmatch->sortme 1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?{>1}/d10^n10?neg0_ =
        
        d3^(1 - n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        d1^(1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) -
        (2*d1^(-1 - n1)*d9^(1 - n9)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        (-1 + n9)) - (d10^(-1 - n10)*d8^(1 - n8)*d9^(1 - n9)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) +
        (d10^(-1 - n10)*d3^(1 - n3)*d9^(1 - n9)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*
        d8^n8*(-1 + n9)) + (d10^(-1 - n10)*d2^(1 - n2)*d9^(1 - n9)*n10)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) -
        (d10^(-1 - n10)*d9^(2 - n9)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*
        (-1 + n9)) + (d4^(1 - n4)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*(-1 + n9)) -
        (d1^(1 - n1)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*
        (-1 + n9)) - (d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*
        d7^n7*d8^n8*(-1 + n9)) + (d2^(1 - n2)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) -
        (d1^(1 - n1)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*
        (-1 + n9)) - (d8^(-1 - n8)*d9^(1 - n9)*n8)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*
        d6^n6*d7^n7*(-1 + n9)) + (d9^(1 - n9)*rat(1 + d - 2*n1 - n5 - n8 - n9, -1 + n9))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8);

* 17: n10     (-)              
        id,only,ifmatch->sortme 1/d1/d2/d3/d4/d5/d6/d7/d8/d9/d10^n10?neg_ =

        (d10^(-1 - n10)*rat(-4, -4 + d - n10))/(d1*d2*d3*d4*d5^2*d6*d7*d9) +
        (d10^(-1 - n10)*rat(-4, -4 + d - n10))/(d1^2*d2*d3*d4*d5*d6*d7*d9) +
        (d10^(-1 - n10)*rat(-4, -4 + d - n10))/(d1*d3*d4*d5*d6*d7*d8^2*d9) +
        (d10^(-1 - n10)*rat(-4, -4 + d - n10))/(d1*d2*d3*d4*d5*d6*d7*d8^2*d9) +
        (d10^(-1 - n10)*rat(-4, -4 + d - n10))/(d1*d2*d3*d5^2*d6*d7*d8*d9) +
        rat(-3, -4 + d - n10)/(d1*d10^n10*d2*d3*d4*d6*d7^2*d8*d9) +
        rat(-3, -4 + d - n10)/(d10^n10*d2*d3^2*d4*d5*d6*d7*d8*d9) +
        rat(1, -4 + d - n10)/(d1*d10^n10*d2*d3*d4*d6*d7*d8*d9^2) +
        rat(1, -4 + d - n10)/(d10^n10*d2*d3*d4*d5*d6*d7*d8*d9^2) +
        rat(1, -4 + d - n10)/(d1*d10^n10*d2*d3*d5*d6*d7^2*d8*d9) +
        rat(1, -4 + d - n10)/(d1*d10^n10*d2*d3^2*d5*d6*d7*d8*d9) +
        rat(1, 4 - d + n10)/(d1*d10^n10*d2*d3*d4*d5*d6*d8*d9^2) +
        rat(1, 4 - d + n10)/(d1*d10^n10*d2*d4*d5*d6*d7*d8*d9^2) +
        rat(1, 4 - d + n10)/(d1*d10^n10*d2*d3^2*d4*d5*d6*d8*d9) +
        rat(1, 4 - d + n10)/(d1*d10^n10*d2*d4*d5*d6*d7^2*d8*d9) +
        rat(3, -4 + d - n10)/(d1*d10^n10*d2*d3*d4*d5*d6*d7^2*d8) +
        rat(3, -4 + d - n10)/(d1*d10^n10*d2*d3^2*d4*d5*d6*d7*d8) +
        (d10^(-1 - n10)*rat(4, -4 + d - n10))/(d1*d3*d4*d5*d6*d7*d8*d9^2) +
        (d10^(-1 - n10)*rat(4, -4 + d - n10))/(d1*d2*d3*d4*d5*d6*d7*d8*d9^2) +
        (d10^(-1 - n10)*rat(4, -4 + d - n10))/(d2*d3*d4*d5*d6*d7*d8^2*d9) +
        (d10^(-1 - n10)*rat(4, -4 + d - n10))/(d1*d2*d3*d4*d5^2*d7*d8*d9) +
        (d10^(-1 - n10)*rat(4, -4 + d - n10))/(d2*d3*d4*d5^2*d6*d7*d8*d9) +
        (d10^(-1 - n10)*rat(4, -4 + d - n10))/(d1^2*d3*d4*d5*d6*d7*d8*d9) +
        (d10^(-1 - n10)*rat(4, -4 + d - n10))/(d1^2*d2*d3*d4*d5*d6*d7*d8*d9) +
        (d10^(-2 - n10)*rat(-4 - 4*n10, -4 + d - n10))/(d1*d3*d4*d5*d6*d7*d8*d9) +
        (d10^(-2 - n10)*rat(-4 - 4*n10, -4 + d - n10))/(d1*d2*d3*d4*d5*d6*d7*d8*d9) +
        (d10^(-2 - n10)*rat(-4 - 4*n10, 4 - d + n10))/(d1*d2*d3*d4*d5*d6*d7*d8) +
        (d10^(-1 - n10)*rat(-2*n10, -4 + d - n10))/(d1*d3*d4*d5*d6*d7*d8*d9) +
        (d10^(-1 - n10)*rat(-2*n10, -4 + d - n10))/(d1*d2*d3*d4*d5*d6*d7*d8*d9) +
        (d10^(-1 - n10)*rat(-n10, -4 + d - n10))/(d1*d2*d3*d4*d5*d6*d8*d9) +
        (d10^(-1 - n10)*rat(-n10, -4 + d - n10))/(d1*d2*d3*d4*d5*d7*d8*d9) +
        (d10^(-1 - n10)*rat(-n10, -4 + d - n10))/(d1*d2*d4*d5*d6*d7*d8*d9) +
        (d10^(-1 - n10)*rat(-n10, -4 + d - n10))/(d2*d3*d4*d5*d6*d7*d8*d9) +
        (d10^(-1 - n10)*rat(n10, -4 + d - n10))/(d1*d2*d3*d4*d6*d7*d8*d9) +
        (d10^(-1 - n10)*rat(n10, -4 + d - n10))/(d1*d2*d3*d5*d6*d7*d8*d9) +
        (d10^(-1 - n10)*rat(2*n10, -4 + d - n10))/(d1*d2*d3*d4*d5*d6*d7*d9) +
        (d10^(-1 - n10)*rat(3*n10, -4 + d - n10))/(d1*d2*d3*d4*d5*d6*d7*d8);

        if(count(d10,1)==0) id,only 1/d1/d2/d3/d4/d5/d6/d7/d8/d9 = PR12/intH;        

        #call zeroH        

        endif;

        goto endrec;         
        la sortme;
        $irep = 0;
        la endrec;
        
        ModuleOption,minimum,$irep;
        .sort:redH-`$repcount++';
        #redefine irep "`$irep'"
#enddo
#endprocedure
*--#] redH :

* 
* Planar 8 line
* 
*--#[ redBMW :
#procedure redBMW
#$repcount = 1;        
#do irep=1,1
        #$irep = 1;                
        if(count(intBMW,1));
* [0011111111]                  
* den: minus
* n7 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?{>1}/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        -(d1^(-1 - n1)*d7^(1 - n7)*d9^(1 - n9)*n1)/(2*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) + 
        (d1^(-1 - n1)*d5^(1 - n5)*d7^(1 - n7)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d1^(-1 - n1)*d4^(1 - n4)*d7^(1 - n7)*n1)/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d7^(1 - n7)*n1)/(2*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d1^(-1 - n1)*d7^(1 - n7)*n1)/(2*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d2^(-1 - n2)*d7^(1 - n7)*d9^(1 - n9)*n2)/(2*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7)) - 
        (d2^(-1 - n2)*d7^(1 - n7)*d8^(1 - n8)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) + 
        (d2^(-1 - n2)*d6^(1 - n6)*d7^(1 - n7)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d8^n8*d9^n9*(-1 + n7)) - 
        (d2^(-1 - n2)*d4^(1 - n4)*d7^(1 - n7)*n2)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d7^(1 - n7)*n2)/(2*d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d1^(1 - n1)*d2^(-1 - n2)*d7^(1 - n7)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d2^(-1 - n2)*d7^(1 - n7)*n2)/(2*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d3^(1 - n3)*d4^(-1 - n4)*d7^(1 - n7)*n4)/(2*d1^n1*d10^n10*d2^n2*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d4^(-1 - n4)*d7^(1 - n7)*n4)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d4^(-1 - n4)*d7^(2 - n7)*n4)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d5^(-1 - n5)*d7^(1 - n7)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d6^(-1 - n6)*d7^(1 - n7)*d8^(1 - n8)*n6)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d9^n9*(-1 + n7)) + 
        (d5^(1 - n5)*d6^(-1 - n6)*d7^(1 - n7)*n6)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n7)) + 
        (d6^(-1 - n6)*d7^(1 - n7)*n6)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d8^n8*d9^n9*(-1 + n7)) + 
        (d7^(1 - n7)*(2 + n1 - n2 - n4 + 2*n5 + n6 - 2*n7))/(2*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d5^(1 - n5)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*(-1 + n7)) - 
        (d7^(2 - n7)*d9^(-1 - n9)*n9)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n7));

* n1 != 0 && n3 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?neg0_/d3^n3?{>1}/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        (d1^(-1 - n1)*d9^(1 - n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) - 
        (d1^(-1 - n1)*d7^(1 - n7))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d4^(1 - n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d3^(1 - n3)*d5^(1 - n5)*(-1 - n1))/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) + 
        (d1^(-2 - n1)*d3^(1 - n3)*(-1 - n1))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) + 
        (d1^(-2 - n1)*d3^(1 - n3)*d4^(1 - n4)*(1 + n1))/(d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) - 
        (2*d1^(-1 - n1)*d3^(1 - n3)*d5^(-1 - n5)*n5)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n3)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n3)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n3)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n3)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n3)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n3)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n3)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*rat(-1 + d - n1 - 2*n5 - n7 - n8, -1 + n3))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

* n9 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?{>1}/d10^n10?pos_ =

        d7^(1 - n7)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9) - 
        d5^(1 - n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d1^(-1 - n1)*d5^(1 - n5)*d9^(1 - n9)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*d9^(1 - n9)*n1)/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d1^(-1 - n1)*d9^(1 - n9)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d2^(-1 - n2)*d8^(1 - n8)*d9^(1 - n9)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) - 
        (d2^(-1 - n2)*d6^(1 - n6)*d9^(1 - n9)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n9)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*d9^(1 - n9)*n2)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d1^(1 - n1)*d2^(-1 - n2)*d9^(1 - n9)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (2*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d6^(-1 - n6)*d8^(1 - n8)*d9^(1 - n9)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*(-1 + n9)) - 
        (d5^(1 - n5)*d6^(-1 - n6)*d9^(1 - n9)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n9)) - 
        (d6^(-1 - n6)*d9^(1 - n9)*n6)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n9)) + 
        (d9^(1 - n9)*rat(1 + d - n1 - 2*n5 - n6 - n9, -1 + n9))/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8);

* n1 != 0 && n10 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?{>1} =
        
        (d1^(-1 - n1)*d9^(1 - n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) + 
        (d1^(-1 - n1)*d8^(1 - n8))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) - 
        (d1^(-1 - n1)*d7^(1 - n7))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9) - 
        (d1^(-1 - n1)*d6^(1 - n6))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d5^(1 - n5))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d4^(1 - n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d10^(1 - n10)*d4^(1 - n4)*(-1 - n1))/(d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d1^(-2 - n1)*d10^(1 - n10)*d5^(1 - n5)*(1 + n1))/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d1^(-2 - n1)*d10^(1 - n10)*(1 + n1))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (2*d1^(-1 - n1)*d10^(1 - n10)*d5^(-1 - n5)*n5)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d8^(-1 - n8)*n8)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d7^(1 - n7)*d9^(-1 - n9)*n9)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d5^(1 - n5)*d9^(-1 - n9)*n9)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d9^(-1 - n9)*n9)/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*rat(1 - d + n1 + 2*n5 + n8 + n9, -1 + n10))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

* n8 != 1

        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?{>1}/d9^n9?pos_/d10^n10?pos_ =

        -(d1^(-1 - n1)*d8^(1 - n8)*d9^(1 - n9)*n1)/(2*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n8)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d8^(1 - n8)*n1)/(2*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d1^(-1 - n1)*d8^(1 - n8)*n1)/(2*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d10^(-1 - n10)*d3^(1 - n3)*d8^(1 - n8)*n10)/(d1^n1*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d10^(-1 - n10)*d8^(2 - n8)*n10)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d2^(-1 - n2)*d8^(1 - n8)*d9^(1 - n9)*n2)/(2*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n8)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d8^(1 - n8)*n2)/(2*d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d1^(1 - n1)*d2^(-1 - n2)*d8^(1 - n8)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d2^(-1 - n2)*d8^(1 - n8)*n2)/(2*d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d2^(-1 - n2)*d8^(2 - n8)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d3^(-1 - n3)*d8^(1 - n8)*n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d4^(-1 - n4)*d7^(1 - n7)*d8^(1 - n8)*n4)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d9^n9*(-1 + n8)) + 
        (d3^(1 - n3)*d4^(-1 - n4)*d8^(1 - n8)*n4)/(2*d1^n1*d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d4^(-1 - n4)*d8^(1 - n8)*n4)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d5^(1 - n5)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d9^n9*(-1 + n8)) - 
        (d6^(-1 - n6)*d8^(1 - n8)*n6)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n8)) - 
        (d6^(-1 - n6)*d8^(2 - n8)*n6)/(2*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n8)) + 
        (d8^(1 - n8)*(2 + n1 - n2 + 2*n3 + n4 - n6 - 2*n8))/(2*d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8));

* n1 != 0 && n4 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?neg0_/d3^n3?pos_/d4^n4?{>1}/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        (d1^(-1 - n1)*d5^(1 - n5))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) - 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d4^(1 - n4)*(-2 - 2*n1))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d4^(1 - n4)*d8^(1 - n8)*n10)/(d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n4)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d4^(1 - n4)*n10)/(d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d4^(1 - n4)*n10)/(d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*d8^(1 - n8)*n2)/(d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n4)) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*n2)/(d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4)*n2)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) - 
        (2*d1^(-1 - n1)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d10^n10*d2^n2*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) - 
        (d1^(-1 - n1)*d4^(1 - n4)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d8^n8*(-1 + n4)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d3^n3*d6^n6*d8^n8*d9^n9*(-1 + n4)) - 
        (d1^(-1 - n1)*d4^(1 - n4)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n4)) - 
        (2*d1^(-1 - n1)*d4^(1 - n4)*d9^(-1 - n9)*n9)/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n4)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*rat(-1 + 2*d - 2*n1 - n10 - n2 - 2*n3 - n4 - n7 - 2*n9, -1 + n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

* n1 != 0 && n6 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?{>1}/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        (d1^(-1 - n1)*d9^(1 - n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) + 
        (d1^(-1 - n1)*d8^(1 - n8))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) - 
        (d1^(-1 - n1)*d7^(1 - n7))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d4^(1 - n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d3^(1 - n3))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d1^(-1 - n1)*d10^(1 - n10))/(d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d6^(1 - n6)*d8^(1 - n8)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d6^(1 - n6)*n10)/(d2^n2*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d6^(1 - n6)*n10)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d6^(1 - n6)*d9^(1 - n9)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n6)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d6^(1 - n6)*d8^(1 - n8)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d2^(-1 - n2)*d6^(1 - n6)*n2)/(d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d2^(-1 - n2)*d6^(1 - n6)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d3^(1 - n3)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) - 
        (d1^(-1 - n1)*d10^(1 - n10)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d6^(1 - n6)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n6)) + 
        (d1^(-1 - n1)*d6^(1 - n6)*(-n10 + n8))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6));

* n1 != 0 && n5 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?{>1}/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        (d1^(-1 - n1)*d9^(1 - n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) - 
        (d1^(-1 - n1)*d7^(1 - n7))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d4^(1 - n4))/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d3^(1 - n3)*d5^(1 - n5)*(-1 - n1))/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-2 - n1)*d5^(1 - n5)*d9^(1 - n9)*(1 + n1))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n5)) + 
        (d1^(-2 - n1)*d5^(1 - n5)*(1 + n1))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d10^(-1 - n10)*d5^(1 - n5)*d8^(1 - n8)*n10)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d3^(1 - n3)*d5^(1 - n5)*n10)/(d2^n2*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d10^(-1 - n10)*d5^(1 - n5)*n10)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d5^(1 - n5)*d8^(1 - n8)*n2)/(d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d5^(1 - n5)*n2)/(d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d2^(-1 - n2)*d5^(1 - n5)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d4^(-1 - n4)*d5^(1 - n5)*n4)/(d10^n10*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d4^(-1 - n4)*d5^(1 - n5)*n4)/(d1^n1*d10^n10*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d4^(-1 - n4)*d5^(2 - n5)*n4)/(d10^n10*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d5^(1 - n5)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*(-1 + n5)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d3^n3*d6^n6*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(1 - n5)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d5^(2 - n5)*d7^(-1 - n7)*n7)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d3^(1 - n3)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n5)) + 
        (d1^(-1 - n1)*d10^(1 - n10)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n5)) - 
        (d1^(-1 - n1)*d5^(1 - n5)*d8^(-1 - n8)*n8)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n5)) + 
        (2*d1^(-1 - n1)*d5^(1 - n5)*d9^(-1 - n9)*n9)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n5)) + 
        (d1^(-1 - n1)*d5^(1 - n5)*rat(1 - d + n1 + n10 + n2 + n4 - n8 + 2*n9, -1 + n5))/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9);

* n2 != 0 && n4 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg_/d3^n3?pos_/d4^n4?{>1}/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        (d2^(-1 - n2)*d6^(1 - n6))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9) - 
        (d2^(-1 - n2)*d5^(1 - n5))/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(1 - n1)*d2^(-1 - n2))/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*d8^(1 - n8)*n1)/(d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n4)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d4^(1 - n4)*n1)/(d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) - 
        (d1^(-1 - n1)*d4^(1 - n4)*n1)/(d10^n10*d2^n2*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*d8^(1 - n8)*n3)/(d1^n1*d10^n10*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n4)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d1^n1*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d2^(-1 - n2)*d3^(-1 - n3)*d4^(1 - n4)*n3)/(d1^n1*d10^n10*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d10^n10*d3^n3*d5^n5*d7^n7*d9^n9*(-1 + n4)) - 
        (d2^(-1 - n2)*d4^(1 - n4)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(d1^n1*d10^n10*d3^n3*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*d6^(-1 - n6)*n6)/(d1^n1*d10^n10*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n4)) + 
        (2*d2^(-1 - n2)*d4^(1 - n4)*d8^(-1 - n8)*n8)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n4)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*rat(-d + n1 + n3 + n6 + 2*n8, -1 + n4))/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

* n2 != 0 && n7 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?{>1}/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        -((d2^(-1 - n2)*d8^(1 - n8))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9)) + 
        (d2^(-1 - n2)*d6^(1 - n6))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9) - 
        (d2^(-1 - n2)*d5^(1 - n5))/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d10^(1 - n10)*d2^(-1 - n2))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(1 - n1)*d2^(-1 - n2))/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        d2^(-1 - n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d7^(1 - n7)*d8^(1 - n8)*n1)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d7^(1 - n7)*n1)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) + 
        (d1^(-1 - n1)*d7^(1 - n7)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d2^(-1 - n2)*d3^(-1 - n3)*d7^(1 - n7)*d8^(1 - n8)*n3)/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(-1 - n3)*d7^(1 - n7)*n3)/(d1^n1*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d2^(-1 - n2)*d3^(-1 - n3)*d7^(1 - n7)*n3)/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (d2^(-1 - n2)*d5^(-1 - n5)*d7^(1 - n7)*d8^(1 - n8)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d9^n9*(-1 + n7)) + 
        (d2^(-1 - n2)*d5^(-1 - n5)*d6^(1 - n6)*d7^(1 - n7)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n7)) - 
        (d2^(-1 - n2)*d5^(-1 - n5)*d7^(1 - n7)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d8^n8*d9^n9*(-1 + n7)) - 
        (2*d2^(-1 - n2)*d7^(1 - n7)*d8^(-1 - n8)*n8)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n7)) + 
        (d2^(-1 - n2)*d7^(1 - n7)*rat(d - n1 - n3 - n5 - 2*n8, -1 + n7))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d8^n8*d9^n9);

* n10 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?{>1} =
        
        -(d8^(1 - n8)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9)) + 
        d3^(1 - n3)/(d1^n1*d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d8^(1 - n8)*n2)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) + 
        (d1^(1 - n1)*d10^(1 - n10)*d2^(-1 - n2)*n2)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*n2)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d6^(-1 - n6)*d8^(1 - n8)*n6)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d5^(1 - n5)*d6^(-1 - n6)*n6)/(d1^n1*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d6^(-1 - n6)*n6)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (2*d10^(1 - n10)*d8^(-1 - n8)*n8)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*rat(1 + d - n10 - n2 - n6 - 2*n8, -1 + n10))/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

* n2 != 0 && n9 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?{>1}/d10^n10?pos_ =
        
        -((d2^(-1 - n2)*d3^(1 - n3))/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9)) + 
        (d10^(1 - n10)*d2^(-1 - n2))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(1 - n1)*d2^(-1 - n2))/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d8^(1 - n8)*d9^(1 - n9)*n1)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d9^(1 - n9)*n1)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d1^(-1 - n1)*d9^(1 - n9)*n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d10^(-1 - n10)*d2^(-1 - n2)*d8^(1 - n8)*d9^(1 - n9)*n10)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) - 
        (d10^(-1 - n10)*d2^(-1 - n2)*d3^(1 - n3)*d9^(1 - n9)*n10)/(d1^n1*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d10^(-1 - n10)*d2^(-1 - n2)*d9^(1 - n9)*n10)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d2^(-1 - n2)*d5^(-1 - n5)*d8^(1 - n8)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*(-1 + n9)) - 
        (d2^(-1 - n2)*d5^(-1 - n5)*d6^(1 - n6)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n9)) + 
        (d2^(-1 - n2)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (2*d2^(-1 - n2)*d8^(-1 - n8)*d9^(1 - n9)*n8)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) + 
        (d2^(-1 - n2)*d9^(1 - n9)*rat(-d + n1 + n10 + n5 + 2*n8, -1 + n9))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8);

* n2 != 0 && n6 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?{>1}/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_ =
        
        (d2^(-1 - n2)*d8^(1 - n8))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) - 
        (d2^(-1 - n2)*d5^(1 - n5))/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d2^(-1 - n2)*d4^(1 - n4))/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        (d2^(-1 - n2)*d4^(-1 - n4)*d5^(1 - n5)*d6^(1 - n6)*n4)/(d1^n1*d10^n10*d3^n3*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d1^(1 - n1)*d2^(-1 - n2)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(d10^n10*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-1 - n2)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(d1^n1*d10^n10*d3^n3*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-1 - n2)*d6^(1 - n6)*(n4 - n5))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*d5^(-1 - n5)*d6^(1 - n6)*n5)/(d1^n1*d10^n10*d3^n3*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d1^(1 - n1)*d2^(-1 - n2)*d5^(-1 - n5)*d6^(1 - n6)*n5)/(d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n6)) - 
        (d2^(-1 - n2)*d5^(-1 - n5)*d6^(1 - n6)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d7^n7*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-1 - n2)*d6^(1 - n6)*d7^(-1 - n7)*d9^(1 - n9)*n7)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d8^n8*(-1 + n6)) - 
        (d2^(-1 - n2)*d5^(1 - n5)*d6^(1 - n6)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d3^n3*d4^n4*d8^n8*d9^n9*(-1 + n6)) + 
        (d2^(-1 - n2)*d4^(1 - n4)*d6^(1 - n6)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d3^n3*d5^n5*d8^n8*d9^n9*(-1 + n6)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*d6^(1 - n6)*d7^(-1 - n7)*n7)/(d1^n1*d10^n10*d4^n4*d5^n5*d8^n8*d9^n9*(-1 + n6));
        
* n2 != 0 && n8 != 1
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?{>1}/d9^n9?pos_/d10^n10?pos_ =
        
        (d1^(1 - n1)*d2^(-1 - n2))/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        d2^(-1 - n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (2*d1^(-1 - n1)*d2^(-1 - n2)*d8^(1 - n8)*n1)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d2^(-1 - n2)*d3^(-1 - n3)*d8^(1 - n8)*d9^(1 - n9)*n3)/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n8)) + 
        (d1^(1 - n1)*d2^(-1 - n2)*d3^(-1 - n3)*d8^(1 - n8)*n3)/(d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d2^(-1 - n2)*d3^(-1 - n3)*d8^(1 - n8)*n3)/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d2^(-1 - n2)*d4^(1 - n4)*d5^(-1 - n5)*d8^(1 - n8)*n5)/(d1^n1*d10^n10*d3^n3*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d1^(1 - n1)*d2^(-1 - n2)*d5^(-1 - n5)*d8^(1 - n8)*n5)/(d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n8)) + 
        (d2^(-1 - n2)*d5^(-1 - n5)*d8^(1 - n8)*n5)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d7^n7*d9^n9*(-1 + n8)) - 
        (d2^(-1 - n2)*d7^(-1 - n7)*d8^(1 - n8)*d9^(1 - n9)*n7)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*(-1 + n8)) + 
        (d2^(-1 - n2)*d5^(1 - n5)*d7^(-1 - n7)*d8^(1 - n8)*n7)/(d1^n1*d10^n10*d3^n3*d4^n4*d6^n6*d9^n9*(-1 + n8)) - 
        (d2^(-1 - n2)*d4^(1 - n4)*d7^(-1 - n7)*d8^(1 - n8)*n7)/(d1^n1*d10^n10*d3^n3*d5^n5*d6^n6*d9^n9*(-1 + n8)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*d7^(-1 - n7)*d8^(1 - n8)*n7)/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d9^n9*(-1 + n8)) + 
        (d2^(-1 - n2)*d8^(1 - n8)*rat(-1 - d + 2*n1 + n3 + n5 + n8, -1 + n8))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9);

* n1 != 0 && n9 != 1
        id,ifmatch->sortme 1/d1^n1?neg_/d2^n2?neg0_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?{>1}/d10^n10?pos_ =
        
        (d1^(-1 - n1)*d3^(1 - n3))/(d10^n10*d2^n2*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) - 
        d1^(-1 - n1)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (d1^(-2 - n1)*d9^(1 - n9)*(-2 - 2*n1))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*d8^(1 - n8)*d9^(1 - n9)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*(-1 + n9)) - 
        (d1^(-1 - n1)*d2^(-1 - n2)*d9^(1 - n9)*n2)/(d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d2^(-1 - n2)*d9^(1 - n9)*n2)/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d10^n10*d2^n2*d3^n3*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d1^(-1 - n1)*d5^(-1 - n5)*d9^(1 - n9)*n5)/(d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) - 
        (d5^(-1 - n5)*d9^(1 - n9)*n5)/(d1^n1*d10^n10*d2^n2*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*(-1 + n9)) + 
        (d1^(-1 - n1)*d6^(-1 - n6)*d8^(1 - n8)*d9^(1 - n9)*n6)/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d7^n7*(-1 + n9)) - 
        (d1^(-1 - n1)*d5^(1 - n5)*d6^(-1 - n6)*d9^(1 - n9)*n6)/(d10^n10*d2^n2*d3^n3*d4^n4*d7^n7*d8^n8*(-1 + n9)) + 
        (d1^(-1 - n1)*d4^(1 - n4)*d6^(-1 - n6)*d9^(1 - n9)*n6)/(d10^n10*d2^n2*d3^n3*d5^n5*d7^n7*d8^n8*(-1 + n9)) - 
        (d1^(-1 - n1)*d2^(1 - n2)*d6^(-1 - n6)*d9^(1 - n9)*n6)/(d10^n10*d3^n3*d4^n4*d5^n5*d7^n7*d8^n8*(-1 + n9)) + 
        (d1^(-1 - n1)*d9^(1 - n9)*rat(-1 + d - 2*n1 - n2 - n5 - n9, -1 + n9))/(d10^n10*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8);

* n10 != 1 && n2 != 0
        id,ifmatch->sortme 1/d1^n1?neg0_/d2^n2?neg_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?{>1} =
        
        (d2^(-1 - n2)*d9^(1 - n9))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8) + 
        (d2^(-1 - n2)*d8^(1 - n8))/(d1^n1*d10^n10*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(d1^n1*d10^n10*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9) + 
        (2*d1^(-1 - n1)*d10^(1 - n10)*d2^(-1 - n2)*n1)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d4^(1 - n4)*d5^(-1 - n5)*n5)/(d1^n1*d3^n3*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d1^(1 - n1)*d10^(1 - n10)*d2^(-1 - n2)*d5^(-1 - n5)*n5)/(d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d5^(-1 - n5)*n5)/(d1^n1*d3^n3*d4^n4*d6^n6*d7^n7*d8^n8*d9^n9*(-1 + n10)) + 
        (d1^(1 - n1)*d10^(1 - n10)*d2^(-1 - n2)*d8^(-1 - n8)*n8)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d8^(-1 - n8)*n8)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d8^(-1 - n8)*n8)/(d1^n1*d2^n2*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d9^n9*(-1 + n10)) - 
        (d10^(1 - n10)*d2^(-1 - n2)*d3^(1 - n3)*d9^(-1 - n9)*n9)/(d1^n1*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d1^(1 - n1)*d10^(1 - n10)*d2^(-1 - n2)*d9^(-1 - n9)*n9)/(d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*d9^(-1 - n9)*n9)/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*(-1 + n10)) + 
        (d10^(1 - n10)*d2^(-1 - n2)*rat(-d + 2*n1 + n5 + n8 + n9, -1 + n10))/(d1^n1*d3^n3*d4^n4*d5^n5*d6^n6*d7^n7*d8^n8*d9^n9);

* n1 == 0 && n10 == 1 && n2 == 0 && n7 == 1 && n8 == 1 && n9 == 1 && n5 != 1 && n6 != 1
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3^n3?pos_/d4^n4?pos_/d5^n5?{>1}/d6^n6?{>1}/d7/d8/d9/d10 =
        
        -1/(3*d10*d3^n3*d4^n4*d5^n5*d6^n6*d7*d9) - 
        (2*d3^(1 - n3))/(3*d10*d4^n4*d5^n5*d6^n6*d7*d8*d9) + 
        2/(3*d3^n3*d4^n4*d5^n5*d6^n6*d7*d8*d9) + 
        (2*d5^(1 - n5))/(3*d10*d3^n3*d4^n4*d6^n6*d8*d9^2*(-1 + n5)) - 
        (2*d5^(2 - n5))/(3*d10*d3^n3*d4^n4*d6^n6*d7*d8*d9^2*(-1 + n5)) + 
        (2*d5^(1 - n5))/(3*d10^2*d3^n3*d4^n4*d6^n6*d7*d9*(-1 + n5)) - 
        (2*d3^(1 - n3)*d5^(1 - n5))/(3*d10*d4^n4*d6^n6*d7*d8^2*d9*(-1 + n5)) + 
        (2*d5^(1 - n5))/(3*d3^n3*d4^n4*d6^n6*d7*d8^2*d9*(-1 + n5)) + 
        (2*d4^(1 - n4)*d5^(1 - n5))/(3*d10*d3^n3*d6^n6*d7^2*d8*d9*(-1 + n5)) - 
        (2*d3^(1 - n3)*d5^(1 - n5))/(3*d10*d4^n4*d6^n6*d7^2*d8*d9*(-1 + n5)) - 
        (2*d3^(1 - n3)*d5^(1 - n5))/(3*d10^2*d4^n4*d6^n6*d7*d8*d9*(-1 + n5)) - 
        (2*d3^(-1 - n3)*d5^(1 - n5)*n3)/(d10*d4^n4*d6^n6*d7*d8*d9*(-1 + n5)) + 
        (2*d4^(-1 - n4)*d5^(1 - n5)*n4)/(3*d10*d3^n3*d6^n6*d8*d9*(-1 + n5)) - 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d5^(1 - n5)*n4)/(3*d10*d6^n6*d7*d8*d9*(-1 + n5)) + 
        (2*d6^(1 - n6))/(3*d10^2*d3^n3*d4^n4*d5^n5*d7*d9*(-1 + n6)) - 
        (2*d3^(1 - n3)*d6^(1 - n6))/(3*d10*d4^n4*d5^n5*d7*d8^2*d9*(-1 + n6)) + 
        (2*d6^(1 - n6))/(3*d3^n3*d4^n4*d5^n5*d7*d8^2*d9*(-1 + n6)) - 
        (2*d3^(1 - n3)*d6^(1 - n6))/(3*d10^2*d4^n4*d5^n5*d7*d8*d9*(-1 + n6)) - 
        (2*d3^(-1 - n3)*d6^(1 - n6)*n3)/(d10*d4^n4*d5^n5*d7*d8*d9*(-1 + n6)) + 
        (d4^(-1 - n4)*d6^(1 - n6)*n4)/(d10*d3^n3*d5^n5*d8*d9*(-1 + n6)) - 
        (d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(d10*d5^n5*d7*d8*d9*(-1 + n6)) - 
        (d4^(-1 - n4)*d6^(1 - n6)*n4)/(d10*d3^n3*d5^n5*d7*d8*d9*(-1 + n6)) + 
        (2*d5^(1 - n5)*d6^(-1 - n6)*n6)/(3*d10*d3^n3*d4^n4*d7*d9*(-1 + n5)) - 
        (2*d5^(2 - n5)*d6^(-1 - n6)*n6)/(3*d10*d3^n3*d4^n4*d7*d8*d9*(-1 + n5)) + 
        (d5^(1 - n5)*rat(3 + 2*d - 6*n3 - 3*n5, 3*(-1 + n5)))/(d10*d3^n3*d4^n4*d6^n6*d7*d8*d9) + 
        (d6^(1 - n6)*rat(-1 + 2*d - 6*n3 - 3*n4 + n6, 3*(-1 + n6)))/(d10*d3^n3*d4^n4*d5^n5*d7*d8*d9);
        
* n1 == 0 && n10 == 1 && n2 == 0 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n5 != 1 && n5 != 2
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3^n3?pos_/d4^n4?pos_/d5^n5?{>2}/d6/d7/d8/d9/d10 =
        
        d5^(1 - n5)/(3*d10*d3^n3*d4^n4*d6*d7^2*d8*(-1 + n5)) + 
        d5^(1 - n5)/(3*d10*d3^n3*d4^n4*d6*d8*d9^2*(-1 + n5)) - 
        d5^(2 - n5)/(3*d10*d3^n3*d4^n4*d6*d7*d8*d9^2*(-1 + n5)) + 
        d5^(1 - n5)/(3*d10*d3^n3*d4^n4*d6^2*d7*d9*(-1 + n5)) + 
        d5^(1 - n5)/(3*d10^2*d3^n3*d4^n4*d6*d7*d9*(-1 + n5)) + 
        d5^(1 - n5)/(3*d10*d3^n3*d4^n4*d7*d8^2*d9*(-1 + n5)) - 
        d5^(2 - n5)/(3*d10*d3^n3*d4^n4*d6^2*d7*d8*d9*(-1 + n5)) - 
        (d3^(1 - n3)*d5^(1 - n5))/(3*d10^2*d4^n4*d6*d7*d8*d9*(-1 + n5)) - 
        (d3^(-1 - n3)*d5^(1 - n5)*n3)/(3*d10*d4^n4*d6*d7*d8*d9*(-1 + n5)) + 
        (d4^(-1 - n4)*d5^(1 - n5)*n4)/(3*d10*d3^n3*d6*d8*d9*(-1 + n5)) - 
        (d3^(1 - n3)*d4^(-1 - n4)*d5^(1 - n5)*n4)/(3*d10*d6*d7*d8*d9*(-1 + n5)) - 
        (d3^(-1 - n3)*d5^(2 - n5)*n3)/(3*d4^n4*d6*d7*d8^2*d9*(-2 + n5)*(-1 + n5)) + 
        (d3^(-1 - n3)*d5^(2 - n5)*n3)/(3*d10*d4^n4*d6*d7*d8^2*d9*(-2 + n5)*(-1 + n5)) - 
        (d3^(-1 - n3)*d4^(1 - n4)*d5^(2 - n5)*n3)/(3*d10*d6*d7^2*d8*d9*(-2 + n5)*(-1 + n5)) + 
        (d3^(-1 - n3)*d5^(2 - n5)*n3)/(3*d10*d4^n4*d6*d7^2*d8*d9*(-2 + n5)*(-1 + n5)) + 
        (d3^(-2 - n3)*d5^(2 - n5)*(2*n3 + 2*n3^2))/(3*d10*d4^n4*d6*d7*d8*d9*(-2 + n5)*(-1 + n5)) + 
        (d5^(2 - n5)*(2 + n3 - n5))/(3*d10*d3^n3*d4^n4*d6*d7*d8^2*d9*(-2 + n5)*(-1 + n5)) + 
        (d5^(2 - n5)*(2 + n3 - n5))/(3*d10*d3^n3*d4^n4*d6*d7^2*d8*d9*(-2 + n5)*(-1 + n5)) + 
        (d3^(-1 - n3)*d5^(2 - n5)*rat(4*n3 - d*n3 + 2*n3^2, 3*(-2 + n5)*(-1 + n5)))/(d10*d4^n4*d6*d7*d8*d9) + 
        (d5^(1 - n5)*rat(3 + d - n3 - 3*n5, 3*(-1 + n5)))/(d10*d3^n3*d4^n4*d6*d7*d8*d9);
        
* n1 == 0 && n10 == 1 && n2 == 0 && n5 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n3 != 1 && n6 != 1
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3^n3?{>1}/d4^n4?pos_/d5/d6^n6?{>1}/d7/d8/d9/d10 =
        
        1/(9*d10*d3^n3*d4^n4*d5*d6^n6*d7*d9) - 
        1/(3*d10*d3^n3*d4^n4*d6^n6*d7*d8*d9) + 
        2/(9*d3^n3*d4^n4*d5*d6^n6*d7*d8*d9) + 
        (4*d3^(1 - n3))/(9*d10^2*d4^n4*d5*d6^n6*d7*d9*(-1 + n3)) - 
        (4*d3^(2 - n3))/(9*d10^2*d4^n4*d5*d6^n6*d7*d8*d9*(-1 + n3)) - 
        (2*d3^(1 - n3)*d4^(-1 - n4)*n4)/(9*d10*d5*d6^n6*d7*d8*(-1 + n3)) - 
        (2*d3^(1 - n3)*d4^(-1 - n4)*n4)/(9*d10*d5*d6^n6*d7*d9*(-1 + n3)) + 
        (4*d3^(1 - n3)*d4^(-1 - n4)*n4)/(9*d10*d5*d6^n6*d8*d9*(-1 + n3)) + 
        (2*d3^(1 - n3)*d4^(-1 - n4)*n4)/(9*d10*d6^n6*d7*d8*d9*(-1 + n3)) + 
        (2*d3^(1 - n3)*d4^(-1 - n4)*n4)/(9*d5*d6^n6*d7*d8*d9*(-1 + n3)) - 
        (4*d3^(2 - n3)*d4^(-1 - n4)*n4)/(9*d10*d5*d6^n6*d7*d8*d9*(-1 + n3)) + 
        (2*d6^(1 - n6))/(9*d10*d3^n3*d4^n4*d5*d7^2*d8*(-1 + n6)) + 
        (2*d6^(1 - n6))/(9*d10*d3^n3*d4^n4*d5*d8*d9^2*(-1 + n6)) - 
        (2*d6^(1 - n6))/(9*d10*d3^n3*d4^n4*d7*d8*d9^2*(-1 + n6)) - 
        (2*d6^(1 - n6))/(9*d10*d3^n3*d4^n4*d7^2*d8*d9*(-1 + n6)) - 
        (2*d6^(1 - n6))/(3*d10*d3^n3*d4^n4*d5^2*d7*d8*d9*(-1 + n6)) + 
        (2*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d3^n3*d5*d7*d8*(-1 + n6)) - 
        (d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d3^n3*d5*d8*d9*(-1 + n6)) - 
        (2*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d3^n3*d7*d8*d9*(-1 + n6)) - 
        (d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d10*d3^n3*d5*d7*d8*d9*(-1 + n6)) + 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d5*d7^2*d8*(-1 + n3)*(-1 + n6)) - 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d5^2*d7*d8*(-1 + n3)*(-1 + n6)) - 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d5*d8*d9^2*(-1 + n3)*(-1 + n6)) + 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d7*d8*d9^2*(-1 + n3)*(-1 + n6)) + 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10^2*d5*d7*d9*(-1 + n3)*(-1 + n6)) - 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d7*d8^2*d9*(-1 + n3)*(-1 + n6)) + 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d6^(2 - n6)*n4)/(9*d10*d5*d7*d8^2*d9*(-1 + n3)*(-1 + n6)) + 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d5^2*d8*d9*(-1 + n3)*(-1 + n6)) - 
        (2*d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10*d7^2*d8*d9*(-1 + n3)*(-1 + n6)) - 
        (2*d3^(2 - n3)*d4^(-1 - n4)*d6^(1 - n6)*n4)/(9*d10^2*d5*d7*d8*d9*(-1 + n3)*(-1 + n6)) + 
        (d3^(2 - n3)*d4^(-2 - n4)*d6^(1 - n6)*(-2*n4 - 2*n4^2))/(9*d10*d5*d7*d8*d9*(-1 + n3)*(-1 + n6)) + 
        (d3^(1 - n3)*d4^(-2 - n4)*d6^(1 - n6)*(2*n4 + 2*n4^2))/(9*d10*d5*d8*d9*(-1 + n3)*(-1 + n6)) + 
        (d3^(1 - n3)*d4^(-1 - n4)*d6^(1 - n6)*(n4 - 3*n3*n4 + 2*n4*n6))/(9*d10*d5*d7*d8*d9*(-1 + n3)*(-1 + n6)) + 
        (d3^(1 - n3)*rat(6 + 2*d - 6*n3 - 2*n4, 9*(-1 + n3)))/(d10*d4^n4*d5*d6^n6*d7*d8*d9) + 
        (d6^(1 - n6)*rat(-3 + 2*d + n4 - 3*n6, 9*(-1 + n6)))/(d10*d3^n3*d4^n4*d5*d7*d8*d9);
        
* n1 == 0 && n10 == 1 && n4 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n2 != 0 && n5 != 1
        if(count(d1,1)==0) id,only,ifmatch->sortme  1/d2^n2?neg_/d3^n3?pos_/d4/d5^n5?{>1}/d6/d7/d8/d9/d10 =
        
        -d2^(-1 - n2)/(2*d10*d3^n3*d4*d5^n5*d6*d7*d9) + 
        d2^(-1 - n2)/(2*d10*d3^n3*d4*d5^n5*d7*d8*d9) - 
        d2^(-1 - n2)/(2*d10*d3^n3*d5^n5*d6*d7*d8*d9) + 
        (d1*d2^(-1 - n2))/(2*d10*d3^n3*d4*d5^n5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*d5^(1 - n5))/(d10*d3^n3*d4*d6*d7*d8^2*d9*(1 - n5)) + 
        d5^(1 - n5)/(2*d10*d2^n2*d3^n3*d4*d6*d8*d9^2*(-1 + n5)) - 
        (d2^(-1 - n2)*d5^(1 - n5))/(2*d3^n3*d4*d6*d7*d8*d9^2*(-1 + n5)) + 
        (d2^(-1 - n2)*d5^(1 - n5))/(2*d10*d3^n3*d4*d6*d7*d8*d9^2*(-1 + n5)) - 
        d5^(2 - n5)/(2*d10*d2^n2*d3^n3*d4*d6*d7*d8*d9^2*(-1 + n5)) - 
        (d2^(-1 - n2)*d5^(1 - n5))/(2*d10*d3^n3*d4*d6^2*d7*d9*(-1 + n5)) + 
        d5^(1 - n5)/(2*d10*d2^n2*d3^n3*d4*d6^2*d7*d9*(-1 + n5)) - 
        (d2^(-1 - n2)*d5^(1 - n5))/(2*d10^2*d3^n3*d4*d6*d7*d9*(-1 + n5)) - 
        (d2^(-1 - n2)*d5^(1 - n5))/(2*d10*d3^n3*d6^2*d7*d8*d9*(-1 + n5)) + 
        (d2^(-1 - n2)*d5^(2 - n5))/(2*d10*d3^n3*d4*d6^2*d7*d8*d9*(-1 + n5)) - 
        d5^(2 - n5)/(2*d10*d2^n2*d3^n3*d4*d6^2*d7*d8*d9*(-1 + n5)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*d5^(1 - n5))/(2*d10^2*d4*d6*d7*d8*d9*(-1 + n5)) - 
        (d2^(-1 - n2)*d5^(1 - n5))/(2*d10^2*d3^n3*d4*d6*d7*d8*d9*(-1 + n5)) + 
        (d2^(-2 - n2)*d5^(1 - n5)*(-1 - n2))/(2*d10*d3^n3*d4*d6*d7*d9*(-1 + n5)) + 
        (d2^(-1 - n2)*d5^(1 - n5)*(-1 + n2))/(2*d10*d3^n3*d4*d6*d7*d8*d9*(-1 + n5)) + 
        (d2^(-1 - n2)*d5^(1 - n5)*n2)/(2*d10*d3^n3*d4*d6*d7*d9*(-1 + n5)) - 
        (d2^(-1 - n2)*d5^(1 - n5)*n2)/(2*d10*d3^n3*d4*d7*d8*d9*(-1 + n5)) + 
        (d2^(-1 - n2)*d5^(1 - n5)*n2)/(2*d10*d3^n3*d6*d7*d8*d9*(-1 + n5)) - 
        (d1*d2^(-1 - n2)*d5^(1 - n5)*n2)/(2*d10*d3^n3*d4*d6*d7*d8*d9*(-1 + n5)) + 
        (d2^(-2 - n2)*d5^(1 - n5)*(1 + n2))/(2*d10*d3^n3*d4*d6*d7*d8*d9*(-1 + n5)) + 
        (d1*d2^(-2 - n2)*d5^(1 - n5)*(1 + n2))/(2*d10*d3^n3*d4*d6*d7*d8*d9*(-1 + n5)) + 
        (d5^(1 - n5)*rat(d - 2*n5, 2*(-1 + n5)))/(d10*d2^n2*d3^n3*d4*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n2 == 0 && n5 == 2 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n3 != 1 && n4 != 1
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3^n3?{>1}/d4^n4?{>1}/d5^2/d6/d7/d8/d9/d10 =
        
        1/(8*d10*d3^n3*d4^n4*d5^2*d6*d7*d8) + 
        1/(3*d10*d3^n3*d4^n4*d5^2*d6*d8*d9) - 
        1/(8*d10*d3^n3*d4^n4*d5*d6*d7*d8*d9) + 
        d3^(1 - n3)/(d10*d4^n4*d5^3*d6*d8*d9*(1 - n3)) + 
        d3^(1 - n3)/(8*d10*d4^n4*d5^2*d6*d7^2*d8*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d4^n4*d5^2*d6^2*d7*d8*(-1 + n3)) + 
        d3^(1 - n3)/(12*d10*d4^n4*d5^3*d6*d7*d8*(-1 + n3)) + 
        d3^(1 - n3)/(3*d10^2*d4^n4*d5^2*d6*d7*d9*(-1 + n3)) + 
        d3^(1 - n3)/(8*d10*d4^n4*d5^2*d7*d8^2*d9*(-1 + n3)) + 
        (5*d3^(1 - n3))/(24*d4^n4*d5^2*d6*d7*d8^2*d9*(-1 + n3)) - 
        (5*d3^(2 - n3))/(24*d10*d4^n4*d5^2*d6*d7*d8^2*d9*(-1 + n3)) - 
        d3^(1 - n3)/(8*d10*d4^n4*d5*d6*d7*d8^2*d9*(-1 + n3)) + 
        d3^(1 - n3)/(24*d10*d4^n4*d5^2*d6^2*d8*d9*(-1 + n3)) - 
        (5*d3^(2 - n3))/(24*d10*d4^n4*d5^2*d6*d7^2*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(8*d10*d4^n4*d5*d6*d7^2*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d4^n4*d5^2*d6^2*d7*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(24*d10*d4^n4*d5^2*d6^2*d7*d8*d9*(-1 + n3)) + 
        (11*d3^(2 - n3))/(12*d10*d4^n4*d5^3*d6*d7*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(3*d10^2*d4^n4*d5^2*d6*d7*d8*d9*(-1 + n3)) + 
        d4^(1 - n4)/(8*d10*d3^n3*d5^2*d6*d7^2*d8*(-1 + n4)) + 
        d4^(1 - n4)/(2*d10*d3^n3*d5^2*d6*d8*d9^2*(-1 + n4)) - 
        d4^(1 - n4)/(2*d10*d3^n3*d5*d6*d7*d8*d9^2*(-1 + n4)) + 
        (5*d4^(1 - n4))/(12*d10*d3^n3*d5^2*d6^2*d7*d9*(-1 + n4)) + 
        d4^(1 - n4)/(3*d10^2*d3^n3*d5^2*d6*d7*d9*(-1 + n4)) + 
        d4^(1 - n4)/(3*d3^n3*d5^2*d6*d7*d8^2*d9*(-1 + n4)) + 
        d4^(1 - n4)/(24*d10*d3^n3*d5^2*d6^2*d8*d9*(-1 + n4)) + 
        (5*d4^(2 - n4))/(24*d10*d3^n3*d5^2*d6*d7^2*d8*d9*(-1 + n4)) - 
        d4^(1 - n4)/(8*d10*d3^n3*d5*d6*d7^2*d8*d9*(-1 + n4)) + 
        (5*d4^(1 - n4))/(24*d3^n3*d5^2*d6^2*d7*d8*d9*(-1 + n4)) - 
        d4^(2 - n4)/(24*d10*d3^n3*d5^2*d6^2*d7*d8*d9*(-1 + n4)) - 
        (5*d4^(1 - n4))/(8*d10*d3^n3*d5*d6^2*d7*d8*d9*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(3*d10^2*d5^2*d6*d7*d8*d9*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5^2*d6^2*d7*d8^2*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5^3*d6*d7*d8^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5^2*d6^2*d7^2*d8*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(6*d10*d5^3*d6*d7^2*d8*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5^3*d6^2*d7*d8*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(12*d5^2*d6*d7*d8^2*d9^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(12*d10*d5^2*d6*d7*d8^2*d9^2*(-1 + n3)*(-1 + n4)) - 
        (7*d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5^2*d6^2*d8*d9^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(12*d10*d5^2*d6*d7^2*d8*d9^2*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(2 - n4))/(12*d10*d5^2*d6*d7^2*d8*d9^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5*d6*d7^2*d8*d9^2*(-1 + n3)*(-1 + n4)) + 
        (7*d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5*d6^2*d7*d8*d9^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(6*d10*d5^3*d6*d7*d8*d9^2*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(2 - n4))/(6*d10*d5^3*d6*d7*d8*d9^2*(-1 + n3)*(-1 + n4)) - 
        (7*d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5^2*d6^3*d7*d9*(-1 + n3)*(-1 + n4)) + 
        (2*d3^(1 - n3)*d4^(1 - n4))/(3*d10*d5^3*d6^2*d7*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10^2*d5^2*d6^2*d7*d9*(-1 + n3)*(-1 + n4)) - 
        (5*d3^(1 - n3)*d4^(1 - n4))/(6*d10^2*d5^3*d6*d7*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5^2*d6^2*d8^2*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5^3*d6*d8^2*d9*(-1 + n3)*(-1 + n4)) - 
        (5*d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5^3*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(4*d5^2*d6^2*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(4*d10*d5^2*d6^2*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(8*d10*d5*d6^2*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(2*d5^3*d6*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(2*d10*d5^3*d6*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5^3*d6^2*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(4*d10*d5^2*d6^2*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(2 - n4))/(4*d10*d5^2*d6^2*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5*d6^2*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(6*d10*d5^3*d6*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(2 - n4))/(6*d10*d5^3*d6*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (7*d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5*d6^3*d7*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(2*d5^3*d6^2*d7*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(2*d10*d5^3*d6^2*d7*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(2 - n3)*d4^(1 - n4))/(24*d10^2*d5^2*d6^2*d7*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (5*d3^(2 - n3)*d4^(1 - n4))/(6*d10^2*d5^3*d6*d7*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(-1 - n3)*d4^(1 - n4)*n3)/(d10*d5^2*d6*d7*d8*d9*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(-1 - n4)*n4)/(3*d10*d5^2*d6*d8*d9*(-1 + n3)) - 
        (d3^(2 - n3)*d4^(-1 - n4)*n4)/(3*d10*d5^2*d6*d7*d8*d9*(-1 + n3)) + 
        (d3^(1 - n3)*d4^(1 - n4)*(2 - 5*n3 + 5*n4))/(24*d10*d5^2*d6*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*rat(19 + 4*d - 29*n3, 24*(-1 + n3)))/(d10*d4^n4*d5^2*d6*d7*d8*d9) + 
        (d3^(1 - n3)*d4^(1 - n4)*rat(16 + d - 8*n3, 24*(-1 + n3)*(-1 + n4)))/(d10*d5^2*d6*d7*d8^2*d9) + 
        (d3^(1 - n3)*d4^(1 - n4)*rat(1 - d + 2*n3, 12*(-1 + n3)*(-1 + n4)))/(d10*d5^2*d6*d7*d8*d9^2) + 
        (d3^(1 - n3)*d4^(1 - n4)*rat(6 - 2*d + 7*n3 - 4*n4, 24*(-1 + n3)*(-1 + n4)))/(d10*d5^2*d6^2*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-52 + 15*d - 24*n3 - 3*n4, 24*(-1 + n4)))/(d10*d3^n3*d5^2*d6*d7*d8*d9) + 
        (d3^(1 - n3)*d4^(1 - n4)*rat(-42 - 4*d + 22*n3 + 7*n4, 12*(-1 + n3)*(-1 + n4)))/(d10*d5^3*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n2 == 0 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n3 != 1 && n4 != 1
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3^n3?{>1}/d4^n4?{>1}/d5/d6/d7/d8/d9/d10 =
        
        1/(8*d10*d3^n3*d4^n4*d5*d6*d7*d8) + 
        1/(3*d10*d3^n3*d4^n4*d5*d6*d8*d9) - 
        1/(8*d10*d3^n3*d4^n4*d6*d7*d8*d9) + 
        d3^(1 - n3)/(8*d10*d4^n4*d5*d6*d7^2*d8*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d4^n4*d5*d6^2*d7*d8*(-1 + n3)) + 
        d3^(1 - n3)/(24*d10*d4^n4*d5^2*d6*d7*d8*(-1 + n3)) + 
        d3^(1 - n3)/(3*d10^2*d4^n4*d5*d6*d7*d9*(-1 + n3)) + 
        d3^(1 - n3)/(8*d10*d4^n4*d5*d7*d8^2*d9*(-1 + n3)) - 
        d3^(1 - n3)/(8*d10*d4^n4*d6*d7*d8^2*d9*(-1 + n3)) + 
        (5*d3^(1 - n3))/(24*d4^n4*d5*d6*d7*d8^2*d9*(-1 + n3)) - 
        (5*d3^(2 - n3))/(24*d10*d4^n4*d5*d6*d7*d8^2*d9*(-1 + n3)) + 
        d3^(1 - n3)/(24*d10*d4^n4*d5*d6^2*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(2*d10*d4^n4*d5^2*d6*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(8*d10*d4^n4*d6*d7^2*d8*d9*(-1 + n3)) - 
        (5*d3^(2 - n3))/(24*d10*d4^n4*d5*d6*d7^2*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d4^n4*d5*d6^2*d7*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(24*d10*d4^n4*d5*d6^2*d7*d8*d9*(-1 + n3)) + 
        (11*d3^(2 - n3))/(24*d10*d4^n4*d5^2*d6*d7*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(3*d10^2*d4^n4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        d4^(1 - n4)/(8*d10*d3^n3*d5*d6*d7^2*d8*(-1 + n4)) + 
        d4^(1 - n4)/(2*d10*d3^n3*d5*d6*d8*d9^2*(-1 + n4)) - 
        d4^(1 - n4)/(2*d10*d3^n3*d6*d7*d8*d9^2*(-1 + n4)) + 
        (5*d4^(1 - n4))/(12*d10*d3^n3*d5*d6^2*d7*d9*(-1 + n4)) + 
        d4^(1 - n4)/(3*d10^2*d3^n3*d5*d6*d7*d9*(-1 + n4)) + 
        d4^(1 - n4)/(3*d3^n3*d5*d6*d7*d8^2*d9*(-1 + n4)) + 
        d4^(1 - n4)/(24*d10*d3^n3*d5*d6^2*d8*d9*(-1 + n4)) - 
        d4^(1 - n4)/(8*d10*d3^n3*d6*d7^2*d8*d9*(-1 + n4)) + 
        (5*d4^(2 - n4))/(24*d10*d3^n3*d5*d6*d7^2*d8*d9*(-1 + n4)) - 
        (5*d4^(1 - n4))/(8*d10*d3^n3*d6^2*d7*d8*d9*(-1 + n4)) + 
        (5*d4^(1 - n4))/(24*d3^n3*d5*d6^2*d7*d8*d9*(-1 + n4)) - 
        d4^(2 - n4)/(24*d10*d3^n3*d5*d6^2*d7*d8*d9*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(3*d10^2*d5*d6*d7*d8*d9*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5*d6^2*d7*d8^2*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5^2*d6*d7*d8^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5*d6^2*d7^2*d8*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5^2*d6*d7^2*d8*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5^2*d6^2*d7*d8*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(12*d5*d6*d7*d8^2*d9^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(12*d10*d5*d6*d7*d8^2*d9^2*(-1 + n3)*(-1 + n4)) - 
        (7*d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5*d6^2*d8*d9^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(12*d10*d6*d7^2*d8*d9^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(12*d10*d5*d6*d7^2*d8*d9^2*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(2 - n4))/(12*d10*d5*d6*d7^2*d8*d9^2*(-1 + n3)*(-1 + n4)) + 
        (7*d3^(1 - n3)*d4^(1 - n4))/(24*d10*d6^2*d7*d8*d9^2*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(12*d10*d5^2*d6*d7*d8*d9^2*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(2 - n4))/(12*d10*d5^2*d6*d7*d8*d9^2*(-1 + n3)*(-1 + n4)) - 
        (7*d3^(1 - n3)*d4^(1 - n4))/(12*d10*d5*d6^3*d7*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(3*d10*d5^2*d6^2*d7*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10^2*d5*d6^2*d7*d9*(-1 + n3)*(-1 + n4)) - 
        (5*d3^(1 - n3)*d4^(1 - n4))/(12*d10^2*d5^2*d6*d7*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5*d6^2*d8^2*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5^2*d6*d8^2*d9*(-1 + n3)*(-1 + n4)) - 
        (5*d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5^2*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(8*d10*d6^2*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(4*d5*d6^2*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(4*d10*d5*d6^2*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(4*d5^2*d6*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(4*d10*d5^2*d6*d7*d8^2*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d5^2*d6^2*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(24*d10*d6^2*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(4*d10*d5*d6^2*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(2 - n4))/(4*d10*d5*d6^2*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(12*d10*d5^2*d6*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(2 - n4))/(12*d10*d5^2*d6*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (7*d3^(1 - n3)*d4^(1 - n4))/(12*d10*d6^3*d7*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(1 - n3)*d4^(1 - n4))/(4*d5^2*d6^2*d7*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(2 - n3)*d4^(1 - n4))/(4*d10*d5^2*d6^2*d7*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(2 - n3)*d4^(1 - n4))/(24*d10^2*d5*d6^2*d7*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (5*d3^(2 - n3)*d4^(1 - n4))/(12*d10^2*d5^2*d6*d7*d8*d9*(-1 + n3)*(-1 + n4)) - 
        (d3^(-1 - n3)*d4^(1 - n4)*n3)/(d10*d5*d6*d7*d8*d9*(-1 + n4)) + 
        (d3^(1 - n3)*d4^(-1 - n4)*n4)/(3*d10*d5*d6*d8*d9*(-1 + n3)) - 
        (d3^(2 - n3)*d4^(-1 - n4)*n4)/(3*d10*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d3^(1 - n3)*d4^(1 - n4)*(-5*n3 + 5*n4))/(24*d10*d5*d6*d7^2*d8*d9*(-1 + n3)*(-1 + n4)) + 
        (d3^(1 - n3)*rat(26 + 4*d - 29*n3, 24*(-1 + n3)))/(d10*d4^n4*d5*d6*d7*d8*d9) + 
        (d3^(1 - n3)*d4^(1 - n4)*rat(12 + d - 8*n3, 24*(-1 + n3)*(-1 + n4)))/(d10*d5*d6*d7*d8^2*d9) + 
        (d3^(1 - n3)*d4^(1 - n4)*rat(-d + 2*n3, 12*(-1 + n3)*(-1 + n4)))/(d10*d5*d6*d7*d8*d9^2) + 
        (d3^(1 - n3)*d4^(1 - n4)*rat(2 - 2*d + 7*n3 - 4*n4, 24*(-1 + n3)*(-1 + n4)))/(d10*d5*d6^2*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-10 + 5*d - 8*n3 - n4, 8*(-1 + n4)))/(d10*d3^n3*d5*d6*d7*d8*d9) + 
        (d3^(1 - n3)*d4^(1 - n4)*rat(-42 - 4*d + 22*n3 + 7*n4, 24*(-1 + n3)*(-1 + n4)))/(d10*d5^2*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n4 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n2 != 0 && n3 != 1
        if(count(d1,1)==0) id,only,ifmatch->sortme 1/d2^n2?neg_/d3^n3?{>1}/d4/d5/d6/d7/d8/d9/d10 =
        
        -d2^(-1 - n2)/(2*d10*d3^n3*d4*d5*d6*d7*d8) - 
        d2^(-1 - n2)/(2*d10*d3^n3*d4*d5*d6*d7*d9) + 
        d2^(-1 - n2)/(2*d3^n3*d4*d5*d6*d7*d8*d9) + 
        (d1*d2^(-1 - n2))/(2*d10*d3^n3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(d10*d4*d5*d6*d7*d8^2*d9*(1 - n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10^2*d4*d5*d6*d7*d8*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4*d5*d6^2*d7*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10^2*d4*d5*d6*d7*d9*(-1 + n3)) + 
        d3^(1 - n3)/(2*d10^2*d2^n2*d4*d5*d6*d7*d9*(-1 + n3)) + 
        d3^(1 - n3)/(2*d10*d2^n2*d4^2*d5*d6*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4^2*d5*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4*d6^2*d7*d8*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4*d5*d6^2*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3))/(2*d10*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(2*d10*d2^n2*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(2 - n3))/(2*d10^2*d4*d5*d6*d7*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(2*d10^2*d2^n2*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(-1 - n2))/(2*d10*d4*d5*d6*d7*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*(-1 + n2))/(2*d10*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5*d6*d7*d8*(-1 + n3)) + 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5*d6*d7*d9*(-1 + n3)) - 
        (d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d4*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (d1*d2^(-1 - n2)*d3^(1 - n3)*n2)/(2*d10*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d2^(-2 - n2)*d3^(1 - n3)*(1 + n2))/(2*d10*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d1*d2^(-2 - n2)*d3^(1 - n3)*(1 + n2))/(2*d10*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d3^(1 - n3)*rat(d - 2*n3, 2*(-1 + n3)))/(d10*d2^n2*d4*d5*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n2 == 0 && n3 == 1 && n5 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n6 != 1 && n6 != 2
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3/d4^n4?pos_/d5/d6^n6?{>2}/d7/d8/d9/d10 =

        5/(3*d10*d3*d4^n4*d5*d6^n6*d7*d9) - 
        5/(3*d10*d3*d4^n4*d6^n6*d7*d8*d9) + 
        (2*d6^(1 - n6))/(3*d10*d3*d4^n4*d5*d7^2*d8*(-1 + n6)) + 
        (4*d6^(1 - n6))/(3*d10*d3*d4^n4*d5*d8*d9^2*(-1 + n6)) - 
        (4*d6^(1 - n6))/(3*d10*d3*d4^n4*d7*d8*d9^2*(-1 + n6)) - 
        (2*d6^(1 - n6))/(3*d10*d3*d4^n4*d5^2*d7*d9*(-1 + n6)) + 
        d6^(1 - n6)/(d10*d3^2*d4^n4*d5*d7*d9*(-1 + n6)) - 
        (4*d6^(1 - n6))/(3*d10^2*d3*d4^n4*d5*d7*d9*(-1 + n6)) - 
        (2*d6^(1 - n6))/(3*d10*d4^n4*d5*d7*d8^2*d9*(-1 + n6)) + 
        (2*d6^(1 - n6))/(3*d3*d4^n4*d5*d7*d8^2*d9*(-1 + n6)) - 
        (2*d6^(1 - n6))/(3*d10*d3*d4^n4*d7^2*d8*d9*(-1 + n6)) + 
        (2*d4^(1 - n4)*d6^(1 - n6))/(3*d10*d3*d5*d7^2*d8*d9*(-1 + n6)) - 
        (2*d6^(1 - n6))/(3*d10*d4^n4*d5*d7^2*d8*d9*(-1 + n6)) - 
        (2*d6^(1 - n6))/(3*d10*d4^n4*d5^2*d7*d8*d9*(-1 + n6)) + 
        (2*d6^(1 - n6))/(3*d3*d4^n4*d5^2*d7*d8*d9*(-1 + n6)) - 
        (10*d6^(1 - n6))/(3*d10*d3*d4^n4*d5^2*d7*d8*d9*(-1 + n6)) + 
        (4*d6^(1 - n6))/(3*d10^2*d4^n4*d5*d7*d8*d9*(-1 + n6)) - 
        d6^(1 - n6)/(d3^2*d4^n4*d5*d7*d8*d9*(-1 + n6)) + 
        d6^(1 - n6)/(d10*d3^2*d4^n4*d5*d7*d8*d9*(-1 + n6)) + 
        (5*d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d10*d3*d5*d7*d8*(-1 + n6)) + 
        (2*d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d10*d3*d5*d7*d9*(-1 + n6)) - 
        (5*d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d10*d3*d5*d8*d9*(-1 + n6)) - 
        (2*d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d10*d3*d7*d8*d9*(-1 + n6)) + 
        (5*d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d10*d5*d7*d8*d9*(-1 + n6)) - 
        (5*d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d3*d5*d7*d8*d9*(-1 + n6)) + 
        (7*d4^(-1 - n4)*d6^(1 - n6)*n4)/(3*d10*d3*d5*d7*d8*d9*(-1 + n6)) - 
        (2*d6^(2 - n6))/(3*d10^2*d3*d4^n4*d5^2*d7*d9*(-2 + n6)*(-1 + n6)) - 
        (2*d6^(2 - n6))/(3*d10*d4^n4*d5^2*d7*d8^2*d9*(-2 + n6)*(-1 + n6)) + 
        (2*d6^(2 - n6))/(3*d3*d4^n4*d5^2*d7*d8^2*d9*(-2 + n6)*(-1 + n6)) - 
        (2*d6^(2 - n6))/(3*d10*d3*d4^n4*d5^2*d7*d8^2*d9*(-2 + n6)*(-1 + n6)) + 
        (2*d6^(2 - n6))/(3*d10^2*d4^n4*d5^2*d7*d8*d9*(-2 + n6)*(-1 + n6)) + 
        (2*d6^(2 - n6))/(3*d10^2*d3*d4^n4*d5^2*d7*d8*d9*(-2 + n6)*(-1 + n6)) + 
        (d4^(-1 - n4)*d6^(2 - n6)*n4)/(d10*d3*d5^2*d7*d8*(-2 + n6)*(-1 + n6)) + 
        (2*d4^(-1 - n4)*d6^(2 - n6)*n4)/(d10*d3*d5*d7*d8*d9^2*(-2 + n6)*(-1 + n6)) + 
        (2*d4^(-1 - n4)*d6^(2 - n6)*n4)/(3*d10^2*d3*d5*d7*d9*(-2 + n6)*(-1 + n6)) + 
        (2*d4^(-1 - n4)*d6^(2 - n6)*n4)/(3*d10*d5*d7*d8^2*d9*(-2 + n6)*(-1 + n6)) - 
        (2*d4^(-1 - n4)*d6^(2 - n6)*n4)/(3*d3*d5*d7*d8^2*d9*(-2 + n6)*(-1 + n6)) + 
        (2*d4^(-1 - n4)*d6^(2 - n6)*n4)/(3*d10*d3*d5*d7*d8^2*d9*(-2 + n6)*(-1 + n6)) - 
        (d4^(-1 - n4)*d6^(2 - n6)*n4)/(d10*d3*d5^2*d8*d9*(-2 + n6)*(-1 + n6)) + 
        (d4^(-1 - n4)*d6^(2 - n6)*n4)/(d10*d3*d5^2*d7*d8*d9*(-2 + n6)*(-1 + n6)) - 
        (2*d4^(-1 - n4)*d6^(2 - n6)*n4)/(3*d10^2*d5*d7*d8*d9*(-2 + n6)*(-1 + n6)) - 
        (2*d4^(-1 - n4)*d6^(2 - n6)*n4)/(3*d10^2*d3*d5*d7*d8*d9*(-2 + n6)*(-1 + n6)) + 
        (d6^(1 - n6)*rat(-2 + d + 4*n4 - 3*n6, 3*(-1 + n6)))/(d10*d3*d4^n4*d5*d7*d8*d9) + 
        (d4^(-1 - n4)*d6^(2 - n6)*rat(n4 - d*n4 + n4*n6, (-2 + n6)*(-1 + n6)))/(d10*d3*d5*d7*d8*d9);

* n10 == 1 && n3 == 1 && n4 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n1 != 0 && n2 != 0
        id,only,ifmatch->sortme 1/d1^n1?neg_/d2^n2?neg_/d3/d4/d5/d6/d7/d8/d9/d10 =
        
        (d2^(-1 - n2)*rat(-3, -4 + d - n1))/(d1^n1*d10*d3*d4*d5*d6*d7*d8*d9^2) + 
        (d2^(-1 - n2)*rat(-3, -4 + d - n1))/(d1^n1*d10^2*d3*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-2, -4 + d - n1))/(d10^2*d2^n2*d3*d4*d5*d6*d7*d9) + 
        (d1^(-1 - n1)*rat(-2, -4 + d - n1))/(d10*d2^n2*d4*d5*d6*d7*d8^2*d9) + 
        (d1^(-1 - n1)*rat(-2, -4 + d - n1))/(d10*d2^n2*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d1^(-1 - n1)*rat(-2, -4 + d - n1))/(d10*d2^n2*d3*d4*d5^2*d6*d8*d9) + 
        (d2^(-1 - n2)*rat(-2, -4 + d - n1))/(d1^n1*d10*d3*d4*d5^2*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-2, -4 + d - n1))/(d1^n1*d10*d3^2*d4*d5*d6*d7*d8*d9) + 
        rat(1, -4 + d - n1)/(d1^n1*d10*d2^n2*d3*d4*d6*d7*d8*d9^2) + 
        (d2^(-1 - n2)*rat(1, -4 + d - n1))/(d1^n1*d10*d3*d4*d5*d6^2*d7*d9) + 
        (d2^(-1 - n2)*rat(1, -4 + d - n1))/(d1^n1*d10^2*d3*d4*d5*d6*d7*d9) + 
        (d2^(-1 - n2)*rat(1, -4 + d - n1))/(d1^n1*d10*d4*d5*d6*d7*d8^2*d9) + 
        (d2^(-1 - n2)*rat(1, -4 + d - n1))/(d1^n1*d10*d3*d4*d5^2*d6*d8*d9) + 
        (d2^(-1 - n2)*rat(1, -4 + d - n1))/(d1^n1*d10*d3*d5*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(1, -4 + d - n1))/(d1^n1*d10*d3*d4*d5^2*d7*d8*d9) + 
        rat(1, -4 + d - n1)/(d1^n1*d10*d2^n2*d3*d4*d6^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(1, -4 + d - n1))/(d1^n1*d10*d3*d5*d6^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(1, 4 - d + n1))/(d1^n1*d10*d3*d4*d5^2*d6*d7*d8) + 
        rat(1, 4 - d + n1)/(d1^n1*d10*d2^n2*d3*d4*d5*d6*d8*d9^2) + 
        (d2^(-1 - n2)*rat(1, 4 - d + n1))/(d1^n1*d3*d4*d5*d6*d7*d8*d9^2) + 
        rat(1, 4 - d + n1)/(d1^n1*d10*d2^n2*d3*d4*d5*d6^2*d7*d9) + 
        (d2^(-1 - n2)*rat(1, 4 - d + n1))/(d1^n1*d10*d3*d4*d5^2*d6*d7*d9) + 
        (d2^(-1 - n2)*rat(1, 4 - d + n1))/(d1^n1*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d2^(-1 - n2)*rat(1, 4 - d + n1))/(d1^n1*d10*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d2^(-1 - n2)*rat(1, 4 - d + n1))/(d1^n1*d10*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(1, 4 - d + n1))/(d1^n1*d10*d3*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(1, 4 - d + n1))/(d1^n1*d10*d3*d4*d6^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(1, 4 - d + n1))/(d1^n1*d10^2*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(2, -4 + d - n1))/(d10*d2^n2*d3*d4*d5^2*d6*d7*d8) + 
        (d1^(-1 - n1)*rat(2, -4 + d - n1))/(d10*d2^n2*d4*d5*d6*d7*d8*d9^2) + 
        (d1^(-1 - n1)*rat(2, -4 + d - n1))/(d10*d2^n2*d3*d4*d5*d6*d7*d8*d9^2) + 
        (d1^(-1 - n1)*rat(2, -4 + d - n1))/(d2^n2*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d1^(-1 - n1)*rat(2, -4 + d - n1))/(d10*d2^n2*d3*d5^2*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(2, -4 + d - n1))/(d10^2*d2^n2*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(2, -4 + d - n1))/(d10^2*d2^n2*d3*d4*d5*d6*d7*d8*d9) + 
        (d1^(-2 - n1)*rat(-2 - 2*n1, 4 - d + n1))/(d10*d2^n2*d3*d4*d5*d6*d7*d8) + 
        (d1^(-1 - n1)*rat(-n1, 4 - d + n1))/(d10*d2^n2*d3*d4*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*rat(-n1, 4 - d + n1))/(d10*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*rat(n1, 4 - d + n1))/(d10*d3*d4*d5*d6*d7*d8) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*rat(n1, 4 - d + n1))/(d10*d3*d4*d5*d6*d7*d9) + 
        (d1^(-1 - n1)*rat(n1, 4 - d + n1))/(d10*d2^n2*d3*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*d2^(-1 - n2)*rat(2*n1, 4 - d + n1))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d1^(-2 - n1)*rat(2 + 2*n1, 4 - d + n1))/(d10*d2^n2*d4*d5*d6*d7*d8*d9) + 
        (d1^(-2 - n1)*rat(2 + 2*n1, 4 - d + n1))/(d10*d2^n2*d3*d4*d5*d6*d7*d8*d9) + 
        (d1^(1 - n1)*d2^(-2 - n2)*rat(-1 - n2, -4 + d - n1))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-2 - n2)*rat(-1 - n2, -4 + d - n1))/(d1^n1*d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-13 + 3*d - 2*n1 - n2, -4 + d - n1))/(d1^n1*d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-n2, -4 + d - n1))/(d1^n1*d10*d3*d4*d5*d6*d7*d9) + 
        (d2^(-1 - n2)*rat(-n2, -4 + d - n1))/(d1^n1*d10*d3*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(n2, -4 + d - n1))/(d1^n1*d10*d3*d4*d5*d7*d8*d9) + 
        (d1^(1 - n1)*d2^(-1 - n2)*rat(n2, -4 + d - n1))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-2 - n2)*rat(1 + n2, -4 + d - n1))/(d1^n1*d10*d3*d4*d5*d6*d7*d9);

* n1 == 0 && n10 == 1 && n2 == 0 && n4 == 1 && n5 == 2 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n3 != 1
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3^n3?{>1}/d4/d5^2/d6/d7/d8/d9/d10 =

        2/(3*d10*d3^n3*d4*d5^2*d6*d7*d8) + 
        1/(3*d10*d3^n3*d4*d5*d6*d8*d9^2) - 
        1/(3*d10*d3^n3*d4*d6*d7*d8*d9^2) + 
        1/(3*d10*d3^n3*d4*d5*d6^2*d7*d9) + 
        1/(3*d10^2*d3^n3*d4*d5*d6*d7*d9) + 
1/(3*d3^n3*d4*d5*d6*d7*d8^2*d9) - 2/(3*d10*d3^n3*d4*d5^2*d6*d8*d9) + 
1/(3*d10*d3^n3*d4^2*d5*d6*d8*d9) + 
1/(3*d10*d3^n3*d5*d6*d7^2*d8*d9) - 1/(3*d10*d3^n3*d4*d6^2*d7*d8*d9) + 2/(3*d10*d3^n3*d5^2*d6*d7*d8*d9) - d3^(2 - n3)/(3*d10*d4*d5^2*d6*d7^2*d8*(-1 + n3)) + d3^(1 - n3)/(3*d10^2*d4*d5^2*d6*d7*d8*(-1 + n3)) + d3^(2 - n3)/(3*d10^2*d4*d5^2*d6*d7*d8*(-1 + n3)) - d3^(2 - n3)/(3*d10*d4*d5^2*d6*d8*d9^2*(-1 + n3)) - d3^(2 - n3)/(3*d10*d4*d5^2*d6*d7*d8*d9^2*(-1 + n3)) + d3^(2 - n3)/(3*d10*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        (2*d3^(1 - n3))/(3*d10^2*d4*d5^2*d6*d7*d9*(-1 + n3)) - 
        (d1*d3^(1 - n3))/(3*d10^2*d4*d5^2*d6*d7*d9*(-1 + n3)) + 
        (2*d3^(2 - n3))/(3*d10^2*d4*d5^2*d6*d7*d9*(-1 + n3)) - d3^(2 - n3)/(3*d10*d4*d5^2*d7*d8^2*d9*(-1 + n3)) + d3^(2 - n3)/(3*d10*d4*d5^2*d6*d7*d8^2*d9*(-1 + n3)) + d3^(1 - n3)/(3*d10*d4^2*d5^2*d6*d8*d9*(-1 + n3)) - 
        (d1*d3^(1 - n3))/(3*d10*d4^2*d5^2*d6*d8*d9*(-1 + n3)) - d3^(1 - n3)/(3*d10^2*d4*d5^2*d6*d8*d9*(-1 + n3)) - d3^(2 - n3)/(3*d10^2*d4*d5^2*d6*d8*d9*(-1 + n3)) - d3^(2 - n3)/(3*d10*d4*d5^2*d6*d7^2*d8*d9*(-1 + n3)) - d3^(1 - n3)/(3*d10^2*d4*d5^2*d7*d8*d9*(-1 + n3)) - d3^(2 - n3)/(3*d10^2*d4*d5^2*d7*d8*d9*(-1 + n3)) + 
        (4*d3^(2 - n3))/(3*d10*d4*d5^3*d6*d7*d8*d9*(-1 + n3)) + d3^(1 - n3)/(3*d10^2*d5^2*d6*d7*d8*d9*(-1 + n3)) + d3^(2 - n3)/(3*d10^2*d5^2*d6*d7*d8*d9*(-1 + n3)) - 
        (2*d3^(2 - n3))/(3*d10*d4^2*d5^2*d6*d7*d8*d9*(-1 + n3)) - d3^(2 - n3)/(3*d10^2*d4*d5^2*d6*d7*d8*d9*(-1 + n3)) - d3^(3 - n3)/(3*d10^2*d4*d5^2*d6*d7*d8*d9*(-1 + n3)) + d3^(2 - n3)/(3*d10*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) + d3^(2 - n3)/(3*d10^2*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d3^(1 - n3)*(2 - n3))/(3*d10*d4^2*d5*d6*d7*d8*d9*(-1 + n3)) + 
        (d3^(1 - n3)*(2 - n3))/(3*d10^2*d4*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (d3^(-1 - n3)*n3)/(d10*d4*d5*d6*d7*d8*d9) + 
        rat(1 + 3*d - 7*n3, 3)/(d10*d3^n3*d4*d5*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(4 + d - 4*n3, 3*(-1 + n3)))/(d10*d4*d5^2*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(2 + d - 3*n3, 3*(-1 + n3)))/(d10*d4*d5*d6*d7*d8^2*d9) + 
        (d3^(1 - n3)*rat(2 + d - 3*n3, 3*(-1 + n3)))/(d10*d4*d5*d6*d7^2*d8*d9) + 
        (d3^(1 - n3)*rat(1 + d - 2*n3, 3*(-1 + n3)))/(d10*d4*d5^2*d6*d8*d9) + 
        (d3^(2 - n3)*rat(2 + d - 2*n3, 3*(-1 + n3)))/(d10*d4*d5*d6*d7*d8^2*d9) + 
        (d3^(2 - n3)*rat(2 + d - 2*n3, 3*(-1 + n3)))/(d10*d4*d5*d6*d7^2*d8*d9) + 
        (d3^(2 - n3)*rat(4 + d - 2*n3, 3*(-1 + n3)))/(d10*d4*d5^2*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(-1 - d + 2*n3, 3*(-1 + n3)))/(d4*d5*d6*d7*d8^2*d9) + 
        (d3^(1 - n3)*rat(-1 - d + 2*n3, 3*(-1 + n3)))/(d10*d5*d6*d7^2*d8*d9) + 
        (d3^(1 - n3)*rat(-1 - d + 2*n3, 3*(-1 + n3)))/(d10*d5^2*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(-d + 2*n3, 3*(-1 + n3)))/(d10*d4*d5^2*d6*d7*d8) + 
        (d3^(1 - n3)*rat(-d - d^2 + 2*n3 + 4*d*n3 - 4*n3^2, 3*(-1 + n3)))/(d10*d4*d5*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n2 == 0 && n4 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n3 != 1 && n3 != 2
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3^n3?{>2}/d4/d5/d6/d7/d8/d9/d10 =
        
        1/(3*d10*d3^n3*d4*d5*d6*d7*d8) - 
        1/(3*d10*d3^n3*d4*d5*d6*d8*d9) + 
        1/(3*d10*d3^n3*d5*d6*d7*d8*d9) + 
        d3^(1 - n3)/(6*d10*d4*d5*d6*d7^2*d8*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10^2*d4*d5*d6*d7*d8*(-1 + n3)) + 
        d3^(2 - n3)/(6*d10^2*d4*d5*d6*d7*d8*(-1 + n3)) - 
        d3^(1 - n3)/(3*d10*d4*d5*d6*d8*d9^2*(-1 + n3)) - 
        d3^(2 - n3)/(6*d10*d4*d5*d6*d8*d9^2*(-1 + n3)) + 
        d3^(1 - n3)/(3*d10*d4*d6*d7*d8*d9^2*(-1 + n3)) + 
        d3^(2 - n3)/(6*d10*d4*d6*d7*d8*d9^2*(-1 + n3)) - 
        d3^(2 - n3)/(6*d10*d4*d5*d6*d7*d8*d9^2*(-1 + n3)) + 
        d3^(1 - n3)/(d10^2*d4*d5*d6*d7*d9*(-1 + n3)) - 
        (d1*d3^(1 - n3))/(6*d10^2*d4*d5*d6*d7*d9*(-1 + n3)) + 
        d3^(2 - n3)/(3*d10^2*d4*d5*d6*d7*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d4*d5*d7*d8^2*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d4*d6*d7*d8^2*d9*(-1 + n3)) - 
        (d1*d3^(1 - n3))/(6*d10*d4^2*d5*d6*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10^2*d4*d5*d6*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(6*d10^2*d4*d5*d6*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10*d4*d6*d7^2*d8*d9*(-1 + n3)) - 
        d3^(1 - n3)/(6*d10^2*d4*d5*d7*d8*d9*(-1 + n3)) - 
        d3^(2 - n3)/(6*d10^2*d4*d5*d7*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10*d4^2*d6*d7*d8*d9*(-1 + n3)) + 
        d3^(2 - n3)/(6*d10*d4^2*d6*d7*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10^2*d4*d6*d7*d8*d9*(-1 + n3)) + 
        d3^(2 - n3)/(6*d10^2*d4*d6*d7*d8*d9*(-1 + n3)) + 
        d3^(1 - n3)/(6*d10^2*d5*d6*d7*d8*d9*(-1 + n3)) + 
        d3^(2 - n3)/(6*d10^2*d5*d6*d7*d8*d9*(-1 + n3)) - 
        (5*d3^(2 - n3))/(6*d10^2*d4*d5*d6*d7*d8*d9*(-1 + n3)) - 
        d3^(3 - n3)/(6*d10^2*d4*d5*d6*d7*d8*d9*(-1 + n3)) + 
        rat(-3, 4*(2 - d + n3))/(d3^n3*d4*d5*d6*d7*d8*d9^2) + 
        (d3^(1 - n3)*rat(-2, 3*(-1 + n3)*(2 - d + n3)))/(d10*d4*d5*d6*d8*d9^3) + 
        (d3^(2 - n3)*rat(-2, 3*(-1 + n3)*(2 - d + n3)))/(d10*d4^3*d5*d6*d8*d9) + 
        (d3^(1 - n3)*rat(-2, 3*(-1 + n3)*(2 - d + n3)))/(d10^2*d5*d6^2*d7*d8*d9) + 
        rat(-1, 3*(2 - d + n3))/(d10^2*d3^n3*d4*d5*d6*d7*d9) + 
        rat(-1, 6*(2 - d + n3))/(d10*d3^n3*d4*d6*d7*d8^2*d9) + 
        rat(-1, 6*(2 - d + n3))/(d10*d3^n3*d5*d6*d7*d8^2*d9) + 
        rat(-1, 6*(2 - d + n3))/(d10*d3^n3*d4*d6*d7^2*d8*d9) + 
        rat(-1, 6*(2 - d + n3))/(d10*d3^n3*d4^2*d6*d7*d8*d9) + 
        rat(-1, 6*(2 - d + n3))/(d10^2*d3^n3*d4*d6*d7*d8*d9) + 
        rat(-1, 6*(2 - d + n3))/(d10^2*d3^n3*d5*d6*d7*d8*d9) + 
        rat(-1, 12*(2 - d + n3))/(d10*d3^n3*d4*d5*d6*d8*d9^2) + 
        (d3^(1 - n3)*rat(-1, 2*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d5^2*d6*d8*d9) + 
        (d3^(1 - n3)*rat(-1, 3*(-1 + n3)*(2 - d + n3)))/(d10*d4*d5*d6^2*d7*d9^2) + 
        (d3^(2 - n3)*rat(-1, 3*(-1 + n3)*(2 - d + n3)))/(d10^2*d4^2*d5*d6*d7*d9) + 
        (d3^(1 - n3)*rat(-1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d5*d6*d7^2*d9) + 
        (d3^(1 - n3)*rat(-1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d5^2*d6*d7*d9) + 
        (d3^(1 - n3)*rat(-1, 6*(-1 + n3)*(2 - d + n3)))/(d10*d4^2*d6*d7*d8^2*d9) + 
        (d3^(1 - n3)*rat(-1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d6*d7*d8^2*d9) + 
        (d3^(1 - n3)*rat(-1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d5*d6*d7*d8^2*d9) + 
        (d3^(1 - n3)*rat(-1, 6*(-1 + n3)*(2 - d + n3)))/(d10*d4^2*d6*d7^2*d8*d9) + 
        (d3^(1 - n3)*rat(-1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d6*d7^2*d8*d9) + 
        rat(1, 3*(2 - d + n3))/(d3^n3*d4*d5*d6*d7*d8^2*d9) + 
        rat(1, 4*(2 - d + n3))/(d10*d3^n3*d4*d6^2*d7*d8*d9) + 
        rat(1, 6*(2 - d + n3))/(d10*d3^n3*d4*d5*d6*d7^2*d8) + 
        rat(1, 6*(2 - d + n3))/(d10*d3^n3*d4^2*d5*d6*d7*d8) + 
        rat(1, 6*(2 - d + n3))/(d10*d3^n3*d4*d5*d6*d8^2*d9) + 
        rat(1, 6*(2 - d + n3))/(d10*d3^n3*d4*d5*d7*d8^2*d9) + 
        rat(1, 6*(2 - d + n3))/(d10^2*d3^n3*d4*d5*d6*d8*d9) + 
        rat(1, 6*(2 - d + n3))/(d10^2*d3^n3*d4*d5*d7*d8*d9) + 
        rat(1, 12*(2 - d + n3))/(d10*d3^n3*d4*d6*d7*d8*d9^2) + 
        (d3^(1 - n3)*rat(1, 2*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d5^2*d6*d7*d8) + 
        (d3^(1 - n3)*rat(1, 3*(-1 + n3)*(2 - d + n3)))/(d10*d4*d6^2*d7*d8*d9^2) + 
        (d3^(1 - n3)*rat(1, 3*(-1 + n3)*(2 - d + n3)))/(d10^2*d4^2*d5*d6*d7*d9) + 
        (d3^(1 - n3)*rat(1, 3*(-1 + n3)*(2 - d + n3)))/(d10*d4*d6^2*d7*d8^2*d9) + 
        (d3^(1 - n3)*rat(1, 3*(-1 + n3)*(2 - d + n3)))/(d10*d4*d6^2*d7^2*d8*d9) + 
        (d3^(3 - n3)*rat(1, 3*(-1 + n3)*(2 - d + n3)))/(d10^2*d4^2*d5*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(1, 6*(-1 + n3)*(2 - d + n3)))/(d10*d4^2*d5*d6*d7^2*d8) + 
        (d3^(1 - n3)*rat(1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d5*d6*d7^2*d8) + 
        (d3^(1 - n3)*rat(1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d5*d6*d8^2*d9) + 
        (d3^(1 - n3)*rat(1, 6*(-1 + n3)*(2 - d + n3)))/(d10*d4^2*d5*d7*d8^2*d9) + 
        (d3^(1 - n3)*rat(1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d5*d7*d8^2*d9) + 
        (d3^(1 - n3)*rat(1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d5*d6*d7^2*d8*d9) + 
        (d3^(1 - n3)*rat(1, 6*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d5^2*d7*d8*d9) + 
        (d3^(1 - n3)*rat(2, 3*(-1 + n3)*(2 - d + n3)))/(d10^2*d4*d5*d6^2*d7*d8) + 
        (d3^(1 - n3)*rat(2, 3*(-1 + n3)*(2 - d + n3)))/(d10*d4*d6*d7*d8*d9^3) + 
        (d3^(1 - n3)*rat(2, 3*(-1 + n3)*(2 - d + n3)))/(d10*d4^3*d5*d6*d8*d9) + 
        (d3^(3 - n3)*rat(2, 3*(-1 + n3)*(2 - d + n3)))/(d10*d4^3*d5*d6*d7*d8*d9) + 
        rat(3, 4*(2 - d + n3))/(d10*d3^n3*d4*d5*d6*d7*d9^2) + 
        (d3^(2 - n3)*rat(9 - 2*d, 3*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(4 - d, 3*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5*d6^2*d7*d9) + 
        (d3^(1 - n3)*rat(-4 + d, 3*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d6^2*d7*d8*d9) + 
        (d3^(1 - n3)*rat(16 - d - 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5*d6^2*d8^2*d9) + 
        (d3^(1 - n3)*rat(16 - d - 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d6*d8^2*d9) + 
        (d3^(1 - n3)*rat(16 - d - 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d6^2*d8*d9) + 
        rat(8 + d - 6*n3, 12*(-4 + d)*(-2 + d - n3))/(d3^n3*d4*d5*d6*d7^2*d8*d9) + 
        rat(8 + d - 6*n3, 12*(-4 + d)*(-2 + d - n3))/(d3^n3*d4^2*d5*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(8 + d - 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d4^2*d5*d6*d7^2*d8*d9) + 
        (d3^(1 - n3)*rat(-8 + 5*d - 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d7*d8^2*d9) + 
        rat(-16 + 7*d - 6*n3, 12*(-4 + d)*(-2 + d - n3))/(d10*d3^n3*d4^2*d5*d6*d8*d9) + 
        (d3^(1 - n3)*rat(14 - 2*d - 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5*d6*d7^2*d9^2) + 
        (d3^(1 - n3)*rat(14 - 2*d - 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d6*d7*d9^2) + 
        (d3^(1 - n3)*rat(14 - 2*d - 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d6*d7^2*d9) + 
        rat(2 + d - 3*n3, 6*(-4 + d)*(-2 + d - n3))/(d3^n3*d4*d5*d6^2*d7*d8*d9) + 
        (d3^(1 - n3)*rat(2 + d - 3*n3, 4*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5*d6^2*d7^2*d8) + 
        (d3^(1 - n3)*rat(2 + d - 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5*d6^2*d7*d8) + 
        (d3^(2 - n3)*rat(3 + d - 3*n3, 6*(-2 + n3)*(-1 + n3)))/(d10*d4*d5*d6*d7^2*d8*d9) + 
        (d3^(1 - n3)*rat(-2 + 2*d - 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5^2*d6*d8*d9) + 
        (d3^(1 - n3)*rat(-10 + 4*d - 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d6*d7^2*d8) + 
        (d3^(1 - n3)*rat(-14 + 5*d - 3*n3, 6*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5*d6*d8*d9) + 
        (d3^(1 - n3)*rat(-22 + 7*d - 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d4^2*d5*d6*d7*d8*d9^2) + 
        rat(16 - 3*d - 2*n3, 4*(-4 + d)*(-2 + d - n3))/(d10*d3^n3*d4*d5*d6^2*d8*d9) + 
        (d3^(1 - n3)*rat(d - 2*n3, 4*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5*d6*d8^2*d9) + 
        (d3^(1 - n3)*rat(1 + d - 2*n3, 6*(-1 + n3)))/(d10*d4*d5*d6*d8*d9) + 
        (d3^(1 - n3)*rat(10 - 2*d - n3, 4*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5^2*d7*d8*d9) + 
        (d3^(2 - n3)*rat(-1 + d - n3, 3*(-2 + n3)*(-1 + n3)))/(d10*d4*d5^2*d6*d7*d8*d9) + 
        (d3^(2 - n3)*rat(-1 + d - n3, 6*(-2 + n3)*(-1 + n3)))/(d10*d4*d6*d7*d8^2*d9) + 
        (d3^(2 - n3)*rat(-1 + d - n3, 6*(-2 + n3)*(-1 + n3)))/(d10*d4*d5*d6*d7*d8^2*d9) + 
        (d3^(2 - n3)*rat(-1 + d - n3, 6*(-2 + n3)*(-1 + n3)))/(d10*d4*d6*d7^2*d8*d9) + 
        (d3^(2 - n3)*rat(1 - d + n3, 6*(-2 + n3)*(-1 + n3)))/(d10*d4*d5*d6*d7^2*d8) + 
        (d3^(2 - n3)*rat(1 - d + n3, 6*(-2 + n3)*(-1 + n3)))/(d10*d4*d5*d7*d8^2*d9) + 
        (d3^(1 - n3)*rat(-10 + 2*d + n3, 4*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5^2*d6*d7*d9) + 
        (d3^(1 - n3)*rat(-1 - d + 2*n3, 6*(-1 + n3)))/(d10*d5*d6*d7*d8*d9) + 
        (d3^(1 - n3)*rat(-d + 2*n3, 6*(-1 + n3)))/(d10*d4*d5*d6*d7*d8) + 
        (d3^(1 - n3)*rat(-d + 2*n3, 4*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d4^2*d5*d6*d7*d8^2*d9) + 
        rat(-16 + 3*d + 2*n3, 4*(-4 + d)*(-2 + d - n3))/(d10*d3^n3*d5*d6^2*d7*d8*d9) + 
        (d3^(1 - n3)*rat(22 - 7*d + 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5*d7*d8*d9^2) + 
        (d3^(1 - n3)*rat(10 - 4*d + 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d6*d8*d9^2) + 
        (d3^(1 - n3)*rat(2 - 2*d + 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5^2*d6*d7*d8) + 
        (d3^(1 - n3)*rat(-2 - d + 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d4^2*d5*d6^2*d7*d8*d9) + 
        (d3^(1 - n3)*rat(-14 + 2*d + 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5*d7^2*d8*d9^2) + 
        (d3^(1 - n3)*rat(-14 + 2*d + 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d7*d8*d9^2) + 
        (d3^(1 - n3)*rat(-14 + 2*d + 3*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d7^2*d8*d9) + 
        rat(16 - 7*d + 6*n3, 12*(-4 + d)*(-2 + d - n3))/(d10*d3^n3*d5*d6*d7^2*d8*d9) + 
        rat(-8 - d + 6*n3, 12*(-4 + d)*(-2 + d - n3))/(d10*d3^n3*d4*d5*d6*d7^2*d9) + 
        rat(-8 - d + 6*n3, 12*(-4 + d)*(-2 + d - n3))/(d10*d3^n3*d4^2*d5*d6*d7*d9) + 
        (d3^(1 - n3)*rat(-8 - d + 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4^2*d5*d6*d7^2*d9) + 
        (d3^(1 - n3)*rat(-8 - d + 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d6^2*d7*d9) + 
        rat(-16 + d + 6*n3, 12*(-4 + d)*(-2 + d - n3))/(d10*d3^n3*d4*d5*d6^2*d7*d9) + 
        (d3^(1 - n3)*rat(-16 + d + 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5*d6^2*d7*d8^2) + 
        (d3^(1 - n3)*rat(-16 + d + 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d6*d7*d8^2) + 
        (d3^(1 - n3)*rat(-16 + d + 6*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5^2*d6^2*d7*d8) + 
        (d3^(1 - n3)*rat(-22 + d + 9*n3, 12*(-4 + d)*(-2 + d - n3)*(-1 + n3)))/(d10*d4*d5*d6^2*d7^2*d9) + 
        (d3^(1 - n3)*rat(24 - 6*d - 2*d^2 - 12*n3 + 11*d*n3 - 8*n3^2, 6*(-1 + n3)*(2 - d + n3)))/(d10*d4*d5*d6*d7*d8*d9) + 
        (d3^(2 - n3)*rat(d - d^2 - 2*n3 + 3*d*n3 - 2*n3^2, 6*(-2 + n3)*(-1 + n3)))/(d10*d4*d5*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n2 == 0 && n3 == 1 && n5 == 1 && n6 == 2 && n7 == 1 && n8 == 1 && n9 == 1
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3/d4^n4?pos_/d5/d6^2/d7/d8/d9/d10 =
        
        -2/(3*d10^2*d3*d4^n4*d5*d6*d7*d9) + 
        2/(3*d10^2*d4^n4*d5*d6*d7*d8*d9) + 
        2/(d10*d3^2*d4^n4*d5*d6*d7*d8*d9) - 
        (d4^(-1 - n4)*n4)/(d10*d3*d5*d6*d8*d9) + 
        (d4^(-1 - n4)*n4)/(d10*d5*d6*d7*d8*d9) + 
        (d4^(-1 - n4)*n4)/(d10*d3*d5*d6*d7*d8*d9) + 
        rat(-4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6^3*d7*d8) + 
        rat(-4, 3*(-3 + d - n4))/(d10^2*d3*d4^n4*d5*d6^2*d7*d8) + 
        rat(-4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6^3*d7*d9) + 
        (d1*rat(-4, 3*(-3 + d - n4)))/(d10*d3*d4^n4*d5*d6^3*d7*d9) + 
        (d4^(1 - n4)*rat(-4, 3*(-3 + d - n4)))/(d10*d3^2*d5*d6^2*d7*d9) + 
        rat(-4, 3*(-3 + d - n4))/(d10*d4^n4*d5*d6^2*d7*d8^2*d9) + 
        (d4^(1 - n4)*rat(-4, 3*(-3 + d - n4)))/(d10*d3*d6^3*d7*d8*d9) + 
        rat(-4, 3*(-3 + d - n4))/(d10*d4^n4*d5*d6^3*d7*d8*d9) + 
        rat(-4, 3*(-3 + d - n4))/(d10^2*d3*d4^n4*d6^2*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-4, 3*(-3 + d - n4)))/(d10^2*d5*d6^2*d7*d8*d9) + 
        rat(-2, -3 + d - n4)/(d10^2*d3*d4^n4*d5*d6^2*d7*d9) + 
        (d4^(1 - n4)*rat(-2, -3 + d - n4))/(d10*d3^2*d5*d6^2*d7*d8*d9) + 
        (d1*rat(-2, 3*(-3 + d - n4)))/(d10*d3*d4^n4*d5*d6^2*d7^2*d8) + 
        rat(-2, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6*d7^2*d8) + 
        rat(-2, 3*(-3 + d - n4))/(d10*d3^2*d4^n4*d5*d6^2*d7*d8) + 
        rat(-2, 3*(-3 + d - n4))/(d10*d3*d4^n4*d6^2*d7*d8*d9^2) + 
        rat(-2, 3*(-3 + d - n4))/(d10*d4^n4*d5*d6^2*d7*d8*d9^2) + 
        (d2*rat(-2, 3*(-3 + d - n4)))/(d10*d4^n4*d5*d6^2*d7*d8*d9^2) + 
        (d4^(1 - n4)*rat(-2, 3*(-3 + d - n4)))/(d10*d3*d5*d6^2*d7^2*d9) + 
        (d4^(1 - n4)*rat(-2, 3*(-3 + d - n4)))/(d10*d3*d5^2*d6^2*d7*d9) + 
        rat(-2, 3*(-3 + d - n4))/(d10*d3*d4^n4*d6^2*d7*d8^2*d9) + 
        rat(-2, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5^2*d6^2*d8*d9) + 
        (d4^(1 - n4)*rat(-2, 3*(-3 + d - n4)))/(d10*d3*d6^2*d7^2*d8*d9) + 
        (d2*rat(-2, 3*(-3 + d - n4)))/(d10*d3*d4^n4*d6^2*d7^2*d8*d9) + 
        (d4^(1 - n4)*rat(-2, 3*(-3 + d - n4)))/(d10*d5*d6^2*d7^2*d8*d9) + 
        rat(-2, 3*(-3 + d - n4))/(d10*d4^n4*d5*d6^2*d7^2*d8*d9) + 
        (d2*rat(-2, 3*(-3 + d - n4)))/(d10*d4^n4*d5*d6^2*d7^2*d8*d9) + 
        rat(-2, 3*(-3 + d - n4))/(d3^2*d4^n4*d5*d6^2*d7*d8*d9) + 
        (d1*rat(-2, 3*(-3 + d - n4)))/(d3^2*d4^n4*d5*d6^2*d7*d8*d9) + 
        (d2*rat(2, 3*(-3 + d - n4)))/(d10*d3*d4^n4*d5*d6^2*d7^2*d8) + 
        rat(2, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5^2*d6^2*d7*d8) + 
        rat(2, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6^2*d8*d9^2) + 
        rat(2, 3*(-3 + d - n4))/(d3*d4^n4*d5*d6^2*d7*d8*d9^2) + 
        (d1*rat(2, 3*(-3 + d - n4)))/(d3*d4^n4*d5*d6^2*d7*d8*d9^2) + 
        rat(2, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5^2*d6^2*d7*d9) + 
        (d1*rat(2, 3*(-3 + d - n4)))/(d10*d3*d4^n4*d5^2*d6^2*d7*d9) + 
        rat(2, 3*(-3 + d - n4))/(d10*d3^2*d4^n4*d5*d6^2*d7*d9) + 
        (d1*rat(2, 3*(-3 + d - n4)))/(d10*d3^2*d4^n4*d5*d6^2*d7*d9) + 
        rat(2, 3*(-3 + d - n4))/(d10*d3^2*d4^n4*d5*d6^2*d8*d9) + 
        (d1*rat(2, 3*(-3 + d - n4)))/(d10*d3*d4^n4*d6^2*d7^2*d8*d9) + 
        (d4^(1 - n4)*rat(2, 3*(-3 + d - n4)))/(d3*d5*d6^2*d7^2*d8*d9) + 
        (d1*rat(2, 3*(-3 + d - n4)))/(d10*d4^n4*d5*d6^2*d7^2*d8*d9) + 
        rat(2, 3*(-3 + d - n4))/(d10*d3*d4^n4*d6*d7^2*d8*d9) + 
        rat(2, 3*(-3 + d - n4))/(d10*d4^n4*d5*d6*d7^2*d8*d9) + 
        rat(2, 3*(-3 + d - n4))/(d10^2*d4^n4*d5*d6^2*d7*d8*d9) + 
        (d4^(1 - n4)*rat(4, 3*(-3 + d - n4)))/(d10*d3*d5*d6^3*d7*d9) + 
        (d4^(1 - n4)*rat(4, 3*(-3 + d - n4)))/(d10^2*d3*d5*d6^2*d7*d9) + 
        rat(4, 3*(-3 + d - n4))/(d3*d4^n4*d5*d6^2*d7*d8^2*d9) + 
        rat(4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6^3*d8*d9) + 
        rat(4, 3*(-3 + d - n4))/(d10^2*d3*d4^n4*d5*d6^2*d8*d9) + 
        (d1*rat(4, 3*(-3 + d - n4)))/(d10*d3*d4^n4*d6^3*d7*d8*d9) + 
        rat(4, 3*(-3 + d - n4))/(d3*d4^n4*d5*d6^3*d7*d8*d9) + 
        (d4^(1 - n4)*rat(4, 3*(-3 + d - n4)))/(d3^2*d5*d6^2*d7*d8*d9) + 
        rat(-15 + 5*d - 9*n4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6^2*d7*d9) + 
        rat(-8 + 4*d - 8*n4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6*d8*d9^2) + 
        rat(-10 + 4*d - 6*n4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6^2*d7*d8) + 
        rat(-6 + 4*d - 6*n4, 3*(-3 + d - n4))/(d10*d4^n4*d5*d6^2*d7*d8*d9) + 
        rat(-4 + 2*d - 4*n4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5^2*d6*d7*d8) + 
        rat(-4 + 2*d - 2*n4, 3*(-3 + d - n4))/(d10*d4^n4*d5*d6*d7*d8^2*d9) + 
        (d4^(1 - n4)*rat(-4 + 2*d - 2*n4, 3*(-3 + d - n4)))/(d10*d3*d5*d6^2*d7*d8*d9) + 
        (d4^(-1 - n4)*rat(-4*n4, 3*(-3 + d - n4)))/(d10*d3*d6^2*d7*d8*d9) + 
        (d2*d4^(-1 - n4)*rat(-2*n4, 3*(-3 + d - n4)))/(d10*d3*d6^2*d7*d8*d9) + 
        (d1*d4^(-1 - n4)*rat(-2*n4, 3*(-3 + d - n4)))/(d10*d5*d6^2*d7*d8*d9) + 
        (d1*d4^(-1 - n4)*rat(2*n4, 3*(-3 + d - n4)))/(d10*d3*d5*d6^2*d8*d9) + 
        (d4^(-1 - n4)*rat(2*n4, 3*(-3 + d - n4)))/(d10*d3*d6*d7*d8*d9) + 
        rat(4 - 2*d + 2*n4, 3*(-3 + d - n4))/(d3*d4^n4*d5*d6*d7*d8^2*d9) + 
        rat(4 - 2*d + 4*n4, -3 + d - n4)/(d10*d3*d4^n4*d5^2*d6*d7*d8*d9) + 
        rat(4 - 2*d + 4*n4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5^2*d6*d8*d9) + 
        rat(3 - d + 5*n4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d6^2*d7*d8*d9) + 
        rat(8 - 4*d + 6*n4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6^2*d8*d9) + 
        rat(8 - 4*d + 6*n4, 3*(-3 + d - n4))/(d3*d4^n4*d5*d6^2*d7*d8*d9) + 
        rat(8 - 4*d + 8*n4, 3*(-3 + d - n4))/(d10*d3*d4^n4*d6*d7*d8*d9^2) + 
        rat(1 - d + 2*n4 + d*n4 - 3*n4^2, 3*(-3 + d - n4))/(d10*d3*d4^n4*d5*d6*d7*d8*d9);
        
* n1 == 0 && n10 == 1 && n2 == 0 && n3 == 1 && n5 == 2 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n4 != 1
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3/d4^n4?{>1}/d5^2/d6/d7/d8/d9/d10 =
        
        (d4^(1 - n4)*rat(-3, (-1 + n4)*(-14 + d + 5*n4)))/(d10^2*d5^2*d6*d7*d8^2*d9) + 
        (d4^(1 - n4)*rat(-3, 2*(-1 + n4)*(-14 + d + 5*n4)))/(d10*d3^2*d5*d6^2*d7*d8*d9) + 
        rat(-1, -2 + d - n4)/(d10*d3*d4^n4*d5^2*d7*d8*d9^2) + 
        (d4^(1 - n4)*rat(-1, (-2 + d - n4)*(-1 + n4)))/(d10*d3*d5^3*d6*d7^2*d8) + 
        (d4^(1 - n4)*rat(-1, (-2 + d - n4)*(-1 + n4)))/(d10*d3*d5^3*d7*d8*d9^2) + 
        (d4^(1 - n4)*rat(-1, (-2 + d - n4)*(-1 + n4)))/(d10*d3*d5^3*d7^2*d8*d9) + 
        (d4^(1 - n4)*rat(-1, 2*(-2 + d - n4)*(-1 + n4)))/(d10*d3^2*d5^2*d6*d7*d9^2) + 
        (d4^(1 - n4)*rat(-1, 2*(-2 + d - n4)*(-1 + n4)))/(d10*d3*d5^2*d7^2*d8*d9^2) + 
        rat(1, -2 + d - n4)/(d3*d4^n4*d5^2*d6*d7*d8*d9^2) + 
        (d4^(1 - n4)*rat(1, (-2 + d - n4)*(-1 + n4)))/(d10*d3*d5^3*d6*d7*d9^2) + 
        (d4^(1 - n4)*rat(1, (-2 + d - n4)*(-1 + n4)))/(d10*d3*d5^3*d6*d8*d9^2) + 
        (d4^(1 - n4)*rat(1, (-2 + d - n4)*(-1 + n4)))/(d10*d3*d5^3*d6*d7^2*d9) + 
        (d4^(1 - n4)*rat(1, 2*(-2 + d - n4)*(-1 + n4)))/(d10*d3*d5^2*d6*d7^2*d9^2) + 
        (d4^(1 - n4)*rat(1, 2*(-2 + d - n4)*(-1 + n4)))/(d3^2*d5^2*d6*d7*d8*d9^2) + 
        (d4^(1 - n4)*rat(3, (-1 + n4)*(-14 + d + 5*n4)))/(d10^2*d3^2*d5^2*d6*d7*d9) + 
        (d4^(1 - n4)*rat(30 + 3*d - 21*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3^2*d5*d6*d7*d8*d9^2) + 
        rat(32 - d - 14*n4, 2*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2))/(d10*d3*d4^n4*d5^2*d6*d7^2*d9) + 
        rat(32 - d - 14*n4, 2*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2))/(d10*d3^2*d4^n4*d5^2*d6*d7*d9) + 
        (d4^(1 - n4)*rat(32 - d - 14*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3^2*d5^2*d6*d7^2*d9) + 
        (d4^(1 - n4)*rat(22 + d - 13*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10^2*d5^2*d6*d7^2*d8*d9) + 
        rat(12 + 3*d - 12*n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10*d4^n4*d5^2*d6*d7^2*d8*d9) + 
        rat(18 - 9*n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10*d4^n4*d5^2*d6*d7*d8*d9^2) + 
        rat(18 - 9*n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10^2*d4^n4*d5^2*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(20 - d - 8*n4, 2*(-1 + n4)*(-14 + d + 5*n4)))/(d10*d3*d5^2*d6*d8*d9^2) + 
        (d4^(1 - n4)*rat(8 + 2*d - 8*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10^2*d5^2*d6*d7*d8*d9^2) + 
        rat(24 - 3*d - 6*n4, -28 + 16*d - d^2 - 4*n4 - 4*d*n4 + 5*n4^2)/(d10*d3*d4^n4*d5^3*d6*d8*d9) + 
        rat(4 + d - 4*n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10*d3*d4^n4*d5^3*d6*d7*d9) + 
        (d4^(1 - n4)*rat(4 + d - 4*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^3*d6*d7*d8^2) + 
        (d4^(1 - n4)*rat(4 + d - 4*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^3*d6^2*d7*d8) + 
        (d4^(1 - n4)*rat(4 + d - 4*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^2*d6^2*d7*d8^2) + 
        (d4^(1 - n4)*rat(4 + d - 4*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3^2*d5^2*d6^2*d7*d9) + 
        (d4^(2 - n4)*rat(4 + d - 4*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3^2*d5^2*d6^2*d7*d8*d9) + 
        rat(6 - 3*n4, -14 + d + 5*n4)/(d10*d3*d4^n4*d5^2*d6*d8*d9) + 
        (d4^(1 - n4)*rat(8 - d - 2*n4, (-1 + n4)*(-14 + d + 5*n4)))/(d10^2*d5^2*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(8 - d - 2*n4, 2*(-1 + n4)*(-14 + d + 5*n4)))/(d10*d3*d5*d6^2*d7*d8*d9) + 
        rat(-16 + 5*d - 2*n4, 2*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2))/(d10*d3*d4^n4*d5^2*d6*d8^2*d9) + 
        (d4^(1 - n4)*rat(-16 + 5*d - 2*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^3*d7*d8^2*d9) + 
        rat(10 - 2*d - n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10*d4^n4*d5^2*d6^2*d7*d8*d9) + 
        (d4^(1 - n4)*rat(10 - 2*d - n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^2*d6^2*d7*d9^2) + 
        (d4^(1 - n4)*rat(10 - 2*d - n4, 2*(-1 + n4)*(-28 + 16*d - d^2 - 4*n4 - 4*d*n4 + 5*n4^2)))/(d10*d3*d5^2*d6^2*d7^2*d8) + 
        (d4^(1 - n4)*rat(10 - 2*d - n4, 2*(-1 + n4)*(-28 + 16*d - d^2 - 4*n4 - 4*d*n4 + 5*n4^2)))/(d3^2*d5^2*d6^2*d7*d8*d9) + 
        rat(-10 + 2*d + n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10*d3*d4^n4*d5^2*d6^2*d8*d9) + 
        (d4^(1 - n4)*rat(-10 + 2*d + n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^2*d6^2*d8*d9^2) + 
        (d4^(1 - n4)*rat(-10 + 2*d + n4, 2*(-1 + n4)*(-28 + 16*d - d^2 - 4*n4 - 4*d*n4 + 5*n4^2)))/(d10*d3*d5^2*d6^2*d7^2*d9) + 
        (d4^(1 - n4)*rat(16 - 5*d + 2*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^3*d6^2*d7*d9) + 
        (d4^(2 - n4)*rat(16 - 5*d + 2*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3^2*d5^2*d6*d7^2*d8*d9) + 
        (d4^(1 - n4)*rat(-8 + d + 2*n4, (-1 + n4)*(-14 + d + 5*n4)))/(d10^2*d3*d5^2*d6*d7*d9) + 
        (d4^(1 - n4)*rat(-8 + d + 2*n4, 2*(-1 + n4)*(-14 + d + 5*n4)))/(d10*d3*d5^2*d6^2*d7*d9) + 
        rat(-6 + 3*n4, -14 + d + 5*n4)/(d10*d4^n4*d5^2*d6*d7*d8*d9) + 
        rat(-4 - d + 4*n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10*d4^n4*d5^2*d6*d7*d8^2*d9) + 
        rat(-4 - d + 4*n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10*d3*d4^n4*d5^3*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-4 - d + 4*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^3*d6*d8^2*d9) + 
        (d4^(1 - n4)*rat(-4 - d + 4*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^3*d6^2*d8*d9) + 
        (d4^(1 - n4)*rat(-4 - d + 4*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3*d5^2*d6^2*d8^2*d9) + 
        (d4^(1 - n4)*rat(-4 - d + 4*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3^2*d5^2*d6^2*d8*d9) + 
        rat(-24 + 3*d + 6*n4, -28 + 16*d - d^2 - 4*n4 - 4*d*n4 + 5*n4^2)/(d10*d3*d4^n4*d5^3*d6*d7*d8) + 
        rat(-24 + 3*d + 6*n4, 2*(-28 + 16*d - d^2 - 4*n4 - 4*d*n4 + 5*n4^2))/(d3*d4^n4*d5^2*d6*d7*d8^2*d9) + 
        (d4^(1 - n4)*rat(-8 - 2*d + 8*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10^2*d3*d5^2*d6*d7*d9^2) + 
        (d4^(1 - n4)*rat(-20 + d + 8*n4, 2*(-1 + n4)*(-14 + d + 5*n4)))/(d10*d3*d5*d6*d7*d8*d9^2) + 
        rat(-18 + 9*n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10*d3*d4^n4*d5^2*d6*d8*d9^2) + 
        rat(-18 + 9*n4, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)/(d10^2*d3*d4^n4*d5^2*d6*d7*d9) + 
        (d4^(1 - n4)*rat(-22 - d + 13*n4, (-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10^2*d3*d5^2*d6*d7^2*d9) + 
        rat(-32 + d + 14*n4, 2*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2))/(d3*d4^n4*d5^2*d6*d7^2*d8*d9) + 
        rat(-32 + d + 14*n4, 2*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2))/(d3^2*d4^n4*d5^2*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-32 + d + 14*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d3^2*d5^2*d6*d7^2*d8*d9) + 
        (d4^(1 - n4)*rat(-30 - 3*d + 21*n4, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3^2*d5^2*d6*d8*d9^2) + 
        rat(-40 - d + 22*n4, 2*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2))/(d10*d3^2*d4^n4*d5^2*d6*d8*d9) + 
        (d4^(1 - n4)*rat(24 - 20*d + d^2 + 16*n4 + 6*d*n4 - 10*n4^2, 2*(-1 + n4)*(-14 + d + 5*n4)))/(d10*d3*d5^2*d6*d7*d8*d9) + 
        (d4^(-1 - n4)*rat(18*n4 - 9*n4^2, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2))/(d10*d5^2*d6*d7*d8*d9) + 
        (d4^(-1 - n4)*rat(-18*n4 + 9*n4^2, 28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2))/(d10*d3*d5^2*d6*d8*d9) + 
        (d4^(1 - n4)*rat(-136 + 66*d - 2*d^2 + 4*n4 - 25*d*n4 + 24*n4^2, 2*(-1 + n4)*(28 - 16*d + d^2 + 4*n4 + 4*d*n4 - 5*n4^2)))/(d10*d3^2*d5^2*d6*d7*d8*d9);

* n1 == 0 && n10 == 1 && n2 == 0 && n3 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n4 != 1
        if((count(d1,1)==0) && (count(d2,1)==0)) id,only,ifmatch->sortme  1/d3/d4^n4?{>1}/d5/d6/d7/d8/d9/d10 =
        
        rat(-3, 2*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5^2*d6*d7*d8) + 
        rat(-3, 2*(-8 + d + 2*n4))/(d3*d4^n4*d5*d6*d7*d8^2*d9) + 
        (d4^(1 - n4)*rat(-3, 2*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d5*d6^2*d7^2*d9) + 
        rat(-1, 2 - d + n4)/(d3*d4^n4*d5*d6*d7*d8*d9^2) + 
        rat(-1, 2*(2 - d + n4))/(d10*d3*d4^n4*d5*d6^2*d8*d9) + 
        rat(-1, 2*(2 - d + n4))/(d10*d4^n4*d5^2*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-1, (-1 + n4)*(2 - d + n4)))/(d10*d3*d6^3*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5*d6^2*d7*d8^2) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5^2*d6*d7*d8^2) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5^2*d6^2*d7*d8) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5*d6*d7^2*d9^2) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5^2*d6*d7*d9^2) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5*d6^2*d8*d9^2) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d6^2*d7*d8*d9^2) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d3^2*d5*d6*d7*d8*d9^2) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5^2*d6*d7^2*d9) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d6^2*d7^2*d8*d9) + 
        (d4^(1 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10^2*d5*d6^2*d7*d8*d9) + 
        (d4^(2 - n4)*rat(-1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3^2*d5*d6^2*d7*d8*d9) + 
        rat(1, 2 - d + n4)/(d10*d3*d4^n4*d5*d7*d8*d9^2) + 
        rat(1, 2 - d + n4)/(d10*d4^n4*d5*d6*d7*d8*d9^2) + 
        rat(1, 2*(2 - d + n4))/(d10*d4^n4*d5*d6^2*d7*d8*d9) + 
        (d4^(1 - n4)*rat(1, (-1 + n4)*(2 - d + n4)))/(d10*d3*d5*d6^2*d7*d9^2) + 
        (d4^(1 - n4)*rat(1, (-1 + n4)*(2 - d + n4)))/(d10*d3*d5*d6^3*d7*d9) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5^2*d6*d7^2*d8) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3^2*d5*d6*d7*d9^2) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5*d7^2*d8*d9^2) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5^2*d7*d8*d9^2) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10^2*d3*d5*d6^2*d7*d9) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5*d6^2*d8^2*d9) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5^2*d6*d8^2*d9) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5^2*d6^2*d8*d9) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3^2*d5*d6^2*d8*d9) + 
        (d4^(1 - n4)*rat(1, 2*(-1 + n4)*(2 - d + n4)))/(d10*d3*d5^2*d7^2*d8*d9) + 
        rat(3, 2*(-8 + d + 2*n4))/(d10*d4^n4*d5*d6*d7*d8^2*d9) + 
        (d4^(1 - n4)*rat(3, (-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d6*d7^2*d8*d9^2) + 
        (d4^(1 - n4)*rat(4 - d, 2*(-2 + d - n4)*(-1 + n4)))/(d10^2*d5*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-4 + d, 2*(-2 + d - n4)*(-1 + n4)))/(d10^2*d3*d5*d6*d7*d9) + 
        (d4^(1 - n4)*rat(48 + 3*d - 30*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3^2*d6*d7*d8*d9^2) + 
        rat(46 + d - 25*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d4^n4*d5*d6*d7^2*d8*d9) + 
        rat(28 + d - 16*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10^2*d4^n4*d5*d6*d7*d8*d9) + 
        rat(12 + 3*d - 12*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d6*d7*d8*d9^2) + 
        rat(12 + 3*d - 12*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d6^2*d7*d8*d9) + 
        rat(26 - d - 11*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5*d6*d7^2*d9) + 
        rat(26 - d - 11*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3^2*d4^n4*d5*d6*d7*d9) + 
        (d4^(1 - n4)*rat(26 - d - 11*n4, (-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10^2*d5*d6*d7^2*d8*d9) + 
        (d4^(1 - n4)*rat(26 - d - 11*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3^2*d5*d6*d7^2*d9) + 
        rat(18 - 9*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5*d6*d7^2*d8) + 
        rat(18 - 9*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3^2*d4^n4*d5*d6*d7*d8) + 
        rat(18 - 9*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5*d7*d8^2*d9) + 
        (d4^(1 - n4)*rat(18 - 9*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3^2*d5*d6*d7^2*d8) + 
        (d4^(1 - n4)*rat(18 - 9*n4, 2*(-1 + n4)*(2 - d + n4)*(-8 + d + 2*n4)))/(d10*d3^2*d6*d7^2*d8*d9) + 
        (d4^(1 - n4)*rat(-4 + 5*d - 8*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d5^2*d7*d8^2*d9) + 
        rat(10 + d - 7*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5^2*d6*d7*d9) + 
        rat(2 + 2*d - 5*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5^2*d6*d8*d9) + 
        (d4^(1 - n4)*rat(2 + 2*d - 5*n4, (-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10^2*d3*d5^2*d6*d7*d9) + 
        (d4^(1 - n4)*rat(2 + 2*d - 5*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d5*d6^2*d7^2*d8) + 
        (d4^(1 - n4)*rat(2 + 2*d - 5*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3^2*d6^2*d7*d8*d9) + 
        (d4^(1 - n4)*rat(2 + 2*d - 5*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d3^2*d5*d6^2*d7*d8*d9) + 
        rat(4 + d - 4*n4, 2*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5*d6*d8*d9) + 
        (d4^(1 - n4)*rat(4 + d - 4*n4, (-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d5*d6*d8*d9^2) + 
        rat(-4 - d + 4*n4, 2*(-8 + d + 2*n4))/(d10*d4^n4*d5*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-4 - d + 4*n4, (-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d6*d7*d8*d9^2) + 
        (d4^(1 - n4)*rat(-2 - 2*d + 5*n4, (-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3^2*d5*d6^2*d7*d9) + 
        (d4^(1 - n4)*rat(-2 - 2*d + 5*n4, (-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10^2*d5^2*d6*d7*d8*d9) + 
        rat(-10 - d + 7*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5^2*d7*d8*d9) + 
        (d4^(1 - n4)*rat(4 - 5*d + 8*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d5^2*d6*d8*d9^2) + 
        (d4^(1 - n4)*rat(4 - 5*d + 8*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d6^2*d7*d8^2*d9) + 
        (d4^(2 - n4)*rat(4 - 5*d + 8*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3^2*d5*d6*d7^2*d8*d9) + 
        rat(-18 + 9*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d6*d7*d8^2*d9) + 
        rat(-18 + 9*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d6*d7^2*d8*d9) + 
        rat(-18 + 9*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3^2*d4^n4*d6*d7*d8*d9) + 
        rat(-26 + d + 11*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d3*d4^n4*d5*d6*d7^2*d8*d9) + 
        rat(-26 + d + 11*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d3^2*d4^n4*d5*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-26 + d + 11*n4, (-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10^2*d3*d5*d6*d7^2*d9) + 
        (d4^(1 - n4)*rat(-26 + d + 11*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d3^2*d5*d6*d7^2*d8*d9) + 
        rat(-12 - 3*d + 12*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5*d6^2*d7*d9) + 
        rat(-28 - d + 16*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3*d4^n4*d5*d6*d8*d9^2) + 
        rat(-28 - d + 16*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10^2*d3*d4^n4*d5*d6*d7*d9) + 
        rat(-50 + 4*d + 17*n4, 2*(-2 + d - n4)*(-8 + d + 2*n4))/(d10*d3^2*d4^n4*d5*d6*d8*d9) + 
        (d4^(1 - n4)*rat(-48 - 3*d + 30*n4, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3^2*d5*d6*d8*d9^2) + 
        (d4^(-1 - n4)*rat(28*n4 + d*n4 - 16*n4^2, 2*(-2 + d - n4)*(-8 + d + 2*n4)))/(d10*d5*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(-24 + 10*d - 2*d^2 + 4*n4 + 3*d*n4 - 4*n4^2, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d6^2*d7*d8*d9) + 
        (d4^(1 - n4)*rat(24 - 14*d + d^2 + 4*n4 + 3*d*n4 - 4*n4^2, 2*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d5*d6*d7*d8*d9) + 
        (d4^(1 - n4)*rat(24 - 10*d + 2*d^2 - 4*n4 - 3*d*n4 + 4*n4^2, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3*d5*d6^2*d7*d9) + 
        (d4^(-1 - n4)*rat(-28*n4 - d*n4 + 16*n4^2, 2*(-2 + d - n4)*(-8 + d + 2*n4)))/(d10*d3*d5*d6*d8*d9) + 
        (d4^(1 - n4)*rat(-88 + 54*d - 2*d^2 - 20*n4 - 19*d*n4 + 24*n4^2, 2*(-2 + d - n4)*(-1 + n4)*(-8 + d + 2*n4)))/(d10*d3^2*d5*d6*d7*d8*d9);
        
* n10 == 1 && n2 == 0 && n3 == 1 && n4 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n1 != 0
        if(count(d2,1)==0) id,only,ifmatch->sortme  1/d1^n1?neg_/d3/d4/d5/d6/d7/d8/d9/d10^n10?pos_ =
        
        (d1^(-1 - n1)*rat(-2, -4 + d - n1))/(d10*d3*d4*d5^2*d6*d8*d9) + 
        (d1^(-1 - n1)*rat(-2, -4 + d - n1))/(d10*d3^2*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(1, -4 + d - n1))/(d10*d3*d4*d5*d6^2*d7*d8) + 
        rat(1, -4 + d - n1)/(d1^n1*d10*d3*d4*d6*d7*d8*d9^2) + 
        (d1^(-1 - n1)*rat(1, -4 + d - n1))/(d10*d4*d5*d6*d7*d8*d9^2) + 
        (d1^(-1 - n1)*rat(1, -4 + d - n1))/(d10*d3*d4*d5*d6*d7*d8*d9^2) + 
        (d1^(-1 - n1)*rat(1, -4 + d - n1))/(d10*d3*d4*d5*d6^2*d7*d9) + 
        (d1^(-1 - n1)*rat(1, -4 + d - n1))/(d3*d4*d5*d6*d7*d8^2*d9) + 
        (d1^(-1 - n1)*rat(1, -4 + d - n1))/(d10*d3*d5*d6*d7^2*d8*d9) + 
        rat(1, -4 + d - n1)/(d1^n1*d10*d3*d4*d6^2*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(1, -4 + d - n1))/(d10*d3*d5*d6^2*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(1, -4 + d - n1))/(d10*d4*d5*d6^2*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(1, -4 + d - n1))/(d10*d3*d4*d5*d6^2*d7*d8*d9) + 
        rat(1, 4 - d + n1)/(d1^n1*d10*d3*d4*d5*d6*d8*d9^2) + 
        rat(1, 4 - d + n1)/(d1^n1*d10*d3*d4*d5*d6^2*d7*d9) + 
        (d1^(-1 - n1)*rat(1, 4 - d + n1))/(d10*d4*d5*d6*d7*d8^2*d9) + 
        (d1^(-1 - n1)*rat(1, 4 - d + n1))/(d10*d3*d4*d5*d6*d7*d8^2*d9) + 
        (d1^(-1 - n1)*rat(1, 4 - d + n1))/(d10*d3*d4*d5*d6^2*d8*d9) + 
        (d1^(-1 - n1)*rat(1, 4 - d + n1))/(d10*d4*d5*d6*d7^2*d8*d9) + 
        (d1^(-1 - n1)*rat(1, 4 - d + n1))/(d10*d3*d4*d5*d6*d7^2*d8*d9) + 
        (d1^(-1 - n1)*rat(1, 4 - d + n1))/(d3*d4*d5*d6^2*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(2, -4 + d - n1))/(d10*d3*d4*d5^2*d6*d7*d8) + 
        (d1^(-1 - n1)*rat(2, -4 + d - n1))/(d10*d3*d5^2*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-5 + d - n1, -4 + d - n1))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(-n1, 4 - d + n1))/(d10*d3*d4*d6*d7*d8*d9) + 
        (d1^(-1 - n1)*rat(n1, 4 - d + n1))/(d10*d3*d5*d6*d7*d8*d9) + 
        (d1^(-2 - n1)*rat(2 + 2*n1, -4 + d - n1))/(d10*d3*d4*d5*d6*d7*d8) + 
        (d1^(-2 - n1)*rat(2 + 2*n1, 4 - d + n1))/(d10*d4*d5*d6*d7*d8*d9) + 
        (d1^(-2 - n1)*rat(2 + 2*n1, 4 - d + n1))/(d10*d3*d4*d5*d6*d7*d8*d9);
        
* n1 == 0 && n10 == 1 && n3 == 1 && n4 == 1 && n5 == 1 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1 && n2 != 0
        if(count(d1,1)==0) id,only,ifmatch->sortme  1/d2^n2?neg_/d3/d4/d5/d6/d7/d8/d9/d10 =
        
        (-2*d2^(-1 - n2))/(d10*d3*d4*d5*d6*d7*d8) - 
        (2*d2^(-1 - n2))/(d10*d3*d4*d5*d6*d7*d9) + 
        d2^(-1 - n2)/(d10*d3*d4*d5*d7*d8*d9) - 
        d2^(-1 - n2)/(d10*d3*d5*d6*d7*d8*d9) + 
        (2*d2^(-1 - n2))/(d3*d4*d5*d6*d7*d8*d9) + 
        (2*d1*d2^(-1 - n2))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        d2^(-1 - n2)/(d10^2*d3*d4*d5*d6*d7*d8*n2) - 
        1/(d10*d2^n2*d3*d4*d5*d6*d8*d9^2*n2) + 
        1/(d10*d2^n2*d3*d4*d6*d7*d8*d9^2*n2) - 
        d2^(-1 - n2)/(d3*d4*d5*d6*d7*d8*d9^2*n2) - 
        1/(d10^2*d2^n2*d3*d4*d5*d6*d7*d9*n2) - 
        d2^(-1 - n2)/(d10*d3*d4*d5*d6*d7*d8^2*d9*n2) - 
        2/(d10*d2^n2*d3*d4^2*d5*d6*d8*d9*n2) + 
        d2^(-1 - n2)/(d10*d3*d5*d6*d7^2*d8*d9*n2) - 
        d2^(-1 - n2)/(d10*d3^2*d5*d6*d7*d8*d9*n2) + 
        2/(d10*d2^n2*d4^2*d5*d6*d7*d8*d9*n2) + 
        1/(d10^2*d2^n2*d4*d5*d6*d7*d8*d9*n2) + 
        (d2^(-1 - n2)*rat(-4, 4 - d + n2))/(d10*d3*d4^2*d5*d6*d7*d9) + 
        (d2^(-1 - n2)*rat(-4, 4 - d + n2))/(d10*d3*d4*d5*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(-4, 4 - d + n2))/(d10*d3*d4^2*d5*d6*d7*d8*d9) + 
        rat(-4, n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4*d5*d6*d8*d9^2) + 
        rat(-4, n2*(4 - d + n2))/(d10^2*d2^n2*d3*d4*d5*d6*d7^2*d9) + 
        rat(-4, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^3*d5*d6*d8*d9) + 
        rat(-4, n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4^2*d5*d6*d8*d9) + 
        (d2^(-1 - n2)*rat(-2, 4 - d + n2))/(d10*d3*d4^2*d5*d6*d7*d8) + 
        (d2^(-1 - n2)*rat(-2, 4 - d + n2))/(d10*d3*d4*d5*d6^2*d7*d9) + 
        (d2^(-1 - n2)*rat(-2, 4 - d + n2))/(d10*d3^2*d4*d5*d6*d7*d9) + 
        (d2^(-1 - n2)*rat(-2, 4 - d + n2))/(d3*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(-2, 4 - d + n2))/(d10*d3*d4*d5*d6^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-2, 4 - d + n2))/(d10^2*d4*d5*d6*d7*d8*d9) + 
        rat(-2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d5*d6*d8*d9^2) + 
        rat(-2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d5*d6^2*d7*d9) + 
        rat(-2, n2*(4 - d + n2))/(d10^2*d2^n2*d3*d4^2*d5*d6*d7*d9) + 
        rat(-2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d6*d7*d8^2*d9) + 
        rat(-2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d6*d7^2*d8*d9) + 
        rat(-2, n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4*d6*d7^2*d8*d9) + 
        rat(-2, n2*(4 - d + n2))/(d2^n2*d3*d4^2*d5*d6*d7^2*d8*d9) + 
        rat(-2, n2*(4 - d + n2))/(d2^n2*d3^2*d4*d5*d6*d7^2*d8*d9) + 
        rat(-2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d5^2*d7*d8*d9) + 
        rat(-2, n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4^2*d6*d7*d8*d9) + 
        rat(-2, n2*(4 - d + n2))/(d2^n2*d3^2*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(2, 4 - d + n2))/(d10*d3^2*d4*d5*d6*d7*d8) + 
        (d2^(-1 - n2)*rat(2, 4 - d + n2))/(d10^2*d3*d4*d5*d6*d7*d9) + 
        (d2^(-1 - n2)*rat(2, 4 - d + n2))/(d10*d3*d4*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(2, 4 - d + n2))/(d10*d3*d4*d6^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(2, 4 - d + n2))/(d3*d4^2*d5*d6*d7*d8*d9) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d5*d6*d7^2*d8) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4*d5*d6*d7^2*d8) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4^2*d5*d6*d7*d8) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d6*d7*d8*d9^2) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d5*d6*d7^2*d9) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4*d5*d6*d7^2*d9) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d5^2*d6*d7*d9) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4^2*d5*d6*d7*d9) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d5*d7*d8^2*d9) + 
        rat(2, n2*(4 - d + n2))/(d10*d2^n2*d3*d4^2*d6^2*d7*d8*d9) + 
        rat(2, n2*(4 - d + n2))/(d10^2*d2^n2*d4^2*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(4, 4 - d + n2))/(d10*d3^2*d4*d5*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(4, 4 - d + n2))/(d10*d3*d4^2*d6*d7*d8*d9) + 
        rat(4, n2*(4 - d + n2))/(d10*d2^n2*d3^2*d4*d6*d7*d8*d9^2) + 
        rat(4, n2*(4 - d + n2))/(d10*d2^n2*d4^2*d5*d6*d7^2*d8*d9) + 
        rat(4, n2*(4 - d + n2))/(d10^2*d2^n2*d4*d5*d6*d7^2*d8*d9) + 
        rat(4, n2*(4 - d + n2))/(d10*d2^n2*d4^3*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-4 + d - 5*n2, n2*(4 - d + n2)))/(d10*d3^2*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(4 - d - 3*n2, n2*(4 - d + n2)))/(d10*d3*d4*d5*d6*d7*d8*d9^2) + 
        (d2^(-1 - n2)*rat(4 - d - 3*n2, n2*(4 - d + n2)))/(d10*d3*d4*d6*d7*d8^2*d9) + 
        (d2^(-1 - n2)*rat(4 - d - 3*n2, n2*(4 - d + n2)))/(d10*d3^2*d4*d5*d6*d8*d9) + 
        (d2^(-1 - n2)*rat(4 - d - 3*n2, n2*(4 - d + n2)))/(d10*d3*d4*d5*d6*d7^2*d8*d9) + 
        (d2^(-1 - n2)*rat(4 - d - 3*n2, n2*(4 - d + n2)))/(d10*d3*d4*d5^2*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(4 - d - 3*n2, n2*(4 - d + n2)))/(d10*d3*d4*d5^2*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-4 + d - 3*n2, n2*(4 - d + n2)))/(d10^2*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-1 - n2)*rat(-30 + 7*d - 3*n2, 4 - d + n2))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-2 - n2)*rat(-2 - 2*n2, 4 - d + n2))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d1*d2^(-2 - n2)*rat(-2 - 2*n2, 4 - d + n2))/(d10*d3*d4*d5*d6*d7*d8*d9) + 
        (d2^(-2 - n2)*rat(2 + 2*n2, 4 - d + n2))/(d10*d3*d4*d5*d6*d7*d9) + 
        (d2^(-1 - n2)*rat(-4 + d + 3*n2, n2*(4 - d + n2)))/(d10*d3*d4*d5^2*d6*d7*d9) + 
        (d2^(-1 - n2)*rat(-4 + d + 3*n2, n2*(4 - d + n2)))/(d10*d3*d4*d5*d7*d8^2*d9) + 
        (d2^(-1 - n2)*rat(-4 + d + 3*n2, n2*(4 - d + n2)))/(d10*d4*d5*d6*d7^2*d8*d9);

* n1 == 0 && n10 == 1 && n2 == 0 && n3 == 1 && n4 == 1 && n5 == 2 && n6 == 1 && n7 == 1 && n8 == 1 && n9 == 1

        if((count(d1,1)==0) && (count(d2,1)==0)) 
        id,only 1/d3/d4/d5^2/d6/d7/d8/d9/d10 = 1/d3^2/d4/d5/d6/d7/d8/d9/d10;
        
        if((count(d1,1)==0) && (count(d2,1)==0)) 
        id,only 1/d3^2/d4/d5/d6/d7/d8/d9/d10 = PR11d/intBMW;
        
        if((count(d1,1)==0) && (count(d2,1)==0)) 
        id,only 1/d3/d4/d5/d6/d7/d8/d9/d10 = PR11/intBMW;

        #call zeroBMW
        
        endif;
        
        goto endrec;         
        la sortme;
        $irep = 0;
        la endrec;
        
        ModuleOption,minimum,$irep;
        .sort:redBMW-`$repcount++';
        #redefine irep "`$irep'"
#enddo
#endprocedure
*--#] redBMW :


* 
* 
*************************************************************************** 
* 
*   New reduction for 2-loop propagators 
*

* all to 125
*--#[ uniqueTK :
#procedure uniqueT2new(TOPO)
        if(count(int`TOPO',1));        
        if((count(tarC1,-1) > 0) && (count(tarC2,-1) > 0) && (count(tarC3,-1) <= 0) && (count(tarC4,-1) <= 0) && (count(tarC5,-1) > 0))
        Multiply intT2/int`TOPO';

        if((count(tarC1,-1) <= 0) && (count(tarC2,-1) <= 0) && (count(tarC3,-1) > 0) && (count(tarC4,-1) > 0) && (count(tarC5,-1) > 0))
        Multiply replace_(tarC4,tarC1, tarC3,tarC2, tarC2,tarC3, tarC1,tarC4)*intT2/int`TOPO';
        endif;        
#endprocedure
*--#] uniqueTK :

*--#[ uniqueTJ :
#procedure uniqueJnew(TOPO)
        if(count(int`TOPO',1));        
        if((count(tarC1,-1) <= 0) && (count(tarC2,-1) > 0) && (count(tarC3,-1) > 0) && (count(tarC4,-1) <= 0) && (count(tarC5,-1) > 0))
        Multiply intJ/int`TOPO';

        if((count(tarC1,-1) > 0) && (count(tarC2,-1) <= 0) && (count(tarC3,-1) <= 0) && (count(tarC4,-1) > 0) && (count(tarC5,-1) > 0))
        Multiply replace_(tarC3,tarC1, tarC4,tarC2, tarC1,tarC3, tarC2,tarC4)*intJ/int`TOPO';
        endif;        
#endprocedure
*--#] uniqueTJ :

*--#[ uniqueTV :
#procedure uniqueVnew(TOPO)
        if(count(int`TOPO',1));        
        if((count(tarC1,-1) <= 0) && (count(tarC2,-1) > 0) && (count(tarC3,-1) > 0) && (count(tarC4,-1) > 0) && (count(tarC5,-1) > 0))
        Multiply intV/int`TOPO';
        if((count(tarC1,-1) > 0) && (count(tarC2,-1) <= 0) && (count(tarC3,-1) > 0) && (count(tarC4,-1) > 0) && (count(tarC5,-1) > 0))
        Multiply replace_(tarC1,tarC2, tarC2,tarC1, tarC3,tarC4, tarC4,tarC3)*intV/int`TOPO';

        if((count(tarC1,-1) > 0) && (count(tarC2,-1) > 0) && (count(tarC3,-1) <= 0) && (count(tarC4,-1) > 0) && (count(tarC5,-1) > 0))
        Multiply replace_(tarC1,tarC3, tarC3,tarC1, tarC2,tarC4, tarC4,tarC2)*intV/int`TOPO';

        if((count(tarC1,-1) > 0) && (count(tarC2,-1) > 0) && (count(tarC3,-1) > 0) && (count(tarC4,-1) <= 0) && (count(tarC5,-1) > 0))
        Multiply replace_(tarC4,tarC1, tarC1,tarC4, tarC3,tarC2, tarC2,tarC3)*intV/int`TOPO';
        endif;
#endprocedure        
*--#] uniqueTV :


#procedure red01111
#$repcount = 1;        
#do irep=1,1
#$irep = 1;                
* do not need partfrac=1        
#$nopf = 1;                        
if(count(intV,1));

*   n1 != 0 && n5 != 1
id,ifmatch->sortme 1/tarC1^n1?neg_/tarC2^n2?pos_/tarC3^n3?pos_/tarC4^n4?pos_/tarC5^n5?{>1} = 
(n4*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC4^(-1 - n4)*tarC5^(1 - n5))/((-1 + n5)*tarC3^n3) + (tarC1^(-1 - n1)*tarC2^(1 - n2))/(tarC3^n3*tarC4^n4*tarC5^n5) + (tarC1^(-1 - n1)*tarC3^(1 - n3)*(1 - pp/(3*tmm)))/(tarC2^n2*tarC4^n4*tarC5^n5) + (tarC1^(-1 - n1)*tarC4^(1 - n4)*(-1 + pp/(3*tmm)))/(tarC2^n2*tarC3^n3*tarC5^n5) + (tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC4^(1 - n4)*tarC5^(1 - n5)*(-(n2/(-1 + n5)) + (n2*pp)/(2*(-1 + n5)*tmm)))/tarC3^n3 + (tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC5^(1 - n5)*((2*n2*pp)/(-1 + n5) - (n2*pp^2)/(2*(-1 + n5)*tmm)))/(tarC3^n3*tarC4^n4) + (n3*pp*tarC1^(-1 - n1)*tarC3^(-1 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC2^n2*tmm) + ((1 + n1)*pp*tarC1^(-2 - n1)*tarC3^(1 - n3)*tarC5^(1 - n5))/(6*(-1 + n5)*tarC2^n2*tarC4^n4*tmm) + ((1 + n1)*pp*tarC1^(-2 - n1)*tarC2^(1 - n2)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC3^n3*tarC4^n4*tmm) - ((1 + n1)*pp^2*tarC1^(-2 - n1)*tarC5^(1 - n5))/(6*(-1 + n5)*tarC2^n2*tarC3^n3*tarC4^n4*tmm) - (n3*pp*tarC1^(-1 - n1)*tarC3^(-1 - n3)*tarC5^(2 - n5))/(3*(-1 + n5)*tarC2^n2*tarC4^n4*tmm) - ((1 + n1)*pp*tarC1^(-2 - n1)*tarC5^(2 - n5))/(3*(-1 + n5)*tarC2^n2*tarC3^n3*tarC4^n4*tmm) + (tarC1^(-1 - n1)*tarC5^(1 - n5)*((n2 - n4)/(-1 + n5) - (pp*rat(1 + 2*d + n1 - 3*n2 - 6*n4, 1))/(6*(-1 + n5)*tmm)))/(tarC2^n2*tarC3^n3*tarC4^n4);



*   n1 != 0 && n3 != 1
id,ifmatch->sortme 1/tarC1^n1?neg_/tarC2^n2?pos_/tarC3^n3?{>1}/tarC4^n4?pos_/tarC5^n5?pos_ = 
-((n4*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3)*tarC4^(-1 - n4))/((-1 + n3)*tarC5^n5)) + (tarC1^(-2 - n1)*tarC3^(2 - n3)*((1 + n1)/(-1 + n3) - (5*(1 + n1)*pp)/(6*(-1 + n3)*tmm)))/(tarC2^n2*tarC4^n4*tarC5^n5) + (tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC3^(1 - n3)*tarC4^(1 - n4)*(n2/(-1 + n3) - (n2*pp)/(2*(-1 + n3)*tmm)))/tarC5^n5 + (tarC1^(-2 - n1)*tarC3^(1 - n3)*((-2*(1 + n1)*pp)/(-1 + n3) + (5*(1 + n1)*pp^2)/(6*(-1 + n3)*tmm)))/(tarC2^n2*tarC4^n4*tarC5^n5) + (tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC3^(1 - n3)*((-2*n2*pp)/(-1 + n3) + (n2*pp^2)/(2*(-1 + n3)*tmm)))/(tarC4^n4*tarC5^n5) + (n5*pp*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC4^(1 - n4)*tarC5^(-1 - n5))/(3*(-1 + n3)*tarC2^n2*tmm) - (n5*pp*tarC1^(-1 - n1)*tarC3^(2 - n3)*tarC5^(-1 - n5))/(3*(-1 + n3)*tarC2^n2*tarC4^n4*tmm) + (2*(1 + n1)*pp*tarC1^(-2 - n1)*tarC3^(1 - n3)*tarC5^(1 - n5))/(3*(-1 + n3)*tarC2^n2*tarC4^n4*tmm) + (2*pp*tarC1^(-1 - n1)*tarC5^(1 - n5))/(3*tarC2^n2*tarC3^n3*tarC4^n4*tmm) - (2*pp*tarC1^(-1 - n1)*tarC4^(1 - n4))/(3*tarC2^n2*tarC3^n3*tarC5^n5*tmm) - (2*(1 + n1)*pp*tarC1^(-2 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3))/(3*(-1 + n3)*tarC4^n4*tarC5^n5*tmm) + (tarC1^(-1 - n1)*tarC3^(1 - n3)*((-2 - n1 - n2 + n3 + n4)/(-1 + n3) + (pp*rat(5 + 4*d - n1 - 3*n2 - 6*n3 - 6*n4, 1))/(6*(-1 + n3)*tmm)))/(tarC2^n2*tarC4^n4*tarC5^n5);



*   n1 != 0 && n2 != 1
id,ifmatch->sortme 1/tarC1^n1?neg_/tarC2^n2?{>1}/tarC3^n3?pos_/tarC4^n4?pos_/tarC5^n5?pos_ = 
-((n5*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC4^(1 - n4)*tarC5^(-1 - n5))/((-1 + n2)*tarC3^n3)) + (n5*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3)*tarC5^(-1 - n5))/((-1 + n2)*tarC4^n4) + (n4*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC4^(-1 - n4)*tarC5^(1 - n5))/((-1 + n2)*tarC3^n3) - (n3*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(-1 - n3)*tarC5^(1 - n5))/((-1 + n2)*tarC4^n4) - ((1 + n1)*tarC1^(-2 - n1)*tarC2^(1 - n2)*tarC5^(1 - n5))/((-1 + n2)*tarC3^n3*tarC4^n4) + (tarC1^(-1 - n1)*tarC5^(1 - n5))/(tarC2^n2*tarC3^n3*tarC4^n4) - (n4*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3)*tarC4^(-1 - n4))/((-1 + n2)*tarC5^n5) + (n3*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(-1 - n3)*tarC4^(1 - n4))/((-1 + n2)*tarC5^n5) - (tarC1^(-1 - n1)*tarC4^(1 - n4))/(2*tarC2^n2*tarC3^n3*tarC5^n5) + ((1 + n1)*tarC1^(-2 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3))/(2*(-1 + n2)*tarC4^n4*tarC5^n5) - ((1 + n1)*pp*tarC1^(-2 - n1)*tarC2^(1 - n2))/(2*(-1 + n2)*tarC3^n3*tarC4^n4*tarC5^n5) + ((-2 - n1 + n2)*tarC1^(-1 - n1)*tarC2^(1 - n2))/(2*(-1 + n2)*tarC3^n3*tarC4^n4*tarC5^n5) + ((1 + n1)*tarC1^(-2 - n1)*tarC2^(2 - n2))/((-1 + n2)*tarC3^n3*tarC4^n4*tarC5^n5) + (pp*tarC1^(-1 - n1))/(2*tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5);



*   n5 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?pos_/tarC3^n3?pos_/tarC4^n4?pos_/tarC5^n5?{>1} = 
(2*n3*tarC3^(-1 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC1^n1*tarC2^n2*tmm) + (n1*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC2^n2*tarC4^n4*tmm) + (2*n1*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC3^n3*tarC4^n4*tmm) - (n1*pp*tarC1^(-1 - n1)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC2^n2*tarC3^n3*tarC4^n4*tmm) - (2*n3*tarC3^(-1 - n3)*tarC5^(2 - n5))/(3*(-1 + n5)*tarC1^n1*tarC2^n2*tarC4^n4*tmm) - (2*n1*tarC1^(-1 - n1)*tarC5^(2 - n5))/(3*(-1 + n5)*tarC2^n2*tarC3^n3*tarC4^n4*tmm) - tarC4^(1 - n4)/(3*tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5*tmm) + tarC3^(1 - n3)/(3*tarC1^n1*tarC2^n2*tarC4^n4*tarC5^n5*tmm) + (tarC5^(1 - n5)*rat(3 + d - n1 - 3*n5, 1))/(3*(-1 + n5)*tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*tmm);



*   n4 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?pos_/tarC3^n3?pos_/tarC4^n4?{>1}/tarC5^n5?pos_ = 
(tarC2^(-1 - n2)*tarC4^(1 - n4)*(-(n2/(-1 + n4)) + (n2*pp)/(2*(-1 + n4)*tmm)))/(tarC1^n1*tarC3^n3*tarC5^n5) + (n5*tarC3^(1 - n3)*tarC4^(1 - n4)*tarC5^(-1 - n5))/(3*(-1 + n4)*tarC1^n1*tarC2^n2*tmm) - (n5*tarC4^(2 - n4)*tarC5^(-1 - n5))/(3*(-1 + n4)*tarC1^n1*tarC2^n2*tarC3^n3*tmm) + (n3*tarC3^(-1 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n4)*tarC1^n1*tarC2^n2*tmm) + (n1*tarC1^(-1 - n1)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n4)*tarC2^n2*tarC3^n3*tmm) - (n1*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC4^(1 - n4))/(6*(-1 + n4)*tarC2^n2*tarC5^n5*tmm) - (n1*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC4^(1 - n4))/(3*(-1 + n4)*tarC3^n3*tarC5^n5*tmm) + (n1*pp*tarC1^(-1 - n1)*tarC4^(1 - n4))/(6*(-1 + n4)*tarC2^n2*tarC3^n3*tarC5^n5*tmm) - (n3*tarC3^(-1 - n3)*tarC4^(2 - n4))/(3*(-1 + n4)*tarC1^n1*tarC2^n2*tarC5^n5*tmm) - (n2*tarC2^(-1 - n2)*tarC4^(2 - n4))/(2*(-1 + n4)*tarC1^n1*tarC3^n3*tarC5^n5*tmm) + (tarC4^(1 - n4)*rat(6 + 2*d + n1 - 3*n2 - 6*n4, 1))/(6*(-1 + n4)*tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5*tmm);



*   n3 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?pos_/tarC3^n3?{>1}/tarC4^n4?pos_/tarC5^n5?pos_ = 
(tarC1^(-1 - n1)*tarC3^(1 - n3)*(-(n1/(-1 + n3)) + (2*n1*pp)/(3*(-1 + n3)*tmm)))/(tarC2^n2*tarC4^n4*tarC5^n5) + (2*n5*tarC3^(1 - n3)*tarC4^(1 - n4)*tarC5^(-1 - n5))/(3*(-1 + n3)*tarC1^n1*tarC2^n2*tmm) - (2*n5*tarC3^(2 - n3)*tarC5^(-1 - n5))/(3*(-1 + n3)*tarC1^n1*tarC2^n2*tarC4^n4*tmm) + (n1*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC5^(1 - n5))/(3*(-1 + n3)*tarC2^n2*tarC4^n4*tmm) + tarC5^(1 - n5)/(3*tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*tmm) - tarC4^(1 - n4)/(3*tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5*tmm) - (n1*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3))/(3*(-1 + n3)*tarC4^n4*tarC5^n5*tmm) - (2*n1*tarC1^(-1 - n1)*tarC3^(2 - n3))/(3*(-1 + n3)*tarC2^n2*tarC4^n4*tarC5^n5*tmm) + (tarC3^(1 - n3)*rat(3 + d - n1 - 3*n3, 1))/(3*(-1 + n3)*tarC1^n1*tarC2^n2*tarC4^n4*tarC5^n5*tmm);



*   n1 == 0 && n2 != 1 && n2 != 2 && n3 == 1 && n4 == 1 && n5 == 1
        if((count(tarC1,-1)) == 0) id,only,ifmatch->dopartfrac 1/tarC2^n2?{>2}/tarC3/tarC4/tarC5 =
(1/(2*pp) + 1/(2*[pp-4*tmm]))/(tarC2^n2*tarC3*tarC5) + (tarC2^(1 - n2)*(-1/(2*(-1 + n2)*pp) - 1/(2*(-1 + n2)*[pp-4*tmm])))/(tarC3^2*tarC4) + (tarC2^(2 - n2)*(-1/(4*(-1 + n2)*pp) + 1/(4*(-1 + n2)*[pp-4*tmm])))/(tarC3^2*tarC4*tarC5) + (tarC2^(1 - n2)*(1/(4*(-1 + n2)*pp) + 3/(4*(-1 + n2)*[pp-4*tmm])))/(tarC3^2*tarC5) + (tarC2^(2 - n2)*(1/((-2 + n2)*(-1 + n2)*pp) - 1/((-2 + n2)*(-1 + n2)*[pp-4*tmm])))/(tarC3^3*tarC4) + (tarC2^(2 - n2)*(1/(2*(-2 + n2)*(-1 + n2)*pp) - 1/(2*(-2 + n2)*(-1 + n2)*[pp-4*tmm])))/(tarC3^2*tarC5^2) + (tarC2^(2 - n2)*(-1/(2*(-2 + n2)*(-1 + n2)*pp) + 1/(2*(-2 + n2)*(-1 + n2)*[pp-4*tmm])))/(tarC3^2*tarC4^2) + (tarC2^(2 - n2)*(-1/(2*(-2 + n2)*(-1 + n2)*pp) + 1/(2*(-2 + n2)*(-1 + n2)*[pp-4*tmm])))/(tarC3*tarC4*tarC5^2) + (tarC2^(2 - n2)*(-(1/((-2 + n2)*(-1 + n2)*pp)) + 1/((-2 + n2)*(-1 + n2)*[pp-4*tmm])))/(tarC3^3*tarC5) + (tarC2^(2 - n2)*((3 - n2)/(2*(-2 + n2)*(-1 + n2)*pp) + (-3 + n2)/(2*(-2 + n2)*(-1 + n2)*[pp-4*tmm])))/(tarC3*tarC4^2*tarC5) + (tarC2^(1 - n2)*((3 - n2)/(2*(-1 + n2)*pp) + (-3 + n2 - 2*rat(-2 + d - n2, 1))/(2*(-1 + n2)*[pp-4*tmm])))/(tarC3*tarC4*tarC5);



*   n2 == 1 && n3 == 1 && n4 == 1 && n5 == 1 && n1 != 0
id,only,ifmatch->sortme 1/tarC1^n1?neg_/tarC2/tarC3/tarC4/tarC5 = 
(2*pp*tarC1^(-1 - n1)*rat(1, -3 + d - n1))/(tarC2*tarC3^2*tarC4) - rat(1, -3 + d - n1)/(tarC1^n1*tarC2*tarC3^2*tarC4) + (2*(1 + n1)*pp*tarC1^(-2 - n1)*rat(1, -3 + d - n1))/(tarC2*tarC3*tarC4) - (n1*tarC1^(-1 - n1)*rat(1, -3 + d - n1))/(tarC2*tarC3*tarC4) + (pp*tarC1^(-1 - n1)*rat(1, -3 + d - n1))/(tarC2*tarC3*tarC5^2) - (2*rat(1, -3 + d - n1))/(tarC1^n1*tarC2*tarC3*tarC5^2) - (pp*tarC1^(-1 - n1)*rat(1, -3 + d - n1))/(tarC2*tarC4*tarC5^2) + (2*rat(1, -3 + d - n1))/(tarC1^n1*tarC2*tarC4*tarC5^2) - (2*pp*tarC1^(-1 - n1)*rat(1, -3 + d - n1))/(tarC2*tarC3^2*tarC5) + rat(1, -3 + d - n1)/(tarC1^n1*tarC2*tarC3^2*tarC5) + (2*n1*tarC1^(-1 - n1)*rat(1, -3 + d - n1))/(tarC2*tarC4*tarC5) - (2*(1 + n1)*pp*tarC1^(-2 - n1)*rat(1, -3 + d - n1))/(tarC3*tarC4*tarC5) + (n1*tarC1^(-1 - n1)*rat(1, -3 + d - n1))/(tarC3*tarC4*tarC5) - (3*tarC1^(-1 - n1)*tmm*rat(1, -3 + d - n1))/(tarC3*tarC4^2*tarC5) + (tarC1^(-1 - n1)*((-3*pp*rat(1, -3 + d - n1))/2 + 3*tmm*rat(1, -3 + d - n1)))/(tarC2^2*tarC3*tarC5) + (tarC1^(-2 - n1)*((-5*(1 + n1)*pp*rat(1, -3 + d - n1))/2 + 3*(1 + n1)*tmm*rat(1, -3 + d - n1)))/(tarC2*tarC4*tarC5) + (tarC1^(-1 - n1)*((3*pp^2*rat(1, -3 + d - n1))/2 - 6*pp*tmm*rat(1, -3 + d - n1)))/(tarC2^2*tarC3*tarC4*tarC5) + (tarC1^(-2 - n1)*((5*(1 + n1)*pp^2*rat(1, -3 + d - n1))/2 - 6*(1 + n1)*pp*tmm*rat(1, -3 + d - n1)))/(tarC2*tarC3*tarC4*tarC5) + (pp*tarC1^(-1 - n1)*rat(-16 + 4*d - 5*n1, -3 + d - n1))/(2*tarC2*tarC3*tarC4*tarC5);



*   n1 == 0 && n2 == 2 && n3 == 1 && n4 == 1 && n5 == 1
if((count(tarC1,-1)) == 0) id,only,ifmatch->dopartfrac 1/tarC1^n1?neg0_/tarC2^2/tarC3/tarC4/tarC5 = 
rat(1, -3 + d)/(2*tarC2^2*tarC3^2*tarC4) + rat(1, -3 + d)/(2*tarC2^2*tarC4*tarC5^2) - rat(1, -3 + d)/(2*tarC2*tarC3*tarC5^2*[pp-4*tmm]) + rat(1, -3 + d)/(2*tarC2*tarC4*tarC5^2*[pp-4*tmm]) - (2*tmm*rat(1, -3 + d))/(tarC2*tarC3^3*tarC4*[pp-4*tmm]) - (2*tmm*rat(1, -3 + d))/(tarC2*tarC3*tarC5^3*[pp-4*tmm]) + (2*tmm*rat(1, -3 + d))/(tarC2*tarC4*tarC5^3*[pp-4*tmm]) + (2*tmm*rat(1, -3 + d))/(tarC2*tarC3^3*tarC5*[pp-4*tmm]) + ((-3*tmm*rat(1, -3 + d))/(4*pp) + (3*tmm*rat(1, -3 + d))/(4*[pp-4*tmm]))/(tarC3*tarC4^2*tarC5^2) + ((-3*tmm*rat(1, -3 + d))/(4*pp) + (3*tmm*rat(1, -3 + d))/(4*[pp-4*tmm]))/(tarC3^2*tarC4*tarC5^2) + ((-3*tmm*rat(1, -3 + d))/(4*pp) + (3*tmm*rat(1, -3 + d))/(4*[pp-4*tmm]))/(tarC3^2*tarC4^2*tarC5) + (-rat(1, -3 + d)/2 + (3*tmm*rat(1, -3 + d))/(4*pp) + (3*tmm*rat(1, -3 + d))/(4*[pp-4*tmm]))/(tarC2^2*tarC3*tarC5^2) + (-rat(1, -3 + d)/2 + (3*tmm*rat(1, -3 + d))/(4*pp) + (3*tmm*rat(1, -3 + d))/(4*[pp-4*tmm]))/(tarC2^2*tarC3^2*tarC5) + ((3*tmm*rat(1, -3 + d))/(4*pp) + (9*tmm*rat(1, -3 + d))/(4*[pp-4*tmm]))/(tarC2*tarC3^2*tarC5^2) - rat(-3 + d, 1)/(tarC2*tarC3*tarC4*tarC5*[pp-4*tmm]) + rat(-7 + 2*d, -3 + d)/(2*tarC2*tarC3^2*tarC4*[pp-4*tmm]) - rat(-7 + 2*d, -3 + d)/(2*tarC2*tarC3^2*tarC5*[pp-4*tmm]);


endif;

* End recursion
goto endrec;         
la dopartfrac;
$nopf = 0;
la sortme;
$irep = 0;
la endrec;
        
ModuleOption,minimum,$irep,$nopf;
.sort:red2l-V-`$repcount++';
#if `$nopf'==0
        #call partfrac        
#endif                
#redefine irep "`$irep'"
#enddo
#endprocedure



#procedure red01101
#$repcount = 1;        
#do irep=1,1
#$irep = 1;
* do not need partfrac=1        
#$nopf = 1;                        
if(count(intJ,1));

*   n4 != 0 && n5 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?pos_/tarC3^n3?pos_/tarC4^n4?neg_/tarC5^n5?{>1} =
(n1*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC4^(-1 - n4)*tarC5^(1 - n5))/((-1 + n5)*tarC2^n2) + (tarC3^(1 - n3)*tarC4^(-1 - n4))/(tarC1^n1*tarC2^n2*tarC5^n5) + (2*n3*tarC3^(-1 - n3)*tarC4^(-1 - n4)*tarC5^(1 - n5)*tmm)/((-1 + n5)*tarC1^n1*tarC2^n2) + (tarC4^(-1 - n4)*tmm)/(tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5) + (tarC1^(-1 - n1)*tarC4^(-1 - n4)*tarC5^(1 - n5)*(-((n1*pp)/(-1 + n5)) + (2*n1*tmm)/(-1 + n5)))/(tarC2^n2*tarC3^n3) - (tarC4^(-1 - n4)*tarC5^(1 - n5)*rat(1 + d - n1 - 2*n3 - n5, 1))/((-1 + n5)*tarC1^n1*tarC2^n2*tarC3^n3);



*   n1 != 0 && n5 != 1
id,ifmatch->sortme 1/tarC1^n1?neg_/tarC2^n2?pos_/tarC3^n3?pos_/tarC4^n4?neg0_/tarC5^n5?{>1} =
(n4*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC4^(-1 - n4)*tarC5^(1 - n5))/((-1 + n5)*tarC3^n3) + (tarC1^(-1 - n1)*tarC2^(1 - n2))/(tarC3^n3*tarC4^n4*tarC5^n5) + (2*n2*tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC5^(1 - n5)*tmm)/((-1 + n5)*tarC3^n3*tarC4^n4) + (tarC1^(-1 - n1)*tmm)/(tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5) + (tarC1^(-1 - n1)*tarC4^(-1 - n4)*tarC5^(1 - n5)*(-((n4*pp)/(-1 + n5)) + (2*n4*tmm)/(-1 + n5)))/(tarC2^n2*tarC3^n3) - (tarC1^(-1 - n1)*tarC5^(1 - n5)*rat(1 + d - 2*n2 - n4 - n5, 1))/((-1 + n5)*tarC2^n2*tarC3^n3*tarC4^n4);



*   n3 != 1 && n4 != 0
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?pos_/tarC3^n3?{>1}/tarC4^n4?neg_/tarC5^n5?pos_ =
(n1*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC4^(-1 - n4)*tarC5^(1 - n5))/((-1 + n3)*tarC2^n2) + (tarC4^(-1 - n4)*tarC5^(1 - n5))/(tarC1^n1*tarC2^n2*tarC3^n3) - (n1*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3)*tarC4^(-1 - n4))/((-1 + n3)*tarC5^n5) + (2*n5*tarC3^(1 - n3)*tarC4^(-1 - n4)*tarC5^(-1 - n5)*tmm)/((-1 + n3)*tarC1^n1*tarC2^n2) + (n1*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC4^(-1 - n4)*tmm)/((-1 + n3)*tarC2^n2*tarC5^n5) + (tarC4^(-1 - n4)*tmm)/(tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5) - (tarC3^(1 - n3)*tarC4^(-1 - n4)*rat(1 + d - n1 - n3 - 2*n5, 1))/((-1 + n3)*tarC1^n1*tarC2^n2*tarC5^n5);



*   n1 != 0 && n3 != 1
id,ifmatch->sortme 1/tarC1^n1?neg_/tarC2^n2?pos_/tarC3^n3?{>1}/tarC4^n4?neg0_/tarC5^n5?pos_ =
-((n4*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3)*tarC4^(-1 - n4))/((-1 + n3)*tarC5^n5)) + (tarC1^(-1 - n1)*(pp - 2*tmm))/(tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5) - (2*n5*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC5^(-1 - n5)*tmm)/((-1 + n3)*tarC2^n2*tarC4^n4) - (2*n2*tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC3^(1 - n3)*tmm)/((-1 + n3)*tarC4^n4*tarC5^n5) - (2*(1 + n1)*tarC1^(-2 - n1)*tarC3^(1 - n3)*tmm)/((-1 + n3)*tarC2^n2*tarC4^n4*tarC5^n5) + (tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC4^(-1 - n4)*((n4*pp)/(-1 + n3) - (2*n4*tmm)/(-1 + n3)))/(tarC2^n2*tarC5^n5) + (tarC1^(-1 - n1)*tarC3^(1 - n3)*rat(-1 + 2*d - 2*n1 - 2*n2 - n3 - n4 - 2*n5, 1))/((-1 + n3)*tarC2^n2*tarC4^n4*tarC5^n5);



*   n2 != 1 && n4 != 0
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?{>1}/tarC3^n3?pos_/tarC4^n4?neg_/tarC5^n5?pos_ =
-((n1*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3)*tarC4^(-1 - n4))/((-1 + n2)*tarC5^n5)) + (tarC4^(-1 - n4)*(pp - 2*tmm))/(tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5) - (2*n5*tarC2^(1 - n2)*tarC4^(-1 - n4)*tarC5^(-1 - n5)*tmm)/((-1 + n2)*tarC1^n1*tarC3^n3) - (2*(1 + n4)*tarC2^(1 - n2)*tarC4^(-2 - n4)*tmm)/((-1 + n2)*tarC1^n1*tarC3^n3*tarC5^n5) - (2*n3*tarC2^(1 - n2)*tarC3^(-1 - n3)*tarC4^(-1 - n4)*tmm)/((-1 + n2)*tarC1^n1*tarC5^n5) + (tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC4^(-1 - n4)*((n1*pp)/(-1 + n2) - (2*n1*tmm)/(-1 + n2)))/(tarC3^n3*tarC5^n5) + (tarC2^(1 - n2)*tarC4^(-1 - n4)*rat(-1 + 2*d - n1 - n2 - 2*n3 - 2*n4 - 2*n5, 1))/((-1 + n2)*tarC1^n1*tarC3^n3*tarC5^n5);



*   n1 != 0 && n2 != 1
id,ifmatch->sortme 1/tarC1^n1?neg_/tarC2^n2?{>1}/tarC3^n3?pos_/tarC4^n4?neg0_/tarC5^n5?pos_ =
(n4*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC4^(-1 - n4)*tarC5^(1 - n5))/((-1 + n2)*tarC3^n3) + (tarC1^(-1 - n1)*tarC5^(1 - n5))/(tarC2^n2*tarC3^n3*tarC4^n4) - (n4*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3)*tarC4^(-1 - n4))/((-1 + n2)*tarC5^n5) + (2*n5*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC5^(-1 - n5)*tmm)/((-1 + n2)*tarC3^n3*tarC4^n4) + (n4*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC4^(-1 - n4)*tmm)/((-1 + n2)*tarC3^n3*tarC5^n5) + (tarC1^(-1 - n1)*tmm)/(tarC2^n2*tarC3^n3*tarC4^n4*tarC5^n5) - (tarC1^(-1 - n1)*tarC2^(1 - n2)*rat(1 + d - n2 - n4 - 2*n5, 1))/((-1 + n2)*tarC3^n3*tarC4^n4*tarC5^n5);



*   n1 == 0 && n4 == 0 && n5 != 1 && n5 != 2
if((count(tarC1,-1) == 0) && (count(tarC4,-1) == 0)) id,ifmatch->sortme 1/tarC2^n2?pos_/tarC3^n3?pos_/tarC5^n5?{>2} =
-((n2*tarC2^(-1 - n2)*tarC5^(1 - n5))/((-1 + n5)*tarC3^n3)) - (n2*n3*tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC5^(2 - n5))/((-2 + n5)*(-1 + n5)) + (tarC3^(-1 - n3)*tarC5^(1 - n5)*((-3*n3)/(2*(-1 + n5)) + (n3*pp)/(2*(-1 + n5)*tmm)))/tarC2^n2 - (n3*tarC2^(1 - n2)*tarC3^(-1 - n3)*tarC5^(1 - n5))/(2*(-1 + n5)*tmm) + (tarC5^(1 - n5)*rat(2 + 2*d - 2*n2 - n3 - 2*n5, 1))/(2*(-1 + n5)*tarC2^n2*tarC3^n3*tmm) + (n3*tarC3^(-1 - n3)*tarC5^(2 - n5)*rat(2 + d - 2*n2 - n5, 1))/(2*(-2 + n5)*(-1 + n5)*tarC2^n2*tmm);



*   n1 == 0 && n4 == 0 && n5 == 2 && n3 != 1
if((count(tarC1,-1) == 0) && (count(tarC4,-1) == 0)) id,only,ifmatch->dopartfrac 1/tarC2^n2?pos_/tarC3^n3?{>1}/tarC5^2 =
tarC2^(1 - n2)/(tarC3^n3*tarC5^2*[pp-3*tmm]) + (2*n2*tarC2^(-1 - n2)*tarC3^(1 - n3)*tmm)/((-1 + n3)*tarC5^2*[pp-3*tmm]) + (2*n3*tarC3^(-1 - n3)*tmm)/(tarC2^n2*tarC5*[pp-3*tmm]) + (2*n2*tarC2^(-1 - n2)*tmm)/(tarC3^n3*tarC5*[pp-3*tmm]) - rat(-1 + 2*d - 2*n2 - 2*n3, 1)/(tarC2^n2*tarC3^n3*tarC5*[pp-3*tmm]) - (tarC3^(1 - n3)*rat(1 + d - 2*n2 - n3, 1))/((-1 + n3)*tarC2^n2*tarC5^2*[pp-3*tmm]);



*   n1 == 0 && n3 != 1 && n3 != 2 && n4 == 0 && n5 == 1
if((count(tarC1,-1) == 0) && (count(tarC4,-1) == 0)) id,only,ifmatch->sortme 1/tarC2^n2?pos_/tarC3^n3?{>2}/tarC5 =
(n2*(1 + n2)*tarC2^(-2 - n2)*tarC3^(2 - n3))/((-2 + n3)*(-1 + n3)*tarC5) + (n2*tarC2^(-1 - n2)*tarC3^(2 - n3)*rat(2 - d + 2*n2, 1))/(2*(-2 + n3)*(-1 + n3)*tarC5*tmm) + (tarC3^(1 - n3)*rat(2 + d - 2*n3, 1))/(2*(-1 + n3)*tarC2^n2*tarC5*tmm);



*   n2 == 1 && n3 == 1 && n5 == 1 && n4 != 0
id,only,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2/tarC3/tarC4^n4?neg_/tarC5 =
(n4*tarC4^(-1 - n4)*rat(1, -3 + d - n4))/(tarC1^n1*tarC2*tarC3) + rat(1, -3 + d - n4)/(tarC1^n1*tarC2^2*tarC3*tarC4^n4) + ((-1 + n1 - n4)*tarC4^(-1 - n4)*rat(1, -3 + d - n4))/(tarC1^n1*tarC2*tarC5) + (2*tarC4^(-1 - n4)*tmm*rat(1, -3 + d - n4))/(tarC1^n1*tarC2*tarC5^2) + (2*tarC4^(-1 - n4)*tmm*rat(1, -3 + d - n4))/(tarC1^n1*tarC3*tarC5^2) + (2*(1 + n4)*tarC1^(1 - n1)*tarC4^(-2 - n4)*tmm*rat(1, -3 + d - n4))/(tarC2*tarC3*tarC5) - ((1 + n4)*pp*tarC4^(-2 - n4)*tmm*rat(1, -3 + d - n4))/(tarC1^n1*tarC2*tarC3*tarC5) + (n1*tarC1^(-1 - n1)*tarC4^(-1 - n4)*tmm*rat(1, -3 + d - n4))/(tarC2*tarC5) + (tarC4^(-1 - n4)*(-(pp*rat(1, -3 + d - n4)) + 2*tmm*rat(1, -3 + d - n4)))/(tarC1^n1*tarC2^2*tarC3) + (tarC4^(-2 - n4)*((1 + n4)*pp*rat(1, -3 + d - n4) - 2*(1 + n4)*tmm*rat(1, -3 + d - n4)))/(tarC1^n1*tarC2*tarC5) + (tarC4^(-2 - n4)*(-((1 + n4)*pp*rat(1, -3 + d - n4)) + 2*(1 + n4)*tmm*rat(1, -3 + d - n4)))/(tarC1^n1*tarC2*tarC3) + (tarC4^(-1 - n4)*(2*pp*tmm*rat(1, -3 + d - n4) - 2*tmm^2*rat(1, -3 + d - n4)))/(tarC1^n1*tarC2*tarC3^2*tarC5) + (tarC4^(-1 - n4)*(-2*pp*tmm*rat(1, -3 + d - n4) + 2*tmm^2*rat(1, -3 + d - n4)))/(tarC1^n1*tarC2*tarC3*tarC5^2) + (tarC1^(-1 - n1)*tarC4^(-1 - n4)*(-(n1*pp*tmm*rat(1, -3 + d - n4)) - 2*n1*tmm^2*rat(1, -3 + d - n4)))/(tarC2*tarC3*tarC5) - (tarC1^(1 - n1)*tarC4^(-1 - n4)*rat(-6 + 2*d - n1 - 2*n4, -3 + d - n4))/(tarC2*tarC3*tarC5) + (tarC4^(-1 - n4)*(-((-1 + n1 - n4)*tmm*rat(1, -3 + d - n4)) + pp*rat(-3 + d - n1 - n4, -3 + d - n4)))/(tarC1^n1*tarC2*tarC3*tarC5);



*   n1 == 0 && n3 == 1 && n4 == 0 && n5 == 2 && n2 != 1
if((count(tarC1,-1) == 0) && (count(tarC4,-1) == 0)) id,ifmatch->dopartfrac 1/tarC2^n2?{>1}/tarC3/tarC5^2 =
(n2*tarC2^(-1 - n2))/(tarC3*tarC5) - (2*tarC2^(1 - n2))/((-1 + n2)*tarC5^3*[pp-tmm]) - (2*tarC2^(2 - n2))/((-1 + n2)*tarC3*tarC5^3*[pp-tmm]) + (tarC2^(2 - n2)*rat(-4 + d, 1))/(2*(-1 + n2)*tarC3*tarC5^2*tmm*[pp-tmm]) - (tarC2^(1 - n2)*rat((-3 + d)*(3*d - 4*n2), 1))/(2*(-1 + n2)*tarC3*tarC5*tmm*[pp-tmm]) + (-rat(d - 2*n2, 1)/(2*tmm) + (-rat(d - 2*n2, 1) + rat(-6 + 3*d - 2*n2, 1))/(2*[pp-tmm]))/(tarC2^n2*tarC3*tarC5) + (tarC2^(1 - n2)*rat(-1 + d - n2, 1))/((-1 + n2)*tarC3*tarC5^2*[pp-tmm]) + (2*tarC2^(1 - n2)*rat(-1 + d - n2, 1))/((-1 + n2)*tarC3^2*tarC5*[pp-tmm]) + (tarC2^(1 - n2)*rat(-1 + d - n2, 1))/((-1 + n2)*tarC5^2*tmm*[pp-tmm]);



*   n1 == 0 && n3 == 2 && n4 == 0 && n5 == 1 && n2 != 1
if((count(tarC1,-1) == 0) && (count(tarC4,-1) == 0)) id,only,ifmatch->dopartfrac 1/tarC2^n2?{>1}/tarC3^2/tarC5 =
(n2*tarC2^(-1 - n2))/(tarC3*tarC5) + (-1/(2*tmm) + 1/(2*[pp-tmm]))/(tarC2^n2*tarC3^2) + (tarC2^(1 - n2)*(1/((-1 + n2)*tmm) + 1/((1 - n2)*[pp-tmm])))/tarC3^3 - (tarC1*tarC2^(1 - n2))/((-1 + n2)*tarC3^3*tmm*[pp-tmm]) + tarC4/(2*tarC2^n2*tarC3^2*tmm*[pp-tmm]) + ((rat(-2 + 3*d - 4*n2, 1) - rat(d - 2*n2, 1))/(2*[pp-tmm]) - rat(d - 2*n2, 1)/(2*tmm))/(tarC2^n2*tarC3*tarC5) + (2*tarC2^(1 - n2)*rat(-1 + d - n2, 1))/((-1 + n2)*tarC3*tarC5^2*[pp-tmm]) + (tarC2^(1 - n2)*rat(-1 + d - n2, 1))/((-1 + n2)*tarC3^2*tarC5*[pp-tmm]) + (tarC2^(1 - n2)*rat(-1 + d - n2, 1))/(2*(-1 + n2)*tarC3^2*tmm*[pp-tmm]) - (tarC2^(1 - n2)*rat((-4 + 3*d - 2*n2)*(-1 + d - n2), 1))/(2*(-1 + n2)*tarC3*tarC5*tmm*[pp-tmm]);



*   n1 == 0 && n2 != 1 && n2 != 2 && n3 == 1 && n4 == 0 && n5 == 1
if((count(tarC1,-1) == 0) && (count(tarC4,-1) == 0)) id,only,ifmatch->dopartfrac 1/tarC2^n2?{>2}/tarC3/tarC5 =
(-1/(2*tmm) + 1/(4*[pp-tmm]) - 7/(4*[pp-9*tmm]))/(tarC2^n2*tarC3) + (tarC2^(1 - n2)*(-1/(4*(-1 + n2)*[pp-tmm]) + 1/(4*(-1 + n2)*[pp-9*tmm])))/tarC5^2 + (tarC2^(1 - n2)*(1/(2*(-1 + n2)*tmm) - 1/(8*(-1 + n2)*[pp-tmm]) + 13/(8*(-1 + n2)*[pp-9*tmm])))/tarC3^2 + (tarC2^(3 - n2)*(-1/(4*(-2 + n2)*(-1 + n2)*[pp-tmm]) + 1/(4*(-2 + n2)*(-1 + n2)*[pp-9*tmm])))/(tarC3^2*tarC5^2) + (tarC4*(1/(4*tmm*[pp-tmm]) + 1/(4*tmm*[pp-9*tmm])))/(tarC2^n2*tarC3) + (tarC1*tarC2^(1 - n2)*(-3/(8*(-1 + n2)*tmm*[pp-tmm]) - 1/(8*(-1 + n2)*tmm*[pp-9*tmm])))/tarC3^2 + (tarC2^(1 - n2)*((-rat(29*d - 32*n2, 1) + 2*rat(7 + 7*d - 11*n2, 1) - rat(2 + d - 2*n2, 1))/(16*(-1 + n2)*[pp-tmm]) + rat(2 + d - 2*n2, 1)/(2*(-1 + n2)*tmm) + (rat(29*d - 32*n2, 1) - 18*rat(7 + 7*d - 11*n2, 1) + 81*rat(2 + d - 2*n2, 1))/(16*(-1 + n2)*[pp-9*tmm])))/(tarC3*tarC5) + (tarC2^(1 - n2)*((rat(2 + 3*d - 5*n2, 1) - 9*rat(d - n2, 1))/(16*(-1 + n2)*tmm*[pp-9*tmm]) + (-rat(2 + 3*d - 5*n2, 1) + rat(d - n2, 1))/(16*(-1 + n2)*tmm*[pp-tmm])))/tarC3 + (tarC2^(2 - n2)*((rat(-3 + 5*d - 4*n2, 1) - 9*rat(d - n2, 1))/(4*(-2 + n2)*(-1 + n2)*[pp-9*tmm]) + (-rat(-3 + 5*d - 4*n2, 1) + rat(d - n2, 1))/(4*(-2 + n2)*(-1 + n2)*[pp-tmm])))/(tarC3*tarC5^2) + (tarC2^(2 - n2)*((rat(-3 + 5*d - 4*n2, 1) - 9*rat(d - n2, 1))/(4*(-2 + n2)*(-1 + n2)*[pp-9*tmm]) + (-rat(-3 + 5*d - 4*n2, 1) + rat(d - n2, 1))/(4*(-2 + n2)*(-1 + n2)*[pp-tmm])))/(tarC3^2*tarC5) + (tarC2^(2 - n2)*((3*rat((-2 + 3*d - 2*n2)*(d - n2), 1))/(8*(-2 + n2)*(-1 + n2)*tmm*[pp-tmm]) + rat((-2 + 3*d - 2*n2)*(d - n2), 1)/(8*(-2 + n2)*(-1 + n2)*tmm*[pp-9*tmm])))/(tarC3*tarC5);



*   n2 == 1 && n3 == 1 && n4 == 0 && n5 == 1 && n1 != 0
if(count(tarC4,-1) == 0) id,only,ifmatch->sortme 1/tarC1^n1?neg_/tarC2/tarC3/tarC5 =
(2*tarC1^(-1 - n1)*tarC4*rat(1, -8 + 3*d - 2*n1))/(tarC2^2*tarC3) + (2*(1 + n1)*tarC1^(-1 - n1)*rat(1, -8 + 3*d - 2*n1))/(tarC2*tarC5) + (2*(1 + n1)*tarC1^(-2 - n1)*tmm*rat(1, -8 + 3*d - 2*n1))/(tarC2*tarC5) + (tarC1^(-1 - n1)*(-2*pp*rat(1, -8 + 3*d - 2*n1) + 4*tmm*rat(1, -8 + 3*d - 2*n1)))/(tarC2^2*tarC3) + (tarC1^(-1 - n1)*(4*pp*tmm*rat(1, -8 + 3*d - 2*n1) - 4*tmm^2*rat(1, -8 + 3*d - 2*n1)))/(tarC2*tarC3^2*tarC5) - rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1))/(tarC1^n1*tarC2*tarC3^2) + ((1 + n1)*tarC1^(-1 - n1)*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1)))/(tarC3*tarC5) - (2*(1 + n1)*tarC1^(-2 - n1)*tarC4*tmm*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1)))/(tarC2*tarC3*tarC5) + (tarC1^(-1 - n1)*(pp*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1)) - 2*tmm*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1))))/(tarC2*tarC3^2) + (tarC1^(-2 - n1)*((1 + n1)*pp*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1)) - 2*(1 + n1)*tmm*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1))))/(tarC2*tarC3) + (tarC1^(-2 - n1)*(-((1 + n1)*pp*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1))) + 2*(1 + n1)*tmm*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1))))/(tarC3*tarC5) + (tarC1^(-1 - n1)*(-2*pp*tmm*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1)) + 2*tmm^2*rat(-2 + d, (-8 + 3*d - 2*n1)*(-3 + d - n1))))/(tarC2^2*tarC3*tarC5) + (2*tarC1^(-1 - n1)*tmm*rat(-4 + d - 2*n1, (-8 + 3*d - 2*n1)*(-3 + d - n1)))/(tarC2*tarC5^2) + (2*tarC1^(-1 - n1)*tmm*rat(-4 + d - 2*n1, (-8 + 3*d - 2*n1)*(-3 + d - n1)))/(tarC3*tarC5^2) + (tarC1^(-1 - n1)*(-2*pp*tmm*rat(-4 + d - 2*n1, (-8 + 3*d - 2*n1)*(-3 + d - n1)) + 2*tmm^2*rat(-4 + d - 2*n1, (-8 + 3*d - 2*n1)*(-3 + d - n1))))/(tarC2*tarC3*tarC5^2) + (tarC1^(-2 - n1)*(-((1 + n1)*pp*tmm*rat(-4 + d - 2*n1, (-8 + 3*d - 2*n1)*(-3 + d - n1))) - 4*(1 + n1)*tmm^2*rat(-3 + d - n1, (-8 + 3*d - 2*n1)*(-3 + d - n1))))/(tarC2*tarC3*tarC5) + (tarC1^(-1 - n1)*(-((1 + n1)*tmm*rat(-8 + 3*d - 2*n1, (-8 + 3*d - 2*n1)*(-3 + d - n1))) + pp*rat((-4 + d - 2*n1)*(-3 + d - n1), (-8 + 3*d - 2*n1)*(-3 + d - n1))))/(tarC2*tarC3*tarC5) - (tarC1^(-1 - n1)*rat(-6 + 2*d - 4*n1 + d*n1, (-8 + 3*d - 2*n1)*(-3 + d - n1)))/(tarC2*tarC3);



*   n1 == 0 && n2 == 1 && n3 == 1 && n4 == 0 && n5 == 2
if((count(tarC1,-1) == 0) && (count(tarC4,-1) == 0)) id,only,ifmatch->sortme 1/tarC2/tarC3/tarC5^2 = 1/(tarC2*tarC3^2*tarC5);



*   n1 == 0 && n2 == 1 && n3 == 2 && n4 == 0 && n5 == 1
if((count(tarC1,-1) == 0) && (count(tarC4,-1) == 0)) id,only,ifmatch->sortme 1/tarC2/tarC3^2/tarC5 = 1/(tarC2^2*tarC3*tarC5);

endif;

* End recursion
goto endrec;
la dopartfrac;
$nopf = 0;
la sortme;
$irep = 0;
la endrec;
        
ModuleOption,minimum,$irep,$nopf;
.sort:red2l-J-`$repcount++';
#if `$nopf'==0
        #call partfrac        
#endif                
#redefine irep "`$irep'"
#enddo
#endprocedure




#procedure red00111
#$repcount = 1;        
#do irep=1,1
#$irep = 1;                
if(count(intT2,1));

*   n2 != 0 && n5 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?neg_/tarC3^n3?pos_/tarC4^n4?pos_/tarC5^n5?{>1} =
(n3*tarC1^(1 - n1)*tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC5^(1 - n5))/((-1 + n5)*tarC4^n4) + (tarC1^(1 - n1)*tarC2^(-1 - n2))/(tarC3^n3*tarC4^n4*tarC5^n5) + (tarC1^(-1 - n1)*tarC5^(1 - n5)*(n1/(-1 + n5) - (n1*pp)/(3*(-1 + n5)*tmm)))/(tarC2^n2*tarC3^n3*tarC4^n4) + (tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC5^(2 - n5)*(-(n1/(-1 + n5)) + (n1*pp)/(3*(-1 + n5)*tmm)))/(tarC3^n3*tarC4^n4) + (tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC3^(1 - n3)*tarC5^(1 - n5)*(-n1/(2*(-1 + n5)) + (n1*pp)/(3*(-1 + n5)*tmm)))/tarC4^n4 + (tarC2^(-2 - n2)*tarC5^(2 - n5)*((1 + n2)/(-1 + n5) - (2*(1 + n2)*pp)/(3*(-1 + n5)*tmm)))/(tarC1^n1*tarC3^n3*tarC4^n4) + (tarC2^(-2 - n2)*tarC4^(1 - n4)*tarC5^(1 - n5)*(-(1 + n2)/(2*(-1 + n5)) + ((1 + n2)*pp)/(3*(-1 + n5)*tmm)))/(tarC1^n1*tarC3^n3) + (tarC1^(1 - n1)*tarC2^(-2 - n2)*tarC5^(1 - n5)*(-((1 + n2)/(-1 + n5)) + (2*(1 + n2)*pp)/(3*(-1 + n5)*tmm)))/(tarC3^n3*tarC4^n4) + (tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5)*(n3/(-1 + n5) - (n3*pp)/(3*(-1 + n5)*tmm)))/tarC1^n1 + (tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC5^(2 - n5)*(-(n3/(-1 + n5)) + (n3*pp)/(3*(-1 + n5)*tmm)))/(tarC1^n1*tarC4^n4) + (tarC2^(-1 - n2)*tarC4^(-1 - n4)*tarC5^(2 - n5)*(n4/(-1 + n5) - (2*n4*pp)/(3*(-1 + n5)*tmm)))/(tarC1^n1*tarC3^n3) + (tarC2^(-1 - n2)*tarC3^(1 - n3)*tarC4^(-1 - n4)*tarC5^(1 - n5)*(-(n4/(-1 + n5)) + (2*n4*pp)/(3*(-1 + n5)*tmm)))/tarC1^n1 + (tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC5^(1 - n5)*((3*n1*pp)/(2*(-1 + n5)) - (n1*pp^2)/(3*(-1 + n5)*tmm)))/(tarC3^n3*tarC4^n4) + (tarC2^(-2 - n2)*tarC5^(1 - n5)*(((1 + n2)*pp)/(2*(-1 + n5)) - ((1 + n2)*pp^2)/(3*(-1 + n5)*tmm)))/(tarC1^n1*tarC3^n3*tarC4^n4) + (tarC2^(-1 - n2)*tarC5^(1 - n5)*((1 + n1 + n2 - 2*n3)/(2*(-1 + n5)) - (pp*rat(1 + d - 2*n1 + n2 - 3*n3, 1))/(3*(-1 + n5)*tmm)))/(tarC1^n1*tarC3^n3*tarC4^n4);



*   n2 != 0 && n4 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?neg_/tarC3^n3?pos_/tarC4^n4?{>1}/tarC5^n5?pos_ =
-((n3*tarC1^(1 - n1)*tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC4^(1 - n4))/((-1 + n4)*tarC5^n5)) + (tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC3^(1 - n3)*tarC4^(1 - n4)*(n1/(-1 + n4) - (2*n1*pp)/(3*(-1 + n4)*tmm)))/tarC5^n5 + (tarC2^(-2 - n2)*tarC4^(2 - n4)*((1 + n2)/(-1 + n4) - (2*(1 + n2)*pp)/(3*(-1 + n4)*tmm)))/(tarC1^n1*tarC3^n3*tarC5^n5) + (tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC4^(1 - n4)*((-2*n1*pp)/(-1 + n4) + (2*n1*pp^2)/(3*(-1 + n4)*tmm)))/(tarC3^n3*tarC5^n5) + (tarC2^(-2 - n2)*tarC4^(1 - n4)*((-2*(1 + n2)*pp)/(-1 + n4) + (2*(1 + n2)*pp^2)/(3*(-1 + n4)*tmm)))/(tarC1^n1*tarC3^n3*tarC5^n5) + (n3*pp*tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n4)*tarC1^n1*tmm) + ((1 + n2)*pp*tarC2^(-2 - n2)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n4)*tarC1^n1*tarC3^n3*tmm) + (n1*pp*tarC1^(-1 - n1)*tarC2^(-1 - n2)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n4)*tarC3^n3*tmm) + (pp*tarC2^(-1 - n2)*tarC5^(1 - n5))/(3*tarC1^n1*tarC3^n3*tarC4^n4*tmm) - ((1 + n2)*pp*tarC1^(1 - n1)*tarC2^(-2 - n2)*tarC4^(1 - n4))/(3*(-1 + n4)*tarC3^n3*tarC5^n5*tmm) - (n1*pp*tarC1^(-1 - n1)*tarC4^(1 - n4))/(3*(-1 + n4)*tarC2^n2*tarC3^n3*tarC5^n5*tmm) - (n3*pp*tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC4^(2 - n4))/(3*(-1 + n4)*tarC1^n1*tarC5^n5*tmm) - (pp*tarC2^(-1 - n2)*tarC3^(1 - n3))/(3*tarC1^n1*tarC4^n4*tarC5^n5*tmm) + (tarC2^(-1 - n2)*tarC4^(1 - n4)*((-2 - n1 - n2 + n3 + n4)/(-1 + n4) + (pp*rat(2 + 2*d - n1 - n2 - 3*n3 - 3*n4, 1))/(3*(-1 + n4)*tmm)))/(tarC1^n1*tarC3^n3*tarC5^n5);



*   n5 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?neg0_/tarC3^n3?pos_/tarC4^n4?pos_/tarC5^n5?{>1} =
(n4*tarC3^(1 - n3)*tarC4^(-1 - n4)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC1^n1*tarC2^n2*tmm) + (n3*tarC3^(-1 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC1^n1*tarC2^n2*tmm) + (n2*tarC2^(-1 - n2)*tarC4^(1 - n4)*tarC5^(1 - n5))/(6*(-1 + n5)*tarC1^n1*tarC3^n3*tmm) + (n1*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC5^(1 - n5))/(6*(-1 + n5)*tarC2^n2*tarC4^n4*tmm) + (n2*tarC1^(1 - n1)*tarC2^(-1 - n2)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC3^n3*tarC4^n4*tmm) - (n2*pp*tarC2^(-1 - n2)*tarC5^(1 - n5))/(6*(-1 + n5)*tarC1^n1*tarC3^n3*tarC4^n4*tmm) + (n1*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC5^(1 - n5))/(3*(-1 + n5)*tarC3^n3*tarC4^n4*tmm) - (n1*pp*tarC1^(-1 - n1)*tarC5^(1 - n5))/(6*(-1 + n5)*tarC2^n2*tarC3^n3*tarC4^n4*tmm) - (n4*tarC4^(-1 - n4)*tarC5^(2 - n5))/(3*(-1 + n5)*tarC1^n1*tarC2^n2*tarC3^n3*tmm) - (n3*tarC3^(-1 - n3)*tarC5^(2 - n5))/(3*(-1 + n5)*tarC1^n1*tarC2^n2*tarC4^n4*tmm) - (n2*tarC2^(-1 - n2)*tarC5^(2 - n5))/(3*(-1 + n5)*tarC1^n1*tarC3^n3*tarC4^n4*tmm) - (n1*tarC1^(-1 - n1)*tarC5^(2 - n5))/(3*(-1 + n5)*tarC2^n2*tarC3^n3*tarC4^n4*tmm) + (tarC5^(1 - n5)*rat(6 + 2*d - n1 - n2 - 6*n5, 1))/(6*(-1 + n5)*tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*tmm);



*   n4 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?neg0_/tarC3^n3?pos_/tarC4^n4?{>1}/tarC5^n5?pos_ =
(tarC2^(-1 - n2)*tarC4^(1 - n4)*(-(n2/(-1 + n4)) + (n2*pp)/(3*(-1 + n4)*tmm)))/(tarC1^n1*tarC3^n3*tarC5^n5) + (2*n3*tarC3^(-1 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n4)*tarC1^n1*tarC2^n2*tmm) - (n2*tarC2^(-1 - n2)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n4)*tarC1^n1*tarC3^n3*tmm) + (2*n1*tarC1^(-1 - n1)*tarC4^(1 - n4)*tarC5^(1 - n5))/(3*(-1 + n4)*tarC2^n2*tarC3^n3*tmm) - tarC5^(1 - n5)/(3*tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*tmm) - (n1*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC4^(1 - n4))/(3*(-1 + n4)*tarC2^n2*tarC5^n5*tmm) + (n2*tarC1^(1 - n1)*tarC2^(-1 - n2)*tarC4^(1 - n4))/(3*(-1 + n4)*tarC3^n3*tarC5^n5*tmm) - (2*n1*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC4^(1 - n4))/(3*(-1 + n4)*tarC3^n3*tarC5^n5*tmm) + (n1*pp*tarC1^(-1 - n1)*tarC4^(1 - n4))/(3*(-1 + n4)*tarC2^n2*tarC3^n3*tarC5^n5*tmm) - (2*n3*tarC3^(-1 - n3)*tarC4^(2 - n4))/(3*(-1 + n4)*tarC1^n1*tarC2^n2*tarC5^n5*tmm) - (n2*tarC2^(-1 - n2)*tarC4^(2 - n4))/(3*(-1 + n4)*tarC1^n1*tarC3^n3*tarC5^n5*tmm) + tarC3^(1 - n3)/(3*tarC1^n1*tarC2^n2*tarC4^n4*tarC5^n5*tmm) + (tarC4^(1 - n4)*rat(3 + d + n1 - 2*n2 - 3*n4, 1))/(3*(-1 + n4)*tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5*tmm);



*   n3 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?neg0_/tarC3^n3?{>1}/tarC4^n4?pos_/tarC5^n5?pos_ =
(tarC1^(-1 - n1)*tarC3^(1 - n3)*(-(n1/(-1 + n3)) + (n1*pp)/(3*(-1 + n3)*tmm)))/(tarC2^n2*tarC4^n4*tarC5^n5) + (2*n4*tarC3^(1 - n3)*tarC4^(-1 - n4)*tarC5^(1 - n5))/(3*(-1 + n3)*tarC1^n1*tarC2^n2*tmm) + (2*n2*tarC2^(-1 - n2)*tarC3^(1 - n3)*tarC5^(1 - n5))/(3*(-1 + n3)*tarC1^n1*tarC4^n4*tmm) - (n1*tarC1^(-1 - n1)*tarC3^(1 - n3)*tarC5^(1 - n5))/(3*(-1 + n3)*tarC2^n2*tarC4^n4*tmm) - tarC5^(1 - n5)/(3*tarC1^n1*tarC2^n2*tarC3^n3*tarC4^n4*tmm) - (2*n4*tarC3^(2 - n3)*tarC4^(-1 - n4))/(3*(-1 + n3)*tarC1^n1*tarC2^n2*tarC5^n5*tmm) - (n2*tarC2^(-1 - n2)*tarC3^(1 - n3)*tarC4^(1 - n4))/(3*(-1 + n3)*tarC1^n1*tarC5^n5*tmm) + tarC4^(1 - n4)/(3*tarC1^n1*tarC2^n2*tarC3^n3*tarC5^n5*tmm) - (2*n2*tarC1^(1 - n1)*tarC2^(-1 - n2)*tarC3^(1 - n3))/(3*(-1 + n3)*tarC4^n4*tarC5^n5*tmm) + (n2*pp*tarC2^(-1 - n2)*tarC3^(1 - n3))/(3*(-1 + n3)*tarC1^n1*tarC4^n4*tarC5^n5*tmm) + (n1*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(1 - n3))/(3*(-1 + n3)*tarC4^n4*tarC5^n5*tmm) - (n1*tarC1^(-1 - n1)*tarC3^(2 - n3))/(3*(-1 + n3)*tarC2^n2*tarC4^n4*tarC5^n5*tmm) + (tarC3^(1 - n3)*rat(3 + d - 2*n1 + n2 - 3*n3, 1))/(3*(-1 + n3)*tarC1^n1*tarC2^n2*tarC4^n4*tarC5^n5*tmm);



*   n5 != 1
id,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?neg0_/tarC3^n3?pos_/tarC4^n4?pos_/tarC5^n5?{>1} =
(n4*tarC4^(-1 - n4)*tarC5^(1 - n5))/((-1 + n5)*tarC1^n1*tarC2^n2*tarC3^n3) - ((1 + n3)*tarC3^(-2 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5))/((-1 + n5)*tarC1^n1*tarC2^n2) + (n2*tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC4^(1 - n4)*tarC5^(1 - n5))/(2*(-1 + n5)*tarC1^n1) + (n2*tarC1^(1 - n1)*tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC5^(1 - n5))/((-1 + n5)*tarC4^n4) - (n2*pp*tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC5^(1 - n5))/(2*(-1 + n5)*tarC1^n1*tarC4^n4) - (n1*tarC1^(-1 - n1)*tarC2^(1 - n2)*tarC3^(-1 - n3)*tarC5^(1 - n5))/((-1 + n5)*tarC4^n4) + (n1*pp*tarC1^(-1 - n1)*tarC3^(-1 - n3)*tarC5^(1 - n5))/(2*(-1 + n5)*tarC2^n2*tarC4^n4) + ((n1 - n2)*tarC3^(-1 - n3)*tarC5^(1 - n5))/(2*(-1 + n5)*tarC1^n1*tarC2^n2*tarC4^n4) - (n1*tarC1^(-1 - n1)*tarC5^(1 - n5))/(2*(-1 + n5)*tarC2^n2*tarC3^n3*tarC4^n4) - (n4*tarC3^(-1 - n3)*tarC4^(-1 - n4)*tarC5^(2 - n5))/((-1 + n5)*tarC1^n1*tarC2^n2) + ((1 + n3)*tarC3^(-2 - n3)*tarC5^(2 - n5))/((-1 + n5)*tarC1^n1*tarC2^n2*tarC4^n4) - (n2*tarC2^(-1 - n2)*tarC3^(-1 - n3)*tarC5^(2 - n5))/((-1 + n5)*tarC1^n1*tarC4^n4) + (n1*tarC1^(-1 - n1)*tarC3^(-1 - n3)*tarC5^(2 - n5))/((-1 + n5)*tarC2^n2*tarC4^n4) + (tarC3^(-1 - n3)*tarC4^(1 - n4))/(tarC1^n1*tarC2^n2*tarC5^n5);



*   n3 == 1 && n4 == 1 && n5 == 1 && n2 != 0
id,only,ifmatch->sortme 1/tarC1^n1?neg0_/tarC2^n2?neg_/tarC3/tarC4/tarC5 =
-((tarC1^(1 - n1)*tarC2^(-1 - n2))/((-1 + n1)*tarC3*tarC4^2)) + (tarC1^(1 - n1)*tarC2^(-1 - n2))/((-1 + n1)*tarC3^2*tarC4) - ((1 + n2)*tarC1^(1 - n1)*tarC2^(-2 - n2))/((-1 + n1)*tarC3*tarC4) + tarC2^(-1 - n2)/(tarC1^n1*tarC3*tarC4) + (tarC1^(1 - n1)*tarC2^(-1 - n2))/((-1 + n1)*tarC3*tarC5^2) - (tarC1^(1 - n1)*tarC2^(-1 - n2))/((-1 + n1)*tarC4*tarC5^2) - (tarC1^(1 - n1)*tarC2^(-1 - n2))/((-1 + n1)*tarC3^2*tarC5) + ((1 + n2)*tarC1^(1 - n1)*tarC2^(-2 - n2))/(2*(-1 + n1)*tarC3*tarC5) + (tarC1^(1 - n1)*tarC2^(-1 - n2))/((-1 + n1)*tarC4^2*tarC5) - tarC2^(-1 - n2)/(2*tarC1^n1*tarC4*tarC5) - ((1 + n2)*pp*tarC1^(1 - n1)*tarC2^(-2 - n2))/(2*(-1 + n1)*tarC3*tarC4*tarC5) + ((1 + n2)*tarC1^(2 - n1)*tarC2^(-2 - n2))/((-1 + n1)*tarC3*tarC4*tarC5) + ((-2 + n1 - n2)*tarC1^(1 - n1)*tarC2^(-1 - n2))/(2*(-1 + n1)*tarC3*tarC4*tarC5) + (pp*tarC2^(-1 - n2))/(2*tarC1^n1*tarC3*tarC4*tarC5);



*   n2 == 0 && n3 == 1 && n4 == 1 && n5 == 1 && n1 != 0
if(count(tarC2,-1) == 0) id,only,ifmatch->sortme 1/tarC1^n1?neg_/tarC3/tarC4/tarC5 =
(2*pp*tarC1^(-1 - n1)*rat(1, -2 + d - n1))/(3*tarC3^2*tarC4) + (2*(1 + n1)*pp*tarC1^(-2 - n1)*rat(1, -2 + d - n1))/(3*tarC3*tarC4) - (2*pp*tarC1^(-1 - n1)*rat(1, -2 + d - n1))/(3*tarC3^2*tarC5) + ((-1 + n1)*tarC1^(-1 - n1)*rat(1, -2 + d - n1))/(tarC4*tarC5) - (2*(1 + n1)*pp*tarC1^(-2 - n1)*tarC2*rat(1, -2 + d - n1))/(3*tarC3*tarC4*tarC5) - (pp*tarC1^(-1 - n1)*tmm*rat(1, -2 + d - n1))/(tarC3*tarC4^2*tarC5) + (tarC1^(-1 - n1)*((2*pp*rat(1, -2 + d - n1))/3 - 2*tmm*rat(1, -2 + d - n1)))/(tarC3*tarC4^2) + (tarC1^(-1 - n1)*((-2*pp*rat(1, -2 + d - n1))/3 + tmm*rat(1, -2 + d - n1)))/(tarC4^2*tarC5) + (tarC1^(-2 - n1)*((-4*(1 + n1)*pp*rat(1, -2 + d - n1))/3 + 2*(1 + n1)*tmm*rat(1, -2 + d - n1)))/(tarC4*tarC5) + (tarC1^(-2 - n1)*((4*(1 + n1)*pp^2*rat(1, -2 + d - n1))/3 - 4*(1 + n1)*pp*tmm*rat(1, -2 + d - n1)))/(tarC3*tarC4*tarC5) + ((-2 + n1)*rat(1, 2 - d + n1))/(n1*tarC1^n1*tarC3*tarC4^2) - ((-2 + n1)*rat(1, 2 - d + n1))/(n1*tarC1^n1*tarC4^2*tarC5) - (4*tmm*rat(1, 2 - d + n1))/(n1*tarC1^n1*tarC3*tarC4^3) + (2*tmm*rat(1, 2 - d + n1))/(n1*tarC1^n1*tarC3^2*tarC4^2) - (2*tmm*rat(1, 2 - d + n1))/(n1*tarC1^n1*tarC4^2*tarC5^2) + (4*tmm*rat(1, 2 - d + n1))/(n1*tarC1^n1*tarC4^3*tarC5) + (pp*tarC1^(-1 - n1)*rat(-11 + 4*d - 5*n1, -2 + d - n1))/(3*tarC3*tarC4*tarC5);

endif;

goto endrec;         
la sortme;
$irep = 0;
la endrec;
        
ModuleOption,minimum,$irep;
.sort:red2l-T2-`$repcount++';
#redefine irep "`$irep'"
#enddo
#endprocedure



* 
* count bad for neg num!!!!!
* 

#procedure to1l(TOPO)
        if(count(int`TOPO',1));        
***************************************************         
*       GxG        
        if(((count(tarC1,-1) > 0) || (count(tarC3,-1) > 0)) && ((count(tarC2,-1) > 0) || (count(tarC4,-1) > 0)) && (count(tarC5,1) >= 0));
        Multiply intGxG/int`TOPO';        

***************************************************         
*       TxG
        elseif((count(tarC1,1) >= 0) && ((count(tarC2,-1) > 0) || (count(tarC4,-1) > 0)) && (count(tarC3,1) >= 0) && (count(tarC5,-1) > 0));
        Multiply intTxG/int`TOPO';        

***************************************************         
*       GxT
        elseif(((count(tarC1,-1) > 0) || (count(tarC3,-1) > 0)) && (count(tarC2,1) >= 0) && (count(tarC4,1) >= 0) && (count(tarC5,-1) > 0));
        Multiply intGxT/int`TOPO';        
        endif;        
        endif;
#endprocedure



* Reduce two-loop tadpole by hand
#procedure red2l

        if(count(intFG,1)) Multiply sdim(0)*intF/intFG;


        #message
        #message reduce F                
        #message
        #$repcount = 1;        
        #do irep=1,1
                #$irep = 1;                
                if(count(intF,1));
                #call redF
                id,only intF/tarC1/tarC2/tarC3/tarC4/tarC5 = intFG*TFI(1,1,1,1,1);                
                endif;
*                 #breakdo
                goto endrec;
                la dopartfrac;
*                 #call partfrac
                $irep = 0;
                la endrec;
                
                ModuleOption,minimum,$irep;
                .sort:redF-noshift-`$repcount++';
                #if `$irep'==0
                        #call partfrac        
                #endif                                
                #redefine irep "`$irep'"
        #enddo
        id sdim(0) = 1;        
        #printtimes
        .sort:red-F-done;

        #call zeroTFI(F)

        #call uniqueVnew(F)
        #call uniqueJnew(F)
        #call uniqueT2new(F)
        #call to1l(F)
*       reduce V
        #message
        #message reduce V                
        #message
        #call red01111
        #printtimes
        .sort:red-V-done;        
        #call zeroTFI(V)        
        #call partfrac
*       rename V
        #call uniqueVnew(V)
        id,only intV/tarC2/tarC3/tarC4/tarC5 = intFG*TFI(0,1,1,1,1); 
        #call uniqueJnew(V)
        #call uniqueT2new(V)
        #call to1l(V)
*       reduce J
        #message
        #message reduce J
        #message        
        #call red01101
        #printtimes
        .sort:red-J-done;        
        #call zeroTFI(J)
        #call partfrac
*       rename J
        #call uniqueJnew(J)
        id,only intJ/tarC2^2/tarC3/tarC5 = intFG*TFI(2,0,0,1,1);
        id,only intJ/tarC2/tarC3/tarC5   = intFG*TFI(1,0,0,1,1);                
        #call uniqueT2new(J)
        #call to1l(J)
*       reduce T2
        id tarC1^n1?pos_ = (tk1.tk1 - tmm)^n1;
        id tarC2^n2?pos_ = (tk2.tk2 - tmm)^n2;
        id tarC3^n3?pos_ = (tk1.tk1 - 2*tk1.tp + pp - tmm)^n3;
        id tarC4^n4?pos_ = (tk2.tk2 - 2*tk2.tp + pp - tmm)^n4;
        id tarC5^n5?pos_ = (tk1.tk1 - 2*tk1.tk2 + tk2.tk2 - tmm)^n5;
        
        
*       all tarC in num con to SP !!!        
        #call spcontract1(T2,1)
        
        
        if((count(tarC1,-1) == 0) && (count(tarC3,-1) == 0)) Multiply replace_(intT2,intTxG);
        if((count(tarC2,-1) == 0) && (count(tarC4,-1) == 0)) Multiply replace_(intT2,intGxT);
        
        
        #call redTad125(T2)
        
        
* Now reduce remaining low powers of (k1.p) and (k2.p)
* eq B.10 in Davydychev, Smirnov, Tausk 
* Nucl.Phys. B410 (1993) 325-342
        
        if(count(intT2,1));
        id tp.tk1^x?pos0_*tp.tk2^y?pos0_/tarC1^n1?pos_/tarC2^n2?pos_/tarC5^n5? = 
        mod_(x+y+1,2)*(1/2)^(x+y)*fac_(x)*fac_(y)*pp^((x+y)/2)*PochhammerINV((x+y)/2,d/2)*
        sum_(j3,0,min_(x,y),
        mod_(x-j3+1,2)*mod_(y-j3+1,2)*tk1.tk1^((x-j3)/2)*tk2.tk2^((y-j3)/2)*(2*tk1.tk2)^j3*invfac_((x-j3)/2)*invfac_((y-j3)/2)*invfac_(j3))/tarC1^n1/tarC2^n2/tarC5^n5;
        endif;
        
        
        #call spcontract1(T2,1)
        
        
        
*       T2=125        
        if((count(tarC3,-1) == 0) && (count(tarC4,-1) == 0) && (count(tarC1,-1) > 0) && (count(tarC2,-1) > 0) && (count(tarC5,-1) == 0)) Multiply replace_(intT2,intTxT);
        if((count(tarC3,-1) == 0) && (count(tarC4,-1) == 0) && (count(tarC1,-1) > 0) && (count(tarC2,-1) == 0) && (count(tarC5,-1) > 0)) Multiply replace_(intT2,intTxT);
        if((count(tarC3,-1) == 0) && (count(tarC4,-1) == 0) && (count(tarC1,-1) == 0) && (count(tarC2,-1) > 0) && (count(tarC5,-1) > 0)) Multiply replace_(intT2,intTxT);

*       T2=345        
        if((count(tarC1,-1) == 0) && (count(tarC2,-1) == 0) && (count(tarC3,-1) > 0) && (count(tarC4,-1) > 0) && (count(tarC5,-1) == 0)) Multiply replace_(intT2,intTxT);
        if((count(tarC1,-1) == 0) && (count(tarC2,-1) == 0) && (count(tarC3,-1) > 0) && (count(tarC4,-1) == 0) && (count(tarC5,-1) > 0)) Multiply replace_(intT2,intTxT);
        if((count(tarC1,-1) == 0) && (count(tarC2,-1) == 0) && (count(tarC3,-1) == 0) && (count(tarC4,-1) > 0) && (count(tarC5,-1) > 0)) Multiply replace_(intT2,intTxT);
        
        if(count(tarC1,-1) == 0) Multiply replace_(intT2,intTxG);
        if(count(tarC2,-1) == 0) Multiply replace_(intT2,intGxT);
        
        if((count(tk1,1) == 0) && (count(tk2,1) == 0)) id intT2/tarC1^n1?pos_/tarC2^n2?pos_/tarC5^n5?pos_ = intFG*T2(0,n1,n2,n5);
        if((count(tk1,1) == 0) && (count(tk2,1) == 0)) id intT2/tarC1^n1?pos_/tarC2^n2?pos_ = intFG*T1(0,n1)*T1(0,n2);
        
        
*       Reduce T2        
        #message
        #message reduce T2                
        #message
        #do i=1,1
                if(count(intFG,1));
                #call redT2
                if(match(T2(0,n1?pos_,n2?pos_,n3?pos_))) redefine i "0";
                endif;
                .sort:red-T2;
        #enddo
        #printtimes

        if(count(intFG,1));
        id T2(dp?,n1?pos_,n2?pos_,0) = T1(dp,n1)*T1(dp,n2);
        id T2(dp?,n1?pos_,0,n3?pos_) = T1(dp,n1)*T1(dp,n3);
        id T2(dp?,0,n2?pos_,n3?pos_) = T1(dp,n2)*T1(dp,n3);
        
        if(match(T2(dp?,n1?,n2?,n3?))) exit "Unreduced T2";
        endif;

        #call partfrac
        
        .sort
        #call zeroTFI(GxG)
        #call zeroTFI(GxT)
        #call zeroTFI(TxG)
        
        
*       Tensor reduction for factorized integrals         
        if(count(intGxG,1));
        #call tensG(tarC1,tarC3,tk1,tp)
        endif;
        .sort:intGxG-1;
        if(count(intGxG,1));
        #call tensG(tarC2,tarC4,tk2,tp)
        id tp.tp^n? = pp^n;        
        Multiply intFG/intGxG;        
        endif;        
        .sort:intGxG-2;
        
        if(count(intTxG,1));
*       First integrate tadpole with tk1,tk2 in denominator:        
        Multiply replace_(tk1,[tk1-tk2] + tk2);
* Mark terms wo den as zero        
        #call tensT1(tarC5,[tk1-tk2],0)
*         Multiply replace_(zero,0);
        endif;
        .sort:intTxG-1;
        if(count(intTxG,1));
*       tk1.tk2 in numerator wo denominator gives zero scale tadpole        
        Multiply replace_([tk1-tk2],tk1 - tk2);        
        #call tensG(tarC2,tarC4,tk2,tp)
        id tp.tp^n? = pp^n;        
        Multiply intFG/intTxG;        
        endif;        
        .sort:intTxG-2;

        if(count(intGxT,1));
        Multiply replace_(tk2, tk1 - [tk1-tk2]);
        #call tensT1(tarC5,[tk1-tk2],0)
        endif;
        .sort:intGxT-1;
        if(count(intGxT,1));
*       tk1.tk2 in numerator wo denominator gives zero scale tadpole        
        Multiply replace_([tk1-tk2],tk1 - tk2);        
        #call tensG(tarC1,tarC3,tk1,tp)
        id tp.tp^n? = pp^n;                
         Multiply intFG/intGxT;        
        endif;        
        .sort:intGxT-2;


*       One-loop tadpoles        
        if(count(intTxT,1) && (count(tarC1,-1) > 0) && (count(tarC2,-1) > 0));
        #call tensT1(tarC1,tk1,0)
        #call tensT1(tarC2,tk2,0)
*         Multiply replace_(zero,0);
        Multiply intFG/intTxT;        
        endif;
        .sort:intTxT-(12);

        if(count(intTxT,1) && (count(tarC1,-1) > 0) && (count(tarC5,-1) > 0));
        Multiply replace_(tk2, tk1 - [tk1-tk2]);
        #call tensT1(tarC1,tk1,0)
        #call tensT1(tarC5,[tk1-tk2],0)
*         Multiply replace_(zero,0);
        Multiply intFG/intTxT;        
        endif;
        .sort:intTxT-(15);

        if(count(intTxT,1) && (count(tarC2,-1) > 0) && (count(tarC5,-1) > 0));
        Multiply replace_(tk1, [tk1-tk2] + tk2);
        #call tensT1(tarC2,tk2,0)
        #call tensT1(tarC5,[tk1-tk2],0)
*         Multiply replace_(zero,0);
        Multiply intFG/intTxT;        
        endif;
        .sort:intTxT-(25);

        #call subpoch
        
*       Reducing one-loop G and T1        
        #message
        #message reduce G                
        #message
        
        #call drrG(FG)

        #call redG(FG)
        id G(0,1,1) = Gx11;
        
        if(count(intFG,1));
        id G(dp?,n1?,0) = T1(dp,n1);
        id G(dp?,0,n2?) = T1(dp,n2);
        endif;
        #printtimes

*       Now reduce one-loop tadpoles        
        #message
        #message drr T1                
        #message
        #do i=1,1
                if(count(intFG,1));
                #call drrT1
                if(match(T1(dp?pos_,n1?pos_))) redefine i "0";
                endif;
                .sort:drr-T1;
        #enddo
        

        #message
        #message reduce T1                
        #message
        #do i=1,1
                if(count(intFG,1));
                #call redT1
                if(match(T1(0,n1?pos_))) redefine i "0";
                endif;
                .sort:red-T1;
        #enddo

        if(count(intFG,1)) id T1(dp?,0) = 0;
        #call partfrac
        #call subpoch

#endprocedure



* 
* Planar 8 line
* 
*--#[ redFGnew :
#procedure redFGnew
* Reduce factorized topology
* F x G, where F and G are reduced wih TARCER separatelly


* [1, 0, 1, 1, 1, 1, 1, 1, 1, 0]

* F = 1,9,4,7,3
*
*
*         1/|\9          _6_
*        /  |  \        /   \
*   ---<    3    >--5--(     )-- 
*        \  |  /        \_ _/
*         4\|/7           8
*
* 
* G = 6,8

*         
* TARCER:
* 
* C1 = k1^2-tm^2
* C2 = k2^2-tm^2
* C3 = (k1-p)^2-tm^2
* C4 = (k2-p)^2-tm^2
* C5 = (k1-k2)^2-tm^2
*         


* Numerator reduction:
*
*  - all sp[k1,k4] and sp[k2,k4]
*    we treat using expression for
*    one-loop tensor integral k4(mu)...*k4(nu)...
*    leading to (p.k1), (p.k2), k1^2, k2^2

        if(count(intFG,1));

        if((count(d2,1) >= 0) && (count(d10,1) >= 0));

*       Numerator to scalar products        
        Multiply replace_(d2, k2.k2 - tmm, d10, k1.k1 - 2*k1.k2 + k2.k2 - 2*k1.k3 + 2*k2.k3 + k3.k3 - tmm);
        Multiply replace_(d1,tarC1, d9,tarC2, d4,tarC3, d7,tarC4, d3,tarC5, d5,[pp-tmm]);        


*       Inverse powers of denominators to scalar products
        id tarC1^n1?pos_=(k1.k1 - tmm)^n1;
        id tarC2^n2?pos_=(k1.k1 - 2*k1.k3 + k3.k3 - tmm)^n2;
        id tarC3^n3?pos_=(k4.k4 - tmm)^n3;
        id tarC4^n4?pos_=(k3.k3 - 2*k3.k4 + k4.k4 - tmm)^n4;
        id tarC5^n5?pos_=(k3.k3 - tmm)^n5;

*       Numerator to TARCER momenta definitions       
        Multiply replace_(k1,tk1,   
                          k2,tk1 - tk4,   
                          k3,tk1 - tk2,   
                          k4,tk1 - tp);

        endif;
        endif;
        .sort:redFG-num;

*       Now reduce with dimensional shift one-loop subdiagram        
        if(count(intFG,1));
        Multiply replace_(tk4,-mtk4, tp,-mtp);        
        #call tensG(d8,d6,mtk4,mtp)
        Multiply replace_(mtk4,-tk4, mtp,-tp);        
        id tp.tp = pp;        
        endif;        
        .sort:intFG-tensG;        

* All to tarC1...tarC5        

        id tk1.tk1^n?pos_ = (tarC1 + tmm)^n;
        id tk2.tk2^n?pos_ = (tarC2 + tmm)^n;
        id tk1.tp^n?pos_  = ((tarC1 - tarC3 + pp)/2)^n;
        id tk2.tp^n?pos_  = ((tarC2 - tarC4 + pp)/2)^n;
        id tk1.tk2^n?pos_ = ((tarC1 + tarC2 - tarC5 + tmm)/2)^n;
        
        #call red2l
        
#endprocedure
*--#] redFGnew :



** 
* *
*  *  One-dimensional recurence relations
* *
**

*--#[ rec1dm0 :
#procedure rec1dm0
#$repcount = 1;

#do irep = 1,1
#$irep = 1;
                
********************************************************************************                 
   id,ifmatch->sortme m0denJ2010X11100(n?{>1}) =

       + theta_( - 4 + n) * (
          + m0denJ0021X10110( - 3 + n)*rat( - 10*n^2*d + 30*n^2 + 9*n*d^2 - 24
         *n*d - 3*n - 9*d^2 + 31*d - 18,144*n^3 - 72*n^2*d - 144*n^2 + 108*n*d
          - 36*n - 36*d + 36)
          + 1/( - 36 + 36*n)*m0denJ2010X11100( - 3 + n)*n*rat(1,1)
          + 1/( - 36 + 36*n)*m0denJ2010X11100( - 3 + n)*rat(-4,1)
          )

       + theta_( - 3 + n) * (
          + m0denJ0001X01110( - 2 + n)*rat( - 2*n*d^2 + 8*n*d - 8*n + 3*d^3 - 
         19*d^2 + 40*d - 28, - 144*n^2 + 72*n*d - 36*d + 36)
          + m0denJ0011X10110( - 2 + n)*rat( - 6*n*d^2 + 34*n*d - 48*n + 9*d^3
          - 72*d^2 + 191*d - 168, - 144*n^2 + 72*n*d - 36*d + 36)
          + m0denJ0021X00110( - 2 + n)*rat(8*n^2*d^2 - 32*n^2*d + 32*n^2 - 12*
         n*d^3 + 56*n*d^2 - 72*n*d + 16*n + 3*d^4 - 14*d^3 + 17*d^2 - 2*d,144*
         n^4 - 144*n^3*d - 144*n^3 + 36*n^2*d^2 + 180*n^2*d - 36*n^2 - 54*n*
         d^2 - 18*n*d + 36*n + 18*d^2 - 18*d)
          + m0denJ0021X10110( - 2 + n)*rat(46*n^2*d - 138*n^2 - 27*n*d^2 + 28*
         n*d + 131*n + 27*d^2 - 60*d - 35,72*n^3 - 36*n^2*d - 72*n^2 + 54*n*d
          - 18*n - 18*d + 18)
          + 1/( - 18 + 18*n)*m0denJ2010X11100( - 2 + n)*n*rat(-7,1)
          + 1/( - 18 + 18*n)*m0denJ2010X11100( - 2 + n)*rat(21,1)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat(2*n*d^2 - 8*n*d + 8*n - 3*d^3 + 19*
         d^2 - 40*d + 28, - 36*n^2 + 18*n*d - 9*d + 9)
          + m0denJ0001X00110( - 1 + n)*rat( - 4*n^2*d^3 + 24*n^2*d^2 - 48*n^2*
         d + 32*n^2 + 12*n*d^4 - 98*n*d^3 + 300*n*d^2 - 408*n*d + 208*n - 3*
         d^5 + 23*d^4 - 60*d^3 + 48*d^2 + 32*d - 48,288*n^4 - 288*n^3*d - 288*
         n^3 + 72*n^2*d^2 + 360*n^2*d - 72*n^2 - 108*n*d^2 - 36*n*d + 72*n + 
         36*d^2 - 36*d)
          + m0denJ0011X10110( - 1 + n)*rat(24*n*d^2 - 136*n*d + 192*n - 27*d^3
          + 228*d^2 - 641*d + 600, - 144*n^2 + 72*n*d - 36*d + 36)
          + m0denJ0011X00110( - 1 + n)*rat( - 12*n^2*d^3 + 92*n^2*d^2 - 232*
         n^2*d + 192*n^2 + 36*n*d^4 - 354*n*d^3 + 1294*n*d^2 - 2084*n*d + 1248
         *n - 9*d^5 + 84*d^4 - 271*d^3 + 296*d^2 + 108*d - 288,288*n^4 - 288*
         n^3*d - 288*n^3 + 72*n^2*d^2 + 360*n^2*d - 72*n^2 - 108*n*d^2 - 36*n*
         d + 72*n + 36*d^2 - 36*d)
          + m0denJ0021X00110( - 1 + n)*rat( - 20*n^2*d^2 + 80*n^2*d - 80*n^2
          + 90*n*d^2 - 400*n*d + 440*n - 15*d^3 + 45*d^2 + 30*d - 120,72*n^4
          - 72*n^3*d - 72*n^3 + 18*n^2*d^2 + 90*n^2*d - 18*n^2 - 27*n*d^2 - 9*
         n*d + 18*n + 9*d^2 - 9*d)
          + m0denJ0021X10110( - 1 + n)*rat( - 178*n^2*d + 534*n^2 + 27*n*d^2
          + 310*n*d - 1075*n - 27*d^2 - 181*d + 688,144*n^3 - 72*n^2*d - 144*
         n^2 + 108*n*d - 36*n - 36*d + 36)
          + 1/( - 36 + 36*n)*m0denJ2010X11100( - 1 + n)*n*rat(49,1)
          + 1/( - 36 + 36*n)*m0denJ2010X11100( - 1 + n)*rat(-98,1)
          )

       + delta_( - 3 + n) * (
          + m0denJ0011X10110(0)*rat( - 18*d^3 + 171*d^2 - 523*d + 520,900*d - 
         6300)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat(3*d^5 - 47*d^4 + 272*d^3 - 744*d^2 + 976*d
          - 496,108*d^2 - 972*d + 2160)
          + m0denJ0011X10110(0)*rat(126*d^3 - 1113*d^2 + 3227*d - 3080,540*d
          - 2700)
          );

   id,ifmatch->sortme m0denJ2010X01100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m0denJ0021X00110( - 2 + n)*rat(4*n*d - 12*n - 3*d^2 + 9*d,18*n^2
          - 9*n*d - 18*n + 9*d)
          + 1/( - 9 + 9*n)*m0denJ2010X01100( - 2 + n)*n*rat(-1,1)
          + 1/( - 9 + 9*n)*m0denJ2010X01100( - 2 + n)*rat(3,1)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X00110( - 1 + n)*rat(2*n*d^2 - 8*n*d + 8*n - 3*d^3 + 18*
         d^2 - 36*d + 24, - 36*n^2 + 18*n*d + 36*n - 18*d)
          + m0denJ0011X00110( - 1 + n)*rat(6*n*d^2 - 34*n*d + 48*n - 9*d^3 + 
         69*d^2 - 174*d + 144, - 36*n^2 + 18*n*d + 36*n - 18*d)
          + m0denJ0021X00110( - 1 + n)*rat( - 20*n*d + 60*n + 40*d - 120,18*
         n^2 - 9*n*d - 18*n + 9*d)
          + 1/( - 9 + 9*n)*m0denJ2010X01100( - 1 + n)*n*rat(10,1)
          + 1/( - 9 + 9*n)*m0denJ2010X01100( - 1 + n)*rat(-20,1)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat(3*d^3 - 22*d^2 + 52*d - 40,18*d - 72)
          );

   id,ifmatch->sortme m0denJ1011X11110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m0denJ0021X10110( - 3 + n)*rat( - 4*n^2 + 6*n*d - 6*n + 3*d - 2,
         864*n^3 - 1296*n^2*d + 2592*n^2 + 648*n*d^2 - 2592*n*d + 2376*n - 108
         *d^3 + 648*d^2 - 1188*d + 648)
          )

       + theta_( - 3 + n) * (
          + m0denJ0001X01110( - 2 + n)*rat( - 2*n*d^2 + 8*n*d - 8*n - d^2 + 4*
         d - 4,864*n^3 - 1296*n^2*d + 2592*n^2 + 648*n*d^2 - 2592*n*d + 2376*n
          - 108*d^3 + 648*d^2 - 1188*d + 648)
          + m0denJ0011X10110( - 2 + n)*rat( - 6*n*d^2 + 34*n*d - 48*n - 3*d^2
          + 17*d - 24,864*n^3 - 1296*n^2*d + 2592*n^2 + 648*n*d^2 - 2592*n*d
          + 2376*n - 108*d^3 + 648*d^2 - 1188*d + 648)
          + m0denJ0021X00110( - 2 + n)*rat(4*n^2*d - 8*n^2 - 4*n*d^2 + 6*n*d
          + 4*n - 2*d^3 + 5*d^2 - 2*d,864*n^4 - 1728*n^3*d + 2592*n^3 + 1296*
         n^2*d^2 - 3888*n^2*d + 2376*n^2 - 432*n*d^3 + 1944*n*d^2 - 2376*n*d
          + 648*n + 54*d^4 - 324*d^3 + 594*d^2 - 324*d)
          + m0denJ0021X10110( - 2 + n)*rat(28*n^2 - 18*n*d - 2*n - 9*d - 8,432
         *n^3 - 648*n^2*d + 1296*n^2 + 324*n*d^2 - 1296*n*d + 1188*n - 54*d^3
          + 324*d^2 - 594*d + 324)
          + m0denJ1011X01110( - 2 + n)*rat( - n*d + 2*n + d^2 - 3*d + 2,48*n^2
          - 48*n*d + 120*n + 12*d^2 - 60*d + 72)
          + m0denJ1011X11110( - 2 + n)*rat( - 2*n + 3*d - 6,24*n - 12*d + 36)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat( - n*d^2 + 4*n*d - 4*n + 18*d^3 - 
         86*d^2 + 128*d - 56,864*n^3 - 1296*n^2*d + 2592*n^2 + 648*n*d^2 - 
         2592*n*d + 2376*n - 108*d^3 + 648*d^2 - 1188*d + 648)
          + m0denJ0001X11110( - 1 + n)*rat(2*n*d - 4*n + d - 2,48*n^2 - 48*n*d
          + 120*n + 12*d^2 - 60*d + 72)
          + m0denJ0001X00110( - 1 + n)*rat(4*n*d^3 - 24*n*d^2 + 48*n*d - 32*n
          + d^4 - 6*d^3 + 12*d^2 - 8*d,864*n^4 - 1728*n^3*d + 2592*n^3 + 1296*
         n^2*d^2 - 3888*n^2*d + 2376*n^2 - 432*n*d^3 + 1944*n*d^2 - 2376*n*d
          + 648*n + 54*d^4 - 324*d^3 + 594*d^2 - 324*d)
          + m0denJ0011X01110( - 1 + n)*rat(n*d^2 - 5*n*d + 6*n - d^3 + 6*d^2
          - 11*d + 6,24*n^3 - 36*n^2*d + 72*n^2 + 18*n*d^2 - 72*n*d + 66*n - 3
         *d^3 + 18*d^2 - 33*d + 18)
          + m0denJ0011X10110( - 1 + n)*rat(54*n^2*d - 144*n^2 - 84*n*d^2 + 341
         *n*d - 312*n + 36*d^3 - 201*d^2 + 325*d - 120,864*n^3 - 1296*n^2*d + 
         2592*n^2 + 648*n*d^2 - 2592*n*d + 2376*n - 108*d^3 + 648*d^2 - 1188*d
          + 648)
          + m0denJ0011X11100( - 1 + n)*rat( - 2*n^2*d + 4*n^2 + 4*n*d^2 - 15*n
         *d + 14*n - 2*d^3 + 11*d^2 - 19*d + 10,96*n^3 - 144*n^2*d + 288*n^2
          + 72*n*d^2 - 288*n*d + 264*n - 12*d^3 + 72*d^2 - 132*d + 72)
          + m0denJ0011X11110( - 1 + n)*rat(2*n*d - 6*n - 3*d^2 + 15*d - 18,24*
         n^2 - 24*n*d + 60*n + 6*d^2 - 30*d + 36)
          + m0denJ0011X00110( - 1 + n)*rat(12*n*d^3 - 92*n*d^2 + 232*n*d - 192
         *n + 3*d^4 - 23*d^3 + 58*d^2 - 48*d,864*n^4 - 1728*n^3*d + 2592*n^3
          + 1296*n^2*d^2 - 3888*n^2*d + 2376*n^2 - 432*n*d^3 + 1944*n*d^2 - 
         2376*n*d + 648*n + 54*d^4 - 324*d^3 + 594*d^2 - 324*d)
          + m0denJ0021X00110( - 1 + n)*rat( - 76*n^2*d + 152*n^2 - 4*n*d^2 + 
         170*n*d - 324*n - 9*d^3 + 57*d^2 - 78*d,864*n^4 - 1728*n^3*d + 2592*
         n^3 + 1296*n^2*d^2 - 3888*n^2*d + 2376*n^2 - 432*n*d^3 + 1944*n*d^2
          - 2376*n*d + 648*n + 54*d^4 - 324*d^3 + 594*d^2 - 324*d)
          + m0denJ0021X10110( - 1 + n)*rat( - 170*n^2 + 144*n*d - 199*n - 45*
         d^2 + 162*d - 57,432*n^3 - 648*n^2*d + 1296*n^2 + 324*n*d^2 - 1296*n*
         d + 1188*n - 54*d^3 + 324*d^2 - 594*d + 324)
          + m0denJ1001X11110( - 1 + n)*rat(d - 3, - 12*n + 6*d - 18)
          + m0denJ1011X01110( - 1 + n)*rat(6*n*d - 12*n - 7*d^2 + 22*d - 16,96
         *n^2 - 96*n*d + 240*n + 24*d^2 - 120*d + 144)
          + m0denJ1011X11110( - 1 + n)*rat(7*n - 7*d + 16,12*n - 6*d + 18)
          )

       + delta_( - 3 + n) * (
          + m0denJ0011X10110(0)*rat( - 42*d^2 + 217*d - 280,540*d^3 - 12960*
         d^2 + 103140*d - 272160)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat( - d^4 - 2*d^3 + 36*d^2 - 88*d + 64,54*d^4
          - 1188*d^3 + 9666*d^2 - 34452*d + 45360)
          + m0denJ0011X10110(0)*rat(42*d^2 - 217*d + 280,108*d^3 - 1944*d^2 + 
         11556*d - 22680)
          + m0denJ1011X01110(0)*rat(d^2 - 5*d + 6,12*d^2 - 156*d + 504)
          + m0denJ1011X11110(0)*rat( - 3*d + 10,12*d - 84)
          );

   id,ifmatch->sortme m0denJ1011X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m0denJ0021X00110( - 2 + n)*rat(4*n^2 - 4*n*d,216*n^3 - 324*n^2*d
          + 324*n^2 + 162*n*d^2 - 324*n*d + 108*n - 27*d^3 + 81*d^2 - 54*d)
          + m0denJ1011X01110( - 2 + n)*rat( - n + d - 1,12*n - 6*d + 12)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat(n*d - 2*n,24*n^2 - 24*n*d + 36*n + 
         6*d^2 - 18*d + 12)
          + m0denJ0001X11110( - 1 + n)*rat(d - 3, - 12*n + 6*d - 12)
          + m0denJ0001X00110( - 1 + n)*rat(2*n*d^2 - 8*n*d + 8*n,216*n^3 - 324
         *n^2*d + 324*n^2 + 162*n*d^2 - 324*n*d + 108*n - 27*d^3 + 81*d^2 - 54
         *d)
          + m0denJ0011X01110( - 1 + n)*rat(n*d - 3*n - d^2 + 4*d - 3,12*n^2 - 
         12*n*d + 18*n + 3*d^2 - 9*d + 6)
          + m0denJ0011X00110( - 1 + n)*rat(6*n*d^2 - 34*n*d + 48*n,216*n^3 - 
         324*n^2*d + 324*n^2 + 162*n*d^2 - 324*n*d + 108*n - 27*d^3 + 81*d^2
          - 54*d)
          + m0denJ0021X00110( - 1 + n)*rat( - 76*n^2 + 36*n*d + 62*n - 9*d^2
          + 9*d,216*n^3 - 324*n^2*d + 324*n^2 + 162*n*d^2 - 324*n*d + 108*n - 
         27*d^3 + 81*d^2 - 54*d)
          + m0denJ1011X01110( - 1 + n)*rat(14*n - 11*d + 16,24*n - 12*d + 24)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat(4*d^2 - 16*d + 16,27*d^3 - 405*d^2 + 1998*
         d - 3240)
          + m0denJ1011X01110(0)*rat( - d + 3,6*d - 36)
          );

   id,ifmatch->sortme m0denJ1010X11100(n?{>1}) =

       + theta_( - 4 + n) * (
          + m0denJ0021X10110( - 3 + n)*rat( - 2*n + 3*d - 2,24*n^2 - 12*n*d + 
         6*d - 6)
          )

       + theta_( - 3 + n) * (
          + m0denJ0001X01110( - 2 + n)*rat(d^2 - 4*d + 4, - 24*n^2 + 12*n*d - 
         6*d + 6)
          + m0denJ0011X10110( - 2 + n)*rat(3*d^2 - 17*d + 24, - 24*n^2 + 12*n*
         d - 6*d + 6)
          + m0denJ0021X00110( - 2 + n)*rat(2*n^2*d - 4*n^2 - 4*n*d^2 + 9*n*d
          - 2*n + d^3 - 2*d^2,24*n^4 - 24*n^3*d - 24*n^3 + 6*n^2*d^2 + 30*n^2*
         d - 6*n^2 - 9*n*d^2 - 3*n*d + 6*n + 3*d^2 - 3*d)
          + m0denJ0021X10110( - 2 + n)*rat(14*n - 9*d - 8,12*n^2 - 6*n*d + 3*d
          - 3)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat( - 2*d^2 + 8*d - 8, - 12*n^2 + 6*n*
         d - 3*d + 3)
          + m0denJ0001X00110( - 1 + n)*rat(4*n*d^3 - 24*n*d^2 + 48*n*d - 32*n
          - d^4 + 5*d^3 - 6*d^2 - 4*d + 8,48*n^4 - 48*n^3*d - 48*n^3 + 12*n^2*
         d^2 + 60*n^2*d - 12*n^2 - 18*n*d^2 - 6*n*d + 12*n + 6*d^2 - 6*d)
          + m0denJ0011X10110( - 1 + n)*rat( - 6*n^2*d + 18*n^2 + 18*n*d^2 - 
         109*n*d + 165*n - 18*d^2 + 112*d - 174,48*n^3 - 24*n^2*d - 48*n^2 + 
         36*n*d - 12*n - 12*d + 12)
          + m0denJ0011X00110( - 1 + n)*rat(12*n*d^3 - 92*n*d^2 + 232*n*d - 192
         *n - 3*d^4 + 20*d^3 - 35*d^2 - 10*d + 48,48*n^4 - 48*n^3*d - 48*n^3
          + 12*n^2*d^2 + 60*n^2*d - 12*n^2 - 18*n*d^2 - 6*n*d + 12*n + 6*d^2
          - 6*d)
          + m0denJ0021X00110( - 1 + n)*rat( - 20*n^2*d + 40*n^2 + 70*n*d - 140
         *n - 10*d^2 + 40,24*n^4 - 24*n^3*d - 24*n^3 + 6*n^2*d^2 + 30*n^2*d - 
         6*n^2 - 9*n*d^2 - 3*n*d + 6*n + 3*d^2 - 3*d)
          + m0denJ0021X10110( - 1 + n)*rat( - 80*n + 9*d + 111,24*n^2 - 12*n*d
          + 6*d - 6)
          + 1/( - 4 + 4*n)*m0denJ1010X11100( - 1 + n)*n*rat(1,1)
          + 1/( - 4 + 4*n)*m0denJ1010X11100( - 1 + n)*rat(-2,1)
          )

       + delta_( - 3 + n) * (
          + m0denJ0011X10110(0)*rat( - 6*d^2 + 31*d - 40,150*d - 1050)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat(d^4 - 13*d^3 + 54*d^2 - 92*d + 56,18*d^2
          - 162*d + 360)
          + m0denJ0011X10110(0)*rat(42*d^2 - 217*d + 280,90*d - 450)
          );

   id,ifmatch->sortme m0denJ1010X01100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m0denJ0021X00110( - 2 + n)*rat(2*n - 2*d,6*n^2 - 3*n*d - 6*n + 3*d
         )
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X00110( - 1 + n)*rat( - d^2 + 4*d - 4, - 6*n^2 + 3*n*d
          + 6*n - 3*d)
          + m0denJ0011X00110( - 1 + n)*rat( - 3*d^2 + 17*d - 24, - 6*n^2 + 3*n
         *d + 6*n - 3*d)
          + m0denJ0021X00110( - 1 + n)*rat( - 20*n + 40,6*n^2 - 3*n*d - 6*n + 
         3*d)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat(d^2 - 4*d + 4,3*d - 12)
          );

   id,ifmatch->sortme m0denJ1001X11110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat(3*n*d^2 - 12*n*d + 12*n - 3*d^3 + 
         15*d^2 - 24*d + 12,32*n^3 - 48*n^2*d + 96*n^2 + 24*n*d^2 - 96*n*d + 
         88*n - 4*d^3 + 24*d^2 - 44*d + 24)
          + m0denJ0001X11110( - 1 + n)*rat(6*n*d - 12*n - 9*d^2 + 36*d - 36,32
         *n^2 - 32*n*d + 80*n + 8*d^2 - 40*d + 48)
          + m0denJ1001X11110( - 1 + n)*rat(n - 2*d + 5,4*n - 2*d + 6)
          );

   id,ifmatch->sortme m0denJ1001X01110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,8*n^2 - 8
         *n*d + 12*n + 2*d^2 - 6*d + 4)
          + m0denJ0001X11110( - 1 + n)*rat( - 2*n*d + 6*n + d - 2,8*n^2 - 4*n*
         d + 4*d - 8)
          + 1/( - 4 + 4*n)*m0denJ1001X01110( - 1 + n)*n*rat(1,1)
          + 1/( - 4 + 4*n)*m0denJ1001X01110( - 1 + n)*rat(-2,1)
          );

   id,ifmatch->sortme m0denJ1000X01100(n?{>1}) = 0;

   id,ifmatch->sortme m0denJ0021X10110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m0denJ0021X10110( - 3 + n)*rat(4*n^2 - 12*n*d + 18*n + 9*d^2 - 27*
         d + 14,144*n^2 - 72*n*d + 36*d - 36)
          )

       + theta_( - 3 + n) * (
          + m0denJ0001X01110( - 2 + n)*rat( - 2*n*d^2 + 8*n*d - 8*n + 3*d^3 - 
         19*d^2 + 40*d - 28, - 144*n^2 + 72*n*d - 36*d + 36)
          + m0denJ0011X10110( - 2 + n)*rat( - 6*n*d^2 + 34*n*d - 48*n + 9*d^3
          - 72*d^2 + 191*d - 168, - 144*n^2 + 72*n*d - 36*d + 36)
          + m0denJ0021X00110( - 2 + n)*rat(8*n^2*d^2 - 32*n^2*d + 32*n^2 - 12*
         n*d^3 + 56*n*d^2 - 72*n*d + 16*n + 3*d^4 - 14*d^3 + 17*d^2 - 2*d,144*
         n^4 - 144*n^3*d - 144*n^3 + 36*n^2*d^2 + 180*n^2*d - 36*n^2 - 54*n*
         d^2 - 18*n*d + 36*n + 18*d^2 - 18*d)
          + m0denJ0021X10110( - 2 + n)*rat( - 28*n^2 + 60*n*d - 82*n - 27*d^2
          + 39*d + 56,72*n^2 - 36*n*d + 18*d - 18)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat(2*n*d^2 - 8*n*d + 8*n - 3*d^3 + 19*
         d^2 - 40*d + 28, - 36*n^2 + 18*n*d - 9*d + 9)
          + m0denJ0001X00110( - 1 + n)*rat( - 4*n^2*d^3 + 24*n^2*d^2 - 48*n^2*
         d + 32*n^2 + 12*n*d^4 - 98*n*d^3 + 300*n*d^2 - 408*n*d + 208*n - 3*
         d^5 + 23*d^4 - 60*d^3 + 48*d^2 + 32*d - 48,288*n^4 - 288*n^3*d - 288*
         n^3 + 72*n^2*d^2 + 360*n^2*d - 72*n^2 - 108*n*d^2 - 36*n*d + 72*n + 
         36*d^2 - 36*d)
          + m0denJ0011X10110( - 1 + n)*rat(24*n*d^2 - 136*n*d + 192*n - 27*d^3
          + 228*d^2 - 641*d + 600, - 144*n^2 + 72*n*d - 36*d + 36)
          + m0denJ0011X00110( - 1 + n)*rat( - 12*n^2*d^3 + 92*n^2*d^2 - 232*
         n^2*d + 192*n^2 + 36*n*d^4 - 354*n*d^3 + 1294*n*d^2 - 2084*n*d + 1248
         *n - 9*d^5 + 84*d^4 - 271*d^3 + 296*d^2 + 108*d - 288,288*n^4 - 288*
         n^3*d - 288*n^3 + 72*n^2*d^2 + 360*n^2*d - 72*n^2 - 108*n*d^2 - 36*n*
         d + 72*n + 36*d^2 - 36*d)
          + m0denJ0021X00110( - 1 + n)*rat( - 20*n^2*d^2 + 80*n^2*d - 80*n^2
          + 90*n*d^2 - 400*n*d + 440*n - 15*d^3 + 45*d^2 + 30*d - 120,72*n^4
          - 72*n^3*d - 72*n^3 + 18*n^2*d^2 + 90*n^2*d - 18*n^2 - 27*n*d^2 - 9*
         n*d + 18*n + 9*d^2 - 9*d)
          + m0denJ0021X10110( - 1 + n)*rat(196*n^2 - 276*n*d + 338*n + 27*d^2
          + 279*d - 786,144*n^2 - 72*n*d + 36*d - 36)
          )

       + delta_( - 3 + n) * (
          + m0denJ0011X10110(0)*rat( - 18*d^3 + 171*d^2 - 523*d + 520,900*d - 
         6300)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat(3*d^5 - 47*d^4 + 272*d^3 - 744*d^2 + 976*d
          - 496,108*d^2 - 972*d + 2160)
          + m0denJ0011X10110(0)*rat(126*d^3 - 1113*d^2 + 3227*d - 3080,540*d
          - 2700)
          );

   id,ifmatch->sortme m0denJ0021X00110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m0denJ0021X00110( - 2 + n)*rat( - 2*n^2 + 5*n*d - 6*n - 3*d^2 + 6*
         d,18*n^2 - 9*n*d - 18*n + 9*d)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X00110( - 1 + n)*rat(2*n*d^2 - 8*n*d + 8*n - 3*d^3 + 18*
         d^2 - 36*d + 24, - 36*n^2 + 18*n*d + 36*n - 18*d)
          + m0denJ0011X00110( - 1 + n)*rat(6*n*d^2 - 34*n*d + 48*n - 9*d^3 + 
         69*d^2 - 174*d + 144, - 36*n^2 + 18*n*d + 36*n - 18*d)
          + m0denJ0021X00110( - 1 + n)*rat(20*n^2 - 30*n*d + 20*n + 60*d - 120
         ,18*n^2 - 9*n*d - 18*n + 9*d)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat(3*d^3 - 22*d^2 + 52*d - 40,18*d - 72)
          );

   id,ifmatch->sortme m0denJ0011X11110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m0denJ0021X10110( - 3 + n)*rat(20*n^2 - 42*n*d + 42*n + 18*d^2 - 
         45*d + 22,576*n^3 - 576*n^2*d + 576*n^2 + 144*n*d^2 - 144*n*d - 144*n
          - 72*d^2 + 216*d - 144)
          )

       + theta_( - 3 + n) * (
          + m0denJ0001X01110( - 2 + n)*rat(10*n*d^2 - 40*n*d + 40*n - 6*d^3 + 
         35*d^2 - 68*d + 44,576*n^3 - 576*n^2*d + 576*n^2 + 144*n*d^2 - 144*n*
         d - 144*n - 72*d^2 + 216*d - 144)
          + m0denJ0011X10110( - 2 + n)*rat(30*n*d^2 - 170*n*d + 240*n - 18*d^3
          + 135*d^2 - 331*d + 264,576*n^3 - 576*n^2*d + 576*n^2 + 144*n*d^2 - 
         144*n*d - 144*n - 72*d^2 + 216*d - 144)
          + m0denJ0021X00110( - 2 + n)*rat( - 20*n^3*d + 40*n^3 + 54*n^2*d^2
          - 142*n^2*d + 68*n^2 - 36*n*d^3 + 123*n*d^2 - 112*n*d + 20*n + 6*d^4
          - 22*d^3 + 19*d^2 + 2*d,576*n^5 - 864*n^4*d + 432*n^3*d^2 + 432*n^3*
         d - 720*n^3 - 72*n^2*d^3 - 432*n^2*d^2 + 720*n^2*d + 108*n*d^3 - 108*
         n*d^2 - 216*n*d + 144*n - 36*d^3 + 108*d^2 - 72*d)
          + m0denJ0021X10110( - 2 + n)*rat( - 140*n^2 + 174*n*d - 74*n - 54*
         d^2 + 51*d + 88,288*n^3 - 288*n^2*d + 288*n^2 + 72*n*d^2 - 72*n*d - 
         72*n - 36*d^2 + 108*d - 72)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat( - 10*n*d^2 + 40*n*d - 40*n + 6*d^3
          - 35*d^2 + 68*d - 44,144*n^3 - 144*n^2*d + 144*n^2 + 36*n*d^2 - 36*n
         *d - 36*n - 18*d^2 + 54*d - 36)
          + m0denJ0001X00110( - 1 + n)*rat( - 20*n^2*d^3 + 120*n^2*d^2 - 240*
         n^2*d + 160*n^2 + 18*n*d^4 - 126*n*d^3 + 324*n*d^2 - 360*n*d + 144*n
          - 3*d^5 + 20*d^4 - 42*d^3 + 12*d^2 + 56*d - 48,576*n^5 - 864*n^4*d
          + 432*n^3*d^2 + 432*n^3*d - 720*n^3 - 72*n^2*d^3 - 432*n^2*d^2 + 720
         *n^2*d + 108*n*d^3 - 108*n*d^2 - 216*n*d + 144*n - 36*d^3 + 108*d^2
          - 72*d)
          + m0denJ0011X01110( - 1 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,16*n^2 - 
         16*n*d + 24*n + 4*d^2 - 12*d + 8)
          + m0denJ0011X10110( - 1 + n)*rat( - 54*n^2*d + 144*n^2 - 66*n*d^2 + 
         509*n*d - 888*n + 54*d^3 - 447*d^2 + 1177*d - 984,576*n^3 - 576*n^2*d
          + 576*n^2 + 144*n*d^2 - 144*n*d - 144*n - 72*d^2 + 216*d - 144)
          + m0denJ0011X11100( - 1 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,32*n^2 - 
         32*n*d + 48*n + 8*d^2 - 24*d + 16)
          + m0denJ0011X11110( - 1 + n)*rat(2*n - 3*d + 6,8*n - 4*d + 8)
          + m0denJ0011X00110( - 1 + n)*rat( - 60*n^2*d^3 + 460*n^2*d^2 - 1160*
         n^2*d + 960*n^2 + 54*n*d^4 - 468*n*d^3 + 1458*n*d^2 - 1908*n*d + 864*
         n - 9*d^5 + 75*d^4 - 202*d^3 + 122*d^2 + 252*d - 288,576*n^5 - 864*
         n^4*d + 432*n^3*d^2 + 432*n^3*d - 720*n^3 - 72*n^2*d^3 - 432*n^2*d^2
          + 720*n^2*d + 108*n*d^3 - 108*n*d^2 - 216*n*d + 144*n - 36*d^3 + 108
         *d^2 - 72*d)
          + m0denJ0021X00110( - 1 + n)*rat(100*n^3*d - 200*n^3 - 70*n^2*d^2 - 
         90*n^2*d + 460*n^2 + 285*n*d^2 - 880*n*d + 620*n - 30*d^3 + 45*d^2 + 
         150*d - 240,288*n^5 - 432*n^4*d + 216*n^3*d^2 + 216*n^3*d - 360*n^3
          - 36*n^2*d^3 - 216*n^2*d^2 + 360*n^2*d + 54*n*d^3 - 54*n*d^2 - 108*n
         *d + 72*n - 18*d^3 + 54*d^2 - 36*d)
          + m0denJ0021X10110( - 1 + n)*rat(490*n^2 - 330*n*d - 115*n + 27*d^2
          + 306*d - 633,288*n^3 - 288*n^2*d + 288*n^2 + 72*n*d^2 - 72*n*d - 72
         *n - 36*d^2 + 108*d - 72)
          )

       + delta_( - 3 + n) * (
          + m0denJ0011X10110(0)*rat(36*d^3 - 432*d^2 + 1511*d - 1640,1800*d^2
          - 27000*d + 100800)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat( - 3*d^5 + 56*d^4 - 374*d^3 + 1140*d^2 - 
         1624*d + 880,108*d^3 - 1620*d^2 + 7992*d - 12960)
          + m0denJ0011X10110(0)*rat( - 252*d^3 + 2604*d^2 - 8407*d + 8680,1080
         *d^2 - 11880*d + 32400)
          );

   id,ifmatch->sortme m0denJ0011X11100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m0denJ0011X11100( - 1 + n)*rat(n - d + 1,4*n - 2*d + 2)
          );

   id,ifmatch->sortme m0denJ0011X10110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m0denJ0021X10110( - 3 + n)*rat( - 2*n + 3*d - 2,24*n^2 - 12*n*d + 
         6*d - 6)
          )

       + theta_( - 3 + n) * (
          + m0denJ0001X01110( - 2 + n)*rat(d^2 - 4*d + 4, - 24*n^2 + 12*n*d - 
         6*d + 6)
          + m0denJ0011X10110( - 2 + n)*rat(3*d^2 - 17*d + 24, - 24*n^2 + 12*n*
         d - 6*d + 6)
          + m0denJ0021X00110( - 2 + n)*rat(2*n^2*d - 4*n^2 - 4*n*d^2 + 9*n*d
          - 2*n + d^3 - 2*d^2,24*n^4 - 24*n^3*d - 24*n^3 + 6*n^2*d^2 + 30*n^2*
         d - 6*n^2 - 9*n*d^2 - 3*n*d + 6*n + 3*d^2 - 3*d)
          + m0denJ0021X10110( - 2 + n)*rat(14*n - 9*d - 8,12*n^2 - 6*n*d + 3*d
          - 3)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat( - 2*d^2 + 8*d - 8, - 12*n^2 + 6*n*
         d - 3*d + 3)
          + m0denJ0001X00110( - 1 + n)*rat(4*n*d^3 - 24*n*d^2 + 48*n*d - 32*n
          - d^4 + 5*d^3 - 6*d^2 - 4*d + 8,48*n^4 - 48*n^3*d - 48*n^3 + 12*n^2*
         d^2 + 60*n^2*d - 12*n^2 - 18*n*d^2 - 6*n*d + 12*n + 6*d^2 - 6*d)
          + m0denJ0011X10110( - 1 + n)*rat(6*n^2 - 6*n*d + 3*n + 9*d^2 - 53*d
          + 84,24*n^2 - 12*n*d + 6*d - 6)
          + m0denJ0011X00110( - 1 + n)*rat(12*n*d^3 - 92*n*d^2 + 232*n*d - 192
         *n - 3*d^4 + 20*d^3 - 35*d^2 - 10*d + 48,48*n^4 - 48*n^3*d - 48*n^3
          + 12*n^2*d^2 + 60*n^2*d - 12*n^2 - 18*n*d^2 - 6*n*d + 12*n + 6*d^2
          - 6*d)
          + m0denJ0021X00110( - 1 + n)*rat( - 20*n^2*d + 40*n^2 + 70*n*d - 140
         *n - 10*d^2 + 40,24*n^4 - 24*n^3*d - 24*n^3 + 6*n^2*d^2 + 30*n^2*d - 
         6*n^2 - 9*n*d^2 - 3*n*d + 6*n + 3*d^2 - 3*d)
          + m0denJ0021X10110( - 1 + n)*rat( - 80*n + 9*d + 111,24*n^2 - 12*n*d
          + 6*d - 6)
          )

       + delta_( - 3 + n) * (
          + m0denJ0011X10110(0)*rat( - 6*d^2 + 31*d - 40,150*d - 1050)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat(d^4 - 13*d^3 + 54*d^2 - 92*d + 56,18*d^2
          - 162*d + 360)
          + m0denJ0011X10110(0)*rat(42*d^2 - 217*d + 280,90*d - 450)
          );

   id,ifmatch->sortme m0denJ0011X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m0denJ0021X00110( - 2 + n)*rat( - 5*n^2 + 8*n*d - 3*n - 3*d^2 + 3*
         d,36*n^3 - 36*n^2*d - 18*n^2 + 9*n*d^2 + 27*n*d - 18*n - 9*d^2 + 9*d)
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X00110( - 1 + n)*rat( - 5*n*d^2 + 20*n*d - 20*n + 3*d^3
          - 15*d^2 + 24*d - 12,72*n^3 - 72*n^2*d - 36*n^2 + 18*n*d^2 + 54*n*d
          - 36*n - 18*d^2 + 18*d)
          + m0denJ0011X01110( - 1 + n)*rat(n - d + 1,4*n - 2*d + 2)
          + m0denJ0011X00110( - 1 + n)*rat( - 15*n*d^2 + 85*n*d - 120*n + 9*
         d^3 - 60*d^2 + 123*d - 72,72*n^3 - 72*n^2*d - 36*n^2 + 18*n*d^2 + 54*
         n*d - 36*n - 18*d^2 + 18*d)
          + m0denJ0021X00110( - 1 + n)*rat(50*n^2 - 30*n*d - 70*n + 60*d - 60,
         36*n^3 - 36*n^2*d - 18*n^2 + 9*n*d^2 + 27*n*d - 18*n - 9*d^2 + 9*d)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat( - 3*d^3 + 25*d^2 - 64*d + 52,18*d^2 - 162
         *d + 360)
          );

   id,ifmatch->sortme m0denJ0011X01100(n?{>1}) = 0;

   id,ifmatch->sortme m0denJ0011X00110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m0denJ0021X00110( - 2 + n)*rat(2*n - 2*d,6*n^2 - 3*n*d - 6*n + 3*d
         )
          )

       + theta_( - 2 + n) * (
          + m0denJ0001X00110( - 1 + n)*rat( - d^2 + 4*d - 4, - 6*n^2 + 3*n*d
          + 6*n - 3*d)
          + m0denJ0011X00110( - 1 + n)*rat( - 3*d^2 + 17*d - 24, - 6*n^2 + 3*n
         *d + 6*n - 3*d)
          + m0denJ0021X00110( - 1 + n)*rat( - 20*n + 40,6*n^2 - 3*n*d - 6*n + 
         3*d)
          )

       + delta_( - 2 + n) * (
          + m0denJ0001X01110(0)*rat(d^2 - 4*d + 4,3*d - 12)
          );

   id,ifmatch->sortme m0denJ0001X11110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,8*n^2 - 8
         *n*d + 12*n + 2*d^2 - 6*d + 4)
          + m0denJ0001X11110( - 1 + n)*rat(2*n - 3*d + 6,8*n - 4*d + 8)
          );

   id,ifmatch->sortme m0denJ0001X01110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m0denJ0001X01110( - 1 + n)*rat(n - d + 1,4*n - 2*d + 2)
          );

   id,ifmatch->sortme m0denJ0001X00110(n?{>1}) = 0;

********************************************************************************                 

goto endrec;

la sortme;
$irep = 0;
la endrec;

ModuleOption,minimum,$irep;

.sort:rec1d-m0-`$repcount++';

#redefine irep "`$irep'"
#enddo
#endprocedure
*--#] rec1dm0 :

*--#[ rec1dm1 :
#procedure rec1dm1
#$repcount = 1;

#do irep = 1,1
#$irep = 1;
                
********************************************************************************             

   id,ifmatch->sortme m1denJ2010X11100(n?{>1}) =

       + theta_( - 4 + n) * (
          + m1denJ0021X10110( - 3 + n)*rat(n*d - 4*n + 5*d - 14,48*n^2 - 48*n*
         d + 48*n + 48*d - 96)
          + 1/( - 24 + 24*n)*m1denJ2010X11100( - 3 + n)*n*rat(-1,1)
          + 1/( - 24 + 24*n)*m1denJ2010X11100( - 3 + n)*rat(4,1)
          )

       + theta_( - 3 + n) * (
          + m1denJ0001X01110( - 2 + n)*rat(d^3 - 5*d^2 + 8*d - 4, - 48*n^2 + 
         48*n*d - 48*n - 48*d + 96)
          + m1denJ0011X10110( - 2 + n)*rat(9*n*d^2 - 51*n*d + 72*n - 12*d^3 + 
         83*d^2 - 181*d + 120, - 48*n^2 + 48*n*d - 48*n - 48*d + 96)
          + m1denJ0021X00110( - 2 + n)*rat(d - 2, - 24*n + 24*d - 48)
          + m1denJ0021X10110( - 2 + n)*rat(11*n*d - 22*n - 18*d^2 + 51*d - 2,
         48*n^2 - 48*n*d + 48*n + 48*d - 96)
          + 1/( - 12 + 12*n)*m1denJ2010X11100( - 2 + n)*n*rat(5,1)
          + 1/( - 12 + 12*n)*m1denJ2010X11100( - 2 + n)*rat(-15,1)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X01110( - 1 + n)*rat( - d^3 + 5*d^2 - 8*d + 4, - 48*n^2
          + 48*n*d - 48*n - 48*d + 96)
          + m1denJ0001X00110( - 1 + n)*rat( - d^3 + 6*d^2 - 12*d + 8, - 48*n^2
          + 48*n*d - 48*n - 48*d + 96)
          + m1denJ0011X10110( - 1 + n)*rat( - 27*n*d^2 + 153*n*d - 216*n + 24*
         d^3 - 169*d^2 + 379*d - 264, - 48*n^2 + 48*n*d - 48*n - 48*d + 96)
          + m1denJ0011X00110( - 1 + n)*rat(3*d^3 - 23*d^2 + 58*d - 48, - 24*
         n^2 + 24*n*d - 24*n - 24*d + 48)
          + m1denJ0021X00110( - 1 + n)*rat( - d + 2, - 3*n + 3*d - 6)
          + m1denJ0021X10110( - 1 + n)*rat( - 31*n*d + 81*n + 27*d^2 - 109*d
          + 96,24*n^2 - 24*n*d + 24*n + 24*d - 48)
          + 1/( - 24 + 24*n)*m1denJ2010X11100( - 1 + n)*n*rat(-13,1)
          + 1/( - 24 + 24*n)*m1denJ2010X11100( - 1 + n)*rat(26,1)
          )

       + delta_( - 3 + n) * (
          + m1denJ0011X10110(0)*rat( - 6*d^2 + 31*d - 40,240*d - 1200)
          )

       + delta_( - 2 + n) * (
          + m1denJ0001X01110(0)*rat(d^3 - 6*d^2 + 12*d - 8,48*d - 192)
          + m1denJ0011X10110(0)*rat( - 24*d^3 + 226*d^2 - 687*d + 680,240*d - 
         960)
          );

   id,ifmatch->sortme m1denJ2010X01100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0021X00110( - 2 + n)*rat( - d + 3,4*n^2 - 4*n*d + 4*n + 4*d
          - 8)
          + 1/( - 8 + 8*n)*m1denJ2010X01100( - 2 + n)*n*rat(1,1)
          + 1/( - 8 + 8*n)*m1denJ2010X01100( - 2 + n)*rat(-3,1)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X00110( - 1 + n)*rat( - d^3 + 4*d^2 - 4*d, - 32*n^2 + 32
         *n*d - 32*n - 32*d + 64)
          + m1denJ0011X00110( - 1 + n)*rat( - 9*n*d^2 + 51*n*d - 72*n + 9*d^3
          - 60*d^2 + 123*d - 72, - 16*n^2 + 16*n*d - 16*n - 16*d + 32)
          + m1denJ0021X00110( - 1 + n)*rat( - 9*n*d + 27*n + 9*d^2 - 35*d + 24
         ,8*n^2 - 8*n*d + 8*n + 8*d - 16)
          + 1/( - 8 + 8*n)*m1denJ2010X01100( - 1 + n)*n*rat(-7,1)
          + 1/( - 8 + 8*n)*m1denJ2010X01100( - 1 + n)*rat(14,1)
          )

       + delta_( - 2 + n) * (
          + m1denJ0001X01110(0)*rat(d^2 - 4*d + 4,16*d - 64)
          );

   id,ifmatch->sortme m1denJ1011X11110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m1denJ1011X11110( - 3 + n)*rat( - 2*n + 3*d - 4,12*n - 12)
          )

       + theta_( - 3 + n) * (
          + m1denJ0001X11110( - 2 + n)*rat(d - 2,6*n - 6)
          + m1denJ1001X11110( - 2 + n)*rat( - d + 3,6*n - 6)
          + m1denJ1011X01110( - 2 + n)*rat( - d + 2,6*n - 6)
          + m1denJ1011X11110( - 2 + n)*rat(2*n - 2*d + 3,3*n - 3)
          + 1/( - 3 + 3*n)*m1denJ0021X10110( - 2 + n)*rat(-1,1)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X11110( - 1 + n)*rat( - d + 2,6*n - 6)
          + m1denJ0011X10110( - 1 + n)*rat(3*d - 8,6*n - 6)
          + m1denJ0011X11100( - 1 + n)*rat( - d + 2,6*n - 6)
          + m1denJ0011X11110( - 1 + n)*rat(2*d - 6,3*n - 3)
          + m1denJ1001X11110( - 1 + n)*rat( - d + 3,6*n - 6)
          + m1denJ1011X01110( - 1 + n)*rat(d - 2,3*n - 3)
          + m1denJ1011X11110( - 1 + n)*rat( - 2*n + d - 8,12*n - 12)
          + 1/( - 3 + 3*n)*m1denJ0021X10110( - 1 + n)*rat(-2,1)
          )

       + delta_( - 3 + n) * (
          + m1denJ1011X11110(0)*rat(3*d - 10,24)
          )

       + delta_( - 2 + n) * (
          + m1denJ0011X01110(0)*rat(5*d^2 - 20*d + 20,96*d - 288)
          + m1denJ0011X10110(0)*rat( - 4*d^2 + 24*d - 35,48*d - 144)
          + m1denJ0011X11110(0)*rat(3*d - 8,16)
          + m1denJ0031X10110(0)*rat(25,48*d - 144)
          + m1denJ1001X11110(0)*rat( - d + 3,6)
          + m1denJ1011X11110(0)*rat( - 3*d + 10,6)
          + m1denJ2011X11110(0)*rat( - 1,2)
          );

   id,ifmatch->sortme m1denJ1011X01110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m1denJ1011X01110( - 3 + n)*rat( - n + d,6*n - 6)
          )

       + theta_( - 3 + n) * (
          + m1denJ0001X01110( - 2 + n)*rat(d - 2,6*n - 6)
          + m1denJ0001X11110( - 2 + n)*rat( - d + 3,6*n - 6)
          + m1denJ1011X01110( - 2 + n)*rat(8*n - 7*d + 4,12*n - 12)
          + 1/( - 3 + 3*n)*m1denJ0021X00110( - 2 + n)*rat(-1,1)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X01110( - 1 + n)*rat( - d + 2,6*n - 6)
          + m1denJ0001X11110( - 1 + n)*rat( - d + 3,6*n - 6)
          + m1denJ0011X01100( - 1 + n)*rat( - d + 2,6*n - 6)
          + m1denJ0011X01110( - 1 + n)*rat(2*d - 6,3*n - 3)
          + m1denJ0011X00110( - 1 + n)*rat(3*d - 8,6*n - 6)
          + m1denJ1011X01110( - 1 + n)*rat( - 2*n + 3*d - 8,12*n - 12)
          + 1/( - 3 + 3*n)*m1denJ0021X00110( - 1 + n)*rat(-2,1)
          )

       + delta_( - 3 + n) * (
          + m1denJ1011X01110(0)*rat(d - 3,12)
          )

       + delta_( - 2 + n) * (
          + m1denJ0001X11110(0)*rat( - 3*d + 8,12)
          + m1denJ0011X01110(0)*rat(d - 2,3)
          + m1denJ1011X01110(0)*rat( - d + 3,2)
          );

   id,ifmatch->sortme m1denJ1010X11100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0011X10110( - 2 + n)*rat( - 2*d + 6,3*n - 3)
          + 1/( - 3 + 3*n)*m1denJ1010X11100( - 2 + n)*n*rat(1,1)
          + 1/( - 3 + 3*n)*m1denJ1010X11100( - 2 + n)*rat(-3,1)
          + 1/( - 1 + n)*m1denJ0021X10110( - 2 + n)*rat(1,1)
          )

       + theta_( - 2 + n) * (
          + m1denJ0011X10110( - 1 + n)*rat(4*d - 12,3*n - 3)
          + m1denJ0011X00110( - 1 + n)*rat(d - 2,3*n - 3)
          + 1/( - 3 + 3*n)*m1denJ1010X11100( - 1 + n)*n*rat(-2,1)
          + 1/( - 3 + 3*n)*m1denJ1010X11100( - 1 + n)*rat(4,1)
          + 1/( - 1 + n)*m1denJ0021X10110( - 1 + n)*rat(-3,1)
          )

       + delta_( - 2 + n) * (
          + m1denJ0011X10110(0)*rat( - 4*d + 10,15)
          );

   id,ifmatch->sortme m1denJ1010X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m1denJ0011X00110( - 1 + n)*rat(3*d - 8,2*n - 2)
          + 1/( - 1 + n)*m1denJ0021X00110( - 1 + n)*rat(-3,1)
          + 1/( - 1 + n)*m1denJ1010X01100( - 1 + n)*n*rat(-1,1)
          + 1/( - 1 + n)*m1denJ1010X01100( - 1 + n)*rat(2,1)
          );

   id,ifmatch->sortme m1denJ1001X11110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ1001X11110( - 2 + n)*rat(n - 2*d + 4,3*n - 3)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X11110( - 1 + n)*rat(d - 2,n - 1)
          + 1/( - 3 + 3*n)*m1denJ1001X11110( - 1 + n)*n*rat(-2,1)
          + 1/( - 3 + 3*n)*m1denJ1001X11110( - 1 + n)*rat(1,1)
          )

       + delta_( - 2 + n) * (
          + m1denJ1001X11110(0)*rat( - 2*d + 6,3)
          );

   id,ifmatch->sortme m1denJ1001X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0001X11110( - 2 + n)*rat( - 3*d + 10,6*n - 6)
          + 1/( - 3 + 3*n)*m1denJ1001X01110( - 2 + n)*n*rat(1,1)
          + 1/( - 3 + 3*n)*m1denJ1001X01110( - 2 + n)*rat(-3,1)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X01110( - 1 + n)*rat(2*d - 4,3*n - 3)
          + m1denJ0001X11110( - 1 + n)*rat(d - 6,6*n - 6)
          + 1/( - 3 + 3*n)*m1denJ1001X01110( - 1 + n)*n*rat(-2,1)
          + 1/( - 3 + 3*n)*m1denJ1001X01110( - 1 + n)*rat(4,1)
          )

       + delta_( - 2 + n) * (
          + m1denJ0001X11110(0)*rat( - 3*d + 8,6)
          );

   id,ifmatch->sortme m1denJ1000X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m1denJ0001X00110( - 1 + n)*rat(d - 2,2*n - 2)
          + 1/( - 1 + n)*m1denJ1000X01100( - 1 + n)*n*rat(-1,1)
          + 1/( - 1 + n)*m1denJ1000X01100( - 1 + n)*rat(2,1)
          );

   id,ifmatch->sortme m1denJ0021X10110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m1denJ0021X10110( - 3 + n)*rat( - 2*n + 3*d - 2,48*n - 48*d + 96)
          )

       + theta_( - 3 + n) * (
          + m1denJ0001X01110( - 2 + n)*rat(d^3 - 5*d^2 + 8*d - 4, - 48*n^2 + 
         48*n*d - 48*n - 48*d + 96)
          + m1denJ0011X10110( - 2 + n)*rat(9*n*d^2 - 51*n*d + 72*n - 12*d^3 + 
         83*d^2 - 181*d + 120, - 48*n^2 + 48*n*d - 48*n - 48*d + 96)
          + m1denJ0021X00110( - 2 + n)*rat(d - 2, - 24*n + 24*d - 48)
          + m1denJ0021X10110( - 2 + n)*rat(20*n^2 - 9*n*d - 42*n - 18*d^2 + 
         111*d - 122,48*n^2 - 48*n*d + 48*n + 48*d - 96)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X01110( - 1 + n)*rat( - d^3 + 5*d^2 - 8*d + 4, - 48*n^2
          + 48*n*d - 48*n - 48*d + 96)
          + m1denJ0001X00110( - 1 + n)*rat( - d^3 + 6*d^2 - 12*d + 8, - 48*n^2
          + 48*n*d - 48*n - 48*d + 96)
          + m1denJ0011X10110( - 1 + n)*rat( - 27*n*d^2 + 153*n*d - 216*n + 24*
         d^3 - 169*d^2 + 379*d - 264, - 48*n^2 + 48*n*d - 48*n - 48*d + 96)
          + m1denJ0011X00110( - 1 + n)*rat(3*d^3 - 23*d^2 + 58*d - 48, - 24*
         n^2 + 24*n*d - 24*n - 24*d + 48)
          + m1denJ0021X00110( - 1 + n)*rat( - d + 2, - 3*n + 3*d - 6)
          + m1denJ0021X10110( - 1 + n)*rat( - 13*n^2 - 18*n*d + 81*n + 27*d^2
          - 135*d + 148,24*n^2 - 24*n*d + 24*n + 24*d - 48)
          )

       + delta_( - 3 + n) * (
          + m1denJ0011X10110(0)*rat( - 6*d^2 + 31*d - 40,240*d - 1200)
          )

       + delta_( - 2 + n) * (
          + m1denJ0001X01110(0)*rat(d^3 - 6*d^2 + 12*d - 8,48*d - 192)
          + m1denJ0011X10110(0)*rat( - 24*d^3 + 226*d^2 - 687*d + 680,240*d - 
         960)
          );

   id,ifmatch->sortme m1denJ0021X00110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0021X00110( - 2 + n)*rat(n - d,8*n - 8*d + 16)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X00110( - 1 + n)*rat( - d^3 + 4*d^2 - 4*d, - 32*n^2 + 32
         *n*d - 32*n - 32*d + 64)
          + m1denJ0011X00110( - 1 + n)*rat( - 9*n*d^2 + 51*n*d - 72*n + 9*d^3
          - 60*d^2 + 123*d - 72, - 16*n^2 + 16*n*d - 16*n - 16*d + 32)
          + m1denJ0021X00110( - 1 + n)*rat( - 7*n^2 - 2*n*d + 27*n + 9*d^2 - 
         49*d + 52,8*n^2 - 8*n*d + 8*n + 8*d - 16)
          )

       + delta_( - 2 + n) * (
          + m1denJ0001X01110(0)*rat(d^2 - 4*d + 4,16*d - 64)
          );

   id,ifmatch->sortme m1denJ0011X11110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0011X11110( - 2 + n)*rat(2*n - 3*d + 4,6*n - 6)
          )

       + theta_( - 2 + n) * (
          + m1denJ0011X01110( - 1 + n)*rat(d - 2,3*n - 3)
          + m1denJ0011X10110( - 1 + n)*rat( - 3*d + 8,6*n - 6)
          + m1denJ0011X11100( - 1 + n)*rat(d - 2,6*n - 6)
          + m1denJ0011X11110( - 1 + n)*rat( - 4*n + d + 2,6*n - 6)
          + 1/( - 3 + 3*n)*m1denJ0021X10110( - 1 + n)*rat(5,1)
          )

       + delta_( - 2 + n) * (
          + m1denJ0011X11110(0)*rat( - 3*d + 8,6)
          );

   id,ifmatch->sortme m1denJ0011X11100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0011X11100( - 2 + n)*rat(n - d,3*n - 3)
          )

       + theta_( - 2 + n) * (
          + m1denJ0011X01100( - 1 + n)*rat(d - 2,3*n - 3)
          + m1denJ0011X11100( - 1 + n)*rat( - 2*n + d + 1,3*n - 3)
          )

       + delta_( - 2 + n) * (
          + m1denJ0011X01110(0)*rat( - d + 2,3)
          );

   id,ifmatch->sortme m1denJ0011X10110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0011X10110( - 2 + n)*rat(n - 2*d + 3,3*n - 3)
          + 1/( - 1 + n)*m1denJ0021X10110( - 2 + n)*rat(1,1)
          )

       + theta_( - 2 + n) * (
          + m1denJ0011X10110( - 1 + n)*rat( - 2*n + 4*d - 8,3*n - 3)
          + m1denJ0011X00110( - 1 + n)*rat(d - 2,3*n - 3)
          + 1/( - 1 + n)*m1denJ0021X10110( - 1 + n)*rat(-3,1)
          )

       + delta_( - 2 + n) * (
          + m1denJ0011X10110(0)*rat( - 4*d + 10,15)
          );

   id,ifmatch->sortme m1denJ0011X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0011X01110( - 2 + n)*rat(n - d,3*n - 3)
          )

       + theta_( - 2 + n) * (
          + m1denJ0011X01100( - 1 + n)*rat(d - 2,6*n - 6)
          + m1denJ0011X01110( - 1 + n)*rat( - 2*n + d + 1,3*n - 3)
          + m1denJ0011X00110( - 1 + n)*rat( - 3*d + 8,6*n - 6)
          + 1/( - 3 + 3*n)*m1denJ0021X00110( - 1 + n)*rat(5,1)
          )

       + delta_( - 2 + n) * (
          + m1denJ0011X01110(0)*rat( - d + 2,3)
          );

   id,ifmatch->sortme m1denJ0011X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m1denJ0011X01100( - 1 + n)*rat( - 2*n + d + 2,2*n - 2)
          );

   id,ifmatch->sortme m1denJ0011X00110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m1denJ0011X00110( - 1 + n)*rat( - 2*n + 3*d - 4,2*n - 2)
          + 1/( - 1 + n)*m1denJ0021X00110( - 1 + n)*rat(-3,1)
          );

   id,ifmatch->sortme m1denJ0001X11110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0001X11110( - 2 + n)*rat(2*n - 3*d + 4,6*n - 6)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X01110( - 1 + n)*rat(2*d - 4,3*n - 3)
          + m1denJ0001X11110( - 1 + n)*rat( - 4*n + d + 2,6*n - 6)
          )

       + delta_( - 2 + n) * (
          + m1denJ0001X11110(0)*rat( - 3*d + 8,6)
          );

   id,ifmatch->sortme m1denJ0001X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m1denJ0001X01110( - 2 + n)*rat(n - d,3*n - 3)
          )

       + theta_( - 2 + n) * (
          + m1denJ0001X01110( - 1 + n)*rat( - 2*n + d + 1,3*n - 3)
          + m1denJ0001X00110( - 1 + n)*rat(d - 2,3*n - 3)
          )

       + delta_( - 2 + n) * (
          + m1denJ0001X01110(0)*rat( - d + 2,3)
          );

   id,ifmatch->sortme m1denJ0001X00110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m1denJ0001X00110( - 1 + n)*rat( - 2*n + d + 2,2*n - 2)
          );    

********************************************************************************                 

goto endrec;

la sortme;
$irep = 0;
la endrec;

ModuleOption,minimum,$irep;

.sort:rec1d-m1-`$repcount++';

#redefine irep "`$irep'"
#enddo
#endprocedure
*--#] rec1dm1 :

*--#[ rec1dm3 :
#procedure rec1dm3
#$repcount = 1;

#do irep = 1,1
#$irep = 1;
                
********************************************************************************                 
   id,ifmatch->sortme m3denJ2010X11100(n?{>1}) =

       + theta_( - 5 + n) * (
          + m3denJ0021X10110( - 4 + n)*rat(3*d - 10,72*n - 72)
          + 1/( - 36 + 36*n)*m3denJ2010X11100( - 4 + n)*n*rat(-1,1)
          + 1/( - 36 + 36*n)*m3denJ2010X11100( - 4 + n)*rat(5,1)
          )

       + theta_( - 4 + n) * (
          + m3denJ0001X01110( - 3 + n)*rat( - d^2 + 4*d - 4,72*n - 72)
          + m3denJ0011X10110( - 3 + n)*rat( - 3*d^2 + 17*d - 24,72*n - 72)
          + m3denJ0021X00110( - 3 + n)*rat( - d + 2,36*n - 36)
          + m3denJ0021X10110( - 3 + n)*rat(9*d - 22,72*n - 72)
          + 1/( - 18 + 18*n)*m3denJ2010X11100( - 3 + n)*n*rat(1,1)
          + 1/( - 18 + 18*n)*m3denJ2010X11100( - 3 + n)*rat(-4,1)
          )

       + theta_( - 3 + n) * (
          + m3denJ0001X01110( - 2 + n)*rat( - d^2 + 4*d - 4,36*n - 36)
          + m3denJ0011X10110( - 2 + n)*rat(3*d^2 - 17*d + 24,72*n - 72)
          + m3denJ0021X00110( - 2 + n)*rat(d - 2,9*n - 9)
          + m3denJ0021X10110( - 2 + n)*rat( - 9*d + 29,18*n - 18)
          + 1/( - 36 + 36*n)*m3denJ2010X11100( - 2 + n)*n*rat(23,1)
          + 1/( - 36 + 36*n)*m3denJ2010X11100( - 2 + n)*rat(-69,1)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X01110( - 1 + n)*rat(d^2 - 4*d + 4,24*n - 24)
          + m3denJ0021X00110( - 1 + n)*rat(d - 2,3*n - 3)
          + m3denJ0021X10110( - 1 + n)*rat( - 3*d + 8,6*n - 6)
          + 1/( - 3 + 3*n)*m3denJ2010X11100( - 1 + n)*n*rat(1,1)
          + 1/( - 3 + 3*n)*m3denJ2010X11100( - 1 + n)*rat(-2,1)
          )

       + delta_( - 4 + n) * (
          + m3denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,1080)
          )

       + delta_( - 3 + n) * (
          + m3denJ0011X10110(0)*rat( - 6*d^2 + 31*d - 40,720)
          )

       + delta_( - 2 + n) * (
          + m3denJ0011X10110(0)*rat( - 2*d^2 + 11*d - 15,12)
          + m3denJ0031X10110(0)*rat(5,12)
          );

   id,ifmatch->sortme m3denJ2010X01100(n?{>1}) =

       + theta_( - 4 + n) * (
          + m3denJ0021X00110( - 3 + n)*rat( - d + 3,36*n - 36)
          + 1/( - 36 + 36*n)*m3denJ2010X01100( - 3 + n)*n*rat(1,1)
          + 1/( - 36 + 36*n)*m3denJ2010X01100( - 3 + n)*rat(-4,1)
          )

       + theta_( - 3 + n) * (
          + m3denJ0001X00110( - 2 + n)*rat(d^2 - 4*d + 4,72*n - 72)
          + m3denJ0011X00110( - 2 + n)*rat(3*d^2 - 17*d + 24,72*n - 72)
          + m3denJ0021X00110( - 2 + n)*rat( - d + 3,6*n - 6)
          + 1/( - 36 + 36*n)*m3denJ2010X01100( - 2 + n)*n*rat(-1,1)
          + 1/( - 36 + 36*n)*m3denJ2010X01100( - 2 + n)*rat(3,1)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X00110( - 1 + n)*rat(d^2 - 4*d + 4,24*n - 24)
          + 1/( - 3 + 3*n)*m3denJ2010X01100( - 1 + n)*n*rat(-2,1)
          + 1/( - 3 + 3*n)*m3denJ2010X01100( - 1 + n)*rat(4,1)
          )

       + delta_( - 3 + n) * (
          + m3denJ0001X01110(0)*rat( - d^2 + 4*d - 4,144)
          )

       + delta_( - 2 + n) * (
          + m3denJ0001X01110(0)*rat( - d^2 + 4*d - 4,24)
          );

   id,ifmatch->sortme m3denJ1011X11110(n?{>1}) =

       + theta_( - 5 + n) * (
          + m3denJ0021X10110( - 4 + n)*rat( - 4*n + 6*d,54*n^2 - 27*n*d + 27*d
          - 54)
          )

       + theta_( - 4 + n) * (
          + m3denJ0001X01110( - 3 + n)*rat(2*d^2 - 8*d + 8, - 54*n^2 + 27*n*d
          - 27*d + 54)
          + m3denJ0011X10110( - 3 + n)*rat(6*d^2 - 34*d + 48, - 54*n^2 + 27*n*
         d - 27*d + 54)
          + m3denJ0021X00110( - 3 + n)*rat( - 4*d + 8,54*n^2 - 27*n*d + 27*d
          - 54)
          + m3denJ0021X10110( - 3 + n)*rat(8*n + 18*d - 76,54*n^2 - 27*n*d + 
         27*d - 54)
          )

       + theta_( - 3 + n) * (
          + m3denJ0001X01110( - 2 + n)*rat(4*d^2 - 16*d + 16, - 54*n^2 + 27*n*
         d - 27*d + 54)
          + m3denJ0001X11110( - 2 + n)*rat( - 2*n*d + 4*n + 3*d^2 - 10*d + 8,
         18*n^2 - 9*n*d + 9*d - 18)
          + m3denJ0011X10110( - 2 + n)*rat( - 18*n*d + 48*n + 42*d^2 - 184*d
          + 192,54*n^2 - 27*n*d + 27*d - 54)
          + m3denJ0011X11100( - 2 + n)*rat(2*n*d - 4*n - 2*d^2 + 4*d,18*n^2 - 
         9*n*d + 9*d - 18)
          + m3denJ0011X11110( - 2 + n)*rat( - 8*n*d + 24*n + 12*d^2 - 52*d + 
         48,18*n^2 - 9*n*d + 9*d - 18)
          + m3denJ0021X00110( - 2 + n)*rat(16*d - 32,54*n^2 - 27*n*d + 27*d - 
         54)
          + m3denJ0021X10110( - 2 + n)*rat(92*n - 126*d + 100,54*n^2 - 27*n*d
          + 27*d - 54)
          + m3denJ1001X11110( - 2 + n)*rat(2*n*d - 6*n - 4*d^2 + 20*d - 24,6*
         n^2 - 3*n*d + 3*d - 6)
          + m3denJ1011X11110( - 2 + n)*rat(2*n - 3*d + 6,6*n - 3*d + 6)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X01110( - 1 + n)*rat(2*d^2 - 8*d + 8, - 18*n^2 + 9*n*d
          - 9*d + 18)
          + m3denJ0001X11110( - 1 + n)*rat( - 10*n*d + 20*n + 23*d^2 - 100*d
          + 108,18*n^2 - 9*n*d + 9*d - 18)
          + m3denJ0011X01100( - 1 + n)*rat( - 2*d^2 + 8*d - 8, - 18*n^2 + 9*n*
         d - 9*d + 18)
          + m3denJ0011X01110( - 1 + n)*rat(8*d^2 - 40*d + 48, - 18*n^2 + 9*n*d
          - 9*d + 18)
          + m3denJ0011X10110( - 1 + n)*rat( - 12*n*d + 32*n + 12*d^2 - 44*d + 
         32,18*n^2 - 9*n*d + 9*d - 18)
          + m3denJ0011X11100( - 1 + n)*rat(4*n*d - 8*n - 6*d^2 + 22*d - 20,18*
         n^2 - 9*n*d + 9*d - 18)
          + m3denJ0011X11110( - 1 + n)*rat( - 16*n*d + 48*n + 20*d^2 - 84*d + 
         72,18*n^2 - 9*n*d + 9*d - 18)
          + m3denJ0011X00110( - 1 + n)*rat(6*d^2 - 28*d + 32, - 18*n^2 + 9*n*d
          - 9*d + 18)
          + m3denJ0021X00110( - 1 + n)*rat(16*d - 32,18*n^2 - 9*n*d + 9*d - 18
         )
          + m3denJ0021X10110( - 1 + n)*rat(28*n - 46*d + 92,18*n^2 - 9*n*d + 9
         *d - 18)
          + m3denJ1001X11110( - 1 + n)*rat(6*n*d - 18*n - 8*d^2 + 36*d - 36,6*
         n^2 - 3*n*d + 3*d - 6)
          + m3denJ1011X01110( - 1 + n)*rat( - 2*d + 4, - 6*n + 3*d - 6)
          + m3denJ1011X11110( - 1 + n)*rat(4*n - 4*d + 4,6*n - 3*d + 6)
          )

       + delta_( - 4 + n) * (
          + m3denJ0011X10110(0)*rat( - 12*d^2 + 62*d - 80,405*d - 4050)
          )

       + delta_( - 3 + n) * (
          + m3denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,135*d - 1080)
          )

       + delta_( - 2 + n) * (
          + m3denJ0001X11110(0)*rat( - 3*d^2 + 14*d - 16,9*d - 54)
          + m3denJ0011X01110(0)*rat(2*d^2 - 8*d + 8,9*d - 54)
          + m3denJ0011X10110(0)*rat(16*d^2 - 96*d + 140,45*d - 270)
          + m3denJ0011X11110(0)*rat( - 12*d^2 + 68*d - 96,9*d - 54)
          + m3denJ0031X10110(0)*rat(-20,9*d - 54)
          + m3denJ1001X11110(0)*rat(4*d^2 - 24*d + 36,3*d - 18)
          + m3denJ1011X11110(0)*rat(3*d - 10,3*d - 18)
          );

   id,ifmatch->sortme m3denJ1011X01110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m3denJ0021X00110( - 3 + n)*rat(4*n - 4*d - 4,54*n^2 - 27*n*d + 27*
         d - 54)
          )

       + theta_( - 3 + n) * (
          + m3denJ0001X01110( - 2 + n)*rat( - 2*n*d + 4*n + 2*d^2 - 4*d,18*n^2
          - 9*n*d + 9*d - 18)
          + m3denJ0001X11110( - 2 + n)*rat(2*n*d - 6*n - 3*d^2 + 13*d - 12,6*
         n^2 - 3*n*d + 3*d - 6)
          + m3denJ0001X00110( - 2 + n)*rat( - 2*d^2 + 8*d - 8, - 54*n^2 + 27*n
         *d - 27*d + 54)
          + m3denJ0011X01110( - 2 + n)*rat( - 8*n*d + 24*n + 8*d^2 - 24*d,18*
         n^2 - 9*n*d + 9*d - 18)
          + m3denJ0011X00110( - 2 + n)*rat( - 6*d^2 + 34*d - 48, - 54*n^2 + 27
         *n*d - 27*d + 54)
          + m3denJ0021X00110( - 2 + n)*rat( - 4*n - 24*d + 84,54*n^2 - 27*n*d
          + 27*d - 54)
          + m3denJ1011X01110( - 2 + n)*rat(2*n - 2*d + 2,6*n - 3*d + 6)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X01110( - 1 + n)*rat( - 10*n*d + 20*n + 14*d^2 - 56*d + 
         56,18*n^2 - 9*n*d + 9*d - 18)
          + m3denJ0001X11110( - 1 + n)*rat(6*n*d - 18*n - 5*d^2 + 19*d - 12,6*
         n^2 - 3*n*d + 3*d - 6)
          + m3denJ0011X01100( - 1 + n)*rat( - 2*n*d + 4*n - 3*d^2 + 20*d - 28,
         18*n^2 - 9*n*d + 9*d - 18)
          + m3denJ0011X01110( - 1 + n)*rat( - 16*n*d + 48*n + 8*d^2 - 16*d - 
         24,18*n^2 - 9*n*d + 9*d - 18)
          + m3denJ0011X00110( - 1 + n)*rat(6*n*d - 16*n + 3*d^2 - 32*d + 64,18
         *n^2 - 9*n*d + 9*d - 18)
          + m3denJ0021X00110( - 1 + n)*rat( - 20*n - 22*d + 124,18*n^2 - 9*n*d
          + 9*d - 18)
          + m3denJ1011X01110( - 1 + n)*rat(4*n - d - 4,6*n - 3*d + 6)
          )

       + delta_( - 3 + n) * (
          + m3denJ0001X01110(0)*rat(d^2 - 4*d + 4,27*d - 216)
          )

       + delta_( - 2 + n) * (
          + m3denJ0001X11110(0)*rat(3*d^2 - 17*d + 24,3*d - 18)
          + m3denJ0011X01110(0)*rat( - 8*d^2 + 40*d - 48,9*d - 54)
          + m3denJ1011X01110(0)*rat(2*d - 6,3*d - 18)
          );

   id,ifmatch->sortme m3denJ1010X11100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m3denJ0011X10110( - 2 + n)*rat( - 2*d + 6,3*n - 3)
          + 1/( - 3 + 3*n)*m3denJ1010X11100( - 2 + n)*n*rat(1,1)
          + 1/( - 3 + 3*n)*m3denJ1010X11100( - 2 + n)*rat(-3,1)
          + 1/( - 1 + n)*m3denJ0021X10110( - 2 + n)*rat(1,1)
          )

       + theta_( - 2 + n) * (
          + m3denJ0011X00110( - 1 + n)*rat(d - 2,3*n - 3)
          + 1/( - 3 + 3*n)*m3denJ1010X11100( - 1 + n)*n*rat(2,1)
          + 1/( - 3 + 3*n)*m3denJ1010X11100( - 1 + n)*rat(-4,1)
          + 1/( - 1 + n)*m3denJ0021X10110( - 1 + n)*rat(-1,1)
          )

       + delta_( - 2 + n) * (
          + m3denJ0011X10110(0)*rat( - 4*d + 10,15)
          );

   id,ifmatch->sortme m3denJ1010X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m3denJ0011X00110( - 1 + n)*rat(3*d - 8,6*n - 6)
          + 1/( - 3 + 3*n)*m3denJ1010X01100( - 1 + n)*n*rat(-1,1)
          + 1/( - 3 + 3*n)*m3denJ1010X01100( - 1 + n)*rat(2,1)
          + 1/( - 1 + n)*m3denJ0021X00110( - 1 + n)*rat(-1,1)
          );

   id,ifmatch->sortme m3denJ1001X11110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m3denJ1001X11110( - 2 + n)*rat(n - 2*d + 4,3*n - 3)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X11110( - 1 + n)*rat(d - 2,n - 1)
          + m3denJ1001X11110( - 1 + n)*rat(2*n - 4*d + 7,3*n - 3)
          )

       + delta_( - 2 + n) * (
          + m3denJ1001X11110(0)*rat( - 2*d + 6,3)
          );

   id,ifmatch->sortme m3denJ1001X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m3denJ0001X11110( - 2 + n)*rat( - 3*d + 10,6*n - 6)
          + 1/( - 3 + 3*n)*m3denJ1001X01110( - 2 + n)*n*rat(1,1)
          + 1/( - 3 + 3*n)*m3denJ1001X01110( - 2 + n)*rat(-3,1)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X01110( - 1 + n)*rat(2*d - 4,3*n - 3)
          + m3denJ0001X11110( - 1 + n)*rat( - 5*d + 14,6*n - 6)
          + 1/( - 3 + 3*n)*m3denJ1001X01110( - 1 + n)*n*rat(2,1)
          + 1/( - 3 + 3*n)*m3denJ1001X01110( - 1 + n)*rat(-4,1)
          )

       + delta_( - 2 + n) * (
          + m3denJ0001X11110(0)*rat( - 3*d + 8,6)
          );

   id,ifmatch->sortme m3denJ1000X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m3denJ0001X00110( - 1 + n)*rat(d - 2,6*n - 6)
          + 1/( - 3 + 3*n)*m3denJ1000X01100( - 1 + n)*n*rat(-1,1)
          + 1/( - 3 + 3*n)*m3denJ1000X01100( - 1 + n)*rat(2,1)
          );

   id,ifmatch->sortme m3denJ0021X10110(n?{>1}) =

       + theta_( - 5 + n) * (
          + m3denJ0021X10110( - 4 + n)*rat( - 2*n + 3*d,72*n - 72)
          )

       + theta_( - 4 + n) * (
          + m3denJ0001X01110( - 3 + n)*rat( - d^2 + 4*d - 4,72*n - 72)
          + m3denJ0011X10110( - 3 + n)*rat( - 3*d^2 + 17*d - 24,72*n - 72)
          + m3denJ0021X00110( - 3 + n)*rat( - d + 2,36*n - 36)
          + m3denJ0021X10110( - 3 + n)*rat(4*n + 9*d - 38,72*n - 72)
          )

       + theta_( - 3 + n) * (
          + m3denJ0001X01110( - 2 + n)*rat( - d^2 + 4*d - 4,36*n - 36)
          + m3denJ0011X10110( - 2 + n)*rat(3*d^2 - 17*d + 24,72*n - 72)
          + m3denJ0021X00110( - 2 + n)*rat(d - 2,9*n - 9)
          + m3denJ0021X10110( - 2 + n)*rat(23*n - 18*d - 11,36*n - 36)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X01110( - 1 + n)*rat(d^2 - 4*d + 4,24*n - 24)
          + m3denJ0021X00110( - 1 + n)*rat(d - 2,3*n - 3)
          + m3denJ0021X10110( - 1 + n)*rat(2*n - 3*d + 4,6*n - 6)
          )

       + delta_( - 4 + n) * (
          + m3denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,1080)
          )

       + delta_( - 3 + n) * (
          + m3denJ0011X10110(0)*rat( - 6*d^2 + 31*d - 40,720)
          )

       + delta_( - 2 + n) * (
          + m3denJ0011X10110(0)*rat( - 2*d^2 + 11*d - 15,12)
          + m3denJ0031X10110(0)*rat(5,12)
          );

   id,ifmatch->sortme m3denJ0021X00110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m3denJ0021X00110( - 3 + n)*rat(n - d - 1,36*n - 36)
          )

       + theta_( - 3 + n) * (
          + m3denJ0001X00110( - 2 + n)*rat(d^2 - 4*d + 4,72*n - 72)
          + m3denJ0011X00110( - 2 + n)*rat(3*d^2 - 17*d + 24,72*n - 72)
          + m3denJ0021X00110( - 2 + n)*rat( - n - 6*d + 21,36*n - 36)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X00110( - 1 + n)*rat(d^2 - 4*d + 4,24*n - 24)
          + 1/( - 3 + 3*n)*m3denJ0021X00110( - 1 + n)*n*rat(-2,1)
          + 1/( - 3 + 3*n)*m3denJ0021X00110( - 1 + n)*rat(4,1)
          )

       + delta_( - 3 + n) * (
          + m3denJ0001X01110(0)*rat( - d^2 + 4*d - 4,144)
          )

       + delta_( - 2 + n) * (
          + m3denJ0001X01110(0)*rat( - d^2 + 4*d - 4,24)
          );

   id,ifmatch->sortme m3denJ0011X11110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m3denJ0011X11110( - 2 + n)*rat(2*n - 3*d + 4,6*n - 6)
          )

       + theta_( - 2 + n) * (
          + m3denJ0011X01110( - 1 + n)*rat(d - 2,3*n - 3)
          + m3denJ0011X10110( - 1 + n)*rat( - 3*d + 8,6*n - 6)
          + m3denJ0011X11100( - 1 + n)*rat(d - 2,6*n - 6)
          + m3denJ0011X11110( - 1 + n)*rat(4*n - 5*d + 6,6*n - 6)
          + 1/( - 3 + 3*n)*m3denJ0021X10110( - 1 + n)*rat(5,1)
          )

       + delta_( - 2 + n) * (
          + m3denJ0011X11110(0)*rat( - 3*d + 8,6)
          );

   id,ifmatch->sortme m3denJ0011X11100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m3denJ0011X11100( - 2 + n)*rat(n - d,3*n - 3)
          )

       + theta_( - 2 + n) * (
          + m3denJ0011X01100( - 1 + n)*rat(d - 2,3*n - 3)
          + m3denJ0011X11100( - 1 + n)*rat(2*n - d - 1,3*n - 3)
          )

       + delta_( - 2 + n) * (
          + m3denJ0011X01110(0)*rat( - d + 2,3)
          );

   id,ifmatch->sortme m3denJ0011X10110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m3denJ0011X10110( - 2 + n)*rat(n - 2*d + 3,3*n - 3)
          + 1/( - 1 + n)*m3denJ0021X10110( - 2 + n)*rat(1,1)
          )

       + theta_( - 2 + n) * (
          + m3denJ0011X00110( - 1 + n)*rat(d - 2,3*n - 3)
          + 1/( - 3 + 3*n)*m3denJ0011X10110( - 1 + n)*n*rat(2,1)
          + 1/( - 3 + 3*n)*m3denJ0011X10110( - 1 + n)*rat(-4,1)
          + 1/( - 1 + n)*m3denJ0021X10110( - 1 + n)*rat(-1,1)
          )

       + delta_( - 2 + n) * (
          + m3denJ0011X10110(0)*rat( - 4*d + 10,15)
          );

   id,ifmatch->sortme m3denJ0011X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m3denJ0011X01110( - 2 + n)*rat(n - d,3*n - 3)
          )

       + theta_( - 2 + n) * (
          + m3denJ0011X01100( - 1 + n)*rat(d - 2,6*n - 6)
          + m3denJ0011X01110( - 1 + n)*rat(2*n - d - 1,3*n - 3)
          + m3denJ0011X00110( - 1 + n)*rat( - 3*d + 8,6*n - 6)
          + 1/( - 3 + 3*n)*m3denJ0021X00110( - 1 + n)*rat(5,1)
          )

       + delta_( - 2 + n) * (
          + m3denJ0011X01110(0)*rat( - d + 2,3)
          );

   id,ifmatch->sortme m3denJ0011X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m3denJ0011X01100( - 1 + n)*rat( - 2*n + d + 2,6*n - 6)
          );

   id,ifmatch->sortme m3denJ0011X00110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m3denJ0011X00110( - 1 + n)*rat( - 2*n + 3*d - 4,6*n - 6)
          + 1/( - 1 + n)*m3denJ0021X00110( - 1 + n)*rat(-1,1)
          );

   id,ifmatch->sortme m3denJ0001X11110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m3denJ0001X11110( - 2 + n)*rat(2*n - 3*d + 4,6*n - 6)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X01110( - 1 + n)*rat(2*d - 4,3*n - 3)
          + m3denJ0001X11110( - 1 + n)*rat(4*n - 5*d + 6,6*n - 6)
          )

       + delta_( - 2 + n) * (
          + m3denJ0001X11110(0)*rat( - 3*d + 8,6)
          );

   id,ifmatch->sortme m3denJ0001X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m3denJ0001X01110( - 2 + n)*rat(n - d,3*n - 3)
          )

       + theta_( - 2 + n) * (
          + m3denJ0001X01110( - 1 + n)*rat(2*n - d - 1,3*n - 3)
          + m3denJ0001X00110( - 1 + n)*rat(d - 2,3*n - 3)
          )

       + delta_( - 2 + n) * (
          + m3denJ0001X01110(0)*rat( - d + 2,3)
          );

   id,ifmatch->sortme m3denJ0001X00110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m3denJ0001X00110( - 1 + n)*rat( - 2*n + d + 2,6*n - 6)
          );

********************************************************************************                 

goto endrec;

la sortme;
$irep = 0;
la endrec;

ModuleOption,minimum,$irep;

.sort:rec1d-m3-`$repcount++';

#redefine irep "`$irep'"
#enddo
#endprocedure
*--#] rec1dm3 :

*--#[ rec1dm4 :
#procedure rec1dm4
#$repcount = 1;

#do irep = 1,1
#$irep = 1;
                
********************************************************************************                 
   id,ifmatch->sortme m4denJ2010X11100(n?{>1}) =

       + theta_( - 4 + n) * (
          + m4denJ0021X00110( - 3 + n)*rat( - n*d + 2*n + d^2 - d - 2,240*n^2
          - 120*n*d - 120*n + 120*d - 120)
          + m4denJ0021X10110( - 3 + n)*rat( - 2*n*d + 7*n - d + 2,120*n^2 - 60
         *n*d - 60*n + 60*d - 60)
          + 1/( - 60 + 60*n)*m4denJ2010X11100( - 3 + n)*n*rat(1,1)
          + 1/( - 60 + 60*n)*m4denJ2010X11100( - 3 + n)*rat(-4,1)
          )

       + theta_( - 3 + n) * (
          + m4denJ0001X01110( - 2 + n)*rat( - d^2 + 4*d - 4, - 120*n + 60*d - 
         60)
          + m4denJ0001X00110( - 2 + n)*rat(d^3 - 6*d^2 + 12*d - 8, - 480*n^2
          + 240*n*d + 240*n - 240*d + 240)
          + m4denJ0011X10110( - 2 + n)*rat( - 3*d^2 + 17*d - 24, - 120*n + 60*
         d - 60)
          + m4denJ0011X00110( - 2 + n)*rat(3*d^3 - 23*d^2 + 58*d - 48, - 480*
         n^2 + 240*n*d + 240*n - 240*d + 240)
          + m4denJ0021X00110( - 2 + n)*rat(n*d - 2*n + 4*d^2 - 19*d + 22,120*
         n^2 - 60*n*d - 60*n + 60*d - 60)
          + m4denJ0021X10110( - 2 + n)*rat( - 8*n*d + 23*n + 6*d - 17,60*n^2
          - 30*n*d - 30*n + 30*d - 30)
          + 1/( - 30 + 30*n)*m4denJ2010X11100( - 2 + n)*n*rat(1,1)
          + 1/( - 30 + 30*n)*m4denJ2010X11100( - 2 + n)*rat(-3,1)
          )

       + theta_( - 2 + n) * (
          + m4denJ0001X01110( - 1 + n)*rat( - d^2 + 4*d - 4, - 30*n + 15*d - 
         15)
          + m4denJ0001X00110( - 1 + n)*rat(d^3 - 6*d^2 + 12*d - 8, - 120*n^2
          + 60*n*d + 60*n - 60*d + 60)
          + m4denJ0011X10110( - 1 + n)*rat( - 3*d^2 + 17*d - 24, - 120*n + 60*
         d - 60)
          + m4denJ0011X00110( - 1 + n)*rat(3*d^3 - 23*d^2 + 58*d - 48, - 480*
         n^2 + 240*n*d + 240*n - 240*d + 240)
          + m4denJ0021X00110( - 1 + n)*rat(15*n*d - 30*n + 7*d^2 - 73*d + 118,
         240*n^2 - 120*n*d - 120*n + 120*d - 120)
          + m4denJ0021X10110( - 1 + n)*rat( - 14*n*d + 27*n + 37*d - 96,120*
         n^2 - 60*n*d - 60*n + 60*d - 60)
          + 1/( - 60 + 60*n)*m4denJ2010X11100( - 1 + n)*n*rat(-23,1)
          + 1/( - 60 + 60*n)*m4denJ2010X11100( - 1 + n)*rat(46,1)
          )

       + delta_( - 3 + n) * (
          + m4denJ0001X01110(0)*rat( - d^3 + 6*d^2 - 12*d + 8,480*d - 3360)
          + m4denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,300*d - 2100)
          )

       + delta_( - 2 + n) * (
          + m4denJ0001X01110(0)*rat( - d^3 + 6*d^2 - 12*d + 8,60*d - 300)
          + m4denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,300*d - 1500)
          );

   id,ifmatch->sortme m4denJ2010X01100(n?{>1}) =

       + theta_( - 4 + n) * (
          + m4denJ0021X00110( - 3 + n)*rat( - d + 3,60*n - 60)
          + 1/( - 60 + 60*n)*m4denJ2010X01100( - 3 + n)*n*rat(1,1)
          + 1/( - 60 + 60*n)*m4denJ2010X01100( - 3 + n)*rat(-4,1)
          )

       + theta_( - 3 + n) * (
          + m4denJ0001X00110( - 2 + n)*rat(d^2 - 4*d + 4,120*n - 120)
          + m4denJ0011X00110( - 2 + n)*rat(3*d^2 - 17*d + 24,120*n - 120)
          + m4denJ0021X00110( - 2 + n)*rat( - 2*d + 6,15*n - 15)
          + 1/( - 30 + 30*n)*m4denJ2010X01100( - 2 + n)*n*rat(1,1)
          + 1/( - 30 + 30*n)*m4denJ2010X01100( - 2 + n)*rat(-3,1)
          )

       + theta_( - 2 + n) * (
          + m4denJ0001X00110( - 1 + n)*rat(d^2 - 4*d + 4,30*n - 30)
          + m4denJ0011X00110( - 1 + n)*rat(3*d^2 - 17*d + 24,120*n - 120)
          + m4denJ0021X00110( - 1 + n)*rat( - 7*d + 21,60*n - 60)
          + 1/( - 60 + 60*n)*m4denJ2010X01100( - 1 + n)*n*rat(-23,1)
          + 1/( - 60 + 60*n)*m4denJ2010X01100( - 1 + n)*rat(46,1)
          )

       + delta_( - 3 + n) * (
          + m4denJ0001X01110(0)*rat( - d^2 + 4*d - 4,240)
          )

       + delta_( - 2 + n) * (
          + m4denJ0001X01110(0)*rat( - d^2 + 4*d - 4,30)
          );

   id,ifmatch->sortme m4denJ1011X11110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m4denJ0021X00110( - 3 + n)*rat(n*d - 2*n - d^2 + d + 2,48*n^3 - 72
         *n^2*d + 72*n^2 + 24*n*d^2 - 72*n - 24*d^2 + 72*d - 48)
          + m4denJ0021X10110( - 3 + n)*rat( - 2*n + 3*d - 2,48*n^2 - 72*n*d + 
         120*n + 24*d^2 - 72*d + 48)
          )

       + theta_( - 3 + n) * (
          + m4denJ0001X01110( - 2 + n)*rat( - d^2 + 4*d - 4,48*n^2 - 72*n*d + 
         120*n + 24*d^2 - 72*d + 48)
          + m4denJ0001X00110( - 2 + n)*rat(d^3 - 6*d^2 + 12*d - 8,96*n^3 - 144
         *n^2*d + 144*n^2 + 48*n*d^2 - 144*n - 48*d^2 + 144*d - 96)
          + m4denJ0011X10110( - 2 + n)*rat( - 3*d^2 + 17*d - 24,48*n^2 - 72*n*
         d + 120*n + 24*d^2 - 72*d + 48)
          + m4denJ0011X00110( - 2 + n)*rat(3*d^3 - 23*d^2 + 58*d - 48,96*n^3
          - 144*n^2*d + 144*n^2 + 48*n*d^2 - 144*n - 48*d^2 + 144*d - 96)
          + m4denJ0021X00110( - 2 + n)*rat( - 2*d^2 + 9*d - 10,12*n^3 - 18*n^2
         *d + 18*n^2 + 6*n*d^2 - 18*n - 6*d^2 + 18*d - 12)
          + m4denJ0021X10110( - 2 + n)*rat( - 2*n + 9*d - 20,24*n^2 - 36*n*d
          + 60*n + 12*d^2 - 36*d + 24)
          + m4denJ1011X01110( - 2 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,16*n^2 - 
         8*n*d - 8*n + 8*d - 8)
          + m4denJ1011X11110( - 2 + n)*rat( - 2*n + 3*d - 6,8*n - 4*d + 4)
          )

       + theta_( - 2 + n) * (
          + m4denJ0001X01110( - 1 + n)*rat( - 5*n^2*d^2 + 20*n^2*d - 20*n^2 + 
         18*n*d^3 - 125*n*d^2 + 284*n*d - 212*n - 9*d^4 + 72*d^3 - 203*d^2 + 
         236*d - 92,96*n^4 - 288*n^3*d + 480*n^3 + 264*n^2*d^2 - 720*n^2*d + 
         360*n^2 - 72*n*d^3 + 120*n*d^2 + 360*n*d - 600*n + 72*d^3 - 384*d^2
          + 648*d - 336)
          + m4denJ0001X11110( - 1 + n)*rat(8*n^3*d - 16*n^3 - 28*n^2*d^2 + 120
         *n^2*d - 128*n^2 + 26*n*d^3 - 150*n*d^2 + 272*n*d - 152*n - 3*d^4 + 4
         *d^3 + 67*d^2 - 220*d + 188,64*n^4 - 192*n^3*d + 320*n^3 + 176*n^2*
         d^2 - 480*n^2*d + 240*n^2 - 48*n*d^3 + 80*n*d^2 + 240*n*d - 400*n + 
         48*d^3 - 256*d^2 + 432*d - 224)
          + m4denJ0001X00110( - 1 + n)*rat(14*n*d^3 - 84*n*d^2 + 168*n*d - 112
         *n - 39*d^4 + 364*d^3 - 1248*d^2 + 1872*d - 1040,768*n^4 - 2304*n^3*d
          + 3840*n^3 + 2112*n^2*d^2 - 5760*n^2*d + 2880*n^2 - 576*n*d^3 + 960*
         n*d^2 + 2880*n*d - 4800*n + 576*d^3 - 3072*d^2 + 5184*d - 2688)
          + m4denJ0011X01100( - 1 + n)*rat( - 2*n*d^2 + 8*n*d - 8*n + d^3 - 2*
         d^2 - 4*d + 8,64*n^3 - 96*n^2*d + 96*n^2 + 32*n*d^2 - 96*n - 32*d^2
          + 96*d - 64)
          + m4denJ0011X01110( - 1 + n)*rat(n*d^2 - 5*n*d + 6*n - d^3 + 6*d^2
          - 11*d + 6,8*n^3 - 12*n^2*d + 12*n^2 + 4*n*d^2 - 12*n - 4*d^2 + 12*d
          - 8)
          + m4denJ0011X10110( - 1 + n)*rat( - 9*n*d + 24*n + 15*d^2 - 67*d + 
         72,48*n^2 - 72*n*d + 120*n + 24*d^2 - 72*d + 48)
          + m4denJ0011X11100( - 1 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,16*n^2 - 
         24*n*d + 40*n + 8*d^2 - 24*d + 16)
          + m4denJ0011X11110( - 1 + n)*rat( - 2*n*d + 6*n + 3*d^2 - 15*d + 18,
         8*n^2 - 12*n*d + 20*n + 4*d^2 - 12*d + 8)
          + m4denJ0011X00110( - 1 + n)*rat(18*n*d^2 - 84*n*d + 96*n - 21*d^3
          + 116*d^2 - 196*d + 96,192*n^3 - 288*n^2*d + 288*n^2 + 96*n*d^2 - 
         288*n - 96*d^2 + 288*d - 192)
          + m4denJ0021X00110( - 1 + n)*rat( - 26*n*d + 52*n + d^2 + 76*d - 156
         ,96*n^3 - 144*n^2*d + 144*n^2 + 48*n*d^2 - 144*n - 48*d^2 + 144*d - 
         96)
          + m4denJ0021X10110( - 1 + n)*rat(11*n - 6*d - 13,24*n^2 - 36*n*d + 
         60*n + 12*d^2 - 36*d + 24)
          + m4denJ1001X11110( - 1 + n)*rat( - d^2 + 6*d - 9,8*n^2 - 16*n*d + 
         32*n + 6*d^2 - 20*d + 14)
          + m4denJ1011X01110( - 1 + n)*rat(2*n*d - 4*n - 5*d^2 + 18*d - 16,32*
         n^2 - 16*n*d - 16*n + 16*d - 16)
          + m4denJ1011X11110( - 1 + n)*rat( - 5*n + 5*d - 8,4*n - 2*d + 2)
          )

       + delta_( - 3 + n) * (
          + m4denJ0001X01110(0)*rat( - d^3 + 6*d^2 - 12*d + 8,96*d^2 - 1152*d
          + 3360)
          + m4denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,120*d^2 - 1440*d + 4200)
          )

       + delta_( - 2 + n) * (
          + m4denJ0001X01110(0)*rat( - d^3 + 6*d^2 - 12*d + 8,12*d^2 - 108*d
          + 240)
          + m4denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,120*d^2 - 1080*d + 2400)
          + m4denJ1011X01110(0)*rat(d^2 - 5*d + 6,8*d - 40)
          + m4denJ1011X11110(0)*rat( - 3*d + 10,4*d - 20)
          );

   id,ifmatch->sortme m4denJ1011X01110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m4denJ0021X00110( - 3 + n)*rat( - n + d + 1,24*n^2 - 12*n*d - 12*n
          + 12*d - 12)
          )

       + theta_( - 3 + n) * (
          + m4denJ0001X00110( - 2 + n)*rat(d^2 - 4*d + 4, - 48*n^2 + 24*n*d + 
         24*n - 24*d + 24)
          + m4denJ0011X00110( - 2 + n)*rat(3*d^2 - 17*d + 24, - 48*n^2 + 24*n*
         d + 24*n - 24*d + 24)
          + m4denJ0021X00110( - 2 + n)*rat( - n + 4*d - 9,12*n^2 - 6*n*d - 6*n
          + 6*d - 6)
          + m4denJ1011X01110( - 2 + n)*rat( - n + d - 1,4*n - 4)
          )

       + theta_( - 2 + n) * (
          + m4denJ0001X01110( - 1 + n)*rat(n^2*d - 2*n^2 - 2*n*d^2 + 9*n*d - 
         10*n + d^3 - 6*d^2 + 11*d - 6,8*n^3 - 12*n^2*d + 12*n^2 + 4*n*d^2 - 
         12*n - 4*d^2 + 12*d - 8)
          + m4denJ0001X11110( - 1 + n)*rat(d^2 - 5*d + 6, - 8*n^2 + 8*n*d - 8*
         n - 8*d + 16)
          + m4denJ0001X00110( - 1 + n)*rat( - 10*n*d^2 + 40*n*d - 40*n + 13*
         d^3 - 90*d^2 + 204*d - 152,192*n^3 - 288*n^2*d + 288*n^2 + 96*n*d^2
          - 288*n - 96*d^2 + 288*d - 192)
          + m4denJ0011X01100( - 1 + n)*rat(2*n*d - 4*n - d^2 + 4,32*n^2 - 16*n
         *d - 16*n + 16*d - 16)
          + m4denJ0011X01110( - 1 + n)*rat( - n*d + 3*n + d^2 - 4*d + 3,4*n^2
          - 2*n*d - 2*n + 2*d - 2)
          + m4denJ0011X00110( - 1 + n)*rat( - 18*n*d + 48*n + 21*d^2 - 74*d + 
         48,96*n^2 - 48*n*d - 48*n + 48*d - 48)
          + m4denJ0021X00110( - 1 + n)*rat(22*n - d - 74,48*n^2 - 24*n*d - 24*
         n + 24*d - 24)
          + m4denJ1011X01110( - 1 + n)*rat( - 10*n + 5*d,8*n - 8)
          )

       + delta_( - 3 + n) * (
          + m4denJ0001X01110(0)*rat( - d^2 + 4*d - 4,48*d - 336)
          )

       + delta_( - 2 + n) * (
          + m4denJ0001X01110(0)*rat( - d^2 + 4*d - 4,6*d - 30)
          + m4denJ1011X01110(0)*rat(d - 3,4)
          );

   id,ifmatch->sortme m4denJ1010X11100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0011X10110( - 1 + n)*rat(3*n*d - 9*n - 2*d + 6,8*n^2 - 4*n*d
          - 4*n + 4*d - 4)
          + m4denJ0011X00110( - 1 + n)*rat(2*n*d - 4*n - 3*d^2 + 10*d - 8,32*
         n^2 - 16*n*d - 16*n + 16*d - 16)
          + m4denJ0021X00110( - 1 + n)*rat(3*d - 6,16*n^2 - 8*n*d - 8*n + 8*d
          - 8)
          + m4denJ0021X10110( - 1 + n)*rat(-3,4*n - 2*d + 2)
          + 1/( - 4 + 4*n)*m4denJ1010X11100( - 1 + n)*n*rat(-1,1)
          + 1/( - 4 + 4*n)*m4denJ1010X11100( - 1 + n)*rat(2,1)
          );

   id,ifmatch->sortme m4denJ1010X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0011X00110( - 1 + n)*rat(3*d - 8,8*n - 8)
          + 1/( - 4 + 4*n)*m4denJ0021X00110( - 1 + n)*rat(-3,1)
          + 1/( - 4 + 4*n)*m4denJ1010X01100( - 1 + n)*n*rat(-1,1)
          + 1/( - 4 + 4*n)*m4denJ1010X01100( - 1 + n)*rat(2,1)
          );

   id,ifmatch->sortme m4denJ1001X11110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0001X01110( - 1 + n)*rat( - 3*n*d^2 + 12*n*d - 12*n + 3*d^3
          - 15*d^2 + 24*d - 12,32*n^3 - 96*n^2*d + 192*n^2 + 88*n*d^2 - 336*n*
         d + 312*n - 24*d^3 + 128*d^2 - 216*d + 112)
          + m4denJ0001X11110( - 1 + n)*rat(6*n*d - 12*n - 9*d^2 + 36*d - 36,32
         *n^2 - 80*n*d + 176*n + 48*d^2 - 208*d + 224)
          + m4denJ0001X00110( - 1 + n)*rat(6*n*d^3 - 36*n*d^2 + 72*n*d - 48*n
          - 3*d^4 + 12*d^3 - 48*d + 48,256*n^4 - 768*n^3*d + 1280*n^3 + 704*
         n^2*d^2 - 1920*n^2*d + 960*n^2 - 192*n*d^3 + 320*n*d^2 + 960*n*d - 
         1600*n + 192*d^3 - 1024*d^2 + 1728*d - 896)
          + m4denJ1001X11110( - 1 + n)*rat( - n + 2*d - 5,4*n - 6*d + 14)
          );

   id,ifmatch->sortme m4denJ1001X01110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0001X01110( - 1 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,8*n^2 - 
         12*n*d + 20*n + 4*d^2 - 12*d + 8)
          + m4denJ0001X11110( - 1 + n)*rat(n*d - 4*n + d - 2,8*n^2 - 8*n*d + 8
         *n + 8*d - 16)
          + m4denJ0001X00110( - 1 + n)*rat( - 2*n*d^2 + 8*n*d - 8*n + d^3 - 2*
         d^2 - 4*d + 8,64*n^3 - 96*n^2*d + 96*n^2 + 32*n*d^2 - 96*n - 32*d^2
          + 96*d - 64)
          + 1/( - 4 + 4*n)*m4denJ1001X01110( - 1 + n)*n*rat(-1,1)
          + 1/( - 4 + 4*n)*m4denJ1001X01110( - 1 + n)*rat(2,1)
          );

   id,ifmatch->sortme m4denJ1000X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0001X00110( - 1 + n)*rat(d - 2,8*n - 8)
          + 1/( - 4 + 4*n)*m4denJ1000X01100( - 1 + n)*n*rat(-1,1)
          + 1/( - 4 + 4*n)*m4denJ1000X01100( - 1 + n)*rat(2,1)
          );

   id,ifmatch->sortme m4denJ0021X10110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m4denJ0021X00110( - 3 + n)*rat( - n*d + 2*n + d^2 - d - 2,240*n^2
          - 120*n*d - 120*n + 120*d - 120)
          + m4denJ0021X10110( - 3 + n)*rat(2*n - 3*d + 2,120*n - 60*d + 60)
          )

       + theta_( - 3 + n) * (
          + m4denJ0001X01110( - 2 + n)*rat( - d^2 + 4*d - 4, - 120*n + 60*d - 
         60)
          + m4denJ0001X00110( - 2 + n)*rat(d^3 - 6*d^2 + 12*d - 8, - 480*n^2
          + 240*n*d + 240*n - 240*d + 240)
          + m4denJ0011X10110( - 2 + n)*rat( - 3*d^2 + 17*d - 24, - 120*n + 60*
         d - 60)
          + m4denJ0011X00110( - 2 + n)*rat(3*d^3 - 23*d^2 + 58*d - 48, - 480*
         n^2 + 240*n*d + 240*n - 240*d + 240)
          + m4denJ0021X00110( - 2 + n)*rat(n*d - 2*n + 4*d^2 - 19*d + 22,120*
         n^2 - 60*n*d - 60*n + 60*d - 60)
          + m4denJ0021X10110( - 2 + n)*rat(2*n - 9*d + 20,60*n - 30*d + 30)
          )

       + theta_( - 2 + n) * (
          + m4denJ0001X01110( - 1 + n)*rat( - d^2 + 4*d - 4, - 30*n + 15*d - 
         15)
          + m4denJ0001X00110( - 1 + n)*rat(d^3 - 6*d^2 + 12*d - 8, - 120*n^2
          + 60*n*d + 60*n - 60*d + 60)
          + m4denJ0011X10110( - 1 + n)*rat( - 3*d^2 + 17*d - 24, - 120*n + 60*
         d - 60)
          + m4denJ0011X00110( - 1 + n)*rat(3*d^3 - 23*d^2 + 58*d - 48, - 480*
         n^2 + 240*n*d + 240*n - 240*d + 240)
          + m4denJ0021X00110( - 1 + n)*rat(15*n*d - 30*n + 7*d^2 - 73*d + 118,
         240*n^2 - 120*n*d - 120*n + 120*d - 120)
          + m4denJ0021X10110( - 1 + n)*rat( - 46*n + 9*d + 50,120*n - 60*d + 
         60)
          )

       + delta_( - 3 + n) * (
          + m4denJ0001X01110(0)*rat( - d^3 + 6*d^2 - 12*d + 8,480*d - 3360)
          + m4denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,300*d - 2100)
          )

       + delta_( - 2 + n) * (
          + m4denJ0001X01110(0)*rat( - d^3 + 6*d^2 - 12*d + 8,60*d - 300)
          + m4denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,300*d - 1500)
          );

   id,ifmatch->sortme m4denJ0021X00110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m4denJ0021X00110( - 3 + n)*rat(n - d - 1,60*n - 60)
          )

       + theta_( - 3 + n) * (
          + m4denJ0001X00110( - 2 + n)*rat(d^2 - 4*d + 4,120*n - 120)
          + m4denJ0011X00110( - 2 + n)*rat(3*d^2 - 17*d + 24,120*n - 120)
          + m4denJ0021X00110( - 2 + n)*rat(n - 4*d + 9,30*n - 30)
          )

       + theta_( - 2 + n) * (
          + m4denJ0001X00110( - 1 + n)*rat(d^2 - 4*d + 4,30*n - 30)
          + m4denJ0011X00110( - 1 + n)*rat(3*d^2 - 17*d + 24,120*n - 120)
          + m4denJ0021X00110( - 1 + n)*rat( - 23*n - 7*d + 67,60*n - 60)
          )

       + delta_( - 3 + n) * (
          + m4denJ0001X01110(0)*rat( - d^2 + 4*d - 4,240)
          )

       + delta_( - 2 + n) * (
          + m4denJ0001X01110(0)*rat( - d^2 + 4*d - 4,30)
          );

   id,ifmatch->sortme m4denJ0011X11110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m4denJ0021X00110( - 3 + n)*rat(n*d - 2*n - d^2 + d + 2,96*n^3 - 
         144*n^2*d + 144*n^2 + 48*n*d^2 - 144*n - 48*d^2 + 144*d - 96)
          + m4denJ0021X10110( - 3 + n)*rat( - 2*n + 3*d - 2,96*n^2 - 144*n*d
          + 240*n + 48*d^2 - 144*d + 96)
          )

       + theta_( - 3 + n) * (
          + m4denJ0001X01110( - 2 + n)*rat( - d^2 + 4*d - 4,96*n^2 - 144*n*d
          + 240*n + 48*d^2 - 144*d + 96)
          + m4denJ0001X00110( - 2 + n)*rat(d^3 - 6*d^2 + 12*d - 8,192*n^3 - 
         288*n^2*d + 288*n^2 + 96*n*d^2 - 288*n - 96*d^2 + 288*d - 192)
          + m4denJ0011X10110( - 2 + n)*rat( - 3*d^2 + 17*d - 24,96*n^2 - 144*n
         *d + 240*n + 48*d^2 - 144*d + 96)
          + m4denJ0011X00110( - 2 + n)*rat(3*d^3 - 23*d^2 + 58*d - 48,192*n^3
          - 288*n^2*d + 288*n^2 + 96*n*d^2 - 288*n - 96*d^2 + 288*d - 192)
          + m4denJ0021X00110( - 2 + n)*rat( - 2*d^2 + 9*d - 10,24*n^3 - 36*n^2
         *d + 36*n^2 + 12*n*d^2 - 36*n - 12*d^2 + 36*d - 24)
          + m4denJ0021X10110( - 2 + n)*rat( - 2*n + 9*d - 20,48*n^2 - 72*n*d
          + 120*n + 24*d^2 - 72*d + 48)
          )

       + theta_( - 2 + n) * (
          + m4denJ0001X01110( - 1 + n)*rat( - d^2 + 4*d - 4,24*n^2 - 36*n*d + 
         60*n + 12*d^2 - 36*d + 24)
          + m4denJ0001X00110( - 1 + n)*rat(d^3 - 6*d^2 + 12*d - 8,48*n^3 - 72*
         n^2*d + 72*n^2 + 24*n*d^2 - 72*n - 24*d^2 + 72*d - 48)
          + m4denJ0011X01100( - 1 + n)*rat( - 2*n*d^2 + 8*n*d - 8*n + d^3 - 2*
         d^2 - 4*d + 8,128*n^3 - 192*n^2*d + 192*n^2 + 64*n*d^2 - 192*n - 64*
         d^2 + 192*d - 128)
          + m4denJ0011X01110( - 1 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,16*n^2 - 
         24*n*d + 40*n + 8*d^2 - 24*d + 16)
          + m4denJ0011X10110( - 1 + n)*rat( - 9*n*d + 24*n + 15*d^2 - 67*d + 
         72,96*n^2 - 144*n*d + 240*n + 48*d^2 - 144*d + 96)
          + m4denJ0011X11100( - 1 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,32*n^2 - 
         48*n*d + 80*n + 16*d^2 - 48*d + 32)
          + m4denJ0011X11110( - 1 + n)*rat( - 2*n + 3*d - 6,8*n - 8*d + 16)
          + m4denJ0011X00110( - 1 + n)*rat(18*n*d^2 - 84*n*d + 96*n - 21*d^3
          + 116*d^2 - 196*d + 96,384*n^3 - 576*n^2*d + 576*n^2 + 192*n*d^2 - 
         576*n - 192*d^2 + 576*d - 384)
          + m4denJ0021X00110( - 1 + n)*rat( - 38*n*d + 76*n + 13*d^2 + 28*d - 
         108,192*n^3 - 288*n^2*d + 288*n^2 + 96*n*d^2 - 288*n - 96*d^2 + 288*d
          - 192)
          + m4denJ0021X10110( - 1 + n)*rat(23*n - 18*d + 11,48*n^2 - 72*n*d + 
         120*n + 24*d^2 - 72*d + 48)
          )

       + delta_( - 3 + n) * (
          + m4denJ0001X01110(0)*rat( - d^3 + 6*d^2 - 12*d + 8,192*d^2 - 2304*d
          + 6720)
          + m4denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,240*d^2 - 2880*d + 8400)
          )

       + delta_( - 2 + n) * (
          + m4denJ0001X01110(0)*rat( - d^3 + 6*d^2 - 12*d + 8,24*d^2 - 216*d
          + 480)
          + m4denJ0011X10110(0)*rat(6*d^2 - 31*d + 40,240*d^2 - 2160*d + 4800)
          );

   id,ifmatch->sortme m4denJ0011X11100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0011X01100( - 1 + n)*rat(2*n*d - 4*n - d^2 + 4,32*n^2 - 16*n
         *d - 16*n + 16*d - 16)
          + m4denJ0011X11100( - 1 + n)*rat( - n + d - 1,4*n - 2*d + 2)
          );

   id,ifmatch->sortme m4denJ0011X10110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0011X10110( - 1 + n)*rat( - n + 2*d - 4,4*n - 2*d + 2)
          + m4denJ0011X00110( - 1 + n)*rat(2*n*d - 4*n - 3*d^2 + 10*d - 8,32*
         n^2 - 16*n*d - 16*n + 16*d - 16)
          + m4denJ0021X00110( - 1 + n)*rat(3*d - 6,16*n^2 - 8*n*d - 8*n + 8*d
          - 8)
          + m4denJ0021X10110( - 1 + n)*rat(-3,4*n - 2*d + 2)
          );

   id,ifmatch->sortme m4denJ0011X01110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m4denJ0021X00110( - 3 + n)*rat( - n + d + 1,48*n^2 - 24*n*d - 24*n
          + 24*d - 24)
          )

       + theta_( - 3 + n) * (
          + m4denJ0001X00110( - 2 + n)*rat(d^2 - 4*d + 4, - 96*n^2 + 48*n*d + 
         48*n - 48*d + 48)
          + m4denJ0011X00110( - 2 + n)*rat(3*d^2 - 17*d + 24, - 96*n^2 + 48*n*
         d + 48*n - 48*d + 48)
          + m4denJ0021X00110( - 2 + n)*rat( - n + 4*d - 9,24*n^2 - 12*n*d - 12
         *n + 12*d - 12)
          )

       + theta_( - 2 + n) * (
          + m4denJ0001X00110( - 1 + n)*rat(d^2 - 4*d + 4, - 24*n^2 + 12*n*d + 
         12*n - 12*d + 12)
          + m4denJ0011X01100( - 1 + n)*rat(2*n*d - 4*n - d^2 + 4,64*n^2 - 32*n
         *d - 32*n + 32*d - 32)
          + m4denJ0011X01110( - 1 + n)*rat( - n + d - 1,4*n - 2*d + 2)
          + m4denJ0011X00110( - 1 + n)*rat( - 18*n*d + 48*n + 21*d^2 - 74*d + 
         48,192*n^2 - 96*n*d - 96*n + 96*d - 96)
          + m4denJ0021X00110( - 1 + n)*rat(46*n - 13*d - 62,96*n^2 - 48*n*d - 
         48*n + 48*d - 48)
          )

       + delta_( - 3 + n) * (
          + m4denJ0001X01110(0)*rat( - d^2 + 4*d - 4,96*d - 672)
          )

       + delta_( - 2 + n) * (
          + m4denJ0001X01110(0)*rat( - d^2 + 4*d - 4,12*d - 60)
          );

   id,ifmatch->sortme m4denJ0011X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0011X01100( - 1 + n)*rat( - 2*n + d + 2,8*n - 8)
          );

   id,ifmatch->sortme m4denJ0011X00110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0011X00110( - 1 + n)*rat( - 2*n + 3*d - 4,8*n - 8)
          + 1/( - 4 + 4*n)*m4denJ0021X00110( - 1 + n)*rat(-3,1)
          );

   id,ifmatch->sortme m4denJ0001X11110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0001X01110( - 1 + n)*rat(n*d - 2*n - d^2 + 3*d - 2,8*n^2 - 
         12*n*d + 20*n + 4*d^2 - 12*d + 8)
          + m4denJ0001X11110( - 1 + n)*rat( - 2*n + 3*d - 6,8*n - 8*d + 16)
          + m4denJ0001X00110( - 1 + n)*rat( - 2*n*d^2 + 8*n*d - 8*n + d^3 - 2*
         d^2 - 4*d + 8,64*n^3 - 96*n^2*d + 96*n^2 + 32*n*d^2 - 96*n - 32*d^2
          + 96*d - 64)
          );

   id,ifmatch->sortme m4denJ0001X01110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0001X01110( - 1 + n)*rat( - n + d - 1,4*n - 2*d + 2)
          + m4denJ0001X00110( - 1 + n)*rat(2*n*d - 4*n - d^2 + 4,32*n^2 - 16*n
         *d - 16*n + 16*d - 16)
          );

   id,ifmatch->sortme m4denJ0001X00110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m4denJ0001X00110( - 1 + n)*rat( - 2*n + d + 2,8*n - 8)
          );

********************************************************************************                 

goto endrec;

la sortme;
$irep = 0;
la endrec;

ModuleOption,minimum,$irep;

.sort:rec1d-m4-`$repcount++';

#redefine irep "`$irep'"
#enddo
#endprocedure
*--#] rec1dm4 :

*--#[ rec1dm9 :
#procedure rec1dm9
#$repcount = 1;

#do irep = 1,1
#$irep = 1;
                
********************************************************************************             

   id,ifmatch->sortme m9denJ2010X11100(n?{>1}) =

       + theta_( - 4 + n) * (
          + m9denJ0021X10110( - 3 + n)*rat(n*d - 4*n + 5*d - 14,720*n^2 - 720*
         n*d + 720*n + 720*d - 1440)
          + 1/( - 360 + 360*n)*m9denJ2010X11100( - 3 + n)*n*rat(-1,1)
          + 1/( - 360 + 360*n)*m9denJ2010X11100( - 3 + n)*rat(4,1)
          )

       + theta_( - 3 + n) * (
          + m9denJ0001X01110( - 2 + n)*rat(d^3 - 5*d^2 + 8*d - 4, - 720*n^2 + 
         720*n*d - 720*n - 720*d + 1440)
          + m9denJ0011X10110( - 2 + n)*rat(3*n*d^2 - 17*n*d + 24*n + 12*d^3 - 
         95*d^2 + 249*d - 216, - 2160*n^2 + 2160*n*d - 2160*n - 2160*d + 4320)
          + m9denJ0021X00110( - 2 + n)*rat(d - 2, - 360*n + 360*d - 720)
          + m9denJ0021X10110( - 2 + n)*rat(19*n*d - 70*n + 6*d^2 + 35*d - 146,
         720*n^2 - 720*n*d + 720*n + 720*d - 1440)
          + 1/( - 180 + 180*n)*m9denJ2010X11100( - 2 + n)*n*rat(-11,1)
          + 1/( - 180 + 180*n)*m9denJ2010X11100( - 2 + n)*rat(33,1)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X01110( - 1 + n)*rat(7*d^3 - 35*d^2 + 56*d - 28, - 720*
         n^2 + 720*n*d - 720*n - 720*d + 1440)
          + m9denJ0001X00110( - 1 + n)*rat( - d^3 + 6*d^2 - 12*d + 8, - 720*
         n^2 + 720*n*d - 720*n - 720*d + 1440)
          + m9denJ0011X10110( - 1 + n)*rat(15*n*d^2 - 85*n*d + 120*n + 72*d^3
          - 555*d^2 + 1409*d - 1176, - 2160*n^2 + 2160*n*d - 2160*n - 2160*d
          + 4320)
          + m9denJ0011X00110( - 1 + n)*rat( - 3*d^3 + 23*d^2 - 58*d + 48, - 
         1080*n^2 + 1080*n*d - 1080*n - 1080*d + 2160)
          + m9denJ0021X00110( - 1 + n)*rat(d - 2, - 45*n + 45*d - 90)
          + m9denJ0021X10110( - 1 + n)*rat(41*n*d - 143*n + 15*d^2 + 31*d - 
         208,360*n^2 - 360*n*d + 360*n + 360*d - 720)
          + 1/( - 360 + 360*n)*m9denJ2010X11100( - 1 + n)*n*rat(-157,1)
          + 1/( - 360 + 360*n)*m9denJ2010X11100( - 1 + n)*rat(314,1)
          )

       + delta_( - 3 + n) * (
          + m9denJ0011X10110(0)*rat( - 6*d^2 + 31*d - 40,3600*d - 18000)
          )

       + delta_( - 2 + n) * (
          + m9denJ0001X01110(0)*rat(d^3 - 6*d^2 + 12*d - 8,720*d - 2880)
          + m9denJ0011X10110(0)*rat(24*d^3 - 394*d^2 + 1555*d - 1800,10800*d
          - 43200)
          );

   id,ifmatch->sortme m9denJ2010X01100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0021X00110( - 2 + n)*rat(d - 3,36*n^2 - 36*n*d + 36*n + 36*d
          - 72)
          + 1/( - 72 + 72*n)*m9denJ2010X01100( - 2 + n)*n*rat(-1,1)
          + 1/( - 72 + 72*n)*m9denJ2010X01100( - 2 + n)*rat(3,1)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X00110( - 1 + n)*rat(d^3 - 4*d^2 + 4*d, - 288*n^2 + 288*
         n*d - 288*n - 288*d + 576)
          + m9denJ0011X00110( - 1 + n)*rat(3*n*d^2 - 17*n*d + 24*n + 9*d^3 - 
         72*d^2 + 191*d - 168, - 432*n^2 + 432*n*d - 432*n - 432*d + 864)
          + m9denJ0021X00110( - 1 + n)*rat(n*d - 3*n + 3*d^2 - d - 24,72*n^2
          - 72*n*d + 72*n + 72*d - 144)
          + 1/( - 72 + 72*n)*m9denJ2010X01100( - 1 + n)*n*rat(-17,1)
          + 1/( - 72 + 72*n)*m9denJ2010X01100( - 1 + n)*rat(34,1)
          )

       + delta_( - 2 + n) * (
          + m9denJ0001X01110(0)*rat( - d^2 + 4*d - 4,144*d - 576)
          );

   id,ifmatch->sortme m9denJ1011X11110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m9denJ1011X11110( - 3 + n)*rat( - 2*n + 3*d - 4,540*n - 540)
          )

       + theta_( - 3 + n) * (
          + m9denJ0001X11110( - 2 + n)*rat(d - 2,270*n - 270)
          + m9denJ1001X11110( - 2 + n)*rat( - d + 3,270*n - 270)
          + m9denJ1011X01110( - 2 + n)*rat( - d + 2,270*n - 270)
          + m9denJ1011X11110( - 2 + n)*rat( - 10*n + 10*d - 9,135*n - 135)
          + 1/( - 135 + 135*n)*m9denJ0021X10110( - 2 + n)*rat(-1,1)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X11110( - 1 + n)*rat(7*d - 14,270*n - 270)
          + m9denJ0011X10110( - 1 + n)*rat(3*d - 8,270*n - 270)
          + m9denJ0011X11100( - 1 + n)*rat( - d + 2,270*n - 270)
          + m9denJ0011X11110( - 1 + n)*rat(2*d - 6,135*n - 135)
          + m9denJ1001X11110( - 1 + n)*rat( - d + 3,30*n - 30)
          + m9denJ1011X01110( - 1 + n)*rat( - d + 2,45*n - 45)
          + m9denJ1011X11110( - 1 + n)*rat( - 86*n + 43*d + 8,180*n - 180)
          + 1/( - 27 + 27*n)*m9denJ0021X10110( - 1 + n)*rat(-2,1)
          )

       + delta_( - 3 + n) * (
          + m9denJ1011X11110(0)*rat(3*d - 10,1080)
          )

       + delta_( - 2 + n) * (
          + m9denJ0011X01110(0)*rat(d^2 - 4*d + 4,864*d - 2592)
          + m9denJ0011X10110(0)*rat( - 4*d^2 + 24*d - 35,2160*d - 6480)
          + m9denJ0011X11110(0)*rat(3*d - 8,720)
          + m9denJ0031X10110(0)*rat(5,432*d - 1296)
          + m9denJ1001X11110(0)*rat( - d + 3,270)
          + m9denJ1011X11110(0)*rat(3*d - 10,90)
          + m9denJ2011X11110(0)*rat( - 1,90)
          );

   id,ifmatch->sortme m9denJ1011X01110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m9denJ1011X01110( - 3 + n)*rat( - n + d,270*n - 270)
          )

       + theta_( - 3 + n) * (
          + m9denJ0001X01110( - 2 + n)*rat(d - 2,270*n - 270)
          + m9denJ0001X11110( - 2 + n)*rat( - d + 3,270*n - 270)
          + m9denJ1011X01110( - 2 + n)*rat( - 8*n + 5*d + 4,108*n - 108)
          + 1/( - 135 + 135*n)*m9denJ0021X00110( - 2 + n)*rat(-1,1)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X01110( - 1 + n)*rat(7*d - 14,270*n - 270)
          + m9denJ0001X11110( - 1 + n)*rat( - d + 3,30*n - 30)
          + m9denJ0011X01100( - 1 + n)*rat( - d + 2,270*n - 270)
          + m9denJ0011X01110( - 1 + n)*rat(2*d - 6,135*n - 135)
          + m9denJ0011X00110( - 1 + n)*rat(3*d - 8,270*n - 270)
          + m9denJ1011X01110( - 1 + n)*rat( - 86*n + 25*d + 72,180*n - 180)
          + 1/( - 27 + 27*n)*m9denJ0021X00110( - 1 + n)*rat(-2,1)
          )

       + delta_( - 3 + n) * (
          + m9denJ1011X01110(0)*rat(d - 3,540)
          )

       + delta_( - 2 + n) * (
          + m9denJ0001X11110(0)*rat( - 3*d + 8,540)
          + m9denJ0011X01110(0)*rat(d - 2,135)
          + m9denJ1011X01110(0)*rat(d - 3,54)
          );

   id,ifmatch->sortme m9denJ1010X11100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0011X10110( - 2 + n)*rat(2*d - 6,45*n - 45)
          + 1/( - 45 + 45*n)*m9denJ1010X11100( - 2 + n)*n*rat(-1,1)
          + 1/( - 45 + 45*n)*m9denJ1010X11100( - 2 + n)*rat(3,1)
          + 1/( - 15 + 15*n)*m9denJ0021X10110( - 2 + n)*rat(-1,1)
          )

       + theta_( - 2 + n) * (
          + m9denJ0011X10110( - 1 + n)*rat(4*d - 12,15*n - 15)
          + m9denJ0011X00110( - 1 + n)*rat( - d + 2,45*n - 45)
          + 1/( - 45 + 45*n)*m9denJ1010X11100( - 1 + n)*n*rat(-14,1)
          + 1/( - 45 + 45*n)*m9denJ1010X11100( - 1 + n)*rat(28,1)
          + 1/( - 3 + 3*n)*m9denJ0021X10110( - 1 + n)*rat(-1,1)
          )

       + delta_( - 2 + n) * (
          + m9denJ0011X10110(0)*rat(4*d - 10,225)
          );

   id,ifmatch->sortme m9denJ1010X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m9denJ0011X00110( - 1 + n)*rat(3*d - 8,18*n - 18)
          + 1/( - 9 + 9*n)*m9denJ1010X01100( - 1 + n)*n*rat(-1,1)
          + 1/( - 9 + 9*n)*m9denJ1010X01100( - 1 + n)*rat(2,1)
          + 1/( - 3 + 3*n)*m9denJ0021X00110( - 1 + n)*rat(-1,1)
          );

   id,ifmatch->sortme m9denJ1001X11110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ1001X11110( - 2 + n)*rat( - n + 2*d - 4,45*n - 45)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X11110( - 1 + n)*rat( - d + 2,15*n - 15)
          + m9denJ1001X11110( - 1 + n)*rat( - 14*n + 16*d - 25,45*n - 45)
          )

       + delta_( - 2 + n) * (
          + m9denJ1001X11110(0)*rat(2*d - 6,45)
          );

   id,ifmatch->sortme m9denJ1001X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0001X11110( - 2 + n)*rat(3*d - 10,90*n - 90)
          + 1/( - 45 + 45*n)*m9denJ1001X01110( - 2 + n)*n*rat(-1,1)
          + 1/( - 45 + 45*n)*m9denJ1001X01110( - 2 + n)*rat(3,1)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X01110( - 1 + n)*rat( - 2*d + 4,45*n - 45)
          + m9denJ0001X11110( - 1 + n)*rat(23*d - 74,90*n - 90)
          + 1/( - 45 + 45*n)*m9denJ1001X01110( - 1 + n)*n*rat(-14,1)
          + 1/( - 45 + 45*n)*m9denJ1001X01110( - 1 + n)*rat(28,1)
          )

       + delta_( - 2 + n) * (
          + m9denJ0001X11110(0)*rat(3*d - 8,90)
          );

   id,ifmatch->sortme m9denJ1000X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m9denJ0001X00110( - 1 + n)*rat(d - 2,18*n - 18)
          + 1/( - 9 + 9*n)*m9denJ1000X01100( - 1 + n)*n*rat(-1,1)
          + 1/( - 9 + 9*n)*m9denJ1000X01100( - 1 + n)*rat(2,1)
          );

   id,ifmatch->sortme m9denJ0021X10110(n?{>1}) =

       + theta_( - 4 + n) * (
          + m9denJ0021X10110( - 3 + n)*rat( - 2*n + 3*d - 2,720*n - 720*d + 
         1440)
          )

       + theta_( - 3 + n) * (
          + m9denJ0001X01110( - 2 + n)*rat(d^3 - 5*d^2 + 8*d - 4, - 720*n^2 + 
         720*n*d - 720*n - 720*d + 1440)
          + m9denJ0011X10110( - 2 + n)*rat(3*n*d^2 - 17*n*d + 24*n + 12*d^3 - 
         95*d^2 + 249*d - 216, - 2160*n^2 + 2160*n*d - 2160*n - 2160*d + 4320)
          + m9denJ0021X00110( - 2 + n)*rat(d - 2, - 360*n + 360*d - 720)
          + m9denJ0021X10110( - 2 + n)*rat( - 44*n^2 + 63*n*d - 26*n + 6*d^2
          - 97*d + 118,720*n^2 - 720*n*d + 720*n + 720*d - 1440)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X01110( - 1 + n)*rat(7*d^3 - 35*d^2 + 56*d - 28, - 720*
         n^2 + 720*n*d - 720*n - 720*d + 1440)
          + m9denJ0001X00110( - 1 + n)*rat( - d^3 + 6*d^2 - 12*d + 8, - 720*
         n^2 + 720*n*d - 720*n - 720*d + 1440)
          + m9denJ0011X10110( - 1 + n)*rat(15*n*d^2 - 85*n*d + 120*n + 72*d^3
          - 555*d^2 + 1409*d - 1176, - 2160*n^2 + 2160*n*d - 2160*n - 2160*d
          + 4320)
          + m9denJ0011X00110( - 1 + n)*rat( - 3*d^3 + 23*d^2 - 58*d + 48, - 
         1080*n^2 + 1080*n*d - 1080*n - 1080*d + 2160)
          + m9denJ0021X00110( - 1 + n)*rat(d - 2, - 45*n + 45*d - 90)
          + m9denJ0021X10110( - 1 + n)*rat( - 157*n^2 + 198*n*d - 143*n + 15*
         d^2 - 283*d + 420,360*n^2 - 360*n*d + 360*n + 360*d - 720)
          )

       + delta_( - 3 + n) * (
          + m9denJ0011X10110(0)*rat( - 6*d^2 + 31*d - 40,3600*d - 18000)
          )

       + delta_( - 2 + n) * (
          + m9denJ0001X01110(0)*rat(d^3 - 6*d^2 + 12*d - 8,720*d - 2880)
          + m9denJ0011X10110(0)*rat(24*d^3 - 394*d^2 + 1555*d - 1800,10800*d
          - 43200)
          );

   id,ifmatch->sortme m9denJ0021X00110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0021X00110( - 2 + n)*rat( - n + d,72*n - 72*d + 144)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X00110( - 1 + n)*rat(d^3 - 4*d^2 + 4*d, - 288*n^2 + 288*
         n*d - 288*n - 288*d + 576)
          + m9denJ0011X00110( - 1 + n)*rat(3*n*d^2 - 17*n*d + 24*n + 9*d^3 - 
         72*d^2 + 191*d - 168, - 432*n^2 + 432*n*d - 432*n - 432*d + 864)
          + m9denJ0021X00110( - 1 + n)*rat( - 17*n^2 + 18*n*d - 3*n + 3*d^2 - 
         35*d + 44,72*n^2 - 72*n*d + 72*n + 72*d - 144)
          )

       + delta_( - 2 + n) * (
          + m9denJ0001X01110(0)*rat( - d^2 + 4*d - 4,144*d - 576)
          );

   id,ifmatch->sortme m9denJ0011X11110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0011X11110( - 2 + n)*rat( - 2*n + 3*d - 4,90*n - 90)
          )

       + theta_( - 2 + n) * (
          + m9denJ0011X01110( - 1 + n)*rat( - d + 2,45*n - 45)
          + m9denJ0011X10110( - 1 + n)*rat(3*d - 8,90*n - 90)
          + m9denJ0011X11100( - 1 + n)*rat( - d + 2,90*n - 90)
          + m9denJ0011X11110( - 1 + n)*rat( - 28*n + 23*d - 18,90*n - 90)
          + 1/( - 9 + 9*n)*m9denJ0021X10110( - 1 + n)*rat(-1,1)
          )

       + delta_( - 2 + n) * (
          + m9denJ0011X11110(0)*rat(3*d - 8,90)
          );

   id,ifmatch->sortme m9denJ0011X11100(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0011X11100( - 2 + n)*rat( - n + d,45*n - 45)
          )

       + theta_( - 2 + n) * (
          + m9denJ0011X01100( - 1 + n)*rat( - d + 2,45*n - 45)
          + m9denJ0011X11100( - 1 + n)*rat( - 14*n + 7*d + 7,45*n - 45)
          )

       + delta_( - 2 + n) * (
          + m9denJ0011X01110(0)*rat(d - 2,45)
          );

   id,ifmatch->sortme m9denJ0011X10110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0011X10110( - 2 + n)*rat( - n + 2*d - 3,45*n - 45)
          + 1/( - 15 + 15*n)*m9denJ0021X10110( - 2 + n)*rat(-1,1)
          )

       + theta_( - 2 + n) * (
          + m9denJ0011X10110( - 1 + n)*rat( - 14*n + 12*d - 8,45*n - 45)
          + m9denJ0011X00110( - 1 + n)*rat( - d + 2,45*n - 45)
          + 1/( - 3 + 3*n)*m9denJ0021X10110( - 1 + n)*rat(-1,1)
          )

       + delta_( - 2 + n) * (
          + m9denJ0011X10110(0)*rat(4*d - 10,225)
          );

   id,ifmatch->sortme m9denJ0011X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0011X01110( - 2 + n)*rat( - n + d,45*n - 45)
          )

       + theta_( - 2 + n) * (
          + m9denJ0011X01100( - 1 + n)*rat( - d + 2,90*n - 90)
          + m9denJ0011X01110( - 1 + n)*rat( - 14*n + 7*d + 7,45*n - 45)
          + m9denJ0011X00110( - 1 + n)*rat(3*d - 8,90*n - 90)
          + 1/( - 9 + 9*n)*m9denJ0021X00110( - 1 + n)*rat(-1,1)
          )

       + delta_( - 2 + n) * (
          + m9denJ0011X01110(0)*rat(d - 2,45)
          );

   id,ifmatch->sortme m9denJ0011X01100(n?{>1}) =

       + theta_( - 2 + n) * (
          + m9denJ0011X01100( - 1 + n)*rat( - 2*n + d + 2,18*n - 18)
          );

   id,ifmatch->sortme m9denJ0011X00110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m9denJ0011X00110( - 1 + n)*rat( - 2*n + 3*d - 4,18*n - 18)
          + 1/( - 3 + 3*n)*m9denJ0021X00110( - 1 + n)*rat(-1,1)
          );

   id,ifmatch->sortme m9denJ0001X11110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0001X11110( - 2 + n)*rat( - 2*n + 3*d - 4,90*n - 90)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X01110( - 1 + n)*rat( - 2*d + 4,45*n - 45)
          + m9denJ0001X11110( - 1 + n)*rat( - 28*n + 23*d - 18,90*n - 90)
          )

       + delta_( - 2 + n) * (
          + m9denJ0001X11110(0)*rat(3*d - 8,90)
          );

   id,ifmatch->sortme m9denJ0001X01110(n?{>1}) =

       + theta_( - 3 + n) * (
          + m9denJ0001X01110( - 2 + n)*rat( - n + d,45*n - 45)
          )

       + theta_( - 2 + n) * (
          + m9denJ0001X01110( - 1 + n)*rat( - 14*n + 7*d + 7,45*n - 45)
          + m9denJ0001X00110( - 1 + n)*rat( - d + 2,45*n - 45)
          )

       + delta_( - 2 + n) * (
          + m9denJ0001X01110(0)*rat(d - 2,45)
          );

   id,ifmatch->sortme m9denJ0001X00110(n?{>1}) =

       + theta_( - 2 + n) * (
          + m9denJ0001X00110( - 1 + n)*rat( - 2*n + d + 2,18*n - 18)
          );

********************************************************************************                 

goto endrec;

la sortme;
$irep = 0;
la endrec;

ModuleOption,minimum,$irep;

.sort:rec1d-m9-`$repcount++';

#redefine irep "`$irep'"
#enddo
#endprocedure
*--#] rec1dm9 :

*--#[ rec1dnum :
#procedure rec1dnum
#$repcount = 1;

#do irep = 1,1
#$irep = 1;
                
********************************************************************************                 
   id,ifmatch->sortme numJ2010X11100(n?pos_) =

       + theta_( - 3 + n) * (
          + numJ0011X10110( - 3 + n)*rat(288 - 204*d + 36*d^2, - 8 + 3*d + 2*n
         )
          + numJ0021X00110( - 3 + n)*rat( - 36 + 18*d, - 8 + 3*d + 2*n)
          + numJ0021X10110( - 3 + n)*rat( - 540 + 144*d + 540*n - 180*n*d, - 8
          + 3*d - 6*n + 3*n*d + 2*n^2)
          - 72/(1 + n)*numJ2010X11100( - 3 + n)
          + 36/(1 + n)*numJ2010X11100( - 3 + n)*n
          )

       + theta_( - 2 + n) * (
          + numJ0001X01110( - 2 + n)*rat( - 16 + 16*d - 4*d^2, - 8 + 3*d + 2*n
         )
          + numJ0011X10110( - 2 + n)*rat( - 168 + 119*d - 21*d^2, - 8 + 3*d + 
         2*n)
          + numJ0021X00110( - 2 + n)*rat(40 - 20*d, - 8 + 3*d + 2*n)
          + numJ0021X10110( - 2 + n)*rat(512 - 138*d - 468*n + 156*n*d, - 8 + 
         3*d - 6*n + 3*n*d + 2*n^2)
          + 49/(1 + n)*numJ2010X11100( - 2 + n)
          - 49/(1 + n)*numJ2010X11100( - 2 + n)*n
          )

       + theta_( - 1 + n) * (
          + numJ0001X01110( - 1 + n)*rat(4 - 4*d + d^2, - 8 + 3*d + 2*n)
          + numJ0011X10110( - 1 + n)*rat(24 - 17*d + 3*d^2, - 8 + 3*d + 2*n)
          + numJ0021X00110( - 1 + n)*rat( - 4 + 2*d, - 8 + 3*d + 2*n)
          + numJ0021X10110( - 1 + n)*rat( - 68 + 18*d + 72*n - 24*n*d, - 8 + 3
         *d - 6*n + 3*n*d + 2*n^2)
          + 14/(1 + n)*numJ2010X11100( - 1 + n)*n
          )

       + delta_( - 2 + n) * (
          + numJ0011X10110(0)*rat(30 - 42*d + 12*d^2, - 20 + 15*d)
          + numJ0031X10110(0)*rat(30, - 4 + 3*d)
          )

       + delta_( - 1 + n) * (
          + numJ0011X10110(0)*rat( - 280 + 217*d - 42*d^2, - 30 + 15*d)
          );

   id,ifmatch->sortme numJ2010X01100(n?pos_) =

       + theta_( - 2 + n) * (
          + numJ0011X00110( - 2 + n)*rat( - 72 + 51*d - 9*d^2, - 4 + 2*d + 2*n
         )
          + numJ0021X00110( - 2 + n)*rat( - 54*n + 18*n*d, - 2 + d - n + n*d
          + n^2)
          + 9/(1 + n)*numJ2010X01100( - 2 + n)
          - 9/(1 + n)*numJ2010X01100( - 2 + n)*n
          )

       + theta_( - 1 + n) * (
          + numJ0001X00110( - 1 + n)*rat(4 - 4*d + d^2, - 4 + 2*d + 2*n)
          + numJ0011X00110( - 1 + n)*rat(24 - 17*d + 3*d^2, - 4 + 2*d + 2*n)
          + numJ0021X00110( - 1 + n)*rat(30*n - 10*n*d, - 2 + d - n + n*d + 
         n^2)
          + 10/(1 + n)*numJ2010X01100( - 1 + n)*n
          );

   id,ifmatch->sortme numJ1011X11110(n?pos_) =

       + theta_( - 2 + n) * (
          + numJ0001X11110( - 2 + n)*rat( - 8 + 4*d, - 10 + 3*d + 2*n)
          + numJ0011X10110( - 2 + n)*rat(16 - 6*d, - 10 + 3*d + 2*n)
          + numJ0011X11100( - 2 + n)*rat( - 4 + 2*d, - 10 + 3*d + 2*n)
          + numJ0011X11110( - 2 + n)*rat(24 - 8*d, - 10 + 3*d + 2*n)
          + numJ0021X10110( - 2 + n)*rat(4, - 10 + 3*d + 2*n)
          + numJ1011X01110( - 2 + n)*rat(12 - 6*d, - 10 + 3*d + 2*n)
          + numJ1011X11110( - 2 + n)*rat(84 - 12*d - 24*n, - 10 + 3*d + 2*n)
          )

       + theta_( - 1 + n) * (
          + numJ0001X11110( - 1 + n)*rat(4 - 2*d, - 10 + 3*d + 2*n)
          + numJ0021X10110( - 1 + n)*rat(4, - 10 + 3*d + 2*n)
          + numJ1001X11110( - 1 + n)*rat( - 6 + 2*d, - 10 + 3*d + 2*n)
          + numJ1011X01110( - 1 + n)*rat( - 4 + 2*d, - 10 + 3*d + 2*n)
          + numJ1011X11110( - 1 + n)*rat( - 60 + 14*d + 14*n, - 10 + 3*d + 2*n
         )
          )

       + delta_( - 1 + n) * (
          + numJ0011X01110(0)*rat(20 - 20*d + 5*d^2,192 - 136*d + 24*d^2)
          + numJ0011X10110(0)*rat( - 35 + 24*d - 4*d^2,96 - 68*d + 12*d^2)
          + numJ0011X11110(0)*rat( - 24 + 9*d, - 32 + 12*d)
          + numJ0031X10110(0)*rat(25,96 - 68*d + 12*d^2)
          + numJ1001X11110(0)*rat(6 - 2*d, - 8 + 3*d)
          + numJ1011X11110(0)*rat(30 - 9*d, - 8 + 3*d)
          + numJ2011X11110(0)*rat(-6, - 8 + 3*d)
          );

   id,ifmatch->sortme numJ1011X01110(n?pos_) =

       + theta_( - 2 + n) * (
          + numJ0001X01110( - 2 + n)*rat( - 4 + 2*d, - 3 + d + n)
          + numJ0011X01100( - 2 + n)*rat( - 2 + d, - 3 + d + n)
          + numJ0011X01110( - 2 + n)*rat(12 - 4*d, - 3 + d + n)
          + numJ0011X00110( - 2 + n)*rat(8 - 3*d, - 3 + d + n)
          + numJ0021X00110( - 2 + n)*rat(2, - 3 + d + n)
          + numJ1011X01110( - 2 + n)*rat(36 - 6*d - 12*n, - 3 + d + n)
          )

       + theta_( - 1 + n) * (
          + numJ0001X01110( - 1 + n)*rat(2 - d, - 3 + d + n)
          + numJ0001X11110( - 1 + n)*rat( - 3 + d, - 3 + d + n)
          + numJ0021X00110( - 1 + n)*rat(2, - 3 + d + n)
          + numJ1011X01110( - 1 + n)*rat( - 44 + 11*d + 14*n, - 6 + 2*d + 2*n)
          )

       + delta_( - 1 + n) * (
          + numJ0001X11110(0)*rat(8 - 3*d, - 4 + 2*d)
          + numJ0011X01110(0)*rat( - 4 + 2*d, - 2 + d)
          + numJ1011X01110(0)*rat(12 - 4*d, - 2 + d)
          );

   id,ifmatch->sortme numJ1010X11100(n?pos_) =

       + theta_( - 3 + n) * (
          + numJ0011X10110( - 3 + n)*rat(864 - 612*d + 108*d^2,40 - 31*d + 6*
         d^2 - 18*n + 7*n*d + 2*n^2)
          + numJ0021X00110( - 3 + n)*rat( - 108 + 54*d,40 - 31*d + 6*d^2 - 18*
         n + 7*n*d + 2*n^2)
          + numJ0021X10110( - 3 + n)*rat(108 - 216*d + 216*n,40 - 31*d + 6*d^2
          - 18*n + 7*n*d + 2*n^2)
          )

       + theta_( - 2 + n) * (
          + numJ0001X01110( - 2 + n)*rat( - 48 + 48*d - 12*d^2,40 - 31*d + 6*
         d^2 - 18*n + 7*n*d + 2*n^2)
          + numJ0011X10110( - 2 + n)*rat( - 504 + 357*d - 63*d^2,40 - 31*d + 6
         *d^2 - 18*n + 7*n*d + 2*n^2)
          + numJ0021X00110( - 2 + n)*rat(120 - 60*d,40 - 31*d + 6*d^2 - 18*n
          + 7*n*d + 2*n^2)
          + numJ0021X10110( - 2 + n)*rat(360 + 27*d - 294*n,40 - 31*d + 6*d^2
          - 18*n + 7*n*d + 2*n^2)
          )

       + theta_( - 1 + n) * (
          + numJ0001X01110( - 1 + n)*rat(12 - 12*d + 3*d^2,40 - 31*d + 6*d^2
          - 18*n + 7*n*d + 2*n^2)
          + numJ0011X10110( - 1 + n)*rat(216 - 153*d + 27*d^2 - 12*n - 5*n*d
          + 3*n*d^2 + 12*n^2 - 4*n^2*d,40 - 31*d + 6*d^2 + 22*n - 24*n*d + 6*n
         *d^2 - 16*n^2 + 7*n^2*d + 2*n^3)
          + numJ0011X00110( - 1 + n)*rat( - 2 + d, - 5 + 2*d + n)
          + numJ0021X00110( - 1 + n)*rat( - 12 + 6*d,40 - 31*d + 6*d^2 - 18*n
          + 7*n*d + 2*n^2)
          + numJ0021X10110( - 1 + n)*rat( - 108 + 18*d + 60*n,40 - 31*d + 6*
         d^2 - 18*n + 7*n*d + 2*n^2)
          + 4/(1 + n)*numJ1010X11100( - 1 + n)*n
          )

       + delta_( - 2 + n) * (
          + numJ0011X10110(0)*rat(90 - 126*d + 36*d^2,60 - 85*d + 30*d^2)
          + numJ0031X10110(0)*rat(90,12 - 17*d + 6*d^2)
          )

       + delta_( - 1 + n) * (
          + numJ0011X10110(0)*rat( - 840 + 651*d - 126*d^2,120 - 120*d + 30*
         d^2)
          );

   id,ifmatch->sortme numJ1010X01100(n?pos_) =

       + theta_( - 2 + n) * (
          + numJ0011X00110( - 2 + n)*rat( - 216 + 153*d - 27*d^2,12 - 12*d + 3
         *d^2 - 10*n + 5*n*d + 2*n^2)
          + numJ0021X00110( - 2 + n)*rat( - 108 + 54*d - 54*n,12 - 12*d + 3*
         d^2 - 10*n + 5*n*d + 2*n^2)
          )

       + theta_( - 1 + n) * (
          + numJ0001X00110( - 1 + n)*rat(12 - 12*d + 3*d^2,12 - 12*d + 3*d^2
          - 10*n + 5*n*d + 2*n^2)
          + numJ0011X00110( - 1 + n)*rat(72 - 51*d + 9*d^2,12 - 12*d + 3*d^2
          - 10*n + 5*n*d + 2*n^2)
          + numJ0021X00110( - 1 + n)*rat(60*n,12 - 12*d + 3*d^2 - 10*n + 5*n*d
          + 2*n^2)
          );

   id,ifmatch->sortme numJ1001X11110(n?pos_) =

       + theta_( - 1 + n) * (
          + numJ0001X11110( - 1 + n)*rat( - 6 + 3*d, - 6 + 2*d + n)
          + numJ1001X11110( - 1 + n)*rat( - 10 + 2*d + 4*n, - 6 + 2*d + n)
          );

   id,ifmatch->sortme numJ1001X01110(n?pos_) =

       + theta_( - 1 + n) * (
          + numJ0001X01110( - 1 + n)*rat( - 8 + 4*d, - 8 + 3*d + 2*n)
          + numJ0001X11110( - 1 + n)*rat( - 16 + 4*d + 24*n - 8*n*d, - 8 + 3*d
          - 6*n + 3*n*d + 2*n^2)
          + 4/(1 + n)*numJ1001X01110( - 1 + n)*n
          );

   id,ifmatch->sortme numJ1000X01100(n?pos_) = 0;

   id,ifmatch->sortme numJ0021X10110(n?pos_) =

       + theta_( - 3 + n) * (
          + numJ0011X10110( - 3 + n)*rat(288 - 204*d + 36*d^2, - 8 + 3*d + 2*n
         )
          + numJ0021X00110( - 3 + n)*rat( - 36 + 18*d, - 8 + 3*d + 2*n)
          + numJ0021X10110( - 3 + n)*rat(36 - 72*d + 72*n, - 8 + 3*d + 2*n)
          )

       + theta_( - 2 + n) * (
          + numJ0001X01110( - 2 + n)*rat( - 16 + 16*d - 4*d^2, - 8 + 3*d + 2*n
         )
          + numJ0011X10110( - 2 + n)*rat( - 168 + 119*d - 21*d^2, - 8 + 3*d + 
         2*n)
          + numJ0021X00110( - 2 + n)*rat(40 - 20*d, - 8 + 3*d + 2*n)
          + numJ0021X10110( - 2 + n)*rat(120 + 9*d - 98*n, - 8 + 3*d + 2*n)
          )

       + theta_( - 1 + n) * (
          + numJ0001X01110( - 1 + n)*rat(4 - 4*d + d^2, - 8 + 3*d + 2*n)
          + numJ0011X10110( - 1 + n)*rat(24 - 17*d + 3*d^2, - 8 + 3*d + 2*n)
          + numJ0021X00110( - 1 + n)*rat( - 4 + 2*d, - 8 + 3*d + 2*n)
          + numJ0021X10110( - 1 + n)*rat( - 68 + 18*d + 28*n, - 8 + 3*d + 2*n)
          )

       + delta_( - 2 + n) * (
          + numJ0011X10110(0)*rat(30 - 42*d + 12*d^2, - 20 + 15*d)
          + numJ0031X10110(0)*rat(30, - 4 + 3*d)
          )

       + delta_( - 1 + n) * (
          + numJ0011X10110(0)*rat( - 280 + 217*d - 42*d^2, - 30 + 15*d)
          );

   id,ifmatch->sortme numJ0021X00110(n?pos_) =

       + theta_( - 2 + n) * (
          + numJ0011X00110( - 2 + n)*rat( - 72 + 51*d - 9*d^2, - 4 + 2*d + 2*n
         )
          + numJ0021X00110( - 2 + n)*rat( - 18 + 9*d - 9*n, - 2 + d + n)
          )

       + theta_( - 1 + n) * (
          + numJ0001X00110( - 1 + n)*rat(4 - 4*d + d^2, - 4 + 2*d + 2*n)
          + numJ0011X00110( - 1 + n)*rat(24 - 17*d + 3*d^2, - 4 + 2*d + 2*n)
          + numJ0021X00110( - 1 + n)*rat(10*n, - 2 + d + n)
          );

   id,ifmatch->sortme numJ0011X11110(n?pos_) =

       + theta_( - 1 + n) * (
          + numJ0011X01110( - 1 + n)*rat( - 4 + 2*d, - 8 + 3*d + 2*n)
          + numJ0011X10110( - 1 + n)*rat(8 - 3*d, - 8 + 3*d + 2*n)
          + numJ0011X11100( - 1 + n)*rat( - 2 + d, - 8 + 3*d + 2*n)
          + numJ0011X11110( - 1 + n)*rat( - 16 + 4*d + 8*n, - 8 + 3*d + 2*n)
          + numJ0021X10110( - 1 + n)*rat(10, - 8 + 3*d + 2*n)
          );

   id,ifmatch->sortme numJ0011X11100(n?pos_) =

       + theta_( - 1 + n) * (
          + numJ0011X01100( - 1 + n)*rat( - 2 + d, - 2 + d + n)
          + numJ0011X11100( - 1 + n)*rat( - 6 + 2*d + 4*n, - 2 + d + n)
          );

   id,ifmatch->sortme numJ0011X10110(n?pos_) =

       + theta_( - 3 + n) * (
          + numJ0011X10110( - 3 + n)*rat(864 - 612*d + 108*d^2,40 - 31*d + 6*
         d^2 - 18*n + 7*n*d + 2*n^2)
          + numJ0021X00110( - 3 + n)*rat( - 108 + 54*d,40 - 31*d + 6*d^2 - 18*
         n + 7*n*d + 2*n^2)
          + numJ0021X10110( - 3 + n)*rat(108 - 216*d + 216*n,40 - 31*d + 6*d^2
          - 18*n + 7*n*d + 2*n^2)
          )

       + theta_( - 2 + n) * (
          + numJ0001X01110( - 2 + n)*rat( - 48 + 48*d - 12*d^2,40 - 31*d + 6*
         d^2 - 18*n + 7*n*d + 2*n^2)
          + numJ0011X10110( - 2 + n)*rat( - 504 + 357*d - 63*d^2,40 - 31*d + 6
         *d^2 - 18*n + 7*n*d + 2*n^2)
          + numJ0021X00110( - 2 + n)*rat(120 - 60*d,40 - 31*d + 6*d^2 - 18*n
          + 7*n*d + 2*n^2)
          + numJ0021X10110( - 2 + n)*rat(360 + 27*d - 294*n,40 - 31*d + 6*d^2
          - 18*n + 7*n*d + 2*n^2)
          )

       + theta_( - 1 + n) * (
          + numJ0001X01110( - 1 + n)*rat(12 - 12*d + 3*d^2,40 - 31*d + 6*d^2
          - 18*n + 7*n*d + 2*n^2)
          + numJ0011X10110( - 1 + n)*rat(216 - 153*d + 27*d^2 - 68*n + 24*n*d
          + 8*n^2,40 - 31*d + 6*d^2 - 18*n + 7*n*d + 2*n^2)
          + numJ0011X00110( - 1 + n)*rat( - 2 + d, - 5 + 2*d + n)
          + numJ0021X00110( - 1 + n)*rat( - 12 + 6*d,40 - 31*d + 6*d^2 - 18*n
          + 7*n*d + 2*n^2)
          + numJ0021X10110( - 1 + n)*rat( - 108 + 18*d + 60*n,40 - 31*d + 6*
         d^2 - 18*n + 7*n*d + 2*n^2)
          )

       + delta_( - 2 + n) * (
          + numJ0011X10110(0)*rat(90 - 126*d + 36*d^2,60 - 85*d + 30*d^2)
          + numJ0031X10110(0)*rat(90,12 - 17*d + 6*d^2)
          )

       + delta_( - 1 + n) * (
          + numJ0011X10110(0)*rat( - 840 + 651*d - 126*d^2,120 - 120*d + 30*
         d^2)
          );

   id,ifmatch->sortme numJ0011X01110(n?pos_) =

       + theta_( - 1 + n) * (
          + numJ0011X01100( - 1 + n)*rat( - 2 + d, - 4 + 2*d + 2*n)
          + numJ0011X01110( - 1 + n)*rat( - 6 + 2*d + 4*n, - 2 + d + n)
          + numJ0011X00110( - 1 + n)*rat(8 - 3*d, - 4 + 2*d + 2*n)
          + numJ0021X00110( - 1 + n)*rat(5, - 2 + d + n)
          );

   id,ifmatch->sortme numJ0011X01100(n?pos_) = 0;

   id,ifmatch->sortme numJ0011X00110(n?pos_) =

       + theta_( - 2 + n) * (
          + numJ0011X00110( - 2 + n)*rat( - 216 + 153*d - 27*d^2,12 - 12*d + 3
         *d^2 - 10*n + 5*n*d + 2*n^2)
          + numJ0021X00110( - 2 + n)*rat( - 108 + 54*d - 54*n,12 - 12*d + 3*
         d^2 - 10*n + 5*n*d + 2*n^2)
          )

       + theta_( - 1 + n) * (
          + numJ0001X00110( - 1 + n)*rat(12 - 12*d + 3*d^2,12 - 12*d + 3*d^2
          - 10*n + 5*n*d + 2*n^2)
          + numJ0011X00110( - 1 + n)*rat(72 - 51*d + 9*d^2,12 - 12*d + 3*d^2
          - 10*n + 5*n*d + 2*n^2)
          + numJ0021X00110( - 1 + n)*rat(60*n,12 - 12*d + 3*d^2 - 10*n + 5*n*d
          + 2*n^2)
          );

   id,ifmatch->sortme numJ0001X11110(n?pos_) =

       + theta_( - 1 + n) * (
          + numJ0001X01110( - 1 + n)*rat( - 8 + 4*d, - 8 + 3*d + 2*n)
          + numJ0001X11110( - 1 + n)*rat( - 16 + 4*d + 8*n, - 8 + 3*d + 2*n)
          );

   id,ifmatch->sortme numJ0001X01110(n?pos_) =

       + theta_( - 1 + n) * (
          + numJ0001X01110( - 1 + n)*rat( - 6 + 2*d + 4*n, - 2 + d + n)
          + numJ0001X00110( - 1 + n)*rat( - 2 + d, - 2 + d + n)
          );

   id,ifmatch->sortme numJ0001X00110(n?pos_) = 0;

********************************************************************************                 

goto endrec;

la sortme;
$irep = 0;
la endrec;

ModuleOption,minimum,$irep;

.sort:rec1d-num-`$repcount++';

#redefine irep "`$irep'"
#enddo
#endprocedure
*--#] rec1dnum :



*--#[ fmft :
#procedure fmft

*       Only if massive denominators present in TERM        
        if((count(d1,1) != 0) || (count(d2,1) != 0) || (count(d3,1) != 0) || (count(d4,1) != 0) || (count(d5,1) != 0) || (count(d6,1) != 0) || (count(d7,1) != 0) || (count(d8,1) != 0) || (count(d9,1) != 0) || (count(d10,1) != 0)) 
        Multiply intFMFT;


*       rewrite scalar products in terms of denominators        
        #call sp2den(FMFT)        

        if(count(intFMFT,1));

        if((count(d1,-1) <= 0) && (count(d10,-1) > 0));
        Multiply intX/intFMFT;

        elseif(match(1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?neg0_/d7^n7?pos_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_));
        Multiply replace_(d5,d3, d9,d4, d3,d5, d10,d6, d4,d9, d6,d10)*intH/intFMFT;

        elseif(match(1/d1^n1?pos_/d2^n2?pos_/d3^n3?pos_/d4^n4?pos_/d5^n5?pos_/d6^n6?pos_/d7^n7?neg0_/d8^n8?pos_/d9^n9?pos_/d10^n10?pos_));
        Multiply replace_(d3,d2, d5,d3, d8,d4, d2,d5, d10,d6, d6,d7, d9,d8, d4,d9, d7,d10)*intH/intFMFT;

        elseif(count(d10,-1) <= 0);
        Multiply intH/intFMFT;
        endif;

        endif;


        if(count(intFMFT,1) > 0);
        exit "We can not map topology with 9 or higher number of lines on H or X";
        endif;

*       First we remap to simpler topologies BMW and FG        
        #call mapBMW(X)
        #call mapFG(X)

*       Only after mapping on simpler we start to reduce top level topologies        
        #message
        #message reduction X started                 
        #message
        #call redX
        #call mapBMW(X)
        #call mapFG(X)

        #call mapBMW(H)
        #call mapFG(H)
        #printtimes


        #message
        #message reduction H started                 
        #message
        #call redH
        #printtimes

        #call mapBMW(H)
        #call mapFG(H)
        #message
        #message reduction BMW started                 
        #message
        #call redBMW
        #call mapFG(BMW)
        #printtimes
        
        if(count(intX,1,intH,1,intBMW,1) > 0);
        exit "High-level topologies X,H,BMW not reduced";
        endif;
        .sort:redX-H-BMW;

*       TODO:        
*       Improved power expand in numerator        
        id num(x?,y?)=x^y;
        .sort:mapFG;


        
        #message
        #message reduction FG started                 
        #message
        #call redFGnew
        #printtimes

*       Numerator reduction is based on topology with (pp^n5)        

*       Integrals with denominator with power n5 now have argument (n5), n5 > 0
*       If there is no denominator integral is ...num(0)       


        if(count(intFG,1));
*       Rewrite in numerater everithing in terms of pp        
        id [pp-tmm]^n5?pos_ = (pp - tmm)^n5; 
*       Set mass to one        
        Multiply replace_(tmm,1);        
        endif;

*         
*       FxG
*         
        id intFG*TFI(1,1,1,1,1)*Gx11/pp^n5?pos_          = int1d*m0denJ1011X11110(n5);
        id intFG*TFI(1,1,1,1,1)*Gx11/[pp-tmm]^n5?pos_    = int1d*m1denJ1011X11110(n5);
        id intFG*TFI(1,1,1,1,1)*Gx11/[pp-3*tmm]^n5?pos_  = int1d*m3denJ1011X11110(n5);
        id intFG*TFI(1,1,1,1,1)*Gx11/[pp-4*tmm]^n5?pos_  = int1d*m4denJ1011X11110(n5);
        id intFG*TFI(1,1,1,1,1)*Gx11/[pp-9*tmm]^n5?pos_  = int1d*m9denJ1011X11110(n5);
        id intFG*TFI(1,1,1,1,1)*Gx11*pp^n5?pos0_         =   int1d*numJ1011X11110(n5);
*         
*       VxG
*         
        id intFG*TFI(0,1,1,1,1)*Gx11/pp^n5?pos_          = int1d*m0denJ0011X11110(n5);
        id intFG*TFI(0,1,1,1,1)*Gx11/[pp-tmm]^n5?pos_    = int1d*m1denJ0011X11110(n5);
        id intFG*TFI(0,1,1,1,1)*Gx11/[pp-3*tmm]^n5?pos_  = int1d*m3denJ0011X11110(n5);
        id intFG*TFI(0,1,1,1,1)*Gx11/[pp-4*tmm]^n5?pos_  = int1d*m4denJ0011X11110(n5);
        id intFG*TFI(0,1,1,1,1)*Gx11/[pp-9*tmm]^n5?pos_  = int1d*m9denJ0011X11110(n5);
        id intFG*TFI(0,1,1,1,1)*Gx11*pp^n5?pos0_         =   int1d*numJ0011X11110(n5);
*         
*       J2xG
*         
        id intFG*TFI(2,0,0,1,1)*Gx11/pp^n5?pos_          = int1d*m0denJ2010X11100(n5);
        id intFG*TFI(2,0,0,1,1)*Gx11/[pp-tmm]^n5?pos_    = int1d*m1denJ2010X11100(n5);
        id intFG*TFI(2,0,0,1,1)*Gx11/[pp-3*tmm]^n5?pos_  = int1d*m3denJ2010X11100(n5);
        id intFG*TFI(2,0,0,1,1)*Gx11/[pp-4*tmm]^n5?pos_  = int1d*m4denJ2010X11100(n5);
        id intFG*TFI(2,0,0,1,1)*Gx11/[pp-9*tmm]^n5?pos_  = int1d*m9denJ2010X11100(n5);
        id intFG*TFI(2,0,0,1,1)*Gx11*pp^n5?pos0_         =   int1d*numJ2010X11100(n5);
*         
*       J1xG        
*         
        id intFG*TFI(1,0,0,1,1)*Gx11/pp^n5?pos_          = int1d*m0denJ1010X11100(n5);
        id intFG*TFI(1,0,0,1,1)*Gx11/[pp-tmm]^n5?pos_    = int1d*m1denJ1010X11100(n5);
        id intFG*TFI(1,0,0,1,1)*Gx11/[pp-3*tmm]^n5?pos_  = int1d*m3denJ1010X11100(n5);
        id intFG*TFI(1,0,0,1,1)*Gx11/[pp-4*tmm]^n5?pos_  = int1d*m4denJ1010X11100(n5);
        id intFG*TFI(1,0,0,1,1)*Gx11/[pp-9*tmm]^n5?pos_  = int1d*m9denJ1010X11100(n5);
        id intFG*TFI(1,0,0,1,1)*Gx11*pp^n5?pos0_         =   int1d*numJ1010X11100(n5);
*         
*       FxT1        
*         
        id intFG*TFI(1,1,1,1,1)*T1x1/pp^n5?pos_          = int1d*m0denJ1011X01110(n5);
        id intFG*TFI(1,1,1,1,1)*T1x1/[pp-tmm]^n5?pos_    = int1d*m1denJ1011X01110(n5);
        id intFG*TFI(1,1,1,1,1)*T1x1/[pp-3*tmm]^n5?pos_  = int1d*m3denJ1011X01110(n5);
        id intFG*TFI(1,1,1,1,1)*T1x1/[pp-4*tmm]^n5?pos_  = int1d*m4denJ1011X01110(n5);
        id intFG*TFI(1,1,1,1,1)*T1x1/[pp-9*tmm]^n5?pos_  = int1d*m9denJ1011X01110(n5);
        id intFG*TFI(1,1,1,1,1)*T1x1*pp^n5?pos0_         =   int1d*numJ1011X01110(n5);
*         
*       VxT1
*         
        id intFG*TFI(0,1,1,1,1)*T1x1/pp^n5?pos_          = int1d*m0denJ0011X01110(n5);
        id intFG*TFI(0,1,1,1,1)*T1x1/[pp-tmm]^n5?pos_    = int1d*m1denJ0011X01110(n5);
        id intFG*TFI(0,1,1,1,1)*T1x1/[pp-3*tmm]^n5?pos_  = int1d*m3denJ0011X01110(n5);
        id intFG*TFI(0,1,1,1,1)*T1x1/[pp-4*tmm]^n5?pos_  = int1d*m4denJ0011X01110(n5);
        id intFG*TFI(0,1,1,1,1)*T1x1/[pp-9*tmm]^n5?pos_  = int1d*m9denJ0011X01110(n5);
        id intFG*TFI(0,1,1,1,1)*T1x1*pp^n5?pos0_         =   int1d*numJ0011X01110(n5);
*         
*       J2xT1        
*         
        id intFG*TFI(2,0,0,1,1)*T1x1/pp^n5?pos_          = int1d*m0denJ2010X01100(n5);
        id intFG*TFI(2,0,0,1,1)*T1x1/[pp-tmm]^n5?pos_    = int1d*m1denJ2010X01100(n5);
        id intFG*TFI(2,0,0,1,1)*T1x1/[pp-3*tmm]^n5?pos_  = int1d*m3denJ2010X01100(n5);
        id intFG*TFI(2,0,0,1,1)*T1x1/[pp-4*tmm]^n5?pos_  = int1d*m4denJ2010X01100(n5);
        id intFG*TFI(2,0,0,1,1)*T1x1/[pp-9*tmm]^n5?pos_  = int1d*m9denJ2010X01100(n5);
        id intFG*TFI(2,0,0,1,1)*T1x1*pp^n5?pos0_         =   int1d*numJ2010X01100(n5);
*         
*       J1xT1        
*         
        id intFG*TFI(1,0,0,1,1)*T1x1/pp^n5?pos_          = int1d*m0denJ1010X01100(n5);
        id intFG*TFI(1,0,0,1,1)*T1x1/[pp-tmm]^n5?pos_    = int1d*m1denJ1010X01100(n5);
        id intFG*TFI(1,0,0,1,1)*T1x1/[pp-3*tmm]^n5?pos_  = int1d*m3denJ1010X01100(n5);
        id intFG*TFI(1,0,0,1,1)*T1x1/[pp-4*tmm]^n5?pos_  = int1d*m4denJ1010X01100(n5);
        id intFG*TFI(1,0,0,1,1)*T1x1/[pp-9*tmm]^n5?pos_  = int1d*m9denJ1010X01100(n5);
        id intFG*TFI(1,0,0,1,1)*T1x1*pp^n5?pos0_         =   int1d*numJ1010X01100(n5);
*         
*       GxGxG 
*         
        id intFG*Gx11^3/pp^n5?pos_                       = int1d*m0denJ1001X11110(n5);
        id intFG*Gx11^3/[pp-tmm]^n5?pos_                 = int1d*m1denJ1001X11110(n5);
        id intFG*Gx11^3/[pp-3*tmm]^n5?pos_               = int1d*m3denJ1001X11110(n5);
        id intFG*Gx11^3/[pp-4*tmm]^n5?pos_               = int1d*m4denJ1001X11110(n5);
        id intFG*Gx11^3/[pp-9*tmm]^n5?pos_               = int1d*m9denJ1001X11110(n5);
        id intFG*Gx11^3*pp^n5?pos0_                      =   int1d*numJ1001X11110(n5);
*         
*       GxGxT1
*         
        id intFG*Gx11^2*T1x1/pp^n5?pos_                  = int1d*m0denJ1001X01110(n5);
        id intFG*Gx11^2*T1x1/[pp-tmm]^n5?pos_            = int1d*m1denJ1001X01110(n5);
        id intFG*Gx11^2*T1x1/[pp-3*tmm]^n5?pos_          = int1d*m3denJ1001X01110(n5);
        id intFG*Gx11^2*T1x1/[pp-4*tmm]^n5?pos_          = int1d*m4denJ1001X01110(n5);
        id intFG*Gx11^2*T1x1/[pp-9*tmm]^n5?pos_          = int1d*m9denJ1001X01110(n5);
        id intFG*Gx11^2*T1x1*pp^n5?pos0_                 =   int1d*numJ1001X01110(n5);
*         
*       GxT1xT1
*         
        id intFG*Gx11*T1x1^2/pp^n5?pos_                  = int1d*m0denJ0001X01110(n5);
        id intFG*Gx11*T1x1^2/[pp-tmm]^n5?pos_            = int1d*m1denJ0001X01110(n5);
        id intFG*Gx11*T1x1^2/[pp-3*tmm]^n5?pos_          = int1d*m3denJ0001X01110(n5);
        id intFG*Gx11*T1x1^2/[pp-4*tmm]^n5?pos_          = int1d*m4denJ0001X01110(n5);
        id intFG*Gx11*T1x1^2/[pp-9*tmm]^n5?pos_          = int1d*m9denJ0001X01110(n5);
        id intFG*Gx11*T1x1^2*pp^n5?pos0_                 =   int1d*numJ0001X01110(n5);   
*         
*       T2xG
*         
        id intFG*T2x111*Gx11/pp^n5?pos_                  = int1d*m0denJ0011X11100(n5);
        id intFG*T2x111*Gx11/[pp-tmm]^n5?pos_            = int1d*m1denJ0011X11100(n5);
        id intFG*T2x111*Gx11/[pp-3*tmm]^n5?pos_          = int1d*m3denJ0011X11100(n5);
        id intFG*T2x111*Gx11/[pp-4*tmm]^n5?pos_          = int1d*m4denJ0011X11100(n5);
        id intFG*T2x111*Gx11/[pp-9*tmm]^n5?pos_          = int1d*m9denJ0011X11100(n5);
        id intFG*T2x111*Gx11*pp^n5?pos0_                 =   int1d*numJ0011X11100(n5);
*         
*       T2xT1
*         
        id intFG*T2x111*T1x1/pp^n5?pos_                  = int1d*m0denJ0011X01100(n5);
        id intFG*T2x111*T1x1/[pp-tmm]^n5?pos_            = int1d*m1denJ0011X01100(n5);
        id intFG*T2x111*T1x1/[pp-3*tmm]^n5?pos_          = int1d*m3denJ0011X01100(n5);
        id intFG*T2x111*T1x1/[pp-4*tmm]^n5?pos_          = int1d*m4denJ0011X01100(n5);
        id intFG*T2x111*T1x1/[pp-9*tmm]^n5?pos_          = int1d*m9denJ0011X01100(n5);
        id intFG*T2x111*T1x1*pp^n5?pos0_                 =   int1d*numJ0011X01100(n5);
*         
*       T1xT1xT1
*         
        id intFG*T1x1^3/pp^n5?pos_                       = int1d*m0denJ1000X01100(n5);
        id intFG*T1x1^3/[pp-tmm]^n5?pos_                 = int1d*m1denJ1000X01100(n5);
        id intFG*T1x1^3/[pp-3*tmm]^n5?pos_               = int1d*m3denJ1000X01100(n5);
        id intFG*T1x1^3/[pp-4*tmm]^n5?pos_               = int1d*m4denJ1000X01100(n5);
        id intFG*T1x1^3/[pp-9*tmm]^n5?pos_               = int1d*m9denJ1000X01100(n5);
        id intFG*T1x1^3*pp^n5?pos0_                      =   int1d*numJ1000X01100(n5);

        .sort:toJ;

        #message
        #message Started 1d reduction        
        #message

        #call rec1dm9
        #call rec1dm4
        #call rec1dm3
        #call rec1dm1
        #call rec1dm0
        #call rec1dnum
        #printtimes
        .sort:rec1d-done;


*       Substituting reduction rules for integrals with [pp-9*tmm]
        id int1d*m9denJ0001X01110(0) = PR1;
        id int1d*m9denJ0001X11110(0) = PR2;
        id int1d*m9denJ0011X01100(0) = 0;
        id int1d*m9denJ0011X01110(0) = PR3;
        id int1d*m9denJ0011X10110(0) = PR4; 
        id int1d*m9denJ0011X11100(0) = PR3;
        id int1d*m9denJ0011X11110(0) = PR7;
        id int1d*m9denJ0021X00110(1) = PR1*rat(-(-2 + d)^2, 16*(-3 + d)) + int1d*m9denJ0001X00110(1)*rat((-2 + d)^2, 16*(-3 + d)) + int1d*m9denJ0011X00110(1)*rat(-8 + 3*d, 24);
        id int1d*m9denJ0021X10110(1) = PR4d*rat(-1, 24*(-3 + d)) + int1d*m9denJ0001X01110(1)*rat((-2 + d)^2, 16*(-3 + d)) + PR4*rat(-(-5 + 2*d)^2, 120*(-3 + d)) + int1d*m9denJ0011X10110(1)*rat(-8 + 3*d, 24);
        id int1d*m9denJ0031X10110(0) = PR4d;
        id int1d*m9denJ1000X01100(0) = 0;
        id int1d*m9denJ1000X01100(1) = int1d*m9denJ0001X00110(1);
        id int1d*m9denJ1001X01110(0) = PR2;
        id int1d*m9denJ1001X11110(0) = PR5;
        id int1d*m9denJ1001X01110(1) = int1d*m9denJ0001X11110(1);
        id int1d*m9denJ1010X01100(0) = PR1; 
        id int1d*m9denJ1010X11100(0) = PR4;
        id int1d*m9denJ1010X01100(1) = int1d*m9denJ0011X00110(1);
        id int1d*m9denJ1010X11100(1) = int1d*m9denJ0011X10110(1);
        id int1d*m9denJ1011X01110(0) = PR6;
        id int1d*m9denJ1011X11110(0) = PR9;
        id int1d*m9denJ2010X01100(1) = PR1*rat(-(-2 + d)^2, 16*(-3 + d)) + int1d*m9denJ0001X00110(1)*rat((-2 + d)^2, 16*(-3 + d)) + int1d*m9denJ0011X00110(1)*rat(-8 + 3*d, 24);
        id int1d*m9denJ2010X11100(1) = PR4d*rat(-1, 24*(-3 + d)) + int1d*m9denJ0001X01110(1)*rat((-2 + d)^2, 16*(-3 + d)) + PR4*rat(-(-5 + 2*d)^2, 120*(-3 + d)) + int1d*m9denJ0011X10110(1)*rat(-8 + 3*d, 24);
        id int1d*m9denJ2011X11110(0) = PR9x;

*       Substituting reduction rules for integrals with [pp-4*tmm]
        id int1d*m4denJ0001X01110(0) = PR1;
        id int1d*m4denJ0001X11110(0) = PR2;
        id int1d*m4denJ0001X01110(1) = PR1*rat(2 - d, 2*(-3 + d)) + int1d*m4denJ0001X00110(1)*rat(-2 + d, 2*(-3 + d));
        id int1d*m4denJ0001X11110(1) = PR2*rat(8 - 3*d, 8*(-3 + d)) + PR1*rat(-(-2 + d)^2, 4*(-3 + d)^2) + int1d*m4denJ0001X00110(1)*rat((-2 + d)^2, 4*(-3 + d)^2);
        id int1d*m4denJ0011X01100(0) = 0;
        id int1d*m4denJ0011X01110(0) = PR3;
        id int1d*m4denJ0011X10110(0) = PR4;
        id int1d*m4denJ0011X11100(0) = PR3;
        id int1d*m4denJ0011X11110(0) = PR7;
        id int1d*m4denJ0011X01110(1) = int1d*m4denJ0021X00110(1)*rat(5, 2*(-3 + d)) + int1d*m4denJ0011X00110(1)*rat(8 - 3*d, 4*(-3 + d)) + PR3*rat(2 - d, 2*(-3 + d)) + int1d*m4denJ0011X01100(1)*rat(-2 + d, 4*(-3 + d));
        id int1d*m4denJ0011X10110(1) = PR4*rat(5 - 2*d, 5*(-3 + d)) + int1d*m4denJ0011X00110(1)*rat(-2 + d, 2*(-3 + d));
        id int1d*m4denJ0011X11100(1) = PR3*rat(2 - d, 2*(-3 + d)) + int1d*m4denJ0011X01100(1)*rat(-2 + d, 2*(-3 + d));
        id int1d*m4denJ0011X11110(1) = PR4d*rat(5, 8*(-3 + d)^2) + PR7*rat(8 - 3*d, 8*(-3 + d)) + int1d*m4denJ0021X00110(1)*rat(5*(-2 + d), 4*(-3 + d)^2) + PR3*rat(-3*(-2 + d)^2, 16*(-3 + d)^2) + int1d*m4denJ0011X01100(1)*rat((-2 + d)^2, 8*(-3 + d)^2) + PR4*rat(-((-7 + 2*d)*(-5 + 2*d)), 40*(-3 + d)^2) + int1d*m4denJ0011X00110(1)*rat(-((-2 + d)*(-8 + 3*d)), 8*(-3 + d)^2);
        id int1d*m4denJ0021X10110(1) = PR4d*rat(1, 2*(-3 + d)) + PR4*rat(5 - 2*d, 10) + int1d*m4denJ0021X00110(1)*rat(-2 + d, 2*(-3 + d));
        id int1d*m4denJ0031X10110(0) = PR4d;
        id int1d*m4denJ1000X01100(0) = 0;
        id int1d*m4denJ1000X01100(1) = int1d*m4denJ0001X00110(1);
        id int1d*m4denJ1001X01110(0) = PR2;
        id int1d*m4denJ1001X11110(0) = PR5;
        id int1d*m4denJ1001X01110(1) = PR2*rat(8 - 3*d, 8*(-3 + d)) + PR1*rat(-(-2 + d)^2, 4*(-3 + d)^2) + int1d*m4denJ0001X00110(1)*rat((-2 + d)^2, 4*(-3 + d)^2);
        id int1d*m4denJ1001X11110(1) = -PR5/3 + PR1*rat(-(-2 + d)^3, 8*(-3 + d)^3) + int1d*m4denJ0001X00110(1)*rat((-2 + d)^3, 8*(-3 + d)^3) + PR2*rat(-((-2 + d)*(-8 + 3*d)), 16*(-3 + d)^2);
        id int1d*m4denJ1010X01100(0) = PR1;
        id int1d*m4denJ1010X11100(0) = PR4;
        id int1d*m4denJ1010X01100(1) = int1d*m4denJ0011X00110(1);
        id int1d*m4denJ1010X11100(1) = PR4*rat(5 - 2*d, 5*(-3 + d)) + int1d*m4denJ0011X00110(1)*rat(-2 + d, 2*(-3 + d));
        id int1d*m4denJ1011X01110(0) = PR6;
        id int1d*m4denJ1011X11110(0) = PR9;
        id int1d*m4denJ1011X11110(1) = -PR5/6 + PR4d*rat(-5, 16*(-3 + d)^2) + PR9x*rat(3, 2*(-3 + d)) + PR9*rat(10 - 3*d, 4*(-3 + d)) + int1d*m4denJ1011X01110(1)*rat(-2 + d, 2*(-3 + d)) + PR3*rat(-(-2 + d)^2, 32*(-3 + d)^2) + PR4*rat((-7 + 2*d)*(-5 + 2*d), 80*(-3 + d)^2) + PR7*rat(-8 + 3*d, 16*(-3 + d));
        id int1d*m4denJ2010X01100(1) = int1d*m4denJ0021X00110(1);
        id int1d*m4denJ2010X11100(1) = PR4d*rat(1, 2*(-3 + d)) + PR4*rat(5 - 2*d, 10) + int1d*m4denJ0021X00110(1)*rat(-2 + d, 2*(-3 + d));
        id int1d*m4denJ2011X11110(0) = PR9x;

*       Substituting reduction rules for integrals with [pp-3*tmm]
        id int1d*m3denJ0001X01110(0) = PR1;
        id int1d*m3denJ0001X11110(0) = PR2;
        id int1d*m3denJ0011X01100(0) = 0;
        id int1d*m3denJ0011X01110(0) = PR3;
        id int1d*m3denJ0011X10110(0) = PR4;
        id int1d*m3denJ0011X11100(0) = PR3;
        id int1d*m3denJ0011X11110(0) = PR7;
        id int1d*m3denJ0031X10110(0) = PR4d;
        id int1d*m3denJ1000X01100(0) = 0;
        id int1d*m3denJ1000X01100(1) = int1d*m3denJ0001X00110(1);
        id int1d*m3denJ1001X01110(0) = PR2;
        id int1d*m3denJ1001X11110(0) = PR5;
        id int1d*m3denJ1001X01110(1) = int1d*m3denJ0001X11110(1);
        id int1d*m3denJ1010X01100(0) = PR1;
        id int1d*m3denJ1010X11100(0) = PR4;
        id int1d*m3denJ1010X01100(1) = int1d*m3denJ0011X00110(1);
        id int1d*m3denJ1010X11100(1) = int1d*m3denJ0011X10110(1);
        id int1d*m3denJ1011X01110(0) = PR6;
        id int1d*m3denJ1011X11110(0) = PR9;
        id int1d*m3denJ1011X01110(1) = int1d*m3denJ0021X00110(1)*rat(-16, 3*(-4 + d)) + PR2*rat(8 - 3*d, 3*(-4 + d)) + int1d*m3denJ0001X11110(1)*rat(-2*(-3 + d), -4 + d) + PR6*rat(-2*(-3 + d), 3*(-4 + d)) + int1d*m3denJ0011X01110(1)*rat(8*(-3 + d), 3*(-4 + d)) + int1d*m3denJ0011X01100(1)*rat(-2*(-2 + d), 3*(-4 + d)) + int1d*m3denJ0001X01110(1)*rat(2*(-2 + d), 3*(-4 + d)) + PR3*rat(4*(-2 + d), 3*(-4 + d)) + int1d*m3denJ0011X00110(1)*rat(2*(-8 + 3*d), 3*(-4 + d));
        id int1d*m3denJ1011X11110(1) = int1d*m3denJ0021X10110(1)*rat(-16, 3*(-4 + d)) + PR9x*rat(-2, -4 + d) + PR4d*rat(25, 12*(-4 + d)*(-3 + d)) + int1d*m3denJ1001X11110(1)*rat(-2*(-3 + d), -4 + d) + PR5*rat(-2*(-3 + d), 3*(-4 + d)) + int1d*m3denJ0011X11110(1)*rat(8*(-3 + d), 3*(-4 + d)) + int1d*m3denJ0011X11100(1)*rat(-2*(-2 + d), 3*(-4 + d)) + int1d*m3denJ0001X11110(1)*rat(2*(-2 + d), 3*(-4 + d)) + PR3*rat(5*(-2 + d)^2, 24*(-4 + d)*(-3 + d)) + PR4*rat(-((-7 + 2*d)*(-5 + 2*d)), 12*(-4 + d)*(-3 + d)) + PR7*rat(-8 + 3*d, 4*(-4 + d)) + int1d*m3denJ0011X10110(1)*rat(2*(-8 + 3*d), 3*(-4 + d));
        id int1d*m3denJ2010X01100(1) = int1d*m3denJ0021X00110(1);
        id int1d*m3denJ2010X11100(1) = int1d*m3denJ0021X10110(1);
        id int1d*m3denJ2011X11110(0) = PR9x;

*       Substituting reduction rules for integrals with [pp-tmm]
        id int1d*m1denJ0001X01110(0) = PR1;
        id int1d*m1denJ0001X11110(0) = PR2;
        id int1d*m1denJ0001X00110(1) = PR1;
        id int1d*m1denJ0001X01110(1) = PR3;
        id int1d*m1denJ0001X11110(1) = PR6;
        id int1d*m1denJ0011X01100(0) = 0;
        id int1d*m1denJ0011X01110(0) = PR3;
        id int1d*m1denJ0011X10110(0) = PR4;
        id int1d*m1denJ0011X11100(0) = PR3;
        id int1d*m1denJ0011X11110(0) = PR7;
        id int1d*m1denJ0011X00110(1) = PR2;
        id int1d*m1denJ0011X01100(1) = PR3;
        id int1d*m1denJ0011X01110(1) = PR6;
        id int1d*m1denJ0011X10110(1) = PR7;
        id int1d*m1denJ0011X11100(1) = PR8;
        id int1d*m1denJ0011X11110(1) = PR10;
        id int1d*m1denJ0021X00110(1) = PR2*rat(-8 + 3*d, 8);
        id int1d*m1denJ0021X10110(1) = PR4d*rat(-5, 8*(-3 + d)) + PR3*rat(-(-2 + d)^2, 16*(-3 + d)) + PR4*rat((-7 + 2*d)*(-5 + 2*d), 40*(-3 + d)) + PR7*rat(-8 + 3*d, 8);
        id int1d*m1denJ0031X10110(0) = PR4d;
        id int1d*m1denJ1000X01100(0) = 0;
        id int1d*m1denJ1000X01100(1) = PR1;
        id int1d*m1denJ1001X01110(0) = PR2;
        id int1d*m1denJ1001X11110(0) = PR5;
        id int1d*m1denJ1001X01110(1) = PR6;
        id int1d*m1denJ1001X11110(1) = PR14;
        id int1d*m1denJ1010X01100(0) = PR1;
        id int1d*m1denJ1010X11100(0) = PR4;
        id int1d*m1denJ1010X01100(1) = PR2;
        id int1d*m1denJ1010X11100(1) = PR7;
        id int1d*m1denJ1011X01110(0) = PR6;
        id int1d*m1denJ1011X11110(0) = PR9;
        id int1d*m1denJ1011X01110(1) = PR13;
        id int1d*m1denJ1011X11110(1) = PR15;
        id int1d*m1denJ2010X01100(1) = PR2*rat(-8 + 3*d, 8);
        id int1d*m1denJ2010X11100(1) = PR4d*rat(-5, 8*(-3 + d)) + PR3*rat(-(-2 + d)^2, 16*(-3 + d)) + PR4*rat((-7 + 2*d)*(-5 + 2*d), 40*(-3 + d)) + PR7*rat(-8 + 3*d, 8);
        id int1d*m1denJ2011X11110(0) = PR9x;

*       Substituting reduction rules for integrals with (pp)
        id int1d*m0denJ0001X01110(0) = PR1;
        id int1d*m0denJ0001X11110(0) = PR2;
        id int1d*m0denJ0001X00110(1) = 0;
        id int1d*m0denJ0001X01110(1) = PR1*rat(-2 + d, 2*(-3 + d));
        id int1d*m0denJ0001X11110(1) = PR1*rat(-(-2 + d)^2, 2*(-4 + d)*(-3 + d)) + PR2*rat(-8 + 3*d, 4*(-4 + d));
        id int1d*m0denJ0011X01100(0) = 0;
        id int1d*m0denJ0011X01110(0) = PR3;
        id int1d*m0denJ0011X10110(0) = PR4;
        id int1d*m0denJ0011X11100(0) = PR3;
        id int1d*m0denJ0011X11110(0) = PR7;
        id int1d*m0denJ0011X01100(1) = 0;
        id int1d*m0denJ0011X01110(1) = int1d*m0denJ0011X00110(1)*rat(8 - 3*d, 6*(-3 + d)) + PR3*rat(-2 + d, 2*(-3 + d));
        id int1d*m0denJ0011X10110(1) = PR4d*rat(-5, -3 + d) + int1d*m0denJ0011X00110(1)*rat(-2 + d, 2) + PR4*rat((-4 + d)*(-5 + 2*d), 5*(-3 + d));
        id int1d*m0denJ0011X11100(1) = PR3*rat(-2 + d, 2*(-3 + d));
        id int1d*m0denJ0011X11110(1) = PR3*rat(-3*(-2 + d)^2, 8*(-4 + d)*(-3 + d)) + PR4d*rat(5*(-7 + 2*d), 4*(-4 + d)*(-3 + d)) + PR7*rat(-8 + 3*d, 4*(-4 + d)) + int1d*m0denJ0011X00110(1)*rat(-((-2 + d)*(-8 + 3*d)), 12*(-3 + d)) + PR4*rat(-((-5 + 2*d)*(99 - 50*d + 6*d^2)), 60*(-4 + d)*(-3 + d));
        id int1d*m0denJ0021X00110(1) = int1d*m0denJ0011X00110(1)*rat(-8 + 3*d, 6);
        id int1d*m0denJ0021X10110(1) = (-5*PR4d)/2 + PR4*rat((-5 + 2*d)*(-13 + 3*d), 30) + int1d*m0denJ0011X00110(1)*rat((-2 + d)*(-8 + 3*d), 12);
        id int1d*m0denJ0031X10110(0) = PR4d;
        id int1d*m0denJ1000X01100(0) = 0;
        id int1d*m0denJ1000X01100(1) = 0;
        id int1d*m0denJ1001X01110(0) = PR2;
        id int1d*m0denJ1001X11110(0) = PR5;
        id int1d*m0denJ1001X01110(1) = PR1*rat(-(-2 + d)^2, 2*(-4 + d)*(-3 + d)) + PR2*rat(-8 + 3*d, 4*(-4 + d));
        id int1d*m0denJ1001X11110(1) = PR5*rat(-3 + d, -5 + d) + PR1*rat(3*(-2 + d)^3, 4*(-5 + d)*(-4 + d)*(-3 + d)) + PR2*rat(-3*(-2 + d)*(-8 + 3*d), 8*(-5 + d)*(-4 + d));
        id int1d*m0denJ1010X01100(0) = PR1;
        id int1d*m0denJ1010X11100(0) = PR4;
        id int1d*m0denJ1010X01100(1) = int1d*m0denJ0011X00110(1);
        id int1d*m0denJ1010X11100(1) = PR4d*rat(-5, -3 + d) + int1d*m0denJ0011X00110(1)*rat(-2 + d, 2) + PR4*rat((-4 + d)*(-5 + 2*d), 5*(-3 + d));
        id int1d*m0denJ1011X01110(0) = PR6;
        id int1d*m0denJ1011X11110(0) = PR9;
        id int1d*m0denJ1011X01110(1) = PR6*rat(2*(-3 + d), 3*(-4 + d)) + PR3*rat(-2*(-2 + d), 3*(-4 + d)) + PR1*rat((-2 + d)^2, 6*(-4 + d)*(-3 + d)) + PR2*rat(-8 + 3*d, 12*(-4 + d));
        id int1d*m0denJ1011X11110(1) = PR9x*rat(1, 2*(-5 + d)) + PR5*rat(-3 + d, 6*(-5 + d)) + PR6*rat(-((-3 + d)*(-2 + d)), 3*(-5 + d)*(-4 + d)) + PR1*rat(-(-2 + d)^3, 4*(-5 + d)*(-4 + d)*(-3 + d)) + PR9*rat(-10 + 3*d, 4*(-5 + d)) + PR2*rat((-2 + d)*(-8 + 3*d), 24*(-5 + d)*(-4 + d)) + PR4d*rat(-5*(-12 + 5*d), 48*(-5 + d)*(-4 + d)*(-3 + d)) + PR7*rat(-((-8 + 3*d)*(-36 + 11*d)), 48*(-5 + d)*(-4 + d)) + PR3*rat((-2 + d)^2*(-180 + 59*d), 96*(-5 + d)*(-4 + d)*(-3 + d)) + PR4*rat(-((-5 + 2*d)*(132 - 95*d + 18*d^2)), 720*(-5 + d)*(-4 + d)*(-3 + d));
        id int1d*m0denJ2010X01100(1) = int1d*m0denJ0011X00110(1)*rat(-8 + 3*d, 6);
        id int1d*m0denJ2010X11100(1) = (-5*PR4d)/2 + PR4*rat((-5 + 2*d)*(-13 + 3*d), 30) + int1d*m0denJ0011X00110(1)*rat((-2 + d)*(-8 + 3*d), 12);
        id int1d*m0denJ2011X11110(0) = PR9x;

*       Integrals with numerators from n5        
        id int1d*numJ0001X00110(0) = 0;
        id int1d*numJ0001X01110(0) = PR1;
        id int1d*numJ0001X11110(0) = PR2;
        id int1d*numJ0011X00110(0) = PR1;
        id int1d*numJ0011X01100(0) = 0;
        id int1d*numJ0011X01110(0) = PR3;
        id int1d*numJ0011X10110(0) = PR4;
        id int1d*numJ0011X11100(0) = PR3;
        id int1d*numJ0011X11110(0) = PR7;
        id int1d*numJ0021X00110(0) = PR1*rat(-2 + d, 2);
        id int1d*numJ0021X10110(0) = PR4*rat(-5 + 2*d, 5);
        id int1d*numJ0031X10110(0) = PR4d;
        id int1d*numJ1000X01100(0) = 0;
        id int1d*numJ1001X01110(0) = PR2;
        id int1d*numJ1001X11110(0) = PR5;
        id int1d*numJ1010X01100(0) = PR1;
        id int1d*numJ1010X11100(0) = PR4;
        id int1d*numJ1011X01110(0) = PR6;
        id int1d*numJ1011X11110(0) = PR9;
        id int1d*numJ2010X01100(0) = PR1*rat(-2 + d, 2);
        id int1d*numJ2010X11100(0) = PR4*rat(-5 + 2*d, 5);
        id int1d*numJ2011X11110(0) = PR9x;

*         
*       9x = [0002111111] -> 9d = [0110211110]
*         
        id PR9x = PR4d*rat(-5, 8*(-3 + d)) + PR9d*rat(-1, 1) + 
        PR5*rat(3 - d, 3) + PR9*rat(-4 + d, 2) + 
        PR3*rat(-(-2 + d)^2, 16*(-3 + d)) + 
        PR4*rat((-7 + 2*d)*(-5 + 2*d), 40*(-3 + d)) + PR7*rat(-8 + 3*d, 8);
        
        .sort:sub-int1d;
*         
*       Now we check that all terms with m^2 != tmm in denominator are eliminated       
*         
        if(count(int1d,1)) Print "Uneliminated integral with aux mass in term %t";
        
        
#endprocedure
*--#] fmft :


*--#[ submi :
#procedure submi
*         
* Expansion of integrals from  
* hep-ph/0411261 Nucl.Phys. B710 (2005) 485-498.        
*         
   id PR0 =
       + PR0ep0

       + ep * (
          Oep(1,PR0)
          )
        ;

   id PR12 =
       + PR12ep0

       + ep * (
          Oep(1,PR12)
          )
        ;

   id PR15 =

       + ep^-2 * (
          + 3/2*z3
          )

       + ep^-1 * (
          - 3/4*z4
          + D6
          + 3/2*z3
          )

       + PR15ep0

       + ep * (
          Oep(1,PR15)
          )
        ;

   id PR11 =

       + ep^-1 * (
          + 5*z5
          )

       + PR11ep0

       + ep * (
          Oep(1,PR11)
          )
        ;

   id PR11d =
       + PR11dep0

       + ep * (
          Oep(1,PR11d)
          )
        ;

   id PR14 =

       + ep^-4 * (
          + 3/4
          )

       + ep^-3 * (
          + 25/4
          )

       + ep^-2 * (
          + 137/4
          + 3/2*z2
          - 81/2*S2
          )

       + ep^-1 * (
          + 363/4
          - 3*T1ep
          - z2
          - 162*S2
          - 27/2*z3
          )

       + PR14ep0

       + ep * (
          Oep(1,PR14)
          )
        ;

   id PR13 =

       + ep^-2 * (
          + 2*z3
          )

       + ep^-1 * (
          + D6
          + 2*z3
          )

       + PR13ep0

       + ep * (
          Oep(1,PR13)
          )
        ;

   id PR10 =

       + ep^-4 * (
          + 1/2
          )

       + ep^-3 * (
          + 53/12
          )

       + ep^-2 * (
          + 51/2
          + z2
          - 27*S2
          )

       + ep^-1 * (
          + 937/12
          - 2*T1ep
          - 1/6*z2
          - 243/2*S2
          - 32/3*z3
          )

       + PR10ep0

       + ep * (
          Oep(1,PR10)
          )
        ;

   id PR9 =

       + ep^-4 * (
          + 1/4
          )

       + ep^-3 * (
          + 7/3
          )

       + ep^-2 * (
          + 169/12
          + 1/2*z2
          - 27/2*S2
          + z3
          )

       + ep^-1 * (
          + 143/3
          - T1ep
          + 1/6*z2
          - 135/2*S2
          + 3/2*z4
          - 4/3*z3
          )

       + PR9ep0

       + ep * (
          + PR9ep1
          )

       + ep^2 * (
          Oep(2,PR9)
          )
        ;

   id PR9d =

       + ep^-4 * (
          + 1/12
          )

       + ep^-3 * (
          + 1/3
          )

       + ep^-2 * (
          + 7/12
          + 1/6*z2
          - 9/2*S2
          )

       + ep^-1 * (
          - 26/3
          - 1/3*T1ep
          - 5/6*z2
          + 27/2*S2
          + 29/9*z3
          )

       + PR9dep0

       + ep * (
          + PR9dep1
          )

       + ep^2 * (
          Oep(2,PR9d)
          )
        ;

   id PR8 =

       + ep^-4 * (
          + 9/4
          )

       + ep^-3 * (
          + 27/2
          )

       + ep^-2 * (
          + 207/4
          + 9/2*z2
          - 81/2*S2
          )

       + ep^-1 * (
          + 189/2
          - 3*T1ep
          + 27/2*z2
          - 243/2*S2
          )

       + PR8ep0

       + ep * (
          Oep(1,PR8)
          )
        ;

   id PR7 =

       + ep^-4 * (
          + 13/8
          )

       + ep^-3 * (
          + 491/48
          )

       + ep^-2 * (
          + 3719/96
          + 13/4*z2
          - 81/4*S2
          )

       + ep^-1 * (
          + 13741/192
          - 3/2*T1ep
          + 329/24*z2
          - 459/8*S2
          - 2/3*z3
          )

       + 381313/10368
          - 3/2*T1ep2
          - 1/4*PR4ep0
          - 5/4*PR4dep0
          - 17/4*T1ep
          + 6805/144*z2
          - 1593/16*S2
          - 81/4*S2*z2
          + 153/16*z4
          + 61/2*z3

       + ep * (
          + PR7ep1
          )

       + ep^2 * (
          Oep(2,PR7)
          )
        ;

   id PR6 =

       + ep^-4 * (
          + 1
          )

       + ep^-3 * (
          + 20/3
          )

       + ep^-2 * (
          + 29
          + 2*z2
          - 27*S2
          )

       + ep^-1 * (
          + 181/3
          - 2*T1ep
          + 13/3*z2
          - 81*S2
          - 16/3*z3
          )

       + PR6ep0

       + ep * (
          Oep(1,PR6)
          )
        ;

   id PR5 =

       + ep^-4 * (
          + 3/2
          )

       + ep^-3 * (
          + 19/2
          )

       + ep^-2 * (
          + 67/2
          + 3*z2
          )

       + ep^-1 * (
          + 127/2
          + 19*z2
          - 5*z3
          )

       + PR5ep0

       + ep * (
          Oep(1,PR5)
          )
        ;

   id PR4 =

       + ep^-4 * (
          + 5/2
          )

       + ep^-3 * (
          + 35/3
          )

       + ep^-2 * (
          + 4565/144
          + 5*z2
          )

       + ep^-1 * (
          + 58345/864
          + 70/3*z2
          - 10/3*z3
          )

       + PR4ep0

       + ep * (
          + PR4ep1
          )

       + ep^2 * (
          Oep(2,PR4)
          )
        ;

   id PR4d =

       + ep^-3 * (
          - 7/6
          )

       + ep^-2 * (
          - 215/48
          )

       + ep^-1 * (
          - 965/96
          - 7/3*z2
          )

       + PR4dep0

       + ep * (
          + PR4dep1
          )

       + ep^2 * (
          Oep(2,PR4d)
          )
        ;

   id PR3 =

       + ep^-4 * (
          + 3/2
          )

       + ep^-3 * (
          + 15/2
          )

       + ep^-2 * (
          + 24
          + 3*z2
          - 27/2*S2
          )

       + ep^-1 * (
          + 81/2
          - T1ep
          + 21/2*z2
          - 27*S2
          - z3
          )

       + 57
          - T1ep2
          - 2*T1ep
          + 57/2*z2
          - 81/2*S2
          - 27/2*S2*z2
          + 51/8*z4
          - 5*z3

       + ep * (
          + PR3ep1
          )

       + ep^2 * (
          Oep(2,PR3)
          )
        ;

   id PR2 =

       + ep^-4 * (
          + 2
          )

       + ep^-3 * (
          + 29/3
          )

       + ep^-2 * (
          + 163/6
          + 4*z2
          )

       + ep^-1 * (
          + 601/12
          + 58/3*z2
          - 8/3*z3
          )

       + 635/24
          + 163/3*z2
          + 12*z4
          + 220/9*z3

       + ep * (
          + PR2ep1
          )

       + ep^2 * (
          Oep(2,PR2)
          )
        ;

   id PR1 =

       + ep^-4 * (
          + 1
          )

       + ep^-3 * (
          + 4
          )

       + ep^-2 * (
          + 10
          + 2*z2
          )

       + ep^-1 * (
          + 20
          + 8*z2
          - 4/3*z3
          )

       + 35
          + 20*z2
          + 6*z4
          - 16/3*z3

       + ep * (
          + PR1ep1
          )

       + ep^2 * (
          Oep(2,PR1)
          )
        ;

        
#endprocedure        
*--#] submi :


*--#[ exp4d :
#procedure exp4d(maxeppow)
* <1>   Expansion near d=4-2*ep      
        Multiply replace_(d,4-2*ep);
        
* <2>   substituting master integrals        
        #call submi
*
*	Expands the PolyRatFun to sufficient powers in ep.
*
.sort:expansion-1;
PolyRatFun;
id	rat(x1?,x2?) = num(x1)*den(x2);
SplitArg,den;
id	den(?a,x1?) = den(x1,?a);
repeat id den(x1?,x2?,x3?,?a) = den(x1,x2+x3,?a);
id	den(x1?,x2?) = den(1,x2/x1)/x1;
id	den(x1?) = 1/x1;
.sort:expansion-2;
id	num(x1?) = x1;
if ( count(ep,1) > `maxeppow' ) discard;
repeat;
	id den(1,x?) = 1-x*den(1,x);
	if ( count(ep,1) > `maxeppow' ) discard;
endrepeat;
.sort:expansion-3;
Symbol ep(:`maxeppow');
        .sort
* Switch back to ep without truncation        
S ep;        
#endprocedure        
*--#] exp4d :

*end-file
