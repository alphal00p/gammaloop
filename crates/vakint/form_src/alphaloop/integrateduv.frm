Auto S cMi1L1, cmi1L2, cmi2L2, alarm;
S eulergamma, log4pi, pi, cl2, sqrt3, MASTER2m22EP, MASTER3m111111FINPIECE, MASTER3m011111FIN, MASTER3m011101EP, MASTER3m001111FIN;
CF uvid;

#if `MAXPOLE' == 3
    #define SPURIOUSPOLE "1"
#else
    #define SPURIOUSPOLE "0"
#endif

#procedure TruncateExpansion(SPURIOUSPOWER)
    argument rat;
        id ep^{`MAXPOLE'+`SELECTEDEPSILONORDER'+`SPURIOUSPOWER'+1} = 0;
    endargument;
#endprocedure

#procedure IntegrateUV1L()
    if ((count(vxs,1) == 1) && match(vxs(k1?vector_,-k1?)));
* reduce the numerator
        repeat id g(k1?,k1?)*uvprop(k1?,n1?) = uvprop(k1, n1-1) + mUV^2 * uvprop(k1, n1);
        id uvprop(k1?,n?) = uvprop(n);

* 1-loop IBP
        id uvprop(n1?{<1}) = 0;
        repeat id uvprop(n1?{>1}) = uvprop(-1 + n1)*rat((2 + (4-2*ep) - 2*n1), 2* (-1 + n1)) / mUV^2;
        id uvprop(1) = uvid(1,1) * rat(1,ep);
    endif;
    id vxs(k1?,-k1?) = 1;
#endprocedure

#procedure IntegrateUV2L()
    id uvprop(k1?,n1?) = uvprop(k1,n1)*tmps(k1,-k1);

    id vxs(k1?,kn1?,k2?,kn2?)*tmps(k1?,kn1?)*tmps(k2?,kn2?) = uvid(2,1,0,0,0)*map(kk1,k1,1)*map(kk1,kn1,-1)*map(kk2,k2,1)*map(kk2,kn2,-1);
    id vxs(k1?,k2?,k3?)*vxs(kn1?,kn2?,kn3?)*tmps(k1?,kn1?)*tmps(k2?,kn2?)*tmps(k3?,kn3?) =
        uvid(2,1,0,0,0)*map(kk1,k1,1)*map(kk1,kn1,-1)*map(kk2,k2,1)*map(kk2,kn2,-1)*map(kk3,k3,1)*map(kk3,kn3,-1);

    id tmps(k1?,k2?) = 1;

    id map(kk1,k?,nn?)*uvprop(k?,n?)*uvid(2,1,n1?,n2?,n3?) = map(kk1,k,nn)*uvid(2,1,n1+n,n2,n3);
    id map(kk2,k?,nn?)*uvprop(k?,n?)*uvid(2,1,n1?,n2?,n3?) = map(kk2,k,nn)*uvid(2,1,n1,n2+n,n3);
    id map(kk3,k?,nn?)*uvprop(k?,n?)*uvid(2,1,n1?,n2?,n3?) = map(kk3,k,nn)*uvid(2,1,n1,n2,n3+n);

* reduce the numerator
    repeat id g(k1?,k1?)*map(kk1,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = (uvid(2,1,n1-1,n2,n3) + mUV^2 * uvid(2,1,n1,n2,n3))*map(kk1,k1,n);
    repeat id g(k1?,k1?)*map(kk2,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = (uvid(2,1,n1,n2-1,n3) + mUV^2 * uvid(2,1,n1,n2,n3))*map(kk2,k1,n);
    repeat id g(k1?,k1?)*map(kk3,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = (uvid(2,1,n1,n2,n3-1) + mUV^2 * uvid(2,1,n1,n2,n3))*map(kk3,k1,n);

    repeat id g(k1?,k2?)*map(kk1,k1?,nn1?)*map(kk2,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*map(kk1,k1,nn1)*map(kk2,k2,nn2)*(
        uvid(2,1,n1,n2,n3-1) - uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2-1,n3) - mUV^2 * uvid(2,1,n1,n2,n3)
    );

    repeat id g(k1?,k2?)*map(kk1,k1?,nn1?)*map(kk3,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*map(kk1,k1,nn1)*map(kk3,k2,nn2)*(
        uvid(2,1,n1,n2-1,n3) - uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2,n3-1) - mUV^2 * uvid(2,1,n1,n2,n3)
    );

    repeat id g(k1?,k2?)*map(kk2,k1?,nn1?)*map(kk3,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*map(kk2,k1,nn1)*map(kk3,k2,nn2)*(
        uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2-1,n3) - uvid(2,1,n1,n2,n3-1) - mUV^2 * uvid(2,1,n1,n2,n3)
    );

    id map(?a) = 1;

    #define nloop "1"
    #do loop=1,1
* zero sectors
        id uvid(2,1,n1?{<1},n2?{<1},n3?{<1}) = 0;
        id uvid(2,1,n1?{<1},n2?{<1},n3?) = 0;
        id uvid(2,1,n1?{<1},n2?,n3?{<1}) = 0;
        id uvid(2,1,n1?,n2?{<1},n3?{<1}) = 0;

* sector mappings
        id uvid(2, 1, n1?{>0}, n2?{<1}, n3?{>0}) = uvid(2, 1, n2, n3, n1);
        id uvid(2, 1, n1?{>0}, n2?{>0}, n3?{<1}) = uvid(2, 1, n3, n2, n1);

* actual reduction
        id ifmatch->end uvid(2, 1, n1?{<0}, n2?{>0}, n3?{>1}) =
            uvid(2, 1, 1 + n1, -1 + n2, n3)*rat(1,1) +
            uvid(2, 1, 1 + n1, n2, -1 + n3)*rat((-3 + (4-2*ep) - 3*n1),2*(-1 + n3)) +
            uvid(2, 1, 2 + n1, -1 + n2, -1 + n3)*rat((-1 - n1),2*(-1 + n3)) +
            uvid(2, 1, 2 + n1, n2, -2 + n3)*rat((1 + n1),2*(-1 + n3)) -
            uvid(2, 1, 2 + n1, n2, -1 + n3)*mUV^2*rat(3*(1 + n1),2*(-1 + n3));

        id ifmatch->end uvid(2, 1, n1?{<0}, n2?{>1}, n3?{>0}) =
            uvid(2, 1, 1 + n1, -1 + n2, n3)*rat((-3 + (4-2*ep) - 3*n1),2*(-1 + n2)) +
            uvid(2, 1, 1 + n1, n2, -1 + n3)*rat(1,1) +
            uvid(2, 1, 2 + n1, -2 + n2, n3)*rat((1 + n1),2*(-1 + n2)) +
            uvid(2, 1, 2 + n1, -1 + n2, -1 + n3)*rat((-1 - n1),2*(-1 + n2)) -
            uvid(2, 1, 2 + n1, -1 + n2, n3)*mUV^2*rat(3*(1 + n1),2*(-1 + n2));

        id ifmatch->end uvid(2, 1, n1?{<1}, n2?{>0}, n3?{>1})  =
            uvid(2, 1, n1, n2, -1 + n3)*mUV^-2*rat(2 + (4-2*ep) - n1 - 2*n3,2*(-1 + n3)) +
            uvid(2, 1, 1 + n1, -1 + n2, -1 + n3)*mUV^-2*rat(n1,2*(-1 + n3)) -
            uvid(2, 1, 1 + n1, n2, -2 + n3)*mUV^2*rat(n1,2*(-1 + n3)) -
            uvid(2, 1, 1 + n1, n2, -1 + n3)*rat(n1,2*(-1 + n3));

        id ifmatch->end uvid(2, 1, n1?{<1}, n2?{>1}, n3?{>0})  =
            uvid(2, 1, n1, -1 + n2, n3)*mUV^-2*rat(2 + (4-2*ep) - n1 - 2*n2,2*(-1 + n2)) -
            uvid(2, 1, 1 + n1, -2 + n2, n3)*mUV^-2*rat(n1,2*(-1 + n2)) +
            uvid(2, 1, 1 + n1, -1 + n2, -1 + n3)*mUV^-2*rat(n1,2*(-1 + n2)) -
            uvid(2, 1, 1 + n1, -1 + n2, n3)*rat(n1,2*(-1 + n2));

        id ifmatch->end uvid(2, 1, n1?{<0}, 1, 1)  =
            uvid(2, 1, 1 + n1, 0, 1)*rat(n1,-2 + (4-2*ep) - n1) -
            uvid(2, 1, 1 + n1, 1, 0)*rat(n1,-2 + (4-2*ep) - n1) +
            uvid(2, 1, 1 + n1, 1, 1)*mUV^2*rat(-3 + (4-2*ep) - 2*n1,-2 + (4-2*ep) - n1) +
            uvid(2, 1, 1 + n1, 2, 0)*mUV^2*rat(2,-2 + (4-2*ep) - n1) -
            uvid(2, 1, 2 + n1, 0, 1)*mUV^2*rat(1 + n1,2 - (4-2*ep) + n1) -
            uvid(2, 1, 2 + n1, 1, 0)*mUV^2*rat(1 + n1,-2 + (4-2*ep) - n1) -
            uvid(2, 1, 2 + n1, 1, 1)*mUV^4*rat(3*(1 + n1),-2 + (4-2*ep) - n1);

        id ifmatch->end uvid(2, 1, n1?{>0}, n2?{>0}, n3?{>1})  =
            uvid(2, 1, -1 + n1, 1 + n2, -1 + n3)*mUV^-2*rat(n2,3*(-1 + n3)) +
            uvid(2, 1, n1, n2, -1 + n3)*mUV^-2*rat((3 + (4-2*ep) - 3*n3),3*(-1 + n3)) -
            uvid(2, 1, n1, 1 + n2, -2 + n3)*mUV^-2*rat(n2,3*(-1 + n3)) +
            uvid(2, 1, 1 + n1, -1 + n2, -1 + n3)*mUV^-2*rat(n1,3*(-1 + n3)) -
            uvid(2, 1, 1 + n1, n2, -2 + n3)*mUV^-2*rat(n1,3*(-1 + n3));

        id ifmatch->end uvid(2, 1, n1?{>0}, n2?{>1}, n3?{>0}) =
            uvid(2, 1, -1 + n1, n2, n3)*mUV^-2*rat(1,3) +
            uvid(2, 1, n1, -1 + n2, n3)*mUV^-2*rat(3 + (4-2*ep) - 3*n2,3*(-1 + n2)) -
            uvid(2, 1, n1, n2, -1 + n3)*mUV^-2*rat(1,3) -
            uvid(2, 1, 1 + n1, -2 + n2, n3)*mUV^-2*rat(2*n1,3*(-1 + n2)) +
            uvid(2, 1, 1 + n1, -1 + n2, -1 + n3)*mUV^-2*rat(2*n1,3*(-1 + n2));

        id ifmatch->end uvid(2, 1, n1?{>1}, n2?{>0}, n3?{>0}) =
            uvid(2, 1, -2 + n1, 1 + n2, n3)*mUV^-2*rat(-2*n2,3*(-1 + n1)) +
            uvid(2, 1, -1 + n1, n2, n3)*mUV^-2*rat((3 + (4-2*ep) - 3*n1),3*(-1 + n1)) +
            uvid(2, 1, -1 + n1, 1 + n2, -1 + n3)*mUV^-2*rat(2*n2,3*(-1 + n1)) +
            uvid(2, 1, n1, -1 + n2, n3)*mUV^-2*rat(1,3) -
            uvid(2, 1, n1, n2, -1 + n3)*mUV^-2*rat(1,3);

        id ifmatch->end uvid(2, 1, n1?{>0}, n2?{>0}, n3?{>1})  =
            uvid(2, 1, n1, 1 + n2, -1 + n3)*rat(n2,-1 + n3) +
            uvid(2, 1, 1 + n1, -1 + n2, n3)*rat(1,1) -
            uvid(2, 1, 1 + n1, 1 + n2, -2 + n3)*rat(n2,-1 + n3) +
            uvid(2, 1, 2 + n1, -1 + n2, -1 + n3)*rat((-1 - n1),-1 + n3) +
            uvid(2, 1, 2 + n1, n2, -2 + n3)*rat((1 + n1),-1 + n3);

        label end;

        id uvid(2,1,0,1,1) = tp2P011;
        id uvid(2,1,1,1,1) = tp2P111;

        if (match(uvid(2, 1, ?a)));
            redefine loop "0";
        endif;

        .sort:reduction-`nloop++';
    #enddo
* MIs are uvid(2,1,0,2,2) and uvid(2,1,2,2,1)
    id tp2P011 = uvid(2,1)*rat(1,ep^2)*mUV^4*rat(4,(4-2*ep)^2 - 4*(4-2*ep) + 4);
    id tp2P111 = mUV^2*rat(9, (4-2*ep)^2 - 5*(4-2*ep) + 6)*(uvid(2,2) - uvid(2,1)*rat(1,ep^2)*rat(-1,3));
#endprocedure

#procedure IntegrateUV3L()
    B+ vxs,uvprop,uvid,g;
    .sort:mapping-3l-start;
    Keep brackets;

    id uvprop(k1?,n1?) = uvprop(k1,n1)*tmps(k1,-k1);

* mercedes
    id ifmatch->matchend3l vxs(kn1?,k2?,k4?)*vxs(kn2?,k3?,k5?)*vxs(kn4?,kn5?,k6?)*vxs(k1?,kn3?,kn6?)*tmps(k1?,kn1?)*tmps(k2?,kn2?)*tmps(k3?,kn3?)*tmps(k4?,kn4?)*tmps(k5?,kn5?)*tmps(k6?,kn6?) =
        uvid(3,1,0,0,0,0,0,0)*map(kk1,k1,1)*map(kk1,kn1,-1)*map(kk2,k2,1)*map(kk2,kn2,-1)*map(kk3,k3,1)*map(kk3,kn3,-1)*map(kk4,k4,1)*map(kk4,kn4,-1)*map(kk5,k5,1)*map(kk5,kn5,-1)*map(kk6,k6,1)*map(kk6,kn6,-1);

* 5-edge
    id ifmatch->matchend3l vxs(kn1?,k2?,k4?)*vxs(kn2?,k3?,k5?)*vxs(k1?,kn3?,kn4?,kn5?)*tmps(k1?,kn1?)*tmps(k2?,kn2?)*tmps(k3?,kn3?)*tmps(k4?,kn4?)*tmps(k5?,kn5?) =
        uvid(3,1,0,0,0,0,0,0)*map(kk1,k1,1)*map(kk1,kn1,-1)*map(kk2,k2,1)*map(kk2,kn2,-1)*map(kk3,k3,1)*map(kk3,kn3,-1)*map(kk4,k4,1)*map(kk4,kn4,-1)*map(kk5,k5,1)*map(kk5,kn5,-1);

* banana
    id ifmatch->matchend3l vxs(kn1?,k3?,k4?,k5?)*vxs(k1?,kn3?,kn4?,kn5?)*tmps(k1?,kn1?)*tmps(k3?,kn3?)*tmps(k4?,kn4?)*tmps(k5?,kn5?) =
        uvid(3,1,0,0,0,0,0,0)*map(kk1,k1,1)*map(kk1,kn1,-1)*map(kk3,k3,1)*map(kk3,kn3,-1)*map(kk4,k4,1)*map(kk4,kn4,-1)*map(kk5,k5,1)*map(kk5,kn5,-1);

* sunrise-bubble
    id ifmatch->matchend3l vxs(kn1?,k2?,k4?)*vxs(k1?,kn2?,kn4?,k3?,kn3?)*tmps(k1?,kn1?)*tmps(k2?,kn2?)*tmps(k3?,kn3?)*tmps(k4?,kn4?) =
        uvid(3,1,0,0,0,0,0,0)*map(kk1,k1,1)*map(kk1,kn1,-1)*map(kk2,k2,1)*map(kk2,kn2,-1)*map(kk3,k3,1)*map(kk3,kn3,-1)*map(kk4,k4,1)*map(kk4,kn4,-1);

* triple-bubble
    id vxs(k1?,k2?,k3?,kn1?,kn2?,kn3?)*tmps(k1?,kn1?)*tmps(k2?,kn2?)*tmps(k3?,kn3?) =
        uvid(3,1,0,0,0,0,0,0)*map(kk1,k1,1)*map(kk1,kn1,-1)*map(kk2,k2,1)*map(kk2,kn2,-1)*map(kk3,k3,1)*map(kk3,kn3,-1);

    label matchend3l;

    id map(kk1,k?,nn?)*uvprop(k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = map(kk1,k,nn)*uvid(3,1,n1+n,n2,n3,n4,n5,n6);
    id map(kk2,k?,nn?)*uvprop(k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = map(kk2,k,nn)*uvid(3,1,n1,n2+n,n3,n4,n5,n6);
    id map(kk3,k?,nn?)*uvprop(k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = map(kk3,k,nn)*uvid(3,1,n1,n2,n3+n,n4,n5,n6);
    id map(kk4,k?,nn?)*uvprop(k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = map(kk4,k,nn)*uvid(3,1,n1,n2,n3,n4+n,n5,n6);
    id map(kk5,k?,nn?)*uvprop(k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = map(kk5,k,nn)*uvid(3,1,n1,n2,n3,n4,n5+n,n6);
    id map(kk6,k?,nn?)*uvprop(k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = map(kk6,k,nn)*uvid(3,1,n1,n2,n3,n4,n5,n6+n);

* reduce the numerator
    #define nloop "1"
    #do i=1,1
        id ifmatch->endmap`i' g(k?,k?)*map(kk1,k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = (uvid(3,1,n1-1,n2,n3,n4,n5,n6) + mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6))*map(kk1,k,n);
        id ifmatch->endmap`i' g(k?,k?)*map(kk2,k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = (uvid(3,1,n1,n2-1,n3,n4,n5,n6) + mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6))*map(kk2,k,n);
        id ifmatch->endmap`i' g(k?,k?)*map(kk3,k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = (uvid(3,1,n1,n2,n3-1,n4,n5,n6) + mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6))*map(kk3,k,n);
        id ifmatch->endmap`i' g(k?,k?)*map(kk4,k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = (uvid(3,1,n1,n2,n3,n4-1,n5,n6) + mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6))*map(kk4,k,n);
        id ifmatch->endmap`i' g(k?,k?)*map(kk5,k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = (uvid(3,1,n1,n2,n3,n4,n5-1,n6) + mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6))*map(kk5,k,n);
        id ifmatch->endmap`i' g(k?,k?)*map(kk6,k?,n?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = (uvid(3,1,n1,n2,n3,n4,n5,n6-1) + mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6))*map(kk6,k,n);

        id ifmatch->endmap`i' g(k1?,k2?)*map(kk1,k1?,nn1?)*map(kk2,k2?,nn2?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn1*nn2*map(kk1,k1,nn1)*map(kk2,k2,nn2)*(
            uvid(3,1,n1,n2,n3,n4-1,n5,n6) - uvid(3,1,n1-1,n2,n3,n4,n5,n6) - uvid(3,1,n1,n2-1,n3,n4,n5,n6) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k1?,k3?)*map(kk1,k1?,nn1?)*map(kk3,k3?,nn3?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn1*nn3*map(kk1,k1,nn1)*map(kk3,k3,nn3)*(
            uvid(3,1,n1,n2,n3,n4,n5,n6-1) - uvid(3,1,n1-1,n2,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3-1,n4,n5,n6) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k1?,k4?)*map(kk1,k1?,nn1?)*map(kk4,k4?,nn4?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn1*nn4*map(kk1,k1,nn1)*map(kk4,k4,nn4)*(
            uvid(3,1,n1,n2-1,n3,n4,n5,n6) - uvid(3,1,n1-1,n2,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3,n4-1,n5,n6) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k1?,k5?)*map(kk1,k1?,nn1?)*map(kk5,k5?,nn5?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn1*nn5*map(kk1,k1,nn1)*map(kk5,k5,nn5)*(
            uvid(3,1,n1,n2,n3,n4-1,n5,n6) + uvid(3,1,n1,n2,n3-1,n4,n5,n6) - uvid(3,1,n1,n2-1,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5,n6-1));
        id ifmatch->endmap`i' g(k1?,k6?)*map(kk1,k1?,nn1?)*map(kk6,k6?,nn6?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn1*nn6*map(kk1,k1,nn1)*map(kk6,k6,nn6)*(
            uvid(3,1,n1,n2,n3-1,n4,n5,n6) - uvid(3,1,n1-1,n2,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5,n6-1) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k2?,k3?)*map(kk2,k2?,nn2?)*map(kk3,k3?,nn3?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn2*nn3*map(kk2,k2,nn2)*map(kk3,k3,nn3)*(
            uvid(3,1,n1,n2,n3,n4,n5-1,n6) - uvid(3,1,n1,n2-1,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3-1,n4,n5,n6) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k2?,k4?)*map(kk2,k2?,nn2?)*map(kk4,k4?,nn4?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = 1/2*nn2*nn4*map(kk2,k2,nn2)*map(kk4,k4,nn4)*(
            uvid(3,1,n1-1,n2,n3,n4,n5,n6) - uvid(3,1,n1,n2-1,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3,n4-1,n5,n6) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k2?,k5?)*map(kk2,k2?,nn2?)*map(kk5,k5?,nn5?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn2*nn5*map(kk2,k2,nn2)*map(kk5,k5,nn5)*(
            uvid(3,1,n1,n2,n3-1,n4,n5,n6) - uvid(3,1,n1,n2-1,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5-1,n6) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k2?,k6?)*map(kk2,k2?,nn2?)*map(kk6,k6?,nn6?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn2*nn6*map(kk2,k2,nn2)*map(kk6,k6,nn6)*(
            uvid(3,1,n1,n2,n3-1,n4,n5,n6) + uvid(3,1,n1,n2,n3,n4-1,n5,n6) - uvid(3,1,n1-1,n2,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5-1,n6));
        id ifmatch->endmap`i' g(k3?,k4?)*map(kk3,k3?,nn3?)*map(kk4,k4?,nn4?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn3*nn4*map(kk3,k3,nn3)*map(kk4,k4,nn4)*(
            uvid(3,1,n1,n2-1,n3,n4,n5,n6) + uvid(3,1,n1,n2,n3,n4,n5,n6-1) - uvid(3,1,n1-1,n2,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5-1,n6));
        id ifmatch->endmap`i' g(k3?,k5?)*map(kk3,k3?,nn3?)*map(kk5,k5?,nn5?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = 1/2*nn3*nn5*map(kk3,k3,nn3)*map(kk5,k5,nn5)*(
            uvid(3,1,n1,n2-1,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3-1,n4,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5-1,n6) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k3?,k6?)*map(kk3,k3?,nn3?)*map(kk6,k6?,nn6?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = 1/2*nn3*nn6*map(kk3,k3,nn3)*map(kk6,k6,nn6)*(
            uvid(3,1,n1-1,n2,n3,n4,n5,n6) - uvid(3,1,n1,n2,n3-1,n4,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5,n6-1) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k4?,k5?)*map(kk4,k4?,nn4?)*map(kk5,k5?,nn5?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = 1/2*nn4*nn5*map(kk4,k4,nn4)*map(kk5,k5,nn5)*(
            uvid(3,1,n1,n2,n3,n4,n5,n6-1) - uvid(3,1,n1,n2,n3,n4-1,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5-1,n6) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k4?,k6?)*map(kk4,k4?,nn4?)*map(kk6,k6?,nn6?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn4*nn6*map(kk4,k4,nn4)*map(kk6,k6,nn6)*(
            uvid(3,1,n1,n2,n3,n4,n5-1,n6) - uvid(3,1,n1,n2,n3,n4-1,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5,n6-1) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));
        id ifmatch->endmap`i' g(k5?,k6?)*map(kk5,k5?,nn5?)*map(kk6,k6?,nn6?)*uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?) = -1/2*nn5*nn6*map(kk5,k5,nn5)*map(kk6,k6,nn6)*(
            uvid(3,1,n1,n2,n3,n4-1,n5,n6) - uvid(3,1,n1,n2,n3,n4,n5-1,n6) - uvid(3,1,n1,n2,n3,n4,n5,n6-1) - mUV^2 * uvid(3,1,n1,n2,n3,n4,n5,n6));

        goto realmapend;
        label endmap`i';
        redefine i "0";
        label realmapend;

        B+ g,map,uvid;
        .sort:num-reduction-3l-`nloop++';
        Keep brackets;
    #enddo

    id map(?a) = 1;

    id tmps(p1?,p2?) = 1;

    B+ uvid;
    .sort:mapping-3l;
    Keep brackets;

    #define nloop "1"
    #do i=1,1
* 0 sectors
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?{<1}, n4?{<1}, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?{<1}, n4?{<1}, n5?{<1}, n6?) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?{<1}, n4?{<1}, n5?, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?{<1}, n4?, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?, n4?{<1}, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?, n3?{<1}, n4?{<1}, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?, n2?{<1}, n3?{<1}, n4?{<1}, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?{<1}, n4?{<1}, n5?, n6?) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?{<1}, n4?, n5?{<1}, n6?) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?{<1}, n4?, n5?, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?, n4?{<1}, n5?{<1}, n6?) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?, n4?{<1}, n5?, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?, n4?, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?, n3?{<1}, n4?{<1}, n5?{<1}, n6?) = 0;
        id uvid(3,1,n1?{<1}, n2?, n3?{<1}, n4?{<1}, n5?, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?, n3?{<1}, n4?, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?, n3?, n4?{<1}, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?, n2?{<1}, n3?{<1}, n4?{<1}, n5?{<1}, n6?) = 0;
        id uvid(3,1,n1?, n2?{<1}, n3?{<1}, n4?{<1}, n5?, n6?{<1}) = 0;
        id uvid(3,1,n1?, n2?{<1}, n3?{<1}, n4?, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?, n2?{<1}, n3?, n4?{<1}, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?, n2?, n3?{<1}, n4?{<1}, n5?{<1}, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?{<1}, n3?{<1}, n4?, n5?, n6?) = 0;
        id uvid(3,1,n1?{<1}, n2?, n3?, n4?{<1}, n5?, n6?{<1}) = 0;
        id uvid(3,1,n1?{<1}, n2?, n3?, n4?{<1}, n5?, n6?{<1}) = 0;
        id uvid(3,1,n1?, n2?{<1}, n3?, n4?{<1}, n5?{<1}, n6?) = 0;
        id uvid(3,1,n1?, n2?, n3?{<1}, n4?, n5?{<1}, n6?{<1}) = 0;

        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{<1}, n5?{<1}, n6?{<1}) = uvid(3,1,n6, n5, n3, n4, n2, n1);
        id uvid(3,1,n1?{<1}, n2?{>0}, n3?{<1}, n4?{>0}, n5?{>0}, n6?{>0}) = uvid(3,1,n3, n1, n2, n6, n4, n5);
        id uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{<1}, n5?{>0}, n6?{>0}) = uvid(3,1,n4, n1, n6, n2, n3, n5);
        id uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{<1}) = uvid(3,1,n6, n1, n4, n3, n2, n5);
        id uvid(3,1,n1?{>0}, n2?{<1}, n3?{<1}, n4?{>0}, n5?{>0}, n6?{>0}) = uvid(3,1,n3, n2, n1, n5, n4, n6);
        id uvid(3,1,n1?{>0}, n2?{<1}, n3?{>0}, n4?{<1}, n5?{>0}, n6?{>0}) = uvid(3,1,n4, n2, n5, n1, n3, n6);
        id uvid(3,1,n1?{>0}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{<1}, n6?{>0}) = uvid(3,1,n5, n2, n4, n3, n1, n6);
        id uvid(3,1,n1?{>0}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{<1}) = uvid(3,1,n6, n5, n3, n4, n2, n1);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{<1}, n4?{<1}, n5?{>0}, n6?{>0}) = uvid(3,1,n4, n6, n1, n5, n3, n2);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{<1}, n4?{>0}, n5?{<1}, n6?{>0}) = uvid(3,1,n5, n3, n6, n2, n1, n4);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{<1}, n4?{>0}, n5?{>0}, n6?{<1}) = uvid(3,1,n6, n3, n5, n1, n2, n4);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{<1}, n5?{<1}, n6?{>0}) = uvid(3,1,n5, n4, n2, n6, n1, n3);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{<1}, n5?{>0}, n6?{<1}) = uvid(3,1,n6, n4, n1, n5, n2, n3);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{<1}, n6?{<1}) = uvid(3,1,n6, n5, n3, n4, n2, n1);
        id uvid(3,1,n1?{>0}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>0}) = uvid(3,1,n2, n5, n4, n3, n6, n1);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{<1}, n4?{>0}, n5?{>0}, n6?{>0}) = uvid(3,1,n3, n6, n5, n1, n4, n2);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{<1}, n5?{>0}, n6?{>0}) = uvid(3,1,n4, n6, n1, n5, n3, n2);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{<1}, n6?{>0}) = uvid(3,1,n5, n6, n3, n4, n1, n2);
        id uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{<1}) = uvid(3,1,n6, n5, n3, n4, n2, n1);

        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{<0}, n5?{>0}, n6?{>1}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3,1 + n4,n5, - 1 + n6)*rat(n1,n6 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,1 + n4,n5, - 1 + n6)*rat( - n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 2 + n6)*rat(n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 1 + n6)*mUV^2*rat( - n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - n1,n6 - 1)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5,n6)*rat(1,1)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 1 + n6)*rat( - 2*ep - n1 - 3*n4 + 1,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,2 + n4, - 1 + n5, - 1 + n6)*rat( - n4 - 1,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,2 + n4,n5, - 2 + n6)*rat(n4 + 1,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,2 + n4,n5, - 1 + n6)*mUV^2*rat( - 3*n4 - 3,2*n6 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{<1}, n3?{>0}, n4?{<1}, n5?{>0}, n6?{>1}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3,1 + n4,n5, - 1 + n6)*rat(n4,n6 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*rat( - n4,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 2 + n6)*rat(n4,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 1 + n6)*mUV^2*rat( - n4,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - 2*ep - 3*n1 - n4 + 1,2*n6 - 2)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - n1 - 1,2*n6 - 2)
            + uvid(3,1,2 + n1,n2,n3,n4,n5, - 2 + n6)*rat(n1 + 1,2*n6 - 2)
            + uvid(3,1,2 + n1,n2,n3,n4,n5, - 1 + n6)*mUV^2*rat( - 3*n1 - 3,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 1 + n6)*rat( - n4,n6 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{<0}, n5?{>1}, n6?{>0}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3,1 + n4, - 1 + n5,n6)*rat(n2,n5 - 1)
            + uvid(3,1,n1,1 + n2, - 1 + n3,1 + n4, - 1 + n5,n6)*rat( - n2,2*n5 - 2)
            + uvid(3,1,n1,1 + n2,n3,1 + n4, - 2 + n5,n6)*rat(n2,2*n5 - 2)
            + uvid(3,1,n1,1 + n2,n3,1 + n4, - 1 + n5,n6)*mUV^2*rat( - n2,2*n5 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4, - 1 + n5,n6)*rat( - n2,n5 - 1)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5,n6)*rat( - 2*ep - n2 - 3*n4 + 1,2*n5 - 2)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 1 + n6)*rat(1,1)
            + uvid(3,1,n1,n2,n3,2 + n4, - 2 + n5,n6)*rat(n4 + 1,2*n5 - 2)
            + uvid(3,1,n1,n2,n3,2 + n4, - 1 + n5, - 1 + n6)*rat( - n4 - 1,2*n5 - 2)
            + uvid(3,1,n1,n2,n3,2 + n4, - 1 + n5,n6)*mUV^2*rat( - 3*n4 - 3,2*n5 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<0}, n3?{>0}, n4?{<1}, n5?{>1}, n6?{>0}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3,1 + n4, - 1 + n5,n6)*rat(n4,n5 - 1)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,n1,1 + n2,n3,1 + n4, - 2 + n5,n6)*rat(n4,2*n5 - 2)
            + uvid(3,1,n1,1 + n2,n3,1 + n4, - 1 + n5, - 1 + n6)*rat( - n4,2*n5 - 2)
            + uvid(3,1,n1,1 + n2,n3,1 + n4, - 1 + n5,n6)*mUV^2*rat( - n4,2*n5 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4, - 1 + n5,n6)*rat( - 2*ep - 3*n2 - n4 + 1,2*n5 - 2)
            + uvid(3,1,n1,2 + n2, - 1 + n3,n4, - 1 + n5,n6)*rat( - n2 - 1,2*n5 - 2)
            + uvid(3,1,n1,2 + n2,n3,n4, - 2 + n5,n6)*rat(n2 + 1,2*n5 - 2)
            + uvid(3,1,n1,2 + n2,n3,n4, - 1 + n5,n6)*mUV^2*rat( - 3*n2 - 3,2*n5 - 2)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5,n6)*rat( - n4,n5 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<0}, n3?{>1}, n4?{<1}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1,1 + n2, - 2 + n3,n4,n5,n6)*rat(n1,2*n3 - 2)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*rat(n1,n3 - 1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - n1,2*n3 - 2)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - n1,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat( - n1,n3 - 1)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5,n6)*rat( - 2*ep - n1 - 3*n2 + 1,2*n3 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4, - 1 + n5,n6)*rat(1,1)
            + uvid(3,1,n1,2 + n2, - 2 + n3,n4,n5,n6)*rat(n2 + 1,2*n3 - 2)
            + uvid(3,1,n1,2 + n2, - 1 + n3,n4, - 1 + n5,n6)*rat( - n2 - 1,2*n3 - 2)
            + uvid(3,1,n1,2 + n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - 3*n2 - 3,2*n3 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{<1}, n3?{>1}, n4?{<1}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1,1 + n2, - 2 + n3,n4,n5,n6)*rat(n2,2*n3 - 2)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*rat(n2,n3 - 1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*rat( - n2,2*n3 - 2)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - n2,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat( - 2*ep - 3*n1 - n2 + 1,2*n3 - 2)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat(1,1)
            + uvid(3,1,2 + n1,n2, - 2 + n3,n4,n5,n6)*rat(n1 + 1,2*n3 - 2)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - n1 - 1,2*n3 - 2)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - 3*n1 - 3,2*n3 - 2)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5,n6)*rat( - n2,n3 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{<1}, n5?{>0}, n6?{>1}) =
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 2 + n6)*mUV^-2*rat( - n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - n1,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(n4,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 2 + n6)*mUV^-2*rat( - n4,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 1 + n6)*rat( - n4,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - 2*ep - n1 - n4 - 2*n6 + 6,2*n6 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{<1}, n5?{>1}, n6?{>0}) =
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(n2,2*n5 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4, - 2 + n5,n6)*mUV^-2*rat( - n2,2*n5 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4, - 1 + n5,n6)*rat( - n2,2*n5 - 2)
            + uvid(3,1,n1,n2,n3,1 + n4, - 2 + n5,n6)*mUV^-2*rat( - n4,2*n5 - 2)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(n4,2*n5 - 2)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5,n6)*rat( - n4,2*n5 - 2)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat( - 2*ep - n2 - n4 - 2*n5 + 6,2*n5 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>1}, n4?{<1}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1,n2, - 2 + n3,n4,n5,n6)*mUV^-2*rat( - n1,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(n1,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat( - n1,2*n3 - 2)
            + uvid(3,1,n1,1 + n2, - 2 + n3,n4,n5,n6)*mUV^-2*rat( - n2,2*n3 - 2)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(n2,2*n3 - 2)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5,n6)*rat( - n2,2*n3 - 2)
            + uvid(3,1,n1,n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat( - 2*ep - n1 - n2 - 2*n3 + 6,2*n3 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, 1, n4?{<0}, 1, 1) =
            + uvid(3,1, - 1 + n1,1 + n2,1,1 + n4,1,1)*mUV^2*rat( - 2*n2,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,1 + n2,0,1 + n4,1,1)*mUV^2*rat(n2,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,1 + n2,0,n4,1,1)*rat(n2,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,1 + n2,1,1 + n4,0,1)*mUV^2*rat( - n2,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,1 + n2,1,1 + n4,1,1)*mUV^4*rat(n2,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,1 + n2,1,n4,0,1)*rat( - n2,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,1 + n2,1,n4,1,1)*mUV^2*rat(n2,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,n2,1,1 + n4,0,1)*rat( - n4,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,n2,1,1 + n4,1,0)*rat(n4,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,n2,1,1 + n4,1,1)*mUV^2*rat(2*ep + n2 + 2*n4 - 1,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,n2,1,1 + n4,2,0)*mUV^2*rat(-2,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,n2,1,2 + n4,0,1)*mUV^2*rat( - n4 - 1,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,n2,1,2 + n4,1,0)*mUV^2*rat(n4 + 1,2*ep + n2 + n4 - 2)
            + uvid(3,1,n1,n2,1,2 + n4,1,1)*mUV^4*rat(3*n4 + 3,2*ep + n2 + n4 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<0}, 1, 0, 1, 1) =
            + uvid(3,1,1 + n1,1 + n2,0,0,1,1)*mUV^2*rat( - n1,2*ep + n1 + n2 - 2)
            + uvid(3,1,1 + n1,1 + n2,1,-1,1,1)*mUV^2*rat( - 2*n1,2*ep + n1 + n2 - 2)
            + uvid(3,1,1 + n1,1 + n2,1,0,1,0)*mUV^2*rat(n1,2*ep + n1 + n2 - 2)
            + uvid(3,1,1 + n1,1 + n2,1,0,1,1)*mUV^4*rat(n1,2*ep + n1 + n2 - 2)
            + uvid(3,1,1 + n1,n2,0,0,1,1)*rat( - n1,2*ep + n1 + n2 - 2)
            + uvid(3,1,1 + n1,n2,1,0,1,0)*rat(n1,2*ep + n1 + n2 - 2)
            + uvid(3,1,1 + n1,n2,1,0,1,1)*mUV^2*rat(n1,2*ep + n1 + n2 - 2)
            + uvid(3,1,n1,1 + n2,0,0,1,1)*rat( - n2,2*ep + n1 + n2 - 2)
            + uvid(3,1,n1,1 + n2,1,0,0,1)*rat(n2,2*ep + n1 + n2 - 2)
            + uvid(3,1,n1,1 + n2,1,0,1,1)*mUV^2*rat(2*ep + n1 + 2*n2 - 1,2*ep + n1 + n2 - 2)
            + uvid(3,1,n1,1 + n2,2,0,0,1)*mUV^2*rat(-2,2*ep + n1 + n2 - 2)
            + uvid(3,1,n1,2 + n2,0,0,1,1)*mUV^2*rat( - n2 - 1,2*ep + n1 + n2 - 2)
            + uvid(3,1,n1,2 + n2,1,0,0,1)*mUV^2*rat(n2 + 1,2*ep + n1 + n2 - 2)
            + uvid(3,1,n1,2 + n2,1,0,1,1)*mUV^4*rat(3*n2 + 3,2*ep + n1 + n2 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, 0, 1, 0, 1, 1)=
            + uvid(3,1,1 + n1,0,0,0,1,1)*rat( - n1,2*ep + n1 - 2)
            + uvid(3,1,1 + n1,0,1,0,1,0)*rat(n1,2*ep + n1 - 2)
            + uvid(3,1,1 + n1,0,1,0,1,1)*mUV^2*rat(2*ep + 2*n1 - 1,2*ep + n1 - 2)
            + uvid(3,1,1 + n1,0,2,0,1,0)*mUV^2*rat(-2,2*ep + n1 - 2)
            + uvid(3,1,2 + n1,0,0,0,1,1)*mUV^2*rat( - n1 - 1,2*ep + n1 - 2)
            + uvid(3,1,2 + n1,0,1,0,1,0)*mUV^2*rat(n1 + 1,2*ep + n1 - 2)
            + uvid(3,1,2 + n1,0,1,0,1,1)*mUV^4*rat(3*n1 + 3,2*ep + n1 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>1}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3,1 + n4,n5, - 1 + n6)*rat(n4,n6 - 1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5, - 1 + n6)*rat(n2,6*n6 - 6)
            + uvid(3,1,1 + n1,1 + n2,n3, - 1 + n4,n5, - 1 + n6)*rat(n2,3*n6 - 3)
            + uvid(3,1,1 + n1,1 + n2,n3,n4, - 1 + n5, - 1 + n6)*rat( - n2,6*n6 - 6)
            + uvid(3,1,1 + n1,1 + n2,n3,n4,n5, - 1 + n6)*mUV^2*rat(n2,6*n6 - 6)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,1 + n5, - 1 + n6)*rat(n5,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*rat( - 2*n4,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 2 + n6)*rat(2*n4,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4,1 + n5, - 2 + n6)*rat( - n5,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - 4*ep - 9*n1 + n2 - 1,6*n6 - 6)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - n1 - 1,2*n6 - 2)
            + uvid(3,1,2 + n1,n2,n3,n4,n5, - 2 + n6)*rat(n1 + 1,2*n6 - 2)
            + uvid(3,1,2 + n1,n2,n3,n4,n5, - 1 + n6)*mUV^2*rat( - 3*n1 - 3,2*n6 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4,n5, - 1 + n6)*rat( - n2,3*n6 - 3)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 1 + n6)*rat( - n4,n6 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<0}, n3?{>0}, n4?{>0}, n5?{>1}, n6?{>0}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3,1 + n4, - 1 + n5,n6)*rat(n4,n5 - 1)
            + uvid(3,1, - 1 + n1,2 + n2,n3,n4, - 1 + n5,n6)*rat( - n2 - 1,3*n5 - 3)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4,n5,n6)*rat(1,3)
            + uvid(3,1,n1,1 + n2,n3,1 + n4, - 2 + n5,n6)*rat(n4,3*n5 - 3)
            + uvid(3,1,n1,1 + n2,n3,1 + n4, - 1 + n5, - 1 + n6)*rat( - n4,3*n5 - 3)
            + uvid(3,1,n1,1 + n2,n3,n4, - 1 + n5,n6)*rat( - 2*ep - 4*n2,3*n5 - 3)
            + uvid(3,1,n1,1 + n2,n3,n4,n5, - 1 + n6)*rat(-1,3)
            + uvid(3,1,n1,2 + n2, - 1 + n3,n4, - 1 + n5,n6)*rat( - n2 - 1,3*n5 - 3)
            + uvid(3,1,n1,2 + n2,n3, - 1 + n4, - 1 + n5,n6)*rat(n2 + 1,3*n5 - 3)
            + uvid(3,1,n1,2 + n2,n3,n4, - 2 + n5,n6)*rat(n2 + 1,3*n5 - 3)
            + uvid(3,1,n1,2 + n2,n3,n4, - 1 + n5,n6)*mUV^2*rat( - 4*n2 - 4,3*n5 - 3)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5,n6)*rat( - n4,n5 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<0}, n3?{>1}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1,1 + n2, - 2 + n3,n4,n5,n6)*rat(n1,2*n3 - 2)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*rat(n1,n3 - 1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - n1,2*n3 - 2)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - n1,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat( - n1,n3 - 1)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5,n6)*rat( - 2*ep - n1 - 3*n2 + 1,2*n3 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4, - 1 + n5,n6)*rat(1,1)
            + uvid(3,1,n1,2 + n2, - 2 + n3,n4,n5,n6)*rat(n2 + 1,2*n3 - 2)
            + uvid(3,1,n1,2 + n2, - 1 + n3,n4, - 1 + n5,n6)*rat( - n2 - 1,2*n3 - 2)
            + uvid(3,1,n1,2 + n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - 3*n2 - 3,2*n3 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{<1}, n3?{>1}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1,1 + n2, - 2 + n3,n4,n5,n6)*rat(n2,2*n3 - 2)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*rat(n2,n3 - 1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*rat( - n2,2*n3 - 2)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - n2,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat( - 2*ep - 3*n1 - n2 + 1,2*n3 - 2)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat(1,1)
            + uvid(3,1,2 + n1,n2, - 2 + n3,n4,n5,n6)*rat(n1 + 1,2*n3 - 2)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - n1 - 1,2*n3 - 2)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - 3*n1 - 3,2*n3 - 2)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5,n6)*rat( - n2,n3 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>1}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - n2,3*n6 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 2 + n6)*mUV^-2*rat( - n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - n1,2*n6 - 2)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(n2,6*n6 - 6)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4,n5, - 1 + n6)*mUV^-2*rat(n2,3*n6 - 3)
            + uvid(3,1,n1,1 + n2,n3,n4, - 1 + n5, - 1 + n6)*mUV^-2*rat( - n2,6*n6 - 6)
            + uvid(3,1,n1,1 + n2,n3,n4,n5, - 1 + n6)*rat(n2,6*n6 - 6)
            + uvid(3,1,n1,n2,n3, - 1 + n4,1 + n5, - 1 + n6)*mUV^-2*rat(n5,3*n6 - 3)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(n4,3*n6 - 3)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 2 + n6)*mUV^-2*rat( - n4,3*n6 - 3)
            + uvid(3,1,n1,n2,n3,n4,1 + n5, - 2 + n6)*mUV^-2*rat( - n5,3*n6 - 3)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - 4*ep - 3*n1 + n2 - 6*n6 + 14,6*n6 - 6)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{>1}, n6?{>0}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat( - n2,3*n5 - 3)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(2*n2,3*n5 - 3)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4, - 1 + n5,n6)*mUV^-2*rat(n2,3*n5 - 3)
            + uvid(3,1,n1,1 + n2,n3,n4, - 2 + n5,n6)*mUV^-2*rat( - 2*n2,3*n5 - 3)
            + uvid(3,1,n1,1 + n2,n3,n4, - 1 + n5,n6)*rat( - n2,3*n5 - 3)
            + uvid(3,1,n1,n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(1,3)
            + uvid(3,1,n1,n2,n3,1 + n4, - 2 + n5,n6)*mUV^-2*rat( - 2*n4,3*n5 - 3)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(2*n4,3*n5 - 3)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat( - 2*ep - n2 - 3*n5 + 7,3*n5 - 3)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat(-1,3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{>1}, n5?{>0}, n6?{>0}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(2*n2,3*n4 - 3)
            + uvid(3,1,n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*mUV^-2*rat( - n2,3*n4 - 3)
            + uvid(3,1,n1,1 + n2,n3, - 2 + n4,n5,n6)*mUV^-2*rat( - 2*n2,3*n4 - 3)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4, - 1 + n5,n6)*mUV^-2*rat(n2,3*n4 - 3)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4,n5,n6)*rat( - n2,3*n4 - 3)
            + uvid(3,1,n1,n2,n3, - 2 + n4,1 + n5,n6)*mUV^-2*rat( - 2*n5,3*n4 - 3)
            + uvid(3,1,n1,n2,n3, - 1 + n4,1 + n5, - 1 + n6)*mUV^-2*rat(2*n5,3*n4 - 3)
            + uvid(3,1,n1,n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat( - 2*ep - n2 - 3*n4 + 7,3*n4 - 3)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat(1,3)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat(-1,3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>1}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1,n2, - 2 + n3,n4,n5,n6)*mUV^-2*rat( - n1,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(n1,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat( - n1,2*n3 - 2)
            + uvid(3,1,n1,1 + n2, - 2 + n3,n4,n5,n6)*mUV^-2*rat( - n2,2*n3 - 2)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(n2,2*n3 - 2)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5,n6)*rat( - n2,2*n3 - 2)
            + uvid(3,1,n1,n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat( - 2*ep - n1 - n2 - 2*n3 + 6,2*n3 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>1}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3,1 + n4,n5, - 1 + n6)*rat( - n2,n6 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,1 + n4,n5, - 1 + n6)*rat(n1,n6 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,1 + n4,n5, - 1 + n6)*rat( - n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 2 + n6)*rat(n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 1 + n6)*mUV^2*rat( - n1,2*n6 - 2)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - n1,n6 - 1)
            + uvid(3,1,n1,1 + n2, - 1 + n3,1 + n4,n5, - 1 + n6)*rat(n2,2*n6 - 2)
            + uvid(3,1,n1,1 + n2,n3,1 + n4, - 1 + n5, - 1 + n6)*rat( - n2,2*n6 - 2)
            + uvid(3,1,n1,1 + n2,n3,1 + n4,n5, - 1 + n6)*mUV^2*rat(n2,2*n6 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4,n5, - 1 + n6)*rat(n2,n6 - 1)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5,n6)*rat(1,1)
            + uvid(3,1,n1,n2,n3,1 + n4,1 + n5, - 2 + n6)*rat( - n5,n6 - 1)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 1 + n6)*rat( - n1 + n2,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,2 + n4, - 1 + n5, - 1 + n6)*rat( - n4 - 1,n6 - 1)
            + uvid(3,1,n1,n2,n3,2 + n4,n5, - 2 + n6)*rat(n4 + 1,n6 - 1)
            + uvid(3,1,n1,n2,n3,n4,1 + n5, - 1 + n6)*rat(n5,n6 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<0}, 1, 1, 1, 1) =
            + uvid(3,1, - 2 + n1,2 + n2,1,1,1,1)*rat(n2 + 1,n1 - 1)
            + uvid(3,1, - 1 + n1,1 + n2,1,0,1,2)*rat(1,n1 - 1)
            + uvid(3,1, - 1 + n1,1 + n2,1,0,2,1)*rat(-1,n1 - 1)
            + uvid(3,1, - 1 + n1,1 + n2,1,1,0,2)*rat(-1,n1 - 1)
            + uvid(3,1, - 1 + n1,1 + n2,1,1,1,1)*rat(n1 - n2 - 2,2*n1 - 2)
            + uvid(3,1, - 1 + n1,1 + n2,1,1,2,0)*rat(1,n1 - 1)
            + uvid(3,1, - 1 + n1,1 + n2,1,2,0,1)*rat(1,n1 - 1)
            + uvid(3,1, - 1 + n1,1 + n2,1,2,1,0)*rat(-1,n1 - 1)
            + uvid(3,1, - 1 + n1,2 + n2,0,1,1,1)*rat( - n2 - 1,2*n1 - 2)
            + uvid(3,1, - 1 + n1,2 + n2,1,0,1,1)*rat( - n2 - 1,n1 - 1)
            + uvid(3,1, - 1 + n1,2 + n2,1,1,0,1)*rat(n2 + 1,2*n1 - 2)
            + uvid(3,1, - 1 + n1,2 + n2,1,1,1,1)*mUV^2*rat( - n2 - 1,2*n1 - 2)
            + uvid(3,1,n1,1 + n2,0,1,1,1)*rat(1,2)
            + uvid(3,1,n1,1 + n2,1,0,1,1)*rat(1,1)
            + uvid(3,1,n1,1 + n2,1,1,1,0)*rat(-1,2)
            + uvid(3,1,n1,1 + n2,1,1,1,1)*mUV^2*rat(1,2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, 0, 1, 1, 1, 1)=
            + uvid(3,1,1 + n1,0,0,1,1,1)*rat( - n1,2*ep + n1 - 2)
            + uvid(3,1,1 + n1,0,1,1,1,0)*rat(n1,2*ep + n1 - 2)
            + uvid(3,1,1 + n1,0,1,1,1,1)*mUV^2*rat(2*ep + 2*n1 - 1,2*ep + n1 - 2)
            + uvid(3,1,1 + n1,0,2,1,1,0)*mUV^2*rat(-2,2*ep + n1 - 2)
            + uvid(3,1,2 + n1,0,0,1,1,1)*mUV^2*rat( - n1 - 1,2*ep + n1 - 2)
            + uvid(3,1,2 + n1,0,1,1,1,0)*mUV^2*rat(n1 + 1,2*ep + n1 - 2)
            + uvid(3,1,2 + n1,0,1,1,1,1)*mUV^4*rat(3*n1 + 3,2*ep + n1 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{<0}, n6?{>1}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,1 + n5, - 1 + n6)*rat( - n1,n6 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,1 + n5, - 1 + n6)*rat(n1,n6 - 1)
            + uvid(3,1,n1,1 + n2,n3,n4,1 + n5, - 1 + n6)*mUV^2*rat( - n2,n6 - 1)
            + uvid(3,1,n1,n2,1 + n3,n4,1 + n5, - 1 + n6)*mUV^2*rat( - n3,n6 - 1)
            + uvid(3,1,n1,n2,n3, - 1 + n4,1 + n5,n6)*rat(1,1)
            + uvid(3,1,n1,n2,n3,1 + n4,1 + n5, - 1 + n6)*mUV^2*rat(n4,n6 - 1)
            + uvid(3,1,n1,n2,n3,n4,1 + n5, - 1 + n6)*rat( - ep - n2 - n3 + n4 - n5 + 1,n6 - 1)
            + uvid(3,1,n1,n2,n3,n4,2 + n5, - 1 + n6)*mUV^2*rat( - n5 - 1,n6 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{<1}, n6?{>1}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,1 + n5, - 1 + n6)*rat( - n5,n6 - 1)
            + uvid(3,1,1 + n1,1 + n2,n3,n4,n5, - 1 + n6)*mUV^2*rat( - n2,n6 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,1 + n5, - 1 + n6)*rat(n5,n6 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,n2,1 + n3,n4,n5, - 1 + n6)*mUV^2*rat(n3,n6 - 1)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 1 + n6)*mUV^2*rat( - n4,n6 - 1)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - ep - n1 - n2 + n3 - n4 + 1,n6 - 1)
            + uvid(3,1,2 + n1,n2,n3,n4,n5, - 1 + n6)*mUV^2*rat( - n1 - 1,n6 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{>1}, n5?{<0}, n6?{>0}) =
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,1 + n5,n6)*rat( - n1,n4 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,1 + n5, - 1 + n6)*rat(n1,n4 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,1 + n5,n6)*mUV^2*rat( - n1,n4 - 1)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4,1 + n5,n6)*mUV^2*rat( - 2*n2,n4 - 1)
            + uvid(3,1,n1,n2,1 + n3, - 1 + n4,1 + n5,n6)*mUV^2*rat( - 2*n3,n4 - 1)
            + uvid(3,1,n1,n2,n3, - 1 + n4,1 + n5,n6)*rat( - 4*ep - n1 - 2*n2 - 2*n3 - n4 - 2*n5 + 7,n4 - 1)
            + uvid(3,1,n1,n2,n3, - 1 + n4,2 + n5,n6)*mUV^2*rat( - 2*n5 - 2,n4 - 1)
            + uvid(3,1,n1,n2,n3,n4,1 + n5, - 1 + n6)*rat(1,1)
            + uvid(3,1,n1,n2,n3,n4,1 + n5,n6)*mUV^2*rat(-1,1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{>0}, n3?{>0}, n4?{>1}, n5?{<1}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,1 + n5,n6)*rat(n5,n4 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,1 + n2,n3, - 1 + n4,n5,n6)*mUV^2*rat(2*n2,n4 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,1 + n5,n6)*rat( - n5,n4 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,1 + n5,n6)*mUV^2*rat(n5,n4 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5,n6)*rat(2*ep + 2*n2 + n4 + n5 - 5,n4 - 1)
            + uvid(3,1,1 + n1,n2,n3,n4,n5,n6)*mUV^2*rat(1,1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>1}, n4?{>0}, n5?{<0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2, - 1 + n3,n4,1 + n5,n6)*rat(n1,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,1 + n5,n6)*rat( - n1,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,1 + n5,n6)*mUV^2*rat(n1,n3 - 1)
            + uvid(3,1,n1, - 1 + n2,n3,n4,1 + n5,n6)*rat(1,1)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,1 + n5,n6)*mUV^2*rat(2*n2,n3 - 1)
            + uvid(3,1,n1,n2, - 1 + n3,n4,1 + n5,n6)*rat(2*ep + n1 + 2*n2 + n3 - 5,n3 - 1)
            + uvid(3,1,n1,n2,n3,n4,1 + n5,n6)*mUV^2*rat(1,1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{>0}, n3?{>1}, n4?{>0}, n5?{<1}, n6?{>0}) =
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - 2*n2,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,1 + n5,n6)*rat( - n5,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,1 + n4,n5,n6)*mUV^2*rat( - 2*n4,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,1 + n5, - 1 + n6)*rat(n5,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,1 + n5,n6)*mUV^2*rat( - n5,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat( - 4*ep - 2*n1 - 2*n2 - n3 - 2*n4 - n5 + 7,n3 - 1)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat(1,1)
            + uvid(3,1,1 + n1,n2,n3,n4,n5,n6)*mUV^2*rat(-1,1)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - 2*n1 - 2,n3 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>1}, n3?{>0}, n4?{>0}, n5?{<0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2, - 1 + n3,n4,1 + n5,n6)*rat(n1,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,1 + n5, - 1 + n6)*rat( - n1,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,1 + n5,n6)*mUV^2*rat(n1,n2 - 1)
            + uvid(3,1,n1, - 1 + n2,1 + n3,n4,1 + n5,n6)*mUV^2*rat(2*n3,n2 - 1)
            + uvid(3,1,n1, - 1 + n2,n3,n4,1 + n5,n6)*rat(2*ep + n1 + n2 + 2*n3 - 5,n2 - 1)
            + uvid(3,1,n1,n2, - 1 + n3,n4,1 + n5,n6)*rat(1,1)
            + uvid(3,1,n1,n2,n3,n4,1 + n5,n6)*mUV^2*rat(1,1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{>1}, n3?{>0}, n4?{>0}, n5?{<1}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,1 + n5,n6)*rat(n5,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,1 + n4,n5,n6)*mUV^2*rat(2*n4,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,1 + n5, - 1 + n6)*rat( - n5,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,1 + n5,n6)*mUV^2*rat(n5,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5,n6)*rat(2*ep + n2 + 2*n4 + n5 - 5,n2 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,n2,n3,n4,n5,n6)*mUV^2*rat(1,1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{<1}, n6?{>1}) =
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - n1,n6 - 1)
            + uvid(3,1,n1,1 + n2,n3,n4,n5, - 1 + n6)*rat( - n2,n6 - 1)
            + uvid(3,1,n1,n2,1 + n3,n4,n5, - 1 + n6)*rat( - n3,n6 - 1)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 1 + n6)*rat( - n4,n6 - 1)
            + uvid(3,1,n1,n2,n3,n4,1 + n5, - 1 + n6)*rat( - n5,n6 - 1)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - 3*ep - n1 - n2 - n3 - n4 - n5 - n6 + 7,n6 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,0, n2?{>0}, n3?{>0}, n4?{>2}, 0, 1)=
            + uvid(3,1,0, - 1 + n2,1 + n3, - 1 + n4,0,1)*mUV^-2*rat( - n3,2*n4 - 2)
            + uvid(3,1,0,1 + n2,1 + n3, - 2 + n4,0,1)*n3*rat( - n2,n4^2 - 3*n4 + 2)
            + uvid(3,1,0,1 + n2,n3, - 1 + n4,0,1)*rat( - n2,n4 - 1)
            + uvid(3,1,0,n2,1 + n3, - 2 + n4,0,1)*mUV^-2*rat( - ep*n3 + 3*n3,n4^2 - 3*n4 + 2)
            + uvid(3,1,0,n2,1 + n3, - 2 + n4,0,1)*mUV^-2*n4*rat( - n3,2*n4^2 - 6*n4 + 4)
            + uvid(3,1,0,n2,1 + n3, - 2 + n4,0,1)*mUV^-2*n3*rat( - n2,n4^2 - 3*n4 + 2)
            + uvid(3,1,0,n2,1 + n3, - 1 + n4,0,0)*mUV^-2*rat(n3,2*n4 - 2)
            + uvid(3,1,0,n2,1 + n3, - 1 + n4,0,1)*rat( - n3,n4 - 1)
            + uvid(3,1,0,n2,n3, - 1 + n4,0,1)*mUV^-2*rat( - 4*ep - 2*n2 - n3 - 2*n4 + 10,2*n4 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,0, n2?{>0}, n3?{>1}, 2, 0, 1) =
            + uvid(3,1,0, - 1 + n2,n3,2,0,1)*mUV^-2*rat(-1,2)
            + uvid(3,1,0,1 + n2, - 1 + n3,2,0,1)*rat( - n2,n3 - 1)
            + uvid(3,1,0,1 + n2,n3,1,0,1)*rat( - n2,1)
            + uvid(3,1,0,n2, - 1 + n3,2,0,1)*mUV^-2*rat( - 2*ep - 2*n2 - n3 + 5,2*n3 - 2)
            + uvid(3,1,0,n2,1 + n3,1,0,1)*rat( - n3,1)
            + uvid(3,1,0,n2,n3,1,0,1)*mUV^-2*rat( - 4*ep - 2*n2 - 2*n3 + 7,2)
            + uvid(3,1,0,n2,n3,2,0,0)*mUV^-2*rat(1,2)
            ;
        id ifmatch->end`i' uvid(3,1,0, n2?{>0}, n3?{>2}, 1, 0, 1) =
            + uvid(3,1,0,1 + n2, - 2 + n3,1,0,1)*mUV^-2*rat(ep*n2 + n2^2 - n2,n3^2 - 3*n3 + 2)
            + uvid(3,1,0,2 + n2, - 2 + n3,1,0,1)*rat(n2^2 + n2,n3^2 - 3*n3 + 2)
            + uvid(3,1,0,n2, - 1 + n3,1,0,1)*mUV^-2*rat( - ep - n3 + 3,n3 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, 1, 1, 1, n5?{<0}, 1) =
            + uvid(3,1,1 + n1,1,1,1,1 + n5,1)*mUV^4*rat( - n1,3*ep + n1 + n5 - 2)
            + uvid(3,1,1 + n1,1,1,1,n5,1)*mUV^2*rat( - n1,3*ep + n1 + n5 - 2)
            + uvid(3,1,n1,0,2,1,1 + n5,1)*mUV^2*rat(-1,3*ep + n1 + n5 - 2)
            + uvid(3,1,n1,1,1,0,1 + n5,2)*mUV^2*rat(-1,3*ep + n1 + n5 - 2)
            + uvid(3,1,n1,1,1,1,1 + n5,1)*mUV^2*rat(ep - n1 + 2*n5 + 1,3*ep + n1 + n5 - 2)
            + uvid(3,1,n1,1,1,1,2 + n5,1)*mUV^4*rat(3*n5 + 3,3*ep + n1 + n5 - 2)
            + uvid(3,1,n1,1,1,2,1 + n5,0)*mUV^2*rat(-1,3*ep + n1 + n5 - 2)
            + uvid(3,1,n1,2,0,1,1 + n5,1)*mUV^2*rat(-1,3*ep + n1 + n5 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,0, n2?{>0}, 1, 2, 0, 1)=
            + uvid(3,1,0, - 1 + n2,1,2,0,1)*mUV^-2*rat( - 2*ep - n2 + 3,2*ep + n2 - 2)
            + uvid(3,1,0,1 + n2,1,1,0,1)*rat( - 3*n2,1)
            + uvid(3,1,0,n2,0,2,0,2)*rat(-1,2*ep + n2 - 2)
            + uvid(3,1,0,n2,1,1,0,1)*mUV^-2*rat( - 9*ep*n2 - 8*ep - 3*n2^2 + 12*n2 - 11,2*ep + n2 - 2)
            + uvid(3,1,0,n2,1,1,0,1)*ep*mUV^-2*rat(24,2*ep + n2 - 2)
            + uvid(3,1,0,n2,1,1,0,1)*ep^2*mUV^-2*rat(-6,2*ep + n2 - 2)
            + uvid(3,1,0,n2,2,2,0,0)*rat(-1,2*ep + n2 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,0, n2?{>0}, 2, 1, 0, 1)=
            + uvid(3,1,0, - 1 + n2,2,1,0,1)*mUV^-2*rat( - 2*ep - n2 + 3,2*ep + n2 - 2)
            + uvid(3,1,0,1 + n2,1,1,0,1)*rat( - 3*n2,1)
            + uvid(3,1,0,n2,1,1,0,1)*mUV^-2*rat( - 9*ep*n2 - 8*ep - 3*n2^2 + 12*n2 - 11,2*ep + n2 - 2)
            + uvid(3,1,0,n2,1,1,0,1)*ep*mUV^-2*rat(24,2*ep + n2 - 2)
            + uvid(3,1,0,n2,1,1,0,1)*ep^2*mUV^-2*rat(-6,2*ep + n2 - 2)
            + uvid(3,1,0,n2,2,0,0,2)*rat(-1,2*ep + n2 - 2)
            + uvid(3,1,0,n2,2,2,0,0)*rat(-1,2*ep + n2 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,0, n2?{>1}, 1, 1, 0, 1)=
            + uvid(3,1,0, - 2 + n2,2,1,0,1)*mUV^-2*rat( - 6*ep - 3*n2 + 12,16*ep*n2 - 16*ep + 8*n2^2 - 32*n2 + 24)
            + uvid(3,1,0, - 1 + n2,1,1,0,1)*mUV^-2*rat( - 22*ep*n2 + 10*ep - 8*n2^2 + 47*n2 - 66,16*ep*n2 - 16*ep + 8*n2^2 - 32*n2 + 24)
            + uvid(3,1,0, - 1 + n2,1,1,0,1)*ep*mUV^-2*rat(6,2*ep*n2 - 2*ep + n2^2 - 4*n2 + 3)
            + uvid(3,1,0, - 1 + n2,1,1,0,1)*ep^2*mUV^-2*rat(-3,4*ep*n2 - 4*ep + 2*n2^2 - 8*n2 + 6)
            + uvid(3,1,0, - 1 + n2,2,0,0,1)*mUV^-2*rat( - 3*n2 + 6,16*ep*n2 - 16*ep + 8*n2^2 - 32*n2 + 24)
            + uvid(3,1,0, - 1 + n2,2,1,0,0)*mUV^-2*rat(3*n2 - 6,16*ep*n2 - 16*ep + 8*n2^2 - 32*n2 + 24)
            + uvid(3,1,0, - 1 + n2,2,2,0,0)*rat(-3,8*ep*n2 - 8*ep + 4*n2^2 - 16*n2 + 12)
            + uvid(3,1,0,n2,0,1,0,2)*rat(-1,8*ep + 4*n2 - 12)
            + uvid(3,1,0,n2,1,0,0,2)*rat(1,8*ep + 4*n2 - 12)
            + uvid(3,1,0,n2,2,0,0,1)*rat(-1,4*ep + 2*n2 - 6)
            + uvid(3,1,0,n2,2,1,0,0)*rat(1,4*ep + 2*n2 - 6)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, 1, 1, 1, 0, 1)=
            + uvid(3,1,1 + n1,0,1,2,0,1)*mUV^2*rat(-1,3*ep + n1 - 2)
            + uvid(3,1,1 + n1,1,0,1,0,2)*mUV^2*rat(-1,3*ep + n1 - 2)
            + uvid(3,1,1 + n1,1,1,1,0,1)*mUV^2*rat(ep + 2*n1 + 1,3*ep + n1 - 2)
            + uvid(3,1,1 + n1,1,2,1,0,0)*mUV^2*rat(-1,3*ep + n1 - 2)
            + uvid(3,1,1 + n1,2,1,0,0,1)*mUV^2*rat(-1,3*ep + n1 - 2)
            + uvid(3,1,2 + n1,1,1,1,0,1)*mUV^4*rat(3*n1 + 3,3*ep + n1 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>1}) =
            + uvid(3,1,1 + n1, - 1 + n2,1 + n3,n4,n5, - 1 + n6)*rat(n3,n6 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,1 + n5, - 1 + n6)*rat( - n5,n6 - 1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - n2,n6 - 1)
            + uvid(3,1,1 + n1,1 + n2,n3,n4, - 1 + n5, - 1 + n6)*rat(n2,n6 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,1 + n5, - 1 + n6)*rat(n5,n6 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,n2,1 + n3,n4, - 1 + n5, - 1 + n6)*rat( - n3,n6 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5,n6)*rat(2,3)
            + uvid(3,1,1 + n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*rat(n4,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 2 + n6)*rat( - n4,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4, - 1 + n5,n6)*rat(-2,3)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - ep - 2*n1,3*n6 - 3)
            + uvid(3,1,2 + n1, - 1 + n2,n3,n4,n5, - 1 + n6)*rat(n1 + 1,3*n6 - 3)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - 2*n1 - 2,3*n6 - 3)
            + uvid(3,1,2 + n1,n2,n3, - 1 + n4,n5, - 1 + n6)*rat( - n1 - 1,3*n6 - 3)
            + uvid(3,1,2 + n1,n2,n3,n4,n5, - 2 + n6)*rat(2*n1 + 2,3*n6 - 3)
            + uvid(3,1,2 + n1,n2,n3,n4,n5, - 1 + n6)*mUV^2*rat( - 2*n1 - 2,3*n6 - 3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{>0}, n3?{>0}, n4?{>1}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2,1 + n3, - 1 + n4,n5,n6)*rat( - n3,n4 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,1 + n5,n6)*rat(n5,n4 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*rat(n2,n4 - 1)
            + uvid(3,1,1 + n1,1 + n2,n3, - 1 + n4, - 1 + n5,n6)*rat( - n2,n4 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,1 + n5,n6)*rat( - n5,n4 - 1)
            + uvid(3,1,1 + n1,n2,1 + n3, - 1 + n4, - 1 + n5,n6)*rat(n3,n4 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 2 + n4,n5,1 + n6)*rat( - n6,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4, - 1 + n5,1 + n6)*rat(n6,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5,n6)*rat( - ep - 2*n1,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4, - 1 + n5,n6)*rat(-2,3)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat(2,3)
            + uvid(3,1,2 + n1, - 1 + n2,n3, - 1 + n4,n5,n6)*rat( - 2*n1 - 2,3*n4 - 3)
            + uvid(3,1,2 + n1,n2, - 1 + n3, - 1 + n4,n5,n6)*rat(n1 + 1,3*n4 - 3)
            + uvid(3,1,2 + n1,n2,n3, - 2 + n4,n5,n6)*rat(2*n1 + 2,3*n4 - 3)
            + uvid(3,1,2 + n1,n2,n3, - 1 + n4,n5, - 1 + n6)*rat( - n1 - 1,3*n4 - 3)
            + uvid(3,1,2 + n1,n2,n3, - 1 + n4,n5,n6)*mUV^2*rat( - 2*n1 - 2,3*n4 - 3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{>0}, n3?{>1}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5,n6)*rat(2,3)
            + uvid(3,1,1 + n1,1 + n2, - 2 + n3,n4,n5,n6)*rat( - n2,3*n3 - 3)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*rat(n2,3*n3 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,1 + n5,n6)*rat( - n5,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,n5,1 + n6)*rat(n6,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,1 + n4, - 1 + n5,n6)*rat(n4,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,1 + n4,n5, - 1 + n6)*rat( - n4,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4, - 1 + n5,1 + n6)*rat( - n6,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,1 + n5, - 1 + n6)*rat(n5,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat( - ep - 2*n1,3*n3 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4, - 1 + n5,n6)*rat(-2,3)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat(1,1)
            + uvid(3,1,2 + n1, - 1 + n2, - 1 + n3,n4,n5,n6)*rat( - n1 - 1,3*n3 - 3)
            + uvid(3,1,2 + n1,n2, - 2 + n3,n4,n5,n6)*rat(2*n1 + 2,3*n3 - 3)
            + uvid(3,1,2 + n1,n2, - 1 + n3, - 1 + n4,n5,n6)*rat(n1 + 1,3*n3 - 3)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - 2*n1 - 2,3*n3 - 3)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5,n6)*mUV^2*rat( - 2*n1 - 2,3*n3 - 3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, n2?{>1}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 2 + n2,1 + n3,n4,n5,n6)*rat( - n3,3*n2 - 3)
            + uvid(3,1,1 + n1, - 1 + n2,1 + n3,n4, - 1 + n5,n6)*rat(n3,3*n2 - 3)
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,1 + n5,n6)*rat(n5,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,n5,1 + n6)*rat( - n6,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,1 + n4, - 1 + n5,n6)*rat( - n4,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,1 + n4,n5, - 1 + n6)*rat(n4,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4, - 1 + n5,1 + n6)*rat(n6,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,1 + n5, - 1 + n6)*rat( - n5,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5,n6)*rat( - ep - 2*n1,3*n2 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat(2,3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,n2,n3,n4, - 1 + n5,n6)*rat(-2,3)
            + uvid(3,1,2 + n1, - 2 + n2,n3,n4,n5,n6)*rat(2*n1 + 2,3*n2 - 3)
            + uvid(3,1,2 + n1, - 1 + n2, - 1 + n3,n4,n5,n6)*rat( - n1 - 1,3*n2 - 3)
            + uvid(3,1,2 + n1, - 1 + n2,n3, - 1 + n4,n5,n6)*rat( - 2*n1 - 2,3*n2 - 3)
            + uvid(3,1,2 + n1, - 1 + n2,n3,n4,n5, - 1 + n6)*rat(n1 + 1,3*n2 - 3)
            + uvid(3,1,2 + n1, - 1 + n2,n3,n4,n5,n6)*mUV^2*rat( - 2*n1 - 2,3*n2 - 3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>1}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - n1,3*n6 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(2*n1,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5, - 1 + n6)*mUV^-2*rat(n1,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 2 + n6)*mUV^-2*rat( - 2*n1,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat( - n1,3*n6 - 3)
            + uvid(3,1,n1,n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(1,3)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(2*n4,3*n6 - 3)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 2 + n6)*mUV^-2*rat( - 2*n4,3*n6 - 3)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat(-1,3)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - 2*ep - n1 - 3*n6 + 7,3*n6 - 3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{>1}, n6?{>0}) =
            + uvid(3,1,1 + n1,n2,n3,n4, - 1 + n5,n6)*rat(n1,3*n5 - 3)
            + uvid(3,1,n1, - 1 + n2,1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(n3,3*n5 - 3)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(n2,3*n5 - 3)
            + uvid(3,1,n1,1 + n2,n3,n4, - 2 + n5,n6)*mUV^-2*rat( - n2,3*n5 - 3)
            + uvid(3,1,n1,n2,1 + n3,n4, - 2 + n5,n6)*mUV^-2*rat( - n3,3*n5 - 3)
            + uvid(3,1,n1,n2,n3, - 1 + n4, - 1 + n5,1 + n6)*mUV^-2*rat(n6,3*n5 - 3)
            + uvid(3,1,n1,n2,n3,1 + n4, - 2 + n5,n6)*mUV^-2*rat( - n4,3*n5 - 3)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(n4,3*n5 - 3)
            + uvid(3,1,n1,n2,n3,n4, - 2 + n5,1 + n6)*mUV^-2*rat( - n6,3*n5 - 3)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat( - ep + n1 - 3*n5 + 5,3*n5 - 3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{>1}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(2*n1,3*n4 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,n5,n6)*mUV^-2*rat( - n1,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 2 + n4,n5,n6)*mUV^-2*rat( - 2*n1,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5, - 1 + n6)*mUV^-2*rat(n1,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5,n6)*rat( - n1,3*n4 - 3)
            + uvid(3,1,n1,n2,n3, - 2 + n4,n5,1 + n6)*mUV^-2*rat( - 2*n6,3*n4 - 3)
            + uvid(3,1,n1,n2,n3, - 1 + n4, - 1 + n5,1 + n6)*mUV^-2*rat(2*n6,3*n4 - 3)
            + uvid(3,1,n1,n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat( - 2*ep - n1 - 3*n4 + 7,3*n4 - 3)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat(-1,3)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat(1,3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>1}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat(n1,3*n3 - 3)
            + uvid(3,1,1 + n1,n2, - 2 + n3,n4,n5,n6)*mUV^-2*rat( - 2*n1,3*n3 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,n5,n6)*mUV^-2*rat( - n1,3*n3 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(2*n1,3*n3 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat( - n1,3*n3 - 3)
            + uvid(3,1,n1, - 1 + n2,n3,n4,n5,n6)*mUV^-2*rat(1,3)
            + uvid(3,1,n1,1 + n2, - 2 + n3,n4,n5,n6)*mUV^-2*rat( - 2*n2,3*n3 - 3)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(2*n2,3*n3 - 3)
            + uvid(3,1,n1,n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat( - 2*ep - n1 - 3*n3 + 7,3*n3 - 3)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat(-1,3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>1}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 2 + n2,n3,n4,n5,n6)*mUV^-2*rat( - 2*n1,3*n2 - 3)
            + uvid(3,1,1 + n1, - 1 + n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat(n1,3*n2 - 3)
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(2*n1,3*n2 - 3)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - n1,3*n2 - 3)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5,n6)*rat( - n1,3*n2 - 3)
            + uvid(3,1,n1, - 2 + n2,1 + n3,n4,n5,n6)*mUV^-2*rat( - 2*n3,3*n2 - 3)
            + uvid(3,1,n1, - 1 + n2,1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(2*n3,3*n2 - 3)
            + uvid(3,1,n1, - 1 + n2,n3,n4,n5,n6)*mUV^-2*rat( - 2*ep - n1 - 3*n2 + 7,3*n2 - 3)
            + uvid(3,1,n1,n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat(1,3)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat(-1,3)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{<0}, 1, 1, 1, 1, 1) =
            + uvid(3,1,1 + n1,0,1,1,1,1)*rat( - 2*n1,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,0,2,1,1,1)*mUV^2*rat(1,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,0,1,1,1)*rat(n1,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,1,0,1,1)*rat(2*n1,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,1,0,1,2)*mUV^2*rat(3,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,1,0,2,1)*mUV^2*rat(-3,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,1,1,0,2)*mUV^2*rat(-3,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,1,1,1,0)*rat( - n1,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,1,1,1,1)*mUV^2*rat(ep + n1,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,1,1,2,0)*mUV^2*rat(3,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,1,2,0,1)*mUV^2*rat(3,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,1,2,1,0)*mUV^2*rat(-3,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,1,2,1,0,1)*mUV^2*rat(-1,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,2,0,1,1,1)*mUV^2*rat(-2,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,2,1,0,1,1)*mUV^2*rat(-3,2*ep + n1 - 1)
            + uvid(3,1,1 + n1,2,1,1,0,1)*mUV^2*rat(2,2*ep + n1 - 1)
            + uvid(3,1,2 + n1,0,1,1,1,1)*mUV^2*rat( - 2*n1 - 2,2*ep + n1 - 1)
            + uvid(3,1,2 + n1,1,0,1,1,1)*mUV^2*rat(n1 + 1,2*ep + n1 - 1)
            + uvid(3,1,2 + n1,1,1,0,1,1)*mUV^2*rat(2*n1 + 2,2*ep + n1 - 1)
            + uvid(3,1,2 + n1,1,1,1,1,0)*mUV^2*rat( - n1 - 1,2*ep + n1 - 1)
            + uvid(3,1,2 + n1,1,1,1,1,1)*mUV^4*rat(2*n1 + 2,2*ep + n1 - 1)
            + uvid(3,1,n1,0,2,1,1,1)*rat(-2,2*ep + n1 - 1)
            + uvid(3,1,n1,1,2,1,0,1)*rat(2,2*ep + n1 - 1)
            + uvid(3,1,n1,2,0,1,1,1)*rat(1,2*ep + n1 - 1)
            + uvid(3,1,n1,2,1,1,0,1)*rat(-1,2*ep + n1 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>1}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat(n2,2*n6 - 2)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - 2*n1,3*n6 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(5*n1,6*n6 - 6)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5, - 1 + n6)*mUV^-2*rat(2*n1,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 2 + n6)*mUV^-2*rat( - 5*n1,6*n6 - 6)
            + uvid(3,1,n1, - 1 + n2,1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(n3,6*n6 - 6)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat( - n2,3*n6 - 3)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4,n5, - 1 + n6)*mUV^-2*rat( - n2,2*n6 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(n2,3*n6 - 3)
            + uvid(3,1,n1,n2,1 + n3,n4, - 1 + n5, - 1 + n6)*mUV^-2*rat( - n3,6*n6 - 6)
            + uvid(3,1,n1,n2,n3, - 1 + n4,1 + n5, - 1 + n6)*mUV^-2*rat( - n5,2*n6 - 2)
            + uvid(3,1,n1,n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(5,6)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(7*n4,6*n6 - 6)
            + uvid(3,1,n1,n2,n3,1 + n4,n5, - 2 + n6)*mUV^-2*rat( - 7*n4,6*n6 - 6)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat(-5,6)
            + uvid(3,1,n1,n2,n3,n4,1 + n5, - 2 + n6)*mUV^-2*rat(n5,2*n6 - 2)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - ep - 2*n6 + 4,2*n6 - 2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{>1}, n6?{>0}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat( - n2,2*n5 - 2)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat(n1,3*n5 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat( - n1,6*n5 - 6)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4, - 1 + n5,n6)*mUV^-2*rat( - n1,3*n5 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(n1,6*n5 - 6)
            + uvid(3,1,n1, - 1 + n2,1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(n3,6*n5 - 6)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(2*n2,3*n5 - 3)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4, - 1 + n5,n6)*mUV^-2*rat(n2,2*n5 - 2)
            + uvid(3,1,n1,1 + n2,n3,n4, - 2 + n5,n6)*mUV^-2*rat( - 2*n2,3*n5 - 3)
            + uvid(3,1,n1,n2,1 + n3,n4, - 2 + n5,n6)*mUV^-2*rat( - n3,6*n5 - 6)
            + uvid(3,1,n1,n2,n3, - 1 + n4, - 1 + n5,1 + n6)*mUV^-2*rat( - n6,6*n5 - 6)
            + uvid(3,1,n1,n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(1,2)
            + uvid(3,1,n1,n2,n3,1 + n4, - 2 + n5,n6)*mUV^-2*rat( - 5*n4,6*n5 - 6)
            + uvid(3,1,n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*mUV^-2*rat(5*n4,6*n5 - 6)
            + uvid(3,1,n1,n2,n3,n4, - 2 + n5,1 + n6)*mUV^-2*rat(n6,6*n5 - 6)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat( - ep - 2*n5 + 4,2*n5 - 2)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat(-1,2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{>1}, n5?{>0}, n6?{>0}) =
            + uvid(3,1, - 1 + n1,1 + n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(n2,2*n4 - 2)
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(n1,3*n4 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,n5,n6)*mUV^-2*rat( - n1,6*n4 - 6)
            + uvid(3,1,1 + n1,n2,n3, - 2 + n4,n5,n6)*mUV^-2*rat( - n1,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5, - 1 + n6)*mUV^-2*rat(n1,6*n4 - 6)
            + uvid(3,1,n1, - 1 + n2,1 + n3, - 1 + n4,n5,n6)*mUV^-2*rat(n3,6*n4 - 6)
            + uvid(3,1,n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*mUV^-2*rat( - n2,3*n4 - 3)
            + uvid(3,1,n1,1 + n2,n3, - 2 + n4,n5,n6)*mUV^-2*rat( - n2,2*n4 - 2)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4, - 1 + n5,n6)*mUV^-2*rat(n2,3*n4 - 3)
            + uvid(3,1,n1,n2,1 + n3, - 1 + n4, - 1 + n5,n6)*mUV^-2*rat( - n3,6*n4 - 6)
            + uvid(3,1,n1,n2,n3, - 2 + n4,1 + n5,n6)*mUV^-2*rat( - n5,2*n4 - 2)
            + uvid(3,1,n1,n2,n3, - 2 + n4,n5,1 + n6)*mUV^-2*rat( - n6,6*n4 - 6)
            + uvid(3,1,n1,n2,n3, - 1 + n4, - 1 + n5,1 + n6)*mUV^-2*rat(n6,6*n4 - 6)
            + uvid(3,1,n1,n2,n3, - 1 + n4,1 + n5, - 1 + n6)*mUV^-2*rat(n5,2*n4 - 2)
            + uvid(3,1,n1,n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat( - ep - 2*n4 + 4,2*n4 - 2)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat(1,6)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat(-1,6)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>0}, n3?{>1}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1, - 1 + n1,1 + n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat(n2,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 2 + n3,n4,n5,n6)*mUV^-2*rat( - n1,2*n3 - 2)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*mUV^-2*rat(n1,2*n3 - 2)
            + uvid(3,1,n1, - 1 + n2,n3,n4,n5,n6)*mUV^-2*rat(1,2)
            + uvid(3,1,n1,1 + n2, - 2 + n3,n4,n5,n6)*mUV^-2*rat( - n2,n3 - 1)
            + uvid(3,1,n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*mUV^-2*rat( - n2,2*n3 - 2)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(n2,n3 - 1)
            + uvid(3,1,n1,n2, - 1 + n3, - 1 + n4,1 + n5,n6)*mUV^-2*rat( - n5,2*n3 - 2)
            + uvid(3,1,n1,n2, - 1 + n3, - 1 + n4,n5,1 + n6)*mUV^-2*rat(n6,2*n3 - 2)
            + uvid(3,1,n1,n2, - 1 + n3,1 + n4, - 1 + n5,n6)*mUV^-2*rat(n4,2*n3 - 2)
            + uvid(3,1,n1,n2, - 1 + n3,1 + n4,n5, - 1 + n6)*mUV^-2*rat( - n4,2*n3 - 2)
            + uvid(3,1,n1,n2, - 1 + n3,n4, - 1 + n5,1 + n6)*mUV^-2*rat( - n6,2*n3 - 2)
            + uvid(3,1,n1,n2, - 1 + n3,n4,1 + n5, - 1 + n6)*mUV^-2*rat(n5,2*n3 - 2)
            + uvid(3,1,n1,n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat( - ep - 2*n3 + 4,2*n3 - 2)
            + uvid(3,1,n1,n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat(-1,2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>1}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1, - 1 + n1,n2,n3,n4,n5,n6)*mUV^-2*rat(1,2)
            + uvid(3,1,1 + n1, - 2 + n2,n3,n4,n5,n6)*mUV^-2*rat( - n1,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat(n1,2*n2 - 2)
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(n1,n2 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat( - n1,2*n2 - 2)
            + uvid(3,1,n1, - 2 + n2,1 + n3,n4,n5,n6)*mUV^-2*rat( - n3,2*n2 - 2)
            + uvid(3,1,n1, - 1 + n2,1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(n3,2*n2 - 2)
            + uvid(3,1,n1, - 1 + n2,n3, - 1 + n4,1 + n5,n6)*mUV^-2*rat( - n5,2*n2 - 2)
            + uvid(3,1,n1, - 1 + n2,n3, - 1 + n4,n5,1 + n6)*mUV^-2*rat(n6,2*n2 - 2)
            + uvid(3,1,n1, - 1 + n2,n3,1 + n4, - 1 + n5,n6)*mUV^-2*rat(n4,2*n2 - 2)
            + uvid(3,1,n1, - 1 + n2,n3,1 + n4,n5, - 1 + n6)*mUV^-2*rat( - n4,2*n2 - 2)
            + uvid(3,1,n1, - 1 + n2,n3,n4, - 1 + n5,1 + n6)*mUV^-2*rat( - n6,2*n2 - 2)
            + uvid(3,1,n1, - 1 + n2,n3,n4,1 + n5, - 1 + n6)*mUV^-2*rat(n5,2*n2 - 2)
            + uvid(3,1,n1, - 1 + n2,n3,n4,n5,n6)*mUV^-2*rat( - ep - 2*n2 + 4,2*n2 - 2)
            + uvid(3,1,n1,n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(-1,2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{>1}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1, - 2 + n1,1 + n2,n3,n4,n5,n6)*mUV^-2*rat( - 3*n2,2*n1 - 2)
            + uvid(3,1, - 1 + n1, - 1 + n2,1 + n3,n4,n5,n6)*mUV^-2*rat( - n3,2*n1 - 2)
            + uvid(3,1, - 1 + n1,1 + n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat(n2,n1 - 1)
            + uvid(3,1, - 1 + n1,1 + n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(3*n2,2*n1 - 2)
            + uvid(3,1, - 1 + n1,1 + n2,n3,n4, - 1 + n5,n6)*mUV^-2*rat( - n2,n1 - 1)
            + uvid(3,1, - 1 + n1,n2,1 + n3,n4, - 1 + n5,n6)*mUV^-2*rat(n3,2*n1 - 2)
            + uvid(3,1, - 1 + n1,n2,n3, - 1 + n4,1 + n5,n6)*mUV^-2*rat(3*n5,2*n1 - 2)
            + uvid(3,1, - 1 + n1,n2,n3, - 1 + n4,n5,1 + n6)*mUV^-2*rat( - 3*n6,2*n1 - 2)
            + uvid(3,1, - 1 + n1,n2,n3,1 + n4, - 1 + n5,n6)*mUV^-2*rat( - 3*n4,2*n1 - 2)
            + uvid(3,1, - 1 + n1,n2,n3,1 + n4,n5, - 1 + n6)*mUV^-2*rat(3*n4,2*n1 - 2)
            + uvid(3,1, - 1 + n1,n2,n3,n4, - 1 + n5,1 + n6)*mUV^-2*rat(3*n6,2*n1 - 2)
            + uvid(3,1, - 1 + n1,n2,n3,n4,1 + n5, - 1 + n6)*mUV^-2*rat( - 3*n5,2*n1 - 2)
            + uvid(3,1, - 1 + n1,n2,n3,n4,n5,n6)*mUV^-2*rat( - ep - 2*n1 + 4,2*n1 - 2)
            + uvid(3,1,n1, - 1 + n2,n3,n4,n5,n6)*mUV^-2*rat(1,1)
            + uvid(3,1,n1,n2, - 1 + n3,n4,n5,n6)*mUV^-2*rat(-1,2)
            + uvid(3,1,n1,n2,n3, - 1 + n4,n5,n6)*mUV^-2*rat(-1,1)
            + uvid(3,1,n1,n2,n3,n4,n5, - 1 + n6)*mUV^-2*rat(1,2)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{>1}) =
            + uvid(3,1,1 + n1, - 1 + n2,1 + n3,n4,n5, - 1 + n6)*rat(4*n3,3*n6 - 3)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,1 + n5, - 1 + n6)*rat( - n5,n6 - 1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - 5*n2,3*n6 - 3)
            + uvid(3,1,1 + n1,1 + n2,n3, - 1 + n4,n5, - 1 + n6)*rat( - n2,n6 - 1)
            + uvid(3,1,1 + n1,1 + n2,n3,n4, - 1 + n5, - 1 + n6)*rat(5*n2,3*n6 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,1 + n5, - 1 + n6)*rat(n5,n6 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,n2,1 + n3,n4, - 1 + n5, - 1 + n6)*rat( - 4*n3,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,1 + n5, - 1 + n6)*rat( - n5,n6 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,n5,n6)*rat(5,3)
            + uvid(3,1,1 + n1,n2,n3,1 + n4, - 1 + n5, - 1 + n6)*rat(4*n4,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,1 + n4,n5, - 2 + n6)*rat( - 4*n4,3*n6 - 3)
            + uvid(3,1,1 + n1,n2,n3,n4, - 1 + n5,n6)*rat(-5,3)
            + uvid(3,1,1 + n1,n2,n3,n4,1 + n5, - 2 + n6)*rat(n5,n6 - 1)
            + uvid(3,1,2 + n1, - 1 + n2,n3,n4,n5, - 1 + n6)*rat( - n1 - 1,3*n6 - 3)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - n1 - 1,3*n6 - 3)
            + uvid(3,1,2 + n1,n2,n3, - 1 + n4,n5, - 1 + n6)*rat(n1 + 1,3*n6 - 3)
            + uvid(3,1,2 + n1,n2,n3,n4,n5, - 2 + n6)*rat(n1 + 1,3*n6 - 3)
            + uvid(3,1,n1,1 + n2,n3,n4,n5, - 1 + n6)*rat(n2,n6 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>0}, n3?{>0}, n4?{>1}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2,1 + n3, - 1 + n4,n5,n6)*rat( - 2*n3,3*n4 - 3)
            + uvid(3,1,1 + n1, - 1 + n2,n3, - 1 + n4,1 + n5,n6)*rat(n5,n4 - 1)
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*rat(n2,3*n4 - 3)
            + uvid(3,1,1 + n1,1 + n2,n3, - 2 + n4,n5,n6)*rat( - n2,n4 - 1)
            + uvid(3,1,1 + n1,1 + n2,n3, - 1 + n4, - 1 + n5,n6)*rat( - n2,3*n4 - 3)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,1 + n5,n6)*rat( - n5,n4 - 1)
            + uvid(3,1,1 + n1,n2,1 + n3, - 1 + n4, - 1 + n5,n6)*rat(2*n3,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 2 + n4,1 + n5,n6)*rat( - n5,n4 - 1)
            + uvid(3,1,1 + n1,n2,n3, - 2 + n4,n5,1 + n6)*rat(2*n6,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4, - 1 + n5,1 + n6)*rat( - 2*n6,3*n4 - 3)
            + uvid(3,1,1 + n1,n2,n3, - 1 + n4,1 + n5, - 1 + n6)*rat(n5,n4 - 1)
            + uvid(3,1,1 + n1,n2,n3,n4, - 1 + n5,n6)*rat(1,3)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat(-1,3)
            + uvid(3,1,2 + n1, - 1 + n2,n3, - 1 + n4,n5,n6)*rat( - 4*n1 - 4,3*n4 - 3)
            + uvid(3,1,2 + n1,n2, - 1 + n3, - 1 + n4,n5,n6)*rat(2*n1 + 2,3*n4 - 3)
            + uvid(3,1,2 + n1,n2,n3, - 2 + n4,n5,n6)*rat(4*n1 + 4,3*n4 - 3)
            + uvid(3,1,2 + n1,n2,n3, - 1 + n4,n5, - 1 + n6)*rat( - 2*n1 - 2,3*n4 - 3)
            + uvid(3,1,n1,1 + n2,n3, - 1 + n4,n5,n6)*rat(n2,n4 - 1)
            ;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>0}, n3?{>1}, n4?{>0}, n5?{>0}, n6?{>0}) =
            + uvid(3,1,1 + n1, - 1 + n2,n3,n4,n5,n6)*rat(1,1)
            + uvid(3,1,1 + n1,1 + n2, - 2 + n3,n4,n5,n6)*rat( - n2,n3 - 1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3, - 1 + n4,n5,n6)*rat( - n2,n3 - 1)
            + uvid(3,1,1 + n1,1 + n2, - 1 + n3,n4, - 1 + n5,n6)*rat(n2,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,1 + n5,n6)*rat( - 2*n5,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3, - 1 + n4,n5,1 + n6)*rat(2*n6,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,1 + n4, - 1 + n5,n6)*rat(2*n4,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,1 + n4,n5, - 1 + n6)*rat( - 2*n4,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4, - 1 + n5,1 + n6)*rat( - 2*n6,n3 - 1)
            + uvid(3,1,1 + n1,n2, - 1 + n3,n4,1 + n5, - 1 + n6)*rat(2*n5,n3 - 1)
            + uvid(3,1,1 + n1,n2,n3,n4, - 1 + n5,n6)*rat(-1,1)
            + uvid(3,1,1 + n1,n2,n3,n4,n5, - 1 + n6)*rat(1,1)
            + uvid(3,1,2 + n1, - 1 + n2, - 1 + n3,n4,n5,n6)*rat( - n1 - 1,n3 - 1)
            + uvid(3,1,2 + n1,n2, - 2 + n3,n4,n5,n6)*rat(n1 + 1,n3 - 1)
            + uvid(3,1,2 + n1,n2, - 1 + n3, - 1 + n4,n5,n6)*rat(n1 + 1,n3 - 1)
            + uvid(3,1,2 + n1,n2, - 1 + n3,n4,n5, - 1 + n6)*rat( - n1 - 1,n3 - 1)
            + uvid(3,1,n1,1 + n2, - 1 + n3,n4,n5,n6)*rat(n2,n3 - 1)
            ;

        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{<1}, n6?{>0}) = uvid(3,1,n5, 0, n4, n1, n3, n6)*(uvid(3,1,-1, 0, 0, 0, 0, 0)- uvid(3,1,0, -1, 0, 0, 0, 0) + uvid(3,1,0, 0, -1, 0, 0, 0) + uvid(3,1,0, 0,0, -1, 0, 0) + uvid(3,1,0, 0, 0, 0, -1, 0) - uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0, 0, 0, 0, 0, 0))^-n2;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{>0}, n6?{<1}) = uvid(3,1,0, n6, n4, n2, n5, n3)*(-uvid(3,1,-1, 0, 0, 0, 0, 0) + uvid(3,1,0, -1, 0, 0, 0, 0) + uvid(3,1,0,0, -1, 0, 0, 0) + uvid(3,1,0, 0, 0, -1, 0, 0) - uvid(3,1,0, 0, 0, 0, -1, 0)+ uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0, 0, 0, 0, 0, 0))^-n1;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{<1}, n4?{<1}, n5?{>0}, n6?{>0}) = uvid(3,1,0, n4, n6, n3, n5, n2)*(-uvid(3,1,-1, 0, 0, 0, 0, 0) +uvid(3,1,0, -1, 0, 0, 0, 0) + uvid(3,1,0, 0, -1, 0, 0, 0) + uvid(3,1,0, 0, 0,-1, 0, 0) - uvid(3,1,0, 0, 0, 0, -1, 0) + uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0, 0, 0, 0, 0, 0))^-n1;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{<1}, n4?{>0}, n5?{<1}, n6?{>0}) = uvid(3,1,n5, n1, n4, 0, n2, n6)*(uvid(3,1,-1, 0, 0, 0, 0, 0) + uvid(3,1,0, -1, 0, 0, 0, 0) - uvid(3,1,0, 0, -1, 0, 0,0) - uvid(3,1,0, 0, 0, -1, 0, 0) + uvid(3,1,0, 0, 0, 0, -1, 0) + uvid(3,1,0,0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0, 0, 0, 0, 0, 0))^-n3;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{<1}, n4?{>0}, n5?{>0}, n6?{<1}) = uvid(3,1,n6, n3, n5, n1, n2, n4);
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{<1}, n5?{<1}, n6?{>0}) = uvid(3,1,n5, n1, n3, 0, n6, n2)*(uvid(3,1,-1, 0, 0, 0, 0,0) + uvid(3,1,0, -1, 0, 0, 0, 0) - uvid(3,1,0, 0, -1, 0, 0, 0) - uvid(3,1,0,0, 0, -1, 0, 0) + uvid(3,1,0, 0, 0, 0, -1, 0) + uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0, 0, 0, 0, 0, 0))^-n4;
        id ifmatch->end`i' uvid(3,1,n1?{<1}, n2?{>0}, n3?{>0}, n4?{>0}, n5?{<1}, n6?{<1}) = uvid(3,1,n5, 0, n3, n1, n4, n2)*(uvid(3,1,-1, 0, 0, 0, 0, 0) - uvid(3,1,0, -1, 0, 0, 0, 0) + uvid(3,1,0,0, -1, 0, 0, 0) + uvid(3,1,0, 0, 0, -1, 0, 0) + uvid(3,1,0, 0, 0, 0, -1, 0)- uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0, 0, 0, 0, 0, 0))^-n6;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{<1}, n3?{<1}, n4?{<1}, n5?{>0}, n6?{>0}) = uvid(3,1,n4, n3, n6, 0, n1, n5)*(uvid(3,1,-1, 0, 0, 0, 0, 0) + uvid(3,1,0, -1, 0, 0, 0, 0) - uvid(3,1,0, 0, -1, 0, 0, 0) - uvid(3,1,0, 0, 0, -1, 0,0) + uvid(3,1,0, 0, 0, 0, -1, 0) + uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0,0, 0, 0, 0, 0))^-n2;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{<1}, n3?{<1}, n4?{>0}, n5?{<1}, n6?{>0}) = uvid(3,1,n5, n3, n6, n2, n1, n4);
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{<1}, n3?{<1}, n4?{>0}, n5?{>0}, n6?{<1}) = uvid(3,1,n6, n2, n4, 0, n1, n5)*(uvid(3,1,-1, 0, 0, 0, 0, 0) + uvid(3,1,0, -1, 0, 0, 0, 0) -uvid(3,1,0, 0, -1, 0, 0, 0) - uvid(3,1,0, 0, 0, -1, 0, 0) + uvid(3,1,0, 0, 0,0, -1, 0) + uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0, 0, 0, 0, 0,0))^-n3;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{<1}, n3?{>0}, n4?{<1}, n5?{>0}, n6?{<1}) = uvid(3,1,n6, n2, n3, 0, n5, n1)*(uvid(3,1,-1, 0, 0, 0, 0, 0) + uvid(3,1,0, -1, 0, 0, 0, 0) - uvid(3,1,0, 0, -1, 0, 0, 0) - uvid(3,1,0, 0, 0, -1, 0,0) + uvid(3,1,0, 0, 0, 0, -1, 0) + uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0,0, 0, 0, 0, 0))^-n4;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{<1}, n3?{>0}, n4?{>0}, n5?{<1}, n6?{<1}) = uvid(3,1,0, n6, n3, n2, n1, n4)*(-uvid(3,1,-1, 0, 0, 0,0, 0) + uvid(3,1,0, -1, 0, 0, 0, 0) + uvid(3,1,0, 0, -1, 0, 0, 0) + uvid(3,1,0, 0, 0, -1, 0, 0) - uvid(3,1,0, 0, 0, 0, -1, 0) + uvid(3,1,0, 0, 0, 0, 0,-1) + mUV^2*uvid(3,1,0, 0, 0, 0, 0, 0))^-n5;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>0}, n3?{<1}, n4?{<1}, n5?{<1}, n6?{>0}) = uvid(3,1,0, n4, n2, n3, n1, n6)*(-uvid(3,1,-1, 0, 0, 0, 0, 0) + uvid(3,1,0, -1, 0, 0, 0, 0) +uvid(3,1,0, 0, -1, 0, 0, 0) + uvid(3,1,0, 0, 0, -1, 0, 0) - uvid(3,1,0, 0, 0,0, -1, 0) + uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0, 0, 0, 0, 0,0))^-n5;
        id ifmatch->end`i' uvid(3,1,n1?{>0}, n2?{>0}, n3?{<1}, n4?{<1}, n5?{>0}, n6?{<1}) = uvid(3,1,n4, n3, n2, 0, n5, n1)*(uvid(3,1,-1, 0, 0, 0, 0, 0) + uvid(3,1,0, -1, 0, 0, 0, 0) - uvid(3,1,0, 0, -1, 0, 0, 0) - uvid(3,1,0, 0, 0, -1, 0,0) + uvid(3,1,0, 0, 0, 0, -1, 0) + uvid(3,1,0, 0, 0, 0, 0, -1) + mUV^2*uvid(3,1,0,0, 0, 0, 0, 0))^-n6;

        goto realend;
        label end`i';
        redefine i "0";
        label realend;
* The two lines below wh
        id D = rat(4-2*ep, 1);
        Multiply replace_(D, 4-2*ep);

* some reduction rules have generated raising operators that have to be merged
        repeat id uvid(3,1,n1?,n2?,n3?,n4?,n5?,n6?)*uvid(3,1,nn1?,nn2?,nn3?,nn4?,nn5?,nn6?) = uvid(3,1,n1+nn1,n2+nn2,n3+nn3,n4+nn4,n5+nn5,n6+nn6);

        B+ uvid;
        .sort:reduction-3l-`nloop++';
        Keep brackets;
    #enddo
    .sort:reduction-3l-end;

    if (match(uvid(3,?a)));
        argument rat;
            if (count(ep,1) == -2);
                Print "Deep spurious pole detected";
                exit "Critical error";
            endif;
        endargument;
    endif;

#endprocedure

#procedure IntegrateUV()
* match the topology, working down in loop count
    #do i = {3L,2L,1L}
        #call IntegrateUV`i'()
    #enddo

    if (count(vxs,1,uvprop,1));
        Print "Unsubstituted UV graph: %t";
        exit "Critical error";
    endif;
#endprocedure

*#procedure Masters()
*    id uvid(1,1) = rat(ep,1) * mUV^2 * (rat(1,ep) + 1 + (1 + 1/12*pi^2)*rat(ep,1) + (1 + 1/12*pi^2 - 1/3*z3)*rat(ep^2,1) + alarmt1m3*rat(ep^3,1));
*    id uvid(2,1) = rat(ep^2,1) * (rat(1,ep^2) + pi^2/6 - 2/3*z3 * rat(ep, 1) + alarmt2m21*rat(ep^2,1));
*    id uvid(2,2) = 2 * ( -8453430797319808/21639331090712481 + 28375860251062315/42167665933747557 * rat(ep, 1) + alarmt2m22*rat(ep^2,1));
*
** mercedes
*    id uvid(3,1,1,1,1,1,1,1) = (2 * z3*rat(1,ep) + 1/55 * (-23459561604811217099/24730222849782180 + 330 * z3) + alarmt3m111111*rat(ep,1));
** 5-edge
*    id uvid(3,1,0,1,1,1,1,1) = mUV^2*(rat(1,ep^3) + rat(17/3,ep^2) + 62425469161433460/3513164752244341 * rat(1,ep)
*        + 693238079765918839/11778746848918154 + alarmt3m011111*rat(ep,1));
** banana
*    id uvid(3,1,0,1,1,1,0,1) = mUV^4*(rat(2,ep^3) + rat(23,3 * ep^2) + (35 + pi^2)*rat(1,2 * ep) + 1/12 * (275 + 23 * pi^2 - 24 * z3)
*        + 531237557421164097/8442111747388969*rat(ep,1) + alarmt3m011101*rat(ep^2,1));
** sunrise-bubble
*    id uvid(3,1,0,0,1,1,1,1) = mUV^4*(rat(3,2 * ep^3) + rat(6,ep^2) + 1287900757431154619/9648502930354515*rat(1,8 * ep)
*        + 650828837758949887/14793500421610002 + alarmt3m001111*rat(ep,1));
** triple-bubble
*    id uvid(3,1,0,0,1,0,1,1) = mUV^6*(rat(1,ep^3) + rat(3,ep^2) + (24 + pi^2)*rat(1,4 * ep) + 1/4 * (40 + 3 * pi^2 - 4 * z3)
*        + 1/480*(7200 + 720*pi^2 + 19*pi^4 - 1440*z3)*rat(ep,1) + alarmt3m001011*rat(ep^2,1));
*#endprocedure

#procedure Masters()

    id uvid(1,1) = rat(ep,1) * mUV^2 * (rat(1,ep) + 1 + (1 + 1/12*pi^2)*rat(ep,1) + (1 + 1/12*pi^2 - 1/3*z3)*rat(ep^2,1) + alarmt1m3*rat(ep^3,1));

    id uvid(2,1) = rat(ep^2,1) * (rat(1,ep^2) + pi^2/6 - 2/3*z3 * rat(ep, 1) + alarmt2m21*rat(ep^2,1));

    id uvid(2,2) = 2 * ( -2/9*sqrt3*cl2 + MASTER2m22EP * rat(ep, 1) + alarmt2m22*rat(ep^2,1));

* mercedes
    id uvid(3,1,1,1,1,1,1,1) = (2 * z3*rat(1,ep) + 1/55 * (MASTER3m111111FINPIECE + 330 * z3) + alarmt3m111111*rat(ep,1));

* 5-edge
    id uvid(3,1,0,1,1,1,1,1) = mUV^2*(rat(1,ep^3) + rat(17/3,ep^2) + (67/3 - 4*sqrt3*cl2 + pi^2/4) * rat(1,ep)
        + MASTER3m011111FIN + alarmt3m011111*rat(ep,1));

* banana
    id uvid(3,1,0,1,1,1,0,1) = mUV^4*(rat(2,ep^3) + rat(23,3 * ep^2) + (35 + pi^2)*rat(1,2 * ep) + 1/12 * (275 + 23 * pi^2 - 24 * z3)
        + MASTER3m011101EP*rat(ep,1) + alarmt3m011101*rat(ep^2,1));

* sunrise-bubble
    id uvid(3,1,0,0,1,1,1,1) = mUV^4*(rat(3,2 * ep^3) + rat(6,ep^2) + (-16*sqrt3*cl2 + 3*(44 + pi^2))*rat(1,8 * ep)
        + MASTER3m001111FIN + alarmt3m001111*rat(ep,1));

* triple-bubble
    id uvid(3,1,0,0,1,0,1,1) = mUV^6*(rat(1,ep^3) + rat(3,ep^2) + (24 + pi^2)*rat(1,4 * ep) + 1/4 * (40 + 3 * pi^2 - 4 * z3)
        + 1/480*(7200 + 720*pi^2 + 19*pi^4 - 1440*z3)*rat(ep,1) + alarmt3m001011*rat(ep^2,1));

#endprocedure

#procedure SubstituteMasters()
    Multiply replace_(D, 4-2*ep);
    id ep^n1? = rat(ep^n1,1);

    B+ uvid;
    .sort:masters-1;
    PolyRatFun rat(expand,ep,{`MAXPOLE'+`SELECTEDEPSILONORDER'+`SPURIOUSPOLE'});
    Keep brackets;

    .sort
    #call TruncateExpansion(`SPURIOUSPOLE')

    B+ uvid;
    .sort:normalization;
    Keep brackets;

    #call Masters()
    .sort:substitution;

    #call TruncateExpansion(0)
    .sort:truncation;
    PolyRatFun rat(expand,ep,{`MAXPOLE'+`SELECTEDEPSILONORDER'});

* fill in constants, they have not been substituted before to have exact cancellation
*    id pi = 30246273033735921/9627687726852338;
*   id z3 = 15408859284534249/12818743641862727;

#endprocedure
