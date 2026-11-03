#procedure SchemeZnostar5()
  #ifndef `TURBOTHRESHOLD'
    #define TURBOTHRESHOLD "4"
  #endif
*--#[ def :
* internals: p5,p6,p7
* externals: p0
* propagators: p1.p1,p2.p2,p3.p3,p4.p4,p5.p5,p6.p6,p7.p7,p8.p8,p0.p2
* replacements:
*   p3 = p0 + p5 + p7,
*   p1 = p0 + p6,
*   p9 = p0,
*   p8 = p5 - p6,
*   p2 = p0 + p6 + p7,
*   p4 = -p5 - p7,
*   p0.p0 = 1,
*   p9.p9 = 1,
* signs: +,+,+,+,*,+,+,+,-
*--#] def : 
*--#[ z2y :
#procedure z2y()
  #call Z2Y(9,8,5)
#endprocedure
*--#] z2y : 
*--#[ turbo :
#procedure turbo()
* depth = 3 (same as 4), threshold = 12
*--#[ n5 >= 2 :
#if `TURBOTHRESHOLD' >= 3
*--#[ nostar5_02_turbo-782643a7c6 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n9), (n2+n3+n4+n8)
* +:         n3,             n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,             n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n4)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n4)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n7+n3+2*n2+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar5_02_turbo-782643a7c6 : 
*--#[ nostar5_02_turbo-8bba96812f : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n6), [(n9)]
* decreases: (-n9), (n1+n2+n6+n7)
* +: n1,                         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3,         n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n6)
+`Z'(n1+1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n6)
+`Z'(n1+1,n2,n3,n4,n5,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(n1,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n8+2*n2+n1+2*ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
;
*--#] nostar5_02_turbo-8bba96812f : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ nostar5_02_turbo-40da96a34f : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1), [(n9)]
* decreases: (-n9), (n1+n3+n4+n5+n8)
* +:                     n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4, n5,         n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n6,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8-1,n9+1)*rat(-1,2)*rat(n6,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n7+n6+2*n5+n1+4*ep,-1+n1)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar5_02_turbo-40da96a34f : 
*--#[ nostar5_02_turbo-cbafaa47ec : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3), [(n9)]
* decreases: (-n9), (n1+n3+n5+n6+n7)
* +:             n4,             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3,     n5, n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4+1,n5-1,n6,n7,n8,n9+1)*rat(1,2)*rat(n4,-1+n3)
+`Z'(n1,n2,n3-1,n4+1,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n4,-1+n3)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n8+2*n5+n4+n3+4*ep,-1+n3)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
;
*--#] nostar5_02_turbo-cbafaa47ec : 
*--#[ nostar5_02_turbo-16995588ec : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 17+0=17 terms (0+0+6+9+2) ~ 3.067627 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n2), [(1+n9), (n9)]
* decreases: (-n9), (n1+n2+n3+n4+n5+n6+n7+n8)
* +: n1,                 n6, n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?{<=-2}) =
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n2)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n8+2*n6+n2+n1+2*ep,-1+n2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9+2)*rat(1,4)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7+1,n8,n9+2)*rat(-1,4)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6+1,n7,n8,n9+2)*rat(-1,4)*rat(n6,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6+1,n7,n8-1,n9+2)*rat(1,4)*rat(n6,-1+n2)
+`Z'(n1+1,n2-1,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(-1,4)*rat(n1,-1+n2)
+`Z'(n1+1,n2-1,n3,n4-1,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(n1,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5-1,n6,n7,n8,n9+2)*rat(-1,4)*rat(n1,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5,n6,n7,n8-1,n9+2)*rat(1,4)*rat(n1,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(-1,4)*rat(-3+n9+n7+n6+2*n5+n1+4*ep,-1+n2)
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
;
*--#] nostar5_02_turbo-16995588ec : 
*--#[ nostar5_02_turbo-054b6a0857 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 13+0=13 terms (0+0+9+0+4) ~ 3.079446 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n2), [(n9)]
* decreases: (-n9), (n1+n2+n3+n4+n5+n6+n7+n8)
* +:                         n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n2)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-3+n9+n8+n7+2*n5+4*ep,-1+n2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
;
*--#] nostar5_02_turbo-054b6a0857 : 
*--#[ nostar5_02_turbo-11da8b5af0 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 19+0=19 terms (0+0+11+6+2) ~ 3.624329 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n7), [(1+n9), (n9)]
* decreases: (-n9)
* +: n1, n2,         n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?{>=2},n8?pos_,n9?{<=-2}) =
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(n5+ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5+ep,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n2,-1+n7)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n7)
+`Z'(n1+1,n2,n3,n4,n5,n6,n7-1,n8,n9+2)*rat(1,4)*rat(n1,-1+n7)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(-3+n9+2*n8+n5+n2+3*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9+2)*rat(-1,4)*rat(n5+ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+2)*rat(1,4)*rat(n5+ep,-1+n7)
+`Z'(n1+1,n2,n3,n4,n5,n6-1,n7-1,n8,n9+2)*rat(-1,4)*rat(n1,-1+n7)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+2)*rat(-1,4)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+2)*rat(-1,4)*rat(-4+n9+n7+2*n6+n5+n1+3*ep,-1+n7)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7-1,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7-1,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n7)
;
*--#] nostar5_02_turbo-11da8b5af0 : 
*--#[ nostar5_02_turbo-2d456ff8f9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 19+0=19 terms (0+0+11+6+2) ~ 3.624329 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n8), [(1+n9), (n9)]
* decreases: (-n9)
* +:     n2, n3,     n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2,     n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?{>=2},n9?{<=-2}) =
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n2,-1+n8)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n5+ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5+ep,-1+n8)
+`Z'(n1,n2,n3+1,n4,n5,n6,n7,n8-1,n9+2)*rat(1,4)*rat(n3,-1+n8)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n8)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(-3+n9+2*n7+n5+n2+3*ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+2)*rat(-1,4)*rat(n5+ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+2)*rat(1,4)*rat(n5+ep,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7,n8-1,n9+2)*rat(-1,4)*rat(n3,-1+n8)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+2)*rat(-1,4)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+2)*rat(-1,4)*rat(-4+n9+n8+n5+2*n4+n3+3*ep,-1+n8)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8-1,n9+2)*rat(-1,4)*rat(1+n9,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8-1,n9+2)*rat(1,4)*rat(1+n9,-1+n8)
;
*--#] nostar5_02_turbo-2d456ff8f9 : 
*--#[ nostar5_02_turbo-3202f6bde5 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+2+6) ~ 3.667707 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n7), [(n9)]
* decreases: (-n9)
* +: n1,             n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,             n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?{>=2},n8?pos_,n9?neg_) =
+`Z'(n1+1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n1,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(n5+ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5+ep,-1+n7)
+`Z'(n1+1,n2,n3,n4,n5,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(n1,-1+n7)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-5+n7+2*n6+n5+n1+3*ep,-1+n7)
;
*--#] nostar5_02_turbo-3202f6bde5 : 
*--#[ nostar5_02_turbo-0324b87989 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+2+6) ~ 3.667707 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n8), [(n9)]
* decreases: (-n9)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,     n4,         n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?{>=2},n9?neg_) =
+`Z'(n1,n2,n3+1,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n3,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n5+ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5+ep,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n8)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(-5+n8+n5+2*n4+n3+3*ep,-1+n8)
;
*--#] nostar5_02_turbo-0324b87989 : 
#endif
#if `TURBOTHRESHOLD' >= 6
*--#[ nostar5_02_turbo-4e3af4a35e : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1), (-1+n6)
* decreases: (n1+n6)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,         n4,     n6
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6-1,n7+1,n8,n9)*rat(n7,-1+n1)*rat(n5+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6-1,n7+1,n8,n9)*rat(-n7,-1+n1)*rat(-1+n7+n5+2*ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-n5-ep,-1+n1)*rat(-1+n7+n5+2*ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n6)
;
*--#] nostar5_02_turbo-4e3af4a35e : 
*--#[ nostar5_02_turbo-3036f59470 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (n3+n4)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:         n3, n4,     n6
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?{>=2},n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2,n3-1,n4-1,n5+1,n6-1,n7,n8+1,n9)*rat(n8,-1+n3)*rat(n5+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(-n8,-1+n3)*rat(-1+n8+n5+2*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5+1,n6,n7,n8,n9)*rat(-n5-ep,-1+n3)*rat(-1+n8+n5+2*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n4)
;
*--#] nostar5_02_turbo-3036f59470 : 
#endif
#if `TURBOTHRESHOLD' >= 7
*--#[ nostar5_02_turbo-21952b836d : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,n9)*rat(-n5-ep,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9)*rat(n5+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-9+n9+2*n8+n7+n5+n4+2*n3+2*n2+5*ep,-1+n4)
;
*--#] nostar5_02_turbo-21952b836d : 
*--#[ nostar5_02_turbo-6089d927f7 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n6)
* decreases: (n1+n6)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,n9)*rat(n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(-n8,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,n9)*rat(-n5-ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9)*rat(n5+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9)*rat(-9+n9+n8+2*n7+n6+n5+2*n2+2*n1+5*ep,-1+n6)
;
*--#] nostar5_02_turbo-6089d927f7 : 
*--#[ nostar5_02_turbo-1b4f9e09a9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 9+0=9 terms (0+6+1+2) ~ 6.236841 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1)
* decreases: (n1+n6)
* +:     n2,                 n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6,     n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(-n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,n9)*rat(n7,-1+n1)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,n9)*rat(-n2,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,n9)*rat(n2,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)*rat(-9+n9+2*n8+n7+2*n6+2*n5+n2+n1+6*ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n1)
;
*--#] nostar5_02_turbo-1b4f9e09a9 : 
*--#[ nostar5_02_turbo-5f3583346c : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 9+0=9 terms (0+6+1+2) ~ 6.236841 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3)
* decreases: (n3+n4)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,n9)*rat(-n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(-n8,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,n9)*rat(n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)*rat(-9+n9+n8+2*n7+2*n5+2*n4+n3+n2+6*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
;
*--#] nostar5_02_turbo-5f3583346c : 
#endif
*--#] n5 >= 2 : 
*--#[ n5 <= 0 :
#if `TURBOTHRESHOLD' >= 3
*--#[ nostar5_02_turbo-98e72ab6ab : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 9+0=9 terms (0+0+4+2+3) ~ 2.348634 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n6)
* decreases: (-n5), (n6+n8)
* +: n1,             n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -:     n2, n3, n4,     n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1+1,n2-1,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(n1,-1+n6)
+`Z'(n1+1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9)*rat(-n1,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-1+n8-n7-n5-n4+n3+n2-ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9)
+`Z'(n1,n2,n3-1,n4,n5+1,n6-1,n7,n8,n9)*rat(8-n9-n8-n7-n6-n5-n4-n3-n2-n1-5*ep,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,n9)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n6)
+`Z'(n1,n2-1,n3,n4,n5+1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5+1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8-1,n9+1)*rat(-1,2)*rat(n9,-1+n6)
;
*--#] nostar5_02_turbo-98e72ab6ab : 
*--#[ nostar5_02_turbo-8e8be19b9f : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 9+0=9 terms (0+0+4+2+3) ~ 2.348634 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (-n5), (n4+n7)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -: n1, n2,     n4,     n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?neg0_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3+1,n4-1,n5+1,n6,n7,n8,n9)*rat(n3,-1+n4)
+`Z'(n1,n2,n3+1,n4-1,n5+1,n6,n7,n8-1,n9)*rat(-n3,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8,n9)*rat(-1-n8+n7-n6-n5+n2+n1-ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9)
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6,n7,n8,n9)*rat(8-n9-n8-n7-n6-n5-n4-n3-n2-n1-5*ep,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,n9)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n9,-1+n4)
;
*--#] nostar5_02_turbo-8e8be19b9f : 
*--#[ nostar5_02_turbo-77935dee09 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 9+0=9 terms (0+1+5+0+2+1) ~ 2.879719 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n7), (-1+n8), [(n3)]
* decreases: (-n5), (n6+n7+n8)
* +:     n2, n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -:             n4,     n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?{>=2},n8?{>=2},n9?neg0_) =
+`Z'(n1,n2,n3+1,n4,n5+1,n6,n7,n8-1,n9)*rat(n3,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5+1,n6,n7,n8-1,n9)*rat(-n3,-1+n8)
+`Z'(n1,n2+1,n3,n4,n5+1,n6,n7-1,n8-1,n9)*rat(n2,-1+n7)*rat(-1+n9+n2+ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9)*rat(1-n8-n5-2*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9)*rat(5-n8-n7-2*n5-2*n4-n3-5*ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+2,n6-1,n7-1,n8,n9)*rat(1+n5+ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+2,n6,n7-1,n8-1,n9)*rat(-1-n5-ep,-1+n7)*rat(-1+n8+n5+2*ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+2)*rat(1,4)*rat(n9,-1+n7)*rat(1+n9,-1+n8)
;
*--#] nostar5_02_turbo-77935dee09 : 
*--#[ nostar5_02_turbo-16c82d9859 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 9+0=9 terms (0+1+5+1+2) ~ 2.939690 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n7)
* decreases: (-n5), (n4+n7)
* +: n1, n2,         n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -:         n3, n4,     n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?{>=2},n8?pos_,n9?neg0_) =
+`Z'(n1+1,n2,n3,n4,n5+1,n6,n7-1,n8,n9)*rat(n1,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5+1,n6,n7-1,n8,n9)*rat(n2,-1+n7)
+`Z'(n1,n2+1,n3,n4,n5+1,n6,n7-1,n8-1,n9)*rat(-n2,-1+n7)
+`Z'(n1+1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9)*rat(-n1,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9)*rat(7-n9-2*n8-n7-2*n6-2*n5-n2-n1-6*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n9,-1+n7)
+`Z'(n1,n2,n3-1,n4,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n9,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n9,-1+n7)
;
*--#] nostar5_02_turbo-16c82d9859 : 
*--#[ nostar5_02_turbo-2e466704b7 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 9+0=9 terms (0+1+5+1+2) ~ 2.939690 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n8)
* decreases: (-n5), (n6+n8)
* +:     n2, n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -: n1,         n4,     n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?pos_,n8?{>=2},n9?neg0_) =
+`Z'(n1,n2,n3+1,n4,n5+1,n6,n7,n8-1,n9)*rat(n3,-1+n8)
+`Z'(n1-1,n2+1,n3,n4,n5+1,n6,n7,n8-1,n9)*rat(n2,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5+1,n6,n7,n8-1,n9)*rat(-n3,-1+n8)
+`Z'(n1,n2+1,n3,n4,n5+1,n6,n7-1,n8-1,n9)*rat(-n2,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9)*rat(7-n9-n8-2*n7-2*n5-2*n4-n3-n2-6*ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n9,-1+n8)
+`Z'(n1-1,n2,n3,n4,n5+1,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n9,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8-1,n9+1)*rat(-1,2)*rat(n9,-1+n8)
;
*--#] nostar5_02_turbo-2e466704b7 : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ nostar5_02_turbo-9f6b9f62ab : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+2+3+1) ~ 3.096633 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n2), [(n9)]
* decreases: (n2), (-n9)
* +: n1,     n3,                     n9
* 0: n1,     n3, n4, n5, n6, n7, n8
* -:     n2,     n4,     n6
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3+1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n3,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n1,-1+n2)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-13+n9+2*n8+2*n7+2*n6+2*n5+2*n4+n3+2*n2+n1+8*ep,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n2)
;
*--#] nostar5_02_turbo-9f6b9f62ab : 
#endif
*--#] n5 <= 0 : 
*--#[ n5 = 1 :
#if `TURBOTHRESHOLD' >= 3
*--#[ nostar5_02_turbo-6dd20d56b6 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n9), (n2+n3+n4+n8)
* +:         n3,             n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,             n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},1,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,1,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,1,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3+1,n4-1,1,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n4)
+`Z'(n1,n2,n3+1,n4-1,1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n4)
+`Z'(n1,n2-1,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,1,n6,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n7+n3+2*n2+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,1,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar5_02_turbo-6dd20d56b6 : 
*--#[ nostar5_02_turbo-cfcfb03191 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n6), [(n9)]
* decreases: (-n9), (n1+n2+n6+n7)
* +: n1,                         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3,         n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,1,n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3,n4,1,n6-1,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,1,n6-1,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n6)
+`Z'(n1+1,n2-1,n3,n4,1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n6)
+`Z'(n1+1,n2,n3,n4,1,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(n1,-1+n6)
+`Z'(n1-1,n2,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n8+2*n2+n1+2*ep,-1+n6)
+`Z'(n1,n2,n3,n4,1,n6,n7-1,n8,n9+1)*rat(-1,2)
;
*--#] nostar5_02_turbo-cfcfb03191 : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ nostar5_02_turbo-f181b693a8 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+2+3+1) ~ 3.096633 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n2), [(n9)]
* decreases: (n2), (-n9)
* +: n1,     n3,                     n9
* 0: n1,     n3, n4, n5, n6, n7, n8
* -:     n2,     n4,     n6
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,1,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3+1,n4,1,n6,n7,n8,n9+1)*rat(-1,2)*rat(n3,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,1,n6,n7,n8,n9+1)*rat(-1,2)*rat(n1,-1+n2)
+`Z'(n1,n2-1,n3+1,n4-1,1,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n2)
+`Z'(n1,n2-1,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)*rat(-11+n9+2*n8+2*n7+2*n6+2*n4+n3+2*n2+n1+8*ep,-1+n2)
+`Z'(n1,n2-1,n3,n4,1,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n2)
;
*--#] nostar5_02_turbo-f181b693a8 : 
*--#[ nostar5_02_turbo-4f1035cc74 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+2+6) ~ 3.667707 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n7), [(n9)]
* decreases: (-n9)
* +: n1,             n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,             n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,1,n6?pos_,n7?{>=2},n8?pos_,n9?neg_) =
+`Z'(n1+1,n2,n3,n4,1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n1,-1+n7)
+`Z'(n1,n2,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,2,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(1+ep,-1+n7)
+`Z'(n1,n2,n3,n4,2,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(1+ep,-1+n7)
+`Z'(n1+1,n2,n3,n4,1,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(n1,-1+n7)
+`Z'(n1,n2-1,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,1,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-4+n7+2*n6+n1+3*ep,-1+n7)
;
*--#] nostar5_02_turbo-4f1035cc74 : 
*--#[ nostar5_02_turbo-0f58529e23 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+2+6) ~ 3.667707 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n8), [(n9)]
* decreases: (-n9)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,     n4,         n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,1,n6?pos_,n7?pos_,n8?{>=2},n9?neg_) =
+`Z'(n1,n2,n3+1,n4,1,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n3,-1+n8)
+`Z'(n1,n2,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,2,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(1+ep,-1+n8)
+`Z'(n1,n2,n3,n4,2,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(1+ep,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n8)
+`Z'(n1,n2-1,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,1,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(-4+n8+2*n4+n3+3*ep,-1+n8)
;
*--#] nostar5_02_turbo-0f58529e23 : 
#endif
#if `TURBOTHRESHOLD' >= 6
*--#[ nostar5_02_turbo-5060334c88 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n1), (-1+n6)
* decreases: (n1+n6)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,         n4,     n6
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,1,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,2,n6-1,n7+1,n8,n9)*rat(n7,-1+n1)*rat(1+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,1,n6-1,n7+1,n8,n9)*rat(-n7,-1+n1)*rat(n7+2*ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,2,n6-1,n7,n8,n9)*rat(-1-ep,-1+n1)*rat(n7+2*ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,1,n6,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n1)
+`Z'(n1,n2,n3,n4,1,n6-1,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n6)
;
*--#] nostar5_02_turbo-5060334c88 : 
*--#[ nostar5_02_turbo-d0c3555946 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (n3+n4)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:         n3, n4,     n6
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?{>=2},1,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2,n3-1,n4-1,2,n6-1,n7,n8+1,n9)*rat(n8,-1+n3)*rat(1+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,1,n6,n7,n8+1,n9)*rat(-n8,-1+n3)*rat(n8+2*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,2,n6,n7,n8,n9)*rat(-1-ep,-1+n3)*rat(n8+2*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,1,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,1,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n4)
;
*--#] nostar5_02_turbo-d0c3555946 : 
#endif
#if `TURBOTHRESHOLD' >= 7
*--#[ nostar5_02_turbo-c4daddb357 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},1,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,1,n6,n7+1,n8,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,1,n6,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1,n2,n3,n4-1,2,n6-1,n7,n8,n9)*rat(-1-ep,-1+n4)
+`Z'(n1,n2,n3,n4-1,2,n6,n7,n8-1,n9)*rat(1+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,1,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,1,n6,n7,n8,n9)*rat(-8+n9+2*n8+n7+n4+2*n3+2*n2+5*ep,-1+n4)
;
*--#] nostar5_02_turbo-c4daddb357 : 
*--#[ nostar5_02_turbo-828c4ac422 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n6)
* decreases: (n1+n6)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,1,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4,1,n6-1,n7,n8+1,n9)*rat(n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,1,n6-1,n7,n8+1,n9)*rat(-n8,-1+n6)
+`Z'(n1,n2,n3,n4-1,2,n6-1,n7,n8,n9)*rat(-1-ep,-1+n6)
+`Z'(n1,n2,n3,n4,2,n6-1,n7-1,n8,n9)*rat(1+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,1,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4,1,n6-1,n7,n8,n9)*rat(-8+n9+n8+2*n7+n6+2*n2+2*n1+5*ep,-1+n6)
;
*--#] nostar5_02_turbo-828c4ac422 : 
#endif
*--#] n5 = 1 : 
#endprocedure
*--#] turbo : 
*--#[ reduce1346a :
#procedure reduce1346a()
*--#[ nostar5_02_x16 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n1), (-1+n6)
* decreases: (n1+n6)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,         n4,     n6
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6-1,n7+1,n8,n9)*rat(n7,-1+n1)*rat(n5+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6-1,n7+1,n8,n9)*rat(-n7,-1+n1)*rat(-1+n7+n5+2*ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-n5-ep,-1+n1)*rat(-1+n7+n5+2*ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n6)
;
*--#] nostar5_02_x16 : 
*--#[ nostar5_02_x34 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (n3+n4)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:         n3, n4,     n6
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?{>=2},n5?,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2,n3-1,n4-1,n5+1,n6-1,n7,n8+1,n9)*rat(n8,-1+n3)*rat(n5+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(-n8,-1+n3)*rat(-1+n8+n5+2*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5+1,n6,n7,n8,n9)*rat(-n5-ep,-1+n3)*rat(-1+n8+n5+2*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n4)
;
*--#] nostar5_02_x34 : 
#endprocedure
*--#] reduce1346a : 
*--#[ reduce9 :
#procedure reduce9()
* #call turbo()
* #call reduce1346a()
*--#[ nostar5_02_x29 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+2+3+1) ~ 3.096633 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n2), [(n9)]
* decreases: (n2), (-n9)
* +: n1,     n3,                     n9
* 0: n1,     n3, n4, n5, n6, n7, n8
* -:     n2,     n4,     n6
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3+1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n3,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n1,-1+n2)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-13+n9+2*n8+2*n7+2*n6+2*n5+2*n4+n3+2*n2+n1+8*ep,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n2)
;
*--#] nostar5_02_x29 : 
*--#[ nostar5_02_x9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+5+0+1) ~ 5.049357 [+,+,+,+,x,+,+,+,-]
* checked poles: (-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep), [(n9)]
* decreases: (-n9)
* +: n1,     n3,                     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,                 n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3+1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1,n2,n3+1,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1+1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n1,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1+1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n1,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n3+2*n2+n1+2*ep,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,2)*rat(1+n9,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
;
*--#] nostar5_02_x9 : 
#endprocedure
*--#] reduce9 : 
*--#[ reduce1346 :
#procedure reduce1346()
  #call reduce1346a()
*--#[ nostar5_02_p1 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,0) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1)
* decreases: (n1+n6)
* +:     n2,                 n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6,     n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,0)*rat(-n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,0)*rat(n7,-1+n1)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,0)*rat(-n2,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,0)*rat(n2,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,0)*rat(-9+2*n8+n7+2*n6+2*n5+n2+n1+6*ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,0)
;
*--#] nostar5_02_p1 : 
*--#[ nostar5_02_p3 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,0) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3)
* decreases: (n3+n4)
* +:     n2,                     n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,0)*rat(-n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,0)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,0)*rat(-n8,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,0)*rat(n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,0)*rat(-9+n8+2*n7+2*n5+2*n4+n3+n2+6*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,0)
;
*--#] nostar5_02_p3 : 
*--#[ nostar5_02_x4 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,0) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?,n6?pos_,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,0)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,0)*rat(n7,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,0)*rat(-n5-ep,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,0)*rat(n5+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,0)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,0)*rat(-9+2*n8+n7+n5+n4+2*n3+2*n2+5*ep,-1+n4)
;
*--#] nostar5_02_x4 : 
*--#[ nostar5_02_x6 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,0) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n6)
* decreases: (n1+n6)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?,n6?{>=2},n7?pos_,n8?pos_,0) =
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,0)*rat(n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,0)*rat(-n8,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,0)*rat(-n5-ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,0)*rat(n5+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,0)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,0)*rat(-9+n8+2*n7+n6+n5+2*n2+2*n1+5*ep,-1+n6)
;
*--#] nostar5_02_x6 : 
*--#[ nostar5_02_x1 : Z(n1,n2,n3,1,n5,1,n7,n8,0)
* Z(n1,n2,n3,1,n5+ep,1,n7,n8,0) gives 6+2=8 terms (0+6)+(0+2) ~ 6.020222 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n1), (-5+n8+n7+n3+n2+n1+2*ep)
* decreases: (n1)
* +:     n2,         n5,     n7
* 0:     n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3,                 n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,1,n5?,1,n7?pos_,n8?pos_,0) =
+`Z'(n1-2,n2,n3,1,n5,1,n7+1,n8,0)*rat(-n7,-1+n1)*rat(-2+n8+n7+n5+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1-1,n2-1,n3,1,n5,1,n7+1,n8,0)*rat(n7,-1+n1)*rat(-2+n8+n7+n5+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1-1,n2,n3,1,n5+1,1,n7,n8-1,0)*rat(n5+ep,-1+n1)*rat(2+n7+n5-n2-n1+ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1-1,n2+1,n3-1,1,n5,1,n7,n8,0)*rat(-n2,-1+n1)*rat(-3+n8+n5+n3+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1-1,n2+1,n3,1,n5,1,n7,n8-1,0)*rat(n2,-1+n1)*rat(-3+n8+n5+n3+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1-1,n2,n3,1,n5,1,n7,n8,0)*rat(1,-5+n8+n7+n3+n2+n1+2*ep)*rat(27-13*n8+2*n8^2-11*n7+3*n7*n8+n7^2-12*n5+4*n5*n8+2*n5*n7+n5^2-7*n3+2*n3*n8+2*n3*n7+2*n3*n5-6*n2+n2*n8+2*n2*n7+2*n2*n5+n2*n3-8*n1+n1*n8+n1*n7+n1*n5+n1*n3+n1*n2+n1^2-42*ep+12*ep*n8+8*ep*n7+8*
  ep*n5+6*ep*n3+6*ep*n2+5*ep*n1+15*ep^2,-1+n1)
+`Y'(n1-1,n2,n3,1,n5+1,0,n7,n8,0)*rat(-n5-ep,-1+n1)*rat(2+n7+n5-n2-n1+ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Y'(n1,n2,n3,1,n5,0,n7,n8,0)
;
*--#] nostar5_02_x1 : 
*--#[ nostar5_02_x3 : Z(n1,n2,n3,1,n5,1,n7,n8,0)
* Z(n1,n2,n3,1,n5+ep,1,n7,n8,0) gives 6+2=8 terms (0+6)+(0+2) ~ 6.020222 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n3), (-5+n8+n7+n3+n2+n1+2*ep)
* decreases: (n3)
* +:     n2,         n5,         n8
* 0: n1, n2,     n4, n5, n6, n7, n8, n9
* -: n1, n2, n3,             n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},1,n5?,1,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2+1,n3-1,1,n5,1,n7,n8,0)*rat(-n2,-1+n3)*rat(-3+n7+n5+n1+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2-1,n3-1,1,n5,1,n7,n8+1,0)*rat(n8,-1+n3)*rat(-2+n8+n7+n5+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2,n3-2,1,n5,1,n7,n8+1,0)*rat(-n8,-1+n3)*rat(-2+n8+n7+n5+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2,n3-1,1,n5+1,1,n7-1,n8,0)*rat(n5+ep,-1+n3)*rat(2+n8+n5-n3-n2+ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2+1,n3-1,1,n5,1,n7-1,n8,0)*rat(n2,-1+n3)*rat(-3+n7+n5+n1+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2,n3-1,1,n5,1,n7,n8,0)*rat(1,-5+n8+n7+n3+n2+n1+2*ep)*rat(27-11*n8+n8^2-13*n7+3*n7*n8+2*n7^2-12*n5+2*n5*n8+4*n5*n7+n5^2-8*n3+n3*n8+n3*n7+n3*n5+n3^2-6*n2+2*n2*n8+n2*n7+2*n2*n5+n2*n3-7*n1+2*n1*n8+2*n1*n7+2*n1*n5+n1*n3+n1*n2-42*ep+8*ep*n8+12*ep*n7+8*
  ep*n5+5*ep*n3+6*ep*n2+6*ep*n1+15*ep^2,-1+n3)
+`Y'(n1,n2,n3-1,0,n5+1,1,n7,n8,0)*rat(-n5-ep,-1+n3)*rat(2+n8+n5-n3-n2+ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Y'(n1,n2,n3,0,n5,1,n7,n8,0)
;
*--#] nostar5_02_x3 : 
#endprocedure
*--#] reduce1346 : 
*--#[ reduce78 :
#procedure reduce78()
*--#[ nostar5_02_x77 : Z(1,n2,1,1,n5,1,n7,n7,0)
* Z(1,n2,1,1,n5+ep,1,n7,n7,0) gives 5+6=11 terms (0+5)+(0+6) ~ 5.046025 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n7), (-3+2*n7+n2+2*ep)
* decreases: (n2+n7+n8)
* +:     n2,         n5
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,                 n7, n8
`id' `Z'(1,n2?pos_,1,1,n5?,1,n7?{>=2},n7?{>=2},0) =
+`Z'(1,n2,1,1,n5+1,1,n7-1,n7-1,0)*rat(-n5-ep,-1+n7)*rat(-6+3*n7+n5+2*n2+5*ep,-3+2*n7+n2+2*ep)
+`Z'(1,n2+1,1,1,n5,1,n7-1,n7-1,0)*rat(-n2,-1+n7)*rat(-6+3*n7+n5+2*n2+5*ep,-3+2*n7+n2+2*ep)
+`Z'(1,n2-1,1,1,n5,1,n7,n7,0)*rat(4-2*n7-n5-n2-3*ep,-3+2*n7+n2+2*ep)
+`Z'(1,n2,1,1,n5,1,n7-1,n7,0)*rat(2-n7-n2-2*ep,-3+2*n7+n2+2*ep)
+`Z'(1,n2,1,1,n5,1,n7,n7-1,0)*rat(5-3*n7-n5-n2-3*ep,-1+n7)*rat(-6+3*n7+n5+2*n2+5*ep,-3+2*n7+n2+2*ep)
+`Y'(0,n2+1,0,1,n5,1,n7,n7,0)*rat(-2*n2,-3+2*n7+n2+2*ep)
+`Y'(0,n2+1,1,1,n5,1,n7,n7-1,0)*rat(n2,-1+n7)*rat(-8+5*n7+n5+2*n2+5*ep,-3+2*n7+n2+2*ep)
+`Y'(1,n2,1,0,n5+1,1,n7,n7-1,0)*rat(n5+ep,-1+n7)*rat(-6+3*n7+n5+2*n2+5*ep,-3+2*n7+n2+2*ep)
+`Y'(0,n2,1,1,n5,1,n7,n7,0)*rat(-9+5*n7+2*n5+2*n2+8*ep,-3+2*n7+n2+2*ep)
+`Y'(1,n2,0,1,n5,1,n7,n7,0)*rat(-6+3*n7+n5+2*n2+5*ep,-3+2*n7+n2+2*ep)
+`Y'(1,n2,1,1,n5,0,n7,n7,0)
;
*--#] nostar5_02_x77 : 
*--#[ nostar5_02_x78 : Z(1,n2,1,1,n5,1,n7,n8,0)
* Z(1,n2,1,1,n5+ep,1,n7,n8,0) gives 3+7=10 terms (0+3)+(0+7) ~ 3.045463 [+,+,+,+,x,+,+,+,-]
* checked poles: (-3+n8+n7+n2+2*ep), (n8-n7)
* decreases: (n2+n7+n8)
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,                 n7, n8
if (match(`Z'(1,n2?,1,1,n5?,1,n7?!{n8?},n8?!{n7?},0)))
`id' `Z'(1,n2?pos_,1,1,n5?,1,n7?pos_,n8?pos_,0) =
+`Z'(1,n2-1,1,1,n5,1,n7,n8,0)*rat(4-n8-n7-n5-n2-3*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2,1,1,n5,1,n7-1,n8,0)*rat(-4+n8+n7+n5+n2+3*ep,n8-n7)*rat(-6+2*n8+n7+n5+2*n2+5*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2,1,1,n5,1,n7,n8-1,0)*rat(4-n8-n7-n5-n2-3*ep,n8-n7)*rat(-6+n8+2*n7+n5+2*n2+5*ep,-3+n8+n7+n2+2*ep)
+`Y'(0,n2+1,0,1,n5,1,n7,n8,0)*rat(-2*n2,-3+n8+n7+n2+2*ep)
+`Y'(0,n2+1,1,1,n5,1,n7,n8-1,0)*rat(2*n2,n8-n7)*rat(-6+2*n8+n7+n5+2*n2+5*ep,-3+n8+n7+n2+2*ep)
+`Y'(1,n2+1,0,1,n5,1,n7-1,n8,0)*rat(-2*n2,n8-n7)*rat(-6+n8+2*n7+n5+2*n2+5*ep,-3+n8+n7+n2+2*ep)
+`Y'(0,n2,1,1,n5,1,n7,n8,0)*rat(-3+2*n8+n5+3*ep,n8-n7)*rat(-6+2*n8+n7+n5+2*n2+5*ep,-3+n8+n7+n2+2*ep)
+`Y'(1,n2,0,1,n5,1,n7,n8,0)*rat(3-2*n7-n5-3*ep,n8-n7)*rat(-6+n8+2*n7+n5+2*n2+5*ep,-3+n8+n7+n2+2*ep)
+`Y'(1,n2,1,0,n5,1,n7,n8,0)*rat(6-n8-2*n7-n5-2*n2-5*ep,n8-n7)
+`Y'(1,n2,1,1,n5,0,n7,n8,0)*rat(-6+2*n8+n7+n5+2*n2+5*ep,n8-n7)
;
*--#] nostar5_02_x78 : 
#endprocedure
*--#] reduce78 : 
*--#[ reduce2 :
#procedure reduce2()
*--#[ nostar5_02_s2 : Z(1,n2,1,1,n5,1,1,1,0)
* Z(1,n2,1,1,n5+ep,1,1,1,0) gives 2+14=16 terms (0+2)+(0+5+9) ~ 2.054901 [+,+,+,+,*,+,+,+,-]
* checked poles: (-1+n2), (-1+n2+2*ep), (-1+n2+ep), (-2+2*n5+n2+6*ep), [(-2+n5+n2+3*ep), (-2+n5+n2+5*ep)]
* decreases: (n2)
* +:                 n5
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -:     n2
`id' `Z'(1,n2?{>=2},1,1,n5?,1,1,1,0) =
+`Z'(1,n2-1,1,1,n5+1,1,1,1,0)*rat(n5+ep,-1+n2)*rat(n5+2*ep,-1+n2+ep)*rat(n5+3*ep,-1+n2+2*ep)*rat(-3+n5+2*n2+5*ep,-2+2*n5+n2+6*ep)
+`Z'(1,n2-1,1,1,n5,1,1,1,0)*rat(-1-n5+n2-ep,-1+n2)*rat(-2+n5+n2+3*ep,-1+n2+ep)*rat(-2+n5+n2+4*ep,-1+n2+2*ep)*rat(-2+n5+n2+5*ep,-2+2*n5+n2+6*ep)
+`Y'(0,n2+1,0,1,n5,1,1,1,0)*rat(n2,-1+n2+ep)*rat(-1+n5+3*ep,-1+n2+2*ep)
+`Y'(0,n2+1,1,1,n5,1,1,0,0)*rat(-n2,-1+n2+ep)*rat(-1+n5+3*ep,-1+n2+2*ep)
+`Y'(0,n2,1,1,n5,1,1,1,0)*rat(-1+n5+3*ep,-1+n2+2*ep)
+`Y'(1,n2,1,1,n5,0,1,1,0)
+`Y'(1,n2,1,1,n5,1,0,1,0)*rat(-2+n5+n2+4*ep,-1+n2+2*ep)*rat(2-n5-n2-3*ep,-1+n2+ep)
+`Y'(0,n2-1,1,1,n5+1,1,1,1,0)*rat(1,2)*rat(n5+ep,-1+n2)*rat(n5+3*ep,-1+n2+ep)*rat(-2+n2+2*ep,-1+n2+2*ep)*rat(-3+n5+2*n2+5*ep,-2+2*n5+n2+6*ep)
+`Y'(0,n2,0,1,n5+1,1,1,1,0)*rat(-n5-ep,-1+n2+ep)*rat(n5+3*ep,-1+n2+2*ep)*rat(-3+n5+2*n2+5*ep,-2+2*n5+n2+6*ep)
+`Y'(1,n2-1,0,1,n5+1,1,1,1,0)*rat(1,2)*rat(n5+ep,-1+n2)*rat(n5+3*ep,-1+n2+ep)*rat(-2+n2+2*ep,-1+n2+2*ep)*rat(-3+n5+2*n2+5*ep,-2+2*n5+n2+6*ep)
+`Y'(1,n2-1,1,0,n5+1,1,1,1,0)*rat(-1,2)*rat(n5+ep,-1+n2)*rat(-2+n2+2*ep,-1+n2+ep)*rat(-3+n5+2*n2+5*ep,-1+n2+2*ep)
+`Y'(1,n2-1,1,1,n5+1,0,1,1,0)*rat(-1,2)*rat(n5+ep,-1+n2)*rat(-2+n2+2*ep,-1+n2+ep)*rat(-3+n5+2*n2+5*ep,-1+n2+2*ep)
+`Y'(1,n2-1,1,1,n5+1,1,0,1,0)*rat(1,2)*rat(n5+ep,-1+n2)*rat(-2+n2+2*ep,-1+n2+ep)*rat(-2+n5+n2+3*ep,-1+n2+2*ep)*rat(-3+n5+2*n2+5*ep,-2+2*n5+n2+6*ep)
+`Y'(1,n2-1,1,1,n5+1,1,1,0,0)*rat(1,2)*rat(n5+ep,-1+n2)*rat(-2+n2+2*ep,-1+n2+ep)*rat(-2+n5+n2+3*ep,-1+n2+2*ep)*rat(-3+n5+2*n2+5*ep,-2+2*n5+n2+6*ep)
+`Y'(1,n2,1,1,n5+1,0,0,1,0)*rat(n5+ep,-1+n2+ep)*rat(-3+n5+2*n2+5*ep,-1+n2+2*ep)
+`Y'(1,n2,1,1,n5+1,1,0,0,0)*rat(-n5-ep,-1+n2+ep)*rat(-2+n5+n2+3*ep,-1+n2+2*ep)*rat(-3+n5+2*n2+5*ep,-2+2*n5+n2+6*ep)
;
*--#] nostar5_02_s2 : 
#endprocedure
*--#] reduce2 : 
*--#[ reduce5 :
#procedure reduce5()
*--#[ nostar5_02_p5 : Z(1,1,1,1,n5,1,1,1,0)
* Z(1,1,1,1,n5+ep,1,1,1,0) gives 1+9=10 terms (0+1)+(0+9) ~ 1.070403 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n5+3*ep), (-1+n5+ep)
* decreases: (n5)
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -:                 n5
`id' `Z'(1,1,1,1,n5?{>=2},1,1,1,0) =
+`Z'(1,1,1,1,n5-1,1,1,1,0)*rat(1-n5-5*ep,-1+n5+ep)
+`Y'(1,1,1,0,n5+1,0,1,1,0)*rat(-2,1)*rat(n5+ep,-1+n5+3*ep)
+`Y'(1,1,1,1,n5+1,0,0,1,0)*rat(-2,1)*rat(n5+ep,-1+n5+ep)*rat(-1+n5+5*ep,-1+n5+3*ep)
+`Y'(1,1,1,1,n5+1,1,0,0,0)*rat(4,1)*rat(n5+ep,-1+n5+ep)
+`Y'(0,1,1,1,n5,1,1,1,0)*rat(1-n5-5*ep,-1+n5+ep)
+`Y'(1,0,1,1,n5,1,1,1,0)*rat(2,1)*rat(-1+n5+3*ep,-1+n5+ep)
+`Y'(1,1,1,0,n5,1,1,1,0)*rat(2,1)*rat(-1+n5+3*ep,-1+n5+ep)
+`Y'(1,1,1,1,n5,0,1,1,0)*rat(2,1)*rat(-1+n5+4*ep,-1+n5+3*ep)
+`Y'(1,1,1,1,n5,1,0,1,0)*rat(4*ep,-1+n5+ep)
+`Y'(1,1,1,1,n5,1,1,0,0)*rat(-1,1)
;
*--#] nostar5_02_p5 : 
*--#[ nostar5_02_m5 : Z(1,1,1,1,n5,1,1,1,0)
* Z(1,1,1,1,n5+ep,1,1,1,0) gives 1+9=10 terms (0+1)+(0+2+7) ~ 1.070403 [+,+,+,+,-,+,+,+,-]
* checked poles: (n5+3*ep), (n5+5*ep)
* decreases: (-n5)
* +:                 n5
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
`id' `Z'(1,1,1,1,n5?neg0_,1,1,1,0) =
+`Z'(1,1,1,1,n5+1,1,1,1,0)*rat(-n5-ep,n5+5*ep)
+`Y'(2,0,1,1,n5+1,1,1,1,0)*rat(1,n5+5*ep)*rat(2,1)*rat(n5+ep,n5+3*ep)
+`Y'(2,1,1,1,n5+1,1,0,1,0)*rat(-4,1)*rat(1,n5+5*ep)
+`Y'(2,1,1,1,n5+1,0,0,1,0)*rat(1,n5+3*ep)*rat(2,1)
+`Y'(0,1,1,1,n5+1,1,1,1,0)*rat(-4*ep,n5+5*ep)
+`Y'(1,0,1,1,n5+1,1,1,1,0)*rat(2*ep,n5+3*ep)
+`Y'(1,1,0,1,n5+1,1,1,1,0)*rat(-n5-ep,n5+5*ep)
+`Y'(1,1,1,0,n5+1,1,1,1,0)*rat(2,1)*rat(n5+3*ep,n5+5*ep)
+`Y'(1,1,1,1,n5+1,0,1,1,0)*rat(2*ep,n5+3*ep)*rat(n5+ep,n5+5*ep)
+`Y'(1,1,1,1,n5+1,1,0,1,0)*rat(3,1)
;
*--#] nostar5_02_m5 : 
#endprocedure
*--#] reduce5 : 
* Currently, the best performance for Znostar5(6,6,6,6,6,6,6,6,-5) is achieved
* by using only rule x29 as a turbo rule. (TU 5 Feb 2016)
* #call StdReduceProc(turbo,9,1346,78,2,5)
  #call StdReduceProc(9,1346,78,2,5)
#endprocedure
* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4
