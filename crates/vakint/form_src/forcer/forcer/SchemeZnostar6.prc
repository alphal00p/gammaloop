#procedure SchemeZnostar6()
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
* signs: +,+,+,+,+,*,+,+,-
*--#] def : 
*--#[ z2y :
#procedure z2y()
  #call Z2Y(9,8,6)
#endprocedure
*--#] z2y : 
*--#[ turbo :
#procedure turbo()
* depth = 3 (same as 4), threshold = 12
*--#[ n6 >= 2 :
#if `TURBOTHRESHOLD' >= 3
*--#[ nostar6_02_turbo-441e28532a : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n9), (n2+n3+n4+n8)
* +:         n3,             n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,             n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n4)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n4)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n7+n3+2*n2+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar6_02_turbo-441e28532a : 
*--#[ nostar6_02_turbo-dca6bc5dd8 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n5), [(n9)]
* decreases: (-n9), (n1+n2+n3+n5+n7+n8)
* +:                         n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3,     n5,     n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n5)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n5)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n5)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(-1,2)*rat(-3+n9+n8+n7+2*n2+2*ep,-1+n5)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar6_02_turbo-dca6bc5dd8 : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ nostar6_02_turbo-907a44e5ca : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1), [(n9)]
* decreases: (-n9), (n1+n3+n4+n5+n8)
* +:                     n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4, n5,         n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n6+ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8-1,n9+1)*rat(-1,2)*rat(n6+ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n7+n6+2*n5+n1+3*ep,-1+n1)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar6_02_turbo-907a44e5ca : 
*--#[ nostar6_02_turbo-445aa60ca9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3), [(n9)]
* decreases: (-n9), (n1+n3+n5+n6+n7)
* +:             n4,             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3,     n5, n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4+1,n5-1,n6,n7,n8,n9+1)*rat(1,2)*rat(n4,-1+n3)
+`Z'(n1,n2,n3-1,n4+1,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n4,-1+n3)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n8+2*n5+n4+n3+2*ep,-1+n3)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
;
*--#] nostar6_02_turbo-445aa60ca9 : 
*--#[ nostar6_02_turbo-63a8df048d : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 17+0=17 terms (0+0+6+9+2) ~ 3.067627 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n2), [(1+n9), (n9)]
* decreases: (-n9), (n1+n2+n3+n4+n5+n6+n7+n8)
* +: n1,                 n6, n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?{<=-2}) =
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n2)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n8+2*n6+n2+n1+4*ep,-1+n2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9+2)*rat(1,4)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7+1,n8,n9+2)*rat(-1,4)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6+1,n7,n8,n9+2)*rat(-1,4)*rat(n6+ep,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6+1,n7,n8-1,n9+2)*rat(1,4)*rat(n6+ep,-1+n2)
+`Z'(n1+1,n2-1,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(-1,4)*rat(n1,-1+n2)
+`Z'(n1+1,n2-1,n3,n4-1,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(n1,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5-1,n6,n7,n8,n9+2)*rat(-1,4)*rat(n1,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5,n6,n7,n8-1,n9+2)*rat(1,4)*rat(n1,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(-1,4)*rat(-3+n9+n7+n6+2*n5+n1+3*ep,-1+n2)
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
;
*--#] nostar6_02_turbo-63a8df048d : 
*--#[ nostar6_02_turbo-4b08b0da41 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 13+0=13 terms (0+0+9+0+4) ~ 3.079446 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n2), [(n9)]
* decreases: (-n9), (n1+n2+n3+n4+n5+n6+n7+n8)
* +:                         n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n2)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-3+n9+n8+n7+2*n5+2*ep,-1+n2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
;
*--#] nostar6_02_turbo-4b08b0da41 : 
*--#[ nostar6_02_turbo-86672079e9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 19+0=19 terms (0+0+11+6+2) ~ 3.624329 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n7), [(1+n9), (n9)]
* decreases: (-n9)
* +: n1, n2,         n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?{>=2},n8?pos_,n9?{<=-2}) =
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(n5,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n2,-1+n7)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n7)
+`Z'(n1+1,n2,n3,n4,n5,n6,n7-1,n8,n9+2)*rat(1,4)*rat(n1,-1+n7)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(-3+n9+2*n8+n5+n2+2*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9+2)*rat(-1,4)*rat(n5,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+2)*rat(1,4)*rat(n5,-1+n7)
+`Z'(n1+1,n2,n3,n4,n5,n6-1,n7-1,n8,n9+2)*rat(-1,4)*rat(n1,-1+n7)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+2)*rat(-1,4)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+2)*rat(-1,4)*rat(-4+n9+n7+2*n6+n5+n1+4*ep,-1+n7)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7-1,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7-1,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n7)
;
*--#] nostar6_02_turbo-86672079e9 : 
*--#[ nostar6_02_turbo-51b7fbf688 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 19+0=19 terms (0+0+11+6+2) ~ 3.624329 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n8), [(1+n9), (n9)]
* decreases: (-n9)
* +:     n2, n3,     n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2,     n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?{>=2},n9?{<=-2}) =
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n2,-1+n8)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n5,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5,-1+n8)
+`Z'(n1,n2,n3+1,n4,n5,n6,n7,n8-1,n9+2)*rat(1,4)*rat(n3,-1+n8)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n8)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(-3+n9+2*n7+n5+n2+2*ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+2)*rat(-1,4)*rat(n5,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+2)*rat(1,4)*rat(n5,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7,n8-1,n9+2)*rat(-1,4)*rat(n3,-1+n8)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+2)*rat(-1,4)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+2)*rat(-1,4)*rat(-4+n9+n8+n5+2*n4+n3+2*ep,-1+n8)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8-1,n9+2)*rat(-1,4)*rat(1+n9,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8-1,n9+2)*rat(1,4)*rat(1+n9,-1+n8)
;
*--#] nostar6_02_turbo-51b7fbf688 : 
*--#[ nostar6_02_turbo-361a1e584d : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 13+0=13 terms (0+1+9+1+2) ~ 3.645976 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n7), [(n9)]
* decreases: (-n9)
* +:     n2,         n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?{>=2},n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(n5,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n2,-1+n7)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n7)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(-3+n9+2*n8+n5+n2+2*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n7)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7-1,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7-1,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n7)
;
*--#] nostar6_02_turbo-361a1e584d : 
*--#[ nostar6_02_turbo-8553c6a1d1 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 8+0=8 terms (0+2+6) ~ 3.667707 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n8), [(n9)]
* decreases: (-n9)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,     n4,         n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?{>=2},n9?neg_) =
+`Z'(n1,n2,n3+1,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n3,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n5,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n8)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(-5+n8+n5+2*n4+n3+2*ep,-1+n8)
;
*--#] nostar6_02_turbo-8553c6a1d1 : 
#endif
#if `TURBOTHRESHOLD' >= 6
*--#[ nostar6_02_turbo-8ee784fdf4 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (n3+n4)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:         n3, n4,     n6
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?{>=2},n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2,n3-1,n4-1,n5+1,n6-1,n7,n8+1,n9)*rat(n5,-1+n3)*rat(n8,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(-n8,-1+n3)*rat(-1+n8+n5+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5+1,n6,n7,n8,n9)*rat(-n5,-1+n3)*rat(-1+n8+n5+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n4)
;
*--#] nostar6_02_turbo-8ee784fdf4 : 
#endif
#if `TURBOTHRESHOLD' >= 7
*--#[ nostar6_02_turbo-bc5d4cfa24 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,n9)*rat(-n5,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9)*rat(n5,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-9+n9+2*n8+n7+n5+n4+2*n3+2*n2+4*ep,-1+n4)
;
*--#] nostar6_02_turbo-bc5d4cfa24 : 
*--#[ nostar6_02_turbo-f02b2eff47 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+6+1+2) ~ 6.236841 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1)
* decreases: (n1+n6)
* +:     n2,                 n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6,     n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
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
*--#] nostar6_02_turbo-f02b2eff47 : 
*--#[ nostar6_02_turbo-61fe0ccf70 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+6+1+2) ~ 6.236841 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3)
* decreases: (n3+n4)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,n9)*rat(-n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(-n8,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,n9)*rat(n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)*rat(-9+n9+n8+2*n7+2*n5+2*n4+n3+n2+4*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
;
*--#] nostar6_02_turbo-61fe0ccf70 : 
#endif
*--#] n6 >= 2 : 
*--#[ n6 <= 0 :
#if `TURBOTHRESHOLD' >= 3
*--#[ nostar6_02_turbo-d95f625a58 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 5+0=5 terms (0+1+4) ~ 2.573482 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n5), (-1+n8)
* decreases: (-n6), (n5+n8)
* +:         n3, n4,     n6
* 0: n1, n2, n3, n4, n5,     n7, n8, n9
* -:                 n5,         n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?neg0_,n7?pos_,n8?{>=2},n9?neg0_) =
+`Z'(n1,n2,n3+1,n4+1,n5-1,n6+1,n7,n8-1,n9)*rat(n3,-1+n5)*rat(n4,-1+n8)
+`Z'(n1,n2,n3,n4+1,n5-1,n6+1,n7,n8-1,n9)*rat(-n4,-1+n5)*rat(-1+n4+n3+ep,-1+n8)
+`Z'(n1,n2,n3+1,n4,n5-1,n6+1,n7,n8-1,n9)*rat(-n3,-1+n5)*rat(-1+n4+n3+ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)*rat(-3+n8+n5+ep,-1+n5)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8-1,n9)*rat(-3+n8+n5+ep,-1+n8)
;
*--#] nostar6_02_turbo-d95f625a58 : 
*--#[ nostar6_02_turbo-4c7ef564ee : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+0+7+0+2) ~ 2.702762 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n1), [(n9)]
* decreases: (-n9), (n1+n3+n4+n5+n8)
* +:                     n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4, n5,         n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?neg0_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n7+n6+2*n5+n1+3*ep,-1+n1)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n6+ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8-1,n9+1)*rat(-1,2)*rat(n6+ep,-1+n1)
;
*--#] nostar6_02_turbo-4c7ef564ee : 
*--#[ nostar6_02_turbo-04b6c5ff26 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+0+7+1+1) ~ 2.760043 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n5), [(n9)]
* decreases: (-n9), (n1+n3+n4+n5+n8)
* +:                     n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5,         n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?neg0_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n5)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n7+n6+n5+2*n1+3*ep,-1+n5)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6+1,n7,n8,n9+1)*rat(-1,2)*rat(n6+ep,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n6+ep,-1+n5)
;
*--#] nostar6_02_turbo-04b6c5ff26 : 
*--#[ nostar6_02_turbo-39d4ad7162 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 6+0=6 terms (0+1+5) ~ 2.806142 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n5)
* decreases: (-n6), (n5+n8)
* +:             n4,     n6, n7
* 0: n1, n2, n3, n4, n5,     n7, n8, n9
* -: n1, n2, n3,     n5,         n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?neg0_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2,n3,n4+1,n5-1,n6+1,n7,n8,n9)*rat(-n4,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7+1,n8,n9)*rat(-n7,-1+n5)
+`Z'(n1,n2-1,n3,n4,n5-1,n6+1,n7+1,n8,n9)*rat(n7,-1+n5)
+`Z'(n1,n2,n3-1,n4+1,n5-1,n6+1,n7,n8,n9)*rat(n4,-1+n5)
+`Z'(n1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)*rat(-9+n9+2*n8+n7+n5+n4+2*n3+2*n2+4*ep,-1+n5)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8-1,n9)
;
*--#] nostar6_02_turbo-39d4ad7162 : 
*--#[ nostar6_02_turbo-b5e4992fcf : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n9), (n2+n3+n4+n8)
* +:         n3,             n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,             n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,n6?neg0_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n4)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n4)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n7+n3+2*n2+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar6_02_turbo-b5e4992fcf : 
*--#[ nostar6_02_turbo-5149eda8a3 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+1+5+1+1+1) ~ 2.924719 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n1), (-1+n8)
* decreases: (-n6)
* +:     n2,             n6, n7,     n9
* 0: n1, n2, n3, n4, n5,     n7, n8, n9
* -: n1,         n4, n5,         n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?neg0_,n7?pos_,n8?{>=2},n9?neg0_) =
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8,n9)
+`Z'(n1-1,n2,n3,n4-1,n5,n6+1,n7+1,n8,n9)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7+1,n8,n9)*rat(-n7,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6+1,n7,n8-1,n9)*rat(n2,-1+n1)*rat(-1+n9+n2+ep,-1+n8)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8,n9)*rat(5-n8-n7-2*n6-2*n5-n1-5*ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8-1,n9)*rat(1-n6-n1-2*ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6+2,n7,n8-1,n9)*rat(1+n6+ep,-1+n8)
+`Z'(n1-1,n2,n3,n4,n5,n6+2,n7,n8-1,n9)*rat(-1-n6-ep,-1+n1)*rat(-1+n6+n1+2*ep,-1+n8)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8-1,n9+2)*rat(1,4)*rat(n9,-1+n1)*rat(1+n9,-1+n8)
;
*--#] nostar6_02_turbo-5149eda8a3 : 
*--#[ nostar6_02_turbo-66567d400c : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+1+5+1+2) ~ 2.939690 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n1)
* decreases: (-n6)
* +:     n2,             n6, n7,     n9
* 0: n1, n2, n3, n4, n5,     n7, n8, n9
* -: n1,     n3, n4, n5,         n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?neg0_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8,n9)
+`Z'(n1-1,n2,n3,n4-1,n5,n6+1,n7+1,n8,n9)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7+1,n8,n9)*rat(-n7,-1+n1)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6+1,n7,n8,n9)*rat(n2,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6+1,n7,n8-1,n9)*rat(-n2,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8,n9)*rat(7-n9-2*n8-n7-2*n6-2*n5-n2-n1-6*ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6+1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n1)
;
*--#] nostar6_02_turbo-66567d400c : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ nostar6_02_turbo-872613192f : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 9+0=9 terms (0+1+6+1+1) ~ 3.109078 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n8)
* decreases: (n5+n8), (-n6+n8), (-n6-n9)
* +:     n2, n3,         n6,         n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,         n4, n5,     n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?neg0_,n7?pos_,n8?{>=2},n9?neg0_) =
+`Z'(n1,n2,n3+1,n4,n5,n6+1,n7,n8-1,n9)*rat(-n3,-1+n8)
+`Z'(n1-1,n2+1,n3,n4,n5,n6+1,n7,n8-1,n9)*rat(-n2,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5,n6+1,n7,n8-1,n9)*rat(n3,-1+n8)
+`Z'(n1,n2+1,n3,n4,n5,n6+1,n7-1,n8-1,n9)*rat(n2,-1+n8)
+`Z'(n1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n9,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8-1,n9)*rat(-9+n9+n8+2*n7+2*n5+2*n4+n3+n2+4*ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8-1,n9+1)*rat(-1,2)*rat(n9,-1+n8)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8-1,n9+1)*rat(-1,2)*rat(n9,-1+n8)
;
*--#] nostar6_02_turbo-872613192f : 
*--#[ nostar6_02_turbo-6f87c81505 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 10+0=10 terms (0+1+6+1+2) ~ 3.131594 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n2), [(n9)]
* decreases: (-n9)
* +:         n3,             n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3, n4, n5,         n8
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?neg0_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3+1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n3,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n7+2*n4+n3+n2+2*ep,-1+n2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
;
*--#] nostar6_02_turbo-6f87c81505 : 
#endif
#if `TURBOTHRESHOLD' >= 5
*--#[ nostar6_02_turbo-b909befe02 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 24+0=24 terms (0+0+21+0+3) ~ 4.603829 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (-n6-n9), (n3+n4-n6), (n2+n3-n6+n8), (n1+n2+n3+n4+n5+n7+n8)
* +:     n2, n3, n4,     n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4, n5,     n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?{>=2},n5?pos_,n6?neg0_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3-1,n4-1,n5,n6+1,n7+1,n8,n9)*rat(n2,-1+n3)*rat(n7,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,n6+1,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6+1,n7,n8,n9)*rat(n2,-1+n3)
+`Z'(n1,n2-1,n3,n4-1,n5,n6+1,n7+1,n8,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6+1,n7,n8,n9)*rat(-n3,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6+1,n7+1,n8,n9)*rat(n7,-1+n3)*rat(-1+n4-n2,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5-1,n6+1,n7+1,n8,n9)*rat(-n7,-1+n3)
+`Z'(n1,n2,n3-1,n4+1,n5-1,n6+1,n7,n8,n9)*rat(n4,-1+n3)
+`Z'(n1,n2,n3-1,n4+1,n5,n6+1,n7-1,n8,n9)*rat(-n4,-1+n3)
+`Z'(n1,n2,n3+1,n4-1,n5,n6+1,n7,n8-1,n9)*rat(n3,-1+n4)
+`Z'(n1,n2+1,n3-1,n4-1,n5,n6+1,n7,n8,n9)*rat(-n2,-1+n3)*rat(-3+n9+n7+n3+2*n2+2*ep,-1+n4)
+`Z'(n1,n2+1,n3-1,n4,n5,n6+1,n7-1,n8,n9)*rat(-n2,-1+n3)
+`Z'(n1,n2+1,n3,n4-1,n5,n6+1,n7,n8-1,n9)*rat(n2,-1+n4)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8,n9)
+`Z'(n1,n2-1,n3,n4,n5,n6+1,n7,n8,n9)*rat(-1,1)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6+1,n7,n8,n9)*rat(-n9-n7+n4-n2,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6+1,n7,n8,n9)*rat(4-n9-n7-n3-3*n2-2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7-1,n8,n9)*rat(-1,1)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8-1,n9)
+`Z'(n1-1,n2,n3-1,n4,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-2,n4,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4-1,n5,n6+1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
;
*--#] nostar6_02_turbo-b909befe02 : 
*--#[ nostar6_02_turbo-cab9ee97f1 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 19+0=19 terms (0+2+14+1+2) ~ 4.961732 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n4), (-1+n7)
* decreases: (-n6)
* +:     n2, n3, n4,     n6, n7,     n9
* 0: n1, n2, n3, n4, n5,     n7, n8, n9
* -: n1, n2, n3, n4, n5,     n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,n6?neg0_,n7?{>=2},n8?pos_,n9?neg0_) =
+`Z'(n1,n2,n3,n4+1,n5,n6+1,n7-1,n8,n9)*rat(n4,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8,n9)
+`Z'(n1-1,n2,n3,n4-1,n5,n6+1,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6+1,n7+1,n8,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6+1,n7,n8,n9)*rat(-n3,-1+n4)
+`Z'(n1,n2,n3-1,n4+1,n5,n6+1,n7-1,n8,n9)*rat(-n4,-1+n7)
+`Z'(n1,n2,n3+1,n4-1,n5,n6+1,n7,n8-1,n9)*rat(n3,-1+n4)
+`Z'(n1,n2+1,n3-1,n4,n5,n6+1,n7-1,n8,n9)*rat(-n2,-1+n7)
+`Z'(n1,n2+1,n3,n4,n5,n6+1,n7-1,n8-1,n9)*rat(n2,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8,n9)
+`Z'(n1,n2-1,n3,n4,n5,n6+1,n7,n8,n9)*rat(-1,1)
+`Z'(n1,n2,n3-1,n4,n5,n6+1,n7,n8,n9)*rat(-1,1)
+`Z'(n1,n2,n3,n4-1,n5,n6+1,n7,n8,n9)*rat(5-n9-n7-n4-n3-2*n2-2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7-1,n8,n9)*rat(5-n7-n4-2*n3-n2-2*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8-1,n9)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7-1,n8,n9+1)*rat(-1,2)*rat(n9,-1+n7)
+`Z'(n1,n2,n3-1,n4,n5,n6+1,n7-1,n8,n9+1)*rat(-1,2)*rat(n9,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5,n6+1,n7-1,n8,n9+1)*rat(1,2)*rat(n9,-1+n7)
;
*--#] nostar6_02_turbo-cab9ee97f1 : 
#endif
#if `TURBOTHRESHOLD' >= 6
*--#[ nostar6_02_turbo-3fa680e8aa : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 22+0=22 terms (0+2+15+2+3) ~ 5.122275 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n2), (-1+n3)
* decreases: (-n6-n9), (n2+n3-n6)
* +:     n2, n3, n4,     n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4, n5,     n7, n8
`id' `Z'(n1?pos_,n2?{>=2},n3?{>=2},n4?pos_,n5?pos_,n6?neg0_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3+1,n4,n5,n6+1,n7,n8,n9)*rat(n3,-1+n2)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8,n9)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6+1,n7,n8,n9)*rat(n2,-1+n3)
+`Z'(n1,n2-1,n3,n4-1,n5,n6+1,n7+1,n8,n9)*rat(-n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6+1,n7+1,n8,n9)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6+1,n7,n8,n9)*rat(-n3,-1+n2)
+`Z'(n1,n2,n3-1,n4+1,n5-1,n6+1,n7,n8,n9)*rat(n4,-1+n3)
+`Z'(n1,n2,n3-1,n4+1,n5,n6+1,n7-1,n8,n9)*rat(-n4,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6+1,n7-1,n8,n9)*rat(-n2,-1+n3)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8,n9)
+`Z'(n1,n2-1,n3,n4,n5,n6+1,n7,n8,n9)*rat(5-n7-2*n4-n3-n2-2*ep,-1+n2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6+1,n7,n8,n9)*rat(5-n9-2*n7-n4-n3-n2-2*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6+1,n7,n8,n9)*rat(-1,1)
+`Z'(n1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7-1,n8,n9)*rat(-1,1)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8-1,n9)
+`Z'(n1,n2-1,n3,n4,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n2)
+`Z'(n1,n2,n3-1,n4,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
+`Z'(n1,n2-1,n3-1,n4,n5,n6+1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n2)
;
*--#] nostar6_02_turbo-3fa680e8aa : 
*--#[ nostar6_02_turbo-ffa5d04433 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 22+0=22 terms (0+3+13+2+4) ~ 5.532391 [+,+,+,+,+,-,+,+,-]
* checked poles: (-1+n2), (-1+n7)
* decreases: (-n6)
* +:     n2, n3, n4,     n6, n7,     n9
* 0: n1, n2, n3, n4, n5,     n7, n8, n9
* -: n1, n2, n3, n4, n5,     n7, n8
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?neg0_,n7?{>=2},n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3+1,n4,n5,n6+1,n7,n8,n9)*rat(n3,-1+n2)
+`Z'(n1,n2,n3,n4+1,n5,n6+1,n7-1,n8,n9)*rat(n4,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8,n9)
+`Z'(n1,n2-1,n3,n4-1,n5,n6+1,n7+1,n8,n9)*rat(-n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6+1,n7+1,n8,n9)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6+1,n7,n8,n9)*rat(-n3,-1+n2)
+`Z'(n1,n2,n3-1,n4+1,n5,n6+1,n7-1,n8,n9)*rat(-n4,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6+1,n7-1,n8,n9)*rat(-n2,-1+n7)
+`Z'(n1,n2+1,n3,n4,n5,n6+1,n7-1,n8-1,n9)*rat(n2,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8,n9)
+`Z'(n1,n2-1,n3,n4,n5,n6+1,n7,n8,n9)*rat(5-n7-2*n4-n3-n2-2*ep,-1+n2)
+`Z'(n1,n2,n3-1,n4,n5,n6+1,n7,n8,n9)*rat(-1,1)
+`Z'(n1,n2,n3,n4-1,n5,n6+1,n7,n8,n9)*rat(-1,1)
+`Z'(n1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7-1,n8,n9)*rat(5-n7-n4-2*n3-n2-2*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7,n8-1,n9)
+`Z'(n1,n2-1,n3,n4,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n2)
+`Z'(n1,n2,n3,n4,n5,n6+1,n7-1,n8,n9+1)*rat(-1,2)*rat(n9,-1+n7)
+`Z'(n1,n2-1,n3-1,n4,n5,n6+1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n2)
+`Z'(n1,n2,n3-1,n4,n5,n6+1,n7-1,n8,n9+1)*rat(-1,2)*rat(n9,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5,n6+1,n7-1,n8,n9+1)*rat(1,2)*rat(n9,-1+n7)
;
*--#] nostar6_02_turbo-ffa5d04433 : 
#endif
*--#] n6 <= 0 : 
*--#[ n6 = 1 :
#if `TURBOTHRESHOLD' >= 3
*--#[ nostar6_02_turbo-90093e9d67 : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n9), (n2+n3+n4+n8)
* +:         n3,             n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,             n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,1,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,1,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,1,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3+1,n4-1,n5,1,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n4)
+`Z'(n1,n2,n3+1,n4-1,n5,1,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n4)
+`Z'(n1,n2-1,n3,n4,n5,1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n7+n3+2*n2+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5,1,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar6_02_turbo-90093e9d67 : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ nostar6_02_turbo-851615c6d6 : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n1), [(n9)]
* decreases: (-n9), (n1+n3+n4+n5+n8)
* +:                     n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4, n5,         n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,1,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,1,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,1,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,2,n7,n8,n9+1)*rat(1,2)*rat(1+ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,2,n7,n8-1,n9+1)*rat(-1,2)*rat(1+ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,1,n7,n8,n9+1)*rat(1,2)*rat(-4+n7+2*n5+n1+3*ep,-1+n1)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,1,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar6_02_turbo-851615c6d6 : 
*--#[ nostar6_02_turbo-99d8cbc4cf : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n5), [(n9)]
* decreases: (-n9), (n1+n2+n3+n5+n7+n8)
* +:                         n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3,     n5,     n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},1,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4,n5-1,1,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n5)
+`Z'(n1,n2-1,n3,n4,n5-1,1,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n5)
+`Z'(n1,n2-1,n3,n4,n5-1,1,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n5)
+`Z'(n1,n2,n3-1,n4,n5-1,1,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5-1,1,n7,n8,n9+1)*rat(-1,2)*rat(-3+n9+n8+n7+2*n2+2*ep,-1+n5)
+`Z'(n1,n2,n3,n4,n5,1,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,1,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] nostar6_02_turbo-99d8cbc4cf : 
*--#[ nostar6_02_turbo-b3e81464d3 : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 10+0=10 terms (0+1+6+1+2) ~ 3.131594 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n2), [(n9)]
* decreases: (-n9)
* +:         n3,             n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3, n4, n5,         n8
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,1,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3+1,n4,n5,1,n7,n8,n9+1)*rat(-1,2)*rat(n3,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5,1,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,1,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3+1,n4-1,n5,1,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,1,n7,n8,n9+1)*rat(1,2)*rat(-5+n7+2*n4+n3+n2+2*ep,-1+n2)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,1,n7,n8-1,n9+1)*rat(-1,2)
+`Z'(n1,n2-1,n3,n4,n5,1,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3-1,n4,n5,1,n7,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5,1,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
;
*--#] nostar6_02_turbo-b3e81464d3 : 
*--#[ nostar6_02_turbo-07b59a339d : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 8+0=8 terms (0+2+6) ~ 3.667707 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n8), [(n9)]
* decreases: (-n9)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,     n4,         n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,1,n7?pos_,n8?{>=2},n9?neg_) =
+`Z'(n1,n2,n3+1,n4,n5,1,n7,n8-1,n9+1)*rat(1,2)*rat(n3,-1+n8)
+`Z'(n1,n2,n3,n4,n5,1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5+1,1,n7,n8-1,n9+1)*rat(-1,2)*rat(n5,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,1,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5,1,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n8)
+`Z'(n1,n2-1,n3,n4,n5,1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,1,n7,n8-1,n9+1)*rat(-1,2)*rat(-5+n8+n5+2*n4+n3+2*ep,-1+n8)
;
*--#] nostar6_02_turbo-07b59a339d : 
#endif
#if `TURBOTHRESHOLD' >= 8
*--#[ nostar6_02_turbo-09e24b754c : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 10+0=10 terms (0+7+1+2) ~ 7.203038 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (n3+n4)
* +:     n2,                 n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5,         n8
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?{>=2},n5?pos_,1,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3-1,n4-1,n5,1,n7+1,n8,n9)*rat(n2,-1+n3)*rat(n7,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,1,n7+1,n8,n9)*rat(n7,-1+n3)*rat(-1+n4-n2,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5-1,1,n7+1,n8,n9)*rat(-n7,-1+n3)
+`Z'(n1,n2+1,n3-1,n4-1,n5,1,n7,n8,n9)*rat(-n2,-1+n3)*rat(-3+n9+n7+n3+2*n2+2*ep,-1+n4)
+`Z'(n1,n2+1,n3,n4-1,n5,1,n7,n8-1,n9)*rat(n2,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9)*rat(-5+n7+2*n4+n3+2*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9)*rat(-1+n4-n2,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-2,n4,n5,1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4-1,n5,1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
;
*--#] nostar6_02_turbo-09e24b754c : 
#endif
#if `TURBOTHRESHOLD' >= 11
*--#[ nostar6_02_turbo-de773b0580 : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 13+0=13 terms (0+10+1+2) ~ 10.144817 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n3), (-1+n5)
* decreases: (n3+n4)
* +:     n2,                 n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5,     n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?{>=2},1,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3-1,n4,n5-1,1,n7+1,n8,n9)*rat(n7,-1+n5)*rat(-n2,-1+n3)
+`Z'(n1,n2+1,n3-2,n4,n5-1,1,n7,n8+1,n9)*rat(n8,-1+n5)*rat(-n2,-1+n3)
+`Z'(n1-1,n2+1,n3-1,n4,n5,1,n7,n8,n9)*rat(-n2,-1+n3)
+`Z'(n1,n2,n3-1,n4-1,n5,1,n7+1,n8,n9)*rat(n7,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,1,n7,n8+1,n9)*rat(n2,-1+n3)*rat(n8,-1+n5)
+`Z'(n1,n2,n3-1,n4,n5-1,1,n7+1,n8,n9)*rat(-n7,-1+n3)*rat(-1+n5-n2,-1+n5)
+`Z'(n1,n2+1,n3-1,n4,n5-1,1,n7,n8,n9)*rat(n2,-1+n3)*rat(-2+n9+n8+n7+2*n2+2*ep,-1+n5)
+`Z'(n1,n2+1,n3-1,n4,n5,1,n7-1,n8,n9)*rat(n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9)*rat(-5+n7+2*n4+n3+n2+2*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-2,n4,n5,1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4-1,n5,1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
;
*--#] nostar6_02_turbo-de773b0580 : 
#endif
#if `TURBOTHRESHOLD' >= 12
*--#[ nostar6_02_turbo-4b761b1286 : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 14+0=14 terms (0+11+1+2) ~ 11.133033 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n1), (-1+n3)
* decreases: (n3+n4)
* +:     n2,             n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5,         n8
`id' `Z'(n1?{>=2},n2?pos_,n3?{>=2},n4?pos_,n5?pos_,1,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3-1,n4-1,n5,1,n7+1,n8,n9)*rat(n2,-1+n1)*rat(n7,-1+n3)
+`Z'(n1-1,n2+1,n3-1,n4,n5-1,1,n7+1,n8,n9)*rat(n7,-1+n3)*rat(-n2,-1+n1)
+`Z'(n1-1,n2+1,n3-1,n4,n5-1,2,n7,n8,n9)*rat(-n2,-1+n1)*rat(1+ep,-1+n3)
+`Z'(n1-1,n2+1,n3-1,n4,n5,2,n7,n8-1,n9)*rat(n2,-1+n1)*rat(1+ep,-1+n3)
+`Z'(n1-1,n2+1,n3-1,n4,n5,1,n7,n8,n9)*rat(-n2,-1+n1)*rat(-4+n7+2*n5+n1+3*ep,-1+n3)
+`Z'(n1,n2,n3-1,n4-1,n5,1,n7+1,n8,n9)*rat(n7,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,1,n7+1,n8,n9)*rat(-n7,-1+n3)
+`Z'(n1,n2+1,n3-1,n4-1,n5,1,n7,n8,n9)*rat(n2,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5-1,1,n7,n8,n9)*rat(-n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9)*rat(-5+n7+2*n4+n3+n2+2*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-2,n4,n5,1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4-1,n5,1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
;
*--#] nostar6_02_turbo-4b761b1286 : 
*--#[ nostar6_02_turbo-a3775befea : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 16+0=16 terms (0+11+2+3) ~ 11.243449 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n4), (-1+n8)
* decreases: (n3+n4)
* +:     n2,         n5,     n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,1,n7?pos_,n8?{>=2},n9?neg0_) =
+`Z'(n1-1,n2+1,n3,n4-1,n5,1,n7+1,n8-1,n9)*rat(n7,-1+n8)*rat(-n2,-1+n4)
+`Z'(n1,n2,n3,n4-2,n5+1,1,n7+1,n8-1,n9)*rat(n7,-1+n8)*rat(-n5,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5,0,n7+1,n8-1,n9+1)*rat(1,2)*rat(n7,-1+n4)*rat(n9,-1+n8)
+`Z'(n1-1,n2,n3,n4-1,n5,1,n7+1,n8,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,1,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5,1,n7+1,n8-1,n9)*rat(n7,-1+n4)*rat(-2+n9+2*n7+n5+n2+2*ep,-1+n8)
+`Z'(n1,n2,n3,n4-1,n5+1,1,n7,n8-1,n9)*rat(n5,-1+n4)*rat(n7,-1+n8)
+`Z'(n1,n2+1,n3-1,n4-1,n5,1,n7,n8,n9)*rat(n2,-1+n4)
+`Z'(n1,n2+1,n3,n4-1,n5,1,n7,n8-1,n9)*rat(-n2,-1+n4)*rat(-1+n8-n7,-1+n8)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9)*rat(-5+n7+n4+2*n3+n2+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5,1,n7+1,n8-1,n9+1)*rat(-1,2)*rat(n7,-1+n4)*rat(n9,-1+n8)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,1,n7+1,n8-1,n9+1)*rat(-1,2)*rat(n7,-1+n4)*rat(n9,-1+n8)
+`Z'(n1,n2,n3-1,n4-1,n5,1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n4)
+`Z'(n1,n2,n3,n4-2,n5,1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n4)
;
*--#] nostar6_02_turbo-a3775befea : 
*--#[ nostar6_02_turbo-f24f333033 : Z(n1,n2,n3,n4,n5,1,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,1+ep,n7,n8,n9) gives 17+0=17 terms (0+11+2+4) ~ 11.258302 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n1), (-1+n4)
* decreases: (n3+n4)
* +:     n2,             n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5,         n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?{>=2},n5?pos_,1,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,n5-1,2,n7+1,n8,n9)*rat(n7,-1+n1)*rat(1+ep,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,2,n7+1,n8-1,n9)*rat(-n7,-1+n1)*rat(1+ep,-1+n4)
+`Z'(n1-1,n2+1,n3-1,n4-1,n5,1,n7+1,n8,n9)*rat(n2,-1+n1)*rat(n7,-1+n4)
+`Z'(n1-1,n2+1,n3,n4-1,n5,1,n7+1,n8-1,n9)*rat(n7,-1+n4)*rat(-n2,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,1,n7+1,n8,n9)*rat(-n7,-1+n1)*rat(-4+n9+2*n8+n2+n1+3*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,1,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5,1,n7+1,n8-1,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2+1,n3-1,n4-1,n5,1,n7,n8,n9)*rat(n2,-1+n4)
+`Z'(n1,n2+1,n3,n4-1,n5,1,n7,n8-1,n9)*rat(-n2,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,1,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9)*rat(-5+n7+n4+2*n3+n2+2*ep,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,1,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n1)*rat(n9,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5,1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n4)
+`Z'(n1-1,n2,n3-1,n4-1,n5,1,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n1)*rat(n9,-1+n4)
+`Z'(n1-1,n2,n3,n4-2,n5,1,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n1)*rat(n9,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n4)
+`Z'(n1,n2,n3,n4-2,n5,1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n4)
;
*--#] nostar6_02_turbo-f24f333033 : 
#endif
*--#] n6 = 1 : 
#endprocedure
*--#] turbo : 
*--#[ reduce9 :
#procedure reduce9
  #call turbo()
*--#[ nostar6_02_x9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,n9) gives 6+0=6 terms (0+5+0+1) ~ 5.049357 [+,+,+,+,+,x,+,+,-]
* checked poles: (-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep), [(n9)]
* decreases: (-n9)
* +: n1,     n3,                     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,                 n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3+1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1,n2,n3+1,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1+1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n1,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1+1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n1,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n3+2*n2+n1+2*ep,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,2)*rat(1+n9,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep)
;
*--#] nostar6_02_x9 : 
#endprocedure
*--#] reduce9 : 
*--#[ reduce34 :
#procedure reduce34
* NOTE: x3 and x4 contain a term with a worse complexity.
*--#[ nostar6_02_x34 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,0) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (n3+n4)
* +:     n2,                 n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?{>=2},n5?pos_,n6?,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2+1,n3-1,n4-1,n5,n6,n7+1,n8,0)*rat(n2,-1+n3)*rat(n7,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7+1,n8,0)*rat(-n7,-1+n3)*rat(-1+n7+n2+ep,-1+n4)
+`Z'(n1,n2+1,n3-1,n4-1,n5,n6,n7,n8,0)*rat(-n2,-1+n3)*rat(-1+n7+n2+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,0)*rat(-3+n4+n3+ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,0)*rat(-3+n4+n3+ep,-1+n4)
;
*--#] nostar6_02_x34 : 
*--#[ nostar6_02_x3 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,0) gives 6+0=6 terms (1|0+5) [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n3)
* decreases: (n3+n4)
* +:     n2,                     n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?pos_,n6?,n7?pos_,n8?pos_,0) =
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,0)*rat(-n8,-1+n3)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,0)*rat(-n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,0)*rat(n8,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,0)*rat(n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,0)*rat(-9+n8+2*n7+2*n5+2*n4+n3+n2+4*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,0)
;
*--#] nostar6_02_x3 : 
*--#[ nostar6_02_x4 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,0) gives 6+0=6 terms (1|0+5) [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,n6?,n7?pos_,n8?pos_,0) =
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,0)*rat(-n5,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,0)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,0)*rat(n7,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,0)*rat(n5,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,0)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,0)*rat(-9+2*n8+n7+n5+n4+2*n3+2*n2+4*ep,-1+n4)
;
*--#] nostar6_02_x4 : 
#endprocedure
*--#] reduce34 : 
*--#[ reduce1 :
#procedure reduce1
*--#[ nostar6_02_p1 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5,n6+ep,n7,n8,0) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1)
* decreases: (n1+n6)
* +:     n2,                 n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6,     n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,0)*rat(-n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,0)*rat(n7,-1+n1)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,0)*rat(-n2,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,0)*rat(n2,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,0)*rat(-9+2*n8+n7+2*n6+2*n5+n2+n1+6*ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,0)
;
*--#] nostar6_02_p1 : 
*--#[ nostar6_02_x1 : Z(n1,n2,1,1,n5,n6,n7,n8,0)
* Z(n1,n2,1,1,n5,n6+ep,n7,n8,0) gives 4+6=10 terms (0+4)+(0+4+2) ~ 4.044187 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n1), (-4+n7+n5+n2+n1+2*ep), [(2*n8-n7+n6-n5-n1+3*ep)]
* decreases: (n1+n5+n7)
* +:                         n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,             n5,     n7
`id' `Z'(n1?{>=2},n2?pos_,1,1,n5?pos_,n6?,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2,1,1,n5-1,n6,n7+1,n8,0)*rat(n7,-1+n1)*rat(-8+n8+2*n7+n6+2*n5+n2+n1+5*ep,-4+n7+n5+n2+n1+2*ep)
+`Z'(n1-1,n2,1,1,n5,n6,n7,n8,0)*rat(-5+n8+n6+n5+n2+n1+3*ep,-1+n1)*rat(-8+n8+2*n7+n6+2*n5+n2+n1+5*ep,-4+n7+n5+n2+n1+2*ep)
+`Z'(n1,n2,1,1,n5-1,n6,n7,n8,0)*rat(-4+n8+n7+n5+n2+2*ep,-4+n7+n5+n2+n1+2*ep)
+`Z'(n1,n2,1,1,n5,n6,n7-1,n8,0)*rat(4-n8-n7-n5-n2-2*ep,-4+n7+n5+n2+n1+2*ep)
+`Y'(n1-1,n2,1,0,n5,n6,n7+1,n8,0)*rat(-n7,-1+n1)*rat(-10+n8+2*n7+n6+2*n5+n2+3*n1+5*ep,-4+n7+n5+n2+n1+2*ep)
+`Y'(n1,n2-1,1,0,n5,n6,n7+1,n8,0)*rat(2*n7,-4+n7+n5+n2+n1+2*ep)
+`Y'(n1,n2,0,1,n5,n6,n7,n8,0)
+`Y'(n1,n2,1,0,n5,n6,n7,n8,0)*rat(-2+n8+2*n2-n1+2*ep,-4+n7+n5+n2+n1+2*ep)
+`Y'(n1-1,n2,0,1,n5,n6,n7,n8,0)*rat(7-n8-n7-n6-n5-n2-n1-5*ep,-1+n1)*rat(-8+n8+2*n7+n6+2*n5+n2+n1+5*ep,-4+n7+n5+n2+n1+2*ep)
+`Y'(n1-1,n2,1,0,n5,n6,n7,n8,0)*rat(-7+n8+n7+n6+n5+n2+n1+5*ep,-1+n1)*rat(-8+n8+2*n7+n6+2*n5+n2+n1+5*ep,-4+n7+n5+n2+n1+2*ep)
;
*--#] nostar6_02_x1 : 
#endprocedure
*--#] reduce1 : 
*--#[ reduce58 :
#procedure reduce58
*--#[ nostar6_02_x5 : Z(1,n2,1,1,n5,n6,n7,n8,0)
* Z(1,n2,1,1,n5,n6+ep,n7,n8,0) gives 6+5=11 terms (0+6)+(0+5) ~ 6.041733 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n5), (-3+n7+n5+n2+2*ep)
* decreases: (n2+n5+n7+n8)
* +:                         n7, n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,         n5,     n7, n8
`id' `Z'(1,n2?pos_,1,1,n5?{>=2},n6?,n7?pos_,n8?pos_,0) =
+`Z'(1,n2-1,1,1,n5-1,n6,n7,n8+1,0)*rat(-n8,-1+n5)*rat(-8+2*n8+n7+n6+2*n5+2*n2+5*ep,-3+n7+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5-1,n6,n7+1,n8,0)*rat(-n7,-1+n5)*rat(-6+2*n8+n7+n6+2*n2+5*ep,-3+n7+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5,n6,n7,n8,0)*rat(2-n7-n2-2*ep,-3+n7+n5+n2+2*ep)
+`Z'(1,n2,1,1,n5-1,n6,n7,n8,0)*rat(-1,-3+n7+n5+n2+2*ep)*rat(35-16*n8+2*n8^2-13*n7+3*n7*n8+n7^2-5*n6+n6*n8+n6*n7-12*n5+2*n5*n8+3*n5*n7+n5*n6+n5^2-24*n2+6*n2*n8+4*n2*n7+2*n2*n6+4*n2*n5+4*n2^2-39*ep+9*ep*n8+7*ep*n7+2*ep*n6+7*ep*n5+14*ep*n2+10*ep^2,-1+n5)
+`Z'(1,n2,1,1,n5,n6,n7-1,n8,0)*rat(8-2*n8-n7-n6-2*n5-2*n2-5*ep,-3+n7+n5+n2+2*ep)
+`Z'(1,n2,1,1,n5,n6,n7,n8-1,0)
+`Y'(0,n2,1,0,n5,n6,n7+1,n8,0)*rat(-2*n7,-3+n7+n5+n2+2*ep)
+`Y'(0,n2,1,1,n5-1,n6,n7+1,n8,0)*rat(n7,-1+n5)*rat(-8+2*n8+n7+n6+2*n5+2*n2+5*ep,-3+n7+n5+n2+2*ep)
+`Y'(1,n2,0,1,n5-1,n6,n7,n8+1,0)*rat(n8,-1+n5)*rat(-8+2*n8+n7+n6+2*n5+2*n2+5*ep,-3+n7+n5+n2+2*ep)
+`Y'(0,n2,1,1,n5,n6,n7,n8,0)*rat(-11+2*n8+2*n7+2*n6+3*n5+2*n2+8*ep,-3+n7+n5+n2+2*ep)
+`Y'(1,n2,1,0,n5,n6,n7,n8,0)*rat(-10+2*n8+2*n7+n6+2*n5+3*n2+7*ep,-3+n7+n5+n2+2*ep)
;
*--#] nostar6_02_x5 : 
*--#[ nostar6_02_x8 : Z(1,n2,1,1,n5,n6,n7,n8,0)
* Z(1,n2,1,1,n5,n6+ep,n7,n8,0) gives 6+5=11 terms (0+6)+(0+5) ~ 6.041733 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n8), (-3+n8+n7+n2+2*ep)
* decreases: (n2+n5+n7+n8)
* +:     n2,         n5
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,         n5,     n7, n8
`id' `Z'(1,n2?pos_,1,1,n5?pos_,n6?,n7?pos_,n8?{>=2},0) =
+`Z'(1,n2,1,1,n5+1,n6,n7-1,n8-1,0)*rat(-n5,-1+n8)*rat(-8+2*n8+2*n7+n6+2*n5+n2+5*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2+1,1,1,n5,n6,n7-1,n8-1,0)*rat(-n2,-1+n8)*rat(-6+2*n7+n6+2*n5+n2+5*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2-1,1,1,n5,n6,n7,n8,0)*rat(8-2*n8-2*n7-n6-2*n5-n2-5*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2,1,1,n5-1,n6,n7,n8,0)
+`Z'(1,n2,1,1,n5,n6,n7-1,n8,0)*rat(2-n7-n2-2*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2,1,1,n5,n6,n7,n8-1,0)*rat(-1,-3+n8+n7+n2+2*ep)*rat(35-12*n8+n8^2-24*n7+4*n7*n8+4*n7^2-5*n6+n6*n8+2*n6*n7-16*n5+2*n5*n8+6*n5*n7+n5*n6+2*n5^2-13*n2+3*n2*n8+4*n2*n7+n2*n6+3*n2*n5+n2^2-39*ep+7*ep*n8+14*ep*n7+2*ep*n6+9*ep*n5+7*ep*n2+10*ep^2,-1+n8)
+`Y'(0,n2+1,0,1,n5,n6,n7,n8,0)*rat(-2*n2,-3+n8+n7+n2+2*ep)
+`Y'(0,n2+1,1,1,n5,n6,n7,n8-1,0)*rat(n2,-1+n8)*rat(-8+2*n8+2*n7+n6+2*n5+n2+5*ep,-3+n8+n7+n2+2*ep)
+`Y'(1,n2,1,0,n5+1,n6,n7,n8-1,0)*rat(n5,-1+n8)*rat(-8+2*n8+2*n7+n6+2*n5+n2+5*ep,-3+n8+n7+n2+2*ep)
+`Y'(0,n2,1,1,n5,n6,n7,n8,0)*rat(-11+3*n8+2*n7+2*n6+2*n5+2*n2+8*ep,-3+n8+n7+n2+2*ep)
+`Y'(1,n2,0,1,n5,n6,n7,n8,0)*rat(-10+2*n8+3*n7+n6+2*n5+2*n2+7*ep,-3+n8+n7+n2+2*ep)
;
*--#] nostar6_02_x8 : 
#endprocedure
*--#] reduce58 : 
*--#[ reduce27 :
#procedure reduce27
*--#[ nostar6_02_x27 : Z(1,n2,1,1,1,n6,n7,1,0)
* Z(1,n2,1,1,1,n6+ep,n7,1,0) gives 2+8=10 terms (0+2)+(0+8) ~ 2.044431 [+,+,+,+,+,x,+,+,-]
* checked poles: (-2+n7+n2+2*ep), (n7-n2)
* decreases: (n2+n7)
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,                 n7
if (match(`Z'(1,n2?!{n7?},1,1,1,n6?,n7?!{n2?},1,0)))
`id' `Z'(1,n2?pos_,1,1,1,n6?,n7?pos_,1,0) =
+`Z'(1,n2-1,1,1,1,n6,n7,1,0)*rat(-2+n6+n2+3*ep,n7-n2)*rat(-4+2*n7+n6+n2+5*ep,-2+n7+n2+2*ep)
+`Z'(1,n2,1,1,1,n6,n7-1,1,0)*rat(2-n7-n6-3*ep,n7-n2)*rat(-4+n7+n6+2*n2+5*ep,-2+n7+n2+2*ep)
+`Y'(0,n2+1,0,1,1,n6,n7,1,0)*rat(-2*n2,-2+n7+n2+2*ep)
+`Y'(1,n2+1,0,1,1,n6,n7-1,1,0)*rat(2*n2,n7-n2)*rat(-4+2*n7+n6+n2+5*ep,-2+n7+n2+2*ep)
+`Y'(1,n2+1,1,1,1,n6,n7-1,0,0)*rat(-2*n2,n7-n2)*rat(-4+n7+n6+2*n2+5*ep,-2+n7+n2+2*ep)
+`Y'(0,n2,1,1,1,n6,n7,1,0)*rat(2-n6-n2-3*ep,n7-n2)*rat(-4+n7+n6+2*n2+5*ep,-2+n7+n2+2*ep)
+`Y'(1,n2,0,1,1,n6,n7,1,0)*rat(-6+3*n7+n6+2*n2+7*ep,-2+n7+n2+2*ep)
+`Y'(1,n2,1,0,1,n6,n7,1,0)*rat(-4+2*n7+n6+n2+5*ep,n7-n2)
+`Y'(1,n2,1,1,0,n6,n7,1,0)*rat(4-n7-n6-2*n2-5*ep,n7-n2)
+`Y'(1,n2,1,1,1,n6,n7,0,0)*rat(2-n6-n2-3*ep,-2+n7+n2+2*ep)
;
*--#] nostar6_02_x27 : 
*--#[ nostar6_02_x22 : Z(1,n2,1,1,1,n6,n2,1,0)
* Z(1,n2,1,1,1,n6+ep,n2,1,0) gives 2+12=14 terms (0+2)+(0+12) ~ 2.051963 [+,+,+,+,+,x,+,+,-]
* checked poles: (-1+n2), (-1+n2+ep), (-2+n6+n2+3*ep), [(-4+n6+3*n2+5*ep), (n2)]
* decreases: (n7)
* +:     n2
* 0: n1, n2, n3, n4, n5, n6,     n8, n9
* -:                         n7
`id' `Z'(1,n2?{>=2},1,1,1,n6?,n2?{>=2},1,0) =
+`Z'(1,n2+1,1,1,1,n6,n2-2,1,0)*rat(1,4)*rat(n2,-1+n2)*rat(n2+ep,-1+n2+ep)*rat(-1+n6+n2+3*ep,-1+n2+ep)*rat(-3+n6+n2+3*ep,-2+n6+n2+3*ep)*rat(-4+n6+3*n2+5*ep,1)
+`Z'(1,n2,1,1,1,n6,n2-1,1,0)*rat(-1,4)*rat(1,-1+n2+ep)^2*rat(1,-2+n6+n2+3*ep)*rat(24-8*n6-56*n2+16*n2*n6+36*n2^2+7*n2^2*n6-8*n2^2*n6^2+n2^2*n6^3+5*n2^3-20*n2^3*n6+5*n2^3*n6^2-12*n2^4+7*n2^4*n6+3*n2^5-80*ep+16*ep*n6+124*ep*n2-ep*n2*n6-8*ep*n2*n6^2+ep*n2*
  n6^3+8*ep*n2^2-72*ep*n2^2*n6+16*ep*n2^2*n6^2-76*ep*n2^3+41*ep*n2^3*n6+26*ep*n2^4+88*ep^2-8*ep^2*n6-29*ep^2*n2-52*ep^2*n2*n6+11*ep^2*n2*n6^2-148*ep^2*n2^2+73*ep^2*n2^2*n6+80*ep^2*n2^3-32*ep^3-84*ep^3*n2+39*ep^3*n2*n6+102*ep^3*n2^2+45*ep^4*n2,-1+n2)
+`Y'(0,n2,1,0,1,n6,n2+1,1,0)*rat(-n2,-1+n2+ep)
+`Y'(0,n2+1,1,0,1,n6,n2,1,0)*rat(-1,2)*rat(n2,-1+n2+ep)*rat(-1+n6+n2+3*ep,-1+n2+ep)*rat(-4+n6+3*n2+5*ep,-2+n6+n2+3*ep)
+`Y'(0,n2+1,1,1,0,n6,n2,1,0)*rat(1,2)*rat(n2,-1+n2+ep)*rat(-3+n6+n2+3*ep,-1+n2+ep)*rat(-4+n6+3*n2+5*ep,-2+n6+n2+3*ep)
+`Y'(0,n2+1,1,1,1,n6,n2-1,1,0)*rat(-1,4)*rat(n2,-1+n2)*rat(-2+n2+ep,-1+n2+ep)*rat(-1+n6+n2+3*ep,-1+n2+ep)*rat(-3+n6+n2+3*ep,-2+n6+n2+3*ep)*rat(-4+n6+3*n2+5*ep,1)
+`Y'(1,n2-1,1,0,1,n6,n2+1,1,0)*rat(n2,-1+n2+ep)
+`Y'(1,n2+1,0,1,1,n6,n2-1,1,0)*rat(1,2)*rat(n2,-1+n2)*rat(-4+n6+3*n2+5*ep,-1+n2+ep)*rat(2-3*n2+n2*n6+n2^2-5*ep+ep*n6+4*ep*n2+3*ep^2,-2+n6+n2+3*ep)
+`Y'(1,n2+1,1,0,1,n6,n2-1,1,0)*rat(1,2)*rat(n2,-1+n2)*rat(n2+ep,-1+n2+ep)*rat(-1+n6+n2+3*ep,-1+n2+ep)*rat(-4+n6+3*n2+5*ep,-2+n6+n2+3*ep)
+`Y'(1,n2+1,1,1,0,n6,n2-1,1,0)*rat(-1,2)*rat(n2,-1+n2)*rat(n2+ep,-1+n2+ep)*rat(-3+n6+n2+3*ep,-1+n2+ep)*rat(-4+n6+3*n2+5*ep,-2+n6+n2+3*ep)
+`Y'(1,n2+1,1,1,1,n6,n2-1,0,0)*rat(-1,2)*rat(n2,-1+n2)*rat(-4+n6+3*n2+5*ep,-1+n2+ep)*rat(-4+n2+n2*n6+n2^2+ep+ep*n6+4*ep*n2+3*ep^2,-2+n6+n2+3*ep)
+`Y'(1,n2,0,1,1,n6,n2,1,0)
+`Y'(1,n2,1,0,1,n6,n2,1,0)*rat(1,2)*rat(1,-1+n2+ep)*rat(1,-2+n6+n2+3*ep)*rat(-20+6*n6+46*n2-10*n2*n6-30*n2^2-n2^2*n6+n2^2*n6^2+n2^3+4*n2^3*n6+3*n2^4+66*ep-12*ep*n6-98*ep*n2+5*ep*n2*n6+ep*n2*n6^2+14*ep*n2^2+12*ep*n2^2*n6+17*ep*n2^3-72*ep^2+6*ep^2*n6+39*
  ep^2*n2+8*ep^2*n2*n6+29*ep^2*n2^2+26*ep^3+15*ep^3*n2,-1+n2+ep)
+`Y'(1,n2,1,1,0,n6,n2,1,0)*rat(-1,2)*rat(1,-1+n2+ep)*rat(1,-2+n6+n2+3*ep)*rat(12-4*n6-24*n2+6*n2*n6+26*n2^2-9*n2^2*n6+n2^2*n6^2-15*n2^3+4*n2^3*n6+3*n2^4-40*ep+8*ep*n6+66*ep*n2-13*ep*n2*n6+ep*n2*n6^2-56*ep*n2^2+12*ep*n2^2*n6+17*ep*n2^3+44*ep^2-4*ep^2*n6-57
  *ep^2*n2+8*ep^2*n2*n6+29*ep^2*n2^2-16*ep^3+15*ep^3*n2,-1+n2+ep)
;
*--#] nostar6_02_x22 : 
#endprocedure
*--#] reduce27 : 
*--#[ reduce6 :
#procedure reduce6
*--#[ nostar6_02_p6 : Z(1,1,1,1,1,n6,1,1,0)
* Z(1,1,1,1,1,n6+ep,1,1,0) gives 1+5=6 terms (0+1)+(0+5) ~ 1.064139 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n6+ep)
* decreases: (n6)
* 0: n1, n2, n3, n4, n5,     n7, n8, n9
* -:                     n6
`id' `Z'(1,1,1,1,1,n6?{>=2},1,1,0) =
+`Z'(1,1,1,1,1,n6-1,1,1,0)*rat(-1+n6+5*ep,-1+n6+ep)
+`Y'(1,0,1,1,1,n6-1,1,2,0)*rat(1,-1+n6+ep)
+`Y'(1,1,0,1,1,n6-1,1,2,0)*rat(-1,-1+n6+ep)
+`Y'(1,1,1,0,2,n6-1,1,1,0)*rat(-1,-1+n6+ep)
+`Y'(1,1,1,1,2,n6-1,0,1,0)*rat(1,-1+n6+ep)
+`Y'(0,1,1,1,1,n6,1,1,0)
;
*--#] nostar6_02_p6 : 
*--#[ nostar6_02_m6 : Z(1,1,1,1,1,n6,1,1,0)
* Z(1,1,1,1,1,n6+ep,1,1,0) gives 1+5=6 terms (0+1)+(0+4+1) ~ 1.064139 [+,+,+,+,+,-,+,+,-]
* checked poles: (n6+5*ep)
* decreases: (-n6)
* +:                     n6
* 0: n1, n2, n3, n4, n5,     n7, n8, n9
`id' `Z'(1,1,1,1,1,n6?neg0_,1,1,0) =
+`Z'(1,1,1,1,1,n6+1,1,1,0)*rat(n6+ep,n6+5*ep)
+`Y'(1,0,2,1,1,n6+1,1,1,0)*rat(-1,n6+5*ep)
+`Y'(1,1,1,2,0,n6+1,1,1,0)*rat(1,n6+5*ep)
+`Y'(1,1,1,2,1,n6+1,0,1,0)*rat(-1,n6+5*ep)
+`Y'(1,1,2,1,1,n6+1,1,0,0)*rat(1,n6+5*ep)
+`Y'(0,1,1,1,1,n6+1,1,1,0)
;
*--#] nostar6_02_m6 : 
#endprocedure
*--#] reduce6 : 
  #call StdReduceProc(turbo,9,34,1,58,27,6)
#endprocedure
* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4
