#procedure SchemeZfastar2()
  #ifndef `TURBOTHRESHOLD'
    #define TURBOTHRESHOLD "4"
  #endif
*--#[ def :
* internals: p1,p6,p7
* externals: p0
* propagators: p1.p1,p2.p2,p3.p3,p4.p4,p5.p5,p6.p6,p7.p7,p1.p3,p0.p2
* replacements:
*   p3 = p1 + p6 + p7,
*   p5 = -p0 + p1,
*   p8 = p0,
*   p2 = p1 + p6,
*   p4 = p0 - p1 - p6 - p7,
*   p0.p0 = 1,
*   p8.p8 = 1,
* signs: +,*,+,+,+,+,+,-,-
*--#] def : 
*--#[ z2y :
#procedure z2y()
  #call Z2Y(9,7,2)
#endprocedure
*--#] z2y : 
*--#[ turbo :
#procedure turbo()
* depth = 3 (same as 4), threshold = 12
*--#[ n2 >= 2 :
#if `TURBOTHRESHOLD' >= 3
*--#[ fastar2_1302_turbo-345b42463b : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+3+3+3) ~ 2.246969 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n4), [(n8)]
* decreases: (-n8), (n1+n3+n4)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3, n4, n5, n6
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(-9+n9+n8+n7+n6+2*n5+n4+n3+n2+n1+5*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)
+`Z'(n1,n2,n3,n4-1,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n4)
;
*--#] fastar2_1302_turbo-345b42463b : 
*--#[ fastar2_1302_turbo-04aff0d6bd : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+3+3+3) ~ 2.246969 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n5), [(n8)]
* decreases: (-n8), (n1+n3+n5)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3, n4, n5,     n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8+1,n9)*rat(1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-1+n5)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-9+n9+n8+n7+n6+n5+2*n4+n3+n2+n1+5*ep,-1+n5)
+`Z'(n1,n2,n3,n4-1,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-1+n5)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n5)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n5)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n5)
;
*--#] fastar2_1302_turbo-04aff0d6bd : 
*--#[ fastar2_1302_turbo-c771e34453 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 6+0=6 terms (0+0+5+0+0+1) ~ 2.268871 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n6), (-1+n7), [(n8)]
* decreases: (-n8), (n1+n3+n6), (n1+n3+n7)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3,         n6, n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?{>=2},n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(-3+n7+n6+ep,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(-3+n7+n6+ep,-1+n6)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7-1,n8+1,n9)*rat(1,2)*rat(-3+n7+n6+n2+2*ep,-1+n6)*rat(-4+n9+n7+n6+n2+3*ep,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7-1,n8+1,n9+2)*rat(1,8)*rat(n9,-1+n6)*rat(1+n9,-1+n7)
;
*--#] fastar2_1302_turbo-c771e34453 : 
*--#[ fastar2_1302_turbo-005f07a760 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+0+3+3+5) ~ 2.318280 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n1), [(n8)]
* decreases: (n1), (-n8)
* +:     n2,                     n8, n9
* 0:     n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3, n4,         n7
`id' `Z'(n1?{>=2},n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n2+ep,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(1,2)*rat(n2+ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-5+n9+n7+n6+n5+n4-n3+n1+2*ep,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n1)
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+2,n9)*rat(-1,2)*rat(1+n8,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n1)
;
*--#] fastar2_1302_turbo-005f07a760 : 
*--#[ fastar2_1302_turbo-b73e924f0e : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+0+3+3+5) ~ 2.318280 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n3), [(n8)]
* decreases: (n3), (-n8)
* +:     n2,                     n8, n9
* 0: n1, n2,     n4, n5, n6, n7,     n9
* -: n1, n2, n3,     n5, n6
`id' `Z'(n1?pos_,n2?{>=2},n3?{>=2},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n2+ep,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(1,2)*rat(n2+ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-5+n9+n7+n6+n5+n4+n3-n1+2*ep,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+2,n9)*rat(-1,2)*rat(1+n8,-1+n3)
+`Z'(n1,n2-1,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n3)
;
*--#] fastar2_1302_turbo-b73e924f0e : 
*--#[ fastar2_1302_turbo-844672029f : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 17+0=17 terms (0+0+2+4+8+1+2) ~ 2.349528 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n7), [(n8), (n9)]
* decreases: (-n9), (n2+n7-n8), (n1+n2+n3+n6+n7)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5, n6, n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg_,n9?neg_) =
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(4-n8-n7-n4-n3-2*ep,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(1,4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(1,4)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(-4+2*n8+n7+n5+2*n3+2*ep,-1+n7)
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(-1,4)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(-2+2*n8+n7-n5+2*n3+2*ep,-1+n7)
+`Z'(n1,n2-1,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(-1,4)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(-11+n9+3*n8+2*n7+2*n6+n5+2*n3+2*n2+2*n1+8*ep,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7-1,n8+2,n9+1)*rat(-1,4)*rat(1+n8,-1+n7)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7-1,n8+2,n9+1)*rat(1,4)*rat(1+n8,-1+n7)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7-1,n8+2,n9+1)*rat(1,4)*rat(1+n8,-1+n7)
;
*--#] fastar2_1302_turbo-844672029f : 
*--#[ fastar2_1302_turbo-4be6e1923c : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 17+0=17 terms (0+0+2+4+8+1+2) ~ 2.349528 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n6), [(n8), (n9)]
* decreases: (-n9), (n2+n6-n8), (n1+n2+n3+n6+n7)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5, n6, n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg_,n9?neg_) =
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(4-n8-n6-n5-n1-2*ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(1,4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(1,4)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(-4+2*n8+n6+n4+2*n1+2*ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(-1,4)
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)
+`Z'(n1,n2-1,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,4)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,4)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(-2+2*n8+n6-n4+2*n1+2*ep,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(-11+n9+3*n8+2*n7+2*n6+n4+2*n3+2*n2+2*n1+8*ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(1,4)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+2,n9+1)*rat(-1,4)*rat(1+n8,-1+n6)
+`Z'(n1-1,n2,n3-1,n4,n5,n6-1,n7,n8+2,n9+1)*rat(1,4)*rat(1+n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5-1,n6-1,n7,n8+2,n9+1)*rat(1,4)*rat(1+n8,-1+n6)
;
*--#] fastar2_1302_turbo-4be6e1923c : 
*--#[ fastar2_1302_turbo-3ba8670d05 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 6+0=6 terms (0+1+4+0+0+0+1) ~ 2.584517 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n5), (-1+n6), [(n9)]
* decreases: (-n9)
* +: n1,                         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3,     n5, n6
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?{>=2},n6?{>=2},n7?pos_,n8?neg0_,n9?neg_) =
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1+1,n2,n3,n4,n5-1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n5)*rat(-1+n8+n1+ep,-1+n6)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(-1,2)*rat(-3+n6+n5+ep,-1+n5)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-3+n6+n5+ep,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5-1,n6-1,n7,n8+2,n9+1)*rat(1,8)*rat(n8,-1+n5)*rat(1+n8,-1+n6)
;
*--#] fastar2_1302_turbo-3ba8670d05 : 
*--#[ fastar2_1302_turbo-c72892f111 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 6+0=6 terms (0+1+4+0+0+0+1) ~ 2.584517 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n4), (-1+n7), [(n9)]
* decreases: (-n9)
* +:         n3,                 n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2,     n4,         n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?{>=2},n8?neg0_,n9?neg_) =
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n3,-1+n4)*rat(-1+n8+n3+ep,-1+n7)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-3+n7+n4+ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-3+n7+n4+ep,-1+n7)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7-1,n8+2,n9+1)*rat(1,8)*rat(n8,-1+n4)*rat(1+n8,-1+n7)
;
*--#] fastar2_1302_turbo-c72892f111 : 
*--#[ fastar2_1302_turbo-831ea847f9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+6+1+2) ~ 2.607858 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n6), [(n8)]
* decreases: (-n8), (n1+n3+n6+n7)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(1,2)*rat(n2+ep,-1+n6)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7-1,n8+1,n9)*rat(-1,2)*rat(n2+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(-5+n9+2*n7+n6+n2+3*ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n6)
;
*--#] fastar2_1302_turbo-831ea847f9 : 
*--#[ fastar2_1302_turbo-123db3b30d : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+6+1+2) ~ 2.607858 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n7), [(n8)]
* decreases: (-n8), (n1+n3+n6+n7)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3,     n5, n6, n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(1,2)*rat(n2+ep,-1+n7)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7-1,n8+1,n9)*rat(-1,2)*rat(n2+ep,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(-5+n9+n7+2*n6+n2+3*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n7)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n7)
;
*--#] fastar2_1302_turbo-123db3b30d : 
*--#[ fastar2_1302_turbo-fc9e506cb8 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+1+5+0+3) ~ 2.895401 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n5), [(n9)]
* decreases: (-n9)
* +: n1,                         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3,     n5, n6, n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1+1,n2-1,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n5)
+`Z'(n1+1,n2,n3,n4,n5-1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n1,-1+n5)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(-1,2)*rat(-5+n8+2*n6+n5+n1+2*ep,-1+n5)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n5)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n5)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(n8,-1+n5)
;
*--#] fastar2_1302_turbo-fc9e506cb8 : 
*--#[ fastar2_1302_turbo-2142c4d47e : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+1+5+0+3) ~ 2.895401 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n9)
* +:         n3,                 n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2,     n4,     n6, n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n4)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n3,-1+n4)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-5+n8+2*n7+n4+n3+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n8,-1+n4)
;
*--#] fastar2_1302_turbo-2142c4d47e : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ fastar2_1302_turbo-b4720c0f1a : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+2+4+1+2) ~ 3.369162 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n6), [(n9)]
* decreases: (-n9)
* +: n1,                         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3, n4, n5, n6
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg0_,n9?neg_) =
+`Z'(n1+1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n6)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1+1,n2,n3,n4,n5-1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n1,-1+n6)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-5+n8+n6+2*n5+n1+2*ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n8,-1+n6)
;
*--#] fastar2_1302_turbo-b4720c0f1a : 
*--#[ fastar2_1302_turbo-9d75b9e4bd : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+2+4+1+2) ~ 3.369162 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n7), [(n9)]
* decreases: (-n9)
* +:         n3,                 n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2,     n4, n5,     n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg0_,n9?neg_) =
+`Z'(n1,n2,n3+1,n4,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n3,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n3,-1+n7)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-5+n8+n7+2*n4+n3+2*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n7)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(n8,-1+n7)
;
*--#] fastar2_1302_turbo-9d75b9e4bd : 
*--#[ fastar2_1302_turbo-d6b4ce48fd : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 7+0=7 terms (0+3+1+2+0+1) ~ 3.487540 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n1), (-1+n5)
* decreases: (n1+n5)
* +:                     n6,     n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5
`id' `Z'(n1?{>=2},n2?{>=2},n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)*rat(-n6,-1+n1)*rat(-1+n6+ep,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)*rat(-3+n8+n5+n1+ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9)*rat(-3+n8+n5+n1+ep,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(n8,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5-1,n6,n7,n8+2,n9)*rat(1,4)*rat(n8,-1+n1)*rat(1+n8,-1+n5)
;
*--#] fastar2_1302_turbo-d6b4ce48fd : 
#endif
#if `TURBOTHRESHOLD' >= 5
*--#[ fastar2_1302_turbo-af365da268 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 4+0=4 terms (0+4) ~ 4.000000 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                         n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2, n3, n4
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-5+n8+n7+n4+2*n3+2*ep,-1+n4)
;
*--#] fastar2_1302_turbo-af365da268 : 
*--#[ fastar2_1302_turbo-e617710486 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 4+0=4 terms (0+4) ~ 4.000000 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n5)
* decreases: (n1+n5)
* +:                         n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3,     n5
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7+1,n8,n9)*rat(n7,-1+n5)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7+1,n8,n9)*rat(-n7,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9)*rat(-9+n9+n8+n7+2*n6+n5+2*n2+2*n1+6*ep,-1+n5)
;
*--#] fastar2_1302_turbo-e617710486 : 
#endif
*--#] n2 >= 2 : 
*--#[ n2 <= 0 :
#if `TURBOTHRESHOLD' >= 2
*--#[ fastar2_1302_turbo-199619afb7 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 6+0=6 terms (0+0+2+0+3+1) ~ 1.786025 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n6), (-1+n7), [(n8)]
* decreases: (-n8), (n1+n3+n6), (n1+n3+n7), (n1-n2+n3-n9)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3,         n6, n7
`id' `Z'(n1?pos_,n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?{>=2},n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(-3+n7+n6+ep,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(-3+n7+n6+ep,-1+n6)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7-1,n8+1,n9)*rat(1,2)*rat(-3+n7+n6+n2+2*ep,-1+n6)*rat(-4+n9+n7+n6+n2+3*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7-1,n8+1,n9+2)*rat(1,8)*rat(n9,-1+n6)*rat(1+n9,-1+n7)
;
*--#] fastar2_1302_turbo-199619afb7 : 
#endif
#if `TURBOTHRESHOLD' >= 3
*--#[ fastar2_1302_turbo-a5c4f828d1 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 7+0=7 terms (0+0+3+1+2+0+1) ~ 2.030855 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n3), (-1+n7)
* decreases: (-n2), (n3+n7)
* +:     n2,     n4,             n8
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3,     n5,     n7
`id' `Z'(n1?pos_,n2?neg0_,n3?{>=2},n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg0_,n9?neg0_) =
+`Z'(n1,n2+1,n3-1,n4+1,n5,n6,n7-1,n8,n9)*rat(-n4,-1+n3)*rat(-1+n8+n4+ep,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7,n8,n9)*rat(-3+n7+n3+ep,-1+n3)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8,n9)*rat(-3+n7+n3+ep,-1+n7)
+`Z'(n1,n2+1,n3-1,n4+1,n5,n6,n7-1,n8+1,n9)*rat(1,2)*rat(n4,-1+n3)*rat(n8,-1+n7)
+`Z'(n1-1,n2+1,n3-1,n4+1,n5,n6,n7-1,n8+1,n9)*rat(1,2)*rat(n4,-1+n3)*rat(n8,-1+n7)
+`Z'(n1,n2+1,n3-1,n4+1,n5-1,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(n4,-1+n3)*rat(n8,-1+n7)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7-1,n8+2,n9)*rat(-1,4)*rat(n8,-1+n3)*rat(1+n8,-1+n7)
;
*--#] fastar2_1302_turbo-a5c4f828d1 : 
*--#[ fastar2_1302_turbo-22ee9df5ca : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 7+0=7 terms (0+0+3+1+2+0+1) ~ 2.030855 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n1), (-1+n6)
* decreases: (-n2), (n1+n6)
* +:     n2,         n5,         n8
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4,     n6
`id' `Z'(n1?{>=2},n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-n5,-1+n1)*rat(-1+n8+n5+ep,-1+n6)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n1)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n6)
+`Z'(n1-1,n2+1,n3,n4,n5+1,n6-1,n7,n8+1,n9)*rat(1,2)*rat(n5,-1+n1)*rat(n8,-1+n6)
+`Z'(n1-1,n2+1,n3-1,n4,n5+1,n6-1,n7,n8+1,n9)*rat(1,2)*rat(n5,-1+n1)*rat(n8,-1+n6)
+`Z'(n1-1,n2+1,n3,n4-1,n5+1,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(n5,-1+n1)*rat(n8,-1+n6)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6-1,n7,n8+2,n9)*rat(-1,4)*rat(n8,-1+n1)*rat(1+n8,-1+n6)
;
*--#] fastar2_1302_turbo-22ee9df5ca : 
*--#[ fastar2_1302_turbo-d413534641 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+3+1+5) ~ 2.157813 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n7), [(n8)]
* decreases: (-n8), (n1+n3+n7)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3,     n5, n6, n7
`id' `Z'(n1?pos_,n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(-4+n9+n7+n6+n2+3*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n7)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(-1+n6-n2-ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7-1,n8+1,n9)*rat(1,2)*rat(-4+n9+n7+n6+n2+3*ep,-1+n7)
;
*--#] fastar2_1302_turbo-d413534641 : 
*--#[ fastar2_1302_turbo-f4a79584a9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+3+1+5) ~ 2.157813 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n6), [(n8)]
* decreases: (-n8), (n1+n3+n6)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(-4+n9+n7+n6+n2+3*ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n6)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n6)
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(-1+n7-n2-ep,-1+n6)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7-1,n8+1,n9)*rat(1,2)*rat(-4+n9+n7+n6+n2+3*ep,-1+n6)
;
*--#] fastar2_1302_turbo-f4a79584a9 : 
*--#[ fastar2_1302_turbo-c245863444 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+0+2+3+6) ~ 2.185564 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n1), [(n8)]
* decreases: (n1), (-n8)
* +:     n2,                     n8, n9
* 0:     n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3, n4,         n7
`id' `Z'(n1?{>=2},n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-5+n9+n7+n6+n5+n4-n3+n1+2*ep,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+2,n9)*rat(-1,2)*rat(1+n8,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n1)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n2+ep,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(1,2)*rat(n2+ep,-1+n1)
;
*--#] fastar2_1302_turbo-c245863444 : 
*--#[ fastar2_1302_turbo-912f883e9c : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+0+2+3+6) ~ 2.185564 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n3), [(n8)]
* decreases: (n3), (-n8)
* +:     n2,                     n8, n9
* 0: n1, n2,     n4, n5, n6, n7,     n9
* -: n1, n2, n3,     n5, n6
`id' `Z'(n1?pos_,n2?neg0_,n3?{>=2},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1,n2-1,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-5+n9+n7+n6+n5+n4+n3-n1+2*ep,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+2,n9)*rat(-1,2)*rat(1+n8,-1+n3)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n2+ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(1,2)*rat(n2+ep,-1+n3)
;
*--#] fastar2_1302_turbo-912f883e9c : 
*--#[ fastar2_1302_turbo-3d1b63ef43 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 8+0=8 terms (0+0+4+1+2+0+1) ~ 2.234417 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n3), (-1+n6)
* decreases: (-n2), (n1+n6)
* +:     n2,     n4,         n7, n8
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3,     n5, n6
`id' `Z'(n1?pos_,n2?neg0_,n3?{>=2},n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7+1,n8,n9)*rat(-n7,-1+n3)*rat(-1+n7+ep,-1+n6)
+`Z'(n1,n2+1,n3-1,n4+1,n5,n6-1,n7,n8,n9)*rat(n4,-1+n3)*rat(-1+n8+n4+ep,-1+n6)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7,n8,n9)*rat(5-n9-2*n7-n6-n3-2*n2-5*ep,-1+n6)
+`Z'(n1,n2+1,n3-1,n4+1,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(n4,-1+n3)*rat(n8,-1+n6)
+`Z'(n1-1,n2+1,n3-1,n4+1,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(n4,-1+n3)*rat(n8,-1+n6)
+`Z'(n1,n2+1,n3-1,n4+1,n5-1,n6-1,n7,n8+1,n9)*rat(1,2)*rat(n4,-1+n3)*rat(n8,-1+n6)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6-1,n7,n8+2,n9)*rat(1,4)*rat(n8,-1+n3)*rat(1+n8,-1+n6)
;
*--#] fastar2_1302_turbo-3d1b63ef43 : 
*--#[ fastar2_1302_turbo-9b3582a320 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 8+0=8 terms (0+0+4+1+2+0+1) ~ 2.234417 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n1), (-1+n7), [(n4)]
* decreases: (-n2), (n3+n7)
* +:     n2,         n5, n6,     n8
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4,         n7
`id' `Z'(n1?{>=2},n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg0_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3,n4,n5,n6+1,n7-1,n8,n9)*rat(-n6,-1+n1)*rat(-1+n6+ep,-1+n7)
+`Z'(n1-1,n2+1,n3,n4,n5+1,n6,n7-1,n8,n9)*rat(n5,-1+n1)*rat(-1+n8+n5+ep,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8,n9)*rat(5-n9-n7-2*n6-2*n2-n1-5*ep,-1+n7)
+`Z'(n1-1,n2+1,n3,n4,n5+1,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(n5,-1+n1)*rat(n8,-1+n7)
+`Z'(n1-1,n2+1,n3-1,n4,n5+1,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(n5,-1+n1)*rat(n8,-1+n7)
+`Z'(n1-1,n2+1,n3,n4-1,n5+1,n6,n7-1,n8+1,n9)*rat(1,2)*rat(n5,-1+n1)*rat(n8,-1+n7)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7-1,n8+2,n9)*rat(1,4)*rat(n8,-1+n1)*rat(1+n8,-1+n7)
;
*--#] fastar2_1302_turbo-9b3582a320 : 
*--#[ fastar2_1302_turbo-8c5a3dfec6 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 4+0=4 terms (0+1+3) ~ 2.311676 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n6)
* decreases: (-n2), (n1+n6)
* +:     n2,         n5
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -: n1,                 n6
`id' `Z'(n1?pos_,n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,n2+1,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-n5,-1+n6)
+`Z'(n1-1,n2+1,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(n5,-1+n6)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7,n8,n9)*rat(-5+n8+n6+n5+2*n1+2*ep,-1+n6)
;
*--#] fastar2_1302_turbo-8c5a3dfec6 : 
*--#[ fastar2_1302_turbo-a591587031 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 4+0=4 terms (0+1+3) ~ 2.311676 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n7)
* decreases: (-n2), (n3+n7)
* +:     n2,     n4
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -:         n3,             n7
`id' `Z'(n1?pos_,n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg0_,n9?neg0_) =
+`Z'(n1,n2+1,n3,n4+1,n5,n6,n7-1,n8,n9)*rat(-n4,-1+n7)
+`Z'(n1,n2+1,n3-1,n4+1,n5,n6,n7-1,n8,n9)*rat(n4,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8,n9)*rat(-5+n8+n7+n4+2*n3+2*ep,-1+n7)
;
*--#] fastar2_1302_turbo-a591587031 : 
*--#[ fastar2_1302_turbo-44e8faaa6e : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+4+3+2) ~ 2.392541 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n4), [(n8)]
* decreases: (-n8), (n1+n3+n4)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3, n4, n5, n6
`id' `Z'(n1?pos_,n2?neg0_,n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n4)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(-9+n9+n8+n7+n6+2*n5+n4+n3+n2+n1+5*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)
+`Z'(n1,n2,n3,n4-1,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n4)
;
*--#] fastar2_1302_turbo-44e8faaa6e : 
*--#[ fastar2_1302_turbo-f5d415cba9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+4+3+2) ~ 2.392541 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n5), [(n8)]
* decreases: (-n8), (n1+n3+n5)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3, n4, n5,     n7
`id' `Z'(n1?pos_,n2?neg0_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8+1,n9)*rat(1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-1+n5)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-9+n9+n8+n7+n6+n5+2*n4+n3+n2+n1+5*ep,-1+n5)
+`Z'(n1,n2,n3,n4-1,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-1+n5)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n5)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n5)
;
*--#] fastar2_1302_turbo-f5d415cba9 : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ fastar2_1302_turbo-e81b8a426d : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 7+0=7 terms (0+3+1+2+0+1) ~ 3.487540 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (n3+n4)
* +:                         n7, n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5
`id' `Z'(n1?pos_,n2?neg0_,n3?{>=2},n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7+1,n8,n9)*rat(-n7,-1+n3)*rat(-1+n7+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)*rat(-3+n8+n4+n3+ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-3+n8+n4+n3+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(1,2)*rat(n8,-1+n3)
+`Z'(n1-1,n2,n3-1,n4-1,n5,n6,n7,n8+2,n9)*rat(1,4)*rat(n8,-1+n3)*rat(1+n8,-1+n4)
;
*--#] fastar2_1302_turbo-e81b8a426d : 
*--#[ fastar2_1302_turbo-174a4cba40 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 7+0=7 terms (0+3+1+2+0+1) ~ 3.487540 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n1), (-1+n5)
* decreases: (n1+n5)
* +:                     n6,     n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5
`id' `Z'(n1?{>=2},n2?neg0_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)*rat(-n6,-1+n1)*rat(-1+n6+ep,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)*rat(-3+n8+n5+n1+ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9)*rat(-3+n8+n5+n1+ep,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(n8,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5-1,n6,n7,n8+2,n9)*rat(1,4)*rat(n8,-1+n1)*rat(1+n8,-1+n5)
;
*--#] fastar2_1302_turbo-174a4cba40 : 
#endif
#if `TURBOTHRESHOLD' >= 5
*--#[ fastar2_1302_turbo-1a4cae91d8 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 18+0=18 terms (0+3+3+6+2+4) ~ 4.179887 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n3)
* decreases: (-n2+n3+n4), (n3+n4+n7)
* +:     n2,                     n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6, n7,     n9
`id' `Z'(n1?pos_,n2?neg0_,n3?{>=2},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7,n8,n9-1)*rat(-14+n9+2*n8+3*n7+2*n6+n4+3*n3+2*n2+n1+8*ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)*rat(-1,2)*rat(-13+n9+2*n8+4*n7+3*n6+n5+n3+2*n2+n1+8*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7,n8,n9)*rat(-1,2)*rat(-3+n9+n6-n5+n3+2*n2+n1+4*ep,-1+n3)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8,n9)*rat(1,2)
+`Z'(n1,n2,n3-2,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(1,4)*rat(n8,-1+n3)
+`Z'(n1,n2+1,n3-2,n4,n5,n6,n7,n8,n9)*rat(5-n8-n7-n4-n3-2*ep,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5-1,n6,n7,n8,n9)*rat(-1,2)*rat(-4+n8+n6+n5+n1+2*ep,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7,n8,n9)*rat(1,2)*rat(-4+n8+n6+n5+n1+2*ep,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,n9)*rat(-5+n8+n7+n4+n3+2*ep,-1+n3)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(1,4)*rat(n8,-1+n3)
+`Z'(n1-1,n2+1,n3-2,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(n8,-1+n3)
+`Z'(n1-1,n2+1,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(1,4)*rat(n8,-1+n3)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,4)*rat(n8,-1+n3)
+`Z'(n1,n2+1,n3-2,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n3)
;
*--#] fastar2_1302_turbo-1a4cae91d8 : 
*--#[ fastar2_1302_turbo-200e85317c : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 18+0=18 terms (0+3+3+6+3+3) ~ 4.187524 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n1)
* decreases: (n1-n2+n5), (n1+n5+n6)
* +:     n2,                     n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6, n7,     n9
`id' `Z'(n1?{>=2},n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8,n9-1)*rat(-6+n8+n6+n5+2*n1+2*ep,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)*rat(-1,2)*rat(n6-n5,-1+n1)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9)*rat(1,2)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8,n9)*rat(1,2)*rat(n6-n5,-1+n1)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7,n8,n9)*rat(1,2)
+`Z'(n1-2,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,4)*rat(n8,-1+n1)
+`Z'(n1-2,n2+1,n3,n4,n5,n6,n7,n8,n9)*rat(5-n8-n6-n5-n1-2*ep,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,4)*rat(n8,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(1,4)*rat(n8,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5-1,n6,n7,n8,n9)*rat(1,2)*rat(-5+n8+n6+n5+n1+2*ep,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6-1,n7,n8,n9)*rat(1,2)*rat(-5+n8+n6+n5+n1+2*ep,-1+n1)
+`Z'(n1-2,n2+1,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,4)*rat(n8,-1+n1)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,4)*rat(n8,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(1,4)*rat(n8,-1+n1)
+`Z'(n1-2,n2+1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(n8,-1+n1)
+`Z'(n1-2,n2+1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,4)*rat(n8,-1+n1)
+`Z'(n1-2,n2+1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,4)*rat(n8,-1+n1)
;
*--#] fastar2_1302_turbo-200e85317c : 
#endif
*--#] n2 <= 0 : 
*--#[ n2 = 1 :
#if `TURBOTHRESHOLD' >= 3
*--#[ fastar2_1302_turbo-38dc06d074 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 6+0=6 terms (0+0+5+0+0+1) ~ 2.268871 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n6), (-1+n7), [(n8)]
* decreases: (-n8), (n1+n3+n6), (n1+n3+n7)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3,         n6, n7
`id' `Z'(n1?pos_,1,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?{>=2},n8?neg_,n9?neg0_) =
+`Z'(n1-1,2,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(-3+n7+n6+ep,-1+n7)
+`Z'(n1,2,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(-3+n7+n6+ep,-1+n6)
+`Z'(n1,2,n3,n4,n5,n6-1,n7-1,n8+1,n9)*rat(1,2)*rat(-2+n7+n6+2*ep,-1+n6)*rat(-3+n9+n7+n6+3*ep,-1+n7)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3,n4,n5,n6-1,n7-1,n8+1,n9+2)*rat(1,8)*rat(n9,-1+n6)*rat(1+n9,-1+n7)
;
*--#] fastar2_1302_turbo-38dc06d074 : 
*--#[ fastar2_1302_turbo-dc7cfcbe65 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 14+0=14 terms (0+0+1+9+1+3) ~ 2.336930 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n1), (-1+n5), [(n8)]
* decreases: (n1), (-n8)
* +:     n2,             n6,     n8, n9
* 0:     n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3, n4, n5,     n7
`id' `Z'(n1?{>=2},1,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n1)
+`Z'(n1-2,1,n3,n4,n5-1,n6+1,n7,n8+1,n9)*rat(-1,2)*rat(n6,-1+n1)*rat(-1+n6+ep,-1+n5)
+`Z'(n1-1,1,n3-1,n4,n5-1,n6+1,n7,n8+1,n9)*rat(-1,2)*rat(n6,-1+n1)*rat(-1+n6+ep,-1+n5)
+`Z'(n1-1,1,n3,n4,n5-1,n6+1,n7-1,n8+1,n9)*rat(1,2)*rat(n6,-1+n1)*rat(-1+n6+ep,-1+n5)
+`Z'(n1-1,2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(1+ep,-1+n1)*rat(-4+n8+n6+n5+n1+2*ep,-1+n5)
+`Z'(n1-1,2,n3,n4,n5-1,n6,n7-1,n8+1,n9)*rat(1,2)*rat(1+ep,-1+n1)*rat(-4+n8+n6+n5+n1+2*ep,-1+n5)
+`Z'(n1-2,1,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-3+n8+n5+n1+ep,-1+n1)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-3+n5+n1+ep,-1+n1)
+`Z'(n1-1,1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(-3+n5+n1+ep,-1+n1)
+`Z'(n1-1,1,n3,n4,n5-1,n6,n7,n8+1,n9)*rat(1,2)*rat(1,-1+n5)*rat(30-4*n9-12*n8+n8*n9+n8^2-8*n7+2*n7*n8-10*n6+n6*n9+2*n6*n8+2*n6*n7+n6^2-12*n5+n5*n9+2*n5*n8+2*n5*n7+2*n5*n6+n5^2-15*n1+n1*n9+3*n1*n8+2*n1*n7+2*n1*n6+3*n1*n5+2*n1^2-30*ep+2*ep*n9+6*ep*n8+4*ep*
  n7+6*ep*n6+6*ep*n5+7*ep*n1+8*ep^2,-1+n1)
+`Z'(n1-1,1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n1)*rat(-4+n8+n6+n5+n1+2*ep,-1+n5)
+`Z'(n1-1,1,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n1)*rat(-4+n8+n6+n5+n1+2*ep,-1+n5)
+`Z'(n1-1,1,n3-1,n4,n5-1,n6,n7,n8+2,n9)*rat(-1,4)*rat(1+n8,-1+n1)*rat(-6+n8+2*n5+2*n1+2*ep,-1+n5)
+`Z'(n1-1,1,n3,n4-1,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n1)*rat(-4+n8+n6+n5+n1+2*ep,-1+n5)
;
*--#] fastar2_1302_turbo-dc7cfcbe65 : 
*--#[ fastar2_1302_turbo-43133af460 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 14+0=14 terms (0+0+1+9+1+3) ~ 2.336930 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n3), (-1+n4), [(n8)]
* decreases: (n3), (-n8)
* +:     n2,                 n7, n8, n9
* 0: n1, n2,     n4, n5, n6, n7,     n9
* -: n1,     n3, n4, n5, n6
`id' `Z'(n1?pos_,1,n3?{>=2},n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n3)
+`Z'(n1-1,1,n3-1,n4-1,n5,n6,n7+1,n8+1,n9)*rat(-1,2)*rat(n7,-1+n3)*rat(-1+n7+ep,-1+n4)
+`Z'(n1-1,2,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(1+ep,-1+n3)*rat(-4+n8+n7+n4+n3+2*ep,-1+n4)
+`Z'(n1,1,n3-2,n4-1,n5,n6,n7+1,n8+1,n9)*rat(-1,2)*rat(n7,-1+n3)*rat(-1+n7+ep,-1+n4)
+`Z'(n1,1,n3-1,n4-1,n5,n6-1,n7+1,n8+1,n9)*rat(1,2)*rat(n7,-1+n3)*rat(-1+n7+ep,-1+n4)
+`Z'(n1,2,n3-1,n4-1,n5,n6-1,n7,n8+1,n9)*rat(1,2)*rat(1+ep,-1+n3)*rat(-4+n8+n7+n4+n3+2*ep,-1+n4)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-3+n4+n3+ep,-1+n3)
+`Z'(n1,1,n3-2,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-3+n8+n4+n3+ep,-1+n3)
+`Z'(n1,1,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(1,-1+n4)*rat(30-4*n9-12*n8+n8*n9+n8^2-10*n7+n7*n9+2*n7*n8+n7^2-8*n6+2*n6*n8+2*n6*n7-12*n4+n4*n9+2*n4*n8+2*n4*n7+2*n4*n6+n4^2-15*n3+n3*n9+3*n3*n8+2*n3*n7+2*n3*n6+3*n3*n4+2*n3^2-30*ep+2*ep*n9+6*ep*n8+6*ep*
  n7+4*ep*n6+6*ep*n4+7*ep*n3+8*ep^2,-1+n3)
+`Z'(n1,1,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-3+n4+n3+ep,-1+n3)
+`Z'(n1,1,n3-1,n4-1,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n3)*rat(-4+n8+n7+n4+n3+2*ep,-1+n4)
+`Z'(n1-1,1,n3-1,n4-1,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n3)*rat(-4+n8+n7+n4+n3+2*ep,-1+n4)
+`Z'(n1-1,1,n3-1,n4-1,n5,n6,n7,n8+2,n9)*rat(-1,4)*rat(1+n8,-1+n3)*rat(-6+n8+2*n4+2*n3+2*ep,-1+n4)
+`Z'(n1,1,n3-1,n4-1,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n3)*rat(-4+n8+n7+n4+n3+2*ep,-1+n4)
;
*--#] fastar2_1302_turbo-43133af460 : 
*--#[ fastar2_1302_turbo-8157639bd4 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+4+3+2) ~ 2.392541 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n4), [(n8)]
* decreases: (-n8), (n1+n3+n4)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3, n4, n5, n6
`id' `Z'(n1?pos_,1,n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1,0,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n4)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-3+n8+n7+n6+n3+n1+3*ep,-1+n4)
+`Z'(n1-1,1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(-8+n9+n8+n7+n6+2*n5+n4+n3+n1+5*ep,-1+n4)
+`Z'(n1,1,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)
+`Z'(n1,1,n3,n4-1,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-3+n8+n7+n6+n3+n1+3*ep,-1+n4)
+`Z'(n1-1,1,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n4)
+`Z'(n1,1,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n4)
;
*--#] fastar2_1302_turbo-8157639bd4 : 
*--#[ fastar2_1302_turbo-e52c64846d : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+4+3+2) ~ 2.392541 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n5), [(n8)]
* decreases: (-n8), (n1+n3+n5)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3, n4, n5,     n7
`id' `Z'(n1?pos_,1,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1,0,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n5)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3,n4,n5-1,n6,n7,n8+1,n9)*rat(1,2)*rat(-3+n8+n7+n6+n3+n1+3*ep,-1+n5)
+`Z'(n1-1,1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)
+`Z'(n1,1,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-8+n9+n8+n7+n6+n5+2*n4+n3+n1+5*ep,-1+n5)
+`Z'(n1,1,n3,n4-1,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-3+n8+n7+n6+n3+n1+3*ep,-1+n5)
+`Z'(n1,1,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n5)
+`Z'(n1,1,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n5)
;
*--#] fastar2_1302_turbo-e52c64846d : 
*--#[ fastar2_1302_turbo-a5b5ddf8f3 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+0+4+3+4) ~ 2.456929 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n1), [(n8)]
* decreases: (n1), (-n8)
* +:     n2,                     n8, n9
* 0:     n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3, n4,         n7
`id' `Z'(n1?{>=2},1,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1-1,0,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n1)
+`Z'(n1-1,2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(1+ep,-1+n1)
+`Z'(n1-1,2,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(1,2)*rat(1+ep,-1+n1)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-5+n9+n7+n6+n5+n4-n3+n1+2*ep,-1+n1)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-7+n9+n8+n7+n6+n5+n4+n3+n1+5*ep,-1+n1)
+`Z'(n1-1,1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(-7+n9+n8+n7+n6+n5+n4+n3+n1+5*ep,-1+n1)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n1)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n9,-1+n1)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+2,n9)*rat(-1,2)*rat(1+n8,-1+n1)
+`Z'(n1-1,1,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n1)
+`Z'(n1-1,1,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n1)
;
*--#] fastar2_1302_turbo-a5b5ddf8f3 : 
*--#[ fastar2_1302_turbo-b2a4bcb800 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+0+4+3+4) ~ 2.456929 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n3), [(n8)]
* decreases: (n3), (-n8)
* +:     n2,                     n8, n9
* 0: n1, n2,     n4, n5, n6, n7,     n9
* -: n1, n2, n3,     n5, n6
`id' `Z'(n1?pos_,1,n3?{>=2},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1-1,2,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(1+ep,-1+n3)
+`Z'(n1,0,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n3)
+`Z'(n1,2,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(1,2)*rat(1+ep,-1+n3)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-5+n9+n7+n6+n5+n4+n3-n1+2*ep,-1+n3)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-7+n9+n8+n7+n6+n5+n4+n3+n1+5*ep,-1+n3)
+`Z'(n1,1,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-7+n9+n8+n7+n6+n5+n4+n3+n1+5*ep,-1+n3)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n3)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+2,n9)*rat(-1,2)*rat(1+n8,-1+n3)
+`Z'(n1,1,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n3)
+`Z'(n1,1,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n3)
;
*--#] fastar2_1302_turbo-b2a4bcb800 : 
*--#[ fastar2_1302_turbo-492a86c214 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+6+1+2) ~ 2.607858 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n7), [(n8)]
* decreases: (-n8), (n1+n3+n6+n7)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3,     n5, n6, n7
`id' `Z'(n1?pos_,1,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg_,n9?neg0_) =
+`Z'(n1-1,2,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(1,2)*rat(1+ep,-1+n7)
+`Z'(n1,2,n3,n4,n5,n6-1,n7-1,n8+1,n9)*rat(-1,2)*rat(1+ep,-1+n7)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)
+`Z'(n1,1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)*rat(-4+n9+n7+2*n6+3*ep,-1+n7)
+`Z'(n1,1,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n7)
+`Z'(n1-1,1,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n7)
+`Z'(n1,1,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n7)
;
*--#] fastar2_1302_turbo-492a86c214 : 
*--#[ fastar2_1302_turbo-d04b909f81 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+6+1+2) ~ 2.607858 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n6), [(n8)]
* decreases: (-n8), (n1+n3+n6+n7)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1,     n3, n4,     n6, n7
`id' `Z'(n1?pos_,1,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1,2,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(1,2)*rat(1+ep,-1+n6)
+`Z'(n1,2,n3,n4,n5,n6-1,n7-1,n8+1,n9)*rat(-1,2)*rat(1+ep,-1+n6)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(1,2)
+`Z'(n1,1,n3,n4,n5,n6-1,n7,n8+1,n9)*rat(-1,2)*rat(-4+n9+2*n7+n6+3*ep,-1+n6)
+`Z'(n1,1,n3,n4,n5,n6,n7-1,n8+1,n9)*rat(-1,2)
+`Z'(n1,1,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n6)
+`Z'(n1,1,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-1+n6)
+`Z'(n1,1,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-1+n6)
;
*--#] fastar2_1302_turbo-d04b909f81 : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ fastar2_1302_turbo-c18a66450b : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 8+0=8 terms (0+2+5+0+0+0+1) ~ 3.475059 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n5), (-1+n6), [(n9)]
* decreases: (-n9)
* +: n1,             n5,         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3,     n5, n6
`id' `Z'(n1?pos_,1,n3?pos_,n4?pos_,n5?{>=2},n6?{>=2},n7?pos_,n8?neg0_,n9?neg_) =
+`Z'(n1,1,n3,n4,n5+1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n5,-1+n6)
+`Z'(n1,1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1-1,1,n3,n4,n5+1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n5,-1+n6)
+`Z'(n1+1,1,n3,n4,n5-1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n5)*rat(-1+n8+n1+ep,-1+n6)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,1,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(-1,2)*rat(-3+n6+n5+ep,-1+n5)
+`Z'(n1,1,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(-2+n8+2*n1+ep,-1+n6)
+`Z'(n1,1,n3-1,n4,n5-1,n6-1,n7,n8+2,n9+1)*rat(1,8)*rat(n8,-1+n5)*rat(1+n8,-1+n6)
;
*--#] fastar2_1302_turbo-c18a66450b : 
*--#[ fastar2_1302_turbo-437c8f0c6c : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 8+0=8 terms (0+2+5+0+0+0+1) ~ 3.475059 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n4), (-1+n7), [(n9)]
* decreases: (-n9)
* +:         n3,     n5,         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4,         n7
`id' `Z'(n1?pos_,1,n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?{>=2},n8?neg0_,n9?neg_) =
+`Z'(n1,1,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n5,-1+n7)
+`Z'(n1,1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1-1,1,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n5,-1+n7)
+`Z'(n1,1,n3+1,n4-1,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n3,-1+n4)*rat(-1+n8+n3+ep,-1+n7)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,1,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-3+n7+n4+ep,-1+n4)
+`Z'(n1,1,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-9+n9+n8+2*n7+2*n6+n5+n4+2*n1+7*ep,-1+n7)
+`Z'(n1-1,1,n3,n4-1,n5,n6,n7-1,n8+2,n9+1)*rat(1,8)*rat(n8,-1+n4)*rat(1+n8,-1+n7)
;
*--#] fastar2_1302_turbo-437c8f0c6c : 
*--#[ fastar2_1302_turbo-cf90c089e9 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 7+0=7 terms (0+3+1+2+0+1) ~ 3.487540 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (n3+n4)
* +:                         n7, n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5
`id' `Z'(n1?pos_,1,n3?{>=2},n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,1,n3-1,n4-1,n5,n6,n7+1,n8,n9)*rat(-n7,-1+n3)*rat(-1+n7+ep,-1+n4)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8,n9)*rat(-3+n8+n4+n3+ep,-1+n3)
+`Z'(n1,1,n3,n4-1,n5,n6,n7,n8,n9)*rat(-3+n8+n4+n3+ep,-1+n4)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n3)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n3)
+`Z'(n1,1,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(1,2)*rat(n8,-1+n3)
+`Z'(n1-1,1,n3-1,n4-1,n5,n6,n7,n8+2,n9)*rat(1,4)*rat(n8,-1+n3)*rat(1+n8,-1+n4)
;
*--#] fastar2_1302_turbo-cf90c089e9 : 
*--#[ fastar2_1302_turbo-7469b3e5c2 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 7+0=7 terms (0+3+1+2+0+1) ~ 3.487540 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n1), (-1+n5)
* decreases: (n1+n5)
* +:                     n6,     n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5
`id' `Z'(n1?{>=2},1,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1-1,1,n3,n4,n5-1,n6+1,n7,n8,n9)*rat(-n6,-1+n1)*rat(-1+n6+ep,-1+n5)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8,n9)*rat(-3+n8+n5+n1+ep,-1+n1)
+`Z'(n1,1,n3,n4,n5-1,n6,n7,n8,n9)*rat(-3+n8+n5+n1+ep,-1+n5)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n1)
+`Z'(n1-1,1,n3-1,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(n8,-1+n1)
+`Z'(n1-1,1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(n8,-1+n1)
+`Z'(n1-1,1,n3-1,n4,n5-1,n6,n7,n8+2,n9)*rat(1,4)*rat(n8,-1+n1)*rat(1+n8,-1+n5)
;
*--#] fastar2_1302_turbo-7469b3e5c2 : 
#endif
#if `TURBOTHRESHOLD' >= 5
*--#[ fastar2_1302_turbo-2b842ded1e : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+3+5+1+2) ~ 4.292055 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n6), [(n9)]
* decreases: (-n9)
* +: n1,             n5,         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4, n5, n6
`id' `Z'(n1?pos_,1,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg0_,n9?neg_) =
+`Z'(n1,1,n3,n4,n5+1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n5,-1+n6)
+`Z'(n1+1,1,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n6)
+`Z'(n1,1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1-1,1,n3,n4,n5+1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n5,-1+n6)
+`Z'(n1+1,1,n3,n4,n5-1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n1,-1+n6)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,1,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,1,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n5-n1,-1+n6)
+`Z'(n1,1,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n6)
+`Z'(n1,1,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n6)
+`Z'(n1,1,n3,n4-1,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n8,-1+n6)
;
*--#] fastar2_1302_turbo-2b842ded1e : 
*--#[ fastar2_1302_turbo-4d60a5d71c : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+3+5+1+2) ~ 4.292055 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n7), [(n9)]
* decreases: (-n9)
* +:         n3,     n5,         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4, n5,     n7
`id' `Z'(n1?pos_,1,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg0_,n9?neg_) =
+`Z'(n1,1,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n5,-1+n7)
+`Z'(n1,1,n3+1,n4,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n3,-1+n7)
+`Z'(n1,1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1-1,1,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n5,-1+n7)
+`Z'(n1,1,n3+1,n4-1,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n3,-1+n7)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,1,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,1,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-11+n9+2*n8+2*n7+2*n6+n5+2*n4+n3+2*n1+8*ep,-1+n7)
+`Z'(n1,1,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n7)
+`Z'(n1-1,1,n3,n4,n5,n6,n7-1,n8+1,n9+1)*rat(1,4)*rat(n8,-1+n7)
+`Z'(n1,1,n3,n4,n5-1,n6,n7-1,n8+1,n9+1)*rat(-1,4)*rat(n8,-1+n7)
;
*--#] fastar2_1302_turbo-4d60a5d71c : 
*--#[ fastar2_1302_turbo-3211babf81 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 8+0=8 terms (0+4+1+2+0+1) ~ 4.367357 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n1), (-1+n4)
* decreases: (n3+n4)
* +:                 n5, n6,     n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4
`id' `Z'(n1?{>=2},1,n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1-1,1,n3,n4-1,n5,n6+1,n7,n8,n9)*rat(n6,-1+n1)*rat(-1+n6+ep,-1+n4)
+`Z'(n1-1,1,n3,n4-1,n5+1,n6,n7,n8,n9)*rat(-n5,-1+n1)*rat(-1+n8+n5+ep,-1+n4)
+`Z'(n1,1,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,1,n3,n4-1,n5,n6,n7,n8,n9)*rat(-9+n9+n8+2*n7+2*n6+n4+2*n3+n1+7*ep,-1+n4)
+`Z'(n1-1,1,n3,n4-1,n5+1,n6,n7,n8+1,n9)*rat(1,2)*rat(n5,-1+n1)*rat(n8,-1+n4)
+`Z'(n1-1,1,n3-1,n4-1,n5+1,n6,n7,n8+1,n9)*rat(1,2)*rat(n5,-1+n1)*rat(n8,-1+n4)
+`Z'(n1-1,1,n3,n4-2,n5+1,n6,n7,n8+1,n9)*rat(-1,2)*rat(n5,-1+n1)*rat(n8,-1+n4)
+`Z'(n1-1,1,n3-1,n4-1,n5,n6,n7,n8+2,n9)*rat(-1,4)*rat(n8,-1+n1)*rat(1+n8,-1+n4)
;
*--#] fastar2_1302_turbo-3211babf81 : 
*--#[ fastar2_1302_turbo-9603a7b449 : Z(n1,1,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,1+ep,n3,n4,n5,n6,n7,n8,n9) gives 8+0=8 terms (0+4+1+2+0+1) ~ 4.367357 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n3), (-1+n5)
* decreases: (n1+n5)
* +:             n4,         n7, n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3,     n5
`id' `Z'(n1?pos_,1,n3?{>=2},n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,1,n3-1,n4,n5-1,n6,n7+1,n8,n9)*rat(n7,-1+n3)*rat(-1+n7+ep,-1+n5)
+`Z'(n1,1,n3-1,n4+1,n5-1,n6,n7,n8,n9)*rat(-n4,-1+n3)*rat(-1+n8+n4+ep,-1+n5)
+`Z'(n1-1,1,n3,n4,n5,n6,n7,n8,n9)
+`Z'(n1,1,n3,n4,n5-1,n6,n7,n8,n9)*rat(-9+n9+n8+2*n7+2*n6+n5+n3+2*n1+7*ep,-1+n5)
+`Z'(n1,1,n3-1,n4+1,n5-1,n6,n7,n8+1,n9)*rat(1,2)*rat(n4,-1+n3)*rat(n8,-1+n5)
+`Z'(n1-1,1,n3-1,n4+1,n5-1,n6,n7,n8+1,n9)*rat(1,2)*rat(n4,-1+n3)*rat(n8,-1+n5)
+`Z'(n1,1,n3-1,n4+1,n5-2,n6,n7,n8+1,n9)*rat(-1,2)*rat(n4,-1+n3)*rat(n8,-1+n5)
+`Z'(n1-1,1,n3-1,n4,n5-1,n6,n7,n8+2,n9)*rat(-1,4)*rat(n8,-1+n3)*rat(1+n8,-1+n5)
;
*--#] fastar2_1302_turbo-9603a7b449 : 
#endif
*--#] n2 = 1 : 
#endprocedure
*--#] turbo : 
*--#[ reduce4567 :
#procedure reduce4567()
* Always applicable. Decrease the complexity. n7, n6, n5, n4 may be increased,
* respectively.
*--#[ fastar2_1302_p4 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 4+0=4 terms (0+4) ~ 4.000000 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                         n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2, n3, n4
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-5+n8+n7+n4+2*n3+2*ep,-1+n4)
;
*--#] fastar2_1302_p4 : 
*--#[ fastar2_1302_p5 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 4+0=4 terms (0+4) ~ 4.000000 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n5)
* decreases: (n1+n5)
* +:                     n6
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2,         n5
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7,n8,n9)*rat(n6,-1+n5)
+`Z'(n1,n2-1,n3,n4,n5-1,n6+1,n7,n8,n9)*rat(-n6,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9)*rat(-5+n8+n6+n5+2*n1+2*ep,-1+n5)
;
*--#] fastar2_1302_p5 : 
*--#[ fastar2_1302_m6 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 4+0=4 terms (0+1+3) ~ 2.311676 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n6)
* decreases: (-n2), (n1+n6)
* +:     n2,         n5
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -: n1,                 n6
`id' `Z'(n1?pos_,n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?neg0_,n9?neg0_) =
+`Z'(n1,n2+1,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-n5,-1+n6)
+`Z'(n1-1,n2+1,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(n5,-1+n6)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7,n8,n9)*rat(-5+n8+n6+n5+2*n1+2*ep,-1+n6)
;
*--#] fastar2_1302_m6 : 
*--#[ fastar2_1302_m7 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 4+0=4 terms (0+1+3) ~ 2.311676 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n7)
* decreases: (-n2), (n3+n7)
* +:     n2,     n4
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -:         n3,             n7
`id' `Z'(n1?pos_,n2?neg0_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?neg0_,n9?neg0_) =
+`Z'(n1,n2+1,n3,n4+1,n5,n6,n7-1,n8,n9)*rat(-n4,-1+n7)
+`Z'(n1,n2+1,n3-1,n4+1,n5,n6,n7-1,n8,n9)*rat(n4,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8,n9)*rat(-5+n8+n7+n4+2*n3+2*ep,-1+n7)
;
*--#] fastar2_1302_m7 : 
#endprocedure
*--#] reduce4567 : 
*--#[ reduce89 :
#procedure reduce89()
  #call turbo()
  #call reduce4567()
* Remove the numerators. n2 may be decreased.
*--#[ fastar2_1302_s8 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+2+4+4+1) ~ 3.528956 [+,*,+,+,+,+,+,-,-]
* checked poles: (-3+n8+n7+n4+n3+2*ep), [(n8)]
* decreases: (-n8)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7,     n9
* -: n1, n2, n3,     n5, n6
`id' `Z'(n1?pos_,n2?,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(-4+n9+n7+n6+n2+3*ep,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(1,2)*rat(n9,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9)*rat(-1,2)*rat(-5+n9+2*n6+2*n5-n3+n1+2*ep,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8+1,n9)*rat(-1,2)*rat(-4+n8+n7+n6+n3+n2+n1+3*ep,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9)*rat(1,2)*rat(-4+n9+n7+n6+n2+3*ep,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n9,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(-1,4)*rat(n9,-3+n8+n7+n4+n3+2*ep)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+2,n9)*rat(1,2)*rat(1+n8,-3+n8+n7+n4+n3+2*ep)
;
*--#] fastar2_1302_s8 : 
*--#[ fastar2_1302_s9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,n8,n9) gives 11+0=11 terms (0+1+5+2+3) ~ 3.030558 [+,*,+,+,+,+,+,-,-]
* checked poles: (-3+n9+n7+n6+n2+3*ep), [(n9)]
* decreases: (-n9)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5
`id' `Z'(n1?pos_,n2?,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?neg0_,n9?neg_) =
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n7+n6-n5-n4+n3+2*n2+n1+4*ep,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,2)*rat(1+n9,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-4+n8+n6+n5+n1+2*ep,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-4+n8+n7+n4+n3+2*ep,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-4+n8+n7+n4+n3+2*ep,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)*rat(-4+n8+n6+n5+n1+2*ep,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(1,4)*rat(n8,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n8,-3+n9+n7+n6+n2+3*ep)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9+1)*rat(-1,4)*rat(n8,-3+n9+n7+n6+n2+3*ep)
;
*--#] fastar2_1302_s9 : 
#endprocedure
*--#] reduce89 : 
*--#[ reduce67 :
#procedure reduce67()
  #call reduce4567()
* Reduce n6 and n7. n2 is increased. n7 and n6 may be increased, respectively,
* so the two rules have to be grouped.
*--#[ fastar2_1302_s6 : Z(n1,n2,n3,n4,n5,n6,n7,0,0)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,0,0) gives 3+0=3 terms (0+3) ~ 3.000000 [+,*,+,+,+,+,+,-,-]
* checked poles: (-1+n6), (-3+n7+n6+n2+3*ep)
* decreases: (n1+n6)
* +:     n2,                 n7
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3,         n6
`id' `Z'(n1?pos_,n2?,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,0,0) =
+`Z'(n1,n2+1,n3-1,n4,n5,n6-1,n7+1,0,0)*rat(n7,-1+n6)*rat(-1+n7+ep,-3+n7+n6+n2+3*ep)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,0,0)*rat(-1+n7+n2+2*ep,-3+n7+n6+n2+3*ep)
+`Z'(n1,n2+1,n3,n4,n5,n6-1,n7,0,0)*rat(1-n7-n2-2*ep,-1+n6)
;
*--#] fastar2_1302_s6 : 
*--#[ fastar2_1302_s7 : Z(n1,n2,n3,n4,n5,n6,n7,0,0)
* Z(n1,n2+ep,n3,n4,n5,n6,n7,0,0) gives 3+0=3 terms (0+3) ~ 3.000000 [+,*,+,+,+,+,+,-,-]
* checked poles: (-1+n7), (-3+n7+n6+n2+3*ep)
* decreases: (n3+n7)
* +:     n2,             n6
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3,             n7
`id' `Z'(n1?pos_,n2?,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},0,0) =
+`Z'(n1-1,n2+1,n3,n4,n5,n6+1,n7-1,0,0)*rat(n6,-1+n7)*rat(-1+n6+ep,-3+n7+n6+n2+3*ep)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7,0,0)*rat(-1+n6+n2+2*ep,-3+n7+n6+n2+3*ep)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,0,0)*rat(1-n6-n2-2*ep,-1+n7)
;
*--#] fastar2_1302_s7 : 
#endprocedure
*--#] reduce67 : 
*--#[ reduce45 :
#procedure reduce45()
* Reduce n4 and n5. No side effect.
*--#[ fastar2_1302_x4 : Z(n1,n2,n3,n4,n5,1,1,0,0)
* Z(n1,n2+ep,n3,n4,n5,1,1,0,0) gives 6+0=6 terms (0+2+4) ~ 3.251678 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n4), (-3+n5+n4+2*ep)
* decreases: (n3+n4)
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5
`id' `Z'(n1?pos_,n2?,n3?pos_,n4?{>=2},n5?pos_,1,1,0,0) =
+`Z'(n1,n2,n3-1,n4,n5,1,1,0,0)
+`Z'(n1,n2,n3,n4-1,n5,1,1,0,0)*rat(-1,-3+n5+n4+2*ep)*rat(-19+2*n5+n5^2+13*n4-n4*n5-2*n4^2+3*n3-2*n3*n5-2*n3*n4+n3^2+3*n2-n2*n5-n2*n4+n2*n3+3*n1-n1*n2-n1^2+24*ep-2*ep*n5-8*ep*n4-2*ep*n2-4*ep*n1-8*ep^2,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,1,1,0,0)*rat(3-n5-n1-2*ep,-1+n4)*rat(-4+n5+n2+n1+4*ep,-3+n5+n4+2*ep)
+`Z'(n1,n2,n3-1,n4-1,n5,1,1,0,0)*rat(-4+n4+n3+2*ep,-1+n4)*rat(-5+n4+n3+n2+4*ep,-3+n5+n4+2*ep)
+`Z'(n1,n2,n3,n4-2,n5,1,1,0,0)*rat(4-n4-n3-2*ep,-1+n4)*rat(-5+n4+n3+n2+4*ep,-3+n5+n4+2*ep)
+`Z'(n1,n2,n3,n4-1,n5-1,1,1,0,0)*rat(-3+n5+n1+2*ep,-1+n4)*rat(-4+n5+n2+n1+4*ep,-3+n5+n4+2*ep)
;
*--#] fastar2_1302_x4 : 
*--#[ fastar2_1302_x5 : Z(n1,n2,n3,n4,n5,1,1,0,0)
* Z(n1,n2+ep,n3,n4,n5,1,1,0,0) gives 6+0=6 terms (0+2+4) ~ 3.251678 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n5), (-3+n5+n4+2*ep), [(-2+n2+3*ep), (-3+n5+n1+ep)]
* decreases: (n1+n5)
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5
`id' `Z'(n1?pos_,n2?,n3?pos_,n4?pos_,n5?{>=2},1,1,0,0) =
+`Z'(n1-1,n2,n3,n4,n5,1,1,0,0)
+`Z'(n1,n2,n3,n4,n5-1,1,1,0,0)*rat(1,-3+n5+n4+2*ep)*rat(19-13*n5+2*n5^2-2*n4+n4*n5-n4^2-3*n3+n3^2-3*n2+n2*n5+n2*n4+n2*n3-3*n1+2*n1*n5+2*n1*n4-n1*n2-n1^2-24*ep+8*ep*n5+2*ep*n4+4*ep*n3+2*ep*n2+8*ep^2,-1+n5)
+`Z'(n1-1,n2,n3,n4,n5-1,1,1,0,0)*rat(-4+n5+n1+2*ep,-1+n5)*rat(-5+n5+n2+n1+4*ep,-3+n5+n4+2*ep)
+`Z'(n1,n2,n3-1,n4,n5-1,1,1,0,0)*rat(3-n4-n3-2*ep,-1+n5)*rat(-4+n4+n3+n2+4*ep,-3+n5+n4+2*ep)
+`Z'(n1,n2,n3,n4-1,n5-1,1,1,0,0)*rat(-3+n4+n3+2*ep,-1+n5)*rat(-4+n4+n3+n2+4*ep,-3+n5+n4+2*ep)
+`Z'(n1,n2,n3,n4,n5-2,1,1,0,0)*rat(4-n5-n1-2*ep,-1+n5)*rat(-5+n5+n2+n1+4*ep,-3+n5+n4+2*ep)
;
*--#] fastar2_1302_x5 : 
#endprocedure
*--#] reduce45 : 
*--#[ reduce13 :
#procedure reduce13()
* Reduce n1 and n3. n2 may be increased.
*--#[ fastar2_1302_x1 : Z(n1,n2,n3,1,1,1,1,0,0)
* Z(n1,n2+ep,n3,1,1,1,1,0,0) gives 3+2=5 terms (0+3)+(0+2) ~ 3.020867 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n1), (-3+n3+n2+n1+3*ep), [(-1+n2+2*ep), (-2+n1+2*ep)]
* decreases: (n1+n3)
* +:     n2
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3
`id' `Z'(n1?{>=2},n2?,n3?pos_,1,1,1,1,0,0) =
+`Z'(n1-1,n2+1,n3-1,1,1,1,1,0,0)*rat(-n2-ep,-1+n1)*rat(-2+n1+2*ep,-3+n3+n2+n1+3*ep)
+`Z'(n1-1,n2,n3,1,1,1,1,0,0)*rat(-2+n1+2*ep,-1+n1)*rat(-3+n2+n1+4*ep,-3+n3+n2+n1+3*ep)
+`Z'(n1,n2,n3-1,1,1,1,1,0,0)*rat(-3+n3+n2+4*ep,-3+n3+n2+n1+3*ep)
+`Y'(n1,n2,n3,0,1,1,1,0,0)*rat(3-n3-n2-4*ep,-1+n1)
+`Y'(n1,n2,n3,1,1,1,0,0,0)*rat(-2+n2+3*ep,-1+n1)*rat(-3+n2+n1+4*ep,-3+n3+n2+n1+3*ep)
;
*--#] fastar2_1302_x1 : 
*--#[ fastar2_1302_x3 : Z(n1,n2,n3,1,1,1,1,0,0)
* Z(n1,n2+ep,n3,1,1,1,1,0,0) gives 3+2=5 terms (0+3)+(0+2) ~ 3.020867 [+,x,+,+,+,+,+,-,-]
* checked poles: (-1+n3), (-3+n3+n2+n1+3*ep), [(-1+n2+2*ep), (-2+n3+2*ep)]
* decreases: (n1+n3)
* +:     n2
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3
`id' `Z'(n1?pos_,n2?,n3?{>=2},1,1,1,1,0,0) =
+`Z'(n1-1,n2+1,n3-1,1,1,1,1,0,0)*rat(-n2-ep,-1+n3)*rat(-2+n3+2*ep,-3+n3+n2+n1+3*ep)
+`Z'(n1-1,n2,n3,1,1,1,1,0,0)*rat(-3+n2+n1+4*ep,-3+n3+n2+n1+3*ep)
+`Z'(n1,n2,n3-1,1,1,1,1,0,0)*rat(-2+n3+2*ep,-1+n3)*rat(-3+n3+n2+4*ep,-3+n3+n2+n1+3*ep)
+`Y'(n1,n2,n3,1,0,1,1,0,0)*rat(3-n2-n1-4*ep,-1+n3)
+`Y'(n1,n2,n3,1,1,0,1,0,0)*rat(-2+n2+3*ep,-1+n3)*rat(-3+n3+n2+4*ep,-3+n3+n2+n1+3*ep)
;
*--#] fastar2_1302_x3 : 
#endprocedure
*--#] reduce13 : 
*--#[ reduce2 :
#procedure reduce2()
* Reduce n2.
*--#[ fastar2_1302_p2 : Z(1,n2,1,1,1,1,1,0,0)
* Z(1,n2+ep,1,1,1,1,1,0,0) gives 1+4=5 terms (0+1)+(0+4) ~ 1.061772 [+,+,+,+,+,+,+,-,-]
* checked poles: (-1+n2+2*ep), (-1+n2+3*ep)
* decreases: (n2)
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
* -:     n2
`id' `Z'(1,n2?{>=2},1,1,1,1,1,0,0) =
+`Z'(1,n2-1,1,1,1,1,1,0,0)*rat(-2+n2+3*ep,-1+n2+2*ep)*rat(-2+n2+4*ep,-1+n2+3*ep)
+`Y'(1,n2-1,1,1,0,2,1,0,0)*rat(-1,-1+n2+2*ep)
+`Y'(0,n2,1,1,1,1,1,0,0)*rat(-1+2*ep,-1+n2+3*ep)
+`Y'(1,n2,0,1,1,1,1,0,0)*rat(-1+2*ep,-1+n2+2*ep)*rat(-2+n2+4*ep,-1+n2+3*ep)
+`Y'(1,n2,1,1,1,1,0,0,0)*rat(1-2*ep,-1+n2+2*ep)*rat(-2+n2+3*ep,-1+n2+3*ep)
;
*--#] fastar2_1302_p2 : 
*--#[ fastar2_1302_m2 : Z(1,n2,1,1,1,1,1,0,0)
* Z(1,n2+ep,1,1,1,1,1,0,0) gives 1+4=5 terms (0+1)+(0+0+4) ~ 1.061772 [+,-,+,+,+,+,+,-,-]
* checked poles: (-1+n2+3*ep), (-1+n2+4*ep)
* decreases: (-n2)
* +:     n2
* 0: n1,     n3, n4, n5, n6, n7, n8, n9
`id' `Z'(1,n2?neg0_,1,1,1,1,1,0,0) =
+`Z'(1,n2+1,1,1,1,1,1,0,0)*rat(n2+2*ep,-1+n2+3*ep)*rat(n2+3*ep,-1+n2+4*ep)
+`Y'(0,n2+1,1,1,1,1,1,0,0)*rat(1-2*ep,-1+n2+3*ep)
+`Y'(1,n2+1,0,1,1,1,1,0,0)*rat(1-2*ep,-1+n2+3*ep)*rat(n2+2*ep,-1+n2+4*ep)
+`Y'(1,n2+1,1,0,1,1,1,0,0)*rat(-1+2*ep,-1+n2+3*ep)*rat(n2+3*ep,-1+n2+4*ep)
+`Y'(1,n2+1,1,1,1,0,1,0,0)*rat(-1+2*ep,-1+n2+4*ep)
;
*--#] fastar2_1302_m2 : 
#endprocedure
*--#] reduce2 : 
  #call StdReduceProc(turbo,89,67,45,13,2)
#endprocedure
* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4
