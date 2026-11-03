#procedure SchemeZlastar5()
  #ifndef `LASTAR5SCHEME'
    #define LASTAR5SCHEME "3"
  #endif
  #ifndef `TURBOTHRESHOLD'
    #define TURBOTHRESHOLD "4"
  #endif
*--#[ def :
* internals: p6,p7,p8
* externals: p0
* propagators: p1.p1,p2.p2,p3.p3,p4.p4,p5.p5,p6.p6,p7.p7,p8.p8,p1.p3
* replacements:
*   p3 = p0 + p6 + p7 + p8,
*   p5 = p6 + p7,
*   p1 = p0 + p6,
*   p9 = p0,
*   p2 = p0 + p6 + p7,
*   p4 = -p6 - p7 - p8,
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
*--#[ lastar5_13_turbo-45f39a35a2 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+0+3+2+1) ~ 2.066852 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3), [(n9)]
* decreases: (n3), (-n9)
* +:     n2,                         n9
* 0: n1, n2,     n4, n5, n6, n7, n8
* -: n1,     n3,         n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n2,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n8+n7+n6+n5+n4+n3-n1+3*ep,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n3)
;
*--#] lastar5_13_turbo-45f39a35a2 : 
*--#[ lastar5_13_turbo-2ddb1221c7 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+0+3+2+1) ~ 2.066852 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1), [(n9)]
* decreases: (n1), (-n9)
* +:     n2,                         n9
* 0:     n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4,             n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n2,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n8+n7+n6+n5+n4-n3+n1+3*ep,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n1)
;
*--#] lastar5_13_turbo-2ddb1221c7 : 
*--#[ lastar5_13_turbo-6ca7a0f36c : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+4+4) ~ 2.401892 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n9), (n1+n3+n4)
* +:                 n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n5+ep,-1+n4)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5+1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n5+ep,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-9+n9+n8+n7+2*n6+n5+n4+n3+n2+n1+5*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n4)
;
*--#] lastar5_13_turbo-6ca7a0f36c : 
*--#[ lastar5_13_turbo-9f3ad2650a : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+4+4) ~ 2.401892 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n6), [(n9)]
* decreases: (-n9), (n1+n3+n6)
* +:                 n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8-1,n9+1)*rat(1,2)*rat(n5+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n6)
+`Z'(n1,n2-1,n3,n4-1,n5+1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n5+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-9+n9+n8+n7+n6+n5+2*n4+n3+n2+n1+5*ep,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n6)
;
*--#] lastar5_13_turbo-9f3ad2650a : 
*--#[ lastar5_13_turbo-34d818f4c9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n8), [(n9)]
* decreases: (-n9), (n1+n3+n7+n8)
* +:     n2,         n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3,         n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?{>=2},n9?neg_) =
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8-1,n9+1)*rat(1,2)*rat(n5+ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(-1,2)*rat(n5+ep,-1+n8)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8-1,n9+1)*rat(-1,2)*rat(n2,-1+n8)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(-5+n8+2*n7+n5+n2+3*ep,-1+n8)
;
*--#] lastar5_13_turbo-34d818f4c9 : 
*--#[ lastar5_13_turbo-8b682f6f45 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n7), [(n9)]
* decreases: (-n9), (n1+n3+n7+n8)
* +:     n2,         n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4,         n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?{>=2},n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n5+ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(-1,2)*rat(n5+ep,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n2,-1+n7)
+`Z'(n1,n2+1,n3,n4,n5,n6,n7-1,n8-1,n9+1)*rat(-1,2)*rat(n2,-1+n7)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-5+2*n8+n7+n5+n2+3*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] lastar5_13_turbo-8b682f6f45 : 
#endif
#if `TURBOTHRESHOLD' >= 5
*--#[ lastar5_13_turbo-3580abd1f4 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 4+0=4 terms (0+4) ~ 4.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                             n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2, n3, n4
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-n8,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(n8,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-5+n9+n8+n4+2*n3+2*ep,-1+n4)
;
*--#] lastar5_13_turbo-3580abd1f4 : 
*--#[ lastar5_13_turbo-595ac180b9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 4+0=4 terms (0+4) ~ 4.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n6)
* decreases: (n1+n6)
* +:                         n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2,             n6
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5,n6-1,n7+1,n8,n9)*rat(n7,-1+n6)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7+1,n8,n9)*rat(-n7,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9)*rat(-5+n9+n7+n6+2*n1+2*ep,-1+n6)
;
*--#] lastar5_13_turbo-595ac180b9 : 
*--#[ lastar5_13_turbo-5d9a35dd58 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 7+0=7 terms (0+4+1+2) ~ 4.358189 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1)
* decreases: (n1+n6)
* +:                         n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,n9)*rat(-n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6-1,n7+1,n8,n9)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)*rat(-5+n9+n7+2*n6+n1+2*ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n1)
;
*--#] lastar5_13_turbo-5d9a35dd58 : 
*--#[ lastar5_13_turbo-cd3093aade : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 7+0=7 terms (0+4+1+2) ~ 4.358189 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3)
* decreases: (n3+n4)
* +:                             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,n9)*rat(-n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)*rat(-5+n9+n8+2*n4+n3+2*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n3)
;
*--#] lastar5_13_turbo-cd3093aade : 
#endif
#if `TURBOTHRESHOLD' >= 7
*--#[ lastar5_13_turbo-073a6f3727 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n2)
* decreases: (n2+n5)
* +:                         n7, n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,     n4, n5, n6
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-n8,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7,n8+1,n9)*rat(n8,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5-1,n6,n7+1,n8,n9)*rat(n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7+1,n8,n9)*rat(-n7,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9)*rat(-5+n8+n7+2*n5+n2+4*ep,-1+n2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9)
;
*--#] lastar5_13_turbo-073a6f3727 : 
#endif
*--#] n5 >= 2 : 
*--#[ n5 <= 0 :
#if `TURBOTHRESHOLD' >= 2
*--#[ lastar5_13_turbo-8bbbda602a : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 4+0=4 terms (0+0+3+0+0+0+1) ~ 1.766880 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n4), (-1+n8)
* decreases: (-n5), (n4+n8)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -: n1,         n4,             n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?neg0_,n6?pos_,n7?pos_,n8?{>=2},n9?neg0_) =
+`Z'(n1,n2,n3+1,n4-1,n5+1,n6,n7,n8-1,n9)*rat(-n3,-1+n4)*rat(-1+n9+n3+ep,-1+n8)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8,n9)*rat(-3+n8+n4+ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9)*rat(-3+n8+n4+ep,-1+n8)
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+2)*rat(-1,4)*rat(n9,-1+n4)*rat(1+n9,-1+n8)
;
*--#] lastar5_13_turbo-8bbbda602a : 
*--#[ lastar5_13_turbo-d390d9f640 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 4+0=4 terms (0+0+3+0+0+0+1) ~ 1.766880 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n6), (-1+n7), [(n2)]
* decreases: (-n5), (n6+n7)
* +: n1,             n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -:         n3,         n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?{>=2},n7?{>=2},n8?pos_,n9?neg0_) =
+`Z'(n1+1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9)*rat(-n1,-1+n6)*rat(-1+n9+n1+ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-3+n7+n6+ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9)*rat(-3+n7+n6+ep,-1+n7)
+`Z'(n1,n2,n3-1,n4,n5+1,n6-1,n7-1,n8,n9+2)*rat(-1,4)*rat(n9,-1+n6)*rat(1+n9,-1+n7)
;
*--#] lastar5_13_turbo-d390d9f640 : 
*--#[ lastar5_13_turbo-0885502fe5 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 11+0=11 terms (0+0+1+2+6+0+2) ~ 1.963449 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n5), (n2-n9), (n4+n8), (n4-n9)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?neg0_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3+1,n4-1,n5+1,n6,n7,n8,n9)*rat(-n3,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9+1)*rat(1,2)
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(-1,2)*rat(-3+n9-n8+2*n7+2*n6+n5-n3+n1+3*ep,-1+n4)
+`Z'(n1-1,n2,n3,n4,n5+1,n6,n7,n8-1,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(-1,2)*rat(-2+n9+n7+n5+n2+3*ep,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(-2+n9+n7+n5+n2+3*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8-1,n9+1)*rat(-1,2)
+`Z'(n1-1,n2,n3-1,n4-1,n5+1,n6,n7,n8,n9+2)*rat(1,2)*rat(1+n9,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+2)*rat(-1,2)*rat(1+n9,-1+n4)
;
*--#] lastar5_13_turbo-0885502fe5 : 
#endif
#if `TURBOTHRESHOLD' >= 3
*--#[ lastar5_13_turbo-b22309e7f8 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+0+3+2+1) ~ 2.066852 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n1), [(n9)]
* decreases: (n1), (-n9)
* +:     n2,                         n9
* 0:     n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4,             n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n2,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n8+n7+n6+n5+n4-n3+n1+3*ep,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n1)
;
*--#] lastar5_13_turbo-b22309e7f8 : 
*--#[ lastar5_13_turbo-879c470934 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+0+3+2+1) ~ 2.066852 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n3), [(n9)]
* decreases: (n3), (-n9)
* +:     n2,                         n9
* 0: n1, n2,     n4, n5, n6, n7, n8
* -: n1,     n3,         n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?neg0_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n2,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n8+n7+n6+n5+n4+n3-n1+3*ep,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-8+n9+n8+n7+n6+n5+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n3)
;
*--#] lastar5_13_turbo-879c470934 : 
*--#[ lastar5_13_turbo-50d4590fde : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 12+0=12 terms (0+0+1+3+8) ~ 2.098995 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n7), [(n9)]
* decreases: (-n9), (n3-n5), (-n5+n7), (n6+n7), (n2+n3+n7+n8)
* +:     n2,         n5,         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3, n4,     n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?{>=2},n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3-1,n4,n5,n6,n7-1,n8+1,n9+1)*rat(-n8,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n2,-1+n7)
+`Z'(n1,n2+1,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(-1,2)*rat(n2,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n5+n4-n2+ep,-1+n7)
+`Z'(n1,n2,n3-1,n4-1,n5+1,n6,n7-1,n8+1,n9+1)*rat(n8,-1+n7)
+`Z'(n1,n2-1,n3,n4,n5+1,n6-1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-4+n8+n7+n5+n2+3*ep,-1+n7)
+`Z'(n1,n2,n3-1,n4,n5+1,n6-1,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-2+n8+n7+2*n5-n4+n2+4*ep,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(-3+n9+n5+n4+2*n3+3*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8-1,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(-4+n8+n7+n5+n2+3*ep,-1+n7)
;
*--#] lastar5_13_turbo-50d4590fde : 
*--#[ lastar5_13_turbo-e40ae59461 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 11+0=11 terms (0+1+0+1+7+0+2) ~ 2.099631 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n8), [(n9)]
* decreases: (-n5), (n4+n8), (n8-n9), (n1+n2+n7+n8)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?pos_,n8?{>=2},n9?neg_) =
+`Z'(n1,n2,n3+1,n4,n5+1,n6,n7,n8-1,n9)*rat(-n3,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9+1)*rat(1,2)*rat(-6+n8+n7+n4+n3+n2+n1+2*ep,-1+n8)
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1-1,n2,n3,n4,n5+1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(-3+n9+2*n7+2*n6+n5-n4-n3+n1+3*ep,-1+n8)
+`Z'(n1,n2-1,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4,n5+1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(-4+n8+n7+n5+n2+3*ep,-1+n8)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8-1,n9+1)*rat(-1,2)*rat(-6+n8+n7+n4+n3+n2+n1+2*ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(-4+n8+n7+n5+n2+3*ep,-1+n8)
+`Z'(n1-1,n2,n3-1,n4,n5+1,n6,n7,n8-1,n9+2)*rat(1,2)*rat(1+n9,-1+n8)
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+2)*rat(-1,2)*rat(1+n9,-1+n8)
;
*--#] lastar5_13_turbo-e40ae59461 : 
*--#[ lastar5_13_turbo-c0d27c20d5 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 7+0=7 terms (0+0+4+0+3) ~ 2.163348 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n6)
* decreases: (-n5), (n6+n7)
* +: n1,             n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -:     n2, n3,         n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1+1,n2-1,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-n1,-1+n6)
+`Z'(n1+1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9)*rat(n1,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-5+n9+2*n7+n6+n1+2*ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9)
+`Z'(n1,n2-1,n3,n4,n5+1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5+1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8-1,n9+1)*rat(1,2)*rat(n9,-1+n6)
;
*--#] lastar5_13_turbo-c0d27c20d5 : 
*--#[ lastar5_13_turbo-137c8a4e70 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 7+0=7 terms (0+0+4+0+3) ~ 2.163348 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (-n5), (n4+n8)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -: n1, n2,     n4,         n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?neg0_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3+1,n4-1,n5+1,n6,n7,n8,n9)*rat(-n3,-1+n4)
+`Z'(n1,n2,n3+1,n4-1,n5+1,n6,n7,n8-1,n9)*rat(n3,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8,n9)*rat(-5+n9+2*n8+n4+n3+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9)
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n9,-1+n4)
;
*--#] lastar5_13_turbo-137c8a4e70 : 
*--#[ lastar5_13_turbo-aae2b09917 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 7+0=7 terms (0+1+3+1+2) ~ 2.513766 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n7)
* decreases: (-n5), (n6+n7)
* +: n1,             n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -:         n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?{>=2},n8?pos_,n9?neg0_) =
+`Z'(n1+1,n2,n3,n4,n5+1,n6,n7-1,n8,n9)*rat(-n1,-1+n7)
+`Z'(n1+1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9)*rat(n1,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9)*rat(-5+n9+n7+2*n6+n1+2*ep,-1+n7)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n9,-1+n7)
+`Z'(n1,n2,n3-1,n4,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n9,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n9,-1+n7)
;
*--#] lastar5_13_turbo-aae2b09917 : 
*--#[ lastar5_13_turbo-39e852b440 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 7+0=7 terms (0+1+3+1+2) ~ 2.513766 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n8)
* decreases: (-n5), (n4+n8)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -: n1,         n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?pos_,n8?{>=2},n9?neg0_) =
+`Z'(n1,n2,n3+1,n4,n5+1,n6,n7,n8-1,n9)*rat(-n3,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5+1,n6,n7,n8-1,n9)*rat(n3,-1+n8)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9)*rat(-5+n9+n8+2*n4+n3+2*ep,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n9,-1+n8)
+`Z'(n1-1,n2,n3,n4,n5+1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n9,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7,n8-1,n9+1)*rat(1,2)*rat(n9,-1+n8)
;
*--#] lastar5_13_turbo-39e852b440 : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ lastar5_13_turbo-0e8601ec30 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 16+0=16 terms (0+2+3+3+6+1+1) ~ 3.386047 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n2), [(n9)]
* decreases: (-n5)
* +: n1,     n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1+1,n2-1,n3,n4,n5+1,n6,n7,n8,n9)*rat(n1,-1+n2)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8,n9)
+`Z'(n1,n2-1,n3+1,n4-1,n5+1,n6,n7,n8,n9)*rat(-n3,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-n1,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5+1,n6,n7,n8,n9)*rat(8-n9-n8-2*n7-2*n6-2*n5-n4-n2-n1-6*ep,-1+n2)
+`Z'(n1-1,n2,n3,n4,n5+1,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2-1,n3,n4,n5+1,n6,n7,n8,n9+1)*rat(1,2)*rat(2+n9+n6+n5+n4-n2-2*n1+ep,-1+n2)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(1,2)
+`Z'(n1-1,n2-1,n3,n4,n5+1,n6,n7,n8,n9+1)*rat(1,2)*rat(-11+2*n9+n8+3*n7+3*n6+2*n5+n4+n2+2*n1+8*ep,-1+n2)
+`Z'(n1,n2-2,n3,n4,n5+1,n6,n7,n8,n9+1)*rat(1,2)*rat(-4+n8+n7+n5+n2+3*ep,-1+n2)
+`Z'(n1,n2-1,n3-1,n4,n5+1,n6,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5+1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-2+n9+n6+n5+n4+3*ep,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5+1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-4+n8+n7+n5+n2+3*ep,-1+n2)
+`Z'(n1-1,n2-1,n3,n4,n5+1,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n2)
+`Z'(n1-1,n2-1,n3-1,n4,n5+1,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n2)
;
*--#] lastar5_13_turbo-0e8601ec30 : 
*--#[ lastar5_13_turbo-3e9547241a : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 11+0=11 terms (0+3+3+1+4) ~ 3.940637 [+,+,+,+,-,+,+,+,-]
* checked poles: (-1+n2)
* decreases: (-n5)
* +: n1,     n3,     n5,             n9
* 0: n1, n2, n3, n4,     n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?neg0_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3+1,n4,n5+1,n6,n7,n8,n9)*rat(n3,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5+1,n6,n7,n8,n9)*rat(n1,-1+n2)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7,n8,n9)
+`Z'(n1,n2-1,n3+1,n4-1,n5+1,n6,n7,n8,n9)*rat(-n3,-1+n2)
+`Z'(n1+1,n2-1,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-n1,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5+1,n6,n7,n8,n9)*rat(11-2*n9-2*n8-2*n7-2*n6-2*n5-2*n4-n3-n2-n1-8*ep,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5+1,n6,n7,n8,n9+1)*rat(n9,-1+n2)
+`Z'(n1-1,n2-1,n3,n4,n5+1,n6,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n2)
+`Z'(n1,n2-1,n3-1,n4,n5+1,n6,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n2)
+`Z'(n1,n2-1,n3,n4-1,n5+1,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5+1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n2)
;
*--#] lastar5_13_turbo-3e9547241a : 
#endif
*--#] n5 <= 0 : 
*--#[ n5 = 1 :
#if `TURBOTHRESHOLD' >= 3
*--#[ lastar5_13_turbo-b09fbbde34 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+0+3+2+1) ~ 2.066852 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n1), [(n9)]
* decreases: (n1), (-n9)
* +:     n2,                         n9
* 0:     n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4,             n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,1,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2+1,n3-1,n4,1,n6,n7,n8,n9+1)*rat(-1,2)*rat(n2,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,1,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n1)
+`Z'(n1-1,n2,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)*rat(-4+n8+n7+n6+n4-n3+n1+3*ep,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,1,n6,n7,n8,n9+1)*rat(1,2)*rat(-7+n9+n8+n7+n6+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,1,n6,n7,n8,n9+1)*rat(-1,2)*rat(-7+n9+n8+n7+n6+n4+n3+n2+n1+5*ep,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,1,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n1)
;
*--#] lastar5_13_turbo-b09fbbde34 : 
*--#[ lastar5_13_turbo-b66c6a0868 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+0+3+2+1) ~ 2.066852 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n3), [(n9)]
* decreases: (n3), (-n9)
* +:     n2,                         n9
* 0: n1, n2,     n4, n5, n6, n7, n8
* -: n1,     n3,         n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,1,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2+1,n3-1,n4,1,n6,n7,n8,n9+1)*rat(-1,2)*rat(n2,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,1,n6,n7,n8,n9+1)*rat(1,2)*rat(-4+n8+n7+n6+n4+n3-n1+3*ep,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,1,n6,n7,n8,n9+1)*rat(1,2)*rat(-7+n9+n8+n7+n6+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1,n2,n3-1,n4,1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-7+n9+n8+n7+n6+n4+n3+n2+n1+5*ep,-1+n3)
+`Z'(n1-1,n2,n3-1,n4,1,n6,n7,n8,n9+2)*rat(-1,2)*rat(1+n9,-1+n3)
;
*--#] lastar5_13_turbo-b66c6a0868 : 
*--#[ lastar5_13_turbo-491a653403 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+4+4) ~ 2.401892 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n9), (n1+n3+n4)
* +:                 n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},1,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3,n4-1,2,n6,n7-1,n8,n9+1)*rat(1,2)*rat(1+ep,-1+n4)
+`Z'(n1-1,n2,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,1,n6,n7,n8,n9+1)*rat(1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,2,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(1+ep,-1+n4)
+`Z'(n1-1,n2,n3,n4-1,1,n6,n7,n8,n9+1)*rat(-1,2)*rat(-8+n9+n8+n7+2*n6+n4+n3+n2+n1+5*ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,1,n6-1,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4-1,1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n4)
;
*--#] lastar5_13_turbo-491a653403 : 
*--#[ lastar5_13_turbo-2c03dbea80 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+4+4) ~ 2.401892 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n6), [(n9)]
* decreases: (-n9), (n1+n3+n6)
* +:                 n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,1,n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3,n4,2,n6-1,n7,n8-1,n9+1)*rat(1,2)*rat(1+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,1,n6-1,n7,n8,n9+1)*rat(1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n6)
+`Z'(n1,n2-1,n3,n4-1,2,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(1+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4-1,1,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3-1,n4,1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-8+n9+n8+n7+n6+2*n4+n3+n2+n1+5*ep,-1+n6)
+`Z'(n1,n2,n3,n4-1,1,n6-1,n7,n8,n9+1)*rat(-1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-1+n6)
;
*--#] lastar5_13_turbo-2c03dbea80 : 
*--#[ lastar5_13_turbo-58d27c683e : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n7), [(n9)]
* decreases: (-n9), (n1+n3+n7+n8)
* +:     n2,         n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4,         n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,1,n6?pos_,n7?{>=2},n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3,n4-1,2,n6,n7-1,n8,n9+1)*rat(1,2)*rat(1+ep,-1+n7)
+`Z'(n1,n2,n3,n4,2,n6,n7-1,n8-1,n9+1)*rat(-1,2)*rat(1+ep,-1+n7)
+`Z'(n1,n2+1,n3-1,n4,1,n6,n7-1,n8,n9+1)*rat(1,2)*rat(n2,-1+n7)
+`Z'(n1,n2+1,n3,n4,1,n6,n7-1,n8-1,n9+1)*rat(-1,2)*rat(n2,-1+n7)
+`Z'(n1-1,n2,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,1,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(-4+2*n8+n7+n2+3*ep,-1+n7)
+`Z'(n1,n2,n3,n4,1,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] lastar5_13_turbo-58d27c683e : 
*--#[ lastar5_13_turbo-d7af3e9e45 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n8), [(n9)]
* decreases: (-n9), (n1+n3+n7+n8)
* +:     n2,         n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3,         n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,1,n6?pos_,n7?pos_,n8?{>=2},n9?neg_) =
+`Z'(n1-1,n2+1,n3,n4,1,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n2,-1+n8)
+`Z'(n1,n2,n3,n4,2,n6-1,n7,n8-1,n9+1)*rat(1,2)*rat(1+ep,-1+n8)
+`Z'(n1,n2,n3,n4,2,n6,n7-1,n8-1,n9+1)*rat(-1,2)*rat(1+ep,-1+n8)
+`Z'(n1,n2+1,n3,n4,1,n6,n7-1,n8-1,n9+1)*rat(-1,2)*rat(n2,-1+n8)
+`Z'(n1-1,n2,n3,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,1,n6,n7-1,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(-4+n8+2*n7+n2+3*ep,-1+n8)
;
*--#] lastar5_13_turbo-d7af3e9e45 : 
#endif
#if `TURBOTHRESHOLD' >= 5
*--#[ lastar5_13_turbo-daa93036aa : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 4+0=4 terms (0+4) ~ 4.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                             n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2, n3, n4
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},1,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4-1,1,n6,n7,n8+1,n9)*rat(-n8,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,1,n6,n7,n8+1,n9)*rat(n8,-1+n4)
+`Z'(n1,n2,n3-1,n4,1,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,1,n6,n7,n8,n9)*rat(-5+n9+n8+n4+2*n3+2*ep,-1+n4)
;
*--#] lastar5_13_turbo-daa93036aa : 
*--#[ lastar5_13_turbo-c57344ff18 : Z(n1,n2,n3,n4,1,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,1+ep,n6,n7,n8,n9) gives 4+0=4 terms (0+4) ~ 4.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n6)
* decreases: (n1+n6)
* +:                         n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2,             n6
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,1,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4,1,n6-1,n7+1,n8,n9)*rat(n7,-1+n6)
+`Z'(n1,n2-1,n3,n4,1,n6-1,n7+1,n8,n9)*rat(-n7,-1+n6)
+`Z'(n1-1,n2,n3,n4,1,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4,1,n6-1,n7,n8,n9)*rat(-5+n9+n7+n6+2*n1+2*ep,-1+n6)
;
*--#] lastar5_13_turbo-c57344ff18 : 
#endif
*--#] n5 = 1 : 
#endprocedure
*--#] turbo : 
*--#[ reduce4 :
#procedure reduce4()
*--#[ lastar5_13_x4 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 4+0=4 terms (0+4) ~ 4.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                             n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2, n3, n4
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7,n8+1,n9)*rat(-n8,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(n8,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-5+n9+n8+n4+2*n3+2*ep,-1+n4)
;
*--#] lastar5_13_x4 : 
#endprocedure
*--#] reduce4 : 
*--#[ reduce9 :
#procedure reduce9()
*--#[ lastar5_13_x9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,n9) gives 6+0=6 terms (0+1+4+0+1) ~ 2.613745 [+,+,+,+,x,+,+,+,-]
* checked poles: (-3+n9+n7+n6+n1+2*ep), [(n9)]
* decreases: (-n9)
* +:                                 n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,             n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-3+n9+n7+n6+n1+2*ep)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-4+n8+n7+n5+n2+3*ep,-3+n9+n7+n6+n1+2*ep)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-5+2*n8+n5+2*n4+n3-n1+3*ep,-3+n9+n7+n6+n1+2*ep)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(-4+n9+n8+n7+n3+n2+n1+2*ep,-3+n9+n7+n6+n1+2*ep)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(-4+n8+n7+n5+n2+3*ep,-3+n9+n7+n6+n1+2*ep)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+2)*rat(1,2)*rat(1+n9,-3+n9+n7+n6+n1+2*ep)
;
*--#] lastar5_13_x9 : 
#endprocedure
*--#] reduce9 : 
*--#[ reduce12368 :
#procedure reduce12368()
*--#[ lastar5_13_x12368 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5+ep,n6,n7,n8,0) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,x,+,+,+,-]
* checked poles: (-2+n8-n6+n3+n2+ep), (-5+n8+n7+n3+n2+n1+2*ep), [(-4+n8+n4+n3+2*ep)]
* decreases: (n1+n2+n3+n6+n8)
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3,         n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?,n6?pos_,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,0)*rat(-4+n7+n6+n1+2*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,0)*rat(-4+n8+n7+n5+n2+3*ep,-2+n8-n6+n3+n2+ep)*rat(-7+2*n8+n7+n5+n4+n3+n2+4*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,0)*rat(-4+n8+n4+n3+2*ep,-2+n8-n6+n3+n2+ep)*rat(-7+2*n8+n7+n5+n4+n3+n2+4*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,0)*rat(4-n7-n6-n1-2*ep,-2+n8-n6+n3+n2+ep)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,0)*rat(4-n8-n4-n3-2*ep,-2+n8-n6+n3+n2+ep)*rat(-4+n8+n7+n5+n2+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
;
*--#] lastar5_13_x12368 : 
#endprocedure
*--#] reduce12368 : 
*--#[ reduce1237 :
#procedure reduce1237()
*--#[ lastar5_13_x1237 : Z(n1,n2,n3,1,n5,n6,n7,n8,0)
* Z(n1,n2,n3,1,n5+ep,n6,n7,n8,0) gives 4+1=5 terms (0+4)+(0+1) ~ 4.011524 [+,+,+,+,x,+,+,+,-]
* checked poles: (-3+n7+n2+n1+ep), (-5+n8+n7+n3+n2+n1+2*ep)
* decreases: (n1+n2+n3+n7)
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3,             n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,1,n5?,n6?pos_,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2,n3,1,n5,n6,n7,n8,0)*rat(-4+n7+n6+n1+2*ep,-3+n7+n2+n1+ep)*rat(-7+n8+2*n7+n6+n5+n2+n1+4*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2-1,n3,1,n5,n6,n7,n8,0)*rat(-4+n8+n7+n5+n2+3*ep,-3+n7+n2+n1+ep)*rat(-7+n8+2*n7+n6+n5+n2+n1+4*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2,n3-1,1,n5,n6,n7,n8,0)*rat(-3+n8+n3+2*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Z'(n1,n2,n3,1,n5,n6,n7-1,n8,0)*rat(4-n7-n6-n1-2*ep,-3+n7+n2+n1+ep)*rat(-4+n8+n7+n5+n2+3*ep,-5+n8+n7+n3+n2+n1+2*ep)
+`Y'(n1,n2,n3,0,n5,n6,n7,n8,0)*rat(3-n8-n3-2*ep,-3+n7+n2+n1+ep)
;
*--#] lastar5_13_x1237 : 
#endprocedure
*--#] reduce1237 : 
*--#[ reduce :
#procedure reduce()
  #call turbo()
* Scheme 2 is simple, but may give spurious poles.
  #if `LASTAR5SCHEME' == 2
    #call reduce9()
    #call reduce12368()
  #elseif `LASTAR5SCHEME' == 3
    #call reduce4()
    #call reduce9()
    #call reduce1237()
  #endif
#endprocedure
*--#] reduce : 
  #call StdReduceProc(turbo,reduce)
#endprocedure
* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4
