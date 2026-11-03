#procedure SchemeZno()
  #ifndef `NOSCHEME'
    #define NOSCHEME "1"
  #endif
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
* signs: +,+,+,+,+,+,+,+,-
*--#] def : 
*--#[ z2y :
#procedure z2y()
  #call Z2Y(9,8)
#endprocedure
*--#] z2y : 
*--#[ turbo :
#procedure turbo()
* depth = 3 (same as 4), threshold = 12
#if `TURBOTHRESHOLD' >= 3
*--#[ no_02_turbo-5c02ccbd91 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n4), [(n9)]
* decreases: (-n9), (n2+n3+n4+n8)
* +:         n3,             n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4,             n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n4)
+`Z'(n1,n2-1,n3+1,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-1+n4)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n4)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n7+n3+2*n2+2*ep,-1+n4)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] no_02_turbo-5c02ccbd91 : 
*--#[ no_02_turbo-8346376252 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 8+0=8 terms (0+0+8) ~ 2.828427 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n6), [(n9)]
* decreases: (-n9), (n1+n2+n6+n7)
* +: n1,                         n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3,         n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,n9+1)*rat(1,2)*rat(n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9+1)*rat(-1,2)*rat(n8,-1+n6)
+`Z'(n1+1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(n1,-1+n6)
+`Z'(n1+1,n2,n3,n4,n5,n6-1,n7-1,n8,n9+1)*rat(-1,2)*rat(n1,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n8+2*n2+n1+2*ep,-1+n6)
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)
;
*--#] no_02_turbo-8346376252 : 
#endif
#if `TURBOTHRESHOLD' >= 4
*--#[ no_02_turbo-86a2e23d78 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n5), [(n9)]
* decreases: (-n9), (n1+n2+n3+n5+n7+n8)
* +:                         n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3,     n5,     n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?{>=2},n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
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
*--#] no_02_turbo-86a2e23d78 : 
*--#[ no_02_turbo-8b509bb5f7 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1), [(n9)]
* decreases: (-n9), (n1+n3+n4+n5+n8)
* +:                     n6, n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3, n4, n5,         n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9+1)*rat(-1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,n9+1)*rat(1,2)*rat(n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6+1,n7,n8,n9+1)*rat(1,2)*rat(n6,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6+1,n7,n8-1,n9+1)*rat(-1,2)*rat(n6,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-5+n7+n6+2*n5+n1+2*ep,-1+n1)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5-1,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)
;
*--#] no_02_turbo-8b509bb5f7 : 
*--#[ no_02_turbo-f13ba82355 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+0+9) ~ 3.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3), [(n9)]
* decreases: (-n9), (n1+n3+n5+n6+n7)
* +:             n4,             n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1,     n3,     n5, n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
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
*--#] no_02_turbo-f13ba82355 : 
*--#[ no_02_turbo-27fdb60717 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 17+0=17 terms (0+0+6+9+2) ~ 3.067627 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n2), [(1+n9), (n9)]
* decreases: (-n9), (n1+n2+n3+n4+n5+n6+n7+n8)
* +: n1,                 n6, n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?{<=-2}) =
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
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(-1,4)*rat(-3+n9+n7+n6+2*n5+n1+2*ep,-1+n2)
+`Z'(n1-1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n2)
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n2)
;
*--#] no_02_turbo-27fdb60717 : 
*--#[ no_02_turbo-4e8aea4872 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 13+0=13 terms (0+0+9+0+4) ~ 3.079446 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n2), [(n9)]
* decreases: (-n9), (n1+n2+n3+n4+n5+n6+n7+n8)
* +:                         n7, n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?{>=2},n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
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
*--#] no_02_turbo-4e8aea4872 : 
*--#[ no_02_turbo-2e112d6d74 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 19+0=19 terms (0+0+11+6+2) ~ 3.624329 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n7), [(1+n9), (n9)]
* decreases: (-n9)
* +: n1, n2,         n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?pos_,n9?{<=-2}) =
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
+`Z'(n1,n2,n3,n4,n5,n6,n7-1,n8,n9+2)*rat(-1,4)*rat(-4+n9+n7+2*n6+n5+n1+2*ep,-1+n7)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7-1,n8,n9+2)*rat(-1,4)*rat(1+n9,-1+n7)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7-1,n8,n9+2)*rat(1,4)*rat(1+n9,-1+n7)
;
*--#] no_02_turbo-2e112d6d74 : 
*--#[ no_02_turbo-035715063a : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 19+0=19 terms (0+0+11+6+2) ~ 3.624329 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n8), [(1+n9), (n9)]
* decreases: (-n9)
* +:     n2, n3,     n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -: n1, n2,     n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>=2},n9?{<=-2}) =
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
*--#] no_02_turbo-035715063a : 
*--#[ no_02_turbo-940d0b3629 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 13+0=13 terms (0+1+9+1+2) ~ 3.645976 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n7), [(n9)]
* decreases: (-n9)
* +:     n2,         n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2, n3, n4, n5, n6, n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?{>=2},n8?pos_,n9?neg_) =
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
*--#] no_02_turbo-940d0b3629 : 
*--#[ no_02_turbo-a796f59059 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 8+0=8 terms (0+2+6) ~ 3.667707 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n8), [(n9)]
* decreases: (-n9)
* +:         n3,     n5,             n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,     n4,         n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?{>=2},n9?neg_) =
+`Z'(n1,n2,n3+1,n4,n5,n6,n7,n8-1,n9+1)*rat(1,2)*rat(n3,-1+n8)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n5,-1+n8)
+`Z'(n1,n2,n3,n4,n5+1,n6,n7-1,n8-1,n9+1)*rat(1,2)*rat(n5,-1+n8)
+`Z'(n1,n2,n3+1,n4-1,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-1+n8)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(-1,2)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(-5+n8+n5+2*n4+n3+2*ep,-1+n8)
;
*--#] no_02_turbo-a796f59059 : 
#endif
#if `TURBOTHRESHOLD' >= 6
*--#[ no_02_turbo-4ac11cabb3 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1), (-1+n6)
* decreases: (n1+n6)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,         n4,     n6
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,n5+1,n6-1,n7+1,n8,n9)*rat(n5,-1+n1)*rat(n7,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6-1,n7+1,n8,n9)*rat(-n7,-1+n1)*rat(-1+n7+n5+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5+1,n6-1,n7,n8,n9)*rat(-n5,-1+n1)*rat(-1+n7+n5+ep,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9)*rat(-3+n6+n1+ep,-1+n6)
;
*--#] no_02_turbo-4ac11cabb3 : 
*--#[ no_02_turbo-901ba44d2b : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 5+0=5 terms (0+5) ~ 5.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3), (-1+n4)
* decreases: (n3+n4)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:         n3, n4,     n6
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2,n3-1,n4-1,n5+1,n6-1,n7,n8+1,n9)*rat(n5,-1+n3)*rat(n8,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5,n6,n7,n8+1,n9)*rat(-n8,-1+n3)*rat(-1+n8+n5+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4-1,n5+1,n6,n7,n8,n9)*rat(-n5,-1+n3)*rat(-1+n8+n5+ep,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-3+n4+n3+ep,-1+n4)
;
*--#] no_02_turbo-901ba44d2b : 
#endif
#if `TURBOTHRESHOLD' >= 7
*--#[ no_02_turbo-a81b261d03 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(n7,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,n9)*rat(-n5,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,n9)*rat(n5,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,n9)*rat(-9+n9+2*n8+n7+n5+n4+2*n3+2*n2+4*ep,-1+n4)
;
*--#] no_02_turbo-a81b261d03 : 
*--#[ no_02_turbo-94275316c3 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n6)
* decreases: (n1+n6)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,n9)*rat(n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,n9)*rat(-n8,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,n9)*rat(-n5,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,n9)*rat(n5,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9)*rat(-9+n9+n8+2*n7+n6+n5+2*n2+2*n1+4*ep,-1+n6)
;
*--#] no_02_turbo-94275316c3 : 
*--#[ no_02_turbo-b2cbcb38bd : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+6+1+2) ~ 6.236841 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1)
* decreases: (n1+n6)
* +:     n2,                 n7,     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6,     n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,n9)*rat(-n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,n9)*rat(n7,-1+n1)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,n9)*rat(-n2,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,n9)*rat(n2,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9)*rat(-9+n9+2*n8+n7+2*n6+2*n5+n2+n1+4*ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,n9)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3-1,n4,n5,n6,n7,n8,n9+1)*rat(-1,2)*rat(n9,-1+n1)
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n9,-1+n1)
;
*--#] no_02_turbo-b2cbcb38bd : 
*--#[ no_02_turbo-b10eadbb5a : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 9+0=9 terms (0+6+1+2) ~ 6.236841 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3)
* decreases: (n3+n4)
* +:     n2,                     n8, n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg0_) =
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
*--#] no_02_turbo-b10eadbb5a : 
#endif
#endprocedure
*--#] turbo : 
*--#[ reduce9 :
#procedure reduce9()
  #call turbo()
*--#[ no_02_9 : Z(n1,n2,n3,n4,n5,n6,n7,n8,n9)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,n9) gives 6+0=6 terms (0+5+0+1) ~ 5.049357 [+,+,+,+,+,+,+,+,-]
* checked poles: (-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep), [(n9)]
* decreases: (-n9)
* +: n1,     n3,                     n9
* 0: n1, n2, n3, n4, n5, n6, n7, n8
* -:     n2,                 n7, n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,n9?neg_) =
+`Z'(n1,n2-1,n3+1,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n3,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
+`Z'(n1,n2,n3+1,n4,n5,n6,n7,n8-1,n9+1)*rat(-1,2)*rat(n3,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
+`Z'(n1+1,n2-1,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(n1,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
+`Z'(n1+1,n2,n3,n4,n5,n6,n7-1,n8,n9+1)*rat(-1,2)*rat(n1,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
+`Z'(n1,n2,n3,n4,n5,n6,n7,n8,n9+1)*rat(1,2)*rat(-3+n9+n3+2*n2+n1+2*ep,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
+`Z'(n1,n2-1,n3,n4,n5,n6,n7,n8,n9+2)*rat(1,2)*rat(1+n9,-7+n9+n8+n7+n6+n5+n4+n3+n2+n1+4*ep)
;
*--#] no_02_9 : 
#endprocedure
*--#] reduce9 : 
*--#[ reduce1346 :
#procedure reduce1346()
*--#[ no_02_1 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,0) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n1)
* decreases: (n1+n6)
* +:     n2,                 n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6,     n8
`id' `Z'(n1?{>=2},n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,0)*rat(-n7,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5-1,n6,n7+1,n8,0)*rat(n7,-1+n1)
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,0)*rat(-n2,-1+n1)
+`Z'(n1-1,n2+1,n3,n4,n5,n6,n7,n8-1,0)*rat(n2,-1+n1)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,0)*rat(-9+2*n8+n7+2*n6+2*n5+n2+n1+4*ep,-1+n1)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,0)
;
*--#] no_02_1 : 
*--#[ no_02_3 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,0) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n3)
* decreases: (n3+n4)
* +:     n2,                     n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1,     n3, n4, n5, n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?{>=2},n4?pos_,n5?pos_,n6?pos_,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2+1,n3-1,n4,n5,n6,n7,n8,0)*rat(-n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5-1,n6,n7,n8+1,0)*rat(n8,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,0)*rat(-n8,-1+n3)
+`Z'(n1,n2+1,n3-1,n4,n5,n6,n7-1,n8,0)*rat(n2,-1+n3)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,0)*rat(-9+n8+2*n7+2*n5+2*n4+n3+n2+4*ep,-1+n3)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,0)
;
*--#] no_02_3 : 
*--#[ no_02_4 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,0) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n4)
* decreases: (n3+n4)
* +:                 n5,     n7
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6,     n8
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?{>=2},n5?pos_,n6?pos_,n7?pos_,n8?pos_,0) =
+`Z'(n1-1,n2,n3,n4-1,n5,n6,n7+1,n8,0)*rat(-n7,-1+n4)
+`Z'(n1,n2-1,n3,n4-1,n5,n6,n7+1,n8,0)*rat(n7,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,0)*rat(-n5,-1+n4)
+`Z'(n1,n2,n3,n4-1,n5+1,n6,n7,n8-1,0)*rat(n5,-1+n4)
+`Z'(n1,n2,n3-1,n4,n5,n6,n7,n8,0)
+`Z'(n1,n2,n3,n4-1,n5,n6,n7,n8,0)*rat(-9+2*n8+n7+n5+n4+2*n3+2*n2+4*ep,-1+n4)
;
*--#] no_02_4 : 
*--#[ no_02_6 : Z(n1,n2,n3,n4,n5,n6,n7,n8,0)
* Z(n1,n2,n3,n4,n5,n6,n7,n8,0) gives 6+0=6 terms (0+6) ~ 6.000000 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n6)
* decreases: (n1+n6)
* +:                 n5,         n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -: n1, n2, n3, n4,     n6, n7
`id' `Z'(n1?pos_,n2?pos_,n3?pos_,n4?pos_,n5?pos_,n6?{>=2},n7?pos_,n8?pos_,0) =
+`Z'(n1,n2-1,n3,n4,n5,n6-1,n7,n8+1,0)*rat(n8,-1+n6)
+`Z'(n1,n2,n3-1,n4,n5,n6-1,n7,n8+1,0)*rat(-n8,-1+n6)
+`Z'(n1,n2,n3,n4-1,n5+1,n6-1,n7,n8,0)*rat(-n5,-1+n6)
+`Z'(n1,n2,n3,n4,n5+1,n6-1,n7-1,n8,0)*rat(n5,-1+n6)
+`Z'(n1-1,n2,n3,n4,n5,n6,n7,n8,0)
+`Z'(n1,n2,n3,n4,n5,n6-1,n7,n8,0)*rat(-9+n8+2*n7+n6+n5+2*n2+2*n1+4*ep,-1+n6)
;
*--#] no_02_6 : 
#endprocedure
*--#] reduce1346 : 
*--#[ reduce2578 :
* Scheme 2 tends to be slightly slower than scheme 1, which may depend on
* the performance for reduced integrals. (TU 2 Nov 2015)
#procedure reduce2578()
#if `NOSCHEME' == 1
*--#[ no_02_2 : Z(1,n2,1,1,n5,1,n7,n8,0)
* Z(1,n2,1,1,n5,1,n7,n8,0) gives 6+5=11 terms (0+6)+(0+5) ~ 6.041733 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n2), (-3+n8+n5+n2+2*ep)
* decreases: (n2+n5+n7+n8)
* +:                         n7, n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,         n5,     n7, n8
`id' `Z'(1,n2?{>=2},1,1,n5?pos_,1,n7?pos_,n8?pos_,0) =
+`Z'(1,n2-1,1,1,n5-1,1,n7,n8+1,0)*rat(-n8,-1+n2)*rat(-5+n8+2*n7+2*n5+4*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5-1,1,n7+1,n8,0)*rat(-n7,-1+n2)*rat(-7+n8+2*n7+2*n5+2*n2+4*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5,1,n7,n8,0)*rat(-1,-3+n8+n5+n2+2*ep)*rat(30-12*n8+n8^2-15*n7+3*n7*n8+2*n7^2-22*n5+4*n5*n8+6*n5*n7+4*n5^2-11*n2+3*n2*n8+2*n2*n7+4*n2*n5+n2^2-32*ep+6*ep*n8+8*ep*n7+12*ep*n5+6*ep*n2+8*ep^2,-1+n2)
+`Z'(1,n2,1,1,n5-1,1,n7,n8,0)*rat(2-n8-n5-2*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2,1,1,n5,1,n7-1,n8,0)
+`Z'(1,n2,1,1,n5,1,n7,n8-1,0)*rat(7-n8-2*n7-2*n5-2*n2-4*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2-1,1,0,n5,1,n7+1,n8,0)*rat(n7,-1+n2)*rat(-7+n8+2*n7+2*n5+2*n2+4*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2-1,1,1,n5,0,n7,n8+1,0)*rat(n8,-1+n2)*rat(-7+n8+2*n7+2*n5+2*n2+4*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2,0,1,n5,0,n7,n8+1,0)*rat(-2*n8,-3+n8+n5+n2+2*ep)
+`Y'(1,n2,0,1,n5,1,n7,n8,0)*rat(-9+2*n8+2*n7+3*n5+2*n2+6*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2,1,1,n5,0,n7,n8,0)*rat(-9+2*n8+2*n7+2*n5+3*n2+6*ep,-3+n8+n5+n2+2*ep)
;
*--#] no_02_2 : 
*--#[ no_02_5 : Z(1,n2,1,1,n5,1,n7,n8,0)
* Z(1,n2,1,1,n5,1,n7,n8,0) gives 6+5=11 terms (0+6)+(0+5) ~ 6.041733 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n5), (-3+n8+n5+n2+2*ep)
* decreases: (n2+n5+n7+n8)
* +:                         n7, n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,         n5,     n7, n8
`id' `Z'(1,n2?pos_,1,1,n5?{>=2},1,n7?pos_,n8?pos_,0) =
+`Z'(1,n2-1,1,1,n5-1,1,n7,n8+1,0)*rat(-n8,-1+n5)*rat(-5+n8+2*n7+2*n2+4*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5-1,1,n7+1,n8,0)*rat(-n7,-1+n5)*rat(-7+n8+2*n7+2*n5+2*n2+4*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5,1,n7,n8,0)*rat(2-n8-n2-2*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2,1,1,n5-1,1,n7,n8,0)*rat(-1,-3+n8+n5+n2+2*ep)*rat(30-12*n8+n8^2-15*n7+3*n7*n8+2*n7^2-11*n5+3*n5*n8+2*n5*n7+n5^2-22*n2+4*n2*n8+6*n2*n7+4*n2*n5+4*n2^2-32*ep+6*ep*n8+8*ep*n7+6*ep*n5+12*ep*n2+8*ep^2,-1+n5)
+`Z'(1,n2,1,1,n5,1,n7-1,n8,0)
+`Z'(1,n2,1,1,n5,1,n7,n8-1,0)*rat(7-n8-2*n7-2*n5-2*n2-4*ep,-3+n8+n5+n2+2*ep)
+`Y'(0,n2,1,1,n5-1,1,n7+1,n8,0)*rat(n7,-1+n5)*rat(-7+n8+2*n7+2*n5+2*n2+4*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2,0,1,n5-1,1,n7,n8+1,0)*rat(n8,-1+n5)*rat(-7+n8+2*n7+2*n5+2*n2+4*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2,0,1,n5,0,n7,n8+1,0)*rat(-2*n8,-3+n8+n5+n2+2*ep)
+`Y'(1,n2,0,1,n5,1,n7,n8,0)*rat(-9+2*n8+2*n7+3*n5+2*n2+6*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2,1,1,n5,0,n7,n8,0)*rat(-9+2*n8+2*n7+2*n5+3*n2+6*ep,-3+n8+n5+n2+2*ep)
;
*--#] no_02_5 : 
*--#[ no_02_7 : Z(1,n2,1,1,n5,1,n7,n8,0)
* Z(1,n2,1,1,n5,1,n7,n8,0) gives 6+5=11 terms (0+6)+(0+5) ~ 6.041733 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n7), (-3+n8+n7+n5+2*ep)
* decreases: (n2+n5+n7+n8)
* +:     n2,         n5
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,         n5,     n7, n8
`id' `Z'(1,n2?pos_,1,1,n5?pos_,1,n7?{>=2},n8?pos_,0) =
+`Z'(1,n2,1,1,n5+1,1,n7-1,n8-1,0)*rat(-n5,-1+n7)*rat(-5+2*n8+n5+2*n2+4*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2+1,1,1,n5,1,n7-1,n8-1,0)*rat(-n2,-1+n7)*rat(-7+2*n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2-1,1,1,n5,1,n7,n8,0)
+`Z'(1,n2,1,1,n5-1,1,n7,n8,0)*rat(7-2*n8-2*n7-n5-2*n2-4*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2,1,1,n5,1,n7-1,n8,0)*rat(-1,-3+n8+n7+n5+2*ep)*rat(30-22*n8+4*n8^2-11*n7+4*n7*n8+n7^2-12*n5+4*n5*n8+3*n5*n7+n5^2-15*n2+6*n2*n8+2*n2*n7+3*n2*n5+2*n2^2-32*ep+12*ep*n8+6*ep*n7+6*ep*n5+8*ep*n2+8*ep^2,-1+n7)
+`Z'(1,n2,1,1,n5,1,n7,n8-1,0)*rat(2-n8-n5-2*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,0,n5+1,0,n7,n8,0)*rat(-2*n5,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,1,n5+1,0,n7-1,n8,0)*rat(n5,-1+n7)*rat(-7+2*n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2+1,0,1,n5,1,n7-1,n8,0)*rat(n2,-1+n7)*rat(-7+2*n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,0,n5,1,n7,n8,0)*rat(-9+3*n8+2*n7+2*n5+2*n2+6*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,1,n5,0,n7,n8,0)*rat(-9+2*n8+3*n7+2*n5+2*n2+6*ep,-3+n8+n7+n5+2*ep)
;
*--#] no_02_7 : 
*--#[ no_02_8 : Z(1,n2,1,1,n5,1,n7,n8,0)
* Z(1,n2,1,1,n5,1,n7,n8,0) gives 6+5=11 terms (0+6)+(0+5) ~ 6.041733 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n8), (-3+n8+n7+n5+2*ep)
* decreases: (n2+n5+n7+n8)
* +:     n2,         n5
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,         n5,     n7, n8
`id' `Z'(1,n2?pos_,1,1,n5?pos_,1,n7?pos_,n8?{>=2},0) =
+`Z'(1,n2,1,1,n5+1,1,n7-1,n8-1,0)*rat(-n5,-1+n8)*rat(-5+2*n7+n5+2*n2+4*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2+1,1,1,n5,1,n7-1,n8-1,0)*rat(-n2,-1+n8)*rat(-7+2*n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2-1,1,1,n5,1,n7,n8,0)
+`Z'(1,n2,1,1,n5-1,1,n7,n8,0)*rat(7-2*n8-2*n7-n5-2*n2-4*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2,1,1,n5,1,n7-1,n8,0)*rat(2-n7-n5-2*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2,1,1,n5,1,n7,n8-1,0)*rat(-1,-3+n8+n7+n5+2*ep)*rat(30-11*n8+n8^2-22*n7+4*n7*n8+4*n7^2-12*n5+3*n5*n8+4*n5*n7+n5^2-15*n2+2*n2*n8+6*n2*n7+3*n2*n5+2*n2^2-32*ep+6*ep*n8+12*ep*n7+6*ep*n5+8*ep*n2+8*ep^2,-1+n8)
+`Y'(0,n2+1,1,1,n5,1,n7,n8-1,0)*rat(n2,-1+n8)*rat(-7+2*n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,0,n5+1,0,n7,n8,0)*rat(-2*n5,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,0,n5+1,1,n7,n8-1,0)*rat(n5,-1+n8)*rat(-7+2*n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,0,n5,1,n7,n8,0)*rat(-9+3*n8+2*n7+2*n5+2*n2+6*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,1,n5,0,n7,n8,0)*rat(-9+2*n8+3*n7+2*n5+2*n2+6*ep,-3+n8+n7+n5+2*ep)
;
*--#] no_02_8 : 
#elseif `NOSCHEME' == 2
*--#[ no_02_2 : Z(1,n2,1,1,n5,1,n7,n8,0)
* Z(1,n2,1,1,n5,1,n7,n8,0) gives 5+6=11 terms (0+5)+(0+6) ~ 5.046025 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n2), (-3+n8+n5+n2+2*ep)
* decreases: (n2+n5+n8)
* +:                         n7, n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,         n5,         n8
`id' `Z'(1,n2?{>=2},1,1,n5?pos_,1,n7?pos_,n8?pos_,0) =
+`Z'(1,n2-1,1,1,n5-1,1,n7,n8+1,0)*rat(-n8,-1+n2)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5-1,1,n7+1,n8,0)*rat(-n7,-1+n2)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5,1,n7,n8,0)*rat(5-n8-n7-2*n5-n2-2*ep,-1+n2)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2,1,1,n5-1,1,n7,n8,0)*rat(2-n8-n5-2*ep,-3+n8+n5+n2+2*ep)
+`Z'(1,n2,1,1,n5,1,n7,n8-1,0)*rat(4-n8-n7-n5-n2-2*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2-1,1,0,n5,1,n7+1,n8,0)*rat(n7,-1+n2)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2-1,1,1,n5,0,n7,n8+1,0)*rat(n8,-1+n2)*rat(-8+2*n8+n7+2*n5+3*n2+4*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2,0,1,n5,0,n7,n8+1,0)*rat(-2*n8,-3+n8+n5+n2+2*ep)
+`Y'(0,n2,1,1,n5,1,n7,n8,0)
+`Y'(1,n2,0,1,n5,1,n7,n8,0)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n5+n2+2*ep)
+`Y'(1,n2,1,1,n5,0,n7,n8,0)*rat(-9+2*n8+2*n7+2*n5+3*n2+6*ep,-3+n8+n5+n2+2*ep)
;
*--#] no_02_2 : 
*--#[ no_02_5 : Z(1,n2,1,1,n5,1,n7,n8,0)
* Z(1,n2,1,1,n5,1,n7,n8,0) gives 5+6=11 terms (0+5)+(0+6) ~ 5.046025 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n5), (-3+n7+n5+n2+2*ep)
* decreases: (n2+n5+n7)
* +:                         n7, n8
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,         n5,     n7
`id' `Z'(1,n2?pos_,1,1,n5?{>=2},1,n7?pos_,n8?pos_,0) =
+`Z'(1,n2-1,1,1,n5-1,1,n7,n8+1,0)*rat(-n8,-1+n5)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n7+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5-1,1,n7+1,n8,0)*rat(-n7,-1+n5)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n7+n5+n2+2*ep)
+`Z'(1,n2-1,1,1,n5,1,n7,n8,0)*rat(2-n7-n2-2*ep,-3+n7+n5+n2+2*ep)
+`Z'(1,n2,1,1,n5-1,1,n7,n8,0)*rat(5-n8-n7-n5-2*n2-2*ep,-1+n5)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n7+n5+n2+2*ep)
+`Z'(1,n2,1,1,n5,1,n7-1,n8,0)*rat(4-n8-n7-n5-n2-2*ep,-3+n7+n5+n2+2*ep)
+`Y'(0,n2,1,0,n5,1,n7+1,n8,0)*rat(-2*n7,-3+n7+n5+n2+2*ep)
+`Y'(0,n2,1,1,n5-1,1,n7+1,n8,0)*rat(n7,-1+n5)*rat(-8+n8+2*n7+3*n5+2*n2+4*ep,-3+n7+n5+n2+2*ep)
+`Y'(1,n2,0,1,n5-1,1,n7,n8+1,0)*rat(n8,-1+n5)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n7+n5+n2+2*ep)
+`Y'(0,n2,1,1,n5,1,n7,n8,0)*rat(-9+2*n8+2*n7+3*n5+2*n2+6*ep,-3+n7+n5+n2+2*ep)
+`Y'(1,n2,1,0,n5,1,n7,n8,0)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n7+n5+n2+2*ep)
+`Y'(1,n2,1,1,n5,0,n7,n8,0)
;
*--#] no_02_5 : 
*--#[ no_02_7 : Z(1,n2,1,1,n5,1,n7,n8,0)
* Z(1,n2,1,1,n5,1,n7,n8,0) gives 5+6=11 terms (0+5)+(0+6) ~ 5.046025 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n7), (-3+n8+n7+n5+2*ep)
* decreases: (n5+n7+n8)
* +:     n2,         n5
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:                 n5,     n7, n8
`id' `Z'(1,n2?pos_,1,1,n5?pos_,1,n7?{>=2},n8?pos_,0) =
+`Z'(1,n2,1,1,n5+1,1,n7-1,n8-1,0)*rat(-n5,-1+n7)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2+1,1,1,n5,1,n7-1,n8-1,0)*rat(-n2,-1+n7)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2,1,1,n5-1,1,n7,n8,0)*rat(4-n8-n7-n5-n2-2*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2,1,1,n5,1,n7-1,n8,0)*rat(5-2*n8-n7-n5-n2-2*ep,-1+n7)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n7+n5+2*ep)
+`Z'(1,n2,1,1,n5,1,n7,n8-1,0)*rat(2-n8-n5-2*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,0,n5+1,0,n7,n8,0)*rat(-2*n5,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,1,n5+1,0,n7-1,n8,0)*rat(n5,-1+n7)*rat(-8+2*n8+3*n7+2*n5+n2+4*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2+1,0,1,n5,1,n7-1,n8,0)*rat(n2,-1+n7)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n7+n5+2*ep)
+`Y'(0,n2,1,1,n5,1,n7,n8,0)
+`Y'(1,n2,1,0,n5,1,n7,n8,0)*rat(-6+2*n8+n7+2*n5+n2+4*ep,-3+n8+n7+n5+2*ep)
+`Y'(1,n2,1,1,n5,0,n7,n8,0)*rat(-9+2*n8+3*n7+2*n5+2*n2+6*ep,-3+n8+n7+n5+2*ep)
;
*--#] no_02_7 : 
*--#[ no_02_8 : Z(1,n2,1,1,n5,1,n7,n8,0)
* Z(1,n2,1,1,n5,1,n7,n8,0) gives 5+6=11 terms (0+5)+(0+6) ~ 5.046025 [+,+,+,+,+,+,+,+,-]
* checked poles: (-1+n8), (-3+n8+n7+n2+2*ep)
* decreases: (n2+n7+n8)
* +:     n2,         n5
* 0: n1, n2, n3, n4, n5, n6, n7, n8, n9
* -:     n2,                 n7, n8
`id' `Z'(1,n2?pos_,1,1,n5?pos_,1,n7?pos_,n8?{>=2},0) =
+`Z'(1,n2,1,1,n5+1,1,n7-1,n8-1,0)*rat(-n5,-1+n8)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2+1,1,1,n5,1,n7-1,n8-1,0)*rat(-n2,-1+n8)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2-1,1,1,n5,1,n7,n8,0)*rat(4-n8-n7-n5-n2-2*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2,1,1,n5,1,n7-1,n8,0)*rat(2-n7-n2-2*ep,-3+n8+n7+n2+2*ep)
+`Z'(1,n2,1,1,n5,1,n7,n8-1,0)*rat(5-n8-2*n7-n5-n2-2*ep,-1+n8)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n2+2*ep)
+`Y'(0,n2+1,0,1,n5,1,n7,n8,0)*rat(-2*n2,-3+n8+n7+n2+2*ep)
+`Y'(0,n2+1,1,1,n5,1,n7,n8-1,0)*rat(n2,-1+n8)*rat(-8+3*n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n2+2*ep)
+`Y'(1,n2,1,0,n5+1,1,n7,n8-1,0)*rat(n5,-1+n8)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n2+2*ep)
+`Y'(0,n2,1,1,n5,1,n7,n8,0)*rat(-9+3*n8+2*n7+2*n5+2*n2+6*ep,-3+n8+n7+n2+2*ep)
+`Y'(1,n2,0,1,n5,1,n7,n8,0)*rat(-6+n8+2*n7+n5+2*n2+4*ep,-3+n8+n7+n2+2*ep)
+`Y'(1,n2,1,1,n5,0,n7,n8,0)
;
*--#] no_02_8 : 
#endif
#endprocedure
*--#] reduce2578 : 
  #call StdReduceProc(turbo,9,1346,2578)
#endprocedure
* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4
