#import "../theme.typ": palette

#set page(width: auto, height: auto, margin: (x: 1mm, y: 1mm), fill: none)
#set text(fill: palette.ink, size: 22pt)

$
  frac(dif sigma, dif cal(O))
  = sum_(Gamma in upright("graphs"))
    integral product_(i = 1)^(n_"loop") dif^3 bold(k)_i
    sum_(c in upright("cuts")) I_(Gamma,c)^upright("LU")
    delta(cal(O)(c, bold(k)_i))
$
