#import "draw-raw-dot.typ": *
#set page(width: auto, height: auto)
#show: draw_dot

```dot
digraph BNL {
  num = "1";
  overall_factor = "1";
  overall_factor_evaluated = "1";
  projector = "vbar(0,spenso::bis(4,hedge(0)))
    * u(1,spenso::bis(4,hedge(1)))
    * v(2,spenso::bis(4,hedge(2)))
    * ubar(3,spenso::bis(4,hedge(3)))
    * ebar(4,spenso::mink(4,hedge(4)))
    * ebar(5,spenso::mink(4,hedge(5)))
    * ebar(6,spenso::mink(4,hedge(6)))";

  edge [label=""]
  ext [style=invis]
  ext  -> 0:0 [dir=back id=0 particle="e+" label="$e^+(p_#edge)$" pos="x:@-left!,y:0" "label-anchor"="east"];
  ext -> 0:1 [id=1 particle="e-" label="$e^-(p_#edge)$" pos="x:@-left!,y:1" "label-anchor"="east"];

  1:2 -> ext [ dir=back id=2 particle="mu+" label="$mu^+(p_#edge)$" pos="x:@+right!,y:20" "label-anchor"="west"];
  1:3 -> ext [id=3 particle="mu-" label="$mu^-(p_#edge)$" pos="x:@+right!,y:10" "label-anchor"="west"];
  2:4 -> ext [id=4 particle="a" label="$gamma(p_#edge)$" pos="x:@+right!,y:0" "label-anchor"="west"];
  3:5 -> ext [id=5 particle="a" label="$gamma(p_#edge)$" pos="x:@+right!,y:-10" "label-anchor"="west"];
  4:6 -> ext [id=6 particle="a" label="$gamma(p_#edge)$" pos="x:@+right!,y:-20" "label-anchor"="west"];

  0 -> 5 [id=7 particle="a"  label="$gamma^star$"];
  6 -> 1 [id=8 particle="a" label="$gamma^star$"];

  5 -> 6 [id=9 particle="t" lmb_id=0];
  6 -> 2 [id=10 particle="t"];
  2 -> 3 [id=11 particle="t" ];
  3 -> 7 [id=12 particle="t"];
  7 -> 4 [id=13 particle="t"];
  4 -> 8 [id=14 particle="t"];
  8 -> 5 [id=15 particle="t" label="$t$"];
  7 -> 8 [id=16 particle="g" lmb_id=1];
}

```
