#import "../src/lib.typ": draw, graph, layout, subgraph

#set page(width: auto, height: auto)
#set text(size: 5pt)

#let tree-label(g) = {
  let value = graph.info(g).global-statements.at("tree", default: none)
  if value == none {
    panic("tree.typ: expected a graph-level `tree` statement")
  }
  str(value).trim("\"")
}

#let half-edge-count(g) = {
  let count = 0
  for edge in graph.edges(g) {
    for endpoint in (edge.source, edge.sink) {
      if endpoint != none {
        count = calc.max(count, endpoint.hedge + 1)
      }
    }
  }
  count
}

#let complement-subgraph(g, selected) = {
  let selected-hedges = subgraph.hedges(selected)
  subgraph.bits(
    g,
    range(half-edge-count(g)).map(i => not selected-hedges.contains(i)),
  )
}

#let tree-depths(g, tree, root: 0) = {
  let tree-hedges = subgraph.hedges(tree)
  let adjacency = (:)
  for edge in graph.edges(g) {
    if edge.source != none and edge.sink != none {
      let in-tree = tree-hedges.contains(edge.source.hedge) or tree-hedges.contains(edge.sink.hedge)
      if in-tree {
        let source = str(edge.source.node)
        let sink = str(edge.sink.node)
        adjacency.insert(source, adjacency.at(source, default: ()) + (edge.sink.node,))
        adjacency.insert(sink, adjacency.at(sink, default: ()) + (edge.source.node,))
      }
    }
  }

  let depths = (:)
  let frontier = (root,)
  depths.insert(str(root), 0)
  while frontier.len() > 0 {
      let next = ()
      for node in frontier {
        for child in adjacency.at(str(node), default: ()) {
        if not (str(child) in depths) {
          depths.insert(str(child), depths.at(str(node)) + 1)
          next.push(child)
        }
      }
    }
    frontier = next
  }
  depths
}

#let tree-edge(edge, tree-hedges) = {
  let source = edge.at("source-half-edge")
  let sink = edge.at("sink-half-edge")
  let source-in-tree = source != none and tree-hedges.contains(source.hedge)
  let sink-in-tree = sink != none and tree-hedges.contains(sink.hedge)
  source-in-tree or sink-in-tree
}

#let edge-depth(edge, depths) = {
  let depth = none
  for half-edge in (edge.at("source-half-edge"), edge.at("sink-half-edge")) {
    if half-edge != none {
      let value = depths.at(str(half-edge.node), default: 0)
      if depth == none or value < depth {
        depth = value
      }
    }
  }
  if depth == none { 0 } else { depth }
}

#let node-label(node) = {
  let value = node.at("label", default: node.name)
  if value == none {
    none
  } else {
    let value = str(value).trim("\"")
    if value == "" {
      value = str(node.name)
    }
    text(weight: "bold")[#value]
  }
}

#let edge-style(depths, tree-hedges) = edge => {
  if tree-edge(edge, tree-hedges) {
    let depth = edge-depth(edge, depths)
    (stroke: rgb("#555555") + 1.7pt / (1 + depth / 3))
  } else {
    (stroke: rgb("#b9bec7") + 0.35pt)
  }
}

#let draw_tree_dot(doc) = {
  show raw: it => if it.at("lang") == "dot" {
    for g in graph.parse(it.text) {
      let tree = subgraph.label(g, tree-label(g))
      let rest = complement-subgraph(g, tree)
      let tree-hedges = subgraph.hedges(tree)
      let g = layout(
        g,
        layout-algo: "tree",
        subgraph: tree,
        layout-roots: (0,),
        tree-dx: 26.0,
        tree-dy: 12.0,
        label-steps: 0,
      )
      let g = layout(
        g,
        layout-algo: "dot",
        subgraph: rest,
        layout-nodes: "fixed",
        tree-dx: 26.0,
        tree-dy: 12.0,
        label-steps: 0,
      )
      let depths = tree-depths(g, tree)
      draw(
        g,
        subgraph: tree,
        unit: 1.0,
        edge-stroke: 0.35pt,
        edge-omega: 1.18,
        source-style: edge-style(depths, tree-hedges),
        sink-style: edge-style(depths, tree-hedges),
        node-fill: none,
        node-stroke: none,
        node-min-radius: 0.01,
        node-label-padding: 0.0,
        node-label: node-label,
        node-label-style: (
          padding: 0.08,
          frame: "rect",
          fill: white,
          stroke: 0.55pt + rgb("#777777"),
        ),
        node-outset: 0.28,
        subgraph-edge-style: (stroke: rgb("#72b7b2") + 1.2pt),
      )
    }
  }
  doc
}

#show: draw_tree_dot

```dot
digraph {
  node	 [shape=circle,height=0.1];
  overlap = "scale";
  layout = "neato";
  tree = "3UDua8U1UeS1uIlguU5mUxjPAMQX5nq1VXmmhhorLyBLdGBJYSqx8ebV2lOp"
  0	 [label = "∏"];
  1	 [label = "S:scalar(0)"];
  2	 [label = "∑"];
  3	 [label = "∏"];
  4	 [label = "S:scalar(1)"];
  5	 [label = "∑"];
  6	 [label = "∏"];
  7	 [label = "S:scalar(2)"];
  8	 [label= "L:gamma()"];
  9	 [label = "∏"];
  10	 [label = "S:scalar(3)"];
  11	 [label = "∑"];
  12	 [label = "∏"];
  13	 [label = "S:4:-1"];
  14	 [label = "T:Q3(16)"];
  15	 [label = "∏"];
  16	 [label = "S:scalar(5)"];
  17	 [label = "T:δ(cind(0))"];
  18	 [label = "∑"];
  19	 [label = "∏"];
  20	 [label = "S:6:-1"];
  21	 [label = "T:Q3(16)"];
  22	 [label = "∏"];
  23	 [label = "S:scalar(7)"];
  24	 [label = "T:δ(cind(0))"];
  25	 [label= "L:gamma()"];
  26	 [label= "L:gamma()"];
  27	 [label= "L:gamma()"];
  28	 [label= "L:gamma()"];
  29	 [label= "L:gamma()"];
  30	 [label = "T:ϵbar(6)"];
  31	 [label = "∏"];
  32	 [label = "S:scalar(8)"];
  33	 [label = "∑"];
  34	 [label = "∏"];
  35	 [label = "S:9:MT"];
  36	 [label= "L:g()"];
  37	 [label = "∏"];
  38	 [label = "S:scalar(10)"];
  39	 [label = "T:()"];
  40	 [label = "T:()"];
  41	 [label = "∑"];
  42	 [label = "∏"];
  43	 [label = "S:11:MT"];
  44	 [label= "L:g()"];
  45	 [label = "∏"];
  46	 [label = "S:scalar(12)"];
  47	 [label = "T:()"];
  48	 [label = "T:()"];
  49	 [label = "T:()"];
  50	 [label= "L:gamma()"];
  51	 [label= "L:gamma()"];
  52	 [label = "∑"];
  53	 [label = "∏"];
  54	 [label = "S:13:MT"];
  55	 [label= "L:g()"];
  56	 [label = "∏"];
  57	 [label = "S:scalar(14)"];
  58	 [label = "T:()"];
  59	 [label = "T:()"];
  60	 [label = "∑"];
  61	 [label = "∏"];
  62	 [label = "S:15:MT"];
  63	 [label= "L:g()"];
  64	 [label = "∏"];
  65	 [label = "S:scalar(16)"];
  66	 [label = "T:()"];
  67	 [label = "T:()"];
  68	 [label = "∑"];
  69	 [label = "∏"];
  70	 [label = "S:17:MT"];
  71	 [label= "L:g()"];
  72	 [label = "∏"];
  73	 [label = "S:scalar(18)"];
  74	 [label = "T:()"];
  75	 [label = "T:()"];
  76	 [label = "∑"];
  77	 [label = "∏"];
  78	 [label = "S:19:MT"];
  79	 [label= "L:g()"];
  80	 [label = "∏"];
  81	 [label = "S:scalar(20)"];
  82	 [label = "T:()"];
  83	 [label = "T:()"];
  84	 [label = "∑"];
  85	 [label = "∏"];
  86	 [label = "S:21:MT"];
  87	 [label= "L:g()"];
  88	 [label = "∏"];
  89	 [label = "S:scalar(22)"];
  90	 [label = "T:()"];
  91	 [label = "T:()"];
  92	 [label = "T:()"];
  93	 [label = "T:()"];
  94	 [label = "T:()"];
  95	 [label = "T:()"];
  96	 [label= "L:gamma()"];
  97	 [label= "L:gamma()"];
  ext0	 [style=invis];
  0:0:s	-> ext0	 [id=0 color="red"];
  97:353:s	-> 0:1:s	 [id=1  color="red:blue;0.5"];
  96:349:s	-> 0:2:s	 [id=2  color="red:blue;0.5"];
  95:346:s	-> 0:3:s	 [id=3  color="red:blue;0.5"];
  94:343:s	-> 0:4:s	 [id=4  color="red:blue;0.5"];
  93:341:s	-> 0:5:s	 [id=5  color="red:blue;0.5"];
  92:339:s	-> 0:6:s	 [id=6  color="red:blue;0.5"];
  84:310:s	-> 0:7:s	 [id=7  color="red:blue;0.5"];
  76:281:s	-> 0:8:s	 [id=8  color="red:blue;0.5"];
  68:252:s	-> 0:9:s	 [id=9  color="red:blue;0.5"];
  60:223:s	-> 0:10:s	 [id=10  color="red:blue;0.5"];
  52:194:s	-> 0:11:s	 [id=11  color="red:blue;0.5"];
  2:15:s	-> 0:12:s	 [id=12  color="red:blue;0.5"];
  1:14:s	-> 0:13:s	 [id=13  color="red:blue;0.5"];
  44:172:s	-> 41:165:s	 [id=14 dir=none  color="red:blue;0.5" label="bis4|h22-0"];
  3:24:s	-> 2:16:s	 [id=15  color="red:blue;0.5"];
  31:117:s	-> 2:17:s	 [id=16  color="red:blue;0.5"];
  2:18:s	-> 60:231:s	 [id=17 dir=none  color="red:blue;0.5" label="bis4|h12-0"];
  50:187:s	-> 2:19:s	 [id=18 dir=none  color="red:blue;0.5" label="bis4|h12-0"];
  5:32:s	-> 2:20:s	 [id=19 dir=none  color="red:blue;0.5" label="bis4|h12-0"];
  2:21:s	-> 84:314:s	 [id=20 dir=none  color="red:blue;0.5" label="bis4|h25-0"];
  51:192:s	-> 2:22:s	 [id=21 dir=none  color="red:blue;0.5" label="bis4|h25-0"];
  5:35:s	-> 2:23:s	 [id=22 dir=none  color="red:blue;0.5" label="bis4|h25-0"];
  46:176:s	-> 45:175:s	 [id=23  color="red:blue;0.5"];
  30:115:s	-> 3:25:s	 [id=24  color="red:blue;0.5"];
  5:29:s	-> 3:26:s	 [id=25  color="red:blue;0.5"];
  4:28:s	-> 3:27:s	 [id=26  color="red:blue;0.5"];
  20:86:s	-> 19:85:s	 [id=27  color="red:blue;0.5"];
  6:41:s	-> 5:30:s	 [id=28  color="red:blue;0.5"];
  9:49:s	-> 5:31:s	 [id=29  color="red:blue;0.5"];
  44:170:s	-> 42:167:s	 [id=30  color="red:blue;0.5"];
  26:100:s	-> 5:33:s	 [id=31 dir=none  color="red:blue;0.5" label="bis4|h12-0"];
  8:46:s	-> 5:34:s	 [id=32 dir=none  color="red:blue;0.5" label="bis4|h12-0"];
  50:189:s	-> 51:193:s	 [id=33 dir=none  color="red:blue;0.5" label="mink4|h24-0"];
  27:105:s	-> 5:36:s	 [id=34 dir=none  color="red:blue;0.5" label="bis4|h25-0"];
  8:47:s	-> 5:37:s	 [id=35 dir=none  color="red:blue;0.5" label="bis4|h25-0"];
  5:38:s	-> 30:116:s	 [id=36 dir=none  color="red:blue;0.5" label="mink4|h6-0"];
  25:98:s	-> 5:39:s	 [id=37 dir=none  color="red:blue;0.5" label="mink4|h6-0"];
  8:48:s	-> 5:40:s	 [id=38 dir=none  color="red:blue;0.5" label="mink4|h6-0"];
  27:104:s	-> 28:109:s	 [id=39 dir=none  color="red:blue;0.5" label="bis4|h14-0"];
  8:45:s	-> 6:42:s	 [id=40  color="red:blue;0.5"];
  7:44:s	-> 6:43:s	 [id=41  color="red:blue;0.5"];
  25:96:s	-> 28:108:s	 [id=42 dir=none  color="red:blue;0.5" label="bis4|h13-0"];
  25:97:s	-> 29:113:s	 [id=43 dir=none  color="red:blue;0.5" label="bis4|h22-0"];
  26:101:s	-> 29:112:s	 [id=44 dir=none  color="red:blue;0.5" label="bis4|h21-0"];
  26:102:s	-> 27:106:s	 [id=45 dir=none  color="red:blue;0.5" label="mink4|h24-0"];
  29:111:s	-> 9:50:s	 [id=46  color="red:blue;0.5"];
  28:107:s	-> 9:51:s	 [id=47  color="red:blue;0.5"];
  27:103:s	-> 9:52:s	 [id=48  color="red:blue;0.5"];
  26:99:s	-> 9:53:s	 [id=49  color="red:blue;0.5"];
  25:95:s	-> 9:54:s	 [id=50  color="red:blue;0.5"];
  18:77:s	-> 9:55:s	 [id=51  color="red:blue;0.5"];
  11:59:s	-> 9:56:s	 [id=52  color="red:blue;0.5"];
  10:58:s	-> 9:57:s	 [id=53  color="red:blue;0.5"];
  16:74:s	-> 15:73:s	 [id=54  color="red:blue;0.5"];
  12:65:s	-> 11:60:s	 [id=55  color="red:blue;0.5"];
  15:71:s	-> 11:61:s	 [id=56  color="red:blue;0.5"];
  11:62:s	-> 29:114:s	 [id=57 dir=none  color="red:blue;0.5" label="mink4|e13-1"];
  17:76:s	-> 11:63:s	 [id=58 dir=none  color="red:blue;0.5" label="mink4|e13-1"];
  14:70:s	-> 11:64:s	 [id=59 dir=none  color="red:blue;0.5" label="mink4|e13-1"];
  17:75:s	-> 15:72:s	 [id=60  color="red:blue;0.5"];
  14:69:s	-> 12:66:s	 [id=61  color="red:blue;0.5"];
  13:68:s	-> 12:67:s	 [id=62  color="red:blue;0.5"];
  23:92:s	-> 22:91:s	 [id=63  color="red:blue;0.5"];
  19:83:s	-> 18:78:s	 [id=64  color="red:blue;0.5"];
  22:89:s	-> 18:79:s	 [id=65  color="red:blue;0.5"];
  18:80:s	-> 28:110:s	 [id=66 dir=none  color="red:blue;0.5" label="mink4|e14-1"];
  24:94:s	-> 18:81:s	 [id=67 dir=none  color="red:blue;0.5" label="mink4|e14-1"];
  21:88:s	-> 18:82:s	 [id=68 dir=none  color="red:blue;0.5" label="mink4|e14-1"];
  24:93:s	-> 22:90:s	 [id=69  color="red:blue;0.5"];
  21:87:s	-> 19:84:s	 [id=70  color="red:blue;0.5"];
  43:169:s	-> 42:168:s	 [id=71  color="red:blue;0.5"];
  51:190:s	-> 31:118:s	 [id=72  color="red:blue;0.5"];
  50:186:s	-> 31:119:s	 [id=73  color="red:blue;0.5"];
  49:183:s	-> 31:120:s	 [id=74  color="red:blue;0.5"];
  41:154:s	-> 31:121:s	 [id=75  color="red:blue;0.5"];
  33:125:s	-> 31:122:s	 [id=76  color="red:blue;0.5"];
  32:124:s	-> 31:123:s	 [id=77  color="red:blue;0.5"];
  39:148:s	-> 37:145:s	 [id=78  color="red:blue;0.5"];
  34:137:s	-> 33:126:s	 [id=79  color="red:blue;0.5"];
  37:144:s	-> 33:127:s	 [id=80  color="red:blue;0.5"];
  40:151:s	-> 33:128:s	 [id=81  color="red:blue;0.5"];
  33:129:s	-> 49:184:s	 [id=82 dir=none  color="red:blue;0.5" label="bis4|h13-0"];
  40:152:s	-> 33:130:s	 [id=83 dir=none  color="red:blue;0.5" label="bis4|h13-0"];
  39:149:s	-> 33:131:s	 [id=84 dir=none  color="red:blue;0.5" label="bis4|h13-0"];
  36:142:s	-> 33:132:s	 [id=85 dir=none  color="red:blue;0.5" label="bis4|h13-0"];
  33:133:s	-> 51:191:s	 [id=86 dir=none  color="red:blue;0.5" label="bis4|h14-0"];
  40:153:s	-> 33:134:s	 [id=87 dir=none  color="red:blue;0.5" label="bis4|h14-0"];
  39:150:s	-> 33:135:s	 [id=88 dir=none  color="red:blue;0.5" label="bis4|h14-0"];
  36:143:s	-> 33:136:s	 [id=89 dir=none  color="red:blue;0.5" label="bis4|h14-0"];
  38:147:s	-> 37:146:s	 [id=90  color="red:blue;0.5"];
  36:141:s	-> 34:138:s	 [id=91  color="red:blue;0.5"];
  35:140:s	-> 34:139:s	 [id=92  color="red:blue;0.5"];
  47:177:s	-> 45:174:s	 [id=93  color="red:blue;0.5"];
  42:166:s	-> 41:155:s	 [id=94  color="red:blue;0.5"];
  45:173:s	-> 41:156:s	 [id=95  color="red:blue;0.5"];
  48:180:s	-> 41:157:s	 [id=96  color="red:blue;0.5"];
  41:158:s	-> 50:188:s	 [id=97 dir=none  color="red:blue;0.5" label="bis4|h21-0"];
  48:181:s	-> 41:159:s	 [id=98 dir=none  color="red:blue;0.5" label="bis4|h21-0"];
  47:178:s	-> 41:160:s	 [id=99 dir=none  color="red:blue;0.5" label="bis4|h21-0"];
  44:171:s	-> 41:161:s	 [id=100 dir=none  color="red:blue;0.5" label="bis4|h21-0"];
  41:162:s	-> 49:185:s	 [id=101 dir=none  color="red:blue;0.5" label="bis4|h22-0"];
  48:182:s	-> 41:163:s	 [id=102 dir=none  color="red:blue;0.5" label="bis4|h22-0"];
  47:179:s	-> 41:164:s	 [id=103 dir=none  color="red:blue;0.5" label="bis4|h22-0"];
  58:217:s	-> 56:214:s	 [id=104  color="red:blue;0.5"];
  53:206:s	-> 52:195:s	 [id=105  color="red:blue;0.5"];
  56:213:s	-> 52:196:s	 [id=106  color="red:blue;0.5"];
  59:220:s	-> 52:197:s	 [id=107  color="red:blue;0.5"];
  52:198:s	-> 94:344:s	 [id=108 dir=none  color="red:blue;0.5" label="bis4|h9-0"];
  59:221:s	-> 52:199:s	 [id=109 dir=none  color="red:blue;0.5" label="bis4|h9-0"];
  58:218:s	-> 52:200:s	 [id=110 dir=none  color="red:blue;0.5" label="bis4|h9-0"];
  55:211:s	-> 52:201:s	 [id=111 dir=none  color="red:blue;0.5" label="bis4|h9-0"];
  52:202:s	-> 95:347:s	 [id=112 dir=none  color="red:blue;0.5" label="bis4|h10-0"];
  59:222:s	-> 52:203:s	 [id=113 dir=none  color="red:blue;0.5" label="bis4|h10-0"];
  58:219:s	-> 52:204:s	 [id=114 dir=none  color="red:blue;0.5" label="bis4|h10-0"];
  55:212:s	-> 52:205:s	 [id=115 dir=none  color="red:blue;0.5" label="bis4|h10-0"];
  57:216:s	-> 56:215:s	 [id=116  color="red:blue;0.5"];
  55:210:s	-> 53:207:s	 [id=117  color="red:blue;0.5"];
  54:209:s	-> 53:208:s	 [id=118  color="red:blue;0.5"];
  62:238:s	-> 61:237:s	 [id=119  color="red:blue;0.5"];
  61:235:s	-> 60:224:s	 [id=120  color="red:blue;0.5"];
  64:242:s	-> 60:225:s	 [id=121  color="red:blue;0.5"];
  67:249:s	-> 60:226:s	 [id=122  color="red:blue;0.5"];
  60:227:s	-> 95:348:s	 [id=123 dir=none  color="red:blue;0.5" label="bis4|h11-0"];
  67:250:s	-> 60:228:s	 [id=124 dir=none  color="red:blue;0.5" label="bis4|h11-0"];
  66:247:s	-> 60:229:s	 [id=125 dir=none  color="red:blue;0.5" label="bis4|h11-0"];
  63:240:s	-> 60:230:s	 [id=126 dir=none  color="red:blue;0.5" label="bis4|h11-0"];
  66:246:s	-> 64:243:s	 [id=127  color="red:blue;0.5"];
  67:251:s	-> 60:232:s	 [id=128 dir=none  color="red:blue;0.5" label="bis4|h12-0"];
  66:248:s	-> 60:233:s	 [id=129 dir=none  color="red:blue;0.5" label="bis4|h12-0"];
  63:241:s	-> 60:234:s	 [id=130 dir=none  color="red:blue;0.5" label="bis4|h12-0"];
  65:245:s	-> 64:244:s	 [id=131  color="red:blue;0.5"];
  63:239:s	-> 61:236:s	 [id=132  color="red:blue;0.5"];
  74:275:s	-> 72:272:s	 [id=133  color="red:blue;0.5"];
  69:264:s	-> 68:253:s	 [id=134  color="red:blue;0.5"];
  72:271:s	-> 68:254:s	 [id=135  color="red:blue;0.5"];
  75:278:s	-> 68:255:s	 [id=136  color="red:blue;0.5"];
  68:256:s	-> 96:350:s	 [id=137 dir=none  color="red:blue;0.5" label="bis4|h15-0"];
  75:279:s	-> 68:257:s	 [id=138 dir=none  color="red:blue;0.5" label="bis4|h15-0"];
  74:276:s	-> 68:258:s	 [id=139 dir=none  color="red:blue;0.5" label="bis4|h15-0"];
  71:269:s	-> 68:259:s	 [id=140 dir=none  color="red:blue;0.5" label="bis4|h15-0"];
  68:260:s	-> 97:354:s	 [id=141 dir=none  color="red:blue;0.5" label="bis4|h16-0"];
  75:280:s	-> 68:261:s	 [id=142 dir=none  color="red:blue;0.5" label="bis4|h16-0"];
  74:277:s	-> 68:262:s	 [id=143 dir=none  color="red:blue;0.5" label="bis4|h16-0"];
  71:270:s	-> 68:263:s	 [id=144 dir=none  color="red:blue;0.5" label="bis4|h16-0"];
  73:274:s	-> 72:273:s	 [id=145  color="red:blue;0.5"];
  71:268:s	-> 69:265:s	 [id=146  color="red:blue;0.5"];
  70:267:s	-> 69:266:s	 [id=147  color="red:blue;0.5"];
  82:304:s	-> 80:301:s	 [id=148  color="red:blue;0.5"];
  77:293:s	-> 76:282:s	 [id=149  color="red:blue;0.5"];
  80:300:s	-> 76:283:s	 [id=150  color="red:blue;0.5"];
  83:307:s	-> 76:284:s	 [id=151  color="red:blue;0.5"];
  76:285:s	-> 97:355:s	 [id=152 dir=none  color="red:blue;0.5" label="bis4|h19-0"];
  83:308:s	-> 76:286:s	 [id=153 dir=none  color="red:blue;0.5" label="bis4|h19-0"];
  82:305:s	-> 76:287:s	 [id=154 dir=none  color="red:blue;0.5" label="bis4|h19-0"];
  79:298:s	-> 76:288:s	 [id=155 dir=none  color="red:blue;0.5" label="bis4|h19-0"];
  76:289:s	-> 94:345:s	 [id=156 dir=none  color="red:blue;0.5" label="bis4|h20-0"];
  83:309:s	-> 76:290:s	 [id=157 dir=none  color="red:blue;0.5" label="bis4|h20-0"];
  82:306:s	-> 76:291:s	 [id=158 dir=none  color="red:blue;0.5" label="bis4|h20-0"];
  79:299:s	-> 76:292:s	 [id=159 dir=none  color="red:blue;0.5" label="bis4|h20-0"];
  81:303:s	-> 80:302:s	 [id=160  color="red:blue;0.5"];
  79:297:s	-> 77:294:s	 [id=161  color="red:blue;0.5"];
  78:296:s	-> 77:295:s	 [id=162  color="red:blue;0.5"];
  86:325:s	-> 85:324:s	 [id=163  color="red:blue;0.5"];
  85:322:s	-> 84:311:s	 [id=164  color="red:blue;0.5"];
  88:329:s	-> 84:312:s	 [id=165  color="red:blue;0.5"];
  91:336:s	-> 84:313:s	 [id=166  color="red:blue;0.5"];
  90:333:s	-> 88:330:s	 [id=167  color="red:blue;0.5"];
  91:337:s	-> 84:315:s	 [id=168 dir=none  color="red:blue;0.5" label="bis4|h25-0"];
  90:334:s	-> 84:316:s	 [id=169 dir=none  color="red:blue;0.5" label="bis4|h25-0"];
  87:327:s	-> 84:317:s	 [id=170 dir=none  color="red:blue;0.5" label="bis4|h25-0"];
  84:318:s	-> 96:351:s	 [id=171 dir=none  color="red:blue;0.5" label="bis4|h26-0"];
  91:338:s	-> 84:319:s	 [id=172 dir=none  color="red:blue;0.5" label="bis4|h26-0"];
  90:335:s	-> 84:320:s	 [id=173 dir=none  color="red:blue;0.5" label="bis4|h26-0"];
  87:328:s	-> 84:321:s	 [id=174 dir=none  color="red:blue;0.5" label="bis4|h26-0"];
  89:332:s	-> 88:331:s	 [id=175  color="red:blue;0.5"];
  87:326:s	-> 85:323:s	 [id=176  color="red:blue;0.5"];
  92:340:s	-> 97:356:s	 [id=177 dir=none  color="red:blue;0.5" label="mink4|h18-0"];
  93:342:s	-> 96:352:s	 [id=178 dir=none  color="red:blue;0.5" label="mink4|h8-0"];
}
```
