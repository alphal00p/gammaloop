#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, cetz,edge,hide
#let mom_arr =(stroke:black+0.3mm,marks:((inherit:"head",rev:false,pos:1,scale:40%),))
#let p = plugin("./linnest.wasm")

// step=0.24
//     beta =9.1
//     k_spring=17.6;
//     g_center=0
//     gamma_dangling=50
//     gamma_ee=0.3
//     gamma_ev=0.01
//     length_scale = 0.2
#let layout(input,scope:(:),columns:(1fr),unit:1,additional_data:(:))={
  let a= p.layout_graph(bytes(input),cbor.encode((steps:200,step:1.2,step_shrink:0.21,temp:0.1,beta:9.1,k_spring:20.,g_center:0,gamma_dangling:50,gamma_ee:0.4,gamma_ev:0.1,length_scale:0.2,)+additional_data));


  let graphs= cbor(a);

  let diags = ();
  for (g,parse) in graphs{

  let noed = ();
  let n =(:);

  for (i,v) in g.nodes.enumerate(){
    n.insert(str(i),v)
    let (x,y) = v.remove("pos")
    let b= v.remove("shift")
    let ev = v.remove("eval",default:"(:)")
    if ev == none{
      ev = "(:)"
    }
    noed.push(node(pos:(x,y),name:label(str(i)),..eval(ev,scope: scope ,mode: "code"),layer:2))
    // v.percent=xpercent
  }


  for (i,e) in g.edges.enumerate(){
    let start = e.data.remove("from")
    let end = e.data.remove("to")
    let data = e.remove("data")

    let ev_sink = data.remove("eval_sink",default:"(:)")
    if ev_sink == none{
      ev_sink = "(:)"
    }

    let ev_source = data.remove("eval_source",default:"(:)")
    if ev_source == none{
      ev_source = "(:)"
    }

    let mev = data.remove("mom_eval",default:"(:)")
    if mev == none{
      mev = "(:)"
    }

    let bend-angle = data.remove("bend")

    let bend = bend-angle.remove("Ok",default:0.)

    let o = e.remove("orientation")

    let ev_sink = eval(ev_sink,scope: scope ,mode: "code")
    let ev_source = eval(ev_source,scope: scope ,mode: "code")

    let snmlab = label("sm"+str(i))
    let (start-node,start-node-pos) =if start != none{
      let nodelab = label(str(start))

      if o == "Reversed"{

      noed.push(edge(vertices:((data.pos.x*unit,data.pos.y*unit),nodelab),bend: bend * 0.5rad,..ev_sink))}else{
        noed.push(edge(vertices:(nodelab,(data.pos.x*unit,data.pos.y*unit)),bend: bend * -0.5rad,..ev_source))

      }

      (nodelab,n.at(str(start)).pos)
    } else{
      let lab = label("exts"+str(i))
      noed.push(node((data.pos.x*unit,data.pos.y*unit),name:lab,outset:-5mm,radius:5mm,fill:none))
      (lab,data.pos)
    }

    let enmlab = label("em"+str(i))
    let (end-node,end-node-pos) =if end != none{
      let nodelab = label(str(end))
      if o == "Reversed"{

        let side = if bend < 0{
          "ri"
        }else{
          "right"
        }

      noed.push(edge(vertices:(nodelab,(data.pos.x*unit,data.pos.y*unit)),bend: bend * 0.5rad,..ev_source,label-side:left))}else{
        noed.push(edge(vertices:((data.pos.x*unit,data.pos.y*unit),nodelab),bend: bend * -0.5rad,..ev_sink))

      }

      (nodelab,n.at(str(end)).pos)
    } else{
      let lab = label("exte"+str(i))
      noed.push(node((data.pos.x*unit,data.pos.y*unit),name:lab,outset:-5mm,radius:5mm,fill:none))
      (lab,data.pos)
    }


    let percentb = 1+calc.abs(bend/calc.pi)

    let a  =  calc.sqrt(calc.pow(start-node-pos.x - end-node-pos.x,2)+calc.pow(start-node-pos.y - end-node-pos.y,2))*2.5mm*percentb*unit

    noed.push(node(pos:(start-node-pos.x*unit,start-node-pos.y*unit),name:snmlab,outset:a))
    noed.push(node(pos:(end-node-pos.x*unit,end-node-pos.y*unit),name:enmlab,outset:a))


    let shift = if bend != 0.{
      bend *0.15
    } else{1} * 1.5mm




    // noed.push(edge(vertices:(snmlab,enmlab),bend:bend * (percentb -1) * -1rad,shift:shift,..mom_arr,..eval(mev,scope: scope ,mode: "code")))

  }


 diags.push( diagram(
 node-shape:circle,
	node-fill: black,
 edge-stroke:0.1em,
	spacing: 2em,
  ..noed,
))

// diags.push([#parse ])
  }
  set align(center)

  grid(
    align:center+horizon,gutter:4em,
    columns: columns,
    ..diags
  )
  }
