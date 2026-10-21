"""Exercise exported installed_ltd_display.py output with Playwright.

Usage: python ltd_display_browser.py DIR [browser-executable] [chromium|webkit]
"""

import html
import json
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright


def main(directory, executable=None, engine="chromium"):
    sources = sorted(directory.glob("ltd-*.html"))
    assert sources, "Run installed_ltd_display.py first"
    with sync_playwright() as p:
        browser = getattr(p, engine).launch(
            executable_path=executable,
            args=["--no-sandbox"] if engine == "chromium" else [],
        )
        page = browser.new_page(
            viewport={"width": 800, "height": 1100}, device_scale_factor=1.5
        )
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        for source in sources:
            page.set_content(source.read_text())
            root = page.locator(".feynkit-ltd-result").first
            assert root.get_attribute("data-ready") == "true"
            assert root.locator("details").count() == 2
            assert root.locator("details").evaluate_all(
                "nodes => nodes.every(e => !e.open)"
            )
            assert root.locator(".lp-graph").count() == 0
            assert root.locator(".lp-equation").is_visible()
            root.locator("[data-residue-graph] > summary").click()
            root.locator(".lp-graph").wait_for(state="visible")
            assert root.locator(".lp-path").count() or source.stem == "ltd-contact"
            root.evaluate("""root => {
              const svg=root.querySelector('.lp-graph');
              const native=root.querySelector('[data-drawing]').content;
              const arrows=[...native.querySelectorAll('[data-linnet-momentum]')];
              if(svg.querySelectorAll('[data-momentum-edge]').length!==arrows.length)throw Error('missing momentum arrows');
              for(const original of arrows) {
                const edge=JSON.parse(original.dataset.linnetDetail).edge;
                const rendered=svg.querySelector(`[data-momentum-edge="${edge}"]`);
                const paths=group=>JSON.stringify([...group.querySelectorAll('path')].map(p=>p.getAttribute('d')));
                if(paths(original)!==paths(rendered))throw Error('momentum geometry differs from Linnet');
                if(rendered.querySelectorAll('path').length!==2)throw Error('missing parallel curve or chevron');
              }
              for(const original of native.querySelectorAll('a[data-linnet-kind="edge"][transform]')) {
                const detail=JSON.parse(original.dataset.linnetDetail);
                if('half-edge' in detail)continue;
                const label=[...svg.querySelectorAll('.lp-math-label')].find(label=>
                  ['q','p'].some(symbol=>label.dataset.labelKey===`${symbol}-${detail.edge}`
                    ||label.dataset.labelKey.startsWith(`${symbol}-${detail.edge}-`)));
                if(!label.getAttribute('transform').startsWith(original.getAttribute('transform')))
                  throw Error('label moved from its native placement onto an edge hit target');
              }
              for(const label of svg.querySelectorAll('.lp-math-label')) {
                const use=label.querySelector('use');
                const glyph=document.getElementById(use.getAttribute('href').slice(1));
                if(!svg.contains(glyph)||!glyph.querySelector('use'))throw Error('missing local Typst label');
                if(!use.getBBox().width||!use.getBBox().height)throw Error('empty Typst label');
              }
            }""")
            # Reconstruct every displayed product, including common factors on
            # ancestors. Raised poles exercise signed sums and repeated energies.
            root.evaluate("""root => {
              const g=JSON.parse(root.querySelector('[data-ltd]').textContent);
              for(const r of g.residues) {
                const input=root.querySelector('[data-sum-input]');
                input.value=r.id;
                input.dispatchEvent(new KeyboardEvent('keydown',{key:'Enter',bubbles:true}));
                const expected=r.terms.flatMap(t=>t.chains.map(chain=>({
                  coefficient:t.coefficient,
                  factors:[...t.energies.map(i=>`energy-${i}`),...chain.map(i=>`surface-${i}`)].sort(),
                })));
                const actual=[...root.querySelectorAll('[data-ltd-term]')].map(leaf=>{
                  const factors=[...leaf.querySelectorAll('[data-factor-id]')].flatMap(e=>Array(Number(e.dataset.factorPower || 1)).fill(e.dataset.factorId));
                  for(let node=leaf.parentElement;node&&!node.matches('.lp-equation');node=node.parentElement)
                    if(node.matches('.hs-factor-node'))
                      factors.push(...[...node.querySelectorAll(':scope > .hs-factored-fraction [data-factor-id]')].flatMap(e=>Array(Number(e.dataset.factorPower || 1)).fill(e.dataset.factorId)));
                  const product={coefficient:leaf.querySelector('[data-coefficient]').dataset.coefficient,factors:factors.sort()};
                  if(JSON.stringify(product)!==JSON.stringify(expected[+leaf.dataset.ltdTerm]))throw Error(`altered residue ${r.id}`);
                  return +leaf.dataset.ltdTerm;
                });
                if(new Set(actual).size!==expected.length||actual.length!==expected.length)throw Error(`lost residue ${r.id} terms`);
              }
            }""")
            assert not errors, (source.name, errors)

        results = {}
        for name in ("triangle", "three-loop", "reordered"):
            page.set_content((directory / f"ltd-{name}.html").read_text())
            page.locator("[data-residue-graph] > summary").click()
            page.locator(".lp-graph").wait_for(state="visible")
            results[name] = page.evaluate("""() => {
              const root=document.querySelector('.feynkit-ltd-result');
              const g=JSON.parse(root.querySelector('[data-ltd]').textContent);
              const carriers=[...root.querySelector('[data-drawing]').content.querySelectorAll('[data-linnet-carrier]')].map(e=>JSON.parse(e.dataset.linnetDetail));
              const choose=id=>{const f=root.querySelector('[data-sum-input]');f.value=id;f.dispatchEvent(new KeyboardEvent('keydown',{key:'Enter',bubbles:true}));};
              const region=()=>JSON.stringify([...root.querySelectorAll('[data-linnet-shade-halfedge]')].map(e=>[e.dataset.linnetShadeHalfedge,e.getAttribute('d')]).sort());
              const poles=()=>JSON.stringify([...root.querySelectorAll('[data-pole-role=cut]')].map(e=>[e.dataset.poleEdge,e.dataset.poleSign]));
              const arrows=()=>JSON.stringify([...root.querySelectorAll('[data-momentum-edge]')].map(e=>[e.dataset.momentumEdge,...[...e.querySelectorAll('path')].map(p=>p.getAttribute('d'))]));
              let pairs=0,cycles=0;
              for(const r of g.residues) {
                choose(r.id);
                const fixedPoles=poles(),fixedArrows=arrows();
                const displayed=[...root.querySelectorAll('.lp-equation [data-surface]')].map(e=>+e.dataset.surface).sort((a,b)=>a-b);
                if(JSON.stringify(displayed)!==JSON.stringify([...r.terms[0].chains[0]].sort((a,b)=>a-b)))throw Error('lost native denominator factors');
                for(const [key,pair] of Object.entries(r.propagator_factors)) {
                  const edge=+key;let fixedRegion;
                  for(let branch=0;branch<2;branch++) {
                    root.querySelector(`.lp-equation [data-factor-edge="${edge}"][data-surface="${pair[branch]}"]`).click();
                    if(fixedRegion&&fixedRegion!==region())throw Error('factor selection changed region');
                    fixedRegion=region();
                    const part=new Set([g.routing[edge].head]);let changed=true;
                    while(changed){changed=false;g.routing.forEach((e,i)=>{if(i!==edge&&!r.cuts.includes(i)&&part.has(e.tail)!==part.has(e.head)){part.add(e.tail);part.add(e.head);changed=true;}});}
                    if(part.has(g.routing[edge].tail))throw Error('selected edge did not separate forest');
                    const expected=carriers.filter(p=>part.has(p.node)).map(p=>p["half-edge"]).sort((a,b)=>a-b);
                    const shaded=[...root.querySelectorAll('[data-linnet-shade-halfedge]')].map(p=>+p.dataset.linnetShadeHalfedge).sort((a,b)=>a-b);
                    if(JSON.stringify(expected)!==JSON.stringify(shaded))throw Error('wrong shaded half-edges');
                    const nodes=[...root.querySelectorAll('[data-linnet-shade-node]')].map(p=>+p.dataset.linnetShadeNode).sort((a,b)=>a-b);
                    if(JSON.stringify([...part].sort((a,b)=>a-b))!==JSON.stringify(nodes))throw Error('wrong shaded vertices');
                    const external=Array(g.routing[0].signature.external_signature.length).fill(0);
                    for(const e of Object.values(g.external_routing))e.external_coefficients.forEach((c,i)=>external[i]+=(+part.has(e.source)-+part.has(e.destination))*c);
                    const shifts=new Map(g.surfaces[pair[branch]].expression.external.map(([i,c])=>[i,+c]));
                    if(external.some((c,i)=>c!==(shifts.get(i)||0)))throw Error('external boundary differs from native surface');
                    const coeffs=new Map(g.surfaces[pair[branch]].expression.internal.map(([i,c])=>[g.edges[i],+c]));
                    for(const marker of root.querySelectorAll('[data-pole-edge]')) {
                      const physical=+marker.dataset.poleEdge,dense=g.edges.indexOf(physical),route=g.routing[dense];
                      const sigma=dense===edge?(branch===0?1:-1):(r.signs[dense]==='Reversed'?-1:1);
                      const incidence=+part.has(route.tail)-+part.has(route.head);
                      if(+marker.dataset.poleSign!==sigma)throw Error('wrong native pole sign');
                      if((coeffs.get(physical)||0)!==incidence*sigma)throw Error('region boundary disagrees with native surface');
                      if(marker.dataset.regionDirection!==(incidence?(sigma>0?'out':'in'):'none'))throw Error('pole marker has wrong region direction');
                      if(+marker.dataset.gapWidth!==(r.cuts.includes(dense)?22:8))throw Error('wrong cut gap');
                      const cutout=root.querySelector(`[data-cutout-edge="${physical}"]`);
                      if(cutout.getAttribute('transform')!==marker.getAttribute('transform'))throw Error('aperture detached from marker');
                      const shade=root.querySelector('[data-region-layer]');
                      if(!shade.parentNode.getAttribute('mask'))throw Error('shading bridges a cut gap');
                      const path=root.querySelector(`.lp-path[data-edge="${physical}"][data-flow=source]`),len=path.getTotalLength();
                      const a=path.getPointAtLength(len),b=path.getPointAtLength(Math.max(0,len-1)),m=marker.getCTM();
                      if(incidence&&(m.a*(a.x-b.x)+m.b*(a.y-b.y))*incidence*sigma<0)throw Error('painted gap points the wrong way');
                    }
                    root.querySelector(`[data-surface-edge="${g.edges[edge]}"]`).dispatchEvent(new MouseEvent('click',{bubbles:true}));
                    if(+root.querySelector('.lp-equation [aria-pressed=true]').dataset.surface!==pair[1-branch])throw Error('tree edge did not switch pair');
                    if(fixedRegion!==region()||fixedPoles!==poles()||fixedArrows!==arrows())throw Error('inspection changed fixed routing');
                  }
                  pairs++;
                }
                for(const cut of r.cuts) {
                  const trigger=root.querySelector(`[data-trace-cut="${cut}"]`);
                  trigger.dispatchEvent(new PointerEvent('pointerover',{bubbles:true,pointerType:'mouse'}));
                  const cycle=root.querySelector(`[data-cycle-arrow="${g.edges[cut]}"]`);
                  if(!cycle||+cycle.dataset.cycleDirection!==(r.signs[cut]==='Reversed'?-1:1))throw Error('cycle does not explain native sign');
                  const mark=root.querySelector(`[data-cut-edge="${g.edges[cut]}"]`);
                  mark.dispatchEvent(new MouseEvent('click',{bubbles:true}));
                  const choices=g.residues.filter(x=>x.cuts.includes(cut)),at=choices.findIndex(x=>x.id===r.id);
                  if(+root.querySelector('[data-sum-input]').value!==choices[(at+1)%choices.length].id)throw Error('cut click did not navigate matching residues');
                  choose(r.id);cycles++;
                }
              }
              return {pairs,cycles};
            }""")
            assert not errors, errors

        source = (directory / "ltd-three-loop.html").read_text()
        page.set_content(
            '<body data-theme="dark"><iframe style="width:736px;border:0" srcdoc="'
            + html.escape(source, quote=True)
            + '"></iframe></body>'
        )
        frame = page.frames[1]
        root = frame.locator(".feynkit-ltd-result")
        frame.wait_for_selector('.feynkit-ltd-result[data-theme="dark"]')
        initial = root.bounding_box()["height"]
        root.locator("[data-residue-graph] > summary").click()
        root.locator(".lp-path").first.wait_for()
        page.wait_for_function(
            "parseFloat(document.querySelector('iframe').style.height)>"
            + str(initial + 100)
        )
        field = root.locator("[data-sum-input]")
        field.fill("10")
        field.press("Escape")
        assert field.input_value() == "0"
        for invalid in ("", "-1", "1.5", "56"):
            field.fill(invalid)
            field.press("Enter")
            assert field.input_value() == "0"
            assert root.locator(".hs-status").inner_text()
        field.fill("10")
        root.get_by_role("button", name="Show tree residue 1", exact=True).click()
        assert field.input_value() == "1", "Unsubmitted input swallowed navigation"
        field.fill("10")
        field.press("Enter")
        assert root.locator("[data-sum-id]").evaluate_all(
            "nodes=>nodes.map(n=>+n.dataset.sumId)"
        ) == [0, 9, 10, 11, 55]
        assert root.locator("[data-sum-input]").bounding_box()["width"] < 32
        root.get_by_role("button", name="Show tree residue 55", exact=True).click()
        assert field.input_value() == "55"
        root.get_by_role("button", name="Show tree residue 0", exact=True).click()
        field.fill("1")
        field.press("Enter")
        assert root.locator("[data-residue-graph]").evaluate("e=>e.open")
        root.locator(".lp-equation [data-surface]").first.press("Enter")
        assert root.locator(".lp-equation [data-surface]").first.evaluate(
            "e=>document.activeElement===e"
        )
        camera = root.locator(".linnet-viewport")
        camera_handle = camera.element_handle()
        original = camera.get_attribute("viewBox")
        camera.press("+")
        zoomed = camera.get_attribute("viewBox")
        assert zoomed != original
        assert root.locator(".linnet-toolbar output").inner_text() == "105%"
        root.locator(".lp-equation [data-surface]").first.click()
        assert camera_handle.evaluate("e=>e.isConnected"), (
            "selection replaced Linnet's camera"
        )
        assert camera.get_attribute("viewBox") == zoomed
        marker = root.locator("[data-surface-edge]")
        before_surface = root.locator(".lp-equation [aria-pressed=true]").get_attribute(
            "data-surface"
        )
        marker.press("Enter")
        assert (
            root.locator(".lp-equation [aria-pressed=true]").get_attribute(
                "data-surface"
            )
            != before_surface
        )
        assert marker.evaluate("e=>e===document.activeElement")
        assert camera.get_attribute("viewBox") == zoomed
        # A drag that starts on a cut must pan, without selecting another tree.
        old_id = field.input_value()
        gap = root.locator("[data-cut-edge] rect").first
        point = gap.bounding_box()
        x, y = point["x"] + point["width"] / 2, point["y"] + point["height"] / 2
        page.mouse.move(x, y)
        page.mouse.down()
        page.mouse.move(x + 25, y + 15, steps=5)
        page.mouse.up()
        assert root.locator("[data-sum-input]").input_value() == old_id
        panned = camera.get_attribute("viewBox")
        assert panned != zoomed
        gap = root.locator("[data-cut-edge] rect").first
        old_id = field.input_value()
        gap.hover()
        target = gap.element_handle()
        page.wait_for_timeout(30)
        assert target.evaluate("e=>e.isConnected"), "hover replaced the click target"
        cut = gap.locator("..").get_attribute("data-pole-edge")
        sigma = gap.locator("..").get_attribute("data-pole-sign")
        key = f"pole-{cut}-{'minus' if sigma == '-1' else 'plus'}"
        assert root.locator(f'.lp-math-label[data-label-key="{key}"]').count() == 1
        gap.click()
        assert root.locator("[data-sum-input]").input_value() != old_id
        assert camera.get_attribute("viewBox") == panned
        assert camera_handle.evaluate("e=>e.isConnected")
        camera.press("0")
        assert camera.get_attribute("viewBox") == original
        root.locator("[data-residue-graph] > summary").click()
        root.locator(".lp-graph").wait_for(state="hidden")
        field.fill("10")
        field.press("Enter")
        root.locator(".lp-equation [data-surface]").first.click()
        assert not root.locator("[data-residue-graph]").evaluate("e=>e.open")
        root.locator("[data-residue-graph] > summary").click()
        root.locator(".lp-graph").wait_for(state="visible")
        assert camera_handle.evaluate("e=>e.isConnected")
        assert camera.get_attribute("viewBox") == original
        assert root.locator("[data-cut-edge]").evaluate_all(
            "nodes=>nodes.map(e=>+e.dataset.cutEdge).sort((a,b)=>a-b)"
        ) == root.evaluate("""root => {
            const g=JSON.parse(root.querySelector('[data-ltd]').textContent);
            return g.residues.find(r=>r.id===10).cuts.map(i=>g.edges[i]).sort((a,b)=>a-b);
        }""")
        page.mouse.move(0, 0)
        page.locator("body").evaluate("e=>e.dataset.theme='light'")
        frame.wait_for_selector('.feynkit-ltd-result[data-theme="light"]')
        root.screenshot(path=str(directory / "ltd-three-loop.png"))
        for width in (736, 360, 320):
            page.locator("iframe").evaluate("(e,w)=>e.style.width=w+'px'", width)
            page.wait_for_timeout(50)
            assert root.evaluate("e=>e.scrollWidth<=e.clientWidth+2"), width
        root.screenshot(path=str(directory / "ltd-mobile.png"))

        page.set_content((directory / "ltd-multiple.html").read_text())
        roots = page.locator(".feynkit-ltd-result")
        assert roots.count() == 2
        roots.first.locator('[data-sum-index="1"]').click()
        assert roots.nth(1).locator("[data-sum-input]").input_value() == "0"
        for root in roots.all():
            root.locator("[data-residue-graph] > summary").click()
            root.locator(".lp-graph").wait_for(state="visible")
        assert page.locator(".lp-graph [id]").evaluate_all(
            "nodes=>new Set(nodes.map(n=>n.id)).size===nodes.length"
        ), "Typeset glyph IDs leaked between display instances"
        assert not errors, errors
        print(json.dumps(results))
        browser.close()


if __name__ == "__main__":
    main(Path(sys.argv[1]), *sys.argv[2:])
