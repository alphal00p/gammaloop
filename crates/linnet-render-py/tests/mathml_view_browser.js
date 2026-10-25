/* Test transport only: measurement and mounting use the production widget ESM;
   every final SVG comes from the retained native producer, never a mock. */
import widget from '/view.js';

const assert = (condition, message) => {
  if (!condition) throw new Error(message);
};
const settle = async doc => {
  await doc.fonts.ready;
  await new Promise(resolve => doc.defaultView.requestAnimationFrame(resolve));
  await new Promise(resolve => doc.defaultView.requestAnimationFrame(resolve));
};
const box = element => {
  const r = element.getBoundingClientRect();
  return {left:r.left, top:r.top, right:r.right, bottom:r.bottom,
    width:r.width, height:r.height, centerY:(r.top+r.bottom)/2};
};

export async function inspect(doc = document, request = null) {
  const faces = await doc.fonts.load('9px "Linnet Offline Math"');
  await settle(doc);
  assert(faces.length > 0 && faces.every(face => face.status === 'loaded'),
    'bundled font must actually load');
  assert(doc.fonts.check('9px "Linnet Offline Math"'), 'bundled font check failed');
  const svg = doc.querySelector('svg[data-linnet-interactive]');
  assert(svg && svg.linnetSelection !== undefined, 'native SVG controls did not initialize');
  assert(svg.querySelectorAll('.linnet-viewport').length === 1, 'viewport must initialize once');
  assert(svg.querySelectorAll('.linnet-toolbar').length === 1, 'toolbar must initialize once');
  const labels = [...svg.querySelectorAll('.linnet-math-label')];
  assert(labels.length === 9, `lost native labels: ${labels.length}`);
  const diagnostic = {labels:[],regressions:[]};
  doc.defaultView.linnetGeometryDiagnostic = diagnostic;
  const rects = labels.map((label, index) => {
    const math = label.querySelector('math');
    const fo = label.closest('foreignObject');
    const bounds = box(fo), measured = box(math);
    const authored = fo.width.baseVal.value;
    const scale = bounds.width / authored;
    // Subpixel bounding rectangles may differ by one fractional pixel. This
    // tolerance is in native points and is retained in the capture evidence.
    const epsilon = Math.max(0.15, scale * 0.35);
    const matrix = fo.getScreenCTM();
    const ctm = matrix && Object.fromEntries(['a','b','c','d','e','f'].map(key=>[key,matrix[key]]));
    const expected = request?.labels.find(source => {
      const probe = doc.createElementNS('http://www.w3.org/1998/Math/MathML','math');
      probe.innerHTML = source.mathml;
      return probe.innerHTML === math.innerHTML;
    });
    const evidence = {index,label_id:expected?.id,text:math.textContent,mathml:math.innerHTML,
      bounds,measured,epsilon,ctm,style:label.getAttribute('style'),
      tokens:[...math.querySelectorAll('mi,mn,mo,mtext,ms')].map(token=>({text:token.textContent,box:box(token)}))};
    diagnostic.labels.push(evidence);
    assert(measured.width > 0 && measured.height > 0, `empty MathML label ${index}`);
    const fit = r => r.left >= bounds.left-epsilon && r.right <= bounds.right+epsilon
      && r.top >= bounds.top-epsilon && r.bottom <= bounds.bottom+epsilon;
    assert(fit(measured), `MathML exceeds foreignObject ${index}: ${JSON.stringify({bounds, measured, epsilon})}`);
    const font = doc.defaultView.getComputedStyle(label);
    assert(font.fontFamily.includes('Linnet Offline Math'), `wrong label font: ${font.fontFamily}`);
    assert(font.fontSize === '9px', `wrong local label font size: ${font.fontSize}`);
    assert(font.fontWeight === '400' && font.fontStyle === 'normal', 'host styles changed label font');
    assert(!/overflow\s*:\s*visible/i.test(fo.getAttribute('style') || ''),
      'overflow style masks invalid label geometry');
    const tokens = [...math.querySelectorAll('mi,mn,mo,mtext,ms')].map(token => {
      const rect = box(token);
      assert(rect.width > 0 && rect.height > 0, `empty token ${token.textContent}`);
      assert(fit(rect), `token exceeds foreignObject ${index}: ${token.textContent} ${JSON.stringify({bounds, rect, epsilon})}`);
      return {text:token.textContent, box:rect};
    });
    const scripts = [];
    for (const node of math.querySelectorAll('msub,msup,msubsup,mmultiscripts')) {
      const base = box(node.children[0]);
      let prescript = false, position = 0;
      for (const child of [...node.children].slice(1)) {
        if (child.localName === 'mprescripts') { prescript=true; position=0; continue; }
        const superscript = node.localName === 'msup' ||
          (node.localName !== 'msub' && position % 2 === 1);
        position += 1;
        if (child.localName === 'none') continue;
        const rect = box(child);
        assert(superscript ? rect.centerY < base.centerY : rect.centerY > base.centerY,
          `script on wrong baseline: ${node.outerHTML}`);
        assert(prescript ? rect.right <= base.left+epsilon : rect.left >= base.right-epsilon,
          `script on wrong side: ${node.outerHTML}`);
        scripts.push({superscript,prescript,base,script:rect});
      }
    }
    return {index, label_id:expected?.id, text:math.textContent, bounds, measured, epsilon, tokens, scripts,ctm,
      fraction:math.querySelector('mfrac') !== null,multiscripts:math.querySelector('mmultiscripts') !== null,
      authored:{width:fo.width.baseVal.value,height:fo.height.baseVal.value}, font:font.fontFamily,font_size:font.fontSize};
  });
  if (request) {
    const cases = [
      {name:'fraction-script-minimum-p2-over-m-plus-epsilon',id:4,count:4,
        matches:label=>label.fraction && label.tokens.map(token=>token.text).join('') === 'p2m+ε'
          && label.scripts.some(script=>script.superscript)},
      {name:'tensor-multiscripts',id:1,count:1,
        matches:label=>label.multiscripts && label.scripts.some(script=>script.prescript)},
      {name:'greek-subscript-and-superscript',id:2,count:1,
        matches:label=>label.scripts.some(script=>script.superscript) && label.scripts.some(script=>!script.superscript)},
    ];
    for (const regression of cases) {
      const captured=rects.filter(label=>label.label_id === regression.id);
      const passed=captured.length === regression.count && captured.every(regression.matches);
      const result={name:regression.name,label_id:regression.id,expected_count:regression.count,
        actual_count:captured.length,passed,foreign_objects:captured.map(label=>label.index)};
      diagnostic.regressions.push(result);
      assert(passed,`targeted MathML regression missing: ${JSON.stringify(result)}`);
    }
  }
  for (let i=0;i<rects.length;i++) for (let j=i+1;j<rects.length;j++) {
    const a=rects[i].bounds,b=rects[j].bounds;
    const w=Math.min(a.right,b.right)-Math.max(a.left,b.left);
    const h=Math.min(a.bottom,b.bottom)-Math.max(a.top,b.top);
    assert(w <= 0.35 || h <= 0.35, `native label boxes overlap: ${i},${j} ${JSON.stringify({a,b,w,h})}`);
  }
  const viewport = svg.querySelector('.linnet-viewport');
  const graph = box(viewport);
  assert(graph.width > 0 && graph.height > 0, 'graph viewport is empty');
  for (const rect of rects) {
    assert(rect.bounds.left >= graph.left-0.75 && rect.bounds.right <= graph.right+0.75
      && rect.bounds.top >= graph.top-0.75 && rect.bounds.bottom <= graph.bottom+0.75,
    `label clipped by final viewport: ${rect.index}`);
  }
  const identities = {};
  for (const kind of ['node','edge','halfedge']) identities[kind] = [...new Set(
    [...svg.querySelectorAll(`[data-linnet-kind="${kind}"]`)].map(node=>Number(node.dataset.linnetId)))].sort((a,b)=>a-b);
  assert(JSON.stringify(identities) === JSON.stringify({node:[0,1,2],edge:[0,1,2,3,4],halfedge:[0,1,2,3,4,5,6,7]}),
    `topology changed: ${JSON.stringify(identities)}`);
  return {font_loaded:true,labels:rects,identities,viewport:graph,regressions:diagnostic.regressions,
    width:doc.documentElement.clientWidth,height:doc.documentElement.clientHeight,
    font_embedded:doc.querySelector('style[data-linnet-font]')?.textContent.includes('data:font/woff2;base64,') === true};
}

export async function mountPair(engine, width) {
  const manifest = await (await fetch('/manifest')).json();
  assert(manifest.native === true, 'missing native producer');
  // Deliberately hostile notebook CSS must not affect the measurement iframe.
  const hostile = document.createElement('style');
  hostile.textContent = 'math{font-size:72px!important;font-family:monospace!important}.linnet-math-label{line-height:4!important;letter-spacing:10px!important}';
  document.head.append(hostile);
  const records = [], failures = [];
  for (const slot of ['a','b']) {
    const name = `${engine}-${width}-${slot}`;
    const request = manifest.requests[name];
    assert(request, `native job missing: ${name}`);
    assert(request.coordinate_system === 'svg-user-units' && request.font_style.fontSize === '9px'
      && !('points_per_css_pixel' in request), 'stale measurement coordinate contract');
    const host = document.createElement('section');
    host.className = 'host'; host.style.width = `${width}px`;
    host.dataset.job = name;
    document.querySelector('#hosts').append(host);
    const fields = {request,rendered:{}}, callbacks = new Map();
    const record = {name,request,host,fields,callbacks,messages:[],events:[],receipt:null};
    host.addEventListener('linnet:rendered', event => {
      const doc=event.detail.document;
      record.events.push({id:event.detail.id, controlsReady:doc.querySelector('svg')?.linnetSelection !== undefined});
    });
    record.update = rendered => {
      fields.rendered=rendered;
      for (const callback of callbacks.get('change:rendered') || []) callback();
    };
    record.model = {
      get:name=>fields[name],
      on:(name,callback)=>{ if (!callbacks.has(name)) callbacks.set(name,new Set()); callbacks.get(name).add(callback); },
      off:(name,callback)=>callbacks.get(name)?.delete(callback),
      send:message=>{
        record.messages.push(message);
        if (message.action === 'measure') record.receipt=message.measurements;
      },
    };
    records.push(record);
  }
  await Promise.all(records.map(async record => { record.dispose=await widget.render({model:record.model,el:record.host}); }));
  for (const record of records) {
    assert(record.receipt?.id === record.request.id, 'receipt lost its original identity');
    assert(record.receipt.metrics.length === record.request.labels.length, 'incomplete browser receipt');
    record.update({id:'stale-other-view',html:'<p>stale</p>'});
    assert(!record.host.querySelector('iframe'), 'stale rendered state was mounted');
  }
  const completionOrder=[];
  const waitReady = record => new Promise((resolve,reject) => {
    if (record.events.length) return resolve();
    const timeout=setTimeout(()=>reject(new Error(`no ready event for ${record.name}`)),15000);
    record.host.addEventListener('linnet:rendered',()=>{clearTimeout(timeout);resolve();},{once:true});
  });
  // Both widgets measure first; native completion arrives in reverse order.
  for (const record of [...records].reverse()) {
    const response=await fetch(`/receipt/${record.name}`,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(record.receipt)});
    const result=await response.json();
    assert(response.ok && result.native === true, `native completion failed: ${JSON.stringify(result)}`);
    assert(result.id === record.request.id, 'backend crossed view identity');
    record.result=result;
    assert(result.request.pending === false, "backend did not publish completion");
    record.fields.request={...record.request,...result.request};
    record.update(result);
    await waitReady(record);
    completionOrder.push(record.name);
  }
  const captures=[];
  for (const record of records) {
    assert(record.events.length === 1 && record.events[0].id === record.request.id && record.events[0].controlsReady,
      `ready event preceded controls or was duplicated: ${JSON.stringify(record.events)}`);
    const doc=record.host.querySelector('iframe').contentDocument;
    captures.push({name:record.name,request_id:record.request.id,events:record.events,
      source_labels:record.request.labels,receipt:record.receipt,geometry:await inspect(doc,record.request)});
  }
  const firstDoc=records[0].host.querySelector('iframe').contentDocument;
  const secondDoc=records[1].host.querySelector('iframe').contentDocument;
  const firstSvg=firstDoc.querySelector('svg'),secondSvg=secondDoc.querySelector('svg');
  firstSvg.querySelector('[data-linnet-kind="node"][data-linnet-id="0"][tabindex="0"]')
    .dispatchEvent(new firstDoc.defaultView.MouseEvent('click',{bubbles:true,shiftKey:true}));
  assert(firstSvg.linnetSelection.nodes.includes(0), 'selection controls are inert');
  assert(secondSvg.linnetSelection.nodes.length === 0, 'one view selected another');
  const viewport=firstSvg.querySelector('.linnet-viewport');
  viewport.focus();
  assert(firstDoc.activeElement === viewport,'graph keyboard focus did not reach viewport');
  const before=viewport.getAttribute('viewBox');
  viewport.dispatchEvent(new firstDoc.defaultView.KeyboardEvent('keydown',{key:'+',bubbles:true}));
  assert(viewport.getAttribute('viewBox') !== before,'zoom keyboard control is inert');
  viewport.dispatchEvent(new firstDoc.defaultView.KeyboardEvent('keydown',{key:'0',bubbles:true}));
  await settle(firstDoc);
  assert(viewport.getAttribute('viewBox') === before,'fit control did not restore camera');
  window.fixtureRecords=records;
  window.fixtureFailures=failures;
  return {captures,completion_order:completionOrder,source_sha256:manifest.source_sha256,
    controls:{selection:true,isolated:true,zoom:true,fit:true}};
}

export async function disposeAndRemount() {
  const records=window.fixtureRecords;
  for (const record of records) {
    record.dispose();
    assert([...record.callbacks.values()].every(callbacks=>callbacks.size === 0), 'widget callbacks leaked');
    assert(!record.host.querySelector('iframe'), 'iframe leaked after dispose');
    assert(record.messages.filter(message=>message.action === 'dispose').length === 1, 'dispose transport missing');
  }
  for (const record of records) {
    const before=record.messages.filter(message=>message.action === 'measure').length;
    const ready=new Promise((resolve,reject)=>{
      const timeout=setTimeout(()=>reject(new Error(`remounted view never ready: ${record.name}`)),15000);
      record.host.addEventListener('linnet:rendered',()=>{clearTimeout(timeout);resolve();},{once:true});
    });
    record.dispose=await widget.render({model:record.model,el:record.host});
    await ready;
    assert(record.messages.filter(message=>message.action === 'measure').length === before,
      'completed notebook view silently remeasured');
    await inspect(record.host.querySelector('iframe').contentDocument,record.request);
    record.dispose();
    assert([...record.callbacks.values()].every(callbacks=>callbacks.size === 0), 'remounted callbacks leaked');
  }
  return {callbacks_released:true,frames_removed:true,completed_view_reused:true};
}

export async function mountSavedImage(name) {
  let timer;
  try {
    return await Promise.race([(async () => {
      const url=new URL(`/saved/${name}`,location.href).href;
      const reference=document.createElement('iframe');
      reference.id='image-reference';
      // Both carriers use the identical saved SVG URL with scripts disabled.
      // The document reference exposes real DOM/font boxes for raster evidence.
      reference.setAttribute('sandbox','allow-same-origin');
      Object.assign(reference.style,{position:'absolute',left:'0',top:'0',border:'0',
        width:'2000px',height:'2000px',display:'block'});
      const loaded=new Promise(resolve=>reference.addEventListener('load',resolve,{once:true}));
      reference.src=url;
      document.body.append(reference);
      await loaded;
      const doc=reference.contentDocument,svg=doc?.querySelector('svg');
      assert(svg,'saved image reference did not load the native SVG');
      const faces=await doc.fonts.load('9px "Linnet Offline Math"');
      await doc.fonts.ready;
      // The reference disallows scripts; WebKit suppresses its rAF callbacks.
      // Its loaded fonts still paint on the active parent's rendering frames.
      await settle(document);
      assert(faces.length > 0 && faces.every(face=>face.status === 'loaded'),
        'saved image reference must load the bundled font');
      const bounds=box(svg);
      assert(bounds.width > 0 && bounds.height > 0 && bounds.width < 2000 && bounds.height < 2000,
        `unexpected saved image dimensions: ${JSON.stringify(bounds)}`);
      const image=document.createElement('img');
      image.id='saved-svg-image';
      image.alt='Native graph thumbnail with MathML labels';
      Object.assign(image.style,{position:'absolute',left:'0',top:'0',display:'block',
        width:`${bounds.width}px`,height:`${bounds.height}px`,visibility:'hidden'});
      image.src=url;
      document.body.append(image);
      await image.decode();
      await settle(document);
      assert(image.complete && image.naturalWidth > 0 && image.naturalHeight > 0,
        'saved SVG failed to decode in an img element');
      const labels=[...svg.querySelectorAll('.linnet-math-label')].map((label,index)=>{
        const rectangle=box(label.closest('foreignObject'));
        return {index,text:label.textContent,left:rectangle.left-bounds.left,
          top:rectangle.top-bounds.top,right:rectangle.right-bounds.left,
          bottom:rectangle.bottom-bounds.top,
          tokens:[...label.querySelectorAll('mi,mn,mo,mtext,ms')]
            .filter(token=>token.textContent.trim()).map((token,tokenIndex)=>{
              const tokenBox=box(token);
              return {label_index:index,index:tokenIndex,tag:token.localName,text:token.textContent,
                left:tokenBox.left-bounds.left,top:tokenBox.top-bounds.top,
                right:tokenBox.right-bounds.left,bottom:tokenBox.bottom-bounds.top};
            })};
      });
      assert(labels.length === 9,'saved image reference lost MathML labels');
      const response=await fetch(url);
      assert(response.ok,'negative control could not read the original saved SVG');
      const source=await response.text();
      const sha256=async text=>[...new Uint8Array(await crypto.subtle.digest('SHA-256',
        new TextEncoder().encode(text)))].map(byte=>byte.toString(16).padStart(2,'0')).join('');
      const controlDocument=new DOMParser().parseFromString(source,'image/svg+xml');
      assert(!controlDocument.querySelector('parsererror'),'native SVG is not valid XML');
      // Preserve every non-MathML element/attribute and all font/style/script
      // text. The control changes only the contents of existing math roots.
      const invariant=doc=>JSON.stringify({
        elements:[...doc.querySelectorAll('*')]
          .filter(node=>node.localName === 'math' || !node.closest('math'))
          .map(node=>[node.namespaceURI,node.localName,[...node.attributes].map(attr=>[attr.namespaceURI,attr.name,attr.value])]),
        styles:[...doc.querySelectorAll('style,script')].map(node=>node.textContent),
      });
      const authored=invariant(controlDocument);
      const mathRoots=[...controlDocument.querySelectorAll('.linnet-math-label math')];
      assert(mathRoots.length === labels.length,'negative control must retain all nine MathML roots');
      for (const math of mathRoots) math.replaceChildren();
      assert(invariant(controlDocument) === authored,'negative control changed authored geometry/fonts/styles');
      const negative=document.createElement('img');
      negative.id='saved-svg-negative-control';
      negative.alt='Negative control: same native SVG with only MathML contents emptied';
      negative.style.cssText=image.style.cssText;
      const negativeSource=new XMLSerializer().serializeToString(controlDocument);
      const negativeBlob=new Blob([negativeSource],{type:'image/svg+xml'});
      const negativeUrl=URL.createObjectURL(negativeBlob);
      assert(new URL(negativeUrl).origin === location.origin,'negative control must remain a same-origin blob');
      negative.src=negativeUrl;
      document.body.append(negative);
      await negative.decode();
      await settle(document);
      assert(negative.complete && negative.naturalWidth === image.naturalWidth
        && negative.naturalHeight === image.naturalHeight,'negative control changed intrinsic SVG image dimensions');
      const result={name,url,reference_font_loaded:true,labels,width:bounds.width,height:bounds.height,
        natural_width:image.naturalWidth,natural_height:image.naturalHeight,
        original_svg_sha256:await sha256(source),
        negative_control:{kind:'same saved SVG clone with only MathML contents emptied',
          math_roots_retained:mathRoots.length,geometry_fonts_styles_unchanged:true,
          altered_native_svg_written:false,url:negativeUrl,origin:new URL(negativeUrl).origin,
          byte_length:negativeBlob.size,sha256:await sha256(negativeSource)}};
      window.savedImageFixture={reference,image,negative,negativeUrl,result};
      return result;
    })(),new Promise((_,reject)=>{
      timer=setTimeout(()=>reject(new Error(`saved SVG img test timed out: ${name}`)),15000);
    })]);
  } finally { clearTimeout(timer); }
}

export function showSavedImage(showNegative=false) {
  const {reference,image,negative}=window.savedImageFixture;
  reference.style.visibility='hidden';
  image.style.visibility=showNegative?'hidden':'visible';
  negative.style.visibility=showNegative?'visible':'hidden';
}

export async function compareSavedImagePixels(referencePng,imagePng,negativePng=null) {
  // Screenshot PNGs are origin-clean even where SVG foreignObject deliberately
  // taints a canvas. Never interpret SVG canvas-read restrictions as lost ink.
  const decode=async encoded=>{
    const image=new Image();
    image.src=`data:image/png;base64,${encoded}`;
    await image.decode();
    const canvas=document.createElement('canvas');
    canvas.width=image.naturalWidth;canvas.height=image.naturalHeight;
    const context=canvas.getContext('2d',{willReadFrequently:true});
    context.drawImage(image,0,0);
    return {width:canvas.width,height:canvas.height,
      pixels:context.getImageData(0,0,canvas.width,canvas.height).data};
  };
  const [reference,image,negative]=await Promise.all([
    decode(referencePng),decode(imagePng),negativePng?decode(negativePng):null]);
  assert(reference.width === image.width && reference.height === image.height,
    `saved image screenshot dimensions differ: ${reference.width}x${reference.height} / ${image.width}x${image.height}`);
  const labels=window.savedImageFixture.result.labels.map(label=>{
    let mismatches=0;
    for (let y=Math.max(0,Math.floor(label.top));y<Math.min(reference.height,Math.ceil(label.bottom));y++) {
      for (let x=Math.max(0,Math.floor(label.left));x<Math.min(reference.width,Math.ceil(label.right));x++) {
        const offset=4*(y*reference.width+x);
        if ([0,1,2,3].some(channel=>reference.pixels[offset+channel] !== image.pixels[offset+channel])) mismatches++;
      }
    }
    return {...label,mismatching_pixels:mismatches,passed:mismatches === 0};
  });
  // Retain strict cross-context comparison as a diagnostic. SVG-document and
  // SVG-image raster contexts need not produce byte-identical antialiasing.
  const result={comparison:'exact screenshot pixels within native MathML label boxes (diagnostic)',
    width:reference.width,height:reference.height,labels,passed:labels.every(label=>label.passed)};
  if (negative) {
    assert(negative.width === image.width && negative.height === image.height,
      'same-carrier negative control screenshot dimensions differ');
    const union=new Uint8Array(image.width*image.height);
    const sourceLabels=window.savedImageFixture.result.labels;
    const sourceTokens=sourceLabels.flatMap(label=>label.tokens);
    assert(sourceTokens.length === 39,'saved image reference lost native MathML tokens');
    for (const label of sourceLabels) {
      for (let y=Math.max(0,Math.floor(label.top));y<Math.min(image.height,Math.ceil(label.bottom));y++) {
        for (let x=Math.max(0,Math.floor(label.left));x<Math.min(image.width,Math.ceil(label.right));x++) union[y*image.width+x]=1;
      }
    }
    const darkening=offset=>negative.pixels[offset]+negative.pixels[offset+1]+negative.pixels[offset+2]
      -image.pixels[offset]-image.pixels[offset+1]-image.pixels[offset+2];
    let outsideChanged=0,outsideAlphaChanged=0,controlNoise=0;
    for (let index=0;index<union.length;index++) {
      if (union[index]) continue;
      const offset=4*index;
      if ([0,1,2,3].some(channel=>image.pixels[offset+channel] !== negative.pixels[offset+channel])) outsideChanged++;
      if (image.pixels[offset+3] !== negative.pixels[offset+3]) outsideAlphaChanged++;
      controlNoise=Math.max(controlNoise,Math.abs(darkening(offset)));
    }
    // Emptying MathML can change vector antialiasing by a few RGB levels even
    // with identical authored geometry. Derive the bound from this same pair,
    // then require every label/token to darken more than any control pixel.
    const ink=rectangle=>{
      let changed=0,darkened=0,aboveNoise=0,peak=0;
      for (let y=Math.max(0,Math.floor(rectangle.top));y<Math.min(image.height,Math.ceil(rectangle.bottom));y++) {
        for (let x=Math.max(0,Math.floor(rectangle.left));x<Math.min(image.width,Math.ceil(rectangle.right));x++) {
          const offset=4*(y*image.width+x),signal=darkening(offset);
          if ([0,1,2,3].some(channel=>image.pixels[offset+channel] !== negative.pixels[offset+channel])) changed++;
          if (signal > 0) darkened++;
          if (signal > controlNoise) aboveNoise++;
          peak=Math.max(peak,signal);
        }
      }
      return {...rectangle,changed_pixels:changed,darkened_ink_pixels:darkened,
        ink_pixels_above_control_noise:aboveNoise,peak_rgb_sum_darkening:peak,passed:aboveNoise > 0};
    };
    const paintedLabels=sourceLabels.map(ink),tokens=sourceTokens.map(ink);
    const sourceUnchanged=window.savedImageFixture.result.negative_control.geometry_fonts_styles_unchanged === true;
    result.negative_control={comparison:'same SVG img carrier, MathML contents present versus emptied',
      labels:paintedLabels,tokens,outside_label_changed_pixels:outsideChanged,
      vector_pixels_identical_outside_labels:outsideChanged === 0,
      outside_alpha_changed_pixels:outsideAlphaChanged,
      control_noise_units:'absolute difference of RGB sums (0..765); alpha excluded',
      control_noise_bound_rgb_sum:controlNoise,source_geometry_fonts_styles_unchanged:sourceUnchanged,
      passed:sourceUnchanged && paintedLabels.every(label=>label.passed) && tokens.every(token=>token.passed)};
  }
  return result;
}
