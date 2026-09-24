(() => {
  'use strict';
  window.addEventListener('message',event => {
    if (event.source===parent && event.data?.type==='spenso-theme'
        && ['light','dark'].includes(event.data.theme)) {
      document.documentElement.style.colorScheme=event.data.theme;
    }
  });
  const data = JSON.parse(document.getElementById('tensor-data').textContent);
  const shape = data.shape;
  const rank = shape.length;
  const size = shape.reduce((a,b) => a*b,1);
  const entries = data.entries.map(([index, [bytes,plain,html]]) => ({index,bytes,plain,html}));
  const stored = new Map(entries.map(cell => [cell.index.join(','),cell]));
  const axisName = axis => `${axis}: ${data.axes[axis]}`;
  const remaining = () => shape.map((_,i) => i).filter(i => i !== state.row && i !== state.column);
  const state = {
    view:rank === 3 && shape[0] <= 4 ? 'atlas' : 'slice',
    display:'grid', expanded:false,
    row:Math.max(0,rank-2), column:rank > 1 ? rank-1 : -1,
    fixed:Array(rank).fill(0), starts:Array(rank).fill(0),
    selected:entries[0]?.index.slice() || Array(rank).fill(0),
  };
  const allBytes = entries.map(cell => cell.bytes);
  if (data.sparse && size > data.stored) allBytes.push(0);
  const low = allBytes.length ? Math.min(...allBytes) : 0;
  const high = allBytes.length ? Math.max(...allBytes) : 0;
  const payload = entries.reduce((sum,cell) => sum+cell.bytes,0);
  const get = index => stored.get(index.join(',')) || (data.sparse && data.complete ? {
    index, bytes:0, plain:data.default[1], html:data.default[2], implicit:true,
  } : null);
  const el = (tag, text, className) => {
    const node = document.createElement(tag);
    if (text !== undefined) node.textContent = text;
    if (className) node.className = className;
    return node;
  };
  const area = name => document.getElementById(name);
  const label = (text, control) => {
    const node = el('label'); node.append(el('span',text));
    if (control.tagName==='SELECT') {
      const wrapper=el('span',undefined,'select-wrap'); wrapper.append(control); node.append(wrapper);
    } else node.append(control);
    return node;
  };
  function select(name, options, value, change) {
    const control = el('select'); control.dataset.control = name;
    for (const [key,text] of options) { const option=el('option',text); option.value=key; control.append(option); }
    control.value = value;
    control.addEventListener('change',() => {
      change(control.value); render();
      area('controls').querySelector(`[data-control="${name}"]`).focus();
    });
    return control;
  }
  function coordinate(axis, value, change) {
    if (shape[axis]<=16) {
      const control=select(`coordinate-${axis}`,Array.from({length:shape[axis]},(_,i)=>[i,i]),value,next=>change(Number(next)));
      control.dataset.axis=axis;
      return control;
    }
    const input = el('input'); input.type='number'; input.min=0; input.max=Math.max(0,shape[axis]-1); input.step=1;
    input.value=value; input.dataset.axis=axis;
    input.addEventListener('change',() => {
      const parsed = input.valueAsNumber;
      if (!Number.isInteger(parsed) || parsed < 0 || parsed >= shape[axis]) { input.value=value; return; }
      change(parsed); render();
      area('controls').querySelector(`[data-axis="${axis}"]`).focus();
    });
    return input;
  }
  function controls() {
    const host=area('controls'); host.replaceChildren();
    const toolbar=el('div',undefined,'control-bar');
    const fields=el('div',undefined,'control-fields');
    fields.id='control-fields'; fields.hidden=!state.expanded;
    host.append(toolbar);
    const views=[['slice','One slice'],['heaviest','Heaviest first']];
    if (rank>2) views.splice(1,0,['atlas','All slices']);
    const display=el('div',undefined,'display-toggle');
    display.setAttribute('role','group'); display.setAttribute('aria-label','Component display');
    for (const [value,title] of [['grid','Memory grid'],['matrix','Matrix']]) {
      const toggle=el('button',title); toggle.type='button'; toggle.dataset.display=value;
      toggle.setAttribute('aria-pressed',String(state.display===value));
      toggle.addEventListener('click',() => {
        state.display=value;
        if (state.view==='heaviest') state.view='slice';
        render();
        area('controls').querySelector(`[data-display="${value}"]`).focus();
      });
      display.append(toggle);
    }
    const summary=el('output',undefined,'control-summary'); summary.hidden=state.expanded;
    summary.setAttribute('aria-live','polite');
    const description=[views.find(([value])=>value===state.view)[1]];
    if (state.view!=='heaviest') {
      description.push(`Rows ${axisName(state.row)}${state.column<0 ? '' : ' · columns '+axisName(state.column)}`);
      remaining().forEach((axis,i)=>description.push(`${data.axes[axis]} ${state.view==='atlas' && i===0 ? '≥' : '='} ${state.fixed[axis]}`));
      [state.row,state.column].filter(axis=>axis>=0 && state.starts[axis]>0)
        .forEach(axis=>description.push(`Start ${axisName(axis)} = ${state.starts[axis]}`));
    }
    summary.textContent=description.join(' · ');
    const expand=el('button',undefined,'controls-expand'); expand.type='button';
    const action=state.expanded ? 'Collapse controls' : 'Expand controls';
    expand.setAttribute('aria-label',action); expand.title=action;
    expand.setAttribute('aria-expanded',String(state.expanded)); expand.setAttribute('aria-controls',fields.id);
    expand.addEventListener('click',()=> {
      state.expanded=!state.expanded; controls();
      area('controls').querySelector('.controls-expand').focus();
    });
    toolbar.append(display,fields,summary,expand);
    fields.append(label('View',select('view',views,state.view,value => {
      state.view=value;
      if (value==='heaviest') state.display='grid';
    })));
    if (state.view==='heaviest') { fields.classList.add('single'); return; }
    if (rank>1) {
      const axes=shape.map((_,i) => [i,axisName(i)]);
      fields.append(label('Rows',select('row',axes,state.row,value => {
        const old=state.row; state.row=Number(value); if (state.column===state.row) state.column=old;
      })));
      fields.append(label('Columns',select('column',axes,state.column,value => {
        const old=state.column; state.column=Number(value); if (state.row===state.column) state.row=old;
      })));
    }
    const free=remaining();
    free.forEach((axis,i) => fields.append(label(
      state.view==='atlas' && i===0 ? `First slice · ${axisName(axis)}` : `Fix ${axisName(axis)}`,
      coordinate(axis,state.fixed[axis],value => { state.fixed[axis]=value; state.selected[axis]=value; })
    )));
    for (const axis of [state.row,state.column].filter(axis => axis>=0 && shape[axis]>8)) {
      fields.append(label(`Start ${axisName(axis)}`,coordinate(axis,state.starts[axis],value => { state.starts[axis]=value; })));
    }
    fields.style.setProperty('--other-fields',fields.children.length-1);
    fields.classList.toggle('single',fields.children.length===1);
    fields.classList.toggle('many',fields.children.length>4);
  }
  function detail() {
    const host=area('detail'); host.replaceChildren();
    if (size===0) { host.textContent='Empty tensor'; return; }
    const cell=get(state.selected);
    const title=el('div',`Component [${state.selected.join(',')}] · ${cell ? cell.implicit ? 'implicit default · no per-cell payload' : cell.bytes+' bytes' : 'not included in preview'}`);
    const formula=el('div',undefined,'formula');
    if (cell) { formula.innerHTML=cell.html; formula.setAttribute('role','img'); formula.setAttribute('aria-label',cell.plain); }
    host.append(title,formula);
    area('plot').querySelectorAll('[data-coordinate]').forEach(button => button.setAttribute('aria-pressed',String(button.dataset.coordinate===state.selected.join(','))));
  }
  function button(index, matrix=false) {
    const cell=get(index);
    const button=el('button',undefined,matrix ? '' : 'tile'); button.type='button';
    button.dataset.coordinate=index.join(',');
    button.setAttribute('aria-label',`Component [${index.join(',')}], ${cell ? cell.implicit ? 'implicit default' : cell.bytes+' bytes' : 'not included in preview'}`);
    button.setAttribute('aria-pressed',String(index.join(',')===state.selected.join(',')));
    if (matrix && cell) button.innerHTML=cell.html;
    else button.textContent=cell ? index.join(',') : '…';
    if (!cell) button.classList.add('missing');
    if (!matrix && cell) {
      const fraction=high===low ? (high===0 ? 0 : 1) : (cell.bytes-low)/(high-low);
      button.style.backgroundColor=`color-mix(in srgb,var(--blue) ${8+54*fraction}%,var(--bg))`;
    }
    button.addEventListener('click',() => { state.selected=index.slice(); detail(); });
    return button;
  }
  function slice(panel, fixed) {
    const rows=Array.from({length:Math.min(8,shape[state.row]-state.starts[state.row])},(_,i)=>i+state.starts[state.row]);
    const columns=state.column<0 ? [0] : Array.from({length:Math.min(8,shape[state.column]-state.starts[state.column])},(_,i)=>i+state.starts[state.column]);
    const fixedLabel=remaining().map(axis=>`${data.axes[axis]} = ${fixed[axis]}`).join(' · ');
    if (fixedLabel) panel.append(el('div',fixedLabel,'panel-heading'));
    panel.append(el('div',`Rows ${axisName(state.row)}${state.column<0 ? '' : ' · columns '+axisName(state.column)}`,'axes'));
    if (state.display==='matrix') {
      const viewport=el('div',undefined,'matrix-viewport'); const table=el('table',undefined,'matrix');
      table.setAttribute('aria-label',`Matrix slice ${fixedLabel}`);
      const body=el('tbody');
      for (const row of rows) {
        const tr=el('tr');
        for (const column of columns) {
          const index=fixed.slice(); index[state.row]=row; if (state.column>=0) index[state.column]=column;
          const td=el('td'); td.append(button(index,true)); tr.append(td);
        }
        body.append(tr);
      }
      table.append(body); viewport.append(table); panel.append(viewport);
    } else {
      const grid=el('div',undefined,'grid'); grid.style.gridTemplateColumns=`24px repeat(${columns.length},minmax(0,1fr))`;
      grid.append(el('span'));
      columns.forEach(column=>grid.append(el('span',column,'index')));
      for (const row of rows) {
        grid.append(el('span',row,'index'));
        for (const column of columns) {
          const index=fixed.slice(); index[state.row]=row; if (state.column>=0) index[state.column]=column;
          grid.append(button(index));
        }
      }
      panel.append(grid);
    }
    if (rows.length < shape[state.row] || (state.column>=0 && columns.length<shape[state.column])) {
      panel.append(el('div',`Showing rows ${rows[0]}–${rows.at(-1)}${state.column<0?'':' · columns '+columns[0]+'–'+columns.at(-1)}`,'axes'));
    }
  }
  function render() {
    area('shape').textContent=`Rank ${rank} · ${shape.join(' × ')} · ${size.toLocaleString()} components${data.sparse ? ' · '+data.stored.toLocaleString()+' stored' : ''}`;
    controls();
    const plot=area('plot'); plot.replaceChildren();
    if (size===0) { detail(); area('status').textContent='No components'; return; }
    if (state.view==='heaviest') {
      const table=el('table',undefined,'heaviest');
      const head=el('thead'); const row=el('tr'); row.append(el('th','Component'),el('th','Payload bytes')); head.append(row); table.append(head);
      const body=el('tbody');
      [...entries].sort((a,b)=>b.bytes-a.bytes).slice(0,16).forEach(cell => {
        const row=el('tr'); const entry=el('td'); entry.append(button(cell.index)); row.append(entry,el('td',cell.bytes)); body.append(row);
      });
      table.append(body); plot.append(table);
    } else {
      const panels=el('div',undefined,state.view==='atlas' ? 'panels' : 'panels single');
      if (state.display==='matrix') panels.classList.add('matrices');
      const facet=remaining()[0];
      const count=state.view==='atlas' && facet!==undefined ? Math.min(4,shape[facet]-state.fixed[facet]) : 1;
      for (let i=0;i<count;i++) {
        const fixed=state.fixed.slice(); if (state.view==='atlas' && facet!==undefined) fixed[facet]+=i;
        const panel=el('section',undefined,'panel'); slice(panel,fixed); panels.append(panel);
      }
      plot.append(panels);
    }
    const legend=area('legend'); legend.replaceChildren();
    if (state.display==='grid' && state.view!=='heaviest') {
      legend.append(el('span',low+' B'),el('span',undefined,'ramp'),el('span',high+' B'),el('span','Deeper shade = larger payload · shared scale'));
    }
    area('status').textContent=[
      `${payload.toLocaleString()} B ${data.complete ? 'stored component payload' : 'in included components'}`,
      ...(data.sparse ? [`${data.default[0]} B shared default`] : []),
      ...(data.complete ? [] : [`Preview includes ${entries.length} of ${data.stored.toLocaleString()} stored components; … means not loaded`]),
      ...(state.view==='heaviest' ? ['16 heaviest included components'] : []),
      'Excludes allocation overhead and shared symbol metadata',
    ].join(' · ');
    detail();
  }
  render();
})();
