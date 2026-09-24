// Runs in the notebook page when an explorer frame loads. WebKit does not
// consistently propagate an iframe's color-scheme to its media queries.
(frame => {
  const key=Symbol.for('spenso.explorer.theme');
  if (!window[key]) {
    const media=matchMedia('(prefers-color-scheme: dark)');
    const sent=new WeakMap();
    const sync=loaded => {
      document.querySelectorAll('iframe[data-spenso-explorer]').forEach(frame => {
        const schemes=getComputedStyle(frame).colorScheme.split(/\s+/);
        const theme=schemes.includes('dark') && !schemes.includes('light') ? 'dark'
          : schemes.includes('light') && !schemes.includes('dark') ? 'light'
          : media.matches ? 'dark' : 'light';
        if (loaded===frame || sent.get(frame)!==theme) {
          frame.contentWindow.postMessage({type:'spenso-theme',theme},'*');
          sent.set(frame,theme);
        }
      });
    };
    let scheduled=false;
    const schedule=() => {
      if (scheduled) return;
      scheduled=true;
      requestAnimationFrame(() => { scheduled=false; sync(); });
    };
    // One observer per page, with no strong references to removed frames.
    new MutationObserver(schedule).observe(document.documentElement,{
      subtree:true, attributes:true,
      attributeFilter:['class','style','data-theme','data-jp-theme-light'],
    });
    media.addEventListener('change',schedule);
    window[key]=sync;
  }
  window[key](frame);
})(this);
