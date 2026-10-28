(() => {
  /* The Marimo adapter emits a hidden container before this document loads.
     Keep its iframe measurable, then reveal only this initialized document.
     Snapshot display never depends on optional font downloads. */
  const frame = window.frameElement;
  if (!frame) return;
  /* A reactive host can navigate the same iframe instead of replacing it.
     This prefix runs before any new drawing markup can paint. */
  frame.style.visibility = 'hidden';
  frame.removeAttribute('data-linnet-ready');
  let shown = false, fallback;
  const reveal = () => {
    if (shown || !frame.isConnected || frame.contentDocument !== document) return;
    shown = true;
    window.clearTimeout(fallback);
    try {
      const style = getComputedStyle(document.body);
      const contents = document.createRange();
      contents.selectNodeContents(document.body);
      const bottom = contents.getBoundingClientRect().bottom + window.scrollY;
      frame.style.height = Math.ceil(bottom + parseFloat(style.paddingBottom) + parseFloat(style.marginBottom)) + 'px';
    } finally {
      /* An initialization or sizing error must not leave a usable static
         snapshot permanently hidden. Exports themselves remain unmodified. */
      frame.style.visibility = 'visible';
      frame.setAttribute('data-linnet-ready', '');
      frame.dispatchEvent(new CustomEvent('linnet:rendered', { bubbles: true }));
    }
  };
  const ready = () => {
    /* Display fallback glyphs while this snapshot's font loads. The separate
       measurement document still requires the exact bundled face. */
    for (const face of document.fonts) {
      if (face.family.replace(/['"]/g, '') === 'Linnet Offline Math') face.display = 'swap';
    }
    /* Inline graph scripts have installed their theme, camera and controls.
       Their queued resize runs before this reveal. Background documents can
       defer animation frames, so retain a finite DOM-ready fallback. */
    requestAnimationFrame(reveal);
    fallback = window.setTimeout(reveal, 100);
  };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', ready, { once: true });
  else ready();
})();
