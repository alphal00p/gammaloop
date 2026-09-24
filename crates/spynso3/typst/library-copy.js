(async button => {
  const code=button.closest('.sl-code').querySelector('pre');
  try {
    await navigator.clipboard.writeText(code.textContent);
    button.textContent='Copied';
  } catch {
    const range=document.createRange(); range.selectNodeContents(code);
    const selection=window.getSelection(); selection.removeAllRanges(); selection.addRange(range);
    button.textContent='Selected';
  }
  setTimeout(() => { button.textContent='Copy'; },1800);
})(this);
