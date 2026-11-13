function outsetGutter(ViewPlugin) {
  return ViewPlugin.fromClass(class {
      constructor(view) { this.measure(view); }
      update(update) { this.measure(update.view); }
      measure(view) {
        view.requestMeasure({key: this,
          read: () => view.dom.querySelector('.cm-gutters')?.getBoundingClientRect().width || 0,
          write: width => {
            const value = `${width}px`;
            if (view.dom.style.getPropertyValue('--notebook-gutter-width') !== value) {
              view.dom.style.setProperty('--notebook-gutter-width', value);
              view.requestMeasure();
            }
          },
        });
      }
    });
}
// One constructor-bound source editor and outset-gutter owner for direct cells and Marimo.
const editorAPIs = new Map();
export async function load(runtime, profile = 'islands') {
  const root = new URL('./', runtime).href;
  const key = profile + root;
  if (!editorAPIs.has(key)) editorAPIs.set(key, (async () => {
    const [markdown, {t: DOMPurify}] = await Promise.all([
      import(new URL(profile === 'full' ? 'marked.esm-B_3iMIUq.js' : 'toast-DvY70eFX.js', root).href),
      import(new URL(profile === 'full' ? 'purify.es-CYRiAeUx.js' : 'purify.es-DgKyQviY.js', root).href),
    ]);
    const marked = profile === 'full' ? markdown.t : markdown.c;
    const renderMarkdown = source => DOMPurify.sanitize(marked.parse(source));
    if (profile === 'full') {
      const [state, completion] = await Promise.all([
        import(new URL('dist-BczecT6V.js', root).href),
        import(new URL('dist-Ccyblc8N.js', root).href),
      ]);
      return {renderMarkdown, outsetGutter: outsetGutter(state.ut), EditorView: state.st, EditorState: state.Wt,
        StateEffect: state.Zt, StateField: state.Qt, ViewPlugin: state.ut,
        hoverTooltip: state.wt, syntaxHighlighting: state.V, classHighlighter: state.K,
        intelligenceField: state.Qt.define({create: () => true, update: () => true}),
        autocompletion: completion.r};
    }
    const [state, completion, python, history] = await Promise.all([
      import(new URL('dist-Kxgjxakw.js', root).href),
      import(new URL('dist-B1vYLK28.js', root).href),
      import(new URL('dist-C2KzV-vf.js', root).href),
      import(new URL('dist-CO1AwwG_.js', root).href),
    ]);
    // Bind extensions to the actual runtime owning this view. Identical bundle
    // bytes at different content URLs still create different JS constructors.
    return {renderMarkdown, outsetGutter: outsetGutter(state.tt), ViewPlugin: state.tt, StateField: state.Rt, syntaxHighlighting: state.F, classHighlighter: state.z,
      highlightActiveLineGutter: state.dt, EditorState: state.At, EditorView: state.$, StateEffect: state.Lt,
      intelligenceField: state.Rt.define({create: () => true, update: () => true}),
      lineNumbers: state.ht, highlightActiveLine: state.ut, keymap: state.mt,
      hoverTooltip: state.pt, autocompletion: completion.r,
      python: python.i, highlightStyle: state.z,
      StreamLanguage: state.f,
      history: history.o, historyKeymap: history.c, completionKeymap: completion.c};
  })());
  return editorAPIs.get(key);
}
// The pinned Marimo bundle supplies Python and the shared stream parser API.
// Typst uses that same CodeMirror instance, including nested comment/raw state.
export function typstLanguage(cm) {
  return cm.StreamLanguage.define({
    name: 'typst',
    startState: () => ({comment: 0, raw: '', quote: false}),
    token(stream, state) {
      if (state.raw) {
        while (!stream.eol()) {
          if (stream.match(state.raw)) { state.raw = ''; break; }
          stream.next();
        }
        return 'string';
      }
      if (!state.comment && (state.quote || stream.peek() === '"')) {
        if (!state.quote) stream.next();
        state.quote = true;
        while (!stream.eol()) {
          const character = stream.next();
          if (character === '\\') stream.next();
          else if (character === '"') { state.quote = false; break; }
        }
        return 'string';
      }
      if (!state.comment && stream.match('/*')) state.comment = 1;
      if (state.comment) {
        while (!stream.eol()) {
          if (stream.match('/*')) ++state.comment;
          else if (stream.match('*/')) { if (!--state.comment) break; }
          else stream.next();
        }
        return 'comment';
      }
      if (stream.eatSpace()) return null;
      if (stream.match('//')) { stream.skipToEnd(); return 'comment'; }
      const raw = stream.match(/`+/);
      if (raw) { state.raw = raw[0]; return 'string'; }
      if (stream.match(/#?(?:let|set|show|import|include|as|if|else|for|in|while|return|break|continue|true|false|none|auto)\b/)) return 'keyword';
      if (stream.match(/(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?(?:pt|mm|cm|in|em|deg|rad|%|fr)?/i)) return 'number';
      if (stream.match(/<[^<>\n]+>|@[\w-]+/)) return 'labelName';
      if (stream.match(/#?[\p{L}_][\p{L}\d_-]*/u)) return /^\s*\(/.test(stream.string.slice(stream.pos)) ? 'function' : 'variableName';
      return /[+*\-\/=<>!&|:]/.test(stream.next()) ? 'operator' : null;
    },
  });
}
