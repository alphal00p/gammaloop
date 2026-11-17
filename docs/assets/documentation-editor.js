import {frameSource} from './code-blocks.js';

// The Markdown boundary owns presentation; worker documentation stays verbatim.
function documentationMarkdown(marked, DOMPurify, python, highlightCode, highlighter) {
  const renderer = new marked.Renderer();
  const sectionBoundary = /\n(?=(?:Parameters|Other Parameters|Returns|Yields|Raises|Warns|Warnings|Notes|See Also|References|Examples|Attributes|Methods)[ \t]*\n[-=]{3,}[ \t]*(?:\n|$))/;
  const cell = (source, language, output = false) => {
    const code = document.createElement('code');
    if (!output && /^(python|py|pycon|python-console)$/.test(language)) {
      highlightCode(source, python.parser.parse(source), highlighter, (text, classes) => {
        const span = document.createElement('span');
        span.className = classes;
        span.textContent = text;
        code.append(span);
      }, () => code.append('\n'));
    } else code.textContent = source;
    const pre = document.createElement('pre');
    pre.append(code);
    if (output) {
      pre.className = 'documentation-example-output';
      pre.setAttribute('aria-label', 'Expected output');
      return pre.outerHTML;
    }
    const container = document.createElement('div');
    container.append(pre);
    frameSource(pre, {compact: true});
    return container.innerHTML;
  };
  renderer.code = ({text, lang = ''}) => {
    const language = lang.split(/\s+/)[0].toLowerCase();
    if (!/^(|python|py|pycon|python-console)$/.test(language) || !/^ *>>> /.test(text)) {
      return cell(text, language);
    }
    // ty can include the following NumPy section inside its example fence.
    // Recover that structural boundary without rewriting the library docstring.
    const section = sectionBoundary.exec(text);
    const suffix = section ? renderMarkdown(text.slice(section.index + 1)) : '';
    if (section) text = text.slice(0, section.index);
    const indent = text.match(/^( *)>>> /)[1];
    let html = '', source = [], output = [], command = [], continuing = false;
    const incomplete = () => {
      let error = false;
      python.parser.parse(command.join('\n')).iterate({enter: node => { if (node.type.isError) error = true; }});
      return error;
    };
    const flush = () => {
      if (source.length) html += cell(source.join('\n'), 'python');
      if (output.length) html += cell(output.join('\n'), '', true);
      source = []; output = [];
    };
    for (const line of text.split('\n')) {
      const local = line.startsWith(indent) ? line.slice(indent.length) : line;
      if (/^>>>(?: |$)/.test(local)) {
        if (output.length) flush();
        command = [local.slice(4)]; continuing = false;
        source.push(local.slice(4));
      } else if (!output.length && source.length && (continuing || incomplete()) && /^(\.\.\. |\.\.\.$)/.test(local)) {
        command.push(local.slice(4));
        continuing = local.length > 4 || incomplete();
        if (!continuing) command = [];
        source.push(local.slice(4));
      } else output.push(local);
    }
    flush();
    return html + suffix;
  };
  // Block tokenizers also run inside lists. Blank lines end unfenced runs;
  // explicit continuation prompts retain blank input lines. Fences stay owned
  // by Marked, so other languages and literal prompts are never rewritten.
  const transcript = source => {
    if (!/^ {0,3}>>> /.test(source)) return;
    const section = sectionBoundary.exec(source);
    const body = section ? source.slice(0, section.index + 1) : source;
    const raw = body.match(/^[^\n]*(?:\n(?![ \t]*\n)[^\n]+)*(?:\n|$)/)[0];
    return {type: 'code', raw, text: raw.replace(/\n$/, ''), lang: 'pycon'};
  };
  const renderMarkdown = source => {
    // These built-in constructors only need the module alias used in their
    // real class example. Keep other context-dependent examples untouched.
    const modelSetup = /Using the setup in the `+Model`+ class example:[ \t]*(?:\r?\n)*/;
    const builtInModel = /^([ \t]*)>>> model = hep\.Model\.(?:phi3|phi4|phi_3_4|standard_model|scalar_qed)\(\)/m;
    if (modelSetup.test(source) && builtInModel.test(source)) {
      source = source.replace(modelSetup, '').replace(builtInModel,
        (line, indent) => `${indent}>>> from symbolica.community import hepkit as hep\n${line}`);
    }
    const tokens = marked.lexer(source, {...marked.defaults, extensions: {
      block: [transcript],
      startBlock: [source => { const match = /(?:^|\n) {0,3}>>> /.exec(source); return match ? match.index : undefined; }],
    }});
    // NumPy parameter sections are plain text in completion responses. Retain
    // their term/description structure instead of collapsing it into a paragraph.
    for (let index = 1; index < tokens.length; ++index) {
      const token = tokens[index], heading = tokens[index - 1];
      if (heading.type !== 'heading' || !/^(Parameters|Other Parameters|Returns|Yields|Attributes)\s*$/.test(heading.text) ||
          token.type !== 'paragraph' || !/^\S[^\n]*\n[ \t]+\S/.test(token.text)) continue;
      let html = '', term = '', description = [];
      const flush = () => {
        if (term) html += `<dt>${marked.parseInline(term)}</dt><dd>${marked.parseInline(description.join('\n'))}</dd>`;
      };
      for (const line of token.text.split('\n')) {
        if (/^[ \t]/.test(line)) description.push(line.trim());
        else { flush(); term = line.trim(); description = []; }
      }
      flush();
      token.type = 'html'; token.text = `<dl class="documentation-parameters">${html}</dl>`;
    }
    return DOMPurify.sanitize(marked.parser(tokens, {renderer}));
  };
  return renderMarkdown;
}

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
    if (profile === 'full') {
      const [state, completion, python] = await Promise.all([
        import(new URL('dist-BczecT6V.js', root).href),
        import(new URL('dist-Ccyblc8N.js', root).href),
        import(new URL('dist-DueXZKz9.js', root).href),
      ]);
      const renderMarkdown = documentationMarkdown(marked, DOMPurify, python.a, state.q, state.K);
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
    const renderMarkdown = documentationMarkdown(marked, DOMPurify, python.a, state.B, state.z);
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
