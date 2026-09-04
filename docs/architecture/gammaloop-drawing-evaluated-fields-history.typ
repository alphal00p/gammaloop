= GammaLoop Evaluated-Field Proposal

#quote(block: true)[
#strong[Status:] Archived proposal; superseded by the direct callback interface

This label/style precedence and interpolation scheme was proposed but is not the current
interface. It is preserved to explain the migration history. In particular, the current code has
no special `*-eval` fields and performs no `{field}` placeholder interpolation. Use
#link("gammaloop-drawing-architecture.typ")[the current drawing architecture] instead.
]

== Label Fields

/ `label-eval`: Always interpolated and evaluated as Typst content for the edge label. If
  present, it takes precedence over the other label fields.

/ `label`: Plain label text by default. In the GammaLoop embedded template, this is used
  only if `display-label` and `label-template` are absent.

/ `display-label`: Preferred drawing label template.

/ `label-template`: Fallback drawing label template.

Label templates interpolate `{field}` placeholders from the edge callback data:

```dot
a -> b [particle="a", id=7, display-label="{particle} edge {eid}"];
```

Use `{{` and `}}` for literal braces. Unknown placeholders are left unchanged.
The DOT `id` field is consumed as Linnest's stable edge id; draw callbacks expose
that value as `eid`.

== Style Fields

/ `source-style`: Drawing-only style fragment for the source half-edge. In plain mode, string
  values are ignored by the generated GammaLoop callback. In eval mode, the
  string is interpolated and evaluated as Typst, and must evaluate to a
  dictionary.

/ `sink-style`: Drawing-only style fragment for the sink half-edge, with the same behavior as
  `source-style`.

/ `source-style-eval`: Always interpolated and evaluated as a Typst dictionary for the source
  half-edge.

/ `sink-style-eval`: Always interpolated and evaluated as a Typst dictionary for the sink
  half-edge.

== Precedence

Source half-edge style is assembled in this order:

+ the generated model style selected by `particle`;
+ `source-style`, evaluated only when `typst-fields` is `"eval"`;
+ `source-style-eval`, always evaluated.

Sink half-edge style is analogous: generated model style, then `sink-style`,
then `sink-style-eval`. Later dictionaries override earlier keys.

Edge labels use the first available value in this order:

+ `label-eval`, always evaluated;
+ `display-label`;
+ `label-template`;
+ `label`;
+ the generated model label.

The normal GammaLoop-generated path does not rely on `*-eval`. These fields are
an escape hatch for manually edited DOT files.

== Interpolation And Templating

String templates are converted to text, outer quotes are trimmed, and then
`{field}` placeholders are replaced from the edge callback dictionary. Escaped
`{{` and `}}` become literal braces. Unknown placeholders remain in the string,
which makes misspelled fields visible in the rendered output or in Typst eval
errors.

Eval fields are interpolated first and evaluated afterward. Their eval scope is
the generated style scope plus the edge callback dictionary, so expressions can
refer to helpers such as `source-stroke`, `sink-stroke`, `wave`, `coil`,
`zigzag`, and to edge fields such as `particle`, `label`, `eid`,
`source-half-edge`, and `sink-half-edge`.

== Eval Mode

The embedded figure template reads `sys.inputs.typst-fields`. The default is
plain mode:

```bash
linnet draw graphs --input typst-fields=plain
```

Plain mode interpolates label templates but does not execute normal render
fields. This keeps generated GammaLoop drawings data-only.

Eval mode is opt-in:

```bash
linnet draw graphs --input typst-fields=eval
```

In eval mode, the known render fields `label`, `display-label`,
`label-template`, `source-style`, and `sink-style` are interpolated and then
passed to Typst `eval`. Explicit `label-eval`, `source-style-eval`, and
`sink-style-eval` fields are evaluated in both modes because their names are
already an opt-in.

Example:

```dot
digraph styled {
  a -> b [
    particle="a",
    id=7,
    display-label="[#text(fill: red)[{particle} edge {eid}]]",
    source-style="(stroke: red + 1.2pt)",
    sink-style="(stroke: blue + 1.2pt)"
  ];
}
```

This requires `--input typst-fields=eval`.
