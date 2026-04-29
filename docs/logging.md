# Logging

This file summarizes the default tracing configuration used by the CLI/API at startup.

## Default directives

Normal session defaults come from `GlobalSettings`:

- `display_directive = "info"`
- `logfile_directive = "off"`

These are the application defaults written into `global_settings.toml` when defaults are shown or persisted.

In practice, `info` is the current convention for normal screen-facing output. Selective internal diagnostics should generally stay on `debug` and be turned on via narrower target/tag filters rather than by raising the whole display logger.

## Startup precedence

Startup applies logging in this order:

1. Base settings come from `GlobalSettings`.
2. Environment overrides replace the base spec:
   - `GL_ALL_LOG_FILTER` overrides both stderr and logfile filters.
   - `GL_DISPLAY_FILTER` overrides only the stderr/display filter.
   - `GL_LOGFILE_FILTER` overrides only the logfile filter.
3. CLI overrides are then applied:
   - `-l/--level` overrides the stderr/display filter for the session.
   - `-L/--logfile-level` overrides the logfile filter for the session.
4. Some CLI modes hard-disable logfile logging for the session:
   - `--read-only-state`
   - `--logfile-level off`

When logfile logging is hard-disabled at boot, later settings changes cannot re-enable it for that session.

## CLI level mapping

`-l/--level` and `-L/--logfile-level` expand to explicit crate directives:

- `off` -> `gammaloop_api=off,gammalooprs=off`
- `error` -> `gammaloop_api=error,gammalooprs=error`
- `warn` -> `gammaloop_api=warn,gammalooprs=warn`
- `info` -> `gammaloop_api=info,gammalooprs=info`
- `debug` -> `gammaloop_api=debug,gammalooprs=debug`
- `trace` -> `gammaloop_api=trace,gammalooprs=trace`

## Empty-spec fallback floors

The filter builders also have hardcoded fallback floors used when a spec is empty:

- display/stderr builder fallback: `off`
- logfile builder fallback: `warn`

This is separate from the normal application defaults above. In other words:

- if startup uses the normal settings path, the defaults are `info` for display and `off` for logfile
- if an empty filter string is explicitly parsed, the display builder falls back to `off` and the logfile builder falls back to `warn`

## Runtime updates

Changing `global.display_directive` or `global.logfile_directive` after startup reloads the active tracing filters. The CLI stderr override from `-l/--level` is treated as a session override and is separate from the persisted global setting.

## Tag-based debug filtering

For selective debug logging, GammaLoop uses its own logging DSL rather than raw `EnvFilter`.

This follows the general idea in Tom Mrazik's tag-based logging note:
- https://mmapped.blog/posts/44-tag-based-logging

The main reason for owning the DSL is that GammaLoop needs composable tag groups and stable semantics around presence, absence, and explicit boolean values.

- Use stable targets for broad subsystems such as `gammalooprs::uv::forest` or `gammalooprs::integrands::process::cross_section`.
- Add boolean event fields as composable tags.
- Filter by field presence with `target[{tag_a,tag_b,!tag_c}]=debug`, or by displayed boolean values with `target[{#tag_a,#!tag_c}]=debug`.

This is preferable to span-name filtering because span-based filters can admit child-library events emitted while the span is active. GammaLoop tag filters only look at the event target and the boolean tag fields on the event itself.

### Supported directive syntax

The supported directive forms are:

- `level`
- `target=level`
- `target[{tag_a,tag_b,!tag_c}]=level`
- `target[{#tag_a,#!tag_c}]=level`
- `target`

Notes:

- Directives are comma-separated.
- `target` is matched as a module-style prefix, so `gammalooprs::uv=debug` matches `gammalooprs::uv` and `gammalooprs::uv::forest`.
- A directive without an explicit `=level` is treated as `=trace`.
- Tag groups are conjunctions: `[{generation,uv}]` means both tags must match.
- `!tag` means the field is absent.
- `#tag` means the field is present with boolean value `true`.
- `#!tag` means the field is present with boolean value `false`.
- `field=value` matches an explicit field value.

### Tag matching model

Tags are represented by event fields.

- `#generation` at the callsite emits `generation = true`.
- In a directive, `generation` means the event must carry a field named `generation`.
- In a directive, `#generation` means the event must carry `generation = true`.
- In a directive, `!inspect` means the event must not carry a field named `inspect`.
- In a directive, `#!inspect` means the event must carry `inspect = false`.
- In a directive, `inspect=true`, `mode=summary`, or `label="soft region"` means the event must carry that value specifically.

This means the callsite controls both:

- whether a log is selectable by a presence/absence query
- whether a log can be selected by an explicit value query such as `inspect=false` or `mode=summary`

### Static vs dynamic filtering

Presence/absence-only directives stay on the lazy metadata path.

- `generation`
- `!inspect`
- `gammalooprs::uv::forest[{generation,uv,!dump}]=debug`

These can be decided from the callsite target and declared field names alone, so the filter can answer `always` or `never` at callsite registration time.

Value-matching directives are dynamic.

- `inspect=false`
- `#!inspect`
- `#generation`
- `mode=summary`
- `label="soft region"`

Those depend on the event instance, so they are evaluated in the per-event pass.

### Precedence

When multiple directives match the same event, GammaLoop prefers the most specific one:

- longer target prefix wins
- then the directive with more tag requirements wins
- then later directives win if the earlier specificity is identical

This lets broad directives such as `gammalooprs=info` coexist with narrow tag-specific directives such as `gammalooprs::uv::forest[{generation,uv,dump}]=debug`.

### Copy-pasteable display tags

Display logs render boolean fields next to the source as a directive tag group:

```text
@gammalooprs::uv::forest[{#generation, #uv, #!inspect}] DEBUG: message
```

Every boolean field is rendered this way, not only the standard tag vocabulary. `true` becomes `#field_name`; `false` becomes `#!field_name`. These boolean fields are removed from the normal display field table.

Use the `full` logging prefix when you want the rendered source to be a valid directive head:

```bash
gammaloop --logging-prefix full -l debug ...
```

Then the text after `@` and before the log level can be copied into a directive by adding `=debug`, for example:

```toml
display_directive = "gammalooprs::uv::forest[{#generation, #uv, #!inspect}]=debug"
```

This copy-paste form assumes the display source is the module target. Enabling `full_line_source` adds file and line information to the source display, which is useful for locating code but is not a valid directive target.

### Sink-only field prefixes

This repo also has sink-routing prefixes for fields:

- `file.<name>`: only rendered into the logfile/json sink
- `display.<name>`: only rendered into the stderr/display sink

Examples:

- `file.integrands = %...` will appear only in the file logger output
- `display.progress = %...` will appear only in the display logger output

These prefixes are a formatting/routing feature, not a separate event kind. They decide where a field is rendered after the event has been emitted. This lets a callsite attach large payloads such as `file.integrands` without dumping them to stderr.

### Pipeline-wide tag set

Prefer a small stable vocabulary that cuts across the whole pipeline.

Primary phase tags:

- `generation`: work that builds or prepares runtime objects before Monte Carlo integration starts
- `integration`: work performed while evaluating samples or managing the adaptive integration loop
- `profile`: UV/IR/profile-style diagnostic scans and their analysis
- `persistence`: reading or writing saved state, manifests, checkpoints, and exported results

Core domain tags:

- `uv`: ultraviolet counterterms, UV forests, UV profiles, or UV-specific generation/evaluation logic
- `ir`: infrared subtraction, IR profiles, or IR-specific generation/evaluation logic
- `subtraction`: threshold, UV, or IR subtraction logic when the main concern is subtraction rather than the specific regime
- `sampling`: channel choice, discrete axes, parameterizations, or sample-generation choices
- `stability`: precision escalation, retries, instability diagnosis, and numerical safety checks
- `observables`: histogramming, event-to-observable projection, and observable snapshot production
- `selectors`: event-selection logic and selector decisions
- `cache`: cache lookup, reuse, invalidation, or cache-debug instrumentation

Common work-unit tags:

- `graph`: the log is about one graph or graph-local data
- `group`: the log is about a graph group or another explicitly grouped aggregate
- `orientation`: the log is about orientation-dependent data or choosing/summing orientations
- `channel`: the log is about multi-channeling or a particular channel
- `cut`: the log is about a cut, raised cut, or cut-local computation
- `event`: the log is about generated/retained event objects or event-group processing
- `sample`: the log is about one evaluation sample, its coordinates, or per-sample intermediate values
- `iteration`: the log is about one adaptive integration iteration or iteration-level summaries
- `term`: the log is about one algebraic term, summand, or term-local contribution

Common purpose tags:

- `solver`: the log is about root-finding, linear solves, fitting, or similar numerical solver state
- `compile`: the log is about evaluator/code generation or compilation-oriented preparation
- `inspect`: the log exposes detailed intermediate state for debugging, rather than a high-level milestone
- `summary`: the log is a compact roll-up rather than a step-by-step trace
- `dump`: the log emits large or structured payloads such as expressions, tables, or serialized views

### Tag boundaries

Use the most general tag that accurately captures the reason you want to turn the log on or off.

- Prefer `generation` over `compile` when the message is a broad generation milestone.
- Add `compile` only when the message is specifically about evaluator construction or compilation-like work.
- Prefer `uv` or `ir` when the regime matters to the user; use `subtraction` when the subtraction mechanism is the real concern.
- Prefer `inspect` for verbose intermediate values that are mainly useful while debugging internals.
- Add `dump` when the payload is large enough that users may want to suppress it separately from lighter debug logs.
- Do not use `graph`, `cut`, or `sample` as substitutes for `generation` or `integration`; they refine phase tags rather than replace them.

### Naming guidance

- Use phase tags first. They answer "when in the pipeline did this happen?"
- Use domain tags second. They answer "what subsystem is this about?"
- Use work-unit tags third. They answer "what object is being worked on?"
- Use purpose tags last for optional refinement.
- Prefer tags that are meaningful across multiple modules.
- Avoid tags that merely restate a function name.
- Avoid tags tied to one internal representation unless that representation is a stable concept in user-facing debugging.

`parametric` is usually not a core pipeline tag. It is acceptable as a local refinement when needed, but should not be treated as part of the primary vocabulary unless it becomes a consistently useful cross-cutting concept.

### Example queries

- all generation-time UV forest logs:
  - `gammalooprs::uv::forest[{generation,uv}]=debug`
- only the large UV forest dumps for orientation-dependent generation logs:
  - `gammalooprs::uv::forest[{generation,uv,orientation,dump}]=debug`
- UV forest term-by-term generation logs:
  - `gammalooprs::uv::forest[{generation,uv,term}]=debug`
- all integration-time subtraction debug in amplitude evaluation:
  - `gammalooprs::integrands::process::amplitude[{integration,subtraction}]=debug`
- per-sample integration inspection in cross-section evaluation:
  - `gammalooprs::integrands::process::cross_section[{integration,sample,inspect}]=debug`
- all integration-time cut solver debug in cross-section evaluation:
  - `gammalooprs::integrands::process::cross_section[{integration,cut,solver}]=debug`
- integration-time debug excluding inspect-tagged logs:
  - `gammalooprs::integrands::process::cross_section[{integration,!inspect}]=debug`
- only events that explicitly set `inspect=true`:
  - `gammalooprs::integrands::process::cross_section[{integration,inspect=true}]=debug`
- only events with a specific textual mode field:
  - `gammalooprs::integrands::process::cross_section[{integration,mode=summary}]=debug`

Example display directive:

```toml
[cli_settings.global]
display_directive = "gammaloop_api=info,gammalooprs=info,symbolica=off,poly::gcd=off,gammalooprs::uv::forest[{generation,uv,orientation,dump}]=debug"
```

## Unsupported `EnvFilter` syntax

GammaLoop no longer treats the directive string as generic `EnvFilter` syntax.

In particular, do not rely on:

- span-name filters such as `target[span_name]=debug`
- generic field filters such as `[{field=value}]`
- `EnvFilter` regex semantics

If those are needed later, they must be added explicitly to the GammaLoop DSL.

## Release-build note

`gammalooprs` is compiled with `tracing` feature `release_max_level_info`, so `debug!` and `trace!` callsites in that crate are compiled out in release builds.
