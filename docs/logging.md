# Logging

This file summarizes the default tracing configuration used by the CLI/API at startup.

## Default directives

Normal session defaults come from `GlobalSettings`:

- `display_directive = "info"`
- `logfile_directive = "off"`

These are the application defaults written into `global_settings.toml` when defaults are shown or persisted.

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

## Release-build note

`gammalooprs` is compiled with `tracing` feature `release_max_level_info`, so `debug!` and `trace!` callsites in that crate are compiled out in release builds.
