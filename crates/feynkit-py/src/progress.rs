//! Shared notebook presentation for caller-owned native progress streams.

use pyo3::{prelude::*, types::PyDict};

struct Indicator {
    context: Py<PyAny>,
    widget: Py<PyAny>,
    total: Option<usize>,
    completed: usize,
}

/// A live Marimo indicator, independent of the native operation's event type.
///
/// Callers retain their own snapshots, throttling, callbacks and cancellation.
/// No notebook dependency is imported in ordinary scripts. Finish the display
/// on both success and error, passing the original error before returning it.
pub struct MarimoProgress {
    status: Py<PyAny>,
    flush: Py<PyAny>,
    active: Option<Indicator>,
}

impl MarimoProgress {
    /// Start a spinner only when an already imported Marimo is running a notebook.
    pub fn detect(py: Python<'_>, title: &str, subtitle: &str) -> PyResult<Option<Self>> {
        let modules = py.import("sys")?.getattr("modules")?;
        let marimo = modules.cast::<PyDict>()?.get_item("marimo")?;
        let Some(marimo) = marimo.filter(|module| !module.is_none()) else {
            return Ok(None);
        };
        if !marimo.call_method0("running_in_notebook")?.is_truthy()? {
            return Ok(None);
        }
        let flush = marimo
            .getattr("output")?
            .getattr("_output")?
            .getattr("flush")?
            .unbind();
        let status = marimo.getattr("status")?;
        let context = status.call_method1("spinner", (title, subtitle))?;
        let widget = context.call_method0("__enter__")?.unbind();
        Ok(Some(Self {
            status: status.unbind(),
            flush,
            active: Some(Indicator {
                context: context.unbind(),
                widget,
                total: None,
                completed: 0,
            }),
        }))
    }

    /// Display absolute stage counts, including resets between stages.
    ///
    /// Unknown or empty work uses a spinner. Positive totals use a bar with no
    /// speculative rate or ETA. Each delivered update is explicitly flushed so
    /// Marimo's own throttle cannot hide a stage's only update.
    pub fn update(
        &mut self,
        py: Python<'_>,
        title: &str,
        subtitle: &str,
        completed: usize,
        total: Option<usize>,
    ) -> PyResult<()> {
        let total = total.filter(|total| *total > 0);
        if self.active.as_ref().map(|state| state.total) != Some(total) {
            if let Some(previous) = self.active.take() {
                previous
                    .context
                    .call_method1(py, "__exit__", (py.None(), py.None(), py.None()))?;
            }
            let kwargs = PyDict::new(py);
            kwargs.set_item("title", title)?;
            kwargs.set_item("subtitle", subtitle)?;
            kwargs.set_item("remove_on_exit", true)?;
            let kind = if let Some(total) = total {
                kwargs.set_item("total", total)?;
                kwargs.set_item("show_rate", false)?;
                kwargs.set_item("show_eta", false)?;
                "progress_bar"
            } else {
                "spinner"
            };
            let context = self.status.call_method(py, kind, (), Some(&kwargs))?;
            let widget = context.call_method0(py, "__enter__")?;
            self.active = Some(Indicator {
                context,
                widget,
                total,
                completed: 0,
            });
        }
        let state = self.active.as_mut().unwrap();
        if total.is_some() {
            let increment = completed as i128 - state.completed as i128;
            state
                .widget
                .call_method1(py, "update", (increment, title, subtitle))?;
        } else {
            state.widget.call_method1(py, "update", (title, subtitle))?;
        }
        state.completed = completed;
        self.flush.call0(py)?;
        Ok(())
    }

    /// Close the active context once; cleanup never replaces an original error.
    pub fn finish(
        &mut self,
        py: Python<'_>,
        error: Option<&PyErr>,
        stopped_title: &str,
    ) -> PyResult<()> {
        let Some(state) = self.active.take() else {
            return Ok(());
        };
        if let Some(error) = error {
            let subtitle = error.to_string();
            let _ = if state.total.is_some() {
                state
                    .widget
                    .call_method1(py, "update", (0, stopped_title, subtitle))
            } else {
                state
                    .widget
                    .call_method1(py, "update", (stopped_title, subtitle))
            };
            let _ = state.context.call_method1(
                py,
                "__exit__",
                (error.get_type(py), error.value(py), error.traceback(py)),
            );
        } else {
            state
                .context
                .call_method1(py, "__exit__", (py.None(), py.None(), py.None()))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use super::*;

    #[test]
    fn shared_display_preserves_counts_and_error_cleanup() {
        Python::initialize();
        Python::attach(|py| {
            let locals = PyDict::new(py);
            let run =
                |code: &str| py.run(&CString::new(code).unwrap(), Some(&locals), Some(&locals));
            run(r#"
from types import SimpleNamespace
from unittest.mock import MagicMock
context = MagicMock()
indicator = context.__enter__.return_value
marimo = SimpleNamespace(
    running_in_notebook=MagicMock(return_value=True),
    status=SimpleNamespace(
        spinner=MagicMock(return_value=context),
        progress_bar=MagicMock(return_value=context)),
    output=SimpleNamespace(_output=SimpleNamespace(flush=MagicMock())))
"#)
            .unwrap();
            let modules = py.import("sys").unwrap().getattr("modules").unwrap();
            let modules = modules.cast::<PyDict>().unwrap();
            let previous = modules.get_item("marimo").unwrap();
            let result = (|| -> PyResult<()> {
                modules.set_item("marimo", py.None())?;
                assert!(MarimoProgress::detect(py, "test", "start")?.is_none());
                modules.set_item("marimo", locals.get_item("marimo")?.unwrap())?;
                run("marimo.running_in_notebook.return_value = False")?;
                assert!(MarimoProgress::detect(py, "test", "start")?.is_none());
                run(
                    "marimo.status.spinner.assert_not_called()\nmarimo.running_in_notebook.return_value = True",
                )?;
                let mut display = MarimoProgress::detect(py, "test", "start")?.unwrap();
                display.update(py, "empty", "zero", 0, Some(0))?;
                display.update(py, "work", "three", 3, Some(5))?;
                display.update(py, "next stage", "one", 1, Some(5))?;
                display.update(py, "unknown", "four", 4, None)?;
                display.finish(py, None, "stopped")?;
                display.finish(py, None, "stopped")?;
                run(r#"
assert [call.args for call in indicator.update.call_args_list] == [
    ('empty', 'zero'), (3, 'work', 'three'),
    (-2, 'next stage', 'one'), ('unknown', 'four')]
assert marimo.output._output.flush.call_count == 4
assert context.__enter__.call_count == context.__exit__.call_count == 3
context.__exit__.assert_called_with(None, None, None)
context.reset_mock()
indicator.update.side_effect = RuntimeError('widget failed')
context.__exit__.side_effect = RuntimeError('cleanup failed')
"#)?;
                let mut display = MarimoProgress::detect(py, "test", "start")?.unwrap();
                let error = pyo3::exceptions::PyKeyboardInterrupt::new_err("interrupted");
                locals.set_item("original_error", error.value(py))?;
                display.finish(py, Some(&error), "stopped")?;
                run(r#"
assert context.__exit__.call_count == 1
assert context.__exit__.call_args.args[:2] == (KeyboardInterrupt, original_error)
context.reset_mock()
indicator.update.side_effect = None
context.__exit__.side_effect = None
marimo.status.progress_bar.side_effect = RuntimeError('cannot open bar')
"#)?;
                let mut display = MarimoProgress::detect(py, "test", "start")?.unwrap();
                let error = display.update(py, "work", "one", 1, Some(2)).unwrap_err();
                display.finish(py, Some(&error), "stopped")?;
                run("assert context.__enter__.call_count == context.__exit__.call_count == 1")?;
                Ok(())
            })();
            if let Some(previous) = previous {
                modules.set_item("marimo", previous).unwrap();
            } else {
                modules.del_item("marimo").unwrap();
            }
            result.unwrap();
        });
    }
}
