//! Shared notebook presentation for caller-owned native progress streams.

use pyo3::prelude::{Py, PyAny, PyAnyMethods, PyDictMethods, PyErr, PyResult, Python};
use pyo3::types::PyDict;

/// A live Marimo spinner, independent of the native operation's event type.
///
/// One spinner serves the whole operation: its title names the current stage and
/// its subtitle carries the counts, so stage changes never replace the element.
/// Callers retain their own snapshots, throttling, callbacks and cancellation.
/// No notebook dependency is imported in ordinary scripts. Finish the display
/// on both success and error, passing the original error before returning it.
pub struct MarimoProgress {
    flush: Py<PyAny>,
    /// The entered spinner context and its widget, until the display finishes.
    active: Option<(Py<PyAny>, Py<PyAny>)>,
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
        let context = marimo
            .getattr("status")?
            .call_method1("spinner", (title, subtitle))?;
        let widget = context.call_method0("__enter__")?.unbind();
        Ok(Some(Self {
            flush,
            active: Some((context.unbind(), widget)),
        }))
    }

    /// Show a stage and its counts.
    ///
    /// Each update is explicitly flushed so Marimo's own throttle cannot hide a
    /// stage's only update.
    pub fn update(&mut self, py: Python<'_>, title: &str, subtitle: &str) -> PyResult<()> {
        let Some((_, widget)) = &self.active else {
            return Ok(());
        };
        widget.call_method1(py, "update", (title, subtitle))?;
        self.flush.call0(py)?;
        Ok(())
    }

    /// Close the spinner once; cleanup never replaces an original error.
    pub fn finish(
        &mut self,
        py: Python<'_>,
        error: Option<&PyErr>,
        stopped_title: &str,
    ) -> PyResult<()> {
        let Some((context, widget)) = self.active.take() else {
            return Ok(());
        };
        if let Some(error) = error {
            let _ = widget.call_method1(py, "update", (stopped_title, error.to_string()));
            let _ = context.call_method1(
                py,
                "__exit__",
                (error.get_type(py), error.value(py), error.traceback(py)),
            );
        } else {
            context.call_method1(py, "__exit__", (py.None(), py.None(), py.None()))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use super::*;
    use pyo3::exceptions::PyKeyboardInterrupt;

    #[test]
    fn one_spinner_shows_every_update_and_error_cleanup() {
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
    status=SimpleNamespace(spinner=MagicMock(return_value=context)),
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
                display.update(py, "work", "3 / 5 processed")?;
                display.update(py, "next stage", "1 / 5 processed")?;
                display.update(py, "unknown", "4 processed")?;
                display.finish(py, None, "stopped")?;
                display.finish(py, None, "stopped")?;
                display.update(py, "late", "ignored")?;
                run(r#"
marimo.status.spinner.assert_called_once_with("test", "start")
assert [call.args for call in indicator.update.call_args_list] == [
    ('work', '3 / 5 processed'), ('next stage', '1 / 5 processed'), ('unknown', '4 processed')]
assert marimo.output._output.flush.call_count == 3
assert context.__enter__.call_count == context.__exit__.call_count == 1
context.__exit__.assert_called_with(None, None, None)
context.reset_mock()
indicator.update.side_effect = RuntimeError('widget failed')
context.__exit__.side_effect = RuntimeError('cleanup failed')
"#)?;
                let mut display = MarimoProgress::detect(py, "test", "start")?.unwrap();
                let error = PyKeyboardInterrupt::new_err("interrupted");
                locals.set_item("original_error", error.value(py))?;
                display.finish(py, Some(&error), "stopped")?;
                run(r#"
assert context.__exit__.call_count == 1
assert context.__exit__.call_args.args[:2] == (KeyboardInterrupt, original_error)
context.reset_mock()
context.__exit__.side_effect = None
"#)?;
                let mut display = MarimoProgress::detect(py, "test", "start")?.unwrap();
                let error = display.update(py, "work", "one").unwrap_err();
                display.finish(py, Some(&error), "stopped")?;
                run(
                    "assert context.__enter__.call_count == context.__exit__.call_count == 1\nindicator.update.side_effect = None",
                )?;
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
