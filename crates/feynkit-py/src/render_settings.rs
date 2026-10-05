//! Typed Python settings for the native diagram renderer.
use feynkit_graph::SceneOptions;
use linnest::svg::Config;
use pyo3::{prelude::*, types::PyDict};
pub use spynso3::display::graph::{PyLayoutSettings, PyStrokeStyle};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
/// Immutable presentation settings for diagrams, subgraphs, and process schematics.
/// None leaves an option to the renderer; it does not force a default override.
/// Constructors and read-only properties expose options to help() and completion.
///
/// Examples
/// --------
/// >>> from symbolica.community.hepkit import RenderSettings
/// >>> settings = RenderSettings(node_radius=5, show_particle=False)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "RenderSettings",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Debug, Default)]
pub struct PyRenderSettings {
    /// Plain-text title above the drawing.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().title
    #[pyo3(get)]
    title: Option<String>,
    /// Native graph layout options; omitted values retain renderer defaults.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().layout
    #[pyo3(get)]
    layout: Option<PyLayoutSettings>,
    /// Finite nonnegative vertex radius in drawing units.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().node_radius
    #[pyo3(get)]
    node_radius: Option<f64>,
    /// Vertex fill as a CSS color.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().node_fill
    #[pyo3(get)]
    node_fill: Option<String>,
    /// Vertex outline paint, width, and dash pattern.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().node_stroke
    #[pyo3(get)]
    node_stroke: Option<PyStrokeStyle>,
    /// Edge paint, width, and dash pattern.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().edge_stroke
    #[pyo3(get)]
    edge_stroke: Option<PyStrokeStyle>,
    /// Show particle labels; enabled by default.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().show_particle
    #[pyo3(get)]
    show_particle: Option<bool>,
    /// Show momentum labels; None follows momenta or the supplied momentum basis.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().show_momentum
    #[pyo3(get)]
    show_momentum: Option<bool>,
    /// Show native edge IDs; disabled by default.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().show_edge_index
    #[pyo3(get)]
    show_edge_index: Option<bool>,
    /// Show native vertex IDs; disabled by default.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().show_node_index
    #[pyo3(get)]
    show_node_index: Option<bool>,
    /// Draw momentum arrows; None follows momenta or the supplied momentum basis.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().momentum_arrows
    #[pyo3(get)]
    momentum_arrows: Option<bool>,
    /// Open sewn initial-state connections in cross sections; enabled by default.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().split_initial_state
    #[pyo3(get)]
    split_initial_state: Option<bool>,
    /// Show both vertex and edge IDs; disabled by default.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> value = RenderSettings().debug
    #[pyo3(get)]
    debug: Option<bool>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyRenderSettings {
    /// Construct immutable overrides; None preserves the renderer's defaults.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> settings = RenderSettings(node_radius=5, show_particle=False)
    ///
    /// Parameters
    /// ----------
    /// title : str or None, optional
    ///     Plain-text title above the drawing.
    /// layout : LayoutSettings or None, optional
    ///     Native graph layout options; omitted values retain renderer defaults.
    /// node_radius : float or None, optional
    ///     Finite nonnegative vertex radius in drawing units.
    /// node_fill : str or None, optional
    ///     Vertex fill as a CSS color.
    /// node_stroke : StrokeStyle or None, optional
    ///     Vertex outline paint, width, and dash pattern.
    /// edge_stroke : StrokeStyle or None, optional
    ///     Edge paint, width, and dash pattern.
    /// show_particle : bool or None, optional
    ///     Show particle labels; enabled by default.
    /// show_momentum : bool or None, optional
    ///     Show momentum labels; None follows momenta or the supplied momentum basis.
    /// show_edge_index : bool or None, optional
    ///     Show native edge IDs; disabled by default.
    /// show_node_index : bool or None, optional
    ///     Show native vertex IDs; disabled by default.
    /// momentum_arrows : bool or None, optional
    ///     Draw momentum arrows; None follows momenta or the supplied momentum basis.
    /// split_initial_state : bool or None, optional
    ///     Open sewn initial-state connections in cross sections; enabled by default.
    /// debug : bool or None, optional
    ///     Show both vertex and edge IDs; disabled by default.
    #[new]
    #[pyo3(signature = (*, title=None, layout=None, node_radius=None, node_fill=None, node_stroke=None, edge_stroke=None, show_particle=None, show_momentum=None, show_edge_index=None, show_node_index=None, momentum_arrows=None, split_initial_state=None, debug=None))]
    #[allow(clippy::too_many_arguments)] // Independently discoverable Python settings.
    fn new(
        title: Option<String>,
        layout: Option<PyLayoutSettings>,
        node_radius: Option<f64>,
        node_fill: Option<String>,
        node_stroke: Option<PyStrokeStyle>,
        edge_stroke: Option<PyStrokeStyle>,
        show_particle: Option<bool>,
        show_momentum: Option<bool>,
        show_edge_index: Option<bool>,
        show_node_index: Option<bool>,
        momentum_arrows: Option<bool>,
        split_initial_state: Option<bool>,
        debug: Option<bool>,
    ) -> PyResult<Self> {
        let settings = Self {
            title,
            layout,
            node_radius,
            node_fill,
            node_stroke,
            edge_stroke,
            show_particle,
            show_momentum,
            show_edge_index,
            show_node_index,
            momentum_arrows,
            split_initial_state,
            debug,
        };
        settings.validate()?;
        Ok(settings)
    }

    /// Inspect the selected overrides without rendering a diagram.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import RenderSettings
    /// >>> text = repr(RenderSettings(node_radius=5, show_particle=False))
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let fields = PyDict::new(py);
        if let Some(value) = &self.title {
            fields.set_item("title", value.clone())?;
        }
        if let Some(value) = &self.layout {
            fields.set_item("layout", value.clone())?;
        }
        if let Some(value) = &self.node_radius {
            fields.set_item("node_radius", *value)?;
        }
        if let Some(value) = &self.node_fill {
            fields.set_item("node_fill", value.clone())?;
        }
        if let Some(value) = &self.node_stroke {
            fields.set_item("node_stroke", value.clone())?;
        }
        if let Some(value) = &self.edge_stroke {
            fields.set_item("edge_stroke", value.clone())?;
        }
        if let Some(value) = &self.show_particle {
            fields.set_item("show_particle", *value)?;
        }
        if let Some(value) = &self.show_momentum {
            fields.set_item("show_momentum", *value)?;
        }
        if let Some(value) = &self.show_edge_index {
            fields.set_item("show_edge_index", *value)?;
        }
        if let Some(value) = &self.show_node_index {
            fields.set_item("show_node_index", *value)?;
        }
        if let Some(value) = &self.momentum_arrows {
            fields.set_item("momentum_arrows", *value)?;
        }
        if let Some(value) = &self.split_initial_state {
            fields.set_item("split_initial_state", *value)?;
        }
        if let Some(value) = &self.debug {
            fields.set_item("debug", *value)?;
        }
        let arguments = fields
            .iter()
            .map(|(key, value)| {
                Ok(format!(
                    "{}={}",
                    key.str()?.to_str()?,
                    value.repr()?.to_str()?
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(format!("RenderSettings({})", arguments.join(", ")))
    }
}
impl PyRenderSettings {
    /// Resolve physics defaults once and retain native layout/style overrides.
    pub(crate) fn resolve(settings: Option<&Self>, momenta: bool) -> (Config, SceneOptions) {
        let empty = Self::default();
        let settings = settings.unwrap_or(&empty);
        let defaults = SceneOptions::default();
        let options = SceneOptions {
            momentum_arrows: settings.momentum_arrows.unwrap_or(momenta),
            show_momentum: settings.show_momentum,
            show_particle: settings.show_particle.unwrap_or(defaults.show_particle),
            show_edge_index: settings
                .show_edge_index
                .or(settings.debug)
                .unwrap_or(defaults.show_edge_index),
            show_node_index: settings
                .show_node_index
                .or(settings.debug)
                .unwrap_or(defaults.show_node_index),
            split_initial_state: settings
                .split_initial_state
                .unwrap_or(defaults.split_initial_state),
            ..defaults
        };
        (settings.config(), options)
    }

    fn validate(&self) -> PyResult<()> {
        self.graph_settings().map(|_| ())
    }

    fn graph_settings(&self) -> PyResult<spynso3::display::graph::PyRenderSettings> {
        spynso3::display::graph::PyRenderSettings::new(
            self.title.clone(),
            self.layout.clone(),
            self.node_radius,
            self.node_fill.clone(),
            self.node_stroke.clone(),
            self.edge_stroke.clone(),
        )
    }

    fn config(&self) -> Config {
        // Constructor validation and frozen fields keep these overrides valid.
        self.graph_settings()
            .expect("validated graph settings")
            .config()
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyLayoutSettings>()?;
    module.add_class::<PyStrokeStyle>()?;
    module.add_class::<PyRenderSettings>()?;
    Ok(())
}
