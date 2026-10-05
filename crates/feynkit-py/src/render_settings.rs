//! Typed Python settings for the native diagram renderer.
use feynkit_graph::SceneOptions;
use linnest::svg::Config;
use linnet::half_edge::layout::impred::ImpredConfig;
use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};
use serde_json::{Map, Value, json};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
/// Immutable native graph layout overrides.
/// None leaves an option to the renderer; it does not force a default override.
/// Constructors and read-only properties expose options to help() and completion.
///
/// Examples
/// --------
/// >>> from symbolica.community.hepkit import LayoutSettings
/// >>> settings = LayoutSettings(impred_steps=100, impred_labels=True)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "LayoutSettings",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Debug, Default)]
pub struct PyLayoutSettings {
    /// Layout algorithm: "impred" preserves the embedding; "dot" uses layered ranks.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().layout_algo
    #[pyo3(get)]
    layout_algo: Option<String>,
    /// Horizontal spacing for layered layouts, in drawing units.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().tree_dx
    #[pyo3(get)]
    tree_dx: Option<f64>,
    /// Vertical spacing for layered layouts, in drawing units.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().tree_dy
    #[pyo3(get)]
    tree_dy: Option<f64>,
    /// Cooling schedule duration in reference steps; the native default is 500.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_steps
    #[pyo3(get)]
    impred_steps: Option<usize>,
    /// Maximum integration stride after settling; the native default is 2.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_step_scale
    #[pyo3(get)]
    impred_step_scale: Option<usize>,
    /// Target spacing in drawing units; the native default is 2.4.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_spacing
    #[pyo3(get)]
    impred_spacing: Option<f64>,
    /// Nonnegative repulsion multiplier; the native default is 2.5.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_repulsion
    #[pyo3(get)]
    impred_repulsion: Option<f64>,
    /// Nonnegative attraction multiplier; the native default is 2.5.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_attraction
    #[pyo3(get)]
    impred_attraction: Option<f64>,
    /// Parallel-edge attraction balance; the native default is 1.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_parallel_balance
    #[pyo3(get)]
    impred_parallel_balance: Option<f64>,
    /// External-leg pull multiplier; the native default is 0.45.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_pull
    #[pyo3(get)]
    impred_pull: Option<f64>,
    /// Balance of external-leg pulls; the native default is 1.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_pull_balance
    #[pyo3(get)]
    impred_pull_balance: Option<f64>,
    /// Maximum intermediate external-route points, from 0 to 3; default 2.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_external_max_points
    #[pyo3(get)]
    impred_external_max_points: Option<usize>,
    /// Subdivision threshold relative to spacing; the native default is 1.5.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_split_length_ratio
    #[pyo3(get)]
    impred_split_length_ratio: Option<f64>,
    /// Positive contraction threshold below the subdivision threshold; default 1.25.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_contract_chord_ratio
    #[pyo3(get)]
    impred_contract_chord_ratio: Option<f64>,
    /// Positive edge clearance in drawing units; the native default is 0.4.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_edge_clearance
    #[pyo3(get)]
    impred_edge_clearance: Option<f64>,
    /// Nonnegative node-edge repulsion multiplier; the native default is 4.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_node_edge_strength
    #[pyo3(get)]
    impred_node_edge_strength: Option<f64>,
    /// Refine the layout around measured labels; disabled by default.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_labels
    #[pyo3(get)]
    impred_labels: Option<bool>,
    /// Rotate toward the external-pull optimum; enabled by default.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> value = LayoutSettings().impred_level
    #[pyo3(get)]
    impred_level: Option<bool>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyLayoutSettings {
    /// Construct immutable overrides; None preserves the renderer's defaults.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> settings = LayoutSettings(impred_steps=100, impred_labels=True)
    ///
    /// Parameters
    /// ----------
    /// layout_algo : str or None, optional
    ///     Layout algorithm: "impred" preserves the embedding; "dot" uses layered ranks.
    /// tree_dx : float or None, optional
    ///     Horizontal spacing for layered layouts, in drawing units.
    /// tree_dy : float or None, optional
    ///     Vertical spacing for layered layouts, in drawing units.
    /// impred_steps : int or None, optional
    ///     Cooling schedule duration in reference steps; the native default is 500.
    /// impred_step_scale : int or None, optional
    ///     Maximum integration stride after settling; the native default is 2.
    /// impred_spacing : float or None, optional
    ///     Target spacing in drawing units; the native default is 2.4.
    /// impred_repulsion : float or None, optional
    ///     Nonnegative repulsion multiplier; the native default is 2.5.
    /// impred_attraction : float or None, optional
    ///     Nonnegative attraction multiplier; the native default is 2.5.
    /// impred_parallel_balance : float or None, optional
    ///     Parallel-edge attraction balance; the native default is 1.
    /// impred_pull : float or None, optional
    ///     External-leg pull multiplier; the native default is 0.45.
    /// impred_pull_balance : float or None, optional
    ///     Balance of external-leg pulls; the native default is 1.
    /// impred_external_max_points : int or None, optional
    ///     Maximum intermediate external-route points, from 0 to 3; default 2.
    /// impred_split_length_ratio : float or None, optional
    ///     Subdivision threshold relative to spacing; the native default is 1.5.
    /// impred_contract_chord_ratio : float or None, optional
    ///     Positive contraction threshold below the subdivision threshold; default 1.25.
    /// impred_edge_clearance : float or None, optional
    ///     Positive edge clearance in drawing units; the native default is 0.4.
    /// impred_node_edge_strength : float or None, optional
    ///     Nonnegative node-edge repulsion multiplier; the native default is 4.
    /// impred_labels : bool or None, optional
    ///     Refine the layout around measured labels; disabled by default.
    /// impred_level : bool or None, optional
    ///     Rotate toward the external-pull optimum; enabled by default.
    #[new]
    #[pyo3(signature = (*, layout_algo=None, tree_dx=None, tree_dy=None, impred_steps=None, impred_step_scale=None, impred_spacing=None, impred_repulsion=None, impred_attraction=None, impred_parallel_balance=None, impred_pull=None, impred_pull_balance=None, impred_external_max_points=None, impred_split_length_ratio=None, impred_contract_chord_ratio=None, impred_edge_clearance=None, impred_node_edge_strength=None, impred_labels=None, impred_level=None))]
    #[allow(clippy::too_many_arguments)] // Independently discoverable Python settings.
    fn new(
        #[gen_stub(override_type(type_repr="typing.Literal['dot', 'impred'] | None", imports=("typing")))]
        layout_algo: Option<String>,
        tree_dx: Option<f64>,
        tree_dy: Option<f64>,
        impred_steps: Option<usize>,
        impred_step_scale: Option<usize>,
        impred_spacing: Option<f64>,
        impred_repulsion: Option<f64>,
        impred_attraction: Option<f64>,
        impred_parallel_balance: Option<f64>,
        impred_pull: Option<f64>,
        impred_pull_balance: Option<f64>,
        impred_external_max_points: Option<usize>,
        impred_split_length_ratio: Option<f64>,
        impred_contract_chord_ratio: Option<f64>,
        impred_edge_clearance: Option<f64>,
        impred_node_edge_strength: Option<f64>,
        impred_labels: Option<bool>,
        impred_level: Option<bool>,
    ) -> PyResult<Self> {
        let settings = Self {
            layout_algo,
            tree_dx,
            tree_dy,
            impred_steps,
            impred_step_scale,
            impred_spacing,
            impred_repulsion,
            impred_attraction,
            impred_parallel_balance,
            impred_pull,
            impred_pull_balance,
            impred_external_max_points,
            impred_split_length_ratio,
            impred_contract_chord_ratio,
            impred_edge_clearance,
            impred_node_edge_strength,
            impred_labels,
            impred_level,
        };
        settings.validate()?;
        Ok(settings)
    }

    /// Inspect the selected overrides without rendering a diagram.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import LayoutSettings
    /// >>> text = repr(LayoutSettings(impred_steps=100, impred_labels=True))
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let fields = PyDict::new(py);
        if let Some(value) = &self.layout_algo {
            fields.set_item("layout_algo", value.clone())?;
        }
        if let Some(value) = &self.tree_dx {
            fields.set_item("tree_dx", *value)?;
        }
        if let Some(value) = &self.tree_dy {
            fields.set_item("tree_dy", *value)?;
        }
        if let Some(value) = &self.impred_steps {
            fields.set_item("impred_steps", *value)?;
        }
        if let Some(value) = &self.impred_step_scale {
            fields.set_item("impred_step_scale", *value)?;
        }
        if let Some(value) = &self.impred_spacing {
            fields.set_item("impred_spacing", *value)?;
        }
        if let Some(value) = &self.impred_repulsion {
            fields.set_item("impred_repulsion", *value)?;
        }
        if let Some(value) = &self.impred_attraction {
            fields.set_item("impred_attraction", *value)?;
        }
        if let Some(value) = &self.impred_parallel_balance {
            fields.set_item("impred_parallel_balance", *value)?;
        }
        if let Some(value) = &self.impred_pull {
            fields.set_item("impred_pull", *value)?;
        }
        if let Some(value) = &self.impred_pull_balance {
            fields.set_item("impred_pull_balance", *value)?;
        }
        if let Some(value) = &self.impred_external_max_points {
            fields.set_item("impred_external_max_points", *value)?;
        }
        if let Some(value) = &self.impred_split_length_ratio {
            fields.set_item("impred_split_length_ratio", *value)?;
        }
        if let Some(value) = &self.impred_contract_chord_ratio {
            fields.set_item("impred_contract_chord_ratio", *value)?;
        }
        if let Some(value) = &self.impred_edge_clearance {
            fields.set_item("impred_edge_clearance", *value)?;
        }
        if let Some(value) = &self.impred_node_edge_strength {
            fields.set_item("impred_node_edge_strength", *value)?;
        }
        if let Some(value) = &self.impred_labels {
            fields.set_item("impred_labels", *value)?;
        }
        if let Some(value) = &self.impred_level {
            fields.set_item("impred_level", *value)?;
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
        Ok(format!("LayoutSettings({})", arguments.join(", ")))
    }
}
/// Immutable stroke overrides shared by vertices and edges.
/// None leaves an option to the renderer; it does not force a default override.
/// Constructors and read-only properties expose options to help() and completion.
///
/// Examples
/// --------
/// >>> from symbolica.community.hepkit import StrokeStyle
/// >>> settings = StrokeStyle(paint="#6f4d85", thickness=1.5)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "StrokeStyle",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Debug, Default)]
pub struct PyStrokeStyle {
    /// Stroke paint as a CSS color, for example "#6f4d85".
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import StrokeStyle
    /// >>> value = StrokeStyle().paint
    #[pyo3(get)]
    paint: Option<String>,
    /// Finite nonnegative stroke width in points.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import StrokeStyle
    /// >>> value = StrokeStyle().thickness
    #[pyo3(get)]
    thickness: Option<f64>,
    /// Line pattern: "solid", "dotted", or "dashed".
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import StrokeStyle
    /// >>> value = StrokeStyle().dash
    #[pyo3(get)]
    dash: Option<String>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyStrokeStyle {
    /// Construct immutable overrides; None preserves the renderer's defaults.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import StrokeStyle
    /// >>> settings = StrokeStyle(paint="#6f4d85", thickness=1.5)
    ///
    /// Parameters
    /// ----------
    /// paint : str or None, optional
    ///     Stroke paint as a CSS color, for example "#6f4d85".
    /// thickness : float or None, optional
    ///     Finite nonnegative stroke width in points.
    /// dash : str or None, optional
    ///     Line pattern: "solid", "dotted", or "dashed".
    #[new]
    #[pyo3(signature = (*, paint=None, thickness=None, dash=None))]
    fn new(
        paint: Option<String>,
        thickness: Option<f64>,
        #[gen_stub(override_type(type_repr="typing.Literal['solid', 'dotted', 'dashed'] | None", imports=("typing")))]
        dash: Option<String>,
    ) -> PyResult<Self> {
        let settings = Self {
            paint,
            thickness,
            dash,
        };
        settings.validate()?;
        Ok(settings)
    }

    /// Inspect the selected overrides without rendering a diagram.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import StrokeStyle
    /// >>> text = repr(StrokeStyle(paint="#6f4d85", thickness=1.5))
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let fields = PyDict::new(py);
        if let Some(value) = &self.paint {
            fields.set_item("paint", value.clone())?;
        }
        if let Some(value) = &self.thickness {
            fields.set_item("thickness", *value)?;
        }
        if let Some(value) = &self.dash {
            fields.set_item("dash", value.clone())?;
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
        Ok(format!("StrokeStyle({})", arguments.join(", ")))
    }
}
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
impl PyLayoutSettings {
    fn options(&self) -> Map<String, Value> {
        let mut options = Map::from_iter([
            ("layout-algo".into(), json!(self.layout_algo)),
            ("tree-dx".into(), json!(self.tree_dx)),
            ("tree-dy".into(), json!(self.tree_dy)),
            ("impred-steps".into(), json!(self.impred_steps)),
            ("impred-step-scale".into(), json!(self.impred_step_scale)),
            ("impred-spacing".into(), json!(self.impred_spacing)),
            ("impred-repulsion".into(), json!(self.impred_repulsion)),
            ("impred-attraction".into(), json!(self.impred_attraction)),
            (
                "impred-parallel-balance".into(),
                json!(self.impred_parallel_balance),
            ),
            ("impred-pull".into(), json!(self.impred_pull)),
            (
                "impred-pull-balance".into(),
                json!(self.impred_pull_balance),
            ),
            (
                "impred-external-max-points".into(),
                json!(self.impred_external_max_points),
            ),
            (
                "impred-split-length-ratio".into(),
                json!(self.impred_split_length_ratio),
            ),
            (
                "impred-contract-chord-ratio".into(),
                json!(self.impred_contract_chord_ratio),
            ),
            (
                "impred-edge-clearance".into(),
                json!(self.impred_edge_clearance),
            ),
            (
                "impred-node-edge-strength".into(),
                json!(self.impred_node_edge_strength),
            ),
            ("impred-labels".into(), json!(self.impred_labels)),
            ("impred-level".into(), json!(self.impred_level)),
        ]);
        options.retain(|_, value| !value.is_null());
        options
    }
}
impl PyStrokeStyle {
    fn options(&self) -> Map<String, Value> {
        let mut options = Map::from_iter([
            ("paint".into(), json!(self.paint)),
            ("thickness".into(), json!(self.thickness)),
            ("dash".into(), json!(self.dash)),
        ]);
        options.retain(|_, value| !value.is_null());
        options
    }
}

impl PyLayoutSettings {
    fn validate(&self) -> PyResult<()> {
        if self
            .layout_algo
            .as_deref()
            .is_some_and(|algo| !matches!(algo, "dot" | "impred"))
        {
            return Err(PyValueError::new_err(
                "layout_algo must be 'dot' or 'impred'",
            ));
        }
        for (name, value) in [("tree_dx", self.tree_dx), ("tree_dy", self.tree_dy)] {
            if value.is_some_and(|value| !value.is_finite() || value <= 0.0) {
                return Err(PyValueError::new_err(format!(
                    "{name} must be finite and positive"
                )));
            }
        }
        let defaults = ImpredConfig::default();
        ImpredConfig {
            steps: self.impred_steps.unwrap_or(defaults.steps),
            step_scale: self.impred_step_scale.unwrap_or(defaults.step_scale),
            target: self.impred_spacing.unwrap_or(defaults.target),
            repulsion: self.impred_repulsion.unwrap_or(defaults.repulsion),
            attraction: self.impred_attraction.unwrap_or(defaults.attraction),
            parallel_attraction_balance: self
                .impred_parallel_balance
                .unwrap_or(defaults.parallel_attraction_balance),
            pull: self.impred_pull.unwrap_or(defaults.pull),
            pull_balance: self.impred_pull_balance.unwrap_or(defaults.pull_balance),
            external_max_points: self
                .impred_external_max_points
                .unwrap_or(defaults.external_max_points),
            split_length_ratio: self
                .impred_split_length_ratio
                .unwrap_or(defaults.split_length_ratio),
            contract_chord_ratio: self
                .impred_contract_chord_ratio
                .unwrap_or(defaults.contract_chord_ratio),
            edge_clearance: self
                .impred_edge_clearance
                .unwrap_or(defaults.edge_clearance),
            node_edge_strength: self
                .impred_node_edge_strength
                .unwrap_or(defaults.node_edge_strength),
            labels: self.impred_labels.unwrap_or(defaults.labels),
            level: self.impred_level.unwrap_or(defaults.level),
            ..defaults
        }
        .validate()
        .map_err(PyValueError::new_err)
    }
}

impl PyStrokeStyle {
    fn validate(&self) -> PyResult<()> {
        if self
            .thickness
            .is_some_and(|width| !width.is_finite() || width < 0.0)
        {
            return Err(PyValueError::new_err(
                "thickness must be finite and nonnegative",
            ));
        }
        if self
            .dash
            .as_deref()
            .is_some_and(|dash| !matches!(dash, "solid" | "dotted" | "dashed"))
        {
            return Err(PyValueError::new_err(
                "dash must be 'solid', 'dotted', or 'dashed'",
            ));
        }
        Ok(())
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
        if self
            .node_radius
            .is_some_and(|radius| !radius.is_finite() || radius < 0.0)
        {
            return Err(PyValueError::new_err(
                "node_radius must be finite and nonnegative",
            ));
        }
        Ok(())
    }

    fn config(&self) -> Config {
        let mut drawing = Map::new();
        if let Some(radius) = self.node_radius {
            drawing.insert("node-radius".into(), json!(radius));
        }
        let mut node_style = Map::new();
        if let Some(fill) = &self.node_fill {
            node_style.insert("fill".into(), json!(fill));
        }
        if let Some(stroke) = &self.node_stroke {
            node_style.insert("stroke".into(), Value::Object(stroke.options()));
        }
        let mut style = Map::new();
        if !node_style.is_empty() {
            style.insert("node-style".into(), Value::Object(node_style));
        }
        if let Some(stroke) = &self.edge_stroke {
            style.insert("edge-style".into(), json!({"stroke": stroke.options()}));
        }
        Config {
            title: self.title.clone(),
            layouts: self
                .layout
                .as_ref()
                .map_or_else(Map::new, PyLayoutSettings::options),
            drawing,
            style,
            template_options: Map::new(),
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyLayoutSettings>()?;
    module.add_class::<PyStrokeStyle>()?;
    module.add_class::<PyRenderSettings>()?;
    Ok(())
}
