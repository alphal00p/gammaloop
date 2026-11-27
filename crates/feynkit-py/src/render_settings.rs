//! Physics presentation options, independent of generic graph rendering.
use feynkit_graph::{HighlightOverlap, HullOptions, MINIATURE_STAMP_WEIGHT, SceneOptions};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Particle and momentum presentation for Feynman diagrams.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> style = hep.DiagramStyle(show_momentum=True)
/// >>> style.show_momentum
/// True
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DiagramStyle",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object,
    get_all
)]
#[derive(Clone, Debug, Default)]
pub struct PyDiagramStyle {
    /// Show particle labels; None preserves the renderer default.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(show_particle=True).show_particle
    /// True
    show_particle: Option<bool>,
    /// Show momentum labels; None follows the render request.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(show_momentum=True).show_momentum
    /// True
    show_momentum: Option<bool>,
    /// Show the physical edge indices.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(show_edge_index=True).show_edge_index
    /// True
    show_edge_index: Option<bool>,
    /// Show the physical vertex indices.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(show_node_index=True).show_node_index
    /// True
    show_node_index: Option<bool>,
    /// Draw momentum arrows beside particle lines.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(momentum_arrows=True).momentum_arrows
    /// True
    momentum_arrows: Option<bool>,
    /// Separate the incoming states of a cross-section drawing.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(split_initial_state=True).split_initial_state
    /// True
    split_initial_state: Option<bool>,
    /// Show both vertex and edge indices unless individually overridden.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(debug=True).debug
    /// True
    debug: Option<bool>,
    /// Use stamp recolouring by default.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(miniature=True).miniature
    /// True
    miniature: bool,
    /// Innermost hull radius in points, or the default.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(hull_radius=4.2).hull_radius
    /// 4.2
    hull_radius: Option<f64>,
    /// Additional hull radius per nesting level, in points.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(hull_step=2.6).hull_step
    /// 2.6
    hull_step: Option<f64>,
    /// Line weight multiplier of recoloured miniatures, or the default.
    ///
    /// Examples
    /// --------
    /// >>> hep.DiagramStyle(miniature=True, stamp_weight=1.6).stamp_weight
    /// 1.6
    stamp_weight: Option<f64>,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramStyle {
    /// Select physics presentation without changing graph layout or drawing settings.
    ///
    /// Examples
    /// --------
    /// >>> style = hep.DiagramStyle(show_particle=False, show_momentum=True)
    ///
    /// Parameters
    /// ----------
    /// show_particle : bool or None, optional
    ///     Show particle labels; None preserves the renderer default.
    /// show_momentum : bool or None, optional
    ///     Show momentum labels; None follows the render request.
    /// show_edge_index : bool or None, optional
    ///     Show the physical edge indices.
    /// show_node_index : bool or None, optional
    ///     Show the physical vertex indices.
    /// momentum_arrows : bool or None, optional
    ///     Draw momentum arrows beside particle lines.
    /// split_initial_state : bool or None, optional
    ///     Separate the incoming states of a cross-section drawing.
    /// debug : bool or None, optional
    ///     Show both vertex and edge indices unless individually overridden.
    /// miniature : bool, optional
    ///     Default to recolouring and increase ghost visibility for small stamps.
    /// hull_radius : float or None, optional
    ///     Innermost hull radius in points, in [0.001, 1000]. Defaults to 4.2.
    ///     Effective radii including nesting and counterterm padding must
    ///     also be in this range. Thin hulls on long regions may require a larger radius.
    /// hull_step : float or None, optional
    ///     Finite nonnegative radius increment per inclusion level. Defaults to 2.6 points.
    /// stamp_weight : float or None, optional
    ///     Positive line weight multiplier of recoloured miniatures. Defaults to
    ///     1.35, which keeps lines legible at about 28 points.
    #[new]
    #[pyo3(signature=(*,show_particle=None,show_momentum=None,show_edge_index=None,show_node_index=None,momentum_arrows=None,split_initial_state=None,debug=None,miniature=false,hull_radius=None,hull_step=None,stamp_weight=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        show_particle: Option<bool>,
        show_momentum: Option<bool>,
        show_edge_index: Option<bool>,
        show_node_index: Option<bool>,
        momentum_arrows: Option<bool>,
        split_initial_state: Option<bool>,
        debug: Option<bool>,
        miniature: bool,
        hull_radius: Option<f64>,
        hull_step: Option<f64>,
        stamp_weight: Option<f64>,
    ) -> PyResult<Self> {
        if stamp_weight.is_some_and(|weight| !(weight.is_finite() && weight > 0.0)) {
            return Err(PyValueError::new_err(
                "stamp_weight must be positive and finite",
            ));
        }
        Ok(Self {
            show_particle,
            show_momentum,
            show_edge_index,
            show_node_index,
            momentum_arrows,
            split_initial_state,
            debug,
            miniature,
            hull_radius,
            hull_step,
            stamp_weight,
        })
    }
    /// Summarize the explicitly selected physics options.
    ///
    /// Examples
    /// --------
    /// >>> repr(hep.DiagramStyle())
    /// 'DiagramStyle()'
    fn __repr__(&self) -> String {
        let mut fields = [
            ("show_particle", self.show_particle),
            ("show_momentum", self.show_momentum),
            ("show_edge_index", self.show_edge_index),
            ("show_node_index", self.show_node_index),
            ("momentum_arrows", self.momentum_arrows),
            ("split_initial_state", self.split_initial_state),
            ("debug", self.debug),
        ]
        .into_iter()
        .filter_map(|(name, value)| {
            value.map(|value| format!("{name}={}", if value { "True" } else { "False" }))
        })
        .collect::<Vec<_>>();
        if self.miniature {
            fields.push("miniature=True".into());
        }
        for (name, value) in [
            ("hull_radius", self.hull_radius),
            ("hull_step", self.hull_step),
            ("stamp_weight", self.stamp_weight),
        ] {
            if let Some(value) = value {
                fields.push(format!("{name}={value:?}"));
            }
        }
        format!("DiagramStyle({})", fields.join(", "))
    }
}
impl PyDiagramStyle {
    pub(crate) fn resolve(settings: Option<&Self>, momenta: bool) -> SceneOptions {
        let empty = Self::default();
        let s = settings.unwrap_or(&empty);
        let defaults = SceneOptions::default();
        SceneOptions {
            momentum_arrows: s.momentum_arrows.unwrap_or(momenta),
            show_momentum: s.show_momentum,
            show_particle: s.show_particle.unwrap_or(defaults.show_particle),
            show_edge_index: s
                .show_edge_index
                .or(s.debug)
                .unwrap_or(defaults.show_edge_index),
            show_node_index: s
                .show_node_index
                .or(s.debug)
                .unwrap_or(defaults.show_node_index),
            split_initial_state: s
                .split_initial_state
                .unwrap_or(defaults.split_initial_state),
            overlap: if s.miniature {
                HighlightOverlap::Recolour
            } else {
                defaults.overlap
            },
            hulls: HullOptions {
                radius: s.hull_radius.unwrap_or(defaults.hulls.radius),
                step: s.hull_step.unwrap_or(defaults.hulls.step),
            },
            stamp_weight: s.stamp_weight.unwrap_or(if s.miniature {
                MINIATURE_STAMP_WEIGHT
            } else {
                defaults.stamp_weight
            }),
            ..defaults
        }
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDiagramStyle>()
}
