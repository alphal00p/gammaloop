//! UV spinneys, woods, forests and counterterms of FeynKit diagrams.
//!
//! The combinatorics and drawings live in `feynkit-uv`; these classes keep the
//! diagram they belong to and convert errors at the Python boundary.

use std::collections::BTreeSet;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::sync::{Arc, OnceLock};

use feynkit_graph::{
    FeynmanDiagram, HighlightOverlap, LoopMomentumBasis, MINIATURE_STAMP_WEIGHT, SceneOptions,
};
use feynkit_uv::{
    Counterterm, CountertermScheme, Expansion, ForestTerm, Forests,
    RenormalizationPrescriptionSettings, STAMP_SIZE, Spinney, SpinneyRelation, StampOptions,
    UnknownCountertermScheme, UvDisplayError, UvGraph, UvRegion, Wood,
};
use linnest::svg::Config;
use linnet::half_edge::NodeIndex;
use linnet::half_edge::subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike};
use linnet_render_py::{
    PyDiagramRender, PyDiagramView, PyRenderSettings, highlight_overlap, typst_export,
};
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::{
    Bound, Py, PyAny, PyAnyMethods, PyModuleMethods, PyRef, PyResult, Python, pyclass, pymethods,
};
use pyo3::types::{PyDict, PyDictMethods, PyList, PyModule};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::display::{self, escape_html};
use crate::error::{self, UvError};
use crate::graph::{PyFeynmanDiagram, PySubgraph};
use crate::render_settings::PyDiagramStyle;

type Basis = LoopMomentumBasis;

/// Woods and their lazily unfolded forests, shared by every object taken from them.
struct WoodState {
    diagram: PyFeynmanDiagram,
    wood: Wood<Basis>,
    /// Nonempty spinneys in the wood's deterministic order.
    nodes: Vec<NodeIndex>,
    forests: OnceLock<Forests<Basis>>,
    expansion: OnceLock<Expansion>,
}

impl WoodState {
    fn new(diagram: PyFeynmanDiagram, wood: Wood<Basis>) -> Arc<Self> {
        let nodes = wood
            .nodes_by_edges(&diagram.inner)
            .into_iter()
            .filter(|node| !wood.spinney(*node).is_empty())
            .collect();
        Arc::new(Self {
            diagram,
            wood,
            nodes,
            forests: OnceLock::new(),
            expansion: OnceLock::new(),
        })
    }

    fn forests(&self) -> &Forests<Basis> {
        self.forests.get_or_init(|| self.wood.clone().unfold())
    }

    fn expansion(&self) -> &Expansion {
        self.expansion
            .get_or_init(|| self.forests().expansion(&self.diagram.inner))
    }

    fn spinney(&self, spinney: &Spinney<Basis>) -> PySpinney {
        PySpinney {
            diagram: self.diagram.clone(),
            inner: spinney.clone(),
        }
    }
}

/// Stamp size and line weight from the Python keywords.
fn stamps(size: f64, weight: f64) -> StampOptions {
    StampOptions { size, weight }
}

/// The drawing configuration of every stamp.
fn stamp_config(config: Option<&PyRenderSettings>) -> PyResult<Config> {
    config.map_or_else(|| Ok(Config::default()), PyRenderSettings::config)
}

/// A static view of a complete stamp-graph SVG with a caption.
fn stamp_graph_view(svg: String, caption: &str) -> PyDiagramView {
    let svg = display::themed_stamp_graph_svg(svg);
    let html = format!(
        "<figure class=\"feynkit-diagram\" style=\"max-width:100%;margin:0\">\
         <div style=\"max-width:100%;overflow-x:auto\">{svg}</div>\
         <figcaption style=\"font-size:.85em;opacity:.75;text-align:center;margin-top:.2rem\">{}</figcaption></figure>",
        escape_html(caption)
    );
    PyDiagramView::from_render(PyDiagramRender::new(svg, html))
}

/// Draw `regions` on `diagram` with the usual physics options.
fn region_view(
    py: Python<'_>,
    diagram: &PyFeynmanDiagram,
    regions: &[UvRegion],
    config: Option<&PyRenderSettings>,
    style: Option<&PyDiagramStyle>,
    overlap: Option<&str>,
) -> PyResult<PyDiagramView> {
    let mut options: SceneOptions = PyDiagramStyle::resolve(style, false);
    if let Some(overlap) = overlap {
        options.overlap = highlight_overlap(overlap)?;
    } else if options.overlap == HighlightOverlap::Auto {
        // A single ordinary layer would keep plain strokes, which carry no
        // hover text; hulls show every spinney's details.
        options.overlap = HighlightOverlap::Hulls;
    }
    options.hulls.validate().map_err(PyValueError::new_err)?;
    diagram.draw_layers(py, config, &options, None, &UvRegion::layers(regions))
}

/// An explicit scheme, or `None` to derive it from the approximations.
fn parse_scheme(scheme: Option<&str>) -> PyResult<Option<CountertermScheme>> {
    scheme
        .map(|name| {
            name.parse()
                .map_err(|error: UnknownCountertermScheme| PyValueError::new_err(error.to_string()))
        })
        .transpose()
}

/// One UV subgraph (spinney) of a Feynman diagram: a union of cycles of its
/// internal lines with a nonnegative degree of divergence.
///
/// Power counting is in four dimensions and uses the diagram's numerator as
/// stored: apply Feynman rules first (``apply_feynman_rules``) for diagrams read
/// from DOT with non-scalar particles. A numerator that vanishes only after index
/// contraction is not detected.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> sunrise = hep.FeynmanDiagram.from_dot(hep.Model.phi4(), '''digraph sunrise {
/// ...     edge [particle=phi]; ext [style=invis];
/// ...     ext -> a; a -> b; a -> b; a -> b; b -> ext; }''')
/// >>> full = hep.Spinney(sunrise.subgraph(edges=[0, 1, 2]))
/// >>> (full.loop_count, full.degree, full.superficial_degree)
/// (2, 2, 2)
/// >>> bubble = hep.Spinney(sunrise.subgraph(edges=[0, 1]))
/// >>> bubble.relation(full)
/// 'contained in'
/// >>> view = full.render()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Spinney",
    module = "symbolica.community.hepkit",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PySpinney {
    diagram: PyFeynmanDiagram,
    inner: Spinney<Basis>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PySpinney {
    /// Classify the internal lines of ``subgraph`` as a spinney.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> spinney = hep.Spinney(sunrise.subgraph(edges=[0, 2]))
    ///
    /// Parameters
    /// ----------
    /// subgraph : Subgraph
    ///     A selection of the diagram's internal lines, each with both halves
    ///     selected; external legs are dropped.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If the selection holds only one half of an internal line, is not a
    ///     nonempty union of cycles, has a superficially convergent connected
    ///     component, or its loop momenta cannot be routed.
    #[new]
    fn new(subgraph: PyRef<'_, PySubgraph>) -> PyResult<Self> {
        let diagram = subgraph.into_super();
        let whole = diagram.whole();
        let inner = whole.inner.as_ref();
        let graph = inner.hedge_graph();
        let lines = inner.uv_domain(&diagram.selection());
        if let Some(half) = lines
            .included_iter()
            .find(|hedge| !lines.includes(&graph.inv(*hedge)))
        {
            return Err(UvError::new_err(format!(
                "the subgraph selects only one half of internal edge {}; select whole lines",
                graph[&half].0
            )));
        }
        let selection = InternalSubGraph::cleaned_filter_optimist(lines, graph);
        if selection.is_empty() || !graph.bridges_of(&selection.filter).is_empty() {
            return Err(UvError::new_err(
                "the subgraph is not a nonempty union of cycles of internal lines",
            ));
        }
        let convergent = || {
            UvError::new_err(
                "the subgraph or one of its connected components is superficially convergent",
            )
        };
        let spinney = Spinney::new(selection, inner, inner.loop_momentum_basis())
            .map_err(error::uv)?
            .ok_or_else(convergent)?;
        if !spinney.has_divergent_factors(inner).map_err(error::uv)? {
            return Err(convergent());
        }
        Ok(Self {
            diagram: whole,
            inner: spinney,
        })
    }

    /// The diagram this spinney belongs to.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> assert full.diagram.id == sunrise.id
    #[getter]
    fn diagram(&self) -> PyFeynmanDiagram {
        self.diagram.clone()
    }

    /// The spinney as a selection of its diagram.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> region = full.subgraph
    #[getter]
    fn subgraph(&self, py: Python<'_>) -> PyResult<Py<PySubgraph>> {
        self.diagram.subgraph_view(py, self.inner.filter().clone())
    }

    /// The internal edge IDs, in increasing order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> full.edges
    /// [0, 1, 2]
    #[getter]
    fn edges(&self) -> Vec<usize> {
        self.inner
            .edges(self.diagram.inner.as_ref())
            .into_iter()
            .map(|edge| edge.0)
            .collect()
    }

    /// The connected components, as selections of the diagram.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> assert len(full.components) == 1
    #[getter]
    fn components(&self, py: Python<'_>) -> PyResult<Vec<Py<PySubgraph>>> {
        self.inner
            .components()
            .iter()
            .map(|component| self.diagram.subgraph_view(py, component.clone()))
            .collect()
    }

    /// The number of independent loops.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> bubble.loop_count
    /// 1
    #[getter]
    fn loop_count(&self) -> usize {
        self.inner.loop_count()
    }

    /// The classifying degree of divergence in four dimensions.
    ///
    /// Loop momenta of the spinney are scaled with its external momenta fixed;
    /// numerator factors carrying only external momenta do not count. Zero is
    /// logarithmic, positive is power divergent.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> bubble.degree
    /// 0
    #[getter]
    fn degree(&self) -> i32 {
        self.inner.dod()
    }

    /// The local superficial degree of divergence in four dimensions.
    ///
    /// Enumeration prunes with the larger of this and ``degree``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> bubble.superficial_degree
    /// 0
    ///
    /// Raises
    /// ------
    /// DiagramError
    ///     If the diagram cannot power count the spinney.
    #[getter]
    fn superficial_degree(&self) -> PyResult<i32> {
        self.inner
            .superficial_degree(self.diagram.inner.as_ref())
            .map_err(error::diagram)
    }

    /// Whether the spinney is UV divergent: every connected component has a
    /// nonnegative classifying degree. Spinneys are only built divergent, so
    /// this holds unless the diagram's power counting changed.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> assert full.is_divergent
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If the diagram cannot power count a component.
    #[getter]
    fn is_divergent(&self) -> PyResult<bool> {
        self.inner
            .has_divergent_factors(self.diagram.inner.as_ref())
            .map_err(error::uv)
    }

    /// The subtraction scheme the renormalization prescription selected:
    /// ``"MUV"``, ``"PolePart"``, ``"OS"`` or ``"IR"``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> full.approximation
    /// 'MUV'
    #[getter]
    fn approximation(&self) -> String {
        self.inner.renormalization_scheme().to_string()
    }

    /// A short name listing the internal edges, as shown on hover.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> full.label
    /// 'γ{0,1,2}'
    #[getter]
    fn label(&self) -> String {
        self.inner.label(self.diagram.inner.as_ref())
    }

    /// How this spinney relates to ``other`` of the same diagram: ``"equal"``,
    /// ``"contains"``, ``"contained in"``, ``"disjoint"`` (no common vertex),
    /// ``"touching"`` (a common vertex but no common line) or ``"overlapping"``
    /// (a common line).
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> full.relation(bubble)
    /// 'contains'
    ///
    /// Parameters
    /// ----------
    /// other : Spinney
    ///     A spinney of the same diagram.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``other`` belongs to another diagram.
    fn relation(&self, other: PyRef<'_, Self>) -> PyResult<String> {
        Ok(self.checked_relation(&other)?.to_string())
    }

    /// Whether this spinney and ``other`` can belong to one forest: nested or
    /// disjoint (no common vertex), not touching or overlapping.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> assert bubble.is_compatible(full)
    ///
    /// Parameters
    /// ----------
    /// other : Spinney
    ///     A spinney of the same diagram.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``other`` belongs to another diagram.
    fn is_compatible(&self, other: PyRef<'_, Self>) -> PyResult<bool> {
        Ok(self.checked_relation(&other)?.is_compatible())
    }

    /// This spinney contracted to a point: a counterterm graph.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> counterterm = bubble.counterterm("MS-bar")
    ///
    /// Parameters
    /// ----------
    /// scheme : {"MS-bar", "on-shell", "unspecified"} or None, optional
    ///     Case-insensitive renormalization scheme. None derives it from the
    ///     approximations: on-shell when every contracted spinney is subtracted
    ///     on shell (``"OS"``), unspecified otherwise.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the scheme is unknown.
    #[pyo3(signature = (scheme=None))]
    fn counterterm(&self, scheme: Option<&str>) -> PyResult<PyCounterterm> {
        PyCounterterm::build(
            self.diagram.clone(),
            [self.inner.clone()],
            parse_scheme(scheme)?,
        )
    }

    /// Draw the diagram with this spinney highlighted, as a hull by default.
    /// Hover shows its edges, loop count, degrees and scheme.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Spinney`` class example:
    ///
    /// >>> view = full.render(style=hep.DiagramStyle(show_particle=False))
    ///
    /// Parameters
    /// ----------
    /// config : RenderSettings or None, optional
    ///     Layout and drawing settings, as in ``FeynmanDiagram.render``.
    /// style : DiagramStyle or None, optional
    ///     Particle and index presentation; ``miniature=True`` draws a stamp.
    /// overlap : {"hulls", "recolour", "casings", "bands", "bands-miniature", "lanes"} or None, optional
    ///     Highlight geometry, as in ``FeynmanDiagram.render``. Defaults to
    ///     hulls, which carry the hover text, or recolouring for miniatures.
    #[pyo3(signature = (*, config=None, style=None, overlap=None))]
    fn render(
        &self,
        py: Python<'_>,
        config: Option<&PyRenderSettings>,
        style: Option<&PyDiagramStyle>,
        #[gen_stub(override_type(type_repr = "typing.Literal['hulls', 'recolour', 'casings', 'bands', 'bands-miniature', 'lanes'] | None", imports = ("typing")))]
        overlap: Option<&str>,
    ) -> PyResult<PyDiagramView> {
        let regions = [self.inner.region(self.diagram.inner.as_ref())];
        region_view(py, &self.diagram, &regions, config, style, overlap)
    }

    /// Display the highlighted diagram in Jupyter and IPython.
    ///
    /// Examples
    /// --------
    /// >>> from IPython.display import display
    /// >>> display(full)
    ///
    /// Parameters
    /// ----------
    /// include : list[str] or None, optional
    ///     MIME types requested by the frontend; None allows all available types.
    /// exclude : list[str] or None, optional
    ///     MIME types omitted from the frontend response.
    #[pyo3(signature = (include=None, exclude=None))]
    fn _repr_mimebundle_(
        &self,
        py: Python<'_>,
        include: Option<Vec<String>>,
        exclude: Option<Vec<String>>,
    ) -> PyResult<Py<PyAny>> {
        self.render(py, None, None, None)?
            ._repr_mimebundle_(py, include, exclude)
    }

    /// Return Marimo MIME data for the highlighted diagram.
    ///
    /// Examples
    /// --------
    /// >>> data = full._mime_()
    fn _mime_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.render(py, None, None, None)?._mime_(py)
    }

    /// Whether both select the same lines of the same diagram.
    ///
    /// Examples
    /// --------
    /// >>> assert full == hep.Spinney(sunrise.subgraph(edges=[0, 1, 2]))
    ///
    /// Parameters
    /// ----------
    /// other : Any
    ///     The object to compare with.
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other.extract::<PyRef<'_, Self>>().is_ok_and(|other| {
            Arc::ptr_eq(&self.diagram.inner, &other.diagram.inner) && self.inner == other.inner
        })
    }

    /// Hash consistent with equality.
    ///
    /// Examples
    /// --------
    /// >>> spinneys = {full, bubble}
    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    /// A one-line summary.
    ///
    /// Examples
    /// --------
    /// >>> repr(bubble)
    /// 'Spinney(γ{0,1}, loops=1, degree=0, MUV)'
    fn __repr__(&self) -> String {
        format!(
            "Spinney({}, loops={}, degree={}, {})",
            self.label(),
            self.inner.loop_count(),
            self.inner.dod(),
            self.inner.renormalization_scheme()
        )
    }
}

impl PySpinney {
    fn checked_relation(&self, other: &Self) -> PyResult<SpinneyRelation> {
        if !Arc::ptr_eq(&self.diagram.inner, &other.diagram.inner) {
            return Err(PyValueError::new_err(
                "spinneys belong to different diagrams",
            ));
        }
        Ok(self
            .inner
            .relation(&other.inner, self.diagram.inner.as_ref()))
    }
}

/// The UV-divergent spinneys of a diagram, ordered by inclusion.
///
/// ``Wood(diagram)`` enumerates the divergent unions of cycles of the internal
/// lines with FeynKit's pruning and default renormalization prescription.
/// Iteration yields spinneys by loop count, then by edge IDs; the empty root of the
/// inclusion order is implicit. Display draws the Hasse diagram of stamps.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> basketball = hep.FeynmanDiagram.from_dot(
/// ...     hep.Model.phi4(), "digraph { edge [particle=phi]; a -> b; a -> b; a -> b; a -> b; }")
/// >>> wood = hep.Wood(basketball)
/// >>> len(wood)
/// 11
/// >>> [spinney.loop_count for spinney in wood].count(1)
/// 6
/// >>> forests = wood.forests()
/// >>> len(forests)
/// 46
/// >>> view = wood.operation_graph()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Wood",
    module = "symbolica.community.hepkit",
    frozen,
    skip_from_py_object
)]
pub struct PyWood {
    state: Arc<WoodState>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyWood {
    /// Enumerate and classify the divergent spinneys of ``diagram``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> wood = hep.Wood(basketball)
    ///
    /// Parameters
    /// ----------
    /// diagram : FeynmanDiagram
    ///     A complete diagram, not a subgraph view.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If power counting, particle lookup or routing fails.
    #[new]
    fn new(py: Python<'_>, diagram: PyRef<'_, PyFeynmanDiagram>) -> PyResult<Self> {
        let diagram = diagram.clone();
        let inner = Arc::clone(diagram.as_shared_diagram()?);
        // Box the error so the detached closure returns a small result.
        let wood = py
            .detach(|| {
                let spinneys = inner
                    .classify_spinneys(
                        &inner.hedge_graph().full_filter(),
                        inner.model(),
                        &RenormalizationPrescriptionSettings::default(),
                        inner.loop_momentum_basis(),
                    )
                    .map_err(Box::new)?;
                Wood::from_spinneys(spinneys, inner.as_ref(), inner.loop_momentum_basis())
                    .map_err(Box::new)
            })
            .map_err(|error| error::uv(*error))?;
        Ok(Self {
            state: WoodState::new(diagram, wood),
        })
    }

    /// Build the wood of user-given spinneys of ``diagram``.
    ///
    /// Each spinney is checked against ``diagram``: its selection width, its
    /// lines, components and loop count, and the loop-momentum routing the
    /// diagram gives it. A spinney classified on a structurally identical diagram
    /// is accepted; one routed differently (for example with a reversed line) or
    /// taken from another topology is rejected. Each spinney is power counted
    /// in ``diagram``, so one convergent here (for example after changing the
    /// numerators) is rejected as ``Spinney(...)`` rejects it. Products of
    /// vertex-disjoint spinneys, connected or not, are added as
    /// ``Wood(diagram)`` has them, so their product forests are kept.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> pair = hep.Spinney(basketball.subgraph(edges=[0, 1]))
    /// >>> triple = hep.Spinney(basketball.subgraph(edges=[0, 1, 2]))
    /// >>> chain = hep.Wood.from_spinneys(basketball, [pair, triple])
    /// >>> len(chain.forests())
    /// 4
    ///
    /// Parameters
    /// ----------
    /// diagram : FeynmanDiagram
    ///     The complete diagram the spinneys must belong to.
    /// spinneys : Sequence[Spinney]
    ///     The spinneys; equal selections are merged.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If a spinney does not belong to ``diagram``; the message names the
    ///     failed check ("its width differs", "its routed sub-basis differs", ...).
    ///     Also if a spinney or one of its components is convergent in ``diagram``.
    #[staticmethod]
    fn from_spinneys(
        diagram: PyRef<'_, PyFeynmanDiagram>,
        spinneys: Vec<PyRef<'_, PySpinney>>,
    ) -> PyResult<Self> {
        let diagram = diagram.clone();
        let inner = Arc::clone(diagram.as_shared_diagram()?);
        let (graph, basis) = (inner.as_ref(), inner.loop_momentum_basis());
        // Equal selections count once, so repeated inputs cost nothing extra.
        let mut seen = BTreeSet::new();
        let given: Vec<_> = spinneys
            .iter()
            .filter(|spinney| seen.insert(spinney.inner.filter().clone()))
            .map(|spinney| spinney.inner.clone())
            .collect();
        // Validate the given spinneys before forming their products.
        Wood::from_spinneys(given.iter().cloned(), graph, basis).map_err(error::uv)?;
        // Power count each spinney in this diagram, as Spinney(...) does: a
        // spinney taken from a diagram with other numerators keeps neither its
        // degree nor its divergence here.
        let mut given = given
            .into_iter()
            .map(|spinney| {
                let label = spinney.label(graph);
                let convergent = || {
                    UvError::new_err(format!(
                        "{label} or one of its connected components is superficially \
                         convergent in this diagram"
                    ))
                };
                let dod = graph
                    .degree_of_divergence(spinney.filter())
                    .map_err(|error| UvError::new_err(error.to_string()))?;
                let rebuilt = Spinney::with_scheme(
                    spinney.subgraph().clone(),
                    graph,
                    basis,
                    spinney.renormalization_scheme(),
                    dod,
                )
                .map_err(error::uv)?
                .ok_or_else(convergent)?;
                if !rebuilt.has_divergent_factors(graph).map_err(error::uv)? {
                    return Err(convergent());
                }
                Ok(rebuilt)
            })
            .collect::<PyResult<Vec<_>>>()?;
        let products = Wood::disjoint_products(&given, graph, basis).map_err(error::uv)?;
        given.extend(products);
        let wood = Wood::from_spinneys(given, graph, basis).map_err(error::uv)?;
        Ok(Self {
            state: WoodState::new(diagram, wood),
        })
    }

    /// The diagram of this wood.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> assert wood.diagram.id == basketball.id
    #[getter]
    fn diagram(&self) -> PyFeynmanDiagram {
        self.state.diagram.clone()
    }

    /// The number of nonempty spinneys.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> len(wood)
    /// 11
    fn __len__(&self) -> usize {
        self.state.nodes.len()
    }

    /// The spinney at ``index`` in iteration order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> largest = wood[-1]
    ///
    /// Parameters
    /// ----------
    /// index : int
    ///     Position, negative values counting from the end.
    ///
    /// Raises
    /// ------
    /// IndexError
    ///     If ``index`` is out of range.
    fn __getitem__(&self, index: isize) -> PyResult<PySpinney> {
        let length = self.state.nodes.len() as isize;
        let position = if index < 0 { index + length } else { index };
        if !(0..length).contains(&position) {
            return Err(PyIndexError::new_err("wood index out of range"));
        }
        let node = self.state.nodes[position as usize];
        Ok(self.state.spinney(self.state.wood.spinney(node)))
    }

    /// Iterate over the spinneys by loop count, then by lines.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> loops = [spinney.loop_count for spinney in wood]
    fn __iter__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let spinneys = self.spinneys();
        Ok(PyList::new(py, spinneys)?
            .call_method0("__iter__")?
            .unbind())
    }

    /// The spinneys, in iteration order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> spinneys = wood.spinneys()
    fn spinneys(&self) -> Vec<PySpinney> {
        self.state
            .nodes
            .iter()
            .map(|node| self.state.spinney(self.state.wood.spinney(*node)))
            .collect()
    }

    /// The cover relation of inclusion between nonempty spinneys: pairs
    /// ``(smaller, larger)`` with no spinney strictly between them.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> len(wood.covers())
    /// 16
    fn covers(&self) -> Vec<(PySpinney, PySpinney)> {
        let wood = &self.state.wood;
        wood.cover_relation()
            .into_iter()
            .filter(|(smaller, _)| !wood.spinney(*smaller).is_empty())
            .map(|(smaller, larger)| {
                (
                    self.state.spinney(wood.spinney(smaller)),
                    self.state.spinney(wood.spinney(larger)),
                )
            })
            .collect()
    }

    /// The forests of the forest formula in ``k`` order, as in ``k_table``: the
    /// empty (unsubtracted) forest first, then by depth.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> forests = wood.forests()
    /// >>> assert forests[0].is_empty
    /// >>> assert [forest.k for forest in forests] == list(range(len(forests)))
    fn forests(&self, py: Python<'_>) -> Vec<PyForest> {
        let state = &self.state;
        py.detach(|| state.expansion());
        state
            .expansion()
            .terms()
            .iter()
            .map(|term| PyForest {
                state: Arc::clone(state),
                node: term.forest,
            })
            .collect()
    }

    /// Draw the Hasse diagram: one stamp per spinney, with the largest on top and
    /// the bare diagram at the bottom. Hover shows each spinney's details.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> svg = wood.render(stamp_weight=1.6).to_svg()
    ///
    /// Parameters
    /// ----------
    /// config : RenderSettings or None, optional
    ///     Layout settings applied to every stamp.
    /// stamp_size : float, optional
    ///     Larger side of the bare stamp, in points. Defaults to 28.
    /// stamp_weight : float, optional
    ///     Line weight multiplier of the stamps. Defaults to 1.35.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If drawing fails or a stamp option is not positive.
    #[pyo3(signature = (*, config=None, stamp_size=STAMP_SIZE, stamp_weight=MINIATURE_STAMP_WEIGHT))]
    fn render(
        &self,
        py: Python<'_>,
        config: Option<&PyRenderSettings>,
        stamp_size: f64,
        stamp_weight: f64,
    ) -> PyResult<PyDiagramView> {
        let config = stamp_config(config)?;
        let state = &self.state;
        let svg = py
            .detach(|| {
                state.wood.hasse_svg(
                    &state.diagram.inner,
                    &config,
                    stamps(stamp_size, stamp_weight),
                )
            })
            .map_err(error::uv_display)?;
        Ok(stamp_graph_view(
            svg,
            &format!(
                "Wood of {}: {} spinney{}",
                state.diagram.inner.name(),
                state.nodes.len(),
                if state.nodes.len() == 1 { "" } else { "s" }
            ),
        ))
    }

    /// Draw the forest operation graph in the style of Linnest's
    /// ``tree-auto-root``. Each box holds its forest's term ``Op[m(S|G) · inner]``:
    /// the stamp shows the approximated spinney ``S`` in blue over the treated
    /// spinneys ``G`` in red, ``inner`` is the nested term for ``G`` and the empty
    /// forest is 𝟙. Boxes are ranked left to right by the size of their cover;
    /// the root is yellow and forests covering every line of every spinney are
    /// light blue. Hover shows each term as ``K_{S|G}``, its stored quantity ``k``
    /// and its sign.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> svg = wood.operation_graph(op="T").to_svg()
    ///
    /// Parameters
    /// ----------
    /// op : str, optional
    ///     Operator symbol, such as ``"K"`` (the default), ``"T"`` or
    ///     ``"T^{3D}"``; ``^{…}`` and ``_{…}`` set scripts.
    /// config : RenderSettings or None, optional
    ///     Layout settings applied to every stamp.
    /// stamp_size : float, optional
    ///     Larger side of the bare stamp, in points. Defaults to 28.
    /// stamp_weight : float, optional
    ///     Line weight multiplier of the stamps. Defaults to 1.35.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If drawing fails or a stamp option is not positive.
    #[pyo3(signature = (*, op="K", config=None, stamp_size=STAMP_SIZE, stamp_weight=MINIATURE_STAMP_WEIGHT))]
    fn operation_graph(
        &self,
        py: Python<'_>,
        op: &str,
        config: Option<&PyRenderSettings>,
        stamp_size: f64,
        stamp_weight: f64,
    ) -> PyResult<PyDiagramView> {
        let svg = self.forest_svg(py, config, |forests, diagram, expansion, config| {
            forests.operation_graph_svg(
                diagram,
                expansion,
                config,
                stamps(stamp_size, stamp_weight),
                op,
            )
        })?;
        let count = self.state.forests().n_forests();
        Ok(stamp_graph_view(
            svg,
            &format!(
                "Forests of {}: {count} forest{}",
                self.state.diagram.inner.name(),
                if count == 1 { "" } else { "s" }
            ),
        ))
    }

    /// Draw the forest formula term by term, one row per depth:
    /// ``m(Γ|∅) − m(Γ|γ)·Op[m(γ|∅)] − ⋯``. Each term carries its cofactor
    /// ``m(Γ|γ)``, the whole diagram in blue with the term's outermost spinneys
    /// in red. Hover shows each term as ``K_{S|G}`` and its stored quantity ``k``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> view = wood.expansion()
    /// >>> typst_source = wood.expansion_typst(op="T")
    ///
    /// Parameters
    /// ----------
    /// op : str, optional
    ///     Operator symbol, as in ``operation_graph``.
    /// config : RenderSettings or None, optional
    ///     Layout settings applied to every stamp.
    /// stamp_size : float, optional
    ///     Larger side of the bare stamp, in points. Defaults to 28.
    /// stamp_weight : float, optional
    ///     Line weight multiplier of the stamps. Defaults to 1.35.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If drawing fails or a stamp option is not positive.
    #[pyo3(signature = (*, op="K", config=None, stamp_size=STAMP_SIZE, stamp_weight=MINIATURE_STAMP_WEIGHT))]
    fn expansion(
        &self,
        py: Python<'_>,
        op: &str,
        config: Option<&PyRenderSettings>,
        stamp_size: f64,
        stamp_weight: f64,
    ) -> PyResult<PyDiagramView> {
        let svg = self.forest_svg(py, config, |forests, diagram, expansion, config| {
            forests.expansion_svg(
                diagram,
                expansion,
                config,
                stamps(stamp_size, stamp_weight),
                op,
            )
        })?;
        Ok(stamp_graph_view(
            svg,
            &format!("Forest formula of {}", self.state.diagram.inner.name()),
        ))
    }

    /// The forest formula of ``expansion`` as a self-contained Typst document
    /// embedding the same SVG.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> from pathlib import Path
    /// >>> Path("forests.typ").write_text(wood.expansion_typst())
    ///
    /// Parameters
    /// ----------
    /// op : str, optional
    ///     Operator symbol, as in ``operation_graph``.
    /// config : RenderSettings or None, optional
    ///     Layout settings applied to every stamp.
    /// stamp_size : float, optional
    ///     Larger side of the bare stamp, in points. Defaults to 28.
    /// stamp_weight : float, optional
    ///     Line weight multiplier of the stamps. Defaults to 1.35.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If drawing fails or a stamp option is not positive.
    #[pyo3(signature = (*, op="K", config=None, stamp_size=STAMP_SIZE, stamp_weight=MINIATURE_STAMP_WEIGHT))]
    fn expansion_typst(
        &self,
        py: Python<'_>,
        op: &str,
        config: Option<&PyRenderSettings>,
        stamp_size: f64,
        stamp_weight: f64,
    ) -> PyResult<String> {
        let svg = self.forest_svg(py, config, |forests, diagram, expansion, config| {
            forests.expansion_svg(
                diagram,
                expansion,
                config,
                stamps(stamp_size, stamp_weight),
                op,
            )
        })?;
        Ok(typst_export::svg_source(&display::themed_stamp_graph_svg(
            svg,
        )))
    }

    /// The stored quantities of the forest formula, one line per term in
    /// topological order: ``k0 = K[∅]``, ``k5 = K_{Γ|a}[k2]`` where the bracket
    /// names the nested term, and ``k3 = k1·k2`` for a product of disjoint
    /// factors. Spinneys are named ``a``, ``b``, … in iteration order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> print(wood.k_table())
    ///
    /// Parameters
    /// ----------
    /// op : str, optional
    ///     Operator symbol, ``"K"`` by default.
    #[pyo3(signature = (op="K"))]
    fn k_table(&self, py: Python<'_>, op: &str) -> String {
        let state = &self.state;
        py.detach(|| state.expansion().k_table(op).join("\n"))
    }

    /// The name of ``spinney`` in the notation of ``k_table``, such as ``"a"``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> wood.name(wood[0])
    /// 'a'
    ///
    /// Parameters
    /// ----------
    /// spinney : Spinney
    ///     A spinney of this wood.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``spinney`` is not a spinney of this wood.
    fn name(&self, spinney: PyRef<'_, PySpinney>) -> PyResult<String> {
        let state = &self.state;
        Arc::ptr_eq(&state.diagram.inner, &spinney.diagram.inner)
            .then(|| {
                state
                    .wood
                    .spinneys()
                    .find(|(_, candidate)| candidate.filter() == spinney.inner.filter())
                    .and_then(|(node, _)| state.expansion().name(node))
                    .map(str::to_owned)
            })
            .flatten()
            .ok_or_else(|| PyValueError::new_err("the spinney is not in this wood"))
    }

    /// A copy of this wood whose spinneys carry the names in ``names`` in the
    /// expansion, ``k_table``, ``repr`` and hover text; the others keep theirs.
    /// By default the whole graph is ``Γ`` and the other spinneys are ``a``,
    /// ``b``, … in iteration order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Wood`` class example:
    ///
    /// >>> named = wood.renamed({wood[0]: "x", "Γ": "G"})
    /// >>> named.name(wood[0])
    /// 'x'
    ///
    /// Parameters
    /// ----------
    /// names : dict
    ///     Maps a spinney of this wood, or its current name, to its new name: a
    ///     nonempty string without whitespace, control characters, ``|``, braces,
    ///     brackets, commas, ``·`` or ``∅``, not of the form ``k<digits>``, and
    ///     ``Γ`` only for the whole graph.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If a key is no spinney or name of this wood, a spinney is keyed twice,
    ///     a name is invalid, or two spinneys would share a name.
    fn renamed(&self, py: Python<'_>, names: &Bound<'_, PyDict>) -> PyResult<Self> {
        let state = &self.state;
        let expansion = py.detach(|| state.expansion().clone());
        let mut renames = Vec::new();
        for (key, name) in names.iter() {
            let node = if let Ok(spinney) = key.extract::<PyRef<'_, PySpinney>>() {
                Arc::ptr_eq(&state.diagram.inner, &spinney.diagram.inner)
                    .then(|| {
                        state
                            .wood
                            .spinneys()
                            .find(|(_, candidate)| candidate.filter() == spinney.inner.filter())
                            .map(|(node, _)| node)
                    })
                    .flatten()
                    .filter(|node| expansion.name(*node).is_some())
                    .ok_or_else(|| {
                        PyValueError::new_err(format!(
                            "{} is not a spinney of this wood",
                            spinney.inner.label(&*spinney.diagram.inner)
                        ))
                    })?
            } else if let Ok(current) = key.extract::<String>() {
                state
                    .nodes
                    .iter()
                    .copied()
                    .find(|node| expansion.name(*node) == Some(current.as_str()))
                    .ok_or_else(|| {
                        PyValueError::new_err(format!(
                            "no spinney of this wood is named {current:?}"
                        ))
                    })?
            } else {
                return Err(PyTypeError::new_err(
                    "names maps spinneys or their current names to new names",
                ));
            };
            if renames.iter().any(|(other, _)| *other == node) {
                return Err(PyValueError::new_err(format!(
                    "spinney {} is given more than one new name",
                    expansion.name(node).unwrap_or("?")
                )));
            }
            renames.push((node, name.extract::<String>()?));
        }
        let mut expansion = expansion;
        expansion
            .rename(renames)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let renamed = WoodState::new(state.diagram.clone(), state.wood.clone());
        let _ = renamed.expansion.set(expansion);
        Ok(Self { state: renamed })
    }

    /// Display the Hasse diagram in Jupyter and IPython.
    ///
    /// Examples
    /// --------
    /// >>> from IPython.display import display
    /// >>> display(wood)
    ///
    /// Parameters
    /// ----------
    /// include : list[str] or None, optional
    ///     MIME types requested by the frontend; None allows all available types.
    /// exclude : list[str] or None, optional
    ///     MIME types omitted from the frontend response.
    #[pyo3(signature = (include=None, exclude=None))]
    fn _repr_mimebundle_(
        &self,
        py: Python<'_>,
        include: Option<Vec<String>>,
        exclude: Option<Vec<String>>,
    ) -> PyResult<Py<PyAny>> {
        self.render(py, None, STAMP_SIZE, MINIATURE_STAMP_WEIGHT)?
            ._repr_mimebundle_(py, include, exclude)
    }

    /// Return Marimo MIME data for the Hasse diagram.
    ///
    /// Examples
    /// --------
    /// >>> data = wood._mime_()
    fn _mime_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.render(py, None, STAMP_SIZE, MINIATURE_STAMP_WEIGHT)?
            ._mime_(py)
    }

    /// A one-line summary.
    ///
    /// Examples
    /// --------
    /// >>> repr(wood)
    /// 'Wood(11 spinneys, max loops 3)'
    fn __repr__(&self) -> String {
        format!(
            "Wood({} spinney{}, max loops {})",
            self.state.nodes.len(),
            if self.state.nodes.len() == 1 { "" } else { "s" },
            self.state.wood.max_loop_count()
        )
    }
}

impl PyWood {
    /// Draw from the unfolded forests without the GIL.
    fn forest_svg(
        &self,
        py: Python<'_>,
        config: Option<&PyRenderSettings>,
        draw: impl FnOnce(
            &Forests<Basis>,
            &FeynmanDiagram,
            &Expansion,
            &Config,
        ) -> Result<String, UvDisplayError>
        + Send,
    ) -> PyResult<String> {
        let config = stamp_config(config)?;
        let state = &self.state;
        py.detach(|| {
            draw(
                state.forests(),
                &state.diagram.inner,
                state.expansion(),
                &config,
            )
        })
        .map_err(error::uv_display)
    }
}

/// One forest of the forest formula: pairwise nested or vertex-disjoint spinneys
/// of a wood, applied as a trace of subtractions in Foata normal form.
///
/// Each Foata level holds commuting subtractions; later levels act on the result
/// of earlier ones. A product of disjoint spinneys lists each factor as its own
/// element. The sign is ``(-1)**depth``. Display stacks the spinneys as
/// hulls on the diagram: blue outside, red innermost, violet in between.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> eye = hep.FeynmanDiagram.from_dot(hep.Model.phi4(), '''digraph eye {
/// ...     edge [particle=phi]; ext [style=invis]; ext -> a; ext -> a;
/// ...     a -> b; a -> c; b -> c; b -> c; b -> ext; c -> ext; }''')
/// >>> forests = hep.Wood(eye).forests()
/// >>> nested = forests[-1]
/// >>> (nested.depth, nested.sign)
/// (2, 1)
/// >>> [[spinney.loop_count for spinney in level] for level in nested.levels]
/// [[1], [2]]
/// >>> view = nested.render()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Forest",
    module = "symbolica.community.hepkit",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PyForest {
    state: Arc<WoodState>,
    node: NodeIndex,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyForest {
    /// The forest's spinneys, in application order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> [spinney.loop_count for spinney in nested.spinneys]
    /// [1, 2]
    #[getter]
    fn spinneys(&self) -> Vec<PySpinney> {
        self.forests()
            .forest_spinneys(self.node)
            .into_iter()
            .map(|spinney| self.state.spinney(spinney))
            .collect()
    }

    /// The Foata levels, first applied first, each listing its spinneys.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> levels = nested.levels
    #[getter]
    fn levels(&self) -> Vec<Vec<PySpinney>> {
        self.forests()
            .forest_levels(self.node)
            .into_iter()
            .map(|level| {
                level
                    .into_iter()
                    .map(|spinney| self.state.spinney(spinney))
                    .collect()
            })
            .collect()
    }

    /// The number of subtractions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> nested.depth
    /// 2
    #[getter]
    fn depth(&self) -> usize {
        self.forests().operation(self.node).depth()
    }

    /// The sign of the forest term, ``(-1)**depth``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> forests[1].sign
    /// -1
    #[getter]
    fn sign(&self) -> i8 {
        self.forests().operation(self.node).sign()
    }

    /// Whether this is the empty, unsubtracted forest.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> assert forests[0].is_empty
    #[getter]
    fn is_empty(&self) -> bool {
        self.forests().operation(self.node).key.is_empty()
    }

    /// Whether the forest is the product of forests of disjoint components of a
    /// disconnected spinney.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> assert not nested.is_disjoint_union
    #[getter]
    fn is_disjoint_union(&self) -> bool {
        self.forests().is_disjoint_union(self.node)
    }

    /// The union of the forest's spinneys, or None for the empty forest.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> nested.cover.loop_count
    /// 2
    #[getter]
    fn cover(&self) -> Option<PySpinney> {
        let forests = self.forests();
        (!forests.operation(self.node).key.is_empty())
            .then(|| self.state.spinney(forests.source_spinney(self.node)))
    }

    /// The forests one subtraction smaller that this forest extends.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> [parent.depth for parent in nested.parents]
    /// [1]
    #[getter]
    fn parents(&self) -> Vec<PyForest> {
        self.neighbours(true)
    }

    /// The forests that extend this forest by one subtraction.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> children = forests[0].children
    #[getter]
    fn children(&self) -> Vec<PyForest> {
        self.neighbours(false)
    }

    /// The forest's spinneys contracted: a counterterm graph.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> counterterm = nested.counterterm("on-shell")
    ///
    /// Parameters
    /// ----------
    /// scheme : {"MS-bar", "on-shell", "unspecified"} or None, optional
    ///     Case-insensitive renormalization scheme. None derives it from the
    ///     approximations: on-shell when every contracted spinney is subtracted
    ///     on shell (``"OS"``), unspecified otherwise.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     For the empty forest.
    /// ValueError
    ///     If the scheme is unknown.
    #[pyo3(signature = (scheme=None))]
    fn counterterm(&self, scheme: Option<&str>) -> PyResult<PyCounterterm> {
        PyCounterterm::build(
            self.state.diagram.clone(),
            self.forests()
                .forest_spinneys(self.node)
                .into_iter()
                .cloned(),
            parse_scheme(scheme)?,
        )
    }

    /// Draw the diagram with the forest's spinneys stacked as regions: hulls at
    /// reading size, recoloured with ``DiagramStyle(miniature=True)``. Hover
    /// shows each spinney's details.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> view = nested.render(overlap="lanes")
    ///
    /// Parameters
    /// ----------
    /// config : RenderSettings or None, optional
    ///     Layout and drawing settings, as in ``FeynmanDiagram.render``.
    /// style : DiagramStyle or None, optional
    ///     Particle and index presentation.
    /// overlap : {"hulls", "recolour", "casings", "bands", "bands-miniature", "lanes"} or None, optional
    ///     Highlight geometry; defaults to hulls, which carry the hover text,
    ///     or recolouring for miniatures.
    #[pyo3(signature = (*, config=None, style=None, overlap=None))]
    fn render(
        &self,
        py: Python<'_>,
        config: Option<&PyRenderSettings>,
        style: Option<&PyDiagramStyle>,
        #[gen_stub(override_type(type_repr = "typing.Literal['hulls', 'recolour', 'casings', 'bands', 'bands-miniature', 'lanes'] | None", imports = ("typing")))]
        overlap: Option<&str>,
    ) -> PyResult<PyDiagramView> {
        let regions = self.forests().regions(self.node, &self.state.diagram.inner);
        region_view(py, &self.state.diagram, &regions, config, style, overlap)
    }

    /// Display the forest in Jupyter and IPython.
    ///
    /// Examples
    /// --------
    /// >>> from IPython.display import display
    /// >>> display(nested)
    ///
    /// Parameters
    /// ----------
    /// include : list[str] or None, optional
    ///     MIME types requested by the frontend; None allows all available types.
    /// exclude : list[str] or None, optional
    ///     MIME types omitted from the frontend response.
    #[pyo3(signature = (include=None, exclude=None))]
    fn _repr_mimebundle_(
        &self,
        py: Python<'_>,
        include: Option<Vec<String>>,
        exclude: Option<Vec<String>>,
    ) -> PyResult<Py<PyAny>> {
        self.render(py, None, None, None)?
            ._repr_mimebundle_(py, include, exclude)
    }

    /// Return Marimo MIME data for the forest.
    ///
    /// Examples
    /// --------
    /// >>> data = nested._mime_()
    fn _mime_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.render(py, None, None, None)?._mime_(py)
    }

    /// Whether both are the same forest of the same wood.
    ///
    /// Examples
    /// --------
    /// >>> assert nested == forests[-1]
    ///
    /// Parameters
    /// ----------
    /// other : Any
    ///     The object to compare with.
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other
            .extract::<PyRef<'_, Self>>()
            .is_ok_and(|other| Arc::ptr_eq(&self.state, &other.state) && self.node == other.node)
    }

    /// Hash consistent with equality.
    ///
    /// Examples
    /// --------
    /// >>> unique = set(forests)
    fn __hash__(&self) -> u64 {
        self.node.0 as u64
    }

    /// The forest's term in text, such as ``K_{Γ|a}``: each outermost spinney
    /// given the spinneys nested directly in it, named as in ``Wood.k_table``.
    /// The empty forest is ``K[∅]``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> nested.notation()
    /// 'K_{Γ|a}'
    ///
    /// Parameters
    /// ----------
    /// op : str, optional
    ///     Operator symbol, ``"K"`` by default.
    #[pyo3(signature = (op="K"))]
    fn notation(&self, op: &str) -> String {
        let expansion = self.state.expansion();
        expansion.notation(self.forest_term(), op)
    }

    /// The index of the forest's stored quantity ``k``, in the topological order
    /// of ``Wood.k_table``; the empty forest is ``k0``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> nested.k
    /// 3
    #[getter]
    fn k(&self) -> usize {
        self.forest_term().k
    }

    /// Draw the forest's term of the forest formula: its sign, the cofactor
    /// ``m(Γ|γ)`` with the outermost spinneys in red, and the factors
    /// ``Op[m(S|G) · inner]``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Forest`` class example:
    ///
    /// >>> view = nested.term(op="T")
    ///
    /// Parameters
    /// ----------
    /// op : str, optional
    ///     Operator symbol, as in ``Wood.operation_graph``.
    /// config : RenderSettings or None, optional
    ///     Layout settings applied to every stamp.
    /// stamp_size : float, optional
    ///     Larger side of the bare stamp, in points. Defaults to 28.
    /// stamp_weight : float, optional
    ///     Line weight multiplier of the stamps. Defaults to 1.35.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If drawing fails or a stamp option is not positive.
    #[pyo3(signature = (*, op="K", config=None, stamp_size=STAMP_SIZE, stamp_weight=MINIATURE_STAMP_WEIGHT))]
    fn term(
        &self,
        py: Python<'_>,
        op: &str,
        config: Option<&PyRenderSettings>,
        stamp_size: f64,
        stamp_weight: f64,
    ) -> PyResult<PyDiagramView> {
        let config = stamp_config(config)?;
        let state = &self.state;
        let node = self.node;
        let svg = py
            .detach(|| {
                state.forests().term_svg(
                    node,
                    &state.diagram.inner,
                    state.expansion(),
                    &config,
                    stamps(stamp_size, stamp_weight),
                    op,
                )
            })
            .map_err(error::uv_display)?;
        Ok(stamp_graph_view(
            svg,
            &format!(
                "{} = k{}",
                self.state.expansion().notation(self.forest_term(), op),
                self.forest_term().k
            ),
        ))
    }

    /// A one-line summary: the term, its ``k`` and its sign.
    ///
    /// Examples
    /// --------
    /// >>> repr(nested)
    /// 'Forest(K_{Γ|a}, k3, depth=2, sign=+1)'
    fn __repr__(&self) -> String {
        let forest = self.forests().operation(self.node);
        format!(
            "Forest({}, k{}, depth={}, sign={:+})",
            self.state.expansion().notation(self.forest_term(), "K"),
            self.forest_term().k,
            forest.depth(),
            forest.sign()
        )
    }
}

impl PyForest {
    fn forests(&self) -> &Forests<Basis> {
        self.state.forests()
    }

    fn forest_term(&self) -> &ForestTerm {
        self.state
            .expansion()
            .term(self.node)
            .expect("every forest of a wood has a term")
    }

    fn neighbours(&self, parents: bool) -> Vec<PyForest> {
        self.forests()
            .edges()
            .filter_map(|edge| {
                let (from, to) = if parents {
                    (edge.child, edge.parent)
                } else {
                    (edge.parent, edge.child)
                };
                (from == self.node).then_some(to)
            })
            .map(|node| PyForest {
                state: Arc::clone(&self.state),
                node,
            })
            .collect()
    }
}

/// A counterterm graph: a diagram with a forest of spinneys contracted, in a
/// renormalization scheme.
///
/// This is the combinatorial object only; counterterm expressions and
/// integrands are not built yet. A disconnected spinney contracts as its
/// connected components, as its product forest does. Display keeps the
/// diagram's layout and fades the contracted spinneys inside hulls, one colour
/// per nesting level. Each outermost spinney carries one circled cross: open
/// for MS-bar, filled on shell and in a dashed ring while the scheme is
/// unspecified. Its hover lists the contracted chain, inner to outer, with the
/// scheme. ``scheme`` is positional-or-keyword here, in ``with_scheme`` and in
/// ``Spinney.counterterm`` and ``Forest.counterterm``.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> sunrise = hep.FeynmanDiagram.from_dot(hep.Model.phi4(), '''digraph sunrise {
/// ...     edge [particle=phi]; ext [style=invis];
/// ...     ext -> a; a -> b; a -> b; a -> b; b -> ext; }''')
/// >>> bubble = hep.Spinney(sunrise.subgraph(edges=[0, 1]))
/// >>> counterterm = hep.Counterterm([bubble], scheme="on-shell")
/// >>> counterterm.scheme
/// 'on-shell'
/// >>> view = counterterm.render()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Counterterm",
    module = "symbolica.community.hepkit",
    frozen,
    skip_from_py_object
)]
pub struct PyCounterterm {
    diagram: PyFeynmanDiagram,
    inner: Counterterm<Basis>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyCounterterm {
    /// Contract a spinney, a forest, or pairwise nested or disjoint spinneys.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Counterterm`` class example:
    ///
    /// >>> local = hep.Counterterm(bubble, scheme="MS-bar")
    ///
    /// Parameters
    /// ----------
    /// contracted : Spinney or Forest or Sequence[Spinney]
    ///     The spinneys to contract, all of one diagram, pairwise nested or
    ///     disjoint (no common vertex). A disconnected spinney stands for its
    ///     connected components.
    /// scheme : {"MS-bar", "on-shell", "unspecified"} or None, optional
    ///     Case-insensitive renormalization scheme. None derives it from the
    ///     approximations: on-shell when every contracted spinney is subtracted
    ///     on shell (``"OS"``), unspecified otherwise.
    ///
    /// Raises
    /// ------
    /// UvError
    ///     If nothing nonempty is contracted, two spinneys overlap or touch, or a
    ///     spinney is selected twice (also as a component of a disconnected one).
    /// ValueError
    ///     If the spinneys belong to different diagrams or the scheme is unknown.
    /// TypeError
    ///     If ``contracted`` has another type.
    #[new]
    #[pyo3(signature = (contracted, scheme=None))]
    fn new(
        #[gen_stub(override_type(type_repr = "Spinney | Forest | typing.Sequence[Spinney]", imports = ("typing")))]
        contracted: &Bound<'_, PyAny>,
        scheme: Option<&str>,
    ) -> PyResult<Self> {
        if let Ok(spinney) = contracted.extract::<PyRef<'_, PySpinney>>() {
            return spinney.counterterm(scheme);
        }
        if let Ok(forest) = contracted.extract::<PyRef<'_, PyForest>>() {
            return forest.counterterm(scheme);
        }
        let scheme = parse_scheme(scheme)?;
        let spinneys: Vec<PyRef<'_, PySpinney>> = contracted.extract().map_err(|_| {
            PyTypeError::new_err("contract a Spinney, a Forest or a sequence of spinneys")
        })?;
        let first = spinneys
            .first()
            .ok_or_else(|| UvError::new_err("a counterterm contracts at least one spinney"))?;
        if spinneys
            .iter()
            .any(|spinney| !Arc::ptr_eq(&spinney.diagram.inner, &first.diagram.inner))
        {
            return Err(PyValueError::new_err(
                "spinneys belong to different diagrams",
            ));
        }
        Self::build(
            first.diagram.clone(),
            spinneys.iter().map(|spinney| spinney.inner.clone()),
            scheme,
        )
    }

    /// The diagram whose spinneys are contracted.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Counterterm`` class example:
    ///
    /// >>> assert counterterm.diagram.id == sunrise.id
    #[getter]
    fn diagram(&self) -> PyFeynmanDiagram {
        self.diagram.clone()
    }

    /// The contracted spinneys, outermost first.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Counterterm`` class example:
    ///
    /// >>> [spinney.label for spinney in counterterm.contracted]
    /// ['γ{0,1}']
    #[getter]
    fn contracted(&self) -> Vec<PySpinney> {
        self.inner
            .contracted()
            .iter()
            .map(|spinney| PySpinney {
                diagram: self.diagram.clone(),
                inner: spinney.clone(),
            })
            .collect()
    }

    /// The renormalization scheme: ``"MS-bar"``, ``"on-shell"`` or
    /// ``"unspecified"``.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Counterterm`` class example:
    ///
    /// >>> counterterm.scheme
    /// 'on-shell'
    #[getter]
    fn scheme(&self) -> String {
        self.inner.scheme().to_string()
    }

    /// The same contraction in another scheme.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Counterterm`` class example:
    ///
    /// >>> counterterm.with_scheme("MS-bar").scheme
    /// 'MS-bar'
    ///
    /// Parameters
    /// ----------
    /// scheme : {"MS-bar", "on-shell", "unspecified"} or None, optional
    ///     Case-insensitive renormalization scheme. None derives it from the
    ///     approximations, as the constructor does.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the scheme is unknown.
    #[pyo3(signature = (scheme=None))]
    fn with_scheme(&self, scheme: Option<&str>) -> PyResult<Self> {
        let scheme = parse_scheme(scheme)?
            .unwrap_or_else(|| CountertermScheme::implied_by(self.inner.contracted()));
        Ok(Self {
            diagram: self.diagram.clone(),
            inner: self.inner.clone().with_scheme(scheme),
        })
    }

    /// Draw the counterterm graph in place: each outermost contracted spinney
    /// fades inside a red hull marked by one circled cross, open for MS-bar,
    /// filled on shell and in a dashed ring while unspecified. Its hover lists the
    /// contracted chain, inner to outer; nested spinneys show as coloured hulls.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Counterterm`` class example:
    ///
    /// >>> stamp = counterterm.render(style=hep.DiagramStyle(miniature=True))
    ///
    /// Parameters
    /// ----------
    /// config : RenderSettings or None, optional
    ///     Layout and drawing settings, as in ``FeynmanDiagram.render``.
    /// style : DiagramStyle or None, optional
    ///     Particle and index presentation; ``miniature=True`` draws a stamp.
    #[pyo3(signature = (*, config=None, style=None))]
    fn render(
        &self,
        py: Python<'_>,
        config: Option<&PyRenderSettings>,
        style: Option<&PyDiagramStyle>,
    ) -> PyResult<PyDiagramView> {
        let regions = self.inner.regions(&self.diagram.inner);
        region_view(py, &self.diagram, &regions, config, style, None)
    }

    /// Display the counterterm graph in Jupyter and IPython.
    ///
    /// Examples
    /// --------
    /// >>> from IPython.display import display
    /// >>> display(counterterm)
    ///
    /// Parameters
    /// ----------
    /// include : list[str] or None, optional
    ///     MIME types requested by the frontend; None allows all available types.
    /// exclude : list[str] or None, optional
    ///     MIME types omitted from the frontend response.
    #[pyo3(signature = (include=None, exclude=None))]
    fn _repr_mimebundle_(
        &self,
        py: Python<'_>,
        include: Option<Vec<String>>,
        exclude: Option<Vec<String>>,
    ) -> PyResult<Py<PyAny>> {
        self.render(py, None, None)?
            ._repr_mimebundle_(py, include, exclude)
    }

    /// Return Marimo MIME data for the counterterm graph.
    ///
    /// Examples
    /// --------
    /// >>> data = counterterm._mime_()
    fn _mime_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.render(py, None, None)?._mime_(py)
    }

    /// A one-line summary.
    ///
    /// Examples
    /// --------
    /// >>> repr(counterterm)
    /// 'Counterterm(γ{0,1}, scheme=on-shell)'
    fn __repr__(&self) -> String {
        let contracted = self
            .inner
            .contracted()
            .iter()
            .map(|spinney| spinney.label(self.diagram.inner.as_ref()))
            .collect::<Vec<_>>();
        format!(
            "Counterterm({}, scheme={})",
            contracted.join(" "),
            self.inner.scheme()
        )
    }

    /// Whether both contract the same spinneys of the same diagram in the same
    /// scheme.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Counterterm`` class example:
    ///
    /// >>> assert counterterm == hep.Counterterm([bubble], "on-shell")
    ///
    /// Parameters
    /// ----------
    /// other : Any
    ///     The object to compare with.
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other.extract::<PyRef<'_, Self>>().is_ok_and(|other| {
            Arc::ptr_eq(&self.diagram.inner, &other.diagram.inner)
                && self.inner.scheme() == other.inner.scheme()
                && self.filters().eq(other.filters())
        })
    }

    /// Hash consistent with equality.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Counterterm`` class example:
    ///
    /// >>> unique = {counterterm, counterterm.with_scheme("MS-bar")}
    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        for filter in self.filters() {
            filter.hash(&mut hasher);
        }
        self.inner.scheme().to_string().hash(&mut hasher);
        hasher.finish()
    }
}

impl PyCounterterm {
    /// The contracted selections, outermost first.
    fn filters(&self) -> impl Iterator<Item = &SuBitGraph> {
        self.inner.contracted().iter().map(Spinney::filter)
    }

    fn build(
        diagram: PyFeynmanDiagram,
        spinneys: impl IntoIterator<Item = Spinney<Basis>>,
        scheme: Option<CountertermScheme>,
    ) -> PyResult<Self> {
        let spinneys: Vec<_> = spinneys.into_iter().collect();
        let scheme = scheme.unwrap_or_else(|| CountertermScheme::implied_by(&spinneys));
        let graph = diagram.inner.as_ref();
        let inner = Counterterm::new(spinneys, scheme, graph, graph.loop_momentum_basis())
            .map_err(error::counterterm)?;
        Ok(Self { diagram, inner })
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySpinney>()?;
    module.add_class::<PyWood>()?;
    module.add_class::<PyForest>()?;
    module.add_class::<PyCounterterm>()
}
