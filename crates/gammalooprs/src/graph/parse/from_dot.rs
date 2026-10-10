//! Import DOT graphs through FeynKit, then apply GammaLoop's own DOT attributes.

use std::{collections::HashSet, fs, path::Path, sync::Arc};

use color_eyre::Result;
use eyre::{Context, eyre};
use feynkit_amplitude::ColorRepresentation;
use feynkit_generator::FeynmanDiagramRulesExt;
use feynkit_graph::{
    DiagramEndpoint, DiagramHalfEdge, DotAttributes, DotDiagram, FeynmanDiagram, symbols,
};
use linnet::half_edge::{NodeIndex, involution::EdgeIndex};
use symbolica::atom::Atom;

use crate::{
    graph::{Graph, GroupId, attribute_warnings::warn_about_unknown_attributes},
    model::Model,
    processes::AmplitudeGraph,
};

use super::{ParseGraph, StripParse};

impl Graph {
    /// Import every graph of a DOT document.
    ///
    /// FeynKit reads compact physics DOT and its own annotated DOT, matches the
    /// model's particles and interactions, and fills in every numerator fragment
    /// and projector the input omits from the model's Feynman rules. GammaLoop then
    /// applies the attributes FeynKit does not interpret: `group_id`,
    /// `is_group_master` and `params` on the graph, `dod` on vertices, and `name`,
    /// `mass`, `dod`, `momtrop_edge_power` and `vakint_edge_power` on edges. An edge's
    /// `mass` replaces the particle mass of its propagator denominator only; a
    /// numerator the model fills in keeps the model's mass.
    ///
    /// The product of all numerator fragments, the projector and the factors must
    /// contract every tensor index but the colour indices of external legs, which the
    /// model's projector leaves open as in generated amplitudes. Explicit `num` and
    /// `projector` attributes therefore name the indices of the model's Feynman rules:
    /// the k-th index of half-edge h is `hedge(h,k)`, counted from one.
    pub fn from_dot(input: &str, model: &Model) -> Result<Vec<Self>> {
        Self::from_dot_with(input, &Arc::new(model.clone()))
    }

    /// Import the graphs of a DOT file, or of every `.dot` file in a directory in
    /// file-name order.
    pub fn from_dot_path(path: impl AsRef<Path>, model: &Model) -> Result<Vec<Self>> {
        let path = path.as_ref();
        // All graphs share one model, whose fingerprint is then computed once.
        let model = Arc::new(model.clone());
        if !path.is_dir() {
            return Self::from_dot_file(path, &model);
        }
        let mut files = fs::read_dir(path)
            .with_context(|| format!("Could not read graph directory '{}'.", path.display()))?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<std::io::Result<Vec<_>>>()?;
        files.retain(|file| file.is_file() && file.extension().is_some_and(|ext| ext == "dot"));
        files.sort();
        if files.is_empty() {
            return Err(eyre!(
                "No .dot files found in directory: {}",
                path.display()
            ));
        }
        files.iter().try_fold(Vec::new(), |mut graphs, file| {
            graphs.extend(Self::from_dot_file(file, &model)?);
            Ok(graphs)
        })
    }

    fn from_dot_file(path: &Path, model: &Arc<Model>) -> Result<Vec<Self>> {
        let input = fs::read_to_string(path)
            .with_context(|| format!("Could not read graph file '{}'.", path.display()))?;
        Self::from_dot_with(&input, model)
            .with_context(|| format!("Could not import graphs from '{}'.", path.display()))
    }

    fn from_dot_with(input: &str, model: &Arc<Model>) -> Result<Vec<Self>> {
        DotDiagram::from_dot_set(Arc::clone(model), input)?
            .iter()
            .map(Self::from_dot_diagram)
            .collect()
    }

    fn from_dot_diagram(dot: &DotDiagram) -> Result<Self> {
        warn_about_unknown_attributes(dot);
        let diagram = dot.apply_feynman_rules()?;
        // Uncontracted indices leave no scalar to evaluate; integrand generation would
        // only find out after expanding the whole network.
        Self::require_contracted_indices(
            &(diagram.numerator()
                * diagram.numerator_prefactor()
                * diagram.projector()
                * diagram.overall_factor()),
            &Self::external_colour_indices(&diagram)?,
        )
        .with_context(|| {
            format!(
                "Graph '{}' does not contract every tensor index, so it has no scalar value. \
                 Explicit `num` and `projector` attributes must name the indices of the \
                 model's Feynman rules, hedge(h,k) for the k-th index of half-edge h, \
                 counted from one",
                diagram.name().unwrap_or("<unnamed>")
            )
        })?;
        let mut parsed = ParseGraph::from_feynkit_diagram(&diagram)?;
        parsed
            .apply_dot_attributes(&dot.attributes)
            .with_context(|| {
                format!(
                    "Invalid GammaLoop attribute in graph '{}'",
                    diagram.name().unwrap_or("<unnamed>")
                )
            })?;
        Self::from_feynkit_parse_graph(&diagram, parsed)
    }

    /// The colour indices of the external legs of `diagram`, in either orientation.
    fn external_colour_indices(diagram: &FeynmanDiagram) -> Result<HashSet<Atom>> {
        let mut indices = HashSet::new();
        for (edge, _, data) in diagram
            .edges()
            .filter(|(_, _, data)| data.external.is_some())
        {
            let particle = diagram.model().particle_by_id(data.particle)?;
            let Some(colour) = ColorRepresentation::from_ufo(particle.color) else {
                continue;
            };
            for endpoint in [DiagramEndpoint::Source, DiagramEndpoint::Target] {
                let Some(hedge) = diagram.half_edge_id(DiagramHalfEdge { edge, endpoint }) else {
                    continue;
                };
                let index = symbols::hedge_index().call((hedge.0, 1));
                for representation in [colour, colour.dual()] {
                    indices.insert(representation.representation().to_symbolic([index.clone()]));
                }
            }
        }
        Ok(indices)
    }
}

impl ParseGraph {
    /// Apply GammaLoop's own DOT attributes to the graph converted from FeynKit.
    fn apply_dot_attributes(&mut self, attributes: &DotAttributes) -> Result<()> {
        let global = &mut self.global_data;
        if let Some(group_id) = attributes.graph.get("group_id") {
            global.group_id = Some(GroupId(group_id.strip_parse().context("group_id")?));
        }
        if let Some(master) = attributes.graph.get("is_group_master") {
            global.is_group_master = master.strip_parse().context("is_group_master")?;
        }
        if let Some(params) = attributes.graph.get("params") {
            global.parameters = params
                .split(';')
                .map(str::trim)
                .filter(|param| !param.is_empty())
                .map(|param| {
                    param
                        .strip_parse::<Atom>()
                        .with_context(|| format!("params entry {param}"))
                })
                .collect::<Result<_>>()?;
        }
        for (vertex, attributes) in attributes.vertices.iter().enumerate() {
            if let Some(dod) = attributes.get("dod") {
                self.graph[NodeIndex(vertex)].dod = Some(dod.strip_parse().context("vertex dod")?);
            }
        }
        for (edge, attributes) in attributes.edges.iter().enumerate() {
            let data = &mut self.graph[EdgeIndex(edge)];
            if let Some(name) = attributes.get("name") {
                data.name = Some(name.clone());
            }
            if let Some(mass) = attributes.get("mass") {
                let mass = mass.strip_parse::<Atom>().context("mass")?;
                data.particle = data.particle.clone().override_mass(mass);
            }
            if let Some(dod) = attributes.get("dod") {
                data.dod = Some(dod.strip_parse().context("edge dod")?);
            }
            if let Some(power) = attributes.get("momtrop_edge_power") {
                data.momtrop_edge_power = Some(power.strip_parse().context("momtrop_edge_power")?);
            }
            if let Some(power) = attributes.get("vakint_edge_power") {
                data.vakint_edge_power = Some(
                    power
                        .parse()
                        .map_err(|_| eyre!("invalid vakint_edge_power '{power}'"))?,
                );
            }
        }
        Ok(())
    }
}

/// Build GammaLoop graphs from DOT text, as in tests and inline fixtures.
pub trait IntoGraph<T> {
    fn into_graph(self, model: &Model) -> Result<T>;
}

impl IntoGraph<Vec<Graph>> for &str {
    fn into_graph(self, model: &Model) -> Result<Vec<Graph>> {
        Graph::from_dot(self, model)
    }
}

/// The single graph of a DOT document, whose full numerator must contract every index.
impl IntoGraph<Graph> for &str {
    fn into_graph(self, model: &Model) -> Result<Graph> {
        let [graph] = <[Graph; 1]>::try_from(Graph::from_dot(self, model)?)
            .map_err(|graphs| eyre!("expected one DOT graph, found {}", graphs.len()))?;
        graph
            .validate_full_numerator_tensor_network()
            .with_context(|| {
                format!(
                    "Failed to validate full numerator tensor network for graph {}",
                    graph.name
                )
            })?;
        Ok(graph)
    }
}

impl IntoGraph<AmplitudeGraph> for &str {
    fn into_graph(self, model: &Model) -> Result<AmplitudeGraph> {
        IntoGraph::<Graph>::into_graph(self, model).map(AmplitudeGraph::new)
    }
}

impl<T> IntoGraph<T> for String
where
    for<'a> &'a str: IntoGraph<T>,
{
    fn into_graph(self, model: &Model) -> Result<T> {
        self.as_str().into_graph(model)
    }
}

/// Build GammaLoop graphs from DOT tokens, against the Standard Model unless a
/// model name or reference follows: `dot!(digraph g { ... }, "scalars")`.
#[macro_export]
macro_rules! dot {
    // ------------------ Internal Rules (Do not call directly) ------------------

    (@internal [$($code:tt)*], $model:literal) => {

        stringify!($($code)*).into_graph(&$crate::utils::load_generic_model($model))
    };

    // Internal rule: End of parsing, with an optional argument.
    // This is matched when the accumulator has collected the code block and we hit a comma.
    (@internal [$($code:tt)*], $model:expr) => {
       stringify!($($code)*).into_graph($model)
    };

    // Internal rule: End of parsing, no optional argument.
    // This is matched when the accumulator has run out of tokens to process.
    (@internal [$($code:tt)*]) => {
        stringify!($($code)*).into_graph(&$crate::utils::load_generic_model("sm"))
    };

    // Internal rule: The "accumulator".
    // It takes the next token ($next), adds it to the $code accumulator,
    // and recursively calls the macro with the rest of the tokens ($($rest)*).
    (@internal [$($code:tt)*] $next:tt $($rest:tt)*) => {
        dot!(@internal [$($code)* $next] $($rest)*)
    };

    // ------------------ Public Entry Point ------------------

    // This is the only rule users should call. It kicks off the process by
    // invoking the internal rules with an empty accumulator `[]`.
    ($($all_tokens:tt)+) => {
        dot!(@internal [] $($all_tokens)+)
    };
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::involution::EdgeIndex;
    use symbolica::atom::{Atom, AtomCore};

    use super::IntoGraph;
    use crate::{
        graph::Graph, initialisation::test_initialise, momentum::sample::ExternalIndex,
        utils::load_generic_model,
    };

    const EE_MUMU: &str = r#"digraph ee_mumu {
        ext [style=invis];
        ext -> a [id=0, particle="e-"];
        ext -> a [id=1, particle="e+"];
        a -> b [id=2, particle="a"];
        b -> ext [id=3, particle="mu-"];
        b -> ext [id=4, particle="mu+"];
    }"#;

    #[test]
    fn plain_dot_takes_numerators_and_projector_from_the_model() {
        test_initialise().unwrap();
        let model = load_generic_model("sm");
        // The single-graph conversion also requires a fully contracted numerator.
        let graph: Graph = EE_MUMU.into_graph(&model).unwrap();
        for (_, _, vertex) in graph.underlying.iter_nodes() {
            assert!(vertex.vertex_rule.is_some());
            assert!(!vertex.num.autogenerated && vertex.num.value != Atom::one());
        }
        let photon = &graph.underlying[EdgeIndex(2)];
        assert_ne!(photon.num.value, Atom::one());
        let projector = graph.global_prefactor.projector.to_canonical_string();
        for spinor in ["::u(", "::vbar(", "::ubar(", "::v("] {
            assert!(projector.contains(spinor), "{projector}");
        }

        // Explicit factors override the model's, in the index identities of its rules.
        let (_, _, vertex) = graph.underlying.iter_nodes().next().unwrap();
        let explicit = |statements: &str| {
            Graph::from_dot(
                &EE_MUMU.replace(
                    "ext [style=invis];",
                    &format!("{statements} ext [style=invis];"),
                ),
                &model,
            )
        };
        let overridden = explicit(&format!(
            r#"projector="{projector}"; a [num="3*{}"];"#,
            vertex.num.value.to_canonical_string()
        ))
        .unwrap()
        .remove(0);
        assert_eq!(
            overridden.global_prefactor.projector,
            graph.global_prefactor.projector
        );
        let (_, _, overridden_vertex) = overridden.underlying.iter_nodes().next().unwrap();
        assert_eq!(overridden_vertex.num.value, &vertex.num.value * 3);
        // `hedge(h)` is no index of the model's rules: the spinors stay uncontracted.
        let Err(error) = explicit(&format!(
            r#"projector="{}";"#,
            projector.replace(",1)", ")")
        )) else {
            panic!("an uncontracted projector must not import");
        };
        assert!(format!("{error:?}").contains("hedge(h,k)"), "{error:?}");
    }

    #[test]
    fn amplitudes_may_leave_only_external_colour_open() {
        test_initialise().unwrap();
        let model = load_generic_model("sm");
        let ddx_aa = |statements: &str| {
            format!(
                r#"digraph ddx_aa {{
                    {statements}
                    ext [style=invis];
                    ext -> a:0 [particle="d"];
                    ext -> b:1 [particle="d~"];
                    a -> b [particle="d"];
                    a:2 -> ext [particle="a"];
                    b:3 -> ext [particle="a"];
                }}"#
            )
        };
        // As in a generated amplitude, the quark colours stay open.
        Graph::from_dot(&ddx_aa(""), &model).unwrap();
        let average = |slot: &str| {
            format!(
                r#"num="spenso::g(spenso::dind(spenso::cof(3,hedge(0{slot}))),spenso::cof(3,hedge(1{slot})))/3";"#
            )
        };
        Graph::from_dot(&ddx_aa(&average(",1")), &model).unwrap();
        // `hedge(h)` names no index of the model's rules.
        let Err(error) = Graph::from_dot(&ddx_aa(&average("")), &model) else {
            panic!("a colour average over hedge(h) indices must not import");
        };
        assert!(format!("{error:?}").contains("hedge(h,k)"), "{error:?}");
    }

    #[test]
    fn amplitude_external_momenta_follow_the_dangling_half_edges() {
        test_initialise().unwrap();
        let model = load_generic_model("sm");
        let graph: Graph = r#"digraph amplitude_order {
                    ext [style=invis];
                    ext -> v:3 [id=0, particle="H"];
                    v:2 -> ext [id=1, particle="H"];
                    v:1 -> ext [id=2, particle="H"];
                    ext -> v:0 [id=3, particle="H"];
                    v [num="1"];
                }"#
        .into_graph(&model)
        .unwrap();
        let order = (0..4)
            .map(|index| graph.loop_momentum_basis.ext_edges[ExternalIndex(index)].0)
            .collect::<Vec<_>>();
        assert_eq!(order, [3, 2, 1, 0]);
    }

    #[test]
    fn cross_section_external_momenta_follow_the_cut_tags() {
        test_initialise().unwrap();
        let model = load_generic_model("scalars");
        let graph: Graph = r#"digraph forward {
                    node [num=1];
                    ext [style=invis];
                    ext -> a [particle="scalar_2", is_cut=5];
                    ext -> a [particle="scalar_0", is_cut=2];
                    a -> b [particle="scalar_1"];
                    a -> b [particle="scalar_1"];
                    b -> ext [particle="scalar_2", is_cut=5];
                    b -> ext [particle="scalar_0", is_cut=2];
                }"#
        .into_graph(&model)
        .unwrap();
        let particles = (0..2)
            .map(|index| {
                let edge = graph.loop_momentum_basis.ext_edges[ExternalIndex(index)];
                model
                    .particle_by_id(graph.underlying[edge].particle())
                    .unwrap()
                    .name
                    .clone()
            })
            .collect::<Vec<_>>();
        assert_eq!(particles, ["scalar_0", "scalar_2"]);
        assert!(!graph.finalized_cuts.is_empty());
    }
}
