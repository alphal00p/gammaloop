use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::{self, Display},
    str::FromStr,
};

use color_eyre::Result;
use eyre::{Context, eyre};
use itertools::Itertools;
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{Cycle, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};
use typed_index_collections::TiVec;

use crate::{
    graph::{Graph, GraphGroup, GroupId, edge::EdgeMass},
    model::{ArcParticle, Model},
};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum GraphGroupSelectionMode {
    #[default]
    MasterGraphs,
    CrossSectionAmplitudeGraphs,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GraphGroupSelectionSpec {
    rules: Vec<GraphGroupSelectionRule>,
    mode: GraphGroupSelectionMode,
}

impl GraphGroupSelectionSpec {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_master_graph_names(graph_names: Vec<String>) -> Self {
        Self::new().with_master_graph_names(graph_names)
    }

    pub fn with_mode(mut self, mode: GraphGroupSelectionMode) -> Self {
        self.mode = mode;
        self
    }

    pub fn mode(&self) -> GraphGroupSelectionMode {
        self.mode
    }

    pub fn with_master_graph_names(self, graph_names: Vec<String>) -> Self {
        self.with_master_graph_names_polarity(SelectionPolarity::With, graph_names)
    }

    pub fn with_master_graph_names_polarity(
        mut self,
        polarity: SelectionPolarity,
        graph_names: Vec<String>,
    ) -> Self {
        if !graph_names.is_empty() {
            self.rules.push(GraphGroupSelectionRule::MasterGraphNames {
                polarity,
                graph_names,
            });
        }
        self
    }

    pub fn with_raised_propagator_signatures(
        mut self,
        polarity: SelectionPolarity,
        scope: RaisedPropagatorScope,
        signatures: Vec<RaisedPropagatorSignature>,
    ) -> Self {
        if !signatures.is_empty() {
            self.rules.push(GraphGroupSelectionRule::RaisedPropagators {
                polarity,
                scope,
                signatures,
            });
        }
        self
    }

    pub fn with_raised_cut_signatures(
        mut self,
        polarity: SelectionPolarity,
        scope: RaisedPropagatorScope,
        signatures: Vec<RaisedPropagatorSignature>,
    ) -> Self {
        if !signatures.is_empty() {
            self.rules.push(GraphGroupSelectionRule::RaisedCuts {
                polarity,
                scope,
                signatures,
            });
        }
        self
    }

    pub fn with_cycle_signatures(
        mut self,
        polarity: SelectionPolarity,
        signatures: Vec<CycleSignature>,
    ) -> Self {
        if !signatures.is_empty() {
            self.rules.push(GraphGroupSelectionRule::Cycles {
                polarity,
                signatures,
            });
        }
        self
    }

    pub fn with_vertex_signatures(
        mut self,
        polarity: SelectionPolarity,
        signatures: Vec<VertexSignature>,
    ) -> Self {
        if !signatures.is_empty() {
            self.rules.push(GraphGroupSelectionRule::Vertices {
                polarity,
                signatures,
            });
        }
        self
    }

    pub fn with_particle_signatures(
        mut self,
        polarity: SelectionPolarity,
        signatures: Vec<ParticleSignature>,
    ) -> Self {
        if !signatures.is_empty() {
            self.rules.push(GraphGroupSelectionRule::Particles {
                polarity,
                signatures,
            });
        }
        self
    }

    pub fn is_empty(&self) -> bool {
        self.rules.is_empty()
    }

    pub fn has_raised_cut_rules(&self) -> bool {
        self.rules
            .iter()
            .any(GraphGroupSelectionRule::uses_cut_analysis)
    }

    pub(crate) fn has_graph_analysis_rules(&self) -> bool {
        self.rules
            .iter()
            .any(GraphGroupSelectionRule::uses_graph_analysis)
    }

    pub fn plan<'a, F>(
        &self,
        graph_group_structure: &TiVec<GroupId, GraphGroup>,
        graph_by_id: F,
    ) -> Result<GraphGroupSelectionPlan>
    where
        F: FnMut(usize) -> Option<&'a Graph>,
    {
        self.plan_with_analysis_contexts(
            graph_group_structure,
            graph_by_id,
            |_master_graph_id, master_graph| {
                Ok(vec![GraphSelectionSubject::whole_graph(master_graph)])
            },
            |_master_graph_id, _master_graph| Ok(Vec::new()),
            "Graph-group selection structural filters have no graph analysis subjects.",
            "Graph-group selection raised-cut filters require cross-section Cutkosky cuts.",
        )
    }

    pub(crate) fn plan_with_analysis_contexts<'a, F, S, C>(
        &self,
        graph_group_structure: &TiVec<GroupId, GraphGroup>,
        mut graph_by_id: F,
        mut analysis_subjects_by_master_id: S,
        mut cut_subjects_by_master_id: C,
        no_structural_subjects_message: &str,
        no_cut_subjects_message: &str,
    ) -> Result<GraphGroupSelectionPlan>
    where
        F: FnMut(usize) -> Option<&'a Graph>,
        S: FnMut(usize, &'a Graph) -> Result<Vec<GraphSelectionSubject<'a>>>,
        C: FnMut(usize, &'a Graph) -> Result<Vec<GraphCutSelectionSubject<'a>>>,
    {
        if self.rules.is_empty() {
            return Err(eyre!("No graph-group selection rules were provided."));
        }
        if graph_group_structure.is_empty() {
            return Err(eyre!("Cannot select graph groups from an empty integrand."));
        }

        let candidates = graph_group_structure
            .iter_enumerated()
            .map(|(group_id, group)| {
                let master_graph_id = group.master();
                let master_graph = graph_by_id(master_graph_id).ok_or_else(|| {
                    eyre!(
                        "Graph group {} refers to missing master graph id {}.",
                        group_id.0,
                        master_graph_id
                    )
                })?;
                let graph_names = group
                    .into_iter()
                    .map(|graph_id| {
                        graph_by_id(graph_id)
                            .map(|graph| graph.name.clone())
                            .ok_or_else(|| {
                                eyre!(
                                    "Graph group {} refers to missing graph id {}.",
                                    group_id.0,
                                    graph_id
                                )
                            })
                    })
                    .collect::<Result<Vec<_>>>()?;
                let analysis_subjects =
                    analysis_subjects_by_master_id(master_graph_id, master_graph)?;
                let cut_subjects = if self.has_raised_cut_rules() {
                    cut_subjects_by_master_id(master_graph_id, master_graph)?
                } else {
                    Vec::new()
                };
                Ok(GraphGroupSelectionCandidate {
                    group_id,
                    master_graph_id,
                    master_graph_name: master_graph.name.clone(),
                    analysis_subjects,
                    cut_subjects,
                    graph_names,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        if self.has_graph_analysis_rules()
            && candidates
                .iter()
                .all(|candidate| candidate.analysis_subjects.is_empty())
        {
            return Err(eyre!(no_structural_subjects_message.to_string()));
        }
        if self.has_raised_cut_rules()
            && candidates
                .iter()
                .all(|candidate| candidate.cut_subjects.is_empty())
        {
            return Err(eyre!(no_cut_subjects_message.to_string()));
        }

        let mut master_name_to_group = BTreeMap::<String, GroupId>::new();
        for candidate in &candidates {
            if let Some(previous) =
                master_name_to_group.insert(candidate.master_graph_name.clone(), candidate.group_id)
            {
                return Err(eyre!(
                    "Master graph name '{}' is ambiguous: it appears in groups {} and {}.",
                    candidate.master_graph_name,
                    previous.0,
                    candidate.group_id.0
                ));
            }
        }

        let authoritative_group_ids =
            self.authoritative_master_graph_group_ids(&candidates, &master_name_to_group)?;
        let non_authoritative_rules = self
            .rules
            .iter()
            .filter(|rule| !rule.is_authoritative_master_graph_name_rule())
            .collect::<Vec<_>>();

        let mut retained_group_ids = if non_authoritative_rules.is_empty() {
            authoritative_group_ids.clone()
        } else {
            candidates
                .iter()
                .map(|candidate| candidate.group_id)
                .collect::<BTreeSet<_>>()
        };
        for rule in non_authoritative_rules {
            let rule_groups = rule.resolve(&candidates, &master_name_to_group)?;
            retained_group_ids = retained_group_ids
                .intersection(&rule_groups)
                .copied()
                .collect::<BTreeSet<_>>();
        }
        retained_group_ids.extend(authoritative_group_ids);

        if retained_group_ids.is_empty() {
            return Err(eyre!(
                "Graph-group selection would remove all graph groups."
            ));
        }

        let retained_group_ids = candidates
            .iter()
            .filter_map(|candidate| {
                retained_group_ids
                    .contains(&candidate.group_id)
                    .then_some(candidate.group_id)
            })
            .collect::<Vec<_>>();
        let old_to_new_group_id = retained_group_ids
            .iter()
            .enumerate()
            .map(|(new_group_id, old_group_id)| (*old_group_id, GroupId(new_group_id)))
            .collect::<BTreeMap<_, _>>();
        let kept_master_graphs = candidates
            .iter()
            .filter(|candidate| old_to_new_group_id.contains_key(&candidate.group_id))
            .map(|candidate| candidate.master_graph_name.clone())
            .collect::<Vec<_>>();
        let removed_master_graphs = candidates
            .iter()
            .filter(|candidate| !old_to_new_group_id.contains_key(&candidate.group_id))
            .map(|candidate| candidate.master_graph_name.clone())
            .collect::<Vec<_>>();
        let removed_graphs = candidates
            .iter()
            .filter(|candidate| !old_to_new_group_id.contains_key(&candidate.group_id))
            .flat_map(|candidate| candidate.graph_names.iter().cloned())
            .collect::<Vec<_>>();

        Ok(GraphGroupSelectionPlan {
            retained_group_ids,
            old_to_new_group_id,
            report: GraphGroupSelectionReport {
                kept_master_graphs,
                removed_master_graphs,
                removed_graphs,
            },
        })
    }

    fn authoritative_master_graph_group_ids<'a>(
        &self,
        candidates: &[GraphGroupSelectionCandidate<'a>],
        master_name_to_group: &BTreeMap<String, GroupId>,
    ) -> Result<BTreeSet<GroupId>> {
        let mut group_ids = BTreeSet::new();
        for rule in self
            .rules
            .iter()
            .filter(|rule| rule.is_authoritative_master_graph_name_rule())
        {
            group_ids.extend(
                rule.resolve_authoritative_master_graph_names(candidates, master_name_to_group)?,
            );
        }
        Ok(group_ids)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GraphGroupSelectionPlan {
    retained_group_ids: Vec<GroupId>,
    old_to_new_group_id: BTreeMap<GroupId, GroupId>,
    report: GraphGroupSelectionReport,
}

impl GraphGroupSelectionPlan {
    pub fn retained_group_ids(&self) -> &[GroupId] {
        &self.retained_group_ids
    }

    pub fn new_group_id_for_old(&self, group_id: GroupId) -> Option<GroupId> {
        self.old_to_new_group_id.get(&group_id).copied()
    }

    pub fn report(&self) -> &GraphGroupSelectionReport {
        &self.report
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GraphGroupSelectionReport {
    pub kept_master_graphs: Vec<String>,
    pub removed_master_graphs: Vec<String>,
    pub removed_graphs: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionPolarity {
    With,
    Without,
}

impl SelectionPolarity {
    fn keep_if_match(self, matched: bool) -> bool {
        match self {
            Self::With => matched,
            Self::Without => !matched,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RaisedPropagatorScope {
    All,
    Massive,
    Massless,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum RaisedPropagatorSignature {
    Exact(Vec<usize>),
    AnyRaising,
}

impl RaisedPropagatorSignature {
    pub const ANY_RAISING_KEYWORD: &'static str = "ANY_RAISING";

    pub fn new(mut multiplicities: Vec<usize>) -> Result<Self> {
        if let Some(invalid) = multiplicities
            .iter()
            .find(|multiplicity| **multiplicity < 2)
        {
            return Err(eyre!(
                "Raised-propagator signatures only accept multiplicities >= 2; found {}.",
                invalid
            ));
        }
        multiplicities.sort_unstable();
        Ok(Self::Exact(multiplicities))
    }

    pub fn canonical(&self) -> String {
        self.to_string()
    }

    fn matches(&self, actual: &Self) -> bool {
        match self {
            Self::Exact(_) => self == actual,
            Self::AnyRaising => {
                matches!(actual, Self::Exact(multiplicities) if !multiplicities.is_empty())
            }
        }
    }
}

impl Display for RaisedPropagatorSignature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Exact(multiplicities) => write!(f, "[{}]", multiplicities.iter().join(",")),
            Self::AnyRaising => write!(f, "{}", Self::ANY_RAISING_KEYWORD),
        }
    }
}

impl FromStr for RaisedPropagatorSignature {
    type Err = eyre::Report;

    fn from_str(raw: &str) -> Result<Self> {
        if raw.trim().eq_ignore_ascii_case(Self::ANY_RAISING_KEYWORD) {
            return Ok(Self::AnyRaising);
        }
        let content = bracket_content(raw, '[', ']')
            .with_context(|| format!("Invalid raised-propagator signature '{raw}'"))?;
        if content.trim().is_empty() {
            return Ok(Self::Exact(Vec::new()));
        }
        let entries = parse_comma_separated_list(content)
            .with_context(|| format!("Invalid raised-propagator signature '{raw}'"))?;
        let multiplicities = entries
            .into_iter()
            .map(|value| {
                value.parse::<usize>().map_err(|_| {
                    eyre!("Invalid raised-propagator multiplicity '{value}' in '{raw}'.")
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Self::new(multiplicities)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CycleSignature(Vec<CycleRequirement>);

impl CycleSignature {
    pub fn new(requirements: Vec<CycleRequirement>) -> Result<Self> {
        if requirements.is_empty() {
            return Err(eyre!(
                "Cycle signatures require at least one cycle requirement."
            ));
        }
        let mut seen = BTreeSet::new();
        for requirement in &requirements {
            if !seen.insert(requirement.clone()) {
                return Err(eyre!(
                    "Duplicate cycle requirement '{}' is ambiguous; specify each required cycle once.",
                    requirement
                ));
            }
        }
        Ok(Self(requirements))
    }

    pub fn parse(raw: &str, model: &Model) -> Result<Self> {
        let content = bracket_content(raw, '[', ']')
            .with_context(|| format!("Invalid cycle signature '{raw}'"))?;
        let requirements = parse_cycle_requirements(content, model)
            .with_context(|| format!("Invalid cycle signature '{raw}'"))?;
        Self::new(requirements)
    }

    pub fn canonical(&self) -> String {
        self.to_string()
    }
}

impl Display for CycleSignature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}]", self.0.iter().join(","))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CycleRequirement(Vec<CycleMatcher>);

impl CycleRequirement {
    pub fn new(mut matchers: Vec<CycleMatcher>) -> Result<Self> {
        if matchers.is_empty() {
            return Err(eyre!("Cycle requirements cannot be empty."));
        }
        matchers.sort();
        matchers.dedup();
        Ok(Self(matchers))
    }
}

impl Display for CycleRequirement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({})", self.0.iter().join(","))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum CycleMatcher {
    Pdg(isize),
    Fermion,
    Ghost,
    Goldstone,
}

impl CycleMatcher {
    fn matches(&self, particle: &ArcParticle) -> bool {
        match self {
            Self::Pdg(pdg) => particle.pdg_code.abs() == *pdg,
            Self::Fermion => particle.is_fermion(),
            Self::Ghost => particle.is_ghost(),
            Self::Goldstone => particle.is_goldstone(),
        }
    }
}

impl Display for CycleMatcher {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Pdg(pdg) => write!(f, "{pdg}"),
            Self::Fermion => f.write_str("fermion"),
            Self::Ghost => f.write_str("ghost"),
            Self::Goldstone => f.write_str("goldstone"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct VertexSignature(BTreeMap<String, usize>);

impl VertexSignature {
    pub fn new(vertex_rule_names: Vec<String>) -> Result<Self> {
        if vertex_rule_names.is_empty() {
            return Err(eyre!("Vertex signatures require at least one vertex rule."));
        }
        let mut counts = BTreeMap::<String, usize>::new();
        for vertex_rule_name in vertex_rule_names {
            if vertex_rule_name.trim().is_empty() {
                return Err(eyre!("Vertex rule names cannot be empty."));
            }
            *counts.entry(vertex_rule_name).or_default() += 1;
        }
        Ok(Self(counts))
    }

    pub fn parse(raw: &str) -> Result<Self> {
        let content = bracket_content(raw, '[', ']')
            .with_context(|| format!("Invalid vertex signature '{raw}'"))?;
        let names = parse_comma_separated_identifiers(content)
            .with_context(|| format!("Invalid vertex signature '{raw}'"))?;
        Self::new(names)
    }

    pub fn vertex_rule_names(&self) -> impl Iterator<Item = &str> {
        self.0.keys().map(String::as_str)
    }

    pub fn canonical(&self) -> String {
        self.to_string()
    }
}

impl Display for VertexSignature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let names = self
            .0
            .iter()
            .flat_map(|(name, count)| std::iter::repeat_n(name, *count))
            .join(",");
        write!(f, "[{names}]")
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ParticleSignature(BTreeSet<isize>);

impl ParticleSignature {
    pub fn new(pdgs: impl IntoIterator<Item = isize>) -> Result<Self> {
        let pdgs = pdgs.into_iter().collect::<BTreeSet<_>>();
        if pdgs.is_empty() {
            return Err(eyre!("Particle signatures require at least one particle."));
        }
        Ok(Self(pdgs))
    }

    pub fn parse(raw: &str, model: &Model) -> Result<Self> {
        let raw = raw.trim();
        let content = if raw.starts_with('[') {
            bracket_content(raw, '[', ']')
        } else {
            bracket_content(raw, '(', ')')
        }
        .with_context(|| format!("Invalid particle signature '{raw}'"))?;
        let tokens = parse_comma_separated_list(content)
            .with_context(|| format!("Invalid particle signature '{raw}'"))?;
        let pdgs = tokens
            .into_iter()
            .map(|token| parse_particle_signature_pdg(&token, model))
            .collect::<Result<Vec<_>>>()?;
        Self::new(pdgs)
    }

    pub fn canonical(&self) -> String {
        self.to_string()
    }
}

impl Display for ParticleSignature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}]", self.0.iter().join(","))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GraphSelectionSignatureInventory {
    pub raised_all: Vec<String>,
    pub raised_massive: Vec<String>,
    pub raised_massless: Vec<String>,
    pub cycles: Vec<String>,
    pub vertices: Vec<String>,
}

impl GraphSelectionSignatureInventory {
    pub fn empty() -> Self {
        Self {
            raised_all: vec!["[]".to_string(), "[2]".to_string()],
            raised_massive: vec!["[]".to_string(), "[2]".to_string()],
            raised_massless: vec!["[]".to_string(), "[2]".to_string()],
            cycles: Vec::new(),
            vertices: Vec::new(),
        }
    }

    pub fn from_master_graphs<'a>(graphs: impl IntoIterator<Item = &'a Graph>) -> Self {
        Self::from_analysis_subjects(graphs.into_iter().map(GraphSelectionSubject::whole_graph))
    }

    pub(crate) fn from_analysis_subjects<'a>(
        subjects: impl IntoIterator<Item = GraphSelectionSubject<'a>>,
    ) -> Self {
        let mut raised_all = BTreeSet::from(["[]".to_string(), "[2]".to_string()]);
        let mut raised_massive = BTreeSet::from(["[]".to_string(), "[2]".to_string()]);
        let mut raised_massless = BTreeSet::from(["[]".to_string(), "[2]".to_string()]);
        let mut cycles = BTreeSet::new();
        let mut vertices = BTreeSet::new();

        for subject in subjects {
            raised_all.insert(
                subject
                    .raised_propagator_signature(RaisedPropagatorScope::All)
                    .canonical(),
            );
            raised_massive.insert(
                subject
                    .raised_propagator_signature(RaisedPropagatorScope::Massive)
                    .canonical(),
            );
            raised_massless.insert(
                subject
                    .raised_propagator_signature(RaisedPropagatorScope::Massless)
                    .canonical(),
            );
            cycles.extend(
                subject
                    .cycle_requirements()
                    .into_iter()
                    .map(|requirement| CycleSignature(vec![requirement]).canonical()),
            );
            if let Some(signature) = subject.vertex_signature() {
                vertices.insert(signature.canonical());
            }
        }

        Self {
            raised_all: raised_all.into_iter().collect(),
            raised_massive: raised_massive.into_iter().collect(),
            raised_massless: raised_massless.into_iter().collect(),
            cycles: cycles.into_iter().collect(),
            vertices: vertices.into_iter().collect(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RaisedCutSignatureInventory {
    pub all: Vec<String>,
    pub massive: Vec<String>,
    pub massless: Vec<String>,
}

impl RaisedCutSignatureInventory {
    pub fn empty() -> Self {
        Self {
            all: vec!["[]".to_string(), "[2]".to_string()],
            massive: vec!["[]".to_string(), "[2]".to_string()],
            massless: vec!["[]".to_string(), "[2]".to_string()],
        }
    }

    pub(crate) fn from_cut_subjects<'a>(
        subjects: impl IntoIterator<Item = GraphCutSelectionSubject<'a>>,
    ) -> Self {
        let mut all = BTreeSet::from(["[]".to_string(), "[2]".to_string()]);
        let mut massive = BTreeSet::from(["[]".to_string(), "[2]".to_string()]);
        let mut massless = BTreeSet::from(["[]".to_string(), "[2]".to_string()]);

        for subject in subjects {
            all.insert(
                subject
                    .raised_cut_signature(RaisedPropagatorScope::All)
                    .canonical(),
            );
            massive.insert(
                subject
                    .raised_cut_signature(RaisedPropagatorScope::Massive)
                    .canonical(),
            );
            massless.insert(
                subject
                    .raised_cut_signature(RaisedPropagatorScope::Massless)
                    .canonical(),
            );
        }

        Self {
            all: all.into_iter().collect(),
            massive: massive.into_iter().collect(),
            massless: massless.into_iter().collect(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum GraphGroupSelectionRule {
    MasterGraphNames {
        polarity: SelectionPolarity,
        graph_names: Vec<String>,
    },
    RaisedPropagators {
        polarity: SelectionPolarity,
        scope: RaisedPropagatorScope,
        signatures: Vec<RaisedPropagatorSignature>,
    },
    RaisedCuts {
        polarity: SelectionPolarity,
        scope: RaisedPropagatorScope,
        signatures: Vec<RaisedPropagatorSignature>,
    },
    Cycles {
        polarity: SelectionPolarity,
        signatures: Vec<CycleSignature>,
    },
    Vertices {
        polarity: SelectionPolarity,
        signatures: Vec<VertexSignature>,
    },
    Particles {
        polarity: SelectionPolarity,
        signatures: Vec<ParticleSignature>,
    },
}

impl GraphGroupSelectionRule {
    fn uses_graph_analysis(&self) -> bool {
        matches!(
            self,
            Self::RaisedPropagators { .. }
                | Self::Cycles { .. }
                | Self::Vertices { .. }
                | Self::Particles { .. }
        )
    }

    fn uses_cut_analysis(&self) -> bool {
        matches!(self, Self::RaisedCuts { .. })
    }

    fn is_authoritative_master_graph_name_rule(&self) -> bool {
        matches!(
            self,
            Self::MasterGraphNames {
                polarity: SelectionPolarity::With,
                ..
            }
        )
    }

    fn resolve_authoritative_master_graph_names(
        &self,
        candidates: &[GraphGroupSelectionCandidate],
        master_name_to_group: &BTreeMap<String, GroupId>,
    ) -> Result<BTreeSet<GroupId>> {
        let Self::MasterGraphNames {
            polarity: SelectionPolarity::With,
            graph_names,
        } = self
        else {
            return Ok(BTreeSet::new());
        };
        Self::resolve_master_graph_name_set(graph_names, candidates, master_name_to_group)
    }

    fn resolve(
        &self,
        candidates: &[GraphGroupSelectionCandidate],
        master_name_to_group: &BTreeMap<String, GroupId>,
    ) -> Result<BTreeSet<GroupId>> {
        match self {
            Self::MasterGraphNames {
                polarity,
                graph_names,
            } => {
                let matched = Self::resolve_master_graph_name_set(
                    graph_names,
                    candidates,
                    master_name_to_group,
                )?;
                Ok(candidates
                    .iter()
                    .filter_map(|candidate| {
                        polarity
                            .keep_if_match(matched.contains(&candidate.group_id))
                            .then_some(candidate.group_id)
                    })
                    .collect())
            }
            Self::RaisedPropagators {
                polarity,
                scope,
                signatures,
            } => matching_candidates(candidates, |candidate| {
                let matched = candidate.analysis_subjects.iter().any(|subject| {
                    let actual = subject.raised_propagator_signature(*scope);
                    signatures
                        .iter()
                        .any(|signature| signature.matches(&actual))
                });
                polarity.keep_if_match(matched)
            }),
            Self::RaisedCuts {
                polarity,
                scope,
                signatures,
            } => matching_candidates(candidates, |candidate| {
                let matched = candidate.cut_subjects.iter().any(|subject| {
                    let actual = subject.raised_cut_signature(*scope);
                    signatures
                        .iter()
                        .any(|signature| signature.matches(&actual))
                });
                polarity.keep_if_match(matched)
            }),
            Self::Cycles {
                polarity,
                signatures,
            } => matching_candidates(candidates, |candidate| {
                let matched = candidate.analysis_subjects.iter().any(|subject| {
                    signatures
                        .iter()
                        .any(|signature| subject.matches_cycle_signature(signature))
                });
                polarity.keep_if_match(matched)
            }),
            Self::Vertices {
                polarity,
                signatures,
            } => matching_candidates(candidates, |candidate| {
                let matched = candidate.analysis_subjects.iter().any(|subject| {
                    let counts = subject.vertex_rule_name_counts();
                    signatures.iter().any(|signature| {
                        signature.0.iter().all(|(name, required_count)| {
                            counts.get(name).copied().unwrap_or_default() >= *required_count
                        })
                    })
                });
                polarity.keep_if_match(matched)
            }),
            Self::Particles {
                polarity,
                signatures,
            } => matching_candidates(candidates, |candidate| {
                let matched = candidate.analysis_subjects.iter().any(|subject| {
                    let pdgs = subject.particle_pdgs();
                    signatures
                        .iter()
                        .any(|signature| signature.0.iter().all(|pdg| pdgs.contains(pdg)))
                });
                polarity.keep_if_match(matched)
            }),
        }
    }

    fn resolve_master_graph_name_set(
        graph_names: &[String],
        candidates: &[GraphGroupSelectionCandidate],
        master_name_to_group: &BTreeMap<String, GroupId>,
    ) -> Result<BTreeSet<GroupId>> {
        if graph_names.is_empty() {
            return Err(eyre!(
                "Graph-name selection requires at least one graph name."
            ));
        }

        let mut seen = BTreeSet::<&str>::new();
        for graph_name in graph_names {
            if !seen.insert(graph_name.as_str()) {
                return Err(eyre!(
                    "Duplicate graph name '{}' in graph-name selection.",
                    graph_name
                ));
            }
        }

        graph_names
            .iter()
            .map(|graph_name| {
                resolve_master_graph_name(graph_name, candidates, master_name_to_group)
                    .with_context(|| {
                        format!(
                            "Available master graphs are: {}",
                            master_name_to_group.keys().join(", ")
                        )
                    })
            })
            .collect::<Result<BTreeSet<_>>>()
    }
}

#[derive(Clone)]
struct GraphGroupSelectionCandidate<'a> {
    group_id: GroupId,
    master_graph_id: usize,
    master_graph_name: String,
    analysis_subjects: Vec<GraphSelectionSubject<'a>>,
    cut_subjects: Vec<GraphCutSelectionSubject<'a>>,
    graph_names: Vec<String>,
}

#[derive(Clone)]
pub(crate) struct GraphSelectionSubject<'a> {
    graph: &'a Graph,
    subgraph: Option<SuBitGraph>,
    raised_edge_policy: RaisedEdgePolicy,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum RaisedEdgePolicy {
    LoopDependentOnly,
    AllInternalInSubject,
}

impl<'a> GraphSelectionSubject<'a> {
    pub(crate) fn whole_graph(graph: &'a Graph) -> Self {
        Self {
            graph,
            subgraph: None,
            raised_edge_policy: RaisedEdgePolicy::LoopDependentOnly,
        }
    }

    #[cfg(test)]
    pub(crate) fn subgraph(graph: &'a Graph, subgraph: SuBitGraph) -> Self {
        Self {
            graph,
            subgraph: Some(subgraph),
            raised_edge_policy: RaisedEdgePolicy::LoopDependentOnly,
        }
    }

    pub(crate) fn cut_side_amplitude_subgraph(graph: &'a Graph, subgraph: SuBitGraph) -> Self {
        Self {
            graph,
            subgraph: Some(subgraph),
            raised_edge_policy: RaisedEdgePolicy::AllInternalInSubject,
        }
    }

    fn raised_edge_subgraph(&self) -> SuBitGraph {
        match &self.subgraph {
            Some(subgraph) => subgraph
                .subtract(&self.graph.external_filter::<SuBitGraph>())
                .subtract(&self.graph.initial_state_cut.left)
                .subtract(&self.graph.initial_state_cut.right),
            None => self.graph.underlying.full_filter(),
        }
    }

    fn internal_edge_subgraph(&self) -> SuBitGraph {
        let mut subgraph = match &self.subgraph {
            Some(subgraph) => subgraph.clone(),
            None => self.graph.underlying.full_filter(),
        }
        .subtract(&self.graph.external_filter::<SuBitGraph>())
        .subtract(&self.graph.initial_state_cut.left)
        .subtract(&self.graph.initial_state_cut.right);

        for (pair, _, edge) in self.graph.underlying.iter_edges() {
            if edge.data.is_dummy {
                subgraph.sub(pair);
            }
        }
        subgraph
    }

    fn vertex_subgraph(&self) -> Option<SuBitGraph> {
        self.subgraph
            .as_ref()
            .map(|_| self.internal_edge_subgraph())
    }
}

#[derive(Clone)]
pub(crate) struct GraphCutSelectionSubject<'a> {
    graph: &'a Graph,
    cut_edges: SuBitGraph,
}

impl<'a> GraphCutSelectionSubject<'a> {
    pub(crate) fn new(graph: &'a Graph, cut_edges: SuBitGraph) -> Self {
        Self { graph, cut_edges }
    }

    fn raised_cut_signature(&self, scope: RaisedPropagatorScope) -> RaisedPropagatorSignature {
        let cut_edge_ids = self
            .graph
            .underlying
            .iter_edges_of(&self.cut_edges)
            .map(|(_, edge_id, _)| edge_id)
            .collect::<BTreeSet<_>>();

        let multiplicities = GraphSelectionSubject::whole_graph(self.graph)
            .raised_edge_groups()
            .into_iter()
            .filter(|group| group.len() > 1)
            .filter(|group| raised_group_matches_scope(self.graph, group, scope))
            .filter(|group| group.iter().any(|edge| cut_edge_ids.contains(edge)))
            .map(|group| group.len())
            .collect::<Vec<_>>();

        RaisedPropagatorSignature::new(multiplicities)
            .expect("group lengths are always valid raised-cut multiplicities")
    }
}

trait GraphSelectionAnalysis {
    fn raised_propagator_signature(
        &self,
        scope: RaisedPropagatorScope,
    ) -> RaisedPropagatorSignature;
    fn cycle_requirements(&self) -> BTreeSet<CycleRequirement>;
    fn cycle_particle_sets(&self) -> Vec<Vec<ArcParticle>>;
    fn matches_cycle_signature(&self, signature: &CycleSignature) -> bool;
    fn vertex_rule_name_counts(&self) -> BTreeMap<String, usize>;
    fn particle_pdgs(&self) -> BTreeSet<isize>;
}

impl GraphSelectionAnalysis for Graph {
    fn raised_propagator_signature(
        &self,
        scope: RaisedPropagatorScope,
    ) -> RaisedPropagatorSignature {
        GraphSelectionSubject::whole_graph(self).raised_propagator_signature(scope)
    }

    fn cycle_requirements(&self) -> BTreeSet<CycleRequirement> {
        GraphSelectionSubject::whole_graph(self).cycle_requirements()
    }

    fn cycle_particle_sets(&self) -> Vec<Vec<ArcParticle>> {
        GraphSelectionSubject::whole_graph(self).cycle_particle_sets()
    }

    fn matches_cycle_signature(&self, signature: &CycleSignature) -> bool {
        GraphSelectionSubject::whole_graph(self).matches_cycle_signature(signature)
    }

    fn vertex_rule_name_counts(&self) -> BTreeMap<String, usize> {
        GraphSelectionSubject::whole_graph(self).vertex_rule_name_counts()
    }

    fn particle_pdgs(&self) -> BTreeSet<isize> {
        GraphSelectionSubject::whole_graph(self).particle_pdgs()
    }
}

impl GraphSelectionAnalysis for GraphSelectionSubject<'_> {
    fn raised_propagator_signature(
        &self,
        scope: RaisedPropagatorScope,
    ) -> RaisedPropagatorSignature {
        let multiplicities = self
            .raised_edge_groups()
            .into_iter()
            .filter(|group| group.len() > 1)
            .filter(|group| raised_group_matches_scope(self.graph, group, scope))
            .map(|group| group.len())
            .collect::<Vec<_>>();
        RaisedPropagatorSignature::new(multiplicities)
            .expect("group lengths are always valid raised-propagator multiplicities")
    }

    fn cycle_requirements(&self) -> BTreeSet<CycleRequirement> {
        self.cycle_particle_sets()
            .into_iter()
            .filter_map(cycle_requirement)
            .collect()
    }

    fn cycle_particle_sets(&self) -> Vec<Vec<ArcParticle>> {
        self.internal_simple_cycles()
            .into_iter()
            .filter_map(|cycle| cycle_particles(self.graph, &cycle))
            .collect()
    }

    fn matches_cycle_signature(&self, signature: &CycleSignature) -> bool {
        let cycle_particle_sets = self.cycle_particle_sets();
        signature.0.iter().all(|requirement| {
            cycle_particle_sets
                .iter()
                .any(|particles| cycle_matches_requirement(particles, requirement))
        })
    }

    fn vertex_rule_name_counts(&self) -> BTreeMap<String, usize> {
        let mut counts = BTreeMap::<String, usize>::new();
        if let Some(subgraph) = self.vertex_subgraph() {
            for (_, _, vertex) in self.graph.underlying.iter_nodes_of(&subgraph) {
                if let Some(vertex_rule) = &vertex.vertex_rule {
                    *counts.entry(vertex_rule.name.to_string()).or_default() += 1;
                }
            }
        } else {
            for (_, _, vertex) in self.graph.underlying.iter_nodes() {
                if let Some(vertex_rule) = &vertex.vertex_rule {
                    *counts.entry(vertex_rule.name.to_string()).or_default() += 1;
                }
            }
        }
        counts
    }

    fn particle_pdgs(&self) -> BTreeSet<isize> {
        self.graph
            .underlying
            .iter_edges_of(&self.internal_edge_subgraph())
            .filter_map(|(_, edge_id, _)| self.graph[edge_id].particle())
            .map(|particle| particle.pdg_code.abs())
            .collect()
    }
}

trait GraphCycleSelectionAnalysis {
    fn internal_simple_cycles(&self) -> Vec<Cycle>;
}

impl GraphSelectionSubject<'_> {
    fn raised_edge_groups(&self) -> Vec<Vec<EdgeIndex>> {
        let mut result = Vec::<Vec<EdgeIndex>>::new();
        let subgraph = self.raised_edge_subgraph();

        for (_, edge_index, edge) in self.graph.underlying.iter_edges_of(&subgraph) {
            let loop_independent = self.graph.loop_momentum_basis.edge_signatures[edge_index]
                .internal
                .iter()
                .all(|sign| sign.is_zero());
            if edge.data.is_dummy
                || (self.raised_edge_policy == RaisedEdgePolicy::LoopDependentOnly
                    && loop_independent)
            {
                continue;
            }

            let group_position = result.iter().position(|group| {
                group.iter().all(|edge| {
                    self.graph
                        .loop_momentum_basis
                        .edges_are_raised(*edge, edge_index)
                        && self.graph[edge_index].mass == self.graph[*edge].mass
                })
            });

            if let Some(pos) = group_position {
                result[pos].push(edge_index);
            } else {
                result.push(vec![edge_index]);
            }
        }

        result.iter_mut().for_each(|group| group.sort());
        result
    }

    fn vertex_signature(&self) -> Option<VertexSignature> {
        let counts = self.vertex_rule_name_counts();
        if counts.is_empty() {
            None
        } else {
            Some(VertexSignature(counts))
        }
    }
}

impl GraphCycleSelectionAnalysis for GraphSelectionSubject<'_> {
    fn internal_simple_cycles(&self) -> Vec<Cycle> {
        let subgraph = self.internal_edge_subgraph();
        if self.graph.underlying.cyclotomatic_number(&subgraph) == 0 {
            return Vec::new();
        }
        let basis = self.graph.underlying.cycle_basis_of(&subgraph).0;
        Cycle::all_sum_powerset_filter_map(&basis, &|mut cycle| {
            if cycle.is_circuit(&self.graph.underlying) {
                cycle.loop_count = Some(1);
                Some(cycle)
            } else {
                None
            }
        })
        .map(|cycles| cycles.into_iter().collect())
        .unwrap_or_default()
    }
}

fn raised_group_matches_scope(
    graph: &Graph,
    group: &[EdgeIndex],
    scope: RaisedPropagatorScope,
) -> bool {
    match scope {
        RaisedPropagatorScope::All => true,
        RaisedPropagatorScope::Massive => group
            .iter()
            .all(|edge| !matches!(graph[*edge].mass, EdgeMass::Zero)),
        RaisedPropagatorScope::Massless => group
            .iter()
            .all(|edge| matches!(graph[*edge].mass, EdgeMass::Zero)),
    }
}

fn cycle_particles(graph: &Graph, cycle: &Cycle) -> Option<Vec<ArcParticle>> {
    let mut particles = BTreeMap::<EdgeIndex, ArcParticle>::new();
    for hedge in cycle.filter.included_iter() {
        let edge_id = graph.underlying[&hedge];
        if graph[edge_id].is_dummy {
            continue;
        }
        let particle = graph[edge_id].particle()?;
        particles.insert(edge_id, particle);
    }
    Some(particles.into_values().collect())
}

fn cycle_requirement(particles: Vec<ArcParticle>) -> Option<CycleRequirement> {
    let matchers = particles
        .into_iter()
        .map(|particle| CycleMatcher::Pdg(particle.pdg_code.abs()))
        .collect::<Vec<_>>();
    CycleRequirement::new(matchers).ok()
}

fn cycle_matches_requirement(particles: &[ArcParticle], requirement: &CycleRequirement) -> bool {
    if particles.is_empty() {
        return false;
    }
    let mut matched_requirements = vec![false; requirement.0.len()];
    for particle in particles {
        let mut particle_matched = false;
        for (index, matcher) in requirement.0.iter().enumerate() {
            if matcher.matches(particle) {
                matched_requirements[index] = true;
                particle_matched = true;
            }
        }
        if !particle_matched {
            return false;
        }
    }
    matched_requirements.into_iter().all(|matched| matched)
}

fn matching_candidates(
    candidates: &[GraphGroupSelectionCandidate],
    predicate: impl Fn(&GraphGroupSelectionCandidate) -> bool,
) -> Result<BTreeSet<GroupId>> {
    Ok(candidates
        .iter()
        .filter(|candidate| predicate(candidate))
        .map(|candidate| candidate.group_id)
        .collect())
}

fn resolve_master_graph_name(
    graph_name: &str,
    candidates: &[GraphGroupSelectionCandidate],
    master_name_to_group: &BTreeMap<String, GroupId>,
) -> Result<GroupId> {
    if let Some(group_id) = master_name_to_group.get(graph_name) {
        return Ok(*group_id);
    }

    let containing_groups = candidates
        .iter()
        .filter(|candidate| candidate.graph_names.iter().any(|name| name == graph_name))
        .collect::<Vec<_>>();
    match containing_groups.as_slice() {
        [] => Err(eyre!("Unknown graph '{}'.", graph_name)),
        [candidate] => Err(eyre!(
            "Graph '{}' is not the master graph of its group; use '{}' instead.",
            graph_name,
            candidate.master_graph_name
        )),
        many => Err(eyre!(
            "Graph name '{}' is ambiguous across graph groups with masters: {}.",
            graph_name,
            many.iter()
                .map(|candidate| format!(
                    "{} (group {}, master graph id {})",
                    candidate.master_graph_name, candidate.group_id.0, candidate.master_graph_id
                ))
                .join(", ")
        )),
    }
}

fn bracket_content(raw: &str, open: char, close: char) -> Result<&str> {
    let raw = raw.trim();
    if !raw.starts_with(open) || !raw.ends_with(close) {
        return Err(eyre!("Expected value enclosed by '{open}' and '{close}'."));
    }
    Ok(&raw[open.len_utf8()..raw.len() - close.len_utf8()])
}

fn parse_comma_separated_identifiers(content: &str) -> Result<Vec<String>> {
    parse_comma_separated_list(content)
}

fn parse_comma_separated_list(content: &str) -> Result<Vec<String>> {
    let trimmed = content.trim();
    if trimmed.is_empty() {
        return Ok(Vec::new());
    }

    let tokens = trimmed.split(',').map(str::trim).collect::<Vec<_>>();
    let last_index = tokens.len().saturating_sub(1);
    tokens
        .into_iter()
        .enumerate()
        .filter_map(|(index, token)| {
            if token.is_empty() {
                if index == last_index {
                    None
                } else {
                    Some(Err(eyre!("Empty entry in comma-separated list.")))
                }
            } else {
                Some(Ok(token.to_string()))
            }
        })
        .collect()
}

fn parse_cycle_requirements(content: &str, model: &Model) -> Result<Vec<CycleRequirement>> {
    let mut requirements = Vec::new();
    let mut rest = content.trim();
    while !rest.is_empty() {
        if !rest.starts_with('(') {
            return Err(eyre!("Expected '(' at '{}'.", rest));
        }
        let close = rest
            .find(')')
            .ok_or_else(|| eyre!("Unclosed cycle requirement in '{content}'."))?;
        let tuple_content = &rest[1..close];
        let matchers = parse_cycle_matchers(tuple_content, model)?;
        requirements.push(CycleRequirement::new(matchers)?);
        rest = rest[close + 1..].trim_start();
        if rest.is_empty() {
            break;
        }
        if !rest.starts_with(',') {
            return Err(eyre!(
                "Expected ',' after cycle requirement in '{content}'."
            ));
        }
        rest = rest[1..].trim_start();
        if rest.is_empty() {
            break;
        }
    }
    Ok(requirements)
}

fn parse_cycle_matchers(content: &str, model: &Model) -> Result<Vec<CycleMatcher>> {
    let trimmed = content.trim();
    if trimmed.is_empty() {
        return Err(eyre!("Cycle requirements cannot be empty."));
    }
    let tokens = trimmed
        .split(',')
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .map(str::to_string)
        .collect::<Vec<_>>();
    if tokens.is_empty() {
        return Err(eyre!("Cycle requirements cannot be empty."));
    }
    tokens
        .into_iter()
        .map(|token| parse_cycle_matcher(&token, model))
        .collect()
}

fn parse_cycle_matcher(token: &str, model: &Model) -> Result<CycleMatcher> {
    match token {
        "fermion" => return Ok(CycleMatcher::Fermion),
        "ghost" => return Ok(CycleMatcher::Ghost),
        "goldstone" => return Ok(CycleMatcher::Goldstone),
        _ => {}
    }
    if let Ok(pdg) = token.parse::<isize>() {
        if pdg < 0 {
            return Err(eyre!(
                "Cycle signatures use particle PDGs only; specify {} instead of anti-particle PDG {}.",
                pdg.abs(),
                pdg
            ));
        }
        model.try_get_particle_from_pdg(pdg)?;
        return Ok(CycleMatcher::Pdg(pdg));
    }

    let particle = model.try_get_particle(token)?;
    if particle.pdg_code < 0 || particle.name.as_str() != token {
        return Err(eyre!(
            "Cycle signatures use particles only; specify '{}' instead of anti-particle '{}'.",
            particle.get_anti_particle(model).name,
            token
        ));
    }
    Ok(CycleMatcher::Pdg(particle.pdg_code))
}

fn parse_particle_signature_pdg(token: &str, model: &Model) -> Result<isize> {
    if let Ok(pdg) = token.parse::<isize>() {
        if pdg == 0 {
            return Err(eyre!("Particle signatures do not accept PDG 0."));
        }
        let abs_pdg = pdg.abs();
        model
            .try_get_particle_from_pdg(abs_pdg)
            .or_else(|_| model.try_get_particle_from_pdg(-abs_pdg))?;
        return Ok(abs_pdg);
    }

    let particle = model.try_get_particle(token)?;
    Ok(particle.pdg_code.abs())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, OnceLock};

    use crate::{
        dot,
        graph::{
            Graph, LoopMomentumBasis,
            edge::EdgeMass,
            parse::{IntoGraph, complete_group_parsing},
        },
        initialisation::test_initialise,
        model::{ArcVertexRule, ColorStructure, ParameterName, Particle, UFOSymbol, VertexRule},
        momentum::signature::LoopExtSignature,
    };

    fn scalar_model() -> &'static Model {
        static MODEL: OnceLock<Model> = OnceLock::new();
        MODEL.get_or_init(|| {
            test_initialise().expect("test initialization should succeed");
            crate::utils::load_generic_model("scalars")
        })
    }

    fn minimal_model() -> Model {
        let mut model = Model::default();
        for particle in [
            test_particle(3, "s", "s~", 2, 3, false, 0),
            test_particle(-3, "s~", "s", 2, -3, false, 0),
            test_particle(4, "c", "c~", 2, 3, false, 0),
            test_particle(11, "e-", "e+", 2, 1, false, 0),
            test_particle(-11, "e+", "e-", 2, 1, false, 0),
            test_particle(21, "g", "g", 3, 8, false, 0),
            test_particle(22, "a", "a", 3, 1, false, 0),
            test_particle(25, "h", "h", 1, 1, true, 0),
            test_particle(82, "ghG", "ghG~", 1, 8, false, 1),
        ] {
            model
                .particle_pdg_to_position
                .insert(particle.pdg_code, model.particles.len());
            model
                .particle_name_to_position
                .insert(particle.name.clone(), model.particles.len());
            model.particles.push(ArcParticle(Arc::new(particle)));
        }
        model
    }

    fn test_particle(
        pdg_code: isize,
        name: &str,
        antiname: &str,
        spin: isize,
        color: isize,
        goldstone: bool,
        ghost_number: isize,
    ) -> Particle {
        Particle {
            pdg_code,
            name: name.into(),
            antiname: antiname.into(),
            spin,
            color,
            mass: ParameterName(UFOSymbol::zero()),
            width: ParameterName(UFOSymbol::zero()),
            texname: name.into(),
            antitexname: antiname.into(),
            charge: 0.0,
            ghost_number,
            lepton_number: 0,
            y_charge: 0,
            goldstone,
        }
    }

    fn raised_test_graph() -> Result<Graph> {
        let mut graph: Graph = dot!(
            digraph raised {
                edge [particle=scalar_1]
                ext [style=invis]
                ext -> A [id=3]
                B -> ext [id=4]
                A -> B [id=0]
                A -> B [id=1]
                A -> B [id=2]
                A -> B [id=5]
                A -> B [id=6]
            },
            scalar_model()
        )?;

        graph.loop_momentum_basis = LoopMomentumBasis {
            tree: graph.underlying.empty_subgraph(),
            loop_edges: vec![
                EdgeIndex::from(0),
                EdgeIndex::from(1),
                EdgeIndex::from(2),
                EdgeIndex::from(5),
                EdgeIndex::from(6),
            ]
            .into(),
            ext_edges: vec![EdgeIndex::from(3), EdgeIndex::from(4)].into(),
            edge_signatures: graph
                .underlying
                .new_edgevec(|_, edge_id, _| match edge_id.0 {
                    0 => LoopExtSignature::from((vec![1], vec![])),
                    1 => LoopExtSignature::from((vec![-1], vec![])),
                    2 => LoopExtSignature::from((vec![0], vec![])),
                    5 => LoopExtSignature::from((vec![1], vec![])),
                    6 => LoopExtSignature::from((vec![-1], vec![])),
                    _ => LoopExtSignature::from((vec![0], vec![])),
                }),
        };

        for edge_id in [0, 1, 2] {
            graph.underlying[EdgeIndex::from(edge_id)].mass = EdgeMass::Zero;
        }
        Ok(graph)
    }

    fn tree_raised_test_graph() -> Result<Graph> {
        let mut graph: Graph = dot!(
            digraph tree_raised {
                edge [particle=scalar_1]
                ext [style=invis]
                ext -> A [id=3]
                B -> ext [id=4]
                A -> B [id=0]
                A -> B [id=1]
                A -> B [id=2]
            },
            scalar_model()
        )?;

        graph.loop_momentum_basis = LoopMomentumBasis {
            tree: graph.underlying.empty_subgraph(),
            loop_edges: vec![EdgeIndex::from(2)].into(),
            ext_edges: vec![EdgeIndex::from(3), EdgeIndex::from(4)].into(),
            edge_signatures: graph
                .underlying
                .new_edgevec(|_, edge_id, _| match edge_id.0 {
                    0 => LoopExtSignature::from((vec![0], vec![1])),
                    1 => LoopExtSignature::from((vec![0], vec![-1])),
                    2 => LoopExtSignature::from((vec![1], vec![])),
                    _ => LoopExtSignature::from((vec![0], vec![])),
                }),
        };

        Ok(graph)
    }

    fn cycle_test_graph() -> Result<Graph> {
        dot!(
            digraph cycle_test {
                node [num="1"]
                A -> B [particle=scalar_0, id=0]
                B -> C [particle=scalar_1, id=1]
                C -> A [particle=scalar_0, id=2]
                B -> D [particle=scalar_2, id=3]
                D -> C [particle=scalar_2, id=4]
            },
            scalar_model()
        )
    }

    fn vertex_rule(name: &str) -> ArcVertexRule {
        ArcVertexRule(Arc::new(VertexRule {
            name: name.into(),
            couplings: Vec::new(),
            lorentz_structures: Vec::new(),
            particles: Vec::new(),
            color_structures: ColorStructure {
                color_structure: Vec::new(),
            },
            dod: 0,
        }))
    }

    fn graph_with_vertex_rules(name: &str, vertex_names: &[&str]) -> Result<Graph> {
        let mut graph: Graph = "digraph vertex_test {
            node [num=\"1\"]
            edge [particle=scalar_1]
            A -> B [id=0]
            B -> C [id=1]
            C -> A [id=2]
        }"
        .into_graph(scalar_model())?;
        graph.name = name.to_string();
        let vertex_rules = vertex_names
            .iter()
            .map(|name| vertex_rule(name))
            .collect::<Vec<_>>();
        for ((_, _, vertex), vertex_rule) in graph
            .underlying
            .iter_nodes_mut()
            .zip(vertex_rules.into_iter().cycle())
        {
            vertex.vertex_rule = Some(vertex_rule);
        }
        Ok(graph)
    }

    fn particle_test_graph(name: &str) -> Result<Graph> {
        let mut graph: Graph = "digraph particle_test {
            node [num=\"1\"]
            A -> B [particle=scalar_0, id=0]
            B -> C [particle=scalar_1, id=1]
            C -> A [particle=scalar_2, id=2]
        }"
        .into_graph(scalar_model())?;
        graph.name = name.to_string();
        Ok(graph)
    }

    fn subgraph_from_edge_ids(graph: &Graph, edge_ids: &[usize]) -> SuBitGraph {
        let wanted = edge_ids.iter().copied().collect::<BTreeSet<_>>();
        let mut subgraph: SuBitGraph = graph.underlying.empty_subgraph();
        for (pair, edge_id, _) in graph.underlying.iter_edges() {
            if wanted.contains(&edge_id.0) {
                subgraph.add(pair);
            }
        }
        subgraph
    }

    fn vertex_test_subjects_by_graph_id(
        graph_id: usize,
        graph: &Graph,
    ) -> Result<Vec<GraphSelectionSubject<'_>>> {
        let edge_ids = if graph_id == 0 { &[0][..] } else { &[1][..] };
        Ok(vec![GraphSelectionSubject::subgraph(
            graph,
            subgraph_from_edge_ids(graph, edge_ids),
        )])
    }

    fn raised_cut_test_subjects_by_graph_id(
        graph_id: usize,
        graph: &Graph,
    ) -> Result<Vec<GraphCutSelectionSubject<'_>>> {
        let edge_ids = if graph_id == 0 { &[0][..] } else { &[2][..] };
        Ok(vec![GraphCutSelectionSubject::new(
            graph,
            subgraph_from_edge_ids(graph, edge_ids),
        )])
    }

    fn raised_any_test_subjects_by_graph_id(
        graph_id: usize,
        graph: &Graph,
    ) -> Result<Vec<GraphSelectionSubject<'_>>> {
        let subject = if graph_id == 0 {
            GraphSelectionSubject::whole_graph(graph)
        } else {
            GraphSelectionSubject::subgraph(graph, subgraph_from_edge_ids(graph, &[0, 2]))
        };
        Ok(vec![subject])
    }

    #[test]
    fn raised_signature_parses_and_canonicalizes() {
        assert_eq!(
            "[2,3,4]",
            RaisedPropagatorSignature::from_str("[3,2,4]")
                .unwrap()
                .canonical()
        );
        assert_eq!(
            "[]",
            RaisedPropagatorSignature::from_str("[]")
                .unwrap()
                .canonical()
        );
        assert_eq!(
            "[2]",
            RaisedPropagatorSignature::from_str("[2,]")
                .unwrap()
                .canonical()
        );
        assert_eq!(
            "ANY_RAISING",
            RaisedPropagatorSignature::from_str("ANY_RAISING")
                .unwrap()
                .canonical()
        );
        assert!(RaisedPropagatorSignature::from_str("[1]").is_err());
    }

    #[test]
    fn any_raising_matches_any_non_empty_raised_signature() -> Result<()> {
        let any = RaisedPropagatorSignature::from_str("ANY_RAISING")?;
        let none = RaisedPropagatorSignature::from_str("[]")?;
        let raised = RaisedPropagatorSignature::from_str("[2]")?;

        assert!(any.matches(&raised));
        assert!(!any.matches(&none));
        assert!(raised.matches(&raised));
        assert!(!raised.matches(&RaisedPropagatorSignature::from_str("[3]")?));

        Ok(())
    }

    #[test]
    fn raised_signature_groups_up_to_sign_and_splits_by_mass_scope() -> Result<()> {
        let graph = raised_test_graph()?;

        assert_eq!(
            "[2,2]",
            graph
                .raised_propagator_signature(RaisedPropagatorScope::All)
                .canonical()
        );
        assert_eq!(
            "[2]",
            graph
                .raised_propagator_signature(RaisedPropagatorScope::Massless)
                .canonical()
        );
        assert_eq!(
            "[2]",
            graph
                .raised_propagator_signature(RaisedPropagatorScope::Massive)
                .canonical()
        );

        Ok(())
    }

    #[test]
    fn subgraph_raised_signature_only_uses_edges_in_view() -> Result<()> {
        let graph = raised_test_graph()?;
        let repeated_side =
            GraphSelectionSubject::subgraph(&graph, subgraph_from_edge_ids(&graph, &[0, 1, 2]));
        let cut_edge_removed_side =
            GraphSelectionSubject::subgraph(&graph, subgraph_from_edge_ids(&graph, &[0, 2]));

        assert_eq!(
            "[2]",
            repeated_side
                .raised_propagator_signature(RaisedPropagatorScope::All)
                .canonical()
        );
        assert_eq!(
            "[]",
            cut_edge_removed_side
                .raised_propagator_signature(RaisedPropagatorScope::All)
                .canonical()
        );

        Ok(())
    }

    #[test]
    fn cut_side_amplitude_raised_signature_includes_tree_like_internal_edges() -> Result<()> {
        let graph = tree_raised_test_graph()?;
        let tree_raised_edges = subgraph_from_edge_ids(&graph, &[0, 1]);
        let normal_subject = GraphSelectionSubject::subgraph(&graph, tree_raised_edges.clone());
        let cut_side_subject =
            GraphSelectionSubject::cut_side_amplitude_subgraph(&graph, tree_raised_edges.clone());

        assert_eq!(
            "[]",
            normal_subject
                .raised_propagator_signature(RaisedPropagatorScope::All)
                .canonical()
        );
        assert_eq!(
            "[2]",
            cut_side_subject
                .raised_propagator_signature(RaisedPropagatorScope::All)
                .canonical()
        );
        assert_eq!(
            "[]",
            GraphCutSelectionSubject::new(&graph, tree_raised_edges)
                .raised_cut_signature(RaisedPropagatorScope::All)
                .canonical()
        );

        Ok(())
    }

    #[test]
    fn raised_cut_signature_counts_touched_raised_groups_once() -> Result<()> {
        let graph = raised_test_graph()?;

        assert_eq!(
            "[2]",
            GraphCutSelectionSubject::new(&graph, subgraph_from_edge_ids(&graph, &[0]))
                .raised_cut_signature(RaisedPropagatorScope::All)
                .canonical()
        );
        assert_eq!(
            "[2]",
            GraphCutSelectionSubject::new(&graph, subgraph_from_edge_ids(&graph, &[0, 1]))
                .raised_cut_signature(RaisedPropagatorScope::All)
                .canonical()
        );
        assert_eq!(
            "[2,2]",
            GraphCutSelectionSubject::new(&graph, subgraph_from_edge_ids(&graph, &[0, 5]))
                .raised_cut_signature(RaisedPropagatorScope::All)
                .canonical()
        );
        assert_eq!(
            "[]",
            GraphCutSelectionSubject::new(&graph, subgraph_from_edge_ids(&graph, &[2]))
                .raised_cut_signature(RaisedPropagatorScope::All)
                .canonical()
        );

        Ok(())
    }

    #[test]
    fn raised_cut_signature_respects_mass_scope() -> Result<()> {
        let graph = raised_test_graph()?;
        let cut_subject =
            GraphCutSelectionSubject::new(&graph, subgraph_from_edge_ids(&graph, &[0, 5]));

        assert_eq!(
            "[2,2]",
            cut_subject
                .raised_cut_signature(RaisedPropagatorScope::All)
                .canonical()
        );
        assert_eq!(
            "[2]",
            cut_subject
                .raised_cut_signature(RaisedPropagatorScope::Massless)
                .canonical()
        );
        assert_eq!(
            "[2]",
            cut_subject
                .raised_cut_signature(RaisedPropagatorScope::Massive)
                .canonical()
        );

        Ok(())
    }

    #[test]
    fn vertex_signature_preserves_multiplicity() {
        let signature = VertexSignature::parse("[V_9,V_6,V_9]").unwrap();
        assert_eq!(signature.0.get("V_9"), Some(&2));
        assert_eq!(signature.0.get("V_6"), Some(&1));
        assert_eq!(signature.canonical(), "[V_6,V_9,V_9]");
        assert_eq!(
            VertexSignature::parse("[V_1,]").unwrap().canonical(),
            "[V_1]"
        );
        assert!(VertexSignature::parse("[]").is_err());
    }

    #[test]
    fn particle_signature_parses_names_pdgs_antiparticles_and_tuples() {
        let model = minimal_model();
        assert_eq!(
            "[11,21]",
            ParticleSignature::parse("[e+,g]", &model)
                .unwrap()
                .canonical()
        );
        assert_eq!(
            "[11,21]",
            ParticleSignature::parse("(e-,-21,)", &model)
                .unwrap()
                .canonical()
        );
        assert!(ParticleSignature::parse("[]", &model).is_err());
        assert!(ParticleSignature::parse("[0]", &model).is_err());
    }

    #[test]
    fn cycle_signature_parses_and_canonicalizes() {
        let model = minimal_model();
        assert_eq!(
            "[(3)]",
            CycleSignature::parse("[(3)]", &model).unwrap().canonical()
        );
        assert_eq!(
            "[(3)]",
            CycleSignature::parse("[(3,),]", &model)
                .unwrap()
                .canonical()
        );
        assert_eq!(
            "[(3,21),(4,22)]",
            CycleSignature::parse("[(3,21), (4,22)]", &model)
                .unwrap()
                .canonical()
        );
        assert_eq!(
            "[(fermion,ghost,goldstone)]",
            CycleSignature::parse("[(fermion,ghost,goldstone)]", &model)
                .unwrap()
                .canonical()
        );
    }

    #[test]
    fn cycle_requirements_include_simple_cycles_beyond_basis() -> Result<()> {
        let graph = cycle_test_graph()?;
        let scalar_0 = scalar_model().try_get_particle("scalar_0")?.pdg_code.abs();
        let scalar_1 = scalar_model().try_get_particle("scalar_1")?.pdg_code.abs();
        let scalar_2 = scalar_model().try_get_particle("scalar_2")?.pdg_code.abs();
        let requirements = graph
            .cycle_requirements()
            .into_iter()
            .map(|requirement| requirement.to_string())
            .collect::<BTreeSet<_>>();

        let basis_cycle_a = format!("({},{})", scalar_0.min(scalar_1), scalar_0.max(scalar_1));
        let basis_cycle_b = format!("({},{})", scalar_1.min(scalar_2), scalar_1.max(scalar_2));
        let combined_cycle = format!("({},{})", scalar_0.min(scalar_2), scalar_0.max(scalar_2));

        assert!(
            requirements.contains(&basis_cycle_a),
            "missing first basis cycle in {requirements:?}"
        );
        assert!(
            requirements.contains(&basis_cycle_b),
            "missing second basis cycle in {requirements:?}"
        );
        assert!(
            requirements.contains(&combined_cycle),
            "missing combined simple cycle in {requirements:?}"
        );

        Ok(())
    }

    #[test]
    fn subgraph_cycles_are_broken_by_removed_edges() -> Result<()> {
        let graph = cycle_test_graph()?;
        let full_subject = GraphSelectionSubject::whole_graph(&graph);
        let cut_side_subject =
            GraphSelectionSubject::subgraph(&graph, subgraph_from_edge_ids(&graph, &[0, 1]));

        assert!(!full_subject.cycle_requirements().is_empty());
        assert!(cut_side_subject.cycle_requirements().is_empty());

        Ok(())
    }

    #[test]
    fn cycle_signature_rejects_antiparticles() {
        let model = minimal_model();
        assert!(CycleSignature::parse("[(-3)]", &model).is_err());
        assert!(CycleSignature::parse("[(e+)]", &model).is_err());
    }

    #[test]
    fn cycle_category_requirement_matches_particles_directly() {
        let model = minimal_model();
        let gluon = model.try_get_particle_from_pdg(21).unwrap();
        let strange = model.try_get_particle_from_pdg(3).unwrap();
        let requirement =
            CycleRequirement::new(vec![CycleMatcher::Pdg(3), CycleMatcher::Pdg(21)]).unwrap();
        assert!(cycle_matches_requirement(
            &[strange.clone(), gluon.clone()],
            &requirement
        ));
        let fermion_requirement = CycleRequirement::new(vec![CycleMatcher::Fermion]).unwrap();
        assert!(cycle_matches_requirement(&[strange], &fermion_requirement));
        assert!(!cycle_matches_requirement(&[gluon], &fermion_requirement));
    }

    #[test]
    fn subgraph_vertex_signature_counts_only_vertices_in_view() -> Result<()> {
        let graph = graph_with_vertex_rules("g0", &["V_A", "V_B", "V_C"])?;
        let subject = GraphSelectionSubject::subgraph(&graph, subgraph_from_edge_ids(&graph, &[0]));
        let counts = subject.vertex_rule_name_counts();

        assert_eq!(counts.values().copied().sum::<usize>(), 2);
        assert_eq!(
            GraphSelectionSubject::whole_graph(&graph)
                .vertex_rule_name_counts()
                .values()
                .copied()
                .sum::<usize>(),
            3
        );

        Ok(())
    }

    #[test]
    fn subgraph_particle_signature_counts_only_particles_in_view() -> Result<()> {
        let graph = particle_test_graph("g0")?;
        let subject = GraphSelectionSubject::subgraph(&graph, subgraph_from_edge_ids(&graph, &[0]));
        let scalar_0 = scalar_model().try_get_particle("scalar_0")?.pdg_code.abs();
        let scalar_1 = scalar_model().try_get_particle("scalar_1")?.pdg_code.abs();

        assert!(subject.particle_pdgs().contains(&scalar_0));
        assert!(!subject.particle_pdgs().contains(&scalar_1));
        assert!(
            GraphSelectionSubject::whole_graph(&graph)
                .particle_pdgs()
                .contains(&scalar_1)
        );

        Ok(())
    }

    #[test]
    fn selection_spec_combines_rules_with_and_and_respects_vertex_multiplicity() -> Result<()> {
        let mut graphs = vec![
            graph_with_vertex_rules("g0", &["V_A", "V_A", "V_B"])?,
            graph_with_vertex_rules("g1", &["V_A", "V_B", "V_C"])?,
            graph_with_vertex_rules("g2", &["V_A", "V_A", "V_C"])?,
        ];
        let graph_group_structure = complete_group_parsing(&mut graphs)?;

        let spec = GraphGroupSelectionSpec::new()
            .with_vertex_signatures(
                SelectionPolarity::With,
                vec![VertexSignature::parse("[V_A,V_A]")?],
            )
            .with_vertex_signatures(
                SelectionPolarity::Without,
                vec![VertexSignature::parse("[V_C]")?],
            );
        let plan = spec.plan(&graph_group_structure, |graph_id| graphs.get(graph_id))?;

        assert_eq!(plan.retained_group_ids(), &[GroupId(0)]);
        assert_eq!(plan.report().kept_master_graphs, vec!["g0".to_string()]);
        assert_eq!(
            plan.report().removed_master_graphs,
            vec!["g1".to_string(), "g2".to_string()]
        );

        Ok(())
    }

    #[test]
    fn graph_name_selection_supports_without_polarity() -> Result<()> {
        let mut graphs = vec![
            graph_with_vertex_rules("g0", &["V_A"])?,
            graph_with_vertex_rules("g1", &["V_B"])?,
            graph_with_vertex_rules("g2", &["V_C"])?,
        ];
        let graph_group_structure = complete_group_parsing(&mut graphs)?;

        let spec = GraphGroupSelectionSpec::new()
            .with_master_graph_names_polarity(SelectionPolarity::Without, vec!["g1".to_string()]);
        let plan = spec.plan(&graph_group_structure, |graph_id| graphs.get(graph_id))?;

        assert_eq!(plan.report().kept_master_graphs, vec!["g0", "g2"]);
        assert_eq!(plan.report().removed_master_graphs, vec!["g1"]);

        Ok(())
    }

    #[test]
    fn with_graph_names_are_authoritative_over_vetoes() -> Result<()> {
        let mut graphs = vec![
            graph_with_vertex_rules("g0", &["V_A"])?,
            graph_with_vertex_rules("g1", &["V_B"])?,
        ];
        let graph_group_structure = complete_group_parsing(&mut graphs)?;

        let spec = GraphGroupSelectionSpec::new()
            .with_master_graph_names(vec!["g1".to_string()])
            .with_vertex_signatures(
                SelectionPolarity::Without,
                vec![VertexSignature::parse("[V_B]")?],
            );
        let plan = spec.plan(&graph_group_structure, |graph_id| graphs.get(graph_id))?;

        assert_eq!(plan.report().kept_master_graphs, vec!["g0", "g1"]);
        assert!(plan.report().removed_master_graphs.is_empty());

        Ok(())
    }

    #[test]
    fn particle_selection_matches_required_sets() -> Result<()> {
        let mut graphs = vec![
            particle_test_graph("g0")?,
            graph_with_vertex_rules("g1", &["V_A"])?,
        ];
        let graph_group_structure = complete_group_parsing(&mut graphs)?;
        let scalar_0 = scalar_model()
            .try_get_particle("scalar_0")?
            .name
            .to_string();
        let scalar_1 = scalar_model()
            .try_get_particle("scalar_1")?
            .name
            .to_string();
        let scalar_2 = scalar_model()
            .try_get_particle("scalar_2")?
            .name
            .to_string();

        let spec = GraphGroupSelectionSpec::new()
            .with_particle_signatures(
                SelectionPolarity::With,
                vec![ParticleSignature::parse(
                    &format!("[{scalar_0},{scalar_1}]"),
                    scalar_model(),
                )?],
            )
            .with_particle_signatures(
                SelectionPolarity::Without,
                vec![ParticleSignature::parse(
                    &format!("[{scalar_2}]"),
                    scalar_model(),
                )?],
            );

        let err = spec
            .plan(&graph_group_structure, |graph_id| graphs.get(graph_id))
            .unwrap_err();
        assert!(
            err.to_string()
                .contains("Graph-group selection would remove all graph groups"),
            "{err:?}"
        );

        let spec = GraphGroupSelectionSpec::new().with_particle_signatures(
            SelectionPolarity::With,
            vec![ParticleSignature::parse(
                &format!("({scalar_0},{scalar_1})"),
                scalar_model(),
            )?],
        );
        let plan = spec.plan(&graph_group_structure, |graph_id| graphs.get(graph_id))?;
        assert_eq!(plan.report().kept_master_graphs, vec!["g0"]);

        Ok(())
    }

    #[test]
    fn structural_selection_matches_any_analysis_subject_and_without_vetoes() -> Result<()> {
        let mut graphs = vec![
            graph_with_vertex_rules("g0", &["V_A", "V_B", "V_C"])?,
            graph_with_vertex_rules("g1", &["V_D", "V_E", "V_F"])?,
        ];
        let graph_group_structure = complete_group_parsing(&mut graphs)?;
        let with_spec = GraphGroupSelectionSpec::new()
            .with_master_graph_names(vec!["g0".to_string(), "g1".to_string()])
            .with_vertex_signatures(
                SelectionPolarity::With,
                vec![VertexSignature::parse("[V_A]")?],
            );
        let without_spec = GraphGroupSelectionSpec::new().with_vertex_signatures(
            SelectionPolarity::Without,
            vec![VertexSignature::parse("[V_A]")?],
        );

        let with_plan = with_spec.plan_with_analysis_contexts(
            &graph_group_structure,
            |graph_id| graphs.get(graph_id),
            vertex_test_subjects_by_graph_id,
            |_graph_id, _graph| Ok(Vec::new()),
            "no subjects",
            "no cuts",
        )?;
        assert_eq!(with_plan.report().kept_master_graphs, vec!["g0", "g1"]);

        let without_plan = without_spec.plan_with_analysis_contexts(
            &graph_group_structure,
            |graph_id| graphs.get(graph_id),
            vertex_test_subjects_by_graph_id,
            |_graph_id, _graph| Ok(Vec::new()),
            "no subjects",
            "no cuts",
        )?;
        assert_eq!(without_plan.report().kept_master_graphs, vec!["g1"]);

        Ok(())
    }

    #[test]
    fn raised_propagator_selection_supports_any_raising_wildcard() -> Result<()> {
        let mut graphs = vec![raised_test_graph()?, raised_test_graph()?];
        graphs[0].name = "g0".to_string();
        graphs[1].name = "g1".to_string();
        let graph_group_structure = complete_group_parsing(&mut graphs)?;
        let any = RaisedPropagatorSignature::from_str("ANY_RAISING")?;
        let with_spec = GraphGroupSelectionSpec::new().with_raised_propagator_signatures(
            SelectionPolarity::With,
            RaisedPropagatorScope::All,
            vec![any.clone()],
        );
        let without_spec = GraphGroupSelectionSpec::new().with_raised_propagator_signatures(
            SelectionPolarity::Without,
            RaisedPropagatorScope::All,
            vec![any],
        );

        let with_plan = with_spec.plan_with_analysis_contexts(
            &graph_group_structure,
            |graph_id| graphs.get(graph_id),
            raised_any_test_subjects_by_graph_id,
            |_graph_id, _graph| Ok(Vec::new()),
            "no graph subjects",
            "no cuts",
        )?;
        assert_eq!(with_plan.report().kept_master_graphs, vec!["g0"]);

        let without_plan = without_spec.plan_with_analysis_contexts(
            &graph_group_structure,
            |graph_id| graphs.get(graph_id),
            raised_any_test_subjects_by_graph_id,
            |_graph_id, _graph| Ok(Vec::new()),
            "no graph subjects",
            "no cuts",
        )?;
        assert_eq!(without_plan.report().kept_master_graphs, vec!["g1"]);

        Ok(())
    }

    #[test]
    fn raised_cut_selection_matches_any_cut_and_without_vetoes() -> Result<()> {
        let mut graphs = vec![raised_test_graph()?, raised_test_graph()?];
        graphs[0].name = "g0".to_string();
        graphs[1].name = "g1".to_string();
        let graph_group_structure = complete_group_parsing(&mut graphs)?;
        let with_spec = GraphGroupSelectionSpec::new().with_raised_cut_signatures(
            SelectionPolarity::With,
            RaisedPropagatorScope::All,
            vec![RaisedPropagatorSignature::from_str("[2]")?],
        );
        let without_spec = GraphGroupSelectionSpec::new().with_raised_cut_signatures(
            SelectionPolarity::Without,
            RaisedPropagatorScope::All,
            vec![RaisedPropagatorSignature::from_str("[2]")?],
        );

        let with_plan = with_spec.plan_with_analysis_contexts(
            &graph_group_structure,
            |graph_id| graphs.get(graph_id),
            |_graph_id, graph| Ok(vec![GraphSelectionSubject::whole_graph(graph)]),
            raised_cut_test_subjects_by_graph_id,
            "no graph subjects",
            "no cuts",
        )?;
        assert_eq!(with_plan.report().kept_master_graphs, vec!["g0"]);

        let without_plan = without_spec.plan_with_analysis_contexts(
            &graph_group_structure,
            |graph_id| graphs.get(graph_id),
            |_graph_id, graph| Ok(vec![GraphSelectionSubject::whole_graph(graph)]),
            raised_cut_test_subjects_by_graph_id,
            "no graph subjects",
            "no cuts",
        )?;
        assert_eq!(without_plan.report().kept_master_graphs, vec!["g1"]);

        Ok(())
    }

    #[test]
    fn raised_cut_selection_supports_any_raising_wildcard() -> Result<()> {
        let mut graphs = vec![raised_test_graph()?, raised_test_graph()?];
        graphs[0].name = "g0".to_string();
        graphs[1].name = "g1".to_string();
        let graph_group_structure = complete_group_parsing(&mut graphs)?;
        let any = RaisedPropagatorSignature::from_str("ANY_RAISING")?;
        let with_spec = GraphGroupSelectionSpec::new().with_raised_cut_signatures(
            SelectionPolarity::With,
            RaisedPropagatorScope::All,
            vec![any.clone()],
        );
        let without_spec = GraphGroupSelectionSpec::new().with_raised_cut_signatures(
            SelectionPolarity::Without,
            RaisedPropagatorScope::All,
            vec![any],
        );

        let with_plan = with_spec.plan_with_analysis_contexts(
            &graph_group_structure,
            |graph_id| graphs.get(graph_id),
            |_graph_id, graph| Ok(vec![GraphSelectionSubject::whole_graph(graph)]),
            raised_cut_test_subjects_by_graph_id,
            "no graph subjects",
            "no cuts",
        )?;
        assert_eq!(with_plan.report().kept_master_graphs, vec!["g0"]);

        let without_plan = without_spec.plan_with_analysis_contexts(
            &graph_group_structure,
            |graph_id| graphs.get(graph_id),
            |_graph_id, graph| Ok(vec![GraphSelectionSubject::whole_graph(graph)]),
            raised_cut_test_subjects_by_graph_id,
            "no graph subjects",
            "no cuts",
        )?;
        assert_eq!(without_plan.report().kept_master_graphs, vec!["g1"]);

        Ok(())
    }
}
