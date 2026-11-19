//! Assemble the ordinary UV preparation and the UV form of the IR correction.

use std::{
    collections::BTreeMap,
    sync::Arc,
    time::{Duration, Instant},
};

use eyre::{Result, WrapErr, ensure, eyre};
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{SuBitGraph, SubSetLike},
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    parse, symbol,
};
use vakint::{Vakint, VakintSettings};

use crate::uv::approx::integrated::Integrated;
use crate::{
    debug_tags,
    graph::Graph,
    model::Model,
    utils::GS,
    uv::{UVgenerationSettings, UltravioletGraph, uv_graph::UVE},
};

use super::{
    ChildMass, CountertermCall, IrCorrection, IrForestTerm, MassiveUvEvaluator, RStarPlan,
    UvChildrenRoute, UvPreparationTerm,
    boundary::{self, BoundaryMode, BoundaryRequest, Built, ClosedRootNormalization},
    expression::{DiracBoundary, DiracKind, IntegrationTimings, RStarIntegrand, RootProjector},
    ir_dirac::IrDiracMode,
    moments::hard_moments,
    physical::PreparedPhysicalUvTerm,
    top_mass::{self, RStarTerms, TopMassTerms},
};

/// The `(coefficient, insertion)` moments of one IR correction term's hard
/// part, whether the closed path computed them, and the cost of today's
/// path (`IR_CLOSED`, `IR_HARD_COST`).
struct TermMoments {
    moments: Vec<(Atom, Atom)>,
    use_closed: bool,
    taylor_seconds: f64,
    taylor_bytes: usize,
    project_seconds: f64,
    project_bytes: usize,
}

/// `GAMMALOOP_IR_HARD_KEYS` (any value): log the canonical identities of the
/// IR hard components (`IR_HARD`, `IR_HARD_PLAN`) and their costs
/// (`IR_HARD_COST`). Diagnostics only.
fn ir_hard_keys() -> bool {
    static KEYS: std::sync::LazyLock<bool> =
        std::sync::LazyLock::new(|| std::env::var_os("GAMMALOOP_IR_HARD_KEYS").is_some());
    *KEYS
}

/// `GAMMALOOP_IR_HARD_CENSUS=1|only` (diagnostics): one `IR_HARD_CENSUS`
/// line per IR correction of every R* plan evaluated; `only` also skips the
/// bare-term export of `pint6_dump` after the census. `Some(only)` when set.
fn ir_hard_census() -> Option<bool> {
    static MODE: std::sync::LazyLock<Option<bool>> = std::sync::LazyLock::new(|| {
        match std::env::var("GAMMALOOP_IR_HARD_CENSUS")
            .as_deref()
            .map(str::trim)
        {
            Ok("1") => Some(false),
            Ok("only") => Some(true),
            _ => None,
        }
    });
    *MODE
}

/// `GAMMALOOP_IR_HARD_DECLARE` (any value): one `IR_HARD_DECL` line per IR
/// hard term on today's path (diagnostics).
fn ir_hard_declare() -> bool {
    static ON: std::sync::LazyLock<bool> =
        std::sync::LazyLock::new(|| std::env::var_os("GAMMALOOP_IR_HARD_DECLARE").is_some());
    *ON
}

/// `GAMMALOOP_IR_HARD_DRYRUN` (any value): plan only. Every R* counterterm is
/// returned as zero after its `IR_HARD` lines; results are meaningless.
fn ir_hard_dryrun() -> bool {
    static DRY: std::sync::LazyLock<bool> =
        std::sync::LazyLock::new(|| std::env::var_os("GAMMALOOP_IR_HARD_DRYRUN").is_some());
    *DRY
}

impl RStarPlan {
    /// Complete every original proper-UV child and prepare its physical
    /// quotient without reducing or integrating the remaining loop tensors.
    ///
    /// The optional linear root projector acts on each complete quotient
    /// numerator exactly once. Its coefficients must be independent of the
    /// remaining loop momenta. Child Laurent products and their correlations
    /// remain inside the terms, including when an independently grouped bare
    /// amplitude vanishes. Group matching terms with [`super::PhysicalUvGroup`]
    /// and apply the final `-K` to the resulting preparation plus bare poles.
    pub fn prepare_physical_uv_terms(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_projector: Option<&RootProjector<'_>>,
    ) -> Result<Vec<PreparedPhysicalUvTerm>> {
        let evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            ChildMass::Zero,
        )?;
        let started = Instant::now();
        let settings: Arc<str> = serde_json::to_string(&uv_settings.vakint)?.into();
        let terms = self
            .ordinary_uv
            .iter()
            .filter(|term| !term.counterterms.is_empty())
            .map(|term| {
                PreparedPhysicalUvTerm::new(
                    evaluation.physical_remainder(term)?,
                    model,
                    Arc::clone(&settings),
                    root_projector,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        tracing::info!(
            uv = true, rstar = true, timing = true,
            stage = "rstar_prepare_physical_quotients",
            graph = %graph.name,
            prepared_quotients = terms.len(),
            proper_uv_seconds = evaluation.proper_uv_elapsed.as_secs_f64(),
            preparation_seconds = started.elapsed().as_secs_f64(),
            "Prepared original-forest physical quotients before integration"
        );
        Ok(terms)
    }

    /// Compute the local MS counterterm of a homogeneous massless graph.
    ///
    /// The result retains its open tensor indices and external edge momenta.
    /// Global graph prefactors and external projectors are not included.
    /// All proper UV maps use the existing massive R operation; only the
    /// remaining single-mass integrals use Vakint's p-integral backend. Include
    /// `forcer` in the evaluation methods for that backend.
    pub fn counterterm(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            ChildMass::Reference,
        )?
        .counterterm()
    }

    /// Compute a root-projected local MS counterterm using the same R* forest.
    ///
    /// The linear projector acts on each complete outer Taylor polynomial,
    /// before its remaining vacuum tensor reduction. It is never applied to
    /// proper UV children or decorated IR quotient counterterms, whose open
    /// tensor boundaries are retained. Its coefficients must be regular in
    /// epsilon and it must not introduce loop momenta. Rational functions of
    /// external momenta are allowed when their denominators are scalar and
    /// independent of every loop momentum. These factors are spectators:
    /// the root Taylor operation has already finished.
    ///
    /// Root normalization factors such as `1/(d-1)` participate in the final
    /// Laurent projection. Scalar marker coefficients may carry several root
    /// projection channels through a single evaluation. Graph-wide prefactors
    /// remain excluded, as for [`Self::counterterm`].
    pub fn projected_counterterm(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_projector: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_projector),
            ChildMass::Reference,
        )?
        .counterterm()
    }

    /// Local export patch: run only the ordinary-UV preparation terms of
    /// [`Self::projected_counterterm`] that carry no proper UV counterterm
    /// (the empty spinney), through the same root projector, outer Taylor
    /// expansion and single-mass integration. Proper UV children and IR
    /// corrections are not computed. Returns the number of terms integrated;
    /// the result values are discarded (meaningful only with requests skipped
    /// through `PINT6_SKIP_MIN_P_LOOPS`, for request export).
    pub fn projected_bare_terms_only(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_projector: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<usize> {
        self.bare_terms_only_with(graph, model, vakint, uv_settings, Some(root_projector))
    }

    /// As [`Self::projected_counterterm`] for a metric root projector whose
    /// momentum-free metric is applied in the root scope instead: every
    /// root-level remainder is multiplied by [`RStarIntegrand::root_closure`]
    /// together with its child insertions, so the remainder algebra contracts
    /// a closed scalar. IR quotient counterterms and proper UV children keep
    /// their open indices; IR corrections are closed after their quotient `K`.
    /// `root_normalization` then supplies only the momentum-dependent rest of
    /// the projector (normalization, rerouting, channel tags) after Taylor.
    pub fn projected_counterterm_closing_root(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            ChildMass::Reference,
        )?
        .closing_root()?
        .counterterm()
    }

    /// As [`Self::projected_counterterm_from_pole_part`], closing the root
    /// metric in the root scope (see [`Self::projected_counterterm_closing_root`]).
    /// `bare_poles` keep their open slots; they are closed where the metric
    /// projector applied before, after their own pole projection.
    pub fn projected_counterterm_from_pole_part_closing_root(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        bare_poles: &Atom,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            ChildMass::Zero,
        )?
        .closing_root()?
        .counterterm_from_pole_part(bare_poles)
    }

    /// Local export patch: as [`Self::projected_bare_terms_only`], closing the
    /// metric root in the root scope (see [`Self::projected_counterterm_closing_root`]).
    pub fn projected_bare_terms_only_closing_root(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<usize> {
        self.bare_terms_only_closing(
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            true,
            None,
        )
    }

    /// [`Self::projected_counterterm_closing_root`] with the closed metric-root
    /// normalisation also given as data (`closed_root`: the callback is
    /// `closed_root.factor · σ` on closed scalars). With
    /// `GAMMALOOP_REMAINDER_BOUNDARY` unset this is
    /// [`Self::projected_counterterm_closing_root`]; otherwise the admitted
    /// root-scope remainders take the remainder-boundary path.
    pub fn projected_counterterm_closing_root_scalar(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
        closed_root: &ClosedRootNormalization,
    ) -> Result<Atom> {
        let mut evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            ChildMass::Reference,
        )?
        .closing_root()?;
        evaluation.closed_root = Some(closed_root);
        evaluation.counterterm()
    }

    /// [`Self::projected_bare_terms_only_closing_root`] with the closed
    /// metric-root normalisation as data (see
    /// [`Self::projected_counterterm_closing_root_scalar`]).
    pub fn projected_bare_terms_only_closing_root_scalar(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
        closed_root: &ClosedRootNormalization,
    ) -> Result<usize> {
        self.bare_terms_only_closing(
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            true,
            Some(closed_root),
        )
    }

    /// Local export patch: as [`Self::projected_bare_terms_only`], without a
    /// root projector, i.e. the bare terms of [`Self::pole_part`] (the route
    /// of insertion surrogates above four loops).
    pub fn bare_terms_only_unprojected(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<usize> {
        self.bare_terms_only_with(graph, model, vakint, uv_settings, None)
    }

    fn bare_terms_only_with<'a>(
        &'a self,
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        uv_settings: &'a UVgenerationSettings,
        root_projector: Option<&'a RootProjector<'a>>,
    ) -> Result<usize> {
        self.bare_terms_only_closing(
            graph,
            model,
            vakint,
            uv_settings,
            root_projector,
            false,
            None,
        )
    }

    fn bare_terms_only_closing<'a>(
        &'a self,
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        uv_settings: &'a UVgenerationSettings,
        root_projector: Option<&'a RootProjector<'a>>,
        close_root: bool,
        closed_root: Option<&'a ClosedRootNormalization>,
    ) -> Result<usize> {
        let mut evaluation =
            // Bare terms only: no proper UV child is computed.
            Evaluation::new_with(
                self,
                graph,
                model,
                vakint,
                uv_settings,
                root_projector,
                ChildMass::Reference,
                false,
            )?;
        if close_root {
            evaluation = evaluation.closing_root()?;
        }
        evaluation.closed_root = closed_root;
        // GAMMALOOP_IR_HARD_CENSUS: the plan's IR corrections (`only`: and
        // nothing else; the bare terms are not integrated).
        if let Some(only) = ir_hard_census() {
            evaluation.log_ir_hard_census();
            if only {
                return Ok(0);
            }
        }
        let mut count = 0;
        for term in &self.ordinary_uv {
            if !term.counterterms.is_empty() || term.remainder.is_scaleless(graph, self.marked_edge)
            {
                continue;
            }
            super::pint6_dump::set_term(format!(
                "graph={} marked_edge={} term=ordinary_uv[empty spinney] remaining_loops={} taylor_order={}",
                graph.name,
                usize::from(self.marked_edge),
                term.remainder.loops,
                self.taylor_order
            ));
            if let Err(error) = evaluation.integrate_ordinary(term) {
                // The counts up to the failure, before the error propagates.
                evaluation.log_boundary_summary("error");
                return Err(error);
            }
            count += 1;
        }
        evaluation.log_boundary_summary("ok");
        Ok(count)
    }

    /// The local MS counterterm of a logarithmically divergent graph whose
    /// open slots are four Minkowski slots of one representation, through its
    /// three scalar metric channels instead of the open rank-four tensor.
    ///
    /// The counterterm is `sum_j K[P_j . X] T_j` over the metric pairings
    /// `T_j` ([`RStarIntegrand::isotropic_boundary`]) with the dual projectors
    /// `P_j = ((d+2) T_j - sum_i T_i) / (d (d-1) (d+2))`: MS acts on the
    /// coefficients of a d-dimensional basis, which these projections are.
    /// The momentum-free polynomial part of `sum_j tau_j P_j` closes every
    /// root-scope remainder (as [`Self::projected_counterterm_closing_root`]);
    /// its rational normalization joins after the Taylor step, inside this
    /// graph's `K`. Nested children and quotients keep their open slots. The
    /// channel markers `tau_j` are then replaced by `T_j`.
    pub fn isotropic_counterterm(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        pairings: &[Atom; 3],
    ) -> Result<Atom> {
        let d = Atom::var(GS.dim);
        let markers = [
            symbol!("gammalooprs::uv_isotropic_channel_1"; Scalar),
            symbol!("gammalooprs::uv_isotropic_channel_2"; Scalar),
            symbol!("gammalooprs::uv_isotropic_channel_3"; Scalar),
        ];
        // sum_j tau_j ((d+2) T_j - sum_i T_i) = sum_j ((d+2) tau_j - sum_i tau_i) T_j.
        let markers_sum = markers
            .iter()
            .fold(Atom::Zero, |sum, marker| sum + Atom::var(*marker));
        let closure =
            markers
                .iter()
                .zip(pairings)
                .fold(Atom::Zero, |closure, (marker, pairing)| {
                    closure + ((&d + 2) * Atom::var(*marker) - &markers_sum) * pairing
                });
        let normalization = Atom::one() / (&d * (&d - 1) * (&d + 2));
        let normalize = |tensor: &Atom| -> Result<Atom> { Ok(tensor * &normalization) };
        let mut evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(&normalize),
            ChildMass::Reference,
        )?;
        evaluation.root_closure = Some(closure);
        evaluation.channel = ("isotropic", "1/(d (d-1) (d+2))");
        let channels = evaluation.counterterm()?;
        Ok(markers
            .iter()
            .zip(pairings)
            .fold(channels, |tensor, (marker, pairing)| {
                tensor.replace(Atom::var(*marker)).with(pairing.clone())
            }))
    }

    /// The local MS counterterm of a graph with one open Dirac line whose
    /// counterterm has a unique Dirac structure
    /// ([`RStarIntegrand::dirac_boundary`]), as a closed trace instead of the
    /// open matrix.
    ///
    /// A degree-0 vertex is `c gamma^mu` with `c = Tr[gamma_mu X^mu] /
    /// (Tr[1] d)`; a degree-1 self-energy is `V_nu gamma^nu` with `V_nu =
    /// Tr[gamma_nu X] / Tr[1]` (linear in the boundary momenta). MS acts on
    /// the coefficient of the basis element, which this projection is. As in
    /// [`Self::isotropic_counterterm`], the momentum-free `gamma` closure
    /// closes every root-scope remainder, the rational normalization joins
    /// after the Taylor step inside this graph's `K`, and nested children
    /// and quotients keep their open slots. The open structure is then
    /// multiplied back.
    pub fn dirac_counterterm(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        boundary: &DiracBoundary,
    ) -> Result<Atom> {
        let normalize = |tensor: &Atom| -> Result<Atom> { Ok(tensor * &boundary.normalization) };
        let mut evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(&normalize),
            ChildMass::Reference,
        )?;
        evaluation.root_closure = Some(boundary.closure.clone());
        evaluation.channel = match boundary.kind {
            DiracKind::Vertex => ("dirac vertex", "1/(Tr[1] d)"),
            DiracKind::SelfEnergy => ("dirac self-energy", "1/Tr[1]"),
        };
        let trace = evaluation.counterterm()?;
        Ok(trace * &boundary.structure)
    }

    /// Reconstruct all poles of the bare, infrared-finite massless propagator
    /// graph at its original off-shell momentum, including subdivergence poles.
    ///
    /// From `Delta_U(G) = -K R'_U G`, this returns
    /// `K G = -Delta_U(G) - K sum_{S != empty} Delta_U(S) G/S`.
    /// Proper UV maps use the existing massive R operation; their contracted
    /// remainders are lower-loop p-integrals, with no Taylor expansion or
    /// auxiliary mass. Vakint must support these lower-loop integrals and the
    /// extra epsilon orders required by the insertion poles. Physical IR
    /// divergences in the input are not reconstructed by this identity.
    ///
    /// As for `counterterm`, graph prefactors and external projectors are
    /// excluded; open tensor indices and external momenta are retained.
    pub fn pole_part(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<Atom> {
        let (_, external_routing) = RStarIntegrand::physical_external_routing(
            graph,
            &graph.external_filter::<SuBitGraph>(),
        )?;
        let mut evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            ChildMass::Reference,
        )?;
        let result = -evaluation.counterterm()? - evaluation.physical_uv_preparation()?;
        evaluation.poles(&super::linearize::linearize(
            result.replace_multiple(&external_routing),
        ))
    }

    /// Recover the original local MS counterterm from independently known bare
    /// poles: `Delta_U(G) = -K [K G + sum_{S != empty} Delta_U(S) G/S]`.
    ///
    /// This retains the original graph's proper UV forest. It therefore also
    /// accepts bare poles obtained by integrating propagator insertions or
    /// replacing them with normalized scalar bubbles: those substitutions need
    /// not preserve the original subdivergences or their nested pole maps.
    ///
    /// `bare_poles` must describe this infrared-finite massless two-point graph,
    /// in the supplied Vakint normalization, with its original open tensor
    /// slots and external momentum labels. Exclude graph-wide prefactors and
    /// external projectors, as in [`Self::pole_part`]. Boundary momenta are
    /// identified using the same physical routing as that method. All proper
    /// UV counterterms belong to the original graph; the root R*
    /// counterterm is not evaluated.
    pub fn counterterm_from_pole_part(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        bare_poles: &Atom,
    ) -> Result<Atom> {
        // Only the physical remainders consume the children, at zero mass.
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            ChildMass::Zero,
        )?
        .counterterm_from_pole_part(bare_poles)
    }

    /// Project the complete physical preparation before its final pole map.
    ///
    /// This is [`Self::counterterm_from_pole_part`] with a linear root projector
    /// whose coefficients are regular in epsilon and contain no loop momenta.
    /// Proper UV children retain their open tensor boundaries. Contractions
    /// and scalar dimensions introduced by the projector are completed before
    /// `K`, including normalization factors such as `1/(d-1)`.
    pub fn projected_counterterm_from_pole_part(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        bare_poles: &Atom,
        root_projector: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_projector),
            ChildMass::Zero,
        )?
        .counterterm_from_pole_part(bare_poles)
    }
}

#[cfg(test)]
impl RStarPlan {
    /// Tests: the outer Taylor polynomial of every ordinary-UV remainder of a
    /// closing-root evaluation (children included), formed exactly as
    /// [`RStarIntegrand::integrate`] forms it, with the remainder's loop
    /// count, its contracted expression and the idenso certificate of that
    /// contraction ([`RStarIntegrand::contracted_certified`]). Nothing is
    /// projected or integrated.
    pub(super) fn closing_root_outer_taylor_polynomials(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<Vec<(usize, Atom, bool, Atom)>> {
        let evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            ChildMass::Reference,
        )?
        .closing_root()?;
        let mut polynomials = Vec::new();
        for term in &self.ordinary_uv {
            if term.remainder.is_scaleless(graph, self.marked_edge) {
                continue;
            }
            let remainder = evaluation.remainder(term)?;
            let (expression, certified) = remainder.contracted_certified()?;
            let polynomial = remainder.outer_taylor(
                expression,
                &evaluation.root_edges,
                usize::try_from(self.taylor_order)?,
                evaluation.reference_mass_squared,
                &evaluation.auxiliary_mass_squared,
                &mut super::expression::IntegrationTimings::new(),
            )?;
            polynomials.push((
                remainder.region.loops,
                expression.clone(),
                *certified,
                polynomial,
            ));
        }
        Ok(polynomials)
    }
}

struct Evaluation<'a> {
    plan: &'a RStarPlan,
    graph: &'a Graph,
    model: &'a Model,
    vakint: &'a Vakint,
    uv_settings: &'a UVgenerationSettings,
    settings: VakintSettings,
    epsilon: Symbol,
    reference_mass_squared: Symbol,
    auxiliary_mass_squared: Atom,
    children: BTreeMap<SuBitGraph, Atom>,
    /// Carriers scaled by the outer Taylor operation: the graph's external
    /// edges, or the full crown of a proper support root.
    root_edges: Vec<EdgeIndex>,
    uv: MassiveUvEvaluator<'a>,
    quotient_counterterms: BTreeMap<(CountertermCall, Atom), Atom>,
    integration_id: usize,
    root_projector: Option<&'a RootProjector<'a>>,
    /// The root metric closure of a closing evaluation (root scope only).
    root_closure: Option<Atom>,
    /// The root projection channel and its normalisation, recorded in the
    /// closure dump headers (`PINT6_CLOSURE_DUMP`).
    channel: (&'static str, &'static str),
    /// The closed metric-root normalisation as data, given by the `_scalar`
    /// entry points; used only when `GAMMALOOP_REMAINDER_BOUNDARY` is set.
    closed_root: Option<&'a ClosedRootNormalization>,
    /// Time on the remainder-boundary path (inert when the mode is off).
    boundary_elapsed: Duration,
    /// Remainder-boundary counts of this evaluation (`BOUNDARY_SUMMARY`).
    summary: boundary::BoundarySummary,
    started: Instant,
    proper_uv_elapsed: Duration,
    quotient_elapsed: Duration,
    integration_elapsed: Duration,
    /// `GAMMALOOP_RSTAR_TOP_MASS_TERMS` of a top-level root (`None` for a
    /// child root): the scope of every outer integration of this evaluation.
    top_mass: Option<TopMassTerms>,
    /// `GAMMALOOP_RSTAR_TERMS` of a top-level root (`All` for a child root).
    rstar_terms: RStarTerms,
    /// Root terms skipped by `rstar_terms` (logged when nonzero).
    skipped_terms: usize,
    /// `GAMMALOOP_IR_DIRAC_CLOSURE`, read when the evaluation starts.
    ir_dirac: IrDiracMode,
    /// The IR Dirac closure of each IR correction (or its decline reason),
    /// computed on first use when `ir_dirac` is on.
    ir_dirac_boundaries: BTreeMap<usize, std::result::Result<DiracBoundary, String>>,
}

impl<'a> Evaluation<'a> {
    fn new(
        plan: &'a RStarPlan,
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        uv_settings: &'a UVgenerationSettings,
        root_projector: Option<&'a RootProjector<'a>>,
        child_mass: ChildMass,
    ) -> Result<Self> {
        Self::new_with(
            plan,
            graph,
            model,
            vakint,
            uv_settings,
            root_projector,
            child_mass,
            true,
        )
    }

    fn new_with(
        plan: &'a RStarPlan,
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        uv_settings: &'a UVgenerationSettings,
        root_projector: Option<&'a RootProjector<'a>>,
        child_mass: ChildMass,
        compute_children: bool,
    ) -> Result<Self> {
        let started = Instant::now();
        // The top-level switches act at a top-level root only; read before
        // any work so that a bad value fails early.
        let (top_mass, rstar_terms) = if plan.root_boundary.is_none() {
            let mode = match top_mass::top_mass_terms()? {
                // The odd orders vanish only for a projected (closed scalar)
                // root channel.
                TopMassTerms::DropOdd if root_projector.is_none() => TopMassTerms::Drop,
                mode => mode,
            };
            // The bare/nonbare split is a split of the projected diagram root
            // only, never of an unprojected propagator-insertion pole part.
            let terms = top_mass::rstar_terms()?;
            (
                Some(mode),
                if root_projector.is_some() {
                    terms
                } else {
                    RStarTerms::All
                },
            )
        } else {
            (None, RStarTerms::All)
        };
        for (_, _, edge) in graph.iter_edges_of(&plan.support) {
            ensure!(
                edge.data.mass_atom(model).is_zero(),
                "R* currently requires a massless input graph"
            );
        }
        let ir_dirac = super::ir_dirac::ir_dirac_mode()?;
        let settings = uv_settings.vakint.true_settings();
        let epsilon = parse!(&settings.epsilon_symbol)
            .get_symbol()
            .ok_or_else(|| eyre!("Vakint epsilon must be a symbol"))?;
        let reference_mass_squared = *super::REFERENCE_MASS_SQUARED;
        let reference_mass = Atom::var(reference_mass_squared).sqrt();
        let children_started = Instant::now();
        let children = if plan.counterterms.is_empty() || !compute_children {
            BTreeMap::new()
        } else {
            plan.evaluate_uv_children_for(
                graph,
                model,
                &reference_mass,
                vakint,
                uv_settings,
                child_mass,
            )
            .wrap_err(match plan.uv_children {
                UvChildrenRoute::Massive => {
                    "while computing the proper UV counterterms with massive R"
                }
                UvChildrenRoute::Rstar { .. } => {
                    "while computing the proper UV counterterms with recursive R*"
                }
            })?
        };
        let root_edges = match &plan.root_boundary {
            None => {
                let root_external = graph.external_filter::<SuBitGraph>();
                graph
                    .iter_edges_of(&root_external)
                    .map(|(_, edge, _)| edge)
                    .collect::<Vec<_>>()
            }
            Some(boundary) => graph
                .iter_edges_of(boundary)
                .map(|(_, edge, _)| edge)
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect(),
        };
        let proper_uv_elapsed = children_started.elapsed();
        let mut uv = MassiveUvEvaluator::new(graph, model, vakint, uv_settings);
        if let Some(cache) = &plan.counterterm_cache {
            uv = uv.with_counterterm_cache(cache.clone());
        }
        Ok(Self {
            plan,
            graph,
            model,
            vakint,
            uv_settings,
            settings,
            epsilon,
            reference_mass_squared,
            auxiliary_mass_squared: Atom::var(symbol!("rstar_mass_squared"; Scalar)),
            children,
            root_edges,
            uv,
            quotient_counterterms: BTreeMap::new(),
            integration_id: 0,
            root_projector,
            root_closure: None,
            channel: if root_projector.is_some() {
                (
                    "projection open-root",
                    "root projector after the outer Taylor step",
                )
            } else {
                ("none", "1")
            },
            closed_root: None,
            boundary_elapsed: Duration::ZERO,
            summary: Default::default(),
            started,
            proper_uv_elapsed,
            quotient_elapsed: Duration::ZERO,
            integration_elapsed: Duration::ZERO,
            top_mass,
            rstar_terms,
            skipped_terms: 0,
            ir_dirac,
            ir_dirac_boundaries: BTreeMap::new(),
        })
    }

    /// `GAMMALOOP_RSTAR_TERMS`: whether this evaluation computes a root term
    /// with (`nonbare_term`) or without proper UV children or IR counterterms.
    fn keeps_term(&mut self, nonbare_term: bool) -> bool {
        let keeps = self.rstar_terms.keeps(nonbare_term);
        if !keeps {
            self.skipped_terms += 1;
        }
        keeps
    }

    /// One line per top-level root evaluation whose switches are not at
    /// their defaults.
    fn log_top_level(&self) {
        if self.rstar_terms != RStarTerms::All
            || self.top_mass.is_some_and(|mode| mode != TopMassTerms::Keep)
        {
            eprintln!(
                "RSTAR_TOP graph={} marked_edge={} taylor_order={} rstar_terms={} top_mass_terms={} skipped_terms={}",
                self.graph.name,
                usize::from(self.plan.marked_edge),
                self.plan.taylor_order,
                self.rstar_terms.as_str(),
                self.top_mass.map_or("child", TopMassTerms::as_str),
                self.skipped_terms
            );
        }
    }

    /// Close the root metric in every root-scope remainder of this evaluation.
    fn closing_root(mut self) -> Result<Self> {
        self.root_closure = Some(RStarIntegrand::root_closure(
            self.graph,
            &self.plan.support,
        )?);
        self.channel = ("projection metric", "1/(i p^2 d_A (d-1))");
        Ok(self)
    }

    fn counterterm(&mut self) -> Result<Atom> {
        let result = self.counterterm_terms();
        if result.is_err() {
            // The counts up to the failure, before the error propagates
            // (log_timing prints them on success).
            self.log_boundary_summary("error");
        }
        result
    }

    fn counterterm_terms(&mut self) -> Result<Atom> {
        if self.graph.compute_dod(&self.plan.support) < 0 {
            self.log_timing();
            return Ok(Atom::Zero);
        }
        let mut prepared = Atom::Zero;
        if ir_hard_keys() {
            self.log_ir_hard_keys();
        }
        if ir_hard_census().is_some() {
            self.log_ir_hard_census();
        }
        if ir_hard_dryrun() {
            self.log_timing();
            return Ok(Atom::Zero);
        }
        for term in &self.plan.ordinary_uv {
            if term
                .remainder
                .is_scaleless(self.graph, self.plan.marked_edge)
            {
                continue;
            }
            // GAMMALOOP_RSTAR_TERMS: the bare term is the empty spinney.
            if !self.keeps_term(!term.counterterms.is_empty()) {
                continue;
            }
            prepared += self.integrate_ordinary(term)?;
        }
        let plan = self.plan;
        for (correction_index, correction) in plan.ir.iter().enumerate() {
            // GAMMALOOP_IR_HARD=child|child-check: the whole correction from
            // its integrated hard part H(γ); a decline takes today's terms.
            let child_mode = super::closed::child_mode();
            if child_mode != super::closed::ChildMode::Off
                && self.rstar_terms.keeps(true)
                && let Some(value) =
                    self.ir_correction_child(correction_index, correction, child_mode)?
            {
                prepared += value;
                continue;
            }
            for (term_index, term) in correction.preparation.iter().enumerate() {
                if term
                    .remainder
                    .is_scaleless(self.graph, self.plan.marked_edge)
                {
                    continue;
                }
                // GAMMALOOP_RSTAR_TERMS: every IR correction is a nonbare term.
                if !self.keeps_term(true) {
                    continue;
                }
                prepared += self.ir_term_today(correction_index, correction, term_index, term)?;
            }
        }
        let result = -self.poles(&prepared)?;
        self.log_top_level();
        self.log_timing();
        Ok(result)
    }

    /// One IR correction term on today's path (closed or projected moments,
    /// then the generic or analytic outer step): its contribution to the
    /// prepared sum.
    fn ir_term_today(
        &mut self,
        correction_index: usize,
        correction: &IrCorrection,
        term_index: usize,
        term: &UvPreparationTerm,
    ) -> Result<Atom> {
        let mut prepared = Atom::Zero;
        // The hard part feeds quotient counterterms (a `K` each): it
        // keeps its open root slots, closed after the quotient below.
        let hard = self.open_remainder(term)?;
        let order = usize::try_from(correction.taylor_order)?;
        // GAMMALOOP_IR_DIRAC_CLOSURE=on|check: the moments of the hard part
        // with its open Dirac line closed into a trace, reassembled.
        let (hard, today) = if self.ir_dirac == IrDiracMode::Off {
            let today =
                self.ir_term_moments(correction_index, correction, term_index, term, &hard, order)?;
            (hard, today)
        } else {
            self.ir_term_dirac(correction_index, correction, term_index, term, hard, order)?
        };
        let TermMoments {
            moments,
            use_closed,
            taylor_seconds,
            taylor_bytes,
            project_seconds,
            project_bytes,
        } = today;
        let downstream_started = Instant::now();
        let quotient_before = self.quotient_elapsed;
        let integration_before = self.integration_elapsed;
        // GAMMALOOP_IR_OUTER (any R* root; closed moments only).
        let outer = top_mass::ir_outer()?;
        let analytic = if outer == top_mass::IrOuter::Analytic || outer == top_mass::IrOuter::Check
        {
            let started = Instant::now();
            let outcome = if !use_closed {
                Err("not_closed".to_owned())
            } else {
                self.ir_outer_analytic(&hard, &correction.forest, &moments)?
            };
            Some((outcome, started.elapsed().as_secs_f64()))
        } else {
            None
        };
        let generic_started = Instant::now();
        let analytic_value = match &analytic {
            Some((Ok((value, _)), _)) if outer == top_mass::IrOuter::Analytic => {
                Some(value.clone())
            }
            _ => None,
        };
        let generic_value = if let Some(value) = analytic_value {
            prepared += value;
            None
        } else {
            let mut generic = Atom::Zero;
            for forest in &correction.forest {
                let mut expression = Atom::Zero;
                for (coefficient, insertion) in &moments {
                    expression += coefficient * self.ir_forest(forest, insertion)?;
                }
                if !expression.is_zero() {
                    if let Some(closure) = &self.root_closure {
                        expression *= closure;
                    }
                    generic += self.integrate(&hard, &expression)?;
                }
            }
            prepared += &generic;
            Some(generic)
        };
        let generic_seconds = generic_started.elapsed().as_secs_f64();
        if outer != top_mass::IrOuter::Unset {
            let (outcome, stats, analytic_seconds) = match &analytic {
                None => ("generic".to_owned(), String::new(), 0.0),
                Some((Ok((_, stats)), seconds)) => ("analytic".to_owned(), stats.clone(), *seconds),
                Some((Err(reason), seconds)) => {
                    (format!("declined:{reason}"), String::new(), *seconds)
                }
            };
            let check = match (&analytic, &generic_value) {
                (Some((Ok((value, _)), _)), Some(generic)) if outer == top_mass::IrOuter::Check => {
                    let difference = self.poles(&(generic - value))?;
                    let equal = difference.is_zero();
                    ensure!(
                        equal,
                        "GAMMALOOP_IR_OUTER=check: the analytic IR outer Taylor differs for {} corr={correction_index} term={term_index}: poles {difference}",
                        self.graph.name
                    );
                    " check=equal"
                }
                _ => "",
            };
            eprintln!(
                "IR_OUTER graph={} corr={correction_index} term={term_index} rem_loops={} order={order} moments={} forests={} outcome={outcome}{stats}{check} analytic_s={analytic_seconds:.3} generic_s={generic_seconds:.3}",
                self.graph.name,
                term.remainder.loops,
                moments.len(),
                correction.forest.len(),
            );
        }
        if ir_hard_keys() {
            eprintln!(
                "IR_HARD_COST graph={} corr={correction_index} term={term_index} rem_loops={} taylor_s={taylor_seconds:.3} taylor_bytes={taylor_bytes} project_s={project_seconds:.3} project_bytes={project_bytes} moments={} quotient_s={:.3} integrate_s={:.3} downstream_s={:.3}",
                self.graph.name,
                term.remainder.loops,
                moments.len(),
                (self.quotient_elapsed - quotient_before).as_secs_f64(),
                (self.integration_elapsed - integration_before).as_secs_f64(),
                downstream_started.elapsed().as_secs_f64(),
            );
        }
        Ok(prepared)
    }

    /// The `(coefficient, insertion)` moments of the hard part `hard` of one
    /// IR correction term at Taylor order `order`: closed
    /// (`GAMMALOOP_IR_HARD=closed|child`), else projected (today's path), with
    /// the `IR_CLOSED` line when the closed path is on.
    fn ir_term_moments(
        &self,
        correction_index: usize,
        correction: &IrCorrection,
        term_index: usize,
        term: &UvPreparationTerm,
        hard: &RStarIntegrand,
        order: usize,
    ) -> Result<TermMoments> {
        // GAMMALOOP_IR_HARD=closed|check (prototype, default off).
        let mode = super::closed::closed_mode();
        let mut closed_stats = super::closed::ClosedStats::default();
        let closed_started = Instant::now();
        let closed = if mode == super::closed::ClosedMode::Off {
            None
        } else {
            Some(
                hard.closed_moments(order, &mut closed_stats)
                    .wrap_err("while computing closed hard moments")?,
            )
        };
        let closed_seconds = closed_started.elapsed().as_secs_f64();
        let use_closed = mode == super::closed::ClosedMode::On && matches!(closed, Some(Ok(_)));
        let (mut taylor_seconds, mut taylor_bytes, mut project_seconds, mut project_bytes) =
            (0.0, 0, 0.0, 0);
        let today = if use_closed {
            None
        } else {
            if ir_hard_declare() {
                self.log_ir_hard_declaration(correction_index, correction, term_index, term, &hard);
            }
            let cost_started = Instant::now();
            let polynomial = hard
                .taylor(order)
                .wrap_err("while Taylor-expanding the prepared hard component")?;
            taylor_seconds = cost_started.elapsed().as_secs_f64();
            taylor_bytes = polynomial.as_view().get_byte_size();
            let project_started = Instant::now();
            let covariant = hard.project_hard(&polynomial)?;
            drop(polynomial);
            project_seconds = project_started.elapsed().as_secs_f64();
            project_bytes = covariant.as_view().get_byte_size();
            Some(hard_moments(&covariant, hard.soft_edges()).wrap_err_with(|| {
                format!(
                    "while separating hard moments: hard={}, remaining={}, soft={:?}, covariant={covariant}",
                    correction.hard.string_label(),
                    term.remainder.remaining.string_label(),
                    hard.soft_edges(),
                )
            })?)
        };
        if let Some(closed) = &closed {
            let (outcome, check) = match (closed, &today) {
                (Err(reason), _) => (format!("declined:{reason}"), String::new()),
                (Ok(_), None) => ("closed".to_owned(), String::new()),
                (Ok(moments), Some(reference)) => {
                    let point = |rows: &[(Atom, Atom)]| {
                        rows.iter()
                            .fold(Atom::Zero, |sum, (coefficient, insertion)| {
                                sum + super::closed::random_point(coefficient)
                                    * super::closed::random_point(insertion)
                            })
                    };
                    let mut ours = moments
                        .iter()
                        .map(|(_, i)| i.to_canonical_string())
                        .collect::<Vec<_>>();
                    let mut theirs = reference
                        .iter()
                        .map(|(_, i)| i.to_canonical_string())
                        .collect::<Vec<_>>();
                    ours.sort();
                    theirs.sort();
                    let equal = (point(moments) - point(reference)).is_zero();
                    (
                        "checked".to_owned(),
                        format!(
                            " equal={equal} insertions_same={} ours={} today={}",
                            ours == theirs,
                            ours.len(),
                            theirs.len()
                        ),
                    )
                }
            };
            eprintln!(
                "IR_CLOSED graph={} corr={correction_index} term={term_index} hard_loops={} rem_loops={} order={order} outcome={outcome}{check} slots={} sectors={} odd={} closures={} seconds={closed_seconds:.3} taylor_s={:.3} polarize_s={:.3} close_s={:.3} taylor_bytes={} closed_bytes={} certified={} uncertified={} dropped_zero={} zero_test_skipped={} today_taylor_s={taylor_seconds:.3} today_project_s={project_seconds:.3} today_taylor_bytes={taylor_bytes} today_project_bytes={project_bytes}",
                self.graph.name,
                self.graph.n_loops(&correction.hard),
                term.remainder.loops,
                closed_stats.slots,
                closed_stats.sectors,
                closed_stats.odd_sectors,
                closed_stats.closures,
                closed_stats.taylor_seconds,
                closed_stats.polarize_seconds,
                closed_stats.close_seconds,
                closed_stats.taylor_bytes,
                closed_stats.closed_bytes,
                closed_stats.certified,
                closed_stats.uncertified,
                closed_stats.dropped_zero,
                closed_stats.zero_test_skipped,
            );
        }
        let moments = match (today, closed) {
            (Some(today), _) => today,
            (None, Some(Ok(closed))) => closed,
            _ => unreachable!("the closed path runs only when it succeeded"),
        };
        Ok(TermMoments {
            moments,
            use_closed,
            taylor_seconds,
            taylor_bytes,
            project_seconds,
            project_bytes,
        })
    }

    /// The IR Dirac closure of an IR correction (cached per correction): the
    /// closure of its hard part's open Dirac line, or the decline reason. An
    /// error while reading the shape declines too (fail closed).
    fn ir_dirac_boundary(
        &mut self,
        correction_index: usize,
        correction: &IrCorrection,
    ) -> std::result::Result<DiracBoundary, String> {
        if let Some(boundary) = self.ir_dirac_boundaries.get(&correction_index) {
            return boundary.clone();
        }
        let boundary =
            super::ir_dirac::boundary(self.graph, &correction.hard, correction.taylor_order)
                .unwrap_or_else(|error| {
                    Err(format!("error:{error}")
                        .chars()
                        .map(|c| if c.is_whitespace() { '_' } else { c })
                        .take(160)
                        .collect())
                });
        self.ir_dirac_boundaries
            .insert(correction_index, boundary.clone());
        boundary
    }

    /// One `IR_DIRAC` line (`GAMMALOOP_IR_DIRAC_CLOSURE=on|check`).
    #[allow(clippy::too_many_arguments)]
    fn log_ir_dirac(
        &self,
        correction_index: usize,
        correction: &IrCorrection,
        term: &str,
        path: &str,
        rem_loops: usize,
        boundary: Option<&DiracBoundary>,
        outcome: &str,
        seconds: f64,
        today_seconds: Option<f64>,
    ) {
        #[cfg(test)]
        super::ir_dirac::count(format!(
            "{path}:{}:{}",
            super::ir_dirac::kind_label(boundary),
            outcome.split(':').take(2).collect::<Vec<_>>().join(":")
        ));
        eprintln!(
            "IR_DIRAC graph={} corr={correction_index} term={term} path={path} hard_loops={} rem_loops={rem_loops} omega={} kind={} outcome={outcome} seconds={seconds:.3}{}",
            self.graph.name,
            self.graph.n_loops(&correction.hard),
            correction.taylor_order,
            super::ir_dirac::kind_label(boundary),
            today_seconds.map_or_else(String::new, |seconds| format!(" today_s={seconds:.3}")),
        );
    }

    /// `GAMMALOOP_IR_DIRAC_CLOSURE=on|check` for one IR correction term on
    /// today's path: the moments of the hard part closed by its Dirac closure
    /// (`open` times the closure, in `open`'s loop basis) and reassembled
    /// ([`super::ir_dirac::reassemble`]), with that closed hard part for the
    /// downstream integration; on a decline, `open` and its moments. In
    /// `check` the open moments are computed too and must agree
    /// (`IR_DIRAC_CHECK`, an error otherwise). One `IR_DIRAC` line.
    fn ir_term_dirac(
        &mut self,
        correction_index: usize,
        correction: &IrCorrection,
        term_index: usize,
        term: &UvPreparationTerm,
        open: RStarIntegrand,
        order: usize,
    ) -> Result<(RStarIntegrand, TermMoments)> {
        let started = Instant::now();
        let label = term_index.to_string();
        let rem_loops = term.remainder.loops;
        let boundary = match self.ir_dirac_boundary(correction_index, correction) {
            Ok(boundary) => boundary,
            Err(reason) => {
                let today = self.ir_term_moments(
                    correction_index,
                    correction,
                    term_index,
                    term,
                    &open,
                    order,
                )?;
                self.log_ir_dirac(
                    correction_index,
                    correction,
                    &label,
                    "today",
                    rem_loops,
                    None,
                    &format!("declined:{reason}"),
                    started.elapsed().as_secs_f64(),
                    None,
                );
                return Ok((open, today));
            }
        };
        let hard = open.times_momentum_free(&boundary.closure)?;
        let closed =
            self.ir_term_moments(correction_index, correction, term_index, term, &hard, order)?;
        let rows = match super::ir_dirac::reassemble(&closed.moments, &boundary) {
            Ok(rows) => rows,
            Err(reason) => {
                drop(hard);
                let today = self.ir_term_moments(
                    correction_index,
                    correction,
                    term_index,
                    term,
                    &open,
                    order,
                )?;
                self.log_ir_dirac(
                    correction_index,
                    correction,
                    &label,
                    "today",
                    rem_loops,
                    Some(&boundary),
                    &format!("declined:{reason}"),
                    started.elapsed().as_secs_f64(),
                    None,
                );
                return Ok((open, today));
            }
        };
        let closed_seconds = started.elapsed().as_secs_f64();
        let today_seconds = if self.ir_dirac == IrDiracMode::Check {
            let open_started = Instant::now();
            let reference =
                self.ir_term_moments(correction_index, correction, term_index, term, &open, order)?;
            let open_seconds = open_started.elapsed().as_secs_f64();
            let agreement =
                super::ir_dirac::agree(&closed.moments, &rows, &reference.moments, &boundary)?;
            self.log_ir_dirac_check(
                correction_index,
                correction,
                &label,
                "today",
                &boundary,
                agreement,
                closed_seconds,
                open_seconds,
            )?;
            Some(open_seconds)
        } else {
            None
        };
        drop(open);
        self.log_ir_dirac(
            correction_index,
            correction,
            &label,
            "today",
            rem_loops,
            Some(&boundary),
            "closed",
            closed_seconds,
            today_seconds,
        );
        Ok((
            hard,
            TermMoments {
                moments: rows,
                ..closed
            },
        ))
    }

    /// One `IR_DIRAC_CHECK` line; an error when the closed moments differ
    /// from the open ones.
    #[allow(clippy::too_many_arguments)]
    fn log_ir_dirac_check(
        &self,
        correction_index: usize,
        correction: &IrCorrection,
        term: &str,
        path: &str,
        boundary: &DiracBoundary,
        agreement: super::ir_dirac::Agreement,
        closed_seconds: f64,
        open_seconds: f64,
    ) -> Result<()> {
        #[cfg(test)]
        super::ir_dirac::count_check(agreement);
        eprintln!(
            "IR_DIRAC_CHECK graph={} corr={correction_index} term={term} path={path} hard_loops={} omega={} kind={} identical={} equal={} structure={} coefficients={} closed_s={closed_seconds:.3} open_s={open_seconds:.3}",
            self.graph.name,
            self.graph.n_loops(&correction.hard),
            correction.taylor_order,
            super::ir_dirac::kind_label(Some(boundary)),
            agreement.identical,
            agreement.equal(),
            agreement.structure,
            agreement.coefficients,
        );
        ensure!(
            agreement.equal(),
            "IR_DIRAC_CHECK: the Dirac-closed IR hard part of {} corr={correction_index} term={term} ({path}) differs from the open one: {agreement:?}",
            self.graph.name
        );
        Ok(())
    }

    /// `GAMMALOOP_IR_HARD=child|child-check` for one IR correction: its
    /// contribution `Σ_F Σ_j T_comp[h_j Z̃_F(G/γ; E_j) closure]` from the
    /// integrated hard part ([`Self::ir_hard_child`]), or `None` when the
    /// child route declines (today's terms follow). In `child-check` today's
    /// terms are computed too, their pole parts must equal the child value's
    /// (an error otherwise), and today's value is used. One `IR_HARD_CHILD`
    /// line per correction.
    fn ir_correction_child(
        &mut self,
        correction_index: usize,
        correction: &IrCorrection,
        mode: super::closed::ChildMode,
    ) -> Result<Option<Atom>> {
        let started = Instant::now();
        // GAMMALOOP_IR_HARD_CHILD_GUARD=1 (opt-in): the closed path's scoped
        // budget around the whole H(γ). Off by default: its RSS checks run
        // inside every vakint integration of H and cost as much system time
        // as the work itself (L4s: 370 s sys against 4 s); the process
        // budget (BETA_PROCESS_RSS_LIMIT_GIB) applies either way.
        static GUARD: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
            matches!(
                std::env::var("GAMMALOOP_IR_HARD_CHILD_GUARD")
                    .as_deref()
                    .map(str::trim),
                Ok("1")
            )
        });
        // GAMMALOOP_IR_HARD_CHILD_SELECT=all|omega1|omega1-spinor (default
        // all): which corrections take the child route; the others stay on
        // the closed path. `omega1`: hard Taylor order ω >= 1 only (order-0
        // hard parts close to a few scalar channels, cheaper than an open
        // tensor integration); `omega1-spinor`: also ω = 0 with spinor slots
        // (which the closed path declines to the projected path).
        static SELECT: std::sync::LazyLock<u8> =
            std::sync::LazyLock::new(|| {
                match std::env::var("GAMMALOOP_IR_HARD_CHILD_SELECT")
                    .as_deref()
                    .map(str::trim)
                {
                    Ok("omega1") => 1,
                    Ok("omega1-spinor") => 2,
                    _ => 0,
                }
            });
        if *SELECT != 0 && correction.taylor_order < 1 {
            let spinor = *SELECT == 2
                && correction
                    .preparation
                    .iter()
                    .find(|term| term.counterterms.is_empty())
                    .and_then(|term| self.open_remainder(term).ok())
                    .and_then(|hard| super::closed::open_slot_list(hard.routed_expression()).ok())
                    .is_some_and(|slots| slots.iter().any(|(rep, _)| !rep.contains("mink")));
            if !spinor {
                eprintln!(
                    "IR_HARD_CHILD graph={} corr={correction_index} hard_loops={} omega={} outcome=skipped:select seconds=0.000",
                    self.graph.name,
                    self.graph.n_loops(&correction.hard),
                    correction.taylor_order,
                );
                return Ok(None);
            }
        }
        let attempt = if mode == super::closed::ChildMode::Check || !*GUARD {
            self.ir_hard_child(correction_index, correction)?
        } else {
            // Fail-closed under the closed path's scoped budget
            // (GAMMALOOP_IR_HARD_BUDGET_GIB over the RSS now): an exhausted
            // budget, an error or a panic declines the correction.
            let limit = super::closed::budget_extra_bytes()
                .and_then(|extra| super::expression::process_rss_now().map(|rss| rss + extra));
            super::closed::guarded(limit, super::expression::memory_release(), || {
                self.ir_hard_child(correction_index, correction)
            })
        };
        let seconds = started.elapsed().as_secs_f64();
        let (value, stats) = match attempt {
            Ok(outcome) => outcome,
            Err(reason) => {
                eprintln!(
                    "IR_HARD_CHILD graph={} corr={correction_index} hard_loops={} omega={} outcome=declined:{reason} seconds={seconds:.3}",
                    self.graph.name,
                    self.graph.n_loops(&correction.hard),
                    correction.taylor_order,
                );
                return Ok(None);
            }
        };
        let (value, check) = if mode == super::closed::ChildMode::Check {
            let today_started = Instant::now();
            let mut today = Atom::Zero;
            for (term_index, term) in correction.preparation.iter().enumerate() {
                if term
                    .remainder
                    .is_scaleless(self.graph, self.plan.marked_edge)
                {
                    continue;
                }
                today += self.ir_term_today(correction_index, correction, term_index, term)?;
            }
            let difference = self.poles(&(&today - &value))?;
            // A massive-route quotient Z̃ (GAMMALOOP_IR_QUOTIENT unset or
            // massive) carries non-local log(mUV) pieces that cancel only in
            // the complete sum; they depend on how the insertion is written
            // (its boundary momenta), so the two routes may differ there
            // alone. Classified, never accepted silently: `equal_mod_log_mUV`.
            let verdict = if difference.is_zero() {
                "equal"
            } else {
                let mass = GS.m_uv_vacuum;
                let without = difference
                    .replace_map(|view, _, out| {
                        if let symbolica::atom::AtomView::Fun(function) = view
                            && function.get_symbol() == Symbol::LOG
                            && view.contains_symbol(mass)
                        {
                            **out = Atom::Zero;
                        }
                    })
                    .expand();
                ensure!(
                    without.is_zero(),
                    "GAMMALOOP_IR_HARD=child-check: the child IR hard part differs for {} corr={correction_index}: poles {difference}",
                    self.graph.name
                );
                "equal_mod_log_mUV"
            };
            let today_zero = self.poles(&today)?.is_zero();
            (
                today,
                format!(
                    " check={verdict} today_poles_zero={today_zero} today_s={:.3}",
                    today_started.elapsed().as_secs_f64()
                ),
            )
        } else {
            (value, String::new())
        };
        eprintln!(
            "IR_HARD_CHILD graph={} corr={correction_index} hard_loops={} omega={} outcome=child{stats}{check} seconds={seconds:.3}",
            self.graph.name,
            self.graph.n_loops(&correction.hard),
            correction.taylor_order,
        );
        Ok(Some(value))
    }

    /// `GAMMALOOP_IR_HARD=child`: the integrated IR hard part of one
    /// correction, `H(γ) = Σ_σ ∫ T_ω[γ/σ · Π Δ(σ_i)]` with the overall `Δ(γ)`
    /// among the σ (`preparation(…, true)`, exactly today's terms), computed as
    /// its own support root: each remainder keeps its open boundary slots, is
    /// truncated to grade `ω` in the boundary momenta (the crown of γ: soft
    /// lines and root legs), contracted in the UV children's edge basis
    /// (streamed when large), Taylor-expanded at the fixed reference mass and
    /// integrated by the UV children's single-mass call
    /// ([`RStarIntegrand::integrate`] → `integrate_polynomial`: `to_vakint` +
    /// `integrate_single_mass`; no `K`, finite parts kept). The Laurent depth
    /// is the remainder's loops plus the quotient's loops, the pole bound of
    /// `Z̃_F(G/γ)`, so that every product with the quotient is exact to
    /// `ε^-1`. The σ terms are summed per remainder loop count `L` (the
    /// analytic outer Taylor acts on `(M²)^{a-εL}`), read off as `Σ_j h_j E_j`
    /// by [`hard_moments`] (the same `(coefficient, insertion)` split as
    /// today's), and passed to [`Self::ir_outer_integrated`].
    fn ir_hard_child(
        &mut self,
        correction_index: usize,
        correction: &IrCorrection,
    ) -> Result<std::result::Result<(Atom, String), String>> {
        if self.epsilon != GS.dim_epsilon {
            return Ok(Err("epsilon_symbol".into()));
        }
        let omega = correction.taylor_order;
        let Ok(degree) = usize::try_from(omega) else {
            return Ok(Err("negative_order".into()));
        };
        let hard_loops = self.graph.n_loops(&correction.hard);
        let quotient_loops = self
            .graph
            .n_loops(&self.plan.support)
            .saturating_sub(hard_loops);
        let edge_basis = super::expression::ir_hard_child_edge_basis(omega);
        let (closures_before, integrated_before) = super::integrate::single_mass_counts();
        // The Laurent depth of H beyond its own poles: the quotient's loops
        // (the pole bound of Z̃_F(G/γ)), raised when the quotient factors
        // actually carry more poles (an inserted pole part outside γ, as in
        // a propagator-insertion surrogate); then H is computed again.
        // GAMMALOOP_IR_DIRAC_CLOSURE=on|check: every remainder of H(γ) closed
        // by the hard part's Dirac closure, the moments reassembled.
        let dirac = if self.ir_dirac == IrDiracMode::Off {
            None
        } else {
            Some(self.ir_dirac_boundary(correction_index, correction))
        };
        let closure = dirac
            .as_ref()
            .and_then(|boundary| boundary.as_ref().ok())
            .cloned();
        let dirac_check = closure.is_some() && self.ir_dirac == IrDiracMode::Check;
        let mut extra = quotient_loops;
        let mut attempts = 0_usize;
        loop {
            attempts += 1;
            let started = Instant::now();
            let mut groups = BTreeMap::<usize, Atom>::new();
            // The open H(γ) per loop count and the times (IR_DIRAC_CHECK).
            let mut open_groups = BTreeMap::<usize, Atom>::new();
            let (mut closed_seconds, mut open_seconds) = (0.0, 0.0);
            let mut representative = None;
            let mut soft = Vec::new();
            let mut terms = 0_usize;
            for (term_index, term) in correction.preparation.iter().enumerate() {
                if term
                    .remainder
                    .is_scaleless(self.graph, self.plan.marked_edge)
                {
                    continue;
                }
                let term_started = Instant::now();
                let (closures, _) = super::integrate::single_mass_counts();
                let remainder = self.open_remainder(term)?;
                soft = remainder.soft_edges().to_vec();
                let (remainder, open) = match &closure {
                    Some(boundary) => {
                        let closed = remainder.times_momentum_free(&boundary.closure)?;
                        (closed, dirac_check.then_some(remainder))
                    }
                    None => (remainder, None),
                };
                let remainder = remainder.with_ir_hard_truncation(&soft, omega)?;
                let routed_bytes = remainder.routed_bytes();
                let value = self
                    .integrate_ir_hard_term(&remainder, &soft, degree, extra, edge_basis)
                    .wrap_err_with(|| {
                        format!(
                            "while integrating the IR hard part corr={correction_index} term={term_index} of {}",
                            self.graph.name
                        )
                    })?;
                drop(remainder);
                eprintln!(
                    "IR_HARD_CHILD_TERM graph={} corr={correction_index} term={term_index} attempt={attempts} hard_loops={hard_loops} rem_loops={} children={} omega={omega} edge={edge_basis} extra_depth={extra} routed_bytes={routed_bytes} h_bytes={} closures={} seconds={:.3}",
                    self.graph.name,
                    term.remainder.loops,
                    term.counterterms.len(),
                    value.as_view().get_byte_size(),
                    super::integrate::single_mass_counts().0 - closures,
                    term_started.elapsed().as_secs_f64(),
                );
                let term_seconds = term_started.elapsed().as_secs_f64();
                closed_seconds += term_seconds;
                let open_term_seconds = match open {
                    Some(open) => {
                        let open_started = Instant::now();
                        let open = open.with_ir_hard_truncation(&soft, omega)?;
                        let open_value = self
                            .integrate_ir_hard_term(&open, &soft, degree, extra, edge_basis)
                            .wrap_err_with(|| {
                                format!(
                                    "while integrating the open IR hard part corr={correction_index} term={term_index} of {} (IR_DIRAC_CHECK)",
                                    self.graph.name
                                )
                            })?;
                        *open_groups
                            .entry(term.remainder.loops)
                            .or_insert(Atom::Zero) += open_value;
                        let seconds = open_started.elapsed().as_secs_f64();
                        open_seconds += seconds;
                        Some(seconds)
                    }
                    None => None,
                };
                if let Some(boundary) = &dirac {
                    let outcome = match boundary {
                        Ok(_) => "closed".to_owned(),
                        Err(reason) => format!("declined:{reason}"),
                    };
                    self.log_ir_dirac(
                        correction_index,
                        correction,
                        &term_index.to_string(),
                        "child",
                        term.remainder.loops,
                        closure.as_ref(),
                        &outcome,
                        term_seconds,
                        open_term_seconds,
                    );
                }
                representative.get_or_insert(term_index);
                *groups.entry(term.remainder.loops).or_insert(Atom::Zero) += value;
                terms += 1;
            }
            let h_seconds = started.elapsed().as_secs_f64();
            let (closures_after, integrated_after) = super::integrate::single_mass_counts();
            let Some(representative) = representative else {
                return Ok(Ok((Atom::Zero, " terms=0".into())));
            };
            let mut rows = Vec::new();
            let mut h_bytes = 0;
            // Per loop count: the closed moments and their reassembly.
            let mut closed_moments =
                BTreeMap::<usize, (Vec<(Atom, Atom)>, Vec<(Atom, Atom)>)>::new();
            for (loops, h) in &groups {
                h_bytes += h.as_view().get_byte_size();
                if h.is_zero() {
                    continue;
                }
                let moments = match hard_moments(h, &soft) {
                    Ok(moments) => moments,
                    Err(error) => {
                        let reason = format!("{error}")
                            .chars()
                            .map(|c| if c.is_whitespace() { '_' } else { c })
                            .take(120)
                            .collect::<String>();
                        return Ok(Err(format!("moments:{reason}")));
                    }
                };
                let moments = match &closure {
                    Some(boundary) => match super::ir_dirac::reassemble(&moments, boundary) {
                        Ok(reassembled) => {
                            closed_moments.insert(*loops, (moments, reassembled.clone()));
                            reassembled
                        }
                        Err(reason) => {
                            // The closed H(γ) is not of the assumed form: the
                            // correction declines to today's terms.
                            self.log_ir_dirac(
                                correction_index,
                                correction,
                                "all",
                                "child",
                                *loops,
                                Some(boundary),
                                &format!("declined:{reason}"),
                                h_seconds,
                                None,
                            );
                            return Ok(Err(format!("dirac_{reason}")));
                        }
                    },
                    None => moments,
                };
                rows.extend(
                    moments
                        .into_iter()
                        .map(|(h, insertion)| (h, insertion, *loops)),
                );
            }
            if let (true, Some(boundary)) = (dirac_check, &closure) {
                let mut all_loops = groups
                    .keys()
                    .copied()
                    .collect::<std::collections::BTreeSet<_>>();
                all_loops.extend(open_groups.keys().copied());
                for loops in all_loops {
                    let (closed_rows, reassembled) =
                        closed_moments.remove(&loops).unwrap_or_default();
                    let open = open_groups.remove(&loops).unwrap_or_default();
                    let open_rows = if open.is_zero() {
                        Vec::new()
                    } else {
                        hard_moments(&open, &soft).wrap_err_with(|| {
                            format!(
                                "IR_DIRAC_CHECK: the open moments of {} corr={correction_index} loops={loops}",
                                self.graph.name
                            )
                        })?
                    };
                    let agreement =
                        super::ir_dirac::agree(&closed_rows, &reassembled, &open_rows, boundary)?;
                    self.log_ir_dirac_check(
                        correction_index,
                        correction,
                        &format!("all_rem{loops}"),
                        "child",
                        boundary,
                        agreement,
                        closed_seconds,
                        open_seconds,
                    )?;
                }
            }
            let hard = self.open_remainder(&correction.preparation[representative])?;
            let outer_started = Instant::now();
            let outcome = self.ir_outer_integrated(&hard, &correction.forest, &rows)?;
            let (value, outer_stats, valuation) = match outcome {
                Ok(outcome) => outcome,
                Err(reason) => return Ok(Err(reason)),
            };
            let needed = usize::try_from(-valuation.min(0)).unwrap_or(0);
            if needed > extra {
                eprintln!(
                    "IR_HARD_CHILD_DEPTH graph={} corr={correction_index} attempt={attempts} extra_depth={extra} needed={needed}",
                    self.graph.name
                );
                if attempts >= 3 {
                    return Ok(Err(format!("depth:{needed}")));
                }
                extra = needed;
                continue;
            }
            return Ok(Ok((
                value,
                format!(
                    " quotient_loops={quotient_loops} extra_depth={extra} attempts={attempts} terms={terms} groups={} rows={} h_bytes={h_bytes} h_s={h_seconds:.3} closures={} integrated_terms={} edge={edge_basis} forests={} outer_s={:.3}{outer_stats}",
                    groups.len(),
                    rows.len(),
                    closures_after - closures_before,
                    integrated_after - integrated_before,
                    correction.forest.len(),
                    outer_started.elapsed().as_secs_f64(),
                ),
            )));
        }
    }

    /// One remainder of an IR hard part `H(γ)` ([`Self::ir_hard_child`]):
    /// the UV children's outer integration (stream, else the whole
    /// contraction) in the IR hard scope, with boundary edges `soft`, degree
    /// `degree` and `quotient_loops` extra Laurent terms, no root projector.
    fn integrate_ir_hard_term(
        &mut self,
        remainder: &RStarIntegrand,
        soft: &[EdgeIndex],
        degree: usize,
        quotient_loops: usize,
        edge_basis: bool,
    ) -> Result<Atom> {
        self.integration_id += 1;
        let id = self.integration_id;
        let mut settings = self.settings.clone();
        settings.number_of_terms_in_epsilon_expansion =
            i64::try_from(remainder.region.loops + quotient_loops)?;
        let started = Instant::now();
        let (graph, reference, auxiliary, vakint, uv_settings) = (
            self.graph,
            self.reference_mass_squared,
            &self.auxiliary_mass_squared,
            self.vakint,
            self.uv_settings,
        );
        let result = super::expression::in_ir_hard_scope(edge_basis, || -> Result<Atom> {
            let integrate = |expression: &Atom| {
                remainder.integrate(
                    expression,
                    graph,
                    soft,
                    degree,
                    reference,
                    auxiliary,
                    vakint,
                    uv_settings,
                    &settings,
                    id,
                    None,
                )
            };
            if super::expression::uv_stream() != super::expression::UvStream::Off {
                let mut sum = crate::uv::approx::streamed::SumAccumulator::default();
                let streamed = remainder.stream_contracted(&mut |piece| {
                    sum.push(integrate(&piece)?);
                    Ok(())
                })?;
                if streamed {
                    return Ok(sum.finish());
                }
            }
            integrate(remainder.expression()?)
        });
        self.integration_elapsed += started.elapsed();
        result
    }

    /// The analytic outer step of [`Self::ir_outer_analytic`] for already
    /// integrated hard moments `rows = (h_j, E_j, L)`: `h_j ∝ (M²)^{a_j-εL}`
    /// (the auxiliary mass, `log(M²)` inside), `E_j` the insertion. For every
    /// forest and row: `z = Z̃_F(E_j) · closure` (the generic path's quotient
    /// call and closure), its outer orders `[λ^b] z`, times
    /// `Σ_k (-1)^k binom(a_j-εL, k)`; then the generic tail (summand
    /// simplification, root projection and regulation, tensor finish) with
    /// the `h_j` substituted last. `Err(reason)` declines the correction.
    fn ir_outer_integrated(
        &mut self,
        hard: &RStarIntegrand,
        forests: &[IrForestTerm],
        rows: &[(Atom, Atom, usize)],
    ) -> Result<std::result::Result<(Atom, String, i64), String>> {
        use symbolica::{function, id::Replacement};
        let projector = self.root_projector;
        let mode = self.top_mass.unwrap_or(TopMassTerms::Keep);
        let degree = usize::try_from(self.plan.taylor_order)?;
        let mass = top_mass::auxiliary_mass_squared();
        let mut powers = Vec::with_capacity(rows.len());
        for (h, _, _) in rows {
            let power = if h.is_zero() {
                None
            } else {
                match top_mass::mass_power(h, mass) {
                    Some(power) => Some(power),
                    None if h.expand().is_zero() => None,
                    None => return Ok(Err("mass_power".into())),
                }
            };
            powers.push(power);
        }
        let started = Instant::now();
        let placeholder = symbolica::symbol!("ir_hard_child_integral");
        let mut assembled = Vec::new();
        let mut uses = 0_usize;
        let mut valuation = 0_i64;
        for forest in forests {
            for (i, (_, insertion, loops)) in rows.iter().enumerate() {
                let Some(power) = powers[i] else {
                    continue;
                };
                let mut z = self.ir_forest(forest, insertion)?;
                if z.is_zero() {
                    continue;
                }
                if let Some(closure) = &self.root_closure {
                    z *= closure;
                }
                let orders = match hard.ir_outer_orders(
                    &z,
                    &self.root_edges,
                    degree,
                    self.reference_mass_squared,
                )? {
                    Ok(orders) => orders,
                    Err(reason) => return Ok(Err(reason)),
                };
                let h = function!(placeholder, Atom::num(i as i64));
                for (order, value) in orders.iter().enumerate() {
                    if value.is_zero() {
                        continue;
                    }
                    let factor = top_mass::compensation_factor(
                        mode,
                        power,
                        *loops,
                        order,
                        degree,
                        GS.dim_epsilon,
                    );
                    if factor.is_zero() {
                        continue;
                    }
                    let regulated = Integrated::dimensionally_regularized(value);
                    let series = regulated
                        .series(
                            GS.dim_epsilon,
                            0,
                            symbolica::poly::series::SeriesDepth::relative(1),
                        )
                        .wrap_err("could not find the Laurent order of an IR quotient factor")?;
                    if !series.is_zero() {
                        let exponent = series.get_trailing_exponent();
                        let Some(pole) = exponent
                            .is_integer()
                            .then(|| exponent.numerator().to_i64())
                            .flatten()
                        else {
                            return Ok(Err("quotient_epsilon_power".into()));
                        };
                        valuation = valuation.min(pole);
                    }
                    assembled.push(&h * factor * value);
                }
                uses += 1;
            }
        }
        let quotient_seconds = started.elapsed().as_secs_f64();
        self.integration_id += 1;
        let id = self.integration_id;
        let value = if assembled.is_empty() {
            Atom::Zero
        } else {
            let sum = Atom::add_many(assembled.iter());
            let simplified = RStarIntegrand::simplify_summands(sum.as_view())
                .wrap_err("while simplifying an IR hard child outer polynomial")?;
            let regulated = match projector {
                Some(projector) => Integrated::dimensionally_regularized(
                    &projector(&simplified)
                        .wrap_err("while projecting an IR hard child outer polynomial")?,
                ),
                None => simplified,
            };
            let finished = RStarIntegrand::finish_tensor(regulated, id)?;
            let replacements = rows
                .iter()
                .enumerate()
                .map(|(i, (h, _, _))| {
                    Replacement::new(
                        function!(placeholder, Atom::num(i as i64)).to_pattern(),
                        h.clone(),
                    )
                })
                .collect::<Vec<_>>();
            finished.replace_multiple(&replacements)
        };
        Ok(Ok((
            value,
            format!(
                " uses={uses} quotient_s={quotient_seconds:.3} powers={:?}",
                powers
                    .iter()
                    .flatten()
                    .collect::<std::collections::BTreeSet<_>>()
            ),
            valuation,
        )))
    }

    /// `GAMMALOOP_IR_OUTER=analytic` for one IR correction term of the
    /// top-level root: `Σ_F Σ_j T_comp[c_j Z̃_F(E_j) closure]` with each
    /// distinct hard coefficient `c_j` integrated once at fixed mass
    /// (`h (M²)^{a-εL}`), the quotient counterterms and the closure computed
    /// exactly as on the generic path, and the compensated outer Taylor
    /// applied analytically ([`top_mass::compensation_factor`]); then the
    /// generic path's summand simplification, root projection, regulation
    /// and tensor finish. `Err(reason)` declines the term (generic path).
    /// The hard integrals are requested deep enough for the poles of the
    /// product: the generic path's depth, raised by the pole order of the
    /// quotient factor.
    #[allow(clippy::type_complexity)]
    fn ir_outer_analytic(
        &mut self,
        hard: &RStarIntegrand,
        forests: &[IrForestTerm],
        moments: &[(Atom, Atom)],
    ) -> Result<std::result::Result<(Atom, String), String>> {
        use symbolica::{function, id::Replacement, poly::series::SeriesDepth};
        // As on the generic path: a projected root is projected and regulated
        // after the Taylor step; an unprojected one (a child root, a
        // propagator-insertion pole part) goes to the tensor finish as is.
        let projector = self.root_projector;
        if self.epsilon != GS.dim_epsilon {
            return Ok(Err("epsilon_symbol".into()));
        }
        let mode = self.top_mass.unwrap_or(TopMassTerms::Keep);
        let degree = usize::try_from(self.plan.taylor_order)?;
        let loops = hard.region.loops;
        let started = Instant::now();
        // The quotient factors' orders, per use of a distinct coefficient.
        let mut distinct = Vec::<Atom>::new();
        let mut index = std::collections::HashMap::<Atom, usize>::new();
        let mut uses = Vec::<(usize, Vec<Atom>)>::new();
        for forest in forests {
            for (coefficient, insertion) in moments {
                if coefficient.is_zero() {
                    continue;
                }
                let mut z = self.ir_forest(forest, insertion)?;
                if z.is_zero() {
                    continue;
                }
                if let Some(closure) = &self.root_closure {
                    z *= closure;
                }
                let orders = match hard.ir_outer_orders(
                    &z,
                    &self.root_edges,
                    degree,
                    self.reference_mass_squared,
                )? {
                    Ok(orders) => orders,
                    Err(reason) => return Ok(Err(reason)),
                };
                if orders.iter().all(Atom::is_zero) {
                    continue;
                }
                let i = *index.entry(coefficient.clone()).or_insert_with(|| {
                    distinct.push(coefficient.clone());
                    distinct.len() - 1
                });
                uses.push((i, orders));
            }
        }
        let quotient_seconds = started.elapsed().as_secs_f64();
        // Depth: the lowest ε power of each coefficient's quotient factors.
        let mut valuation = vec![0_i64; distinct.len()];
        for (i, orders) in &uses {
            for order in orders.iter().filter(|order| !order.is_zero()) {
                let regulated = Integrated::dimensionally_regularized(order);
                let series = regulated
                    .series(GS.dim_epsilon, 0, SeriesDepth::relative(1))
                    .wrap_err("could not find the Laurent order of an IR quotient factor")?;
                if series.is_zero() {
                    continue;
                }
                let exponent = series.get_trailing_exponent();
                let Some(power) = exponent
                    .is_integer()
                    .then(|| exponent.numerator().to_i64())
                    .flatten()
                else {
                    return Ok(Err("quotient_epsilon_power".into()));
                };
                valuation[*i] = valuation[*i].min(power);
            }
        }
        let hard_started = Instant::now();
        let mass = top_mass::auxiliary_mass_squared();
        let mut integrals = Vec::with_capacity(distinct.len());
        for (i, coefficient) in distinct.iter().enumerate() {
            self.integration_id += 1;
            let mut settings = self.settings.clone();
            settings.number_of_terms_in_epsilon_expansion =
                i64::try_from(loops)? - valuation[i].min(0);
            let integral = match hard.integrate_hard_coefficient(
                coefficient,
                self.graph,
                &self.root_edges,
                self.reference_mass_squared,
                &self.auxiliary_mass_squared,
                self.vakint,
                self.uv_settings,
                &settings,
                self.integration_id,
                projector.is_some(),
            )? {
                Ok(integral) => integral,
                Err(reason) => return Ok(Err(reason)),
            };
            let power = if integral.is_zero() {
                0
            } else {
                match top_mass::mass_power(&integral, mass) {
                    Some(power) => power,
                    None => return Ok(Err("mass_power".into())),
                }
            };
            integrals.push((integral, power));
        }
        let hard_seconds = hard_started.elapsed().as_secs_f64();
        self.integration_elapsed += hard_started.elapsed();
        let assemble_started = Instant::now();
        let placeholder = symbolica::symbol!("top_mass_ir_hard_integral");
        let mut assembled = Vec::new();
        for (i, orders) in &uses {
            let (integral, power) = &integrals[*i];
            if integral.is_zero() {
                continue;
            }
            let h = function!(placeholder, Atom::num(*i as i64));
            for (order, value) in orders.iter().enumerate() {
                if value.is_zero() {
                    continue;
                }
                let factor = top_mass::compensation_factor(
                    mode,
                    *power,
                    loops,
                    order,
                    degree,
                    GS.dim_epsilon,
                );
                if factor.is_zero() {
                    continue;
                }
                assembled.push(&h * factor * value);
            }
        }
        self.integration_id += 1;
        let id = self.integration_id;
        let value = if assembled.is_empty() {
            Atom::Zero
        } else {
            let sum = Atom::add_many(assembled.iter());
            let simplified = RStarIntegrand::simplify_summands(sum.as_view())
                .wrap_err("while simplifying an analytic IR outer polynomial")?;
            let regulated = match projector {
                Some(projector) => Integrated::dimensionally_regularized(
                    &projector(&simplified)
                        .wrap_err("while projecting an analytic IR outer polynomial")?,
                ),
                None => simplified,
            };
            let finished = RStarIntegrand::finish_tensor(regulated, id)?;
            let replacements = integrals
                .iter()
                .enumerate()
                .map(|(i, (integral, _))| {
                    Replacement::new(
                        function!(placeholder, Atom::num(i as i64)).to_pattern(),
                        integral.clone(),
                    )
                })
                .collect::<Vec<_>>();
            finished.replace_multiple(&replacements)
        };
        let stats = format!(
            " distinct={} uses={} quotient_s={quotient_seconds:.3} hard_s={hard_seconds:.3} assemble_s={:.3} powers={:?}",
            distinct.len(),
            uses.len(),
            assemble_started.elapsed().as_secs_f64(),
            integrals
                .iter()
                .map(|(_, power)| *power)
                .collect::<std::collections::BTreeSet<_>>(),
        );
        Ok(Ok((value, stats)))
    }

    /// `GAMMALOOP_IR_HARD_KEYS`: one `IR_HARD` line per IR preparation term of
    /// this evaluation, with the canonical identity of its hard component
    /// (marked line coloured; `key_nomark` without it) and of the term (the
    /// contracted spinney coloured as well). Plan-only: no Taylor operation.
    /// `GAMMALOOP_IR_HARD_CENSUS=1|only` (diagnostics, plan level): one
    /// `IR_HARD_CENSUS` line per IR correction of this evaluation's plan,
    /// from the insertion-free (σ = ∅) hard remainder: its open boundary
    /// slots by representation, the boundary vectors (crown edges; root legs
    /// among them; those present in the routed product), ω and the loops.
    fn log_ir_hard_census(&self) {
        let scope = if self.plan.root_boundary.is_some() {
            "child"
        } else {
            "root"
        };
        let external = self
            .graph
            .iter_edges_of(&self.graph.external_filter::<SuBitGraph>())
            .map(|(_, edge, _)| usize::from(edge))
            .collect::<std::collections::BTreeSet<_>>();
        for (correction_index, correction) in self.plan.ir.iter().enumerate() {
            let bare = correction
                .preparation
                .iter()
                .find(|term| term.counterterms.is_empty());
            let fields = match bare.map(|term| self.open_remainder(term)) {
                None => "slots=na reason=no_bare_term".to_owned(),
                Some(Err(error)) => format!(
                    "slots=na reason={}",
                    format!("{error}").replace(char::is_whitespace, "_")
                ),
                Some(Ok(hard)) => {
                    let routed = hard.routed_expression();
                    let soft = hard
                        .soft_edges()
                        .iter()
                        .map(|edge| usize::from(*edge))
                        .collect::<Vec<_>>();
                    let present = soft
                        .iter()
                        .filter(|edge| {
                            let mut found = false;
                            routed.visitor(&mut |part| {
                                if let symbolica::atom::AtomView::Fun(momentum) = part
                                    && momentum.get_symbol() == GS.emr_mom
                                    && usize::try_from(momentum.get(0)).ok() == Some(**edge)
                                {
                                    found = true;
                                }
                                !found
                            });
                            found
                        })
                        .count();
                    let roots = soft.iter().filter(|edge| external.contains(edge)).count();
                    let (slots, reps) = match super::closed::open_slot_list(routed) {
                        Ok(list) => {
                            let mut by_rep = BTreeMap::<String, usize>::new();
                            for (rep, _) in &list {
                                *by_rep.entry(rep.clone()).or_default() += 1;
                            }
                            (
                                list.len().to_string(),
                                by_rep
                                    .iter()
                                    .map(|(rep, n)| format!("{rep}:{n}"))
                                    .collect::<Vec<_>>()
                                    .join(","),
                            )
                        }
                        Err(reason) => ("na".to_owned(), reason.replace(char::is_whitespace, "_")),
                    };
                    let mink = super::closed::open_slot_list(routed)
                        .map(|list| {
                            list.iter()
                                .filter(|(rep, _)| rep.contains("mink"))
                                .count()
                                .to_string()
                        })
                        .unwrap_or_else(|_| "na".to_owned());
                    let spinor = super::closed::open_slot_list(routed)
                        .map(|list| {
                            list.iter()
                                .filter(|(rep, _)| rep.contains("bis"))
                                .count()
                                .to_string()
                        })
                        .unwrap_or_else(|_| "na".to_owned());
                    format!(
                        "slots={slots} mink_slots={mink} spinor_slots={spinor} slot_reps={reps} softs={} roots={roots} soft_present={present} soft_edges={soft:?}",
                        soft.len()
                    )
                }
            };
            eprintln!(
                "IR_HARD_CENSUS graph={} scope={scope} support_loops={} marked={} corr={correction_index} hard_loops={} omega={} terms={} live_terms={} forests={} {fields}",
                self.graph.name,
                self.graph.n_loops(&self.plan.support),
                usize::from(self.plan.marked_edge),
                self.graph.n_loops(&correction.hard),
                correction.taylor_order,
                correction.preparation.len(),
                correction
                    .preparation
                    .iter()
                    .filter(|term| !term
                        .remainder
                        .is_scaleless(self.graph, self.plan.marked_edge))
                    .count(),
                correction.forest.len(),
            );
        }
    }

    /// `GAMMALOOP_IR_HARD_DECLARE=1` (diagnostics): the boundary declaration
    /// of one IR hard term on today's path, printed just before its hard
    /// Taylor step, with the number of the next `GAMMALOOP_ALGEBRA_DUMP` file
    /// (its `hard_routed-<n>` dump when that stage is dumped).
    fn log_ir_hard_declaration(
        &self,
        correction_index: usize,
        correction: &IrCorrection,
        term_index: usize,
        term: &UvPreparationTerm,
        hard: &RStarIntegrand,
    ) {
        let routed = hard.routed_expression();
        let external = self
            .graph
            .iter_edges_of(&self.graph.external_filter::<SuBitGraph>())
            .map(|(_, edge, _)| usize::from(edge))
            .collect::<std::collections::BTreeSet<_>>();
        let slots = match super::closed::open_slot_list(routed) {
            Ok(list) => list
                .iter()
                .map(|(rep, slot)| format!("{rep}|{slot}"))
                .collect::<Vec<_>>()
                .join(";"),
            Err(reason) => format!("error:{reason}"),
        };
        let softs = hard
            .soft_edges()
            .iter()
            .map(|edge| {
                let edge = usize::from(*edge);
                format!(
                    "{edge}:{}",
                    if external.contains(&edge) {
                        "root"
                    } else {
                        "soft"
                    }
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let degree = usize::try_from(correction.taylor_order).unwrap_or(0);
        let sectors = match super::closed::taylor_sectors(routed, hard.soft_edges(), degree) {
            Ok(sectors) => sectors
                .iter()
                .map(|(edges, nonzero)| {
                    format!(
                        "[{}]{}",
                        edges
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(","),
                        if *nonzero { "" } else { "=0" }
                    )
                })
                .collect::<Vec<_>>()
                .join(";"),
            Err(error) => format!(
                "error:{}",
                format!("{error}").replace(char::is_whitespace, "_")
            ),
        };
        let (loop_edges, routes) = hard.route_table();
        eprintln!(
            "IR_HARD_DECL graph={} next_dump={} corr={correction_index} term={term_index} hard_loops={} rem_loops={} children={:?} omega={} marked={} reference_mass={} slots={slots} softs={softs} sectors={sectors} loop_edges={loop_edges:?} routes={} routed_bytes={}",
            self.graph.name,
            crate::uv::approx::streamed::dump_next(),
            self.graph.n_loops(&correction.hard),
            term.remainder.loops,
            term.counterterms
                .iter()
                .map(|call| call.argument.support.string_label())
                .collect::<Vec<_>>(),
            correction.taylor_order,
            usize::from(self.plan.marked_edge),
            self.reference_mass_squared,
            routes
                .iter()
                .map(|(edge, momentum)| format!("{edge}:{momentum}"))
                .collect::<Vec<_>>()
                .join(";"),
            routed.as_view().get_byte_size(),
        );
    }

    fn log_ir_hard_keys(&self) {
        let scope = if self.plan.root_boundary.is_some() {
            "child"
        } else {
            "root"
        };
        let support_loops = self.graph.n_loops(&self.plan.support);
        let empty = std::collections::BTreeSet::new();
        let digest = |support: &SuBitGraph,
                      marked: Option<EdgeIndex>,
                      decorated: &std::collections::BTreeSet<usize>| {
            crate::uv::counterterm_identity::CountertermIdentity::for_rstar_hard(
                self.graph, support, self.model, marked, decorated,
            )
            .map(|identity| format!("{:016x}", identity.key_digest()))
            .unwrap_or_else(|error| format!("err:{}", format!("{error}").replace(' ', "_")))
        };
        let support_key = digest(&self.plan.support, Some(self.plan.marked_edge), &empty);
        eprintln!(
            "IR_HARD_PLAN graph={} scope={scope} support_loops={support_loops} support_key={support_key} marked={} ordinary={} ordinary_live={} corrections={} taylor_order={}",
            self.graph.name,
            usize::from(self.plan.marked_edge),
            self.plan.ordinary_uv.len(),
            self.plan
                .ordinary_uv
                .iter()
                .filter(|term| !term
                    .remainder
                    .is_scaleless(self.graph, self.plan.marked_edge))
                .count(),
            self.plan.ir.len(),
            self.plan.taylor_order,
        );
        for (correction_index, correction) in self.plan.ir.iter().enumerate() {
            let boundary = self.graph.dummy_less_full_crown(&correction.hard);
            let mut legs = Vec::new();
            let mut external_legs = 0;
            for (pair, edge, data) in self.graph.iter_edges_of(&boundary) {
                if !pair.is_paired() {
                    external_legs += 1;
                }
                legs.push(format!(
                    "{}:{}",
                    usize::from(edge),
                    data.data
                        .particle()
                        .map(|particle| particle.index() as i64)
                        .unwrap_or(-1)
                ));
            }
            let key = digest(&correction.hard, Some(self.plan.marked_edge), &empty);
            let key_nomark = digest(&correction.hard, None, &empty);
            let live = correction
                .preparation
                .iter()
                .filter(|term| {
                    !term
                        .remainder
                        .is_scaleless(self.graph, self.plan.marked_edge)
                })
                .count();
            for (term_index, term) in correction.preparation.iter().enumerate() {
                let scaleless = term
                    .remainder
                    .is_scaleless(self.graph, self.plan.marked_edge);
                let contracted = term
                    .remainder
                    .insertions
                    .iter()
                    .flat_map(|slot| {
                        self.graph
                            .iter_edges_of(&slot.support)
                            .filter(|(pair, _, _)| pair.is_paired())
                            .map(|(_, edge, _)| usize::from(edge))
                            .collect::<Vec<_>>()
                    })
                    .collect::<std::collections::BTreeSet<_>>();
                let term_key = if scaleless {
                    "-".to_owned()
                } else {
                    digest(&correction.hard, Some(self.plan.marked_edge), &contracted)
                };
                eprintln!(
                    "IR_HARD graph={} scope={scope} support_loops={support_loops} marked={} corr={correction_index} hard={} hard_loops={} order={} legs={} ext_legs={external_legs} leg_particles={} key={key} key_nomark={key_nomark} terms={} live={live} forests={} forest_steps={} term={term_index} rem_loops={} children={} scaleless={scaleless} term_key={term_key}",
                    self.graph.name,
                    usize::from(self.plan.marked_edge),
                    correction.hard.string_label(),
                    self.graph.n_loops(&correction.hard),
                    correction.taylor_order,
                    legs.len(),
                    legs.join(","),
                    correction.preparation.len(),
                    correction.forest.len(),
                    correction
                        .forest
                        .iter()
                        .map(|forest| forest.steps.len())
                        .sum::<usize>(),
                    term.remainder.loops,
                    term.counterterms.len(),
                );
            }
        }
    }

    /// The remainder-boundary summary of this evaluation, when the mode is
    /// not off (`BOUNDARY_SUMMARY`, DESIGN 8).
    fn log_boundary_summary(&self, outcome: &str) {
        if !matches!(
            boundary::remainder_boundary_mode(),
            Ok(BoundaryMode::On | BoundaryMode::Check)
        ) {
            return;
        }
        let root = if self.root_closure.is_some() && self.root_projector.is_some() {
            if self.plan.root_boundary.is_some() {
                "child"
            } else {
                "closed"
            }
        } else if self.plan.root_boundary.is_some() || self.root_projector.is_none() {
            "child"
        } else {
            "open"
        };
        self.summary.log(&self.graph.name, root, outcome);
    }

    fn log_timing(&self) {
        self.log_boundary_summary("ok");
        let elapsed = self.started.elapsed();
        let preparation = elapsed.saturating_sub(
            self.proper_uv_elapsed + self.quotient_elapsed + self.integration_elapsed,
        );
        if super::routing::log_level() >= 1 {
            eprintln!(
                "RR_PHASE graph={} child={} total={:.3} proper_uv={:.3} ir_quotient={:.3} outer_integration={:.3} preparation={:.3} outer_integrations={}",
                self.graph.name,
                self.plan.root_boundary.is_some(),
                elapsed.as_secs_f64(),
                self.proper_uv_elapsed.as_secs_f64(),
                self.quotient_elapsed.as_secs_f64(),
                self.integration_elapsed.as_secs_f64(),
                preparation.as_secs_f64(),
                self.integration_id,
            );
        }
        tracing::info!(
            uv = true,
            rstar = true,
            timing = true,
            stage = "rstar_counterterm_phase_timing",
            graph = %self.graph.name,
            root_projected = self.root_projector.is_some(),
            total_ms = elapsed.as_millis() as u64,
            proper_uv_ms = self.proper_uv_elapsed.as_millis() as u64,
            ir_quotient_ms = self.quotient_elapsed.as_millis() as u64,
            outer_integration_ms = self.integration_elapsed.as_millis() as u64,
            preparation_ms = preparation.as_millis() as u64,
            quotient_calls = self.quotient_counterterms.len(),
            outer_integrations = self.integration_id,
            "R* counterterm phase timings"
        );
    }

    fn poles(&self, expression: &Atom) -> Result<Atom> {
        Ok(expression.series(self.epsilon, 0, -1)?.to_atom().expand())
    }

    fn counterterm_from_pole_part(&mut self, bare_poles: &Atom) -> Result<Atom> {
        let (_, external_routing) = RStarIntegrand::physical_external_routing(
            self.graph,
            &self.graph.external_filter::<SuBitGraph>(),
        )?;
        // Closed bare poles stand where the metric projector contracted them.
        let bare_poles = match &self.root_closure {
            Some(closure) => bare_poles * closure,
            None => bare_poles.clone(),
        };
        // GAMMALOOP_RSTAR_TERMS: here the bare term is the supplied bare poles.
        let bare_poles = if self.keeps_term(false) {
            bare_poles
        } else {
            Atom::Zero
        };
        let physical = if self.rstar_terms.keeps(true) {
            self.physical_uv_preparation()?
        } else {
            self.skipped_terms += 1;
            Atom::Zero
        };
        let mut prepared = bare_poles + physical;
        if let Some(projector) = self.root_projector {
            prepared = projector(&prepared)?;
        }
        // The independent bare-pole provider and the root projector can both
        // introduce tensor contractions with evanescent scalar coefficients.
        // Complete them before K, just as for an ordinary integrated UV node.
        let prepared = RStarIntegrand::finish_tensor(
            super::linearize::linearize(prepared.replace_multiple(&external_routing)),
            self.integration_id + 1,
        )?;
        let result = -self.poles(&prepared)?;
        self.log_timing();
        Ok(result)
    }

    /// The physical, nonempty proper-UV sum, before the root pole projection.
    /// Child counterterms remain in their original scopes and multiply their
    /// quotients before the scalar integral's Laurent depth is determined.
    fn physical_uv_preparation(&mut self) -> Result<Atom> {
        let mut result = Atom::Zero;
        for term in self
            .plan
            .ordinary_uv
            .iter()
            .filter(|term| !term.counterterms.is_empty())
        {
            let remainder = self.physical_remainder(term)?;
            self.integration_id += 1;
            let mut settings = self.settings.clone();
            settings.number_of_terms_in_epsilon_expansion = i64::try_from(term.remainder.loops)?;
            let started = Instant::now();
            result += remainder.integrate_physical(
                self.vakint,
                self.uv_settings,
                &settings,
                self.integration_id,
            )?;
            self.integration_elapsed += started.elapsed();
        }
        Ok(result)
    }

    fn physical_remainder(&self, term: &UvPreparationTerm) -> Result<RStarIntegrand> {
        let insertions = self
            .insertions(term)?
            .into_iter()
            .map(|insertion| insertion.replace(self.reference_mass_squared).with(0))
            .collect::<Vec<_>>();
        RStarIntegrand::new_with_closure(
            self.graph,
            self.model,
            &term.remainder,
            &insertions,
            self.plan.marked_edge,
            Atom::Zero,
            self.root_closure.as_ref(),
        )
    }

    fn insertions(&self, term: &UvPreparationTerm) -> Result<Vec<Atom>> {
        term.counterterms
            .iter()
            .map(|child| {
                self.children
                    .get(&child.argument.support)
                    .cloned()
                    .ok_or_else(|| eyre!("missing a completed proper UV counterterm"))
            })
            .collect()
    }

    /// A root-scope remainder, closed when this evaluation closes the root.
    fn remainder(&self, term: &UvPreparationTerm) -> Result<RStarIntegrand> {
        self.remainder_with(term, self.root_closure.as_ref())
            .and_then(|remainder| {
                remainder.with_root_truncation(
                    self.plan.root_boundary.is_some(),
                    &self.root_edges,
                    self.plan.taylor_order,
                )
            })
    }

    /// A remainder that keeps its open root slots (an IR hard part).
    fn open_remainder(&self, term: &UvPreparationTerm) -> Result<RStarIntegrand> {
        self.remainder_with(term, None)
    }

    fn remainder_with(
        &self,
        term: &UvPreparationTerm,
        closure: Option<&Atom>,
    ) -> Result<RStarIntegrand> {
        let insertions = self.insertions(term)?;
        RStarIntegrand::new_with_closure(
            self.graph,
            self.model,
            &term.remainder,
            &insertions,
            self.plan.marked_edge,
            Atom::var(self.reference_mass_squared),
            closure,
        )
    }

    fn quotient(&mut self, call: &CountertermCall, insertion: &Atom) -> Result<Atom> {
        let key = (call.clone(), insertion.clone());
        if let Some(value) = self.quotient_counterterms.get(&key) {
            return Ok(value.clone());
        }
        let started = Instant::now();
        // GAMMALOOP_IR_QUOTIENT: the massive R operation (default) or R* on
        // the same decorated graph.
        let mode = top_mass::ir_quotient()?;
        let massive = |this: &Self| -> Result<Atom> {
            this.uv
                .counterterm(call, std::slice::from_ref(insertion))
                .wrap_err_with(|| {
                    format!(
                        "while applying massive R to the decorated quotient {}",
                        call.argument.remaining.string_label()
                    )
                })
        };
        let use_rstar = match mode {
            top_mass::IrQuotient::Rstar => true,
            top_mass::IrQuotient::Auto => call.argument.loops > 4,
            _ => false,
        };
        let (result, route) = match mode {
            top_mass::IrQuotient::Check => {
                let reference = massive(self)?;
                let massive_seconds = started.elapsed().as_secs_f64();
                let rstar_started = Instant::now();
                let outcome = match self.quotient_rstar(call, insertion)? {
                    Ok(value) => {
                        let difference = (&reference - &value).expand();
                        let without_muv =
                            difference.replace(GS.m_uv_vacuum).with(Atom::Zero).expand();
                        let sample = if without_muv.is_zero() {
                            String::new()
                        } else {
                            let text = format!(
                                " degree={:?} insertion={} difference={}",
                                insertion_degree(insertion),
                                insertion,
                                without_muv
                            );
                            text.chars().take(1500).collect()
                        };
                        format!(
                            "equal={} equal_without_muv={} massive_has_muv={} massive_s={massive_seconds:.3} rstar_s={:.3}{sample}",
                            difference.is_zero(),
                            without_muv.is_zero(),
                            reference.contains_symbol(GS.m_uv_vacuum),
                            rstar_started.elapsed().as_secs_f64()
                        )
                    }
                    Err(reason) => format!("declined:{reason}"),
                };
                eprintln!(
                    "IR_QUOTIENT_CHECK graph={} loops={} support={} {outcome}",
                    self.graph.name,
                    call.argument.loops,
                    call.argument.remaining.string_label()
                );
                (reference, "massive")
            }
            _ if use_rstar => match self.quotient_rstar(call, insertion)? {
                Ok(value) => (value, "rstar"),
                Err(reason) => {
                    eprintln!(
                        "IR_QUOTIENT graph={} loops={} declined={reason}",
                        self.graph.name, call.argument.loops
                    );
                    (massive(self)?, "massive_declined")
                }
            },
            _ => (massive(self)?, "massive"),
        };
        if mode != top_mass::IrQuotient::Unset && mode != top_mass::IrQuotient::Check {
            eprintln!(
                "IR_QUOTIENT graph={} loops={} route={route} seconds={:.3}",
                self.graph.name,
                call.argument.loops,
                started.elapsed().as_secs_f64()
            );
        }
        self.quotient_elapsed += started.elapsed();
        self.quotient_counterterms.insert(key, result.clone());
        Ok(result)
    }

    /// `GAMMALOOP_IR_QUOTIENT=rstar|auto|check`: the decorated quotient
    /// counterterm by R* on the decorated graph the massive route builds
    /// (the contracted hard part is one vertex carrying the insertion
    /// polynomial), through [`RStarPlan::rstar_child_counterterm`] (automatic
    /// marked line, recursive R* children, Δ-cache under the R* route's key).
    /// `Err(reason)` when the quotient is not an admissible R* support root.
    fn quotient_rstar(
        &self,
        call: &CountertermCall,
        insertion: &Atom,
    ) -> Result<std::result::Result<Atom, String>> {
        if call.argument.loops == 0 || insertion.is_zero() {
            return Ok(Ok(Atom::Zero));
        }
        // The decorated vertex carries the insertion's own momentum degree:
        // R* power-counts UV and IR subgraphs with it (the massive forest
        // recomputes it from the polynomial; the structural bound would
        // misstate the IR degrees).
        let Some(degree) = insertion_degree(insertion) else {
            return Ok(Err("inhomogeneous_insertion".into()));
        };
        let graph = self.uv.decorated_graph_with_dod(
            &call.argument,
            std::slice::from_ref(insertion),
            Some(&[degree]),
        )?;
        let support = &call.argument.remaining;
        if graph.compute_dod(support) < 0 {
            return Ok(Ok(Atom::Zero));
        }
        if !RStarPlan::admits_support_root(&graph, support) {
            return Ok(Err("support_structure".into()));
        }
        let cache =
            self.plan.counterterm_cache.clone().unwrap_or_else(|| {
                Arc::new(crate::uv::counterterm_cache::CountertermCache::default())
            });
        let route = match self.plan.uv_children {
            UvChildrenRoute::Massive => UvChildrenRoute::Rstar { min_loops: 1 },
            route => route,
        };
        let delta = RStarPlan::rstar_child_counterterm(
            &cache,
            &graph,
            support,
            self.model,
            self.vakint,
            self.uv_settings,
            route,
            self.plan.uv_child_closure,
        )
        .wrap_err_with(|| {
            format!(
                "while applying R* to the decorated quotient {}",
                call.argument.remaining.string_label()
            )
        })?;
        Ok(Ok(delta))
    }

    fn ir_forest(&mut self, forest: &IrForestTerm, insertion: &Atom) -> Result<Atom> {
        let mut rows = vec![(Atom::num(forest.sign as i64), insertion.clone())];
        for (index, call) in forest.steps.iter().enumerate() {
            let mut next = Vec::new();
            for (coefficient, polynomial) in rows {
                let counterterm = self.quotient(call, &polynomial)?;
                if index + 1 == forest.steps.len() {
                    next.push((coefficient, counterterm));
                    continue;
                }
                // The c-operation product is coefficientwise in epsilon.
                // Retain the complete tensor polynomial inside the next UV
                // call: its metric contractions can depend on the dimension.
                for (power, tensor_polynomial) in counterterm.series(self.epsilon, 0, -1)?.terms() {
                    if !tensor_polynomial.is_zero() {
                        next.push((
                            &coefficient * Atom::var(self.epsilon).pow(power),
                            tensor_polynomial.clone(),
                        ));
                    }
                }
            }
            rows = next;
        }
        Ok(rows
            .into_iter()
            .fold(Atom::Zero, |sum, (coefficient, polynomial)| {
                sum + coefficient * polynomial
            }))
    }

    /// Remainder-boundary admission (DESIGN 10a): every reason is silent and
    /// counted, never an error, except the vakint gates under STRICT.
    fn boundary_admits(&self, remainder: &RStarIntegrand) -> std::result::Result<(), String> {
        if self.closed_root.is_none() {
            return Err("no_closed_root".into());
        }
        if self.root_closure.is_none() || self.root_projector.is_none() {
            return Err("open_root".into());
        }
        if self.plan.root_boundary.is_some() {
            return Err("proper_support_root".into());
        }
        if remainder.region.loops < 2 {
            return Err("one_loop_tree_cut".into());
        }
        if !RStarIntegrand::plain_remainder_contraction() {
            return Err("edge_contraction".into());
        }
        if [
            "direct_lowering",
            "tensor_factors",
            "topology_groups",
            "term_rewrites",
            "factor_maps",
            "dot_notation",
        ]
        .iter()
        .any(|name| crate::uv::approx::streamed::legacy(name))
        {
            return Err("algebra_legacy".into());
        }
        if self.settings.verify_numerator_identification {
            return Err("strict_split_branch".into());
        }
        boundary::vakint_gates()
            .and_then(|()| boundary::settings_gate(&self.settings))
            .map_err(|reason| format!("vakint_gates:{reason}"))
    }

    /// One root-scope ordinary remainder (the ordinary loop and the census
    /// bare-term loop): today's [`Self::integrate`] with the mode off; the
    /// remainder-boundary path for an admitted remainder otherwise (DESIGN
    /// 6a, 7).
    fn integrate_ordinary(&mut self, term: &UvPreparationTerm) -> Result<Atom> {
        let mode = boundary::remainder_boundary_mode()?;
        let remainder = self.remainder(term)?;
        if mode == BoundaryMode::Off {
            if let Some(result) = self.integrate_streamed(&remainder)? {
                return Ok(result);
            }
            return self.integrate(&remainder, remainder.expression()?);
        }
        self.summary.remainders += 1;
        let strict = boundary::boundary_strict();
        if let Err(reason) = self.boundary_admits(&remainder) {
            if strict && reason.starts_with("vakint_gates") {
                return Err(eyre!(
                    "BOUNDARY strict: the vakint gates reject an admissible remainder of {} ({reason})",
                    self.graph.name
                ));
            }
            *self.summary.not_admitted.entry(reason).or_default() += 1;
            if let Some(result) = self.integrate_streamed(&remainder)? {
                return Ok(result);
            }
            return self.integrate(&remainder, remainder.expression()?);
        }
        let errors = self.summary.errors;
        let result = self.integrate_admitted(&remainder, mode, strict);
        if result.is_err() && self.summary.errors == errors {
            // A hard error of an admitted remainder (the build, the
            // evaluation, or today's path in check mode or as a fallback).
            self.summary.errors += 1;
        }
        result
    }

    /// One admitted remainder (DESIGN 6a, 7): the build, then a decline
    /// (fallback, or an error under STRICT or when large), the evaluation
    /// (mode on) or the comparison with today's path (check).
    fn integrate_admitted(
        &mut self,
        remainder: &RStarIntegrand,
        mode: BoundaryMode,
        strict: bool,
    ) -> Result<Atom> {
        let closed_root = self
            .closed_root
            .expect("admission requires the closed-root normalisation");
        let id = self.integration_id + 1;
        let started = Instant::now();
        let mut settings = self.settings.clone();
        settings.number_of_terms_in_epsilon_expansion = i64::try_from(remainder.region.loops)?;
        let degree = usize::try_from(self.plan.taylor_order)?;
        let built = boundary::build_boundary(
            remainder,
            &BoundaryRequest {
                graph_name: &self.graph.name,
                marked_edge: self.plan.marked_edge,
                integration_id: id,
                root_edges: &self.root_edges,
                degree,
                reference_mass_squared: self.reference_mass_squared,
                auxiliary_mass_squared: &self.auxiliary_mass_squared,
                vakint: self.vakint,
                settings: &settings,
                additional_normalization: &self.uv_settings.vakint.additional_normalization,
                epsilon: self.epsilon,
                closed_root,
                mass_terms: self.top_mass.unwrap_or(TopMassTerms::Keep),
            },
        )
        .wrap_err_with(|| {
            format!(
                "while building the remainder boundary {id} of {}",
                self.graph.name
            )
        })?;
        let built = match built {
            Built::Declined(decline, stats) => {
                let large = boundary::large_reason(
                    &stats,
                    remainder.region.loops,
                    boundary::census_dump_active(),
                );
                eprintln!(
                    "BOUNDARY_FALLBACK graph={} id={id} loops={} stage={} reason={} contracted={} p_summands={} p_bytes={} large={}",
                    self.graph.name,
                    remainder.region.loops,
                    decline.stage,
                    decline.reason,
                    stats.contracted,
                    stats.p_summands,
                    stats.p_bytes,
                    large.is_some()
                );
                self.summary.declined(&decline);
                self.summary.contracted(&stats);
                self.boundary_elapsed += started.elapsed();
                self.summary.seconds_boundary += started.elapsed().as_secs_f64();
                if let Some(why) = large {
                    self.summary.errors += 1;
                    return Err(eyre!(
                        "BOUNDARY_DECLINED_LARGE graph={} id={id} loops={}: {}:{} ({why}, {} the contraction; no fallback to today's path)",
                        self.graph.name,
                        remainder.region.loops,
                        decline.stage,
                        decline.reason,
                        if stats.contracted { "after" } else { "before" }
                    ));
                }
                if strict {
                    self.summary.errors += 1;
                    return Err(eyre!(
                        "BOUNDARY strict: the remainder boundary {id} of {} declined at {}: {}",
                        self.graph.name,
                        decline.stage,
                        decline.reason
                    ));
                }
                return self.integrate(remainder, remainder.expression()?);
            }
            built => built,
        };
        match mode {
            BoundaryMode::Check => self.check_boundary(remainder, built, &settings, id, started),
            _ => self.evaluate_boundary(remainder, built, &settings, id, started),
        }
    }

    /// Mode on: the value of a built boundary, finished as today's outer
    /// integration (one integration id, kind `single_mass`).
    fn evaluate_boundary(
        &mut self,
        remainder: &RStarIntegrand,
        built: Built,
        settings: &VakintSettings,
        id: usize,
        started: Instant,
    ) -> Result<Atom> {
        let loops = remainder.region.loops;
        let (value, mut stats, built) = match built {
            Built::Zero(stats) => {
                self.summary.zero += 1;
                (Atom::Zero, stats, None)
            }
            Built::Boundary(built, stats) => (Atom::Zero, stats, Some(built)),
            Built::Declined(..) => unreachable!("declines are handled by the caller"),
        };
        self.summary.contracted(&stats);
        let mut value = value;
        let evaluate_started = Instant::now();
        let mut dump_seconds = 0.0;
        if let Some(built) = &built {
            self.summary.boundary += 1;
            let p_loops = built.loops.saturating_sub(1);
            let skipped =
                super::pint6_dump::settings().is_some_and(|dump| p_loops >= dump.skip_min_p_loops);
            if let Some(dump) = super::pint6_dump::settings() {
                // The census record, before the evaluation (DESIGN 6b).
                let dump_started = Instant::now();
                self.dump_boundary(dump, built, &stats, settings, !skipped, None)
                    .wrap_err_with(|| {
                        format!(
                            "while dumping the remainder boundary {id} of {}",
                            self.graph.name
                        )
                    })?;
                dump_seconds = dump_started.elapsed().as_secs_f64();
            }
            if !skipped {
                // A large remainder runs vakint strict: a pattern decline is
                // an error, never an `ordinary()` expansion (DESIGN 10c).
                let options = vakint::p_integral::BoundaryOptions {
                    strict: boundary::boundary_strict() || boundary::is_large(&stats, loops),
                };
                let (evaluated, reads, hits) =
                    boundary::evaluate_boundary(&built.boundary, self.vakint, settings, options)
                        .wrap_err_with(|| {
                            format!(
                                "while evaluating the remainder boundary {id} of {}",
                                self.graph.name
                            )
                        })?;
                value = evaluated;
                stats.block_reads = reads;
                stats.block_hits = hits;
            }
        }
        let evaluate_seconds = evaluate_started.elapsed().as_secs_f64() - dump_seconds;
        if super::pint6_dump::settings().is_some() {
            stats.seconds.push(("dump", dump_seconds));
        }
        stats.seconds.push(("evaluate", evaluate_seconds));
        stats.rss_gb[3] = boundary::rss_gb();
        let finish_started = Instant::now();
        let result = RStarIntegrand::finish_tensor(value, id)?;
        let finish_seconds = finish_started.elapsed().as_secs_f64();
        stats.seconds.push(("finish", finish_seconds));
        self.integration_id = id;
        let stage = |name: &str| {
            stats
                .seconds
                .iter()
                .filter(|(stage, _)| *stage == name)
                .map(|(_, seconds)| seconds)
                .sum::<f64>()
        };
        let taylor = stage("den") + stage("collect");
        let conversion = [
            "frames", "cut", "routing", "contract", "project", "assemble", "verify",
        ]
        .iter()
        .map(|name| stage(name))
        .sum::<f64>();
        IntegrationTimings::boundary(
            started,
            taylor,
            conversion,
            evaluate_seconds,
            finish_seconds,
        )
        .log(self.graph, id, loops, "single_mass");
        boundary::log_stats(&self.graph.name, id, &stats);
        boundary::record_stats(&self.graph.name, id, &stats);
        let elapsed = started.elapsed();
        self.integration_elapsed += elapsed;
        self.boundary_elapsed += elapsed;
        self.summary.seconds_boundary += elapsed.as_secs_f64();
        Ok(result)
    }

    /// Mode check (DESIGN 7): the boundary's value, then today's path under
    /// capture (its closures, value and root callback), compared per vacuum
    /// integral; one `BOUNDARY_CHECK` line. Returns today's value; any
    /// MISMATCH is an error.
    fn check_boundary(
        &mut self,
        remainder: &RStarIntegrand,
        built: Built,
        settings: &VakintSettings,
        id: usize,
        started: Instant,
    ) -> Result<Atom> {
        let (built, mut stats) = match built {
            Built::Zero(stats) => {
                self.summary.zero += 1;
                (None, stats)
            }
            Built::Boundary(built, stats) => {
                self.summary.boundary += 1;
                (Some(built), stats)
            }
            Built::Declined(..) => unreachable!("declines are handled by the caller"),
        };
        self.summary.contracted(&stats);
        let closed_root = self
            .closed_root
            .expect("admission requires the closed-root normalisation");
        let p_loops = remainder.region.loops.saturating_sub(1);
        let skipped =
            super::pint6_dump::settings().is_some_and(|dump| p_loops >= dump.skip_min_p_loops);
        let evaluate_started = Instant::now();
        let new_value = match &built {
            Some(_) if skipped => None,
            Some(built) => {
                let options = vakint::p_integral::BoundaryOptions {
                    strict: boundary::boundary_strict()
                        || boundary::is_large(&stats, remainder.region.loops),
                };
                let (value, reads, hits) =
                    boundary::evaluate_boundary(&built.boundary, self.vakint, settings, options)
                        .wrap_err_with(|| {
                            format!(
                                "while evaluating the remainder boundary {id} of {}",
                                self.graph.name
                            )
                        })?;
                stats.block_reads = reads;
                stats.block_hits = hits;
                Some(value)
            }
            None => Some(Atom::Zero),
        };
        stats
            .seconds
            .push(("evaluate", evaluate_started.elapsed().as_secs_f64()));
        stats.rss_gb[3] = boundary::rss_gb();
        let seconds_new = started.elapsed().as_secs_f64();
        self.summary.seconds_boundary += seconds_new;
        self.boundary_elapsed += started.elapsed();
        // Today's path, captured.
        let old_started = Instant::now();
        let ((today, callback), capture) = super::integrate::capturing_closures(|| {
            boundary::with_callback_check(closed_root, || -> Result<Atom> {
                let expression = remainder.expression()?;
                self.integrate(remainder, expression)
            })
        });
        let today = today?;
        let seconds_old = old_started.elapsed().as_secs_f64();
        let check_started = Instant::now();
        let old_value = if capture.closures.iter().all(|closure| closure.evaluated) {
            capture.value.as_ref()
        } else {
            None
        };
        let label = format!("graph={} id={id}", self.graph.name);
        let report = boundary::compare_pairs(
            built.as_ref(),
            &capture.closures,
            new_value.as_ref(),
            old_value,
            boundary::check_pair_max_bytes(),
            &label,
        );
        let seconds_check = check_started.elapsed().as_secs_f64();
        if let (Some(dump), Some(built)) = (super::pint6_dump::settings(), built.as_ref()) {
            // The boundary record beside today's v1 records, paired with them
            // (DESIGN 6b `ref_closure`); written before a mismatch is raised.
            let dump_started = Instant::now();
            self.dump_boundary(
                dump,
                built,
                &stats,
                settings,
                !skipped,
                Some(&report.ref_closures),
            )
            .wrap_err_with(|| {
                format!(
                    "while dumping the remainder boundary {id} of {}",
                    self.graph.name
                )
            })?;
            stats
                .seconds
                .push(("dump", dump_started.elapsed().as_secs_f64()));
        }
        boundary::log_stats(&self.graph.name, id, &stats);
        boundary::record_stats(&self.graph.name, id, &stats);
        eprintln!(
            "BOUNDARY_CHECK graph={} id={id} loops={} patterns={} zero_patterns={} classes={} old_closures={} pairs={} input={} value={} callback={callback} seconds_new={seconds_new:.3} seconds_old={seconds_old:.3} seconds_check={seconds_check:.3}",
            self.graph.name,
            remainder.region.loops,
            stats.patterns,
            stats.zero_patterns,
            stats.classes,
            report.old_closures,
            report.pairs,
            report.input,
            report.value,
        );
        for line in &report.failing {
            eprintln!("{line}");
        }
        self.summary.checks += 1;
        if report.input.starts_with("skipped") {
            self.summary.checks_skipped += 1;
        }
        if report.mismatch(callback) {
            self.summary.errors += 1;
            return Err(eyre!(
                "BOUNDARY_CHECK mismatch graph={} id={id}: input={} value={} callback={callback}",
                self.graph.name,
                report.input,
                report.value
            ));
        }
        Ok(today)
    }

    /// `GAMMALOOP_UV_STREAM=on|auto`: integrate the remainder's edge-basis
    /// contraction batch by batch ([`RStarIntegrand::stream_contracted`]),
    /// adding up only the integrated batches; all batches share one integration
    /// id (it labels the restored dummy indices of [`RStarIntegrand::integrate`]).
    /// `None` when the stream does not apply; the caller takes the whole path.
    fn integrate_streamed(&mut self, remainder: &RStarIntegrand) -> Result<Option<Atom>> {
        if super::expression::uv_stream() == super::expression::UvStream::Off {
            return Ok(None);
        }
        let id = self.integration_id + 1;
        let mut settings = self.settings.clone();
        settings.number_of_terms_in_epsilon_expansion = i64::try_from(remainder.region.loops)?;
        let started = Instant::now();
        let channel = super::pint6_dump::ClosureChannel {
            channel: self.channel.0,
            taylor_order: i64::from(self.plan.taylor_order),
            normalisation: self.channel.1,
            uv_child_truncate: remainder.truncation_dump_label(),
        };
        let degree = usize::try_from(self.plan.taylor_order)?;
        let (graph, root_edges, reference, auxiliary, vakint, uv_settings, projector) = (
            self.graph,
            &self.root_edges,
            self.reference_mass_squared,
            &self.auxiliary_mass_squared,
            self.vakint,
            self.uv_settings,
            self.root_projector,
        );
        let top = self.top_mass;
        let mut sum = crate::uv::approx::streamed::SumAccumulator::default();
        let streamed = remainder.stream_contracted(&mut |piece| {
            let integrated = super::pint6_dump::with_channel(channel.clone(), || {
                top_mass::in_root_scope(top, || {
                    remainder.integrate(
                        &piece,
                        graph,
                        root_edges,
                        degree,
                        reference,
                        auxiliary,
                        vakint,
                        uv_settings,
                        &settings,
                        id,
                        projector,
                    )
                })
            })?;
            sum.push(integrated);
            Ok(())
        })?;
        if !streamed {
            return Ok(None);
        }
        self.integration_id = id;
        self.integration_elapsed += started.elapsed();
        self.summary.seconds_today += started.elapsed().as_secs_f64();
        Ok(Some(sum.finish()))
    }

    fn integrate(&mut self, remainder: &RStarIntegrand, expression: &Atom) -> Result<Atom> {
        self.integration_id += 1;
        // All insertion poles and hard projector coefficients are already in
        // this expression. Only summation and the root -K follow, so request
        // poles of the complete product. The integration helper increases the
        // scalar integral's depth by the coefficient's negative valuation.
        let mut settings = self.settings.clone();
        settings.number_of_terms_in_epsilon_expansion = i64::try_from(remainder.region.loops)?;
        debug_tags!(#uv, #integrated, #rstar;
            integration_id = self.integration_id,
            loops = remainder.region.loops,
            log.expression = expression,
            "Integrating a compensated R* remainder"
        );
        let started = Instant::now();
        let channel = super::pint6_dump::ClosureChannel {
            channel: self.channel.0,
            taylor_order: i64::from(self.plan.taylor_order),
            normalisation: self.channel.1,
            uv_child_truncate: remainder.truncation_dump_label(),
        };
        let result = super::pint6_dump::with_channel(channel, || {
            top_mass::in_root_scope(self.top_mass, || {
                remainder.integrate(
                    expression,
                    self.graph,
                    &self.root_edges,
                    usize::try_from(self.plan.taylor_order)?,
                    self.reference_mass_squared,
                    &self.auxiliary_mass_squared,
                    self.vakint,
                    self.uv_settings,
                    &settings,
                    self.integration_id,
                    self.root_projector,
                )
            })
        });
        self.integration_elapsed += started.elapsed();
        self.summary.seconds_today += started.elapsed().as_secs_f64();
        result
    }

    /// The census record of a built boundary (DESIGN 6b). `references`: the
    /// v1 records of today's path (check mode); in mode on the v1 records
    /// are written from the boundary first when `PINT6_BOUNDARY_V1=1` (6c),
    /// with this outer integration's channel, as today's path writes them.
    fn dump_boundary(
        &self,
        dump: &super::pint6_dump::DumpSettings,
        built: &boundary::BuiltBoundary,
        stats: &boundary::BoundaryStats,
        settings: &VakintSettings,
        evaluated: bool,
        references: Option<&BTreeMap<usize, usize>>,
    ) -> Result<()> {
        let references = match references {
            Some(references) => references.clone(),
            None if super::pint6_dump::boundary_v1() => {
                let channel = super::pint6_dump::ClosureChannel {
                    channel: self.channel.0,
                    taylor_order: i64::from(self.plan.taylor_order),
                    normalisation: self.channel.1,
                    uv_child_truncate: super::expression::untruncated_dump_label(),
                };
                super::pint6_dump::with_channel(channel, || {
                    super::pint6_dump::write_v1_from_boundary(dump, settings, built)
                })?
            }
            None => BTreeMap::new(),
        };
        let counts = boundary::stats_json(&self.graph.name, built.integration_id, stats);
        let written = super::pint6_dump::write_boundary(
            dump,
            settings,
            built,
            &counts,
            evaluated,
            &references,
        )?;
        eprintln!(
            "BOUNDARY_DUMP graph={} id={} sequence={} header={} block_files={} text_bytes={} ref_closures={} rss_gb={:.2}",
            self.graph.name,
            built.integration_id,
            written.sequence,
            written.header.display(),
            written.block_files,
            written.text_bytes,
            references.len(),
            boundary::rss_gb()
        );
        if super::pint6_dump::boundary_verify() {
            let started = Instant::now();
            super::pint6_dump::verify_boundary(&built.boundary, &written.header)?;
            eprintln!(
                "BOUNDARY_DUMP_VERIFY graph={} id={} sequence={} identical seconds={:.3} rss_gb={:.2}",
                self.graph.name,
                built.integration_id,
                written.sequence,
                started.elapsed().as_secs_f64(),
                boundary::rss_gb()
            );
        }
        Ok(())
    }
}

/// The momentum degree of a homogeneous insertion polynomial: the common
/// number of edge momenta (`Q`) in its expanded summands; `None` when the
/// summands differ in degree.
fn insertion_degree(insertion: &Atom) -> Option<i32> {
    use symbolica::atom::AtomView;
    fn count(view: AtomView<'_>) -> Option<i64> {
        match view {
            AtomView::Fun(function) if function.get_symbol() == GS.emr_mom => Some(1),
            AtomView::Fun(function) => {
                let mut total = 0;
                for argument in function.iter() {
                    total += count(argument)?;
                }
                Some(total)
            }
            AtomView::Mul(product) => {
                let mut total = 0;
                for factor in product.iter() {
                    total += count(factor)?;
                }
                Some(total)
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let inner = count(base)?;
                if inner == 0 {
                    Some(0)
                } else {
                    Some(inner * i64::try_from(exponent).ok()?)
                }
            }
            AtomView::Add(sum) => {
                let mut degrees = sum
                    .iter()
                    .map(count)
                    .collect::<Option<std::collections::BTreeSet<_>>>()?;
                degrees.remove(&-1);
                (degrees.len() <= 1).then(|| degrees.into_iter().next().unwrap_or(0))
            }
            _ => Some(0),
        }
    }
    let expanded = insertion.expand();
    count(expanded.as_view()).and_then(|degree| i32::try_from(degree).ok())
}
