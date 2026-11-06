//! Assemble the ordinary UV preparation and the UV form of the IR correction.

use std::collections::BTreeMap;

use eyre::{Result, WrapErr, ensure, eyre};
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubSetLike};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    parse, symbol,
};
use vakint::{Vakint, VakintSettings};

use crate::{
    debug_tags,
    graph::Graph,
    model::Model,
    uv::{UVgenerationSettings, UltravioletGraph, uv_graph::UVE},
};

use super::{
    CountertermCall, IrForestTerm, MassiveUvEvaluator, RStarPlan, UvPreparationTerm,
    expression::RStarIntegrand, moments::hard_moments,
};

impl RStarPlan {
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
        Evaluation::new(self, graph, model, vakint, uv_settings)?.counterterm()
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
        let mut evaluation = Evaluation::new(self, graph, model, vakint, uv_settings)?;
        let mut result = -evaluation.counterterm()?;
        for term in self
            .ordinary_uv
            .iter()
            .filter(|term| !term.counterterms.is_empty())
        {
            let insertions = evaluation
                .insertions(term)?
                .into_iter()
                .map(|insertion| insertion.replace(evaluation.reference_mass_squared).with(0))
                .collect::<Vec<_>>();
            let remainder = RStarIntegrand::new(
                graph,
                model,
                &term.remainder,
                &insertions,
                self.marked_edge,
                Atom::Zero,
            )?;
            evaluation.integration_id += 1;
            let mut settings = evaluation.settings.clone();
            settings.number_of_terms_in_epsilon_expansion = i64::try_from(term.remainder.loops)?;
            result -= remainder.integrate_physical(
                vakint,
                uv_settings,
                &settings,
                evaluation.integration_id,
            )?;
        }
        evaluation.poles(&result.replace_multiple(&external_routing))
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
    uv: MassiveUvEvaluator<'a>,
    quotient_counterterms: BTreeMap<(CountertermCall, Atom), Atom>,
    integration_id: usize,
}

impl<'a> Evaluation<'a> {
    fn new(
        plan: &'a RStarPlan,
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        uv_settings: &'a UVgenerationSettings,
    ) -> Result<Self> {
        for (_, _, edge) in graph.iter_edges_of(&plan.support) {
            ensure!(
                edge.data.mass_atom(model).is_zero(),
                "R* currently requires a massless input graph"
            );
        }
        let settings = uv_settings.vakint.true_settings();
        let epsilon = parse!(&settings.epsilon_symbol)
            .get_symbol()
            .ok_or_else(|| eyre!("Vakint epsilon must be a symbol"))?;
        let reference_mass_squared = symbol!("rstar_reference_mass_squared"; Scalar);
        let reference_mass = Atom::var(reference_mass_squared).sqrt();
        let children = if plan.counterterms.is_empty() {
            BTreeMap::new()
        } else {
            plan.evaluate_uv_children(graph, model, &reference_mass, vakint, uv_settings)
                .wrap_err("while computing the proper UV counterterms with massive R")?
        };
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
            uv: MassiveUvEvaluator::new(graph, model, vakint, uv_settings),
            quotient_counterterms: BTreeMap::new(),
            integration_id: 0,
        })
    }

    fn counterterm(&mut self) -> Result<Atom> {
        if self.graph.compute_dod(&self.plan.support) < 0 {
            return Ok(Atom::Zero);
        }
        let mut prepared = Atom::Zero;
        for term in &self.plan.ordinary_uv {
            if self.scaleless(term) {
                continue;
            }
            let remainder = self.remainder(term)?;
            prepared += self.integrate(&remainder, &remainder.expression)?;
        }
        for correction in &self.plan.ir {
            for term in &correction.preparation {
                if self.scaleless(term) {
                    continue;
                }
                let hard = self.remainder(term)?;
                let polynomial = hard
                    .taylor(usize::try_from(correction.taylor_order)?)
                    .wrap_err("while Taylor-expanding the prepared hard component")?;
                let covariant = hard.project_hard(&polynomial)?;
                let moments = hard_moments(&covariant, hard.soft_edges())?;
                for forest in &correction.forest {
                    let mut expression = Atom::Zero;
                    for (coefficient, insertion) in &moments {
                        expression += coefficient * self.ir_forest(forest, insertion)?;
                    }
                    if !expression.is_zero() {
                        prepared += self.integrate(&hard, &expression)?;
                    }
                }
            }
        }
        Ok(-self.poles(&prepared)?)
    }

    fn poles(&self, expression: &Atom) -> Result<Atom> {
        Ok(expression.series(self.epsilon, 0, -1)?.to_atom().expand())
    }

    fn scaleless(&self, term: &UvPreparationTerm) -> bool {
        // A local mixed contribution can retain the reference mass after its
        // marked edge is contracted. Only unintegrated massless vacuum loops
        // vanish here, not a local polynomial or a UV counterterm operation.
        term.remainder.loops > 0
            && !term
                .remainder
                .remaining
                .includes(&self.graph[&self.plan.marked_edge].1)
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

    fn remainder(&self, term: &UvPreparationTerm) -> Result<RStarIntegrand> {
        let insertions = self.insertions(term)?;
        RStarIntegrand::new(
            self.graph,
            self.model,
            &term.remainder,
            &insertions,
            self.plan.marked_edge,
            Atom::var(self.reference_mass_squared),
        )
    }

    fn quotient(&mut self, call: &CountertermCall, insertion: &Atom) -> Result<Atom> {
        let key = (call.clone(), insertion.clone());
        if let Some(value) = self.quotient_counterterms.get(&key) {
            return Ok(value.clone());
        }
        let result = self
            .uv
            .counterterm(call, std::slice::from_ref(insertion))
            .wrap_err_with(|| {
                format!(
                    "while applying massive R to the decorated quotient {}",
                    call.argument.remaining.string_label()
                )
            })?;
        self.quotient_counterterms.insert(key, result.clone());
        Ok(result)
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
        remainder.integrate(
            expression,
            self.graph,
            usize::try_from(self.plan.taylor_order)?,
            self.reference_mass_squared,
            &self.auxiliary_mass_squared,
            self.vakint,
            self.uv_settings,
            &settings,
            self.integration_id,
        )
    }
}
