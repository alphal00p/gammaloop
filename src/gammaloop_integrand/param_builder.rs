use std::{borrow::Cow, fmt::Display, ops::Deref};

use bincode_trait_derive::{Decode, Encode};
use idenso::color::CS;
use itertools::Itertools;
use linnet::half_edge::involution::{HedgePair, Orientation};
use log::debug;
use spenso::{
    algebra::{algebraic_traits::RefOne, complex::Complex},
    iterators::IteratableTensor,
    network::ExecutionResult,
    structure::concrete_index::ExpandedIndex,
    tensors::parametric::AtomViewOrConcrete,
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    domains::rational::Rational,
    evaluate::FunctionMap,
    id::Replacement,
    symbol,
};
use tabled::settings::Style;

use crate::{
    cff::expression::GraphOrientation,
    graph::{FeynmanGraph, Graph},
    model::Model,
    momentum::{Helicity, PolType},
    momentum_sample::{ExternalFourMomenta, MomentumSample},
    numerator::ParsingNet,
    utils::{f128, FloatLike, ToCoefficient, F, GS, TENSORLIB},
    GammaLoopContext,
};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ParamValuePairs<T: FloatLike> {
    pub values: Vec<Complex<F<T>>>,
    pub params: Vec<Atom>,
}

impl<T: FloatLike> ParamValuePairs<T> {
    pub fn validate(&self) {
        assert_eq!(
            self.values.len(),
            self.params.len(),
            "Number of values and parameters must match"
        );
    }

    pub fn map_values<U: FloatLike>(
        &self,
        mut map: impl FnMut(&Complex<F<T>>) -> Complex<F<U>>,
    ) -> ParamValuePairs<U> {
        let values = self.values.iter().map(&mut map).collect();
        ParamValuePairs {
            values,
            params: self.params.clone(),
        }
    }
    pub fn replacement(&self) -> Vec<Replacement> {
        let mut replacements = Vec::new();
        for (p, v) in self.params.iter().zip_eq(self.values.iter()) {
            println!("{p}");
            let replacement = Replacement::new(
                p.clone().to_pattern(),
                Atom::num(v.clone().to_coefficient()),
            );
            replacements.push(replacement);
        }
        replacements
    }

    pub fn extract_and_fill(&mut self, values: &mut Vec<Complex<F<T>>>) {
        self.values = values.split_off(self.params.len())
    }

    pub fn default_from_symbol(symbol: Symbol) -> Self {
        Self {
            values: vec![Complex::new_re(F::from_f64(0.0))],
            params: vec![Atom::var(symbol)],
        }
    }

    pub fn set_zeros(&mut self) {
        self.values = vec![Complex::new_re(F::from_f64(0.0)); self.params.len()];
    }
}

impl<T: FloatLike> Default for ParamValuePairs<T> {
    fn default() -> Self {
        Self {
            values: Vec::new(),
            params: Vec::new(),
        }
    }
}

// impl<T:FloatLike> ParamValuePairs<T>{
//     pub fn map_values<U:FloatLike>(&self,map:impl FnMut(&Complex<F<T>>)->Complex<F<U>>)->ParamValuePairs<U>{

//     }
// }

impl<T: FloatLike> ParamValuePairs<T>
where
    T::Higher: FloatLike,
    T::Lower: FloatLike,
{
    fn higher(&self) -> ParamValuePairs<T::Higher> {
        ParamValuePairs {
            values: self
                .values
                .iter()
                .map(|v| v.map_ref(|v| v.higher()))
                .collect(),
            params: self.params.clone(),
        }
    }
    fn lower(&self) -> ParamValuePairs<T::Lower> {
        ParamValuePairs {
            values: self
                .values
                .iter()
                .map(|v| v.map_ref(|v| v.lower()))
                .collect(),
            params: self.params.clone(),
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ParamBuilder<T: FloatLike = f64> {
    // values: Vec<Complex<F<T>>
    m_uv: ParamValuePairs<T>,
    idenso_vars: ParamValuePairs<T>,
    mu_r_sq: ParamValuePairs<T>,
    orientations: ParamValuePairs<T>,
    pub model_parameters: ParamValuePairs<T>,
    pub external_energies: ParamValuePairs<T>,
    pub external_spatial: ParamValuePairs<T>,
    polarizations: ParamValuePairs<T>,
    emr_spatial: ParamValuePairs<T>,
    tstar: ParamValuePairs<T>,
    h_function: ParamValuePairs<T>,
    esurface_derivative: ParamValuePairs<T>,
    uv_damp_plus: ParamValuePairs<T>,
    uv_damp_minus: ParamValuePairs<T>,
    radius: ParamValuePairs<T>,
    radius_star: ParamValuePairs<T>,
    pub reps: Vec<(Atom, Atom)>,
    // pub eager_const_map: HashMap<Atom, Complex<F<T>>>,
    // pub eager_function_map: HashMap<Symbol, EvaluationFn<Atom, Complex<F<T>>>>,
    // pub eager_fn_map:
    pub fn_map: FunctionMap,
}

pub struct ThresholdParams<T: FloatLike> {
    pub radius: F<T>,
    pub radius_star: F<T>,
    pub esurface_derivative: F<T>,
    pub uv_damp_plus: F<T>,
    pub uv_damp_minus: F<T>,
    pub h_function: F<T>,
}

pub trait UpdateAndGetParams<T: FloatLike> {
    fn update_emr_and_get_params<'a>(
        &'a mut self,
        sample: &'a MomentumSample<T>,
        graph: &'a Graph,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<T>>,
    ) -> Cow<'a, Vec<Complex<F<T>>>>;
}

impl UpdateAndGetParams<f64> for ParamBuilder<f64> {
    fn update_emr_and_get_params<'a>(
        &'a mut self,
        sample: &'a MomentumSample<f64>,
        graph: &'a Graph,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<f64>>,
    ) -> Cow<'a, Vec<Complex<F<f64>>>> {
        let emr_spatial: Vec<_> = graph
            .iter_edges()
            .flat_map(|(pair, edge_id, _)| {
                if let HedgePair::Paired { .. } = pair {
                    let emr_vec = graph.loop_momentum_basis.edge_signatures[edge_id]
                        .compute_three_momentum_from_four(
                            sample.loop_moms(),
                            sample.external_moms(),
                        );
                    // println!("{edge_id}:{emr_vec}");
                    vec![
                        Complex::new_re(emr_vec.px),
                        Complex::new_re(emr_vec.py),
                        Complex::new_re(emr_vec.pz),
                    ]
                } else {
                    vec![]
                }
            })
            .collect();

        // parse!("s").evaluator(fn_map, params, optimization_settings).unwrap().

        self.external_spatial_value(sample);
        self.external_energies_value(sample);
        self.emr_spatial.values = emr_spatial;
        self.polarizations_values(graph, sample.external_moms(), helicities);

        if let Some(threshold_params) = threshold_params {
            self.threshold_params(threshold_params);
        }

        debug!("params:\n {}", self);

        Cow::Owned(self.clone().build_values()) // ideally borrows a single vec
    }
}

impl UpdateAndGetParams<f128> for ParamBuilder<f64> {
    fn update_emr_and_get_params(
        &mut self,
        sample: &MomentumSample<f128>,
        graph: &Graph,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<f128>>,
    ) -> Cow<Vec<Complex<F<f128>>>> {
        let emr_spatial: Vec<_> = graph
            .iter_edges()
            .flat_map(|(pair, edge_id, _)| {
                if let HedgePair::Paired { .. } = pair {
                    let emr_vec = graph.loop_momentum_basis.edge_signatures[edge_id]
                        .compute_three_momentum_from_four(
                            sample.loop_moms(),
                            sample.external_moms(),
                        );
                    vec![
                        Complex::new_re(emr_vec.px),
                        Complex::new_re(emr_vec.py),
                        Complex::new_re(emr_vec.pz),
                    ]
                } else {
                    vec![]
                }
            })
            .collect();

        let mut new_param_builder = self.higher();
        new_param_builder.emr_spatial.values = emr_spatial;
        new_param_builder.external_spatial_value(sample);
        new_param_builder.external_energies_value(sample);
        new_param_builder.polarizations_values(graph, sample.external_moms(), helicities);

        if let Some(threshold_params) = threshold_params {
            new_param_builder.threshold_params(threshold_params);
        }

        // println!("ParamBuilder before eval f128:\n{}", self);

        Cow::Owned(new_param_builder.build_values())
    }
}

impl<T: FloatLike> ParamBuilder<T>
where
    T::Higher: FloatLike,
    T::Lower: FloatLike,
{
    fn higher(&self) -> ParamBuilder<T::Higher> {
        ParamBuilder {
            fn_map: self.fn_map.clone(),
            reps: self.reps.clone(),
            m_uv: self.m_uv.higher(),
            mu_r_sq: self.mu_r_sq.higher(),
            idenso_vars: self.idenso_vars.higher(),
            orientations: self.orientations.higher(),
            model_parameters: self.model_parameters.higher(),
            external_energies: self.external_energies.higher(),
            external_spatial: self.external_spatial.higher(),
            polarizations: self.polarizations.higher(),
            emr_spatial: self.emr_spatial.higher(),
            tstar: self.tstar.higher(),
            h_function: self.h_function.higher(),
            esurface_derivative: self.esurface_derivative.higher(),
            uv_damp_plus: self.uv_damp_plus.higher(),
            uv_damp_minus: self.uv_damp_minus.higher(),
            radius: self.radius.higher(),
            radius_star: self.radius_star.higher(),
        }
    }
    fn lower(&self) -> ParamBuilder<T::Lower> {
        ParamBuilder {
            orientations: self.orientations.lower(),
            fn_map: self.fn_map.clone(),
            reps: self.reps.clone(),
            m_uv: self.m_uv.lower(),
            idenso_vars: self.idenso_vars.lower(),
            mu_r_sq: self.mu_r_sq.lower(),
            model_parameters: self.model_parameters.lower(),
            external_energies: self.external_energies.lower(),
            external_spatial: self.external_spatial.lower(),
            polarizations: self.polarizations.lower(),
            emr_spatial: self.emr_spatial.lower(),
            tstar: self.tstar.lower(),
            h_function: self.h_function.lower(),
            esurface_derivative: self.esurface_derivative.lower(),
            uv_damp_plus: self.uv_damp_plus.lower(),
            uv_damp_minus: self.uv_damp_minus.lower(),
            radius: self.radius.lower(),
            radius_star: self.radius_star.lower(),
        }
    }
}
impl<T: FloatLike> ParamBuilder<T> {
    pub fn validate(&self) {
        debug!("Validating mu_r_sq");
        self.mu_r_sq.validate();
        debug!("Validating model_parameters");
        self.model_parameters.validate();
        debug!("Validating external_energies");
        self.external_energies.validate();
        debug!("Validating external_spatial");
        self.external_spatial.validate();
        debug!("Validating polarizations");
        self.polarizations.validate();
        debug!("Validating emr_spatial");
        self.emr_spatial.validate();
        debug!("Validating tstar");
        self.tstar.validate();
        debug!("Validating h_function");
        self.h_function.validate();
        debug!("Validating derivative_at_tstar");
        self.esurface_derivative.validate();
        debug!("Validating uv_damp_plus");
        self.uv_damp_plus.validate();
        debug!("Validating uv_damp_minus");
        self.uv_damp_minus.validate();
        debug!("Validating radius");
        self.radius.validate();
        debug!("Validating radius_star");
        self.radius_star.validate();
    }

    pub fn fill_in_values(&mut self, mut values: Vec<Complex<F<T>>>) {
        self.m_uv.extract_and_fill(&mut values);
        self.mu_r_sq.extract_and_fill(&mut values);
        self.model_parameters.extract_and_fill(&mut values);
        self.external_energies.extract_and_fill(&mut values);
        self.external_spatial.extract_and_fill(&mut values);
        self.polarizations.extract_and_fill(&mut values);
        self.emr_spatial.extract_and_fill(&mut values);
        self.tstar.extract_and_fill(&mut values);
        self.h_function.extract_and_fill(&mut values);
        self.esurface_derivative.extract_and_fill(&mut values);
        self.uv_damp_plus.extract_and_fill(&mut values);
        self.uv_damp_minus.extract_and_fill(&mut values);
        self.radius.extract_and_fill(&mut values);
        self.radius_star.extract_and_fill(&mut values);
    }

    pub fn replace_non_emr(&self, atom: impl AtomCore) -> Atom {
        let reps = self
            .reps
            .iter()
            .map(|(a, b)| Replacement::new(a.clone().to_pattern(), b.clone()))
            .collect_vec();
        let mut a = atom.replace_multiple(&reps);

        for r in [
            &self.m_uv,
            &self.mu_r_sq,
            &self.model_parameters,
            &self.external_energies,
            &self.external_spatial,
            // &self.polarizations,
            // &self.emr_spatial,
            &self.tstar,
            &self.h_function,
            &self.esurface_derivative,
            &self.uv_damp_plus,
            &self.uv_damp_minus,
            &self.radius,
            &self.radius_star,
        ]
        .into_iter()
        {
            a = a.replace_multiple(&r.replacement());
        }
        a
    }

    pub fn add_tagged_function(
        &mut self,
        name: Symbol,
        tags: Vec<Atom>,
        rename: String,
        args: Vec<Symbol>,
        body: Atom,
    ) -> Result<(), &str> {
        self.reps.push((
            FunctionBuilder::new(name)
                .add_args(&tags)
                .add_args(&args)
                .finish(),
            body.clone(),
        ));

        self.fn_map
            .add_tagged_function(name, tags, rename, args, body)

        // body.evaluate(coeff_map, const_map, function_map)
    }

    pub fn add_function(
        &mut self,
        name: Symbol,
        rename: String,
        args: Vec<Symbol>,
        body: Atom,
    ) -> Result<(), &str> {
        self.reps.push((
            FunctionBuilder::new(name)
                // .add_args()
                .add_args(&args)
                .finish(),
            body.clone(),
        ));
        self.fn_map.add_function(name, rename, args, body)
    }

    pub fn add_constant(&mut self, key: Atom, value: symbolica::domains::float::Complex<Rational>) {
        self.fn_map.add_constant(key, value)
    }

    pub(crate) fn build_values(self) -> Vec<Complex<F<T>>> {
        let mut values = Vec::with_capacity(100);
        for value in self.into_iter() {
            values.extend(value.values);
        }
        values
    }
    pub(crate) fn build_params(self) -> Vec<Atom> {
        let mut values = Vec::with_capacity(100);
        for value in self.into_iter() {
            values.extend(value.params);
        }
        values
    }

    pub(crate) fn new_empty() -> Self {
        Self {
            fn_map: FunctionMap::default(),
            idenso_vars: ParamValuePairs::default(),
            orientations: ParamValuePairs::default(),
            m_uv: ParamValuePairs::default(),
            mu_r_sq: ParamValuePairs::default(),
            model_parameters: ParamValuePairs::default(),
            external_energies: ParamValuePairs::default(),
            polarizations: ParamValuePairs::default(),
            external_spatial: ParamValuePairs::default(),
            emr_spatial: ParamValuePairs::default(),
            tstar: ParamValuePairs::default(),
            h_function: ParamValuePairs::default(),
            esurface_derivative: ParamValuePairs::default(),
            uv_damp_plus: ParamValuePairs::default(),
            uv_damp_minus: ParamValuePairs::default(),
            radius: ParamValuePairs::default(),
            radius_star: ParamValuePairs::default(),
            reps: Vec::new(),
        }
    }

    pub(crate) fn idenso_vars(&mut self) {
        self.idenso_vars.params = [CS.tr, CS.nc].into_iter().map(Atom::var).collect();

        self.idenso_vars.values = vec![
            Complex::new_re(F(T::from_f64(0.5))),
            Complex::new_re(F(T::from_f64(3.))),
        ];
    }

    pub(crate) fn new(graph: &Graph, model: &Model) -> Self {
        let mut new = Self::new_empty();

        new.m_uv = ParamValuePairs::default_from_symbol(GS.m_uv);

        new.mu_r_sq = ParamValuePairs::default_from_symbol(GS.mu_r_sq);
        new.tstar = ParamValuePairs::default_from_symbol(GS.rescale_star);
        new.radius = ParamValuePairs::default_from_symbol(GS.radius);
        new.radius_star = ParamValuePairs::default_from_symbol(GS.radius_star);
        new.h_function = ParamValuePairs::default_from_symbol(GS.hfunction);
        new.esurface_derivative = ParamValuePairs::default_from_symbol(GS.deta);
        new.uv_damp_plus = ParamValuePairs::default_from_symbol(GS.uv_damp_plus);
        new.uv_damp_minus = ParamValuePairs::default_from_symbol(GS.uv_damp_minus);
        new.idenso_vars();

        new.update_model_data(model);
        new.external_energies_atom(graph);
        new.orientation_params(graph);
        new.polarization_params(graph);
        new.external_spatial_atom(graph);
        new.emr_spatial_atom(graph);

        let pi_rational = Rational::from(std::f64::consts::PI);

        for (_, e, _) in graph.iter_edges() {
            new.add_tagged_function(
                GS.ose,
                vec![Atom::num(e.0 as i64)],
                format!("OSE{e}"),
                vec![],
                graph.explicit_ose_atom(e),
            )
            .unwrap();
        }
        new.add_constant(Atom::PI.into(), pi_rational.into());

        new
    }

    pub(crate) fn m_uv_value(&mut self, m_uv: Complex<F<T>>) {
        self.m_uv.values = vec![m_uv];
    }

    pub(crate) fn mu_r_sq_value(&mut self, mu_r_sq: Complex<F<T>>) {
        self.mu_r_sq.values = vec![mu_r_sq];
    }

    pub(crate) fn update_model_data(&mut self, model: &Model) {
        self.model_parameters.params = model.generate_params();
        self.model_parameters.values = model.generate_values();
    }

    pub(crate) fn external_energies_atom(&mut self, graph: &Graph) {
        self.external_energies.params = graph.get_external_energy_atoms();
        self.external_energies.set_zeros();
    }

    pub(crate) fn add_external_four_mom(&mut self, ext: &ExternalFourMomenta<F<T>>) {
        self.external_energies.values = ext
            .iter()
            .map(|x| Complex::new_re(x.temporal.value.clone()))
            .collect_vec();
        self.external_spatial.values = ext
            .iter()
            .flat_map(|x| x.spatial.clone().into_iter().map(|c| Complex::new_re(c)))
            .collect_vec();
    }

    pub(crate) fn polarization_params(&mut self, graph: &Graph) {
        // println!("{}");c
        let mut pols = graph.global_prefactor.polarizations();
        pols.sort_by(|a, b| a.0.cmp(&b.0));
        let mut params = Vec::new();

        for (p, a) in pols {
            // println!("{a}");
            match ParsingNet::try_from_view(a.as_view(), TENSORLIB.read().unwrap().deref())
                .unwrap()
                .result_tensor(TENSORLIB.read().unwrap().deref())
                .unwrap()
            {
                ExecutionResult::One => {}
                ExecutionResult::Zero => {}
                ExecutionResult::Val(a) => {
                    for (_, val) in a.iter_flat() {
                        let AtomViewOrConcrete::Atom(a) = val else {
                            panic!("SHOULD BE ATOMVIEW")
                        };

                        params.push(a.to_owned());
                    }
                }
            }
        }

        self.polarizations.params = params;
        self.polarizations.set_zeros();

        // self.polarizations.params = graph.generate_polarization_params();
    }

    pub(crate) fn orientation_params(&mut self, graph: &Graph) {
        let mut params = Vec::new();

        for (p, i, d) in graph.iter_edges() {
            params.push(GS.sign(i));
            params.push(GS.sign_theta(GS.sign(i)));
            params.push(GS.sign_theta(-GS.sign(i)));
        }

        self.orientations.params = params;
        self.orientations.set_zeros();
    }

    pub(crate) fn orientation_value<O: GraphOrientation>(&mut self, orientation: &O) {
        let mut values = Vec::new();
        let zero: Complex<F<T>> = Complex::new_re(F(T::from_f64(0.)));
        let one = zero.ref_one();
        let minusone = -(one.clone());

        for (_, i) in orientation.orientation() {
            match i {
                Orientation::Default => {
                    values.push(one.clone());
                    values.push(one.clone());
                    values.push(zero.clone());
                }
                Orientation::Reversed => {
                    values.push(minusone.clone());
                    values.push(zero.clone());
                    values.push(one.clone());
                }
                Orientation::Undirected => {
                    values.push(zero.clone());
                    values.push(zero.clone());
                    values.push(zero.clone());
                }
            }
        }

        self.orientations.values = values;
    }

    pub(crate) fn polarizations_values(
        &mut self,
        graph: &Graph,
        ext: &ExternalFourMomenta<F<T>>,
        helicities: &[Helicity],
    ) {
        let mut pols = graph.global_prefactor.polarizations();
        pols.sort_by(|a, b| a.0.cmp(&b.0));

        let mut vals = Vec::new();

        for (p, a) in pols {
            let extid = graph.loop_momentum_basis.ext_from(p.eid).unwrap();
            let hel = p.hel.unwrap_or(helicities[extid.0]);
            // println!("MOM:{}", ext[extid]);
            // println!("Pol {},{},{:?}", extid, hel, p.pol_type);
            let pol = match p.pol_type {
                PolType::Epsilon => ext[extid].pol(hel),
                PolType::EpsilonBar => ext[extid].pol(hel).bar(),
                PolType::Scalar => {
                    continue;
                }
                PolType::U => ext[extid].u(hel.try_into().unwrap()),
                PolType::V => ext[extid].v(hel.try_into().unwrap()),
                PolType::UBar => ext[extid].u(hel.try_into().unwrap()).bar(),
                PolType::VBar => ext[extid].v(hel.try_into().unwrap()).bar(),
            };

            for (_, val) in pol.tensor.iter_flat() {
                // println!("{val}");
                vals.push(val.clone());
            }
        }

        self.polarizations.values = vals;

        // self.polarizations.params = graph.generate_polarization_params();
    }

    pub(crate) fn external_energies_value(&mut self, momentum_sample: &MomentumSample<T>) {
        self.external_energies.values = momentum_sample
            .external_moms()
            .iter()
            .map(|x| Complex::new_re(x.temporal.value.clone()))
            .collect_vec();
    }

    pub(crate) fn external_spatial_atom(&mut self, graph: &Graph) {
        self.external_spatial.params = graph
            .iter_edges()
            .flat_map(|(pair, edge_id, _)| {
                if let HedgePair::Unpaired { .. } = pair {
                    vec![
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                    ]
                } else {
                    vec![]
                }
            })
            .collect();
    }

    pub(crate) fn external_spatial_value(&mut self, momentum_sample: &MomentumSample<T>) {
        self.external_spatial.values = momentum_sample
            .external_moms()
            .iter()
            .flat_map(|x| x.spatial.clone().into_iter().map(|c| Complex::new_re(c)))
            .collect_vec();
    }

    pub(crate) fn emr_spatial_atom(&mut self, graph: &Graph) {
        self.emr_spatial.params = graph
            .iter_edges()
            .flat_map(|(pair, edge_id, _)| {
                if let HedgePair::Paired { .. } = pair {
                    vec![
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                    ]
                } else {
                    vec![]
                }
            })
            .collect();
    }

    pub(crate) fn emr_spatial_value(&mut self, emr_spatial: Vec<Complex<F<T>>>) {
        self.emr_spatial.values = emr_spatial;
    }

    pub(crate) fn tstar_value(&mut self, tstar: Complex<F<T>>) {
        self.tstar.values = vec![tstar];
    }

    pub(crate) fn h_function_value(&mut self, h_function: Complex<F<T>>) {
        self.h_function.values = vec![h_function];
    }

    pub(crate) fn derivative_at_tstar_value(&mut self, derivative_at_tstar: Complex<F<T>>) {
        self.esurface_derivative.values = vec![derivative_at_tstar];
    }

    pub(crate) fn radius_star_atom(&mut self, radius_star: Atom) {
        self.radius_star.params = vec![radius_star];
    }

    pub(crate) fn threshold_params(&mut self, threshold_params: &ThresholdParams<T>) {
        self.radius.values = vec![Complex::new_re(threshold_params.radius.clone())];
        self.radius_star.values = vec![Complex::new_re(threshold_params.radius_star.clone())];
        self.esurface_derivative.values = vec![Complex::new_re(
            threshold_params.esurface_derivative.clone(),
        )];
        self.uv_damp_plus.values = vec![Complex::new_re(threshold_params.uv_damp_plus.clone())];
        self.uv_damp_minus.values = vec![Complex::new_re(threshold_params.uv_damp_minus.clone())];
        self.h_function.values = vec![Complex::new_re(threshold_params.h_function.clone())];
    }
}

impl<T: FloatLike> Display for ParamBuilder<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = tabled::builder::Builder::new();

        for (lhs, rhs) in &self.reps {
            table.push_record(vec![lhs.to_string(), rhs.to_string()]);
        }

        for i in self {
            if i.values.is_empty() {
                for p in i.params.iter() {
                    table.push_record(vec![p.to_string(), "N/A".to_string()]);
                }
            } else {
                for (v, p) in i.values.iter().zip(i.params.iter()) {
                    table.push_record(vec![p.to_string(), v.to_string()]);
                }
            }
        }

        table.build().with(Style::rounded()).to_string().fmt(f)
    }
}

impl<T: FloatLike> IntoIterator for ParamBuilder<T> {
    type Item = ParamValuePairs<T>;
    type IntoIter = std::array::IntoIter<Self::Item, 16>;

    fn into_iter(self) -> Self::IntoIter {
        [
            self.m_uv,
            self.mu_r_sq,
            self.idenso_vars,
            self.model_parameters,
            self.external_energies,
            self.external_spatial,
            self.polarizations,
            self.emr_spatial,
            self.tstar,
            self.h_function,
            self.esurface_derivative,
            self.uv_damp_plus,
            self.uv_damp_minus,
            self.radius,
            self.radius_star,
            self.orientations,
        ]
        .into_iter()
    }
}

impl<'a, T: FloatLike> IntoIterator for &'a ParamBuilder<T> {
    type Item = &'a ParamValuePairs<T>;
    type IntoIter = std::array::IntoIter<Self::Item, 16>;

    fn into_iter(self) -> Self::IntoIter {
        [
            &self.m_uv,
            &self.mu_r_sq,
            &self.idenso_vars,
            &self.model_parameters,
            &self.external_energies,
            &self.external_spatial,
            &self.polarizations,
            &self.emr_spatial,
            &self.tstar,
            &self.h_function,
            &self.esurface_derivative,
            &self.uv_damp_plus,
            &self.uv_damp_minus,
            &self.radius,
            &self.radius_star,
            &self.orientations,
        ]
        .into_iter()
    }
}

#[test]
fn evaltest() {
    use symbolica::evaluate::{FunctionMap, OptimizationSettings};
    use symbolica::{atom::AtomCore, parse};
    let expr = parse!("delta_sigma(1,1,1,1)");
    let args: Vec<_> = ["x1", "x2", "x3", "x4"]
        .into_iter()
        .map(|a| symbol!(a))
        .collect();
    let body = parse!("theta(x1*x(1))*theta(x2*x(2))*theta(x3*x(3))*theta(x4*x(4))");

    let mut fn_map = FunctionMap::new();

    let params: Vec<Atom> = vec!["theta(x(1))", "theta(x(2))", "theta(x(3))", "theta(x(4))"]
        .iter()
        .map(|a| parse!(a))
        .collect();
    fn_map.add_function(symbol!("delta_sigma"), "delta_sigma".into(), args, body);
    let optimization_settings = OptimizationSettings::default();
    let mut evaluator = expr
        .evaluator(&fn_map, &params, optimization_settings)
        .unwrap()
        .map_coeff(&|x| x.to_real().unwrap().to_f64());
    assert_eq!(evaluator.evaluate_single(&[1.0, 1.0, 4.0, 4.0]), 8.0);
}
