use std::{
    borrow::Cow,
    fmt::Display,
    ops::{Deref, Range},
};

use bincode::{Decode, Encode};
use idenso::color::CS;
use libc::ENOPROTOOPT;
use linnet::half_edge::involution::{EdgeIndex, HedgePair, Orientation};
use log::debug;
use ringbuffer::{ConstGenericRingBuffer, RingBuffer};
use spenso::{
    algebra::{algebraic_traits::RefOne, complex::Complex},
    iterators::IteratableTensor,
    network::ExecutionResult,
    tensors::parametric::AtomViewOrConcrete,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, FunctionBuilder, Symbol},
    domains::rational::Rational,
    evaluate::FunctionMap,
    id::Replacement,
};
use tabled::{settings::Style, Table};

use crate::{
    cff::expression::GraphOrientation,
    graph::Graph,
    model::Model,
    momentum::{Helicity, PolType},
    momentum_sample::{ExternalFourMomenta, MomentumSample},
    numerator::ParsingNet,
    status_debug, status_info, status_warn,
    utils::{
        f128, symbolica_ext::LOGPRINTOPTS, tracing::StatusRenderable, FloatLike,
        PrecisionUpgradable, F, GS, TENSORLIB,
    },
    GammaLoopContext,
};

/// Check if debug cache mode is enabled
#[inline]
fn is_debug_cache_enabled() -> bool {
    std::env::var("GAMMALOOP_DEBUG_CACHE")
        .map(|v| v == "1" || v.to_lowercase() == "true")
        .unwrap_or(false)
        || cfg!(debug_assertions)
}

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ParamValuePairs {
    pub value_range: Range<usize>,
    pub params: Vec<Atom>,
}

impl ParamValuePairs {
    pub fn validate(&self) {
        assert_eq!(
            self.value_range.len(),
            self.params.len(),
            "Number of values and parameters must match"
        );
    }

    pub fn default_from_symbol(param: Symbol) -> Self {
        Self {
            value_range: 0..1,
            params: vec![Atom::var(param)],
        }
    }

    // pub fn replacement(&self) -> Vec<Replacement> {
    //     let mut replacements = Vec::new();
    //     for (p, v) in self.params.iter().zip_eq(self.values.iter()) {
    //         println!("{p}");
    //         let replacement = Replacement::new(
    //             p.clone().to_pattern(),
    //             Atom::num(v.clone().to_coefficient()),
    //         );
    //         replacements.push(replacement);
    //     }
    //     replacements
    // }
}

impl<'a, A: Into<AtomOrView<'a>>> FromIterator<A> for ParamValuePairs {
    fn from_iter<T: IntoIterator<Item = A>>(iter: T) -> Self {
        let params: Vec<Atom> = iter.into_iter().map(|a| a.into().into_owned()).collect();
        let len = params.len();
        Self {
            value_range: 0..len,
            params,
        }
    }
}

impl Default for ParamValuePairs {
    fn default() -> Self {
        Self {
            value_range: Range::default(),
            params: Vec::new(),
        }
    }
}

pub trait SplitPolarizations {
    fn polarizations(&self) -> Vec<Atom>;
}

pub trait ParamBuilderGraph {
    fn get_external_energy_atoms(&self) -> Vec<Atom>;
    fn iter_edge_ids(&self) -> impl Iterator<Item = EdgeIndex> + '_;
    fn external_spatial_params(&self) -> Vec<Atom>;
    fn emr_spatial_params(&self) -> Vec<Atom>;
    fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom;
    fn get_ose_replacements(&self) -> Vec<Replacement>;
}

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GammaLoopPairs {
    m_uv: ParamValuePairs,
    idenso_vars: ParamValuePairs,
    mu_r_sq: ParamValuePairs,
    orientations: ParamValuePairs,
    pub model_parameters: ParamValuePairs,
    external_energies: ParamValuePairs,
    external_spatial: ParamValuePairs,
    pub polarizations: ParamValuePairs,
    emr_spatial: ParamValuePairs,
    tstar: ParamValuePairs,
    h_function: ParamValuePairs,
    esurface_derivative: ParamValuePairs,
    uv_damp_plus: ParamValuePairs,
    uv_damp_minus: ParamValuePairs,
    radius: ParamValuePairs,
    radius_star: ParamValuePairs,
}

impl IntoIterator for GammaLoopPairs {
    type Item = ParamValuePairs;
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

impl<'a> IntoIterator for &'a GammaLoopPairs {
    type Item = &'a ParamValuePairs;
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

impl<'a> IntoIterator for &'a mut GammaLoopPairs {
    type Item = &'a mut ParamValuePairs;
    type IntoIter = std::array::IntoIter<Self::Item, 16>;

    fn into_iter(self) -> Self::IntoIter {
        [
            &mut self.m_uv,
            &mut self.mu_r_sq,
            &mut self.idenso_vars,
            &mut self.model_parameters,
            &mut self.external_energies,
            &mut self.external_spatial,
            &mut self.polarizations,
            &mut self.emr_spatial,
            &mut self.tstar,
            &mut self.h_function,
            &mut self.esurface_derivative,
            &mut self.uv_damp_plus,
            &mut self.uv_damp_minus,
            &mut self.radius,
            &mut self.radius_star,
            &mut self.orientations,
        ]
        .into_iter()
    }
}

impl GammaLoopPairs {
    pub fn update_ranges(&mut self) -> usize {
        let mut start = 0;
        for pair in self {
            let len = pair.params.len();
            pair.value_range = start..start + len;
            start += len;
        }
        start
    }

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

    pub(crate) fn new<G: ParamBuilderGraph + SplitPolarizations>(
        model: &Model,
        graph: &G,
    ) -> (Self, usize) {
        let mut pairs: GammaLoopPairs = Default::default();

        pairs.m_uv = ParamValuePairs::default_from_symbol(GS.m_uv);

        pairs.mu_r_sq = ParamValuePairs::default_from_symbol(GS.mu_r_sq);
        pairs.tstar = ParamValuePairs::default_from_symbol(GS.rescale_star);
        pairs.radius = ParamValuePairs::default_from_symbol(GS.radius);
        pairs.radius_star = ParamValuePairs::default_from_symbol(GS.radius_star);
        pairs.h_function = ParamValuePairs::default_from_symbol(GS.hfunction);
        pairs.esurface_derivative = ParamValuePairs::default_from_symbol(GS.deta);
        pairs.uv_damp_plus = ParamValuePairs::default_from_symbol(GS.uv_damp_plus);
        pairs.uv_damp_minus = ParamValuePairs::default_from_symbol(GS.uv_damp_minus);

        pairs.idenso_vars();

        pairs.update_model(model);
        pairs.external_energies(graph);
        pairs.orientations(graph);
        pairs.polarizations(graph);
        pairs.external_spatial(graph);
        pairs.emr_spatial(graph);

        let len = pairs.update_ranges();
        (pairs, len)
    }

    pub(crate) fn idenso_vars(&mut self) {
        self.idenso_vars = [CS.tr, CS.nc].into_iter().collect();
        self.update_ranges();
    }

    pub(crate) fn update_model(&mut self, model: &Model) {
        self.model_parameters = model.generate_params().into_iter().collect();
    }

    pub(crate) fn external_spatial<G: ParamBuilderGraph>(&mut self, graph: &G) {
        self.external_spatial.params = graph.external_spatial_params();
    }

    pub(crate) fn emr_spatial<G: ParamBuilderGraph>(&mut self, graph: &G) {
        self.emr_spatial.params = graph.emr_spatial_params();
    }

    pub(crate) fn polarizations<G: ParamBuilderGraph + SplitPolarizations>(&mut self, graph: &G) {
        let mut params = Vec::new();

        for a in graph.polarizations() {
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

        self.polarizations = params.into_iter().collect();
    }

    pub(crate) fn orientations<G: ParamBuilderGraph>(&mut self, graph: &G) {
        let mut params = Vec::new();

        for i in graph.iter_edge_ids() {
            params.push(GS.sign(i));
            params.push(GS.sign_theta(GS.sign(i)));
            params.push(GS.sign_theta(-GS.sign(i)));
        }

        self.orientations.params = params;
    }

    pub(crate) fn external_energies<G: ParamBuilderGraph>(&mut self, graph: &G) {
        self.external_energies = graph.get_external_energy_atoms().into_iter().collect();
    }

    pub(crate) fn add_external_four_mom_impl<T: FloatLike>(
        &self,
        ext: &ExternalFourMomenta<F<T>>,
        values: &mut [Complex<F<T>>],
    ) {
        let mut e_start = self.external_energies.value_range.start;
        let mut s_start = self.external_spatial.value_range.start;

        for e in ext {
            values[e_start] = Complex::new_re(e.temporal.value.clone());
            e_start += 1;
            for c in &e.spatial {
                values[s_start] = Complex::new_re(c.clone());
                s_start += 1;
            }
        }

        debug_assert_eq!(e_start, self.external_energies.value_range.end);
        debug_assert_eq!(s_start, self.external_spatial.value_range.end);
    }

    pub(crate) fn polarizations_values<T: FloatLike>(
        &self,
        graph: &Graph,
        ext: &ExternalFourMomenta<F<T>>,
        helicities: &[Helicity],
        // values: &mut [Complex<F<T>>],
    ) -> Vec<Complex<F<T>>> {
        // let p_start = self.polarizations.value_range.start;
        // let mut p_shift = 0;
        let mut values = Vec::with_capacity(self.polarizations.value_range.len());

        for (p, _) in &graph.polarizations {
            let extid = graph.loop_momentum_basis.ext_from(p.eid).unwrap();
            let hel = p.hel.unwrap_or(helicities[extid.0]);
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

            for val in pol.tensor.data.into_iter() {
                // info!("{}:{}", self.polarizations.params[p_shift], val);
                values.push(val);
            }
        }
        values
    }

    pub(crate) fn threshold_params<T: FloatLike>(
        &self,
        threshold_params: &ThresholdParams<T>,
        values: &mut [Complex<F<T>>],
    ) {
        values[self.radius.value_range.start] = Complex::new_re(threshold_params.radius.clone());
        values[self.radius_star.value_range.start] =
            Complex::new_re(threshold_params.radius_star.clone());
        values[self.esurface_derivative.value_range.start] =
            Complex::new_re(threshold_params.esurface_derivative.clone());

        values[self.uv_damp_plus.value_range.start] =
            Complex::new_re(threshold_params.uv_damp_plus.clone());
        values[self.uv_damp_minus.value_range.start] =
            Complex::new_re(threshold_params.uv_damp_minus.clone());
        values[self.h_function.value_range.start] =
            Complex::new_re(threshold_params.h_function.clone());
    }
}

impl Default for GammaLoopPairs {
    fn default() -> Self {
        Self {
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
        }
    }
}

#[derive(Clone)]
pub struct ParamCache<T: FloatLike> {
    cache: ConstGenericRingBuffer<Vec<Complex<F<T>>>, 64>,
    len: usize,
}

impl<T: FloatLike> ParamCache<T> {
    pub fn get(&self, key: usize) -> Option<&Vec<Complex<F<T>>>> {
        if self.len <= key || key + 64 < self.len {
            None
        } else if self.len < 64 {
            Some(&self.cache[key])
        } else {
            let cachekey = (64 + key) - self.len;
            Some(&self.cache[cachekey])
        }
    }

    pub fn checked_push(&mut self, key: usize, value: Vec<Complex<F<T>>>) {
        debug_assert_eq!(self.len, key);
        self.len += 1;
        self.cache.enqueue(value);
    }
}

impl<T: FloatLike> Default for ParamCache<T> {
    fn default() -> Self {
        Self {
            cache: ConstGenericRingBuffer::new(),
            len: 0,
        }
    }
}

impl<T: FloatLike + Encode> Encode for ParamCache<T> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        _encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        Ok(())
    }
}

impl<C, T: FloatLike + Decode<C>> Decode<C> for ParamCache<T> {
    fn decode<D: bincode::de::Decoder>(
        _decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        Ok(Self::default())
    }
}

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ParamBuilder<T: FloatLike = f64> {
    pub values: Vec<Complex<F<T>>>,
    pub pairs: GammaLoopPairs,
    pub polarization_cache: ParamCache<T>,

    pub reps: Vec<(Atom, Atom)>,
    // pub eager_const_map: HashMap<Atom, Complex<F<T>>>,
    // pub eager_function_map: HashMap<Symbol, EvaluationFn<Atom, Complex<F<T>>>>,
    // pub eager_fn_map:
    pub fn_map: FunctionMap,
}

impl<T: FloatLike> Default for ParamBuilder<T> {
    fn default() -> Self {
        Self::new_empty()
    }
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
        cache: bool,
        sample: &'a MomentumSample<T>,
        graph: &'a Graph,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<T>>,
    ) -> Cow<'a, Vec<Complex<F<T>>>>;
}

impl UpdateAndGetParams<f64> for ParamBuilder<f64> {
    fn update_emr_and_get_params<'a>(
        &'a mut self,
        cache: bool,
        sample: &'a MomentumSample<f64>,
        graph: &'a Graph,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<f64>>,
    ) -> Cow<'a, Vec<Complex<F<f64>>>> {
        let emr_start = self.pairs.emr_spatial.value_range.start;
        let mut shift = 0;
        graph.iter_edges().for_each(|(pair, edge_id, _)| {
            if let HedgePair::Paired { .. } = pair {
                let emr_vec = graph.loop_momentum_basis.edge_signatures[edge_id]
                    .compute_three_momentum_from_four(sample.loop_moms(), sample.external_moms());
                // info!("px:{}", self.pairs.emr_spatial.params[shift]);
                //
                self.values[emr_start + shift] = Complex::new_re(emr_vec.px);
                shift += 1;
                // info!("py:{}", self.pairs.emr_spatial.params[shift]);
                self.values[emr_start + shift] = Complex::new_re(emr_vec.py);
                shift += 1;
                // info!("pz:{}", self.pairs.emr_spatial.params[shift]);
                self.values[emr_start + shift] = Complex::new_re(emr_vec.pz);
                shift += 1;
            }
        });

        // parse!("s").evaluator(fn_map, params, optimization_settings).unwrap().

        self.add_external_four_mom(sample.external_moms());
        self.polarizations_values(cache, graph, sample, helicities);

        if let Some(threshold_params) = threshold_params {
            self.threshold_params(threshold_params);
        }

        Cow::Borrowed(&self.values)
    }
}

impl UpdateAndGetParams<f128> for ParamBuilder<f64> {
    fn update_emr_and_get_params(
        &mut self,
        _cache: bool,
        sample: &MomentumSample<f128>,
        graph: &Graph,
        helicities: &[Helicity],
        threshold_params: Option<&ThresholdParams<f128>>,
    ) -> Cow<Vec<Complex<F<f128>>>> {
        let mut emr_start = self.pairs.emr_spatial.value_range.start;
        let mut values = self.higher();
        graph.iter_edges().for_each(|(pair, edge_id, _)| {
            if let HedgePair::Paired { .. } = pair {
                let emr_vec = graph.loop_momentum_basis.edge_signatures[edge_id]
                    .compute_three_momentum_from_four(sample.loop_moms(), sample.external_moms());
                // println!("{edge_id}:{emr_vec}");
                values[emr_start] = Complex::new_re(emr_vec.px);
                emr_start += 1;
                values[emr_start] = Complex::new_re(emr_vec.py);
                emr_start += 1;
                values[emr_start] = Complex::new_re(emr_vec.pz);
                emr_start += 1;
            }
        });

        self.pairs
            .add_external_four_mom_impl(sample.external_moms(), &mut values);
        values[self.pairs.polarizations.value_range.clone()].clone_from_slice(
            &self
                .pairs
                .polarizations_values(graph, sample.external_moms(), helicities),
        );

        if let Some(threshold_params) = threshold_params {
            self.pairs.threshold_params(threshold_params, &mut values);
        }

        Cow::Owned(values)
    }
}

impl<T: FloatLike> ParamBuilder<T>
where
    T::Higher: FloatLike,
    T::Lower: FloatLike,
{
    fn higher(&self) -> Vec<Complex<F<T::Higher>>> {
        self.values.iter().map(|v| v.higher()).collect()
    }
    // fn lower(&self) -> Vec<Complex<F<T::Lower>>> {
    //     self.values.iter().map(|v| v.lower()).collect()
    // }
}
impl<T: FloatLike> ParamBuilder<T> {
    pub fn model_values(&self) -> &[Complex<F<T>>] {
        let range = self.pairs.model_parameters.value_range.clone();
        &self.values[range]
    }
    pub fn validate(&self) {
        self.pairs.validate();
    }

    pub fn add_tagged_function(
        &mut self,
        name: Symbol,
        tags: Vec<Atom>,
        rename: String,
        args: Vec<Symbol>,
        body: Atom,
    ) -> Result<(), String> {
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
    ) -> Result<(), String> {
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

    pub(crate) fn new_empty() -> Self {
        Self {
            polarization_cache: ParamCache::default(),

            fn_map: FunctionMap::default(),
            pairs: GammaLoopPairs::default(),
            values: Vec::new(),

            reps: Vec::new(),
        }
    }

    pub(crate) fn new<G: ParamBuilderGraph + SplitPolarizations>(graph: &G, model: &Model) -> Self {
        let (pairs, len) = GammaLoopPairs::new(model, graph);

        let mut new = Self {
            pairs,
            ..Default::default()
        };

        let pi_rational = Rational::from(std::f64::consts::PI);

        for e in graph.iter_edge_ids() {
            new.add_tagged_function(
                GS.ose,
                vec![Atom::num(e.0 as i64)],
                format!("OSE{e}"),
                vec![],
                graph.explicit_ose_atom(e),
            )
            .unwrap();
        }
        new.add_constant(GS.pi.into(), pi_rational.into());

        new.values = vec![Complex::new_re(F(T::from_f64(0.))); len];
        new.update_model_values(model);
        new.update_idenso_values();
        new
    }

    #[inline]
    pub(crate) fn m_uv_value(&mut self, m_uv: Complex<F<T>>) {
        debug_assert!(self.pairs.m_uv.value_range.len() == 1);
        self.values[self.pairs.m_uv.value_range.start] = m_uv;
    }

    pub(crate) fn mu_r_sq_value(&mut self, mu_r_sq: Complex<F<T>>) {
        debug_assert!(self.pairs.mu_r_sq.value_range.len() == 1);
        self.values[self.pairs.mu_r_sq.value_range.start] = mu_r_sq;
    }

    pub(crate) fn update_idenso_values(&mut self) {
        let tr_value = Complex::new_re(F(T::from_f64(0.5)));
        let nc_value = Complex::new_re(F(T::from_f64(3.)));

        debug_assert!(self.pairs.idenso_vars.value_range.len() == 2);
        self.values[self.pairs.idenso_vars.value_range.start] = tr_value;
        self.values[self.pairs.idenso_vars.value_range.start + 1] = nc_value;
    }

    pub(crate) fn update_model_values(&mut self, model: &Model) {
        let mut pos = self.pairs.model_parameters.value_range.start;
        for cpl in model.couplings.values().filter(|c| c.value.is_some()) {
            if let Some(value) = cpl.value {
                self.values[pos] = value.map(F::from_f64);
                pos += 1;
            }
        }
        for param in model.parameters.values().filter(|p| p.value.is_some()) {
            if let Some(value) = param.value {
                let value = Complex::new(F::<T>::from_ff64(value.re), F::<T>::from_ff64(value.im));
                self.values[pos] = value;
                pos += 1;
            }
        }

        debug_assert_eq!(pos, self.pairs.model_parameters.value_range.end);
    }

    pub(crate) fn add_external_four_mom(&mut self, ext: &ExternalFourMomenta<F<T>>) {
        self.pairs.add_external_four_mom_impl(ext, &mut self.values);
    }

    pub(crate) fn orientation_value<O: GraphOrientation>(&mut self, orientation: &O) {
        let zero: Complex<F<T>> = Complex::new_re(F(T::from_f64(0.)));
        let one = zero.ref_one();
        let minusone = -(one.clone());

        let mut o_start = self.pairs.orientations.value_range.start;

        for (_, i) in orientation.orientation() {
            match i {
                Orientation::Default => {
                    self.values[o_start] = one.clone();
                    o_start += 1;
                    self.values[o_start] = one.clone();
                    o_start += 1;
                    self.values[o_start] = zero.clone();
                    o_start += 1;
                }
                Orientation::Reversed => {
                    self.values[o_start] = minusone.clone();
                    o_start += 1;
                    self.values[o_start] = zero.clone();
                    o_start += 1;
                    self.values[o_start] = one.clone();
                    o_start += 1;
                }
                Orientation::Undirected => {
                    self.values[o_start] = zero.clone();
                    o_start += 1;
                    self.values[o_start] = zero.clone();
                    o_start += 1;
                    self.values[o_start] = zero.clone();
                    o_start += 1;
                }
            }
        }
    }

    pub(crate) fn polarizations_values(
        &mut self,
        cache: bool,
        graph: &Graph,
        sample: &MomentumSample<T>,
        helicities: &[Helicity],
    ) {
        let cache_id = sample.sample.external_mom_cache_id;
        let base_cache_id = sample.sample.external_mom_base_cache_id;

        // Validate cache consistency
        if cache_id < base_cache_id {
            status_warn!(
                "WARNING: Cache ID inconsistency detected! current={}, base={}",
                cache_id,
                base_cache_id
            );
        }

        // Compute expected polarizations for debug search
        let expected_pols = if is_debug_cache_enabled() {
            Some(
                self.pairs
                    .polarizations_values(graph, sample.external_moms(), helicities),
            )
        } else {
            None
        };

        // Try to get from cache first
        let pols = if let Some(v) = self.polarization_cache.get(cache_id) {
            status_debug!(
                "Cache HIT for external_cache_id={} (base={})",
                cache_id,
                base_cache_id
            );

            // Validate cached values in debug mode
            if let Some(ref expected) = expected_pols {
                debug_assert_eq!(
                    v, expected,
                    "Cached polarizations don't match computed values for cache_id={}",
                    cache_id
                );
            }

            if !cache {
                self.values[self.pairs.polarizations.value_range.clone()].clone_from_slice(v);
                return;
            }
            v
        } else {
            status_debug!(
                "Cache MISS for external_cache_id={} (base={})",
                cache_id,
                base_cache_id
            );

            // DEBUG: Search entire cache for matching polarizations to detect missed hits
            if is_debug_cache_enabled() && cache {
                if let Some(ref expected) = expected_pols {
                    self.debug_search_cache_for_matching_polarizations(
                        expected,
                        cache_id,
                        base_cache_id,
                        sample.external_moms(),
                    );
                }
            }

            if !cache {
                let computed_pols = expected_pols.unwrap_or_else(|| {
                    self.pairs
                        .polarizations_values(graph, sample.external_moms(), helicities)
                });
                self.values[self.pairs.polarizations.value_range.clone()]
                    .clone_from_slice(&computed_pols);
                return;
            }

            // Compute and cache new polarizations
            let computed_pols = expected_pols.unwrap_or_else(|| {
                self.pairs
                    .polarizations_values(graph, sample.external_moms(), helicities)
            });

            self.polarization_cache
                .checked_push(cache_id, computed_pols);
            status_debug!("Cached new polarizations for cache_id={}", cache_id);

            self.polarization_cache.get(cache_id).unwrap()
        };

        self.values[self.pairs.polarizations.value_range.clone()].clone_from_slice(pols);
    }

    /// Debug function to search cache for matching polarizations to detect missed cache hits
    fn debug_search_cache_for_matching_polarizations(
        &self,
        expected_pols: &[Complex<F<T>>],
        current_cache_id: usize,
        base_cache_id: usize,
        external_moms: &ExternalFourMomenta<F<T>>,
    ) {
        // Only run this expensive search if debug mode is enabled
        if !is_debug_cache_enabled() {
            return;
        }

        let tolerance = F::<T>::from_f64(1e-12);
        let mut found_matches = Vec::new();

        // Search through all possible cache IDs that could contain our polarizations
        let max_search_id = std::cmp::max(current_cache_id + 10, self.polarization_cache.len);

        for search_id in 0..max_search_id {
            if let Some(cached_pols) = self.polarization_cache.get(search_id) {
                if cached_pols.len() == expected_pols.len() {
                    let matches =
                        cached_pols
                            .iter()
                            .zip(expected_pols.iter())
                            .all(|(cached, expected)| {
                                (&cached.re - &expected.re).abs() < tolerance
                                    && (&cached.im - &expected.im).abs() < tolerance
                            });

                    if matches {
                        found_matches.push(search_id);
                    }
                }
            }
        }

        if !found_matches.is_empty() {
            status_warn!(
                "🔍 DEBUG CACHE SEARCH: Found {} matching polarization(s) in cache but using cache_id={}!",
                found_matches.len(),
                current_cache_id
            );
            status_warn!(
                "   ⚠️  MISSED CACHE HITS: Found identical polarizations at cache_id(s): {:?}",
                found_matches
            );
            status_warn!(
                "   📊 Current: cache_id={}, base_cache_id={}",
                current_cache_id,
                base_cache_id
            );

            // Provide actionable debugging information
            if found_matches.contains(&base_cache_id) {
                status_warn!(
                    "   💡 HINT: Base cache_id {} has matching polarizations. You should call revert_to_base_external_cache_id()!",
                    base_cache_id
                );
            }

            // Check if any found matches are close to base_cache_id (indicating related configurations)
            let related_ids: Vec<_> = found_matches
                .iter()
                .filter(|&&id| id >= base_cache_id && id <= base_cache_id + 20)
                .collect();

            if !related_ids.is_empty() {
                status_warn!(
                    "   🔗 RELATED IDs: Cache IDs {:?} are related to base_id {} and contain identical polarizations",
                    related_ids,
                    base_cache_id
                );
            }

            // Log external momentum info for debugging
            status_debug!(
                "   🔍 External momenta: [{:.6}, {:.6}, {:.6}, {:.6}] (first external)",
                external_moms
                    .first()
                    .map_or(&F::from_f64(0.0), |m| &m.temporal.value),
                external_moms
                    .first()
                    .map_or(&F::from_f64(0.0), |m| &m.spatial.px),
                external_moms
                    .first()
                    .map_or(&F::from_f64(0.0), |m| &m.spatial.py),
                external_moms
                    .first()
                    .map_or(&F::from_f64(0.0), |m| &m.spatial.pz),
            );
        } else {
            // Only log this in very verbose debug mode
            if std::env::var("GAMMALOOP_DEBUG_CACHE_VERBOSE").is_ok() {
                status_debug!(
                    "🔍 DEBUG CACHE SEARCH: No matching polarizations found in cache (searched {} entries). This is a genuine cache miss.",
                    max_search_id
                );
            }
        }
    }

    /// Public method to enable/disable debug cache mode at runtime
    pub fn set_debug_cache_mode(enabled: bool) {
        std::env::set_var("GAMMALOOP_DEBUG_CACHE", if enabled { "1" } else { "0" });
        status_info!(
            "Debug cache mode {}: Set GAMMALOOP_DEBUG_CACHE={}",
            if enabled { "enabled" } else { "disabled" },
            if enabled { "1" } else { "0" }
        );
    }

    /// Check current debug cache mode status
    pub fn get_debug_cache_status() -> DebugCacheStatus {
        let env_enabled = std::env::var("GAMMALOOP_DEBUG_CACHE")
            .map(|v| v == "1" || v.to_lowercase() == "true")
            .unwrap_or(false);
        let verbose_enabled = std::env::var("GAMMALOOP_DEBUG_CACHE_VERBOSE").is_ok();
        let debug_assertions = cfg!(debug_assertions);

        DebugCacheStatus {
            environment_enabled: env_enabled,
            debug_assertions,
            verbose_enabled,
            effective_enabled: env_enabled || debug_assertions,
        }
    }

    /// Debug function that's available in release mode for critical cache analysis
    #[allow(dead_code)]
    fn analyze_cache_efficiency(&self, sample: &MomentumSample<T>) -> CacheAnalysisResult<T> {
        let mut analysis = CacheAnalysisResult::<T> {
            _phantom: std::marker::PhantomData,
            total_entries: self.polarization_cache.len,
            current_cache_id: sample.sample.external_mom_cache_id,
            base_cache_id: sample.sample.external_mom_base_cache_id,
            entries_since_base: sample
                .sample
                .external_mom_cache_id
                .saturating_sub(sample.sample.external_mom_base_cache_id),
            potential_reuses: 0,
            cache_span: 0,
        };

        // Calculate cache span (range of used cache IDs)
        if self.polarization_cache.len > 0 {
            analysis.cache_span = std::cmp::min(64, self.polarization_cache.len);
        }

        // Estimate potential reuses based on pattern
        if analysis.entries_since_base > 0 {
            // If we're far from base, there might be reuse opportunities
            analysis.potential_reuses = analysis.entries_since_base.saturating_sub(1);
        }

        analysis
    }

    pub(crate) fn threshold_params(&mut self, threshold_params: &ThresholdParams<T>) {
        self.pairs
            .threshold_params(threshold_params, &mut self.values);
    }

    pub fn table(&self) -> Table {
        let mut table = tabled::builder::Builder::new();

        for (lhs, rhs) in &self.reps {
            table.push_record(vec![
                lhs.printer(LOGPRINTOPTS).to_string(),
                rhs.printer(LOGPRINTOPTS).to_string(),
            ]);
        }

        for i in &self.pairs {
            for (p, v) in i.params.iter().zip(i.value_range.clone().into_iter()) {
                table.push_record(vec![
                    p.printer(LOGPRINTOPTS).to_string().to_string(),
                    self.values[v].to_string(),
                    format!("{}", v),
                ]);
            }
        }

        table.build()
    }
}

impl<T: FloatLike> StatusRenderable for ParamBuilder<T> {
    fn status_json(&self) -> serde_json::Value {
        let mut map = serde_json::Map::new();

        for (lhs, rhs) in &self.reps {
            map.insert(lhs.to_canonical_string(), rhs.to_canonical_string().into());
        }

        for pairs in &self.pairs {
            for (param, value) in pairs
                .params
                .iter()
                .zip(pairs.value_range.clone().into_iter())
            {
                map.insert(
                    param.to_canonical_string(),
                    serde_json::to_value(&self.values[value]).unwrap(),
                );
            }
        }

        serde_json::Value::Object(map)
    }
    fn status_pretty(&self) -> Cow<'_, str> {
        format!("\n{}", self).into()
    }
}

impl<T: FloatLike> Display for ParamBuilder<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.table().with(Style::rounded()).to_string().fmt(f)
    }
}

/// Debug analysis result for cache efficiency
#[derive(Debug, Clone)]
pub struct CacheAnalysisResult<T: FloatLike> {
    _phantom: std::marker::PhantomData<T>,
    pub total_entries: usize,
    pub current_cache_id: usize,
    pub base_cache_id: usize,
    pub entries_since_base: usize,
    pub potential_reuses: usize,
    pub cache_span: usize,
}

impl<T: FloatLike> std::fmt::Display for CacheAnalysisResult<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Cache Analysis: {} entries, current={}, base={}, span={}, potential_reuses={}",
            self.total_entries,
            self.current_cache_id,
            self.base_cache_id,
            self.cache_span,
            self.potential_reuses
        )
    }
}

/// Status of debug cache mode
#[derive(Debug, Clone)]
pub struct DebugCacheStatus {
    pub environment_enabled: bool,
    pub debug_assertions: bool,
    pub verbose_enabled: bool,
    pub effective_enabled: bool,
}

impl std::fmt::Display for DebugCacheStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Debug Cache Status: {} (env={}, debug_assertions={}, verbose={})",
            if self.effective_enabled {
                "ENABLED"
            } else {
                "DISABLED"
            },
            self.environment_enabled,
            self.debug_assertions,
            self.verbose_enabled
        )
    }
}

#[test]
fn evaltest() {
    use symbolica::evaluate::{FunctionMap, OptimizationSettings};
    use symbolica::{atom::AtomCore, parse, symbol};
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
    let _ = fn_map.add_function(symbol!("delta_sigma"), "delta_sigma".into(), args, body);
    let optimization_settings = OptimizationSettings::default();
    let mut evaluator = expr
        .evaluator(&fn_map, &params, optimization_settings)
        .unwrap()
        .map_coeff(&|x| x.to_real().unwrap().to_f64());
    assert_eq!(evaluator.evaluate_single(&[1.0, 1.0, 4.0, 4.0]), 8.0);
}
