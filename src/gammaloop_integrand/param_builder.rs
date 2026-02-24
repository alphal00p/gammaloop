use std::{
    borrow::Cow,
    fmt::Display,
    ops::{Deref, Range},
};

use bincode::{Decode, Encode};
use idenso::color::CS;

use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};
use ringbuffer::{ConstGenericRingBuffer, RingBuffer};
use spenso::{
    algebra::{
        algebraic_traits::{RefOne, RefZero},
        complex::Complex,
    },
    iterators::IteratableTensor,
    network::{ExecutionResult, parsing::ParseSettings},
    structure::concrete_index::ExpandedIndex,
    tensors::parametric::AtomViewOrConcrete,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, FunctionBuilder, Indeterminate, Symbol},
    domains::rational::Rational,
    evaluate::FunctionMap,
    id::Replacement,
    parse_lit, symbol,
};
use tabled::{Table, settings::Style};
use tracing::debug;
use tracing::warn;

use crate::{
    GammaLoopContext,
    cff::expression::GraphOrientation,
    gammaloop_integrand::{
        amplitude::export::atom_to_bytes_for_mode,
        evaluators::{InputParams, SliceMut},
    },
    graph::{Graph, LoopMomentumBasis},
    model::Model,
    momentum::{Helicity, PolType},
    momentum_sample::{ExternalFourMomenta, MomentumSample},
    numerator::ParsingNet,
    processes::StandaloneExportSettings,
    utils::{
        ArbPrec, F, FloatLike, GS, PrecisionUpgradable, TENSORLIB, f128,
        hyperdual_utils::DualOrNot, symbolica_ext::LOGPRINTOPTS, tracing::StatusRenderable,
    },
};

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
#[derive(Default)]
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

pub trait SplitPolarizations {
    fn polarizations(&self) -> Vec<Atom>;
}

pub trait ParamBuilderGraph {
    fn get_external_energy_atoms(&self) -> Vec<Atom>;
    fn iter_edge_ids(&self) -> impl Iterator<Item = EdgeIndex> + '_;
    fn external_spatial_params(&self) -> Vec<Atom>;
    fn loop_mom_params(&self, lmb: &LoopMomentumBasis) -> Vec<Atom>;
    fn explicit_ose_atom(&self, edge: EdgeIndex) -> Atom;
    fn get_ose_replacements(&self) -> Vec<Replacement>;
}

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
#[derive(Default)]
pub struct GammaLoopPairs {
    m_uv: ParamValuePairs,
    idenso_vars: ParamValuePairs,
    mu_r_sq: ParamValuePairs,
    orientations: ParamValuePairs,
    override_if: ParamValuePairs,
    pub model_parameters: ParamValuePairs,
    external_energies: ParamValuePairs,
    external_spatial: ParamValuePairs,
    pub polarizations: ParamValuePairs,
    loop_moms_spatial: ParamValuePairs,
    tstar: ParamValuePairs,
    h_function_lu_cut: ParamValuePairs,
    h_function_left_th: ParamValuePairs,
    h_function_right_th: ParamValuePairs,
    esurface_derivative_lu_cut: ParamValuePairs,
    esurface_derivative_left_th: ParamValuePairs,
    esurface_derivative_right_th: ParamValuePairs,
    uv_damp_plus_left: ParamValuePairs,
    uv_damp_minus_left: ParamValuePairs,
    radius_left: ParamValuePairs,
    radius_star_left: ParamValuePairs,
    uv_damp_plus_right: ParamValuePairs,
    uv_damp_minus_right: ParamValuePairs,
    radius_right: ParamValuePairs,
    radius_star_right: ParamValuePairs,
    pub additional_params: ParamValuePairs,
}

impl IntoIterator for GammaLoopPairs {
    type Item = ParamValuePairs;
    type IntoIter = std::array::IntoIter<Self::Item, 26>;

    fn into_iter(self) -> Self::IntoIter {
        [
            self.m_uv,
            self.mu_r_sq,
            self.idenso_vars,
            self.model_parameters,
            self.external_energies,
            self.external_spatial,
            self.polarizations,
            self.loop_moms_spatial,
            self.tstar,
            self.h_function_lu_cut,
            self.h_function_left_th,
            self.h_function_right_th,
            self.esurface_derivative_lu_cut,
            self.esurface_derivative_left_th,
            self.esurface_derivative_right_th,
            self.uv_damp_plus_left,
            self.uv_damp_minus_left,
            self.radius_left,
            self.radius_star_left,
            self.uv_damp_plus_right,
            self.uv_damp_minus_right,
            self.radius_right,
            self.radius_star_right,
            self.orientations,
            self.override_if,
            self.additional_params,
        ]
        .into_iter()
    }
}

impl<'a> IntoIterator for &'a GammaLoopPairs {
    type Item = &'a ParamValuePairs;
    type IntoIter = std::array::IntoIter<Self::Item, 26>;

    fn into_iter(self) -> Self::IntoIter {
        [
            &self.m_uv,
            &self.mu_r_sq,
            &self.idenso_vars,
            &self.model_parameters,
            &self.external_energies,
            &self.external_spatial,
            &self.polarizations,
            &self.loop_moms_spatial,
            &self.tstar,
            &self.h_function_lu_cut,
            &self.h_function_left_th,
            &self.h_function_right_th,
            &self.esurface_derivative_lu_cut,
            &self.esurface_derivative_left_th,
            &self.esurface_derivative_right_th,
            &self.uv_damp_plus_left,
            &self.uv_damp_minus_left,
            &self.radius_left,
            &self.radius_star_left,
            &self.uv_damp_plus_right,
            &self.uv_damp_minus_right,
            &self.radius_right,
            &self.radius_star_right,
            &self.orientations,
            &self.override_if,
            &self.additional_params,
        ]
        .into_iter()
    }
}

impl<'a> IntoIterator for &'a mut GammaLoopPairs {
    type Item = &'a mut ParamValuePairs;
    type IntoIter = std::array::IntoIter<Self::Item, 26>;

    fn into_iter(self) -> Self::IntoIter {
        [
            &mut self.m_uv,
            &mut self.mu_r_sq,
            &mut self.idenso_vars,
            &mut self.model_parameters,
            &mut self.external_energies,
            &mut self.external_spatial,
            &mut self.polarizations,
            &mut self.loop_moms_spatial,
            &mut self.tstar,
            &mut self.h_function_lu_cut,
            &mut self.h_function_left_th,
            &mut self.h_function_right_th,
            &mut self.esurface_derivative_lu_cut,
            &mut self.esurface_derivative_left_th,
            &mut self.esurface_derivative_right_th,
            &mut self.uv_damp_plus_left,
            &mut self.uv_damp_minus_left,
            &mut self.radius_left,
            &mut self.radius_star_left,
            &mut self.uv_damp_plus_right,
            &mut self.uv_damp_minus_right,
            &mut self.radius_right,
            &mut self.radius_star_right,
            &mut self.orientations,
            &mut self.override_if,
            &mut self.additional_params,
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
        debug!("Validating orientations");
        self.orientations.validate();
        self.override_if.validate();
        debug!("Validating emr_spatial");
        self.loop_moms_spatial.validate();
        debug!("Validating tstar");
        self.tstar.validate();
        debug!("Validating h_function");
        self.h_function_lu_cut.validate();
        self.h_function_left_th.validate();
        self.h_function_right_th.validate();
        debug!("Validating derivative_at_tstar");
        self.esurface_derivative_lu_cut.validate();
        self.esurface_derivative_left_th.validate();
        self.esurface_derivative_right_th.validate();
        debug!("Validating uv_damp_plus");
        self.uv_damp_plus_left.validate();
        self.uv_damp_plus_right.validate();
        debug!("Validating uv_damp_minus");
        self.uv_damp_minus_left.validate();
        self.uv_damp_minus_right.validate();
        debug!("Validating radius");
        self.radius_left.validate();
        self.radius_right.validate();
        debug!("Validating radius_star");
        self.radius_star_left.validate();
        self.radius_star_right.validate();
    }

    pub(crate) fn new<
        'a,
        A: Into<AtomOrView<'a>>,
        G: ParamBuilderGraph + SplitPolarizations,
        T: IntoIterator<Item = A>,
    >(
        model: &Model,
        graph: &G,
        lmb: &LoopMomentumBasis,
        additional_params: T,
    ) -> (Self, usize) {
        let mut pairs = GammaLoopPairs {
            m_uv: ParamValuePairs::default_from_symbol(GS.m_uv),
            mu_r_sq: ParamValuePairs::default_from_symbol(GS.mu_r_sq),
            tstar: ParamValuePairs::default_from_symbol(GS.rescale_star),
            radius_left: ParamValuePairs::default_from_symbol(GS.radius_left),
            radius_right: ParamValuePairs::default_from_symbol(GS.radius_right),
            radius_star_left: ParamValuePairs::default_from_symbol(GS.radius_star_left),
            radius_star_right: ParamValuePairs::default_from_symbol(GS.radius_star_right),
            h_function_lu_cut: ParamValuePairs::default_from_symbol(GS.hfunction_lu_cut),
            h_function_left_th: ParamValuePairs::default_from_symbol(GS.hfunction_left_th),
            h_function_right_th: ParamValuePairs::default_from_symbol(GS.hfunction_right_th),
            esurface_derivative_lu_cut: ParamValuePairs::default_from_symbol(GS.deta_lu_cut),
            esurface_derivative_left_th: ParamValuePairs::default_from_symbol(GS.deta_left_th),
            esurface_derivative_right_th: ParamValuePairs::default_from_symbol(GS.deta_right_th),
            uv_damp_plus_left: ParamValuePairs::default_from_symbol(GS.uv_damp_plus_left),
            uv_damp_plus_right: ParamValuePairs::default_from_symbol(GS.uv_damp_plus_right),
            uv_damp_minus_left: ParamValuePairs::default_from_symbol(GS.uv_damp_minus_left),
            uv_damp_minus_right: ParamValuePairs::default_from_symbol(GS.uv_damp_minus_right),
            additional_params: additional_params.into_iter().collect(),
            ..Default::default()
        };

        pairs.idenso_vars();

        pairs.update_model(model);
        pairs.external_energies(graph);
        pairs.orientations(graph);
        pairs.polarizations(graph);
        pairs.external_spatial(graph);
        pairs.loop_moms_spatial(graph, lmb);

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

    pub(crate) fn loop_moms_spatial<G: ParamBuilderGraph>(
        &mut self,
        graph: &G,
        lmb: &LoopMomentumBasis,
    ) {
        self.loop_moms_spatial.params = graph.loop_mom_params(lmb);
    }

    pub(crate) fn polarizations<G: ParamBuilderGraph + SplitPolarizations>(&mut self, graph: &G) {
        let mut params = Vec::new();

        for a in graph.polarizations() {
            match ParsingNet::try_from_view(
                a.as_view(),
                TENSORLIB.read().unwrap().deref(),
                &ParseSettings::default(),
            )
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
        }

        self.override_if.params = vec![Atom::var(GS.override_if)];
        self.orientations.params = params;
    }

    pub(crate) fn external_energies<G: ParamBuilderGraph>(&mut self, graph: &G) {
        self.external_energies = graph.get_external_energy_atoms().into_iter().collect();
    }

    pub(crate) fn add_external_four_mom_impl<T: FloatLike>(
        &self,
        ext: &ExternalFourMomenta<F<T>>,
        values: &mut [Complex<F<T>>],
        multiplicative_offset: usize,
    ) {
        let mut e_start = self.external_energies.value_range.start * multiplicative_offset;
        let mut s_start = self.external_spatial.value_range.start * multiplicative_offset;

        for e in ext {
            values[e_start] = Complex::new_re(e.temporal.value.clone());
            e_start += multiplicative_offset;
            for c in &e.spatial {
                values[s_start] = Complex::new_re(c.clone());
                s_start += multiplicative_offset;
            }
        }

        debug_assert_eq!(
            e_start,
            self.external_energies.value_range.end * multiplicative_offset
        );
        debug_assert_eq!(
            s_start,
            self.external_spatial.value_range.end * multiplicative_offset
        );
    }

    pub(crate) fn add_additional_params<T: FloatLike>(
        &self,
        additional_params: &[F<T>],
        values: &mut [Complex<F<T>>],
        multiplicative_offset: usize,
    ) {
        let mut start = self.additional_params.value_range.start * multiplicative_offset;

        for p in additional_params {
            values[start] = Complex::new_re(p.clone());
            start += multiplicative_offset;
        }

        debug_assert_eq!(
            start,
            self.additional_params.value_range.end * multiplicative_offset,
            "Not filled up additional params"
        );
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

    pub(crate) fn left_threshold_params<T: FloatLike>(
        &self,
        threshold_params: &ThresholdParams<T>,
        values: &mut [Complex<F<T>>],
    ) {
        values[self.radius_left.value_range.start] =
            Complex::new_re(threshold_params.radius.clone());
        values[self.radius_star_left.value_range.start] =
            Complex::new_re(threshold_params.radius_star.clone());
        values[self.esurface_derivative_left_th.value_range.start] =
            Complex::new_re(threshold_params.esurface_derivative.clone());

        values[self.uv_damp_plus_left.value_range.start] =
            Complex::new_re(threshold_params.uv_damp_plus.clone());
        values[self.uv_damp_minus_left.value_range.start] =
            Complex::new_re(threshold_params.uv_damp_minus.clone());
        values[self.h_function_left_th.value_range.start] =
            Complex::new_re(threshold_params.h_function.clone());
    }

    pub(crate) fn right_threshold_params<T: FloatLike>(
        &self,
        threshold_params: &ThresholdParams<T>,
        values: &mut [Complex<F<T>>],
    ) {
        values[self.radius_right.value_range.start] =
            Complex::new_re(threshold_params.radius.clone());
        values[self.radius_star_right.value_range.start] =
            Complex::new_re(threshold_params.radius_star.clone());
        values[self.esurface_derivative_right_th.value_range.start] =
            Complex::new_re(threshold_params.esurface_derivative.clone());

        values[self.uv_damp_plus_right.value_range.start] =
            Complex::new_re(threshold_params.uv_damp_plus.clone());
        values[self.uv_damp_minus_right.value_range.start] =
            Complex::new_re(threshold_params.uv_damp_minus.clone());
        values[self.h_function_right_th.value_range.start] =
            Complex::new_re(threshold_params.h_function.clone());
    }

    pub(crate) fn lu_params<T: FloatLike>(
        &self,
        lu_params: &LUParams<T>,
        values: &mut [Complex<F<T>>],
    ) {
        match (&lu_params.tstar, &lu_params.h_function) {
            (DualOrNot::Dual(tstar), DualOrNot::Dual(h_function)) => {
                let multiplicative_offset = tstar.values.len();

                values[self.tstar.value_range.start * multiplicative_offset..]
                    .iter_mut()
                    .zip(tstar.values.iter())
                    .for_each(|(v, t)| *v = Complex::new_re(t.clone()));

                values[self.h_function_lu_cut.value_range.start * multiplicative_offset..]
                    .iter_mut()
                    .zip(h_function.values.iter())
                    .for_each(|(v, h)| *v = Complex::new_re(h.clone()));
            }
            (DualOrNot::NonDual(tstar), DualOrNot::NonDual(h_function)) => {
                values[self.tstar.value_range.start] = Complex::new_re(tstar.clone());
                values[self.h_function_lu_cut.value_range.start] =
                    Complex::new_re(h_function.clone());
            }
            _ => {
                unreachable!("LU params must both be dual or non-dual");
            }
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
    pub values: Vec<Vec<Complex<F<T>>>>,
    pub pairs: GammaLoopPairs,
    pub polarization_cache: ParamCache<T>,

    pub reps: Vec<FnMapEntry>,
    // pub eager_const_map: HashMap<Atom, Complex<F<T>>>,
    // pub eager_function_map: HashMap<Symbol, EvaluationFn<Atom, Complex<F<T>>>>,
    // pub eager_fn_map:
    pub fn_map: FunctionMap,
}

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode, Debug)]
#[trait_decode(trait = GammaLoopContext)]
pub struct FnMapEntry {
    pub lhs: Atom,
    pub rhs: Atom,
    pub args: Vec<Indeterminate>,
    pub tags: Vec<Atom>,
}

impl FnMapEntry {
    pub fn to_bytes(
        &self,
        settings: &StandaloneExportSettings,
    ) -> Result<(Vec<u8>, Vec<u8>, Vec<Vec<u8>>, Vec<Vec<u8>>)> {
        Ok((
            atom_to_bytes_for_mode(&self.lhs, settings.mode)?,
            atom_to_bytes_for_mode(&self.rhs, settings.mode)?,
            self.tags
                .iter()
                .map(|t| atom_to_bytes_for_mode(t, settings.mode))
                .collect::<Result<Vec<_>>>()?,
            self.args
                .iter()
                .map(|t| atom_to_bytes_for_mode(&Atom::from(t.clone()), settings.mode))
                .collect::<Result<Vec<_>>>()?,
        ))
    }
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

pub struct LUParams<T: FloatLike> {
    pub tstar: DualOrNot<F<T>>,
    pub h_function: DualOrNot<F<T>>,
}

pub trait UpdateAndGetParams<T: FloatLike> {
    #[allow(clippy::too_many_arguments)]
    fn update_emr_and_get_params<'a>(
        &'a mut self,
        cache: (bool, bool),
        sample: &'a MomentumSample<T>,
        graph: &'a Graph,
        helicities: &[Helicity],
        additional_params: &[F<T>],
        left_threshold_params: Option<&ThresholdParams<T>>,
        right_threshold_params: Option<&ThresholdParams<T>>,
        lu_params: Option<&LUParams<T>>,
    ) -> InputParams<'a, T>;
}

impl UpdateAndGetParams<f64> for ParamBuilder<f64> {
    #[allow(clippy::too_many_arguments)]
    fn update_emr_and_get_params<'a>(
        &'a mut self,
        cache: (bool, bool),
        sample: &'a MomentumSample<f64>,
        graph: &'a Graph,
        helicities: &[Helicity],
        additional_params: &[F<f64>],
        left_threshold_params: Option<&ThresholdParams<f64>>,
        right_threshold_params: Option<&ThresholdParams<f64>>,
        lu_params: Option<&LUParams<f64>>,
    ) -> InputParams<'a, f64> {
        let multiplicative_offset = if let Some(dual_loops) = &sample.sample.dual_loop_moms {
            dual_loops.first().unwrap().px.values.len()
        } else {
            1
        };

        let value_index = multiplicative_offset - 1;

        let loop_moms_start =
            self.pairs.loop_moms_spatial.value_range.start * multiplicative_offset;

        let flattened_loop_momenta = if let Some(dual_loop_moms) = &sample.sample.dual_loop_moms {
            dual_loop_moms
                .iter()
                .flat_map(|mom| {
                    vec![
                        mom.px.values.clone(),
                        mom.py.values.clone(),
                        mom.pz.values.clone(),
                    ]
                    .into_iter()
                    .flatten()
                })
                .collect_vec()
        } else {
            sample
                .loop_moms()
                .iter()
                .flat_map(|mom| vec![mom.px.clone(), mom.py.clone(), mom.pz.clone()].into_iter())
                .collect_vec()
        };

        self.values[value_index][loop_moms_start..]
            .iter_mut()
            .zip(flattened_loop_momenta)
            .for_each(|(value, loop_mom_component)| *value = Complex::new_re(loop_mom_component));

        self.pairs.add_additional_params(
            additional_params,
            &mut self.values[value_index],
            multiplicative_offset,
        );

        self.add_external_four_mom(sample.external_moms(), multiplicative_offset);

        // parse!("s").evaluator(fn_map, params, optimization_settings).unwrap().

        self.polarizations_values(cache, graph, sample, helicities, multiplicative_offset);

        if let Some(threshold_params) = left_threshold_params {
            self.left_threshold_params(threshold_params);
        }

        if let Some(threshold_params) = right_threshold_params {
            self.right_threshold_params(threshold_params);
        }

        if let Some(lu_params) = lu_params {
            self.pairs
                .lu_params(lu_params, &mut self.values[value_index]);
        }

        InputParams {
            values: SliceMut::Borrowed(&mut self.values[value_index]),
            multiplicative_offset,
            override_pos: self.pairs.override_if.value_range.start,
            orientations_start: self.pairs.orientations.value_range.start,
        }
    }
}

impl UpdateAndGetParams<f128> for ParamBuilder<f64> {
    #[allow(clippy::too_many_arguments)]
    fn update_emr_and_get_params<'a>(
        &'a mut self,
        _cache: (bool, bool),
        sample: &MomentumSample<f128>,
        graph: &Graph,
        helicities: &[Helicity],
        additional_params: &[F<f128>],
        left_threshold_params: Option<&ThresholdParams<f128>>,
        right_threshold_params: Option<&ThresholdParams<f128>>,
        lu_params: Option<&LUParams<f128>>,
    ) -> InputParams<'a, f128> {
        let multiplicative_offset = if let Some(dual_loops) = &sample.sample.dual_loop_moms {
            dual_loops.first().unwrap().px.values.len()
        } else {
            1
        };

        let value_index = multiplicative_offset - 1;

        let loop_mom_start = self.pairs.loop_moms_spatial.value_range.start * multiplicative_offset;

        let mut values = self.higher(value_index);

        let flattened_loop_momenta = if let Some(dual_loop_moms) = &sample.sample.dual_loop_moms {
            dual_loop_moms
                .iter()
                .flat_map(|mom| {
                    vec![
                        mom.px.values.clone(),
                        mom.py.values.clone(),
                        mom.pz.values.clone(),
                    ]
                    .into_iter()
                    .flatten()
                })
                .collect_vec()
        } else {
            sample
                .loop_moms()
                .clone()
                .into_iter()
                .flat_map(|mom| vec![mom.px, mom.py, mom.pz])
                .collect_vec()
        };

        values[loop_mom_start..]
            .iter_mut()
            .zip(flattened_loop_momenta)
            .for_each(|(value, loop_mom_component)| *value = Complex::new_re(loop_mom_component));

        self.pairs
            .add_additional_params(additional_params, &mut values, multiplicative_offset);

        self.pairs.add_external_four_mom_impl(
            sample.external_moms(),
            &mut values,
            multiplicative_offset,
        );

        let polarization_values =
            self.pairs
                .polarizations_values(graph, sample.external_moms(), helicities);

        let mut polarization_start =
            self.pairs.polarizations.value_range.start * multiplicative_offset;

        for val in polarization_values {
            values[polarization_start] = val;
            polarization_start += multiplicative_offset;
        }

        if let Some(threshold_params) = left_threshold_params {
            self.pairs
                .left_threshold_params(threshold_params, &mut values);
        }

        if let Some(threshold_params) = right_threshold_params {
            self.pairs
                .right_threshold_params(threshold_params, &mut values);
        }

        if let Some(lu_params) = lu_params {
            self.pairs.lu_params(lu_params, &mut values);
        }

        InputParams {
            values: SliceMut::Owned(values),
            multiplicative_offset,
            override_pos: self.pairs.override_if.value_range.start,
            orientations_start: self.pairs.orientations.value_range.start,
        }
    }
}

impl UpdateAndGetParams<ArbPrec> for ParamBuilder<f64> {
    #[allow(clippy::too_many_arguments)]
    fn update_emr_and_get_params<'a>(
        &'a mut self,
        _cache: (bool, bool),
        sample: &MomentumSample<ArbPrec>,
        graph: &Graph,
        helicities: &[Helicity],
        additional_params: &[F<ArbPrec>],
        left_threshold_params: Option<&ThresholdParams<ArbPrec>>,
        right_threshold_params: Option<&ThresholdParams<ArbPrec>>,
        lu_params: Option<&LUParams<ArbPrec>>,
    ) -> InputParams<'a, ArbPrec> {
        let multiplicative_offset = if let Some(dual_loops) = &sample.sample.dual_loop_moms {
            dual_loops.first().unwrap().px.values.len()
        } else {
            1
        };

        let value_index = multiplicative_offset - 1;

        let loop_mom_start = self.pairs.loop_moms_spatial.value_range.start * multiplicative_offset;

        let mut values: Vec<Complex<F<ArbPrec>>> = self.values[value_index]
            .iter()
            .map(|v| v.higher().higher())
            .collect();

        let flattened_loop_momenta = if let Some(dual_loop_moms) = &sample.sample.dual_loop_moms {
            dual_loop_moms
                .iter()
                .flat_map(|mom| {
                    vec![
                        mom.px.values.clone(),
                        mom.py.values.clone(),
                        mom.pz.values.clone(),
                    ]
                    .into_iter()
                    .flatten()
                })
                .collect_vec()
        } else {
            sample
                .loop_moms()
                .clone()
                .into_iter()
                .flat_map(|mom| vec![mom.px, mom.py, mom.pz])
                .collect_vec()
        };

        values[loop_mom_start..]
            .iter_mut()
            .zip(flattened_loop_momenta)
            .for_each(|(value, loop_mom_component)| *value = Complex::new_re(loop_mom_component));

        self.pairs
            .add_additional_params(additional_params, &mut values, multiplicative_offset);

        self.pairs.add_external_four_mom_impl(
            sample.external_moms(),
            &mut values,
            multiplicative_offset,
        );

        let polarization_values =
            self.pairs
                .polarizations_values(graph, sample.external_moms(), helicities);

        let mut polarization_start =
            self.pairs.polarizations.value_range.start * multiplicative_offset;

        for val in polarization_values {
            values[polarization_start] = val;
            polarization_start += multiplicative_offset;
        }

        if let Some(threshold_params) = left_threshold_params {
            self.pairs
                .left_threshold_params(threshold_params, &mut values);
        }

        if let Some(threshold_params) = right_threshold_params {
            self.pairs
                .right_threshold_params(threshold_params, &mut values);
        }

        if let Some(lu_params) = lu_params {
            self.pairs.lu_params(lu_params, &mut values);
        }

        InputParams {
            values: SliceMut::Owned(values),
            multiplicative_offset,
            override_pos: self.pairs.override_if.value_range.start,
            orientations_start: self.pairs.orientations.value_range.start,
        }
    }
}

impl<T: FloatLike> ParamBuilder<T>
where
    T::Higher: FloatLike,
    T::Lower: FloatLike,
{
    fn higher(&self, value_index: usize) -> Vec<Complex<F<T::Higher>>> {
        self.values[value_index]
            .iter()
            .map(|v| v.higher())
            .collect()
    }
    // fn lower(&self) -> Vec<Complex<F<T::Lower>>> {
    //     self.values.iter().map(|v| v.lower()).collect()
    // }
}
impl<T: FloatLike> ParamBuilder<T> {
    pub fn model_values(&self) -> &[Complex<F<T>>] {
        let range = self.pairs.model_parameters.value_range.clone();
        &self.values[0][range]
    }
    pub fn validate(&self) {
        self.pairs.validate();
    }

    pub fn add_tagged_function<A: Into<Indeterminate> + Clone>(
        &mut self,
        name: Symbol,
        tags: Vec<Atom>,
        rename: String,
        args: Vec<A>,
        body: Atom,
    ) -> Result<(), String> {
        let args = args.into_iter().map(|a| a.into()).collect_vec();
        let atom_args = args.iter().map(|a| Atom::from(a.clone())).collect_vec();

        self.reps.push(FnMapEntry {
            lhs: FunctionBuilder::new(name)
                .add_args(&tags)
                .add_args(&atom_args)
                .finish(),
            rhs: body.clone(),
            tags: tags.clone(),
            args: args.clone().into_iter().map(|a| a.into()).collect(),
        });

        self.fn_map
            .add_tagged_function(name, tags, rename, args, body)

        // body.evaluate(coeff_map, const_map, function_map)
    }

    pub fn add_function<A: Into<Indeterminate> + Clone>(
        &mut self,
        name: Symbol,
        rename: String,
        args: Vec<A>,
        body: Atom,
    ) -> Result<(), String> {
        let args = args.into_iter().map(|a| a.into()).collect_vec();
        let atom_args = args.iter().map(|a| Atom::from(a.clone())).collect_vec();
        self.reps.push(FnMapEntry {
            lhs: FunctionBuilder::new(name).add_args(&atom_args).finish(),
            rhs: body.clone(),
            tags: vec![],
            args: args.clone().into_iter().map(|a| a.into()).collect(),
        });

        self.fn_map.add_function(name, rename, args, body)
    }

    pub fn initialize_t_derivatives(&mut self, num_derivatives: usize) {
        let mut higher_order_values = (1..=num_derivatives)
            .map(|derivative_order| {
                self.values
                    .first()
                    .unwrap()
                    .clone()
                    .into_iter()
                    .map(move |value| {
                        let mut constant_dual_values = vec![value.clone()];
                        for _ in 0..derivative_order {
                            constant_dual_values.push(value.ref_zero());
                        }
                        constant_dual_values
                    })
                    .flatten()
                    .collect()
            })
            .collect();
        self.values.append(&mut higher_order_values);
        println!("Initialized {} derivative values", self.values.len() - 1);

        for values in self.values.iter() {
            println!("values at order: ");
            println!("len: {}", values.len());
        }
    }

    pub fn add_constant(&mut self, key: Atom, value: symbolica::domains::float::Complex<Rational>) {
        self.reps.push(FnMapEntry {
            lhs: key.clone(),
            rhs: Atom::num(value.clone()),
            tags: vec![],
            args: vec![],
        });

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

    pub(crate) fn new<
        'a,
        A: Into<AtomOrView<'a>>,
        P: IntoIterator<Item = A>,
        G: ParamBuilderGraph + SplitPolarizations,
    >(
        graph: &G,
        model: &Model,
        lmb: &LoopMomentumBasis,
        additional_params: P,
    ) -> Self {
        let (pairs, len) = GammaLoopPairs::new(model, graph, lmb, additional_params);

        let mut new = Self {
            pairs,
            ..Default::default()
        };

        let pi_rational = Rational::try_from(std::f64::consts::PI).unwrap();

        for e in graph.iter_edge_ids() {
            new.add_tagged_function::<Symbol>(
                GS.ose,
                vec![Atom::num(e.0 as i64)],
                format!("OSE{e}"),
                vec![],
                graph.explicit_ose_atom(e),
            )
            .unwrap();
        }

        for (edge_id, signature) in lmb.edge_signatures.iter() {
            if !lmb.loop_edges.contains(&edge_id) {
                let start = if signature.internal.iter().any(|sign| sign.is_sign()) {
                    //has loop mom->is a non-tree edge -> energy is OSE-> no need for Q(0) rep
                    1
                } else {
                    0
                };
                for i in start..4 {
                    new.add_tagged_function::<Symbol>(
                        GS.emr_mom,
                        vec![
                            Atom::num(edge_id.0 as i64),
                            Atom::from(ExpandedIndex::from_iter([i])),
                        ],
                        format!("Q({edge_id}, {i})"),
                        vec![],
                        lmb.loop_atom(
                            edge_id,
                            GS.emr_mom,
                            &[Atom::from(ExpandedIndex::from_iter([i]))],
                            true,
                        ) + lmb.ext_atom(
                            edge_id,
                            GS.emr_mom,
                            &[Atom::from(ExpandedIndex::from_iter([i]))],
                            true,
                        ),
                    )
                    .unwrap();
                }
            }
        }

        new.add_function(
            GS.broadcasting_sqrt,
            "sqrt".to_string(),
            vec![symbol!("x")],
            parse_lit!(sqrt(x)),
        )
        .unwrap();

        // new.fn_map.add_conditional(GS.orientation_if);

        new.add_constant(GS.pi.into(), pi_rational.into());

        new.values = vec![vec![Complex::new_re(F(T::from_f64(0.))); len]];
        new.update_model_values(model);
        new.update_idenso_values();

        //println!("self: {}", new);
        //panic!();
        new
    }

    #[inline]
    pub(crate) fn m_uv_value(&mut self, m_uv: Complex<F<T>>) {
        debug_assert!(self.pairs.m_uv.value_range.len() == 1);

        for (index, values) in self.values.iter_mut().enumerate() {
            let multiplicative_offset = index + 1;
            values[self.pairs.m_uv.value_range.start * multiplicative_offset] = m_uv.clone();
        }
    }

    pub(crate) fn mu_r_sq_value(&mut self, mu_r_sq: Complex<F<T>>) {
        debug_assert!(self.pairs.mu_r_sq.value_range.len() == 1);

        for (index, values) in self.values.iter_mut().enumerate() {
            let multiplicative_offset = index + 1;
            values[self.pairs.mu_r_sq.value_range.start * multiplicative_offset] = mu_r_sq.clone();
        }
    }

    pub(crate) fn update_idenso_values(&mut self) {
        let tr_value = Complex::new_re(F(T::from_f64(0.5)));
        let nc_value = Complex::new_re(F(T::from_f64(3.)));
        debug_assert!(self.pairs.idenso_vars.value_range.len() == 2);

        for (index, values) in self.values.iter_mut().enumerate() {
            let multiplicative_offset = index + 1;
            values[self.pairs.idenso_vars.value_range.start * multiplicative_offset] =
                tr_value.clone();
            values[self.pairs.idenso_vars.value_range.start * multiplicative_offset
                + multiplicative_offset] = nc_value.clone();
        }
    }

    pub(crate) fn update_model_values(&mut self, model: &Model) {
        for (value_index, values) in self.values.iter_mut().enumerate() {
            let multiplicative_offset = value_index + 1;
            let mut pos = self.pairs.model_parameters.value_range.start * multiplicative_offset;
            let value_index = multiplicative_offset - 1;
            for cpl in model.couplings.values().filter(|c| c.value.is_some()) {
                if let Some(value) = cpl.value {
                    values[pos] = value.map(F::from_f64);
                    pos += multiplicative_offset;
                }
            }
            for param in model.parameters.values().filter(|p| p.value.is_some()) {
                if let Some(value) = param.value {
                    let value =
                        Complex::new(F::<T>::from_ff64(value.re), F::<T>::from_ff64(value.im));
                    values[pos] = value.clone();
                    pos += multiplicative_offset;
                }
            }

            debug_assert_eq!(
                pos,
                self.pairs.model_parameters.value_range.end * multiplicative_offset
            );
        }
    }

    pub(crate) fn add_external_four_mom(
        &mut self,
        ext: &ExternalFourMomenta<F<T>>,
        multiplicative_offset: usize,
    ) {
        let value_index = multiplicative_offset - 1;
        self.pairs.add_external_four_mom_impl(
            ext,
            &mut self.values[value_index],
            multiplicative_offset,
        );
    }

    pub(crate) fn add_external_four_mom_all_derivatives(
        &mut self,
        ext: &ExternalFourMomenta<F<T>>,
    ) {
        for value_index in 0..self.values.len() {
            let multiplicative_offset = value_index + 1;
            self.pairs.add_external_four_mom_impl(
                ext,
                &mut self.values[value_index],
                multiplicative_offset,
            );
        }
    }

    pub(crate) fn orientation_value<O: GraphOrientation>(
        &mut self,
        orientation: &O,
        multiplicative_offset: usize,
    ) {
        let zero: Complex<F<T>> = Complex::new_re(F(T::from_f64(0.)));
        let one = zero.ref_one();
        let minusone = -(one.clone());

        let mut o_start = self.pairs.orientations.value_range.start * multiplicative_offset;
        let value_index = multiplicative_offset - 1;

        for (_, i) in orientation.orientation() {
            match i {
                Orientation::Default => {
                    self.values[value_index][o_start] = one.clone();
                    o_start += multiplicative_offset;
                }
                Orientation::Reversed => {
                    self.values[value_index][o_start] = minusone.clone();
                    o_start += multiplicative_offset;
                }
                Orientation::Undirected => {
                    self.values[value_index][o_start] = zero.clone();
                    o_start += multiplicative_offset;
                }
            }
        }
    }

    pub(crate) fn set_override_if(&mut self, over_ride: bool, multiplicative_offset: usize) {
        let zero: Complex<F<T>> = Complex::new_re(F(T::from_f64(0.)));
        let one = zero.ref_one();

        let o_start = self.pairs.override_if.value_range.start * multiplicative_offset;
        let value_index = multiplicative_offset - 1;

        if over_ride {
            self.values[value_index][o_start] = zero;
        } else {
            self.values[value_index][o_start] = one;
        }
    }

    pub(crate) fn polarizations_values(
        &mut self,
        (cache, debug_cache): (bool, bool),
        graph: &Graph,
        sample: &MomentumSample<T>,
        helicities: &[Helicity],
        multiplicative_offset: usize,
    ) {
        let cache_id = sample.sample.external_mom_cache_id;
        let base_cache_id = sample.sample.external_mom_base_cache_id;
        let value_index = multiplicative_offset - 1;

        // Validate cache consistency
        if cache_id < base_cache_id {
            warn!(
                "WARNING: Cache ID inconsistency detected! current={}, base={}",
                cache_id, base_cache_id
            );
        }

        // Compute expected polarizations for debug search
        let expected_pols = if debug_cache {
            Some(
                self.pairs
                    .polarizations_values(graph, sample.external_moms(), helicities),
            )
        } else {
            None
        };

        // Try to get from cache first
        let pols = if let Some(v) = self.polarization_cache.get(cache_id) {
            debug!(
                "Cache HIT for external_cache_id={} (base={})",
                cache_id, base_cache_id
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
                let mut polarization_start =
                    self.pairs.polarizations.value_range.start * multiplicative_offset;

                for pol_value in v {
                    self.values[value_index][polarization_start] = pol_value.clone();
                    polarization_start += multiplicative_offset;
                }
                return;
            }
            v
        } else {
            debug!(
                "Cache MISS for external_cache_id={} (base={})",
                cache_id, base_cache_id
            );

            // DEBUG: Search entire cache for matching polarizations to detect missed hits
            if debug_cache
                && cache
                && let Some(ref expected) = expected_pols
            {
                self.debug_search_cache_for_matching_polarizations(
                    expected,
                    cache_id,
                    base_cache_id,
                    sample.external_moms(),
                );
            }

            if !cache {
                let computed_pols = expected_pols.unwrap_or_else(|| {
                    self.pairs
                        .polarizations_values(graph, sample.external_moms(), helicities)
                });

                let mut polarization_start =
                    self.pairs.polarizations.value_range.start * multiplicative_offset;

                for pol_value in &computed_pols {
                    self.values[value_index][polarization_start] = pol_value.clone();
                    polarization_start += multiplicative_offset;
                }

                return;
            }

            // Compute and cache new polarizations
            let computed_pols = expected_pols.unwrap_or_else(|| {
                self.pairs
                    .polarizations_values(graph, sample.external_moms(), helicities)
            });

            self.polarization_cache
                .checked_push(cache_id, computed_pols);
            debug!("Cached new polarizations for cache_id={}", cache_id);

            self.polarization_cache.get(cache_id).unwrap()
        };

        let mut polarization_start =
            self.pairs.polarizations.value_range.start * multiplicative_offset;

        for pol_value in pols {
            self.values[value_index][polarization_start] = pol_value.clone();
            polarization_start += multiplicative_offset;
        }
    }

    /// Debug function to search cache for matching polarizations to detect missed cache hits
    fn debug_search_cache_for_matching_polarizations(
        &self,
        expected_pols: &[Complex<F<T>>],
        current_cache_id: usize,
        base_cache_id: usize,
        external_moms: &ExternalFourMomenta<F<T>>,
    ) {
        let tolerance = F::<T>::from_f64(1e-12);
        let mut found_matches = Vec::new();

        // Search through all possible cache IDs that could contain our polarizations
        let max_search_id = std::cmp::max(current_cache_id + 10, self.polarization_cache.len);

        for search_id in 0..max_search_id {
            if let Some(cached_pols) = self.polarization_cache.get(search_id)
                && cached_pols.len() == expected_pols.len()
            {
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

        if !found_matches.is_empty() {
            warn!(
                "🔍 DEBUG CACHE SEARCH: Found {} matching polarization(s) in cache but using cache_id={}!",
                found_matches.len(),
                current_cache_id
            );
            warn!(
                "   ⚠️  MISSED CACHE HITS: Found identical polarizations at cache_id(s): {:?}",
                found_matches
            );
            warn!(
                "   📊 Current: cache_id={}, base_cache_id={}",
                current_cache_id, base_cache_id
            );

            // Provide actionable debugging information
            if found_matches.contains(&base_cache_id) {
                warn!(
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
                warn!(
                    "   🔗 RELATED IDs: Cache IDs {:?} are related to base_id {} and contain identical polarizations",
                    related_ids, base_cache_id
                );
            }

            // Log external momentum info for debugging
            debug!(
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
                debug!(
                    "🔍 DEBUG CACHE SEARCH: No matching polarizations found in cache (searched {} entries). This is a genuine cache miss.",
                    max_search_id
                );
            }
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

    pub(crate) fn left_threshold_params(&mut self, threshold_params: &ThresholdParams<T>) {
        self.pairs
            .left_threshold_params(threshold_params, &mut self.values[0]);
    }

    pub(crate) fn right_threshold_params(&mut self, threshold_params: &ThresholdParams<T>) {
        self.pairs
            .right_threshold_params(threshold_params, &mut self.values[0]);
    }

    pub fn table(&self) -> Table {
        let mut table = tabled::builder::Builder::new();

        for FnMapEntry { lhs, rhs, .. } in &self.reps {
            table.push_record(vec![
                lhs.printer(LOGPRINTOPTS).to_string(),
                rhs.printer(LOGPRINTOPTS).to_string(),
            ]);
        }

        for i in &self.pairs {
            for (p, v) in i.params.iter().zip(i.value_range.clone()) {
                table.push_record(vec![
                    p.printer(LOGPRINTOPTS).to_string().to_string(),
                    self.values[0][v].to_string(),
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

        for FnMapEntry { lhs, rhs, .. } in &self.reps {
            map.insert(lhs.to_canonical_string(), rhs.to_canonical_string().into());
        }

        for pairs in &self.pairs {
            for (param, value) in pairs.params.iter().zip(pairs.value_range.clone()) {
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
