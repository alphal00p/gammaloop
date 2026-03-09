use std::{borrow::Cow, fmt::Display, path::Path, sync::Arc};

use ahash::{AHashMap, AHashSet, HashMap};
use eyre::eyre;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate, KeyLookup, Symbol},
    coefficient::ConvertToRing,
    domains::{
        EuclideanDomain, InternalOrdering,
        factorized_rational_polynomial::{
            FactorizedRationalPolynomial, FromNumeratorAndFactorizedDenominator,
        },
        float::{Complex as SymComplex, FloatLike, Real, SingleFloat},
        rational::Rational,
        rational_polynomial::{FromNumeratorAndDenominator, RationalPolynomial},
    },
    evaluate::{
        CompileOptions, CompiledCode, CompiledNumber, EvalTree, EvaluationFn, ExportNumber,
        ExportSettings, ExportedCode, ExpressionEvaluator, FunctionMap, OptimizationSettings,
    },
    id::{BorrowReplacement, Context, Pattern},
    poly::{
        PolyVariable, PositiveExponent, factor::Factorize, gcd::PolynomialGCD,
        polynomial::MultivariatePolynomial,
    },
    symbol,
    utils::{BorrowedOrOwned, Settable},
};

use crate::{
    algebra::{
        complex::{Complex, symbolica_traits::CompiledComplexEvaluatorSpenso},
        upgrading_arithmetic::TrySmallestUpgrade,
    },
    iterators::IteratableTensor,
    shadowing::{
        ShadowMapping, Shadowable,
        symbolica_utils::{IntoArgs, IntoSymbol},
    },
    structure::{
        HasName, TensorStructure, ToSymbolic,
        slot::{AbsInd, IsAbstractSlot, ParseableAind},
    },
    tensors::{
        complex::RealOrComplexTensor,
        data::DataTensor,
        parametric::{
            AtomViewOrConcrete, CompiledEvalTensor, EvalTensor, EvalTreeTensor, MixedTensor,
            ParamTensor,
            atomcore::{PatternReplacement, ReplaceBuilderGeneric, TensorAtomMaps, TensorAtomOps},
        },
    },
};

use super::{
    ExecutionResult, Network, TensorNetworkError,
    store::{NetworkStore, TensorScalarStore, TensorScalarStoreMapping},
};

impl<'a> From<ExecutionResult<Cow<'a, Atom>>> for Atom {
    fn from(value: ExecutionResult<Cow<'a, Atom>>) -> Self {
        match value {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(v) => v.into_owned(),
        }
    }
}

impl<
    Store: TensorScalarStore<Tensor = MixedTensor<T, S>, Scalar: AtomCore>,
    T,
    S,
    K,
    FK,
    Aind: AbsInd,
> Network<Store, K, FK, Aind>
where
    S: TensorStructure + Clone,
    T: Clone,
{
    pub fn replace<'b, P: Into<BorrowedOrOwned<'b, Pattern>>>(
        &self,
        pattern: P,
    ) -> ReplaceBuilderGeneric<'b, &'_ Self, Self> {
        ReplaceBuilderGeneric::new(self, pattern)
    }

    #[allow(clippy::type_complexity)]
    pub fn replace_multiple<R: BorrowReplacement>(
        &self,
        replacements: &[R],
    ) -> Network<Store::Store<MixedTensor<T, S>, Atom>, K, FK, Aind>
    where
        K: Clone,
        FK: Clone,
    {
        // println!("Replacing");
        self.map_ref(
            |a| a.replace_multiple(replacements),
            |a| a.replace_multiple(replacements),
        )
    }

    pub fn generate_params(&mut self) -> AHashSet<Atom>
    where
        K: Clone,
        FK: Clone,
    {
        let mut params = AHashSet::new();
        for n in self.iter_tensors().filter(|t| t.is_parametric()) {
            for (_, a) in n.iter_flat() {
                if let AtomViewOrConcrete::Atom(a) = a {
                    params.insert(a.to_owned());
                }
            }
        }
        params
    }
}

impl<
    Store: TensorScalarStore<Tensor = T, Scalar = S>,
    T,
    S,
    K: Clone,
    FK: Clone,
    Aind: AbsInd + ParseableAind,
> Network<Store, K, FK, Aind>
where
    T: Shadowable + HasName<Name = Symbol, Args: IntoArgs>,
    T::Structure: Clone + ToSymbolic + TensorStructure,
    <T::Structure as TensorStructure>::Slot: IsAbstractSlot<Aind = Aind>,
{
    #[allow(clippy::type_complexity)]
    pub fn sym_shadow(
        &mut self,
        tensor_name: &str,
        scalar_name: &str,
    ) -> Network<Store::Store<ParamTensor<T::Structure>, Atom>, K, FK, Aind> {
        self.map_ref_mut_enumerate(
            |(i, _)| Atom::var(symbol!(format!("{}{}", scalar_name, i))),
            |(i, t)| {
                t.set_name(symbol!(format!("{}{}", tensor_name, i)));

                let node = t.expanded_shadow().unwrap();
                ParamTensor::<T::Structure>::param(node.clone().into())
            },
        )
    }
}

impl<
    Store: TensorScalarStore<Tensor = T, Scalar = S>,
    T,
    S,
    K: Clone,
    FK: Clone,
    Aind: AbsInd + ParseableAind,
> Network<Store, K, FK, Aind>
where
    T: HasName<Name: IntoSymbol, Args: IntoArgs> + TensorStructure,
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    pub fn append_map<U>(&self, fn_map: &mut FunctionMap<U>)
    where
        T: ShadowMapping<U>,
        T::Structure: Clone + ToSymbolic + TensorStructure,
        <T::Structure as TensorStructure>::Slot: IsAbstractSlot<Aind = Aind>,
        S: Clone,
    {
        for n in self.iter_tensors() {
            n.expanded_append_map(fn_map)
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn shadow(&self) -> Network<Store::Store<ParamTensor<T::Structure>, S>, K, FK, Aind>
    where
        T: Shadowable,
        T::Structure: Clone + ToSymbolic + TensorStructure<Slot = T::Slot>,
        S: Clone,
        K: Clone,
        FK: Clone,
    {
        self.map_ref(Clone::clone, |t| {
            let node = t.expanded_shadow().unwrap();
            ParamTensor::<T::Structure>::param(node.clone().into())
        })
    }
}

impl<Store: TensorScalarStore<Tensor = T, Scalar = S>, T, S, FK: Clone, K: Clone, Aind: AbsInd>
    Network<Store, K, FK, Aind>
where
    T: HasName + TensorStructure,
{
    pub fn name(&mut self, name: T::Name)
    where
        T::Name: From<std::string::String> + Display,
    {
        for (id, t) in self.iter_tensors_mut().enumerate() {
            t.set_name(format!("{}{}", name, id).into());
        }
    }

    pub fn namesym(&mut self, name: &str)
    where
        T::Name: IntoSymbol,
    {
        for (id, t) in self.iter_tensors_mut().enumerate() {
            t.set_name(IntoSymbol::from_str(&format!("{}{}", name, id)));
        }
    }
}

impl<
    Store: TensorScalarStore<Tensor = P, Scalar = Atom>,
    P: TensorAtomMaps,
    K: Clone,
    FK: Clone,
    Aind: AbsInd,
> TensorAtomMaps for Network<Store, K, FK, Aind>
where
    P: Clone + TensorStructure,
{
    type ContainerData<T> = Network<Store::Store<P::ContainerData<T>, T>, K, FK, Aind>;
    type AtomContainer = Network<Store::Store<P::AtomContainer, Atom>, K, FK, Aind>;
    type Ref<'a>
        = &'a Self
    where
        Self: 'a;

    // fn replace_all<R: symbolica::id::BorrowPatternOrMap>(
    //     &self,
    //     pattern: &Pattern,
    //     rhs: R,
    //     conditions: Option<&Condition<PatternRestriction>>,
    //     settings: Option<&MatchSettings>,
    // ) -> Self {
    //     let rhs = rhs.borrow();
    //     let graph = self
    //         .graph
    //         .map_nodes_ref(|(_, a)| a.replace_all(pattern, rhs, conditions, settings));
    //     TensorNetwork {
    //         graph,
    //         scalar: self
    //             .scalar
    //             .as_ref()
    //             .map(|a| a.replace_all(pattern, rhs, conditions, settings)),
    //     }
    // }
    fn replace<'b, Pat: Into<symbolica::utils::BorrowedOrOwned<'b, Pattern>>>(
        &self,
        pattern: Pat,
    ) -> crate::tensors::parametric::atomcore::ReplaceBuilderGeneric<'b, Self::Ref<'_>, Self>
    where
        Self: Sized,
    {
        ReplaceBuilderGeneric::new(self, pattern)
    }

    fn apart(&self, x: Symbol) -> Self::AtomContainer {
        self.map_ref(|a| a.apart(x), |a| a.apart(x))
    }

    fn cancel(&self) -> Self::AtomContainer {
        self.map_ref(|a| a.cancel(), |a| a.cancel())
    }

    fn expand(&self) -> Self::AtomContainer {
        self.map_ref(|a| a.expand(), |a| a.expand())
    }

    fn factor(&self) -> Self::AtomContainer {
        self.map_ref(|a| a.factor(), |a| a.factor())
    }

    fn nsolve<
        'a,
        N: SingleFloat + Real + PartialOrd,
        V: Into<BorrowedOrOwned<'a, Indeterminate>>,
    >(
        &self,
        x: V,
        init: N,
        prec: N,
        max_iterations: usize,
    ) -> std::result::Result<Self::ContainerData<N>, std::string::String> {
        let binding = x.into();
        let x = binding.borrow();
        self.map_ref_result(
            |a| a.nsolve(x, init.clone(), prec.clone(), max_iterations),
            |a| a.nsolve(x, init.clone(), prec.clone(), max_iterations),
        )
    }

    fn series<T: AtomCore>(
        &self,
        x: Symbol,
        expansion_point: T,
        depth: Rational,
        depth_is_absolute: bool,
    ) -> std::result::Result<
        Self::ContainerData<symbolica::poly::series::Series<symbolica::domains::atom::AtomField>>,
        String,
    > {
        self.map_ref_result(
            |a| {
                a.series(
                    x,
                    expansion_point.as_atom_view(),
                    depth.clone(),
                    depth_is_absolute,
                )
            },
            |a| {
                a.series(
                    x,
                    expansion_point.as_atom_view(),
                    depth.clone(),
                    depth_is_absolute,
                )
            },
        )
    }

    fn evaluate<A: AtomCore + KeyLookup, T: Real, F: Fn(&Rational) -> T + Copy>(
        &self,
        coeff_map: F,
        const_map: &HashMap<A, T>,
        function_map: &HashMap<Symbol, EvaluationFn<A, T>>,
        // cache: &mut HashMap<AtomView<'b>, T>,
    ) -> std::result::Result<Self::ContainerData<T>, std::string::String> {
        self.map_ref_result(
            |a| a.evaluate(coeff_map, const_map, function_map),
            |a| a.evaluate(coeff_map, const_map, function_map),
        )
    }

    fn together(&self) -> Self::AtomContainer {
        self.map_ref(AtomCore::together, TensorAtomMaps::together)
    }

    fn expand_in<T: AtomCore>(&self, var: T) -> Self::AtomContainer {
        let var = var.as_atom_view();
        self.map_ref(|a| a.expand_in(var), |a| a.expand_in(var))
    }

    fn map_terms(
        &self,
        f: impl Fn(AtomView) -> Atom + Send + Sync + Clone,
        n_cores: usize,
    ) -> Self::AtomContainer {
        self.map_ref(
            |a| a.map_terms(f.clone(), n_cores),
            |a| a.map_terms(f.clone(), n_cores),
        )
    }

    fn zero_test(
        &self,
        iterations: usize,
        tolerance: f64,
    ) -> Self::ContainerData<symbolica::id::ConditionResult> {
        self.map_ref(
            |a| a.zero_test(iterations, tolerance),
            |a| a.zero_test(iterations, tolerance),
        )
    }

    fn derivative(&self, x: Symbol) -> Self::AtomContainer {
        self.map_ref(|a| a.derivative(x), |a| a.derivative(x))
    }

    fn expand_num(&self) -> Self::AtomContainer {
        self.map_ref(AtomCore::expand_num, TensorAtomMaps::expand_num)
    }

    fn to_pattern(&self) -> Self::ContainerData<Pattern> {
        self.map_ref(AtomCore::to_pattern, TensorAtomMaps::to_pattern)
    }

    fn coefficient<T: AtomCore>(&self, x: T) -> Self::AtomContainer {
        let x = x.as_atom_view();
        self.map_ref(|a| a.coefficient(x), |a| a.coefficient(x))
    }

    fn collect_num(&self) -> Self::AtomContainer {
        self.map_ref(AtomCore::collect_num, TensorAtomMaps::collect_num)
    }

    fn replace_map<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(
        &self,
        m: &F,
    ) -> Self::AtomContainer {
        self.map_ref(|a| a.replace_map(m), |a| a.replace_map(m))
    }

    fn to_polynomial<R: EuclideanDomain + ConvertToRing, E: symbolica::poly::Exponent>(
        &self,
        field: &R,
        var_map: Option<Arc<Vec<PolyVariable>>>,
    ) -> Self::ContainerData<MultivariatePolynomial<R, E>> {
        self.map_ref(
            |a| a.to_polynomial(field, var_map.clone()),
            |a| a.to_polynomial(field, var_map.clone()),
        )
    }

    fn expand_via_poly<E: symbolica::poly::Exponent, T: AtomCore>(
        &self,
        var: Option<T>,
    ) -> Self::AtomContainer {
        let var = var.as_ref().map(|a| a.as_atom_view());

        self.map_ref(
            |a| a.expand_via_poly::<E, AtomView>(var),
            |a| a.expand_via_poly::<E, AtomView>(var),
        )
    }

    fn expand_in_symbol(&self, var: Symbol) -> Self::AtomContainer {
        self.map_ref(|a| a.expand_in_symbol(var), |a| a.expand_in_symbol(var))
    }

    fn map_coefficient<
        F: Fn(symbolica::coefficient::CoefficientView) -> symbolica::coefficient::Coefficient + Copy,
    >(
        &self,
        f: F,
    ) -> Self::AtomContainer {
        self.map_ref(|a| a.map_coefficient(f), |a| a.map_coefficient(f))
    }

    fn replace_map_mut<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(&mut self, m: &F) {
        self.iter_tensors_mut().for_each(|a| a.replace_map_mut(m));

        for a in self.iter_scalars_mut() {
            a.replace_map_mut(m);
        }
    }

    fn to_canonical_string(&self) -> Self::ContainerData<std::string::String> {
        self.map_ref(
            AtomCore::to_canonical_string,
            TensorAtomMaps::to_canonical_string,
        )
    }

    fn replace_multiple<T: BorrowReplacement>(&self, replacements: &[T]) -> Self::AtomContainer {
        self.map_ref(
            |a| a.replace_multiple(replacements),
            |a| a.replace_multiple(replacements),
        )
    }

    fn set_coefficient_ring(&self, vars: &Arc<Vec<PolyVariable>>) -> Self::AtomContainer {
        self.map_ref(
            |a| a.set_coefficient_ring(vars),
            |a| a.set_coefficient_ring(vars),
        )
    }

    fn to_float(&self, f: u32) -> Self::AtomContainer {
        self.map_ref(|a| a.to_float(f), |a| a.to_float(f))
    }

    fn map_terms_single_core(&self, f: impl Fn(AtomView) -> Atom + Clone) -> Self::AtomContainer {
        self.map_ref(
            |a| a.map_terms_single_core(f.clone()),
            |a| a.map_terms_single_core(f.clone()),
        )
    }

    fn to_polynomial_in_vars<E: symbolica::poly::Exponent>(
        &self,
        var_map: &Arc<Vec<PolyVariable>>,
    ) -> Self::ContainerData<MultivariatePolynomial<symbolica::domains::atom::AtomField, E>> {
        self.map_ref(
            |a| a.to_polynomial_in_vars(var_map),
            |a| a.to_polynomial_in_vars(var_map),
        )
    }

    fn to_rational_polynomial<
        R: EuclideanDomain + ConvertToRing,
        RO: EuclideanDomain + PolynomialGCD<E>,
        E: PositiveExponent,
    >(
        &self,
        field: &R,
        out_field: &RO,
        var_map: Option<Arc<Vec<PolyVariable>>>,
    ) -> Self::ContainerData<RationalPolynomial<RO, E>>
    where
        RationalPolynomial<RO, E>:
            FromNumeratorAndDenominator<R, RO, E> + FromNumeratorAndDenominator<RO, RO, E>,
    {
        self.map_ref(
            |a| a.to_rational_polynomial(field, out_field, var_map.clone()),
            |a| a.to_rational_polynomial(field, out_field, var_map.clone()),
        )
    }

    fn rationalize(&self, relative_error: &Rational) -> Self::AtomContainer {
        self.map_ref(
            |a| a.rationalize(relative_error),
            |a| a.rationalize(relative_error),
        )
    }

    fn replace_multiple_mut<T: BorrowReplacement>(&mut self, replacements: &[T]) {
        self.iter_tensors_mut()
            .for_each(|a| a.replace_multiple_mut(replacements));
        for a in self.iter_scalars_mut() {
            a.replace_multiple_mut(replacements);
        }
    }

    fn replace_multiple_repeat<T: BorrowReplacement>(
        &self,
        replacements: &[T],
    ) -> Self::AtomContainer {
        self.map_ref(
            |a| a.replace_multiple_repeat(replacements),
            |a| a.replace_multiple_repeat(replacements),
        )
    }

    fn replace_multiple_repeat_mut<T: BorrowReplacement>(&mut self, replacements: &[T]) {
        self.iter_tensors_mut()
            .for_each(|a| a.replace_multiple_repeat_mut(replacements));
        for a in self.iter_scalars_mut() {
            a.replace_multiple_repeat_mut(replacements);
        }
    }

    fn to_factorized_rational_polynomial<
        R: EuclideanDomain + ConvertToRing,
        RO: EuclideanDomain + PolynomialGCD<E>,
        E: PositiveExponent,
    >(
        &self,
        field: &R,
        out_field: &RO,
        var_map: Option<Arc<Vec<PolyVariable>>>,
    ) -> Self::ContainerData<FactorizedRationalPolynomial<RO, E>>
    where
        FactorizedRationalPolynomial<RO, E>: FromNumeratorAndFactorizedDenominator<R, RO, E>
            + FromNumeratorAndFactorizedDenominator<RO, RO, E>,
        MultivariatePolynomial<RO, E>: Factorize,
    {
        self.map_ref(
            |a| a.to_factorized_rational_polynomial(field, out_field, var_map.clone()),
            |a| a.to_factorized_rational_polynomial(field, out_field, var_map.clone()),
        )
    }
}

pub trait Evaluate<Str: TensorScalarStore, K, FK, S, T, Aind> {
    #[allow(clippy::type_complexity)]
    fn evaluate(&mut self, params: &[T]) -> Network<Str::Store<DataTensor<T, S>, T>, K, FK, Aind>
    where
        S: TensorStructure + Clone;

    // fn evaluate_float(
    //     &mut self,
    //     params: &[f64],
    // ) -> Network<Str::Store<DataTensor<f64, S>, f64>, K>where S:TensorStructure+Clone;;

    // fn evaluate_complex(
    //     &mut self,
    //     params: &[SymComplex<f64>],
    // ) -> Network<Str::Store<DataTensor<SymComplex<f64>, S>, SymComplex<f64>>, K>where S:TensorStructure+Clone;;
}

impl<
    Store: TensorScalarStore<Tensor = ParamTensor<S>, Scalar = Sc>,
    S: TensorStructure + Clone,
    Sc: AtomCore,
    K: Clone,
    FK: Clone,
    Aind: Clone + AbsInd,
> Network<Store, K, FK, Aind>
{
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn eval_tree(
        &self,
        fn_map: &FunctionMap,
        params: &[Atom],
    ) -> Result<
        Network<
            Store::Store<EvalTreeTensor<SymComplex<Rational>, S>, EvalTree<SymComplex<Rational>>>,
            K,
            FK,
            Aind,
        >,
        String,
    >
    where
        S: TensorStructure,
    {
        self.map_ref_result(
            |a| a.to_evaluation_tree(fn_map, params),
            |a| a.to_evaluation_tree(fn_map, params),
        )
    }

    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn evaluate_direct<A: AtomCore + KeyLookup, D, F: Fn(&Rational) -> D + Copy>(
        &self,
        coeff_map: F,
        const_map: &AHashMap<A, D>,
        function_map: &HashMap<Symbol, EvaluationFn<A, D>>,
    ) -> Result<Network<Store::Store<DataTensor<D, S>, D>, K, FK, Aind>, String>
    where
        D: Clone
            + symbolica::domains::float::Real
            + for<'c> std::convert::From<&'c symbolica::domains::rational::Rational>,
    {
        self.map_ref_result(
            |a| a.evaluate(coeff_map, const_map, function_map),
            |a| a.evaluate(coeff_map, const_map, function_map),
        )
    }
}

impl<T: Real, S: TensorStructure, K: Clone, FK: Clone, Aind: AbsInd>
    Evaluate<NetworkStore<EvalTreeTensor<T, S>, EvalTree<T>>, K, FK, S, T, Aind>
    for Network<NetworkStore<EvalTreeTensor<T, S>, EvalTree<T>>, K, FK, Aind>
{
    fn evaluate(&mut self, params: &[T]) -> Network<NetworkStore<DataTensor<T, S>, T>, K, FK, Aind>
    where
        S: TensorStructure + Clone,
    {
        self.map_ref_mut(
            |a| {
                let mut out = [params[0].zero()];
                a.evaluate(params, &mut out);
                let [o] = out;
                o
            },
            |a| a.evaluate(params),
        )
    }
}

impl<
    Store: TensorScalarStore<
            Tensor = EvalTreeTensor<SymComplex<Rational>, S>,
            Scalar = EvalTree<SymComplex<Rational>>,
        >,
    S: TensorStructure,
    K: Clone,
    FK: Clone,
    Aind: AbsInd,
> Network<Store, K, FK, Aind>
{
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn linearize(
        self,
        settings: &OptimizationSettings,
    ) -> Network<
        Store::Store<
            EvalTensor<ExpressionEvaluator<SymComplex<Rational>>, S>,
            ExpressionEvaluator<SymComplex<Rational>>,
        >,
        K,
        FK,
        Aind,
    >
    where
        S: Clone,
    {
        self.map(|mut a| a.linearize(settings), |a| a.linearize(settings))
    }
}

impl<
    Store: TensorScalarStore<Tensor = EvalTreeTensor<T, S>, Scalar = EvalTree<T>>,
    T,
    S: TensorStructure,
    K: Clone,
    FK: Clone,
    Aind: AbsInd,
> Network<Store, K, FK, Aind>
{
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn map_coeff<T2, F: Fn(&T) -> T2>(
        &self,
        f: &F,
    ) -> Network<Store::Store<EvalTreeTensor<T2, S>, EvalTree<T2>>, K, FK, Aind>
    where
        T: Clone + PartialEq,
        S: Clone,
    {
        self.map_ref(|a| a.map_coeff(f), |a| a.map_coeff(f))
        // self.map_data_ref(|x| x.map_coeff(f))
    }

    pub fn common_subexpression_elimination(&mut self)
    where
        T: std::fmt::Debug + std::hash::Hash + Eq + InternalOrdering + Clone + Default,
        S: Clone,
    {
        self.iter_tensors_mut()
            .for_each(|a| a.common_subexpression_elimination());
        for a in self.iter_scalars_mut() {
            a.common_subexpression_elimination()
        }
    }
}

impl<T: Real, S: TensorStructure + Clone, K: Clone, FK: Clone, Aind: AbsInd>
    Evaluate<
        NetworkStore<EvalTensor<ExpressionEvaluator<T>, S>, ExpressionEvaluator<T>>,
        K,
        FK,
        S,
        T,
        Aind,
    >
    for Network<
        NetworkStore<EvalTensor<ExpressionEvaluator<T>, S>, ExpressionEvaluator<T>>,
        K,
        FK,
        Aind,
    >
{
    fn evaluate(&mut self, params: &[T]) -> Network<NetworkStore<DataTensor<T, S>, T>, K, FK, Aind>
    where
        S: TensorStructure + Clone,
        T: Real,
    {
        self.map_ref_mut(
            |a| {
                let mut out = [params[0].zero()];
                a.evaluate(params, &mut out);
                let [o] = out;
                o
            },
            |a| a.evaluate(params),
        )
    }
}

impl<
    Store: TensorScalarStore<
            Tensor = EvalTensor<ExpressionEvaluator<T>, S>,
            Scalar = ExpressionEvaluator<T>,
        >,
    T,
    S: TensorStructure,
    FK: Clone + Display,
    K: Clone + Display,
    Aind: AbsInd,
> Network<Store, K, FK, Aind>
{
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn export_cpp<F: CompiledNumber>(
        &self,
        path: impl AsRef<Path>,
        function_name: &str,
        settings: ExportSettings,
    ) -> Result<
        Network<Store::Store<EvalTensor<ExportedCode<F>, S>, ExportedCode<F>>, K, FK, Aind>,
        TensorNetworkError<K, FK>,
    >
    where
        T: FloatLike + ExportNumber + SingleFloat,
        S: Clone,
    {
        // TODO @Lucien with the new export_cpp you are now able to put these different functions in the same file!

        self.map_ref_result_enumerate(
            |(_, a)| {
                a.export_cpp(path.as_ref(), function_name, settings.clone())
                    .map_err(|a| TensorNetworkError::InOut(a))
            },
            |(_, a)| {
                a.export_cpp(path.as_ref(), function_name, settings.clone())
                    .map_err(|a| TensorNetworkError::InOut(a))
            },
        )
    }
}

impl<
    F: CompiledNumber,
    Store: TensorScalarStore<Tensor = EvalTensor<ExportedCode<F>, S>, Scalar = ExportedCode<F>>,
    S: TensorStructure,
    K: Clone + Display,
    FK: Clone + Display,
    Aind: AbsInd,
> Network<Store, K, FK, Aind>
{
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn compile(
        &self,
        out: &str,
        options: CompileOptions,
    ) -> Result<
        Network<Store::Store<EvalTensor<CompiledCode<F>, S>, CompiledCode<F>>, K, FK, Aind>,
        TensorNetworkError<K, FK>,
    >
    where
        S: Clone,
    {
        self.map_ref_result(
            |a| {
                a.compile(out, options.clone())
                    .map_err(|a| TensorNetworkError::InOut(a))
            },
            |a| {
                a.compile(out, options.clone())
                    .map_err(|a| TensorNetworkError::InOut(a))
            },
        )
    }

    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn compile_and_load(
        &self,
        out: &str,
        options: CompileOptions,
    ) -> Result<
        Network<Store::Store<EvalTensor<F::Evaluator, S>, F::Evaluator>, K, FK, Aind>,
        TensorNetworkError<K, FK>,
    >
    where
        S: Clone,
    {
        self.map_ref_result(
            |a| {
                a.compile(out, options.clone())?
                    .load()
                    .map_err(|a| TensorNetworkError::Other(eyre!("Loading error:{a}")))
            },
            |a| {
                a.compile(out, options.clone())?
                    .load()
                    .map_err(|a| TensorNetworkError::Other(eyre!("Loading error:{a}")))
            },
        )
    }
}

impl<S: TensorStructure, K: Clone, FK: Clone, Aind: AbsInd>
    Evaluate<
        NetworkStore<CompiledEvalTensor<S>, CompiledComplexEvaluatorSpenso>,
        K,
        FK,
        S,
        Complex<f64>,
        Aind,
    >
    for Network<NetworkStore<CompiledEvalTensor<S>, CompiledComplexEvaluatorSpenso>, K, FK, Aind>
{
    fn evaluate(
        &mut self,
        params: &[Complex<f64>],
    ) -> Network<NetworkStore<DataTensor<Complex<f64>, S>, Complex<f64>>, K, FK, Aind>
    where
        S: TensorStructure + Clone,
    {
        self.map_ref_mut(
            |a| {
                let mut out = [Complex::default()];
                a.evaluate(params, &mut out);
                let [o] = out;
                o
            },
            |a| a.evaluate(params),
        )
    }
}

impl<
    Store: TensorScalarStore<
            Tensor = EvalTreeTensor<SymComplex<Rational>, S>,
            Scalar = EvalTree<SymComplex<Rational>>,
        >,
    S: Clone + TensorStructure,
    K: Clone,
    FK: Clone,
    Aind: AbsInd,
> Network<Store, K, FK, Aind>
{
    pub fn horner_scheme(&mut self) {
        self.iter_tensors_mut().for_each(|a| a.horner_scheme());
        for a in self.iter_scalars_mut() {
            a.horner_scheme()
        }
    }
}

impl<
    Store: TensorScalarStore<Tensor = MixedTensor<T, S>, Scalar = Sc>,
    T,
    S,
    Sc: AtomCore,
    K,
    FK,
    Aind: AbsInd,
> Network<Store, K, FK, Aind>
where
    S: Clone + TensorStructure + std::fmt::Debug,
    T: Clone,
    K: Clone,
    FK: Clone,
{
    pub fn evaluate_real<A: AtomCore + KeyLookup, F: Fn(&Rational) -> T + Copy>(
        &mut self,
        coeff_map: F,
        const_map: &AHashMap<A, T>,
        function_map: &HashMap<Symbol, EvaluationFn<A, T>>,
    ) where
        T: Real + for<'c> From<&'c Rational>,
    {
        for t in self.iter_tensors_mut() {
            t.evaluate_real(coeff_map, const_map, function_map);
        }
    }

    pub fn evaluate_complex<A: AtomCore + KeyLookup, F: Fn(&Rational) -> SymComplex<T> + Copy>(
        &mut self,
        coeff_map: F,
        const_map: &AHashMap<A, SymComplex<T>>,
        function_map: &HashMap<Symbol, EvaluationFn<A, SymComplex<T>>>,
    ) where
        T: Real + for<'c> From<&'c Rational>,
        SymComplex<T>: Real + for<'c> From<&'c Rational>,
    {
        for t in self.iter_tensors_mut() {
            t.evaluate_complex(coeff_map, const_map, function_map);
        }
    }

    pub fn to_fully_parametric(self) -> Network<Store::Store<ParamTensor<S>, Sc>, K, FK, Aind>
    where
        T: TrySmallestUpgrade<Atom, LCM = Atom>,
        Complex<T>: TrySmallestUpgrade<Atom, LCM = Atom>,
    {
        self.map(
            |s| s,
            |t| match t {
                MixedTensor::Concrete(RealOrComplexTensor::Real(t)) => {
                    ParamTensor::composite(t.try_upgrade::<Atom>().unwrap().into_owned())
                }
                MixedTensor::Concrete(RealOrComplexTensor::Complex(t)) => {
                    ParamTensor::composite(t.try_upgrade::<Atom>().unwrap().into_owned())
                }
                MixedTensor::Param(t) => t.clone(),
            },
        )
    }
}
