use std::{marker::PhantomData, sync::Arc};

use ahash::{HashMap, HashSet, HashSetExt};
use dyn_clone::DynClone;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate, KeyLookup, Symbol},
    coefficient::{Coefficient, CoefficientView, ConvertToRing},
    domains::{
        EuclideanDomain, InternalOrdering,
        atom::AtomField,
        factorized_rational_polynomial::{
            FactorizedRationalPolynomial, FromNumeratorAndFactorizedDenominator,
        },
        float::{Complex as SymComplex, Real, SingleFloat},
        integer::Z,
        rational::Rational,
        rational_polynomial::{
            FromNumeratorAndDenominator, RationalPolynomial, RationalPolynomialField,
        },
    },
    evaluate::{EvaluationFn, ExpressionEvaluator, FunctionMap, OptimizationSettings},
    id::{
        BorrowReplacement, Condition, ConditionResult, Context, MatchMap, MatchSettings, Pattern,
        PatternRestriction, ReplaceBuilder,
    },
    poly::{
        Exponent, PolyVariable, PositiveExponent, factor::Factorize, gcd::PolynomialGCD,
        polynomial::MultivariatePolynomial, series::Series,
    },
    solve::SolveError,
    state::RecycledAtom,
    tensors::matrix::Matrix,
    utils::{BorrowedOrOwned, Settable},
};

use crate::{
    iterators::IteratableTensor,
    structure::{HasStructure, TensorStructure},
    tensors::data::{DataTensor, DenseTensor, SparseTensor, StorageTensor},
};

use super::{EvalTensor, EvalTreeTensor, MixedTensor, ParamOrConcrete, ParamTensor};

pub trait PatternReplacement {
    fn replace_multiple_repeat<T: BorrowReplacement>(&self, replacements: &[T]) -> Self;
    fn replace_multiple_mut<T: BorrowReplacement>(&mut self, replacements: &[T]);
    fn replace_multiple_repeat_mut<T: BorrowReplacement>(&mut self, replacements: &[T]);
    fn replace_map_mut<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(&mut self, m: &F);
}

impl PatternReplacement for Atom {
    fn replace_multiple_mut<T: BorrowReplacement>(&mut self, replacements: &[T]) {
        *self = self.as_atom_view().replace_multiple(replacements);
    }

    fn replace_multiple_repeat<T: BorrowReplacement>(&self, replacements: &[T]) -> Atom {
        let mut out = self.clone();
        let mut out_mut = out.clone();

        while out.replace_multiple_into(replacements, &mut out_mut) {
            if out == out_mut {
                break;
            }
            std::mem::swap(&mut out, &mut out_mut)
        }
        out
    }

    fn replace_multiple_repeat_mut<T: BorrowReplacement>(&mut self, replacements: &[T]) {
        let atom = self.replace_multiple_repeat(replacements);
        *self = atom;
    }

    fn replace_map_mut<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(&mut self, m: &F) {
        *self = self.replace_map(m);
    }
}

pub trait ClonableAtomMap: Fn(AtomView, &mut Atom) + DynClone {}

pub trait TensorAtomMapsTest {
    type ContainerData<T>;
}

impl<S: StorageTensor<Data = Atom>> TensorAtomMapsTest for S {
    type ContainerData<T> = <S as StorageTensor>::ContainerData<T>;
}

/// Construct a replacement by specifying the pattern and finishing it with the right-hand side
/// using [ReplaceBuilder::with], [ReplaceBuilder::with_into], or [ReplaceBuilder::iter].
#[derive(Debug, Clone)]
pub struct ReplaceBuilderGeneric<'b, R, T> {
    target: R,
    pattern: BorrowedOrOwned<'b, Pattern>,
    conditions: Option<BorrowedOrOwned<'b, Condition<PatternRestriction>>>,
    settings: MatchSettings,
    repeat: bool,
    _marker: PhantomData<T>,
}

pub trait ReplaceWithBuilder: Sized {
    type Ref<'a>: Copy
    where
        Self: 'a;
    type RefMut<'a>
    where
        Self: 'a;
    fn with<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
        with: C,
    ) -> Self;

    fn with_into<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
        with: C,
        into: &mut Self,
    ) -> bool;

    fn with_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
        replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
        rhs: M,
    ) -> Self;

    fn recycled(ref_self: Self::Ref<'_>) -> Self;
    // fn fill_in<'a>(&'a mut self, ref_self: Self::Ref<'a>);

    fn ref_mut(&mut self) -> Self::RefMut<'_>;
    fn reference(&self) -> Self::Ref<'_>;

    fn with_owned<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: ReplaceBuilderGeneric<'c, Self, Self>,
        with: C,
    ) -> Self;

    fn with_owned_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
        replacement: ReplaceBuilderGeneric<'c, Self, Self>,
        rhs: M,
    ) -> Self;

    fn with_mut<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: &'a mut ReplaceBuilderGeneric<'c, Self::RefMut<'a>, Self>,
        with: C,
    );

    fn with_mut_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
        replacement: &'a mut ReplaceBuilderGeneric<'c, Self::RefMut<'a>, Self>,
        rhs: M,
    );
}

impl<'b, R, Phantom> ReplaceBuilderGeneric<'b, R, Phantom> {
    pub fn new<T: Into<BorrowedOrOwned<'b, Pattern>>>(target: R, replacement: T) -> Self {
        ReplaceBuilderGeneric {
            target,
            pattern: replacement.into(),
            conditions: None,
            settings: MatchSettings::default(),
            repeat: false,
            _marker: PhantomData,
        }
    }

    pub fn with_target<RR, P2>(&self, target: RR) -> ReplaceBuilderGeneric<'b, RR, P2> {
        ReplaceBuilderGeneric {
            target,
            pattern: self.pattern.clone(),
            conditions: self.conditions.clone(),
            settings: self.settings.clone(),
            repeat: self.repeat,
            _marker: PhantomData,
        }
    }

    pub fn update_symbolica_builder<'c, 'a>(
        &'c self,
        builder: ReplaceBuilder<'a, 'b>,
    ) -> ReplaceBuilder<'a, 'b>
    where
        'c: 'b,
    {
        let builder = if self.repeat {
            builder.repeat()
        } else {
            builder
        };

        if let Some(conditions) = &self.conditions {
            builder.when(conditions.borrow())
        } else {
            builder
        }
        .non_greedy_wildcards(self.settings.non_greedy_wildcards.clone())
        .allow_new_wildcards_on_rhs(self.settings.allow_new_wildcards_on_rhs)
        .level_range(self.settings.level_range)
        .level_is_tree_depth(self.settings.level_is_tree_depth)
        .rhs_cache_size(self.settings.rhs_cache_size)
    }

    /// Specifies wildcards that try to match as little as possible.
    pub fn non_greedy_wildcards(mut self, non_greedy_wildcards: Vec<Symbol>) -> Self {
        self.settings.non_greedy_wildcards = non_greedy_wildcards;
        self
    }
    /// Specifies the `[min,max]` level at which the pattern is allowed to match.
    /// The first level is 0 and the level is increased when entering a function, or going one level deeper in the expression tree,
    /// depending on `level_is_tree_depth`.
    pub fn level_range(mut self, level_range: (usize, Option<usize>)) -> Self {
        self.settings.level_range = level_range;
        self
    }
    /// Determine whether a level reflects the expression tree depth or the function depth.
    pub fn level_is_tree_depth(mut self, level_is_tree_depth: bool) -> Self {
        self.settings.level_is_tree_depth = level_is_tree_depth;
        self
    }
    /// Allow wildcards on the right-hand side that do not appear in the pattern.
    pub fn allow_new_wildcards_on_rhs(mut self, allow: bool) -> Self {
        self.settings.allow_new_wildcards_on_rhs = allow;
        self
    }
    /// The maximum size of the cache for the right-hand side of a replacement.
    /// This can be used to prevent expensive recomputations.
    pub fn rhs_cache_size(mut self, rhs_cache_size: usize) -> Self {
        self.settings.rhs_cache_size = rhs_cache_size;
        self
    }

    /// Add a condition to the replacement.
    pub fn when<C: Into<BorrowedOrOwned<'b, Condition<PatternRestriction>>>>(
        mut self,
        conditions: C,
    ) -> Self {
        self.conditions = Some(conditions.into());
        self
    }

    /// Repeat the replacement until no more matches are found.
    pub fn repeat(mut self) -> Self {
        self.repeat = true;
        self
    }
}

impl<'b, 'a, Phantom> From<ReplaceBuilderGeneric<'b, AtomView<'a>, Phantom>>
    for ReplaceBuilder<'a, 'b>
{
    fn from(builder: ReplaceBuilderGeneric<'b, AtomView<'a>, Phantom>) -> Self {
        ReplaceBuilder::new(builder.target, builder.pattern)
    }
}

impl<'a, 'b, R: ReplaceWithBuilder> ReplaceBuilderGeneric<'b, R::Ref<'a>, R> {
    /// Execute the replacement by specifying the right-hand side.
    ///
    /// To use a map as a right-hand side, use [ReplaceBuilder::with_map].
    pub fn with<'c: 'a, C: Into<BorrowedOrOwned<'c, Pattern>>>(&'a self, rhs: C) -> R
    where
        'b: 'c,
    {
        // let rhs = ReplaceWith::Pattern(rhs.into());
        R::with(self, rhs)
    }

    /// Execute the replacement by specifying the right-hand side and writing the result in `out`.
    pub fn with_into<'c: 'a, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        &'a self,
        rhs: C,
        out: &mut R,
    ) -> bool
    where
        'b: 'c,
    {
        R::with_into(self, rhs, out)
    }

    /// Execute the replacement by specifying the right-hand side as a map on the matched wildcards.
    ///
    pub fn with_map<'c: 'a, M: MatchMap + 'static + Clone>(&'a self, rhs: M) -> R {
        R::with_map(self, rhs)
    }
}

pub trait TensorAtomMaps {
    type ContainerData<T>;
    type AtomContainer;
    type Ref<'a>
    where
        Self: 'a;

    fn replace<'b, P: Into<BorrowedOrOwned<'b, Pattern>>>(
        &self,
        pattern: P,
    ) -> ReplaceBuilderGeneric<'b, Self::Ref<'_>, Self>
    where
        Self: Sized;

    /// Replace all occurrences of the patterns, where replacements are tested in the order that they are given.
    fn replace_multiple<T: BorrowReplacement>(&self, replacements: &[T]) -> Self::AtomContainer;
    fn replace_multiple_repeat<T: BorrowReplacement>(
        &self,
        replacements: &[T],
    ) -> Self::AtomContainer;
    fn replace_multiple_mut<T: BorrowReplacement>(&mut self, replacements: &[T]);
    fn replace_multiple_repeat_mut<T: BorrowReplacement>(&mut self, replacements: &[T]);
    fn replace_map_mut<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(&mut self, m: &F);

    /// Collect terms involving the literal occurrence of `x`.
    fn coefficient<T: AtomCore>(&self, x: T) -> Self::AtomContainer;

    /// Write the expression over a common denominator.
    fn together(&self) -> Self::AtomContainer;

    /// Write the expression as a sum of terms with minimal denominators.
    fn apart(&self, x: Symbol) -> Self::AtomContainer;

    /// Cancel all common factors between numerators and denominators.
    /// Any non-canceling parts of the expression will not be rewritten.
    fn cancel(&self) -> Self::AtomContainer;

    /// Factor the expression over the rationals.
    fn factor(&self) -> Self::AtomContainer;

    /// Collect numerical factors by removing the numerical content from additions.
    /// For example, `-2*x + 4*x^2 + 6*x^3` will be transformed into `-2*(x - 2*x^2 - 3*x^3)`.
    ///
    /// The first argument of the addition is normalized to a positive quantity.
    fn collect_num(&self) -> Self::AtomContainer;

    /// Expand an expression. The function [AtomCore::expand_via_poly] may be faster.
    fn expand(&self) -> Self::AtomContainer;

    /// Expand the expression by converting it to a polynomial, optionally
    /// only in the indeterminate `var`. The parameter `E` should be a numerical type
    /// that fits the largest exponent in the expanded expression. Often,
    /// `u8` or `u16` is sufficient.
    fn expand_via_poly<E: Exponent, T: AtomCore>(&self, var: Option<T>) -> Self::AtomContainer;

    /// Expand an expression in the variable `var`. The function [AtomCore::expand_via_poly] may be faster.
    fn expand_in<T: AtomCore>(&self, var: T) -> Self::AtomContainer;

    /// Expand an expression in the variable `var`.
    fn expand_in_symbol(&self, var: Symbol) -> Self::AtomContainer;

    // /// Expand an expression, returning `true` iff the expression changed.
    // fn expand_into(&self, var: Option<AtomView>, out: &mut Atom) -> bool {
    //     self.as_atom_view().expand_into(var, out)
    // }

    /// Distribute numbers in the expression, for example:
    /// `2*(x+y)` -> `2*x+2*y`.
    fn expand_num(&self) -> Self::AtomContainer;

    /// Take a derivative of the expression with respect to `x`.
    fn derivative(&self, x: Symbol) -> Self::AtomContainer;

    // /// Take a derivative of the expression with respect to `x` and
    // /// write the result in `out`.
    // /// Returns `true` if the derivative is non-zero.
    // fn derivative_into(&self, x: Symbol, out: &mut Atom) -> bool {
    //     self.as_atom_view().derivative_into(x, out)
    // }

    /// Series expand in `x` around `expansion_point` to depth `depth`.
    fn series<T: AtomCore>(
        &self,
        x: Symbol,
        expansion_point: T,
        depth: Rational,
        depth_is_absolute: bool,
    ) -> Result<Self::ContainerData<Series<AtomField>>, String>;

    /// Find the root of a function in `x` numerically over the reals using Newton's method.
    fn nsolve<'a, N: SingleFloat + Real + PartialOrd, V: Into<BorrowedOrOwned<'a, Indeterminate>>>(
        &self,
        x: V,
        init: N,
        prec: N,
        max_iterations: usize,
    ) -> Result<Self::ContainerData<N>, String>;

    /// Evaluate a (nested) expression a single time.
    /// For repeated evaluations, use [Self::evaluator()] and convert
    /// to an optimized version or generate a compiled version of your expression.
    ///
    /// All variables and all user functions in the expression must occur in the map.
    fn evaluate<A: AtomCore + KeyLookup, T: Real, F: Fn(&Rational) -> T + Copy>(
        &self,
        coeff_map: F,
        const_map: &HashMap<A, T>,
        function_map: &HashMap<Symbol, EvaluationFn<A, T>>,
    ) -> Result<Self::ContainerData<T>, String>;

    /// Check if the expression could be 0, using (potentially) numerical sampling with
    /// a given tolerance and number of iterations.
    fn zero_test(&self, iterations: usize, tolerance: f64) -> Self::ContainerData<ConditionResult>;

    /// Set the coefficient ring to the multivariate rational polynomial with `vars` variables.
    fn set_coefficient_ring(&self, vars: &Arc<Vec<PolyVariable>>) -> Self::AtomContainer;

    /// Convert all coefficients to floats with a given precision `decimal_prec``.
    /// The precision of floating point coefficients in the input will be truncated to `decimal_prec`.
    fn to_float(&self, decimal_prec: u32) -> Self::AtomContainer;

    // /// Convert all coefficients to floats with a given precision `decimal_prec``.
    // /// The precision of floating point coefficients in the input will be truncated to `decimal_prec`.
    // fn coefficients_to_float_into(&self, decimal_prec: u32, out: &mut Atom) {
    //     self.as_atom_view()
    //         .coefficients_to_float_into(decimal_prec, out);
    // }

    /// Map all coefficients using a given function.
    fn map_coefficient<F: Fn(CoefficientView) -> Coefficient + Copy>(
        &self,
        f: F,
    ) -> Self::AtomContainer;

    // /// Map all coefficients using a given function.
    // fn map_coefficient_into<F: Fn(CoefficientView) -> Coefficient + Copy>(
    //     &self,
    //     f: F,
    //     out: &mut Atom,
    // ) {
    //     self.as_atom_view().map_coefficient_into(f, out);
    // }

    /// Map all floating point and rational coefficients to the best rational approximation
    /// in the interval `[self*(1-relative_error),self*(1+relative_error)]`.
    fn rationalize(&self, relative_error: &Rational) -> Self::AtomContainer;

    /// Convert the atom to a polynomial, optionally in the variable ordering
    /// specified by `var_map`. If new variables are encountered, they are
    /// added to the variable map. Similarly, non-polynomial parts are automatically
    /// defined as a new independent variable in the polynomial.
    fn to_polynomial<R: EuclideanDomain + ConvertToRing, E: Exponent>(
        &self,
        field: &R,
        var_map: Option<Arc<Vec<PolyVariable>>>,
    ) -> Self::ContainerData<MultivariatePolynomial<R, E>>;

    /// Convert the atom to a polynomial in specific variables.
    /// All other parts will be collected into the coefficient, which
    /// is a general expression.
    ///
    /// This routine does not perform expansions.
    fn to_polynomial_in_vars<E: Exponent>(
        &self,
        var_map: &Arc<Vec<PolyVariable>>,
    ) -> Self::ContainerData<MultivariatePolynomial<AtomField, E>>;

    /// Convert the atom to a rational polynomial, optionally in the variable ordering
    /// specified by `var_map`. If new variables are encountered, they are
    /// added to the variable map. Similarly, non-rational polynomial parts are automatically
    /// defined as a new independent variable in the rational polynomial.
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
            FromNumeratorAndDenominator<R, RO, E> + FromNumeratorAndDenominator<RO, RO, E>;
    /// Convert the atom to a rational polynomial with factorized denominators, optionally in the variable ordering
    /// specified by `var_map`. If new variables are encountered, they are
    /// added to the variable map. Similarly, non-rational polynomial parts are automatically
    /// defined as a new independent variable in the rational polynomial.
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
        MultivariatePolynomial<RO, E>: Factorize;

    /// Print the atom in a form that is unique and independent of any implementation details.
    ///
    /// Anti-symmetric functions are not supported.
    fn to_canonical_string(&self) -> Self::ContainerData<String>;

    /// Map the function `f` over all terms.
    fn map_terms_single_core(&self, f: impl Fn(AtomView) -> Atom + Clone) -> Self::AtomContainer;

    /// Map the function `f` over all terms, using parallel execution with `n_cores` cores.
    fn map_terms(
        &self,
        f: impl Fn(AtomView) -> Atom + Send + Sync + Clone,
        n_cores: usize,
    ) -> Self::AtomContainer;

    fn to_pattern(&self) -> Self::ContainerData<Pattern>;

    fn replace_map<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(
        &self,
        m: &F,
    ) -> Self::AtomContainer;

    // fn replace_iter<'a>(
    //     &'a self,
    //     pattern: &'a Pattern,
    //     rhs: &'a PatternOrMap,
    //     conditions: Option<&'a Condition<PatternRestriction>>,
    //     settings: Option<&'a MatchSettings>,
    // ) -> Self::ContainerData<ReplaceIterator<'a, 'a>>;
}

// impl<T: ReplaceWithBuilder + Clone + PartialEq, K: Clone, Aind: AbsInd> ReplaceWithBuilder
//     for Network<NetworkStore<T, Atom>, K, Aind>
// // where
// //     for<'a> T::Ref<'a>: Clone + PartialEq,
// //     for<'a> T::RefMut<'a>: Clone + PartialEq,
// {
//     type Ref<'a>
//         = &'a Self
//     where
//         Self: 'a;
//     type RefMut<'a>
//         = &'a mut Self
//     where
//         Self: 'a;

//     fn ref_mut(&mut self) -> Self::RefMut<'_> {
//         self
//     }

//     fn recycled(ref_self: Self::Ref<'_>) -> Self {
//         Network::zero()
//     }

//     fn with<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
//         replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
//         with: C,
//     ) -> Self {
//         let new_store = replacement.target.store.map_ref(
//             |a| ReplaceBuilder::from(replacement.with_target(a.as_atom_view())).with(with),
//             |t| {
//                 replacement
//                     .with_target::<T, T::Ref<'a>>(t.clone())
//                     .with(with)
//             },
//         );

//         Network {
//             graph: replacement.target.graph.clone(),
//             store: new_store,
//         }
//     }
// }

impl<T: Clone + PartialEq, S: TensorStructure + Clone + PartialEq> ReplaceWithBuilder
    for MixedTensor<T, S>
{
    type Ref<'a>
        = &'a Self
    where
        Self: 'a;
    type RefMut<'a>
        = &'a mut Self
    where
        Self: 'a;

    fn ref_mut(&mut self) -> Self::RefMut<'_> {
        self
    }

    fn recycled(ref_self: Self::Ref<'_>) -> Self {
        match ref_self {
            ParamOrConcrete::Param(a) => ParamOrConcrete::param(
                <DataTensor<Atom, S> as ReplaceWithBuilder>::recycled(a.tensor.reference()),
            ),
            _ => ref_self.clone(),
        }
    }
    fn reference(&self) -> Self::Ref<'_> {
        self
    }

    fn with<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
        with: C,
    ) -> Self {
        let with: BorrowedOrOwned<'c, Pattern> = with.into();

        let with_borrowed = with.borrow();
        match replacement.target {
            ParamOrConcrete::Param(a) => ParamOrConcrete::Param(a.map_data_ref_self(|a| {
                let rep =
                    replacement.update_symbolica_builder(a.replace(replacement.pattern.borrow()));
                rep.with(with_borrowed)
            })),
            _ => replacement.target.clone(),
        }
    }

    fn with_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
        replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
        rhs: M,
    ) -> Self {
        match replacement.target {
            ParamOrConcrete::Param(a) => ParamOrConcrete::Param(a.map_data_ref_self(|a| {
                let rep =
                    replacement.update_symbolica_builder(a.replace(replacement.pattern.borrow()));
                rep.with_map(dyn_clone::clone_box(&rhs))
            })),
            _ => replacement.target.clone(),
        }
    }

    fn with_into<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
        with: C,
        into: &mut Self,
    ) -> bool {
        let with: BorrowedOrOwned<'c, Pattern> = with.into();

        let with_borrowed = with.borrow();
        *into = match replacement.target {
            ParamOrConcrete::Param(a) => ParamOrConcrete::Param(a.map_data_ref_self(|a| {
                let rep =
                    replacement.update_symbolica_builder(a.replace(replacement.pattern.borrow()));
                rep.with(with_borrowed)
            })),
            _ => replacement.target.clone(),
        };
        replacement.target == into
    }

    fn with_mut<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: &'a mut ReplaceBuilderGeneric<'c, Self::RefMut<'a>, Self>,
        with: C,
    ) {
        let with: BorrowedOrOwned<'c, Pattern> = with.into();

        let with_borrowed = with.borrow();
        // let pattern = replacement.pattern.borrow().clone();
        let pattern_borrowed = replacement.pattern.borrow();
        if let ParamOrConcrete::Param(a) = &mut replacement.target {
            a.map_data_mut(|a| {
                let builder = a.replace(pattern_borrowed);

                let builder = if replacement.repeat {
                    builder.repeat()
                } else {
                    builder
                };

                let builder = if let Some(conditions) = &replacement.conditions {
                    builder.when(conditions.borrow())
                } else {
                    builder
                }
                .non_greedy_wildcards(replacement.settings.non_greedy_wildcards.clone())
                .allow_new_wildcards_on_rhs(replacement.settings.allow_new_wildcards_on_rhs)
                .level_range(replacement.settings.level_range)
                .level_is_tree_depth(replacement.settings.level_is_tree_depth)
                .rhs_cache_size(replacement.settings.rhs_cache_size);
                *a = builder.with(with_borrowed);
            })
        }
    }

    fn with_mut_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
        replacement: &'a mut ReplaceBuilderGeneric<'c, Self::RefMut<'a>, Self>,
        rhs: M,
    ) {
        if let ParamOrConcrete::Param(a) = &mut replacement.target {
            a.map_data_mut(|a| {
                let builder = a.replace(replacement.pattern.borrow());

                let builder = if replacement.repeat {
                    builder.repeat()
                } else {
                    builder
                };

                let builder = if let Some(conditions) = &replacement.conditions {
                    builder.when(conditions.borrow())
                } else {
                    builder
                }
                .non_greedy_wildcards(replacement.settings.non_greedy_wildcards.clone())
                .allow_new_wildcards_on_rhs(replacement.settings.allow_new_wildcards_on_rhs)
                .level_range(replacement.settings.level_range)
                .level_is_tree_depth(replacement.settings.level_is_tree_depth)
                .rhs_cache_size(replacement.settings.rhs_cache_size);
                *a = builder.with_map(dyn_clone::clone_box(&rhs));
            })
        }
    }

    fn with_owned<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: ReplaceBuilderGeneric<'c, Self, Self>,
        with: C,
    ) -> Self {
        let with: BorrowedOrOwned<'c, Pattern> = with.into();

        let with_borrowed = with.borrow();
        match replacement.target {
            ParamOrConcrete::Param(a) => ParamOrConcrete::Param(a.map_data_ref_self(|a| {
                let builder = a.replace(replacement.pattern.borrow());

                let builder = if replacement.repeat {
                    builder.repeat()
                } else {
                    builder
                };

                let builder = if let Some(conditions) = &replacement.conditions {
                    builder.when(conditions.borrow())
                } else {
                    builder
                }
                .non_greedy_wildcards(replacement.settings.non_greedy_wildcards.clone())
                .allow_new_wildcards_on_rhs(replacement.settings.allow_new_wildcards_on_rhs)
                .level_range(replacement.settings.level_range)
                .level_is_tree_depth(replacement.settings.level_is_tree_depth)
                .rhs_cache_size(replacement.settings.rhs_cache_size);
                builder.with(with_borrowed)
            })),
            a => a.clone(),
        }
    }

    fn with_owned_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
        replacement: ReplaceBuilderGeneric<'c, Self, Self>,
        rhs: M,
    ) -> Self {
        match replacement.target {
            ParamOrConcrete::Param(a) => ParamOrConcrete::Param(a.map_data_self(|a| {
                let builder = a.replace(replacement.pattern.borrow());

                let builder = if replacement.repeat {
                    builder.repeat()
                } else {
                    builder
                };

                let builder = if let Some(conditions) = &replacement.conditions {
                    builder.when(conditions.borrow())
                } else {
                    builder
                }
                .non_greedy_wildcards(replacement.settings.non_greedy_wildcards.clone())
                .allow_new_wildcards_on_rhs(replacement.settings.allow_new_wildcards_on_rhs)
                .level_range(replacement.settings.level_range)
                .level_is_tree_depth(replacement.settings.level_is_tree_depth)
                .rhs_cache_size(replacement.settings.rhs_cache_size);
                builder.with_map(dyn_clone::clone_box(&rhs))
            })),
            a => a,
        }
    }
}
impl<S: StorageTensor<Data = Atom> + Clone + PartialEq> ReplaceWithBuilder for S {
    type Ref<'a>
        = &'a Self
    where
        Self: 'a;
    type RefMut<'a>
        = &'a mut Self
    where
        Self: 'a;

    fn ref_mut(&mut self) -> Self::RefMut<'_> {
        self
    }

    fn recycled(ref_self: Self::Ref<'_>) -> Self {
        ref_self.map_data_ref_self(|_| RecycledAtom::new().into_inner())
    }
    fn reference(&self) -> Self::Ref<'_> {
        self
    }

    fn with<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
        with: C,
    ) -> Self {
        let with: BorrowedOrOwned<'c, Pattern> = with.into();

        let with_borrowed = with.borrow();
        replacement.target.map_data_ref_self(|a| {
            let rep = replacement.update_symbolica_builder(a.replace(replacement.pattern.borrow()));
            rep.with(with_borrowed)
        })
    }

    fn with_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
        replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
        rhs: M,
    ) -> Self {
        // let rhs = Box::new(rhs);
        // dyn_clone::clone_box(&rhs);
        replacement.target.map_data_ref_self(|a| {
            let rep = replacement.update_symbolica_builder(a.replace(replacement.pattern.borrow()));
            rep.with_map(dyn_clone::clone_box(&rhs))
        })
    }

    fn with_into<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: &'a ReplaceBuilderGeneric<'c, Self::Ref<'a>, Self>,
        with: C,
        into: &mut Self,
    ) -> bool {
        let with: BorrowedOrOwned<'c, Pattern> = with.into();

        let with_borrowed = with.borrow();
        *into = replacement.target.map_data_ref_self(move |expr| {
            let rep =
                replacement.update_symbolica_builder(expr.replace(replacement.pattern.borrow()));
            rep.with(with_borrowed)
        });
        replacement.target == into
    }

    fn with_mut<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: &'a mut ReplaceBuilderGeneric<'c, Self::RefMut<'a>, Self>,
        with: C,
    ) {
        let with: BorrowedOrOwned<'c, Pattern> = with.into();

        let with_borrowed = with.borrow();
        // let pattern = replacement.pattern.borrow().clone();
        let pattern_borrowed = replacement.pattern.borrow();
        replacement.target.map_data_mut(|a| {
            let builder = a.replace(pattern_borrowed);

            let builder = if replacement.repeat {
                builder.repeat()
            } else {
                builder
            };

            let builder = if let Some(conditions) = &replacement.conditions {
                builder.when(conditions.borrow())
            } else {
                builder
            }
            .non_greedy_wildcards(replacement.settings.non_greedy_wildcards.clone())
            .allow_new_wildcards_on_rhs(replacement.settings.allow_new_wildcards_on_rhs)
            .level_range(replacement.settings.level_range)
            .level_is_tree_depth(replacement.settings.level_is_tree_depth)
            .rhs_cache_size(replacement.settings.rhs_cache_size);
            *a = builder.with(with_borrowed);
        })
    }

    fn with_mut_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
        replacement: &'a mut ReplaceBuilderGeneric<'c, Self::RefMut<'a>, Self>,
        rhs: M,
    ) {
        replacement.target.map_data_mut(|a| {
            let builder = a.replace(replacement.pattern.borrow());

            let builder = if replacement.repeat {
                builder.repeat()
            } else {
                builder
            };

            let builder = if let Some(conditions) = &replacement.conditions {
                builder.when(conditions.borrow())
            } else {
                builder
            }
            .non_greedy_wildcards(replacement.settings.non_greedy_wildcards.clone())
            .allow_new_wildcards_on_rhs(replacement.settings.allow_new_wildcards_on_rhs)
            .level_range(replacement.settings.level_range)
            .level_is_tree_depth(replacement.settings.level_is_tree_depth)
            .rhs_cache_size(replacement.settings.rhs_cache_size);
            *a = builder.with_map(dyn_clone::clone_box(&rhs));
        })
    }

    fn with_owned<'a, 'c, C: Into<BorrowedOrOwned<'c, Pattern>>>(
        replacement: ReplaceBuilderGeneric<'c, Self, Self>,
        with: C,
    ) -> Self {
        let with: BorrowedOrOwned<'c, Pattern> = with.into();

        let with_borrowed = with.borrow();
        replacement.target.map_data_self(|a| {
            let builder = a.replace(replacement.pattern.borrow());

            let builder = if replacement.repeat {
                builder.repeat()
            } else {
                builder
            };

            let builder = if let Some(conditions) = &replacement.conditions {
                builder.when(conditions.borrow())
            } else {
                builder
            }
            .non_greedy_wildcards(replacement.settings.non_greedy_wildcards.clone())
            .allow_new_wildcards_on_rhs(replacement.settings.allow_new_wildcards_on_rhs)
            .level_range(replacement.settings.level_range)
            .level_is_tree_depth(replacement.settings.level_is_tree_depth)
            .rhs_cache_size(replacement.settings.rhs_cache_size);
            builder.with(with_borrowed)
        })
    }

    fn with_owned_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
        replacement: ReplaceBuilderGeneric<'c, Self, Self>,
        rhs: M,
    ) -> Self {
        replacement.target.map_data_self(|a| {
            let builder = a.replace(replacement.pattern.borrow());

            let builder = if replacement.repeat {
                builder.repeat()
            } else {
                builder
            };

            let builder = if let Some(conditions) = &replacement.conditions {
                builder.when(conditions.borrow())
            } else {
                builder
            }
            .non_greedy_wildcards(replacement.settings.non_greedy_wildcards.clone())
            .allow_new_wildcards_on_rhs(replacement.settings.allow_new_wildcards_on_rhs)
            .level_range(replacement.settings.level_range)
            .level_is_tree_depth(replacement.settings.level_is_tree_depth)
            .rhs_cache_size(replacement.settings.rhs_cache_size);
            builder.with_map(dyn_clone::clone_box(&rhs))
        })
    }

    // fn with_owned_map<'a, 'c: 'a, M: MatchMap + 'static + Clone>(
    //     replacement: ReplaceBuilderGeneric<'c, Self, Self>,
    //     rhs: M,
    // ) -> Self {

    // }

    // fn replace_mut<'a, 'b>(
    //     replacement: &'a mut ReplaceBuilderGeneric<'b, Self::RefMut<'a>, Self>,
    //     with: &ReplaceWith<'a>,
    // ) {
    //     replacement.target.map_data_mut(|expr_ref| {
    //         let rep = replacement
    //             .update_symbolica_builder(expr_ref.replace(replacement.pattern.borrow()));

    //         match with {
    //             ReplaceWith::Pattern(p) => {
    //                 *expr_ref = rep.with(p.borrow());
    //             }
    //             ReplaceWith::Map(m) => {
    //                 *expr_ref = rep.with_map(m.as_ref());
    //             }
    //         }
    //     });
    // }

    // fn replace_owned<'a, 'b>(
    //     replacement: ReplaceBuilderGeneric<'b, Self, Self>,
    //     with: &ReplaceWith<'a>,
    // ) -> Self {
    //     replacement.target.map_data_self(|expr_ref| {
    //         let rep = replacement
    //             .update_symbolica_builder(expr_ref.replace(replacement.pattern.borrow()));
    //         match with {
    //             ReplaceWith::Pattern(p) => rep.with(p.borrow()),
    //             ReplaceWith::Map(m) => rep.with_map(m.as_ref()),
    //         }
    //     })
    // }

    // fn replace<'a, 'b>(
    //     replacement: &'a ReplaceBuilderGeneric<'b, Self::Ref<'a>, Self>,
    //     with: &ReplaceWith<'a>,
    // ) -> Self {
    //     replacement.target.map_data_ref_self(|expr| {
    //         let rep = replacement.update_symbolica_builder(expr.replace(&replacement.pattern));
    //         rep.with(with)
    //     })
    // }

    // fn replace_into<'a, 'b>(
    //     replacement: &'a ReplaceBuilderGeneric<'b, Self::Ref<'a>, Self>,
    //     with: &ReplaceWith<'a>,
    //     into: &mut Self,
    // ) -> bool {
    //     let mut modified = false;
    //     *into = replacement.target.map_data_ref_self(|expr| {
    //         let rep = replacement.update_symbolica_builder(expr.replace(&replacement.pattern));
    //         let out = rep.with(with);
    //         if out != expr {
    //             modified = true;
    //         }
    //         out
    //     });
    //     modified
    // }
}

impl<S: StorageTensor<Data = Atom>> TensorAtomMaps for S {
    type AtomContainer = S;
    type ContainerData<T> = <S as StorageTensor>::ContainerData<T>;
    type Ref<'a>
        = &'a Self
    where
        Self: 'a;

    fn replace<'b, P: Into<BorrowedOrOwned<'b, Pattern>>>(
        &self,
        pattern: P,
    ) -> ReplaceBuilderGeneric<'b, Self::Ref<'_>, Self> {
        ReplaceBuilderGeneric::new(self, pattern)
    }

    /// Replace all occurrences of the patterns, where replacements are tested in the order that they are given.
    fn replace_multiple<T: BorrowReplacement>(&self, replacements: &[T]) -> Self {
        self.map_data_ref_self(|a| a.replace_multiple(replacements))
    }

    fn replace_multiple_repeat<T: BorrowReplacement>(&self, replacements: &[T]) -> Self {
        self.map_data_ref_self(|a| a.replace_multiple_repeat(replacements))
    }
    fn replace_multiple_mut<T: BorrowReplacement>(&mut self, replacements: &[T]) {
        self.map_data_mut(|a| a.replace_multiple_mut(replacements));
    }
    fn replace_multiple_repeat_mut<T: BorrowReplacement>(&mut self, replacements: &[T]) {
        self.map_data_mut(|a| a.replace_multiple_repeat_mut(replacements));
    }
    fn replace_map_mut<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(&mut self, m: &F) {
        self.map_data_mut(|a| a.replace_map_mut(m));
    }

    /// Collect terms involving the literal occurrence of `x`.
    fn coefficient<T: AtomCore>(&self, x: T) -> Self {
        let x = x.as_atom_view();
        self.map_data_ref_self(|a| a.coefficient(x))
    }

    /// Write the expression over a common denominator.
    fn together(&self) -> Self {
        self.map_data_ref_self(|a| a.together())
    }

    /// Write the expression as a sum of terms with minimal denominators.
    fn apart(&self, x: Symbol) -> Self {
        self.map_data_ref_self(|a| a.apart(x))
    }

    /// Cancel all common factors between numerators and denominators.
    /// Any non-canceling parts of the expression will not be rewritten.
    fn cancel(&self) -> Self {
        self.map_data_ref_self(|a| a.cancel())
    }

    /// Factor the expression over the rationals.
    fn factor(&self) -> Self {
        self.map_data_ref_self(|a| a.factor())
    }

    /// Collect numerical factors by removing the numerical content from additions.
    /// For example, `-2*x + 4*x^2 + 6*x^3` will be transformed into `-2*(x - 2*x^2 - 3*x^3)`.
    ///
    /// The first argument of the addition is normalized to a positive quantity.
    fn collect_num(&self) -> Self {
        self.map_data_ref_self(|a| a.collect_num())
    }

    /// Expand an expression. The function [AtomCore::expand_via_poly] may be faster.
    fn expand(&self) -> Self {
        self.map_data_ref_self(|a| a.expand())
    }

    /// Expand the expression by converting it to a polynomial, optionally
    /// only in the indeterminate `var`. The parameter `E` should be a numerical type
    /// that fits the largest exponent in the expanded expression. Often,
    /// `u8` or `u16` is sufficient.
    fn expand_via_poly<E: Exponent, T: AtomCore>(&self, var: Option<T>) -> Self {
        let var = var.as_ref().map(|v| v.as_atom_view());
        self.map_data_ref_self(|a| a.expand_via_poly::<E, AtomView>(var))
    }

    /// Expand an expression in the variable `var`. The function [AtomCore::expand_via_poly] may be faster.
    fn expand_in<T: AtomCore>(&self, var: T) -> Self {
        let var = var.as_atom_view();
        self.map_data_ref_self(|a| a.expand_in(var))
    }

    /// Expand an expression in the variable `var`.
    fn expand_in_symbol(&self, var: Symbol) -> Self {
        self.map_data_ref_self(|a| a.expand_in_symbol(var))
    }

    // /// Expand an expression, returning `true` iff the expression changed.
    // fn expand_into(&self, var: Option<AtomView>, out: &mut Atom) -> bool {
    //     self.as_atom_view().expand_into(var, out)
    // }

    /// Distribute numbers in the expression, for example:
    /// `2*(x+y)` -> `2*x+2*y`.
    fn expand_num(&self) -> Self {
        self.map_data_ref_self(|a| a.expand_num())
    }

    /// Take a derivative of the expression with respect to `x`.
    fn derivative(&self, x: Symbol) -> Self {
        self.map_data_ref_self(|a| a.derivative(x))
    }

    // /// Take a derivative of the expression with respect to `x` and
    // /// write the result in `out`.
    // /// Returns `true` if the derivative is non-zero.
    // fn derivative_into(&self, x: Symbol, out: &mut Atom) -> bool {
    //     self.as_atom_view().derivative_into(x, out)
    // }

    /// Series expand in `x` around `expansion_point` to depth `depth`.
    fn series<T: AtomCore>(
        &self,
        x: Symbol,
        expansion_point: T,
        depth: Rational,
        depth_is_absolute: bool,
    ) -> Result<Self::ContainerData<Series<AtomField>>, String> {
        let expansion_point = expansion_point.as_atom_view();
        self.map_data_ref_result(|a| a.series(x, expansion_point, depth.clone(), depth_is_absolute))
    }

    /// Find the root of a function in `x` numerically over the reals using Newton's method.
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
    ) -> Result<Self::ContainerData<N>, String> {
        let binding = x.into();
        let x = binding.borrow();
        self.map_data_ref_result(|a| {
            a.nsolve::<N, Indeterminate>(x.clone(), init.clone(), prec.clone(), max_iterations)
        })
    }

    /// Evaluate a (nested) expression a single time.
    /// For repeated evaluations, use [Self::evaluator()] and convert
    /// to an optimized version or generate a compiled version of your expression.
    ///
    /// All variables and all user functions in the expression must occur in the map.
    fn evaluate<A: AtomCore + KeyLookup, T: Real, F: Fn(&Rational) -> T + Copy>(
        &self,
        coeff_map: F,
        const_map: &HashMap<A, T>,
        function_map: &HashMap<Symbol, EvaluationFn<A, T>>,
    ) -> Result<Self::ContainerData<T>, String> {
        self.map_data_ref_result(|a| a.evaluate(coeff_map, const_map, function_map))
    }

    /// Check if the expression could be 0, using (potentially) numerical sampling with
    /// a given tolerance and number of iterations.
    fn zero_test(&self, iterations: usize, tolerance: f64) -> Self::ContainerData<ConditionResult> {
        self.map_data_ref(|a| a.zero_test(iterations, tolerance))
    }

    /// Set the coefficient ring to the multivariate rational polynomial with `vars` variables.
    fn set_coefficient_ring(&self, vars: &Arc<Vec<PolyVariable>>) -> Self {
        self.map_data_ref_self(|a| a.set_coefficient_ring(vars))
    }

    /// Convert all coefficients to floats with a given precision `decimal_prec``.
    /// The precision of floating point coefficients in the input will be truncated to `decimal_prec`.
    fn to_float(&self, decimal_prec: u32) -> Self {
        self.map_data_ref_self(|a| a.to_float(decimal_prec))
    }

    // /// Convert all coefficients to floats with a given precision `decimal_prec``.
    // /// The precision of floating point coefficients in the input will be truncated to `decimal_prec`.
    // fn coefficients_to_float_into(&self, decimal_prec: u32, out: &mut Atom) {
    //     self.as_atom_view()
    //         .coefficients_to_float_into(decimal_prec, out);
    // }

    /// Map all coefficients using a given function.
    fn map_coefficient<F: Fn(CoefficientView) -> Coefficient + Copy>(&self, f: F) -> Self {
        self.map_data_ref_self(|a| a.map_coefficient(f))
    }

    // /// Map all coefficients using a given function.
    // fn map_coefficient_into<F: Fn(CoefficientView) -> Coefficient + Copy>(
    //     &self,
    //     f: F,
    //     out: &mut Atom,
    // ) {
    //     self.as_atom_view().map_coefficient_into(f, out);
    // }

    /// Map all floating point and rational coefficients to the best rational approximation
    /// in the interval `[self*(1-relative_error),self*(1+relative_error)]`.
    fn rationalize(&self, relative_error: &Rational) -> Self {
        self.map_data_ref_self(|a| a.rationalize(relative_error))
    }

    /// Convert the atom to a polynomial, optionally in the variable ordering
    /// specified by `var_map`. If new variables are encountered, they are
    /// added to the variable map. Similarly, non-polynomial parts are automatically
    /// defined as a new independent variable in the polynomial.
    fn to_polynomial<R: EuclideanDomain + ConvertToRing, E: Exponent>(
        &self,
        field: &R,
        var_map: Option<Arc<Vec<PolyVariable>>>,
    ) -> Self::ContainerData<MultivariatePolynomial<R, E>> {
        self.map_data_ref(|a| a.to_polynomial(field, var_map.clone()))
    }

    /// Convert the atom to a polynomial in specific variables.
    /// All other parts will be collected into the coefficient, which
    /// is a general expression.
    ///
    /// This routine does not perform expansions.
    fn to_polynomial_in_vars<E: Exponent>(
        &self,
        var_map: &Arc<Vec<PolyVariable>>,
    ) -> Self::ContainerData<MultivariatePolynomial<AtomField, E>> {
        self.map_data_ref(|a| a.to_polynomial_in_vars(var_map))
    }

    /// Convert the atom to a rational polynomial, optionally in the variable ordering
    /// specified by `var_map`. If new variables are encountered, they are
    /// added to the variable map. Similarly, non-rational polynomial parts are automatically
    /// defined as a new independent variable in the rational polynomial.
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
        self.map_data_ref(|a| a.to_rational_polynomial(field, out_field, var_map.clone()))
    }

    /// Convert the atom to a rational polynomial with factorized denominators, optionally in the variable ordering
    /// specified by `var_map`. If new variables are encountered, they are
    /// added to the variable map. Similarly, non-rational polynomial parts are automatically
    /// defined as a new independent variable in the rational polynomial.
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
        self.map_data_ref(|a| {
            a.to_factorized_rational_polynomial(field, out_field, var_map.clone())
        })
    }

    // /// Construct a printer for the atom with special options.
    // fn printer<'a>(&'a self, opts: PrintOptions) -> Self::ContainerData<AtomPrinter<'a>> {
    //     self.map_data_ref(|a| a.printer(opts))
    // }

    /// Print the atom in a form that is unique and independent of any implementation details.
    ///
    /// Anti-symmetric functions are not supported.
    fn to_canonical_string(&self) -> Self::ContainerData<String> {
        self.map_data_ref(|a| a.to_canonical_string())
    }

    /// Map the function `f` over all terms.
    fn map_terms_single_core(&self, f: impl Fn(AtomView) -> Atom + Clone) -> Self {
        self.map_data_ref_self(|a| a.map_terms_single_core(f.clone()))
    }

    /// Map the function `f` over all terms, using parallel execution with `n_cores` cores.
    fn map_terms(
        &self,
        f: impl Fn(AtomView) -> Atom + Send + Sync + Clone,
        n_cores: usize,
    ) -> Self {
        self.map_data_ref_self(|a| a.map_terms(f.clone(), n_cores))
    }

    fn to_pattern(&self) -> Self::ContainerData<Pattern> {
        self.map_data_ref(|a| a.to_pattern())
    }

    fn replace_map<F: Fn(AtomView, &Context, &mut Settable<'_, Atom>)>(
        &self,
        m: &F,
    ) -> Self::AtomContainer {
        self.map_data_ref_self(|a| a.replace_map(m))
    }
    // fn replace_map<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(&self, m: &F) -> Self {
    //     self.map_data_ref_self(|a| a.replace_map(m))
    // }

    // /// Return an iterator that replaces the pattern in the target once.
    // fn replace_iter<'a>(
    //     &'a self,
    //     pattern: &'a Pattern,
    //     rhs: &'a PatternOrMap,
    //     conditions: Option<&'a Condition<PatternRestriction>>,
    //     settings: Option<&'a MatchSettings>,
    // ) -> Self::ContainerData<ReplaceIterator<'a, 'a>> {
    //     self.map_data_ref(|a| a.replace_iter(pattern, rhs, conditions, settings))
    // }
}

impl<S: TensorStructure> DenseTensor<Atom, S> {
    /// Solve a non-linear system numerically over the reals using Newton's method.
    pub fn nsolve_system<
        N: SingleFloat + Real + PartialOrd + InternalOrdering + Eq + std::hash::Hash,
    >(
        &self,
        vars: &[Indeterminate],
        init: &[N],
        prec: N,
        max_iterations: usize,
    ) -> Result<Vec<N>, String> {
        <Atom as AtomCore>::nsolve_system::<N, Atom>(&self.data, vars, init, prec, max_iterations)
    }

    /// Solve a system that is linear in `vars`, if possible.
    /// Each expression in `system` is understood to yield 0.
    pub fn solve_linear_system<E: PositiveExponent, T: AtomCore>(
        &self,
        vars: &[T],
    ) -> Result<Vec<Atom>, SolveError> {
        <Atom as AtomCore>::solve_linear_system::<E, Atom, T>(&self.data, vars)
    }
    /// Convert a system of linear equations to a matrix representation, returning the matrix
    /// and the right-hand side.
    #[allow(clippy::type_complexity)]
    pub fn system_to_matrix<E: PositiveExponent, T: AtomCore>(
        &self,
        vars: &[T],
    ) -> Result<
        (
            Matrix<RationalPolynomialField<Z, E>>,
            Matrix<RationalPolynomialField<Z, E>>,
        ),
        String,
    > {
        <Atom as AtomCore>::system_to_matrix::<E, Atom, T>(&self.data, vars)
    }
}

pub trait TensorAtomOps: HasStructure {
    /// Collect terms involving the same power of `x` in `xs`, where `xs` is a list of indeterminates.
    /// Return the list of key-coefficient pairs
    fn coefficient_list<E: Exponent, T: AtomCore>(&self, xs: &[T]) -> Vec<(Atom, Atom)>;

    /// Convert nested expressions to a tree suitable for repeated evaluations with
    /// different values for `params`.
    /// All variables and all user functions in the expression must occur in the map.
    fn to_evaluation_tree(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
    ) -> Result<EvalTreeTensor<SymComplex<Rational>, Self::Structure>, String>;

    /// Create an efficient evaluator for a (nested) expression.
    /// All free parameters must appear in `params` and all other variables
    /// and user functions in the expression must occur in the function map.
    /// The function map may have nested expressions.
    fn evaluator(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
        optimization_settings: &OptimizationSettings,
    ) -> Result<EvalTensor<ExpressionEvaluator<SymComplex<Rational>>, Self::Structure>, String>;

    /// Get all symbols in the expression, optionally including function symbols.
    fn get_all_symbols(&self, include_function_symbols: bool) -> HashSet<Symbol>;

    /// Get all variables and functions in the expression.
    fn get_all_indeterminates(&self, enter_functions: bool) -> HashSet<AtomView<'_>>;

    /// Returns true iff `self` contains the symbol `s`.
    fn contains_symbol(&self, s: Symbol) -> bool;

    /// Returns true iff `self` contains `a` literally.
    fn contains<T: AtomCore>(&self, s: T) -> bool;

    /// Check if the expression can be considered a polynomial in some variables, including
    /// redefinitions. For example `f(x)+y` is considered a polynomial in `f(x)` and `y`, whereas
    /// `f(x)+x` is not a polynomial.
    ///
    /// Rational powers or powers in variables are not rewritten, e.g. `x^(2y)` is not considered
    /// polynomial in `x^y`.
    fn is_polynomial(
        &self,
        allow_not_expanded: bool,
        allow_negative_powers: bool,
    ) -> Option<HashSet<AtomView<'_>>>;
}

impl<S: TensorStructure + Clone> TensorAtomOps for DenseTensor<Atom, S> {
    fn get_all_indeterminates(&self, enter_functions: bool) -> HashSet<AtomView<'_>> {
        let mut indeterminates = HashSet::new();

        for (_, a) in self.iter_flat() {
            indeterminates.extend(a.get_all_indeterminates(enter_functions));
        }

        indeterminates
    }

    fn contains_symbol(&self, s: Symbol) -> bool {
        self.iter_flat().any(|(_, a)| a.contains_symbol(s))
    }

    fn get_all_symbols(&self, include_function_symbols: bool) -> HashSet<Symbol> {
        let mut all_symbols = HashSet::new();
        for (_, a) in self.iter_flat() {
            all_symbols.extend(a.get_all_symbols(include_function_symbols));
        }
        all_symbols
    }

    fn coefficient_list<E: Exponent, T: AtomCore>(&self, xs: &[T]) -> Vec<(Atom, Atom)> {
        let mut list = vec![];

        for (_, a) in self.iter_flat() {
            list.extend(a.coefficient_list::<E>(xs))
        }
        list
    }

    fn contains<T: AtomCore>(&self, s: T) -> bool {
        let s = s.as_atom_view();
        self.iter_flat().any(|(_, a)| a.contains(s))
    }

    fn to_evaluation_tree(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
    ) -> Result<EvalTreeTensor<SymComplex<Rational>, Self::Structure>, String> {
        let atomviews: Vec<AtomView> = self.data.iter().map(|a| a.as_view()).collect();
        let eval = AtomView::to_eval_tree_multiple(&atomviews, fn_map, params)?;

        Ok(EvalTreeTensor {
            eval,
            indexmap: None,
            structure: self.structure.clone(),
        })
    }

    fn evaluator(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
        optimization_settings: &OptimizationSettings,
    ) -> Result<EvalTensor<ExpressionEvaluator<SymComplex<Rational>>, Self::Structure>, String>
    {
        let mut tree = self.to_evaluation_tree(fn_map, params)?;

        Ok(tree.optimize(optimization_settings))
    }

    fn is_polynomial(
        &self,
        allow_not_expanded: bool,
        allow_negative_powers: bool,
    ) -> Option<HashSet<AtomView<'_>>> {
        let mut is_polynomial = true;
        let mut set = HashSet::new();
        for (_, v) in self.iter_flat() {
            if let Some(v) = v.is_polynomial(allow_not_expanded, allow_negative_powers) {
                set.extend(v);
            } else {
                is_polynomial = false;
                break;
            }
        }

        if is_polynomial { Some(set) } else { None }
    }
}

impl<S: TensorStructure + Clone> TensorAtomOps for SparseTensor<Atom, S> {
    fn get_all_indeterminates(&self, enter_functions: bool) -> HashSet<AtomView<'_>> {
        let mut indeterminates = HashSet::new();

        for (_, a) in self.iter_flat() {
            indeterminates.extend(a.get_all_indeterminates(enter_functions));
        }

        indeterminates
    }

    fn contains_symbol(&self, s: Symbol) -> bool {
        self.iter_flat().any(|(_, a)| a.contains_symbol(s))
    }

    fn get_all_symbols(&self, include_function_symbols: bool) -> HashSet<Symbol> {
        let mut all_symbols = HashSet::new();
        for (_, a) in self.iter_flat() {
            all_symbols.extend(a.get_all_symbols(include_function_symbols));
        }
        all_symbols
    }

    fn coefficient_list<E: Exponent, T: AtomCore>(&self, xs: &[T]) -> Vec<(Atom, Atom)> {
        let mut list = vec![];

        for (_, a) in self.iter_flat() {
            list.extend(a.coefficient_list::<E>(xs))
        }
        list
    }

    fn contains<T: AtomCore>(&self, s: T) -> bool {
        let s = s.as_atom_view();
        self.iter_flat().any(|(_, a)| a.contains(s))
    }

    fn to_evaluation_tree(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
    ) -> Result<EvalTreeTensor<SymComplex<Rational>, Self::Structure>, String> {
        let atomviews: Vec<AtomView> = self.iter_flat().map(|(_, a)| a.as_view()).collect();
        let eval = AtomView::to_eval_tree_multiple(&atomviews, fn_map, params)?;

        Ok(EvalTreeTensor {
            eval,
            indexmap: None,
            structure: self.structure.clone(),
        })
    }

    fn evaluator(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
        optimization_settings: &OptimizationSettings,
    ) -> Result<EvalTensor<ExpressionEvaluator<SymComplex<Rational>>, Self::Structure>, String>
    {
        let mut tree = self.to_evaluation_tree(fn_map, params)?;

        Ok(tree.optimize(optimization_settings))
    }

    fn is_polynomial(
        &self,
        allow_not_expanded: bool,
        allow_negative_powers: bool,
    ) -> Option<HashSet<AtomView<'_>>> {
        let mut is_polynomial = true;
        let mut set = HashSet::new();
        for (_, v) in self.iter_flat() {
            if let Some(v) = v.is_polynomial(allow_not_expanded, allow_negative_powers) {
                set.extend(v);
            } else {
                is_polynomial = false;
                break;
            }
        }

        if is_polynomial { Some(set) } else { None }
    }
}

impl<S: TensorStructure + Clone> TensorAtomOps for DataTensor<Atom, S> {
    fn contains<T: AtomCore>(&self, s: T) -> bool {
        let s = s.as_atom_view();
        match self {
            DataTensor::Dense(d) => d.contains(s),
            DataTensor::Sparse(d) => d.contains(s),
        }
    }

    fn evaluator(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
        optimization_settings: &OptimizationSettings,
    ) -> Result<EvalTensor<ExpressionEvaluator<SymComplex<Rational>>, Self::Structure>, String>
    {
        match self {
            DataTensor::Dense(d) => d.evaluator(fn_map, params, optimization_settings),
            DataTensor::Sparse(s) => s.evaluator(fn_map, params, optimization_settings),
        }
    }

    fn is_polynomial(
        &self,
        allow_not_expanded: bool,
        allow_negative_powers: bool,
    ) -> Option<HashSet<AtomView<'_>>> {
        match self {
            DataTensor::Dense(d) => d.is_polynomial(allow_not_expanded, allow_negative_powers),
            DataTensor::Sparse(s) => s.is_polynomial(allow_not_expanded, allow_negative_powers),
        }
    }

    fn contains_symbol(&self, s: Symbol) -> bool {
        match self {
            DataTensor::Dense(d) => d.contains_symbol(s),
            DataTensor::Sparse(d) => d.contains_symbol(s),
        }
    }

    fn get_all_symbols(&self, include_function_symbols: bool) -> HashSet<Symbol> {
        match self {
            DataTensor::Dense(d) => d.get_all_symbols(include_function_symbols),
            DataTensor::Sparse(s) => s.get_all_symbols(include_function_symbols),
        }
    }

    fn coefficient_list<E: Exponent, T: AtomCore>(&self, xs: &[T]) -> Vec<(Atom, Atom)> {
        match self {
            DataTensor::Dense(d) => d.coefficient_list::<E, T>(xs),
            DataTensor::Sparse(s) => s.coefficient_list::<E, T>(xs),
        }
    }

    fn to_evaluation_tree(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
    ) -> Result<EvalTreeTensor<SymComplex<Rational>, Self::Structure>, String> {
        match self {
            DataTensor::Dense(d) => d.to_evaluation_tree(fn_map, params),
            DataTensor::Sparse(s) => s.to_evaluation_tree(fn_map, params),
        }
    }

    fn get_all_indeterminates(&self, enter_functions: bool) -> HashSet<AtomView<'_>> {
        match self {
            DataTensor::Dense(d) => d.get_all_indeterminates(enter_functions),
            DataTensor::Sparse(s) => s.get_all_indeterminates(enter_functions),
        }
    }
}

impl<S: TensorStructure + Clone> TensorAtomOps for ParamTensor<S> {
    fn contains<T: AtomCore>(&self, s: T) -> bool {
        self.tensor.contains(s)
    }

    fn evaluator(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
        optimization_settings: &OptimizationSettings,
    ) -> Result<EvalTensor<ExpressionEvaluator<SymComplex<Rational>>, Self::Structure>, String>
    {
        self.tensor.evaluator(fn_map, params, optimization_settings)
    }

    fn is_polynomial(
        &self,
        allow_not_expanded: bool,
        allow_negative_powers: bool,
    ) -> Option<HashSet<AtomView<'_>>> {
        self.tensor
            .is_polynomial(allow_not_expanded, allow_negative_powers)
    }

    fn contains_symbol(&self, s: Symbol) -> bool {
        self.tensor.contains_symbol(s)
    }

    fn get_all_symbols(&self, include_function_symbols: bool) -> HashSet<Symbol> {
        self.tensor.get_all_symbols(include_function_symbols)
    }

    fn coefficient_list<E: Exponent, T: AtomCore>(&self, xs: &[T]) -> Vec<(Atom, Atom)> {
        self.tensor.coefficient_list::<E, T>(xs)
    }

    fn to_evaluation_tree(
        &self,
        fn_map: &FunctionMap<SymComplex<Rational>>,
        params: &[Atom],
    ) -> Result<EvalTreeTensor<SymComplex<Rational>, Self::Structure>, String> {
        self.tensor.to_evaluation_tree(fn_map, params)
    }

    fn get_all_indeterminates(&self, enter_functions: bool) -> HashSet<AtomView<'_>> {
        self.tensor.get_all_indeterminates(enter_functions)
    }
}
