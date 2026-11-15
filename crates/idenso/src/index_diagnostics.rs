//! Borrowed diagnostics for contractions left in a symbolic expression.

use std::sync::LazyLock;

use spenso::structure::{
    representation::{LibraryRep, Lorentz, Minkowski},
    slot::{SlotMatch, SlotMatcher, SlotView},
};
use symbolica::atom::{AtomView, Symbol};

static LORENTZ_SYMBOLS: LazyLock<[Symbol; 2]> = LazyLock::new(|| {
    [
        LibraryRep::from(Minkowski {}).symbol(),
        LibraryRep::from(Lorentz {}).symbol(),
    ]
});

#[derive(Default)]
pub(crate) struct LorentzOccurrences<'a> {
    slots: Vec<SlotView<'a>>,
    repeated: bool,
}

impl<'a> LorentzOccurrences<'a> {
    pub(crate) fn has_repeated(expression: AtomView<'a>) -> bool {
        let representations = *LORENTZ_SYMBOLS;
        Self::inspect(expression, &mut SlotMatcher::default(), &representations).repeated
    }

    fn inspect(
        expression: AtomView<'a>,
        matcher: &mut SlotMatcher,
        representations: &[Symbol; 2],
    ) -> Self {
        match matcher.classify(expression) {
            SlotMatch::Explicit(slot) => {
                return Self {
                    slots: if representations.contains(&slot.representation().head())
                        && !slot.is_concrete_index()
                    {
                        vec![slot]
                    } else {
                        Vec::new()
                    },
                    repeated: false,
                };
            }
            // Dimensions and index payloads are metadata, not more tensor slots.
            SlotMatch::Opaque => return Self::default(),
            SlotMatch::Other => {}
        }
        match expression {
            AtomView::Add(sum) => Self::combine(sum.iter(), true, matcher, representations),
            AtomView::Mul(product) => {
                Self::combine(product.iter(), false, matcher, representations)
            }
            AtomView::Fun(function) => {
                Self::combine(function.iter(), false, matcher, representations)
            }
            AtomView::Pow(power) => {
                let mut base = Self::inspect(power.get_base(), matcher, representations);
                if i64::try_from(power.get_exp())
                    .is_ok_and(|exponent| !(-1..=1).contains(&exponent))
                {
                    // A squared tensor component repeats its free abstract port.
                    // Already internal dummies belong to each independent copy.
                    base.repeated |= !base.slots.is_empty();
                }
                // Powers do not create an implicit index scope. Even an inverse
                // or fractional power retaining a free port must be checked
                // against the surrounding term. Explicit scopes remain part
                // of the exact index payload compared below.
                // Tensor-valued exponents are invalid, but if one reaches this
                // diagnostic its surviving indices must not escape the check.
                let exponent = Self::inspect(power.get_exp(), matcher, representations);
                Self::combine_occurrences([base, exponent].into_iter(), false)
            }
            _ => Self::default(),
        }
    }

    fn combine(
        children: impl Iterator<Item = AtomView<'a>>,
        alternatives: bool,
        matcher: &mut SlotMatcher,
        representations: &[Symbol; 2],
    ) -> Self {
        Self::combine_occurrences(
            children.map(|child| Self::inspect(child, matcher, representations)),
            alternatives,
        )
    }

    fn combine_occurrences(children: impl Iterator<Item = Self>, alternatives: bool) -> Self {
        let mut result = Self::default();
        for child in children {
            if child.repeated {
                result.repeated = true;
                return result;
            }
            for slot in child.slots {
                if !alternatives
                    && result.slots.iter().any(|previous| {
                        previous.index() == slot.index()
                            && previous.dimension() == slot.dimension()
                            && previous.representation().head() == slot.representation().head()
                    })
                {
                    result.repeated = true;
                    return result;
                }
                // A sum contributes alternatives, not simultaneous occurrences.
                // Their union still detects a contraction with an outer factor.
                if !result.slots.iter().any(|previous| {
                    previous.index() == slot.index()
                        && previous.dimension() == slot.dimension()
                        && previous.representation().head() == slot.representation().head()
                }) {
                    result.slots.push(slot);
                }
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::{abstract_index::AbstractIndex, representation::RepName},
    };
    use symbolica::{
        atom::{Atom, AtomCore},
        function, symbol,
    };

    use super::*;
    use crate::{
        IndexTooling,
        dirac::{AGS, GammaSimplifySettings},
        representations::Bispinor,
        tensor::{AlgebraSettings, SymbolicTensor},
        test_support::test_initialize,
    };

    #[test]
    fn repeated_lorentz_indices_survive_unknown_tensors_and_separate_open_chains() {
        let r = test_initialize();
        let mu = r.mink_d.to_symbolic([Atom::num(101)]);
        let unknown = SPENSO_TAG.tensor_symbol("idenso::unresolved_lorentz_tensor");
        assert!(function!(unknown, &mu, &mu).has_repeated_lorentz_indices());
        let spin = Bispinor {}.new_rep(4);
        let ports: [_; 4] = std::array::from_fn(|i| spin.to_symbolic([Atom::num(i)]));
        let source = function!(AGS.gamma, &ports[0], &ports[1], &mu)
            * function!(AGS.gamma, &ports[2], &ports[3], &mu);
        let reduced = SymbolicTensor::infer(source)
            .unwrap()
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::canonical()),
                ..Default::default()
            })
            .unwrap();
        assert!(reduced.expression().has_repeated_lorentz_indices());
    }

    #[test]
    fn repeated_lorentz_indices_do_not_pair_independent_sum_alternatives() {
        let r = test_initialize();
        let mu = r.mink_d.to_symbolic([Atom::num(102)]);
        let first = SPENSO_TAG.tensor_symbol("idenso::unresolved_lorentz_first");
        let second = SPENSO_TAG.tensor_symbol("idenso::unresolved_lorentz_second");
        let first = function!(first, &mu);
        let second = function!(second, &mu);
        assert!(!(&first + &second).has_repeated_lorentz_indices());
        assert!((&first * (&first + &second)).has_repeated_lorentz_indices());
        assert!(first.pow(2).has_repeated_lorentz_indices());
        assert!(first.pow(-2).has_repeated_lorentz_indices());
        assert!((first.pow(-1) * &second).has_repeated_lorentz_indices());
        assert!((first.pow(Atom::num(1) / Atom::num(2)) * &second).has_repeated_lorentz_indices());
        assert!(first.pow(&second).has_repeated_lorentz_indices());
        assert!(
            Atom::num(2)
                .pow(&first * &second)
                .has_repeated_lorentz_indices()
        );
    }

    #[test]
    fn repeated_lorentz_indices_respect_typed_scoped_and_fixed_components() {
        let r = test_initialize();
        let unknown = SPENSO_TAG.tensor_symbol("idenso::unresolved_lorentz_metadata");
        let slot = |dimension: Atom, index: Atom| Minkowski {}.to_symbolic([dimension, index]);
        let scope = symbol!("idenso::unresolved_lorentz_scope");
        let index = AbstractIndex::Normal(103);
        let plain = r.mink_d.to_symbolic([Atom::from(index)]);
        let scoped = r.mink_d.to_symbolic([Atom::from(index.scoped(scope))]);
        assert!(!function!(unknown, &plain, &scoped).has_repeated_lorentz_indices());
        let dimension = Atom::var(symbol!("idenso::unresolved_lorentz_other_dimension"));
        assert!(
            !function!(unknown, &plain, slot(dimension, Atom::from(index)))
                .has_repeated_lorentz_indices()
        );
        let concrete = function!(spenso::structure::abstract_index::AIND_SYMBOLS.cind, 1);
        let concrete = r.mink_d.to_symbolic([concrete]);
        assert!(!function!(unknown, &concrete, &concrete).has_repeated_lorentz_indices());
        let spin = r.bis4.to_symbolic([Atom::from(index)]);
        assert!(!function!(unknown, &plain, &spin).has_repeated_lorentz_indices());
    }

    #[test]
    fn repeated_lorentz_indices_include_equal_and_complementary_variances() {
        use spenso::structure::abstract_index::AIND_SYMBOLS;

        test_initialize();
        let unknown = SPENSO_TAG.tensor_symbol("idenso::unresolved_lorentz_variance");
        let base = Lorentz {}.new_rep(4).to_symbolic([Atom::num(105)]);
        let upper = function!(AIND_SYMBOLS.uind, &base);
        let lower = function!(AIND_SYMBOLS.dind, &base);
        // Equal variances are malformed Einstein notation, but are still
        // residual repeated Lorentz labels and deserve the same diagnostic.
        for (first, second) in [
            (&base, &base),
            (&upper, &upper),
            (&lower, &lower),
            (&upper, &lower),
        ] {
            assert!(function!(unknown, first, second).has_repeated_lorentz_indices());
        }
        let upper_tensor = function!(unknown, &upper);
        let lower_tensor = function!(unknown, &lower);
        assert!(upper_tensor.pow(2).has_repeated_lorentz_indices());
        assert!(lower_tensor.pow(-2).has_repeated_lorentz_indices());
        assert!(!(&upper_tensor + &lower_tensor).has_repeated_lorentz_indices());
        let other = Lorentz {}.new_rep(4).to_symbolic([Atom::num(106)]);
        assert!(!function!(unknown, &upper, other).has_repeated_lorentz_indices());
        let other_representation = Minkowski {}.new_rep(4).to_symbolic([Atom::num(105)]);
        assert!(!function!(unknown, &upper, other_representation).has_repeated_lorentz_indices());
    }

    #[test]
    fn repeated_lorentz_indices_are_checked_after_same_chain_clifford_reduction() {
        let r = test_initialize();
        let mu = r.mink_d.to_symbolic([Atom::num(104)]);
        let spin = Bispinor {}.new_rep(4);
        let a = spin.to_symbolic([Atom::num(1)]);
        let b = spin.to_symbolic([Atom::num(2)]);
        let source = spenso::chain!(
            a,
            b,
            function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, &mu),
            function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, &mu)
        );
        assert!(source.has_repeated_lorentz_indices());
        let reduced = SymbolicTensor::infer(source)
            .unwrap()
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::canonical()),
                ..Default::default()
            })
            .unwrap();
        assert!(!reduced.expression().has_repeated_lorentz_indices());
    }
}
