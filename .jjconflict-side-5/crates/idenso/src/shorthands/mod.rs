use spenso::{
    network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};
use symbolica::atom::{Atom, AtomView};

use crate::{
    NetworkToolingError,
    tensor::{SymbolicNetExt, SymbolicNetParse},
};

pub mod chain;
pub mod metric;
pub mod schoonschip;

pub trait UndoShorthands {
    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;

    fn undo_all<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;
    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;
    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;
    fn undo_trace<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;
}

impl UndoShorthands for Atom {
    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_schoonschip::<Aind>()
    }

    fn undo_all<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_all::<Aind>()
    }

    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_dots::<Aind>()
    }

    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_chain::<Aind>()
    }

    fn undo_trace<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_trace::<Aind>()
    }
}

impl<'a> UndoShorthands for AtomView<'a> {
    fn undo_all<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::full(),
                    trace: true,
                    chain: true,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }
    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::none(),
                    trace: false,
                    chain: true,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }

    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: true,
                        expand_inside_chains: false,
                        expand_schoonship: false,
                    },
                    trace: false,
                    chain: false,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }

    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: false,
                        expand_inside_chains: true,
                        expand_schoonship: true,
                    },
                    trace: false,
                    chain: false,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }

    fn undo_trace<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::none(),
                    trace: true,
                    chain: false,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }
}

#[cfg(test)]
mod tests {
    use spenso::{
        chain, chain_factor, dot, mink, network::tags::SPENSO_TAG,
        structure::abstract_index::AbstractIndex, vector,
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder},
        symbol,
    };
    use symbolica_utils::AtomPrintExt;

    use crate::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};
    use crate::test_support::test_initialize;

    use super::*;

    #[test]
    fn undo_dots_positive_integer_powers_roundtrip() {
        let _ = test_initialize();
        let product = dot!(
            vector!(powered_dot_p, mink!(4)),
            vector!(powered_dot_q, mink!(4))
        );

        for exponent in [2, 3] {
            let expression = product.clone().pow(Atom::num(exponent));
            let expanded = expression.undo_dots::<AbstractIndex>().unwrap();

            assert_eq!(
                expanded
                    .schoonschip_net::<AbstractIndex>()
                    .unwrap()
                    .metric_shorthand_to_dot(),
                expression,
                "independent contractions must survive power {exponent}"
            );
        }
    }

    #[test]
    fn undo_dots_powers_preserve_factorized_sums() {
        let _ = test_initialize();
        let product = dot!(vector!(sum_p, mink!(4)), vector!(sum_q, mink!(4)));
        let coefficient = (Atom::var(symbol!("dot_coefficient")) + Atom::num(1)).pow(5);
        let expression = coefficient * (product + Atom::num(1)).pow(2);
        let expanded = expression.undo_dots::<AbstractIndex>().unwrap();

        let AtomView::Mul(factors) = expanded.as_view() else {
            panic!("the product of sums must stay factorized: {expanded}");
        };
        assert_eq!(
            factors
                .iter()
                .filter(|factor| matches!(factor, AtomView::Add(_)))
                .count(),
            2,
        );
        assert_eq!(
            expanded
                .schoonschip_net::<AbstractIndex>()
                .unwrap()
                .metric_shorthand_to_dot(),
            expression,
        );
    }

    #[test]
    fn undo_dots_powers_reserve_existing_indices() {
        let _ = test_initialize();
        let spectator = vector!(
            reserved_dot_r,
            mink!(4, AbstractIndex::new_dummy_at(1_000_000).to_atom())
        );
        let product = dot!(
            vector!(reserved_dot_p, mink!(4)),
            vector!(reserved_dot_q, mink!(4))
        );
        let expression = spectator * product.pow(2);
        let expanded = expression.undo_dots::<AbstractIndex>().unwrap();

        assert_eq!(
            expanded
                .schoonschip_net::<AbstractIndex>()
                .unwrap()
                .metric_shorthand_to_dot(),
            expression,
        );
    }

    #[test]
    fn undo_dots_powers_preserve_self_dual_external_ports() {
        let _ = test_initialize();
        let product = dot!(vector!(open_dot_p, mink!(4)), vector!(open_dot_q, mink!(4)));
        let external = vector!(open_dot_r, mink!(4, mu));
        let norm = dot!(vector!(open_dot_r, mink!(4)), vector!(open_dot_r, mink!(4)));

        for exponent in [2, 3] {
            // Keep the external port and lowered dot in the same powered parser group.
            let expression = spenso::bracket!(&external, &product).pow(exponent);
            let expanded = expression.undo_dots::<AbstractIndex>().unwrap();
            let expected =
                product.pow(exponent) * norm.pow(exponent / 2) * external.pow(exponent % 2);
            assert_eq!(
                expanded
                    .schoonschip_net::<AbstractIndex>()
                    .unwrap()
                    .metric_shorthand_to_dot(),
                expected,
            );
        }
    }

    #[test]
    fn undo_schoonschip_across_chain() {
        let _ = test_initialize();
        let expr = chain!(
            mink!(4, i),
            mink!(4, j),
            chain_factor!(
                schoonschip_only_factor,
                in,
                out,
                (vector!(compact_p, mink!(4)))
            )
        );

        insta::assert_snapshot!( expr.undo_schoonschip::<AbstractIndex>().unwrap().to_bare_ordered_string(), @"chain(mink(4,i),mink(4,j),schoonschip_only_factor(in,out,mink(4,d_1000000)))*compact_p(mink(4,d_1000000))");
    }

    #[test]
    fn network_tooling_malformed_shorthand_returns_a_parse_error() {
        let _ = test_initialize();
        let expression = FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(Atom::var(symbol!("malformed_dot_operand")))
            .finish();

        assert!(matches!(
            expression
                .undo_dots::<AbstractIndex>()
                .expect_err("malformed dot notation should return an error"),
            NetworkToolingError::Parse { reason } if reason.contains("Invalid dot function")
        ));
    }
}
