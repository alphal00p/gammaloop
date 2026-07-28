use spenso::{
    network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};
use symbolica::atom::{Atom, AtomView};

use crate::tensor::{SymbolicNetExt, SymbolicNetParse};

pub mod chain;
pub mod metric;
pub mod schoonschip;

pub trait UndoShorthands {
    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom;

    fn undo_all<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom;
    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom;
    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom;
    fn undo_trace<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom;
}

impl UndoShorthands for Atom {
    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
        self.as_view().undo_schoonschip::<Aind>()
    }

    fn undo_all<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
        self.as_view().undo_all::<Aind>()
    }

    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
        self.as_view().undo_dots::<Aind>()
    }

    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
        self.as_view().undo_chain::<Aind>()
    }

    fn undo_trace<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
        self.as_view().undo_trace::<Aind>()
    }
}

impl<'a> UndoShorthands for AtomView<'a> {
    fn undo_all<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::full(),
                    trace: true,
                    chain: true,
                },
                ..Default::default()
            })
            .unwrap();

        net.simple_execute::<()>()
    }
    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::none(),
                    trace: false,
                    chain: true,
                },
                ..Default::default()
            })
            .unwrap();

        net.simple_execute::<()>()
    }

    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
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
            .unwrap();

        net.simple_execute::<()>()
    }

    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
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
            .unwrap();

        net.simple_execute::<()>()
    }

    fn undo_trace<Aind: AbsInd + DummyAind + ParseableAind>(&self) -> Atom {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::none(),
                    trace: true,
                    chain: false,
                },
                ..Default::default()
            })
            .unwrap();

        net.simple_execute::<()>()
    }
}

#[cfg(test)]
mod tests {
    use spenso::{
        chain, chain_factor, mink, shadowing::symbolica_utils::AtomCoreExt,
        structure::abstract_index::AbstractIndex, vector,
    };

    use crate::test_support::test_initialize;

    use super::*;

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

        insta::assert_snapshot!( expr.undo_schoonschip::<AbstractIndex>().to_bare_ordered_string(), @"chain(mink(4,i),mink(4,j),schoonschip_only_factor(in,out,mink(4,d_1000000)))*compact_p(mink(4,d_1000000))");
    }
}
