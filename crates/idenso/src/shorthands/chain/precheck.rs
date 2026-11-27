//! Necessary conditions for the chain endpoint matcher. Accepted functions
//! still go through the rule, which owns placeholder restrictions and binding.
use spenso::network::library::symbolic::ETS;
use spenso::network::tags::SPENSO_TAG as T;
use spenso::structure::representation::{LibraryRep, RepName};
use spenso::structure::slot::SlotMatcher;
use symbolica::atom::representation::FunView;
use symbolica::atom::{AtomCore, AtomView, Symbol};

pub(super) fn rule_head(head: Symbol) -> bool {
    head != T.chain && head != T.trace && head != ETS.metric && !T.is_scalar_metadata(head)
}

/// At distinct argument positions, the matcher needs a start followed by its
/// dual of the same dimension. Permutable heads can also bind the reverse order.
fn has_endpoint_pair(
    function: FunView<'_>,
    representation: LibraryRep,
    slots: &mut SlotMatcher,
) -> bool {
    if function.get_nargs() < 2 {
        return false;
    }
    let head = function.get_symbol();
    let permutable = head.is_symmetric() || head.is_antisymmetric() || head.is_cyclesymmetric();
    let dual = representation.dual();
    let mut starts = Vec::new();
    let mut ends = Vec::new();
    for argument in function.iter() {
        let Some((rep, dimension)) = slots.port_representation(argument) else {
            continue;
        };
        if (rep == dual && starts.contains(&dimension))
            || (permutable && rep == representation && ends.contains(&dimension))
        {
            return true;
        }
        if rep == representation {
            starts.push(dimension);
        } else if permutable && rep == dual {
            ends.push(dimension);
        }
    }
    false
}

/// The endpoint rule has max_level(0) and can only match this function.
pub(super) fn rule_may_match(
    function: FunView<'_>,
    representation: LibraryRep,
    slots: &mut SlotMatcher,
) -> bool {
    rule_head(function.get_symbol()) && has_endpoint_pair(function, representation, slots)
}

/// Inner-product selection searches both operands at every level. Existing
/// chain heads are eligible in addition to heads accepted by the endpoint rule.
pub(super) fn selection_may_match(
    operand: AtomView<'_>,
    representation: LibraryRep,
    slots: &mut SlotMatcher,
) -> bool {
    let mut found = false;
    operand.visitor(&mut |node| {
        if !found && let AtomView::Fun(function) = node {
            let head = function.get_symbol();
            found = (head == T.chain || rule_head(head))
                && has_endpoint_pair(function, representation, slots);
        }
        !found
    });
    found
}

#[cfg(test)]
mod tests;
