//! Transport scalarized numerators with an already validated parent witness.

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;
use symbolica::id::Replacement;

use crate::Vakint;
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;

pub(super) fn route_to_parent(
    numerator: AtomView,
    loop_momenta: &[Atom],
    coordinates: Option<&[Atom]>,
) -> Atom {
    let Some(coordinates) = coordinates else {
        return numerator.to_owned();
    };
    // Input-to-canonical routing has already been applied once by VakintTerm.
    // Compose only the stored canonical-to-parent witness, simultaneously at
    // component level, so no dot(sum, sum) reaches the scalar lowerer. External
    // and scalar spectators are untouched.
    let component = Atom::var(vk_symbol!("rustred_parent_component_"));
    let components = loop_momenta
        .iter()
        .enumerate()
        .map(|(axis, source)| {
            Replacement::new(
                source.to_pattern(),
                function!(S.k, Atom::num(axis + 1), &component).to_pattern(),
            )
            .allow_new_wildcards_on_rhs(true)
        })
        .collect::<Vec<_>>();
    let routing = coordinates
        .iter()
        .enumerate()
        .map(|(axis, target)| {
            Replacement::new(
                function!(S.k, Atom::num(axis + 1), &component).to_pattern(),
                target.replace_multiple(&components).to_pattern(),
            )
        })
        .collect::<Vec<_>>();
    Vakint::convert_from_dot_notation(numerator)
        .replace_multiple(&routing)
        .expand()
}
