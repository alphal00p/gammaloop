use symbolica::atom::Atom;

use crate::family::IntegralFamily;
use crate::masters::{MasterBasis, OneLoopMasters};
use crate::reduce::reduce;

pub fn amplitude(family: &IntegralFamily) -> Atom {
    let basis = OneLoopMasters;
    reduce(family)
        .terms
        .iter()
        .fold(Atom::Zero, |acc, (coeff, master)| {
            acc + coeff * basis.symbol(master)
        })
}
