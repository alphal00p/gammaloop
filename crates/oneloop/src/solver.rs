use symbolica::atom::{Atom, AtomCore};

use crate::error::OneLoopError;

pub struct RationalSolver;

impl RationalSolver {
    pub fn solve(&self, system: &[Atom], vars: &[Atom]) -> Result<Vec<Atom>, OneLoopError> {
        Atom::solve_linear_system::<u8, _, _>(system, vars)
            .map_err(|e| OneLoopError::SolverFailed(e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::RationalSolver;
    use symbolica::atom::Atom;
    use symbolica::symbol;

    #[test]
    fn solves_a_two_by_two_system() {
        crate::ensure_symbolica_license();
        let x = Atom::var(symbol!("oneloop::x"));
        let y = Atom::var(symbol!("oneloop::y"));
        // x + y - 3 = 0,  x - y - 1 = 0   ->   x = 2, y = 1
        let system = [&x + &y - Atom::num(3), &x - &y - Atom::num(1)];
        let sol = RationalSolver
            .solve(&system, &[x, y])
            .expect("a determined 2x2 system has a unique solution");
        assert_eq!(sol, vec![Atom::num(2), Atom::num(1)]);
    }
}
