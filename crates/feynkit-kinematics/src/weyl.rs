//! Exact Weyl factorization for complex null momenta.

use symbolica::atom::Atom;

use crate::{FourMomentum, SymbolicKinematicsError};

impl FourMomentum<Atom> {
    /// The Weyl matrix `[[E+pz, px-i py], [px+i py, E-pz]]`.
    ///
    /// Its determinant is the mostly-minus invariant `p.p`. Components are
    /// exact Symbolica expressions, suitable for precision-generic evaluators.
    pub fn weyl_matrix(&self) -> [[Atom; 2]; 2] {
        let [energy, px, py, pz] = self.components();
        let i_py = Atom::i() * py;
        [[energy + pz, px - &i_py], [px + i_py, energy - pz]]
    }

    /// Factor a nonzero null momentum into angle and square Weyl components.
    ///
    /// Returns `[angle, square]`, where `weyl_matrix()[i][j]` equals
    /// `angle[i] * square[j]` on `p.p = 0`. The angle is the selected column;
    /// the square is the selected row divided by its pivot. This fixes the
    /// reciprocal rescaling freedom without introducing a square-root branch.
    ///
    /// The caller must select a nonzero entry `[row, column]` at evaluation
    /// time and supply a null momentum. Selection can use primal values while
    /// numerical evaluation retains all precision and Dual components. A chart
    /// is valid only where its pivot remains nonzero. This method rejects an
    /// out-of-range or identically zero pivot; it does not infer assumptions
    /// about symbolic components or test the on-shell relation.
    pub fn null_weyl_spinors(
        &self,
        pivot: [usize; 2],
    ) -> Result<[[Atom; 2]; 2], SymbolicKinematicsError> {
        let [row, column] = pivot;
        if row >= 2 || column >= 2 {
            return Err(SymbolicKinematicsError::InvalidWeylPivot(pivot));
        }
        let matrix = self.weyl_matrix();
        let denominator = &matrix[row][column];
        if denominator.is_zero() {
            return Err(SymbolicKinematicsError::ZeroWeylPivot(pivot));
        }
        Ok([
            [matrix[0][column].clone(), matrix[1][column].clone()],
            [&matrix[row][0] / denominator, &matrix[row][1] / denominator],
        ])
    }

    /// Construct the complex null momentum whose Weyl matrix is `angle ⊗ square`.
    ///
    /// Mixing the angle of one null momentum with the square of another gives
    /// a null vector orthogonal to both. No conjugation is applied: the two
    /// spinors are independent complex components. This is the convention
    /// `p = 1/2 <angle | gamma | square]` with `(E, px, py, pz)` ordering.
    pub fn from_weyl_spinors(angle: &[Atom; 2], square: &[Atom; 2]) -> Self {
        let diagonal = [&angle[0] * &square[0], &angle[1] * &square[1]];
        let off_diagonal = [&angle[0] * &square[1], &angle[1] * &square[0]];
        Self::from_args(
            (&diagonal[0] + &diagonal[1]) / 2,
            (&off_diagonal[0] + &off_diagonal[1]) / 2,
            Atom::i() * (&off_diagonal[0] - &off_diagonal[1]) / 2,
            (&diagonal[0] - &diagonal[1]) / 2,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::{atom::AtomCore, parse};

    fn assert_zero(value: Atom) {
        assert!(value.together().is_zero(), "{value}");
    }

    #[test]
    fn determinant_and_all_charts_preserve_the_null_relation() {
        let momentum = FourMomentum::from_args(parse!("e"), parse!("x"), parse!("y"), parse!("z"));
        let matrix = momentum.weyl_matrix();
        assert_zero(
            &matrix[0][0] * &matrix[1][1] - &matrix[0][1] * &matrix[1][0] - momentum.mass_squared(),
        );
        let on_shell = parse!("x^2+y^2+z^2");
        for (row, entries) in matrix.iter().enumerate() {
            for (column, pivot) in entries.iter().enumerate() {
                let [angle, square] = momentum.null_weyl_spinors([row, column]).unwrap();
                let reconstructed = FourMomentum::from_weyl_spinors(&angle, &square);
                for (actual, expected) in reconstructed
                    .components()
                    .into_iter()
                    .zip(momentum.components())
                {
                    let difference = ((actual - expected) * pivot).together().expand();
                    assert_zero(difference.replace(parse!("e^2")).with(on_shell.clone()));
                }
            }
        }
    }

    #[test]
    fn cross_spinors_form_the_transverse_null_pair() {
        let a = [parse!("a0"), parse!("a1")];
        let b = [parse!("b0"), parse!("b1")];
        let c = [parse!("c0"), parse!("c1")];
        let d = [parse!("d0"), parse!("d1")];
        let first = FourMomentum::from_weyl_spinors(&a, &b);
        let second = FourMomentum::from_weyl_spinors(&c, &d);
        let transverse = [
            FourMomentum::from_weyl_spinors(&a, &d),
            FourMomentum::from_weyl_spinors(&c, &b),
        ];
        for vector in [&first, &second, &transverse[0], &transverse[1]] {
            assert_zero(vector.mass_squared());
        }
        for vector in &transverse {
            assert_zero(vector.dot(&first));
            assert_zero(vector.dot(&second));
        }
        assert_zero(transverse[0].dot(&transverse[1]) + first.dot(&second));
    }

    #[test]
    fn complex_axis_charts_and_invalid_pivots() {
        // Each rank-one matrix with a single nonzero entry exercises a chart
        // that cannot be replaced by one universally nonzero light-cone entry.
        for row in 0..2 {
            for column in 0..2 {
                let mut angle = [Atom::num(0), Atom::num(0)];
                let mut square = [Atom::num(0), Atom::num(0)];
                angle[row] = Atom::i();
                square[column] = Atom::num(2);
                let momentum = FourMomentum::from_weyl_spinors(&angle, &square);
                let [a, b] = momentum.null_weyl_spinors([row, column]).unwrap();
                let reconstructed = FourMomentum::from_weyl_spinors(&a, &b);
                for (actual, expected) in reconstructed
                    .components()
                    .into_iter()
                    .zip(momentum.components())
                {
                    assert_zero(actual - expected);
                }
                assert!(matches!(
                    momentum.null_weyl_spinors([1 - row, column]),
                    Err(SymbolicKinematicsError::ZeroWeylPivot(_))
                ));
                assert!(matches!(
                    momentum.null_weyl_spinors([row, 2]),
                    Err(SymbolicKinematicsError::InvalidWeylPivot(_))
                ));
            }
        }
    }
}
