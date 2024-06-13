use color_eyre::Report;
use smallvec::SmallVec;
use std::ops::{Add, Index, IndexMut, Mul, Sub};
use symbolica::domains::float::{NumericalFloatComparison, NumericalFloatLike, Real};

use crate::utils::{FloatLike, F};

#[derive(Debug, Clone, Default)]
/// square symmetric matrix for use in the tropical sampling algorithm
pub struct SquareMatrix<T> {
    data: SmallVec<[T; 36]>, // this allows us to store the L matrix on the stack for 6 loops and below
    dim: usize,
}

impl<T> Index<(usize, usize)> for SquareMatrix<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[index.0 * self.dim + index.1]
    }
}

impl<T> IndexMut<(usize, usize)> for SquareMatrix<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.data[index.0 * self.dim + index.1]
    }
}

impl<'a, T: FloatLike> Mul<&SquareMatrix<F<T>>> for &'a SquareMatrix<F<T>> {
    type Output = SquareMatrix<F<T>>;

    fn mul(self, rhs: &SquareMatrix<F<T>>) -> Self::Output {
        let mut result = self.new_zeros(self.dim);

        for row in 0..self.dim {
            for col in 0..self.dim {
                for k in 0..self.dim {
                    result[(row, col)] += &self[(row, k)] * &rhs[(k, col)];
                }
            }
        }

        result
    }
}

impl<T: FloatLike> Mul<F<T>> for SquareMatrix<F<T>> {
    type Output = SquareMatrix<F<T>>;

    fn mul(self, rhs: F<T>) -> Self::Output {
        let mut result = self.new_zeros(self.dim);

        for row in 0..self.dim {
            for col in 0..self.dim {
                result[(row, col)] = &self[(row, col)] * &rhs;
            }
        }

        result
    }
}

impl<'a, T: FloatLike> Add<&SquareMatrix<F<T>>> for &'a SquareMatrix<F<T>> {
    type Output = SquareMatrix<F<T>>;

    fn add(self, rhs: &SquareMatrix<F<T>>) -> Self::Output {
        let mut result = self.new_zeros(self.dim);

        for row in 0..self.dim {
            for col in 0..self.dim {
                result[(row, col)] = &self[(row, col)] + &rhs[(row, col)];
            }
        }

        result
    }
}

impl<'a, T: FloatLike> Sub<&SquareMatrix<F<T>>> for &'a SquareMatrix<F<T>> {
    type Output = SquareMatrix<F<T>>;

    fn sub(self, rhs: &SquareMatrix<F<T>>) -> Self::Output {
        let mut result = self.new_zeros(self.dim);

        for row in 0..self.dim {
            for col in 0..self.dim {
                result[(row, col)] = &self[(row, col)] - &rhs[(row, col)];
            }
        }

        result
    }
}

impl SquareMatrix<F<f64>> {
    pub fn zeros(dim: usize) -> Self {
        let zero = F(0.0);
        Self {
            data: SmallVec::from_elem(zero, dim * dim),
            dim,
        }
    }
}

impl<T: FloatLike> SquareMatrix<F<T>> {
    #[must_use]
    pub fn new_zeros(&self, dim: usize) -> Self {
        Self::zeros_from_scalar(&self.data[0], dim)
    }

    #[must_use]
    pub fn zeros_from_scalar(scalar: &F<T>, dim: usize) -> Self {
        let zero = scalar.zero();
        Self::fill(zero, dim)
    }

    pub fn fill(scalar: F<T>, dim: usize) -> Self {
        Self {
            data: SmallVec::from_elem(scalar, dim * dim),
            dim,
        }
    }
    /// Performs operations on a matrix for tropical sampling
    /// # Errors
    /// Returns an error if the matrix is not invertible
    pub fn decompose_for_tropical(&self) -> Result<DecompositionResult<F<T>>, Report> {
        let zero = self.data[0].zero();
        let one = zero.one();
        // start cholesky decomposition

        let mut q: SquareMatrix<F<T>> = self.new_zeros(self.dim);

        for i in 0..self.dim {
            let mut diagonal_entry_squared = self[(i, i)].clone();
            for j in 0..i {
                diagonal_entry_squared -= q[(i, j)].square();
            }

            let diagonal_entry = diagonal_entry_squared.sqrt();
            q[(i, i)] = diagonal_entry.clone();

            for j in i + 1..self.dim {
                let mut entry = self[(i, j)].clone();
                for k in 0..i {
                    entry -= q[(i, k)].square();
                }
                q[(j, i)] = entry / &diagonal_entry;
            }
        }
        // end cholesky decomposition

        // compute the determinant of Q and store the invesrses of the diagonal elements
        let mut det_q = one.clone();
        let mut inverse_diagonal_entries = SmallVec::<[F<T>; 6]>::new();

        for i in 0..self.dim {
            let q_ii = q[(i, i)].clone();
            det_q *= &q_ii;
            inverse_diagonal_entries.push(q_ii.inv());
        }

        (0..self.dim).fold(one.clone(), |acc, i| acc * &q[(i, i)]);
        let determinant = det_q.square();

        if det_q.is_zero() {
            return Err(Report::msg("Determinant of Q is zero!"));
        }

        // the matrix N is defined through Q = D(I + N)
        let mut n_matrix: SquareMatrix<F<T>> = self.new_zeros(self.dim);
        for row in 1..self.dim {
            let inverse_diagonal_element = inverse_diagonal_entries[row].clone();
            for col in 0..row {
                n_matrix[(row, col)] = inverse_diagonal_element.clone() * &q[(row, col)];
            }
        }

        let max_non_zero_power_of_n = self.dim - 1;

        // this algorithm is unoptizmied, will optimize later if this works
        let mut powers_of_n = SmallVec::<[SquareMatrix<F<T>>; 5]>::new();
        powers_of_n.push(n_matrix);

        for _ in 1..max_non_zero_power_of_n {
            let last_power_of_n = powers_of_n
                .last()
                .unwrap_or_else(|| unreachable!("Never empty due to push before"));
            let first_power_of_n = powers_of_n
                .first()
                .unwrap_or_else(|| unreachable!("Never empty due to push before"));
            powers_of_n.push(last_power_of_n * first_power_of_n);
        }

        let n_sum =
            powers_of_n
                .iter()
                .enumerate()
                .fold(self.new_zeros(self.dim), |acc, (i, mat)| {
                    if i % 2 == 0 {
                        &acc - mat
                    } else {
                        &acc + mat
                    }
                });

        let mut inverse_q = n_sum;
        for row in 0..self.dim {
            inverse_q[(row, row)] += &one;
            for col in 0..self.dim {
                inverse_q[(row, col)] *= &inverse_diagonal_entries[col];
            }
        }

        let mut q_transposed_inverse = self.new_zeros(self.dim);
        for row in 0..self.dim {
            for col in 0..self.dim {
                q_transposed_inverse[(row, col)] = inverse_q[(col, row)].clone();
            }
        }

        let inverse = &q_transposed_inverse * &inverse_q;

        Ok(DecompositionResult {
            determinant,
            inverse,
            q_transposed_inverse,
        })
    }

    pub fn get_dim(&self) -> usize {
        self.dim
    }
}

#[derive(Debug)]
pub struct DecompositionResult<T> {
    pub determinant: T,
    pub inverse: SquareMatrix<T>,
    pub q_transposed_inverse: SquareMatrix<T>,
}

#[cfg(test)]
mod tests {
    use crate::utils::assert_approx_eq;
     //FIXME: add f128 back

    use super::*;
    const EPSILON: F<f64> = F(1e-12);

    #[test]
    fn test_decompose_for_tropical_2x2() {
        let mut test_matrix = SquareMatrix::zeros(2);

        test_matrix[(0, 0)] = F(2.0);
        test_matrix[(1, 1)] = F(4.0);
        test_matrix[(0, 1)] = F(1.0);
        test_matrix[(1, 0)] = F(1.0);

        let decomposition_result = test_matrix.decompose_for_tropical().unwrap();

        assert_approx_eq(&decomposition_result.determinant, &7.0.into(), &EPSILON);

        let identity = &test_matrix * &decomposition_result.inverse;

        assert_approx_eq(&identity[(0, 0)], &1.0.into(), &EPSILON);
        assert_approx_eq(&identity[(0, 1)], &0.0.into(), &EPSILON);
        assert_approx_eq(&identity[(1, 0)], &0.0.into(), &EPSILON);
        assert_approx_eq(&identity[(1, 1)], &1.0.into(), &EPSILON);
    }

    #[test]
    fn test_decompose_for_tropical_f128_2x2() {
        let mut test_matrix = SquareMatrix::zeros(2);

        test_matrix[(0, 0)] = 2.0.into();
        test_matrix[(1, 1)] = 4.0.into();
        test_matrix[(0, 1)] = 1.0.into();
        test_matrix[(1, 0)] = 1.0.into();

        let decomposition_result = test_matrix.decompose_for_tropical().unwrap();

        assert_approx_eq(&decomposition_result.determinant, &7.0.into(), &EPSILON);

        let identity = &test_matrix * &decomposition_result.inverse;

        assert_approx_eq(&identity[(0, 0)], &1.0.into(), &EPSILON);
        assert_approx_eq(&identity[(0, 1)], &0.0.into(), &EPSILON);
        assert_approx_eq(&identity[(1, 0)], &0.0.into(), &EPSILON);
        assert_approx_eq(&identity[(1, 1)], &1.0.into(), &EPSILON);
    }

    #[test]
    fn test_decompose_for_tropical_4x4() {
        let mut wilson_matrix = SquareMatrix::zeros(4);

        wilson_matrix[(0, 0)] = 5.0.into();
        wilson_matrix[(1, 1)] = 10.0.into();
        wilson_matrix[(2, 2)] = 10.0.into();
        wilson_matrix[(3, 3)] = 10.0.into();
        wilson_matrix[(0, 1)] = 7.0.into();
        wilson_matrix[(1, 0)] = 7.0.into();
        wilson_matrix[(0, 2)] = 6.0.into();
        wilson_matrix[(2, 0)] = 6.0.into();
        wilson_matrix[(0, 3)] = 5.0.into();
        wilson_matrix[(3, 0)] = 5.0.into();
        wilson_matrix[(1, 2)] = 8.0.into();
        wilson_matrix[(2, 1)] = 8.0.into();
        wilson_matrix[(1, 3)] = 7.0.into();
        wilson_matrix[(3, 1)] = 7.0.into();
        wilson_matrix[(2, 3)] = 9.0.into();
        wilson_matrix[(3, 2)] = 9.0.into();

        let decomposition_result = wilson_matrix.decompose_for_tropical().unwrap();

        let identity = &wilson_matrix * &decomposition_result.inverse;

        let one = 1.0.into();
        let zero = 0.0.into();
        for row in 0..4 {
            for col in 0..4 {
                assert_approx_eq(
                    &identity[(row, col)],
                    if row == col { &one } else { &zero },
                    &EPSILON,
                );
            }
        }

        let nalgebra_dmatrix = nalgebra::DMatrix::from_row_slice(
            4,
            4,
            &[
                5.0, 7.0, 6.0, 5.0, 7.0, 10.0, 8.0, 7.0, 6.0, 8.0, 10.0, 9.0, 5.0, 7.0, 9.0, 10.0,
            ],
        );

        let determinant = nalgebra_dmatrix.determinant();

        assert_approx_eq(
            &decomposition_result.determinant,
            &determinant.into(),
            &EPSILON,
        );
    }

    #[test]
    fn test_decompose_for_tropical_4x4_f128() {
        let mut wilson_matrix = SquareMatrix::zeros(4);

        wilson_matrix[(0, 0)] = 5.0.into();
        wilson_matrix[(1, 1)] = 10.0.into();
        wilson_matrix[(2, 2)] = 10.0.into();
        wilson_matrix[(3, 3)] = 10.0.into();
        wilson_matrix[(0, 1)] = 7.0.into();
        wilson_matrix[(1, 0)] = 7.0.into();
        wilson_matrix[(0, 2)] = 6.0.into();
        wilson_matrix[(2, 0)] = 6.0.into();
        wilson_matrix[(0, 3)] = 5.0.into();
        wilson_matrix[(3, 0)] = 5.0.into();
        wilson_matrix[(1, 2)] = 8.0.into();
        wilson_matrix[(2, 1)] = 8.0.into();
        wilson_matrix[(1, 3)] = 7.0.into();
        wilson_matrix[(3, 1)] = 7.0.into();
        wilson_matrix[(2, 3)] = 9.0.into();
        wilson_matrix[(3, 2)] = 9.0.into();

        let decomposition_result = wilson_matrix.decompose_for_tropical().unwrap();
        let identity = &wilson_matrix * &decomposition_result.inverse;

        let one = 1.0.into();
        let zero = 0.0.into();
        for row in 0..4 {
            for col in 0..4 {
                assert_approx_eq(
                    &identity[(row, col)],
                    if row == col { &one } else { &zero },
                    &EPSILON,
                );
            }
        }

        let nalgebra_dmatrix = nalgebra::DMatrix::from_row_slice(
            4,
            4,
            &[
                5.0, 7.0, 6.0, 5.0, 7.0, 10.0, 8.0, 7.0, 6.0, 8.0, 10.0, 9.0, 5.0, 7.0, 9.0, 10.0,
            ],
        );

        let determinant = nalgebra_dmatrix.determinant();

        assert_approx_eq(
            &decomposition_result.determinant,
            &determinant.into(),
            &EPSILON,
        );
    }
}
