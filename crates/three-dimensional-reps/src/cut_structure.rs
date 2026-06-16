use itertools::Itertools;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::utils::{determinant_i32_is_nonzero, determinant_i32_signum};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ContourClosure {
    Above,
    Below,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Residue {
    pub basis: Vec<usize>,
    pub sigmas: Vec<i32>,
    pub sign: i32,
    pub spanning_tree_index: usize,
}

#[derive(Debug, Error)]
pub enum CutStructureError {
    #[error("empty loop-line signature list")]
    EmptySignatures,
    #[error("contour closure has length {actual}, expected {expected}")]
    InvalidContourClosure { actual: usize, expected: usize },
}

pub type Result<T> = std::result::Result<T, CutStructureError>;

pub fn ltd_residues(
    loop_line_signatures: &[Vec<i32>],
    contour_closure: &[ContourClosure],
) -> Result<Vec<Residue>> {
    if loop_line_signatures.is_empty() {
        return Err(CutStructureError::EmptySignatures);
    }
    let n_loops = loop_line_signatures[0].len();
    if contour_closure.len() != n_loops {
        return Err(CutStructureError::InvalidContourClosure {
            actual: contour_closure.len(),
            expected: n_loops,
        });
    }

    let spanning_trees = spanning_tree_generators(loop_line_signatures);
    let mut residues = Vec::new();
    for sigmas in (0..n_loops).map(|_| [1, -1]).multi_cartesian_product() {
        for (tree_index, tree) in spanning_trees.iter().enumerate() {
            residues.extend(tree.residues_per_tree(&sigmas, contour_closure, tree_index));
        }
    }
    residues.sort_by_key(|residue| residue.spanning_tree_index);
    Ok(residues)
}

fn spanning_tree_generators(loop_line_signatures: &[Vec<i32>]) -> Vec<SpanningTreeGenerator> {
    let n_loops = loop_line_signatures[0].len();
    (0..loop_line_signatures.len())
        .combinations(n_loops)
        .filter_map(|basis| {
            let matrix = basis
                .iter()
                .map(|edge_id| loop_line_signatures[*edge_id].clone())
                .collect::<Vec<_>>();
            if !determinant_i32_is_nonzero(&matrix) {
                None
            } else {
                Some(SpanningTreeGenerator { matrix, basis })
            }
        })
        .collect()
}

#[derive(Debug, Clone)]
struct SpanningTreeGenerator {
    matrix: Vec<Vec<i32>>,
    basis: Vec<usize>,
}

impl SpanningTreeGenerator {
    fn residues_per_tree(
        &self,
        sigmas: &[i32],
        contour_closure: &[ContourClosure],
        tree_index: usize,
    ) -> Vec<Residue> {
        let residue_generators = self.residue_generators();
        let mut residues = residue_generators
            .into_iter()
            .filter_map(|generator| {
                generator
                    .residue_sign(sigmas, contour_closure)
                    .map(|sign| Residue {
                        basis: self.basis.clone(),
                        sigmas: sigmas.to_vec(),
                        sign,
                        spanning_tree_index: tree_index,
                    })
            })
            .collect::<Vec<_>>();

        let mut index = 0;
        while index < residues.len() {
            let opposite =
                residues
                    .iter()
                    .enumerate()
                    .skip(index + 1)
                    .find_map(|(other_index, other)| {
                        (other.sign == -residues[index].sign
                            && other.basis == residues[index].basis
                            && other.sigmas == residues[index].sigmas)
                            .then_some(other_index)
                    });
            if let Some(other_index) = opposite {
                residues.remove(other_index);
                residues.remove(index);
            } else {
                index += 1;
            }
        }
        residues
    }

    fn residue_generators(&self) -> Vec<ResidueGenerator> {
        let n_loops = self.matrix.len();
        (0..n_loops)
            .permutations(n_loops)
            .filter_map(|permutation| {
                let permuted = permutation
                    .iter()
                    .map(|row_index| self.matrix[*row_index].clone())
                    .collect::<Vec<_>>();
                let allowed = (0..n_loops).all(|rank| {
                    let sub = leading_submatrix(&permuted, rank + 1);
                    determinant_i32_is_nonzero(&sub)
                });
                allowed.then_some(ResidueGenerator {
                    matrix: permuted,
                    permutation,
                })
            })
            .collect()
    }
}

#[derive(Debug, Clone)]
struct ResidueGenerator {
    matrix: Vec<Vec<i32>>,
    permutation: Vec<usize>,
}

impl ResidueGenerator {
    fn residue_sign(&self, sigmas: &[i32], contour_closure: &[ContourClosure]) -> Option<i32> {
        let contour_sign = if contour_closure
            .iter()
            .filter(|closure| matches!(closure, ContourClosure::Below))
            .count()
            % 2
            == 0
        {
            1
        } else {
            -1
        };
        let sigma_sign = sigmas.iter().product::<i32>();
        let det_sign = determinant_i32_signum(&self.matrix);
        let heaviside_sign = self.heaviside_sign(sigmas, contour_closure)?;
        Some(contour_sign * sigma_sign * det_sign * heaviside_sign)
    }

    fn heaviside_sign(&self, sigmas: &[i32], contour_closure: &[ContourClosure]) -> Option<i32> {
        let n_loops = self.matrix.len();
        for (rank, contour) in contour_closure.iter().enumerate().take(n_loops) {
            let sub = leading_submatrix(&self.matrix, rank + 1);
            let det_sub_sign = determinant_i32_signum(&sub);
            let mut values = vec![None; n_loops];

            if rank == 0 {
                let basis_index = self.permutation[rank];
                values[basis_index] = Some(det_sub_sign * sigmas[basis_index]);
            } else {
                for row_to_remove in 0..=rank {
                    let basis_index = self.permutation[row_to_remove];
                    let minor = leading_minor(&self.matrix, rank, row_to_remove);
                    let cofactor_sign = if (row_to_remove + rank) % 2 == 0 {
                        1
                    } else {
                        -1
                    };
                    let numerator_sign = cofactor_sign * determinant_i32_signum(&minor);
                    let value_sign = numerator_sign * det_sub_sign * sigmas[basis_index];
                    values[basis_index] = Some(value_sign);
                }
            }

            if matches!(contour, ContourClosure::Above) {
                for value in values.iter_mut().flatten() {
                    *value *= -1;
                }
            }

            let present = values.iter().flatten().copied().collect::<Vec<_>>();
            if present.is_empty() || present.iter().all(|value| *value <= 0) {
                return None;
            }
            if present.iter().all(|value| *value >= 0) {
                continue;
            }
            return None;
        }
        Some(1)
    }
}

fn leading_submatrix(matrix: &[Vec<i32>], size: usize) -> Vec<Vec<i32>> {
    (0..size)
        .map(|row| (0..size).map(|column| matrix[row][column]).collect())
        .collect()
}

fn leading_minor(matrix: &[Vec<i32>], rank: usize, row_to_remove: usize) -> Vec<Vec<i32>> {
    (0..=rank)
        .filter(|row| *row != row_to_remove)
        .map(|row| (0..rank).map(|column| matrix[row][column]).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_loop_residues_match_two_cut_signs() {
        let residues = ltd_residues(&[vec![1]], &[ContourClosure::Below]).unwrap();
        assert_eq!(residues.len(), 1);
        assert_eq!(residues[0].basis, vec![0]);
        assert_eq!(residues[0].sigmas, vec![1]);
    }

    #[test]
    fn two_loop_triangle_has_ltd_residues() {
        let residues = ltd_residues(
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
            &[ContourClosure::Below, ContourClosure::Below],
        )
        .unwrap();
        assert!(!residues.is_empty());
        assert!(residues.iter().all(|residue| residue.basis.len() == 2));
        assert!(residues.iter().all(|residue| residue.sign != 0));
    }
}
