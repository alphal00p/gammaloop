//! Loop-momentum basis change between a user routing and a matched topology.
//!
//! A topology match identifies every input propagator with a canonical one and
//! records whether its orientation is reversed. The input loop momenta are then
//! fixed by the linear system `q_in(e) = s(e) q_canonical(e)` over all matched
//! propagators. Solving it exactly is the only way to rewrite a numerator: an
//! input loop momentum need not be a propagator momentum, and one that is may
//! carry a sign. The solution must be consistent with every matched propagator
//! (up to the sign of a line, which only enters squared) and measure preserving
//! (`|det| = 1`).

use std::collections::BTreeMap;

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::domains::rational::Rational;

use crate::VakintError;

fn invalid(detail: impl std::fmt::Display) -> VakintError {
    VakintError::InvalidMomentumExpression(format!("loop momentum routing: {detail}"))
}

/// A linear combination of momentum heads such as `k(2)`, keyed by the head.
type LinearMomentum = BTreeMap<Atom, Rational>;

/// Parse an expanded linear combination `sum_i c_i f_i(n_i)` with rational `c_i`.
fn linear_momentum(momentum: AtomView) -> Result<LinearMomentum, VakintError> {
    let expanded = momentum.expand();
    let mut result = LinearMomentum::new();
    let mut add_term = |term: AtomView| -> Result<(), VakintError> {
        let (coefficient, head) = match term {
            AtomView::Fun(_) => (Rational::one(), term.to_owned()),
            AtomView::Mul(product) => {
                let mut coefficient = None;
                let mut head = None;
                for factor in product.iter() {
                    match factor {
                        AtomView::Num(_) if coefficient.is_none() => {
                            coefficient = Some(Rational::try_from(factor).map_err(|_| {
                                invalid(format!("non-rational coefficient in {momentum}"))
                            })?)
                        }
                        AtomView::Fun(_) if head.is_none() => head = Some(factor.to_owned()),
                        _ => return Err(invalid(format!("{momentum} is not linear"))),
                    }
                }
                match (coefficient, head) {
                    (Some(coefficient), Some(head)) => (coefficient, head),
                    _ => return Err(invalid(format!("{momentum} is not linear"))),
                }
            }
            _ => return Err(invalid(format!("{momentum} is not linear"))),
        };
        let entry = result.entry(head).or_insert_with(Rational::zero);
        *entry = &*entry + &coefficient;
        Ok(())
    };
    match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().try_for_each(&mut add_term)?,
        view if view.is_zero() => {}
        view => add_term(view)?,
    }
    result.retain(|_, coefficient| !coefficient.is_zero());
    Ok(result)
}

/// Express each input loop momentum head in the canonical loop momenta.
///
/// `equations` pairs every matched input propagator momentum with the
/// canonical propagator momentum it was identified with, multiplied by `-1`
/// for an edge whose written orientation is reversed. Returns
/// `(input head, canonical linear combination)`.
///
/// A propagator only enters through its square, so the change of variables
/// `k_in = M k_canonical` is correct whenever `|det M| = 1` and it maps every
/// input line momentum onto plus or minus its canonical image. The written
/// orientations are tried first; if they are inconsistent (a line written
/// against the direction its momentum flows), the signs of an independent set
/// of lines are searched as well. An input for which no sign choice works does
/// not route the matched graph and is rejected.
pub(crate) fn input_loop_momenta_in_canonical_basis(
    equations: &[(Atom, Atom)],
) -> Result<Vec<(Atom, Atom)>, VakintError> {
    let parsed = equations
        .iter()
        .map(|(input, canonical)| {
            Ok((
                linear_momentum(input.as_view())?,
                linear_momentum(canonical.as_view())?,
            ))
        })
        .collect::<Result<Vec<_>, VakintError>>()?;
    let heads = |side: fn(&(LinearMomentum, LinearMomentum)) -> &LinearMomentum| {
        parsed
            .iter()
            .flat_map(|row| side(row).keys().cloned())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>()
    };
    let unknowns = heads(|row| &row.0);
    let targets = heads(|row| &row.1);
    if unknowns.len() != targets.len() {
        return Err(invalid(format!(
            "the input routing uses {} loop momenta but the matched topology has {}",
            unknowns.len(),
            targets.len()
        )));
    }
    let dense = |momentum: &LinearMomentum, basis: &[Atom]| {
        basis
            .iter()
            .map(|head| momentum.get(head).cloned().unwrap_or_else(Rational::zero))
            .collect::<Vec<_>>()
    };
    let inputs = parsed
        .iter()
        .map(|(input, _)| dense(input, &unknowns))
        .collect::<Vec<_>>();
    let images = parsed
        .iter()
        .map(|(_, canonical)| dense(canonical, &targets))
        .collect::<Vec<_>>();
    let basis = independent_rows(&inputs).ok_or_else(|| {
        invalid("the propagator momenta do not determine every input loop momentum")
    })?;
    let lines_agree = |solution: &[Vec<Rational>]| {
        inputs.iter().zip(&images).all(|(input, image)| {
            let mapped = (0..targets.len())
                .map(|column| {
                    input
                        .iter()
                        .zip(solution)
                        .fold(Rational::zero(), |sum, (coefficient, row)| {
                            &sum + &(coefficient * &row[column])
                        })
                })
                .collect::<Vec<_>>();
            mapped == *image || mapped.iter().zip(image).all(|(a, b)| *a == -b.clone())
        })
    };
    // Sign mask 0 keeps the written orientations of the basis lines.
    for signs in 0..1usize << basis.len() {
        let matrix = basis.iter().map(|&row| inputs[row].clone()).collect();
        let rhs = basis
            .iter()
            .enumerate()
            .map(|(position, &row)| {
                images[row]
                    .iter()
                    .map(|entry| {
                        if signs & (1 << position) == 0 {
                            entry.clone()
                        } else {
                            -entry.clone()
                        }
                    })
                    .collect()
            })
            .collect();
        let solution =
            solve_square(matrix, rhs).expect("independent rows form an invertible matrix");
        if determinant(solution.clone()).abs().is_one() && lines_agree(&solution) {
            return Ok(unknowns
                .into_iter()
                .zip(solution)
                .map(|(head, coefficients)| {
                    let image = targets
                        .iter()
                        .zip(coefficients)
                        .filter(|(_, coefficient)| !coefficient.is_zero())
                        .fold(Atom::Zero, |image, (target, coefficient)| {
                            image + Atom::num(coefficient) * target
                        });
                    (head, image)
                })
                .collect());
        }
    }
    Err(invalid(
        "no unimodular change of loop momenta maps the input propagator momenta onto the \
         matched topology (momentum conservation is violated or the basis is not unimodular)",
    ))
}

/// Indices of the first rows that are linearly independent, if they span
/// the full column space.
fn independent_rows(rows: &[Vec<Rational>]) -> Option<Vec<usize>> {
    let columns = rows.first().map_or(0, Vec::len);
    let mut echelon: Vec<(usize, Vec<Rational>)> = Vec::new();
    let mut selected = Vec::new();
    for (index, row) in rows.iter().enumerate() {
        let mut reduced = row.clone();
        for (pivot, pivot_row) in &echelon {
            if reduced[*pivot].is_zero() {
                continue;
            }
            let factor = &reduced[*pivot] / &pivot_row[*pivot];
            subtract_multiple(&mut reduced, pivot_row, &factor);
        }
        if let Some(pivot) = (0..columns).find(|&column| !reduced[column].is_zero()) {
            echelon.push((pivot, reduced));
            selected.push(index);
            if selected.len() == columns {
                return Some(selected);
            }
        }
    }
    None
}

/// `target -= factor * source`, entry by entry.
fn subtract_multiple(target: &mut [Rational], source: &[Rational], factor: &Rational) {
    for (entry, subtrahend) in target.iter_mut().zip(source) {
        *entry = &*entry - &(factor * subtrahend);
    }
}

/// Solve `matrix * x = rhs` for an invertible square `matrix`.
fn solve_square(
    mut matrix: Vec<Vec<Rational>>,
    mut rhs: Vec<Vec<Rational>>,
) -> Option<Vec<Vec<Rational>>> {
    let n = matrix.len();
    for column in 0..n {
        let found = (column..n).find(|&row| !matrix[row][column].is_zero())?;
        matrix.swap(column, found);
        rhs.swap(column, found);
        let pivot = matrix[column][column].clone();
        for entry in matrix[column].iter_mut().chain(rhs[column].iter_mut()) {
            *entry = &*entry / &pivot;
        }
        for row in 0..n {
            if row == column || matrix[row][column].is_zero() {
                continue;
            }
            let factor = matrix[row][column].clone();
            let (pivot_row, pivot_rhs) = (matrix[column].clone(), rhs[column].clone());
            subtract_multiple(&mut matrix[row], &pivot_row, &factor);
            subtract_multiple(&mut rhs[row], &pivot_rhs, &factor);
        }
    }
    Some(rhs)
}

fn determinant(mut matrix: Vec<Vec<Rational>>) -> Rational {
    let n = matrix.len();
    let mut determinant = Rational::one();
    for column in 0..n {
        let Some(found) = (column..n).find(|&row| !matrix[row][column].is_zero()) else {
            return Rational::zero();
        };
        if found != column {
            matrix.swap(found, column);
            determinant = -determinant;
        }
        let pivot = matrix[column][column].clone();
        determinant = &determinant * &pivot;
        for row in column + 1..n {
            if matrix[row][column].is_zero() {
                continue;
            }
            let factor = &matrix[row][column] / &pivot;
            let pivot_row = matrix[column].clone();
            subtract_multiple(&mut matrix[row], &pivot_row, &factor);
        }
    }
    determinant
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;

    fn solve(equations: &[(&str, &str)]) -> Result<Vec<(String, Atom)>, VakintError> {
        let equations = equations
            .iter()
            .map(|(input, canonical)| (vk_parse!(*input).unwrap(), vk_parse!(*canonical).unwrap()))
            .collect::<Vec<_>>();
        Ok(input_loop_momenta_in_canonical_basis(&equations)?
            .into_iter()
            .map(|(head, image)| (head.to_canonical_string(), image))
            .collect())
    }

    fn assert_images(equations: &[(&str, &str)], expected: &[&str]) {
        let images = solve(equations).unwrap();
        assert_eq!(images.len(), expected.len());
        for ((head, image), expected) in images.iter().zip(expected) {
            assert_eq!(
                (image - vk_parse!(*expected).unwrap()).expand(),
                Atom::Zero,
                "{head} -> {image}, expected {expected}"
            );
        }
    }

    #[test]
    fn reversed_single_momentum_keeps_its_sign() {
        // Sunset: edge 2 written as -k(2) on the reversed edge.
        assert_images(
            &[
                ("k(1)", "k(1)"),
                ("-k(2)", "-k(2)"),
                ("k(1)+k(2)", "k(1)+k(2)"),
            ],
            &["k(1)", "k(2)"],
        );
        // The same input matched with edge 2 flipped against canonical k(2).
        assert_images(
            &[
                ("k(1)", "k(1)"),
                ("-k(2)", "-(k(2))"),
                ("k(1)+k(2)", "k(1)+k(2)"),
            ],
            &["k(1)", "k(2)"],
        );
    }

    #[test]
    fn loop_momenta_that_are_not_propagator_momenta_are_solved() {
        // props (k1+k2, k2, k1+2k2) identified with canonical (k1, k2, k1+k2).
        assert_images(
            &[
                ("k(1)+k(2)", "k(1)"),
                ("k(2)", "k(2)"),
                ("k(1)+2*k(2)", "k(1)+k(2)"),
            ],
            &["k(1)-k(2)", "k(2)"],
        );
    }

    #[test]
    fn a_line_written_against_its_momentum_flow_only_needs_its_square() {
        // All three sunset lines written 1 -> 2 although k(1)+k(2) flows back:
        // the propagators are unchanged, so k(1), k(2) keep their meaning.
        assert_images(
            &[
                ("k(1)", "k(1)"),
                ("k(2)", "k(2)"),
                ("k(1)+k(2)", "-k(1)-k(2)"),
            ],
            &["k(1)", "k(2)"],
        );
        // k(1), k(2), k(1)-k(2) is the sunset with one single-momentum line
        // reversed; the first sign choice that works reverses k(1).
        assert_images(
            &[
                ("k(1)", "k(1)"),
                ("k(2)", "k(2)"),
                ("k(1)-k(2)", "k(1)+k(2)"),
            ],
            &["-k(1)", "k(2)"],
        );
    }

    #[test]
    fn routings_of_another_graph_and_non_unimodular_bases_are_rejected() {
        // No choice of signs makes k(1)+2k(2) a sunset line.
        assert!(
            solve(&[
                ("k(1)", "k(1)"),
                ("k(2)", "k(2)"),
                ("k(1)+2*k(2)", "k(1)+k(2)")
            ])
            .is_err()
        );
        // (k1+k2, k1-k2) spans only an index-2 sublattice.
        assert!(
            solve(&[
                ("k(1)+k(2)", "k(1)"),
                ("k(1)-k(2)", "k(2)"),
                ("2*k(1)", "k(1)+k(2)")
            ])
            .is_err()
        );
        // Nonlinear momentum.
        assert!(solve(&[("k(1)^2", "k(1)")]).is_err());
    }
}
