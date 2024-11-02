use core::panic;

use crate::{
    gammaloop_integrand::DefaultSample,
    graph::{BareGraph, EdgeType, Graph, LoopExtSignature, LoopMomentumBasis},
    momentum::{Energy, FourMomentum, Polarization, Signature, ThreeMomentum},
    numerator::{Evaluate, Evaluators, Numerator, NumeratorState},
    utils::{FloatLike, F},
    Settings,
};
use bincode::{Decode, Encode};
use itertools::Itertools;
use log::debug;
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use spenso::{
    arithmetic::ScalarMul, complex::Complex, data::DataTensor, structure::ScalarTensor,
    upgrading_arithmetic::FallibleAdd,
};
use symbolica::domains::float::{NumericalFloatLike, Real};

#[derive(Debug, Clone, Copy, PartialEq)]
enum ContourClosure {
    Above,
    Below,
}

struct CutStructureGenerator {
    loop_line_signatures: Vec<Signature>,
    n_loops: usize,
    n_loop_lines: usize,
}

impl CutStructureGenerator {
    fn new(loop_line_signatures: Vec<Signature>) -> Self {
        let n_loops = loop_line_signatures[0].len();
        let n_loop_lines = loop_line_signatures.len();

        Self {
            loop_line_signatures,
            n_loops,
            n_loop_lines,
        }
    }

    fn generate_structure(
        &self,
        contour_closure: &[ContourClosure],
        simplify: bool,
    ) -> Vec<Vec<f64>> {
        if self.n_loops != contour_closure.len() {
            panic!("number of loops and contour closures do not match")
        }

        let residue_elements = self.get_residues(contour_closure, simplify);
        let mut cut_structure = vec![];

        for residue_element in residue_elements.iter() {
            let mut cut_vec_iter = residue_element.sigmas.clone().unwrap().into_iter();
            let residue_basis = residue_element.basis.clone().unwrap();
            let mut new_element = vec![];
            for i in 0..self.n_loop_lines {
                if residue_basis.contains(&i) {
                    new_element.push(cut_vec_iter.next().unwrap());
                } else {
                    new_element.push(0.);
                }
            }

            cut_structure.push(new_element);
        }

        cut_structure
    }

    fn get_residues(&self, contour_closure: &[ContourClosure], simplify: bool) -> Vec<Residue> {
        if self.n_loops != contour_closure.len() {
            panic!("number of loops and contour closures do not match")
        }

        let mut residue_elements = vec![];

        let spanning_trees = self.get_spanning_tree_generators();

        for sigmas in vec![[1., -1.]; self.n_loops]
            .into_iter()
            .multi_cartesian_product()
        {
            for (tree_nr, spanning_tree) in spanning_trees.iter().enumerate() {
                let residues_per_tree =
                    spanning_tree.get_residues_per_tree(&sigmas, contour_closure, simplify);

                for residue in residues_per_tree.iter() {
                    let mut new_residue = residue.clone();
                    new_residue.tree_nr = Some(tree_nr);
                    residue_elements.push(new_residue);
                }
            }
        }

        residue_elements.sort_by(|a, b| a.tree_nr.cmp(&b.tree_nr));
        residue_elements
    }

    fn get_spanning_tree_generators(&self) -> Vec<SpanningTreeGenerator> {
        let mut spanning_tree_generators = vec![];

        let reference_signature_matrices = self
            .loop_line_signatures
            .clone()
            .into_iter()
            .combinations(self.n_loops);

        let bases = (0..self.n_loop_lines).combinations(self.n_loops);

        for (basis, reference_signature_matrix) in bases.zip(reference_signature_matrices) {
            let reference_dmatrix = DMatrix::from_vec(
                self.n_loops,
                self.n_loops,
                reference_signature_matrix
                    .iter()
                    .flat_map(|row| row.into_iter().map(|s| (s as i8) as f64).collect_vec())
                    .collect_vec(),
            );

            if reference_dmatrix.determinant() != 0. {
                spanning_tree_generators.push(SpanningTreeGenerator::new(
                    reference_signature_matrix,
                    basis,
                ))
            } else {
                continue;
            }
        }

        spanning_tree_generators
    }
}

struct SpanningTreeGenerator {
    n_loops: usize,
    reference_signature_matrix: Vec<Signature>,
    basis: Vec<usize>,
}

impl SpanningTreeGenerator {
    fn new(reference_signature_matrix: Vec<Signature>, basis: Vec<usize>) -> Self {
        Self {
            n_loops: reference_signature_matrix[0].len(),
            reference_signature_matrix,
            basis,
        }
    }

    fn get_residues_per_tree(
        &self,
        sigmas: &[f64],
        contour_closure: &[ContourClosure],
        simplify: bool,
    ) -> Vec<Residue> {
        let mut residues_per_tree = vec![];

        let residue_generators = self.get_residue_generators();

        for residue_generator in residue_generators.iter() {
            let residue = residue_generator.get_residue(sigmas, contour_closure, simplify);

            if simplify {
                if residue.sign != 0. {
                    residues_per_tree.push(residue);
                }
            } else {
                residues_per_tree.push(residue);
            }
        }

        let mut residues_per_tree_with_cancelling_removed = vec![];
        for residue_element in residues_per_tree.iter() {
            let cancelling_residue = Residue {
                sign: -residue_element.sign,
                heavisides_arguments: residue_element.heavisides_arguments.clone(),
                sigmas: residue_element.sigmas.clone(),
                basis: residue_element.basis.clone(),
                tree_nr: residue_element.tree_nr,
            };

            if !residues_per_tree.contains(&cancelling_residue) {
                residues_per_tree_with_cancelling_removed.push(residue_element.clone());
            }
        }

        for residue_per_tree in residues_per_tree_with_cancelling_removed.iter() {
            if residue_per_tree.heavisides_arguments.is_empty() {
                continue;
            } else {
                panic!("residue left in heaviside")
            }
        }

        for residue_per_tree in residues_per_tree_with_cancelling_removed.iter_mut() {
            residue_per_tree.sigmas = Some(sigmas.to_vec());
            residue_per_tree.basis = Some(self.basis.clone());
        }

        residues_per_tree_with_cancelling_removed
    }

    fn get_residue_generators(&self) -> Vec<ResidueGenerator> {
        let mut residue_generators = vec![];

        let permuted_signature_matrices = self
            .reference_signature_matrix
            .iter()
            .permutations(self.n_loops);

        let permutations = (0..self.n_loops).permutations(self.n_loops);

        for (permutated_signature_matrix, permutation) in
            permuted_signature_matrices.zip(permutations)
        {
            let mut allowed = true;
            for r in 0..self.n_loops {
                let sub_matrix = permutated_signature_matrix
                    .iter()
                    .take(r + 1)
                    .map(|row| {
                        row.into_iter()
                            .take(r + 1)
                            .map(|s| (s as i8) as f64)
                            .collect_vec()
                    })
                    .collect_vec();

                let sub_matrix = DMatrix::from_vec(r + 1, r + 1, sub_matrix.concat());

                allowed = sub_matrix.determinant() != 0.;
                if !allowed {
                    break;
                }
            }

            if allowed {
                let new_dmatrix = DMatrix::from_row_slice(
                    self.n_loops,
                    self.n_loops,
                    &permutated_signature_matrix
                        .iter()
                        .flat_map(|row| row.into_iter().map(|s| (s as i8) as f64).collect_vec())
                        .collect_vec(),
                );

                residue_generators.push(ResidueGenerator {
                    n_loops: self.n_loops,
                    signature_matrix: new_dmatrix,
                    permutation,
                })
            }
        }

        residue_generators
    }
}

#[derive(Debug)]
struct ResidueGenerator {
    n_loops: usize,
    signature_matrix: DMatrix<f64>,
    permutation: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
struct Residue {
    sign: f64,
    heavisides_arguments: Vec<Vec<Option<f64>>>,
    sigmas: Option<Vec<f64>>,
    basis: Option<Vec<usize>>,
    tree_nr: Option<usize>,
}

impl ResidueGenerator {
    fn get_residue(
        &self,
        sigmas: &[f64],
        contour_closure: &[ContourClosure],
        simplify: bool,
    ) -> Residue {
        let sign = self.get_sign(sigmas, contour_closure);
        let (heavisides, heaviside_sign) = self.get_heavisides(sigmas, contour_closure, simplify);
        if heaviside_sign {
            Residue {
                sign,
                heavisides_arguments: heavisides,
                sigmas: None,
                basis: None,
                tree_nr: None,
            }
        } else {
            Residue {
                sign: 0.,
                heavisides_arguments: vec![],
                sigmas: None,
                basis: None,
                tree_nr: None,
            }
        }
    }

    fn get_sign(&self, sigmas: &[f64], contour_closure: &[ContourClosure]) -> f64 {
        let countour_sign = (-1_f64).powi(
            contour_closure
                .iter()
                .filter(|c| **c == ContourClosure::Below)
                .count() as i32,
        );

        let cut_sign = sigmas.iter().product::<f64>();
        let det_sign = self.signature_matrix.determinant();

        countour_sign * cut_sign * det_sign
    }

    fn get_heavisides(
        &self,
        sigmas: &[f64],
        countour_closure: &[ContourClosure],
        simplify: bool,
    ) -> (Vec<Vec<Option<f64>>>, bool) {
        let mut heavisides = vec![];

        for (r, closure) in countour_closure.iter().enumerate() {
            let mut heaviside = Heaviside::new(self.n_loops);

            let mut sub_matrix = DMatrix::from_element(r + 1, r + 1, 0.);
            for i in 0..=r {
                for j in 0..=r {
                    sub_matrix[(i, j)] = self.signature_matrix[(i, j)];
                }
            }

            if r == 0 {
                let basis_index = self.permutation[r];

                let heaviside_value = 1. / sub_matrix.determinant() * sigmas[basis_index];
                heaviside.add(basis_index, heaviside_value);
            } else {
                for m in 0..=r {
                    let basis_index = self.permutation[m];
                    let subsub_matrix = sub_matrix.clone().remove_column(r).remove_row(m);

                    let mut heaviside_value = 1. / sub_matrix.determinant();
                    heaviside_value *= (-1_f64).powf((m + r) as f64) * subsub_matrix.determinant();
                    heaviside_value *= sigmas[basis_index];
                    heaviside.add(basis_index, heaviside_value);
                }
            }

            if *closure == ContourClosure::Above {
                heaviside.invert_sign();
            }

            if simplify {
                heaviside.simplify();
                if !heaviside.sign {
                    return (vec![], false);
                }
            }

            if !heaviside.arguments.is_empty() {
                heavisides.push(heaviside.arguments);
            }
        }

        (heavisides, true)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Heaviside {
    arguments: Vec<Option<f64>>,
    sign: bool,
}

impl Heaviside {
    fn new(n_loops: usize) -> Self {
        Self {
            arguments: vec![None; n_loops],
            sign: true,
        }
    }

    fn add(&mut self, index: usize, value: f64) {
        self.arguments[index] = Some(value);
    }

    fn invert_sign(&mut self) {
        for arg in self.arguments.iter_mut().flatten() {
            *arg = arg.neg();
        }
    }

    fn simplify(&mut self) {
        if self.arguments.iter().all(|v| match v {
            Some(v) => *v <= 0.,
            None => true,
        }) {
            self.arguments.clear();
            self.sign = false;
        } else if self.arguments.iter().all(|v| match v {
            Some(v) => *v >= 0.,
            None => true,
        }) {
            self.arguments.clear();
            self.sign = true;
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
struct LTDTerm {
    associated_lmb: Vec<(usize, f64)>,
    signature_of_lmb: Vec<LoopExtSignature>,
}

impl LTDTerm {
    #[allow(clippy::type_complexity)]
    fn evaluate<T: FloatLike>(
        &self,
        emr: (&[ThreeMomentum<F<T>>], Option<&[ThreeMomentum<F<T>>]>),
        external_moms: (&[FourMomentum<F<T>>], Option<&[FourMomentum<F<T>>]>),
        polarizations: (
            &[Polarization<Complex<F<T>>>],
            Option<&[Polarization<Complex<F<T>>>]>,
        ),
        graph: &BareGraph,
        num: &mut Numerator<Evaluators>,
        setting: &Settings,
    ) -> DataTensor<Complex<F<T>>> {
        // compute on shell energies of the momenta in associated_lmb

        let zero = emr.0.iter().next().unwrap().px.zero();
        let one = zero.one();
        let two = zero.from_i64(2);

        let edge_momenta_of_associated_lmb = self
            .associated_lmb
            .iter()
            .map(|(i, s)| {
                let mut momentum = emr.0[*i].clone().into_on_shell_four_momentum(
                    graph.edges[*i]
                        .particle
                        .mass
                        .value
                        .map(|m| F::<T>::from_ff64(m.re)),
                );

                momentum.temporal *= Energy::new(F::<T>::from_f64(*s));
                momentum
            })
            .collect_vec();

        let edge_momenta_of_associated_lmb_rot = if let Some(rot_emr) = emr.1 {
            self.associated_lmb
                .iter()
                .map(|(i, s)| {
                    let mut momentum = rot_emr[*i].clone().into_on_shell_four_momentum(
                        graph.edges[*i]
                            .particle
                            .mass
                            .value
                            .map(|m| F::<T>::from_ff64(m.re)),
                    );

                    momentum.temporal *= Energy::new(F::<T>::from_f64(*s));
                    momentum
                })
                .collect_vec()
        } else {
            edge_momenta_of_associated_lmb.clone()
        };

        // iterate over remaining propagators
        let mut inv_res = one.clone();
        let mut energy_product = one.clone();

        let mut ltd_emr = vec![];

        for (index, edge) in graph.edges.iter().enumerate()
        // .filter(|(index, e)| {
        // e.edge_type == EdgeType::Virtual && self.associated_lmb.iter().all(|(i, _)| i != index)
        // })
        {
            if let Some(i) = self.associated_lmb.iter().position(|(i, _)| i == &index) {
                if setting.stability.rotate_numerator {
                    ltd_emr.push(edge_momenta_of_associated_lmb_rot[i].clone());
                } else {
                    ltd_emr.push(edge_momenta_of_associated_lmb[i].clone());
                }
            } else {
                match edge.edge_type {
                    EdgeType::Virtual => {
                        let momentum = self.signature_of_lmb[index].compute_momentum(
                            &edge_momenta_of_associated_lmb_rot,
                            external_moms.1.unwrap_or(external_moms.0),
                        );

                        if setting.stability.rotate_numerator {
                            ltd_emr.push(momentum.clone());
                        } else {
                            ltd_emr.push(self.signature_of_lmb[index].compute_momentum(
                                &edge_momenta_of_associated_lmb,
                                external_moms.0,
                            ));
                        }

                        match edge.particle.mass.value {
                            Some(mass) => {
                                inv_res *=
                                    momentum.clone().square() - F::<T>::from_ff64(mass.re).square();
                                energy_product *= (momentum.spatial.norm_squared()
                                    + F::<T>::from_ff64(mass.re).square())
                                .sqrt()
                                    * &two;
                            }
                            None => {
                                inv_res *= momentum.clone().square();
                                energy_product *= momentum.spatial.norm() * &two;
                            }
                        }
                    }
                    _ => {
                        if setting.stability.rotate_numerator {
                            ltd_emr.push(self.signature_of_lmb[index].compute_momentum(
                                &edge_momenta_of_associated_lmb_rot,
                                external_moms.1.unwrap_or(external_moms.0),
                            ));
                        } else {
                            ltd_emr.push(self.signature_of_lmb[index].compute_momentum(
                                &edge_momenta_of_associated_lmb,
                                external_moms.0,
                            ));
                        }
                    }
                }
            }
        }

        let polarizations = if setting.stability.rotate_numerator {
            polarizations.1.unwrap_or(polarizations.0)
        } else {
            polarizations.0
        };

        let num = num.evaluate_single(&ltd_emr, polarizations, None, setting);

        num.scalar_mul(&(inv_res.inv() * energy_product)).unwrap()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct LTDExpression {
    terms: Vec<LTDTerm>,
}

impl LTDExpression {
    pub fn evaluate<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        graph: &BareGraph,
        num: &mut Numerator<Evaluators>,
        setting: &Settings,
    ) -> DataTensor<Complex<F<T>>> {
        let zero = sample.zero();
        let possibly_rotated_emr = sample
            .rotated_sample
            .as_ref()
            .map(|s| graph.compute_emr(&s.loop_moms, &s.external_moms));
        let unrotated_emr =
            graph.compute_emr(&sample.sample.loop_moms, &sample.sample.external_moms);

        let external_moms = sample.external_mom_pair();

        self.terms
            .iter()
            .map(|term| {
                term.evaluate(
                    (&unrotated_emr, possibly_rotated_emr.as_deref()),
                    external_moms,
                    sample.polarizations_pair(),
                    graph,
                    num,
                    setting,
                )
            })
            .reduce(|acc, e| acc.add_fallible(&e).unwrap())
            .unwrap_or(DataTensor::new_scalar(Complex::new_re(zero)))
    }

    pub fn evaluate_in_lmb<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        graph: &BareGraph,
        lmb: &LoopMomentumBasis,
        num: &mut Numerator<Evaluators>,
        setting: &Settings,
    ) -> DataTensor<Complex<F<T>>> {
        let zero = sample.zero();
        let possibly_rotated_emr = sample
            .rotated_sample
            .as_ref()
            .map(|s| graph.compute_emr_in_lmb(&s.loop_moms, &s.external_moms, lmb));
        let unrotated_emr =
            graph.compute_emr_in_lmb(&sample.sample.loop_moms, &sample.sample.external_moms, lmb);
        let external_moms = sample.external_mom_pair();
        self.terms
            .iter()
            .map(|term| {
                term.evaluate(
                    (&unrotated_emr, possibly_rotated_emr.as_deref()),
                    external_moms,
                    sample.polarizations_pair(),
                    graph,
                    num,
                    setting,
                )
            })
            .reduce(|acc, e| acc.add_fallible(&e).unwrap())
            .unwrap_or(DataTensor::new_scalar(Complex::new_re(zero)))
    }
}

pub fn generate_ltd_expression<S: NumeratorState>(graph: &mut Graph<S>) -> LTDExpression {
    debug!(
        "generating ltd expression for graph: {:?}",
        graph.bare_graph.name
    );

    let loop_line_signatures = graph
        .bare_graph
        .get_virtual_edges_iterator()
        .map(|(index, _e)| {
            graph.bare_graph.loop_momentum_basis.edge_signatures[index]
                .internal
                .clone()
        })
        .collect_vec();

    let loop_number = loop_line_signatures[0].len();

    let position_map = graph
        .bare_graph
        .get_virtual_edges_iterator()
        .map(|(index, _e)| index)
        .collect_vec();

    let cut_structure_generator = CutStructureGenerator::new(loop_line_signatures);
    let countour_closure =
        vec![ContourClosure::Above; graph.bare_graph.loop_momentum_basis.basis.len()];
    let cut_structure = cut_structure_generator.generate_structure(&countour_closure, true);

    graph.generate_loop_momentum_bases_if_not_exists();
    debug!(
        "number of spanning trees: {}",
        graph
            .derived_data
            .as_ref()
            .unwrap()
            .loop_momentum_bases
            .as_ref()
            .unwrap()
            .len()
    );

    let mut ltd_terms = Vec::with_capacity(cut_structure.len());

    for cut_signature in cut_structure.iter() {
        let mut associated_lmb = Vec::with_capacity(loop_number);
        for (index, sign) in cut_signature.iter().enumerate() {
            if *sign != 0. {
                associated_lmb.push((position_map[index], *sign));
            }
        }
        // associated_lmb.sort_by(|a, b| a.0.cmp(&b.0));

        let mut found_signature = false;

        for loop_momentum_basis in graph
            .derived_data
            .as_ref()
            .unwrap()
            .loop_momentum_bases
            .as_ref()
            .unwrap()
            .iter()
        {
            if loop_momentum_basis
                .basis
                .iter()
                .zip(associated_lmb.iter())
                .all(|(a, b)| *a == b.0)
            {
                ltd_terms.push(LTDTerm {
                    associated_lmb: associated_lmb.clone(),
                    signature_of_lmb: loop_momentum_basis.edge_signatures.clone(),
                });
                found_signature = true;
                break;
            }
        }

        if !found_signature {
            panic!(
                "cut structure has no equivalent in the lmb: cut_structure: {:?}. associated_lmb: {:?}. all_lmbs: {:?}",
                cut_signature,
                associated_lmb  ,
                graph.derived_data.as_ref().unwrap().loop_momentum_bases.as_ref().unwrap(),
            )
        }
    }

    LTDExpression { terms: ltd_terms }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_heavisides() {
        let mut test_heaviside = Heaviside::new(3);
        test_heaviside.simplify();
        assert_eq!(test_heaviside.arguments, vec![]);
        assert!(!test_heaviside.sign);

        test_heaviside = Heaviside::new(3);
        test_heaviside.add(0, 1.);
        test_heaviside.add(1, 1.);
        test_heaviside.add(2, 1.);
        test_heaviside.simplify();
        assert_eq!(test_heaviside.arguments, vec![]);
        assert!(test_heaviside.sign);

        test_heaviside = Heaviside::new(3);
        test_heaviside.add(0, 1.);
        test_heaviside.add(1, 1.);
        test_heaviside.add(2, 1.);
        test_heaviside.invert_sign();
        test_heaviside.simplify();
        assert_eq!(test_heaviside.arguments, vec![]);
        assert!(!test_heaviside.sign);

        test_heaviside = Heaviside::new(3);
        test_heaviside.add(0, 1.);
        test_heaviside.add(1, 1.);
        test_heaviside.add(2, -1.);
        test_heaviside.simplify();
        assert_ne!(test_heaviside.arguments, vec![]);
    }

    #[test]
    fn test_residue_generator() {
        // test sign

        let test_permutation = vec![0, 1, 2];
        let test_signature_matrix = [vec![-1, 0, 1], vec![0, 1, 0], vec![1, 0, 0]];
        let test_sigmas = vec![1., -1., 1.];

        let test_residue_generator = ResidueGenerator {
            n_loops: 3,
            signature_matrix: DMatrix::from_row_slice(
                3,
                3,
                &test_signature_matrix
                    .iter()
                    .flat_map(|row| row.iter().map(|s| *s as f64))
                    .collect_vec(),
            ),
            permutation: test_permutation,
        };

        let test_contour_closure = vec![
            ContourClosure::Below,
            ContourClosure::Below,
            ContourClosure::Above,
        ];

        assert_eq!(
            test_residue_generator.get_sign(&test_sigmas, &test_contour_closure),
            1.
        );

        let test_contour_closure = vec![
            ContourClosure::Below,
            ContourClosure::Above,
            ContourClosure::Above,
        ];

        assert_eq!(
            test_residue_generator.get_sign(&test_sigmas, &test_contour_closure),
            -1.
        );

        let residue = test_residue_generator.get_residue(&test_sigmas, &test_contour_closure, true);
        assert_eq!(residue.sign, 0.0);
        for arg in residue.heavisides_arguments.iter() {
            assert_eq!(*arg, vec![]);
        }
    }

    #[test]
    fn test_spanning_tree_generator() {
        let contour_closure = [
            ContourClosure::Below,
            ContourClosure::Below,
            ContourClosure::Above,
        ];

        let basis = vec![0, 1, 2];

        let test_sigmas_1 = [1., -1., 1.];
        let test_sigmas_2 = [-1., 1., -1.];

        let reference_signature_matrix = vec![
            vec![-1, 0, 1].into(),
            vec![0, 1, 0].into(),
            vec![1, 0, 0].into(),
        ];
        let spanning_tree_generator = SpanningTreeGenerator::new(reference_signature_matrix, basis);

        let residue_generators = spanning_tree_generator.get_residue_generators();
        assert_eq!(residue_generators.len(), 2);

        assert_eq!(residue_generators[0].permutation, vec![0, 1, 2]);
        assert_eq!(residue_generators[1].permutation, vec![2, 1, 0]);

        assert_eq!(
            residue_generators[0].signature_matrix,
            DMatrix::from_row_slice(3, 3, &[-1., 0., 1., 0., 1., 0., 1., 0., 0.])
        );

        assert_eq!(
            residue_generators[1].signature_matrix,
            DMatrix::from_row_slice(3, 3, &[1., 0., 0., 0., 1., 0., -1., 0., 1.])
        );

        let residues_per_tree =
            spanning_tree_generator.get_residues_per_tree(&test_sigmas_1, &contour_closure, true);
        println!("residues: {:?}", residues_per_tree);
        assert!(residues_per_tree.is_empty());

        let residues_per_tree =
            spanning_tree_generator.get_residues_per_tree(&test_sigmas_2, &contour_closure, true);

        assert_eq!(residues_per_tree.len(), 1);
    }

    #[test]
    fn test_ltd() {
        let loop_line_signatures = vec![
            vec![1, 0, 0].into(),
            vec![0, 1, 0].into(),
            vec![0, 0, 1].into(),
            vec![1, -1, 0].into(),
            vec![-1, 0, 1].into(),
            vec![0, 1, -1].into(),
        ];

        let cut_structure_generator = CutStructureGenerator::new(loop_line_signatures);
        let countour_closure = [
            ContourClosure::Below,
            ContourClosure::Below,
            ContourClosure::Above,
        ];

        let cut_structure = cut_structure_generator.generate_structure(&countour_closure, true);
        assert_eq!(cut_structure.len(), 16);
        assert_eq!(cut_structure[0], vec![1., 1., -1., 0., 0., 0.]);
        assert_eq!(cut_structure[1], vec![-1., 1., 0., 0., -1., 0.]);
        assert_eq!(cut_structure[2], vec![1., -1., 0., 0., 0., 1.]);
        assert_eq!(cut_structure[3], vec![1., 0., -1., -1., 0., 0.]);
        assert_eq!(cut_structure[4], vec![1., 0., -1., 0., 0., 1.]);
        assert_eq!(cut_structure[5], vec![-1., 0., 0., -1., -1., 0.]);
        assert_eq!(cut_structure[6], vec![-1., 0., 0., 1., 0., 1.]);
        assert_eq!(cut_structure[7], vec![-1., 0., 0., 0., -1., 1.]);
        assert_eq!(cut_structure[8], vec![0., 1., -1., 1., 0., 0.]);
        assert_eq!(cut_structure[9], vec![0., 1., -1., 0., -1., 0.]);
        assert_eq!(cut_structure[10], vec![0., -1., 0., -1., -1., 0.]);
        assert_eq!(cut_structure[11], vec![0., -1., 0., 1., 0., 1.]);
        assert_eq!(cut_structure[12], vec![0., -1., 0., 0., -1., 1.]);
        assert_eq!(cut_structure[13], vec![0., 0., -1., -1., -1., 0.]);
        assert_eq!(cut_structure[14], vec![0., 0., -1., 1., 0., 1.]);
        assert_eq!(cut_structure[15], vec![0., 0., -1., 0., -1., 1.]);
    }
}
