// Gamma chain example

use std::{ops::Neg, time::Instant};
use wide::f64x4;

use _gammaloop::tensor::{
    ufo_spin_tensors::{gamma, sigma},
    AbstractIndex, ContractableWithDense, DenseTensor, Expr,
    Representation::Lorentz,
    Representation::{self, Euclidean},
    TensorStructure, VecSlotExtension,
};

use num::complex::Complex64;
use num::traits::{Num, ToPrimitive};
use symbolica::{
    poly::{
        evaluate::{BorrowedHornerScheme, InstructionSetPrinter},
        polynomial::MultivariatePolynomial,
    },
    representations::{Atom, Identifier},
    rings::rational::RationalField,
    state::{State, Workspace},
};

fn pslash(indices: (usize, usize), p: &[Complex64; 4]) -> DenseTensor<Complex64> {
    let minkindex = indices.0 + indices.1;

    let p = mink_four_vector(minkindex, p);
    gamma(minkindex, indices).contract_with_dense(&p).unwrap()
}

#[allow(dead_code)]
fn mink_four_vector<T: std::clone::Clone>(index: usize, p: &[T; 4]) -> DenseTensor<T> {
    DenseTensor::from_data(p, TensorStructure::from_idxsing(&[index], &[Lorentz(4)])).unwrap()
}

fn labeled_mink_four_vector(
    label: &str,
    index: AbstractIndex,
    ws: &Workspace,
    state: &mut State,
) -> DenseTensor<Atom> {
    let structure = TensorStructure::from_idxsing(&[index], &[Lorentz(4)]);

    DenseTensor::symbolic_labels(label, structure, ws, state)
    // .try_into()
    // .unwrap_or_else(|v: Vec<_>| panic!("Expected a Vec of length 4 but it was {}", v.len()))
}
#[allow(dead_code)]
fn numbered_labeled_mink_four_vector<'a>(
    label: Identifier,
    number: usize,
    index: AbstractIndex,
    state: &'a State,
    ws: &'a Workspace,
) -> DenseTensor<Expr<'a>> {
    let structure = TensorStructure::from_idxsing(&[index], &[Lorentz(4)]);
    DenseTensor::numbered_labeled_builder(number, label, structure, ws, state)
}

fn eucl_four_vector<T>(index: usize, p: [T; 4]) -> DenseTensor<T>
where
    T: Num + std::default::Default + std::clone::Clone,
{
    DenseTensor::from_data(&p, TensorStructure::from_idxsing(&[index], &[Euclidean(4)])).unwrap()
}

fn labeled_eucl_four_vector(
    label: &str,
    index: AbstractIndex,
    ws: &Workspace,
    state: &mut State,
) -> DenseTensor<Atom> {
    let structure = TensorStructure::from_idxsing(&[index], &[Euclidean(4)]);
    DenseTensor::symbolic_labels(label, structure, ws, state)
}

#[allow(dead_code)]

fn numbered_labeled_eucl_four_vector<'a>(
    label: Identifier,
    number: usize,
    index: AbstractIndex,
    state: &'a State,
    ws: &'a Workspace,
) -> DenseTensor<Expr<'a>> {
    let structure = TensorStructure::from_idxsing(&[index], &[Euclidean(4)]);
    DenseTensor::numbered_labeled_builder(number, label, structure, ws, state)
}

fn benchmark_chain(
    minkindices: &[i32],
    vbar: [Complex64; 4],
    u: [Complex64; 4],
) -> DenseTensor<Complex64> {
    let mut i = 0;
    let mut contracting_index = 0;
    let mut result = eucl_four_vector(contracting_index, vbar);
    for m in minkindices {
        if *m > 0 {
            let p = [
                Complex64::new(1.0 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.1 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.2 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.3 + 0.01 * i.to_f64().unwrap(), 0.0),
            ];
            i += 1;
            let pslash = pslash((contracting_index, contracting_index + 1), &p);
            result = pslash.contract_with_dense(&result).unwrap();
        } else {
            result = gamma(
                usize::try_from(m.neg()).unwrap(),
                (contracting_index, contracting_index + 1),
            )
            .contract_with_dense(&result)
            .unwrap();
        }
        contracting_index += 1;
    }
    result
        .contract_with_dense(&eucl_four_vector(contracting_index, u))
        .unwrap()
}
#[allow(dead_code)]
fn symbolic_chain_function<'a>(
    minkindices: &[i32],
    ws: &'a Workspace,
    state: &'a mut State,
) -> DenseTensor<Expr<'a>> {
    let mut contracting_index = 0;
    let id = state.get_or_insert_fn("p", None).unwrap();
    let vbar = labeled_eucl_four_vector("vbar", contracting_index, ws, state).builder(state, ws);
    let mut result = vbar;

    // for (i, m) in minkindices.iter().filter(|&x| *x > 0).enumerate() {
    //     let p = labeled_mink_four_vector(&format!("p{}", i), ws, state);
    //     internal_mom.push(p);
    // }

    let mut i = 0;

    for m in minkindices {
        if *m > 0 {
            let minkindex = 2 * contracting_index;
            let p = numbered_labeled_mink_four_vector(id, i, minkindex, state, ws);
            let pslash = gamma(minkindex, (contracting_index, contracting_index + 1))
                .to_symbolic(ws, state)
                .builder(state, ws)
                .contract_with_dense(&p)
                .unwrap();
            result = pslash.contract_with_dense(&result).unwrap();

            i += 1;
        } else {
            result = gamma(
                usize::try_from(m.neg()).unwrap(),
                (contracting_index, contracting_index + 1),
            )
            .to_symbolic_builder(ws, state)
            .contract_with_dense(&result)
            .unwrap();
        }
        contracting_index += 1;
    }
    let result = result.finish();

    let u = labeled_eucl_four_vector("u", contracting_index, ws, state);
    result
        .builder(state, ws)
        .contract_with_dense(&u.builder(state, ws))
        .unwrap()
}

fn symbolic_chain<'a>(
    minkindices: &[i32],
    ws: &'a Workspace,
    state: &'a mut State,
) -> DenseTensor<Expr<'a>> {
    let mut contracting_index = 0;
    let vbar = labeled_eucl_four_vector("vbar", contracting_index, ws, state);
    let mut result = vbar;

    // for (i, m) in minkindices.iter().filter(|&x| *x > 0).enumerate() {
    //     let p = labeled_mink_four_vector(&format!("p{}", i), ws, state);
    //     internal_mom.push(p);
    // }

    let mut i = 0;

    for m in minkindices {
        if *m > 0 {
            let minkindex = 2 * contracting_index;
            let p = labeled_mink_four_vector(&format!("p{}", i), minkindex, ws, state);
            let pslash = gamma(minkindex, (contracting_index, contracting_index + 1))
                .to_symbolic(ws, state)
                .builder(state, ws)
                .contract_with_dense(&p.builder(state, ws))
                .unwrap();
            result = pslash
                .contract_with_dense(&result.builder(state, ws))
                .unwrap()
                .finish();

            i += 1;
        } else {
            result = gamma(
                usize::try_from(m.neg()).unwrap(),
                (contracting_index, contracting_index + 1),
            )
            .to_symbolic_builder(ws, state)
            .contract_with_dense(&result.builder(state, ws))
            .unwrap()
            .finish();
        }
        contracting_index += 1;
    }
    let u = labeled_eucl_four_vector("u", contracting_index, ws, state);

    result
        .builder(state, ws)
        .contract_with_dense(&u.builder(state, ws))
        .unwrap()
}

#[allow(dead_code)]
fn identity(indices: (usize, usize), signature: Representation) -> DenseTensor<Complex64> {
    let structure = TensorStructure::from_idxsing(&[indices.0, indices.1], &[signature, signature]);
    let mut identity = DenseTensor::default(structure);
    for i in 0..signature.into() {
        identity.set(&[i, i], Complex64::new(1.0, 0.0));
    }
    identity
}

#[allow(unused_variables)]
fn main() {
    let start = Instant::now();

    let s = sigma((1, 2), (3, 3));
    let duration = start.elapsed();
    println!("{:?} in {:?}", s, duration);
    // let p = [Complex64::new(1.0, 0.0); 4];
    // let p1 = pslash((1, 2), p);

    // let trg = gamma(1, (2, 2)).internal_contract();

    let one = Complex64::new(1.0, 0.0);
    let zero = Complex64::new(0.0, 0.0);

    // let gammamu = gamma(1, (2, 3))
    //     .contract_with_dense(&identity((3, 4), Euclidean(4)))
    //     .unwrap();

    // // println!("{:?}", gammamu);

    // let gammunu = gamma(1, (4, 5))
    //     .contract_with_dense(&identity((5, 6), Euclidean(4)))
    //     .unwrap();

    // println!("{:?}", gammamu.contract_with_dense(&gammunu).unwrap());

    let vbar = [one * 3.0, one * 3.1, one * 3.2, one * 3.3];
    let u = [one * 4.0, one * 4.1, one * 4.2, one * 4.3];
    // let a = eucl_four_vector(1, vbar);
    // let b = eucl_four_vector(1, u);
    // println!("{:?}", a.contract_with_dense(&b));

    // let p11 = pslash((1, 2), [one, zero, zero, zero])
    //     .contract_with_dense(&pslash((2, 1), [one, zero, zero, zero]));

    // println!("P {:?}", p11);

    let spacings: [i32; 2] = [4, 5];
    let mut start = 1;
    let mut ranges = Vec::new();

    for &spacing in spacings.iter() {
        ranges.push((start..start + spacing).chain(std::iter::once(-1)));
        start += spacing;
    }

    let vec: Vec<i32> = ranges.into_iter().flatten().collect();

    let start = Instant::now();
    let chain = benchmark_chain(&vec, vbar, u);
    let duration = start.elapsed();

    println!("{:?} in {:?}", chain, duration);

    let mut state = State::new();
    let ws = Workspace::new();

    let start = Instant::now();
    let chain = symbolic_chain(&vec, &ws, &mut state);
    let duration = start.elapsed();

    // println!("{:?} in {:?}", chain, duration);
    let mut out = ws.new_atom();
    let s = chain.finish().data.remove(0);

    s.as_view().expand(&ws, &state, &mut out);

    println!("{}", out.printer(&state));

    let poly: MultivariatePolynomial<_, u8> = out
        .as_view()
        .to_polynomial(&RationalField::new(), None)
        .unwrap();

    let (h, _ops, scheme) = poly.optimize_horner_scheme(4000);
    let mut i = h.to_instr(poly.nvars);

    println!(
        "Number of operations={}, with scheme={:?}",
        BorrowedHornerScheme::from(&h).op_count_cse(),
        scheme,
    );

    i.fuse_operations();

    for _ in 0..100_000 {
        if !i.common_pair_elimination() {
            break;
        }
        i.fuse_operations();
    }

    let op_count = i.op_count();
    let o = i.to_output(poly.var_map.as_ref().unwrap().to_vec(), true);
    let o_f64 = o.convert::<f64>();

    println!("Writing output to evaluate.cpp");
    std::fs::write(
        "evaluate.cpp",
        format!(
            "{}",
            InstructionSetPrinter {
                instr: &o,
                state: &state,
                mode: symbolica::poly::evaluate::InstructionSetMode::CPP(
                    symbolica::poly::evaluate::InstructionSetModeCPPSettings {
                        write_header_and_test: true,
                        always_pass_output_array: false,
                    }
                )
            }
        ),
    )
    .unwrap();

    let mut evaluator = o_f64.evaluator();

    let start = Instant::now();
    assert!(!evaluator
        .evaluate(&(0..poly.nvars).map(|x| x as f64 + 1.).collect::<Vec<_>>())
        .is_empty());
    let duration = start.elapsed();

    println!("Final number of operations={}", op_count);
    println!("Evaluation = {:?}", duration);

    // evaluate with simd
    let o_f64x4 = o.convert::<f64x4>();
    let mut evaluator = o_f64x4.evaluator();

    println!(
        "Evaluation with simd = {:?}",
        evaluator.evaluate(
            &(0..poly.nvars)
                .map(|x| f64x4::new([x as f64 + 1., x as f64 + 2., x as f64 + 3., x as f64 + 4.]))
                .collect::<Vec<_>>()
        )[0]
    );
}
