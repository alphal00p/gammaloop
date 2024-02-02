// Gamma chain example

use num::Complex;
use smartstring::alias::String;
use std::{
    collections::HashMap,
    fmt::{format, Debug},
    ops::{AddAssign, DivAssign, Mul, MulAssign, Neg, RemAssign, SubAssign},
    time::Instant,
};

use _gammaloop::tensor::{
    mixed_tensor::{MixedTensor, MixedTensors},
    ufo_spin_tensors::{
        euclidean_four_vector, euclidean_four_vector_sym, gamma, gammasym, mink_four_vector,
        mink_four_vector_sym, param_euclidean_four_vector, param_mink_four_vector,
    },
    AbstractIndex, Contract, DenseTensor, Expr, HasTensorData, HasTensorStructure, IntoId,
    NumTensors,
    Representation::Euclidean,
    Representation::Lorentz,
    SparseTensor, TensorNetwork, TensorSkeleton,
};

use num::complex::Complex64;
use num::traits::{Num, ToPrimitive};
use symbolica::{
    printer::{AtomPrinter, PrintOptions},
    representations::{default::Linear, Atom, Identifier},
    state::{State, Workspace},
};

#[allow(dead_code)]
fn pslash<T>(indices: (usize, usize), p: &[Complex<T>; 4]) -> DenseTensor<Complex<T>>
where
    T: Num
        + std::default::Default
        + Copy
        + Neg<Output = T>
        + Debug
        + AddAssign
        + SubAssign
        + MulAssign
        + DivAssign
        + RemAssign,
{
    let minkindex = indices.0 + indices.1;

    let p: DenseTensor<num::Complex<T>> = mink_four_vector(minkindex, p);
    p.contract(&gamma(minkindex, indices)).unwrap()
}

#[allow(dead_code)]
fn labeled_mink_four_vector(
    label: &str,
    index: AbstractIndex,
    ws: &Workspace,
    state: &mut State,
) -> DenseTensor<Atom> {
    let structure = TensorSkeleton::from_idxsing(&[(index, Lorentz(4))], "p".into());

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
    let structure = TensorSkeleton::from_idxsing(&[(index, Lorentz(4))], "p".into());
    DenseTensor::numbered_labeled_builder(number, label, structure, ws, state)
}

#[allow(dead_code)]
fn labeled_eucl_four_vector(
    label: &str,
    index: AbstractIndex,
    ws: &Workspace,
    state: &mut State,
) -> DenseTensor<Atom> {
    let structure = TensorSkeleton::from_idxsing(&[(index, Euclidean(4))], "x".into());
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
    let structure = TensorSkeleton::from_idxsing(&[(index, Euclidean(4))], "x".into());
    DenseTensor::numbered_labeled_builder(number, label, structure, ws, state)
}

#[allow(dead_code)]
fn benchmark_chain(
    minkindices: &[i32],
    vbar: [Complex64; 4],
    u: [Complex64; 4],
) -> DenseTensor<Complex64> {
    let mut i = 0;
    let mut contracting_index = 0;
    let mut result = euclidean_four_vector(contracting_index, &vbar);
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
            result = pslash.contract(&result).unwrap();
        } else {
            result = gamma::<f64>(
                usize::try_from(m.neg()).unwrap(),
                (contracting_index, contracting_index + 1),
            )
            .contract(&result)
            .unwrap();
        }
        contracting_index += 1;
    }
    result
        .contract(&euclidean_four_vector(contracting_index, &u))
        .unwrap()
}
// #[allow(dead_code)]
// fn symbolic_chain_function<'a>(
//     minkindices: &[i32],
//     ws: &'a Workspace,
//     state: &'a mut State,
// ) -> DenseTensor<Expr<'a>> {
//     let mut contracting_index = 0;
//     let id = state.get_or_insert_fn("p", None).unwrap();
//     let vbar = labeled_eucl_four_vector("vbar", contracting_index, ws, state).builder(state, ws);
//     let mut result = vbar;

#[allow(dead_code)]
fn gamma_trace<T>(minkindices: &[i32]) -> SparseTensor<Complex<T>>
where
    T: Num
        + std::default::Default
        + Copy
        + Neg<Output = T>
        + Debug
        + AddAssign
        + SubAssign
        + MulAssign
        + DivAssign
        + RemAssign,
{
    let mink = minkindices[0];
    let mut result = gamma(usize::try_from(mink).unwrap(), (0, 1));
    let mut contracting_index = 1;
    for m in minkindices[1..].iter() {
        let ui = contracting_index;
        contracting_index += 1;

        let uj = if contracting_index < minkindices.len() {
            contracting_index
        } else {
            0
        };

        if *m > 0 {
            let gamma: SparseTensor<Complex<T>> = gamma(usize::try_from(*m).unwrap(), (ui, uj));
            result = gamma.contract(&result).unwrap();
        } else {
            result = gamma::<T>(usize::try_from(m.neg()).unwrap(), (ui, uj))
                .contract(&result)
                .unwrap();
        }
    }
    result
}

#[allow(dead_code)]
fn gamma_chain<T>(minkindices: &[i32]) -> SparseTensor<Complex<T>>
where
    T: Num
        + Copy
        + Neg<Output = T>
        + Debug
        + AddAssign
        + SubAssign
        + MulAssign
        + DivAssign
        + RemAssign
        + Default,
{
    let mink = minkindices[0];
    let mut result = gamma(usize::try_from(mink).unwrap(), (0, 1));
    let mut contracting_index = 1;
    for m in minkindices[1..].iter() {
        let ui = contracting_index;
        contracting_index += 1;

        let uj = contracting_index;

        if *m > 0 {
            let gamma: SparseTensor<Complex<T>> = gamma(usize::try_from(*m).unwrap(), (ui, uj));
            result = gamma.contract(&result).unwrap();
        } else {
            result = gamma::<T>(usize::try_from(m.neg()).unwrap() + 10000, (ui, uj))
                .contract(&result)
                .unwrap();
        }
    }
    result
}

#[allow(dead_code)]
fn gamma_net(
    minkindices: &[i32],
    vbar: [Complex64; 4],
    u: [Complex64; 4],
    state: &mut State,
) -> TensorNetwork<NumTensors<Identifier>> {
    let mut i = 0;
    let mut contracting_index = 0;
    let mut result: Vec<NumTensors<Identifier>> =
        vec![euclidean_four_vector_sym(contracting_index, &vbar, state).into()];
    for m in minkindices {
        let ui = contracting_index;
        contracting_index += 1;
        let uj = contracting_index;
        if *m > 0 {
            let p = [
                Complex64::new(1.0 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.1 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.2 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.3 + 0.01 * i.to_f64().unwrap(), 0.0),
            ];
            i += 1;
            result.push(mink_four_vector_sym(usize::try_from(*m).unwrap(), &p, state).into());
            result.push(gammasym(usize::try_from(*m).unwrap(), (ui, uj), state).into());
        } else {
            result
                .push(gammasym(usize::try_from(m.neg()).unwrap() + 10000, (ui, uj), state).into());
        }
    }
    result.push(euclidean_four_vector_sym(contracting_index, &u, state).into());
    TensorNetwork::new(result)
}

#[allow(dead_code)]
fn defered_chain(
    minkindices: &[i32],
    gamma_chain: &SparseTensor<Complex<f64>>,
    vbar: [Complex64; 4],
    u: [Complex64; 4],
) -> DenseTensor<Complex<f64>> {
    let mut result = euclidean_four_vector(0, &vbar);
    result = gamma_chain.contract(&result).unwrap();
    result = result
        .contract(&euclidean_four_vector(minkindices.len(), &u))
        .unwrap();

    let mut i = 0;
    for m in minkindices {
        if *m > 0 {
            let p = [
                Complex64::new(1.0 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.1 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.2 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.3 + 0.01 * i.to_f64().unwrap(), 0.0),
            ];
            i += 1;
            let pmu = mink_four_vector(usize::try_from(*m).unwrap(), &p);
            result = pmu.contract(&result).unwrap();
        }
    }
    result
}

#[allow(dead_code)]
fn gamma_net_param(
    minkindices: &[i32],
    state: &mut State,
    ws: &Workspace,
) -> TensorNetwork<MixedTensor<Identifier>> {
    let mut i = 0;
    let mut contracting_index = 0;
    let mut result: Vec<MixedTensor<Identifier>> =
        vec![
            param_euclidean_four_vector(contracting_index, "vbar".into_id(state), state, ws).into(),
        ];
    for m in minkindices {
        let ui = contracting_index;
        contracting_index += 1;
        let uj = contracting_index;
        if *m > 0 {
            let pname = format!("p{}", i).into_id(state);
            i += 1;
            result.push(
                param_mink_four_vector(usize::try_from(*m).unwrap(), pname, state, ws).into(),
            );
            result.push(gammasym(usize::try_from(*m).unwrap(), (ui, uj), state).into());
        } else {
            result
                .push(gammasym(usize::try_from(m.neg()).unwrap() + 10000, (ui, uj), state).into());
        }
    }
    result
        .push(param_euclidean_four_vector(contracting_index, "u".into_id(state), state, ws).into());
    TensorNetwork::new(result)
}

fn dump_c_with_func(levels: Vec<Vec<(Identifier, Vec<Atom>)>>) {}
fn dump_c(levels: Vec<Vec<HashMap<Identifier, Vec<Atom>>>>) {}

// #[allow(dead_code)]
// fn symbolic_chain_function<'a>(
//     minkindices: &[i32],
//     ws: &'a Workspace,
//     state: &'a mut State,
// ) -> DenseTensor<Expr<'a>> {
//     let mut contracting_index = 0;
//     let id = state.get_or_insert_fn("p", None).unwrap();
//     let vbar = labeled_eucl_four_vector("vbar", contracting_index, ws, state).builder(state, ws);
//     let mut result = vbar;

//     // for (i, m) in minkindices.iter().filter(|&x| *x > 0).enumerate() {
//     //     let p = labeled_mink_four_vector(&format!("p{}", i), ws, state);
//     //     internal_mom.push(p);
//     // }

//     let mut i = 0;

//     for m in minkindices {
//         if *m > 0 {
//             let minkindex = 2 * contracting_index;
//             let p = numbered_labeled_mink_four_vector(id, i, minkindex, state, ws);
//             let pslash = gamma(minkindex, (contracting_index, contracting_index + 1))
//                 .to_symbolic(ws, state)
//                 .builder(state, ws)
//                 .contract(&p)
//                 .unwrap();
//             result = pslash.contract(&result).unwrap();

//             i += 1;
//         } else {
//             result = gamma(
//                 usize::try_from(m.neg()).unwrap(),
//                 (contracting_index, contracting_index + 1),
//             )
//             .to_symbolic_builder(ws, state)
//             .contract(&result)
//             .unwrap();
//         }
//         contracting_index += 1;
//     }
//     let result = result.finish();

//     let u = labeled_eucl_four_vector("u", contracting_index, ws, state);
//     result
//         .builder(state, ws)
//         .contract(&u.builder(state, ws))
//         .unwrap()
// }

// fn symbolic_chain<'a>(
//     minkindices: &[i32],
//     ws: &'a Workspace,
//     state: &'a mut State,
// ) -> DenseTensor<Expr<'a>> {
//     let mut contracting_index = 0;
//     let vbar = labeled_eucl_four_vector("vbar", contracting_index, ws, state);
//     let mut result = vbar;

//     // for (i, m) in minkindices.iter().filter(|&x| *x > 0).enumerate() {
//     //     let p = labeled_mink_four_vector(&format!("p{}", i), ws, state);
//     //     internal_mom.push(p);
//     // }

//     let mut i = 0;

//     for m in minkindices {
//         if *m > 0 {
//             let minkindex = 2 * contracting_index;
//             let p = labeled_mink_four_vector(&format!("p{}", i), minkindex, ws, state);
//             let pslash = gamma(minkindex, (contracting_index, contracting_index + 1))
//                 .to_symbolic(ws, state)
//                 .builder(state, ws)
//                 .contract(&p.builder(state, ws))
//                 .unwrap();
//             result = pslash
//                 .contract(&result.builder(state, ws))
//                 .unwrap()
//                 .finish();

//             i += 1;
//         } else {
//             result = gamma(
//                 usize::try_from(m.neg()).unwrap(),
//                 (contracting_index, contracting_index + 1),
//             )
//             .to_symbolic_builder(ws, state)
//             .contract(&result.builder(state, ws))
//             .unwrap()
//             .finish();
//         }
//         contracting_index += 1;
//     }
//     let u = labeled_eucl_four_vector("u", contracting_index, ws, state);

//     result
//         .builder(state, ws)
//         .contract(&u.builder(state, ws))
//         .unwrap()
// }

#[allow(unused_variables)]
fn main() {
    let one = Complex64::new(1.0, 0.0);
    let zero = Complex64::new(0.0, 0.0);

    let vbar = [one * 3.0, one * 3.1, one * 3.2, one * 3.3];
    let u = [one * 4.0, one * 4.1, one * 4.2, one * 4.3];

    let spacings: [i32; 2] = [20, 24];
    let mut start = 1;
    let mut ranges = Vec::new();

    for &spacing in spacings.iter() {
        ranges.push((start..start + spacing).chain(std::iter::once(-1)));
        start += spacing;
    }

    let vec: Vec<i32> = ranges.into_iter().flatten().collect();

    println!("{:?}", vec);

    // // let vec = (1..=3).collect::<Vec<_>>();
    // // let start = Instant::now();
    // // let chain = gamma_chain(&vec); //, vbar, u);
    // // let duration = start.elapsed();

    // // println!(
    // //     "Gamma chain {:?} gammas, size {:?}  in {:?}",
    // //     vec.len(),
    // //     chain.size(),
    // //     duration
    // // );

    // // let start = Instant::now();
    // // let chain = defered_chain(&vec, &chain, vbar, u);
    // // let duration = start.elapsed();

    // // println!(
    // //     "Defered pslash with {:?} gammas, size {:?}  in {:?}, gives {:?}",
    // //     vec.len(),
    // //     chain.size(),
    // //     duration,
    // //     chain.data,
    // // );

    let startfull = Instant::now();
    let mut state = State::new();

    let mut chain = gamma_net(&vec, vbar, u, &mut state);
    let start = Instant::now();
    chain.contract();
    let duration = start.elapsed();
    let durationfull = startfull.elapsed();

    println!(
        "Gamma net with {} gammas, fully numeric, takes {:?} for the contraction, and {:?} with initialization",
        vec.len(),
        duration,
        durationfull
    );

    // [Complex { re: 5.3418526123694e16, im: -136854212797684.0 }] for 20, 24
    println!(
        "Result: {:?}",
        chain.result().try_as_complex().unwrap().data()
    );

    // let mut chain = gamma_net(&vec, vbar, u);
    let ws: Workspace<Linear> = Workspace::new();

    let atom = Atom::parse("A+P", &mut state, &ws).unwrap();

    let printops = PrintOptions {
        terms_on_new_line: false,
        color_top_level_sum: false,
        color_builtin_functions: false,
        print_finite_field: false,
        explicit_rational_polynomial: false,
        multiplication_operator: '*',
        square_brackets_for_function: false,
        number_thousands_separator: None,
        num_exp_as_superscript: false,
        latex: false,
    };
    let print = AtomPrinter::new_with_options(atom.as_view(), printops, &state);

    let satom = format!("{}", print);
    let natom = Atom::parse(&satom, &mut state, &ws).unwrap();

    println!("Print {}", natom.printer(&state));
    // // println!("{}", chain.dot());
    // let start = Instant::now();
    // chain.contract_sym(&state, &ws);
    // let duration = start.elapsed();
    // println!(
    //     "Benchmark net {:?} gammas, size in {:?}",
    //     vec.len(),
    //     duration,
    // );

    // println!("{:?}", chain.result().is_scalar());

    // println!("{}", chain.result().structure());

    // let mut chain_param = gamma_net_param(&vec, &mut state, &ws);

    // println!("{}", chain_param.dot());
    // let start = Instant::now();
    // chain_param.contract_sym_depth(5, &state, &ws);
    // let duration = start.elapsed();
    // // println!(
    // //     "Benchmark net param {:?} gammas, size in {:?}",
    // //     vec.len(),
    // //     duration,
    // // );
    // println!("{}", chain_param.dot());

    // println!("{:?}", chain_param.result().is_scalar());

    // println!("{}", chain_param.result().structure());

    let mut chain_param = gamma_net_param(&vec, &mut state, &ws);

    println!("{}", chain_param.dotsym(&state));
    let params: Vec<Atom> = chain_param
        .clone()
        .to_symbolic_tensor_vec()
        .into_iter()
        .flat_map(|x| x.data())
        .collect();

    let paramstr = params
        .iter()
        .map(|a| {
            format!(
                "{}",
                AtomPrinter::new_with_options(a.as_view(), printops, &state)
            )
        })
        .collect::<Vec<_>>();

    serde_yaml::to_writer(std::fs::File::create("params.yaml").unwrap(), &paramstr).unwrap();

    let paramstr: Vec<String> =
        serde_yaml::from_reader(std::fs::File::open("params.yaml").unwrap()).unwrap();

    let params: Vec<Atom> = paramstr
        .iter()
        .map(|x| Atom::parse(x, &mut state, &ws).unwrap())
        .collect();

    for p in params {
        print!("{}", p.printer(&state));
    }

    chain_param.contract_sym_depth(9, &state, &ws);

    let mut shadow = chain_param.symbolic_shadow("S", &mut state, &ws);
    println!("{}", chain_param.dotsym(&state));
    let a = chain_param
        .clone()
        .to_symbolic_tensor_vec()
        .into_iter()
        .map(|x| {
            (
                state.get_name(*x.global_name().unwrap()).clone(),
                x.data()
                    .into_iter()
                    .map(|x| {
                        format!(
                            "{}",
                            AtomPrinter::new_with_options(x.as_view(), printops, &state)
                        )
                        .into()
                    })
                    .collect(),
            )
        })
        .collect();

    let amap = chain_param
        .to_symbolic_tensor_vec()
        .into_iter()
        .map(|x| x.symhashmap(*x.global_name().unwrap(), &mut state, &ws))
        .collect::<Vec<_>>();

    let amapstr = amap
        .iter()
        .map(|x| {
            let mut a = HashMap::new();
            for (k, v) in x.iter() {
                a.insert(
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(k.as_view(), printops, &state)
                    )
                    .into(),
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(v.as_view(), printops, &state)
                    )
                    .into(),
                );
            }
            a
        })
        .collect::<Vec<_>>();

    println!("{}", shadow.dotsym(&state));
    let start = Instant::now();
    shadow.contract_sym_depth(10, &state, &ws);
    let duration = start.elapsed();

    // println!(
    //     "Shadow net param {:?} gammas, size in {:?}",
    //     vec.len(),
    //     duration,
    // );

    let mut shadow2 = shadow.symbolic_shadow("T", &mut state, &ws);
    println!("{}", shadow.dotsym(&state));
    let b: Vec<(String, Vec<String>)> = shadow
        .clone()
        .to_symbolic_tensor_vec()
        .into_iter()
        .map(|x| {
            (
                state.get_name(*x.global_name().unwrap()).clone(),
                x.data()
                    .into_iter()
                    .map(|x| {
                        format!(
                            "{}",
                            AtomPrinter::new_with_options(x.as_view(), printops, &state)
                        )
                        .into()
                    })
                    .collect(),
            )
        })
        .collect();

    let bmap = shadow
        .to_symbolic_tensor_vec()
        .into_iter()
        .map(|x| x.symhashmap(*x.global_name().unwrap(), &mut state, &ws))
        .collect::<Vec<_>>();

    let bmapstr = bmap
        .iter()
        .map(|x| {
            let mut a = HashMap::new();
            for (k, v) in x.iter() {
                a.insert(
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(k.as_view(), printops, &state)
                    )
                    .into(),
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(v.as_view(), printops, &state)
                    )
                    .into(),
                );
            }
            a
        })
        .collect::<Vec<_>>();

    println!("{}", shadow2.dotsym(&state));
    let start = Instant::now();
    shadow2.contract_sym_depth(10, &state, &ws);
    let duration = start.elapsed();

    // println!(
    //     "Shadow2 net param {:?} gammas, size in {:?}",
    //     vec.len(),
    //     duration,
    // );

    shadow2.namesym("U", &mut state);
    println!("{}", shadow2.dotsym(&state));

    let c: Vec<(String, Vec<Atom>)> = shadow2
        .clone()
        .to_symbolic_tensor_vec()
        .into_iter()
        .map(|x| (state.get_name(*x.global_name().unwrap()).clone(), x.data()))
        .collect();

    let e: Vec<(String, Vec<String>)> = shadow2
        .clone()
        .to_symbolic_tensor_vec()
        .into_iter()
        .map(|x| {
            (
                state.get_name(*x.global_name().unwrap()).clone(),
                x.data()
                    .into_iter()
                    .map(|x| {
                        format!(
                            "{}",
                            AtomPrinter::new_with_options(x.as_view(), printops, &state)
                        )
                        .into()
                    })
                    .collect(),
            )
        })
        .collect();

    let cmap = shadow2
        .to_symbolic_tensor_vec()
        .into_iter()
        .map(|x| x.symhashmap(*x.global_name().unwrap(), &mut state, &ws))
        .collect::<Vec<_>>();

    let cmapstr = cmap
        .iter()
        .map(|x| {
            let mut a = HashMap::new();
            for (k, v) in x.iter() {
                a.insert(
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(k.as_view(), printops, &state)
                    )
                    .into(),
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(v.as_view(), printops, &state)
                    )
                    .into(),
                );
            }
            a
        })
        .collect::<Vec<_>>();

    let d = e
        .iter()
        .map(|(s, v)| {
            (
                s,
                v.iter()
                    .map(|x| {
                        println!("Hi: {}", x);
                        Atom::parse(x, &mut state, &ws).unwrap()
                    })
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Vec<_>>();

    println!("{:?}", e);

    // for (s, v) in d.iter() {
    //     for x in v.iter() {
    //         println!("{}", x.printer(&state));
    //     }
    // }

    for (s, v) in c.iter() {
        for x in v.iter() {
            println!("{}", x.printer(&state));
        }
    }

    let out: Vec<Vec<(String, Vec<String>)>> = vec![a, b, e];
    let outmap: Vec<Vec<HashMap<String, String>>> = vec![amapstr, bmapstr, cmapstr];

    serde_yaml::to_writer(std::fs::File::create("outmap.yaml").unwrap(), &outmap).unwrap();

    serde_yaml::to_writer(std::fs::File::create("out.yaml").unwrap(), &out).unwrap();

    let from_file: Vec<Vec<(String, Vec<String>)>> =
        serde_yaml::from_reader(std::fs::File::open("out.yaml").unwrap()).unwrap();

    let from_file_map: Vec<Vec<HashMap<String, String>>> =
        serde_yaml::from_reader(std::fs::File::open("outmap.yaml").unwrap()).unwrap();

    let levelsmap = from_file_map
        .iter()
        .map(|x| {
            x.iter()
                .map(|x| {
                    let mut a = HashMap::new();
                    for (k, v) in x.iter() {
                        a.insert(
                            Atom::parse(k, &mut state, &ws).unwrap(),
                            Atom::parse(v, &mut state, &ws).unwrap(),
                        );
                    }
                    a
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let levels: Vec<Vec<(Identifier, Vec<Atom>)>> = from_file
        .iter()
        .map(|x| {
            x.iter()
                .map(|(s, v)| {
                    (
                        state.get_or_insert_fn(s, None).unwrap(),
                        v.iter()
                            .map(|x| Atom::parse(x, &mut state, &ws).unwrap())
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    dump_c_with_func(levels);

    // println!("{:?}", shadow.result().is_scalar());

    // println!("{:?}", shadow.result());

    // println!("{:?}", chain.result());
    // for (i, c) in chain.iter() {
    //     if *c == Complex::<i8>::new(0, 0) {
    //         print!("hi")
    //     }
    //     if *c == Complex::<i8>::new(-0, 0) {
    //         print!("hello")
    //     }
    //     if *c == Complex::<i8>::new(-0, -0) {
    //         print!("hello")
    //     }
    //     if *c == Complex::<i8>::new(0, -0) {
    //         print!("hello")
    //     }
    //     print!("{}", c.re);
    // }

    let start = Instant::now();
    // let chain = benchmark_chain(&vec, vbar, u);
    let duration = start.elapsed();

    // println!("{:?} in {:?}", chain, duration);

    // let start = Instant::now();
    // // let chain = symbolic_chain(&vec, &ws, &mut state);
    // let duration = start.elapsed();

    // // println!("{:?} in {:?}", chain, duration);
    // let mut out = ws.new_atom();
    // let s = chain.finish().data.remove(0);

    // s.as_view().expand(&ws, &state, &mut out);

    // println!("{}", out.printer(&state));

    // let poly: MultivariatePolynomial<_, u8> = out
    //     .as_view()
    //     .to_polynomial(&RationalField::new(), None)
    //     .unwrap();

    // let (h, _ops, scheme) = poly.optimize_horner_scheme(4000);
    // let mut i = h.to_instr(poly.nvars);

    // println!(
    //     "Number of operations={}, with scheme={:?}",
    //     BorrowedHornerScheme::from(&h).op_count_cse(),
    //     scheme,
    // );

    // i.fuse_operations();

    // for _ in 0..100_000 {
    //     if !i.common_pair_elimination() {
    //         break;
    //     }
    //     i.fuse_operations();
    // }

    // let op_count = i.op_count();
    // let o = i.to_output(poly.var_map.as_ref().unwrap().to_vec(), true);
    // let o_f64 = o.convert::<f64>();

    // println!("Writing output to evaluate.cpp");
    // std::fs::write(
    //     "evaluate.cpp",
    //     format!(
    //         "{}",
    //         InstructionSetPrinter {
    //             instr: &o,
    //             state: &state,
    //             mode: symbolica::poly::evaluate::InstructionSetMode::CPP(
    //                 symbolica::poly::evaluate::InstructionSetModeCPPSettings {
    //                     write_header_and_test: true,
    //                     always_pass_output_array: false,
    //                 }
    //             )
    //         }
    //     ),
    // )
    // .unwrap();

    // let mut evaluator = o_f64.evaluator();

    // let start = Instant::now();
    // assert!(!evaluator
    //     .evaluate(&(0..poly.nvars).map(|x| x as f64 + 1.).collect::<Vec<_>>())
    //     .is_empty());
    // let duration = start.elapsed();

    // println!("Final number of operations={}", op_count);
    // println!("Evaluation = {:?}", duration);

    // // evaluate with simd
    // let o_f64x4 = o.convert::<f64x4>();
    // let mut evaluator = o_f64x4.evaluator();

    // println!(
    //     "Evaluation with simd = {:?}",
    //     evaluator.evaluate(
    //         &(0..poly.nvars)
    //             .map(|x| f64x4::new([x as f64 + 1., x as f64 + 2., x as f64 + 3., x as f64 + 4.]))
    //             .collect::<Vec<_>>()
    //     )[0]
    // );
}
