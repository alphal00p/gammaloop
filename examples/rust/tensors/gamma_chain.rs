// Gamma chain example

use std::{ops::Neg, time::Instant};

use _gammaloop::tensor::{
    parametric::MixedTensor,
    ufo::{
        euclidean_four_vector, euclidean_four_vector_sym, gammasym, mink_four_vector,
        mink_four_vector_sym, param_euclidean_four_vector, param_mink_four_vector,
    },
    AbstractIndex, Contract, DenseTensor, FallibleMul, HasTensorData, HistoryStructure, IntoId,
    NumTensor, SparseTensor, TensorNetwork,
};

use num::traits::ToPrimitive;
use symbolica::domains::float::Complex;
use symbolica::representations::{Atom, Symbol};

// #[allow(dead_code)]
// fn gamma_trace<T>(minkindices: &[i32]) -> SparseTensor<Complex<T>>
// where
//     T: TrySmallestUpgrade<T, LCM = T>
//         + One
//         + Zero
//         + Default
//         + Copy
//         + Neg<Output = T>
//         + FallibleAddAssign<T>
//         + FallibleSubAssign<T>,
//     for<'a, 'b> &'a T: FallibleMul<&'b T, Output = T> + TrySmallestUpgrade<&'b T, LCM = T>,
// {
//     let mink = minkindices[0];
//     let mut result = gamma(usize::try_from(mink).unwrap().into(), (0.into(), 1.into()));
//     let mut contracting_index = 1;
//     for m in minkindices[1..].iter() {
//         let ui = contracting_index;
//         contracting_index += 1;

//         let uj = if contracting_index < minkindices.len() {
//             contracting_index
//         } else {
//             0
//         };

//         if *m > 0 {
//             let gamma: SparseTensor<Complex<T>> =
//                 gamma(usize::try_from(*m).unwrap().into(), (ui.into(), uj.into()));
//             result = gamma.contract(&result).unwrap();
//         } else {
//             result = gamma::<T>(
//                 usize::try_from(m.neg()).unwrap().into(),
//                 (ui.into(), uj.into()),
//             )
//             .contract(&result)
//             .unwrap();
//         }
//     }
//     result
// }

// #[allow(dead_code)]
// fn gamma_chain<T>(minkindices: &[i32]) -> SparseTensor<Complex<T>>
// where
//     T: Num
//         + Copy
//         + Neg<Output = T>
//         + Debug
//         + AddAssign
//         + SubAssign
//         + MulAssign
//         + DivAssign
//         + RemAssign
//         + Default,
// {
//     let mink = minkindices[0];
//     let mut result = gamma(usize::try_from(mink).unwrap().into(), (0.into(), 1.into()));
//     let mut contracting_index = 1;
//     for m in minkindices[1..].iter() {
//         let ui = contracting_index;
//         contracting_index += 1;

//         let uj = contracting_index;

//         if *m > 0 {
//             let gamma: SparseTensor<Complex<T>> =
//                 gamma(usize::try_from(*m).unwrap().into(), (ui.into(), uj.into()));
//             result = gamma.contract(&result).unwrap();
//         } else {
//             result = gamma::<T>(
//                 AbstractIndex::from(usize::try_from(m.neg()).unwrap() + 10000),
//                 (ui.into(), uj.into()),
//             )
//             .contract(&result)
//             .unwrap();
//         }
//     }
//     result
// }

#[allow(dead_code)]
fn gamma_net(
    minkindices: &[i32],
    vbar: [Complex<f64>; 4],
    u: [Complex<f64>; 4],
) -> TensorNetwork<NumTensor<HistoryStructure<Symbol>>> {
    let mut i = 0;
    let mut contracting_index: AbstractIndex = 0.into();
    let mut result: Vec<NumTensor<HistoryStructure<Symbol>>> =
        vec![euclidean_four_vector_sym(contracting_index, &vbar).into()];
    for m in minkindices {
        let ui = contracting_index;
        contracting_index += 1.into();
        let uj = contracting_index;
        if *m > 0 {
            let p = [
                Complex::<f64>::new(1.0 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex::<f64>::new(1.1 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex::<f64>::new(1.2 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex::<f64>::new(1.3 + 0.01 * i.to_f64().unwrap(), 0.0),
            ];
            i += 1;
            result.push(mink_four_vector_sym(usize::try_from(*m).unwrap().into(), &p).into());
            result.push(gammasym(usize::try_from(*m).unwrap().into(), (ui, uj)).into());
        } else {
            result.push(
                gammasym(
                    AbstractIndex::from(usize::try_from(m.neg()).unwrap() + 10000),
                    (ui, uj),
                )
                .into(),
            );
        }
    }
    result.push(euclidean_four_vector_sym(contracting_index, &u).into());
    TensorNetwork::from(result)
}

#[allow(dead_code)]
fn defered_chain(
    minkindices: &[i32],
    gamma_chain: &SparseTensor<Complex<f64>>,
    vbar: [Complex<f64>; 4],
    u: [Complex<f64>; 4],
) -> DenseTensor<Complex<f64>> {
    let mut result = euclidean_four_vector(0.into(), &vbar);
    result = gamma_chain.contract(&result).unwrap();
    result = result
        .contract(&euclidean_four_vector(minkindices.len().into(), &u))
        .unwrap();

    let mut i = 0;
    for m in minkindices {
        if *m > 0 {
            let p = [
                Complex::<f64>::new(1.0 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex::<f64>::new(1.1 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex::<f64>::new(1.2 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex::<f64>::new(1.3 + 0.01 * i.to_f64().unwrap(), 0.0),
            ];
            i += 1;
            let pmu = mink_four_vector(usize::try_from(*m).unwrap().into(), &p);
            result = pmu.contract(&result).unwrap();
        }
    }
    result
}

#[allow(dead_code)]
fn gamma_net_param(minkindices: &[i32]) -> TensorNetwork<MixedTensor<HistoryStructure<Symbol>>> {
    let mut i = 0;
    let mut contracting_index: AbstractIndex = 0.into();
    let mut result: Vec<MixedTensor<HistoryStructure<Symbol>>> =
        vec![param_euclidean_four_vector(contracting_index, "vbar".into_id()).into()];
    for m in minkindices {
        let ui = contracting_index;
        contracting_index += 1.into();
        let uj = contracting_index;
        if *m > 0 {
            let pname = format!("p{}", i).into_id();
            i += 1;
            result.push(param_mink_four_vector(usize::try_from(*m).unwrap().into(), pname).into());
            result.push(gammasym(usize::try_from(*m).unwrap().into(), (ui, uj)).into());
        } else {
            result.push(
                gammasym(
                    AbstractIndex::from(usize::try_from(m.neg()).unwrap() + 10000),
                    (ui, uj),
                )
                .into(),
            );
        }
    }
    result.push(param_euclidean_four_vector(contracting_index, "u".into_id()).into());
    TensorNetwork::from(result)
}

#[allow(dead_code)]
fn dump_c_with_func(_levels: Vec<Vec<(Symbol, Vec<Atom>)>>) {}

#[allow(unused_variables)]
fn main() {
    let one = Complex::<f64>::new(1.0, 0.0);
    let zero = Complex::<f64>::new(0.0, 0.0);

    let vbar = [
        one.mul_fallible(3.0).unwrap(),
        one.mul_fallible(3.1).unwrap(),
        one.mul_fallible(3.2).unwrap(),
        one.mul_fallible(3.3).unwrap(),
    ];
    let u = [
        one.mul_fallible(4.0).unwrap(),
        one.mul_fallible(4.1).unwrap(),
        one.mul_fallible(4.2).unwrap(),
        one.mul_fallible(4.3).unwrap(),
    ];

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

    let mut chain = gamma_net(&vec, vbar, u);
    println!("{}", chain.graph.edges.len());
    println!("{}", chain.graph.nodes.len());
    println!("{}", chain.graph.involution.len());
    println!("{}", chain.graph.neighbors.len());

    println!("{}", chain.dot());
    let start = Instant::now();
    chain.contract();
    let duration = start.elapsed();
    let durationfull = startfull.elapsed();

    println!("{}", chain.dot());

    println!(
        "Gamma net with {} gammas, fully numeric, takes {:?} for the contraction, and {:?} with initialization",
        vec.len(),
        duration,
        durationfull
    );

    // [Complex { re: 5.341852612369398e16, im: -136854212797686.44 }]
    // [Complex { re: 5.3418526123694e16, im: -136854212797684.0 }] for 20, 24
    println!(
        "Result: {:?}",
        chain.result().try_as_complex().unwrap().data()
    );

    // // let mut chain = gamma_net(&vec, vbar, u);
    // let ws: Workspace<Linear> = Workspace::new();

    // let atom = Atom::parse("A+P", &mut  &ws).unwrap();

    // let printops = PrintOptions {
    //     terms_on_new_line: false,
    //     color_top_level_sum: false,
    //     color_builtin_functions: false,
    //     print_finite_field: false,
    //     explicit_rational_polynomial: false,
    //     multiplication_operator: '*',
    //     square_brackets_for_function: false,
    //     number_thousands_separator: None,
    //     num_exp_as_superscript: false,
    //     latex: false,
    // };
    // let print = AtomPrinter::new_with_options(atom.as_view(), printops, &state);

    // let satom = format!("{}", print);
    // let natom = Atom::parse(&satom, &mut  &ws).unwrap();

    // println!("Print {}", natom.printer(&state));
    // // // println!("{}", chain.dot());
    // // let start = Instant::now();
    // // chain.contract_sym(& &ws);
    // // let duration = start.elapsed();
    // // println!(
    // //     "Benchmark net {:?} gammas, size in {:?}",
    // //     vec.len(),
    // //     duration,
    // // );

    // // println!("{:?}", chain.result().is_scalar());

    // // println!("{}", chain.result().structure());

    // // let mut chain_param = gamma_net_param(&vec, &mut  &ws);

    // // println!("{}", chain_param.dot());
    // // let start = Instant::now();
    // // chain_param.contract_sym_depth(5, & &ws);
    // // let duration = start.elapsed();
    // // // println!(
    // // //     "Benchmark net param {:?} gammas, size in {:?}",
    // // //     vec.len(),
    // // //     duration,
    // // // );
    // // println!("{}", chain_param.dot());

    // // println!("{:?}", chain_param.result().is_scalar());

    // // println!("{}", chain_param.result().structure());

    // let mut chain_param = gamma_net_param(&vec, &mut  &ws);

    // println!("{}", chain_param.dotsym(&state));
    // let params: Vec<Atom> = chain_param
    //     .clone()
    //     .to_symbolic_tensor_vec()
    //     .into_iter()
    //     .flat_map(|x| x.data())
    //     .collect();

    // let paramstr = params
    //     .iter()
    //     .map(|a| {
    //         format!(
    //             "{}",
    //             AtomPrinter::new_with_options(a.as_view(), printops, &state)
    //         )
    //     })
    //     .collect::<Vec<_>>();

    // serde_yaml::to_writer(std::fs::File::create("params.yaml").unwrap(), &paramstr).unwrap();

    // let paramstr: Vec<String> =
    //     serde_yaml::from_reader(std::fs::File::open("params.yaml").unwrap()).unwrap();

    // let params: Vec<Atom> = paramstr
    //     .iter()
    //     .map(|x| Atom::parse(x, &mut  &ws).unwrap())
    //     .collect();

    // for p in params {
    //     print!("{}", p.printer(&state));
    // }

    // chain_param.contract_sym_depth(9, & &ws);

    // let mut shadow = chain_param.symbolic_shadow("S", &mut  &ws);
    // // println!("{}", chain_param.dotsym(&state));
    // let a = chain_param
    //     .clone()
    //     .to_symbolic_tensor_vec()
    //     .into_iter()
    //     .map(|x| {
    //         (
    //             state.get_name(*x.name().unwrap()).clone(),
    //             x.data()
    //                 .into_iter()
    //                 .map(|x| {
    //                     format!(
    //                         "{}",
    //                         AtomPrinter::new_with_options(x.as_view(), printops, &state)
    //                     )
    //                     .into()
    //                 })
    //                 .collect(),
    //         )
    //     })
    //     .collect();

    // let amap = chain_param
    //     .to_symbolic_tensor_vec()
    //     .into_iter()
    //     .map(|x| x.symhashmap(*x.name().unwrap(), &mut  &ws))
    //     .collect::<Vec<_>>();

    // let amapstr = amap
    //     .iter()
    //     .map(|x| {
    //         let mut a = HashMap::new();
    //         for (k, v) in x.iter() {
    //             a.insert(
    //                 format!(
    //                     "{}",
    //                     AtomPrinter::new_with_options(k.as_view(), printops, &state)
    //                 )
    //                 .into(),
    //                 format!(
    //                     "{}",
    //                     AtomPrinter::new_with_options(v.as_view(), printops, &state)
    //                 )
    //                 .into(),
    //             );
    //         }
    //         a
    //     })
    //     .collect::<Vec<_>>();

    // println!("{}", shadow.dotsym(&state));
    // let start = Instant::now();
    // shadow.contract_sym_depth(10, & &ws);
    // let duration = start.elapsed();

    // // println!(
    // //     "Shadow net param {:?} gammas, size in {:?}",
    // //     vec.len(),
    // //     duration,
    // // );

    // let mut shadow2 = shadow.symbolic_shadow("T", &mut  &ws);
    // println!("{}", shadow.dotsym(&state));
    // let b: Vec<(String, Vec<String>)> = shadow
    //     .clone()
    //     .to_symbolic_tensor_vec()
    //     .into_iter()
    //     .map(|x| {
    //         (
    //             state.get_name(*x.name().unwrap()).clone(),
    //             x.data()
    //                 .into_iter()
    //                 .map(|x| {
    //                     format!(
    //                         "{}",
    //                         AtomPrinter::new_with_options(x.as_view(), printops, &state)
    //                     )
    //                     .into()
    //                 })
    //                 .collect(),
    //         )
    //     })
    //     .collect();

    // let bmap = shadow
    //     .to_symbolic_tensor_vec()
    //     .into_iter()
    //     .map(|x| x.symhashmap(*x.name().unwrap(), &mut  &ws))
    //     .collect::<Vec<_>>();

    // let bmapstr = bmap
    //     .iter()
    //     .map(|x| {
    //         let mut a = HashMap::new();
    //         for (k, v) in x.iter() {
    //             a.insert(
    //                 format!(
    //                     "{}",
    //                     AtomPrinter::new_with_options(k.as_view(), printops, &state)
    //                 )
    //                 .into(),
    //                 format!(
    //                     "{}",
    //                     AtomPrinter::new_with_options(v.as_view(), printops, &state)
    //                 )
    //                 .into(),
    //             );
    //         }
    //         a
    //     })
    //     .collect::<Vec<_>>();

    // println!("{}", shadow2.dotsym(&state));
    // let start = Instant::now();
    // shadow2.contract_sym_depth(10, & &ws);
    // let duration = start.elapsed();

    // // println!(
    // //     "Shadow2 net param {:?} gammas, size in {:?}",
    // //     vec.len(),
    // //     duration,
    // // );

    // shadow2.namesym("U", &mut state);
    // println!("{}", shadow2.dotsym(&state));

    // let c: Vec<(String, Vec<Atom>)> = shadow2
    //     .clone()
    //     .to_symbolic_tensor_vec()
    //     .into_iter()
    //     .map(|x| (state.get_name(*x.name().unwrap()).clone(), x.data()))
    //     .collect();

    // let e: Vec<(String, Vec<String>)> = shadow2
    //     .clone()
    //     .to_symbolic_tensor_vec()
    //     .into_iter()
    //     .map(|x| {
    //         (
    //             state.get_name(*x.name().unwrap()).clone(),
    //             x.data()
    //                 .into_iter()
    //                 .map(|x| {
    //                     format!(
    //                         "{}",
    //                         AtomPrinter::new_with_options(x.as_view(), printops, &state)
    //                     )
    //                     .into()
    //                 })
    //                 .collect(),
    //         )
    //     })
    //     .collect();

    // let cmap = shadow2
    //     .to_symbolic_tensor_vec()
    //     .into_iter()
    //     .map(|x| x.symhashmap(*x.name().unwrap(), &mut  &ws))
    //     .collect::<Vec<_>>();

    // let cmapstr = cmap
    //     .iter()
    //     .map(|x| {
    //         let mut a = HashMap::new();
    //         for (k, v) in x.iter() {
    //             a.insert(
    //                 format!(
    //                     "{}",
    //                     AtomPrinter::new_with_options(k.as_view(), printops, &state)
    //                 )
    //                 .into(),
    //                 format!(
    //                     "{}",
    //                     AtomPrinter::new_with_options(v.as_view(), printops, &state)
    //                 )
    //                 .into(),
    //             );
    //         }
    //         a
    //     })
    //     .collect::<Vec<_>>();

    // let d = e
    //     .iter()
    //     .map(|(s, v)| {
    //         (
    //             s,
    //             v.iter()
    //                 .map(|x| {
    //                     println!("Hi: {}", x);
    //                     Atom::parse(x, &mut  &ws).unwrap()
    //                 })
    //                 .collect::<Vec<_>>(),
    //         )
    //     })
    //     .collect::<Vec<_>>();

    // println!("{:?}", e);

    // // for (s, v) in d.iter() {
    // //     for x in v.iter() {
    // //         println!("{}", x.printer(&state));
    // //     }
    // // }

    // for (s, v) in c.iter() {
    //     for x in v.iter() {
    //         println!("{}", x.printer(&state));
    //     }
    // }

    // let out: Vec<Vec<(String, Vec<String>)>> = vec![a, b, e];
    // let outmap: Vec<Vec<HashMap<String, String>>> = vec![amapstr, bmapstr, cmapstr];

    // serde_yaml::to_writer(std::fs::File::create("outmap.yaml").unwrap(), &outmap).unwrap();

    // serde_yaml::to_writer(std::fs::File::create("out.yaml").unwrap(), &out).unwrap();

    // let from_file: Vec<Vec<(String, Vec<String>)>> =
    //     serde_yaml::from_reader(std::fs::File::open("out.yaml").unwrap()).unwrap();

    // let from_file_map: Vec<Vec<HashMap<String, String>>> =
    //     serde_yaml::from_reader(std::fs::File::open("outmap.yaml").unwrap()).unwrap();

    // let levelsmap = from_file_map
    //     .iter()
    //     .map(|x| {
    //         x.iter()
    //             .map(|x| {
    //                 let mut a = HashMap::new();
    //                 for (k, v) in x.iter() {
    //                     a.insert(
    //                         Atom::parse(k, &mut  &ws).unwrap(),
    //                         Atom::parse(v, &mut  &ws).unwrap(),
    //                     );
    //                 }
    //                 a
    //             })
    //             .collect::<Vec<_>>()
    //     })
    //     .collect::<Vec<_>>();

    // let levels: Vec<Vec<(Symbol, Vec<Atom>)>> = from_file
    //     .iter()
    //     .map(|x| {
    //         x.iter()
    //             .map(|(s, v)| {
    //                 (
    //                     state.get_or_insert_fn(s, None).unwrap(),
    //                     v.iter()
    //                         .map(|x| Atom::parse(x, &mut  &ws).unwrap())
    //                         .collect::<Vec<_>>(),
    //                 )
    //             })
    //             .collect::<Vec<_>>()
    //     })
    //     .collect::<Vec<_>>();

    // dump_c_with_func(levels);

    // // println!("{:?}", shadow.result().is_scalar());

    // // println!("{:?}", shadow.result());

    // // println!("{:?}", chain.result());
    // // for (i, c) in chain.iter() {
    // //     if *c == Complex::<i8>::new(0, 0) {
    // //         print!("hi")
    // //     }
    // //     if *c == Complex::<i8>::new(-0, 0) {
    // //         print!("hello")
    // //     }
    // //     if *c == Complex::<i8>::new(-0, -0) {
    // //         print!("hello")
    // //     }
    // //     if *c == Complex::<i8>::new(0, -0) {
    // //         print!("hello")
    // //     }
    // //     print!("{}", c.re);
    // // }

    // let start = Instant::now();
    // // let chain = benchmark_chain(&vec, vbar, u);
    // let duration = start.elapsed();

    // // println!("{:?} in {:?}", chain, duration);

    // // let start = Instant::now();
    // // // let chain = symbolic_chain(&vec, &ws, &mut state);
    // // let duration = start.elapsed();

    // // // println!("{:?} in {:?}", chain, duration);
    // // let mut out = ws.new_atom();
    // // let s = chain.finish().data.remove(0);

    // // s.as_view().expand(&ws, & &mut out);

    // // println!("{}", out.printer(&state));

    // // let poly: MultivariatePolynomial<_, u8> = out
    // //     .as_view()
    // //     .to_polynomial(&RationalField::new(), None)
    // //     .unwrap();

    // // let (h, _ops, scheme) = poly.optimize_horner_scheme(4000);
    // // let mut i = h.to_instr(poly.nvars);

    // // println!(
    // //     "Number of operations={}, with scheme={:?}",
    // //     BorrowedHornerScheme::from(&h).op_count_cse(),
    // //     scheme,
    // // );

    // // i.fuse_operations();

    // // for _ in 0..100_000 {
    // //     if !i.common_pair_elimination() {
    // //         break;
    // //     }
    // //     i.fuse_operations();
    // // }

    // // let op_count = i.op_count();
    // // let o = i.to_output(poly.var_map.as_ref().unwrap().to_vec(), true);
    // // let o_f64 = o.convert::<f64>();

    // // println!("Writing output to evaluate.cpp");
    // // std::fs::write(
    // //     "evaluate.cpp",
    // //     format!(
    // //         "{}",
    // //         InstructionSetPrinter {
    // //             instr: &o,
    // //             state: &
    // //             mode: symbolica::poly::evaluate::InstructionSetMode::CPP(
    // //                 symbolica::poly::evaluate::InstructionSetModeCPPSettings {
    // //                     write_header_and_test: true,
    // //                     always_pass_output_array: false,
    // //                 }
    // //             )
    // //         }
    // //     ),
    // // )
    // // .unwrap();

    // // let mut evaluator = o_f64.evaluator();

    // // let start = Instant::now();
    // // assert!(!evaluator
    // //     .evaluate(&(0..poly.nvars).map(|x| x as f64 + 1.).collect::<Vec<_>>())
    // //     .is_empty());
    // // let duration = start.elapsed();

    // // println!("Final number of operations={}", op_count);
    // // println!("Evaluation = {:?}", duration);

    // // // evaluate with simd
    // // let o_f64x4 = o.convert::<f64x4>();
    // // let mut evaluator = o_f64x4.evaluator();

    // // println!(
    // //     "Evaluation with simd = {:?}",
    // //     evaluator.evaluate(
    // //         &(0..poly.nvars)
    // //             .map(|x| f64x4::new([x as f64 + 1., x as f64 + 2., x as f64 + 3., x as f64 + 4.]))
    // //             .collect::<Vec<_>>()
    // //     )[0]
    // // );
}
