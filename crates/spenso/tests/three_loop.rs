// use std::{fs::File, io::BufReader};

// use ahash::{AHashMap, HashMap, HashMapExt};

// use approx::{assert_relative_eq, RelativeEq};
// use common::WEYLIB;
// use spenso::{
//     complex::Complex,
//     network::Network,
//     parametric::{atomcore::TensorAtomOps, MixedTensor},
//     structure::{HasStructure, SmartShadowStructure},
//     symbolic::SymbolicTensor,
//     symbolica_utils::SerializableAtom,
// };
// use symbolica::{
//     atom::{Atom, AtomCore, AtomView, Symbol},
//     domains::rational::Rational,
//     evaluate::{CompileOptions, FunctionMap, InlineASM},
//     id::Replacement,
//     parse, parse_lit, symbol,
// };

// mod common;

// use symbolica::domains::float::Complex as SymComplex;

// fn main() {
//     // let atom = parse_lit!(
//     //     -64 / 729 * G
//     //         ^ 4 * ee
//     //         ^ 6
//     // *(MT*id(bis(4,47),bis(4,135))+Q(15,mink(4,149)))*gamma(mink(4,149),bis(4,47),bis(4,135)))
//     // *(MT*id(bis(4,83),bis(4,46))+Q(6,mink(4,138)))*gamma(mink(4,138),bis(4,83),bis(4,46))
//     // *(MT*id(bis(4,88),bis(4,82))+gamma(mink(4,140),bis(4,88),bis(4,82)))*Q(7,mink(4,140))
//     // *(MT*id(bis(4,96),bis(4,142))+gamma(mink(4,141),bis(4,96),bis(4,142)))*Q(8,mink(4,141))
//     // *(MT*id(bis(4,103),bis(4,95))+gamma(mink(4,143),bis(4,103),bis(4,95)))*Q(9,mink(4,143))
//     // *(MT*id(bis(4,110),bis(4,102))+gamma(mink(4,144),bis(4,110),bis(4,102)))*Q(10,mink(4,144))
//     // *(MT*id(bis(4,117),bis(4,109))+gamma(mink(4,145),bis(4,117),bis(4,109)))*Q(11,mink(4,145))
//     // *(MT*id(bis(4,122),bis(4,116))+gamma(mink(4,146),bis(4,122),bis(4,116)))*Q(12,mink(4,146))
//     // *(MT*id(bis(4,129),bis(4,123))+gamma(mink(4,147),bis(4,129),bis(4,123)))*Q(13,mink(4,147))
//     // *(MT*id(bis(4,134),bis(4,130))+gamma(mink(4,148),bis(4,134),bis(4,130)))*Q(14,mink(4,148))
//     // *gamma(mink(4,45),bis(4,47),bis(4,46))*gamma(mink(4,81),bis(4,83),bis(4,82))*gamma(mink(4,87),bis(4,88),bis(4,142))*gamma(mink(4,94),bis(4,96),bis(4,95))
//     // *gamma(mink(4,101),bis(4,103),bis(4,102))*gamma(mink(4,108),bis(4,110),bis(4,109))*gamma(mink(4,115),bis(4,117),bis(4,116))*gamma(mink(4,121),bis(4,122),bis(4,123))
//     // *gamma(mink(4,128),bis(4,129),bis(4,130))*gamma(mink(4,133),bis(4,134),bis(4,135))*g(mink(4,121),mink(4,87))*g(mink(4,133),mink(4,101))).unwrap();
//     // *ϵ(0,mink(4,45))*ϵ(1,mink(4,81))*ϵbar(2,mink(4,94))*ϵbar(3,mink(4,108))*ϵbar(4,lord(4,115))*ϵbar(5,mink(4,128))
//     // *id(coaf(3,46),cof(3,47))*id(coaf(3,82),cof(3,83))*id(coaf(3,95),cof(3,96))*id(coaf(3,109),cof(3,110))*id(coaf(3,116),cof(3,117))*id(coaf(3,130),cof(3,129))",
//     // *T(coad(8,87),cof(3,88),coaf(3,46))*T(coad(8,101),cof(3,103),coaf(3,102))*T(coad(8,121),cof(3,122),coaf(3,123))*T(coad(8,133),cof(3,134),coaf(3,135))",
//     // let sym_tensor: SymbolicTensor = Atom::Zero.try_into().unwrap();

//     // let mut network = sym_tensor
//     //     .to_network(WEYLIB.read().as_ref().unwrap())
//     //     .unwrap();

//     // let file = File::open("./examples/data.json").unwrap();
//     // let reader = BufReader::new(file);

//     // let data_string_map: AHashMap<String, Complex<f64>> = serde_json::from_reader(reader).unwrap();

//     // let file = File::open("./examples/const.json").unwrap();
//     // let reader = BufReader::new(file);

//     // let const_string_map: AHashMap<String, Complex<f64>> = serde_json::from_reader(reader).unwrap();

//     // let data_atom_map: (Vec<Atom>, Vec<Complex<f64>>) = data_string_map
//     //     .into_iter()
//     //     .map(|(k, v)| (parse!(&k).unwrap(), v))
//     //     .unzip();

//     // let mut const_atom_map: AHashMap<Symbol, Complex<f64>> = const_string_map
//     //     .into_iter()
//     //     .map(|(k, v)| (symbol!(&k), v))
//     //     .collect();

//     // const_atom_map.insert(Atom::I, Complex::i());

//     // let mut const_map: AHashMap<AtomView<'_>, symbolica::domains::float::Complex<f64>> =
//     //     data_atom_map
//     //         .0
//     //         .iter()
//     //         .zip(data_atom_map.1.iter())
//     //         .map(|(k, v)| (k.as_view(), (*v).into()))
//     //         .collect();

//     // let mut constvec = AHashMap::new();

//     // for (k, v) in const_atom_map.iter() {
//     //     constvec.insert(Atom::var(*k), *v);
//     // }
//     // for (k, &v) in constvec.iter() {
//     //     const_map.insert(k.as_view(), v.into());
//     // }

//     // let mut replacements = vec![];
//     // let mut fn_map: FunctionMap<Rational> = FunctionMap::new();

//     // for (k, v) in const_atom_map.iter() {
//     //     let name_re = Atom::var(symbol!(k.to_string() + "_re"));
//     //     let name_im = Atom::var(symbol!(k.to_string() + "_im"));
//     //     let i = Atom::var(Atom::I);
//     //     let pat = &name_re + i * &name_im;
//     //     replacements.push(Replacement::new(
//     //         Atom::var(*k).to_pattern(),
//     //         pat.to_pattern(),
//     //     ));

//     //     fn_map.add_constant(name_re, Rational::from(v.re));
//     //     fn_map.add_constant(name_im, Rational::from(v.im));
//     // }

//     // let mut params = data_atom_map.0.clone();
//     // params.push(Atom::var(Atom::I));

//     // let mut truth_net = network.clone();

//     // let function_map = HashMap::new();
//     // truth_net.evaluate_complex(|i| i.into(), &const_map, &function_map);
//     // truth_net.contract().unwrap();
//     // let truth = truth_net
//     //     .result()
//     //     .unwrap()
//     //     .0
//     //     .scalar()
//     //     .unwrap()
//     //     .try_into_concrete()
//     //     .unwrap()
//     //     .try_into_complex()
//     //     .unwrap();
//     // let mut postcontracted = network.clone();
//     // postcontracted.contract().unwrap();
//     // assert_relative_eq!(
//     //     truth,
//     //     postcontracted
//     //         .result()
//     //         .unwrap()
//     //         .0
//     //         .scalar()
//     //         .unwrap()
//     //         .try_into_concrete()
//     //         .unwrap()
//     //         .try_into_complex()
//     //         .unwrap(),
//     //     epsilon = 0.1
//     // );

//     // let mut counting_network: Network<MixedTensor<_, SmartShadowStructure<_, _>>, Atom> =
//     //     network.clone().cast();

//     // // (&replacements);
//     // let mut values: Vec<SymComplex<f64>> = data_atom_map.1.iter().map(|c| (*c).into()).collect();
//     // values.push(SymComplex::from(Complex::i()));

//     // let mut postcontracted_eval_tree_tensor = counting_network
//     //     .clone()
//     //     .to_fully_parametric()
//     //     .eval_tree(&fn_map, &params)
//     //     .unwrap();

//     // postcontracted_eval_tree_tensor.horner_scheme();
//     // // postcontracted_eval_tree_tensor.common_pair_elimination();
//     // postcontracted_eval_tree_tensor.common_subexpression_elimination();

//     // let mut mapped_postcontracted_eval_tree_tensor =
//     //     postcontracted_eval_tree_tensor.map_coeff::<SymComplex<f64>, _>(&|r| r.into());

//     // let mut out = mapped_postcontracted_eval_tree_tensor.evaluate(&values);
//     // out.contract().unwrap();
//     // assert_relative_eq!(
//     //     truth,
//     //     out.result().unwrap().0.scalar().unwrap().into(),
//     //     epsilon = 0.1
//     // );

//     // let mut levels: Levels<_, _, Atom> = counting_network.clone().into();
//     // let mut levels2 = levels.clone();

//     // let mut eval_tree_leveled_tensor = levels
//     //     .contract(1, &mut fn_map)
//     //     .to_evaluation_tree(&fn_map, &params)
//     //     .unwrap();

//     // eval_tree_leveled_tensor.horner_scheme();
//     // eval_tree_leveled_tensor.common_subexpression_elimination();
//     // // evaluator_tensor.common_pair_elimination();

//     // let mut eval_tree_leveled_tensor_depth2 = levels2
//     //     .contract(2, &mut fn_map)
//     //     .to_evaluation_tree(&fn_map, &params)
//     //     .unwrap();

//     // eval_tree_leveled_tensor_depth2.horner_scheme();
//     // eval_tree_leveled_tensor_depth2.common_subexpression_elimination();
//     // // eval_tree_leveled_tensor_depth2.common_pair_elimination();
//     // // evaluator_tensor.evaluate(&values);
//     // let mut neet = eval_tree_leveled_tensor.map_coeff::<SymComplex<f64>, _>(&|r| r.into());

//     // let mut neet2 = eval_tree_leveled_tensor_depth2.map_coeff::<SymComplex<f64>, _>(&|r| r.into());

//     // let out = neet.evaluate(&values).scalar().unwrap();
//     // assert!(truth.relative_eq(&out.into(), 0.1, 1.));

//     // let out = neet2.evaluate(&values).scalar().unwrap();
//     // assert!(truth.relative_eq(&out.into(), 0.1, 1.));
//     // network.contract().unwrap();

//     // let mut precontracted = network.clone();
//     // precontracted.evaluate_complex(|i| i.into(), &const_map, &function_map);
//     // assert!(truth.relative_eq(
//     //     &precontracted
//     //         .result()
//     //         .unwrap()
//     //         .0
//     //         .scalar()
//     //         .unwrap()
//     //         .try_into_concrete()
//     //         .unwrap()
//     //         .try_into_complex()
//     //         .unwrap(),
//     //     0.1,
//     //     1.
//     // ));

//     // let mut contracted_counting_network = counting_network.clone();
//     // contracted_counting_network.contract().unwrap();

//     // let mut precontracted_eval_tree_net = contracted_counting_network
//     //     .clone()
//     //     .to_fully_parametric()
//     //     .eval_tree(&fn_map, &params)
//     //     .unwrap();
//     // precontracted_eval_tree_net.horner_scheme();
//     // precontracted_eval_tree_net.common_subexpression_elimination();
//     // // precontracted_eval_tree_net.common_pair_elimination();

//     // let mut mapped_precontracted_eval_tree_net =
//     //     precontracted_eval_tree_net.map_coeff::<SymComplex<f64>, _>(&|r| r.into());

//     // let out = mapped_precontracted_eval_tree_net.evaluate(&values);
//     // assert!(truth.relative_eq(&out.result().unwrap().0.scalar().unwrap().into(), 0.1, 1.));

//     // let mut mapped_precontracted_eval_net = mapped_precontracted_eval_tree_net.linearize(Some(1));

//     // let out = mapped_precontracted_eval_net.evaluate(&values);
//     // assert!(truth.relative_eq(&out.result().unwrap().0.scalar().unwrap().into(), 0.1, 1.));

//     // let mut neeet = precontracted_eval_tree_net
//     //     .map_coeff::<f64, _>(&|r| r.into())
//     //     .linearize(None)
//     //     .export_cpp(
//     //         "nested_evaluation_asm",
//     //         "nested_evaluation_asm",
//     //         true,
//     //         InlineASM::X64,
//     //     )
//     //     .unwrap()
//     //     .compile_and_load("nested_evaluation_asm", CompileOptions::default())
//     //     .unwrap();

//     // let out = neeet.evaluate_complex(&values);
//     // assert!(truth.relative_eq(&(out.result().unwrap().0.scalar().unwrap()).into(), 0.1, 1.),);
// }
