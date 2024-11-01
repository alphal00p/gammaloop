use itertools::Itertools;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use ref_ops::RefNeg;
use std::fs::File;
use std::io::BufReader;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use symbolica::atom::{Atom, AtomView};
use symbolica::domains::float::{Complex, Real};
use symbolica::evaluate::{CompileOptions, ExpressionEvaluator, FunctionMap};
use symbolica::id::{Pattern, Replacement};
use symbolica::poly::Variable;
use symbolica::state::{State, Workspace};
use symbolica::{self, fun, symb};

fn main() {
    let name = "sixphotons_1L";

    fn to_shadowed_poly_impl(
        poly: &symbolica::poly::polynomial::MultivariatePolynomial<
            symbolica::domains::atom::AtomField,
            u8,
        >,
        workspace: &Workspace,
        reps: Arc<Mutex<Vec<Atom>>>,
    ) -> Atom {
        if poly.is_zero() {
            return Atom::new_num(0);
        }

        let mut add = Atom::new_num(0);
        let coef = symb!("coef");
        let shift = reps.as_ref().lock().unwrap().len();

        let mut mul_h = workspace.new_atom().into_inner();
        let mut var_h = workspace.new_atom();
        let mut num_h = workspace.new_atom();
        let mut pow_h = workspace.new_atom();

        for (i, monomial) in poly.into_iter().enumerate() {
            mul_h = Atom::new_num(1);
            for (var_id, &pow) in poly.variables.iter().zip(monomial.exponents) {
                if pow > 0 {
                    match var_id {
                        Variable::Symbol(v) => {
                            var_h.to_var(*v);
                        }
                        Variable::Temporary(_) => {
                            unreachable!("Temporary variable in expression")
                        }
                        Variable::Function(_, a) | Variable::Other(a) => {
                            var_h.set_from_view(&a.as_view());
                        }
                    }

                    if pow > 0 {
                        num_h.to_num((pow as i64).into());
                        pow_h.to_pow(var_h.as_view(), num_h.as_view());
                        mul_h = mul_h * pow_h.as_view();
                    } else {
                        mul_h = mul_h * var_h.as_view();
                    }
                }
            }

            reps.lock()
                .as_mut()
                .unwrap()
                .push(monomial.coefficient.clone());

            mul_h = mul_h * fun!(coef, Atom::new_num((i as usize + shift) as i64));
            add = add + mul_h.as_view();
        }

        add
    }
    fn shadow_poly(
        poly: symbolica::poly::polynomial::MultivariatePolynomial<
            symbolica::domains::atom::AtomField,
            u8,
        >,
        reps: Arc<Mutex<Vec<Atom>>>,
    ) -> Atom {
        Workspace::get_local().with(|ws| to_shadowed_poly_impl(&poly, ws, reps))
    }
    let f = File::open(&format!("examples/data/{}.dat", name)).unwrap();
    let mut reader = brotli::Decompressor::new(BufReader::new(f), 4096);

    let a = Atom::import(&mut reader, None).unwrap();

    println!("{}", a);

    let var_map: Arc<Vec<Variable>> = Arc::new(
        [
            "Q(0,cind(0))",
            "Q(1,cind(0))",
            "Q(2,cind(0))",
            "Q(3,cind(0))",
            "Q(4,cind(0))",
            "Q(5,cind(0))",
            "Q(6,cind(0))",
            "Q(7,cind(0))",
            "Q(8,cind(0))",
            "Q(9,cind(0))",
            "Q(10,cind(0))",
            "Q(11,cind(0))",
        ]
        .into_iter()
        .map(|a| Variable::from(Atom::parse(a).unwrap()))
        .collect(),
    );

    let poly = a.as_view().to_polynomial_in_vars::<u8>(&var_map);
    let reps = Arc::new(Mutex::new(Vec::new()));

    let shadowed = shadow_poly(poly, reps.clone());

    println!("{shadowed}");
    let mut fn_map_shadowed = FunctionMap::new();
    let mut fn_map = FunctionMap::new();

    let coefs = Arc::try_unwrap(reps).unwrap().into_inner().unwrap();

    fn_map.add_constant(Atom::parse("MT").unwrap(), 1.into());

    fn_map_shadowed.add_constant(Atom::parse("MT").unwrap(), 1.into());
    fn_map.add_constant(Atom::parse("ee").unwrap(), 2.into());
    fn_map_shadowed.add_constant(Atom::parse("ee").unwrap(), 2.into());
    for (v, k) in coefs.iter().enumerate() {
        fn_map_shadowed
            .add_tagged_function(
                symb!("coef"),
                vec![Atom::new_num(v as i64).into()],
                format!("coef{v}"),
                vec![],
                k.as_view(),
            )
            .unwrap();
    }

    let q = symb!("Q");
    let cind = symb!("cind");
    let cinds = [
        fun!(cind, Atom::new_num(0)),
        fun!(cind, Atom::new_num(1)),
        fun!(cind, Atom::new_num(2)),
        fun!(cind, Atom::new_num(3)),
    ];

    let eps = symb!("ϵ");
    let eps_inds: Vec<i32> = vec![0, 2];
    let epsb = symb!("ϵbar");

    let epsb_inds: Vec<i32> = vec![5, 7, 9, 11];
    let u = symb!("u");

    let u_inds: Vec<i32> = vec![];
    let ubar = symb!("ubar");
    let ubar_inds: Vec<i32> = vec![];

    let mut params = (0..12)
        .map(|i| {
            cinds
                .iter()
                .map(move |c| fun!(q, Atom::new_num(i), c.as_view()))
        })
        .flatten()
        .collect_vec();

    params.extend(
        eps_inds
            .iter()
            .map(|&i| cinds.iter().map(move |c| fun!(eps, Atom::new_num(i), c)))
            .flatten()
            .chain(
                epsb_inds
                    .iter()
                    .map(|&i| cinds.iter().map(move |c| fun!(epsb, Atom::new_num(i), c)))
                    .flatten(),
            )
            .chain(
                u_inds
                    .iter()
                    .map(|&i| cinds.iter().map(move |c| fun!(u, Atom::new_num(i), c)))
                    .flatten(),
            )
            .chain(
                ubar_inds
                    .iter()
                    .map(|&i| cinds.iter().map(move |c| fun!(ubar, Atom::new_num(i), c)))
                    .flatten(),
            ),
    );

    params.push(Atom::new_var(State::I));

    let reps = (0..12)
        .map(|i| {
            (
                fun!(q, Atom::new_num(i), cinds[0]).into_pattern(),
                fun!(q, Atom::new_num(i), cinds[0])
                    .ref_neg()
                    .into_pattern()
                    .into(),
            )
        })
        .collect_vec();

    let signatures = vec![
        vec![0, 1, 2],
        vec![2, 3, 6, 7],
        vec![0, 6, 7, 9],
        vec![0, 1, 2, 3],
        vec![3, 4, 7],
        vec![0, 1, 2, 3, 4, 5],
        vec![0, 1, 2, 3, 4, 5, 6, 7],
        vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        vec![0, 1, 5, 6, 7, 8, 9, 10, 11],
        vec![0, 1, 2, 4, 5, 6, 7, 8, 11],
        vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        vec![0, 1, 2, 6, 7, 8, 9, 11],
        vec![0, 1, 3, 4, 5, 6, 7],
        vec![0, 1, 2, 5, 6, 7, 8, 11],
        vec![0, 1, 2, 3, 4, 11],
    ];

    let mut polys = Vec::new();
    let mut atoms = Vec::new();

    for s in signatures {
        let reps: Vec<_> = s
            .iter()
            .map(|&i| Replacement::new(&reps[i].0, &reps[i].1))
            .collect();

        polys.push(shadowed.replace_all_multiple(&reps));
        atoms.push(a.replace_all_multiple(&reps));
    }

    let poly_views: Vec<_> = polys.iter().map(|p| p.as_view()).collect();
    let atom_views: Vec<_> = atoms.iter().map(|a| a.as_view()).collect();

    let mut poly_views_iter = poly_views.iter();
    let poly_view = poly_views_iter.next().unwrap();
    let time = Instant::now();
    let mut iterative_poly_eval = poly_view
        .to_evaluation_tree(&fn_map_shadowed, &params)
        .unwrap();

    let elapsed = time.elapsed();
    println!("Poly: To eval single {elapsed:?}");
    let time = Instant::now();
    iterative_poly_eval.common_subexpression_elimination();
    let elapsed = time.elapsed();
    println!("Poly: CSE single {elapsed:?}");
    let time = Instant::now();
    iterative_poly_eval.horner_scheme();
    let elapsed = time.elapsed();
    println!("Poly: horner single {elapsed:?}");
    let time = Instant::now();
    let mut iterative_poly_eval_lin = iterative_poly_eval.linearize(Some(2));
    let elapsed = time.elapsed();
    println!("Poly: Linearize single {elapsed:?}");
    let time = Instant::now();

    for i in poly_views_iter {
        let time = Instant::now();
        let mut poly_eval = i.to_evaluation_tree(&fn_map_shadowed, &params).unwrap();
        poly_eval.common_subexpression_elimination();
        poly_eval.horner_scheme();
        let poly_eval_lin = poly_eval.linearize(Some(2));
        iterative_poly_eval_lin
            .merge(poly_eval_lin, Some(2))
            .unwrap();
        let elapsed = time.elapsed();
        println!("Poly: merge {elapsed:?}");
    }

    let elapsed = time.elapsed();
    println!("Poly: Push optimize {elapsed:?}");

    let mut flat_views_iter = atom_views.iter();
    let flat_view = flat_views_iter.next().unwrap();
    let time = Instant::now();
    let mut iterative_flat_eval = flat_view
        .to_evaluation_tree(&fn_map_shadowed, &params)
        .unwrap();

    let elapsed = time.elapsed();
    println!("To eval single {elapsed:?}");
    let time = Instant::now();
    iterative_flat_eval.common_subexpression_elimination();
    let elapsed = time.elapsed();
    println!("CSE single {elapsed:?}");
    let time = Instant::now();
    iterative_flat_eval.horner_scheme();
    let elapsed = time.elapsed();
    println!("horner single {elapsed:?}");
    let time = Instant::now();
    let mut iterative_flat_eval_lin = iterative_flat_eval.linearize(Some(2));
    let elapsed = time.elapsed();
    println!("Linearize single {elapsed:?}");
    let time = Instant::now();

    for i in flat_views_iter {
        let time = Instant::now();
        let mut flat_eval = i.to_evaluation_tree(&fn_map_shadowed, &params).unwrap();
        flat_eval.common_subexpression_elimination();
        flat_eval.horner_scheme();
        let flat_eval_lin = flat_eval.linearize(Some(2));
        iterative_flat_eval_lin
            .merge(flat_eval_lin, Some(2))
            .unwrap();
        let elapsed = time.elapsed();
        println!("merge {elapsed:?}");
    }

    let elapsed = time.elapsed();
    println!("Push optimize {elapsed:?}");

    let time = Instant::now();
    let mut atom_eval = AtomView::to_eval_tree_multiple(&atom_views, &fn_map, &params).unwrap();
    let elapsed = time.elapsed();
    println!("To eval {elapsed:?}");

    let time = Instant::now();
    atom_eval.common_subexpression_elimination();
    let elapsed = time.elapsed();
    println!("CSE {elapsed:?}");
    let time = Instant::now();
    atom_eval.horner_scheme();
    let elapsed = time.elapsed();
    println!("horner {elapsed:?}");

    let time = Instant::now();
    let mut lin_eval: ExpressionEvaluator<symbolica::domains::float::Complex<f64>> = atom_eval
        .clone()
        .linearize(Some(2))
        .map_coeff(&|a| a.into());

    let elapsed = time.elapsed();

    let mut lin_eval_comp = atom_eval
        .clone()
        .linearize(Some(2))
        .export_cpp(
            &format!("{}.cpp", name),
            name,
            true,
            symbolica::evaluate::InlineASM::None,
        )
        .unwrap()
        .compile(&format!("{}.so", name), CompileOptions::default())
        .unwrap()
        .load()
        .unwrap();
    let mut lin_eval_comp_asm = atom_eval
        .linearize(Some(2))
        .export_cpp(
            &format!("{}_asm.cpp", name),
            name,
            true,
            symbolica::evaluate::InlineASM::X64,
        )
        .unwrap()
        .compile(&format!("{}_asm.so", name), CompileOptions::default())
        .unwrap()
        .load()
        .unwrap();
    println!("Linearized {elapsed:?}");
    let time = Instant::now();

    let mut poly_eval =
        AtomView::to_eval_tree_multiple(&poly_views, &fn_map_shadowed, &params).unwrap();
    let elapsed = time.elapsed();
    println!("Poly: To eval {elapsed:?}");

    let time = Instant::now();
    poly_eval.common_subexpression_elimination();
    let elapsed = time.elapsed();
    println!("Poly: CSE {elapsed:?}");
    let time = Instant::now();
    poly_eval.horner_scheme();
    let elapsed = time.elapsed();
    println!("Poly: horner {elapsed:?}");
    let time = Instant::now();
    let mut lin_poly_eval: ExpressionEvaluator<symbolica::domains::float::Complex<f64>> = poly_eval
        .clone()
        .linearize(Some(2))
        .map_coeff(&|a| a.into());

    let elapsed = time.elapsed();

    let mut lin_poly_eval_comp = poly_eval
        .clone()
        .linearize(Some(2))
        .export_cpp(
            &format!("{}-poly.cpp", name),
            name,
            true,
            symbolica::evaluate::InlineASM::None,
        )
        .unwrap()
        .compile(&format!("{}-poly.so", name), CompileOptions::default())
        .unwrap()
        .load()
        .unwrap();

    let mut lin_poly_eval_comp_asm = poly_eval
        .clone()
        .linearize(Some(2))
        .export_cpp(
            &format!("{}-poly_asm.cpp", name),
            name,
            true,
            symbolica::evaluate::InlineASM::X64,
        )
        .unwrap()
        .compile(&format!("{}-poly_asm.so", name), CompileOptions::default())
        .unwrap()
        .load()
        .unwrap();

    println!("Poly: Linearized {elapsed:?}");

    let mut out_poly = vec![Complex::new(0.0, 0.0); poly_views.len()];
    let mut rng = SmallRng::seed_from_u64(4);
    let params = (0..params.len())
        .map(|i| Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
        .collect::<Vec<_>>();

    let time = Instant::now();
    lin_poly_eval.evaluate(&params, &mut out_poly);
    let elapsed = time.elapsed();

    println!("Poly: eval {elapsed:?}");

    let mut out_poly_comp = out_poly.clone();
    let time = Instant::now();
    lin_poly_eval_comp.evaluate(&params, &mut out_poly_comp);
    let elapsed = time.elapsed();

    println!("Poly: compiled eval {elapsed:?}");

    let mut out_poly_comp_asm = out_poly.clone();
    let time = Instant::now();
    lin_poly_eval_comp_asm.evaluate(&params, &mut out_poly_comp_asm);
    let elapsed = time.elapsed();

    println!("Poly: compiled eval {elapsed:?}");

    let mut out = vec![Complex::new(0.0, 0.0); poly_views.len()];
    let time = Instant::now();
    lin_eval.evaluate(&params, &mut out);
    let elapsed = time.elapsed();

    println!("eval {elapsed:?}");

    let mut out_comp = out_poly.clone();
    let time = Instant::now();
    lin_eval_comp.evaluate(&params, &mut out_comp);
    let elapsed = time.elapsed();

    println!("compiled eval {elapsed:?}");

    let mut out_comp_asm = out_poly.clone();
    let time = Instant::now();
    lin_eval_comp_asm.evaluate(&params, &mut out_comp_asm);
    let elapsed = time.elapsed();

    println!("compiled eval {elapsed:?}");

    let coef = symb!("coef");

    let coefs_syms: Vec<_> = (0..coefs.len())
        .map(|i| fun!(coef, Atom::new_num(i as i64)).into_pattern())
        .collect();

    let coefs_reps: Vec<_> = coefs.iter().map(|a| a.into_pattern().into()).collect();

    let mut reps: Vec<_> = coefs_syms
        .iter()
        .zip(coefs_reps.iter())
        .map(|(p, rhs)| Replacement::new(p, rhs))
        .collect();

    let pats = [
        (
            &Pattern::parse("ϵbar(n_, cind(x_))").unwrap(),
            &Pattern::parse("n_+4+x_").unwrap().into(),
        ),
        (
            &Pattern::parse("ϵ(n_, cind(x_))").unwrap(),
            &Pattern::parse("n_+4+x_").unwrap().into(),
        ),
        (
            &Pattern::parse("Q(n_, cind(x_))").unwrap(),
            &Pattern::parse("n_*x_").unwrap().into(),
        ),
    ];

    let other_reps: Vec<_> = pats
        .iter()
        .map(|(p, rhs)| Replacement::new(p, rhs))
        .collect();

    assert_eq!(
        Atom::new_num(0),
        (&atoms[0] - &polys[0].replace_all_multiple(&reps)).replace_all_multiple(&other_reps)
    );

    for (a, b) in out.iter().zip(out_poly.iter()) {
        assert!((a - b).norm().re < 1e-10, "{:?} {:?}", a, b);
    }

    for (a, b) in out_comp.iter().zip(out_poly_comp.iter()) {
        assert!((a - b).norm().re < 1e-10, "{:?} {:?}", a, b);
    }

    for (a, b) in out_comp_asm.iter().zip(out_poly_comp_asm.iter()) {
        assert!((a - b).norm().re < 1e-10, "{:?} {:?}", a, b);
    }
}

fn time_step<F, R>(label: &str, func: F) -> R
where
    F: FnOnce() -> R,
{
    let start = Instant::now();
    let result = func();
    println!("{} took {:?}", label, start.elapsed());
    result
}

// fn optimize_view(
//     view: AtomView,
//     fn_map: &FunctionMap,
//     params: &[Atom],
//     label: &str,
// ) -> ExpressionEvaluator<Complex<f64>> {

// }

// fn optimize_and_linearize_with_steps(
//     views: &[AtomView],
//     fn_map: &FunctionMap,
//     params: &[Atom],
//     label: &str,
// ) -> ExpressionEvaluator<Complex<f64>> {
//     let mut iter = views.iter();

//     // Initialize with the first item
//     let mut eval_tree = time_step(&format!("{} eval tree single", label), || {
//         iter.next()
//             .unwrap()
//             .to_evaluation_tree(fn_map, params)
//             .unwrap()
//     });
//     time_step(&format!("{} cse single", label), || {
//         eval_tree.common_subexpression_elimination()
//     });
//     time_step(&format!("{} horner single", label), || {
//         eval_tree.horner_scheme()
//     });

//     let mut lin_eval: ExpressionEvaluator<Complex<f64>> =
//         time_step(&format!("{} linearize single", label), || {
//             eval_tree.linearize(Some(2)).map_coeff(&Complex::from)
//         });

//     // Merge remaining items with granular timing
//     for (index, view) in iter.enumerate() {
//         let mut eval = time_step(&format!("to_evaluation_tree for view {}", index), || {
//             view.to_evaluation_tree(fn_map, params).unwrap()
//         });
//         time_step(&format!("CSE for view {}", index), || {
//             eval.common_subexpression_elimination()
//         });
//         time_step(&format!("Horner scheme for view {}", index), || {
//             eval.horner_scheme()
//         });
//         let lin = time_step(&format!("linearize view {}", index), || {
//             eval.linearize(Some(2)).map_coeff(Into::into)
//         });
//         time_step(&format!("merge view {}", index), || {
//             lin_eval.merge(lin, Some(2)).unwrap()
//         });
//     }

//     lin_eval
// }
