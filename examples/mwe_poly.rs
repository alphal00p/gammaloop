use itertools::Itertools;
use ref_ops::RefNeg;
use std::fs::File;
use std::io::BufReader;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use symbolica::atom::{Atom, AtomView};
use symbolica::evaluate::FunctionMap;
use symbolica::id::Replacement;
use symbolica::poly::Variable;
use symbolica::state::{State, Workspace};
use symbolica::{self, fun, symb};

fn main() {
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
    let f = File::open("examples/data/6photons_1L.dat").unwrap();
    let mut reader = brotli::Decompressor::new(BufReader::new(f), 4096);

    let a = Atom::import(&mut reader, None).unwrap();

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
    let mut fn_map_shadowed = FunctionMap::new();
    let mut fn_map = FunctionMap::new();

    let reps = Arc::try_unwrap(reps).unwrap().into_inner().unwrap();

    fn_map.add_constant(Atom::parse("MT").unwrap(), 1.into());

    fn_map_shadowed.add_constant(Atom::parse("MT").unwrap(), 1.into());
    fn_map.add_constant(Atom::parse("ee").unwrap(), 2.into());
    fn_map_shadowed.add_constant(Atom::parse("ee").unwrap(), 2.into());
    for (v, k) in reps.iter().enumerate() {
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
    let epsb = symb!("ϵbar");

    let mut params = (0..12)
        .map(|i| {
            cinds
                .iter()
                .map(move |c| fun!(q, Atom::new_num(i), c.as_view()))
        })
        .flatten()
        .collect_vec();

    params.extend(
        [0, 2]
            .iter()
            .map(|&i| cinds.iter().map(move |c| fun!(eps, Atom::new_num(i), c)))
            .flatten()
            .chain(
                [5, 7, 9, 11]
                    .iter()
                    .map(|&i| cinds.iter().map(move |c| fun!(epsb, Atom::new_num(i), c)))
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
    atom_eval.linearize(Some(2));

    let elapsed = time.elapsed();
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
    poly_eval.linearize(Some(2));

    let elapsed = time.elapsed();
    println!("Poly: Linearized {elapsed:?}");
}
