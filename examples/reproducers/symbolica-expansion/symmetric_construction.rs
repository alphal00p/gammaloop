use std::{
    hint::black_box,
    sync::{Arc, Mutex},
    time::Instant,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, EvaluationInfo, FunctionBuilder, Symbol},
    domains::float::{Complex, Float},
    parser::ParseSettings,
};

fn parse(source: &str) -> Atom {
    Atom::parse(source, "symmetry_probe", ParseSettings::symbolica()).unwrap()
}

fn build(head: Symbol, args: &[Atom]) -> Atom {
    args.iter()
        .fold(FunctionBuilder::new(head), |builder, arg| {
            builder.add_arg(arg)
        })
        .finish()
}

fn record(name: &str, value: &Atom) {
    let bytes = value
        .as_view()
        .get_data()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    println!("GATE {name} {} {bytes}", value.to_canonical_string());
}

fn gates() {
    let symmetric = symbolica::symbol!("symmetry_probe::f"; Symmetric);
    let antisymmetric = symbolica::symbol!("symmetry_probe::a"; Antisymmetric);
    let cyclic = symbolica::symbol!("symmetry_probe::c"; Cyclesymmetric);
    let linear = symbolica::symbol!("symmetry_probe::l"; Symmetric, Linear);
    let mut args: Vec<_> = [
        "-1",
        "3/7",
        "x",
        "x^2",
        "x*y",
        "x+y",
        "p(mink(4))",
        "q(a,b)",
        "w_",
        "1.25",
        "1.25`8",
        "1.25`30",
    ]
    .iter()
    .map(|value| parse(value))
    .collect();
    // Function arguments use the inherent comparator, not the raw-byte Ord API.
    assert!(args.iter().any(|a| args.iter().any(|b| {
        let (a, b) = (a.as_view(), b.as_view());
        AtomView::cmp(&a, &b) != <AtomView<'_> as Ord>::cmp(&a, &b)
    })));
    println!("CHECK distinct_atom_orderings true");
    args.sort_by(|a, b| a.as_view().cmp(&b.as_view()));
    let reference = build(symmetric, &args);
    for shift in 0..args.len() {
        let mut rotated = args.clone();
        rotated.rotate_left(shift);
        assert_eq!(build(symmetric, &rotated), reference);
        rotated.reverse();
        assert_eq!(build(symmetric, &rotated), reference);
    }
    for (name, head) in [
        ("symmetric", symmetric),
        ("antisymmetric", antisymmetric),
        ("cyclic", cyclic),
        ("linear", linear),
    ] {
        for source in [
            vec![],
            vec![parse("x")],
            vec![parse("x"), parse("x")],
            vec![parse("y"), parse("x")],
            vec![parse("x+y"), parse("2*x")],
            args.clone(),
        ] {
            record(name, &build(head, &source));
        }
    }
    for arity in [19, 20, 21, 128, 257] {
        let values: Vec<_> = (0..arity).map(|i| args[i % args.len()].clone()).collect();
        let expected = build(symmetric, &values);
        let mut reversed = values.clone();
        reversed.reverse();
        assert_eq!(build(symmetric, &reversed), expected);
        record(&format!("arity_{arity}"), &expected);
    }
    assert_eq!(
        build(antisymmetric, &[parse("x"), parse("x")]),
        Atom::num(0)
    );
    assert_eq!(
        build(antisymmetric, &[parse("y"), parse("x")]),
        -build(antisymmetric, &[parse("x"), parse("y")])
    );

    let calls = Arc::new(Mutex::new(Vec::<Atom>::new()));
    let observed = Arc::clone(&calls);
    let callback = symbolica::symbol!("symmetry_probe::callback"; Symmetric; norm = move |view, out| {
        observed.lock().unwrap().push(view.to_owned());
        let AtomView::Fun(fun) = view else { unreachable!() };
        if fun.get_nargs() == 2 && fun.iter().all(|arg| arg == parse("x").as_view()) {
            **out = Atom::num(17);
        }
    });
    let mut sorted = vec![parse("p(mink(4))"), parse("q(mink(4))")];
    sorted.sort_by(|a, b| a.as_view().cmp(&b.as_view()));
    let canonical = build(callback, &sorted);
    calls.lock().unwrap().clear();
    let same = build(callback, &sorted);
    sorted.reverse();
    let reversed = build(callback, &sorted);
    assert_eq!(same, reversed);
    assert_eq!(*calls.lock().unwrap(), vec![canonical.clone(), canonical]);
    assert_eq!(build(callback, &[parse("x"), parse("x")]), Atom::num(17));
    println!("CHECK callback_order_count_and_replacement true");

    let transcript = Arc::new(Mutex::new(Vec::<String>::new()));
    let normalizations = Arc::clone(&transcript);
    let evaluations = Arc::clone(&transcript);
    let evaluated = symbolica::symbol!("symmetry_probe::evaluated"; Symmetric;
        norm = move |view, _| normalizations.lock().unwrap().push(format!("normalize:{:?}", view.get_data())),
        eval = EvaluationInfo::new().register(move |args: &[Complex<Float>]| {
            evaluations.lock().unwrap().push(format!("evaluate:{args:?}"));
            args[0].clone()
        })
    );
    let small = parse("1.25`8");
    let large = parse("2.5`30");
    let first = build(evaluated, &[large.clone(), small.clone()]);
    let second = build(evaluated, &[small.clone(), large]);
    assert_eq!(first, small);
    assert_eq!(first, second);
    let observed = transcript.lock().unwrap();
    assert_eq!(observed.len(), 4);
    assert_eq!(observed[0], observed[2]);
    assert_eq!(observed[1], observed[3]);
    println!("CHECK numeric_callback_transcript {observed:?}");
    record("evaluated", &first);
}

fn main() {
    gates();
    if std::env::args().nth(1).as_deref() != Some("bench") {
        return;
    }
    let head = symbolica::symbol!("symmetry_probe::metric"; Symmetric);
    let count: usize = std::env::var("SYMMETRY_PROBE_ITERATIONS")
        .ok()
        .map(|value| value.parse().unwrap())
        .unwrap_or(50_000);
    for (name, arity, reverse, wide) in [
        ("ordered_pair", 2, false, false),
        ("reversed_pair", 2, true, false),
        ("ordered_eight", 8, false, false),
        ("reversed_eight", 8, true, false),
        ("reversed_thirty_two", 32, true, false),
        ("wide_pair", 2, false, true),
    ] {
        let mut args: Vec<_> = (0..arity)
            .map(|i| {
                if wide {
                    parse(&format!(
                        "p{i}({})",
                        (0..64)
                            .map(|j| format!("x{j}*y{j}"))
                            .collect::<Vec<_>>()
                            .join("+")
                    ))
                } else {
                    parse(&format!("p{i}(mink(4))"))
                }
            })
            .collect();
        args.sort_by(|a, b| a.as_view().cmp(&b.as_view()));
        if reverse {
            args.reverse();
        }
        let expected = build(head, &args);
        let mut samples = Vec::new();
        for _ in 0..5 {
            let start = Instant::now();
            for _ in 0..count {
                drop(black_box(build(black_box(head), black_box(&args))));
            }
            samples.push(start.elapsed().as_nanos());
            assert_eq!(build(head, &args), expected);
        }
        println!("BENCH {name} {count} {samples:?}");
    }
}
