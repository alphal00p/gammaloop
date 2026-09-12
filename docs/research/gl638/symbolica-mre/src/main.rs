use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    domains::float::Complex,
    parse,
};

fn main() {
    // Use Symbolica's normal license initialization; no application/OEM setup.
    let params: Vec<_> = [
        "f", "dL", "dR", "rL", "sL", "uL", "vL", "hL", "rR", "sR", "uR", "vR", "hR",
    ]
    .into_iter()
    .map(|s| parse!(s))
    .collect();
    let common = parse!("f/(dL*dR*rL^2*rR^2)");
    let left = parse!("uL/(rL-sL)+vL/(-rL-sL)");
    let right = parse!("uR/(rR-sR)+vR/(-rR-sR)");
    let pi = Atom::var(Symbol::PI);
    let expressions = [
        &common * &left * &right,
        Atom::i() * &pi * &common * &left * parse!("hR"),
        -Atom::i() * &pi * &common * parse!("hL") * &right,
        &pi * &pi * &common * parse!("hL*hR"),
    ];
    let mut input: Vec<_> = [1., 1., 1., 2., 1., 1., 1., 0.3, 3., 1., 1., 1., 0.4]
        .into_iter()
        .map(|x| Complex::new(x, 0.0))
        .collect();
    input[0] = Complex::new(0.0, 1.0);
    let pi = std::f64::consts::PI;
    let expected = [
        Complex::new(0.0, 1.0 / 216.0),
        Complex::new(-pi / 135.0, 0.0),
        Complex::new(pi / 480.0, 0.0),
        Complex::new(0.0, pi * pi / 300.0),
    ];
    let mut all_match = true;
    for expanded in [false, true] {
        let expressions: Vec<_> = expressions
            .iter()
            .map(|a| if expanded { a.expand() } else { a.clone() })
            .collect();
        let evaluators: Vec<_> = expressions
            .iter()
            .map(|a| a.evaluator(&params).build().unwrap())
            .collect();
        let separate: Vec<_> = evaluators
            .iter()
            .cloned()
            .map(|e| {
                e.map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()))
                    .evaluate_single(&input)
            })
            .collect();
        let mut merged = evaluators[0].clone();
        for evaluator in evaluators.iter().skip(1) {
            merged.merge(evaluator.clone(), None).unwrap();
        }
        let mut merged_output = [Complex::new(0.0, 0.0); 4];
        merged
            .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()))
            .evaluate(&input, &mut merged_output);
        let mut joint_output = [Complex::new(0.0, 0.0); 4];
        Atom::evaluator_multiple(&expressions, &params)
            .build()
            .unwrap()
            .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()))
            .evaluate(&input, &mut joint_output);
        println!("expanded={expanded}");
        println!("expected: {expected:?}");
        for (route, output) in [
            ("separate", separate.as_slice()),
            ("joint", joint_output.as_slice()),
            ("merged", merged_output.as_slice()),
        ] {
            println!("{route}: {output:?}");
            all_match &= output.iter().zip(&expected).all(|(actual, expected)| {
                (actual.re - expected.re).abs() < 1e-12 && (actual.im - expected.im).abs() < 1e-12
            });
        }
    }
    assert!(
        all_match,
        "evaluator output differs from the analytical value"
    );
}
