use idenso::{
    dirac::GammaSimplifySettings,
    gamma, gamma5,
    tensor::{SymbolicTensor, aliases::AliasInterfaces},
};
use serde_json::json;
use spenso::{
    chain, mink, p, q,
    structure::{dimension::Dimension, partial::PartialStructure},
    trace,
};
use std::{hint::black_box, path::Path, sync::Arc, time::Instant};
use symbolica::{
    atom::{AliasedAtom, Atom},
    symbol,
};

type Tensor = SymbolicTensor<PartialStructure>;
type Aliased = Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>;
enum Input {
    Plain(Tensor),
    Aliased(Aliased),
}
struct Case {
    name: String,
    input: Input,
}
fn trace_word(n: usize, dim: &Dimension, axial: bool, paired: bool, suffix: usize) -> Atom {
    let mut word = if axial { vec![gamma5!()] } else { Vec::new() };
    for i in 0..n {
        let mu = if paired {
            if i % 2 == 0 {
                p!(mink!(dim.clone()))
            } else {
                q!(mink!(dim.clone()))
            }
        } else {
            mink!(
                dim.clone(),
                symbol!(&format!("trusted_benchmark::mu_{suffix}_{i}"))
            )
        };
        word.push(gamma!(mu));
    }
    trace!(idenso::bis!(4);word)
}
fn cases() -> Vec<Case> {
    let four: Dimension = 4usize.into();
    let d: Dimension = symbol!("trusted_benchmark::D").into();
    let mut cases = Vec::new();
    for (name, n, dim, axial, paired) in [
        ("trace4_free8", 8, &four, false, false),
        ("trace4_free12", 12, &four, false, false),
        ("trace4_axial12", 12, &four, true, false),
        ("tracen_free8", 8, &d, false, false),
        ("tracen_free10", 10, &d, false, false),
        ("tracen_alternating12", 12, &d, false, true),
    ] {
        cases.push(Case {
            name: name.into(),
            input: Input::Plain(Tensor::infer(trace_word(n, dim, axial, paired, 0)).unwrap()),
        });
    }
    let factors = [0, 1, 2, 3, 1, 0]
        .into_iter()
        .map(|i| gamma!(mink!(4, symbol!(&format!("trusted_benchmark::open_mu{i}")))))
        .collect::<Vec<_>>();
    let open = chain!(idenso::bis!(4,symbol!("trusted_benchmark::spin_left")),idenso::bis!(4,symbol!("trusted_benchmark::spin_right"));factors);
    cases.push(Case {
        name: "open_repeated6".into(),
        input: Input::Plain(Tensor::infer(open).unwrap()),
    });
    for count in [1, 16, 64] {
        let mut definitions = Vec::new();
        let mut terms = Vec::new();
        for i in 0..=count {
            let coefficient = Atom::var(symbol!(&format!("trusted_benchmark::c{i}")));
            let body = Tensor::infer(coefficient * trace_word(4, &four, false, true, i)).unwrap();
            let handle = body.alias_handle().unwrap();
            if i < count {
                terms.push(handle.expression().clone());
            }
            definitions.push((handle, body));
        }
        let root = Tensor::infer(Atom::add_many(terms))
            .unwrap()
            .with_aliases(definitions)
            .unwrap();
        cases.push(Case {
            name: format!("aliased_domains_{count}"),
            input: Input::Aliased(Arc::new(root)),
        });
    }
    cases
}
#[inline(never)]
fn run(case: &Case, mode: &str) -> Aliased {
    match &case.input {
        Input::Plain(input) if mode == "admission" => {
            Tensor::infer(black_box(input.expression()).clone())
                .unwrap()
                .simplify_gamma(GammaSimplifySettings::default())
                .unwrap()
        }
        Input::Plain(input) => black_box(input)
            .simplify_gamma(GammaSimplifySettings::default())
            .unwrap(),
        Input::Aliased(input) => black_box(input)
            .simplify_gamma(GammaSimplifySettings::default())
            .unwrap(),
    }
}
fn rendered(result: &Aliased) -> String {
    let value = result.resolved().unwrap();
    format!(
        "{}\n{:?}\n",
        value.expression().to_plain_string(),
        value.structure()
    )
}
fn main() {
    idenso::representations::initialize();
    let args: Vec<_> = std::env::args().collect();
    let output = Path::new(&args[1]);
    std::fs::create_dir_all(output).unwrap();
    let action = args.get(2).map(String::as_str).unwrap_or("measure");
    let reference = args.get(3).map(Path::new);
    for case in cases() {
        let modes: &[&str] = match &case.input {
            Input::Plain(_) => &["typed", "admission"],
            Input::Aliased(_) => &["aliased"],
        };
        for mode in modes {
            let result = run(&case, mode);
            let text = rendered(&result);
            let repeated = result
                .simplify_gamma(GammaSimplifySettings::default())
                .unwrap();
            assert_eq!(text, rendered(&repeated), "idempotence: {}", case.name);
            if let Some(reference) = reference {
                assert_eq!(
                    text,
                    std::fs::read_to_string(reference.join(format!("{}-{mode}.output", case.name)))
                        .unwrap(),
                    "baseline exact output: {}",
                    case.name
                );
            }
            std::fs::write(output.join(format!("{}-{mode}.output", case.name)), &text).unwrap();
            println!(
                "{}",
                json!({"kind":"check","case":case.name,"mode":mode,"exact_reference":reference.is_some(),"idempotent":true,"resolved_bytes":text.len(),"alias_count":result.aliases().unwrap().len()})
            );
            if action == "check" {
                continue;
            }
            let timer = Instant::now();
            drop(black_box(run(&case, mode)));
            let ns = timer.elapsed().as_nanos().max(1);
            let repetitions = (10_000_000 / ns).clamp(1, 2000) as usize;
            for sample in 0..5 {
                let timer = Instant::now();
                for _ in 0..repetitions {
                    drop(black_box(run(&case, mode)));
                }
                println!(
                    "{}",
                    json!({"kind":"timing","case":case.name,"mode":mode,"sample":sample,"iterations":repetitions,"ns":timer.elapsed().as_nanos()as f64/repetitions as f64})
                );
            }
        }
    }
}
