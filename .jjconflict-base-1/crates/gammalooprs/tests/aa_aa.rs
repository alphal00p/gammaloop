use std::time::Instant;

use gammalooprs::LogMessage;
use gammalooprs::graph::parse::IntoGraph;
use gammalooprs::integrands::process::evaluators::EvaluatorStack;
use gammalooprs::processes::{ContractionMode, EvaluatorSettings, ExecutionMode};
use gammalooprs::utils::{F, GS};
use gammalooprs::{dot, graph::Graph, uv::UltravioletGraph};
use idenso::dirac::GammaSimplifier;
use idenso::shorthands::metric::MetricSimplifier;
use idenso::shorthands::schoonschip::Schoonschip;
use linnet::half_edge::subgraph::SubSetOps;
use rand::Rng;
use spenso::algebra::complex::Complex;
use spenso::shadowing::TensorCollectExt;
use spenso::shadowing::symbolica_utils::LogPrint;
use spenso::structure::concrete_index::ExpandedIndex;
use spenso::structure::representation::Minkowski;
use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::function;
use tabled::{
    builder::Builder,
    settings::{Alignment, Modify, Style, object::Columns},
};
#[test]
fn aaa() {
    let graph: Graph = dot!(digraph GL0{
                num = "1";
                overall_factor = "NumeratorDependentGrouping(0,1,(AutG(1))^(-1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(1,1,(AutG(1))^(-1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))";
                overall_factor_evaluated = "-2";
                projector = "ϵ(0,spenso::mink(4,hedge(0)))*ϵ(1,spenso::mink(4,hedge(1)))*ϵbar(2,spenso::mink(4,hedge(2)))*ϵbar(3,spenso::mink(4,hedge(3)))";


            exte0	 [style=invis];
            exte0	-> 3:0	 [id=0 dir=none particle="a" pin="x:@-left"];
            exte1	 [style=invis];
            exte1	-> 1:1	 [id=1 particle="a" pin="x:@-left"];
            exte2	 [style=invis];
            0:2	-> exte2	 [id=2 dir=none particle="a" pin="x:@+right"];
            exte3	 [style=invis];
            2:3	-> exte3	 [id=3 dir=none  particle="a" pin="x:@+right"];
            0:4	-> 2:5	 [id=4   particle="t"];

            3:6	-> a [id=5  particle="t"];
            a -> 0:7	 [particle="t"];
            2:8	-> b [particle="t"];
            a->b[particle = "g"]

            b->1:9	 [id=6  particle="t"];
            1:10	-> 3:11	 [id=7  particle="t"];
            }
    ).unwrap();

    let num = graph
        .numerator(&graph.full_filter(), &graph.empty_subgraph())
        .state
        .expr
        * &graph.global_prefactor.projector;
    let loop_edges = graph.full_filter().subtract(&graph.tree_edges);
    let mut reps = vec![];
    for (_, e, _) in graph.iter_edges_of(&loop_edges) {
        let rep = GS.split_mom_pattern_simple(e);
        println!("{}", rep);
        reps.push(rep);
    }
    // num = num.replace_multiple(&reps);
    // num.collect_symbol()
    let timn = Instant::now();

    println!("{}", num.log_display());
    let num = num.simplify_metrics();
    let elapsed = timn.elapsed();
    println!(
        "Metrics simplified: {},terms in  {:?}",
        num.log_print(Some(120)),
        elapsed
    );

    let timn = Instant::now();
    let simplified = num
        .simplify_gamma()
        .collect_rep(Minkowski {}.into())
        .to_dots();
    let elapsed = timn.elapsed();

    println!(
        "Gamma simplified: {}, {} terms in  {:?}",
        simplified.log_print(Some(120)),
        simplified.nterms(),
        elapsed
    );

    let expanded_simplified = simplified.expand();
    println!(
        "Gamma simplified: {}, {} terms",
        expanded_simplified.log_print(Some(120)),
        expanded_simplified.nterms(),
    );

    let evaluator_settings = EvaluatorSettings {
        iterative_orientation_optimization: false,
        horner_iterations: 10,
        spenso_execution_mode: (ExecutionMode::Sequential, ContractionMode::MinResultRank),
        store_atom: true,
        verbose: false,
        ..EvaluatorSettings::default()
    };
    let mut evaluator_param_builder = graph.param_builder.clone();

    // println!("Hello");
    for (edge_id, signature) in graph.loop_momentum_basis.edge_signatures.iter() {
        if signature.internal.iter().any(|sign| sign.is_sign()) {
            if graph.loop_momentum_basis.loop_edges.contains(&edge_id) {
                evaluator_param_builder
                    .pairs
                    .external_energies
                    .params
                    .push(function!(
                        GS.emr_mom,
                        Atom::num(edge_id.0 as i64),
                        Atom::from(ExpandedIndex::from_iter([0]))
                    ));

                evaluator_param_builder.pairs.update_ranges();
            } else {
                //has loop mom->is a non-tree edge -> energy is OSE-> no need for Q(0) rep
                let i = 0;

                evaluator_param_builder
                    .add_tagged_function::<Symbol>(
                        GS.emr_mom,
                        vec![
                            Atom::num(edge_id.0 as i64),
                            Atom::from(ExpandedIndex::from_iter([i])),
                        ],
                        format!("Q({edge_id}, {i})"),
                        vec![],
                        graph.loop_momentum_basis.loop_atom(
                            edge_id,
                            GS.emr_mom,
                            &[Atom::from(ExpandedIndex::from_iter([i]))],
                            true,
                        ) + graph.loop_momentum_basis.ext_atom(
                            edge_id,
                            GS.emr_mom,
                            &[Atom::from(ExpandedIndex::from_iter([i]))],
                            true,
                        ),
                    )
                    .unwrap();
            }
        }
    }

    // for r in &evaluator_param_builder.reps {
    //     println!("{}", r.replacement())
    // }

    // for r in &evaluator_param_builder.pairs {
    //     for param in &r.params {
    //         println!("{}", param)
    //     }
    // }
    // println!("Hi");
    let mut size_table_builder = Builder::new();
    size_table_builder.push_record([
        "expression".to_string(),
        "size (bytes)".to_string(),
        "number of multiplications".to_string(),
        "number of additions".to_string(),
        "spenso time".to_string(),
        "symbolica time".to_string(),
        "average time per evaluation".to_string(),
        "re".to_string(),
        "im".to_string(),
    ]);
    let params: Vec<Atom> = (&evaluator_param_builder.pairs)
        .into_iter()
        .flat_map(|p| p.params.clone())
        .collect();
    let mut rng = rand::rng();
    let values = params
        .iter()
        .map(|_| {
            Complex::new(
                F(rng.random_range(-1.0..1.0)),
                F(rng.random_range(-1.0..1.0)),
            )
        })
        .collect::<Vec<_>>();

    for (name, expr) in [
        ("Concretized", num),
        ("Simplified", simplified),
        ("Expanded", expanded_simplified),
    ] {
        println!("Creating evaluator for {}", name);
        let (mut evaluator, timings) = EvaluatorStack::new_with_timings(
            std::slice::from_ref(&expr),
            &evaluator_param_builder,
            &[],
            None,
            &evaluator_settings,
        )
        .unwrap();
        let operation_count = evaluator.single_parametric.f64_eager.count_operations();

        let mut out = vec![Complex::new(F(0.0), F(0.0))];

        println!("{}:", name);
        for inst in evaluator
            .single_parametric
            .f64_eager
            .export_instructions()
            .instructions
        {
            println!("{}", inst);
        }

        let time = std::time::Instant::now();
        for _ in 0..10000 {
            evaluator
                .single_parametric
                .f64_eager
                .evaluate(&values, &mut out);
        }
        let elapsed = time.elapsed();

        size_table_builder.push_record([
            name.to_string(),
            evaluator.single_parametric.exprs.as_ref().unwrap()[0]
                .as_view()
                .get_byte_size()
                .to_string(),
            operation_count.multiplications.to_string(),
            operation_count.additions.to_string(),
            format!("{:?}", timings.spenso_time),
            format!("{:?}", timings.symbolica_time),
            format!("{:?}", elapsed / 10000),
            out[0].re.to_string(),
            out[0].im.to_string(),
        ]);
    }
    let mut size_table = size_table_builder.build();
    size_table.with(Style::rounded());
    size_table.with(Modify::new(Columns::new(1..)).with(Alignment::right()));
    println!("{size_table}");
}
