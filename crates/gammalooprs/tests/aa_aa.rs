use std::ops::Deref;
use std::time::Instant;

use gammalooprs::graph::parse::IntoGraph;
use gammalooprs::integrands::process::{
    evaluators::EvaluatorStack, param_builder::ParamBuilderGraph,
};
use gammalooprs::numerator::symbolica_ext::AtomCoreExt;
use gammalooprs::processes::EvaluatorSettings;
use gammalooprs::utils::symbolica_ext::{CallSymbol, LogPrint};
use gammalooprs::utils::{FUN_LIB, GS, TENSORLIB};
use gammalooprs::{dot, graph::Graph, uv::UltravioletGraph};
use idenso::dirac::GammaSimplifier;
use idenso::shorthands::metric::MetricSimplifier;
use idenso::shorthands::schoonschip::Schoonschip;
use linnet::half_edge::subgraph::SubSetOps;
use spenso::shadowing::TensorCollectExt;
use spenso::structure::representation::Minkowski;
use spenso::{
    network::{MinResultRank, Sequential},
    structure::concrete_index::ExpandedIndex,
};
use symbolica::atom::{Atom, AtomCore, FunctionBuilder, Symbol};
use symbolica::id::Replacement;
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

    let mut num = graph
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

    let evaluator_settings = EvaluatorSettings {
        iterative_orientation_optimization: false,
        verbose: true,
        ..EvaluatorSettings::default()
    };
    let mut evaluator_param_builder = graph.param_builder.clone();

    println!("Hello");
    for (edge_id, signature) in graph.loop_momentum_basis.edge_signatures.iter() {
        // if !lmb.loop_edges.contains(&edge_id) {
        if signature.internal.iter().any(|sign| sign.is_sign()) {
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

    println!("Hi");
    let mut size_table_builder = Builder::new();
    size_table_builder.push_record([
        "expression".to_string(),
        "size (bytes)".to_string(),
        "number of multiplications".to_string(),
        "number of additions".to_string(),
        "spenso time".to_string(),
        "symbolica time".to_string(),
    ]);

    for (name, expr) in [("Concretized", num), ("Simplified", simplified)] {
        println!("Creating evaluator for {}", name);
        let (evaluator, timings) = EvaluatorStack::new_with_timings(
            std::slice::from_ref(&expr),
            &evaluator_param_builder,
            &[],
            None,
            &evaluator_settings,
        )
        .unwrap();
        let (additions, multiplications) = evaluator.single_parametric.f64_eager.count_operations();
        println!("{}", name);
        for inst in evaluator
            .single_parametric
            .f64_eager
            .export_instructions()
            .0
        {
            println!("{}", inst);
        }
        size_table_builder.push_record([
            name.to_string(),
            expr.as_view().get_byte_size().to_string(),
            multiplications.to_string(),
            additions.to_string(),
            format!("{:?}", timings.spenso_time),
            format!("{:?}", timings.symbolica_time),
        ]);
    }
    let mut size_table = size_table_builder.build();
    size_table.with(Style::rounded());
    size_table.with(Modify::new(Columns::new(1..)).with(Alignment::right()));
    println!("{size_table}");
}
