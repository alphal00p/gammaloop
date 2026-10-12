use linnet::half_edge::generate::{
    Countable, GenerationError, Legs, Plan, Query, Structure, Stub, Theory,
};

#[test]
fn explanations_snapshot_the_lowered_request_and_infeasibility() {
    let photon = Stub::undirected('a');
    let electron = Stub::incoming('e');
    let qed = Theory::builder()
        .vertex([electron.clone(), electron.flipped(), photon.clone()], [])
        .build()
        .unwrap();
    let legs = Legs::new().leg(1u8, photon.clone()).leg(2, photon);
    let plan = Plan::new(
        &qed,
        &Query::new(&legs)
            .loops(2..=2)
            .structure(Structure::one_pi()),
    )
    .unwrap();
    assert_eq!(
        plan.explain().to_string(),
        include_str!("generation_reports/qed.txt")
    );

    let p = Stub::undirected('p');
    let weighted = Theory::builder()
        .weights(["g"])
        .vertex([p.clone(), p.clone(), p.clone()], [2])
        .vertex([p.clone(), p.clone(), p.clone()], [-1])
        .build()
        .unwrap();
    let plan = Plan::new(
        &weighted,
        &Query::<_, u8>::new(&Legs::new())
            .loops(3..=3)
            .weight("g", -2..=2)
            .structure(Structure {
                self_loops: true,
                zero_flow_edges: true,
                ..Structure::default()
            }),
    )
    .unwrap();
    assert_eq!(
        plan.explain().to_string(),
        include_str!("generation_reports/signed_vacuum.txt")
    );

    let scalar = Theory::builder()
        .vertex([p.clone(), p.clone(), p.clone()], [])
        .build()
        .unwrap();
    let identical = Legs::new().identical(0u8, p.clone(), 3);
    let plan = Plan::new(
        &scalar,
        &Query::new(&identical)
            .loops(0..=0)
            .count(Countable::Vertices(vec![0]), 1..=1),
    )
    .unwrap();
    assert_eq!(
        plan.explain().to_string(),
        include_str!("generation_reports/identical_legs.txt")
    );

    let two = Theory::builder()
        .vertex([p.clone(), p.clone()], [])
        .build()
        .unwrap();
    let legs = Legs::new().leg(1u8, p.clone()).leg(2, p.clone());
    let plan = Plan::new(
        &two,
        &Query::new(&legs)
            .structure(Structure {
                max_vertices: Some(3),
                ..Structure::default()
            })
            .prune_skeletons(|_| true)
            .prune_dressings(|_| true)
            .select(|_| true),
    )
    .unwrap();
    assert_eq!(
        plan.explain().to_string(),
        include_str!("generation_reports/finite_vertex_bound.txt")
    );

    let cubic = Theory::builder()
        .weights(["g"])
        .vertex([p.clone(), p.clone(), p], [1])
        .build()
        .unwrap();
    let report = match Plan::new(
        &cubic,
        &Query::new(&identical).loops(0..=0).weight("g", 2..=2),
    ) {
        Err(GenerationError::Infeasible(report)) => report,
        _ => panic!("one cubic vertex cannot carry coupling power two"),
    };
    assert_eq!(
        report.to_string(),
        include_str!("generation_reports/infeasible_order.txt")
    );
}
