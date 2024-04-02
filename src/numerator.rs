use symbolica::{id::Pattern, representations::Atom};

use crate::{graph::Graph, model::Model};

pub fn generate_numerator(graph: &Graph, model: &Model) -> Atom {
    let _ = model;
    let vatoms: Vec<Atom> = graph
        .vertices
        .iter()
        .flat_map(|v| v.contracted_vertex_rule(graph))
        .collect();
    // graph
    //     .edges
    //     .iter()
    //     .filter(|e| e.edge_type == EdgeType::Virtual)
    //     .map(|e| e.particle);
    let eatoms: Vec<Atom> = graph.edges.iter().map(|e| e.numerator(graph)).collect();
    let mut builder = Atom::new_num(1);

    for v in &vatoms {
        builder = builder * v;
    }

    for e in &eatoms {
        builder = builder * e;
    }

    let identity = Pattern::parse("Identity(i_,j_)").unwrap();
    builder = identity.replace_all(
        builder.as_view(),
        &Pattern::parse("id(euc(4,i_),euc(4,j_))").unwrap(),
        None,
        None,
    );

    let identity_lor = Pattern::parse("IdentityL(mu_,nu_)").unwrap();
    builder = identity_lor.replace_all(
        builder.as_view(),
        &Pattern::parse("id(lor(4,mu_),lor(4,nu_))").unwrap(),
        None,
        None,
    );

    let gamma = Pattern::parse("Gamma(mu_,i_,j_)").unwrap();
    builder = gamma.replace_all(
        builder.as_view(),
        &Pattern::parse("γ(lor(4,mu_),bis(4,i_),bis(4,j_))").unwrap(),
        None,
        None,
    );

    let gamma5 = Pattern::parse("Gamma5(i_,j_)").unwrap();
    builder = gamma5.replace_all(
        builder.as_view(),
        &Pattern::parse("γ5(bis(4,i_),bis(4,j_))").unwrap(),
        None,
        None,
    );

    let proj_m = Pattern::parse("ProjM(i_,j_)").unwrap();
    builder = proj_m.replace_all(
        builder.as_view(),
        &Pattern::parse("ProjM(bis(4,i_),bis(4,j_))").unwrap(),
        None,
        None,
    );

    let proj_p = Pattern::parse("ProjP(i_,j_)").unwrap();
    builder = proj_p.replace_all(
        builder.as_view(),
        &Pattern::parse("ProjP(bis(4,i_),bis(4,j_))").unwrap(),
        None,
        None,
    );

    let sigma = Pattern::parse("Sigma(mu_,nu_,i_,j_)").unwrap();
    builder = sigma.replace_all(
        builder.as_view(),
        &Pattern::parse("σ(lor(4,mu_),lor(4,nu_),bis(4,i_),bis(4,j_))").unwrap(),
        None,
        None,
    );

    let charge_conj = Pattern::parse("C(i_,j_)").unwrap();
    builder = charge_conj.replace_all(
        builder.as_view(),
        &Pattern::parse("C(bis(4,i_),bis(4,j_))").unwrap(),
        None,
        None,
    );

    let metric = Pattern::parse("Metric(mu_,nu_)").unwrap();
    builder = metric.replace_all(
        builder.as_view(),
        &Pattern::parse("Metric(lor(4,mu_),lor(4,nu_))").unwrap(),
        None,
        None,
    );

    builder
}

#[cfg(test)]
mod tests;
