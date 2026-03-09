use super::*;
use crate::{
    data::SparseOrDense,
    shadowing::test::EXPLICIT_TENSOR_MAP,
    structure::{
        abstract_index::{ABSTRACTIND, AIND_SYMBOLS},
        representation::{Minkowski, RepName},
        NamedStructure, SmartShadowStructure,
    },
    symbolic::SymbolicTensor,
    symbolica_utils::SerializableSymbol,
    tensor_library::{ShadowedStructure, ETS},
    upgrading_arithmetic::FallibleSub,
};
use constcat::concat;
use symbolica::parse;

#[test]
fn other_network() {
    let mut net: Network<
        DataTensor<
            Complex<Rational>,
            SmartShadowStructure<SerializableSymbol, Vec<SerializableAtom>>,
        >,
        SerializableAtom,
    > = Network::new();

    net.contract().unwrap();
}

#[test]
fn parsing_id() {
    let expr = "γ(aind(mink(4,192),euc(4,105),euc(4,175)))";
    let atom = parse!(expr).unwrap();

    let sym_tensor: SymbolicTensor = atom.try_into().unwrap();

    let mut network = sym_tensor
        .to_network(&EXPLICIT_TENSOR_MAP.read().unwrap())
        .unwrap();

    // for (id, n) in &network.graph.nodes {
    //     println!("Node id: {:?}", id);
    //     println!("Structure {}", n.structure())
    // }
    // println!("{}", network.dot());
    network.contract().unwrap();
    // println!("Network res: {}", network.result().unwrap().0);
}

#[test]
fn parsing_single_contract() {
    let expr = "Q(15,mink(4,192))*γ(aind(mink(4,192),euc(4,105),euc(4,175)))";
    let atom = parse!(expr).unwrap();

    let sym_tensor: SymbolicTensor = atom.try_into().unwrap();

    let mut network = sym_tensor
        .to_network(&EXPLICIT_TENSOR_MAP.read().unwrap())
        .unwrap();

    println!("{}", network.dot());
    network.contract().unwrap();
    println!("Network res: {}", network.result().unwrap().0);
}

#[test]
fn parsing_addition_and_mul() {
    let expr = "(MT*id(aind(euc(4,105),euc(4,175)))+Q(15,aind(mink(4,192)))*γ(aind(mink(4,192),euc(4,105),euc(4,175))))";
    let atom = parse!(expr).unwrap();

    let sym_tensor: SymbolicTensor = atom.try_into().unwrap();

    let mut network = sym_tensor
        .to_network(&EXPLICIT_TENSOR_MAP.read().unwrap())
        .unwrap();

    for (id, n) in &network.graph.nodes {
        println!("Node id: {:?}", id);
        println!("Structure {}", n.structure())
    }

    network.contract().unwrap();
    println!("Network res: {}", network.result().unwrap().0);
}

#[test]
fn pslash_parse() {
    let _ = AIND_SYMBOLS.dind;

    let expr = "Q(15,dind(lor(4,75257)))   *γ(lor(4,75257),euc(4,1),euc(4,18))";
    let atom = parse!(expr).unwrap();

    let sym_tensor: SymbolicTensor = atom.try_into().unwrap();

    let network = sym_tensor
        .to_network(&EXPLICIT_TENSOR_MAP.read().unwrap())
        .unwrap();

    println!("{}", network.dot());
}

#[test]
fn three_loop_photon_parse() {
    // let _ = ETS.gamma;

    let expr = concat!(
        "-64/729*ee^6*G^4",
        "*(MT*id(euc(4,1),euc(4,18)))", //+Q(15,aind(loru(4,75257)))    *γ(aind(loru(4,75257),euc(4,1),euc(4,18))))",
        "*(MT*id(euc(4,3),euc(4,0)))", //+Q(6,aind(loru(4,17)))        *γ(aind(loru(4,17),euc(4,3),euc(4,0))))",
        "*(MT*id(euc(4,5),euc(4,2))   )", //+Q(7,aind(loru(4,35)))        *γ(aind(loru(4,35),euc(4,5),euc(4,2))))",
        "*(MT*id(euc(4,7),euc(4,4))   )", //+Q(8,aind(loru(4,89)))        *γ(aind(loru(4,89),euc(4,7),euc(4,4))))",
        "*(MT*id(euc(4,9),euc(4,6))   )", //+Q(9,aind(loru(4,233)))       *γ(aind(loru(4,233),euc(4,9),euc(4,6))))",
        "*(MT*id(euc(4,11),euc(4,8))  )", //+Q(10,aind(loru(4,611)))      *γ(aind(loru(4,611),euc(4,11),euc(4,8))))",
        "*(MT*id(euc(4,13),euc(4,10)) )", //+Q(11,aind(loru(4,1601)))    *γ(aind(loru(4,1601),euc(4,13),euc(4,10))))",
        "*(MT*id(euc(4,15),euc(4,12)) )", //+Q(12,aind(loru(4,4193)))    *γ(aind(loru(4,4193),euc(4,15),euc(4,12))))",
        "*(MT*id(euc(4,17),euc(4,14)) )", //+Q(13,aind(loru(4,10979)))   *γ(aind(loru(4,10979),euc(4,17),euc(4,14))))",
        "*(MT*id(euc(4,19),euc(4,16)) )", //+Q(14,aind(loru(4,28745)))   *γ(aind(loru(4,28745),euc(4,19),euc(4,16))))",
        "*Metric(mink(4,13),mink(4,8))",
        "*Metric(mink(4,15),mink(4,10))",
        // "*T(coad(8,9),cof(3,8),coaf(3,7))",
        // "*T(coad(8,14),cof(3,13),coaf(3,12))",
        // "*T(coad(8,21),cof(3,20),coaf(3,19))",
        // "*T(coad(8,26),cof(3,25),coaf(3,24))",
        // "*id(coaf(3,3),cof(3,4))*id(coaf(3,4),cof(3,24))*id(coaf(3,5),cof(3,6))*id(coaf(3,6),cof(3,3))",
        // "*id(coaf(3,8),cof(3,5))*id(coaf(3,10),cof(3,11))*id(coaf(3,11),cof(3,7))*id(coaf(3,13),cof(3,10))",
        // "*id(coaf(3,15),cof(3,16))*id(coaf(3,16),cof(3,12))*id(coaf(3,17),cof(3,18))*id(coaf(3,18),cof(3,15))",
        // "*id(coaf(3,20),cof(3,17))*id(coaf(3,22),cof(3,23))*id(coaf(3,23),cof(3,19))*id(coaf(3,25),cof(3,22))*id(coad(8,21),coad(8,9))*id(coad(8,26),coad(8,14))",
        "*γ(mink(4,6),euc(4,1),euc(4,0))",
        "*γ(mink(4,7),euc(4,3),euc(4,2))",
        "*γ(mink(4,8),euc(4,5),euc(4,4))",
        "*γ(mink(4,9),euc(4,7),euc(4,6))",
        "*γ(mink(4,10),euc(4,9),euc(4,8))",
        "*γ(mink(4,11),euc(4,11),euc(4,10))",
        "*γ(mink(4,12),euc(4,13),euc(4,12))",
        "*γ(mink(4,13),euc(4,15),euc(4,14))",
        "*γ(mink(4,14),euc(4,17),euc(4,16))",
        "*γ(mink(4,15),euc(4,19),euc(4,18))",
        "*ϵ(0,mink(4,6))",
        "*ϵ(1,mink(4,7))",
        "*ϵbar(2,mink(4,14))",
        "*ϵbar(3,mink(4,12))",
        "*ϵbar(4,mink(4,11))",
        "*ϵbar(5,mink(4,9))"
    );

    let atom = parse!(expr).unwrap();

    let sym_tensor: SymbolicTensor = atom.try_into().unwrap();

    let _network = sym_tensor
        .to_network(&EXPLICIT_TENSOR_MAP.read().unwrap())
        .unwrap();

    // println!("{}", network.rich_graph().dot());
}

// fn g(i: usize,) -> Atom {
//     let mink = Minkowski::rep(4);
//     let euc = Bispinor::rep(4);

//     function!(
//         ETS.gamma,
//         mink.new_slot(mu).to_atom(),
//         bis.new_slot(i).to_atom(),
//         bis.new_slot(j).to_atom()
//     )
// }
// fn g(mu: usize, nu: usize) -> Atom {
//     let mink = Minkowski::rep(4);

//     function!(
//         ETS.metric,
//         mink.new_slot(mu).to_atom(),
//         mink.new_slot(nu).to_atom()
//     )
// }

fn g_concrete(mu: usize, nu: usize) -> RealOrComplexTensor<f64, ShadowedStructure> {
    let _ = AIND_SYMBOLS.dind;
    let mink = LibraryRep::from(Minkowski {}).rep(4);

    NamedStructure::<_, (), LibraryRep>::from_iter([mink.slot(mu), mink.slot(nu)], ETS.metric, None)
        .to_shell()
        .to_explicit(&EXPLICIT_TENSOR_MAP.read().unwrap())
        .unwrap()
        .try_into_concrete()
        .unwrap()
}
#[test]
fn sparse_dense_addition() {
    let a = g_concrete(1, 3)
        .contract(&g_concrete(2, 4))
        .unwrap()
        .add_fallible(&g_concrete(1, 2).contract(&g_concrete(3, 4)).unwrap())
        .unwrap(); // * (g(1, 4) * g(2, 3));

    let b = g_concrete(1, 3)
        .to_dense()
        .contract(&g_concrete(2, 4).to_dense())
        .unwrap()
        .add_fallible(
            &g_concrete(1, 2)
                .to_dense()
                .contract(&g_concrete(3, 4).to_dense())
                .unwrap(),
        )
        .unwrap();

    println!("{}", a.sub_fallible(&b).unwrap().to_sparse());

    // println!("{}", b.to_sparse());
}

// #[test]
// fn expanded_vs_not() {
//     let a = (g(1, 3) * g(2, 4) + g(1, 2) * g(3, 4)); // * (g(1, 4) * g(2, 3));

//     let sym_tensor: SymbolicTensor = a.try_into().unwrap();
//     let mut network = sym_tensor.to_network::<PhysReps>().unwrap();
//     network.contract();
//     let res = network.result_tensor_smart().unwrap();

//     println!("{:?}", res.structure());
//     println!("{}", res);
// }
