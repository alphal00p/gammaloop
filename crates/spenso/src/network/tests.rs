use crate::network::library::DummyKey;

use super::TensorNetworkError;

#[test]
fn display() {
    let a = TensorNetworkError::<i8, DummyKey>::Infallible;

    println!("{a}")
}

#[cfg(feature = "shadowing")]
#[test]
fn scalar_alias_refs_resolve_to_original_atom() {
    use symbolica::{atom::Atom, symbol};

    use super::{
        Network, NetworkLeaf, NetworkNode,
        store::{NetworkStore, TensorScalarStore},
        tags::scalar_store_alias,
    };

    let original = Atom::var(symbol!("x"));
    let mut net: Network<NetworkStore<(), Atom>, i8, i8> = Network::from_scalar(original.clone());
    let aliases = net.alias_scalar_refs(|_, _| true);

    assert_eq!(aliases.aliases_created(), 1);
    let (node, _, _) = net.graph.result().unwrap();
    let NetworkNode::Leaf(NetworkLeaf::Scalar(scalar)) = node else {
        panic!("expected scalar result node");
    };

    assert_eq!(net.store.get_scalar_ref(*scalar), &scalar_store_alias(0));
    assert_eq!(
        net.resolve_scalar_aliases(&aliases, scalar_store_alias(0)),
        original
    );
}
