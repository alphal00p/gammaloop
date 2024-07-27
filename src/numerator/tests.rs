use ahash::AHashMap;
use symbolica::atom::Atom;

use super::Numerator;

#[test]
fn color_3L_physical_photon() {
    let color_expr= Atom::parse(concat!(
        // "T(aind(coad(8,9),cof(3,8),coaf(3,7)))*T(aind(coad(8,14),cof(3,13),coaf(3,12)))*T(aind(coad(8,21),cof(3,20),coaf(3,19)))*T(aind(coad(8,26),cof(3,25),coaf(3,24)))*",
    // "id(aind(coaf(3,3),cof(3,4)))*id(aind(coaf(3,5),cof(3,6)))*id(aind(coaf(3,10),cof(3,11)))*id(aind(coaf(3,15),cof(3,16)))*id(aind(coaf(3,17),cof(3,18)))*id(aind(coaf(3,22),cof(3,23)))*id(aind(cof(3,4),coaf(3,24)))*id(aind(cof(3,6),coaf(3,3)))*id(aind(cof(3,8),coaf(3,5)))*id(aind(cof(3,11),coaf(3,7)))*id(aind(cof(3,13),coaf(3,10)))*id(aind(cof(3,16),coaf(3,12)))*id(aind(cof(3,18),coaf(3,15)))*id(aind(cof(3,20),coaf(3,17)))*id(aind(cof(3,23),coaf(3,19)))*id(aind(cof(3,25),coaf(3,22)))*id(aind(coad(8,21),coad(8,9)))*id(aind(coad(8,26),coad(8,14)))",
"id(aind(coaf(3,3),cof(3,4)))*id(aind(coaf(3,5),cof(3,6)))*id(aind(coaf(3,10),cof(3,11)))*id(aind(coaf(3,15),cof(3,16)))*id(aind(coaf(3,17),cof(3,18)))*id(aind(coaf(3,22),cof(3,23)))*id(aind(cof(3,4),coaf(3,24)))*id(aind(cof(3,6),coaf(3,3)))*id(aind(cof(3,8),coaf(3,5)))*id(aind(cof(3,11),coaf(3,7)))*id(aind(cof(3,13),coaf(3,10)))*id(aind(cof(3,16),coaf(3,12)))*id(aind(cof(3,18),coaf(3,15)))*id(aind(cof(3,20),coaf(3,17)))*id(aind(cof(3,23),coaf(3,19)))*id(aind(cof(3,25),coaf(3,22)))*id(aind(coad(8,21),coad(8,9)))*id(aind(coad(8,26),coad(8,14)))")).unwrap();

    let mut color_num = Numerator {
        expression: color_expr,
        network: None,
        const_map: AHashMap::new(),
    };

    color_num.process_color_simple();

    println!("{}", color_num.expression)
}
