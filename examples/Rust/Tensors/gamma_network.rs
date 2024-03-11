use std::ops::Neg;

use _gammaloop::tensor::{
    ufo::{euclidean_four_vector, gamma, mink_four_vector},
    AbstractIndex, ContractionCountStructure, FallibleMul, NumTensor, TensorNetwork,
};
use num::ToPrimitive;

use symbolica::domains::float::Complex;

fn gamma_net_num(
    minkindices: &[i32],
    vbar: [Complex<f64>; 4],
    u: [Complex<f64>; 4],
) -> TensorNetwork<NumTensor<ContractionCountStructure>> {
    let mut i: i32 = 0;
    let mut contracting_index = 0.into();
    let mut result: Vec<NumTensor<ContractionCountStructure>> =
        vec![euclidean_four_vector(contracting_index, &vbar).into()];
    for m in minkindices {
        let ui = contracting_index;
        contracting_index += 1.into();
        let uj = contracting_index;
        if *m > 0 {
            let p = [
                Complex::<f64>::new(1.0 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex::<f64>::new(1.1 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex::<f64>::new(1.2 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex::<f64>::new(1.3 + 0.01 * i.to_f64().unwrap(), 0.0),
            ];
            i += 1;
            result.push(mink_four_vector(usize::try_from(*m).unwrap().into(), &p).into());
            result.push(gamma(usize::try_from(*m).unwrap().into(), (ui, uj)).into());
        } else {
            result.push(
                gamma(
                    AbstractIndex::from(usize::try_from(m.neg()).unwrap() + 10000),
                    (ui, uj),
                )
                .into(),
            );
        }
    }
    result.push(euclidean_four_vector(contracting_index, &u).into());
    TensorNetwork::from(result)
}

fn main() {
    let one = Complex::<f64>::new(1.0, 0.0);

    let vbar = [
        one.mul_fallible(3.0).unwrap(),
        one.mul_fallible(3.1).unwrap(),
        one.mul_fallible(3.2).unwrap(),
        one.mul_fallible(3.3).unwrap(),
    ];
    let u = [
        one.mul_fallible(4.0).unwrap(),
        one.mul_fallible(4.1).unwrap(),
        one.mul_fallible(4.2).unwrap(),
        one.mul_fallible(4.3).unwrap(),
    ];

    let spacings: [i32; 2] = [20, 24];
    let mut start = 1;
    let mut ranges = Vec::new();

    for &spacing in spacings.iter() {
        ranges.push((start..start + spacing).chain(std::iter::once(-1)));
        start += spacing;
    }

    let vec: Vec<i32> = ranges.into_iter().flatten().collect();

    let mut net = gamma_net_num(&vec, vbar, u);
    net.contract_algo(|tn| tn.edge_to_min_degree_node_with_depth(1));

    println!("{}", net.dot());
    // assert_eq!(
    //     Complex {
    //         re: 5.341852612369398e16,
    //         im: -136854212797686.44
    //     },
    //     net.result().try_as_complex().unwrap().data()[0]
    // );
}
