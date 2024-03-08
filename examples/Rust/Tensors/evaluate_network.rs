use std::ops::Neg;

use _gammaloop::tensor::{
    ufo::{euclidean_four_vector, gamma, mink_four_vector},
    AbstractIndex, DenseTensor, FallibleMul, HasTensorData, MixedTensor, MixedTensors,
    NamedStructure, NumTensor, Representation, SetTensorData, Shadowable, Slot, SparseTensor,
    TensorNetwork, TensorStructure, VecStructure,
};
use ahash::{AHashMap, HashMap, HashMapExt};
use num::ToPrimitive;
use rand::{distributions::Uniform, Rng, SeedableRng};
use rand_xoshiro::Xoroshiro64Star;
use symbolica::{
    domains::float::Complex,
    representations::{Atom, AtomView},
    state::State,
};
use yaml_rust::yaml::Hash;

fn gamma_net_param(
    minkindices: &[i32],
    vbar: [Complex<f64>; 4],
    u: [Complex<f64>; 4],
) -> TensorNetwork<MixedTensor> {
    let mut i: i32 = 0;
    let mut contracting_index = 0.into();
    let mut result: Vec<MixedTensor> = vec![euclidean_four_vector(contracting_index, &vbar).into()];
    for m in minkindices {
        let ui = contracting_index;
        contracting_index += 1.into();
        let uj = contracting_index;
        if *m > 0 {
            let p: VecStructure = vec![Slot::from((
                usize::try_from(*m).unwrap().into(),
                Representation::Lorentz(4.into()),
            ))]
            .into();
            i += 1;
            let pid = State::get_global_state()
                .write()
                .unwrap()
                .get_or_insert_fn(&format!("p{}", i), None)
                .unwrap();

            result.push(p.shadow_with(pid).into());

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

fn test_tensor<D, S>(structure: S, seed: u64, range: Option<(D, D)>) -> SparseTensor<D, S>
where
    S: TensorStructure,
    D: rand::distributions::uniform::SampleUniform,
    Uniform<D>: Copy,

    rand::distributions::Standard: rand::distributions::Distribution<D>,
{
    let mut rng: Xoroshiro64Star = Xoroshiro64Star::seed_from_u64(seed);

    let mut tensor = SparseTensor::empty(structure);

    let density = tensor.size();

    if let Some((low, high)) = range {
        let multipliable = Uniform::new(low, high);
        for _ in 0..density {
            tensor
                .set_flat(rng.gen_range(0..tensor.size()), rng.sample(multipliable))
                .unwrap();
        }
    } else {
        for _ in 0..density {
            tensor
                .set_flat(rng.gen_range(0..tensor.size()), rng.gen())
                .unwrap();
        }
    }

    tensor
}

fn const_map_gen<'a, 'b>(params: &'a [MixedTensor], const_map: &mut HashMap<AtomView<'b>, f64>)
where
    'a: 'b,
{
    for (i, p) in params.iter().enumerate() {
        let pdata = test_tensor(p.structure().clone(), i as u64, Some((1.0, 10.0))).to_dense();
        p.try_as_symbolic()
            .unwrap()
            .try_as_dense()
            .unwrap()
            .append_const_map(&pdata, const_map);
    }
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
    let spacings: [i32; 2] = [2, 4];
    let mut start = 1;
    let mut ranges = Vec::new();

    for &spacing in spacings.iter() {
        ranges.push((start..start + spacing).chain(std::iter::once(-1)));
        start += spacing;
    }

    let vec: Vec<i32> = ranges.into_iter().flatten().collect();

    let mut net = gamma_net_param(&vec, vbar, u);
    net.generate_params();
    let mut const_map = AHashMap::new();
    let params = net.params.clone();
    const_map_gen(&params, &mut const_map);

    net.evaluate_float(&const_map);
    // println!("{:#?}", net.graph.nodes);
    net.contract();
    println!("{:?}", net.result().try_as_complex().unwrap().data()[0]);
}
