use std::{fmt::Debug, ops::Neg};

use _gammaloop::tensor::{
    ufo::{euclidean_four_vector, gamma},
    AbstractIndex, ContractionCountStructure, FallibleMul, MixedTensor, Representation,
    SetTensorData, Slot, SparseTensor, TensorNetwork, TensorStructure,
};

use ahash::{AHashMap, HashMap};
use criterion::{criterion_group, criterion_main, Criterion};

use rand::{distributions::Uniform, Rng, SeedableRng};
use rand_xoshiro::Xoroshiro64Star;
use symbolica::{
    domains::float::Complex,
    representations::{Atom, AtomView},
    state::State,
};

fn indices(n: i32, m: i32) -> Vec<i32> {
    let spacings: [i32; 2] = [n, m];
    let mut start = 1;
    let mut ranges = Vec::new();

    for &spacing in spacings.iter() {
        ranges.push((start..start + spacing).chain(std::iter::once(-1)));
        start += spacing;
    }

    ranges.into_iter().flatten().collect()
}

fn gamma_net_param(
    minkindices: &[i32],
    vbar: [Complex<f64>; 4],
    u: [Complex<f64>; 4],
) -> TensorNetwork<MixedTensor<ContractionCountStructure>> {
    let mut i: i32 = 0;
    let mut contracting_index = 0.into();
    let mut result: Vec<MixedTensor<ContractionCountStructure>> =
        vec![euclidean_four_vector(contracting_index, &vbar).into()];
    for m in minkindices {
        let ui = contracting_index;
        contracting_index += 1.into();
        let uj = contracting_index;
        if *m > 0 {
            let p: ContractionCountStructure = vec![Slot::from((
                usize::try_from(*m).unwrap().into(),
                Representation::Lorentz(4.into()),
            ))]
            .into_iter()
            .collect();
            i += 1;
            let pid = State::get_symbol(&format!("p{}", i));

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

fn test_tensor<S>(structure: S) -> SparseTensor<symbolica::domains::float::Complex<f64>, S>
where
    S: TensorStructure,
{
    let mut rng: Xoroshiro64Star = Xoroshiro64Star::from_entropy();

    let mut tensor = SparseTensor::empty(structure);

    let density = tensor.size();

    let multipliable = Uniform::new(1., 10.);

    for _ in 0..density {
        tensor
            .set_flat(
                rng.gen_range(0..tensor.size()),
                Complex::<f64>::new(rng.sample(multipliable), rng.sample(multipliable)),
            )
            .unwrap();
    }

    tensor
}

fn const_map_gen<'a, 'b, I>(
    params: &'a [MixedTensor<I>],
    const_map: &mut HashMap<AtomView<'b>, symbolica::domains::float::Complex<f64>>,
) where
    'a: 'b,
    I: TensorStructure + Clone + Debug,
{
    #[allow(clippy::unused_enumerate_index)]
    for (_i, p) in params.iter().enumerate() {
        let pdata = test_tensor(p.structure().clone()).to_dense();
        p.try_as_symbolic()
            .unwrap()
            .try_as_dense()
            .unwrap()
            .append_const_map(&pdata, const_map);
    }
}
fn criterion_benchmark(c: &mut Criterion) {
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
    let minkindices = indices(20, 24);

    let mut net = gamma_net_param(&minkindices, vbar, u);
    net.generate_params();
    let params = net.params.clone();
    println!("{:?}", params.len());
    net.contract_algo(|tn| tn.edge_to_min_degree_node_with_depth(2));
    let mut const_map = AHashMap::new();

    let i = Atom::new_var(State::I);
    const_map.insert(i.as_view(), Complex::<f64>::new(0., 1.));

    let mut group = c.benchmark_group("evaluate_net");

    group.bench_function("Evaluate_net", |b| {
        b.iter_batched(
            || net.clone(),
            |mut net| {
                const_map_gen(&params, &mut const_map);

                net.evaluate_complex(&const_map);

                net.contract();
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
