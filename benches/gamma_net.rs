use std::ops::Neg;

use _gammaloop::tensor::{
    ufo::{
        euclidean_four_vector, euclidean_four_vector_sym, gamma, gammasym, mink_four_vector,
        mink_four_vector_sym,
    },
    AbstractIndex, FallibleMul, HistoryStructure, NumTensor, TensorNetwork,
};

use criterion::{criterion_group, criterion_main, Criterion};
use num::ToPrimitive;
use symbolica::domains::float::Complex;

use symbolica::representations::Symbol;
fn gamma_net_sym(
    minkindices: &[i32],
    vbar: [Complex<f64>; 4],
    u: [Complex<f64>; 4],
) -> TensorNetwork<NumTensor<HistoryStructure<Symbol>>> {
    let mut i = 0;
    let mut contracting_index = 0.into();
    let mut result: Vec<NumTensor<HistoryStructure<Symbol>>> =
        vec![euclidean_four_vector_sym(contracting_index, &vbar).into()];
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
            result.push(mink_four_vector_sym(usize::try_from(*m).unwrap().into(), &p).into());
            result.push(gammasym(usize::try_from(*m).unwrap().into(), (ui, uj)).into());
        } else {
            result.push(
                gammasym(
                    AbstractIndex::from(usize::try_from(m.neg()).unwrap() + 10000),
                    (ui, uj),
                )
                .into(),
            );
        }
    }
    result.push(euclidean_four_vector_sym(contracting_index, &u).into());
    TensorNetwork::from(result)
}

fn gamma_net(
    minkindices: &[i32],
    vbar: [Complex<f64>; 4],
    u: [Complex<f64>; 4],
) -> TensorNetwork<NumTensor> {
    let mut i = 0;
    let mut contracting_index = 0.into();
    let mut result: Vec<NumTensor> = vec![euclidean_four_vector(contracting_index, &vbar).into()];
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

fn criterion_benchmark(c: &mut Criterion) {
    let one = Complex::<f64>::new(1.0, 0.0);
    let _zero = Complex::<f64>::new(0.0, 0.0);

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

    let netsym = gamma_net_sym(&minkindices, vbar, u);
    let net = gamma_net(&minkindices, vbar, u);

    let mut group = c.benchmark_group("gamma_net");

    group.bench_function("gamma_net_contract_sym", |b| {
        b.iter_batched(
            || netsym.clone(),
            |mut netsym| {
                netsym.contract();
            },
            criterion::BatchSize::SmallInput,
        )
    });

    group.bench_function("gamma_net_contraction", |b| {
        b.iter_batched(
            || net.clone(),
            |mut net| {
                net.contract();
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
