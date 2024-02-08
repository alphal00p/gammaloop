use std::ops::Neg;

use _gammaloop::tensor::{
    ufo_spin_tensors::{
        euclidean_four_vector, euclidean_four_vector_sym, gamma, gammasym, mink_four_vector,
        mink_four_vector_sym,
    },
    NumTensors, TensorNetwork,
};

use criterion::{criterion_group, criterion_main, Criterion};
use num::{complex::Complex64, ToPrimitive};
use smartstring::alias::String;
use symbolica::{representations::Identifier, state::State};
fn gamma_net_sym(
    minkindices: &[i32],
    vbar: [Complex64; 4],
    u: [Complex64; 4],
    state: &mut State,
) -> TensorNetwork<NumTensors<Identifier>> {
    let mut i = 0;
    let mut contracting_index = 0;
    let mut result: Vec<NumTensors<Identifier>> =
        vec![euclidean_four_vector_sym(contracting_index, &vbar, state).into()];
    for m in minkindices {
        let ui = contracting_index;
        contracting_index += 1;
        let uj = contracting_index;
        if *m > 0 {
            let p = [
                Complex64::new(1.0 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.1 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.2 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.3 + 0.01 * i.to_f64().unwrap(), 0.0),
            ];
            i += 1;
            result.push(mink_four_vector_sym(usize::try_from(*m).unwrap(), &p, state).into());
            result.push(gammasym(usize::try_from(*m).unwrap(), (ui, uj), state).into());
        } else {
            result
                .push(gammasym(usize::try_from(m.neg()).unwrap() + 10000, (ui, uj), state).into());
        }
    }
    result.push(euclidean_four_vector_sym(contracting_index, &u, state).into());
    TensorNetwork::new(result)
}

fn gamma_net(
    minkindices: &[i32],
    vbar: [Complex64; 4],
    u: [Complex64; 4],
) -> TensorNetwork<NumTensors<String>> {
    let mut i = 0;
    let mut contracting_index = 0;
    let mut result: Vec<NumTensors<String>> =
        vec![euclidean_four_vector(contracting_index, &vbar).into()];
    for m in minkindices {
        let ui = contracting_index;
        contracting_index += 1;
        let uj = contracting_index;
        if *m > 0 {
            let p = [
                Complex64::new(1.0 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.1 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.2 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.3 + 0.01 * i.to_f64().unwrap(), 0.0),
            ];
            i += 1;
            result.push(mink_four_vector(usize::try_from(*m).unwrap(), &p).into());
            result.push(gamma(usize::try_from(*m).unwrap(), (ui, uj)).into());
        } else {
            result.push(gamma(usize::try_from(m.neg()).unwrap() + 10000, (ui, uj)).into());
        }
    }
    result.push(euclidean_four_vector(contracting_index, &u).into());
    TensorNetwork::new(result)
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
    let mut state = State::new();
    let one = Complex64::new(1.0, 0.0);
    let _zero = Complex64::new(0.0, 0.0);

    let vbar = [one * 3.0, one * 3.1, one * 3.2, one * 3.3];
    let u = [one * 4.0, one * 4.1, one * 4.2, one * 4.3];

    let minkindices = indices(20, 24);

    let netsym = gamma_net_sym(&minkindices, vbar, u, &mut state);
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
