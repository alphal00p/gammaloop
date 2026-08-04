use std::{cell::Cell, collections::HashMap, env, hint::black_box};

use criterion::{BatchSize, BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use spenso::{
    contraction::Contract,
    network::FastTensorSum,
    structure::{
        OrderedStructure,
        concrete_index::FlatIndex,
        representation::{Euclidean, RepName},
    },
    symbolic_parallelism::{
        SymbolicParallelism, set_symbolica_rayon_enabled, symbolica_rayon_enabled,
    },
    tensors::{
        data::{DataTensor, SparseTensor},
        parametric::ParamTensor,
    },
};
use symbolica::{atom::Atom, parser::ParseSettings, wrap_input};

type Structure = OrderedStructure<Euclidean>;
type Tensor = ParamTensor<Structure>;

const POLICIES: [(&str, SymbolicParallelism); 3] = [
    ("serial", SymbolicParallelism::Serial),
    ("auto", SymbolicParallelism::Auto),
    ("parallel", SymbolicParallelism::Parallel),
];
const FIXTURE_POOL_SIZE: usize = 8;

struct BuiltTensor {
    tensor: Tensor,
    atom_bytes: usize,
    max_atom_bytes: usize,
}

struct FastSumFixture {
    inputs: Vec<Tensor>,
    input_entries: usize,
    atom_bytes: usize,
    max_atom_bytes: usize,
}

#[derive(Clone, Copy)]
enum GroupedContractKind {
    SingleAxis,
    MultiAxis,
}

impl GroupedContractKind {
    fn name(self) -> &'static str {
        match self {
            Self::SingleAxis => "single_axis",
            Self::MultiAxis => "multi_axis",
        }
    }
}

struct ContractFixture {
    left: Tensor,
    right: Tensor,
    input_entries: usize,
    output_entries: usize,
    products: usize,
    atom_bytes: usize,
    max_atom_bytes: usize,
}

fn generated_atom(seed: &str, terms: usize) -> Atom {
    let expression = (0..terms)
        .map(|term| format!("x_{seed}_{term}*y_{seed}_{term}"))
        .collect::<Vec<_>>()
        .join("+");
    Atom::parse_with_default_namespace(wrap_input!(&expression), ParseSettings::default())
        .expect("generated benchmark Atom should parse")
}

fn sparse_tensor(
    structure: Structure,
    entry_count: usize,
    atom_terms: usize,
    tensor_name: &str,
) -> BuiltTensor {
    let mut atom_bytes = 0;
    let mut max_atom_bytes = 0;
    let elements = (0..entry_count)
        .map(|index| {
            let atom = generated_atom(&format!("{tensor_name}_{index}"), atom_terms);
            let bytes = atom.as_view().get_byte_size();
            atom_bytes += bytes;
            max_atom_bytes = max_atom_bytes.max(bytes);
            (FlatIndex::from(index), atom)
        })
        .collect::<HashMap<_, _>>();

    BuiltTensor {
        tensor: ParamTensor::composite(DataTensor::Sparse(SparseTensor {
            elements,
            zero: Atom::Zero,
            structure,
        })),
        atom_bytes,
        max_atom_bytes,
    }
}

fn fast_sum_fixture(
    entries: usize,
    input_count: usize,
    atom_terms: usize,
    fixture_index: usize,
) -> FastSumFixture {
    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(entries, 1)]).structure;
    let mut atom_bytes = 0;
    let mut max_atom_bytes = 0;
    let inputs = (0..input_count)
        .map(|tensor_index| {
            let built = sparse_tensor(
                structure.clone(),
                entries,
                atom_terms,
                &format!("fixture_{fixture_index:02}_sum_{tensor_index}"),
            );
            atom_bytes += built.atom_bytes;
            max_atom_bytes = max_atom_bytes.max(built.max_atom_bytes);
            built.tensor
        })
        .collect();

    FastSumFixture {
        inputs,
        input_entries: entries * input_count,
        atom_bytes,
        max_atom_bytes,
    }
}

fn contract_fixture(
    kind: GroupedContractKind,
    dimension: usize,
    atom_terms: usize,
    fixture_index: usize,
) -> ContractFixture {
    let (left_structure, right_structure, input_entries, output_entries, products) = match kind {
        GroupedContractKind::SingleAxis => (
            OrderedStructure::new(vec![
                Euclidean {}.new_slot(dimension, 1),
                Euclidean {}.new_slot(dimension, 2),
            ])
            .structure,
            OrderedStructure::new(vec![
                Euclidean {}.new_slot(dimension, 2),
                Euclidean {}.new_slot(dimension, 3),
            ])
            .structure,
            dimension.pow(2),
            dimension.pow(2),
            dimension.pow(3),
        ),
        GroupedContractKind::MultiAxis => (
            OrderedStructure::new(vec![
                Euclidean {}.new_slot(dimension, 1),
                Euclidean {}.new_slot(dimension, 2),
                Euclidean {}.new_slot(dimension, 3),
            ])
            .structure,
            OrderedStructure::new(vec![
                Euclidean {}.new_slot(dimension, 2),
                Euclidean {}.new_slot(dimension, 3),
                Euclidean {}.new_slot(dimension, 4),
            ])
            .structure,
            dimension.pow(3),
            dimension.pow(2),
            dimension.pow(4),
        ),
    };

    // Fully populated sparse storage reliably exercises Spenso's grouped-Atom
    // contraction paths while keeping the fixture shape easy to reproduce.
    let left = sparse_tensor(
        left_structure,
        input_entries,
        atom_terms,
        &format!("fixture_{fixture_index:02}_contract_left"),
    );
    let right = sparse_tensor(
        right_structure,
        input_entries,
        atom_terms,
        &format!("fixture_{fixture_index:02}_contract_right"),
    );

    ContractFixture {
        left: left.tensor,
        right: right.tensor,
        input_entries: input_entries * 2,
        output_entries,
        products,
        atom_bytes: left.atom_bytes + right.atom_bytes,
        max_atom_bytes: left.max_atom_bytes.max(right.max_atom_bytes),
    }
}

fn fast_sum(inputs: &[&Tensor]) -> Tensor {
    <Tensor as FastTensorSum>::fast_tensor_sum(black_box(inputs), None)
        .expect("sparse symbolic tensors should use FastTensorSum")
}

fn grouped_contract(fixture: &ContractFixture) -> Tensor {
    black_box(&fixture.left)
        .contract(black_box(&fixture.right))
        .expect("grouped sparse contraction should succeed")
}

fn assert_fast_sum_modes_match(inputs: &[&Tensor], case: &str) {
    set_symbolica_rayon_enabled(SymbolicParallelism::Serial);
    let serial = fast_sum(inputs);
    set_symbolica_rayon_enabled(SymbolicParallelism::Auto);
    let automatic = fast_sum(inputs);
    set_symbolica_rayon_enabled(SymbolicParallelism::Parallel);
    let parallel = fast_sum(inputs);
    assert_eq!(serial, automatic, "fast-sum Auto differs for {case}");
    assert_eq!(serial, parallel, "fast-sum modes differ for {case}");
}

fn assert_contract_modes_match(fixture: &ContractFixture, case: &str) {
    set_symbolica_rayon_enabled(SymbolicParallelism::Serial);
    let serial = grouped_contract(fixture);
    set_symbolica_rayon_enabled(SymbolicParallelism::Auto);
    let automatic = grouped_contract(fixture);
    set_symbolica_rayon_enabled(SymbolicParallelism::Parallel);
    let parallel = grouped_contract(fixture);
    assert_eq!(serial, automatic, "contraction Auto differs for {case}");
    assert_eq!(serial, parallel, "contraction modes differ for {case}");
}

fn benchmark_fast_sums(criterion: &mut Criterion, rayon_threads: usize) {
    let mut group = criterion.benchmark_group(format!(
        "symbolic_parallelism/fast_sparse_sum/rayon_threads_{rayon_threads}"
    ));

    for (entries, atom_terms) in [(4, 1), (64, 1), (256, 1), (4, 32), (64, 8), (256, 8)] {
        let input_count = 4;
        let fixtures = (0..FIXTURE_POOL_SIZE)
            .map(|fixture_index| fast_sum_fixture(entries, input_count, atom_terms, fixture_index))
            .collect::<Vec<_>>();
        let input_pool = fixtures
            .iter()
            .map(|fixture| fixture.inputs.iter().collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let fixture = &fixtures[0];
        let case = format!(
            "entries_{entries}_inputs_{input_count}_atom_terms_{atom_terms}_input_bytes_{}_max_atom_bytes_{}_fixture_pool_{FIXTURE_POOL_SIZE}",
            fixture.atom_bytes, fixture.max_atom_bytes,
        );

        for inputs in &input_pool {
            assert_fast_sum_modes_match(inputs, &case);
        }
        group.throughput(Throughput::Elements(fixture.input_entries as u64));
        for (mode, policy) in POLICIES {
            group.bench_with_input(
                BenchmarkId::new(mode, &case),
                &input_pool,
                |bench, input_pool| {
                    set_symbolica_rayon_enabled(policy);
                    let next_fixture = Cell::new(0);
                    bench.iter_batched(
                        || {
                            let index = next_fixture.get();
                            next_fixture.set((index + 1) % input_pool.len());
                            &input_pool[index]
                        },
                        |inputs| black_box(fast_sum(black_box(inputs))),
                        BatchSize::SmallInput,
                    );
                },
            );
        }
    }

    group.finish();
}

fn benchmark_grouped_contractions(criterion: &mut Criterion, rayon_threads: usize) {
    for (kind, cases) in [
        (
            GroupedContractKind::SingleAxis,
            &[(2, 1), (8, 1), (8, 8)][..],
        ),
        (GroupedContractKind::MultiAxis, &[(2, 1), (4, 1)][..]),
    ] {
        let mut group = criterion.benchmark_group(format!(
            "symbolic_parallelism/grouped_contract_{}/rayon_threads_{rayon_threads}",
            kind.name(),
        ));

        for &(dimension, atom_terms) in cases {
            let fixtures = (0..FIXTURE_POOL_SIZE)
                .map(|fixture_index| contract_fixture(kind, dimension, atom_terms, fixture_index))
                .collect::<Vec<_>>();
            let fixture = &fixtures[0];
            let case = format!(
                "dim_{dimension}_input_entries_{}_output_entries_{}_products_{}_atom_terms_{atom_terms}_input_bytes_{}_max_atom_bytes_{}_fixture_pool_{FIXTURE_POOL_SIZE}",
                fixture.input_entries,
                fixture.output_entries,
                fixture.products,
                fixture.atom_bytes,
                fixture.max_atom_bytes,
            );

            for fixture in &fixtures {
                assert_contract_modes_match(fixture, &case);
            }
            group.throughput(Throughput::Elements(fixture.products as u64));
            for (mode, policy) in POLICIES {
                group.bench_with_input(
                    BenchmarkId::new(mode, &case),
                    &fixtures,
                    |bench, fixtures| {
                        set_symbolica_rayon_enabled(policy);
                        let next_fixture = Cell::new(0);
                        bench.iter_batched(
                            || {
                                let index = next_fixture.get();
                                next_fixture.set((index + 1) % fixtures.len());
                                &fixtures[index]
                            },
                            |fixture| black_box(grouped_contract(black_box(fixture))),
                            BatchSize::SmallInput,
                        );
                    },
                );
            }
        }

        group.finish();
    }
}

fn symbolic_parallelism(criterion: &mut Criterion) {
    set_symbolica_rayon_enabled(SymbolicParallelism::Auto);
    assert!(
        symbolica_rayon_enabled(),
        "the symbolic_parallelism benchmark requires a licensed Symbolica because it deliberately forces Parallel mode"
    );

    // Broadcast both initializes the global pool before any measurement and
    // confirms the worker count encoded into every benchmark group name.
    let workers = rayon::broadcast(|context| context.index());
    let rayon_threads = rayon::current_num_threads();
    assert_eq!(workers.len(), rayon_threads);
    eprintln!(
        "symbolic_parallelism benchmark metadata: rayon_threads={rayon_threads} RAYON_NUM_THREADS={} available_parallelism={}",
        env::var("RAYON_NUM_THREADS").unwrap_or_else(|_| "default".to_owned()),
        std::thread::available_parallelism().map_or(1, usize::from),
    );

    benchmark_fast_sums(criterion, rayon_threads);
    benchmark_grouped_contractions(criterion, rayon_threads);
    set_symbolica_rayon_enabled(SymbolicParallelism::Auto);
}

criterion_group!(benches, symbolic_parallelism);
criterion_main!(benches);
