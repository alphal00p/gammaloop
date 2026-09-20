#import "../../shared.typ": source-link

#let quickstart-rust = [
= Using FeynKit from Rust

Use the source checkout corresponding to this manual. From its root, the maintained facade
example accepts a normalized model file:

// docs-example: syntax
```sh
nix develop
cargo run -p feynkit --example generate -- crates/feynkit-model/tests/fixtures/sm.json
```

That example generates electron–positron to muon–antimuon tree diagrams. For a smaller first
program, use the scalar fixture and the same owning APIs:

// docs-example: compile feynkit-rust-quickstart
```rust
use feynkit::{GenerationOptions, Generator, Model, Process};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model = Model::from_path(
        "crates/feynkit-model/tests/fixtures/scalars_2p_3p.json",
    )?;
    let generator = Generator::new(model);
    let process = Process::amplitude(["scalar_0"], ["scalar_0", "scalar_0"]);
    let options = GenerationOptions::default().max_vertices(3);
    let result = generator.generate(&process, &options)?;

    assert!(!result.diagrams.is_empty());
    for diagram in result.diagrams {
        diagram.validate()?;
        assert_eq!(diagram.loop_count(), 0);
        println!("{}", diagram.to_dot()?);
    }
    Ok(())
}
```

The program's model path is relative to the repository root. Success means a nonempty set of
validated tree diagrams with deterministic names and DOT output. Symbolica operations require
the license configuration appropriate to your environment.

For an application outside the workspace, depend on `feynkit` from the same checkout or pin the
repository revision you tested. The default facade features include model, graph, generator,
CFF, kinematics, and tensor reduction; enable `ufo` only when the application also supplies an
attached Python interpreter. Focused crates are available when fewer dependencies are needed.
See #source-link("crates/feynkit/Cargo.toml", label: "the feature manifest") and the
#link("reference/interfaces/")[Rust component map].
]
