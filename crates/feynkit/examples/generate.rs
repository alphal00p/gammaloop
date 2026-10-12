use feynkit::{
    Model, ModelProcessExt,
    generator::{Loops, Naming, NumeratorOptions, Run, TopologyRules},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model_path = std::env::args()
        .nth(1)
        .ok_or_else(|| std::io::Error::other("usage: generate <model.json>"))?;
    let model = Model::from_path(model_path)?;
    let process = model.process(["e-", "e+"], ["mu-", "mu+"])?;
    let run = Run::default();
    let result = process
        .amplitude(Loops::total(0..=0))
        .with_topology(TopologyRules::new().with_self_loops(false))
        .generate(model.into(), &run)?
        .build(
            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
            &run,
        )?;

    for diagram in result.iter() {
        println!("{}", diagram.to_dot()?);
    }
    Ok(())
}
