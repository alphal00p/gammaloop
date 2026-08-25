use feynkit::{GenerationOptions, Generator, Model, Process};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model_path = std::env::args()
        .nth(1)
        .ok_or_else(|| std::io::Error::other("usage: generate <model.json>"))?;
    let model = Model::from_path(model_path)?;
    let process = Process::amplitude(["e-", "e+"], ["mu-", "mu+"]);
    let result = Generator::new(model).generate(&process, &GenerationOptions::default())?;

    for diagram in result.diagrams {
        println!("{}", diagram.to_dot()?);
    }
    Ok(())
}
