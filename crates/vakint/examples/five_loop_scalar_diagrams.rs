//! Generate reusable vacuum DOT fixtures through the public model pipeline.

use feynkit_generator::{GenerationOptions, Process};
use feynkit_model::Model;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let directory = std::env::args_os()
        .nth(1)
        .ok_or("usage: five_loop_scalar_diagrams OUTPUT_DIRECTORY")?;
    let directory = std::path::Path::new(&directory);
    if directory.exists() {
        return Err("output directory must not already exist".into());
    }
    let options = GenerationOptions::default()
        .with_loop_count(5, 5)?
        .threads(1)
        .allow_self_loops(true)
        .max_vertices(8);
    let generated = Process::new(Vec::<&str>::new(), Vec::<&str>::new())
        .generate_diagrams(Model::phi_3_4(), &options)?;
    std::fs::create_dir_all(directory)?;
    for (index, diagram) in generated.diagrams.iter().enumerate() {
        std::fs::write(
            directory.join(format!("vacuum_{index:04}.dot")),
            diagram.to_dot()?,
        )?;
    }
    println!(
        "{} generated five-loop scalar vacuum graphs; graph counts are not integral-class counts",
        generated.diagrams.len()
    );
    Ok(())
}
