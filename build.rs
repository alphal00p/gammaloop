use vergen_gitcl::{Emitter, GitclBuilder};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(not(test))]
    {
        let git = GitclBuilder::default().branch(true).build()?;
        Emitter::default().add_instructions(&git)?.emit()?;
    }

    #[cfg(target_os = "macos")]
    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        println!("cargo:rustc-link-lib=gcc_s");
        println!("cargo:rustc-link-lib=gcc");
        println!("cargo:rustc-link-lib=gfortran");
    }
    Ok(())
}
