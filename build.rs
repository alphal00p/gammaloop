use vergen::Emitter;
use vergen_gitcl::GitclBuilder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // if cfg!(feature = "python_api") {
    //    pyo3_build_config::add_extension_module_link_args();
    // }
    // if cfg!(feature = "fjcore") {
    //     println!("cargo:rustc-link-search=./python/gammaloop/dependencies/fjcore");
    //     println!("cargo:rustc-link-lib=stdc++");
    // }

    #[cfg(not(test))]
    {
        let git = GitclBuilder::default().branch(true).build()?;
        Emitter::default().add_instructions(&git)?.emit()?;
    }

    // Help Cargo know when to rerun
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/refs/heads");

    #[cfg(target_os = "macos")]
    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        println!("cargo:rustc-link-lib=gcc_s");
        println!("cargo:rustc-link-lib=gcc");
        println!("cargo:rustc-link-lib=gfortran");
    }
    Ok(())
}
