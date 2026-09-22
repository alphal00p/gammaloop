fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("cargo:rerun-if-env-changed=EXTRA_MACOS_LIBS_FOR_GNU_GCC");

    #[cfg(target_os = "macos")]
    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        println!("cargo:rustc-link-lib=gcc_s");
        println!("cargo:rustc-link-lib=gcc");
        println!("cargo:rustc-link-lib=gfortran");
    }
    Ok(())
}
