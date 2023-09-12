fn main() {
    if cfg!(feature = "python_api") {
        pyo3_build_config::add_extension_module_link_args();
    }
    println!("cargo:rustc-link-search=./dependencies/fjcore");
    println!("cargo:rustc-link-lib=stdc++");
    #[cfg(target_os = "macos")]
    {
        println!("cargo:rustc-link-lib=gcc");
        println!("cargo:rustc-link-lib=gfortran");
    }
}
