fn main() {
    // if cfg!(feature = "fjcore") {
    //     println!("cargo:rustc-link-search=./python/gammaloop/dependencies/fjcore");
    //     println!("cargo:rustc-link-lib=stdc++");
    // }
    #[cfg(target_os = "macos")]
    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        println!("cargo:rustc-link-lib=gcc_s");
        println!("cargo:rustc-link-lib=gcc");
        println!("cargo:rustc-link-lib=gfortran");
    }
}
