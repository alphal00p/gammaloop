use std::path::Path;

fn main() {

    if cfg!(target_os = "macos") && std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        let lib_dir = "/opt/local/lib/libgcc";
        let lib_name = "gcc_s.1.1";
        let lib_path = format!("{}/lib{}.dylib", lib_dir, lib_name);
        if Path::new(&lib_path).exists() {
            println!("cargo:rustc-link-search=native={}", lib_dir);
            println!("cargo:rustc-link-lib={}", lib_name);
            println!("cargo:rustc-link-arg=-l{}", lib_name);
        }
    }

}
