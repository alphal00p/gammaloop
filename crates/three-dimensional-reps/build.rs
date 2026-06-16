fn main() {
    emit_macos_gcc_runtime_link();
}

fn emit_macos_gcc_runtime_link() {
    println!("cargo:rerun-if-env-changed=EXTRA_MACOS_LIBS_FOR_GNU_GCC");

    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_none() {
        return;
    }
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("macos") {
        return;
    }

    let lib_dir = "/opt/local/lib/libgcc";
    println!("cargo:rustc-link-search=native={lib_dir}");
    println!("cargo:rustc-link-lib=dylib=gcc_s.1.1");
    println!("cargo:rustc-link-arg=-lgcc_s.1.1");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{lib_dir}");
}
