use std::path::PathBuf;

fn main() {
    let manifest_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").expect("manifest dir"));
    let fjcore_dir = manifest_dir.join("resources").join("fjcore");
    let target = std::env::var("TARGET").expect("target triple");

    println!(
        "cargo:rerun-if-changed={}",
        fjcore_dir.join("fjcore.cc").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        fjcore_dir.join("fjcore.hh").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        fjcore_dir.join("wrapper.cc").display()
    );
    println!("cargo:rustc-link-lib={}", cpp_stdlib_for_target(&target));

    cc::Build::new()
        .cpp(true)
        .files([fjcore_dir.join("fjcore.cc"), fjcore_dir.join("wrapper.cc")])
        .include(&fjcore_dir)
        .flag_if_supported("-std=c++14")
        .warnings(false)
        .compile("fjcore_test");
}

fn cpp_stdlib_for_target(target: &str) -> &'static str {
    if target.contains("apple") || target.contains("freebsd") || target.contains("openbsd") {
        "c++"
    } else {
        "stdc++"
    }
}
