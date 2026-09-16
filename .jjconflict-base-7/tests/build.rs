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
    cc::Build::new()
        .cpp(true)
        .files([fjcore_dir.join("fjcore.cc"), fjcore_dir.join("wrapper.cc")])
        .include(&fjcore_dir)
        // `cc-rs` treats any stderr output during `flag_if_supported` probes as a failure.
        // The Nix clang wrapper warns about Cargo's Darwin target spelling, so probing
        // drops the flag and clang falls back to C++17, where `std::auto_ptr` is removed.
        .std("c++14")
        .warnings(false)
        .compile("fjcore_test");

    // Keep the C++ standard library after the static archive on the linker line.
    println!("cargo:rustc-link-lib={}", cpp_stdlib_for_target(&target));
}

fn cpp_stdlib_for_target(target: &str) -> &'static str {
    if target.contains("apple") || target.contains("freebsd") || target.contains("openbsd") {
        "c++"
    } else {
        "stdc++"
    }
}
