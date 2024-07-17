#!/usr/bin/env python3
import subprocess


def revert_revision():
    # pass
    # Pin the revision of the branch currently working for gammaLoop
    print(subprocess.check_output(
        ['git', 'reset', '--hard', '13d6bfffd5c6b40064739bdf4bdcfc96beb4f324']).decode('ascii').strip())


def patch_lib_rs():
    with open('./src/lib.rs', 'r', encoding='utf8') as f_in:
        symbolica_lib_rs = f_in.read()
    if 'GAMMALOOP_USER' not in symbolica_lib_rs:
        with open('./src/lib.rs', 'w', encoding='utf8') as f_out:
            f_out.write(symbolica_lib_rs.replace(
                '.or(env::var("SYMBOLICA_LICENSE").ok());',
                """.or(env::var("SYMBOLICA_LICENSE").ok());
            
        /* START OF GAMMALOOP CUSTOM LICENSE MODIFICATION */
        /* DISCLAIMER
        |
        | Allow the user to set the GAMMALOOP_USER environment variable to bypass the license
        | This is a special measure unique to the gammaLoop fork of Symbolica and is not allowed
        | to be used outside of the context of gammaLoop.
        |
        */
        if let Some(gammaloop_key) = key.clone() {
            if gammaloop_key == "GAMMALOOP_USER" {
                LICENSED.store(true, Relaxed);                
                return Ok(());
            }
        }
        /* END OF GAMMALOOP CUSTOM LICENSE  MODIFICATION */"""))


def patch_build_rs():
    with open('./build.rs', 'r', encoding='utf8') as f_in:
        symbolica_build_rs = f_in.read()
    if 'pyo3_build_config::add_extension_module_link_args();' not in symbolica_build_rs:
        with open('./build.rs', 'w', encoding='utf8') as f_out:
            f_out.write(symbolica_build_rs.replace(
                'fn main() {',
                """fn main() {
    pyo3_build_config::add_extension_module_link_args();"""))


def patch_cargo_toml():
    with open('./Cargo.toml', 'r', encoding='utf8') as f_in:
        cargo_toml = f_in.read()
    if 'pyo3-build-config = "*"' not in cargo_toml:
        with open('./Cargo.toml', 'w', encoding='utf8') as f_out:
            f_out.write('\n'.join([
                cargo_toml.replace(
                    'crate-type = ["lib"]', 'crate-type = ["lib", "cdylib"]'),
                '',
                '[build-dependencies]',
                'pyo3-build-config = "*"'
            ]))


if __name__ == '__main__':
    revert_revision()
    patch_lib_rs()
    patch_build_rs()
    patch_cargo_toml()
    print("Symbolica successfully patched for use in GammaLoop")
