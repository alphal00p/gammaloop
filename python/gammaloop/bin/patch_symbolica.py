#!/usr/bin/env python3
import subprocess
import os
import re

GL_PATH = os.path.abspath(os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir, os.path.pardir))

DEFAULT_SYMBOLICA_REVISION = 'latest'
# DEFAULT_SYMBOLICA_REVISION = 'eeed330c824e6ad8f4bc78472864e6d8121b1a12'


def revert_revision():
    target_symbolica_revision = os.environ.get(
        'SYMBOLICA_REVISION_HASH', DEFAULT_SYMBOLICA_REVISION)
    # Pin the revision of the branch currently working for gammaLoop
    if target_symbolica_revision != 'latest':
        print(subprocess.check_output(
            ['git', 'reset', '--hard', target_symbolica_revision]).decode('ascii').strip())


def patch_lib_rs():
    with open('./src/lib.rs', 'r', encoding='utf8') as f_in:
        symbolica_lib_rs = f_in.read()
    if 'GAMMALOOP_USER' not in symbolica_lib_rs:
        symbolica_lib_rs = "#![allow(warnings)]\n" + symbolica_lib_rs

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
    if cfg!(target_os = "macos") {
        pyo3_build_config::add_extension_module_link_args();
        if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
            println!("cargo:rustc-link-lib=gcc_s");
        }
    }"""))


# Something like this can be useful on some user's systems with this
# error when loading the pyo3 symbolica library
#   Symbol not found: ___emutls_get_address
# This will however not be called by default
def patch_build_linking():
    with open('./build.rs', 'r', encoding='utf8') as f_in:
        symbolica_build_rs = f_in.read()
    if 'gcc_s' not in symbolica_build_rs:
        with open('./build.rs', 'w', encoding='utf8') as f_out:
            f_out.write(symbolica_build_rs.replace(
                'fn main() {',
                """fn main() {
    println!("cargo:rustc-link-lib=gcc_s");
    println!("cargo:rustc-link-search=/opt/local/lib/libgcc");"""))


def get_symbolica_version_in_gammaloop_cargo_toml():
    try:
        if not os.path.isfile(os.path.join(GL_PATH, 'Cargo.toml')):
            return None
        with open(os.path.join(GL_PATH, 'Cargo.toml'), 'r') as f_in:
            symbolica_version = re.findall(
                r'symbolica\s*=\s*\{\s*version\s*=\s*\"(?P<version>.*)\"\s*\}', f_in.read())[0]
        if symbolica_version == '*':
            raise BaseException(
                "Symbolica version is specified in Cargo.toml must be pinned to a specific version")
        return symbolica_version
    except Exception as e:
        raise BaseException(
            "Could not identify symbolica version specified in Cargo.toml") from e


def patch_cargo_toml():
    with open('./Cargo.toml', 'r', encoding='utf8') as f_in:
        cargo_toml = f_in.read()

    requested_symbolica_version = get_symbolica_version_in_gammaloop_cargo_toml()
    current_version_number = None
    modified_cargo_toml = []
    in_package_group = False
    symbolica_version_modified = None
    for line in cargo_toml.split('\n'):
        if not in_package_group:
            if '[package]' in line:
                in_package_group = True
            modified_cargo_toml.append(line)
        else:
            if line.strip().startswith('version'):
                try:
                    current_version_number = re.findall(
                        r'version\s*=\s*\"(?P<version>.*)\"', line)[0]
                    if requested_symbolica_version is not None and current_version_number != requested_symbolica_version:
                        symbolica_version_modified = True
                        modified_cargo_toml.append(
                            f'version = "{requested_symbolica_version}"')
                    else:
                        symbolica_version_modified = False
                        modified_cargo_toml.append(line)
                except Exception as e:
                    raise BaseException(
                        "Could not identify symbolica version specified in symbolica's Cargo.toml") from e
            else:
                if line.strip().startswith('['):
                    in_package_group = False
                modified_cargo_toml.append(line)

    if current_version_number is None:
        raise BaseException(
            "Could not find version field in symbolica's Cargo.toml")
    elif symbolica_version_modified:
        print(f"Patched symbolica's Cargo.toml to specify version {requested_symbolica_version} requested by gammaLoop (it was {current_version_number})")  # nopep8
        cargo_toml = '\n'.join(modified_cargo_toml)
        with open('./Cargo.toml', 'w', encoding='utf8') as f_out:
            f_out.write(cargo_toml)

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
    # patch_build_linking()
    patch_cargo_toml()
    print("Symbolica successfully patched for use in GammaLoop")
