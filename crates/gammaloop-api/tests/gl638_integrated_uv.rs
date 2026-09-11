//! Standalone reproducer for the GL638 integrated-UV empty-residue failure.
#![cfg(feature = "cli")]

use std::{fs, process::Command};

#[test]
#[ignore = "known bug: GL638 integrated UV generation returns an empty direct-3D residue family"]
fn gl638_single_orientation_integrated_uv_generation() -> std::io::Result<()> {
    let scratch = tempfile::tempdir()?;
    // Preserve the physical SM numerator, original edge IDs/LMB and process cuts.
    // Threshold metadata, integration, state export and post-generation model
    // changes are unnecessary: the error occurs while generating the UV terms.
    fs::write(
        scratch.path().join("GL638.dot"),
        r#"digraph GL638{
    num = "1𝑖";
    overall_factor = "NumeratorDependentGrouping(638,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(639,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))";
    overall_factor_evaluated = "-4";
    projector = "u(1,spenso::bis(4,hedge(1)))*ubar(1,spenso::bis(4,hedge(29)))*v(0,spenso::bis(4,hedge(28)))*vbar(0,spenso::bis(4,hedge(0)))";

    0 [int_id="V_141"];
    1 [int_id="V_141"];
    2 [int_id="V_134"];
    3 [int_id="V_134"];
    4 [int_id="V_137"];
    5 [int_id="V_137"];
    6 [int_id="V_137"];
    7 [int_id="V_137"];
    8 [int_id="V_98"];
    9 [int_id="V_98"];
    exte0     [style=invis];
    exte0    -> 9:0     [id=0 is_cut="28" particle="e+"];
    exte1     [style=invis];
    exte1    -> 9:1     [id=1 is_cut="29" particle="e-"];
    0:2    -> 1:3     [id=2 particle="H"];
    0:4    -> 6:5     [id=3 lmb_id="0" particle="t"];
    7:6    -> 0:7     [id=4 lmb_id="1" particle="t"];
    3:8    -> 1:9     [id=5 particle="t"];
    1:10    -> 7:11     [id=6 particle="t"];
    2:12    -> 5:13     [id=7 lmb_id="2" particle="t"];
    6:14    -> 2:15     [id=8 particle="t"];
    2:16    -> 9:17     [id=9 particle="a"];
    4:18    -> 3:19     [id=10 lmb_id="3" particle="t"];
    3:20    -> 8:21     [id=11 particle="a"];
    5:22    -> 4:23     [id=12 particle="t"];
    4:24    -> 7:25     [id=13 particle="g"];
    5:26    -> 6:27     [id=14 particle="g"];
    exte15     [style=invis];
    8:28    -> exte15     [id=15 is_cut="28" particle="e+"];
    exte16     [style=invis];
    8:29    -> exte16     [id=16 is_cut="29" particle="e-"];
}"#,
    )?;
    fs::write(
        scratch.path().join("run.toml"),
        r#"
[cli_settings.state]
folder = "state"

[cli_settings.global.n_cores]
generate = 1

[cli_settings.global.generation.orientation_pattern]
pat = "(+,+,+,+,-,-,-,+,-,0,+,0,+,-,+)"

[cli_settings.global.generation.evaluator]
do_algebra = false
iterative_orientation_optimization = false
store_atom = true
compile = false
summed = false
summed_function_map = false

[cli_settings.global.generation.uv]
softct = false
inner_products = true
generate_integrated = true
subtract_uv = true

[cli_settings.global.generation.threshold_subtraction]
enable_thresholds = true
check_esurface_at_generation = true
esurface_existence_threshold = 1.0e-7
skip_thresholds_that_are_cuts = true

[default_runtime_settings.kinematics]
e_cm = 1000.0

[default_runtime_settings.kinematics.externals]
type = "constant"

[default_runtime_settings.kinematics.externals.data]
momenta = [[500.0, 0.0, 0.0, 500.0], [500.0, 0.0, 0.0, -500.0]]
helicities = ["summed_averaged", "summed_averaged"]
"#,
    )?;

    let output = Command::new(env!("CARGO_BIN_EXE_gammaloop"))
        .current_dir(scratch.path())
        .args([
            "-n",
            "-l",
            "info",
            "-L",
            "off",
            "-p",
            "none",
            "run.toml",
            "run",
            "-c",
            "import model sm-default.json; \
             import graphs GL638.dot --process-spec 'e+ e- > t t~ h | e+ e- g t t~ h d d~ ghG ghG~ a QCD^2==4 QED^2==6 [{{4}} QCD=2]' -p epem_a_tth -i NNLO; \
             generate existing -p epem_a_tth -i NNLO",
        ])
        .output()?;
    // Successful generation is the regression expectation. This currently exposes
    // the error in DirectResidueBranches instead of treating that error as success.
    assert!(
        output.status.success(),
        "GL638 single-orientation integrated UV generation failed ({}).\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    Ok(())
}
