use super::utils::*;
use super::*;

mod important {
    use super::*;

    #[test]
    #[serial]
    fn cross_section_fermion_spin_sum_matches_explicit_incoming_helicity_average() -> Result<()> {
        let test_name = "cross_section_fermion_spin_sum_matches_explicit_incoming_helicity_average";
        let mut cli = setup_epem_tth_spin_sum_cli(test_name)?;
        let point = nontrivial_xspace_point_for(&mut cli, "epem_a_tth_sum", "LO")?;

        let summed = inspect_xspace_process(&mut cli, "epem_a_tth_sum", "LO", &point)?;
        let explicit_average = average_over_generated_helicity_processes(
            &mut cli,
            &[
                "epem_a_tth_pp",
                "epem_a_tth_pm",
                "epem_a_tth_mp",
                "epem_a_tth_mm",
            ],
            "LO",
            &point,
        )?;

        assert_complex_approx_eq(
            summed,
            explicit_average,
            "fermion incoming summed_averaged inspect should match the average over explicit helicity configurations",
        );

        clean_test(&cli.cli_settings.state.folder);
        Ok(())
    }
}

#[test]
#[serial]
fn cross_section_vector_spin_sum_matches_explicit_incoming_helicity_average() -> Result<()> {
    let test_name = "cross_section_vector_spin_sum_matches_explicit_incoming_helicity_average";
    let mut cli = setup_aa_ttbar_spin_sum_cli(test_name)?;
    let point = nontrivial_xspace_point_for(&mut cli, "aa_ttbar_sum", "LO")?;

    let summed = inspect_xspace_process(&mut cli, "aa_ttbar_sum", "LO", &point)?;
    let explicit_average = average_over_generated_helicity_processes(
        &mut cli,
        &["aa_ttbar_pp", "aa_ttbar_pm", "aa_ttbar_mp", "aa_ttbar_mm"],
        "LO",
        &point,
    )?;

    assert_complex_approx_eq(
        summed,
        explicit_average,
        "vector incoming summed_averaged inspect should match the average over explicit helicity configurations in light-like axial gauge",
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
