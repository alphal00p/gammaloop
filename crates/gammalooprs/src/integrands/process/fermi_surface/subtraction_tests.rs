use color_eyre::Result;
use symbolica::atom::AtomCore;

use crate::{
    dot,
    graph::parse::IntoGraph,
    initialisation::test_initialise,
    integrands::process::param_builder::{ParamBuilderGraph, ThermalDistributionReplacement},
    numerator::symbolica_ext::NumeratorAtomExt,
    processes::AmplitudeGraph,
    settings::global::{GenerationSettings, MediumMode},
    utils::{GS, symbols::ThermalDistributionLimit},
};

use super::FermiSurfaceSector;

#[test]
fn fermi_sectors_retain_vacuum_subtraction_and_both_uv_contributions() -> Result<()> {
    test_initialise()?;
    // The repeated fermion pole supplies N'(E_d). A proper scalar tadpole
    // counterterm leaves that cycle in its thermal cograph, so both local and
    // integrated UV contributions must reach the same Fermi-sector boundary.
    // Every numerator is the constant on its own edge or vertex.
    let original: AmplitudeGraph = dot!(digraph fermi_cycle_with_uv_tadpole {
        node [num=1]
        edge [num=1]
        A -> B [id=0 particle="d" lmb_id=0]
        B -> A [id=1 particle="d"]
        A -> A [id=2 particle="H" lmb_id=1]
    })?;
    let mut generated = Vec::new();
    for (vacuum_subtraction, subtract_uv, generate_integrated) in [
        (false, false, false),
        (true, false, false),
        (true, true, false),
        (true, true, true),
    ] {
        let mut settings = GenerationSettings::default();
        settings.medium.mode = MediumMode::ZeroTemperatureEquilibrium;
        settings.medium.vacuum_subtraction = vacuum_subtraction;
        settings.threshold_subtraction.enable_thresholds = false;
        settings.explicit_orientation_sum_only = true;
        settings.uv.subtract_uv = subtract_uv;
        settings.uv.generate_integrated = generate_integrated;
        settings.uv.softct = false;
        let mut amplitude = original.clone();
        amplitude.generate_cff(&settings)?;
        amplitude.build_integrands(&settings, crate::utils::vakint()?)?;

        // Exercise the actual evaluator input before resolving any retained
        // numerator definitions. Sector collection keeps those bodies opaque.
        let expression = &amplitude.derived_data.all_mighty_integrand;
        let (bulk, sectors) = FermiSurfaceSector::extract(&amplitude.graph, expression)?;
        assert!(
            !sectors.is_empty(),
            "the generated cycle must have Fermi support"
        );
        assert!(!FermiSurfaceSector::is_present(&bulk)?);
        let reconstructed = sectors.iter().fold(bulk, |sum, sector| {
            sum + sector
                .product
                .factors()
                .iter()
                .zip(&sector.orientations)
                .fold(
                    sector.coefficient.clone(),
                    |coefficient, (factor, orientation)| {
                        coefficient
                            * GS.thermal_distribution(
                                factor.edge_id.0 as i64,
                                factor.derivative_order as i64,
                                0,
                                factor.sign,
                                *orientation,
                            )
                    },
                )
        });
        let difference = reconstructed - expression.unwrap_function(GS.thermal_weight_wrapper);
        assert!(
            difference.expand().is_zero(),
            "sector extraction changed UV={subtract_uv}, integrated={generate_integrated}: {difference}"
        );

        // These are scalar diagnostic copies. Only now resolve the definitions
        // to compare separately generated UV settings with distinct scope tags.
        generated.push(amplitude.derived_data.resolved_integrand()?);
    }

    let raw = &generated[0];
    let vacuum = original.graph.make_thermal_distributions_explicit(
        raw,
        ThermalDistributionLimit::Vacuum,
        original.graph.iter_edge_ids(),
        ThermalDistributionReplacement::All,
    )?;
    let difference = (&generated[1] - (raw - vacuum)).unwrap_function(GS.thermal_weight_wrapper);
    assert!(
        difference.expand().is_zero(),
        "vacuum subtraction must act on the complete distribution weight: {difference}"
    );

    for (name, contribution) in [
        ("local UV", &generated[2] - &generated[1]),
        ("integrated UV", &generated[3] - &generated[2]),
    ] {
        let (_, sectors) = FermiSurfaceSector::extract(&original.graph, &contribution)?;
        assert!(
            sectors
                .iter()
                .any(|sector| !sector.coefficient.expand().is_zero()),
            "the {name} contribution must retain the untouched fermion cycle's Fermi delta"
        );
    }
    Ok(())
}
