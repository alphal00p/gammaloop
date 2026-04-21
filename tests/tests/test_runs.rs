#![allow(dead_code)]
#![allow(unused_variables)]

pub(crate) use color_eyre::Result;
pub(crate) use colored::{ColoredString, Colorize};
pub(crate) use gammaloop_api::{
    commands::{
        Profile, Renormalize,
        evaluate_samples::{EvaluateSamples, evaluate_sample},
        inspect::Inspect,
        integrate::Integrate,
        profile::{InfraRedProfile, UltraVioletProfile},
    },
    state::{ProcessRef, SyncSettings},
};
pub(crate) use gammaloop_integration_tests::{
    clean_test, example_run_card, get_example_cli, get_test_cli, get_tests_workspace_path,
    run_commands, setup_sm_differential_lu_cli,
};
pub(crate) use gammalooprs::utils as gloop_utils;
pub(crate) use gammalooprs::{
    integrands::HasIntegrand,
    model::UFOSymbol,
    settings::runtime::{IntegralEstimate, IntegrationResult as RuntimeIntegrationResult},
    utils::{ApproxEq, F},
};
pub(crate) use insta::assert_snapshot;
pub(crate) use itertools::Itertools;
pub(crate) use momtrop::assert_approx_eq;
pub(crate) use ndarray::arr2;
pub(crate) use serde_json::Value as JsonValue;
pub(crate) use serial_test::serial;
pub(crate) use spenso::algebra::complex::Complex;
pub(crate) use std::path::{Path, PathBuf};
pub(crate) use std::{env, fmt::Write, time::Duration};
pub(crate) use symbolica::symbol;
pub(crate) use tabled::{Table, Tabled, settings::Style};
pub(crate) use tracing::info;

#[path = "test_runs/aa_aa.rs"]
mod aa_aa;
#[path = "test_runs/differential.rs"]
mod differential;
#[path = "test_runs/events.rs"]
mod events;
#[path = "test_runs/examples.rs"]
mod examples;
#[path = "test_runs/inspect.rs"]
mod inspect;
#[path = "test_runs/integrations.rs"]
mod integrations;
#[path = "test_runs/multi_integrand.rs"]
mod multi_integrand;
#[path = "test_runs/photonic.rs"]
mod photonic;
#[path = "test_runs/profile_bulk.rs"]
mod profile_bulk;
#[path = "test_runs/smoke.rs"]
mod smoke;
#[path = "test_runs/spin_sums.rs"]
mod spin_sums;
#[path = "test_runs/test_integrated_uv_cts.rs"]
mod test_integrated_uv_cts;
#[path = "test_runs/utils.rs"]
mod utils;
