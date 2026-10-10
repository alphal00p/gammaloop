// Keep the NLO acceptances in one harness to share executable optimization and linking.

#[path = "test_nlo_acceptance/epem_a_ddx.rs"]
mod epem_a_ddx;

#[path = "test_nlo_acceptance/epem_a_ttx.rs"]
mod epem_a_ttx;

#[path = "test_nlo_acceptance/gamma_star_ddx.rs"]
mod gamma_star_ddx;

#[path = "test_nlo_acceptance/gamma_star_ttx.rs"]
mod gamma_star_ttx;
