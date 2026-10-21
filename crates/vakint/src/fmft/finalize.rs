//! Pure Symbolica finalization of an expression already in the FMFT PR basis.
//!
//! This is not an integral reducer: the caller supplies a reduced expression
//! with its classical mass dimension restored. No executable is consulted.

use colored::Colorize;
use log::debug;
use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::domains::rational::Rational;
use symbolica::function;

use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{FMFTOptions, Vakint, VakintError};

use super::FMFT;

impl FMFT {
    /// Finalize native scalar output in the existing Minkowski FMFT PR basis.
    ///
    /// Native and historical FORM callers share exact normalization and Laurent
    /// expansion, then the checked precision-tracking numerical boundary.
    /// No four-loop artifact is registered by exposing this internal boundary.
    #[allow(dead_code)] // Used by the optional experimental reduction bridge.
    pub(crate) fn finalize_native_reduced_masters(
        &self,
        evaluated_integral: Atom,
        muv_sq_atom: &Atom,
        options: &FMFTOptions,
    ) -> Result<Atom, VakintError> {
        crate::citations::USED_CITATIONS.record(crate::citations::CitationSource::Fmft);
        // The native reducer can return coefficients in d or the configured
        // epsilon variable. Normalize these before the Laurent expansion so
        // poles such as 1/(d-4) cannot hide required unknown master orders.
        let evaluated_integral = self.process_fmft_form_output(evaluated_integral)?;
        self.finalize_master_expression(evaluated_integral, 4, muv_sq_atom, options)
    }

    pub(super) fn finalize_master_expression(
        &self,
        evaluated_integral: Atom,
        loop_count: i64,
        muv_sq_atom: &Atom,
        options: &FMFTOptions,
    ) -> Result<Atom, VakintError> {
        let staged =
            self.stage_master_expression(evaluated_integral, loop_count, muv_sq_atom, options)?;
        if options.expand_masters && options.susbstitute_masters {
            crate::precision::materialize(
                staged.as_view(),
                &self.settings,
                &std::collections::HashMap::new(),
                None,
            )
        } else {
            Ok(staged)
        }
    }

    pub(super) fn stage_master_expression(
        &self,
        mut evaluated_integral: Atom,
        loop_count: i64,
        muv_sq_atom: &Atom,
        options: &FMFTOptions,
    ) -> Result<Atom, VakintError> {
        let settings = &self.settings;
        crate::citations::USED_CITATIONS.record(crate::citations::CitationSource::Fmft);
        let normalization = vk_parse!(
            format!(
                "((𝑖*(𝜋^((4-2*{eps})/2)))\
                *(exp(-EulerGamma))^({eps})\
                *(exp(-logmUVmu-log_mu_sq))^({eps}))^{loop_count}",
                eps = settings.epsilon_symbol,
            )
            .as_str()
        )
        .unwrap();

        // Adjust normalization factor.
        let complete_normalization = (normalization
            * settings
                .get_staged_integral_normalization_factor_atom()?
                .replace(S.n_loops.to_pattern())
                .with(Atom::num(loop_count).to_pattern()))
        .replace(Atom::var(vk_symbol!(settings.epsilon_symbol.as_str())).to_pattern())
        .with(vk_parse!("ep").unwrap().to_pattern());
        evaluated_integral *= complete_normalization;

        if options.expand_masters {
            let expansion_depth = settings.number_of_terms_in_epsilon_expansion - loop_count - 1;
            debug!(
                "{}: Expanding master integrals through {}^{} ...",
                "FMFT".green(),
                settings.epsilon_symbol,
                expansion_depth
            );
            evaluated_integral = self.expand_masters(evaluated_integral.as_view())?;
            evaluated_integral = evaluated_integral
                .series(
                    vk_symbol!("ep"),
                    Atom::Zero.as_view(),
                    Rational::from(expansion_depth),
                )
                .map_err(|error| VakintError::SymbolicaError(error.to_string()))?
                .to_atom();

            // Sanity check: coefficients can expose otherwise higher unknown
            // master orders through spurious epsilon poles.
            Self::reject_unknown_orders(&evaluated_integral)?;
            if options.susbstitute_masters {
                // Distribute exact Laurent coefficients before replacing PR
                // coefficients by approximate table values. Otherwise an
                // exact cancellation can become a tiny floating pole. This
                // is algebraic normalization, not a numerical zero threshold.
                evaluated_integral = evaluated_integral.expand();
                debug!(
                    "{}: Substituting master coefficients and period constants...",
                    "FMFT".green()
                );
                evaluated_integral = self.substitute_masters(evaluated_integral.as_view())?;
                evaluated_integral = evaluated_integral.expand();
                evaluated_integral = crate::hpl::normalize_constants(evaluated_integral.as_view())?;
                // Sanity check, including tails inside the numerical tables.
                Self::reject_unknown_orders(&evaluated_integral)?;
            }
        }

        evaluated_integral = evaluated_integral
            .replace(vk_parse!("ep").unwrap().to_pattern())
            .with(Atom::var(vk_symbol!(settings.epsilon_symbol.as_str())).to_pattern());
        if !settings.use_dot_product_notation {
            evaluated_integral = Vakint::convert_from_dot_notation(evaluated_integral.as_view());
        }
        let log_muv_mu_sq = function!(
            Symbol::LOG,
            muv_sq_atom.clone() / Atom::var(vk_symbol!(settings.mu_r_sq_symbol.as_str()))
        );
        let log_mu_sq = function!(
            Symbol::LOG,
            Atom::var(vk_symbol!(settings.mu_r_sq_symbol.as_str()))
        );
        evaluated_integral = evaluated_integral
            .replace(vk_parse!("logmUVmu").unwrap().to_pattern())
            .with(log_muv_mu_sq.to_pattern())
            .replace(vk_parse!("log_mu_sq").unwrap().to_pattern())
            .with(log_mu_sq.to_pattern());

        // println!(
        //     "evaluated_integral: {}",
        //     evaluated_integral.to_canonical_string()
        // );
        Ok(evaluated_integral)
    }

    fn reject_unknown_orders(expression: &Atom) -> Result<(), VakintError> {
        if let Some(matched) = expression
            .pattern_match(&vk_parse!("Oep(x_,y_)").unwrap().to_pattern(), None, None)
            .next()
        {
            return Err(VakintError::FMFTError(format!(
                "FMFT expansion yielded terms beyond expansion depth supported: Oep({},{})",
                matched.get(&S.x_).unwrap(),
                matched.get(&S.y_).unwrap(),
            )));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
