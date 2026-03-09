use std::collections::{BTreeMap, HashMap};

use crate::fmft_numerics::{MASTERS_EXPANSION, MASTERS_NUMERIC_SUBSTITUTIONS};
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::utils::{self, set_precision_in_polynomial_atom, undress_vakint_symbols};
use crate::utils::{get_full_name, set_precision_in_float_atom};
use crate::{
    fmft_numerics::{ADDITIONAL_CONSTANTS, POLY_GAMMA_SUBSTITUTIONS},
    gt_condition,
};
use colored::Colorize;
use log::debug;
use regex::Regex;
use string_template_plus::{Render, RenderOptions, Template};
use symbolica::atom::Symbol;
use symbolica::printer::{AtomPrinter, PrintOptions};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::{integer::Integer, rational::Rational},
    function,
    id::Condition,
};

use crate::{
    FMFTOptions, TEMPLATES, get_integer_from_atom, number_condition, symbol_condition, symbols::S,
};

use crate::{ReplacementRules, Vakint, VakintError, VakintSettings};

pub struct FMFT {
    settings: VakintSettings,
}

impl FMFT {
    pub fn with_settings(settings: VakintSettings) -> Self {
        FMFT { settings }
    }

    pub fn process_fmft_form_output(
        &self,
        processed_form_output: Atom,
    ) -> Result<Atom, VakintError> {
        let mut res = processed_form_output.to_owned();
        res = res
            .replace(vk_parse!("d").unwrap().to_pattern())
            .with(vk_parse!("4-2*ep").unwrap().to_pattern());
        // Temporarily work with the variable "ep" instead of the UTF-8 symbol for epsilon
        res = res
            .replace(Atom::var(vk_symbol!(self.settings.epsilon_symbol.as_str())).to_pattern())
            .with(vk_parse!("ep").unwrap().to_pattern());
        res = res
            .replace(
                function!(S.vkdot, function!(S.p, S.id1_a), function!(S.p, S.id2_a)).to_pattern(),
            )
            .when(
                &(Condition::from((S.id1_, number_condition()))
                    & Condition::from((S.id2_, number_condition()))),
            )
            .with(function!(S.dot, function!(S.p, S.id1_a), function!(S.p, S.id2_a)).to_pattern());
        Ok(res)
    }

    pub fn substitute_gam_functions(&self, result: AtomView) -> Atom {
        let mut res = result.to_owned();
        res = res
            .replace(vk_parse!("Gam(x_,y_)").unwrap().to_pattern())
            .with(
                vk_parse!("exp(ep*y_*EulerGamma)*Gamma(x_+ep*y_)")
                    .unwrap()
                    .to_pattern(),
            );
        res = res
            .replace(vk_parse!("iGam(x_,y_)").unwrap().to_pattern())
            .with(
                vk_parse!("exp(-ep*y_*EulerGamma)/Gamma(x_+ep*y_)")
                    .unwrap()
                    .to_pattern(),
            );
        res
    }

    pub fn expand_masters(&self, result: AtomView) -> Result<Atom, VakintError> {
        let mut r = result.to_owned();
        r.repeat_map(Box::new(move |av: AtomView| {
            let mut res = av.to_owned();
            res = self.substitute_gam_functions(res.as_view());
            for (src, (trgt, restriction)) in MASTERS_EXPANSION.iter() {
                res = res
                    .replace(src.to_pattern())
                    .when(restriction)
                    .with(trgt.to_pattern());
            }
            res
        }));
        if let Some(m) = r
            .pattern_match(&vk_parse!("GammaArgs(x_,y_)").unwrap().into(), None, None)
            .next()
        {
            return Err(VakintError::FMFTError(format!(
                "FMFT result contains a Gamma function whose numerical evaluation is not implemented in vakint: Gamma({}+{}*{})",
                m.get(&S.x_).unwrap(),
                m.get(&S.y_).unwrap(),
                self.settings.epsilon_symbol
            )));
        }
        Ok(r)
    }

    pub fn substitute_masters(&self, result: AtomView) -> Result<Atom, VakintError> {
        let processed_constants = MASTERS_NUMERIC_SUBSTITUTIONS
            .iter()
            .map(|(src, (trgt, condition))| {
                (
                    src,
                    (
                        set_precision_in_polynomial_atom(
                            trgt.as_view(),
                            vk_symbol!("ep"),
                            &self.settings,
                        ),
                        condition.clone(),
                    ),
                )
            })
            .collect::<Vec<_>>();
        // for (a, (b, _c)) in processed_constants.clone() {
        //     println!("{} -> {}", a, b);
        // }
        let mut r = result.to_owned();
        r.repeat_map(Box::new(move |av: AtomView| {
            let mut res = av.to_owned();
            for (src, (trgt, matching_condition)) in processed_constants.iter() {
                res = res
                    .replace(src.to_pattern())
                    .when(matching_condition)
                    .with(trgt.to_pattern());
            }
            res
        }));
        // println!("DONE! {}", r);
        Ok(r)
    }

    pub fn substitute_poly_gamma(&self, result: AtomView) -> Result<Atom, VakintError> {
        let processed_constants = POLY_GAMMA_SUBSTITUTIONS
            .iter()
            .map(|(src, trgt)| {
                (
                    src,
                    set_precision_in_float_atom(trgt.as_view(), &self.settings),
                )
            })
            .collect::<Vec<_>>();
        let mut r = result.to_owned();

        r = r
            .replace(vk_parse!("Gamma(n_)").unwrap().to_pattern())
            .when(Condition::from((S.n_, gt_condition(0))))
            .with_map(move |match_in| {
                let n = get_integer_from_atom(match_in.get(S.n_).unwrap().to_atom().as_view())
                    .unwrap() as u32;
                Atom::num(Integer::factorial(n - 1))
            });
        r.repeat_map(Box::new(move |av: AtomView| {
            let mut res = av.to_owned();
            for (src, trgt) in processed_constants.iter() {
                res = res.replace(src.to_pattern()).with(trgt.to_pattern());
            }
            res
        }));

        if let Some(m) = r
            .pattern_match(&vk_parse!("PolyGamma(x_,y_)").unwrap().into(), None, None)
            .next()
        {
            return Err(VakintError::FMFTError(format!(
                "FMFT result contains a PolyGamma function whose numerical evaluation is not implemented in vakint: PolyGamma({},{})",
                m.get(&S.x_).unwrap(),
                m.get(&S.y_).unwrap()
            )));
        }
        Ok(r)
    }

    pub fn substitute_additional_constants(&self, result: AtomView) -> Result<Atom, VakintError> {
        let processed_constants = ADDITIONAL_CONSTANTS
            .iter()
            .map(|(src, trgt)| {
                (
                    src,
                    set_precision_in_float_atom(trgt.as_view(), &self.settings),
                )
            })
            .collect::<Vec<_>>();
        let mut r = result.to_owned();
        r.repeat_map(Box::new(move |av: AtomView| {
            let mut res = av.to_owned();
            for (src, trgt) in processed_constants.iter() {
                res = res.replace(src.to_pattern()).with(trgt.to_pattern());
            }
            res
        }));

        Ok(r)
    }
}

impl Vakint {
    pub fn fmft_evaluate(
        &self,
        settings: &VakintSettings,
        input_numerator: AtomView,
        integral_specs: &ReplacementRules,
        options: &FMFTOptions,
    ) -> Result<Atom, VakintError> {
        let integral = integral_specs.canonical_topology.get_integral();

        debug!(
            "Processing the following integral with {}:\n{}",
            "FMFT".green(),
            integral
        );

        let err = Err(VakintError::InvalidGenericExpression(format!(
            "Could not find the shorthand integral name in the FMFT expression: {}",
            integral
                .short_expression
                .as_ref()
                .map(|a| a.to_canonical_string())
                .unwrap_or("None".to_string())
        )));

        let fmft = FMFT::with_settings(settings.clone());

        let integral_name = if let Some(short_expression) = integral.short_expression.as_ref() {
            if let Some(m) = short_expression
                .pattern_match(
                    &function!(S.fun_, S.any_a___).to_pattern(),
                    Some(&Condition::from((S.fun_, symbol_condition()))),
                    None,
                )
                .next()
            {
                AtomPrinter::new_with_options(
                    m.get(&S.fun_).unwrap().as_atom_view(),
                    PrintOptions::file_no_namespace(),
                )
                .to_string()
            } else {
                return err;
            }
        } else {
            return err;
        };

        let (muv_atom, muv_sq_atom) =
            Vakint::identify_uv_mass_symbols(integral.canonical_expression.as_ref().unwrap())?;

        // Here we map the propagators in the correct order for the definition of the topology in FMFT
        let vakint_to_fmft_edge_map = match integral_name.as_str().split("_pinch_").next().unwrap()
        {
            "I4L_H" => vec![2, 3, 4, 5, 6, 7, 8, 9, 1],
            "I4L_X" => vec![2, 3, 4, 5, 6, 7, 8, 9, 10],
            "I4L_BMW" => vec![3, 4, 5, 6, 7, 8, 9, 10],
            "I4L_FG" => vec![1, 3, 4, 5, 6, 7, 8, 9],
            _ => {
                return Err(VakintError::InvalidGenericExpression(format!(
                    "Integral {} is not supported by FMFT.",
                    integral_name
                )));
            }
        };

        let mut numerator = Vakint::convert_to_dot_notation(input_numerator);

        if utils::could_match(
            &function!(S.dot, function!(S.p, S.id1_a), function!(S.k, S.id2_a)).to_pattern(),
            numerator.as_view(),
        ) {
            return Err(VakintError::InvalidNumerator(format!(
                "Make sure the numerator has been tensor-reduced before being processed by FMFT : {}",
                numerator
            )));
        }

        // Now map all exterior dot products into a special function `vkdot` so that it does not interfere with FMFT
        // Do not forget to normalize by the dimensionality, which is muv^2 in this case.
        numerator = numerator
            .replace(
                function!(S.dot, function!(S.p, S.id1_a), function!(S.p, S.id2_a)).to_pattern(),
            )
            .when(
                &(Condition::from((S.id1_, number_condition()))
                    & Condition::from((S.id2_, number_condition()))),
            )
            .with(
                function!(S.vkdot, function!(S.p, S.id1_a), function!(S.p, S.id2_a)).to_pattern(),
            );

        // And finally map all interior products into the form p<i>.p<i> expected by fmft, where <i> is the edge id carrying momentum k<i>.
        let momenta: HashMap<usize, Atom> = integral_specs.get_propagator_property_list("q_");
        let mut lmb_prop_indices = vec![];
        for i_loop in 1..=integral.n_loops {
            if let Some((i_edge, _)) = momenta.iter().find(|(_i_prop, k)| {
                k.as_view() == function!(S.k, Atom::num(i_loop as i64)).as_view()
            }) {
                lmb_prop_indices.push(*i_edge as isize);
            } else {
                lmb_prop_indices.push(-1);
            }
        }

        let vakint_to_fmft_edge_map_copy = vakint_to_fmft_edge_map.clone();
        let muv_sq_atom_clone = muv_sq_atom.clone();
        numerator = numerator.replace(
                function!(S.dot, function!(S.k, S.id1_a), function!(S.k, S.id2_a))
            .to_pattern()).when(&(Condition::from((S.id1_, number_condition()))
            & Condition::from((S.id2_, number_condition())))).with_map(
                move |match_in| {
                    let id1 = lmb_prop_indices[(get_integer_from_atom(match_in.get(S.id1_).unwrap().to_atom().as_view()).unwrap()-1) as usize];
                    if id1 < 0 {
                        panic!(
                            "Could not find LMB edge for momentum k({}) in a topology supported by FMFT and used in numerator.",
                            get_integer_from_atom(match_in.get(S.id1_).unwrap().to_atom().as_view()).unwrap()
                        );
                    }
                    let id2 = lmb_prop_indices[(get_integer_from_atom(match_in.get(S.id2_).unwrap().to_atom().as_view()).unwrap()-1) as usize];
                    if id2 < 0 {
                        panic!(
                            "Could not find LMB edge for momentum k({}) in a topology supported by FMFT and used in numerator.",
                            get_integer_from_atom(match_in.get(S.id2_).unwrap().to_atom().as_view()).unwrap()
                        );
                    }
                    let i_edge1 =
                        vakint_to_fmft_edge_map_copy[(id1 as usize) - 1];
                    let i_edge2 =
                        vakint_to_fmft_edge_map_copy[(id2 as usize) - 1];
                    // in FMFT, the loop momenta dot products need to be written p<i>.p<i>
                    // The outter dot will be converted to an inner dot in the next step
                    // Again we must normalize by the dimensionality muv^2
                    vk_parse!(format!("dot(p{},p{})", i_edge1, i_edge2).as_str()).unwrap()
                        * (muv_sq_atom_clone.clone())
                }
            );

        // Substitute eps by (4-d)/2
        numerator = numerator
            .replace(vk_parse!(&settings.epsilon_symbol).unwrap().to_pattern())
            .with(vk_parse!("(4-d)/2").unwrap().to_pattern());

        let muv_symbol = match (muv_atom.as_view(), muv_sq_atom.as_view()) {
            (AtomView::Var(s), _) => s.get_symbol(),
            (_, AtomView::Var(s2)) => s2.get_symbol(),
            _ => {
                return Err(VakintError::MalformedGraph(
                    "Could not identify the UV mass symbol in the integral.".into(),
                ));
            }
        };

        let (form_header_additions, expression_str, indices) =
            self.sanitize_user_expressions(settings, numerator.as_view(), true, &[muv_symbol])?;

        // let expression_str =
        //     &AtomPrinter::new_with_options(numerator.as_view(), PrintOptions::file_no_namespace())
        //         .to_string();

        let dot_produce_replacer =
            Regex::new(r"dot\((?<vecA>[\w|\d]+),(?<vecB>[\w|\d]+)\)").unwrap();
        let numerator_string = dot_produce_replacer
            .replace_all(&expression_str, "($vecA.$vecB)")
            .to_string();

        let powers = (1..=integral.n_props)
            .map(|i_prop| {
                integral_specs
                    .canonical_expression_substitutions
                    .get(&function!(S.pow, Atom::num(i_prop as i64)))
                    .map(|a| a.try_into().unwrap())
                    .unwrap_or(0)
            })
            .collect::<Vec<_>>();

        let integral_string = powers
            .iter()
            .zip(vakint_to_fmft_edge_map)
            .map(|(pwr, fmft_edge_index)| format!("d{}^{}", fmft_edge_index, -pwr))
            .collect::<Vec<_>>()
            .join("*");

        // Replace functions with 1 and get all remaining symbols
        // let mut numerator_additional_symbols = input_numerator
        //     .replace(vk_parse!("f_(args__)").unwrap().to_pattern())
        //     .with(vk_parse!("1").unwrap().to_pattern())
        //     .get_all_symbols(false);
        // let eps_symbol = vk_symbol!(settings.epsilon_symbol.clone());
        // numerator_additional_symbols.retain(|&s| s != eps_symbol);
        let numerator_additional_symbols: std::collections::HashSet<
            symbolica::atom::Symbol,
            ahash::RandomState,
        > = std::collections::HashSet::default();

        let template = Template::parse_template(TEMPLATES.get("run_fmft.txt").unwrap()).unwrap();
        let mut vars: HashMap<String, String> = HashMap::new();
        vars.insert("numerator".into(), numerator_string);
        vars.insert("integral".into(), integral_string);
        // match (muv_atom.as_view(), muv_sq_atom.as_view()) {
        //     (AtomView::Var(s), _) => vars.insert("symbols".into(), s.get_symbol().to_string()),
        //     (_, AtomView::Var(s2)) => vars.insert("symbols".into(), s2.get_symbol().to_string()),
        //     _ => {
        //         return Err(VakintError::MalformedGraph(
        //             "Could not identify the UV mass symbol in the integral.".into(),
        //         ));
        //     }
        // };
        vars.insert("symbols".into(), "NOPREDEDUBEDVAKINTSYMBOL".into());

        let mut additional_symbols_str = if !numerator_additional_symbols.is_empty() {
            format!(
                "Auto S {};\n",
                numerator_additional_symbols
                    .iter()
                    .map(|item| item.get_stripped_name())
                    .collect::<Vec<_>>()
                    .join(", "),
            )
        } else {
            "".into()
        };
        additional_symbols_str = format!(
            "{}\nCF {};",
            additional_symbols_str,
            undress_vakint_symbols(&get_full_name(&S.g))
        );
        if !form_header_additions.is_empty() {
            if !additional_symbols_str.is_empty() {
                additional_symbols_str =
                    format!("{}\n{}", form_header_additions, additional_symbols_str);
            } else {
                additional_symbols_str = form_header_additions;
            }
        }
        vars.insert("additional_symbols".into(), additional_symbols_str);

        let rendered = template
            .render(&RenderOptions {
                variables: vars,
                ..Default::default()
            })
            .unwrap();

        let form_result = self.run_form(
            settings,
            &["fmft.frm".into()],
            ("run_fmft.frm".into(), rendered),
            vec![],
            settings.clean_tmp_dir,
            settings.temporary_directory.clone(),
        )?;

        let processed_form_result =
            self.process_form_output(settings, form_result, indices, BTreeMap::new())?;
        let mut evaluated_integral = fmft.process_fmft_form_output(processed_form_result)?;
        debug!(
            "{}: raw result from FORM:\n{}",
            "FMFT".green(),
            evaluated_integral
        );

        // Restore dimensionality now

        // Offset dimensionality by the denominators
        let muv_sq_dimension = 2 * (integral.n_loops as i64) - powers.iter().sum::<i64>();

        evaluated_integral *= vk_parse!(
            format!(
                "({})^{}",
                muv_sq_atom.to_canonical_string(),
                muv_sq_dimension
            )
            .as_str()
        )
        .unwrap();

        let fmft_normalization_correction = vk_parse!(
            format!(
                "(
                (𝑖*(𝜋^((4-2*{eps})/2)))\
              * (exp(-EulerGamma))^({eps})\
              * (exp(-logmUVmu-log_mu_sq))^({eps})\
             )^{n_loops}",
                eps = settings.epsilon_symbol,
                n_loops = integral.n_loops
            )
            .as_str()
        )
        .unwrap();

        // Adjust normalization factor
        let mut complete_normalization = fmft_normalization_correction
            * settings
                .get_integral_normalization_factor_atom()?
                .replace(S.n_loops.to_pattern())
                .with(Atom::num(integral.n_loops as i64).to_pattern());
        complete_normalization = complete_normalization
            .replace(Atom::var(vk_symbol!(settings.epsilon_symbol.as_str())).to_pattern())
            .with(vk_parse!("ep").unwrap().to_pattern());

        evaluated_integral *= complete_normalization;

        if options.expand_masters {
            let expansion_depth =
                settings.number_of_terms_in_epsilon_expansion - (integral.n_loops as i64) - 1;
            debug!(
                "{}: Expanding master integrals with terms up to and including {}^{} ...",
                "FMFT".green(),
                settings.epsilon_symbol,
                expansion_depth
            );
            evaluated_integral = fmft.expand_masters(evaluated_integral.as_view())?;

            debug!(
                "{}: Series expansion of the result up to and including terms of order {}^{} ...",
                "FMFT".green(),
                settings.epsilon_symbol,
                expansion_depth
            );
            evaluated_integral = match evaluated_integral.series(
                vk_symbol!("ep"),
                Atom::Zero.as_view(),
                Rational::from(expansion_depth),
                true,
            ) {
                Ok(a) => a,
                Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
            }
            .to_atom();

            // Sanity check
            if let Some(m) = evaluated_integral
                .pattern_match(&vk_parse!("Oep(x_,y_)").unwrap().to_pattern(), None, None)
                .next()
            {
                return Err(VakintError::FMFTError(format!(
                    "FMFT expansion yielded terms beyond expansion depth supported: Oep({},{})",
                    m.get(&S.x_).unwrap(),
                    m.get(&S.y_).unwrap(),
                )));
            }

            if options.susbstitute_masters {
                debug!(
                    "{}: Substituting master integrals coefficient with their numerical evaluations...",
                    "FMFT".green()
                );
                evaluated_integral = fmft.substitute_masters(evaluated_integral.as_view())?;
                debug!(
                    "{}: Substituting PolyGamma and period constants...",
                    "FMFT".green()
                );
                evaluated_integral = evaluated_integral.expand();
                evaluated_integral = fmft.substitute_poly_gamma(evaluated_integral.as_view())?;
                evaluated_integral =
                    fmft.substitute_additional_constants(evaluated_integral.as_view())?;
                // Sanity check
                if let Some(m) = evaluated_integral
                    .pattern_match(&vk_parse!("Oep(x_,y_)").unwrap().to_pattern(), None, None)
                    .next()
                {
                    return Err(VakintError::FMFTError(format!(
                        "FMFT expansion yielded terms beyond expansion depth supported: Oep({},{})",
                        m.get(&S.x_).unwrap(),
                        m.get(&S.y_).unwrap(),
                    )));
                }
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
            muv_sq_atom / Atom::var(vk_symbol!(settings.mu_r_sq_symbol.as_str()))
        );

        let log_mu_sq = function!(
            Symbol::LOG,
            Atom::var(vk_symbol!(settings.mu_r_sq_symbol.as_str()))
        );

        evaluated_integral = evaluated_integral
            .replace(vk_parse!("logmUVmu").unwrap().to_pattern())
            .with((log_muv_mu_sq).to_pattern());
        evaluated_integral = evaluated_integral
            .replace(vk_parse!("log_mu_sq").unwrap().to_pattern())
            .with((log_mu_sq).to_pattern());

        // println!(
        //     "evaluated_integral: {}",
        //     evaluated_integral.to_canonical_string()
        // );

        Ok(evaluated_integral)
    }
}
