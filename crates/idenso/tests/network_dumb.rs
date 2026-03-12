use std::time::Instant;

use idenso::rep_symbols::RS;
use idenso::{metric::MetricSimplifier, representations::initialize};
use spenso::network::library::symbolic::ETS;
use spenso::{
    network::parsing::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::atom::FunctionBuilder;
use symbolica::{
    LicenseManager,
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::{Match, MatchSettings, MatchStack, Pattern, Replacement},
    parse, symbol,
    transformer::Transformer,
};

use symbolica::{parse_lit, tag};

fn run_network_informed() {
    initialize();
    let d = symbol!("d"; Symmetric, Linear); // metric

    let mink = Minkowski {}.new_rep(4);
    let rep = tag!("rep");

    let k0 = symbol!("k0");
    let muw1_ = symbol!(
        "muw1_",
        tags = [&SPENSO_TAG.representation, &SPENSO_TAG.self_dual]
    );
    // let muw1_ = symbol!("muw1_", tags = [&rep]);
    let k1_ = symbol!("k1_");
    let x___ = symbol!("x___");

    symbol!("mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "mu10", "mu11";
        tags=["spenso::index"]);

    let index = tag!("index");

    let mink = symbol!(
        "mink",
        norm = move |ain, out| {
            let AtomView::Fun(f) = ain else {
                return;
            };

            // println!("Normalizing representation function: {:#}", ain);

            let name = f.get_symbol();

            if f.get_nargs() == 3 {
                let mut args = f.iter();
                let dim = args.next().unwrap();
                let a = args.next().unwrap();
                let b = args.next().unwrap();

                // A rep is linear in each index slot;
                match a {
                    AtomView::Add(a) => {
                        let builder = FunctionBuilder::new(name);

                        let mut new_out = Atom::Zero;
                        for t in a.iter() {
                            new_out += builder.clone().add_args(&[dim, t, b]).finish();
                        }
                        **out = new_out;
                        return;
                    }

                    AtomView::Mul(a) => {
                        let builder = FunctionBuilder::new(name);

                        let mut new_out = Atom::num(1);
                        let mut coef = Atom::num(1);
                        for t in a.iter() {
                            if let AtomView::Num(_) = t {
                                coef *= t;
                            } else {
                                new_out *= t;
                            }
                        }

                        if coef != Atom::num(1) {
                            **out = coef * builder.add_args(&[dim, new_out.as_view(), b]).finish();
                            return;
                        }
                    }

                    _ => {}
                }

                match b {
                    AtomView::Add(b) => {
                        let builder = FunctionBuilder::new(name);

                        let mut new_out = Atom::Zero;
                        for t in b.iter() {
                            new_out += builder.clone().add_args(&[dim, a, t]).finish();
                        }
                        **out = new_out;
                        return;
                    }

                    AtomView::Mul(b) => {
                        let builder = FunctionBuilder::new(name);

                        let mut new_out = Atom::num(1);
                        let mut coef = Atom::num(1);
                        for t in b.iter() {
                            if let AtomView::Num(_) = t {
                                coef *= t;
                            } else {
                                new_out *= t;
                            }
                        }

                        if coef != Atom::num(1) {
                            **out = coef * builder.add_args(&[dim, a, new_out.as_view()]).finish();
                            return;
                        }
                    }

                    AtomView::Pow(_) => {
                        panic!("Powers are not supported in index slots {ain}");
                    }
                    _ => {}
                }

                let a_build = match a {
                    AtomView::Fun(f) => {
                        let sym = f.get_symbol();
                        if sym.get_wildcard_level() > 0 {
                            return;
                        }

                        if sym.has_tag(&index) {
                            None
                        } else {
                            Some((sym, f.iter().collect::<Vec<_>>()))
                        }
                    }
                    AtomView::Var(v) => {
                        let sym = v.get_symbol();
                        if sym.get_wildcard_level() > 0 {
                            return;
                        }
                        if sym.has_tag(&index) {
                            None
                        } else {
                            Some((sym, vec![]))
                        }
                    }
                    _ => None,
                };

                let b_build = match b {
                    AtomView::Fun(f) => {
                        let sym = f.get_symbol();
                        if sym.get_wildcard_level() > 0 {
                            return;
                        }
                        if sym.has_tag(&index) {
                            None
                        } else {
                            Some((sym, f.iter().collect::<Vec<_>>()))
                        }
                    }
                    AtomView::Var(v) => {
                        let sym = v.get_symbol();
                        if sym.get_wildcard_level() > 0 {
                            return;
                        }
                        // sym
                        if sym.has_tag(&index) {
                            None
                        } else {
                            Some((sym, vec![]))
                        }
                    }
                    _ => None,
                };

                match (a_build, b_build) {
                    (Some((a, args)), Some((b, args_b))) => {
                        // println!(
                        //     "Both are functions or variables, so this is a dot product;"
                        // );
                        let builder = FunctionBuilder::new(d);
                        let rep = FunctionBuilder::new(name).add_arg(dim).finish();

                        let a = if args.is_empty() {
                            Atom::var(a)
                        } else {
                            FunctionBuilder::new(a).add_args(&args).finish()
                        };

                        let b = if args_b.is_empty() {
                            Atom::var(b)
                        } else {
                            FunctionBuilder::new(b).add_args(&args_b).finish()
                        };

                        let atom = builder.add_args(&[rep, a, b]).finish();
                        // println!("{atom}");
                        **out = atom;
                    }
                    (None, Some((b, args))) => {
                        // println!("a is an index and b is a tensor, we just append to b");
                        let index = FunctionBuilder::new(name).add_arg(dim).add_arg(a).finish();
                        // println!("index: {index}");
                        let atom = FunctionBuilder::new(b)
                            .add_args(&args)
                            .add_arg(index)
                            .finish();
                        // println!("atom:{atom}");
                        **out = atom;
                    }
                    (Some((a, args)), None) => {
                        // println!(
                        //     "b is a dual index (in second position) and a is a tensor, we just append to a"
                        // );
                        let index = FunctionBuilder::new(name).add_arg(dim).add_arg(b).finish();
                        // println!("index: {index}");
                        let atom = FunctionBuilder::new(a)
                            .add_args(&args)
                            .add_arg(index)
                            .finish();
                        // println!("atom: {atom}");
                        **out = atom;
                    }
                    (None, None) => {
                        // println!("both are indices so this is a metric shorthand");
                        let index_builder = FunctionBuilder::new(name).add_arg(dim);

                        let atom = FunctionBuilder::new(d)
                            .add_args(&[
                                index_builder.clone().add_arg(a).finish(),
                                index_builder.add_arg(b).finish(),
                            ])
                            .finish();

                        // println!("atom:{atom}");
                        **out = atom;
                    }
                }
            }
            // else {
            //     // println!("Does not have 3 arguments: {}", f.get_nargs());
            // }
        },
        tags = [&rep]
    );

    let mut r = parse!(
        "vx(1,-k(0), k(0)-k(1), k(1), k(10), mu1, mu8)
    vx(2,-k(1), k(2), k(1)-k(2), mu1, mu2, mu9)
    vx(3,-k(2), k(3), k(2)-k(3), mu2, mu3, mu10)
    vx(4,-k(3), k(4), k(3)-k(4), mu3, mu4, mu11)
    vx(5,-k(4), k(0), k(4)-k(0), mu4, k(20), mu5)
    vx(6,-k(4)+k(0), -k(3)+k(4), k(3)-k(0), mu5, mu11, mu6)
    vx(7,-k(3)+k(0), -k(2)+k(3), k(2)-k(0), mu6, mu10, mu7)
    vx(8,-k(2)+k(0), -k(1)+k(2), k(1)-k(0), mu7, mu9, mu8)"
    );

    let gluon_rule = parse!(
        "(- spenso::mink(4,mu1_, mu3_) * spenso::mink(4,k1_, mu2_)
                + spenso::mink(4,mu1_, mu2_) * spenso::mink(4,k1_, mu3_)
                + spenso::mink(4,mu2_, mu3_) * spenso::mink(4,k2_, mu1_)
                - spenso::mink(4,mu1_, mu2_) * spenso::mink(4,k2_, mu3_)
                - spenso::mink(4,mu2_, mu3_) * spenso::mink(4,k3_, mu1_)
                + spenso::mink(4,mu1_, mu3_) * spenso::mink(4,k3_, mu2_)
                )"
    )
    .to_pattern();

    let rhs_subs = Pattern::Transformer(Box::new((
        Some(gluon_rule),
        vec![Transformer::Map(Box::new(move |input, _state, out| {
            let mut exp = input.expand();
            println!("Expanded");
            exp = exp.simplify_metrics();
            println!("Simplified");
            exp = exp.to_dots();
            // println!("To dots", exp);

            *out = exp;
            Ok(())
        }))],
    )));

    let rhs_subs = {
        let gluon_rule = parse!(
            "(- spenso::mink(4,mu1_, mu3_) * spenso::mink(4,k1_, mu2_)
                    + spenso::mink(4,mu1_, mu2_) * spenso::mink(4,k1_, mu3_)
                    + spenso::mink(4,mu2_, mu3_) * spenso::mink(4,k2_, mu1_)
                    - spenso::mink(4,mu1_, mu2_) * spenso::mink(4,k2_, mu3_)
                    - spenso::mink(4,mu2_, mu3_) * spenso::mink(4,k3_, mu1_)
                    + spenso::mink(4,mu1_, mu3_) * spenso::mink(4,k3_, mu2_)
                    )"
        )
        .to_pattern();

        let rhs = move |m: &MatchStack| {
            let mut input = gluon_rule.replace_wildcards_with_matches(m);

            input = input.expand();
            // input = input.simplify_metrics();
            // input = input.to_dots();
            input
        };

        rhs
    };

    // r = parse_lit!(
    //     vx(3, -k(2), k(3), k(2) - k(3), mu2, mu3, mu10)
    //         * vx(4, -k(3), k(4), k(3) - k(4), mu3, mu4, mu11)
    //         * vx(2, -k(1), k(2), k(1) - k(2), k(0), mu2, mu9)
    //         * vx(5, -k(4), k(0), -k(0) + k(4), mu4, k(20), mu5)
    //         * vx(6, k(0) - k(4), -k(3) + k(4), -k(0) + k(3), mu5, mu11, mu6)
    //         * vx(7, k(0) - k(3), -k(2) + k(3), -k(0) + k(2), mu6, mu10, mu7)
    //         * vx(8, k(0) - k(2), -k(1) + k(2), -k(0) + k(1), mu7, mu9, k(10))
    //         + vx(3, -k(2), k(3), k(2) - k(3), mu2, mu3, mu10)
    //             * vx(4, -k(3), k(4), k(3) - k(4), mu3, mu4, mu11)
    //             * vx(2, -k(1), k(2), k(1) - k(2), k(1), mu2, mu9)
    //             * vx(5, -k(4), k(0), -k(0) + k(4), mu4, k(20), mu5)
    //             * vx(6, k(0) - k(4), -k(3) + k(4), -k(0) + k(3), mu5, mu11, mu6)
    //             * vx(7, k(0) - k(3), -k(2) + k(3), -k(0) + k(2), mu6, mu10, mu7)
    //             * vx(8, k(0) - k(2), -k(1) + k(2), -k(0) + k(1), mu7, mu9, k(10))
    //         - 2 * vx(3, -k(2), k(3), k(2) - k(3), mu2, mu3, mu10)
    //             * vx(4, -k(3), k(4), k(3) - k(4), mu3, mu4, mu11)
    //             * vx(2, -k(1), k(2), k(1) - k(2), k(10), mu2, mu9)
    //             * vx(5, -k(4), k(0), -k(0) + k(4), mu4, k(20), mu5)
    //             * vx(6, k(0) - k(4), -k(3) + k(4), -k(0) + k(3), mu5, mu11, mu6)
    //             * vx(7, k(0) - k(3), -k(2) + k(3), -k(0) + k(2), mu6, mu10, mu7)
    //             * vx(8, k(0) - k(2), -k(1) + k(2), -k(0) + k(1), mu7, mu9, k(0))
    //         + vx(3, -k(2), k(3), k(2) - k(3), mu2, mu3, mu10)
    //             * vx(4, -k(3), k(4), k(3) - k(4), mu3, mu4, mu11)
    //             * vx(2, -k(1), k(2), k(1) - k(2), k(10), mu2, mu9)
    //             * vx(5, -k(4), k(0), -k(0) + k(4), mu4, k(20), mu5)
    //             * vx(6, k(0) - k(4), -k(3) + k(4), -k(0) + k(3), mu5, mu11, mu6)
    //             * vx(7, k(0) - k(3), -k(2) + k(3), -k(0) + k(2), mu6, mu10, mu7)
    //             * vx(8, k(0) - k(2), -k(1) + k(2), -k(0) + k(1), mu7, mu9, k(1))
    //         + spenso::g(spenso::mink(4), k(0), k(10))
    //             * vx(3, -k(2), k(3), k(2) - k(3), mu2, mu3, mu10)
    //             * vx(4, -k(3), k(4), k(3) - k(4), mu3, mu4, mu11)
    //             * vx(5, -k(4), k(0), -k(0) + k(4), mu4, k(20), mu5)
    //             * vx(2, -k(1), k(2), k(1) - k(2), mink(4, mu8), mu2, mu9)
    //             * vx(6, k(0) - k(4), -k(3) + k(4), -k(0) + k(3), mu5, mu11, mu6)
    //             * vx(7, k(0) - k(3), -k(2) + k(3), -k(0) + k(2), mu6, mu10, mu7)
    //             * vx(8, k(0) - k(2), -k(1) + k(2), -k(0) + k(1), mu7, mu9, mu8)
    //         - 2 * spenso::g(spenso::mink(4), k(1), k(10))
    //             * vx(3, -k(2), k(3), k(2) - k(3), mu2, mu3, mu10)
    //             * vx(4, -k(3), k(4), k(3) - k(4), mu3, mu4, mu11)
    //             * vx(5, -k(4), k(0), -k(0) + k(4), mu4, k(20), mu5)
    //             * vx(2, -k(1), k(2), k(1) - k(2), mink(4, mu8), mu2, mu9)
    //             * vx(6, k(0) - k(4), -k(3) + k(4), -k(0) + k(3), mu5, mu11, mu6)
    //             * vx(7, k(0) - k(3), -k(2) + k(3), -k(0) + k(2), mu6, mu10, mu7)
    //             * vx(8, k(0) - k(2), -k(1) + k(2), -k(0) + k(1), mu7, mu9, mu8)
    // );
    for i in 0..8 {
        // r = r.to_dots();
        println!("Hi");
        let mut init = Instant::now();
        r = r
            .replace(parse!(format!(
                "vx({}, k1_, k2_, k3_, mu1_, mu2_, mu3_)",
                i + 1
            )))
            .level_range((0, Some(0)))
            .rhs_cache_size(1000)
            .with_map(rhs_subs.clone());

        println!("RHS time: {:?}", init.elapsed());
        println!("{}", r.nterms());
        init = Instant::now();
        r = r.expand();
        println!("Expansion time: {:?}", init.elapsed());
        println!("{}", r.nterms());
        init = Instant::now();
        // g(mink(4,mu), mink(4,nu)) <-> d(mu,nu)
        // g(k,mink(4),p)<-> mink(4,k,p)
        // mink(4,k+p,l)->mink(4,k,l)+mink(4,p,l)->g(mink(4),k,l)+g(mink(4),p,l)
        // g(mu_,nu_)*f_(..,mu_,..)->f_(..,nu_,..)
        // r = r.simplify_metrics();
        r = r
            .replace(parse!("spenso::g(a_(d_,j_),muw1_(d_,i_))*vx(x___,i_,y___)"))
            .repeat()
            .level_range((0, Some(0)))
            .with(parse!("vx(x___, j_, y___)")); //.simplify_metrics();
        // println!("Metric simplification time: {:?}", init.elapsed());
        println!("{}", r.nterms());
        init = Instant::now();
        r = r
            .replace(parse!("k(a_, muw1_(d_,i_))*vx(x___,i_,y___)"))
            .level_range((0, Some(0)))
            .repeat()
            .with(parse!("vx(x___,k(a_),y___)"));
        println!("Schoonify: {:?}", init.elapsed());
        println!("{}", r.nterms());
        // println!("{:>}", r);
    }
}

#[test]
fn network_informed() {
    run_network_informed();
}
