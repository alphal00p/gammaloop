//! Native drawings and the written forest formula of woods, forests and
//! counterterms.
#![cfg(feature = "svg")]

use std::collections::BTreeSet;

use feynkit_graph::CountertermScheme as Glyph;
use feynkit_graph::{DiagramError, FeynmanDiagram, HighlightStyle, LoopMomentumBasis};
use feynkit_model::Model;
use feynkit_uv::{
    ApproximationType, Counterterm, CountertermScheme, FULL_FILL, ROOT_FILL, RenameError,
    RenormalizationPrescriptionSettings, StampOptions, UvDisplayError, UvGraph, Wood,
};
use linnest::svg::Config;

fn wood_with(
    model: Model,
    dot: &str,
    prescription: &RenormalizationPrescriptionSettings,
) -> (FeynmanDiagram, Wood<LoopMomentumBasis>) {
    let diagram = FeynmanDiagram::from_dot(model, dot).unwrap();
    let spinneys = diagram
        .classify_spinneys(
            &diagram.hedge_graph().full_filter(),
            diagram.model(),
            prescription,
            diagram.loop_momentum_basis(),
        )
        .unwrap();
    let wood = Wood::from_spinneys(spinneys, &diagram, diagram.loop_momentum_basis()).unwrap();
    (diagram, wood)
}

fn build(model: Model, dot: &str) -> (FeynmanDiagram, Wood<LoopMomentumBasis>) {
    wood_with(model, dot, &RenormalizationPrescriptionSettings::default())
}

const BASKETBALL: &str =
    "digraph basketball { edge [particle=phi]; a -> b; a -> b; a -> b; a -> b; }";
/// The vacuum sunrise of the reference BPHZ notes: three overlapping bubbles
/// a, b, c in the whole graph Γ.
const SUNRISE: &str = "digraph sunrise { edge [particle=phi]; a -> b; a -> b; a -> b; }";
/// A cubic stand-in for the dotted sunrise: a bubble a nested in the whole graph
/// Γ, with the dot's two-point vertex replaced by an external leg.
const DOTTED: &str = "digraph dotted { edge [particle=phi]; ext [style=invis]; \
     a -> b; a -> b; a -> c; c -> b; ext -> c; }";

#[test]
fn the_sunrise_expansion_follows_the_written_forest_formula() {
    let (sunrise, wood) = build(Model::phi3(), SUNRISE);
    let forests = wood.unfold();
    let expansion = forests.expansion(&sunrise);
    // m(Γ|∅) − m(Γ|Γ)·K[m(Γ|∅)] − m(Γ|a)·K[m(a|∅)] − … + m(Γ|Γ)·K[m(Γ|a)·K[m(a|∅)]] + …
    let written: Vec<_> = expansion
        .terms()
        .iter()
        .map(|term| (term.k, term.sign, expansion.notation(term, "K")))
        .collect();
    assert_eq!(
        written,
        [
            (0, 1, "K[∅]".to_owned()),
            (1, -1, "K_{Γ}".to_owned()),
            (2, -1, "K_{a}".to_owned()),
            (3, -1, "K_{b}".to_owned()),
            (4, -1, "K_{c}".to_owned()),
            (5, 1, "K_{Γ|a}".to_owned()),
            (6, 1, "K_{Γ|b}".to_owned()),
            (7, 1, "K_{Γ|c}".to_owned()),
        ]
    );
    // The nested term of K_{Γ|a} is K_{a}, and its cofactor treats Γ.
    let nested = &expansion.terms()[5];
    assert_eq!(nested.outermost.len(), 1);
    assert_eq!(expansion.name(nested.outermost[0]), Some("Γ"));
    let inner = &nested.factors[0].inner;
    assert_eq!(inner.len(), 1);
    assert!(inner[0].inner.is_empty() && inner[0].given.is_empty());
    assert_eq!(expansion.name(inner[0].spinney), Some("a"));

    let svg = forests
        .expansion_svg(
            &sunrise,
            &expansion,
            &Config::default(),
            StampOptions::default(),
            "K",
        )
        .unwrap();
    assert_eq!(svg.matches("data-feynkit-term=").count(), 8);
    // Cofactor and factor stamps: 1 + 2·4 + 3·3 uses for the eight terms.
    assert_eq!(svg.matches("<use ").count(), 1 + 2 * 4 + 3 * 3);
    assert_eq!(svg.matches(">−<").count(), 4);
    assert!(svg.contains("<title>K_{Γ|a} = k5 · depth 2, sign +</title>"));
    // Stamps are cropped to their ink: no symbol keeps Linnest's blank border.
    let symbols: Vec<_> = svg.split("<symbol ").skip(1).collect();
    assert!(!symbols.is_empty());
    for symbol in symbols {
        let header = &symbol[..symbol.find('>').unwrap()];
        assert!(!header.contains("viewBox=\"0 0 "), "{header}");
    }
}

#[test]
fn spinneys_take_user_names_and_reject_unknown_or_duplicate_ones() {
    let (sunrise, wood) = build(Model::phi3(), SUNRISE);
    let forests = wood.unfold();
    let mut expansion = forests.expansion(&sunrise);
    let node = |k: usize| expansion.terms()[k].factors[0].spinney;
    let (whole, a, b) = (node(1), node(2), node(3));
    expansion
        .rename([(whole, "G".to_owned()), (a, "a1".to_owned())])
        .unwrap();
    assert_eq!(expansion.notation(&expansion.terms()[5], "K"), "K_{G|a1}");
    assert_eq!(
        expansion.rename([(b, "a1".to_owned())]),
        Err(RenameError::Duplicate("a1".to_owned()))
    );
    assert_eq!(
        expansion.rename([(b, "x|y".to_owned())]),
        Err(RenameError::InvalidName("x|y".to_owned()))
    );
    // Forest nodes are no spinneys of the wood.
    let forest = expansion.terms()[7].forest;
    assert_eq!(
        expansion.rename([(forest, "z".to_owned())]),
        Err(RenameError::UnknownSpinney(forest))
    );
    // A failed rename leaves the names unchanged.
    assert_eq!(expansion.name(b), Some("b"));
    // Names that read as notation, Γ for another spinney, and two names for one
    // spinney are rejected.
    for name in ["∅", "x·y", "k1", "k12", "a b", "x\u{1}y", "[x]"] {
        assert_eq!(
            expansion.rename([(b, name.to_owned())]),
            Err(RenameError::InvalidName(name.to_owned())),
            "{name:?}"
        );
    }
    assert_eq!(
        expansion.rename([(b, "Γ".to_owned())]),
        Err(RenameError::ReservedName)
    );
    assert_eq!(
        expansion.rename([(b, "x".to_owned()), (b, "y".to_owned())]),
        Err(RenameError::RenamedTwice(b))
    );
    // "k" alone and a k with letters are ordinary names.
    expansion
        .rename([(b, "k".to_owned()), (a, "kx".to_owned())])
        .unwrap();
}

#[test]
fn expansions_group_wrapped_lines_by_depth() {
    let (basketball, wood) = build(Model::phi4(), BASKETBALL);
    let forests = wood.unfold();
    let expansion = forests.expansion(&basketball);
    let svg = forests
        .expansion_svg(
            &basketball,
            &expansion,
            &Config::default(),
            StampOptions::default(),
            "K",
        )
        .unwrap();
    // Each line names its depth; depths 0 to 3 appear in order, and only the
    // first line of each depth is not a continuation.
    let lines: Vec<(usize, bool)> = svg
        .split("<g data-feynkit-depth=\"")
        .skip(1)
        .map(|rest| {
            let header = &rest[..rest.find('>').unwrap()];
            let depth = header[..header.find('"').unwrap()].parse().unwrap();
            (depth, header.contains("data-feynkit-continued=\"true\""))
        })
        .collect();
    let depths: Vec<usize> = lines.iter().map(|(depth, _)| *depth).collect();
    assert!(depths.windows(2).all(|pair| pair[0] <= pair[1]));
    assert_eq!(depths.iter().copied().collect::<BTreeSet<_>>().len(), 4);
    for (index, (depth, continued)) in lines.iter().enumerate() {
        let starts = index == 0 || lines[index - 1].0 != *depth;
        assert_eq!(*continued, !starts, "line {index} of depth {depth}");
    }
    // Wide depths wrap, so there are more lines than depths.
    assert!(lines.len() > 4);
    assert_eq!(svg.matches("data-feynkit-term=").count(), 46);
}

#[test]
fn branching_counterterm_hover_groups_siblings_under_their_parent() {
    // Three a–b lines (0–2), two c–d lines (3–4), connectors a–c (5) and b–d (6).
    let (graph, wood) = build(
        Model::phi4(),
        "digraph branch { edge [particle=phi]; ext [style=invis]; \
             a -> b [id=0]; a -> b [id=1]; a -> b [id=2]; c -> d [id=3]; c -> d [id=4]; \
             a -> c [id=5]; b -> d [id=6]; ext -> c; d -> ext; }",
    );
    let pick = |label: &str| {
        wood.spinneys()
            .map(|(_, spinney)| spinney)
            .find(|spinney| spinney.label(&graph) == label)
            .unwrap_or_else(|| panic!("{label} is a spinney"))
            .clone()
    };
    let counterterm = Counterterm::new(
        [
            pick("γ{0,1}"),
            pick("γ{0,1,2}"),
            pick("γ{3,4}"),
            pick("γ{0,1,2,3,4,5,6}"),
        ],
        CountertermScheme::MsBar,
        &graph,
        graph.loop_momentum_basis(),
    )
    .unwrap();
    let regions = counterterm.regions(&graph);
    let crossed: Vec<_> = regions
        .iter()
        .filter(|region| region.details.contracted)
        .collect();
    assert_eq!(crossed.len(), 1);
    assert_eq!(
        crossed[0].details.label,
        "(γ{0,1} ⊂ γ{0,1,2}, γ{3,4}) ⊂ γ{0,1,2,3,4,5,6}"
    );
}

#[test]
fn the_scene_api_rejects_a_bad_stamp_weight() {
    let (sunrise, wood) = build(Model::phi3(), SUNRISE);
    let (_, spinney) = wood
        .spinneys()
        .find(|(_, spinney)| !spinney.is_empty())
        .unwrap();
    let regions = [spinney.region(&sunrise)];
    for weight in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let options = feynkit_graph::SceneOptions {
            overlap: feynkit_graph::HighlightOverlap::Recolour,
            stamp_weight: weight,
            ..Default::default()
        };
        assert!(matches!(
            sunrise.to_scene_with_highlights(
                &feynkit_uv::UvRegion::layers(&regions),
                None,
                &options
            ),
            Err(DiagramError::SceneConfiguration(_))
        ));
    }
}

#[test]
fn the_dotted_sunrise_numbers_its_stored_quantities_topologically() {
    let (dotted, wood) = build(Model::phi3(), DOTTED);
    let forests = wood.unfold();
    let expansion = forests.expansion(&dotted);
    // As in the expected table: K[∅] = k0, K(Γ) = k1, K(a) = k2, K(Γ,a,k2) = k3.
    assert_eq!(
        expansion.k_table("K"),
        ["k0 = K[∅]", "k1 = K_{Γ}", "k2 = K_{a}", "k3 = K_{Γ|a}[k2]"]
    );
    let signs: Vec<_> = expansion.terms().iter().map(|term| term.sign).collect();
    assert_eq!(signs, [1, -1, -1, 1]);
    assert_eq!(expansion.k_table("T^{3D}")[3], "k3 = T^{3D}_{Γ|a}[k2]");
}

#[test]
fn products_of_disjoint_factors_multiply_their_stored_quantities() {
    let (spectacles, wood) = build(
        Model::phi3(),
        "digraph spectacles { edge [particle=phi]; a -> a; a -> b; b -> b; }",
    );
    let forests = wood.unfold();
    let expansion = forests.expansion(&spectacles);
    assert_eq!(
        expansion.k_table("K"),
        ["k0 = K[∅]", "k1 = K_{a}", "k2 = K_{b}", "k3 = k1·k2"]
    );
    let product = &expansion.terms()[3];
    assert_eq!(expansion.notation(product, "K"), "K_{a}K_{b}");
    assert_eq!(product.outermost.len(), 2);
}

#[test]
fn operation_graphs_write_each_forest_as_its_term() {
    let (basketball, wood) = build(Model::phi4(), BASKETBALL);
    let forests = wood.unfold();
    let svg = forests
        .operation_graph_svg(
            &basketball,
            &forests.expansion(&basketball),
            &Config::default(),
            StampOptions::default(),
            "K",
        )
        .unwrap();
    assert_eq!(svg.matches("data-feynkit-forest=").count(), 46);
    assert_eq!(svg.matches("data-feynkit-k=").count(), 46);
    // Each distinct stamp is defined once.
    let ids: BTreeSet<_> = svg
        .split("<symbol id=\"")
        .skip(1)
        .map(|rest| rest.split('"').next().unwrap())
        .collect();
    assert_eq!(ids.len(), svg.matches("<symbol ").count());
    // One stamp per factor, no cofactors: one use per operation of every forest.
    let operations: usize = forests.operations().map(|(_, forest)| forest.depth()).sum();
    assert_eq!(svg.matches("<use ").count(), operations);
    assert!(svg.contains(">𝟙<"));
    assert!(!svg.contains("<a "));
    assert_eq!(svg.matches(ROOT_FILL).count(), 1);
    assert!(svg.contains(FULL_FILL));
    // Any operator symbol; scripts are set as superscripts.
    let local = forests
        .operation_graph_svg(
            &basketball,
            &forests.expansion(&basketball),
            &Config::default(),
            StampOptions::default(),
            "T^{3D}",
        )
        .unwrap();
    assert!(local.contains(r#"baseline-shift="0.45em">3D<"#));
}

#[test]
fn stamp_weight_scales_every_stamp_line() {
    let (basketball, wood) = build(Model::phi4(), BASKETBALL);
    let widths = |weight: f64| {
        let svg = wood
            .hasse_svg(
                &basketball,
                &Config::default(),
                StampOptions {
                    weight,
                    ..StampOptions::default()
                },
            )
            .unwrap();
        svg.split("stroke-width=\"")
            .skip(1)
            .filter_map(|rest| rest.split('"').next()?.parse::<f64>().ok())
            .fold(0.0, f64::max)
    };
    let (light, heavy) = (widths(1.0), widths(2.0));
    assert!(heavy > 1.9 * light, "{light} {heavy}");
}

#[test]
fn hasse_diagrams_rank_spinneys_and_reject_bad_stamp_options() {
    let (basketball, wood) = build(Model::phi4(), BASKETBALL);
    let svg = wood
        .hasse_svg(&basketball, &Config::default(), StampOptions::default())
        .unwrap();
    assert_eq!(svg.matches("data-feynkit-spinney=").count(), 12);
    assert!(svg.contains(
        "<title>γ{0,1,2,3} · 3 loops, degree 4 (superficial 4), MUV approximation</title>"
    ));
    for bad in [0.0, -1.0, f64::NAN] {
        let size = StampOptions {
            size: bad,
            ..StampOptions::default()
        };
        assert!(matches!(
            wood.hasse_svg(&basketball, &Config::default(), size),
            Err(UvDisplayError::StampSize(_))
        ));
        let weight = StampOptions {
            weight: bad,
            ..StampOptions::default()
        };
        assert!(matches!(
            wood.hasse_svg(&basketball, &Config::default(), weight),
            Err(UvDisplayError::StampWeight(_))
        ));
    }
}

#[test]
fn counterterms_cross_each_outermost_spinney_once() {
    let (basketball, wood) = build(Model::phi4(), BASKETBALL);
    let forests = wood.unfold();
    let (node, _) = forests
        .operations()
        .find(|(_, forest)| forest.depth() == 3)
        .unwrap();
    let regions = forests.regions(node, &basketball);
    let styles: Vec<_> = regions.iter().map(|region| region.style).collect();
    assert_eq!(
        styles,
        [
            HighlightStyle::ActedOn,
            HighlightStyle::Middle,
            HighlightStyle::Treated
        ]
    );
    let counterterm =
        Counterterm::of_forest::<DiagramError>(&forests, node, CountertermScheme::OnShell).unwrap();
    let contracted = counterterm.regions(&basketball);
    // One crossed region, the outermost; the chain is listed inner to outer.
    let crossed: Vec<_> = contracted
        .iter()
        .filter(|region| region.details.contracted)
        .collect();
    assert_eq!(crossed.len(), 1);
    assert_eq!(crossed[0].style, HighlightStyle::Treated);
    assert_eq!(crossed[0].details.scheme, Glyph::OnShell);
    assert_eq!(crossed[0].details.label.matches(" ⊂ ").count(), 2);
    assert!(crossed[0].details.label.ends_with("γ{0,1,2,3}"));
    // The nested ones are hulls coloured by depth.
    let nested: Vec<_> = contracted[1..].iter().map(|region| region.style).collect();
    assert_eq!(nested, [HighlightStyle::Middle, HighlightStyle::ActedOn]);

    // Vertex-disjoint outermost spinneys each get their own cross.
    let (spectacles, wood) = build(
        Model::phi3(),
        "digraph spectacles { edge [particle=phi]; a -> a; a -> b; b -> b; }",
    );
    let forests = wood.unfold();
    let (product, _) = forests
        .operations()
        .find(|(node, _)| forests.is_disjoint_union(*node))
        .unwrap();
    let counterterm =
        Counterterm::of_forest::<DiagramError>(&forests, product, CountertermScheme::MsBar)
            .unwrap();
    let regions = counterterm.regions(&spectacles);
    assert_eq!(regions.len(), 2);
    assert!(regions.iter().all(|region| region.details.contracted));
}

#[test]
fn the_default_scheme_follows_on_shell_approximations_only() {
    // Subtract the tadpoles on shell; the default prescription keeps MUV.
    let dot = "digraph spectacles { edge [particle=phi]; a -> a; a -> b; b -> b; }";
    let (_, muv) = build(Model::phi3(), dot);
    let tadpoles: Vec<_> = muv
        .spinneys()
        .map(|(_, spinney)| spinney)
        .filter(|spinney| spinney.n_components() == 1 && !spinney.is_empty())
        .collect();
    assert_eq!(
        CountertermScheme::implied_by(tadpoles.iter().copied()),
        CountertermScheme::Unspecified
    );
    let prescription = RenormalizationPrescriptionSettings {
        log_divergent: ApproximationType::OS,
        massive_power_divergent: ApproximationType::OS,
        massless_power_divergent: ApproximationType::OS,
        ..RenormalizationPrescriptionSettings::default()
    };
    let (_, on_shell) = wood_with(Model::phi3(), dot, &prescription);
    let tadpoles: Vec<_> = on_shell
        .spinneys()
        .map(|(_, spinney)| spinney)
        .filter(|spinney| spinney.n_components() == 1 && !spinney.is_empty())
        .collect();
    assert!(
        tadpoles
            .iter()
            .all(|spinney| spinney.renormalization_scheme() == ApproximationType::OS),
        "{:?}",
        tadpoles
            .iter()
            .map(|spinney| spinney.renormalization_scheme())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        CountertermScheme::implied_by(tadpoles.iter().copied()),
        CountertermScheme::OnShell
    );
}

#[test]
fn product_forests_draw_and_contract_each_factor() {
    let (pair, wood) = build(
        Model::phi4(),
        "digraph pair { edge [particle=phi]; ext [style=invis]; \
         a -> b [id=0]; a -> b [id=1]; b -> c [id=2]; c -> d [id=3]; c -> d [id=4]; \
         ext -> a; ext -> a; b -> ext; c -> ext; d -> ext; d -> ext; }",
    );
    let forests = wood.unfold();
    let (product, _) = forests
        .operations()
        .find(|(node, _)| forests.is_disjoint_union(*node))
        .unwrap();
    // Both bubbles are outermost: one paint, one region each, no union region.
    let regions = forests.regions(product, &pair);
    let mut labels: Vec<_> = regions.iter().map(|r| r.details.label.as_str()).collect();
    labels.sort_unstable();
    assert_eq!(labels, ["γ{0,1}", "γ{3,4}"]);
    assert!(regions.iter().all(|r| r.style == HighlightStyle::ActedOn));
    let counterterm =
        Counterterm::of_forest::<DiagramError>(&forests, product, CountertermScheme::MsBar)
            .unwrap();
    let mut contracted: Vec<_> = counterterm
        .contracted()
        .iter()
        .map(|s| s.label(&pair))
        .collect();
    contracted.sort_unstable();
    assert_eq!(contracted, ["γ{0,1}", "γ{3,4}"]);
    // The operation graph multiplies two single-bubble factors.
    let svg = forests
        .operation_graph_svg(
            &pair,
            &forests.expansion(&pair),
            &Config::default(),
            StampOptions::default(),
            "K",
        )
        .unwrap();
    assert!(svg.contains("<title>K_{a}K_{b} = k3 · depth 2, sign +</title>"));
}
