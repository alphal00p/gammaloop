# Documentation image provenance

The portraits below are size-optimized WebP derivatives of public professional profile images. Every SVG used by the website is generated from editable Typst under `docs/assets/typst/`, with graph layout shared from GammaLoop's `assets/embedded/drawing/templates/layout-core.typ`; run `nix develop --command just docs-svg-assets` after changing them.

`STIXTwoMath-Regular.woff2` is the unmodified STIX Two Math 2.13 webfont from the [official STIX Fonts release](https://github.com/stipub/stixfonts/releases/tag/v2.13). It is distributed under the SIL Open Font License 1.1 in `STIX-Two-OFL.txt` and is bundled so MathML renders consistently without relying on a system font.

| Local asset | Subject | Public source |
| --- | --- | --- |
| `people/valentin.webp` | Valentin Hirschi | [University of Bern profile](https://www.itp.unibe.ch/about_us/people/people/index_eng.html?id=248) · [original portrait](https://itpcenter.itp.unibe.ch/ajax/proxy.php/people/Hirschi_Valentin.jpg) |
| `people/mathijs.webp` | Mathijs Fraaije | [University of Bern profile](https://www.itp.unibe.ch/about_us/people/people/index_eng.html?id=250) · [original portrait](https://itpcenter.itp.unibe.ch/ajax/proxy.php/people/Fraaije_Mathijs.jpg) |
| `people/lucien.webp` | Lucien Huber | [University of Bern profile](https://www.itp.unibe.ch/about_us/people/people/index_eng.html?id=260) · [original portrait](https://itpcenter.itp.unibe.ch/ajax/proxy.php/people/Huber_Lucien.jpg) |
| `people/kaapo.webp` | Kaapo Seppänen | [University of Bern profile](https://www.itp.unibe.ch/about_us/people/people/index_eng.html?id=290) · [original portrait](https://itpcenter.itp.unibe.ch/ajax/proxy.php/people/Sepp%C3%A4nen_Kaapo.jpg) |
| `people/zeno.webp` | Zeno Capatti | [University of Bern profile](https://www.itp.unibe.ch/about_us/people/people/index_eng.html?id=262) · [original portrait](https://itpcenter.itp.unibe.ch/ajax/proxy.php/people/Capatti_Zeno.jpg) |
| `people/ben.webp` | Ben Ruijl | [Symbolica profile](https://symbolica.io/about.html) · [original portrait](https://symbolica.io/ben.jpg) |
| `spensologo.svg` | Spenso | Native Typst reconstruction of the [canonical logo at the pinned upstream revision](https://github.com/alphal00p/spenso/blob/c052f22dc98a18535e114eede758674befc2758f/spensologo.svg) |

The graph sources read the real generated-process and test-resource DOT files directly and pass them to Linnest from Typst; no Graphviz, `just draw`, or generated template bundle is part of the website build. Four sources enable Linnest momentum arrows explicitly, and the shared edge style owns the website-specific line weights and light/dark palette.

The portraits and referenced project design are public professional/project material. No explicit public redistribution license was found for the portrait sources; preserve this provenance and obtain permission where required before reusing them elsewhere.
