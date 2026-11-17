= Website font provenance

The approved browser families are original Atkinson Hyperlegible for reading,
Literata for headings, Inter for controls, and IBM Plex Mono for code.
All sources below come from the Google Fonts repository at commit
`5bb28733cec313b63fbd415b574b5ae35c8aa8b8`, retrieved 2026-11-10.
Atkinson comes from `ofl/atkinsonhyperlegible`, never the Next family.

== Selected faces

- Atkinson Hyperlegible: 400 and 700, upright and italic. Intermediate CSS
  weights resolve to the nearest real face; no nonexistent 500/600 face is declared.
- Literata: one upright variable face, weight 200–900, optical size 7–72.
- Inter: one upright variable face, weight 100–900, optical size 14–32.
- IBM Plex Mono: upright 400, 500, 600, 700; italic 400 and 700.

The twelve WOFF2 files total 1093448 bytes. These are complete upstream fonts,
without language or glyph subsetting, preserving punctuation, accents and every
upstream character. The files retain their names, outlines, metrics and variable
axes. The WOFF2 encoding preserves font data; FontTools removes obsolete DSIG
signature tables as required when repackaging. No duplicate TTF or
static version of a variable face is distributed. Literata and Inter italics are
not shipped because their roles use upright headings and controls; prose and
code have genuine regular/bold italics.

`fonts.css` declares exact family names, accurate weight/style descriptors,
asset-relative URLs and `font-display: swap`. The consuming styles supply system
sans, serif or monospace fallback stacks. Unsupported characters use the normal
browser fallback cascade; these text faces do not replace the mathematics font.
There is no external font service request and no `local()` alias that could
select an unrelated installed font.

== License

Each family is distributed under the SIL Open Font License 1.1. The accompanying
`fonts/atkinson-OFL.txt`, `fonts/literata-OFL.txt`, `fonts/inter-OFL.txt` and `fonts/plex-OFL.txt` are the
unaltered upstream license and copyright notices. Retain these files when copying
the font directory. Source hashes for the notices appear below.

== Reproduction and verification

Download the source paths below from the pinned commit. Validate their SHA-256
and Git blob SHA-1 (SHA-1 of the Git blob header followed by the original bytes).
Using FontTools 4.63.0 with Brotli 1.2.0 from the repository-pinned nixpkgs,
load each TTF with `TTFont(path, recalcTimestamp=False)`, assign
`font.flavor = "woff2"`, and call `font.save(output)`.
No instancing, glyph subsetting or hint removal is performed.

Validation reopened every WOFF2 and compared character maps, horizontal metrics,
all name records, OS/2 weight/style metadata, table sets (excluding removed
DSIG signatures), glyph outlines, hint programs and variable-axis values
against the source. ASCII code advances are equal in every IBM Plex Mono face.
SHA-256 values below identify both downloaded source and shipped bytes.

== Files and hashes

=== `AtkinsonHyperlegible-Regular.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/atkinsonhyperlegible/AtkinsonHyperlegible-Regular.ttf")[Pinned upstream source]

Source SHA-256: `7fb917c89019896d0b52ee84b7cbb3304c18cb90b19a62f5e32712bd23e97669`.

Git blob SHA-1: `b3cc4200fcf09d600eac74002ab4ea65b404a536`.

WOFF2 SHA-256: `2df4ba17804bc7a36f123127966075d8427bff2df58d0d76820c1130bb1a4150`.

Size: 23196 bytes; style normal; weight 400; 342 mapped characters.

=== `AtkinsonHyperlegible-Italic.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/atkinsonhyperlegible/AtkinsonHyperlegible-Italic.ttf")[Pinned upstream source]

Source SHA-256: `021beda4d3c6edfc78872e436d74009f9a1bcb294331908fe5747c61d3dcc514`.

Git blob SHA-1: `9ee68183aed5bf84af3bd4eac52e3c723a02d5a5`.

WOFF2 SHA-256: `60bedf0954ffa7e7800eddd529b06c7558d09283667c74fef60170d8d241089e`.

Size: 24908 bytes; style italic; weight 400; 342 mapped characters.

=== `AtkinsonHyperlegible-Bold.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/atkinsonhyperlegible/AtkinsonHyperlegible-Bold.ttf")[Pinned upstream source]

Source SHA-256: `5a3b0c8cc8ca545155150b4512a1fa248298df121c50d6557e651e61fbdab92f`.

Git blob SHA-1: `ca3f613847de62b8da83d2dab500ec391c8fe477`.

WOFF2 SHA-256: `da8fce41a04f8498fbf79076f92d304b12e70c76f71b143c5dcfb6536c93c075`.

Size: 23776 bytes; style normal; weight 700; 342 mapped characters.

=== `AtkinsonHyperlegible-BoldItalic.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/atkinsonhyperlegible/AtkinsonHyperlegible-BoldItalic.ttf")[Pinned upstream source]

Source SHA-256: `7a46f1b0ef2f64b61065c1a6b5a11fcb9b019b44062d88ba121a083cc39e86bf`.

Git blob SHA-1: `2804fdb8f06c328afb96409d7498c6f5f287b79c`.

WOFF2 SHA-256: `f27e95143bc12e8f40f955cf86c74e6afef4441abb89d271282ea0f9e93c1fd5`.

Size: 25160 bytes; style italic; weight 700; 342 mapped characters.

=== `atkinson-OFL.txt`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/atkinsonhyperlegible/OFL.txt")[Pinned upstream source]

Source SHA-256: `f32d22b3908fcad2c86a74000614ec22e6a7f66ea7e867e616026a27aebdc143`.

Git blob SHA-1: `2befc89b0ee81d5a37e6ca6503bb0db099feb70f`.

=== `Literata-Variable.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/literata/Literata%5Bopsz%2Cwght%5D.ttf")[Pinned upstream source]

Source SHA-256: `b41138c9373112f32abb589cc22e8674b06ed4048b0c513be922bdd26f274440`.

Git blob SHA-1: `0a81b4984a4a4eb101e0ddc28babc85bb137efef`.

WOFF2 SHA-256: `230e5b1de4f4ace4d49c9fda9c8cf0851acb5425dce9067d89ff88aed93f325d`.

Size: 397212 bytes; style normal; weight 200 900; 1163 mapped characters.

=== `literata-OFL.txt`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/literata/OFL.txt")[Pinned upstream source]

Source SHA-256: `8742963604cd89dc81437811a850018fc03b2bfad686d7422c8235967c87614e`.

Git blob SHA-1: `9d3b4e3b1e9068fbda799f57833ff831d05cc271`.

=== `Inter-Variable.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/inter/Inter%5Bopsz%2Cwght%5D.ttf")[Pinned upstream source]

Source SHA-256: `29160a80ff49ddcab2c97711247e08b1fab27a484a329ce8b813d820dc559031`.

Git blob SHA-1: `047c92f6e2212473dc436020afed689527076d44`.

WOFF2 SHA-256: `13b4cddc57045d40b411008e61588c8f81aa0bd72fe884624831a0ff9c4d159a`.

Size: 350460 bytes; style normal; weight 100 900; 2849 mapped characters.

=== `inter-OFL.txt`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/inter/OFL.txt")[Pinned upstream source]

Source SHA-256: `5b9321a4298cfeb6b34354164a1c3afc3db114569984c502b9b35d988fd58c57`.

Git blob SHA-1: `21f6aff961064c2e429f570995e446bcdd555422`.

=== `IBMPlexMono-Regular.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/ibmplexmono/IBMPlexMono-Regular.ttf")[Pinned upstream source]

Source SHA-256: `6a3412f058c7d8dfd9170c41e85ade48e5156ecb89356110ca57a0a27734af46`.

Git blob SHA-1: `0c9770d5183ba60dc4350d3e011b782a320761ae`.

WOFF2 SHA-256: `13b16cd770a73f60304510221c3e7186d396ba6ca5d64c52b5c6fb983d0e900a`.

Size: 39176 bytes; style normal; weight 400; 930 mapped characters.

=== `IBMPlexMono-Medium.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/ibmplexmono/IBMPlexMono-Medium.ttf")[Pinned upstream source]

Source SHA-256: `a9b4c49bb299e05b5f6c481e7fb5e78943d2793249a0c8874ab574a2d1ea6755`.

Git blob SHA-1: `33c546f68b6007a45a3f10e845523abb2db25399`.

WOFF2 SHA-256: `0498bcb48294756b5ac9062290e1639e907b1ec69b5dc2e8cfdec3fab58e0bc8`.

Size: 40080 bytes; style normal; weight 500; 930 mapped characters.

=== `IBMPlexMono-SemiBold.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/ibmplexmono/IBMPlexMono-SemiBold.ttf")[Pinned upstream source]

Source SHA-256: `d3c38e55c78f5b0f28009fddba4834ec503278936a5986032424c9bd2d23aa46`.

Git blob SHA-1: `8f8998a139d8f61d3b96bdd3b512fca501a1f6b9`.

WOFF2 SHA-256: `5148e1fccb46851e8a7dcc489ff7d3d5e54717ee72e2865c2fa225e3fc4744e8`.

Size: 41136 bytes; style normal; weight 600; 930 mapped characters.

=== `IBMPlexMono-Bold.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/ibmplexmono/IBMPlexMono-Bold.ttf")[Pinned upstream source]

Source SHA-256: `ac27abd6450a64dd94467580a02fe6235156d5b92f2926ebbc8e7489df64e0be`.

Git blob SHA-1: `17b406d221e7f27179344707bd543fb8470728ac`.

WOFF2 SHA-256: `9f9e5597981645d576b209a25021c5221e8644a4580581f44e8b4ef954d95d03`.

Size: 40348 bytes; style normal; weight 700; 930 mapped characters.

=== `IBMPlexMono-Italic.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/ibmplexmono/IBMPlexMono-Italic.ttf")[Pinned upstream source]

Source SHA-256: `3362fc791b0652193328b862c1c5f23a789bc7288b1617fa63302f88689a2a34`.

Git blob SHA-1: `d7f69134614ebdcaa10ebd36b0be1a44a035e157`.

WOFF2 SHA-256: `ee9d5440c380188b64e98ebbfe6b6aa45f3299573e39f86f97456290320772b8`.

Size: 43716 bytes; style italic; weight 400; 930 mapped characters.

=== `IBMPlexMono-BoldItalic.woff2`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/ibmplexmono/IBMPlexMono-BoldItalic.ttf")[Pinned upstream source]

Source SHA-256: `af4e05a761e98c1adf064c48a6352c9bec1a6ad70982cd2a544149323391f98e`.

Git blob SHA-1: `f0291778c969fca0c255fe7e06c90eb97b34770d`.

WOFF2 SHA-256: `ea6bad6fa3077ee1a65ba470d5c7bda931a52f91f1f5a26b800ccd5401462d3d`.

Size: 44280 bytes; style italic; weight 700; 930 mapped characters.

=== `plex-OFL.txt`

#link("https://raw.githubusercontent.com/google/fonts/5bb28733cec313b63fbd415b574b5ae35c8aa8b8/ofl/ibmplexmono/OFL.txt")[Pinned upstream source]

Source SHA-256: `7e6b2818edbd8f6a01ae80641cc8f16a51080d08fb4e532be3a0b6f74adb07da`.

Git blob SHA-1: `c35c4c618fab33da8695177b3a6cefe0810b7b28`.

