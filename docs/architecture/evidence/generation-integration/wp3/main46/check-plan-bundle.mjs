#!/usr/bin/env node

import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";

const expected = {
  sourceSha256: "e1ff7645596b88d02de20dc21316c890888746a5768443e4e5d0e33bf506ee85",
  wrapperSha256: "3a9deca513708fd960c6b91612311a4799b7b18ecd5de0425f35b17a17ef3103",
  argvSha256: "ec9ec57c1e20ea65a20b86d85f6c3b0b61de86d1159e18239f9ff6f898ff4700",
  currentHtmlSha256: "626f89df70bd7f65a64a7c7ba8fcfda83f07ab31f7f93eaee70d61866fd2eacd",
  staleHtmlSha256: "fd5f11c55b247eadad8518bace3d9d35b6b823fb253cd59d5cfe65172f309539",
  paragraph: "Source46 ordered integration 2026-10-29: independently approved color-cycle preparation eca80902 / code eb173840 is reconciled onto approved MAIN44–45. Current Idenso owners remain authoritative; only reducing adjoint-cycle identities publish a rewrite, with sixteen networks and 17,484 fresh registration-order permutations retained. Combined evidence is recorded in evidence/generation-integration/wp3/main46/receipt.json. The generic Nix Cargo-target owner now includes immediate child-directory main.rs targets; the corrected Git-pinned Idenso archive realizes and passes 192 profile tests plus 2 separately selected default-filter controls. The inherited TensorReduction declared-layout consumer remains failed and separately owned. A NEW independent Opus A/B review is required before the WP3 rebuilt-tip exit matrix; replay completion alone does not complete WP3.",
};

function sha256(path) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

function decode(text) {
  return text
    .replace(/<[^>]*>/g, "")
    .replaceAll("&lt;", "<")
    .replaceAll("&gt;", ">")
    .replaceAll("&amp;", "&")
    .replace(/&#(\d+);/g, (_, value) => String.fromCodePoint(Number(value)))
    .replace(/&#x([0-9a-f]+);/gi, (_, value) => String.fromCodePoint(Number.parseInt(value, 16)))
    .replace(/\s+/g, " ")
    .trim();
}

function source46Paragraph(path) {
  const html = readFileSync(path, "utf8");
  const paragraphs = [...html.matchAll(/<p>([\s\S]*?)<\/p>/g)].map((match) => decode(match[1]));
  const matches = paragraphs.filter((paragraph) => paragraph.startsWith("Source46 ordered integration 2026-10-29:"));
  if (matches.length !== 1) throw new Error(`expected one source46 paragraph, found ${matches.length}`);
  return matches[0];
}

const [source, wrapper, argv, currentHtml, staleHtml] = process.argv.slice(2);
if (!staleHtml) throw new Error("usage: check-plan-bundle.mjs SOURCE WRAPPER ARGV CURRENT_HTML STALE_HTML");

for (const [label, path, wanted] of [
  ["source", source, expected.sourceSha256],
  ["wrapper", wrapper, expected.wrapperSha256],
  ["argv", argv, expected.argvSha256],
  ["current HTML", currentHtml, expected.currentHtmlSha256],
  ["stale HTML", staleHtml, expected.staleHtmlSha256],
]) {
  const actual = sha256(path);
  if (actual !== wanted) throw new Error(`${label} sha256 mismatch: ${actual}`);
}

const currentParagraph = source46Paragraph(currentHtml);
const staleParagraph = source46Paragraph(staleHtml);
if (currentParagraph !== expected.paragraph) throw new Error("current source46 paragraph is not the frozen expected text");
if (staleParagraph === expected.paragraph) throw new Error("stale HTML unexpectedly matches current text");
if (!staleParagraph.includes("passes192 profile tests plus2")) throw new Error("stale control lacks the known pre-fix text");

console.log(JSON.stringify({
  status: "PASS",
  currentParagraph,
  staleRejected: true,
  hashes: expected,
}, null, 2));
