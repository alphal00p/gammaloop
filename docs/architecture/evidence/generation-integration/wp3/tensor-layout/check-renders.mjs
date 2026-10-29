#!/usr/bin/env node

import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";

const expected = [
  {
    sourceSha256: "566a90a1648a4e570baa562cb73ece8dff5d7b1001f5ccd4e4d6f80e4c3269d0",
    htmlSha256: "b376397502c87e3b25e2e2db7c772f59f73f6357f981aed4ca6c598abd5d20fc",
    paragraph: "Fresh combined validation records 203 affected tensor/HEP passes with 3 ignored, 17 tensor doctests, exact in-tree and held-out materializer probes, 16 Idenso owner controls, 192 color controls plus 2 default-filter controls, and bounded Python, graph and Vakint consumer checks. Source44 physical captures, source46 frame/color proofs and the realized source46 Idenso archive are reused only after exact source and derivation invariance checks.",
  },
  {
    sourceSha256: "54cef3a7c7911d71a89133c71e3f7db6ca71e933273307d62c68d0245527d669",
    htmlSha256: "1c5132ad651f1fde0ba4f9288f070808ef71e9333fe63e925258637416fb71fc",
    paragraph: "TensorReduction declared-layout integration 2026-10-29: approved preparation 97269ac1 is reconciled byte-for-byte onto approved MAIN46. Checked positional sums retain their established correspondence; other sums are restored through the reducer-owned source carrier, while generic Idenso arithmetic remains canonical. Fresh bounded consumer, held-out, color, Python and Vakint validation is recorded in evidence/generation-integration/wp3/tensor-layout/receipt.json. A NEW independent Opus integration review is required before the rebuilt-tip WP3 exit matrix; this integration does not complete WP3.",
  },
  {
    sourceSha256: "cf321633ad121ba884711329d159cba42a54dccf72f0d836dce67ac8a3837c7e",
    htmlSha256: "04e790d43b2bc267a61fdb15395fa0e1e916139dcf63a72c0f1397f3ae0834a2",
    paragraph: "The raw Vakint adapter uses TensorReducer::admit to admit each converted selected numerator once as a SymbolicTensor<PartialStructure>, then calls TensorReducer::reduce with that retained value. Rejected source expressions keep the reducer’s precise index/power error when its existing parser can refine the native admission error; validation is never bypassed. TensorReductionTerm keeps its tensor body in the same carrier; its separate Atom coefficient is the exact scalar orbit or projector weight. Foreign-representation spectators remain tensors during selected Lorentz reduction and composition. TensorReduction also retains the source interface for an empty or identically vanishing result, including odd-rank moments. Its borrowed expression() and consuming into_expression() use generic typed addition for compatibility and algebra. Positional sums retain the established logical correspondence owned by that checked addition. Named sums, whose generic arithmetic intentionally publishes canonical axes, are validated against the reducer-owned source carrier to restore its declared logical layout. This is a caller-owned boundary; arbitrary positional rewrites still require their own identity witness. Only the adapter’s return to Vakint’s raw projector/index ABI decodes the completed result; raw VakintTerm storage remains that backend boundary.",
  },
];

function sha256(path) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

function decode(text) {
  return text.replace(/<[^>]*>/g, "")
    .replaceAll("&lt;", "<").replaceAll("&gt;", ">").replaceAll("&amp;", "&")
    .replace(/&#(\d+);/g, (_, value) => String.fromCodePoint(Number(value)))
    .replace(/&#x([0-9a-f]+);/gi, (_, value) => String.fromCodePoint(Number.parseInt(value, 16)))
    .replace(/\s+/g, " ").trim();
}

function paragraphs(path) {
  const html = readFileSync(path, "utf8");
  return [...html.matchAll(/<p>([\s\S]*?)<\/p>/g)].map((match) => decode(match[1]));
}

const args = process.argv.slice(2);
if (args.length !== 6) throw new Error("usage: check-renders.mjs SOURCE HTML SOURCE HTML SOURCE HTML");

for (let index = 0; index < expected.length; index += 1) {
  const source = args[index * 2];
  const html = args[index * 2 + 1];
  if (sha256(source) !== expected[index].sourceSha256) throw new Error(`source ${index} hash mismatch`);
  if (sha256(html) !== expected[index].htmlSha256) throw new Error(`HTML ${index} hash mismatch`);
  if (!paragraphs(html).includes(expected[index].paragraph)) throw new Error(`HTML ${index} lacks exact frozen paragraph`);
}

console.log(JSON.stringify({ status: "PASS", expected }, null, 2));
