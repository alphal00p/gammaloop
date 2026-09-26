#!/usr/bin/env node
// Render the front-page GammaLoop tour (docs/assets/showcase.js|css) to a video, or write a
// standalone preview page. The stage is a pure function of time, so frames are captured
// deterministically with Playwright and piped to ffmpeg together with the soundtrack
// (docs/assets/showcase-hard-boiled.mp3).
//
//   node scripts/render-showcase-video.mjs                 target/showcase/gammaloop-showcase.mp4
//   node scripts/render-showcase-video.mjs --preview       target/showcase/index.html (needs fonts online)
//   node scripts/render-showcase-video.mjs --standalone    one self-contained HTML page with embedded media
//
// Options: --out DIR, --fps 30, --from 0, --to END, --width 1920, --height 1080, --poster 5.6,
//          --assets DIR, --ffmpeg PATH, --no-fonts, --no-music.
// Requirements: Node 20+, the `playwright` package with a Chromium download, ffmpeg with
// libx264 and aac, and curl for the one-time Google Fonts download (cached under --out).
import { spawn, spawnSync } from 'node:child_process';
import { once } from 'node:events';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const args = process.argv.slice(2);
const opt = (name, fallback) => {
  const i = args.indexOf(`--${name}`);
  return i >= 0 && i + 1 < args.length ? args[i + 1] : fallback;
};
const flag = (name) => args.includes(`--${name}`);

const assets = path.resolve(opt('assets', path.join(here, '..', 'docs', 'assets')));
const out = path.resolve(opt('out', path.join(here, '..', 'target', 'showcase')));
const fps = Number(opt('fps', 30));
const width = Number(opt('width', 1920));
const height = Number(opt('height', 1080));
fs.mkdirSync(out, { recursive: true });

const css = fs.readFileSync(path.join(assets, 'showcase.css'), 'utf8');
const js = fs.readFileSync(path.join(assets, 'showcase.js'), 'utf8');
const soundtrackName = 'showcase-hard-boiled.mp3';
const soundtrack = path.join(assets, soundtrackName);

const FAMILIES =
  'family=IBM+Plex+Mono:wght@400;500;600&family=Inter:wght@400;500;600;700&family=Newsreader:ital,wght@0,400;0,600;1,400&display=swap';
const GOOGLE_FONTS_CSS = `https://fonts.googleapis.com/css2?${FAMILIES}`;
const USER_AGENT = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/128.0 Safari/537.36';

function curl(url, destination) {
  const result = spawnSync('curl', ['-sS', '-L', '-m', '60', '-A', USER_AGENT, '-o', destination, url], {
    stdio: ['ignore', 'ignore', 'pipe'],
  });
  if (result.status !== 0) throw new Error(`curl failed for ${url}: ${result.stderr}`);
}

// Download the Google Fonts CSS and its woff2 files so the offline render uses the site faces.
function localFontCss() {
  const directory = path.join(out, 'fonts');
  fs.mkdirSync(directory, { recursive: true });
  const cached = path.join(directory, 'fonts.css');
  if (fs.existsSync(cached)) return fs.readFileSync(cached, 'utf8');
  try {
    const temporary = path.join(directory, 'google.css');
    curl(GOOGLE_FONTS_CSS, temporary);
    let text = fs.readFileSync(temporary, 'utf8');
    const urls = [...new Set([...text.matchAll(/url\((https:\/\/fonts\.gstatic\.com[^)]+)\)/g)].map((m) => m[1]))];
    for (const url of urls) {
      const name = url.split('/').slice(-2).join('-');
      const destination = path.join(directory, name);
      if (!fs.existsSync(destination)) curl(url, destination);
      text = text.split(url).join(`fonts/${name}`);
    }
    fs.writeFileSync(cached, text);
    return text;
  } catch (error) {
    console.error(`font download failed (${error.message}); falling back to the Google Fonts link`);
    return null;
  }
}

const embedFonts = (fontCss) => {
  let embedded = fontCss;
  for (const name of new Set([...fontCss.matchAll(/fonts\/([^)]+\.woff2)/g)].map((m) => m[1]))) {
    const data = fs.readFileSync(path.join(out, 'fonts', name)).toString('base64');
    embedded = embedded.split(`fonts/${name}`).join(`data:font/woff2;base64,${data}`);
  }
  return embedded;
};

const pageHtml = (fontCss, soundtrackSrc) => `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GammaLoop Showcase</title>
${fontCss ? `<style>${fontCss}</style>` : `<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n<link rel="stylesheet" href="${GOOGLE_FONTS_CSS}">`}
<style>${css}</style>
<style>
html, body { margin: 0; min-height: 100%; background: #211a23; }
.showcase-host { max-width: 1920px; margin: 0 auto; }
</style>
</head>
<body>
<div class="showcase-host" data-showcase data-autoplay data-poster="5.6" data-soundtrack="${soundtrackSrc}"></div>
<script>${js}</script>
</body>
</html>`;

const fontCss = flag('no-fonts') ? null : localFontCss();
const indexFile = path.join(out, 'index.html');
fs.copyFileSync(soundtrack, path.join(out, soundtrackName));
fs.writeFileSync(indexFile, pageHtml(fontCss, soundtrackName));
if (flag('standalone')) {
  const file = path.join(out, 'gammaloop-showcase.html');
  const audio = `data:audio/mpeg;base64,${fs.readFileSync(soundtrack).toString('base64')}`;
  fs.writeFileSync(file, pageHtml(fontCss ? embedFonts(fontCss) : null, audio));
  console.log(file);
  process.exit(0);
}
if (flag('preview')) {
  console.log(indexFile);
  process.exit(0);
}

async function loadPlaywright() {
  try {
    return await import('playwright');
  } catch {
    const root = spawnSync('npm', ['root', '-g'], { encoding: 'utf8' }).stdout.trim();
    return import(path.join(root, 'playwright', 'index.mjs'));
  }
}

function findFfmpeg() {
  const explicit = opt('ffmpeg', process.env.FFMPEG);
  if (explicit) return explicit;
  const which = spawnSync('sh', ['-c', 'command -v ffmpeg'], { encoding: 'utf8' });
  if (which.status === 0 && which.stdout.trim()) return which.stdout.trim();
  const python = spawnSync('python3', ['-c', 'import imageio_ffmpeg;print(imageio_ffmpeg.get_ffmpeg_exe())'], { encoding: 'utf8' });
  if (python.status === 0 && python.stdout.trim()) return python.stdout.trim();
  throw new Error('ffmpeg not found; pass --ffmpeg PATH or set FFMPEG');
}

const { chromium } = await loadPlaywright();
const ffmpegBin = findFfmpeg();
const browser = await chromium.launch();
const page = await browser.newPage({ viewport: { width, height }, deviceScaleFactor: 1 });
await page.goto(`file://${indexFile}?render=1`);
await page.evaluate(async () => {
  const faces = ['400', '500', '600'].map((w) => `${w} 22px "IBM Plex Mono"`)
    .concat(['400', '600', '700'].map((w) => `${w} 22px "Inter"`))
    .concat(['400 22px "Newsreader"', '600 22px "Newsreader"', 'italic 400 22px "Newsreader"']);
  await Promise.all(faces.map((face) => document.fonts.load(face).catch(() => null)));
  await document.fonts.ready;
  await window.showcase.ready;
});
const duration = await page.evaluate(() => window.showcase.duration);
const from = Number(opt('from', 0));
const to = Number(opt('to', duration));
const total = Math.round((to - from) * fps);

await page.evaluate((t) => window.showcase.seek(t), Number(opt('poster', 5.6)));
fs.writeFileSync(path.join(out, 'poster.png'), await page.screenshot({ type: 'png' }));

const videoFile = path.join(out, opt('name', 'gammaloop-showcase.mp4'));
const withMusic = !flag('no-music');
const ffmpegArgs = ['-y', '-loglevel', 'error', '-f', 'image2pipe', '-framerate', String(fps), '-i', 'pipe:0'];
if (withMusic) ffmpegArgs.push('-ss', String(from), '-i', soundtrack);
ffmpegArgs.push('-c:v', 'libx264', '-preset', 'medium', '-crf', '17', '-pix_fmt', 'yuv420p', '-movflags', '+faststart');
if (withMusic) ffmpegArgs.push('-c:a', 'aac', '-b:a', '160k', '-shortest');
ffmpegArgs.push(videoFile);
const ffmpeg = spawn(ffmpegBin, ffmpegArgs, { stdio: ['pipe', 'inherit', 'inherit'] });
const started = Date.now();
for (let i = 0; i < total; i++) {
  const t = from + i / fps;
  await page.evaluate((tt) => window.showcase.seek(tt), t);
  const frame = await page.screenshot({ type: 'png' });
  if (!ffmpeg.stdin.write(frame)) await once(ffmpeg.stdin, 'drain');
  if (i % Math.max(1, Math.round(total / 20)) === 0) {
    console.error(`frame ${i}/${total} t=${t.toFixed(2)}s elapsed ${((Date.now() - started) / 1000).toFixed(0)}s`);
  }
}
ffmpeg.stdin.end();
await once(ffmpeg, 'close');
await browser.close();
console.log(videoFile);
