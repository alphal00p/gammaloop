"""Execute every supported documentation example with pinned browser workers."""

import argparse
import functools
import http.server
import json
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--bundle", type=Path, required=True)
parser.add_argument("--runtime", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--chromium", type=Path, required=True)
args = parser.parse_args()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
runtime_manifest = json.loads((args.runtime / "runtime-manifest.json").read_text())
runtime_path = "assets/notebook-runtime/" + runtime_manifest["identity"]
for name, path in [
    ("bundle", args.bundle),
    (runtime_path, args.runtime),
    ("workers", args.bundle),
]:
    link = output / name
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise ValueError(f"Test staging path already exists: {link}")
    link.symlink_to(path.resolve(), target_is_directory=True)
for name, digest in runtime_manifest["shared_files"].items():
    link = output / "assets/notebook-runtime-files" / digest / Path(name).name
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.is_symlink():
        link.unlink()
    link.symlink_to((args.runtime / name).resolve())
(output / "runtime-config.json").write_text(
    json.dumps(
        {
            "runtime": runtime_path + "/",
            "wheel": "../../notebook-runtime-files/"
            + runtime_manifest["shared_files"][
                "wheels/symbolica-3.0.1-cp310-abi3-pyemscripten_2026_0_wasm32.whl"
            ]
            + "/symbolica-3.0.1-cp310-abi3-pyemscripten_2026_0_wasm32.whl",
            "stdlib": "../../notebook-runtime-files/"
            + runtime_manifest["shared_files"]["python_stdlib.zip"]
            + "/python_stdlib.zip",
        }
    )
)


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


server = http.server.ThreadingHTTPServer(
    ("127.0.0.1", 0), functools.partial(Quiet, directory=str(output))
)
threading.Thread(target=server.serve_forever, daemon=True).start()
try:
    with sync_playwright() as p:
        browser = p.chromium.launch(
            executable_path=str(args.chromium), headless=True, args=["--no-sandbox"]
        )
        context = browser.new_context()
        denied = []

        def route_request(route):
            if route.request.url.startswith(f"http://127.0.0.1:{server.server_port}/"):
                route.continue_()
            else:
                denied.append(route.request.url)
                route.abort()

        context.route("**/*", route_request)
        page = context.new_page()
        page.goto(f"http://127.0.0.1:{server.server_port}/")
        page.on("console", lambda msg: print(msg.text, flush=True))
        result = page.evaluate("""async () => {
          const root = new URL('./',location.href);
          const config = {...await(await fetch('runtime-config.json')).json(),compilerRoot:new URL('bundle/',root).href,
            files:await(await fetch('bundle/python-files.json')).json(),fonts:await(await fetch('bundle/fonts.json')).json()};
          config.runtime = new URL(config.runtime,root).href;
          const measure = (await import(new URL('bundle/graph-measure.js',root))).measure;
          const cells = await(await fetch('bundle/cells.json')).json();
          cells.push(...await(await fetch('bundle/checks.json')).json());
          const sourceFor = (cell, seen = new Set()) => {
            if(seen.has(cell.id)) return []; seen.add(cell.id);
            return [...(cell.depends||[]).flatMap(id=>sourceFor(cells.find(c=>c.id===id),seen)),
              ...(cell.setup?[[cell.id+':setup',cell.setup]]:[]),[cell.id,cell.source]];
          };
          const results=[];
          for (const language of ['python','typst']) {
            const worker = new Worker(new URL(`workers/documentation-${language}-worker.js`,root),{type:'module'});
            const pending = new Map();
            await new Promise((resolve,reject)=>{
              worker.onerror=e=>reject(Error(e.message));
              worker.onmessage=async ({data})=>{
                if(data.type==='measure') {
                  try { worker.postMessage({type:'measurement',id:data.id,receipt:await measure(data.request)}); }
                  catch(e) {worker.postMessage({type:'measurement',id:data.id,error:String(e)});}
                } else if(data.type==='ready')resolve();
                else if(data.type==='error')reject(Error(data.message));
                else {const done=pending.get(data.id);pending.delete(data.id);done?.(data);}
              };
              worker.postMessage({type:'init',config});
            });
            for (const cell of cells.filter(c=>c.language===language&&!c.reason)) {
              const request = language==='python'?{sources:sourceFor(cell),session:cell.file,packages:cell.packages||[]}:
                {source:(cell.prefix||'')+(cell.execution||cell.source),path:cell.path};
              const result = await new Promise(resolve=>{pending.set(cell.id,resolve);worker.postMessage({type:'run',id:cell.id,...request});});
              results.push({...result,outputs:result.success?result.outputs.map(o=>({mimes:Object.keys(o)})):result.outputs});console.log(cell.id,result.success,result.success?'':JSON.stringify(result.outputs).slice(-2000));
            }
            if(language==='typst') {
              for(const [id,source] of [['invalid','#let ='],['recovery','Hello recovery']]) {
                const result=await new Promise(resolve=>{pending.set(id,resolve);worker.postMessage({type:'run',id,source});});
                results.push(result);console.log(id,result.success);
              }
            }
            worker.terminate();
          }
          return results;
        }""")
        (output / "runtime-results.json").write_text(json.dumps(result, indent=2))
        if denied:
            raise AssertionError(f"External runtime requests: {denied}")
        failures = [
            r for r in result if r.get("success", False) != (r["id"] != "invalid")
        ]
        print(
            "RESULT",
            len(result) - len(failures),
            "/",
            len(result),
            "expected outcomes",
            flush=True,
        )
        if failures:
            raise AssertionError(
                f"Documentation execution failures: {[r['id'] for r in failures]}"
            )
        browser.close()
finally:
    server.shutdown()
    server.server_close()
