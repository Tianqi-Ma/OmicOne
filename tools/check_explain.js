#!/usr/bin/env node
// Run every step animation through a mock canvas, in both languages, over two
// full loops: no exception, no non-finite coordinate, and every registered
// step has a scene. Run after any change to inst/app/www/explain*.js:
//
//   node tools/check_explain.js            (from the repository root)
//
// Exit 0 = all scenes fine, 1 = a problem (printed). Visual review: open
// docs/explain-gallery.html in a browser.
"use strict";
const fs = require("fs");
const path = require("path");

const root = path.resolve(__dirname, "..");
const www = path.join(root, "inst", "app", "www");
const files = fs.readdirSync(www).filter((f) => /^explain(-[a-z]+)?\.js$/.test(f))
  .sort((a, b) => (a === "explain.js" ? -1 : b === "explain.js" ? 1 : a.localeCompare(b)));

let nonFinite = 0;
let calls = 0;
let where = "";
function mockCtx() {
  const target = {};
  return new Proxy(target, {
    get(t, k) {
      if (k in t) return t[k];
      if (k === "measureText") return (s) => ({ width: String(s).length * 6 });
      return (...args) => {
        calls++;
        if (args.some((a) => typeof a === "number" && !Number.isFinite(a))) {
          nonFinite++;
          if (nonFinite <= 5) console.log("non-finite argument to", String(k), "in", where);
        }
      };
    },
    set(t, k, v) { t[k] = v; return true; }
  });
}

global.window = global;
global.addEventListener = () => {};
global.document = {
  readyState: "complete", visibilityState: "visible",
  documentElement: { getAttribute: () => "light" },
  querySelectorAll: () => [], addEventListener() {}, body: { contains: () => true }
};
global.getComputedStyle = () => ({ getPropertyValue: () => "" });
global.requestAnimationFrame = () => 0;
global.matchMedia = () => ({ matches: false });

for (const f of files) {
  // the global is renamed per package by tools/sync_mirror.sh
  eval(fs.readFileSync(path.join(www, f), "utf8"));
}
const X = global.OmicOneExplain || global.scStudioExplain;
if (!X) { console.log("explain.js did not register its global"); process.exit(1); }

// step keys of this repo, read from the registries
const keys = [];
for (const f of ["steps_sc.R", "steps.R"]) {
  const p = path.join(root, "R", f);
  if (!fs.existsSync(p)) continue;
  const src = fs.readFileSync(p, "utf8");
  for (const m of src.matchAll(/list\(v = "([a-z_]+)"[^)]*?ui = mod_(?!placeholder)/g)) keys.push(m[1]);
}

let errors = 0;
for (const lang of ["en", "zh"]) {
  global.OmicOneLang = global.scStudioLang = lang;
  for (const name of X.names()) {
    const canvas = { clientWidth: 560, width: 0, height: 0, getContext: () => mockCtx() };
    const stages = X.stagesOf(name) || [];
    const period = stages.length ? stages[stages.length - 1][1] : 10;
    try {
      for (let t = 0; t <= period * 2; t += 0.1) {
        where = `${name} (${lang}, t=${t.toFixed(1)})`;
        X.drawAt(canvas, name, t % period);
      }
    } catch (e) {
      errors++;
      console.log("ERROR", where, e.message);
    }
  }
}
const missing = keys.filter((k) => !X.has(k));
console.log(`${X.names().length} scenes; ${calls} draw calls; ` +
            `${errors} errors; ${nonFinite} non-finite; missing: ${missing.join(", ") || "none"}`);
process.exit(errors || nonFinite || missing.length ? 1 : 0);
