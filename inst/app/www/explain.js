/* OmicOne step animations — engine.
 *
 * Every pipeline step has a short looping "scene" that shows, in pictures, what
 * the step does. A scene is drawn on <canvas class="omicone-explain"
 * data-scene="<step key>"> (rendered by preview_plot_ui() in the empty preview
 * area before the first run, and by explain_scene() elsewhere).
 *
 * No dependencies and nothing fetched: the app must work offline. Scenes are
 * drawn in a fixed 560 x 260 logical space and scaled to the canvas, so a scene
 * is written once and looks the same at any size. Text follows the EN / 中
 * switch, colours follow the light / dark theme, a canvas only animates while
 * it is on screen, and prefers-reduced-motion shows a single still frame.
 *
 * Scenes are registered by explain-sc.js / explain-wes.js:
 *   OmicOneExplain.register("qc", {
 *     period: 8,                                   // seconds per loop
 *     stages: [[t0, t1, "English", "中文"], ...],   // captions, top-left
 *     still: 6.5,                                  // reduced-motion frame
 *     init: function (R, H) { return state; },     // R = seeded random()
 *     draw: function (g, s, t, H) { ... }          // g = drawing helper
 *   });
 */
(function () {
  "use strict";

  var W = 560, HH = 260, TAU = Math.PI * 2;
  var FONT = 'Inter, system-ui, -apple-system, "Segoe UI", "PingFang SC", "Microsoft YaHei", sans-serif';
  var MONO = 'ui-monospace, SFMono-Regular, Menlo, Consolas, "Courier New", monospace';

  // ---- small maths ---------------------------------------------------------
  function clamp(x, a, b) { return x < a ? a : x > b ? b : x; }
  function lerp(a, b, p) { return a + (b - a) * p; }
  function easeIO(p) { return p < 0.5 ? 4 * p * p * p : 1 - Math.pow(-2 * p + 2, 3) / 2; }
  function easeOut(p) { return 1 - Math.pow(1 - p, 3); }
  // progress of t through [a, b], eased (pass false for linear)
  function seg(t, a, b, ease) {
    var p = clamp((t - a) / (b - a), 0, 1);
    if (ease === false) return p;
    return (ease || easeIO)(p);
  }
  function rng(seed) {                       // mulberry32: deterministic per scene
    var a = seed >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function gauss(R) {
    var u = 0, v = 0;
    while (u === 0) u = R();
    while (v === 0) v = R();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(TAU * v);
  }
  function hash(s) {
    var h = 2166136261;
    for (var i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); }
    return h >>> 0;
  }
  // "#rrggbb", "#rgb" or "rgb(r,g,b)" (mix() returns the latter, and mixes get mixed again)
  function hexRgb(c) {
    c = String(c).trim();
    var m = /^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i.exec(c);
    if (m) return [parseInt(m[1], 16), parseInt(m[2], 16), parseInt(m[3], 16)];
    m = /^#([0-9a-f])([0-9a-f])([0-9a-f])$/i.exec(c);    // short form, e.g. a theme's #fff
    if (m) return [parseInt(m[1] + m[1], 16), parseInt(m[2] + m[2], 16), parseInt(m[3] + m[3], 16)];
    m = /^rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)/i.exec(c);
    if (m) return [+m[1], +m[2], +m[3]];
    return [128, 128, 128];
  }
  // mix two hex colours (p = 0 -> c1, 1 -> c2)
  function mix(c1, c2, p) {
    var a = hexRgb(c1), b = hexRgb(c2);
    return "rgb(" + Math.round(lerp(a[0], b[0], p)) + "," + Math.round(lerp(a[1], b[1], p)) +
      "," + Math.round(lerp(a[2], b[2], p)) + ")";
  }
  var VIRIDIS = ["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"];
  function ramp(stops, p) {
    p = clamp(p, 0, 1) * (stops.length - 1);
    var i = Math.min(stops.length - 2, Math.floor(p));
    return mix(stops[i], stops[i + 1], p - i);
  }
  // a Gaussian blob of points
  function blob(R, n, cx, cy, sx, sy, extra) {
    var out = [];
    for (var i = 0; i < n; i++) {
      var p = { x: cx + gauss(R) * sx, y: cy + gauss(R) * (sy == null ? sx : sy) };
      if (extra) for (var k in extra) if (Object.prototype.hasOwnProperty.call(extra, k)) p[k] = extra[k];
      out.push(p);
    }
    return out;
  }

  // ---- language & theme ----------------------------------------------------
  function lang() { return window.OmicOneLang === "zh" ? "zh" : "en"; }
  function L(en, zh) { return lang() === "zh" ? (zh == null ? en : zh) : en; }
  function readTheme(el) {
    var cs = window.getComputedStyle(el);
    function v(n, d) { var x = cs.getPropertyValue(n); x = x && x.trim(); return x || d; }
    var dark = document.documentElement.getAttribute("data-bs-theme") === "dark";
    return {
      dark: dark,
      text: v("--sc-text", "#1c2530"), muted: v("--sc-muted", "#66727f"),
      border: v("--sc-border", "#e3e7ed"), card: v("--sc-card", "#ffffff"),
      panel: v("--sc-panel", "#f7f9fb"), accent: v("--sc-accent", "#2f81c7"),
      accent2: v("--sc-accent2", "#2669a8"), ok: v("--sc-ok", "#2f9e6d"),
      warn: v("--sc-warn", "#b7791f"),
      faint: dark ? "#3a4654" : "#d5dbe3",      // empty tiles, neutral marks
      ink: dark ? "#e6edf3" : "#1c2530"
    };
  }

  // ---- drawing helper --------------------------------------------------------
  function Gfx(ctx, th) { this.c = ctx; this.th = th; this.A = 1; this.W = W; this.H = HH; }
  var G = Gfx.prototype;
  G.alpha = function (a) { this.c.globalAlpha = clamp(this.A * (a == null ? 1 : a), 0, 1); };
  G.dot = function (x, y, r, col, a) {
    if (!(r > 0)) return;
    var c = this.c; this.alpha(a); c.fillStyle = col;
    c.beginPath(); c.arc(x, y, r, 0, TAU); c.fill();
  };
  G.ring = function (x, y, r, col, w, a) {
    if (!(r > 0)) return;
    var c = this.c; this.alpha(a); c.strokeStyle = col; c.lineWidth = w || 1;
    c.beginPath(); c.arc(x, y, r, 0, TAU); c.stroke();
  };
  G.arc = function (x, y, r, a0, a1, col, w, a) {
    var c = this.c; this.alpha(a); c.strokeStyle = col; c.lineWidth = w || 1;
    c.beginPath(); c.arc(x, y, r, a0, a1); c.stroke();
  };
  G.line = function (x1, y1, x2, y2, col, w, a, dash) {
    var c = this.c; this.alpha(a); c.strokeStyle = col; c.lineWidth = w || 1;
    c.setLineDash(dash || []); c.beginPath(); c.moveTo(x1, y1); c.lineTo(x2, y2); c.stroke();
    c.setLineDash([]);
  };
  G.poly = function (pts, col, w, a, dash) {
    if (!pts || pts.length < 2) return;
    var c = this.c; this.alpha(a); c.strokeStyle = col; c.lineWidth = w || 1;
    c.lineJoin = "round"; c.lineCap = "round";
    c.setLineDash(dash || []); c.beginPath(); c.moveTo(pts[0][0], pts[0][1]);
    for (var i = 1; i < pts.length; i++) c.lineTo(pts[i][0], pts[i][1]);
    c.stroke(); c.setLineDash([]);
  };
  G.fillPoly = function (pts, col, a) {
    if (!pts || pts.length < 3) return;
    var c = this.c; this.alpha(a); c.fillStyle = col; c.beginPath(); c.moveTo(pts[0][0], pts[0][1]);
    for (var i = 1; i < pts.length; i++) c.lineTo(pts[i][0], pts[i][1]);
    c.closePath(); c.fill();
  };
  G.curve = function (x1, y1, cx, cy, x2, y2, col, w, a) {      // quadratic
    var c = this.c; this.alpha(a); c.strokeStyle = col; c.lineWidth = w || 1;
    c.beginPath(); c.moveTo(x1, y1); c.quadraticCurveTo(cx, cy, x2, y2); c.stroke();
  };
  function roundRect(c, x, y, w, h, r) {
    r = Math.max(0, Math.min(r || 0, Math.abs(w) / 2, Math.abs(h) / 2));
    c.beginPath();
    c.moveTo(x + r, y); c.lineTo(x + w - r, y); c.quadraticCurveTo(x + w, y, x + w, y + r);
    c.lineTo(x + w, y + h - r); c.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
    c.lineTo(x + r, y + h); c.quadraticCurveTo(x, y + h, x, y + h - r);
    c.lineTo(x, y + r); c.quadraticCurveTo(x, y, x + r, y); c.closePath();
  }
  G.rect = function (x, y, w, h, col, a, r) {
    if (!(w > 0) || !(h > 0)) return;
    var c = this.c; this.alpha(a); c.fillStyle = col; roundRect(c, x, y, w, h, r); c.fill();
  };
  G.box = function (x, y, w, h, col, lw, a, r, dash) {
    var c = this.c; this.alpha(a); c.strokeStyle = col; c.lineWidth = lw || 1;
    c.setLineDash(dash || []); roundRect(c, x, y, w, h, r); c.stroke(); c.setLineDash([]);
  };
  G.text = function (s, x, y, o) {
    o = o || {};
    var c = this.c; this.alpha(o.a);
    c.fillStyle = o.col || this.th.text;
    c.font = (o.italic ? "italic " : "") + (o.bold ? "600 " : "") + (o.size || 11) + "px " + (o.mono ? MONO : FONT);
    c.textAlign = o.align || "left"; c.textBaseline = o.base || "middle";
    if (o.rot) {
      c.save(); c.translate(x, y); c.rotate(o.rot); c.fillText(s, 0, 0); c.restore();
    } else {
      c.fillText(s, x, y);
    }
  };
  G.t = function (en, zh, x, y, o) { this.text(L(en, zh), x, y, o); };
  G.measure = function (s, size, bold) {
    this.c.font = (bold ? "600 " : "") + (size || 11) + "px " + FONT;
    return this.c.measureText(s).width;
  };
  G.arrow = function (x1, y1, x2, y2, col, w, a, head) {
    var dx = x2 - x1, dy = y2 - y1, len = Math.sqrt(dx * dx + dy * dy);
    if (len < 0.5) return;
    var hs = head == null ? Math.min(8, len * 0.45) : head;
    var ux = dx / len, uy = dy / len;
    this.line(x1, y1, x2 - ux * hs * 0.6, y2 - uy * hs * 0.6, col, w, a);
    this.fillPoly([[x2, y2], [x2 - ux * hs - uy * hs * 0.5, y2 - uy * hs + ux * hs * 0.5],
                   [x2 - ux * hs + uy * hs * 0.5, y2 - uy * hs - ux * hs * 0.5]], col, a);
  };
  // L-shaped axes; (x0, y0) is the origin (bottom-left)
  G.axes = function (x0, y0, w, h, xl, yl, a) {
    var th = this.th;
    this.line(x0, y0, x0 + w, y0, th.muted, 1, a == null ? 0.7 : a);
    this.line(x0, y0, x0, y0 - h, th.muted, 1, a == null ? 0.7 : a);
    if (xl) this.text(xl, x0 + w / 2, y0 + 13, { col: th.muted, size: 10, align: "center", a: a });
    if (yl) this.text(yl, x0 - 12, y0 - h / 2, { col: th.muted, size: 10, align: "center",
                                                 rot: -Math.PI / 2, a: a });
  };
  G.chip = function (s, x, y, col, a, o) {
    o = o || {};
    var size = o.size || 10.5, w = this.measure(s, size, o.bold) + 14, h = size + 9;
    var x0 = o.align === "left" ? x : x - w / 2;
    this.rect(x0, y - h / 2, w, h, o.bg || col, (a == null ? 1 : a) * (o.bg ? 1 : 0.14), h / 2);
    this.box(x0, y - h / 2, w, h, col, 1, (a == null ? 1 : a) * 0.6, h / 2);
    this.text(s, x0 + w / 2, y + 0.5, { col: o.fg || col, size: size, align: "center",
                                        bold: o.bold, a: a });
    return w;
  };
  G.note = function (en, zh, a) {            // footnote line, bottom-right
    this.text(L(en, zh), W - 12, HH - 10, { col: this.th.muted, size: 9.5, align: "right",
                                             italic: true, a: a });
  };

  // ---- stage caption ---------------------------------------------------------
  function drawStages(g, stages, t) {
    if (!stages) return;
    for (var i = 0; i < stages.length; i++) {
      var s = stages[i];
      if (t < s[0] - 0.05 || t > s[1] + 0.05) continue;
      var a = Math.min(seg(t, s[0], s[0] + 0.35, false), 1 - seg(t, s[1] - 0.3, s[1], false));
      if (a <= 0) continue;
      var x = 14, y = 15;
      g.dot(x + 7, y, 7.5, g.th.accent, a);
      g.text(String(i + 1), x + 7, y + 0.5, { col: "#ffffff", size: 9.5, bold: true,
                                               align: "center", a: a });
      g.text(L(s[2], s[3]), x + 20, y, { col: g.th.text, size: 11.5, bold: true, a: a });
    }
  }

  // ---- registry & runtime ----------------------------------------------------
  var SCENES = {};
  var H = { clamp: clamp, lerp: lerp, seg: seg, easeIO: easeIO, easeOut: easeOut,
            gauss: gauss, blob: blob, mix: mix, ramp: ramp, VIRIDIS: VIRIDIS,
            L: L, lang: lang, TAU: TAU, W: W, H: HH,
            PAL: ["#2f81c7", "#e4572e", "#3fb37f", "#b5179e", "#f4a261", "#4361ee",
                  "#2a9d8f", "#9c6ade", "#e63946", "#ffca3a"],
            KEPT: "#3b6ea5", FLAG: "#c1476b" };

  var items = [];
  var reduced = !!(window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches);
  var raf = null;

  function attach(el) {
    if (el.__omExplain) return;
    var name = el.getAttribute("data-scene");
    var sc = SCENES[name];
    var guide = el.closest ? el.closest(".omicone-preview-guide, .omicone-explain-wrap") : null;
    if (!sc) { el.style.display = "none"; return; }
    if (guide) guide.classList.add("has-scene");
    var it = { el: el, sc: sc, name: name, ctx: el.getContext("2d"), state: null,
               t0: null, visible: false, th: null, thAge: 0, w: 0, still: false };
    it.state = sc.init ? sc.init(rng(hash(name)), H) : {};
    el.__omExplain = it;
    items.push(it);
    if ("IntersectionObserver" in window) {
      if (!attach.io) {
        attach.io = new IntersectionObserver(function (entries) {
          entries.forEach(function (e) {
            var x = e.target.__omExplain;
            if (x) x.visible = e.isIntersecting;
          });
          kick();
        }, { threshold: 0.05 });
      }
      attach.io.observe(el);
    } else {
      it.visible = true;
    }
    kick();
  }

  function sizeCanvas(it) {
    var cw = it.el.clientWidth;
    if (!cw) return false;
    var dpr = Math.min(window.devicePixelRatio || 1, 2.5);
    if (cw !== it.w || it.dpr !== dpr) {
      it.w = cw;
      it.dpr = dpr;
      it.el.width = Math.round(cw * dpr);
      it.el.height = Math.round(cw * HH / W * dpr);
    }
    var k = it.el.width / W;
    it.ctx.setTransform(k, 0, 0, k, 0, 0);
    return true;
  }

  function drawItem(it, now) {
    if (!sizeCanvas(it)) return;
    if (!it.th || now - it.thAge > 800) { it.th = readTheme(it.el); it.thAge = now; }
    var sc = it.sc, P = sc.period || 8;
    if (it.t0 === null) it.t0 = now;
    var t = reduced ? (sc.still != null ? sc.still : P * 0.8) : ((now - it.t0) / 1000) % P;
    var g = new Gfx(it.ctx, it.th);
    it.ctx.clearRect(0, 0, W, HH);
    g.A = reduced ? 1 : Math.min(seg(t, 0, 0.35, false), 1 - seg(t, P - 0.45, P, false));
    try {
      sc.draw(g, it.state, t, H);
      drawStages(g, sc.stages, t);
    } catch (e) {
      if (window.console) console.warn("explain scene '" + it.name + "':", e);
      it.broken = true;
    }
    it.ctx.globalAlpha = 1;
  }

  function frame(now) {
    raf = null;
    var any = false;
    for (var i = items.length - 1; i >= 0; i--) {
      var it = items[i];
      if (!document.body.contains(it.el)) {          // UI re-rendered: drop it
        if (attach.io) attach.io.unobserve(it.el);
        items.splice(i, 1);
        continue;
      }
      if (!it.visible || it.broken || it.el.offsetParent === null) continue;
      any = true;
      if (reduced && it.drawnAt && now - it.drawnAt < 1000) continue;
      drawItem(it, now);
      it.drawnAt = now;
    }
    if (any && document.visibilityState !== "hidden") raf = window.requestAnimationFrame(frame);
  }
  function kick() { if (raf === null) raf = window.requestAnimationFrame(frame); }

  function scan(root) {
    var list = (root || document).querySelectorAll
      ? (root || document).querySelectorAll("canvas.omicone-explain") : [];
    for (var i = 0; i < list.length; i++) attach(list[i]);
  }

  // canvases arrive with renderUI(), so watch the DOM
  function start() {
    scan(document);
    if ("MutationObserver" in window) {
      new MutationObserver(function (muts) {
        for (var i = 0; i < muts.length; i++) {
          var added = muts[i].addedNodes;
          for (var j = 0; j < added.length; j++) {
            var n = added[j];
            if (n.nodeType !== 1) continue;
            if (n.matches && n.matches("canvas.omicone-explain")) attach(n);
            else if (n.querySelector && n.querySelector("canvas.omicone-explain")) scan(n);
          }
        }
        kick();
      }).observe(document.body, { childList: true, subtree: true });
    }
    // switching tabs/steps changes visibility without a DOM insert
    document.addEventListener("shown.bs.tab", kick);
    document.addEventListener("click", function () { setTimeout(kick, 50); });
    document.addEventListener("visibilitychange", kick);
    window.addEventListener("resize", kick);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();

  window.OmicOneExplain = {
    register: function (name, scene) { SCENES[name] = scene; scan(document); },
    has: function (name) { return !!SCENES[name]; },
    names: function () { return Object.keys(SCENES); },
    helpers: H,
    rescan: function () { scan(document); kick(); },
    stagesOf: function (name) { return SCENES[name] ? SCENES[name].stages || null : null; },
    // drawing one frame at time t: used by the gallery and the tests
    drawAt: function (canvas, name, t) {
      var sc = SCENES[name];
      if (!sc) return false;
      var it = { el: canvas, ctx: canvas.getContext("2d"), w: 0 };
      if (!sizeCanvas(it)) { canvas.width = W; canvas.height = HH; it.ctx.setTransform(1, 0, 0, 1, 0, 0); }
      var g = new Gfx(it.ctx, readTheme(canvas));
      it.ctx.clearRect(0, 0, W, HH);
      sc.draw(g, sc.init ? sc.init(rng(hash(name)), H) : {}, t, H);
      drawStages(g, sc.stages, t);
      return true;
    }
  };
})();
