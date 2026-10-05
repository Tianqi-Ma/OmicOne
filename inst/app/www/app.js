/* OmicOne client-side helpers: language switch + startup animation */
(function () {
  "use strict";

  // ---------------------------------------------------------------------------
  // Language switch — only ONE language visible at a time.
  // Each translatable element has class "i18n" and data-en / data-zh holding the
  // two versions (may be HTML). We swap innerHTML to the active language.
  // ---------------------------------------------------------------------------
  window.OmicOneLang = "en";

  function applyLang(root) {
    var lang = window.OmicOneLang;
    var scope = root || document;
    var nodes = scope.querySelectorAll(".i18n");
    nodes.forEach(function (el) {
      var val = el.getAttribute("data-" + lang);
      if (val === null) val = el.getAttribute("data-en");
      if (val !== null && el.innerHTML !== val) el.innerHTML = val;
    });
  }

  window.OmicOneSetLang = function (lang) {
    window.OmicOneLang = (lang === "zh") ? "zh" : "en";
    document.body.setAttribute("data-lang", window.OmicOneLang);
    document.documentElement.setAttribute("lang", window.OmicOneLang === "zh" ? "zh-CN" : "en");
    applyLang(document);
    // reflect active state on the segmented control
    document.querySelectorAll(".omicone-lang-btn").forEach(function (b) {
      b.classList.toggle("active", b.getAttribute("data-lang") === window.OmicOneLang);
    });
    try { window.localStorage.setItem("omicone-lang", window.OmicOneLang); } catch (e) {}
    sendLang();
  };

  // Let the server know too (input$lang), e.g. for text drawn inside plots.
  function sendLang() {
    if (window.Shiny && typeof Shiny.setInputValue === "function") {
      Shiny.setInputValue("lang", window.OmicOneLang);
    }
  }
  document.addEventListener("shiny:connected", sendLang);

  // Start in the language chosen last time; on a first visit, follow the
  // browser's language.
  function initialLang() {
    try {
      var saved = window.localStorage.getItem("omicone-lang");
      if (saved === "en" || saved === "zh") return saved;
    } catch (e) {}
    var nav = (navigator.languages && navigator.languages[0]) || navigator.language || "en";
    return /^zh/i.test(nav) ? "zh" : "en";
  }
  function initLang() { window.OmicOneSetLang(initialLang()); }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initLang);
  } else {
    initLang();
  }

  // Re-apply language to any freshly rendered (renderUI) content.
  document.addEventListener("shiny:idle", function () { applyLang(document); });
  document.addEventListener("shiny:value", function (e) {
    if (e.target) setTimeout(function () { applyLang(e.target); }, 0);
  });

  // ---------------------------------------------------------------------------
  // Startup animation — scattered "cells" fly in and coalesce into coloured
  // UMAP-like clusters, then the logo fades in; the splash fades out on load.
  // ---------------------------------------------------------------------------
  function runSplash() {
    var splash = document.getElementById("omicone-splash");
    if (!splash) return;
    var canvas = document.getElementById("omicone-splash-canvas");
    if (!canvas || !canvas.getContext) { fadeOut(splash); return; }
    var ctx = canvas.getContext("2d");
    var W, H;
    function resize() {
      W = canvas.width = splash.clientWidth;
      H = canvas.height = splash.clientHeight;
    }
    resize();
    window.addEventListener("resize", resize);

    // brand palette (blue-cyan research theme + accents)
    var pal = ["#2f81c7", "#4f9fd0", "#3fb37f", "#e4572e", "#b5179e",
               "#f4a261", "#9c6ade", "#2a9d8f"];
    var K = 6;                     // number of UMAP-like clusters
    var N = 340;                   // number of cells
    var cx = W / 2, cy = H / 2 - 20;
    var R = Math.min(W, H) * 0.26;
    var centers = [];
    for (var k = 0; k < K; k++) {
      var a = (k / K) * Math.PI * 2;
      centers.push({ x: cx + Math.cos(a) * R, y: cy + Math.sin(a) * R,
                     col: pal[k % pal.length] });
    }
    var cells = [];
    for (var i = 0; i < N; i++) {
      var c = centers[i % K];
      cells.push({
        sx: Math.random() * W, sy: Math.random() * H,             // start (scattered)
        tx: c.x + (Math.random() - 0.5) * R * 0.7,                // target (in cluster)
        ty: c.y + (Math.random() - 0.5) * R * 0.7,
        col: c.col, r: 1.5 + Math.random() * 2.2
      });
    }

    var start = null, DUR = 1600;
    function ease(t) { return 1 - Math.pow(1 - t, 3); }          // easeOutCubic
    function frame(ts) {
      if (start === null) start = ts;
      var p = Math.min(1, (ts - start) / DUR);
      var e = ease(p);
      ctx.clearRect(0, 0, W, H);
      for (var i = 0; i < cells.length; i++) {
        var c = cells[i];
        var x = c.sx + (c.tx - c.sx) * e;
        var y = c.sy + (c.ty - c.sy) * e;
        ctx.globalAlpha = 0.35 + 0.55 * e;
        ctx.fillStyle = c.col;
        ctx.beginPath();
        ctx.arc(x, y, c.r, 0, Math.PI * 2);
        ctx.fill();
      }
      ctx.globalAlpha = 1;
      if (p < 1) {
        requestAnimationFrame(frame);
      } else {
        var logo = document.getElementById("omicone-splash-logo");
        if (logo) logo.classList.add("show");
        setTimeout(function () { fadeOut(splash); }, 750);
      }
    }
    requestAnimationFrame(frame);
    // safety: never let the splash trap the user
    setTimeout(function () { fadeOut(splash); }, 5000);
  }

  function fadeOut(splash) {
    if (!splash || splash.classList.contains("hide")) return;
    splash.classList.add("hide");
    setTimeout(function () { if (splash && splash.parentNode) splash.style.display = "none"; }, 600);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", runSplash);
  } else {
    runSplash();
  }

  // ---------------------------------------------------------------------------
  // Landing-page omics cards — each has a themed canvas animation that plays on
  // hover (mouseenter starts a requestAnimationFrame loop; mouseleave stops it).
  // ---------------------------------------------------------------------------
  var PAL = ["#2f81c7", "#4f9fd0", "#3fb37f", "#e4572e", "#b5179e",
             "#f4a261", "#9c6ade", "#2a9d8f"];

  function cardAnim(canvas, kind) {
    var ctx = canvas.getContext && canvas.getContext("2d");
    if (!ctx) return { stop: function () {} };
    var W, H, raf = null, t0 = null, running = true;
    function resize() { W = canvas.width = canvas.clientWidth; H = canvas.height = canvas.clientHeight; }
    resize();
    var pts = [];
    function seed() {
      pts = [];
      var n = 90;
      for (var i = 0; i < n; i++) {
        pts.push({ x: Math.random() * W, y: Math.random() * H,
                   c: PAL[i % PAL.length], r: 1.5 + Math.random() * 2,
                   a: Math.random() * Math.PI * 2, sp: 0.4 + Math.random() * 1.2,
                   col: i % 5, row: i % 6 });
      }
    }
    seed();
    function frame(ts) {
      if (!running) return;
      if (t0 === null) t0 = ts;
      var e = (ts - t0) / 1000;
      ctx.clearRect(0, 0, W, H);
      var cx = W / 2, cy = H / 2, i, p, x, y;
      if (kind === "sc") {                     // cells drift into cluster blobs
        for (i = 0; i < pts.length; i++) {
          p = pts[i];
          var ang = (p.col / 5) * Math.PI * 2 + e * 0.3;
          var tx = cx + Math.cos(ang) * W * 0.22, ty = cy + Math.sin(ang) * H * 0.22;
          p.x += (tx - p.x) * 0.03; p.y += (ty - p.y) * 0.03;
          dot(ctx, p.x + Math.sin(e + i) * 3, p.y + Math.cos(e + i) * 3, p.r, p.c);
        }
      } else if (kind === "bulk") {            // volcano: points spread outward
        for (i = 0; i < pts.length; i++) {
          p = pts[i];
          var side = (i % 2 === 0) ? 1 : -1;
          x = cx + side * ((i * 7 % (W / 2)) * (0.5 + 0.5 * Math.sin(e)));
          y = H - ((i * 13) % H) * (0.5 + 0.5 * Math.abs(Math.sin(e * 0.8 + i)));
          dot(ctx, x, y, 2.2, side > 0 ? "#e4572e" : "#2f81c7");
        }
      } else if (kind === "wes") {             // oncoplot grid fills tile by tile
        var cols = 12, rows = 7, cw = W / cols, ch = H / rows, k = Math.floor(e * 18) % (cols * rows);
        for (i = 0; i <= k; i++) {
          var cc = i % cols, rr = Math.floor(i / cols) % rows;
          ctx.fillStyle = PAL[(cc + rr) % PAL.length];
          ctx.globalAlpha = 0.85;
          ctx.fillRect(cc * cw + 1, rr * ch + 1, cw - 2, ch - 2);
        }
        ctx.globalAlpha = 1;
      } else if (kind === "spatial") {         // hex spots light up in a gradient
        var r = Math.min(W, H) / 16;
        for (i = 0; i < pts.length; i++) {
          var gx = cx + ((i % 9) - 4) * r * 1.7 + ((Math.floor(i / 9) % 2) ? r * 0.85 : 0);
          var gy = cy + (Math.floor(i / 9) - 3) * r * 1.5;
          var lit = 0.5 + 0.5 * Math.sin(e * 1.5 - (gx + gy) * 0.02);
          hex(ctx, gx, gy, r * 0.75, PAL[0], lit);
        }
      } else {                                  // integration: streams merge to node
        for (i = 0; i < pts.length; i++) {
          p = pts[i];
          var lane = i % 4;
          var sx = (lane < 2 ? 0 : W), sy = (lane % 2 === 0 ? 0 : H);
          var prog = (e * 0.5 + i / pts.length) % 1;
          x = sx + (cx - sx) * prog; y = sy + (cy - sy) * prog;
          dot(ctx, x, y, 2, PAL[lane]);
        }
        dot(ctx, cx, cy, 6 + 2 * Math.sin(e * 3), "#4f9fd0");
      }
      raf = requestAnimationFrame(frame);
    }
    function dot(c, x, y, r, col) { c.fillStyle = col; c.beginPath(); c.arc(x, y, r, 0, Math.PI * 2); c.fill(); }
    function hex(c, x, y, r, col, a) {
      c.globalAlpha = a; c.fillStyle = col; c.beginPath();
      for (var j = 0; j < 6; j++) { var ang = Math.PI / 3 * j + Math.PI / 6;
        var px = x + r * Math.cos(ang), py = y + r * Math.sin(ang);
        j === 0 ? c.moveTo(px, py) : c.lineTo(px, py); }
      c.closePath(); c.fill(); c.globalAlpha = 1;
    }
    raf = requestAnimationFrame(frame);
    return { stop: function () { running = false; if (raf) cancelAnimationFrame(raf); if (ctx) ctx.clearRect(0, 0, W, H); } };
  }

  function wireCards() {
    document.querySelectorAll(".omicone-omcard").forEach(function (card) {
      if (card.dataset.wired) return;
      card.dataset.wired = "1";
      var canvas = card.querySelector(".omicone-omcanvas");
      var kind = card.getAttribute("data-anim");
      var anim = null;
      card.addEventListener("mouseenter", function () { if (canvas && !anim) anim = cardAnim(canvas, kind); });
      card.addEventListener("mouseleave", function () { if (anim) { anim.stop(); anim = null; } });
    });
  }
  // Cards are rendered dynamically (renderUI); wire them whenever the DOM settles.
  document.addEventListener("shiny:idle", wireCards);
  document.addEventListener("shiny:value", function () { setTimeout(wireCards, 30); });
})();
