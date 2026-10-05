/* OmicOne step animations — WES / somatic mutation scenes (one per step key).
 * Drawing space: 560 x 260; captions occupy the top 30 px. See explain.js. */
(function () {
  "use strict";
  var X = window.OmicOneExplain;
  if (!X) return;
  var H = X.helpers, seg = H.seg, lerp = H.lerp, PAL = H.PAL, L = H.L;

  // maftools' default variant-classification colours, so the scenes match the plots
  var VC = { Missense: "#33A02C", Nonsense: "#E31A1C", Frame_Shift: "#1F78B4",
             Splice: "#FB9A99", In_Frame: "#CAB2D6", Multi_Hit: "#555555" };
  var VCN = ["Missense", "Nonsense", "Frame_Shift", "Splice", "In_Frame"];
  function vcPick(R) {
    var u = R();
    return u < 0.6 ? "Missense" : u < 0.75 ? "Nonsense" : u < 0.88 ? "Frame_Shift" : u < 0.95 ? "Splice" : "In_Frame";
  }
  // COSMIC SBS colours; C>G is drawn in the theme ink so it shows in dark mode
  function sbsCols(th) { return ["#1EBFF0", th.ink, "#E62725", "#CBCACB", "#A1CF64", "#EDC8C5"]; }
  var SBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"];

  function km(R, n, rate, horizon) {
    var times = [];
    for (var i = 0; i < n; i++) times.push(Math.min(horizon, -Math.log(1 - R()) / rate));
    times.sort(function (a, b) { return a - b; });
    var s = 1, pts = [[0, 1]], cens = [];
    times.forEach(function (tt, i) {
      if (tt >= horizon || R() < 0.25) { cens.push([tt, s]); return; }
      pts.push([tt, s]); s *= 1 - 1 / (n - i); pts.push([tt, s]);
    });
    pts.push([horizon, s]);
    return { pts: pts, cens: cens };
  }
  function drawKM(g, k, x0, y0, w, h, horizon, col, p) {
    var out = [], lim = horizon * p, last = 1;
    for (var i = 0; i < k.pts.length; i++) {
      var q = k.pts[i];
      if (q[0] > lim) { out.push([x0 + lim / horizon * w, y0 - last * h]); break; }
      out.push([x0 + q[0] / horizon * w, y0 - q[1] * h]);
      last = q[1];
    }
    g.poly(out, col, 2.2, 1);
    k.cens.forEach(function (c) {
      if (c[0] > lim) return;
      g.line(x0 + c[0] / horizon * w, y0 - c[1] * h - 3.5, x0 + c[0] / horizon * w, y0 - c[1] * h + 3.5, col, 1.4, 0.9);
    });
  }

  // ---- 1. import MAF ---------------------------------------------------------
  X.register("wes_import", {
    period: 8.5,
    stages: [[0, 2.6, "A MAF file lists one mutation per row", "MAF 文件每行记录一个突变"],
             [2.6, 6, "Rows are gathered per sample and gene", "按样本与基因汇总"],
             [6, 8.5, "The basis of the cohort's mutation landscape", "得到队列突变全景的基础"]],
    still: 7.2,
    init: function (R) {
      var genes = ["TP53", "NOTCH1", "PIK3CA", "KMT2D"], rows = [];
      for (var i = 0; i < 9; i++) rows.push({ s: Math.floor(R() * 6), g: Math.floor(R() * 4), v: vcPick(R) });
      return { genes: genes, rows: rows };
    },
    draw: function (g, s, t) {
      var th = g.th, tx = 26, ty = 52;
      g.text("Tumor_Sample_Barcode   Hugo_Symbol   Variant_Classification", tx, ty - 10,
             { col: th.muted, size: 8.5, mono: true, a: seg(t, 0, 0.6) });
      s.rows.forEach(function (r, i) {
        var a = seg(t, 0.2 + i * 0.2, 0.6 + i * 0.2), fly = seg(t, 2.8 + i * 0.25, 3.8 + i * 0.25);
        var y = ty + 6 + i * 19;
        g.text("S" + (r.s + 1), tx, y, { size: 9.5, mono: true, a: a * (1 - fly * 0.5) });
        g.text(s.genes[r.g], tx + 70, y, { size: 9.5, mono: true, a: a * (1 - fly * 0.5) });
        g.rect(tx + 140, y - 6, 60, 12, VC[r.v], a * 0.85 * (1 - fly * 0.5), 3);
        g.text(r.v.replace("_", " "), tx + 170, y, { col: "#ffffff", size: 8, align: "center", a: a * (1 - fly * 0.5) });
        if (fly > 0 && fly < 1) {
          var gx = 340 + r.s * 30 + 14, gy = 70 + r.g * 34 + 14;
          g.rect(lerp(tx + 140, gx - 12, fly), lerp(y - 6, gy - 12, fly), lerp(60, 24, fly), lerp(12, 24, fly), VC[r.v], 0.85, 3);
        }
      });
      var ga = seg(t, 2.6, 3.2);
      s.genes.forEach(function (gn, i) { g.text(gn, 334, 70 + i * 34 + 14, { col: th.muted, size: 10, align: "right", a: ga }); });
      for (var c = 0; c < 6; c++) g.text("S" + (c + 1), 340 + c * 30 + 14, 62, { col: th.muted, size: 9, align: "center", a: ga });
      for (var i = 0; i < 4; i++) for (c = 0; c < 6; c++) g.box(342 + c * 30, 72 + i * 34, 24, 24, th.faint, 1, ga, 3);
      s.rows.forEach(function (r, i) {
        if (seg(t, 2.8 + i * 0.25, 3.8 + i * 0.25) < 1) return;
        g.rect(342 + r.s * 30, 72 + r.g * 34, 24, 24, VC[r.v], 0.9, 3);
      });
      g.note("Silent variants are kept but not counted as non-synonymous", "同义突变会保留，但不计入非同义突变", seg(t, 6.2, 7));
    }
  });

  // ---- 2. cohort summary -----------------------------------------------------
  X.register("wes_summary", {
    period: 8.5,
    stages: [[0, 3, "Which kinds of variants are there?", "有哪些类型的变异？"],
             [3, 6, "How many per sample?", "每个样本有多少？"],
             [6, 8.5, "Check outliers before going further", "继续之前先检查离群样本"]],
    still: 7.4,
    init: function (R) {
      var samples = [];
      for (var i = 0; i < 22; i++) samples.push(i === 3 ? 160 : 18 + Math.floor(R() * 50));
      samples.sort(function (a, b) { return b - a; });
      return { samples: samples };
    },
    draw: function (g, s, t) {
      var th = g.th, ba = seg(t, 0.2, 1.8);
      var val = [1, 0.3, 0.22, 0.14, 0.08];
      VCN.forEach(function (v, i) {
        var y = 58 + i * 34;
        g.text(v.replace("_", " "), 110, y + 8, { col: th.muted, size: 10, align: "right", a: ba });
        g.rect(116, y, 140 * val[i] * ba, 16, VC[v], 0.9, 3);
      });
      var x0 = 300, y0 = 222, w = 240, h = 170, sa = seg(t, 3, 4.6, false);
      g.axes(x0, y0, w, h, L("samples", "样本"), L("variants", "变异数"), seg(t, 2.8, 3.4));
      s.samples.forEach(function (n, i) {
        if (i >= sa * s.samples.length) return;
        var bh = Math.min(h - 10, n) * 0.95, bx = x0 + 4 + i * 10.6;
        g.rect(bx, y0 - bh, 8, bh * 0.62, VC.Missense, 0.9);
        g.rect(bx, y0 - bh * 0.38 - bh * 0.62, 8, bh * 0.38, VC.Frame_Shift, 0.9);
      });
      var med = 40 * 0.95, ma = seg(t, 4.6, 5.2);
      g.line(x0, y0 - med, x0 + w, y0 - med, th.ink, 1, ma * 0.7, [4, 3]);
      g.t("median", "中位数", x0 + w, y0 - med - 8, { col: th.muted, size: 9.5, align: "right", a: ma });
      var oa = seg(t, 6.2, 7);
      g.ring(x0 + 8, y0 - 150, 14, H.FLAG, 1.6, oa);
      g.t("hypermutated? (POLE / MMR-d is real biology)", "超突变？（POLE / MMR 缺陷是真实生物学）",
          x0 + 26, y0 - 150, { col: H.FLAG, size: 9.5, a: oa });
    }
  });

  // ---- 3. oncoplot -----------------------------------------------------------
  X.register("wes_onco", {
    period: 10.5,
    stages: [[0, 3.2, "Each tile: a mutated gene in one sample", "每个方块：某个样本中一个突变的基因"],
             [3.2, 6.8, "Sort genes by frequency, samples by pattern", "基因按频率排序，样本按突变模式排序"],
             [6.8, 10.5, "Bars: burden per sample, % mutated per gene", "柱：每个样本的负荷，每个基因的突变比例"]],
    still: 9.4,
    init: function (R) {
      var genes = ["TP53", "NOTCH1", "PIK3CA", "KMT2D", "FAT1", "CDKN2A"], p = [0.75, 0.4, 0.28, 0.22, 0.18, 0.12];
      var n = 20, M = [];
      for (var g = 0; g < genes.length; g++) {
        M.push([]);
        for (var s = 0; s < n; s++) M[g].push(R() < p[g] ? (R() < 0.08 ? "Multi_Hit" : vcPick(R)) : null);
      }
      var ord = genes.map(function (gn, i) { return i; }).sort(function (a, b) {
        return M[b].filter(Boolean).length - M[a].filter(Boolean).length;
      });
      genes = ord.map(function (i) { return genes[i]; });
      M = ord.map(function (i) { return M[i]; });
      var idx = [];
      for (s = 0; s < n; s++) idx.push(s);
      var start = idx.slice().sort(function () { return R() - 0.5; });
      var sorted = idx.slice().sort(function (a, b) {
        for (var g2 = 0; g2 < genes.length; g2++) {
          var d = (M[g2][b] ? 1 : 0) - (M[g2][a] ? 1 : 0);
          if (d) return d;
        }
        return a - b;
      });
      var tmb = idx.map(function () { return 0.2 + R() * 0.8; });
      return { genes: genes, M: M, n: n, start: start, sorted: sorted, tmb: tmb };
    },
    draw: function (g, s, t) {
      var th = g.th, x0 = 90, y0 = 70, cw = 18.5, rh = 24, so = seg(t, 3.6, 6);
      s.genes.forEach(function (gn, i) {
        g.text(gn, x0 - 6, y0 + i * rh + rh / 2, { col: th.text, size: 10, align: "right", italic: true, a: seg(t, 0, 0.6) });
      });
      for (var c = 0; c < s.n; c++) {
        var from = s.start.indexOf(c), to = s.sorted.indexOf(c), col = lerp(from, to, so);
        var x = x0 + col * cw;
        for (var r = 0; r < s.genes.length; r++) {
          var v = s.M[r][c], y = y0 + r * rh, a = seg(t, 0.2 + (from + r) * 0.06, 0.6 + (from + r) * 0.06);
          g.rect(x + 1, y + 1, cw - 2, rh - 2, th.faint, a * 0.6, 2);
          if (v) g.rect(x + 1, y + 1, cw - 2, rh - 2, VC[v], a * 0.95, 2);
        }
        var ta = seg(t, 7, 8);
        g.rect(x + 2, y0 - 6 - s.tmb[c] * 26 * ta, cw - 4, s.tmb[c] * 26 * ta, th.muted, 0.7, 1);
      }
      var fa = seg(t, 7.4, 8.4);
      s.genes.forEach(function (gn, i) {
        var f = s.M[i].filter(Boolean).length / s.n;
        g.rect(x0 + s.n * cw + 10, y0 + i * rh + 4, 60 * f * fa, rh - 8, th.accent, 0.8, 2);
        g.text(Math.round(f * 100) + "%", x0 + s.n * cw + 14 + 60 * f * fa, y0 + i * rh + rh / 2, { col: th.muted, size: 9.5, a: fa });
      });
    }
  });

  // ---- 4. Ti/Tv & VAF --------------------------------------------------------
  X.register("wes_titv", {
    period: 9.5,
    stages: [[0, 3, "Six classes of single-base substitution", "六类单碱基替换"],
             [3, 6.2, "Transitions (C>T, T>C) vs transversions", "转换（C>T、T>C）与颠换"],
             [6.2, 9.5, "VAF: the fraction of reads carrying the variant", "VAF：携带该变异的读段比例"]],
    still: 8.4,
    init: function (R) {
      var vaf = [[0.42, 0.06], [0.2, 0.08], [0.16, 0.07], [0.3, 0.1]].map(function (m) {
        var v = []; for (var i = 0; i < 14; i++) v.push(H.clamp(m[0] + H.gauss(R) * m[1], 0.02, 0.95)); v.sort(); return v;
      });
      return { vaf: vaf };
    },
    draw: function (g, s, t) {
      var th = g.th, cols = sbsCols(th), val = [0.16, 0.08, 0.42, 0.06, 0.2, 0.08], ba = seg(t, 0.2, 1.8);
      var x0 = 50, y0 = 215;
      g.axes(x0, y0, 220, 170, "", L("fraction", "比例"), seg(t, 0, 0.6));
      SBS.forEach(function (lab, i) {
        var bh = val[i] * 340 * ba, x = x0 + 10 + i * 35;
        var ti = i === 2 || i === 4, dim = seg(t, 3.2, 4) * (ti ? 0 : 0.55);
        g.rect(x, y0 - bh, 26, bh, cols[i], 0.9 - dim, 2);
        g.text(lab, x + 13, y0 + 10, { col: th.muted, size: 9.5, align: "center", a: ba });
      });
      var tl = seg(t, 3.4, 4.2);
      g.t("Ti", "转换", x0 + 10 + 2 * 35 + 13, y0 - 0.42 * 340 - 10, { col: cols[2], size: 10.5, bold: true, align: "center", a: tl });
      g.t("Ti", "转换", x0 + 10 + 4 * 35 + 13, y0 - 0.2 * 340 - 10, { col: cols[4], size: 10.5, bold: true, align: "center", a: tl });
      g.note("C>T and C>A also come from FFPE / oxidative damage", "C>T、C>A 也可能来自 FFPE 或氧化损伤", tl);
      var va = seg(t, 6.4, 7.4), vx = 330, vy = 215, vw = 200, vh = 170;
      if (va > 0) {
        g.axes(vx, vy, vw, vh, "", "VAF", va);
        ["TP53", "NOTCH1", "PIK3CA", "KMT2D"].forEach(function (gn, k) {
          var cx = vx + 26 + k * 48;
          s.vaf[k].forEach(function (v, i) { g.dot(cx + ((i * 7) % 13 - 6), vy - v * vh * va, 2.4, PAL[k], 0.75); });
          var med = s.vaf[k][7];
          g.line(cx - 12, vy - med * vh * va, cx + 12, vy - med * vh * va, th.ink, 2, va);
          g.text(gn, cx, vy + 10, { col: th.muted, size: 9, align: "center", italic: true, a: va });
        });
      }
    }
  });

  // ---- 5. TMB ---------------------------------------------------------------
  X.register("wes_tmb", {
    period: 9.5,
    stages: [[0, 3, "Count non-synonymous mutations in the captured region", "统计捕获区域内的非同义突变数"],
             [3, 6, "Divide by the size of that region in megabases", "除以该区域的大小（Mb）"],
             [6, 9.5, "Compare samples on a log scale", "在对数坐标上比较各样本"]],
    still: 8.4,
    init: function (R) {
      var ticks = [], samples = [];
      for (var i = 0; i < 42; i++) ticks.push(R());
      for (i = 0; i < 40; i++) samples.push({ v: Math.exp(H.gauss(R) * 0.9 + Math.log(2.5)) * (i === 5 ? 15 : 1), j: H.gauss(R) });
      return { ticks: ticks, samples: samples };
    },
    draw: function (g, s, t) {
      var th = g.th, sx = 60, sw = 440, sy = 70, ea = seg(t, 0.1, 0.8) * (1 - seg(t, 5.8, 6.4));
      g.rect(sx, sy, sw, 16, th.faint, ea * 0.7, 4);
      g.t("captured exome (e.g. 35.8 Mb)", "捕获的外显子区域（如 35.8 Mb）", sx, sy - 12, { col: th.muted, size: 10, a: ea });
      var n = Math.round(seg(t, 0.6, 2.8, false) * s.ticks.length);
      for (var i = 0; i < n; i++) g.line(sx + s.ticks[i] * sw, sy - 3, sx + s.ticks[i] * sw, sy + 19, H.FLAG, 1.6, ea);
      g.text(String(n), sx + sw + 8, sy + 8, { col: H.FLAG, size: 12, bold: true, a: ea });
      var fa = seg(t, 3.1, 3.9) * (1 - seg(t, 5.8, 6.4));
      g.t("42 mutations ÷ 35.8 Mb = 1.2 mut/Mb", "42 个突变 ÷ 35.8 Mb = 1.2 mut/Mb", 280, 150,
          { size: 16, bold: true, align: "center", col: th.accent2, a: fa });
      g.t("the denominator must match your capture kit", "分母必须与所用捕获试剂盒一致", 280, 178,
          { col: th.muted, size: 10.5, align: "center", a: fa });
      var pa = seg(t, 6.2, 7.2), x0 = 90, y0 = 222, h = 175;
      if (pa > 0) {
        g.axes(x0, y0, 420, h, L("samples, sorted", "样本（已排序）"), "", pa);
        g.text("mut/Mb (log)", x0 - 34, y0 - h / 2, { col: th.muted, size: 10, align: "center", rot: -Math.PI / 2, a: pa });
        var lg = function (v) { return y0 - (Math.log10(v) + 1) / 3 * h; };
        [0.1, 1, 10, 100].forEach(function (v) { g.text(String(v), x0 - 4, lg(v), { col: th.muted, size: 8.5, align: "right", a: pa }); });
        var srt = s.samples.slice().sort(function (a, b) { return a.v - b.v; });
        srt.forEach(function (p, i) {
          var hi = p.v >= 10, hl = seg(t, 7.6, 8.4);
          g.dot(x0 + 10 + i * 10, lg(p.v), hi ? 3.6 : 3, hi ? H.mix(th.accent, H.FLAG, hl) : th.accent, pa * 0.85);
        });
        g.line(x0, lg(10), x0 + 420, lg(10), th.warn, 1.4, pa, [5, 4]);
        g.t("10 mut/Mb (panel-based cut-off, not calibrated for WES)", "10 mut/Mb（基于 panel 的阈值，未针对 WES 校准）",
            x0 + 6, lg(10) - 9, { col: th.warn, size: 9.5, a: pa });
      }
    }
  });

  // ---- 6. lollipop ------------------------------------------------------------
  X.register("wes_lolli", {
    period: 9.5,
    stages: [[0, 2.6, "The protein, with its domains", "蛋白及其结构域"],
             [2.6, 6.2, "Each lollipop: mutations at one residue", "每根棒棒糖：一个氨基酸位点上的突变"],
             [6.2, 9.5, "Tall stacks are hotspots", "高高的堆叠就是热点"]],
    still: 8.4,
    init: function (R) {
      var muts = [[175, 8, "Missense"], [248, 10, "Missense"], [273, 9, "Missense"], [220, 4, "Missense"],
                  [245, 5, "Missense"], [196, 2, "Nonsense"], [213, 2, "Nonsense"], [306, 2, "Nonsense"]];
      for (var i = 0; i < 9; i++) muts.push([20 + Math.floor(R() * 370), 1, vcPick(R)]);
      return { muts: muts };
    },
    draw: function (g, s, t) {
      var th = g.th, x0 = 60, w = 450, y = 196, len = 393, a = seg(t, 0.1, 1);
      function px(aa) { return x0 + aa / len * w; }
      g.rect(x0, y, w * a, 10, th.faint, 1, 3);
      [[6, 30, "TAD", PAL[4]], [95, 289, "DNA-binding", PAL[0]], [318, 358, "Tet", PAL[2]]].forEach(function (d) {
        var da = seg(t, 0.8, 1.6);
        g.rect(px(d[0]), y - 4, px(d[1]) - px(d[0]), 18, d[3], da * 0.85, 4);
        g.text(d[2], (px(d[0]) + px(d[1])) / 2, y + 5, { col: "#ffffff", size: 9, bold: true, align: "center", a: da });
      });
      g.text("0", x0, y + 24, { col: th.muted, size: 9, align: "center", a: a });
      g.text("393 aa", x0 + w, y + 24, { col: th.muted, size: 9, align: "center", a: a });
      g.text("TP53", x0, 46, { size: 12, bold: true, italic: true, a: a });
      s.muts.forEach(function (m, i) {
        var p = seg(t, 2.8 + i * 0.12, 3.6 + i * 0.12);
        if (p <= 0) return;
        var x = px(m[0]), top = y - 8 - m[1] * 13 * p;
        g.line(x, y - 4, x, top, th.muted, 1, p);
        g.dot(x, top, 4.2, VC[m[2]], p);
      });
      var la = seg(t, 6.4, 7.2);
      [["R175H", 175, 8], ["R248Q", 248, 10], ["R273H", 273, 9]].forEach(function (h) {
        g.text(h[0], px(h[1]), y - 18 - h[2] * 13, { col: th.text, size: 9.5, bold: true, align: "center", a: la });
      });
      g.note("Positions depend on the transcript the MAF was annotated with", "位置取决于 MAF 注释所用的转录本", la);
    }
  });

  // ---- 7. drivers & interactions ---------------------------------------------
  X.register("wes_driver", {
    period: 9.5,
    stages: [[0, 3, "Drivers tend to hit the same positions", "驱动突变倾向于命中相同位点"],
             [3, 6, "Passengers scatter along the protein", "乘客突变沿蛋白随机分布"],
             [6, 9.5, "Score positional clustering, control the FDR", "评估位置聚集程度，并控制 FDR"]],
    still: 8.4,
    init: function (R) {
      var a = [], b = [], genes = [];
      for (var i = 0; i < 14; i++) a.push(i < 10 ? 0.45 + H.gauss(R) * 0.012 : R());
      for (i = 0; i < 14; i++) b.push(R());
      for (i = 0; i < 16; i++) genes.push({ x: R() * 0.5, y: R() * 1.1, s: 2 + R() * 3 });
      genes.push({ x: 0.75, y: 3.4, s: 7, lab: "GENE A" });
      return { a: a, b: b, genes: genes };
    },
    draw: function (g, s, t) {
      var th = g.th, x0 = 40, w = 300;
      function stack(arr, y, p, col) {
        var bins = {};
        arr.forEach(function (u, i) {
          var k = Math.round(u * 60), n = bins[k] = (bins[k] || 0) + 1;
          g.dot(x0 + u * w, y - 8 - (n - 1) * 8 * p, 3.4, col, seg(p * 14, i, i + 1));
        });
      }
      var aa = seg(t, 0.2, 2.6, false), ba = seg(t, 3.2, 5.4, false);
      g.rect(x0, 96, w, 9, th.faint, seg(t, 0, 0.5), 3);
      g.text("GENE A", x0, 60, { size: 11, bold: true, italic: true, a: seg(t, 0, 0.5) });
      stack(s.a, 96, aa, H.FLAG);
      g.t("clustered → driver-like", "聚集 → 疑似驱动", x0 + w, 60, { col: H.FLAG, size: 10, align: "right", a: seg(t, 2, 2.6) });
      g.rect(x0, 196, w, 9, th.faint, seg(t, 3, 3.5), 3);
      g.text("GENE B", x0, 160, { size: 11, bold: true, italic: true, a: seg(t, 3, 3.5) });
      stack(s.b, 196, ba, th.muted);
      g.t("scattered → passenger-like", "分散 → 疑似乘客", x0 + w, 160, { col: th.muted, size: 10, align: "right", a: seg(t, 5, 5.6) });
      var sa = seg(t, 6.2, 7.2), sx = 395, sy = 215, sw = 140, sh = 165;
      if (sa > 0) {
        g.axes(sx, sy, sw, sh, L("fraction in clusters", "聚集比例"), "−log10 FDR", sa);
        g.line(sx, sy - 1 / 4 * sh, sx + sw, sy - 1 / 4 * sh, th.warn, 1.2, sa, [4, 3]);
        s.genes.forEach(function (q) {
          var sig = q.y > 1;
          g.dot(sx + q.x * sw, sy - q.y / 4 * sh, q.s, sig ? H.FLAG : th.muted, sa * 0.75);
          if (q.lab) g.text(q.lab, sx + q.x * sw - 8, sy - q.y / 4 * sh - 12, { size: 9, bold: true, align: "center", a: sa });
        });
      }
      g.note("Scattered truncating drivers (tumour suppressors) are missed", "散在截短型驱动（抑癌基因）会被漏检", seg(t, 7, 7.8));
    }
  });

  // ---- 8. mutational signatures ----------------------------------------------
  X.register("wes_sig", {
    period: 10.5,
    stages: [[0, 3, "96 trinucleotide contexts form the spectrum", "96 种三核苷酸背景构成突变谱"],
             [3, 6.8, "NMF splits it into a few signatures", "NMF 把它分解为少数几个特征"],
             [6.8, 10.5, "Match each one to the COSMIC catalogue", "与 COSMIC 目录逐一比对"]],
    still: 9.4,
    init: function (R) {
      function sig(fn) { var v = []; for (var i = 0; i < 96; i++) v.push(Math.max(0.01, fn(Math.floor(i / 16), i % 16) + R() * 0.04)); return v; }
      var A = sig(function (cls, ctx) { return cls === 2 && ctx % 4 === 2 ? 0.9 : 0.05; });          // clock-like CpG C>T
      var B = sig(function (cls, ctx) { return (cls === 2 || cls === 1) && Math.floor(ctx / 4) === 3 ? 0.8 : 0.04; }); // APOBEC TCN
      var C = sig(function (cls) { return cls === 0 ? 0.55 + R() * 0.3 : 0.06; });                   // tobacco C>A
      var w = [0.45, 0.3, 0.25], mixv = [];
      for (var i = 0; i < 96; i++) mixv.push(w[0] * A[i] + w[1] * B[i] + w[2] * C[i]);
      return { A: A, B: B, C: C, w: w, mix: mixv };
    },
    draw: function (g, s, t) {
      var th = g.th, cols = sbsCols(th);
      function spec(v, x0, y0, w, h, a, scale) {
        var bw = w / 96, mx = scale || Math.max.apply(null, v);
        for (var i = 0; i < 96; i++) g.rect(x0 + i * bw, y0 - v[i] / mx * h, Math.max(0.8, bw - 0.6), v[i] / mx * h, cols[Math.floor(i / 16)], a);
        for (var k = 0; k < 6; k++) g.rect(x0 + k * 16 * bw, y0 + 2, 16 * bw - 1, 3, cols[k], a);
      }
      var ma = seg(t, 0.2, 1.6);
      spec(s.mix, 40, 120, 480, 72, ma);
      SBS.forEach(function (lab, k) { g.text(lab, 40 + (k + 0.5) * 80, 132, { col: th.muted, size: 9, align: "center", a: ma }); });
      var da = seg(t, 3.2, 4.4), names = [["SBS1 · clock-like", "SBS1 · 时钟样"], ["SBS2/13 · APOBEC", "SBS2/13 · APOBEC"],
                                           ["SBS4 · tobacco", "SBS4 · 烟草"]];
      [s.A, s.B, s.C].forEach(function (v, k) {
        var x0 = 40 + k * 168, p = seg(t, 3.4 + k * 0.4, 4.4 + k * 0.4);
        spec(v, x0, 222, 144, 46, p);
        g.text((s.w[k] * 100).toFixed(0) + "%", x0 + 72, 160, { col: th.text, size: 11, bold: true, align: "center", a: p });
        if (k < 2) g.text("+", x0 + 156, 200, { col: th.muted, size: 16, align: "center", a: da });
        var ca = seg(t, 7 + k * 0.4, 7.8 + k * 0.4);
        g.text(L(names[k][0], names[k][1]) + "  cos 0.9", x0 + 72, 244, { col: th.accent2, size: 9.5, align: "center", a: ca });
      });
      g.t("Small cohorts: refit known signatures instead of extracting new ones", "小队列：用已知特征拟合，而不是提取新特征",
          548, 40, { col: th.muted, size: 9.5, italic: true, align: "right", a: seg(t, 8.4, 9.2) });
    }
  });

  // ---- 9. clinical / pathway / drug -------------------------------------------
  X.register("wes_clin", {
    period: 9.5,
    stages: [[0, 3, "Split samples by a clinical feature", "按临床特征划分样本"],
             [3, 6.2, "Compare each gene's mutation rate (Fisher test)", "逐个基因比较突变率（Fisher 检验）"],
             [6.2, 9.5, "Keep genes that pass the FDR", "保留通过 FDR 的基因"]],
    still: 8.4,
    init: function (R) {
      var sm = [];
      for (var i = 0; i < 24; i++) sm.push({ g: i % 2, j: Math.floor(i / 2), x: R(), y: R() });
      return { sm: sm };
    },
    draw: function (g, s, t) {
      var th = g.th, sp = seg(t, 1, 2.4), GC = [PAL[0], PAL[4]];
      s.sm.forEach(function (p, i) {
        var x = lerp(60 + p.x * 160, 52 + p.g * 112 + (p.j % 4) * 18, sp), y = lerp(60 + p.y * 150, 74 + Math.floor(p.j / 4) * 26, sp);
        g.dot(x, y, 5, GC[p.g], seg(t, 0.1, 0.8) * (1 - seg(t, 3, 3.6)));
      });
      g.t("Stage I–II", "I–II 期", 80, 50, { col: GC[0], size: 10.5, bold: true, align: "center", a: sp * (1 - seg(t, 3, 3.6)) });
      g.t("Stage III–IV", "III–IV 期", 190, 50, { col: GC[1], size: 10.5, bold: true, align: "center", a: sp * (1 - seg(t, 3, 3.6)) });
      var ba = seg(t, 3.2, 4.4), genes = ["TP53", "NOTCH1", "PIK3CA", "FAT1", "NFE2L2"], r = [[0.8, 0.85], [0.2, 0.52], [0.25, 0.28], [0.15, 0.19], [0.1, 0.12]];
      genes.forEach(function (gn, i) {
        var y = 56 + i * 36, sig = i === 1, q = seg(t, 6.4, 7.2);
        g.text(gn, 300, y + 10, { col: th.text, size: 10, align: "right", italic: true, a: ba * (sig || q < 1 ? 1 : 1 - 0.5 * q) });
        g.rect(308, y, 200 * r[i][0] * ba, 9, GC[0], 0.85 * (sig ? 1 : 1 - 0.55 * q), 2);
        g.rect(308, y + 11, 200 * r[i][1] * ba, 9, GC[1], 0.85 * (sig ? 1 : 1 - 0.55 * q), 2);
        if (sig) g.text("q < 0.05", 312 + 200 * r[i][1] + 4, y + 10, { col: H.FLAG, size: 10, bold: true, a: q });
      });
      g.t("% samples mutated", "突变样本比例", 408, 240, { col: th.muted, size: 9.5, align: "center", a: ba });
    }
  });

  // ---- 10. cohort comparison --------------------------------------------------
  X.register("wes_compare", {
    period: 9.5,
    stages: [[0, 3, "Two cohorts, the same genes", "两个队列，同一批基因"],
             [3, 6.2, "Odds ratio of being mutated, with 95% CI", "突变的比值比及 95% 置信区间"],
             [6.2, 9.5, "Only FDR-significant genes are coloured", "只有 FDR 显著的基因被着色"]],
    still: 8.4,
    draw: function (g, s, t) {
      var th = g.th, x0 = 150, x1 = 510, cx = (x0 + x1) / 2, y0 = 62;
      function lx(or) { return cx + Math.log2(or) / 4 * (x1 - x0) / 2; }
      var rows = [["TP53", 1.1, 0.8, 1.5, false], ["NOTCH1", 3.2, 1.8, 5.6, true], ["PIK3CA", 0.9, 0.5, 1.6, false],
                  ["KMT2D", 0.3, 0.15, 0.6, true], ["FAT1", 1.6, 0.7, 3.4, false], ["CDKN2A", 0.7, 0.3, 1.5, false]];
      var a = seg(t, 0.1, 1);
      g.line(cx, y0 - 14, cx, y0 + rows.length * 28, th.muted, 1.2, a, [4, 3]);
      [0.25, 0.5, 1, 2, 4].forEach(function (v) { g.text(String(v), lx(v), y0 + rows.length * 28 + 12, { col: th.muted, size: 9, align: "center", a: a }); });
      g.t("OR (log scale)", "OR（对数刻度）", cx, y0 + rows.length * 28 + 26, { col: th.muted, size: 9.5, align: "center", a: a });
      g.t("← higher in cohort B", "← 队列 B 更高", x0, 44, { col: PAL[4], size: 10, a: a });
      g.t("higher in cohort A →", "队列 A 更高 →", x1, 44, { col: PAL[0], size: 10, align: "right", a: a });
      rows.forEach(function (r, i) {
        var y = y0 + i * 28 + 8, p = seg(t, 3.2 + i * 0.25, 4.2 + i * 0.25), c = seg(t, 6.4, 7.2);
        g.text(r[0], x0 - 12, y, { col: th.text, size: 10, align: "right", italic: true, a: a });
        var col = r[4] ? H.mix(th.muted, r[1] > 1 ? PAL[0] : PAL[4], c) : th.muted;
        g.line(lerp(cx, lx(r[2]), p), y, lerp(cx, lx(r[3]), p), y, col, 2, p);
        g.rect(lerp(cx, lx(r[1]), p) - 5, y - 5, 10, 10, col, p, 2);
      });
    }
  });

  // ---- 11. mutation vs survival -----------------------------------------------
  X.register("wes_surv", {
    period: 10.5,
    stages: [[0, 3, "Mutant: carries a mutation in the gene", "突变型：该基因携带突变"],
             [3, 5.8, "Wild-type includes sequenced samples with no mutation at all", "野生型包括已测序但完全没有突变记录的样本"],
             [5.8, 10.5, "Compare survival (Kaplan–Meier, log-rank)", "比较生存（Kaplan–Meier，log-rank）"]],
    still: 9.4,
    init: function (R) {
      var s = [];
      for (var i = 0; i < 24; i++) s.push({ m: i < 7, x: R(), y: R() });
      return { s: s, mut: km(R, 30, 1 / 10, 60), wt: km(R, 60, 1 / 26, 60) };
    },
    draw: function (g, s, t) {
      var th = g.th, sp = seg(t, 1, 2.4), add = seg(t, 3.2, 4.6);
      s.s.forEach(function (p, i) {
        var j = p.m ? i : i - 7;
        var tx = p.m ? 60 + (j % 4) * 20 : 160 + (j % 5) * 20, ty = 70 + Math.floor(j / (p.m ? 4 : 5)) * 22;
        g.dot(lerp(50 + p.x * 200, tx, sp), lerp(60 + p.y * 150, ty, sp), 5.5, p.m ? H.FLAG : PAL[0], seg(t, 0.1, 0.8) * (1 - seg(t, 5.8, 6.4)));
      });
      for (var k = 0; k < 4; k++) {
        g.dot(lerp(250, 160 + ((17 + k) % 5) * 20, add), lerp(220, 70 + Math.floor((17 + k) / 5) * 22, add), 5.5,
              H.mix(th.muted, PAL[0], add), seg(t, 3, 3.5) * (1 - seg(t, 5.8, 6.4)));
      }
      var la = seg(t, 1.6, 2.4) * (1 - seg(t, 5.8, 6.4));
      g.t("mutant", "突变型", 90, 52, { col: H.FLAG, size: 10.5, bold: true, align: "center", a: la });
      g.t("wild-type", "野生型", 200, 52, { col: PAL[0], size: 10.5, bold: true, align: "center", a: la });
      g.t("no mutation in the MAF → WT", "MAF 中无突变记录 → 野生型", 250, 232, { col: th.muted, size: 9.5, align: "center", a: seg(t, 3.2, 3.8) * (1 - seg(t, 5.8, 6.4)) });
      var kp = seg(t, 6.2, 9, false), kx = 90, ky = 215, kw = 400, kh = 165;
      if (kp > 0) {
        g.axes(kx, ky, kw, kh, L("months", "月"), L("survival", "生存率"), seg(t, 6, 6.6));
        drawKM(g, s.wt, kx, ky, kw, kh, 60, PAL[0], kp);
        drawKM(g, s.mut, kx, ky, kw, kh, 60, H.FLAG, kp);
        g.t("log-rank test · HR (95% CI)", "log-rank 检验 · HR（95% CI）", kx + kw, ky - kh + 6, { col: th.muted, size: 10, align: "right", a: seg(t, 8.6, 9.2) });
      }
    }
  });

  // ---- 12. heterogeneity -------------------------------------------------------
  X.register("wes_hetero", {
    period: 9.5,
    stages: [[0, 3, "Each variant has an allele frequency (VAF)", "每个变异都有一个等位基因频率（VAF）"],
             [3, 6.2, "VAFs pile up into clusters", "VAF 聚集成若干簇"],
             [6.2, 9.5, "A main clone plus subclones", "一个主克隆加上若干亚克隆"]],
    still: 8.4,
    init: function (R) {
      var v = [];
      for (var i = 0; i < 70; i++) v.push({ x: i < 44 ? H.clamp(0.4 + H.gauss(R) * 0.05, 0.05, 0.95) : H.clamp(0.14 + H.gauss(R) * 0.035, 0.02, 0.95), c: i < 44 ? 0 : 1, d: R() });
      return { v: v };
    },
    draw: function (g, s, t) {
      var th = g.th, x0 = 60, w = 450, y0 = 215, cl = seg(t, 6.4, 7.2);
      g.axes(x0, y0, w, 165, "", "", seg(t, 0, 0.6));
      g.text("VAF", x0 + w + 8, y0, { col: th.muted, size: 10, a: seg(t, 0, 0.6) });
      [0, 0.25, 0.5, 0.75, 1].forEach(function (u) { g.text(String(u), x0 + u * w, y0 + 12, { col: th.muted, size: 9, align: "center", a: seg(t, 0, 0.6) }); });
      var bins = {};
      s.v.forEach(function (p) {
        var k = Math.round(p.x * 60), n = bins[k] = (bins[k] || 0) + 1;
        var land = seg(t, 0.3 + p.d * 2.4, 0.9 + p.d * 2.4);
        g.dot(x0 + p.x * w, lerp(40, y0 - 5 - (n - 1) * 8, land), 3.4, cl > 0 ? H.mix(th.muted, PAL[p.c ? 4 : 0], cl) : th.muted, land * 0.85);
      });
      var dn = seg(t, 3.2, 5, false), pts = [];
      for (var u = 0; u <= dn; u += 0.005) {
        var d = 0;
        s.v.forEach(function (p) { d += Math.exp(-Math.pow((u - p.x) / 0.03, 2) / 2); });
        pts.push([x0 + u * w, y0 - 70 - d * 3.2]);
      }
      g.poly(pts, th.ink, 1.8, 0.8);
      g.t("main clone", "主克隆", x0 + 0.4 * w, 52, { col: PAL[0], size: 10.5, bold: true, align: "center", a: cl });
      g.t("subclone", "亚克隆", x0 + 0.14 * w, 78, { col: PAL[4], size: 10.5, bold: true, align: "center", a: cl });
      g.note("The clonal peak sits near purity / 2, not at 0.5", "克隆峰位于纯度 / 2 附近，而不是 0.5", cl);
    }
  });
  // ---- variant filters ----------------------------------------------------------
  X.register("wes_filter", {
    period: 10,
    stages: [[0, 2.6, "Every call from the caller enters", "检测工具给出的每个突变都进入过滤"],
             [2.6, 6.6, "Each filter drops what it cannot trust", "每一道过滤去掉不可信的突变"],
             [6.6, 10, "Report what each filter removed", "报告每一道过滤去掉了多少"]],
    still: 9,
    init: function (R) {
      var fails = [0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3], calls = [];
      for (var i = 0; i < 40; i++) {
        calls.push({ y: 74 + R() * 92, d: R() * 1.6, vc: vcPick(R),
                     fail: i < fails.length ? fails[i] : -1, drop: 20 + R() * 30 });
      }
      calls.sort(function (a, b) { return a.d - b.d; });
      return { calls: calls, gates: [170, 260, 350, 440] };
    },
    draw: function (g, s, t) {
      var th = g.th, G = s.gates, lab = [["PASS", "PASS"], ["depth ≥ 10", "深度 ≥ 10"], ["VAF ≥ 5%", "VAF ≥ 5%"],
                                         ["gnomAD ≤ 0.1%", "gnomAD ≤ 0.1%"]];
      var summ = seg(t, 6.6, 7.4);
      G.forEach(function (x, k) {
        var a = seg(t, 0.3 + k * 0.25, 0.9 + k * 0.25) * (1 - summ * 0.6);
        g.line(x, 62, x, 176, th.accent, 2, a, [4, 4]);
        g.text(L(lab[k][0], lab[k][1]), x, 52, { col: th.accent2, size: 10, bold: true, align: "center", a: a });
      });
      var removed = [0, 0, 0, 0], kept = 0;
      s.calls.forEach(function (c) {
        var p = seg(t, 0.6 + c.d, 5.6 + c.d, false), x = lerp(40, 530, p);
        if (c.fail >= 0 && x >= G[c.fail]) {
          var q = H.clamp((x - G[c.fail]) / 40, 0, 1);
          removed[c.fail]++;
          g.dot(G[c.fail] + 6, c.y + q * c.drop, 3.2, H.FLAG, (1 - q * 0.75) * (1 - summ));
          return;
        }
        if (p >= 1) kept++;
        g.dot(x, c.y, 3.2, VC[c.vc], 0.9 * (1 - summ));
      });
      G.forEach(function (x, k) {
        if (removed[k]) g.text("−" + removed[k], x + 12, 196, { col: H.FLAG, size: 11, bold: true, a: 1 - summ });
      });
      if (summ > 0) {
        var left = [40, 37, 33, 30, 28], xs = [80, 215, 305, 395, 485], base = 236;
        left.forEach(function (n, k) {
          var hh = n / 40 * 150 * seg(t, 6.8 + k * 0.15, 7.6 + k * 0.15);
          g.rect(xs[k] - 16, base - hh, 32, hh, k === 4 ? th.ok : th.accent, 0.85, 3);
          g.text(String(n), xs[k], base - hh - 9, { col: th.text, size: 10.5, bold: true, align: "center", a: summ });
        });
        g.t("in", "输入", xs[0], 249, { col: th.muted, size: 9.5, align: "center", a: summ });
        g.t("kept", "保留", xs[4], 249, { col: th.ok, size: 9.5, bold: true, align: "center", a: summ });
      }
    }
  });

  // ---- TMB vs outcome -----------------------------------------------------------
  X.register("wes_tmbclin", {
    period: 10,
    stages: [[0, 3, "Each patient's TMB, on a log scale", "每位患者的 TMB，对数坐标"],
             [3, 6.4, "Split by outcome: does TMB shift?", "按结局分开：TMB 是否偏移？"],
             [6.4, 10, "One number per doubling, plus the ROC curve", "以每翻一倍衡量，再看 ROC 曲线"]],
    still: 9,
    init: function (R) {
      var pts = [];
      for (var i = 0; i < 46; i++) {
        var resp = R() < 0.4, v = H.gauss(R) * 0.75 + (resp ? 0.5 : -0.1);
        pts.push({ v: v, resp: resp, j: H.gauss(R) });
      }
      var pos = pts.filter(function (p) { return p.resp; }), neg = pts.filter(function (p) { return !p.resp; });
      var thr = pts.map(function (p) { return p.v; }).sort(function (a, b) { return b - a; }), roc = [[0, 0]];
      thr.forEach(function (c) {
        roc.push([neg.filter(function (p) { return p.v >= c; }).length / neg.length,
                  pos.filter(function (p) { return p.v >= c; }).length / pos.length]);
      });
      return { pts: pts, roc: roc };
    },
    draw: function (g, s, t) {
      var th = g.th, x0 = 50, x1 = 330, sp = seg(t, 3.2, 4.4), roc = seg(t, 6.6, 8.6, false);
      var xs = function (v) { return lerp(x0, x1, H.clamp((v + 2) / 4, 0, 1)); };
      var a0 = seg(t, 0.2, 1.2);
      g.line(x0, 214, x1, 214, th.muted, 1, a0 * 0.8);
      [["0.1", -1.6], ["1", 0], ["10", 1.6]].forEach(function (tk) {
        g.text(tk[0], lerp(x0, x1, (tk[1] + 2) / 4), 226, { col: th.muted, size: 9.5, align: "center", a: a0 });
      });
      g.t("TMB (mut/Mb, log)", "TMB（mut/Mb，对数）", (x0 + x1) / 2, 242, { col: th.muted, size: 10, align: "center", a: a0 });
      var ym = [0, 0], nm = [0, 0];
      s.pts.forEach(function (p, i) {
        var y = lerp(140 + p.j * 18, (p.resp ? 92 : 172) + p.j * 9, sp);
        var col = sp > 0 ? H.mix("#8a96a3", p.resp ? th.accent : "#8a96a3", sp) : th.accent;
        g.dot(xs(p.v), y, 3.3, col, seg(t, 0.3 + i * 0.03, 0.9 + i * 0.03) * 0.85);
        ym[p.resp ? 1 : 0] += p.v; nm[p.resp ? 1 : 0]++;
      });
      if (sp > 0) {
        g.t("responders", "应答", x0, 66, { col: th.accent2, size: 10, bold: true, a: sp });
        g.t("others", "其他", x0, 196, { col: th.muted, size: 10, bold: true, a: sp });
        [[1, 92], [0, 172]].forEach(function (r) {
          var mx = xs(ym[r[0]] / nm[r[0]]);
          g.line(mx, r[1] - 22, mx, r[1] + 22, r[0] ? th.accent : th.muted, 2.2, seg(t, 4.6, 5.4));
        });
      }
      if (roc > 0) {
        var rx = 380, ry = 222, rw = 150, rh = 150, pts = [];
        g.axes(rx, ry, rw, rh, L("1 − specificity", "1 − 特异度"), L("sensitivity", "灵敏度"), seg(t, 6.4, 7));
        g.line(rx, ry, rx + rw, ry - rh, th.faint, 1.2, seg(t, 6.4, 7), [4, 4]);
        var k = Math.max(1, Math.round(roc * (s.roc.length - 1)));
        for (var i = 0; i <= k; i++) pts.push([rx + s.roc[i][0] * rw, ry - s.roc[i][1] * rh]);
        g.poly(pts, th.accent, 2.2, 1);
        g.t("OR per doubling, AUC", "每翻一倍的 OR、AUC", rx + rw, 48, { col: th.accent2, size: 10.5, bold: true, align: "right", a: seg(t, 8.4, 9) });
      }
    }
  });

  // ---- WES report -----------------------------------------------------------------
  X.register("wes_report", {
    period: 8.5,
    stages: [[0, 2.8, "Every WES step you ran is logged", "你运行的每个 WES 步骤都被记录"],
             [2.8, 5.8, "Parameters and code are written up", "参数与代码被整理成文"],
             [5.8, 8.5, "Download the report and the R script", "下载报告与 R 脚本"]],
    still: 7.2,
    init: function (R) {
      var lines = [];
      for (var i = 0; i < 9; i++) lines.push(0.5 + R() * 0.45);
      return { lines: lines };
    },
    draw: function (g, s, t) {
      var th = g.th, dx = 210, dy = 40, dw = 170, dh = 196;
      ["Import", "Filters", "Oncoplot", "TMB", "Survival"].forEach(function (st, i) {
        var p = seg(t, 0.3 + i * 0.4, 1.1 + i * 0.4);
        g.chip(st, lerp(80, dx + 20, seg(t, 2.6 + i * 0.12, 3.2 + i * 0.12)), 70 + i * 30,
               PAL[i], p * (1 - seg(t, 2.8 + i * 0.12, 3.3 + i * 0.12)), { size: 10, align: "left" });
      });
      var da = seg(t, 2.6, 3.2);
      g.rect(dx, dy, dw, dh, th.card, da, 6);
      g.box(dx, dy, dw, dh, th.border, 1.2, da, 6);
      g.rect(dx + 14, dy + 14, 90, 8, th.ink, da * 0.8, 3);
      [["Methods", "方法"], ["Parameters", "参数"], ["Code", "代码"]].forEach(function (hd, k) {
        var y = dy + 42 + k * 50, p = seg(t, 3.2 + k * 0.7, 3.8 + k * 0.7);
        g.text(L(hd[0], hd[1]), dx + 14, y, { col: th.accent2, size: 10, bold: true, a: p });
        for (var j = 0; j < 2; j++) {
          var w = (dw - 28) * s.lines[k * 3 + j] * seg(t, 3.4 + k * 0.7 + j * 0.2, 4 + k * 0.7 + j * 0.2, false);
          g.rect(dx + 14, y + 12 + j * 10, w, 4, k === 2 ? th.ok : th.faint, 0.9, 2);
        }
      });
      var dl = seg(t, 6, 6.8);
      g.chip("report.html", 460, 110, th.accent, dl, { size: 10.5 });
      g.chip("analysis.R", 460, 145, th.ok, dl, { size: 10.5 });
      g.arrow(dx + dw + 8, 128, dx + dw + 8 + 30 * dl, 128, th.accent, 2, dl);
    }
  });
})();
