/* OmicOne step animations — single-cell scenes (one per step key).
 * Drawing space: 560 x 260; captions occupy the top 30 px. See explain.js. */
(function () {
  "use strict";
  var X = window.OmicOneExplain;
  if (!X) return;
  var H = X.helpers, seg = H.seg, lerp = H.lerp, PAL = H.PAL, L = H.L;

  function clusters(R, spec) {                 // spec: [[cx, cy, n, sd], ...]
    var pts = [];
    spec.forEach(function (c, k) {
      H.blob(R, c[2], c[0], c[1], c[3], c[4]).forEach(function (p, i) {
        p.k = k; p.i = pts.length; p.r = 2.2 + R() * 1.2; p.d = R(); pts.push(p);
      });
    });
    return pts;
  }
  function nearest(pts, k) {                   // k nearest neighbour index lists
    return pts.map(function (p, i) {
      var d = pts.map(function (q, j) {
        return { j: j, d: (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) };
      });
      d.sort(function (a, b) { return a.d - b.d; });
      return d.slice(1, k + 1).map(function (o) { return o.j; });
    });
  }
  function kmCurve(R, n, rate, horizon) {      // simple exponential KM step points
    var times = [];
    for (var i = 0; i < n; i++) times.push(Math.min(horizon, -Math.log(1 - R()) / rate));
    times.sort(function (a, b) { return a - b; });
    var s = 1, pts = [[0, 1]], cens = [];
    times.forEach(function (tt, i) {
      var atRisk = n - i;
      if (tt >= horizon || R() < 0.25) { cens.push([tt, s]); return; }
      pts.push([tt, s]); s = s * (1 - 1 / atRisk); pts.push([tt, s]);
    });
    pts.push([horizon, s]);
    return { pts: pts, cens: cens };
  }
  function drawKM(g, km, x0, y0, w, h, horizon, col, p) {
    var out = [], lim = horizon * p;
    for (var i = 0; i < km.pts.length; i++) {
      var q = km.pts[i];
      if (q[0] > lim) { out.push([x0 + lim / horizon * w, y0 - (out.length ? (y0 - out[out.length - 1][1]) / h : 1) * h]); break; }
      out.push([x0 + q[0] / horizon * w, y0 - q[1] * h]);
    }
    g.poly(out, col, 2.2, 1);
    km.cens.forEach(function (c) {
      if (c[0] > lim) return;
      var cx = x0 + c[0] / horizon * w, cy = y0 - c[1] * h;
      g.line(cx, cy - 3.5, cx, cy + 3.5, col, 1.4, 0.9);
    });
  }

  // ---- 1. import -------------------------------------------------------------
  X.register("import", {
    period: 8,
    stages: [[0, 2.2, "Read the count file (.h5, .rds or a table)", "读取计数文件（.h5、.rds 或表格）"],
             [2.2, 5, "Genes × cells count matrix: most entries are 0", "基因 × 细胞的计数矩阵：大多数元素为 0"],
             [5, 8, "Each column becomes one cell", "每一列就是一个细胞"]],
    init: function (R) {
      var m = [];
      for (var i = 0; i < 8; i++) { m.push([]); for (var j = 0; j < 12; j++) m[i].push(R() < 0.33 ? 1 + Math.floor(R() * 9) : 0); }
      var cell = [];
      for (j = 0; j < 12; j++) cell.push({ k: j % 3, tx: 455 + (j % 3) * 30 + H.gauss(R) * 9, ty: 95 + Math.floor(j / 3) * 28 + H.gauss(R) * 6 });
      return { m: m, cell: cell };
    },
    draw: function (g, s, t) {
      var th = g.th, a1 = seg(t, 0.1, 0.9);
      // file icon
      var fx = 34, fy = 78, fw = 70, fh = 92;
      g.rect(fx, fy + (1 - a1) * 10, fw, fh, th.panel, a1, 8);
      g.box(fx, fy + (1 - a1) * 10, fw, fh, th.muted, 1.2, a1, 8);
      g.fillPoly([[fx + fw - 18, fy + (1 - a1) * 10], [fx + fw, fy + 18 + (1 - a1) * 10], [fx + fw - 18, fy + 18 + (1 - a1) * 10]], th.faint, a1);
      for (var r = 0; r < 4; r++) g.line(fx + 12, fy + 34 + r * 12, fx + fw - 12, fy + 34 + r * 12, th.muted, 1.4, a1 * 0.6);
      g.text("counts", fx + fw / 2, fy + fh + 14, { col: th.muted, size: 10, align: "center", a: a1 });
      // arrow
      var a2 = seg(t, 0.8, 1.8);
      g.arrow(116, 124, 116 + 44 * a2, 124, th.accent, 2, a2);
      // matrix
      var x0 = 222, y0 = 60, cs = 15;
      for (var i = 0; i < 8; i++) {
        g.text("GENE" + (i + 1), x0 - 6, y0 + i * cs + cs / 2, { col: th.muted, size: 9, align: "right", a: seg(t, 1.9, 2.4) });
        for (var j = 0; j < 12; j++) {
          var aj = seg(t, 2 + j * 0.18, 2.4 + j * 0.18);
          var v = s.m[i][j], x = x0 + j * cs, y = y0 + i * cs;
          if (v) g.rect(x + 1, y + 1, cs - 2, cs - 2, th.accent, aj * (0.25 + v / 12), 2);
          else g.box(x + 1.5, y + 1.5, cs - 3, cs - 3, th.faint, 1, aj, 2);
        }
      }
      g.text(L("cells →", "细胞 →"), x0 + 6 * cs, y0 - 9, { col: th.muted, size: 10, align: "center", a: seg(t, 2, 2.6) });
      // columns fly out as cells
      for (j = 0; j < 12; j++) {
        var p = seg(t, 5 + j * 0.16, 5.9 + j * 0.16);
        if (p <= 0) continue;
        var cx = x0 + j * cs + cs / 2;
        if (p < 1) g.box(x0 + j * cs, y0, cs, 8 * cs, PAL[s.cell[j].k], 1.4, 1 - p, 2);
        var c = s.cell[j];
        g.dot(lerp(cx, c.tx, p), lerp(y0 + 4 * cs, c.ty, p), lerp(2, 5.5, p), PAL[c.k], 0.85);
      }
    }
  });

  // ---- 2. QC -------------------------------------------------------------------
  X.register("qc", {
    period: 8.5,
    stages: [[0, 2.6, "Each dot is a cell: depth vs mitochondrial %", "每个点是一个细胞：测序深度与线粒体比例"],
             [2.6, 5.2, "MAD thresholds adapt to this dataset", "MAD 阈值随本数据集自适应"],
             [5.2, 8.5, "Outliers are flagged and removed", "离群细胞被标记并移除"]],
    still: 6.2,
    init: function (R) {
      var pts = [];
      for (var i = 0; i < 130; i++) pts.push({ x: H.clamp(0.55 + H.gauss(R) * 0.11, 0.25, 0.85), y: Math.abs(0.16 + H.gauss(R) * 0.06) });
      for (i = 0; i < 9; i++) pts.push({ x: 0.35 + R() * 0.4, y: 0.55 + R() * 0.38 });
      for (i = 0; i < 6; i++) pts.push({ x: 0.03 + R() * 0.13, y: 0.1 + R() * 0.4 });
      for (i = 0; i < 3; i++) pts.push({ x: 0.9 + R() * 0.07, y: 0.12 + R() * 0.15 });
      pts.forEach(function (p, k) { p.bad = p.x < 0.2 || p.x > 0.88 || p.y > 0.45; p.d = R(); });
      return { pts: pts, kept: pts.filter(function (p) { return !p.bad; }).length };
    },
    draw: function (g, s, t) {
      var th = g.th, x0 = 70, y0 = 222, w = 420, h = 172;
      g.axes(x0, y0, w, h, L("UMI count (log)", "UMI 数（对数）"), L("mito %", "线粒体 %"), seg(t, 0, 0.6));
      var lo = x0 + 0.2 * w, hi = x0 + 0.88 * w, mt = y0 - 0.45 * h, ab = seg(t, 2.7, 3.8);
      if (ab > 0) {
        g.rect(lo, mt, (hi - lo) * ab, y0 - mt, th.ok, 0.07);
        g.line(lo, y0, lo, y0 - h * ab, th.accent, 1.4, 0.9, [5, 4]);
        g.line(hi, y0, hi, y0 - h * ab, th.accent, 1.4, 0.9, [5, 4]);
        g.line(x0, mt, x0 + w * ab, mt, th.accent, 1.4, 0.9, [5, 4]);
        g.t("median ± k·MAD", "中位数 ± k·MAD", hi - 4, y0 - h + 8, { col: th.accent2, size: 10, align: "right", a: ab });
      }
      var flag = seg(t, 4.2, 5.0), gone = seg(t, 5.6, 6.8);
      s.pts.forEach(function (p) {
        var a = seg(t, 0.2 + p.d * 1.8, 0.5 + p.d * 1.8);
        var x = x0 + p.x * w, y = y0 - p.y * h, r = 2.6;
        var col = H.KEPT;
        if (p.bad) {
          col = H.mix(H.KEPT, H.FLAG, flag);
          a *= 1 - gone;
          y += gone * 18;
          r *= 1 + 0.4 * flag;
        }
        g.dot(x, y, r, col, a * 0.8);
      });
      var c = seg(t, 5.8, 6.6);
      g.t("kept " + s.kept + " / " + s.pts.length, "保留 " + s.kept + " / " + s.pts.length,
          x0 + w, 42, { col: th.ok, size: 12, bold: true, align: "right", a: c });
    }
  });

  // ---- 3. doublets ----------------------------------------------------------
  X.register("doublet", {
    period: 9,
    stages: [[0, 2.6, "Two cells caught in one droplet look like a hybrid", "两个细胞落进同一液滴，看起来像混合体"],
             [2.6, 5.8, "Simulate doublets by mixing random cell pairs", "随机配对细胞，模拟双细胞"],
             [5.8, 9, "Cells resembling the simulations are flagged", "与模拟双细胞相似的细胞被标记"]],
    still: 7.4,
    init: function (R) {
      var A = H.blob(R, 45, 160, 145, 24, 22), B = H.blob(R, 45, 400, 145, 24, 22);
      var real = H.blob(R, 6, 280, 145, 14, 16), pairs = [];
      for (var i = 0; i < 9; i++) pairs.push([A[Math.floor(R() * A.length)], B[Math.floor(R() * B.length)], R()]);
      return { A: A, B: B, real: real, pairs: pairs };
    },
    draw: function (g, s, t) {
      var th = g.th, cA = PAL[0], cB = PAL[2];
      var a = seg(t, 0.1, 1.2);
      s.A.forEach(function (p) { g.dot(p.x, p.y, 2.8, cA, a * 0.75); });
      s.B.forEach(function (p) { g.dot(p.x, p.y, 2.8, cB, a * 0.75); });
      g.t("cell type A", "细胞类型 A", 160, 214, { col: cA, size: 10.5, align: "center", a: a });
      g.t("cell type B", "细胞类型 B", 400, 214, { col: cB, size: 10.5, align: "center", a: a });
      var flag = seg(t, 6, 6.8);
      s.real.forEach(function (p, i) {
        var ar = seg(t, 0.8 + i * 0.12, 1.4 + i * 0.12);
        g.box(p.x - 7.5, p.y - 6, 15, 12, th.muted, 1, ar * 0.6 * (1 - flag), 6);
        g.dot(p.x - 2.4, p.y, 2.7, cA, ar * 0.9);
        g.dot(p.x + 2.4, p.y, 2.7, cB, ar * 0.9);
        if (flag > 0) g.ring(p.x, p.y, 7 + 2 * Math.sin(t * 6 + i), H.FLAG, 1.8, flag);
      });
      s.pairs.forEach(function (q, i) {
        var p = seg(t, 2.8 + i * 0.22, 3.8 + i * 0.22), out = 1 - seg(t, 6, 6.8);
        if (p <= 0) return;
        var mx = (q[0].x + q[1].x) / 2, my = (q[0].y + q[1].y) / 2 + (q[2] - 0.5) * 30;
        g.line(q[0].x, q[0].y, lerp(q[0].x, mx, p), lerp(q[0].y, my, p), cA, 0.8, 0.45 * out);
        g.line(q[1].x, q[1].y, lerp(q[1].x, mx, p), lerp(q[1].y, my, p), cB, 0.8, 0.45 * out);
        if (p >= 1) g.fillPoly([[mx, my - 4.5], [mx + 4.5, my], [mx, my + 4.5], [mx - 4.5, my]], "#9c6ade", 0.9 * out);
      });
      g.t("◆ simulated doublet", "◆ 模拟双细胞", 280, 52, { col: "#9c6ade", size: 10.5, align: "center", a: seg(t, 3, 3.6) * (1 - seg(t, 6, 6.8)) });
      g.t("flagged: high doublet score", "被标记：双细胞得分高", 280, 52, { col: H.FLAG, size: 10.5, align: "center", bold: true, a: flag });
      g.note("Run per sample: doublets only form within one capture", "按样本运行：双细胞只会在同一次捕获内形成", seg(t, 6.2, 7));
    }
  });

  // ---- 4. normalize ---------------------------------------------------------
  X.register("normalize", {
    period: 8.5,
    stages: [[0, 2.6, "Raw counts depend on sequencing depth", "原始计数受测序深度影响"],
             [2.6, 5.4, "Scale every cell to the same total (10,000)", "把每个细胞缩放到相同总量（10,000）"],
             [5.4, 8.5, "Then log-transform: log1p(x)", "再取对数：log1p(x)"]],
    init: function (R) {
      var tot = [150, 62, 110, 172, 86, 128], cells = [];
      tot.forEach(function (T) {
        var f = [0.42, 0.28, 0.19, 0.11].map(function (x) { return x * (0.85 + R() * 0.3); });
        var sum = f.reduce(function (a, b) { return a + b; }, 0);
        cells.push({ T: T, f: f.map(function (x) { return x / sum; }) });
      });
      return { cells: cells };
    },
    draw: function (g, s, t) {
      var th = g.th, base = 222, eq = 120;
      g.axes(70, base, 430, 180, "",
             t < 5.4 ? L("counts", "计数") : "log1p(CP10k)", 0.8);
      var grow = seg(t, 0.2, 1.8), norm = seg(t, 2.8, 4.3), lg = seg(t, 5.6, 7);
      if (norm > 0) {
        g.line(70, base - eq, 500, base - eq, th.accent, 1.2, norm * (1 - lg), [5, 4]);
        g.text("10,000", 502, base - eq, { col: th.accent2, size: 10, a: norm * (1 - lg) });
      }
      s.cells.forEach(function (c, i) {
        var x = 100 + i * 66, total = lerp(c.T, eq, norm) * grow, y = base;
        var logs = c.f.map(function (f) { return Math.log(1 + f * 30); });
        var lsum = logs.reduce(function (a, b) { return a + b; }, 0);
        c.f.forEach(function (f, k) {
          var hLin = total * f, hLog = 92 * logs[k] / lsum * grow;
          var hh = lerp(hLin, hLog, lg);
          g.rect(x, y - hh, 36, hh - 1, PAL[k], 0.85, 2);
          y -= hh;
        });
        g.text(L("cell ", "细胞 ") + (i + 1), x + 18, base + 10, { col: th.muted, size: 9.5, align: "center", a: grow });
      });
      ["GENE A", "GENE B", "GENE C", "GENE D"].forEach(function (nm, k) {
        g.rect(452, 40 + k * 14, 9, 9, PAL[k], seg(t, 0.4, 1.2));
        g.text(nm, 466, 45 + k * 14, { col: th.muted, size: 9.5, a: seg(t, 0.4, 1.2) });
      });
    }
  });

  // ---- 5. features / PCA ----------------------------------------------------
  X.register("reduce", {
    period: 10,
    stages: [[0, 3.2, "Cells live in thousands of gene dimensions", "细胞处在成千上万个基因维度中"],
             [3.2, 6.2, "PC1 follows the direction of largest variation", "PC1 沿着变异最大的方向"],
             [6.2, 10, "Keep the top PCs: a compact summary", "保留前几个主成分：紧凑的概要"]],
    still: 8.5,
    init: function (R) {
      var pts = [];
      for (var i = 0; i < 140; i++) {
        var k = i % 3;
        pts.push({ u: (k - 1) * 2.1 + H.gauss(R) * 0.75, v: H.gauss(R) * 0.9, w: H.gauss(R) * 0.35, k: k });
      }
      // fixed rotation from PC space into "gene space"
      var a = 0.7, b = -0.5;
      var R1 = [[Math.cos(a), 0, Math.sin(a)], [0, 1, 0], [-Math.sin(a), 0, Math.cos(a)]];
      var R2 = [[1, 0, 0], [0, Math.cos(b), -Math.sin(b)], [0, Math.sin(b), Math.cos(b)]];
      var M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      for (var r = 0; r < 3; r++) for (var c = 0; c < 3; c++) for (var q = 0; q < 3; q++) M[r][c] += R2[r][q] * R1[q][c];
      return { pts: pts, M: M };
    },
    draw: function (g, s, t) {
      var th = g.th, flat = seg(t, 6.3, 8), spin = (t < 6.3 ? t : 6.3) * 0.55 * (1 - flat);
      var cx = 230, cy = 145, sc = 26;
      function rot(u, v, w) {                         // PC space -> screen
        var M = s.M, p = [u, v, w * (1 - flat)];
        var x = 0, y = 0, z = 0;
        for (var c = 0; c < 3; c++) {
          var m = H.lerp(M[0][c], c === 0 ? 1 : 0, flat), n = H.lerp(M[1][c], c === 1 ? 1 : 0, flat),
              o = H.lerp(M[2][c], c === 2 ? 1 : 0, flat);
          x += m * p[c]; y += n * p[c]; z += o * p[c];
        }
        var xs = x * Math.cos(spin) + z * Math.sin(spin), zs = -x * Math.sin(spin) + z * Math.cos(spin);
        var f = 9 / (9 + zs);
        return [cx + xs * sc * f, cy - y * sc * f, zs];
      }
      var o = rot(0, 0, 0), ax = seg(t, 0.2, 1.2) * (1 - flat) * (1 - 0.75 * seg(t, 3.4, 4.4));
      [[3.6, 0, 0], [0, 3, 0], [0, 0, 3]].forEach(function (d, i) {
        var e = rot(d[0], d[1], d[2]);
        g.line(o[0], o[1], e[0], e[1], th.faint, 1, ax);
        g.text(["gene 1", "gene 2", "gene 3"][i], e[0] + 4, e[1], { col: th.muted, size: 9, a: ax });
      });
      var proj = s.pts.map(function (p) { var q = rot(p.u, p.v, p.w); return { x: q[0], y: q[1], z: q[2], k: p.k }; });
      proj.sort(function (a, b) { return b.z - a.z; });
      proj.forEach(function (p) { g.dot(p.x, p.y, 2.7 - p.z * 0.12, PAL[p.k], seg(t, 0.2, 1.4) * 0.8); });
      var pa = seg(t, 3.4, 4.4);
      if (pa > 0) {
        var e1 = rot(4.2 * pa, 0, 0), e2 = rot(0, 2.2 * seg(t, 4.2, 5), 0), n1 = rot(-4.2 * pa, 0, 0);
        g.arrow(n1[0], n1[1], e1[0], e1[1], th.accent, 2.4, 1);
        g.text("PC1", e1[0] + 6, e1[1], { col: th.accent, size: 11, bold: true, a: pa });
        g.arrow(o[0], o[1], e2[0], e2[1], th.ok, 2.2, seg(t, 4.2, 5));
        g.text("PC2", e2[0] + 6, e2[1] - 4, { col: th.ok, size: 11, bold: true, a: seg(t, 4.2, 5) });
      }
      // elbow plot
      var ea = seg(t, 6.8, 7.8), var_ = [1, 0.52, 0.3, 0.17, 0.11, 0.08, 0.06, 0.05];
      if (ea > 0) {
        g.axes(410, 210, 128, 130, "PC", L("variance", "方差"), ea);
        var_.forEach(function (v, i) {
          g.rect(416 + i * 15, 210 - v * 118 * ea, 11, v * 118 * ea, i < 3 ? th.accent : th.faint, 0.9, 2);
        });
        g.t("keep", "保留", 438, 70, { col: th.accent2, size: 10, align: "center", a: ea });
      }
    }
  });

  // ---- 6. integrate ---------------------------------------------------------
  X.register("integrate", {
    period: 9.5,
    stages: [[0, 3, "Colour = sample: one cell type splits by batch", "颜色 = 样本：同一细胞类型因批次分开"],
             [3, 6.2, "Harmony aligns the samples", "Harmony 对齐各个样本"],
             [6.2, 9.5, "Colour = cell type: the biology is kept", "颜色 = 细胞类型：生物学差异被保留"]],
    still: 8.2,
    init: function (R) {
      var types = [[160, 110], [300, 185], [420, 105]], pts = [];
      types.forEach(function (c, k) {
        for (var b = 0; b < 2; b++) {
          H.blob(R, 26, c[0], c[1], 15, 13).forEach(function (p) { p.k = k; p.b = b; pts.push(p); });
        }
      });
      return { pts: pts };
    },
    draw: function (g, s, t) {
      var th = g.th, al = seg(t, 3.2, 5.6), rc = seg(t, 6.4, 7.4), BC = ["#2f81c7", "#f4a261"];
      var TC = [PAL[2], PAL[3], PAL[1]];
      s.pts.forEach(function (p) {
        var off = (p.b ? 1 : -1) * (1 - al);
        var x = p.x + off * 44, y = p.y + off * 20;
        g.dot(x, y, 2.8, H.mix(BC[p.b], TC[p.k], rc), seg(t, 0.1, 1.2) * 0.8);
      });
      var la = seg(t, 0.4, 1.2);
      if (rc < 1) {
        g.dot(400, 230, 4, BC[0], la * (1 - rc)); g.t("sample 1", "样本 1", 408, 230, { col: th.muted, size: 10, a: la * (1 - rc) });
        g.dot(470, 230, 4, BC[1], la * (1 - rc)); g.t("sample 2", "样本 2", 478, 230, { col: th.muted, size: 10, a: la * (1 - rc) });
      }
      [["T cells", "T 细胞"], ["B cells", "B 细胞"], ["Mono", "单核"]].forEach(function (n, k) {
        g.dot(330 + k * 72, 230, 4, TC[k], rc); g.t(n[0], n[1], 338 + k * 72, 230, { col: th.muted, size: 10, a: rc });
      });
    }
  });

  // ---- 7. cluster -----------------------------------------------------------
  X.register("cluster", {
    period: 9.5,
    stages: [[0, 3, "Connect each cell to its nearest neighbours", "把每个细胞连到最近的邻居"],
             [3, 6.2, "Densely linked communities become clusters", "连接紧密的群落成为聚类"],
             [6.2, 9.5, "Higher resolution splits clusters further", "分辨率越高，聚类分得越细"]],
    still: 8.4,
    init: function (R) {
      var pts = clusters(R, [[140, 115, 26, 17], [245, 185, 26, 17], [335, 100, 26, 17], [430, 150, 18, 13], [478, 182, 16, 12]]);
      pts.forEach(function (p) { if (p.k === 4) { p.sub = 1; p.k = 3; } });
      var nn = nearest(pts, 3), edges = [];
      nn.forEach(function (l, i) { l.forEach(function (j) { if (i < j || nn[j].indexOf(i) < 0) edges.push([i, j, R()]); }); });
      return { pts: pts, edges: edges };
    },
    draw: function (g, s, t) {
      var th = g.th, col = seg(t, 3.2, 5.2), split = seg(t, 6.5, 7.5);
      s.edges.forEach(function (e) {
        var a = seg(t, 0.4 + e[2] * 2, 0.8 + e[2] * 2);
        var p = s.pts[e[0]], q = s.pts[e[1]], same = p.k === q.k;
        g.line(p.x, p.y, q.x, q.y, same ? th.muted : th.faint, 0.8, a * (same ? 0.55 : 0.55 * (1 - col)));
      });
      s.pts.forEach(function (p) {
        var c = PAL[p.k];
        if (p.sub) c = H.mix(PAL[3], PAL[4], split);
        g.dot(p.x, p.y, 3, H.mix(th.dark ? "#8b98a5" : "#9aa5b1", c, seg(t, 3.2 + p.d * 1.6, 3.8 + p.d * 1.6)), 0.9);
      });
      g.text("resolution " + (split > 0.5 ? "1.0" : "0.5"), 540, 42,
             { col: th.accent2, size: 11, bold: true, align: "right", a: seg(t, 6.2, 6.8) });
    }
  });

  // ---- 8. embed (UMAP) ------------------------------------------------------
  X.register("embed", {
    period: 9.5,
    stages: [[0, 2.8, "Start: cells scattered, neighbours known from the graph", "起点：细胞散落，邻居关系来自图"],
             [2.8, 6.2, "Pull neighbours together, push the rest apart", "拉近邻居，推开其余细胞"],
             [6.2, 9.5, "A 2-D map of the neighbourhood structure", "一张保留邻域结构的二维图"]],
    still: 8.4,
    init: function (R) {
      var pts = clusters(R, [[130, 120, 24, 15], [250, 190, 24, 15], [360, 105, 24, 15], [460, 175, 24, 15]]);
      pts.forEach(function (p) { p.sx = 50 + R() * 460; p.sy = 45 + R() * 190; p.sw = (R() - 0.5) * 60; });
      return { pts: pts, nn: nearest(pts, 2) };
    },
    draw: function (g, s, t) {
      var th = g.th;
      var pos = s.pts.map(function (p) {
        var q = seg(t, 2.9 + p.d * 0.8, 5.4 + p.d * 0.8);
        return [lerp(p.sx, p.x, q) + Math.sin(q * Math.PI) * p.sw, lerp(p.sy, p.y, q) - Math.sin(q * Math.PI) * p.sw * 0.4];
      });
      var la = seg(t, 0.5, 1.5) * (1 - seg(t, 6.2, 7));
      s.nn.forEach(function (l, i) {
        if (i % 2) return;
        l.forEach(function (j) { g.line(pos[i][0], pos[i][1], pos[j][0], pos[j][1], PAL[s.pts[i].k], 0.8, la * 0.5); });
      });
      s.pts.forEach(function (p, i) { g.dot(pos[i][0], pos[i][1], 2.9, PAL[p.k], seg(t, 0.1, 1) * 0.85); });
      var lab = seg(t, 6.4, 7.2);
      [[130, 92], [250, 162], [360, 77], [460, 147]].forEach(function (c, k) {
        g.text(String(k), c[0], c[1], { col: PAL[k], size: 13, bold: true, align: "center", a: lab });
      });
      g.note("Distances between clusters are not quantitative", "簇与簇之间的距离不具定量意义", lab);
      g.text("UMAP 1", 290, 248, { col: th.muted, size: 9.5, align: "center", a: lab });
    }
  });

  // ---- 9. markers -----------------------------------------------------------
  X.register("markers", {
    period: 9.5,
    stages: [[0, 3, "Compare one cluster with all the other cells", "把一个聚类与其余所有细胞比较"],
             [3, 6.2, "High here, low elsewhere: a marker gene", "这里高、别处低：标志基因"],
             [6.2, 9.5, "Dot plot: size = % expressing, colour = mean", "点图：大小 = 表达比例，颜色 = 平均表达"]],
    still: 8.4,
    init: function (R) {
      return { pts: clusters(R, [[70, 110, 18, 12], [150, 175, 18, 12], [200, 95, 18, 12], [95, 200, 14, 10]]) };
    },
    draw: function (g, s, t) {
      var th = g.th, hi = 2, focus = seg(t, 0.8, 1.8);
      s.pts.forEach(function (p) {
        g.dot(p.x, p.y, 2.7, PAL[p.k], seg(t, 0.1, 0.9) * (p.k === hi ? 0.95 : lerp(0.85, 0.22, focus)));
      });
      g.ring(200, 95, 38, PAL[hi], 1.4, focus * 0.7);
      g.text(L("cluster 2", "簇 2"), 200, 145, { col: PAL[hi], size: 10.5, bold: true, align: "center", a: focus });
      g.t("vs rest", "对比其余", 120, 236, { col: th.muted, size: 10, align: "center", a: focus });
      var vio = seg(t, 3.2, 4.4) * (1 - seg(t, 6.2, 6.8));
      if (vio > 0) {
        g.axes(290, 222, 240, 170, "", "MS4A1", vio);
        [0.12, 0.18, 0.9, 0.1].forEach(function (m, k) {
          var cx = 320 + k * 58, top = 222 - m * 160 * vio, pts = [], pts2 = [];
          for (var y = 0; y <= 1; y += 0.1) {
            var wv = Math.sin(y * Math.PI) * (k === hi ? 18 : 12) * (0.6 + 0.4 * Math.cos((y - 0.5) * 3));
            pts.push([cx - wv, lerp(222, top, y)]); pts2.unshift([cx + wv, lerp(222, top, y)]);
          }
          g.fillPoly(pts.concat(pts2), PAL[k], 0.6 * vio);
          g.text(String(k), cx, 232, { col: th.muted, size: 10, align: "center", a: vio });
        });
      }
      var dp = seg(t, 6.5, 7.4);
      if (dp > 0) {
        var genes = ["CD3E", "LYZ", "MS4A1", "NKG7"];
        genes.forEach(function (gn, r) {
          g.text(gn, 320, 70 + r * 38, { col: th.muted, size: 10, align: "right", a: dp });
          for (var c = 0; c < 4; c++) {
            var on = (r === 0 && c === 0) || (r === 1 && c === 1) || (r === 2 && c === hi) || (r === 3 && c === 3);
            g.dot(355 + c * 50, 70 + r * 38, (on ? 12 : 3.5) * dp, on ? th.accent : th.faint, on ? 0.9 : 0.8);
          }
        });
        for (var c2 = 0; c2 < 4; c2++) g.text(String(c2), 355 + c2 * 50, 222, { col: th.muted, size: 10, align: "center", a: dp });
      }
    }
  });

  // ---- 10. annotate ---------------------------------------------------------
  X.register("annotate", {
    period: 9.5,
    stages: [[0, 2.8, "So far, clusters are only numbers", "此时聚类只是编号"],
             [2.8, 6, "Known marker genes point to identities", "已知的标志基因指向细胞身份"],
             [6, 9.5, "Name each cluster", "为每个聚类命名"]],
    still: 8.4,
    init: function (R) {
      return { pts: clusters(R, [[140, 110, 26, 16], [270, 185, 26, 16], [360, 95, 26, 16], [455, 175, 22, 14]]) };
    },
    draw: function (g, s, t) {
      var th = g.th, C = [[140, 110], [270, 185], [360, 95], [455, 175]];
      s.pts.forEach(function (p) { g.dot(p.x, p.y, 2.9, PAL[p.k], seg(t, 0.1, 1) * 0.8); });
      var gene = ["CD3E", "MS4A1", "LYZ", "NKG7"], name = [["T cells", "T 细胞"], ["B cells", "B 细胞"],
                  ["Monocytes", "单核细胞"], ["NK cells", "NK 细胞"]];
      var nm = seg(t, 6.2, 7);
      C.forEach(function (c, k) {
        g.text(String(k), c[0], c[1] - 40, { col: PAL[k], size: 14, bold: true, align: "center", a: seg(t, 0.4, 1) * (1 - nm) });
        var p = seg(t, 3 + k * 0.35, 3.8 + k * 0.35);
        if (p > 0) g.chip(gene[k], lerp(560, c[0] + 44, p), c[1] + 26, PAL[k], p * (1 - nm * 0.7), { size: 9.5 });
        g.text(L(name[k][0], name[k][1]), c[0], c[1] - 40, { col: PAL[k], size: 12.5, bold: true, align: "center", a: nm });
      });
    }
  });

  // ---- 11. enrichment (GSEA) ------------------------------------------------
  X.register("enrichment", {
    period: 9.5,
    stages: [[0, 2.8, "Rank all tested genes by fold change", "按倍数变化给所有被检验的基因排序"],
             [2.8, 6.4, "Walk down the list: step up at gene-set members", "沿列表向下：遇到通路基因就上升"],
             [6.4, 9.5, "The peak is the enrichment score", "峰值就是富集分数"]],
    still: 8.4,
    init: function (R) {
      var n = 90, hit = [];
      for (var i = 0; i < n; i++) hit.push(R() < (i < 30 ? 0.42 : 0.08));
      var nh = hit.filter(Boolean).length, es = [0], cur = 0, best = 0, bi = 0;
      hit.forEach(function (h, i) { cur += h ? 1 / nh : -1 / (n - nh); es.push(cur); if (cur > best) { best = cur; bi = i + 1; } });
      return { n: n, hit: hit, es: es, best: best, bi: bi };
    },
    draw: function (g, s, t) {
      var th = g.th, x0 = 60, x1 = 520, w = (x1 - x0) / s.n;
      var sa = seg(t, 0.1, 2);
      for (var i = 0; i < s.n; i++) {
        if (i / s.n > sa) break;
        var v = 1 - i / (s.n - 1);
        g.rect(x0 + i * w, 222, w + 0.4, 16, v > 0.5 ? H.mix("#f4f6f9", "#e4572e", (v - 0.5) * 2) : H.mix("#2f81c7", "#f4f6f9", v * 2), 0.95);
      }
      g.t("up in cluster", "在该簇上调", x0, 248, { col: th.muted, size: 9.5, a: sa });
      g.t("down", "下调", x1, 248, { col: th.muted, size: 9.5, align: "right", a: sa });
      var ta = seg(t, 1.2, 2.4);
      s.hit.forEach(function (h, i) { if (h) g.line(x0 + (i + 0.5) * w, 200, x0 + (i + 0.5) * w, 215, th.text, 1.2, ta); });
      var walk = seg(t, 3, 6, false), y0 = 150, sc = 105, pts = [];
      g.line(x0, y0, x1, y0, th.faint, 1, seg(t, 2.8, 3.2));
      g.text("ES", x0 - 8, 60, { col: th.muted, size: 10, align: "right", a: seg(t, 2.8, 3.2) });
      for (i = 0; i <= Math.round(walk * s.n); i++) pts.push([x0 + i * w, y0 - s.es[i] * sc]);
      g.poly(pts, PAL[1], 2.2, 1);
      if (walk > 0 && walk < 1) g.line(x0 + walk * s.n * w, 200, x0 + walk * s.n * w, 240, th.accent, 1.5, 0.9);
      var pk = seg(t, 6.6, 7.4);
      if (pk > 0) {
        var px = x0 + s.bi * w, py = y0 - s.best * sc;
        g.line(px, py, px, y0, th.accent, 1.2, pk, [4, 3]);
        g.dot(px, py, 4.5, th.accent, pk);
        g.text(L("enrichment score", "富集分数"), px + 8, py - 4, { col: th.accent2, size: 11, bold: true, a: pk });
      }
    }
  });

  // ---- 12. trajectory -------------------------------------------------------
  X.register("trajectory", {
    period: 9.5,
    stages: [[0, 2.8, "Cells capture snapshots of a continuous process", "细胞记录了一个连续过程的快照"],
             [2.8, 6, "Fit a path through the cells from a root", "从起点出发，拟合穿过细胞的路径"],
             [6, 9.5, "Pseudotime: each cell's position along the path", "拟时序：细胞在路径上的位置"]],
    still: 8.4,
    init: function (R) {
      function path(u, br) {                    // u in [0,1]
        if (u < 0.45) return [70 + u / 0.45 * 200, 145];
        var v = (u - 0.45) / 0.55;
        return [270 + v * 230, 145 + (br ? 1 : -1) * Math.pow(v, 0.9) * 75];
      }
      var pts = [];
      for (var i = 0; i < 130; i++) {
        var u = R(), br = R() < 0.5, p = path(u, br);
        pts.push({ x: p[0] + H.gauss(R) * 9, y: p[1] + H.gauss(R) * 9, u: u });
      }
      return { pts: pts, path: path };
    },
    draw: function (g, s, t) {
      var th = g.th, fit = seg(t, 3, 5.4, false), col = seg(t, 6.2, 8.2, false);
      s.pts.forEach(function (p) {
        var c = p.u <= col ? H.ramp(H.VIRIDIS, p.u) : (th.dark ? "#6b7785" : "#a7b0ba");
        g.dot(p.x, p.y, 2.7, c, seg(t, 0.1, 1.2) * 0.85);
      });
      [0, 1].forEach(function (br) {
        var pts = [];
        for (var u = 0; u <= fit + 1e-6; u += 0.02) pts.push(s.path(u, br));
        g.poly(pts, th.ink, 2.6, 0.85);
      });
      var ra = seg(t, 3, 3.6);
      g.dot(70, 145, 6, th.accent, ra);
      g.t("root", "起点", 70, 128, { col: th.accent2, size: 10.5, bold: true, align: "center", a: ra });
      var cb = seg(t, 6.2, 7);
      for (var i = 0; i < 30; i++) g.rect(390 + i * 4, 40, 4.2, 8, H.ramp(H.VIRIDIS, i / 29), cb);
      g.t("early", "早", 388, 56, { col: th.muted, size: 9.5, a: cb });
      g.t("late", "晚", 512, 56, { col: th.muted, size: 9.5, align: "right", a: cb });
    }
  });

  // ---- 13. RNA velocity -----------------------------------------------------
  X.register("velocity", {
    period: 9.5,
    stages: [[0, 3, "Unspliced RNA is new; spliced RNA is mature", "未剪接 RNA 是新转录的，剪接 RNA 是成熟的"],
             [3, 6.2, "Their balance says whether a gene is turning on or off", "两者的比例说明基因在开启还是关闭"],
             [6.2, 9.5, "Arrows show where each cell is heading", "箭头显示每个细胞的去向"]],
    still: 8.6,
    init: function (R) {
      function arc(u) { return [50 + u * 290, 205 - Math.sin(u * 2.4) * 120 - u * 20]; }
      var cells = [];
      for (var i = 0; i < 60; i++) { var u = R(), p = arc(u); cells.push({ x: p[0] + H.gauss(R) * 8, y: p[1] + H.gauss(R) * 8, u: u }); }
      return { arc: arc, cells: cells };
    },
    draw: function (g, s, t) {
      var th = g.th, ca = seg(t, 0.1, 1);
      // phase portrait inset
      var ix = 380, iy = 215, iw = 150, ih = 150, pa = seg(t, 0.5, 1.8);
      g.axes(ix, iy, iw, ih, L("spliced", "剪接"), L("unspliced", "未剪接"), pa);
      g.line(ix, iy, ix + iw, iy - ih * 0.8, th.faint, 1, pa, [4, 3]);
      var up = [], down = [];
      for (var k = 0; k <= 20; k++) {
        var q = k / 20;
        up.push([ix + iw * q * 0.95, iy - ih * (0.9 * Math.sqrt(q))]);
        down.push([ix + iw * 0.95 * (1 - q), iy - ih * 0.9 * (1 - q) * (1 - q) * 0.75 - 2]);
      }
      g.poly(up, PAL[1], 2, pa * 0.9);
      g.poly(down, PAL[0], 2, pa * 0.9);
      g.t("turning on", "开启中", ix + 30, iy - ih + 18, { col: PAL[1], size: 9.5, a: seg(t, 3.2, 4) });
      g.t("turning off", "关闭中", ix + iw - 5, iy - 30, { col: PAL[0], size: 9.5, align: "right", a: seg(t, 3.2, 4) });
      // cells on the arc, then arrows and flowing streams
      var ar = seg(t, 6.4, 7.4);
      s.cells.forEach(function (c) {
        g.dot(c.x, c.y, 2.8, H.ramp(H.VIRIDIS, c.u * 0.9), ca * 0.85);
        if (ar > 0) {
          var a = s.arc(Math.min(1, c.u + 0.04)), b = s.arc(c.u), dx = a[0] - b[0], dy = a[1] - b[1], n = Math.sqrt(dx * dx + dy * dy) || 1;
          g.arrow(c.x, c.y, c.x + dx / n * 15 * ar, c.y + dy / n * 15 * ar, th.ink, 1.1, 0.75, 4.5);
        }
      });
      if (ar > 0) {
        for (var f = 0; f < 8; f++) {
          var u = ((t * 0.12 + f / 8) % 1), p = s.arc(u);
          g.dot(p[0], p[1] - 14, 2.2, th.accent, ar * Math.sin(u * Math.PI));
        }
      }
    }
  });

  // ---- 14. dynamic features -------------------------------------------------
  X.register("dynamic", {
    period: 9.5,
    stages: [[0, 2.8, "Order the cells by pseudotime", "按拟时序排列细胞"],
             [2.8, 6.2, "Fit each gene's expression along it", "拟合每个基因沿拟时序的表达"],
             [6.2, 9.5, "Genes that change are dynamic features", "随时间变化的基因就是动态特征"]],
    still: 8.4,
    init: function (R) {
      var cells = [];
      for (var i = 0; i < 60; i++) cells.push({ u: R(), j: H.gauss(R) });
      return { cells: cells };
    },
    draw: function (g, s, t) {
      var th = g.th, x0 = 70, y0 = 165, w = 440, h = 120;
      g.axes(x0, y0, w, h, "", L("expression", "表达"), seg(t, 0, 0.6));
      s.cells.forEach(function (c) {
        g.dot(x0 + c.u * w, y0 + 10 + c.j * 2, 2.6, H.ramp(H.VIRIDIS, c.u), seg(t, 0.2 + c.u, 0.6 + c.u));
      });
      g.t("pseudotime →", "拟时序 →", x0 + w, y0 + 22, { col: th.muted, size: 10, align: "right", a: seg(t, 0.4, 1) });
      var f = [function (u) { return 1 / (1 + Math.exp(-(u - 0.5) * 10)); },
               function (u) { return 1 - 1 / (1 + Math.exp(-(u - 0.4) * 9)); },
               function (u) { return Math.exp(-Math.pow((u - 0.55) / 0.15, 2)); }];
      var dr = seg(t, 3, 5.6, false);
      f.forEach(function (fn, k) {
        var pts = [];
        for (var u = 0; u <= dr + 1e-6; u += 0.01) pts.push([x0 + u * w, y0 - 8 - fn(u) * (h - 16)]);
        g.poly(pts, PAL[k], 2.3, 1);
      });
      var hm = seg(t, 6.4, 7.6, false), genes = [2, 0, 1, 2, 0, 1];
      for (var r = 0; r < 6; r++) {
        if (r >= hm * 6) break;
        for (var cI = 0; cI < 40; cI++) {
          var u2 = cI / 39, sh = [0.08, -0.06, 0.05, -0.1, 0.12, -0.03][r];
          var v = f[genes[r]](H.clamp(u2 + sh, 0, 1));
          g.rect(x0 + cI * (w / 40), 196 + r * 9, w / 40 + 0.3, 8.4, H.mix(th.dark ? "#1b2530" : "#f1f4f7", th.accent, v), 1);
        }
      }
    }
  });

  // ---- 15. cell cycle & signatures ------------------------------------------
  X.register("cellcycle", {
    period: 9.5,
    stages: [[0, 3, "Cells cycle through G1, S and G2/M", "细胞在 G1、S、G2/M 之间循环"],
             [3, 6.2, "Score S-phase and G2/M gene sets in every cell", "为每个细胞计算 S 期与 G2/M 期基因集得分"],
             [6.2, 9.5, "The two scores assign a phase", "两个得分决定周期时相"]],
    still: 8.4,
    init: function (R) {
      var cells = [];
      for (var i = 0; i < 90; i++) {
        var ph = R() < 0.55 ? 0 : (R() < 0.5 ? 1 : 2);
        var sx = ph === 1 ? 0.55 + H.gauss(R) * 0.13 : -0.15 + H.gauss(R) * 0.13;
        var sy = ph === 2 ? 0.55 + H.gauss(R) * 0.13 : -0.15 + H.gauss(R) * 0.13;
        cells.push({ ph: ph, sx: sx, sy: sy, a0: R() * H.TAU, v: 0.5 + R() * 0.5 });
      }
      return { cells: cells };
    },
    draw: function (g, s, t) {
      var th = g.th, cx = 135, cy = 140, r = 70, C = [PAL[0], PAL[4], PAL[2]];
      var a = seg(t, 0.1, 1);
      var arcs = [[-Math.PI / 2, Math.PI * 0.5, "G1"], [Math.PI * 0.5, Math.PI, "S"], [Math.PI, Math.PI * 1.5, "G2/M"]];
      arcs.forEach(function (ac, k) {
        g.arc(cx, cy, r, ac[0], ac[1], C[k], 9, a * 0.7);
        var mid = (ac[0] + ac[1]) / 2;
        g.text(ac[2], cx + Math.cos(mid) * (r + 22), cy + Math.sin(mid) * (r + 22), { col: C[k], size: 11.5, bold: true, align: "center", a: a });
      });
      for (var i = 0; i < 14; i++) {
        var ang = -Math.PI / 2 + ((t * 0.5 + i / 14) % 1) * H.TAU;
        g.dot(cx + Math.cos(ang) * r, cy + Math.sin(ang) * r, 3.5, th.ink, a * 0.8);
      }
      var x0 = 330, y0 = 225, w = 200, h = 175, sa = seg(t, 3.2, 4.4);
      g.axes(x0, y0, w, h, "S score", "G2/M score", sa);
      var cl = seg(t, 6.4, 7.4);
      if (cl > 0) {
        g.line(x0 + w * 0.42, y0, x0 + w * 0.42, y0 - h, th.faint, 1, cl, [4, 3]);
        g.line(x0, y0 - h * 0.42, x0 + w, y0 - h * 0.42, th.faint, 1, cl, [4, 3]);
      }
      s.cells.forEach(function (c, i) {
        var p = seg(t, 3.4 + (i % 30) * 0.05, 4.2 + (i % 30) * 0.05);
        g.dot(x0 + (c.sx + 0.4) / 1.3 * w, y0 - (c.sy + 0.4) / 1.3 * h, 2.6,
              H.mix(th.dark ? "#6b7785" : "#a7b0ba", C[c.ph], cl), p * 0.85);
      });
    }
  });

  // ---- 16. cell-cell communication ------------------------------------------
  X.register("cellcomm", {
    period: 9.5,
    stages: [[0, 2.6, "Cell types send and receive signals", "细胞类型之间收发信号"],
             [2.6, 6.2, "A ligand on the sender meets a receptor on the receiver", "发送方的配体与接收方的受体结合"],
             [6.2, 9.5, "Edge width = how many significant pairs", "连线粗细 = 显著配体–受体对的数量"]],
    still: 8.4,
    draw: function (g, s, t) {
      var th = g.th, cx = 280, cy = 145, R0 = 88;
      var names = [["T", "T"], ["B", "B"], ["Mono", "单核"], ["NK", "NK"], ["Fibro", "成纤维"]];
      var pos = names.map(function (n, i) { var a = -Math.PI / 2 + i / 5 * H.TAU; return [cx + Math.cos(a) * R0 * 1.5, cy + Math.sin(a) * R0]; });
      var edges = [[4, 2, 0.95], [2, 0, 0.8], [0, 1, 0.5], [2, 3, 0.4], [3, 0, 0.3]];
      edges.forEach(function (e, k) {
        var p = seg(t, 2.8 + k * 0.5, 3.6 + k * 0.5, false);
        if (p <= 0) return;
        var a = pos[e[0]], b = pos[e[1]], mx = (a[0] + b[0]) / 2 + (cy - (a[1] + b[1]) / 2) * 0.35,
            my = (a[1] + b[1]) / 2 - (cx - (a[0] + b[0]) / 2) * 0.35;
        var wdt = lerp(1.6, 1.6 + e[2] * 7, seg(t, 6.4, 7.4));
        g.curve(a[0], a[1], mx, my, b[0], b[1], PAL[e[0]], wdt, 0.55 * Math.min(1, p * 2));
        for (var f = 0; f < 3; f++) {
          var u = (t * 0.45 + f / 3 + k * 0.17) % 1;
          var x = (1 - u) * (1 - u) * a[0] + 2 * (1 - u) * u * mx + u * u * b[0], y = (1 - u) * (1 - u) * a[1] + 2 * (1 - u) * u * my + u * u * b[1];
          g.dot(x, y, 2.4, PAL[e[0]], Math.min(1, p * 2) * 0.9);
        }
      });
      pos.forEach(function (p, i) {
        var a = seg(t, 0.2 + i * 0.15, 0.8 + i * 0.15);
        g.dot(p[0], p[1], 17, PAL[i], a * 0.9);
        g.text(L(names[i][0], names[i][1]), p[0], p[1], { col: "#ffffff", size: 10, bold: true, align: "center", a: a });
      });
      g.chip("CXCL12 → CXCR4", 470, 236, PAL[4], seg(t, 3, 3.8) * (1 - seg(t, 6.2, 6.8)), { size: 10 });
    }
  });

  // ---- 17. malignancy / CNV -------------------------------------------------
  X.register("malignancy", {
    period: 9.5,
    stages: [[0, 3, "Average expression along each chromosome", "沿每条染色体平均基因表达"],
             [3, 6.2, "Compare with reference (normal) cells", "与参照（正常）细胞比较"],
             [6.2, 9.5, "Large gains and losses mark malignant cells", "大片段的扩增与缺失提示恶性细胞"]],
    still: 8.6,
    init: function (R) {
      var sizes = [249, 243, 198, 191, 181, 171, 159, 145, 138, 134, 135, 133, 115, 107, 102, 90, 83, 80, 59, 64, 47, 51];
      var tot = sizes.reduce(function (a, b) { return a + b; }, 0), x = 60, chr = [];
      sizes.forEach(function (sz, i) { var w = sz / tot * 420; chr.push({ x: x, w: w, n: i + 1 }); x += w; });
      var rows = [];
      for (var r = 0; r < 14; r++) {
        var tum = [1, 2, 4, 5, 7, 8, 10, 11, 12].indexOf(r) >= 0, v = [];
        for (var c = 0; c < 22; c++) {
          var base = H.gauss(R) * 0.12;
          if (tum && (c === 2 || c === 7)) base += 0.8;
          if (tum && (c === 8 || c === 16)) base -= 0.75;
          v.push(base);
        }
        rows.push({ tum: tum, v: v });
      }
      var order = rows.map(function (r, i) { return i; }).sort(function (a, b) { return (rows[a].tum ? 1 : 0) - (rows[b].tum ? 1 : 0) || a - b; });
      rows.forEach(function (r, i) { r.final = order.indexOf(i); });
      return { chr: chr, rows: rows };
    },
    draw: function (g, s, t) {
      var th = g.th, y0 = 50, rh = 12, sweep = seg(t, 0.2, 2.6, false), colr = seg(t, 3.2, 4.6), sort = seg(t, 6.4, 7.6);
      var neutral = th.dark ? "#2a3440" : "#eef1f4";
      s.chr.forEach(function (c, i) {
        if (i % 2 === 0 || i < 9) g.text(String(c.n), c.x + c.w / 2, y0 - 8, { col: th.muted, size: 8.5, align: "center", a: seg(t, 0, 0.6) });
      });
      s.rows.forEach(function (r, i) {
        var y = y0 + lerp(i, r.final, sort) * rh;
        s.chr.forEach(function (c, j) {
          if (c.x + c.w > 60 + 420 * sweep + 0.5) return;
          var v = r.v[j] * colr, col = v > 0 ? H.mix(neutral, "#d6403a", H.clamp(v, 0, 1)) : H.mix(neutral, "#3473c4", H.clamp(-v, 0, 1));
          g.rect(c.x + 0.4, y + 0.5, c.w - 0.8, rh - 1, col, 1);
        });
        var lab = seg(t, 6.8, 7.6);
        if (lab > 0) g.text(r.tum ? L("malignant", "恶性") : L("normal", "正常"), 488, y + rh / 2,
                            { col: r.tum ? H.FLAG : th.ok, size: 9, a: lab });
      });
      var ref = seg(t, 3.2, 3.8) * (1 - sort);
      s.rows.forEach(function (r, i) { if (!r.tum) g.box(58, y0 + i * rh, 424, rh, th.ok, 1, ref * 0.7, 2); });
      g.t("reference cells", "参照细胞", 58, 226, { col: th.ok, size: 10, align: "left", a: ref });
      g.t("gain", "扩增", 300, 236, { col: "#d6403a", size: 10, bold: true, a: colr });
      g.t("loss", "缺失", 340, 236, { col: "#3473c4", size: 10, bold: true, a: colr });
      g.note("Pick reference cells from the same patient", "参照细胞应来自同一位患者", seg(t, 3.4, 4.2));
    }
  });

  // ---- 18. clinical & survival ----------------------------------------------
  X.register("clinical", {
    period: 10.5,
    stages: [[0, 3.2, "Summarise the cells of each patient into one value", "把每位患者的细胞汇总成一个数值"],
             [3.2, 6, "Split patients at the median of that value", "按该数值的中位数把患者分组"],
             [6, 10.5, "Compare survival between the groups", "比较两组的生存"]],
    still: 9.2,
    init: function (R) {
      var pats = [];
      for (var i = 0; i < 8; i++) {
        var frac = 0.15 + R() * 0.7, cells = [];
        for (var c = 0; c < 10; c++) cells.push(R() < frac);
        pats.push({ frac: cells.filter(Boolean).length / 10, cells: cells });
      }
      var sorted = pats.slice().sort(function (a, b) { return a.frac - b.frac; });
      pats.forEach(function (p) { p.rank = sorted.indexOf(p); p.hi = p.rank >= 4; });
      return { pats: pats, lo: kmCurve(R, 40, 1 / 30, 60), hi: kmCurve(R, 40, 1 / 13, 60) };
    },
    draw: function (g, s, t) {
      var th = g.th, agg = seg(t, 1.6, 2.8), sp = seg(t, 3.4, 4.6), x0 = 30;
      s.pats.forEach(function (p, i) {
        var x = x0 + 18 + lerp(i, p.rank, sp) * 28;
        g.text("P" + (i + 1), x0 + 18 + i * 28, 228, { col: th.muted, size: 9, align: "center", a: seg(t, 0.2, 0.8) * (1 - sp) });
        p.cells.forEach(function (on, c) {
          var y = 70 + c * 14;
          g.dot(x0 + 18 + i * 28, lerp(y, 160, agg), 4 * (1 - agg) + 0.5, on ? PAL[3] : th.faint, seg(t, 0.2, 1) * (1 - agg));
        });
        var bh = p.frac * 110 * agg;
        g.rect(x - 9, 215 - bh, 18, bh, sp > 0 ? H.mix(PAL[3], p.hi ? PAL[1] : PAL[0], sp) : PAL[3], 0.85, 3);
      });
      var md = seg(t, 3.6, 4.4);
      g.line(x0 + 4, 215 - 0.5 * (s.pats[3].frac + s.pats[4].frac) * 110, x0 + 240, 215 - 0.5 * (s.pats[3].frac + s.pats[4].frac) * 110, th.ink, 1, md * 0.7, [4, 3]);
      g.t("median", "中位数", x0 + 4, 207 - 0.5 * (s.pats[3].frac + s.pats[4].frac) * 110, { col: th.muted, size: 9.5, a: md });
      g.t("fraction of cell type", "细胞类型占比", x0 + 128, 50, { col: PAL[3], size: 10, align: "center", a: agg });
      var km = seg(t, 6.2, 9, false), kx = 320, ky = 210, kw = 210, kh = 160;
      if (km > 0) {
        g.axes(kx, ky, kw, kh, L("months", "月"), L("survival", "生存率"), seg(t, 6, 6.6));
        drawKM(g, s.lo, kx, ky, kw, kh, 60, PAL[0], km);
        drawKM(g, s.hi, kx, ky, kw, kh, 60, PAL[1], km);
        g.t("low", "低", kx + kw + 4, ky - 0.62 * kh, { col: PAL[0], size: 10, bold: true, a: seg(t, 8.6, 9.2) });
        g.t("high", "高", kx + kw + 4, ky - 0.18 * kh, { col: PAL[1], size: 10, bold: true, a: seg(t, 8.6, 9.2) });
      }
      g.note("One row per patient; cells are never the unit of a survival test", "每位患者一行；生存检验的单位永远不是细胞", seg(t, 6.4, 7.2));
    }
  });

  // ---- 19. visualize --------------------------------------------------------
  X.register("viz", {
    period: 9.5,
    stages: [[0, 3, "An embedding coloured by cluster", "按聚类着色的降维图"],
             [3, 6.2, "The same map, coloured by one gene", "同一张图，按一个基因的表达着色"],
             [6.2, 9.5, "Or the same cells as violins per cluster", "或者按聚类画成小提琴图"]],
    still: 8.4,
    init: function (R) {
      var pts = clusters(R, [[150, 110, 24, 16], [260, 185, 24, 16], [370, 105, 24, 16], [460, 180, 20, 14]]);
      pts.forEach(function (p) {
        p.e = H.clamp((p.k === 2 ? 0.75 : 0.12) + H.gauss(R) * 0.15, 0, 1);
        p.vx = 140 + p.k * 100 + H.gauss(R) * (6 + 12 * Math.sin(p.e * Math.PI));
        p.vy = 215 - p.e * 160;
      });
      return { pts: pts };
    },
    draw: function (g, s, t) {
      var th = g.th, fe = seg(t, 3.2, 4.2), vi = seg(t, 6.4, 7.8);
      s.pts.forEach(function (p) {
        var c = H.mix(PAL[p.k], H.ramp(H.VIRIDIS, p.e), fe);
        g.dot(lerp(p.x, p.vx, vi), lerp(p.y, p.vy, vi), 2.7, vi > 0.5 ? PAL[p.k] : c, seg(t, 0.1, 1) * 0.85);
      });
      g.chip("MS4A1", 500, 42, th.accent, fe * (1 - vi), { size: 10 });
      if (vi > 0) {
        g.axes(110, 222, 410, 180, "", "MS4A1", vi);
        for (var k = 0; k < 4; k++) g.text(String(k), 140 + k * 100, 233, { col: th.muted, size: 10, align: "center", a: vi });
      }
    }
  });

  // ---- 20. report -----------------------------------------------------------
  X.register("report", {
    period: 8.5,
    stages: [[0, 2.8, "Every step you ran is logged", "你运行的每一步都被记录"],
             [2.8, 5.8, "Parameters and code are written up", "参数与代码被整理成文"],
             [5.8, 8.5, "Download it as HTML or Markdown", "下载为 HTML 或 Markdown"]],
    still: 7.2,
    init: function (R) {
      var lines = [];
      for (var i = 0; i < 9; i++) lines.push(0.5 + R() * 0.45);
      return { lines: lines };
    },
    draw: function (g, s, t) {
      var th = g.th, dx = 210, dy = 40, dw = 170, dh = 196;
      var steps = ["QC", "Normalize", "PCA", "Cluster", "UMAP"];
      steps.forEach(function (st, i) {
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
      g.chip("report.md", 460, 145, th.muted, dl, { size: 10.5 });
      g.arrow(dx + dw + 8, 128, dx + dw + 8 + 30 * dl, 128, th.accent, 2, dl);
    }
  });

  // ---- 21. export -----------------------------------------------------------
  X.register("export", {
    period: 8.5,
    stages: [[0, 2.6, "The analysed object", "分析完成的对象"],
             [2.6, 5.6, "Save the object, the tables and a script", "保存对象、表格与脚本"],
             [5.6, 8.5, "Run the script to reproduce the analysis", "运行脚本即可复现分析"]],
    still: 7.2,
    draw: function (g, s, t) {
      var th = g.th, ox = 150, oy = 140, a = seg(t, 0.1, 1);
      for (var i = 0; i < 3; i++) g.rect(ox - 44 + i * 6, oy - 44 + i * 6, 76, 76, [PAL[0], PAL[2], PAL[3]][i], a * 0.35 + i * 0.15, 10);
      g.text("obj", ox, oy + 2, { col: "#ffffff", size: 14, bold: true, align: "center", a: a });
      var files = [".rds", ".csv", ".h5ad", "analysis.R"];
      files.forEach(function (f, k) {
        var p = seg(t, 2.8 + k * 0.45, 3.8 + k * 0.45);
        if (p <= 0) return;
        var tx = 400, ty = 70 + k * 40, mx = (ox + tx) / 2, my = Math.min(oy, ty) - 50;
        var x = (1 - p) * (1 - p) * ox + 2 * (1 - p) * p * mx + p * p * tx, y = (1 - p) * (1 - p) * oy + 2 * (1 - p) * p * my + p * p * ty;
        g.chip(f, x, y, k === 3 ? th.ok : th.accent, p, { size: 10.5, bold: k === 3 });
      });
      var rr = seg(t, 5.8, 7);
      if (rr > 0) {
        g.curve(400, 210, 290, 250, ox + 10, oy + 44, th.ok, 1.8, rr);
        g.arrow(ox + 24, oy + 52, ox + 12, oy + 44, th.ok, 1.8, rr, 7);
        g.t("re-run → same result", "重新运行 → 相同结果", 300, 238, { col: th.ok, size: 10.5, align: "center", a: rr });
      }
    }
  });
})();
