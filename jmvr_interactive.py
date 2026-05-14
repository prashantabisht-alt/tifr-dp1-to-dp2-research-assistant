"""
Build a single self-contained HTML file with embedded precomputed JMVR data.
Drops the result at /TIFR_DP2/jmvr_interactive.html so you can open it
directly in any browser. Sliders for gamma and epsilon. Two lattices side by side.
"""
import json
import os

here = os.path.dirname(__file__)
data_path = os.path.join(os.path.dirname(here), os.pardir, "widget_data_final.json")
data_path = os.path.abspath(data_path)
out_path = os.path.join(here, "jmvr_interactive.html")

with open(data_path) as f:
    data = json.load(f)

# Compact JSON for embedding
data_json = json.dumps(data, separators=(",", ":"))

html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>JMVR active random walker — interactive</title>
<style>
:root {
  --bg: #fafaf8;
  --panel: #ffffff;
  --border: rgba(0,0,0,0.12);
  --text: #1a1a1a;
  --muted: #6a6a6a;
  --accent: #5b3aa8;
}
* { box-sizing: border-box; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  margin: 0; padding: 24px;
  background: var(--bg); color: var(--text);
}
h1 { font-size: 22px; font-weight: 500; margin: 0 0 4px; }
.sub { color: var(--muted); font-size: 14px; margin-bottom: 24px; }
.controls {
  display: grid; grid-template-columns: auto 1fr auto;
  gap: 10px 18px; align-items: center;
  background: var(--panel); padding: 16px 20px;
  border: 1px solid var(--border); border-radius: 10px;
  max-width: 720px; margin-bottom: 24px;
}
.controls label { font-size: 14px; color: var(--muted); }
.controls input[type=range] { width: 100%; }
.controls .val { font-size: 14px; font-weight: 500; min-width: 64px; text-align: right; font-variant-numeric: tabular-nums; }
.panels { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; max-width: 900px; }
.panel {
  background: var(--panel); border: 1px solid var(--border);
  border-radius: 10px; padding: 16px;
}
.panel h2 { font-size: 16px; font-weight: 500; margin: 0 0 4px; }
.panel .descr { font-size: 13px; color: var(--muted); margin-bottom: 12px; }
canvas { width: 100%; height: auto; image-rendering: pixelated; image-rendering: crisp-edges; display: block; }
.footer { color: var(--muted); font-size: 12px; margin-top: 18px; max-width: 900px; }
</style>
</head>
<body>

<h1>JMVR active random walker</h1>
<div class="sub">Single walker. P(x, y, t = """ + str(int(data["t"])) + """). Sheared k-grid + corrected c<sub>3</sub> sign.</div>

<div class="controls">
  <label for="g">γ (rotation rate)</label>
  <input id="g" type="range" min="0" max="5" step="1" value="2">
  <span class="val" id="gv">0.020</span>

  <label for="e">ε (chirality)</label>
  <input id="e" type="range" min="0" max="5" step="1" value="3">
  <span class="val" id="ev">—</span>
</div>

<div class="panels">
  <div class="panel">
    <h2>Square lattice</h2>
    <div class="descr">4 directors. Sites at integer (n₁, n₂). Walker rolled to centre.</div>
    <canvas id="sq" width="16" height="16"></canvas>
  </div>
  <div class="panel">
    <h2>Triangular lattice</h2>
    <div class="descr">6 directors. Sheared lattice basis. Walker rolled to centre.</div>
    <canvas id="tr" width="16" height="16"></canvas>
  </div>
</div>

<div class="footer">
  Pre-computed on an L×L torus at fixed t. Slider snaps to grid values:
  γ ∈ {0.005, 0.01, 0.02, 0.05, 0.10, 0.20}, ε swept across each lattice's physical range.
  Probability rendered with the same black→purple→red→yellow palette used in Confinement-2021 Fig 11.
</div>

<script>
const DATA = __DATA__;
const L = DATA.L;
const PALETTE = [[5,5,5],[50,16,95],[180,43,214],[237,28,36],[255,140,0],[255,242,0]];
const STOPS = [0.00, 0.22, 0.42, 0.60, 0.78, 1.00];

function colorAt(t) {
  t = Math.max(0, Math.min(1, t));
  for (let i = 1; i < STOPS.length; i++) {
    if (t <= STOPS[i]) {
      const u = (t - STOPS[i-1]) / (STOPS[i] - STOPS[i-1]);
      const a = PALETTE[i-1], b = PALETTE[i];
      return [a[0]+u*(b[0]-a[0]), a[1]+u*(b[1]-a[1]), a[2]+u*(b[2]-a[2])];
    }
  }
  return PALETTE[PALETTE.length-1];
}

function renderHeatmap(canvasId, flat, vmax) {
  const c = document.getElementById(canvasId);
  c.width = L; c.height = L;
  const ctx = c.getContext('2d');
  const img = ctx.createImageData(L, L);
  for (let n2 = 0; n2 < L; n2++) {
    for (let n1 = 0; n1 < L; n1++) {
      const p = flat[n2 * L + n1];
      const t = p / vmax;
      const col = colorAt(t);
      const idx = 4 * (n2 * L + n1);
      img.data[idx] = col[0]; img.data[idx+1] = col[1]; img.data[idx+2] = col[2]; img.data[idx+3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
  // Scale up via CSS for crisp pixels (already pixelated via image-rendering)
  c.style.width = '320px'; c.style.height = '320px';
}

function vmaxOf(flat) { let m = 0; for (const v of flat) if (v > m) m = v; return m; }

function update() {
  const gi = +document.getElementById('g').value;
  const ei = +document.getElementById('e').value;
  const gamma = DATA.gammas[gi];
  const epsSq = DATA.eps_square[ei];
  const epsTr = DATA.eps_triangular[ei];
  document.getElementById('gv').textContent = gamma.toFixed(3);
  document.getElementById('ev').textContent = 'sq ' + epsSq.toFixed(2) + ' | tr ' + epsTr.toFixed(2);
  const sqFlat = DATA.square[gi][ei];
  const trFlat = DATA.triangular[gi][ei];
  renderHeatmap('sq', sqFlat, vmaxOf(sqFlat));
  renderHeatmap('tr', trFlat, vmaxOf(trFlat));
}
document.getElementById('g').addEventListener('input', update);
document.getElementById('e').addEventListener('input', update);
update();
</script>
</body>
</html>
"""

html = html.replace("__DATA__", data_json)
with open(out_path, "w") as f:
    f.write(html)

print(f"Saved: {out_path}")
print(f"Size: {os.path.getsize(out_path):,} bytes")
