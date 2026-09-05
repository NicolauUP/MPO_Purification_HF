"""Render the fixed-cap QTT-MPO SP2 SCF scaling control.

The measurements are timing-only SCF solves: no full real-space observable
extraction, checkpoint write, or energy post-processing is included.
"""

from pathlib import Path


# From job 45358745 (128^2) and array job 45359269 (256^2--2048^2).
side = [128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]
qtt_levels = [14, 16, 18, 20, 22, 24, 26, 28, 30]
# Mean over SCF iterations 2:end. The first iteration includes initialization
# and is deliberately excluded from this steady-state scaling figure.
steady_iteration_s = [9.325, 10.329, 11.226, 12.587, 12.421,
                      12.446, 11.296, 11.688, 11.345]
purification_s = [5.738, 6.509, 6.908, 7.496, 7.710,
                  7.547, 6.829, 7.175, 7.116]
mean_field_s = [1.626, 1.659, 1.981, 2.576, 2.272,
                2.529, 2.130, 2.114, 1.989]

output = Path(__file__).resolve().parent / "qtt_mpo_sp2_scf_scaling.svg"

# Dependency-free SVG renderer: intentional, so the figure can be reproduced on
# the login machine without installing a plotting package.
width, height = 900, 540
left, right, top, bottom = 110, 45, 75, 100
xmin, xmax, ymin, ymax = 13.5, 30.5, 0.0, 14.0
sx = lambda x: left + (x - xmin) / (xmax - xmin) * (width - left - right)
sy = lambda y: height - bottom - (y - ymin) / (ymax - ymin) * (height - top - bottom)

series = [
    ("full SCF iteration", steady_iteration_s, "#1769aa", "circle"),
    ("SP2 purification", purification_s, "#d95f02", "square"),
    ("Hartree/Fock update", mean_field_s, "#2a9d8f", "triangle"),
]
parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<style>text{{font-family:Arial,sans-serif;fill:#222}} .tick{{font-size:14px}} .label{{font-size:17px}} .title{{font-size:20px;font-weight:bold}}</style>
<text x="{width/2}" y="32" text-anchor="middle" class="title">Fixed-χmax = 256 QTT-MPO SP2 SCF scaling (V₂ = 0, U = 1)</text>''']

for y in range(0, 15, 2):
    yy = sy(y)
    parts.append(f'<line x1="{left}" y1="{yy:.1f}" x2="{width-right}" y2="{yy:.1f}" stroke="#d9d9d9"/>')
    parts.append(f'<text x="{left-12}" y="{yy+5:.1f}" text-anchor="end" class="tick">{y}</text>')
for x in qtt_levels:
    xx = sx(x)
    parts.append(f'<line x1="{xx:.1f}" y1="{height-bottom}" x2="{xx:.1f}" y2="{height-bottom+6}" stroke="#222"/>')
    parts.append(f'<text x="{xx:.1f}" y="{height-bottom+26}" text-anchor="middle" class="tick">{x}</text>')

parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#222" stroke-width="1.5"/>')
parts.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#222" stroke-width="1.5"/>')

for label, values, colour, marker in series:
    points = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in zip(qtt_levels, values))
    parts.append(f'<polyline points="{points}" fill="none" stroke="{colour}" stroke-width="3"/>')
    for x, y in zip(qtt_levels, values):
        xx, yy = sx(x), sy(y)
        if marker == "circle":
            parts.append(f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="6" fill="{colour}"/>')
        elif marker == "square":
            parts.append(f'<rect x="{xx-5:.1f}" y="{yy-5:.1f}" width="10" height="10" fill="{colour}"/>')
        else:
            parts.append(f'<path d="M {xx:.1f} {yy-6:.1f} L {xx-6:.1f} {yy+5:.1f} L {xx+6:.1f} {yy+5:.1f} Z" fill="{colour}"/>')

for x, n, t in zip(qtt_levels, side, steady_iteration_s):
    label = "1.07B sites" if n == 32768 else f"{n*n:,} sites"
    parts.append(f'<text x="{sx(x):.1f}" y="{sy(t)-12:.1f}" text-anchor="middle" class="tick">{label}</text>')

legend_x, legend_y = 130, 92
for i, (label, _, colour, _) in enumerate(series):
    y = legend_y + i*27
    parts.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x+28}" y2="{y}" stroke="{colour}" stroke-width="3"/>')
    parts.append(f'<text x="{legend_x+38}" y="{y+5}" class="tick">{label}</text>')

parts.append(f'<text x="{width/2}" y="{height-28}" text-anchor="middle" class="label">Total QTT length L_QTT = 2 log₂(L_side)</text>')
parts.append(f'<text x="27" y="{height/2}" text-anchor="middle" class="label" transform="rotate(-90 27 {height/2})">steady-state wall time per SCF iteration (s)</text>')
parts.append('</svg>')
output.write_text("\n".join(parts))
print(output)
