"""
plot_figure1AB.py
  Figure A — Handedness 与SD的 LDSC cor network plot
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

 # BASE = r''

# ───────────────────────────────────────────────────────────────────────
# LDSC rg（ldsc_handedness_disease_rg.csv）
ldsc_rg = {
    'SCZ':  {'rg':  0.1301, 'p': 0.000378},
    'ADHD': {'rg': -0.0542, 'p': 0.2293},
    'PTSD': {'rg':  0.0348, 'p': 0.3861},
}

# h²（LDSC log，SE）
h2 = {
    'LH':   {'h2': 0.0106, 'se': 0.0010},
    'SCZ':  {'h2': 0.1008, 'se': 0.0035},
    'ADHD': {'h2': 0.0947, 'se': 0.0045},
    'PTSD': {'h2': 0.0282, 'se': 0.0011},
}

deep = {'LH': '#4DA6D0', 'SCZ': '#8B4DAB', 'ADHD': '#D94040', 'PTSD': '#34A86B'}
lite = {'SCZ': '#CBAADF', 'ADHD': '#F0AAAA', 'PTSD': '#96D4B3'}

# ═══════════════════════════════════════════════════════════════════════════
# Figure A 
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(5.5, 5.5))
ax.set_xlim(-1.95, 1.95)
ax.set_ylim(-1.85, 1.80)
ax.set_aspect('equal')
ax.axis('off')

pos = {
    'LH':   np.array([ 0.00,  0.05]),
    'SCZ':  np.array([-1.15,  1.00]),
    'ADHD': np.array([ 1.15,  1.00]),
    'PTSD': np.array([ 0.00, -1.20]),
}

h2_min = min(v['h2'] for v in h2.values())
h2_max = max(v['h2'] for v in h2.values())
r_min, r_max = 0.22, 0.38
node_r = {
    k: r_min + (v['h2'] - h2_min) / (h2_max - h2_min) * (r_max - r_min)
    for k, v in h2.items()
}

lw_ref = 8.0
rg_ref = 0.15

for dis, vals in ldsc_rg.items():
    rg  = vals['rg']
    p   = vals['p']
    sig = p < 0.05

    lw = max(1.5, lw_ref * abs(rg) / rg_ref)
    ec = deep[dis] if sig else lite[dis]

    p1 = pos['LH']; p2 = pos[dis]
    dv = p2 - p1
    dn = dv / np.linalg.norm(dv)

    s = p1 + dn * node_r['LH']
    e = p2 - dn * node_r[dis]

    ax.plot([s[0], e[0]], [s[1], e[1]],
            color=ec, lw=lw, solid_capstyle='round', zorder=1)

    pcw  = np.array([ dn[1], -dn[0]])   
    pccw = np.array([-dn[1],  dn[0]])   

    if dis == 'SCZ':
        pu = pcw    
    elif dis == 'ADHD':
        pu = pccw   
    else:           
        pu = pccw

    mid  = (s + e) / 2
    lpos = mid + pu * 0.30

    rg_s = f'$r_g$ = {rg:.3f}'
    p_s = f'P = {p:.4e}' if p < 0.0001 else f'P = {p:.4f}'
    tc   = deep[dis] if sig else '#999999'
    fw   = 'bold' if sig else 'normal'

    ax.text(lpos[0], lpos[1], f'{rg_s}\n{p_s}',
            ha='center', va='center', fontsize=8.5, color=tc, fontweight=fw,
            zorder=5,
            bbox=dict(boxstyle='round,pad=0.18', fc='white', ec='none', alpha=0.80))

for nm, xy in pos.items():
    circ = plt.Circle(xy, node_r[nm], color=deep[nm], zorder=3,
                      ec='white', linewidth=1.5)
    ax.add_patch(circ)
    ax.text(xy[0], xy[1], nm, ha='center', va='center',
            fontsize=11, fontweight='bold', color='white', zorder=4)

ax.set_title('A', loc='left', fontsize=16, fontweight='bold', pad=4)

fig.tight_layout(pad=0.3)
out_a = f'{BASE}/figure1A_network.png'
out_a_pdf = f'{BASE}/figure1A_network.pdf'
fig.savefig(out_a, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out_a_pdf, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f"Saved: figure1A_network.png and figure1A_network.pdf")

# ═══════════════════════════════════════════════════════════════════════════
# Figure B — h²_SNP  (Moved to R script: plot_figure1B_h2.R)
# ═══════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════

