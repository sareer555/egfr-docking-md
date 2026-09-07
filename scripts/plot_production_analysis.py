#!/usr/bin/env python3
"""
4-panel summary plot for the EGFR-curcumin 30 ps production MD run.
Mirrors the Project 2 (curcumin-GO-MD-GROMACS) production_analysis.png style.
All inputs are PBC-corrected (production_pbc.xtc based).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load_xvg(path, cols=None):
    data = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            data.append([float(x) for x in parts])
    arr = np.array(data)
    if cols is not None:
        return arr[:, cols]
    return arr

# --- Load data ---
mindist = load_xvg("mindist.xvg")          # time(ps), min dist (nm)
prot_com = load_xvg("protein_com.xvg")     # time, x, y, z (nm)
curc_com = load_xvg("curcumin_com.xvg")    # time, x, y, z (nm)
rmsd = load_xvg("curcumin_rmsd.xvg")       # time(ps), rmsd (nm)
rg = load_xvg("curcumin_rg.xvg")           # time(ps), Rg, RgX, RgY, RgZ (nm)

t_mindist, d_mindist = mindist[:, 0], mindist[:, 1] * 10   # nm -> A
t_com = prot_com[:, 0]
com_dist = np.linalg.norm(prot_com[:, 1:4] - curc_com[:, 1:4], axis=1)  # nm
t_rmsd, d_rmsd = rmsd[:, 0], rmsd[:, 1] * 10               # nm -> A
t_rg, d_rg = rg[:, 0], rg[:, 1] * 10                       # nm -> A

NAVY = "#1F3B57"
TEAL = "#2E8B8B"
RUST = "#B5651D"
PLUM = "#6A4C93"
GRID = "#E3E7EB"

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
fig.suptitle("EGFR + Curcumin — 30 ps Production MD (PBC-corrected)", fontsize=14, fontweight="bold", color=NAVY)

def style(ax, title, ylabel):
    ax.set_title(title, fontsize=11, fontweight="bold", color=NAVY)
    ax.set_xlabel("Time (ps)", fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

# Panel 1: closest contact distance
ax = axes[0, 0]
ax.plot(t_mindist, d_mindist, color=TEAL, linewidth=1.6)
ax.fill_between(t_mindist, d_mindist, alpha=0.12, color=TEAL)
style(ax, "Closest Contact Distance (Curcumin ↔ Protein)", "Distance (Å)")

# Panel 2: COM distance
ax = axes[0, 1]
ax.plot(t_com, com_dist, color=NAVY, linewidth=1.6)
ax.fill_between(t_com, com_dist, alpha=0.12, color=NAVY)
style(ax, "Center-of-Mass Distance", "Distance (nm)")

# Panel 3: curcumin RMSD
ax = axes[1, 0]
ax.plot(t_rmsd, d_rmsd, color=RUST, linewidth=1.6)
ax.fill_between(t_rmsd, d_rmsd, alpha=0.12, color=RUST)
style(ax, "Curcumin Internal RMSD (vs. docked pose)", "RMSD (Å)")

# Panel 4: curcumin radius of gyration
ax = axes[1, 1]
ax.plot(t_rg, d_rg, color=PLUM, linewidth=1.6)
ax.fill_between(t_rg, d_rg, alpha=0.12, color=PLUM)
style(ax, "Curcumin Radius of Gyration", "Rg (Å)")

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("production_analysis.png", dpi=160, facecolor="white")
print("saved production_analysis.png")
print(f"mindist: {d_mindist.min():.2f}-{d_mindist.max():.2f} A, mean {d_mindist.mean():.2f}+/-{d_mindist.std():.2f}")
print(f"com_dist: {com_dist.min():.3f}-{com_dist.max():.3f} nm, mean {com_dist.mean():.3f}+/-{com_dist.std():.3f}")
print(f"rmsd: {d_rmsd.min():.2f}-{d_rmsd.max():.2f} A, mean {d_rmsd.mean():.2f}+/-{d_rmsd.std():.2f}")
print(f"rg: {d_rg.min():.2f}-{d_rg.max():.2f} A, mean {d_rg.mean():.2f}+/-{d_rg.std():.2f}")
