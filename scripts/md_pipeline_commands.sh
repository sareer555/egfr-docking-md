#!/usr/bin/env bash
# EGFR (1M17) + Curcumin — GROMACS protein-ligand MD pipeline
# Force field: AMBER99SB-ILDN (protein), TIP3P (water), custom curcumin.itp (ligand)
set -e

# --- 1. Fix missing loop (residues 965-976 disordered in the 1M17 crystal) ---
# Insert a TER record between residue 964 and 977 so pdb2gmx treats it as a
# genuine chain break instead of bonding across the gap (avoids a "Long Bond" warning).

# --- 2. Protein topology (real validated force field, unlike Project 2's by-analogy fallback) ---
gmx pdb2gmx -f receptor_protein_fixed.pdb -o protein_processed.gro \
    -ff amber99sb-ildn -water tip3p -ignh

# --- 3. Combine protein + curcumin's docked pose into one system ---
python3 build_complex_gro.py     # protein_processed.gro (5032 atoms) + curcumin (47 atoms, reordered/repositioned) -> complex.gro

# curcumin.itp: renamed custom "H" atomtype -> "H_curc" to avoid colliding with
# AMBER99SB-ILDN's own "H" atomtype (a silent global-scope conflict — see README).

# --- 4. Box + solvation + neutralization ---
gmx editconf -f complex.gro -o complex_boxed.gro -c -d 1.2 -bt dodecahedron
gmx solvate -cp complex_boxed.gro -cs spc216.gro -o complex_solv.gro -p topol.top
gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s ions.tpr -o complex_ions.gro -p topol.top -pname NA -nname CL -neutral

# --- 5. Custom index group (avoids atom overlap between default T-coupling groups) ---
gmx make_ndx -f complex_ions.gro -o index.ndx <<EOF
1 | 13
q
EOF
# group 21 = "Protein_CURC" (5079 atoms, no ions)

# --- 6. Energy minimization ---
gmx grompp -f minim.mdp -c complex_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em
# Result: converged in 1355 steps, Fmax 983.9

# --- 7. NVT equilibration (position-restrained, 10 ps) ---
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt
# Result: 297.8 K

# --- 8. NPT equilibration (10 ps) ---
# npt.mdp includes refcoord_scaling = com (required with position restraints + pressure coupling)
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt
# Result: 300.1 K, -36.8 bar

# --- 9. Production (30 ps, reduced scope vs. Project 2's 500 ps — see README) ---
gmx grompp -f production.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o production.tpr
gmx mdrun -deffnm production -cpt 1

# --- 10. PBC correction (required before ANY spatial-distance analysis) ---
gmx trjconv -f production.xtc -s production.tpr -n index.ndx \
    -pbc mol -center -o production_pbc.xtc
# (choose "Protein" to center on, "System" to output)

# --- 11. Analysis ---
gmx mindist  -f production_pbc.xtc -s production.tpr -n index.ndx -o mindist.xvg          # closest contact
gmx traj     -f production_pbc.xtc -s production.tpr -n index.ndx -com -ox protein_com.xvg  # (group: Protein)
gmx traj     -f production_pbc.xtc -s production.tpr -n index.ndx -com -ox curcumin_com.xvg # (group: CURC)
gmx rms      -f production_pbc.xtc -s production.tpr -n index.ndx -o curcumin_rmsd.xvg      # (group: CURC vs CURC)
gmx gyrate   -f production_pbc.xtc -s production.tpr -n index.ndx -o curcumin_rg.xvg        # (group: CURC)

python3 plot_production_analysis.py   # -> figures/production_analysis.png
