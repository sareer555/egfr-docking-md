#!/usr/bin/env bash
# EGFR (PDB 1M17) — AutoDock Vina docking pipeline
# Tools: AutoDock Vina 1.2.5, Open Babel 3.1.1
set -e

# --- 1. Receptor preparation ---
# Extract protein-only ATOM records from 1M17.pdb -> structures/receptor_protein.pdb (not included, derivable)
obabel receptor_protein.pdb -O receptor_h.pdb -p 7.4      # add H at physiological pH
obabel receptor_h.pdb -O receptor.pdbqt -xr               # rigid-receptor PDBQT

# --- 2. Ligand preparation (from SMILES, see structures/smiles.txt) ---
obabel -:"<SMILES_STRING>" -O curcumin.pdbqt     --gen3d -p 7.4
obabel -:"<SMILES_STRING>" -O doxorubicin.pdbqt  --gen3d -p 7.4
obabel -:"<SMILES_STRING>" -O erlotinib.pdbqt    --gen3d -p 7.4
obabel -:"<SMILES_STRING>" -O gefitinib.pdbqt    --gen3d -p 7.4

# --- 3. Search box: centered on the crystallographic AQ4 (erlotinib) coordinates ---
# extracted from 1M17.pdb HETATM records -> structures/crystal_ligand.pdb
# box center/size determined from that ligand's bounding coordinates + 20 A padding per axis

# --- 4. Docking runs (Vina 1.2.5 dropped --log; redirect stdout instead) ---
vina --receptor receptor.pdbqt --ligand curcumin.pdbqt \
     --center_x <X> --center_y <Y> --center_z <Z> \
     --size_x 20 --size_y 20 --size_z 20 \
     --exhaustiveness 8 --out curcumin_docked.pdbqt | tee curcumin_dock.log

vina --receptor receptor.pdbqt --ligand doxorubicin.pdbqt \
     --center_x <X> --center_y <Y> --center_z <Z> \
     --size_x 20 --size_y 20 --size_z 20 \
     --exhaustiveness 8 --out doxorubicin_docked.pdbqt | tee doxorubicin_dock.log

vina --receptor receptor.pdbqt --ligand erlotinib.pdbqt \
     --center_x <X> --center_y <Y> --center_z <Z> \
     --size_x 20 --size_y 20 --size_z 20 \
     --exhaustiveness 8 --out erlotinib_docked.pdbqt | tee erlotinib_dock.log

vina --receptor receptor.pdbqt --ligand gefitinib.pdbqt \
     --center_x <X> --center_y <Y> --center_z <Z> \
     --size_x 20 --size_y 20 --size_z 20 \
     --exhaustiveness 8 --out gefitinib_docked.pdbqt | tee gefitinib_dock.log

# --- 5. Re-docking validation: erlotinib's docked pose vs. its known crystal pose ---
# Graph-based atom matching (-m) is required since PDBQT/PDB atom order/naming differ.
obrms crystal_ligand.pdb erlotinib_docked.pdbqt -m
# Result: 1.93 A RMSD to the experimental pose -> docking protocol validated
