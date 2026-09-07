from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

ref = Chem.MolFromXYZFile("/home/claude/curcumin_dft.xyz")
rdDetermineBonds.DetermineBonds(ref, charge=0)
docked = Chem.MolFromPDBFile("curcumin_pose_allH.pdb", removeHs=False)
match = docked.GetSubstructMatch(ref)
conf = docked.GetConformer()
positions = [conf.GetAtomPosition(match[i]) for i in range(len(match))]

names = []
with open("/home/claude/curcumin.itp") as f:
    in_atoms = False
    for line in f:
        if line.strip().startswith("[ atoms ]"):
            in_atoms = True; continue
        if in_atoms:
            if line.strip().startswith("["): break
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                names.append(parts[4])

# Read the protein .gro file GROMACS just built
with open("protein_processed.gro") as f:
    lines = f.readlines()

title = lines[0]
natoms_protein = int(lines[1].strip())
protein_atom_lines = lines[2:2+natoms_protein]
box_line = lines[2+natoms_protein]

# Build curcumin's .gro atom lines (PDB is in Angstrom -> .gro needs nm, so divide by 10)
resnum = 313  # next residue number after the protein's 312
curc_lines = []
for i, (nm, pos) in enumerate(zip(names, positions), start=1):
    atomnum = natoms_protein + i
    x, y, z = pos.x/10.0, pos.y/10.0, pos.z/10.0
    curc_lines.append(f"{resnum:5d}{'CURC':<5s}{nm:>5s}{atomnum:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")

total_atoms = natoms_protein + len(names)
with open("complex.gro", "w") as f:
    f.write(title)
    f.write(f"{total_atoms:5d}\n")
    f.writelines(protein_atom_lines)
    f.writelines(curc_lines)
    f.write(box_line)

print(f"protein atoms: {natoms_protein}, curcumin atoms: {len(names)}, total: {total_atoms}")
print("wrote complex.gro")
