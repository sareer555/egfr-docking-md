from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds

# Reference = the exact atom order used in curcumin.itp (built from the Project 1 DFT geometry)
ref = Chem.MolFromXYZFile("/home/claude/curcumin_dft.xyz")
rdDetermineBonds.DetermineBonds(ref, charge=0)

# Query = our docked pose, now with all 47 atoms restored
docked = Chem.MolFromPDBFile("curcumin_pose_allH.pdb", removeHs=False)

print("ref atoms:", ref.GetNumAtoms(), "docked atoms:", docked.GetNumAtoms())

match = docked.GetSubstructMatch(ref)
if not match:
    # try the other direction
    match = ref.GetSubstructMatch(docked)
    print("Matched ref-as-query, len:", len(match))
else:
    print("Matched docked-as-query, len:", len(match))

print(match)

# Reorder docked's 3D coordinates into ref (itp) atom order
conf = docked.GetConformer()
positions_in_ref_order = [conf.GetAtomPosition(match[i]) for i in range(len(match))]

# Pull the itp's atom names (C1, O2, C3 ... in exact order) straight from curcumin.itp
names = []
with open("/home/claude/curcumin.itp") as f:
    in_atoms = False
    for line in f:
        if line.strip().startswith("[ atoms ]"):
            in_atoms = True
            continue
        if in_atoms:
            if line.strip().startswith("["):
                break
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                names.append(parts[4])

print("itp atom names, first 10:", names[:10])
assert len(names) == 47

# Write a PDB with itp atom order/names but the DOCKED 3D positions (in Angstrom -> PDB uses Angstrom already)
with open("curcumin_docked_itporder.pdb", "w") as f:
    for i, (nm, pos) in enumerate(zip(names, positions_in_ref_order), start=1):
        f.write(f"HETATM{i:5d}  {nm:<3s} CURC A   1    {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00          {nm[0]:>2s}\n")
    f.write("END\n")
print("wrote curcumin_docked_itporder.pdb")
