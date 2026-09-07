from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

ref = Chem.MolFromXYZFile("/home/claude/curcumin_dft.xyz")
rdDetermineBonds.DetermineBonds(ref, charge=0)
docked = Chem.MolFromPDBFile("curcumin_pose_allH.pdb", removeHs=False)

match = docked.GetSubstructMatch(ref)

# Check 1: element identity must match at every mapped position
mismatches = 0
for i in range(ref.GetNumAtoms()):
    e_ref = ref.GetAtomWithIdx(i).GetSymbol()
    e_docked = docked.GetAtomWithIdx(match[i]).GetSymbol()
    if e_ref != e_docked:
        mismatches += 1
        print(f"MISMATCH at ref atom {i}: ref={e_ref} vs docked={e_docked}")
print(f"Element check: {mismatches} mismatches out of {ref.GetNumAtoms()} atoms")

# Check 2: bond lengths in the REORDERED docked coordinates, using the itp's own bond list
conf = docked.GetConformer()
positions = [conf.GetAtomPosition(match[i]) for i in range(len(match))]

bonds = []
with open("/home/claude/curcumin.itp") as f:
    in_bonds = False
    for line in f:
        if line.strip().startswith("[ bonds ]"):
            in_bonds = True
            continue
        if in_bonds:
            if line.strip().startswith("[") or not line.strip():
                if line.strip().startswith("["):
                    break
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                bonds.append((int(parts[0])-1, int(parts[1])-1))

import math
bad = 0
lens = []
for a,b in bonds:
    p1, p2 = positions[a], positions[b]
    d = math.dist((p1.x,p1.y,p1.z),(p2.x,p2.y,p2.z))
    lens.append(d)
    if d < 0.8 or d > 1.8:
        bad += 1
        print(f"SUSPECT bond {a+1}-{b+1}: {d:.3f} A")
print(f"Bond length check: {len(bonds)} bonds, {bad} outside normal 0.8-1.8 A range")
print(f"Bond length range: {min(lens):.3f} - {max(lens):.3f} A")
