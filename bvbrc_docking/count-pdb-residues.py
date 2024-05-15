#
# Count residues in the given PDB file
#

import MDAnalysis as mda
from openbabel import pybel
import sys
import tempfile
import os

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} pdb-file", file=sys.stderr)
    sys.exit(1)

pdb_file = sys.argv[1]

if pdb_file.endswith("cif.gz") or pdb_file.endswith("cif"):
    mol = next(pybel.readfile("cif", pdb_file))

    temp_pdb = tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False)
    temp_pdb.write(mol.write(format="pdb"))
    temp_pdb.close()

    mda_u = mda.Universe(temp_pdb.name)
    os.unlink(temp_pdb.name)
else:
    mda_u = mda.Universe(pdb_file)

protein = mda_u.select_atoms("protein")
print(protein.n_residues)

