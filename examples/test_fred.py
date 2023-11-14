import os
import sys

sys.path.append("../")
from bvbrc_docking import fred_docking

receptor_pdb = "/lambda_stor/homes/heng.ma/Research/BVBRC/docking/structures/1AH5.pdb"
drugs_dbs = (
    "/lambda_stor/homes/heng.ma/Research/BVBRC/docking/oe_fred/example/renin/all.smi"
)

docking = fred_docking(receptor_pdb, drugs_dbs)
docking.run()
