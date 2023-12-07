import glob
import os
from typing import Optional

import MDAnalysis as mda
from openeye import oechem

from bvbrc_docking.utils import run_and_save


def oe_convert(input_file, output_file):
    ifs = oechem.oemolistream()
    ifs.open(input_file)

    ofs = oechem.oemolostream()
    ofs.open(output_file)
    for mol in ifs.GetOEMols():
        oechem.OEWriteMolecule(ofs, mol)


def pdb_split(pdb_file):
    fp = open(pdb_file, "r")
    file_name, file_ext = os.path.splitext(pdb_file)
    i = 0
    output_pdbs = []
    for line in fp.readlines():
        if line.startswith("COMPND"):
            # lig_name = line.split()[1]
            file_output = f"{file_name}_{i}{file_ext}"
            output_pdbs.append(file_output)
            i += 1
            fo = open(file_output, "w")
            fo.write(line)
        elif line.startswith("END"):
            fo.write(line)
            fo.close()
        else:
            fo.write(line)

    return output_pdbs


class fred_dock(object):
    """_summary_

    Parameters
    ----------
    receptor_pdb : string
        Path name for the receptor pdb file
    drug_dbs : string
        Path name for the drug database smile file
    """

    def __init__(self, receptor_pdb, drug_dbs, n_cpus: int = 1, fred_path="", **kwargs):
        self.receptor_pdb = os.path.abspath(receptor_pdb)
        self.drug_dbs = os.path.abspath(drug_dbs)
        self.n_cpus = n_cpus
        self.fred_path = fred_path
        self.label = os.path.basename(receptor_pdb).split(".")[0]
        self.hostdir = os.path.abspath("./")
        self.run_dir = f"{self.hostdir}/run_{self.label}"
        os.makedirs(self.run_dir)

    def prepare_receptor(self):
        if self.receptor_pdb.endswith("oedu"):
            self.oe_receptor = self.receptor_pdb
        else:
            spruce_cmd = f"spruce -in {self.receptor_pdb}"
            process = run_and_save(spruce_cmd)
            process.wait()

            spruce_out = glob.glob(f"{self.run_dir}/{self.label.upper()}*.oedu")
            if spruce_out == []:
                raise BaseException("spruce run failed...")

            self.oe_receptor = f"{self.run_dir}/{self.label}.oedu"
            MKreceptor_cmd = f"receptorindu -in {spruce_out[0]} -out {self.oe_receptor}"
            process = run_and_save(MKreceptor_cmd)
            process.wait()

    def prepare_lig(self):
        if self.drug_dbs.endswith("oeb.gz"):
            self.oe_dbs = self.drug_dbs
        else:
            self.oe_dbs = f"{self.run_dir}/{self.label}.oeb.gz"
            MKlig_cmd = f"oeomega classic -in {self.drug_dbs} -out {self.oe_dbs}"
            process = run_and_save(MKlig_cmd)
            process.wait()

    def run_fred(self):
        self.oe_docked = f"{self.run_dir}/{self.label}_docked.oeb.gz"
        fred_exec = "fred"
        if self.n_cpus > 1:
            fred_exec += f" -mpi_np {self.n_cpus}"
        fred_cmd = (
            f"{fred_exec} -receptor {self.oe_receptor} "
            f"-dbase {self.oe_dbs} -docked_molecule_file {self.oe_docked}"
        )
        process = run_and_save(fred_cmd)
        process.wait()

    def prepare_report(self):
        report_cmd = (
            f"docking_report -docked_poses {self.oe_docked} "
            f"-receptor {self.oe_receptor} -report_file {self.run_dir}/{self.label}.pdf"
        )
        process = run_and_save(report_cmd)
        process.wait()

    def prepare_output(self):
        output_pdb = f"{self.run_dir}/{self.label}.pdb"
        oe_convert(self.oe_docked, output_pdb)
        lig_pdbs = pdb_split(output_pdb)

        pro_u = mda.Universe(self.receptor_pdb)
        proteins = pro_u.select_atoms("protein")
        for lig_pdb in lig_pdbs:
            lig_u = mda.Universe(lig_pdb)
            comp_u = mda.Merge(proteins.atoms, lig_u.atoms)
            comp_u.atoms.write(lig_pdb)

    def run(self, n_cpus=1):
        os.chdir(self.run_dir)

        self.prepare_receptor()
        self.prepare_lig()

        self.run_fred(n_cpus=n_cpus)

        self.prepare_report()
        if self.receptor_pdb.endswith("pdb"):
            self.prepare_output()

        os.chdir(self.hostdir)
