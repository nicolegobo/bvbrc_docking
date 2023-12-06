import logging
import os
import shlex
import subprocess

import MDAnalysis as mda

three_to_one = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "MSE": "M",  # MSE this is almost the same AA as MET. The sulfur is just replaced by Selen
    "PHE": "F",
    "PRO": "P",
    "PYL": "O",
    "SER": "S",
    "SEC": "U",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "ASX": "B",
    "GLX": "Z",
    "XAA": "X",
    "XLE": "J",
}


def build_logger(debug=0):
    logger_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=logger_level, format="%(asctime)s %(message)s")
    logger = logging.getLogger(__name__)
    return logger


def run_and_save(cmd, cwd=None, output_file=None):
    process = subprocess.Popen(
        cmd,
        shell=True,
        cwd=cwd,
        stdout=output_file,
        stderr=subprocess.STDOUT,
    )
    return process


def pdb2seq(pdb_file):
    mda_u = mda.Universe(pdb_file)
    protein = mda_u.select_atoms("protein")
    seq = ""
    for res in protein.residues:
        if res.resname in three_to_one:
            seq += three_to_one[res.resname]
        else:
            seq += "-"
    return seq


def get_pdblabel(pdb_file) -> str:
    return os.path.basename(pdb_file)[:-4]


def comb_pdb(prot_pdb, lig_pdb, comp_pdb=None):
    prot_u = mda.Universe(prot_pdb)
    lig_u = mda.Universe(lig_pdb)
    merged = mda.Merge(prot_u.atoms, lig_u.atoms)
    if comp_pdb is None:
        comp_pdb = f"{os.path.dirname(lig_pdb)}/{get_pdblabel(prot_pdb)}_{get_pdblabel(lig_pdb)}.pdb"
    merged.atoms.write(comp_pdb)
    return comp_pdb


def sdf2pdb(sdf_file, pdb_file=None):
    if pdb_file is None:
        pdb_file = sdf_file[:-3] + "pdb"
    cmd = f"obabel -isdf {shlex.quote(sdf_file)} -opdb > {shlex.quote(pdb_file)}"
    p = run_and_save(cmd)
    p.wait()
    return pdb_file


def clean_pdb(pdb_file, output_pdb=None):
    mda_u = mda.Universe(pdb_file)
    protein = mda_u.select_atoms("protein")
    protein.write(output_pdb)
    return output_pdb
