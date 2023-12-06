import glob
import os

import pandas as pd

from .utils import clean_pdb, comb_pdb, run_and_save, sdf2pdb


class diff_dock(object):
    def __init__(
        self, pdb_files, lig_dbs, work_dir, output_path, top_n: int = 1
    ) -> None:
        if isinstance(pdb_files, str):
            self.pdb_files = [pdb_files]
        else:
            self.pdb_files = pdb_files
        self.lig_dbs = lig_dbs
        self.work_dir = work_dir
        self.output_path = os.path.abspath(output_path)

        self.top_n = top_n
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

    def prepare_inputs(self):
        with open(self.lig_dbs, "r") as fp:
            smiless = [line.split()[0] for line in fp]

        self.all_runs = f"{self.output_path}/all.csv"
        with open(self.all_runs, "w") as out:
            out.write("protein_path,ligand\n")
            for pdb_file in self.pdb_files:
                output_pdb = os.path.join(self.output_path, os.path.basename(pdb_file))
                pdb_file = clean_pdb(pdb_file, output_pdb)
                for smiles in smiless:
                    out.write(f"{pdb_file},{smiles}\n")

    def run(self):
        self.prepare_inputs()
        self.get_esm_embeddings()
        self.run_docking()
        self.post_process()

    def get_esm_embeddings(self):
        cmd_getSeq = (
            f"python {self.work_dir}/datasets/esm_embedding_preparation.py "
            f"--protein_ligand_csv {self.all_runs} "
            f"--out_file data/prepared_for_esm.fasta"
        )
        proc = run_and_save(cmd_getSeq, cwd=self.work_dir)
        proc.wait()

        cmd_esm = (
            f"python {self.work_dir}/esm/scripts/extract.py "
            f"esm2_t33_650M_UR50D data/prepared_for_esm.fasta data/esm2_output "
            f"--repr_layers 33 --include per_tok --truncation_seq_length 30000"
        )
        proc = run_and_save(cmd_esm, cwd=self.work_dir)
        proc.wait()

    def run_docking(self):
        cmd_diffdock = (
            f"python -m inference "
            f"--protein_ligand_csv {self.all_runs} "
            f"--out_dir {self.output_path} "
            f"--inference_steps 20 --samples_per_complex 40 --batch_size 6"
        )
        proc = run_and_save(cmd_diffdock, cwd=self.work_dir)
        proc.wait()

    def post_process(self):
        # result_paths = glob.glob(f"{self.output_path}/index*")
        input_df = pd.read_csv(self.all_runs)
        output_df = []
        for i, row in input_df.iterrows():
            result_path = f"{self.output_path}/index{i}_{row['protein_path'].replace('/', '-')}____{row['ligand']}"
            print(result_path)
            for j in range(self.top_n):
                sdf = glob.glob(
                    f"{glob.escape(result_path)}/rank{j+1}_confidence-*.sdf"
                )[0]
                score = float(os.path.basename(sdf).split("-")[1][:-4])
                local_dict = row.to_dict()
                local_dict["lig_sdf"] = sdf
                local_dict["score"] = score
                local_dict["comp_pdb"] = comb_pdb(
                    local_dict["protein_path"], sdf2pdb(sdf)
                )
                output_df.append(local_dict)

        output_df = pd.DataFrame(output_df)
        output_df.to_csv(f"{self.output_path}/result.csv")
        return output_df
