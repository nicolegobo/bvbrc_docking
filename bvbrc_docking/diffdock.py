import time
import glob
import os
import sys

import numpy as np
import pandas as pd

from bvbrc_docking.utils import clean_pdb, comb_pdb, run_and_save, sdf2pdb


class diff_dock(object):
    def __init__(
        self, receptor_pdb, drug_dbs, diffdock_dir, output_dir, top_n: int = 1, **kwargs
    ) -> None:
        self.receptor_pdb = receptor_pdb
        self.label = os.path.basename(receptor_pdb).split(".")[0]
        self.drug_dbs = drug_dbs
        self.diffdock_dir = diffdock_dir
        self.output_dir = os.path.abspath(output_dir)
        os.makedirs(self.output_dir, exist_ok=True)

        self.run_dir = f"{self.output_dir}/run_{self.label}"
        os.makedirs(self.run_dir)

        self.top_n = top_n

    def prepare_inputs(self):
        with open(self.drug_dbs, "r") as fp:
            smiless = [line.split()[0] for line in fp]

        self.all_runs = f"{self.run_dir}/all.csv"
        with open(self.all_runs, "w") as out:
            out.write("protein_path,ligand\n")
            output_pdb = os.path.join(self.run_dir, os.path.basename(self.receptor_pdb))
            pdb_file = clean_pdb(self.receptor_pdb, output_pdb)
            for smiles in smiless:
                out.write(f"{pdb_file},{smiles}\n")

    def run(self):
        log_file = f"{self.run_dir}/diffdock_log"
        self.log_handle = open(log_file, "w")

        t1 = time.time()
        self.prepare_inputs()
        t2 = time.time()
        self.get_esm_embeddings()
        t3 = time.time()
        self.run_docking()
        t4 = time.time()
        self.post_process()
        t5 = time.time()
        #
        # Summarize individual docks
        #

        time_prep = t2 - t1
        time_embed = t3 - t2
        time_dock = t4 - t3
        time_post = t5 - t4
        time_total = t5 - t1
        with np.printoptions(precision=1, suppress=True):
            dockinfo = np.load(f"{self.run_dir}/run_times.npy")
            print(f"dock times: {np.array2string(dockinfo)}")
            total_dock = np.sum(dockinfo)
            print(f"prep={time_prep:.1f} embed={time_embed:.1f} dock={time_dock:.1f} dock_gpu={total_dock:.1f} post={time_post:.1f} total={time_total:.1f}")

        self.log_handle.close()
        return { "prep": time_prep,
                 "embed": time_embed,
                 "dock": time_dock,
                 "dock_gpu_total": total_dock,
                 "dock_gpu_per_ligand": numpy.asarray(dockinfo).tolist(),
                 "total": time_total
                 };

    def get_esm_embeddings(self):
        cmd_getSeq = (
            f"python {self.diffdock_dir}/datasets/esm_embedding_preparation.py "
            f"--protein_ligand_csv {self.all_runs} "
            f"--out_file {self.run_dir}/prepared_for_esm.fasta"
        )
        proc = run_and_save(cmd_getSeq, cwd=self.run_dir, output_file=self.log_handle)

        model_def = os.getenv("BVDOC_ESM_MODEL")
        if model_def is None:
            model_def = "esm2_t33_650M_UR50D"
        elif not os.path.exists(model_def):
            print(f"Model path defined as {model_def} but that path does not exist", file=sys.stderr)
            sys.exit(1)
        else:
            print(f"Using model path {model_def}")

        cmd_esm = (
            f"python {self.diffdock_dir}/esm/scripts/extract.py "
            f"{model_def} {self.run_dir}/prepared_for_esm.fasta {self.run_dir}/esm2_output "
            f"--repr_layers 33 --include per_tok --truncation_seq_length 30000"
        )
        proc = run_and_save(cmd_esm, cwd=self.run_dir, output_file=self.log_handle)

    def run_docking(self):
        cmd_diffdock = (
            f"python {self.diffdock_dir}/inference.py "
            f"--protein_ligand_csv {self.all_runs} "
            f"--out_dir {self.run_dir} --esm_embeddings_path {self.run_dir}/esm2_output "
            f"--cache_path {self.run_dir}/cache "
            f"--model_dir {self.diffdock_dir}/workdir/paper_score_model "
            f"--confidence_model_dir {self.diffdock_dir}/workdir/paper_confidence_model "
            f"--inference_steps 20 --samples_per_complex 40 --batch_size 6"
        )
        proc = run_and_save(cmd_diffdock, cwd=self.run_dir, output_file=self.log_handle)

    def post_process(self):
        # result_paths = glob.glob(f"{self.run_dir}/index*")
        input_df = pd.read_csv(self.all_runs)
        output_df = []
        for i, row in input_df.iterrows():
            result_path = f"{self.run_dir}/index{i}_{row['protein_path'].replace('/', '-')}____{row['ligand']}"
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
        output_df.to_csv(f"{self.run_dir}/result.csv")
        return output_df
