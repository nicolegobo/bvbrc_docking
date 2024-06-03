import imp
import re

import pandas as pd
import parsl
from parsl.app.app import python_app
from parsl.channels import LocalChannel
from parsl.config import Config
from parsl.configs.local_threads import config
from parsl.executors import ThreadPoolExecutor
from parsl.providers import LocalProvider
from tqdm import tqdm
from utils import process_line

local_config = Config(
    app_cache=True,
    checkpoint_files=None,
    checkpoint_mode=None,
    checkpoint_period=None,
    executors=(
        ThreadPoolExecutor(
            label="threads",
            max_threads=2,
            storage_access=None,
            thread_name_prefix="",
            working_dir=None,
        ),
    ),
    garbage_collect=True,
    initialize_logging=True,
    internal_tasks_max_threads=10,
    max_idletime=120.0,
    monitoring=None,
    retries=1,
    retry_handler=None,
    run_dir="runinfo",
    strategy="simple",
    usage_tracking=False,
)
parsl.load(local_config)

process_line_parsl = python_app(process_line)


if __name__ == "__main__":
    smi_file = "/nfs/lambda_stor_01/homes/heng.ma/Research/BVBRC/docking/uc4-base.smi"

    smi_df = []
    with open(smi_file, "r") as fp:
        for line in tqdm(fp.readlines()):
            mol = process_line_parsl(line)
            smi_df.append(mol)

    fail_runs = "failed.smi"
    df = []
    with open(fail_runs, "w") as fp:
        for i, mol in enumerate(smi_df):
            try:
                df.append(mol.result())
            except Exception as e:
                print(f"Failed run {mol.task_record['args']}, due to {e}", file=fp)
    # print(smi_df)
    df = pd.DataFrame(df)
    df.to_csv("uc4_base_lambda1.csv")
