import logging
import os
import sys
import warnings
from argparse import ArgumentParser

from pydantic import Field

from bvbrc_docking.config import DockConfig
from bvbrc_docking.utils import BaseModel, mkdir_validator


class WorkflowConfig(BaseModel):
    output_dir: str = "./"
    dock: DockConfig = Field(..., description="Docking config.")
    # compute_settings: ComputeSettings

    def configure_logging(self) -> None:
        """Set up logging."""
        logging.basicConfig(
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            level=logging.INFO,
            handlers=[
                logging.FileHandler(f"{self.output_dir}/runtime.log"),
                logging.StreamHandler(sys.stdout),
            ],
        )

    # validators
    _output_dir_mkdir = mkdir_validator("output_dir")
    # _protein_interactions_dir_exists = path_validator("protein_interactions_dir")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--config")
    parser.add_argument(
        "-n",
        "--name",
        choices=["diffdock", "diffdock_1_1", "fred"],
        default="diffdock_1_1",
    )
    parser.add_argument("-r", "--receptor-pdb")
    parser.add_argument("-d", "--drug-dbs")
    parser.add_argument(
        "-D", "--diffdock-dir", default=os.getenv("BVDOCK_DIFFDOCK_DIR")
    )
    parser.add_argument("-t", "--top-n", default=3)
    parser.add_argument("output_dir")
    # parser.add_argument("-o", "--output_dir", default="./xout")

    args = parser.parse_args()

    if args.config is None:
        cfg = WorkflowConfig.from_args(args)
    else:
        cfg = WorkflowConfig.from_yaml(args.config)

    cfg.configure_logging()
    logging.info(f"Starting Docking run for {cfg.dock.receptor_pdb}")

    if cfg.dock.name == "fred":
        from bvbrc_docking.fred import fred_dock as docking
    elif cfg.dock.name == "diffdock":
        from bvbrc_docking.diffdock import diff_dock as docking
    elif cfg.dock.name == "diffdock_1_1":
        from bvbrc_docking.diffdock_1_1 import diff_dock as docking

    dock = docking(**cfg.dock.dict())

    warnings.filterwarnings("ignore", ".*guess.*", UserWarning)
    warnings.filterwarnings("ignore", ".*Found no information.*", UserWarning)
    warnings.filterwarnings("ignore", ".*Unit cell dimensions not found.*", UserWarning)
    warnings.filterwarnings("ignore", ".*Found missing chainIDs.*", UserWarning)

    dock.run()
    logging.info(f"Finished Docking run for {cfg.dock.receptor_pdb}")
    #    for w in wlist:
    # print(f"message={w.message} cat={w.category} fn={w.filename} line={w.lineno}")
