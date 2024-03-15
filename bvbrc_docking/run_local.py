import logging
import sys
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
    parser.add_argument("-c", "--config", required=True)
    args = parser.parse_args()
    cfg = WorkflowConfig.from_yaml(args.config)

    if cfg.dock.name == "fred":
        from bvbrc_docking.fred import fred_dock as docking
    elif cfg.dock.name == "diffdock":
        from bvbrc_docking.diffdock import diff_dock as docking
    elif cfg.dock.name == "diffdock_1_1":
        from bvbrc_docking.diffdock_1_1 import diff_dock as docking

    dock = docking(**cfg.dock.dict())
    dock.run()
