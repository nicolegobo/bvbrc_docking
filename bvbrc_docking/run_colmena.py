import logging
import sys
from argparse import ArgumentParser
from functools import partial, update_wrapper
from pathlib import Path
from typing import List

from colmena.models import Result
from colmena.queue.python import PipeQueues
from colmena.task_server import ParslTaskServer
from colmena.thinker import BaseThinker, agent, result_processor
from proxystore.store import register_store
from proxystore.store.file import FileStore
from pydantic import Field

from bvbrc_docking.config import DockConfig
from bvbrc_docking.parsl import ComputeSettingsTypes
from bvbrc_docking.utils import BaseModel, mkdir_validator, path_validator


def run_docking(pdb_file, output_dir, dock_type: str = "fred", **kwargs):
    if "fred" in dock_type:
        from bvbrc_docking.fred import fred_dock as docking
    elif "diffdock" in dock_type:
        from bvbrc_docking.diffdock import diff_dock as docking
    dock = docking(receptor_pdb=pdb_file, output_dir=output_dir, **kwargs)
    dock.run()


class ResultLogger:
    def __init__(self, result_dir: Path) -> None:
        result_dir.mkdir(exist_ok=True)
        self.result_dir = result_dir

    def log(self, result: Result, topic: str) -> None:
        """Write a JSON result per line of the output file."""
        with open(self.result_dir / f"{topic}.json", "a") as f:
            print(result.json(exclude={"inputs", "value"}), file=f)


class Thinker(BaseThinker):
    def __init__(
        self,
        proteins: List[Path],
        result_logger: ResultLogger,
        n_workers: int = 1,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.proteins = proteins
        self.n_workers = n_workers
        self.task_idx = 0
        self.result_logger = result_logger

    def submit_docking_task(self):
        if self.task_idx >= len(self.proteins):
            self.done.set()
            return

        pdb_file = self.proteins[self.task_idx]

        self.queues.send_inputs(
            pdb_file,
            method="run_docking",
            topic="docking",
            keep_inputs=False,
        )
        self.task_idx += 1

    @agent(startup=True)
    def start_tasks(self) -> None:
        for i in range(self.n_workers):
            self.submit_docking_task()

    @result_processor(topic="docking")
    def process_esm_result(self, result: Result):
        self.result_logger.log(result, "docking")
        if not result.success:
            logging.warning(f"Bad inference result: {result.json()}")

        # The old task is finished, start a new one
        self.submit_docking_task()


class WorkflowConfig(BaseModel):
    input_dir: Path = Field(..., description="Path to protein pdbs")
    output_dir: Path = Field(..., description="Path to the workflow output dir")
    n_workers: int = Field(..., description="Number of workers")
    dock: DockConfig = Field(..., description="Docking config.")
    compute_settings: ComputeSettingsTypes

    # validators
    _output_dir_mkdir = mkdir_validator("output_dir")
    _protein_interactions_dir_exists = path_validator("input_dir")

    def configure_logging(self) -> None:
        """Set up logging."""
        logging.basicConfig(
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            level=logging.INFO,
            handlers=[
                logging.FileHandler(self.output_dir / "runtime.log"),
                logging.StreamHandler(sys.stdout),
            ],
        )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", required=True)
    args = parser.parse_args()
    cfg = WorkflowConfig.from_yaml(args.config)
    cfg.configure_logging()

    # Make the proxy store
    store = FileStore(name="file", store_dir=str(cfg.output_dir / "proxy-store"))
    register_store(store)

    # Make the queues
    queues = PipeQueues(
        serialization_method="pickle",
        topics=["docking"],
        proxystore_name="file",
        proxystore_threshold=10000,
    )

    # Define the parsl configuration (this can be done using the get_config
    # for common use cases or by defining your own configuration.)
    parsl_config = cfg.compute_settings.get_config(f"{cfg.output_dir}/run-info")

    dock_output = cfg.output_dir / "dock_result"
    dock_output.mkdir()

    my_run_docking = partial(
        run_docking,
        output_dir=dock_output,
        dock_type=cfg.dock.name,
        **cfg.dock.dict(),
    )
    update_wrapper(my_run_docking, run_docking)

    doer = ParslTaskServer([my_run_docking], queues, parsl_config)

    # Create the result logger
    result_logger = ResultLogger(cfg.output_dir / "result")

    logging.info("Loading proteins")
    proteins = sorted(list(cfg.input_dir.glob("1A*.pdb")))
    thinker = Thinker(
        queue=queues,
        proteins=proteins,
        n_workers=cfg.n_workers,
        result_logger=result_logger,
    )
    logging.info("Created the task server and task generator")

    try:
        # Launch the servers
        doer.start()
        thinker.start()
        logging.info("Launched the servers")

        # Wait for the task generator to complete
        thinker.join()
        logging.info("Task generator has completed")
    finally:
        # Send the kill signal to the task server
        queues.send_kill_signal()

        # Wait for the task server to complete
        doer.join()

        # Clean up proxy store
        store.close()

    logging.info("Workflow has completed")
