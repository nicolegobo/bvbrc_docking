import logging
import subprocess


def build_logger(debug=0):
    logger_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=logger_level, format="%(asctime)s %(message)s")
    logger = logging.getLogger(__name__)
    return logger


def run_and_save(cmd, output_file=None):
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=output_file,
        stderr=subprocess.STDOUT,
    )
    return process
