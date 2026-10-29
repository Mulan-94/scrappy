# scrappy/logger/config.py
import logging
from pathlib import Path

import rich.logging

PACKAGE_NAME = "scrappy"

LOG_FORMATTER = logging.Formatter(
    datefmt="%d.%m.%Y %H:%M:%S", fmt="%(asctime)s -- %(levelname)s : %(message)s"
)

WLOG_FORMATTER = logging.Formatter(
    datefmt="%d.%m.%Y %H:%M:%S",
    fmt="%(asctime)s -- %(name)s : %(levelname)s - %(message)s",
)


def setup_rich_handler():
    # https://rich.readthedocs.io/en/latest/reference/logging.html#logging
    handler = rich.logging.RichHandler(
        rich_tracebacks=True,
        show_time=True,
        show_path=True,
        show_level=True,
        tracebacks_code_width=85,
        markup=True,
    )
    return handler


def setup_filehandler(filename):
    l_handler = logging.FileHandler(filename, mode="w")
    l_handler.setLevel(logging.WARNING)
    l_handler.setFormatter(WLOG_FORMATTER)
    return l_handler


def setup_streamhandler():
    s_handler = logging.StreamHandler()
    s_handler.setLevel(logging.INFO)
    s_handler.setFormatter(LOG_FORMATTER)
    return s_handler


def configure_logger(out_dir=None, level=None):
    logger = logging.getLogger(PACKAGE_NAME)  # always "scrappy"

    # TODO: Allow level and outdir to be set a module level
    # guard against duplicate handlers on repeated calls
    if logger.handlers:
        return logger

    # file handler
    out_dir = Path(out_dir) if out_dir else Path(".")
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.addHandler(setup_filehandler(out_dir / "scrappy.log"))

    # stream handler
    # logger.addHandler(setup_streamhandler())

    # rich handler
    logger.addHandler(setup_rich_handler())

    # set level
    if level is not None:
        level = level.upper() if isinstance(level, str) else int(level)
        logger.setLevel(level)
        for handler in logger.handlers:
            handler.setLevel(level)
    else:
        logger.setLevel(logging.INFO)

    return logger
