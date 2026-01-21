import logging
import os
import warnings

LOG_FORMATTER = logging.Formatter(
        datefmt='%d.%m.%Y %H:%M:%S',
        fmt="%(asctime)s -- %(module)s : %(message)s")

WLOG_FORMATTER = logging.Formatter(
        datefmt='%d.%m.%Y %H:%M:%S',
        fmt="%(asctime)s -- %(module)s : %(levelname)s - %(message)s")

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
    # ignore overflow errors, assume these to be mostly flagged data
    warnings.simplefilter("ignore")

    if out_dir is None:
        out_dir = "."
    else:
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
    
    l_handler = setup_filehandler(
        os.path.join(out_dir, "xcrapping.log"))

    s_handler = setup_streamhandler()

    logger = logging.getLogger(__name__)

    logger.addHandler(l_handler)
    logger.addHandler(s_handler)

    if level is not None:
        level = level.upper() if isinstance(level, str) else int(level)
        logger.setLevel(level)
        for handler in logger.handlers:
            handler.setLevel(level)
    else:
        logger.setLevel(logging.INFO)

    
    return logger