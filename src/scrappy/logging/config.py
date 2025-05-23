import logging
import logging.config
import os
import warnings
from pathlib import Path

from dataclasses import dataclass

@dataclass
class SnitchSettings:
    LOG_LEVEL: str = "info"
    LOG_DIR: str = "."
    LOG_FILE_NAME: str = "scraplog"

logger_vars = SnitchSettings()


# CONSOLE_LOG_FORMATTER = dict(
#         class="logging.Formatter",
#         datefmt='%H:%M:%S %d.%m.%Y',
#         fmt="%(asctime)s : %(message)s")

# FILE_LOG_FORMATTER = dict(
#     class = logging.Formatter,
#     datefmt='%H:%M:%S %d.%m.%Y',
#     fmt="%(asctime)s : %(levelname)s - %(message)s")


DEFAULT_CONFIG = {
    'version': 1,
    'formatters': {
        'standard': {
            'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
        },
        'consolelog':{
            "class": "logging.Formatter",
            "datefmt": '%H:%M:%S %d.%m.%Y',
            "fmt": "%(asctime)s : %(message)s"
        },
        'filelog':{
            "class" :  "logging.Formatter",
            "datefmt": '%H:%M:%S %d.%m.%Y',
            "fmt": "%(asctime)s : %(levelname)s - %(message)s"
        }
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'level': 'INFO',
            'formatter': "consolelog",
            'stream': 'ext://sys.stdout',
        },
        'file': {
            'class': 'logging.handlers.RotatingFileHandler',
            'level': 'WARNING',
            'formatter': "filelog",
            'filename': Path(os.environ.get(
                logger_vars.LOG_DIR, '.')) / f"{logger_vars.LOG_FILE_NAME}.log",
            'maxBytes': 10485760,  # 10MB
            'backupCount': 5,
        }
    },
    'loggers': {
        'scrappy': {
            'level': os.environ.get('LOG_LEVEL', 'INFO'),
            'handlers': ['console', 'file'],
            'propagate': False
        }
    }
}

def setup_logging():
    """Configure logging for the package."""
    logging.config.dictConfig(DEFAULT_CONFIG)
    return logging.getLogger('scrappy')