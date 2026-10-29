import importlib.metadata
import logging

logging.getLogger("scrappy")

# from scrappy.logging.config import configure_logger

__version__ = importlib.metadata.version("scrappy")
