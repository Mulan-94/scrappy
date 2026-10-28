import logging
import os
import shutil
import sys
from contextlib import contextmanager
from importlib import util
from pathlib import Path
from time import perf_counter

import numpy as np

snitch = logging.getLogger(__name__)


@contextmanager
def super_write(path, mode="r", **kwargs):
    path = Path(path)
    with open(path, mode, **kwargs) as file:
        yield file
    if "a" in mode or "w" in mode:
        print(f"Saved at: {path.resolve()}")


def np_savetxt(path, data, **kwargs):
    path = Path(path)
    np.savetxt(path, data, **kwargs)
    print(f"Saved at: {path.resolve()}")


def np_savez(path, **kwargs):
    path = Path(path)
    np.savez(path, **kwargs)
    print(f"Saved at: {path.resolve()}")


def savefig(fig, path, **kwargs):
    path = Path(path)
    fig.savefig(path, **kwargs)
    print(f"Saved at: {path.resolve()}")


def example_context(fname, mode="r"):
    manager = open(fname, mode)
    fobj = manager.__enter__()
    try:
        data = fobj.read()
    finally:
        manager.__exit__()
    return data


def append_extension(name: str, extension: str):
    extensions = "fits html json npz png reg txt".split()
    extensions = {ext: f".{ext}" for ext in extensions}
    if extension not in name:
        name += extensions[extension]
    return name


def make_out_dir(dir_name: str, delete_if_exists: bool = False):
    if os.path.isdir(dir_name) and delete_if_exists:
        snitch.info(f"Deleting '{dir_name}' because it already exists")
        shutil.rmtree(dir_name)

    if not os.path.isdir(dir_name):
        snitch.info(f"Creating directory: {dir_name}")
        os.makedirs(dir_name)
    return os.path.abspath(dir_name)


def fullpath(*args):  # TODO: perhaps change the name of this function?
    return os.path.join(*args)


def read_sorted_filnames(fname: Path) -> list[Path]:
    with open(fname) as post:
        items = post.readlines()
        items = [_.replace("\n", "") for _ in items]
    return items


def timer(func):
    def wrapper(*args, **kwargs):
        start = perf_counter()
        result = func(*args, **kwargs)
        print(f"'{func.__name__}' run in: {perf_counter() - start:.2f} sec")
        return result

    return wrapper


def read_npz(filename: Path):
    with np.load(filename, allow_pickle=True) as dat:
        datas = dict(dat)
    return datas


def make_importable():
    import os
    import sys

    PATH = set(sys.path)
    PROJECT_ROOT = os.path.abspath(
        os.path.join(os.path.dirname(__file__), os.pardir)
    )
    if not PATH.issuperset(PROJECT_ROOT):
        sys.path.append(PROJECT_ROOT)


def does_specified_file_exist(*args):
    not_found = []
    for fname in args:
        if fname is None:
            continue
        if not os.path.isfile(fname):
            not_found.append(fname)

    if len(not_found) > 0:
        snitch.error("The following files/dirs are required but were not found:")
        for fname in not_found:
            snitch.error(f"-> {fname}")
        snitch.error("Please ensure that they exist if intended for use")
        sys.exit()
    else:
        return True


def lazy_loader(name):
    # Check if already in cache
    if name in sys.modules:
        return sys.modules[name]

    spec = util.find_spec(name)
    if spec is None:
        snitch.error(f"Module {name} was not found. Ensure it is installed")
        raise ImportError(f"Required module '{name}' missing")

    # 2. Create a lazy loader based on the existing loader
    loader = util.LazyLoader(spec.loader)
    spec.loader = loader

    # 3. Create the module object from the spec
    module = util.module_from_spec(spec)

    # 4. Insert into sys.modules and execute the loader
    sys.modules[name] = module
    loader.exec_module(module)

    snitch.debug(f"Lazily importing heaviness: {name}")

    return module


# from https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
class dicto(dict):
    """Use dot operator to access to dictionary attributes"""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
