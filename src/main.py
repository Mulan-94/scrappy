import argparse
import configparser
import os
from shutil import copy

from argparse import RawTextHelpFormatter

__version__ = "0.1.0"

def get_console_scripts():
    config = configparser.ConfigParser()
    config.read_file(open(os.path.join(os.path.dirname(__file__), "..", r"setup.cfg")))

    scripts = dict([
        tuple(_.split("=")) for _ in 
        config.get("options.entry_points", "console_scripts").split("\n") if len(_)>0])
    
    return list(scripts.keys())


def parser():
    parser = argparse.ArgumentParser(add_help=False,  formatter_class=RawTextHelpFormatter)
    parser.add_argument('-v', '--version', action='version',
        version='%(prog)s' + f" {__version__}", 
        help="Show %(prog)s version number and exit") 
    parser.add_argument('-h', "--help", action="help",
    help=f"""
This is a collection of python scripts that perform some RM related
tasks, and geared towards (but not limited to) generating files to be used
with polarvis: https://github.com/Mulan-94/polarvis\n""" + \
"The following command line tools are available and help accesible via: \033[1m" + \
" ".join([f'{_.strip()} -h\n' for _ in get_console_scripts()]) + "\033[0m"
        )

    parser.add_argument("-i", "--initialize", action="store_true", dest="initialize",
    help="Set up files required to run the showrunner. It will not run without those files")

    parser.add_argument("-f", "--force", action="store_true", dest="overwrite",
    help="""Force overwrite of the files required for showrunner in case they 
    already exist in the current directory."""
    )

    parser.add_argument("-p", "--pica", action="store_true", dest="pica",
    help="Initialise some pica files specifically")
    return parser


def initialize(pica=False, overwrite=False):
    op_dir = os.path.join(os.path.dirname(__file__),  "..", "posta")
    work_dir = os.path.abspath(os.path.curdir)

    pica_files = ["pica_box_region.reg", "pica_noise_region.reg"]

    default_files = ["showrunner.sh", "env-vars", "bk_plots.yml"]

    if pica:
        default_files.extend(pica_files)
    
    
    print("\n")
    for fil in default_files:
        if not os.path.exists(os.path.join(work_dir, fil)) or overwrite:
            print(f"{fil:20} : not found. Copying to {work_dir}")
            copy(os.path.join(op_dir, fil), os.path.join(work_dir, fil))
        else:
            print(f"""{fil:20} : already exists, skipping. 
            Use the '-f' option to enforce an overwrite,""")

def console():
    opts = parser().parse_args()
    if opts.initialize:
        initialize(opts.pica, overwrite=opts.overwrite)
    
    return


if __name__ == "__main__":
    console()