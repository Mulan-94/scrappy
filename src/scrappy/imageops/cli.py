import argparse
import typer

from types import SimpleNamespace
from typing import Optional, Annotated, List
from pathlib import Path

sm_app = typer.Typer()
# plot_beam_app = typer.Typer()
# new_fit_app = typer.Typer()




def simple_mask_parser():
    parse = argparse.ArgumentParser()
    parse.add_argument("iname",
        help="Input image from which to generate mask")
    parse.add_argument("-o", "--outname", dest="oname",
        default=[], metavar="", type=str, nargs="+",
        help="Where to dump the output(s). We can iterate over multiple reg files. See '-rb' ")
    parse.add_argument("-nodata", "--nodata", action='store_true',
        help="Exclude the data and just make a mask with the previded region mask")
    parse.add_argument("-above", "--above", dest="above", default=None,
        metavar="", type=float,
        help="Include everything above this value. ie. > above. Aka the lower limit")
    parse.add_argument("-below", "--below", dest="below", default=None,
        metavar="", type=float,
        help="Include everything below this value. i.e < below. Aka the upper limit")
    parse.add_argument("-rb", "--region-bound", dest="regname", default=[],
        metavar="", type=str, nargs="+",
        help="DS9 region file(s) within which to make our mask")
    parse.add_argument("-er", "--exclude-region", dest="ex_regname", default=None,
        metavar="", type=str,
        help="DS9 region file(s) containing the regions that should be excluded")
    return parse


@sm_app.command()
def simple_mask_cli(
    iname: Annotated[Path,
        typer.Argument(help="Input image from which to generate mask")],
    oname: Annotated[Optional[str],
        typer.Option("-o", "--outname", show_default=False,
        help="""Where to dump the output(s)""")] = None,
    nodata: Annotated[Optional[bool],
        typer.Option("-nodata", "--nodata",
        help="Exclude the data and just make a mask with the provided region mask")] = None,
    above: Annotated[Optional[float],
        typer.Option("-above", "--above", show_default=False,
        help="Include everything above this value. ie. > above. Aka the lower limit")] = None,
    below: Annotated[Optional[float],
        typer.Option("-below", "--below", show_default=False,
        help="Include everything below this value. i.e < below. Aka the upper limit")] = None,
    regname: Annotated[Optional[List[Path]],
        typer.Option("-rb", "--region-bound", show_default=False,
        help="""DS9 region file(s) within which to make our mask. Use this
            option option multiple times to specify files. E.g.
            -rb file1.reg -rb file2.reg ...""")] = None,
    ex_regname: Annotated[Optional[List[Path]],
        typer.Option("-er", "--exclude-region", show_default=False,
        help="""DS9 region file(s) containing the regions that should be excluded.
            See -rb for instructions on multiple files.""")] = None
):
    """Generate a simple FITS mask using the given input image 'iname' as a base"""
    opts = SimpleNamespace(**locals())

    from scrappy.imageops.simple_mask import main
    main(opts)
    return




def simple_mask_console():
    sm_app()



def plot_bmaj_bmin_parser():
    parser = argparse.ArgumentParser()


    parser.add_argument("idir", metavar="", type=str,
        help="Directory where to find the input images"
    )

    parser.add_argument("-o", "--output", dest="output", metavar="",
        type=str, default=".",
        help="Directory where to dump outputs. Default is current directory."
    )
    parser.add_argument("-search", "--search", type=str,
        help="Regex of what string to search in the directory",
        default="*[0-9q]-I-image.fits")

    parser.add_argument("-t", "--threshold", dest="threshold", metavar="",
        type=float, default=0.5,
        help="""Channels with WSUM/WSUM.max() below this value will be excluded.
        Only used if -s is active"""
    )
    parser.add_argument("-sc", "--select-channels", dest="sel", default=None,
        type=str, help="""Custom channel selections. This can be a file
        containing those channels, or a comma separated list
        e.g. '0,1,2,3'"""
        )
    parser.add_argument("-as", "--auto-select", action="store_true", dest="auto_select",
        help="Try and suggest valid channels for selection. Use in conjuction with '-t'.")

    parser.add_argument("-pre", "--prefix", type=str, default="00",
        help="Prefix to append on the output images")

    parser.add_argument("-ec", "--exclude-channels", dest="exclude", type=str,
        help="""Channels to be excluded during autoselection.
        This can be a file containing those channels, or a comma separated list
        e.g. '0,1,2,3'.""")
    return parser


def plot_bmaj_bmin_console():
    opts = plot_bmaj_bmin_parser().parse_args()
    from scrappy.imageops.plot_bmaj_bmin import main
    main(opts)


def new_fit_parser():
    parser = argparse.ArgumentParser(
        description="Fit a simple taylor polynomial on Stokes I cube, generating a model I")
    parser.add_argument("cube", type=str,
        help="Name of the input I cube to be modelled")
    parser.add_argument("freqs", type=str,
        help="Text file containing respective frequencies. The freqs should be in Hz")

    parser.add_argument("-mask", "--mask-name", dest="mask", type=str, default=None,
        help="Region within which to perform fits. This is a FITS MASK")

    parser.add_argument("-ex", "--exclude-channels", dest="exclude", type=str, default=None,
        help="""Channels that should be excluded while fitting the image cube.
        This should be a text file containing those channel numbers. Data
        corresponding to the specified channels will be set to zero.""")

    parser.add_argument("-u", "--unstack", dest="unstack", action="store_true",
        help="""Unstack the output model cube as single channelised images. If -ex was specified
        those channels will also be excluded.""")

    parser.add_argument("-deg", "--degree", dest="deg", type=int, default=2,
        help="Degree of the taylor polynimial to be fitted.")

    parser.add_argument("-o", "--output-dir", dest="odir", type=str, default="fitness",
        help="Directory where to dump the output files.")

    return parser


def new_fit_console():
    opts = new_fit_parser().parse_args()
    from scrappy.imageops.new_fit import main
    main(opts)
