import argparse
from pathlib import Path
from types import SimpleNamespace
from typing import Annotated, List, Optional

import typer

sm_app = typer.Typer(context_settings={"help_option_names": ["-h", "--help"]})
plot_beam_app = typer.Typer(context_settings={"help_option_names": ["-h", "--help"]})

# new_fit_app = typer.Typer()


@sm_app.command()
def simple_mask_cli(
    iname: Annotated[
        Path, typer.Argument(help="Input image from which to generate mask")
    ],
    oname: Annotated[
        Optional[str],
        typer.Option(
            "-o",
            "--outname",
            show_default=False,
            help="""Where to dump the output(s)""",
        ),
    ] = None,
    above: Annotated[
        Optional[float],
        typer.Option(
            "-above",
            "--above",
            show_default=False,
            help="Include everything above this value. ie. > above. Aka the lower limit",
        ),
    ] = None,
    below: Annotated[
        Optional[float],
        typer.Option(
            "-below",
            "--below",
            show_default=False,
            help="Include everything below this value. i.e < below. Aka the upper limit",
        ),
    ] = None,
    ex_regname: Annotated[
        Optional[List[Path]],
        typer.Option(
            "-er",
            "--exclude-region",
            show_default=False,
            help="""DS9 region file(s) containing the regions that should be excluded.
            See -rb for instructions on multiple files.""",
        ),
    ] = None,
    nodata: Annotated[
        Optional[bool],
        typer.Option(
            "-nodata",
            "--nodata",
            help="Exclude the data and just make a mask with the provided region mask",
        ),
    ] = None,
    regname: Annotated[
        Optional[List[Path]],
        typer.Option(
            "-rb",
            "--region-bound",
            show_default=False,
            help="""DS9 region file(s) within which to make our mask. Use this
            option option multiple times to specify files. E.g.
            -rb file1.reg -rb file2.reg ...""",
        ),
    ] = None,
):
    """
    Generate a simple FITS mask using the given input image 'iname' as a base
    """
    opts = SimpleNamespace(**locals())

    from scrappy.imageops.simple_mask import main

    main(opts)
    return


def simple_mask_console():
    sm_app()


@plot_beam_app.command()
def plot_bmaj_bmin_cli(
    idir: Annotated[
        Path,
        typer.Option(
            "--idir",
            rich_help_panel="Input/Output",
            help="Directory where to find the input images",
        ),
    ],
    output: Annotated[
        str,
        typer.Option(
            "-o",
            "--output",
            rich_help_panel="Input/Output",
            help="Directory where to dump outputs. Default is current directory.",
        ),
    ],
    auto_select: Annotated[
        Optional[bool],
        typer.Option(
            "-as",
            "--auto-select",
            help="Try and suggest valid channels for selection. Use in conjuction with '-t'.",
        ),
    ] = False,
    exclude: Annotated[
        Optional[str],
        typer.Option(
            "-ec",
            "--exclude-channels",
            help="""Channels to be excluded during autoselection.
            This can be a file containing those channels, or a comma separated list
            e.g. '0,1,2,3'.""",
        ),
    ] = None,
    prefix: Annotated[
        Optional[str],
        typer.Option(
            "-pre",
            "--prefix",
            show_default=True,
            help="Prefix to append on the output images",
        ),
    ] = "00",
    search_str: Annotated[
        Optional[str],
        typer.Option(
            "-search",
            "--search",
            show_default=True,
            help="Regex of what string to search in the directory",
        ),
    ] = "*[0-9q]-I-image.fits",
    sel: Annotated[
        Optional[str],
        typer.Option(
            "-sc",
            "--select-channels",
            help="""Custom channel selections. This can be a file
        containing those channels, or a comma separated list
        e.g. '0,1,2,3'""",
        ),
    ] = None,
    threshold: Annotated[
        Optional[float],
        typer.Option(
            "-t",
            "--threshold",
            show_default=True,
            help="""Channels with WSUM/WSUM.max() below this value will be excluded.
            Only used if -s is active""",
        ),
    ] = 0.5,
):
    """
    Make a simple plot of of BMAJ and BMIN for the various
    specified frequency channels
    """

    opts = SimpleNamespace(**locals())
    from scrappy.imageops.plot_bmaj_bmin import main

    main(opts)
    return


def plot_bmaj_bmin_console():
    plot_beam_app()


def new_fit_parser():
    parser = argparse.ArgumentParser(
        description="Fit a simple taylor polynomial on Stokes I cube, generating a model I"
    )
    parser.add_argument(
        "cube", type=str, help="Name of the input I cube to be modelled"
    )
    parser.add_argument(
        "freqs",
        type=str,
        help="Text file containing respective frequencies. The freqs should be in Hz",
    )

    parser.add_argument(
        "-mask",
        "--mask-name",
        dest="mask",
        type=str,
        default=None,
        help="Region within which to perform fits. This is a FITS MASK",
    )

    parser.add_argument(
        "-ex",
        "--exclude-channels",
        dest="exclude",
        type=str,
        default=None,
        help="""Channels that should be excluded while fitting the image cube.
        This should be a text file containing those channel numbers. Data
        corresponding to the specified channels will be set to zero.""",
    )

    parser.add_argument(
        "-u",
        "--unstack",
        dest="unstack",
        action="store_true",
        help="""Unstack the output model cube as single channelised images. If -ex was specified
        those channels will also be excluded.""",
    )

    parser.add_argument(
        "-deg",
        "--degree",
        dest="deg",
        type=int,
        default=2,
        help="Degree of the taylor polynimial to be fitted.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        dest="odir",
        type=str,
        default="fitness",
        help="Directory where to dump the output files.",
    )

    return parser


def new_fit_console():
    opts = new_fit_parser().parse_args()
    from scrappy.imageops.new_fit import main

    main(opts)
