from pathlib import Path
from types import SimpleNamespace
from typing import Annotated, List, Optional

import typer

from scrappy.logging.config import configure_logger

snitch = configure_logger()

help_set = {"context_settings": {"help_option_names": ["-h", "--help"]}}

sm_app = typer.Typer(**help_set)
plot_beam_app = typer.Typer(**help_set)
new_fit_app = typer.Typer(**help_set)

# TODO: configure logger under each *_cli() to enable outdir perhaps
# TODO: do it just once for each module


def validate_inputs(files: list[Path | None]):
    """
    Check that specified input files actually exist before the program starts
    """
    errors = []
    for file in files:
        if file is not None and not file.exists():
            errors.append(f"Input file not found: {file} ")
    if errors:
        for _ in errors:
            snitch.error(_)
        raise typer.Exit(1)


@sm_app.command()
def simple_mask_cli(
    iname: Annotated[
        Path, typer.Argument(help="Input image from which to generate mask")
    ],
    outname: Annotated[
        Optional[str],
        typer.Option(
            "-o",
            "--outname",
            show_default=False,
            help="""Name of output mask""",
        ),
    ] = "sc-output-mask.fits",
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
    outdir: Annotated[
        str,
        typer.Option(
            "-o",
            "--output",
            rich_help_panel="Input/Output",
            help="Directory where to dump outputs. Default is current directory.",
        ),
    ] = "sc-beam-plot-io",
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
    validate_inputs([idir])

    opts = SimpleNamespace(**locals())
    from scrappy.imageops.plot_bmaj_bmin import main

    main(opts)
    return


def plot_bmaj_bmin_console():
    plot_beam_app()


@new_fit_app.command()
def new_fit_cli(
    cube: Annotated[
        Path,
        typer.Option(
            "-cube",
            "--cube",
            rich_help_panel="Required",
            help="Input Stokes I cube to be modelled",
        ),
    ],
    freqs: Annotated[
        Path,
        typer.Option(
            "-freqs",
            "--freqs",
            rich_help_panel="Required",
            help="""A text file containing respective frequencies.
            The frequencies should be in Hz""",
        ),
    ],
    deg: Annotated[
        int,
        typer.Option(
            "-deg", "--degree", help="Degree of the taylor polynimial to be fitted."
        ),
    ] = 2,
    exclude: Annotated[
        str,
        typer.Option(
            "-ex",
            "--exclude-channels",
            help="""A list of channels that should be excluded during fitting.
            Can be a:
            1. Comma-separated list
            2. A range specified by a colon e.g. 23:56
            3. A text file containing channels to be excluded, line by line.
            Data corresponding to the specified channel will be set to Nan""",
        ),
    ] = None,
    mask: Annotated[
        Path,
        typer.Option(
            "-mask", "--mask", help="FITS mask within which to perform fitting"
        ),
    ] = None,
    outdir: Annotated[
        Path,
        typer.Option(
            "-o", "--output-dir", help="Directory in which to dump the output files."
        ),
    ] = "sc-fit-io",
    unstack: Annotated[
        bool,
        typer.Option(
            "-unstack",
            "--unstack",
            help="""Unstack the output model cube as single channelised images.
            If -ex was specified those channels will also be excluded.""",
        ),
    ] = False,
):
    """
    Fit a simple taylor polynomial on Stokes I cube, generating a model I
    """
    validate_inputs([cube, freqs, mask])
    opts = SimpleNamespace(**locals())
    from scrappy.imageops.new_fit import main

    main(opts)
    return


def new_fit_console():
    new_fit_app()
