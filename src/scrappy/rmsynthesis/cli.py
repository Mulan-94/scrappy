from pathlib import Path
from types import SimpleNamespace
from typing import Annotated

import typer

from scrappy.logging.config import configure_logger

snitch = configure_logger()

help_set = {
    "context_settings": {"help_option_names": ["-h", "--help"]},
    "rich_markup_mode": "markdown",
}
rm_maps_app = typer.Typer(**help_set)
los_rm_app = typer.Typer(**help_set)


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


@rm_maps_app.command()
def rm_maps_cli(
    icube: Annotated[
        Path,
        typer.Option(
            "-i", "--icube", rich_help_panel="Required", help="Stokes I cube (fits)"
        ),
    ],
    qcube: Annotated[
        Path,
        typer.Option(
            "-q", "--qcube", rich_help_panel="Required", help="Stokes Q cube (fits)"
        ),
    ],
    ucube: Annotated[
        Path,
        typer.Option(
            "-u", "--ucube", rich_help_panel="Required", help="Stokes U cube (fits)"
        ),
    ],
    freqs: Annotated[
        Path,
        typer.Option(
            "-f",
            "--frequencies",
            rich_help_panel="Required",
            help="A text file containing the cubes frequencies",
        ),
    ],
    max_depth: Annotated[
        int,
        typer.Option(
            "-md",
            "--max-depth",
            rich_help_panel="RM-Synthesis",
            help="Maximum Faraday depth to fit for.",
        ),
    ],
    debug: Annotated[
        bool, typer.Option("-debug", "--debug", help="Enable debug mode")
    ] = False,
    mask: Annotated[
        Path,
        typer.Option(
            "-mask", "--mask", help="FITS mask image to use for processing"
        ),
    ] = None,
    numProcessor: Annotated[
        int,
        typer.Option("-ncore", "--ncore", help="Number of computing cores to use."),
    ] = 20,
    outdir: Annotated[
        Path,
        typer.Option(
            "-o", "--outdir", help="Where to deposit outputs from this script."
        ),
    ] = "rmmap-io",
    depth_step: Annotated[
        int,
        typer.Option(
            "--depth-step", rich_help_panel="RM-Synthesis", help="Depth stepping."
        ),
    ] = 1,
    niters: Annotated[
        int,
        typer.Option(
            "-niters",
            "--niters",
            rich_help_panel="RM-Synthesis",
            help="Number of clean iterations.",
        ),
    ] = 1000,
    no_derotate: Annotated[
        bool,
        typer.Option(
            "-nd",
            "--no-derotate",
            rich_help_panel="RM-Synthesis",
            help="Use this switch to NOT derotate the RMTF by the mean lambda squared.",
        ),
    ] = False,
    snr: Annotated[
        float,
        typer.Option(
            "-snr",
            "--snr-threshold",
            help="Signal to NOise ratio below which data is excluded",
        ),
    ] = 10,
):
    """
    Performs linear least squares fitting to Q and U image cubes.

    ### Outputs
    1. Fractional polarization (fpol.fits)
    2. Linear polarization (lpol.fits)
    3. Rotation measure (rm.fits)
    4. Polarization angle map (pa.fits)
    5. Intrinsic polarization angle (po.fits)
    6. Magnetic field LIC map (mag-field.png)
    7. SNR (snr.fits)
    """
    in_files = [icube, qcube, ucube, freqs, mask]
    validate_inputs(in_files)
    opts = SimpleNamespace(**locals())
    from scrappy.rmsynthesis.rmmap_x2 import main

    main(opts)
    return


def rm_maps_console():
    rm_maps_app()


@los_rm_app.command()
def per_los_rm_static_cli(
    indir: Annotated[
        Path,
        typer.Option(
            "-id",
            "--input-dir",
            help="Directory containing the various LoS data files (reg files)",
        ),
    ],
    max_fdepth: Annotated[
        int,
        typer.Option(
            "-md",
            "--max-fdepth",
            rich_help_panel="RM-synthesis",
            help="Maximum Faraday depth.",
        ),
    ],
    outdir: Annotated[
        Path,
        typer.Option(
            "-o", "--output-dir", help="Where to dump the output plots if available"
        ),
    ] = "sc-rm-io",
    depth_step: Annotated[
        int,
        typer.Option(
            "-ds",
            "--depth-step",
            rich_help_panel="RM-synthesis",
            help="Faraday depth step.",
        ),
    ] = 10,
    niters: Annotated[
        int,
        typer.Option(
            "-iters",
            "--iters",
            rich_help_panel="RM-synthesis",
            help="Number of RM clean iterations",
        ),
    ] = 100,
    plot: Annotated[
        bool,
        typer.Option("-nop", "--no-plot", help="Switch plotting off for this data"),
    ] = False,
    no_derotate: Annotated[
        bool,
        typer.Option(
            "-nd",
            "--no-derotate",
            rich_help_panel="RM-synthesis",
            help="Use this switch to NOT derotate the RMTF by the mean lambda squared.",
        ),
    ] = True,
    debug: Annotated[
        bool,
        typer.Option(
            "-debug", "--debug", help="Enable debug mode, will run in serial mode"
        ),
    ] = False,
):
    """
    This script takes in I, Q and U data and does the process
    of RM-Synthesis and RM-CLEAN. It gives the following outputs

    ### Pickle files containing:
    1. The dirty FDF (Faraday Dispersion Funtion / Faraday spectrum)
    2. The cleaned FDF
    3. Faraday depths used. These numpy dict keys used are:
        - 'depths'
        - 'fdirty'
        - 'fclean'
        - 'rmtf'

    ### Plots of:
    1. The dirty and clean FDF and position angle vs wavelength sq and its
    linear squares fit
    2. The data's RMSF
    """

    in_files = [indir]
    validate_inputs(in_files)
    opts = SimpleNamespace(**locals())
    from scrappy.rmsynthesis.rm_synthesis import main

    main(opts)
    return


def per_los_rm_static_console():
    los_rm_app()
