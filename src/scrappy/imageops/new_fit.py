import matplotlib
import numpy as np

matplotlib.use("Agg")
import sys
from copy import deepcopy
from pathlib import Path

from scrappy.imageops.cli import new_fit_app
from scrappy.logging.config import configure_logger
from scrappy.utils.astroimage import AstroImage
from scrappy.utils.genutils import save_fits, selector

snitch = configure_logger()
# snitch = logging.getLogger(__name__)

OUTDIR = None


def polynomial_fit(x, y_raw, degree):
    """
    Perform a polynomial fit on the input data
    Not that neither should contain Nan values. Zeroes are also discouraged
    Unless they form part of the raw data.

    parameters
    ----------
    x: np.array
        the independent variable (i.e. frequency for this specific script)
    y_raw: np.array
        the raw data being fitted
    degree: int
        Degree of poynomial fit

    Returns
    -------
    y_fit: np.array
        The fitted model to y
    coeff: tuple
        Coefficients to the fitted polynomial
    """
    if hasattr(y_raw, "mask"):
        coeffs = np.ma.polyfit(x, y_raw, degree)
    else:
        coeffs = np.polyfit(x, y_raw, degree)
    polynom = np.poly1d(coeffs)
    y_fit = polynom(x)
    return y_fit, tuple(coeffs.data)


def split_cube(oname, out_cube, in_cube, exclude=None):
    """
    Generate single channelised images rather than a stacked cube
    """
    snitch.info("Unstacking / splitting the output cube")
    clean_hdr = in_cube.get_clean_header()

    for chan in range(out_cube.shape[0]):
        if (exclude is not None) and (chan in exclude):
            continue

        clean_hdr.update(
            {
                "BMAJ": in_cube.header.get(f"BMAJ{chan+1}"),
                "BMIN": in_cube.header.get(f"BMIN{chan+1}"),
                "BPA": in_cube.header.get(f"BPA{chan+1}"),
            }
        )

        _oname = oname.with_stem(oname.stem + f"-{chan:04}.fits")
        snitch.info(f"Writing:    {_oname}")

        save_fits(_oname, out_cube[chan], header=clean_hdr)


def fit_the_data(freq, data, deg=2):
    nf, nx, ny = data.shape
    data = data.reshape(nf, -1)
    snitch.info(f"Fitting data shape with {( nf, nx, ny)}")
    if hasattr(data, "mask"):
        snitch.info("Using masked fitting")
        coeffs = np.ma.polyfit(freq, data, deg=deg)
    else:
        coeffs = np.polyfit(freq, data, deg=deg)
    vander = np.polynomial.polynomial.polyvander(freq, deg)  # (F, degree+1)
    fitted = vander @ coeffs  # (F, Nx*Ny)

    residuals = data - fitted

    return (
        fitted.reshape(nf, nx, ny),
        coeffs.reshape(coeffs.shape[0], nx, ny),
        residuals.reshape(nf, nx, ny),
    )


def main(opts):
    global OUTDIR

    cube = opts.cube
    mask = opts.mask
    freq_fname = opts.freqs
    deg = opts.deg
    OUTDIR = opts.outdir

    if not freq_fname.exists():
        snitch.error(f"Frequency file {freq_fname.absolute()} was not found.")
        sys.exit(1)
    if not cube.exists():
        snitch.error(f"Cube file {cube.absolute()} was not found.")
        sys.exit(1)
    if mask and not mask.exists():
        snitch.error(f"Mask file {mask.absolute()} was not found.")
        sys.exit(1)

    if opts.exclude:
        if ":" in opts.exclude or "," in opts.exclude:
            excludes = selector(opts.exclude)
        else:
            snitch.info("Reading frequency selection as file")
            opts.exclude = Path(opts.exclude)
            if opts.exclude.exists():
                excludes = np.loadtxt(opts.exclude).astype(int)
            else:
                snitch.error(
                    f"The specified Exclusion file: {opts.exclude.absolute()} was not found."
                )
                sys.exit(1)
    else:
        excludes = None

    snitch.info(f"Creating/Overwriting directory {OUTDIR.absolute()}")
    OUTDIR.mkdir(parents=True, exist_ok=True)

    image = AstroImage(cube)
    image.frequency = np.loadtxt(freq_fname) / 1e9

    if excludes is not None:
        image.data[excludes, ...] = np.nan
        image.frequency[excludes] = np.nan

    if mask:
        mask_im = AstroImage(mask)
        mask_data = np.broadcast_to(~mask_im.data.astype(bool), image.data.shape)
        image.data = np.ma.array(data=image.data, mask=mask_data)

    model_data, coeffs, residuals = fit_the_data(
        image.frequency, image.data, deg=deg
    )

    model_hdr = deepcopy(image.header)

    model_hdr.remove("HISTORY", remove_all=True)
    clean_hdr = image.get_clean_header()

    if opts.unstack:
        split_cube(
            OUTDIR.joinpath(f"sc-ifit-i-model-deg-{deg}.fits"),
            model_data,
            image,
            exclude=excludes,
        )
    else:
        # start saving files
        save_fits(
            OUTDIR.joinpath(f"sc-ifit-i-model-deg-{deg}.fits"),
            model_data,
            header=model_hdr,
            overwrite=True,
        )
    save_fits(
        OUTDIR.joinpath(f"sc-ifit-i-residual-deg-{deg}.fits"),
        residuals,
        header=clean_hdr,
        overwrite=True,
    )

    for _ in range(len(coeffs)):
        save_fits(
            OUTDIR.joinpath(f"sc-ifit-coefficient-{_}.fits"),
            coeffs[_],
            header=clean_hdr,
            overwrite=True,
        )

    snitch.info("---- Done. Bye :D ! ----")


def console():
    """A console run entry point for setup.cfg"""
    new_fit_app()
    snitch.info("Bye :D !")


if __name__ == "__main__":
    console()
