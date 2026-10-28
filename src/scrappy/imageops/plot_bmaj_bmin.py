import os
from glob import glob
from types import SimpleNamespace

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from natsort import natsorted

matplotlib.use("Agg")

from operator import itemgetter

from scrappy.logging import logging, setup_streamhandler
from scrappy.utils.astroimage import AstroImage
from scrappy.utils.genutils import (
    append_extension,
    fullpath,
    make_out_dir,
    np_savetxt,
    np_savez,
    savefig,
    super_write,
)

snitch = logging.getLogger(__name__)
snitch.addHandler(setup_streamhandler())
snitch.setLevel("INFO")

OUTDIR = None


def get_input_files(indir: str, search_str: str) -> tuple[str]:
    images = natsorted(glob(os.path.join(indir, search_str)))
    snitch.info(f"Found {len(images)} FITS files meeting search criteria")
    return images


def filter_images_and_get_beam_dims(
    indir,
    search_str,
    threshold,
    selection: str = None,
    exclusion: str = None,
    autoselect: bool = False,
):
    images = get_input_files(indir, search_str)
    dims = get_image_beam_dims(indir, search_str)
    excludes = ()
    includes = list(range(len(images)))

    if selection is not None:
        includes = [int(_) for _ in selection.split(",")]
        snitch.info(f"You've chosen to include channels: {includes}")
    if exclusion is not None:
        excludes = [int(_) for _ in exclusion.split(",")]
        snitch.info(f"You've chosen to exclude channels: {excludes}")

    if (selection and exclusion) is None and autoselect:
        # then auto selection happens
        snitch.info(f"Autoselecting based on wsums and threshold {threshold}")
        wsums = np.round(dims.wsum / dims.wsum.max(), 2)
        excludes = np.where(wsums < threshold)[0].data.tolist()
        includes = np.where(wsums > threshold)[0].data.tolist()

    if len(excludes) > 0:
        _ex = list(itemgetter(*excludes)(images))
        snitch.info("Excluding the following images:")
        for _ in _ex:
            snitch.info(f"--> {_}")
        with super_write(fullpath(OUTDIR, "excluded-images.txt"), "w") as file:
            file.writelines([f"{_}\n" for _ in _ex])

    snitch.info("Including the following images:")
    images = list(itemgetter(*includes)(images))
    for _ in images:
        snitch.info(f"--> {_}")
        with super_write(fullpath(OUTDIR, "included-images.txt"), "w") as file:
            file.writelines([f"{_}\n" for _ in _ex])

    dims.norm_wsums = wsums
    dims.mask = np.zeros(dims.wsum.size, dtype=bool)
    dims.mask[excludes] = True
    return dims


def get_image_beam_dims(indir: str, search_str: str) -> SimpleNamespace:
    images = get_input_files(indir, search_str)

    # getting beam infomation
    dims = {"chan": [], "bmaj": [], "bmin": [], "bpa": [], "wsum": [], "freqs": []}
    for image in images:
        image = AstroImage(image)
        im_dims = image.get_beam_dims()
        dims["chan"].append(im_dims["chan"])
        dims["bmaj"].append(im_dims["bmaj"])
        dims["bmin"].append(im_dims["bmin"])
        dims["bpa"].append(im_dims["bpa"])
        dims["wsum"].append(image.wsum)
        dims["freqs"].append(image.frequency)

    # convert all to numpy array
    for key, value in dims.items():
        dims[key] = np.array(value)
    dims = SimpleNamespace(**dims)
    dims.mask = np.zeros(dims.wsum.size, dtype=bool)
    return dims


def plot_beams(data: SimpleNamespace, oname: str) -> None:
    oname = append_extension(oname, "png")
    snitch.info(f"Saving file at: {oname}")

    plt.close("all")
    fig, ax = plt.subplots(figsize=(16, 9))

    ax.plot(
        data.chan[~data.mask],
        data.bmaj[~data.mask],
        "o",
        markersize=8,
        color="tab:blue",
        alpha=0.7,
    )
    ax.set_xlabel("Channel numbers")
    ax.set_ylabel("BMAJ [Deg]")
    ax.minorticks_on()
    ax.grid(True, axis="x")
    ax.tick_params(axis="y", labelcolor="tab:blue")

    ax2 = ax.twinx()
    ax2.plot(
        data.chan[~data.mask],
        data.bmin[~data.mask],
        "o",
        markersize=8,
        color="tab:orange",
        alpha=0.7,
    )
    ax2.minorticks_on()

    ax2.set_ylabel("BMIN [Deg]")
    ax2.tick_params(axis="y", labelcolor="tab:orange")

    if data.mask.sum() > 0:
        ax.plot(
            data.chan[data.mask],
            data.bmaj[data.mask],
            "kx",
            alpha=0.8,
            markersize=8,
            label="Excluded",
        )
        ax2.plot(
            data.chan[data.mask],
            data.bmin[data.mask],
            "kx",
            alpha=0.8,
            markersize=8,
            label="Excluded",
        )

    if hasattr(data, "freqs"):
        ax3 = ax.twiny()
        ax3.plot(data.freqs[~data.mask] / 1e9, data.bmin[~data.mask], "o", alpha=0)
        ax3.set_xlabel("Freqs [GHz]")

    fig.tight_layout()
    savefig(fig, oname)
    return


def main(opts):
    global OUTDIR
    OUTDIR = "beam-plot-io" if opts.output is None else opts.output

    make_out_dir(OUTDIR)

    select = True if (opts.sel or opts.exclude or opts.auto_select) else False

    if select:
        beam_dims = filter_images_and_get_beam_dims(
            opts.idir,
            opts.search_str,
            opts.threshold,
            selection=opts.sel,
            exclusion=opts.exclude,
            autoselect=opts.auto_select,
        )
    else:
        beam_dims = get_image_beam_dims(opts.idir, opts.search_str)

    plot_beams(beam_dims, fullpath(OUTDIR, "beam-plots.png"))

    # saving some auxiliary files
    snitch.info("saving wsums of selected channels to text file")
    np_savetxt(fullpath(OUTDIR, "wsums.txt"), beam_dims.wsum[~beam_dims.mask])

    snitch.info("Saving selected frequencies")
    np_savetxt(fullpath(OUTDIR, "frequencies.txt"), beam_dims.freqs[~beam_dims.mask])

    snitch.info("Saving beam info for selected channels")
    np_savez(
        fullpath(OUTDIR, "beams"),
        bmajs=beam_dims.bmaj[~beam_dims.mask],
        bpas=beam_dims.bpa[~beam_dims.mask],
        bmins=beam_dims.bmin[~beam_dims.mask],
        freqs=beam_dims.freqs[~beam_dims.mask],
    )

    # largest beam size available
    with super_write(fullpath(OUTDIR, "beam-dims.txt"), "w") as file:
        file.write(f"{beam_dims.bmaj[~beam_dims.mask][0]}\n")
        file.write(f"{beam_dims.bmin[~beam_dims.mask][0]}\n")
        file.write(f"{beam_dims.bpa[~beam_dims.mask][0]}\n")

    return


def console():
    """A console run entry point for setup.cfg"""
    from scrappy.imageops.cli import plot_beam_app

    plot_beam_app()
    snitch.info("Bye :D !")


if __name__ == "__main__":
    console()
