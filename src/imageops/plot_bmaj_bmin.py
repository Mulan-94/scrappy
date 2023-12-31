#############
## Lexy plotting BMAJ bmin BLABAL
#######################


import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

from astropy.io import fits
from glob import glob
from natsort import natsorted


from utils.logger import logging, LOG_FORMATTER, setup_streamhandler
from utils.genutils import make_out_dir
snitch = logging.getLogger(__name__)
snitch.addHandler(setup_streamhandler())
snitch.setLevel("INFO")

def read_fits(imname):
    outs = {}
    hdr = fits.getheader(imname)
    for _ in "BMAJ BMIN BPA CRVAL3".split():
            outs[_] = hdr[_]
    snitch.info(f"Reading         :{imname}")
    return outs


def read_cube_fits(imname, channels=None):
    outs = []

    with fits.open(imname) as hdul:
        if channels is not None:
            snitch.info(f"Selecting {len(ps.channels)} channels out of {hdul[0].header['NAXIS3']}")
            snitch.info(f"Channels: {channels}")
            channels = np.array(channels)
            channels += 1
        else:
            channels = range(1, hdul[0].header["NAXIS3"]+1)
        for freq_id in channels:
            outs.append((
                hdul[0].header[f"BMAJ{freq_id}"],
                hdul[0].header[f"BMIN{freq_id}"],
                hdul[0].header[f"BPA{freq_id}"],
                freq_id))
    return outs


def plotme(freqs, bmajs, bmins, outname, title="All freqs"):
    chans = np.arange(freqs.size)

    snitch.info("Plotting beam dimensions")
    
    fig, ax = plt.subplots(figsize=(16,9), ncols=1)
    
    ax.plot(chans, bmajs, "o", markersize=8, color="tab:blue", alpha=0.7)
    ax.set_xlabel("Channel numbers")
    ax.set_ylabel("BMAJ [Deg]")
    ax.minorticks_on()
    ax.grid(True, axis="x")
    ax.tick_params(axis='y', labelcolor="tab:blue")



    ax2 = ax.twinx()
    ax2.plot(chans, bmins, "o", markersize=8, color="tab:orange", alpha=0.7)
    ax2.minorticks_on()

    ax2.set_ylabel("BMIN [Deg]")
    ax2.tick_params(axis='y', labelcolor="tab:orange")

    fig.tight_layout()
    snitch.info(f"Saving file at: {outname}")
    fig.savefig(outname)

    return


def get_params(pairs):
    snitch.info("Populating beam parameters for all images")
    bmajs = np.array([_["BMAJ"] for _ in pairs])
    bmins = np.array([_["BMIN"] for _ in pairs])
    bpas = np.array([_["BPA"] for _ in pairs])
    freqs = np.array([_["CRVAL3"] for _ in pairs])
    return bmajs, bmins, bpas, freqs



def get_and_plot_beam_info(indir=None, search="*[0-9]-I-image.fits", dump=".",
    prefix=None, oname="beam_vs_freq.png"):
    """
    read file containing filenames for the images to be processed
    these are stored in e.g eds.txt
    if a file containing these file names is not avaliable, 
    just use the files as they are
    """

    data = natsorted(glob(os.path.join(indir, search)))

    pairs = [read_fits(dat) for dat in data]
    bmajs, bmins, bpas, freqs = get_params(pairs)


    # save this beam information into beams file
    np.savez(os.path.join(dump, "beams" if prefix is None else f"{prefix}-beams"),
        bpas=bpas, bmajs=bmajs, bmins=bmins,
        freqs=freqs)

    oname = os.path.join(dump, oname if prefix is None else f"{prefix}-{oname}")
    plotme(freqs, bmajs, bmins, oname, title="All freqs")

    return


# def single_cube_file(cube_name, oname=None, channels=None):
#     ## in the case of multiple data cubes
#     pairs = read_cube_fits(cube_name, channels=channels)
#     bmajs, bmins, bpas, freqs = get_params(pairs)
#     if oname is None:
#         oname = f"bmaj-bmin-vs-freq-{cube_name}.png"
#     plotme(freqs, bmajs, bmins, oname, title="All freqs")



#------------------------------------------============
# Some Basic automated channel selection using wsums
#------------------------------------------============
def read_wsums(image):
    hdr = fits.getheader(image)
    return hdr["WSCVWSUM"]


def channel_selection(folder, dump, threshold=0.5, autoselect=True, sel=None,
    exclude=None):
    """
    folder: str
        The directory containing your intended images
    threshold: float
        A value between zero and one. We look at the individual images wsum
        which is normalised by the maximum available. This threshold is checked
        against the normalised wsum. Values below the threshold are ignored.
    """
    snitch.info("Starting channel selection")

    images = natsorted(glob(os.path.join(folder, "*[0-9]-I-image.fits")))
    owsums = np.zeros(len(images))
    for i, image in enumerate(images):
        owsums[i] = read_wsums(image)

    ONAME = os.path.join(dump, "orig-wsums.txt")
    snitch.info(f"Saving WSUMS to: {ONAME}")
    np.savetxt(ONAME, owsums)

    if exclude:
        exclude = exclude.split(",")
        if len(exclude) > 1 or exclude[0].isnumeric():
            exc = np.array(exclude).astype(int)
        else:
            exclude = exclude[0]
            exc = np.loadtxt(exclude).astype(int)
        snitch.info(f"You've chosen to exclude channels: {', '.join(exc.astype(str))}")
        owsums[exc] = 0

    
    if autoselect:
        wsums = np.round(owsums/owsums.max(), 2)
        not_sel, = np.where(np.ma.masked_less_equal(wsums, threshold).mask==True)
        ONAME = os.path.join(dump, "not-selected-channels.txt")
        snitch.info(f"Saving channel not selected to: {ONAME}")
        with open(ONAME, "w") as file:
            file.writelines([f"{_}".zfill(4) + "\n" for _ in not_sel])

        sel, = np.where(np.ma.masked_less_equal(wsums, threshold).mask==False)
        snitch.info(f"{sel.size} of {wsums.size} channels selected.")
        ONAME = os.path.join(dump, "selected-channels.txt")
        snitch.info(f"Saving channel sel to: {ONAME}")
        with open(ONAME, "w") as file:
            file.writelines([f"{_}".zfill(4) + "\n" for _ in sel])

        ONAME = os.path.join(dump, "wsums.txt")
        snitch.info(f"Saving WSUMS of selected channels to: {ONAME}")
        np.savetxt(ONAME, owsums[sel])
       

    else:
        sel = sel.split(",")
        if len(sel) > 1 or sel[0].isnumeric():
            sel = np.array(sel).astype(int)
        else:
            sel = sel[0]
            snitch.info(f"Reading custom selected channel file: {sel}")
            sel = np.loadtxt(sel).astype(int)
        snitch.info(f"You've selected channels: {', '.join(sel.astype(str))}")
    
        not_sel = np.zeros(len(images))
        not_sel[sel] = 1
        not_sel, = np.where(not_sel==0)

        ONAME = os.path.join(dump, "wsums.txt")
        snitch.info(f"Saving WSUMS of selected channels to: {ONAME}")
        np.savetxt(ONAME, owsums[sel])

        ONAME = os.path.join(dump, "selected-channels.txt")
        snitch.info(f"Saving channel sel to: {ONAME}")
        with open(ONAME, "w") as file:
            file.writelines([f"{_}".zfill(4) + "\n" for _ in sel])
       

        ONAME = os.path.join(dump, "not-selected-channels.txt")
        snitch.info(f"Saving channel not selected to: {ONAME}")
        with open(ONAME, "w") as file:
            file.writelines([f"{_}".zfill(4) + "\n" for _ in not_sel])

    
    return sel, not_sel


def read_and_plot_beams2(folder, dump=".", prefix=None, beam_file="beams.npz",
    threshold=0.5, sel=None, autoselect=True, exclude=None):
    """
    folder: str
        Where the input images are
    dump: str
        Where to dump the outputs
    """
    
    sel, not_sel = channel_selection(folder, dump, threshold=threshold,
        sel=sel, autoselect=autoselect, exclude=exclude)
    
    chans = [f"{_}".zfill(4) for _ in sel]
    
    bm = dict(np.load(os.path.join(dump, beam_file)))
    bma, bmi, bpa, freqs = bm["bmajs"], bm["bmins"], bm["bpas"], bm["freqs"]

    # write the largest beam available
    ONAME = os.path.join(dump, "beam-dims.txt")
    with open(ONAME, "w") as file:
        file.write(f"{bma[sel[0]]}\n")
        file.write(f"{bmi[sel[0]]}\n")
        file.write(f"{bpa[sel[0]]}\n")


    
    ONAME = os.path.join(dump, "frequencies.txt")
    snitch.info(f"Saving selected freqs : {ONAME}")

    np.savetxt(ONAME, freqs[sel])


    chans = np.arange(bma.size)

    plt.close("all")
    fig, ax = plt.subplots(figsize=(16, 9), ncols=1, nrows=1)

    ax.plot(chans[sel], bma[sel], "o", alpha=0.7, markersize=8, color="tab:blue")
    ax.plot(chans[not_sel], bma[not_sel], "kx", alpha=0.8, markersize=8, label="Excluded")
    ax.minorticks_on()
    ax.grid(True, axis="x")
    ax.set_ylabel("BMAJ [Deg]", color="tab:blue")
    ax.set_xlabel("Channel Number")
    ax.tick_params(axis='y', labelcolor="tab:blue")

    ax2 = ax.twinx()


    ax2.plot(chans[sel], bmi[sel], "o", alpha=0.7, markersize=8, color="tab:orange")
    ax2.plot(chans[not_sel], bmi[not_sel], "kx", alpha=0.8, markersize=8, label="Excluded")
    ax2.minorticks_on()
    ax2.set_ylabel("BMIN [Deg]")
    ax2.set_xlabel("Channel Number")
    ax2.tick_params(axis='y', labelcolor="tab:orange")
    ax.legend()

    ax3 = ax.twiny()
    ax3.plot(freqs[sel]/1e9, bmi[sel], "o", alpha=0)
    ax3.set_xlabel("Freqs [GHz]")

    
    if prefix is None:
        ONAME = os.path.join(dump, "selected_beam_vs_freq.png")
    else:
        ONAME = os.path.join(dump, f"{prefix}-selected_beam_vs_freq.png")
    fig.savefig(ONAME, bbox_inches="tight", dpi=300)

    return
    


#------------------------------------------


def parser():
    ps = argparse.ArgumentParser()

    
    ps.add_argument("idir", metavar="", type=str,
        help="Directory where to find the input images"
    )

    ps.add_argument("-o", "--output", dest="output", metavar="",
        type=str, default=".", 
        help="Directory where to dump outputs. Default is current directory."
    )
    ps.add_argument("-search", "--search", type=str,
        help="Regex of what string to search in the directory",
        default="*[0-9q]-I-image.fits")

    ps.add_argument("-t", "--threshold", dest="threshold", metavar="",
        type=float, default=0.5, 
        help="""Channels with WSUM/WSUM.max() below this value will be excluded.
        Only used if -s is active"""
    )
    ps.add_argument("-sc", "--select-channels", dest="sel", default=None,
        type=str, help="""Custom channel selections. This can be a file 
        containing those channels, or a comma separated list
        e.g. '0,1,2,3'"""
        )
    ps.add_argument("-as", "--auto-select", action="store_true", dest="auto_select", 
        help="Try and suggest valid channels for selection. Use in conjuction with '-t'.")

    ps.add_argument("-pre", "--prefix", type=str, default="00", 
        help="Prefix to append on the output images")

    ps.add_argument("-ec", "--exclude-channels", dest="exclude", type=str,
        help="""Channels to be excluded during autoselection. 
        This can be a file containing those channels, or a comma separated list
        e.g. '0,1,2,3'.""")
    return ps


def main():
    ps = parser().parse_args()

    if ps.output:
        make_out_dir(ps.output)

    # get info about beams and plot them
    get_and_plot_beam_info(ps.idir, search=ps.search, dump=ps.output, prefix=ps.prefix)
       
    if ps.auto_select or ps.sel:
        read_and_plot_beams2(ps.idir, dump=ps.output, beam_file=f"{ps.prefix}-beams.npz",
            threshold=ps.threshold, sel=ps.sel, autoselect=ps.auto_select,
            exclude=ps.exclude)


def console():
    """A console run entry point for setup.cfg"""
    main()
    snitch.info("Bye :D !")
    

if __name__ == "__main__":
    console()
