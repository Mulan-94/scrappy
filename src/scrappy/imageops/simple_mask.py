import os
import matplotlib
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from itertools import zip_longest
from astropy.wcs import FITSFixedWarning
import warnings

from astropy.io import fits
from astropy.coordinates import SkyCoord, FK5
from astropy.wcs import WCS
from regions import (Regions, LineSkyRegion, RectangleSkyRegion)
from scipy.ndimage import rotate, label

from mpl_toolkits.axes_grid1 import make_axes_locatable

from scrappy.logging import logging, LOG_FORMATTER, setup_streamhandler

from scrappy.utils.genutils import timer, append_extension
matplotlib.use("Agg")


snitch = logging.getLogger(__name__)
snitch.addHandler(setup_streamhandler())
snitch.setLevel("INFO")


def get_wcs(name:str, pixels=False):
    """
    Return the image wcs object or the image dimensions

    name:
        Image from which to get the wcs stuff
    pixels: bool
        If specified, will return the image dimensions in pixels
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(fits.getheader(name))
    if wcs.naxis > 2:
        dropped = wcs.naxis - 2
        # remove extra and superflous axes. These become problematic
        for _ in range(dropped):
            try:
                # sometimes crazy things happen
                wcs = wcs.dropaxis(-1)
            except IndexError:
                continue
 
    if pixels:
        # get the image dimensions
        return wcs.pixel_shape
    else:
        return wcs


def world_to_pixel_coords(ra, dec, wcs_ref):
    """
    Convert world coordinates to pixel coordinates.
    The assumed reference is FK5
    ra: float
        Right ascension in degrees
    dec: float
        Declination in degrees
    wcs_ref:
        Image to use for WCS information

    Returns
    -------
        x and y pixel information
    """
    if isinstance(wcs_ref, str):
        wcs = get_wcs(wcs_ref)
    else:
        wcs = wcs_ref
    world_coord = FK5(ra=ra, dec=dec)
    skies = SkyCoord(world_coord)
    x, y = skies.to_pixel(wcs)
    return int(np.round(x,0)), int(np.round(y,0))

#######################################################


def merotate2(img, angle, pivot):
    """
    From:
    https://stackoverflow.com/questions/25458442/rotate-a-2d-image-around-specified-origin-in-python
    """
    padx = [img.shape[1] - pivot[0], pivot[0]]
    pady = [img.shape[0] - pivot[1], pivot[1]]
    paded = np.pad(img, [pady, padx], 'constant')
    rotated = rotate(paded, angle, reshape=False, output="uint8")
    return rotated[pady[0]:-pady[1], padx[0]:-padx[1]]


def gen_reg_mask(reg, wcs, shape):
    """
    reg: Region
        Region being worked on
    wcs:
        Wcs info for transformation
    shape: tuple
        Shape of the resulting image
    """
    reg = reg.to_pixel(wcs)
    mask = reg.to_mask()
    mask = mask.to_image(shape)
    return mask.astype("uint8")



def get_region_mask(wcs: astropy.wcs.wcs.WCS, reg: regions.Region, fill_val: int=1, default_arr: str="zeros"):
    """
    fname: str
        Name of the input image. To be used for wcs reference
    data: np.ndarray
        Numpy array containining image data
    reg: :obj:`Regions`
        Region being currently processed
    fill_value: int
        What to the chosen valid region area with. Default is 1. Meaning that
        within the bounds of this region, the value there will be 1.
    default_arr: str
        What the buffering array will be composed of. Intention is to be:
        'zeros' of 'ones' and nothing else! 
    """


    x_npix, y_npix = wcs.pixel_shape
    bitmask = getattr(np, default_arr)((x_npix,y_npix), dtype="uint8")

    if isinstance(reg, LineSkyRegion):
        snitch.info("""We don't do lines here :), skipping""")
        return bitmask

    if isinstance(reg, RectangleSkyRegion):
        # using a custom masker for rectangular regions because it gives
        # better results!!! Do NOT change!
        cx, cy = world_to_pixel_coords(reg.center.ra, reg.center.dec, wcs_ref=wcs)

        reg = reg.to_pixel(wcs)
        
        w, h = np.ceil(np.array((int(reg.width), int(reg.height)))//2).astype(int)
                
        minx, maxx = cx-w, cx+w+1
        miny, maxy = cy-h, cy+h+1

        xs = np.ma.masked_outside(np.arange(minx, maxx), 0, x_npix).compressed()
        ys = np.ma.masked_outside(np.arange(miny, maxy), 0, y_npix).compressed()

        # notice how we pass y before x, returns x before y
        # Leave it this way!!!!
        mesh = tuple(np.meshgrid(ys, xs))
        bitmask[mesh] = fill_val
        
        if hasattr(reg, "angle"):
            pivot = cx, cy
            bitmask = merotate2(bitmask, -reg.angle.value, pivot=pivot)
    else:
        bitmask = gen_reg_mask(reg, wcs, (x_npix, y_npix))

    return bitmask



def read_regions(regs: List[Path]):
    oregs = Regions()
    for reg in regs:
        oregs.extend(Regions.read(reg, format="ds9").regions)
    snitch.info(f"Found {len(oregs)} regions")
    return oregs


def plot_mask(data, oname):
    snitch.info(f"Writing output mask into {oname}")

    vmax = np.percentile(data, 97.5)
    vmin = np.percentile(data, 2.5)
    _mask, features = label(data)
    cmap = plt.get_cmap('viridis', features)

    fig, ax = plt.subplots(figsize=(16, 9))
    cs = ax.imshow(_mask, origin="lower", cmap=cmap, vmin=0, vmax=features)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(cs, cax=cax, label="Tag number")
    fig.tight_layout()
    fig.savefig(oname, dpi=400)
    return


def make_mask(fname, outname, above=None, below=None, regnames=None, ex_regnames=None, nodata=False):
    """
    Make simple mask

    fname: str
        Input image for reference
    outname: str
        Output image name for refrence
    above: float
        Values above which to mask
    below: float
        Values below which to mask

    Returns
    -------
    Mask: np.ndarray
        A numpy array containing the mask
    """
    data = fits.getdata(fname).squeeze()
    
    wcs = get_wcs(fname)
    hdr = wcs.to_header()
    if "HISTORY" in hdr:
        del hdr["HISTORY"]



    # cater for empty masks
    if nodata:
        data = np.ma.masked_array(data=np.zeros_like(data))
        data = np.ma.masked_array(data=data, mask=np.ma.make_mask(data))
    elif (above is None and below is None):
        # convert all masks to true, create an empty mask. ie everything is masked
        data = np.ma.masked_array(data=data, mask=np.ma.make_mask(data))
    else:
        if above is not None:
            # mask everything BELOW this value. Everything ABOVE it isn't masked
            data = np.ma.masked_less_equal(data, above)
        if below is not None:
            # mask everything ABOVE this value. Everything BELOW it isn't masked
            data = np.ma.masked_greater_equal(data, below)

    # here unwanted regions are set to one. But in our mask, we want the wanted
    # regions set to 1, so we invert the mask
    data = np.ma.masked_invalid(data)
    data_mask = ~data.mask
    data_mask = data_mask.astype("uint8")


    if regnames is not None:
        reg_mask = np.zeros(data.shape, dtype="uint8")
        regs = read_regions(regnames)

        for reg in regs:
            reg_mask = np.bitwise_or(reg_mask, get_region_mask(wcs, reg))

        if reg_mask.sum() == 0:
            snitch.info("Invalid region(s). We're ignoring this")
        else:
            data_mask = reg_mask * data_mask
        data_mask = reg_mask * data_mask
  
    if ex_regnames is not None:
        reg_mask = np.zeros(data.shape, dtype="uint8")
        regs = read_regions(ex_regnames)

        for reg in regs:
            reg_mask = np.bitwise_or(reg_mask, get_region_mask(wcs, reg))

        # invert this mask to contain 1s where
        reg_mask = np.logical_not(reg_mask).astype("uint8")

        if reg_mask.sum() == 0:
            snitch.info("Invalid region(s). We're ignoring this")
        else:
            data_mask = reg_mask * data_mask

    plot_mask(data_mask, f"{outname}-overlay-final.png")

    outname = append_extension(outname, "fits")

    snitch.info(f"Writing output mask into {outname}")

    if data_mask.max() == 1:
        data_mask = data_mask.astype("uint8")
    fits.writeto(outname, data_mask, header=hdr, overwrite=True)

    return


def main(opts):

    EMPTY = ""
    
    if opts.oname is None:
        opts.oname = (EMPTY)
    if opts.regname is None:
        opts.regname = (EMPTY)
    
    if opts.oname is None:
        oname = f"output-mask-{i}.fits"
    
    make_mask(opts.iname, opts.oname, above=opts.above, below=opts.below,
                regnames=opts.regname, ex_regnames=opts.ex_regname, nodata=opts.nodata)

    # for i,( oname, regname) in enumerate(zip_longest(opts.oname, opts.regname)):
    #     if oname is None:
    #         oname = f"output-mask-{i}.fits"
    #     make_mask(opts.iname, oname, above=opts.above, below=opts.below,
    #         regname=regname, ex_regname=opts.ex_regname, nodata=opts.nodata)
    snitch.info('------- Done ------')


def console():
    """A console run entry point for setup.cfg"""
    from scrappy.imageops.cli import sm_app
    sm_app()
    snitch.info("Bye :D !")
    

if __name__ == "__main__":
    console()

