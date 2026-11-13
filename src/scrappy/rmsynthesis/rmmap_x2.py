#! /usr/bin/env python

"""
Courtesy of Lerato
Edited by Lexy
"""

from itertools import product
from pathlib import Path

import dask.array as da
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from pydantic import validate_call
from scipy import signal

from scrappy.logging.config import configure_logger
from scrappy.utils.astroimage import AstroImage
from scrappy.utils.genutils import save_fits
from scrappy.utils.rmmath import (
    SPEED_OF_LIGHT,
    RMSynthesis,
    frac_polzn,
)

matplotlib.use("Agg")

snitch = configure_logger(level="DEBUG")


OUTDIR = None

GET_FRAC_POL = False


@validate_call
def make_new_2d_fits(fname: Path, hdr):
    nx, ny = hdr["NAXIS2"], hdr["NAXIS1"]
    data = np.empty((nx, ny))
    save_fits(fname, data=data, header=hdr, overwrite=True)
    return fname


@validate_call
def debug_diagnostic_plots(*datas: tuple, xaxis=None, fname: Path = "comparison"):
    fname = fname.with_suffix(".png")
    plt.close("all")

    for i, data in enumerate(datas):
        plt.plot(xaxis, data, label=f"p-{i}")
    print(f"Plotting --> {fname}")
    plt.legend()
    plt.savefig(fname)


# def read_freq_file(filename):
#     _, ext = os.path.splitext(filename)

#     reader = {
#         ".txt": np.loadtxt,
#         ".dat": np.loadtxt,
#         ".npy": read_npz,
#         ".npz": read_npz,
#     }

#     try:
#         freqs = reader[ext.lower()](filename)

#         if isinstance(freqs, dict):
#             snitch.info("The following keys were found in frequency file")
#             snitch.info("{}".format(bold_blue(list(freqs.keys()).join(", "))))
#             snitch.info(f"Using {freqs.get(freqs.keys()[0])}")
#             freqs = freqs.get(freqs.keys()[0])
#         elif isinstance(freqs, tuple):
#             freqs = np.array(freqs)

#     except ValueError:
#         snitch.error(f"Unable to reaa {red_bold(filename)}. Please check on it")
#         sys.exit()
#     return freqs


def apply_mask(filename, data):
    mask = fits.getdata(filename)
    mask.dtype = bool
    mask = ~mask
    mask = np.broadcast_to(mask, shape=data.shape)

    if isinstance(data, np.ndarray):
        data = np.ma.MaskedArray(data=data, mask=mask, fill_value=0)
    elif isinstance(data, da.Array):
        data = da.ma.masked_array(data=data, mask=mask, fill_value=0)
        # we assume that its a dask array

    # return data.filled()? # a little slow but can be explored if necessary
    return data


# def chunked_cube_read(fname, xchunk=1024, ychunk=1024):
#     def get_image_slices(xiters, yiters):
#         xs, ys = 0, 0
#         sections = []
#         for xiter in range(xiters):
#             for yiter in range(yiters):
#                 xstop = xs + xchunk
#                 ystop = ys + ychunk
#                 # :, slice(xs, xstop), slice(yx, ystop)
#                 snitch.debug(
#                     f"This is slice {xiter, yiter}: as ':, slice{xs, xs+xchunk}, slice{ys, ys+xchunk}'"
#                 )
#                 sections.append(
#                     (slice(None), slice(xs, xs + xchunk), slice(ys, ys + ychunk))
#                 )
#                 ys += ychunk
#             xs += xchunk
#             ys = 0
#         return sections

#     GB = 2**30
#     MB = 2**20

#     # data structure is freq x xaxis x yaxis
#     # for RM Synthesis , freq axis MUST BE FULL, DO NOT CHUNK!!
#     xchunk, ychunk = 1024, 1024

#     hdr = fits.getheader(fname)
#     nx, ny, nfreq = tuple(map(hdr.get, "NAXIS1 NAXIS2 NAXIS3".split()))

#     xiters = ceil(nx / xchunk)
#     yiters = ceil(ny / ychunk)

#     sections = get_image_slices(xiters, yiters)

#     with fits.open(fname, mmap=True, mode="denywrite") as hdul:
#         hdu = hdul[0]

#         hdu.data = hdu.data.squeeze()

#         data_size = hdu.size / GB
#         snitch.info(f"This data is {bold_yellow(data_size)} GB in size")
#         snitch.info(f"The data is chunked into {len(sections)} sections")

#         for i, section in enumerate(sections):
#             snitch.debug(f"Yielding section {i}: {section}")
#             yield hdu.section[section]

#     return


# def main2(opts):
#     """
#     1. Validate that the input files exits
#     2. Read the frequency file
#     3. Read the actual data files
#     4. If mask is available, read that mask
#     5. Now mask the data (just one, the rest ops with MA will follow suite)
#     6. Concoct the required datas, fpol, lpol and whatever
#     -- Maybe write out the FITS data templates at this point?
#     7. Perform RM-synthesis on that data
#     8
#     """

#     validate_required_input_files(args)
#     frequencies = read_freq_file(opts.freq)
#     wavelength_sq = simple_lsq(frequencies)

#     # ----------------------------------------------------------------------
#     # testing chunked reading

#     tqdata = chunked_cube_read(opts.qcube)
#     tudata = chunked_cube_read(opts.ucube)
#     for sect, (qdata, udata) in enumerate(zip(tqdata, tudata), 1):
#         pdata = linear_polzn(qdata, udata)
#         sumn = rm_synthesis(
#             wavelength_sq,
#             pdata,
#             phi_max=opts.max_depth,
#             phi_step=opts.depth_step,
#             niter=opts.niters,
#             gain=0.1,
#             plot=False,
#             derotate=False,
#         )
#         print("-----------")

#     # ----------------------------------------------------------------------

#     # read the data files in
#     read = "dask"  # for testing
#     qdata = read_image_cube(opts.qcube, reader=read)
#     udata = read_image_cube(opts.ucube, reader=read)

#     if opts.mask:
#         # mask this one thing and every other operation will be masked
#         qdata = apply_mask(opts.mask, qdata)

#     if opts.icube:
#         idata = read_image_cube(opts.icube)
#         fp_data = frac_polzn(idata, qdata, udata)

#     pdata = linear_polzn(qdata, udata)

#     output_files = prep_output_files(opts.qcube, OUTDIR)

#     # RM-synthesis now
#     # Options
#     # 1. process the entire dataset
#     # 2. Manually chunk the data (maybe 1k x 1k x 4k)? and process iteratively
#     # 3. Use dask chunking

#     snr = opts.snr


# ---------------------------- End of New Functions ------------------------


def faraday_to_lambda(lam2, phi_range, pol_lam2):
    """
    Computes Faraday Spectra using RM-Synthesis
    as defined by Brentjens and de Bruyn (2005)

    Returns:
    fdata: shape(phi_range.size)

    """

    N = len(lam2)
    l20 = lam2.mean()
    fdata = np.zeros([len(phi_range)], dtype=complex)

    for k, phi in enumerate(phi_range):
        kernel = np.exp(-2j * (lam2 - l20) * phi)
        fdata[k] = pow(N, -1) * np.sum(pol_lam2 * kernel)
    return fdata


def rm_clean(lam2, phi_range, fspectrum, niter=500, gain=0.1, derotate=True):
    # Lexy, copy the value and not the reference
    fspectrum = fspectrum.copy()

    if derotate:
        fwhm = 3.8 / abs(lam2[0] - lam2[-1])
    else:
        # use the peak of the real component for the gaussian as the CLEAN beam
        # according to rudnick and cotton 2023
        fwhm = 2 / (lam2[-1] + lam2[0])

    sigma = fwhm / 2.35482
    Gauss = np.exp(-0.5 * (phi_range / sigma) ** 2)

    # I am padding here to avoid edge effects.

    pad = abs(phi_range[-1]) * 2
    dpad = abs(phi_range[0] - phi_range[1])
    phi_pad = np.arange(-pad, pad, dpad)
    dshift = int(pad / (2.0 * dpad))

    rmsf_fixed = faraday_to_lambda(lam2, phi_pad, 1)
    debug_diagnostic_plots(rmsf_fixed, xaxis=phi_pad, fname="orig-fixed-rmsf")
    components = np.zeros([len(phi_range)], dtype=complex)

    for _n in range(niter):
        temp = np.zeros([len(phi_range)], dtype=complex)
        f_amp = np.absolute(fspectrum)
        ind = np.where(f_amp == f_amp.max())[0]
        f_comp = fspectrum[ind[0]] * gain
        temp[ind[0]] = f_comp
        components += temp

        # # replacing
        dirac = np.zeros(len(phi_range))
        dirac[ind[0]] = 1
        rmsf = signal.convolve(rmsf_fixed, dirac, mode="same")
        rmsf_center = rmsf[dshift : -dshift + 1]

        fspectrum -= f_comp * rmsf_center

    Fres = fspectrum
    fclean = signal.convolve(components, Gauss, mode="same") + Fres

    return fclean, components, fwhm


def rm_synthesis(
    lam,
    pdata,
    phi_max=500,
    phi_step=5,
    niter=1000,
    gain=0.1,
    plot=False,
    derotate=True,
):
    phi_range = np.arange(-phi_max, phi_max + phi_step, phi_step)
    # this ensures that the middle value is zero.
    fdirty = faraday_to_lambda(lam, phi_range, pdata)

    fclean, fcomp, fwhm = rm_clean(
        lam, phi_range, fdirty, niter=niter, gain=gain, derotate=derotate
    )

    return phi_range, fclean, fdirty, fwhm


# -==========================================================================
# some RM noise stuff
# -==========================================================================
def rmtools_median_abs_dev(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:
    median(abs(a - median(a))) / c
    c = 0.6745 is the constant to convert from MAD to std
    """

    a = np.ma.masked_where(a != a, a)
    if a.ndim == 1:
        d = np.ma.median(a)
        m = np.ma.median(np.ma.fabs(a - d) / c)
    else:
        d = np.ma.median(a, axis=axis)
        if axis > 0:
            aswp = np.ma.swapaxes(a, 0, axis)
        else:
            aswp = a
        m = np.ma.median(np.ma.fabs(aswp - d) / c, axis=0)

    return m


def median_abs_dev(arr):
    """median Absolute Deviation"""
    mad = np.nanmedian(np.abs(arr - np.nanmean(arr)))
    return mad


def mean_abs_dev(arr):
    """Mean Absolute Deviation (MAD) statistic"""
    mad = np.nansum(np.abs(arr - np.nanmean(arr))) / arr.size
    return mad


def root_mean_square(arr):
    rms = np.sqrt(np.nanmean(np.square(arr)))
    return rms


# def peak_snr(fdf, rmsf_fwhm, phi_range):
#     """
#     With individual noise for Q and U because they are not the same

#     1. Get the complex FDF
#     2. Get its amplitudeqlq
#     3. get index of at peak amplitude
#     4. Get peak amplitude
#     5. mask the THE WIDTH of the FDF peak
#     6. select all the unmasked, and remaining values of the absolute FDF
#     7. if fraction of data remaining is less than 30%, ie if more 70pc of the data is flagged:
#             - MAD
#             - RMS of the unflagged fdf
#         otherwise, using the flagged data:
#             - MAD
#             - RMS to generate the noise

#         This will give the fdf noise
#     8. snr is therefore peak fdf to this noise
#     """
#     fdf_amp = np.abs(fdf)

#     max_peak = np.max(fdf_amp)
#     max_peak_loc = np.nanargmax(fdf_amp[1:-1])

#     # find the phi range where the main rmsf peak is located
#     phi_delta = np.nanmin(np.diff(phi_range))

#     rmsf_fwhm_channel = np.ceil(rmsf_fwhm / phi_delta)
#     rmsf_width_start = int(max(0, max_peak_loc - rmsf_fwhm_channel * 2))
#     rmsf_width_end = int(max(0, max_peak_loc + rmsf_fwhm_channel * 2))

#     flagged_fdf_amp = np.copy(fdf_amp)
#     flagged_fdf_amp[rmsf_width_start:rmsf_width_end] = np.nan
#     flagged_fdf_amp = np.ma.masked_invalid(flagged_fdf_amp).compressed()

#     if flagged_fdf_amp.size / fdf_amp.size < 0.3:
#         rms = root_mean_square(fdf_amp)
#         mad = median_abs_dev(fdf_amp)
#         # mad = rmtools_median_abs_dev(fdf_amp)
#     else:
#         rms = root_mean_square(flagged_fdf_amp)
#         mad = median_abs_dev(flagged_fdf_amp)
#         # mad = rmtools_median_abs_dev(flagged_fdf_amp)

#     # signal to noise here
#     snr = max_peak / mad

#     return snr


# -==========================================================================


# def call_fitting(
#     args, wavelengths=None, max_depth=500, niters=100, phi_step=1, derotate=True
# ):
#     """
#     Returns
#     -------
#     p0_map
#     PA_map
#     RM_map
#     fp_map

#     peak_rm_amp: p0_map
#         Amplitude of the clean peak RM
#     pol_angle: PA_map
#         Polarisation angle at clean peak
#     peak_depth
#         Faraday depth at clean RM peak
#     cleaned_amp
#         Amplitude of cleaned the entire Farady  spectra
#     peak_fpol
#         Fpol at the center freq, or the peak at all freqs


#     pdata: is the total linearly polarised intensity
#     fpdata: is the fractional polarisation
#     """
#     x, y = args

#     snitch.info(f"Processing pixel {x} x {y}")
#     # all freqs single pixel
#     single_pxl_poln = pdata[:, x, y]
#     ind_nans = ~np.isnan(abs(single_pxl_poln))
#     single_pxl_poln = single_pxl_poln[ind_nans]
#     wave_sq = wavelengths[ind_nans]

#     # RM-synth and RM-CLEAN on this single pixel. We pass in the complex polarised data
#     phi_range, fcleaned, fdirty, fwhm = rm_synthesis(
#         wave_sq,
#         single_pxl_poln,
#         phi_max=max_depth,
#         phi_step=phi_step,
#         niter=niters,
#         gain=0.1,
#         plot=False,
#         derotate=derotate,
#     )

#     # 6
#     snr_clean = peak_snr(fcleaned, fwhm, phi_range)
#     # snr_dirty = peak_snr(fdirty, fwhm, phi_range)
#     # return nans for pixel if snr is less than 10
#     # if snr_clean < 10:
#     #     return [np.nan]*6 + [snr_clean]

#     # 3 is fcleaned(complex)
#     abs_fclean = np.abs(fcleaned)

#     # get the peak index of the clean faraday spectra
#     peak_idx = np.where(abs_fclean == abs_fclean.max())[0][0]

#     # get the depth at which RM is peak. Is this the estimated RM?
#     # 2
#     peak_depth = phi_range[peak_idx]

#     # get the peak complex fclean
#     peak_complex_rm = fcleaned[peak_idx]

#     # this gives the estimated intrinsic fractional polarisation??
#     # 0
#     peak_rm_amp = np.abs(peak_complex_rm)

#     # get the polarisation angle at the peak fday spectrum (intrinsic fpol?)
#     # of that component. Note that this is from the cleaned FDF
#     # 1
#     pol_angle = 0.5 * np.arctan2(peak_complex_rm.imag, peak_complex_rm.real)

#     # Trying to get the frac pol at max polarised intensity for this x,y pixel
#     # nfp = np.abs(fp_data[:, x, y])
#     nfp = fp_data[:, x, y]

#     # get location of peak polzn intensity
#     # 5
#     peak_lpol = np.abs(single_pxl_poln).max()
#     # max_pol_idx = np.where(np.abs(single_pxl_poln)==peak_lpol)[0][0]
#     # snitch.info(f"x: {y}, y:{x}, max polzn intensity @ chan: {max_pol_idx}")

#     # take fpol at the highest frequency
#     # 4
#     peak_fpol = nfp[-1]

#     return (
#         peak_rm_amp,
#         pol_angle,
#         peak_depth,
#         abs(fcleaned),
#         peak_fpol,
#         peak_lpol,
#         snr_clean,
#     )


# def read_mask(fits):
#     """get pixel coordinates for a mask"""
#     maskdata, mhdr = read_data(fits, freq=False)
#     xpix, ypix = np.where(maskdata > 0)

#     return xpix, ypix


# def depolarisation_ratio(idata, qdata, udata, e_file, mask=None, output="depol-map"):
#     """
#     Calculate the depolarization ratio

#     Parameters
#     ----------
#     (i|q|u)data: np.ndarray
#         Data cubes containing the stokes
#     e_file: str
#         Name of the file containing errors. This is found in the file generated
#         from scrappy's 'sc-los' tool and is a numpy pickle file.
#     mask: str
#         Name of mask to use while calculating depolarisation ratio
#     output: str
#         Name of output depolarisation ratio map

#     Returns
#     -------
#     Maps of depolarization ratio and its error

#     """

#     def depolzn_err(fp_lo, fer_lo, fp_hi, fer_hi):
#         return np.sqrt(np.square(fer_lo / fp_lo) + np.square(fer_hi / fp_hi))

#     def read_npz(fname):
#         return dict(np.load(fname))

#     ihdr = fits.getheader(idata)
#     idata = fits.getdata(idata).squeeze()
#     qdata = fits.getdata(qdata).squeeze()
#     udata = fits.getdata(udata).squeeze()

#     sel = [0, -1]

#     npz = read_npz(e_file)
#     ierr, qerr, uerr = npz["I_err"][sel], npz["Q_err"][sel], npz["U_err"][sel]

#     idata, qdata, udata = idata[sel], qdata[sel], udata[sel]

#     if mask:
#         mask = fits.getdata(mask).squeeze()
#         idata = np.ma.masked_equal(idata * mask, 0).filled()
#         qdata = np.ma.masked_equal(qdata * mask, 0).filled()
#         udata = np.ma.masked_equal(udata * mask, 0).filled()

#     fpol = frac_polzn(idata, qdata, udata)

#     fperr_lo = frac_polzn_error(
#         idata[0], qdata[0], udata[0], ierr[0], qerr[0], uerr[0]
#     )
#     fperr_hi = frac_polzn_error(
#         idata[-1], qdata[-1], udata[-1], ierr[-1], qerr[-1], uerr[-1]
#     )

#     fp_lo, fp_hi = fpol[0], fpol[-1]

#     depoln = fp_lo / fp_hi
#     depoln_err = depolzn_err(fp_lo, fperr_lo, fp_hi, fperr_hi)
#     depoln_err_ratio = depoln_err / depoln

#     depol_hdr = modify_fits_header(ihdr, ctype="depol", unit="unit")

#     # depolarization
#     fits.writeto(output + ".fits", depoln, depol_hdr, overwrite=True)

#     # depolarization
#     fits.writeto(output + "-err.fits", depoln_err, depol_hdr, overwrite=True)

#     fits.writeto(
#         output + "-err-ratio.fits", depoln_err_ratio, depol_hdr, overwrite=True
#     )

#     # write out fpol error at the highest freq (shortest wavelength)
#     fperr_hdr = modify_fits_header(ihdr, ctype="fpol_err", unit="unit")
#     fits.writeto(output + "-fpol-err.fits", fperr_hi, fperr_hdr, overwrite=True)

#     snitch.info("Depolarization maps ready")
#     return


# def depol_shim():
#     snitch.info("Generating depolarization maps")
#     snitch.info(f"Running: depolarisation_ratio({', '.join(sys.argv[1:])})")
#     depolarisation_ratio(*sys.argv[1:])
#     return


def test_orig_rm(lam2, pdata, pixel: tuple):
    phi_range = np.arange(-200, 200 + 1)
    one_pix = pdata[:, *pixel]

    fspec_dirty = faraday_to_lambda(lam2, phi_range, one_pix)

    fspec_clean, components, fwhm = rm_clean(
        lam2, phi_range, fspec_dirty, niter=500, gain=0.1, derotate=True
    )
    return fspec_dirty, fspec_clean, components, fwhm


def main(opts):
    global GET_FRAC_POL, OUTDIR

    OUTDIR = opts.outdir

    OUTDIR.mkdir(parents=True, exist_ok=True)

    # if opts.snr:
    #     snr = opts.snr

    phi_range = np.arange(-opts.max_depth, opts.max_depth + 1, opts.depth_step)

    freqs = np.loadtxt(opts.freqs)
    waves = SPEED_OF_LIGHT / freqs
    wave_sq = waves**2

    qcube = AstroImage(opts.qcube)
    ucube = AstroImage(opts.ucube)
    qcube.frequency = freqs
    qcube.wavelength = waves
    ucube.frequency = freqs
    ucube.wavelength = waves

    if opts.icube:
        GET_FRAC_POL = True
        icube = AstroImage(opts.icube, chunked=False)
        icube.frequency = freqs
        icube.wavelength = waves

    if opts.mask:
        mask = AstroImage(opts.mask)
        mask.data = np.broadcast_to(~mask.data.astype(bool), qcube.data.shape)
        # qcube.data = np.ma.array(data=qcube.data, mask=mask.data, fill_value=0)

    clean_hdr = qcube.get_clean_2d_header()

    chunks = 5
    if GET_FRAC_POL:
        fpol_data = frac_polzn(icube.data, qcube.data, ucube.data)[
            :, :chunks, :chunks
        ]

    # lpol_data = linear_polzn(qcube.data, ucube.data)[:, :chunks, :chunks]

    # Debuggging
    # for testing
    pixels = list(product(range(chunks), repeat=2))
    dirties = []
    cleans = []

    for pixel in pixels:
        print(f"processing pixel for original algorithm: {pixel}")
        fspec_dirty, fspec_clean, components, fwhm = test_orig_rm(
            wave_sq, fpol_data, pixel
        )
        dirties.append(fspec_dirty)
        cleans.append(fspec_clean)

    fdirty = RMSynthesis.rm_synthesis(wave_sq, phi_range, fpol_data)

    fclean, components, fwhm = RMSynthesis.rm_clean(
        wave_sq, phi_range, fdirty, niter=500, gain=0.1, derotate=True, threshold=0
    )

    for _i, pixel in enumerate(pixels):
        print(f"Plotting pixel for original algorithm: {pixel}")
        debug_diagnostic_plots(
            dirties[_i],
            fdirty[:, *pixel],
            xaxis=phi_range,
            fname=f"fdirty-orig-vs-new-algo-{pixel[0], pixel[1]}",
        )
        debug_diagnostic_plots(
            cleans[_i],
            fclean[:, *pixel],
            xaxis=phi_range,
            fname=f"fclean-orig-vs-new-algo-px-{pixel[0], pixel[1]}",
        )

    print("----")
    # call_fitting(_v, wavelengths=wavelengths,
    #             max_depth=args.max_depth,
    #                 niters=args.niters, phi_step=args.depth_step,
    #                 derotate=args.no_derotate)

    # 0. Mask the data
    # 1. Get the fractional polarization
    # 2. Get the linear polarization
    # 3. Generate the:
    #    - RM map
    #    - p0 map
    #    - pA map
    #    - SNR map

    # create the new files
    new_fits_files = {
        ("p", "inherent-fractional-pol"): {
            "CTYPE3": "p",
            "CUNIT3": "ratio",
            "comment": "p0 at the peak of FDF",
        },
        ("pa", "pol-angle"): {
            "CTYPE3": "PA",
            "CUNIT3": "radians",
            "comment": "Polarization angle at the of peak FDF",
        },
        ("rm", "rotation-measure"): {
            "CTYPE3": "RM",
            "CUNIT3": "rad/m^2",
            "comment": "Faraday depth at the peak of FDF",
        },
        ("fpol", "fractional-pol"): {
            "CTYPE3": "FPOL",
            "CUNIT3": "x100%",
            "comment": "Fractional polarization at the central frequency",
        },
        ("lpol", "linear-pol"): {
            "CTYPE3": "LPOL",
            "CUNIT3": "unit",
            "comment": "Linear polarization at the central frequency",
        },
        ("snr", "signal-to-noise"): {
            "CTYPE3": "SNR",
            "CUNIT3": "unit",
            "comment": "Signal to noise at peak",
        },
    }

    ofiles = {}
    for (ctype, fname), hkeys in new_fits_files.items():
        _oname = OUTDIR.joinpath(f"{fname}-{ctype}").with_suffix(".fits")
        # hkeys = hkeys.update({"CUNIT3": hkeys["unit"]})
        clean_hdr.update(hkeys)
        ofiles[fname] = make_new_2d_fits(_oname, clean_hdr)

    print("Testing iteratives")


def console():
    """A console run entry point for setup.cfg"""
    main()
    # main2()
    snitch.info("---- Done. Bye :D ! ----")


if __name__ == "__main__":
    console()

    """
    Running

    time python pica_rm-x2.py -q Q-image-cubes.fits -u U-image-cubes.fits \
        -i I-image-cubes.fits -ncore 120 -o turbo \
        -f Frequencies-PicA-Masked.txt -mask true_mask_572.fits

    python pica_rm-x2.py -q Q-image-cubes.fits -u U-image-cubes.fits \
        -i I-image-cubes.fits -ncore 120 -o from_rick/outputs/with-NHS-data-mask \
        -f Frequencies-PicA-Masked.txt -mask masks/nhs-mask-572.fits

    python pica_rm-x2.py -q Q-image-cubes.fits -u U-image-cubes.fits \
        -i I-image-cubes.fits -ncore 120 -o from_rick/outputs/with-ricks-data-mask \
        -f Frequencies-PicA-Masked.txt -mask from_rick/masks/ricks-data-mask.fits

    """
