"""
Some useful equations from Brentjens 2005
Useful for some of the usual suspect calculations


Also contained a class for performing RM synthesis

In this code:

RM: rotation measure
SQ: squared
SQRT: sqaure root
"""

from dataclasses import dataclass
from typing import Iterable, Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import c

from scrappy.logging.config import configure_logger

matplotlib.use("Agg")

snitch = configure_logger(level="DEBUG")

# Get the value from this quantity
SPEED_OF_LIGHT = c.value

type FloatOrFloats = Union[float, Iterable[float]]


def simple_lsq(x):
    # a simple conversion from frequency to lambda squared
    return (SPEED_OF_LIGHT / x) ** 2


@dataclass
class RMResult:
    phi_range: np.ndarray
    fclean: np.ndarray
    fdirty: np.ndarray
    fwhm: float
    peak_phi: float
    peak_amp: float


class RMSynthesis:
    """
    A group of methods used to perform rotation measure synthesis
    on a 3-d data cube.
    """

    @classmethod
    def debug_plot_shift(cls, rmsf_shifts, xaxis=None, skip=50, prefix=None):
        if prefix is None:
            prefix = "rmsf-shift-phi"

        for i in range(0, len(rmsf_shifts), skip):
            plt.close("all")
            plt.plot(xaxis, rmsf_shifts[i])
            print("plotting rmsf-shift: ", i)
            plt.savefig(f"{prefix}-{i}.png")
        return

    @classmethod
    def rm_synthesis(
        cls, wave_sq: np.array, phi_range: np.array, p_lambda: np.ndarray | float
    ) -> np.ndarray:
        """
        Computes Faraday Spectra using RM-Synthesis
        as defined by Brentjens and de Bruyn (2005)
        A linear algebraic implementation

        Parameters
        ----------
        wave_sq: 1-d array
            wavelength squared
        phi_range: 1-d array
            Faraday depth range for which we are searching for a signal
        p_lambda: 3-d array
            Polarised/(fractional polarised) intensity data cube of shape (freq, ypix, xpix)

        Returns
        -------
        fdata : 3-d array
            A dirty Faraday spectrum of shape(phi_range.size, ypix, xpix)
        """

        n_waves = len(wave_sq)
        wave_sq0 = wave_sq.mean()

        # Shape: (len(phi_range), len(wave_sq))
        kernels = np.exp(-2j * np.outer(phi_range, wave_sq - wave_sq0))

        if isinstance(p_lambda, np.ndarray) and p_lambda.ndim > 2:
            _, nx, ny = p_lambda.shape
            p_lambda = p_lambda.reshape(n_waves, -1)

            fdata = (kernels @ p_lambda / n_waves).reshape(phi_range.size, nx, ny)
        else:
            fdata = kernels @ p_lambda / n_waves

        return fdata

    # @classmethod
    # def _defunct_rm_clean0(cls, wave_sq, phi_range, dirty_fspectrum, niter=500, gain=0.1,
    #             derotate=True, threshold=None):

    #     fspec_dirty = dirty_fspectrum.copy()  # don't modify input in place

    #     # Restoring beam
    #     if derotate:
    #         fwhm = 3.8 / abs(wave_sq[0] - wave_sq[-1])
    #     else:
    #         fwhm = 2 / (wave_sq[-1] + wave_sq[0])
    #     sigma = fwhm / 2.35482
    #     gaussian_norm = np.exp(-0.5 * (phi_range / sigma) ** 2)

    #     # Auto threshold
    #     if threshold is None:
    #         threshold = 1 * np.std(np.abs(fspec_dirty))

    #     # # Precompute RMSF shifts
    #     pad = abs(phi_range[-1]) * 2
    #     dpad = abs(phi_range[0] - phi_range[1])
    #     phi_pad = np.arange(-pad, pad+1, dpad)
    #     dshift = int(pad/(2.0 * dpad))

    #     rmsf_fixed = cls.rm_synthesis(wave_sq, phi_pad,
    #         np.ones(len(wave_sq), dtype=complex))

    #     nx, ny = dirty_fspectrum.shape[1:]
    #     ix, iy = np.indices((nx, ny))
    #     components = np.zeros_like(dirty_fspectrum, dtype=complex)

    #     for n in range(niter):
    #         f_amp = np.abs(fspec_dirty)
    #         peak = f_amp.max()

    #         # if peak < threshold:
    #         #     print(f"Converged after {n} iterations")
    #         #     break

    #         ind = np.argmax(f_amp, axis=0)              # (5, 5)
    #         f_comp = fspec_dirty[ind, ix, iy] * gain    # (5, 5)
    #         components[ind, ix, iy] += f_comp

    #         # Replicate original exactly, per pixel
    #         for xi in range(nx):
    #             for yi in range(ny):
    #                 dirac = np.zeros(len(phi_range))
    #                 dirac[ind[xi, yi]] = 1
    #                 rmsf = signal.convolve(rmsf_fixed, dirac, mode='same')
    #                 rmsf = rmsf[dshift:-dshift+1]
    #                 fspec_dirty[:, xi, yi] -= f_comp[xi, yi] * rmsf

    #     fspec_clean = fspec_dirty + (
    #         convolve1d(components.real, gaussian_norm, axis=0) +
    #         1j * convolve1d(components.imag, gaussian_norm, axis=0))

    #     return fspec_clean, components, fwhm

    @classmethod
    def rm_clean(
        cls,
        lam2,
        phi_range,
        dirty_fspectrum,
        niter=500,
        gain=0.1,
        derotate=True,
        threshold=0.0,
    ):
        """
        lam2:             (n_chan,)
        phi_range:        (n_phi,)  uniformly spaced
        fspec_dirty:  (n_phi,) or (n_phi, nx, ny)
        returns:          fclean, components (same shape as input), fwhm
        """
        fspec_dirty = dirty_fspectrum.copy()  # don't modify input in place

        shape = fspec_dirty.shape
        n_phi = shape[0]

        # Work on a copy, flattened to (n_phi, n_pixels)
        fres = np.array(fspec_dirty, dtype=complex).reshape(n_phi, -1)
        pix = np.arange(fres.shape[1])

        # Restoring beam width
        if derotate:
            fwhm = 3.8 / abs(lam2[0] - lam2[-1])
        else:
            fwhm = 2 / (lam2[-1] + lam2[0])
        sigma = fwhm / 2.35482

        # RMSF on offsets -(n_phi-1) ... +(n_phi-1): exactly 2*n_phi - 1 points
        dphi = phi_range[1] - phi_range[0]
        phi_pad = dphi * np.arange(-(n_phi - 1), n_phi)
        rmsf_fixed = cls.rm_synthesis(
            lam2, phi_pad, np.ones(len(lam2), dtype=complex)
        )
        center = n_phi - 1  # index of the RMSF peak

        # rmsf_cols[:, i] = RMSF centred on phi_range[i], sampled on phi_range
        rmsf_cols = np.array(
            [rmsf_fixed[center - i : center - i + n_phi] for i in range(n_phi)]
        ).T.copy()  # (n_phi, n_phi)

        # cls. debug_plot_shift(rmsf_cols, xaxis=phi_range, skip=50, prefix="rm-shifts-new")
        components = np.zeros_like(fres)

        if threshold is None:
            threshold = 1 * np.std(np.abs(fspec_dirty))

        for n in range(niter):
            f_amp = np.abs(fres)
            ind = np.argmax(f_amp, axis=0)  # peak index per pixel
            active = f_amp[ind, pix] > threshold  # pixels still being cleaned

            if not active.any():
                print(f"Converged after {n} iterations")
                break

            f_comp = np.where(active, fres[ind, pix] * gain, 0)
            components[ind, pix] += f_comp
            fres -= rmsf_cols[:, ind] * f_comp  # (n_phi, n_pixels)

        # Restore: peak-1 Gaussian centred on each component, plus residuals
        restore = np.exp(
            -0.5 * ((phi_range[:, None] - phi_range[None, :]) / sigma) ** 2
        )
        fclean = restore @ components + fres

        return fclean.reshape(shape), components.reshape(shape), fwhm

        @classmethod
        def rmsynth_pipeline(
            cls,
            wave_sq,
            pol_data,
            phi_max=500,
            phi_step=5,
            niter=1000,
            gain=0.1,
            derotate=True,
            threshold=None,
        ):
            # Validate
            assert len(wave_sq) == len(
                pol_data
            ), "wave_sq and pol_data must have same length"
            assert np.iscomplexobj(pol_data), "pol_data must be complex (Q + iU)"

            # Build phi axis (guaranteed zero at centre)
            phi_range = np.arange(-phi_max, phi_max + phi_step, phi_step)

            # RM synthesis
            fdirty = cls.rm_synthesis(wave_sq, phi_range, pol_data)

            # RMCLEAN
            fclean, fcomp, fwhm = cls.rmclean(
                wave_sq,
                phi_range,
                fdirty,
                niter=niter,
                gain=gain,
                derotate=derotate,
                threshold=threshold,
            )

            # Peak
            peak_ind = np.argmax(np.abs(fclean))

            return RMResult(
                phi_range=phi_range,
                fclean=fclean,
                fdirty=fdirty,
                fwhm=fwhm,
                peak_phi=phi_range[peak_ind],
                peak_amp=np.abs(fclean[peak_ind]),
            )


def chan_width_wavelengths(center_freq: float, chan_width: float) -> float:
    """Equation 35
    Channel width in wavelengths
    """
    cww = ((2 * SPEED_OF_LIGHT**2 * chan_width) / center_freq**3) * (
        1 + 0.5 * (chan_width / center_freq) ** 2
    )
    return cww


def freq_sqrt(lambda_sq: FloatOrFloats) -> FloatOrFloats:
    lambda_sqrt = np.sqrt(lambda_sq)
    # frequency to wavelength
    freq_ghz = SPEED_OF_LIGHT / lambda_sqrt
    snitch.debug("Returning frequency in GHz")
    return freq_ghz


def lambda_sqd(center_freq: float, chan_width: float) -> float:
    """Equation 34
    Calculate the center frequency in lambda squared for a single channel
    """
    lsq = (SPEED_OF_LIGHT**2 / center_freq**2) * (
        1 + (0.75 * (chan_width / center_freq) ** 2)
    )
    return lsq


def lambda_sqd_err(
    n_chans: int, lambda_sq_0: float, lambda_sq: Iterable[float]
) -> float:
    """Equation 53
    Error in the lambda squared
    """
    lambda_sum = (lambda_sq**2 - lambda_sq_0**2).sum()

    assert (
        n_chans == lambda_sq.size
    ), f"Provided value of n_chans {n_chans} != size of wavelengths squared provided: {lambda_sq.size}"
    err_lambda = np.sqrt(lambda_sum / (n_chans - 1))
    return err_lambda


def linear_polzn(stokes_Q: FloatOrFloats, stokes_U: FloatOrFloats) -> FloatOrFloats:
    """complex linear polarised intensity"""
    linear = stokes_Q + 1j * stokes_U
    return linear


def linear_polzn_error(
    stokes_Q: FloatOrFloats,
    stokes_U: FloatOrFloats,
    stokes_Q_err: FloatOrFloats,
    stokes_U_err: FloatOrFloats,
) -> FloatOrFloats:
    """Appendix A.7
    Error in total linear polarised intensity
    """
    lpol = np.abs(linear_polzn(stokes_Q, stokes_U))
    err_sq = (stokes_Q**2 * stokes_Q_err**2) / lpol**2 + (
        stokes_U**2 * stokes_U_err**2
    ) / lpol**2
    return np.sqrt(err_sq)


def polzn_angle(
    stokes_Q: FloatOrFloats, stokes_U: FloatOrFloats, deg=False
) -> FloatOrFloats:
    """Angle of polarisation"""
    angle = 0.5 * np.arctan2(stokes_U, stokes_Q)
    if deg:
        angle = np.rad2deg(angle)
    return angle


def polzn_angle_error(
    stokes_Q: FloatOrFloats,
    stokes_U: FloatOrFloats,
    stokes_Q_err: FloatOrFloats,
    stokes_U_err: FloatOrFloats,
) -> FloatOrFloats:
    """Equation 54 and 55"""
    lpol = np.abs(linear_polzn(stokes_Q, stokes_U))
    err_sq = ((stokes_U**2 * stokes_Q_err**2) + (stokes_Q**2 * stokes_U_err**2)) / (
        4 * lpol**4
    )
    return np.sqrt(err_sq)


def rm_error(pangle_err: float, lam_sq_err: float, n_chans: int):
    """Equation 52 and 53"""
    err_rm = pangle_err / (lam_sq_err * np.sqrt(n_chans - 2))
    return err_rm


def frac_polzn(
    stokes_I: FloatOrFloats, stokes_Q: FloatOrFloats, stokes_U: FloatOrFloats
) -> FloatOrFloats:
    """
    Note that the input stokes_Q, stokes_U are not fractional at this point
    """
    # linear fractional polarisation
    frac_pol = np.abs(linear_polzn(stokes_Q / stokes_I, stokes_U / stokes_I))
    return frac_pol


def frac_polzn_error(
    stokes_I: FloatOrFloats,
    stokes_Q: FloatOrFloats,
    stokes_U: FloatOrFloats,
    stokes_I_err: FloatOrFloats,
    stokes_Q_err: FloatOrFloats,
    stokes_U_err: FloatOrFloats,
) -> FloatOrFloats:
    """
    Equation
    """
    fpol = frac_polzn(stokes_I, stokes_Q, stokes_U)
    p = np.abs(linear_polzn(stokes_Q / stokes_I, stokes_U / stokes_I))
    p_err = linear_polzn_error(stokes_Q, stokes_U, stokes_Q_err, stokes_U_err)
    res = np.abs(fpol) * np.sqrt(
        np.square(p_err / p) + np.square(stokes_I_err / stokes_I)
    )
    return res


# TODO: NEW implementation. Check that it works conventionally! Delete previous
def new_frac_polzn_error(
    stokes_I: FloatOrFloats,
    stokes_Q: FloatOrFloats,
    stokes_U: FloatOrFloats,
    stokes_I_err: FloatOrFloats,
    stokes_Q_err: FloatOrFloats,
    stokes_U_err: FloatOrFloats,
) -> FloatOrFloats:
    """
    From error propagation, GPT and ClaudeAI
    """
    p = frac_polzn(stokes_I, stokes_Q, stokes_U)

    numerator_QU = (stokes_Q**2) * (stokes_Q_err**2) + (stokes_U**2) * (
        stokes_U_err**2
    )

    # Calculate the error squared
    sigma_p_sq = (numerator_QU) / (stokes_I**2 * (stokes_Q**2 + stokes_U**2)) + (
        (p * stokes_I_err) / stokes_I
    ) ** 2

    return np.sqrt(sigma_p_sq)


######################################
# some checks here
######################################


def polarised_snr(
    stokes_Q: FloatOrFloats,
    stokes_U: FloatOrFloats,
    stokes_Q_err: FloatOrFloats,
    stokes_U_err: FloatOrFloats,
) -> FloatOrFloats:
    """
    From Feain 2012: Section 5
    Signal-to-noise ratio of the polarised signal
    Taken to be the ratio between the:
    1. Total polarised intensity
    2. Standard error of (1)
    """
    lpol = np.abs(linear_polzn(stokes_Q, stokes_U))
    lpol_err = linear_polzn_error(stokes_Q, stokes_U, stokes_Q_err, stokes_U_err)
    return lpol / lpol_err


######################################
# RM SYNTHESIS some values
######################################


def max_faraday_depth(center_freq: float, chan_width: float) -> float:
    """
    Equation 63
    Approximate max faraday depth at which one has more than 50% sensitivity
    """
    cw_wave = chan_width_wavelengths(center_freq, chan_width)
    mfd = 3**0.5 / cw_wave
    return mfd


def max_faraday_scale(lambda_sq_min: float) -> float:
    """
    Equation 62
    Scale in depth space where sensitivty is at 50%
    """
    mscale = np.pi / lambda_sq_min
    return mscale


def rmtf_resolution(lambda_sq_min: float, lambda_sq_max: float) -> float:
    """
    Equation 61
    Approximate FWHM of the RMSF

    But changed by schnitzeler to 3.8/lambda^2
    We use schnitzelers version
    """
    # msens = (2 * 3**0.5)/(lambda_sq_max - lambda_sq_min)
    msens = 3.8 / (lambda_sq_max - lambda_sq_min)
    return msens


def max_chan_width(max_rm, max_rot, obs_freq) -> float:
    """
    Get the maximum allowable channel width before bandwidth depolzn
    takes over. Sebokolodi 2020 Eq 4

    max_rm: float
        Maximum known RM of source / maximum fdepth possible
    max_rot: float
        Maximum allowable rotation per channel in degrees
    """
    max_rot = np.deg2rad(max_rot)
    wave_sq = simple_lsq(obs_freq)
    max_cw = (max_rot * obs_freq) / (2 * wave_sq * max_rm)
    print(f"{max_cw / 1e6:.2f} MHz")
    return max_cw
