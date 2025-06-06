"""
Some useful equations from Brentjens 2005
Useful for some of the usual suspect calculations


In this code:

RM: rotation measure
SQ: squared
SQRT: sqaure root
"""
import numpy as np
from astropy.constants import c
from scrappy.logging.config import setup_logging


# Get the value from this quantity
SPEED_OF_LIGHT = c.value
snitch = setup_logging()

def chan_width_wavelengths(center_freq: float, chan_width: float):
    """Equation 35
    Channel width in wavelengths
    """
    cww = ((2 * SPEED_OF_LIGHT**2 * chan_width)/center_freq**3) \
        * (1 + 0.5 * (chan_width/center_freq)**2)
    return cww

def freq_sqrt(lambda_sq):
    lambda_sqrt = np.sqrt(lambda_sq)
    # frequency to wavelength
    freq_ghz = (SPEED_OF_LIGHT/lambda_sqrt)
    snitch.debug("Returning frequency in GHz")
    return freq_ghz


simple_lsq = lambda x: (SPEED_OF_LIGHT/x)**2


def lambda_sqd(center_freq: float, chan_width: float):
    """Equation 34
    Calculate the center frequency in lambda squared for a single channel
    """
    lsq = (SPEED_OF_LIGHT**2 / center_freq**2) * \
        (1 + ( 0.75 * (chan_width/center_freq)**2))
    return lsq

def lambda_sqd_err(n_chans: int, lambda_sq_0: float, lambda_sq: tuple):
    """Equation 53
    Error in the lambda squared
    """
    lambda_sum = (lambda_sq**2 - lambda_sq_0**2).sum()

    assert n_chans == lambda_sq.size, f""""Provided value of n_chans {n_chans} != size of wavelengths squared providied: {lambda_sq.size}"""
    err_lambda = np.sqrt(lambda_sum/(n_chans - 1))
    return err_lambda

def linear_polzn(stokes_Q: float, stokes_U: float):
    """complex linear polarised intensity"""
    linear = stokes_Q + 1j*stokes_U
    return linear

def linear_polzn_error(stokes_Q: float, stokes_U: float, stokes_Q_err: float, stokes_U_err: float):
    """Appendix A.7
    Error in total linear polarised intensity
    """
    lpol = np.abs(linear_polzn(stokes_Q, stokes_U))
    err_sq = (stokes_Q**2 * stokes_Q_err**2)/ lpol**2 + \
             (stokes_U**2 * stokes_U_err**2)/lpol**2
    return np.sqrt(err_sq)

def polzn_angle(stokes_Q: float, stokes_U: float, deg=False):
    """Angle of polarisation"""
    angle = 0.5 * np.arctan2(stokes_U,stokes_Q)
    if deg:
        angle = np.rad2deg(angle)
    return angle

def polzn_angle_error(stokes_Q: float, stokes_U: float, stokes_Q_err: float, stokes_U_err: float):
    """Equation 54 and 55"""
    lpol = np.abs(linear_polzn(stokes_Q, stokes_U))
    err_sq = ((stokes_U**2 * stokes_Q_err**2) + (stokes_Q**2 * stokes_U_err**2)) \
            / (4 * lpol**4)
    return np.sqrt(err_sq)

def rm_error(pangle_err: float, lam_sq_err: float, n_chans: int):
    """Equation 52 and 53"""
    err_rm = pangle_err / (lam_sq_err * np.sqrt(n_chans - 2))
    return err_rm

def frac_polzn(stokes_I: float, stokes_Q: float, stokes_U: float):
    """
    Note that the input stokes_Q, stokes_U are not fractional at this point
    """
    # linear fractional polarisation
    frac_pol = np.abs(linear_polzn(stokes_Q/stokes_I, stokes_U/stokes_I))
    return frac_pol



def frac_polzn_error(stokes_I, stokes_Q, stokes_U, stokes_I_err, stokes_Q_err, stokes_U_err):
    """
    Equation 
    """
    fpol = frac_polzn(stokes_I, stokes_Q, stokes_U)
    p = np.abs(linear_polzn(stokes_Q/stokes_I, stokes_U/stokes_I))
    p_err = linear_polzn_error(stokes_Q, stokes_U, stokes_Q_err, stokes_U_err)
    res = np.abs(fpol) * np.sqrt(np.square((p_err/p)) + np.square((stokes_I_err/stokes_I)))
    return res


# TODO: NEW implementation. Check that it works conventionally! Delete previous
def new_frac_polzn_error(stokes_I, stokes_Q, stokes_U, stokes_I_err, stokes_Q_err, stokes_U_err):
    """
    From error propagation, GPT and ClaudeAI
    """
    p = frac_polzn(stokes_I, stokes_Q, stokes_U)

    numerator_QU = (stokes_Q**2) * (stokes_Q_err**2) + (stokes_U**2) * (stokes_U_err**2)
    
    # Calculate the error squared
    sigma_p_sq = (numerator_QU) / (stokes_I**2 * (stokes_Q**2 + stokes_U**2)) + ((p * stokes_I_err) / stokes_I)**2

    return np.sqrt(sigma_p_sq)



######################################
# some checks here
######################################


def polarised_snr(stokes_Q: float, stokes_U: float, stokes_Q_err: float, stokes_U_err: float):
    """
    From Feain 2012: Section 5
    Signal-to-noise ratio of the polarised signal
    Taken to be the ratio between the:
    1. Total polarised intensity
    2. Standard error of (1)
    """
    lpol = np.abs(linear_polzn(stokes_Q, stokes_U))
    lpol_err = linear_polzn_error(stokes_Q, stokes_U, stokes_Q_err, stokes_U_err)
    return lpol/lpol_err



######################################
# RM SYNTHESIS some values
######################################


def max_faraday_depth(center_freq: float, chan_width: float):
    """
    Equation 63
    Approximate max faraday depth at which one has more than 50% sensitivity
    """
    cw_wave = chan_width_wavelengths(center_freq, chan_width)
    mfd = 3**0.5 / cw_wave
    return mfd

def max_faraday_scale(lambda_sq_min: float):
    """
    Equation 62
    Scale in depth space where sensitivty is at 50%
    """
    mscale = np.pi / lambda_sq_min
    return mscale

def rmtf_resolution(lambda_sq_min: float, lambda_sq_max: float):
    """
    Equation 61
    Approximate FWHM of the RMSF

    But changed by schnitzeler to 3.8/lambda^2
    We use schnitzelers version
    """
    # msens = (2 * 3**0.5)/(lambda_sq_max - lambda_sq_min)
    msens = 3.8/(lambda_sq_max - lambda_sq_min)
    return msens


def max_chan_width(max_rm, max_rot, obs_freq):
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
    max_cw = (max_rot * obs_freq) /(2 * wave_sq * max_rm)
    print(f"{max_cw/1e6:.2f} MHz")
    return max_cw

