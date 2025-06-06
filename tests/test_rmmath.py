import pytest
import math

import scrappy.utils.rmmath as rm
import numpy as np

@pytest.mark.parametrize(
    "central_freq_of_chan, chan_width_freq", [
        (1e9, 1e6)
    ])
def test_chan_width_wavelengths(central_freq_of_chan, chan_width_freq):
    result = rm.chan_width_wavelengths(central_freq_of_chan, chan_width_freq)

    chan_width_in_wavelength_sq = (central_freq_of_chan**-5) * \
        ((central_freq_of_chan**2 * 2 * rm.SPEED_OF_LIGHT**2 * chan_width_freq) \
            + (rm.SPEED_OF_LIGHT**2 * chan_width_freq**3))

    assert chan_width_in_wavelength_sq == pytest.approx(result)



@pytest.mark.parametrize("lambda_sq", [(0.09, 0.009, 120)])
def test_freq_sqrt(lambda_sq):
    freq_orig = rm.freq_sqrt(lambda_sq)
    frequency = np.sqrt((rm.SPEED_OF_LIGHT**2/lambda_sq))
    # check that we get the frequency
    assert freq_orig == pytest.approx(frequency)

    # testing conversion to lambda sqrd
    lambda_sq_test = rm.simple_lsq(freq_orig)
    assert  lambda_sq == pytest.approx(lambda_sq_test)



@pytest.mark.parametrize(
    "center_freq, chan_width", [
        (1e9, 1e6)
    ])
def test_lambda_sq(center_freq, chan_width):
    lambda_sq_orig = rm.lambda_sqd(center_freq, chan_width)
    lambda_sq_test = (rm.SPEED_OF_LIGHT**2 * ( (4 * center_freq**2) + 3*(chan_width**2))) / (4*center_freq**4)
    assert  lambda_sq_orig == pytest.approx(lambda_sq_test)



@pytest.mark.parametrize(
    "n_chans, lambda_sq_0, lambda_sq",
    [
        (4096, 0.09, np.linspace(0.09, 0.5, 4096))
    ])
def test_lambda_sqd_err(n_chans, lambda_sq_0, lambda_sq):
    # getting variance
    sigma_sq_test = np.sum(lambda_sq**2 - lambda_sq_0**2) / (n_chans - 1)
    sigma_test = np.sqrt(sigma_sq_test)
    sigma_orig = rm.lambda_sqd_err(n_chans, lambda_sq_0, lambda_sq)
    sigma_sq_orig = sigma_orig**2

    assert sigma_orig == pytest.approx(sigma_test)
    assert sigma_sq_orig == pytest.approx(sigma_sq_test)


@pytest.mark.parametrize(
    "stokes_Q, stokes_U, stokes_Q_err, stokes_U_err",[
        (0.5, 0.3, 2e-4, 2.05e-4)
    ]
)
def test_linear_polzn_error(stokes_Q, stokes_U, stokes_Q_err, stokes_U_err):
    # if sigma_q ~ sigma_u == sigma, then sigma_p^2 == sigma^2
    sigma_p = rm.linear_polzn_error(stokes_Q, stokes_U, stokes_Q_err, stokes_U_err)
    assert pytest.approx(sigma_p, rel=1e4) == stokes_Q_err




@pytest.mark.parametrize(
    "stokes_Q, stokes_U, deg",[
        (0.5, 0.8, False),
        (0.5, 0.8, True)
    ]
)
def test_polzn_angle(stokes_Q, stokes_U, deg):
    angle = rm.polzn_angle(stokes_Q, stokes_U, deg)
    if deg:
        angle = np.deg2rad(angle)
    assert pytest.approx(np.tan(2*angle)) == (stokes_U/stokes_Q)



@pytest.mark.parametrize(
    "stokes_Q, stokes_U, stokes_Q_err, stokes_U_err",[
        (0.5, 0.3, 2e-4, 2.05e-4)
    ]
)
def test_polzn_angle_error(stokes_Q, stokes_U, stokes_Q_err, stokes_U_err):
    # equations say if sigma_q ~ sigma_u = sigma (testing using eq 54)
    sigma_angle_orig = rm.polzn_angle_error(stokes_Q, stokes_U, stokes_Q_err, stokes_U_err)
    sigma_angle_test = 0.5 * stokes_Q_err/(np.abs(stokes_Q + stokes_U*1j))
    assert sigma_angle_test == pytest.approx(sigma_angle_orig, rel=1e5)



@pytest.mark.math
@pytest.mark.parametrize(
    "stokes_I, stokes_Q, stokes_U, stokes_I_err, stokes_Q_err, stokes_U_err",[
        (6.5, 0.5, 0.3, 2e-2, 2e-4, 2.05e-4)
    ]
)
def test_frac_polzn_error(stokes_I, stokes_Q, stokes_U, stokes_I_err, stokes_Q_err, stokes_U_err):
    fpol_sigma_test_sq = ((1 / (stokes_I**2 * (stokes_Q**2 + stokes_U**2))) * (stokes_Q**2 * stokes_Q_err**2 + stokes_U**2 * stokes_U_err**2)) + ((rm.frac_polzn(stokes_I,stokes_Q,stokes_U)/stokes_I)*stokes_I_err)**2
    
    fpol_sigma_orig = rm.new_frac_polzn_error(stokes_I, stokes_Q, stokes_U, stokes_I_err, stokes_Q_err, stokes_U_err)
    
    assert fpol_sigma_orig == pytest.approx(np.sqrt(fpol_sigma_test_sq))