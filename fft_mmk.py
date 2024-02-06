import numpy as np
from scipy.fft import fft, fftfreq, fftshift

def fft_output(signal, dt):

    """
    Calculates the power spectrum of a signal using the FFT algorithm.

    Parameters:
    signal (array-like): The input signal.
    dt (float): The time interval between samples.

    Returns:
    periods (ndarray): The periods corresponding to the frequencies.
    power_spec (ndarray): The power spectrum calculated as squares of the absolute value of the Fourier coefficients.
    """

    # Check inputs
    signal = np.asarray(signal)
    if signal.ndim != 1:
        raise ValueError("Input signal must be 1D.")
        
    # Take fourier transform
    fft_coeffs = fft(signal)

    # Generate frequency axis
    xf = fftfreq(len(signal), dt)

    # isolate indices of negative frequencies
    index = np.where(xf < 0)

    # create new array with coefficients corresponding to negative frequencies ONLY
    c_negative = fft_coeffs[index]
    print("c_negative ", c_negative, "\n")

    # create new array with negative frequencies ONLY, and flip to give positive values for plotting purposes
    freq_negative = -xf[index]
    print("freq_negative ", freq_negative, "\n")

    # generate periods from frequencies
    periods = 1/freq_negative
    print("periods ", periods, "\n")

    # calculate power spectrum 
    power_spec = np.abs(c_negative)**2
    
    return periods, power_spec

