import numpy as np
import xarray as xr

def calculate_spectrum_kx_ky(ds):
    # Requires dataset to be on a regular lat lon grid
    
    # Perform fft
    fft_result = np.fft.fft2(ds)

    # Create spectral grid
    dx = ds.lon[1] - ds.lon[0]
    dy = ds.lat[1] - ds.lat[0]
    kx = np.fft.fftfreq(len(ds.lon), d=dx.values)
    ky = np.fft.fftfreq(len(ds.lat), d=dy.values)

    # Create spectrum xr DataArray
    ds_spectrum = xr.DataArray(fft_result,
                               dims=['ky','kx'],
                               coords=({'ky':ky, 'kx':kx}))

    return ds_spectrum

def save_spectrum(ds_spectrum, path, filename):
    ds = xr.Dataset({'real': ds_spectrum.real, 'imag': ds_spectrum.imag}).to_netcdf(path + filename)
    