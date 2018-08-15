import numpy as np
import math
from bisect import bisect_left
import astropy.convolution
from scipy import special
import flux_model as fm
import velocity_model as vm
import read_params as rp
import tools
from astropy import constants as const


def observed_lambda(geom, models, instr, v):
    """
	Calculate the doppler shift due to cosmological redshift and to the proper
	motion of each element in the galaxy. It returns the redshifted position
	of the central wavelength and the corresponding wavelegth dispersion
	(broadening), including the effect of the spectral resolution.

    Parameters
    ----------
    v : float, array
		can be a single value or an array, e.g., a velocity map

    Returns
    -------
	lambda_obs : float or array with the same shape of v
		central position of the redshifted emission lines
	sigma : float or array with the same shape of v
		1-sigma broadening of the gaussian emission lines redshifted

    """

    # Central wavelength of the redshifted emission lines at each pixel
    lambda_obs = models.l0 * (1+geom.z) * (1 + v/const.c.to('km/s').value)
    # Broadening of the detected emission lines
    sigma_lambda = models.l0 * (1+geom.z) * models.disp / const.c.to('km/s').value
    sigma_lambda = sigma_lambda * np.ones(lambda_obs.shape)
    # Effect of the spectral resolution
    sigma_res = models.l0 * (1+geom.z) / instr.R / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    # Combined effects
    sigma = np.sqrt(sigma_lambda**2 + sigma_res**2)
    return lambda_obs, sigma


def psf_footprint(geom, instr):
    """
	Perform the spatial convolution of the psf with a top-hat function having
	the shape of the fiber. This is to compute in advance the fraction of the
	flux at each pixel that will actually make it to be collected by the fiber
	once the spatial convolution is taken into account.

	If the PSF is none, or if it is cero, then all the light in the footprint
	of the fiber is considered, and nothing else

    Returns
    -------
	psf_mask : array
		footprint of the fiber psf-convolved, with the same shape of the input
		maps

    """

    # Identify those pixels which are actually covered by the fiber
    y, x = np.indices((instr.size, instr.size))
    # center the fiber mask in the center of the map
    R = np.sqrt((x - (instr.size/2 - 0.5))**2 + (y - (instr.size/2 - 0.5))**2)
    mask = (R <= instr.fiber)
    if (instr.psf==0) | (instr.psf==None):
        return mask
	# We can have more than one convolution scheme here. Currently there are
	# two schemes implemented in tools.py
    psf_mask = tools.spatial_convolution_fft(mask,instr.psf)
    return psf_mask


def make_spectrum(instr, flux, lambda_obs, sigma, chunk = 50):
    """
	This is the main function, it takes every pixel in the models and compute
	the center/width of the associated emission line

    Parameters
    ----------
	flux  : array
			map of the flux that  will be used to weight the contribution of
			each individual spectrum to the total. This map should have been
			already weighted by the effect of the PSF.
	lambda_obs : array
			central positions of the redshifted emission lines
    chunk : int
            Number of pixels to be processed at once is chunk**2. The larger
			the faster, but this number will be limited by the available RAM 

    """

	# vector to store the flux of the final spectrum 
    spectrum = np.zeros(instr.channels.size)
	# position of the first element of each chunk of pixels
    index = list(np.arange(0,instr.size,chunk))
	# to correctly consider the final position of the last chunk
    index.append(instr.size)
    for i in range(len(index)-1):
        for j in range(len(index)-1):
			# We will take a rectangular slice from the model maps. npix2 is the
			# number of pixels in that region
            npix2 = (index[i+1] - index[i]) * (index[j+1] - index[j])
			# number of channels in the spectrum
            n_ch = instr.channels.size
			# Replicates the centers of the channels into a rectangular grid
            ch_center = np.tile(instr.channels, (npix2,1))
			# constant width of the spectral channels, replicated to grid shape
            ch_width = (instr.channels[1] - instr.channels[0]) * np.ones(ch_center.shape)
			# array/slice of centers of the individual emission lines
            mu_array = np.transpose(np.tile(np.ravel(lambda_obs[index[i]:index[i+1],index[j]:index[j+1]]), (n_ch, 1)))
            # array/slice of individual wavelength dispersions
            sigma_arr = np.transpose(np.tile(np.ravel(sigma[index[i]:index[i+1],index[j]:index[j+1]]), (n_ch, 1)))
			# array/slice of the flux intensity map
            Ins = np.transpose(np.tile(np.ravel(flux[index[i]:index[i+1],index[j]:index[j+1]]), (n_ch, 1)))

            line_Ha = int_gaussian(ch_center, ch_width, mu_array, sigma_arr) * Ins
			# Next line was a test, to see what happens if we take the central
			# position of the flux inside a channel instead of actually
			# integrating the gaussian shape inside the limits of the channel
			# No difference was visible, though it was a rather shallow test
			# line_Ha     = (1./(sigma_arr * np.sqrt(2 * np.pi))) * np.exp(-(mu_array-ch_center)**2/2./sigma_arr**2) * Ins

            spectrum = spectrum + line_Ha.sum(axis=0)

    return spectrum


def int_gaussian(x,dx,mu,sigma):
    """
    Compute the integral of a normalized gaussian inside some limits.
    The error functions (special.erf) come from the analytical solution of the
    integral.

    Parameters
    ----------
    x : float, array
        center of the interval where the gaussian is to be integrated.
    dx : float, array
        width of the integration interval.
    mu: float, array
        mean of the gaussian.
    sigma: float, array
        standard deviation of the gaussian.

    """

    # This functions come from the analytic solution
    A = special.erf((x+dx/2-mu)/np.sqrt(2)/sigma)
    B = special.erf((x-dx/2-mu)/np.sqrt(2)/sigma)
    return np.abs((A-B)/2)




