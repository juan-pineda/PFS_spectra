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

def sky_coord_to_galactic(xcen, ycen, pos_angl, incl, im_size=(240, 240)):
    """
	Compute the coordinates in the galactic plane associated to each pixel in
	the plane of the sky

    Parameters
    ----------
    xcen : float
		position of the center of the models in pixels
    ycen : float
		position of the center in the models in pixels
    pos_angl : float
		position angle of the major axis in degrees
	incl : float
		inclination of the disk in degrees
	im_size: tuple
		shape of the model arrays to be created
    
    Returns
	-------
	r : array
		real radius in pixels from the center, of a projected pixel in the sky
	ctheta : array
		cosine of the azimuthal angle in the plane of the disk, associated to
		each píxel in the projection on plane of the sky

    """

    y, x = np.indices(im_size)
    den = (y - ycen) * math.cos(math.radians(pos_angl)) - (x - xcen) * math.sin(math.radians(pos_angl))
    num = - (x - xcen) * math.cos(math.radians(pos_angl)) - (y - ycen) * math.sin(math.radians(pos_angl))
    r = (den ** 2 + (num / math.cos(math.radians(incl))) ** 2) ** 0.5
    tpsi = num * 1.

    tpsi[np.where(den != 0)] /= den[np.where(den != 0)]  # to avoid a NaN at the center
    den2 = math.cos(math.radians(incl)) ** 2 + tpsi ** 2
    sg = np.sign(den)  # signe
    ctheta = sg * (math.cos(math.radians(incl)) ** 2 / den2) ** 0.5  # azimuth in galaxy plane

    return [r, ctheta]


def next_odd(x):
    """
	Compute next odd integer larger or equal to a given number

    """

    x = np.ceil(x)
    if x % 2 == 1:
        return x
    else:
        return x+1

def spatial_convolution(array,psf):
    """
	Spatial convolution of an array with a gaussian kernel of given sigma (psf)    

    """
	# 20 is because I am getting 10 times sigma in each direction
    side = next_odd(20*psf)
    psf = astropy.convolution.Gaussian2DKernel(psf, x_size=side, y_size=side)
    new_array = astropy.convolution.convolve(array,psf)
    return new_array


def spatial_convolution_fft(array,psf):
    side = next_odd(20*psf)
    psf = astropy.convolution.Gaussian2DKernel(psf, x_size=side, y_size=side)
    new_array = astropy.convolution.convolve_fft(array, psf, psf_pad=True, normalize_kernel=np.sum, allow_huge=True)
    return new_array



