import numpy as np
import tools

# This file contains the models for the intensity maps.

# Each block corresponds to a 1D light profile. They in general depend on 3
# parameters: a central brigthness (I0), a characteristic scale length (rdf),
# and a truncation radius (rt). 

# The last block evaluates the 2D intensity map using one of the radial profiles
# and the specific parameters passed in the configfile.


def flat_disk_intensity(r, I0, rdf, rt):
	"""
	Model with constant surface brigthness 

	Parameters
	----------
	r : float, array
		radial positions where the model is to be computed
	I0 : float 
		nuisance parameter in this case
	rdf : float
		nuisance parameter in this case
	rt : float
		truncation radius

	Returns
	-------
	Array with the same shape of r, containing the model intensity profile/map
	
	"""

	flux[(r <= rtrunc)] = I0
	flux[(r > rtrunc)] = 0.

	return flux


def exponential_disk_intensity(r, I0, rdf, rt):
	"""
	Surface brigthness profile of an exponential disc
	r, rdf, and rt must be in the same units
	
	Parameters
	----------
	r : array
		radial positions where the model is to be computed
	I0 : float 
		central surface brigthness in arbitrary units
	rdf : float
		scale-length of the disc
	rt : float
		truncation radius
	
	Returns
	-------
	Array with the same shape of r, containing the model intensity profile/map
	
	"""

	flux = I0 * np.exp(-r/rdf)
	flux[(r > rt)] = 0.

	return flux


	
	
	##--- compute the 2D velocity model with the specific desired geometry ---##

def make_fluxmap(geom, models, instr):
	"""
	Creates a 2D model intensity map using the parameters in the configuration
	file  via the objects *geom*, *models*, *instr*.
	
	For a detailed description of the configuration file and the associated
	objects, see the file read_params.py.

	"""

	# Compute the maps of the real radius and the cosine of azymutal angle theta,
	# corresponding to every pixel in the plane of the disc
	r,ctheta = tools.sky_coord_to_galactic(geom.x0, geom.y0, geom.pa, geom.incl,
			im_size=(instr.size, instr.size))
	# Map of the flux at any given radius
	if models.fm == "exp":
		flux = exponential_disk_intensity(r, models.i0, models.rdf, models.rt)
	elif models.fm == "flat":
		flux = flat_disk_intensity(r, models.i0, models.rdf, models.rt)

	# this is optional as we are using arbitrary units by now
	flux = flux / np.sin(np.deg2rad(geom.incl))
	
	return flux


