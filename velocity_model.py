import numpy as np
import tools
from scipy.special import i0, i1, k0, k1

# This file contain the models for the velocity fields.

# Each block corresponds to a 1D radial velocity profile. They in general
# depend on 2 parameters: a characteristic scale length (rd) and maximum
# rotation velocity (vmax). 

# The last block evaluates the 2D velocity map using one of the radial profiles
# and the specific parameters passed in the configfile.


def exponential_velocity(r, rd, vmax):
    """
    Velocity function for an exponential profile
	r and rd must be in the same units

	Parameters
	----------
	r : array
		radial positions where the model is to be computed
    rd : float
		radius at which the maximum velocity is reached
    vmax : float
		Maximum velocity of the model

	Returns
	-------
	Array with the same shape of r, containing the model velocity curve/map

    """

	# disk scale length
    rd2 = rd / 2.15   
    vr = np.zeros(np.shape(r))
	# to prevent any problem in the center
    q = np.where(r != 0)      
    vr[q] = r[q] / rd2 * vmax / 0.88 * np.sqrt(i0(0.5 * r[q] / rd2) * k0(0.5 * r[q] / rd2) - i1(0.5 * r[q] / rd2) * k1(0.5 * r[q] / rd2))
    return vr


def flat_velocity(r, rt, vmax):
    """
    Velocity function for a flat profile. It is composed of a linearly rising
	inner part, and an horizontal part after a given radius
	r and rd must be in the same units

    Parameters
    ----------
    r : array
        radial positions where the model is to be computed
    rt : float
        radius at which the maximum velocity is reached
    vmax : float
        Maximum velocity of the model

    Returns
    -------
    Array with the same shape of r, containing the model velocity curve/map

    """

    vr = np.zeros(np.shape(r))
    vr[(r <= rt)] = vmax*r[(r <= rt)]/rt
    vr[(r > rt)] = vmax

    return vr


def arctan_velocity(r, rd, vmax):
    """
    Velocity function for an arctan profile
	r and rd must be in the same units

    Parameters
    ----------
    r : float, array
        radial positions where the model is to be computed
    rt : float
        radius at which the maximum velocity is reached
    vmax : float
        Maximum velocity of the model

    Returns
    -------
    Array with the same shape of r, containing the model velocity curve/map


    """

    return 2*vmax/np.pi*np.arctan(2*r/rd)




  ##--- compute the 2D velocity model with the specific desired geometry ---##

def make_velmap(geom, models, instr):
    """
    Creates a model velocity map using the parameters in the configuration file
    via the objects *geom*, *models*, *instr*.

    For a detailed description of the configuration file and the associated
    objects, see the file read_params.py.
    For a list of the available models and the required parameters see the file
    flux_model.py.

    """

    # Compute the maps of the real radius and the cosine of azymutal angle theta,
    # corresponding to every pixel in the plane of the disc
    r,ctheta = tools.sky_coord_to_galactic(geom.x0, geom.y0, geom.pa, geom.incl,
	           im_size=(instr.size, instr.size))
    # Map of the circular velocity at any given radius
    if models.vm == "exp":
        vr = exponential_velocity(r, models.rdv, models.vmax)
    elif models.vm == "flat":
        vr = flat_velocity(r, models.rdv, models.vmax)
    elif models.vm == "arctan":
        vr = arctan_velocity(r, models.rdv, models.vmax)

    # Get the line-of-sight component of the circular velocity at each pixel
    v = vr * np.sin(np.radians(geom.incl)) * ctheta + geom.vsys

    return v


