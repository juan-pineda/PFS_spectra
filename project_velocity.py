import math
import numpy as np
from pylab import *
from scipy import interpolate

def sky_coord_to_galactic(xcen, ycen, pos_angl, incl, im_size=(240, 240)):
    """
    Convert position from Sky coordinates to Galactic coordinates

    :param float xcen: position of the center in arcsec
    :param float ycen: position of the center in arcsec
    :param float pos_angl: position angle of the major axis degree
    :param float incl: inclination of the disk in degree
    :param ndarray im_size: maximum radius of the scene (arcsec),
                          im_size should be larger than the slit length + seeing (Default im_size=100)
    :param float res: resolution of the high resolution data (arcsec),
                      res should be at least n x pixel size (Default res=0.04)
    :return ndarray: [r, theta]
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

# RC is the rotation curve model. It has to have 2 columns. First column is radius,
# and second one the velocities. To avoid errors, the radial range must begin at 0m
# and go beyond the maximum radius you expect in your map
def make_velmap(RC, Vsys, xcen, ycen, pos_angl, incl, im_size=(240, 240)):
	r, ctheta= sky_coord_to_galactic(xcen, ycen, pos_angl, incl, im_size)
	spl = interpolate.interp1d(RC[:,0], RC[:,1])
	velmap = spl(r) * np.sin(np.radians(incl)) * ctheta + Vsys
	return velmap

#####--- EXAMPLE ---#####

# We will use a map of size (240, 240), and an inclination of 70 degrees
# This will put a maximum radius of approx. 481.28736. Look:
r, ctheta = sky_coord_to_galactic(120, 120, 30, 70, im_size=(240, 240))[:]
np.max(r)
#Out[61]: 481.2873645452771


# Now we define the rotation curve:
rad = np.arange(0,490,1)
Vcirc = np.arctan(rad/100)
# and concatenate it into a single 2-column array:
RC = np.concatenate([[rad],[Vcirc]]).T

# That's it, now we can create the model velocity map, with the desired projected
# geometry, adding an arbitrary systemic velocity, for instance Vsys=500

velmap = make_velmap(RC,0,120, 120, 30, 70, im_size=(240, 240))

# Take a look on the result
imshow(velmap)
colorbar()
contour(velmap,10,colors="white")
show()




