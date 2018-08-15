import sys
sys.path.append("../")
import read_params as rp
from pylab import *
import operations as op
from astropy.io import fits
import glob

# In this example, a fits file is provided containing a model intensity map,
# a model velocity map, and the resulting mock spectrum.
# The idea is to create the mock spectrum with PFS_spectra to compare the
# results. Then, as the maps are provided all the parameters usually passed
# to PFS_spectra to create the models are nuisance in this example.

# For this specific example the galaxy is placed at z=1.4, and the provided
# spectrum consists of 64 channels separated by 16 km s^-1 and it is centered
# in the center of the wavelength/velocity range.

# As those values do not correspond to any of the PFS modes, a new mode named
# 'test_z_1.4' was added to the read_params.py file. As PFS_spectra works in
# wavelength, not velocity, the equivalent channels were computed by hand:
#  
# center of the spectrum: 6562.8 * (1+z) = 15750.72 [A]
# width of the channels: delta_lambda = 6562.8 * (1+z) * delta_v / c = 0.8400384 [A]


# read in the provided fits file
x = glob.glob('*.fits')[0]

# this is the structure of the file
hdul = fits.open(x)
spec = hdul[0].data
flux = hdul[1].data
velmap = hdul[2].data


# Read in the needed parameters
geom, models, instr, envi = rp.get_allinput('params.config')

# get the size of the maps
instr.size = velmap.shape[0]

# Map of flux, normalized according to the effect of the spatial resolution
psf_mask    = op.psf_footprint(geom, instr)
flux        = flux * psf_mask

# Maps of velocity, central wavelength of the corresponding emission, and dispersion
lambda_obs, sigma	= op.observed_lambda(geom, models, instr, velmap)

# Calculate the final integrated spectrum
spectrum	= op.make_spectrum(instr, flux, lambda_obs, sigma, chunk = 100)

plot(spec/np.max(spec),"or",zorder=10, label= "provided spectrum")
plot(spectrum/np.max(spectrum),'-k',zorder=9, label="PFS_spectra")
xlabel("wavelenght [A]")
ylabel("Flux [arbitrary units]")
legend(frameon=False, loc="upper right")
show()


