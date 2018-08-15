import sys
sys.path.append("../")
import flux_model as fm
import velocity_model as vm
import read_params as rp
from pylab import *
import operations as op

# in this example we use z=0.676 in order to shift the Ha rest frame wavelength
# to the center of the spectral coverage of the instrument in the nir mode

# Read in the needed parameters
geom, models, instr, envi = rp.get_allinput('params.config')

# Map of flux, normalized according to the effect of the spatial resolution
psf_mask = op.psf_footprint(geom, instr)
flux = fm.make_fluxmap(geom, models, instr)
flux = flux * psf_mask

# Maps of velocity, central wavelength of the corresponding emission, and dispersion
velmap = vm.make_velmap(geom, models, instr)
lambda_obs, sigma = op.observed_lambda(geom, models, instr, velmap)

# Calculate the final integrated spectrum
spectrum	= op.make_spectrum(instr, flux, lambda_obs, sigma, chunk = 100)

# check the result
plot(instr.channels,spectrum,'-',color="Gray")
xlabel("wavelenght [A]")
ylabel("Flux [arbitrary units]")
xlim((10990,11010))
show()


