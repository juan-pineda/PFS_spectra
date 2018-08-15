import os
import numpy as np
import math
import configparser
from astropy.cosmology import Planck13 as cosmo


class GeometryObj():
	"""
	Group the main parameters related to the geometrical orientation
	adopted for the mock observations. They include:

	incl : degrees
		inclination
	pa : degreees
		position angle
	x0 : pixels?
		position of the center
	y0 :  pixels ?
		position of the center
	vsys : km s^-1
		Systemic velocity of the galaxy
	z : redshift 

	"""

	def __init__(self):
		pass

	def parse_input(self, ConfigFile):
		"""
		Parse the geom. parameters from the configuration file to the object

		Parameters
		----------
		ConfigFile : string
			path to the configuration file (filename)

		Self consistent calculations of the cosmological luminosity distance}
		and the angular diameter distance are also performed and added as
		attributes of the object

		"""

		run_config = configparser.SafeConfigParser({}, allow_no_value=True)
		run_config.read(ConfigFile)
		self.incl = run_config.getfloat('geometry','incl')
		self.pa = run_config.getfloat('geometry','pa')
		self.x0 = run_config.getfloat('geometry','x0')
		self.y0 = run_config.getfloat('geometry','y0')
		self.vsys = run_config.getfloat('geometry','Vsys')
		self.z = run_config.getfloat('geometry','z')
		self.dl = cosmo.luminosity_distance(self.z).to('Mpc').value
		# note that the following could have been computed from the angular
		# diameter distance as well, it is equivalent
		self.dtheta = cosmo.kpc_proper_per_arcmin(self.z).to('kpc / arcsec').value


class ModelsObj():
	"""
	Group the main parameters describing the models selected to represent the
	2D flux and velocity maps. They include:

	fm : string
		short name of the chosen flux model
	i0 : float
		central brightness
	rdf : float
		characteristic scale-length of the model light profile
	rt :  float
		truncation radius
	vm : string
		short name of the chosen velocity model
	rdv :  float
		characteristic scale-length of the model velocity profile
	l0 : float
		central wavelngth of the emission line at rest frame
	disp : float
		velocity dispersion in km s^-1

	"""

	def __init__(self):
		pass

	def parse_input(self, ConfigFile):
        """
        Parse the geom. parameters from the configuration file to the object

        Parameters
        ----------
        ConfigFile : string
            path to the configuration file (filename)

        """

		run_config = configparser.SafeConfigParser({}, allow_no_value=True)
		run_config.read(ConfigFile)
		self.fm		= run_config.get('models','fm')
		self.i0   	= run_config.getfloat('models','I0')
		self.rdf	= run_config.getfloat('models','rdf')
		self.rt     = run_config.getfloat('models','rt')
		self.vm		= run_config.get('models','vm')
		self.vmax   = run_config.getfloat('models','vmax')
		self.rdv	= run_config.getfloat('models','rdv')
		self.l0		= run_config.getfloat('models','lambda0')
		self.disp	= run_config.getfloat('models','v_disp')


class InstrumentObj():
    """
    Group the main parameters describing the instrumental set up to be mocked
    They include:

    mode : string
		operation mode of the PFS instrument (e.g, blue, red, nir). Determines
		the wavelength range and number of channels
	psf : float
		spatial resolution at FWHM given in arcsec
	pixsize : float
		spatial sampling of the models in kpc
	fiber : float
		diameter of the fiber in arcsec

	"""

	def __init__(self):
		pass

	def parse_input(self, ConfigFile):
        """
        Parse the instr. parameters from the configuration file to the object
        
        Parameters
        ----------
        ConfigFile : string
            path to the configuration file (filename)
        
        """

		run_config = configparser.SafeConfigParser({}, allow_no_value=True)
		run_config.read(ConfigFile)
		self.mode	= run_config.get('instrument','mode')
		self.psf	= run_config.getfloat('instrument','psf')
		self.pixsize= run_config.getfloat('instrument','pixsize')
		self.fiber	= run_config.getfloat('instrument','fiber')

		if self.mode == 'blue':
			self.channels = np.arange(3800,6500.7,0.7)
			self.R	= 2300
		elif self.mode == 'red_low':
			self.channels = np.arange(6300,9700.9,0.9)
			self.R  = 3000
		elif self.mode == 'red_mid':
			self.channels = np.arange(7100,8850.4,0.4)
			self.R  = 5000
		elif self.mode == 'nir':
			self.channels = np.arange(9400,12600.8,0.8)
			self.R  = 4300
		# for running with much lower number of channelsm centered
		# by hand at l0*(1+z) for a specific configuration
		elif self.mode == 'test_z_0.5':
			self.channels = (np.arange(0,64,1)-30.5)*0.525024 + 9844.2
			self.R  = 4200
		# for running with much lower number of channels, centered
		# by hand at l0*(1+z) for a specific configuration
		elif self.mode == 'test_z_1.4':
			self.channels = (np.arange(0,64,1)-30.5)*0.8400384 + 15750.72
			self.R  = 4200



class EnvironmentObj():
	"""
	This class, to be developed, aims at managing information about the
	filesystem, like the path for i/o, etc.
	"""
	
	def __init__(self):
		pass
	
	def parse_input(self, ConfigFile):
		run_config = configparser.SafeConfigParser({}, allow_no_value=True)
		run_config.read(ConfigFile)
		self.path	= run_config.get('environment','path')
		self.filename	= run_config.get('environment','filename')


def get_allinput(ConfigFile):
	"""
	Read all attributes from the ConfigFile and adjust them for self-consistency
	Currently there are no default values for the missing parameters, so be
	careful handling the ConfigFile parameters

	Parameters
	----------
	ConfigFile : string
		location of the configuration file.

	"""

	# parse the input parameters from ConfigFile
	if(not os.path.isfile(ConfigFile)):
		print('// ' + ConfigFile + ' not found')
		return

	# create the main objects
	geom	= GeometryObj()
	geom.parse_input(ConfigFile)
	models	= ModelsObj()
	models.parse_input(ConfigFile)
	instr	= InstrumentObj()
	instr.parse_input(ConfigFile)
	envi	= EnvironmentObj()
	envi.parse_input(ConfigFile)
	# Change all scales to pixels units
	models.rdf	= models.rdf / instr.pixsize
	models.rdv  = models.rdv / instr.pixsize
	models.rt   = models.rt / instr.pixsize
    instr.fiber = instr.fiber * geom.dtheta / instr.pixsize
    instr.psf   = instr.psf * geom.dtheta / instr.pixsize
    geom.x0     = (geom.x0 * geom.dtheta / instr.pixsize)
    geom.y0     = (geom.y0 * geom.dtheta / instr.pixsize)
	# convert the psf from FWHM to 1-sigma std dev.
	instr.psf	= instr.psf / (2.0 * math.sqrt(2.0 * math.log(2.0)))
	# adapt the size of the models to cover at least 10 times psf-sigma in each
	# direction, beyond the size of the fiber
	instr.size	= int(math.ceil((instr.fiber + 20 * instr.psf)/2) * 2)
	# refer the misscentering with respect to the center of the maps
	geom.x0     = geom.x0 + (instr.size-1)/2
	geom.y0     = geom.y0 + (instr.size-1)/2

	return geom, models, instr, envi






