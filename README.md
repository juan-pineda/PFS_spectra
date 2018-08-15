# PFS_spectra
Python library for mocking spectral observations of galaxies

Originally conceived to mimic the PFS instrument at SUBARU,
**PFS_spectra** can be easily adapted to any other configuration.   


## Contributors:  
Juan Carlos Basto Pineda  
Jérémy Dumoulin


# How to use PFS_spectra

To have a quick first idea of the way to use PFS_spectra just open the file
`example.params.config` provided in the main distribution.

As you see there, this file is used to feed PFS_spectra all the necessary
parameters regarding the instrumental and observational set up to be mimicked.

The params.config file is divided in three/four sections.

1. `[geometry]` aggregates the parameters related to the geometry of the
mock observation.
2. `[models]` aggregates the parameters needed to define the models for the
intensity and velocity maps.
3. `[instrument]` aggregates the parameters relevant to define the exact
instrumental configuration to be mocked.
4. `[environment]` is projected (but not used yet) to store the variables
related to the filesystem/environment.

Once you put all the required parameters in the **ConfigFile**, the way to pass
them to PFS_spectra is with the functions stored in the file **read_params.py**.
The main function there is `get_allinput`. As you see, it receives as input the
ConfigFile, reads in the variables stored therein, and create three objects
which will have those values as their attributes. Each object correspond to
one of the blocks described before, as identified by the names of the classes
and of the objects themselves:

```
geom    = GeometryObj()
geom.parse_input(ConfigFile)
models  = ModelsObj()
models.parse_input(ConfigFile)
instr   = InstrumentObj()
instr.parse_input(ConfigFile)
```

For instance, across the code the geometrical parameters can be invoked at
any time through calls of the instances `geom.incl`, `geom.pa`, `geom.z`, etc. To get
the exact name of all the instances please refer to the `read_params.py` file.

OK, so far PFS_spectra has read the input parameters and created the objects
`geom`, `models`, and `instr` enclosing all the relevant information. Now, how to use
it to mock the observations?

The first step is to create the mock velocity and intensity maps. To do that,
you simply need to use the function `make_velmap(geom, models, instr)` from the
**velocity_models.py** module, and the function `make_fluxmap(geom, models, instr)`
from the **flux_model.py** module. Those functions make use of the routine
`sky_coord_to_galactic` in the module **tools.py**, to project the coordinates in the
plane of the disc to the coordinates in the plane of the sky.

With the maps at hand, the main routines for the next step are in the file
**operations.py**. This includes a function to parametrize the emission line coming
out of each pixel of the maps, shifted as necessary according to the cosmology
and to the proper line-of-sight velocities. A function to record the footprint
of the fiber convolved with the psf, which allows to weight the contribution 
of the flux of individual pixels
in the final spectrum, and the main function for the creation of mock spectra:
`make_spectrum(instr, flux, lambda_obs, sigma, chunk = 50)`.

To understand the order in which these functions should be concatenated in a
clean and efficient way, please check the example **test_0.py**, provided in the
repository. A second example is given in **test_1.py**. After that you are invited
to navigate through the code to get to know better the details.


# Set up

As a final remark, just a friendly reminder for python beginners:  

In order for
**PFS_spectra** to work, you need to make its routines accessible to python from
wherever you are working. This can be achieved en either of 3 ways:

1. by copying all necessary files in the working directory
2. by adding the location of the repository to the `sys.path` variable at the
   beginning of your script with:

```
import sys
sys.path.append("path_to_PFS_files")
```

3. By adding the location of the repository permanently to your `PYTHONPATH`
   variable, adding at your **.bashrc** file something like:

```
export PYTHONPATH="${PYTHONPATH}:/home/juan/Desktop/PFS_spectra/"
```

Then you can import the necessary modules, doing for example:

```
import read_params as rp
import operations as op
```




