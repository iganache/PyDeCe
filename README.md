# PyDeCe
A python code for modeling the dense endmember of pyroclastic density currents (PDCs) generated either by impulsive column collapse or sustained fountaining eruptions. Dense, particle rich PDC is modeled as solid-fluid mixture driven by gravity analogous to the granular flow models of Iverson and Denlinger (2001). Flow movement over real topography is realized by using a digital elevation model (DEM) file as one of the model inputs. Other model inputs include simulation time, flow density and viscosity, x and y coordinates (or longitude and latitude) of the source, among others, which are input to the model either using a config file or via command line arguments.

See [usage instructions](https://github.com/iganache/PyDeCe/wiki) for more information on the model and execution.

Requires [NumPy](https://numpy.org/doc/stable/contents.html), [SciPy](https://www.scipy.org/), and [Rasterio](https://rasterio.readthedocs.io/en/latest/).
