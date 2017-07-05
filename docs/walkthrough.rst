Comprehensive Walkthrough
=========================
The following examples should form a comprehensive walkthrough of downloading the package, getting EnviroOrg output data into the right format for importing, generating a CrossTable instance, calculating and saving a summary table, and plotting various van Krevelen plots.

For detailed information on class attributes, methods, and parameters, consult the `Package Reference Documentation` or use the ``help()`` command within Python.

Quick Guide
-----------

Basic runthrough::

	#import modules
	import matplotlib.pyplot as plt
	import numpy as np
	import pandas as pd
	import fouriertransform as ft

	#generate string to path containing EnviroOrg data
	dir_name = 'path_name'

	#make the CrossTable instance
	ct = ft.CrossTable.from_eo(
		dir_name,
		file_names = 'all', #can replace with list of names
		rescale = 'fraction')

	#make the summary table
	sum_tab = ct.generate_summary()

	#save summary to csv
	sum_tab.to_csv('~/Desktop/summary_table.csv')

	#plot formula intensities and classes for a single sample
	fig, ax = plt.subplots(2,1)

	#make van krevelen plot of intensity
	ax[0] = ct.plot_sample_vk(
		'sample_name',
		ax = ax[0],
		plot_type = 'intensity',
		log = True,
		edgecolor = 'w',
		s = 40,
		cmap = 'YlGnBu')

	#make van krevelen plot of classes
	ax[1] = ct.plot_sample_vk(
		'sample_name',
		ax = ax[1]
		plot_type = 'class',
		edgecolor = 'w',
		s = 40)

	#plot the presence-absense difference between two samples
	fig, ax = plt.subplots(1,1)

	ax = ct.plot_difference_vk(
		'sample_name1',
		'sample_name2',
		ax = ax,
		edgecolor = 'w',
		s = 40)

	#plot the Spearman correlations with some environmental parameter
	fig, ax = plt.subplots(1,1)

	#make van krevelen plot
	ax = ct.plot_correlation_vk(
		env_param, #values of other parameter (e.g. d13C, D14C, etc.)
		ax = ax,
		corr_type = 'Spearman',
		f = 1, #fraction of samples that formula must be in
		alpha = 0.05, #significance cutoff
		edgecolor = 'w',
		s = 40,
		cmap = 'coolwarm',
		vmin = -1,
		vmax = 1)

Downloading the package
-----------------------

Using the ``pip`` package manager
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``fouriertransform`` and the associated dependencies can be downloaded directly from the command line using ``pip``::

	$ pip install fouriertransform

You can check that your installed version is up to date with the latest release by doing::

	$ pip freeze


Downloading from source
~~~~~~~~~~~~~~~~~~~~~~~
Alternatively, ``fouriertransform`` source code can be downloaded directly from `my github repo <http://github.com/FluvialSeds/fouriertransform>`_. Or, if you have git installed::

	$ git clone git://github.com/FluvialSeds/fouriertransform.git

And keep up-to-date with the latest version by doing::

	$ git pull

from within the fouriertransform directory.


Dependencies
~~~~~~~~~~~~
The following packages are required to run ``fouriertransform``:

* `python <http://www.python.org>`_ >= 2.7, including Python 3.x

* `matplotlib <http://matplotlib.org>`_ >= 1.5.2

* `numpy <http://www.numpy.org>`_ >= 1.11.1

* `pandas <http://pandas.pydata.org>`_ >= 0.18.1

* `scipy <http://www.scipy.org>`_ >= 0.18.0

If downloading using ``pip``, these dependencies (except python) are installed
automatically.

Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~
The following packages are not required but are highly recommended:

* `ipython <http://www.ipython.org>`_ >= 4.1.1

Additionally, if you are new to the Python environment or programming using the command line, consider using a Python integrated development environment (IDE) such as:

* `wingware <http://wingware.com>`_

* `Enthought Canopy <https://store.enthought.com/downloads/#default>`_

* `Anaconda <https://www.continuum.io/downloads>`_

* `Spyder <https://github.com/spyder-ide/spyder>`_

Python IDEs provide a "MATLAB-like" environment as well as package management. This option should look familiar for users coming from a MATLAB or RStudio background.

Detailed Walkthrough
--------------------
Coming soon!



