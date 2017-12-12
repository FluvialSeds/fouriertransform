'''
This module contains functions for filtering sample data
'''

from __future__ import(
	division,
	print_function,
	)

__docformat__ = 'restructuredtext en'
__all__ = [
	]

import numpy as np
import pandas as pd
import warnings

from scipy import stats

#import exceptions
from .exceptions import(
	# DimError,
	# FileError,
	# FormulaError,
	# LengthError,
	# SampleError,
	)

#define function to remove abundance outliers
def _rem_abund_outliers(df_in,
	Dmz = 10,
	niter = 10,
	thresh = 2.5):
	'''
	Function to remove outliers based on high relative abundance (e.g.
	contaminants, highly ionizable compounds, etc.). This function uses the
	derivative of a plot of m/z vs. rel. abund. and removes peaks that
	produce a spike in this function.

	Parameters
	----------
	df_in : pd.DataFrame
		Dataframe of input data, containing 'Theor_mz' and 'rel_abund' columns.

	Dmz : int
		Bin size, in integer m/z units. Defaults to `10`.

	niter : int
		Maximum number of times to iterate peak removal. Defaults to `10`.

	thresh : float or int
		Threshold for derivative changes. Peaks producing a change in the
		derivative function above this threshold will be removed. Defaults to
		`1`. (i.e. 1% rel. abund.)

	Returns
	-------

	Raises
	------

	'''

	#check input data format and shape

	#check thresh and niter types

	#rescale to ensure largest peak = 100

	#extract columns
	mz = df_in['Theor_mz']
	ra = df_in['Rel_abund']

	#calculate max and min
	mz_min = np.floor(mz.min())
	mz_max = np.ceil(mz.max())

	#make bins and calc. ave. m/z in each bin
	bins = np.arange(mz_min, mz_max + Dmz, Dmz)
	mz_av = bins[:-1] + Dmz/2

	#sort each peak into an integer bin and calculate the max for each bin
	x = pd.cut(mz, bins)
	ra_max = ra.groupby(x).max() #max value in each bin
	ra_max_ind = ra.groupby(x).idxmax() #index of max value in each bin

	#calculate derivative of ra w.r.t. mz and store in series
	dradmz = np.gradient(ra_max)/Dmz
	dradmz = pd.Series(dradmz, index = ra_max.index)

	#determine formula of largest peak within bins > thresh
	i = dradmz[dradmz > thresh].index
	fliers = ra_max_ind[i].values








