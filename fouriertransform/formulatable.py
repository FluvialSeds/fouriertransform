'''
This module contains the FormulaTable class.
'''

from __future__ import(
	division,
	print_function,
	)

__docformat__ = 'restructuredtext en'
__all__ = ['FormulaTable']

import numpy as np
import pandas as pd
import warnings


#import exceptions
from .exceptions import(
	)

#import helper functions
from .formulatable_helper import(
	)

class FormulaTable(objet):
	__doc__ = '''
	Class for generating formula tables for FT-ICR MS sample sets, generating
	summary tables and statistics, and correlating formula abundances with
	environmental measurements.

	Parameters
	----------

	Warnings
	--------

	Notes
	-----

	See Also
	--------

	Examples
	--------

	**Attributes**

	'''

	def __init__(self):


	#define classmethod to generate instance from EnviroOrg output files
	@classmethod
	def from_eo(
		names,
		dir_path,
		rescale = None):
		'''
		Classmethod to generate a ``FormulaTable`` instance by inputting 
		EnviroOrg files generated for each sample from a given sample set.

		Parameters
		----------
		names : str or list
			Either a list of strings containing the filenames to be imported or
			the string 'all'. If 'all', method will automatically import all 
			files within the provided directory.

		dir_path : str
			String containing the (absolute) path pointing to the directory
			containing files to be imported.

		rescale : str or None
			Tells the method if and how to rescale peak intensities for cross-
			sample comparisons. Either `None`, 'fraction', or 'max_peak'. If
			'fraction', scales each sample such that the sum of intensities 
			for a given sample is equal to unity. If 'max_peak', scales such
			that the most intense peak detected in any sample is equal to 100
			and all other peaks are scaled accordingly. Defaults to `None`.

		Warnings
		--------

		Raises
		------

		Notes
		-----

		See Also
		--------

		References
		----------
		'''

	#define method for generating a summary table
	def generate_summary():
		'''
		Method for generating a summary table of the formula table, including
		relative compound class and category abundances (both by formula
		number and by MS intensity), total formulas detected, and average m/z
		(separated into classes and categories).

		Parameters
		----------

		Returns
		-------

		Warnings
		--------

		Raises
		------

		Notes
		-----

		See Also
		--------

		References
		----------
		'''


	#define method for plotting a sample van Krevelen
	def sample_vankrevelen(name, ax = None, type = 'class'):
		'''
		Method for generating a van Krevelen plot for a given sample within
		the sample set. Can color-code points either according to compound
		class (i.e. CHO, CHON, CHOS) or according to peak intensity.

		Parameters
		----------

		Returns
		-------

		Warnings
		--------

		Raises
		------

		Notes
		-----

		See Also
		--------

		References
		----------
		'''

	#define method for plotting a sample van Krevelen
	def correlation_vankrevelen(
		env_param, 
		ax = None, 
		fraction_present = 1.0,
		corr_coeff = 'Spearman'):
		'''
		Method for generating a van Krevelen plot to compare individual
		formula relative abundances against environmental parameters of
		interest. Can handle peaks that are only present in a subset of
		samples and can correlate using a range of statistical methods.

		Parameters
		----------

		Returns
		-------

		Warnings
		--------

		Raises
		------

		Notes
		-----

		See Also
		--------

		References
		----------
		'''








