'''
This module contains the CrossTable class.
'''

from __future__ import(
	division,
	print_function,
	)

__docformat__ = 'restructuredtext en'
__all__ = ['CrossTable']

import numpy as np
import pandas as pd
import warnings


#import exceptions
from .exceptions import(
	FormulaError,
	)

#import helper functions
from .crosstable_helper import *


class CrossTable(object):
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

	References
	----------
	Koch and Dittmar (2006), Rapid Comm. Mass. Spec., 20, 926-932.


	See Also
	--------

	Examples
	--------

	**Attributes**

	AImod : pd.Series
		Series of modified AI value for each formula, with index as formula
		composition. Length `nF`.

	cmpd_cat : pd.Series
		Series of category of each formula (e.g. aliphatic, condensed
		aromatic high oxygen content, etc.), with index as formula 
		composition. Length `nF`.

	cmpd_class : pd.Series
		Series of class of each formula (i.e. CHO, CHON, CHOS, CHOP, CHONS), 
		with index as formula composition. Length `nF`.

	formulae : list
		List of strings containing each molecular formula, in the format:
		C(1-45)H(1-92)O(1-25)N(0-4)S(0-2). Length `nF`.

	nF : int
		Number of total detected formulae in dataset.

	nS : int
		Number of samples in dataset.

	sam_names : list
		List of strings containing each sample name. Length `nS`.

	intensities : pd.DataFrame
		DataFrame containing the MS intensities of each formula within each 
		sample. Shape [`nF` x `nS`].

	_chem_comp : pd.DataFrame
		Protected attribute. Dataframe of number of each atom contained in
		each formula (delete this text later, for my bookkeeping only.)

	'''

	def __init__(self, intensities, formulae = None, sam_names = None):

		#check that intensities and formulae are in right format
		ints, forms, sams = _check_int(
			intensities, 
			formulae = formulae, 
			sam_names = sam_names)

		nF, nS = np.shape(ints)
		_check_forms(forms)

		#store results
		self.intensities = ints
		self.sam_names = sams
		self.formulae = forms
		self.nF = nF
		self.nS = nS

		#generate a chemical composition table
		self._chem_comp = _gen_chem_comp(forms)

		#calculate compound category, class, and AImod
		self.AImod = _calc_AImod(self)
		self.cmpd_class = _calc_class(self)
		self.cmpd_cat = _calc_category(self)

	#define classmethod to generate instance from EnviroOrg output files
	@classmethod
	def from_eo(
		cls,
		dir_path,
		file_names = 'all',
		rescale = None):
		'''
		Classmethod to generate a ``CrossTable`` instance by inputting 
		EnviroOrg files generated for each sample from a given sample set.

		Parameters
		----------
		dir_path : str
			String containing the (absolute) path pointing to the directory
			containing files to be imported.

		file_names : str or list
			Either a list of strings containing the filenames to be imported or
			the string 'all'. If 'all', method will automatically import all 
			files within the provided directory. Defaults to 'all'.

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

		#combine all files into a single dataframe
		intensities = _combine_EO_samples(dir_path, file_names)

		#return CrossTable instance
		return cls(intensities, formulae = None, sam_names = None)



	#define method for generating a summary table
	def generate_summary():
		'''
		Method for generating a summary table of the formula table, including
		relative compound class and category abundances (both by formula
		number and by MS intensity), total formulae detected, and average m/z
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








