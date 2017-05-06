'''
This module contains the CrossTable class.
'''

from __future__ import(
	division,
	print_function,
	)

__docformat__ = 'restructuredtext en'
__all__ = ['CrossTable']

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings


#import exceptions
from .exceptions import(
	FormulaError,
	SampleError
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

	Santi-Temkiv et al. (2013), PLoS One, doi:10.1371/journal.pone.0053550.

	EnviroOrg_ software, National High Magnetic Field Laboratory,
	Tallahassee, FL.

	.. _EnviroOrg: https://nationalmaglab.org/user-facilities/icr/icr-software

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

	cmpd_mass : pd.Series
		Series of the mass of each formula, with index as formula composition.
		Length `nF`.

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

	def __init__(
		self, 
		intensities, 
		formulae = None,
		rescale = None, 
		sam_names = None):

		#check that intensities and formulae are in right format
		ints, forms, sams = _check_int(
			intensities, 
			formulae = formulae, 
			sam_names = sam_names)

		nF, nS = np.shape(ints)
		_check_forms(forms)

		#rescale intensities if necessary
		if rescale in ['fraction', 'Fraction']:
			m = ints.sum()
			ints = ints/m

		elif rescale in ['max_peak', 'Max_peak', 'Max_Peak']:
			m = ints.max().max()
			ints = ints*100/m

		elif rescale is not None:
			raise ValueError(
				'Rescale value of %r is not recognized. Must be "fraction",'
				' "max_peak", or "None".' % rescale)

		#store results
		self.intensities = ints
		self.sam_names = sams
		self.formulae = forms
		self.nF = nF
		self.nS = nS

		#generate a chemical composition table
		self._chem_comp = _gen_chem_comp(forms)

		#calculate compound mass, category, class, and AImod
		self.cmpd_mass = _calc_mass(self)
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

		Raises
		------
		FileError
			If the directory path does not exist.

		FileError
			If the inputted file does not exist or is not in .csv format.

		ValueError
			If the value for `rescale` is not 'fraction', 'max_peak', or 
			`None`.

		Notes
		-----
		This function requires that inputted files are in .csv format. File
		names will be used as sample names in resulting cross table (omitting
		the .csv extension), so using shortened filenames is encouraged for
		brevity.

		References
		----------
		EnviroOrg_ software, National High Magnetic Field Laboratory,
		Tallahassee, FL.

		.. _EnviroOrg: https://nationalmaglab.org/user-facilities/icr/icr-software
		'''

		#combine all files into a single dataframe
		intensities = _combine_EO_samples(dir_path, file_names)

		#return CrossTable instance
		return cls(
			intensities, 
			formulae = None, 
			rescale = rescale,
			sam_names = None)

	#define method for generating a summary table
	def generate_summary(self):
		'''
		Method for generating a summary table of the ``CrossTable`` instance, 
		including relative compound class and category abundances (both by 
		formula number and by MS intensity), total formulae detected, and 
		average m/z (both by formula number and by MS intensity).

		Returns
		-------
		sum_df : pd.DataFrame
			Dataframe of the summary statistics, including total formulae,
			average mass, percent compound class, and percent compound
			category (all but total formulae calculate both by peak number
			and by relative intensity) for each sample in the
			``CrossTable`` instance.

		References
		----------
		Santi-Temkiv et al. (2013), PLoS One, doi:10.1371/journal.pone.0053550.
		'''

		sum_df = _gen_sum_tab(self)

		return sum_df

	#define method for plotting a sample van Krevelen
	def plot_sample_vk(
		self,
		sam_name, 
		ax = None, 
		plot_type = 'class', 
		log = True,
		**kwargs):
		'''
		Method for generating a van Krevelen plot for a given sample within
		the sample set. Can color-code points either according to compound
		class (i.e. CHO, CHON, CHOS, CHOP, CHONS) or according to peak 
		intensity.

		Parameters
		----------
		sam_name : str
			String of the sample name to be plotted. Must be contained within
			the ``CrossTable`` instance.

		ax : None or plt.axis
			Axis to plot on. If `None`, creates an axis.

		plot_type : str
			Type of plot to make, either 'class' for plotting color-coded
			compound classes or 'intensity' for plotting peak intensities.

		log : boolean
			For intensity plots, tells the method whether or not to perform a
			log transform on intensities. If there are many low-abundance
			compounds, log transform will better show their presence.

		Returns
		-------
		ax : plt.axis
			Axis containing the van Krevelen plot of interest.

		Raises
		------
		SampleError
			If sam_name is not present within the ``CrossTable`` instance.

		ValueError
			If `type` is not 'class' or 'intensity'.

		See Also
		--------
		plot_difference_vk
			Function for plotting presence/absense differences between two
			samples.

		plot_correlation_vk
			Function for plotting correlation coefficients for individual
			chemical formulae relative to a given environmental controlling
			parameter.
		'''

		#check for errors
		if sam_name not in self.sam_names:
			raise SampleError(
				'Sample name %r is not contained within the CrossTable!' 
				% sam_name)

		if plot_type not in ['class','Class','intensity','Intensity']:
			raise ValueError(
				'Plot type %r not recognized, must be "class" or "intensity"'
				% plot_type)

		#make axis if necessary
		if ax is None:
			fig, ax = plt.subplots(1,1)

		#set axis labels
		ax.set_ylabel('H/C')
		ax.set_xlabel('O/C')

		#calculate H/C and O/C for formulae contained within the sample
		ints = self.intensities[self.intensities[sam_name] > 0][sam_name]
		ind = ints.index

		HC = self._chem_comp.loc[ind,'H']/self._chem_comp.loc[ind,'C']
		OC = self._chem_comp.loc[ind,'O']/self._chem_comp.loc[ind,'C']

		#if plot_type is 'intensities', log if necessary and plot
		if plot_type in ['intensity','Intensity']:
			#log transform if necessary
			if log:
				c = np.log10(ints)
				lab = r'$log_{10}$ peak intensity'
			else:
				c = ints
				lab = 'peak intensity'

			#sort by ascending intensity
			ind_sort = np.argsort(c)

			#plot results
			vk = ax.scatter(
				OC[ind_sort], 
				HC[ind_sort], 
				c = c[ind_sort].values, 
				**kwargs)

			#add colorbar
			cbar = fig.colorbar(vk, label = lab)

		#if plot type is 'class', extract classes and plot
		else:
			#extract all unique classes
			ccls = self.cmpd_class.unique()
			classes = self.cmpd_class[ind]

			#make a list of colors, extending to all possible classes
			colors = [
				[.071, 0, .40], #CHO
				[.431, .082, .035], #CHON
				[.769, .455, 0], #CHOS
				[0, .561, .569], #CHOP
				[.133, .376, .035], #CHONS
				]

			#loop through and plot
			for i, cc in enumerate(ccls):
				#extract indices
				ci = classes[classes == cc].index

				#plot
				ax.scatter(
					OC[ci],
					HC[ci],
					c = colors[i],
					label = cc,
					**kwargs)

			#add legend
			ax.legend(
				loc = 'best',
				framealpha = 1,
				scatterpoints = 1,
				markerscale = 2,
				fontsize = 'small')

		return ax

	#define method for plotting a sample van Krevelen
	def plot_difference_vk(
		self,
		sam_name1, 
		sam_name2, 
		ax = None,
		**kwargs):
		'''
		Method for generating a van Krevelen plot of the difference between
		two samples within a given sample set. Can color-code points either 
		according to compound class (i.e. CHO, CHON, CHOS) or according to 
		peak intensity.

		Parameters
		----------
		sam_name1 : str
			String of the name of the sample containing all formulae of 
			interest. (i.e. formulae present in this sample but not in 
			`sam_name2` will be plotted). Must be contained within the 
			``CrossTable`` instance.

		sam_name2 : str
			String of the name of the sample to be subtracted from
			`sam_name1`. Formulae absent from this sample but present in
			`sam_name1` will be plotted. Must be contained within the 
			``CrossTable`` instance.

		ax : None or plt.axis
			Axis to plot on. If `None`, creates an axis.

		Returns
		-------
		ax : plt.axis
			Axis containing the van Krevelen plot of interest.

		Raises
		------
		SampleError
			If sam_name is not present within the ``CrossTable`` instance.

		See Also
		--------
		plot_correlation_vk
			Function for plotting correlation coefficients for individual
			chemical formulae relative to a given environmental controlling
			parameter.

		plot_sample_vk
			Function for plotting a van Krevelen diagram for a single sample,
			color coded either by peak intensity or by compound class.
		'''

		#check that sample names exist
		if sam_name1 not in self.sam_names:
			raise SampleError(
				'Sample name %r is not contained within the CrossTable!' 
				% sam_name1)

		elif sam_name2 not in self.sam_names:
			raise SampleError(
				'Sample name %r is not contained within the CrossTable!' 
				% sam_name2)


		#make axis if necessary
		if ax is None:
			fig, ax = plt.subplots(1,1)

		#set axis labels
		ax.set_ylabel('H/C')
		ax.set_xlabel('O/C')

		#calculate index in sam_name1 but missing from sam_name2
		ind = self.intensities[
			(self.intensities[sam_name1] > 0) & 
			(self.intensities[sam_name2] == 0)
			].index

		#calculate H/C and O/C for formulae contained within the sample
		HC = self._chem_comp.loc[ind,'H']/self._chem_comp.loc[ind,'C']
		OC = self._chem_comp.loc[ind,'O']/self._chem_comp.loc[ind,'C']

		#extract all unique classes
		ccls = self.cmpd_class.unique()
		classes = self.cmpd_class[ind]

		#make a list of colors, extending to all possible classes
		colors = [
			[.071, 0, .40], #CHO
			[.431, .082, .035], #CHON
			[.769, .455, 0], #CHOS
			[0, .561, .569], #CHOP
			[.133, .376, .035], #CHONS
			]

		#loop through and plot
		for i, cc in enumerate(ccls):
			#extract indices
			ci = classes[classes == cc].index

			#plot
			ax.scatter(
				OC[ci],
				HC[ci],
				c = colors[i],
				label = cc,
				**kwargs)

		#add legend
		ax.legend(
			loc = 'best',
			framealpha = 1,
			scatterpoints = 1,
			markerscale = 2,
			fontsize = 'small')

		return ax

	#define method for plotting a sample van Krevelen
	def plot_correlation_vk(
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








