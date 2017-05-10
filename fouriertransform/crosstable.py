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
	LengthError,
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

	Legendre and Legendre (2012), Numerical Ecology, Elsevier, 2nd Ed.

	LaRowe and van Cappellen (2011), Geochim. Cosmochim. Ac., 67, 2030-2042.

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

	NOSC : pd.Series
		Series of the nominal oxidation state of carbon (NOSC) for each formula,
		with index as formula composition. Length `nS`.

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
			ints = 100*ints/m

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

		#calculate compound mass, category, class, NOSC, and AImod
		self.cmpd_mass = _calc_mass(self)
		self.AImod = _calc_AImod(self)
		self.NOSC = _calc_NOSC(self)
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

		#make axis if necessary
		if ax is None:
			fig, ax = plt.subplots(1,1)

		#set axis labels and title
		ax.set_title(sam_name + ' ' + plot_type)
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
		elif plot_type in ['class', 'Class']:
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

		#raise error if type is not class of intensity
		else:
			raise ValueError(
				'Plot type %r not recognized, must be "class" or "intensity"'
				% plot_type)

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

		#set axis labels and title
		ax.set_title(sam_name1 + ' - ' + sam_name2 + ' compound classes')
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
		self,
		env_param, 
		ax = None, 
		alpha = 0.05,
		f = 1.0,
		corr_type = 'Pearson',
		**kwargs):
		'''
		Method for generating a van Krevelen plot to compare individual
		formula relative abundances against environmental parameters of 
		interest. Can handle peaks that are only present in a subset of
		samples and can correlate using a range of statistical methods.

		Parameters
		----------
		env_param : list, array, or pd.Series
			List, array, or series of values for a particular environmental 
			variable. If list or array, assumes values are in the same order 
			as sample names in the ``CrossTable`` instance. Length `nS`. If
			Series with sample names as index, only considers samples present
			in the index for correlation (i.e. can handle missing samples).
			Samples with NaN values will also be dropped.

		ax : None or plt.axis
			Axis to plot on. If `None`, creates an axis. Defaults to `None`.

		alpha : float
			The significance value to use for retaining statistically
			significant formulae for plotting, must be between 0 and 1. 
			Defaults to `0.05`.

		f : float
			The fraction of total samples in which a formula must be present
			in order to be considered for correlation, ranging between 0.0
			and 1.0. If some samples are missing env_param data, then f is the
			fraction of retained samples in which a formula must be present.
			Defaults to 1.

		corr_type : str
			String saying which statistical method to use for correlation.
			Currently accepts:

				Pearson, \n
				Spearman \n

		Returns
		-------
		ax : plt.axis
			Axis containing the van Krevelen plot of interest.

		Raises
		------
		TypeError
			If `env_param` is not list, np.ndarray, or pd.Series.

		LengthError
			If `env_param` is not of length `nS` (if list or array only).

		ValueError
			If `fraction_present` is not between 0 and 1.

		LengthError
			If number of samples is less than 3.

		ValueError
			If `corr_type` is not 'Pearson' or 'Spearman'.

		TypeError
			If `intensities` is not ``pd.DataFrame``.

		TypeError
			If `env_param` is not ``pd.Series``.

		TypeError
			If `alpha` is not float between 0 and 1.

		SampleError
			If sample names in `env_param` and `intensities` do not match.


		See Also
		--------
		plot_difference_vk
			Function for plotting presence/absense differences between two
			samples.

		plot_sample_vk
			Function for plotting a van Krevelen diagram for a single sample,
			color coded either by peak intensity or by compound class.

		References
		----------
		Legendre and Legendre (2012), Numerical Ecology, Elsevier, 2nd Ed.
		
		'''

		#ensure env_param is the right type
		stype = type(env_param).__name__
		l = len(env_param)

		if stype == 'Series':
			#drop missing values if they exist
			env_param.dropna(inplace = True)

		elif stype in ['list', 'ndarray']:
			#ensure right length
			if l != self.nS:
				raise LengthError(
					'when list or array, `env_param` must contain the same'
					' number of samples as the ``CrossTable`` object.')

			env_param = pd.Series(env_param, index = self.sam_names)

		else:
			raise TypeError(
				'env_param is type %r. Must be list, np.ndarray, or pd.Series'
				% s)	

		#ensure fraction present between 0 and 1
		if type(f) not in [float, int] or f < 0 or f > 1:
			raise ValueError(
				'f must be int or float between 0 and 1!')

		#make axis if necessary
		if ax is None:
			fig, ax = plt.subplots(1,1)

		#set axis labels and title
		title = 'Environmental correlations (' + corr_type + \
			', f = ' + str(f) + \
			', alpha = ' + str(alpha) + ')'
		
		ax.set_title(title)
		ax.set_ylabel('H/C')
		ax.set_xlabel('O/C')

		#retain only intensities from samples that exist in env_param
		# i.e. drop missing samples
		sams = env_param.index
		ints = self.intensities[sams]
		nS = len(sams)

		#extract indices of peaks present in >= f fraction of samples
		ind = ints[(ints > 0).sum(axis = 1) >= f * nS].index

		#store number of retained formulae
		nRet = len(ind)

		#calculate correlations and only keep significantly correlated forms
		ind_sig, rhos, pvals = _calc_corr(
			ints.ix[ind], 
			env_param, 
			alpha = alpha,
			corr_type = corr_type)

		nSig = len(ind_sig)

		#retain significant formulae and sort by ascending rho
		c = rhos[ind_sig]
		lab = corr_type + 'corr. coef.'
		ind_sort = np.argsort(c)

		#calculate H/C and O/C
		HC = self._chem_comp.loc[ind_sig,'H']/self._chem_comp.loc[ind_sig,'C']
		OC = self._chem_comp.loc[ind_sig,'O']/self._chem_comp.loc[ind_sig,'C']

		#plot results
		vk = ax.scatter(
			OC[ind_sort], 
			HC[ind_sort], 
			c = c[ind_sort].values, 
			**kwargs)

		#add colorbar
		cbar = fig.colorbar(vk, label = lab)

		#calculate summary statistics and plot as text
		mu = rhos[ind_sig].abs().mean()
		sig = rhos[ind_sig].abs().std()

		pct_sig = 100*nSig / nRet
		l1 = r'$n_{sig}$ = %.0f (%.1f %% of total forms); ' % (nSig, pct_sig)
		l2 = r'$\mu_{\| \rho \|}$ = %.2f; ' % mu
		l3 = r'$\sigma_{\| \rho \|}$ = %.2f' % sig
		text = l1 + l2 + l3

		#make text
		ax.text(0.05, 0.95, text,
			transform = ax.transAxes,
			horizontalalignment = 'left',
			verticalalignment = 'top')

		return ax
