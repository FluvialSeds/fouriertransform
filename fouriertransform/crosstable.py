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
	intensities : 2D array-like
		Either 2D np.ndarray or pd.Dataframe of all formula intensities for
		all samples.

	drop_int_above : none, int, or float
		Tells the method whether to drop formulae that represent greater
		than some amount of the total intensity -- i.e. fliers. The value
		of `drop_int_above` is the percent of total intensity (i.e. the
		relative abundance) value above which peaks will be dropped.
		Defaults to `None`.

	drop_high_OC : Boolean
		If `True`, formulae with O/C ratio above 1.0 will be dropped.
		Defaults to `True`.

	drop_high_HC : Boolean
		If `True`, formulae with H/C ratio above 2.0 will be dropped.
		Defaults to `True`.

	formulae : list
		List of formula names. If `None`, pulls formula names from 
		`intensities` index values. Defaults to `None`.

	rescale : `None` or str
		String telling the method how to rescale peak intensities, either
		'fraction', 'max_peak', or None. If 'fraction', scales intensities such
		that the sum for *each sample* is equal to 100. If 'max_peak', scales
		intensities such that the largest peak *in the entire sample set* is
		equal to 100.

	sam_names : `None` or list
		List of sample names. If `None`, pulls sample names from column names
		of `intensities`. Defaults to `None`.

	Warnings
	--------
	UserWarning
		If `sam_names` = `None` and `intensities` does not contain sample
		names, then names each sample according to column number.

	Notes
	-----
	When calculating correlation van Krevelen plots, `CrossTable` scaling
	should be set to 'fraction'. This way, it is the relative abundance within
	a given sample that is correlated with environmental parameters.

	References
	----------
	EnviroOrg_ software, National High Magnetic Field Laboratory,
	Tallahassee, FL.

	Koch and Dittmar (2006), Rapid Comm. Mass. Spec., 20, 926-932.

	LaRowe and van Cappellen (2011), Geochim. Cosmochim. Ac., 67, 2030-2042.

	Legendre and Legendre (2012), Numerical Ecology, Elsevier, 2nd Ed.

	Santi-Temkiv et al. (2013), PLoS One, doi:10.1371/journal.pone.0053550.

	.. _EnviroOrg: https://nationalmaglab.org/user-facilities/icr/icr-software

	Examples
	--------
	Generating an arbitrary bare-bones cross table containing made-up
	formulae::

		#import modules
		import fouriertransform as ft
		import numpy as np
		import pandas as pd

		#generate arbitrary data
		intensities = np.ones([2,2])
		formulae = ['C1H1O1N0S0P0', 'C2H2O1N0S1P0']
		sam_names = ['sample_1','sample_2']

		#create crosstable
		ct = ft.CrossTable(
			intensities,
			formulae = formulae,
			rescale = None,
			sam_names = sam_names)

	Importing all EnviroOrg files contained within a directory and combining
	to for a cross table::

		#input directory name
		dir_name = '~/Desktop/EO_directory'

		#create crosstable
		ct = ft.CrossTable.from_eo(
			dir_name,
			file_names = 'all', #can replace with list of names for subset
			rescale = 'fraction',
			drop_int_above = 5, #drops fliers above 5 pct total intensity
			drop_high_HC = True, #drops peaks with H/C > 2
			drop_high_OC = True) #drops peaks with O/C > 1

	Generating a summary table for the samples contained within the
	``CrossTable`` instance `ct`::

		#make the summary table
		sum_tab = ct.generate_summary()

		#save as csv
		sum_tab.to_csv('~/Desktop/summary_table.csv')

	Plotting the cross table data in van Krevelen space. You can plot either
	intensities or class for a single sample, the difference between two
	samples, or correlations between relative intensities across all samples
	and some other environmental parameter (e.g. bulk carbon isotopes). For
	all plotting methods, matplotlib scatterplot keyword arguments can be
	used.

	Plotting a single sample's intensity. Here, I'm taking the log of
	intensities to see a wider dynamic range::

		#import additional modules
		import matplotlib.pyplot as plt

		#extract sample name
		s = 'sample_1'

		#make axis
		fig, ax = plt.subplots(1,1)

		#make van krevelen plot
		ax = ct.plot_sample_vk(
			s,
			ax = ax,
			plot_type = 'intensity',
			log = True,
			edgecolor = 'w',
			s = 40,
			cmap = 'YlGnBu')

	Plotting compounds by class for a single sample::

		#make axis
		fig, ax = plt.subplots(1,1)

		#make van krevelen plot
		ax = ct.plot_sample_vk(
			s,
			ax = ax
			plot_type = 'class',
			edgecolor = 'w',
			s = 40)

	Plotting the difference between two samples. Formulae that are present
	in 'sample_1' but not present in 'sample_2' will be plotted and color
	coded according to class::

		#extract second sample name
		s2 = 'sample_2'

		#make axis
		fig, ax = plt.subplots(1,1)

		#make van krevelen plot
		ax = ct.plot_difference_vk(
			s,
			s2,
			ax = ax,
			edgecolor = 'w',
			s = 40)

	Plotting correlations with some independently measured environmental
	parameter, here called `env_param`. You can either plot correlations as
	Pearson or Spearman coefficients, and you can decide what fraction of
	the total samples a formula must be contained in in order to be considered
	for correlation. You can also determine the significance cutoff (alpha)
	and can only consider certain compound classes (i.e. CHO, CHON, CHOS)
	rather than all compounds. Resulting plot only retains statistically
	significant formulae. Resulting correlation coefficients and p-values for
	statistically significant formulae are also stored in 'stats'.
	An example::

		#make axis
		fig, ax = plt.subplots(1,2)

		#make van krevelen plot for all formulae
		ax[0], stats = ct.plot_correlation_vk(
			env_param,
			ax = ax[0],
			corr_type = 'Spearman',
			f = 1,
			alpha = 0.05,
			cmpd_cls = 'all',
			edgecolor = 'w',
			s = 40,
			cmap = 'coolwarm',
			vmin = -1,
			vmax = 1)	

		#make van krevelen plot for CHON formulae only
		ax[1], stats = ct.plot_correlation_vk(
			env_param,
			ax = ax[1],
			corr_type = 'Spearman',
			f = 1,
			alpha = 0.05,
			cmpd_cls = 'CHON',
			edgecolor = 'w',
			s = 40,
			cmap = 'coolwarm',
			vmin = -1,
			vmax = 1)	

	Any figure can then be saved to the drive according to::

		fig.savefig('~/Desktop/van_krevelen_figure.pdf')

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
		C(1-99)H(1-99)O(1-99)N(0-9)S(0-9)P(0-9). Length `nF`.

	intensities : pd.DataFrame
		DataFrame containing the MS intensities of each formula within each 
		sample. Shape [`nF` x `nS`].

	nF : int
		Number of total detected formulae in dataset.

	NOSC : pd.Series
		Series of the nominal oxidation state of carbon (NOSC) for each formula,
		with index as formula composition. Length `nS`.

	nS : int
		Number of samples in dataset.

	sam_names : list
		List of strings containing each sample name. Length `nS`.
	'''

	def __init__(
		self, 
		intensities,
		drop_int_above = None, #new v0.0.3
		drop_high_OC = True, #new v0.0.3
		drop_high_HC = True, #new v0.0.3
		formulae = None,
		rescale = None, 
		sam_names = None):

		#check that intensities and formulae are in right format
		ints_all, forms_all, sams = _check_int(
			intensities,
			formulae = formulae, 
			sam_names = sam_names)

		#new in v0.0.3: 
		#drop high H/C, high O/C, and fliers
		comps, forms, ints  = _drop_forms(
			forms_all,
			ints_all,
			drop_int_above = drop_int_above,
			drop_high_HC = drop_high_HC,
			drop_high_OC = drop_high_OC)

		nF, nS = np.shape(ints)
		_check_forms(forms)

		#rescale intensities if necessary
		if rescale in ['fraction', 'Fraction']:
			m = ints.sum()
			ints = 100*ints/m

		elif rescale in ['max_peak', 'Max_peak', 'Max_Peak']:
			m = ints.max()
			ints = ints*100/m

		elif rescale is not None:
			raise ValueError(
				'Rescale value of %r is not recognized. Must be "fraction",'
				' "max_peak", or "None".' % rescale)

		#generate a chemical composition table
		#comps = _gen_chem_comp(forms)

		#store results
		self.intensities = ints
		self.sam_names = sams
		self.formulae = forms
		self.nF = nF
		self.nS = nS
		self._chem_comp = comps

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
		drop_int_above = None,
		drop_high_OC = True,
		drop_high_HC = True,
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

		drop_int_above : none, int, or float
			Tells the method whether to drop formulae that represent greater
			than some amount of the total intensity -- i.e. fliers. The value
			of `drop_int_above` is the percent of total intensity (i.e. the
			relative abundance) value above which peaks will be dropped.
			Defaults to `None`.

		drop_high_OC : Boolean
			If `True`, formulae with O/C ratio above 1.0 will be dropped.
			Defaults to `True`.

		drop_high_HC : Boolean
			If `True`, formulae with H/C ratio above 2.0 will be dropped.
			Defaults to `True`.

		file_names : str or list
			Either a list of strings containing the filenames to be imported or
			the string 'all'. If 'all', method will automatically import all 
			files within the provided directory. Defaults to 'all'.

		rescale : str or None
			Tells the method if and how to rescale peak intensities for cross-
			sample comparisons. Either `None`, 'fraction', or 'max_peak'. If
			'fraction', scales each sample such that the sum of intensities 
			for a given sample is equal to unity. If 'max_peak', scales such
			that the most intense peak detected in each sample is equal to 100
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

		ValueError
			If `drop_int_above` is not `None`, float, or int.

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
		intensities = _combine_EO_samples(
			dir_path, 
			file_names)

		#return CrossTable instance
		return cls(
			intensities,
			drop_int_above = drop_int_above,
			drop_high_OC = drop_high_OC,
			drop_high_HC = drop_high_HC,
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
		plot_type = 'class',
		int_denom = 'max',
		shared = False,
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

		plot_type : str
			Type of plot to make, either 'class' for plotting color-coded
			compound classes or 'intensity' for plotting peak intensities.
			If 'intensity', plotted values are the difference in relative
			intensity between `sam_name1` and `sam_name2`, divided by a
			normalizing denominator (defined below):

					x = (r_1 - r_2) / denom

			i.e. `sam_name1` is the "final" sample and `sam_name2` is the
			"initila" sample.

		int_denom : str
			Normalizing denominator for 'intensity' difference plots,
			either 'max' or 'sum'. If 'max', divides by the maximum relative
			intensity in either `sam_name` or `sam_name2`. If 'sum', divides by
			the sum in both samples.

		shared : boolean
			For 'intensity' difference plots, if `True`, only considers
			formulae that are present in both samples. Defaults to `False`.

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

		#if plot type is 'class', extract differences and plot by class
		if plot_type in ['class', 'Class']:

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

		#if plot_type is 'intensities', calculate the difference in intensity
		# and plot
		elif plot_type in ['intensity','Intensity']:

			#calculate index in either sample

			#if shared, only formulae that are in both
			if shared:
				ind = self.intensities[
					(self.intensities[sam_name1] > 0) & 
					(self.intensities[sam_name2] > 0)
					].index				

			else:
				ind = self.intensities[
					(self.intensities[sam_name1] > 0) | 
					(self.intensities[sam_name2] > 0)
					].index

			#calculate H/C and O/C for formulae contained within the sample
			HC = self._chem_comp.loc[ind,'H']/self._chem_comp.loc[ind,'C']
			OC = self._chem_comp.loc[ind,'O']/self._chem_comp.loc[ind,'C']

			#calculate the relative intensities within each sample
			ints = self.intensities.loc[ind, [sam_name1, sam_name2]]

			if int_denom in ['sum','Sum']:
				#calculate the normalized difference
				c = (ints[sam_name1] - ints[sam_name2]) / ints.sum(axis = 1)
				lab = 'peak intensity difference (sum normalized)'

			elif int_denom in ['max','Max']:
				#calculate the normalized difference
				c = (ints[sam_name1] - ints[sam_name2]) / ints.max(axis = 1)
				lab = 'peak intensity difference (max normalized)'

			#raise error if denom type is not max or sum
			else:
				raise ValueError(
					'int_denom %r not recognized, must be "max" or "sum"'
					% int_denom)

			#sort by ascending intensity
			ind_sort = np.argsort(np.abs(c))

			#plot results
			vk = ax.scatter(
				OC[ind_sort], 
				HC[ind_sort], 
				c = c[ind_sort].values, 
				**kwargs)

			#add colorbar
			cbar = fig.colorbar(vk, label = lab)

		#raise error if type is not class of intensity
		else:
			raise ValueError(
				'Plot type %r not recognized, must be "class" or "intensity"'
				% plot_type)

		return ax

	#define method for plotting a sample van Krevelen
	def plot_correlation_vk(
		self,
		env_param,
		alpha = 0.05,
		ax = None,
		cmpd_cls = 'all',
		corr_type = 'Pearson',
		f = 1.0,
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

		alpha : float
			The significance value to use for retaining statistically
			significant formulae for plotting, must be between 0 and 1. 
			Defaults to `0.05`.

		ax : None or plt.axis
			Axis to plot on. If `None`, creates an axis. Defaults to `None`.

		cmpd_cls: str
			String saying which compound class to use for correlations.
			Must be one of:

				all, \n
				CHO, \n
				CHON, \n
				CHOS \n
			
			Defaults to 'all'.

		corr_type : str
			String saying which statistical method to use for correlation.
			Currently accepts:

				Pearson, \n
				Spearman \n

			Defaults to 'Pearson'.

		f : float
			The fraction of total samples in which a formula must be present
			in order to be considered for correlation, ranging between 0.0
			and 1.0. If some samples are missing env_param data, then f is the
			fraction of retained samples in which a formula must be present.
			Defaults to 1.

		Returns
		-------
		ax : plt.axis
			Axis containing the van Krevelen plot of interest.

		stats : df.DataFrame
			Dataframe of resulting r (for Pearson) or rho (for Spearman)
			values (collectively termed 'corr') and corresponding p-values for
			all formulae that are included in the plot.

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

		ValueError
			If `cmpd_cls` is not 'all', 'CHO', 'CHON', or 'CHOS'.

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

		#ensure cmpd_cls is the right string
		if cmpd_cls not in ['all','CHO','CHON','CHOS']:
			raise ValueError(
				'cmpd_cls value of %r not recognized. Must be one of "all", '
				' "CHO", "CHON", or "CHOS"' % cmpd_cls)

		#make axis if necessary
		if ax is None:
			fig, ax = plt.subplots(1,1)

		#set axis labels and title
		title = 'Env. corr. (' + corr_type + \
			', f = ' + str(f) + \
			', alpha = ' + str(alpha) + \
			', cmpds = ' + cmpd_cls + ')'
		
		ax.set_title(title)
		ax.set_ylabel('H/C')
		ax.set_xlabel('O/C')

		#retain only intensities from samples that exist in env_param
		# i.e. drop missing samples
		sams = env_param.index
		ints = self.intensities[sams]
		nS = len(sams)

		#extract indices of peaks present in >= f fraction of samples
		i = ints[(ints > 0).sum(axis = 1) >= f * nS].index

		#new in v.0.0.4!
		#only retain compounds within the selected class
		classes = self.cmpd_class[i]

		if cmpd_cls is not 'all':
			ind = classes[classes == cmpd_cls].index
		else:
			ind = i

		#store number of retained formulae
		nRet = len(ind)

		#calculate correlations and only keep significantly correlated forms
		ind_sig, rhos, pvals = _calc_corr(
			ints.loc[ind,:], 
			env_param, 
			alpha = alpha,
			corr_type = corr_type)

		nSig = len(ind_sig)

		#retain significant formulae and sort by ascending rho
		c = rhos[ind_sig]
		lab = corr_type + 'corr. coef.'

		#new in v.0.0.4!
		#plots highest ABSOLUTE rho values on top!
		ind_sort = np.abs(c).argsort()

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
		cbar = ax.figure.colorbar(vk, label = lab)

		#calculate summary statistics and plot as text
		mu = rhos[ind_sig].abs().mean()
		sig = rhos[ind_sig].abs().std()

		pct_sig = 100*nSig / nRet
		l1 = r'$n_{sig}$ = %.0f (%.1f %%); ' % (nSig, pct_sig)
		l2 = r'$\mu_{\| \rho \|}$ = %.2f; ' % mu
		l3 = r'$\sigma_{\| \rho \|}$ = %.2f' % sig
		text = l1 + l2 + l3

		#make text
		ax.text(0.05, 0.10, text,
			transform = ax.transAxes,
			horizontalalignment = 'left',
			verticalalignment = 'top')

		#set ranges
		ax.set_ylim([0,2])
		ax.set_xlim([0,1])

		#new in version 0.0.4!
		#make statistics dataframe
		stats = pd.DataFrame(
			index = ind_sig,
			columns = ['corr','pvals'])

		stats['corr'] = rhos[ind_sig]
		stats['pvals'] = pvals[ind_sig]

		return ax, stats
