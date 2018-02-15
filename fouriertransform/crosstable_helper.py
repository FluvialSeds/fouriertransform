'''
This module contains helper functions for the CrossTable class.
'''

from __future__ import(
	division,
	print_function,
	)

__docformat__ = 'restructuredtext en'
__all__ = [
	'_calc_AImod',
	'_calc_category',
	'_calc_class',
	'_calc_mass',
	'_calc_NOSC',
	'_calc_pct',
	'_calc_corr',
	'_check_forms',
	'_check_int',
	'_combine_EO_samples',
	'_drop_forms',
	'_gen_chem_comp',
	'_gen_sum_tab',
	'_input_EO_sample',
	]

import numpy as np
import os
import pandas as pd
import re
import warnings

from scipy import stats

#import exceptions
from .exceptions import(
	DimError,
	FileError,
	FormulaError,
	LengthError,
	SampleError,
	)

#define function to calculate modified aromaticity index
def _calc_AImod(ct):
	'''
	Calculates the modified armoaticity index for each formula contained
	within a ``CrossTable`` instance.

	Parameters
	----------
	ct : ft.CrossTable
		``CrossTable`` instance containing the formulae of interest.

	Returns
	-------
	AImod : pd.Series
		Series of the resulting AImod values. Length `nF`.

	References
	----------
	Koch and Dittmar (2006), Rapid Comm. Mass. Spec., 20, 926-932.
	'''

	#extract chemical compositions, name df for shorthand
	df = ct._chem_comp

	#UPDATED TO INCLUDE NITROGEN!
	#calculate DBEai
	#DBEai = 1 + df['C'] - 0.5*df['O'] - df['S'] - 0.5*df['H'] - 0.5*df['N']
	DBEai = 1 + df['C'] - df['S'] - 0.5*(df['O'] + df['H'] + df['N'])

	#calculate Cai
	Cai = df['C'] - 0.5*df['O'] - df['N'] - df['S'] - df['P']

	#calculate AImod and name it
	AImod = DBEai / Cai
	AImod.name = 'AImod'

	#make zero where negative
	AImod[AImod < 0] = 0

	return AImod

#define function to calculate formula compound category
def _calc_category(ct):
	'''
	Calculates the compound category of each formula contained within a
	``CrossTable`` instance as:

		Aliphatic, high oxygen content \n
		Aliphatic, low oxygen content \n
		Condensed Aromatic \n
		Polyphenolic, high oxygen content \n
		Polyphenolic, low oxygen content \n
		Unsaturated and phenolic, high oxygen content \n
		Unsaturated and phenolic, low oxygen content \n
		Peptide-like.

	Parameters
	----------
	ct : ft.CrossTable
		``CrossTable`` instance containing the formulae of interest.

	Returns
	-------
	ccat : pd.Series
		Series of the resulting categories. Length `nF`.

	References
	----------
	Santi-Temkiv et al. (2013), PLoS One, doi:10.1371/journal.pone.0053550.

	Notes
	-----
	This is the only location where compound categories are defined. All other
	functions are dependent on the defintions provided here. If categories
	need to be changed in the future (potentially user-defined categories
	option?), this is the only place to do so!
	'''

	#extract chemical compositions and AImod and combine
	ccs = ct._chem_comp
	ais = ct.AImod
	df = pd.concat([ccs, ais], axis = 1)

	#extract list of formulae stored in CrossTable instance
	ccat = pd.Series(index = df.index, name = 'cmpd_cat')

	#calculate indices for each category

	# #aliphatic, high OC
	# ah_ind = df[
	# 	(df['H'] >= 1.5*df['C']) & 
	# 	(df['O'] >= 0.5*df['C']) &
	# 	(df['N'] == 0)].index

	# #aliphatic, low OC
	# al_ind = df[
	# 	(df['H'] >= 1.5*df['C']) & 
	# 	(df['O'] < 0.5*df['C']) &
	# 	(df['N'] == 0)].index

	# #condensed aromatic
	# ca_ind = df[(
	# 	df['AImod'] >= 0.67)].index

	# #polyphenolic, high OC
	# ph_ind = df[
	# 	(df['AImod'] < 0.67) & 
	# 	(df['AImod'] >= 0.5) & 
	# 	(df['O'] >= 0.5*df['C'])].index

	# #polyphenolic, low OC
	# pl_ind = df[
	# 	(df['AImod'] < 0.67) & 
	# 	(df['AImod'] >= 0.5) & 
	# 	(df['O'] < 0.5*df['C'])].index

	# #unsaturatied and phenolic, high OC
	# uh_ind = df[
	# 	(df['AImod'] < 0.5) & 
	# 	(df['H'] < 1.5*df['C']) & 
	# 	(df['O'] >= 0.5*df['C'])].index

	# #unsaturatied and phenolic, low OC
	# ul_ind = df[
	# 	(df['AImod'] < 0.5) & 
	# 	(df['H'] < 1.5*df['C']) & 
	# 	(df['O'] < 0.5*df['C'])].index

	# #peptide-like
	# pe_ind = df[
	# 	(df['H'] >= 1.5*df['C']) & 
	# 	(df['N'] > 0)].index

	#aliphatic, high OC
	ah_ind = df[
		(df['H'] >= 1.5*df['C']) & 
		(df['O'] >= 0.5*df['C']) &
		(df['O'] <= 0.9*df['C']) &
		(df['N'] == 0)].index

	#aliphatic, low OC
	al_ind = df[
		(df['H'] >= 1.5*df['C']) & 
		(df['O'] < 0.5*df['C']) &
		(df['N'] == 0)].index

	#condensed aromatic
	ca_ind = df[
		(df['AImod'] > 0.66) &
		(df['N'] + df['S'] == 0) & 
		(df['O'] <= 0.9*df['C'])].index

	#condensed aromatic with S or N
	cx_ind = df[
		(df['AImod'] > 0.66) &
		(df['N'] + df['S'] > 0) & 
		(df['O'] <= 0.9*df['C'])].index

	#polyphenolic, high OC
	ph_ind = df[
		(df['AImod'] <= 0.66) & 
		(df['AImod'] > 0.5) & 
		(df['O'] >= 0.5*df['C']) &
		(df['O'] <= 0.9*df['C'])].index

	#polyphenolic, low OC
	pl_ind = df[
		(df['AImod'] <= 0.66) & 
		(df['AImod'] > 0.5) & 
		(df['O'] < 0.5*df['C'])].index

	#unsaturatied and phenolic, high OC
	uh_ind = df[
		(df['AImod'] <= 0.5) & 
		(df['H'] < 1.5*df['C']) & 
		(df['O'] >= 0.5*df['C']) &
		(df['O'] <= 0.9*df['C'])].index

	#unsaturatied and phenolic, low OC
	ul_ind = df[
		(df['AImod'] <= 0.5) & 
		(df['H'] < 1.5*df['C']) & 
		(df['O'] < 0.5*df['C'])].index

	#peptide-like
	pe_ind = df[
		(df['H'] >= 1.5*df['C']) & 
		(df['N'] > 0) &
		(df['O'] <= 0.9*df['C'])].index

	#sugars
	su_ind = df[
		(df['O'] > 0.9*df['C']) &
		(df['N'] + df['S'] == 0)].index

	#sugars with S or N
	sx_ind = df[
		(df['O'] > 0.9*df['C']) &
		(df['N'] + df['S'] > 0)].index

	#set values for compounds
	ccat[ah_ind] = 'aliphatic_highOC'
	ccat[al_ind] = 'aliphatic_lowOC'

	ccat[ca_ind] = 'condensed_aromatic'
	ccat[cx_ind] = 'condensed_aromatic_X'
	
	ccat[ph_ind] = 'polyphenolic_highOC'
	ccat[pl_ind] = 'polyphenolic_lowOC'

	ccat[uh_ind] = 'unsaturated_phenolic_highOC'	
	ccat[ul_ind] = 'unsaturated_phenolic_lowOC'

	ccat[pe_ind] = 'peptide_like'

	ccat[su_ind] = 'sugars'
	ccat[sx_ind] = 'sugars_X'

	return ccat

#define function to calculate formula class
def _calc_class(ct):
	'''
	Calculates the compound class of each formula contained within a
	``CrossTable`` instance as CHO, CHON, CHOS, CHOP, or CHONS.

	Parameters
	----------
	ct : ft.CrossTable
		``CrossTable`` instance containing the formulae of interest.

	Returns
	-------
	cclass : pd.Series
		Series of the resulting classes. Length `nF`.
	'''

	#extract chemical compositions, name df for shorthand
	df = ct._chem_comp

	#extract list of formulae stored in CrossTable instance
	cclass = pd.Series(index = df.index, name = 'cmpd_class')

	#calculate indices using chemical composition logic
	cho_ind = df[(df['N'] == 0) & (df['S'] == 0) & (df['P'] == 0)].index
	chon_ind = df[(df['N'] != 0) & (df['S'] == 0)].index
	chos_ind = df[(df['N'] == 0) & (df['S'] != 0)].index
	chop_ind = df[(df['N'] == 0) & (df['S'] == 0) & (df['P'] != 0)].index
	chons_ind = df[(df['N'] != 0) & (df['S'] != 0)].index

	#set values for compound class
	cclass[cho_ind] = 'CHO'
	cclass[chon_ind] = 'CHON'
	cclass[chos_ind] = 'CHOS'
	cclass[chop_ind] = 'CHOP'
	cclass[chons_ind] = 'CHONS'

	return cclass

#define a function to calculate pearson correlations and retain sig. forms.
def _calc_corr(intensities, env_param, alpha = 0.05, corr_type = 'Pearson'):
	'''
	Calculates the correlation coefficients and retains statistically
	significant compounds.

	Parameters
	----------
	intensities : pd.DataFrame
		Dataframe of peak intensities.

	env_param : pd.Series
		Series containing the values of a given environmental parameter to
		correlate with.

	alpha : float
		Significance alpha value for retaining formulae.

	corr_type : str
		The type or correlation to calculate, either 'Pearson' or 'Spearman'.

	Returns
	-------
	ind_sig : pd.Index
		Index of the statistically significant formulae.

	rhos : pd.Series
		Series of the correlation coefficient values

	pvals : pd.Series
		Series of the correlation significance values.

	Raises
	------
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
	'''

	#check formats and values
	itype = type(intensities).__name__
	etype = type(env_param).__name__

	if itype != 'DataFrame':
		raise TypeError(
			'intensities is currently type %r. Must be pd.DataFrame!' % itype)

	elif etype != 'Series':
		raise TypeError(
			'env_param is currently type %r. Must be pd.Series!' % itype)
	
	elif type(alpha) is not float or alpha > 1 or alpha < 0:
		raise TypeError(
			'alpha must be float between 0 and 1')

	elif (env_param.index != intensities.columns).any():
		raise SampleError(
			'sample names in env_param must match those in intensities.')

	#extract constants and raise error of <= 2 samples
	nS = len(env_param)

	if nS < 3:
		raise LengthError(
			'Must contain >= 3 samples to calculate correlations! Currently'
			' trying to correlate for %r samples' % nS)

	#rank transform if Spearman, rename to X and Y
	if corr_type in ['spearman', 'Spearman']:
		X = intensities.rank(axis = 1)
		Y = env_param.rank()

	elif corr_type in ['pearson', 'Pearson']:
		X = intensities
		Y = env_param

	else:
		raise ValueError(
			'Correlation of type %r not recognized. Must be "Pearson" or'
			' "Spearman"' % corr_type)

	#calculate rho
	num = np.mean(X*Y, axis = 1) - np.mean(X, axis = 1)*np.mean(Y)
	denom = (np.mean(X**2, axis = 1) - np.mean(X, axis = 1)**2)**0.5 * \
		(np.mean(Y**2) - np.mean(Y)**2)**0.5

	#store rho values
	rhos = num / denom

	#calculate t scores
	t = rhos / ((1 - rhos**2)/(nS - 2))**0.5
	p = 2*(1-stats.t.cdf(np.abs(t),nS-2))

	pvals = pd.Series(p, index = rhos.index)

	#calculate indices of statistically significant formulae
	ind_sig = pvals[pvals <= alpha].index

	return ind_sig, rhos, pvals

#define a function to calculate the mass of each compound
def _calc_mass(ct):
	'''
	Calculates the compound masses.

	Parameters
	----------
	ct : ft.CrossTable
		``CrossTable`` instance containing the formulae of interest.

	Returns
	-------
	cmass : pd.Series
		Series of the resulting masses. Length `nF`.
	'''

	#extract chemical compositions, name df for shorthand
	df = ct._chem_comp

	cmass = 12.01*df['C'] + \
		1.01*df['H'] + \
		15.99*df['O'] + \
		14.01*df['N'] + \
		32.07*df['S'] + \
		30.97*df['P']

	return cmass

#define function to calculate modified aromaticity index
def _calc_NOSC(ct):
	'''
	Calculates the nominal oxidation state of carbon for each formula contained
	within a ``CrossTable`` instance.

	Parameters
	----------
	ct : ft.CrossTable
		``CrossTable`` instance containing the formulae of interest.

	Returns
	-------
	NOSC : pd.Series
		Series of the resulting NOSC values. Length `nF`.

	References
	----------
	LaRowe and van Cappellen (2011), Geochim. Cosmochim. Ac., 67, 2030-2042.
	'''

	#extract chemical compositions, name df for shorthand
	df = ct._chem_comp

	#calculate NOSC (assume z is zero)
	z = 0
	num = 4*df['C'] + df['H'] - 2*df['O'] - 3*df['N'] - 2*df['S'] + 5*df['P'] - z
	NOSC = 4 - num/df['C']

	#name it
	NOSC.name = 'NOSC'

	return NOSC

#define a function to calculate percentages for summary table
def _calc_pct(ct, cols, weights):
	'''
	Calculates percentages for generating a summary table.

	Parameters
	----------
	ct : ft.CrossTable
		``CrossTable`` instance containing the formulae of interest.

	cols : list
		List of column names to be calculated (required input for indexing
		purposes. Columns get mis-aligned if omitted).

	weights : str
		String of weights to use for %RA calculates, either compound
		categories or classes. Inputed from the _gen_sum_tab function.

	Returns
	-------
	pcts : pd.DataFrame
		Resulting dataframe of calculated percentages for each sample.

	Raises
	------
	AssertionError
		If percentages do not add up to 100.

	'''

	#call intensites df for shorthand
	df = ct.intensities

	#calculate totals
	tot_N = np.sum(df > 0)
	tot_int = np.sum(df)

	#extract class / category names
	names = weights.unique()

	#make empty dataframe of results
	pcts = pd.DataFrame(index = ct.sam_names, columns = cols)

	#loop through each class / category and store results
	for n in names:

		#extract compounds within class/category
		ind = np.where(weights == n)[0]
		cforms = df.ix[ind]

		#calculate class/category formula numbers and intensities
		c_N = np.sum(cforms > 0)
		c_int = np.sum(cforms)

		#calculate percentages
		c_pct_N = 100*c_N/tot_N
		c_pct_RA = 100*c_int/tot_int

		#store results
		pcts[n + '_%N'] = c_pct_N
		pcts[n + '_%RA'] = c_pct_RA

	#assert that it all adds up to 200 (sum of N-based and RA-based)
	assert (pcts.sum(axis = 1) - 200 < 1e-6).all(), \
		'Calculated percentages do not add up within 1 part per million!' \
		' Something is not assigned a compound class / category!'

	return pcts

#define a function to check formulae format
def _check_forms(formulae):
	'''
	Function to check that formulae are in the right format for inputting
	and generating a chemical composition table.

	Parameters
	----------
	formulae : list
		List of strings containing formula chemical compositions.

	Raises
	------
	FormulaError
		If formula is not in the right format.
	'''

	#check that all atom aassignments exist in all samples
	atoms = ['C', 'H', 'O', 'N', 'S', 'P']

	if not all([set(atoms) < set(f) for f in formulae]):
		raise FormulaError(
			'Some inputted formulae are missing C, H, O, N, S, or P'
			' assignment!')

	#check that all atom letters are followed by a number
		#loop through each atom, extract number, and store
	for atom in atoms:
		#make regex string
		s = '((?<=' + atom + ')\d{1,2})'

		#check that all formulae have a number for each atom
		if any([re.search(s, f) is None for f in formulae]):
			raise FormulaError(
				'Some inputted formulae are missing %r value!' % atom)

#define a function to check intensity and sample name format
def _check_int(
	intensities,
	formulae = None, 
	sam_names = None):
	'''
	Checks that inputted intensity table and sample names are in the correct
	format.

	Parameters
	----------
	intensities : 2D array-like
		Either 2D np.ndarray or pd.Dataframe of all formula intensities for
		all samples.

	formulae : list
		List of formula names. If `None`, pulls formula names from 
		`intensities` index values. Defaults to `None`.

	sam_names : `None` or list
		List of sample names. If `None`, pulls sample names from column names
		of `intensities`. Defaults to `None`.

	Returns
	-------
	ints : pd.DataFrame
		DataFrame of intensities, with index as formula names and columns
		as sample names.

	forms : list
		Array of formula names. Length `nF`.

	sams : list
		Array of sampe names. Length `nS`.

	Warnings
	--------
	UserWarning
		If `sam_names` = `None` and `intensities` does not contain sample
		names, then names each sample according to column number.

	Raises
	------
	DimError
		If `intensities` is not two dimensional.

	FormulaError
		If no formula names exist, either in `formulae` or `intensities` 
		index.

	LengthError
		If `intensities` is not shape [len(`formulae`) x len(`sam_names`)]

	TypeError
		If `intensities` is not np.ndarray or pd.DataFrame.

	TypeError
		If `intensities` data points are not float.

	'''

	#ensure intensities in the right format
	inttype = type(intensities).__name__
	dim = intensities.ndim

	if inttype not in ['ndarray', 'DataFrame']:
		raise TypeError(
			'`intensities` is type: %r. Must be either np.ndarray or'
			' pd.DataFrame!' % inttype)

	elif dim != 2:
		raise DimError(
			'`intensities` dimensionality is %r. Must be 2!' % dim)

	#if array, make into dataframe and add index + column names
	nF, nS = np.shape(intensities)

	if inttype == 'ndarray':
		#raise formula error if necessary
		if formulae is None:
			raise FormulaError(
				'No formula names exist, either as the `formulae` variable'
				' or as the index of the `intensities` variable.')

		elif len(formulae) != nF:
			raise LengthError(
				'`formulae` and `intensities` are not the same length!')

		#if no sample names, warn and store
		if sam_names is None:
			warnings.warn(
				'No sample names exist, using column numbers as sample names')

			ints = pd.DataFrame(
				intensities, 
				index = formulae)

		elif len(sam_names) != nS:
			raise LengthError(
				'`sam_names` and `intensities` are not the same length!')

		else:
			ints = pd.DataFrame(
				intensities, 
				index = formulae, 
				columns = sam_names)

	else:
		ints = intensities

	#check data type
	if not all(ints.dtypes == float):
		raise TypeError(
			'intensity values are not float. Check datatype!')

	#pull formula names if necessary
	if formulae is None:
		forms = ints.index.values

	elif not all(ints.index == formulae):
		raise FormulaError(
			'formula list in `intensities` index and list of inputted'
			' formulae do not match!')

	else:
		forms = formulae

	#pull sample names if necessary
	if sam_names is None:
		sams = ints.columns.values

	elif not all(ints.columns == sam_names):
		raise SampleError(
			'sample name list in `intensities` columns and list of inputted'
			' samples do not match!')

	else:
		sams = sam_names

	return ints, forms, sams

#define function to combine EnviroOrg samples into single dataframe
def _combine_EO_samples(
	dir_path, 
	file_names):
	'''
	Generates an intensities table for multiple EnviroOrg files.

	Parameters
	----------
	dir_path : str
		String containing the path pointing to the directory containing 
		files to be imported.

	file_names : str or list
		Either a list of strings containing the filenames to be imported or
		the string 'all'. If 'all', method will automatically import all 
		files within the provided directory. Defaults to 'all'.
	
	Returns
	-------
	intensities : pd.DataFrame
		DataFrame containing the MS intensities of each formula within each 
		sample. Shape [`nF` x `nS`].

	Raises
	------
	FileError
		If the directory path does not exist.
	'''

	nS = len(file_names)

	#check directory path
	exists = os.path.isdir(dir_path)

	if not exists:
		raise FileError(
			'Attempting to search for files in directory: %r but it does not'
			' exist!' % dir_path)

	#if file names is all, import all names
	if file_names == 'all':
		file_names = next(os.walk(dir_path))[2]

		#only keep .csv files!
		file_names = [f for f in file_names if '.csv' in f]

	#combine dir_name and name to make absolute paths
	path_names = [dir_path + '/' + s for s in file_names]

	#make empty dataframe
	intensities = pd.DataFrame()

	#loop through, import and store
	for p, f in zip(path_names, file_names):
		#drop .csv from sample name for nomenclature brevity
		s = re.sub('.csv$', '', f)
		
		#import
		s_int = _input_EO_sample(p, s, norm = False)

		#store
		intensities = pd.concat([intensities, s_int], axis = 1)

	#fill missing data with zeros
	intensities.fillna(0, inplace = True)

	return intensities

#define function to drop high HC, high OC, or flier peaks
def _drop_forms(
	forms_all,
	ints_all,
	drop_int_above = None,
	drop_high_HC = True,
	drop_high_OC = True):
	'''
	Drops high HC, high OC, and flier formulae

	Parameters
	----------
	forms_all : list
		List of formula names. Length `nF + nDropped`.

	ints_all : pd.DataFrame
		DataFrame of intensities, with index as formula names and columns
		as sample names.

	Returns
	-------
	chem_comp : pd.DataFrame
		Dataframe of chemical compositions after dropping, split into columns 
		for `C`, `H`, `O`, `N`, `S`, and `P`. Shape [`nF` x 6].

	forms : list
		List of formula names after dropping. Length `nF`.

	ints : pd.DataFrame
		DataFrame of intensities after dropping, with index as formula names
		and columns as sample names.

	Raises
	------
	TypeError 
		If `drop_int_above` is not `float` or `int`.

	ValueError
		If `drop_int_above` is not between 0 and 100 (percentage of tot. int.)
	'''

	#generate chemical composition dataframe
	comps_all = _gen_chem_comp(forms_all)

	#drop compounds with O/C > 1
	if drop_high_OC:

		#find indices of high OC compounds
		hOC_ind = comps_all[comps_all['O'] > 1*comps_all['C']].index

	else:
		hOC_ind = pd.Index([])

	#drop compounds with H/C > 2
	if drop_high_HC:

		#find indices of high HC compounds
		hHC_ind = comps_all[comps_all['H'] > 2*comps_all['C']].index

	else:
		hHC_ind = pd.Index([])

	#drop fliers
	if drop_int_above is None:
		#make empty index
		hInt_ind = pd.Index([])

	elif type(drop_int_above) not in [float, int]:
		raise TypeError(
			'drop_int_above must be float or int')

	elif drop_int_above > 100 or drop_int_above <= 0:
		raise ValuerErro(
			'drop_int_above must be between 0 and 100 (percent of total int.)')

	else:
		#rescale to ensure percent intensity
		ints_all_pct = 100*ints_all/ints_all.sum()

		#find indices with high fractional intensity
		max_int = ints_all_pct.max(axis = 1)
		hInt_ind = max_int[max_int > drop_int_above].index

	#union dropped indices
	dr_ind = hInt_ind.union(hHC_ind.union(hOC_ind))

	#drop all flagged formulae from comps, forms, and ints
	comps = comps_all.drop(dr_ind)
	forms = comps.index.values
	ints = ints_all.drop(dr_ind)

	return comps, forms, ints

#define function to generate chemical composition dataframe
def _gen_chem_comp(formulae):
	'''
	Generates a dataframe of chemical compositions using a string of formulae
	names.

	Parameters
	----------
	formulae : list
		List of strings containing each molecular formula, in the format:
		C(1-99)H(0-99)O(0-99)N(0-9)S(0-9)P(0-9). Length `nF`.

	Returns
	-------
	chem_comp : pd.DataFrame
		Dataframe of chemical compositions, split into columns for `C`, `H`,
		`O`, `N`, `S`, and `P`. Shape [`nF` x 6].
	'''

	#make empty dataframe
	atoms = ['C', 'H', 'O', 'N', 'S', 'P']
	chem_comp = pd.DataFrame(index = formulae, columns = atoms)

	#loop through each atom, extract number, and store
	for atom in atoms:
		#make regex string
		s = '((?<=' + atom + ')\d{1,2})'

		#extract and store
		ccstr = chem_comp.index.str
		chem_comp[atom] = ccstr.extract(s, expand = False).astype(int)


		#chem_comp[atom] = chem_comp.index.str.extract(s, expand = False).astype(int)

	#double check that all assignments are within bounds
	mins = [1, 0, 0, 0, 0, 0]
	maxs = [99, 99, 99, 9, 9, 9]

	for atom, mi, ma in zip(atoms, mins, maxs):
		#raise errors if outside bounds
		mif = chem_comp[atom].min()
		maf = chem_comp[atom].max()

		if mif < mi:
			raise FormulaError(
				'Detected formula with %r = %r. Minimum is %r' 
				% (atom, mif, mi)
				)

		if maf > ma:
			raise FormulaError(
				'Detected formula with %r = %r. Maximum is %r' 
				% (atom, maf, ma)
				)

	return chem_comp

#define function to summarize CrossTable (i.e. calculate percentages)
def _gen_sum_tab(ct):
	'''
	Generates a summary table for a given `CrossTable` instance.

	Parameters
	----------
	ct : ft.CrossTable
		``CrossTable`` instance containing the formulae of interest.

	Returns
	-------
	sum_df : pd.DataFrame
		DataFrame containing the summary information for the sample set
		contained within a given ``CrossTable`` isntance.
	'''

	#define constants
	nF = ct.nF
	nS = ct.nS
	tot_N = np.sum(ct.intensities > 0)
	tot_int = np.sum(ct.intensities)

	#extract categories and classes
	ccats = ct.cmpd_cat.unique()
	ccls = ct.cmpd_class.unique()

	#append %N or for column names
	ccats_N = [c + '_%N' for c in ccats]
	ccats_RA = [c + '_%RA' for c in ccats]

	ccls_N = [c + '_%N' for c in ccls]
	ccls_RA = [c + '_%RA' for c in ccls]

	#make list of additional variables to summarize
	mets = [
		'Tot_N_formulae',
		'AveMass_N',
		'AveMass_RA',
		'NOSC_N',
		'NOSC_RA',
		'AImod_N',
		'AImod_RA'
		]

	#concatenate everything into columns list
	cls_mets = ccls_N + ccls_RA
	cat_mets = ccats_N + ccats_RA
	
	cols = mets + cls_mets + cat_mets

	#make empty dataframe
	sum_df = pd.DataFrame(index = ct.sam_names, columns = cols)

	#populate matrix
	exists = ct.intensities > 0

	sum_df['Tot_N_formulae'] = tot_N

	sum_df['AveMass_N'] = \
		np.sum(exists.multiply(ct.cmpd_mass, axis = 0))/tot_N

	sum_df['AveMass_RA'] = \
		np.sum(ct.intensities.multiply(ct.cmpd_mass, axis = 0))/tot_int

	sum_df['NOSC_N'] = \
		np.sum(exists.multiply(ct.NOSC, axis = 0))/tot_N

	sum_df['NOSC_RA'] = \
		np.sum(ct.intensities.multiply(ct.NOSC, axis = 0))/tot_int

	sum_df['AImod_N'] = \
		np.sum(exists.multiply(ct.AImod, axis = 0))/tot_N

	sum_df['AImod_RA'] = \
		np.sum(ct.intensities.multiply(ct.AImod, axis = 0))/tot_int

	sum_df[cls_mets] = _calc_pct(ct, cls_mets, ct.cmpd_class)
	sum_df[cat_mets] = _calc_pct(ct, cat_mets, ct.cmpd_cat)

	return sum_df

#define function to load EnviroOrg sample output file
def _input_EO_sample(file, sam_name, norm = False):
	'''
	Function to input a single sample from an EnviroOrg output file and to
	return a pd.DataFrame of the data.

	Parameters
	----------
	file : str
		String to the (absolute) path of the file containing the EnviroOrg
		data for a given sample to be summarized.

	sam_name : str
		String of the name of a given file, to be used for cross table.

	norm : boolean
		Tells the function whether or not to re-normalize the relative
		intensities such that they sum to unity. Defaults to `False`.

	Returns
	-------
	sam_int : pd.Series
		Series of the intensities of each formula for a given sample.

	Raises
	------
	FileError
		If the inputted file does not exist or is not in .csv format.
	'''

	#check that file exists and is a file, and check extension is .csv
	exists = os.path.isfile(file)

	if not exists:
		raise FileError(
			'Looking for file named: %r but this file does not exist!' % file)

	elif '.csv' not in file:
		raise FileError(
			'Inputted file is not a .csv!')

	#pre-allocate column names
	fvars = [
		'Exp_mz',
		'Recal_mz',
		'Theor_mz',
		'Error',
		sam_name,
		'S2N',
		'DBE',
		'HC',
		'OC']

	forms = ['C',
		'nC',
		'H',
		'nH',
		'N',
		'nN',
		'O',
		'nO',
		'S',
		'nS',
		'P',
		'nP']

	cols = fvars + forms

	#import file as dataframe (skip first 2 rows -- are garbage in EO output)
	df = pd.read_csv(file,
		skiprows = 2,
		index_col = 0,
		header = None)

	#drop final column (is all Nan in EO output) and name columns
	df = df.loc[:,:21]
	df.columns = cols

	#generate formula names and rename index
	#need to re-order column names to get in CHONSP format (EO output is off)
	fr = ['C','nC','H','nH','O','nO','N','nN','S','nS','P','nP']

	form_names = df[fr].apply(lambda x: ''.join(x.astype('str')), axis = 1)
	df.index = form_names
	df.index.name = 'Formula'

	#check for dupilicates and, if they exist, drop those rows
	df = df[~df.index.duplicated(keep=False)]

	#only retain necessary columns
	sam_int = df[sam_name]

	#rescale if necessary
	if norm:
		sam_int = sam_int/sam_int.sum()

	return sam_int
