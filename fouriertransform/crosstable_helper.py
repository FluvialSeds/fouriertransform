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
	'_check_forms',
	'_check_int',
	'_combine_EO_samples',
	'_gen_chem_comp',
	'_input_EO_sample',
	]

import numpy as np
import os
import pandas as pd
import re
import warnings

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

	#calculate DBEai
	DBEai = 1 + df['C'] - 0.5*df['O'] - df['S'] - 0.5*df['H']

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
	'''

	#extract chemical compositions and AImod and combine
	ccs = ct._chem_comp
	ais = ct.AImod
	df = pd.concat([ccs, ais], axis = 1)

	#extract list of formulae stored in CrossTable instance
	ccat = pd.Series(index = df.index, name = 'cmpd_cat')

	#calculate indices for each category

	#aliphatic, high OC
	ah_ind = df[
		(df['H'] >= 1.5*df['C']) & 
		(df['O'] >= 0.5*df['C']) &
		(df['N'] == 0)].index

	#aliphatic, low OC
	al_ind = df[
		(df['H'] >= 1.5*df['C']) & 
		(df['O'] < 0.5*df['C']) &
		(df['N'] == 0)].index

	#condensed aromatic
	ca_ind = df[(
		df['AImod'] >= 0.67)].index

	#polyphenolic, high OC
	ph_ind = df[
		(df['AImod'] < 0.67) & 
		(df['AImod'] >= 0.5) & 
		(df['O'] >= 0.5*df['C'])].index

	#polyphenolic, low OC
	pl_ind = df[
		(df['AImod'] < 0.67) & 
		(df['AImod'] >= 0.5) & 
		(df['O'] < 0.5*df['C'])].index

	#unsaturatied and phenolic, high OC
	uh_ind = df[
		(df['AImod'] < 0.5) & 
		(df['H'] < 1.5*df['C']) & 
		(df['O'] >= 0.5*df['C'])].index

	#unsaturatied and phenolic, low OC
	ul_ind = df[
		(df['AImod'] < 0.5) & 
		(df['H'] < 1.5*df['C']) & 
		(df['O'] < 0.5*df['C'])].index

	#peptide-like
	pe_ind = df[
		(df['H'] >= 1.5*df['C']) & 
		(df['N'] > 0)].index

	#set values for compounds
	ccat[ah_ind] = 'aliphatic_highOC'
	ccat[al_ind] = 'aliphatic_lowOC'

	ccat[ca_ind] = 'condensed_aromatic'
	
	ccat[ph_ind] = 'polyphenolic_highOC'
	ccat[pl_ind] = 'polyphenolic_lowOC'

	ccat[uh_ind] = 'unsaturated_phenolic_highOC'	
	ccat[ul_ind] = 'unsaturated_phenolic_lowOC'

	ccat[pe_ind] = 'peptide_like'

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
def _check_int(intensities, formulae = None, sam_names = None):
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
		List of formula names. Length `nF`.

	sams : list
		List of sampe names. Length `nS`.

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
def _combine_EO_samples(dir_path, file_names):
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

	#combine dir_name and name to make absolute paths
	path_names = [dir_path + '/' + s for s in file_names]

	#make empty dataframe
	intensities = pd.DataFrame()

	#loop through, import and store
	for p, f in zip(path_names, file_names):
		#drop .csv from sample name for nomenclature brevity
		s = re.sub('.csv$', '', f)
		
		#import
		s_int = _input_EO_sample(p, f, norm = False)

		#store
		intensities = pd.concat([intensities, s_int], axis = 1)

	#fill missing data with zeros and rescale so biggest peak = 100
	intensities.fillna(0, inplace = True)

	m = intensities.max().max()
	intensities = intensities*100/m

	return intensities

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
		chem_comp[atom] = chem_comp.index.str.extract(s).astype(int)

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

	#only retain necessary columns
	sam_int = df[sam_name]

	#rescale if necessary
	if norm:
		sam_int = sam_int/sam_int.sum()

	return sam_int












