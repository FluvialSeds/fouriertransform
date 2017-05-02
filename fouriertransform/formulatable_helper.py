'''
This module contains helper functions for the FormulaTable class.
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
	'_gen_chem_comp',
	]

import numpy as np
import pandas as pd
import warnings

#import exceptions
from .exceptions import(
	FormulaError,
	)

#define function to calculate modified aromaticity index
def _calc_AImod(formtab):
	'''
	Calculates the modified armoaticity index for each formula contained
	within a ``FormulaTable`` instance.

	Parameters
	----------
	formtab : ft.FormulaTable
		``FormulaTable`` instance containing the formulae of interest.

	Returns
	-------
	AImod : pd.Series
		Series of the resulting AImod values. Length `nF`.

	References
	----------
	Koch and Dittmar (2006), Rapid Comm. Mass. Spec., 20, 926-932.
	'''

	#extract chemical compositions, name df for shorthand
	df = formtab._chem_comp

	#calculate DBEai
	DBEai = 1 + df['C'] - 0.5*df['O'] - df['S'] - 0.5*df['H']

	#calculate Cai
	Cai = df['C'] - 0.5*df['O'] - df['N'] - df['S']

	#calculate AImod and name it
	AImod = DBEai / Cai
	AImod.name = 'AImod'

	#make zero where negative
	AImod[AImod < 0] = 0

	return AImod

#define function to calculate formula compound category
def _calc_category(formtab):
	'''
	Calculates the compound category of each formula contained within a
	``FormulaTable`` instance as:

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
	formtab : ft.FormulaTable
		``FormulaTable`` instance containing the formulae of interest.

	Returns
	-------
	ccat : pd.Series
		Series of the resulting categories. Length `nF`.

	References
	----------
	Santi-Temkiv et al. (2013), PLoS One, doi:10.1371/journal.pone.0053550.
	'''

	#extract chemical compositions and AImod and combine
	ccs = formtab._chem_comp
	ais = formtab.AImod
	df = pd.concat([ccs, ais], axis = 1)

	#extract list of formulae stored in FormulaTable instance
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
def _calc_class(formtab):
	'''
	Calculates the compound class of each formula contained within a
	``FormulaTable`` instance as CHO, CHON, CHOS, or CHONS.

	Parameters
	----------
	formtab : ft.FormulaTable
		``FormulaTable`` instance containing the formulae of interest.

	Returns
	-------
	cclass : pd.Series
		Series of the resulting classes. Length `nF`.
	'''

	#extract chemical compositions, name df for shorthand
	df = formtab._chem_comp

	#extract list of formulae stored in FormulaTable instance
	cclass = pd.Series(index = df.index, name = 'cmpd_class')

	#calculate indices using chemical composition logic
	cho_ind = df[(df['N'] == 0) & (df['S'] == 0)].index
	chon_ind = df[(df['N'] != 0) & (df['S'] == 0)].index
	chos_ind = df[(df['N'] == 0) & (df['S'] != 0)].index
	chons_ind = df[(df['N'] != 0) & (df['S'] != 0)].index

	#set values for compound class
	cclass[cho_ind] = 'CHO'
	cclass[chon_ind] = 'CHON'
	cclass[chos_ind] = 'CHOS'
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
	
	Warnings
	--------
	UserWarning
		If inputted formulae contain phosphorus (which is not assigned in this
		software), warns that P will be dropped from analysis.

	Raises
	------
	FormulaError
		If formula is not in the right format.
	'''

	#check that all atom aassignments exist in all samples
	atoms = ['C', 'H', 'O', 'N', 'S']

	if not all([set(atoms) < set(f) for f in formulae]):
		raise FormulaError(
			'Some inputted formulae are missing C, H, O, N, or S assignment!')

	#check that all atom letters are followed by a number
		#loop through each atom, extract number, and store
	for atom in atoms:
		#make regex string
		s = '((?<=' + atom + ')\d{1,2})'

		#check that all formulae have a number for each atom
		if any([re.search(s, f) is None for f in formulae]):
			raise FormulaError(
				'Some inputted formulae are missing %r value!' % atom)

	#warn if phosphorus exists
	if  any(['P' < set(f) for f in formulae]):
		warnings.warn(
			'Some inputted formulae have been assigned phosphorus. P'
			' assignment is dubious and will be dropped from analysis!'
			' (will still remain in sample name)')

#define a function to check intensity and sample name format
def _check_int(inteisities, sam_names = None)

	return ints, sams


#define function to generate chemical composition dataframe
def _gen_chem_comp(formulae):
	'''
	Generates a dataframe of chemical compositions using a string of formulae
	names.

	Parameters
	----------
	formulae : list
		List of strings containing each molecular formula, in the format:
		C(1-45)H(1-92)O(1-25)N(0-4)S(0-2). Length `nF`.

	Returns
	-------
	chem_comp : pd.DataFrame
		Dataframe of chemical compositions, split into columns for `C`, `H`,
		`O`, `N`, and `S`. Shape [`nF` x 5].
	'''

	#make empty dataframe
	atoms = ['C', 'H', 'O', 'N', 'S']
	chem_comp = pd.DataFrame(index = formulae, columns = atoms)

	#loop through each atom, extract number, and store
	for atom in atoms:
		#make regex string
		s = '((?<=' + atom + ')\d{1,2})'

		#extract and store
		chem_comp[atom] = chem_comp.index.str.extract(s).astype(int)

	#double check that all assignments are within bounds
	mins = [1, 1, 1, 0, 0]
	maxs = [45, 92, 25, 4, 2]

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





