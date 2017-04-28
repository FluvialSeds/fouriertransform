'''
This module contains helper functions for the FormulaTable class.
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

#import exceptions
from .exceptions import(
	)

#define function to calculate formula class
def _calc_class(formtab):
	'''
	Calculates the compound class of each formula contained within a
	`formulatable` instance as CHO, CHON, CHOS, or CHONS.

	Parameters
	----------
	formtab : ft.FormulaTable
		``FormulaTable`` instance containing the formulae of interest.

	Returns
	-------
	classes : pd.Series
		Series of the resulting classes
	'''

	#extract list of formulae stored in FormulaTable instance
	classes = pd.Series(index = formtab.formulae, name = 'cmpd_class')

	#calculate indices for each category using regex
	cho_ind = classes.index.str.contains('N0S0')
	chon_ind = classes.index.str.contains('N[1-4]S0')
	chos_ind = classes.index.str.contains('N0S[12]')
	chons_ind = classes.index.str.contains('N[1-4]S[12')

	#set values for compound class
	cla[cho_ind] = 'CHO'
	cla[chon_ind] = 'CHON'
	cla[chos_ind] = 'CHOS'
	cla[chons_ind] = 'CHONS'

	return classes