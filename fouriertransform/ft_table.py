'''
Script to generate summary table from raw enviro org files
'''

# QUESTIONS FOR DAVID:

# CONDENSED AROMATIC CURRENTLY INCLUDES O/C > 1 COMPOUNDS!
# NO CHONS? WHAT ABOUT CYSTEINE, METHIONINE, ETC.?
# WHAT ABOUT PHOSPHORUS CONTAINING COMPOUNDS?
# COMPOUND CATEGORY IDENTIFICATION ALGORITHM?



import numpy as np
import os
import pandas as pd

#path to folder containing all files to be compiled
path = 'drake_enviro_org'
rename = False

#make a list of all file names in that folder
files = next(os.walk(path))[2]

#extract only .csv.xls EnviroOrg files
eo_files = [f for f in files if '.csv.xls' in f]

fvars = [
	'Noise',
	'Exp mz',
	'Recal mz',
	'Theor mz',
	'Error',
	'Rel Abund',
	'S2N',
	'DBE',
	'HC',
	'OC',
	]

forms = [
	'C',
	'nC',
	'H',
	'nH',
	'N',
	'nN',
	'O',
	'nO',
	'S',
	'nS',
	]


df = pd.read_csv(path + '/' + eo_files[5],
	skiprows = 2,
	header = None)

#drop P columns and Null column
df = df.loc[:,:19]
df.columns = fvars + forms

#calculate formula names
names = df[forms].apply(lambda x: ''.join(x.astype('str')), axis = 1)


def input_sample(file, sam_id = None, norm = False):
	'''
	Function to generate a summary file for a given sample.

	Parameters
	----------
	file : str
		String to the path of the file containing the EnviroOrg data for a
		given sample to be summarized.

	sam_id : str or `None`
		String of the (possibly shorthand) sample ID to be used for the
		name of the relative abundance column. If `None`, defaults to
		file path string.

	norm : boolean
		Tells the function whether or not to re-normalize the relative
		intensities such that they sum to unity. Defaults to `False`.

	Returns
	-------
	df : pd.DataFrame
		DataFrame of the summary statistics for a given sample.
	'''

	#perform checks

	#TO ADD CHECKS LATER!!!

	#check sample id string
	if sam_id is None:
		sam_id = file

	#pre-allocate column names for imported data. These are the columns of
	# the default EnviroOrg spreadsheet
	fvars = [
		'Noise',
		'Exp_mz',
		'Recal_mz',
		'Theor_mz',
		'Error',
		sam_id,
		'S2N',
		'DBE',
		'HC',
		'OC',
		]

	#molecular formula columns -- 'cC', 'hH', etc. is simply a column of that
	# letter, while 'C', 'H', etc. is the number of each atom
	forms = [
		'cC',
		'C',
		'hH',
		'H',
		'nN',
		'N',
		'oO',
		'O',
		'sS',
		'S',
		]

	#import data into dataframe
	df = pd.read_csv(file,
		skiprows = 2,
		header = None)

	# TO CHECK WITH DAVID: HERE I'M NOT CONSIDERING PHOSPHORUS!!!

	#drop P columns and Null column
	df = df.loc[:,:19]
	df.columns = fvars + forms

	#calculate formula names
	names = df[forms].apply(lambda x: ''.join(x.astype('str')), axis = 1)

	#reset and name index
	df.index = names
	df.index.name = 'Formula'

	#drop unnecessary columns
	df.drop(
		[
			'Noise',
			'Exp_mz',
			'Recal_mz',
			'Error',
			'S2N',
			'cC',
			'hH',
			'nN',
			'oO',
			'sS'
			],
		axis = 1,
		inplace = True
		)

	if norm:
		df[sam_id] = df[sam_id]/df[sam_id].max()


	return df


def summarize(sams, path = None):
	'''
	Inputs multiple samples into single summary table.

	Parameters
	----------
	sams : list
		List of strings pointing to the files for each sample to be
		included in summary table.

	path : `None` or str
		String pointing to the directory containing files. If files are in
		working directory, use `None`. Defaults to `None`.

	Returns
	-------
	sum_tab : pd.DataFrame
		Dataframe of the resulting summary table.
	'''

	n = len(sams)

	if path is not None:
		#append path to each file name
		files = [path + '/' + s for s in sams]

	else:
		files = sams

	#import first file
	df = input_sample(files[0], sams[0])

	#loop through remaining files and add results to df
	for f, s in zip(files[1:], sams[1:]):
		#import file
		df2 = input_sample(f, s)

		#append
		df = df.combine_first(df2)

	#fill missing values with zero
	df.fillna(0, inplace = True)

	#rescale so largest peak is 100
	m = df[sams].max().max()
	df[sams] = df[sams]*100/m

	#add additional summary columns
	df['NC'] = df['N']/df['C']
	df['SC'] = df['S']/df['C']
	df['AImod'] = calc_AImod(df)
	df['Category'] = calc_cat(df)
	df['Class'] = calc_class(df)

	return df













def calc_AImod(df):
	'''
	Calculates the modified Aromaticity Index following Koch and Dittmar 2006.

	Parameters
	----------
	df : pd.DataFrame
		Dataframe containing the formulas to calculate AI

	Returns
	-------
	AImod : pd.Series
		Series of the resulting AImod values
	'''

	#calculate DBEai
	DBEai = 1 + df['C'] - 0.5*df['O'] - df['S'] - 0.5*df['H']

	#calculate Cai
	Cai = df['C'] - 0.5*df['O'] - df['N'] - df['S']

	#calculate AImod
	AImod = DBEai / Cai

	#make zero where negative
	AImod[AImod < 0] = 0

	return AImod

def calc_cat(df):
	'''
	Calculates the compound category of each formula according to H/C, O/C,
	and AI values.

	Parameters
	----------
	df : pd.DataFrame
		Dataframe containing the formulas to calculate categories

	Returns
	-------
	cats : pd.Series
		Series of the resulting categories
	'''

	cat = pd.Series(index = df.index)

	#calculate indices for each category

	#aliphatic, low OC
	al_ind = df[
		(df['HC'] >= 1.5) & 
		(df['OC'] < 0.5) &
		(df['N'] == 0)].index

	#aliphatic, high OC
	ah_ind = df[
		(df['HC'] >= 1.5) & 
		(df['OC'] >= 0.5) &
		(df['N'] == 0)].index

	#condensed aromatic
	ca_ind = df[(
		df['AImod'] >= 0.67)].index

	#polyphenolic, low OC
	pl_ind = df[
		(df['AImod'] < 0.67) & 
		(df['AImod'] >= 0.5) & 
		(df['OC'] < 0.5)].index

	#polyphenolic, high OC
	ph_ind = df[
		(df['AImod'] < 0.67) & 
		(df['AImod'] >= 0.5) & 
		(df['OC'] >= 0.5)].index

	#unsaturatied and phenolic, low OC
	ul_ind = df[
		(df['AImod'] < 0.5) & 
		(df['HC'] < 1.5) & 
		(df['OC'] < 0.5)].index

	#unsaturatied and phenolic, high OC
	uh_ind = df[
		(df['AImod'] < 0.5) & 
		(df['HC'] < 1.5) & 
		(df['OC'] >= 0.5)].index

	#peptide-like
	pe_ind = df[
		(df['HC'] >= 1.5) & 
		(df['N'] > 0)].index

	# #sugars
	# su_ind = df[
	# 	(df['HC'] > 2.0) |
	# 	(df['OC'] > 0.9)].index

	#set values for compounds
	cat[ca_ind] = 'condensed_aromatic'
	
	cat[pl_ind] = 'polyphenolic_lowOC'
	cat[ph_ind] = 'polyphenolic_highOC'
	
	cat[ul_ind] = 'unsaturated_phenolic_lowOC'
	cat[uh_ind] = 'unsaturated_phenolic_highOC'

	cat[al_ind] = 'aliphatic_lowOC'
	cat[ah_ind] = 'aliphatic_highOC'
	
	cat[pe_ind] = 'peptide_like'
	# cat[su_ind] = 'sugar'

	return cat


def calc_class(df):
	'''
	Calculates the compound class of each formula as CHO, CHON, CHOS, CHONS.

	Parameters
	----------
	df : pd.DataFrame
		Dataframe containing the formulas to calculate classes

	Returns
	-------
	cla : pd.Series
		Series of the resulting classes
	'''

	cla = pd.Series(index = df.index)

	#calculate indices for each category
	cho_ind = df[(df['N'] == 0) & (df['S'] == 0)].index
	chon_ind = df[(df['N'] != 0) & (df['S'] == 0)].index
	chos_ind = df[(df['N'] == 0) & (df['S'] != 0)].index
	chons_ind = df[(df['N'] != 0) & (df['S'] != 0)].index

	#set values for Class
	cla[cho_ind] = 'CHO'
	cla[chon_ind] = 'CHON'
	cla[chos_ind] = 'CHOS'
	cla[chons_ind] = 'CHONS'

	return cla





