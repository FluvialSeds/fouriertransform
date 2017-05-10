#import packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

#function for extracting only peaks contained within a particular sample
def extract_sample(df, name, info_col_max = 16, normalize = True):
	'''
	Function to extract information for a particular sample from the
	combined dataframe.

	Parameters
	----------
	df : pd.DataFrame
		dataframe containing all samples and all detected peaks for the
		dataset.

	name : str
		String of the name of the sample to be extracted. Name must be a
		column within `df`.

	info_col_max : int
		Last column containing molecular information (i.e. before intensities
		for each sample in the dataset begin). Defaults to 16 (i.e. the value
		for the Ganges headwater dataset).

	normalize : boolean
		Tells the function to normalize intensities or not. If True,
		intensitites will range from [0, 1]. Defaults to `True`.
	
	Returns
	-------
	res : pd.DataFrame
		dataframe for a given sample, sorted by increasing intensity.
	'''

	#keep only existing peaks and extract necessary info
	df_nz = df.loc[df[name] > 0]
	info = df_nz.iloc[:,0:info_col_max]
	ri = df_nz[name]

	#normalize if necessary
	if normalize:
		ri = ri/np.max(ri)

	#rename intensity to be consistent across samples
	ri.name = 'intensity'

	#concatenate info and intensity
	res = pd.concat([info, ri], axis = 1)

	#sort result
	res.sort_values(by = 'intensity', inplace = True)

	return res

#extract common peaks
def extract_common_peaks(df, names, frac_cutoff = 0.33):
	'''
	Function to extract only the peaks (i.e. rows) that are common to multiple
	samples, as defined by `frac_cutoff`.

	Parameters
	----------
	df : pd.DataFrame
		dataframe containing all samples and all detected peaks for the
		dataset.

	names : list
		List of strings of the sample names to be considered. All strings in
		`names` must be contained in the column names of `df`.

	frac_cutoff : float
		Fraction of total samples in which a peak must appear in order to be
		retained. Defaults to 0.33 (i.e. must be in greater than 1/3 of all
		samples.)


	Returns
	-------
	comm_df : pd.DataFrame
		dataframe retaining only common peaks.
	'''

	#set values and make intensity only dataframe
	df_int = df[names]
	n, m = np.shape(df_int) #n = total peaks, m = samples
	n_cut = np.int(np.ceil(frac_cutoff * m))

	#get the number of samples that each peak is contained in
	n_sams = np.sum(df_int > 0.0, axis = 1)

	#find indicies where n_sams is above cutoff
	inds = np.where(n_sams >= n_cut)[0]

	#only keep those rows
	comm_df = df.ix[inds]

	return comm_df

#function for plotting van krevelen
def plot_vk_intensity(res, ax, s = 15, ecolor = 'w', cmap = 'YlGnBu', log = True, vmin = 0, vmax = 1):
	'''
	Function for plotting van krevelen intensities for a given sample.

	Parameters
	----------
	res : pd.DataFrame
		Sample info resulting from `extract_sample` function.

	ax : plt.axis
		axis object to be plotted on

	Returns
	-------
	ax : plt.axis
		Populated figure axes
	'''

	#plot relative intensity (color coded by log intensity)
	if log:
		#log transform colors
		c = np.log10(res['intensity'].values)
	else:
		c = res['intensity'].values

	ax.scatter(
		res['OC'],
		res['HC'],
		s = s,
		edgecolor = ecolor,
		linewidth = 0.5,
		c = c,
		vmin = vmin,
		vmax = vmax,
		cmap = cmap)

	return ax

def plot_vk_difference(res1, res2, ax, s = 15, ecolor = 'w', plot_S = False):
	'''
	Function to plot a vk of the difference between two samples
	'''

	#find peaks in res1 but not in res2
	r1 = res1['intensity']
	r2 = res2['intensity']
	r1.name = 'res1'
	r2.name = 'res2'

	com = pd.concat([r1,r2],axis=1)
	x = com[np.isnan(com['res2'])]

	#extract rows in res1 only
	res1only = res1.ix[x.index]

	#extract N and S containing compounds
	cho = res1only[(res1only['N'] == 0) & (res1only['S'] == 0)]
	chon = res1only[res1only['N'] > 0]
	chos = res1only[res1only['S'] > 0]

	#plot categorical
	ax.scatter(
		cho['OC'],
		cho['HC'],
		s = s,
		c = [0/255, 58/255, 112/255],
		linewidth = 0.5,
		marker = 'o',
		edgecolor = ecolor,
		label = 'CHO',
		zorder = 2)

	ax.scatter(
		chon['OC'],
		chon['HC'],
		s = s,
		c = [255/255, 189/255, 23/255],
		linewidth = 0.5,
		marker = 'o',
		edgecolor = ecolor,
		label = 'CHON',
		zorder = 1)

	if plot_S:
		ax.scatter(
			chos['OC'],
			chos['HC'],
			s = s,
			c = [164/255, 16/255, 52/255],
			linewidth = 0.5,
			marker = 'o',
			edgecolor = ecolor,
			label = 'CHOS',
			zorder = 3)

	return ax

#define function to calculate spearman rho values
def plot_vk_spearman(env, form, names, param, ax, s = 15, ecolor = 'w', cmap = 'coolwarm'):
	'''
	Calculates the spearman correlation between the `param` column of `env` and
	the formulas in `form` and plots in van krevelen space.
	'''

	f = form.loc[:,names]

	#calculate rho and p-values for each formula, combine into dataframe
	rhos, pvals = stats.spearmanr(f.T, env[param])
	r = rhos[-1,:-1]
	p = pvals[-1,:-1]

	df = pd.DataFrame([r,p],
		index = ['rho','pval'],
		columns = form.index).T

	corr_df = pd.concat([form[['HC','OC']], df], axis = 1)

	#retain only statistically significant formulas
	sig_corr = corr_df[corr_df['pval'] <= 0.05]

	#sort formulas by ascending rho
	sig_corr_sort = sig_corr.sort_values('rho')
	n = len(sig_corr_sort)

	#plot
	sc = ax.scatter(
		sig_corr_sort['OC'],
		sig_corr_sort['HC'],
		s = s,
		edgecolor = ecolor,
		linewidth = 0.5,
		c = sig_corr_sort['rho'],
		vmin = -1.0,
		vmax = 1.0,
		cmap = cmap)

	#write n on the figure and scale colorbar
	# ax.set_clim([-1.0, 1.0])

	mu = sig_corr_sort['rho'].abs().mean()
	sig = sig_corr_sort['rho'].abs().std()

	ax.text(0.05, 0.05,r'n = %.0f; mu = %.2f; sig = %.2f' % (n, mu, sig),
		transform = ax.transAxes,
		fontsize = 8,
		horizontalalignment = 'left',
		verticalalignment = 'bottom')

	return ax, sc

#import data
all_data = pd.read_csv('../input_datasets/kolyma_isotopes.csv',index_col = 0)
all_forms = pd.read_csv('../input_datasets/kolyma_forms.csv',index_col = 0)
names = all_data.index.values

###################
## PLOT FIGURE 1 ##
###################

# make figure
fig1, ax1 = plt.subplots(3,1,figsize=(3.74,9.05),sharex=True,sharey=True)

# set limits and labels
ax1[0].set_xlim([0.0,1.0])
ax1[0].set_ylim([0.0,2.0])
ax1[1].set_ylabel('H/C')
ax1[2].set_xlabel('O/C')

plt.tight_layout()

#extract samples but do not normalize
iwt = extract_sample(all_forms,'IWT',normalize=False)
kol = extract_sample(all_forms,'KOL',normalize=False)

#normalize both samples by the biggest peak
denom = np.max([iwt['intensity'].max(),kol['intensity'].max()])
iwt['intensity'] = iwt['intensity']/denom
kol['intensity'] = kol['intensity']/denom

#calculate vmin
vmin = np.log10(np.min([iwt['intensity'].min(), kol['intensity'].min()]))

# Figure 1a: IWT sample
ax1[0] = plot_vk_intensity(iwt, ax1[0], vmin = vmin, vmax = 0)

# Figure 1b: KOL sample
ax1[1] = plot_vk_intensity(kol, ax1[1], vmin = vmin, vmax = 0)

# Figure 1c: IWT - KOL difference
ax1[2] = plot_vk_difference(iwt, kol, ax1[2])

#save figure
fig1.savefig(
	'../output_figures/Kolyma_intensity_vks.pdf',
	dpi = 300,
	transparent = True,
	frameon = False)

###################
## PLOT FIGURE 2 ##
###################

#extract peaks common to all stream samples
str_forms = extract_common_peaks(all_forms, names, frac_cutoff = 1.0)

#rescale formula intensities
str_forms[names] = str_forms[names]/str_forms[names].sum()

#make formula and environmental parameter dataframe
env_params = all_data.loc[:,['d13C','D14C']]

#make figure
fig2, ax2 = plt.subplots(1, 1, figsize = (3.74, 3.74))

# set limits and labels
ax2.set_xlim([0.0,1.0])
ax2.set_ylim([0.0,2.0])
ax2.set_ylabel('H/C')
ax2.set_xlabel('O/C')

plt.tight_layout()

#plot
ax2, _ = plot_vk_spearman(env_params, str_forms, names, 'D14C', ax2)

#save figure
fig2.savefig(
	'../output_figures/Kolyma_Spearman_vks.pdf',
	dpi = 300,
	transparent = True,
	frameon = False)