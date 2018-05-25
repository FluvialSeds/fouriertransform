'''
Initializes the fouriertransform package.

fouriertransform was created by:

	Jordon D. Hemingway 
	Harvard University
	jordon_hemingway@fas.harvard.edu

source code can be found at:
	
	https://github.com/FluvialSeds/fouriertransform

documentation can be found at:

	http://fouriertransform.readthedocs.io

Version 0.0.8 is current as of 25 Mar 2018.

'''

from __future__ import(
	division,
	print_function,
	)

__version__ = '0.0.8'

__docformat__ = 'restructuredtext en'


#import formulatable classes
from .crosstable import(
	CrossTable,
	)
