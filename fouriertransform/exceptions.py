'''
This module contains all the exceptions for the fouriertransform package.
'''

#define a core exception class for subclassing
class ftException(Exception):
	'''
	Root Exception class for the fouriertransform package. Do not call 
	directly.
	'''
	pass


class FormulaError(ftException):
	'''
	Detected formula is outside of the assignment bounds.
	'''
	pass 