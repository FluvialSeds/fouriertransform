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


class DimError(ftException):
	'''
	If object dimensionality is not what it should be.
	'''
	pass


class FileError(ftException):
	'''
	If attempting to import a file that does not exist.
	'''
	pass

class FormulaError(ftException):
	'''
	Detected formula is outside of the assignment bounds.
	'''
	pass


class LengthError(ftException):
	'''
	If two things that should be the same length are not.
	'''
	pass


class SampleError(ftException):
	'''
	Attempting to work with a sample that isn't in the current sample set.
	'''
	pass

