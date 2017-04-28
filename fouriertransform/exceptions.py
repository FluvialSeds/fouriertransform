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


class ExampleError(ftException):
	'''
	Array-like object is not in the right form (e.g. strings).
	'''
	pass 