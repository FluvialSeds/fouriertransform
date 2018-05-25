from setuptools import setup

def readme():
	with open('README.rst') as f:
		return f.read()

setup(name='fouriertransform',
	version='0.0.8',
	description='FT-ICR MS peak assignment and data analysis',
	long_description=readme(),
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'License :: Free for non-commercial use',
		'Programming Language :: Python',
		'Programming Language :: Python :: 3',
		'Programming Language :: Python :: 2',
		'Programming Language :: Python :: 3.5',
		'Programming Language :: Python :: 2.7',
		'Topic :: Scientific/Engineering'
	],
	url='https://github.com/FluvialSeds/fouriertransform',
	# download_url='https://github.com/FluvialSeds/fouriertransform/tarball/0.0.8',
	download_url='https://github.com/FluvialSeds/fouriertransform/archive/0.0.8.tar.gz',
	keywords=[
		'geochemistry',
		'FT-ICR MS',
		'formula assignment',
		'high-resolution mass spectrometry',
		'carbon cycle'
	],
	author='Jordon D. Hemingway',
	author_email='jordon_hemingway@fas.harvard.edu',
	license='GNU GPL Version 3',
	packages=['fouriertransform'],
	install_requires=[
		'matplotlib',
		'numpy',
		'pandas',
		'scipy'
	],
	test_suite='nose.collector',
	tests_require=['nose'],
	include_package_data=True,
	zip_safe=False)