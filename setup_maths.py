from distutils.core import setup, Extension

module1 = Extension('maths',
	sources = ['mathsmodule.c'])

setup (name = 'Maths',
	version = '1.0',
	author = 'Daniel F. Simola',
	author_email = 'dfsimola@fastmail.fm',
	description = 'Implements math and statistics functions',
	ext_modules = [module1])