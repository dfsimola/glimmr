#! /usr/bin/env python

from distutils.core import setup, Extension
from glimmrAccessories import *

version = """\nsetup.py, version %s

 GlimmrHet: genomic locus identification using multiply maped reads and a diploid (heterozygous) reference genome
 
 Requires Unix/Mac OS X/CYGWIN with Python 2.5+, bowtie v1 or v2 
 (preferably 64-bit for memory allocation > 4GB)
 
 Daniel F. Simola, PhD (simola@upenn.edu)
 Laboratory of Shelley L. Berger, PhD
 University of Pennsylvania
 October 2015
 
 Copyright (c) 2015, Daniel F. Simola and Shelley L. Berger, University of
 Pennsylvania. You may not use this file except in compliance with the terms of our
 License, GNU GENERAL PUBLIC LICENSE v3.0, a copy of which is located in the LICENSE file.
"""%(GLIMMR_VERSION)

"""%s"""%(version)

module1 = Extension('snipercore',
	sources = ['snipercore.c'])

setup (name = 'snipercore',
	version = '1.0',
	author = 'Daniel F. Simola',
	author_email = 'simola@upenn.edu',
	description = 'Implements likelihood computation in C for Glimmr',
	ext_modules = [module1])