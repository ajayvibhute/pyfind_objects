#!/bin/env python
from find_objects import find_objects as fo


#change the name
infile="/Users/avibhute/NRAO/pyfind_objects/find_objects/input.txt"

#create object of image file
imgf=fo.ImageFile()

#init the input
imgf.input(infile)

#process image
imgf.process()

#compute spectral index
imgf.compute_spectral_index()

#save source catalog
imgf.save_source_catalog()

#plot image
imgf.plot_image()

#plot spectra
imgf.plot_spectra()

#save spectra
imgf.save_spectra()
