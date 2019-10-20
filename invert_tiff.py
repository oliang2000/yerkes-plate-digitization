import os
import math
import numpy as np
from astropy import wcs
from astropy.io import fits
import csv
import sys
from array import array

print("\n---Invert image of a tif file---")
filename = input("Enter file name: ")
hdu1 = fits.open(filename + ".fits")
data = hdu1[0].data #[0]
header = hdu1[0].header
sizex = np.int(header['NAXIS1'])
sizey = np.int(header['NAXIS2'])
bits  = np.int(header['BITPIX'])
max = math.pow(2, bits)

ix = 0
jy = 0
for jy in range(0, sizey):
    for ix in range(0, sizex):
        T = data[jy,ix]
        data[jy,ix] = (max - 1) - T

fits.writeto(filename+"_inv.fits", data, header, overwrite=True)
print("Finished.")
