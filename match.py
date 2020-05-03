#imports
import pandas as pd 
import seaborn as sns
import math
import matplotlib.pyplot as plt
import numpy as np
import util
import os


gaia_brightness = 21	#max brightness for gaia query
plate_year = 1916		#year plate was taken (to calculate proper motion)
arcsecpp = 1.39			#arcsec per pixel of plate (from astrometry.net)
cut_percent = 50		#percentage of data not being cut
polycor_order = 5		#order of correction for wiggles

font = {#'family' : 'normal',
        'weight' : 'medium',
        'size'   : 12.5}

plt.rc('font', **font)

file = input("Enter file name: ").strip()

#Invert tiff file
#invert_tiff(file)

#SExtractor
#functions here!

#Import extracted data(APT and get GAIA data
df_apt, df_gaia = util.process_file(file = file, gaia_brightness = gaia_brightness, plate_year = plate_year)

#Graph extracted and GAIA 
util.graph_data_distr(file, df_apt, df_gaia)

#Matching & create graphs
df = util.match_two_tables(df_gaia, df_apt, file)
util.analyze_data(df, file, arcsecpp, cut_percent, graph = False) 

#Correct by wiggles
df_apt = util.correct_scanner_wiggle(file, df_apt, df, polycor_order, cut_percent)

#match again
df = util.match_two_tables(df_gaia, df_apt, file, again = True)
util.analyze_data(df, file, arcsecpp, cut_percent) 











