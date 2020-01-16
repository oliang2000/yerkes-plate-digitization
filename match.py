#imports
import pandas as pd 
import seaborn as sns
import math
import matplotlib.pyplot as plt
import numpy as np
import util
import os


gaia_brightness = 18	#max brightness for gaia query
plate_year = 1901		#year plate was taken
arcsecpp = 1.59			#arcsec per pixel of plate (from astrometry.net)
percent_remove = 50		#percentage of data removed


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
util.analyze_data(df, file, arcsecpp, percent_remove, graph = False) 

#Correct by wiggles
df_apt = util.correct_scanner_wiggle(df_apt, df, 7)

#match again
df = util.match_two_tables(df_gaia, df_apt, file, True)
util.analyze_data(df, file, arcsecpp, percent_remove) 











