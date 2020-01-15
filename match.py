#imports
import pandas as pd 
import seaborn as sns
import math
import matplotlib.pyplot as plt
import numpy as np
import util
import os

file = input("Enter file name: ").strip()

#Invert tiff file
#invert_tiff(file)

#SExtractor

#Import APT data and get GAIA data
df_apt, df_gaia = util.process_file(file, gaia_brightness = 20)

#Graph APT and GAIA 
util.graph_matching(file, df_apt, df_gaia)

#Matching & create graphs
df = util.match_two_tables(df_gaia, df_apt, file)
util.analyze_data(df, file, 1.59, 68) #arcsecpp











