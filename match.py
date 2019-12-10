#imports
import pandas as pd 
import seaborn as sns
import math
import matplotlib.pyplot as plt
import numpy as np
import util
import os


#Import APT data and get GAIA data
file = input("Enter file name: ").strip()
df_apt, df_gaia = util.process_file(file)
#Graph APT and GAIA 
util.graph_matching(file, df_apt, df_gaia)

#Matching & create graphs
df = util.match_two_tables(df_gaia, df_apt, file)
util.analyze_data(df, file)

#df_filtered = df.query('diff<' + input("Cut difference at? "))
#util.p_scatter(file, df, df_filtered, 'phot_bp_mean_mag', 'Magnitude', lr2 = True)
#util.p_scatter(file, df, df_filtered, 'ra', 'del_ra')
#util.p_scatter(file, df, df_filtered, 'dec', 'del_dec')










