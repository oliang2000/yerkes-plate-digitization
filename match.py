#imports
import pandas as pd 
import seaborn as sns
import math
import matplotlib.pyplot as plt
import numpy as np
import utility
import os
my_path = os.path.abspath(__file__)
my_path = my_path.replace("match.py", "")

#basic info of csv file from APT
file = input("Enter file name: ")
df_apt = pd.read_csv(my_path + file + '/' + file + '.csv') 
df_apt = df_apt[['CentroidRA', 'CentroidDec','Magnitude', 'MagUncertainty']]
min_ra = df_apt['CentroidRA'].min()
max_ra = df_apt['CentroidRA'].max()
min_dec = df_apt['CentroidDec'].min()
max_dec = df_apt['CentroidDec'].max()
print("\n------------IMAGE INFORMATION------------")
print(" RA:" ,[min_ra, max_ra], "\n Dec:", [min_dec, max_dec])
print("-----------------------------------------")

#Get data from gaia
getdata = utility.run_next_step("Fetch data from GAIA")
if getdata:
    df_gaia = utility.get_gaia_data(file, min_ra, max_ra, min_dec, max_dec, 18) 
    utility.graph_matching(file, df_apt, df_gaia, (max_ra-min_ra)/(max_dec-min_dec))
elif not getdata:
    df_gaia = pd.read_csv(my_path + file + '/' + file+'_gaia.csv')
    utility.graph_matching(file, df_apt, df_gaia, (max_ra-min_ra)/(max_dec-min_dec))
#Matching
if utility.run_next_step("Proceed matching"):
    df = utility.match_two_tables(df_gaia, df_apt, file)
    df_filtered = df.query('diff<' + input("Cut difference at? "))
    utility.p_scatter(file, df, df_filtered, 'phot_bp_mean_mag', 'Magnitude', lr2 = True)









