#imports
import pandas as pd 
import seaborn as sns
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
import statsmodels.api as sm
import utility
#basic info of csv file from APT
file = 'ca397_600'
filename = file+'.csv'
df_apt = pd.read_csv(filename) 
df_apt = df_apt[['CentroidRA', 'CentroidDec','Magnitude', 'MagUncertainty']]
min_ra = df_apt['CentroidRA'].min()
max_ra = df_apt['CentroidRA'].max()
min_dec = df_apt['CentroidDec'].min()
max_dec = df_apt['CentroidDec'].max()
print("\n------------IMAGE INFORMATION------------")
print("ra:" ,[min_ra, max_ra])
print("dec:", [min_dec, max_dec])


if utility.run_next_step("Proceed fetching data from gaia"):
	dt_gaia = utility.get_gaia_data (min_ra, max_ra, min_dec, max_dec, 15)
	gaia_name = file+"_gaia.csv"
	export_csv = dt_gaia.to_csv (gaia_name, index = None, header=True)
