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

#df_apt = df_apt.rename(columns={"CentroidRA": "ra", "CentroidDec": "dec"})


#if utility.run_next_step("Proceed fetching data from GAIA"):
#	df_gaia = utility.get_gaia_data (min_ra, max_ra, min_dec, max_dec, 15)
#	gaia_name = file+"_gaia.csv"
#	export_csv = df_gaia.to_csv (gaia_name, index = None, header=True)
#sns.regplot('ra', 'dec', data=df_gaia, fit_reg = False)
#sns.regplot('CentroidRA', 'CentroidDec', data=df_apt, fit_reg = False).set_title("Position match")
#fig.savefig(file+'_figs.png')

##a working oneeee
#df_gaia = pd.read_csv("ca397_600_gaia.csv") 
#sns_plot = sns.regplot('ra', 'dec', data=df_gaia, fit_reg = False)
#fig = sns_plot.get_figure()
#fig.savefig(file+'_figs.png')

#anorther working one
#df_gaia = pd.read_csv("ca397_600_gaia.csv") 
#plt.scatter(df_apt["CentroidRA"], df_apt["CentroidDec"])
#plt.scatter(df_gaia["ra"], df_gaia["dec"])
#plt.savefig(file+'_figs.png')

df_gaia = pd.read_csv("ca397_600_gaia.csv") 
plt = utility.p_scatter(df_gaia, df_apt, "ra", "CentroidRA")
plt.savefig(file+'_figs.png')









