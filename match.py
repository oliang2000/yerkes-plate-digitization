#imports
import pandas as pd 
import seaborn as sns
import math
import matplotlib.pyplot as plt
import numpy as np
import utility

#basic info of csv file from APT
file = 'config4_r'
filename = file+'.csv'
df_apt = pd.read_csv(filename) 
df_apt = df_apt[['CentroidRA', 'CentroidDec','Magnitude', 'MagUncertainty']]
min_ra = df_apt['CentroidRA'].min()
max_ra = df_apt['CentroidRA'].max()
min_dec = df_apt['CentroidDec'].min()
max_dec = df_apt['CentroidDec'].max()
print("\n------------IMAGE INFORMATION------------")
print("ra:" ,[min_ra, max_ra], "\n", "dec:", [min_dec, max_dec])
print("------------------------------------------")

#Matching
rp = utility.run_next_step("Fetch data from GAIA")
if rp:
    df_gaia = utility.get_gaia_data (min_ra, max_ra, min_dec, max_dec, 15)
    gaia_name = file+"_gaia.csv"
    export_csv = df_gaia.to_csv (gaia_name, index = None, header=True) 
    plt.scatter(df_apt["CentroidRA"], df_apt["CentroidDec"], alpha = 0.4, s = 3)
    plt.scatter(df_gaia["ra"], df_gaia["dec"], alpha = 0.4, s = 3)
    plt.savefig(file+'_figs.png')
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.legend(bbox_to_anchor=(0, 1), loc='lower left', ncol = 2, fontsize='small')
    title = plt.suptitle("Plate #"+file, fontsize=12)
    plt.savefig(file+'_match.png')
elif not rp:
    df_gaia = pd.read_csv(file+'_gaia.csv')
if utility.run_next_step("Proceed matching"):
    df = utility.match_two_tables(df_gaia, df_apt, file)
    df_filtered = df.query('diff<' + input("Cut difference at? "))
    utility.p_scatter(file, df, df_filtered, 'phot_bp_mean_mag', 'Magnitude', lr2 = True)









