import pandas as pd
import scipy
from scipy.spatial import distance
import numpy as np
import matplotlib.pyplot as plt
import os
import statsmodels.api as sm
from statsmodels.graphics.regressionplots import abline_plot

pnum = "R3170"

#generate queries
apt = pd.read_csv("{0}/{0}.csv".format(pnum))#data from APT
min_ra, max_ra = min(apt['CentroidRA']), max(apt['CentroidRA'])
min_dec, max_dec = min(apt['CentroidDec']), max(apt['CentroidDec'])
query_cal = "select ra, dec, u, g, r, i, z, modelmagerr_u, modelmagerr_g, modelmagerr_r, modelmagerr_i, modelmagerr_z\
from star where \
(ra between {} and {}) and (dec between {} and {}) and modelmag_g < 21.0 and clean = 1".format(min_ra, max_ra, min_dec, max_dec)
query_qso = "select p.ra, p.dec, s.z, snMedian, p.u, p.g, p.r, p.i, p.z, p.modelmagerr_u, \
p.modelmagerr_g, p.modelmagerr_r, p.modelmagerr_i, p.modelmagerr_z from photoobj as p, specobj as s\
where p.specobjid = s.specobjid and snMedian > 7.0 and s.class = 'QSO' and \
(p.ra between {} and {}) and (p.dec between {} and {})".format(min_ra, max_ra, min_dec, max_dec)
print('\n' + '-------------------Plate {}-----------------'.format(pnum) + '\n')
print('Query data from SDSS: http://skyserver.sdss.org/dr16/en/tools/search/sql.aspx' + '\n')
print('SDSS Query for Calibration stars (please save the file as {}_cal.csv):'.format(pnum), '\n' + query_cal + '\n' )
print('SDSS Query for Quasars (please save the file as {}_qso.csv):'.format(pnum), '\n' + query_qso + '\n')


if input("Proceed(y/n): ").strip() == 'y':
    #import tables
    cal = pd.read_csv("{0}/{0}_cal.csv".format(pnum), skiprows=1)#calibration stars from SDSS
    qso = pd.read_csv("{0}/{0}_qso.csv".format(pnum), skiprows=1)#QSOs from sdss


    #match tables
    diff_lim = 0.001 #cutoff distance for mismatches
    if input("Match tables(y/n): ").strip() == 'y':
        df = pd.DataFrame(columns=['ra', 'dec', 'u', 'g','r','i', 'z', 'CentroidRA','CentroidDec','Magnitude','diff'])
        dist = scipy.spatial.distance.cdist(apt[['CentroidRA', 'CentroidDec']], cal[['ra', 'dec']])
        min_dist = np.argmin(dist, axis=1)
        m = 0
        while m < len(apt):
            n = min_dist[m]
            df = df.append({'ra': cal['ra'][n], 'dec': cal['dec'][n], 'u': cal['u'][n], 'g': cal['g'][n], \
                            'r': cal['r'][n], 'i': cal['i'][n], 'z': cal['z'][n],\
                            'CentroidRA': apt['CentroidRA'][m],'CentroidDec': apt['CentroidDec'][m],\
                            'Magnitude': apt['Magnitude'][m],\
                            'diff': dist[m][n]
                           }, ignore_index=True)
            m += 1
        df = df.sort_values('diff').drop_duplicates('ra')
        df = df.sort_values('diff').drop_duplicates('CentroidRA')
        df = df.reset_index(drop=True) 
        df.to_csv("{0}/{0}_match.csv".format(pnum), index = None, header=True)

    else: 
        df = pd.read_csv("{0}/{0}_match.csv".format(pnum))

    df_clean = df[df['diff'] <= diff_lim].copy()
    df_clean['QSO'] = (df_clean.ra.isin(qso.ra) & df_clean.dec.isin(qso.dec))
    df_clean['pg'] = df_clean['g'] + 0.3*(df_clean['u'] - df_clean['g']) - 0.45
    df_qso = df_clean[df_clean['QSO'] == True].copy()

    #plot some figures of the match
    plt.figure(figsize = (5, 4))
    plt.hist(df_clean['diff'])
    plt.title("Distance between APT and SDSS match, R3107")
    plt.savefig("{0}/{0}_dist.png".format(pnum), bbox_inches="tight")

    plt.figure(figsize = (5, 4))
    plt.scatter(df_clean['Magnitude'], df_clean['pg'], color = 'black', s = 5, alpha = 0.7)
    plt.scatter(df_qso['Magnitude'], df_qso['pg'], color = 'red', s = 10)
    plt.title("Magnitude for Plate R3107")
    plt.xlabel("APT mag")
    plt.ylabel("pg")
    plt.savefig("{0}/{0}_mag.png".format(pnum), bbox_inches="tight")

    plt.figure(figsize = (5, 4))
    plt.scatter(df_clean['ra'], df_clean['dec'], color = 'black', s = 5, alpha = 0.7)
    plt.scatter(df_qso['ra'], df_qso['dec'], color = 'red', s = 10)
    plt.title("Position for Plate R3107")
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.savefig("{0}/{0}_coor.png".format(pnum), bbox_inches="tight")


    #measure QSOs and save
    n = 0
    for ind in df_clean[df_clean['QSO'] == True].index:
        n += 1
        qso1_ra, qso1_dec = df_clean.iloc[ind]['ra'], df_clean.iloc[ind]['dec']
        qso1_pg, qso1_mag = df_clean.iloc[ind]['pg'], df_clean.iloc[ind]['Magnitude']
        df_cal = df_clean[(df_clean['ra'] <= qso1_ra + 0.2)&(df_clean['ra'] >= qso1_ra - 0.2)\
            &(df_clean['dec'] <= qso1_dec + 0.2)&(df_clean['dec'] >= qso1_dec - 0.2)
            &(df_clean['pg'] <= qso1_pg + 1)&(df_clean['pg'] >= qso1_pg - 2)].copy()
        df_cal = df_cal.dropna(subset=['pg', 'Magnitude'])
        #create txt file
        print("\nQSO{} has {} calibration stars".format(n, len(df_cal)))   

        os.makedirs('{0}/qso{1}'.format(pnum, n))

        plt.figure(figsize = (5, 4))
        plt.scatter(df_cal['Magnitude'], df_cal['pg'], color = 'black', s = 5, alpha = 0.7)
        plt.scatter([df_clean.iloc[ind]['Magnitude']], [df_clean.iloc[0]['pg']], color = 'red', s = 10)
        plt.title("QSO{} and Calibration star for Plate R3107".format(n))
        plt.xlabel("APT mag")
        plt.ylabel("pg")
        plt.savefig("{0}/qso{1}/{0}_qso_calibration{1}.png".format(pnum, n), bbox_inches="tight")

        #measure, write error + result to bigger folder


        

