#todo
#Linear fit pg
#Correct for proper motion
#Remove outliers

PLATE_YEAR = {'R3060': 19150309, 'R3115': 19150602, 'R3095': 19150416, 'R3170': 19150717, 'R3289': 19151104}
import pandas as pd
import scipy
from scipy.spatial import distance
import numpy as np
import matplotlib.pyplot as plt
import os
import statsmodels.api as sm
from statsmodels.graphics.regressionplots import abline_plot
import csv

pnum = input("Enter plate number: ").strip()  #"R3170"
pdate = PLATE_YEAR[pnum]
year_plate = 1915
year_sdss = 2019

#generate queries
apt = pd.read_csv("{0}/{0}.csv".format(pnum))#data from APT
min_ra, max_ra = min(apt['CentroidRA']), max(apt['CentroidRA'])
min_dec, max_dec = min(apt['CentroidDec']), max(apt['CentroidDec'])
query_cal = "select ra, dec, u, g, r, i, z, modelmagerr_u, modelmagerr_g, modelmagerr_r, modelmagerr_i, modelmagerr_z, mjd \
from star where (ra between {:.2f} and {:.2f}) and (dec between {:.2f} and {:.2f}) \
and modelmag_g < 21.0 and clean = 1".format(min_ra, max_ra, min_dec, max_dec)
#star: http://skyserver.sdss.org/dr7/en/help/browser/browser.asp?n=Star&t=V
query_qso = "select p.ra, p.dec, s.z, snMedian, p.u, p.g, p.r, p.i, p.z, p.modelmagerr_u, \
p.modelmagerr_g, p.modelmagerr_r, p.modelmagerr_i, p.modelmagerr_z from photoobj as p, specobj as s \
where p.specobjid = s.specobjid and snMedian > 7.0 and s.class = 'QSO' and \
(p.ra between {:.2f} and {:.2f}) and (p.dec between {:.2f} and {:.2f})".format(min_ra, max_ra, min_dec, max_dec)
print('\n' + '-------------------Plate {}-----------------'.format(pnum) + '\n')
print('Query data from SDSS: http://skyserver.sdss.org/dr16/en/tools/search/sql.aspx' + '\n')
print('SDSS Query for Calibration stars (please save the file as {}_cal.csv):'.format(pnum), '\n' + query_cal + '\n' )
print('SDSS Query for Quasars (please save the file as {}_qso.csv):'.format(pnum), '\n' + query_qso + '\n')


#analyze data
if input("Proceed(y/n): ").strip() == 'y':
    #import tables
    cal = pd.read_csv("{0}/{0}_cal.csv".format(pnum), skiprows = 1)#calibration stars from SDSS
    qso = pd.read_csv("{0}/{0}_qso.csv".format(pnum), skiprows = 1)#QSOs from sdss

    #match tables
    diff_lim = 0.001 #cutoff distance for mismatches, in deg
    if input("Match tables(y/n): ").strip() == 'y':
        df = pd.DataFrame(columns=['ra', 'dec', 'u', 'g','r','i', 'z', 'CentroidRA','CentroidDec','Magnitude','diff'])
        dist = scipy.spatial.distance.cdist(apt[['CentroidRA', 'CentroidDec']], cal[['ra', 'dec']]) #deg
        min_dist = np.argmin(dist, axis=1)
        m = 0
        while m < len(apt):
            n = min_dist[m]
            df = df.append({'ra': cal['ra'][n], 'dec': cal['dec'][n], 'u': cal['u'][n], 'g': cal['g'][n], \
                            'r': cal['r'][n], 'i': cal['i'][n], 'z': cal['z'][n], 'mjd': cal['mjd'][n], \
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
    df_clean = df_clean.dropna(subset=['pg', 'Magnitude'])
    df_clean.reset_index(inplace=True)
    meanra, meandec = np.mean(df.ra), np.mean(df.dec)
    df_clean['radial_dist'] = np.sqrt((df_clean.ra - meanra)**2 + (df_clean.dec - meandec)**2)

    #correct for influence of radial dist on brightness (improve this later by combining stuff?)
    X = df_clean['Magnitude']
    X = sm.add_constant(X)
    y = df_clean['pg']
    model = sm.OLS(y, X).fit()
    df_clean['pred_pg'] = model.predict(X)

    X = df_clean['radial_dist']
    X = sm.add_constant(X)
    y = df_clean['pg'] - df_clean['pred_pg']
    model = sm.OLS(y, X).fit()
    df_clean['pg_cor_by_rad'] = model.predict(X)

    df_qso = df_clean[df_clean['QSO'] == True].copy()

    #plot some figures of the match
    plt.figure(figsize = (5, 4))
    plt.hist(df_clean['diff'])
    plt.title("Distance between APT and SDSS match, {}".format(pnum))
    plt.savefig("{0}/{0}_dist.png".format(pnum), bbox_inches="tight")

    plt.figure(figsize = (5, 4))
    plt.scatter(df_clean['Magnitude'], df_clean['pg'], color = 'black', s = 5, alpha = 0.7)
    plt.scatter(df_qso['Magnitude'], df_qso['pg'], color = 'red', s = 10)
    plt.title("Magnitude for Plate {}".format(pnum))
    plt.xlabel("APT mag")
    plt.ylabel("pg")
    plt.savefig("{0}/{0}_mag.png".format(pnum), bbox_inches="tight")

    plt.figure(figsize = (5, 4))
    plt.scatter(df_clean['ra'], df_clean['dec'], color = 'black', s = 5, alpha = 0.7)
    plt.scatter(df_qso['ra'], df_qso['dec'], color = 'red', s = 10)
    plt.title("Position for Plate {}".format(pnum))
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.savefig("{0}/{0}_coor.png".format(pnum), bbox_inches="tight")


    #measure QSOs and save
    writedata = input("Write data into csv(y/n): ")
    if len(df_clean[df_clean['QSO'] == True]) == 0:
        print("No qualifying quasars found.")
    else:
        print("{} quasars found.".format(len(df_clean[df_clean['QSO'] == True])))
        n = 0
        for ind in df_clean[df_clean['QSO'] == True].index:
            box_size = 0.2 #radius box in which calibration stars are selected
            n += 1
            qso_loc = df_clean.iloc[ind]
            qso_ra, qso_dec, qso_pg, qso_mag = qso_loc['ra'], qso_loc['dec'], qso_loc['pg'], qso_loc['Magnitude']

            qso_diff = qso_loc['diff']
            qso_radial_dist = qso_loc['radial_dist']

            df_cal = df_clean[(df_clean['ra'] <= qso_ra + box_size)&(df_clean['ra'] >= qso_ra - box_size)\
                &(df_clean['dec'] <= qso_dec + box_size)&(df_clean['dec'] >= qso_dec - box_size)
                &(df_clean['pg'] <= qso_pg + 1)&(df_clean['pg'] >= qso_pg - 2)].copy()


            #linear fit
            X = df_cal['Magnitude']
            X = sm.add_constant(X)
            y = df_cal['pg']
            model = sm.OLS(y, X).fit()
            df_cal['pg_predictions'] = model.predict(X) + df_cal['pg_cor_by_rad'] #new


            #plots
            qsopath = '{0}/qso{1}'.format(pnum, n)
            if not os.path.exists(qsopath):
                os.makedirs(qsopath)

            fig = abline_plot(model_results = model, color = 'black')
            ax = fig.axes[0]
            ax.scatter(df_cal['Magnitude'], df_cal['pg'])
            ax.scatter([qso_mag], [qso_pg], color = 'red')
            plt.title("Linear fit for calibration stars, qso{0}, {1}".format(n, pnum))
            plt.xlabel("APT mag")
            plt.ylabel("pg")
            plt.savefig("{0}/qso{1}/{0}_qso{1}_calibration.png".format(pnum, n), bbox_inches="tight")
            plt.clf()

            df_cal['res'] = df_cal['pg'] - df_cal['pg_predictions']
            plt.hist(df_cal['res'])
            plt.title("Histogram for residuals, qso{0}, {1}".format(n, pnum))
            plt.savefig("{0}/qso{1}/{0}_qso{1}_residual_hist.png".format(pnum, n), bbox_inches="tight")
            plt.clf()

            plt.scatter(df_cal['Magnitude'], df_cal['res'])
            plt.title("Magnitude vs. residuals, qso{0}, {1}".format(n, pnum))
            plt.savefig("{0}/qso{1}/{0}_qso{1}_residual_plot{1}.png".format(pnum, n), bbox_inches="tight")
            plt.clf()

            #create txt file
            with open("{0}/{0}.txt".format(pnum), 'a') as f:
                f.write("\n------------------QSO {}-------------------\n".format(n))
                f.write("\nQSO{} has {} calibration stars".format(n, len(df_cal)) + "\n")
                low_16, high_84 = (np.quantile(df_cal['res'], 0.16), np.quantile(df_cal['res'], 0.84))

                
                qso_mag_measured = df_cal[df_cal['QSO'] == True]['pg_predictions'].values[0]

                f.write("Measured pg from plate({}): {:.4f}".format(year_plate, qso_mag_measured) + "\n")
                f.write("Measured pg from SDSS({}): {:.4f}".format(year_sdss, qso_pg) + "\n")
                f.write("The difference between two measurements are: {:.4f}".format(qso_pg - qso_mag_measured) + "\n")
                f.write("Error for this value (sigma value for residuals): {:.4f}".format((high_84 - low_16)/2) + "\n")

            #write into table
            if writedata == "y":
                with open('qso.csv', 'a', newline='') as csvfile: 
                    f = csv.writer(csvfile) 
                    f.writerow([qso_ra, qso_dec, df_clean.iloc[ind]['mjd'], list(qso[qso['ra'] == qso_ra]['z'])[0], qso_pg, \
                        pnum, pdate, qso_mag_measured, qso_mag_measured - qso_pg, "{:.4f}".format((high_84 - low_16)/2), \
                        qso_diff, qso_radial_dist])
 


        

