from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
import pandas as pd
import scipy
import numpy as np
import seaborn as sns
import statsmodels.api as sm
from astropy import wcs
from astropy.io import fits
import csv
import sys
from array import array
import os
import math

MY_PATH = os.path.abspath(__file__).replace("util.py", "")

def invert_tiff(filename):
    '''
    Inverts a tiff file.
    '''
    hdu1 = fits.open(filename + ".fits")
    data = hdu1[0].data #[0]
    header = hdu1[0].header
    sizex = np.int(header['NAXIS1'])
    sizey = np.int(header['NAXIS2'])
    bits  = np.int(header['BITPIX'])
    max = math.pow(2, bits)

    ix = 0
    jy = 0
    for jy in range(0, sizey):
        for ix in range(0, sizex):
            T = data[jy,ix]
            data[jy,ix] = (max - 1) - T

    fits.writeto(filename + ".fits", data, header, overwrite=True)

def process_file(file, gaia_brightness = 20):
    '''
    Reads table from apt with selected columns, outputs the table
      prints information about the coordinates.
    '''
    df_apt = pd.read_csv(MY_PATH + file + '/' + file + '.csv')
    df_apt = df_apt[['CentroidRA', 'CentroidDec','Magnitude', 'MagUncertainty']]
    min_ra = df_apt['CentroidRA'].min()
    max_ra = df_apt['CentroidRA'].max()
    min_dec = df_apt['CentroidDec'].min()
    max_dec = df_apt['CentroidDec'].max()
    min_brightness = df_apt['Magnitude'].min()
    max_brightness = df_apt['Magnitude'].max()

    file1 = open(MY_PATH + file + '/' + file + "_stats.txt","w") 
    file1.write("\n------------------IMAGE INFO-------------------\n")
    file1.write(" RA: " + str(min_ra) + ", " + str(max_ra) + "\n")
    file1.write(" Dec: " + str(min_dec) + ", " + str(max_dec) + "\n")
    file1.write(" Mag: " + str(min_brightness) + ", " + str(max_brightness) + "\n")
    file1.close() 

    print("\n------------------IMAGE INFO-------------------")
    print(" RA:" ,[min_ra, max_ra], "\n Dec:", [min_dec, max_dec])
    print(" Mag:" ,[max_brightness, min_brightness])
    print("-----------------------------------------------")

    gaia_name = file+"_gaia.csv"
    if os.path.isfile(MY_PATH + file + '/' + gaia_name):
        print ("GAIA data exists.")
        df_gaia = pd.read_csv(MY_PATH + file + '/' + file+'_gaia.csv')
    else:
        df_gaia = get_gaia_data(file, min_ra, max_ra, min_dec, max_dec, gaia_brightness)
    return (df_apt, df_gaia)



def get_gaia_data (file, min_ra, max_ra, min_dec, max_dec, brightness):
    '''
    Fetches star data from gaia inside the box given and upper limit of brightness.
    '''
    print("Getting data from GAIA...")
    query = "SELECT ALL gaia_source.source_id,gaia_source.ra,gaia_source.dec,\
    gaia_source.pmra,gaia_source.pmdec,gaia_source.phot_bp_mean_mag FROM gaiadr2.gaia_source WHERE CONTAINS\
    (POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS'," + \
    str((min_ra + max_ra)/2) + "," + str((min_dec + max_dec)/2) + "," + \
    str(max_ra - min_ra) + "," + str(max_dec - min_dec) + "))=1 AND \
    (gaiadr2.gaia_source.phot_bp_mean_mag<=" + str(brightness) + ")" 
    job = Gaia.launch_job_async(query)
    r = job.get_results()
    df_gaia = r.to_pandas()
    gaia_name = file+"_gaia.csv"
    export_csv = df_gaia.to_csv(MY_PATH + file + '/' + gaia_name, index = None, header=True)
    print("File saved as", gaia_name)
    return df_gaia
   


def match_two_tables(gaia, apt, file):
    if os.path.isfile(MY_PATH + file + '/' + file+'_match.csv'):
        print ("Matching file data exists.")
        df = pd.read_csv(MY_PATH + file + '/' + file+'_match.csv')
    else:
        print("Start matching GAIA and APT data...")
        df = pd.DataFrame(columns=['ra','dec','phot_bp_mean_mag','pmra','pmdec',\
                                   'CentroidRA','CentroidDec','Magnitude','diff'])
        dist = scipy.spatial.distance.cdist(apt[['CentroidRA', 'CentroidDec']], gaia[['ra', 'dec']])
        min_dist = np.argmin(dist, axis=1)

        m = 0
        while m < len(apt):
            n = min_dist[m]
            df = df.append({'ra': gaia['ra'][n], 'dec': gaia['dec'][n],\
                            'phot_bp_mean_mag': gaia['phot_bp_mean_mag'][n],\
                            'pmra': gaia['pmra'][n], 'pmdec': gaia['pmdec'][n],\
                            'CentroidRA': apt['CentroidRA'][m],'CentroidDec': apt['CentroidDec'][m],\
                            'Magnitude': apt['Magnitude'][m],\
                            'diff': dist[m][n]
                           }, ignore_index=True)
            m+=1
        df = df.sort_values('diff').drop_duplicates('ra')
        df = df.sort_values('diff').drop_duplicates('CentroidRA')
        df = df.reset_index(drop=True) 
        df['del_ra'] = df.apply(lambda row: (row.CentroidRA - row.ra) * 3600, axis = 1) #degree -> arcsec
        df['del_dec'] = df.apply(lambda row: (row.CentroidDec - row.dec) * 3600, axis = 1)
        df['err_ra'] = df['del_ra'] + df['pmra']/10
        df['err_dec'] = df['del_dec'] + df['pmdec']/10
        #gaia(pmra, pmdec): mas/yr - milli-second of arc per year (1 milliarcsecond = 10^-3 arcsecond)
        #apt(del_ra, del_dec): arcs per ~100 years

        export_csv = df.to_csv(MY_PATH + file + '/' + file + "_match.csv", index = None, header=True)

        print("Matching finished.")

    file1 = open(MY_PATH + file + '/' + file + "_stats.txt", "a") 
    file1.write("--------------------MATCH----------------------\n")
    file1.write("APT: " + str(len(apt)) + "\n")
    file1.write("GAIA: " + str(len(gaia)) + "\n")
    file1.write("MATCHED: " + str(len(df)) + "\n")
    file1.close() 
    return df



def analyze_data(df, file, dpp, cut_percentage):

    #testtt
    df.hist(column='diff', bins = 15)
    plt.savefig(MY_PATH + file + '/hist1_' + file + '.png')

    #get data with 68% in RA and Dec
    len_before = len(df)
    df = df[np.abs(df["diff"]) <= np.percentile(np.abs(df["diff"]), cut_percentage)] #
    len_after = len(df)

    #figures
    df.hist(column='diff', bins = 15)
    plt.savefig(MY_PATH + file + '/hist2_' + file + '.png')

    df.plot.scatter(x = "ra", y = "del_ra", c = "phot_bp_mean_mag", s = 3, colormap='viridis')
    plt.xticks(rotation=45)
    plt.savefig(MY_PATH + file + '/delra_' + file + '.png', dpi = 300)

    df.plot.scatter(x = "dec", y = "del_dec", c = "phot_bp_mean_mag", s = 3, colormap='viridis')
    plt.xticks(rotation=45)
    plt.savefig(MY_PATH + file + '/deldec_' + file + '.png', dpi = 300)

    df.plot.scatter(x = "phot_bp_mean_mag", y = "Magnitude", c = "diff", s = 3, colormap='viridis')
    plt.savefig(MY_PATH + file + '/mag_' + file + '.png', dpi = 300)

    plt.clf()
    plt.quiver(df['ra'], df['dec'], df['err_ra'], df['err_dec'])
    plt.savefig(MY_PATH + file + '/vector_' + file + '.png', dpi = 300)

    #write txt with data info
    file1 = open(MY_PATH + file + '/' + file + "_stats.txt", "a") 
    file1.write("---------------------RMS-----------------------\n")
    file1.write("Removed " + ("%.3f" % (1 - len_after/len_before)) + "% of data\n")
    file1.write(get_rms(df, "del_ra", "arcsec", dpp, "pixel"))
    file1.write(get_rms(df, "del_dec", "arcsec", dpp, "pixel"))
    file1.close() 



def get_rms(df, col, unit = "", conv_f = None, conv_unit = None):
    '''
    Calculates rms for a column of a datatable, returns a string
      in the format "colname: data"
    '''
    rms = ((df[col] - df[col].mean()) ** 2).mean() ** .5
    if conv_f == None: 
        return (col + ": " + ("%.3f" % rms) + " "+ unit + "\n")
    else:
        conv = rms/conv_f
        return (col + ": " + ("%.3f" % rms) + " "+ unit + \
            "(" + ("%.3f" % conv) + " " + conv_unit + ")" + "\n")
        


def graph_matching(file, apt, gaia):
    apt.plot.scatter(x = "CentroidRA", y = "CentroidDec", s = 3)
    plt.savefig(MY_PATH + file + '/apt_' + file + '.png')
    gaia.plot.scatter(x = "ra", y = "dec", c = 'red', s = 3)
    plt.savefig(MY_PATH + file + '/gaia_' + file + '.png')


def p_scatter(file, df1, df2, x, y, xlim =[0,0] , ylim=[0,0], lr1 = False, lr2 = False):
    '''
    This function uses matplotlib.pyplot to graph two functions and fix axis (calculated if not given).
    df2 should be the smaller subset.
    '''
    plt.clf()
    axes = plt.gca()
    plt.scatter(df1[x], df1[y], alpha = 0.4, s = 4)
    plt.scatter(df2[x], df2[y], alpha = 0.4, s = 4)
    if xlim == [0,0]:
        xlim = [df1[x].min()-(df1[x].max()-df1[x].min())*0.05, df1[x].max() + (df1[x].max()-df1[x].min())*0.05]
    if ylim == [0,0]:
        ylim = [df1[y].min()-(df1[y].max()-df1[y].min())*0.05, df1[y].max() + (df1[y].max()-df1[y].min())*0.05]
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    plt.xlabel(x)
    plt.ylabel(y)
    if lr1 == True:
        model = sm.OLS(df1[y], sm.add_constant(df1[x])).fit()
        pred = df1[x]*model.params[x]+model.params['const']
        plt.plot(df1[x],pred,'#2678b2')
    if lr2 == True:
        model = sm.OLS(df2[y], sm.add_constant(df2[x])).fit()
        pred = df2[x]*model.params[x]+model.params['const']
        plt.plot(df2[x],pred,'#fd7f28')
    title = plt.suptitle(file + ': ' + x + ' vs. ' + y, fontsize=12)
    plt.savefig(MY_PATH + file + '/' + file + '_'+ x + '_' + y + '.png')
    plt.close()


#distance in space