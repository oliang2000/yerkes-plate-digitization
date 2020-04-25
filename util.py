#make extra file for functions generating figures
#add SDSS query

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
    Inverts a tiff file (with single layer!)
    Inputs:
        filename(str): path of the file
    Outputs:
        inverted .fits file
    '''
    hdu1 = fits.open(filename + ".fits")
    data = hdu1[0].data #[0] if many layer
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


def process_file(file, gaia_brightness, plate_year):
    '''
    Reads table from apt with selected columns, extract Gaia data according to that.
    Inputs:
        file(str): filename of information extracted from APT
        gaia_brightness(float): upper limit for brightness of extracted sources in Gaia
        plate_year(int): year the plate was taken, to correct for proper motion
    Outputs:
        .txt file containing basic info about the plate
    Returns
        a tuple containing two pandas.dataFrame (df_apt, df_gaia)
    '''
    
    assert os.path.isfile(MY_PATH + file + '/' + file + '.csv'), "File does not exist!"
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
        df_gaia = get_gaia_data(file, min_ra, max_ra, min_dec, max_dec, gaia_brightness, plate_year)
    return (df_apt, df_gaia)



def get_gaia_data(file, min_ra, max_ra, min_dec, max_dec, brightness, plate_year):
    '''
    Fetches source data from gaia by given parameters, and correct for proper motion.
    Inputs:
        file(str): name of the plate, for saving the data as local file
        min_ra, max_ra, min_dec, max_dec(float): constraints on coordinates of the extracted sources
        brightness(float): constraint on the brightness of the extracted sources
        plate_year(int): year the plate was taken, to correct for proper motion
    Outputs:
        .csv file containing the extracted Gaia data
    Returns:
        pandas.DataFrame with data extracted from Gaia
    '''
    print("Getting data from GAIA...")

    query = "SELECT ALL gaia_source.source_id, gaia_source.ra, gaia_source.dec,\
    gaia_source.pmra, gaia_source.pmdec, gaia_source.phot_g_mean_mag, gaia_source.phot_bp_mean_mag\
    FROM gaiadr2.gaia_source WHERE CONTAINS\
    (POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS'," + \
    str((min_ra + max_ra)/2) + "," + str((min_dec + max_dec)/2) + "," + \
    str(max_ra - min_ra) + "," + str(max_dec - min_dec) + "))=1 AND \
    (gaiadr2.gaia_source.phot_bp_mean_mag<=" + str(brightness) + ")" 

    job = Gaia.launch_job_async(query)
    r = job.get_results()
    df_gaia = r.to_pandas()
    #corrected to supposed coordinates @ time of plate
    #gaia(pmra, pmdec) units: mas/yr - milli-second of arc per year (1 milliarcsecond = 10^-3 arcsecond)
    df_gaia['pmra'].fillna(0, inplace = True)
    df_gaia['pmdec'].fillna(0, inplace = True)
    df_gaia['ra_cor'] = df_gaia['ra'] - (2015 - plate_year) * df_gaia['pmra'] / 1000 / 3600
    df_gaia['dec_cor'] = df_gaia['dec'] - (2015 - plate_year) * df_gaia['pmdec'] / 1000 / 3600
    #save file
    gaia_name = file + "_gaia.csv"
    export_csv = df_gaia.to_csv(MY_PATH + file + '/' + gaia_name, index = None, header=True)
    print("File saved as", gaia_name)
    return df_gaia
   


def match_two_tables(gaia, apt, file, again = False):
    '''
    Matches Gaia and APT data.
    Inputs:
        gaia, apt(pandas.DataFrame): Gaia and APT data of the sources
        file(str): name of the plate
        again(boolean): if True, the matching happens again despite the existence of a matched .csv
    Outputs:
        a .csv file with matched data
        write in .txt file with more information on number of matches
    Returns:
        pandas.DataFrame with matched sources
    '''
    if os.path.isfile(MY_PATH + file + '/' + file+'_match.csv') and not again:
        print ("Matching file data exists.")
        df = pd.read_csv(MY_PATH + file + '/' + file+'_match.csv')
    else:
        print("Start matching GAIA and APT data...")
        df = pd.DataFrame(columns=['ra', 'dec', 'source_id', 'ra_cor','dec_cor','phot_bp_mean_mag', 'phot_g_mean_mag'
            ,'pmra','pmdec', 'CentroidRA','CentroidDec','Magnitude','diff'])
        dist = scipy.spatial.distance.cdist(apt[['CentroidRA', 'CentroidDec']], gaia[['ra_cor', 'dec_cor']])
        min_dist = np.argmin(dist, axis=1)

        m = 0
        while m < len(apt):
            n = min_dist[m]
            df = df.append({'ra': gaia['ra'][n], 'dec': gaia['dec'][n], 'source_id': gaia['source_id'][n],\
                            'ra_cor': gaia['ra_cor'][n], 'dec_cor': gaia['dec_cor'][n],\
                            'phot_bp_mean_mag': gaia['phot_bp_mean_mag'][n],\
                            'phot_g_mean_mag': gaia['phot_g_mean_mag'][n],\
                            'pmra': gaia['pmra'][n], 'pmdec': gaia['pmdec'][n],\
                            'CentroidRA': apt['CentroidRA'][m],'CentroidDec': apt['CentroidDec'][m],\
                            'Magnitude': apt['Magnitude'][m],\
                            'diff': dist[m][n]
                           }, ignore_index=True)
            m += 1
        df = df.sort_values('diff').drop_duplicates('ra_cor')
        df = df.sort_values('diff').drop_duplicates('CentroidRA')
        df = df.reset_index(drop=True) 
        df['del_ra'] = (df['CentroidRA'] - df['ra_cor']) * 3600 #degree -> arcsec
        df['del_dec'] = (df['CentroidDec'] - df['dec_cor']) * 3600 
        export_csv = df.to_csv(MY_PATH + file + '/' + file + "_match.csv", index = None, header=True)
        print("Matching finished.")
    file1 = open(MY_PATH + file + '/' + file + "_stats.txt", "a") 
    file1.write("--------------------MATCH----------------------\n")
    file1.write("APT: " + str(len(apt)) + "\n")
    file1.write("GAIA: " + str(len(gaia)) + "\n")
    file1.write("MATCHED: " + str(len(df)) + "\n")
    file1.close() 
    return df


def correct_scanner_wiggle(file, apt, df, order, cut_percentage):
    '''
    Correct CentroidRA and CentroidDec(APT) by polyfitting wiggle.
    Inputs:
        file(str): name of the file 
        apt, df(pandas.Dataframe): APT and matched file on the sources
        order(int): order of the polynomial to fit the curve
        cut_percentage(int): percentage of data being cut before fitting polynomial
    Outputs:
        two .png showing the polyfit of the wiggles
    Returns:
        pandas.DataFrame with corrected APT data
    '''
    df = df[np.abs(df["diff"]) <= np.percentile(np.abs(df["diff"]), 100 - cut_percentage)]

    z = np.polyfit(df['CentroidDec'], df['del_dec'], order)
    p = np.poly1d(z)
    apt['CentroidDec'] = apt['CentroidDec'] - p(apt['CentroidDec'])/3600
    plt.clf()
    plt.scatter(df['CentroidDec'], p(df['CentroidDec']), s = 1)
    plt.scatter(df['CentroidDec'], df['del_dec'], s = 1, c = df['diff'], alpha = 0.4)
    plt.xticks(rotation=45)
    plt.savefig(MY_PATH + file + '/deldec_cor' + file + '.png', dpi = 300, bbox_inches="tight")

    z = np.polyfit(df['CentroidRA'], df['del_ra'], order)
    p = np.poly1d(z)
    apt['CentroidRA'] = apt['CentroidRA'] - p(apt['CentroidRA'])/3600
    plt.clf()
    plt.scatter(df['CentroidRA'], p(df['CentroidRA']), s = 2)
    plt.scatter(df['CentroidRA'], df['del_ra'], s = 1, c = df['diff'], alpha = 0.4)
    plt.xticks(rotation=45)
    plt.savefig(MY_PATH + file + '/delra_cor' + file + '.png', dpi = 300, bbox_inches="tight")

    print("Corrected for scanner wiggles.")

    return apt


def analyze_data(df, file, dpp, cut_percentage, graph = True):
    if graph == True:
        df.hist(column='diff', bins = 15)
        plt.savefig(MY_PATH + file + '/hist1_' + file + '.png', bbox_inches="tight")



        df.hist(column='phot_bp_mean_mag', bins = 30)
        plt.savefig(MY_PATH + file + '/brightness_' + file + '.png', bbox_inches="tight")






        #get data with 68% in RA and Dec
        len_before = len(df)
        df = df[np.abs(df["diff"]) <= np.percentile(np.abs(df["diff"]), 100 - cut_percentage)] #
        len_after = len(df)
        #figures
        df.hist(column='diff', bins = 15)
        plt.xticks(rotation=45)
        plt.savefig(MY_PATH + file + '/hist2_' + file + '.png', bbox_inches="tight")

        df.plot.scatter(x = "ra_cor", y = "del_ra", c = "phot_bp_mean_mag", s = 3, colormap='viridis')
        plt.xticks(rotation=45)
        plt.savefig(MY_PATH + file + '/delra_' + file + '.png', dpi = 300, bbox_inches="tight")

        df.plot.scatter(x = "dec_cor", y = "del_dec", c = "phot_bp_mean_mag", s = 3, colormap='viridis')
        plt.xticks(rotation=45)
        plt.savefig(MY_PATH + file + '/deldec_' + file + '.png', dpi = 300, bbox_inches="tight")

        df.plot.scatter(x = "Magnitude", y = "phot_bp_mean_mag", c = "diff", s = 3, colormap='viridis')
        plt.savefig(MY_PATH + file + '/mag_' + file + '.png', dpi = 300, bbox_inches="tight")

        plt.clf()
        plt.quiver( - df['ra_cor'], - df['dec_cor'], df['del_ra'], df['del_dec'])
        plt.xticks(rotation=45)
        plt.savefig(MY_PATH + file + '/vector_' + file + '.png', dpi = 300, bbox_inches="tight")
    else:
        len_before = len(df)
        df = df[np.abs(df["diff"]) <= np.percentile(np.abs(df["diff"]), 100 - cut_percentage)] #
        len_after = len(df)
    #write txt with data info
    file1 = open(MY_PATH + file + '/' + file + "_stats.txt", "a") 
    file1.write("---------------------RMS-----------------------\n")
    file1.write("Removed " + str(cut_percentage) + "% of data\n")
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
        


def graph_data_distr(file, apt, gaia):
    '''
    Plot the position of the sources in Gaia and APT
    #try overlap this!!
    '''
    apt.plot.scatter(x = "CentroidRA", y = "CentroidDec", s = 3)
    plt.savefig(MY_PATH + file + '/apt_' + file + '.png', bbox_inches="tight")
    gaia.plot.scatter(x = "ra_cor", y = "dec_cor", c = 'red', s = 3)
    plt.savefig(MY_PATH + file + '/gaia_' + file + '.png', bbox_inches="tight")


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



#SDSS query
# SELECT ALL
#    p.objid,p.ra,p.dec,p.u,p.g,p.r,p.i,p.z,
#    p.run, p.rerun, p.camcol, p.field,
#    s.specobjid, s.class, s.z as redshift,
#    s.plate, s.mjd, s.fiberid
# FROM PhotoObj AS p
#    JOIN SpecObj AS s ON s.bestobjid = p.objid
# WHERE 
#    p.ra BETWEEN  348.029335 AND 349.898045
#    AND p.dec BETWEEN 5.750175 AND 7.636144
#distance in space