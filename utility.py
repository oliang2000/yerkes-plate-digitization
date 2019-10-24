from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
import pandas as pd
import scipy
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import os
my_path = os.path.abspath(__file__)
my_path = my_path.replace("utility.py", "")

def get_gaia_data (file, min_ra, max_ra, min_dec, max_dec, brightness):
    '''
    Fetches star data from gaia inside the box given and upper limit of brightness.
    Inputs:
        min_ra: minimum value of RA
        max_ra: maximum value of RA
        min_dec: minimum value of Dec
        max_dec: maximum value of Dec
    Outputs: pandas datatable
    '''
    query = "SELECT ALL gaia_source.source_id,gaia_source.ra,gaia_source.dec,\
    gaia_source.pmra,gaia_source.pmdec,gaia_source.phot_bp_mean_mag FROM gaiadr2.gaia_source WHERE CONTAINS\
    (POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS'," + \
    str((min_ra + max_ra)/2) + "," + str((min_dec + max_dec)/2) + "," + \
    str(max_ra - min_ra) + "," + str(max_dec - min_dec) + "))=1 AND \
    (gaiadr2.gaia_source.phot_bp_mean_mag<=" + str(brightness) + ")" 
    job = Gaia.launch_job_async(query)
    r = job.get_results()
    df = r.to_pandas()
    gaia_name = file+"_gaia.csv"
    export_csv = df.to_csv(my_path + file + '/' + gaia_name, index = None, header=True)
    return df

def run_next_step(question):
    '''
    Ask whether user would like to proceed to next step.
    '''
    inpt = question + " (y/n)? "
    g = input(inpt)
    if g == "y":
        return True
    elif g == "n":
        return False
    else:
        print("Please enter y/n.")
        return run_next_step(question)
        
def graph_matching(file, apt, gaia, ratio):
    fig = plt.figure(figsize=(60,20/ratio))
    title = fig.suptitle("Plate #"+file, fontsize = 25)
    fig.subplots_adjust(top=0.85, wspace=0.3, hspace=0.5)
    ax1 = fig.add_subplot(131)
    plt.scatter(apt["CentroidRA"], apt["CentroidDec"], s = 3)
    plt.xlabel("RA")
    plt.ylabel("Dec")
    ax1.title.set_text('APT')
    ax2 = fig.add_subplot(132)
    plt.scatter(gaia["ra"], gaia["dec"], s = 3, c = '#fd7f28')
    plt.xlabel("RA")
    plt.ylabel("Dec")
    ax2.title.set_text('GAIA')
    ax3 = fig.add_subplot(133)
    plt.scatter(apt["CentroidRA"], apt["CentroidDec"], alpha = 0.4, s = 3)
    plt.scatter(gaia["ra"], gaia["dec"], alpha = 0.4, s = 3)
    plt.xlabel("RA")
    plt.ylabel("Dec")
    ax3.title.set_text('APT vs. GAIA')
    fig.savefig(my_path + file + '/match_' + file + '.png')
    plt.close()
        
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
    plt.savefig(my_path + file + '/' + file + '_'+ x + '_' + y + '.png')
    plt.close()

def match_two_tables(gaia, apt, file):
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
    df['del_ra'] = df.apply(lambda row: row.ra - row.CentroidRA, axis = 1)
    df['del_dec'] = df.apply(lambda row: row.dec - row.CentroidDec, axis = 1)
    df['del_mag'] = df.apply(lambda row: row.Magnitude - row.phot_bp_mean_mag, axis = 1)
    export_csv = df.to_csv(my_path + file + '/' + file + "_match.csv", index = None, header=True)
    plt.hist(df['diff'], bins = 20)
    plt.savefig(my_path + file + '/' +file+'_hist.png')
    plt.close()
    return df
