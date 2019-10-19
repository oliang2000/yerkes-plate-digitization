from astroquery.gaia import Gaia
import matplotlib.pyplot as plt

def get_gaia_data (min_ra, max_ra, min_dec, max_dec, brightness):
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
    str((min_ra + max_ra)/2) + "," + str(min_dec+max_dec/2) + "," + \
    str(max_ra - min_ra) + "," + str(max_dec - min_dec) + "))=1 AND \
    (gaiadr2.gaia_source.phot_bp_mean_mag<=" + str(brightness) + ")" 
    job = Gaia.launch_job_async(query)
    r = job.get_results() 
    return r.to_pandas()

def run_next_step(question):
    '''
    Ask whether user would like to proceed to next step.
    '''
    inpt = question + " (y/n)?"
    g = input(inpt)
    if g == "y":
        return True
    elif g == "n":
        return False
    else:
        print("Please enter y/n.")
        return run_next_step(question)
        
def p_scatter(df1, df2, x, y, xlim =[0,0] , ylim=[0,0], lr1 = False, lr2 = False):
    '''
    This function uses matplotlib.pyplot to graph two functions and fix axis (calculated if not given).
    df2 should be the smaller subset.
    '''
    axes = plt.gca()
    plt.scatter(df1[x], df1[y])
    plt.scatter(df2[x], df2[y])
    if xlim == [0,0]:
        xlim = [df[x].min()-(df[x].max()-df[x].min())*0.05, df[x].max() + (df[x].max()-df[x].min())*0.05]
    if ylim == [0,0]:
        ylim = [df[y].min()-(df[y].max()-df[y].min())*0.05, df[y].max() + (df[y].max()-df[y].min())*0.05]
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

def match_two_tables(df1, df2):
    df = pd.DataFrame(columns=['ra','dec','phot_bp_mean_mag','pmra','pmdec',\
                               'CentroidRA','CentroidDec','Magnitude','diff'])
    dist = scipy.spatial.distance.cdist(df_apt[['CentroidRA', 'CentroidDec']], df_gaia[['ra', 'dec']])
    min_dist = np.argmin(dist, axis=1)

    m = 0
    while m < len(df_apt):
        n = min_dist[m]
        df = df.append({'ra': df_gaia['ra'][n], 'dec': df_gaia['dec'][n],\
                        'phot_bp_mean_mag': df_gaia['phot_bp_mean_mag'][n],\
                        'pmra': df_gaia['pmra'][n], 'pmdec': df_gaia['pmdec'][n],\
                        'CentroidRA': df_apt['CentroidRA'][m],'CentroidDec': df_apt['CentroidDec'][m],\
                        'Magnitude': df_apt['Magnitude'][m],\
                        'diff': dist[m][n]
                       }, ignore_index=True)
        m+=1
    sns.distplot(df['diff']).set_title('Distribution of difference in match')
