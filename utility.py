from astroquery.gaia import Gaia

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
