#extract tables
import pandas as pd
cal = pd.read_csv("starmatch/R3170/R3170_cal.csv", skiprows=1)#calibration stars from SDSS
qso = pd.read_csv("starmatch/R3170/R3170_QSO.csv", skiprows=1)#QSOs from SDSS
apt = pd.read_csv("starmatch/R3170/R3170.csv")#data from APT


#match tables
df = pd.DataFrame(columns=['ra', 'dec', 'u', 'g','r','i', 'z', 'CentroidRA','CentroidDec','Magnitude','diff'])
dist = scipy.spatial.distance.cdist(apt[['CentroidRA', 'CentroidDec']], cal[['ra', 'dec']])
min_dist = np.argmin(dist, axis=1)

m = 0
while m < len(apt):
    n = min_dist[m]
    df = df.append({'ra': cal['ra'][n], 'dec': cal['dec'][n], 'u': cal['u'][n], 'g': cal['g'][n], \
                    'r': cal['r'][n], 'i': cal['i'][n], 'z': cal['z'][n], 
                    
                    'CentroidRA': apt['CentroidRA'][m],'CentroidDec': apt['CentroidDec'][m],\
                    'Magnitude': apt['Magnitude'][m],\
                    'diff': dist[m][n]
                   }, ignore_index=True)
    m += 1
df = df.sort_values('diff').drop_duplicates('ra')
df = df.sort_values('diff').drop_duplicates('CentroidRA')
df = df.reset_index(drop=True) 

df['QSO'] = (df.ra.isin(qso.ra) & df.dec.isin(qso.dec))


#identify QSOs
#measure QSOs and save


