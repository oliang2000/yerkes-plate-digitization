## Yerkes Plate Digitization

:stars: **convert.py** (unfinished):  
Convert a plate image to fits file.

:stars: **extract.py** (unfinished):  
For a giving .fits file of a plate, generate a .csv file of star data extracted from plate, including:   
- location (RA, Dec)
- magnitude

:stars: **match.py & utility.py**:  
For a given .csv file (of star data extracted from plate), fetches star data of matching field from GAIA and compares data.  
- Saves both .csv files of the fetched GAIA data and matching data
- Outputs .png and .txt of the comparison of GAIA and plate data, accuracy calculated from treating GAIA as "truth:
  - photometry: deviation in mag (comparing to that of GAIA)
  - astrometry: deviation in RA and Dec (in both degrees and pixels)
  - Calculation methods:  *unfinished*

To use:  
1. Create a folder inside the same folder containing match.py & utility.py, this folder should be named the plate name, which will be used in created .csv files and figures in this program
2. The folder should contain .csv file with extraction data from the plate containing information about stars' RA, Dec, and magnitude. The default settings reads table from APT
3. Run *python match.py* in terminal, inside the folder containing match.py & utility.py
4. Proceed according to terminal instruction, figures will be produced in the folder with .csv file

**Packages**  


**References**  
[1] https://www.astromatic.net/software/sextractor   
[2] wwwwww
