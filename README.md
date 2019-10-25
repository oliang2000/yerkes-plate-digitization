## Yerkes Plate Digitization Data Analysis
:stars:

**match.py & utility.py**:  
For a given .csv file (of star data extracted from plate), fetches star data of matching field from GAIA and compares data. Saves both .csv files of the fetched GAIA data and matching data, and .png outputs of comparison.

How to use:  
1. Create a folder inside the same folder containing match.py & utility.py, this folder should be named the plate name, which will be used in created .csv files and figures in this program.
2. The folder should contain .csv file with extraction data from the plate containing information about stars' RA, Dec, and magnitude. The default settings reads table from APT.
3. Run python match.py inside the folder in terminal. 
