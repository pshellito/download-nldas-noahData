# download-nldas-noahData

## Overview
The Matlab function, 'getNldasNoahData.m' downloads [NLDAS](http://ldas.gsfc.nasa.gov/nldas/) Noah data (soil moisture and other values) at any number of point locations and saves them as individual text files that contain the corresponding soil moisture timeseries. If a requested location is not in the NLDAS domain, the script will issue a warning and download the closest NLDAS pixel to the requested location. Soil moisture units are in kg/m2. To convert to volumetric (m3/m3), divide by 1000 (density of water), and then by the depth of the soil layer in meters.

## Requirements

You must first install [nctoolbox](https://github.com/nctoolbox/nctoolbox). Follow the instructions on that page to install.

Matlab R2008a+. You can verify the version of Matlab by typing:

      version

Register with URS and create a file titled .netrc in your home directory. Place the following in that file: (for details see Procedure 2 on the [GES DISC help page](https://disc.sci.gsfc.nasa.gov/recipes/?q=recipes/How-to-Download-Data-Files-from-HTTP-Service-with-wget))

      machine urs.earthdata.nasa.gov login <uid> password <password>

Java version 7 or higher. You can verify the version of Java used by Matlab by typing:

      version('-java')

The version returned should start with 'Java 1.7.' If it doesn't, you can try updating the Matlab JVM: http://www.mathworks.com/support/solutions/en/data/1-1812J/

## Demo

The script, 'callGetNldasNoahData.m' will
* Read the site name, latitude, and longitude of two arbitrary locations found in the input file, 'inFile.txt.'
* Define the range of dates over which to download. This is initially set to only 2 days for testing, which should take about a minute to run.
* Place output in a directory it creates called 'nldasNoahData.' Each location's simulation data will be in a space-delimited .txt file named after the site name. Details about the site's location are provided in a header in each file.

The above options can be changed in the 'callGetNldasNoahData.m' script.

Using a CU internet connection and a 2015 Macbook Pro, it takes about 4 hours to download one year of data.

## Contact
If you have any questions or edits, please [contact me](mailto:peter.shellito@colorado.edu) or create a new pull request.
