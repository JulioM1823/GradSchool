INFO ON STARS

The star data files are:
cd_34d241.txt
hr1544.txt
hr3454.txt
ltt1020.txt

Col l = the wavelength in angstroms
Col 2 = the flux density 1.e16*f_lambda 

The f_lambda are multiplied by the constant 1.e16.  You will want to
convert to actual f_lambda by dividing these values by 1.e16 Do not
worry about the other columns.

INFO ON VEGA

The spectrum is in the file:
alpha_lyr.txt

Col l = the wavelength in angstroms  
Col 2 = the flux density f_lambda 

The f_lambda are in proper flux density units. No conversion is required.

INFO ON FILTER DATA

The Johnson-Cousins filter response curves are in the files:
Vband-response.txt
Bband-response.txt

Col 1 = wavelength in angstroms 
Col 2 = relative transmission at wavelength lambda


