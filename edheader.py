from astropy.io import fits
data, header = fits.getdata("Cal_science1/combined_science_10-12_600.000.fit", header=True)
hdu_number = 0 # HDU means header data unit
header['CALSTAT'] = "BD"
header['TELESCOP'] = "Barber20"
header['INSTRUME'] = "10-C Spectrograph"
header['OBJECT'] = "59 Cyg (SAO 50335) 01:50, 02:01, 02:12"
header['GRATING'] = "1200"
header['GRATING'] = "26.0 degrees"
header['DATE'] = ""

fits.writeto('Cal_science1/output.fit', data, header, overwrite=True)