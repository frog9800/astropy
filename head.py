from astropy.io import fits

#hdul = fits.open('addheadsamp/combined_dark_600.000.fit')

del hdul[1].data
del hdul[1].header
del hdul[2].data
del hdul[2].header

hdul.info()
