from astropy.io import fits

hdul = fits.open('addheadsamp/combined_dark_600.000.fit')


hdul.info()
