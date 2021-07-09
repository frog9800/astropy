from astropy.io import fits

hdul = fits.open('addheadsamp/new_combined_dark_600.000.fit')


hdul.info()
