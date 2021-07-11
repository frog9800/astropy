from astropy.io import fits

hdul = fits.open('Cal_science1/combined_science_10-12_600.000.fit')


hdul.info()
