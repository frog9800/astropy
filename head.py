from astropy.io import fits

hdul = fits.open('Cal_science1/output.fit')


hdul.info()
