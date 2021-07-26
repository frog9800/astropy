from astropy.io import fits
import numpy as np
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
from astropy import stats

with fits.open('Cal_science1/output.fit') as hdul:  # open a FITS file
    data = hdul[0].data  # assume the first extension is an image
print(data[487:510, 0:1535])   # get the pixel value at x=5, y=2

# get values of the subsection from x=11 to 20, y=31 to 40 (inclusive)
#print(stats.sigma_clip(data[487:510, 0:1535], sigma=2, maxiters=5))
hdul1= data[487:510, 0:1535]

hdul1.writeto('Cal_science1/spec.fits', overwrite=True)
#y, x = np.mgrid[0:1023, 0:1535]
#plt.imshow(data[487:510, 0:1535])

#plt.show()
