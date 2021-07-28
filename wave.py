from astropy.io import fits
import numpy as np
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
from astropy import stats
from astropy.io import ascii
from astropy.table import Table

with fits.open('Cal_science1/output.fit') as hdul:  # open a FITS file
    data = hdul[0].data  # assume the primary hdu is an image

# get values of the subsection from x=1 to 1536, y=499 to 500 (inclusive)
hdp = data[498:499, 0:1535]
print(hdp)
pp = hdp.sum(axis=0)
print(np.shape(pp))

num = np.arange(start=1, stop=1536, step=1)
inf = Table()
inf['x'] = np.array(num, dtype=np.float64)
inf['y'] = np.array(pp, dtype=np.float64)
ascii.write(inf, 'values.txt', overwrite=True)



