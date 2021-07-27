from astropy.io import fits
import numpy as np
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
from astropy import stats
from astropy.io import ascii
from astropy.table import Table

with fits.open('Cal_science1/output.fit') as hdul:  # open a FITS file
    data = hdul[0].data  # assume the first extension is an image
#print(data[488:513, 0:1535])   # get the pixel value at x=5, y=2

# get values of the subsection from x=11 to 20, y=31 to 40 (inclusive)
#print(stats.sigma_clip(data[487:510, 0:1535], sigma=2, maxiters=5))
hdp = data[493:504, 0:1535]
pp = hdp.sum(axis=0)

num = np.arange(start=1, stop=1536, step=1)
inf = Table()
inf['x'] = np.array(num, dtype=np.float64)
inf['y'] = np.array(pp, dtype=np.float64)
ascii.write(inf, 'values.txt', overwrite=True)

#plt.imshow(data[487:510, 0:1535])

#plt.show()
