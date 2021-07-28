from astropy.io import fits
import numpy as np
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
from astropy import stats
from astropy.io import ascii
from astropy.table import Table

with fits.open('Cal_science1/output.fit') as hdul:  # open a FITS file
    data = hdul[0].data  # assume that primary hdu is an image

sky = data[0, 0:1535]
avg = np.mean(sky)
stv = np.std(sky)*3
sample = avg + stv
print(sample)

spec = sample < data

line = data[np.all(spec, axis=1), :]

print(np.size(line,0))

leng= np.size(line,0)

mid1= leng/2
mid2= mid1 - 1

print(np.int_(mid1))
print(np.int_(mid2))

consec = line[[np.int_(mid2),np.int_(mid1)]]

print(consec)

adding = consec.sum(axis=0)

print(np.shape(adding))

num = np.arange(start=1, stop=1537, step=1)
inf = Table()
inf['x'] = np.array(num, dtype=np.float64)
inf['y'] = np.array(adding, dtype=np.float64)
ascii.write(inf, 'Cal_science1/algorithm.txt', overwrite=True)

#print(np.random.normal(song, 3, size=(2, 1536)))

