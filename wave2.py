from astropy.io import fits
import numpy as np
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
from astropy import stats
from astropy.io import ascii
from astropy.table import Table

with fits.open('Cal_science1/output.fit') as hdul:  # open a FITS file
    data = hdul[0].data  # assume that primary hdu is an image

sky = data[0:50, 0:1535] #edge of background
avg = np.mean(sky) # background average
cut = np.mean(data, axis=1)
print(np.size(cut, 0))
print(cut)
print(avg)
print(np.std(sky))
stv = np.std(sky)*2.5 # 3x standard deviation
sample = avg + stv
print(sample)


detect = np.argwhere(sample < cut) # coordinate(order of rows) of the spectrum

min = (np.amin(detect)) # bottom of the spectrum
max = (np.amax(detect)) # top of the spectrum
print(min)
print(max)


#spec = sample < data # spectrum condition. spectrum must be more than average + 3x standard deviation of backgrouind pixel

#line = data[np.all(spec, axis=1), :] # slicing arrays by given condition.

#print(np.size(line, 0))

#leng= np.size(line, 0) # size of sliced arrays = total number of rows of the spectrum

# Assign two consecutive rows.  half of the total number of spectrum rows = midpoint.
#mid1 = leng/2
#mid2 = mid1 - 1

#make sure to state them as integer.


#spectrum coordinate from top and bottom
spec = data[np.int_(min):np.int_(max)]

print(np.size(spec, 0))


#Averaging out all spectrum rows.
extract = spec.mean(axis=0)

print(np.shape(extract))

# Writing to 1D extraction file

num = np.arange(start=1, stop=1537, step=1)
inf = Table()
inf['x'] = np.array(num, dtype=np.float64)
inf['y'] = np.array(extract, dtype=np.float64)
ascii.write(inf, 'Cal_science1/algorithm.txt', overwrite=True)



