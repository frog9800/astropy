from astropy.io import fits
import numpy as np
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt

from astropy import stats
from astropy.io import ascii
from astropy.table import Table

sigma1 = 3
sigma2 = 4
sigma3 = 5
sigma4 = 6
sigma5 = 7

with fits.open('Cal_science1/output.fit') as hdul1:  # open a FITS file
    data1 = hdul1[0].data  # assume that primary hdu is an image


sky1 = data1[0:50, 0:1535] #edge of background1


avg1 = np.mean(sky1) # background average1
cut1 = np.mean(data1, axis=1)

stv1 = np.std(sky1)*sigma1
stv2 = np.std(sky1)*sigma2
stv3 = np.std(sky1)*sigma3
stv4 = np.std(sky1)*sigma4
stv5 = np.std(sky1)*sigma5

sample1 = avg1 + stv1
sample2 = avg1 + stv2
sample3 = avg1 + stv3
sample4 = avg1 + stv4
sample5 = avg1 + stv5




detect1 = np.argwhere(sample1 < cut1) # coordinate(order of rows) of the spectrum1
detect2 = np.argwhere(sample2 < cut1)
detect3 = np.argwhere(sample3 < cut1)
detect4 = np.argwhere(sample4 < cut1)
detect5 = np.argwhere(sample5 < cut1)
print(np.size(detect2, 0))

#spectrum coordinate from top and bottom
spec1 = data1[np.int_((np.amin(detect1))):np.int_((np.amax(detect1)))]
spec2 = data1[np.int_((np.amin(detect2))):np.int_((np.amax(detect2)))]
spec3 = data1[np.int_((np.amin(detect3))):np.int_((np.amax(detect3)))]
spec4 = data1[np.int_((np.amin(detect4))):np.int_((np.amax(detect4)))]
spec5 = data1[np.int_((np.amin(detect5))):np.int_((np.amax(detect5)))]
extract1 = spec1.mean(axis=0)
extract2 = spec2.mean(axis=0)
extract3 = spec3.mean(axis=0)
extract4 = spec4.mean(axis=0)
extract5 = spec5.mean(axis=0)

skybot = np.amax(detect2) + 50
skytop = skybot + np.size(spec2, 0) + 1

print(np.size(spec2, 0))

sky2 = data1[np.int_(skybot):np.int_(skytop)]
print(skybot)
print(skytop)
print(np.size(sky2, 0))
print(np.mean(sky2))

skysam = sky2.sum(axis=0)

skyspec2 = np.subtract(spec2.mean(axis=0), skysam)


x = np.arange(start=1, stop=1537, step=1)
y1 = spec2.mean(axis=0)
y2 = skyspec2
y3 = spec3.mean(axis=0)
y4 = spec4.mean(axis=0)
y5 = spec5.mean(axis=0)



plt.title("Extractions of combined_science_10-12_600.000.fit")
plt.xlabel("Column")
plt.ylabel("Average adu")
plt.plot(x, y1, label = "Sigma =" + str(sigma2) + "(raw)")
plt.plot(x, y2, label = "Sigma =" + str(sigma2) + "(sky subtracted)")
#plt.plot(x, y3, label = "Sigma =" + str(sigma3))
#plt.plot(x, y4, label = "Sigma =" + str(sigma4))
#plt.plot(x, y5, label = "Sigma =" + str(sigma5))
plt.legend(bbox_to_anchor=(1.05, 1), loc='best', borderaxespad=0.)
#plt.yticks(np.arange(start=0.5, stop=2, step=0.25))
plt.show()

