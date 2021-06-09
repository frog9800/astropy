import numpy as np

# Set up matplotlib
import matplotlib.pyplot as plt

from astropy.io import fits

hdu_list = fits.open("SampleFITS/201102hg-16.fit")
image_data = hdu_list[0].data

plt.imshow(image_data, cmap='gray')
plt.colorbar()
plt.show()
plt.savefig('samplefit.png')