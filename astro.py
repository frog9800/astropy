import numpy as np

# Set up matplotlib
import matplotlib.pyplot as plt

from astropy.io import fits

image_file = "Cal_science1/combined_science_10-12_600.000.fit"
hdu_list = fits.open(image_file)
image_data = fits.getdata(image_file)

plt.imshow(image_data, cmap='gray')
plt.colorbar()
#plt.savefig('samplefit.png')
plt.show()