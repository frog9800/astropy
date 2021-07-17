from astropy.io import fits
from astropy import units as u
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
quantity_support()  # for getting units on the axes below

f = fits.open('Cal_science1/combined_science_10-12_600.000.fit')
# The spectrum is in the second HDU of this file.
specdata = f[0].data
f.close()

from specutils import Spectrum1D
lamb = 10**specdata['loglam'] * u.AA
flux = specdata['flux'] * 10**-17 * u.Unit('erg cm-2 s-1 AA-1')
spec = Spectrum1D(spectral_axis=lamb, flux=flux)

f, ax = plt.subplots()
ax.step(spec.spectral_axis, spec.flux)