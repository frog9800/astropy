import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import units as u

from astropy.nddata import CCDData, StdDevUncertainty
from specutils import Spectrum1D
from ccdproc import trim_image, Combiner

# need to get import to work in notebook w/o global package install
import sys
sys.path.append('..')
import apextract as ap
import fluxcal as fc
import flatfield as ft

sci = 'Cal_science1/combined_science_10-12_600.000.fit'
sciimg = CCDData.read(sci, unit=u.adu)
# BASIC REDUCTION:
# subtract BIAS, divide FLAT, ExpTime, put in units of ADU/s


sciimg.unit = sciimg.unit / u.s

# trim off bias section
calimg = trim_image(calimg, fits_section=calimg.header['DATASEC'])
# Now remove FLAT from ilum section
calimg.data[ilum,:] = calimg.data[ilum,:] / FLAT

plt.figure(figsize=(6,3))
plt.imshow(sciimg, origin='lower', aspect='auto')
plt.clim(np.percentile(sciimg, (5, 98)))


# Do the trace for both the Science image and the Flux Calibration images
sci_tr = ap.trace(sciimg, display=True, nbins=25)


# Extraction of the spectrum along the trace for both the Science and Flux Cal images
sci_ex = ap.extract(sciimg, sci_tr, display=True, apwidth=10, skysep=5, skywidth=5)
