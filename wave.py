import sys
import numpy as np
from astropy.io import fits
from aspired import image_reduction
from aspired import spectral_reduction
import plotly.io as pio

science = 'Cal_science1/output.fit'
science_frame = fits.open(science)



sci = spectral_reduction.TwoDSpec(
    science_frame,
    readnoise=0,
    cosmicray=False,
    gain=1,
    seeing=1,
    verbose=False
)

sci.ap_trace()
sci.ap_extract(
    apwidth=20,
    skysep=3,
    skywidth=5,
    skydeg=0,
    optimal=True,
    display=True,
    renderer='jpg',
    filename='sciextract')

