import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import units as u
import ccdproc as ccdp
from astropy.nddata import CCDData, StdDevUncertainty
from specutils import Spectrum1D
from ccdproc import trim_image, Combiner





from specreduce.extract import BoxcarExtract

sci = ccdp.ImageFileCollection('Cal_science1') # spectrum of the one target star, Gl 669A
sciimg = CCDData.read(sci, unit=u.adu)
BoxcarExtract(sci)
