from pathlib import Path
from astropy.nddata import CCDData
from astropy.io import fits
from ccdproc import ImageFileCollection
import ccdproc as ccdp

data_directory = 'addheadsamp'


im_collection = ImageFileCollection(data_directory)

print(im_collection.summary['file', 'imagetyp', 'exptime'])