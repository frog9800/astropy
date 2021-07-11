from pathlib import Path
from astropy.nddata import CCDData
from astropy.io import fits
from ccdproc import ImageFileCollection
import ccdproc as ccdp

data_directory = 'Cal_science1'
data, header = fits.getdata("Cal_science1/combined_science_10-12_600.000.fit", header=True)

im_collection = ImageFileCollection(data_directory)

print(im_collection.summary['object','file', 'imagetyp', 'exptime'])

hdu_number = 0 # HDU means header data unit
print(fits.getheader('Cal_science1/output.fit', hdu_number))