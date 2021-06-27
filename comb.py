from pathlib import Path

from astropy.nddata import CCDData
from ccdproc import ImageFileCollection
import ccdproc as ccdp

ex1_path_raw = Path('sample_combine')

ex1_images_raw = ImageFileCollection(ex1_path_raw)

for ccd, file_name in ex1_images_raw.ccds(imagetyp='Dark Frame',  # Just get the dark frames
                                         ccd_kwargs={'unit': 'adu'},  # CCDData requires a unit for the image if
                                         # it is not in the header
                                         return_fname=True  # Provide the file name too.
                                         ):
    # Save the result
    modified= Path('addheadsamp')
    ccd.write(modified / file_name)

