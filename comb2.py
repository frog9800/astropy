from pathlib import Path
import os
import shutil

from astropy.nddata import CCDData
from astropy.stats import mad_std
from astropy.io import fits
import ccdproc as ccdp
import matplotlib.pyplot as plt
import numpy as np

calibrated_path = Path('addheadsamp')
reduced_images = ccdp.ImageFileCollection(calibrated_path)

darks = reduced_images.summary['imagetyp'] == 'Dark Frame'
dark_times = set(reduced_images.summary['exptime'][darks])
print(dark_times)

for exp_time in sorted(dark_times):
    calibrated_darks = reduced_images.files_filtered(imagetyp='Dark Frame', exptime=exp_time,
                                                     include_path=True)
    combined_dark = ccdp.combine(calibrated_darks,
                                 method='median'
                                 )

    combined_dark.meta['combined'] = True
    dark_file_name = 'new_combined_dark_{:6.3f}.fit'.format(exp_time)
    ccdp.fits_ccddata_writer(combined_dark, dark_file_name, hdu_mask=None, hdu_uncertainty = None
    )

shutil.move(dark_file_name, calibrated_path) # moving a file to the original directory

