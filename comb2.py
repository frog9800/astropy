from pathlib import Path
import os

from astropy.nddata import CCDData
from astropy.stats import mad_std

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

    dark_file_name = 'combined_dark_{:6.3f}.fit'.format(exp_time)
    combined_dark.write(calibrated_path / dark_file_name)