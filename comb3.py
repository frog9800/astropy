from pathlib import Path
import os
from glob import glob
import shutil

from astropy.nddata import CCDData
from astropy.stats import mad_std

import ccdproc as ccdp
import matplotlib.pyplot as plt
import numpy as np

Sci_path1 = Path('Cal_science1')
Sci_path2 = Path('Cal_science2')
sci1_images = ccdp.ImageFileCollection(Sci_path1)
sci2_images = ccdp.ImageFileCollection(Sci_path2)

lights = sci1_images.summary['imagetyp'] == 'Light Frame'
light_times = set(sci1_images.summary['exptime'][lights])

lights_two = sci2_images.summary['imagetyp'] == 'Light Frame'
#light_times_two = set(sci2_images.summary['exptime'][lights])


for exp_time in sorted(light_times):
    calibrated_science = sci1_images.files_filtered(imagetyp='Light Frame', exptime=exp_time,
                                                     include_path=True)
    combined_science = ccdp.combine(calibrated_science,
                                 method='average'
                                 )

    combined_science.meta['combined'] = True

    science_file_name = 'combined_science_10-12_{:6.3f}.fit'.format(exp_time)
    ccdp.fits_ccddata_writer(combined_science, science_file_name, hdu_mask=None, hdu_uncertainty=None)
    shutil.move(science_file_name, Sci_path1)

