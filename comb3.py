from pathlib import Path
import os
from glob import glob

from astropy.nddata import CCDData
from astropy.stats import mad_std

import ccdproc as ccdp
import matplotlib.pyplot as plt
import numpy as np

Sci_path1 = Path('Cal_science1')
Sci_path2 = Path('Cal_science2')
sci1_images = ccdp.ImageFileCollection(Sci_path1)
sci2_images = ccdp.ImageFileCollection(Sci_path2)

lights = firstreduced_images.summary['imagetyp'] == 'Light Frame'
light_times = set(firstreduced_images.summary['exptime'][lights])


for exp_time in sorted(light_times):
    calibrated_science = reduced_images.files_filtered(imagetyp='Light Frame', exptime=exp_time,
                                                     include_path=True)
    combined_science = ccdp.combine(calibrated_darks,
                                 method='average'
                                 )

    combined_science.meta['combined'] = True

    science_file_name = 'combined_science_10-12{:6.3f}.fit'.format(exp_time)
    combined_dark.write(calibrated_path / science_file_name)