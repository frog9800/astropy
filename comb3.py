from pathlib import Path
import os
from glob import glob

from astropy.nddata import CCDData
from astropy.stats import mad_std

import ccdproc as ccdp
import matplotlib.pyplot as plt
import numpy as np

files = 'addheadsamp/calibrated201103hr8047-10.fit','addheadsamp/calibrated201103hr8047-11.fit','aaddheadsamp/calibrated201103hr8047-12.fit'
#calibrated_path = Path('addheadsamp')
firstreduced_images = ccdp.ImageFileCollection(files)

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