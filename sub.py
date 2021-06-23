import matplotlib.pyplot as plt
import numpy as np

from astropy.modeling import models
from astropy import units as u
from astropy import nddata
from astropy.io import fits

import ccdproc
import bottleneck as bn
from msumastro import ImageFileCollection, TableTree




imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std()) #...which is useful for determining scale when displaying an image.

data_dir = 'SampleFITS'

images = ImageFileCollection(data_dir, keywords='*')


def sigma_clip_median(image_list, clip_baseline_func=med_over_images):
    """
    Combine a list of images using median

    This function does several steps:


    1. sigma clip image using a median of the unclipped stack as the baseline
    2. combine the images on the list using median

    ** It modifies the images in the input list. **
    """

    combo = ccdproc.Combiner(image_list)
    combo.sigma_clipping(func=clip_baseline_func)
    return combo

exposures = [15, 30]
master_darks = {}
combiners = {}
for exposure in exposures:
    # make list of darks with this exposure time
    a_list = []
    for dark, fname in images.hdus(imagetyp='dark', exptime=exposure, return_fname=True):
        meta = dark.header
        meta['filename'] = fname
        a_list.append(ccdproc.CCDData(data=dark.data, meta=meta, unit="adu"))

        # get the exposure time as it appears in the fits file for use as a dictionary key
        exposure_time_in_fits_file = a_list[0].header['exptime']

        # make a combiner for sigma clipping and median combine
        a_combiner = overscan_trim_and_sigma_clip_median(a_list)
        combiners[exposure_time_in_fits_file] = a_combiner
        master_darks[exposure_time_in_fits_file] = a_combiner.median_combine(median_func=bn_median)

        # set the exposure time in the master -- right now combiner doesn't attempt to combine meta
        master_darks[exposure_time_in_fits_file].header['exptime'] = exposure_time_in_fits_file
        print
        "For exposure {} seconds there are {} bad pixels in the master.".format(exposure_time_in_fits_file,
                                                                                master_darks[
                                                                                    exposure_time_in_fits_file].mask.sum())
d_min, d_max, d_mean, d_std = imstats(np.asarray(master_darks[30.0]))
plt.figure(figsize=(15, 15))
plt.imshow(master_darks[30.0], vmax=d_mean + 4*d_std, vmin=d_mean - 4*d_std)
plt.show()