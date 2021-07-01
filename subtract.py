from pathlib import Path
from matplotlib import pyplot as plt
import ccdproc as ccdp
from astropy import units as u
import numpy as np
from astropy.nddata import CCDData
from astropy.io import fits

reduced_path = Path('addheadsamp')
ifc_reduced = ccdp.ImageFileCollection(reduced_path)

science_imagetyp = 'Light Frame'
dark_imagetyp = 'Dark Frame'
exposure = 'exptime'

pathraw = Path('sample_combine')
ifc_raw = ccdp.ImageFileCollection(pathraw)

combo_calibs = ifc_reduced.summary[ifc_reduced.summary['combined'].filled(False).astype('bool')]

combined_darks = {ccd.header[exposure]: ccd for ccd in ifc_reduced.ccds(imagetyp='Dark Frame', combined=True)}


def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=0.5):
    """
    Find the nearest exposure time of a dark frame to the exposure time of the image,
    raising an error if the difference in exposure time is more than tolerance.

    Parameters
    ----------

    image : astropy.nddata.CCDData
        Image for which a matching dark is needed.

    dark_exposure_times : list
        Exposure times for which there are darks.

    tolerance : float or ``None``, optional
        Maximum difference, in seconds, between the image and the closest dark. Set
        to ``None`` to skip the tolerance test.

    Returns
    -------

    float
        Closest dark exposure time to the image.
    """

    dark_exposures = np.array(list(dark_exposure_times))
    idx = np.argmin(np.abs(dark_exposures - image.header['exptime']))
    closest_dark_exposure = dark_exposures[idx]

    if (tolerance is not None and
            np.abs(image.header['exptime'] - closest_dark_exposure) > tolerance):
        raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                           'time {}.'.format(closest_dark_exposure, a_flat.header['exptime']))

    return closest_dark_exposure

# These two lists are created so that we have copies of the raw and calibrated images
all_reds = []
light_ccds = []
for light, file_name in ifc_raw.ccds(imagetyp=science_imagetyp, return_fname=True, ccd_kwargs=dict(unit='adu')):
    light_ccds.append(light)

    closest_dark = find_nearest_dark_exposure(light, combined_darks.keys())
    reduced = ccdp.subtract_dark(light, combined_darks[closest_dark],
                                 exposure_time=exposure, exposure_unit=u.second
                                 )
    all_reds.append(reduced)
    calibrated_name = 'calibrated_' + file_name
    reduced.write(reduced_path / calibrated_name )

