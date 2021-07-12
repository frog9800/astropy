from pathlib import Path
import shutil
from astropy.nddata import CCDData
from ccdproc import ImageFileCollection
import ccdproc as ccdp
from astropy import units as u
import numpy as np
from astropy.io import fits

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

calibrated_path = Path('addheadsamp')
reduced_images = ccdp.ImageFileCollection(calibrated_path)

darks = reduced_images.summary['imagetyp'] == 'Dark Frame'
dark_times = set(reduced_images.summary['exptime'][darks])


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

data, header = fits.getdata("Cal_science1/combined_science_10-12_600.000.fit", header=True)
hdu_number = 0 # HDU means header data unit
header['CALSTAT'] = "BD"
header['TELESCOP'] = "Barber20"
header['INSTRUME'] = "10-C Spectrograph"
header['OBJECT'] = "59 Cyg (SAO 50335) 01:50, 02:01, 02:12"
header['GRATING'] = "1200"
header['GRATING'] = "26.0 degrees"

fits.writeto('Cal_science1/output.fit', data, header, overwrite=True)