from ccdproc import ImageFileCollection
data_directory = 'SampleFITS/'
im_collection = ImageFileCollection(data_directory)
print(im_collection.summary['file', 'imagetyp', 'exptime'])