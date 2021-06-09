import os
from glob import glob
import pyfits
import shutil
from matplotlib import pyplot as plt

files = glob('SampleFITS/*.fit')

for file_name in files:
    data = pyfits.getdata(file_name)
    name = file_name[len('SampleFITS/'):-len('.fit')] # Remove directory name and .fits from the file name

    plt.imshow(data, cmap='gray')
    plt.colorbar()
    plt.title(name)
    plt.savefig('result_images/' + name + '.png')
    plt.show()

