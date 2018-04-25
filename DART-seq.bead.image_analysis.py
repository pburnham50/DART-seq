###
#  Authors: P. Burnham, I. De Vlaminck (2018) 
#  
#  File: DART-seq.bead.image_analysis.py
#
#  Purpose: Produce intensity measurements from DART-seq beads tagged 
#  with fluorescent hybridization oligos. Images should be in .tif 
#  format.
###

### Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import skimage
import sys
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter
from skimage import data, io, filters, measure
from skimage import img_as_float
from skimage.morphology import reconstruction
from skimage.feature import register_translation
from skimage.draw import ellipse
from skimage.measure import label, regionprops
from skimage.transform import rotate
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb
from skimage.morphology import watershed
from skimage.feature import peak_local_max

### Settings/paths/
pathin = "image_files"
pathout = "bead_stats/"
channel = "cy5"

# Read in two images. Image 1 should be the bright field 
# while Image 2 is the fluorescence image.

dil = sys.argv[1] # indicates dilution of toeholds
rep = sys.argv[2] # indicates image replicate number

imag1 = skimage.io.imread(pathin + 'bf' + '.' + dil + '.' + rep + '.tif')  # bright-field
imag2 = skimage.io.imread(pathin + channel + '.' + dil + '.' + rep + '.tif')  # cy5 channel

# Shift image 1 with respect to imag2
shift_vectors = skimage.feature.register_translation(255-imag1, imag2)
tform = skimage.transform.SimilarityTransform(translation=shift_vectors[0])
warped = skimage.transform.warp(imag1, tform)

thresh1 = threshold_otsu(warped)
bw1 = closing(warped < thresh1, square(5))

# remove artifacts connected to image border
cleared1 = clear_border(bw1)

# label image regions
label_image1 = skimage.measure.label(cleared1)

# use watershed algorithm for segmentation of the beads
distance1 = ndi.distance_transform_edt(label_image1)
local_maxi1 = peak_local_max(distance1, indices=False, footprint=np.ones((3, 3)),
                            labels=cleared1)
markers1 = ndi.label(local_maxi1)[0]
labels1 = watershed(-distance1, markers1, mask=cleared1)

## use the image1 mask and the image 2 intensity to get the 
## properties of the beads

abc1 =regionprops(labels1, intensity_image=imag2)

# create an array with various statistics
## use intensity to make sure not counting empty patches
emp = 1 - 1*skimage.util.img_as_bool(imag2)

arr = np.empty((0,2), int)

for i in range(0,len(abc1)):
    if ((abc1[i].area > 5000) & (abc1[i].mean_intensity > np.mean(emp*imag2))) :
        arr = np.append(arr, np.array([[abc1[i].mean_intensity , abc1[i].area]]), axis=0)

final_array = np.column_stack((arr,np.repeat(dil,arr.shape[0]),np.repeat(rep,arr.shape[0])))

exp_file_name = pathout + channel + '.' + dil + '.' + rep +".bead_details.csv"
np.savetxt(exp_file_name, final_array.astype(float), delimiter=",")
