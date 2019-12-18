# art-testing
testing data generation scripts for Varian adaptive project

Applies a synthetic deformation to an image to generate a synthetic image. This
script is capable of the following actions:
  1) deform a CT scan to produce a new CT scan
  2) deform a CT scan to produce a new CBCT scan
  3) deform a CBCT scan to produce a new CBCT scan

## Installation Instructions

### plastimatch

art-testing uses plastimatch for the heavy lifting. Please install plastimatch on your system first:

1. Download and install plastimatch by [building from source](http://plastimatch.org/building_plastimatch.html) or by [installing the Windows binary](http://plastimatch.org/windows_installation.html).
2. If you did not install plastimatch into a standard path, update the variable `plastimatch_executable` in `make_synthetic_image.py`.

### SimpleITK and numpy

art-testing needs SimpleITK and numpy for image manipulation.

SimpleITK can be installed with pip or using a package manager such as [Anaconda](https://www.anaconda.com/distribution/) (I recommend Anaconda). Installation instructions for SimpleITK can be found [here](https://simpleitk.readthedocs.io/en/master/Documentation/docs/source/installation.html).

Numpy can be installed similarly (by default with Anaconda or with pip).

That's it! You're now ready to use the script `make_synthetic_image.py`

## Usage Instructions

### Input and Output

Input CT images, CBCT images and structure file must follow the following structure:
<pre>
#   #   #   #   #   #   #   #   #   #   #
# Input:
#   Pt1
#       Origin (e.g Ptname)/
#         -CT
#             -CT.dcms
#         -RS
#             -RS.dcm
#         -CBCT
#             -CBCT.dcms
#
# Output:
#   Pt1
#       Origin/
#           -CT
#           -RS
#           -CBCT
#           -input.nii # a warped of CT folder
#       Deform1 (or Fx1)/
#           -gaussian
#               -dicom
#               -vf.mha
#               -warped.nii
#           -CT
#               CT.dcms
#               RS.dcm

#           -CBCT
#             -CBCT.dcms
#             -RS??? No need
#
#       Deform2 (or Fx2)/
#           -gaussian
#               -dicom
#               -vf.mha
#               -warped.nii
#           -CT
#               CT.dcms
#               RS.dcm

#           -CBCT
#             -CBCT.dcms
#             -RS??? No need
</pre>
The Gaussian information needed for the deformation can be entered in a .csv file. An example is given in Gaussian.csv.
Note about the Gaussian information:
- center: x y z in mm. It follows plastimatch coordinate
- magnitude: x y z in mm. It can be negative, which means the vector will point in the reverse direction
- standard deviation: in mm.
- nGaussian: the order of Gaussian components. There can be many Gaussian components.
```
