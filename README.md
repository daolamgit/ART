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

The simplest command requires only a path to a set of DICOM image files and an output directory. This will use default parameters for the deformation.

```
python make_synthetic_image.py --input /path/to/dicom --output /my/output/directory
```

The `--input` flag should point to a directory containing one DICOM study (a set of dcm slices files for a single CT scan, for example).

The `--output` flag should point to a directory. If it doesn't exist, the script will create it. Within that directory, the script will create three files: `input.nii`, which is the input DICOM image in Nifti format, `bspline.txt`, which is the synthetic bspline transform file, and `warped.nii`, the warped image with the bspline transform applied to the input image. In addition, a directory `dicom` will be created within the output directory to hold the DICOM format of the warped image.

A more complicated command can be run to select various aspects of the deformation:

```
python make_synthetic_image.py --input /path/to/dicom --output /my/output/directory -d --maxdispl 20 --numspots 200 --gridspacing 30 30 20
```

The flag `-d` controls whether the deformation is diffeomorphic or not. If this flag is used, then maxdispl flag is ignored and random displacements are selected up to the grid spacing / 2.7. If not used, then synthetic deformations are randomly selected up to the maxdispl value.

The flag `--numspots` controls the number of control points with random deformations. `--gridspacing` controls the bspline grid spacing in x (LR), y (AP), z (SI) directions.
