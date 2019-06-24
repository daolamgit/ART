#!/usr/bin/env python
from __future__ import print_function
"""
Applies a synthetic deformation to an image to generate a synthetic image. This
script is capable of the following actions:
  1) deform a CT scan to produce a new CT scan
  2) deform a CT scan to produce a new CBCT scan
  3) deform a CBCT scan to produce a new CBCT scan
"""

""" make_synthetic_image.py
    (c) Washington University, 2019
    Authors: G. Hugo
    Revision Log:
    0.1    initial build     6/1/2019
"""

from argparse import ArgumentParser
import sys
import hashlib, os
from os import makedirs
from os.path import join, exists
import subprocess
import numpy as np
import SimpleITK as sitk

def np_to_str(s, n):
    return s + ' = ' + ' '.join(map(str,n)) + '\n'

class BSplinePM(object):
    """Plastimatch bspline file.

    # Arguments
        origin: image origin, list or tuple of length 3
        spacing: image spacing, in mm, list or tuple of length 3
        dimension: image dimension, list or tuple of length 3
        vox_per_rgn: grid spacing, in voxels
        coeffs: 4D matrix of coefficients. Stored as [cd,sz,sy,sx], where sx is
            the spatial x direction, sy is spatial y, etc. and cd is the depth
            layer, which is the control point values in the x (d=0), y (d=1),
            and z (d=2) directions.
    """

    def __init__(self,
                 origin=None,
                 spacing=None,
                 dimension=None,
                 vox_per_rgn=None,
                 coeffs=None):
        super(BSplinePM, self).__init__()
        if origin is not None:
            self.origin = np.asarray(origin)
        else:
            self.origin = None
        if spacing is not None:
            self.spacing = np.asarray(spacing)
        else:
            self.spacing = None
        if dimension is not None:
            self.dimension = np.asarray(dimension)
            self.roi_dim = self.dimension
        else:
            self.dimension = None
            self.roi_dim = None
        if vox_per_rgn is not None:
            self.vox_per_rgn = np.asarray(vox_per_rgn)
        else:
            self.vox_per_rgn = None

        # set up internal items
        self.dc = np.ndarray.flatten(np.identity(3))
        self.roi_offset = np.zeros(3)


        self.grid_spacing = None
        self.number_of_regions = np.zeros(3)
        self.number_of_control_points = np.zeros(3)

        # self.coeffs[d,z,y,x]
        if coeffs is not None:
            if len(coeffs.shape) == 4:
                self._setup(coeffs)
            else:
                print('coeffs must be a 4D numpy array')
                return()

        # plastimatch bspline.txt format:
        # MGH_GPUIT_BSP <experimental>
        # img_origin = -249.511719 -457.511719 51.500000
        # img_spacing = 0.976562 0.976562 3.000000
        # img_dim = 512 512 130
        # roi_offset = 0 0 0
        # roi_dim = 512 512 130
        # vox_per_rgn = 30 30 10
        # direction_cosines = 1.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000
        # -0.00167561590205878019
        # ...
        # coeffs stored as [x1,x2,..xN, y1,y2,..yN,z1,z2,..zN]

    def _setup(self, coeffs=None):
        if self.vox_per_rgn is not None and self.spacing is not None:
            self.grid_spacing = np.multiply(self.vox_per_rgn, self.spacing)

        if self.roi_dim is not None and self.vox_per_rgn is not None:
            self.number_of_regions = 1 + (self.roi_dim - 1) / self.vox_per_rgn

        if self.number_of_regions.any() > 0:
            self.number_of_control_points = 3 + self.number_of_regions

        if coeffs is not None:
            self.coeffs = np.asarray(coeffs)
        else:
            self.coeffs = None


    def write_file(self, filename):
        pass
        with open(filename, "wb") as f:
            f.write('MGH_GPUIT_BSP <experimental>\n')

            f.write(np_to_str('img_origin', self.origin))
            f.write(np_to_str('img_spacing', self.spacing))
            f.write(np_to_str('img_dim', self.dimension.astype(int)))
            f.write(np_to_str('roi_offset', self.roi_offset.astype(int)))
            f.write(np_to_str('roi_dim', self.roi_dim.astype(int)))
            f.write(np_to_str('vox_per_rgn', self.vox_per_rgn.astype(int)))
            f.write(np_to_str('direction_cosines', self.dc))

            for d in range(self.coeffs.shape[0]):
                for z in range(self.coeffs.shape[1]):
                    for y in range(self.coeffs.shape[2]):
                        for x in range(self.coeffs.shape[3]):
                            f.write(self.coeffs[d,z,y,x].astype(str) + '\n')


    def read_file(self, filename):
        with open(filename) as f:
            # check first line
            line = f.readline()
            if not line.startswith('MGH_GPUIT_BSP'):
                print('file: ' + filename + ' is not a plastimatch bspline file.')
                return None

            while(True):
                line = f.readline()
                split_line = line.split('=')
                if len(split_line) == 2: # assignment
                    if split_line[0].strip() == 'img_origin':
                        self.origin = np.asarray([float(i) for i in split_line[1].split()])
                    elif split_line[0].strip() == 'img_spacing':
                        self.spacing = np.asarray([float(i) for i in split_line[1].split()])
                    elif split_line[0].strip() == 'img_dim':
                        self.dimension = np.asarray([int(i) for i in split_line[1].split()])
                    elif split_line[0].strip() == 'roi_offset':
                        self.roi_offset = np.asarray([float(i) for i in split_line[1].split()])
                    elif split_line[0].strip() == 'roi_dim':
                        self.roi_dim = np.asarray([int(i) for i in split_line[1].split()])
                    elif split_line[0].strip() == 'vox_per_rgn':
                        self.vox_per_rgn = np.asarray([int(i) for i in split_line[1].split()])
                    elif split_line[0].strip() == 'direction_cosines':
                        self.dc = np.asarray([float(i) for i in split_line[1].split()])
                    else:
                        print('unknown parameter: ' + split_line[0].strip())

                elif len(split_line) == 1: # should be start of coefficents
                    self._setup() # regenerate internals based on read parameters
                    self.coeffs = np.ndarray((3, self.number_of_control_points[2], self.number_of_control_points[1], self.number_of_control_points[0]))
                    self.coeffs[0,0,0,0] = float(split_line[0].strip())
                    for d in range(self.coeffs.shape[0]):
                        for z in range(self.coeffs.shape[1]):
                            for y in range(self.coeffs.shape[2]):
                                for x in range(self.coeffs.shape[3]):
                                    if x == 0 and y == 0 and z == 0 and d == 0:
                                        self.coeffs[0,0,0,0] = float(split_line[0].strip())
                                    else:
                                        line = f.readline()
                                        self.coeffs[d,z,y,x] = float(line.strip())
                    break
                else:
                    print('error parsing string')
                    return None


    def initialize_from_image(self, origin, spacing, dim, vox_per_rgn):
        self.__init__(origin, spacing, dim, vox_per_rgn)
        self._setup()


    def make_coeffs(self, max_displacement=1.0, number_of_spots=1, constrain=False):
        # select random displacement (single voxel) and location based on spacing
        # displacement can be no bigger than gridSize / 2.479 (Choi and Lee).
        if(constrain):
            K=2.479
        else:
            K=max_displacement

        displacement_cap = self.grid_spacing / K

        self.coeffs = np.zeros((3, self.number_of_control_points[2], self.number_of_control_points[1], self.number_of_control_points[0]))

        for spot in range(number_of_spots):
            rand_displ = (2*displacement_cap) * np.random.random(3) + -displacement_cap
            rand_idx = [0,0,0]
            for idx in range(3):
                rand_idx[idx] = np.random.randint(low=0, high=self.number_of_control_points[idx])
            for idx in range(3):
                self.coeffs[idx, rand_idx[2], rand_idx[1], rand_idx[0]] = rand_displ[idx]


class DeformImage(object):
    """Class to import and synthetically deform a DICOM image. Parses a
       directory to find DICOM files and reads these in. Based on arguments,
       either generates a CT (same volume) or CBCT (truncated volume) image.
       Images can be deformed either by a b-spline transform or by a displacement
       field (vector field image).

    # Arguments
        path: path to directory containing DICOM image files.
        constrain: constrain to a diffeomorphic deformation.
        number_of_spots: number of control points to displace.
        max_displacement: if not constrained, maximum displacement, in mm.
    """

    def __init__(self,
                 store_path,
                 read_path,
                 constrain=False,
                 number_of_spots=1,
                 max_displacement=1.0,
                 ):
        super(DeformImage, self).__init__()
        self.path = read_path           # location of DICOM image
        self.store_path = store_path    # output location
        self.constrain = constrain      # diffeomorphic transform?
        self.number_of_spots = number_of_spots      # number of control points to deform 
        self.max_displacement = max_displacement    # max displacement in mm
        self.plastimatch_executable = 'plastimatch' # UPDATE IF NECESSARY
        self.bspline = BSplinePM()
        self.origin = None              # image domain info for building bspline
        self.spacing = None             #
        self.dimension = None           #
        self.image = None               # image for thresholding location

    def read_image(self):
        # generate image with plastimatch
        call_params = [self.plastimatch_executable]
        call_params.append('convert')

        call_params.append('--input')
        call_params.append(self.path)

        call_params.append('--output-img')
        call_params.append(join(self.store_path,'input.nii'))

        subprocess.call(call_params)

        # read in image with SimpleITK
        img = sitk.ReadImage(join(self.store_path,'input.nii'))
        self.origin = np.asarray(img.GetOrigin())
        self.spacing = np.asarray(img.GetSpacing())
        self.dimension = np.asarray(img.GetSize())


    def make_bspline(self, vox_per_rgn=[30,30,30]):
        # make bspline
        if len(vox_per_rgn) == 1:
            vox_per_rgn = np.tile(vox_per_rgn, 3)
        if len(vox_per_rgn) != 3:
            print('vox_per_rgn must be a numpy array of length 3.')
            return
        else:
            vox_per_rgn = np.asarray(vox_per_rgn)

        self.bspline.initialize_from_image(self.origin, self.spacing, self.dimension, vox_per_rgn)
        self.bspline.make_coeffs(max_displacement=self.max_displacement, number_of_spots=self.number_of_spots, constrain=self.constrain)
        self.bspline.write_file(filename=join(self.store_path,'bspline.txt'))


    def deform_image(self, vox_per_rgn):
        # make bspline
        self.make_bspline(vox_per_rgn)

        # deform image
        call_params = [self.plastimatch_executable]
        call_params.append('warp')

        call_params.append('--input')
        call_params.append(join(self.store_path,'input.nii'))

        call_params.append('--output-img')
        call_params.append(join(self.store_path,'warped.nii'))

        call_params.append('--xf')
        call_params.append(join(self.store_path,'bspline.txt'))

        subprocess.call(call_params)


    def write_image(self):
        call_params = [self.plastimatch_executable]
        call_params.append('convert')

        call_params.append('--input')
        call_params.append(join(self.store_path,'warped.nii'))

        call_params.append('--output-dicom')
        call_params.append(join(self.store_path,'dicom'))

        subprocess.call(call_params)

if (__name__ == '__main__'):

    # argument parsing
    parser = ArgumentParser()
    parser.add_argument('-b', '--bspline',  help='bspline file')
    parser.add_argument('-i', '--input', help='path to DICOM input image to be deformed', required=True)
    parser.add_argument('-d', '--diffeomorphic', dest='diffeo', help='Constrain displacement to be diffeomorphic', action='store_true')
    parser.add_argument('-m', '--maxdispl', help='Maximum displacement, in image coordinates', type=float)
    parser.add_argument('-n', '--numspots', help='Number of gaussian spot displacements', type=int)
    parser.add_argument('-g', '--gridspacing', nargs='+', help='grid spacing in voxels', type=int)
    parser.add_argument('-o', '--output', help='output directory', required=True)
    args = parser.parse_args()

    # testing
    #bspline = BSplinePM()
    #bspline.read_file(args.bspline)
    #bspline.write_file(args.output)

    # create output directories
    if not exists(args.output):
        try:
            makedirs(args.output)
        except:
            print("Error:", sys.exc_info()[1])
            exit()
    if not exists(join(args.output, 'dicom')):
        try:
            makedirs(join(args.output, 'dicom'))
        except:
            print("Error:", sys.exc_info()[1])
            exit()

    deformer = DeformImage(read_path=args.input, store_path=args.output, constrain=args.diffeo, number_of_spots=args.numspots, max_displacement=args.maxdispl)
    deformer.read_image()
    deformer.deform_image(args.gridspacing)
    deformer.write_image()

    # plastimatch convert --input CBCT --output-img cbct.nii
    # plastimatch warp cbct.nii
