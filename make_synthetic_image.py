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
    0.1    initial build     06/01/2019
    0.2    remove bones      08/13/2019
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

class RigidPM(object):
    """Plastimatch rigid 6 dof file.

    # Arguments
        rotation: rotation angles, in radians (1x3 tuple)
        translation: translation, in mm (1x3 tuple)
    """

    def __init__(self, rotation=(0,0,0), translation=(0,0,0), center_rotation=(0,0,0)):
        super(RigidPM, self).__init__()
        if len(rotation) != 3:
            self.rotation = np.asarray((0,0,0))
        else:
            self.rotation = np.asarray(rotation)
        if len(translation) != 3:
            self.translation = np.asarray((0,0,0))
        else:
            self.translation = np.asarray(translation)
        if len(center_rotation) != 3:
            self.center_rotation = np.asarray((0,0,0))
        else:
            self.center_rotation = np.asarray(center_rotation)


    def set_rotation(self, rotation):
        self.rotation = np.asarray(rotation)


    def set_translation(self, translation):
        self.translation = np.asarray(translation)


    def set_center_rotation(self, center_rotation):
        self.center_rotation = np.asarray(center_rotation)


    def write_file(self, filename):
        with open(filename, "wb") as f:
            f.write('#Insight Transform File V1.0\n')
            f.write('#Transform 0\n')
            f.write('Transform: VersorRigid3DTransform_double_3_3\n')

            params = np.concatenate((self.rotation, self.translation))
            param_str = 'Parameters: ' + ' '.join(map(str,params)) + '\n'
            f.write(param_str)

            fixed_params_str = 'FixedParameters: ' + ' '.join(map(str,self.center_rotation)) + '\n'
            f.write(fixed_params_str)


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
            self.number_of_regions = 1 + (self.roi_dim.astype(int) - 1) / self.vox_per_rgn.astype(int)

        if self.number_of_regions.any() > 0:
            self.number_of_control_points = 3 + self.number_of_regions.astype(int)

        if coeffs is not None:
            self.coeffs = np.asarray(coeffs)
        else:
            self.coeffs = None


    def write_file(self, filename):
        with open(filename, "w") as f:
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


    def _region_from_cp_index(self, idx):
        # returns a tuple of indices (start_x, start_y, start_z, end_x, end_y, end_z) defining the influence region for control point at index idx
        st = [0,0,0]
        en = [0,0,0]
        for i in range(3):
            # CP in image run from 0 to N-3, due to support needed outside of image
            # but really care about tiles 1 and 2 (out of (0,1,2,3)), so exclude first and last tiles
            st[i] = (idx[i] + 1) * self.vox_per_rgn[i]
            # Each CP covers 4 knot points, but only care about middle 2
            en[i] = (idx[i] + 3) * self.vox_per_rgn[i]
        return st + en


    def _histogram_from_image_region(self, region, image, low_threshold=300, high_threshold=700):
        # returns the histogram of the self.image within region

        # only need 3 bins, air, soft tissue, and bone
        bins = [-1050, low_threshold, high_threshold, 5000]

        # extract region
        image = image[region[2]:region[5], region[1]:region[4], region[0]:region[3]]

        # histogram of 3 bins
        h, bins = np.histogram(image, bins, density=False)

        return h


    def _is_good_grid_point_location(self, idx, image, low_threshold=300, high_threshold=700, percentage_threshold=90.0):
        # use a histogram to determine if the region influenced by the grid point contains enough soft tissue

        # compute histogram in grid point influence region
        region = self._region_from_cp_index(idx)
        h = self._histogram_from_image_region(region, image, low_threshold, high_threshold)
        if len(h) != 3:
            print('Error: histogram of wrong size: ', h)
            exit()
        return h[1] > percentage_threshold


    def make_coeffs(self, max_displacement=1.0, number_of_spots=1, constrain=False, remove_bones=False, image=None):
        # select random displacement (single voxel) and location based on spacing
        # displacement can be no bigger than gridSize / 2.479 (Choi and Lee).
        if(constrain):
            K=2.479
            displacement_cap = self.grid_spacing / K
        else:
            displacement_cap = max_displacement

        self.coeffs = np.zeros((3, self.number_of_control_points[2], self.number_of_control_points[1], self.number_of_control_points[0]))

        count = 0
        if image is not None:
            image = sitk.GetArrayFromImage(image)
        for spot in range(number_of_spots):
            rand_displ = (2*displacement_cap) * np.random.random(3) + -displacement_cap
            rand_idx = [0,0,0]

            if remove_bones:
                # limit grid point to locations with soft tissue (HU -700 to 700)
                found_soft_tissue_location = False
                while not found_soft_tissue_location:
                    for idx in range(3):
                        # CP in image run from 0 to N-3, due to support needed outside of image
                        rand_idx[idx] = np.random.randint(low=0, high=self.number_of_control_points[idx]-3)
                    count += 1
                    if self._is_good_grid_point_location(rand_idx, image):
                        found_soft_tissue_location = True
                        print('found good point: ', rand_idx, ' with displacement: ', rand_displ)
                    else:
                        print('point not within soft tissue, retrying...')
            else:
                # allow grid points anywhere
                for idx in range(3):
                    # CP in image run from 0 to N-3, due to support needed outside of image
                    rand_idx[idx] = np.random.randint(low=0, high=self.number_of_control_points[idx]-3)
                count += 1
            for idx in range(3):
                self.coeffs[idx, rand_idx[2], rand_idx[1], rand_idx[0]] = rand_displ[idx]

        print('Sampled ', count, ' points to obtain ', number_of_spots, ' grid points')


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
        remove_bones: try to constrain grid points to <300 HU
    """

    def __init__(self,
                 store_path,
                 read_path,
                 constrain=False,
                 number_of_spots=1,
                 max_displacement=1.0,
                 remove_bones=False
                 ):
        super(DeformImage, self).__init__()
        self.path = read_path           # location of DICOM image
        self.store_path = store_path    # output location
        self.constrain = constrain      # diffeomorphic transform?
        self.number_of_spots = number_of_spots      # number of control points to deform
        self.max_displacement = max_displacement    # max displacement in mm
        self.remove_bones = remove_bones   # no grid points in bones
        self.plastimatch_executable = 'plastimatch' # UPDATE IF NECESSARY
        self.bspline = BSplinePM()
        self.rigid = RigidPM()
        self.origin = None              # image domain info for building bspline
        self.spacing = None             # image domain info for building bspline
        self.dimension = None           # image domain info for building bspline
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
        self.image = sitk.ReadImage(join(self.store_path,'input.nii'))
        self.origin = np.asarray(self.image.GetOrigin())
        self.spacing = np.asarray(self.image.GetSpacing())
        self.dimension = np.asarray(self.image.GetSize())


    def make_rigid(self, translation, rotation, center_rotation=None):
        self.rigid.set_rotation(rotation)
        self.rigid.set_translation(translation)
        if center_rotation:
            self.rigid.set_center_rotation(center_rotation)
        self.rigid.write_file(filename=join(self.store_path,'rigid.txt'))


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
        self.bspline.make_coeffs(max_displacement=self.max_displacement, number_of_spots=self.number_of_spots, constrain=self.constrain, remove_bones=self.remove_bones, image=self.image)
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


    def deform_image_rigid(self):
        # deform image
        call_params = [self.plastimatch_executable]
        call_params.append('warp')

        call_params.append('--input')
        call_params.append(join(self.store_path,'input.nii'))

        call_params.append('--output-img')
        call_params.append(join(self.store_path,'warped.nii'))

        call_params.append('--xf')
        call_params.append(join(self.store_path,'rigid.txt'))

        call_params.append('--default-value')
        call_params.append('-1000') # assume air

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
    #parser.add_argument('-b', '--bspline',  help='bspline file')
    parser.add_argument('-i', '--input', help='path to DICOM input image to be deformed', required=True)
    parser.add_argument('-d', '--diffeomorphic', dest='diffeo', help='Constrain displacement to be diffeomorphic', action='store_true')
    parser.add_argument('-m', '--maxdispl', help='Maximum displacement, in image coordinates', type=float)
    parser.add_argument('-n', '--numspots', help='Number of gaussian spot displacements', type=int)
    parser.add_argument('-g', '--gridspacing', nargs='+', help='grid spacing in voxels', type=int)
    parser.add_argument('-r', '--rigid', nargs='+', help='rigid transform parameters (3 rotation angles in radians, 3 translation in mm)', type=float)
    parser.add_argument('-c', '--cor', nargs='+', help='center of rotation (3 coordinates x,y,z)', type=float)
    parser.add_argument('-o', '--output', help='output directory', required=True)
    parser.add_argument('-s', '--bones', help='try to constrain deformations to soft tissue', action='store_true')
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

    deformer = DeformImage(read_path=args.input, store_path=args.output, constrain=args.diffeo, number_of_spots=args.numspots, max_displacement=args.maxdispl, remove_bones=args.bones)
    deformer.read_image()
    if args.rigid:
        if len(args.rigid) == 6:
            print('rotation:', args.rigid[0:3])
            print('translation: ', args.rigid[3:6])
            if args.cor and len(args.cor) == 3:
                deformer.make_rigid(rotation=args.rigid[0:3], translation=args.rigid[3:6], center_rotation=args.cor)
            else:
                deformer.make_rigid(rotation=args.rigid[0:3], translation=args.rigid[3:6])
            deformer.deform_image_rigid()
        else:
            print("Error: Rigid parameters must be of length 6 (3 rotation angles, 3 translation)")
            exit()
    else:
        deformer.deform_image(args.gridspacing)
    deformer.write_image()

    # plastimatch convert --input CBCT --output-img cbct.nii
    # plastimatch warp cbct.nii
