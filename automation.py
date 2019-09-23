#!python3

from art_testing.make_synthetic_image import DeformImage

from ct2cbct import CT, ct2cbct

from locations import isDicom, find_image_location, find_structure_location

from roi_deform import deform_roi

import glob
from os.path import join, isdir, basename
from os import listdir, makedirs, replace
from shutil import copy, copytree, move
import numpy as np
import pydicom as dicom
#
minImages = 100
input_folder = '/home/radonc/Projects/Halcyon/Data/test_folder'
output_folder = '/home/radonc/Projects/Halcyon/Data/Ready'
patient_list = sorted(listdir(input_folder))

# bspline distort params:
# isSpline = 1
maxdispl = 20
numspot = 20
gridspacing = [15, 15, 15]
diffeo = True
remove_bones = True

# rigid
PI = 3.14159
a_scale = 1  # degre
t_scale = 0
rigid = [0 * a_scale * PI / 180, 0 * a_scale * PI / 180, 1 * a_scale * PI / 180, t_scale * 5, t_scale * 10,
         t_scale * 2]  # rot - tran
cor = [0, 0, 0]

rigid_folder_name = 'Deformed_Rigid'
spline_folder_name = 'Deformed_Spline'
info_folder_name = 'Deformed Info'

RIGID_PARAM_FILE = 'rigid.txt'
SPLINE_PARAM_FILE = 'bspline.txt'


deform_tmp_folder = 'tmp'

plastimatch_prefix = 'image'
CT_origin_prefix = 'CT'
CT_deformed_prefix = 'CT.Deformed'
CBCT_deformed_prefix = 'CT.CBCT.Deformed'

def calculate_cor(image_location = None, abs_cor = [0,0,0]):

    def load_images(path = None, prefix='CT'):
        dcms = glob.glob(join(path, prefix + '*'))
        # dcms = glob.glob(os.path.join(pt_path, '0*.dcm'))
        if len(dcms) == 0:
            print("No images")
            exit(0)
        slices = [dicom.read_file(dcm) for dcm in dcms]
        slices.sort(key=lambda x: float(x.ImagePositionPatient[2]))

        return slices

    slices = load_images( image_location, prefix = 'CT')
    pt_origin = slices[0].ImagePositionPatient

    # #is it true????
    # pt_origin[0] = -pt_origin[0]
    # pt_origin[1] = -pt_origin[1]

    px_spacing = list( slices[0].PixelSpacing) + [slices[0].SliceThickness]
    dims = list( slices[0].pixel_array.shape) + [len( slices)]
    cor_from_center = np.array( dims) * np.array( px_spacing) /2 + pt_origin
        #[0, 0, 0]

    return cor_from_center 


def deform( image_location = './', output_location='./', isSpline = 1):
    deformer = DeformImage(read_path=image_location, store_path=output_location,
                           constrain=diffeo,
                           number_of_spots=numspot, max_displacement=maxdispl, remove_bones = remove_bones)
    deformer.read_image()

    if isSpline:
        deformer.deform_image(gridspacing)
    else:
        # rigid
        # cor_from_center = calculate_cor(image_location, cor)
        # forget about cor
        deformer.make_rigid(rotation=rigid[0:3], translation=rigid[3:6], center_rotation = cor)
        deformer.deform_image_rigid()

    deformer.write_image()


def structure_and_rename( image_location = None, structure_location = None,
                          output_location = None, param_folder = None, isSpline = 1):

    def move_deformed_files( ):#output_location  = None, param_folder  = None, isSpline = 1):
        '''
        move all the file (not the image in dicom folder,  in the tmp folder which is in out_location
        into param_folder
        :param output_location: tmp
        :param param_folder: Spline or Rigig
        :param isSpline: No need now
        :return:
        '''

        #rewrite move deformed file e.g. bspline
        destination_deformed_info_folder = join(output_location, info_folder_name, param_folder)
        if not isdir( destination_deformed_info_folder):
            makedirs( destination_deformed_info_folder)
        plastimatch_files = glob.glob( join( output_location, '*.*')) #not working if plastimatch doesn't have extension
        for file in plastimatch_files:
            filename = basename( file)
            move( file, join( output_location, info_folder_name, param_folder, filename))



        #move image in the dicom folder
        destination_deformed_ct_folder = join( output_location, 'CT', param_folder)
        if not isdir( destination_deformed_ct_folder):
            #create
            makedirs( destination_deformed_ct_folder)

        #move CT
        images = listdir( join( output_location, 'dicom') )
        for image in images:
            source = join( output_location, 'dicom', image)
            dest = join ( destination_deformed_ct_folder, image.replace(plastimatch_prefix, CBCT_deformed_prefix))
            move( source, dest)

 

    def copy_original():
        '''
        Move/rename original CT into outout CT folder
        Move Structrur ifle int on the smae  folder
        :return:
        '''
        if not isdir( join( output_location, 'CT')):
            makedirs( join( output_location, 'CT', 'CT'))
            # makedirs(join(output_location, info_folder_name))
            #copy original to output folder
        images = listdir( join( image_location))
        for image in images:
            if image.startswith('RP') or image.startswith('RS') or image.startswith('RD'):
                continue
            source = join( image_location, image)
            if image.startswith(CT_origin_prefix):
                destt = join( output_location, 'CT', 'CT', image)
            else:
                destt = join(output_location, 'CT', 'CT', 'CT.' + image)
            copy( source, destt)

        #copy RS to the same folder
        destination_ct_folder = join( output_location, 'CT', 'CT')
        rs_name = basename( structure_location)
        if not rs_name.startswith('RS'):
            rs_name = 'RS.'+ rs_name
        if not rs_name.endswith('.dcm'):
            rs_name = rs_name + '.dcm'
        copy( structure_location, join( destination_ct_folder, rs_name))

    copy_original()
    move_deformed_files()

  

def create_cbct( output_location, cbct_location, isSpline, typeofheader = 'CT'):
    if not isdir(join(output_location, 'CBCT')):
        makedirs(join(output_location, 'CBCT'))

        #copy from CBCT_location
        cbct_profiles = listdir( join( cbct_location))
        profile = cbct_profiles[0]
        for profile in cbct_profiles:
            #rename to .CT
            if not isdir(join( output_location, 'CBCT', profile)):
                makedirs( join( output_location, 'CBCT', profile))
            files = listdir( join( cbct_location, profile))
            cbct_origin = join( cbct_location, profile)
            if len(files) < minImages: #assumpt 50 slice in cbct
                files = listdir(join(cbct_location, profile, 'DICOM'))
                cbct_origin = join(cbct_location, profile, 'DICOM')
            for file in files:
                source = join( cbct_origin, file)
                destt = join( output_location, 'CBCT', profile, 'CT.'+file)
                copy( source, destt)

    cbct_profiles = sorted( listdir(join(cbct_location)))
    profile = cbct_profiles[-1]


    if isSpline:
        deform_folder_name = spline_folder_name
    else:
        deform_folder_name = rigid_folder_name

    ct_location_folder = join( output_location, 'CT', deform_folder_name)
    # cbct_location_folder = join( cbct_location, profile)
    cbct_location_folder = join(output_location, 'CBCT', profile)

    deformed_cbct_location_folder = join( output_location, 'CBCT', deform_folder_name)

    ct = CT( ct_location_folder, image_prefix='CT', rs_prefix = 'RS')

    cbct = CT( cbct_location_folder, image_prefix = 'CT', rs_prefix = None)
    ct2cbct( ct, cbct, deformed_cbct_location_folder, output_prefix = 'CT.Deformed.CBCT', typeofheader = typeofheader)

if __name__ == '__main__':

    # deformer = DeformImage(read_path=args.input, store_path=args.output, constrain=args.diffeo,
    #                        number_of_spots=args.numspots, max_displacement=args.maxdispl)

    subfolders_RS = ['Planning','RS', 'DICOM']
    # subfolders_CT = ['Planning', 'CT', 'DICOM']
    subfolders_CT = ['Planning',  'DICOM']
    prefix_ct = ''
    prefix_rs = ''

    for patient in patient_list:

        image_location = find_image_location(input_folder, patient, prefix=prefix_ct, subfolders= subfolders_CT)
        structure_location = find_structure_location( input_folder, patient, prefix=prefix_rs , subfolders= subfolders_RS)

        #debug
        print( image_location, structure_location)

        output_location = join( output_folder, patient)




        for isSpline in [0,1]:
            if isSpline:
                param_folder = spline_folder_name
                deform_param_file = join( output_location, SPLINE_PARAM_FILE)
            else:
                param_folder = rigid_folder_name
                deform_param_file = join( output_location, RIGID_PARAM_FILE)

            #do the common job
            #deform image
            deform(image_location=image_location, output_location=output_location, isSpline=isSpline)

            #deform ROI
            #rs, reference_ct, output, deform_param_file,
            deform_roi( structure_location, output_location, output_location, deform_param_file )

            structure_and_rename(image_location, structure_location, output_location, param_folder= param_folder)

            cbct_location = join( input_folder,patient, 'Treatment')

            # uncomment this if we want to create CBCT
            # create_cbct( output_location, cbct_location, isSpline, typeofheader='CT')