# to create a CBCT from a CT and a CBCT
#python3 environment

import glob
import os
import pydicom as dicom
import numpy as np

import matplotlib.pyplot as plt
import random

from argparse import ArgumentParser
from skimage.transform import resize

from deform import read_contours_from_rs, create_rs, translate

from datetime import date

# from locations import isDicom, find_image_location, find_structure_location

from shutil import copy

from os.path import join, basename

class CT():
    def __init__(self, path = '.', image_prefix = 'CT', rs_prefix='rs'):
        self.images = self.load_images( path = path, prefix = image_prefix)
        self.rs     = self.load_rs( path = path, prefix=rs_prefix)


    def load_images(self, path = None, prefix='image'):
        dcms = glob.glob(os.path.join(path, prefix + '*'))
        # dcms = glob.glob(os.path.join(pt_path, '0*.dcm'))
        if len(dcms) <= 10:
            # print("No images")
            # exit(0)
            return None

        slices = [dicom.read_file(dcm) for dcm in dcms if not (dcm.endswith('xml') or basename(dcm).startswith('RS'))]
        slices.sort(key=lambda x: float(x.ImagePositionPatient[2]))

        return slices

    def load_rs(self, path = None, prefix = 'rs'):
        '''
        only search rs in the same folder or inside folder
        :param prefix:
        :return:
        '''
        if prefix == None or prefix == '': #don't have RS
            return None

        rs_file = find_rs( path = path, prefix= prefix)
        if rs_file is None or rs_file == '':
            return None

        rsrs = dicom.read_file( rs_file)
        return rsrs

    # def find_rs(selfs, path = None, prefix = 'rs'):
    #     #prefix and post fix
    #     rs = glob.glob(os.path.join(path, prefix + '*.dcm'))
    #     if len( rs) == 1:
    #         return rs[0]
    #
    #     #prefix only
    #     rs = glob.glob(os.path.join(path, prefix))
    #     if len( rs) == 1:
    #         return rs[0]
    #
    #     #go out of the image folder and look for it
    #     rs = glob.glob(os.path.join(path, '../', prefix + '*.dcm'))
    #     if len(rs) == 1:
    #         return rs[0]
    #
    #     # prefix only
    #     rs = glob.glob(os.path.join(path, '../', prefix, '*')) #use prefix as folder
    #     if len(rs) == 1:
    #         return rs[0]
    #
    #     print( "Warning only: No RS file for {}. Check the RS name or location".format(path))
    #     return None
    #     #exit(0)

def find_rs( path, prefix = 'rs'):
    # prefix and post fix
    rs = glob.glob(os.path.join(path, prefix + '*.dcm'))
    if len(rs) == 1:
        return rs[0]

    # prefix only
    rs = glob.glob(os.path.join(path, prefix))
    if len(rs) == 1:
        return rs[0]

    # go out of the image folder and look for it
    rs = glob.glob(os.path.join(path, '../', prefix + '*.dcm'))
    if len(rs) == 1:
        return rs[0]

    # prefix only
    rs = glob.glob(os.path.join(path, '../', prefix, '*'))  # use prefix as folder
    if len(rs) == 1:
        return rs[0]

    print("Warning only: No RS file for {}. Check the RS name or location".format(path))
    return None
    # exit(0)

def resample( ct, cbct):
    '''
    resample ct to cbct
    :param ct:
    :param cbct:
    :return:
    '''
    #cbct is the standard, the end goal
    ct_volume = np.stack([ s.pixel_array for s in ct.images], axis = 0).astype( np.float32)

    # cbct_volume =
    NEW_SIZE = cbct.images[0].Rows #desired size
    PIXEL_SPACING = cbct.images[0].PixelSpacing[0]
    inplane_scale = ct.images[0].PixelSpacing[0] / PIXEL_SPACING
    inplane_size = int( np.rint( inplane_scale * ct.images[0].Rows))

    SLICE_THICKNESS = cbct.images[0].SliceThickness
    zplane_scale = ct.images[0].SliceThickness / SLICE_THICKNESS
    zplane_size = int( np.rint( zplane_scale * len( ct.images)))

    if inplane_size == cbct.images[0].Rows and zplane_scale == 1:
        return ct_volume

    Volume = resize( ct_volume,
                     (zplane_size, inplane_size, inplane_size), mode= 'constant')

    if inplane_scale > 1:
        #crop
        return crop( Volume, NEW_SIZE)
    else:
        #pad
        return pad( Volume, NEW_SIZE)

def resample_z( ct, cbct):

    SLICE_THICKNESS = cbct.images[0].SliceThickness
    zplane_scale = ct.images[0].SliceThickness / SLICE_THICKNESS
    zplane_size = int(np.rint(zplane_scale * len(ct.images)))
    zs = np.array([slice.ImagePositionPatient[2] for slice in ct.images])

    zs_new = resize( zs, zplane_size, mode='constant')
    return zs_new

def crop(array, new_size):
    '''
    Crop a Volume or Mask in xy axis to fit to new image size
    :param array: 3D array, new_size is just for inplane size
    :return:
    '''
    #crop = int((self.inplane_scale['size'] - new_size) / 2)
    crop    = int(( array.shape[2] - new_size) / 2)
    array = array[:, crop : -crop , crop : -crop]
    return array

def pad(array, new_size):
    #pad = int((new_size - self.inplane_scale['size']) / 2)
    pad = int(( new_size - array.shape[2]) / 2)
    array = np.pad(array, ((0, 0), (pad, pad), (pad, pad)),'constant')
    return array

def ct2cbct( ct = None, cbct = None, path = None, output_prefix='CT.', typeofheader = 'CBCT'):
    '''
    A wrapper function fo 2 sub functions: cbct_cbctheader an cbct_ctheader
    :param ct:
    :param cbct:
    :param path: output path of the created CBCT
    :param output_prefix:
    :param typeofheader: 1: cbctheader, 2: ctheader
    :return:
    '''

    if not os.path.exists(path):
        os.mkdir(path)

    if ct.images is None:
        print( "No CT image")
        return 0

    if typeofheader == 'CBCT':
        return cbct_cbctheader( ct = ct, cbct = cbct, path = path)
    else:
        return cbct_ctheader( ct = ct, cbct = cbct, path = path)


def cbct_cbctheader(ct=None, cbct = None, path = None, output_prefix = 'CT.'):
    '''
    Cut the ct input into cbct by using cbct header, start at the Preview structure if available
    else start at the mid slice of PTV
    :param ct:
    :param cbct:
    :param path: path to deform CBCT output: e.g. ptname/CBCT/Deformed_Rigid
    :param output_prefix:
    :return:
    '''
    if not cbct:
        print( "No cbct for cbctfomcbct")
        return

    mid_slice = find_mid_slice(ct, cbct)
    nSlices = find_number_of_slices(ct, cbct)

    #cut use previous working header
    cbct_deformed = write_header_cbct(ct, cbct, mid_slice, nSlices, path, output_prefix)

    #create new RS file from CT RS by shifting contour point from CT coordinate to CBCT point coordinate
    # create_RS_from_CBCT( ct, cbct, mid_slice, nSlices, path, output_name = 'RS.CBCT.dcm')

    return cbct_deformed


def create_RS_from_CBCT( ct, cbct, mid_slice, nSlices, path, output_name = 'RS.CBCT.dcm'):
    #compute translation
    #(x, y, z) = compute_translation()
    xyz_cbct = cbct.images[0].ImagePositionPatient

    xyz_ct = ct.images[ mid_slice - nSlices//2].ImagePositionPatient
    (x, y, z) = np.array( xyz_cbct) - np.array(xyz_ct)

    if ct.rs:
        (contours, structure) =  read_contours_from_rs( ct.rs)
    #new_contours:
    for key, value in contours.items():
        #new_contour = translate( contours, x, y, z)
        contours[ key] = translate( contours[key], y, x, z)
    create_rs( join( path, output_name), structure, contours)
    # dicom.dcmwrite(join( path, output_name), structure)

def cbct_ctheader(ct = None, cbct = None, path = None, output_prefix='CT.'): #change back cbct_ctheader when tested
    '''
    Cut the ct input into cbct without cbct header, start at the Preview structure if available
    else start at the mid slice of PTV
    :param ct:
    :param cbct: not nessary, maybe just for number of slices
    :param path:
    :param output_prefix:
    :return:
    '''
    mid_slice = find_mid_slice( ct, cbct)
    nSlices = find_number_of_slices( ct, cbct)

    cbct_deformed = write_header_ct(ct,  mid_slice, nSlices, path, output_prefix)

    #copy the RS file from CT to CBCT
    dicom.dcmwrite( join( path, 'RS.dcm'), ct.rs)

    return cbct_deformed

def write_header_ct(ct, mid_slice, nSlices, path, output_prefix):
    offset = mid_slice -nSlices//2
    new_cbct_images = []
    for i in range( nSlices):
        new_cbct_image = ct.images[i + offset]
        new_cbct_images.append( new_cbct_image)
        dicom.dcmwrite( os.path.join( path, output_prefix+ str(i) +'.dcm' ), new_cbct_image)

    return new_cbct_images

def write_header_cbct( ct, cbct, mid_slice, nSlices_wrong, path, output_prefix):
    # resample ct to cbct resolution
    new_volume = resample(ct, cbct)
    # new_z       = resample_z( ct, cbct)

    new_cbct_images = cbct.images.copy()

    if not os.path.exists(path):
        os.mkdir(path)

    # CT_len = len( ct.images)
    #find the start slice in resample CT volume and cut from theree
    #due to resample
    new_mid_slice = mid_slice * len( new_volume) //len( ct.images)

    #number of slices in syth CBCT must be equal to CBCT
    nSlices = len( cbct.images)
    offset = new_mid_slice - nSlices// 2

    # copy from ct to cbct
    startUID = random.randint(1, 1000)
    incUID = random.randint(1, 1000)

    #get current date
    today = date.today()
    dd = today.strftime('%Y%m%d')

    for i in range( nSlices):
        new_cbct_image = new_cbct_images[i]

        # new_cbct_image.pixel_array[:] = ct.images[i + offset].pixel_array[:]
        new_cbct_image.pixel_array[:] = new_volume[i + offset]
        new_cbct_image.PixelData = new_cbct_image.pixel_array.tobytes()

        # change UID
        # the matter here is the '2.' to make Aria think this is a new study of new patient
        # Aria still can locate to the exist patient but also have the option of create a new one
        new_cbct_image.SOPInstanceUID = str(startUID) + new_cbct_image.SOPInstanceUID
        new_cbct_image.StudyInstanceUID = str(startUID + incUID) + new_cbct_image.StudyInstanceUID
        new_cbct_image.SeriesInstanceUID = str(startUID) + new_cbct_image.SeriesInstanceUID
        # new_cbct_image.PatientID = "Pydicom2"

        #change the Date of creation so that it's display on Eclipse
        new_cbct_image.SeriesDate = dd
        new_cbct_image.AcquisitionDate = dd
        new_cbct_image.ContentDate = dd


        # change slope and intercept
        new_cbct_image.RescaleIntercept = ct.images[0].RescaleIntercept
        new_cbct_image.RescaleSlope = ct.images[0].RescaleSlope
        new_cbct_image.HighBit = ct.images[0].HighBit


        # write
        dicom.dcmwrite(os.path.join(path, output_prefix + str(i) + '.dcm'), new_cbct_image)

    return new_cbct_images


def write_header_cbct_wrong( ct, cbct, mid_slice, nSlices, path, output_prefix):
    # resample ct to cbct resolution
    new_volume = resample(ct, cbct)
    # new_z       = resample_z( ct, cbct)

    new_cbct_images = cbct.images.copy()

    if not os.path.exists(path):
        os.mkdir(path)

    # CT_len = len( ct.images)
    #find the start slice in resample CT volume and cut from theree
    #due to resample
    new_mid_slice = mid_slice * len( new_volume) //len( ct.images)
    offset = new_mid_slice - nSlices// 2

    # copy from ct to cbct
    startUID = random.randint(1, 7)
    incUID = random.randint(1, 2)

    #get current date
    today = date.today()
    dd = today.strftime('%Y%m%d')

    for i in range( nSlices):
        new_cbct_image = new_cbct_images[i]

        # new_cbct_image.pixel_array[:] = ct.images[i + offset].pixel_array[:]
        new_cbct_image.pixel_array[:] = new_volume[i + offset]
        new_cbct_image.PixelData = new_cbct_image.pixel_array.tobytes()

        # change UID
        # the matter here is the '2.' to make Aria think this is a new study of new patient
        # Aria still can locate to the exist patient but also have the option of create a new one
        new_cbct_image.SOPInstanceUID = str(startUID) + new_cbct_image.SOPInstanceUID
        new_cbct_image.StudyInstanceUID = str(startUID + incUID) + new_cbct_image.StudyInstanceUID
        new_cbct_image.SeriesInstanceUID = str(startUID) + new_cbct_image.SeriesInstanceUID
        # new_cbct_image.PatientID = "Pydicom2"

        #change the Date of creation so that it's display on Eclipse
        new_cbct_image.SeriesDate = dd
        new_cbct_image.AcquisitionDate = dd
        new_cbct_image.ContentDate = dd


        # change slope and intercept
        new_cbct_image.RescaleIntercept = ct.images[0].RescaleIntercept
        new_cbct_image.RescaleSlope = ct.images[0].RescaleSlope
        new_cbct_image.HighBit = ct.images[0].HighBit

        # write
        dicom.dcmwrite(os.path.join(path, output_prefix + str(i) + '.dcm'), new_cbct_image)

    return new_cbct_images


def find_number_of_slices( ct, cbct):
    if cbct.images is None:
        nSlices = len( ct.images) //3
    else:
        #need to consider z scale of ct relative to cbct
        scale = cbct.images[0].SliceThickness / ct.images[0].SliceThickness
        nSlices = int( len(cbct.images) * scale)
    return nSlices

def find_mid_slice(ct = None, cbct = None):
    '''
    Finding mid slice in ct by one of the methods, depending on the data available:
    1. if ct rs Preview available: the preview plane is the mid slice
    2. if rs available: middle of PTV is mid slice
    3. else: mid of ct
    for number slice:
    1. If cbct availalbe: no slices = len(cbct)
    2. else: half of len(ct)
    :param ct:
    :param cbct:
    :return: mid slice and number of slices to cut
    '''
    if not ct.rs:
        mid_slice =  len( ct.images) //2
        return mid_slice

    #read ct.rs to decide if Preview available
    (contours, structure) = read_contours_from_rs( ct.rs)
    z_mid = find_mid_slice_z( contours = contours, iso_center_name= 'PREVIEW')
    # z_mid = 322.3
    if z_mid > 0: # i.e. PREVIEW contour inform available
        zs = np.array([ slice.ImagePositionPatient[2] for slice in ct.images ])
        index_mid = np.searchsorted( zs, z_mid, side='left')
        return index_mid

    else:
        mid_slice = len(ct.images) // 2
        return mid_slice

def find_mid_slice_z( contours = None, iso_center_name = 'PREVIEW'):
    '''
    if Preview
    :param contours:
    :return:
    '''
    # if iso_center_name in
    # if not iso_center_name in contours.keys():
    #     print( "There is no {} iso center structure".format( iso_center_name))
    #     return 0 #no isocenter
    for key in contours.keys():
        if iso_center_name.upper() in key.upper():
            iso_contour_name = key
            break
    else:
        print( "There is no {} iso center structure".format( iso_center_name))
        return 0 #no isocenter

    contour = contours[iso_contour_name]['contour']
    if len( contour) != 1:
        print("Ambiguous contours for isocenter. Select the mid slice in PREVIEW contour")
        return contour[len(contour)//2][2]

    return contour[0][2] #the z


if __name__ == '__main__':
    # ct_path = '/home/radonc/Projects/Halcyon/Data/Ready3/Neo0819liv002/CT/Deformed_Rigid'
    # cbct_path = '/home/radonc/Projects/Halcyon/Data/Ready3/Neo0819liv002/CBCT/FX1_PLFBH'
    # output_path = '/home/radonc/Projects/Halcyon/Data/Ready3/Neo0819liv002/CBCT/Deformed_Rigid Buggy'

    ct_path = '/home/radonc/Projects/Halcyon/Data/ExportedFromHalcyon/PTW/CT_warp'
    cbct_path = '/home/radonc/Projects/Halcyon/Data/ExportedFromHalcyon/PTW/CBCT'
    output_path = '/home/radonc/Projects/Halcyon/Data/ExportedFromHalcyon/PTW/CBCT_warp1'

    ct_path = '/home/radonc/Projects/Halcyon/Data/Emulator/NeoProst002/Fx1/plastimatch_dicom'
    cbct_path = '/home/radonc/Projects/Halcyon/Data/Emulator/NeoProst002/Origin/CBCT'
    output_path = '/home/radonc/Projects/Halcyon/Data/Emulator/NeoProst002/Fx1/CBCT direct'

    image_prefix_ct = 'image'
    rs_prefix_ct = 'rs'
    image_prefix_cbct = 'CT'
    rs_prefix_cbct = ''
    output_prefix = 'CT'
    typeofheader = 'CBCT'

    # parser = ArgumentParser()
    # parser.add_argument( '-c','--ct', help='input CT folder')
    # parser.add_argument('-b','--cbct', help = 'input CBCT, leave blank if no CBCT available')
    # parser.add_argument('-o', '--output', help = 'output CBCT')
    # parser.add_argument('--output_prefix', default='CT', help ='Prefix of the deformed CBCT file name')
    # parser.add_argument('--typeofheader', default='CT', help='which header to use to create deformed CBCT header')
    #
    # parser.add_argument('--image_prefix_ct', default='CT', help = 'CT image prefix')
    # parser.add_argument('--image_prefix_cbct', default='CT', help ='CBCT image prefix')
    #
    # parser.add_argument('--rs_prefix_ct', default='RS', help = 'CT RS prefix')
    # parser.add_argument('--rs_prefix_cbct', default='RS', help ='CBCT RS prefix')
    #
    # args = parser.parse_args()
    # #
    # ct_path = args.ct
    # cbct_path = args.cbct
    # output_path = args.output
    # output_prefix = args.output_prefix
    # typeofheader = args.typeofheader
    #
    # image_prefix_ct = args.image_prefix_ct
    # image_prefix_cbct = args.image_prefix_cbct
    # rs_prefix_ct = args.rs_prefix_ct
    # rs_prefix_cbct = args.rs_prefix_cbct

    ct = CT( path = ct_path, image_prefix= image_prefix_ct, rs_prefix=rs_prefix_ct)

    cbct = CT( path = cbct_path, image_prefix= image_prefix_cbct, rs_prefix=rs_prefix_cbct)

    ct2cbct(ct, cbct, output_path, output_prefix=output_prefix, typeofheader=typeofheader)

