# # Prologue
# This version is an combination of the advantage of v2 and deform_ct_plastimatch
# It is supposes to reduce the disadvantage of v2:
# Features:
# Gaussian deformed at maybe multiple location
# Clean code of organizing Gaussian params
# Remove the rot and bspline


from ct2cbct import CT, find_rs, ct2cbct

import pydicom as dicom
import os, random
from datetime import date
from os.path import join, exists
import glob
import subprocess

from fix_ct_number import write_header_ct2ct
import numpy as np

import pandas as pd

import SimpleITK as sitk

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
#
# Middle:
#     - Gaussian params
#     - Pt names, Pt IDs

#some CONST_NAME
#due to conflict with module name
CT_CT = 'CT'
RS_CT = 'RS'
CBCT_CT = 'CBCT'

INPUT_NII = 'input.nii'
ORIGINAL = 'Origin'
GAUSSIAN = 'gaussian'

WARP_NII = 'warped.nii'
VF_SUM = 'vf_sum.mha'

PLASTIMATCH_IMAGE_PREFIX = 'image'
PLASTIMATCH_RS_PREFIX = 'rtss'
PLASTIMATCH_DICOM = 'plastimatch_dicom'

ORIGINAL_CT_PREFIX = 'CT'
ORIGINAL_CBCT_PREFIX = 'CT'

class Gaussian():
    def __init__(self, mag, std, center):
        '''

        :param mag: x
        :param std: x
        :param center: (x,y,z)
        '''
        signs = np.random.choice([-1, 1], size=3)
        Gauss_mag = mag * signs
        self.Gauss_mag = " ".join( map( str, Gauss_mag))
        self.Gauss_std = str(std)
        Gauss_center = " ".join( map( str, center))
        self.Gauss_center = Gauss_center #a string like "-10 120 710"

class Patient():
    def __init__(self, name, id):
        self.PatientName = name
        self.PatientID = id

class PT_Info():
    #df dataframe with following fields:
    #PT, PT_id, PT_path, Site, Fx, nGaussian, Center, Magnitude, Std, Note
    def __init__(self, df):
        pt_path = df['PT_path'].iloc[0]
        nii_path = join(pt_path, ORIGINAL, INPUT_NII)
        fx = df['Fx'].iloc[0]
        Fx_folder = 'Fx' + str(fx)
        Fx_path = join( pt_path, Fx_folder)
        origin_ct_path = join( pt_path, ORIGINAL, CT_CT)
        origin_cbct_path = join( pt_path, ORIGINAL, CBCT_CT)

        plastimatch_dicom_path = join( Fx_path, PLASTIMATCH_DICOM)

        deformed_ct_path = join( Fx_path, CT_CT)
        deformed_cbct_path = join( Fx_path, CBCT_CT)

        vf_path = join( Fx_path, GAUSSIAN, VF_SUM)

        PT_id = df['PT_id'].iloc[0]
        PT_name = df['PT'].iloc[0]

        self.nii_path = nii_path
        self.pt_path = pt_path
        self.origin_ct_path = origin_ct_path
        self.origin_cbct_path = origin_cbct_path
        # self.Fx_folder = Fx_folder
        self.Fx_path = Fx_path
        self.plastimatch_dicom_path = plastimatch_dicom_path
        self.PatientID = PT_id #to match the legacy code
        self.PatientName = PT_name #to match the legacy

        self.deformed_ct_path = deformed_ct_path
        self.deformed_cbct_path = deformed_cbct_path

        self.vf_path = vf_path


def gen_vf_one(df, i):
    #gnerate one by one vector field
    pt_path = df['PT_path'].iloc[i]
    fx = df['Fx'].iloc[i]
    Fx_folder = 'Fx' + str(fx)
    nGaussians = df['nGaussian'].iloc[i]

    Center = df['Center'].iloc[i]
    Magnitude = df['Magnitude'].iloc[i]
    Std = df['Std'].iloc[i]

    nii_path = join(pt_path, ORIGINAL, INPUT_NII)

    output_loc = join( pt_path, Fx_folder, GAUSSIAN)
    output_name = 'vf' + str(nGaussians) +'.mha'

    call_params = ["plastimatch"]
    call_params.append("synth-vf")
    call_params.append("--fixed")
    call_params.append(nii_path)
    call_params.append("--xf-gauss")
    call_params.append("--output")
    call_params.append(join(output_loc, output_name))

    #gaussian params
    call_params.append("--gauss-mag")
    call_params.append( Magnitude)
    call_params.append("--gauss-std")
    call_params.append( Std)
    call_params.append("--gauss-center")
    call_params.append( Center)

    subprocess.call(call_params)

def gen_vf_ma( nii_path, gaussian_param, pt):

    #lcoation of guassian folder
    output_loc = join( os.path.dirname(nii_path),'..', pt.PatientName, 'gaussian')
    if not os.path.isdir( output_loc):
        os.makedirs( output_loc)

    call_params = ["plastimatch"]
    call_params.append("synth-vf")
    call_params.append("--fixed")
    call_params.append( nii_path)
    call_params.append( "--xf-gauss")
    call_params.append("--output")
    call_params.append( join( output_loc, output_name))

    #gaussian params
    call_params.append("--gauss-mag")
    call_params.append( gaussian_param.Gauss_mag)
    call_params.append("--gauss-std")
    call_params.append( gaussian_param.Gauss_std)
    call_params.append("--gauss-center")
    call_params.append( gaussian_param.Gauss_center)

    subprocess.call(call_params)

def warp( df):
    '''
    apply vector field to the input nii
    :return:
    '''

    pt_info = PT_Info(df)
    nii_path = pt_info.nii_path
    output_loc = join( pt_info.Fx_path)

    call_params = ['plastimatch']
    call_params.append("warp")
    call_params.append("--input")
    call_params.append(nii_path)
    call_params.append("--output-img")
    call_params.append(join( output_loc, WARP_NII))
    call_params.append( "--xf")
    call_params.append( join( output_loc, GAUSSIAN, VF_SUM))


    subprocess.call(call_params)

def convert(df):
    '''
    convert from warped.nii to CT folder
    '''
    pt_info = PT_Info(df)
    output_loc = join( pt_info.Fx_path)

    call_params = ['plastimatch']
    call_params.append("convert")
    call_params.append("--input")
    call_params.append(join( output_loc, WARP_NII))
    call_params.append( "--output-dicom")

    call_params.append(join(output_loc, PLASTIMATCH_DICOM))

    subprocess.call(call_params)

def write_ct_header(df):
    #bring the anonymiyze in plastimatch_dicom to CT
    #with PT Name and correct CT number that can be planned

    pt_info = PT_Info( df)

    plan_ct_path = pt_info.origin_ct_path

    deform_ct_path = pt_info.plastimatch_dicom_path

    plan_ct = CT(path = plan_ct_path, image_prefix=ORIGINAL_CT_PREFIX, rs_prefix=RS_CT)
    deform_ct = CT(path = deform_ct_path, image_prefix=PLASTIMATCH_IMAGE_PREFIX, rs_prefix=PLASTIMATCH_RS_PREFIX)

    write_header_ct2ct(deform_ct, plan_ct, pt_info.deformed_ct_path, pt_info)

def write_cbct_header(df):
    #cut ct to cbct with kick from another cbct
    pt_info = PT_Info( df)

    deform_ct_path = pt_info.deformed_ct_path #can't used this one, so confused!
    deform_ct_path = pt_info.plastimatch_dicom_path
    cbct_path = pt_info.origin_cbct_path
    deform_cbct_path = pt_info.deformed_cbct_path

    deform_ct =     CT( path = deform_ct_path, image_prefix = PLASTIMATCH_IMAGE_PREFIX, rs_prefix= PLASTIMATCH_RS_PREFIX)
    cbct      = CT( path = cbct_path, image_prefix = ORIGINAL_CBCT_PREFIX, rs_prefix= None)

    ct2cbct(deform_ct, cbct, deform_cbct_path, output_prefix= CBCT_CT, typeofheader=CBCT_CT)

def warp_rs( df ):
    '''
    apply vector to Structure file
    :param df:
    :return:
    '''
    pt_info = PT_Info( df)

    plan_ct_path = pt_info.origin_ct_path
    deform_ct_path =  pt_info.deformed_ct_path

    #find RS
    RS = find_rs( plan_ct_path, prefix = RS_CT)
    if RS is None:
        print( "NO RS structure! Stop")
        return

    # vf = join( os.path.dirname(nii_path),'..', pt.PatientName, 'gaussian','vf.mha')
    vf = pt_info.vf_path

    call_params = ['plastimatch']
    call_params.append("warp")
    call_params.append("--input")
    call_params.append(RS)
    call_params.append('--referenced-ct')
    call_params.append( deform_ct_path)
    call_params.append('--output-dicom')
    call_params.append(deform_ct_path)

    call_params.append("--xf")
    call_params.append(vf)

    subprocess.call(call_params)


def write_header( nii_path, pt):

    output_loc = join(os.path.dirname(nii_path), '..', pt.PatientName)

    # plan_ct_path = '/home/radonc/Projects/Halcyon/Data/E2E/plan_dicom'
    plan_ct_path = join( os.path.dirname(nii_path), plan_ct_folder_name)
    # deform_ct_path = '/home/radonc/Projects/Halcyon/Data/E2E/deform_dicom'
    deform_ct_path = join( output_loc, 'gaussian', 'dicom')

    plan_ct = CT(path = plan_ct_path, image_prefix='CT', rs_prefix='RS')
    deform_ct = CT(path = deform_ct_path, image_prefix='image', rs_prefix='rtss')

    write_header_ct2ct(deform_ct, plan_ct, join( output_loc,'new_ct'), pt)

def gaussian_sampling(mag_range, std_range, x_range, y_range, z_range):
    '''
    Sampling a mag, std inside the prostate
    :param mag_range:
    :param std_range:
    :param x_range:
    :param y_range:
    :param z_range:
    :return:
    '''

    r = np.random.random(1)[0]
    mag = (mag_range[1] - mag_range[0])*r + mag_range[0]
    r = np.random.random(1)[0]
    std = (std_range[1] - std_range[0])*r + std_range[0]
    r = np.random.random(1)[0]
    x = (x_range[1] - x_range[0])*r + x_range[0]
    r = np.random.random(1)[0]
    y = (y_range[1] - y_range[0])*r + y_range[0]
    r = np.random.random(1)[0]
    z = (z_range[1] - z_range[0])*r + z_range[0]

    return (mag, std, x, y, z)


def gen_nii(pt_path):
    #check if nii exist, generate if not
    pt_info = PT_Info(df)

    nii_path = pt_info.nii_path

    if exists(nii_path):
        return

    #plastimatch convert --input dir --output-img input.nii
    ct_path = pt_info.origin_ct_path

    call_params = ['plastimatch']
    call_params.append("warp")
    call_params.append("--input")
    call_params.append( ct_path)
    call_params.append( '--output-img')
    call_params.append( nii_path)

    subprocess.call(call_params)

def gen_vf( df):
    '''Generat Gaussian vector field by generate the component Gaussian and add them up using SITK'''

    pt_info = PT_Info(df)
    output_loc = join( pt_info.Fx_path, GAUSSIAN)

    output_name = 'vf' + '_sum.mha'
    if exists( join( output_loc, output_name)):
        return

    #generate component Gaussian
    for i in range(len(df)):
        gen_vf_one( df, i)

    # #use siltk to sum up all vfs
    nGaussian = df['nGaussian'].iloc[0]
    output_name = 'vf' + str(nGaussian) + '.mha'
    vf = sitk.ReadImage( join( output_loc, output_name))
    for i in range(1, len(df)):
        nGaussian = df['nGaussian'].iloc[i]
        output_name = 'vf' + str(nGaussian) + '.mha'
        vf1 = sitk.ReadImage( join( output_loc, output_name))
        vf = vf + vf1

    # #write vf to mha
    output_name = VF_SUM
    sitk.WriteImage( vf, join( output_loc, output_name))

def generate_deform(df):
    #df dataframe that can have multiple Gaussian
    #let's work with the first one

    #1. generate .nii
    #2. generate vector field
        #Prerequisite:
            # a.Gaussian object
                # Prerequisite:
                    # i. Center of gaussian
                    # ii. Magnitude
                    # iii. Std

        #Postrequisite:
            # a.Join multiple vector field

    #3. warp nii with vf to dicom
    #4. Change CT number as well some info
    #5. Warp_RS with new_warped CT folder and
    #6. Cut the CT to CBCT
        #Prerequisite:
            # a. CBCT header
            # b. Preview ROI to cut
            # b. The previous code

    #gnerate nii
    gen_nii(df)
    #
    gen_vf(df)
    #
    warp(df)
    #
    convert(df)
    #
    write_ct_header(df) #just wrap of write_header from other function
    #
    warp_rs(df)

    write_cbct_header(df) #just wrap

if __name__ == '__main__':
    params_file = './PTandGaussians_good_v3.csv'

    dataframe = pd.read_csv(params_file, delimiter='\t')
    # df = dataframe.groupby(['PT', 'Fx'])
    PTs = dataframe['PT'].unique()
    Fxs = dataframe['Fx'].unique()

    for pt in PTs:
        for fx in Fxs:
            df = dataframe[(dataframe['PT'] == pt) & (dataframe['Fx'] == fx)]
            if len(df):
                generate_deform(df)

