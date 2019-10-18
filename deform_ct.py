from ct2cbct import CT

import pydicom as dicom
import os, random
from datetime import date
from os.path import join
import glob
import subprocess

from fix_ct_number import write_header_ct2ct
import numpy as np

import pandas as pd

    # Explanation:
    #
    # 1. plastimatch synth-vf --fixed input.nii --xf-gauss --output gaussian/vf.mha --gauss-mag 10 --gauss-std 20 --gauss-center "-10 120 710"
    # 2. plastimatch warp --input input.nii --output-img gaussian/warped.nii --xf gaussian/vf.mha
    # 3. plastimatch convert --input gaussian/warped.nii --output-dicom gaussian/dicom --patient-id 10081008 --patient-name="E2E^ETHOS"
    # 4. write header ct2ct to bring the pixel_value from deform_dicom ct_dicom
    # 5. plastimatch warp --input plan_dicom/RS.PELCIRSWU19.Original071719.dcm --referenced-ct new_dicom/ --output-dicom gaussian/dicom --xf gaussian/vf.mha --patient-id 10081008 --patient-name E2E^ETHOS
    #
    #
    # Input:
    #     Origin/
    #         -Origin plan_dicom (folder)
    #             -CTs
    #             -RS
    #         -input.nii
    #         -....
    #
    # Output:
    #     Origin/
    #     Pt name1/
    #         -gaussian
    #             -dicom
    #                 -CTs (deformed)
    #             -vf.mha
    #             -warped.nii
    #         -new_plan_dicom
    #             -CTs
    #             -RS
    #
    #     Pt name2/
    #         -gaussian
    #             -dicom
    #                 -CTs (deformed)
    #             -vf.mha
    #             -warped.nii
    #         -new_plan_dicom
    #             -CTs
    #             -RS
    #
    # Middle:
    #     - Gaussian params
    #     - Pt names, Pt IDs


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
    call_params.append( join( output_loc, 'vf.mha'))

    #gaussian params
    call_params.append("--gauss-mag")
    call_params.append( gaussian_param.Gauss_mag)
    call_params.append("--gauss-std")
    call_params.append( gaussian_param.Gauss_std)
    call_params.append("--gauss-center")
    call_params.append( gaussian_param.Gauss_center)

    subprocess.call(call_params)


def warp( nii_path, pt):

    output_loc = join(os.path.dirname(nii_path), '..', pt.PatientName, 'gaussian')

    call_params = ['plastimatch']
    call_params.append("warp")
    call_params.append("--input")
    call_params.append(nii_path)
    call_params.append("--output-img")
    call_params.append(join( output_loc, 'warped.nii'))
    call_params.append( "--xf")
    call_params.append( join( output_loc,"vf.mha"))


    subprocess.call(call_params)

def convert( nii_path,pt):
    output_loc = join(os.path.dirname(nii_path), '..', pt.PatientName, 'gaussian')

    call_params = ['plastimatch']
    call_params.append("convert")
    call_params.append("--input")
    call_params.append(join( output_loc, 'warped.nii'))
    call_params.append( "--output-dicom")
    call_params.append(join(output_loc, 'dicom'))

    subprocess.call(call_params)

def warp_rs( nii_path, pt):

    plan_ct_path = join(os.path.dirname(nii_path), plan_ct_folder_name)
    deform_ct_path = join(os.path.dirname(nii_path), '..', pt.PatientName, 'new_ct')

    RSs = glob.glob( join( plan_ct_path, 'RS.*' ))
    if len(RSs) != 1:
        print( "NO RS structure! Stop")
    else:
        RS = RSs[0]

    vf = join( os.path.dirname(nii_path),'..', pt.PatientName, 'gaussian','vf.mha')

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

if __name__ == '__main__':
    origin_path = '/home/radonc/Projects/Halcyon/Data/CT_Deform/E3E'
    plan_ct_folder_name = 'PelvisWithBowel'

    pts = [
        ["E2E^JCPELVIS6", "10292019"],
        ["E2E^JCPELVIS7", "10302019"],
        ["E2E^JCPELVIS8", "10312019"],
        ["E2E^JCPELVIS9", "10322019"],
        ["E2E^JCPELVIS10", "10332019"],
        ["E2E^JCPELVIS11", "10342019"],
        ["E2E^JCPELVIS12", "10352019"],
        ["E2E^JCPELVIS13", "10362019"],
    ]

    # for random sampling
    mag_range = [10, 20]
    std_range = [15, 35]
    x_range = [-21.1, 1.2]
    y_range = [115.4, 126]
    z_range = [698., 716.]

    #fix random
    np.random.seed(1)

    niis = sorted( glob.glob(join( origin_path,'*.nii')))

    params = []

    for i in range( 0, len(pts)):
        pt = Patient( pts[i][0], pts[i][1])
        # gaussian = Gaussian( gaussian_params[i][0], gaussian_params[i][1], gaussian_params[i][2])
        mag, std, x, y, z = gaussian_sampling( mag_range, std_range, x_range, y_range, z_range)
        # params.append( [mag, std, x, y, z])

        # xyz_string = " ".join( [str(x), str(y), str(z)])
        gaussian = Gaussian( mag, std, (x, y, z))
        params.append( [gaussian.Gauss_mag, gaussian.Gauss_std, gaussian.Gauss_center])
        gen_vf_ma( niis[0], gaussian, pt)
        warp( niis[0], pt)
        convert( niis[0], pt)

        #
        write_header( niis[0], pt)

        warp_rs(niis[0], pt)

    df = pd.DataFrame( params)
    print( df)
    df.to_csv('Gaussian_params.txt', index=False, header=['mag', 'std', 'center'])