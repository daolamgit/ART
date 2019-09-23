import subprocess

import shutil

import glob

from os.path import join

PLASTIMATCH_ROI_NAME = 'rtss'
PLASTIMATCH_IMAGE_NAME = 'dicom'
ROI_FILE_NAME = 'RS.dcm'

def deform_roi( rs, reference_ct, output, deform_param_file, exe_name = 'plastimatch'):
    '''
    plastimatch warp --input RS.dcm --referenced-ct dicom --output-dicom /my/output/directory --xf bspline.txt
    :param rs:
    :param reference_ct:
    :param output:
    :param deform_param_path:
    :param exe_name:
    :return:
    '''
    if not rs:
        return
    call_params = [exe_name]

    call_params.append('warp')

    call_params.append('--input')
    call_params.append(rs)

    call_params.append('--referenced-ct')
    call_params.append( join( reference_ct, PLASTIMATCH_IMAGE_NAME))

    call_params.append('--output-dicom')
    call_params.append(output)

    call_params.append('--xf')
    call_params.append(deform_param_file)

    subprocess.call( call_params)

    #rename the output
    # a = 1
    rs_output = glob.glob( join( output, PLASTIMATCH_ROI_NAME + '*'))[0]
    image_outpput = join( output, PLASTIMATCH_IMAGE_NAME, ROI_FILE_NAME)
    shutil.move( rs_output, image_outpput)


    #fixing UID in new RS so that it's load in eclipse