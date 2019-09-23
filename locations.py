import glob
from os.path import join, isdir, basename
from os import listdir, makedirs, replace
from shutil import copy, copytree, move
import numpy as np
import pydicom as dicom

minImages = 100

def isDicom(file_path, modality = 'RTSTRUCT'):
    '''
    use exception to test if the
    :param file_path:
    :return:
    '''
    try:
        dcm = dicom.dcmread( file_path)
        if dcm.Modality.upper() == modality:
            return 1
        else:
            return 0
    except:
        return 0

def find_image_location(patients_folder, patient, prefix = 'CT', subfolders = None):
    '''
    CT can be exported in current folder
    or one lower folder in the hierarchical
    :param patient_path:
    :return:
    '''
    # find dicom images location
    # either in Planning or Planning/Dicom
    image_location = join( patients_folder, patient)
    dcms = glob.glob(join(image_location, prefix, '*'))
    if len( dcms) > minImages:
        return image_location

    #if image in a hierachical subfolders:
    for subfolder in subfolders:
        image_location = join( image_location, subfolder)
        dcms = glob.glob(join(image_location, prefix, '*'))
        if len(dcms) > minImages:
            return image_location

    else:
        print("There is no CT images")
        exit( 1)

def find_structure_location( patients_folder, patient, prefix = 'RS', subfolders = None):
    '''
    :param patient_path:
    :return:
    '''
    structure_location = join( patients_folder, patient)
    dcms = glob.glob(join(structure_location, prefix, '*'))
    if len(dcms) ==1:
        if isDicom( dcms[0], modality='RTSTRUCT'):
            return structure_location

    #if image in a hierachical subfolders:
    for subfolder in subfolders:
        structure_location = join( structure_location, subfolder)
        dcms = glob.glob(join(structure_location, prefix, '*'))
        if isDicom( dcms[0], modality='RTSTRUCT'):
            return dcms[0]

    else:
        print( 'There is no Structure file')
        return None
