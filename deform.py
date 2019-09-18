#python3.6
import os
import pydicom
import glob

import numpy as np


from skimage.draw import polygon

import matplotlib.pyplot as plt

#simple first
#Read contour data to contours
#distort some contour in contours
#write contours in to rs

#complicate later

rs_path = '../Data/CIRS/RS.1.2.246.352.205.5027836166574254264.11682058419526807215.dcm'
new_rs_path = '../Data/CIRS/RS_new.dcm'

distort_contour_names = ['Bladder', 'Rectum']

distort_params = {'Bladder': {'Translation': (10, 10)},

                  'Rectum': {'Zoom': 1.5}}

ROI_ORDER   = ['Bladder', 'BODY', 'Femur_L', 'Femur_R', 'PTV_7920', 'Rectum']

def read_contours_from_rs(structure):
    contours = {}

    # read all the contours, could be a function itself
    for i in range(len(structure.ROIContourSequence)):
        contour = {}
        contour['name'] = structure.StructureSetROISequence[i].ROIName
        # if not contour['name'] in ROI_ORDER:
        #     continue
        contour['number'] = structure.ROIContourSequence[i].ReferencedROINumber
        assert contour['number'] == structure.StructureSetROISequence[i].ROINumber
        try:
            s = structure.ROIContourSequence[i].ContourSequence
        except:
            print("There is no contour for: ", contour['name'])
            continue
        contour['contour'] = [s.ContourData for s in structure.ROIContourSequence[i].ContourSequence]

        contours[contour['name']] = (contour)

    # only extract ROI_ORDER
    # we need this clear out the name resolution
    return (contours, structure)

def read_contours(rs_path):
    '''
    read and return structure and contours info
    :param rs_path:
    :return:
    '''
    structure = pydicom.dcmread(rs_path)

    return read_contours_from_rs( structure)

def distort_contour( contour, method_names):
    '''
    distorting one contour, with method
    :param contour:
    :method_name: 'Bladder': {'Translation': (x, y'),
                                'Zoom': z,
                                'Rotation': angle}
    :return: distorted contour
    '''
    distorted_contour = contour.copy()
    for key in method_names.keys():
        if key == 'Translation':
            (x, y) = method_names[key]
            distorted_contour = translate( distorted_contour, x, y)

        if key == 'Zoom':
            z = method_names[key]
            distorted_contour = zoom(distorted_contour, z)
        if key == 'Rotate':
            angle = method_names[key]
            distorted_contour = rotate(distorted_contour, angle)

    # distorted_contour = {}
    return distorted_contour


def translate( contour, x, y, z = 0):
    contour_distorted = contour.copy()
    for i, slice_contour in enumerate(contour['contour']):
        slice = np.reshape( slice_contour, (-1, 3))
        slice += np.array([x, y, z])
        slice_distort = np.round( slice.flatten(), 4)

        contour_distorted['contour'][i] = slice_distort
    return contour_distorted

def zoom( contour, z):
    '''
    This piece of code amazes me how it works
    :param contour:
    :param z: zoom factor
    :return:
    '''
    contour_distorted = contour.copy()
    for i, slice_contour in enumerate(contour['contour']):
        slice = np.reshape( slice_contour, (-1, 3))

        #zoom
        p = slice[:, 0:2]
        center = (p.min( axis =0) + p.max(axis = 0)) * .5
        q = center *(1 - z) + p* z #which is center + (p-center)* z

        slice[:, 0:2] = q
        slice_distort = np.round( slice.flatten(), 4) #round to 4 digit
        contour_distorted['contour'][i] = slice_distort

    return contour_distorted

def rotate( contour, angle):
    '''
    Haven't done any thing
    :param contour:
    :param angle:
    :return:
    '''
    # the right tthing is
    p = slice[:, 0:2]
    center = (p.min + p.max) *.5
    rot_mat = np.array( [[np.cos( angle), -np.sin(angle)],
                         np.sin(angle), np.cos(angle)])
    q =  (p - center) * rot_mat + center
    return contour

def distorting_contours(contours, distort_params):

    distorted_contours = contours.copy()
    for key, val in distort_params.items():
        contour = contours[key]
        distorted_contour = distort_contour( contour, val)

        distorted_contours[key] = distorted_contour
    return distorted_contours

def create_rs_very_limited(new_rs_file, structure, contours):
    #update structure
    #update structure name
    structure.StructureSetLabel = "BladderandRectum"
    structure.StructureSetName  = 'BladderandRectum'
    structure.SOPInstanceUID = '1.2.826.0.1.3680043.8.498.67063073255766457385690235251075316643'
    distort_contour_names = distort_params.keys()
    for i in range(len(structure.ROIContourSequence)):
        contour_name = structure.StructureSetROISequence[i].ROIName
        if contour_name in distort_contour_names:
            S = structure.ROIContourSequence[i].ContourSequence
            for slice_index in range( len(contours[contour_name]['contour'])):
                S[slice_index].ContourData = list( contours[contour_name]['contour'][slice_index])


    pydicom.dcmwrite( new_rs_file, structure)


def create_rs(new_rs_file, structure, contours):
    #update structure
    #update structure name
    distort_contour_names = contours.keys()
    for i in range(len(structure.ROIContourSequence)):
        contour_name = structure.StructureSetROISequence[i].ROIName
        if contour_name in distort_contour_names:
            S = structure.ROIContourSequence[i].ContourSequence
            for slice_index in range( len(contours[contour_name]['contour'])):
                S[slice_index].ContourData = list( contours[contour_name]['contour'][slice_index])


    pydicom.dcmwrite( new_rs_file, structure)


if __name__ == '__main__':
    (contours, structure) = read_contours( rs_path)


    distorted_contours = distorting_contours( contours, distort_params )

    create_rs( new_rs_path, structure, distorted_contours)