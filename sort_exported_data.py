import pydicom as dicom
import numpy as np
from os.path import join
import glob
import matplotlib.pyplot as plt


def read_data_slice_only(ct_loc):
    dcms = glob.glob(join(pt_loc, ct_loc, 'CT*'))
    print(len(dcms))
    slices = [dicom.dcmread(dcm) for dcm in dcms if not dcm.endswith('xml')]
    slices.sort(key=lambda x: (x.Rows))
    return slices


def sorted_mixed_ct(slices):
    '''
    Input: mixed bag of CTs
    return a list of list
    '''
    slices_sorted = []
    nRows = 0
    k = -1  # empty list
    for i in range(len(slices)):
        if slices[i].Rows != nRows:
            nRows = slices[i].Rows
            slices_sorted.append([])
            k += 1
            slices_sorted[k].append(slices[i])

        else:
            slices_sorted[k].append(slices[i])

    return slices_sorted


def sorted_ct_images(slices):
    '''
    Input: one CT
    :param slices:
    :return: sorted CT and iamges
    '''
    slices.sort(key=lambda x: float(x.ImagePositionPatient[2]))
    images = np.stack(s.pixel_array * slices[0].RescaleSlope + slices[0].RescaleIntercept for s in slices)
    sopuids = [slice.SOPInstanceUID for slice in slices]
    return slices, images, sopuids


def sorted_umbrella(ct_loc):
    slices = read_data_slice_only(ct_loc)
    slices = sorted_mixed_ct(slices)
    images = [[] for i in range(len(slices))]
    sopuids = [[] for i in range(len(slices))]
    for k in range(len(slices)):
        slices[k], images[k], sopuids[k] = sorted_ct_images(slices[k])

    return slices, images, sopuids

pt_loc = '../Data/ExportedFromHalcyon/Aug19/NeoFpel002/Sessions/'
ct_loc = 'Session_1'

Slices, Images, Sopuids = sorted_umbrella( ct_loc)

rts = glob.glob( join( pt_loc, ct_loc, 'RTSTR*'))

for i in range( len( rts)):
    rs = dicom.dcmread( rts[i])
    for k in range( len( Sopuids)):
        if rs.ROIContourSequence[0].ContourSequence[0].ContourImageSequence[0].ReferencedSOPInstanceUID in Sopuids[k]:
            print(rts[i], ' belongs to CT ', k)
            #Start copy this RS and CT in Sopuids to new folder
            
            break

