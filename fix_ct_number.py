from ct2cbct import CT

import pydicom as dicom
import os, random
from datetime import date

def explanation():
    '''
    1. plastimatch synth-vf --fixed input.nii --xf-gauss --output gaussian/vf.mha --gauss-mag 10 --gauss-std 20 --gauss-center "-10 120 710"
    2. plastimatch warp --input input.nii --output-img gaussian/warped.nii --xf gaussian/vf.mha
    3. plastimatch convert --input gaussian/warped.nii --output-dicom gaussian/dicom --patient-id 10081008 --patient-name="E2E^ETHOS"
    4. write header ct2ct to bring the pixel_value from deform_dicom ct_dicom
    5. plastimatch warp --input plan_dicom/RS.PELCIRSWU19.Original071719.dcm --referenced-ct new_dicom/ --output-dicom gaussian/dicom --xf gaussian/vf.mha --patient-id 10081008 --patient-name E2E^ETHOS


    Input:
    Origin/
        -Origin plan_dicom (folder)
            -CTs
            -RS
        -input1.nii
        -input2.nii
        -....

    Output:
    Origin/
    Pt name1/
        -gaussian
            -dicom
                -CTs (deformed)
            -vf.mha
            -warped.nii
        -new_plan_dicom
            -CTs
            -RS

    Pt name2/
        -gaussian
            -dicom
                -CTs (deformed)
            -vf.mha
            -warped.nii
        -new_plan_dicom
            -CTs
            -RS

    Middle:
    Pt names, Pt IDs
    '''
    pass

class Patient():
    def __init__(self, name, id):
        self.PatientName = name
        self.PatientID = id


def write_header_ct2ct( ct, cbct, path, pt, output_prefix = 'CT.'):
    '''
    basically just move the pixel value deform_ct to plan_ct
    May change some UID to crate a new one
    :param plan_ct:
    :param deform_ct:
    :return:
    '''
    # new_volume = resample(ct, cbct)
    new_volume = ct.images
    # new_z       = resample_z( ct, cbct)

    new_cbct_images = cbct.images.copy()

    if not os.path.exists(path):
        os.mkdir(path)

    # CT_len = len( ct.images)
    #find the start slice in resample CT volume and cut from theree
    #due to resample
    # new_mid_slice = mid_slice * len( new_volume) //len( ct.images)

    #number of slices in syth CBCT must be equal to CBCT
    nSlices = len( cbct.images)
    # offset = new_mid_slice - nSlices// 2
    offset = 0
    #
    # copy from ct to cbct
    startUID = random.randint(1, 1000)
    incUID = random.randint(1, 1000)

    #get current date
    today = date.today()
    dd = today.strftime('%Y%m%d')

    for i in range( nSlices):
        new_cbct_image = new_cbct_images[i]

        # new_cbct_image.pixel_array[:] = ct.images[i + offset].pixel_array[:]
        new_cbct_image.pixel_array[:] = \
            (ct.images[i + offset].pixel_array[:] * ct.images[0].RescaleSlope\
            + ct.images[0].RescaleIntercept - new_cbct_image.RescaleIntercept) / new_cbct_image.RescaleSlope

        new_cbct_image.PixelData = new_cbct_image.pixel_array.tobytes()

        # change UID
        # the matter here is the '2.' to make Aria think this is a new study of new patient
        # Aria still can locate to the exist patient but also have the option of create a new one
        new_cbct_image.SOPInstanceUID = str(startUID) + new_cbct_image.SOPInstanceUID
        new_cbct_image.StudyInstanceUID = str(startUID + incUID) + new_cbct_image.StudyInstanceUID
        new_cbct_image.SeriesInstanceUID = str(startUID) + new_cbct_image.SeriesInstanceUID

        # new_cbct_image.ReferencedSOPInstanceUID = str(startUID) + new_cbct_image.ReferencedSOPInstanceUID
        new_cbct_image.FrameOfReferenceUID = str(startUID) + new_cbct_image.FrameOfReferenceUID
        # new_cbct_image.PatientID = "Pydicom2"

        #change pt name and id
        new_cbct_image.PatientName = pt.PatientName
        new_cbct_image.PatientID = pt.PatientID

        #change the Date of creation so that it's display on Eclipse
        new_cbct_image.SeriesDate = dd
        new_cbct_image.AcquisitionDate = dd
        new_cbct_image.ContentDate = dd


        # change slope and intercept
        # new_cbct_image.RescaleIntercept = ct.images[0].RescaleIntercept
        # new_cbct_image.RescaleSlope = ct.images[0].RescaleSlope
        # new_cbct_image.HighBit = ct.images[0].HighBit


        # write
        dicom.dcmwrite(os.path.join(path, output_prefix + str(i) + '.dcm'), new_cbct_image)

    return new_cbct_images

if __name__ == '__main__':
    pt = Patient("E3E^ETHOS", '09182019')

    plan_ct_path = '/home/radonc/Projects/Halcyon/Data/E2E/plan_dicom'
    deform_ct_path = '/home/radonc/Projects/Halcyon/Data/E2E/deform_dicom'
    plan_ct = CT(path=plan_ct_path, image_prefix='CT', rs_prefix='RS')
    deform_ct = CT(path=deform_ct_path, image_prefix='image', rs_prefix='rtss')

    write_header_ct2ct( deform_ct, plan_ct, 'new_ct', pt)