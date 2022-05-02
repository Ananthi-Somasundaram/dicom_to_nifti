import os
import sys

import pydicom
import nibabel as nib

from typing import Tuple


def read_dicom_files(input_dicom_dir):
    """
    Read the DICOM dataset from the input directory. Set force parameter to 
    True if DICOM dataset not stored in accordance with Standard File Format.
    """

    # TODO: Add more unwanted suffixes as required!
    unwanted_suffixes: Tuple[str] = (
        '.db','.xml','.voi','.prj','.csv','.bmp','.nii','.pdf')

    try:
        dcm_dataset = [pydicom.dcmread(os.path.join(input_dicom_dir, filename), 
                                       force = False) \
                        for filename in os.listdir(input_dicom_dir) \
                        if not filename.endswith(unwanted_suffixes)]
        print("Valid DICOM: " + input_dicom_dir)
    except:
        sys.exit("Not a valid DICOM: " + input_dicom_dir)
        
    return dcm_dataset


def get_dicom_header_data(dcm_set, output_filename, output_dir):
    """
    Dump the DICOM header into a text file with the same filename as NIFTI 
    Set anonymize_flag = "Y" to not add patient data tags in the text file
    """
    
    anonymize_flag = "Y"
    
    patient_tags = ['PatientID', 'PatientName', 'PatientBirthDate']
      
    output_txt_filename = output_filename + ".txt"
    
    with open(os.path.join(output_dir,output_txt_filename), "w") as file:
        for header_tag in dcm_set[0].iterall():
            if anonymize_flag == "Y":
                if header_tag not in [dcm_set[0].data_element(tag) for tag in patient_tags]:
                    file.write(str(header_tag) + '\n')
            elif anonymize_flag == "N":
                file.write(str(header_tag) + '\n')    
                
                
def create_nifti_output_path(output_dir: str, frame_number, scan_duration_in_sec, study_name: str) -> str:
    """Create a full path for NIfTI image storage"""
    nifti_output_filename: str = \
        f"s{frame_number}_{scan_duration_in_sec}s_{study_name}.nii"
    
    return os.path.join(output_dir, nifti_output_filename)


def save_nifti_image_to_file(nifti_image: nib.Nifti1Image, output_path: str) -> None:
    """Save the NIfTI image to a designated output path"""
    
    # TODO: re-incorporate dicom header data storage
    # get_dicom_header_data(dcm_set, output_filename, output_dir)
    
    nib.save(nifti_image, output_path)  
    print(f"NIfTI file successfully saved at {output_path}")
