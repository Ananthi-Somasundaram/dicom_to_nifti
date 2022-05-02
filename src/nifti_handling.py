import sys
import math

import numpy as np
import nibabel as nib

from numpy import ndarray

from dicom_handling import (
    get_img_volume,
    compute_image_volume_array,
    get_pixel_dimensions_in_mm,
    get_rotation_matrix,
    split_dynamic_scan_by_frames)
from dicom_validation import is_static_or_whole_body_scan, is_dynamic_scan
from file_handling import create_nifti_output_path, save_nifti_image_to_file


def convert_lps_to_ras(img_volume_array_in_LPS: ndarray) -> ndarray:
    """Convert from LPS (DICOM) to RAS  format"""
    return np.flipud(np.fliplr(img_volume_array_in_LPS))


def compute_image_center(voxel_size, img_volume_array):
    """Compute the center of the image."""
    # Center of image. Shape -1 or +1?
    return (voxel_size * img_volume_array.shape) / 2


def create_affine_matrix(voxel_size, center) -> ndarray:
    """Create a affine matrix for use in NIfTI file creation"""
    affine_matrix: ndarray = np.diag(
        [voxel_size[0], voxel_size[1], voxel_size[2], 1])
    affine_matrix[:3,3] = np.array([-center[0], -center[1], -center[2]])

    return affine_matrix


def create_nifti_image(img_volume_array, affine_matrix) -> nib.Nifti1Image:
    """Apply the affine matrix to the volume array to create NIfTI image"""
    return nib.Nifti1Image(img_volume_array, affine_matrix)


def calculate_quaternions(R): 
    """
    R: Rotation Matrix
    TODO: Docstring
    """
    
    quaternion_a = 0.5 * math.sqrt(1 + \
                        R[0,0] + \
                        R[1,1] + \
                        R[2,2]) # not stored
    
    if quaternion_a != 0: # simplest case
        quaternion_b = 0.25 * (R[2,1] - R[1,2]) / quaternion_a
        quaternion_c = 0.25 * (R[0,2] - R[2,0]) / quaternion_a
        quaternion_d = 0.25 * (R[1,0] - R[0,1]) / quaternion_a
        
    return quaternion_b, quaternion_c, quaternion_d


def set_nifti_image_headers(
        nifti_image: nib.Nifti1Image, 
        img_volume_array, 
        rotation_matrix, 
        xyz_unit: str) -> nib.Nifti1Image:
    """Set a series of headers for the NIfTI image"""
    
    
    nifti_image.header.set_xyzt_units(xyz=xyz_unit)
    nifti_image.header.set_data_offset(352)
    
    nifti_image.header["qform_code"] = 1
    quatern_b, quatern_c, quatern_d = calculate_quaternions(rotation_matrix)
    nifti_image.header["quatern_b"] = quatern_b
    nifti_image.header["quatern_c"] = quatern_c
    nifti_image.header["quatern_d"] = quatern_d
    
    nifti_image.header["sform_code"] = 2 
    nifti_image.header["intent_code"] = 0 # None
    nifti_image.header["intent_name"] = "PET"
    nifti_image.header["cal_max"] = np.max(img_volume_array)
    nifti_image.header["cal_min"] = np.min(img_volume_array)
    
    return nifti_image


def convert_to_nifti(dcm_dataset) -> nib.Nifti1Image:
    """
    1. Get the volume array in LPS orientation of DICOM
    1. Convert the volume array to RAS orientation of NIfTI
    2. Apply slope and intercept to the volume array
    3. Apply the affine matrix to the volume array to create NIfTI image
    4. Set the NIfTI image header
    5. Save the NIfTI file
    """  
    img_volume_array_in_LPS: ndarray = get_img_volume(dcm_dataset)

    img_volume_array_in_RAS: ndarray = convert_lps_to_ras(
        img_volume_array_in_LPS)
    
    img_volume_array = compute_image_volume_array(
        dcm_dataset, img_volume_array_in_RAS)
    
    voxel_size: ndarray = get_pixel_dimensions_in_mm(dcm_dataset)
    center = compute_image_center(voxel_size, img_volume_array)
    affine_matrix: ndarray = create_affine_matrix(voxel_size, center)
    
    nifti_image: nib.Nifti1Image = create_nifti_image(\
        img_volume_array, affine_matrix)
    
    nifti_image.set_data_dtype(np.float32)
    
    xyz_unit = "mm"
    rotation_matrix = get_rotation_matrix(dcm_dataset)
    
    nifti_image: nib.Nifti1Image = set_nifti_image_headers(
        nifti_image, 
        img_volume_array, 
        rotation_matrix, 
        xyz_unit)
    
    return nifti_image

    
def convert_to_nifti_based_on_scan_series_type(
        sorted_dicom_dataset, study_name: str, output_dir: str) -> None:
    """Call the DICOM to NIfTI functionality, depending on scan series type"""
    frame_number = 1
    nr_of_slices = sorted_dicom_dataset[0].NumberOfSlices               
    scan_series_type = sorted_dicom_dataset[0].SeriesType[0]        
        
    if is_static_or_whole_body_scan(scan_series_type):             
        
        scan_duration_in_sec = int(sorted_dicom_dataset[0].ActualFrameDuration / 1000)                    
        nifti_image: nib.Nifti1Image = convert_to_nifti(sorted_dicom_dataset)  
        
        nifti_output_path: str = create_nifti_output_path(
            output_dir, frame_number, scan_duration_in_sec, study_name)
        
        save_nifti_image_to_file(nifti_image, nifti_output_path)
    
    elif is_dynamic_scan(scan_series_type):
        
        dcm_dataset_split_by_frames: list = split_dynamic_scan_by_frames(
            sorted_dicom_dataset, nr_of_slices)
        
        for dcm_dataset in dcm_dataset_split_by_frames:
                            
            scan_duration_in_sec = int(dcm_dataset[0].ActualFrameDuration / 1000)
            nifti_image: nib.Nifti1Image = convert_to_nifti(dcm_dataset) 
            
            nifti_output_path: str = create_nifti_output_path(
                output_dir, frame_number, scan_duration_in_sec, study_name)
            
            save_nifti_image_to_file(nifti_image, nifti_output_path)
            
            frame_number += 1        
                            
    else:
        sys.exit("Invalid scan series type")
