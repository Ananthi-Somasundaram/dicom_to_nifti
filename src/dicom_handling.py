import numpy as np

from datetime import datetime
from typing import Tuple, List
from numpy import ndarray


patient_study = 1
phantom_study = 0
convert_to_suv_flag = 1


def get_x_direction_cosine(dicom_image) -> ndarray:
    """
    Get x-direction cosine of first row, first column w.r.t. the patient
    (0020,0037) Image Orientation Patient Tag
    """
    return np.array(list(map(float,dicom_image.ImageOrientationPatient[:3])))
    

def get_y_direction_cosine(dicom_image) -> ndarray:
    """
    Get y-direction cosine of first row, first column w.r.t. the patient
    (0020,0037) Image Orientation Patient Tag
    """
    return np.array(list(map(float,dicom_image.ImageOrientationPatient[3:])))
    

def get_z_direction_cosine(x_dir_cos: ndarray, y_dir_cos: ndarray) -> ndarray:
    """
    For vectors a and b, a × b vector is perpendicular to both a and b 
    and thus normal to the plane containing them
    """
    return np.cross(x_dir_cos, y_dir_cos)


def get_direction_cosines(dcm_dataset) -> Tuple:
    """Get the x, y, and z directional cosines"""
    x_dir_cos = get_x_direction_cosine(dcm_dataset[0])
    y_dir_cos = get_y_direction_cosine(dcm_dataset[0])
    z_dir_cos = get_z_direction_cosine(x_dir_cos, y_dir_cos)
    
    return (x_dir_cos, y_dir_cos, z_dir_cos)


def get_image_position_coordinates(dcm_dataset) -> List[float]:
    """
    Get the x, y, z coordinates of the upper left hand corner of the image, 
    (center of the first voxel transmitted) in mm
    (0054,0081): Number of slices per frame Tag
    """
    image_position_list: List[int] = [dcm_dataset[i].ImagePositionPatient._list \
            for i in range(0, dcm_dataset[0].NumberOfSlices)]
        
    return [[float(coordinate) for coordinate in image_position] \
        for image_position in image_position_list]


def get_z_coordinates(dcm_dataset):
    """Get the rotation matrix for the DICOM dataset"""
    x_dir_cos, y_dir_cos, z_dir_cos = get_direction_cosines(dcm_dataset)
    
    image_position_patient: List[float] = \
        get_image_position_coordinates(dcm_dataset)

    # Vector dot product or scalar inner product to get z coordinates
    z_coordinates: ndarray = np.dot(image_position_patient, z_dir_cos)
    
    return z_coordinates


def get_rotation_matrix(dcm_dataset):
    """Get the rotation matrix for the DICOM dataset"""
    x_dir_cos, y_dir_cos, z_dir_cos = get_direction_cosines(dcm_dataset)
    
    return np.array([x_dir_cos, y_dir_cos, z_dir_cos])


def sort_dicom_files(dcm_dataset):    
    """Sort the DICOM dataset into volume based on ImageIndex."""
    sorted_dcm_dataset = dcm_dataset.copy()
    sorted_dcm_dataset.sort(key = lambda x: int(x.ImageIndex))
    
    return sorted_dcm_dataset


def reorder_dicom_files(dcm_dataset):
    """Sort the DICOM dataset into volume based on the z_coordinates."""
    z_coordinates = get_z_coordinates(dcm_dataset)
    
    reordered_dcm_dataset = [x for _,x in sorted(zip(z_coordinates, dcm_dataset))]
    
    return reordered_dcm_dataset


def get_voxel_size_z(dcm_dataset):
    """
    The absolute difference between two consecutive z_coordinates gives the 
    voxel size along the z-direction, in mm. 
    """
    z_coordinates = get_z_coordinates(dcm_dataset) 
    
    voxel_size_z = abs(z_coordinates[0] - z_coordinates[1])
    
    return voxel_size_z


def get_img_volume(dcm_dataset):    
    """
    1. Get the voxel intensity array of each slice
    2. Transpose the array to account for the change in array order from 
       row-major to column-major
    3. Stack it in the z-direction
    """
    reordered_dcm_dataset = reorder_dicom_files(dcm_dataset)
    
    img_volume_array = np.transpose(reordered_dcm_dataset[0].pixel_array)
    
    for dcm_slice in reordered_dcm_dataset[1:]:
        img_slice_array = np.transpose(dcm_slice.pixel_array)
        img_volume_array = np.dstack((img_volume_array, img_slice_array))
        
    return img_volume_array


def calculate_decay_constant(dcm_dataset):
    
    half_life_of_radionuclide_in_hr = (dcm_dataset[0].RadiopharmaceuticalInformationSequence[0].RadionuclideHalfLife)/(60*60) # (0018, 1075) 

    decay_constant = np.log(2) / half_life_of_radionuclide_in_hr
    
    return decay_constant


def calculate_decay_factor(dcm_dataset):
    
    radiopharmaceutical_start_date_time = dcm_dataset[0].RadiopharmaceuticalInformationSequence[0].RadiopharmaceuticalStartDateTime # 20151103.062600 # injection date time: (0018, 1078) using the same time base as Series Time

    series_date = dcm_dataset[0].SeriesDate # (0008, 0021)
    
    series_time = dcm_dataset[0].SeriesTime # (0008, 0031)
    
    series_date_time = series_date + series_time # 20151109.013628 # scan start date time 
    
    fmt = '%Y%m%d%H%M%S.000000' # Change as per date time format in DICOM file # 20151109.013628 or # 20151109013628.000000
    time_stamp_1 = datetime.strptime(radiopharmaceutical_start_date_time, fmt)
    time_stamp_2 = datetime.strptime(series_date_time, fmt)
    
    time_difference = time_stamp_2 - time_stamp_1
    time_elapsed_in_hr = (time_difference.total_seconds() / (60*60))
    
    decay_constant = calculate_decay_constant(dcm_dataset)

    decay_factor = np.exp(-decay_constant * time_elapsed_in_hr)
    
    return decay_factor


def get_decay_corrected_activity_in_Bq_at_scan_start(dcm_dataset):
    
    injected_activity_in_Bq = float(dcm_dataset[0].RadiopharmaceuticalInformationSequence[0].RadionuclideTotalDose) # (0018, 1074)
    
    decay_factor = calculate_decay_factor(dcm_dataset)
    
    decay_corrected_activity_in_Bq_at_scan_start = injected_activity_in_Bq * decay_factor
    
    return decay_corrected_activity_in_Bq_at_scan_start


def convert_to_suv(dcm_dataset, img_volume_array):
    
    decay_corrected_activity_in_Bq_at_scan_start = get_decay_corrected_activity_in_Bq_at_scan_start(dcm_dataset)
    
    if patient_study == 1:
        
        patient_weight_in_g = (dcm_dataset[0].PatientWeight)*1000
    
        suv_scaling_factor = 1/(decay_corrected_activity_in_Bq_at_scan_start/patient_weight_in_g)
    
    elif phantom_study == 1:
        
        background_mean_ac_in_bqml = 1491.4 # get from ACCURATE
        
        suv_scaling_factor = 1/background_mean_ac_in_bqml
    
    suv_scaled_img_volume_array = img_volume_array * suv_scaling_factor
    
    return suv_scaled_img_volume_array
    

def get_pixel_dimensions_in_mm(dcm_dataset) -> ndarray:
    """
    Get the physical dimensions of the pixel.
    (0028,0030) - PixelSpacing Tag: Physical distance between the center of 
    each pixel. Row Spacing and Column Spacing in mm.
    """
    voxel_size_z = get_voxel_size_z(dcm_dataset)
    pixel_spacing = dcm_dataset[0].PixelSpacing
    
    return np.array(
        [float(pixel_spacing[0]), 
         float(pixel_spacing[1]), 
         voxel_size_z])


def get_intercept(dcm_dataset):
    """
    Get the intercept for the dataset.
    (0028,1052): Intercept Tag
    """
    return dcm_dataset[0].RescaleIntercept


def rescale_counts_units_to_bqml_units(dcm_dataset):
    """
    Rescale required if COUNTS units as all subsequent computations expect BQML 
    (7053,1009) Philips Activity Concentration Scale Factor Tag: Multiplying 
    stored pixel values by Rescale Slope then this factor results in Bq/ml.
    """
    return dcm_dataset[0].RescaleSlope * dcm_dataset[0][0x7053,0x1009].value


def is_in_counts_units(dcm_dataset) -> bool:
    """Check if the DICOM dataset is in COUNTS units"""
    return dcm_dataset[0].Units == "CNTS"


def get_slope(dcm_dataset):
    """
    Get the slope for the dataset
    (0028,1053): Slope Tag
    """
    
    if is_in_counts_units(dcm_dataset):
        return rescale_counts_units_to_bqml_units(dcm_dataset)
    else:
        return dcm_dataset[0].RescaleSlope


def compute_image_volume_array(dcm_dataset, img_volume_array_in_RAS):
    """
    Get the slope and intercept to compute the image volume array, change image 
    volume array datatype from int to float, apply slope and intercept
    """
    intercept = get_intercept(dcm_dataset)
    slope = get_slope(dcm_dataset)
    
    img_volume_array = (img_volume_array_in_RAS.astype(float) * slope) \
        + intercept

    if convert_to_suv_flag == 1:
        img_volume_array = convert_to_suv(dcm_dataset, img_volume_array)
    
    return img_volume_array
    

def split_dynamic_scan_by_frames(sorted_dcm_dataset: list, nr_of_slices: int) -> list:
    """Split a dynamic scan into separate frames"""
    return [
        # TODO: replace x with more descriptive name
        sorted_dcm_dataset[x : x + nr_of_slices] \
        for x in range(0, len(sorted_dcm_dataset), nr_of_slices)
        ]
