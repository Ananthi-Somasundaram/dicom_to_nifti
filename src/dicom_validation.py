import sys

def is_pet_image(image) -> bool:
    """
    Check if the image modality is PET, as indicated by the DICOM tag.
    (0008,0060) Modality Tag: Type of equipment that acquired the data
    """
    return image.Modality == 'PT' \
        and image.SOPClassUID == '1.2.840.10008.5.1.4.1.1.128'

def decay_and_attenuation_correction_applied(image) -> bool:
    """
    Check if decay and attenuation correction are applied to the PET data.
    (0028,0051) Corrected Image Tag: Corrections that applied to the image.
    """
    return 'DECY' and 'ATTN' in image.CorrectedImage


def images_decay_corrected_to_scan_start(image) -> bool:
    """
    Check if decay correction is applied to the image at scan start
    (0054,1102) Decay Correction Tag: The real-world event to which images in 
    the Series were decay corrected.
    """
    return image.DecayCorrection == 'START' 


def units_either_bqml_or_counts(image) -> bool:
    """
    Check if image units are within the set of allowed units. 
    (0054,1001) Units Tag: Pixel value units - BQML, CNTS, etc.
    """
    return image.Units == 'BQML' or 'CNTS'


def is_gated_scan(image) -> bool:
    """Check if the scan type is GATED, this is currently unsupported"""
    return image.SeriesType[0] == "GATED"  


def validate_dicom_image(image) -> None:
    """Perform a series of checks to validate the image"""
    
    if not is_pet_image(image):
        sys.exit("Not a PET image. Aborting")
        
    if not decay_and_attenuation_correction_applied(image): 
        sys.exit("Decay and attenuation correction not applied. Aborting")
        
    if not images_decay_corrected_to_scan_start(image):
        sys.exit("Images decay not corrected to scan start. Aborting")
        
    if not units_either_bqml_or_counts(image): 
        sys.exit("Units not either BQML or CNTS. Aborting.")
        
    if is_gated_scan(image):
        sys.exit("GATED scans not supported yet")

        
def is_static_or_whole_body_scan(scan_series_type: str) -> bool:
    """Check whether the scan is static or whole body"""
    return (scan_series_type == 'STATIC' or scan_series_type == 'WHOLE BODY')


def is_dynamic_scan(scan_series_type: str) -> bool:
    """
    Check if the scan series is dynamic, indicating need for extra processing
    """
    return scan_series_type == 'DYNAMIC'
