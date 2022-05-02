# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 19:58:57 2019

@author: SomasundarA
"""

import os
import tkinter as tk

from file_handling import read_dicom_files
from dicom_validation import validate_dicom_image
from dicom_handling import sort_dicom_files
from nifti_handling import convert_to_nifti_based_on_scan_series_type


def dcm2nii_main(input_dicom_dir, output_dir):
    
    study_name: str = os.path.basename(input_dicom_dir)
    
    dicom_dataset = read_dicom_files(input_dicom_dir)
    
    sorted_dicom_dataset = sort_dicom_files(dicom_dataset)
    
    validate_dicom_image(sorted_dicom_dataset[0])
    
    convert_to_nifti_based_on_scan_series_type(
        sorted_dicom_dataset, study_name, output_dir)
    
 
# Get the input DICOM directory through User Interface
# TODO: check how to make the window appear in front in stead of minimized
root = tk.Tk() # Create a blank window
root.withdraw() # We don't want a full GUI, so keep root window from appearing
input_dicom_dir = os.path.abspath(tk.filedialog.askdirectory(
    title="Please select the input DICOM folder"))

## Get the output directory for NIfTI file(s) through User Interface
# TODO: check how to make the window appear in front in stead of minimized
root.deiconify() # make the window visible again
root.withdraw()
output_dir = os.path.abspath(tk.filedialog.askdirectory(
    title="Please select a folder to save output NIfTI file(s)"))
root.destroy() # destroy the root window along with all other tkinter widgets

dcm2nii_main(input_dicom_dir, output_dir)
