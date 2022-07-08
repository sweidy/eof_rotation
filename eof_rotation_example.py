"""
Calculating OMI with rotation

Use this file to generate the rotated OMI as described in Weidman, S., Kleiner, N., and Kuang, Z., 2022. 
A rotation procedure to improve seasonally varying Empirical Orthogonal Function bases for MJO indices, 
submitted to *Geophysical Research Letters*. (preprint can be found at https://doi.org/10.1002/essoar.10511626.1). 

The complete OMI algorithm is described in Kiladis, G.N., J. Dias, K.H. Straub, M.C. Wheeler, S.N. Tulich, K. Kikuchi, 
K.M. Weickmann, and M.J. Ventrice, 2014: A Comparison of OLR and Circulation-Based Indices for Tracking the MJO. 
Mon. Wea. Rev., 142, 1697–1715, https://doi.org/10.1175/MWR-D-13-00301.1

General OMI functions from mjoindices python package published by Hoffmann, C.G., Kiladis, G.N., Gehne, M. 
and von Savigny, C., 2021. A Python Package to Calculate the OLR-Based Index of the Madden- Julian-Oscillation 
(OMI) in Climate Science and Weather Forecasting. Journal of Open Research Software, 9(1), p.9. 
DOI: http://doi.org/10.5334/jors.331 which can be found on github in at https://github.com/cghoffmann/mjoindices.

Rotation algorithm functions found at https://github.com/sweidy/eof_rotation at eof_rotation.py. 
Additional file rotation_plotting_tools.py in the same folder is used for organizing and plotting the PCs. 
"""

from pathlib import Path
import os.path
import inspect

import mjoindices.olr_handling as olr
import mjoindices.omi.omi_calculator as omi
import mjoindices.empirical_orthogonal_functions as eof
import mjoindices.principal_components as pc
import mjoindices.evaluation_tools as eval_tools
import mjoindices.tools as tools
import mjoindices.omi.wheeler_kiladis_mjo_filter as wkfilter
import mjoindices.eof_rotation as omir

import numpy as np

# create filename paths for OLR, EOFs, and PCs

OLR_path = Path(os.path.abspath('')) / 'observed_data' / 'olr.day.mean.nc'
eofs_dir = Path(os.path.abspath('')) / 'EOFs'
pc_dir = Path(os.path.abspath('')) / 'PCs'

# load in OLR dataset
raw_olr = olr.load_noaa_interpolated_olr(OLR_path)

# Restrict the OLR dataset to the dates used in the paper. This can be changed to use a longer time period, 
# but EOF results will change slightly.
time_start = np.datetime64('1979-01-01')
time_end = np.datetime64('2012-12-31')

shorter_olr = olr.restrict_time_coverage(raw_olr, time_start, time_end)

# Interpolate the dataset to the same grid as used in the original OMI calcluation 
# (not required, but testing has only been done using the interpolated grid).

interpolated_olr = olr.interpolate_spatial_grid_to_original(shorter_olr)

# Calculate the EOFs using the rotation algorithm. This may take an hour or more, especially on a personal computer.
"""
    :param olrdata: The OLR dataset, from which OMI should be calculated. Note that OLR values are assumed to be given
        in positive values. The spatial grid of the OLR datasets defines also the spatial grid of the complete OMI
        calculation.
    :param implementation: eofs_package (py:mod:eofs) or internal method
    :param sign_doy1reference: Switch signs of EOF to align with 2014 paper (observations) so all days are consistent.
    :param strict_leap_year_treatment: default for no leap years or non-standard calendar should be False`.

"""  

rot_eofs = omir.calc_eofs_from_olr_with_rotation(interpolated_olr, 
                                 implementation='eofs_package', 
                                 sign_doy1reference=True,
                                 strict_leap_year_treatment=False)


# Save EOFs in your defined EOF folder. 
rot_eofs.save_all_eofs_to_npzfile(eofs_dir / 'EOFs_rotated.npz')

# Similarly, calculate the EOFs from the original OMI algorithm. 
# Change interpolate_eofs to True to directly reproduce the method from Kiladis et al., 2014. 
# This will interpolate the EOFs during the beginning of November. This also may take up to an hour or more.

norot_eofs = omi.calc_eofs_from_olr(interpolated_olr,
                        implementation="eofs_package",
                        sign_doy1reference=True,
                        interpolate_eofs=False,
                        strict_leap_year_treatment=False)


norot_eofs.save_all_eofs_to_npzfile(eofs_dir / 'EOFs_unrotated.npz')

# If you have already calculated the EOFs, uncomment and reload them from your EOF folder.
# eof_rot = eofs_dir / 'EOFs_rotated.npz'
# rot_eofs = eof.restore_all_eofs_from_npzfile(eof_rot)

# eof_norot = eofs_dir / 'EOFs_unrotated.npz'
# norot_eofs = eof.restore_all_eofs_from_npzfile(eof_norot)

# Now use the EOFs to calculate the principal components. 
# The time variable is restricted to the same years as the EOF calculation. 
# If you removed leap years earlier, the PCs will also use no leap years. 
pcs_rot = omi.calculate_pcs_from_olr(raw_olr,
                                 rot_eofs,
                                 time_start,
                                 time_end,
                                 use_quick_temporal_filter=False)

pcs_norot = omi.calculate_pcs_from_olr(raw_olr,
                                 norot_eofs,
                                 time_start,
                                 time_end,
                                 use_quick_temporal_filter=False)

pcs_rot.save_pcs_to_txt_file(pc_dir / 'PCs_rotated.txt')
pcs_norot.save_pcs_to_txt_file(pc_dir / 'PCs_unrotated.txt')

# If you have already calculated PCs, uncomment and reload them here.
pc_rot_path = pc_dir / 'PCs_rotated.txt'
pcs_rot = pc.load_pcs_from_txt_file(pc_rot_path)

pc_norot_path = pc_dir / 'PCs_unrotated.txt'
pcs_norot = pc.load_pcs_from_txt_file(pc_norot_path)