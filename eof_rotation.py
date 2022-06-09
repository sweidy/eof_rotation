"""
This code is an add-on to the OMI calculation package published in Hoffmann, C.G., Kiladis, G.N., Gehne, M. and von Savigny, C., 2021. 
A Python Package to Calculate the OLR-Based Index of the Madden- Julian-Oscillation (OMI) in Climate Science and Weather Forecasting. 
Journal of Open Research Software, 9(1), p.9. DOI: http://doi.org/10.5334/jors.331
which can be found on github in at https://github.com/cghoffmann/mjoindices. 

The complete OMI algorithm is described in Kiladis, G.N., J. Dias, K.H. Straub, M.C. Wheeler, S.N. Tulich, K. Kikuchi, K.M.
Weickmann, and M.J. Ventrice, 2014: A Comparison of OLR and Circulation-Based Indices for Tracking the MJO.
Mon. Wea. Rev., 142, 1697–1715, https://doi.org/10.1175/MWR-D-13-00301.1

The function :meth:'calc_eofs_from_olr_with_rotation' included in this file can be used instead of the function 
:meth:'calc_eofs_from_olr' in the original mjoindices package. The new method includes a projection and rotation postprocessing
step that reduces noise in the original EOF calculation. 

The parameter :param no_leap: is currently not implemented in the original mjoindices package, but this is in progress. It may be
deleted for datasets that include leap years. 

Adaptation written by Sarah Weidman, 2022. Contact: sweidman@g.harvard.edu

"""

from pathlib import Path
from typing import Tuple
import os.path
import inspect

import numpy as np
import warnings
import importlib
import scipy.interpolate
import scipy.optimize as optimize
import scipy.linalg as linalg
import pandas as pd

import mjoindices.empirical_orthogonal_functions as eof
import mjoindices.olr_handling as olr
import mjoindices.omi.omi_calculator as omi
import mjoindices.principal_components as pc
import mjoindices.omi.wheeler_kiladis_mjo_filter as wkfilter
import mjoindices.omi.quick_temporal_filter as qfilter
import mjoindices.tools as tools

eofs_spec = importlib.util.find_spec("eofs")
eofs_package_available = eofs_spec is not None
if eofs_package_available:
    import eofs.standard as eofs_package


def angle_between_eofs(reference: eof.EOFData, target=eof.EOFData):
    """
    Calculates angle between two EOF vectors to determine their "closeness."
    theta = arccos(t . r / (||r||*||t||)), 

    :param reference: The reference-EOFs. This is usually the EOF pair of the previous or "first" DOY.
    :param target: The EOF that you want to find the angle with

    :return: A tuple of the  the angles between the reference and target EOFs for both EOF1 and EOF2
    """

    angle1 = tools.angle_btwn_vectors(reference.eof1vector, target.eof1vector)
    angle2 = tools.angle_btwn_vectors(reference.eof2vector, target.eof2vector)

    return (angle1, angle2)


################### EOF Calculation with rotation

def calc_eofs_from_olr_with_rotation(olrdata: olr.OLRData, implementation: str = "internal", sign_doy1reference: bool = True,
                            strict_leap_year_treatment: bool = False, no_leap: bool = False) -> eof.EOFDataForAllDOYs:
    """
    Computes EOFs in a similar way to the function above, but rotates EOFs such that they are continuous and non-degenerate. 
    Rotation algorithm described in :meth:'post_process_rotation'

    This function executes consistently the preprocessing (filtering), the actual EOF analysis, and the postprocessing.

    :param olrdata: The OLR dataset, from which OMI should be calculated. Note that OLR values are assumed to be given
        in positive values. The spatial grid of the OLR datasets defines also the spatial grid of the complete OMI
        calculation.
    :param implementation: See :meth:`calc_eofs_from_preprocessed_olr` in the mjoindices package.
    :param sign_doy1reference: See :meth:`correct_spontaneous_sign_changes_in_eof_series` in the mjoindices package.
    :param strict_leap_year_treatment: See description in :meth:`mjoindices.tools.find_doy_ranges_in_dates` 
    in the mjoindices package.
    :param no_leap: If true, will act as if there are no leap years in the dataset (365 days consistently)

    :return:
    """

    # preprocess OLR data
    preprocessed_olr = omi.preprocess_olr(olrdata)
    # calculate EOFs from raw data
    raw_eofs = omi.calc_eofs_from_preprocessed_olr(preprocessed_olr, implementation=implementation, 
                                                strict_leap_year_treatment=strict_leap_year_treatment, no_leap=no_leap) 
    # postprocess data
    result = post_process_rotation(raw_eofs, sign_doy1reference=sign_doy1reference, 
                                    strict_leap_year_treatment=strict_leap_year_treatment, no_leap=no_leap)

    return result


def post_process_rotation(eofdata: eof.EOFDataForAllDOYs, sign_doy1reference: bool = True,
                            strict_leap_year_treatment: bool = False, no_leap: bool = False) -> eof.EOFDataForAllDOYs:
    """
    Post processes a series of EOF pairs for all DOYs.

    Postprocessing includes an alignment of EOF signs and a rotation algorithm that rotates the EOFs
    in three steps:
    1. Projects EOFs at DOY = n-1 onto EOF space for DOY = n. This is done to reduce spurious oscillations
    between EOFs on sequential days
    2. Rotate the projected EOFs by 1/366 (or 1/365) per day to ensure continuity across January to December
    3. Renormalize the EOFs to have a length of 1 (this is a small adjustment to account for small numerical
    errors).

    See documentation of  the methods :meth:`correct_spontaneous_sign_changes_in_eof_series` in the original 
    mjoindices package for EOF sign flipping

    Note that it is recommended to use the function :meth:`calc_eofs_from_olr_with_rotation` to cover the complete algorithm.

    :param eofdata: The EOF series, which should be post processed.
    :param sign_doy1reference: See description of :meth:`correct_spontaneous_sign_changes_in_eof_series` 
    in the mjoindices package.
    :param no_leap: if True, assumes all years have 365 days

    :return: the postprocessed series of EOFs
    """
    pp_eofs = omi.correct_spontaneous_sign_changes_in_eof_series(eofdata, doy1reference=sign_doy1reference, no_leap=no_leap)
    rot_eofs = rotate_eofs(pp_eofs, no_leap=no_leap)
    norm_eofs = normalize_eofs(rot_eofs, no_leap=no_leap)
    
    return norm_eofs


def rotate_eofs(orig_eofs: eof.EOFDataForAllDOYs, no_leap: bool = False) -> eof.EOFDataForAllDOYs:
    """
    Rotate EOFs at each DOY to 1) align with the EOFs of the previous day and 2) be continuous across December to
    January boundary. Described more in detail in :meth:'post_process_rotation'

    :param orig_eofs: calculated EOFs, signs have been changed via spontaneous_sign_changes
    :param no_leap: if True, assumes all years have 365 days

    :return: set of rotated EOFs
    """

    delta = calculate_angle_from_discontinuity(orig_eofs, no_leap)

    print('Rotating by ', delta)

    return rotate_each_eof_by_delta(orig_eofs, delta, no_leap)


def rotation_matrix(delta):
    """
    Return 2d rotation matrix for corresponding delta
    """
    return np.array([[np.cos(delta), -np.sin(delta)],[np.sin(delta), np.cos(delta)]])


def calculate_angle_from_discontinuity(orig_eofs: eof.EOFDataForAllDOYs, no_leap: bool):
    """
    Project the matrix to align with previous day's EOFs and calculate the resulting
    discontinuity between January 1 and December 31. Divide by number of days in year to 
    result in delta for rotation matrix. 

    :param orig_eofs: calculated EOFs, signs have been changed via spontaneous_sign_changes
    :param no_leap: if True, assumes all years have 365 days

    :return: float of (negative) average angular discontinuity between EOF1 and EOF2 on the 
    first and last day of year, divided by the length of the year.
    """

    list_of_doys = tools.doy_list(no_leap)
    doy1 = orig_eofs.eofdata_for_doy(1)

    ndoys = len(list_of_doys)
    
    # set DOY1 initialization
    rots = np.array([doy1.eof1vector, doy1.eof2vector])

    # project onto previous day
    for d in list_of_doys:
        if d+1 > ndoys: # for last day in cycle, return to January 1
            doyn = orig_eofs.eofdata_for_doy(1)
        else:
            doyn = orig_eofs.eofdata_for_doy(d+1)

        B = np.array([doyn.eof1vector, doyn.eof2vector]).T 
        A = np.array([rots[0,:], rots[1,:]]).T
    
        rots = np.matmul(np.matmul(B, B.T),A).T
    
    # calculate discontinuity between Jan 1 and Jan 1 at end of rotation cycle
    discont = tools.angle_btwn_vectors(doy1.eof1vector, rots[0,:])

    return -discont/ndoys


def rotate_each_eof_by_delta(orig_eofs: eof.EOFDataForAllDOYs, 
                                delta: float, no_leap: bool) -> eof.EOFDataForAllDOYs:
    """
    Use delta calculated by optimization function to rotate original EOFs by delta.
    First projects EOFs from DOY n-1 onto EOF space for DOY n, then rotates projected
    EOFs by small angle delta. 

    :param orig_eofs: calculated EOFs, signs have been changed via spontaneous_sign_changes
    :param delta: scalar by which to rotate EOFs calculated from discontinuity
    :param no_leap: True if each year has 365 days, False if leap years included in dataset

    :returns: new EOFdata with rotated EOFs.  
    """

    R = rotation_matrix(delta)

    doy1 = orig_eofs.eofdata_for_doy(1)
    list_of_doys = tools.doy_list(no_leap)
    eofdata_rotated = []
    eofdata_rotated.append(doy1) # first doy is unchanged

    # set DOY1 initialization
    rots = np.array([doy1.eof1vector, doy1.eof2vector])

    # project onto previous day and rotate 
    for d in list_of_doys[1:]:
        doyn = orig_eofs.eofdata_for_doy(d)

        B = np.array([doyn.eof1vector, doyn.eof2vector]).T 
        A = np.array([rots[0,:], rots[1,:]]).T
    
        rots = np.matmul(np.matmul(np.matmul(B, B.T),A),R).T

        # create new EOFData variable for rotated EOFs
        eofdata_rotated.append(eof.EOFData(doyn.lat, doyn.long, 
                                np.squeeze(rots[0,:]), 
                                np.squeeze(rots[1,:]),
                                explained_variances=doyn.explained_variances,
                                eigenvalues=doyn.eigenvalues,
                                no_observations=doyn.no_observations))

    return eof.EOFDataForAllDOYs(eofdata_rotated, no_leap)


def normalize_eofs(orig_eofs: eof.EOFDataForAllDOYs, no_leap: bool) -> eof.EOFDataForAllDOYs:
    """
    :param eofdata: The rotated EOF series
    :param no_leap: if True, assumes all years have 365 days

    :return: normalize the EOFs to have length 1
    """
    list_of_doys = tools.doy_list(no_leap)

    eofdata_normalized = []

    for d in list_of_doys:

        doyn = orig_eofs.eofdata_for_doy(d)
        eof1_norm = doyn.eof1vector/np.linalg.norm(doyn.eof1vector)
        eof2_norm = doyn.eof2vector/np.linalg.norm(doyn.eof2vector) 

       # create new EOFData variable for rotated EOFs
        eofdata_normalized.append(eof.EOFData(doyn.lat, doyn.long, 
                                            np.squeeze(eof1_norm), 
                                            np.squeeze(eof2_norm),
                                            explained_variances=doyn.explained_variances,
                                            eigenvalues=doyn.eigenvalues,
                                            no_observations=doyn.no_observations)) 

    return eof.EOFDataForAllDOYs(eofdata_normalized, no_leap) 
