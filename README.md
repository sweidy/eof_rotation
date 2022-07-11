# eof_rotation

This code is an add-on to the OMI calculation package published in Hoffmann, C.G., Kiladis, G.N., Gehne, M. and von Savigny, C., 2021. 
A Python Package to Calculate the OLR-Based Index of the Madden- Julian-Oscillation (OMI) in Climate Science and Weather Forecasting. 
Journal of Open Research Software, 9(1), p.9. DOI: http://doi.org/10.5334/jors.331
which can be found on github in at https://github.com/cghoffmann/mjoindices. 

The complete OMI algorithm is described in Kiladis, G.N., J. Dias, K.H. Straub, M.C. Wheeler, S.N. Tulich, K. Kikuchi, K.M.
Weickmann, and M.J. Ventrice, 2014: A Comparison of OLR and Circulation-Based Indices for Tracking the MJO.
Mon. Wea. Rev., 142, 1697–1715, https://doi.org/10.1175/MWR-D-13-00301.1

A description of the rotation algorithm can be found in the preprint of Sarah Weidman, Ned Kleiner, and Zhiming Kuang. Submitted. “A rotation procedure to improve seasonally varying Empirical Orthogonal Function bases for MJO indices.” Geophysical Research Letters, https://doi.org/10.1002/essoar.10511626.1

The function :meth:'calc_eofs_from_olr_with_rotation' included in this file can be used instead of the function 
:meth:'calc_eofs_from_olr' in the original mjoindices package. The new method includes a projection and rotation postprocessing
step that reduces noise in the original EOF calculation. 

Postprocessing includes an alignment of EOF signs and a rotation algorithm that rotates the EOFs in three steps:

    1. Projects EOFs at DOY = n-1 onto EOF space for DOY = n. This is done to reduce spurious oscillations
    between EOFs on sequential days
    2. Rotate the projected EOFs by 1/366 (or 1/365) per day to ensure continuity across January to December
    3. Renormalize the EOFs to have a length of 1 (this is a small adjustment to account for small numerical
    errors).
    
Figures from the EOF rotation paper can be reproduced using the eof_rotation_example.ipynb file, along with the functions in rotation_plotting_tools.py. An example of how to implement the rotation algorithm using the above packages can be found in either eof_rotation_example.ipynb or eof_rotation_example.py. 

Adaptation written by Sarah Weidman, 2022. Contact: sweidman@g.harvard.edu
