import pandas as pd
import numpy as np

import mjoindices.empirical_orthogonal_functions as eof




def organize_pc_dataframe(pc_path) -> pd.DataFrame:
    """
    Take in PC files for as calculated by calc_OMI_PCs.py
    Restructure the data into a dataframe with dates split into YY, MM, DD
    Calculate the OMI amplitude (RMS of PC1 and PC2)

    :param years: which years to take PCs from (outdated, now all years in same file)
    :param pc_path: path that leads to PC text file / filename

    :returns: Pandas dataframe of PC components
    """

    pcs = pd.read_csv(pc_path, sep=',', header=0)

    # calculate the amplitde for each day (RMS of PC1 and PC2)
    pcs['Amplitude'] = np.sqrt(pcs.PC1**2 + pcs.PC2**2)

    # expand date into YR, MM, DD variables
    yy = [None]*len(pcs.Date)
    mm = [None]*len(pcs.Date)
    dd = [None]*len(pcs.Date)

    for i in range(len(pcs.Date)):
        yy[i],mm[i],dd[i] = map(int, pcs.Date[i].split('-'))

    pcs['Year'] = yy
    pcs['Month'] = mm
    pcs['Day'] = dd

    return pcs



def calc_mean_eof(eofdata, start_doy, end_doy):
    
    """
    Calculate mean EOF at each gridpoint for the time period specified by start_doy and end_doy.
    Used to create EOF anomaly for Figure 3. 
    """
    
    doy1 = eofdata.eofdata_for_doy(1)
    eof_map_grid = np.empty([doy1.eof1map.shape[0], doy1.eof1map.shape[1], 2, end_doy-start_doy])
    
    doy_range = [d for d in range(start_doy, end_doy)]
    
    # select EOFs within DOY range
    for i, d in enumerate(doy_range):
        
        doyn = eofdata.eofdata_for_doy(d)
        eof_map_grid[:,:,0,i] = doyn.eof1map
        eof_map_grid[:,:,1,i] = doyn.eof2map
        
    # calculate average
    eof_map_avg = np.mean(eof_map_grid, axis=3)
    
    return eof_map_avg