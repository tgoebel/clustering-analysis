# -*- coding: utf-8 -*-
# python3.6
'''
Created on April 10th, 2019

    - function required for clustering analysis based on nearest-neighbor distances
    
    - NND_eta - eq. 1 for NND in Zaliaping & Ben-Zion 2013

@author: tgoebel
'''
from __future__ import division
import numpy as np




def NND_eta( eqCat, dConst, **kwargs):
    """
        - NND_eta - eq. 1 for NND in Zaliaping & Ben-Zion 2013
    search for 'parent event' i.e. earthquake that occurred closest in space-time-magnitude domain 
                                   but prior to the current event
        here: [jC]          - are the child events and we try to find the one closest parent, occurring earlier in time
              [sel_tau_par] - are the potential parent evnets that occurred before [jC], we select the closest in time

    Parameters
    ----------
    catalog     - catalog.data['Time'], 'X', 'Y', 'Depth', 'MAG'
                - time, cartesian coordinates (X,Y, Depth), magnitude
    dConst      =  'Mc', 'b', 'D' - completeness , b-value, fractal dimension
    kwargs      - rmax (default: = 200) - maximum space window (for faster computation)
                - tmax (default: =  10) - maximum time window (for faster computation)
                - correct_co_located = True, add gaussian uncertainty to avoid product going to zero for co-located earthquakes
                - haversine = True - use haversine distance at surface instead of 3D cartesian distance
                - M0 - reference magnitude, default: M0 = 0
    Returns
    -------
    - {  'aNND'       : aNND,     - nearest neighbor space-time magnitude distance
         'aEqID_p'    : np.array  - ID of the parent event
         'aEqID_c'    : np.array  - ID of the child  event
        } 

    see: Clustering Analysis of Seismicity and Aftershock Identification, Zaliapin, I. (2008)

    """
    #-------------------------------set args and kwargs----------------------------------------------- 
    rmax = 500 # in km
    tmax = 20 # in years
    M0 = 0
    if 'M0' in kwargs.keys() and kwargs['M0'] is not None:
        M0 = kwargs['M0']
    if 'rmax' in kwargs.keys() and kwargs['rmax'] is not None:
        rmax = kwargs['rmax']
    if 'tmax' in kwargs.keys() and kwargs['tmax'] is not None:
        tmax = kwargs['tmax']
    #-----------------------------add small uncertainty to X in case events are colocated-------------------------- 
    if 'correct_co_located' in kwargs.keys() and kwargs['correct_co_located'] == True:
        vUncer = np.random.randn( eqCat.size())*1e-12 
        eqCat.data['X']    += vUncer
        eqCat.data['Time'] += vUncer
    #------------------------------------------------------------------------------         
    #mEta     = np.zeros( (eqCat.size(), eqCat.size()), dtype = float)
    #aNND     = np.ones(  eqCat.size())
    aNND     = np.zeros(  eqCat.size())
    vID_p    = np.zeros( eqCat.size())
    vID_c    = np.zeros( eqCat.size())
    deltaMag = (eqCat.data['Mag'] - M0)
 
    for jC in range( eqCat.size()):
        print 'event %i of %i'%( jC+1, eqCat.size())
        # interevent times: take events that happend before t_i 
        #           child             - parent                > 0 
        tau         =  eqCat.data['Time'][jC] - eqCat.data['Time']
        sel_tau_par = tau > 0
        if sel_tau_par.sum() > 0:

            vcurr_ID = np.arange( eqCat.size(), dtype = int)[sel_tau_par]
            vR = np.sqrt( (eqCat.data['X'][jC] - eqCat.data['X'][vcurr_ID])**2 + (eqCat.data['Y'][jC] - eqCat.data['Y'][vcurr_ID])**2 )
            #  haversine distance
            # = projUtils.haversine( eqCat.data['Lon'][jC], eqCat.data['Lat'][jC],eqCat.data['Lon'][curr_vID], eqCat.data['Lat'][curr_vID] ) 
            sel_r_par = vR < rmax
            if sel_r_par.sum() > 0:
                vcurr_ID = vcurr_ID[sel_r_par]
                curr_Eta = tau[vcurr_ID]* (vR[sel_r_par]**dConst['D']) *( 10**(-dConst['b']*deltaMag[vcurr_ID]))
                sel_min  = curr_Eta == curr_Eta.min()
                aNND[jC]    = curr_Eta[sel_min][0]
                vID_p[jC]   = eqCat.data['N'][vcurr_ID][sel_min][0]
                vID_c[jC]   = eqCat.data['N'][jC]
                #print 'parent', eqCat.data['N'][vcurr_ID][sel_min][0],  'offspring', eqCat.data['N'][jC]
                #print 'parent', eqCat.data['Time'][vcurr_ID][sel_min][0],  'offspring', eqCat.data['Time'][jC]

                if sel_min.sum() > 1:
                    print aNND[jC], curr_Eta[sel_min], eqCat.data['N'][vcurr_ID][sel_min]
                    print eqCat.data['Lon'][vcurr_ID][sel_min], eqCat.data['Lat'][vcurr_ID][sel_min]
                    print eqCat.data['X'][vcurr_ID][sel_min], eqCat.data['Y'][vcurr_ID][sel_min]
        # else:
        #     aNND[jC]    = 1e-6
        #     vID_p[jC]   = eqCat.data['N'][jC]#[sel_min][0]
        #     vID_c[jC]   = eqCat.data['N'][jC]
    #mEta[mEta<0] = 0
    # remove events with aNND < 0; i.e. event at the beginning with no preceding parent
    return {  'aNND' : aNND[aNND>0], 'aEqID_p' : vID_p[aNND>0], 'aEqID_c' : vID_c[aNND>0]}
