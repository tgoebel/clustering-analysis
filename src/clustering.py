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

#===============================================================================
#                          my modules
#===============================================================================
import data_utils

#===============================================================================
# 
#===============================================================================
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
    dConst      - {'Mc':float, 'b':float, 'D':float} # parameter dictionary
                   - completeness , b-value, fractal dimension
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
        vUncer = np.random.randn( eqCat.size())*1e-10
        eqCat.data['X']    += vUncer
        eqCat.data['Time'] += abs( vUncer)#time has to stay positive otherwise parent-offspring gets switched
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


def rescaled_t_r(catChild, catPar, dConst, **kwargs):
    """
    - compute rescaled time and distance

    Parameters
    ----------
    catChild, catPar - objects of type SeisCatDic containing parent and child events
    dConst      =  'b', 'D' -  b-value, fractal dimension
    kwargs       = 3D_distance = True default : False i.e. 2D Euclidean distance 
    Returns
    -------
    - { 'T' : np.array();  'R' : np.array(} }


    see: Clustering Analysis of Seismicity and Aftershock Identification, Zaliapin, I. (2008)

    """

    #-------------------------------set args and kwargs-----------------------------------------------
    M0 = 0
    if 'M0' in kwargs.keys() and kwargs['M0'] is not None:
        M0 = kwargs['M0']
    #-----------------------------add small uncertainty to X in case events are colocated-------------------------- 
    if 'correct_co_located' in kwargs.keys() and kwargs['correct_co_located'] == True:
        vUncer = np.random.randn( catChild.size())*1e-10
        catChild.data['X']    += vUncer
        catChild.data['Time'] += abs( vUncer) # time stays positive so parent offspring is not switched
    #------------------------------------------------------------------------------         
    #vMagCorr = 10**(-0.5*dConst['b']*(catPar.data['MAG']-M0) )
    vMagCorr = 10**(-0.5*dConst['b']*(catPar.data['Mag']-M0) )
    if '3D_distance' in kwargs.keys() and kwargs['3d_distance'] == True:
        vR       = np.sqrt( (catChild.data['X']-catPar.data['X'])**2 + (catChild.data['Y']-catPar.data['Y'])**2+ (catChild.data['Z']-catPar.data['Z'])**2 )**dConst['D']*vMagCorr
    else:
        vR       = np.sqrt( (catChild.data['X']-catPar.data['X'])**2 + (catChild.data['Y']-catPar.data['Y'])**2 )**dConst['D']*vMagCorr
    return   vR,  (catChild.data['Time']-catPar.data['Time'])*vMagCorr

def assembleClusters_NND( vID_child, vID_parent, vSim, simThreshold, **kwargs):
    """
    !! assuming parent and off-spring is connected via unique measurement (e.g. nearest-neighbor distance)
    - create clusters of event pairs based on some similarity criteria
            e.g. a) based on cross-correlation coefficients between pairs
                 b) based on space-time-magnitude distance
    - main input are pairs of connected events separated in parent and offspring
                  (one parent can have many children, but child has only one parent)
    1) find initial singles beyond threshold
    2) find pairs below threshold and assemble clusters
        - take all event pairs with values below (eta_0) or above (CCC),
           --> pairs beyond the threshold do not have to be considered
            if offspring meets similarity criteria:
                - go through each pair and find cluster for child event by searching if the
                  corresponding ID is already in any of the previous clusters
                -  attach to existing cluster or create new cluster
    3) - check if several offspring are connected to same parent and if 
         different clusters have to be combined in case of ID repition 
         --> this is implemented as a while loop
    4) - remove potential multiple IDs from clusters

    :Input    - simThreshold = similarity parameter
              - vID_parent   - event IDs
              - vID_child
              - vSimValues   - all similarity values
              kwargs['useLargerEvents'] = False, 

    :Return  dClust - python dictionary that contains all clusters labeled numerically
                     from '0' - not clustered
                          '1' - '[nCLmax]' - clustered events
                      each dictionary column contains IDs of children [first row] and parents [second row]
    """
    dNND = { 'aEqID_c' : vID_child, 
             'aEqID_p' : vID_parent,
             'aNND'    : vSim}
    # remove identical parents and off-spring if eq is in catalog several times
    sel = abs(vID_child -vID_parent) > 0
    dNND= data_utils.selDicAll(dNND, sel)

    #==================================1=============================================
    #                  initial  selection of events beyond threshold (single event)
    #================================================================================
    ### events without trigger
    if 'useLargerEvents' in kwargs.keys() and kwargs['useLargerEvents'] == True:
        print 'assuming threshold (%s) is a MINIMUM, select similarity values ABOVE this threshold'%( simThreshold)
        sel_single     = dNND['aNND'] <= simThreshold
        # remove independent events
        dNND_trig = data_utils.selectDataRange( dNND, simThreshold, None, 'aNND')
    else:
        print 'assuming threshold (%s) is a MAXIMUM, select similarity values BELOW this threshold'%( simThreshold)
        sel_single     = dNND['aNND'] >= simThreshold
        # remove independent events
        dNND_trig = data_utils.selDicAll( dNND, np.logical_not( sel_single))
    # preliminary single selection with eta > eta_0, may contain cluster events
    vID_single  = dNND['aEqID_c'][sel_single] # could be singles or parents but not offspring
    print 'intial N events not clustered', sel_single.sum(),
    print 'initial N events clustered', dNND_trig['aEqID_c'].shape[0], dNND_trig['aEqID_p'].shape[0], 'N-tot', vID_parent.shape[0]
    #==================================2=============================================
    #                      find clustered events
    #================================================================================
    # initiate vectors and dic during first run
    curr_child_ID     = dNND_trig['aEqID_c'][0]
    curr_par_ID       = dNND_trig['aEqID_p'][0]
    v_pastEqIDs = np.array(  [curr_child_ID, curr_par_ID] )
    v_pastClIDs = np.array(  [1, 1] )
    # dClust['0'] = singles
    dClust = {  '1'     : np.array( [[curr_child_ID],
                                     [curr_par_ID  ] ])}
    # for each child find the corresponding parent ID
    # if child or parent ID are already part of a cluster append to this cluster
    nCl = 2
    for iEv in xrange(1, dNND_trig['aEqID_p'].shape[0]):
        #print 'nPair', iEv+1, 'out of', len( dNND_trig['aEqID_p']), 'iCl', nCl
        curr_child_ID     = dNND_trig['aEqID_c'][iEv]
        curr_par_ID       = dNND_trig['aEqID_p'][iEv]
        # check if parent or child are part of previous cluster
        sel_child = curr_child_ID == v_pastEqIDs
        sel_par   = curr_par_ID   == v_pastEqIDs

        if sel_par.sum() > 0 or sel_child.sum() > 0:
            # find which cluster event pair belongs to
            if sel_par.sum() and sel_child.sum(): # both already part of a cluster
                curr_cl_ID1 = v_pastClIDs[sel_par][0]
                curr_cl_ID2 = v_pastClIDs[sel_child][0]
                # merge clusters and add IDs
                dClust[str(curr_cl_ID1)] =    np.hstack( (dClust[str(curr_cl_ID1)],
                                                             np.array([[curr_child_ID], [curr_par_ID  ] ])
                                                             ))
                dClust[str(curr_cl_ID1)] = np.hstack( (dClust[str(curr_cl_ID1)], dClust[str(curr_cl_ID2)]))
                # add new events but previous cluster ID
                v_pastEqIDs = np.append(  v_pastEqIDs, np.array([curr_child_ID, curr_par_ID] ) )
                v_pastClIDs = np.append(  v_pastClIDs, np.array([   curr_cl_ID1, curr_cl_ID1] ) )
                # remove second cluster ID from dClust
                dClust.pop( str(curr_cl_ID2))
                # remove from past eq IDs and past cl IDs
                sel = curr_cl_ID2 != v_pastClIDs
                v_pastEqIDs = v_pastEqIDs[sel]
                v_pastClIDs = v_pastClIDs[sel]
            else: # only one is part of a cluster
                if sel_par.sum() > 0: # parent already part of a cluster
                    curr_cl_ID = v_pastClIDs[sel_par][0]
                else:# child already part of a cluster
                    curr_cl_ID = v_pastClIDs[sel_child][0]
                dClust[str(curr_cl_ID)] =    np.hstack( (dClust[str(curr_cl_ID)],
                                                             np.array([[curr_child_ID], [curr_par_ID  ] ])
                                                             ))
                v_pastEqIDs = np.append(  v_pastEqIDs, np.array([curr_child_ID, curr_par_ID] ) )
                v_pastClIDs = np.append(  v_pastClIDs, np.array([   curr_cl_ID, curr_cl_ID        ] ) )
        else: # start a new cluster
            dClust[str(nCl)] =    np.array( [[curr_child_ID],
                                             [curr_par_ID  ] ])
            v_pastEqIDs = np.append(  v_pastEqIDs, np.array([curr_child_ID, curr_par_ID] ) )
            v_pastClIDs = np.append(  v_pastClIDs, np.array([          nCl, nCl        ] ) )
            nCl += 1
    # check if children have same parent
    nTotChild = 0
    #=================================3==========================================================================
    #                 remove events from singles if in cluster, remove multiple IDs
    #============================================================================================================
    # create vector of triggered eqIDs and count triggered events
    vID_Trig_all = np.array([])
    vclID_allEv  = np.array([], dtype = int)
    for tag in sorted( dClust.keys()):
        #print 'iCl', tag, 'nEv in cluster', np.unique( dClust[tag].flatten()).shape[0]
        #print dClust[tag][0]
        aID_flat_uni = np.unique( dClust[tag].flatten())
        #nTotTrig   += aID_flat_uni.shape[0]
        vID_Trig_all = np.append( vID_Trig_all, aID_flat_uni )
        vclID_allEv  = np.append( vclID_allEv, np.ones( aID_flat_uni.shape[0], dtype = int)*int(tag))
        # remove multiple ID entries --> possible since pairs are always appeneded
        dClust[tag] = aID_flat_uni
        nTotChild  += dClust[tag].shape[0]-1
    #====================================4========================================================================
    #                       check for events in more than one cluster
    #============================================================================================================
    aIDs, aCounts = np.unique( vID_Trig_all, return_counts=True)
    selDouble = aCounts > 1
    iD = 1
    while selDouble.sum() > 0:
        print '%i run to remove doubles'%(iD)
        for ID in np.unique( aIDs[selDouble]):
            selCl = ID == vID_Trig_all
            aClID = np.unique( vclID_allEv[selCl])
            for iCl in xrange( len( aClID)-1):
                print 'iCl with same events', str( aClID[0]), str( aClID[iCl+1]),
                print ID
                #A# merge clusters that have same events
                dClust[str(aClID[0])] = np.unique( np.hstack( (dClust[str(int( aClID[0]))], dClust[str( int(aClID[iCl+1]))])))
                #B# remove cluster IDs from dictionary
                dClust.pop( str( int( aClID[iCl+1])))
            #C# remove double event and corresponding clID from:
            # vID_Trig_all
            sel_rem = ID != vID_Trig_all
            vID_Trig_all = vID_Trig_all[sel_rem]
            # and vclID_allEv
            vclID_allEv  = vclID_allEv[sel_rem]
            # leave one  event with new clID, i.e. clId of first cluster that contains ID
            vclID_allEv  = np.append( vclID_allEv, aClID[0])
            vID_Trig_all= np.append( vID_Trig_all, ID)
        aIDs, aCounts = np.unique( vID_Trig_all, return_counts=True)
        selDouble = aCounts > 1
        iD+=1
    # find events within initial single selection (eta > eta_0)
    # which are actually part of clustered events
    sel_single = np.ones( vID_single.shape[0], dtype = int) > 0
    iS = 0
    for ID_single in  vID_single:
        sel = ID_single == vID_Trig_all
        if sel.sum() > 0: # remove this event from singles
            sel_single[iS] = False
        iS += 1
    vID_single = vID_single[sel_single]
    print 'Ntot in cluster ', len( vID_Trig_all), 'N-parent',len(dClust.keys()),  'N singles', vID_single.shape[0], 'Ntot. offspring (includes doubles)', nTotChild
    print 'Ntot in cat.', vSim.shape[0], 'N-trig + N-ind', len( vID_Trig_all)+vID_single.shape[0]
    dClust[str(0)] =  vID_single
    return dClust
