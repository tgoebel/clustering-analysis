# -*- coding: utf-8 -*-
# python3.7
'''
Created on April 10th, 2019

    - function required for clustering analysis based on nearest-neighbor distances
    
    - NND_eta - eq. 1 for NND in Zaliaping & Ben-Zion 2013

@author: tgoebel
'''
import numpy as np
import matplotlib.pyplot as plt
import warnings
#===============================================================================
#                          my modules
#===============================================================================
import src.data_utils as data_utils

#===============================================================================
# 
#===============================================================================
def NND_eta( eqCat, dConst, verbose = False, **kwargs):
    """
        - NND_eta - eq. 1 for NND in Zaliapin & Ben-Zion 2013
    search for 'parent event' i.e. earthquake that occurred closest in space-time-magnitude domain 
                                   but prior to the current event
        here: [jC]          - are the off spring events and we try to find the closest parent, occurring earlier in time
              [sel_tau_par] - are the potential parent events that occurred before [jC], we select the closest in time

    Parameters
    ----------
    catalog     - catalog.data['Time'], 'Lon', 'Lat' (or 'X', 'Y',) 'Depth', 'MAG'
                - time, cartesian coordinates (X,Y, Depth), magnitude
    dConst      - {'Mc':float, 'b':float, 'D':float} #  dictionary with statistical seismicity parameters
                   - completeness , b-value, fractal dimension
    kwargs      - rmax (default: = 500) - maximum space window (for faster computation)
                - tmax (default: =  20) - maximum time window (for faster computation)
                - correct_co_located = True, add gaussian uncertainty to avoid product going to zero for co-located earthquakes
                - haversine = True - use haversine distance at surface instead of 3D cartesian distance
                - M0 - reference magnitude, default: M0 = 0
    Returns
    -------
    - {  'aNND'       : aNND,     - nearest neighbor space-time magnitude distance
         'aEqID_p'    : np.array  - ID of the parent event
         'aEqID_c'    : np.array  - ID of the child  event
         'Time'       : np.array  - origin time of offspring
        } 

    see: Clustering Analysis of Seismicity and Aftershock Identification, Zaliapin, I. (2008)

    """
    #-------------------------------set args and kwargs----------------------------------------------- 
    rmax = 500 # in km
    tmax = 20 # in years
    M0 = 0 # reference mag
    if 'M0' in kwargs.keys() and kwargs['M0'] is not None:
        M0 = kwargs['M0']
    if 'rmax' in kwargs.keys() and kwargs['rmax'] is not None:
        rmax = kwargs['rmax']
    if 'tmax' in kwargs.keys() and kwargs['tmax'] is not None:
        tmax = kwargs['tmax']
    #-----------------------------add small uncertainty to X in case events are colocated-------------------------- 
    if 'correct_co_located' in kwargs.keys() and kwargs['correct_co_located'] == True:
        vUncer = np.random.randn( eqCat.size())*1e-10
        eqCat.data['Lon']    += vUncer
    #------------------------------------------------------------------------------
    aNND     = np.zeros( eqCat.size())
    vID_p    = np.zeros( eqCat.size())
    vID_c    = np.zeros( eqCat.size())
    a_M_MS_ref= (eqCat.data['Mag'] - M0)# mainshock mag with respect to reference
 
    for jC in range( eqCat.size()):
        if verbose == True:
            print( f"event {jC+1:d} of {eqCat.size():d}", end= "\r")
        # interevent times: take events that happend before t_i 
        #           child             - parent                > 0 
        tau         =  eqCat.data['Time'][jC] - eqCat.data['Time']
        sel_tau_par = tau > 0
        if sel_tau_par.sum() > 0:

            vcurr_ID = np.arange( eqCat.size(), dtype = int)[sel_tau_par]
            # if cartesian coordinates are available
            if 'X' in eqCat.data.keys() and 'Y' in eqCat.data.keys():
                vR = np.sqrt( (eqCat.data['X'][jC] - eqCat.data['X'][vcurr_ID])**2 + (eqCat.data['Y'][jC] - eqCat.data['Y'][vcurr_ID])**2 )
            else:
                #  haversine distance
                vR = haversine( eqCat.data['Lon'][jC], eqCat.data['Lat'][jC],eqCat.data['Lon'][vcurr_ID], eqCat.data['Lat'][vcurr_ID] )
            sel_r_par = vR < rmax
            if sel_r_par.sum() > 0:
                vcurr_ID = vcurr_ID[sel_r_par]
                curr_Eta = tau[vcurr_ID]* (vR[sel_r_par]**dConst['D']) *( 10**(-dConst['b']*a_M_MS_ref[vcurr_ID]))
                sel_min  = curr_Eta == curr_Eta.min()
                aNND[jC]    = curr_Eta[sel_min][0]
                vID_p[jC]   = eqCat.data['N'][vcurr_ID][sel_min][0]
                vID_c[jC]   = eqCat.data['N'][jC]
                #print( 'parent', eqCat.data['N'][vcurr_ID][sel_min][0],  'offspring', eqCat.data['N'][jC]
                #print( 'parent', eqCat.data['Time'][vcurr_ID][sel_min][0],  'offspring', eqCat.data['Time'][jC]

                if sel_min.sum() > 1:
                    print( aNND[jC], curr_Eta[sel_min], eqCat.data['N'][vcurr_ID][sel_min])
                    print( eqCat.data['Lon'][vcurr_ID][sel_min], eqCat.data['Lat'][vcurr_ID][sel_min])
    sel2 = aNND > 0
    if np.logical_not(sel2).sum() > 0:
        print( f"{np.logical_not(sel2).sum()} %i events with NND=0 ")
        #raise ValueError, error_str
    #  remove events with aNND < 0; i.e. event at the beginning with no preceding parent
    return {  'aNND' : aNND[sel2], 'aEqID_p' : vID_p[sel2], 'aEqID_c' : vID_c[sel2], 'Time' : eqCat.data['Time'][sel2]}
    #return {  'aNND' : aNND, 'aEqID_p' : vID_p, 'aEqID_c' : vID_c, 'Time' : eqCat.data['Time'][1::]}


def rFromTau( dt, b, D, eta_0, M_MS ):
    """
        - compute maximum distance R for events in cluster
          based on interevent time, eta_0 and D (fractal dimension)
    :INPUT
          dt    - array or float
               interevent times (dt relative to MS or first event in family)
          b     - Gutenberg-Richter b-value
          D     - fractal dimension, usually D~1.6
          eta_0 - empiricallly determined separation line between clustered and background
                  mode
          M_MS  - mainshock magnitude (here we assume only one triggering generation)
    :return:
    """
    return ( -eta_0/dt * 10**( b*M_MS))**(1/D)*1e-3

def rescaled_t_r(catChild, catPar, dConst, **kwargs):
    """
    - compute rescaled time and distance

    Parameters
    ----------
    catChild, catPar - objects of type SeisCatDic containing parent and child events
    dConst      =  'b', 'D' -  b-value, fractal dimension
    kwargs       = distance_3D = True default : False i.e. 2D Euclidean distance

    Returns
    -------
    - a_R, a_tau


    see: Clustering Analysis of Seismicity and Aftershock Identification, Zaliapin, I. (2008)

    """
    #-------------------------------set args and kwargs-----------------------------------------------
    M0 = 0
    if 'M0' in kwargs.keys() and kwargs['M0'] is not None:
        M0 = kwargs['M0']
    #-----------------------------add small uncertainty to X in case events are colocated-------------------------- 
    if 'correct_co_located' in kwargs.keys() and kwargs['correct_co_located'] == True:
        vUncer = np.random.randn( catChild.size())*1e-10
        catChild.data['Lon']    += vUncer
    #------------------------------------------------------------------------------         
    #vMagCorr = 10**(-0.5*dConst['b']*(catPar.data['MAG']-M0) )
    vMagCorr = 10**(-0.5*dConst['b']*(catPar.data['Mag']-M0) )
    # if cartesian coordinates are available
    if 'X' in catChild.data.keys() and 'X' in catPar.data.keys():
        a_R = np.sqrt((catChild.data['X'] - catPar.data['X']) ** 2 + (catChild.data['Y'] - catPar.data['Y']) ** 2) ** \
              dConst['D'] * vMagCorr

    else:
        a_R = haversine(catChild.data['Lon'], catChild.data['Lat'],
                   catPar.data['Lon'],   catPar.data['Lat'])**dConst['D']*vMagCorr

    a_dt = catChild.data['Time']-catPar.data['Time']#interevent times
    a_tau = (a_dt)*vMagCorr
    sel2 = a_tau < 0
    if sel2.sum() > 0:
        #print( catChild.data['N'][sel2])
        #print( catPar.data['N'][sel2])
        error_str = '%i parents occurred after offspring, check order of origin time in catChild, catPar'%(sel2.sum())
        raise( ValueError( error_str))
    return a_R, a_tau


def compileClust( dNND, simThreshold, verbose = True,  **kwargs):
    """
    assuming parent and off-spring is connected via unique measurement (e.g. nearest-neighbor distance)
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
         different clusters have to be combined in case of ID repetition
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
    # dNND = { 'aEqID_c' : vID_child,
    #          'aEqID_p' : vID_parent,
    #          'aNND'    : vSim}
    # remove identical parents and off-spring if eq is in catalog several times
    sel = abs(dNND['aEqID_c']-dNND['aEqID_p']) > 0
    dNND= data_utils.selDicAll(dNND, sel)

    # check that dNND is sorted by time
    if 'Time' in dNND.keys():
        i_sort = np.argsort( dNND['Time'])
        dNND   = data_utils.selDicAll(dNND, i_sort)
    else:
        error_str = "'Time' key missing, add offspring origin time to dNND"
        raise ValueError( error_str)
    #==================================1=============================================
    #                  initial  selection of events beyond threshold (single event)
    #================================================================================
    ### events without trigger
    if 'useLargerEvents' in kwargs.keys() and kwargs['useLargerEvents'] == True:
        print( 'assuming threshold (%s) is a MINIMUM, select similarity values ABOVE this threshold'%( simThreshold))
        sel_single     = dNND['aNND'] <= simThreshold
        # remove independent events
        dNND_trig = data_utils.selectDataRange( dNND, simThreshold, None, 'aNND')
    else:
        print( 'assuming threshold (%s) is a MAXIMUM, select similarity values BELOW this threshold'%( simThreshold))
        sel_single     = dNND['aNND'] >= simThreshold
        # remove independent events
        dNND_trig = data_utils.selDicAll( dNND, np.logical_not( sel_single))
    # preliminary single selection with eta > eta_0, may contain cluster events
    vID_single  = dNND['aEqID_c'][sel_single] # could be singles or parents but not offspring
    sel_first = np.in1d( dNND['aEqID_p'][0], vID_single)
    if dNND['aNND'][0] > simThreshold and sel_first.sum() == 0:
        vID_single = np.append(  dNND['aEqID_p'][0], vID_single)

    if verbose == True:
        print( f"---------compileClust - initial numbers:------")
        print(  f"No. singles: {vID_single.shape[0]}"),
        print(  f"No. triggered: {dNND_trig['aEqID_c'].shape[0]}, {dNND_trig['aEqID_p'].shape[0]},"
                f"No. tot. {dNND_trig['aEqID_p'].shape[0]} {sel_single.sum()+dNND_trig['aEqID_c'].shape[0]}")
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
    for iEv in range(1, dNND_trig['aEqID_p'].shape[0]):
        #print( 'nPair', iEv+1, 'out of', len( dNND_trig['aEqID_p']), 'iCl', nCl
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
        #print( 'iCl', tag, 'nEv in cluster', np.unique( dClust[tag].flatten()).shape[0]
        #print( dClust[tag][0]
        aID_flat_uni = np.unique( dClust[tag].flatten())
        #nTotTrig   += aID_flat_uni.shape[0]
        vID_Trig_all = np.append( vID_Trig_all, aID_flat_uni )
        vclID_allEv  = np.append( vclID_allEv, np.ones( aID_flat_uni.shape[0], dtype = int)*int(tag))
        # remove multiple ID entries --> possible since pairs are always appeneded
        dClust[tag] = aID_flat_uni
        nTotChild  += dClust[tag].shape[0]-1
    #====================================4========================================================================
    #                       check for events in more than one cluster, merge clusters
    #============================================================================================================
    # sel_same = np.in1d( vID_Trig_all, np.array([ 3049419,  9020431,  9172305,  9173365, 15332137]))
    # print( "events in trig_all before double remove: ", sel_same.sum(), vID_Trig_all[sel_same])
    aIDs, aCounts = np.unique( vID_Trig_all, return_counts=True)
    selDouble = aCounts > 1
    if verbose == True:
        print( f"N event IDs in more than one cluster: {selDouble.sum()}")
    i_run = 1
    while selDouble.sum() > 0:
        if verbose == True:
            print( '%i. run to remove doubles'%(i_run))
        for ID in np.unique( aIDs[selDouble]):
            selCl = ID == vID_Trig_all
            aClID = np.unique( vclID_allEv[selCl])
            for iCl in range( len( aClID)-1):
                if verbose == True:
                    print( 'iCl with same events', str( aClID[0]), str( aClID[iCl+1]), 'evID: ', int(ID))
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
            vID_Trig_all = np.append( vID_Trig_all, ID)
        aIDs, aCounts = np.unique( vID_Trig_all, return_counts=True)
        selDouble = aCounts > 1
        i_run += 1
    # find events within initial single selection (eta > eta_0)
    # which are actually part of clustered events
    sel_single = np.ones( vID_single.shape[0], dtype = int) > 0
    iS = 0
    for ID_single in  vID_single:
        sel = ID_single == vID_Trig_all
        if sel.sum() > 0: # remove this event from singles
            sel_single[iS] = False
        iS += 1
    if verbose == True:
        print("initial singles now parents - remove from dClust['0']: ",np.array([~sel_single]).sum())

    vID_single = vID_single[sel_single]
    if verbose == True:
        print( "---------------final result--------------------------")
        print(  f" Ntot in cluster: {len( vID_Trig_all)}, N-parent(=N-clust): {len(dClust.keys())},"
            f"No. singles: {vID_single.shape[0]}, Ntot. offspring (includes doubles):  {nTotChild}")
        print( "trig. fraction: ", round((len( vID_Trig_all)-len(dClust.keys()))/dNND['aNND'].shape[0],2), "frac.MS: ", round( len(dClust.keys())/dNND['aNND'].shape[0],2), "single: ", round((vID_single.shape[0]/dNND['aNND'].shape[0]),2))
        print(  'Ntot in cat.', dNND['aNND'].shape[0]+1, 'N-trig + N-ind', len( vID_Trig_all)+vID_single.shape[0])

    dClust[str(0)] =  vID_single
    return dClust

def addClID2cat( seisCat, dClust, test_plot = False, **kwargs):
    """
    - add new column (i.e. dictionary tag='famID') for seisCat
        that specifies which cluster each event belongs to
    - !note that if offspring generation should be recorded run:
      clustering.offspring_gen() first
      and use output dictionary as input for this fct.

    :param dClust:  python dictionary
                    each dic. element specified by key is a vector of evIDs
                    or three row matrix with evID, iGen and average leaf depth


    :param seisCat:
    :return: seisCat (with new tags:
                        'famID' - record family links between events
                        optional:
                        (note that 'clID' is commonly used for waveform-based relocations)
                        'iGen'  - record offspring generation within family)
                        'LD'    - average lead depth for each cluster

    """
    # sort original catalog to get ID of first event
    seisCat.sortCatalog( 'Time')

    # first row is clusterID and second row event ID from catalog
    nRows  = 2
    b_add_iGen = False
    if len( dClust['0'].shape) > 1:
        b_add_iGen = True
        # additional rows for trig. generation and average lead depth
        nRows = 4
    mClust = np.zeros([nRows, seisCat.size()])
    nGen = 0
    nFam = 0
    nEv = 0
    i   = 0
    for sCl in dClust.keys():
        iCl = int(sCl)
        # print( f"------iCl: {iCl}, nEv: {nEv}--------, evID={dClust[sCl]}")
        #earthquake event IDs
        if b_add_iGen == False:
            nEv = dClust[sCl].shape[0]
            mClust[1, i:i + nEv] = dClust[sCl]
        else:
            nEv = dClust[sCl].shape[1]
            mClust[1, i:i + nEv] = dClust[sCl][0]
        # family IDS
        mClust[0,i:i+nEv] = np.ones(nEv)*iCl
        nFam += len( dClust[sCl])
        if b_add_iGen == True:
            nGen += len( dClust[sCl][1])
            # trig generation
            mClust[2,i:i+nEv] = dClust[sCl][1]
            # ave. lead depth
            mClust[3, i:i + nEv] = dClust[sCl][2]
        i += nEv
    #---------include first event in catalog as single-------------
    selFirst = seisCat.data['N'][0] == mClust[1]
    if selFirst.sum() == 0:
        ID_first = int( seisCat.data['N'][0]) #[~selUni][0])
        print( 'first ev. in catalog -ID:', ID_first, int( seisCat.data['N'][0]), 'last ev. in mClust', mClust[1,-1], 'should=0')
        mClust[1] = np.hstack( (ID_first, mClust[1,0:-1]))# not needed if catalog is sorted by Time
    #sel_same = np.in1d( mClust[1], seisCat.data['N'])
    # check that every event ID is represented only once
    __, aID, aN_uni = np.unique( mClust[1], return_counts = True, return_index=True)
    sel = aN_uni > 1
    if sel.sum() > 0:
        error_str = f"ev. ID represented more than once: {mClust[1][aID[sel]]}, 'N-repeats: ', {aN_uni[sel]}"
        raise ValueError( error_str)
    #--sort both cluster ID matrix and cat with respect to IDs
    sortSel = mClust[1].argsort()
    mClust = mClust.T[sortSel].T
    seisCat.sortCatalog('N') #--otherwise clIDs get assigned to wrong event

    if test_plot == True:
        plt.figure()
        plt.subplot( 211)
        plt.plot( mClust[1], mClust[1]-seisCat.data['N'], 'ko')
        plt.xlabel( 'Event ID in Clust')
        plt.ylabel( 'Diff. Events IDs (0)')
        plt.subplot( 212)
        plt.plot(mClust[1],  mClust[0], 'ko')
        plt.xlabel('Event ID in Clust')
        plt.ylabel('Cluster ID')
        #plt.plot( plt.gca().get_xlim(), plt.gca().get_xlim(), 'r--')
        plt.show()

    seisCat.data['famID'] = np.int32( mClust[0])
    if b_add_iGen == True:
        seisCat.data['iGen'] = np.int16(mClust[2])
    seisCat.sortCatalog( 'Time')
    return seisCat

def offspring_gen( dClust, dNND, f_eta_0, **kwargs):
    """
    - trace back triggering chain chronologically and assign trig generation
    - start with parent generation, then add end leafs
            a) identify all parents within cluster
            b) sort by time
            c) assign the same iGen to offspring of the same parent (hierarchical)
            compute average leaf depth:
                <d> = 1/n sum( d_i) = ave. depth across end leafs
    __________________________________
    input:  seisCat   = object SeismicityCatalog
                        used to get origin times of offspring events
            dNND =
            'aEqID_c'  - unique event IDs of offspring
            'aEqID_p ' - events IDs of parents, these are paired to a_ID_child so order matters
                          parents can have many offspring, so repeats are possible here
            'Time'     - offspring origin time from catalog, in case IDs are not chronological

            dClust - '[famID]' = np.array([ offSpringIDs])
    ----------------------------------
    return:
            dGen    - python dictionary
                    'famID' : np.array([3, N])
                    # dGen[famID][0] = evIDs
                    # dGen[famID][1] = trig generation
                    # dGen[famID][2] = ave. leaf depth
                             - average lead depth (same number for entire cluster)
    """
    #=========================1========================================
    #            count generations of offspring events
    #==================================================================
    dGen = {}
    l_famID = list( dClust.keys())
    # singles are all 0 generation
    dGen['0'] = np.zeros( (3, len( dClust['0'])))
    # set ev IDs in new dic
    dGen['0'][0] = dClust['0']

    # ave LD = 1
    dGen['0'][2] = np.ones( len( dClust['0']))

    # ignore singles below
    l_famID.remove( '0')
    for famID in l_famID:
        ###find ori. time for each child
        sel_chi_t  = np.in1d(  dNND['aEqID_c'], dClust[famID])
        # filter for parent - child NND < eta_0
        sel_chi_t2 = dNND['aNND'][sel_chi_t] < f_eta_0
        curr_iPar  = dNND['aEqID_p'][sel_chi_t][sel_chi_t2]
        curr_iChi  = dNND['aEqID_c'][sel_chi_t][sel_chi_t2]
        curr_tChi  = dNND['Time'][np.in1d( dNND['aEqID_c'],curr_iChi)]

        ##sort cluster IDs with respect to offspring time
        sel_sort  = np.argsort( curr_tChi)
        first_ID  = curr_iChi[sel_sort][0]

        # get unique parents and sort by time
        uni_curr_iPar = np.unique(curr_iPar)
        uni_par_times = curr_tChi[np.in1d(curr_iChi, uni_curr_iPar)]
        uni_curr_iPar = curr_iChi[np.in1d(curr_iChi, uni_curr_iPar)]
        sort_uni_par  = np.argsort( uni_par_times)
        uni_curr_iPar = uni_curr_iPar[sort_uni_par]
        # check if parent of first pair needs to be added
        if np.isin( curr_iPar[0], uni_curr_iPar).sum() == 0:
            uni_curr_iPar = np.hstack(( curr_iPar[0], uni_curr_iPar))
        # add end leafs (offspring that are not parents)
        sel_endLeaf   = ~np.in1d(curr_iChi, curr_iPar)
        uni_curr_iPar = np.hstack((uni_curr_iPar, curr_iChi[sel_endLeaf]))
        #----------initiate new vectors------------------------
        uni_iGen_pastPar = np.zeros( len(curr_tChi)+1)
        uni_id_pastPar   = np.zeros( len(curr_tChi)+1)
        ## assign chronological triggering generation
        curr_iGen        = np.zeros( len(curr_tChi)+1)
        iGen          = 0
        for iPar in range( len(uni_curr_iPar)):
            # check if current parent is offspring of other parent
            pastPar = curr_iPar[uni_curr_iPar[iPar] == curr_iChi]
            if len( pastPar) > 0:
                sel_pastPar = pastPar == uni_id_pastPar
            else:
                sel_pastPar = np.array([False])
            if sel_pastPar.sum() > 0:
                # add 1 to previous parent triggering generation
                curr_iGen[iPar]         = uni_iGen_pastPar[sel_pastPar][0]+1
                uni_iGen_pastPar[iPar]  = uni_iGen_pastPar[sel_pastPar][0]+1
                uni_id_pastPar[iPar]    = uni_curr_iPar[iPar]
            else:
                curr_iGen[iPar]         = iGen
                uni_iGen_pastPar[iPar]  = iGen
                uni_id_pastPar[iPar]    = uni_curr_iPar[iPar]
                iGen += 1  # assign new trig generation
        # save evID, trigger generation in dictionary
        dGen[famID]    = np.zeros( (3, len(uni_id_pastPar)))

        dGen[famID][0] = uni_id_pastPar
        dGen[famID][1] = curr_iGen
        # =========================3========================================
        #             compute ave. leaf depth
        # ==================================================================
        sel_endLeaf = ~np.in1d(curr_iChi, curr_iPar)
        dGen[famID][2] = np.ones( len(dClust[famID]))*curr_iGen[1::][sel_endLeaf].mean()
    return dGen

def offspring_gen_test( dClust, dNND, f_eta_0, **kwargs):
    #=========================1========================================
    #               add origin times from seisCat to dNND
    #==================================================================
    # sort  dNND and seisCat by offspring ID!!- seisCat.data['Time'] is added to dNND
    # sortSel = np.argsort( dNND['aEqID_c'])
    # for tag in list(dNND.keys()):
    #     dNND[tag] = dNND[tag][sortSel]
    # seisCat.sortCatalog('Time')
    # firstEvID = seisCat.data['N'][0]
    # seisCat.sortCatalog('N')
    # # add offspring origin time to dNND
    # sel = firstEvID == seisCat.data['N']
    # dNND['at_c'] = seisCat.data['Time'][~sel]
    # check that dNND is sorted by time
    # if 'Time' in dNND.keys():
    #     i_sort = np.argsort( dNND['Time'])
    #     dNND   = data_utils.selDicAll(dNND, i_sort)
    # else:
    #     error_str = "'Time' key missing, add offspring origin time to dNND"
    #     raise ValueError( error_str)
    #=========================2========================================
    #            count generations of offspring events
    #==================================================================
    dGen = {}
    l_famID = list( dClust.keys())
    # singles are all 0 generation
    dGen['0'] = np.zeros( (3, len( dClust['0'])))
    # set ev IDs in new dic
    dGen['0'][0] = dClust['0']
    # ave LD = 1
    dGen['0'][2] = np.ones( len( dClust['0']))
    # ignore singles below
    l_famID.remove( '0')
    for famID in l_famID:
        ###find ori. time for each child
        sel_chi_t  = np.in1d(  dNND['aEqID_c'], dClust[famID])
        # filter for parent - child NND < eta_0
        sel_chi_t2 = dNND['aNND'][sel_chi_t] < f_eta_0
        curr_iPar  = dNND['aEqID_p'][sel_chi_t][sel_chi_t2]
        curr_iChi  = dNND['aEqID_c'][sel_chi_t][sel_chi_t2]
        curr_tChi  = dNND['Time'][np.in1d( dNND['aEqID_c'],curr_iChi)]
        # curr_tChi = np.zeros( len( curr_iPar))
        # for iP in range( len( curr_iChi)):
        #     curr_tChi[iP] = dNND['at_c'][dNND['aEqID_c']==curr_iChi[iP]]

        ##sort cluster IDs with respect to offspring time
        sel_sort  = np.argsort( curr_tChi)

        curr_tChi = curr_tChi[sel_sort]
        curr_iChi = curr_iChi[sel_sort]
        curr_iPar = curr_iPar[sel_sort]

        # parent IDs have one less element than complete cluster (first event has no parent)
        #sel_sort =  np.hstack((0, sel_sort+1))
        # make sure dClust[famID] = dNND['aEqID_c']+
        firstID = dClust[famID][0]

        uni_curr_iPar = np.unique( curr_iPar)
        # sort unique parents by time
        print( uni_curr_iPar)
        uni_par_times = curr_tChi[np.in1d(curr_iChi, uni_curr_iPar)]
        uni_curr_iPar = curr_iChi[np.in1d(curr_iChi, uni_curr_iPar)]
        sort_uni_par  = np.argsort( uni_par_times)
        uni_curr_iPar = uni_curr_iPar[sort_uni_par]
        # check if parent of first pair needs to be added
        if np.isin( curr_iPar[0], uni_curr_iPar).sum() == 0:
            uni_curr_iPar = np.hstack(( curr_iPar[0], uni_curr_iPar))
        # add end leafs (offspring that are not parents)
        sel_endLeaf   = ~np.in1d(curr_iChi, curr_iPar)
        uni_curr_iPar = np.hstack((uni_curr_iPar, curr_iChi[sel_endLeaf]))
        #----------initiate new vectors------------------------
        uni_iGen_pastPar = np.zeros( len(curr_tChi)+1)
        uni_id_pastPar   = np.zeros( len(curr_tChi)+1)
        ## assign chronological triggering generation
        curr_iGen        = np.zeros( len(curr_tChi)+1)
        iGen          = 0
        for iPar in range( len(uni_curr_iPar)):
            # # assign trig gen starting from oldest parent
            # sel_hier_par = curr_iPar == uni_curr_iPar[iPar]
            # # print("current parent: ", uni_curr_iPar[iPar], "offspring: ", curr_iChi[sel_hier_par])
            # # print( "past parents", uni_id_pastPar)
            # # print( "past parent iGen", uni_iGen_pastPar)
            # check if current parent is offspring of other parent
            pastPar = curr_iPar[uni_curr_iPar[iPar] == curr_iChi]
            #print( uni_curr_iPar[iPar], pastPar, uni_id_pastPar)
            if len( pastPar) > 0:
                sel_pastPar = pastPar == uni_id_pastPar
            else:
                sel_pastPar = np.array([False])
            if sel_pastPar.sum() > 0:
                print(uni_curr_iPar[iPar], "past parent: ", pastPar, "trig gen: ", uni_iGen_pastPar[sel_pastPar][0]+1)
                # add 1 to previous parent triggering generation
                curr_iGen[iPar]         = uni_iGen_pastPar[sel_pastPar][0]+1
                uni_iGen_pastPar[iPar]  = uni_iGen_pastPar[sel_pastPar][0]+1
                uni_id_pastPar[iPar]    = uni_curr_iPar[iPar]
            else:
                print( "trig gen: ", iGen)
                curr_iGen[iPar]         = iGen
                uni_iGen_pastPar[iPar]  = iGen
                uni_id_pastPar[iPar]    = uni_curr_iPar[iPar]
                iGen += 1  # assign new trig generation

        print( uni_id_pastPar)
        print( curr_iGen)
        # save evID, trigger generation in dictionary
        dGen[famID]    = np.zeros( (3, len(uni_id_pastPar)))
        dGen[famID][0] = uni_id_pastPar
        dGen[famID][1] = curr_iGen
        # =========================3========================================
        #             compute ave. lead depth
        # ==================================================================
        # end leafs = events without offspring, curr_iPar is already < eta_0)
        #print( len( curr_iGen), len( curr_iChi), len( uni_curr_iPar))
        sel_endLeaf = ~np.in1d(curr_iChi, curr_iPar)
        print( 'end leaf off. ID', curr_iChi[sel_endLeaf])
        print( ' leaf depths ', curr_iGen[1::][sel_endLeaf])
        print( 'mean leaf depth: ', curr_iGen[1::][sel_endLeaf].mean())
        dGen[famID][2] = np.ones( len(dClust[famID]))*curr_iGen[1::][sel_endLeaf].mean()
    # recall data structure:
    # dGen[famID][0] = evIDs
    # dGen[famID][1] = trig generation
    # dGen[famID][2] = ave. leaf depth
    return dGen
#=================================================================================
#                      create random catalogs
#=================================================================================
# create uniform times
def rand_rate_uni( N, tmin, tmax, **kwargs):
    """  draw N random numbers out of a Poisson distribution defined by mu, between tmin and tmax,

    kwargs: - random uniform variable between min and max

    return: vector of N origin times between tmin and tmax """
    return np.random.uniform( tmin, tmax, size = N)


# ------------------------------------------------------------------------------------------
def haversine(lon1, lat1, lon2, lat2, **kwargs):
    """
    haversine formula implementation
    https://en.wikipedia.org/wiki/Great-circle_distance
    great circle distance between two points
    :input   lon1, lat1
             lon2, lat2

    		  gR - Earth radius (global variable)
    :output  distance - great circle distance in kilometer
    """
    i_radius = 6371
    # convert to radians
    lon1 = lon1 * np.pi / 180
    lon2 = lon2 * np.pi / 180
    lat1 = lat1 * np.pi / 180
    lat2 = lat2 * np.pi / 180
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = i_radius * c
    return distance

# ==================================4==============================================================
#                       T-R density plots
# =================================================================================================
def plot_R_T( a_T, a_R, f_eta_0, **kwargs):
    """
        - plot rescaled distance over rescaled time
        Parameters:
    dPar = {'binx': .1, 'biny': .1,  # used for density and gaussian smoothing
            'sigma': None,  # if None: default = n**(-1./(d+4)),
            'Tmin': -8, 'Tmax': 0,
            'Rmin': -5, 'Rmax': 3,
            'cmap': plt.cm.RdYlGn_r}
    Use kwargs['dPar'] = python dictionary
                        'binx', 'biny', etc. to overwrite
                        defaults for specific or all parameters
    :param kwargs:
    :return: fig - figure handle - use fig.axes to get list of corresponding axes
    """
    dPar = {'binx': .1, 'biny': .1,  # used for density and gaussian smoothing
            'sigma': None,  # if None: default = n**(-1./(d+4)),
            'Tmin': -8, 'Tmax': 0,
            'Rmin': -5, 'Rmax': 3,
            'cmap': plt.cm.RdYlGn_r}
    if 'dPar' in kwargs.keys() and kwargs['dPar'] is not None:
        for tag in kwargs['dPar'].keys():
            print( f"overwrite plot_R_T param: {tag}={kwargs['dPar'][tag]}")
            dPar[tag] = kwargs['dPar'][tag]
    a_Tbin = np.arange(dPar['Tmin'], dPar['Tmax'] + 2 * dPar['binx'], dPar['binx'])
    a_Rbin = np.arange(dPar['Rmin'], dPar['Rmax'] + 2 * dPar['biny'], dPar['biny'])
    sel = a_T > 0
    XX, YY, ZZ = data_utils.density_2D(np.log10(a_T[sel]), np.log10(a_R[sel]), a_Tbin, a_Rbin, sigma=dPar['sigma'])

    fig = plt.figure( figsize=(7, 9))
    ax = plt.subplot(111)
    ax.set_title('Nearest Neighbor Pairs in R-T')
    # ------------------------------------------------------------------------------
    normZZ = ZZ * (dPar['binx'] * dPar['biny'] * len(a_R))
    plot1 = ax.pcolormesh(XX, YY, normZZ, cmap=dPar['cmap'])
    cbar = plt.colorbar(plot1, orientation='horizontal', shrink=.5, aspect=20, )
    # ax.plot(  np.log10( a_T), np.log10( a_R), 'wo', ms = 1.5, alpha = .2)
    # plot eta_0 to divide clustered and background mode
    ax.plot([dPar['Tmin'], dPar['Tmax']], -np.array([dPar['Tmin'], dPar['Tmax']]) + f_eta_0, '-', lw=1.5, color='w')
    ax.plot([dPar['Tmin'], dPar['Tmax']], -np.array([dPar['Tmin'], dPar['Tmax']]) + f_eta_0, '--', lw=1.5, color='.5')
    # -----------------------labels and legends-------------------------------------------------------
    # cbar.set_label( 'Event Pair Density [#ev./dRdT]')
    cbar.set_label('Number of Event Pairs', labelpad=-60)
    ax.set_xlabel('Rescaled Time')
    ax.set_ylabel('Rescaled Distance')
    ax.set_xlim(dPar['Tmin'], dPar['Tmax'])
    ax.set_ylim(dPar['Rmin'], dPar['Rmax'])
    # fig.axes
    return fig