'''
Created on Oct 07, 2019

    1) count number of events after largest mag event (MS) within each family
    2) count singles as events with 0 aftershocks

@author: tgoebel
'''
import matplotlib as mpl
mpl.use( 'Agg')
import numpy as np

import os
import matplotlib.pyplot as plt

#----------------------my modules-------------------------------------------------------- 
import src.data_utils as data_utils
#import src.clustering as clustering
from src.EqCat import *

eqCat   = EqCat() # original catalog
eqCatMc = EqCat() # this catalog will be modified with each Mc iteration


#=================================1==============================================
#                            dir, file, params
#================================================================================
data_dir   = 'data'
plot_dir   = 'plots'
file_in    = 'hs_1981_2011_all.mat'
clust_file = file_in.replace( 'all.mat', 'clusters.mat')

#=================================1==============================================
#          dir, file, params
#================================================================================

dPar  = {
            'a_Mc'         :  np.array([  3.0, 4.0]), # , 3.0, 4.0]), #3.0,4.0]),

            'alpha'       :  1, #exponent for test plot
            #=================plotting==============
            'plotFormat' : 'png',
            }

#=================================2==============================================
#                            load data, select events
#================================================================================
eqCat.loadMatBin(  os.path.join( data_dir, file_in))
print( 'total no. of events', eqCat.size())
eqCat.selectEvents( dPar['a_Mc'][0], None, 'Mag')
#eqCat.selectEvents( tmin, tmax, 'Time')
print( 'no. of events after initial selection', eqCat.size())


iMc = 0
for f_Mc in dPar['a_Mc']:
    # load file with IDs of events within family
    clust_file = file_in.replace( 'all.mat', 'Mc_%.1f_clusters.mat'%( f_Mc))
    dClust = data_utils.loadmat( os.path.join( data_dir,clust_file), )

    # cut below current completeness
    eqCatMc.copy( eqCat)
    eqCatMc.selectEvents( f_Mc, None, 'Mag')
    n_aboveMc = eqCatMc.size()
    print( 'current catalog size: ',eqCatMc.size())

    #=================================1==========================================================================
    #                     singles are counted as MS with 0 AS
    #============================================================================================================
    print( 'total number of clusters', len(  dClust.keys()), 'no. of BG events', dClust['0'].shape[0])
    a_ID_single  = dClust['0']

    # IDs of BG events
    a_iSel       = np.zeros( eqCatMc.size(), dtype = int)
    a_mag_single = np.zeros( len( a_ID_single))
    a_N_AS_single= np.zeros( len( a_ID_single))
    a_N_FS_single= np.zeros( len( a_ID_single))
    for i in range( a_ID_single.shape[0]):
        # event ID may be in catalog more than once
        sel_ev          = eqCatMc.data['N'] == a_ID_single[i]
        a_mag_single[i] = eqCatMc.data['Mag'][sel_ev][0]
        a_iSel[sel_ev] = 1#catalog.data['N'][catalog.data['N']==aEqID[i]][0]
        if sel_ev.sum() != 1:
            error_str = 'more than event found', eqCatMc.data['N'][sel_ev]
            raise( ValueError( error_str))
    ### remove singles from catalog
    eqCatMc.selDicAll( np.logical_not(a_iSel))
    print( 'remaining events', eqCatMc.size(), 'BG events', len( a_mag_single))
    dClust.pop('0') # remove singles
    #=================================2==========================================================================
    #                   get MAGs of MS with aftershocks, count aftershocks
    #============================================================================================================
    a_N_FS    = np.zeros( len( dClust.keys()), dtype = int)
    a_N_AS    = np.zeros( len( dClust.keys()), dtype = int)
    a_MS_mag  = np.zeros( len( dClust.keys()))
    a_MS_ID   = np.zeros( len( dClust.keys()), dtype = int)
    iCl = 0
    for sCl in dClust.keys():
        aEqID = dClust[sCl]# np.unique( dClust[sCl].flatten()) unique is not needed anymore, createCluster has been fixed
        print( 'cl: ', iCl+1,'out of: ', len( dClust.keys()), 'no. of ev. in cl.', len( aEqID), len( np.unique( dClust[sCl])))
        # find MS mag and magnitude of entire family
        atmp_MAG = np.zeros( len( aEqID))
        atmp_Time= np.zeros( len( aEqID))
        a_iSel   = np.zeros( eqCatMc.size(), dtype = int)
        # for each family find: event mag. and origin time
        for iM in range( len( aEqID)):
            sel_ev        = eqCatMc.data['N'] == aEqID[iM]
            if sel_ev.sum() != 1:
                error_str = 'more/less than event found', eqCatMc.data['N'][sel_ev], aEqID[iM]
                raise(  ValueError, error_str)
            atmp_MAG[iM]  = eqCatMc.data['Mag'][sel_ev][0]
            atmp_Time[iM] = eqCatMc.data['Time'][sel_ev][0]
            a_iSel[sel_ev] = 1
        # remove events from catalog
        #catalog.selDicAll( np.logical_not(a_iSel))
        #----------------------------mainshock-------------------------------------------------- 
        selMS     = atmp_MAG == atmp_MAG.max()
        f_tMS     = atmp_Time[selMS][0]
        i_ID_MS   = aEqID[selMS]

        #print( 'tMS', tMS, v_currEqID.shape[0], 'MAG', curr_cat.data['MAG'][selMS][0]
        #----------------------------aftershock-------------------------------------------------- 
        selAS     = atmp_Time > f_tMS
        selFS     = atmp_Time < f_tMS
        #print( 'no. of aftershocks', selAS.sum()
        # save number of aftershocks for each MS mag
        a_MS_mag[iCl] = atmp_MAG[selMS][0]#, dPar['magRound'])
        a_N_AS[iCl]   = selAS.sum()
        a_N_FS[iCl]   = selFS.sum()
        a_MS_ID[iCl]  = int( i_ID_MS[0])
        iCl += 1

    #=================================3==========================================================================
    #                  compare MS+single+FS+AS to original number of events in catalog
    #============================================================================================================
    # combine single without AS with mainshocks that do have aftershocks
    a_N_FS    = np.append( a_N_FS, a_N_FS_single)
    a_N_AS    = np.append( a_N_AS, a_N_AS_single)
    a_MS_mag  = np.append( a_MS_mag, a_mag_single)
    a_MS_ID   = np.append( a_MS_ID, a_ID_single)
    print( 'tot ev. in catalog', n_aboveMc,'tot events in families',a_N_FS.sum() + a_N_AS.sum() + a_MS_mag.shape[0])
    #print( 'N BG', a_mag_single.shape[0], 'FS', a_N_FS_single.sum(), 'AS', a_N_AS_single.sum(), 'MS (MS+BG)', a_MS_mag.shape[0]

    #=================================4==========================================================================
    #                    save to ASCII text
    #============================================================================================================
    file_out = '%s/%s_Nas_MS_Mc_%.1f.txt'%(data_dir, file_in.split('.')[0], f_Mc)#, dPar['magRound'])
    np.savetxt( file_out, np.array([a_MS_mag, a_N_AS, a_N_FS, a_MS_ID]).T, fmt='%10.3f%10i%10i%14i',
                header = 'MAG          N-AS          N-FS        MS-ID; note N_AS=0 highlights singles or FS only')
    iMc += 1










