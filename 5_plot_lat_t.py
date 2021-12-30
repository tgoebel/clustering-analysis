'''
Created on May 16, 2019

- plot cluster families with eta <= eta_0 
- plot lat and time (dec. year)

@author: tgoebel - Thomas Goebel University of Memphis
'''
import matplotlib as mpl
#mpl.use( 'Agg') # uncomment for interactive plotting

import os
import numpy as np
import matplotlib.pyplot as plt

#------------------------------my modules-------------------------------------- 
import src.data_utils as dataIO
#import src.clustering as clustering
from src.EqCat import EqCat

eqCat   = EqCat() # original catalog
eqCatMc = EqCat() # this catalog will be modified with each Mc iteration
catChild=  EqCat()
catParent= EqCat()

#=================================1==============================================
#                            dir, file, params
#================================================================================
data_dir = 'data'
plot_dir = 'plots'
file_in  = 'hs_1981_2011_all.mat'

dPar  = {   'a_Mc'        :  np.array([4.0]), #np.array( [2.0, 2.5, 3.0, 3.5]),
            #separate clustered and background
            'eta_0'       : -5.0, # run 2_eta_0.py and
                                  # if file exists: default = load this value from ASCII file
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
    eta_0_file = '%s/%s_Mc_%.1f_eta_0.txt'%(data_dir, file_in, f_Mc)
    # load eta_0 value
    if os.path.isfile( eta_0_file):
        print( 'load eta_0 from file'),
        f_eta_0 = np.loadtxt( eta_0_file, dtype = float)
        print( f_eta_0)
    else:
        print( 'could not find eta_0 file', eta_0_file, 'use value from dPar', dPar['eta_0'])
        f_eta_0 = dPar['eta_0']
    # cut below current completeness
    eqCatMc.copy( eqCat)
    eqCatMc.selectEvents( f_Mc, None, 'Mag')
    print( 'current catalog size: ',eqCatMc.size())
    # load nearest neighbor distances
    NND_file = '%s_NND_Mc_%.1f.mat'%(os.path.basename( file_in).split('.')[0], f_Mc)
    dNND = dataIO.loadmat( os.path.join( data_dir, NND_file))
    print( dNND.keys())
    dNND['aNND'] = np.log10( dNND['aNND'])
    #==================================3=============================================
    #                          "declustering" step
    #================================================================================  
    #catChild, catPar = create_parent_child_cat( projCat, dNND)
    catChild.copy( eqCat)
    catParent.copy( eqCat)
    catChild.selEventsFromID( dNND['aEqID_c'], repeats = True)
    catParent.selEventsFromID( dNND['aEqID_p'], repeats = True)
    print( 'tot. ev', eqCatMc.size(), 'parents', np.unique( catParent.data['N']).shape[0], 'children', np.unique( catChild.data['N']).shape[0])
    #==================================4=============================================
    #                          spanning tree
    #================================================================================
    plt.figure( 1)
    ax = plt.subplot(111)  
    for iEv in range( catParent.size()):
        print( 'MS', int( catParent.data['N'][iEv]), catParent.data['Time'][iEv], eqCatMc.data['Time'][iEv])

        if dNND['aNND'][iEv] < dPar['eta_0']:#triggered cluster
            ax.plot( [catParent.data['Time'][iEv]], [catParent.data['Lat'][iEv]], 'ro', ms = 12, alpha = .2)
            ax.plot( [catParent.data['Time'][iEv],catChild.data['Time'][iEv]],
                      [catParent.data['Lat'][iEv], catChild.data['Lat'][iEv]], 'k-', marker = 'o', ms = 4, mew =1, mfc = 'none')
        else: # independent events
            ax.plot( [catChild.data['Time'][iEv]], [catChild.data['Lat'][iEv]], 'bo', ms = 5, alpha = .6)
    
    #ax.set_xlim( 2009, 2017)
    #=================================3==============================================
    #                           save results
    #================================================================================

    plt.figure(1)
    #plt.savefig( '%s/%s_spanningTree_Mc_%.1f.png'%(plot_dir, file_in.split('.')[0], f_Mc))
    ## save main shock catalog
    plt.show()
    plt.clf()


    iMc += 1
