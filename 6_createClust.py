# python3.7
'''    
        Created on Oct 7th, 2019

     key step in the analysis during which complete families of triggered events are assembled

     1) select all event pairs with NND <= eta_0
     2) place each event into a new cluster (family) with unique cluster ID
        or append to existing cluster if parent or offspring event are found in list
        of previously created clusters
    
    
    Input:   NND_file = str()
             eta_0    = float()
             
      
    Output:
            dictionary with all event families 
            dic['0'] = singles
            all other cluster are integer-strings followed by the associated eqID number
              e.g.
              {   '0': np.array([[1243245,4253455343]]),
                  '1': np.array([[5235,43455343,3456,56652,54]]),
                  '2':  ....}
            Note that:
             1)   output is saved as matlab binary but file cannot
                  be read with matlab because variable names are integers
             2) cluster '0' are singles
             3) mainshocks can be defined as largest  event or
                first event in a family

@author: Thomas Goebel - University of Memphis
'''
import matplotlib as mpl
#mpl.use( 'Agg')
import matplotlib.pyplot as plt

import numpy as np
import os, scipy.io

#------------------------------my modules-------------------------------------- 
import src.data_utils as dataIO
import src.clustering as clustering
from src.EqCat import *

eqCat   = EqCat() # original catalog
eqCatMc = EqCat() # this catalog will be modified with each Mc iteration

#=================================1==============================================
#                            dir, file, params
#================================================================================
data_dir   = 'data'
plot_dir   = 'plots'
file_in    = 'hs_1981_2011_all.mat'


dPar  = {   'a_Mc'        :  np.array([3.0]), #3.0, 4.0]), #np.array( [2.0, 2.5, 3.0, 3.5]),
            #separate clustered and background
            # set to None or False to use value from file,requires results from: 2_eta_0.py
            'eta_0'       : None, #-5.0,
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
    clust_file = file_in.replace( 'all.mat', 'Mc_%.1f_clusters.mat'%( f_Mc))
    eta_0_file = '%s/%s_Mc_%.1f_eta_0.txt'%(data_dir, file_in, f_Mc)
    if os.path.isfile( eta_0_file):
        print( 'load eta_0 from file'),
        f_eta_0 = np.loadtxt( eta_0_file, dtype = float)
        print( 'eta_0',f_eta_0)
    else:
        f_eta_0 = -5
        print( 'could not find eta_0 file', eta_0_file, 'use value: ', f_eta_0)


    # cut below current completeness
    eqCatMc.copy( eqCat)
    eqCatMc.selectEvents( f_Mc, None, 'Mag')
    print( 'current catalog size: ',eqCatMc.size())
    # load nearest neighbor distances
    NND_file = '%s_NND_Mc_%.1f.mat' % (file_in.split('.')[0], f_Mc)
    dNND = dataIO.loadmat( os.path.join( data_dir, NND_file))
    dNND['aNND'] = np.log10( dNND['aNND'])
 
    #==================================3=============================================
    #                      assemble clusters
    #================================================================================
    print( 'similarity threshold', dPar['eta_0'])
    # clustering according to eta_0 similarity criteria
    dClust = clustering.compileClust( dNND, f_eta_0, useLargerEvents = False)
    #=================================4==========================================================================
    #                           save results
    #============================================================================================================
    scipy.io.savemat( os.path.join( data_dir,clust_file), dClust, do_compression=True)











