# python2.7
'''    
        Created on May 16th, 2019

     key step in the analysis during which complete families of triggered events are assembled
     family, cluster, trigger chain are used inter-changeably
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
              {  '0' : np.array([[1243245,4253455343]]),
                  '1': np.array([[5235,43455343,3456,56652,54]]),
                  '2':  ....}
      

@author: Thomas Goebel University of Memphis
'''
import matplotlib as mpl
#mpl.use( 'Agg')
import matplotlib.pyplot as plt

import numpy as np
import os

#------------------------------my modules-------------------------------------- 
import src.data_utils as dataIO
import src.clustering as clustering
from src.EqCat import *

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

dPar  = {   'a_Mc'        :  np.array([3.0, 4.0]), #np.array( [2.0, 2.5, 3.0, 3.5]),
            #separate clustered and background
            'eta_0'       : -5.0,
            'testPlot'    : True,
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
    # cut below current completeness
    eqCatMc.copy( eqCat)
    eqCatMc.selectEvents( f_Mc, None, 'Mag')
    print( 'current catalog size: ',eqCatMc.size())
    # load nearest neighbor distances
    NND_file = '%s_NND_Mc_%.1f.mat'%(os.path.basename( file_in).split('.')[0], f_Mc)
    dNND = dataIO.loadmat( os.path.join( data_dir, NND_file))
    dNND['aNND'] = np.log10( dNND['aNND'])
 
    #==================================3=============================================
    #                      assemble clusters
    #================================================================================
    print( 'similarity threshold', dPar['eta_0'])
    # clustering according to eta_0 similarity criteria
    dClust = clustering.assembleClusters_NND( dNND['aEqID_c'], dNND['aEqID_p'], dNND['aNND'], dPar['eta_0'], useLargerEvents = False)
    for key in sorted( dClust.keys()):
        print key, len( dClust[key])
    print( sgr)
    #=================================4==========================================================================
    #                           save results
    #============================================================================================================
    os.chdir( gPar.dir['data'])
    print gPar.file['clID']
    dicIO.writeBin( gPar.file['clID'], dClust)#'%s_clIDs_Mc_%.1f.dic'%( dPar['catName'], dPar['Mc']), dClust)
    iMc += 1

    #=================================5==========================================================================
    #                           test plots
    #============================================================================================================
    #TODO










