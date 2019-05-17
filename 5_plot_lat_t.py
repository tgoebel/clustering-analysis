'''
Created on May 16, 2019

- plot cluster families with eta <= eta_0 
- plot lat and time (dec. year)

@author: tgoebel
'''
import matplotlib as mpl
mpl.use( 'Agg') # uncomment for interactive plotting

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
#------------------------------my modules-------------------------------------- 
import data_utils
import src.clustering as clustering
from src.EqCat import *

eqCat   = EqCat() # original catalog
eqCatMc = EqCat() # this catalog wil be modfied with each Mc iteration
catChild=  EqCat()
catParent= EqCat()

#=================================1==============================================
#                            dir, file, params
#================================================================================
dir_in = '%s/data/quakeData/SCSN/relocated'%( os.path.expanduser( '~'))
file_in= 'hs_1981_2011_all.mat'

#file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar  = {   'a_Mc'        :  np.array([3.0, 4.0]), #np.array( [2.0, 2.5, 3.0, 3.5]),

            #separate clustered and background
            'eta_0'       : -5.0,}


#=================================2==============================================
#                            load data, select events
#================================================================================
eqCat.loadMatBin(  os.path.join( dir_in, file_in))
print( 'total no. of events', eqCat.size())
eqCat.selectEvents( dPar['a_Mc'][0], None, 'Mag')
#eqCat.selectEvents( tmin, tmax, 'Time')
print( 'no. of events after initial selection', eqCat.size())

iMc = 0
for dPar['Mc'] in dPar['aMc']:
    # cut below current completeness
    eqCatMc.copy( eqCat)
    eqCatMc.selectEvents( f_Mc, None, 'Mag')
    print( 'current catalog size: ',eqCatMc.size())
        
    #=================================1==============================================
    #                           to cartesian coordinates
    #================================================================================
    # two ways to do the distance comp: 1 project into equal distance azimuthal , comp Cartersian distance in 3D
    #                                   2 get surface distance from lon, lat (haversine), use pythagoras to include depth
    projCat.toCart_coordinates( projection = 'aeqd')
   
    #==================================2=============================================
    #                       compute rescale time and distance
    #================================================================================  

    NND_file = '%s_NND_Mc_%.1f.mat'%(dPar['catName'], dPar['Mc'])
    os.chdir( gPar.dir['data'])
    if os.path.isfile( NND_file):
        import scipy.io
        dNND = scipy.io.loadmat( NND_file)
        dicIO.filterDicColumns( dNND)        
    else:
        dNND  = cluster.NND_eta(  projCat, dConst)    
    print dNND.keys()
    dNND['aNND'] = np.log10( dNND['aNND'])
    #==================================3=============================================
    #                         declustering
    #================================================================================  
    #catChild, catPar = create_parent_child_cat( projCat, dNND)
    catChild = MultiCatalog(type='GPS').selEventsFromID( dNND['aEqID_c'], projCat, repeats = True)
    catPar   = MultiCatalog(type='GPS').selEventsFromID( dNND['aEqID_p'], projCat, repeats = True)

    plb.figure( 1)
    ax = plb.subplot(111)
    #==================================4=============================================
    #                          spanning tree
    #================================================================================  
    for iEv in xrange( catPar.size()):
        print 'tot. ev', projCat.size(), 'parents', np.unique( catPar.data['N']).shape[0], 'children', np.unique( catChild.data['N']).shape[0]
        print 'MS', int( catPar.data['N'][iEv]), catPar.data['Time'][iEv], projCat.data['Time'][iEv]

        if dNND['aNND'][iEv] < dPar['eta_0']:#triggered cluster
            ax.plot( [catPar.data['Time'][iEv]], [catPar.data['Lat'][iEv]], 'ro', ms = 12, alpha = .2)
            ax.plot( [catPar.data['Time'][iEv],catChild.data['Time'][iEv]],
                      [catPar.data['Lat'][iEv], catChild.data['Lat'][iEv]], 'k-', marker = 'o', ms = 4, mew =1, mfc = 'none')
        else: # independent events
            ax.plot( [catChild.data['Time'][iEv]], [catChild.data['Lat'][iEv]], 'bo', ms = 5, alpha = .6)
    
    #ax.set_xlim( 2009, 2017)
    #=================================3==============================================
    #                           save results
    #================================================================================
    os.chdir( gPar.dir['plots'])
    plb.figure(1)
    plb.savefig( '%s_spanningTree_Mc_%.1f.png'%( dPar['catName'], dPar['Mc']))
    ## save main shock catalog
    plb.show()
    plb.clf()


    iMc += 1










