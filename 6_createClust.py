'''
        Created on August 16, 2016

    This is a key step in the clustering analysis during which the full chains
    of triggered events are assembled based an NND and the criteria (eta_0) used
    to differentiate background events (singles) and event families.
    
    For each event pair with NND <= eta_0 a new families will be created or 

- create dictionary with all families for which '0' are singles
  and all other cluster are integer-strings followed by the associated eqID number
  e.g.
  {  '0' : np.array([[1243245,4253455343]]),
      '1': np.array([[5235,43455343,3456,56652,54]]),
      '2':}
      
- family identification needs eta_0 from 2_compute_eta_0.py and eta for each NND pair from 3_compute_eta.py

-Note: clusterAnalysis.assempleCluster_NND takes R0 = -1 and onlyR0 = -1 for OK
                     R0      - use both eta_0 and R0 to separate families from singles
                     onlyR0  - uses only distance criterion to separate families and singles

@author: tgoebel Thomas Goebel University of Memphis
'''
import matplotlib as mpl
#mpl.use( 'Agg')

import numpy as np
import os
import pylab as plb

import scipy.io 
import my_stats
import fractalAnalysis as fractal
import dictionaryIO as dicIO
import clusterAnalysis as cluster

from SmoothingUtils import *
from SeisCatDic import *
from InstantPlot import *
from GloParamsOU import *

smooth_util = SmoothingUtils()
catOri      = SeismicCatalog( type ='GPS')

#=================================1==============================================
#          dir, file, params
#================================================================================
dPar  = {   'catName'     : 'reloc', # 'ANSS', # reloc
            #-------event selection-----------------
            'minMag'    : None, 'maxMag' : 9,
            'depth_max' : 30,
            'tmin'      : 2009, 'tmax': 2018.5,

            'R_0'       : None, # -1, #log10R0 or None 
            'onlyR0'    : False, #None, #-1,
            #'eta_0'       : -4.6, #=5,
            'vMc'         :  np.array( [1.75, 2.0, 2.5, 3.0, 3.5] ),

          }



gPar = GloParamsOU( catName = dPar['catName'])

#================================================================================
#                            load data
#================================================================================
if dPar['catName'] == 'reloc':
    os.chdir( gPar.dir['OGS'])
    catOri.readCatalogFromMatBin( gPar.file['reloc'])
else:
    os.chdir( gPar.dir['%s'%(dPar['catName'])])
    catOri.readCatDic( gPar.file['cat_%s'%(dPar['catName'])])
print dPar['catName'], 'catalog size', catOri.size()
catOri.selectEvents( None,           dPar['depth_max'], 'Depth')
catOri.selectEvents( dPar['tmin'],   dPar['tmax'],      'Time')
catOri.selectEvents( dPar['minMag'], dPar['maxMag'],    'MAG')
print 'catalog size', catOri.size()


iMc = 0
for dPar['Mc'] in dPar['vMc']:
    gPar = GloParamsOU( catName = dPar['catName'], Mc = dPar['Mc'])
    print '------------cat. name and Mc', dPar['Mc'], dPar['catName'], '---------------'
 
    #=================================0==============================================
    #                            load data
    #================================================================================
    os.chdir( gPar.dir['data'])
    mData = np.loadtxt( gPar.file['b_D_Mc']).T
    dConst = { 'b' :  round( mData[0], 1), 
              'Mc' :  dPar['Mc'],
               'D' :  round( mData[2], 1),}

    # preprocessing
    projCat = SeismicCatalog( type ='GPS')
    projCat = MultiCatalog( type = 'GPS').copy( catOri)
    projCat.selectEvents( dConst['Mc'], None, 'MAG')
    print 'catalog size after pre-processing', projCat.size()
    
    os.chdir( gPar.dir['data'])
    if os.path.isfile( gPar.file['eta_0']):
        dEta_0 = scipy.io.loadmat( gPar.file['eta_0'])
        dPar['eta_0'] = dEta_0['eta_0'].mean()
        print 'eta 0', dPar['eta_0']
    else:
        error_str = 'run 2_compute_eta_0, eta 0 file not found'
        raise ValueError, error_str
 
    #==================================2=============================================
    #                        load NNDs for all event pairs
    #================================================================================  
    os.chdir( gPar.dir['data'])
    if os.path.isfile( gPar.file['NND']):
        import scipy.io
        dNND = scipy.io.loadmat( gPar.file['NND'])
        dicIO.filterDicColumns( dNND)        
    else:
        print 'start recomputing NND'
        error_str = 'NND file not found, run 3_compute_eta.py'
        raise ValueError, error_str
        #dNND  = cluster.NND_eta(  projCat, dConst)    
    dNND['aNND'] = np.log10( dNND['aNND'])
 
    #==================================3=============================================
    #                           assemble similarity clusters
    #================================================================================
    print 'similarity threshold', dPar['eta_0']
    # clustering according to eta_0 similarity criteria
    if 'onlyR0' in dPar.keys() and dPar['onlyR0'] is not None and dPar['onlyR0'] != False:
        dClust = cluster.assembleClusters_NND( dNND['aEqID_c'], dNND['aEqID_p'], dNND['aNND'], dPar['eta_0'], 
                                           onlyR0 = dPar['onlyR0'], aR = dNND['aR'],  useLargerEvents = False)
         
    else: 
        dClust = cluster.assembleClusters_NND( dNND['aEqID_c'], dNND['aEqID_p'], dNND['aNND'], dPar['eta_0'], 
                                           R0 = dPar['R_0'], aR = dNND['aR'],  useLargerEvents = False)
 
    #=================================4==========================================================================
    #                           save results
    #============================================================================================================
    os.chdir( gPar.dir['data'])
    print os.getcwd()
    ## save main shock catalog
    if 'R_0' in dPar.keys() and dPar['R_0'] is not None:
        print gPar.file['clID'].replace( '.dic', 'R0_%s.dic'%(dPar['R_0']))
        dicIO.writeBin( gPar.file['clID'].replace( '.dic', 'R0_%s.dic'%(dPar['R_0'])), dClust)
    elif 'onlyR0' in dPar.keys() and dPar['onlyR0'] is not None:
        print gPar.file['clID'].replace( '.dic', 'onlyR0_%s.dic'%(dPar['onlyR0']))
        dicIO.writeBin( gPar.file['clID'].replace( '.dic', 'onlyR0_%s.dic'%(dPar['onlyR0'])), dClust)#'%s_clIDs_Mc_%.1f.dic'%( dPar['catName'], dPar['Mc']), dClust)
    else:
        print gPar.file['clID']
        dicIO.writeBin( gPar.file['clID'], dClust)#'%s_clIDs_Mc_%.1f.dic'%( dPar['catName'], dPar['Mc']), dClust)
     


    iMc += 1










