#!python2.3
'''
Created on April 10th,  2019

- compute nearest-neighbor distance (NND) between all event pairs (see equ. 1 in  Zaliapin & Ben-Zion, 2013)
- test -plot histogram of NNDs: Figure 4c in Zaliapin & Ben-Zion, 2013

output: 'data/%s_NND_Mc_%.1f.mat'%(dPar['catName'], dPar['Mc'])
        which is a python dictionary with: 
            {     'aNND'       : aNND,     - nearest neighbor space-time magnitude distance
                 'aEqID_p'    : np.array  - ID of the parent event
                 'aEqID_c'    : np.array  - ID of the child  event
            }

TODO:
    - constrain Mc, b and D independently through statistical analysis of the actual data

@author: tgoebel
'''
import time
#------------------------------------------------------------------------------
import matplotlib as mpl
#mpl.use( 'Agg') # turn off interactive plot
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import os
 
#------------------------------my modules-------------------------------------- 
import src.clustering as clustering
from src.EqCat import EqCat

eqCat = EqCat( )

#=================================1==============================================
#                            dir, file, params
#================================================================================
dir_in = 'data'
file_in= 'hs_1981_2011_all.mat'

#file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar  = {   'aMc'         :  np.array([3.0]), #3.0, 4.0]), #np.array( [2.0, 2.5, 3.0, 3.5]),
            # fractal dimension and b for eq. (1)
            'D'           : 1.6, # TODO: - these values should be contrained independently
            'b'           : 1.0, # use: https://github.com/tgoebel/magnitude-distribution for b-value
            #=================plotting==============
            'eta_binsize' :  .3,
            'xmin' : -13, 'xmax' : 0,
          }

#=================================2==============================================
#                            load data, select events
#================================================================================
eqCat.loadMatBin(  os.path.join( dir_in, file_in))
print( 'total no. of events', eqCat.size())
eqCat.selectEvents( dPar['aMc'][0], None, 'Mag')
#eqCat.selectEvents( tmin, tmax, 'Time')
print( 'no. of events after initial selection', eqCat.size())
#=================================1==============================================
#                           to cartesian coordinates
#================================================================================
# two ways to do the distance comp: 1 project into equal distance azimuthal , comp Cartersian distance in 3D
#                                   2 get surface distance from lon, lat (haversine), use pythagoras to include depth
eqCat.toCart_coordinates( projection = 'eqdc')

for f_Mc in dPar['aMc']:
    print( '-------------- current Mc:', f_Mc, '---------------------')
    # select magnitude range
    eqCat.selectEvents( f_Mc, None, 'Mag')
    print( 'catalog size after MAG selection', eqCat.size())
    # this dictionary is used in module: clustering
    dConst = {'Mc' : f_Mc,
               'b' : dPar['b'],
               'D' : dPar['D']}
    #==================================2=============================================
    #                       compute space-time-magnitude distance, histogram
    #================================================================================
    eqCat.data['Z'] = eqCat.data['Depth']
    print('depth range: ', eqCat.data['Z'].min(), eqCat.data['Z'].max())
    dCluster = clustering.NND_eta( eqCat, dConst,    distance_3D = False,
                                                    correct_co_located = True, verbose= True)
    ###histogram
    aBins       = np.arange( -13, 1, dPar['eta_binsize'], dtype = float)
    aHist, aBins = np.histogram( np.log10( dCluster['aNND'][dCluster['aNND']>0]), aBins)
    aBins = aBins[0:-1] + dPar['eta_binsize']*.5
    # correct for binsize
    aHist = aHist/dPar['eta_binsize']
    # to pdf (prob. density)
    aHist /= eqCat.size()
    #=================================3==============================================
    #                            save results
    #================================================================================
    import scipy.io
    NND_file = 'data/%s_NND_Mc_%.1f.mat'%( file_in.split('.')[0], f_Mc)
    print( 'save file', NND_file)
    scipy.io.savemat( NND_file, dCluster, do_compression  = True)
    
    #=================================4==============================================
    #                          plot histogram
    #================================================================================
    # load eta_0 value - only for plotting purposes
    eta_0_file = '%s/%s_Mc_%.1f_eta_0.txt'%(dir_in, file_in, f_Mc)
    if os.path.isfile( eta_0_file):
        print( 'load eta_0 from file'),
        f_eta_0 = np.loadtxt( eta_0_file, dtype = float)
        print( 'eta_0',f_eta_0)
    else:
        f_eta_0 = -5
        print( 'could not find eta_0 file', eta_0_file, 'use value: ', f_eta_0)

    fig, ax = plt.subplots()
    #ax.plot( vBin, vHist, 'ko')
    ax.bar( aBins, aHist, width =.8*dPar['eta_binsize'], align = 'edge', color = '.5', label = 'Mc = %.1f'%( f_Mc))
    ax.plot( [f_eta_0, f_eta_0], ax.get_ylim(), 'w-',  lw = 2, label = '$N_\mathrm{tot}$=%i'%( eqCat.size()))
    ax.plot( [f_eta_0, f_eta_0], ax.get_ylim(), 'r--', lw = 2, label = '$N_\mathrm{cl}$=%i'%( dCluster['aNND'][dCluster['aNND']<1e-5].shape[0]))

    ax.legend( loc = 'upper left')
    ax.set_xlabel( 'NND, log$_{10} \eta$')
    ax.set_ylabel( 'Number of Events')
    ax.grid( 'on')
    ax.set_xlim( dPar['xmin'], dPar['xmax'])
    plt.show()

    plotFile = 'plots/%s_NND_hist_Mc_%.1f.png'%( file_in.split('.')[0], f_Mc)
    print( 'save plot', plotFile)
    #plt.savefig( plotFile)
    plt.clf()
    








