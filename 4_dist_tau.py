'''
Created on August 16, 2016

- compute inter-event times and distance and normalize by magnitude
- create colormap of event-pair density (i.e. NND events) for
    corresponding rescaled event times and distances

@author: tgoebel
'''
import matplotlib as mpl
mpl.use( 'Agg') # uncomment for interactive plotting

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
#------------------------------my modules-------------------------------------- 
import src.data_utils as data_utils
import src.clustering as clustering
from src.EqCat import *

eqCat   = EqCat() # original catalog
eqCatMc = EqCat() # this catalog wil be modfied with each Mc iteration
catChild=  EqCat()
catParent= EqCat()
np.random.seed( 123456)
#=================================1==============================================
#                            dir, file, params
#================================================================================
data_dir = 'data'
file_in  = 'hs_1981_2011_all.mat'

#file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar  = {   'a_Mc'        :  np.array([4.0]), #np.array( [2.0, 2.5, 3.0, 3.5]),
            # fractal dimension and b for eq. (1)
            'D'           : 1.6, # TODO: - these values should be constrained independently
            'b'           : 1.0,
 
            #===================smoothing parameters=============
            'binx' : .1, 'biny' : .1,# used for density and gaussian smoothing
            'sigma'   : None, #if None: default = n**(-1./(d+4)),
            #=================plotting==============
            'eta_0'       : -5.0, # run: 2_eta_0.py and
                                  # if eta-0-file exists: default = load this value from ASCII file
            #'xmin' : -13, 'xmax' : 0,
            'Tmin' :  -8, 'Tmax' : 0,
            'Rmin' :  -5, 'Rmax' : 3,
            'cmap' : plt.cm.RdYlGn_r, }

#=================================2==============================================
#                            load data, select events
#================================================================================
eqCat.loadMatBin(  os.path.join( data_dir, file_in))

eqCat.selectEvents( dPar['a_Mc'][0], None, 'Mag')
#eqCat.selectEvents( tmin, tmax, 'Time')
print( 'no. of events after initial selection', eqCat.size())

# two ways to do the distance comp: 1 project into equal distance azimuthal , comp Cartersian distance in 3D
#                                   2 get surface distance from lon, lat (haversine), use pythagoras to include depth
eqCat.toCart_coordinates( projection = 'eqdc')

for i in range( dPar['a_Mc'].shape[0]):
    f_Mc =  dPar['a_Mc'][i]
    eta_0_file = '%s/%s_Mc_%.1f_eta_0.txt'%(data_dir, file_in, f_Mc)
    # load eta_0 value
    if os.path.isfile( eta_0_file):
        print( 'load eta_0 from file'),
        f_eta_0 = np.loadtxt( eta_0_file, dtype = float)
    else:
        print( 'could not find eta_0 file', eta_0_file, 'use value from dPar', dPar['eta_0'])
        f_eta_0 = dPar['eta_0']
    # cut below current completeness
    eqCatMc.copy( eqCat)
    eqCatMc.selectEvents( f_Mc, None, 'Mag')
    print( 'current catalog size: ',eqCatMc.size())
    
    # load nearest neighbor distances 
    NND_file = 'data/%s_NND_Mc_%.1f.mat'%( file_in.split('.')[0], f_Mc)
    dNND   = data_utils.loadmat(NND_file) #,  struct_as_record=True)

    #==================================3=============================================
    #                       compute re-scaled interevent times and distances
    #================================================================================ 
    catChild.copy( eqCatMc)
    catParent.copy( eqCatMc)
    catChild.selEventsFromID(    dNND['aEqID_c'], repeats = True)
    catParent.selEventsFromID(   dNND['aEqID_p'], repeats = True)
    print( 'size of offspring catalog', catChild.size(), 'size of parent cat', catParent.size())  

    # note that dictionary dPar here has to include 'b','D' and 'Mc'
    a_R, a_T = clustering.rescaled_t_r( catChild, catParent, {'b':dPar['b'], 'D':dPar['D'], 'Mc':f_Mc}, correct_co_located = True)
    RT_file = 'data/%s_RT_Mc_%.1f.mat'%( file_in.split('.')[0], f_Mc)
    scipy.io.savemat( RT_file, {'R' : a_R, 'T': a_T}, do_compression  = True)
    #==================================4==============================================================
    #                       T-R density plots
    #=================================================================================================
    a_Tbin = np.arange( dPar['Tmin'], dPar['Tmax']+2*dPar['binx'], dPar['binx'])
    a_Rbin = np.arange( dPar['Rmin'], dPar['Rmax']+2*dPar['biny'], dPar['biny'])
    XX, YY, ZZ = data_utils.density_2D( np.log10( a_T), np.log10( a_R), a_Tbin, a_Rbin, sigma = dPar['sigma'])
    
    fig1 = plt.figure(1, figsize= (8,10))
    ax = plt.subplot(111)
    ax.set_title( 'Nearest Neighbor Pairs in R-T')
    #------------------------------------------------------------------------------ 
    normZZ = ZZ*( dPar['binx']*dPar['biny']*eqCatMc.size())
    plot1 = ax.pcolormesh( XX, YY, normZZ, cmap=dPar['cmap'])
    cbar  = plt.colorbar(plot1, orientation = 'horizontal', shrink = .5, aspect = 20,)
    #ax.plot(  np.log10( a_T), np.log10( a_R), 'wo', ms = 1.5, alpha = .2)
    # plot eta_0 to divide clustered and background mode
    ax.plot( [dPar['Tmin'], dPar['Tmax']],  -np.array([dPar['Tmin'], dPar['Tmax']])+f_eta_0, '-', lw = 1.5, color = 'w' )
    ax.plot( [dPar['Tmin'], dPar['Tmax']],  -np.array([dPar['Tmin'], dPar['Tmax']])+f_eta_0,'--', lw = 1.5, color = '.5' )
    #-----------------------labels and legends-------------------------------------------------------
    #cbar.set_label( 'Event Pair Density [#ev./dRdT]') 
    cbar.set_label( 'Number of Event Pairs',labelpad=-40)
    ax.set_xlabel( 'Rescaled Time')
    ax.set_ylabel( 'Rescaled Distance') 
    ax.set_xlim( dPar['Tmin'], dPar['Tmax'])
    ax.set_ylim( dPar['Rmin'], dPar['Rmax'])

    #=================================5==============================================
    #                           save results
    #================================================================================
    print( 'plot saved in: ','plots/T_R_%s_Mc_%.1f.png'%( file_in.split('.')[0], f_Mc))
    fig1.savefig( 'plots/T_R_%s_Mc_%.1f.png'%( file_in.split('.')[0], f_Mc))
    plt.show()
     
    plt.clf()








