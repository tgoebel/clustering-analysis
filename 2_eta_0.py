'''
Created on October 9, 2019

    1) shuffle event magnitudes and times
    2) computer NND values for all event pairs
    3) use 1st percentile as estimate for eta_0
        -> note that eta_0 is required to separate clustered
           and background events in the following analyses steps
    4) eta_0 is saved in /data directori
       file_out = [file_in]_Mc_[mc]_eta_0.txt
@author: tgoebel
'''
#------------------------------------------------------------------------------
import matplotlib as mpl
#mpl.use( 'Agg') # turn off interactive plot
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import os

#------------------------------my modules--------------------------------------
import src.clustering as clustering
import src.data_utils as data_utils
from   src.EqCat import *

eqCat   = EqCat( ) # original cat
ranCat  = EqCat()  # randomized, Poissonian catalog
eqCatMc = EqCat()  # catalog above completeness
np.random.seed( 123456)
#=================================1==============================================
#                            dir, file, params
#================================================================================
dir_in = 'data'
file_in= 'hs_1981_2011_all.mat'

#file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar  = {   'aMc'         :  np.array([3.0, 4.0]), #np.array( [2.0, 2.5, 3.0, 3.5]),
            # fractal dimension and b for eq. (1)
            'D'           : 1.6, # TODO: - these values should be constrained based on the data
            'b'           : 1.0, # use: https://github.com/tgoebel/magnitude-distribution for b-value

            # number of bootstraps for randomized catalogs
            'nBoot' : 100,
            #=================plotting==============
            'eta_binsize' :  .3,

            'cmin' : 1, 
            'xmin' : -13, 'xmax' : 0,
            ## R-T plot
            'binx' : .1, 'biny' : .1,# used for density and gaussian smoothing
            'sigma'   : None, #if None: default = n**(-1./(d+4)), or set Gaussian bandwidth
            'Tmin' :  -8, 'Tmax' : 0,
            'Rmin' :  -5, 'Rmax' : 3,
            'cmap'        : plt.cm.RdYlGn_r,
            'showPlot'    :  False,
          }

#================================================================================
#                      load data, event selection
#================================================================================
eqCat.loadMatBin(  os.path.join( dir_in, file_in))
print( 'total no. of events', eqCat.size())
eqCat.selectEvents( dPar['aMc'][0], None, 'Mag')
#eqCat.selectEvents( tmin, tmax, 'Time')
print( 'no. of events after initial selection', eqCat.size())
# project to equi-distant coordiante system for cartesian distances
eqCat.toCart_coordinates( projection = 'aeqd')
for f_Mc in dPar['aMc']:
    print( '-------------- current Mc:', f_Mc, '---------------------')
    # select magnitude range
    eqCatMc.copy( eqCat)
    eqCatMc.selectEvents( f_Mc, None, 'Mag')
    print( 'catalog size after MAG selection', eqCat.size())
    # this dictionary is used in module: clustering
    dConst = {'Mc' : f_Mc,
               'b' : dPar['b'],
               'D' : dPar['D']}

    #=============================2===================================================
    #                    randomize catalog
    #=================================================================================
    a_Eta_0 = np.zeros( dPar['nBoot'])
    for i_Bs in xrange( dPar['nBoot']):

        ranCat.copy( eqCatMc)
        ranCat.data['X']     = np.random.uniform( eqCatMc.data['X'].min(), eqCatMc.data['X'].max(), size = eqCatMc.size())
        ranCat.data['Y']     = np.random.uniform( eqCatMc.data['Y'].min(), eqCatMc.data['Y'].max(), size = eqCatMc.size())
        ranCat.data['Time']  = clustering.rand_rate_uni( eqCatMc.size(), eqCatMc.data['Time'].min(), eqCatMc.data['Time'].max())
        ranCat.sortCatalog( 'Time')
        #==================================3=============================================
        #                   compute space-time-magnitude distance, histogram
        #================================================================================
        dNND = clustering.NND_eta( ranCat, dConst,  M0 = 0,   correct_co_located = True,
                                   verbose = False)
        sel = dNND['aNND'] < 0
        a_Eta_0[i_Bs] = round( np.percentile( np.log10(dNND['aNND']), 1), 5)
        print( 'nBoot', i_Bs+1,'out of', dPar['nBoot'], 'eta 0 - 1st', np.percentile( np.log10(dNND['aNND']), 1))
        if dPar['showPlot'] == True: # plots to check if everything is working
            #=================================4==============================================
            #                          plot NND histogram
            #================================================================================
            plt.figure( 1, figsize = (10,5))
            ax = plt.axes( [.12, .12, .83, .83])
            ax.hist( np.log10( dNND['aNND']), np.arange( dPar['xmin'], dPar['xmax'], dPar['eta_binsize']),
                            color = '.5', label = 'Mc = %.1f'%( f_Mc), align = 'mid', rwidth=.9)
            ax.plot( [-5, -5], ax.get_ylim(), 'w-',  lw = 2, )
            ax.plot( [-5, -5], ax.get_ylim(), 'k--', lw = 2, )
            ax.plot( [a_Eta_0[i_Bs], a_Eta_0[i_Bs]], ax.get_ylim(), 'w-',  lw = 2, label = '$N_\mathrm{tot}$=%i'%( ranCat.size()))
            ax.plot( [a_Eta_0[i_Bs], a_Eta_0[i_Bs]], ax.get_ylim(), 'r--', lw = 2, label = '$N_\mathrm{cl}$=%i'%( dNND['aNND'][dNND['aNND']<1e-5].shape[0]))

            ax.legend( loc = 'upper left')
            ax.set_xlabel( 'NND, log$_{10} \eta$')
            ax.set_ylabel( 'Number of Events')
            ax.grid( 'on')
            ax.set_xlim( dPar['xmin'], dPar['xmax'])


            #==================================4==============================================================
            #                           T-R density plot
            #=================================================================================================
            catChild = EqCat()
            catParent= EqCat()
            catChild.copy(  ranCat)
            catParent.copy( ranCat)

            catChild.selEventsFromID(    dNND['aEqID_c'], repeats = True)
            catParent.selEventsFromID(   dNND['aEqID_p'], repeats = True)
            print( catChild.size(), catParent.size(), eqCatMc.size())
            a_R, a_T = clustering.rescaled_t_r( catChild, catParent, dConst, correct_co_located = True)

            a_Tbin = np.arange( dPar['Tmin'], dPar['Tmax']+2*dPar['binx'], dPar['binx'])
            a_Rbin = np.arange( dPar['Rmin'], dPar['Rmax']+2*dPar['biny'], dPar['biny'])
            a_log_T = np.log10( a_T)
            a_log_R = np.log10( a_R)
            XX, YY, ZZ = data_utils.density_2D( a_log_T, a_log_R, a_Tbin, a_Rbin, sigma = dPar['sigma'])

            plt.figure(2, figsize= (8,10))
            ax = plt.subplot(111)
            ax.set_title( 'Nearest Neighbor Pairs in R-T')
            #------------------------------------------------------------------------------
            normZZ = ZZ*( dPar['binx']*dPar['biny']*eqCatMc.size())
            plot1 = ax.pcolormesh( XX, YY, normZZ, cmap=dPar['cmap'])
            cbar  = plt.colorbar(plot1, orientation = 'horizontal', shrink = .5, aspect = 20,)
            #ax.plot(  np.log10( a_T), np.log10( a_R), 'wo', ms = 1.5, alpha = .2)
            # plot eta_0 to divide clustered and background mode
            ax.plot( [dPar['Tmin'], dPar['Tmax']],  -np.array([dPar['Tmin'], dPar['Tmax']])+a_Eta_0[i_Bs], '-', lw = 1.5, color = 'w' )
            ax.plot( [dPar['Tmin'], dPar['Tmax']],  -np.array([dPar['Tmin'], dPar['Tmax']])+a_Eta_0[i_Bs],'--', lw = 1.5, color = '.5' )
            #-----------------------labels and legends-------------------------------------------------------
            #cbar.set_label( 'Event Pair Density [#ev./dRdT]')
            cbar.set_label( 'Number of Event Pairs',labelpad=-40)
            ax.set_xlabel( 'Rescaled Time')
            ax.set_ylabel( 'Rescaled Distance')
            ax.set_xlim( dPar['Tmin'], dPar['Tmax'])
            ax.set_ylim( dPar['Rmin'], dPar['Rmax'])

            plt.show()
    #=================================3==============================================
    #                            save results
    #================================================================================
    f_eta_0 = a_Eta_0.mean()
    print 'medium eta_0', a_Eta_0.mean()
    file_out = '%s/%s_Mc_%.1f_eta_0.txt'%(dir_in, file_in, f_Mc)
    np.savetxt( file_out, np.array([f_eta_0]), fmt = '%10.3f', header='eta_0')
    print 'save results', file_out
    scipy.io.savemat(file_out.replace('txt','mat'),
                     {'eta_0': f_eta_0, 'eta_BS' : a_Eta_0,}, do_compression=True)














