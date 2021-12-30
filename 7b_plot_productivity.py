'''
Created on August 16, 2016

- plot No. of events for each mainshock
- average number of events per mainshock
- fit N ave and Mms --> alpha exponent

@author: tgoebel
'''
import matplotlib as mpl
#mpl.use( 'Agg')

import matplotlib.pyplot as plt
import numpy as np
import os


#=================================1==============================================
#          dir, file, params
#================================================================================
data_dir   = 'data'
plot_dir   = 'plots'
file_in    = 'hs_1981_2011_all.mat'

dPar  = {
            'magRound'    : 1, # for binning
            'a_Mc'        : np.array([2.5, 3.0, 4.0]),

            #=================plotting==============
            'alpha'       : 1.0, # for plotting demonstration
            'xmin' :  2,  'xmax' : 8,
            'ymin' : 0.1, 'ymax' : 1e4,
            'plotFormat' : 'png',
            }

#================================================================================
#                            load file with no. of aftershocks
#================================================================================

iMc = 0
for f_Mc in dPar['a_Mc']:
    file_prod = '%s/%s_Nas_MS_Mc_%.1f.txt'%(data_dir, file_in.split('.')[0], f_Mc)#, dPar['magRound'])

    m_N_as = np.loadtxt( file_prod).T
    print( 'total no. of mainshock', m_N_as[0].shape[0])
    print( 'total no. of AS', m_N_as[1].sum())
    print( 'total no. of FS', m_N_as[2].sum())
    #=================================2==========================================================================
    #                           count ave. no. of aftershocks per MS magnitude
    #============================================================================================================
    aMag_round= np.around( m_N_as[0], dPar['magRound'])
    aMag_bin  = np.array( sorted(np.unique( aMag_round)))
    aAveNo_AS = np.ones( len( aMag_bin))*np.nan
    aNo_Fam   = np.zeros( len( aMag_bin)) # total number of families within mag bin
    aNo_AS20  = np.zeros( len( aMag_bin))
    aNo_AS80  = np.zeros( len( aMag_bin))

    i = 0
    for curr_mag in aMag_bin:
        selMag       = curr_mag == aMag_round
        aAveNo_AS[i] = m_N_as[1][selMag].mean()
        if selMag.sum() > 0:
            aNo_AS20[i]  = np.percentile( m_N_as[1][selMag], 20)
            aNo_AS80[i]  = np.percentile( m_N_as[1][selMag], 80)
        aNo_Fam[i]   = selMag.sum()
        print( curr_mag, 'mean N-AS', round(aAveNo_AS[i],2),  aNo_AS20[i],aNo_AS80[i], 'no. of fam', aNo_Fam[i])

        i += 1

    #=================================3==========================================================================
    #                           plot productivity law
    #============================================================================================================
    plt.figure(1, figsize=(8,6))
    ax = plt.axes([.14,.12,.78,.83])#pPlot.createFigureSquare(1)
    ax.semilogy( m_N_as[0],   m_N_as[1],     'o',  ms = 6, mew =0, mfc = '.7', alpha = .2 )
    #ax.errorbar( aMag_bin,  aAveNo_AS, yerr=[np.zeros(aMag_bin.shape[0]), aNo_AS80-aAveNo_AS],
    #             fmt = 'o', ecolor = 'k', elinewidth=.7,capsize=2.5, mec = 'k', ms = 8, mew = 1, mfc = 'w')
    ax.errorbar( aMag_bin,  aAveNo_AS, yerr=[aAveNo_AS-aNo_AS20, aNo_AS80-aAveNo_AS],
                 fmt = 'o', ecolor = 'k', elinewidth=.7,capsize=2.5, mec = 'k', ms = 8, mew = 1, mfc = 'w')

    #-------------------------exponential - estimate-----------------------------------------------------
    mag_fit    = aMag_bin[10] # force fit through this point
    f_no_AS_pl = aAveNo_AS[aMag_bin == mag_fit]
    preFac     = np.log10( f_no_AS_pl) - dPar['alpha']*mag_fit
    a_N_hat    = 10**( dPar['alpha']*aMag_bin + preFac)
    ax.semilogy( aMag_bin, a_N_hat, 'w-')
    ax.semilogy( aMag_bin, a_N_hat, '-', color = 'r', lw = 2, label = 'exp = %.1f'%( np.round( dPar['alpha'],1)))

    #-------------------------------labels, limits etc.-----------------------------------------------
    ax.set_xlim( dPar['xmin'], dPar['xmax'])
    ax.set_ylim( dPar['ymin'], dPar['ymax'])
    ax.set_xlabel( 'Mainshock Magnitude')
    ax.set_ylabel( 'Number of Aftershocks')
    ax.legend( loc = 'upper left', frameon = False)

    plt.savefig( '%s/%s_Mc_%.1f_ASprod.%s'%(plot_dir, file_in.split('.')[0], f_Mc, dPar['plotFormat']))
    
    plt.show()
    plt.clf()

    iMc += 1










