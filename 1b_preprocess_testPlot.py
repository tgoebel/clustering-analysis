#!python2.7
'''
Created on March 28th,  2019

- event selection
- plot earthquake catalog

TODO: - implement geo-referenced plotting with Basemap


@author: tgoebel
'''
#------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import os, sys
 
#------------------------------my modules-------------------------------------- 
import src.dataIO_utils as dataUtils
from src.EqCat import *

eqCat = EqCat( )

#print dir( dataUtils)
#=================================1==============================================
#                            dir, file, params
#================================================================================
dir_in = '%s/data/quakeData/SCSN/relocated'%( os.path.expanduser( '~'))
file_in= 'hs_1981_2011_all.mat'
#xmin, xmax = -122, -114
#ymin, ymax = 34, 38
Mmin, Mmax = 3, None
tmin, tmax = 1990, 2012


#=================================2==============================================
#                            load data, select events
#================================================================================
os.chdir( dir_in)
eqCat.loadMatBin(  file_in)
print 'total no. of events', eqCat.size()
eqCat.selectEvents( Mmin, Mmax, 'Mag')
eqCat.selectEvents( tmin, tmax, 'Time')
print 'no. of events after initial selection', eqCat.size()
#=================================3==============================================
#                          test plot TODO: use basemap
#================================================================================
plt.figure()
plt.scatter( eqCat.data['Lon'], eqCat.data['Lat'], s = np.exp( eqCat.data['Mag']), c = eqCat.data['Mag'], linewidth = 0)
plt.savefig( file_in.replace( 'mat', 'png'))







