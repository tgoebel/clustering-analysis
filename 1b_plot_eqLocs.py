#!python3.7
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
import os

os.environ["PROJ_LIB"] = f"{os.environ['HOME']}/opt/anaconda3/share/proj"
from mpl_toolkits.basemap import Basemap
#------------------------------my modules-------------------------------------- 

from src.EqCat import EqCat

eqCat = EqCat( )

#print( dir( dataUtils)
#=================================1==============================================
#                            dir, file, params
#================================================================================
dir_in = 'data'
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
print(  'total no. of events', eqCat.size())
eqCat.selectEvents( Mmin, Mmax, 'Mag')
eqCat.selectEvents( tmin, tmax, 'Time')
print( 'no. of events after initial selection', eqCat.size())
#=================================3==============================================
#                          test plot T
#================================================================================
projection = 'cyl'
xmin,xmax = eqCat.data['Lon'].min(), eqCat.data['Lon'].max()
ymin,ymax = eqCat.data['Lat'].min(), eqCat.data['Lat'].max()

# setup equi distance basemap.
m = Basemap( llcrnrlat  =  ymin,urcrnrlat  =  ymax,
             llcrnrlon  =  xmin,urcrnrlon  =  xmax,
             projection = projection,lat_0=(ymin+ymax)*.5,lon_0=(xmin+xmax)*.5,
             resolution = 'l')
m.drawstates( linewidth = 1)
m.drawcoastlines( linewidth= 2)
a_x, a_y = m( eqCat.data['Lon'], eqCat.data['Lat'])
m.plot( a_x, a_y, 'ko', ms = 1)
sel7 = eqCat.data['Mag'] >=7
m.plot( a_x[sel7], a_y[sel7], 'ro', ms = 8, mew= 1.5, mfc = 'none')


m.drawmeridians( np.linspace( int(xmin), xmax, 4),labels=[False,False,False,True],
                 fontsize = 12, fmt = '%.1f')
m.drawparallels( np.linspace( int(ymin), ymax, 4),labels=[True,False,False,False],
                 fontsize = 12, fmt = '%.2f')

plt.savefig( file_in.replace( 'mat', 'png'))
plt.show()






