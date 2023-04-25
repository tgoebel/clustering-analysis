#!python3.7
'''
Created on March 28th,  2019

- load Hauksson, Shearer 2011 eq catalog from scec data center, alterntive catalogs
- save as .mat binary for fast data I/O
- note that the original catalog is not provided and has to be downloaded from the web:
  https://scedc.caltech.edu/eq-catalogs/altcatalogs.html

@author: tgoebel
'''
#------------------------------------------------------------------------------
import os

#------------------------------my modules-------------------------------------- 
from src.EqCat import EqCat

eqCat = EqCat( )

#=================================1==============================================
#                            dir, file, params
#================================================================================
# change to local dir where eq. catalogs are saved
# the original catalog can be found here: https://scedc.caltech.edu/research-tools/altcatalogs.html
dir_in = 'data'
file_in= 'hs_1981_2011_all.txt'

#=================================2==============================================
#                            load data
#================================================================================
import numpy as np
# 0-5 (datetime), 6(ID), 7 (lat), 8 (lon), 9 (depth), 10 (mag)
mData = np.loadtxt( f"{dir_in}/{file_in}", usecols=(0,1,2,3,4,5,6,7,8,9, 10)).T
print( mData.shape)

eqCat.loadEqCat( f"{dir_in}/{file_in}", 'HS_reloc')

print( 'total no. of events: ', eqCat.size())
print( sorted( eqCat.data.keys()))
#=================================3==============================================
#                     test plot and save to .mat binary
#================================================================================
eqCat.saveMatBin( file_in.replace( 'txt', 'mat'))
newEqCat = EqCat( )
newEqCat.loadMatBin( file_in.replace( 'txt', 'mat'))
print( newEqCat.size())
print( sorted( newEqCat.data.keys()))







