#!python2.7
'''
Created on March 28th,  2019

- load Hauksson, Shearer 2011 eq catalog from scec data center, alterntive catalogs
- save as .mat binary for fast data I/O

@author: tgoebel
'''
#------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import os


#------------------------------my modules-------------------------------------- 
from src.EqCat import *

eqCat = EqCat( )

#=================================1==============================================
#                            dir, file, params
#================================================================================
dir_in = 'data'
file_in= 'hs_1981_2011_all.txt'
#file_in = 'test_HS.txt'

#=================================2==============================================
#                            load data
#================================================================================
os.chdir( dir_in)
eqCat.loadEqCat(  file_in, 'HS_reloc', removeColumn=[24,25,26])

print( 'total no. of events: ', eqCat.size())
print( sorted( eqCat.data.keys()))
#=================================3==============================================
#                     test plot and save to .mat binary
#================================================================================
eqCat.saveMatBin( file_in.replace( 'txt', 'mat'))
newEqCat = EqCat( )
newEqCat.loadMatBin( file_in.replace( 'txt', 'mat'))
print newEqCat.size()
print sorted( newEqCat.data.keys())







