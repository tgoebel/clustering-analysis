#!usr/bin/python2.7
# -*- coding: utf-8 -*-
"""

    helper functions for easier file handling (mainly ASCII)
    and data I/O, density estimates, 2D Gaussian smoothing etc.


@author tgoebel - UC Santa Cruz
"""
import os
import numpy as np
import scipy.io
#================================================================================
#                           data I/O
#================================================================================   
def removeColumn( file_in, lCol):
    """
    remove all columns specified in lCol
    1) create duplicate file called 'dummy_file.txt' in cwd
    2) remove column using awk
    3) return file_name of dublicate
    """
    # example syntax to remove three columns
    #os.system( "awk '{\$24=""; \$25=""; \$26=""; print(}' in_file.txt > out_file.txt")
    lStr = []
    for col in lCol:
        lStr.append( "$%s=\"\"; "%( col))
    tmp_file    = 'dummy_file.txt'
    command_str = "awk '{ %s print(}' %s > %s"%( ''.join( lStr), file_in, tmp_file)
    os.system( command_str)           
    return tmp_file

def loadmat(filename, verbose = False):
    '''
    this function should be called instead of directly calling scipy.io.loadmat 
    which is used within the method
        (1) - filters dictionary tags 
        (2) - properly recovers python dictionaries
               from mat files. check dic tag which are still mat-objects
        (3) - correct arrays of the form: np.array([[  1, 2, 3]]) to np.array([  1, 2, 3]), squeeze_me=True
        (4) - can handle 'nested' variables in matlab where variable contain several structures
    
    '''
    data = scipy.io.loadmat(filename, struct_as_record=True, squeeze_me=True)
    data = _check_keys(data)
    for tag in list( data.keys()):
        if tag[0] == '_':
            if verbose == True:
                print( 'remove', tag, data[tag])
            data.pop( tag)
    return data

def _check_keys( dData):
    '''
    checks if entries in dictionary are mat-objects. If yes
    to dict is called to change them to nested dictionaries
    '''
    for key in dData:
        if isinstance(dData[key], scipy.io.matlab.mio5_params.mat_struct):
            dData[key] = _todict(dData[key])
    return dData        

def _todict( matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dData = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
            dData[strg] = _todict(elem)
        else:
            dData[strg] = elem
    return dData
#================================================================================
#                           density estimates and smoothing
#================================================================================   
def density_2D( x, y, x_bin, y_bin, **kwargs):
    """
        2D, smoothed event density for point cloud with coordinates x,y
        uses method: scipy.stats.kde.gaussian_kde
        see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
    :input     x,y            - dataset
               x_bin, y_bin   - binned x and y vectors
    
    
        kwargs['sigma'] - specify gaussian smoothing kernel ('bw_method' in scipy.stats.kde)
                          default: =  n**( -1./(d+3)) adapted scott rule for slightly tighter bandiwdth
                        -  'scott' 
                              sigma = n**( -1./(d+4)), d- number of dimensions, n - number of data points
                        - 'silverman'
                              sigma = (n * (d + 2) / 4.)**(-1. / (d + 4))
                        - float( )   = set Gaussian Bandwidth directlty
                                        
                           
    return XX, YY, ZZ - 2D binned x and y coordinates and density for each cell
    """
    from scipy.stats import kde
    n,d = x.shape[0],2
    sigma = n**( -1./(d+2.5))
    if 'sigma' in kwargs.keys() and kwargs['sigma'] is not None:
        sigma = kwargs['sigma']
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    fct_Gauss2D = kde.gaussian_kde( np.array([x,y]), bw_method = sigma)
    # meshgrid of x and y coordinates
    XX,YY   = np.meshgrid( x_bin, y_bin)
    ZZ = fct_Gauss2D( np.vstack([XX.flatten(), YY.flatten()])).reshape( XX.shape)
    dx, dy = x_bin[1] - x_bin[0], y_bin[1] - y_bin[0]
    # check if integral is ~ zero, better: use midepoint method
    print( 'check if integral ~1', round(ZZ.sum()*( dx*dy),3)) #ZZ[ZZ>0].mean()*(XX.max()-XX.min())*(YY.max()-YY.min()))
    return XX-.5*dx, YY-.5*dy, ZZ

#================================================================================
#                          dictionary processing
#================================================================================ 
def copyDic( dic):
    """ create a copy of dic"""
    import copy
    dCopy = {}
    for tag in dic.keys():
        dCopy[tag] = copy.copy( dic[tag])
    return dCopy

def selectDataRange(dicOri, min, max, tag, **kwargs):
    """
    select data within given range, set min = None or max =None for only lower or upper bound
    """
    dic = copyDic(dicOri)
    if 'includeBoundaryEvents' in kwargs.keys() and kwargs['includeBoundaryEvents'] == True:
        if min == None or max == None:
            error_str = 'both boundaries have to be set to include boundary events'
            raise( ValueError( error_str))
        else:
            sel = np.logical_and( dic[tag] >= float(min), dic[tag] <= float(max ) )          
    if max == None:
        sel = dic[tag] > float(min)
    elif min == None:
        sel = dic[tag] < max
    else:
        sel = np.logical_and( dic[tag] > float(min), dic[tag] < float(max) )
    sel = np.arange( dic[tag].shape[0], dtype = int )[sel]
    if 'returnSel' in kwargs.keys() and kwargs['returnSel'] == True:
        return sel
    else:        
        return selDicAll(dic, sel, **kwargs) 


def selDicAll(dic, curr_sel, **kwargs):
    """apply boolean vector to entire data
    e.g. for sorting or cutting ... """
    newDic = {}
    for tag, vector in dic.items():
        newDic[tag] = dic[tag][curr_sel]
    return newDic



