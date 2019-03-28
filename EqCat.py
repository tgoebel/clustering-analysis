#!usr/bin/python2.7
"""seismic catalog analysis class earthquake catalogs
-  data is stored in dictionary which can be extended without difficulties
as long as all vectors have the same length

- basic functionalities are focused on catalog I/O
  and initial processing (space, time, magnitude window selection) 

"""
from __future__ import division # include true division 
import os
import numpy as np
import scipy.io #to writer and read mat bin




class EqCat:
    """

    (1) 
    EqCat.data - type python dictionary
    e.g.:
    self.data = {       'N'          : , #event number 
                        'Time'       : np.array([]), # indecimal years
                        'Lon'        : np.array([]), #or lon
                        'Lat'        : np.array([]), #or lat
                        'Depth'      : np.array([]), #or depth
                        'MAG'        : np.array([]),

                           }  
    """
    def __init__( self, **kwargs ):
        """initiate data dictionary

        """
        self.data           = {}

        """ input use kwargs to go from cartesian to GPS coordinates,
        tags can be accessed: sLoc1 - sLoc3 , last one is depth or self.sLoc3 """
#         if 'type' in kwargs.keys() and kwargs['type'] == 'GPS':
#             self.sLoc1, self.sLoc2, self.sLoc3 = 'Lon', 'Lat', 'Depth'
#         elif 'type' in kwargs.keys() and kwargs['type'] == 'Cart':
#             self.sLoc1, self.sLoc2, self.sLoc3 = 'X','Y','Z'
#         
        self.sLoc1, self.sLoc2, self.sLoc3 = 'Lon', 'Lat', 'Depth'
    #===========================================================================
    #                         import routines
    #===========================================================================

    def loadEqCat(self, file_in, catalogType, **kwargs):
        """ check what type of catalog and call correct function that handles import
        input: - file         - catalog filename
               - catalogType  = 'hs_reloc', focMech ... etc.
                              = 'WaldhauserReloc' - Waldhauser's selection of repeaters at Hayward
                              = 'hypoDD' - ID, lat, long, depth, x, y, z, x-error,y-error,z-error, yr, month, day, hour, minute, second, magnitude
               - kwargs['header']       - what is the line number of header info. of columns -> used for dic tags
               - kwargs['removeColumn'] - specify columns to be removed from original file prior to loading the file
                                           uses 'awk'
        
        TODO: --> biggest time sink is checking the date-time for every earthquake and converting it to dec. year --> vectorizing should help                                   
        
        return: create eqCat object with self.data = {'Time', 'Lon', 'Lat', 'Depth', 'Mag'}, 
                which are the standard dictionary tags 
        """
        import date_time_utils as dt_utils
        #----------------check kwargs---------------------------------
        if 'header' in kwargs.keys() and kwargs['header'] is not None:
            header = kwargs['header']
        else:
            header = None
        if 'removeColumn' in kwargs.keys() and kwargs['removeColumn'] is not None:
            import dataIO_utils
            # remove columns and change file_name to copy of original file to keep the original
            file_in = dataIO_utils.removeColumn( file_in, kwargs['removeColumn'])
        #-----------choose import routine------------------------------
        if catalogType == 'HS_reloc':
            if header is None:
                header = 1
            #TODO: get dic tag from file header
            headList = ['YR', 'MO', 'DY', 'HR', 'MN','SC', 'N', 'Lat','Lon','Depth', 'Mag', 'nPick', 'distSta', 'rms', 'd/n', 'rMeth', 'clID', 'nEvInCl',  'nlnk','err_h','err_z','rel_err_H', 'rel_err_Z']
            self.data = {}           
            mData = np.loadtxt( file_in, skiprows = header)
            #mData = mData.T
            print 'no of columns', mData[0].shape[0]
            print 'no. of earthquakes', mData[:,0].shape[0]
            for l in xrange( mData[0].shape[0] ):
                self.data[headList[l]] = mData[:,l]
                
        elif catalogType == 'WaldhauserReloc':
            mData = np.loadtxt( file_in)
            mData = mData.T
            # DATE        TIME         LAT          LON         DEP      EH1     EH2    AZ    EZ    MAG       ID
            #1984  1  1  1 19 11.320    36.08787   -120.22890   10.964   0.028   0.015  88   0.028  1.8    1109386
            vStrHeader = ['Time', 'Lat', 'Lon',  'Depth', 'EH1',   'EH2', 'AZ', 'EZ', 'MAG',   'N' ]
            # compute decimal year
            vTime = np.array([])
            for i in xrange( mData[0].shape[0]):
                vTime  = np.append( vTime, toDecimalYear([mData[0][i],mData[1][i],mData[2][i],mData[3][i],mData[4][i],mData[5][i]]))
            self.data['Time'] = vTime

            for i in xrange( len(vStrHeader)-1):
                self.data[vStrHeader[i+1]] = mData[i+6]
#         
#         elif catalogType == 'HS_reloc':
#             if headerLine is None:
#                 headerLine = 4
#             headList = ['YR', 'MO', 'DY', 'HR', 'MN','SC', 'N', 'Lat','Lon','Depth', 'MAG', 'nPick', 'distSta', 'rms', 'd/n',  'clID', 'nEvInCl',  'nlnk','err_h','err_z','rel_err_H', 'rel_err_Z','type', 'lMeth', 'poly']
#             self.data = {}
#             for tag in headList[0:-3]:
#                 if tag == 'id':
#                     self.data[tag] = np.array([], dtype = int)
#                 else:
#                     self.data[tag] = np.array([])
#             file_obj = open( file, 'r')
#             file_obj.next()
#             l = 0
#             for line in file_obj:
#                 row = line.split()
#                 if len(row) > 2:
#                     n = 0
#                     for tag in headList[0:-3]:
#                         if tag == 'id':
#                             self.data[tag] = np.append( self.data[tag],   int( row[n]) )
#                         else:
#                             self.data[tag] = np.append( self.data[tag],   float( row[n]) )
#                         
#                         n = n + 1
#                     print self.data['YR'][l],self.data['MN'][l],self.data['DY'][l]
#                     l = l + 1
#             file_obj.close()
#             


        #convert date to decimal year
        self.data['Time'] = np.array([])
        for i in xrange( self.data['Mag'].shape[0] ):
            print i+1, 'out of', self.data['Mag'].shape[0]
            YR, MO, DY, HR, MN, SC = dt_utils.checkDateTime( [self.data['YR'][i], self.data['MO'][i],self.data['DY'][i], self.data['HR'][i],self.data['MN'][i],self.data['SC'][i]])
            self.data['Time'] = np.append( self.data['Time'], 
                                           dt_utils.datetime2decYr( [YR, MO, DY, HR, MN, SC]))
        self.data.pop( 'YR')
        self.data.pop( 'MO')
        self.data.pop( 'DY')
        self.data.pop( 'HR')
        self.data.pop( 'MN')
        self.data.pop( 'SC')
        #sort catalog chronologically
        self.sortCatalog('Time')

        ##clean up
        if 'removeColumn' in kwargs.keys() and kwargs['removeColumn'] is not None:
            print "remove: %s, than hit: y   "%( file_in)
            removeFile = raw_input( ' ')
            print removeFile
            if os.path.isfile( file_in) and removeFile == 'y':
                os.system( "rm %s"%( file_in))

    #======================================2==========================================
    #                            basic processing and catalog event selection
    #=================================================================================
    def size(self):
        if 'Time' in self.data.keys():
            return len( self.data['Time'])
        else:
            return None
    
    
    def selectEvents(self, min, max, tag, **kwargs):
        """
        returns events with time, coordinates, rel.Magnitude that corresponds to a certain time frame
        -cut catalog includes lower bound (min) but excludes upper bound (max)
        input:  min, max = window of events
                min      - can be set to string for columns that contain strings, e.g. type, magType  etc.
                if min is not a string:
                min = None, select only events below max
                max = None, select only events above min
                tag can be 'Time' or magnitude , location, Mw... depending on dictionary
        kwargs: includeBoundaryEvents = True; include events with times equal to min and max otherwise
                                              include only lower boundary (min event)  
                returnSel             = returns IDs of selected events (type np.array([], int))
        
        example: selectEvents( 3, 5, 'Mag', includeBoundaryEvents = True) - all events between 3 and 5 including M=3 and M=5 events
                 selectEvents( 3, None, 'Mag')  - everything above M=3 excluding M=3 events
                 selectEvents( 4, None, 'Mag') and then selectEvents( 'w', None, 'MagType') - all Mws above Mw = 4
                
        """
        if 'includeBoundaryEvents' in kwargs.keys() and kwargs['includeBoundaryEvents'] == True:
            if min == None or max == None:
                error_str = 'both boundaries have to be set to include boundary events'
                raise ValueError, error_str
            else:
                sel = np.logical_and( self.data[tag] >= float(min), self.data[tag] <= float(max ) )          
        else:
            if isinstance( min, str ):
                #str columns, e.g magType ..
                sel = [i for i, x in enumerate( self.data[tag] ) if x == min]  
            elif isinstance( min, (int, float) ) or min == None:
                if max == None:
                    sel = self.data[tag] >= float(min)
                elif min == None:
                    sel = self.data[tag] < max
                else:
                    sel = np.logical_and( self.data[tag] >= float(min), self.data[tag] < float(max) )
            else:
                error_str = 'unknown input min = %s'%(min)
                raise ValueError, error_str
        #sel = np.arange( self.size(), dtype = int )[sel]
        if 'returnSel' in kwargs.keys() and kwargs['returnSel'] == True:
            return sel
        else:        
            self.selDicAll( sel) 

    def sortCatalog(self, tag, **kwargs):
        """sort catalog according to tag (string) e.g. Time, mag, avePol....
        kwargs: beginWithBiggest = True , sort beginning with Biggest value
                returnSel        = return boolean """
        #get boolean vector for sorting
        vSortBool  = self.data[tag].ravel().argsort()
        if 'beginWithBiggest' in kwargs.keys() and kwargs['beginWithBiggest'] == True:
            if 'returnSel' in kwargs.keys() and kwargs['returnSel'] == True:
                return vSortBool
            else:
                self.selDicAll( vSortBool[::-1])
        else:
            if 'returnSel' in kwargs.keys() and kwargs['returnSel'] == True:
                return vSortBool
            else:
                self.selDicAll( vSortBool)

    def selDicAll(self, sel):
        """apply boolean vector to entire data
        e.g. for sorting or cutting ... """
        for tag, vector in self.data.items(): #loop through all entries (tag = vector name, vector = entries)
            # for NND analysis first event is missing (orphan), so sel.shape = vector.shape - 1
            #if sel.shape[0] != vector.shape[0]:
            #    print tag, 'does not have the right dimension: %i %i'%(vector.shape[0], sel.shape[0])
            #else:
            self.data[tag] = self.data[tag][sel]
            
    #======================================3==========================================
    #                            .mat binary load save
    #=================================================================================
    def check_keys(self, ):
        '''
        checks if entries in dictionary are mat-objects. If yes
        to dict is called to change them to nested dictionaries
        '''
        for key in self.data:
            if isinstance(self.data[key], scipy.io.matlab.mio5_params.mat_struct):
                self.data[key] = _todict( self.data[key])

    def saveMatBin(self, file):
        """save dic to bin file"""
        #scipy.io.savemat(file, self.data, appendmat=False, format = '4', oned_as = 'row' ,  do_compression = True)
        scipy.io.savemat(file, self.data, appendmat=True, format = '5',do_compression = True )
    

    def loadMatBin(self, filename):
        '''
        this function should be called instead of direct scipy.io.loadmat
        as it helps with additional non-variable tags in python dictionaries from .mat files


        --> can handle 'nested' variables in matlab where variable contain several structures
        '''
        
        self.data = scipy.io.loadmat(filename,struct_as_record=False, squeeze_me=True)
        self.check_keys( )
        for tag in self.data.keys():
            if tag[0] == '_':
                #print 'remove', tag, self.data[tag]
                self.data.pop( tag)
            #else:
            #    print tag, self.data[tag].shape[0]











        
        