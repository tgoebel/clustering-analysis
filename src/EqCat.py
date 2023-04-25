#!usr/bin/python3.7
"""seismic catalog analysis class earthquake catalogs
-  data is stored in dictionary which can be extended without difficulties
as long as all vectors have the same length

- basic functionalities are focused on catalog I/O
  and initial processing (space, time, magnitude window selection) 

"""
import os
import numpy as np
import scipy.io #to writer and read mat bin
# the next line sets the path to PROJ LIB, should be found automatically for conda install
#-----------------my modules-----------------------------------------
#import ClusteringAnalysis.src.datetime_utils as dateTime
import src.datetime_utils as dateTime

#--------------------------------------------------------------------
class EqCat:
    """

    (1) 
    EqCat.data - type python dictionary
    e.g.:
    self.data = {       'N'          : , #event number 
                        'Time'       : np.array([]), # in decimal years
                        'Lon'        : np.array([]), #or lon
                        'Lat'        : np.array([]), #or lat
                        'Depth'      : np.array([]), #or depth
                        'Mag'        : np.array([]),

                           }  
    """
    def __init__( self, **kwargs ):
        """initiate data dictionary

        """
        self.data           = {}

        self.methods = [method_name for method_name in dir(self)
             if callable(getattr(self, method_name)) and method_name[0] != '_']

        """ input use kwargs to go from cartesian to GPS coordinates,
        tags can be accessed: sLoc1 - sLoc3 , last one is depth or self.sLoc3 """
#         if 'type' in kwargs.keys() and kwargs['type'] == 'GPS':
#             self.sLoc1, self.sLoc2, self.sLoc3 = 'Lon', 'Lat', 'Depth'
#         elif 'type' in kwargs.keys() and kwargs['type'] == 'Cart':
#             self.sLoc1, self.sLoc2, self.sLoc3 = 'X','Y','Z'
#         
        self.sLoc1, self.sLoc2, self.sLoc3 = 'Lon', 'Lat', 'Depth'
        self.sID = 'N'
        
    def copy(self, catalog ):
        """ deep copy of catalog object"""
        import copy
        try:
            for tag, vector in catalog.data.items():
                self.data[tag] = copy.copy( catalog.data[tag])
        except:
            for tag, vector in catalog.items():
                self.data[tag] = copy.copy( catalog[tag])

    #===========================================================================
    #                         import routines
    #===========================================================================
    def loadEqCat(self, file_in, catalogType, verbose=False, **kwargs):
        """ check what type of catalog and call correct function that handles import
        input: - file         - catalog filename
               - catalogType  = 'hs_reloc', focMech ... etc.
                              = 'WaldhauserReloc' - Waldhauser's selection of repeaters at Hayward
                              = 'hypoDD' - ID, lat, long, depth, x, y, z, x-error,y-error,z-error, yr, month, day, hour, minute, second, magnitude
               - kwargs['header']       - what is the line number of header info. of columns -> used for dic tags
               - kwargs['removeColumn'] - specify columns to be removed from original file prior to loading the file
                                           uses 'awk'
                                        - required since loadtxt assume all table entries to be floats
        
        TODO: --> biggest time sink is checking the date-time for every earthquake and converting it to dec. year --> vectorizing should help                                   
        
        return: create eqCat object with self.data = {'Time', 'Lon', 'Lat', 'Depth', 'Mag'}, 
                which are the standard dictionary tags 
        """
        #----------------check kwargs---------------------------------
        if 'header' in kwargs.keys() and kwargs['header'] is not None:
            header = kwargs['header']
        else:
            header = None
        if 'removeColumn' in kwargs.keys() and kwargs['removeColumn'] is not None:
            import src.data_utils as data_utils
            # remove columns and change file_name to copy of original file to keep the original
            file_in = data_utils.removeColumn( file_in, kwargs['removeColumn'])
        #-----------choose import routine------------------------------
        if catalogType == 'HS_reloc':
            if header is None:
                header = 1
            #TODO: get dic tag from file header
            headList = ['YR', 'MO', 'DY', 'HR', 'MN','SC', 'N', 'Lat','Lon','Depth', 'Mag', 'nPick', 'distSta', 'rms', 'd/n', 'rMeth', 'clID', 'nEvInCl',  'nlnk','err_h','err_z','rel_err_H', 'rel_err_Z']
            self.data = {}
            # 0-5 (datetime), 6(ID), 7 (lat), 8 (lon), 9 (depth), 10 (mag)
            mData = np.loadtxt(f"{file_in}", usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
            print( 'no of columns', mData[0].shape[0])
            print( 'no. of earthquakes', mData[:,0].shape[0])
            for l in range( mData[0].shape[0] ):
                self.data[headList[l]] = mData[:,l]
                
        elif catalogType == 'USGS':
            # 'time', 'latitude', 'longitude', 'depth', 'mag', 'magType', 'nst', 'gap', 'dmin', 'rms', 'net', 'id'
            #    0         1          2           3       4       5         6       7       8       9     10    11

            ###1###Date-time
            mDateTime     = np.genfromtxt( file_in, delimiter=(4,1,2,1,2,1,2,1,2,1,4),
                                       skip_header=1, usecols=(0,2,4,6,8,10)).T
            headDate = ['YR', 'MO', 'DY', 'HR', 'MN', 'SC']
            for i in range( len(headDate)):
                self.data[headDate[i]] = mDateTime[i]
            ###2### ID
            #mID = np.loadtxt( file_in, delimiter=',', skiprows=1, usecols=(10,11), dtype = str).T
            #self.data['ID'] = np.array([ int(mID[1,i].strip( mID[0,i])) for i in range( mID.shape[1])], dtype = int)
            self.data['ID']  = np.arange( len( self.data['YR']))
            ###3### location, magnitude, gap etc.
            header = ['Lat', 'Lon', 'Depth', 'Mag']#, 'Nst', 'Gap', 'Dmin', 'rms']
            mData = np.loadtxt( file_in, delimiter=',', skiprows=1,
                                usecols=(1,2,3,4),#,6,7,8,9),
                                dtype = float).T
            for i in range( len(header)):
                self.data[header[i]] = mData[i]

        elif catalogType == 'Kilauea':
            mData = np.loadtxt( file_in).T
            # :TODO convert np.array to python dictionary

        #convert date to decimal year
        self.data['Time'] = np.array([])
        for i in range( self.data['Mag'].shape[0] ):
            if verbose == True:
                print( i+1, 'out of', self.data['Mag'].shape[0])
            YR, MO, DY, HR, MN, SC = dateTime.checkDateTime( [self.data['YR'][i], self.data['MO'][i],self.data['DY'][i], self.data['HR'][i],self.data['MN'][i],self.data['SC'][i]])
            self.data['Time'] = np.append( self.data['Time'], 
                                           dateTime.dateTime2decYr( [YR, MO, DY, HR, MN, SC]))
        #sort catalog chronologically
        self.sortCatalog('Time')

        ##clean up
        if 'removeColumn' in kwargs.keys() and kwargs['removeColumn'] is not None:
            print( "delete: %s, than hit: y"%( file_in))
            removeFile = input( ' ')
            print( removeFile)
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
                raise( ValueError( error_str))
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
                raise( ValueError( error_str))
        #sel = np.arange( self.size(), dtype = int )[sel]
        if 'returnSel' in kwargs.keys() and kwargs['returnSel'] == True:
            return sel
        else:        
            self.selDicAll( sel) 

    def sortCatalog(self, tag, **kwargs):
        """sort catalog according to tag (string) e.g. Time, Mag, ....
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
            #    print( tag, 'does not have the right dimension: %i %i'%(vector.shape[0], sel.shape[0])
            #else:
            self.data[tag] = self.data[tag][sel]

    def selEventsFromID(self, a_ID, **kwargs):
        """ select events specified by list of IDs (self.data['N'])
        -----------------input

        kwargs:  repeats = True , if eqIDs are repeated keep them in catalog and maintain the same order
                default  = False, every earthquake is only ones in catalog, for several events with same ID keep only the first event               
        
        ----------------return:
        eq catalog that corresponds to vEqID """
        Nev= len( a_ID)     
        repeats = False
        if 'repeats' in kwargs.keys() and kwargs['repeats'] == True:
             a_sel  = np.ones(   Nev, dtype = int)
             v_i    = np.arange( self.size(), dtype = int)
             i = 0
             for currID in a_ID: # put one at location of ID match
                 sel_curr_ev = self.data['N']==int(currID)
                 if sel_curr_ev.sum() > 0:
                     a_sel[i] = int( v_i[sel_curr_ev][0])
                 i += 1
        else:
            a_sel = np.in1d( self.data['N'], a_ID, assume_unique=True)
        self.selDicAll( a_sel)

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
                self.data[key] = self.todict( self.data[key])

    def todict(self, matobj):
        '''
        A recursive function which constructs from matobjects nested dictionaries
        '''
        dData = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
                dData[strg] = self.todict(elem)
            else:
                dData[strg] = elem
        return dData

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
        l_tags = list( self.data.keys())
        for tag in l_tags:
            if tag[0] == '_':
                #print( 'remove', tag, self.data[tag]
                self.data.pop( tag, None)
            #else:
            #    print( tag, self.data[tag].shape[0]

    #======================================4==========================================
    #                            projections, rotations etc.
    #=================================================================================
    def toCart_coordinates(self, **kwargs):
        """
        :input
        **kwargs['projection']  =   'aeqd' - (default) azimuthal equidistant
                                    'eqdc' - equi distant conical projection
                                    'cyl'  - cynlidrical equidistant - not working
                'returnProjection' : True  - return basemap object
        use equidistant projection to convert lon, lat to X, Y coordinates
        :output catalog attributes:   - self.data['X'], self.data['Y'], self.data['Depth'] in km
                return True or basemap object, m
        
        """
        os.environ["PROJ_LIB"] = f"{os.environ['HOME']}/opt/anaconda3/share/proj"# adjust, comment out as needed
        from mpl_toolkits.basemap import Basemap
        projection = 'aeqd'
        if 'projection' in kwargs.keys() and kwargs['projection'] is not None:
            projection = kwargs['projection']
        from mpl_toolkits.basemap import Basemap
        xmin,xmax = self.data['Lon'].min(), self.data['Lon'].max()
        ymin,ymax = self.data['Lat'].min(), self.data['Lat'].max()

        # setup equi distance basemap.
        m = Basemap( llcrnrlat  =  ymin,urcrnrlat  =  ymax,
                     llcrnrlon  =  xmin,urcrnrlon  =  xmax,
                     projection = projection,lat_0=(ymin+ymax)*.5,lon_0=(xmin+xmax)*.5,
                     resolution = 'l')

        self.data['X'], self.data['Y'] = m( self.data['Lon'], self.data['Lat'])
        if projection == 'cyl':
            pass
        else:
            self.data['X'] *= 1e-3
            self.data['Y'] *= 1e-3
        if 'returnProjection' in kwargs.keys() and kwargs['returnProjection'] == True:
            return m
        else:
            return True

    #======================================5==========================================
    #                           shuffling, random catalog
    #=================================================================================
    def randomize_cat(self):
        """
        - create a randomized catalog with same average rate, no. of events and
          spatial extent as the initial catalog

        :return: - random Poissonian catalog, uniform spatial distribution
        """
        ## randomize event times