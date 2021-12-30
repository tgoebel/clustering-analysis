#!usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
convert year month day hour min sec to decimal year and vs. verse

@author tgoebel - UC Santa Cruz
"""
from __future__ import division
import numpy as np

import time, datetime, calendar
from datetime import datetime as dt


def mo_to_sec( value):
    return value*(aveDyYr()/12)*24*3600

def sec_to_mo( value):
    return value/((aveDyMo())*24*3600)

def dy_to_sec( value):
    return value*24*3600

def sec_to_dy( value):
    return value/(24*3600)


def aveDyYr():
    """ how many days in a year"""
    return 365 + 1/4 - 1/100 + 1/400

def aveDyMo(): 
    """ how many days in a month """
    return aveDyYr()/12

def checkDateTime( dateTime):
    """ check that hour != 24, MN != 60, SC != 60 """
    YR, MO, DY, HR, MN, SC = int(dateTime[0]), int(dateTime[1]), int(dateTime[2]), int(dateTime[3]),int(dateTime[4]), float(dateTime[5])
    if isinstance( YR, (float, int)):
        if SC < 0:
            SC = 0
        elif SC - 60 >= 0:
            MN += int((SC/60))
            SC -= 60*int( (SC/60))
        if MN < 0:
            MN = 0
        elif MN - 60 >= 0:
            HR += int((MN/60))
            MN -= 60.*int( (MN/60.))
        if HR < 0:
            HR = 0
        elif HR - 24 >= 0:
            HR = 23
            MN = 59
            SC = 59.999
    elif isinstance( YR, (np.ndarray)):
        #set all values below zero to zero
        sel = SC < 0
        SC[sel] = 0
        sel = MN < 0
        MN[sel] = 0
        sel = HR < 0
        HR[sel] = 0        
        #set 60 to zero and 24 to 23.59.59.99
        sel = abs(SC - 60) < 1e-6
        SC[sel] = 0
        MN[sel] = MN[sel] + 1 
        sel = 60 - MN < 1e-6
        MN[sel] = 0
        HR[sel] = HR[sel] + 1   
        sel = 24 - HR < 1e-6
        HR[sel] = 23
        MN[sel] = 59
        SC[sel] = 59.99
    return YR, MO, DY, HR, MN, SC


#------------------------------------------------------------------------------ 
#                        date-time conversions
#------------------------------------------------------------------------------
def dateTime2decYr( datetime_in, **kwargs ):
    """
    input: datetime_in = array containing time columns year - second
                   out = date in decimal year
    """
    try:
        o_dt = datetime.datetime( int( datetime_in[0] ), int( datetime_in[1] ), int( datetime_in[2] ), int( datetime_in[3] ), int( datetime_in[4] ), int( round( datetime_in[5])-1e-3))
    except:
        error_msg = "datetime array not valid - %s; check if date and time is correct, e.g. no SC > 60.." % datetime_in
        raise( ValueError, error_msg)
    time_sc = o_dt.hour*3600 + o_dt.minute*60 + o_dt.second
    # get no. of day within current year between 0 to 364 and ad time in seconds
    dayOfYear_seconds = ( o_dt.timetuple().tm_yday - 1 ) * 86400.0 + time_sc
    if calendar.isleap( o_dt.year):
        year_fraction = dayOfYear_seconds / ( 86400.0 * 366 )
    else:
        year_fraction = dayOfYear_seconds / ( 86400.0 * 365 )
    # dec year = current year + day_time (in dec year)
    return o_dt.year + year_fraction

def decYr2datetime( decimalYear ):
    """
    convert decimal year to year/month/day... 
    """
    year = np.floor( decimalYear)
    rest = decimalYear-year
    
    if year%4 == 0: # leap year
	    ndays = 366    
	    feb = 29
    else:
	    ndays = 365
	    feb = 28
    decDay = rest * ndays 

    if decDay >= 0 and decDay <= 31:
	    month = 1
	    day  = np.ceil( decDay )
	    rest = (decDay) -np.floor( decDay )
    elif decDay >= 0 and decDay <= 31+feb:
	    month = 2
	    day = np.ceil( decDay-  31 )
	    rest = 1 -(day - (decDay - 31 ))
    elif decDay >= 31+feb and decDay <= 2*31+feb:
	    month = 3
	    day = np.ceil( decDay- (31+feb ))
	    rest = 1 -(day - (decDay -(31+feb )))
    elif decDay >= 2*31+feb and decDay <= 3*31+feb-1:
	    month = 4
	    day = np.ceil( decDay- (2*31+feb))
	    rest = 1 -(day - (decDay -(2*31+feb)))
    elif decDay >= 3*31+feb-1 and decDay <= 4*31+feb-1:
	    month = 5
	    day = np.ceil( decDay -(3*31+feb-1) )
	    rest = 1 -(day - (decDay -(3*31+feb-1)))
    elif decDay >= 4*31+feb-1 and decDay <= 5*31+feb-2:
	    month = 6
	    day = np.ceil( decDay-(4*31+feb-1))
	    rest = 1 -(day - (decDay -(4*31+feb-1)))
    elif decDay >= 5*31+feb-2 and decDay <= 6*31+feb-2:
	    month = 7
	    day = np.ceil( decDay-(5*31+feb-2) )
	    rest = 1 -(day - (decDay -(5*31+feb-2)))
    elif decDay >= 6*31+feb-2 and decDay <= 7*31+feb-2:
	    month = 8
	    day = np.ceil( decDay -(6*31+feb-2))
	    rest = 1 -(day - (decDay -(6*31+feb-2)))
    elif decDay >= 7*31+feb-2 and decDay <= 8*31+feb-3:
	    month = 9
	    day = np.ceil( decDay -(7*31+feb-2) )
	    rest = 1 -(day - (decDay -(7*31+feb-2)))
    elif decDay >= 8*31+feb-3 and decDay <= 9*31+feb-3:
	    month = 10
	    day = np.ceil( decDay -(8*31+feb-3))
	    rest = 1 -(day - (decDay -(8*31+feb-3)))
    elif decDay >= 9*31+feb-3 and decDay <= 10*31+feb-4:
	    month = 11
	    day = np.ceil( decDay -(9*31+feb-3))
	    rest = 1 -(day - (decDay -(9*31+feb-3)))
    elif decDay >= 10*31+feb-4 and decDay <= 11*31+feb-4:
	    month = 12
	    day = np.ceil( decDay -(10*31+feb-4))
	    rest = 1 -(day - (decDay -(10*31+feb-4)))
    else:
	    print( 'wrong input decimal year')
    hour   = np.floor( rest * 24 )
    rest   = 24*rest-hour
    minute = np.floor( rest * 60 )
    rest   = 60*rest-minute
    second =  rest * 60     
    if  day == 0: # for int decimal years
        day = 1
    try:  
        return [int(year[0]), int(month), int(day[0]), int(hour[0]), int(minute[0]), second[0]]
    except:
        return [int(year), int(month), int(day), int(hour), int(minute), second]
