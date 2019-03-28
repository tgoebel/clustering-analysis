#!usr/bin/python2.7
# -*- coding: utf-8 -*-
"""

    helper functions for easier file handling (mainly ASCII)
    and data I/O


@author tgoebel - UC Santa Cruz
"""

import os

def removeColumn( file_in, lCol):
    """
    remove all columns specified in lCol
    1) create duplicate file called 'dummy_file.txt' in cwd
    2) remove column using awk
    3) return file_name of dublicate
    """
    # example syntax to remove three columns
    #os.system( "awk '{\$24=""; \$25=""; \$26=""; print}' in_file.txt > out_file.txt")
    lStr = []
    for col in lCol:
        lStr.append( "$%s=\"\"; "%( col))
    tmp_file    = 'dummy_file.txt'
    command_str = "awk '{ %s print}' %s > %s"%( ''.join( lStr), file_in, tmp_file)
    os.system( command_str)           
    return tmp_file