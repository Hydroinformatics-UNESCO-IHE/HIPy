# -*- coding: utf-8 -*-
"""
Data Save module
"""
import cPickle
import csv

def spkl(var, ID):
    '''
    Function used to save data from interpolation results into native python
    format (PKL)
    '''
    with open(ID +'.pkl','w') as DataFile:
        cPickle.dump(var,DataFile)
    print ID + ' Pickled'
    print ''
    return
    
def scsv(var, ID):
    '''
    Function used to save data from interpolation results into comma separated
    value (CSV)
    ''' 
    with open(ID +'.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                                quotechar=',', 
                                quoting=csv.QUOTE_MINIMAL)        
        for j in xrange(len(var)):
            spamwriter.writerow([var[j]])    
    print ID + ' CSV-ed'
    print ''    
    return