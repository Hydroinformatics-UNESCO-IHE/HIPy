#!/usr/bin/env python
# -*- coding: utf-8 -*-
r'''
Input / Output of time series:
=============================

Reads a txt or mat file that contains info for the project.
In the case of a text file the default format is

    yyyy-mm-dd hh:mm    variable

where the delimiter is a tab.

Gonzalo Andrés PEÑA CASTELLANOS
UNESCO-IHE - MWI-SE 2010-2012
e-mail: goanpeca@gmail.com
'''
import os
import gzip
import cPickle
import datetime

import numpy
import scipy
import matplotlib.dates

def main():
    pass

def saveprojectmat(file_name, project_info=None, result_NSRP=None):
    r"""
    Function to save the different dictionaries in a matlab style file (*.mat)
    for later retrieval or use of data
    
    Parameters
    ----------
    file_name : String
        Path to file
    project_info : list of strings
        general information of the project. [project_info, comments]
    results_NRSP : list of list containing the parameters.
        [['lambda'], ['beta'], ['mu_c'], ['eta'], ['alpha'], ['theta']        
          
    Returns
    -------
    
    """
    ## Save the parameters of the NSRP inside a dictionary
    names = ['lambda_', 'beta', 'mu_c', 'eta', 'alpha', 'theta']        
    if result_NSRP==None:
        result_NSRP=[[],[],[],[],[],[]]
    dic_para = dict(zip(names, result_NSRP))

    ## define the dictionary structure for a project file
    names = ['project_name', 'comments']        
    if project_info==None:
        project_info = ['untitled', '']    
    dic_project_info = dict(zip(names, project_info))

    ## save the final file
    names = ['project_info', 'NSRP_parameters']
    project = [dic_project_info, dic_para]
    dic_project = dict(zip(names, project))
        
    savemat(file_name, dic_project)        

def savemat(filename, object_, oned_as='column'):
    r"""
    Parameters
    ----------
    
    Returns
    -------
    """
    # Pending Include some if functionality to take version into account
    np_ver = [int(i) for i in numpy.version.version.split('.')]
    sp_ver = [int(i) for i in scipy.version.version.split('.')]
    np_req = [1, 3, 0] 
    sp_req = [0, 7, 1]

    scipy.io.savemat(filename, object_, oned_as='column', do_compression=True)
    # elif :
    #     pass
    #sp.io.savemat(filename, object_, oned_as=oned_as) ## does not work with scipy i have in ubuntu 0.7.0
    # elif :
    #     pass
    #sp.io.savemat(filename, object_)  
    return True


def loadmat(file_name):
    r"""
    
    """
    dic = scipy.io.loadmat(file_name, struct_as_record=True)
    return dic

def loadtxt(file_name, variable_name, date_format='%Y-%m-%d %H:%M', 
            row_start=0, missing_value=numpy.NaN, filling_value=numpy.NaN,
            delimiter='\t', save_as=None):
    r"""
    Parameters
    ----------
    file_name :
        -
    variable_name :
        -
    date_format : 
        '%Y-%m-%d %H:%M'
    row_start : 
        0
    missing_value : 
        np.NaN
    filling_value : 
        np.NaN
    delimiter : 
        '\t'    
            
    Returns
    -------
    Dictionary containing the variable and the date, day, year, month,
    hour minute and day of year.
    """    
    f = open(file_name)   
    lines = [l for l in f.readlines() if l.strip()] # Clean&remove empty lines
    f.close()

    size = len(lines) - row_start
    
    dtype_int = numpy.int  
    
    ## Data holders for date, variable and the broken date
    data = numpy.zeros(size)
    date = numpy.zeros(size)
    year = numpy.zeros(size, dtype=dtype_int)
    month = numpy.zeros(size, dtype=dtype_int)
    day = numpy.zeros(size, dtype=dtype_int)
    hour = numpy.zeros(size, dtype=dtype_int)
    minute = numpy.zeros(size, dtype=dtype_int)
    doy = numpy.zeros(size, dtype=dtype_int)
    
    for i in range(size):
        columns = lines[i].strip().split(delimiter)
        datevalue = datetime.datetime.strptime(columns[0], date_format )
        
        date[i] = numpy.float(matplotlib.dates.date2num(datevalue))       

        ## Check for missing values and fill them accordingly
        if columns[1] == missing_value:
            data[i] = filling_value
        else:
            data[i] = numpy.float(columns[1])

        ## Break components of date
        year[i] = datevalue.year
        month[i] = datevalue.month
        day[i] = datevalue.day
        hour[i] = datevalue.hour
        minute[i] = datevalue.minute
        doy[i] = datevalue.timetuple().tm_yday

        #print year[i], month[i], day[i]

#    dic = {'date' : np.transpose(np.array([date])), 
#           variable_name : np.transpose(np.array([data])),
#           'year' : np.transpose(np.array([year])),
#           'month' : np.transpose(np.array([month])),
#           'day' : np.transpose(np.array([day])),           
#           'hour' : np.transpose(np.array([hour])),           
#           'minute' : np.transpose(np.array([minute])),
#           'doy' : np.transpose(np.array([doy]))}

    dic = {'date' : date, 
           variable_name : data,
           'year' : year,
           'month' : month,
           'day' : day,           
           'hour' : hour,           
           'minute' : minute,
           'doy' : doy}
              
    if save_as <> None:
        save(dic, save_as)
        
    return dic
    
def save(object_, filename, protocol=-1):
    """Saves a compressed object to disk.
    
    Parameters
    ----------
    object_ : Python object
        Any Python object that is picklable
    filename : String
        File to save the file
    protocol : Int
        Defines which protocol to use, 0 text, 1 binary
        -1 means use highest protocol available.
    
    Returns
    -------
    None
    """
    dirbase = os.path.dirname(filename)
    
    if dirbase <> '':
        if not(os.path.isdir(dirbase)):
            os.makedirs(dirbase)        
    
    if os.path.exists(filename):
        os.remove(filename)

    file_ = gzip.GzipFile(filename, 'wb')
    file_.write(cPickle.dumps(object_, protocol))
    file_.close()

def load(filename):
    """Loads a compressed object from disk
    
    Parameters
    ----------
    filename : String
    Name of the file to uncompress and unpickle    
    
    Returns
    -------
    Python object
    """
    file_ = gzip.GzipFile(filename, 'rb')
    buffer_ = ""
    while True:
            data = file_.read()
            if data == "":
                    break
            buffer_ += data
    object_ = cPickle.loads(buffer_)
    file_.close()
    return object_        
    
def parallel(func, input_arg):
    r'''Convenience function to make use of parallel calculation using 
    python-pp
    
    Parameters
    ----------
    func : 
        TODO
    input_arg :
        TODO

    Returns
    -------
    '''
    import pp    

    jobs = []
    j = 1

    # Tuple containing the modules that will be used in the excecuted code
    input_modules = ('gzip','cPickle','datetime',
                     'numpy','scipy','matplotlib.dates')

    #--------------------------------------------------------------------------   
    # Creates jobserver with automatically detected number of workers
    #--------------------------------------------------------------------------
    ppservers = ()
    job_server = pp.Server(ppservers=ppservers, secret='')
    workers = job_server.get_ncpus()
    
    print("\nStarting pp with " + str(workers) + " workers")

    for i in input_arg:
        jobs.append(job_server.submit(eval(func), i, (), input_modules))
    
    print('and ' + str(len(jobs)) + ' jobs!\n')

    for job in jobs:
        print('Job: ' + str(j))
        res = job()
        print('\n')
        log(res)
        j += 1
        
    job_server.print_stats()

def log(message, show=True):
    r'''Very simple logger to keep track of running scripts
    
    Parameters
    ----------
    message : String
        TODO
    show : Boolean
        TODO
    
    Returns
    -------
    '''
    script_name = __file__
    log_name = script_name.replace('.py', '_log.txt')
        
    if os.path.exists(log_name):
        f = open(log_name, 'a')
    else:
        f = open(log_name, 'w')

    msg = str('--------------------------\n' + str(datetime.datetime.now()) + 
              '\t\n' + str(message) + '\n.\n')
    
    f.write(msg)
    f.close
    
    if show:
        print msg
    
    return msg
        
if __name__ == '__main__':
    main()