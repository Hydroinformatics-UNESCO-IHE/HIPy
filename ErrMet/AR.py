# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:13:24 2014

@author: chaco3

Arma Model
"""

import scipy.stats
import GetData
import numpy
import scipy.optimize
Consult = 'TESTDATA'

def HarmonicTrend(Consult,k):
    
    ## to get data from database
    Data = GetData(Consult)
    
    # k  is number of harmonics
    init_p = numpy.zeros(3*k)    
    x_values = numpy.array(range(len(Data)))
    
    #Function to fit
    def ArFun (p):
        y = []
        for i in xrange(len(p)/3):
            y = p(3*i)*(numpy.sin(x_values+p(3*i+1)) + numpy.cos(x_values+p(3*i+2))) 
        return y
    
    def ArFit(p):
        HarmSerie = HarFun(p)
        Erro = numpy.square(HarmSerie-Data)    
        RMSE = numpy.sqrt(1.*numpy.sum(Erro)/(len(HarmSerie)))
        return RMSE
    ## Make exponential regression

    Opt_p = scipy.optimize.fmin(HarFit,init_p)
    
    ## detrended series
    DetSer = Data - HarFun(Opt_p)
    return DetSer