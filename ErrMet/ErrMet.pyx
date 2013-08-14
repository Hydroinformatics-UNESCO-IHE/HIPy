# -*- coding: utf-8 -*-
"""
Performance Functions
x - calculated value
y - recorded value
q - Quality tag (0-1)
"""
import numpy
# constants import 

# Or something similar,.... like importing a module of contants 
#or something
ERROR_CODE = 9999

# Example
def is_valid(x, y, *q):
    '''
    x 
    y
    q
    '''
    if len(x) != len(y):
        print 'Calculated and recorded series do not match'
        return False
    return True
    
def MAE2(x,y,*q):
    if is_valid(x,y,q): 
        cdef float F = (1./sum(q))*numpy.abs(numpy.sum(x-y))
    else: 
        return -9999
    return F
    

def DataValid(x,y,*q):
    if len(x) != len(y):
        print 'Calculated and recorded series do not match'
        return -9999
    try:
        q
    except NameError:
        q = numpy.ones(len(x))
    try:
        numpy.sum(x)
    except ValueError:
        print 'Calculated data might contain non-numerical values'
        return -9999
    try:
        numpy.sum(y)
    except ValueError:
        print 'Calculated data might contain non-numerical values'
        return -9999
    
    x = numpy.array(x)
    y = numpy.array(y)    
    return x,y,q

def RMSE(x,y,*q):
    x,y,q = DataValid(x,y,*q)
    if x == -9999:
        return -9999    
    Erro = numpy.square(x-y)*q    
    cdef float F = numpy.sqrt(1.*numpy.sum(Erro)/(numpy.sum(q))
    return F
    
def Bias(x,y,*q):
    x,y,q = DataValid(x,y,*q)
    if x == -9999:
        return -9999
    cdef float F = (1./sum(q))*numpy.sum(x-y)
    return F

def MAE(x,y,*q):
    x,y,q = DataValid(x,y,*q)
    if x == -9999:
        return -9999
    cdef float F = (1./sum(q))*numpy.abs(numpy.sum(x-y))
    return F
    
def MSE(x,y,*q):
    x,y,q = DataValid(x,y,*q)
    if x == -9999:
        return -9999    
    Erro = numpy.square(x-y)*q    
    cdef float F = numpy.sum(Erro)/numpy.sum(q)
    return F   
    
def PercVol(x,y,*q):
    x,y,q = DataValid(x,y,*q)
    if x == -9999:
        return -9999
    cdef float sy = numpy.sum(y)
    cdef float F = (numpy.sum(x)-sy)/sy
    return F

def NSE(x,y,*q):
    x,y,q = DataValid(x,y,*q)
    if x == -9999:
        return -9999
    a = numpy.sum(numpy.power(x-y,2)*q)
    b = numpy.sum(numpy.power(y-numpy.average(y),2)*q)
    cdef float F = 1.0 - a/b
    return F
