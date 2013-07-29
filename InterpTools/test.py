# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 17:27:42 2013

@author: chaco3
"""

import pyximport
pyximport.install()
import IDW

import numpy


a = numpy.array([[1,2],[3,4]])

x = [1,2]
y = [1,2]
val = [5,3]
xt = [1,4,5,9]
yt = [3,4,6,8]

k = IDW.Interp(x,y,val,xt,yt)
print k