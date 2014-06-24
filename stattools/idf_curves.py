#!/usr/bin/env python
# -*- coding: utf-8 -*-
r'''Intensity Frequency Duration Curves:
====================================

This module contains functions to calculate the IDF curves.
The results can be adjusted also to fit a general purpose equation in the form
i = (c * T^m)/(Td^e + f)


TO DO : 
    - This could be handled as a class
    - Make idf Class and wrap funcionality inside


Gonzalo Andrés PEÑA CASTELLANOS
UNESCO-IHE - MWI-SE 2010-2012
e-mail: goanpeca@gmail.com

'''
import os
import numpy
import scipy.optimize
import scipy.stats
#from scipy import optimize, stats

import matplotlib.pyplot

import wg_tool
import wg_io

def main():
    pass
    iori = r'D:\Dropbox\module_14\data\02_preprocessed\02_pickled\2000_tokyo.data'
    idir = r'D:\Dropbox\module_14\data\03_processed\01_present\04_gen\tokyo_agg_1_6_12_48_group_m_weight_boot_stat_CVh_Rlh.para - Copy' 
    idir = r'D:\Dropbox\module_14\data\03_processed\01_present\04_gen\tokyo_agg_1_6_12_48_group_m_weight_boot_stat_CVh_Rlh.para' 
    duration = [1, 2, 6, 12, 24, 48, 96]
    return_period = [50]
    plot_pos = 'w'
    accumulate = False
    
#    print ev1_multiple(idir, duration, return_period, plot_pos, accumulate)
    

    idirs = {'a1b': r'D:\Dropbox\module_14\data\03_processed\01_present\04_gen\tokyo_agg_1_6_12_48_group_m_weight_boot_stat_CVh_Rlh.para',
            'a2': r'D:\Dropbox\module_14\data\03_processed\01_present\04_gen\tokyo_agg_1_6_12_48_group_m_weight_boot_stat_CVh_Rlh.para',
            'b1': r'D:\Dropbox\module_14\data\03_processed\01_present\04_gen\tokyo_agg_1_6_12_48_group_m_weight_boot_stat_CVh_Rlh.para'}
        
    odir = ''
    duration = [1,2,6,12,24]    
    idf_scenarios(iori, idirs, odir, duration, return_period)

def idf_scenarios(iori, idirs, odir, duration, return_period):
    r'''
    '''
    pass
    
    rp = 0
    dic_ori = ev1(iori, duration, return_period)
    idf_data_ori = dic_ori['idf_data']
    adjust_coef_ori = adjust_idf(idf_data_ori, duration, return_period)    
    range_agg_plot = numpy.linspace(1, duration[-1], 1000)
    ori_adjusted = evaluate_idf(adjust_coef_ori[0], range_agg_plot, return_period)[rp]
   
    scenario_adjusted = {}
    
    x = range_agg_plot 
    
    for k in idirs:
        ifile = idirs[k]
        dic_temp = ev1_multiple(ifile, duration, return_period)
        idf_data_temp = dic_temp['idf_data']
        adjust_coef_temp = adjust_idf(idf_data_temp, duration, return_period)    
        scenario_adjusted[k] = evaluate_idf(adjust_coef_temp[0], range_agg_plot, return_period)[rp] 

        if 'sres' in k:
            year_max = dic_temp['year_max']
            year_min = dic_temp['year_min']
        
        scenario_adjusted['year_max'] = year_max
        scenario_adjusted['year_min'] = year_min

    scenario_adjusted['ori'] = ori_adjusted
    scenario_adjusted['x'] = x


    
#    matplotlib.pyplot.plot(x, scenario_adjusted['a1b'])    
#    matplotlib.pyplot.plot(x, scenario_adjusted['a2'])    
#    matplotlib.pyplot.plot(x, scenario_adjusted['b1'])    
#    matplotlib.pyplot.plot(x, scenario_adjusted['ori'])        
#    matplotlib.pyplot.show()
    
    return scenario_adjusted
    
    
#    idf_data = dic['idf_data']
#    rain_adjusted = wg_idf.adjust_idf(idf_data, aggregation_adjust, return_period)    
#    dic_2 = wg_idf.ev1(ifile, aggregation_compare, return_period)
#    X = rain_adjusted[0]
#    idf_data_compare = dic_2['idf_data']
#        
#    max_1 = numpy.max(aggregation_adjust)
#    max_2 = numpy.max(aggregation_compare)
#    max_agg = max(max_1, max_2)
#    range_agg_plot = numpy.linspace(1, max_agg, 1000)
#    idf_adjusted = wg_idf.evaluate_idf(rain_adjusted[0], range_agg_plot, return_period)   
#
#    year_max = dic['year_max']
#    year_min = dic['year_min']
#    period = str(year_min) + '-' + str(year_max)

    

def ev1_multiple(idir, duration, return_period, plot_pos='w', 
                 accumulate=False):
    r'''Function to calculate the extreme value distribution of rainfall 
    for different rainfall duration and return periods. This can be then used
    to produce intensity-duration-frequency curves.
    
    Parameters
    ----------
    rain : 
        numpy array
    year : 
        numpy array
    duration : 
        list
    return_period :
        
    plot_pos : 
        P = (i-b)/(n + 1 - 2*b), where n is the number of items, i is the index
        for the sorted list and b is given 
        
        'w' Weibul        0.0
        'c' Chegodayev    0.30
        't' Turkey        0.33333
        'b' Blom          0.375  
        'g' Gringorten    0.44
        'h' Hazen         0.5
    
    
    Returns
    -------
    A dictionary idf_data
    
    '''
    ifiles = os.listdir(idir)
        
    plot_position = {'w': 0.0, 
                     'c': 0.30,
                     't': 0.333333333333333333,
                     'b': 0.375,
                     'g': 0.440,
                     'h': 0.5}

    if str(plot_pos) in plot_position:
        b = plot_position[plot_pos]
    else:
        print('Selected ploting position incorrect. Using Weibul by default')
        b = plot_position[plot_pos]
        # Plot position: (i-b)/(n + 1 - 2*b))

    # Get basic data from file one
    ifile = os.path.join(idir, ifiles[0])
    dic = wg_io.load(ifile)
    year = dic['year']
    year_max = numpy.int(numpy.amax(year))
    year_min = numpy.int(numpy.amin(year)) 
    
    rain_max_list = []
    
    for i in ifiles:
        print i
        ifile = os.path.join(idir, i)
        dic = wg_io.load(ifile)
        rain = dic['rainfall'] 
        year = dic['year']
        
        dic_max = wg_tool.get_max(rain, year, duration, accumulate=accumulate)
        rain_max = dic_max['rain']
        rain_max_list.append(rain_max)

    # now create median
    size = dic_max['size']
    prob = [(i-b)/(size + 1.0-2.0*b) for i in range(1, size + 1)]
    prob = numpy.array(prob)
    ln_ = -numpy.log(-numpy.log(prob))      
    
    coefficient = []    
    for step in range(len(duration)):
        rain_median = scipy.median(numpy.array(rain_max_list), axis=0)
        (m, b, r, tt, stderr) = scipy.stats.linregress(ln_, rain_median[step])
        coefficient.append([m, b, r])
    
    ## Create the array to use in ploting later
    idf = []
    r_p = numpy.array(return_period)*1.0
    ln_calc = -numpy.log(-numpy.log(1.0 - (1.0 / r_p)))
        
    ## mbr stands for Slope(m) Y-intercept(b) Rcoefficient(r)
    for mbr in coefficient:
        i_calc = ln_calc*mbr[0] + mbr[1]  
        idf.append(i_calc) 
                                
    idf_data = numpy.transpose(numpy.array(idf))

    return {'coefficients': coefficient, ## [slope, y-int, R]
            'rain': rain_max_list ,
            'probability': prob,
            'ln': ln_,
            'idf_data': idf_data,
            'duration': duration,
            'year_max': year_max,
            'year_min': year_min
            }

#    
#    
#    x = -numpy.log(-numpy.log(prob))  
#    x_adj = numpy.arange(-1.5,3,0.1)
#    
#    for step in range(len(duration)):
#        fig = matplotlib.pyplot.figure()
#        fig.suptitle('Rainfall duration: ' + str(duration[step]))        
#        ax = fig.add_subplot(111)
#        for t in rain_max_list:
#            ax.plot(x, t[step], 'k.')
#
#        y_median = scipy.median(numpy.array(rain_max_list), axis=0)
#        y_mean = numpy.mean(numpy.array(rain_max_list), axis=0)
#        
#        (m, b, r, tt, stderr) = scipy.stats.linregress(x, y_median[step])
#        y_adj = m*x_adj + b        
#        ax.plot(x_adj, y_adj, 'g-')
#        ax.plot(x, y_median[step], 'b-')
#        ax.plot(x, y_mean[step], 'r-')
#    
#        matplotlib.pyplot.show()
#

      
#    coefficient = []
#    rain_max_list = []
#    probability = []
#    ln = []
#    #year = numpy.transpose(year)   
#    
#    for step in duration:
#        rain_year_max = []   
#        for y in range(year_min, (year_max+1)):           
#            rain_per_year = rain[(year==y)]
#
###!         ##This is to comply with the rolling window function in case there
#            ## are years with no data.
#            if rain_per_year.shape[0] > step: 
#                rain_year_m_ave = movmean(rain_per_year, step, 1, accumulate)    
#                temp_calc = numpy.nanmax(rain_year_m_ave)
#            
###!         ## If the moving average produces NaN when calculating the
#            ## nanmax ignores the NaN. 
#            ## Is this appropriate in every case?
#                if not numpy.isnan(temp_calc): 
#                    rain_year_max.append(temp_calc)
#        
#        size = len(rain_year_max)
#        sorting =  numpy.sort(rain_year_max)
#
#        print size
#
#        ## Only do this once
#        if step == duration[0]:
#            # ploting position probability
#            prob = [(i-b)/(size+1.0-2.0*b) for i in range(1, size + 1)]
#            prob = numpy.array(prob)
#            #fr = 1 - prob
#            ## Extreme value 1  TO UPDATE
#            
#            ln_ = -numpy.log(-numpy.log(prob))  
#            probability.append(prob)
#            ln.append(ln_)
#            
#        (m, b, r, tt, stderr) = scipy.stats.linregress(ln_, sorting)
#        
#        coefficient.append([m, b, r])
#        rain_max_list.append(sorting)
#        
#    ## Create the array to use in ploting later
#    idf = []
#    r_p = numpy.array(return_period)*1.0
#    ln_calc = -numpy.log(-numpy.log(1.0 - (1.0 / r_p)))
#        
#    ## mbr stands for Slope(m) Y-intercept(b) Rcoefficient(r)
#    for mbr in coefficient:
#        i_calc = ln_calc*mbr[0] + mbr[1]  
#        idf.append(i_calc) 
#                                
#    idf_data = numpy.transpose(numpy.array(idf))
#
#    return {'coefficients': coefficient, ## [slope, y-int, R]
#            'rain': rain_max_list ,
#            'probability': probability,
#            'ln': ln,
#            'idf_data': idf_data,
#            'duration': duration,
#            'year_max': year_max,
#            'year_min': year_min            
#            }


def adjust_idf(idf_data, aggregation, return_period):
    r"""
    Adjusts the calculated ev1 data  to the following general equation
    using a  
        i = (c * T^m)/(Td^e + f), where
        
        i : Rainfall Intensity
        T : the return period
        Td : Duration of storm
        c, m, e, f : Shape coefficients   
   
    
    Parameters
    ----------
    idf_data :
        adas
    aggregation :
        asdasd
    return_period : 
        asdasd
        
    Returns
    -------
    """
    ## Initial conditions   
    ##             c  m  e  f
    Xo = numpy.array([1, 1, 1, 1], dtype=numpy.float)

    ## [Min, Max]   c              m            e             f
    bounds = numpy.array([[None, None], [None, None], [None, None], [None, None]])
    
    opt_res = scipy.optimize.fmin_l_bfgs_b(objective_function, Xo, 
                args=(idf_data, aggregation, return_period), 
                approx_grad=True, bounds=bounds, maxfun=3000)    
                
    return opt_res
                   
def objective_function(Xo, idf_data, duration, return_period):
    r"""
    i = (c * T^m)/(Td^e + f)
    Evaluates    
    
    Parameters
    ----------
    X :
        adas
    aggregation :
        asdasd
    return_period : 
        asdasd
        
    Returns
    -------
    """       
    (c, m, e, f) = Xo        
        
    T = numpy.transpose(numpy.array([return_period]))
    Td = numpy.array(duration)
   
    pre_T = T**m
    pre_i = c / ((Td**e) + f)
    i = (numpy.ones((T.shape[0], Td.shape[0])) * pre_i) * pre_T
      
    error =  numpy.sum(((1 - idf_data/i)**2 + (1 - i/idf_data)**2))     
    return error

def evaluate_idf(X, duration, return_period):
    r"""
    Evaluates    
    
    Parameters
    ----------
    X :
        adas
    duration :
        asdasd
    return_period : 
        asdasd
        
    Returns
    -------
    
    
    """
    ## i = (c * T^m)/(Td^e + f)       
    (c, m, e, f) = X    
        
    T = numpy.transpose(numpy.array([return_period]))
    Td = numpy.array(duration)

    pre_T = T**m
    pre_i = c / ((Td**e) + f)
    i = (numpy.ones((T.shape[0], Td.shape[0])) * pre_i) * pre_T
    
    return i
    
def ev1(ifile, duration, return_period, plot_pos='w', accumulate=False):
    r'''Function to calculate the extreme value distribution of rainfall 
    for different rainfall duration and return periods. This can be then used
    to produce intensity-duration-frequency curves.
    
    Parameters
    ----------
    rain : 
        numpy array
    year : 
        numpy array
    duration : 
        list
    return_period :
        
    plot_pos : 
        P = (i-b)/(n + 1 - 2*b), where n is the number of items, i is the index
        for the sorted list and b is given 
        
        'w' Weibul        0.0
        'c' Chegodayev    0.30
        't' Turkey        0.33333
        'b' Blom          0.375  
        'g' Gringorten    0.44
        'h' Hazen         0.5
    
    
    Returns
    -------
    A dictionary idf_data
    
    '''
    dic = wg_io.load(ifile)
    rain = dic['rainfall'] 
    year = dic['year']
    
    plot_position = {'w': 0.0, 
                     'c': 0.30,
                     't': 0.333333333333333333,
                     'b': 0.375,
                     'g': 0.440,
                     'h': 0.5}

    if str(plot_pos) in plot_position:
        b = plot_position[plot_pos]
    else:
        print('Selected ploting position incorrect. Using Weibul by default')
        b = plot_position[plot_pos]

        # Plot position: (i-b)/(n + 1 - 2*b))

    # Common use functions    
    #movmean = tool.movmean
    movmean = wg_tool.movnanmean
    
    year_max = numpy.int(numpy.amax(year))
    year_min = numpy.int(numpy.amin(year)) 
      
    coefficient = []
    rain_max_list = []
    probability = []
    ln = []
    #year = numpy.transpose(year)   
    
    for step in duration:
        rain_year_max = []   
        for y in range(year_min, (year_max+1)):           
            rain_per_year = rain[(year==y)]

##!         ##This is to comply with the rolling window function in case there
            ## are years with no data.
            if rain_per_year.shape[0] > step: 
                rain_year_m_ave = movmean(rain_per_year, step, 1, accumulate)    
                temp_calc = numpy.nanmax(rain_year_m_ave)
            
##!         ## If the moving average produces NaN when calculating the
            ## nanmax ignores the NaN. 
            ## Is this appropriate in every case?
                if not numpy.isnan(temp_calc): 
                    rain_year_max.append(temp_calc)
        
        size = len(rain_year_max)
        sorting =  numpy.sort(rain_year_max)

        print size

        ## Only do this once
        if step == duration[0]:
            # ploting position probability
            prob = [(i-b)/(size+1.0-2.0*b) for i in range(1, size + 1)]
            prob = numpy.array(prob)
            #fr = 1 - prob
            ## Extreme value 1  TO UPDATE
            
            ln_ = -numpy.log(-numpy.log(prob))  
            probability.append(prob)
            ln.append(ln_)
            
        (m, b, r, tt, stderr) = scipy.stats.linregress(ln_, sorting)
        
        coefficient.append([m, b, r])
        rain_max_list.append(sorting)
        
    ## Create the array to use in ploting later
    idf = []
    r_p = numpy.array(return_period)*1.0
    ln_calc = -numpy.log(-numpy.log(1.0 - (1.0 / r_p)))
        
    ## mbr stands for Slope(m) Y-intercept(b) Rcoefficient(r)
    for mbr in coefficient:
        i_calc = ln_calc*mbr[0] + mbr[1]  
        idf.append(i_calc) 
                                
    idf_data = numpy.transpose(numpy.array(idf))

    return {'coefficients': coefficient, ## [slope, y-int, R]
            'rain': rain_max_list ,
            'probability': probability,
            'ln': ln,
            'idf_data': idf_data,
            'duration': duration,
            'year_max': year_max,
            'year_min': year_min            
            }

if __name__ == '__main__':
    main()