#!/usr/bin/env python
# -*- coding: utf-8 -*-
r'''Stochastic Rainfall Generator:
==============================

Gonzalo Andrés PEÑA CASTELLANOS
UNESCO-IHE - MWI-SE 2010-2012
e-mail: goanpeca@gmail.com

Based on AWE-GEN v1.0:
http://www-personal.umich.edu/~ivanov/HYDROWIT/Models.html
----------------------------------------------------------
Department of Civil and Environmental Engineering, University of Michigan, 
Ann Arbor, Michigan, USA.

Department of Civil and Environmental Engineering, Massachusetts Institute 
of Technology, Cambridge, Massachusetts, USA.

Dipartimento di Ingegneria Civile e Ambientale, Università degli Studi di 
Firenze,  Firenze, Italy.

Institut für Geoökologie, TU Braunschweig, Langer Kamp, Braunschweig, 
Germany.
'''
import datetime
import calendar 

import numpy
import scipy.special
import scipy.integrate
import scipy.optimize
import matplotlib.dates

import wg_tool
import wg_downscale
import wg_io

__version__ = "0.1"
__author__ = "Gonzalo Andrés PEÑA CASTELLANOS"
__license__ = "GNU"

def main():  
    pass

def compute_NSRP(lambda_, beta, mu_c, eta, alpha, theta, year_number=30, 
                 storm=None,DEBUG=False):
    r'''Function to calculate the rainfall for a given set of previously
    calculated NSRP parameters.
    
    Parameters
    ----------
    n_storm : Integer 
        Number of storms
    lambda_ : 
        [1/h]
    beta : 
        [1/h]
    mu_c : 
        per storm
    eta : 
        [1/h]        
    theta : 
        []
        
    Returns
    -------
    rain_h : 
        Precipitation in one hour [mm]
    rain : 
        Precipitation [mm]  
    rain_d : 
        [mm]    
        
    References
    ----------
    .. [1] Weisstein, Eric W. "Gamma Distribution." From MathWorld--A
       Wolfram Web Resource.
       http://mathworld.wolfram.com/GammaDistribution.html
    .. [2] Wikipedia, "Gamma-distribution",
       http://en.wikipedia.org/wiki/Gamma-distribution

    Examples
    --------
    >>>hmm!
    '''
    ## Define common use functions
    zeros = numpy.zeros    
#    array = numpy.array
#    transpose = numpy.transpose
    cumsum = numpy.cumsum
    round_ = numpy.round
    max_ = numpy.max
    floor = numpy.floor
    sum_ = numpy.sum

    ## Define common use functions: Random number generation   
    poissrnd = numpy.random.poisson
    geornd = numpy.random.geometric
    exprnd = numpy.random.exponential
    gamrnd = numpy.random.gamma

    ## number of years needed in model
    if storm is None:
        years = (year_number * 365.25 * 24)
        i = 1    
        flag = True    
        while flag:
            storm_max =  cumsum(poissrnd(1/lambda_, i))[-1]
            if years < storm_max: 
                flag = False
            else:
                i += 1
        storm = int(i * 2)
        if DEBUG: print 'Number of storms calculated=', storm
    
    
    ## List for results
    ## Storms are made of cells and all storms are events that are later
    ## transformed to precipitation (raininfall)
    (rain, rain_h, rain_d, event) = ([], [], [], [])  
    (storm_ori_arr, storm_dur_arr) = (zeros(storm), zeros(storm))

    for s in numpy.arange(storm):
        ## Origin time of the storm in hours and then in minutes
        storm_ori =  poissrnd(1/lambda_) 
        storm_ori_min = round_(60 * storm_ori)        
        
        ## Number of cells in the storm
        ## -1 is to make the distribution give 0 values
        ## default geometric is different in MATLAB and Numpy
        ## see http://en.wikipedia.org/wiki/Geometric_distribution
        cell =  geornd(1/(1+mu_c)) - 1      

        if cell != 0:
            ## Time of cell origin with respect to stom origin in hours
            cell_ori = exprnd(1/beta, cell) 
            ## Time of cell duration in hours
            cell_dur =  exprnd(1/eta, cell) 
            ## Intensity of cell in milimeter/hour
            cell_int = gamrnd(alpha, theta, cell) 

            ## Transform time vector from hours to minutes
            cell_ori_min = round_(60 * cell_ori)
            cell_dur_min = round_(60 * cell_dur)

            ## Calculate storm duration in minutes and create empty i holder            
            storm_dur_min = max_(cell_ori_min + cell_dur_min)
            storm_int = zeros(storm_dur_min)

            if DEBUG:
                print '\nOrigin of storm [h]: ', storm_ori
                print 'Number of cells: ', cell
                print 'Cell duration [h]: ', cell_dur
                print 'Cell delta from origin [h]:' , cell_ori
                print 'Rainfall intensity in cell [mm/h]:' , cell_int        
                print 'Time of storm origin [minutes]: ', storm_ori_min
                print 'Storm duration [minutes]      : ', storm_dur_min

            cell_start = cell_ori_min
            cell_end = cell_ori_min + cell_dur_min
            c_f = (cell_dur * 60) / cell_dur_min # correction factor
            ## Now sum all cells along the same time to obtain total storm i
            for c in range(cell):
                storm_int[cell_start[c]:cell_end[c]] += cell_int[c] * c_f[c]
        else:
            storm_int = []
            storm_dur_min = 0
            
        storm_ori_arr[s]= storm_ori_min
        storm_dur_arr[s] = storm_dur_min
        event.append(storm_int)

    storm_ori_cum = cumsum(storm_ori_arr)
    
    ## Total duration time of synthetic rainfall simulation in minutes
    sim_dur = max_(round_(storm_ori_cum + storm_dur_arr))
    rain = zeros(sim_dur)

    if DEBUG: print 'sim duration in months=', sim_dur/60/24/365*12

    storm_start = round_(storm_ori_cum)
    storm_end = storm_start + round_(storm_dur_arr)

    for s in range(storm):
        if event[s] <> []:
            rain[storm_start[s]:storm_end[s]] += event[s]

    n = sim_dur

    ## Rainfall
    rain = rain # [mm/hour]  in minutes

    ## Rainfall intensity vector in hour 
    fr = 60
    m = floor(n / fr)
    rain_h = rain[:m * fr].reshape(m, fr)
    rain_h = sum_(rain_h, 1) / 60 ## [mm/hr]
    
    ## Rainfall intensity vector in days 
    fr = 60 * 24
    m = floor(n / fr)
    rain_d = rain[:m*fr].reshape(m, fr)
    rain_d = sum_(rain_d, 1) / 60 ## [mm/day] 
        
    #result = [transpose(array([rain_h]))] #,
              #transpose(array([rain])),
              #transpose(array([rain_d]))]
    return rain_h         
    #return result

## Inlcude manual override of parameters

def parameter_NSRP(rain_grouped, ts, aggregation, lag, stats=None, 
                   weights=None, change_factors=None, DEBUG=False, 
                   MY_METHOD=True):
    r'''Estimate rainfall parameters for the NSRP model.
  
    Parameters
    ----------
    rain : array
        A list with numpy arrays with rainfall data aggregated (per month usually)
    ts : 
        time step of the input time series [hours]
    ts_new : 
        list with the aggregation periods.
    lag : 
        lag for the calculation of the autocorrelation coefficient.
    stats : List
        Define the statistics to use in the model. The list should contain 
        numbers from 0 to 5, where [0=Eh, 1=VARh, 2=CVh, 3=Rlh, 4=SKh, 5=FFh]
    weights: 
        a list containing the weight factor applied to each statistic and
        aggregation time. Default is none so a weight of 1 is applied to
        everything
    change_factors :
        a dictionary that holds two items, one is the ready to use change factors
        for the stats  being used
        the second item uis the change factor for the mean only
    
    Returns
    -------
    lambda : 
        [1/h]
    beta :  
        [1/h]
    mu_c :  
        per storm
    eta :  
        [1/h]
    alpha :  
        [-]
    theta :  
        [mm/h]
    '''
    DEBUG = False

    ## Define common use functions
    zeros = numpy.zeros
    array = numpy.array
    ones = numpy.ones
    fmin = scipy.optimize.fmin_l_bfgs_b
#    fmin = optimize.fmin_bfgs
#    fmin = optimize.fmin
#    fmin = optimize.fmin_tnc
    
    ## Empty lists for holding calculated parameters
    w_ = []    
    lambda_ = []
    beta = [] 
    mu_c = [] 
    eta = [] 
    alpha =  [] 
    theta = []
    Qest = []
    
    group_len= len(rain_grouped)
    agg_len = len(aggregation)
    
    ## Empty numpy arrays to hold the statistics and the function evaluation
    Xf = zeros((group_len, 5))
    Q_temp = zeros(group_len)
    
    ## initial conditions   
    ##         lambda   beta    mu_c    eta     alpha
    Xo = array([0.0114, 0.0996, 5.5431, 1.4655, 1.5242])
   
   
    ## Define the statistics to use in the model, this way it can be more inter
    ## actively decided which ones to use
    #  0    1     2    3    4    5
    stats_names = ['Eh', 'VARh', 'CVh', 'Rlh', 'SKh', 'FFh']
    if stats == None:
        stats = [2, 3, 4, 5]  
    len_stats = len(stats)
   
   
    print stats   
    stats_names_used = [stats_names[i] for i in stats]
    if DEBUG: print 'Statistics used: ', stats_names_used
   
    for group in range(group_len):        
        fd_temp = []
        EP = []
        boot_results = []
        
        print 'group: %s' % group 

        for ap in aggregation:           
            (Eh, VARh, CVh, Rlh, SKh, FFh, Fddh, Fwwh) = wg_tool.sample_properties(
              rain_grouped[group], ts, ap, lag, DEBUG=DEBUG, 
              MY_METHOD=MY_METHOD)
              
            stats_list = [Eh, VARh, CVh, Rlh, SKh, FFh]
            stats_used = [stats_list[i] for i in stats]
            
            #fd_temp = fd_temp  + [CVh, Rlh, SKh, FFh]   # Original Value 
            fd_temp = fd_temp  + stats_used
            EP.append(Eh)
            
            
            # IF and only if bootstraping is activated generate some bootstrap
            # errors sampling to get the coefficient of variation FAIL!         
            if weights == 'boot':
                boot_g = wg_tool.aggregate(rain_grouped[group], ts, ap, accumulate=True)                
                boot_verbose = False
                boot_iter = 1500
                boot_dic = wg_tool.bootstrap(boot_g, boot_iter, verbose=boot_verbose)
                boot_results.append(boot_dic['CV'])
                
        fd_ = array(fd_temp)
        EP = array(EP)
        
        # Define change factors
        if change_factors <> None:
            all_stats = change_factors['cf_stats'][group]            
            ep_stats = change_factors['cf_mean'][group]            
#            return {'cf_stats': all_stats, 
#            'cf_mean': ep_stats}            
            fd = fd_ * all_stats
            EP = EP * ep_stats
        else:
            fd = fd_

        # Define weights, one or using bootstrapping coefficient of variation
        if weights == None:
            weights_val = ones(len_stats * agg_len)/(len_stats * agg_len * 1.0)
        elif weights == 'boot':
            temp = [[i]*len_stats for i in boot_results]
            weights_val = []
            for a in temp:
                weights_val += a
            weights_val = 1/numpy.array(weights_val)
            weight_sum = numpy.sum(weights_val)
            weights_val = (weights_val/weight_sum)
            #weights_val = (weights_val/weight_sum)
        else:
            pass
        
        ## [Min, Max]   lambda        beta        mu_c     eta       alpha
        bounds = array([[.0001, .05], [.01, .99], [1, 80], [.5, 30], [.1, 20]])
        #bounds = array([[.0001, None], [.01, None], [1, None], [.5, None], [.1, None]])
        #bounds = array([[None, None], [None, None], [None, None], [None, None], [None, None]])

        ## fmin_bfgs_b   GOOD
        (Xf[group], Q_temp[group], d) = fmin(objective_function, Xo, 
                                             args=(fd, EP, weights_val, 
                                                   aggregation, stats), 
                                             approx_grad=True, 
                                             bounds=bounds, maxfun=3000) 
    
#        # fmin_bfgs  FASTER BUT LESS ACC
#        (Xf[group], Q_temp[group], aa, bb, func_calls, grad_calls, 
#         warnflag) = fmin(objective_function, Xo, None, 
#                   args=(fd, EP, weights, aggregation, stats), full_output=1)      
                   
        ## fmin      BAD
#        (Xf[group], Q_temp[group], iteration, func_calls, 
#         warnflag) = fmin(objective_function, Xo,  
#                   args=(fd, EP, weights, aggregation), full_output=1)     

        ## fmin_tnc VERY BAD
#        (Xf[group], Q_temp[group], rc) = fmin(objective_function, Xo, None, 
#                   args=(fd, EP, weights, aggregation), approx_grad=True,
#                    maxfun=500, bounds= bounds)        
                   
        w_.append(weights_val)
        lambda_.append(Xf[group][0])
        beta.append(Xf[group,1])
        mu_c.append(Xf[group,2])
        eta.append(Xf[group,3])
        alpha.append(Xf[group,4]) 
        theta.append(eta[group] * EP[0] / (alpha[group] * lambda_[group] * 
                     mu_c[group] * 1))
       
        Qest.append(Q_temp[group])

        ## Update initial conditions to find solution faster (hopefully!)
        Xo = Xf[group]
    

    print Q_temp
     
    if DEBUG:
        print 'lan=', ' '.join(['%5.5f' % i for i in lambda_]), ';'
        print 'bet=', ' '.join(['%5.5f' % i for i in beta]), ';'
        print 'muc=', ' '.join(['%5.5f' % i for i in mu_c]), ';'
        print 'eta=', ' '.join(['%5.5f' % i for i in eta]), ';'
        print 'alp=', ' '.join(['%5.5f' % i for i in alpha]), ';'
        print 'tet=', ' '.join(['%5.5f' % i for i in theta]), ';'
        print 'Qest:   ', ' '.join(['%5.5f' % i for i in Qest]), ';'
        print '' 

    dic = {'lambda_': lambda_,
           'beta': beta,
           'mu_c': mu_c,
           'eta': eta,
           'alpha': alpha,
           'theta': theta,
           'Qest': Qest,
           'weights': w_,
           'aggregation': aggregation,
           'statistics': stats_names_used,
           'errorsum': numpy.sum(Qest)}

    return dic
    
    #return [lambda_, beta, mu_c, eta, alpha, theta, Qest]
     
def objective_function(Xo, fd, EP, weights, aggregation, stats):
    r'''
    
    Parameters
    ----------
    Xo: 
        List or NumPy array with the initial values of the NSRP parameters 
    fd: 
        NumPy 1D array containing the statistics use for each aggregation
        period. [E1_ag1, E2_ag1,...E1_agn, E2_agn]
    weights : 
        NumPy 1D array of weights applied to the fd  
    aggregation : 
        List or NumPy array containing the aggregation periods
    
    Returns
    -------

    '''
    ## Define common use functions
    zeros = numpy.zeros
    isnan = numpy.isnan

    ## Initial conditions       
    (lambda_ , beta, mu_c, eta, alpha) = Xo      
    
    
    aggregation = numpy.array(aggregation)
    len_stats = len(stats)
    ## Make uniform, len or shape
    fm = zeros(aggregation.shape[0] * len_stats)

    for i in range(aggregation.shape[0]):
        (Eh, VARh, CVh, Rlh, SKh, FFh, Fddh, Fwwh) = moment_estimation(
          aggregation[i], lambda_, beta, mu_c, eta, alpha, EP[i])
        
        index = len_stats * i 
        
        stats_list = [Eh, VARh, CVh, Rlh, SKh, FFh]
        stats_used = [stats_list[k] for k in stats]

        for j in range(len_stats):
            fm[j + index] = stats_used[j]

    ## WHAT IS THIS DOING?            
    if isnan(numpy.sum(fm)):                        
        fm[isnan(fm)] = 1e+15 
    
    error_value = numpy.sum(weights * ((1 - fd/fm)**2 + (1 - fm/fd)**2)) 
    #print error_value
    return error_value

def moment_estimation(ts, lambda_, beta ,mu_c, eta, alpha , EP): ## Rename EP
    r'''Model Moments Estimation 
    Ref. Cowpertwait 1996 1998 2002  
    
    Parameters
    ----------
    ts : 
        time scale (or time step) interval of properties [h]
    lambda_ : 
        Rate arrival storm [1/h] 
    beta : 
        Rate arrival cell [1/h] 
    eta : 
        cell lifetime [1/h]
    mu_c : 
        per storm [-]
    alpha :
        
    theta :
        Rainfall intensity [mm/h] 
    EP :
        mean precipitation 1 hour aggregation [mm/h]
    
  
    Returns
    -------
    Eh = 
    VARh =
    CVh =
    Rlh =
    SKh =
    FFh =
    Fddh =
    Fwwh 
    '''
    ## Constants
    DEBUG = False
    
    ## Define common use functions
    exp = numpy.exp
    gamma = scipy.special.gamma
    sqrt = numpy.sqrt
    trapz = scipy.integrate.trapz
  
    ## Definition of internal scope functions
    def Gxxh(lag, lambda_, eta, ts, mu_c, EX2, mux, beta, ECC_1):
        r"""Calculate the second moment NSRP
        """
        ## γh,l = COV {Yi_h, Yi+l_h } = λ * η − 3 * A(h, l) * [2 * μ_c * E{X2} + [E{X}]**2 * β**2 * E{C2 − C} / (β**2 − η**2)] − λ * [E{X}]**2 * B(h, l) * E{C2 − C} / [β * (β**2 − η**2)]
        return (lambda_ * (eta**-3) * A(ts, lag, eta) * (2 * mu_c * EX2 + ((mux**2) * (beta**2) * (ECC_1)) / (beta**2 - eta**2)) - (lambda_ * (mux**2) * B(ts, lag, beta) * (ECC_1)) / (beta * (beta**2 - eta**2)))

    def A(ts, lag, eta):
        r"""Coefficient A(h, l) to calculate the second moment of the NSRP

        Source: 
        - Cowpertwait (1998)
        - AWE-GEN-v1.0 Technical Reference - pp 51 - eq. 43
        """
        if lag == 0:
            ## A(h, l) = h * η + exp(−η * h) − 1 , if l = 0
            AA = ts * eta + exp(-eta * ts) - 1
        else:
            ## A(h, l) = 0.5 * (1 − exp(−η * h))**2 * exp(−η * h * (l − 1)), if l > 0
            AA = 0.5 * ((1 - exp(-eta * ts))**2) * exp(-eta * ts * (lag - 1))
        return AA

    def B(ts, lag, beta):
        r"""Coefficient B(h, l) to calculate the second moment of the NSRP

        Source:     
        - Cowpertwait (1998)
        - AWE-GEN-v1.0 Technical Reference - pp 51 - eq. 44
        """
        if lag == 0:
            ## B(h, l) = h * β + exp(−β * h) − 1 , if l = 0
            BB = ts * beta + exp(-beta * ts) - 1
        else:
            ## B(h, l) = 0.5 * (1 − exp(−β * h))**2 * exp(−β * h * (l − 1)) , if l > 0
            BB = 0.5 * ((1 - exp(-beta * ts))**2) * exp(-beta * ts * (lag - 1))
        return BB

    def f_1(eta, beta, ts):
        r"""Coefficient f(η, β, h) to calculate the third moment of the NSRP
    
        Source:     
        - Cowpertwait (1998) 
        - AWE-GEN-v1.0 Technical Reference - pp 51 - eq. 46
        """
        ## f(η, β, h) = − 2 * η**3 * β**2 * exp(−η * h) − 2 * η**3 * β**2 * exp(−β * h)        + η**2 * β**3 * exp(−2 * η * h) + 2 * η**4 * β * exp(−η * h)                + 2 * η**4 * β * exp(−β * h) + 2 * η**3 * β**2 * exp(−(η + β) * h)                     − 2 * η**4 * β * exp(−(η + β) * h) − 8 * η**3 * β**3 * h)                 + 11 * η**2 * β**3 − 2 * η**4 * β + 2 * η**3 * β**2                + 4 * η * β**5 * h + 4 * η**5 * β* h               − 7 * β**5 − 4 * η**5 + 8 * β**5 exp(−η * h) − β**5 * exp(−2 * η * h)                    − 2 * h * η**3 * β**3 * exp(−η * h) − 12 * η**2 * β**3 * exp(−η * h)                  + 2 * h * η * β**5 * exp(−η * h) + 4 * η**5 * exp(−β * h)
        return (-2 * eta**3 * beta**2 * exp(-eta * ts) -2 * eta**3 * beta**2 * exp(-beta * ts) + eta**2 * beta**3 * exp(-2 * eta * ts)+ 2 * eta**4 * beta * exp(-eta * ts) + 2 * eta**4 * beta * exp(-beta * ts) + 2 * eta**3 * beta**2 * exp(-(eta + beta) * ts) - 2 * eta**4 * beta * exp(-(eta + beta) * ts) - 8 * eta**3 * beta**3 * ts + 11 * eta**2 * beta**3 - 2 * eta**4 * beta + 2 * eta**3 * beta**2 + 4 * eta * beta**5 * ts  + 4 * eta**5 * beta * ts - 7 * beta**5 - 4 * eta**5 + 8 * beta**5 * exp(-eta * ts) - beta**5 * exp(-2 * eta * ts) - 2 * ts * eta**3 * beta**3 * exp(-eta * ts) - 12 * eta**2 * beta**3 * exp(-eta * ts) + 2 * ts * eta * beta**5 * exp(-eta * ts) + 4 * eta**5 * exp(-beta * ts))

    def g_1(eta, beta, ts):
        r"""Coefficient g(η, β, h) to calculate the third moment of the NSRP
    
        Source:     
        - Cowpertwait (1998)
        - AWE-GEN-v1.0 Technical Reference - pp 51 - eq. 47
        """
        ## g(η, β, h) = 12 * η**5 * β exp(−β * h) + 9 * η**4 * β**2         + 12 * η * β**5 * exp(−η * h)         + 9 * η**2 * β**4 + 12 * η**3 * β**3 * exp(−(η + β) * h)                 − η**2 * β**4 * exp(−2 * η * h) − 12 * η**3 * β**3 * exp(−β * h)                  − 9 * η**5 * β − 9 * η * β**5 − 3 * η * β**5 * exp(−2 * η * h)                  − η**4 * β**2 * exp(−2 * β * h) − 12* η**3 * β**3 * exp(−η * h)                   + 6 * η**5 * β**2 * h − 10 * η**3 * β**4 * h            + 6 * η**2 * β**5 * h −10 * η**4 * β**3 * h + 4 * η * β**6 * h                   − 8 * β**2 * η**4 * exp(−β * h) + 4 * β * η**6 * h                + 12 * β**3 * η**3  − 8 * β**4 * η**2 * exp(−η * h) − 6 * η**6               − 6 * β**6 − 2 * η**6 exp(−2 * β* h) − 2 * β**6 * exp(−2 * η * h)                   + 8 * η**6 * exp(−β * h) + 8 * β**6 * exp(−η * h) − 3 * β * η**5 * exp(−2 * β * h)
        return (12 * eta**5 * beta * exp(-beta * ts) + 9 * eta**4 * beta**2 + 12 * eta * beta**5 * exp(-eta * ts) + 9 * eta**2 * beta**4 + 12 * eta**3 * beta**3 * exp(-(eta + beta) * ts) - eta**2 * beta**4 * exp(-2 * eta * ts) - 12 * eta**3 * beta**3 * exp(-beta * ts) -9 * eta**5 * beta - 9 * eta * beta**5 - 3 * eta * beta**5 * exp(-2 * eta * ts) - eta**4 * beta**2 * exp(-2 * beta * ts) - 12 * eta**3 * beta**3 * exp(-eta * ts) + 6 * eta**5 * beta**2 * ts -10 * eta**3 * beta**4 * ts + 6 * eta**2 * beta**5 * ts -10 * eta**4 * beta**3 * ts + 4 * eta * beta**6 * ts - 8 * beta**2 * eta**4 * exp(-beta * ts) + 4 * beta * eta**6 * ts + 12 * beta**3 * eta**3 - 8 * beta**4 * eta**2 * exp(-eta * ts) - 6 * eta**6 - 6 * beta**6 - 2 * eta**6 * exp(-2 * beta * ts) - 2 * beta**6 * exp(-2 * eta * ts) + 8 * eta**6 * exp(-beta * ts) + 8 * beta**6 * exp(-eta * ts)-3 * beta * eta**5 * exp(-2 * beta * ts))
     
    def intp(ts, beta, eta, mu_c):
        r"""
        """
        t = numpy.arange(0, 90.01, 0.01)
        ## ph(t) = exp(−β * (t + h)) + 1 − (η * exp(−β * t) − β * exp(−η * t) / (η − β) ] * exp[ −μ_c * β * (exp(−β * t) − exp(−η * t)) / (η − β) − μ_c * exp(−β*t) + μ_c * exp(−β * (t + h)) ]
        pht = ((exp(-beta * (t + ts)) + 1 -(eta * exp(-beta * t) - beta * exp(-eta * t)) / (eta - beta)) * (exp((-mu_c * beta * (exp(-beta * t) - exp(-eta * t)) / (eta - beta)) - mu_c * (exp(-beta * t)) + mu_c * exp(-beta * (t + ts)))))
        return trapz((1 - pht),t)   

    ## Main body of the function
    ## AWE-GEN v1.0 Technical Reference - pp. 7 - eq. 7
    theta = eta * EP / (alpha * lambda_ * mu_c * ts) ## for gamma distribution 
    
    ## WEIBULL 
    # Bw = 1 / alpha; - shape parameter weibull
    # Aw = theta**-Bw; - scale parameter  weibull
    # mux =  Aw * gamma(1 + Bw**-1)      mux = (theta**alp) * (gamma(1 + alpha))
    # EX2 = (Aw**2) * gamma(1 + 2 * Bw**-1) EX2 = (theta**(2 * alpha)) * gamma(1 + 2 * alp)
    # EX3 = (Aw**3) * gamma(1 + 3 * Bw**-1) EX3 = (theta**(3 * alpha)) * gamma(1 + 3 * alp)
  
    ## GAMMA / AWE-GEN Technical Reference pp. 7  
    ## E{Xn} = theta**n * gamma(alpha + n) / gamma(alpha) 
    ## E{X} = μ_x = mu_x
    mu_x = theta * alpha
    EX2 = (theta ** 2) * gamma(alpha + 2) / gamma(alpha)
    EX3 = (theta ** 3) * gamma(alpha + 3) / gamma(alpha)
    
    ## EXPONENTIAL 
    # mu_x = theta 
    # EX2 = theta**2 # theta / 2 
    # EX3 = 2
    # POISSON C  
    # ECC_1 = mu_c**2
    # ECC_12 = mu_c**3
    # POISSON C-1
    # mu_c = mu_c
    # ECC_1 = mu_c**2 - 1 
    
    ## GEOMETRIC - AWE-GEN v1.0 Technical Reference - pp. 7
    ## E{C2 −C} = 2 * μ_c * (μ_c − 1)
    ECC_1 = 2 * mu_c * (mu_c - 1) 
    ## E{(C2 − C)(C − 2)} = 6 * μ_c * (μ_c − 1)**2
    ECC_12 = 6 * mu_c * (mu_c - 1)**2
    ## AWE-GEN Technical Reference Equation (45)
    ## ξh = E{[Yh − E{Yh}]3} = 6 * λ * μc * E{X3} * (η * h − 2 + η * h * exp(−η * h) + 2 * exp(−η * h)) / η**4 + 3 * λ * E{X} * E{X2} * E{C(C − 1)} * f(η, β, h) / [2 * η * 4 * β * (β**2 − η**2)**2] + λ * E{X}**3 * E{(C2 − C)(C − 2)} * g(η, β, h) / [2 * η**4 * β * (η**2 − β**2) * (η − β) * (2 * β + η) * (β + 2 * η)]
    M3h = (6 * lambda_ * mu_c * EX3 * (eta * ts -2 + eta * ts * exp(-eta * ts) + 2 * exp(-eta * ts)) / (eta**4) + 3 * lambda_ * mu_x * EX2 * (ECC_1) * f_1(eta, beta, ts) / (2 * eta**4 * beta * ((beta**2 - eta**2)**2)) + lambda_ * (mu_x**3) * (ECC_12) * g_1(eta, beta, ts) / (2 * eta**4 * beta * (eta**2 - beta**2) * (eta - beta) * (2 * beta + eta) * (beta + 2 * eta)))
   
    ## μh = E{Yh} = λ * μ_c * E{X} * h / η ,
    Eh = lambda_ * mu_c * mu_x * ts / eta
    VARh = Gxxh(0, lambda_, eta, ts, mu_c, EX2, mu_x, beta, ECC_1)
    CVh = sqrt(VARh) / Eh
    Rh =  Gxxh(1, lambda_, eta, ts, mu_c, EX2, mu_x, beta, ECC_1) / VARh
    SKh = M3h / ((sqrt(VARh))**3)
  
    ## Probability of a dry interval
    ## AWE-GEN v1.0 Technical Reference - pp. 52 - eq. 48
    ## Φ(h) = exp(−λ h + λ β−1μ−1c [1 − e(−μc+μce􀀀
h)]
    FFh = (exp(-lambda_ * ts + (lambda_ * (beta**-1)) * (mu_c**-1) * (1 - exp(-mu_c + mu_c * exp(-beta * ts))) - lambda_ * intp(ts, beta, eta, mu_c)))
    FF2h = (exp(-lambda_ * (2 * ts) + (lambda_ * (beta**-1)) * ((mu_c)**-1) * (1 - exp(-mu_c + mu_c * exp(-beta * (2 * ts)))) - lambda_ * intp((2 * ts), beta, eta, mu_c)))
    Fddh = FF2h / FFh 
    Fwwh = (1 - 2 * FFh + FF2h) / (1 - FFh)
    
    if DEBUG:
        print 'theta:', theta
        print 'mu_x:', mu_x 
        print 'EX2:', EX2
        print 'EX3:', EX3
        print 'ECC_1:', ECC_1
        print 'ECC_12:', ECC_12
        print 'M3h:', M3h
        print 'Eh:', Eh
        print 'VARh:', VARh
        print 'CVh:', CVh
        print 'Rh:', Rh
        print 'SKh:', SKh
        print 'FFh:', FFh
        print 'FF2h:', FF2h
        print 'Fddh:', Fddh
        print 'Fwwh:', Fwwh
        print ''
        
    return [Eh, VARh, CVh, Rh, SKh, FFh, Fddh, Fwwh]

def run_para(ifile, ofile, ts, aggregation, lag, group_type='m', stats=None, 
             weights=None, MY_METHOD=True, idir_cf=None, aggregation_cf=None):
    r'''
    '''
    if idir_cf == None:
        change_factors = None
    else:
        # Calculate change factor
        #aggregation_cf = [24, 48, 72, 96]
        change_factors = wg_downscale.change_factors(idir_cf, ifile, aggregation, aggregation_cf)
        
    # load data
    dic = wg_io.load(ifile)
    rain_data = dic['rainfall']

    # Make the grouping for the parameter
    if group_type == 'm':
        group_month = dic['month']        
        rain_grouped = wg_tool.group_by(rain_data, group_month)
        #             Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec 
        day_normal = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        day_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        
    elif group_type == 'y':
        rain_grouped = [rain_data]
        day_normal = [365]
        day_leap = [366]
    else:
        group_tail = True
        doy = dic['doy']
        group_doy = wg_tool.window_doy(doy, group_type, group_tail=group_tail)    
        rain_grouped = wg_tool.group_by(rain_data, group_doy)
        
        tail = (365.0/group_type - numpy.floor(365.0/group_type))*group_type
        
        day_normal = [group_type] * numpy.floor(365.0/group_type) 
        day_leap = [group_type] * numpy.floor(365.0/group_type) 
        
        if  group_tail:
            day_normal[-1] = numpy.int(numpy.round(group_type + tail,0))
            day_leap[-1] = numpy.int(numpy.round(group_type + tail + 1,0))
        else:
            day_normal.append(numpy.int(numpy.round(tail, 0)))
            day_leap.append(numpy.int(numpy.round(tail + 1, 0)))
        
        #10 + [(365-14*10)]
    
    dic = parameter_NSRP(rain_grouped, ts, aggregation, lag,  stats=stats, 
                         weights=weights, change_factors=change_factors, DEBUG=True, 
                         MY_METHOD=MY_METHOD)

    dic['day_normal'] = day_normal
    dic['day_leap'] = day_leap
                                 
    # Here is to save
    print 'Saving: ' + str(ofile)
    wg_io.save(dic, ofile)

    return dic

def run_gen(ifile, ofile, start_year, end_year):
    r'''
    
    Parameters
    ----------
    ifile : String
        TODO
    save_as : String
        TODO
    start_year : Integer
        Inclusive
    end_year : Integer
        Inclusive

    Returns
    -------
    
    '''
    dic = wg_io.load(ifile)
    day_normal = dic['day_normal']
    day_leap = dic['day_leap']
    list_lambda_ = dic['lambda_']
    list_beta = dic['beta']
    list_mu_c = dic['mu_c']
    list_eta = dic['eta']
    list_alpha = dic['alpha'] 
    list_theta = dic['theta']

    size = len(day_normal)
    period = (end_year - start_year)

    day_location = [0]    

    # Create a vector containing the exact amount of hourly time steps needed
    # to fill the period 
    for i in range(period):
        if calendar.isleap(start_year + i): 
            day_location = day_location + day_leap
        else:
            day_location = day_location + day_normal

    hour_index = numpy.array(day_location) * 24    
    hour_index_cum = numpy.cumsum(hour_index)
    data_size = numpy.sum(hour_index)
    hour_grouped = [ [hour_index_cum[i], hour_index_cum[i+1]] for i in range(size*period)]

    data = numpy.zeros(data_size)
            
    for i in range(size):
        day_size = day_leap[i]
        year_size = numpy.int(numpy.round((day_size/365.25) * period, 0))

        lambda_ = list_lambda_[i]
        beta = list_beta[i]
        mu_c = list_mu_c[i]
        eta = list_eta[i]
        alpha = list_alpha[i]
        theta = list_theta[i]

        temp_data = compute_NSRP(lambda_, beta, mu_c, eta, alpha, theta, year_number=year_size, storm=None, DEBUG=False)

        t_i_start = 0
        for p in range(period):            
            list_index = (size*p + i)
            i_start = hour_grouped[list_index][0]
            i_end = hour_grouped[list_index][1]
            delta = i_end - i_start
                        
            t_i_end = t_i_start + delta             
            data[i_start:i_end] = temp_data[t_i_start:t_i_end]
            t_i_start =  t_i_end
            
    # Create time vectors 
    dtype_int = numpy.int  
    
    ## Data holders for date, variable and the broken date
    date = numpy.zeros(data_size)
    year = numpy.zeros(data_size, dtype=dtype_int)
    month = numpy.zeros(data_size, dtype=dtype_int)
    day = numpy.zeros(data_size, dtype=dtype_int)
    hour = numpy.zeros(data_size, dtype=dtype_int)
    minute = numpy.zeros(data_size, dtype=dtype_int)
    doy = numpy.zeros(data_size, dtype=dtype_int)


    date_format='%Y-%m-%d %H:%M'
    date_start_string = str(start_year) + '-01-01 00:00'
    date_start = datetime.datetime.strptime(date_start_string, date_format )


    for t in range(data_size):
        h = datetime.timedelta(hours=t)
        datevalue = date_start + h
        date[t] = numpy.float(matplotlib.dates.date2num(datevalue))
        year[t] = datevalue.year
        month[t] = datevalue.month
        day[t] = datevalue.day
        hour[t] = datevalue.hour
        minute[t] = datevalue.minute
        doy[t] = datevalue.timetuple().tm_yday
    
    comment = ''
    dic = {'date' : date, 
           'rainfall': data,
           'year' : year,
           'month' : month,
           'day' : day,           
           'hour' : hour,           
           'minute' : minute,
           'doy' : doy,
           'comment': comment}
    
    print '\nSaving as: ' + ofile
    wg_io.save(dic, ofile)
  
if __name__ == '__main__':
    main()\00
