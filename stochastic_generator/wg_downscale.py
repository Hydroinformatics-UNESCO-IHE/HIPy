#!/usr/bin/env python
# -*- coding: utf-8 -*-
r'''
'''
import os

import numpy
import scipy.stats.mstats
import scipy.optimize
import matplotlib.pyplot

import wg_io
import wg_tool

def main():
    pass
    idir= r'D:\Dropbox\module_14\data\03_processed\gcm_mcmc\kochi\2081-2100\sresa2'
    ifile = r'D:\Dropbox\module_14\data\02_preprocessed\02_pickled\2000_kochi.data'
    extend_aggregation = [1, 2, 4, 6, 24, 48, 72, 96]
    extend_aggregation = [1, 6, 24, 72]
    ofile = 'test.cf'
    aggregation = [24, 48, 72, 96]
    change_factors(idir, ifile, extend_aggregation, aggregation)




def change_factors(idir, ifile, extend_aggregation, aggregation, stats=None):
    r'''
    '''    
    dic = calculate_change_factor(idir, ifile)
    fmin = scipy.optimize.fmin_l_bfgs_b
    fmin = scipy.optimize.fmin        

    # Changed factors original data
    dsf = dic['dsf']
    mean = dic['mean']
    var = dic['var']
    skew = dic['skew']

    cv_ori = []
    correl_ori = []
    skew_ori = []
    mean_ori = []
    var_ori = []
    dsf_ori = []

    # Find the original 
    # -----------------
    dic = wg_io.load(ifile)
    ts = 1
    lag = 1
    accumulate = True
    group_type = 'm'
    dic_stat = wg_tool.series_statistics(dic, ts, extend_aggregation, lag=lag, 
                                  group_type=group_type, accumulate=accumulate)

    # [Ey, VARy, CVy, SKy, Rly, DSFy]
    dicdic = dic_stat['period']
    for g in range(12):
        mean_index = 0
        mean_ori.append(dicdic[mean_index][g])

        var_index = 1
        var_ori.append(dicdic[var_index][g])

        cv_index = 2
        cv_ori.append(dicdic[cv_index][g])

        skew_index = 3
        skew_ori.append(dicdic[skew_index][g])

        correl_index = 4
        correl_ori.append(dicdic[correl_index][g])        

        dsf_index = 5
        dsf_ori.append(dicdic[dsf_index][g])        



    # Find indexes to use extension values for bellow 24 hours and calculates 
    # values for over 24 hours    
    n_group = len(dsf)    
    pre_index = []
    post_index = []    

    for i in range(len(extend_aggregation)):
        if not(extend_aggregation[i] in aggregation):
            pre_index.append(i)
    
    for i in range(len(aggregation)):    
        if (aggregation[i] in extend_aggregation):
            post_index.append(i)

    h = numpy.array(aggregation)
    h_ext = numpy.array(extend_aggregation)
        
    # Extend mean: E(h) = A*h. Is linear and cuts in 0. Either one value of future
    # is used, or a regression is made...
    # ---------------------------------------------------------------------- 
    mean_new = []
    for g in range(n_group):
        mean_h = mean[g]        
        Xo = [1]   # A
        Xf = fmin(mean_min, Xo, args=(h, mean_h), disp=0)
        A = Xf[0]
        mean_adj = A * h_ext
        
        mean_res = [mean_adj[i] for i in pre_index] + [mean_h[j] for j in post_index]
        mean_new.append(numpy.array(mean_res))

#        matplotlib.pyplot.plot(h_ext, mean_adj, 'o--')
#        matplotlib.pyplot.plot(h, mean_h)
#        matplotlib.pyplot.show()
        
    # Extend dry spell fraction: dsf(h) = A*(e^h*B)
    # ---------------------------------------------
    # TODO: include the different groupings
    dsf_new = []
    for g in range(n_group):
        dsf_h = dsf[g]        
        Xo = [-1]   # A
        
        if False:
            # Only use 24 72        
            hh = [h[0], h[1]]
            dsf_hh = [dsf_h[0], dsf_h[1]]
        else:
            hh = h
            dsf_hh = dsf_h

        #print fmin(dsf_min, Xo, args=(h, dsf_h), approx_grad=True, bounds=bounds, maxfun=3000)     
        Xf = fmin(dsf_min, Xo, args=(hh, dsf_hh), disp=0)
        A = Xf[0]

        
        dsf_adj = numpy.exp(A * h_ext)
        
        dsf_res = [dsf_adj[i] for i in pre_index] + [dsf_h[j] for j in post_index]
        dsf_new.append(numpy.array(dsf_res))

#        matplotlib.pyplot.plot(h_ext, dsf_adj, 'o--')
#        matplotlib.pyplot.plot(h, dsf_h)
#        matplotlib.pyplot.show()

    # Skewness.... hmmmm difficult so far so only observed is downscaled
    # ---------------------------------------------
    # TODO: include the different groupings
    skew_new = []
    for g in range(n_group):
        skew_h = skew[g]
        skew_adj = skew_ori[g]
        skew_res = [skew_adj[i] for i in pre_index] + [skew_h[j] for j in post_index]
        skew_new.append(numpy.array(skew_res))

    # Autocorrelation.... very difficult indeed not taken into account
    # ---------------------------------------------    
    correl_new = correl_ori

    # Variance
    # -------------------------------------------------
    var_new = []
    for g in range(n_group):
        var_h = var[g]

        Xo = [0.5, 1, 0.8]
        Xf = fmin(var_min, Xo, args=(aggregation, var_h), disp=0)

        eps =Xf[0]     ## [h]
        Sig2i = Xf[1]  ## [mm^2]
        alp = Xf[2]    ## []

        var_adj = variance_t_downscaling(extend_aggregation, eps, Sig2i, alp)
        var_res = [var_adj[i] for i in pre_index] + [var_h[j] for j in post_index]
        var_new.append(numpy.array(var_res))

    # Coefficient of variance
    #-----------------------
    cv_new = [(var_new[i])**0.5 / mean_new[i] for i in range(len(mean_new))]
    
    dic = {'aggregation': extend_aggregation,
           'dsf': dsf_new,
           'mean': mean_new,
           'var': var_new,
           'skew': skew_new,
           'correl': correl_new,
           'cv': cv_new}
           
    mean_cf = numpy.array(mean_new) / numpy.array(mean_ori)
    var_cf = numpy.array(var_new) / numpy.array(var_ori)
    cv_cf = numpy.array(cv_new) / numpy.array(cv_ori)
    skew_cf = numpy.array(skew_new) / numpy.array(skew_ori)
    correl_cf = numpy.array(correl_new) / numpy.array(correl_ori)
    dsf_cf = numpy.array(dsf_new) / numpy.array(dsf_ori)

    ## Define the statistics to use in the model, this way it can be more inter
    ## actively decided which ones to use
    #  0    1     2    3    4    5
    stats_names = ['Eh', 'VARh', 'CVh', 'Rlh', 'SKh', 'FFh']
    if stats == None:
        stats = [2, 3, 4, 5]  

    len_stats = len(stats)
  
    stats_names_used = [stats_names[i] for i in stats]
    print stats_names_used


    all_stats = []
    ep_stats = []
    # Now arrange so it is easy to use in the parameter generator
    for g in range(n_group):
        temp_stats = []
        temp_ep = []
        for ag in range(len(extend_aggregation)):
            stats_list = [mean_cf[g][ag], var_cf[g][ag], cv_cf[g][ag], 
                         correl_cf[g][ag], skew_cf[g][ag], dsf_cf[g][ag]]
                    
            stats_used = [stats_list[i] for i in stats]                    
                    
            temp_stats = temp_stats + stats_used
            temp_ep = temp_ep + [mean_cf[g][ag]]
        
        all_stats.append(numpy.array(temp_stats))
        ep_stats.append(numpy.array(temp_ep))
        #stats_list.append(Eh, VARh, CVh, Rlh, SKh, FFh]

    print all_stats
    print ep_stats
    return {'cf_stats': all_stats, 
            'cf_mean': ep_stats}
           
    # Plot test
#    for g in range(n_group):
#        for key in dic:
#            x = dic['aggregation']
#            if not(key == 'aggregation'):
#                y = dic[key][g]
#                
#                fig = matplotlib.pyplot.figure()
##                fig.subplots_adjust(left=0.10, bottom=0.10, wspace=0.40, hspace=0.20) 
#                fig.suptitle(key)
#                ax = fig.add_subplot(111)
#                ax.set_axisbelow(True)
#                ax.plot(x, y, 'b-', label='period statistics')
#                matplotlib.pyplot.show()

    #wg_io.save(dic, ofile)





def extend_fine(idir, ifile, extend_aggregation, aggregation, ofile):
    r'''
    '''    
    dic = calculate_change_factor(idir, ifile)
    fmin = scipy.optimize.fmin_l_bfgs_b
    fmin = scipy.optimize.fmin        

    # Changed factors original data
    dsf = dic['dsf']
    mean = dic['mean']
    var = dic['var']
    skew = dic['skew']
    correl_ori = []
    skew_ori = []

    # Find the original 
    # -----------------
    dic = wg_io.load(ifile)
    ts = 1
    lag = 1
    accumulate = True
    group_type = 'm'
    dic_stat = wg_tool.series_statistics(dic, ts, extend_aggregation, lag=lag, 
                                  group_type=group_type, accumulate=accumulate)

    # [Ey, VARy, CVy, SKy, Rly, DSFy]
    dicdic = dic_stat['period']
    for g in range(12):
        skew_index = 3
        skew_ori.append(dicdic[skew_index][g])

        correl_index = 4
        correl_ori.append(dicdic[correl_index][g])        

    # Find indexes to use extension values for bellow 24 hours and calculates 
    # values for over 24 hours    
    n_group = len(dsf)    
    pre_index = []
    post_index = []    

    for i in range(len(extend_aggregation)):
        if not(extend_aggregation[i] in aggregation):
            pre_index.append(i)
    
    for i in range(len(aggregation)):    
        if (aggregation[i] in extend_aggregation):
            post_index.append(i)

    h = numpy.array(aggregation)
    h_ext = numpy.array(extend_aggregation)
        
    # Extend mean: E(h) = A*h. Is linear and cuts in 0. Either one value of future
    # is used, or a regression is made...
    # ---------------------------------------------------------------------- 
    mean_new = []
    for g in range(n_group):
        mean_h = mean[g]        
        Xo = [1]   # A
        Xf = fmin(mean_min, Xo, args=(h, mean_h), disp=0)
        A = Xf[0]
        mean_adj = A * h_ext
        
        mean_res = [mean_adj[i] for i in pre_index] + [mean_h[j] for j in post_index]
        mean_new.append(numpy.array(mean_res))

#        matplotlib.pyplot.plot(h_ext, mean_adj, 'o--')
#        matplotlib.pyplot.plot(h, mean_h)
#        matplotlib.pyplot.show()
        
    # Extend dry spell fraction: dsf(h) = A*(e^h*B)
    # ---------------------------------------------
    # TODO: include the different groupings
    dsf_new = []
    for g in range(n_group):
        dsf_h = dsf[g]        
        Xo = [-1]   # A
        
        if False:
            # Only use 24 72        
            hh = [h[0], h[1]]
            dsf_hh = [dsf_h[0], dsf_h[1]]
        else:
            hh = h
            dsf_hh = dsf_h

        #print fmin(dsf_min, Xo, args=(h, dsf_h), approx_grad=True, bounds=bounds, maxfun=3000)     
        Xf = fmin(dsf_min, Xo, args=(hh, dsf_hh), disp=0)
        A = Xf[0]

        
        dsf_adj = numpy.exp(A * h_ext)
        
        dsf_res = [dsf_adj[i] for i in pre_index] + [dsf_h[j] for j in post_index]
        dsf_new.append(numpy.array(dsf_res))

#        matplotlib.pyplot.plot(h_ext, dsf_adj, 'o--')
#        matplotlib.pyplot.plot(h, dsf_h)
#        matplotlib.pyplot.show()

    # Skewness.... hmmmm difficult so far so only observed is downscaled
    # ---------------------------------------------
    # TODO: include the different groupings
    skew_new = []
    for g in range(n_group):
        skew_h = skew[g]
        skew_adj = skew_ori[g]
        skew_res = [skew_adj[i] for i in pre_index] + [skew_h[j] for j in post_index]
        skew_new.append(numpy.array(skew_res))

    # Autocorrelation.... very difficult indeed not taken into account
    # ---------------------------------------------    
    correl_new = correl_ori

    # Variance
    # -------------------------------------------------
    var_new = []
    for g in range(n_group):
        var_h = var[g]

        Xo = [0.5, 1, 0.8]
        Xf = fmin(var_min, Xo, args=(aggregation, var_h), disp=0)

        eps =Xf[0]     ## [h]
        Sig2i = Xf[1]  ## [mm^2]
        alp = Xf[2]    ## []

        var_adj = variance_t_downscaling(extend_aggregation, eps, Sig2i, alp)
        var_res = [var_adj[i] for i in pre_index] + [var_h[j] for j in post_index]
        var_new.append(numpy.array(var_res))

    # Coefficient of variance
    #-----------------------
    cv_new = [(var_new[i])**0.5 / mean_new[i] for i in range(len(mean_new))]
    
    dic = {'aggregation': extend_aggregation,
           'dsf': dsf_new,
           'mean': mean_new,
           'var': var_new,
           'skew': skew_new,
           'correl': correl_new,
           'cv': cv_new}
           
    # Plot test
#    for g in range(n_group):
#        for key in dic:
#            x = dic['aggregation']
#            if not(key == 'aggregation'):
#                y = dic[key][g]
#                
#                fig = matplotlib.pyplot.figure()
##                fig.subplots_adjust(left=0.10, bottom=0.10, wspace=0.40, hspace=0.20) 
#                fig.suptitle(key)
#                ax = fig.add_subplot(111)
#                ax.set_axisbelow(True)
#                ax.plot(x, y, 'b-', label='period statistics')
#                matplotlib.pyplot.show()

    wg_io.save(dic, ofile)

    return dic

def variance_t_downscaling(T, eps, Sig2i, alp):
    Sig2 = numpy.zeros(len(T)) 
    exp = numpy.exp
    
    for j in range(len(T)):
        if T[j] <= eps:
            Sig2[j] = 2*Sig2i*(eps/alp)*((eps/alp)*(exp(-alp*T[j]/eps)-1)+T[j])
        else:
            Sig2[j] = 2*Sig2i*( ((eps**alp)*exp(-alp)/((1-alp)*(2-alp)))*(T[j]**(2-alp)) + (eps/alp)*(1- (exp(-alp))/(1-alp))*T[j] + ((eps/alp)**2)*(exp(-alp)-1) + 2*((eps**2)*exp(-alp))/(alp*(2-alp))  )

    return Sig2

def scale_var(VARh, T):
    fmin = scipy.optimize.fmin
    Sig2_obs = VARh    
    Xo = [0.5, 1, 0.8]
    
    Xf = fmin(var_min, Xo, args=(T, Sig2_obs), disp=0)
    
    eps =Xf[0]     ## [h]
    Sig2i = Xf[1]  ## [mm^2]
    alp = Xf[2]    ## []
    
    return [eps, Sig2i, alp]

def var_min(Xo, T, Sig2_obs):
    r'''
    '''
    exp = numpy.exp   
    
    Sig2 = numpy.ones(len(T))
    
    eps = Xo[0]
    Sig2i = Xo[1]
    alp = Xo[2]

    for j in range(len(T)):
        if T[j] <= eps:
            Sig2[j] = 2*Sig2i*(eps/alp)*((eps/alp)*(exp(-alp*T[j]/eps)-1) + T[j])
        else:
            Sig2[j] = 2*Sig2i*( ((eps**alp)*exp(-alp)/((1-alp)*(2-alp)))*(T[j]**(2-alp)) + (eps/alp)*(1- (exp(-alp))/(1-alp))*T[j] + ((eps/alp)**2)*(exp(-alp)-1) + 2*((eps**2)*exp(-alp))/(alp*(2-alp))  )            

    if numpy.isnan(numpy.sum(Sig2)):
        Sig2[numpy.isnan(Sig2)] = 1e+15
    
    error = numpy.sum(((1 - Sig2/Sig2_obs)**2 + (1 - Sig2_obs/Sig2)**2))
#    print error
    return error
        
def dsf_min(Xo, h, dsf):
    r'''
    '''
    # equation dsf(h) = (e^h*A)
    h = numpy.array(h)
    dsf = numpy.array(dsf)
    
    A = Xo[0]
  
    error = numpy.sum(((numpy.exp(A * h))/dsf - 1)**2 + (dsf/(numpy.exp(A * h)) - 1)**2)
    return error
    
def mean_min(Xo, h, mean):
    r'''
    '''
    # equation mean(h) = A*h
    h = numpy.array(h)
    mean = numpy.array(mean)
    
    A = Xo[0]
    
    #error = numpy.sum(((A * h) - mean)**2)
    error = numpy.sum(((A * h)/mean -1)**2 + (mean/(A * h) -1)**2)    
    return error
    
def calculate_change_factor(idir, ifile, prob=[0.10, 0.90], aggregation=[24,48,72,96]):
    r'''
    '''
    ifiles = sorted(os.listdir(idir))
    dic = {}
    dic_dsf= {}
            
    for i in ifiles:
        
        fname = os.path.join(idir, i)
        d = wg_io.load(fname)        
        mu = numpy.array(d[0])
        nu = numpy.array(d[1])
        
        mu_dsf = mu/(1 - mu)
        nu_dsf = nu/(1 - nu)

#        print 'CF: ', numpy.mean(nu/mu)
#        print '\n'

        prob_up = prob[1]
        prob_down = prob[0]
        axis = 0    
        alpha = 0.0
        beta = 1.0

        up = scipy.stats.mstats.mquantiles(nu/mu, prob=prob_up, alphap=alpha, betap=beta, axis=axis)
        down = scipy.stats.mstats.mquantiles(nu/mu, prob=prob_down, alphap=alpha, betap=beta, axis=axis)        
        dic[i] = [numpy.mean(nu/mu), down, up]
        
        up_dsf = scipy.stats.mstats.mquantiles(nu_dsf/mu_dsf, prob=prob_up, alphap=alpha, betap=beta, axis=axis)
        down_dsf = scipy.stats.mstats.mquantiles(nu_dsf/mu_dsf, prob=prob_down, alphap=alpha, betap=beta, axis=axis)        
        dic_dsf[i] = [numpy.mean(nu_dsf/mu_dsf), down_dsf, up_dsf]

#    aggregation = [24]
    
    # Preprocess
    #size_in = size_cm[0]/2.54, size_cm[1]/2.54       

    data_mean = {}
    data_var = {}
    data_skew = {}
    data_dsf = {}
    
    # TODO: include the other type of groupings!
    for g in range(12):
        temp_mean = []
        temp_var = []
        temp_skew = []
        temp_dsf = []

        for k in dic:
            for agg in aggregation:
                if ('agg_' + str(agg) in k) and ('stat_Eh' in k) and ('g_' + str(g) + '_' in k):                
                    temp_mean.append([k, dic[k][0], dic[k][1], dic[k][2]])

                if ('agg_' + str(agg) in k) and ('stat_VARh' in k) and ('g_' + str(g) + '_' in k):                
                    temp_var.append([k, dic[k][0], dic[k][1], dic[k][2]])

                if ('agg_' + str(agg) in k) and ('stat_SKh' in k) and ('g_' + str(g) + '_' in k):                
                    temp_skew.append([k, dic[k][0], dic[k][1], dic[k][2]])

                if ('agg_' + str(agg) in k) and ('stat_FFh' in k) and ('g_' + str(g) + '_' in k):                
                    temp_dsf.append([k, dic[k][0], dic[k][1], dic[k][2]])
                    #temp_dsf.append([k, dic_dsf[k][0], dic_dsf[k][1], dic_dsf[k][2]])

        data_mean[str(g)] = sorted(temp_mean)
        data_var[str(g)] = sorted(temp_var)
        data_skew[str(g)] = sorted(temp_skew)
        data_dsf[str(g)] = sorted(temp_dsf)

            
    dic = wg_io.load(ifile)
    ts = 1
    lag = 1
    accumulate = True
    group_type = 'm'
    dic_stat = wg_tool.series_statistics(dic, ts, aggregation, lag=lag, 
                                  group_type=group_type, accumulate=accumulate)
    
    temp_mean = []
    temp_cv = []
    temp_var = []
    temp_dsf = []
    temp_skew = []
    temp_correl = []
                                  
    # [Ey, VARy, CVy, SKy, Rly, DSFy]
    dicdic = dic_stat['period']
    for g in range(12):
        mean_index = 0
        cf_mean = data_mean[str(g)]
        cf_mean = [i[1] for i in cf_mean]
        mean = dicdic[mean_index][g]
        new_mean = numpy.array(mean) * numpy.array(cf_mean)
        temp_mean.append(new_mean)
    
        var_index = 1
        cf_var = data_var[str(g)]
        cf_var = [i[1] for i in cf_var]        
        var = dicdic[var_index][g]
        new_var = numpy.array(var) * numpy.array(cf_var)
        temp_var.append(new_var)

        new_cv = (new_var**0.5)/new_mean
        temp_cv.append(new_cv)
    
        skew_index = 3
        cf_skew = data_skew[str(g)]
        cf_skew = [i[1] for i in cf_skew]
        skew = dicdic[skew_index][g]
        new_skew = numpy.array(skew) * numpy.array(cf_skew)
        temp_skew.append(new_skew)
                
        dsf_index = 5
        cf_dsf = data_dsf[str(g)]
        cf_dsf = [i[1] for i in cf_dsf]        
        dsf = dicdic[dsf_index][g]
        new_dsf = numpy.array(dsf) * numpy.array(cf_dsf)
        temp_dsf.append(new_dsf)
        
        
    ret_dic = {'mean': temp_mean,
               'var': temp_var,
               'skew': temp_skew,
               'cv': temp_cv,
               'dsf': temp_dsf,
               'aggregation': aggregation}        
        
    return ret_dic
    
if __name__ == '__main__':
    main()