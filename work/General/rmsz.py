#!/usr/bin/env python

def rmsz(data, mu, sigma):
    '''
    Compute the Root-Mean-Square of the z-scorw for a set of data
    using z_i = {x_i - mu / sigma}

    mu and sigma are given to the function and not calculated from the data
    
    Arguments:
    :param data: list of real positive numbers
    :param mu: mean of the population
    :param sigma: the standard deviation of the population
    :returns: (sum[(z_i)^2] / length[data])^0.5
    '''
    if len(data) == 0:
        raise ValueError('Called rmsz() with no data !')
    
    if sigma == 0:
        raise ValueError('Called rmsz() with zero sigma value')
    
    
    z_2 = [pow((x-mu)/sigma,2) for x in data if x>0]
    return pow(sum(z_2)/len(z_2),0.5)

    