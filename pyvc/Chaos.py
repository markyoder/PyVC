#!/usr/bin/env python
import numpy as np
import hashlib
import math
import datetime
from scipy import stats

def permutationEntropy(data, word_length):
    it = np.nditer(data, flags=['c_index'])
    permutation_count = {}
    while not it.finished and it.index + word_length <= data.size:
        #try:
        perm = np.argsort(data[it.index:it.index + word_length])
        #except IndexError:
        #    print "test"
        #    break
        
        #print perm, it.index, data.size
        
        try:
            permutation_count[ hashlib.sha1(perm).hexdigest() ] += 1.0
        except KeyError:
            permutation_count[ hashlib.sha1(perm).hexdigest() ] = 1.0
           
        
        it.iternext()
    
    H = 0
    
    total_words = float(data.size - word_length + 1)
    
    for perm_num in permutation_count.itervalues():
        i = perm_num/total_words
        #print perm_num, total_words, i, math.log(i,2)
        H -= i * math.log(i,2) 
    
    return H
    #print data, word_length

def corrCoeff(x, y):
# x and y are the two 1D arrays to be correlated
# time is an array in the same shape as x and y
    #delta_t = time[-1] - time[0]
    
    #mean_x = (1.0/delta_t.days) * np.sum( x )
    #mean_y = (1.0/delta_t.days) * np.sum( y )
    
    #var_x = (1.0/delta_t.days) * np.sum( (x - mean_x)**2.0 )
    #var_y = (1.0/delta_t.days) * np.sum( (y - mean_y)**2.0 )
    
    std_x = np.std( x )
    std_y = np.std( y )
    
    mean_x = np.mean( x )
    mean_y = np.mean( y )
    
    #print mean_x, mean_y, std_x, std_y, len(x)
    
    #print np.sum( (x - mean_x) * (y - mean_y) ), (1.0/(std_x * std_y)),
    
    return (1.0/(std_x * std_y)) * (1.0/float(len(x))) * np.sum( (x - mean_x) * (y - mean_y) ), np.corrcoef(x,y)[0,1]
    

    
    