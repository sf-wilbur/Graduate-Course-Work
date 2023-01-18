import numpy as np
import math
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt


def gauss_func( z, mu, sig):
    A = 1/(sig*np.sqrt(2*math.pi))
    B = (z-mu)**2
    C= 2*sig**2
    return A*np.exp(-B/C)



# with open('depths1.txt') as D1: #this reads in the depths1.txt file and assisngs it to the variable D1
#     depths1= D1.readlines()
# with open('depths2.txt') as D2:#this reads in the depths2.txt file and assisngs it to the variable D2
#     depths2 = D2.readlines()
#
# D1_array =np.array(depths1).astype(np.float)
# D2_array=np.array(depths2).astype(np.float)
#
# def function(D):
#     dmean = np.mean(D)
#     dstd =  np.std(D)
#     print ("mean:", dmean)
#     print ("std:", dstd)
#     return;
#
# if __name__ == '__main__':
#
#     D = D1_array
#
#     function(D)

