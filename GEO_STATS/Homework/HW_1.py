
#import the appropriate  modules to compute equations, statistics, and populate graphs
from scipy import stats
from scipy.stats import norm


import numpy as np
import math
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt


def gauss_func( z, mu, sig):   #define the gaussian function I'll use later in the script, (z = input data, mu = average value, and sig = standard deviation
    A = 1/(sig*np.sqrt(2*math.pi))    # necessary equations used in the equation following the return command
    B = (z-mu)**2
    C= 2*sig**2
    return A*np.exp(-B/C)   #return command that provides output for gaussian function


if __name__ == '__main__':   #this makes the script run the function above in the current working directory


    with open("depths1.txt") as D1: #this reads in the depths1.txt file and assisngs it to the variable D1
        depths1= D1.readlines()
    with open("depths2.txt") as D2:#this reads in the depths2.txt file and assisngs it to the variable D2
        depths2 = D2.readlines()

D1_array =np.array(depths1).astype(np.float)    #this assigns a variable to to the numpy command that makes depths1 an array
D2_array =np.array(depths2).astype(np.float)    #this assigns a variable to the numpy command that makes depths2 an array

D1_mean = np.mean(D1_array)                     # these two lines of code compute the mean for the Depths1 and Depths2 data
D2_mean = np.mean(D2_array)

D1_std =np.std(D1_array)                        #copmute the standard deviation for the Depths1 and 2 data
D2_std = np.std(D2_array)

D1_med = np.median(D1_array)                    #compute the median for the data using np.median
D2_med = np.median(D2_array)
D1_mode = stats.mode(D1_array)                  #compute the mode uscing scipy.stats.mode
D2_mode = stats.mode(D2_array)

D1_IQR =sp.stats.iqr(D1_array)                  #compute the inter quartile range  for the data tables
D2_IQR=sp.stats.iqr(D2_array)
D1_skew = sp.stats.skew(D1_array)               #compute the skew of each data file
D2_skew = sp.stats.skew(D2_array)
D1_k = sp.stats.kurtosis(D1_array)              #computes the kurtosis for each data file
D2_k = sp.stats.kurtosis(D2_array)
print('The mean is %s, median is %s, standard deviation is %s, mode is %s, IQR is %s, skew is %s, and kurtosis is %s for the depths1 data set ' %(D1_mean, D1_med, D1_std, D1_mode, D1_IQR, D1_skew, D1_k))
print('The mean is %s, median is %s, standard deviation is %s, mode is %s, IQR is %s, skew is %s, and kurtosis is %s for the depths2 data set ' %(D2_mean, D2_med, D2_std, D2_mode, D2_IQR, D2_skew, D2_k))
plt.figure(1)                                   # this creates figure 1 that will be given the input  to create boxplot using matplotlib
plt.title("Notched Box Plot")
plt.boxplot(D1_array, notch = True)             #notch = true is what creates the notch in the figure
plt.show()

plt.figure(2)
plt.boxplot(D2_array, notch = True)            # this does the samwe as the code above except for the depths2 array
plt.title("Notched Box Plot")
plt.show()

plt.figure(3)
plt.hist(D1_array, bins =30)                    #This creates a histrogram with counts on the y axis for the depths1 array with 30 bins
plt.title('D1_Hist')
plt.show()

plt.figure(4)
plt.hist(D2_array, bins= 30)                    #This creates a histrogram with counts on the y axis for the depths2 array with 30 bins
plt.title('D2_Hist')
plt.show()


fig, axs = plt.subplots(2, 2)                       #This line of code sets up a subplot with 2 columns and 2 rows
axs[0,0].boxplot(D1_array, notch = True)               #Axis[0,0] is the position of the first plot in the subplopt
axs[0,0].set_ylim([0,400])                          #this line sets the y limit range to be 400
axs[0,0].set_title("Notched Box Plot_Depths1")



axs[0,1].boxplot(D2_array, notch = True)                    #this sets the position of the boxplot to being in the first row second column
axs[0,1].set_ylim([0,400])                                  # set y limit
axs[0,1].set_title("Notched Box Plot_Depths2")


axs[1,0].hist(D1_array, bins =30)                       #This sets the position of the histrogram to row 2 column 1
axs[1,0].set_xlim([0,400])                                 #sets x axis limits
axs[1,0].set_xlabel("Depths1_Histogram")

axs[1,1].hist(D2_array, bins =30)                   # sets the position for the histogram in row 2 column 2
axs[1,1].set_xlim([0,400])                          #sets x axis limit
axs[1,1].set_xlabel("Depths2_Histogram")
plt.show()


#relative histogram function
plt.figure(6)
hist, bins = np.histogram(D1_array, bins =30, density = True)  #This creates a relative density histogram by setting density= true
widths = np.diff(bins)                                          #find difference between bin widths to get center of bins
normed_value= 1                                                 #set the normalized value
hist *= normed_value                                            #multiplies the stores hist variable by the normed_value and stores hist again
plt.bar(bins[:-1], hist, widths)                                #plots the normalzied histogram using the bar function
#still need to do this for D2_array

#gauss function for D1_array
x = gauss_func(sorted(D1_array), D1_mean, D1_std)  #calls the gaussian function written at the beginning of the script
plt.plot(sorted(D1_array),x, 'red')                 #plots a gaussian curve ontop of my hist function in red
plt.show()
#Reassign Variables:
mu = D1_mean                    #reassigns variables to the mean and standard deviation values for my depths1 array
sig = D1_std
prob = norm(mu, sig).cdf(mu+20) - norm(mu, sig).cdf(mu-20)   #computes the probability of finding value within 20cm of the mean (mu)



print('the probability of the value being less than or greater than 20cm of the mean for depths1 is %s' % prob)

plus_20 =norm(mu, sig).cdf(100000) - norm(mu, sig).cdf(mu+20)
print('the probability of the value being  greater than 20cm of the mean for depths1 is %s' % plus_20)
minus_20 = norm(mu, sig).cdf(mu-20) - norm(mu, sig).cdf(-1000000)
print('the probability of the value being less than 20cm of the mean for depths1 is %s' % minus_20)





plt.figure(7)
hist_2, bins_2 = np.histogram(D2_array, bins =30, density = True)
widths_2 = np.diff(bins_2)
normed_value_2= 1
hist_2 *= normed_value_2
plt.bar(bins_2[:-1], hist_2, widths_2)
plt.xlim(0,600)
x_2 = gauss_func(sorted(D2_array), D2_mean, D2_std)
plt.plot(sorted(D2_array),x_2, 'red')
plt.xlim(0,600)
plt.show()
#Reassign Variables:
mu2 = D2_mean       #reassigns variables to the mean and standard deviation values for my depths2 array
sig2 = D2_std
prob_2 = norm(mu2, sig2).cdf(mu2+20) - norm(mu2, sig2).cdf(mu2-20)    #computes the probability of finding value within 20cm of the mean (mu)

#Reassign Variable

prob_2 = norm(mu2, sig2).cdf(mu2+20) - norm(mu2, sig2).cdf(mu2-20)
print('the probability of the value being less than or greater than 20cm of the mean for depths2 is %s' % prob_2)
plus_20_d2 = 1 - norm(mu2, sig2).cdf(mu2+20)
print('the probability of the value being  greater than 20cm of the mean for depths2 is %s' % plus_20_d2)
minus_20_d2 = norm(mu2, sig2).cdf(mu2-20)
print('the probability of the value being less than 20cm of the mean for depths2 is %s' % minus_20_d2)

