import scipy.io as sp
import pandas as pd
from datetime import datetime, timedelta
import scipy as spy
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from numpy import diff
vol_data = sp.loadmat('./CO2_MaunaLoa.mat')
t_year = vol_data['ts']
CO2 = vol_data['CO2']
C02 = CO2.ravel()
t_year = t_year.ravel()



def date_num(dtyear):

    year = int(dtyear)
    time = dtyear
    rem = time-year

    base = datetime(year, 1, 1) #create the pandas datetime value
    result = base + timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem) #convert the remainder to the appropriate datetime
    print(result)
    return result


if __name__ == '__main__':
    result = []
    for a in range(len(t_year)):
        t = date_num(t_year[a])
        result.append([t])
    df = pd.DataFrame(result, columns=['Time'])

    df['year'] = pd.DatetimeIndex(df['Time']).year
    df['month'] = pd.DatetimeIndex(df['Time']).month


# Problem 1 Plot Time in correct series without spurious data points
medfilt_co2 = signal.medfilt(C02) # filter the Co2 data using a median filter
fig = plt.figure(1)
plt.plot(df['Time'], C02, "-b", label = "Unfiltered Data") # plot the original data
plt.plot(df['Time'], medfilt_co2,"-g", label = "MedFilt Data") #plot the data with median filter
plt.legend()
plt.title('C02 Muana Loa vs Time')
plt.xlabel('time')
plt.ylabel('Co2 Atmoshperic PPM')
plt.show()
#problem 2 ########
N = 30; #the number of samples
x = np.linspace(0,30,31) #This creates
f = 1/60 #define the frequency

h = np.sin(2*np.pi*f*x)
plt.plot(x,h)
plt.show()
norm_h = h/np.sum(h)
conv_data= np.convolve(medfilt_co2, norm_h, mode = 'same')
conv_test = conv_data[int(len(h)/2):len(conv_data) - int(len(h)/2)]
time_test = df.Time.values[int(len(h)/2):len(conv_data) - int(len(h)/2)]


fig, (ax1, ax2, ax3) = plt.subplots(3)
ax1.plot(df['Time'], C02, "-b", label = "Unfiltered Data")
ax1.plot(df['Time'], medfilt_co2,"-g", label = "MedFilt Data")
ax1.legend(loc='lower left')
# fig.legend((l1, l2), ['Unfiltered', 'MedFilt Data'])
ax1.set_xlabel('Time')
ax1.set_ylabel('Co2 ppm')
ax2.plot(time_test, conv_test)
ax2.set_xlabel('Time')
ax2.set_ylabel('Co2 ppm')





######Problem 3 #### # ignore my varibale names I mixed up the high pass and low pass
h_neg = -np.sin(2*np.pi*f*x)


norm = np.sum(h_neg) # set up norm of h_neg
norm_h = h_neg/norm # create the norm variable
norm_h = -norm_h # make it negative to reflect the filter we are using
midval = norm_h[int(len(norm_h)/2)] # save the midval point of the normalized array
low_h = np.where(norm_h == midval, 1+midval, norm_h) #replace the midval with 1 so that the integral is zero

low_pass = np.convolve(medfilt_co2, low_h, mode = 'same') #convolve the two functions
low_test = low_pass[int(len(low_h)/2):len(low_pass) - int(len(low_h)/2)]
lowtime_test = df.Time.values[int(len(low_h)/2):len(low_pass) - int(len(low_h)/2)]

ax3.plot(lowtime_test, low_test)
ax3.set_xlabel('Time')
ax3.set_ylabel('Co2 ppm')
plt.show()



######################### ANSWERS TO QUESTIONS #############################################

# 1. The median filter is removing the spurious data points, in this case the noise. As the moving window filter slides along the data it takes the median value of values within that window and applies it to the data
#    this smooths out the input signal in hand removing the outliers. This is a NON_LINEAR filter. It affects the peaks and troughs by smoothing out values around them but managing to maintian sharp changes if they are isloated in the window size.
#    If you new the median filter being used you would still not be able to retrieve the origninal data since it is a non-linear operation

# 2. If you combined the high pass and low pass filters you would see a signalthat is smoothed out from removing the high frequencies initally but it would maintina the high frequiencies occuring within the low frequency periods
#    so for the data from parts 2 and three it would be a relatively flat smooth wave with tiny frequent oscillations


# 3. The high-pass filter removes the DC offset as it detrends the amplitude of the wavefrom by taking the mean of the highest peak amplitude in this case
#    but maintains the frequency and signal of the input