import scipy.io as sp

from scipy import signal
from scipy.io import wavfile
from scipy.io.wavfile import write
import matplotlib.pyplot as plt
import numpy as np
from numpy import diff
seis = sp.loadmat('./Erebus_seismogram.mat')

xdata = seis['data']

hdr = seis['hdr'] #signal info
sr = hdr['sps'] #sample rate
sr = np.ravel(sr)
conv = hdr['atod'] # conversion factor Volts/count
sensor = 3200 #sensor sensitivity in v/m/s



data_dtrend = xdata - np.mean(xdata)
new_data = data_dtrend*conv
data_mps = new_data/sensor # gets data in m/s
data_mmps = data_mps*1000 # convert to mm/s
t = np.linspace(data_mmps[0],len(data_mmps)/sr, len(data_mmps))
t = t.ravel()
fig, (ax1, ax2, ax3) = plt.subplots(3)

ax1.plot(t, data_mmps)
ax1.set_title('Problem 1')

# set labels

plt.setp(ax1, ylabel='Velocity mm/s')

############################# PROBLEM 2##########################

###################### Problem 2 #############################

differ = np.diff(data_mmps.ravel())*sr
t = np.linspace(data_mmps[0],(len(data_mmps)-1)/sr, (len(data_mmps)-1))
t = t.ravel()
ax2.plot(t,differ)
ax2.set_title('Problem 2')

plt.setp(ax2, ylabel='Acceleration (mm/s^2) ')
######## PROBLEM 3##########################################

cum_sum = np.cumsum(data_mmps/sr)
t = np.linspace(data_mmps[0],len(data_mmps)/sr, len(data_mmps))
t = t.ravel()
ax3.plot(t, cum_sum)
ax3.set_title('Problem 3')


plt.setp(ax3, ylabel='Displacement mm')
plt.show()