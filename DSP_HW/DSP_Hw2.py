import scipy.io as sp
from scipy import signal
from scipy.io import wavfile
from scipy.io.wavfile import write
import matplotlib.pyplot as plot
import numpy as np
sound = sp.loadmat('./whistle.mat')
#region Problem 1
## Problem one
N = 10000; #the number of samples
T = 10 #seconds duration
t = np.linspace(0,T, N+1) # set up time variable to go from 0 to 10 seconds and include 10000 samples
f = 1.9 #define the frequency
# Answer Nyquist = 1/2*Sample rate so 1000 samples/secon*1/2 = 500 samples/sec is nyquist frequency
x = np.sin(2*np.pi*f*t)
# figure(1)
plot.plot(t,x)
plot.show()

#Problem 2
#downsample the signal by a factor of fifty
TD = t[::50]
x_new = np.sin(2*np.pi*f*TD)
plot.plot(t, x, 'g.-', TD, x_new, 'r.-')
plot.show()
# Answer Nyquist = 1/2*Sample rate so (1000/50)samples/second*1/2 = 10 samples/sec is nyquist frequency and the frequency of the singal would be 20Hz
TD = t[::500]
x_new = np.sin(2*np.pi*f*TD)
plot.plot(t, x, 'g.-', TD, x_new, 'b.-')
plot.show()
# Answer Nyquist = 1/2*Sample rate so (1000/500)samples/second*1/2 = 0.5 samples/sec is nyquist frequency and the frequency of the singal would be 2Hz


TD = t[::476]
x_new = np.sin(2*np.pi*f*TD)
plot.plot(t, x, 'g.-', TD, x_new, 'k.-')
plot.show()
# Answer Nyquist = 1/2*Sample rate so (1000/476)samples/second*1/2 = 1.05 samples/sec is nyquist frequency and the frequency of the singal would be 2.1Hz

#Problem 5
#the Nyquist Frequency is about 1 Hz and original frequency is 1.9 the new sample rate found by changing everything by a factor of 500 gives 2Hz
# by subtracting the new sample rate by the nyquist frequency we can find the aliasing frequency 2-1.9= 0.1 Hz
#in otherword by changing the sampling rate by a factor of 500hz you will alias some of the observable frequencies at 0.1 Hz


#Part 2 ######################################################
#Problem 1
Y=sound['Y']
Y = Y.ravel()
Fs = sound['Fs']
Fs = Fs.ravel()[0]
Y_new = np.linspace(Y[0],len(Y)/Fs,len(Y))
plot.plot(Y_new, Y)
plot.show() #plot the sound wave
write("sound_wave1.wav", Fs, Y) # Play this file in quicktime

# The Nyquist Frequency here is 22,050Hz

#Probelm 2 Decimate sample rate by 12
downsample = 12;
Ydec = signal.decimate(Y.ravel(), 12)
t = np.arange(0,(len(Ydec)/Fs)/12, 1/Fs/12)
plot.plot(t, Ydec)
plot.show()
fsdown = Fs/12
fsdown = fsdown.astype(np.int16)
write("sound_wave2_dec.wav", fsdown, Ydec)
# The new sample rate is changed by a factor of twelve so so the nyquist frequency is 1837.5 Hz
# Problem 3
Ydown = Y[::12] #downsample by a factor of twelve
t_down = np.linspace(0,len(Ydown)/(Fs/12),len(Ydown)) #this creates time array and we divide by 12 to create a sample rate changed by a factor of twelve
plot.plot(t_down, Ydown)
plot.show()
write("sound_wave3.wav", fsdown, Ydown)

#Problem 4:

#The reason the whistle sounds are different is because the frequency of each whistle is being changed by the rate at which they are samopled.
#The decimation step is changing the sample rate where as down sampling is just removing certain samples and inevitably aliases the high frequencies in the signal making the sound lower.
# the decimate function doesn't remove the high frequencies in problem two so that the signal was not aliased

#