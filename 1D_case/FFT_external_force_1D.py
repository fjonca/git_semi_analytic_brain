"""
Semi-analytic solving of a hyperelastic bar impact (1D)

This code tests Fourier decomposition of external force to identify 
low-frequency part (quasi-static) and vibrating one (dynamic) 

Author: François Jonca
"""

#%% Libraries import
##from scipy.fft import fft, ifft, fftfreq
##from scipy.signal import find_peaks
##import matplotlib.pyplot as plt
##import numpy as np



#%% Fast-Fourier transform of f(t)
# sample period
T = tf/n
# 
fFreq = fft(f); #print(size(fFreq))
xf = fftfreq(n, T)[:n//2] #Hz

#Plot spectrum
figure(2)
plot(xf, 2.0/n * abs(fFreq[0:n//2]),label = 'External force spectrum')
legend(loc='best')
grid()

#Compare with analytic spectrum

#%% Extraction of 0-Hz component
filter0Hz = zeros(size(fFreq)); filter0Hz[0] = 1
staticFFreq = fFreq*filter0Hz

figure(2)
plot(xf, 2.0/n * abs(staticFFreq[0:n//2]),'r--', label = '0-Hz component')
legend()

zeroHzComponent = ifft(staticFFreq); zeroHzComponent = zeroHzComponent.real

figure(1)
plot(t,zeroHzComponent,'r--',label = '0-Hz component')
legend()

#%% Extraction of vibrating component : 0-Hz component is suppressed
# dynamicFFreq = np.zeros(np.size(fFreq)); dynamicFFreq[1:-1] = fFreq[1:-1];
filterDynamic = ones(size(fFreq)); filterDynamic[0] = 0
dynamicFFreq = fFreq*filterDynamic
#Plot
figure(2)
plot(xf, 2.0/n * abs(dynamicFFreq[0:n//2]),'b--', label = 'Vibrating component')
legend()

vibratingComponent = ifft(dynamicFFreq); vibratingComponent = vibratingComponent.real
#Plot
figure(1)
plot(t,vibratingComponent,'b--',label = 'Vibrating component')
legend()

#%% Other way : extract before and after the first zero-energy value
# spectrum = fFreq
# index = find_peaks(-2.0/n* np.abs(spectrum)); print(index[0][0])

# plt.figure(2)
# plt.plot(xf[index[0][0]], 2.0/n* np.abs(spectrum[index[0][0]]),'ro')

# #Part before first zero-energy frequency
# filterPart1 = np.zeros(np.size(fFreq)); filterPart1[0:index[0][0]+1] = 1
# part1 = fFreq*filterPart1
# #Part after first zero-energy frequency
# filterPart2 = np.zeros(np.size(fFreq)); filterPart2[index[0][0]+1:-1] = 1
# part2 = fFreq*filterPart2
# # part2 = spectrum; part2[0:index[0][0]+1] = 0

# plt.figure(3)
# plt.plot(xf, 2.0/n * np.abs(fFreq[0:n//2]),label = 'External force spectrum')
# plt.plot(xf, 2.0/n * np.abs(part1[0:n//2]),'k--',label = 'Part 1')
# plt.plot(xf, 2.0/n * np.abs(part2[0:n//2]),'r--',label = 'Part 2')
# plt.legend(loc='best')
# plt.grid()

# #Inverse transform
# staticComponent2 = ifft(part1); staticComponent2 = staticComponent2.real
# vibratingComponent2 = ifft(part2); vibratingComponent2 = vibratingComponent2.real

# plt.figure(1)
# plt.plot(t,staticComponent2,'k--',label = 'Part 1')
# plt.plot(t,vibratingComponent2,'g--',label = 'Part 2')

# plt.legend(loc = 'best')

#%% Show plots
