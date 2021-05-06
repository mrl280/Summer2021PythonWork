import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
import numpy as np

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams['figure.figsize'] = 20, 10

#create pulse mask
time = np.arange(10000)
amplitude = np.ones(10000)

amplitude[1000:2000] = 0

amplitude[3000:5000] = 0

amplitude[6000:10000] = 0

plt.title('Pulse Mask')
plt.plot(time,amplitude)
plt.show()

plt.title('Pulse Mask Spectrum')
plt.plot(np.abs(np.fft.fft(amplitude)))
plt.show()

#mix the mask with the TX frequency
amplitude = amplitude*np.exp(1.0j*2*np.pi*time/100)

plt.title('Pulses')
plt.plot(time,amplitude)
plt.show()

plt.title('Pulse Spectrum')
plt.plot(np.abs(np.fft.fft(amplitude)))
plt.show()

#add a Doppler shift and time delay to the transmitted pulses

#time delay of 100 ms
amplitude = np.roll(amplitude,1000)

#Doppler shift of 1 Hz
amplitude = amplitude*np.exp(1.0j*2*np.pi*time/10000)

#sample the waveform
sample_points= 4*np.arange(2499)

sample_values = np.real(amplitude[sample_points])

#25 points is the quarter wavelength of the center frequency of 100 Hz
sample_values_complex = np.real(amplitude[sample_points[:-10]+25])

plt.title('Sampled Signal')
plt.plot(time,np.real(amplitude))
plt.scatter(sample_points,(sample_values))
plt.show()

#example of complex sampling
plt.title('Sampled Complex Signal')
plt.plot(time,np.real(amplitude))
plt.scatter(sample_points,(sample_values))
plt.scatter(sample_points[:-10]+25,(sample_values_complex))
plt.show()

plt.title('Sampled Signal Spectrum')
plt.plot(np.abs(np.fft.fft(sample_values)))
plt.show()

#mix to baseband
sample_values = sample_values*np.exp(-1.0j*2*np.pi*sample_points/100)

plt.title('Mixed Sampled Signal Spectrum')
plt.plot(np.abs(np.fft.fft(sample_values)))
plt.show()

#low pass filter
#typically a low pass filter is performed here to prevent unwanted signals from being aliased in
sample_values = np.convolve(sample_values, np.ones(20)/20, mode='full')

plt.title('Mixed, Filtered, Sampled Signal Spectrum')
plt.plot(np.abs(np.fft.fft(sample_values)))
plt.show()

#decimate
dec_sample_values = sample_values[np.arange(100)*25]

plt.title('Decimated Sampled Signal')
plt.scatter(np.arange(100)*25,np.real(dec_sample_values))
plt.scatter(np.arange(100)*25,np.imag(dec_sample_values))
plt.show()

plt.title('Decimated Sampled Signal Spectrum')
plt.plot(np.abs(np.fft.fft(dec_sample_values)))
plt.show()

#create rawACF lags
lag_x = [0,2,3,5]
lag_values = np.zeros(4,dtype=np.complex64)
lag_values[0]=1.0
lag_values[1]=dec_sample_values[15]*np.conj(dec_sample_values[35])
lag_values[2]=dec_sample_values[35]*np.conj(dec_sample_values[65])
lag_values[3]=dec_sample_values[15]*np.conj(dec_sample_values[65])

plt.title('rawACF Example (phase)')
plt.scatter(lag_x,np.angle(lag_values))
plt.show()
