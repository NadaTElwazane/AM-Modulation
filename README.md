# DSB-SC Modulation using MATLAB
## Background
Double Sideband modulation is the easiest and most direct type of analog modulation. In
this scheme, the modulated signal is obtained using a direct multiplication of the
modulating signal (i.e. the message) by a cosine carrier. If the carrier term is omitted,
the modulation is termed double sideband suppressed carrier (DSB-SC). DSB-TC has a
significant advantage in the receiver design (i.e. the envelop detector). However, the DSB-TC loses to
the other variant (i.e. the SC) in terms of power efficiency.
## Procedure
- Use Matlab to read the attached audio file ("eric.wav") which has a sampling frequency Fs= 48 KHz, plot the spectrum of this signal.
- Use an ideal filter to remove all frequencies greater than 4KHZ.
- Obtain the filtered signal in both frequency and time domain.
- Sound the filtered audio signal (make sure that there is only a small error in the filtered signal)
- Generate DSB-SC modulated signal using a carrier signal with carrier frequency Fc=100 KHz, (Note : you have to increase the sampling frequency of filtered signal to be (Fs=5*Fc) before modulation process)
- Demodulate the modulated signal using a coherent detection receiver with SNR=0 ,10, 30 dB.
- For the three SNR cases; sound the received signal and plot it in both time and frequency domain, what conclusion you make of that? (Note: to sound signal after demodulation process you have to decrease the sampling frequency again).
- Repeat the coherent detection with frequency error, Fc=100.1 KHz instead of 100 KHz at receiver with SNR=10dB.
- Repeat the coherent detection with phase error=200
at receiver with SNR=10db.

## Contributors
- Nada Alaa Mohamed
- [Nada Elwazane](https://github.com/NadaTElwazane)