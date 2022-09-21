clear
% Read attached audio file with a sampling frequency of 48 kHz 
wav_file = 'eric.wav';
[yt,Fs] = audioread(wav_file);
% resample to 48 kHz
yt = resample(yt,48000,Fs);
Fs=48000;
%---------------------------FILTERING-------------------------------------%
% plot the audio file in time domain and in frequency domain
yf = fftshift(fft(yt));
x=-Fs/2:Fs/length(yt):Fs/2-Fs/length(yt); % frequency axis
filter_4k=ideal_filter(x,4000,1); % filter for 4000 Hz
figure(1)
plot(x,filter_4k) 
yf=yf.*filter_4k'; % filter the audio file
yt=ifft(ifftshift(yf)); % inverse fft to get the filtered audio file
time=0:1/Fs:length(yt)/Fs-1/Fs;
%--------------PLOT THE FILTERED AUDIO FILE AND PLAY SOUND----------------%
figure(2)
subplot(2,1,1)
plot(time,yt);   % plot the filtered audio file in time domain
title('Filtered audio file in time domain') 
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,1,2)
plot(x,abs(yf))    % filtered audio file in frequency domain
title('Filtered audio file in frequency domain') 
xlabel('Frequency (Hz)')
ylabel('Amplitude')
sound(ifft(ifftshift(yf)),Fs)
pause(length(yt)/Fs)
%--------------------------MODULATION---------------------------------------%
Fc=100000;
% resample the filtered audio file to 5*Fc
Fs_new=5*Fc;
yt1=resample(yt,Fs_new,Fs);
% DSB-SC Modulation 
%-----------------------Plot the modulated signal-----------------%
time=0:1/Fs_new:length(yt1)/Fs_new-1/Fs_new; %time for new sampling frequency
y_mod_t=yt1.*cos(2*pi*Fc*time'); % Modulated signal w/ no noise /with time transposed
length_of_filtered=length(yt1);
figure('Name','Second Modulated Signal','NumberTitle','off');
FvecIn_mod=linspace(-Fs_new/2,Fs_new/2,length_of_filtered);
subplot(2,1,1);
plot(FvecIn_mod,fftshift(abs(fft(y_mod_t)))/length_of_filtered);
title('Amplitude Spectrum of modulated signal(t)')
xlabel('Frequency(Hz)')
ylabel('Amplitude Spectrum(f)')
subplot(2,1,2);
plot(time',real(y_mod_t));
title('modulated signal(t)')
xlabel('time(s)')
ylabel('Amplitude')
%--------------------------DEMODULATION------------------------------------%
noise=[0,10,30];
for(i=1:3)
time=0:1/Fs_new:length(yt1)/Fs_new-1/Fs_new; %time for new sampling frequency
%demodulate the signal using coherent demodulation with SNR=0
y_noise_mod_t=awgn(y_mod_t,noise(i),'measured'); % add noise to the modulated signal
y_noise_demod_t=y_noise_mod_t.*2.*cos(2*pi*Fc*time'); % demodulate the signal
y_noise_demod_f=fftshift(fft(y_noise_demod_t)); % fft of the demodulated signal
%---------------------COHERENT FILTER------------------------------------------%
x=-Fs_new/2:Fs_new/length(y_noise_demod_t):Fs_new/2-Fs_new/length(y_noise_demod_t); %x-axis for frequency domain
filter_coherent=ideal_filter(x,Fs/2,1); % create filter for coherent demodulation
y_demod_noise_filter_f=y_noise_demod_f.*filter_coherent'; % filter the demodulated signal
y_noise_demod_t=ifft(ifftshift(y_demod_noise_filter_f)); %demodulated signal in time domain
y_noise_demod_t=resample(y_noise_demod_t,Fs,Fs_new); % resample the demodulated signal to Fs (original sampling frequency) 
y_demod=fftshift(fft(y_noise_demod_t)); % fft of the demodulated signal
%--------------------------PLOT THE RECEIVED SIGNAL----------------------%
% plot the demodulated signal in time domain
time=0:1/Fs:length(y_noise_demod_t)/Fs-1/Fs; 
x=-Fs/2:Fs/length(y_noise_demod_t):Fs/2-Fs/length(y_noise_demod_t); % frequency axis
figure
subplot(2,1,1)
plot(time,y_noise_demod_t) % plot the demodulated signal in time domain
title(sprintf('Demodulated with SNR=%d db\n in time domain',noise(i)))
xlabel('Frequency (Hz)')
ylabel('Amplitude')
subplot(2,1,2)
plot(x,abs(y_demod)) % plot the demodulated signal in frequency domain   
title(sprintf('Demodulated with SNR=%d db\n in frequency domain',noise(i)))
xlabel('Frequency (Hz)')
ylabel('Amplitude')
%-------------------9-------PLAY THE RECEIVED SIGNAL----------------------%
%sound(y_noise_demod_t,Fs)
%pause(length(y_noise_demod_t)/Fs)
end

%----------%---------------Freq ERROR Fc=100.1kHZ and snr =10db------------%---------%

%--------------------------MODULATION---------------------------------------%
Fc=100000;
Fc_new=100100;
% resample the filtered audio file to 5*Fc
%Fs_new=5*Fc_new;
Fs_new=5*Fc;
yt2=resample(yt,Fs_new,Fs);
% DSB-SC Modulation
%-----------------------Plot the modulated signal-----------------%
time=0:1/Fs_new:length(yt1)/Fs_new-1/Fs_new; %time for new sampling frequency
y_mod_t=yt2.*cos(2*pi*Fc*time'); % Modulated signal w/ no noise /with time transposed
length_of_filtered=length(yt2);
figure('Name','Third Modulated Signal','NumberTitle','off');
FvecIn_mod=linspace(-Fs_new/2,Fs_new/2,length_of_filtered);
subplot(2,1,1);
plot(FvecIn_mod,fftshift(abs(fft(y_mod_t)))/length_of_filtered);
title('Amplitude Spectrum of modulated signal(t)')
xlabel('Frequency(Hz)')
ylabel('Amplitude Spectrum(f)')
subplot(2,1,2);
plot(time',real(y_mod_t));
title('modulated signal(t)')
xlabel('time(s)')
ylabel('Amplitude')

%--------------------------DEMODULATION------------------------------------%
time=0:1/Fs_new:length(yt2)/Fs_new-1/Fs_new; %time for new sampling frequency
%demodulate the signal using coherent demodulation with SNR=0
y_noise_mod_t=awgn(y_mod_t,10,'measured'); % add noise to the modulated signal
y_noise_demod_t=y_noise_mod_t.*2.*cos(2*pi*Fc_new*time'); % demodulate the signal
y_noise_demod_f=fftshift(fft(y_noise_demod_t)); % fft of the demodulated signal
%---------------------COHERENT FILTER------------------------------------------%
x=-Fs_new/2:Fs_new/length(y_noise_demod_t):Fs_new/2-Fs_new/length(y_noise_demod_t); %x-axis for frequency domain
filter_coherent=ideal_filter(x,Fs/2,1); % create filter for coherent demodulation
y_demod_noise_filter_f=y_noise_demod_f.*filter_coherent'; % filter the demodulated signal
y_noise_demod_t=ifft(ifftshift(y_demod_noise_filter_f)); %demodulated signal in time domain
y_noise_demod_t=resample(y_noise_demod_t,Fs,Fs_new); % resample the demodulated signal to Fs (original sampling frequency) 
y_demod=fftshift(fft(y_noise_demod_t)); % fft of the demodulated signal
%--------------------------PLOT THE RECEIVED SIGNAL----------------------%
% plot the demodulated signal in time domain
time=0:1/Fs:length(y_noise_demod_t)/Fs-1/Fs; 
x=-Fs/2:Fs/length(y_noise_demod_t):Fs/2-Fs/length(y_noise_demod_t); % frequency axis
figure
subplot(2,1,1)
plot(time,y_noise_demod_t) % plot the demodulated signal in time domain
title(sprintf('Demodulated with SNR=%d db\n in time domain and Freqency error',10))
xlabel('Frequency (Hz)')
ylabel('Amplitude')
subplot(2,1,2)
plot(x,abs(y_demod)) % plot the demodulated signal in frequency domain   
title(sprintf('Demodulated with SNR=%d db\n in frequency domain and Freqency error',10))
xlabel('Frequency (Hz)')
ylabel('Amplitude')
%-------------------9-------PLAY THE RECEIVED SIGNAL----------------------%
sound(y_noise_demod_t,Fs)
pause(length(y_noise_demod_t)/Fs)

%----------%---------------Phase ERROR 20 deg and snr =10db------------%---------%


%--------------------------MODULATION---------------------------------------%
Fc=100000;
%Fc_new=100100;
% resample the filtered audio file to 5*Fc
%Fs_new=5*Fc_new;
Fs_new=5*Fc;
yt2=resample(yt,Fs_new,Fs);
% DSB-SC Modulation
time=0:1/Fs_new:length(yt2)/Fs_new-1/Fs_new; %time for new sampling frequency
y_mod_t=yt2.*cos(2*pi*Fc*time'); % Modulated signal w/ no noise /with time transposed

%--------------------------DEMODULATION------------------------------------%
time=0:1/Fs_new:length(yt2)/Fs_new-1/Fs_new; %time for new sampling frequency
%demodulate the signal using coherent demodulation with SNR=0
y_noise_mod_t=awgn(y_mod_t,10,'measured'); % add noise to the modulated signal
y_noise_demod_t=y_noise_mod_t.*2.*cos(2*pi*Fc*time' +((pi*20)/180)); % demodulate the signal
y_noise_demod_f=fftshift(fft(y_noise_demod_t)); % fft of the demodulated signal
%---------------------COHERENT FILTER------------------------------------------%
x=-Fs_new/2:Fs_new/length(y_noise_demod_t):Fs_new/2-Fs_new/length(y_noise_demod_t); %x-axis for frequency domain
%filter_coherent=ideal_filter(x,Fs/2,1); % create filter for coherent demodulation
filter_coherent=ideal_filter(x,4000,1); % create filter for coherent demodulation
y_demod_noise_filter_f=y_noise_demod_f.*filter_coherent'; % filter the demodulated signal
y_noise_demod_t=ifft(ifftshift(y_demod_noise_filter_f)); %demodulated signal in time domain
y_noise_demod_t=resample(y_noise_demod_t,Fs,Fs_new); % resample the demodulated signal to Fs (original sampling frequency) 
y_demod=fftshift(fft(y_noise_demod_t)); % fft of the demodulated signal
%--------------------------PLOT THE RECEIVED SIGNAL----------------------%
% plot the demodulated signal in time domain
time=0:1/Fs:length(y_noise_demod_t)/Fs-1/Fs; 
x=-Fs/2:Fs/length(y_noise_demod_t):Fs/2-Fs/length(y_noise_demod_t); % frequency axis
figure
subplot(2,1,1)
plot(time,y_noise_demod_t) % plot the demodulated signal in time domain
title(sprintf('Demodulated with SNR=%d db\n in time domain and phase error 20 deg',10))
xlabel('Frequency (Hz)')
ylabel('Amplitude')
subplot(2,1,2)
plot(x,abs(y_demod)) % plot the demodulated signal in frequency domain   
title(sprintf('Demodulated with SNR=%d db\n in frequency domain and phase error 20 deg',10))
xlabel('Frequency (Hz)')
ylabel('Amplitude')
%-------------------9-------PLAY THE RECEIVED SIGNAL----------------------%
sound(y_noise_demod_t,Fs)
pause(length(y_noise_demod_t)/Fs)