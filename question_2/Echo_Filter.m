info = audioinfo("dataset_6/goodmorning.wav");
[y,Fs] = audioread("dataset_6/goodmorning.wav"); % reading audio file with 
% Fs = 8000 Hz

t = 0:seconds(1/Fs):seconds(info.Duration);
t = t(1:end-1);

% xlabel('Time')
% ylabel('Audio Signal')
sound(y,Fs) % original sound
Y = fft(y);
f = linspace(0,Fs,length(Y));

[Rmm,lags] = xcorr(y,'unbiased'); % autocorrelation of signal Y%

Rmm = Rmm(lags>0);
lags = lags(lags>0);

figure
plot(lags/Fs,Rmm)
xlabel('Lag (s)')


[~,dl] = findpeaks(Rmm,lags,'MinPeakHeight',0.001126); % find peaks in autocorrelation function
dl
alpha =  0.2863; % setting alpha from calculations
mtNew = filter(1,[1 zeros(1,dl(4)-1) alpha],y); % making a filter using given coefficients


hold on
% plot(t,y)
% subplot(2,1,2)
% plot(t,mtNew)

pause(3);
soundsc(mtNew,Fs) % filtered sound


[FRmm, lagsFin] = xcorr(mtNew, 'unbiased'); % autocorrelation of filtered sound


plot(lagsFin/Fs,FRmm)
ratio_val = FRmm(dl(4)) / Rmm(dl(4))
hold off
filename = 'FinalOutputAudio.wav'; % writin
audiowrite(filename,mtNew, Fs); % writing filtered file



