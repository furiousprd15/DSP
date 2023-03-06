%Prewarped frequency conversion from digital to analog(\tan\left(\frac{W\pi}{2}\right)\cdot\frac{2}{a})

fp1=1.70816137093; fp2=	2.34169913223; fs1=1.45308505601; fs2=3.26370337426;

% Angular analog freuency for passband and stopband
Wp1 = 2*pi*fp1; Wp2 = 2*pi*fp2; Ws1 = 2*pi*fs1; Ws2 = 2*pi*fs2;


%Attenuation at passband and stopband levels(Ripple at passband taken to be 0.5dB)
Ap = 0.5;As = 40;
% Absolute gain of transfer function at passband and stopband
Gp = 10^(-Ap/20); 
Gs = 10^(-As/20); 

%Finding Passband and Stopband critical frequencies for prototype lowpass
%filter
Wp = Wp2-Wp1; W0 = sqrt(Wp1*Wp2);
W1 = Ws1 - W0^2/Ws1; W2 = Ws2 - W0^2/Ws2;
Ws = min(abs([W1,W2]));

%Using the custom built low pass analog filter script file for calculationg
%the coeffient matrix of our analog filter
[N,B,A,za,pa] = lpa(Wp,Ws,Ap,As);


%Frequencxy sample points for finding gain vs frequency plot
f = linspace(0,5,1001); s = 1i*2*pi*f;
% Transformation from lowpass to bandpass
s=s+ W0^2./s; 
H = 1;
% Each row of A matrix is of the format = A(i,1) + A(i,2)*s +
% A(i,3)*s.^2.So we loop over the entire rows of A(Denominator) and
% B(Numerator),multiply them to find the transfer function at each frequency sample point.

for i=1:size(B,1)
H = H .* (B(i,1) + B(i,2)*s + B(i,3)*s.^2) ./ (A(i,1) + A(i,2)*s + A(i,3)*s.^2);
end

%Magnitude plot,logMagnitude plot and phase plot of analog bandpass filter
y1 = 20*log10(abs(H));
y2 =angle(H);
y3 = abs(H);

plot(f,y1,'r')
plot(f,y2,'r')
plot(f,y3,'r')


% Checking whether the filter satisfies the bandpass requirement
index = find(f ==1.45 );   
yDesired = y1(index);

%System Transfer function
sys1 = tf(flip(B(1,:)),flip(A(1,:)));
sys2 = tf(flip(B(2,:)),flip(A(2,:)));
sys3 = tf(flip(B(3,:)),flip(A(3,:)));
sys_  = sys1*sys2*sys3;
%Impulse Response,and group delay
[num,den] = tfdata(sys_,'v')
impz(num,den,10);
grpdelay(num,den,10);

[numd,dend] = bilinear(num,den,1);
freqz(numd,dend,1024);


%       8.381e-05 s^4 + 0.0238 s^2 + 0.9441
%  ------------------------------------------------------
%  0.008381 s^4 + 0.03933 s^3 + 0.2332 s^2 + 0.5733 s + 1
 
%Continuous-time transfer function.


%{

Faulty Bilinear transform 
w_d = pi*linspace(0,1,1001);
z_v = exp(1i * w_d);
s = 2*(-1+z_v)./(1+z_v);

H = 1;
for i=1:size(B,1)
H = H .* (B(i,1) + B(i,2)*s + B(i,3)*s.^2) ./ (A(i,1) + A(i,2)*s + A(i,3)*s.^2);
end
plot(w_d/pi,20*log10(abs(H)),'r')


%}


