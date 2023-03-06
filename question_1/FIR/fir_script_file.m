% filter specifications

% pass-band w
wd_p1 = 0.45*pi; % in rad.
wd_p2 = 0.55*pi;
A_p = 0.5; % in dB

% stop-band w
wd_s1 = 0.40*pi;
wd_s2 = 0.65*pi;
A_s1 = 40;
A_s2 = 50;

% choose a window
K = 8; %for hanning
%K = 4; %for rectangular

% find N --- keep it odd
w_t = min((wd_s2 - wd_p2), (wd_p1- wd_s1)); % transition length

N = ceil((2*K*pi)/w_t); % N points, order of filter

if mod(N,2) == 0 % make it odd
    N = N+1;
end

%N = 13;

% find alpha
% symmetric about it
alpha = (N-1)/2; % matlab indexes from 1 instead of 0

% ideal brick-wall filter specs for window
w_l = wd_p1;
w_h = wd_p2;



% find desired hd[n] for n from 1 to N
hd_coeff = ones(1,N);

for n = 0:N-1
    if n == alpha
        hd_coeff(n+1) = (w_h-w_l)/pi;
    else
        hd_coeff(n+1) = sin(w_h*(n-alpha))/(pi*(n-alpha)) - sin(w_l*(n-alpha))/(pi*(n-alpha));
    end
end

% find window function for given N
window_coeff = ones(1,N);
for n = 0:(N-1)
    window_coeff(n+1) = (0.5*(1-cos((2*pi*n)/(N-1))));
end

% find h_coeff
h_coeff = ones(1,N);
for n = 0:N-1
    h_coeff(n+1) = hd_coeff(n+1)*window_coeff(n+1);
end

% now take Z-Transform
w = linspace(0,pi,1000);
H = zeros(1,1000);H

% fill H(i) for each w(i)
for i = 1:1000
    for n = 0:N-1
        H(i) = H(i) + h_coeff(n+1)*exp(-j*n*w(i)); 
    end
end


plot
subplot(2,1,1);
plot(w/pi, (abs(H))); % magnitude plot
title("Magnitude Plot");
subplot(2,1,2);
plot(w/pi, mag2db(abs(H))); % log magnitude plot
title("Log Magnitude Plot");



