%Script File for Low pass analog prototype filter:
% Wp,Ws = passband and stopband frequencies in rad/sec
% Ap,As = passband and stopband attenuations in dB
% N = filter order
% B,A = (L+1)x3 matrices of numerator and denominator coefficients, N =
% 2*L+r
%r = mod(N,2); can onlu be 0, or 1
% is 1 in case of odd numbered functions as no. of real valued poles in
% that case is 1 while in case of even chareceristic functions the no. of
% real valued poles is 0 as all of them are conjugate

%%%%  --------------Used functions:---------------------
%landen Landen transformation
%cde,acde cd elliptic function and its inverse
%sne,asne sn elliptic function and its inverse
%cne,dne cn and dn elliptic functions (for real arguments)
%ellipk complete elliptic integral K(k), 
%ellipdeg exact solution of degree equation (k from N, k1)



function [N,B,A,za,pa] = lpa(Wp, Ws, Ap, As)
%in case no. of arguments is less,returns error
if nargin==0, help lpa; return; end

%Ripple values based on attenuation
ep = sqrt(10^(Ap/10)-1); es = sqrt(10^(As/10)-1);
%k, k1 are known as the selectivity and discrimination parameters
% A narrow transition width would imply that k~< 1, whereas a deep stopbandor a flat passband would imply that k1<< 1. Thus, for most practical desired specifications, we
%will have k1 << k ~< 1.
k = Wp/Ws; k1 = ep/es;
%Finding the minimum order of the filer from the jacobian integrals of the
%selectivity and discrimination parameters and taking the largest integer
%close to N and finally modifying k to fit the new value of approximated N by uisng ellipdeg 
[K,Kp] = ellipk(k);
[K1,K1p] = ellipk(k1);
Nexact = (K1p/K1)/(Kp/K);
N = ceil(Nexact);
k = ellipdeg(N,k1); 

% The transfer function of an elliptic lowpass analog filter is constructed from its zeros and poles {zai, pai} in the second-order factored form:
% e L is the number of analog second-order sections, related to the filter order by N = 2L+r.
%Again, the notation [F]r means that the factor F is present if r = 1 and absent if r = 0. The
% quantity H0 is the gain at Ω = 0
r = mod(N,2); L = (N-r)/2; i = (1:L)'; u=(2*i-1)/N;

% finding the poles and zeros of jacobian polynomials
v0 = -1i * asne(1i/ep, k1) / N;
za = Wp * 1i ./(k*cde(u,k));
pa = Wp * 1i * cde(u-1i*v0, k);
pa0 = Wp * 1i * sne(1i*v0, k);

%Slicing elements of the coeffecients matrioces A and B to find the lowpass
%equivakent transfer function
B = [ones(L+1,1), zeros(L+1,2)];
A = [ones(L,1), -2*real(1./pa), abs(1./pa).^2];
A = [[1, -r*real(1/pa0), 0]; A];
B(2:L+1,:) = [ones(L,1), -2*real(1./za), abs(1./za).^2];
Gp = 10^(-Ap/20);
B(1,1) = Gp^(1-r);
end