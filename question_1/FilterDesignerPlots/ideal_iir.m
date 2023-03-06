Wp = [ 0.45 0.55];
Ws = [0.4 0.65];
Rp = 0.5;
Rs = 40;

[n,w] = ellipord(Wp,Ws,Rp,Rs)

[z,p,k] = ellip(n,Rp,Rs,Wp)
sos = zp2sos(z,p,k);

[h,w] =freqz(sos,'half',2001);

plot(w/pi,20*log10(abs(h)))
ax = gca;
ax.YLim = [-100 20];
ax.XTick = 0:.5:1;
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
