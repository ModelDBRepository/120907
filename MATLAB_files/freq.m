function [freq power] = freq(dt,V,p)
% plot frequency power analysis for V at time t in ms
% dt = the time interval between points given in V
% V = the discrete time series we want to analyze
% p = 1 to plot, 0 to just return discrete frequency and power
% freq = frequencies analyzed
% power = power at corresponding frequencies in freq

Y=fft(V);
N=length(Y);
power=(Y(1:ceil(N/2)).*conj(Y(1:ceil(N/2))))/N;
freq=((1:ceil(N/2))/N)*(1000/dt); % factor of 1000 to translate
                                  % from 1/ms to Hz

if p==1
	plot(freq,power)
	set(gca,'FontSize',14)
	xlabel('frequency in Hz','FontSize',14)
	ylabel('power','FontSize',14)
end
