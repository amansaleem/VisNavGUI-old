function [f, pow, ampSpectrum, powSpectrum] = powerSpectrum(X, sampFreq, plot_flag, colour)

X(isnan(X)) = [];
X = single(X);
L = length(X);

samp_interval = 1./sampFreq;
if nargin<4
    colour = 'b';
end
if nargin<3
    plot_flag = 1;
end

NyLimit = (1/samp_interval)/2;

F = linspace(0,1,round(L/2))*NyLimit;

Y = fft(X)/L;

ampSpectrum = abs(Y(1:round(L/2)));
powerSpectrum = Y(1:round(L/2)).*conj(Y(1:round(L/2)));

if plot_flag
    figure;
    hold on;
    
    subplot(211)
    plot(F,ampSpectrum,'color', colour);
    title('Amplitude spectrum')
    xlabel('Frequency (Hz)')
    hold on;
    
    subplot(212)
    plot(F,powerSpectrum,'color', colour);
    title('Power spectrum')
    xlabel('Frequency (Hz)')
    hold on;

end

data = X;

WL = (1/samp_interval)*10;
nO = round(WL*0.08);
[pow,f] = pwelch(data,WL,nO,[],(1/samp_interval));

% range = find(f>0.5 & f<3000);
% 
% f = f(range);
% pow = pow(range);