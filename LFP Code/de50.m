function [sigOut, interference]=de50(sig, Fs, f0, NW)


if nargin<3
    % removing 50 Hz interference
    f0=50; % [Hz]
end

if nargin<4
    NW=2.5;
end

if nargin<2
    Fs=1000;
end

if nargin<1
    N=2500;
    t=(0:N-1)'/Fs;
    sig=sin(2*pi*f0*t+26)+0.2*randn(N,1);
else
    N=length(sig);
    t=(0:N-1)'/Fs;
end

df=Fs/N;
f=(0:N-1)*df;
fsymmetric=fftshift(f-Fs*(f>=Fs/2));
% plot(f, abs((fft(sig))))
% xlabel('f [Hz]')
% title('fft of EEG signal');

[v]=dpss(length(sig), NW);

sigRep=repmat(sig(:), 1, NW*2);

X=fft(v.*sigRep);
V=fft(v);

[minimum indFiftyHz]=min(abs(f-f0));
indZeroHz=1;

xi=X(indFiftyHz, :)';
vi=V(indZeroHz,:)';

C=inv(vi'*vi)*vi'*xi;
A=abs(C)*2;
phi=phase(C);

% this is statistical analysis
% Here, if Ff0 is larger than Fthreshold, then the sine interference is
% significantly different from noise
Ff0=(2*NW-1)*(abs(C))^2*(norm(vi))^2/sum((abs(xi-C*vi)).^2);
alpha=1/length(sig);
b=4*NW-2;
Fthreshold=b*(1-alpha^(2/b))/2/alpha^(2/b);


interference=A*cos(2*pi*f0*t-phi);
interference=reshape(interference, size(sig));

sigOut=sig-interference;

if nargout==0
    figure
    subplot(2,1,1)
plot(t, sig, 'b', t, interference, 'r:', t, sigOut, 'k', 'LineWidth', 2);
legend('original', 'interference', 'clean');
subplot(2,1,2)
semilogy(fsymmetric, abs(fftshift(fft(sig))).^2, 'b', fsymmetric, abs(fftshift(fft(interference))).^2, 'r:', fsymmetric, abs(fftshift(fft(sigOut))).^2, 'c:');
end
