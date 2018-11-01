function out = MKnotch(in, chunkSize, Fs, f0)

if nargin<4
    f0 = 50;
end
if nargin<3
    Fs = 1000;
end
if nargin<2
    chunkSize = 3;
end

chunkBins = chunkSize*Fs;

out = zeros(size(in));
idx = 1;
while idx<length(in)-chunkBins
    out(idx:idx+chunkBins) = de50(in(idx:idx+chunkBins));
    idx = idx+chunkBins+1;
end
out(idx:end) = de50(in(idx:end));
