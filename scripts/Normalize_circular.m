function vec_out = Normalize_circular(vec_in,numBins)
m = mod(circ_mean(vec_in(:)/numBins*2*pi,[],1)+2*pi,2*pi)*numBins/(2*pi);
sd = circ_std(vec_in(:)/numBins*2*pi,[],[],1)*numBins/(2*pi);
if sd == 0
    sd = 1;
end
vec_out = (vec_in-m)/sd;
end