function [filteredLFP] = LFPfilter(input, low, high, sampRate)
% Usage:     [filteredLFP] = LFPfilter(input, low, high, sampRate)
% Function to filter LFP with 'sampRate' using a elliptic filter in a band from 
% 'low' to 'high' 

[b,a] = ellip(2, 0.1, 40, [low high]/(sampRate/2));
input(isnan(input)) = 0;
filteredLFP = filtfilt(b, a, input);