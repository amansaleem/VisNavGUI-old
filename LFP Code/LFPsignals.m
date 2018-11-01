function LFPsig = LFPsignals(input, sampRate, low, high)
% Usage: LFPsig = LFPsignals(input, sampRate, low, high)
% Get the instantaneous phase, power and frequency of the LFP recording in
% the low-high band range

LFPsig.sampRate = sampRate;
% if nargin>2
% Filter the LFP
LFPsig.filtSig = LFPfilter(input, low, high, sampRate);
LFPsig.band = [low high];
% end

hilTrans = hilbert(LFPsig.filtSig);

LFPsig.hill = (hilTrans);
% The instantaneous phase in the band
LFPsig.phase = angle(hilTrans);

% Inst. freq is the difference in the phase between two consecutive time
% pts.
LFPsig.instFreq = [0 diff(LFPsig.phase)*sampRate/(2*pi)];
LFPsig.instFreq( abs(LFPsig.instFreq) > high ) = NaN;
LFPsig.instFreq( abs(LFPsig.instFreq) < low ) = NaN;

% The power
LFPsig.power = abs(hilTrans);