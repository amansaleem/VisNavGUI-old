function Xout = getBayesVar(X,latcorrection,sampleRate,Tsmth_win,smthtype,Fcircular,maxX)
X = circshift(X,[round(latcorrection/(1000/sampleRate)) 0]);
if Fcircular
    Xtemp = unwrap(X/maxX*2*pi)*maxX/(2*pi);
    Xtemp = smthInTime(Xtemp, sampleRate, Tsmth_win, 'same', [], smthtype);
    Xout = mod(Xtemp,maxX);
else
    Xout = smthInTime(X, sampleRate, Tsmth_win, 'same', [], smthtype);
end
end