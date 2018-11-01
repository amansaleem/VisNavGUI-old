function Xmodel = VisDistXmodel(Xactual,delta,alpha,latcorrection,es,Tsmth_win, wintype,numBins)
sampleRate = mean(1./es.sampleSize);

Xactual = circshift(Xactual,[round(latcorrection/(1000/sampleRate)) 0]);

Xtemp = unwrap(Xactual/max(Xactual)*2*pi)*max(Xactual)/(2*pi);
Xtemp = smthInTime(Xtemp, sampleRate, Tsmth_win, 'same', [], wintype);
Xtemp = mod(Xtemp,max(Xactual));

Nmemory = 1000;
Taumemory = delta;
if Taumemory > 0
    NTaumemory = round(Taumemory/(1000/sampleRate));
    Xmemory = [0 alpha.^((1:Nmemory)-1)];
    
    trialID = unique(es.trialID);
    Xvisual = zeros(size(Xtemp));
    Xvisual(1) = Xtemp(1);
    Rundist = zeros(size(Xtemp));
    maxX = max(Xtemp);
    for trial = 1:numel(trialID)
        trialidx = min(find(es.trialID == trialID(trial)),numel(Xtemp));
        Xtrial = Xtemp(trialidx);
        Xtrial = unwrap(Xtrial/maxX*2*pi)*maxX/(2*pi);
        gaintrial = es.gain(trialidx)/mode(es.gain);
        Xtemptrial = zeros(size(Xtrial));
        for k = 1:NTaumemory
%             Xtrial(k:NTaumemory:end) = unwrap(Xtrial(k:NTaumemory:end)/maxX*2*pi)*maxX/(2*pi);
            Xconv = conv(Xmemory,Xtrial(k:NTaumemory:end)-Xtrial(k),'full');
            Xtemptrial(k:NTaumemory:end) = Xtrial(k) + Xconv(1:numel(Xtrial(k:NTaumemory:end))) * (1-alpha);
        end
        Xvisual(trialidx) = Xtemptrial;
        
        Rundisttemptrial = zeros(size(Xtrial));
        for k = 1:NTaumemory
            Rundisttrial = [diff(Xtrial(k:NTaumemory:end));0]./(gaintrial(k:NTaumemory:end));
            Sconv = conv(Xmemory,Rundisttrial,'full');
            Rundisttemptrial(k:NTaumemory:end) = Sconv(1:numel(Rundisttrial));
        end
        Rundist(trialidx) = Rundisttemptrial;
    end
    Xmodel = mod(Xvisual + Rundist,maxX);
elseif Taumemory == 0
    Xmodel = Xtemp;
end
if numBins > 0
    [Xmodel, ~] = normalise1var(Xmodel, numBins);
end
end