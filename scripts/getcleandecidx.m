function [goodidx, Xpred, Err] = getcleandecidx(post,X,th,maxtol,trialID)
    Prange = size(post,2);
    mattemp = post;
    
    [maxval, ~] = max(mattemp,[],2);
    maxvalmat = repmat(maxval,1,Prange);
    mattemp(mattemp < (1-maxtol)*maxvalmat) = 0;
    
    
    postfilt = zeros(Prange,1);
    postfilt(1:Prange) = (0:(Prange-1))'/(Prange)*(2*pi)-pi/2;
    a = sum(mattemp*cos(postfilt),2);b = sum(mattemp*sin(postfilt),2);
    pred0 = atan2(-a,b)+pi;
    pred0 = pred0/(2*pi)*Prange +1;
    
%     [~, pred0] = max(post,[],2);  
    Xpred = pred0;
    XX = X;
    predplus = pred0 + Prange;
    predminus = pred0 - Prange;
    pred0 = pred0.* (abs(XX-pred0)<abs(XX-predplus) & abs(XX-pred0)<abs(XX-predminus)) + predplus.* (abs(XX-predplus)<abs(XX-pred0) & abs(XX-predplus)<abs(XX-predminus)) + predminus.* (abs(XX-predminus)<abs(XX-pred0) & abs(XX-predminus)<abs(XX-predplus));
%     Xpred = pred0;
    
    Err = (pred0 - XX);%mod(abs(pred0 - XX),floor(Prange/2));
    goodidx = (mod(abs(Err),floor(Prange/2)) < th);%mod(abs(pred0 - XX),floor(Prange/2));% 
    
    goodidx = maxval >= th;
    
%     goodidx = maxval > th;
    
%     [maxpost,pred0] = max(post,[],2);
%     goodidx = maxpost > 2;
    
%     Xpred = pred0;
    
%     goodidx = (mod(min([abs(pred0 - XX) abs(predplus - XX) abs(predminus - XX)],[],2),Prange) < th);
    
%     trialnum = unique(trialID);
%     trialEV = zeros(1,numel(trialnum));
%     for itrial = 1:numel(trialnum)
%         idx = find(trialID == trialnum(itrial));
%         trialEV(itrial) = sum(mod(abs(pred0(idx) - XX(idx)),floor(Prange/2)).^2)/(numel(idx)*50^2);
%     end
%     goodtrials = find(trialEV <= 0.2);
%     goodidx = ismember(trialID,goodtrials);% & (mod(abs(pred0 - XX),floor(Prange/2)) < th);
    
%     th = 2;
%     [maxPost,~] = max(post,[],2);
%     goodidx = maxPost > th;
end