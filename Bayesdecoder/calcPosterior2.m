function [Posterior, nonNormPosterior] = calcPosterior2(respModel, Y, Tbinsize, baseline, errmodel)
% Function to calculate the posterior probability of a stimulus set (Posterior), given the
% population firing rate (Y) and response model (respModel)
% Adapted from code from Daniel Bendor (Bendor & Wilson, Nat Neuro, 2013)
% 
% Usage: [Posterior] = calcPosterior(respModel,Y)
% 
% Aman Saleem
% October 2013
T       = size(Y,1);        %no.of time bins
numCells= size(Y,2);        %no.of Cells
M       = size(respModel,3);%no.of response model bins
if min(respModel(:)) < (max(respModel(:))/(100*M))
    respModel = respModel + (max(respModel(:))/(100*M));
end
% respModel(sum(isnan(respModel'))>1,:)=0;
ProdFields = log(ones(T,M));
sumFields = zeros(1,M);
ProdProb = ones(T,M);

Xprobmodel = zeros(M);
X = 1:M;
sigma = floor(M/10);
gauss = exp(-0.5*(X - floor(M/2)).^2/sigma^2);
gauss = gauss/sum(gauss);
for xx = 1:M
    Xprobmodel(xx,:) = circshift(gauss, xx-floor(M/2),2);
end

for icell = 1:numCells
    spkCount = max(0,round(Y(:,icell)))+1;  %round(Y(:,icell))+1;  
    spkCount(spkCount > size(respModel,2)) = size(respModel,2);
    ProdProb = ProdProb.*squeeze(respModel(icell,spkCount,:));
end
% [~, Xpred] = max(ProdProb(1,:),[],2);
% for tt = 2:T
%     ProdProb(tt,:) = ProdProb(tt,:).*Xprobmodel(Xpred,:);
%     [~, Xpred] = max(ProdProb(tt,:),[],2);
% end

Ratemod = ones(size(Tbinsize));
nonNormPosterior = ProdProb;
nonNormPosterior = nonNormPosterior;
Posterior = nonNormPosterior./(sum(nonNormPosterior,2)*ones(1,M));

Posterior = Posterior/baseline;
end