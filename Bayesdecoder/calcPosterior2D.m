function [Posterior, nonNormPosterior] = calcPosterior2D(respModel, Z, Tbinsize, baseline, errmodel)
% Function to calculate the posterior probability of a stimulus set (Posterior), given the
% population firing rate (Y) and response model (respModel)
% Adapted from code from Daniel Bendor (Bendor & Wilson, Nat Neuro, 2013)
% 
% Usage: [Posterior] = calcPosterior(respModel,Y)
% 
% Aman Saleem
% October 2013
T       = size(Z,1);        %no.of time bins
numCells= size(Z,2);        %no.of Cells
My       = size(respModel,2);
Mx       = size(respModel,3);
M = My*Mx;
if min(respModel(:)) < (max(respModel(:))/(100*M))
    respModel = respModel + (max(respModel(:))/(100*M));
end
ProdFields = log(ones(T,My,Mx));
sumFields = zeros(1,My,Mx);
for icell = 1:numCells
    if sum(sum(isnan(respModel(icell,:,:)),2),3)>1 % | max(respModel(icell,:)*60)<2
        continue
    end
    spkCount = Z(:,icell);
    ProdFields = ProdFields + log((repmat(respModel(icell,:,:),[T 1 1])).^(repmat(spkCount,[1 My Mx])));
    sumFields = sumFields + respModel(icell,:,:);
end
% sumFields  = nansum(respModel,1);
Ratemod = ones(size(Tbinsize));%sum(Y,2)/(sum(Y(:))/size(Y,1));%
tau = Tbinsize;%ones(T,1);%
nonNormPosterior = ProdFields + log((exp(-repmat(Ratemod.*tau,[1 My Mx]).*repmat(sumFields,[T 1 1]))));%ProdFields + log((exp(-(Ratemod.*tau)*sumFields)));
meannonnormPost = mean((nonNormPosterior(:)));
nonNormPosterior = exp(nonNormPosterior - meannonnormPost);
Posterior = nonNormPosterior./repmat(sum(sum(nonNormPosterior,2),3),[1 size(nonNormPosterior,2) size(nonNormPosterior,3)]);

Posterior = Posterior/baseline;
Posterior(sum(Z,2) < 1,:,:) = NaN;
end