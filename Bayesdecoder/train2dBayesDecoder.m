function [obj, prediction, V, fPosterior, fnPosterior] = trainBayesDecoder(obj, V, Y, smth_win, subset)
% program to implement a 1D linear decoder (using ridge regression) given training, cv and test
% sets , delay
% Usage: [Performance, Prediction, Model] = linearDecoder(Y, V, CVO, train_mean)
% Inputs: Y - firing rate
%         V - Variable being coded/decoded
%         CVO(optional) - crossvalidation object
%         train_mean(optional) - mean over the training set, special cases only
%
% Outputs:Performance - mean fraction of explained variance predicted by linear ridge decoder
%         Prediction  - prediction of the variable values by the decoder
%         Model       - The underlying model parameters (Useful for decoding novel stimuli)
%             WV - weight array (actual model)
%             Performance - same as performance above, but for individual iterations
%             train_mean  - needed for decoding (especially to calculate EV of decoded quantity)
%
% Aman Saleem
% Apr 2014

if isempty(obj.kfold)
    obj.kfold = 5;
end

if isempty(obj.CVO)
    obj.CVO = crossValPartition(ones(1,length(V)),obj.kfold);
end

if obj.kfold == 1
    obj.CVO.kfold = 1;
    obj.CVO.train{1}= ones(1,length(V));
    obj.CVO.cv{1}   = obj.CVO.train{1};
    obj.CVO.test{1} = obj.CVO.train{1};
end

if isempty(obj.train_mean)
    calTrainMean = 1;
else
    calTrainMean = 0;
end

if sum(isnan(V))>0 | sum(isnan(Y(:)))>0
    display('WARNING!!! Nans in the data: making a temp fix');
    t = ones(size(V));
    t(isnan(V)) = 0;
    t(isnan(sum(Y,2))) = 0;
    V = V(t>0);
    Y = Y(t>0,:);
end

for icell = 1:size(Y, 2);
    Y(:,icell) = smthInTime(Y(:,icell), obj.sampleRate, smth_win, 'same',subset);
end

obj.performance = zeros(1,obj.CVO.kfold);
Prediction = [];
fPosterior = [];
fnPosterior = [];

% Bayes decoder needs discretization
[V(:,1) obj.bins]   = normalise1var(V(:,1), obj.numBins);
[V(:,2) obj.binsB]  = normalise1var(V(:,2), obj.numBinsB);

pos_on_map = obj.numBins*(V(:,2)-1) + V(:,1);

% Get the 2D maps
allMap = twoDimMap;
allMap.CVO = obj.CVO;
allMap.kfold = obj.kfold;
allMap.bins  = obj.bins; 
allMap.numBins = obj.numBins;
allMap.numBinsB = obj.numBinsB;
allMap = allMap.trainSpikeMap(V, Y, smth_win);

for iter = 1:obj.CVO.kfold
    
    Ytrain  = Y(obj.CVO.train{iter},:);
    Ycv     = Y(obj.CVO.cv{iter},:);
    Ytest   = Y(obj.CVO.test{iter},:);
    
    Vtrain  = pos_on_map(obj.CVO.train{iter},:);
    Vcv     = pos_on_map(obj.CVO.cv{iter},:);
    Vtest   = pos_on_map(obj.CVO.test{iter},:);    

    if calTrainMean
        train_mean = mean(Vtrain(:),1);
    else
        train_mean = obj.train_mean;
    end
    
    %% the main section of bayes decoder
    % Getting the 1D map for each neuron (calculating the place fields)
    % ...the main training component
    respModel = zeros(size(Y,2),obj.numBins*obj.numBinsB);
    for icell = 1:size(Y,2)
        obj.model.trained(iter).respModel_orig(icell,:,:) = allMap.model.tuning(icell).respModel(iter,:,:);
        obj.model.trained(iter).respModel(icell,:,:) = allMap.model.tuning(icell).respModel(iter,:,:);
        obj.model.EV(iter,icell) = allMap.model.EV(iter,icell);
        tempModel(iter,icell,:) = reshape(obj.model.trained(iter).respModel(icell,:,:),[],1);
        minRate(iter,icell) = min(min(obj.model.trained(iter).respModel_orig(icell,:)));
        maxRate(iter,icell) = max(max(obj.model.trained(iter).respModel_orig(icell,:)));
    end
    obj.model.EV(iter,obj.model.EV(iter,:)<0) = 0;
    goodCells = find(obj.model.EV(iter,:)>=0);
%     goodCells = find(obj.model.EV(iter,:)>=0 & (maxRate(iter,:)>2*minRate(iter,:)));
    
    % Calculate the performance of the response model
    [Posterior] = calcPosterior(obj.model.trained(iter).respModel(goodCells,:), Ytest(:,goodCells)); % Get the posterior estimates
    baseline    = 1./(obj.numBins*obj.numBinsB);
    nPosterior = Posterior./repmat(median(Posterior),size(Posterior,1),1);
    nPosterior = nPosterior./(sum(nPosterior,2)*ones(1,size(nPosterior,2)));
   
    % Converting to log likelihood
    Posterior = log2(Posterior/baseline);
    nPosterior = log2(nPosterior/baseline);
    
    [maxInTime] = max(Posterior,[],2);
    [minInTime] = min(Posterior,[],2);
    % Relative performance
    for t = 1:length(Vtest)
        probPeak(t,1) = Posterior(t,Vtest(t));
    end
%     obj.relPerformance(iter) = nanmean((maxInTime-probPeak)./(maxInTime-minInTime));
    obj.relPerformance(iter) = (nanmean(probPeak));
    obj.confidence(iter)     = nanmean(maxInTime);
        
    
    [~, VFit] = max(Posterior,[],2);
    VFit = VFit';
    
    % Normal performance
    Perf(iter) = 0;%1 - ((nansum((VFit' - Vtest).^2)))./((nansum((Vtest - train_mean).^2)));
    
    %
    obj.model.train_mean(iter) = train_mean;
    %%
    fPosterior = [fPosterior' Posterior']';
    fnPosterior = [fnPosterior' nPosterior']';
    
    Prediction = [Prediction VFit];
end
Perf(Perf<0) = 0;
obj.performance = Perf;
obj.meanPerformance = nanmean(Perf);

obj.model.meanModel = squeeze(nanmedian(tempModel,1));

obj.model.bestModel = zeros(size(obj.model.meanModel));
[~,ibestPerf] = max(obj.relPerformance);

for icell = 1:size(Y,2)
        obj.model.bestModel(icell,:) = obj.model.trained(ibestPerf).respModel(icell,:);
end
prediction  = Prediction';