function [obj, prediction, V, fPosterior, fnonNormPosterior] = trainBayesDecoder2D(obj, Vx, Vy, Z, T, Xsmth_win, Ysmth_win, subset, NiteRand)
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
% Oct 2013

if isempty(obj.kfold)
    obj.kfold = 5;
end

if isempty(obj.CVO)
    obj.CVO = crossValPartition(ones(1,length(Vx)),obj.kfold);
end

if obj.kfold == 1
    obj.CVO.kfold = 1;
    obj.CVO.train{1}= ones(1,length(Vx));
    obj.CVO.cv{1}   = obj.CVO.train{1};
    obj.CVO.test{1} = obj.CVO.train{1};
end

if isempty(obj.train_mean)
    calTrainMean = 1;
else
    calTrainMean = 0;
end

if sum(isnan(Vx))>0 || sum(isnan(Vy))>0 || sum(isnan(Z(:)))>0
    display('WARNING!!! Nans in the data: making a temp fix');
    t = ones(size(Vx));
    t(isnan(Vx)) = 0;
    t(isnan(Vy)) = 0;
    t(isnan(sum(Z,2))) = 0;
    Vx = Vx(t>0);
    Vy = Vy(t>0);
    Z = Z(t>0,:);
end
V = [Vx(:) Vy(:)];

obj.performance = zeros(1,obj.CVO.kfold);
Prediction = [];
% Bayes decoder needs discretization
% [V, obj.bins] = normalise1var(V, obj.numBins);
fPosterior = [];
fnonNormPosterior = [];

% Get the 1D maps
allMap = TtwoDimMap;
allMap.CVO = obj.CVO;
allMap.kfold = obj.kfold;
allMap.binsX  = obj.binsX; 
allMap.numBinsX = obj.numBinsX;
allMap.FcircularX = obj.Fcircular;
allMap.binsY  = obj.binsY; 
allMap.numBinsY = obj.numBinsY;
allMap.FcircularY = false;%obj.FcircularY;
allMap.FcomputeMarg = false;
allMap.FcomputePos = false;


allMap = allMap.trainSpikeMap(Vx, Vy, Z, T, Xsmth_win, Ysmth_win, obj.FoptiSmooth);

for iter = 1:obj.CVO.kfold
    
    Ztrain  = Z(obj.CVO.train{iter},:);
    Zcv     = Z(obj.CVO.cv{iter},:);
    Ztest   = Z(obj.CVO.test{iter},:);    

    Vxtrain  = Vx(obj.CVO.train{iter},:);
    Vxcv     = Vx(obj.CVO.cv{iter},:);
    Vxtest   = Vx(obj.CVO.test{iter},:);    
    
    Vytrain  = Vy(obj.CVO.train{iter},:);
    Vycv     = Vy(obj.CVO.cv{iter},:);
    Vytest   = Vy(obj.CVO.test{iter},:);    
    
    %% the main section of bayes decoder
    % Getting the 1D map for each neuron (calculating the place fields)
    % ...the main training component
    for icell = 1:size(Z,2)
%         [model(icell)] = get1Dmap(Y(:,icell), V', obj.numBins, obj.bins, obj.CVO, iter, obj.sampleRate, obj.smth_win);
        obj.model.trained(iter).respModel_orig(icell,:,:) = allMap.model.tuning(icell).respModel(iter,:,:);
        %         obj.model.trained(iter).respModel(icell,:) = model(icell).tuning./sum(model(icell).tuning);
        obj.model.trained(iter).respModel(icell,:,:) = allMap.model.tuning(icell).respModel(iter,:,:);
        obj.model.EV(iter,icell) = allMap.model.EV(iter,icell);
%         obj.model.L(iter,icell) = allMap.model.L(iter,icell);
%         obj.model.Q(iter,icell) = allMap.model.Q(iter,icell);
        tempModel(iter,icell,:,:) = obj.model.trained(iter).respModel(icell,:,:);
        minRate(iter,icell) = min(min(obj.model.trained(iter).respModel_orig(icell,:,:)));
        maxRate(iter,icell) = max(max(obj.model.trained(iter).respModel_orig(icell,:,:)));
    end
    obj.model.EV(iter,obj.model.EV(iter,:)<0) = 0;    
end

fPosterior = zeros(size(Z,1),obj.numBinsY,obj.numBinsX);
fnonNormPosterior = zeros(size(Z,1),obj.numBinsY,obj.numBinsX);
Prediction = zeros(size(Z,1),2);
for iter = 1:obj.CVO.kfold
    Ztrain  = Z(obj.CVO.train{iter},:);
    Zcv     = Z(obj.CVO.cv{iter},:);
    Ztest   = Z(obj.CVO.test{iter},:);%Y(obj.CVO.test{iter},:); 
    
    Ttest = T(obj.CVO.test{iter});

    Vxtrain  = Vx(obj.CVO.train{iter},:);
    Vxcv     = Vx(obj.CVO.cv{iter},:);
    Vxtest   = Vx(obj.CVO.test{iter},:);     
    
    Vytrain  = Vy(obj.CVO.train{iter},:);
    Vycv     = Vy(obj.CVO.cv{iter},:);
    Vytest   = Vy(obj.CVO.test{iter},:);     
    if calTrainMean
        train_meanx = mean(Vxtrain,1);
        train_meany = mean(Vytrain,1);
    else
        train_meanx = obj.train_mean;
    end
    goodCells = true(1,size(Z, 2));
    % Calculate the performance of the response model
    
    [Posterior, nonNormPosterior] = calcPosterior2D(obj.model.trained(iter).respModel(goodCells,:,:), Ztest(:,goodCells), Ttest, 1./(obj.numBinsX*obj.numBinsY)); % Get the posterior estimates
    
    [maxInTime] = max(max(Posterior,[],2),[],3);    
    % Relative performance
    for t = 1:length(Vxtest)
        probPeak(t,1) = Posterior(t,Vytest(t),Vxtest(t));
    end
%     obj.relPerformance(iter) = nanmean((maxInTime-probPeak)./(maxInTime-minInTime));
    obj.relPerformance(iter) = nanmean(probPeak);
    obj.confidence(iter)     = nanmean(maxInTime);
        
    [~, VxFit] = max(squeeze(sum(Posterior,2)),[],2);
    VxFit = VxFit';
    [~, VyFit] = max(squeeze(sum(Posterior,3)),[],2);
    VyFit = VyFit';
    
    % Normal performance
    Perfx(iter) = 1 - ((nansum((VxFit' - Vxtest).^2)))./((nansum((Vxtest - train_meanx).^2)));
    Perfy(iter) = 1 - ((nansum((VyFit' - Vytest).^2)))./((nansum((Vytest - train_meany).^2)));
    
    %
    obj.model.train_mean(iter) = train_meanx;
    %%
    fPosterior(obj.CVO.test{iter},:,:) = Posterior;    
    fnonNormPosterior(obj.CVO.test{iter},:,:) = nonNormPosterior;
    
    Prediction(obj.CVO.test{iter},:) = [VxFit(:) VyFit(:)];
end
Perfx(Perfx<0) = 0;
Perfy(Perfy<0) = 0;
obj.performance = [Perfx(:) Perfy(:)];
obj.meanPerformance = nanmean(obj.performance,1);

obj.model.meanModel = squeeze((nanmedian(tempModel,1)));
% obj.model.meanModel = reshape(obj.model.meanModel,size(obj.model.meanModel,2),size(obj.model.meanModel,3));

prediction  = Prediction;%Prediction';