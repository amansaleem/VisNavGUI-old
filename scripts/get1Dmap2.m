function [model, spred] = get1Dmap2(zua, variable, scaling, nGrid, bins, CVO, iter, sampFreq, Xsmth_win, FoptiSmooth, Fcircular)
% function to fit and test a 1D model, find the optimal smoothing
% function too
if nargin < 10
    FoptiSmooth = true;
end
if nargin < 11
    Fcircular = true;
end
X_train = variable(CVO.train{iter}>0);
X_test  = variable(CVO.test{iter}>0);
X_cv    = variable(CVO.cv{iter}>0);

% for icell = 1:numCells
if nargin<8
    sampFreq = 200;
end
if nargin < 9
    Xsmth_win = 10;
    smth = 0;
else 
    smth = 1;
end

% if smth
%     zua_train = smthInTime(zua(CVO.train{iter}>0), sampFreq, win);
%     zua_test  = smthInTime(zua(CVO.test{iter}>0), sampFreq, win);
%     zua_cv    = smthInTime(zua(CVO.cv{iter}>0), sampFreq, win);
% else
    zua_train0 = zua(CVO.train{iter}>0);
    zua_test0  = zua(CVO.test{iter}>0);
    zua_cv0    = zua(CVO.cv{iter}>0);
% end
test0 = zua(CVO.test{iter});

n1 = nGrid;
grids = 1:n1;

maxNspk = round(max([zua_train0; zua_test0; zua_cv0]));
nbXbins = max(X_train);
model.tuning = zeros(maxNspk+1,nbXbins);

for i = 0:maxNspk
    zua_train = double(round(zua_train0) == i);
    zua_cv = double(round(zua_cv0) == i);

    % get spike count map
    scMap = full(sparse(X_train, 1, (zua_train), nGrid, 1));
    % get occupancy map
    occMap = full(sparse(X_train, 1, 1, nGrid, 1));
    if Fcircular
        scMap = [scMap;scMap;scMap];
        occMap = [occMap;occMap;occMap];
    end
    
    if FoptiSmooth
        EV = zeros(1,length(grids));
        % special_smooth_1d(input, win, bins, nGrid)
        % get the FR map by smoothing
        parfor w=1:length(grids)
            FRMap = special_smooth_1d(scMap, 1./grids(w), bins, n1)...
                ./special_smooth_1d(occMap, 1./grids(w), bins, n1);
            if Fcircular
                FRMap = FRMap(nGrid+1:2*nGrid);
            end
            %     FRMap = scMap ./ occMap;
            %     FRMap(isnan(FRMap)) = 0;
            %     FRMap = special_smooth_1d(FRMap, 1./grids(w), bins, n1);
            
            pred  = (FRMap(X_cv));
            if smth
                %             spred    = smthInTime(pred, sampFreq, w);
                spred = pred;
                %try this next before trying constant window
                %             spred    = smthInTime(pred, sampFreq, win, [], [], 'boxcar');
            else
                spred = pred;
            end
            EV(w) = calCrossValExpVar(zua_train, zua_cv, spred);
        end
        EV = 100*EV;%round(100*EV);
        [~, idx] = max(EV);
    else
        idx = floor(numel(grids)/Xsmth_win);%numel(grids);%25;%numel(grids);%100;
    end
    
    model.swin = grids(idx);
    model.bins = bins;
    tuning = special_smooth_1d(scMap, 1./model.swin, bins, n1)...
        ./special_smooth_1d(occMap, 1./model.swin, bins, n1);
    
%     if sum(tuning) > 0
%      tuning = tuning./sum(tuning);
%     end
    if Fcircular
        model.tuning(i+1,:) = tuning(nGrid+1:2*nGrid);
    else
        model.tuning(i+1,:) = tuning;
    end
end

% for x = 1:size(model.tuning,2)
%     if sum(model.tuning(:,x)) > 0
%         model.tuning(:,x) = model.tuning(:,x)./sum(model.tuning(:,x));
%     end
% end


pred  = zeros(size(X_test))';%(model.tuning(X_test));

if smth
%     spred    = smthInTime(pred, sampFreq, win, [], [], 'boxcar');
    spred = pred;
else
    spred = pred;
end

[model.EV model.corr model.L model.Q model.train_mean] = calCrossValExpVar(zua_train0, zua_test0, spred, test0, pred);

P_i = 1./length(model.tuning); %occMap./sum(occMap(:));% 
R_mean = nanmean(model.tuning(:)); % nanmean(X_train);
R_i = model.tuning(1,:); % scMap;

model.skaggs = sum( P_i .* (R_i/R_mean) .* log2(R_i/R_mean));

end