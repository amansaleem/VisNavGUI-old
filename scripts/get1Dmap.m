function [model, spred] = get1Dmap(zua, variable, scaling, nGrid, bins, CVO, iter, sampFreq, win, FoptiSmooth,Fcircular)
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
    win = 10;
    smth = 0;
else 
    smth = 1;
end

% if smth
%     zua_train = smthInTime(zua(CVO.train{iter}>0), sampFreq, win);
%     zua_test  = smthInTime(zua(CVO.test{iter}>0), sampFreq, win);
%     zua_cv    = smthInTime(zua(CVO.cv{iter}>0), sampFreq, win);
% else
    zua_train = zua(CVO.train{iter}>0).*scaling(CVO.train{iter}>0);
    zua_test  = zua(CVO.test{iter}>0).*scaling(CVO.test{iter}>0);
    zua_cv    = zua(CVO.cv{iter}>0).*scaling(CVO.cv{iter}>0);
% end
test = zua(CVO.test{iter});

% get spike count map
scMap = full(sparse(X_train, 1, (zua_train), nGrid, 1));
% get occupancy map
occMap = full(sparse(X_train, 1, 1, nGrid, 1));

if Fcircular
    scMap = [scMap;scMap;scMap];
    occMap = [occMap;occMap;occMap];
end
    
n1 = nGrid;
grids = 1:n1;

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
    idx = floor(numel(grids)/win);%/4);%25;%numel(grids);%100;    
end

model.swin = grids(idx);
model.bins = bins;

model.tuning = special_smooth_1d(scMap, 1./model.swin, bins, n1)...
    ./special_smooth_1d(occMap, 1./model.swin, bins, n1);
if Fcircular
    model.tuning = model.tuning(nGrid+1:2*nGrid);
end

% model.tuning = scMap ./ occMap;
% model.tuning(isnan(model.tuning)) = 0;
% model.tuning = special_smooth_1d(model.tuning, 1./model.swin, bins, n1);
    
pred  = (model.tuning(X_test));

if smth
%     spred    = smthInTime(pred, sampFreq, win, [], [], 'boxcar');
    spred = pred;
else
    spred = pred;
end

[model.EV model.corr model.L model.Q model.train_mean] = calCrossValExpVar(zua_train, zua_test, spred, test, pred);


P_i = 1./length(model.tuning); %occMap./sum(occMap(:));% 
R_mean = nanmean(model.tuning); % nanmean(X_train);
R_i = model.tuning; % scMap;

model.skaggs = sum( P_i .* (R_i/R_mean) .* log2(R_i/R_mean));

end