classdef TtwoDimMap < TspikeMap
    %     Create a spikeMap class object of type '2D'. The decoder can be trained using the
    %     function 'trainDecoder' and used to decode with the function
    %     'predictVariable'
    %     [obj, prediction, X] = trainSpikeMap(obj, X, Y, Z, Xsmth_win, Ysmth_win);
    %     
    % Aman Saleem
    % Jan 2014
    properties
        kfold;       % number of partitions made to train the decoder
        train_mean;  % mean response used for training
        CVO;         % cross-validation structure
        binsX;        % the values of the bin limits
        binsY;        % the values of the bin limits
        numBinsX;     % number of bins by which the variable is discretized
        numBinsY;     % number of bins by which the variable is discretized
        FcircularX;   % if true, circular smoothing
        FcircularY;   % if true, circular smoothing
        qthreshold;
        FcomputePos;
        FcomputeMarg;
        FcomputeCorr
        Fdiscarditer
    end
    
    methods
        function obj = TtwoDimMap(varargin)
            pnames = {'dimensionality' 'variable' 'variable_range'...
                'Xsmth_win' 'Ysmth_win' 'model'...
                'kfold' 'performance'...
                'train_mean' 'CVO' 'sampleRate' 'XnumBins' 'FcircularX' 'YnumBins' 'FcircularY' 'qthreshold' 'FcomputePos' 'FcomputeMarg' 'FcomputeCorr' 'Fdiscarditer'};
            dflts  = {'2D' 'PxS' []...
                1 1 []...
                20 []...
                [] [] 60 100 true 18 false 1 true true true true};
            [obj.dimensionality, obj.variable, obj.variable_range, ...
                obj.Xsmth_win, obj.Ysmth_win,obj.model,...
                obj.kfold, obj.performance,...
                obj.train_mean, obj.CVO, obj.sampleRate, obj.numBinsX, obj.FcircularX, obj.numBinsY, obj.FcircularY, obj.qthreshold,...
                obj.FcomputePos, obj.FcomputeMarg, obj.FcomputeCorr, obj.Fdiscarditer] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        
        function [obj, Prediction, X] = trainSpikeMap(obj, X, Y, Z, T, Xsmth_win, Ysmth_win, FoptiSmooth)
            if nargin<6
                Xsmth_win = obj.Xsmth_win;
                Ysmth_win = obj.Ysmth_win;
            else
                obj.Xsmth_win = Xsmth_win;
                obj.Ysmth_win = Ysmth_win;
            end
            if nargin < 8
                FoptiSmooth = false;
            end
                        
            if isempty(obj.kfold)
                obj.kfold = 5;
            end
            
            if isempty(obj.train_mean)
                calTrainMean = 1;
            else
                calTrainMean = 0;
            end
            if ~FoptiSmooth
                disp(['fixed spatial window  = ' num2str(Xsmth_win) '%']);
            end
            
            if size(X,1) ~= size(Z,1)
                error('wrong dimension for X in TtwoDimMap');
            end
            if size(Y,1) ~= size(Z,1)
                error('wrong dimension for Y in TtwoDimMap');
            end
            
            if sum(isnan(X(:)))>0 || sum(isnan(Y(:)))>0
                display('WARNING!!! Nans in the data: making a temp fix');
                t = ones(size(X));
                t(isnan(sum(X,2))) = 0;
                t(isnan(sum(Y,2))) = 0;
                X = X(t>0);
                Y = Y(t>0);
                Z = Z(t>0,:);
            end
            
            % if sum(X<0)>0
            %     display('WARNING!!! -ve variable data present: ignorning them');
            %     t = ones(size(X));
            %     t(X<0) = 0;
            %     X = X(t>0);
            %     Y = Y(t>0,:);
            % end
            if isempty(obj.CVO)
                obj.CVO = crossValPartition(ones(1,size(X,1)),obj.kfold);
                %     obj.CVO = crossValPartition_MD(obj.CVO, ones(1,length(X)), 1:length(X), obj.kfold);
            end
            
            if obj.kfold == 1
                obj.CVO.kfold = 1;
                obj.CVO.train{1}= ones(1,size(X,1));
                obj.CVO.cv{1}   = obj.CVO.train{1};
                obj.CVO.test{1} = obj.CVO.train{1};
            end
            
            obj.performance = zeros(1,obj.CVO.kfold);
            Prediction = [];
            
            % Bayes decoder needs discretization
%             [X, obj.binsX] = normalise1var(X, obj.numBinsX);
%             [Y, obj.binsY] = normalise1var(Y, obj.numBinsY);
            
            Ncells = size(Z,2);
            pEV = NaN(obj.CVO.kfold,Ncells);
            pL = NaN(obj.CVO.kfold,Ncells);
            pQ = NaN(obj.CVO.kfold,Ncells);
            pskaggs = NaN(obj.CVO.kfold,Ncells);
            ptrain_mean = NaN(obj.CVO.kfold,Ncells);
            pswinX = NaN(obj.CVO.kfold,Ncells);
            pswinY = NaN(obj.CVO.kfold,Ncells);
            for icell = 1:Ncells
                Tuning(icell).respModel = zeros(obj.CVO.kfold,obj.numBinsY,obj.numBinsX);
                Tuning(icell).meanrespModel = [];
                Tuning(icell).SErespModel = [];
                if obj.FcomputeMarg
                    Tuning(icell).respModelX = zeros(obj.CVO.kfold,obj.numBinsX);
                    Tuning(icell).respModelY = zeros(obj.CVO.kfold,obj.numBinsY);
                    Tuning(icell).respModelXNorm = zeros(obj.CVO.kfold,obj.numBinsX);
                    Tuning(icell).respModelYNorm = zeros(obj.CVO.kfold,obj.numBinsY);
                    
                    Tuning(icell).meanrespModelX = [];
                    Tuning(icell).meanrespModelY = [];
                    Tuning(icell).meanrespModelXNorm = [];
                    Tuning(icell).meanrespModelYNorm = [];
                    
                    Tuning(icell).SErespModelX = [];
                    Tuning(icell).SErespModelY = [];
                    Tuning(icell).SErespModelXNorm = [];
                    Tuning(icell).SErespModelYNorm = [];
                end
                if obj.FcomputeCorr
                    Tuning(icell).corrModelX = zeros(obj.CVO.kfold,obj.numBinsY,obj.numBinsX);
                    Tuning(icell).corrModelY = zeros(obj.CVO.kfold,obj.numBinsY,obj.numBinsX);
                    Tuning(icell).corrModelXpos = zeros(obj.CVO.kfold,obj.numBinsY);
                    Tuning(icell).corrModelXmax = zeros(obj.CVO.kfold,obj.numBinsY);
                    Tuning(icell).corrModelYpos = zeros(obj.CVO.kfold,obj.numBinsX);
                    Tuning(icell).corrModelYmax = zeros(obj.CVO.kfold,obj.numBinsX);
                    
                    Tuning(icell).corrModelXposNorm = zeros(obj.CVO.kfold,obj.numBinsY);
                    Tuning(icell).corrModelXmaxNorm = zeros(obj.CVO.kfold,obj.numBinsY);
                    Tuning(icell).corrModelYposNorm = zeros(obj.CVO.kfold,obj.numBinsX);
                    Tuning(icell).corrModelYmaxNorm = zeros(obj.CVO.kfold,obj.numBinsX);
                    
                    Tuning(icell).meancorrModelX = [];
                    Tuning(icell).meancorrModelY = [];
                    Tuning(icell).meancorrModelXpos = [];
                    Tuning(icell).meancorrModelXmax = [];
                    Tuning(icell).meancorrModelYpos = [];
                    Tuning(icell).meancorrModelYmax = [];
                    Tuning(icell).meancorrModelXposNorm = [];
                    Tuning(icell).meancorrModelXmaxNorm = [];
                    Tuning(icell).meancorrModelYposNorm = [];
                    Tuning(icell).meancorrModelYmaxNorm = [];
                    
                    Tuning(icell).SEcorrModelX = [];
                    Tuning(icell).SEcorrModelY = [];
                    Tuning(icell).SEcorrModelXpos = [];
                    Tuning(icell).SEcorrModelXmax = [];
                    Tuning(icell).SEcorrModelYpos = [];
                    Tuning(icell).SEcorrModelYmax = [];
                    Tuning(icell).SEcorrModelXposNorm = [];
                    Tuning(icell).SEcorrModelXmaxNorm = [];
                    Tuning(icell).SEcorrModelYposNorm = [];
                    Tuning(icell).SEcorrModelYmaxNorm = [];
                end
                if obj.FcomputePos
                    Tuning(icell).respModelXpos = zeros(obj.CVO.kfold,obj.numBinsY);
                    Tuning(icell).respModelYpos = zeros(obj.CVO.kfold,obj.numBinsX);
                    Tuning(icell).respModelXmax = zeros(obj.CVO.kfold,obj.numBinsY);
                    Tuning(icell).respModelYmax = zeros(obj.CVO.kfold,obj.numBinsX);
                    
                    Tuning(icell).respModelXposNorm = zeros(obj.CVO.kfold,obj.numBinsY);
                    Tuning(icell).respModelYposNorm = zeros(obj.CVO.kfold,obj.numBinsX);
                    Tuning(icell).respModelXmaxNorm = zeros(obj.CVO.kfold,obj.numBinsY);
                    Tuning(icell).respModelYmaxNorm = zeros(obj.CVO.kfold,obj.numBinsX);
                    
                    Tuning(icell).meanrespModelXpos = [];
                    Tuning(icell).meanrespModelYpos = [];
                    Tuning(icell).meanrespModelXmax = [];
                    Tuning(icell).meanrespModelYmax = [];
                    Tuning(icell).meanrespModelXposNorm = [];
                    Tuning(icell).meanrespModelYposNorm = [];
                    Tuning(icell).meanrespModelXmaxNorm = [];
                    Tuning(icell).meanrespModelYmaxNorm = [];
                    
                    Tuning(icell).SErespModelXpos = [];
                    Tuning(icell).SErespModelYpos = [];
                    Tuning(icell).SErespModelXmax = [];
                    Tuning(icell).SErespModelYmax = [];
                    Tuning(icell).SErespModelXposNorm = [];
                    Tuning(icell).SErespModelYposNorm = [];
                    Tuning(icell).SErespModelXmaxNorm = [];
                    Tuning(icell).SErespModelYmaxNorm = [];
                    
                    Tuning(icell).fieldidxX = [];
                    Tuning(icell).fieldidxY = [];
                    
                    Tuning(icell).respModelslopeXY = [];
                    Tuning(icell).respModelphi0XY = [];
                    Tuning(icell).respModelrhoXY = [];
                    Tuning(icell).meanrespModelslopeXY = [];
                    Tuning(icell).meanrespModelphi0XY = [];
                    Tuning(icell).meanrespModelrhoXY = [];
                    Tuning(icell).SErespModelslopeXY = [];
                    Tuning(icell).SErespModelphi0XY = [];
                    Tuning(icell).SErespModelrhoXY = [];
                end
            end
            
            
            pnumBinsX = obj.numBinsX;
            pbinsX = obj.binsX;
            pnumBinsY = obj.numBinsY;
            pbinsY = obj.binsY;
            pCVO = obj.CVO;
            pkfold = obj.CVO.kfold;
            psampleRate = obj.sampleRate;
            pXsmth_win = obj.Xsmth_win;
            pYsmth_win = obj.Ysmth_win;
            pFcircularX = obj.FcircularX;
            pFcircularY = obj.FcircularY;
            pFcomputeMarg = obj.FcomputeMarg;
            pFcomputePos = obj.FcomputePos;
            pFcomputeCorr = obj.FcomputeCorr;
            pqthreshold = obj.qthreshold;
            pFdiscarditer = obj.Fdiscarditer;
            
            for iter = 1:obj.CVO.kfold        
                parfor icell = 1:Ncells
                    [model, ~] = get2Dmap(Z(:,icell), X(:,min(end,icell))', Y(:,min(end,icell))', 1./T, pnumBinsX, pbinsX, pnumBinsY, pbinsY, pCVO, iter, psampleRate, pXsmth_win, pYsmth_win, FoptiSmooth, pFcircularX, pFcircularY);
                    Tuning(icell).respModel(iter,:,:) = model.tuning;
                    if pFcomputeMarg
                        Tuning(icell).respModelX(iter,:) = nanmean(model.tuning,1);
                        Tuning(icell).respModelY(iter,:) = nanmean(model.tuning,2);
                        Tuning(icell).respModelXNorm(iter,:) = (Tuning(icell).respModelX(iter,:)-nanmean(Tuning(icell).respModelX(iter,:)))/nanstd(Tuning(icell).respModelX(iter,:));
                        Tuning(icell).respModelYNorm(iter,:) = (Tuning(icell).respModelY(iter,:)-nanmean(Tuning(icell).respModelY(iter,:)))/nanstd(Tuning(icell).respModelY(iter,:));
                    end
                    pEV(iter,icell) = model.EV;
                    pL(iter,icell) = model.L;
                    pQ(iter,icell) = model.Q;
                    pskaggs(iter,icell) = model.skaggs;
                    ptrain_mean(iter,icell) = model.train_mean;
                    pswinX(iter,icell) = model.swinX;
                    pswinY(iter,icell) = model.swinY;
                end
            end
            
            
            Xmesh = repmat(1:obj.numBinsX,[obj.numBinsY 1]);
            Ymesh = repmat((1:obj.numBinsY)',[1 obj.numBinsX]);
            Xmesh_centered = repmat(-floor(obj.numBinsX/2)+1:floor(obj.numBinsX/2),[obj.numBinsY 1]);
            Ymesh_centered = repmat((-floor(obj.numBinsY/2)+1:floor(obj.numBinsY/2))',[1 obj.numBinsX]);
            outcorrXrange = [1:floor(obj.numBinsX/4) floor(obj.numBinsX*3/4)+1:obj.numBinsX];
            
            parfor icell = 1:Ncells
                Tuning(icell).meanrespModel = squeeze(mean(Tuning(icell).respModel,1));
                if pFcomputeMarg || pFcomputePos
                    Tuning(icell).meanrespModelX = nanmean(squeeze(Tuning(icell).meanrespModel),1);
                    Tuning(icell).meanrespModelY = nanmean(squeeze(Tuning(icell).meanrespModel),2);
                    Tuning(icell).meanrespModelXNorm = (Tuning(icell).meanrespModelX-nanmean(Tuning(icell).meanrespModelX))/nanstd(Tuning(icell).meanrespModelX);
                    Tuning(icell).meanrespModelYNorm = (Tuning(icell).meanrespModelY-nanmean(Tuning(icell).meanrespModelY))/nanstd(Tuning(icell).meanrespModelY);
                end
                if pFcomputeCorr
                    Tuning(icell).meancorrModelX = zeros(pnumBinsY,pnumBinsX);
                    Xmap_norm = Tuning(icell).meanrespModelX;
                    Xmap_norm = Xmap_norm - nanmean(Xmap_norm);
                    Xmap_norm = Xmap_norm./(nanmean(Xmap_norm.^2)).^0.5;
                    map_normX = Tuning(icell).meanrespModel;
                    map_normX = (map_normX - repmat(nanmean(map_normX,2),[1,size(map_normX,2)]));
                    map_normX = map_normX./(repmat(nanmean(map_normX.^2,2),[1,size(map_normX,2)])).^0.5;
                    map_normX(:,[1:floor(pnumBinsX/4) floor(pnumBinsX*3/4)+1:pnumBinsX]);
                    ishift = 0;
                    for xshift = -floor(pnumBinsX/2)+1:floor(pnumBinsX/2)
                        ishift = ishift + 1;
                        Tuning(icell).meancorrModelX(:,ishift) = (map_normX*circshift(Xmap_norm,xshift,2)')/pnumBinsX;
                    end
                    corrmap = Tuning(icell).meancorrModelX;
                    corrmap(:,outcorrXrange) = 0;
                    if pFcircularX
                        Tuning(icell).meancorrModelXpos = getCircularAverage(corrmap',0,1);
                        Tuning(icell).meancorrModelXmax = getCircularAverage(corrmap',0,0.01,0.05);
                        Tuning(icell).meancorrModelXposNorm = Normalize_circular(Tuning(icell).meancorrModelXpos,pnumBinsX);
                        Tuning(icell).meancorrModelXmaxNorm = Normalize_circular(Tuning(icell).meancorrModelXmax,pnumBinsX);
                    else
                        Tuning(icell).meancorrModelXpos = sum(corrmap.*Xmesh_centered,2)./sum(corrmap,2);
                        [~, Tuning(icell).meancorrModelXmax] = max(corrmap,[],2);
                        Tuning(icell).meancorrModelXposNorm = Normalize(Tuning(icell).meancorrModelXpos);
                        Tuning(icell).meancorrModelXmaxNorm = Normalize(Tuning(icell).meancorrModelXmax);
                    end
                    
                    Tuning(icell).meancorrModelY = zeros(pnumBinsY,pnumBinsX);
                    Ymap_norm = Tuning(icell).meanrespModelY;
                    Ymap_norm = Ymap_norm - nanmean(Ymap_norm);
                    Ymap_norm = Ymap_norm./(nanmean(Ymap_norm.^2)).^0.5;
                    map_normY = Tuning(icell).meanrespModel;
                    map_normY = (map_normY - repmat(nanmean(map_normY,1),[size(map_normY,1) 1]));
                    map_normY = map_normY./(repmat(nanmean(map_normY.^2,1),[size(map_normY,1) 1])).^0.5;
                    ishift = 0;
                    for yshift = -floor(pnumBinsY/2)+1:floor(pnumBinsY/2)
                        ishift = ishift + 1;
                        Tuning(icell).meancorrModelY(ishift,:) = (map_normY'*circshift(Ymap_norm,yshift,1))/pnumBinsY;
                    end
                    corrmap = Tuning(icell).meancorrModelY;
                    corrmap(:,outcorrXrange) = 0;
                    if pFcircularY
                        Tuning(icell).meancorrModelYpos = getCircularAverage(corrmap,0,1);
                        Tuning(icell).meancorrModelYmax = getCircularAverage(corrmap,0,0.01,0.05);
                        Tuning(icell).meancorrModelYposNorm = Normalize_circular(Tuning(icell).meancorrModelYpos,pnumBinsY);
                        Tuning(icell).meancorrModelYmaxNorm = Normalize_circular(Tuning(icell).meancorrModelYmax,pnumBinsY);
                    else
                        Tuning(icell).meancorrModelYpos = sum(corrmap.*Ymesh_centered,1)./sum(corrmap,1);
                        [~, Tuning(icell).meancorrModelYmax] = max(corrmap,[],1);
                        Tuning(icell).meancorrModelYposNorm = Normalize(Tuning(icell).meancorrModelYpos);
                        Tuning(icell).meancorrModelYmaxNorm = Normalize(Tuning(icell).meancorrModelYmax);
                    end
                    
                    for iter = 1:pkfold
                        Xmap_norm = Tuning(icell).respModelX(iter,:);
                        Xmap_norm = Xmap_norm - nanmean(Xmap_norm);
                        Xmap_norm = Xmap_norm./(nanmean(Xmap_norm.^2)).^0.5;
                        map_normX = squeeze(Tuning(icell).respModel(iter,:,:));
                        map_normX = (map_normX - repmat(nanmean(map_normX,2),[1,size(map_normX,2)]));
                        map_normX = map_normX./(repmat(nanmean(map_normX.^2,2),[1,size(map_normX,2)])).^0.5;
                        ishift = 0;
                        for xshift = -floor(pnumBinsX/2)+1:floor(pnumBinsX/2)
                            ishift = ishift + 1;
                            Tuning(icell).corrModelX(iter,:,ishift) = (map_normX*circshift(Xmap_norm,xshift,2)')/pnumBinsX;
                        end
                        corrmap = Tuning(icell).corrModelX(iter,:,:);
                        corrmap(1,:,outcorrXrange) = 0;
                        if pFcircularX
                            Tuning(icell).corrModelXpos(iter,:) = getCircularAverage(squeeze(corrmap)',0,1);
                            Tuning(icell).corrModelXmax(iter,:) = getCircularAverage(squeeze(corrmap)',0,0.01,0.05);
                            Tuning(icell).corrModelXposNorm(iter,:) = Normalize_circular(Tuning(icell).corrModelXpos(iter,:),pnumBinsX);
                            Tuning(icell).corrModelXmaxNorm(iter,:) = Normalize_circular(Tuning(icell).corrModelXmax(iter,:),pnumBinsX);
                        else
                            Tuning(icell).corrModelXpos(iter,:) = sum(squeeze(corrmap).*Xmesh_centered,2)./sum(squeeze(corrmap),2);
                            [~, Tuning(icell).corrModelXmax(iter,:)] = max(squeeze(corrmap),[],2);
                            Tuning(icell).corrModelXposNorm(iter,:) = Normalize(Tuning(icell).corrModelXpos(iter,:));
                            Tuning(icell).corrModelXmaxNorm(iter,:) = Normalize(Tuning(icell).corrModelXmax(iter,:));
                        end
                        
                        Ymap_norm = (Tuning(icell).respModelY(iter,:))';
                        Ymap_norm = Ymap_norm - nanmean(Ymap_norm);
                        Ymap_norm = Ymap_norm./(nanmean(Ymap_norm.^2)).^0.5;
                        map_normY = squeeze(Tuning(icell).respModel(iter,:,:));
                        map_normY = (map_normY - repmat(nanmean(map_normY,2),[1,size(map_normY,2)]));
                        map_normY = map_normY./(repmat(nanmean(map_normY.^2,2),[1,size(map_normY,2)])).^0.5;
                        ishift = 0;
                        for yshift = -floor(pnumBinsY/2)+1:floor(pnumBinsY/2)
                            ishift = ishift + 1;
                            Tuning(icell).corrModelY(iter,ishift,:) = (map_normY'*circshift(Ymap_norm,yshift,1))/pnumBinsX;
                        end
                        corrmap = Tuning(icell).corrModelY(iter,:,:);
                        corrmap(1,:,outcorrXrange) = 0;
                        if pFcircularY
                            Tuning(icell).corrModelYpos(iter,:) = getCircularAverage(squeeze(corrmap),0,1);
                            Tuning(icell).corrModelYmax(iter,:) = getCircularAverage(squeeze(corrmap),0,0.01,0.05);
                            Tuning(icell).corrModelYposNorm(iter,:) = Normalize_circular(Tuning(icell).corrModelYpos(iter,:),pnumBinsY);
                            Tuning(icell).corrModelYmaxNorm(iter,:) = Normalize_circular(Tuning(icell).corrModelYmax(iter,:),pnumBinsY);
                        else
                            Tuning(icell).corrModelYpos(iter,:) = sum(squeeze(corrmap).*Ymesh_centered,1)./sum(squeeze(corrmap),1);
                            [~, Tuning(icell).corrModelYmax(iter,:)] = max(squeeze(corrmap),[],1);
                            Tuning(icell).corrModelYposNorm(iter,:) = Normalize(Tuning(icell).corrModelYpos(iter,:));
                            Tuning(icell).corrModelYmaxNorm(iter,:) = Normalize(Tuning(icell).corrModelYmax(iter,:));
                        end
                        
                    end
                end
                if pFcomputePos
                    [q_idxX,fieldX] = findfield(Tuning(icell).meanrespModelX,pqthreshold);
                    Tuning(icell).fieldidxX = fieldX;
                    outfieldX = find(~ismember(1:pnumBinsX,fieldX));
                    [q_idxY,fieldY] = findfield(Tuning(icell).meanrespModelY,pqthreshold);
                    Tuning(icell).fieldidxY = fieldY;
                    outfieldY = find(~ismember(1:pnumBinsY,fieldY));
                    [~,fieldXmax] = max(Tuning(icell).meanrespModelX);
                    
                    tuningX0 = squeeze(Tuning(icell).meanrespModel);
                    tuningX = tuningX0;
                    if ~isempty(fieldX)
                        tuningX(:,outfieldX) = 0;%repmat(min(tuningX(:,fieldX),[],2),[1 numel(outfieldX)]);
                    end
                    if pFcircularX
                        Tuning(icell).meanrespModelXpos = getCircularAverage(tuningX',0,1);
                        Tuning(icell).meanrespModelXmax = getCircularAverage(tuningX0',0,0.01,0.05);
                        Tuning(icell).meanrespModelXposNorm = Normalize_circular(Tuning(icell).meanrespModelXpos,pnumBinsX);
                        Tuning(icell).meanrespModelXmaxNorm = Normalize_circular(Tuning(icell).meanrespModelXmax,pnumBinsX);
                    else
                        Tuning(icell).meanrespModelXpos = sum(tuningX.*Xmesh,2)./sum(tuningX,2);
                        [~, Tuning(icell).meanrespModelXmax] = max(tuningX0,[],2);
                        Tuning(icell).meanrespModelXposNorm = Normalize(Tuning(icell).meanrespModelXpos);
                        Tuning(icell).meanrespModelXmaxNorm = Normalize(Tuning(icell).meanrespModelXmax);
                    end
                    
                    tuningY0 = squeeze(Tuning(icell).meanrespModel);
                    tuningY = tuningY0;
                    if ~isempty(fieldX)
                        tuningY(:,outfieldX) = 0;%repmat(min(tuningY(fieldY,:),[],1),[numel(outfieldY) 1]);
                    end
                    if pFcircularY
                        Tuning(icell).meanrespModelYpos = getCircularAverage(tuningY,0,1);
                        Tuning(icell).meanrespModelYmax = getCircularAverage(tuningY0,0,0.01,0.05);
                        Tuning(icell).meanrespModelYposNorm = Normalize_circular(Tuning(icell).meanrespModelYpos,pnumBinsY);
                        Tuning(icell).meanrespModelYmaxNorm = Normalize_circular(Tuning(icell).meanrespModelYmax,pnumBinsY);
                        
                        tuningY_centered = circshift(tuningY,-fieldXmax+floor(pnumBinsX/2),2);
                        [slope,phi0,rho] = circularlinearfit(Ymesh(:)/pnumBinsY*2*pi,Xmesh_centered(:)/numel(fieldX),tuningY_centered(:));
                        Tuning(icell).meanrespModelslopeXY = slope;
                        Tuning(icell).meanrespModelphi0XY = phi0*pnumBinsY/(2*pi);
                        Tuning(icell).meanrespModelrhoXY = rho;                        
                    else
                        Tuning(icell).meanrespModelYpos = sum(tuningY.*Ymesh,1)./sum(tuningY,1);
                        [~, Tuning(icell).meanrespModelYmax] = max(tuningY0,[],1);
                        Tuning(icell).meanrespModelYposNorm = Normalize(Tuning(icell).meanrespModelYpos);
                        Tuning(icell).meanrespModelYmaxNorm = Normalize(Tuning(icell).meanrespModelYmax);
                    end
                    
                    for iter = 1:pkfold
                        tuningX0 = squeeze(Tuning(icell).respModel(iter,:,:));
                        tuningX = tuningX0;
                        if ~isempty(fieldX)
                            tuningX(:,outfieldX) = 0;%repmat(min(tuningX(:,fieldX),[],2),[1 numel(outfieldX)]);
                        end
                        
                        if pFcircularX
                            Tuning(icell).respModelXpos(iter,:) = getCircularAverage(tuningX',0,1);
                            Tuning(icell).respModelXmax(iter,:) = getCircularAverage(tuningX0',0,0.01,0.05);
                            Tuning(icell).respModelXposNorm(iter,:) = Normalize_circular(Tuning(icell).respModelXpos(iter,:),pnumBinsX);
                            Tuning(icell).respModelXmaxNorm(iter,:) = Normalize_circular(Tuning(icell).respModelXmax(iter,:),pnumBinsX);
                        else
                            Tuning(icell).respModelXpos(iter,:) = sum(tuningX.*Xmesh,1)./sum(tuningX,2);
                            [~, Tuning(icell).respModelXmax(iter,:)] = max(tuningX0,[],2);
                            Tuning(icell).respModelXposNorm(iter,:) = Normalize(Tuning(icell).respModelXpos(iter,:));
                            Tuning(icell).respModelXmaxNorm(iter,:) = Normalize(Tuning(icell).respModelXmax(iter,:));
                        end
                        
                        tuningY0 = squeeze(Tuning(icell).respModel(iter,:,:));
                        tuningY = tuningY0;
                        if ~isempty(fieldX)
                            tuningY(:,outfieldX) = 0;%repmat(min(tuningY(fieldY,:),[],1),[numel(outfieldY) 1]);
                        end
                        if pFcircularY
                            Tuning(icell).respModelYpos(iter,:) = getCircularAverage(tuningY,0,1);
                            Tuning(icell).respModelYmax(iter,:) = getCircularAverage(tuningY0,0,0.01,0.05);
                            Tuning(icell).respModelYposNorm(iter,:) = Normalize_circular(Tuning(icell).respModelYpos(iter,:),pnumBinsY);
                            Tuning(icell).respModelYmaxNorm(iter,:) = Normalize_circular(Tuning(icell).respModelYmax(iter,:),pnumBinsY);
                            
                            tuningY_centered = circshift(tuningY,-fieldXmax+floor(pnumBinsX/2),2);
                            [slope,phi0,rho] = circularlinearfit(Ymesh(:)/pnumBinsY*2*pi,Xmesh_centered(:)/numel(fieldX),tuningY_centered(:));
                            Tuning(icell).respModelslopeXY(iter) = slope;
                            Tuning(icell).respModelphi0XY(iter) = phi0*pnumBinsY/(2*pi);
                            Tuning(icell).respModelrhoXY(iter) = rho;
                        else
                            Tuning(icell).respModelYpos(iter,:) = sum(tuningY.*Ymesh,1)./sum(tuningY,1);
                            [~, Tuning(icell).respModelYmax(iter,:)] = max(tuningY0,[],1);
                            Tuning(icell).respModelYposNorm(iter,:) = Normalize(Tuning(icell).respModelYpos(iter,:));
                            Tuning(icell).respModelYmaxNorm(iter,:) = Normalize(Tuning(icell).respModelYmax(iter,:));
                        end
                    end
                end
                stdresp = 0;
                stdrespX = 0;
                stdrespY = 0;
                stdrespXNorm = 0;
                stdrespYNorm = 0;
                stdrespXpos = 0;
                stdrespXmax = 0;
                stdrespYpos = 0;
                stdrespYmax = 0;
                stdrespXposNorm = 0;
                stdrespXmaxNorm = 0;
                stdrespYposNorm = 0;
                stdrespYmaxNorm = 0;
                stdslopeXY = 0;
                stdphi0XY = 0;
                stdrhoXY = 0;
                stdcorrX = 0;
                stdcorrXpos = 0;
                stdcorrXmax = 0;
                stdcorrY = 0;
                stdcorrYpos = 0;
                stdcorrYmax = 0;
                stdcorrXposNorm = 0;
                stdcorrXmaxNorm = 0;
                stdcorrYposNorm = 0;
                stdcorrYmaxNorm = 0;
                for i = 1:pkfold
                    r_iter = squeeze(Tuning(icell).respModel(i,:,:));
                    r_ave = Tuning(icell).meanrespModel;
                    stdresp = stdresp + (pkfold - 1)/pkfold*((r_iter - mean(r_iter(:))) - (r_ave - mean(r_ave(:)))).^2;
                    if pFcomputeMarg
                        r_iter = squeeze(Tuning(icell).respModelX(i,:));
                        r_ave = Tuning(icell).meanrespModelX;
                        r_iter = r_iter(:);
                        r_ave = r_ave(:);
                        stdrespX = stdrespX + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                        r_iter = squeeze(Tuning(icell).respModelY(i,:));
                        r_ave = Tuning(icell).meanrespModelY;
                        r_iter = r_iter(:);
                        r_ave = r_ave(:);
                        stdrespY = stdrespY + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                        
                        r_iter = squeeze(Tuning(icell).respModelXNorm(i,:));
                        r_ave = Tuning(icell).meanrespModelXNorm;
                        r_iter = r_iter(:);
                        r_ave = r_ave(:);
                        stdrespXNorm = stdrespXNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                        r_iter = squeeze(Tuning(icell).respModelYNorm(i,:));
                        r_ave = Tuning(icell).meanrespModelYNorm;
                        r_iter = r_iter(:);
                        r_ave = r_ave(:);
                        stdrespYNorm = stdrespYNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                    end
                    if pFcomputeCorr
                        r_iter = squeeze(Tuning(icell).corrModelX(i,:,:));
                        r_ave = Tuning(icell).meancorrModelX;
                        stdcorrX = stdcorrX + (pkfold - 1)/pkfold*((r_iter - mean(r_iter(:))) - (r_ave - mean(r_ave(:)))).^2;
                        r_iter = squeeze(Tuning(icell).corrModelY(i,:,:));
                        r_ave = Tuning(icell).meancorrModelY;
                        stdcorrY = stdcorrY + (pkfold - 1)/pkfold*((r_iter - mean(r_iter(:))) - (r_ave - mean(r_ave(:)))).^2;
                        if pFcircularX
                            r_iter = (Tuning(icell).corrModelXpos(i,:)/pnumBinsX*2*pi);
                            r_ave = Tuning(icell).meancorrModelXpos/pnumBinsX*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(Tuning(icell).respModelXpos(i,:)' - Tuning(icell).meanrespModelXpos);
                            stdcorrXpos = stdcorrXpos + (pkfold - 1)/pkfold*iterdist.^2;
                            r_iter = (Tuning(icell).corrModelXmax(i,:)/pnumBinsX*2*pi);
                            r_ave = Tuning(icell).meancorrModelXmax/pnumBinsX*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(Tuning(icell).respModelXmax(i,:)' - Tuning(icell).meanrespModelXmax);
                            stdcorrXmax = stdcorrXmax + (pkfold - 1)/pkfold*iterdist.^2;
                            
                            r_iter = (Tuning(icell).corrModelXposNorm(i,:)/pnumBinsX*2*pi);
                            r_ave = Tuning(icell).meancorrModelXposNorm/pnumBinsX*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(Tuning(icell).respModelXpos(i,:)' - Tuning(icell).meanrespModelXpos);
                            stdcorrXposNorm = stdcorrXposNorm + (pkfold - 1)/pkfold*iterdist.^2;
                            r_iter = (Tuning(icell).corrModelXmaxNorm(i,:)/pnumBinsX*2*pi);
                            r_ave = Tuning(icell).meancorrModelXmaxNorm/pnumBinsX*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(Tuning(icell).respModelXmax(i,:)' - Tuning(icell).meanrespModelXmax);
                            stdcorrXmaxNorm = stdcorrXmaxNorm + (pkfold - 1)/pkfold*iterdist.^2;
                        else
                            r_iter = Tuning(icell).corrModelXpos(i,:)';
                            r_ave = Tuning(icell).meancorrModelXpos;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdcorrXpos = stdcorrXpos + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            r_iter = Tuning(icell).corrModelXmax(i,:)';
                            r_ave = Tuning(icell).meancorrModelXmax;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdcorrXmax = stdcorrXmax + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            
                            r_iter = Tuning(icell).corrModelXposNorm(i,:)';
                            r_ave = Tuning(icell).meancorrModelXposNorm;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdcorrXposNorm = stdcorrXposNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            r_iter = Tuning(icell).corrModelXmaxNorm(i,:)';
                            r_ave = Tuning(icell).meancorrModelXmaxNorm;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdcorrXmaxNorm = stdcorrXmaxNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                        end
                        if pFcircularY
                            r_iter = (Tuning(icell).corrModelYpos(i,:)/pnumBinsY*2*pi);
                            r_ave = Tuning(icell).meancorrModelYpos/pnumBinsY*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(Tuning(icell).respModelXpos(i,:)' - Tuning(icell).meanrespModelXpos);
                            stdcorrYpos = stdcorrYpos + (pkfold - 1)/pkfold*iterdist.^2;
                            r_iter = (Tuning(icell).corrModelYmax(i,:)/pnumBinsY*2*pi);
                            r_ave = Tuning(icell).meancorrModelYmax/pnumBinsY*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(Tuning(icell).respModelXmax(i,:)' - Tuning(icell).meanrespModelXmax);
                            stdcorrYmax = stdcorrYmax + (pkfold - 1)/pkfold*iterdist.^2;
                            
                            r_iter = (Tuning(icell).corrModelYposNorm(i,:)/pnumBinsY*2*pi);
                            r_ave = Tuning(icell).meancorrModelYposNorm/pnumBinsY*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(Tuning(icell).respModelXpos(i,:)' - Tuning(icell).meanrespModelXpos);
                            stdcorrYposNorm = stdcorrYposNorm + (pkfold - 1)/pkfold*iterdist.^2;
                            r_iter = (Tuning(icell).corrModelYmaxNorm(i,:)/pnumBinsY*2*pi);
                            r_ave = Tuning(icell).meancorrModelYmaxNorm/pnumBinsY*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(Tuning(icell).respModelXmax(i,:)' - Tuning(icell).meanrespModelXmax);
                            stdcorrYmaxNorm = stdcorrYmaxNorm + (pkfold - 1)/pkfold*iterdist.^2;
                        else
                            r_iter = Tuning(icell).corrModelYpos(i,:)';
                            r_ave = Tuning(icell).meancorrModelYpos;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdcorrYpos = stdcorrYpos + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            r_iter = Tuning(icell).corrModelYmax(i,:)';
                            r_ave = Tuning(icell).meancorrModelYmax;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdcorrYmax = stdcorrYmax + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            
                            r_iter = Tuning(icell).corrModelYposNorm(i,:)';
                            r_ave = Tuning(icell).meancorrModelYposNorm;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdcorrYposNorm = stdcorrYposNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            r_iter = Tuning(icell).corrModelYmaxNorm(i,:)';
                            r_ave = Tuning(icell).meancorrModelYmaxNorm;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdcorrYmaxNorm = stdcorrYmaxNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                        end
                    end
                    if pFcomputePos
                        if pFcircularX
                            r_iter = (Tuning(icell).respModelXpos(i,:)/pnumBinsX*2*pi);
                            r_ave = Tuning(icell).meanrespModelXpos/pnumBinsX*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(Tuning(icell).respModelXpos(i,:)' - Tuning(icell).meanrespModelXpos);
                            stdrespXpos = stdrespXpos + (pkfold - 1)/pkfold*iterdist.^2;
                            r_iter = (Tuning(icell).respModelXmax(i,:)/pnumBinsX*2*pi);
                            r_ave = Tuning(icell).meanrespModelXmax/pnumBinsX*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(Tuning(icell).respModelXmax(i,:)' - Tuning(icell).meanrespModelXmax);
                            stdrespXmax = stdrespXmax + (pkfold - 1)/pkfold*iterdist.^2;
                            
                            r_iter = (Tuning(icell).respModelXposNorm(i,:)/pnumBinsX*2*pi);
                            r_ave = Tuning(icell).meanrespModelXposNorm/pnumBinsX*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(Tuning(icell).respModelXpos(i,:)' - Tuning(icell).meanrespModelXpos);
                            stdrespXposNorm = stdrespXposNorm + (pkfold - 1)/pkfold*iterdist.^2;
                            r_iter = (Tuning(icell).respModelXmaxNorm(i,:)/pnumBinsX*2*pi);
                            r_ave = Tuning(icell).meanrespModelXmaxNorm/pnumBinsX*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(Tuning(icell).respModelXmax(i,:)' - Tuning(icell).meanrespModelXmax);
                            stdrespXmaxNorm = stdrespXmaxNorm + (pkfold - 1)/pkfold*iterdist.^2;
                        else
                            r_iter = Tuning(icell).respModelXpos(i,:)';
                            r_ave = Tuning(icell).meanrespModelXpos;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespXpos = stdrespXpos + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            r_iter = Tuning(icell).respModelXmax(i,:)';
                            r_ave = Tuning(icell).meanrespModelXmax;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespXmax = stdrespXmax + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            
                            r_iter = Tuning(icell).respModelXposNorm(i,:)';
                            r_ave = Tuning(icell).meanrespModelXposNorm;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespXposNorm = stdrespXposNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            r_iter = Tuning(icell).respModelXmaxNorm(i,:)';
                            r_ave = Tuning(icell).meanrespModelXmaxNorm;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespXmaxNorm = stdrespXmaxNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                        end
                        if pFcircularY
                            r_iter = (Tuning(icell).respModelYpos(i,:)/pnumBinsY*2*pi);
                            r_ave = Tuning(icell).meanrespModelYpos/pnumBinsY*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(Tuning(icell).respModelXpos(i,:)' - Tuning(icell).meanrespModelXpos);
                            stdrespYpos = stdrespYpos + (pkfold - 1)/pkfold*iterdist.^2;
                            r_iter = (Tuning(icell).respModelYmax(i,:)/pnumBinsY*2*pi);
                            r_ave = Tuning(icell).meanrespModelYmax/pnumBinsY*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(Tuning(icell).respModelXmax(i,:)' - Tuning(icell).meanrespModelXmax);
                            stdrespYmax = stdrespYmax + (pkfold - 1)/pkfold*iterdist.^2;
                            
                            r_iter = (Tuning(icell).respModelYposNorm(i,:)/pnumBinsY*2*pi);
                            r_ave = Tuning(icell).meanrespModelYposNorm/pnumBinsY*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(Tuning(icell).respModelXpos(i,:)' - Tuning(icell).meanrespModelXpos);
                            stdrespYposNorm = stdrespYposNorm + (pkfold - 1)/pkfold*iterdist.^2;
                            r_iter = (Tuning(icell).respModelYmaxNorm(i,:)/pnumBinsY*2*pi);
                            r_ave = Tuning(icell).meanrespModelYmaxNorm/pnumBinsY*2*pi;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(Tuning(icell).respModelXmax(i,:)' - Tuning(icell).meanrespModelXmax);
                            stdrespYmaxNorm = stdrespYmaxNorm + (pkfold - 1)/pkfold*iterdist.^2;
                            
                            stdslopeXY = stdslopeXY + (pkfold - 1)/pkfold*(Tuning(icell).respModelslopeXY(i) - Tuning(icell).meanrespModelslopeXY).^2;
                            stdphi0XY = stdphi0XY + (pkfold - 1)/pkfold*(Tuning(icell).respModelphi0XY(i) - Tuning(icell).meanrespModelphi0XY).^2;
                            stdrhoXY = stdrhoXY + (pkfold - 1)/pkfold*(Tuning(icell).respModelrhoXY(i) - Tuning(icell).meanrespModelrhoXY).^2;
                        else
                            r_iter = Tuning(icell).respModelYpos(i,:)';
                            r_ave = Tuning(icell).meanrespModelYpos;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespYpos = stdrespYpos + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            r_iter = Tuning(icell).respModelYmax(i,:)';
                            r_ave = Tuning(icell).meanrespModelYmax;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespYmax = stdrespYmax + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            
                            r_iter = Tuning(icell).respModelYposNorm(i,:)';
                            r_ave = Tuning(icell).meanrespModelYposNorm;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespYposNorm = stdrespYposNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            r_iter = Tuning(icell).respModelYmaxNorm(i,:)';
                            r_ave = Tuning(icell).meanrespModelYmaxNorm;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespYmaxNorm = stdrespYmaxNorm + (pkfold - 1)/pkfold*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                        end
                    end
                end
                Tuning(icell).SErespModel = stdresp.^0.5;
                if pFcomputeMarg
                    Tuning(icell).SErespModelX = stdrespX.^0.5;
                    Tuning(icell).SErespModelY = stdrespY.^0.5;
                    
                    Tuning(icell).SErespModelXNorm = stdrespXNorm.^0.5;
                    Tuning(icell).SErespModelYNorm = stdrespYNorm.^0.5;
                end
                if pFcomputeCorr
                    Tuning(icell).SEcorrModelX = stdcorrX.^0.5;
                    Tuning(icell).SEcorrModelY = stdcorrY.^0.5;
                    Tuning(icell).SEcorrModelXpos = stdcorrXpos.^0.5;
                    Tuning(icell).SEcorrModelYpos = stdcorrYpos.^0.5;
                    Tuning(icell).SEcorrModelXmax = stdcorrXmax.^0.5;
                    Tuning(icell).SEcorrModelYmax = stdcorrYmax.^0.5;
                    
                    Tuning(icell).SEcorrModelXposNorm = stdcorrXposNorm.^0.5;
                    Tuning(icell).SEcorrModelYposNorm = stdcorrYposNorm.^0.5;
                    Tuning(icell).SEcorrModelXmaxNorm = stdcorrXmaxNorm.^0.5;
                    Tuning(icell).SEcorrModelYmaxNorm = stdcorrYmaxNorm.^0.5;
                end
                if pFcomputePos
                    Tuning(icell).SErespModelXpos = stdrespXpos.^0.5;
                    Tuning(icell).SErespModelYpos = stdrespYpos.^0.5;
                    Tuning(icell).SErespModelXmax = stdrespXmax.^0.5;
                    Tuning(icell).SErespModelYmax = stdrespYmax.^0.5;
                    
                    Tuning(icell).SErespModelXposNorm = stdrespXposNorm.^0.5;
                    Tuning(icell).SErespModelYposNorm = stdrespYposNorm.^0.5;
                    Tuning(icell).SErespModelXmaxNorm = stdrespXmaxNorm.^0.5;
                    Tuning(icell).SErespModelYmaxNorm = stdrespYmaxNorm.^0.5;
                    
                    if pFcircularY
                        Tuning(icell).SErespModelslopeXY = stdslopeXY.^0.5;
                        Tuning(icell).SErespModelphi0XY = stdphi0XY.^0.5;
                        Tuning(icell).SErespModelrhoXY = stdrhoXY.^0.5;
                    end
                end
                
                %we empty those to save memory. if necessary, SE of these 
                %components can be recovered from the estimates of the 
                %full map
                if pkfold > 2 && pFdiscarditer
                    Tuning(icell).respModel = [];
                    if pFcomputeMarg
                        Tuning(icell).respModelX = [];
                        Tuning(icell).respModelY = [];
                        Tuning(icell).respModelXNorm = [];
                        Tuning(icell).respModelYNorm = [];
                    end
                    if pFcomputeCorr
                        Tuning(icell).corrModelX = [];
                        Tuning(icell).corrModelY = [];
                        Tuning(icell).corrModelXpos = [];
                        Tuning(icell).corrModelXmax = [];
                        Tuning(icell).corrModelYpos = [];
                        Tuning(icell).corrModelYmax = [];
                        Tuning(icell).corrModelXposNorm = [];
                        Tuning(icell).corrModelXmaxNorm = [];
                        Tuning(icell).corrModelYposNorm = [];
                        Tuning(icell).corrModelYmaxNorm = [];
                    end
                    if pFcomputePos
                        Tuning(icell).respModelXpos = [];
                        Tuning(icell).respModelXmax = [];
                        Tuning(icell).respModelYpos = [];
                        Tuning(icell).respModelYmax = [];
                        Tuning(icell).respModelXposNorm = [];
                        Tuning(icell).respModelXmaxNorm = [];
                        Tuning(icell).respModelYposNorm = [];
                        Tuning(icell).respModelYmaxNorm = [];
                        
                        if pFcircularY
                            Tuning(icell).respModelslopeXY = [];
                            Tuning(icell).respModelphi0XY = [];
                            Tuning(icell).respModelrhoXY = [];
                        end
                    end
                end
            end
            
            
            obj.model.tuning = Tuning;
            obj.model.EV = pEV;
            obj.model.L = pL;
            obj.model.Q = pQ;
            obj.model.skaggs = pskaggs;
            obj.model.train_mean = ptrain_mean;
            obj.model.swinX = pswinX;
            obj.model.swinY = pswinY;
            
        end
        
        function [Y, EV] = testMap(obj, X, Yorig)
            minX = min(obj.bins);
            maxX = max(obj.bins);
            X(X<minX) = NaN;
            X(X>maxX) = NaN;
            
            input = normalise1var([minX X' maxX]', obj.numBins);
            input = input(2:end-1);
            
            Y = zeros(length(input),length(obj.model.tuning));
            t = ~isnan(input);
            Y(~t) = NaN;
            for icell = 1:length(obj.model.tuning)
               [~,bestModel] = max(obj.model.EV(:,icell));
               
               tuning = obj.model.tuning(icell).respModel(bestModel,:);
               train  = obj.model.train_mean(bestModel,icell);
               
               Y(t,icell) = tuning(input(t));
               
               EV(icell)  = calCrossValExpVar(train, Yorig(t,icell), Y(t,icell), Yorig(t,icell), Y(t,icell));
            end
        end
    end
end

