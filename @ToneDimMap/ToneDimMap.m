classdef ToneDimMap < TspikeMap
    %     Create a spikeMap class object of type '1D'. The decoder can be trained using the
    %     function 'trainDecoder' and used to decode with the function
    %     'predictVariable'
    %     [obj, prediction, X] = trainSpikeMap(obj, X, Y, Xsmth_win);
    %     
    % Aman Saleem
    % Jan 2014
    properties
        kfold;       % number of partitions made to train the decoder
        train_mean;  % mean response used for training
        CVO;         % cross-validation structure
        bins;        % the values of the bin limits
        numBins;     % number of bins by which the variable is discretized
        Fcircular;   % if true, circular smoothing
        qthreshold;
        FcomputePos;
        Fdiscarditer;
    end
    
    methods
        function obj = ToneDimMap(varargin)
            
            pnames = {'dimensionality' 'variable' 'variable_range'...
                'Xsmth_win' 'model'...
                'kfold' 'performance'...
                'train_mean' 'CVO' 'sampleRate' 'numBins' 'Fcircular' 'qthreshold' 'FcomputePos' 'Fdiscarditer'};
            dflts  = {'1D' 'P' []...
                1 []...
                20 []...
                [] [] 60 100 true 1 true true};
            [obj.dimensionality, obj.variable, obj.variable_range, ...
                obj.Xsmth_win, obj.model,...
                obj.kfold, obj.performance,...
                obj.train_mean, obj.CVO, obj.sampleRate, obj.numBins, obj.Fcircular, obj.qthreshold, obj.FcomputePos, obj.Fdiscarditer] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        
        function [obj, Prediction, X] = trainSpikeMap(obj, X, Y, T, Xsmth_win, FoptiSmooth)
            if nargin<5
                Xsmth_win = obj.Xsmth_win;
            else
                obj.Xsmth_win = Xsmth_win;
            end
            if nargin < 6
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
            
            if size(X,1) ~= size(Y,1)
                error('wrong dimension for X in ToneDimMap');
            end
            
            if sum(isnan(X(:)))>0 || sum(isnan(Y(:)))>0
                display('WARNING!!! Nans in the data: making a temp fix');
                t = ones(size(X));
                t(isnan(sum(X,2))) = 0;
                t(isnan(sum(Y,2))) = 0;
                X = X(t>0,:);
                Y = Y(t>0,:);
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
%             [X, obj.bins] = normalise1var(X, obj.numBins);
            
            Ncells = size(Y,2);
            pEV = zeros(obj.CVO.kfold,Ncells);
            pL = zeros(obj.CVO.kfold,Ncells);
            pQ = zeros(obj.CVO.kfold,Ncells);
            pskaggs = zeros(obj.CVO.kfold,Ncells);
            ptrain_mean = zeros(obj.CVO.kfold,Ncells);
            pswin = zeros(obj.CVO.kfold,Ncells);
            for icell = 1:Ncells
                Tuning(icell).respModel = zeros(obj.CVO.kfold,obj.numBins);
                
                Tuning(icell).meanrespModel = [];
                Tuning(icell).SErespModel = [];
                
                if obj.FcomputePos
                    Tuning(icell).respModelXpos = zeros(obj.CVO.kfold,1);
                    Tuning(icell).respModelXmax = zeros(obj.CVO.kfold,1);
                    
                    Tuning(icell).meanrespModelXpos = [];
                    Tuning(icell).meanrespModelXmax = [];
                    Tuning(icell).SErespModelXpos = [];
                    Tuning(icell).SErespModelXmax = [];
                    
                    Tuning(icell).fieldidx = [];
                end
            end
            
            pnumBins = obj.numBins;
            pbins = obj.bins;
            pCVO = obj.CVO;
            pkfold = obj.CVO.kfold;
            psampleRate = obj.sampleRate;
            pXsmth_win = obj.Xsmth_win;
            pFcircular = obj.Fcircular;
            pFcomputePos = obj.FcomputePos;
            pqthreshold = obj.qthreshold;
            pFdiscarditer = obj.Fdiscarditer;
            
            for iter = 1:obj.CVO.kfold
                parfor icell = 1:Ncells
                    [model, ~] = get1Dmap(Y(:,icell), X(:,min(icell,end))', 1./T, pnumBins, pbins, pCVO, iter, psampleRate, pXsmth_win, FoptiSmooth, pFcircular);
                    Tuning(icell).respModel(iter,:) = model.tuning;
                    pEV(iter,icell) = model.EV;
                    pL(iter,icell) = model.L;
                    pQ(iter,icell) = model.Q;
                    pskaggs(iter,icell) = model.skaggs;
                    ptrain_mean(iter,icell) = model.train_mean;
                    pswin(iter,icell) = model.swin;
                end
            end
            
            Xmesh = 1:obj.numBins;
            parfor icell = 1:Ncells
                Tuning(icell).meanrespModel = squeeze(mean(Tuning(icell).respModel,1));
                
                if pFcomputePos
                    [~,fieldX] = findfield(Tuning(icell).meanrespModel,pqthreshold);
                    Tuning(icell).fieldidx = fieldX;
                    outfieldX = find(~ismember(1:pnumBins,fieldX));
                    for iter = 1:pkfold
                        tuningX0 = squeeze(Tuning(icell).respModel(iter,:));
                        tuningX = tuningX0;
                        if ~isempty(fieldX)
                            tuningX(outfieldX) = 0;%min(tuningX(fieldX));
                        end
                        if pFcircular
                            Tuning(icell).respModelXpos(iter) = getCircularAverage(tuningX',0,1);
                            Tuning(icell).respModelXmax(iter) = getCircularAverage(tuningX0',0,0.01,0.05);
                        else
                            Tuning(icell).respModelXpos(iter) = sum(tuningX.*Xmesh)./sum(tuningX);
                            [~, Tuning(icell).respModelXmax(iter)] = max(tuningX0);
                        end
                    end
                    
                    tuningX0 = squeeze(Tuning(icell).meanrespModel);
                    tuningX = tuningX0;
                    if ~isempty(fieldX)
                        tuningX(outfieldX) = 0;%min(tuningX(fieldX));
                    end
                    if pFcircular
                        Tuning(icell).meanrespModelXpos = getCircularAverage(tuningX',0,1);%mod(circ_mean(Tuning(icell).respModelXpos/pnumBins*2*pi,[],1)+2*pi,2*pi)*pnumBins/(2*pi);
                        Tuning(icell).meanrespModelXmax = getCircularAverage(tuningX0',0,0.01,0.05);%mod(circ_mean(Tuning(icell).respModelXmax/pnumBins*2*pi,[],1)+2*pi,2*pi)*pnumBins/(2*pi);
                    else
                        Tuning(icell).meanrespModelXpos = sum(tuningX.*Xmesh)./sum(tuningX);%squeeze(nanmean(Tuning(icell).respModelXpos,1));
                        [~, Tuning(icell).meanrespModelXmax] = max(tuningX0);%squeeze(nanmean(Tuning(icell).respModelXmax,1));
                    end
                end
                
                stdresp = 0;
                stdrespXpos = 0;
                stdrespXmax = 0;
                for i = 1:pkfold
                    r_iter = Tuning(icell).respModel(i,:);
                    r_ave = Tuning(icell).meanrespModel;
                    stdresp = stdresp + (pkfold - 1)/pkfold*((r_iter - mean(r_iter(:))) - (r_ave - mean(r_ave(:)))).^2;
                    if pFcomputePos
                        if pFcircular
                            iterdist = circ_dist((Tuning(icell).respModelXpos(i)/pnumBins*2*pi),Tuning(icell).meanrespModelXpos/pnumBins*2*pi)*pnumBins/(2*pi);
                            stdrespXpos = stdrespXpos + (pkfold - 1)/pkfold*iterdist.^2;
                            iterdist = circ_dist((Tuning(icell).respModelXmax(i)/pnumBins*2*pi),Tuning(icell).meanrespModelXmax/pnumBins*2*pi)*pnumBins/(2*pi);
                            stdrespXmax = stdrespXmax + (pkfold - 1)/pkfold*iterdist.^2;
                        else
                            stdrespXpos = stdrespXpos + (pkfold - 1)/pkfold*(Tuning(icell).respModelXpos(i) - Tuning(icell).meanrespModelXpos).^2;
                            stdrespXmax = stdrespXmax + (pkfold - 1)/pkfold*(Tuning(icell).respModelXmax(i) - Tuning(icell).meanrespModelXmax).^2;
                        end
                    end
                end
                Tuning(icell).SErespModel = stdresp.^0.5;
                if pFcomputePos
                    Tuning(icell).SErespModelXpos = stdrespXpos.^0.5;
                    Tuning(icell).SErespModelXmax = stdrespXmax.^0.5;
                end
                if pkfold > 2 && pFdiscarditer
                    Tuning(icell).respModel = [];
                    if pFcomputePos
                        Tuning(icell).respModelXpos = [];
                        Tuning(icell).respModelXmax = [];
                    end
                end
            end
            obj.model.tuning = Tuning;
            obj.model.EV = pEV;
            obj.model.L = pL;
            obj.model.Q = pQ;
            obj.model.skaggs = pskaggs;
            obj.model.train_mean = ptrain_mean;
            obj.model.swin = pswin;
        end
        
        function [obj, Prediction, X] = trainSpikeMap2(obj, X, Y, T, Xsmth_win, FoptiSmooth)
            if nargin<5
                Xsmth_win = obj.Xsmth_win;
            else
                obj.Xsmth_win = Xsmth_win;
            end
            if nargin < 6
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
            
            if size(X,1) ~= size(Y,1)
                error('wrong dimension for X in ToneDimMap');
            end
            
            if sum(isnan(X))>0 || sum(isnan(Y(:)))>0
                display('WARNING!!! Nans in the data: making a temp fix');
                t = ones(size(X));
                t(isnan(X)) = 0;
                t(isnan(sum(Y,2))) = 0;
                X = X(t>0);
                Y = Y(t>0,:);
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
%             [X, obj.bins] = normalise1var(X, obj.numBins);
            
            pred_temp = cell(obj.CVO.kfold,1);
            
            for iter = 1:obj.CVO.kfold
                
                %     Ytrain  = Y(obj.CVO.train{iter},:);
                %     Ycv     = Y(obj.CVO.cv{iter},:);
                %     Ytest   = Y(obj.CVO.test{iter},:);
                %
                %     Vtrain  = X(obj.CVO.train{iter},:);
                %     Vcv     = X(obj.CVO.cv{iter},:);
                %     Vtest   = X(obj.CVO.test{iter},:);
                
                %     if calTrainMean
                %         train_mean = mean(Vtrain,1);
                %     else
                %         train_mean = obj.train_mean;
                %     end
                
                % Getting the 1D map for each neuron (calculating the place fields)
                % ...the main training component
                %     respModel = zeros(size(Y,2),obj.numBins);
                for icell = 1:size(Y,2)
                    [model, pred] = get1Dmap2(Y(:,icell), X', 1./T, obj.numBins, obj.bins, obj.CVO, iter, obj.sampleRate, obj.Xsmth_win, FoptiSmooth, obj.Fcircular);
                    obj.model.tuning(icell).respModel(iter,:,:) = model.tuning;
                    obj.model.EV(iter,icell) = model.EV;
                    obj.model.L(iter,icell) = model.L;
                    obj.model.Q(iter,icell) = model.Q;
                    obj.model.skaggs(iter,icell) = model.skaggs;
                    obj.model.train_mean(iter,icell) = model.train_mean;
                    obj.model.swin(iter,icell) = model.swin;
                    pred_temp{iter}(icell,:) = pred;                                                            
                end
                obj.model.EV(iter,obj.model.EV(iter,:)<0) = 0;                
            end
            for icell = 1:size(Y,2)
                obj.model.tuning(icell).meanrespModel = squeeze(mean(obj.model.tuning(icell).respModel,1));
                stdresp = 0;
                for i = 1:obj.CVO.kfold
%                     stdresp = stdresp + (obj.CVO.kfold - 1)/obj.CVO.kfold*(obj.model.tuning(icell).respModel(i,:) - obj.model.tuning(icell).meanrespModel).^2;
                end
                obj.model.tuning(icell).SErespModel = stdresp.^0.5;
            end
%             for icell = 1:size(Y,2)
%                 fitobject = fit(obj.bins',mean(obj.model.tuning(icell).respModel, 1)',@(A,k,p,x)(A*exp(k*cos(x/100*2*pi - p)-1)), 'StartPoint', [1 1 0] );
%                 vonmisesfit = feval(fitobject,obj.bins);
%                 [~ , obj.model.position(icell)] = max(vonmisesfit);
%             end
            
            if size(pred_temp{1},1)>1
                for iter = 1:size(pred_temp,1)
                    Prediction = [Prediction pred_temp{iter}];
                end
            else
                for iter = 1:size(pred_temp,1)
                    Prediction = [Prediction pred_temp{iter}];
                end
                Prediction = Prediction';
            end
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

