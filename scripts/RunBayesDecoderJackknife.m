function obj = RunBayesDecoderJackknife(obj, varname, predictorname, varargin)
if nargin < 4
    %dialog to ask for training and test sets
    if numel(obj.SubsetVal.contrast) > 1
        [contsel,ok] = listdlg('ListString',strsplit(num2str(obj.SubsetVal.contrast)), 'Name', 'Contrast for training', 'InitialValue', 2);
        obj.Bayes.TrainContrast = contsel;
    else
        obj.Bayes.TrainContrast = 1;
    end
    if numel(obj.SubsetVal.gain) > 1
        [gainsel,ok] = listdlg('ListString',strsplit(num2str(obj.SubsetVal.gain)), 'Name', 'Gain for training', 'InitialValue', 2);
        obj.Bayes.TrainGain = gainsel;
    else
        obj.Bayes.TrainGain = 1;
    end
    if numel(obj.SubsetVal.roomlength) > 1
        [roomsel,ok] = listdlg('ListString',strsplit(num2str(obj.SubsetVal.roomlength)), 'Name', 'Room length for training', 'InitialValue', 2);
        obj.Bayes.TrainRoomlength = roomsel;
    else
        obj.Bayes.TrainRoomlength = 1;
    end
    if numel(obj.SubsetVal.outcome) > 1
        [outcomesel,ok] = listdlg('ListString',strsplit(num2str(obj.SubsetVal.outcome)), 'Name', 'Outcome for training', 'InitialValue', 2);
        obj.Bayes.TrainOutcome = outcomesel;%3;%
    else
        obj.Bayes.TrainOutcome = 1;%3;%
    end
    prompt = {'Window size (ms)';'How many theta phase bins';'theta channel#';'How many runs';'Error threshold';'# folds'};
    dlg_title = 'Decoder Parameters';
    num_lines = 1;
    def = {'50';'6';'37';'1';'30';'1'};
    nruns = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(nruns)
        obj.Bayes.smth_win = str2num(nruns{1});
        obj.Bayes.nthetaphsbins = str2num(nruns{2});
        obj.Bayes.thetaChannel = str2num(nruns{3});
        obj.Bayes.nRuns = str2num(nruns{4});
        obj.Bayes.goodidx_th = str2num(nruns{5});
        obj.Bayes.Nfold = str2num(nruns{6});
    else
        obj.Bayes.smth_win = 50;
        obj.Bayes.nthetaphsbins = 6;
        obj.Bayes.thetaChannel = 37;
        obj.Bayes.nRuns = 1;
        obj.Bayes.goodidx_th = inf;
        obj.Bayes.Nfold = 0;
    end
    
    obj.Bayes.type = 'mean';   
else
    pnames = {'train_contrast' 'train_gain' 'train_roomlength' 'train_outcome' 'type' 'smth_win' 'nthetaphsbins' 'thetaChannel' 'nruns' 'error_th' 'Nfold'};
    dflts  = {obj.SubsetVal.contrast obj.SubsetVal.gain obj.SubsetVal.roomlength obj.SubsetVal.outcome 'mean' 50 6 1 inf 0};
    [obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome , obj.Bayes.type, obj.Bayes.smth_win, obj.Bayes.nthetaphsbins, obj.Bayes.thetaChannel, obj.Bayes.nRuns, obj.Bayes.goodidx_th, obj.Bayes.Nfold] = internal.stats.parseArgs(pnames,dflts,varargin{:});
end

X = obj.data.es.(varname);
Y = obj.data.es.(predictorname);
nbins = (floor(max(obj.data.es.traj-1)/50) + 1)*50;

obj.Bayes.Posterior = cell(obj.Bayes.nthetaphsbins,1);
obj.Bayes.nPosterior = cell(obj.Bayes.nthetaphsbins,1);
obj.Bayes.PosError = cell(obj.Bayes.nthetaphsbins,1);
obj.Bayes.DistError = cell(obj.Bayes.nthetaphsbins,1);

obj.Bayes.PosteriorSE = cell(obj.Bayes.nthetaphsbins,1);
obj.Bayes.nPosteriorSE = cell(obj.Bayes.nthetaphsbins,1);

obj.Bayes.Posterior0 = zeros(numel(X),nbins);
obj.Bayes.nPosterior0 = zeros(numel(X),nbins);
obj.Bayes.PosError0 = zeros(numel(X),nbins);
obj.Bayes.DistError0 = zeros(numel(X),nbins);

obj.Bayes.Posterior0SE = zeros(numel(X),nbins);
obj.Bayes.nPosterior0SE = zeros(numel(X),nbins);

for phs = 1:obj.Bayes.nthetaphsbins
    for irun = 1:obj.Bayes.nRuns
        idx = obj.getSubsets(obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome);
        if irun > 1
            obj.Bayes.goodidx = getcleandecidx(obj.Bayes.Posterior0, obj.Bayes.X,obj.Bayes.goodidx_th);
        else
            obj.Bayes.goodidx = true(size(idx));
        end
        
        phs0 = (phs-1)*360/obj.Bayes.nthetaphsbins;
        phsval = mod(obj.data.es.LFPphase(:,obj.Bayes.thetaChannel)-180,360);
        phsidx = ismember(phsval,mod(phs0:(phs0 + 360/obj.Bayes.nthetaphsbins),360));
        idx0 = idx & obj.Bayes.goodidx & phsidx;
        
        numBins = (floor(max(obj.data.es.traj-1)/50) + 1)*50;
        
        obj.Bayes.prediction{phs} = zeros(numel(X),1);
        obj.Bayes.X = zeros(numel(X),1);
        obj.Bayes.Posterior{phs} = zeros(numel(X),numBins);
        obj.Bayes.nPosterior{phs} = zeros(numel(X),numBins);
        obj.Bayes.PosError{phs} = zeros(numel(X),numBins);
        obj.Bayes.DistError{phs} = zeros(numel(X),numBins);
        
        obj.Bayes.predictionSE{phs} = zeros(numel(X),1);
        obj.Bayes.PosteriorSE{phs} = zeros(numel(X),numBins);
        obj.Bayes.nPosteriorSE{phs} = zeros(numel(X),numBins);
        
        predictionSqrSum = zeros(numel(X),1);
        PosteriorSqrSum = zeros(numel(X),numBins);
        nPosteriorSqrSum = zeros(numel(X),numBins);
        
        for kite = 1:obj.Bayes.Nfold
            disp(['Bayesian decoder: run #' num2str(irun) ' out of ' num2str(obj.Bayes.nRuns) ' ; phs#' num2str(phs) ' out of ' num2str(obj.Bayes.nthetaphsbins) ' ; fold #' num2str(kite) ' out of ' num2str(obj.Bayes.Nfold)  ]);
            idx = idx0;
            if obj.Bayes.Nfold > 1
                idx0true = find(idx0);
                if kite < obj.Bayes.Nfold
                    idxfold = ((kite-1)*floor(numel(idx0true)/obj.Bayes.Nfold)+1):(kite*floor(numel(idx0true)/obj.Bayes.Nfold));
                else
                    idxfold = ((kite-1)*floor(numel(idx0true)/obj.Bayes.Nfold)+1):numel(idx0true);
                end
                idx(idx0true(idxfold)) = false;
            end
            
            obj.Bayes.decOder = TbayesDecoder;
            obj.Bayes.decOder.numBins = numBins;
            obj.Bayes.decOder.kfold = obj.Bayes.Nfold;
            obj.Bayes.decOder.FoptiSmooth = false;
            
            [obj.Bayes.decOder, pred, Xdec, Posterior, nPosterior] = obj.Bayes.decOder.trainDecoder(X(idx), Y(idx,:), obj.Bayes.smth_win);
            
            IDX = find(idx);
            IDX = IDX(1:min(end,numel(pred)));
            
%             obj.Bayes.prediction{phs}(IDX,:) = pred;
%             obj.Bayes.X(IDX,:) = Xdec(1:min(end,numel(pred)));
%             obj.Bayes.Posterior{phs}(IDX,:) = Posterior;
%             obj.Bayes.nPosterior{phs}(IDX,:) = nPosterior;
            
%             obj.Bayes.randPosteriorMean{phs}(IDX,:) = RandPosteriormean;
%             obj.Bayes.randPosteriorStd{phs}(IDX,:) = RandPosteriorstd;
            
            [pred, Xdec, Posterior, nPosterior] = obj.Bayes.decOder.predictBayesDecoder(X(~idx & idx0,:), Y(~idx & idx0,:), obj.Bayes.smth_win, 'All');
            
            obj.Bayes.prediction{phs}(~idx & idx0,:) = obj.Bayes.prediction{phs}(~idx & idx0,:) + sum(pred,2)/obj.Bayes.Nfold;
            obj.Bayes.X(~idx & idx0,:) = Xdec(1:min(end,size(pred,1)));
            obj.Bayes.Posterior{phs}(~idx & idx0,:) = obj.Bayes.Posterior{phs}(~idx & idx0,:) + sum(Posterior,3)/obj.Bayes.Nfold;
            obj.Bayes.nPosterior{phs}(~idx & idx0,:) = obj.Bayes.nPosterior{phs}(~idx & idx0,:) + size(nPosterior,3)/obj.Bayes.Nfold;
            
            predictionSqrSum(~idx & idx0,:) = predictionSqrSum(~idx & idx0,:) + sum(pred.^2,2)/obj.Bayes.Nfold;
            PosteriorSqrSum(~idx & idx0,:) = PosteriorSqrSum(~idx & idx0,:) + sum(Posterior.^2,3)/obj.Bayes.Nfold;
            nPosteriorSqrSum(~idx & idx0,:) = nPosteriorSqrSum(~idx & idx0,:) + sum(nPosterior.^2,3)/obj.Bayes.Nfold;
            
            [pred, Xdec, Posterior, nPosterior] = obj.Bayes.decOder.predictBayesDecoder(X(~idx0,:), Y(~idx0,:), obj.Bayes.smth_win, obj.Bayes.type);
            
            obj.Bayes.prediction{phs}(~idx0,:) = obj.Bayes.prediction{phs}(~idx0,:) + pred/obj.Bayes.Nfold;
            obj.Bayes.X(~idx0,:) = Xdec(1:min(end,size(pred,1)));
            obj.Bayes.Posterior{phs}(~idx0,:) = obj.Bayes.Posterior{phs}(~idx0,:) + Posterior/obj.Bayes.Nfold;
            obj.Bayes.nPosterior{phs}(~idx0,:) = obj.Bayes.nPosterior{phs}(~idx0,:) + nPosterior/obj.Bayes.Nfold;
            
            predictionSqrSum(~idx0,:) = predictionSqrSum(~idx0,:) + (pred.^2)/obj.Bayes.Nfold;
            PosteriorSqrSum(~idx0,:) = PosteriorSqrSum(~idx0,:) + (Posterior.^2)/obj.Bayes.Nfold;
            nPosteriorSqrSum(~idx0,:) = nPosteriorSqrSum(~idx0,:) + (nPosterior.^2)/obj.Bayes.Nfold;
            
%             obj.Bayes.randPosteriorMean{phs}(~idx,:) = RandPosteriormean;
%             obj.Bayes.randPosteriorStd{phs}(~idx,:) = RandPosteriorstd;
            
            % baseline = 1./obj.Bayes.decOder.numBins;
            % obj.Bayes.Posterior = medfilt1(exp(obj.Bayes.Posterior*log(2))*baseline,round(50*60/1000));
            % obj.Bayes.Posterior = obj.Bayes.Posterior./(sum(obj.Bayes.Posterior,2)*ones(1,size(obj.Bayes.Posterior,2)));
            % obj.Bayes.Posterior = log2(obj.Bayes.Posterior/baseline);
            
            Prange = size(obj.Bayes.Posterior{phs},2);
            Xrange = max(obj.Bayes.X);
            for i = 1:Xrange
                obj.Bayes.PosError{phs}(obj.Bayes.X == i,:) = circshift(obj.Bayes.Posterior{phs}(obj.Bayes.X == i,:),floor(Prange/2)-i,2);
            end
            
            Prange = size(obj.Bayes.Posterior{phs},2);
            obj.Bayes.D = floor(obj.data.es.traj ./ (obj.data.es.gain/mode(obj.data.es.gain)));
            Xrange = max(obj.Bayes.D);
            for i = 1:Xrange
                obj.Bayes.DistError{phs}(obj.Bayes.D == i,:) = circshift(obj.Bayes.Posterior{phs}(obj.Bayes.D == i,:),floor(Prange/2)-i,2);
            end
        end
        obj.Bayes.predictionSE{phs} = sqrt((predictionSqrSum - obj.Bayes.prediction{phs}.^2)/(obj.Bayes.Nfold/(obj.Bayes.Nfold-1)));
        obj.Bayes.PosteriorSE{phs} = sqrt((PosteriorSqrSum - obj.Bayes.Posterior{phs}.^2)/(obj.Bayes.Nfold/(obj.Bayes.Nfold-1)));
        obj.Bayes.nPosteriorSE{phs} = sqrt((nPosteriorSqrSum - obj.Bayes.nPosterior{phs}.^2)/(obj.Bayes.Nfold/(obj.Bayes.Nfold-1)));

        obj.Bayes.Posterior0(phsidx,:) = obj.Bayes.Posterior{phs}(phsidx,:);
        obj.Bayes.nPosterior0(phsidx,:) = obj.Bayes.nPosterior{phs}(phsidx,:);
        obj.Bayes.PosError0(phsidx,:) = obj.Bayes.PosError{phs}(phsidx,:);
        obj.Bayes.DistError0(phsidx,:) = obj.Bayes.DistError{phs}(phsidx,:);
        
        obj.Bayes.Posterior0SE(phsidx,:) = obj.Bayes.PosteriorSE{phs}(phsidx,:);
        obj.Bayes.nPosterior0SE(phsidx,:) = obj.Bayes.nPosteriorSE{phs}(phsidx,:);
        
%         obj.Bayes.randPosterior0Mean(phsidx,:) = obj.Bayes.randPosteriorMean{phs}(phsidx,:);
%         obj.Bayes.randPosterior0Std(phsidx,:) = obj.Bayes.randPosteriorStd{phs}(phsidx,:);
    end
end
end