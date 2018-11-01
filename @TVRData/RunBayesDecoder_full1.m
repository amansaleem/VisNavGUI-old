function obj = RunBayesDecoder(obj, varname, predictorname, varargin)
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
    prompt = {'Window size (ms)';'How many theta phase bins';'theta channel#';'How many runs';'Error threshold';'# folds';'Optimal smoothing?'};
    dlg_title = 'Decoder Parameters';
    num_lines = 1;
    def = {'50';'6';'37';'1';'30';'5';'1'};
    nruns = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(nruns)
        obj.Bayes.smth_win = str2num(nruns{1});
        obj.Bayes.nthetaphsbins = str2num(nruns{2});
        obj.Bayes.thetaChannel = str2num(nruns{3});
        obj.Bayes.nRuns = str2num(nruns{4});
        obj.Bayes.goodidx_th = str2num(nruns{5});
        obj.Bayes.kfold = str2num(nruns{6});
        obj.Bayes.FoptiSmooth = boolean(str2num(nruns{7}));
    else
        obj.Bayes.smth_win = 50;
        obj.Bayes.nthetaphsbins = 6;
        obj.Bayes.thetaChannel = 37;
        obj.Bayes.nRuns = 1;
        obj.Bayes.goodidx_th = inf;
        obj.Bayes.kfold = 5;
        obj.Bayes.FoptiSmooth = true;
    end
    
    obj.Bayes.type = 'mean';   
else
    pnames = {'train_contrast' 'train_gain' 'train_roomlength' 'train_outcome' 'type' 'smth_win' 'nthetaphsbins' 'thetaChannel' 'nruns' 'error_th' 'kfold' 'FoptiSmooth'};
    dflts  = {obj.SubsetVal.contrast obj.SubsetVal.gain obj.SubsetVal.roomlength obj.SubsetVal.outcome 'mean' 50 6 1 inf 0};
    [obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome , obj.Bayes.type, obj.Bayes.smth_win, obj.Bayes.nthetaphsbins, obj.Bayes.thetaChannel, obj.Bayes.nRuns, obj.Bayes.goodidx_th, obj.Bayes.kfold, obj.Bayes.FoptiSmooth] = internal.stats.parseArgs(pnames,dflts,varargin{:});
end
tic

Nmodel = 1;
obj.Bayes.Posterior = cell(obj.Bayes.nthetaphsbins,1);
obj.Bayes.nPosterior = cell(obj.Bayes.nthetaphsbins,1);
obj.Bayes.PosError = cell(obj.Bayes.nthetaphsbins,1);
obj.Bayes.DistError = cell(obj.Bayes.nthetaphsbins,1);
obj.Bayes.prediction = cell(obj.Bayes.nthetaphsbins,1);

ntimebins = numel(obj.data.es.(varname));
% if Nmodel > 1
    nbins = max(floor(obj.data.es.(varname)./(obj.data.es.gain/mode(obj.data.es.gain))))+1;
% else
%     nbins = max(floor(obj.data.es.(varname)))+1;
% end

obj.Bayes.X = cell(1,Nmodel+1);
obj.Bayes.Posterior0 = zeros(ntimebins,nbins);
obj.Bayes.nPosterior0 = zeros(ntimebins,nbins);
obj.Bayes.PosError0 = zeros(ntimebins,nbins);
obj.Bayes.DistError0 = zeros(ntimebins,nbins);
obj.Bayes.bestModel0 = zeros(ntimebins,nbins);
obj.Bayes.prediction0 = zeros(1,nbins);  
obj.Bayes.modelError = cell(obj.Bayes.nthetaphsbins,1);
obj.Bayes.bestModel = cell(obj.Bayes.nthetaphsbins,1);
    
for phs = 1:obj.Bayes.nthetaphsbins
    for irun = 1:obj.Bayes.nRuns
        obj.Bayes.modelError{phs} = zeros(numel(obj.data.es.(varname)),Nmodel);
        obj.Bayes.bestModel{phs} = zeros(numel(obj.data.es.(varname)),1);
        
        predictionm = cell(Nmodel,1);
        Xm = cell(Nmodel,1);
        PosteriorM = cell(Nmodel,1);
        nPosteriorM = cell(Nmodel,1);
        PosErrorM = cell(Nmodel,1);
        DistErrorM = cell(Nmodel,1);
        
        for xmodel = 1:Nmodel
            X = obj.data.es.(varname);
            Y = obj.data.es.(predictorname);
            xmax = (floor(max(X-1)/50) + 1)*50;
            
%             if Nmodel > 1
                x0 = 50;%(xmodel-1)*floor(xmax/Nmodel);
                X = ((X - x0) ./ (obj.data.es.gain/mode(obj.data.es.gain)));
                X(X < 0) = X(X < 0)+xmax./(obj.data.es.gain(X < 0)/mode(obj.data.es.gain));
%             end
            nbins = max(floor(X))+1;
            
            if xmodel == 1
                obj.Bayes.prediction{phs} = zeros(numel(X),1);
                obj.Bayes.X = zeros(numel(X),1);
                obj.Bayes.Posterior{phs} = zeros(numel(X),nbins);
                obj.Bayes.nPosterior{phs} = zeros(numel(X),nbins);
                obj.Bayes.PosError{phs} = zeros(numel(X),nbins);
                obj.Bayes.DistError{phs} = zeros(numel(X),nbins);
            end
            
            disp(['Bayesian decoder: run #' num2str(irun) ' out of ' num2str(obj.Bayes.nRuns) ' ; phs#' num2str(phs) ' out of ' num2str(obj.Bayes.nthetaphsbins) ' ; model#' num2str(xmodel) ' out of ' num2str(Nmodel)]);
            
            obj.Bayes.decOder = TbayesDecoder;
            obj.Bayes.decOder.numBins =nbins;
            obj.Bayes.decOder.kfold = obj.Bayes.kfold;
            obj.Bayes.decOder.FoptiSmooth = obj.Bayes.FoptiSmooth;
            
            idx = obj.getSubsets(obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome);            
            
            if irun > 1
                obj.Bayes.goodidx = getcleandecidx(obj.Bayes.Posterior0, obj.Bayes.X,obj.Bayes.goodidx_th);
            else
                obj.Bayes.goodidx = true(size(idx));
            end
            
            phs0 = (phs-1)*360/obj.Bayes.nthetaphsbins;
            phsval = mod(obj.data.es.LFPphase(:,obj.Bayes.thetaChannel)-180,360);
            phsidx = ismember(phsval,mod(phs0:(phs0 + 360/obj.Bayes.nthetaphsbins),360));
            idx = idx & obj.Bayes.goodidx & phsidx;
            
            predictionm{xmodel} = zeros(numel(X),1);
            Xm{xmodel} = zeros(numel(X),1);
            PosteriorM{xmodel} = zeros(numel(X),obj.Bayes.decOder.numBins);
            nPosteriorM{xmodel} = zeros(numel(X),obj.Bayes.decOder.numBins);
            PosErrorM{xmodel} = zeros(numel(X),obj.Bayes.decOder.numBins);
            DistErrorM{xmodel} = zeros(numel(X),obj.Bayes.decOder.numBins);
                                    
            [obj.Bayes.decOder, pred, Xdec, Posterior, nPosterior] = obj.Bayes.decOder.trainDecoder(X(idx), Y(idx,:), obj.Bayes.smth_win);
            
            IDX = find(idx);
            IDX = IDX(1:min(end,numel(pred)));
            
            obj.Bayes.modelError{phs}(IDX,xmodel) = mod(abs(pred - Xdec(1:min(end,numel(pred)))),floor(nbins/2));
            
            predictionm{xmodel}(IDX,:) = pred;
            Xm{xmodel}(IDX,:) = Xdec(1:min(end,numel(pred)));
            PosteriorM{xmodel}(IDX,:) = Posterior;
            nPosteriorM{xmodel}(IDX,:) = nPosterior;
            
            [pred, Xdec, Posterior, nPosterior] = obj.Bayes.decOder.predictBayesDecoder(X(~idx,:), Y(~idx,:), obj.Bayes.smth_win, obj.Bayes.type);
            
            obj.Bayes.modelError{phs}(~idx,xmodel) = mod(abs(pred - Xdec(1:min(end,numel(pred)))),floor(nbins/2));
            
            predictionm{xmodel}(~idx,:) = pred;
            Xm{xmodel}(~idx,:) = Xdec(1:min(end,numel(pred)));
            PosteriorM{xmodel}(~idx,:) = Posterior;
            nPosteriorM{xmodel}(~idx,:) = nPosterior;
            
            %         baseline = 1./obj.Bayes.decOder.numBins;
            %         h = histogram(obj.Bayes.X,obj.Bayes.decOder.numBins);
            %         Prior = (h.Values(:)')/sum(h.Values);
            %         Prior = smooth(medfilt1([Prior Prior Prior],5),5);
            %         Prior = Prior(101:200);
            %         Prior = Prior/sum(Prior);
            %         obj.Bayes.Posterior{phs} = obj.Bayes.Posterior{phs}*baseline.*repmat(smooth(medfilt1(Prior,3),3)',size(obj.Bayes.Posterior{phs},1),1);
            %         obj.Bayes.Posterior{phs} = obj.Bayes.Posterior{phs}./repmat(sum(obj.Bayes.Posterior{phs},2),1,100);
            %         obj.Bayes.Posterior{phs} = obj.Bayes.Posterior{phs}./baseline;
            
            % baseline = 1./obj.Bayes.decOder.numBins;
            % obj.Bayes.Posterior = medfilt1(exp(obj.Bayes.Posterior*log(2))*baseline,round(50*60/1000));
            % obj.Bayes.Posterior = obj.Bayes.Posterior./(sum(obj.Bayes.Posterior,2)*ones(1,size(obj.Bayes.Posterior,2)));
            % obj.Bayes.Posterior = log2(obj.Bayes.Posterior/baseline);
            
            Prange = size(PosteriorM{xmodel},2);
            Xrange = max(Xm{xmodel});
            for i = 1:Xrange
                PosErrorM{xmodel}(Xm{xmodel} == i,:) = circshift(PosteriorM{xmodel}(Xm{xmodel} == i,:),floor(Prange/2)-i,2);
            end
            
            Prange = size(PosteriorM{xmodel},2);
            obj.Bayes.D = floor(obj.data.es.traj ./ (obj.data.es.gain/mode(obj.data.es.gain)));
            Xrange = max(obj.Bayes.D);
            for i = 1:Xrange
                DistErrorM{xmodel}(obj.Bayes.D == i,:) = circshift(PosteriorM{xmodel}(obj.Bayes.D == i,:),floor(Prange/2)-i,2);
            end            
        end
        [~, bestmodel] = min(obj.Bayes.modelError{phs}, [],2);
        obj.Bayes.bestModel{phs}(phsidx) = (double(bestmodel(phsidx))-1)*floor(xmax/Nmodel);
        for xmodel = 1:Nmodel
            x0 = 50;%(xmodel-1)*floor(xmax/Nmodel);
            obj.Bayes.prediction{phs}(phsidx & bestmodel == xmodel,:) = mod(x0 - 1 + predictionm{xmodel}(phsidx & bestmodel == xmodel),xmax);
            obj.Bayes.Posterior{phs}(phsidx & bestmodel == xmodel,:) = circshift(PosteriorM{xmodel}(phsidx & bestmodel == xmodel,:),x0,2);
            obj.Bayes.nPosterior{phs}(phsidx & bestmodel == xmodel,:) = circshift(nPosteriorM{xmodel}(phsidx & bestmodel == xmodel,:),x0,2);
            obj.Bayes.PosError{phs}(phsidx & bestmodel == xmodel,:) = PosErrorM{xmodel}(phsidx & bestmodel == xmodel,:);
            obj.Bayes.DistError{phs}(phsidx & bestmodel == xmodel,:) = DistErrorM{xmodel}(phsidx & bestmodel == xmodel,:);
        end
        obj.Bayes.Posterior0(phsidx,:) = obj.Bayes.Posterior{phs}(phsidx,:);
        obj.Bayes.nPosterior0(phsidx,:) = obj.Bayes.nPosterior{phs}(phsidx,:);
        obj.Bayes.PosError0(phsidx,:) = obj.Bayes.PosError{phs}(phsidx,:);
        obj.Bayes.DistError0(phsidx,:) = obj.Bayes.DistError{phs}(phsidx,:);
        obj.Bayes.bestModel0(phsidx) = obj.Bayes.bestModel{phs}(phsidx);
        nbins = max(floor(obj.data.es.(varname)))+1;
        [obj.Bayes.X,~] = normalise1var(obj.data.es.(varname), nbins);
    end
end

toc

end