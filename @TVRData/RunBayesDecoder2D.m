function obj = RunBayesDecoder2D(obj, varname, predictorname, varargin)
if nargin < 2
    strlistvarname = {'trajPercent','distTrav'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'Contrast for training', 'InitialValue', 1);
    obj.Bayes.varname = strlistvarname{varnamesel};
else
    obj.Bayes.varname = varname;
end
if nargin < 3
    strlistpredictor = {'spikeTrain'};
    [predictorsel,ok] = listdlg('ListString',strlistpredictor, 'Name', 'Contrast for training', 'InitialValue', 1);
    obj.Bayes.predictorname = strlistpredictor{predictorsel};
else
    obj.Bayes.predictorname = predictorname;
end
if nargin < 4
    [nspdbins, smth_spd] = getOptSpdParams(obj);
    Latcorrection = getOptLatParams(obj);
    
    %dialog to ask for training and test sets
    if numel(obj.SubsetVal.contrast) > 1
        [contsel,ok] = listdlg('ListString',strsplit(num2str(obj.SubsetVal.contrast)), 'Name', 'Contrast for training', 'InitialValue', find(obj.SubsetVal.contrast == mode(obj.data.es.contrast)));
        obj.Bayes.TrainContrast = contsel;
    else
        obj.Bayes.TrainContrast = 1;
    end
    if numel(obj.SubsetVal.gain) > 1
        [gainsel,ok] = listdlg('ListString',strsplit(num2str(obj.SubsetVal.gain)), 'Name', 'Gain for training', 'InitialValue', find(obj.SubsetVal.gain == mode(obj.data.es.gain)));
        obj.Bayes.TrainGain = gainsel;
    else
        obj.Bayes.TrainGain = 1;
    end
    if numel(obj.SubsetVal.roomlength) > 1
        [roomsel,ok] = listdlg('ListString',strsplit(num2str(obj.SubsetVal.roomlength)), 'Name', 'Room length for training', 'InitialValue', find(obj.SubsetVal.roomlength == mode(obj.data.es.roomLength)));
        obj.Bayes.TrainRoomlength = roomsel;
    else
        obj.Bayes.TrainRoomlength = 1;
    end
    if numel(obj.SubsetVal.outcome) > 1
        [outcomesel,ok] = listdlg('ListString',strsplit(num2str(obj.SubsetVal.outcome)), 'Name', 'Outcome for training', 'InitialValue', find(obj.SubsetVal.outcome == 2));
        obj.Bayes.TrainOutcome = outcomesel;%3;%
    else
        obj.Bayes.TrainOutcome = 1;%3;%
    end
    prompt = {'Use lookup table';'Window size (ms)';'Spatial smth (%)';'# of X bins';'How many theta phase bins';'theta channel#';'How many speed bins';'How many eye bins';'Speed smooth win. (ms)';'Latency shift (ms)';'alpha';'delta (ms)';'How many runs';'Error threshold';'# folds';'Optimal smoothing?'};
    dlg_title = 'Decoder Parameters';
    num_lines = 1;
    def = {'0';'150';'4';'100';'1';'34';'3';'3';'150';'0';'0';'0';'1';'30';'20';'0'};
    nruns = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(nruns)
        obj.Bayes.Flookuptable = str2num(nruns{1});
        obj.Bayes.Tsmth_win = str2num(nruns{2});
        obj.Bayes.Xsmth_win = str2num(nruns{3});
        obj.Bayes.numBinsX = str2num(nruns{4});
        obj.Bayes.nthetaphsbins = str2num(nruns{5});
        obj.Bayes.thetaChannel = str2num(nruns{6});
        obj.Bayes.nspdbins = str2num(nruns{7});
        obj.Bayes.neyebins = str2num(nruns{8});
        obj.Bayes.smth_spd = str2num(nruns{9});
        obj.Bayes.latcorrection = str2num(nruns{10});
        obj.Bayes.alpha = str2num(nruns{11});
        obj.Bayes.delta = str2num(nruns{12});
        obj.Bayes.nRuns = str2num(nruns{13});
        obj.Bayes.goodidx_th = str2num(nruns{14});
        obj.Bayes.kfold = str2num(nruns{15});
        obj.Bayes.FoptiSmooth = boolean(str2num(nruns{16}));
    else
        obj.Bayes.Tsmth_win = 50;
        obj.Bayes.Flookuptable = false;
        obj.Bayes.Xsmth_win = 1;
        obj.Bayes.numBinsX = 100;
        obj.Bayes.nthetaphsbins = 6;
        obj.Bayes.thetaChannel = 34;
        obj.Bayes.nspdbins = 3;
        obj.Bayes.neyebins = 3;
        obj.Bayes.smth_spd = 800;
        obj.Bayes.latcorrection = 0;
        obj.Bayes.alpha = 0;
        obj.Bayes.delta = 0;
        obj.Bayes.nRuns = 1;
        obj.Bayes.goodidx_th = inf;
        obj.Bayes.kfold = 20;
        obj.Bayes.FoptiSmooth = true;
    end
    
    prompt = {'speed threshold';'good clusters';'unsorted clusters';'MUAs';'max. Spike rate';'Significance thresh. (zscore)';'min Spatial selectivity index'};
    dlg_title = 'Decoder Parameters';
    num_lines = 1;
    def = {'5';'1';'1';'1';'+inf';'-inf';'-inf'};
    nruns = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(nruns)
        obj.Bayes.speed_th = str2num(nruns{1});
        obj.Bayes.FGoodcluster = logical(str2num(nruns{2}));
        obj.Bayes.FUnsortedcluster = logical(str2num(nruns{3}));
        obj.Bayes.FMUAcluster = logical(str2num(nruns{4}));
        obj.Bayes.maxRate = str2num(nruns{5});
        obj.Bayes.zth = str2num(nruns{6});
        obj.Bayes.SSImin = str2num(nruns{7});
    else
        obj.Bayes.speed_th = obj.SpeedThresh;
        obj.Bayes.FGoodcluster = 1;
        obj.Bayes.FUnsortedcluster = 1;
        obj.Bayes.FMUAcluster = 1;
        obj.Bayes.maxRate = 8;
        obj.Bayes.zth = 1.96;
        obj.Bayes.SSImin = 0;
    end
    
    obj.Bayes.type = 'mean';
    obj.Bayes.ProbeID = [1 2];
else
    pnames = {'train_contrast' 'train_gain' 'train_roomlength' 'train_outcome' 'type' 'Flookuptable' 'Tsmth_win' 'Xsmth_win' 'Ysmth_win' 'numBinsX' 'numBinsY' 'nthetaphsbins' 'thetaChannel' 'nspdbins' 'neyebins' 'smth_spd' 'latcorrection' 'alpha' 'delta' 'nruns' 'error_th' 'kfold' 'FoptiSmooth' 'speed_th' 'FGoodcluster' 'FUnsortedcluster' 'FMUAcluster' 'maxRate' 'zth' 'SSImin' 'ProbeID'};
    dflts  = {obj.SubsetVal.contrast obj.SubsetVal.gain obj.SubsetVal.roomlength obj.SubsetVal.outcome 'mean' false 150 4 3 100 10 6 34 3 3 150 0 0 0 1 inf 20 0 5 1 1 1 8 1.96 0 [1 2]};
    [obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome , obj.Bayes.type, obj.Bayes.Flookuptable, obj.Bayes.Tsmth_win, obj.Bayes.Xsmth_win, obj.Bayes.Ysmth_win, obj.Bayes.numBinsX, obj.Bayes.numBinsY, obj.Bayes.nthetaphsbins, obj.Bayes.thetaChannel, obj.Bayes.nspdbins, obj.Bayes.neyebins, obj.Bayes.smth_spd,...
     obj.Bayes.latcorrection, obj.Bayes.alpha, obj.Bayes.delta, obj.Bayes.nRuns, obj.Bayes.goodidx_th, obj.Bayes.kfold, obj.Bayes.FoptiSmooth, obj.Bayes.speed_th, obj.Bayes.FGoodcluster, obj.Bayes.FUnsortedcluster, obj.Bayes.FMUAcluster, obj.Bayes.maxRate, obj.Bayes.zth, obj.Bayes.SSImin, obj.Bayes.ProbeID] = internal.stats.parseArgs(pnames,dflts,varargin{:});     
end
obj.Bayes.sampleRate = mean(1./obj.data.es.sampleSize);
tic

cbase = find(obj.SubsetVal.contrast == mode(obj.data.es.contrast));
gbase = find(obj.SubsetVal.gain == mode(obj.data.es.gain));
rbase = find(obj.SubsetVal.roomlength == mode(obj.data.es.roomLength));
obase = find(obj.SubsetVal.outcome == 2);

SSIlim = [obj.Bayes.SSImin inf];
Ratelim = [0 obj.Bayes.maxRate];
SEfactor = 1;%1/3;
SEmaxlim = size(obj.CellInfo.field{cbase, gbase, rbase, obase},2)*SEfactor;
SSIrangeidx = (obj.CellInfo.SSI{cbase, gbase, rbase, obase} >= SSIlim(1) &...
               obj.CellInfo.SSI{cbase, gbase, rbase, obase} <= SSIlim(2));
raterangeidx = (obj.CellInfo.fieldMin{cbase, gbase, rbase, obase} >= Ratelim(1) &...
                obj.CellInfo.fieldMin{cbase, gbase, rbase, obase} <= Ratelim(2));
goodunits = (double(obj.Bayes.FGoodcluster) * double(obj.CellInfo.Goodcluster) | double(obj.Bayes.FUnsortedcluster) * double(obj.CellInfo.Unsortedcluster) | double(obj.Bayes.FMUAcluster) * double(obj.CellInfo.MUAcluster));
stableplacefields = (obj.CellInfo.fieldPosSE{cbase, gbase, rbase, obase} < SEmaxlim);
CA1cells = obj.CellInfo.Probe == 1;
CA1interneurons = obj.CellInfo.Finterneuron;
[maxZ,~] = max(obj.CellInfo.fieldZ{cbase, gbase, rbase, obase},[],2);
signiplacefields = true(size(stableplacefields));%(abs(maxZ') > obj.Bayes.zth);

distri_Th = 0.99;
SinfoperSpike = obj.CellInfo.SpatialInfoPerSpike{cbase, gbase, rbase, obase};
% SinfoperSpikeZscore = (SinfoperSpike - quantile(obj.CellInfo.SpatialInfoPerSpikeRef{cbase, gbase, rbase, obase},distri_Th,2)')./quantile(obj.CellInfo.SpatialInfoPerSpikeRef{cbase, gbase, rbase, obase},distri_Th,2)';
% signiSinfoperSpike = SinfoperSpikeZscore>1;
signiSinfoperSpike = true(size(SinfoperSpike));
stableplacefields = true(size(stableplacefields));


nbProbe = numel(unique(obj.CellInfo.Probe));
if numel(obj.Bayes.nspdbins) < nbProbe
    obj.Bayes.nspdbins = repmat(obj.Bayes.nspdbins,[nbProbe 1]);
end
if numel(obj.Bayes.neyebins) < nbProbe
    obj.Bayes.neyebins = repmat(obj.Bayes.neyebins,[nbProbe 1]);
end
if numel(obj.Bayes.smth_spd) < nbProbe
    obj.Bayes.smth_spd = repmat(obj.Bayes.smth_spd,[nbProbe 1]);
end
if numel(obj.Bayes.latcorrection) < nbProbe
    obj.Bayes.latcorrection = repmat(obj.Bayes.latcorrection,[nbProbe 1]);
end
if numel(obj.Bayes.alpha) < nbProbe
    obj.Bayes.alpha = repmat(obj.Bayes.alpha,[nbProbe 1]);
end
if numel(obj.Bayes.delta) < nbProbe
    obj.Bayes.delta = repmat(obj.Bayes.delta,[nbProbe 1]);
end
obj.Bayes.Posterior = cell(nbProbe, max(obj.Bayes.nspdbins), max(obj.Bayes.nthetaphsbins));
obj.Bayes.PosError = cell(nbProbe, max(obj.Bayes.nspdbins), max(obj.Bayes.nthetaphsbins));

obj.Bayes.LFPphase = cell(nbProbe,1);
obj.Bayes.LFPphase = cell(nbProbe,1);
% obj.Bayes.LFPphase = obj.data.es.LFPbinphase;%obj.data.es.LFPphase(:,min(end,obj.Bayes.thetaChannel));

nbins = obj.Bayes.numBinsX;
smthtype = 'boxcar_centered';
%temp
wintype = 'boxcarsum_centered';%'gaussian_1';%'boxcarsum';
obj.Bayes.X0 = cell(1,2);
obj.Bayes.Posterior0 = cell(nbProbe,2);
obj.Bayes.PosError0 = cell(nbProbe,2);

for iprobe = 1:nbProbe
    if iprobe == 1 %CA1
        obj.Bayes.DecCellidx{iprobe} = signiSinfoperSpike & raterangeidx & stableplacefields & signiplacefields & SSIrangeidx & goodunits & CA1cells & ~CA1interneurons;
    elseif iprobe == 2 %V1
        obj.Bayes.DecCellidx{iprobe} =   signiSinfoperSpike & stableplacefields & signiplacefields & SSIrangeidx & goodunits & ~CA1cells;
    end
end

for iprobe = 1:nbProbe
    if ((iprobe==1 && sum(obj.CellInfo.Probe == 1)>0) || (iprobe==2 && sum(obj.CellInfo.Probe == 2) > 0)) && ismember(iprobe,obj.Bayes.ProbeID)
        sampleRate = obj.Bayes.sampleRate;
        
        if strcmp(obj.Bayes.varname,'trajPercent')
            maxX = 100;
            Fcircular = true;
        else
            maxX = max(obj.data.es.(obj.Bayes.varname));
            Fcircular = false;
        end
        X = getBayesVar(obj.data.es.(obj.Bayes.varname),obj.Bayes.latcorrection(iprobe),sampleRate,obj.Bayes.Tsmth_win,smthtype,Fcircular,maxX);
        [X, ~] = normalise1var(X, obj.Bayes.numBinsX);
        obj.Bayes.X = X(:);
        
        idxref = obj.getSubsets(obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome,obj.Bayes.speed_th);
        obj.Bayes.varnameXy = 'ballspeed';
        obj.Bayes.Ysmth_win = 1;
        Xy = getBayesVar(obj.data.es.(obj.Bayes.varnameXy),obj.Bayes.latcorrection(iprobe),sampleRate,obj.Bayes.Tsmth_win,smthtype,false,maxX);
%         [Xy, ~] = normalise1var(Xy, nSpeedranges);
        
        Xyprofile = fast1Dmap(X(idxref),Xy(idxref),1,1,1/(obj.Bayes.Xsmth_win/obj.Bayes.numBinsX),Fcircular);
        Xy = Xy./(Xyprofile(X));
        Xy(Xy < 0.5) = 0.4;
        Xy(Xy > 1.5) = 1.6;
        Ybinsize = 0.1;
        nSpeedranges = floor((max(Xy(idxref))-min(Xy(idxref)))/Ybinsize)+1;
        [Xy, ~] = normalise1var(Xy, nSpeedranges);
    
%         nSpeedranges = 10;
%         spdquantilelim = zeros(obj.Bayes.numBinsX,2);        
%         Xybinned = NaN(size(Xy));
%         if nSpeedranges > 1
%             for spd = 1:nSpeedranges
%                 for xx = 1:obj.Bayes.numBinsX
%                     spdquantilelim(xx,1) = quantile(Xy(idxref & X == xx),max(0,(spd-1)/nSpeedranges));
%                     spdquantilelim(xx,2) = quantile(Xy(idxref & X == xx),min(1,(spd)/nSpeedranges));
%                 end
%                 Xybinned(Xy >= spdquantilelim(X,1) & Xy < spdquantilelim(X,2)) = spd;
%                 if spd == 1
%                     Xybinned(Xy <= spdquantilelim(X,1)) = spd;
%                 end
%                 if spd == nSpeedranges
%                     Xybinned(Xy >= spdquantilelim(X,2)) = spd;
%                 end
%             end        
%         else
%             Xybinned = ones(size(Xy));
%         end
%         Xy = Xybinned;

        obj.Bayes.X = [X(:) Xy];        
        Y = obj.data.es.(obj.Bayes.predictorname);%circshift(obj.data.es.(obj.Bayes.predictorname),[-round(obj.Bayes.latcorrection(iprobe)/(1000/sampleRate)) 0]);%
        
        win = round(obj.Bayes.Tsmth_win/(1000/sampleRate));
        T = obj.data.es.sampleSize.*win;
        
        obj.Bayes.paramvarname = 'ballspeed';%'trajspeed';%'ballspeed';%'eyeXpos';%'ballspeed';%'eyeXpos';%'rand';%'ballspeed';%
        disp(['will use ' obj.Bayes.paramvarname]);
        if strcmp(obj.Bayes.varname,'distTrav')
            X = X / 2;
            X = min(X,100);
            speed = NaN(size(obj.data.es.ballspeed));
            speed(~isnan(obj.data.es.ballspeed)) = smthInTime(obj.data.es.ballspeed(~isnan(obj.data.es.ballspeed)), sampleRate, obj.Bayes.smth_spd(iprobe), 'same', [], 'boxcar');
        elseif strcmp(obj.Bayes.paramvarname,'rand')
            speed = NaN(size(obj.data.es.trajspeed));
            speed(~isnan(obj.data.es.trajspeed)) = smthInTime(randn(sum(~isnan(obj.data.es.trajspeed)),1), sampleRate, obj.Bayes.smth_spd(iprobe), 'same', [], 'boxcar');%'boxcar');
        else
            speed = NaN(size(obj.data.es.(obj.Bayes.paramvarname)));
            speed(~isnan(obj.data.es.(obj.Bayes.paramvarname))) = smthInTime(obj.data.es.(obj.Bayes.paramvarname)(~isnan(obj.data.es.(obj.Bayes.paramvarname))), sampleRate, obj.Bayes.smth_spd(iprobe), 'same', [], smthtype);%'boxcar');
        end
        
        eyeX = NaN(size(obj.data.es.eyeXpos));
        eyeX(~isnan(obj.data.es.eyeXpos)) = smthInTime(obj.data.es.eyeXpos(~isnan(obj.data.es.eyeXpos)), sampleRate, obj.Bayes.smth_spd(iprobe), 'same', [], smthtype);%'boxcar');
        
        if iprobe == 1 %CA1
            if isfield(obj.data.es,'LFPphase')
                smthLFPphase = smthInTime(unwrap(obj.data.es.LFPphase(:,min(end,obj.Bayes.thetaChannel))/360*2*pi)*360/(2*pi), sampleRate, obj.Bayes.Tsmth_win, 'same', [], smthtype);
            else
                smthLFPphase = zeros(size(obj.data.es.traj));
            end
        elseif iprobe == 2 %V1
            if isfield(obj.data.es,'LFPphaseV1')
                smthLFPphase = smthInTime(unwrap(obj.data.es.LFPphaseV1(:,min(end,obj.Bayes.thetaChannel))/360*2*pi)*360/(2*pi), sampleRate, obj.Bayes.Tsmth_win, 'same', [], smthtype);
            elseif isfield(obj.data.es,'LFPphase')
                smthLFPphase = smthInTime(unwrap(obj.data.es.LFPphase(:,min(end,obj.Bayes.thetaChannel))/360*2*pi)*360/(2*pi), sampleRate, obj.Bayes.Tsmth_win, 'same', [], smthtype);
            else
                smthLFPphase = zeros(size(obj.data.es.traj));
            end
        end
        Y = Y(:,obj.Bayes.DecCellidx{iprobe});%obj.data.es.CellInfo.MUAcluster obj.data.es.CellInfo.Goodcluster
        disp([num2str(sum(obj.Bayes.DecCellidx{iprobe})) ' cells used for decoding']);
        
        nbCells = size(Y,2);
        Yfilt = zeros(size(Y));
        for icell = 1:nbCells
            Yfilt(:,icell) = smthInTime(Y(:,icell), sampleRate, obj.Bayes.Tsmth_win, 'same', [], wintype);%, 'boxcarsum');
        end
        Y = Yfilt;
        T = smthInTime(T, sampleRate, obj.Bayes.Tsmth_win, 'same', [], wintype);%, 'boxcarsum');
        
        
        obj.Bayes.LFPphase{iprobe} = mod(round(smthLFPphase),360);
        obj.Bayes.LFPphase{iprobe} = obj.Bayes.LFPphase{iprobe}(:);
        
        
        
        itraj = floor(X/(max(round(X))/obj.Bayes.numBinsX))+1;
        ntrajbins = max(itraj);
        spdquantilelim = zeros(ntrajbins,2);        
        obj.Bayes.Spdbin = NaN(size(obj.data.es.trajspeed));
        if obj.Bayes.nspdbins(iprobe) > 1
            for spd = 1:obj.Bayes.nspdbins(iprobe)
                for xx = 1:ntrajbins
                    spdquantilelim(xx,1) = quantile(speed(idxref & itraj == xx),max(0,(spd-1)/obj.Bayes.nspdbins(iprobe)));
                    spdquantilelim(xx,2) = quantile(speed(idxref & itraj == xx),min(1,(spd)/obj.Bayes.nspdbins(iprobe)));
                end
                obj.Bayes.Spdbin(speed >= spdquantilelim(itraj,1) & speed < spdquantilelim(itraj,2)) = spd;
                if spd == 1
                    obj.Bayes.Spdbin(speed <= spdquantilelim(itraj,1)) = spd;
                end
                if spd == obj.Bayes.nspdbins(iprobe)
                    obj.Bayes.Spdbin(speed >= spdquantilelim(itraj,2)) = spd;
                end
            end        
        else
            obj.Bayes.Spdbin = ones(size(obj.Bayes.Spdbin));
        end        
        
        
        eyequantilelim = zeros(ntrajbins,2);        
        obj.Bayes.Eyebin = NaN(size(obj.data.es.eyeXpos));
        if max(eyeX(~isnan(eyeX))) == 0
            obj.Bayes.neyebins(iprobe) = 1;
        end
        if obj.Bayes.neyebins(iprobe) > 1
            for eye = 1:obj.Bayes.neyebins(iprobe)
                for xx = 1:ntrajbins
                    eyequantilelim(xx,1) = quantile(eyeX(idxref & itraj == xx),max(0,(eye-1)/obj.Bayes.neyebins(iprobe)));
                    eyequantilelim(xx,2) = quantile(eyeX(idxref & itraj == xx),min(1,(eye)/obj.Bayes.neyebins(iprobe)));
                end
                obj.Bayes.Eyebin(eyeX >= eyequantilelim(itraj,1) & eyeX < eyequantilelim(itraj,2)) = eye;
                if eye == 1
                    obj.Bayes.Eyebin(eyeX <= eyequantilelim(itraj,1)) = eye;
                end
                if eye == obj.Bayes.nspdbins(iprobe)
                    obj.Bayes.Eyebin(eyeX >= eyequantilelim(itraj,2)) = eye;
                end
            end        
        else
            obj.Bayes.Eyebin = ones(size(obj.Bayes.Eyebin));
        end
        
        
        obj.Bayes.Posterior0{iprobe,1} = zeros(numel(X),obj.Bayes.numBinsX);
        obj.Bayes.PosError0{iprobe,1} = zeros(numel(X),obj.Bayes.numBinsX);
        obj.Bayes.Posterior0{iprobe,2} = zeros(numel(X),nSpeedranges);
        obj.Bayes.PosError0{iprobe,2} = zeros(numel(X),nSpeedranges);
        
        for irun = 1:obj.Bayes.nRuns
            for phs = 1:obj.Bayes.nthetaphsbins
                idxref = obj.getSubsets(obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome,obj.Bayes.speed_th);
                if irun > 1
                    obj.Bayes.Spdbin = obj.Bayes.bestquantile;
                else
                    obj.Bayes.bestquantile = obj.Bayes.Spdbin;
                    obj.Bayes.goodidx = true(size(idxref));
                end
                
                phs0 = (phs-1)*360/obj.Bayes.nthetaphsbins;
                phsval = obj.Bayes.LFPphase{iprobe};%mod(obj.Bayes.LFPphase{iprobe}-180,360);
                phsidx = ismember(phsval,mod(phs0:(phs0 + 360/obj.Bayes.nthetaphsbins)-1,360));
                
                idx = idxref & obj.Bayes.goodidx & phsidx;
                
                obj.Bayes.prediction{iprobe,phs} = zeros(numel(X),2);
                obj.Bayes.X = zeros(numel(X),2);
                obj.Bayes.Posterior{iprobe,phs} = NaN(numel(X),nSpeedranges,obj.Bayes.numBinsX);
                
                disp(['Time Smth win = ' num2str(obj.Bayes.Tsmth_win) ' ms ; Spatial Smth = ' num2str(obj.Bayes.Xsmth_win) '%']);
                nspdbins = max(obj.Bayes.Spdbin);
                neyebins = max(obj.Bayes.Eyebin);
                for spd = 1:nspdbins
                    for eye = 1:neyebins
                        disp(['Bayesian decoder: run #' num2str(irun) ' out of ' num2str(obj.Bayes.nRuns) ' ; phs#' num2str(phs) ' out of ' num2str(obj.Bayes.nthetaphsbins) ' ; spd#' num2str(spd) ' out of ' num2str(obj.Bayes.nspdbins(iprobe)) ' ; eye#' num2str(eye) ' out of ' num2str(obj.Bayes.neyebins(iprobe))]);
                        disp(['Speed thresh. ' num2str(obj.Bayes.speed_th) ' cm/s ; Smth Spd win ' num2str(obj.Bayes.smth_spd(iprobe)) 'ms ; Latency correction ' num2str(obj.Bayes.latcorrection(iprobe)) ' ms; alpha ' num2str(obj.Bayes.alpha(iprobe))]);
                        spdidx = (obj.Bayes.Spdbin == spd);%true(size(idx));%
                        
                        eyeidx = (obj.Bayes.Eyebin == eye);
                        spdidx = spdidx & eyeidx;
                        
                        obj.Bayes.decOder = TbayesDecoder2D;
                        obj.Bayes.decOder.numBinsX = obj.Bayes.numBinsX;
                        obj.Bayes.decOder.numBinsY = nSpeedranges;
                        obj.Bayes.decOder.Fcircular = obj.data.es.CircularMaze;
                        obj.Bayes.decOder.Flookuptable = obj.Bayes.Flookuptable;
                        obj.Bayes.decOder.kfold = obj.Bayes.kfold;
                        obj.Bayes.decOder.FoptiSmooth = obj.Bayes.FoptiSmooth;
                        
                        [obj.Bayes.decOder, pred, Xdec, Posterior, ~, ~] = obj.Bayes.decOder.trainDecoder(X(idx & spdidx), Xy(idx & spdidx), Y(idx & spdidx ,:), T(idx & spdidx ), obj.Bayes.Xsmth_win, obj.Bayes.Ysmth_win);
                        
                        
                        obj.Bayes.prediction{iprobe,phs}(idx & spdidx ,:) = pred;
                        obj.Bayes.X(idx & spdidx ,1) = Xdec(1:min(end,size(pred,1)),1);
                        obj.Bayes.X(idx & spdidx ,2) = Xdec(1:min(end,size(pred,1)),2);
                        obj.Bayes.Posterior{iprobe,phs}(idx & spdidx ,:,:) = Posterior;%nonNormPosterior;%Posterior;
                        
                        [pred, Xdec, Posterior, ~, ~] = obj.Bayes.decOder.predictBayesDecoder(X(~idx & spdidx ,:), Xy(~idx & spdidx ,:), Y(~idx & spdidx ,:), T(~idx & spdidx ), obj.Bayes.type);
                        
                        obj.Bayes.prediction{iprobe,phs}(~idx & spdidx ,:) = pred;
                        obj.Bayes.X(~idx & spdidx ,1) = Xdec(1:min(end,size(pred,1)),1);
                        obj.Bayes.X(~idx & spdidx ,2) = Xdec(1:min(end,size(pred,1)),2);
                        obj.Bayes.Posterior{iprobe,phs}(~idx & spdidx ,:,:) = Posterior;%nonNormPosterior;%Posterior;%
                    end
                end
                obj.Bayes.X0{1} = obj.Bayes.X(:,1);
                obj.Bayes.X0{2} = obj.Bayes.X(:,2);
                obj.Bayes.Posterior0{iprobe,1,phs} = squeeze(nanmean(obj.Bayes.Posterior{iprobe,phs},2));
                obj.Bayes.Posterior0{iprobe,2,phs} = squeeze(nanmean(obj.Bayes.Posterior{iprobe,phs},3));
                PrangeX = size(obj.Bayes.Posterior0{iprobe,1,phs},2);
                PrangeY = size(obj.Bayes.Posterior0{iprobe,2,phs},2);
                Xrange = max(obj.Bayes.X0{1});
                Yrange = max(obj.Bayes.X0{2});
                for xx = 1:Xrange
                    obj.Bayes.PosError0{iprobe,1,phs}(obj.Bayes.X0{1} == xx,:) = circshift(obj.Bayes.Posterior0{iprobe,1,phs}(obj.Bayes.X0{1} == xx,:),floor(PrangeX/2)-xx+1,2);
                end
                for yy = 1:Yrange
                    obj.Bayes.PosError0{iprobe,2,phs}(obj.Bayes.X0{2} == yy,:) = circshift(obj.Bayes.Posterior0{iprobe,2,phs}(obj.Bayes.X0{2} == yy,:),floor(PrangeY/2)-yy+1,2);
                end
            end
            obj.Bayes.time = ((1:size(obj.Bayes.X,1))./obj.Bayes.sampleRate)';
        end
        
        %emptying all those to save memory
        obj.Bayes.decOder = [];
        obj.Bayes.Posterior = [];
        obj.Bayes.PosError = [];
        obj.Bayes.prediction = [];
        obj.Bayes.X = [];
        toc
    else
        for spd = 1:obj.Bayes.nspdbins(iprobe)
            for phs = 1:obj.Bayes.nthetaphsbins
                obj.Bayes.Posterior{iprobe,phs} = NaN(size(obj.Bayes.Posterior{1,phs}));
                obj.Bayes.PosError{iprobe,phs} = NaN(size(obj.Bayes.PosError{1,phs}));
                obj.Bayes.prediction{iprobe,phs} = NaN(size(obj.Bayes.prediction{1,phs}));
            end
        end
        obj.Bayes.LFPphase{iprobe} = NaN(size(obj.Bayes.LFPphase{1}));
        obj.Bayes.Posterior0{iprobe,1} = NaN(size(obj.Bayes.Posterior0{1,1}));
        obj.Bayes.PosError0{iprobe,1} = NaN(size(obj.Bayes.PosError0{1,1}));
        obj.Bayes.Posterior0{iprobe,2} = NaN(size(obj.Bayes.Posterior0{1,1}));
        obj.Bayes.PosError0{iprobe,2} = NaN(size(obj.Bayes.PosError0{1,1}));
    end
end
end