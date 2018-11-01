function obj = RunBayesDecoder(obj, varname, predictorname, varargin)
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
        obj.Bayes.numBins = str2num(nruns{4});
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
        obj.Bayes.numBins = 100;
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
    pnames = {'train_contrast' 'train_gain' 'train_roomlength' 'train_outcome' 'type' 'Flookuptable' 'Tsmth_win' 'Xsmth_win' 'numBins' 'nthetaphsbins' 'thetaChannel' 'nspdbins' 'neyebins' 'smth_spd' 'latcorrection' 'alpha' 'delta' 'nruns' 'error_th' 'kfold' 'FoptiSmooth' 'speed_th' 'FGoodcluster' 'FUnsortedcluster' 'FMUAcluster' 'maxRate' 'zth' 'SSImin' 'ProbeID'};
    dflts  = {obj.SubsetVal.contrast obj.SubsetVal.gain obj.SubsetVal.roomlength obj.SubsetVal.outcome 'mean' false 150 4 100 6 34 3 3 150 0 0 0 1 inf 20 0 5 1 1 1 8 1.96 0 [1 2]};
    [obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome , obj.Bayes.type, obj.Bayes.Flookuptable, obj.Bayes.Tsmth_win, obj.Bayes.Xsmth_win, obj.Bayes.numBins, obj.Bayes.nthetaphsbins, obj.Bayes.thetaChannel, obj.Bayes.nspdbins, obj.Bayes.neyebins, obj.Bayes.smth_spd,...
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
obj.Bayes.DistError = cell(nbProbe, max(obj.Bayes.nspdbins),max(obj.Bayes.nthetaphsbins));
obj.Bayes.prediction = cell(nbProbe, max(obj.Bayes.nspdbins),max(obj.Bayes.nthetaphsbins));

obj.Bayes.LFPphase = cell(nbProbe,1);
obj.Bayes.LFPphase = cell(nbProbe,1);
% obj.Bayes.LFPphase = obj.data.es.LFPbinphase;%obj.data.es.LFPphase(:,min(end,obj.Bayes.thetaChannel));

nbins = obj.Bayes.numBins;
smthtype = 'boxcar_centered';
%temp
wintype = 'boxcarsum_centered';%'gaussian_1';%'boxcarsum';
obj.Bayes.Posterior0 = cell(nbProbe,1);
obj.Bayes.PosError0 = cell(nbProbe,1);
obj.Bayes.DistError0 = cell(nbProbe,1);
obj.Bayes.DecCellidx = cell(nbProbe,1);

for iprobe = 1:nbProbe
    if iprobe == 1 %CA1
        obj.Bayes.DecCellidx{iprobe} = signiSinfoperSpike & raterangeidx & stableplacefields & signiplacefields & SSIrangeidx & goodunits & CA1cells;
    elseif iprobe == 2 %V1
        obj.Bayes.DecCellidx{iprobe} =   signiSinfoperSpike & stableplacefields & signiplacefields & SSIrangeidx & goodunits & ~CA1cells;
    end
end

for iprobe = 1:nbProbe
    if ((iprobe==1 && sum(obj.CellInfo.Probe == 1)>0) || (iprobe==2 && sum(obj.CellInfo.Probe == 2) > 0)) && ismember(iprobe,obj.Bayes.ProbeID)
        sampleRate = obj.Bayes.sampleRate;
%         X = obj.data.es.(obj.Bayes.varname);

%         Xtemp = unwrap(obj.data.es.(obj.Bayes.varname)/max(obj.data.es.(obj.Bayes.varname))*2*pi)*max(obj.data.es.(obj.Bayes.varname))/(2*pi);
%         Xtemp = smthInTime(Xtemp, sampleRate, obj.Bayes.Tsmth_win, 'same', [], 'boxcar_centered');
%         Xtemp = mod(Xtemp,max(obj.data.es.(obj.Bayes.varname)));
        
%         X = mod(Xtemp,max(obj.data.es.(obj.Bayes.varname)));
%         lat = round(obj.Bayes.latcorrection(iprobe)/(1000/sampleRate));
%         Rundist = [Xtemp(lat+1:end) - Xtemp(1:end-lat); zeros(lat,1)]./(obj.data.es.gain/mode(obj.data.es.gain));
        
        %temp
        %%%in process: introduce dependence over past time series
        X = VisDistXmodel(obj.data.es.(obj.Bayes.varname),obj.Bayes.delta(iprobe),obj.Bayes.alpha(iprobe),obj.Bayes.latcorrection(iprobe),obj.data.es,obj.Bayes.Tsmth_win,'boxcar_centered',obj.Bayes.numBins);
%         X = KalmanX2model(obj.data.es.(obj.Bayes.varname),1,1,0,obj.data.es,obj.Bayes.Tsmth_win,'boxcar_centered',obj.Bayes.numBins);
        obj.Bayes.X0 = X;
%         X = VisDistXmodel(obj.data.es.(obj.Bayes.varname),obj.Bayes.latcorrection(iprobe),obj.Bayes.alpha(iprobe),obj.data.es,obj.Bayes.Tsmth_win,'boxcar_centered',obj.Bayes.numBins);
        
        %temp
%         X = smthInTime(Xtemp, sampleRate, obj.Bayes.Tsmth_win, 'same', [], smthtype);%
%         X = mod(X,max(obj.data.es.(obj.Bayes.varname)));%obj.data.es.(obj.Bayes.varname);%
        
        Y = obj.data.es.(obj.Bayes.predictorname);%circshift(obj.data.es.(obj.Bayes.predictorname),[-round(obj.Bayes.latcorrection(iprobe)/(1000/sampleRate)) 0]);%
        
        win = round(obj.Bayes.Tsmth_win/(1000/sampleRate));
        T = obj.data.es.sampleSize.*win;
        
        %temp
%         obj.Bayes.smth_spd(iprobe) = obj.Bayes.Tsmth_win;%NTaumemory*1/sampleRate*1000;%obj.Bayes.latcorrection(iprobe);
        
        idxref = obj.getSubsets(obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome,obj.Bayes.speed_th);
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
                        
%             Xunwraptraj = unwrap(X/max(X)*2*pi)*max(X)/(2*pi);
%             gaintrial = obj.data.es.gain./mode(obj.data.es.gain);
%             speed = [diff(Xunwraptraj);0]./gaintrial;
%             speed(~isnan(speed)) = smthInTime(speed(~isnan(speed)), sampleRate, obj.Bayes.smth_spd(iprobe), 'same', [], smthtype);%'boxcar');
            
%             speed(~isnan(obj.data.es.(obj.Bayes.paramvarname))) = smthInTime([0;diff(obj.data.es.(obj.Bayes.paramvarname)(~isnan(obj.data.es.(obj.Bayes.paramvarname))))], sampleRate, obj.Bayes.smth_spd(iprobe), 'same', [], smthtype);%'boxcar');
            %temp
%             speed = NaN(size(obj.data.es.(obj.Bayes.paramvarname)));
%             speed(~isnan(obj.data.es.(obj.Bayes.paramvarname))) = smthInTime(obj.data.es.(obj.Bayes.paramvarname)(~isnan(obj.data.es.(obj.Bayes.paramvarname))), sampleRate, obj.Bayes.smth_spd(iprobe), 'same', [], 'boxcar');%'boxcar');
%             speedtemp = speed(~isnan(obj.data.es.(obj.Bayes.paramvarname)));
%             for k = 1:NTaumemory
%                 Sconv = conv(Xmemory,speedtemp(k:NTaumemory:end),'full');
%                 speedtemp(k:NTaumemory:end) = Sconv(1:numel(speedtemp(k:NTaumemory:end)));
%             end
%             speed(~isnan(obj.data.es.(obj.Bayes.paramvarname))) = speedtemp;
        end
        
        eyeX = NaN(size(obj.data.es.eyeXpos));
        eyeX(~isnan(obj.data.es.eyeXpos)) = smthInTime(obj.data.es.eyeXpos(~isnan(obj.data.es.eyeXpos)), sampleRate, obj.Bayes.smth_spd(iprobe), 'same', [], smthtype);%'boxcar');
        
%         speed =  Rundist;
        
        %temp
        
        if iprobe == 1 %CA1
            obj.Bayes.DecCellidx{iprobe} = signiSinfoperSpike & raterangeidx & stableplacefields & signiplacefields & SSIrangeidx & goodunits & CA1cells;%true(size(signiSinfoperSpike));%% & SSIrangeidx;% & goodunits & stableplacefields & signiplacefields;%& SSIrangeidx;%raterangeidx &
            if isfield(obj.data.es,'LFPphase')
                smthLFPphase = smthInTime(unwrap(obj.data.es.LFPphase(:,min(end,obj.Bayes.thetaChannel))/360*2*pi)*360/(2*pi), sampleRate, obj.Bayes.Tsmth_win, 'same', [], smthtype);
            else
                smthLFPphase = zeros(size(obj.data.es.traj));
            end
        elseif iprobe == 2 %V1
            obj.Bayes.DecCellidx{iprobe} =   signiSinfoperSpike & stableplacefields & signiplacefields & SSIrangeidx & goodunits & ~CA1cells;% & SSIrangeidx;% & goodunits & stableplacefields & signiplacefields;%& SSIrangeidx;%raterangeidx &
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
%         Y = double(round(Y) > 0);
        T = smthInTime(T, sampleRate, obj.Bayes.Tsmth_win, 'same', [], wintype);%, 'boxcarsum');
        
        
        obj.Bayes.LFPphase{iprobe} = mod(round(smthLFPphase),360);
        obj.Bayes.LFPphase{iprobe} = obj.Bayes.LFPphase{iprobe}(:);
%         obj.Bayes.LFPphase{iprobe} = circshift(obj.Bayes.LFPphase{iprobe},[-round(obj.Bayes.latcorrection(iprobe)/(1000/sampleRate)) 0]);%

        % obj.Bayes.LFPphase{iprobe} = obj.data.es.LFPbinphase;%obj.data.es.LFPphase{iprobe}(:,min(end,obj.Bayes.thetaChannel));
        
        
        
        itraj = floor(X/(max(round(X))/obj.Bayes.numBins))+1;
        ntrajbins = max(itraj);
        spdquantilelim = zeros(ntrajbins,2);        
        obj.Bayes.Spdbin = NaN(size(obj.data.es.trajspeed));
        if obj.Bayes.nspdbins(iprobe) > 1
            for spd = 1:obj.Bayes.nspdbins(iprobe)
                for xx = 1:ntrajbins
                    spdquantilelim(xx,1) = quantile(speed(idxref & itraj == xx),max(0,(spd-1)/obj.Bayes.nspdbins(iprobe)));
                    spdquantilelim(xx,2) = quantile(speed(idxref & itraj == xx),min(1,(spd)/obj.Bayes.nspdbins(iprobe)));

                    %                     spdquantilelim(xx,1) = quantile(visualspeed(obj.Subset.all & obj.data.es.outcome == 2 & itraj == xx),max(0,(spd-2)/obj.Bayes.nspdbins(iprobe)));
                    %                     spdquantilelim(xx,2) = quantile(visualspeed(obj.Subset.all & obj.data.es.outcome == 2 & itraj == xx),min(1,(spd+2)/obj.Bayes.nspdbins(iprobe)));

                    %                     spdquantilelim(xx,1) = quantile(obj.data.es.smthTrajSpd(obj.Subset.all & itraj == xx),(spd-1)/obj.Bayes.nspdbins);
                    %                     spdquantilelim(xx,2) = quantile(obj.data.es.smthTrajSpd(obj.Subset.all & itraj == xx)/obj.Bayes.nspdbins);
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

                    %                     spdquantilelim(xx,1) = quantile(visualspeed(obj.Subset.all & obj.data.es.outcome == 2 & itraj == xx),max(0,(spd-2)/obj.Bayes.nspdbins(iprobe)));
                    %                     spdquantilelim(xx,2) = quantile(visualspeed(obj.Subset.all & obj.data.es.outcome == 2 & itraj == xx),min(1,(spd+2)/obj.Bayes.nspdbins(iprobe)));

                    %                     spdquantilelim(xx,1) = quantile(obj.data.es.smthTrajSpd(obj.Subset.all & itraj == xx),(spd-1)/obj.Bayes.nspdbins);
                    %                     spdquantilelim(xx,2) = quantile(obj.data.es.smthTrajSpd(obj.Subset.all & itraj == xx)/obj.Bayes.nspdbins);
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
        
        nbins = obj.Bayes.numBins;        
        
        obj.Bayes.Posterior0{iprobe} = zeros(numel(X),nbins);
        obj.Bayes.PosError0{iprobe} = zeros(numel(X),nbins);
        obj.Bayes.DistError0{iprobe} = zeros(numel(X),nbins);
        
        for irun = 1:obj.Bayes.nRuns
            for phs = 1:obj.Bayes.nthetaphsbins
                    idxref = obj.getSubsets(obj.Bayes.TrainContrast, obj.Bayes.TrainGain, obj.Bayes.TrainRoomlength, obj.Bayes.TrainOutcome,obj.Bayes.speed_th);
                    if irun > 1
                        obj.Bayes.Spdbin = obj.Bayes.bestquantile;
                        
%                         maxtol = 1;
%                         [obj.Bayes.goodidx, Xpred] = getcleandecidx(obj.Bayes.Posterior0{iprobe}, obj.Bayes.X,obj.Bayes.goodidx_th,maxtol,obj.data.es.trialID);
%                         
% %                         Xpredunwrapped = (unwrap(Xpred/max(Xpred)*2*pi)*max(Xpred)/(2*pi));
% %                         Xpreddiff = [0;diff(Xpredunwrapped)];
% %                         Xpreddiff(~isnan(Xpreddiff)) = smthInTime(Xpreddiff(~isnan(Xpreddiff)), sampleRate, obj.Bayes.smth_spd(iprobe), 'same', [], smthtype);
% %                         Xpred = Xpreddiff;
%                         
%                         obj.Bayes.goodidx = true(size(idxref));
%                         for spd = 1:obj.Bayes.nspdbins(iprobe)
%                             for xx = 1:ntrajbins
%                                 spdquantilelim(xx,1) = quantile(Xpred(idxref & itraj == xx),max(0,(spd-1)/obj.Bayes.nspdbins(iprobe)));
%                                 spdquantilelim(xx,2) = quantile(Xpred(idxref & itraj == xx),min(1,(spd)/obj.Bayes.nspdbins(iprobe)));
%                                 
%                                 %                     spdquantilelim(xx,1) = quantile(visualspeed(obj.Subset.all & obj.data.es.outcome == 2 & itraj == xx),max(0,(spd-2)/obj.Bayes.nspdbins(iprobe)));
%                                 %                     spdquantilelim(xx,2) = quantile(visualspeed(obj.Subset.all & obj.data.es.outcome == 2 & itraj == xx),min(1,(spd+2)/obj.Bayes.nspdbins(iprobe)));
%                                 
%                                 %                     spdquantilelim(xx,1) = quantile(obj.data.es.smthTrajSpd(obj.Subset.all & itraj == xx),(spd-1)/obj.Bayes.nspdbins);
%                                 %                     spdquantilelim(xx,2) = quantile(obj.data.es.smthTrajSpd(obj.Subset.all & itraj == xx)/obj.Bayes.nspdbins);
%                             end
%                             obj.Bayes.Spdbin(Xpred >= spdquantilelim(itraj,1) & Xpred < spdquantilelim(itraj,2)) = spd;
%                             if spd == 1
%                                 obj.Bayes.Spdbin(Xpred < spdquantilelim(itraj,1)) = spd;
%                             end
%                             if spd == obj.Bayes.nspdbins(iprobe)
%                                 obj.Bayes.Spdbin(Xpred >= spdquantilelim(itraj,2)) = spd;
%                             end
%                         end
                    else
%                         obj.Bayes.Spdbin = ones(size(obj.Bayes.Spdbin));
                        obj.Bayes.bestquantile = obj.Bayes.Spdbin;
                        obj.Bayes.goodidx = true(size(idxref));
                    end
                    
                    phs0 = (phs-1)*360/obj.Bayes.nthetaphsbins;
                    phsval = obj.Bayes.LFPphase{iprobe};%mod(obj.Bayes.LFPphase{iprobe}-180,360);
                    phsidx = ismember(phsval,mod(phs0:(phs0 + 360/obj.Bayes.nthetaphsbins)-1,360));
                    
                    idx = idxref & obj.Bayes.goodidx & phsidx;
                    
                    obj.Bayes.prediction{iprobe,phs} = zeros(numel(X),1);
                    obj.Bayes.X = zeros(numel(X),1);
                    obj.Bayes.Posterior{iprobe,phs} = NaN(numel(X),obj.Bayes.numBins);
                    obj.Bayes.PosError{iprobe,phs} = NaN(numel(X),obj.Bayes.numBins);
                    obj.Bayes.DistError{iprobe,phs} = NaN(numel(X),obj.Bayes.numBins);
                    
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
                        
                        obj.Bayes.decOder = TbayesDecoder;
                        obj.Bayes.decOder.numBins = obj.Bayes.numBins;
                        obj.Bayes.decOder.Fcircular = obj.data.es.CircularMaze;
                        obj.Bayes.decOder.Flookuptable = obj.Bayes.Flookuptable;
                        obj.Bayes.decOder.kfold = obj.Bayes.kfold;
                        obj.Bayes.decOder.FoptiSmooth = obj.Bayes.FoptiSmooth;
                        
                        [obj.Bayes.decOder, pred, Xdec, Posterior, ~, nonNormPosterior] = obj.Bayes.decOder.trainDecoder(X(idx & spdidx), Y(idx & spdidx ,:), Yfilt(idx & spdidx ,:), T(idx & spdidx ), obj.Bayes.Xsmth_win);
                        
                        
                        obj.Bayes.prediction{iprobe,phs}(idx & spdidx ,:) = pred;
                        obj.Bayes.X(idx & spdidx ,:) = Xdec(1:min(end,numel(pred)));
                        obj.Bayes.Posterior{iprobe,phs}(idx & spdidx ,:) = Posterior;%nonNormPosterior;%Posterior;
                        
                        %                     [pred, Xdec, Posterior, ~, nonNormPosterior] = obj.Bayes.decOder.predictBayesDecoder(X(~idx | ~spdidx), Yfilt(~idx | ~spdidx,:), T(~idx | ~spdidx), obj.Bayes.Tsmth_win, obj.Bayes.type);
                        %                     if irun == 1
                        %                         testidx = find(~idx | ~spdidx);
                        %                         postatXdec = zeros(size(Posterior,1),1);
                        %                         postatXdecNew = zeros(size(Posterior,1),1);
                        %                         for i = 1:max(Xdec)
                        %                             ttidx = find(Xdec == i);
                        %                             postatXdecNew(ttidx) = Posterior(ttidx,i);
                        %                             postatXdec(ttidx) = obj.Bayes.Posterior{iprobe,phs}(testidx(ttidx),i);
                        %                         end
                        %                         goodspdidx = postatXdecNew>postatXdec;
                        %                         obj.Bayes.bestquantile(testidx(goodspdidx)) = spd;%./obj.Bayes.Spdbin(testidx(goodspdidx));
                        %                         obj.Bayes.prediction{iprobe,phs}(testidx(goodspdidx),:) = pred(goodspdidx);
                        %                         obj.Bayes.X(~idx | ~spdidx,:) = Xdec(1:min(end,numel(pred)));
                        %                         obj.Bayes.Posterior{iprobe,phs}(testidx(goodspdidx),:) = Posterior(goodspdidx,:);%nonNormPosterior;%Posterior;%
                        %                     else
                        
                        [pred, Xdec, Posterior, ~, nonNormPosterior] = obj.Bayes.decOder.predictBayesDecoder(X(~idx & spdidx ,:), Yfilt(~idx & spdidx ,:), T(~idx & spdidx ), obj.Bayes.type);
                        
                        obj.Bayes.prediction{iprobe,phs}(~idx & spdidx ,:) = pred;
                        obj.Bayes.X(~idx & spdidx ,:) = Xdec(1:min(end,numel(pred)));
                        obj.Bayes.Posterior{iprobe,phs}(~idx & spdidx ,:) = Posterior;%nonNormPosterior;%Posterior;%
                        %                     end
                    end
                end
                
                Prange = size(obj.Bayes.Posterior{iprobe,phs},2);
                Xrange = max(obj.Bayes.X);
                for i = 1:Xrange
                    obj.Bayes.PosError{iprobe,phs}(obj.Bayes.X == i,:) = circshift(obj.Bayes.Posterior{iprobe,phs}(obj.Bayes.X == i,:),floor(Prange/2)-i+1,2);
                end

                Prange = size(obj.Bayes.Posterior{iprobe,phs},2);
                dist = floor(obj.data.es.traj ./ (obj.data.es.gain/mode(obj.data.es.gain)));
                [obj.Bayes.D,~] = normalise1var(dist, round(obj.Bayes.numBins/min(obj.data.es.gain/mode(obj.data.es.gain))));
                Xrange = max(obj.Bayes.D);
                for i = 1:Xrange
                    obj.Bayes.DistError{iprobe,phs}(obj.Bayes.D == i,:) = circshift(obj.Bayes.Posterior{iprobe,phs}(obj.Bayes.D == i,:),floor(Prange/2)-i,2);
                end
                obj.Bayes.Posterior0{iprobe}(phsidx,:) = obj.Bayes.Posterior0{iprobe}(phsidx,:) + obj.Bayes.Posterior{iprobe,phs}(phsidx,:);
                obj.Bayes.PosError0{iprobe}(phsidx,:) = obj.Bayes.PosError0{iprobe}(phsidx,:) + obj.Bayes.PosError{iprobe,phs}(phsidx,:);
                obj.Bayes.DistError0{iprobe}(phsidx,:) = obj.Bayes.DistError0{iprobe}(phsidx,:) + obj.Bayes.DistError{iprobe,phs}(phsidx,:);                
            end
            
            %     obj.Bayes.Posterior00 = obj.Bayes.Posterior0;
            %
            obj.Bayes.X(obj.Bayes.X == 0) = max(obj.Bayes.X);
            
            %     phsval = obj.Bayes.LFPphase{iprobe};
            %     idxphschange = [1 ; find(diff(phsval)< -180)+1 ; numel(phsval)+1];
            %     for i = 1:numel(idxphschange)-1
            %         idxmerge = idxphschange(i):idxphschange(i+1)-1;
            %         bestidx = find(Err(idxmerge)-mean(Err(idxmerge)) == min(Err(idxmerge)-mean(Err(idxmerge))),1,'first');
            %         obj.Bayes.Posterior0(idxmerge,:) = repmat(obj.Bayes.Posterior0(idxmerge(bestidx),:),numel(idxmerge),1);
            %     end
            
            %     obj.Bayes.nthetaphsbins = 12;
            %     phsval = mod(obj.Bayes.LFPphase{iprobe},360);
            %     idxphschange = [1 ; find(diff(phsval)< -180)+1 ; numel(phsval)+1];
            %     Prange = size(obj.Bayes.PosError0,2);
            %     Ntimebins = size(obj.Bayes.PosError0,1);
            %     phsbins = floor(obj.Bayes.LFPphase{iprobe}/(360/obj.Bayes.nthetaphsbins));
            %     phsbins = mod(phsbins,obj.Bayes.nthetaphsbins)+1;
            %     thetaErrorVec = zeros(numel(idxphschange)-1,Prange*obj.Bayes.nthetaphsbins);
            %     thetatrainidx = false(numel(idxphschange)-1,1);
            %     thetaphs = [];
            %     fulltimetheta = zeros(numel(idxphschange)-1,obj.Bayes.nthetaphsbins);
            %     for i = 1:numel(idxphschange)-1
            %         idxmerge = max(1,idxphschange(i)-1):min(Ntimebins,idxphschange(i+1));
            %         thetatrainidx(i) = sum(idx(idxmerge)) > 0;
            %         thetaphsbins = phsbins(idxmerge);
            %         thetaphsbins(1) = thetaphsbins(1) - obj.Bayes.nthetaphsbins;
            %         thetaphsbins(end) = thetaphsbins(end) + obj.Bayes.nthetaphsbins;
            %         thetaphsbinsval = unique(thetaphsbins);
            %         thetaError = zeros(numel(thetaphsbinsval),Prange);
            %         timetheta = zeros(numel(thetaphsbinsval),1);
            %         for ph = 1:numel(thetaphsbinsval)
            %             thetaError(ph,:) = mean(obj.Bayes.PosError0(idxmerge(thetaphsbins==thetaphsbinsval(ph)),:),1);
            %             timetheta(ph) = idxmerge(find(thetaphsbins==thetaphsbinsval(ph),1,'first'));
            %         end
            %         thetaErrorFull = zeros(obj.Bayes.nthetaphsbins,Prange);
            %         for xx = 1:Prange
            %             thetaErrorFull(:,xx) = interp1(thetaphsbinsval,thetaError(:,xx),1:obj.Bayes.nthetaphsbins);
            %         end
            %
            %         timetheta = sort(timetheta,'ascend');
            %         fulltimetheta(i,:) =  interp1(thetaphsbinsval,timetheta,1:obj.Bayes.nthetaphsbins);
            %
            %         thetaphs = [thetaphs 1:obj.Bayes.nthetaphsbins];
            %         thetaErrorVec(i,:) = reshape(thetaErrorFull,1,Prange*obj.Bayes.nthetaphsbins);
            %     end
            %     thetaphs = thetaphs(:);
            %     [coeff,score,latent] = pca(thetaErrorVec,'Centered',false);
            %     thetaErrorFull = reshape(score(:,1:20)*coeff(:,1:20)',size(score,1)*size(score,2)/Prange,Prange);
            %     fulltimetheta = fulltimetheta(:);
            %     posError0 = zeros(size(obj.Bayes.PosError0));
            %     for xx = 1:Prange
            %         posError0(:,xx) = interp1(fulltimetheta,thetaErrorFull(:,xx),1:size(posError0,1));
            %     end
            %     obj.Bayes.PosError0 = posError0;
            %     for i = 1:Xrange
            %         obj.Bayes.Posterior0(obj.Bayes.X == i,:) = circshift(posError0(obj.Bayes.X == i,:),-(floor(Prange/2)-i),2);
            %     end
            %     obj.Bayes.nthetaphsbins = 1;
            
            %     maxtol = 1;
            %     [goodidx, ~, Err] = getcleandecidx(obj.Bayes.Posterior0, obj.Bayes.X,obj.Bayes.goodidx_th,maxtol,obj.data.es.trialID);
            %     if irun == 1
            %         theta = LFPsignals(Err', 60, 2, 16);
            %         obj.Bayes.LFPphase{iprobe} = mod(round(theta.phase'/(2*pi)*360),360);
            %     end
            
            obj.Bayes.time = ((1:numel(obj.Bayes.X))./obj.Bayes.sampleRate)';
            
            obj.Bayes.X(obj.Bayes.X == 0) = max(obj.Bayes.X);
            maxtol = 0.1;
            [goodidx, ~] = getcleandecidx(obj.Bayes.Posterior0{iprobe}, obj.Bayes.X,obj.Bayes.goodidx_th,maxtol,obj.data.es.trialID);
            gainratio = obj.data.es.gain/mode(obj.data.es.gain);
            nanidx = find(isnan(gainratio));
            for i = 1:numel(nanidx)
                gainratio(nanidx(i)) = gainratio(nanidx(i)-1);
            end
            gainratioval = unique(gainratio(~isnan(gainratio)));
            obj.Bayes.PosteriorXtrain = zeros(size(obj.Bayes.Posterior0{iprobe},2));
            obj.Bayes.PosteriorXdist = cell(1,numel(gainratioval));
            for g = 1:numel(gainratioval)
                obj.Bayes.PosteriorXdist{g} = zeros(size(obj.Bayes.Posterior0{iprobe},2));
            end
            for i = 1:size(obj.Bayes.Posterior0{iprobe},2)
                obj.Bayes.PosteriorXtrain(:,i) = mean(obj.Bayes.Posterior0{iprobe}(obj.Bayes.X == i  & idxref & goodidx,:),1)';
            end
            for i = 1:size(obj.Bayes.Posterior0{iprobe},2)
                for g = 1:numel(gainratioval)
                    obj.Bayes.PosteriorXdist{g}(:,i) = circshift(obj.Bayes.PosteriorXtrain(i,:),round(i/gainratioval(g)-i),2);
                end
            end
            
            obj.Bayes.Kroom = zeros(size(obj.Bayes.Posterior0{iprobe},1),1);
            obj.Bayes.Kdist = zeros(size(obj.Bayes.Posterior0{iprobe},1),1);
            for i = 1:size(obj.Bayes.Posterior0{iprobe},1)
                matrefX = obj.Bayes.PosteriorXtrain(:,obj.Bayes.X(i));
                matrefX = (matrefX-min(matrefX))/(max(matrefX)-min(matrefX));
                matrefD = obj.Bayes.PosteriorXtrain(:,mod(round(obj.Bayes.X(i)/gainratio(i)),obj.Bayes.numBins)+1);%obj.Bayes.PosteriorXdist{gainratioval == gainratio(i)}(:,obj.Bayes.X(i));
                matrefD = (matrefD-min(matrefD))/(max(matrefD)-min(matrefD));
                post0 = obj.Bayes.Posterior0{iprobe}(i,:);
                post0 = (post0-min(post0))/(max(post0)-min(post0));
                post0Dist = circshift(obj.Bayes.Posterior0{iprobe}(i,:),-round(obj.Bayes.X(i)/gainratio(i)-obj.Bayes.X(i)),2);%obj.Bayes.Posterior0{iprobe}(i,:)
                obj.Bayes.Kroom(i) = (post0*matrefX)/sqrt(sum(post0.^2)*sum(matrefX.^2));
                obj.Bayes.Kdist(i) = (post0*matrefD)/sqrt(sum(post0.^2)*sum(matrefD.^2));%((post0Dist)*(matrefX))/sqrt(sum((post0Dist).^2)*sum((matrefX).^2));%(post0tt*matrefD)/sqrt(sum(post0tt.^2)*sum(matrefD.^2));
                %         obj.Bayes.Kroom(i) = ((obj.Bayes.Posterior0{iprobe}(i,:)-mean(obj.Bayes.Posterior0{iprobe}(i,:)))*(matrefX-mean(matrefX)))/sqrt(sum((obj.Bayes.Posterior0{iprobe}(i,:)-mean(obj.Bayes.Posterior0{iprobe}(i,:))).^2)*sum((matrefX-mean(matrefX)).^2));
                %         obj.Bayes.Kdist(i) = ((post0Dist-mean(post0Dist))*(matrefX-mean(matrefX)))/sqrt(sum((post0Dist-mean(post0Dist)).^2)*sum((matrefX-mean(matrefX)).^2));%(post0tt*matrefD)/sqrt(sum(post0tt.^2)*sum(matrefD.^2));
                
                %         corrcoeff = obj.Bayes.Posterior0{iprobe}(i,:)*obj.Bayes.PosteriorXtrain;
                %         corrcoeff = corrcoeff./sqrt((sum(obj.Bayes.Posterior0{iprobe}(i,:).^2).*sum(obj.Bayes.PosteriorXtrain.^2,1)));
                %         [~,bestX] = max(corrcoeff);
                %         obj.Bayes.Posterior0{iprobe}(i,:) = obj.Bayes.PosteriorXtrain(:,bestX);
            end
            
            %     KroomX = zeros(1, Xrange);
            %     KroomXref = zeros(1, Xrange);
            %     KdistX = zeros(1, Xrange);
            %     KdistXref = zeros(1, Xrange);
            %     for i = 1:Xrange
            %         KroomXref(i) = nanmean(Kroom(X == i  & tidxref & phsidx & goodidx));
            %         KdistXref(i) = nanmean(Kdist(X == i  & tidxref & phsidx & goodidx));
            %         KdistX(i) = nanmean(Kdist(X == i  & tidx & phsidx & goodidx));
            %         KroomX(i) = nanmean(Kroom(X == i  & tidx & phsidx & goodidx));
            %     end
            
            %     if obj.Bayes.nthetaphsbins > 1
            %         baseline = 1./EXP.Bayes.decOder.numBins;
            %         thetasmthwin = 9;
            %         obj.Bayes.Posterior0{iprobe} = conv2(obj.Bayes.Posterior0{iprobe},1/thetasmthwin*ones(thetasmthwin,1),'same');
            %         obj.Bayes.Posterior0{iprobe} = obj.Bayes.Posterior0{iprobe}./repmat(sum(obj.Bayes.Posterior0{iprobe},2),1,100);
            %         obj.Bayes.Posterior0{iprobe} = obj.Bayes.Posterior0{iprobe}./baseline;
            %
            %         obj.Bayes.nPosterior0{iprobe} = conv2(obj.Bayes.nPosterior0{iprobe},1/thetasmthwin*ones(thetasmthwin,1),'same');
            %         obj.Bayes.nPosterior0{iprobe} = obj.Bayes.nPosterior0{iprobe}./repmat(sum(obj.Bayes.nPosterior0{iprobe},2),1,100);
            %         obj.Bayes.nPosterior0{iprobe} = obj.Bayes.nPosterior0{iprobe}./baseline;
            %
            %         obj.Bayes.PosError0 = conv2(obj.Bayes.PosError0,1/thetasmthwin*ones(thetasmthwin,1),'same');
            %         obj.Bayes.PosError0 = obj.Bayes.PosError0./repmat(sum(obj.Bayes.PosError0,2),1,100);
            %         obj.Bayes.PosError0 = obj.Bayes.PosError0./baseline;
            %
            %         obj.Bayes.DistError0 = conv2(obj.Bayes.DistError0,1/thetasmthwin*ones(thetasmthwin,1),'same');
            %         obj.Bayes.DistError0 = obj.Bayes.DistError0./repmat(sum(obj.Bayes.DistError0,2),1,100);
            %         obj.Bayes.DistError0 = obj.Bayes.DistError0./baseline;
            %     end
            
            %     goodidx = getcleandecidx(obj.Bayes.Posterior0{iprobe}, obj.Bayes.X,obj.Bayes.goodidx_th);
            %     Xrange = max(obj.Bayes.X);
            %     Prange = size(obj.Bayes.Posterior0{iprobe},2);
            %     mat = zeros(Prange, Xrange);
            %     for i = 1:Xrange
            %         mat(:,i) = mean(obj.Bayes.Posterior0{iprobe}(obj.Bayes.X == i  & idx & goodidx,:),1)';
            %     end
            %     mattemp = conv2(repmat(mat,3,3), 1/9*ones(3), 'same');
            %     mat = mattemp(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
            %
            %     postfilt = zeros(Prange,1);
            %     postfilt(1:Prange) = (0:(Prange-1))'/(Prange)*(2*pi)-pi/2;
            %
            %     maxtol = 0.001;
            %     [maxval, ~] = max(mat,[],1);
            %     for i = 1:Xrange
            %         mat(mat(:,i) < (1-maxtol)*maxval(i),i) = 0;
            %     end
            %
            %     a = sum((mat'*cos(postfilt)),2);b = sum((mat'*sin(postfilt)),2);
            %     predave = atan2(-a,b)+pi;%./sum(exp(log(2)*Post(tidx,:)),2);
            %     predave = predave/(2*pi)*Prange+1;
            %     predplus = predave + Prange;
            %     predminus = predave - Prange;
            %     Xtemp = (1:Xrange)';
            %     predave = predave.* (abs(Xtemp-predave)<abs(Xtemp-predplus) & abs(Xtemp-predave)<abs(Xtemp-predminus)) + predplus.* (abs(Xtemp-predplus)<abs(Xtemp-predave) & abs(Xtemp-predplus)<abs(Xtemp-predminus)) + predminus.* (abs(Xtemp-predminus)<abs(Xtemp-predave) & abs(Xtemp-predminus)<abs(Xtemp-predplus));
            %     predave = smooth(predave,3);
            %     diffX = diff(X);
            %     resolX = 0.001;
            %     Xvalues = 0:resolX:Xrange;
            %     predave = [predave;predave(1)+Xrange];
            %     predave = interp1(0:Xrange,predave,Xvalues);
            %     X = predave(floor(X/resolX)+1)';
        end
        obj.Bayes.decOder = [];
        toc
    else
        for spd = 1:obj.Bayes.nspdbins(iprobe)
            for phs = 1:obj.Bayes.nthetaphsbins
                obj.Bayes.Posterior{iprobe,phs} = NaN(size(obj.Bayes.Posterior{1,phs}));
                obj.Bayes.PosError{iprobe,phs} = NaN(size(obj.Bayes.PosError{1,phs}));
                obj.Bayes.DistError{iprobe,phs} = NaN(size(obj.Bayes.DistError{1,phs}));
                obj.Bayes.prediction{iprobe,phs} = NaN(size(obj.Bayes.prediction{1,phs}));
            end
        end
        obj.Bayes.LFPphase{iprobe} = NaN(size(obj.Bayes.LFPphase{1}));
        obj.Bayes.Posterior0{iprobe} = NaN(size(obj.Bayes.Posterior0{1}));
        obj.Bayes.PosError0{iprobe} = NaN(size(obj.Bayes.PosError0{1}));
        obj.Bayes.DistError0{iprobe} = NaN(size(obj.Bayes.DistError0{1}));
    end
end
end