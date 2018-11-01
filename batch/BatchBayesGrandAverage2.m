function res = BatchBayesGrandAverage2(batch2p,params)
SetDirs;
if nargin < 1
    strlistvarname = {'2p data','electrophys data'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
    if ok && varnamesel == 1
        batch2p = true;
    elseif ok
        batch2p = false;
    end
end
if batch2p
    expt = getExperimentList2p;
    datadir = 'D:\DATA\batch\2p';%DIRS.data2p;
    res.nSpdbins = 5;
    res.nVisbins = 5;
    res.nEyebins = 1;
    res.nPhsbins = 18;
else
    expt = getExperimentList;
    datadir = 'D:\DATA\batch';%DIRS.multichanspikes;
    res.nSpdbins = 5;%res.nspeedbins;
    res.nVisbins = 5;%
    res.nEyebins = 5;%res.neyebins;
    res.nPhsbins = 18;%res.nspeedbins;
end

% %default values
% params.SpeedThreshold = 5;
% params.nphsbins = 1;
% params.nspeedbins = 5;
% params.neyebins = 1;
% params.Tsmthwin = 15;
% params.Xsmthwin = 4;
% params.Tsmthwin_dec = 250;
% params.cellstr = 'goodonly';
% params.DecAveType = 'DecPost';%'DecPost' or 'DecError'

if nargin < 2
    prompt = {'Speed Threshold';'nphsbins';'nspdbins';'neyebins';'Window size (ms)';'Spatial smth (%)';'Decoding window size (ms)';'suffix';'error type'};
    dlg_title = 'Parameters';
    num_lines = 1;
    def = {'5';'1';'3';'1';'15';'4';'50';'goodonly';'DecPost'};
    nruns = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(nruns)
        res.SpeedThreshold = str2num(nruns{1});
        res.nphsbins = str2num(nruns{2});
        res.nspeedbins = str2num(nruns{3});
        res.neyebins = str2num(nruns{4});
        res.Tsmthwin = str2num(nruns{5});
        res.Xsmthwin = str2num(nruns{6});
        res.Tsmthwin_dec = str2num(nruns{7});
        res.cellstr = nruns{8};
        res.DecAveType = nruns{9};%'DecPost' or 'DecError'
    end
else
    res.SpeedThreshold = params.SpeedThreshold;
    res.nphsbins = params.nphsbins;
    res.nspeedbins = params.nspeedbins;
    res.neyebins = params.neyebins;
    res.Tsmthwin = params.Tsmthwin;
    res.Xsmthwin = params.Xsmthwin;
    res.Tsmthwin_dec = params.Tsmthwin_dec;
    res.cellstr = params.cellstr;
    res.DecAveType = params.DecAveType;%'DecPost' or 'DecError'
end
savedfilename_res = ['D:\DATA\batch\All\res_' res.DecAveType '_Twin' num2str(res.Tsmthwin) '_Xwin' num2str(res.Xsmthwin) '_spdth' num2str(res.SpeedThreshold)...
                     '_Decwin' num2str(res.Tsmthwin_dec) '_' num2str(res.nspeedbins) 'speedbins_' num2str(res.neyebins) 'eyebins_' num2str(res.nphsbins) 'thetabins_' res.cellstr '_spdquantiles_all' '.mat'];
savedfilename_popres_sessions = ['D:\DATA\batch\All\popres_sessions' res.DecAveType '_Twin' num2str(res.Tsmthwin) '_Xwin' num2str(res.Xsmthwin) '_spdth' num2str(res.SpeedThreshold)...
                     '_Decwin' num2str(res.Tsmthwin_dec) '_' num2str(res.nspeedbins) 'speedbins_' num2str(res.neyebins) 'eyebins_' num2str(res.nphsbins) 'thetabins_' res.cellstr '_spdquantiles_all' '.mat'];
savedfilename_popres_trials = ['D:\DATA\batch\All\popres_trials' res.DecAveType '_Twin' num2str(res.Tsmthwin) '_Xwin' num2str(res.Xsmthwin) '_spdth' num2str(res.SpeedThreshold)...
                     '_Decwin' num2str(res.Tsmthwin_dec) '_' num2str(res.nspeedbins) 'speedbins_' num2str(res.neyebins) 'eyebins_' num2str(res.nphsbins) 'thetabins_' res.cellstr '_spdquantiles_all' '.mat'];
                
% res.Tsmthwin = 50;%250;%250;%250;%150;%300;%40;%120;%50
% res.Xsmthwin = 2;%1;%
% res.SpeedThreshold = 5;
% res.nspeedbins = 3;%5;
% res.neyebins = 3;%5;
% res.nphsbins = 0;%1;%
% res.cellstr = 'goodonly';%'All_50bins';%'goodonly';%'goodonly_unwrapped';%'goodonly';%'All';%
if batch2p
%     res.Tsmthwin = 250;
%     res.SpeedThreshold = 1;
%     res.nspeedbins = 3;
%     res.neyebins = 1;%5;
%     res.nphsbins = 0;%1;%
%     res.cellstr = 'goodonly_spktimes';
    filesuffix_EXP = ['Twin' num2str(res.Tsmthwin) '_' 'Xwin' num2str(res.Xsmthwin) '_' 'spdth' num2str(res.SpeedThreshold) '_' num2str(res.nspeedbins) 'speedbins' '_' num2str(res.neyebins) 'eyebins' '_' num2str(res.nphsbins) 'thetabins' '_' res.cellstr];
else
    filesuffix_EXP = ['Twin' num2str(res.Tsmthwin) '_' 'Xwin' num2str(res.Xsmthwin) '_' 'spdth' num2str(res.SpeedThreshold) '_' 'Decwin' num2str(res.Tsmthwin_dec) '_' num2str(res.nspeedbins) 'speedbins' '_' num2str(res.neyebins) 'eyebins' '_' num2str(res.nphsbins) 'thetabins' '_' res.cellstr '_spdquantiles'];
end
% filesuffix_EXP = ['Twin' num2str(res.Tsmthwin) '_' 'Xwin' num2str(res.Xsmthwin) '_' 'spdth' num2str(res.SpeedThreshold) '_' num2str(res.nspeedbins) 'speedbins' '_' num2str(res.nphsbins) 'thetabins' '_' res.cellstr];
disp(filesuffix_EXP);
disp(res.DecAveType);

res.maxtol_ave = 1;
res.maxtol_max = 0.1;
res.amp_th_ave = 0;%1
res.amp_th_max = 0;%1
res.decPostmin = 1.5;
kfold = 20;
contval = [0.2:0.05:0.9];
outcomeVal = 2;
RLval = 1;

lambdaSmooth = 2;
corrmaxlag = 180;
nanimal = numel(expt);

for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        if batch2p
            dDIRname = [datadir  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];
        else
            dDIRname = [datadir  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];
        end
        if exist([dDIRname filesep 'EXP_' filesuffix_EXP '.mat'],'file')
            S = load([dDIRname filesep 'EXP_' filesuffix_EXP '.mat']);
            EXP = TVRData;
            EXP.Copyobj(S.EXP);
%             Slat = load([datadir filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'decoderParams_latency.mat']);
            
            res.gainVal{ianimal,iseries} = EXP.SubsetVal.gain;
            res.roomlengthVal{ianimal,iseries} = EXP.SubsetVal.roomlength;
            res.outcomeVal{ianimal,iseries} = EXP.SubsetVal.outcome;
            res.gain{ianimal,iseries} = EXP.data.es.gain;
            res.contrast{ianimal,iseries} = EXP.data.es.contrast;
            res.outcome{ianimal,iseries} = EXP.data.es.outcome;
            res.trialID{ianimal,iseries} = EXP.data.es.trialID;
            res.blanks{ianimal,iseries} = EXP.data.es.blanks;
            res.afterblanks{ianimal,iseries} = EXP.data.es.afterblanks;
            res.trialgainchange{ianimal,iseries} = EXP.data.es.trialgainchange;
            
            res.traj{ianimal,iseries} = EXP.data.es.traj;
            res.trajspeed{ianimal,iseries} = EXP.data.es.trajspeed;
            res.ballspeed{ianimal,iseries} = EXP.data.es.ballspeed;
            
            res.lick{ianimal,iseries} = EXP.data.es.lick;
            res.firstgoodlick{ianimal,iseries} = EXP.data.es.firstgoodlick;
            res.firstrewlick{ianimal,iseries} = EXP.data.es.firstrewlick;
            res.goodlick{ianimal,iseries} = EXP.data.es.goodlick;
            res.preRewlick{ianimal,iseries} = EXP.data.es.preRewlick;
            res.postRewlick{ianimal,iseries} = EXP.data.es.postRewlick;
            res.passivelick{ianimal,iseries} = EXP.data.es.postRewlick;
            res.firstbadlick{ianimal,iseries} = EXP.data.es.postRewlick;
            res.badlick{ianimal,iseries} = EXP.data.es.postRewlick;
            
            res.pupilSize{ianimal,iseries} = EXP.data.es.pupilSize;
            res.eyeYpos{ianimal,iseries} = EXP.data.es.eyeYpos;
            res.eyeXpos{ianimal,iseries} = EXP.data.es.eyeXpos;
            
            eyepos = NaN(size(EXP.data.es.eyeXpos));
            eyepos(~isnan(EXP.data.es.eyeXpos)) = smthInTime(EXP.data.es.eyeXpos(~isnan(EXP.data.es.eyeXpos)), EXP.Bayes.sampleRate, EXP.Bayes.Tsmth_win(1), 'same', [], 'boxcar_centered');
            res.eyepos{ianimal,iseries} = eyepos;
            if isfield(EXP.Bayes,'Eyebin')
                res.Eyebin{ianimal,iseries} = EXP.Bayes.Eyebin;
            else
                res.Eyebin{ianimal,iseries} = NaN(size(eyepos));
            end
            
            runspeed = NaN(size(EXP.data.es.ballspeed));
            runspeed(~isnan(EXP.data.es.ballspeed)) = smthInTime(EXP.data.es.ballspeed(~isnan(EXP.data.es.ballspeed)), EXP.Bayes.sampleRate, EXP.Bayes.Tsmth_win(1), 'same', [], 'boxcar_centered');
            res.runSpeed{ianimal,iseries} = runspeed;
            visspeed = NaN(size(EXP.data.es.trajspeed));
            visspeed(~isnan(EXP.data.es.trajspeed)) = smthInTime(EXP.data.es.trajspeed(~isnan(EXP.data.es.trajspeed)), EXP.Bayes.sampleRate, EXP.Bayes.Tsmth_win(1), 'same', [], 'boxcar_centered');            
            res.visSpeed{ianimal,iseries} = visspeed;
            res.Spdbin{ianimal,iseries} = EXP.Bayes.Spdbin;
            
            res.X{ianimal,iseries} = EXP.Bayes.Xsmth0{1};%EXP.Bayes.X0{1};
            res.firstgoodlick{ianimal,iseries} = EXP.data.es.firstgoodlick;
            res.contrastVal{ianimal,iseries} = EXP.SubsetVal.contrast;            
            
            for iprobe = 1:numel(EXP.Bayes.LFPphase)
                res.BayesLFPphase{ianimal,iseries,iprobe} = EXP.Bayes.LFPphase{iprobe};
                res.Phsbin{ianimal,iseries,iprobe} = NaN(size(res.BayesLFPphase{ianimal,iseries,1}));
                if res.nPhsbins > 1
                    for iphsbin = 1:res.nPhsbins
                        res.Phsbin{ianimal,iseries,iprobe}(mod(res.BayesLFPphase{ianimal,iseries,iprobe},360) >= (iphsbin-1)*360/res.nPhsbins & mod(res.BayesLFPphase{ianimal,iseries,iprobe},360) < iphsbin*360/res.nPhsbins) = iphsbin;
                    end
                else
                    res.Phsbin{ianimal,iseries,iprobe} = ones(size(res.Phsbin{ianimal,iseries}));
                end
            end
            
            itraj = res.X{ianimal,iseries};
            ntrajbins = max(itraj);
                       
            contref = find(ismember(res.contrastVal{ianimal,iseries}, contval));
            RLref = find(ismember(res.roomlengthVal{ianimal,iseries}, [1]));
            outcomeref = find(ismember(res.outcomeVal{ianimal,iseries}, [2]));
            gainref = find(ismember(res.gainVal{ianimal,iseries}, mode(res.gain{ianimal,iseries}(~isnan(res.gain{ianimal,iseries})))));%gainref = 2;
            idxref = EXP.getSubsets(contref, gainref, RLref, outcomeref, EXP.Bayes.speed_th);%res.tidx{ianimal,iseries,1,contref,gainref,RLref,outcomeref};
            
            if res.nSpdbins >= 3
                res.Spdbin2{ianimal,iseries} = NaN(size(res.runSpeed{ianimal,iseries}));
                speed = res.runSpeed{ianimal,iseries};
                speedprofile = NaN(EXP.Bayes.numBins,1);
                for xx = 1:EXP.Bayes.numBins
                    speedprofile(xx) = median(speed(idxref & res.X{ianimal,iseries}==xx));
                end
                %speedprofile = fast1Dmap(res.X{ianimal,iseries}(idxref),res.runSpeed{ianimal,iseries}(idxref),1,1,1/(EXP.Bayes.Xsmth_win/EXP.Bayes.numBins),EXP.data.es.CircularMaze);
                speed = (speed - speedprofile(res.X{ianimal,iseries}));%./speedprofile(res.X{ianimal,iseries});
                speedrange = 10;%from -5 cm/s to +5 cm/s around the median speed
                if res.nSpdbins > 3
                    binsize = speedrange/(res.nSpdbins-2 - 1);
                else
                    binsize = 0;
                end
                minbinval = -floor(speedrange/2) - binsize/2;
                maxbinval = floor(speedrange/2) + binsize/2;
                speed(speed < minbinval) = minbinval-1;
                speed(speed > maxbinval) = maxbinval+1;
                [res.Spdbin2{ianimal,iseries}(speed >= minbinval & speed <= maxbinval), ~] = normalise1var(speed(speed >= minbinval & speed <= maxbinval), res.nSpdbins-2,[],[minbinval maxbinval]);
                res.Spdbin2{ianimal,iseries}(speed >= minbinval & speed <= maxbinval) = res.Spdbin2{ianimal,iseries}(speed >= minbinval & speed <= maxbinval) + 1;
                res.Spdbin2{ianimal,iseries}(speed < minbinval) = 1;
                res.Spdbin2{ianimal,iseries}(speed > maxbinval) = res.nSpdbins;
            else
                res.Spdbin2{ianimal,iseries} = ones(size(res.runSpeed{ianimal,iseries}));
                res.nSpdbins = 1;
            end
            
%             spdquantilelim = zeros(ntrajbins,2);
%             res.Spdbin2{ianimal,iseries} = NaN(size(res.runSpeed{ianimal,iseries}));
%             if res.nSpdbins > 1
%                 for spd = 1:res.nSpdbins
%                     for xx = 1:ntrajbins
%                         spdquantilelim(xx,1) = quantile(res.runSpeed{ianimal,iseries}(idxref & itraj == xx),max(0,(spd-1)/res.nSpdbins));
%                         spdquantilelim(xx,2) = quantile(res.runSpeed{ianimal,iseries}(idxref & itraj == xx),min(1,(spd)/res.nSpdbins));
%                     end
%                     res.Spdbin2{ianimal,iseries}(res.runSpeed{ianimal,iseries} >= spdquantilelim(itraj,1) & res.runSpeed{ianimal,iseries} < spdquantilelim(itraj,2)) = spd;
%                     if spd == 1
%                         res.Spdbin2{ianimal,iseries}(res.runSpeed{ianimal,iseries} <= spdquantilelim(itraj,1)) = spd;
%                     end
%                     if spd == res.nSpdbins
%                         res.Spdbin2{ianimal,iseries}(res.runSpeed{ianimal,iseries} >= spdquantilelim(itraj,2)) = spd;
%                     end
%                 end
%             else
%                 res.Spdbin2{ianimal,iseries} = ones(size(res.Spdbin2{ianimal,iseries}));
%             end

            if res.nVisbins >= 3
                res.Visbin2{ianimal,iseries} = NaN(size(res.visSpeed{ianimal,iseries}));
                visspeed = res.visSpeed{ianimal,iseries};
                visspeedprofile = NaN(EXP.Bayes.numBins,1);
                for xx = 1:EXP.Bayes.numBins
                    visspeedprofile(xx) = median(visspeed(idxref & res.X{ianimal,iseries}==xx));
                end
                %visspeedprofile = fast1Dmap(res.X{ianimal,iseries}(idxref),res.visSpeed{ianimal,iseries}(idxref),1,1,1/(EXP.Bayes.Xsmth_win/EXP.Bayes.numBins),EXP.data.es.CircularMaze);
                visspeed = (visspeed - visspeedprofile(res.X{ianimal,iseries}));%./visspeedprofile(res.X{ianimal,iseries});
                visspeedrange = 10;%
                if res.nVisbins > 3
                    binsize = visspeedrange/(res.nVisbins-2 - 1);
                else
                    binsize = 0;
                end
                minbinval = -floor(visspeedrange/2) - binsize/2;
                maxbinval = floor(visspeedrange/2) + binsize/2;
                visspeed(visspeed < minbinval) = minbinval-1;
                visspeed(visspeed > maxbinval) = maxbinval+1;
                [res.Visbin2{ianimal,iseries}(visspeed >= minbinval & visspeed <= maxbinval), ~] = normalise1var(visspeed(visspeed >= minbinval & visspeed <= maxbinval), res.nVisbins-2,[],[minbinval maxbinval]);
                res.Visbin2{ianimal,iseries}(visspeed >= minbinval & visspeed <= maxbinval) = res.Visbin2{ianimal,iseries}(visspeed >= minbinval & visspeed <= maxbinval) + 1;
                res.Visbin2{ianimal,iseries}(visspeed < minbinval) = 1;
                res.Visbin2{ianimal,iseries}(visspeed > maxbinval) = res.nVisbins;
            else
                res.Visbin2{ianimal,iseries} = ones(size(res.visSpeed{ianimal,iseries}));
                res.nVisbins = 1;
            end
            
%             visquantilelim = zeros(ntrajbins,2);
%             res.Visbin2{ianimal,iseries} = NaN(size(res.trajspeed{ianimal,iseries}));
%             if res.nVisbins > 1
%                 for ivis = 1:res.nVisbins
%                     for xx = 1:ntrajbins
%                         visquantilelim(xx,1) = quantile(res.visSpeed{ianimal,iseries}(idxref & itraj == xx),max(0,(ivis-1)/res.nVisbins));
%                         visquantilelim(xx,2) = quantile(res.visSpeed{ianimal,iseries}(idxref & itraj == xx),min(1,(ivis)/res.nVisbins));
%                     end
%                     res.Visbin2{ianimal,iseries}(res.visSpeed{ianimal,iseries} >= visquantilelim(itraj,1) & res.visSpeed{ianimal,iseries} < visquantilelim(itraj,2)) = ivis;
%                     if ivis == 1
%                         res.Visbin2{ianimal,iseries}(res.visSpeed{ianimal,iseries} <= visquantilelim(itraj,1)) = ivis;
%                     end
%                     if ivis == res.nVisbins
%                         res.Visbin2{ianimal,iseries}(res.visSpeed{ianimal,iseries} >= visquantilelim(itraj,2)) = ivis;
%                     end
%                 end
%             else
%                 res.Visbin2{ianimal,iseries} = ones(size(res.Visbin2{ianimal,iseries}));
%             end
            
            eyequantilelim = zeros(ntrajbins,2);
            res.Eyebin2{ianimal,iseries} = NaN(size(res.eyeXpos{ianimal,iseries}));
            if res.nEyebins > 1
                for ieye = 1:res.nEyebins
                    for xx = 1:ntrajbins
                        eyequantilelim(xx,1) = quantile(res.eyepos{ianimal,iseries}(idxref & itraj == xx),max(0,(ieye-1)/res.nEyebins));
                        eyequantilelim(xx,2) = quantile(res.eyepos{ianimal,iseries}(idxref & itraj == xx),min(1,(ieye)/res.nEyebins));
                    end
                    res.Eyebin2{ianimal,iseries}(res.eyepos{ianimal,iseries} >= eyequantilelim(itraj,1) & res.eyepos{ianimal,iseries} < eyequantilelim(itraj,2)) = ieye;
                    if ieye == 1
                        res.Eyebin2{ianimal,iseries}(res.eyepos{ianimal,iseries} <= eyequantilelim(itraj,1)) = ieye;
                    end
                    if ieye == res.nEyebins
                        res.Eyebin2{ianimal,iseries}(res.eyepos{ianimal,iseries} >= eyequantilelim(itraj,2)) = ieye;
                    end
                end
            else
                res.Eyebin2{ianimal,iseries} = ones(size(res.Eyebin2{ianimal,iseries}));
            end
            
            
            cont_list = find(ismember(res.contrastVal{ianimal,iseries}, contval));
            RL_list = find(ismember(res.roomlengthVal{ianimal,iseries}, RLval));
            outcome_list = find(ismember(res.outcomeVal{ianimal,iseries}, outcomeVal));
            
            for iprobe = 1:size(EXP.Bayes.Posterior0,1)
                res.DecCells{ianimal,iseries,iprobe} = EXP.Bayes.DecCellidx{iprobe};
                
                nTimes = size(EXP.Bayes.Posterior0{iprobe,1},1);
                Prange = EXP.Bayes.numBins;%size(EXP.Bayes.Posterior0{iprobe,1},2);
                Xrange = EXP.Bayes.numBins;%max(floor(EXP.Bayes.X0{1}));
                baseline = 1/Prange;
                goodPostidx =  ~isnan(sum(EXP.Bayes.Posterior0{iprobe,1},2));
                
%                 if strcmp(res.DecAveType,'DecError')
%                     EXP.Bayes.MaxDecodedPosition0{iprobe,1} = NaN(size(EXP.Bayes.Posterior0{iprobe,1},1),1);
%                     Post = EXP.Bayes.Posterior0{iprobe,1};
%                     Post(Post<=res.decPostmin) = 0;
%                     for tt = 1:size(EXP.Bayes.Posterior0{iprobe,1},1)
%                         [pks,loc] = findpeaks(repmat(Post(tt,:),[1 3]));
%                         loc = loc - Prange;
%                         loc(pks<=res.decPostmin | loc < 1 | loc > Prange) = [];
%                         if ~isempty(loc)
%                             [~,locminidx] = min(abs(Prange/(2*pi)*circ_dist(2*pi/Prange*loc,2*pi/Prange*EXP.Bayes.Xsmth0{1}(tt))));
%                             %             EXP.Bayes.MaxDecodedPosition0{iprobe,1}(tt) = loc(locminidx);
%                             [outfieldidx,fieldix,amp_th_i] = findfield(EXP.Bayes.Posterior0{iprobe,1}(tt,:),1,loc(locminidx),false);
%                             post_tt = EXP.Bayes.Posterior0{iprobe,1}(tt,:);
%                             post_tt(outfieldidx) = 1;
%                             EXP.Bayes.MaxDecodedPosition0{iprobe,1}(tt) = getCircularAverage(post_tt',0,1);%loc(locminidx);
%                             EXP.Bayes.DecodingError0{iprobe,1}(tt) = Prange/(2*pi)*circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe,1}(tt),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tt));
%                         else
%                             EXP.Bayes.MaxDecodedPosition0{iprobe,1}(tt) = NaN;
%                             EXP.Bayes.DecodingError0{iprobe,1}(tt) = NaN;
%                         end
%                     end
%                     goodPostidx = ~isnan(EXP.Bayes.MaxDecodedPosition0{iprobe,1});
%                 end
                EXP.Bayes.MaxDecodedPosition0{iprobe,1}(~goodPostidx) = NaN;
                EXP.Bayes.MeanDecodedPosition0{iprobe,1}(~goodPostidx) = NaN;
                EXP.Bayes.DecodingError0{iprobe,1}(~goodPostidx) = NaN;
                EXP.Bayes.DecodingBias0{iprobe,1}(~goodPostidx) = NaN;
                
%                 EXP.Bayes.PosBias0{iprobe,1} = zeros(size(EXP.Bayes.Posterior0{iprobe,1}));
%                 [~,maxpos] = max(EXP.Bayes.Posterior0{iprobe,1},[],2);
%                 for tt = 1:size(EXP.Bayes.Posterior0{iprobe,1},1)
%                     EXP.Bayes.PosBias0{iprobe,1}(tt,:) = circshift(EXP.Bayes.Posterior0{iprobe,1}(tt,:),floor(Prange/2) - maxpos(tt));
%                 end
                
                Xtidx = false(nTimes,Xrange);
                for i = 1:Xrange
                    Xtidx(:,i) = res.X{ianimal,iseries} == i;
                end
                
                goodtidx = true(nTimes,1);%sum(EXP.Bayes.PosError0{iprobe,1}(:,floor(Prange/4):floor(3*Prange/4)),2)>floor(Prange/2);
                
                for g = [2 1 3]
                    res.tidx{ianimal,iseries,iprobe,g} = false(nTimes,1);
                    for cont = 1:numel(cont_list)
                        for r = 1:numel(RL_list)
                            for o = 1:numel(outcome_list)
                                tidx = EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                res.tidx{ianimal,iseries,iprobe,g}(tidx) = true;
                            end
                        end
                    end
                    
                    res.tidx0{ianimal,iseries,iprobe,g} = res.tidx{ianimal,iseries,iprobe,g};
                    res.tidx{ianimal,iseries,iprobe,g} = res.tidx{ianimal,iseries,iprobe,g} & goodtidx;
                    
                    if sum(res.tidx{ianimal,iseries,iprobe,g}) > 0
                        if strcmp(res.DecAveType,'DecPost')
                            res.PostXSum{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange);
                            res.ErrXSum{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange);
%                             res.BiasXSum{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange);
                        end
                        res.XSum{ianimal,iseries,iprobe,g} = zeros(1, Xrange);

                        res.maxdecErr{ianimal,iseries,iprobe,g} = NaN(1, Xrange);
                        res.meandecErr{ianimal,iseries,iprobe,g} = NaN(1, Xrange);
                        for i = 1:Xrange
                            res.maxdecErr{ianimal,iseries,iprobe,g}(i) = Prange/(2*pi)*sqrt(circ_mean(circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe,1}(Xtidx(:,i) & res.tidx{ianimal,iseries,iprobe,g} & goodPostidx),2*pi/Prange*res.X{ianimal,iseries}(Xtidx(:,i) & res.tidx{ianimal,iseries,iprobe,g} & goodPostidx)).^2));
                            res.meandecErr{ianimal,iseries,iprobe,g}(i) = Prange/(2*pi)*sqrt(circ_mean(circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{iprobe,1}(Xtidx(:,i) & res.tidx{ianimal,iseries,iprobe,g} & goodPostidx),2*pi/Prange*res.X{ianimal,iseries}(Xtidx(:,i) & res.tidx{ianimal,iseries,iprobe,g} & goodPostidx)).^2));
                        end
                        
                        for i = 1:Xrange
                            if strcmp(res.DecAveType,'DecPost')
                                res.PostXSum{ianimal,iseries,iprobe,g}(:,i) = nansum(EXP.Bayes.Posterior0{iprobe,1}(Xtidx(:,i) & res.tidx{ianimal,iseries,iprobe,g} & goodPostidx,:),1)';
                                res.ErrXSum{ianimal,iseries,iprobe,g}(:,i) = nansum(EXP.Bayes.PosError0{iprobe,1}(Xtidx(:,i) & res.tidx{ianimal,iseries,iprobe,g} & goodPostidx,:),1)';
%                                 res.BiasXSum{ianimal,iseries,iprobe,g}(:,i) = nansum(EXP.Bayes.PosBias0{iprobe,1}(Xtidx(:,i) & res.tidx{ianimal,iseries,iprobe,g} & goodPostidx,:),1)';
                            end
                            res.XSum{ianimal,iseries,iprobe,g}(i) = sum(Xtidx(:,i) & res.tidx{ianimal,iseries,iprobe,g} & goodPostidx);
                        end
                        X = res.X{ianimal,iseries};
                        if strcmp(res.DecAveType,'DecError')
                            Y = EXP.Bayes.MaxDecodedPosition0{iprobe,1};
                            [res.PostXSum{ianimal,iseries,iprobe,g}, ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & goodPostidx)],...
                                                                                                  [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                            res.PostXSum{ianimal,iseries,iprobe,g} = res.PostXSum{ianimal,iseries,iprobe,g}/baseline;
                            
                            Y = EXP.Bayes.DecodingError0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
                            [res.ErrXSum{ianimal,iseries,iprobe,g}, ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & goodPostidx)],...
                                                                                                  [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                            res.ErrXSum{ianimal,iseries,iprobe,g} = res.ErrXSum{ianimal,iseries,iprobe,g}/baseline;
                            
%                             Y = EXP.Bayes.DecodingBias0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
%                             [res.BiasXSum{ianimal,iseries,iprobe,g}, ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & goodPostidx)],...
%                                                                                                    [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
%                             res.BiasXSum{ianimal,iseries,iprobe,g} = res.BiasXSum{ianimal,iseries,iprobe,g}/baseline;
                        end
                        
                        idx = find(res.tidx{ianimal,iseries,iprobe,g});
                        [CVO] = crossValPartition(1:numel(idx), kfold);
                        for kiter = 1:CVO.kfold
                            if strcmp(res.DecAveType,'DecPost')
                                res.PostXSumCVO{ianimal,iseries,iprobe,g,kiter} = zeros(Prange, Xrange);
                                res.ErrXSumCVO{ianimal,iseries,iprobe,g,kiter} = zeros(Prange, Xrange);
%                                 res.BiasXSumCVO{ianimal,iseries,iprobe,g,kiter} = zeros(Prange, Xrange);
                            end
                            res.XSumCVO{ianimal,iseries,iprobe,g,kiter} = zeros(1, Xrange);
                            tidxCVO = false(size(res.tidx{ianimal,iseries,iprobe,g}));
                            tidxCVO(idx(CVO.train{kiter})) = true;
                            for i = 1:Xrange
                                if strcmp(res.DecAveType,'DecPost')
                                    res.PostXSumCVO{ianimal,iseries,iprobe,g,kiter}(:,i) = nansum(EXP.Bayes.Posterior0{iprobe,1}(Xtidx(:,i) & tidxCVO & goodPostidx,:),1)';
                                    res.ErrXSumCVO{ianimal,iseries,iprobe,g,kiter}(:,i) = nansum(EXP.Bayes.PosError0{iprobe,1}(Xtidx(:,i) & tidxCVO & goodPostidx,:),1)';
%                                     res.BiasXSumCVO{ianimal,iseries,iprobe,g,kiter}(:,i) = nansum(EXP.Bayes.PosBias0{iprobe,1}(Xtidx(:,i) & tidxCVO & goodPostidx,:),1)';
                                end
                                res.XSumCVO{ianimal,iseries,iprobe,g,kiter}(i) = sum(res.X{ianimal,iseries} == i & tidxCVO & goodPostidx);
                            end
                            
                            X = res.X{ianimal,iseries};
                            if strcmp(res.DecAveType,'DecError')
                                Y = EXP.Bayes.MaxDecodedPosition0{iprobe,1};
                                [res.PostXSumCVO{ianimal,iseries,iprobe,g,kiter}, ~, ~, ~] = smoothhist2D_AS([X(tidxCVO & goodPostidx) Y(tidxCVO & goodPostidx)],...
                                                                                                               [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                                res.PostXSumCVO{ianimal,iseries,iprobe,g,kiter} = res.PostXSumCVO{ianimal,iseries,iprobe,g,kiter}/baseline;
                                
                                Y = EXP.Bayes.DecodingError0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
                                [res.ErrXSumCVO{ianimal,iseries,iprobe,g,kiter}, ~, ~, ~] = smoothhist2D_AS([X(tidxCVO & goodPostidx) Y(tidxCVO & goodPostidx)],...
                                                                                                               [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                                res.ErrXSumCVO{ianimal,iseries,iprobe,g,kiter} = res.ErrXSumCVO{ianimal,iseries,iprobe,g,kiter}/baseline;
                                
%                                 Y = EXP.Bayes.DecodingBias0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
%                                 [res.BiasXSumCVO{ianimal,iseries,iprobe,g,kiter}, ~, ~, ~] = smoothhist2D_AS([X(tidxCVO & goodPostidx) Y(tidxCVO & goodPostidx)],...
%                                                                                                                 [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
%                                 res.BiasXSumCVO{ianimal,iseries,iprobe,g,kiter} = res.BiasXSumCVO{ianimal,iseries,iprobe,g,kiter}/baseline;
                            end
                        end


%                         res.tidxSpdEye{ianimal,iseries,iprobe,g} = false(nTimes, nSpdbins, res.nEyebins);
%                         res.PostXSumSpdEye{ianimal,iseries,iprobe,g} = zeros(nSpdbins, res.nEyebins, Prange, Xrange);
%                         res.ErrXSumSpdEye{ianimal,iseries,iprobe,g} = zeros(nSpdbins, res.nEyebins, Prange, Xrange);
%                         res.XSumSpdEye{ianimal,iseries,iprobe,g} = zeros(nSpdbins, res.nEyebins, Xrange);
%                         for ispd = 1:nSpdbins
%                             spdidx = res.Spdbin2{ianimal,iseries} == ispd;
%                             for ieye = 1:res.nEyebins
%                                 eyeidx = res.Eyebin2{ianimal,iseries} == ieye;
%                                 res.tidxSpdEye{ianimal,iseries,iprobe,g}(:,ispd,ieye) = res.tidx{ianimal,iseries,iprobe,g} & spdidx & eyeidx;
%                                 for i = 1:Xrange
%                                     res.PostXSumSpdEye{ianimal,iseries,iprobe,g}(ispd,ieye,:,i) = nansum(EXP.Bayes.Posterior0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & spdidx & eyeidx,:),1)';
%                                     res.ErrXSumSpdEye{ianimal,iseries,iprobe,g}(ispd,ieye,:,i) = circshift(res.PostXSumSpdEye{ianimal,iseries,iprobe,g}(ispd,ieye,:,i),floor(Prange/2) - i + 1,3);
%                                     res.XSumSpdEye{ianimal,iseries,iprobe,g}(ispd,ieye,i) = sum(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & spdidx & eyeidx & goodPostidx);
%                                 end
%                             end
%                         end
                        
                        res.tidxSpd{ianimal,iseries,iprobe,g} = false(nTimes, res.nSpdbins);
                        res.PostXSumSpd{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange, res.nSpdbins);
                        res.ErrXSumSpd{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange, res.nSpdbins);
%                         res.BiasXSumSpd{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange, res.nSpdbins);
                        res.XSumSpd{ianimal,iseries,iprobe,g} = zeros(1,Xrange, res.nSpdbins);
                        for ispd = 1:res.nSpdbins
                            spdidx = res.Spdbin2{ianimal,iseries} == ispd;
                            res.tidxSpd{ianimal,iseries,iprobe,g}(:,ispd) = res.tidx{ianimal,iseries,iprobe,g} & spdidx;
                            for i = 1:Xrange
                                if strcmp(res.DecAveType,'DecPost')
                                    res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,i,ispd) = nansum(EXP.Bayes.Posterior0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx,:),1)';
                                    res.ErrXSumSpd{ianimal,iseries,iprobe,g}(:,i,ispd) = nansum(EXP.Bayes.PosError0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx,:),1)';
%                                     res.BiasXSumSpd{ianimal,iseries,iprobe,g}(:,i,ispd) = nansum(EXP.Bayes.PosBias0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx,:),1)';
                                end
                                res.XSumSpd{ianimal,iseries,iprobe,g}(1,i,ispd) = sum(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx);
                            end
                            X = res.X{ianimal,iseries};
                            if strcmp(res.DecAveType,'DecError')
                                Y = EXP.Bayes.MaxDecodedPosition0{iprobe,1};
                                [res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,ispd), ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx)],...
                                                                                                                   [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                                res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,ispd) = res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,ispd)/baseline;
                                
                                Y = EXP.Bayes.DecodingError0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
                                [res.ErrXSumSpd{ianimal,iseries,iprobe,g}(:,:,ispd), ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx)],...
                                                                                                                   [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                                res.ErrXSumSpd{ianimal,iseries,iprobe,g}(:,:,ispd) = res.ErrXSumSpd{ianimal,iseries,iprobe,g}(:,:,ispd)/baseline;
                                
%                                 Y = EXP.Bayes.DecodingBias0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
%                                 [res.BiasXSumSpd{ianimal,iseries,iprobe,g}(:,:,ispd), ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & spdidx & goodPostidx)],...
%                                                                                                                     [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
%                                 res.BiasXSumSpd{ianimal,iseries,iprobe,g}(:,:,ispd) = res.BiasXSumSpd{ianimal,iseries,iprobe,g}(:,:,ispd)/baseline;
                            end
                        end
                        
                        res.tidxEye{ianimal,iseries,iprobe,g} = false(nTimes, res.nEyebins);
                        res.PostXSumEye{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange, res.nEyebins);
                        res.ErrXSumEye{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange, res.nEyebins);
                        res.BiasXSumEye{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange, res.nEyebins);
                        res.XSumEye{ianimal,iseries,iprobe,g} = zeros(1,Xrange, res.nEyebins);
                        for ieye = 1:res.nEyebins
                            eyeidx = res.Eyebin2{ianimal,iseries} == ieye;
                            res.tidxEye{ianimal,iseries,iprobe,g}(:,ieye) = res.tidx{ianimal,iseries,iprobe,g} & eyeidx;
                            for i = 1:Xrange
                                if strcmp(res.DecAveType,'DecPost')
                                    res.PostXSumEye{ianimal,iseries,iprobe,g}(:,i,ieye) = nansum(EXP.Bayes.Posterior0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx,:),1)';
                                    res.ErrXSumEye{ianimal,iseries,iprobe,g}(:,i,ieye) = nansum(EXP.Bayes.PosError0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx,:),1)';
%                                     res.BiasXSumEye{ianimal,iseries,iprobe,g}(:,i,ieye) = nansum(EXP.Bayes.PosBias0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx,:),1)';
                                end
                                res.XSumEye{ianimal,iseries,iprobe,g}(1,i,ieye) = sum(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx);
                            end
                            X = res.X{ianimal,iseries};
                            if strcmp(res.DecAveType,'DecError')
                                Y = EXP.Bayes.MaxDecodedPosition0{iprobe,1};
                                [res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,ieye), ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx)],...
                                                                                                                   [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                                res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,ieye) = res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,ieye)/baseline;
                                
                                Y = EXP.Bayes.DecodingError0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
                                [res.ErrXSumEye{ianimal,iseries,iprobe,g}(:,:,ieye), ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx)],...
                                                                                                                   [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                                res.ErrXSumEye{ianimal,iseries,iprobe,g}(:,:,ieye) = res.ErrXSumEye{ianimal,iseries,iprobe,g}(:,:,ieye)/baseline;
                                
%                                 Y = EXP.Bayes.DecodingBias0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
%                                 [res.BiasXSumEye{ianimal,iseries,iprobe,g}(:,:,ieye), ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & eyeidx & goodPostidx)],...
%                                                                                                                     [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
%                                 res.BiasXSumEye{ianimal,iseries,iprobe,g}(:,:,ieye) = res.BiasXSumEye{ianimal,iseries,iprobe,g}(:,:,ieye)/baseline;
                            end
                        end
                        
                        res.tidxVis{ianimal,iseries,iprobe,g} = false(nTimes, res.nVisbins);
                        res.PostXSumVis{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange, res.nVisbins);
                        res.ErrXSumVis{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange, res.nVisbins);
%                         res.BiasXSumVis{ianimal,iseries,iprobe,g} = zeros(Prange, Xrange, res.nVisbins);
                        res.XSumVis{ianimal,iseries,iprobe,g} = zeros(1,Xrange, res.nVisbins);
                        for ivis = 1:res.nVisbins
                            visidx = res.Visbin2{ianimal,iseries} == ivis;
                            res.tidxVis{ianimal,iseries,iprobe,g}(:,ivis) = res.tidx{ianimal,iseries,iprobe,g} & visidx;
                            for i = 1:Xrange
                                if strcmp(res.DecAveType,'DecPost')
                                    res.PostXSumVis{ianimal,iseries,iprobe,g}(:,i,ivis) = nansum(EXP.Bayes.Posterior0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx,:),1)';
                                    res.ErrXSumVis{ianimal,iseries,iprobe,g}(:,i,ivis) = nansum(EXP.Bayes.PosError0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx,:),1)';
%                                     res.BiasXSumVis{ianimal,iseries,iprobe,g}(:,i,ivis) = nansum(EXP.Bayes.PosBias0{iprobe,1}(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx,:),1)';
                                end
                                res.XSumVis{ianimal,iseries,iprobe,g}(1,i,ivis) = sum(Xtidx(:,i) &  res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx);
                            end
                            X = res.X{ianimal,iseries};
                            if strcmp(res.DecAveType,'DecError')
                                Y = EXP.Bayes.MaxDecodedPosition0{iprobe,1};
                                [res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,ivis), ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx)],...
                                                                                                                   [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                                res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,ivis) = res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,ivis)/baseline;
                                
                                Y = EXP.Bayes.DecodingError0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
                                [res.ErrXSumVis{ianimal,iseries,iprobe,g}(:,:,ivis), ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx)],...
                                                                                                                   [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                                res.ErrXSumVis{ianimal,iseries,iprobe,g}(:,:,ivis) = res.ErrXSumVis{ianimal,iseries,iprobe,g}(:,:,ivis)/baseline;
                                
%                                 Y = EXP.Bayes.DecodingBias0{iprobe,1}+floor(Prange/2);Y(Y==0) = Prange;
%                                 [res.BiasXSumVis{ianimal,iseries,iprobe,g}(:,:,ivis), ~, ~, ~] = smoothhist2D_AS([X(res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx) Y(res.tidx{ianimal,iseries,iprobe,g} & visidx & goodPostidx)],...
%                                                                                                                     [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
%                                 res.BiasXSumVis{ianimal,iseries,iprobe,g}(:,:,ivis) = res.BiasXSumVis{ianimal,iseries,iprobe,g}(:,:,ivis)/baseline;
                            end
                        end
                        
                        if res.Tsmthwin_dec <= 50
                            res.nSpdPhsbins = 3;
                            XrangePhs = 20;
                            XPhsbinsize = Xrange/XrangePhs;
                            SpdPhsbinsize = res.nSpdbins/res.nSpdPhsbins;
                            
                            res.thetaperiod{ianimal,iseries,iprobe,g} = zeros(size(res.tidx{ianimal,iseries,iprobe,g}));
                            minidx0 = find([0;abs(diff(mod(res.BayesLFPphase{ianimal,iseries,iprobe},360)))>180]);
                            thetaperiod0 = diff(minidx0);
                            thetaperiod0 = thetaperiod0(:);
                            thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
                            res.thetaperiod{ianimal,iseries,iprobe,g}(minidx0) = thetaperiod0;
                            
                            nbtimebins = 18;
                            idx = find(res.tidx0{ianimal,iseries,iprobe,g});
                            minidx = minidx0(ismember(minidx0,idx));
                            thetaperiod =  thetaperiod0(ismember(minidx0,idx));
                            POdec = zeros(res.nPhsbins,Prange,numel(minidx));
%                             POdecbias = zeros(res.nPhsbins,Prange,numel(minidx));
                            POxpos = zeros(res.nPhsbins,numel(minidx));
                            POspdbin = zeros(res.nPhsbins,numel(minidx));
                            POdecErrxpos = zeros(res.nPhsbins,numel(minidx));
%                             POdecBiasxpos = zeros(res.nPhsbins,numel(minidx));
                            POphs = zeros(res.nPhsbins,numel(minidx));
                            POtheta = 0;
                            for tt = 1:numel(minidx)
                                if minidx(tt) > 1 && minidx(tt)+nbtimebins-1 < size(EXP.Bayes.PosError0{iprobe,1},1)
                                    if strcmp(res.DecAveType,'DecPost')
                                        POmap = EXP.Bayes.PosError0{iprobe}(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),:);
                                        POdec(:,:,tt) = interp1((0:nbtimebins-1)/thetaperiod(tt),POmap,0:1/res.nPhsbins:(1-(1/res.nPhsbins)));
%                                         POmapbias = EXP.Bayes.PosBias0{iprobe}(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),:);
%                                         POdecbias(:,:,tt) = interp1((0:nbtimebins-1)/thetaperiod(tt),POmapbias,0:1/res.nPhsbins:(1-(1/res.nPhsbins)));
                                    end
                                    POthetavec = EXP.data.es.LFPtheta(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),min(end,EXP.Bayes.thetaChannel));
                                    POtheta = POtheta + (interp1((0:nbtimebins-1)/thetaperiod(tt),POthetavec,0:1/res.nPhsbins:(1-(1/res.nPhsbins))))';
                                    POx = EXP.Bayes.X0{1}(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end));
                                    POx = unwrap(POx/(Xrange)*2*pi)*Xrange/(2*pi);
                                    POxpos(:,tt) = mod(interp1((0:nbtimebins-1)/thetaperiod(tt),POx,0:1/res.nPhsbins:(1-(1/res.nPhsbins))),Xrange);
                                    POspd = res.Spdbin2{ianimal,iseries}(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),:);
                                    POspdbin(:,tt) = interp1((0:nbtimebins-1)/thetaperiod(tt),POspd,0:1/res.nPhsbins:(1-(1/res.nPhsbins)),'nearest');
                                    
                                    POphs(:,tt) = 1:res.nPhsbins;
                                    if strcmp(res.DecAveType,'DecError')
                                        POdecErrx = EXP.Bayes.DecodingError0{iprobe,1}(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end))+floor(Prange/2);
                                        POdecErrx = unwrap(POdecErrx/(Xrange)*2*pi)*Xrange/(2*pi);
                                        POdecErrxpos(:,tt) = mod(interp1((0:nbtimebins-1)/thetaperiod(tt),POdecErrx,0:1/res.nPhsbins:(1-(1/res.nPhsbins))),Xrange);
                                        POdecErrxpos(:,tt) = floor(POdecErrxpos(:,tt))+1;
                                        POdecErrxpos(POdecErrxpos(:,tt)>Prange,tt) = Prange;
%                                         POdecBiasx = EXP.Bayes.DecodingBias0{iprobe,1}(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end))+floor(Prange/2);
%                                         POdecBiasx = unwrap(POdecBiasx/(Xrange)*2*pi)*Xrange/(2*pi);
%                                         POdecBiasxpos(:,tt) = mod(interp1((0:nbtimebins-1)/thetaperiod(tt),POdecBiasx,0:1/res.nPhsbins:(1-(1/res.nPhsbins))),Xrange);
%                                         POdecBiasxpos(:,tt) = floor(POdecBiasxpos(:,tt))+1;
%                                         POdecBiasxpos(POdecBiasxpos(:,tt)>Prange,tt) = Prange;
                                    end
                                end
                            end
                            
                            res.ErrXSumPhsAll{ianimal,iseries,iprobe,g} = zeros(Prange, res.nPhsbins, XrangePhs);
%                             res.BiasXSumPhsAll{ianimal,iseries,iprobe,g} = zeros(Prange, res.nPhsbins, XrangePhs);
                            
                            res.XSumPhsAll{ianimal,iseries,iprobe,g} = zeros(Prange, res.nPhsbins, XrangePhs);
                            for i = 1:XrangePhs
                                for iphs = 1:res.nPhsbins
                                    ttidx = find(ismember(min(floor((POxpos(iphs,:)-1)/XPhsbinsize)+1,XrangePhs),i));
                                    if strcmp(res.DecAveType,'DecPost')
                                        res.ErrXSumPhsAll{ianimal,iseries,iprobe,g}(:,iphs,i)= nansum(POdec(iphs,:,ttidx),3)';
%                                         res.BiasXSumPhsAll{ianimal,iseries,iprobe,g}(:,iphs,i)= nansum(POdecbias(iphs,:,ttidx),3)';
                                    end
                                    res.XSumPhsAll{ianimal,iseries,iprobe,g}(:,iphs,i) = nansum(~isnan(POdec(iphs,:,ttidx)),3)';
                                end
                                if strcmp(res.DecAveType,'DecError')
                                    [res.ErrXSumPhsAll{ianimal,iseries,iprobe,g}(:,:,i), ~, ~, ~] = smoothhist2D_AS([POphs(ismember(min(floor((POxpos-1)/XPhsbinsize)+1,XrangePhs),i) & ~isnan(POdecErrxpos)) POdecErrxpos(ismember(min(floor((POxpos-1)/XPhsbinsize)+1,XrangePhs),i) & ~isnan(POdecErrxpos))],...
                                                                                                                       [NaN NaN], [res.nPhsbins Prange], 1:res.nPhsbins, 1:Prange, true, true);
                                    res.ErrXSumPhsAll{ianimal,iseries,iprobe,g}(:,:,i) = res.ErrXSumPhsAll{ianimal,iseries,iprobe,g}(:,:,i)/baseline;
                                    
%                                     [res.BiasXSumPhsAll{ianimal,iseries,iprobe,g}(:,:,i), ~, ~, ~] = smoothhist2D_AS([POphs(ismember(min(floor((POxpos-1)/XPhsbinsize)+1,XrangePhs),i) & ~isnan(POdecBiasxpos)) POdecBiasxpos(ismember(min(floor((POxpos-1)/XPhsbinsize)+1,XrangePhs),i) & ~isnan(POdecBiasxpos))],...
%                                                                                                                         [NaN NaN], [res.nPhsbins Prange], 1:res.nPhsbins, 1:Prange, true, true);
%                                     res.BiasXSumPhsAll{ianimal,iseries,iprobe,g}(:,:,i) = res.BiasXSumPhsAll{ianimal,iseries,iprobe,g}(:,:,i)/baseline;
                                end
                            end
                            
                            [CVO] = crossValPartition(1:size(POdec,3), kfold);
                            for kiter = 1:CVO.kfold
                                tidxCVO = false(1,size(POdec,3));
                                tidxCVO(CVO.train{kiter}) = true;
                                res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter} = zeros(Prange, res.nPhsbins, XrangePhs);
%                                 res.BiasXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter} = zeros(Prange, res.nPhsbins, XrangePhs);
                                
                                res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter} = zeros(Prange, res.nPhsbins, XrangePhs);
                                POphs_iter = POphs(:,tidxCVO);
                                POxpos_iter = POxpos(:,tidxCVO);
                                POdecErrxpos_iter = POdecErrxpos(:,tidxCVO);
%                                 POdecBiasxpos_iter = POdecBiasxpos(:,tidxCVO);
                                for i = 1:XrangePhs
                                    for iphs = 1:res.nPhsbins
                                        ttidx = find(ismember(min(floor((POxpos(iphs,:)-1)/XPhsbinsize)+1,XrangePhs),i) & tidxCVO);
                                        if strcmp(res.DecAveType,'DecPost')
                                            res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,iphs,i) = nansum(POdec(iphs,:,ttidx),3)';
%                                             res.BiasXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,iphs,i) = nansum(POdecbias(iphs,:,ttidx),3)';
                                        end
                                        res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,iphs,i) = nansum(~isnan(POdec(iphs,:,ttidx)),3)';
                                    end
                                    if strcmp(res.DecAveType,'DecError')
                                        [res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,:,i), ~, ~, ~] = smoothhist2D_AS([POphs_iter(ismember(min(floor((POxpos_iter-1)/XPhsbinsize)+1,XrangePhs),i) & ~isnan(POdecErrxpos_iter)) POdecErrxpos_iter(ismember(min(floor((POxpos_iter-1)/XPhsbinsize)+1,XrangePhs),i) & ~isnan(POdecErrxpos_iter))],...
                                                                                                                                    [NaN NaN], [res.nPhsbins Prange], 1:res.nPhsbins, 1:Prange, true, true);
                                        res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,:,i) = res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,:,i)/baseline;
                                        
%                                         [res.BiasXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,:,i), ~, ~, ~] = smoothhist2D_AS([POphs_iter(ismember(min(floor((POxpos_iter-1)/XPhsbinsize)+1,XrangePhs),i) & ~isnan(POdecBiasxpos_iter)) POdecBiasxpos_iter(ismember(min(floor((POxpos_iter-1)/XPhsbinsize)+1,XrangePhs),i) & ~isnan(POdecBiasxpos_iter))],...
%                                                                                                                                      [NaN NaN], [res.nPhsbins Prange], 1:res.nPhsbins, 1:Prange, true, true);
%                                         res.BiasXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,:,i) = res.BiasXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,:,i)/baseline;
                                    end
                                end
                            end
                            
                            for ispd = 1:res.nSpdPhsbins
                                res.ErrXSumPhs{ianimal,iseries,iprobe,ispd,g} = zeros(Prange, res.nPhsbins, XrangePhs);
%                                 res.BiasXSumPhs{ianimal,iseries,iprobe,ispd,g} = zeros(Prange, res.nPhsbins, XrangePhs);
                                
                                res.XSumPhs{ianimal,iseries,iprobe,ispd,g} = zeros(Prange, res.nPhsbins, XrangePhs);
                                for i = 1:XrangePhs
                                    for iphs = 1:res.nPhsbins
                                        ttidx = find(ismember(min(floor(POxpos(iphs,:)/XPhsbinsize)+1,XrangePhs),i) & ismember(floor((POspdbin(iphs,:)-1)/SpdPhsbinsize)+1,ispd));
                                        if strcmp(res.DecAveType,'DecPost')
                                            res.ErrXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,iphs,i) = nansum(POdec(iphs,:,ttidx),3);
%                                             res.BiasXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,iphs,i) = nansum(POdecbias(iphs,:,ttidx),3);
                                        end
                                        res.XSumPhs{ianimal,iseries,iprobe,ispd,g}(:,iphs,i) = nansum(~isnan(POdec(iphs,:,ttidx)),3);
                                    end
                                    if strcmp(res.DecAveType,'DecError')
                                        [res.ErrXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,:,i), ~, ~, ~] = smoothhist2D_AS([POphs(ismember(min(floor((POxpos-1)/XPhsbinsize)+1,XrangePhs),i) & ismember(floor((POspdbin-1)/SpdPhsbinsize)+1,ispd) & ~isnan(POdecErrxpos)) POdecErrxpos(ismember(min(floor((POxpos-1)/XPhsbinsize)+1,XrangePhs),i) & ismember(floor((POspdbin-1)/SpdPhsbinsize)+1,ispd) & ~isnan(POdecErrxpos))],...
                                                                                                                             [NaN NaN], [res.nPhsbins Prange], 1:res.nPhsbins, 1:Prange, true, true);
                                        res.ErrXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,:,i) = res.ErrXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,:,i)/baseline;
                                        
%                                         [res.BiasXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,:,i), ~, ~, ~] = smoothhist2D_AS([POphs(ismember(min(floor((POxpos-1)/XPhsbinsize)+1,XrangePhs),i) & ismember(floor((POspdbin-1)/SpdPhsbinsize)+1,ispd) & ~isnan(POdecBiasxpos)) POdecBiasxpos(ismember(min(floor((POxpos-1)/XPhsbinsize)+1,XrangePhs),i) & ismember(floor((POspdbin-1)/SpdPhsbinsize)+1,ispd) & ~isnan(POdecBiasxpos))],...
%                                                                                                                               [NaN NaN], [res.nPhsbins Prange], 1:res.nPhsbins, 1:Prange, true, true);
%                                         res.BiasXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,:,i) = res.BiasXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,:,i)/baseline;
                                    end
                                end
                            end
                            POdec = [];
                            POxpos = [];
                            POspdbin = [];
                            POdecErrxpos = [];
                            POdecBiasxpos = [];
                            POphs = [];
                            
                            
%                             res.nSpdPhsbins = 3;
%                             phsbinsize = 360/res.nPhsbins;
%                             XrangePhs = 20;
%                             XPhsbinsize = Xrange/XrangePhs;
%                             X0phsbinned = min(floor((EXP.Bayes.X0{1}-1)/XPhsbinsize)+1,XrangePhs);
%                             Thetaphsbinned = min(floor(res.BayesLFPphase{ianimal,iseries,1}/phsbinsize) + 1,res.nPhsbins);
%                             SpdPhsbinsize = res.nSpdbins/res.nSpdPhsbins;
%                             SpdPhsbinned = min(floor((res.Spdbin2{ianimal,iseries}-1)/SpdPhsbinsize) + 1,res.nSpdPhsbins);
%                             
%                             Xtidxphs = false(nTimes,XrangePhs);
%                             for i = 1:XrangePhs
%                                 Xtidxphs(:,i) = ismember(X0phsbinned,i);
%                             end
%                             Thetatidxphs = false(nTimes,res.nPhsbins);
%                             for iphs = 1:res.nPhsbins
%                                 Thetatidxphs(:,iphs) = ismember(Thetaphsbinned,iphs);
%                             end
%                             Spdtidxphs = false(nTimes,res.nSpdPhsbins);
%                             for ispd = 1:res.nSpdPhsbins
%                                 Spdtidxphs(:,ispd) = ismember(SpdPhsbinned,ispd);
%                             end
%                             
%                             res.ErrXSumPhsAll{ianimal,iseries,iprobe,g} = zeros(Prange, XrangePhs, res.nPhsbins);
%                             res.XSumPhsAll{ianimal,iseries,iprobe,g} = zeros(Prange, XrangePhs, res.nPhsbins);
%                             for i = 1:XrangePhs
%                                 for iphs = 1:res.nPhsbins
%                                     ttidx = find(res.tidx{ianimal,iseries,iprobe,g} & Xtidxphs(:,i) & Thetatidxphs(:,iphs));
%                                     res.ErrXSumPhsAll{ianimal,iseries,iprobe,g}(:,i,iphs) = nansum(EXP.Bayes.PosError0{iprobe,1}(ttidx,:),1);
%                                     res.XSumPhsAll{ianimal,iseries,iprobe,g}(:,i,iphs) = nansum(~isnan(EXP.Bayes.PosError0{iprobe,1}(ttidx,:)),1);
%                                 end
%                             end
%                             
%                             idx = find(res.tidx{ianimal,iseries,iprobe,g});
%                             [CVO] = crossValPartition(1:numel(idx), kfold);
%                             for kiter = 1:CVO.kfold
%                                 res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter} = zeros(Prange, XrangePhs, res.nPhsbins);
%                                 res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter} = zeros(Prange, XrangePhs, res.nPhsbins);
%                                 tidxCVO = false(size(res.tidx{ianimal,iseries,iprobe,g}));
%                                 tidxCVO(idx(CVO.train{kiter})) = true;
%                                 for i = 1:XrangePhs
%                                     for iphs = 1:res.nPhsbins
%                                         ttidx = find(res.tidx{ianimal,iseries,iprobe,g} & Xtidxphs(:,i) & Thetatidxphs(:,iphs) & tidxCVO);
%                                         res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,i,iphs) = nansum(EXP.Bayes.PosError0{iprobe,1}(ttidx,:),1);
%                                         res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter}(:,i,iphs) = nansum(~isnan(EXP.Bayes.PosError0{iprobe,1}(ttidx,:)),1);
%                                     end
%                                 end
%                             end
%                             
%                             for ispd = 1:res.nSpdPhsbins
%                                 res.ErrXSumPhs{ianimal,iseries,ispd,iprobe,g} = zeros(Prange, XrangePhs, res.nPhsbins);
%                                 res.XSumPhs{ianimal,iseries,ispd,iprobe,g} = zeros(Prange, XrangePhs, res.nPhsbins);
%                                 for i = 1:XrangePhs
%                                     for iphs = 1:res.nPhsbins
%                                         ttidx = find(res.tidx{ianimal,iseries,iprobe,g} & Xtidxphs(:,i) & Thetatidxphs(:,iphs) & Spdtidxphs(:,ispd));
%                                         res.ErrXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,i,iphs) = nansum(EXP.Bayes.PosError0{iprobe,1}(ttidx,:),1);
%                                         res.XSumPhs{ianimal,iseries,iprobe,ispd,g}(:,i,iphs) = nansum(~isnan(EXP.Bayes.PosError0{iprobe,1}(ttidx,:)),1);
%                                     end
%                                 end
%                             end
% 
%                             minidx0 = find([0;abs(diff(mod(res.BayesLFPphase{ianimal,iseries,1},360)))>180]);
%                             thetaperiod0 = diff(minidx0);
%                             thetaperiod0 = thetaperiod0(:);
%                             thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
%                             res.thetaperiod{ianimal,iseries,iprobe,g}(minidx0) = thetaperiod0;




%                             SpdPhsbinsize = 1;
%                             res.nSpdPhsbins = floor(res.nSpdbins / SpdPhsbinsize);
%                             XPhsbinsize = 5;
%                             XrangePhs = floor(Xrange / XPhsbinsize);
%                             res.tidxPhs{ianimal,iseries,iprobe,g} = false(nTimes, res.nPhsbins);
%                             for ispd = 1:res.nSpdPhsbins
%                                 res.ErrXSumPhs{ianimal,iseries,ispd,iprobe,g} = zeros(Prange, XrangePhs, res.nPhsbins);
%                                 res.XSumPhs{ianimal,iseries,ispd,iprobe,g} = zeros(Prange, XrangePhs, res.nPhsbins);
%                             end
%                             res.thetaperiod{ianimal,iseries,iprobe,g} = zeros(1,nTimes);
%                             
%                             nbtimebins = 18;
%                             minidx0 = find([0;abs(diff(mod(res.BayesLFPphase{ianimal,iseries,1},360)))>180]);
%                             thetaperiod0 = diff(minidx0);
%                             thetaperiod0 = thetaperiod0(:);
%                             thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
%                             res.thetaperiod{ianimal,iseries,iprobe,g}(minidx0) = thetaperiod0;
%                             
%                             idx = find(res.tidx{ianimal,iseries,iprobe,g});
%                             minidx = minidx0(ismember(minidx0,idx));
%                             thetaperiod =  thetaperiod0(ismember(minidx0,idx));
%                             POdec = zeros(res.nPhsbins,Prange,numel(minidx));
%                             POxpos = zeros(res.nPhsbins,numel(minidx));
%                             POspdbin = zeros(res.nPhsbins,numel(minidx));
%                             POtheta = 0;
%                             goodidx = single(~isnan(sum(EXP.Bayes.Posterior0{iprobe,1},2)));
%                             for tt = 1:numel(minidx)
%                                 if minidx(tt) > 1 && minidx(tt)+nbtimebins-1 < size(EXP.Bayes.PosError0{iprobe},1) %&& sum(es.badlick(max(1,minidx(tt)):min(minidx(tt)+20*nbphsbins,end))) > 0 %&& sum(es.goodlick(max(1,minidx(tt)-24*nbphsbins):min(minidx(tt),end))) > 0
%                                     POmap = EXP.Bayes.PosError0{iprobe}(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),:);
%                                     POdec(:,:,tt) = interp1((0:nbtimebins-1)/thetaperiod(tt),POmap,0:1/res.nPhsbins:(1-(1/res.nPhsbins)));%POmap;
%                                     POthetavec = EXP.data.es.LFPtheta(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),min(end,EXP.Bayes.thetaChannel));
%                                     POtheta = POtheta + (interp1((0:nbtimebins-1)/thetaperiod(tt),POthetavec,0:1/res.nPhsbins:(1-(1/res.nPhsbins))))';
%                                     POx = EXP.Bayes.X0{1}(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),:);
%                                     POx = unwrap(POx/(Xrange)*2*pi)*Xrange/(2*pi);
%                                     POxpos(:,tt) = mod(interp1((0:nbtimebins-1)/thetaperiod(tt),POx,0:1/res.nPhsbins:(1-(1/res.nPhsbins))),Xrange);%POmap;
%                                     POspd = res.Spdbin2{ianimal,iseries}(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),:);
%                                     POspdbin(:,tt) = interp1((0:nbtimebins-1)/thetaperiod(tt),POspd,0:1/res.nPhsbins:(1-(1/res.nPhsbins)),'nearest');
%                                 end
%                             end
%                             POdecall = POdec;
%                             for ispd = 1:res.nSpdPhsbins
%                                 res.ErrXSumPhsAll{ianimal,iseries,iprobe,ispd,g} = nansum(POdec,3)';
%                                 res.XSumPhsAll{ianimal,iseries,iprobe,ispd,g} = nansum(~isnan(POdec),3)';
%                             end
%                             
%                             allcycles = 1:size(POdecall,3);
%                             for kiter = 1:kfold
%                                 idxtest = ((kiter-1)*floor(numel(allcycles)/kfold)+1):(kiter)*floor(numel(allcycles)/kfold);
%                                 idxtrain = allcycles(~ismember(allcycles,idxtest));
%                                 res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter} = nansum(POdecall(:,:,idxtrain),3)';
%                                 res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter} = nansum(~isnan(POdecall(:,:,idxtrain)),3)';
%                             end
%                             
%                             for ispd = 1:res.nSpdPhsbins
%                                 for i = 1:XrangePhs
%                                     for iphs = 1:res.nPhsbins
%                                         ttidx = find(ismember(min(floor(POxpos(iphs,:)/XPhsbinsize)+1,XrangePhs),i) & ismember(floor((POspdbin(iphs,:)-1)/SpdPhsbinsize)+1,ispd));
%                                         res.ErrXSumPhs{ianimal,iseries,iprobe,ispd,g}(:,i,iphs) = nansum(POdec(iphs,:,ttidx),3);
%                                         res.XSumPhs{ianimal,iseries,iprobe,ispd,g}(:,i,iphs) = nansum(~isnan(POdec(iphs,:,ttidx)),3);
%                                     end
%                                 end
%                             end
%                             POdec = [];
%                             POxpos = [];
%                             POspdbin = [];
                        end
                    end
                end
            end
        end
    end
end
save(savedfilename_res, 'res','-v7.3');
popres_sessions = PopBayesAnalysis2(res,false,'sessions',false);
save(savedfilename_popres_sessions, 'popres_sessions','-v7.3');
popres_trials = PopBayesAnalysis2(res,false,'trials',false);
save(savedfilename_popres_trials, 'popres_trials','-v7.3');
end