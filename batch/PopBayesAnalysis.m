function  popres = PopBayesAnalysis(res,batch2p,FshowsessionFig,spdbins,nSpdbins,phsbins,nPhsbins)
if nargin < 2
    batch2p = [];
end
if nargin < 3
    FshowsessionFig = [];
end
if nargin < 4
    nSpdbins = 5;
    spdbins = 1:5;
end
if isempty(nSpdbins)
    nSpdbins = 5;
    spdbins = 1:5;
end
if nargin < 6
    phsbins = 1;
    nPhsbins = 1;
end
cl{1} = 'c';
cl{2} = 'k';
cl{3} = 'm';

if isempty(FshowsessionFig)
    figdirname = 'C:\Users\experiment\Desktop\figure_temp';
    prompt = {'Save single session fig?';'figure directory'};
    dlg_title = 'Parameters';
    num_lines = 1;
    def = {'0';figdirname};
    nruns = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(nruns)
        FshowsessionFig = logical(str2num(nruns{1}));
        figdirname = nruns{2};
        if ~isdir(figdirname)
            mkdir(figdirname)
        end
    end
end

lambdaSmooth = 1;
if isempty(batch2p)
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

    strlistvarname = {'V1medial','V1lateral','PPC', 'AL'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'area', 'InitialValue', 1);
    area_str = cell2mat(strlistvarname(varnamesel));
    nProbe = 1;
    popres.probestr{1} = area_str;
    popres.probestr{2} = 'none';
elseif ~batch2p
    expt = getExperimentList;
    area_str = 'CA1V1';
    nProbe = 2;
    popres.probestr{1} = 'CA1';
    popres.probestr{2} = 'V1';
end
contval = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];%[0.2 0.3 0.4];%[0.8 0.9];%

maxtol_post = 1;maxtol_meanpost = 1;%0.1;
maxtol_decmax = 1;maxtol_meandecmax = 1;
maxtol_decave = 1;maxtol_meandecave = 1;%0.1;%1;%
amp_th_post = 0;amp_th_meanpost = 0;%0;
amp_th_decmax = 0;amp_th_meandecmax = 0;%0;
amp_th_decave = 0;amp_th_meandecave = 0;%0;
popres.dx = 0.1;

nanimal = numel(expt);
ngain = 3;
popres.PostXAll = cell(nProbe,ngain);
popres.Xsum = cell(nProbe,ngain);
popres.X = cell(nProbe,ngain);
popres.Xpred_ave = cell(nProbe,ngain);
popres.Xpred_max = cell(nProbe,ngain);
popres.DistriXaveAll = cell(nProbe,ngain);
popres.DistriXmaxAll = cell(nProbe,ngain);
popres.speed = cell(nProbe,ngain);
popres.eyeXpos = cell(nProbe,ngain);
popres.licks = cell(nProbe,ngain);

popres.s_PostXAll = cell(nProbe,ngain);
popres.s_Xsum = cell(nProbe,ngain);
popres.s_X = cell(nProbe,ngain);
popres.s_trialID = cell(nProbe,ngain);
popres.s_Xpred_ave = cell(nProbe,ngain);
popres.s_Xpred_max = cell(nProbe,ngain);
popres.s_DistriXaveAll = cell(nProbe,ngain);
popres.s_DistriXmaxAll = cell(nProbe,ngain);
popres.s_meanErrX_ave = cell(nProbe,ngain);
popres.s_meanErrX_max = cell(nProbe,ngain);
popres.s_meanErrX_aveSeg = cell(nProbe,ngain);
popres.s_speed = cell(nProbe,ngain);
popres.s_licks = cell(nProbe,ngain);
popres.s_meanlickX = cell(nProbe,ngain);
popres.s_correcttrials = cell(nProbe,ngain);
popres.s_incorrecttrials = cell(nProbe,ngain);
popres.performance = cell(nProbe,ngain);
for iprobe = 1:nProbe
    for g = 1:ngain
        popres.PostXAll{iprobe,g} = 0;
        popres.Xsum{iprobe,g} = 0;
        popres.X{iprobe,g} = [];
        popres.Xpred_ave{iprobe,g} = [];
        popres.Xpred_max{iprobe,g} = [];
        popres.speed{iprobe,g} = [];
        popres.eyeXpos{iprobe,g} = [];
        popres.licks{iprobe,g} = [];
        for ianimal = 1:nanimal
            for iseries = 1:numel(expt(ianimal).series)
                if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                    if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
                        popres.s_PostXAll{ianimal,iseries,iprobe,g} = 0;
                        popres.s_Xsum{ianimal,iseries,iprobe,g} = 0;
                        popres.s_X{ianimal,iseries,iprobe,g} = [];
                        popres.s_trialID{ianimal,iseries,iprobe,g} = [];
                        popres.s_Xpred_ave{ianimal,iseries,iprobe,g} = [];
                        popres.s_Xpred_max{ianimal,iseries,iprobe,g} = [];
                        popres.s_meanErrX_ave{ianimal,iseries,iprobe,g} = [];
                        popres.s_meanErrX_aveSeg{ianimal,iseries,iprobe,g} = [];
                        popres.s_meanErrX_max{ianimal,iseries,iprobe,g} = [];
                        popres.s_speed{ianimal,iseries,iprobe,g} = [];
                        popres.s_licks{ianimal,iseries,iprobe,g} = [];
                        popres.s_meanlickX{ianimal,iseries,iprobe,g} = [];
                        popres.s_correcttrials{ianimal,iseries,iprobe,g} = [];
                        popres.s_incorrecttrials{ianimal,iseries,iprobe,g} = [];
                    end
                end
            end
            popres.performance{ianimal,iprobe,g} = [];
        end
    end
end

for ianimal = 1:nanimal
    for iseries = 1:numel(expt(ianimal).series)
        if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
            contref = find(ismember(res.contrastVal{ianimal,iseries}, [0.5 0.6 0.7]));
            RLref = find(ismember(res.roomlengthVal{ianimal,iseries}, [1]));
            outcomeref = find(ismember(res.outcomeVal{ianimal,iseries}, [2]));
            gainref = 2;
            if ~isempty(res.ballspeed{ianimal,iseries})
                speed = NaN(size(res.ballspeed{ianimal,iseries}));
                sampleRate = 60;
                if isfield(res,'Tsmthwin')
                    SpdSmthWin = res.Tsmthwin;
                else
                    SpdSmthWin = 150;
                end
                speed(~isnan(res.ballspeed{ianimal,iseries})) = smthInTime(res.ballspeed{ianimal,iseries}(~isnan(res.ballspeed{ianimal,iseries})), sampleRate, SpdSmthWin, 'same', [], 'boxcar_centered');
%                 speed(~isnan(res.trajspeed{ianimal,iseries})) = smthInTime(res.trajspeed{ianimal,iseries}(~isnan(res.trajspeed{ianimal,iseries})), sampleRate, SpdSmthWin, 'same', [], 'boxcar_centered');
                
                itraj = res.X{ianimal,iseries};% floor(X/(max(round(X))/obj.Bayes.numBins))+1;
                
                ntrajbins = max(itraj);
                spdquantilelim = zeros(ntrajbins,2);
                res.Spdbin2{ianimal,iseries} = NaN(size(res.ballspeed{ianimal,iseries}));
                try
                    idxref = res.tidx{ianimal,iseries,1,contref,gainref,RLref,outcomeref};
                catch
                    keyboard
                end
                if nSpdbins > 1
                    for spd = 1:nSpdbins
                        for xx = 1:ntrajbins
                            spdquantilelim(xx,1) = quantile(speed(idxref & itraj == xx),max(0,(spd-1)/nSpdbins));
                            spdquantilelim(xx,2) = quantile(speed(idxref & itraj == xx),min(1,(spd)/nSpdbins));
                        end
                        res.Spdbin2{ianimal,iseries}(speed >= spdquantilelim(itraj,1) & speed < spdquantilelim(itraj,2)) = spd;
                        if spd == 1
                            res.Spdbin2{ianimal,iseries}(speed <= spdquantilelim(itraj,1)) = spd;
                        end
                        if spd == nSpdbins
                            res.Spdbin2{ianimal,iseries}(speed >= spdquantilelim(itraj,2)) = spd;
                        end
                    end
                else
                    res.Spdbin2{ianimal,iseries} = ones(size(res.Spdbin2{ianimal,iseries}));
                end
            end
        end
    end
end

for ianimal = 1:nanimal
    for iseries = 1:numel(expt(ianimal).series)
        if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
            contref = find(ismember(res.contrastVal{ianimal,iseries}, [0.5 0.6 0.7]));
            RLref = find(ismember(res.roomlengthVal{ianimal,iseries}, [1]));
            outcomeref = find(ismember(res.outcomeVal{ianimal,iseries}, [2]));
            gainref = 2;
            if isfield(res,'LFPphase')
                if ~isempty(res.LFPphase{ianimal,iseries})
                    phs = NaN(size(res.LFPphase{ianimal,iseries}));
                    sampleRate = 60;
                    if isfield(res,'Tsmthwin')
                        SpdSmthWin = res.Tsmthwin;
                    else
                        SpdSmthWin = 300;%150;
                    end
                    phs = res.LFPphase{ianimal,iseries};
                    %                 phs(~isnan(res.LFPphase{ianimal,iseries})) = smthInTime(res.LFPphase{ianimal,iseries}(~isnan(res.LFPphase{ianimal,iseries})), sampleRate, SpdSmthWin, 'same', [], 'boxcar_centered');
                    
                    res.Phsbin{ianimal,iseries} = NaN(size(res.LFPphase{ianimal,iseries}));
                    if nPhsbins > 1
                        for iphsbin = 1:nPhsbins
                            res.Phsbin{ianimal,iseries}(mod(phs,360) >= (iphsbin-1)*360/nPhsbins & mod(phs,360) < iphsbin*360/nPhsbins) = iphsbin;
                        end
                    else
                        res.Phsbin{ianimal,iseries} = ones(size(res.Phsbin{ianimal,iseries}));
                    end
                end
            else
                res.Phsbin{ianimal,iseries} = ones(size(res.traj{ianimal,iseries}));
            end
        end
    end
end

if FshowsessionFig
    f_distriave = figure('Name','DistriAve');
    f_distrimax = figure('Name','DistriMax');
    f_post = figure('Name','Post');
end

outcomeVal = 2;%0;%[3];%[4];%[0 1 3 4]
for iprobe = 1:nProbe
    for ianimal = 1:nanimal
        disp(ianimal)
        nseries = 0;
        for iseries = 1:numel(expt(ianimal).series)
            if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
                    nseries = nseries + 1;
                end
            end
        end
        iseries_valid = 0;
        for iseries = 1:numel(expt(ianimal).series)
            if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
                    iseries_valid = iseries_valid + 1;
                    cont_list = find(ismember(res.contrastVal{ianimal,iseries}, contval));
                    RL_list = find(ismember(res.roomlengthVal{ianimal,iseries}, [1]));
                    outcome_list = find(ismember(res.outcomeVal{ianimal,iseries}, [outcomeVal]));
                    %                 speed(~isnan(res.trajspeed{ianimal,iseries})) = smthInTime(res.trajspeed{ianimal,iseries}(~isnan(res.trajspeed{ianimal,iseries})), sampleRate, SpdSmthWin, 'same', [], 'boxcar_centered');
                    speed(~isnan(res.ballspeed{ianimal,iseries})) = smthInTime(res.ballspeed{ianimal,iseries}(~isnan(res.ballspeed{ianimal,iseries})), sampleRate, SpdSmthWin, 'same', [], 'boxcar_centered');
                    
                    Xrange = max(res.X{ianimal,iseries}(~isnan(res.X{ianimal,iseries})));
                    
%                     res.meanErrtrial{ianimal,iseries,iprobe} = NaN(size(res.ballspeed{ianimal,iseries}));
%                     XErr = res.Xpred_ave{ianimal,iseries,iprobe} - res.X{ianimal,iseries};
%                     XErr(XErr > 50) = XErr(XErr > 50) - 100;
%                     XErr(XErr < -50) = XErr(XErr < -50) + 100;
%                     trialIDs = unique(res.trialID{ianimal,iseries});
%                     for tt = 1:numel(trialIDs)
%                         res.meanErrtrial{ianimal,iseries,iprobe}(res.trialID{ianimal,iseries} == trialIDs(tt)) = nanmean(abs(XErr(res.trialID{ianimal,iseries} == trialIDs(tt))));
%                     end
                    
                    for g = [2 1 3]
                        if ((iprobe == 1 && expt(ianimal).goodCA1dec{g}{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1dec{g}{iseries} == 1))
                            for cont = 1:numel(cont_list)
                                for r = 1:numel(RL_list)
                                    for o = 1:numel(outcome_list)
                                        if ~isempty(res.PostXSum{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)})
                                            tidx = res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)} & ismember(res.Spdbin2{ianimal,iseries},spdbins) & ismember(res.Phsbin{ianimal,iseries},phsbins);% & ismember(res.trialgainchange{ianimal,iseries}, [1:5]);
%                                             tidx = tidx & res.meanErrtrial{ianimal,iseries,iprobe} < 20;
                                            
                                            popres.speed{iprobe,g} = [popres.speed{iprobe,g} ; speed(tidx)];
                                            popres.eyeXpos{iprobe,g} = [popres.eyeXpos{iprobe,g} ; res.eyeXpos{ianimal,iseries}(tidx)];
                                            popres.licks{iprobe,g} = [popres.licks{iprobe,g} ; res.firstgoodlick{ianimal,iseries}(tidx)];
                                            
                                            popres.Xpred_ave{iprobe,g} = [popres.Xpred_ave{iprobe,g} ; res.Xpred_ave{ianimal,iseries,iprobe}(tidx)];
                                            popres.Xpred_max{iprobe,g} = [popres.Xpred_max{iprobe,g} ; res.Xpred_max{ianimal,iseries,iprobe}(tidx)];
                                            popres.X{iprobe,g} = [popres.X{iprobe,g} ; res.X{ianimal,iseries}(tidx)];%res.traj{ianimal,iseries}(tidx)];%
                                            popres.PostXAll{iprobe,g} = popres.PostXAll{iprobe,g} + res.PostXSum{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)};
                                            popres.Xsum{iprobe,g} = popres.Xsum{iprobe,g} + res.XSum{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)};
                                            
                                            popres.s_speed{ianimal,iseries,iprobe,g} = [popres.s_speed{ianimal,iseries,iprobe,g} ; speed(tidx)];
                                            popres.s_licks{ianimal,iseries,iprobe,g} = [popres.s_licks{ianimal,iseries,iprobe,g} ; res.firstgoodlick{ianimal,iseries}(tidx)];
                                            
                                            popres.s_Xpred_ave{ianimal,iseries,iprobe,g} = [popres.s_Xpred_ave{ianimal,iseries,iprobe,g} ; res.Xpred_ave{ianimal,iseries,iprobe}(tidx)];
                                            popres.s_Xpred_max{ianimal,iseries,iprobe,g} = [popres.s_Xpred_max{ianimal,iseries,iprobe,g} ; res.Xpred_max{ianimal,iseries,iprobe}(tidx)];
                                            popres.s_X{ianimal,iseries,iprobe,g} = [popres.s_X{ianimal,iseries,iprobe,g} ; res.X{ianimal,iseries}(tidx)];%res.traj{ianimal,iseries}(tidx)];%
                                            popres.s_PostXAll{ianimal,iseries,iprobe,g} = popres.s_PostXAll{ianimal,iseries,iprobe,g} + res.PostXSum{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)};
                                            popres.s_Xsum{ianimal,iseries,iprobe,g} = popres.s_Xsum{ianimal,iseries,iprobe,g} + res.XSum{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)};
                                            popres.s_trialID{ianimal,iseries,iprobe,g} = [popres.s_trialID{ianimal,iseries,iprobe,g} ; res.trialID{ianimal,iseries}(tidx)];
                                            
                                            if sum(ismember(res.outcomeVal{ianimal,iseries}, [3]))>0
                                                earlytrials = res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),ismember(res.outcomeVal{ianimal,iseries}, [3])};
                                            else
                                                earlytrials = false(size(tidx));
                                            end
                                            if sum(ismember(res.outcomeVal{ianimal,iseries}, [4]))>0
                                                latetrials = res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),ismember(res.outcomeVal{ianimal,iseries}, [4])};
                                            else
                                                latetrials = false(size(tidx));
                                            end
                                            tidx_incorrect = (earlytrials | latetrials) & ismember(res.Spdbin2{ianimal,iseries},spdbins);
                                            popres.s_correcttrials{ianimal,iseries,iprobe,g} = [popres.s_correcttrials{ianimal,iseries,iprobe,g} ; unique(res.trialID{ianimal,iseries}(tidx))];
                                            popres.s_incorrecttrials{ianimal,iseries,iprobe,g} = [popres.s_incorrecttrials{ianimal,iseries,iprobe,g} ; unique(res.trialID{ianimal,iseries}(tidx_incorrect))];
                                        end
                                    end
                                end
                            end
                        end
                        
                        if numel(unique([popres.s_correcttrials{ianimal,iseries,iprobe,g};popres.s_incorrecttrials{ianimal,iseries,iprobe,g}])) > 0
                            perf = numel(popres.s_correcttrials{ianimal,iseries,iprobe,g})...
                                /numel(([popres.s_correcttrials{ianimal,iseries,iprobe,g};popres.s_incorrecttrials{ianimal,iseries,iprobe,g}]));
                            %                         perf = sum(~ismember(popres.s_correcttrials{ianimal,iseries,iprobe,g},popres.s_incorrecttrials{ianimal,iseries,iprobe,g}))...
                            %                                                                          /numel(unique([popres.s_correcttrials{ianimal,iseries,iprobe,g};popres.s_incorrecttrials{ianimal,iseries,iprobe,g}]));
                            popres.performance{ianimal,iprobe,g} = [popres.performance{ianimal,iprobe,g} perf];
                        end
                        
                        dxlick = 1;
                        Xsmthlick = 100/4;
                        maxtol_licks = 1;
                        amp_th_licks = 0.5;
                        [lickmap,x,~] = fast1Dmap(popres.s_X{ianimal,iseries,iprobe,g},popres.s_licks{ianimal,iseries,iprobe,g},dxlick,1,Xsmthlick,1);
                        if ~isempty(lickmap)
                            lickmap = circshift(lickmap/sum(lickmap)*numel(x)/dxlick,round(numel(x)/dxlick/2));
                            if ianimal == 5 && iseries == 7 && g == 3
                                lickmap(1:30) = 0;
                            end
                            [~,imax] = max(lickmap);
                            lickpos = getCircularAverage(lickmap,amp_th_licks,maxtol_licks);
                            if ~isnan(lickpos)
                                lickpos = x(floor(lickpos)) + rem(lickpos,1)*(x(2) - x(1));
                            end
                            popres.s_meanlickX{ianimal,iseries,iprobe,g} = lickpos;
                        else
                            popres.s_meanlickX{ianimal,iseries,iprobe,g} = NaN;
                        end
%                         if g == 2 && popres.s_meanlickX{ianimal,iseries,iprobe,g}<popres.s_meanlickX{ianimal,iseries,iprobe,1}
%                             keyboard
%                         end
                        
                        s_xpredave = popres.s_Xpred_ave{ianimal,iseries,iprobe,g};
%                         trialIDs = unique(popres.s_trialID{ianimal,iseries,iprobe,g});
%                         Fave = NaN(Xrange,Xrange,numel(trialIDs));
%                         for tt = 1:numel(trialIDs)
%                             [Fave(:,:,tt), ~, ~, ~] = smoothhist2D_corrected([popres.s_X{ianimal,iseries,iprobe,g}(popres.s_trialID{ianimal,iseries,iprobe,g}==trialIDs(tt) & ~isnan(s_xpredave)) s_xpredave(popres.s_trialID{ianimal,iseries,iprobe,g}==trialIDs(tt) & ~isnan(s_xpredave))], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
%                         end
%                         Fave = nanmean(Fave,3);
                        
                        [Fave, ~, ~, ~] = smoothhist2D_corrected([popres.s_X{ianimal,iseries,iprobe,g}(~isnan(s_xpredave)) s_xpredave(~isnan(s_xpredave))], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
                        [xorig,yorig] = meshgrid(1:3*Xrange);
                        [xinterp,yinterp] = meshgrid(popres.dx:popres.dx:3*Xrange);
                        Fave_interp = interp2(xorig,yorig,repmat(Fave,[3 3]),xinterp,yinterp);
                        Fave_interp = Fave_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
                        Fave = Fave_interp;
                        
                        popres.s_DistriXaveAll{ianimal,iseries,iprobe,g} = Fave;  
                        popres.s_DistriXavepred{ianimal,iseries,iprobe,g} = getCircularAverage(popres.s_DistriXaveAll{ianimal,iseries,iprobe,g},amp_th_meandecmax,maxtol_meandecmax);
                        
                        xerr = s_xpredave - popres.s_X{ianimal,iseries,iprobe,g};
                        xerr(xerr > floor(Xrange/2)) = xerr(xerr > floor(Xrange/2)) - Xrange;
                        xerr(xerr < -floor(Xrange/2)) = xerr(xerr < -floor(Xrange/2)) + Xrange;
                        xerr = xerr + floor(Xrange/2);
                        [Fave, ~, ~, ~] = smoothhist2D_corrected([popres.s_X{ianimal,iseries,iprobe,g}(~isnan(s_xpredave)) xerr(~isnan(s_xpredave))], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
                        [xorig,yorig] = meshgrid(1:3*Xrange);
                        [xinterp,yinterp] = meshgrid(popres.dx:popres.dx:3*Xrange);
                        Fave_interp = interp2(xorig,yorig,repmat(Fave,[3 3]),xinterp,yinterp);
                        Fave_interp = Fave_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
                        Fave = Fave_interp;
                        
                        popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g} = Fave;
                        popres.s_DistriaveErrpred{ianimal,iseries,iprobe,g} = getCircularAverage(popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g},amp_th_meandecave,maxtol_meandecave);
                        
%                         for i = 1:size(popres.s_DistriXaveAll{ianimal,iseries,iprobe,g},2)
%                             popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g}(:,i) = circshift(popres.s_DistriXaveAll{ianimal,iseries,iprobe,g}(:,i),-i + round(size(popres.s_DistriXaveAll{ianimal,iseries,iprobe,g},1)/2));
%                         end
                        
                        %                     popres.s_meanErrX_ave{ianimal,iseries,iprobe,g} = getCircularAverage(nanmean(popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g},2),amp_th,maxtol);
                        x_interp = popres.dx:popres.dx:Xrange;
                        MeanXerr{g} = nanmean(popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g},2);
                        if sum(isnan(MeanXerr{g})) <= 2 && sum(MeanXerr{g}>amp_th_decave) > 0
                            [~, imax] = max(MeanXerr{g});
                            popres.s_meanErrX_ave{ianimal,iseries,iprobe,g} =  x_interp(min(numel(x_interp),round(getCircularAverage(MeanXerr{g},amp_th_decave,maxtol_decave))));%nanmean(getCircularAverage(popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g},amp_th_decave,maxtol_decave));%
                        else
                            popres.s_meanErrX_ave{ianimal,iseries,iprobe,g} = NaN;
                        end
                        
                        for xx = 1:5
%                             if xx == 1
%                                 xseg = [95:100 1:5 40:50 60:70];
%                             else
%                                 xseg = [15:25 45:55 75:85];
%                             end
%                             MeanXerr{g} =
%                             nanmean(popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g}(:,xseg),2);
                        lickidx = find(popres.s_licks{ianimal,iseries,iprobe,g});
                        timelicksidx = true(size(popres.s_licks{ianimal,iseries,iprobe,g}));
                        %to comment out if one wants to remove time points
                        %after the first licks
                        for l = 1:numel(lickidx)
                            timelicksidx(lickidx(l):lickidx(l)+find(popres.s_X{ianimal,iseries,iprobe,g}(lickidx(l)+1:end)>10,1,'first')) = false;
                        end
                        s_xpredave = popres.s_Xpred_ave{ianimal,iseries,iprobe,g};
                        [Fave, ~, ~, ~] = smoothhist2D_corrected([popres.s_X{ianimal,iseries,iprobe,g}(~isnan(s_xpredave) & timelicksidx) s_xpredave(~isnan(s_xpredave) & timelicksidx)], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
                        
                        [xorig,yorig] = meshgrid(1:3*Xrange);
                        [xinterp,yinterp] = meshgrid(popres.dx:popres.dx:3*Xrange);
                        Fave_interp = interp2(xorig,yorig,repmat(Fave,[3 3]),xinterp,yinterp);
                        Fave_interp = Fave_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
                        Fave = Fave_interp;
                        
                        for i = 1:size(popres.s_DistriXaveAll{ianimal,iseries,iprobe,g},2)
                            Fave(:,i) = circshift(Fave(:,i),-i + floor(size(popres.s_DistriXaveAll{ianimal,iseries,iprobe,g},2)/2));
                        end
                        ampth_lickErr = amp_th_decave;
                        x_interp = popres.dx:popres.dx:Xrange;
%                         MeanXerr{g} = nanmean(popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g}(:,((xx-1)*20+1):end),2);
                        MeanXerr{g} = nanmean(Fave(:,((xx-1)*20+1):end),2);
                        if sum(isnan(MeanXerr{g})) <= 2 && sum(MeanXerr{g}>ampth_lickErr) > 0
                            [~, imax] = max(MeanXerr{g});
                            try
                            popres.s_meanErrX_aveSeg{ianimal,iseries,iprobe,g}(xx) =  x_interp(min(numel(x_interp),round(getCircularAverage(MeanXerr{g},ampth_lickErr,maxtol_decave))));%nanmean(getCircularAverage(popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g},amp_th_decave,maxtol_decave));%
                            catch
                                keyboard
                            end
                        else
                            popres.s_meanErrX_aveSeg{ianimal,iseries,iprobe,g}(xx) = NaN;
                        end
                        end
                        
                        if FshowsessionFig
                            figure(f_distriave)
                            subplot(4,nseries,iseries_valid)
                            hold on;
                            plot(MeanXerr{g},cl{g});
                            hold on;plot([popres.s_meanErrX_ave{ianimal,iseries,iprobe,g} popres.s_meanErrX_ave{ianimal,iseries,iprobe,g}],[0 3],cl{g})
                            set(gca,'PlotBoxAspectRatio', [1 1 1]);
                            title(num2str(expt(ianimal).series{iseries}));
                            subplot(4,nseries,g*nseries+iseries_valid)
                            imagesc(popres.s_DistriXaveAll{ianimal,iseries,iprobe,g})
                            hold on;plot([0 Xrange],[0 Xrange],'k');
                            colormap(parula);%colormap(RedWhiteBlue);%
                            set(gca,'Clim',[0 2],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
                        end
                        
                        s_xpredmax = popres.s_Xpred_max{ianimal,iseries,iprobe,g};
                        [Fmax, ~, ~, ~] = smoothhist2D_corrected([popres.s_X{ianimal,iseries,iprobe,g}(~isnan(s_xpredmax)) s_xpredmax(~isnan(s_xpredmax))], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
                        
                        [xorig,yorig] = meshgrid(1:3*Xrange);
                        [xinterp,yinterp] = meshgrid(popres.dx:popres.dx:3*Xrange);
                        Fmax_interp = interp2(xorig,yorig,repmat(Fmax,[3 3]),xinterp,yinterp);
                        Fmax_interp = Fmax_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
                        Fmax = Fmax_interp;
                        
                        popres.s_DistriXmaxAll{ianimal,iseries,iprobe,g} = Fmax;
                        popres.s_DistriXmaxpred{ianimal,iseries,iprobe,g} = getCircularAverage(popres.s_DistriXmaxAll{ianimal,iseries,iprobe,g},amp_th_meandecmax,maxtol_meandecmax);
                        
                        xerr = s_xpredmax - popres.s_X{ianimal,iseries,iprobe,g};
                        xerr(xerr > floor(Xrange/2)) = xerr(xerr > floor(Xrange/2)) - Xrange;
                        xerr(xerr < -floor(Xrange/2)) = xerr(xerr < -floor(Xrange/2)) + Xrange;
                        xerr = xerr + floor(Xrange/2);
                        [Fave, ~, ~, ~] = smoothhist2D_corrected([popres.s_X{ianimal,iseries,iprobe,g}(~isnan(s_xpredmax)) xerr(~isnan(s_xpredmax))], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
                        [xorig,yorig] = meshgrid(1:3*Xrange);
                        [xinterp,yinterp] = meshgrid(popres.dx:popres.dx:3*Xrange);
                        Fave_interp = interp2(xorig,yorig,repmat(Fave,[3 3]),xinterp,yinterp);
                        Fave_interp = Fave_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
                        Fave = Fave_interp;
                        
                        popres.s_DistrimaxErrAll{ianimal,iseries,iprobe,g} = Fave;
                        popres.s_DistrimaxErrpred{ianimal,iseries,iprobe,g} = getCircularAverage(popres.s_DistrimaxErrAll{ianimal,iseries,iprobe,g},amp_th_meandecave,maxtol_meandecave);
%                         for i = 1:size(popres.s_DistriXmaxAll{ianimal,iseries,iprobe,g},2)
%                             popres.s_DistrimaxErrAll{ianimal,iseries,iprobe,g}(:,i) = circshift(popres.s_DistriXmaxAll{ianimal,iseries,iprobe,g}(:,i),-i + floor(size(popres.s_DistriXmaxAll{ianimal,iseries,iprobe,g},2)/2));
%                         end
                        %                     popres.s_meanErrX_max{ianimal,iseries,iprobe,g} = getCircularAverage(nanmean(popres.s_DistrimaxErrAll{ianimal,iseries,iprobe,g},2),amp_th,maxtol);
                        x_interp = popres.dx:popres.dx:Xrange;
                        MeanXerr{g} = nanmean(popres.s_DistrimaxErrAll{ianimal,iseries,iprobe,g},2);
                        if sum(isnan(MeanXerr{g})) <= 2 && sum(MeanXerr{g}>amp_th_decmax) > 0
                            [~, imax] = max(MeanXerr{g});
                            popres.s_meanErrX_max{ianimal,iseries,iprobe,g} =  x_interp(min(numel(x_interp),round(getCircularAverage(MeanXerr{g},amp_th_decmax,maxtol_decmax))));%nanmean(getCircularAverage(popres.s_DistrimaxErrAll{ianimal,iseries,iprobe,g},amp_th_decmax,maxtol_decmax));%
                        else
                            popres.s_meanErrX_max{ianimal,iseries,iprobe,g} = NaN;
                        end
                        if FshowsessionFig
                            figure(f_distrimax)
                            subplot(4,nseries,iseries_valid)
                            hold on;
                            plot(MeanXerr{g},cl{g});
                            hold on;plot([popres.s_meanErrX_max{ianimal,iseries,iprobe,g} popres.s_meanErrX_max{ianimal,iseries,iprobe,g}],[0 3],cl{g})
                            set(gca,'PlotBoxAspectRatio', [1 1 1]);
                            title(num2str(expt(ianimal).series{iseries}));
                            subplot(4,nseries,g*nseries+iseries_valid)
                            imagesc(popres.s_DistriXmaxAll{ianimal,iseries,iprobe,g})
                            hold on;plot([0 Xrange],[0 Xrange],'k');
                            colormap(parula);%colormap(RedWhiteBlue);%
                            set(gca,'Clim',[0 2],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
                        end
                        
                        if sum(popres.s_PostXAll{ianimal,iseries,iprobe,g}(:)) == 0
                            popres.s_PostXAll{ianimal,iseries,iprobe,g} = NaN(Xrange);
                        end
                        popres.s_PostXAll{ianimal,iseries,iprobe,g} = smooth2D(popres.s_PostXAll{ianimal,iseries,iprobe,g}./repmat(popres.s_Xsum{ianimal,iseries,iprobe,g},[Xrange 1]), lambdaSmooth);
                        [xorig,yorig] = meshgrid(1:3*Xrange);
                        [xinterp,yinterp] = meshgrid(popres.dx:popres.dx:3*Xrange);
                        Post_interp = interp2(xorig,yorig,repmat(popres.s_PostXAll{ianimal,iseries,iprobe,g},[3 3]),xinterp,yinterp);
                        
                        Post_interp = Post_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
                        popres.s_PostXAll{ianimal,iseries,iprobe,g} = Post_interp;
                        
                        for i = 1:size(popres.s_PostXAll{ianimal,iseries,iprobe,g},2)
                            popres.s_PostErrAll{ianimal,iseries,iprobe,g}(:,i) = circshift(popres.s_PostXAll{ianimal,iseries,iprobe,g}(:,i),-i + floor(size(popres.s_PostXAll{ianimal,iseries,iprobe,g},2)/2));
                        end
                        %                     popres.s_meanErrX_Post{ianimal,iseries,iprobe,g} = getCircularAverage(nanmean(popres.s_PostErrAll{ianimal,iseries,iprobe,g},2),amp_th,maxtol);
                        x_interp = popres.dx:popres.dx:Xrange;
                        MeanXerr{g} = nanmean(popres.s_PostErrAll{ianimal,iseries,iprobe,g},2);
                        if sum(isnan(MeanXerr{g})) <= 2 && sum(MeanXerr{g}>amp_th_post) > 0
                            [~, imax] = max(MeanXerr{g});
                            popres.s_meanErrX_Post{ianimal,iseries,iprobe,g} = x_interp(min(numel(x_interp),round(getCircularAverage(MeanXerr{g},amp_th_post,maxtol_post))));%nanmean(getCircularAverage(popres.s_PostErrAll{ianimal,iseries,iprobe,g},amp_th_post,maxtol_post));%
                        else
                            popres.s_meanErrX_Post{ianimal,iseries,iprobe,g} = NaN;
                        end
                        if FshowsessionFig
                            figure(f_post)
                            subplot(4,nseries,iseries_valid)
                            hold on;
                            plot(MeanXerr{g},cl{g});
                            hold on;plot([popres.s_meanErrX_Post{ianimal,iseries,iprobe,g}  popres.s_meanErrX_Post{ianimal,iseries,iprobe,g}],[0 3],cl{g})
                            set(gca,'PlotBoxAspectRatio', [1 1 1]);
                            title(num2str(expt(ianimal).series{iseries}));
                            subplot(4,nseries,g*nseries+iseries_valid)
                            imagesc(popres.s_PostXAll{ianimal,iseries,iprobe,g})
                            hold on;plot([0 Xrange],[0 Xrange],'k');
                            colormap(parula);%colormap(RedWhiteBlue);%
                            set(gca,'Clim',[0 2],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
                        end
                    end
                    %                 pause
                    %                 delete(f)
                end
            end
        end
        if FshowsessionFig
            figure(f_distrimax)
            savefig2pdf([figdirname filesep popres.probestr{iprobe} '-' expt(ianimal).animal 'distri-Max' '.pdf']);
            savefig2pdf([figdirname filesep popres.probestr{iprobe} '-' expt(ianimal).animal 'distri-Max']);
            figure(f_distriave)
            savefig2pdf([figdirname filesep popres.probestr{iprobe} '-' expt(ianimal).animal 'distri-Ave' '.pdf']);
            savefig2pdf([figdirname filesep popres.probestr{iprobe} '-' expt(ianimal).animal 'distri-Ave']);
            figure(f_post)
            savefig2pdf([figdirname filesep popres.probestr{iprobe} '-' expt(ianimal).animal 'post' '.pdf']);
            savefig2pdf([figdirname filesep popres.probestr{iprobe} '-' expt(ianimal).animal 'post']);
            clf(f_post);
            clf(f_distrimax);
            clf(f_distriave);
        end
    end
    
    
    %to average across sessions instead of across trials
%     iseries_ave = zeros(1,3);
%     for ianimal = 1:nanimal
%         for iseries = 1:numel(expt(ianimal).series)
%             if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
%                 for g = 1:3
%                     if ~isempty(popres.s_DistriXaveAll{ianimal,iseries,iprobe,g}) && sum(isnan(popres.s_DistriXaveAll{ianimal,iseries,iprobe,g}(:))) == 0
%                         iseries_ave(g) = iseries_ave(g) + 1;
%                         if isempty(popres.DistriXaveAll{iprobe,g})
%                             popres.DistriXaveAll{iprobe,g} = popres.s_DistriXaveAll{ianimal,iseries,iprobe,g};
%                             popres.DistriXmaxAll{iprobe,g} = popres.s_DistriXmaxAll{ianimal,iseries,iprobe,g};
%                             popres.PostXAll{iprobe,g} = popres.s_PostXAll{ianimal,iseries,iprobe,g};
%                             
%                             popres.DistriaveErrAll{iprobe,g} = popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g};
%                             popres.DistrimaxErrAll{iprobe,g} = popres.s_DistrimaxErrAll{ianimal,iseries,iprobe,g};
%                             popres.PostErrAll{iprobe,g} = popres.s_PostErrAll{ianimal,iseries,iprobe,g};
%                         else
%                             popres.DistriXaveAll{iprobe,g} = popres.DistriXaveAll{iprobe,g} + popres.s_DistriXaveAll{ianimal,iseries,iprobe,g};
%                             popres.DistriXmaxAll{iprobe,g} = popres.DistriXmaxAll{iprobe,g} + popres.s_DistriXmaxAll{ianimal,iseries,iprobe,g};
%                             popres.PostXAll{iprobe,g} = popres.PostXAll{iprobe,g} + popres.s_PostXAll{ianimal,iseries,iprobe,g};
%                             
%                             popres.DistriaveErrAll{iprobe,g} = popres.DistriaveErrAll{iprobe,g} + popres.s_DistriaveErrAll{ianimal,iseries,iprobe,g};
%                             popres.DistrimaxErrAll{iprobe,g} = popres.DistrimaxErrAll{iprobe,g} + popres.s_DistrimaxErrAll{ianimal,iseries,iprobe,g};
%                             popres.PostErrAll{iprobe,g} = popres.PostErrAll{iprobe,g} + popres.s_PostErrAll{ianimal,iseries,iprobe,g};
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     for g = [2 1 3]
%         popres.DistriXaveAll{iprobe,g} = popres.DistriXaveAll{iprobe,g}/iseries_ave(g);
%         popres.DistriXmaxAll{iprobe,g} = popres.DistriXmaxAll{iprobe,g}/iseries_ave(g);
%         popres.PostXAll{iprobe,g} = popres.PostXAll{iprobe,g}/iseries_ave(g);
%         
%         popres.DistriaveErrAll{iprobe,g} = popres.DistriaveErrAll{iprobe,g}/iseries_ave(g);
%         popres.DistrimaxErrAll{iprobe,g} = popres.DistrimaxErrAll{iprobe,g}/iseries_ave(g);
%         popres.PostErrAll{iprobe,g} = popres.PostErrAll{iprobe,g}/iseries_ave(g);
%         
%         popres.DistriXavepred{iprobe,g} = getCircularAverage(popres.DistriXaveAll{iprobe,g},amp_th_meandecave,maxtol_meandecave);
%         popres.DistriaveErrpred{iprobe,g} = getCircularAverage(popres.DistriaveErrAll{iprobe,g},amp_th_meandecave,maxtol_meandecave);
%         
%         popres.DistriXmaxpred{iprobe,g} = getCircularAverage(popres.DistriXmaxAll{iprobe,g},amp_th_meandecmax,maxtol_meandecmax);
%         popres.DistrimaxErrpred{iprobe,g} = getCircularAverage(popres.DistrimaxErrAll{iprobe,g},amp_th_meandecmax,maxtol_meandecmax);
%         
%         popres.PostXpred{iprobe,g} = getCircularAverage(popres.PostXAll{iprobe,g},amp_th_meanpost,maxtol_meanpost);
%         popres.PostErrpred{iprobe,g} = getCircularAverage(popres.PostErrAll{iprobe,g},amp_th_meanpost,maxtol_meanpost);
%         
%         [popres.PostErrpred_reg{iprobe,g},popres.xreg_interp] = Xpredcorrection(popres.PostXpred{iprobe,g},popres.PostXpred{iprobe,2},size(popres.PostXpred{iprobe,2},1),popres.dx);
%         [popres.DistrimaxErrpred_reg{iprobe,g},popres.xreg_interp] = Xpredcorrection(popres.DistriXmaxpred{iprobe,g},popres.DistriXmaxpred{iprobe,2},size(popres.DistriXmaxpred{iprobe,2},1),popres.dx);
%         [popres.DistriaveErrpred_reg{iprobe,g},popres.xreg_interp] = Xpredcorrection(popres.DistriXavepred{iprobe,g},popres.DistriXavepred{iprobe,2},size(popres.DistriXavepred{iprobe,2},1),popres.dx);
%     end
    
    kfold = 3;%30;
    for g = [2 1 3]
        [xorig,yorig] = meshgrid(1:3*Xrange);
        [xinterp,yinterp] = meshgrid(popres.dx:popres.dx:3*Xrange);
        
        xpredave = (popres.Xpred_ave{iprobe,g});
        [Fave, ~, ~, ~] = smoothhist2D_corrected([popres.X{iprobe,g}(~isnan(xpredave)) xpredave(~isnan(xpredave))], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
        Fave_interp = interp2(xorig,yorig,repmat(Fave,[3 3]),xinterp,yinterp);
        Fave_interp = Fave_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
        Fave = Fave_interp;
        popres.DistriXaveAll{iprobe,g} = Fave;
        
        xpredmax = (popres.Xpred_max{iprobe,g});
        [Fmax, ~, ~, ~] = smoothhist2D_corrected([popres.X{iprobe,g}(~isnan(xpredmax)) xpredmax(~isnan(xpredmax))], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
        Fmax_interp = interp2(xorig,yorig,repmat(Fmax,[3 3]),xinterp,yinterp);
        Fmax_interp = Fmax_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
        Fmax = Fmax_interp;
        popres.DistriXmaxAll{iprobe,g} = Fmax;
        
        for i = 1:size(popres.PostXAll{iprobe,g},2)
            popres.PostErrAll{iprobe,g}(:,i) = circshift(popres.PostXAll{iprobe,g}(:,i),-i + floor(size(popres.PostXAll{iprobe,g},2)/2));
        end
        popres.PostErrAll{iprobe,g} = smooth2D(popres.PostErrAll{iprobe,g}./repmat(popres.Xsum{iprobe,g},[Xrange 1]), lambdaSmooth);
        Post_interp = interp2(xorig,yorig,repmat(popres.PostErrAll{iprobe,g},[3 3]),xinterp,yinterp);
        Post_interp = Post_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
        popres.PostErrAll{iprobe,g} = Post_interp;
        
        popres.PostXAll{iprobe,g} = smooth2D(popres.PostXAll{iprobe,g}./repmat(popres.Xsum{iprobe,g},[Xrange 1]), lambdaSmooth);
        Post_interp = interp2(xorig,yorig,repmat(popres.PostXAll{iprobe,g},[3 3]),xinterp,yinterp);
        Post_interp = Post_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
        popres.PostXAll{iprobe,g} = Post_interp;
        
        
        popres.PostXpred{iprobe,g} = getCircularAverage(popres.PostXAll{iprobe,g},amp_th_meanpost,maxtol_meanpost);
        popres.PostErrpred{iprobe,g} = getCircularAverage(popres.PostErrAll{iprobe,g},amp_th_meanpost,maxtol_meanpost);
        
%         for i = 1:size(popres.PostXAll{iprobe,g},2)
%             popres.PostErrAll{iprobe,g}(:,i) = circshift(popres.PostXAll{iprobe,g}(:,i),-i + floor(size(popres.PostXAll{iprobe,g},2)/2));
%         end
%         popres.PostErrpred{iprobe,g} = popres.PostXpred{iprobe,g} - (1:size(popres.PostXAll{iprobe,g},2))';
%         popres.PostErrpred{iprobe,g}(popres.PostErrpred{iprobe,g}>floor(size(popres.PostXAll{iprobe,g},2)/2)) = popres.PostErrpred{iprobe,g}(popres.PostErrpred{iprobe,g}>floor(size(popres.PostXAll{iprobe,g},2)/2))-size(popres.PostXAll{iprobe,g},2);
%         popres.PostErrpred{iprobe,g}(popres.PostErrpred{iprobe,g}<-floor(size(popres.PostXAll{iprobe,g},2)/2)) = popres.PostErrpred{iprobe,g}(popres.PostErrpred{iprobe,g}<-floor(size(popres.PostXAll{iprobe,g},2)/2)) + size(popres.PostXAll{iprobe,g},2);
        
        
        xpredmax = (popres.Xpred_max{iprobe,g});
        xerr = xpredmax - popres.X{iprobe,g};
        xerr(xerr > floor(Xrange/2)) = xerr(xerr > floor(Xrange/2)) - Xrange;
        xerr(xerr < -floor(Xrange/2)) = xerr(xerr < -floor(Xrange/2)) + Xrange;
        xerr = xerr + floor(Xrange/2);
        [Fave, ~, ~, ~] = smoothhist2D_corrected([popres.X{iprobe,g}(~isnan(xpredave)) xerr(~isnan(xpredave))], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
        Fave_interp = interp2(xorig,yorig,repmat(Fave,[3 3]),xinterp,yinterp);
        Fave_interp = Fave_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
        Fave = Fave_interp;
        popres.DistrimaxErrAll{iprobe,g} = Fave;
        popres.DistrimaxErrpred{iprobe,g} = getCircularAverage(popres.DistrimaxErrAll{iprobe,g},amp_th_meandecmax,maxtol_meandecmax);
        
        popres.DistriXmaxpred{iprobe,g} = getCircularAverage(popres.DistriXmaxAll{iprobe,g},amp_th_meandecmax,maxtol_meandecmax);
%         for i = 1:size(popres.DistriXmaxAll{iprobe,g},2)
%             popres.DistrimaxErrAll{iprobe,g}(:,i) = circshift(popres.DistriXmaxAll{iprobe,g}(:,i),-i + floor(size(popres.DistriXmaxAll{iprobe,g},2)/2));
%         end
%         popres.DistrimaxErrpred{iprobe,g} = popres.DistriXmaxpred{iprobe,g} - (1:size(popres.DistriXmaxAll{iprobe,g},2))';
%         popres.DistrimaxErrpred{iprobe,g}(popres.DistrimaxErrpred{iprobe,g}>floor(size(popres.DistriXmaxAll{iprobe,g},2)/2)) = popres.DistrimaxErrpred{iprobe,g}(popres.DistrimaxErrpred{iprobe,g}>floor(size(popres.DistriXmaxAll{iprobe,g},2)/2))-size(popres.DistriXmaxAll{iprobe,g},2);
%         popres.DistrimaxErrpred{iprobe,g}(popres.DistrimaxErrpred{iprobe,g}<-floor(size(popres.DistriXmaxAll{iprobe,g},2)/2)) = popres.DistrimaxErrpred{iprobe,g}(popres.DistrimaxErrpred{iprobe,g}<-floor(size(popres.DistriXmaxAll{iprobe,g},2)/2)) + size(popres.DistriXmaxAll{iprobe,g},2);
        
        xpredave = (popres.Xpred_ave{iprobe,g});
        xerr = xpredave - popres.X{iprobe,g};
        xerr(xerr > floor(Xrange/2)) = xerr(xerr > floor(Xrange/2)) - Xrange;
        xerr(xerr < -floor(Xrange/2)) = xerr(xerr < -floor(Xrange/2)) + Xrange;
        xerr = xerr + floor(Xrange/2);
        [Fave, ~, ~, ~] = smoothhist2D_corrected([popres.X{iprobe,g}(~isnan(xpredave)) xerr(~isnan(xpredave))], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
        Fave_interp = interp2(xorig,yorig,repmat(Fave,[3 3]),xinterp,yinterp);
        Fave_interp = Fave_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
        Fave = Fave_interp;
        popres.DistriaveErrAll{iprobe,g} = Fave;
        popres.DistriaveErrpred{iprobe,g} = getCircularAverage(popres.DistriaveErrAll{iprobe,g},amp_th_meandecave,maxtol_meandecave);
%         diffX = [smooth(diff(repmat(popres.DistriaveErrpred{iprobe,g},[3 1])),10/popres.dx);0]/popres.dx;
%         popres.DistriaveErrpred_diff{iprobe,g} = diffX(round(numel(diffX)/3)+1:round(2*numel(diffX)/3));
        popres.DistriaveErrpred_rel{iprobe,g} = popres.DistriaveErrpred{iprobe,g} - popres.DistriaveErrpred{iprobe,2};
%         popres.DistriaveErrpred_diffrel{iprobe,g} = popres.DistriaveErrpred_diff{iprobe,g} - popres.DistriaveErrpred_diff{iprobe,2};
        
        popres.DistriXavepred{iprobe,g} = getCircularAverage(popres.DistriXaveAll{iprobe,g},amp_th_meandecave,maxtol_meandecave);
        
        [~,isort] = sort(popres.DistriXavepred{iprobe,g},'ascend');
        x_i = [popres.DistriXavepred{iprobe,g}(isort);popres.DistriXavepred{iprobe,g}(isort(end))+popres.DistriXavepred{iprobe,g}(isort);2*popres.DistriXavepred{iprobe,g}(isort(end))+popres.DistriXavepred{iprobe,g}(isort)];
        x_o = 0:(3*max(popres.DistriXavepred{iprobe,2})-1);
        err_i = repmat(popres.DistriaveErrpred{iprobe,g}(isort),[3 1]);
        err = interp1(x_i,err_i,x_o);
        err = err(floor(numel(err)/3)+1:floor(2*numel(err)/3));
        Xvis = (1:numel(popres.DistriXavepred{iprobe,g}))';
        Xvis_sort = Xvis(isort);
        Xvis_sort(Xvis_sort>Xvis + floor(numel(popres.DistriXavepred{iprobe,g})/2)) = Xvis_sort(Xvis_sort>Xvis + floor(numel(popres.DistriXavepred{iprobe,g})/2)) - numel(popres.DistriXavepred{iprobe,g});
        Xvis_sort(Xvis_sort<Xvis - floor(numel(popres.DistriXavepred{iprobe,g})/2)) = Xvis_sort(Xvis_sort<Xvis - floor(numel(popres.DistriXavepred{iprobe,g})/2)) + numel(popres.DistriXavepred{iprobe,g});
        Xvis_sort = [Xvis_sort;Xvis_sort(end) + Xvis_sort; 2*Xvis_sort(end) + Xvis_sort];
        Xvis_i = interp1(x_i,Xvis_sort,x_o);
        Xvis_i = Xvis_i(floor(numel(Xvis_i)/3)+1:floor(2*numel(Xvis_i)/3)) - Xvis_i(floor(numel(Xvis_i)/3));
        Xvis_i(Xvis_i < 0) = Xvis_i(Xvis_i < 0) + numel(popres.DistriXavepred{iprobe,g});
        Xvis_i(Xvis_i > numel(popres.DistriXavepred{iprobe,g})) = Xvis_i(Xvis_i > numel(popres.DistriXavepred{iprobe,g})) - numel(popres.DistriXavepred{iprobe,g});
        popres.DistriaveErrpred_dec{iprobe,g} = err(:);
        
        [~,isort] = sort(popres.DistriXavepred{iprobe,2},'ascend');
        x_i = [popres.DistriXavepred{iprobe,2}(isort);popres.DistriXavepred{iprobe,2}(isort(end))+popres.DistriXavepred{iprobe,2}(isort);2*popres.DistriXavepred{iprobe,2}(isort(end))+popres.DistriXavepred{iprobe,2}(isort)];
        x_o = 0:(3*max(popres.DistriXavepred{iprobe,2})-1);
        err_i = repmat(popres.DistriaveErrpred{iprobe,2}(isort),[3 1]);
        err_ref = interp1(x_i,err_i,x_o);
        err_ref = err_ref(floor(numel(err_ref)/3)+1:floor(2*numel(err_ref)/3));
        [~,isortV] = unique(Xvis_i,'sorted');
        err_diff = interp1(Xvis_i(isortV),err(isortV) - err_ref(isortV),0:(Xrange/popres.dx - 1));
        err_diff = smooth(err_diff,round(1/popres.dx));
        popres.DistriaveErrpred_reldec{iprobe,g} = err_diff(:);%err(:) - err_ref(:);
        
        nx = numel(popres.DistriaveErrpred{iprobe,g},1);
        popres.DistriaveErrpreddiffmat{iprobe,g} = repmat(popres.DistriaveErrpred{iprobe,g},[1 nx]) - repmat(popres.DistriaveErrpred{iprobe,g}',[nx 1]);
        popres.DistriaveErrpreddiffmat_rel{iprobe,g} = repmat(popres.DistriaveErrpred_rel{iprobe,g},[1 nx]) - repmat(popres.DistriaveErrpred_rel{iprobe,g}',[nx 1]);
        popres.DistriaveErrpreddiffmat_dec{iprobe,g} = repmat(popres.DistriaveErrpred_dec{iprobe,g},[1 nx]) - repmat(popres.DistriaveErrpred_dec{iprobe,g}',[nx 1]);
        popres.DistriaveErrpreddiffmat_reldec{iprobe,g} = repmat(popres.DistriaveErrpred_reldec{iprobe,g},[1 nx]) - repmat(popres.DistriaveErrpred_reldec{iprobe,g}',[nx 1]);
        
%         [~,isort] = sort(popres.DistriXavepred{iprobe,g},'ascend');
%         x_i = [popres.DistriXavepred{iprobe,g}(isort);popres.DistriXavepred{iprobe,g}(isort(end))+popres.DistriXavepred{iprobe,g}(isort);2*popres.DistriXavepred{iprobe,g}(isort(end))+popres.DistriXavepred{iprobe,g}(isort)];
%         x_o = 0:(3*max(popres.DistriXavepred{iprobe,2})-1);
%         err_i = repmat(popres.DistriaveErrpred_diff{iprobe,g}(isort),[3 1]);
%         err = interp1(x_i,err_i,x_o,'spline');
%         err = err(floor(numel(err)/3)+1:floor(2*numel(err)/3));
%         popres.DistriaveErrpred_diffdec{iprobe,g} = err(:);
%         
%         [~,isort] = sort(popres.DistriXavepred{iprobe,2},'ascend');
%         x_i = [popres.DistriXavepred{iprobe,2}(isort);popres.DistriXavepred{iprobe,2}(isort(end))+popres.DistriXavepred{iprobe,2}(isort);2*popres.DistriXavepred{iprobe,2}(isort(end))+popres.DistriXavepred{iprobe,2}(isort)];
%         x_o = 0:(3*max(popres.DistriXavepred{iprobe,2})-1);
%         err_i = repmat(popres.DistriaveErrpred_diff{iprobe,2}(isort),[3 1]);
%         err_ref = interp1(x_i,err_i,x_o,'spline');
%         err_ref = err_ref(floor(numel(err_ref)/3)+1:floor(2*numel(err_ref)/3));
%         popres.DistriaveErrpred_diffreldec{iprobe,g} = err(:) - err_ref(:);
        
%         for i = 1:size(popres.DistriXaveAll{iprobe,g},2)
%             popres.DistriaveErrAll{iprobe,g}(:,i) = circshift(popres.DistriXaveAll{iprobe,g}(:,i),-i + floor(size(popres.DistriXaveAll{iprobe,g},2)/2));
%         end
%         popres.DistriaveErrpred{iprobe,g} = popres.DistriXavepred{iprobe,g} - (1:size(popres.DistriXaveAll{iprobe,g},2))';
%         popres.DistriaveErrpred{iprobe,g}(popres.DistriaveErrpred{iprobe,g}>floor(size(popres.DistriXaveAll{iprobe,g},2)/2)) = popres.DistriaveErrpred{iprobe,g}(popres.DistriaveErrpred{iprobe,g}>floor(size(popres.DistriXaveAll{iprobe,g},2)/2))-size(popres.DistriXaveAll{iprobe,g},2);
%         popres.DistriaveErrpred{iprobe,g}(popres.DistriaveErrpred{iprobe,g}<-floor(size(popres.DistriXaveAll{iprobe,g},2)/2)) = popres.DistriaveErrpred{iprobe,g}(popres.DistriaveErrpred{iprobe,g}<-floor(size(popres.DistriXaveAll{iprobe,g},2)/2)) + size(popres.DistriXaveAll{iprobe,g},2);

        [CVO] = crossValPartition(1:numel(popres.Xpred_ave{iprobe,g}), kfold);
        popres.DistriXaveAll_std{iprobe,g} = 0;
        popres.DistriXavepred_std{iprobe,g} = 0;
        popres.DistriaveErrAll_std{iprobe,g} = 0;
        popres.DistriaveErrpred_std{iprobe,g} = 0;
        popres.DistriaveErrpred_relstd{iprobe,g} = 0;
        popres.DistriaveErrpred_reldecstd{iprobe,g} = 0;
        popres.DistriaveErrpred_decstd{iprobe,g} = 0;
        popres.DistriaveErrpred_diffstd{iprobe,g} = 0;
        popres.DistriaveErrpred_diffrelstd{iprobe,g} = 0;
        popres.DistriaveErrpred_diffdecstd{iprobe,g} = 0;
        popres.DistriaveErrpred_diffreldecstd{iprobe,g} = 0;
        popres.DistriaveErrpreddiffmat_std{iprobe,g} = 0;
        popres.DistriaveErrpreddiffmat_relstd{iprobe,g} = 0;
        popres.DistriaveErrpreddiffmat_decstd{iprobe,g} = 0;
        popres.DistriaveErrpreddiffmat_reldecstd{iprobe,g} = 0;
        for kiter = 1:CVO.kfold
            idxCVO = (CVO.train{kiter})';
            xpredave = (popres.Xpred_ave{iprobe,g});
            [Fave, ~, ~, ~] = smoothhist2D_corrected([popres.X{iprobe,g}(~isnan(xpredave) & idxCVO) xpredave(~isnan(xpredave) & idxCVO)], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
            Fave_interp = interp2(xorig,yorig,repmat(Fave,[3 3]),xinterp,yinterp);
            Fave_interp = Fave_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
            Fave = Fave_interp;
            popres.DistriXaveAll_std{iprobe,g} = popres.DistriXaveAll_std{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(Fave-popres.DistriXaveAll{iprobe,g}).^2;
            avepred_iter = getCircularAverage(Fave,amp_th_meandecave,maxtol_meandecave);
            popres.DistriXavepred_std{iprobe,g} = popres.DistriXavepred_std{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(avepred_iter - popres.DistriXavepred{iprobe,g}).^2;
            
            xpredave = (popres.Xpred_ave{iprobe,g});
            xerr = xpredave - popres.X{iprobe,g};
            xerr(xerr > floor(Xrange/2)) = xerr(xerr > floor(Xrange/2)) - Xrange;
            xerr(xerr < -floor(Xrange/2)) = xerr(xerr < -floor(Xrange/2)) + Xrange;
            xerr = xerr + floor(Xrange/2);
            [Fave, ~, ~, ~] = smoothhist2D_corrected([popres.X{iprobe,g}(~isnan(xpredave) & idxCVO) xerr(~isnan(xpredave) & idxCVO)], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
            Fave_interp = interp2(xorig,yorig,repmat(Fave,[3 3]),xinterp,yinterp);
            Fave_interp = Fave_interp((size(xinterp,1)/3+1):(2*(size(xinterp,1)/3)),(size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)));
            Fave = Fave_interp;
            popres.DistriaveErrAll_std{iprobe,g} = popres.DistriaveErrAll_std{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(Fave - popres.DistriaveErrAll{iprobe,g}).^2;
            aveerrpred_iter = getCircularAverage(Fave,amp_th_meandecave,maxtol_meandecave);
%             diffX = [smooth(diff(repmat(aveerrpred_iter,[3 1])),10/popres.dx);0]/popres.dx;
%             aveerrpred_diffiter = diffX(round(numel(diffX)/3)+1:round(2*numel(diffX)/3));
            aveerrpred_reliter = aveerrpred_iter - popres.DistriaveErrpred{iprobe,2};
%             aveerrpred_diffreliter = aveerrpred_diffiter - popres.DistriaveErrpred_diff{iprobe,2};
            
            popres.DistriaveErrpred_std{iprobe,g} = popres.DistriaveErrpred_std{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpred_iter - popres.DistriaveErrpred{iprobe,g}).^2;
            popres.DistriaveErrpred_relstd{iprobe,g} = popres.DistriaveErrpred_relstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpred_reliter - popres.DistriaveErrpred_rel{iprobe,g}).^2;
%             popres.DistriaveErrpred_diffstd{iprobe,g} = popres.DistriaveErrpred_diffstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpred_diffiter - popres.DistriaveErrpred_diff{iprobe,g}).^2;
%             popres.DistriaveErrpred_diffrelstd{iprobe,g} = popres.DistriaveErrpred_diffrelstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpred_diffreliter - popres.DistriaveErrpred_diffrel{iprobe,g}).^2;
            
            
            [~,isort] = sort(avepred_iter,'ascend');
            x_i = [avepred_iter(isort);avepred_iter(isort(end))+avepred_iter(isort);2*avepred_iter(isort(end))+avepred_iter(isort)];
            x_o = 0:(3*max(popres.DistriXavepred{iprobe,2})-1);
            err_i = repmat(aveerrpred_iter(isort),[3 1]);
            err = interp1(x_i,err_i,x_o);
            err = err(floor(numel(err)/3)+1:floor(2*numel(err)/3));
            Xvis = (1:numel(avepred_iter))';
            Xvis_sort = Xvis(isort);
            Xvis_sort(Xvis_sort>Xvis + floor(numel(avepred_iter)/2)) = Xvis_sort(Xvis_sort>Xvis + floor(numel(avepred_iter)/2)) - numel(avepred_iter);
            Xvis_sort(Xvis_sort<Xvis - floor(numel(avepred_iter)/2)) = Xvis_sort(Xvis_sort<Xvis - floor(numel(avepred_iter)/2)) + numel(avepred_iter);
            Xvis_sort = [Xvis_sort;Xvis_sort(end) + Xvis_sort; 2*Xvis_sort(end) + Xvis_sort];
            Xvis_i = interp1(x_i,Xvis_sort,x_o);
            Xvis_i = Xvis_i(floor(numel(Xvis_i)/3)+1:floor(2*numel(Xvis_i)/3)) - Xvis_i(floor(numel(Xvis_i)/3));
            Xvis_i(Xvis_i < 0) = Xvis_i(Xvis_i < 0) + numel(popres.DistriXavepred{iprobe,g});
            Xvis_i(Xvis_i > numel(avepred_iter)) = Xvis_i(Xvis_i > numel(avepred_iter)) - numel(avepred_iter);
            aveerrpred_deciter = err(:);
            popres.DistriaveErrpred_decstd{iprobe,g} = popres.DistriaveErrpred_decstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpred_deciter - popres.DistriaveErrpred_dec{iprobe,g}).^2;
            
            [~,isort] = sort(popres.DistriXavepred{iprobe,2},'ascend');
            x_i = [popres.DistriXavepred{iprobe,2}(isort);popres.DistriXavepred{iprobe,2}(isort(end))+popres.DistriXavepred{iprobe,2}(isort);2*popres.DistriXavepred{iprobe,2}(isort(end))+popres.DistriXavepred{iprobe,2}(isort)];
            x_o = 0:(3*max(popres.DistriXavepred{iprobe,2})-1);
            err_i = repmat(popres.DistriaveErrpred{iprobe,2}(isort),[3 1]);
            err_ref = interp1(x_i,err_i,x_o);
            err_ref = err_ref(floor(numel(err_ref)/3)+1:floor(2*numel(err_ref)/3));
            [~,isortV] = unique(Xvis_i,'sorted');
            aveerrpred_reldeciter = interp1(Xvis_i(isortV),err(isortV) - err_ref(isortV),0:(Xrange/popres.dx - 1));
            aveerrpred_reldeciter = smooth(aveerrpred_reldeciter,round(1/popres.dx));
%             aveerrpred_reldeciter = err(:) - err_ref(:);
            popres.DistriaveErrpred_reldecstd{iprobe,g} = popres.DistriaveErrpred_reldecstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpred_reldeciter - popres.DistriaveErrpred_reldec{iprobe,g}).^2;
            
            nx = numel(aveerrpred_iter,1);
            aveerrpreddiffmat_iter = repmat(aveerrpred_iter,[1 nx]) - repmat(aveerrpred_iter',[nx 1]);
            aveerrpreddiffmat_reliter = repmat(aveerrpred_reliter,[1 nx]) - repmat(aveerrpred_reliter',[nx 1]);
            aveerrpreddiffmat_deciter = repmat(aveerrpred_deciter,[1 nx]) - repmat(aveerrpred_deciter',[nx 1]);
            aveerrpreddiffmat_reldeciter = repmat(aveerrpred_reldeciter,[1 nx]) - repmat(aveerrpred_reldeciter',[nx 1]);
            
            popres.DistriaveErrpreddiffmat_std{iprobe,g} = popres.DistriaveErrpreddiffmat_std{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpreddiffmat_iter - popres.DistriaveErrpreddiffmat{iprobe,g}).^2;
            popres.DistriaveErrpreddiffmat_relstd{iprobe,g} = popres.DistriaveErrpreddiffmat_relstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpreddiffmat_reliter - popres.DistriaveErrpreddiffmat_rel{iprobe,g}).^2;
            popres.DistriaveErrpreddiffmat_decstd{iprobe,g} = popres.DistriaveErrpreddiffmat_decstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpreddiffmat_deciter - popres.DistriaveErrpreddiffmat_dec{iprobe,g}).^2;
            popres.DistriaveErrpreddiffmat_reldecstd{iprobe,g} = popres.DistriaveErrpreddiffmat_reldecstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpreddiffmat_reldeciter - popres.DistriaveErrpreddiffmat_reldec{iprobe,g}).^2;
            
%             [~,isort] = sort(avepred_iter,'ascend');
%             x_i = [avepred_iter(isort);avepred_iter(isort(end))+avepred_iter(isort);2*avepred_iter(isort(end))+avepred_iter(isort)];
%             x_o = 0:(3*max(popres.DistriXavepred{iprobe,2})-1);
%             err_i = repmat(aveerrpred_diffiter(isort),[3 1]);
%             err = interp1(x_i,err_i,x_o,'spline');
%             err = err(floor(numel(err)/3)+1:floor(2*numel(err)/3));
%             popres.DistriaveErrpred_diffdecstd{iprobe,g} = popres.DistriaveErrpred_diffdecstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(err(:) - popres.DistriaveErrpred_diffdec{iprobe,g}).^2;
%             
%             [~,isort] = sort(popres.DistriXavepred{iprobe,2},'ascend');
%             x_i = [popres.DistriXavepred{iprobe,2}(isort);popres.DistriXavepred{iprobe,2}(isort(end))+popres.DistriXavepred{iprobe,2}(isort);2*popres.DistriXavepred{iprobe,2}(isort(end))+popres.DistriXavepred{iprobe,2}(isort)];
%             x_o = 0:(3*max(popres.DistriXavepred{iprobe,2})-1);
%             err_i = repmat(popres.DistriaveErrpred_diff{iprobe,2}(isort),[3 1]);
%             err_ref = interp1(x_i,err_i,x_o,'spline');
%             err_ref = err_ref(floor(numel(err_ref)/3)+1:floor(2*numel(err_ref)/3));
%             aveerrpred_diffreldeciter = err(:) - err_ref(:);
%             popres.DistriaveErrpred_diffreldecstd{iprobe,g} = popres.DistriaveErrpred_diffreldecstd{iprobe,g} + (CVO.kfold - 1)/CVO.kfold*(aveerrpred_diffreldeciter - popres.DistriaveErrpred_diffreldec{iprobe,g}).^2;
%             
%             stdresp = stdresp + (CVO.kfold - 1)/CVO.kfold*(obj.model.tuning(icell).respModel(i,:) - obj.model.tuning(icell).meanrespModel).^2;
        end
    end
end

for iprobe = 1:nProbe
    for g = 1:3
%         popres.PostXpred{iprobe,g} = getCircularAverage(popres.PostXAll{iprobe,g},amp_th_meanpost,maxtol_meanpost);
%         for i = 1:size(popres.PostXAll{iprobe,g},2)
%             popres.PostErrAll{iprobe,g}(:,i) = circshift(popres.PostXAll{iprobe,g}(:,i),-i + floor(size(popres.PostXAll{iprobe,g},2)/2));
%         end
%         popres.PostErrpred{iprobe,g} = popres.PostXpred{iprobe,g} - (1:size(popres.PostXAll{iprobe,g},2))';
%         popres.PostErrpred{iprobe,g}(popres.PostErrpred{iprobe,g}>floor(size(popres.PostXAll{iprobe,g},2)/2)) = popres.PostErrpred{iprobe,g}(popres.PostErrpred{iprobe,g}>floor(size(popres.PostXAll{iprobe,g},2)/2))-size(popres.PostXAll{iprobe,g},2);
%         popres.PostErrpred{iprobe,g}(popres.PostErrpred{iprobe,g}<-floor(size(popres.PostXAll{iprobe,g},2)/2)) = popres.PostErrpred{iprobe,g}(popres.PostErrpred{iprobe,g}<-floor(size(popres.PostXAll{iprobe,g},2)/2)) + size(popres.PostXAll{iprobe,g},2);
%         
%         popres.DistriXmaxpred{iprobe,g} = getCircularAverage(popres.DistriXmaxAll{iprobe,g},amp_th_meandecmax,maxtol_meandecmax);
%         for i = 1:size(popres.DistriXmaxAll{iprobe,g},2)
%             popres.DistrimaxErrAll{iprobe,g}(:,i) = circshift(popres.DistriXmaxAll{iprobe,g}(:,i),-i + floor(size(popres.DistriXmaxAll{iprobe,g},2)/2));
%         end
%         popres.DistrimaxErrpred{iprobe,g} = popres.DistriXmaxpred{iprobe,g} - (1:size(popres.DistriXmaxAll{iprobe,g},2))';
%         popres.DistrimaxErrpred{iprobe,g}(popres.DistrimaxErrpred{iprobe,g}>floor(size(popres.DistriXmaxAll{iprobe,g},2)/2)) = popres.DistrimaxErrpred{iprobe,g}(popres.DistrimaxErrpred{iprobe,g}>floor(size(popres.DistriXmaxAll{iprobe,g},2)/2))-size(popres.DistriXmaxAll{iprobe,g},2);
%         popres.DistrimaxErrpred{iprobe,g}(popres.DistrimaxErrpred{iprobe,g}<-floor(size(popres.DistriXmaxAll{iprobe,g},2)/2)) = popres.DistrimaxErrpred{iprobe,g}(popres.DistrimaxErrpred{iprobe,g}<-floor(size(popres.DistriXmaxAll{iprobe,g},2)/2)) + size(popres.DistriXmaxAll{iprobe,g},2);
%         
%         popres.DistriXavepred{iprobe,g} = getCircularAverage(popres.DistriXaveAll{iprobe,g},amp_th_meandecave,maxtol_meandecave);
%         for i = 1:size(popres.DistriXaveAll{iprobe,g},2)
%             popres.DistriaveErrAll{iprobe,g}(:,i) = circshift(popres.DistriXaveAll{iprobe,g}(:,i),-i + floor(size(popres.DistriXaveAll{iprobe,g},2)/2));
%         end
%         popres.DistriaveErrpred{iprobe,g} = popres.DistriXavepred{iprobe,g} - (1:size(popres.DistriXaveAll{iprobe,g},2))';
%         popres.DistriaveErrpred{iprobe,g}(popres.DistriaveErrpred{iprobe,g}>floor(size(popres.DistriXaveAll{iprobe,g},2)/2)) = popres.DistriaveErrpred{iprobe,g}(popres.DistriaveErrpred{iprobe,g}>floor(size(popres.DistriXaveAll{iprobe,g},2)/2))-size(popres.DistriXaveAll{iprobe,g},2);
%         popres.DistriaveErrpred{iprobe,g}(popres.DistriaveErrpred{iprobe,g}<-floor(size(popres.DistriXaveAll{iprobe,g},2)/2)) = popres.DistriaveErrpred{iprobe,g}(popres.DistriaveErrpred{iprobe,g}<-floor(size(popres.DistriXaveAll{iprobe,g},2)/2)) + size(popres.DistriXaveAll{iprobe,g},2);
    end 
    
    for g = 1:3
        [popres.PostErrpred_reg{iprobe,g},popres.xreg_interp] = Xpredcorrection(popres.PostXpred{iprobe,g},popres.PostXpred{iprobe,2},size(popres.PostXpred{iprobe,2},1),popres.dx);
        [popres.DistrimaxErrpred_reg{iprobe,g},popres.xreg_interp] = Xpredcorrection(popres.DistriXmaxpred{iprobe,g},popres.DistriXmaxpred{iprobe,2},size(popres.DistriXmaxpred{iprobe,2},1),popres.dx);
        [popres.DistriaveErrpred_reg{iprobe,g},popres.xreg_interp] = Xpredcorrection(popres.DistriXavepred{iprobe,g},popres.DistriXavepred{iprobe,2},size(popres.DistriXavepred{iprobe,2},1),popres.dx);
    end
end

% popres.CA1V1corr = cell(3,1);
% popres.CA1V1corr_rand = cell(3,1);
% popres.CA1V1corr_randX = cell(3,1);
% popres.CA1V1corr_randS = cell(3,1);
% popres.CA1V1corr_randXS = cell(3,1);
% 
% popres.CA1V12Dcorr = cell(3,1);
% popres.CA1V12Dcorr_rand = cell(3,1);
% popres.CA1V12Dcorr_randX = cell(3,1);
% popres.CA1V12Dcorr_randS = cell(3,1);
% popres.CA1V12Dcorr_randXS = cell(3,1);
% 
% popres.CA1V12Dcorr_randsqr = cell(3,1);
% popres.CA1V12Dcorr_randXsqr = cell(3,1);
% popres.CA1V12Dcorr_randSsqr = cell(3,1);
% popres.CA1V12Dcorr_randXSsqr = cell(3,1);
% 
% popres.CA1V11Dcorr = cell(3,1);
% popres.CA1V11Dcorr_rand = cell(3,1);
% popres.CA1V11Dcorr_randX = cell(3,1);
% popres.CA1V11Dcorr_randS = cell(3,1);
% popres.CA1V11Dcorr_randXS = cell(3,1);
% Nsession{g} = cell(3,1);
% for g = 1:3
%     popres.CA1V1corr_rand{g} = 0;
%     popres.CA1V1corr_randX{g} = 0;
%     popres.CA1V1corr_randS{g} = 0;
%     popres.CA1V1corr_randXS{g} = 0;
%     
%     popres.CA1V1corr_randsqr{g} = 0;
%     popres.CA1V1corr_randXsqr{g} = 0;
%     popres.CA1V1corr_randSsqr{g} = 0;
%     popres.CA1V1corr_randXSsqr{g} = 0;
% 
%     popres.CA1V12Dcorr{g} = 0;
%     popres.CA1V12Dcorr_rand{g} = 0;
%     popres.CA1V12Dcorr_randX{g} = 0;
%     popres.CA1V12Dcorr_randS{g} = 0;
%     popres.CA1V12Dcorr_randXS{g} = 0;
%     
%     popres.CA1V12Dcorr_randsqr{g} = 0;
%     popres.CA1V12Dcorr_randXsqr{g} = 0;
%     popres.CA1V12Dcorr_randSsqr{g} = 0;
%     popres.CA1V12Dcorr_randXSsqr{g} = 0;
%     
%     popres.CA1V11Dcorr{g} = 0;
%     popres.CA1V11Dcorr_rand{g} = 0;
%     popres.CA1V11Dcorr_randX{g} = 0;
%     popres.CA1V11Dcorr_randS{g} = 0;
%     popres.CA1V11Dcorr_randXS{g} = 0;
%     
%     popres.CA1V11Dcorr_randsqr{g} = 0;
%     popres.CA1V11Dcorr_randXsqr{g} = 0;
%     popres.CA1V11Dcorr_randSsqr{g} = 0;
%     popres.CA1V11Dcorr_randXSsqr{g} = 0;
%     
%     popres.CA1V1mscoh{g} = 0;
%     popres.CA1V1mscoh_rand{g} = 0;
%     popres.CA1V1mscoh_randX{g} = 0;
%     popres.CA1V1mscoh_randS{g} = 0;
%     popres.CA1V1mscoh_randXS{g} = 0;
%     Nsession{g} = 0;
% end
% outvalcorr = 2;%[0 1 2 3 4 5];%[0 1 2 3 4];%
% for g = 1:3
%     Nsession{g} = 0;
% end
% for ianimal = 1:nanimal
%     for iseries = 1:numel(expt(ianimal).series)
%         if (expt(ianimal).goodCA1{iseries} == 1) && (expt(ianimal).goodV1{iseries} == 1)
%             
%             speeds = NaN(size(res.ballspeed{ianimal,iseries}));
%             speeds(~isnan(res.ballspeed{ianimal,iseries})) = smthInTime(res.ballspeed{ianimal,iseries}(~isnan(res.ballspeed{ianimal,iseries})), sampleRate, SpdSmthWin, 'same', [], 'boxcar_centered');
%             
%             XErrCA1 = res.Xpred_ave{ianimal,iseries,1};% - res.X{ianimal,iseries};%res.Xpred_max{ianimal,iseries,1} - res.X{ianimal,iseries};
% %             XErrCA1(XErrCA1>floor(Xrange/2)) = XErrCA1(XErrCA1>floor(Xrange/2)) - Xrange;
% %             XErrCA1(XErrCA1<-floor(Xrange/2)) = XErrCA1(XErrCA1<-floor(Xrange/2)) + Xrange;
%             
%             XErrV1 = res.Xpred_ave{ianimal,iseries,2};% - res.X{ianimal,iseries};%res.Xpred_max{ianimal,iseries,2} - res.X{ianimal,iseries};
% %             XErrV1(XErrV1>floor(Xrange/2)) = XErrV1(XErrV1>floor(Xrange/2)) - Xrange;
% %             XErrV1(XErrV1<-floor(Xrange/2)) = XErrV1(XErrV1<-floor(Xrange/2)) + Xrange;
%             
%             for g = [2 1 3]
%                 cont_list = find(ismember(res.contrastVal{ianimal,iseries}, contval));
%                 RL_list = find(ismember(res.roomlengthVal{ianimal,iseries}, [1]));
%                 outcome_list = find(ismember(res.outcomeVal{ianimal,iseries}, outvalcorr));
%                 tidx = false(size(res.tidx{ianimal,iseries,1,1,2,1,3}));
%                 for cont = 1:numel(cont_list)
%                     for r = 1:numel(RL_list)
%                         for o = 1:numel(outcome_list)
%                             if ~isempty(res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)})
%                                 tidx = tidx | (res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)} & ismember(res.Phsbin{ianimal,iseries},phsbins));% & res.X{ianimal,iseries} > 0 & res.X{ianimal,iseries} <= 10);
%                             end
%                         end
%                     end
%                 end
%                 tidx0 = tidx & ismember(res.Spdbin2{ianimal,iseries},spdbins) & ~isnan(XErrCA1) & ~isnan(XErrV1);
% %                 [Fave_CA1, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(~isnan(XErrCA1) & tidx0) XErrCA1(~isnan(XErrCA1) & tidx0)], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                 [Fave_V1, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(~isnan(XErrV1) & tidx0) XErrV1(~isnan(XErrV1) & tidx0)], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                 Fave_CA1err = getCircularAverage(Fave_CA1,amp_th_meandecave,maxtol_meandecave) - floor(Xrange/2);
% %                 Fave_V1err = getCircularAverage(Fave_V1,amp_th_meandecave,maxtol_meandecave) - floor(Xrange/2);
% 
%                 for ispd = 1:numel(spdbins)
%                 for xx = 1:Xrange
%                     tidx = tidx0 & (res.X{ianimal,iseries} == xx) & (res.Spdbin2{ianimal,iseries} == spdbins(ispd));
% %                     XErrCA1(tidx) = XErrCA1(tidx) - Fave_CA1err(xx);
% %                     XErrCA1(tidx & XErrCA1>floor(Xrange/2)) = XErrCA1(tidx & XErrCA1>floor(Xrange/2)) - Xrange;
% %                     XErrCA1(tidx & XErrCA1<-floor(Xrange/2)) = XErrCA1(tidx & XErrCA1<-floor(Xrange/2)) + Xrange;
% %                     XErrV1(tidx) = XErrV1(tidx) - Fave_V1err(xx);
% %                     XErrV1(tidx & XErrV1>floor(Xrange/2)) = XErrV1(tidx & XErrV1>floor(Xrange/2)) - Xrange;
% %                     XErrV1(tidx & XErrV1<-floor(Xrange/2)) = XErrV1(tidx & XErrV1<-floor(Xrange/2)) + Xrange;
% %                 end
% 
% %                 tidx = tidx0;
%                             
%                 if sum(tidx)>0
%                     Nsession{g} = Nsession{g} + 1;
%                     idx = find(tidx);
% %                     XErrV1CA1 = XErrV1(idx) - XErrCA1(idx);
% %                     XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
% %                     XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
% %                     CA1V1noisecorr = sum((XErrCA1(idx)+XErrV1CA1 - nanmean((XErrCA1(idx)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx)+XErrV1CA1 - nanmean((XErrCA1(idx)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
% % 
%                     XerrV1ang = (XErrV1(idx)/Xrange)*2*pi;
%                     a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
%                     XerrV1ang_mean = atan2(b,a);
%                     XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
%                     a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
%                     XerrCA1ang_mean = atan2(b,a);
% %                     [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
%                     CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;
% 
% %                     CA1V1noisecorr = sum(XErrV1(idx).*XErrCA1(idx))/(sum(XErrV1(idx).^2)*sum(XErrCA1(idx).^2))^0.5;
%                     popres.CA1V1corr{g} = [popres.CA1V1corr{g} CA1V1noisecorr];
%                     
%                     Xwincorr = (-1:1);%(-180:180);
%                     xcorrtemp = zeros(numel(Xwincorr),1);
%                     for tt_corr = Xwincorr
%                         XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
%                         XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
%                         xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
%                     end
%                     popres.CA1V11Dcorr{g} = popres.CA1V11Dcorr{g} + xcorrtemp;%xcorr(XErrCA1(idx),XErrV1(idx),3*60,'coeff');
%                     %                     popres.CA1V1mscoh{g} = popres.CA1V1mscoh{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx),512,128));
%                     
%                     
% %                     [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                     [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
%                     [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx)], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
% %                     Fave = Fave/sum(Fave(:));
%                     popres.CA1V12Dcorr{g} = popres.CA1V12Dcorr{g} + (Fave./(mean(Fave,2)*mean(Fave,1)));%Fave-(mean(Fave,2)*mean(Fave,1));%Fave./repmat(sum(Fave,1),[Xrange 1]);%
%                 else
%                     popres.CA1V1corr{g} = [popres.CA1V1corr{g} NaN];
%                 end
%                 end
%                 end
%             end
%         end
%     end
% end
% for g = 1:3
%     popres.CA1V12Dcorr{g} = popres.CA1V12Dcorr{g}/Nsession{g};
%     popres.CA1V11Dcorr{g} = popres.CA1V11Dcorr{g}/Nsession{g};
% end

% Nranditer = 1;%
% for kiter = 1:Nranditer
%     tempCA1V1corr_rand = cell(3,1);
%     tempCA1V1corr_randX = cell(3,1);
%     tempCA1V1corr_randS = cell(3,1);
%     tempCA1V1corr_randXS = cell(3,1);
%     for g = 1:3
%         tempCA1V12Dcorr{g} = 0;
%         tempCA1V12Dcorr_rand{g} = 0;
%         tempCA1V12Dcorr_randX{g} = 0;
%         tempCA1V12Dcorr_randS{g} = 0;
%         tempCA1V12Dcorr_randXS{g} = 0;
%         
%         tempCA1V11Dcorr{g} = 0;
%         tempCA1V11Dcorr_rand{g} = 0;
%         tempCA1V11Dcorr_randX{g} = 0;
%         tempCA1V11Dcorr_randS{g} = 0;
%         tempCA1V11Dcorr_randXS{g} = 0;
%         
%         tempCA1V1mscoh{g} = 0;
%         tempCA1V1mscoh_rand{g} = 0;
%         tempCA1V1mscoh_randX{g} = 0;
%         tempCA1V1mscoh_randS{g} = 0;
%         tempCA1V1mscoh_randXS{g} = 0;
%         Nsession{g} = 0;
%     end
%     for ianimal = 1:nanimal
%         for iseries = 1:numel(expt(ianimal).series)
%             if (expt(ianimal).goodCA1{iseries} == 1) && (expt(ianimal).goodV1{iseries} == 1)
%                 
%                 speeds = NaN(size(res.ballspeed{ianimal,iseries}));
%                 speeds(~isnan(res.ballspeed{ianimal,iseries})) = smthInTime(res.ballspeed{ianimal,iseries}(~isnan(res.ballspeed{ianimal,iseries})), sampleRate, SpdSmthWin, 'same', [], 'boxcar_centered');
%                 
%                 XErrCA1 = res.Xpred_ave{ianimal,iseries,1};% - res.X{ianimal,iseries};%res.Xpred_max{ianimal,iseries,1} - res.X{ianimal,iseries};
% %                 XErrCA1(XErrCA1>floor(Xrange/2)) = XErrCA1(XErrCA1>floor(Xrange/2)) - Xrange;
% %                 XErrCA1(XErrCA1<-floor(Xrange/2)) = XErrCA1(XErrCA1<-floor(Xrange/2)) + Xrange;
%                 
%                 XErrV1 = res.Xpred_ave{ianimal,iseries,2};% - res.X{ianimal,iseries};%res.Xpred_max{ianimal,iseries,2} - res.X{ianimal,iseries};
% %                 XErrV1(XErrV1>floor(Xrange/2)) = XErrV1(XErrV1>floor(Xrange/2)) - Xrange;
% %                 XErrV1(XErrV1<-floor(Xrange/2)) = XErrV1(XErrV1<-floor(Xrange/2)) + Xrange;
%                 
%                 %             posErrCA1 = 2*randn(size(XErrCA1));
%                 %             posErrV1 = posErrCA1;
%                 %             posErrCA1 = posErrCA1 + randn(size(posErrCA1))*3;
%                 %             posErrV1 = posErrV1 + randn(size(posErrV1))*3;
%                 %
%                 %             XErrCA1 = posErrCA1 - nanmean(posErrCA1);
%                 %             XErrV1 = posErrV1 - nanmean(posErrV1);
%                 %             Xshift = [5 0 -5];
%                 
%                 for g = 1:3
%                     cont_list = find(ismember(res.contrastVal{ianimal,iseries}, contval));
%                     RL_list = find(ismember(res.roomlengthVal{ianimal,iseries}, [1]));
%                     outcome_list = find(ismember(res.outcomeVal{ianimal,iseries}, outvalcorr));
%                     tidx = false(size(res.tidx{ianimal,iseries,1,1,2,1,3}));
%                     for cont = 1:numel(cont_list)
%                         for RF = 1:numel(RL_list)
%                             for o = 1:numel(outcome_list)
%                                 if ~isempty(res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)})
%                                     tidx = tidx | (res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)} & ismember(res.Phsbin{ianimal,iseries},phsbins));% & res.X{ianimal,iseries} > 0 & res.X{ianimal,iseries} <= 10);
%                                 end
%                             end
%                         end
%                     end
%                     tidx0 = tidx & ismember(res.Spdbin2{ianimal,iseries},spdbins) & ~isnan(XErrCA1) & ~isnan(XErrV1);
% %                     [Fave_CA1, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(~isnan(XErrCA1) & tidx0) XErrCA1(~isnan(XErrCA1) & tidx0)], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                     [Fave_V1, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(~isnan(XErrV1) & tidx0) XErrV1(~isnan(XErrV1) & tidx0)], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                     Fave_CA1err = getCircularAverage(Fave_CA1,amp_th_meandecave,maxtol_meandecave) - floor(Xrange/2);
% %                     Fave_V1err = getCircularAverage(Fave_V1,amp_th_meandecave,maxtol_meandecave) - floor(Xrange/2);
%                     for ispd = 1:numel(spdbins)    
%                     for xx = 1:Xrange
%                         tidx = tidx0 & (res.X{ianimal,iseries} == xx) & (res.Spdbin2{ianimal,iseries} == spdbins(ispd));
% %                         XErrCA1(tidx) = XErrCA1(tidx) - Fave_CA1err(xx);
% %                         XErrCA1(tidx & XErrCA1>floor(Xrange/2)) = XErrCA1(tidx & XErrCA1>floor(Xrange/2)) - Xrange;
% %                         XErrCA1(tidx & XErrCA1<-floor(Xrange/2)) = XErrCA1(tidx & XErrCA1<-floor(Xrange/2)) + Xrange;
% %                         XErrV1(tidx) = XErrV1(tidx) - Fave_V1err(xx);
% %                         XErrV1(tidx & XErrV1>floor(Xrange/2)) = XErrV1(tidx & XErrV1>floor(Xrange/2)) - Xrange;
% %                         XErrV1(tidx & XErrV1<-floor(Xrange/2)) = XErrV1(tidx & XErrV1<-floor(Xrange/2)) + Xrange;
% %                     end
% 
% %                     tidx = tidx0;
%                     
%                     if sum(tidx)>0
%                         Nsession{g} = Nsession{g} + 1;
%                         idx = find(tidx);
%                         idx_rand = idx(randperm(numel(idx)));
% %                         XErrV1CA1 = XErrV1(idx_rand) - XErrCA1(idx_rand);
% %                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
% %                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
% %                         CA1V1noisecorr = sum((XErrCA1(idx_rand)+XErrV1CA1 - nanmean((XErrCA1(idx_rand)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx_rand)+XErrV1CA1 - nanmean((XErrCA1(idx_rand)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
%                         
%                         XerrV1ang = (XErrV1(idx_rand)/Xrange)*2*pi;
%                         a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
%                         XerrV1ang_mean = atan2(b,a);
%                         XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
%                         a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
%                         XerrCA1ang_mean = atan2(b,a);
% %                         [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
%                         CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;
%                         
% %                         CA1V1noisecorr = sum(XErrV1(idx_rand).*XErrCA1(idx))/(sum(XErrV1(idx_rand).^2)*sum(XErrCA1(idx).^2))^0.5;
% 
%                         Xwincorr = (-1:1);%(-180:180);
%                         xcorrtemp = zeros(numel(Xwincorr),1);
%                         for tt_corr = Xwincorr
%                             XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
%                             XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
%                             xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
%                         end
%                         tempCA1V11Dcorr_rand{g} = tempCA1V11Dcorr_rand{g} + xcorrtemp;%tempCA1V11Dcorr_rand{g} = tempCA1V11Dcorr_rand{g} + xcorr(XErrCA1(idx),XErrV1(idx_rand),3*60,'coeff');
%                         %                         popres.CA1V1mscoh_rand{g} = popres.CA1V1mscoh_rand{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx_rand),512,128));
%                         
% %                         XErrV1CA1 = res.Xpred_ave{ianimal,iseries,2}(idx_rand) - res.Xpred_ave{ianimal,iseries,1}(idx);
% %                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
% %                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
% %                         [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% 
% %                         [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx_rand)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
%                         [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx_rand)], lambdaSmooth, [Xrange Xrange], 1:Xrange, 1:Xrange, true, true);
%                         tempCA1V12Dcorr_rand{g} = tempCA1V12Dcorr_rand{g} + (Fave - mean(Fave,2)*mean(Fave,1));
%                         tempCA1V1corr_rand{g} = [tempCA1V1corr_rand{g} CA1V1noisecorr];
%                         
%                         
% %                         idx = find(tidx);
% %                         idx_randX = zeros(size(idx));
% %                         for ix = 1:max(res.X{ianimal,iseries})
% %                             idxtemp = idx(res.X{ianimal,iseries}(idx) == ix);
% %                             idx_randX(res.X{ianimal,iseries}(idx) == ix) = idxtemp(randperm(numel(idxtemp)));
% %                         end
% % %                         XErrV1CA1 = XErrV1(idx_randX) - XErrCA1(idx_randX);
% % %                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
% % %                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
% % %                         CA1V1noisecorr = sum((XErrCA1(idx_randX)+XErrV1CA1 - nanmean((XErrCA1(idx_randX)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx_randX)+XErrV1CA1 - nanmean((XErrCA1(idx_randX)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
% %                         
% %                         XerrV1ang = (XErrV1(idx_randX)/Xrange)*2*pi;
% %                         a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
% %                         XerrV1ang_mean = atan2(b,a);
% %                         XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
% %                         a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
% %                         XerrCA1ang_mean = atan2(b,a);
% % %                         [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
% %                         CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;
% %                         
% % %                         CA1V1noisecorr = sum(XErrV1(idx_randX).*XErrCA1(idx))/(sum(XErrV1(idx_randX).^2)*sum(XErrCA1(idx).^2))^0.5;
% %                         Xwincorr = (-180:180);
% %                         xcorrtemp = zeros(numel(Xwincorr),1);
% %                         for tt_corr = Xwincorr
% %                             XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
% %                             XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
% %                             xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
% %                         end
% %                         tempCA1V11Dcorr_randX{g} = tempCA1V11Dcorr_randX{g} + xcorrtemp;%tempCA1V11Dcorr_randX{g} = tempCA1V11Dcorr_randX{g} + xcorr(XErrCA1(idx),XErrV1(idx_randX),3*60,'coeff');
% %                         %                         popres.CA1V1mscoh_randX{g} = popres.CA1V1mscoh_randX{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx_randX),512,128));
% %                         
% % %                         XErrV1CA1 = res.Xpred_ave{ianimal,iseries,2}(idx_randX) - res.Xpred_ave{ianimal,iseries,1}(idx);
% % %                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
% % %                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
% % %                         [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                         [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx_randX)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                         tempCA1V12Dcorr_randX{g} = tempCA1V12Dcorr_randX{g} + Fave;
% %                         tempCA1V1corr_randX{g} = [tempCA1V1corr_randX{g} CA1V1noisecorr];
% %                         
% %                         idx = find(tidx);
% %                         idx_randSpeed = zeros(size(idx));
% %                         Speedrange = [min(speeds(idx)):(max(speeds(idx))-min(speeds(idx)))/10:max(speeds(idx))];
% %                         for ispeed = 1:10
% %                             idxtemp = idx(speeds(idx) >= Speedrange(ispeed) & speeds(idx) <= Speedrange(ispeed+1));
% %                             idx_randSpeed(speeds(idx) >= Speedrange(ispeed) & speeds(idx) <= Speedrange(ispeed+1)) = idxtemp(randperm(numel(idxtemp)));
% %                         end
% % %                         XErrV1CA1 = XErrV1(idx_randSpeed) - XErrCA1(idx_randSpeed);
% % %                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
% % %                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
% % %                         CA1V1noisecorr = sum((XErrCA1(idx_randSpeed)+XErrV1CA1 - nanmean((XErrCA1(idx_randSpeed)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx_randSpeed)+XErrV1CA1 - nanmean((XErrCA1(idx_randSpeed)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
% %                         
% %                         XerrV1ang = (XErrV1(idx_randSpeed)/Xrange)*2*pi;
% %                         a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
% %                         XerrV1ang_mean = atan2(b,a);
% %                         XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
% %                         a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
% %                         XerrCA1ang_mean = atan2(b,a);
% % %                         [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
% %                         CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;
% %                         
% % %                         CA1V1noisecorr = sum(XErrV1(idx_randSpeed).*XErrCA1(idx))/(sum(XErrV1(idx_randSpeed).^2)*sum(XErrCA1(idx).^2))^0.5;
% %                         Xwincorr = (-180:180);
% %                         xcorrtemp = zeros(numel(Xwincorr),1);
% %                         for tt_corr = Xwincorr
% %                             XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
% %                             XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
% %                             xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
% %                         end
% %                         tempCA1V11Dcorr_randS{g} = tempCA1V11Dcorr_randS{g} + xcorrtemp;%tempCA1V11Dcorr_randS{g} = tempCA1V11Dcorr_randS{g} + xcorr(XErrCA1(idx),XErrV1(idx_randSpeed),3*60,'coeff');
% %                         %                         popres.CA1V1mscoh_randS{g} = popres.CA1V1mscoh_randS{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx_randSpeed),512,128));
% %                         
% % %                         XErrV1CA1 = res.Xpred_ave{ianimal,iseries,2}(idx_randSpeed) - res.Xpred_ave{ianimal,iseries,1}(idx);
% % %                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
% % %                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
% % %                         [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                         [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx_randSpeed)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                         tempCA1V12Dcorr_randS{g} = tempCA1V12Dcorr_randS{g} + Fave;
% %                         tempCA1V1corr_randS{g} = [tempCA1V1corr_randS{g} CA1V1noisecorr];
% %                         
% %                         idx = find(tidx);
% %                         idx_randXS = zeros(size(idx));
% %                         for ispeed = 1:nSpdbins
% %                             for ix = min(res.X{ianimal,iseries}(idx)):max(res.X{ianimal,iseries}(idx))
% %                                 idxtemp = idx(res.X{ianimal,iseries}(idx) == ix & res.Spdbin2{ianimal,iseries}(idx) == ispeed);
% %                                 idx_randXS(res.X{ianimal,iseries}(idx) == ix & res.Spdbin2{ianimal,iseries}(idx) == ispeed) = idxtemp(randperm(numel(idxtemp)));
% %                             end
% %                         end
% % %                         XErrV1CA1 = XErrV1(idx_randXS) - XErrCA1(idx_randXS);
% % %                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
% % %                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
% % %                         CA1V1noisecorr = sum((XErrCA1(idx_randXS)+XErrV1CA1 - nanmean((XErrCA1(idx_randXS)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx_randXS)+XErrV1CA1 - nanmean((XErrCA1(idx_randXS)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
% %                         
% %                         XerrV1ang = (XErrV1(idx_randXS)/Xrange)*2*pi;
% %                         a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
% %                         XerrV1ang_mean = atan2(b,a);
% %                         XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
% %                         a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
% %                         XerrCA1ang_mean = atan2(b,a);
% % %                         [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
% %                         CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;
% %                         
% % %                         CA1V1noisecorr = sum(XErrV1(idx_randXS).*XErrCA1(idx))/(sum(XErrV1(idx_randXS).^2)*sum(XErrCA1(idx).^2))^0.5;
% %                         Xwincorr = (-180:180);
% %                         xcorrtemp = zeros(numel(Xwincorr),1);
% %                         for tt_corr = Xwincorr
% %                             XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
% %                             XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
% %                             xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
% %                         end
% %                         tempCA1V11Dcorr_randXS{g} = tempCA1V11Dcorr_randXS{g} + xcorrtemp;%tempCA1V11Dcorr_randXS{g} = tempCA1V11Dcorr_randXS{g} + xcorr(XErrCA1(idx),XErrV1(idx_randXS),3*60,'coeff');
% %                         %                         popres.CA1V1mscoh_randXS{g} = popres.CA1V1mscoh_randXS{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx_randXS),512,128));
% %                         
% % %                         XErrV1CA1 = res.Xpred_ave{ianimal,iseries,2}(idx_randXS) - res.Xpred_ave{ianimal,iseries,1}(idx);
% % %                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
% % %                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
% % %                         [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                         [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx_randXS)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
% %                         tempCA1V12Dcorr_randXS{g} = tempCA1V12Dcorr_randXS{g} + Fave;
% %                         tempCA1V1corr_randXS{g} = [tempCA1V1corr_randXS{g} CA1V1noisecorr];
%                     else
%                         tempCA1V1corr_rand{g} = [tempCA1V1corr_rand{g} NaN];
% %                         tempCA1V1corr_randX{g} = [tempCA1V1corr_randX{g} NaN];
% %                         tempCA1V1corr_randS{g} = [tempCA1V1corr_randS{g} NaN];
% %                         tempCA1V1corr_randXS{g} = [tempCA1V1corr_randXS{g} NaN];
%                     end
%                     end
%                     end
%                 end
%             end
%         end
%     end
%     for g = 1:3
%         tempCA1V12Dcorr_rand{g} = tempCA1V12Dcorr_rand{g}/Nsession{g};
% %         tempCA1V12Dcorr_randX{g} = tempCA1V12Dcorr_randX{g}/Nsession{g};
% %         tempCA1V12Dcorr_randS{g} = tempCA1V12Dcorr_randS{g}/Nsession{g};
% %         tempCA1V12Dcorr_randXS{g} = tempCA1V12Dcorr_randXS{g}/Nsession{g};
%         
%         tempCA1V11Dcorr_rand{g} = tempCA1V11Dcorr_rand{g}/Nsession{g};
% %         tempCA1V11Dcorr_randX{g} = tempCA1V11Dcorr_randX{g}/Nsession{g};
% %         tempCA1V11Dcorr_randS{g} = tempCA1V11Dcorr_randS{g}/Nsession{g};
% %         tempCA1V11Dcorr_randXS{g} = tempCA1V11Dcorr_randXS{g}/Nsession{g};
%         
%         popres.CA1V12Dcorr_rand{g} = popres.CA1V12Dcorr_rand{g} + tempCA1V12Dcorr_rand{g}/Nranditer;
% %         popres.CA1V12Dcorr_randX{g} = popres.CA1V12Dcorr_randX{g} + tempCA1V12Dcorr_randX{g}/Nranditer;
% %         popres.CA1V12Dcorr_randS{g} = popres.CA1V12Dcorr_randS{g} + tempCA1V12Dcorr_randS{g}/Nranditer;
% %         popres.CA1V12Dcorr_randXS{g} = popres.CA1V12Dcorr_randXS{g} + tempCA1V12Dcorr_randXS{g}/Nranditer;
%         
%         popres.CA1V12Dcorr_randsqr{g} = popres.CA1V12Dcorr_randsqr{g} + tempCA1V12Dcorr_rand{g}.^2/Nranditer;
% %         popres.CA1V12Dcorr_randXsqr{g} = popres.CA1V12Dcorr_randXsqr{g} + tempCA1V12Dcorr_randX{g}.^2/Nranditer;
% %         popres.CA1V12Dcorr_randSsqr{g} = popres.CA1V12Dcorr_randSsqr{g} + tempCA1V12Dcorr_randS{g}.^2/Nranditer;
% %         popres.CA1V12Dcorr_randXSsqr{g} = popres.CA1V12Dcorr_randXSsqr{g} + tempCA1V12Dcorr_randXS{g}.^2/Nranditer;
%         
%         popres.CA1V11Dcorr_rand{g} = popres.CA1V11Dcorr_rand{g} + tempCA1V11Dcorr_rand{g}/Nranditer;
% %         popres.CA1V11Dcorr_randX{g} = popres.CA1V11Dcorr_randX{g} + tempCA1V11Dcorr_randX{g}/Nranditer;
% %         popres.CA1V11Dcorr_randS{g} = popres.CA1V11Dcorr_randS{g} + tempCA1V11Dcorr_randS{g}/Nranditer;
% %         popres.CA1V11Dcorr_randXS{g} = popres.CA1V11Dcorr_randXS{g} + tempCA1V11Dcorr_randXS{g}/Nranditer;
%         
%         popres.CA1V11Dcorr_randsqr{g} = popres.CA1V11Dcorr_randsqr{g} + tempCA1V11Dcorr_rand{g}.^2/Nranditer;
% %         popres.CA1V11Dcorr_randXsqr{g} = popres.CA1V11Dcorr_randXsqr{g} + tempCA1V11Dcorr_randX{g}.^2/Nranditer;
% %         popres.CA1V11Dcorr_randSsqr{g} = popres.CA1V11Dcorr_randSsqr{g} + tempCA1V11Dcorr_randS{g}.^2/Nranditer;
% %         popres.CA1V11Dcorr_randXSsqr{g} = popres.CA1V11Dcorr_randXSsqr{g} + tempCA1V11Dcorr_randXS{g}.^2/Nranditer;
%         
%         popres.CA1V1corr_rand{g} = popres.CA1V1corr_rand{g} + tempCA1V1corr_rand{g}/Nranditer;
% %         popres.CA1V1corr_randX{g} = popres.CA1V1corr_randX{g} + tempCA1V1corr_randX{g}/Nranditer;
% %         popres.CA1V1corr_randS{g} = popres.CA1V1corr_randS{g} + tempCA1V1corr_randS{g}/Nranditer;
% %         popres.CA1V1corr_randXS{g} = popres.CA1V1corr_randXS{g} + tempCA1V1corr_randXS{g}/Nranditer;
%         
%         popres.CA1V1corr_randsqr{g} = popres.CA1V1corr_randsqr{g} + tempCA1V1corr_rand{g}.^2/Nranditer;
% %         popres.CA1V1corr_randXsqr{g} = popres.CA1V1corr_randXsqr{g} + tempCA1V1corr_randX{g}.^2/Nranditer;
% %         popres.CA1V1corr_randSsqr{g} = popres.CA1V1corr_randSsqr{g} + tempCA1V1corr_randS{g}.^2/Nranditer;
% %         popres.CA1V1corr_randXSsqr{g} = popres.CA1V1corr_randXSsqr{g} + tempCA1V1corr_randXS{g}.^2/Nranditer;
%         
%         %     popres.CA1V1mscoh{g} = popres.CA1V1mscoh{g}/Nsession{g};
%         %     popres.CA1V1mscoh_rand{g} = popres.CA1V1mscoh_rand{g}/Nsession{g}/Nranditer;
%         %     popres.CA1V1mscoh_randX{g} = popres.CA1V1mscoh_randX{g}/Nsession{g}/Nranditer;
%         %     popres.CA1V1mscoh_randS{g} = popres.CA1V1mscoh_randS{g}/Nsession{g}/Nranditer;
%         %     popres.CA1V1mscoh_randXS{g} = popres.CA1V1mscoh_randXS{g}/Nsession{g}/Nranditer;
%     end
% end



% maxtol = 1;
% nanimal = size(res.PostX,1);
% res = getresSEM(res,maxtol);
% for iprobe = 1:2
%     for g = 1:3
%         res.XpredAll{iprobe,g} = getCircularAverage(res.PostXAll{iprobe,g},maxtol);%getCircularAverage(res.DistriXAll{iprobe,g},maxtol);
%         res.MeanXErrAll{iprobe,g} = nanmean(res.PostErrXAll{iprobe,g},2);
%     end
% end
% 
% % for ianimal = 1:nanimal
% %     figure;
% %     for iprobe = 1:2
% %         for g = 1:3
% %             if sum(res.PostX{ianimal,iprobe,g}(:))>0
% %                 subplot(2,3,(iprobe-1)*3 + g);
% %                 imagesc(res.PostX{ianimal,iprobe,g});%imagesc(res.DistriXAll{iprobe,g});%
% %                 set(gca,'Ydir','normal','Clim',[0.5 1.5],'PlotBoxAspectRatio', [1 1 1])
% %                 colormap(RedWhiteBlue)
% %                 hold on;plot([0 100],[0 100],'w')
% %                 Xpred = getCircularAverage(res.PostX{ianimal,iprobe,g},maxtol);
% %                 hold on;plot(Xpred,'k');
% %             end
% %         end
% %     end
% % end
% c{1} = 'c';
% c{2} = 'k';
% c{3} = 'm';
% titlestr{1} = 'low';
% titlestr{2} = 'med';
% titlestr{3} = 'high';
% probestr{1} = 'CA1';
% probestr{2} = 'V1';
% for iprobe = 1:2
%     figure('name',['Decoding : ' probestr{iprobe}]);
%     for g = 1:3
%         subplot(3,3,g);
%         imagesc(res.PostXAll{iprobe,g});%imagesc(res.DistriXAll{iprobe,g});%
%         set(gca,'Ydir','normal','Clim',[0.5 1.5],'PlotBoxAspectRatio', [1 1 1])
%         colormap(RedWhiteBlue)
%         hold on;plot([0 100],[0 100],'w')
%         hold on;plot(res.XpredAll{iprobe,g},'k');
%         title([probestr{iprobe} ' ' titlestr{g}]);
%         
%         subplot(3,3,4);
%         XErr = res.XpredAll{iprobe,g} - (1:100)';
%         XErr(XErr>50) = XErr(XErr>50) - 100;
%         XErr(XErr<-50) = XErr(XErr<-50) + 100;
%         hold on;ciplot(gca,1:100,XErr,res.XpredAllsem{iprobe,g}(:)',0.2,c{g});%plot(Xpred(:) - (1:100)');
% %         hold on;plot(XErr,c{g});
%         set(gca,'Ylim',[-20 20]);
%         
%         subplot(3,3,5);
%         [~,imax] = max(res.MeanXErrAll{iprobe,2});
%         imax = imax-1;
% %         hold on;ciplot(gca,-50:49,circshift(res.MeanXErrAll{iprobe,g},50-imax),circshift(res.MeanXErrAllsem{iprobe,g},50-imax),0.2,c{g});
%         hold on;plot(-50:49,circshift(res.MeanXErrAll{iprobe,g},50-imax),c{g});
%         set(gca,'Xlim',[-20 20]);
%         [~,imax] = max(circshift(res.MeanXErrAll{iprobe,g},50-imax));
%         hold on;plot([imax-51 imax-51],[0 4],c{g});
%         set(gca,'Xlim',[-20 20]);
%     end
%     shift_low = [];
%     shift_high = [];
%     OptLat_low = [];
%     OptLat_high = [];
%     for ianimal = 1:size(res.MeanXErrMaxPos,1)
%         shift_low = [shift_low res.MeanXErrMaxPos{ianimal,iprobe,1}-res.MeanXErrMaxPos{ianimal,iprobe,2}];
%         shift_high = [shift_high res.MeanXErrMaxPos{ianimal,iprobe,3}-res.MeanXErrMaxPos{ianimal,iprobe,2}];
%         OptLat_low = [OptLat_low res.Latopt{ianimal,iprobe,1}];
%         OptLat_high = [OptLat_high res.Latopt{ianimal,iprobe,3}];
%     end
%     subplot(3,3,6);
%     hold on;
%     hlow = histogram(shift_low,-5:0.5:5,'EdgeColor','none','FaceColor',c{1});
%     hhigh = histogram(shift_high,-5:0.5:5,'EdgeColor','none','FaceColor',c{3});
%     set(gca,'Xlim',[-5 5],'Ylim',[0 18],'PlotBoxAspectRatio', [1 1 1]);
%     hold on;plot([0 0],[0 18],'k');
%     
%     subplot(3,3,8);
%     hold on;
%     hlow = histogram([OptLat_low OptLat_high],-150:100:750,'EdgeColor','none','FaceColor',c{1});
%     hhigh = histogram(OptLat_high,-150:100:7500,'EdgeColor','none','FaceColor',c{3});
%     set(gca,'Xlim',[-200 1000]);
%     
%     subplot(3,3,7);
%     hold on;
%     meanVSrep = nanmean(cellprop.globalXposrep(cellprop.Probe == iprobe & cellprop.Goodcluster ,16:40)+cellprop.globalYposrep(cellprop.Probe == iprobe & cellprop.Goodcluster ,16:40),1);
%     xtimes = (-4:20)*50;
%     stim = zeros(size(meanVSrep));
%     stim(5:end) = sin(xtimes(5:end)/1000*2*pi);
%     plot(xtimes,meanVSrep);
%     plot(xtimes,nanmean(meanVSrep)+(max(meanVSrep)-nanmean(meanVSrep))*stim);
%     set(gca,'Xlim',[xtimes(1) xtimes(end)]);
%     
%     if iprobe == 2
%         V1gainshift = cell(1,3);
%         CA1gainshift = cell(1,3);
%         subplot(3,3,9);
%         for ianimal = 1:size(res.MeanXErrMaxPos,1)
%             CA1V1series = ~isnan(res.MeanXErrMaxPos{ianimal,1,2}) & ~isnan(res.MeanXErrMaxPos{ianimal,2,2});
%             if sum(CA1V1series)>0
%                 for g = [1 3]
%                     hold on;scatter(res.MeanXErrMaxPos{ianimal,1,g}(CA1V1series) - res.MeanXErrMaxPos{ianimal,1,2}(CA1V1series),res.MeanXErrMaxPos{ianimal,2,g}(CA1V1series) - res.MeanXErrMaxPos{ianimal,2,2}(CA1V1series),...
%                         'MarkerEdgeColor',c{g},'MarkerFaceColor',c{g},'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.2);
%                     V1gainshift{g} = [V1gainshift{g} res.MeanXErrMaxPos{ianimal,2,g}(CA1V1series) - res.MeanXErrMaxPos{ianimal,2,2}(CA1V1series)];
%                     CA1gainshift{g} = [CA1gainshift{g} res.MeanXErrMaxPos{ianimal,1,g}(CA1V1series) - res.MeanXErrMaxPos{ianimal,1,2}(CA1V1series)];
%                 end
%             end
%         end
%         CA1nanidx = isnan(CA1gainshift{g});
%         V1nanidx = isnan(V1gainshift{g});
%         CA1gainshift{g} = CA1gainshift{g}(~CA1nanidx & ~V1nanidx);
%         V1gainshift{g} = V1gainshift{g}(~CA1nanidx & ~V1nanidx);
%         for g = [1 3]
%             [r,p] = corr(CA1gainshift{g}',V1gainshift{g}');
%             rho{g} = r; pval{g} = p;
%         end
%         hold on;plot([-5 5],[-5 5],'k');
%         title(['rlow = ' num2str(rho{1}) ' rhigh = ' num2str(rho{3})]);
%         set(gca,'Xlim',[-5 5],'Ylim',[-5 5],'XAxisLocation','origin','YAxisLocation','origin','PlotBoxAspectRatio', [1 1 1]);
%     end
% end
% 
% figure('Name','CA1V1 correlation')
% Xxcorr = ((1:numel(res.CA1V1corrAll{2})) - floor(numel(res.CA1V1corrAll{2})/2)+1)*1000/60;
% for g = 1:3
%     subplot(1,3,g);
%     hold on;plot(Xxcorr,res.CA1V1corrAll{g});
%     hold on;plot(Xxcorr,res.CA1V1corrAll_rand{g});
%     hold on;plot(Xxcorr,res.CA1V1corrAll_randX{g});
%     hold on;plot(Xxcorr,res.CA1V1corrAll_randS{g});
%     set(gca,'Xlim',[Xxcorr(1) Xxcorr(end)],'Ylim', [-0.02 0.07]);
%     title(titlestr{g});
% end


end

%to look at speed slices, run:
% nSpdbins = 3;
% meanErrAll = cell(2,3);
% popres_Rspd = [];
% for spd = 1:nSpdbins
% disp(spd)
% popres_Rspd{spd} = PopBayesAnalysis(res,false,false,spd,nSpdbins);
% n = 0;
% for ianimal = 1:10
% for iseries = 1:size(popres_Rspd{spd}.s_meanErrX_ave,2)
% n=n+1;
% for iprobe = 1:2
% for g = 1:3
% if ~isempty(popres_Rspd{spd}.s_meanErrX_ave{ianimal,iseries,iprobe,g})
% meanErrAll{iprobe,g}(spd,n) = popres_Rspd{spd}.s_meanErrX_ave{ianimal,iseries,iprobe,g};
% else
% meanErrAll{iprobe,g}(spd,n) = NaN;
% end
% end
% end
% end
% end
% end
% figure
% cl{1} = 'c';
% cl{3} = 'm';
% titlestr{1} = 'CA1';
% titlestr{2} = 'V1';
% for iprobe = 1:2
% for g=[1 3]
% for spd = 1:nSpdbins
% subplot(2,1,iprobe)
% hold on;bar(spd,nanmean(meanErrAll{iprobe,g}(spd,:)-meanErrAll{iprobe,2}(spd,:)),'EdgeColor',cl{g},'FaceColor',cl{g},'EdgeAlpha',0.3,'FaceAlpha',0.2);
% hold on;errorbar(spd,nanmean(meanErrAll{iprobe,g}(spd,:)-meanErrAll{iprobe,2}(spd,:)),nanstd(meanErrAll{iprobe,g}(spd,:)-meanErrAll{iprobe,2}(spd,:))/sqrt(sum(~isnan(meanErrAll{iprobe,g}(spd,:)-meanErrAll{iprobe,2}(spd,:)))),'Color',cl{g});
% set(gca,'Xlim',[0.5 5.5],'PlotBoxAspectRatio', [1 1 1])
% title(titlestr{iprobe});
% end
% end
% end

%to look at phase precession, run the following
% nbphsbins = 6;
% for iphs = 1:nbphsbins
%     popres_theta{iphs} = PopBayesAnalysis(res,[],[],[],iphs,nbphsbins);
% end
% 
% meanErrAll = [];
% n = 0;
% for ianimal = 1:10
% for iseries = 1:size(popres.s_meanErrX_ave,2)
% n=n+1;
% for iprobe = 1:2
% for g = 1:3
% for iphs = 1:6
% if ~isempty(popres_theta.s_meanErrX_ave{ianimal,iseries,iprobe,g})
% meanErrAll{iprobe,g}(iphs,n) = popres_theta{iphs}.s_meanErrX_ave{ianimal,iseries,iprobe,g}-50;
% else
% meanErrAll{iprobe,g}(iphs,n) = NaN;
% end
% end
% end
% [~,imax] = max(meanErrAll{1,2}(:,n));
% if ~isempty(imax)
% for g = 1:3
% meanErrAll{iprobe,g}(:,n) = circshift(meanErrAll{iprobe,g}(:,n),-imax);
% end
% end
% end
% end
% end
% figure
% cl{1} = 'c';
% cl{2} = 'k';
% cl{3} = 'm';
% titlestr{1} = 'CA1';
% titlestr{2} = 'V1';
% for iprobe = 1:2
% for iphs = 1:6
% for g=[1 3]
% subplot(2,2,iprobe+2)
% hold on;bar(iphs,nanmean(meanErrAll{iprobe,g}(iphs,:)-meanErrAll{iprobe,2}(iphs,:)),'EdgeColor',cl{g},'FaceColor',cl{g},'EdgeAlpha',0.3,'FaceAlpha',0.2);
% hold on;errorbar(iphs,nanmean(meanErrAll{iprobe,g}(iphs,:)-meanErrAll{iprobe,2}(iphs,:)),nanstd(meanErrAll{iprobe,g}(iphs,:)-meanErrAll{iprobe,2}(iphs,:))/sqrt(sum(~isnan(meanErrAll{iprobe,g}(iphs,:)-meanErrAll{iprobe,2}(iphs,:)))),'Color',cl{g});
% set(gca,'Xlim',[0.5 6.5],'PlotBoxAspectRatio', [1 1 1])
% xlabel('theta phs bin')
% ylabel('decoding error')
% title([titlestr{iprobe} ' lo/hi gain']);
% end
% subplot(2,2,iprobe)
% hold on;bar(iphs,nanmean(meanErrAll{iprobe,2}(iphs,:)),'EdgeColor',cl{2},'FaceColor',cl{2},'EdgeAlpha',0.3,'FaceAlpha',0.2);
% hold on;errorbar(iphs,nanmean(meanErrAll{iprobe,2}(iphs,:)),nanstd(meanErrAll{iprobe,2}(iphs,:))/sqrt(sum(~isnan(meanErrAll{iprobe,2}(iphs,:)))),'Color',cl{2});
% set(gca,'Xlim',[0.5 6.5],'PlotBoxAspectRatio', [1 1 1])
% xlabel('theta phs bin')
% ylabel('decoding error')
% title([titlestr{iprobe} ' med gain']);
% end
% end

function predave = getCircularAverageXdec(mat,maxtol)
    Prange = size(mat,1);
    Xrange = size(mat,2);
    postfilt = zeros(Prange,1);
    postfilt(1:Prange) = (0:(Prange-1))'/(Prange)*(2*pi)-pi/2;

    [maxval, ~] = max(mat,[],1);
    predmax = zeros(size(maxval));
    for i = 1:Xrange
        mat(mat(:,i) < (1-maxtol)*maxval(i),i) = 0;
    end

    a = sum((mat'*cos(postfilt)),2);b = sum((mat'*sin(postfilt)),2);
    predave = atan2(-a,b)+pi;%./sum(exp(log(2)*Post(tidx,:)),2);
    predave = predave/(2*pi)*Prange+1;
    predplus = predave + Prange;
    predminus = predave - Prange;
    xx = (1:Xrange)';
    predave = predave.* (abs(xx-predave)<abs(xx-predplus) & abs(xx-predave)<abs(xx-predminus)) + predplus.* (abs(xx-predplus)<abs(xx-predave) & abs(xx-predplus)<abs(xx-predminus)) + predminus.* (abs(xx-predminus)<abs(xx-predave) & abs(xx-predminus)<abs(xx-predplus));
    predave = smooth(predave,3);
end

function res = getresSEM(res,maxtol)
nanimal = size(res.PostX,1);
for iprobe = 1:2
    for g = 1:3
        for iperm = 0:nanimal
            PostXperm = 0;
            PostErrXperm = 0;
            DistriXperm = 0;
            DistriErrXperm = 0;
            nAll = 0;
            for ianimal = 1:nanimal%[1 2 3 4 5 8 9]%
                if sum(res.PostX{ianimal,iprobe,g}(:))>0 && ianimal ~= iperm
                    nAll = nAll + 1;
                    PostXperm = PostXperm + res.PostX{ianimal,iprobe,g};
                    PostErrXperm = PostErrXperm + res.PostErrX{ianimal,iprobe,g};
                    DistriXperm = DistriXperm + res.DistriX{ianimal,iprobe,g};
                    DistriErrXperm = DistriErrXperm + res.DistriErrX{ianimal,iprobe,g};
                end
            end
            PostXperm = PostXperm/nAll;
            Xpredperm = getCircularAverage(PostXperm,maxtol);
            XErr = Xpredperm - (1:100)';
            XErr(XErr>50) = XErr(XErr>50) - 100;
            XErr(XErr<-50) = XErr(XErr<-50) + 100;
            
            PostErrXperm = PostErrXperm/nAll;
            MeanXErrperm = nanmean(PostErrXperm,2);
            DistriXperm = DistriXperm/nAll;
            DistriErrXperm = DistriErrXperm/nAll;
            
            if iperm == 0
                res.PostXAll{iprobe,g} = PostXperm;%res.PostXAll{iprobe,g} = PostXperm;
                res.PostErrXAll{iprobe,g} = PostErrXperm;
                res.DistriXAll{iprobe,g} = DistriXperm;
                res.DistriErrXAll{iprobe,g} = DistriErrXperm;
                
                res.XpredAll{iprobe,g} = getCircularAverage(res.PostXAll{iprobe,g},maxtol);
                XErrAll{iprobe,g} = XErr;
                res.XpredAll{iprobe,g} = res.XpredAll{iprobe,g};
                res.MeanXErrAll{iprobe,g} = nanmean(res.PostErrXAll{iprobe,g},2);
                
                res.XpredAllsem{iprobe,g} = 0;
                res.MeanXErrAllsem{iprobe,g} = 0;
            else
                res.XpredAllsem{iprobe,g} = res.XpredAllsem{iprobe,g} + (XErr - XErrAll{iprobe,g}).^2;
                res.MeanXErrAllsem{iprobe,g} = res.MeanXErrAllsem{iprobe,g} + (MeanXErrperm - res.MeanXErrAll{iprobe,g}).^2;
            end
        end
        res.XpredAllsem{iprobe,g} = ((res.XpredAllsem{iprobe,g}*(nAll - 1)/nAll).^0.5)/nAll^0.5;
        res.MeanXErrAllsem{iprobe,g} = ((res.MeanXErrAllsem{iprobe,g}*(nAll - 1)/nAll).^0.5)/nAll^0.5;
    end
end


end

function mat_out = smooth2D(mat,lambdaSmooth)
mat(isnan(mat)) = 0;
G = smooth1D(repmat(mat,3,3),lambdaSmooth);
H = smooth1D(G',lambdaSmooth)';
mat_out = H(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
end

function [predave_out,x] = Xpredcorrection(predave,predave_ref,Xrange,dx)
predave_ref(predave_ref - (1:Xrange)' > floor(Xrange/2)) = predave_ref(predave_ref - (1:Xrange)' > floor(Xrange/2)) - Xrange;
predave_ref(predave_ref - (1:Xrange)' < -floor(Xrange/2)) = predave_ref(predave_ref - (1:Xrange)' < -floor(Xrange/2)) + Xrange;
predave(predave - (1:Xrange)' > floor(Xrange/2)) = predave(predave - (1:Xrange)' > floor(Xrange/2)) - Xrange;
predave(predave - (1:Xrange)' < -floor(Xrange/2)) = predave(predave - (1:Xrange)' < -floor(Xrange/2)) + Xrange;

predave_interp = interp1(linspace(0,numel(predave),numel(predave)+1), [predave(1);predave], 0:dx:(Xrange-dx));predave_interp(isnan(predave_interp)) = 0;
predaveref_interp = interp1(linspace(0,numel(predave_ref),numel(predave_ref)+1), [predave_ref(1);predave_ref], 0:dx:(Xrange-dx));predaveref_interp(isnan(predaveref_interp)) = 0;
predave_out = zeros(size(predave_interp));
x = 0:dx:(Xrange-dx);
predaveref_interp = [predaveref_interp(1:end-1)-Xrange predaveref_interp predaveref_interp(2:end)+Xrange];
xrep = [x(1:end-1)-Xrange x x(2:end)+Xrange];
for i = 1:numel(x)
    try
    idxmatch = min(find(abs(predaveref_interp-predave_interp(i)) <= min(abs(predaveref_interp-predave_interp(i)))));
    idxmatch = idxmatch(abs(idxmatch-(numel(x)-1) - i) == min(abs(idxmatch-(numel(x)-1) - i)));
    predave_out(i) = xrep(round(idxmatch)) - x(i);
    catch
        keyboard
    end
end
% predave_out(predave_out>50) = predave_out(predave_out>50) - 100;
% predave_out(predave_out<-50) = predave_out(predave_out<-50) + 100;
end