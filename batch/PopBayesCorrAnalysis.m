function  popres = PopBayesCorrAnalysis(res,popres,batch2p,FshowsessionFig,spdbins,nSpdbins,phsbins,nPhsbins)
if nargin < 3
    batch2p = [];
end
if nargin < 4
    FshowsessionFig = [];
end
if nargin < 5
    nSpdbins = 5;
    spdbins = 1:5;
end
if isempty(nSpdbins)
    nSpdbins = 5;
    spdbins = 1:5;
end
if nargin < 7
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

lambdaSmooth = 2;
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

maxtol_post = 1;maxtol_meanpost = 0.1;%0.1;
maxtol_decmax = 1;maxtol_meandecmax = 0.1;
maxtol_decave = 1;maxtol_meandecave = 0.1;%0.1;%1;%
amp_th_post = 1;amp_th_meanpost = 1;%0;
amp_th_decmax = 1;amp_th_meandecmax = 1;%0;
amp_th_decave = 1;amp_th_meandecave = 1;%0;
popres.dx = 0.1;


nanimal = numel(expt);
ngain = 3;

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

Xrange = 100;
popres.CA1V1corr = cell(3,Xrange);
popres.CA1V1corr_rand = cell(3,Xrange);
popres.CA1V1corr_randX = cell(3,Xrange);
popres.CA1V1corr_randS = cell(3,Xrange);
popres.CA1V1corr_randXS = cell(3,Xrange);

popres.CA1V12Dcorr = cell(3,Xrange);
popres.CA1V12Dcorr_rand = cell(3,Xrange);
popres.CA1V12Dcorr_randX = cell(3,Xrange);
popres.CA1V12Dcorr_randS = cell(3,Xrange);
popres.CA1V12Dcorr_randXS = cell(3,Xrange);

popres.CA1V12Dcorr_randsqr = cell(3,Xrange);
popres.CA1V12Dcorr_randXsqr = cell(3,Xrange);
popres.CA1V12Dcorr_randSsqr = cell(3,Xrange);
popres.CA1V12Dcorr_randXSsqr = cell(3,Xrange);

popres.CA1V11Dcorr = cell(3,Xrange);
popres.CA1V11Dcorr_rand = cell(3,Xrange);
popres.CA1V11Dcorr_randX = cell(3,Xrange);
popres.CA1V11Dcorr_randS = cell(3,Xrange);
popres.CA1V11Dcorr_randXS = cell(3,Xrange);
Nsession = cell(3,Xrange);

for g = 1:3
    for xx = 1:Xrange
    popres.CA1V1corr_rand{g,xx} = 0;
    popres.CA1V1corr_randX{g,xx} = 0;
    popres.CA1V1corr_randS{g,xx} = 0;
    popres.CA1V1corr_randXS{g,xx} = 0;
    
    popres.CA1V1corr_randsqr{g,xx} = 0;
    popres.CA1V1corr_randXsqr{g,xx} = 0;
    popres.CA1V1corr_randSsqr{g,xx} = 0;
    popres.CA1V1corr_randXSsqr{g,xx} = 0;

    popres.CA1V12Dcorr{g,xx} = 0;
    popres.CA1V12Dcorr_rand{g,xx} = 0;
    popres.CA1V12Dcorr_randX{g,xx} = 0;
    popres.CA1V12Dcorr_randS{g,xx} = 0;
    popres.CA1V12Dcorr_randXS{g,xx} = 0;
    
    popres.CA1V12Dcorr_randsqr{g,xx} = 0;
    popres.CA1V12Dcorr_randXsqr{g,xx} = 0;
    popres.CA1V12Dcorr_randSsqr{g,xx} = 0;
    popres.CA1V12Dcorr_randXSsqr{g,xx} = 0;
    
    popres.CA1V11Dcorr{g,xx} = 0;
    popres.CA1V11Dcorr_rand{g,xx} = 0;
    popres.CA1V11Dcorr_randX{g,xx} = 0;
    popres.CA1V11Dcorr_randS{g,xx} = 0;
    popres.CA1V11Dcorr_randXS{g,xx} = 0;
    
    popres.CA1V11Dcorr_randsqr{g,xx} = 0;
    popres.CA1V11Dcorr_randXsqr{g,xx} = 0;
    popres.CA1V11Dcorr_randSsqr{g,xx} = 0;
    popres.CA1V11Dcorr_randXSsqr{g,xx} = 0;
    
    popres.CA1V1mscoh{g,xx} = 0;
    popres.CA1V1mscoh_rand{g,xx} = 0;
    popres.CA1V1mscoh_randX{g,xx} = 0;
    popres.CA1V1mscoh_randS{g,xx} = 0;
    popres.CA1V1mscoh_randXS{g,xx} = 0;
    Nsession{g,xx} = 0;
    end
end
outvalcorr = 2;%[0 1 2 3 4 5];%[0 1 2 3 4];%
iprobe = 1;
for ianimal = 1:nanimal
    for iseries = 1:numel(expt(ianimal).series)
        if (expt(ianimal).goodCA1{iseries} == 1) && (expt(ianimal).goodV1{iseries} == 1)
            Xrange = max(res.X{ianimal,iseries}(~isnan(res.X{ianimal,iseries})));
            
            speeds = NaN(size(res.ballspeed{ianimal,iseries}));
            speeds(~isnan(res.ballspeed{ianimal,iseries})) = smthInTime(res.ballspeed{ianimal,iseries}(~isnan(res.ballspeed{ianimal,iseries})), sampleRate, SpdSmthWin, 'same', [], 'boxcar_centered');
            
            XErrCA1 = res.Xpred_ave{ianimal,iseries,1} - res.X{ianimal,iseries};%res.Xpred_max{ianimal,iseries,1} - res.X{ianimal,iseries};
            XErrCA1(XErrCA1>floor(Xrange/2)) = XErrCA1(XErrCA1>floor(Xrange/2)) - Xrange;
            XErrCA1(XErrCA1<-floor(Xrange/2)) = XErrCA1(XErrCA1<-floor(Xrange/2)) + Xrange;
            
            XErrV1 = res.Xpred_ave{ianimal,iseries,2} - res.X{ianimal,iseries};%res.Xpred_max{ianimal,iseries,2} - res.X{ianimal,iseries};
            XErrV1(XErrV1>floor(Xrange/2)) = XErrV1(XErrV1>floor(Xrange/2)) - Xrange;
            XErrV1(XErrV1<-floor(Xrange/2)) = XErrV1(XErrV1<-floor(Xrange/2)) + Xrange;
            
            XErrCA10 = XErrCA1;
            XErrV10 = XErrV1;
            for g = [2 1 3]
                cont_list = find(ismember(res.contrastVal{ianimal,iseries}, contval));
                RL_list = find(ismember(res.roomlengthVal{ianimal,iseries}, [1]));
                outcome_list = find(ismember(res.outcomeVal{ianimal,iseries}, outvalcorr));
                tidx = false(size(res.tidx{ianimal,iseries,1,1,2,1,3}));
                for cont = 1:numel(cont_list)
                    for r = 1:numel(RL_list)
                        for o = 1:numel(outcome_list)
                            if ~isempty(res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)})
                                tidx = tidx | (res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)} & ismember(res.Phsbin{ianimal,iseries},phsbins));% & res.X{ianimal,iseries} > 0 & res.X{ianimal,iseries} <= 10);
                            end
                        end
                    end
                end
                XErrCA1 = XErrCA10;
                XErrV1 = XErrV10;
                tidx0 = tidx & ismember(res.Spdbin2{ianimal,iseries},spdbins) & ~isnan(XErrCA10) & ~isnan(XErrV10);
%                 [Fave_CA1, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(~isnan(XErrCA1) & tidx0) XErrCA1(~isnan(XErrCA1) & tidx0)], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
%                 [Fave_V1, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(~isnan(XErrV1) & tidx0) XErrV1(~isnan(XErrV1) & tidx0)], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
%                 Fave_CA1err = getCircularAverage(Fave_CA1,amp_th_meandecave,maxtol_meandecave) - floor(Xrange/2);
%                 Fave_V1err = getCircularAverage(Fave_V1,amp_th_meandecave,maxtol_meandecave) - floor(Xrange/2);
%               
                
                for xx = 1:Xrange
                    xwin = mod(xx + (-1:1) -1,Xrange)+1;
                    tidx = tidx0 & ismember(res.X{ianimal,iseries},xwin);
                    XErrCA1(tidx) = XErrCA10(tidx);%- mean(Fave_CA1err);
                    XErrCA1(tidx & XErrCA1>floor(Xrange/2)) = XErrCA1(tidx & XErrCA1>floor(Xrange/2)) - Xrange;
                    XErrCA1(tidx & XErrCA1<-floor(Xrange/2)) = XErrCA1(tidx & XErrCA1<-floor(Xrange/2)) + Xrange;
                    XErrV1(tidx) = XErrV10(tidx);% - mean(Fave_V1err);
                    XErrV1(tidx & XErrV1>floor(Xrange/2)) = XErrV1(tidx & XErrV1>floor(Xrange/2)) - Xrange;
                    XErrV1(tidx & XErrV1<-floor(Xrange/2)) = XErrV1(tidx & XErrV1<-floor(Xrange/2)) + Xrange;
%                 end

%                 XErrCA1 = XErrCA10;
%                 XErrV1 = circshift(XErrV10,xx-floor(Xrange/2));
%                 tidx = tidx0 & ~isnan(XErrCA1) & ~isnan(XErrV1);
                
                            
                if sum(tidx)>0
                    Nsession{g,xx} = Nsession{g,xx} + 1;
                    idx = find(tidx);
%                     XErrV1CA1 = XErrV1(idx) - XErrCA1(idx);
%                     XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
%                     XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
%                     CA1V1noisecorr = sum((XErrCA1(idx)+XErrV1CA1 - nanmean((XErrCA1(idx)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx)+XErrV1CA1 - nanmean((XErrCA1(idx)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
% 
                    XerrV1ang = (XErrV1(idx)/Xrange)*2*pi;
                    a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
                    XerrV1ang_mean = atan2(b,a);
                    XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
                    a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
                    XerrCA1ang_mean = atan2(b,a);
%                     [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
                    CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;

%                     CA1V1noisecorr = sum(XErrV1(idx).*XErrCA1(idx))/(sum(XErrV1(idx).^2)*sum(XErrCA1(idx).^2))^0.5;
                    popres.CA1V1corr{g,xx} = [popres.CA1V1corr{g,xx} CA1V1noisecorr];
                    
                    Xwincorr = (-180:180);
                    xcorrtemp = zeros(numel(Xwincorr),1);
                    for tt_corr = Xwincorr
                        XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
                        XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
                        xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
                    end
                    popres.CA1V11Dcorr{g,xx} = popres.CA1V11Dcorr{g,xx} + xcorrtemp;%xcorr(XErrCA1(idx),XErrV1(idx),3*60,'coeff');
                    %                     popres.CA1V1mscoh{g} = popres.CA1V1mscoh{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx),512,128));
                    
                    
%                     [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                    [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                    popres.CA1V12Dcorr{g,xx} = popres.CA1V12Dcorr{g,xx} + Fave*numel(idx);
                else
                    popres.CA1V1corr{g,xx} = [popres.CA1V1corr{g,xx} NaN];
                end
            end
            end
        end
    end
end
for g = 1:3
    for xx = 1:Xrange
    popres.CA1V12Dcorr{g,xx} = popres.CA1V12Dcorr{g,xx}/sum(popres.CA1V12Dcorr{g,xx}(:))/(1/numel(popres.CA1V12Dcorr{g,xx}(:)));%Nsession{g};
    popres.CA1V11Dcorr{g,xx} = popres.CA1V11Dcorr{g,xx}/Nsession{g,xx};
    end
end

Nranditer = 30;%
for kiter = 1:Nranditer
    tempCA1V1corr_rand = cell(3,Xrange);
    tempCA1V1corr_randX = cell(3,Xrange);
    tempCA1V1corr_randS = cell(3,Xrange);
    tempCA1V1corr_randXS = cell(3,Xrange);
    for g = 1:3
        for xx = 1:Xrange
        tempCA1V12Dcorr{g,xx} = 0;
        tempCA1V12Dcorr_rand{g,xx} = 0;
        tempCA1V12Dcorr_randX{g,xx} = 0;
        tempCA1V12Dcorr_randS{g,xx} = 0;
        tempCA1V12Dcorr_randXS{g,xx} = 0;
        
        tempCA1V11Dcorr{g,xx} = 0;
        tempCA1V11Dcorr_rand{g,xx} = 0;
        tempCA1V11Dcorr_randX{g,xx} = 0;
        tempCA1V11Dcorr_randS{g,xx} = 0;
        tempCA1V11Dcorr_randXS{g,xx} = 0;
        
        tempCA1V1mscoh{g,xx} = 0;
        tempCA1V1mscoh_rand{g,xx} = 0;
        tempCA1V1mscoh_randX{g,xx} = 0;
        tempCA1V1mscoh_randS{g,xx} = 0;
        tempCA1V1mscoh_randXS{g,xx} = 0;
        Nsession{g,xx} = 0;
        end
    end
    for ianimal = 1:nanimal
        for iseries = 1:numel(expt(ianimal).series)
            if (expt(ianimal).goodCA1{iseries} == 1) && (expt(ianimal).goodV1{iseries} == 1)
                
                speeds = NaN(size(res.ballspeed{ianimal,iseries}));
                speeds(~isnan(res.ballspeed{ianimal,iseries})) = smthInTime(res.ballspeed{ianimal,iseries}(~isnan(res.ballspeed{ianimal,iseries})), sampleRate, SpdSmthWin, 'same', [], 'boxcar_centered');
                
                XErrCA1 = res.Xpred_ave{ianimal,iseries,1} - res.X{ianimal,iseries};%res.Xpred_max{ianimal,iseries,1} - res.X{ianimal,iseries};
                XErrCA1(XErrCA1>floor(Xrange/2)) = XErrCA1(XErrCA1>floor(Xrange/2)) - Xrange;
                XErrCA1(XErrCA1<-floor(Xrange/2)) = XErrCA1(XErrCA1<-floor(Xrange/2)) + Xrange;
                
                XErrV1 = res.Xpred_ave{ianimal,iseries,2} - res.X{ianimal,iseries};%res.Xpred_max{ianimal,iseries,2} - res.X{ianimal,iseries};
                XErrV1(XErrV1>floor(Xrange/2)) = XErrV1(XErrV1>floor(Xrange/2)) - Xrange;
                XErrV1(XErrV1<-floor(Xrange/2)) = XErrV1(XErrV1<-floor(Xrange/2)) + Xrange;
                
                %             posErrCA1 = 2*randn(size(XErrCA1));
                %             posErrV1 = posErrCA1;
                %             posErrCA1 = posErrCA1 + randn(size(posErrCA1))*3;
                %             posErrV1 = posErrV1 + randn(size(posErrV1))*3;
                %
                %             XErrCA1 = posErrCA1 - nanmean(posErrCA1);
                %             XErrV1 = posErrV1 - nanmean(posErrV1);
                %             Xshift = [5 0 -5];
                
                XErrCA10 = XErrCA1;
                XErrV10 = XErrV1;
                for g = 1:3
                    cont_list = find(ismember(res.contrastVal{ianimal,iseries}, contval));
                    RL_list = find(ismember(res.roomlengthVal{ianimal,iseries}, [1]));
                    outcome_list = find(ismember(res.outcomeVal{ianimal,iseries}, outvalcorr));
                    tidx = false(size(res.tidx{ianimal,iseries,1,1,2,1,3}));
                    for cont = 1:numel(cont_list)
                        for r = 1:numel(RL_list)
                            for o = 1:numel(outcome_list)
                                if ~isempty(res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)})
                                    tidx = tidx | (res.tidx{ianimal,iseries,iprobe,cont_list(cont),g,RL_list(r),outcome_list(o)} & ismember(res.Phsbin{ianimal,iseries},phsbins));% & res.X{ianimal,iseries} > 0 & res.X{ianimal,iseries} <= 10);
                                end
                            end
                        end
                    end
                    XErrCA1 = XErrCA10;
                    XErrV1 = XErrV10;
                    tidx0 = tidx & ismember(res.Spdbin2{ianimal,iseries},spdbins) & ~isnan(XErrCA1) & ~isnan(XErrV1);
%                     [Fave_CA1, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(~isnan(XErrCA1) & tidx0) XErrCA1(~isnan(XErrCA1) & tidx0)], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
%                     [Fave_V1, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(~isnan(XErrV1) & tidx0) XErrV1(~isnan(XErrV1) & tidx0)], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
%                     Fave_CA1err = getCircularAverage(Fave_CA1,amp_th_meandecave,maxtol_meandecave) - floor(Xrange/2);
%                     Fave_V1err = getCircularAverage(Fave_V1,amp_th_meandecave,maxtol_meandecave) - floor(Xrange/2);
%                     
                    
                    for xx = 1:Xrange
                        xwin = mod(xx + (-1:1) -1,Xrange)+1;
                        tidx = tidx0 & ismember(res.X{ianimal,iseries},xwin);
                        XErrCA1(tidx) = XErrCA10(tidx);% - mean(Fave_CA1err);
                        XErrCA1(tidx & XErrCA1>floor(Xrange/2)) = XErrCA1(tidx & XErrCA1>floor(Xrange/2)) - Xrange;
                        XErrCA1(tidx & XErrCA1<-floor(Xrange/2)) = XErrCA1(tidx & XErrCA1<-floor(Xrange/2)) + Xrange;
                        XErrV1(tidx) = XErrV10(tidx);% - mean(Fave_V1err);
                        XErrV1(tidx & XErrV1>floor(Xrange/2)) = XErrV1(tidx & XErrV1>floor(Xrange/2)) - Xrange;
                        XErrV1(tidx & XErrV1<-floor(Xrange/2)) = XErrV1(tidx & XErrV1<-floor(Xrange/2)) + Xrange;
%                     end

%                     XErrCA1 = XErrCA10;
%                     XErrV1 = circshift(XErrV10,xx-floor(Xrange/2));
%                     tidx = tidx0 & ~isnan(XErrCA1) & ~isnan(XErrV1);
                    
                    if sum(tidx)>0
                        Nsession{g,xx} = Nsession{g,xx} + 1;
                        idx = find(tidx);
                        idx_rand = idx(randperm(numel(idx)));
%                         XErrV1CA1 = XErrV1(idx_rand) - XErrCA1(idx_rand);
%                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
%                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
%                         CA1V1noisecorr = sum((XErrCA1(idx_rand)+XErrV1CA1 - nanmean((XErrCA1(idx_rand)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx_rand)+XErrV1CA1 - nanmean((XErrCA1(idx_rand)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
                        
                        XerrV1ang = (XErrV1(idx_rand)/Xrange)*2*pi;
                        a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
                        XerrV1ang_mean = atan2(b,a);
                        XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
                        a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
                        XerrCA1ang_mean = atan2(b,a);
%                         [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
                        CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;
                        
%                         CA1V1noisecorr = sum(XErrV1(idx_rand).*XErrCA1(idx))/(sum(XErrV1(idx_rand).^2)*sum(XErrCA1(idx).^2))^0.5;

                        Xwincorr = (-180:180);
                        xcorrtemp = zeros(numel(Xwincorr),1);
                        for tt_corr = Xwincorr
                            XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
                            XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
                            xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
                        end
                        tempCA1V11Dcorr_rand{g,xx} = tempCA1V11Dcorr_rand{g,xx} + xcorrtemp;%tempCA1V11Dcorr_rand{g} = tempCA1V11Dcorr_rand{g} + xcorr(XErrCA1(idx),XErrV1(idx_rand),3*60,'coeff');
                        %                         popres.CA1V1mscoh_rand{g} = popres.CA1V1mscoh_rand{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx_rand),512,128));
                        
%                         XErrV1CA1 = res.Xpred_ave{ianimal,iseries,2}(idx_rand) - res.Xpred_ave{ianimal,iseries,1}(idx);
%                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
%                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
%                         [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                        [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx_rand)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                        tempCA1V12Dcorr_rand{g,xx} = tempCA1V12Dcorr_rand{g,xx} + Fave*numel(idx);
                        tempCA1V1corr_rand{g,xx} = [tempCA1V1corr_rand{g,xx} CA1V1noisecorr];
                        
                        
                        idx = find(tidx);
                        idx_randX = zeros(size(idx));
                        for ix = 1:max(res.X{ianimal,iseries})
                            idxtemp = idx(res.X{ianimal,iseries}(idx) == ix);
                            idx_randX(res.X{ianimal,iseries}(idx) == ix) = idxtemp(randperm(numel(idxtemp)));
                        end
%                         XErrV1CA1 = XErrV1(idx_randX) - XErrCA1(idx_randX);
%                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
%                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
%                         CA1V1noisecorr = sum((XErrCA1(idx_randX)+XErrV1CA1 - nanmean((XErrCA1(idx_randX)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx_randX)+XErrV1CA1 - nanmean((XErrCA1(idx_randX)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
                        
                        XerrV1ang = (XErrV1(idx_randX)/Xrange)*2*pi;
                        a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
                        XerrV1ang_mean = atan2(b,a);
                        XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
                        a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
                        XerrCA1ang_mean = atan2(b,a);
%                         [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
                        CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;
                        
%                         CA1V1noisecorr = sum(XErrV1(idx_randX).*XErrCA1(idx))/(sum(XErrV1(idx_randX).^2)*sum(XErrCA1(idx).^2))^0.5;
                        Xwincorr = (-180:180);
                        xcorrtemp = zeros(numel(Xwincorr),1);
                        for tt_corr = Xwincorr
                            XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
                            XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
                            xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
                        end
                        tempCA1V11Dcorr_randX{g,xx} = tempCA1V11Dcorr_randX{g,xx} + xcorrtemp;%tempCA1V11Dcorr_randX{g} = tempCA1V11Dcorr_randX{g} + xcorr(XErrCA1(idx),XErrV1(idx_randX),3*60,'coeff');
                        %                         popres.CA1V1mscoh_randX{g} = popres.CA1V1mscoh_randX{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx_randX),512,128));
                        
%                         XErrV1CA1 = res.Xpred_ave{ianimal,iseries,2}(idx_randX) - res.Xpred_ave{ianimal,iseries,1}(idx);
%                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
%                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
%                         [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                        [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx_randX)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                        tempCA1V12Dcorr_randX{g,xx} = tempCA1V12Dcorr_randX{g,xx} + Fave*numel(idx);
                        tempCA1V1corr_randX{g,xx} = [tempCA1V1corr_randX{g,xx} CA1V1noisecorr];
                        
                        idx = find(tidx);
                        idx_randSpeed = zeros(size(idx));
                        for ispeed = 1:nSpdbins
                            idxtemp = idx(res.Spdbin2{ianimal,iseries}(idx) == ispeed);
                            idx_randSpeed(res.Spdbin2{ianimal,iseries}(idx) == ispeed) = idxtemp(randperm(numel(idxtemp)));
                        end
%                         XErrV1CA1 = XErrV1(idx_randSpeed) - XErrCA1(idx_randSpeed);
%                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
%                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
%                         CA1V1noisecorr = sum((XErrCA1(idx_randSpeed)+XErrV1CA1 - nanmean((XErrCA1(idx_randSpeed)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx_randSpeed)+XErrV1CA1 - nanmean((XErrCA1(idx_randSpeed)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
                        
                        XerrV1ang = (XErrV1(idx_randSpeed)/Xrange)*2*pi;
                        a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
                        XerrV1ang_mean = atan2(b,a);
                        XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
                        a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
                        XerrCA1ang_mean = atan2(b,a);
%                         [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
                        CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;
                        
%                         CA1V1noisecorr = sum(XErrV1(idx_randSpeed).*XErrCA1(idx))/(sum(XErrV1(idx_randSpeed).^2)*sum(XErrCA1(idx).^2))^0.5;
                        Xwincorr = (-180:180);
                        xcorrtemp = zeros(numel(Xwincorr),1);
                        for tt_corr = Xwincorr
                            XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
                            XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
                            xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
                        end
                        tempCA1V11Dcorr_randS{g,xx} = tempCA1V11Dcorr_randS{g,xx} + xcorrtemp;%tempCA1V11Dcorr_randS{g} = tempCA1V11Dcorr_randS{g} + xcorr(XErrCA1(idx),XErrV1(idx_randSpeed),3*60,'coeff');
                        %                         popres.CA1V1mscoh_randS{g} = popres.CA1V1mscoh_randS{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx_randSpeed),512,128));
                        
%                         XErrV1CA1 = res.Xpred_ave{ianimal,iseries,2}(idx_randSpeed) - res.Xpred_ave{ianimal,iseries,1}(idx);
%                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
%                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
%                         [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                        [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx_randSpeed)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                        tempCA1V12Dcorr_randS{g,xx} = tempCA1V12Dcorr_randS{g,xx} + Fave*numel(idx);
                        tempCA1V1corr_randS{g,xx} = [tempCA1V1corr_randS{g,xx} CA1V1noisecorr];
                        
                        idx = find(tidx);
                        idx_randXS = zeros(size(idx));
                        for ispeed = 1:nSpdbins
                            for ix = min(res.X{ianimal,iseries}(idx)):max(res.X{ianimal,iseries}(idx))
                                idxtemp = idx(res.X{ianimal,iseries}(idx) == ix & res.Spdbin2{ianimal,iseries}(idx) == ispeed);
                                idx_randXS(res.X{ianimal,iseries}(idx) == ix & res.Spdbin2{ianimal,iseries}(idx) == ispeed) = idxtemp(randperm(numel(idxtemp)));
                            end
                        end
%                         XErrV1CA1 = XErrV1(idx_randXS) - XErrCA1(idx_randXS);
%                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
%                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
%                         CA1V1noisecorr = sum((XErrCA1(idx_randXS)+XErrV1CA1 - nanmean((XErrCA1(idx_randXS)+XErrV1CA1))).*(XErrCA1(idx)-nanmean(XErrCA1(idx))))/(sum((XErrCA1(idx_randXS)+XErrV1CA1 - nanmean((XErrCA1(idx_randXS)+XErrV1CA1))).^2)*sum((XErrCA1(idx)-nanmean(XErrCA1(idx))).^2))^0.5;
                        
                        XerrV1ang = (XErrV1(idx_randXS)/Xrange)*2*pi;
                        a = sum(cos(XerrV1ang));b = sum(sin(XerrV1ang));
                        XerrV1ang_mean = atan2(b,a);
                        XerrCA1ang = (XErrCA1(idx)/Xrange)*2*pi;
                        a = sum(cos(XerrCA1ang));b = sum(sin(XerrCA1ang));
                        XerrCA1ang_mean = atan2(b,a);
%                         [CA1V1noisecorr, ~] = circ_corrcc(XerrV1ang, XerrCA1ang);
                        CA1V1noisecorr = sum(sin(XerrV1ang - XerrV1ang_mean).*sin(XerrCA1ang - XerrCA1ang_mean))/(sum(sin(XerrV1ang - XerrV1ang_mean).^2)*sum(sin(XerrCA1ang - XerrCA1ang_mean).^2))^0.5;
                        
%                         CA1V1noisecorr = sum(XErrV1(idx_randXS).*XErrCA1(idx))/(sum(XErrV1(idx_randXS).^2)*sum(XErrCA1(idx).^2))^0.5;
                        Xwincorr = (-180:180);
                        xcorrtemp = zeros(numel(Xwincorr),1);
                        for tt_corr = Xwincorr
                            XerrV1angtemp = circshift(XerrV1ang - XerrV1ang_mean, -tt_corr);
                            XerrCA1angtemp = XerrCA1ang - XerrCA1ang_mean;
                            xcorrtemp(tt_corr - Xwincorr(1) + 1) = sum(sin(XerrV1angtemp).*sin(XerrCA1angtemp))/(sum(sin(XerrV1angtemp).^2)*sum(sin(XerrCA1angtemp).^2))^0.5;
                        end
                        tempCA1V11Dcorr_randXS{g,xx} = tempCA1V11Dcorr_randXS{g,xx} + xcorrtemp;%tempCA1V11Dcorr_randXS{g} = tempCA1V11Dcorr_randXS{g} + xcorr(XErrCA1(idx),XErrV1(idx_randXS),3*60,'coeff');
                        %                         popres.CA1V1mscoh_randXS{g} = popres.CA1V1mscoh_randXS{g} + abs(mscohere(XErrCA1(idx),XErrV1(idx_randXS),512,128));
                        
%                         XErrV1CA1 = res.Xpred_ave{ianimal,iseries,2}(idx_randXS) - res.Xpred_ave{ianimal,iseries,1}(idx);
%                         XErrV1CA1(XErrV1CA1>floor(Xrange/2)) = XErrV1CA1(XErrV1CA1>floor(Xrange/2)) - Xrange;
%                         XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) = XErrV1CA1(XErrV1CA1<-floor(Xrange/2)) + Xrange;
%                         [Fave, ~, ~, ~] = smoothhist2D_corrected([res.X{ianimal,iseries}(idx) XErrV1CA1], lambdaSmooth, [Xrange Xrange], 1:Xrange, (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                        [Fave, ~, ~, ~] = smoothhist2D_corrected2([XErrCA1(idx) XErrV1(idx_randXS)], lambdaSmooth, [Xrange Xrange], (-floor(Xrange/2)+1):floor(Xrange/2), (-floor(Xrange/2)+1):floor(Xrange/2), true, true);
                        tempCA1V12Dcorr_randXS{g,xx} = tempCA1V12Dcorr_randXS{g,xx} + Fave*numel(idx);
                        tempCA1V1corr_randXS{g,xx} = [tempCA1V1corr_randXS{g,xx} CA1V1noisecorr];
                    else
                        tempCA1V1corr_rand{g,xx} = [tempCA1V1corr_rand{g,xx} NaN];
                        tempCA1V1corr_randX{g,xx} = [tempCA1V1corr_randX{g,xx} NaN];
                        tempCA1V1corr_randS{g,xx} = [tempCA1V1corr_randS{g,xx} NaN];
                        tempCA1V1corr_randXS{g,xx} = [tempCA1V1corr_randXS{g,xx} NaN];
                    end
                end
                end
            end
        end
    end
    for g = 1:3
        for xx = 1:Xrange
        tempCA1V12Dcorr_rand{g,xx} = tempCA1V12Dcorr_rand{g,xx}/sum(tempCA1V12Dcorr_rand{g,xx}(:))/(1/numel(tempCA1V12Dcorr_rand{g,xx}(:)));%/Nsession{g,xx};
        tempCA1V12Dcorr_randX{g,xx} = tempCA1V12Dcorr_randX{g,xx}/sum(tempCA1V12Dcorr_randX{g,xx}(:))/(1/numel(tempCA1V12Dcorr_randX{g,xx}(:)));%/Nsession{g,xx};
        tempCA1V12Dcorr_randS{g,xx} = tempCA1V12Dcorr_randS{g,xx}/sum(tempCA1V12Dcorr_randS{g,xx}(:))/(1/numel(tempCA1V12Dcorr_randS{g,xx}(:)));%/Nsession{g,xx};
        tempCA1V12Dcorr_randXS{g,xx} = tempCA1V12Dcorr_randXS{g,xx}/sum(tempCA1V12Dcorr_randXS{g,xx}(:))/(1/numel(tempCA1V12Dcorr_randXS{g,xx}(:)));%/Nsession{g,xx};
        
        tempCA1V11Dcorr_rand{g,xx} = tempCA1V11Dcorr_rand{g,xx}/Nsession{g,xx};
        tempCA1V11Dcorr_randX{g,xx} = tempCA1V11Dcorr_randX{g,xx}/Nsession{g,xx};
        tempCA1V11Dcorr_randS{g,xx} = tempCA1V11Dcorr_randS{g,xx}/Nsession{g,xx};
        tempCA1V11Dcorr_randXS{g,xx} = tempCA1V11Dcorr_randXS{g,xx}/Nsession{g,xx};
        
        popres.CA1V12Dcorr_rand{g,xx} = popres.CA1V12Dcorr_rand{g,xx} + tempCA1V12Dcorr_rand{g,xx}/Nranditer;
        popres.CA1V12Dcorr_randX{g,xx} = popres.CA1V12Dcorr_randX{g,xx} + tempCA1V12Dcorr_randX{g,xx}/Nranditer;
        popres.CA1V12Dcorr_randS{g,xx} = popres.CA1V12Dcorr_randS{g,xx} + tempCA1V12Dcorr_randS{g,xx}/Nranditer;
        popres.CA1V12Dcorr_randXS{g,xx} = popres.CA1V12Dcorr_randXS{g,xx} + tempCA1V12Dcorr_randXS{g,xx}/Nranditer;
        
        popres.CA1V12Dcorr_randsqr{g,xx} = popres.CA1V12Dcorr_randsqr{g,xx} + tempCA1V12Dcorr_rand{g,xx}.^2/Nranditer;
        popres.CA1V12Dcorr_randXsqr{g,xx} = popres.CA1V12Dcorr_randXsqr{g,xx} + tempCA1V12Dcorr_randX{g,xx}.^2/Nranditer;
        popres.CA1V12Dcorr_randSsqr{g,xx} = popres.CA1V12Dcorr_randSsqr{g,xx} + tempCA1V12Dcorr_randS{g,xx}.^2/Nranditer;
        popres.CA1V12Dcorr_randXSsqr{g,xx} = popres.CA1V12Dcorr_randXSsqr{g,xx} + tempCA1V12Dcorr_randXS{g,xx}.^2/Nranditer;
        
        popres.CA1V11Dcorr_rand{g,xx} = popres.CA1V11Dcorr_rand{g,xx} + tempCA1V11Dcorr_rand{g,xx}/Nranditer;
        popres.CA1V11Dcorr_randX{g,xx} = popres.CA1V11Dcorr_randX{g,xx} + tempCA1V11Dcorr_randX{g,xx}/Nranditer;
        popres.CA1V11Dcorr_randS{g,xx} = popres.CA1V11Dcorr_randS{g,xx} + tempCA1V11Dcorr_randS{g,xx}/Nranditer;
        popres.CA1V11Dcorr_randXS{g,xx} = popres.CA1V11Dcorr_randXS{g,xx} + tempCA1V11Dcorr_randXS{g,xx}/Nranditer;
        
        popres.CA1V11Dcorr_randsqr{g,xx} = popres.CA1V11Dcorr_randsqr{g,xx} + tempCA1V11Dcorr_rand{g,xx}.^2/Nranditer;
        popres.CA1V11Dcorr_randXsqr{g,xx} = popres.CA1V11Dcorr_randXsqr{g,xx} + tempCA1V11Dcorr_randX{g,xx}.^2/Nranditer;
        popres.CA1V11Dcorr_randSsqr{g,xx} = popres.CA1V11Dcorr_randSsqr{g,xx} + tempCA1V11Dcorr_randS{g,xx}.^2/Nranditer;
        popres.CA1V11Dcorr_randXSsqr{g,xx} = popres.CA1V11Dcorr_randXSsqr{g,xx} + tempCA1V11Dcorr_randXS{g,xx}.^2/Nranditer;
        
        popres.CA1V1corr_rand{g,xx} = popres.CA1V1corr_rand{g,xx} + tempCA1V1corr_rand{g,xx}/Nranditer;
        popres.CA1V1corr_randX{g,xx} = popres.CA1V1corr_randX{g,xx} + tempCA1V1corr_randX{g,xx}/Nranditer;
        popres.CA1V1corr_randS{g,xx} = popres.CA1V1corr_randS{g,xx} + tempCA1V1corr_randS{g,xx}/Nranditer;
        popres.CA1V1corr_randXS{g,xx} = popres.CA1V1corr_randXS{g,xx} + tempCA1V1corr_randXS{g,xx}/Nranditer;
        
        popres.CA1V1corr_randsqr{g,xx} = popres.CA1V1corr_randsqr{g,xx} + tempCA1V1corr_rand{g,xx}.^2/Nranditer;
        popres.CA1V1corr_randXsqr{g,xx} = popres.CA1V1corr_randXsqr{g,xx} + tempCA1V1corr_randX{g,xx}.^2/Nranditer;
        popres.CA1V1corr_randSsqr{g,xx} = popres.CA1V1corr_randSsqr{g,xx} + tempCA1V1corr_randS{g,xx}.^2/Nranditer;
        popres.CA1V1corr_randXSsqr{g,xx} = popres.CA1V1corr_randXSsqr{g,xx} + tempCA1V1corr_randXS{g,xx}.^2/Nranditer;
        
        %     popres.CA1V1mscoh{g,xx} = popres.CA1V1mscoh{g,xx}/Nsession{g,xx};
        %     popres.CA1V1mscoh_rand{g,xx} = popres.CA1V1mscoh_rand{g,xx}/Nsession{g,xx}/Nranditer;
        %     popres.CA1V1mscoh_randX{g,xx} = popres.CA1V1mscoh_randX{g,xx}/Nsession{g,xx}/Nranditer;
        %     popres.CA1V1mscoh_randS{g,xx} = popres.CA1V1mscoh_randS{g,xx}/Nsession{g,xx}/Nranditer;
        %     popres.CA1V1mscoh_randXS{g,xx} = popres.CA1V1mscoh_randXS{g,xx}/Nsession{g,xx}/Nranditer;
        end
    end
end



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