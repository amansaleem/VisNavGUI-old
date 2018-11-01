function popresCorr = BatchJointCorrelationDec2(batch2p,Fcoherence)
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
    datadir = DIRS.data2p;
else
    expt = getExperimentList;
    datadir = 'D:\DATA\batch';%DIRS.multichanspikes;
end
                
popresCorr.Tsmthwin = 50;%250;%250;%150;%300;%40;%120;%50
popresCorr.Xsmthwin = 4;%2%1;%
popresCorr.SpdSmthWin = popresCorr.Tsmthwin;
popresCorr.SpeedThreshold = 5;
popresCorr.nspeedbins = 3;
popresCorr.neyebins = 1;
popresCorr.nthetaphsbins = 0;%1;%
popresCorr.cellstr = 'goodonly';%'All_50bins';%'goodonly';%'goodonly_unwrapped';%'goodonly';%'All';%
filesuffix_EXP = ['Twin' num2str(popresCorr.Tsmthwin) '_' 'Xwin' num2str(popresCorr.Xsmthwin) '_' 'spdth' num2str(popresCorr.SpeedThreshold) '_' num2str(popresCorr.nspeedbins) 'speedbins' '_' num2str(popresCorr.neyebins) 'eyebins' '_' num2str(popresCorr.nthetaphsbins) 'thetabins' '_' popresCorr.cellstr];
disp(filesuffix_EXP);

popresCorr.sampleRate = 60;
popresCorr.nSpdbins = 3;%1;%
popresCorr.nEyebins = 3;%3;%3;
popresCorr.nXbins = 20;%100;

lambdaSmooth = 2;
corrmaxlag = 180;
nanimal = numel(expt);

contval = [0.1:0.05:0.9];%[0.2 0.3 0.4];%[0.8 0.9];%
outvalcorr = 2;%[0 1 2 3 4 5];%[0 1 2 3 4];%

for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        dDIRname = [datadir  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];
        if exist([dDIRname filesep 'EXP_' filesuffix_EXP '.mat'],'file')
            if (expt(ianimal).goodCA1{iseries} == 1) && (expt(ianimal).goodV1{iseries} == 1)
                S = load([dDIRname filesep 'EXP_' filesuffix_EXP '.mat']);
                EXP = TVRData;
                EXP.Copyobj(S.EXP);
                
                cont_list = find(ismember(EXP.SubsetVal.contrast, contval));
                RL_list = find(ismember(EXP.SubsetVal.roomlength, [1]));
                outcome_list = find(ismember(EXP.SubsetVal.outcome, outvalcorr));
                Xrange = max(floor(EXP.Bayes.X));%1;%
                speeds = NaN(size(EXP.data.es.ballspeed));
                speeds(~isnan(EXP.data.es.ballspeed)) = smthInTime(EXP.data.es.ballspeed(~isnan(EXP.data.es.ballspeed)), popresCorr.sampleRate, popresCorr.SpdSmthWin, 'same', [], 'boxcar_centered');
                eyeX = NaN(size(EXP.data.es.eyeXpos));
                eyeX(~isnan(EXP.data.es.eyeXpos)) = smthInTime(EXP.data.es.eyeXpos(~isnan(EXP.data.es.eyeXpos)), popresCorr.sampleRate, popresCorr.SpdSmthWin, 'same', [], 'boxcar_centered');
                
                Speedbinned = cell(1,4);
                Eyebinned = cell(1,4);
                
                gooddecCA1 = false(1,3);
                gooddecV1 = false(1,3);
                for g = [2 1 3]
                    gooddecCA1(g) = expt(ianimal).goodCA1dec{g}{iseries};
                    gooddecV1(g) = expt(ianimal).goodV1dec{g}{iseries};
                end
                gooddecCA1(4) = gooddecCA1(2);
                gooddecV1(4) = gooddecV1(2);
                
                goodPostidx = cell(2,1);
                for iprobe = 1:2
                    goodPostidx{iprobe} =  ~isnan(sum(EXP.Bayes.Posterior0{iprobe,1},2));
                    EXP.Bayes.DecodingError0{iprobe,1}(~goodPostidx{iprobe}) = NaN;
                    EXP.Bayes.DecodingBias0{iprobe,1}(~goodPostidx{iprobe}) = NaN;
                end
                Prange = size(EXP.Bayes.Posterior0{1,1},2);
                
                PosError0{1} = EXP.Bayes.DecodingError0{1};
                PosError0{2} = EXP.Bayes.DecodingError0{2};

                tvec = 1:size(PosError0{2},1);
                PosError0{2} = interp1(tvec(goodPostidx{2}),PosError0{2}(goodPostidx{2}),1:size(PosError0{2},1));
                PosError0{1} = interp1(tvec(goodPostidx{1}),PosError0{1}(goodPostidx{1}),1:size(PosError0{1},1));
                PosError0{2} = PosError0{2}(:);
                PosError0{1} = PosError0{1}(:);
                    
                for g = [2 1 3 4]%2%
                    if gooddecCA1(g) && gooddecV1(g)
                        tidx = false(size(EXP.Bayes.X));
                        for cont = 1:numel(cont_list)
                            for r = 1:numel(RL_list)
                                for o = 1:numel(outcome_list)
                                    if ~isempty(EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th))
                                        tidx = tidx | EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                    end
                                end
                            end
                        end
                        
                        idxref = tidx;
                        
                        itraj = EXP.Bayes.X;
                        ntrajbins = max(itraj);
                        spdquantilelim = zeros(ntrajbins,2);
                        Speedbinned{g} = NaN(size(EXP.data.es.ballspeed));
                        if popresCorr.nSpdbins > 1
                            for spd = 1:popresCorr.nSpdbins
                                for xx = 1:ntrajbins
                                    spdquantilelim(xx,1) = quantile(speeds(idxref & itraj == xx),max(0,(spd-1)/popresCorr.nSpdbins));
                                    spdquantilelim(xx,2) = quantile(speeds(idxref & itraj == xx),min(1,(spd)/popresCorr.nSpdbins));
                                end
                                Speedbinned{g}(speeds >= spdquantilelim(itraj,1) & speeds < spdquantilelim(itraj,2)) = spd;
                                if spd == 1
                                    Speedbinned{g}(speeds <= spdquantilelim(itraj,1)) = spd;
                                end
                                if spd == popresCorr.nSpdbins
                                    Speedbinned{g}(speeds >= spdquantilelim(itraj,2)) = spd;
                                end
                            end
                        else
                            Speedbinned{g} = ones(size(EXP.data.es.ballspeed));
                        end
                        
                        eyequantilelim = zeros(ntrajbins,2);
                        Eyebinned{g} = NaN(size(EXP.data.es.eyeXpos));
                        if popresCorr.nEyebins > 1
                            for ieye = 1:popresCorr.nEyebins
                                for xx = 1:ntrajbins
                                    eyequantilelim(xx,1) = quantile(eyeX(idxref & itraj == xx),max(0,(ieye-1)/popresCorr.nEyebins));
                                    eyequantilelim(xx,2) = quantile(eyeX(idxref & itraj == xx),min(1,(ieye)/popresCorr.nEyebins));
                                end
                                Eyebinned{g}(eyeX >= eyequantilelim(itraj,1) & eyeX < eyequantilelim(itraj,2)) = ieye;
                                if ieye == 1
                                    Eyebinned{g}(eyeX <= eyequantilelim(itraj,1)) = ieye;
                                end
                                if ieye == popresCorr.nEyebins
                                    Eyebinned{g}(eyeX >= eyequantilelim(itraj,2)) = ieye;
                                end
                            end
                        else
                            Eyebinned{g} = ones(size(EXP.data.es.eyeXpos));
                        end
                        
                        Xbinned = floor(((EXP.Bayes.X-1)/Xrange)*popresCorr.nXbins)+1;%
                        %                         fq = [0.1 4];
                        %                         for xdec = 1:size(PosError0{1},2)
                        %                             vec = PosError0{1}(:,xdec);
                        %                             vec = LFPfilter(vec, fq(1), fq(2), 60);
                        %                             PosError0{1}(:,xdec) = vec;%
                        %                             vec = PosError0{2}(:,xdec);
                        %                             vec = LFPfilter(vec, fq(1), fq(2), 60);
                        %                             PosError0{2}(:,xdec) = vec;%
                        %                         end
                    end
                end
                    
                for g = [2 1 3 4]%2%
                    if gooddecCA1(g) && gooddecV1(g)
                        tidx = false(size(EXP.Bayes.X));
                        for cont = 1:numel(cont_list)
                            for r = 1:numel(RL_list)
                                for o = 1:numel(outcome_list)
                                    if g < 4
                                        if ~isempty(EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th))
                                            tidx = tidx | EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                        end
                                    else
                                        if ~isempty(EXP.getSubsets(cont_list(cont),[1 2 3],RL_list(r),outcome_list(o), EXP.Bayes.speed_th))
                                            tidx = tidx | EXP.getSubsets(cont_list(cont),[1 2 3],RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                        end
                                    end
                                end
                            end
                        end
                        
                        
                        nXbins = popresCorr.nXbins;
                        nSpdbins = popresCorr.nSpdbins;
                        nEyebins = popresCorr.nEyebins;
                        s_nDataPoints = zeros(nXbins,nSpdbins,nEyebins);
                        s_CA1V1Joint = NaN(nXbins,nSpdbins,nEyebins,Prange,Prange);
                        s_CA1V1Joint_shuffled = NaN(nXbins,nSpdbins,nEyebins,Prange,Prange);
                        
                        s_CA1V1Jointdiag = NaN(nXbins,nSpdbins,nEyebins,Prange);
                        s_CA1V1Jointdiag_shuffled = NaN(nXbins,nSpdbins,nEyebins,Prange);
                        
                        xsmth = 0;
                        for ieye = 1:nEyebins
                            for spd = 1:nSpdbins
                                for xx = 1:nXbins
                                    tidx_corr = tidx & ~isnan(sum(PosError0{1},2)) & ~isnan(sum(PosError0{2},2)) & ismember(Xbinned,mod((xx-xsmth:xx+xsmth)-1,popresCorr.nXbins)+1) & Speedbinned{g} == spd & Eyebinned{g} == ieye;
                                    idx_corrCA1 = find(tidx_corr);
                                    idx_corrV1 = find(tidx_corr);
                                    s_nDataPoints(xx,spd,ieye) = numel(idx_corrCA1);
                                    s_CA1V1Joint(xx,spd,ieye,:,:) = smoothhist2D_corrected([PosError0{1}(idx_corrCA1)+floor(Prange/2) PosError0{2}(idx_corrV1)+floor(Prange/2)],[2 2], [Prange Prange], 1:Prange, 1:Prange, true, true,false);
                                    shuffledidxV1 = randperm(numel(idx_corrV1));
                                    s_CA1V1Joint_shuffled(xx,spd,ieye,:,:) = smoothhist2D_corrected([PosError0{1}(idx_corrCA1)+floor(Prange/2) PosError0{2}(idx_corrV1(shuffledidxV1))+floor(Prange/2)],[2 2], [Xrange Prange], 1:Xrange, 1:Prange, true, true,false);
                                    
                                    s_CA1V1Jointdiag(xx,spd,ieye,:) = diag(squeeze(s_CA1V1Joint(xx,spd,ieye,:,:)));
                                    s_CA1V1Jointdiag_shuffled(xx,spd,ieye,:) = diag(squeeze(s_CA1V1Joint_shuffled(xx,spd,ieye,:,:)));
                                end
                            end
                        end
                        popresCorr.s_CA1V1Cov{ianimal,iseries,g} = squeeze(nanmean(nanmean(nanmean(s_CA1V1Joint,1),2),3));
                        popresCorr.s_CA1V1Cov_shuffled{ianimal,iseries,g} = squeeze(nanmean(nanmean(nanmean(s_CA1V1Joint_shuffled,1),2),3));
                        
                        popresCorr.s_CA1V1Covdiag{ianimal,iseries,g} = squeeze(nanmean(nanmean(nanmean(s_CA1V1Jointdiag,1),2),3));
                        popresCorr.s_CA1V1Covdiag_shuffled{ianimal,iseries,g} = squeeze(nanmean(nanmean(nanmean(s_CA1V1Jointdiag_shuffled,1),2),3));
                    end
                end
            
                if Fcoherence
                    params.Fs = 60;
                    params.fpass = [0 20];
                    params.tapers = [1 1];
                    params.trialave = 0;%1;
                    params.err = 0;
                    nphsbins = 18;
                    movingwin = [1 0.1];%[1 0.1];
                    ErrPost_V1 = PosError0{2};%EXP.Bayes.PosError0{2};% - repmat(mean(PosError0{2}(tidx,:),1),[sum(tidx) 1]);
                    ErrPost_CA1 = PosError0{1};%EXP.Bayes.PosError0{1};%-repmat(mean(PosError0{1}(tidx,:),1),[sum(tidx) 1]);
                    
                    %                         ErrPost_V1 = getCircularAverage(ErrPost_V1',0,1);
                    %                         ErrPost_CA1 = getCircularAverage(ErrPost_CA1',0,1);
                    
                    %                         [~,~,S12,S1,S2,t,f] = cohgramc(ErrPost_CA1,ErrPost_V1,movingwin,params);
                    
                    params.complex = true;
                    [S1,~,~]=mtspecgramc(ErrPost_CA1,movingwin,params);
                    [S2,t,f]=mtspecgramc(ErrPost_V1,movingwin,params);
                    
                    popresCorr.t = t;
                    popresCorr.f = f;
                    %
                    for g = [2 1 3 4]
                        if gooddecCA1(g) && gooddecV1(g)
                            tidx = false(size(EXP.Bayes.X));
                            for cont = 1:numel(cont_list)
                                for r = 1:numel(RL_list)
                                    for o = 1:numel(outcome_list)
                                        if g < 4
                                            if ~isempty(EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th))
                                                tidx = tidx | EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                            end
                                        else
                                            if ~isempty(EXP.getSubsets(cont_list(cont),[1 2 3],RL_list(r),outcome_list(o), EXP.Bayes.speed_th))
                                                tidx = tidx | EXP.getSubsets(cont_list(cont),[1 2 3],RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                            end
                                        end
                                    end
                                end
                            end
                            tidx_spec = tidx;
                            tidx_spec = interp1((1:numel(tidx_spec))*(1/params.Fs),single(tidx_spec),popresCorr.t,'nearest');
                            tidx_spec = tidx_spec > 0.5;
                            
                            S12 = zeros(size(S1,3),size(S1,3),numel(popresCorr.f));
                            for ifq = 1:numel(popresCorr.f)
                                for itapers = 1:size(S1,4)
                                    for x1 = 1:size(S1,3)
                                        for x2 = 1:size(S1,3)
                                            S12(x1,x2,ifq) = S12(x1,x2,ifq) + squeeze(sum(squeeze(conj(S2(tidx_spec,ifq,x1,itapers))).*squeeze(S1(tidx_spec,ifq,x2,itapers)),1))'/sum(tidx_spec)/size(S2,4);
                                        end
                                    end
                                end
                            end
                            S10 = zeros(size(S1,3),numel(popresCorr.f));
                            for ifq = 1:numel(popresCorr.f)
                                for itapers = 1:size(S1,4)
                                    S10(:,ifq) = S10(:,ifq) + squeeze(sum(squeeze(conj(S1(tidx_spec,ifq,:,itapers))).*squeeze(S1(tidx_spec,ifq,:,itapers)),1))'/sum(tidx_spec)/size(S1,4);
                                end
                            end
                            S20 = zeros(size(S2,3),numel(popresCorr.f));
                            for ifq = 1:numel(popresCorr.f)
                                for itapers = 1:size(S1,4)
                                    S20(:,ifq) = S20(:,ifq) + squeeze(sum(squeeze(conj(S2(tidx_spec,ifq,:,itapers))).*squeeze(S2(tidx_spec,ifq,:,itapers)),1))'/sum(tidx_spec)/size(S2,4);
                                end
                            end
                            
                            popresCorr.s_CA1Spec{ianimal,iseries,g} = squeeze(mean(abs(mean(S1(tidx_spec,:,:,:),4)),1));
                            popresCorr.s_V1Spec{ianimal,iseries,g} = squeeze(mean(abs(mean(S2(tidx_spec,:,:,:),4)),1));
                            
                            popresCorr.s_CA1V1Coh{ianimal,iseries,g} = zeros(size(S12));
                            popresCorr.s_CA1V1Phi{ianimal,iseries,g} = zeros(size(S12));
                            popresCorr.s_CA1V1CovSpec{ianimal,iseries,g} = zeros(size(S12));
                            for ifq = 1:numel(popresCorr.f)
                                popresCorr.s_CA1V1Coh{ianimal,iseries,g}(:,:,ifq) = abs(S12(:,:,ifq)./sqrt(S20(:,ifq)*S10(:,ifq)'));
                                popresCorr.s_CA1V1Phi{ianimal,iseries,g}(:,:,ifq) = angle(S12(:,:,ifq)./sqrt(S20(:,ifq)*S10(:,ifq)'));
                                
                                S1temp = squeeze(abs(mean(S1(tidx_spec,ifq,:,:),4)));
                                S1temp = S1temp - repmat(mean(S1temp,1),[size(S1temp,1) 1]);
                                S2temp = squeeze(abs(mean(S2(tidx_spec,ifq,:,:),4)));
                                S2temp = S2temp - repmat(mean(S2temp,1),[size(S2temp,1) 1]);
                                popresCorr.s_CA1V1CovSpec{ianimal,iseries,g}(:,:,ifq) = (S2temp'*S1temp)./sqrt(sum(S2temp.^2,1)'*sum(S1temp.^2,1));
                            end
                        end
                    end
                    
                    Trandshift = 1000;
                    Nrand = 0;%20;
                    if Nrand > 0
                        Cohrand_iter = zeros(Nrand,4,numel(popresCorr.f));
                        Phirand_iter = zeros(Nrand,4,numel(popresCorr.f));
                        for kiter = 1:Nrand
                            [~,~,S12,S1,S2,t,~]=cohgramc(circshift(ErrPost_CA1,kiter*Trandshift,1),ErrPost_V1,movingwin,params);
                            for g = [2 1 3 4]%2%
                                if gooddecCA1(g) && gooddecV1(g)
                                    tidx = false(size(EXP.Bayes.X));
                                    for cont = 1:numel(cont_list)
                                        for r = 1:numel(RL_list)
                                            for o = 1:numel(outcome_list)
                                                if g < 4
                                                    if ~isempty(EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th))
                                                        tidx = tidx | EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                                    end
                                                else
                                                    tidx = true(size(tidx));
                                                end
                                            end
                                        end
                                    end
                                    tidx_spec = interp1((1:numel(tidx))*(1/params.Fs),single(tidx),t,'nearest');
                                    tidx_spec = tidx_spec > 0.5;
                                    Cohrand_iter(kiter,g,:) = squeeze(abs(nanmean(nanmean(S12(tidx_spec,:,:),3),1)./sqrt(nanmean(nanmean(S1(tidx_spec,:,:),3),1).*nanmean(nanmean(S2(tidx_spec,:,:),3),1))));
                                    Phirand_iter(kiter,g,:) = squeeze(angle(nanmean(nanmean(S12(tidx_spec,:,:),3),1)./sqrt(nanmean(nanmean(S1(tidx_spec,:,:),3),1).*nanmean(nanmean(S2(tidx_spec,:,:),3),1))));
                                    
                                    %                                         Cohrand_iter(kiter,g,:) = squeeze(nanmean(Coh(tidx_spec,:,:),1));%squeeze(abs(nanmean(S12(tidx_spec,:,:),1)./sqrt(nanmean(S1(tidx_spec,:,:),1).*nanmean(S2(tidx_spec,:,:),1))))';
                                    %                                         Phirand_iter(kiter,g,:) = squeeze(circ_mean(phi(tidx_spec,:,:),[],1));
                                end
                            end
                        end
                        for g = [2 1 3 4]
                            popresCorr.s_CA1V1Coh_rand{ianimal,iseries,g} = squeeze(nanmean(Cohrand_iter(:,g,:),1));
                            popresCorr.s_CA1V1Coh_randstd{ianimal,iseries,g} = squeeze(nanstd(Cohrand_iter(:,g,:),0,1));
                            popresCorr.s_CA1V1Phi_rand{ianimal,iseries,g} = squeeze(circ_mean(Phirand_iter(:,g,:),[],1));
                            popresCorr.s_CA1V1Phi_randstd{ianimal,iseries,g} = squeeze(circ_std(Phirand_iter(:,g,:),[],[],1));
                        end
                    else
                        for g = [2 1 3 4]
                            popresCorr.s_CA1V1Coh_rand{ianimal,iseries,g} = 0;
                            popresCorr.s_CA1V1Coh_randstd{ianimal,iseries,g} = 0;
                            popresCorr.s_CA1V1Phi_rand{ianimal,iseries,g} = 0;
                            popresCorr.s_CA1V1Phi_randstd{ianimal,iseries,g} = 0;
                        end
                    end
                end
            end
        end
    end
end

%now we do the grand averages by averaging across sessions
for g = [2 1 3 4]
    count = 0;
    popresCorr.CA1V1Cov{g} = 0;
    popresCorr.CA1V1Cov_shuffled{g} = 0;
    popresCorr.CA1V1Covdiag{g} = 0;
    popresCorr.CA1V1Covdiag_shuffled{g} = 0;
    if Fcoherence
        popresCorr.CA1V1Coh{g} = 0;
        popresCorr.CA1Spec{g} = 0;
        popresCorr.V1Spec{g} = 0;
        popresCorr.CA1V1CovSpec{g} = 0;
        popresCorr.CA1V1Coh_rand{g} = 0;
        popresCorr.CA1V1Coh_randstd{g} = 0;
    end
    for ianimal = 1:size(popresCorr.s_CA1V1Cov,1)
        for iseries = 1:size(popresCorr.s_CA1V1Cov,2)
            if ~isempty(popresCorr.s_CA1V1Cov{ianimal,iseries,g})
               if sum(isnan(popresCorr.s_CA1V1Cov{ianimal,iseries,g}(:))) == 0
                    popresCorr.CA1V1Cov{g} = popresCorr.CA1V1Cov{g} + popresCorr.s_CA1V1Cov{ianimal,iseries,g};
                    popresCorr.CA1V1Cov_shuffled{g} = popresCorr.CA1V1Cov_shuffled{g} + popresCorr.s_CA1V1Cov_shuffled{ianimal,iseries,g};
                    popresCorr.CA1V1Covdiag{g} = popresCorr.CA1V1Covdiag{g} + popresCorr.s_CA1V1Covdiag{ianimal,iseries,g};
                    popresCorr.CA1V1Covdiag_shuffled{g} = popresCorr.CA1V1Covdiag_shuffled{g} + popresCorr.s_CA1V1Covdiag_shuffled{ianimal,iseries,g};
                    
                    if Fcoherence
                        popresCorr.CA1V1Coh{g} = popresCorr.CA1V1Coh{g} + popresCorr.s_CA1V1Coh{ianimal,iseries,g};
                        popresCorr.CA1Spec{g} = popresCorr.CA1Spec{g} + popresCorr.s_CA1Spec{ianimal,iseries,g};
                        popresCorr.V1Spec{g} = popresCorr.V1Spec{g} + popresCorr.s_V1Spec{ianimal,iseries,g};
                        popresCorr.CA1V1CovSpec{g} = popresCorr.CA1V1CovSpec{g} + popresCorr.s_CA1V1CovSpec{ianimal,iseries,g};
                        
                        popresCorr.CA1V1Coh_rand{g} = popresCorr.CA1V1Coh_rand{g} + popresCorr.s_CA1V1Coh_rand{ianimal,iseries,g};
                        popresCorr.CA1V1Coh_randstd{g} = popresCorr.CA1V1Coh_randstd{g} + popresCorr.s_CA1V1Coh_randstd{ianimal,iseries,g}; 
                    end
                    count = count + 1;
                end
            end
        end
    end
    popresCorr.CA1V1Cov{g} = popresCorr.CA1V1Cov{g}./count;
    popresCorr.CA1V1Cov_shuffled{g} = popresCorr.CA1V1Cov_shuffled{g}./count;
    popresCorr.CA1V1Covdiag{g} = popresCorr.CA1V1Covdiag{g}./count;
    popresCorr.CA1V1Covdiag_shuffled{g} = popresCorr.CA1V1Covdiag_shuffled{g}./count;
    if Fcoherence
        popresCorr.CA1V1Coh{g} = popresCorr.CA1V1Coh{g}./count;
        popresCorr.CA1Spec{g} = popresCorr.CA1Spec{g}./count;
        popresCorr.V1Spec{g} = popresCorr.V1Spec{g}./count;
        popresCorr.CA1V1CovSpec{g} = popresCorr.CA1V1CovSpec{g}./count;
        
        popresCorr.CA1V1Coh_rand{g} = popresCorr.CA1V1Coh_rand{g}./count;
        popresCorr.CA1V1Coh_randstd{g} = popresCorr.CA1V1Coh_randstd{g}./count;
    end
end
end



%use the following to visualize the data until there's a GUI or fig
%associated
% for g = [2 1 3]
% for tlag = 1:1
% count = 0;
% CovMat_Ave{tlag,g} = 0;
% CovMat_Shuffled{tlag,g} = 0;
% for ianimal = 1:size(popresCorr.s_CA1V1Cov,1)
% for iseries = 1:size(popresCorr.s_CA1V1Cov,2)
% if ~isempty(popresCorr.s_CA1V1Cov{ianimal,iseries,tlag,g})
% %This will average across sessions than across all trials
% mat = popresCorr.s_CA1V1Joint{ianimal,iseries,tlag,g}./repmat(popresCorr.s_nDataPoints{ianimal,iseries,tlag,g},[1 1 100 100]);
% if sum(isnan(mat)) == 0
% CovMat_Ave{tlag,g} = CovMat_Ave{tlag,g} + mat;
% CovMat_Shuffled{tlag,g} = CovMat_Shuffled{tlag,g} + popresCorr.s_CA1V1Joint_shuffled{ianimal,iseries,tlag,g}./repmat(popresCorr.s_nDataPoints{ianimal,iseries,tlag,g},[1 1 100 100]);
% count = count+ 1;%popresCorr.s_nDataPoints{ianimal,iseries,tlag,g};
% end
% end
% end
% end
% CovMat_Ave{tlag,g} = CovMat_Ave{tlag,g}/count;%./repmat(count,[1 1 100 100]);
% CovMat_Shuffled{tlag,g} = CovMat_Shuffled{tlag,g}/count;%./repmat(count,[1 1 100 100]);
% end
% end
% 
% figure;
% mat1 = [];
% mat2 = [];
% for g = 1:3
% mat1{g} = 0;
% mat2{g} = 0;
% end
% for xx = 1:100
% for g = 1:3
% mattemp1 = CovMat_Ave{1,g}(xx,:,:,:);
% mattemp2 = CovMat_Shuffled{1,g}(xx,:,:,:);
% mattemp1 = squeeze(nanmean(nanmean(mattemp1,1),2));
% mattemp2 = squeeze(nanmean(nanmean(mattemp2,1),2));
% mat1{g} = (mat1{g}*(xx-1) + mattemp1)/xx;
% mat2{g} = (mat2{g}*(xx-1) + mattemp2)/xx;
% subplot(2,3,g)
% imagesc(mat1{g});
% set(gca,'Clim',[-1 1],'Ydir','normal');
% subplot(2,3,g+3)
% imagesc(mat2{g});
% set(gca,'Clim',[-1 1],'Ydir','normal');
% xlabel(num2str(xx))
% end
% pause
% end


%or alternatively
% for g = [2 1 3]
% count = 0;
% crosscount = 0;
% CovMat_Ave{g} = 0;
% CovLagMat_Ave{g} = 0;
% CrossMat_Ave{g} = 0;
% CovMat_Shuffled{g} = 0;
% CovLagMat_Shuffled{g} = 0;
% CrossMat_Shuffled{g} = 0;
% CovMat_stdCA1{g} = 0;
% CovMat_stdV1{g} = 0;
% for ianimal = 1:size(popresCorr.s_CA1V1Cov,1)
% for iseries = 1:size(popresCorr.s_CA1V1Cov,2)
% if ~isempty(popresCorr.s_CA1V1Cov{ianimal,iseries,g})
% mat = popresCorr.s_CA1V1Joint{ianimal,iseries,g};
% matshf = popresCorr.s_CA1V1Joint_shuffled{ianimal,iseries,g};
% matcross = zeros(size(mat));
% matcross_shf = zeros(size(mat));
% for xx = 1:size(mat,1)
% for spd = 1:size(mat,2)
% matcross(xx,spd,:,:) = squeeze(mat(xx,spd,:,:))./(squeeze(popresCorr.s_V1marginal{ianimal,iseries,g}(xx,spd,:))*squeeze(popresCorr.s_CA1marginal{ianimal,iseries,g}(xx,spd,:))').^0.5;
% matcross_shf(xx,spd,:,:) = squeeze(matshf(xx,spd,:,:))./(squeeze(popresCorr.s_V1marginal_shuffled{ianimal,iseries,g}(xx,spd,:))*squeeze(popresCorr.s_CA1marginal_shuffled{ianimal,iseries,g}(xx,spd,:))').^0.5;
% end
% end
% if sum(isnan(mat)) == 0
% CovMat_Ave{g} = CovMat_Ave{g} + mat;
% CrossMat_Ave{g} = CrossMat_Ave{g} + matcross;
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g} + matcross_shf;
% CovLagMat_Ave{g} = CovLagMat_Ave{g} + popresCorr_250.s_CA1V1Jointdiag{ianimal,iseries,g};
% crosscount = crosscount + 1;
% CovMat_stdCA1{g} = CovMat_stdCA1{g} + popresCorr_250.s_CA1marginal{ianimal,iseries,g};
% CovMat_stdV1{g} = CovMat_stdV1{g} + popresCorr_250.s_V1marginal{ianimal,iseries,g};
% CovMat_Shuffled{g} = CovMat_Shuffled{g} + popresCorr_250.s_CA1V1Joint_shuffled{ianimal,iseries,g};
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g} + popresCorr_250.s_CA1V1Jointdiag_shuffled{ianimal,iseries,g};
% count = count+ popresCorr_250.s_nDataPoints{ianimal,iseries,g};
% end
% end
% end
% end
% CovMat_Ave{g} = CovMat_Ave{g}./repmat(count,[1 1 100 100]);
% CovLagMat_Ave{g} = CovLagMat_Ave{g}./repmat(count,[1 1 241 100]);
% CrossMat_Ave{g} = CrossMat_Ave{g}/crosscount;
% CovMat_Shuffled{g} = CovMat_Shuffled{g}./repmat(count,[1 1 100 100]);
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g}./repmat(count,[1 1 241 100]);
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g}/crosscount;
% CovMat_stdCA1{g} = CovMat_stdCA1{g}./repmat(count,[1 1 100]);
% CovMat_stdV1{g} = CovMat_stdV1{g}./repmat(count,[1 1 100]);
% end
% 
% 
% same as above but averaging across sessions
% for g = [2 1 3]
% count = 0;
% crosscount = 0;
% CovMat_Ave{g} = 0;
% CovLagMat_Ave{g} = 0;
% CrossMat_Ave{g} = 0;
% CovMat_Shuffled{g} = 0;
% CovLagMat_Shuffled{g} = 0;
% CrossMat_Shuffled{g} = 0;
% CovMat_stdCA1{g} = 0;
% CovMat_stdV1{g} = 0;
% for ianimal = 1:size(popresCorr.s_CA1V1Cov,1)
% for iseries = 1:size(popresCorr.s_CA1V1Cov,2)
% if ~isempty(popresCorr.s_CA1V1Cov{ianimal,iseries,g})
% mat = popresCorr.s_CA1V1Joint{ianimal,iseries,g}./repmat(popresCorr.s_nDataPoints{ianimal,iseries,g},[1 1 1 100 100]);
% matshf = popresCorr.s_CA1V1Joint_shuffled{ianimal,iseries,g}./repmat(popresCorr.s_nDataPoints{ianimal,iseries,g},[1 1 1 100 100]);
% matcross = zeros(size(mat));
% matcross_shf = zeros(size(mat));
% for xx = 1:size(mat,1)
% for spd = 1:size(mat,2)
% for ieye = 1:size(mat,3)
% matcross(xx,spd,ieye,:,:) = squeeze(popresCorr.s_CA1V1Joint{ianimal,iseries,g}(xx,spd,ieye,:,:))./(squeeze(popresCorr.s_V1marginal{ianimal,iseries,g}(xx,spd,ieye,:))*squeeze(popresCorr.s_CA1marginal{ianimal,iseries,g}(xx,spd,ieye,:))').^0.5;
% matcross_shf(xx,spd,ieye,:,:) = squeeze(popresCorr.s_CA1V1Joint_shuffled{ianimal,iseries,g}(xx,spd,ieye,:,:))./(squeeze(popresCorr.s_V1marginal_shuffled{ianimal,iseries,g}(xx,spd,ieye,:))*squeeze(popresCorr.s_CA1marginal_shuffled{ianimal,iseries,g}(xx,spd,ieye,:))').^0.5;
% end
% end
% end
% if sum(isnan(mat)) == 0
% CovMat_Ave{g} = CovMat_Ave{g} + mat;
% CrossMat_Ave{g} = CrossMat_Ave{g} + matcross;
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g} + matcross_shf;
% CovLagMat_Ave{g} = CovLagMat_Ave{g} + popresCorr.s_CA1V1Jointdiag{ianimal,iseries,g}./repmat(popresCorr.s_nDataPoints{ianimal,iseries,g},[1 1 1 1 100]);
% crosscount = crosscount + 1;
% CovMat_stdCA1{g} = CovMat_stdCA1{g} + popresCorr.s_CA1marginal{ianimal,iseries,g};
% CovMat_stdV1{g} = CovMat_stdV1{g} + popresCorr.s_V1marginal{ianimal,iseries,g};
% CovMat_Shuffled{g} = CovMat_Shuffled{g} + matshf;
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g} + popresCorr.s_CA1V1Jointdiag_shuffled{ianimal,iseries,g}./repmat(popresCorr.s_nDataPoints{ianimal,iseries,g},[1 1 1 1 100]);
% count = count+ 1;%popresCorr_250.s_nDataPoints{ianimal,iseries,g};
% end
% end
% end
% end
% CovMat_Ave{g} = CovMat_Ave{g}./repmat(count,[100 3 3 100 100]);
% CovLagMat_Ave{g} = CovLagMat_Ave{g}./repmat(count,[100 3 3 241 100]);
% CrossMat_Ave{g} = CrossMat_Ave{g}/crosscount;
% CovMat_Shuffled{g} = CovMat_Shuffled{g}./repmat(count,[100 3 3 100 100]);
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g}./repmat(count,[100 3 3 241 100]);
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g}/crosscount;
% CovMat_stdCA1{g} = CovMat_stdCA1{g}./repmat(count,[100 3 3 100]);
% CovMat_stdV1{g} = CovMat_stdV1{g}./repmat(count,[100 3 3 100]);
% end