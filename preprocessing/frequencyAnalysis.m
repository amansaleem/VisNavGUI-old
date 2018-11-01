function res = frequencyAnalysis(resCA1, resV1)
for k = 1:2
    res(k).Spec = [];
    res(k).thetafq = [6 9];
end

fieldname = 'LFPfilt';%'CSDfilt';
params.Fs = 1000;
params.fpass = [0 80];
params.tapers = [3 5];
params.trialave = 0;
params.err = 0;
movingwin = [5 1];

if ~isempty(resCA1)
    CA1chref = resCA1.CA1chref;
    LFPCA1 = resCA1.LFPfilt;%resCA1.LFPfilt(:,CA1chref);
end

nfold = 20;
speed_th = 5;outcomeval = 2;gainval = [0.4 0.5 0.6];
    
spkV1 = [];
spkCA1 = [];
if ~isempty(resV1)
    [lfp_S1_V1,~,~] = mtspecgramc(resV1.(fieldname),movingwin,params);
    for ichan = 1:numel(resV1.chans)
        spkV1(ichan).times = resV1.chans(ichan).spiketimes*0.001;
    end
    spkV1 = spkV1(resV1.Fgoodunit);
    [~,~,spk_S12_V1,lfp_S1_CA1_V1,spk_S2_V1,t,f]=cohgramcpt(repmat(LFPCA1,[1 numel(spkV1)]),spkV1,movingwin,params);
end

if ~isempty(resCA1)
    for ichan = 1:numel(resCA1.chans)
        spkCA1(ichan).times = resCA1.chans(ichan).spiketimes*0.001;
    end
    spkCA1 = spkCA1(resCA1.Fgoodunit);
    [~,~,spk_S12_CA1,lfp_S1_CA1,spk_S2_CA1,t,f]=cohgramcpt(repmat(LFPCA1,[1 numel(spkCA1)]),spkCA1,movingwin,params);
end
res(1).t = t;
res(1).f = f;
res(2).t = t;
res(2).f = f;

samplingdt = 0.001;
if ~isempty(resV1)
    res(2).lfp_SpecAll = squeeze(abs(nanmean(nanmean(lfp_S1_V1(:,:,:),1),3)));
    res(2).spk_SpecChAll = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
    res(2).spk_SpecChAll(:,resV1.Fgoodunit) = squeeze(abs(nanmean(spk_S2_V1(:,:,:),1)));
    res(2).spk_SpecAll = squeeze(abs(nanmean(nanmean(spk_S2_V1(:,:,:),1),3)));
    res(2).spk_CohSpecAll = squeeze(abs(nanmean(nanmean(spk_S12_V1(:,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1_V1(:,:,:),1),3).*nanmean(nanmean(spk_S2_V1(:,:,:),1),3))));
    res(2).spk_PhsCohSpecAll = squeeze(angle(nanmean(nanmean(spk_S12_V1(:,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1_V1(:,:,:),1),3).*nanmean(nanmean(spk_S2_V1(:,:,:),1),3))));
    res(2).spk_CohSpecChAll = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
    res(2).spk_CohSpecChAll(:,resV1.Fgoodunit) = squeeze(abs(nanmean(spk_S12_V1(:,:,:),1)./sqrt(nanmean(lfp_S1_CA1_V1(:,:,:),1).*nanmean(spk_S2_V1(:,:,:),1))));
    res(2).spk_PhsCohSpecChAll = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
    res(2).spk_PhsCohSpecChAll(:,resV1.Fgoodunit) = squeeze(angle(nanmean(spk_S12_V1(:,:,:),1)./sqrt(nanmean(lfp_S1_CA1_V1(:,:,:),1).*nanmean(spk_S2_V1(:,:,:),1))));
    
    tidx = true(size(lfp_S1_V1,1),1);
    idx = find(tidx);
    [CVO] = crossValPartition(1:numel(idx), nfold, false, floor(60/movingwin(2)));
    lfp_Spec2_std = 0;
    spk_SpecCh2_std = 0;
    spk_Spec2_std = 0;
    spk_CohSpec2_std = 0;
    spk_PhsCohSpec2_std = 0;
    spk_CohSpecCh2_std = 0;
    spk_PhsCohSpecCh2_std = 0;
    for kiter = 1:nfold
        tidxCVO = false(size(tidx));
        tidxCVO(idx(CVO.train{kiter})) = true;
        
        res(2).lfp_SpecAll_CVO{kiter} = squeeze(abs(nanmean(nanmean(lfp_S1_V1(tidxCVO,:,:),1),3)));
        lfp_Spec2_iter = res(2).lfp_SpecAll_CVO{kiter};
        lfp_Spec2_std = lfp_Spec2_std + (nfold - 1)/nfold*(lfp_Spec2_iter - res(2).lfp_SpecAll).^2;
        
        res(2).spk_SpecChAll_CVO{kiter} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
        res(2).spk_SpecChAll_CVO{kiter}(:,resV1.Fgoodunit) = squeeze(abs(nanmean(spk_S2_V1(tidxCVO,:,:),1)));
        spk_SpecCh2_iter = res(2).spk_SpecChAll_CVO{kiter}(:,resV1.Fgoodunit);
        spk_SpecCh2_std = spk_SpecCh2_std + (nfold - 1)/nfold*(spk_SpecCh2_iter - res(2).spk_SpecChAll(:,resV1.Fgoodunit)).^2;
        
        res(2).spk_SpecAll_CVO{kiter} = squeeze(abs(nanmean(nanmean(spk_S2_V1(tidxCVO,:,:),1),3)));
        spk_Spec2_iter = res(2).spk_SpecAll_CVO{kiter};
        spk_Spec2_std = spk_Spec2_std + (nfold - 1)/nfold*(spk_Spec2_iter - res(2).spk_SpecAll).^2;
        
        res(2).spk_CohSpecAll_CVO{kiter} = squeeze(abs(nanmean(nanmean(spk_S12_V1(tidxCVO,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1_V1(tidxCVO,:,:),1),3).*nanmean(nanmean(spk_S2_V1(tidxCVO,:,:),1),3))));
        spk_CohSpec2_iter = res(2).spk_CohSpecAll_CVO{kiter};
        spk_CohSpec2_std = spk_CohSpec2_std + (nfold - 1)/nfold*(spk_CohSpec2_iter - res(2).spk_CohSpecAll).^2;
        
        res(2).spk_PhsCohSpecAll_CVO{kiter} = squeeze(angle(nanmean(nanmean(spk_S12_V1(tidxCVO,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1_V1(tidxCVO,:,:),1),3).*nanmean(nanmean(spk_S2_V1(tidxCVO,:,:),1),3))));
        spk_PhsCohSpec2_iter = res(2).spk_PhsCohSpecAll_CVO{kiter};
        spk_PhsCohSpec2_std = spk_PhsCohSpec2_std + (nfold - 1)/nfold*circ_dist(spk_PhsCohSpec2_iter,res(2).spk_PhsCohSpecAll).^2;
        
        res(2).spk_CohSpecChAll_CVO{kiter} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
        res(2).spk_CohSpecChAll_CVO{kiter}(:,resV1.Fgoodunit) = squeeze(abs(nanmean(spk_S12_V1(tidxCVO,:,:),1)./sqrt(nanmean(lfp_S1_CA1_V1(tidxCVO,:,:),1).*nanmean(spk_S2_V1(tidxCVO,:,:),1))));
        spk_CohSpecCh2_iter = res(2).spk_CohSpecChAll_CVO{kiter}(:,resV1.Fgoodunit);
        spk_CohSpecCh2_std = spk_CohSpecCh2_std + (nfold - 1)/nfold*(spk_CohSpecCh2_iter - res(2).spk_CohSpecChAll(:,resV1.Fgoodunit)).^2;
        
        res(2).spk_PhsCohSpecChAll_CVO{kiter} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
        res(2).spk_PhsCohSpecChAll_CVO{kiter}(:,resV1.Fgoodunit) = squeeze(angle(nanmean(spk_S12_V1(tidxCVO,:,:),1)./sqrt(nanmean(lfp_S1_CA1_V1(tidxCVO,:,:),1).*nanmean(spk_S2_V1(tidxCVO,:,:),1))));
        spk_PhsCohSpecCh2_iter = res(2).spk_PhsCohSpecChAll_CVO{kiter}(:,resV1.Fgoodunit);
        spk_PhsCohSpecCh2_std = spk_PhsCohSpecCh2_std + (nfold - 1)/nfold*circ_dist(spk_PhsCohSpecCh2_iter,res(2).spk_PhsCohSpecChAll(:,resV1.Fgoodunit)).^2;
    end
    res(2).lfp_SpecAll_SE = sqrt(lfp_Spec2_std);
    res(2).spk_SpecChAll_SE = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
    res(2).spk_SpecChAll_SE(:,resV1.Fgoodunit) = sqrt(spk_SpecCh2_std);
    res(2).spk_SpecAll_SE = sqrt(spk_Spec2_std);
    res(2).spk_CohSpecAll_SE = sqrt(spk_CohSpec2_std);
    res(2).spk_PhsCohSpecAll_SE = sqrt(spk_PhsCohSpec2_std);
    res(2).spk_CohSpecChAll_SE = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
    res(2).spk_CohSpecChAll_SE(:,resV1.Fgoodunit) = sqrt(spk_CohSpecCh2_std);
    res(2).spk_PhsCohSpecChAll_SE = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
    res(2).spk_PhsCohSpecChAll_SE(:,resV1.Fgoodunit) = sqrt(spk_PhsCohSpecCh2_std);
    
    if ~isempty(resV1.gain)
        res(2).traj = interp1((1:numel(resV1.traj))*samplingdt,resV1.traj,res(2).t,'nearest');
        res(2).smthBallSpd = interp1((1:numel(resV1.smthBallSpd))*samplingdt,resV1.smthBallSpd,res(2).t,'nearest');
        res(2).outcome = interp1((1:numel(resV1.outcome))*samplingdt,resV1.outcome,res(2).t,'nearest');
        res(2).gain = interp1((1:numel(resV1.gain))*samplingdt,resV1.gain,res(2).t,'nearest');
        res(2).blanks = interp1((1:numel(resV1.blanks))*samplingdt,resV1.blanks,res(2).t,'nearest');
        
        for g = [2 1 3]
            tidx = ismember(res(2).outcome,outcomeval) & ismember(res(2).gain,gainval(g)) & res(2).smthBallSpd > speed_th & ~res(2).blanks;
            res(2).meanBallSpd{g} = nanmean(res(2).smthBallSpd(tidx));
        end
        for g = [2 1 3]
            tidx = ismember(res(2).outcome,outcomeval) & ismember(res(2).gain,gainval(g)) & res(2).smthBallSpd > speed_th & ~res(2).blanks;
            
            res(2).lfp_Spec{g} = squeeze(abs(nanmean(nanmean(lfp_S1_V1(tidx,:,:),1),3)));
            res(2).spk_SpecCh{g} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
            res(2).spk_SpecCh{g}(:,resV1.Fgoodunit) = squeeze(abs(nanmean(spk_S2_V1(tidx,:,:),1)));
            res(2).spk_Spec{g} = squeeze(abs(nanmean(nanmean(spk_S2_V1(tidx,:,:),1),3)));
            res(2).spk_CohSpec{g} = squeeze(abs(nanmean(nanmean(spk_S12_V1(tidx,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1_V1(tidx,:,:),1),3).*nanmean(nanmean(spk_S2_V1(tidx,:,:),1),3))));
            res(2).spk_PhsCohSpec{g} = squeeze(angle(nanmean(nanmean(spk_S12_V1(tidx,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1_V1(tidx,:,:),1),3).*nanmean(nanmean(spk_S2_V1(tidx,:,:),1),3))));
            res(2).spk_CohSpecCh{g} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
            res(2).spk_CohSpecCh{g}(:,resV1.Fgoodunit) = squeeze(abs(nanmean(spk_S12_V1(tidx,:,:),1)./sqrt(nanmean(lfp_S1_CA1_V1(tidx,:,:),1).*nanmean(spk_S2_V1(tidx,:,:),1))));
            res(2).spk_PhsCohSpecCh{g} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
            res(2).spk_PhsCohSpecCh{g}(:,resV1.Fgoodunit) = squeeze(angle(nanmean(spk_S12_V1(tidx,:,:),1)./sqrt(nanmean(lfp_S1_CA1_V1(tidx,:,:),1).*nanmean(spk_S2_V1(tidx,:,:),1))));
            
            idx = find(tidx);
            [CVO] = crossValPartition(1:numel(idx), nfold, false, floor(60/movingwin(2)));
            lfp_Spec2_std = 0;
            spk_SpecCh2_std = 0;
            spk_Spec2_std = 0;
            spk_CohSpec2_std = 0;
            spk_PhsCohSpec2_std = 0;
            spk_CohSpecCh2_std = 0;
            spk_PhsCohSpecCh2_std = 0;
            for kiter = 1:nfold
                tidxCVO = false(size(tidx));
                tidxCVO(idx(CVO.train{kiter})) = true;
                
                res(2).lfp_Spec_CVO{kiter,g} = squeeze(abs(nanmean(nanmean(lfp_S1_V1(tidxCVO,:,:),1),3)));
                lfp_Spec2_iter = res(2).lfp_Spec_CVO{kiter,g};
                lfp_Spec2_std = lfp_Spec2_std + (nfold - 1)/nfold*(lfp_Spec2_iter - res(2).lfp_Spec{g}).^2;
                
                res(2).spk_SpecCh_CVO{kiter,g} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
                res(2).spk_SpecCh_CVO{kiter,g}(:,resV1.Fgoodunit) = squeeze(abs(nanmean(spk_S2_V1(tidxCVO,:,:),1)));
                spk_SpecCh2_iter = squeeze(abs(nanmean(spk_S2_V1(tidxCVO,:,:),1)));
                spk_SpecCh2_std = spk_SpecCh2_std + (nfold - 1)/nfold*(spk_SpecCh2_iter - res(2).spk_SpecCh{g}(:,resV1.Fgoodunit)).^2;
                
                res(2).spk_Spec_CVO{kiter,g} = squeeze(abs(nanmean(nanmean(spk_S2_V1(tidxCVO,:,:),1),3)));
                spk_Spec2_iter = res(2).spk_Spec_CVO{kiter,g};
                spk_Spec2_std = spk_Spec2_std + (nfold - 1)/nfold*(spk_Spec2_iter - res(2).spk_Spec{g}).^2;
                
                res(2).spk_CohSpec_CVO{kiter,g} = squeeze(abs(nanmean(nanmean(spk_S12_V1(tidxCVO,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1_V1(tidxCVO,:,:),1),3).*nanmean(nanmean(spk_S2_V1(tidxCVO,:,:),1),3))));
                spk_CohSpec2_iter = res(2).spk_CohSpec_CVO{kiter,g};
                spk_CohSpec2_std = spk_CohSpec2_std + (nfold - 1)/nfold*(spk_CohSpec2_iter - res(2).spk_CohSpec{g}).^2;
                
                res(2).spk_PhsCohSpec_CVO{kiter,g} = squeeze(angle(nanmean(nanmean(spk_S12_V1(tidxCVO,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1_V1(tidxCVO,:,:),1),3).*nanmean(nanmean(spk_S2_V1(tidxCVO,:,:),1),3))));
                spk_PhsCohSpec2_iter = res(2).spk_PhsCohSpec_CVO{kiter,g};
                spk_PhsCohSpec2_std = spk_PhsCohSpec2_std + (nfold - 1)/nfold*circ_dist(spk_PhsCohSpec2_iter,res(2).spk_PhsCohSpec{g}).^2;
                
                res(2).spk_CohSpecCh_CVO{kiter,g} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
                res(2).spk_CohSpecCh_CVO{kiter,g}(:,resV1.Fgoodunit) = squeeze(abs(nanmean(spk_S12_V1(tidxCVO,:,:),1)./sqrt(nanmean(lfp_S1_CA1_V1(tidxCVO,:,:),1).*nanmean(spk_S2_V1(tidxCVO,:,:),1))));
                spk_CohSpecCh2_iter = res(2).spk_CohSpecCh_CVO{kiter,g}(:,resV1.Fgoodunit);
                spk_CohSpecCh2_std = spk_CohSpecCh2_std + (nfold - 1)/nfold*(spk_CohSpecCh2_iter - res(2).spk_CohSpecCh{g}(:,resV1.Fgoodunit)).^2;
                
                res(2).spk_PhsCohSpecCh_CVO{kiter,g} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
                res(2).spk_PhsCohSpecCh_CVO{kiter,g}(:,resV1.Fgoodunit) = squeeze(angle(nanmean(spk_S12_V1(tidxCVO,:,:),1)./sqrt(nanmean(lfp_S1_CA1_V1(tidxCVO,:,:),1).*nanmean(spk_S2_V1(tidxCVO,:,:),1))));
                spk_PhsCohSpecCh2_iter = res(2).spk_PhsCohSpecCh_CVO{kiter,g}(:,resV1.Fgoodunit);
                spk_PhsCohSpecCh2_std = spk_PhsCohSpecCh2_std + (nfold - 1)/nfold*circ_dist(spk_PhsCohSpecCh2_iter,res(2).spk_PhsCohSpecCh{g}(:,resV1.Fgoodunit)).^2;
            end
            
            res(2).lfp_Spec_SE{g} = sqrt(lfp_Spec2_std);
            res(2).spk_SpecCh_SE{g} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
            res(2).spk_SpecCh_SE{g}(:,resV1.Fgoodunit) = sqrt(spk_SpecCh2_std);
            res(2).spk_Spec_SE{g} = sqrt(spk_Spec2_std);
            res(2).spk_CohSpec_SE{g} = sqrt(spk_CohSpec2_std);
            res(2).spk_PhsCohSpec_SE{g} = sqrt(spk_PhsCohSpec2_std);
            res(2).spk_CohSpecCh_SE{g} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
            res(2).spk_CohSpecCh_SE{g}(:,resV1.Fgoodunit) = sqrt(spk_CohSpecCh2_std);
            res(2).spk_PhsCohSpecCh_SE{g} = NaN(size(spk_S12_V1,2),numel(resV1.Fgoodunit));
            res(2).spk_PhsCohSpecCh_SE{g}(:,resV1.Fgoodunit) = sqrt(spk_PhsCohSpecCh2_std);
        end
    end
end
if ~isempty(resCA1)
    res(1).lfp_SpecAll = squeeze(abs(nanmean(nanmean(lfp_S1_CA1(:,:,:),1),3)));
    res(1).spk_SpecChAll = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
    res(1).spk_SpecChAll(:,resCA1.Fgoodunit) = squeeze(abs(nanmean(spk_S2_CA1(:,:,:),1)));
    res(1).spk_SpecAll = squeeze(abs(nanmean(nanmean(spk_S2_CA1(:,:,:),1),3)));
    res(1).spk_CohSpecAll = squeeze(abs(nanmean(nanmean(spk_S12_CA1(:,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1(:,:,:),1),3).*nanmean(nanmean(spk_S2_CA1(:,:,:),1),3))));
    res(1).spk_PhsCohSpecAll = squeeze(angle(nanmean(nanmean(spk_S12_CA1(:,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1(:,:,:),1),3).*nanmean(nanmean(spk_S2_CA1(:,:,:),1),3))));
    res(1).spk_CohSpecChAll = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
    res(1).spk_CohSpecChAll(:,resCA1.Fgoodunit) = squeeze(abs(nanmean(spk_S12_CA1(:,:,:),1)./sqrt(nanmean(lfp_S1_CA1(:,:,:),1).*nanmean(spk_S2_CA1(:,:,:),1))));
    res(1).spk_PhsCohSpecChAll = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
    res(1).spk_PhsCohSpecChAll(:,resCA1.Fgoodunit) = squeeze(angle(nanmean(spk_S12_CA1(:,:,:),1)./sqrt(nanmean(lfp_S1_CA1(:,:,:),1).*nanmean(spk_S2_CA1(:,:,:),1))));
    
    tidx = true(size(lfp_S1_CA1,1),1);
    idx = find(tidx);
    [CVO] = crossValPartition(1:numel(idx), nfold, false, floor(60/movingwin(2)));
    lfp_Spec1_std = 0;
    spk_SpecCh1_std = 0;
    spk_Spec1_std = 0;
    spk_CohSpec1_std = 0;
    spk_PhsCohSpec1_std = 0;
    spk_CohSpecCh1_std = 0;
    spk_PhsCohSpecCh1_std = 0;
    for kiter = 1:nfold
        tidxCVO = false(size(tidx));
        tidxCVO(idx(CVO.train{kiter})) = true;
        
        res(1).lfp_SpecAll_CVO{kiter} = squeeze(abs(nanmean(nanmean(lfp_S1_CA1(tidxCVO,:,:),1),3)));
        lfp_Spec1_iter = res(1).lfp_SpecAll_CVO{kiter};
        lfp_Spec1_std = lfp_Spec1_std + (nfold - 1)/nfold*(lfp_Spec1_iter - res(1).lfp_SpecAll).^2;
        
        res(1).spk_SpecChAll_CVO{kiter} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
        res(1).spk_SpecChAll_CVO{kiter}(:,resCA1.Fgoodunit) = squeeze(abs(nanmean(spk_S2_CA1(tidxCVO,:,:),1)));
        spk_SpecCh1_iter = squeeze(abs(nanmean(spk_S2_CA1(tidxCVO,:,:),1)));
        spk_SpecCh1_std = spk_SpecCh1_std + (nfold - 1)/nfold*(spk_SpecCh1_iter - res(1).spk_SpecChAll(:,resCA1.Fgoodunit)).^2;
        
        res(1).spk_SpecAll_CVO{kiter} = squeeze(abs(nanmean(nanmean(spk_S2_CA1(tidxCVO,:,:),1),3)));
        spk_Spec1_iter = res(1).spk_SpecAll_CVO{kiter};
        spk_Spec1_std = spk_Spec1_std + (nfold - 1)/nfold*(spk_Spec1_iter - res(1).spk_SpecAll).^2;
        
        res(1).spk_CohSpecAll_CVO{kiter} = squeeze(abs(nanmean(nanmean(spk_S12_CA1(tidxCVO,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1(tidxCVO,:,:),1),3).*nanmean(nanmean(spk_S2_CA1(tidxCVO,:,:),1),3))));
        spk_CohSpec1_iter = res(1).spk_CohSpecAll_CVO{kiter};
        spk_CohSpec1_std = spk_CohSpec1_std + (nfold - 1)/nfold*(spk_CohSpec1_iter - res(1).spk_CohSpecAll).^2;
        
        res(1).spk_PhsCohSpecAll_CVO{kiter} = squeeze(angle(nanmean(nanmean(spk_S12_CA1(tidxCVO,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1(tidxCVO,:,:),1),3).*nanmean(nanmean(spk_S2_CA1(tidxCVO,:,:),1),3))));
        spk_PhsCohSpec1_iter = res(1).spk_PhsCohSpecAll_CVO{kiter};
        spk_PhsCohSpec1_std = spk_PhsCohSpec1_std + (nfold - 1)/nfold*circ_dist(spk_PhsCohSpec1_iter,res(1).spk_PhsCohSpecAll).^2;
        
        res(1).spk_CohSpecChAll_CVO{kiter} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
        res(1).spk_CohSpecChAll_CVO{kiter}(:,resCA1.Fgoodunit) = squeeze(abs(nanmean(spk_S12_CA1(tidxCVO,:,:),1)./sqrt(nanmean(lfp_S1_CA1(tidxCVO,:,:),1).*nanmean(spk_S2_CA1(tidxCVO,:,:),1))));
        spk_CohSpecCh1_iter = res(1).spk_CohSpecChAll_CVO{kiter}(:,resCA1.Fgoodunit);
        spk_CohSpecCh1_std = spk_CohSpecCh1_std + (nfold - 1)/nfold*(spk_CohSpecCh1_iter - res(1).spk_CohSpecChAll(:,resCA1.Fgoodunit)).^2;
        
        res(1).spk_PhsCohSpecChAll_CVO{kiter} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
        res(1).spk_PhsCohSpecChAll_CVO{kiter}(:,resCA1.Fgoodunit) = squeeze(angle(nanmean(spk_S12_CA1(tidxCVO,:,:),1)./sqrt(nanmean(lfp_S1_CA1(tidxCVO,:,:),1).*nanmean(spk_S2_CA1(tidxCVO,:,:),1))));
        spk_PhsCohSpecCh1_iter = res(1).spk_PhsCohSpecChAll(:,resCA1.Fgoodunit);
        spk_PhsCohSpecCh1_std = spk_PhsCohSpecCh1_std + (nfold - 1)/nfold*circ_dist(spk_PhsCohSpecCh1_iter,res(1).spk_PhsCohSpecChAll(:,resCA1.Fgoodunit)).^2;
    end
    res(1).lfp_SpecAll_SE = sqrt(lfp_Spec1_std);
    res(1).spk_SpecChAll_SE = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
    res(1).spk_SpecChAll_SE(:,resCA1.Fgoodunit) = sqrt(spk_SpecCh1_std);
    res(1).spk_SpecAll_SE = sqrt(spk_Spec1_std);
    res(1).spk_CohSpecAll_SE = sqrt(spk_CohSpec1_std);
    res(1).spk_PhsCohSpecAll_SE = sqrt(spk_PhsCohSpec1_std);
    res(1).spk_CohSpecChAll_SE = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
    res(1).spk_CohSpecChAll_SE(:,resCA1.Fgoodunit) = sqrt(spk_CohSpecCh1_std);
    res(1).spk_PhsCohSpecChAll_SE = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
    res(1).spk_PhsCohSpecChAll_SE(:,resCA1.Fgoodunit) = sqrt(spk_PhsCohSpecCh1_std);
    
    if ~isempty(resCA1.gain)
        res(1).traj = interp1((1:numel(resCA1.traj))*samplingdt,resCA1.traj,res(1).t,'nearest');
        res(1).smthBallSpd = interp1((1:numel(resCA1.smthBallSpd))*samplingdt,resCA1.smthBallSpd,res(1).t,'nearest');
        res(1).outcome = interp1((1:numel(resCA1.outcome))*samplingdt,resCA1.outcome,res(1).t,'nearest');
        res(1).gain = interp1((1:numel(resCA1.gain))*samplingdt,resCA1.gain,res(1).t,'nearest');
        res(1).blanks = interp1((1:numel(resCA1.blanks))*samplingdt,resCA1.blanks,res(1).t,'nearest');
        
        for g = [2 1 3]
            tidx = ismember(res(1).outcome,outcomeval) & ismember(res(1).gain,gainval(g)) & res(1).smthBallSpd > speed_th & ~res(1).blanks;
            res(1).meanBallSpd{g} = nanmean(res(1).smthBallSpd(tidx));
        end
        for g = [2 1 3]
            tidx = ismember(res(1).outcome,outcomeval) & ismember(res(1).gain,gainval(g)) & res(1).smthBallSpd > speed_th & ~res(1).blanks;
            
            res(1).lfp_Spec{g} = squeeze(abs(nanmean(nanmean(lfp_S1_CA1(tidx,:,:),1),3)));
            res(1).spk_SpecCh{g} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
            res(1).spk_SpecCh{g}(:,resCA1.Fgoodunit) = squeeze(abs(nanmean(spk_S2_CA1(tidx,:,:),1)));
            res(1).spk_Spec{g} = squeeze(abs(nanmean(nanmean(spk_S2_CA1(tidx,:,:),1),3)));
            res(1).spk_CohSpec{g} = squeeze(abs(nanmean(nanmean(spk_S12_CA1(tidx,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1(tidx,:,:),1),3).*nanmean(nanmean(spk_S2_CA1(tidx,:,:),1),3))));
            res(1).spk_PhsCohSpec{g} = squeeze(angle(nanmean(nanmean(spk_S12_CA1(tidx,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1(tidx,:,:),1),3).*nanmean(nanmean(spk_S2_CA1(tidx,:,:),1),3))));
            res(1).spk_CohSpecCh{g} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
            res(1).spk_CohSpecCh{g}(:,resCA1.Fgoodunit) = squeeze(abs(nanmean(spk_S12_CA1(tidx,:,:),1)./sqrt(nanmean(lfp_S1_CA1(tidx,:,:),1).*nanmean(spk_S2_CA1(tidx,:,:),1))));
            res(1).spk_PhsCohSpecCh{g} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
            res(1).spk_PhsCohSpecCh{g}(:,resCA1.Fgoodunit) = squeeze(angle(nanmean(spk_S12_CA1(tidx,:,:),1)./sqrt(nanmean(lfp_S1_CA1(tidx,:,:),1).*nanmean(spk_S2_CA1(tidx,:,:),1))));
            
            idx = find(tidx);
            [CVO] = crossValPartition(1:numel(idx), nfold, false, floor(60/movingwin(2)));
            lfp_Spec1_std = 0;
            spk_SpecCh1_std = 0;
            spk_Spec1_std = 0;
            spk_CohSpec1_std = 0;
            spk_PhsCohSpec1_std = 0;
            spk_CohSpecCh1_std = 0;
            spk_PhsCohSpecCh1_std = 0;
            for kiter = 1:nfold
                tidxCVO = false(size(tidx));
                tidxCVO(idx(CVO.train{kiter})) = true;
                
                lfp_Spec1_iter = squeeze(abs(nanmean(nanmean(lfp_S1_CA1(tidxCVO,:,:),1),3)));
                lfp_Spec1_std = lfp_Spec1_std + (nfold - 1)/nfold*(lfp_Spec1_iter - res(1).lfp_Spec{g}).^2;
                
                res(1).spk_SpecCh_CVO{kiter,g} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
                res(1).spk_SpecCh_CVO{kiter,g}(:,resCA1.Fgoodunit) = squeeze(abs(nanmean(spk_S2_CA1(tidxCVO,:,:),1)));
                spk_SpecCh1_iter = res(1).spk_SpecCh_CVO{kiter,g}(:,resCA1.Fgoodunit);
                spk_SpecCh1_std = spk_SpecCh1_std + (nfold - 1)/nfold*(spk_SpecCh1_iter - res(1).spk_SpecCh{g}(:,resCA1.Fgoodunit)).^2;
                
                res(1).spk_Spec_CVO{kiter,g} = squeeze(abs(nanmean(nanmean(spk_S2_CA1(tidxCVO,:,:),1),3)));
                spk_Spec1_iter = res(1).spk_Spec_CVO{kiter,g};
                spk_Spec1_std = spk_Spec1_std + (nfold - 1)/nfold*(spk_Spec1_iter - res(1).spk_Spec{g}).^2;
                
                res(1).spk_CohSpec_CVO{kiter,g} = squeeze(abs(nanmean(nanmean(spk_S12_CA1(tidxCVO,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1(tidxCVO,:,:),1),3).*nanmean(nanmean(spk_S2_CA1(tidxCVO,:,:),1),3))));
                spk_CohSpec1_iter = res(1).spk_CohSpec_CVO{kiter,g};
                spk_CohSpec1_std = spk_CohSpec1_std + (nfold - 1)/nfold*(spk_CohSpec1_iter - res(1).spk_CohSpec{g}).^2;
                
                res(1).spk_PhsCohSpec_CVO{kiter,g} = squeeze(angle(nanmean(nanmean(spk_S12_CA1(tidxCVO,:,:),1),3)./sqrt(nanmean(nanmean(lfp_S1_CA1(tidxCVO,:,:),1),3).*nanmean(nanmean(spk_S2_CA1(tidxCVO,:,:),1),3))));
                spk_PhsCohSpec1_iter = res(1).spk_PhsCohSpec_CVO{kiter,g};
                spk_PhsCohSpec1_std = spk_PhsCohSpec1_std + (nfold - 1)/nfold*circ_dist(spk_PhsCohSpec1_iter,res(1).spk_PhsCohSpec{g}).^2;
                
                res(1).spk_CohSpecCh_CVO{kiter,g} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
                res(1).spk_CohSpecCh_CVO{kiter,g}(:,resCA1.Fgoodunit) = squeeze(abs(nanmean(spk_S12_CA1(tidxCVO,:,:),1)./sqrt(nanmean(lfp_S1_CA1(tidxCVO,:,:),1).*nanmean(spk_S2_CA1(tidxCVO,:,:),1))));
                spk_CohSpecCh1_iter = res(1).spk_CohSpecCh_CVO{kiter,g}(:,resCA1.Fgoodunit);
                spk_CohSpecCh1_std = spk_CohSpecCh1_std + (nfold - 1)/nfold*(spk_CohSpecCh1_iter - res(1).spk_CohSpecCh{g}(:,resCA1.Fgoodunit)).^2;
                
                res(1).spk_PhsCohSpecCh_CVO{kiter,g} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
                res(1).spk_PhsCohSpecCh_CVO{kiter,g}(:,resCA1.Fgoodunit) = squeeze(angle(nanmean(spk_S12_CA1(tidxCVO,:,:),1)./sqrt(nanmean(lfp_S1_CA1(tidxCVO,:,:),1).*nanmean(spk_S2_CA1(tidxCVO,:,:),1))));
                spk_PhsCohSpecCh1_iter = res(1).spk_PhsCohSpecCh_CVO{kiter,g}(:,resCA1.Fgoodunit);
                spk_PhsCohSpecCh1_std = spk_PhsCohSpecCh1_std + (nfold - 1)/nfold*circ_dist(spk_PhsCohSpecCh1_iter,res(1).spk_PhsCohSpecCh{g}(:,resCA1.Fgoodunit)).^2;
            end
            
            res(1).lfp_Spec_SE{g} = sqrt(lfp_Spec1_std);
            res(1).spk_SpecCh{g} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
            res(1).spk_SpecCh_SE{g}(:,resCA1.Fgoodunit) = sqrt(spk_SpecCh1_std);
            res(1).spk_Spec_SE{g} = sqrt(spk_Spec1_std);
            res(1).spk_CohSpec_SE{g} = sqrt(spk_CohSpec1_std);
            res(1).spk_PhsCohSpec_SE{g} = sqrt(spk_PhsCohSpec1_std);
            res(1).spk_CohSpecCh_SE{g} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
            res(1).spk_CohSpecCh_SE{g}(:,resCA1.Fgoodunit) = sqrt(spk_CohSpecCh1_std);
            res(1).spk_PhsCohSpecCh_SE{g} = NaN(size(spk_S12_CA1,2),numel(resCA1.Fgoodunit));
            res(1).spk_PhsCohSpecCh_SE{g}(:,resCA1.Fgoodunit) = sqrt(spk_PhsCohSpecCh1_std);
        end
    end
end
if isfield(resV1, 'LFPfilt_VS')
    for iprotocol = 1:numel(resCA1.LFPfilt_VS)
        res(1).StimType{iprotocol} = resCA1.StimType{iprotocol};
        res(2).StimType{iprotocol} = resV1.StimType{iprotocol};
        res(1).repVS{iprotocol} = squeeze(mean(resCA1.LFPfilt_VS{iprotocol}(:,:,:),1));
        res(2).repVS{iprotocol} = squeeze(mean(resV1.LFPfilt_VS{iprotocol}(:,:,:),1));
    end
else
    res(1).repVS = [];
    res(2).repVS = [];
end
end