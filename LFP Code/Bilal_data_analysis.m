for num = 1:11 %15 12-15 are aneasthetised
    
    switch num
        case 1
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M120713_9_11_g_reversal_filtered_lfp_psc.mat')
            name = 'M120713 9 11';
        case 2
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M111004_12_13_g_reversal_filtered_lfp_psc.mat')
            name = 'M111004 12 13';
        case 3
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M120713_15_16_g_reversal_filtered_lfp_psc.mat')
            name = 'M120713 15 16';
        case 4
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M120802_5_6_g_reversal_filtered_lfp_psc.mat')
            name = 'M120802 5 6';
        case 5
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M111004_6_7_g_reversal_filtered_lfp_psc.mat')
            name = 'M111004 6 7';
        case 6
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M111003_13_14_g_reversal_filtered_lfp_psc.mat')
            name = 'M111003 13 14';
        case 7
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M110810_4_7_g_reversal_filtered_lfp_psc.mat')
            name = 'M110810 4 7';
        case 8
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M110804_15_17_g_reversal_filtered_lfp_psc.mat')
            name = 'M110804 15 17';
        case 9
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M110926_12_14_g_reversal_filtered_lfp_psc.mat')
            name = 'M110926 12 14';

% % %         case 12
% % %             load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M110927_6_7_g_reversal_filtered_lfp_psc.mat')
% % %             name = 'M110927 6 7';
        case 11
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M110927_10_11_g_reversal_filtered_lfp_psc.mat')
            name = 'M110927 10 11';
        case 10
            load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Vclamp_LFP_PSC\M120727_14_15_g_reversal_filtered_lfp_psc.mat')
            name = 'M120727 14 15';
        
%         case 10
%             load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Anesth_Vclamp_LFP_PSC\M120709_10_14_g_reversal_filtered_lfp_psc.mat')
%             name = 'aneas M111005_10_11';
%         case 11
%             load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Anesth_Vclamp_LFP_PSC\M111005_10_11_g_reversal_filtered_lfp_psc.mat')
%             name = 'aneas M111005_10_11';
%         case 12
%             load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Anesth_Vclamp_LFP_PSC\M110908_13_14_g_reversal_filtered_lfp_psc.mat')
%             name = 'aneas M110908_13_14';
%         case 13
%             load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Anesth_Vclamp_LFP_PSC\M110824_16_17_g_reversal_filtered_lfp_psc.mat')
%             name = 'aneas M110824_16_17';
%         case 14
%             load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Anesth_Vclamp_LFP_PSC\M110909_10_11_g_reversal_filtered_lfp_psc.mat')
%             name = 'aneas M110909_10_11';
%         case 15
%             load('C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Anesth_Vclamp_LFP_PSC\M110824_9_11_g_reversal_filtered_lfp_psc.mat')
%             name = 'aneas M110824_9_11';
%         
    end
    
    params.Fs = 1000;
    params.tapers = [3 2];
    
    params.fpass = [20 100];
    movingwin=[0.33 0.05]; %floor(size(temp.epsc.filtered{1},2)/200)/100
    
    clear C S1 S2 phi S12 epsc_ds elfp_ds
    for istim = 1:length(temp.lfp_e.filtered)
        epsc = (temp.epsc.filtered{istim}');
        elfp = (temp.lfp_e.filtered{istim}');
        trialIdx = 1;
        nBadTrials = 0;
        for itrial = 1:size(epsc,2);
            if ~isempty(temp.lfp_e.tracking)
                if sum((temp.lfp_e.tracking(:,1) - istim)==0 & (temp.lfp_e.tracking(:,2) - itrial)==0)>0
                    nBadTrials = nBadTrials + 1;
                    continue
                end
                if (sum(isnan(epsc(:,itrial)))>0 | sum(isnan(elfp(:,itrial)))>0)
                    nBadTrials = nBadTrials + 1;
                    continue
                end
            end
            epsc_ds(:,trialIdx) = decimate(epsc(:,itrial),20);
            elfp_ds(:,trialIdx) = decimate(elfp(:,itrial),20);
            trialIdx = trialIdx + 1;
        end
        %     [Ct(istim,:,:),phit(istim,:,:),S12t(istim,:,:),S1t(istim,:,:),S2t(istim,:,:),t,f]=cohgramc(elfp_ds(1:350,:),epsc_ds(1:350,:),movingwin,params);
        [Ct,phit,S12t,S1t,S2t,t,f]=cohgramc(elfp_ds,epsc_ds,movingwin,params);
        if nBadTrials>0
            Ct(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
            phit(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
            S12t(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
            S1t(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
            S2t(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
        end
        if length(t)<2
            C(istim,:,:) = Ct;
            phi(istim,:,:) = phit;
            S12(istim,:,:) = S12t;
            S1(istim,:,:) = S1t;
            S2(istim,:,:) = S2t;
        else
            C(istim,:,:) = squeeze(nanmean(Ct,1));
            phi(istim,:,:) = squeeze(nanmean(phit,1));
            S12(istim,:,:) = squeeze(nanmean(S12t,1));
            S1(istim,:,:) = squeeze(nanmean(S1t,1));
            S2(istim,:,:) = squeeze(nanmean(S2t,1));
        end
    end
    clear iC iS1 iS2 iphi iS12 ipsc_ds ilfp_ds
    for istim = 1:length(temp.lfp_i.filtered)
        ipsc = (temp.ipsc.filtered{istim}');
        ilfp = (temp.lfp_i.filtered{istim}');
        trialIdx = 1;
        nBadTrials = 0;
        for itrial = 1:size(ipsc,2);
            if ~isempty(temp.lfp_i.tracking)
                if sum((temp.lfp_i.tracking(:,1) - istim)==0 & (temp.lfp_i.tracking(:,2) - itrial)==0)>0
                    nBadTrials = nBadTrials + 1;
                    continue
                end
                if (sum(isnan(ipsc(:,itrial)))>0 | sum(isnan(ilfp(:,itrial)))>0)
                    nBadTrials = nBadTrials + 1;
                    continue
                end
            end
            ipsc_ds(:,trialIdx) = decimate(ipsc(:,itrial),20);
            ilfp_ds(:,trialIdx) = decimate(ilfp(:,itrial),20);
            trialIdx = trialIdx + 1;
        end
        %     [Ct(istim,:,:),phit(istim,:,:),S12t(istim,:,:),S1t(istim,:,:),S2t(istim,:,:),t,f]=cohgramc(elfp_ds(1:350,:),epsc_ds(1:350,:),movingwin,params);
        [iCt,iphit,iS12t,iS1t,iS2t,t,f]=cohgramc(ilfp_ds,ipsc_ds,movingwin,params);
        if nBadTrials>0
            iCt(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
            iphit(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
            iS12t(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
            iS1t(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
            iS2t(:,:,trialIdx:trialIdx+nBadTrials-1) = nan;
        end
        if length(t)<2
            iC(istim,:,:) = iCt;
            iphi(istim,:,:) = iphit;
            iS12(istim,:,:) = iS12t;
            iS1(istim,:,:) = iS1t;
            iS2(istim,:,:) = iS2t;
        else
            iC(istim,:,:) = squeeze(nanmean(iCt,1));
            iphi(istim,:,:) = squeeze(nanmean(iphit,1));
            iS12(istim,:,:) = squeeze(nanmean(iS12t,1));
            iS1(istim,:,:) = squeeze(nanmean(iS1t,1));
            iS2(istim,:,:) = squeeze(nanmean(iS2t,1));
        end
        
%         for itrial = 1:size(Ct,3)
%             f2 = fit(f', solo_el(num,:)', 'exp1');
%             psolo_el(num,:) = f2(f)-solo_el(num,:)';
%             f2 = fit(f', solo_ep(num,:)', 'exp1');
%             psolo_ep(num,:) = f2(f)-solo_ep(num,:)';
%             f2 = fit(f', solo_el_spont(num,:)', 'exp1');
%             psolo_el_spont(num,:) = f2(f)-solo_el_spont(num,:)';
%             f2 = fit(f', solo_ep_spont(num,:)', 'exp1');
%             psolo_ep_spont(num,:) = f2(f)-solo_ep_spont(num,:)';
%         end
%         for itrial = 1:size(iCt,3)
%             f2 = fit(f', solo_il(num,:)', 'exp1');
%             psolo_il(num,:) = f2(f)-solo_il(num,:)';
%             f2 = fit(f', solo_ip(num,:)', 'exp1');
%             psolo_ip(num,:) = f2(f)-solo_ip(num,:)';
%             f2 = fit(f', solo_il_spont(num,:)', 'exp1');
%             psolo_il_spont(num,:) = f2(f)-solo_il_spont(num,:)';
%             f2 = fit(f', solo_ip_spont(num,:)', 'exp1');
%             psolo_ip_spont(num,:) = f2(f)-solo_ip_spont(num,:)';
%         end
    end
    %%
    solo_el(num,:) = squeeze(nanmean(nanmean(S1(1:end-1,:,:),3),1));
    solo_ep(num,:) = squeeze(nanmean(nanmean(S2(1:end-1,:,:),3),1));
    solo_il(num,:) = squeeze(nanmean(nanmean(iS1(1:end-1,:,:),3),1));
    solo_ip(num,:) = squeeze(nanmean(nanmean(iS2(1:end-1,:,:),3),1));
    
    solo_el_spont(num,:) = squeeze(nanmean(nanmean(S1(end,:,:),3),1));
    solo_ep_spont(num,:) = squeeze(nanmean(nanmean(S2(end,:,:),3),1));
    solo_il_spont(num,:) = squeeze(nanmean(nanmean(iS1(end,:,:),3),1));
    solo_ip_spont(num,:) = squeeze(nanmean(nanmean(iS2(end,:,:),3),1));
    
    try
        all_el_spont(num,:,:) = squeeze(S1(end,:,:));
        all_ep_spont(num,:,:) = squeeze(S2(end,:,:));
    catch
        if size(all_el_spont,3)>size(S1,3)
            all_el_spont(num,:,1:size(S1,3)) = squeeze(S1(end,:,:));
            all_ep_spont(num,:,1:size(S1,3)) = squeeze(S2(end,:,:));
            all_el_spont(num,:,(size(S1,3)+1):end) = nan;
            all_ep_spont(num,:,(size(S1,3)+1):end) = nan;
        else size(all_el_spont,3)<size(S1,3)
            all_el_spont(:,:,end+1:size(S1,3)) = nan;
            all_ep_spont(:,:,end+1:size(S1,3)) = nan;
            all_el_spont(num,:,:) = squeeze(S1(end,:,:));
            all_ep_spont(num,:,:) = squeeze(S2(end,:,:));
        end
    end
    try
        %         if size(all_il_spont,3)==size(iS1,3) | num<1
        all_il_spont(num,:,:) = squeeze(iS1(end,:,:));
        all_ip_spont(num,:,:) = squeeze(iS2(end,:,:));
    catch
        if size(all_il_spont,3)>size(iS1,3)
            all_il_spont(num,:,1:size(iS1,3)) = squeeze(iS1(end,:,:));
            all_ip_spont(num,:,1:size(iS1,3)) = squeeze(iS2(end,:,:));
            all_il_spont(num,:,(size(iS1,3)+1):end) = nan;
            all_ip_spont(num,:,(size(iS1,3)+1):end) = nan;
        else size(all_il_spont,3)<size(iS1,3)
            all_il_spont(:,:,end+1:size(iS1,3)) = nan;
            all_ip_spont(:,:,end+1:size(iS1,3)) = nan;
            all_il_spont(num,:,:) = squeeze(iS1(end,:,:));
            all_ip_spont(num,:,:) = squeeze(iS2(end,:,:));
        end
    end
    joint_e(num,:) = squeeze(nanmean(nanmean(S12(1:end-1,:,:),3),1));
    coh_e(num,:) = squeeze(nanmean(nanmean(C(1:end-1,:,:),3),1));
    
    joint_i(num,:) = squeeze(nanmean(nanmean(iS12(1:end-1,:),1)));
    coh_i(num,:) = squeeze(nanmean(nanmean(iC(1:end-1,:),1)));
    
    joint_e_spont(num,:) = squeeze(nanmean(nanmean(S12(end,:,:),3),1));
    coh_e_spont(num,:) = squeeze(nanmean(nanmean(C(end,:,:),3),1));
    
    joint_i_spont(num,:) = squeeze(nanmean(nanmean(iS12(end,:,:),3),1));
    coh_i_spont(num,:) = squeeze(nanmean(nanmean(iC(end,:,:),3),1));
    
    t_f = (f<55 | f>70);
    f2 = fit(f(t_f)', solo_el(num,(t_f))', 'exp2');
    psolo_el(num,:) = 100*(f2(f)-solo_el(num,:)')./f2(f);
    f2 = fit(f(t_f)', solo_ep(num,(t_f))', 'exp2');
    psolo_ep(num,:) = 100*(f2(f)-solo_ep(num,:)')./f2(f);
    f2 = fit(f(t_f)', solo_el_spont(num,(t_f))', 'exp2');
    psolo_el_spont(num,:) = 100*(f2(f)-solo_el_spont(num,:)')./f2(f);
    f2 = fit(f(t_f)', solo_ep_spont(num,(t_f))', 'exp2');
    psolo_ep_spont(num,:) = 100*(f2(f)-solo_ep_spont(num,:)')./f2(f);
    
    f2 = fit(f(t_f)', solo_il(num,(t_f))', 'exp2');
    psolo_il(num,:) = 100*(f2(f)-solo_il(num,:)')./f2(f);
    f2 = fit(f(t_f)', solo_ip(num,(t_f))', 'exp2');
    psolo_ip(num,:) = 100*(f2(f)-solo_ip(num,:)')./f2(f);
    f2 = fit(f(t_f)', solo_il_spont(num,(t_f))', 'exp2');
    psolo_il_spont(num,:) = 100*(f2(f)-solo_il_spont(num,:)')./f2(f);
    f2 = fit(f(t_f)', solo_ip_spont(num,(t_f))', 'exp2');
    psolo_ip_spont(num,:) = 100*(f2(f)-solo_ip_spont(num,:)')./f2(f);
    
    %%
    figure(num+30);
    set(gcf, 'Position', [9          49        1082        1068]);
    subplot(231)
    plot(f, squeeze(nanmean(nanmean(log(S1(end,:,:)),3),1)),'ko-', 'linewidth', 1); hold on
    plot(f, squeeze(nanmean(nanmean(log(S1(1:end-1,:,:)),3),1)),'ko-', 'linewidth', 1.5);
    errorarea_as(f, squeeze(nanmean(nanmean(log(S1(end,:,:)),3),1)),nansem(squeeze(nanmean(log(S1(end,:,:)),1))'),'k')
    errorarea_as(f, squeeze(nanmean(nanmean(log(S1(1:end-1,:,:)),3),1)),nansem(squeeze(nanmean(log(S1(1:end-1,:,:)),1))'),'k')
    hold off;
    title(['LFP-E power: ' name]);
    
    subplot(232)
    plot((f), squeeze(nanmean(nanmean(log(S2(end,:,:)),3),1)),'ro-', 'linewidth', 1); hold on
    plot((f), squeeze(nanmean(nanmean(log(S2(1:end-1,:,:)),3))),'ro-', 'linewidth', 1.5);
    errorarea_as((f), squeeze(nanmean(nanmean(log(S2(end,:,:)),3),1)),nansem(squeeze(nanmean(log(S1(end,:,:)),1))'),[0.5 0 0])
    errorarea_as((f), squeeze(nanmean(nanmean(log(S2(1:end-1,:,:)),3),1)),nansem(squeeze(nanmean(log(S1(1:end-1,:,:)),1))'),'r')
    hold off;
    title(['EPSC power, tapers: ' num2str(params.tapers)])
    
    subplot(233)
    plot(f, squeeze(nanmean(nanmean(real(C(end,:,:)),3),1)),'ko-', 'linewidth', 1); hold on
    plot(f, squeeze(nanmean(nanmean(real(C(1:end-1,:,:)),3),1)),'ko-', 'linewidth', 1.5);
    errorarea_as(f, squeeze(nanmean(nanmean(real(C(end,:,:)),3),1)),nansem(squeeze(nanmean(real(C(end,:,:)),1))'),[0.5 0.5 0.5])
    errorarea_as(f, squeeze(nanmean(nanmean(real(C(1:end-1,:,:)),3),1)),nansem(squeeze(nanmean(real(C(1:end-1,:,:)),1))'),'k')
    title('Coherence power')
    
    % subplot(233)
    % plot(f, squeeze(nanmean(nanmean(real(S12(end,:,:)),3),1)),'ko-', 'linewidth', 1); hold on
    % plot(f, squeeze(nanmean(nanmean(real(S12(1:end-1,:,:)),3),1)),'ko-', 'linewidth', 1.5); hold off
    % title('Joint power')
    
    subplot(234)
    plot(f, squeeze(nanmean(nanmean(log(iS1(end,:,:)),3),1)),'ko-', 'linewidth', 1); hold on
    plot(f, squeeze(nanmean(nanmean(log(iS1(1:end-1,:,:)),3),1)),'ko-', 'linewidth', 1.5);
    errorarea_as(f, squeeze(nanmean(nanmean(log(iS1(end,:,:)),3),1)),nansem(squeeze(nanmean(log(iS1(end,:,:)),1))'),'k')
    errorarea_as(f, squeeze(nanmean(nanmean(log(iS1(1:end-1,:,:)),3),1)),nansem(squeeze(nanmean(log(iS1(1:end-1,:,:)),1))'),'k')
    hold off;
    title('LFP-I power')
    
    subplot(235)
    plot(f, squeeze(nanmean(nanmean(log(iS2(end,:,:)),3),1)),'bo-', 'linewidth', 1); hold on
    plot(f, squeeze(nanmean(nanmean(log(iS2(1:end-1,:,:)),3))),'bo-', 'linewidth', 1.5);
    errorarea_as(f, squeeze(nanmean(nanmean(log(iS2(end,:,:)),3),1)),nansem(squeeze(nanmean(log(iS1(end,:,:)),1))'),[0 0 0.5])
    errorarea_as(f, squeeze(nanmean(nanmean(log(iS2(1:end-1,:,:)),3),1)),nansem(squeeze(nanmean(log(iS1(1:end-1,:,:)),1))'),'b')
    hold off;
    title('IPSC power')
    
    subplot(236)
    plot(f, squeeze(nanmean(nanmean(real(iC(end,:,:)),3),1)),'ko-', 'linewidth', 1); hold on
    plot(f, squeeze(nanmean(nanmean(real(iC(1:end-1,:,:)),3),1)),'ko-', 'linewidth', 1.5)
    errorarea_as(f, squeeze(nanmean(nanmean(real(iC(end,:,:)),3),1)),nansem(squeeze(nanmean(real(iC(end,:,:)),1))'),[0.5 0.5 0.5])
    errorarea_as(f, squeeze(nanmean(nanmean(real(iC(1:end-1,:,:)),3),1)),nansem(squeeze(nanmean(real(iC(1:end-1,:,:)),1))'),'k')
    title('Coherence power')
    
    % subplot(237)
    % plot(f, squeeze(nanmean(nanmean(real(iS12(end,:,:)),3),1)),'ko-', 'linewidth', 1); hold on
    % plot(f, squeeze(nanmean(nanmean(real(iS12(1:end-1,:,:)),3),1)),'ko-', 'linewidth', 1.5); hold off
    % title('Joint power')
    
    for n = 1:6
        subplot(2,3,n)
        axis tight;
        grid on;
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('Freq (Hz)')
        ylabel('log Power')
    end
    drawnow; % pause
    
    %%
    figure(num+50);
    set(gcf, 'Position', [969    49   944   488]);
    subplot(221)
    plot(f, -psolo_el_spont(num,:),'ko-', 'linewidth', 1); hold on
    plot(f, -psolo_el(num,:),'ko-', 'linewidth', 1.5);
    hold off;
    title(['LFP-E power: ' name]);
    axis tight
    ylims = ylim;
    if sum(abs(ylims)>100)>0
        set(gca, 'YLim',[-100 100]);
    end
    
    subplot(222)
    plot(f, -psolo_ep_spont(num,:),'ro-', 'linewidth', 1); hold on
    plot(f, -psolo_ep(num,:),'ro-', 'linewidth', 1.5);
    hold off;
    title(['EPSC power, tapers: ' num2str(params.tapers)])
    axis tight
    ylims = ylim;
    if sum(abs(ylims)>100)>0
        set(gca, 'YLim',[-100 100]);
    end
    elims = get(gca,'ylim');
    
    subplot(223)
    plot(f, -psolo_el_spont(num,:),'ko-', 'linewidth', 1); hold on
    plot(f, -psolo_el(num,:),'ko-', 'linewidth', 1.5);
    hold off;
    title('LFP-I power')
    axis tight
    ylims = ylim;
    if sum(abs(ylims)>300)>0
        set(gca, 'YLim',[-300 300]);
    end
    
    subplot(224)
    plot(f, -psolo_ip_spont(num,:),'bo-', 'linewidth', 1); hold on
    plot(f, -psolo_ip(num,:),'bo-', 'linewidth', 1.5);
    hold off;
    title('IPSC power')
    axis tight
    ylims = ylim;
    if sum(abs(ylims)>300)>0
        set(gca, 'YLim',[-300 300]);
    end
    ilims = get(gca,'ylim');
    
    for n = 1:4
        subplot(2,2,n)
        
        grid on;
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('Freq (Hz)')
        ylabel('Power - exp fit')
    end
    ielims = [min([ilims elims]) max([ilims elims])];
    subplot(222)
    set(gca, 'ylim', ielims)
    subplot(224)
    set(gca, 'ylim', ielims)
    
    drawnow; 
%     pause
    
end

%%
frange = f>55 & f<70;
fval = f(frange);
[gammaPow_el, peakPos] = max(-psolo_el(:,frange),[],2);
peakVal_el = fval(peakPos');
[gammaPow_ep, peakPos] = max(-psolo_ep(:,frange),[],2);
peakVal_ep = fval(peakPos');
[gammaPow_ep_spont, peakPos] = max(-psolo_ep_spont(:,frange),[],2);
peakVal_ep_spont = fval(peakPos');
[gammaPow_el_spont, peakPos] = max(-psolo_el_spont(:,frange),[],2);
peakVal_el_spont = fval(peakPos');

[gammaPow_il, peakPos] = max(-psolo_il(:,frange),[],2);
peakVal_il = fval(peakPos');
[gammaPow_ip, peakPos] = max(-psolo_ip(:,frange),[],2);
peakVal_ip = fval(peakPos');
[gammaPow_ip_spont, peakPos] = max(-psolo_ip_spont(:,frange),[],2);
peakVal_ip_spont = fval(peakPos');
[gammaPow_il_spont, peakPos] = max(-psolo_il_spont(:,frange),[],2);
peakVal_il_spont = fval(peakPos');