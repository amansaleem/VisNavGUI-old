list = [...
    110519 3;
    110529 3;
    110530 3;
    110720 8;
    110725 3;
    110802 4;
    110803 3;
    110830 3;
    110831 3;
    120713 3;
    120727 8;
    120801 4];
%     120712 5;
for num = [1:length(list)-1]
    
    if ismac
        fileName = ['/Users/aman/Dropbox/Work/Data/E-I-LFP-fromBH/Awake_Cclamp_Vm/' 'M' num2str(list(num,1)) '_' num2str(list(num,2)) '_filtered'];
    else
        fileName = ['C:\Users\aman\Dropbox\Work\Data\E-I-LFP-fromBH\Awake_Cclamp_Vm\' 'M' num2str(list(num,1)) '_' num2str(list(num,2)) '_filtered'];
    end
    load(fileName);
    name = ['M' num2str(list(num,1)) ' ' num2str(list(num,2))];
    
    params.Fs = 1000;
    params.tapers = [3 2];
    
    params.fpass = [20 100];
    movingwin=[0.33 0.05]; %floor(size(temp.epsc.filtered{1},2)/200)/100
    
    clear iC iS1 iS2 iphi iS12 vm_ds lfp_ds vm lfp
    for istim = 1:length(temp.lfp.filtered)
        vm = (temp.vm.filtered{istim}');
        lfp = (temp.lfp.filtered{istim}');
        trialIdx = 1;
        for itrial = 1:size(vm,2);
            if ~isempty(temp.lfp.tracking)
                if sum((temp.lfp.tracking(:,1) - istim)==0 & (temp.lfp.tracking(:,2) - itrial)==0)>0
                    continue
                end
                if (sum(isnan(vm(:,itrial)))>0 | sum(isnan(lfp(:,itrial)))>0)
                    continue
                end
            end
            vm_ds(:,trialIdx) = decimate(vm(:,itrial),20);
            lfp_ds(:,trialIdx) = decimate(lfp(:,itrial),20);
            trialIdx = trialIdx + 1;
            
        end
        %     [iC(istim,:,:),iphi(istim,:,:),iS12(istim,:,:),iS1(istim,:,:),iS2(istim,:,:),t,f]=cohgramc(ilfp_ds(1:350,:),ipsc_ds(1:350,:),movingwin,params);
        [Ct,phit,S12t,S1t,S2t,t,f]=cohgramc(lfp_ds,vm_ds,movingwin,params);
        if length(t)<2
            iC(istim,:) = squeeze(mean(Ct,3));
            iphi(istim,:) = squeeze(mean(phit,3));
            iS12(istim,:) = squeeze(mean(S12t,3));
            iS1(istim,:) = squeeze(mean(S1t,3));
            iS2(istim,:) = squeeze(mean(S2t,3));
        else
            iC(istim,:) = squeeze(mean(mean(Ct,1),3));
            iphi(istim,:) = squeeze(mean(mean(phit,1),3));
            iS12(istim,:) = squeeze(mean(mean(S12t,1),3));
            iS1(istim,:) = squeeze(mean(mean(S1t,1),3));
            iS2(istim,:) = squeeze(mean(mean(S2t,1),3));
        end
    end
    
    %%
    solo_lfp(num,:) = squeeze(mean(iS1(1:end-1,:),1));
    solo_vm(num,:) = squeeze(mean(iS2(1:end-1,:),1));
    joint_vl(num,:) = squeeze(mean(iS12(1:end-1,:),1));
    coh_vl(num,:) = squeeze(mean(iC(1:end-1,:),1));
    
    joint_vl_spont(num,:) = squeeze(mean(mean(iS12(end,:,:),3),1));
    coh_vl_spont(num,:) = squeeze(mean(mean(iC(end,:,:),3),1));
    solo_lfp_spont(num,:) = squeeze(mean(mean(iS1(end,:,:),3),1));
    solo_vm_spont(num,:) = squeeze(mean(mean(iS2(end,:,:),3),1));
    
    f2 = fit(f', solo_lfp(num,:)', 'exp2');
    psolo_lfp(num,:) = 100*(-f2(f)+solo_lfp(num,:)')./f2(f);
    f2 = fit(f', solo_vm(num,:)', 'exp2');
    psolo_vm(num,:) = 100*(-f2(f)+solo_vm(num,:)')./f2(f);
    f2 = fit(f', solo_lfp_spont(num,:)', 'exp2');
    psolo_lfp_spont(num,:) = 100*(-f2(f)+solo_lfp_spont(num,:)')./f2(f);
    f2 = fit(f', solo_vm_spont(num,:)', 'exp2');
    psolo_vm_spont(num,:) = 100*(-f2(f)+solo_vm_spont(num,:)')./f2(f);
    
    %%
    figure(num);
    set(gcf, 'Position', [1   41  1920   1084]);
    subplot(231)
    plot(f, log(squeeze(mean(mean(iS1(end,:,:),3),1))),'ko-', 'linewidth', 1); hold on
    plot(f, log(squeeze(mean(mean(iS1(1:end-1,:,:),3),1))),'ko-', 'linewidth', 1.5); hold off
    title(['LFP power: ' name]);
    
    subplot(232)
    plot(f, log(f.*squeeze(mean(mean(iS2(end,:,:),3),1))),'ko-', 'linewidth', 1); hold on
    plot(f, log(f.*squeeze(mean(mean(iS2(1:end-1,:,:),3),1))),'ko-', 'linewidth', 1.5); hold off
    title(['V_m power, tapers: ' num2str(params.tapers)])
    
    subplot(233)
    plot(f, squeeze(mean(mean(real(iC(end,:,:)),3),1)),'ko-', 'linewidth', 1); hold on
    plot(f, squeeze(mean(mean(real(iC(1:end-1,:,:)),3),1)),'ko-', 'linewidth', 1.5); hold off
    title('Coherence')
    
    for n = 1:3
        subplot(2,3,n)
        axis tight;
        grid on;
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('Freq (Hz)')
        ylabel('log Power')
    end
    subplot(224)
    plot(f, psolo_vm_spont(num,:),'ko-', 'linewidth', 1); hold on
    plot(f, psolo_vm(num,:),'ko-', 'linewidth', 1.5); hold off
        axis tight;
        ylims = ylim;
    if sum(abs(ylims)>300)>0
        set(gca, 'YLim',[-300 300]);
    end
    grid on;
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('Freq (Hz)')
        ylabel('log Power')
    
        subplot(223)
    plot(f, psolo_lfp_spont(num,:),'ko-', 'linewidth', 1); hold on
    plot(f, psolo_lfp(num,:),'ko-', 'linewidth', 1.5); hold off
        axis tight;
    ylims = ylim;
    if sum(abs(ylims)>300)>0
        set(gca, 'YLim',[-300 300]);
    end
        grid on;
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('Freq (Hz)')
        ylabel('log Power')
    
        
    drawnow; % pause
end

%%
frange = f>52 & f<68;
fval = f(frange);
[gammaPow_lfp, peakPos] = max(psolo_lfp(:,frange),[],2);
peakVal_lfp = fval(peakPos');
[gammaPow_vm, peakPos] = max(psolo_vm(:,frange),[],2);
peakVal_vm = fval(peakPos');
[gammaPow_vm_spont, peakPos] = max(psolo_vm_spont(:,frange),[],2);
peakVal_vm_spont = fval(peakPos');
[gammaPow_lfp_spont, peakPos] = max(psolo_lfp_spont(:,frange),[],2);
peakVal_lfp_spont = fval(peakPos');