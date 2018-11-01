%% 
% load Gamma_contrast_data_all
close all;
list_all(1).a = 'es_816_825_6';
list_all(2).a = 'es_815_819_6';
list_all(3).a = 'es_815_823_13';
list_all(4).a = 'es_719_831_6';
list_all(5).a = 'es_804_1_7';
list_all(6).a = 'es_804_3_6';
list_all(7).a = 'es_804_2_6';
clear frequency_pow 
alpha = 0.1;

bLow = 20;
bHigh= 45;

bLow2 = 70;
bHigh2= 90;

nLow = 57;
nHigh= 67;

for expIdx = 1:length(list_all)
    es_LFP = eval(list_all(expIdx).a);
    switch expIdx
        case 1
            [stimTimes,p] = createStimMatrix...
                ('M110816_BALL', 825, 6);
        case 2
            [stimTimes,p] = createStimMatrix...
                ('M110815_BALL', 819, 6);
        case 3
            [stimTimes,p] = createStimMatrix...
                ('M110815_BALL', 823, 13);
        case 4
            [stimTimes,p] = createStimMatrix...
                ('M110719_BALL', 831, 6);
        case 5
            [stimTimes,p] = createStimMatrix...
                ('M110804_BALL', 1, 7);
        case 6
            [stimTimes,p] = createStimMatrix...
                ('M110804_BALL', 2, 6);
        case 7
            [stimTimes,p] = createStimMatrix...
                ('M110804_BALL', 3, 6);
    end    
    [out] = getStimLFP(es_LFP, p, stimTimes, 1); 
    % es_LFP = es;
    figure;
    set(gcf, 'Name', list_all(expIdx).a)
    f = es_LFP.freq>20 & es_LFP.freq<90;
    list = [1+6:8; 9+6:16; 17+6:24; 25+6:32; 33+6:40];
    times = (-0.9:0.2:2.4)';
    
    subplot(2,3,1)
    imagesc(times, es_LFP.freq(f), squeeze(log(nanmean(out(41,:,f),1)))');
%     imagesc(times, es_LFP.freq(f),  - squeeze(log(nanmean(out(list(5),:,f),1)))'))
    istim = 1;
    frequency_pow(expIdx,istim, :,:) = (squeeze(log(nanmean(out(41,:,:),1)))');
    frequency_pow_diff(expIdx,istim, :,:) = sq(frequency_pow(expIdx,istim, :,:))...
        - repmat(mean(sq(frequency_pow(expIdx,istim,:,times>=-0.5 & times<0)),2),1,length(times));
    imagesc(times, es_LFP.freq(f), squeeze(frequency_pow_diff(expIdx,istim,f,:)));
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    title('0')
    subplot(2,3,2)
    imagesc(times, es_LFP.freq(f), squeeze(log(nanmean(out(list(1),:,f),1)))');
    istim = 2;
    frequency_pow(expIdx,istim, :,:) = (squeeze(log(nanmean(out(list(1),:,:),1)))');
    frequency_pow_diff(expIdx,istim, :,:) = sq(frequency_pow(expIdx,istim, :,:))...
        - repmat(mean(sq(frequency_pow(expIdx,istim,:,times>=-0.5 & times<0)),2),1,length(times));
    imagesc(times, es_LFP.freq(f), squeeze(frequency_pow_diff(expIdx,istim,f,:)));
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    title('10')
    subplot(2,3,3)
    imagesc(times, es_LFP.freq(f), squeeze(log(nanmean(out(list(2),:,f),1)))');
    istim = 3;
    frequency_pow(expIdx,istim, :,:) = (squeeze(log(nanmean(out(list(2),:,:),1)))');
    frequency_pow_diff(expIdx,istim, :,:) = sq(frequency_pow(expIdx,istim, :,:))...
        - repmat(mean(sq(frequency_pow(expIdx,istim,:,times>=-0.5 & times<0)),2),1,length(times));
    imagesc(times, es_LFP.freq(f), squeeze(frequency_pow_diff(expIdx,istim,f,:)));
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    title('25')
    subplot(2,3,4)
    imagesc(times, es_LFP.freq(f), squeeze(log(nanmean(out(list(3),:,f),1)))');
    istim = 4;
    frequency_pow(expIdx,istim, :,:) = (squeeze(log(nanmean(out(list(3),:,:),1)))');
    frequency_pow_diff(expIdx,istim, :,:) = sq(frequency_pow(expIdx,istim, :,:))...
        - repmat(mean(sq(frequency_pow(expIdx,istim,:,times>=-0.5 & times<0)),2),1,length(times));
    imagesc(times, es_LFP.freq(f), squeeze(frequency_pow_diff(expIdx,istim,f,:)));
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    title('50')
    subplot(2,3,5)
    imagesc(times, es_LFP.freq(f), squeeze(log(nanmean(out(list(4),:,f),1)))');
    istim = 5;
    frequency_pow(expIdx,istim, :,:) = (squeeze(log(nanmean(out(list(4),:,:),1)))');
    frequency_pow_diff(expIdx,istim, :,:) = sq(frequency_pow(expIdx,istim, :,:))...
        - repmat(mean(sq(frequency_pow(expIdx,istim,:,times>=-0.5 & times<0)),2),1,length(times));
    imagesc(times, es_LFP.freq(f), squeeze(frequency_pow_diff(expIdx,istim,f,:)));
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    title('75')
    subplot(2,3,6)
    imagesc(times, es_LFP.freq(f), squeeze(log(nanmean(out(list(5),:,f),1)))');
    istim = 6;
    frequency_pow(expIdx,istim, :,:) = (squeeze(log(nanmean(out(list(5),:,:),1)))');
    frequency_pow_diff(expIdx,istim, :,:) = sq(frequency_pow(expIdx,istim, :,:))...
        - repmat(mean(sq(frequency_pow(expIdx,istim,:,times>=-0.5 & times<0)),2),1,length(times));
    imagesc(times, es_LFP.freq(f), squeeze(frequency_pow(expIdx,istim,f,:)));
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    title('100')
    for n2 = 1:6
        subplot(2,3,n2)
%         set(gca, 'Clim', [-0.4 0.4]); RedWhiteBlue;
        line([0 0], ylim, 'linestyle','--','color','k')
        line([1 1], ylim, 'linestyle','--','color','k')
    end
    figure(101);
    subplot(3,3,expIdx)
%     imagesc(times, es_LFP.freq(f), -(squeeze(log(nanmean(out(41,:,f),1)))' - squeeze(log(nanmean(out(list(5),:,f),1)))'))
    imagesc(times, es_LFP.freq(f), squeeze(frequency_pow_diff(expIdx,6,f,:)))
    set(gca, 'Clim', [-0.4 0.4]); RedWhiteBlue;
    line([0 0], ylim, 'linestyle','--','color','k')
    line([1 1], ylim, 'linestyle','--','color','k')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    %     title(list_all(n).a
    title(['Experiment :' num2str(expIdx)]);
    
    figure(102);
    subplot(3,3,expIdx)
%     imagesc(times, es_LFP.freq(f), -(squeeze(log(nanmean(out(41,:,f),1)))' - squeeze(log(nanmean(out(list(5),:,f),1)))'))
    imagesc(times, es_LFP.freq(f), squeeze(frequency_pow_diff(expIdx,1,f,:)))
    set(gca, 'Clim', [-0.4 0.4]); RedWhiteBlue;
    line([0 0], ylim, 'linestyle','--','color','k')
    line([1 1], ylim, 'linestyle','--','color','k')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    %     title(list_all(n).a
    title(['Experiment :' num2str(expIdx)]);
end
figure(101);
subplot(3,3,9)
imagesc(times, es_LFP.freq(f), squeeze(nanmean(frequency_pow_diff(:,6,f,:),1)))
% set(gca, 'Clim', [-0.4 0.4]); RedWhiteBlue;
line([0 0], ylim, 'linestyle','--','color','k')
line([1 1], ylim, 'linestyle','--','color','k')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
%     title(list_all(n).a
title(['Average']);

figure(102);
subplot(3,3,9)
imagesc(times, es_LFP.freq(f), squeeze(nanmean(frequency_pow_diff(:,1,f,:),1)))
set(gca, 'Clim', [-0.4 0.4]); RedWhiteBlue;
line([0 0], ylim, 'linestyle','--','color','k')
line([1 1], ylim, 'linestyle','--','color','k')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
%     title(list_all(n).a
title(['Average']);

clear frequency_diff 
% for expIdx = 1:7
%     for n = 2:6
%         frequency_diff(expIdx,n-1,:,:) = frequency_pow(expIdx,n, :,:) - frequency_pow(expIdx,1, :,:);
%     end
% end
% 
% figure(101);
% subplot(3,2,6)
% imagesc(times, es_LFP.freq(f), (squeeze(nanmean(nanmean(frequency_diff(:,:,f,:),2),1))));
% line([0 0], ylim, 'linestyle','--','color','k')
% line([1 1], ylim, 'linestyle','--','color','k')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('Average across expt');
% xlabel('Time from stimulus onset');
% colorbar;

%%
f_b = (es_LFP.freq>=bLow & es_LFP.freq<=bHigh) | (es_LFP.freq>=bLow2 & es_LFP.freq<=bHigh2);
f_n = es_LFP.freq>=nLow & es_LFP.freq<=nHigh;
broad_band_pow  = squeeze(nanmean(frequency_pow(:,:,f_b,:),3));
narrow_band_pow = squeeze(nanmean(frequency_pow(:,:,f_n,:),3));

pre_freq_pow_average = squeeze((nanmean(frequency_pow(:,:,:,times>-0.5 & times<0),4)));

pre_broad_band_pow = squeeze(nanmean(pre_freq_pow_average(:,:,f_b),3));
pre_narrow_band_pow = squeeze(nanmean(pre_freq_pow_average(:,:,f_n),3));

broad_freq_pow_diff  = broad_band_pow - repmat(pre_broad_band_pow,[1,1,size(broad_band_pow,3)]);
narrow_freq_pow_diff = narrow_band_pow - repmat(pre_narrow_band_pow,[1,1,size(narrow_band_pow,3)]);

colour = [0.5 0.5 0.5;
        0.00 0 1.00;
        0.25 0 0.75;
        0.50 0 0.50;
        0.75 0 0.25;
        1.00 0 0.00];
figure(500)
subplot(131)
% f = es_LFP.freq>5 & es_LFP.freq<90;
ax = gca;
hold on;
set(ax, 'ColorOrder', colour);
plot(times, sq(nanmean(broad_freq_pow_diff,1))','linewidth',2);
% hold off;
for iStim = 1:size(broad_freq_pow_diff,2)
    X = times;
    Y = sq(nanmean(broad_freq_pow_diff(:,iStim,:),1))';
    E = sq(nansem(broad_freq_pow_diff(:,iStim,:),[],1))';
    switch iStim
        case 1
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 2
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 3
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 4
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 5
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 6
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
    end    
            
    hold on;
end

set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis tight
line(xlim, [0 0], 'linestyle','--','color','k')
xlabel('Time from stimulus onset (s)')
legend('0% contrast', '10%','25%','50%', '75%', '100%')
ylabel('\Delta log power')
title('Broad band gamma power (30-40Hz)')
line([0 0], ylim, 'linestyle','--','color','k')
line([1 1], ylim, 'linestyle','--','color','k')

subplot(132)
% f = es_LFP.freq>5 & es_LFP.freq<90;
% plot(times, sq(nanmean(narrow_freq_pow_diff,1))','.-');
for iStim = 1:size(broad_freq_pow_diff,2)
    X = times';
    Y = sq(nanmean(narrow_freq_pow_diff(:,iStim,:),1))';
    E = sq(nansem(narrow_freq_pow_diff(:,iStim,:),[],1))';
    switch iStim
        case 1
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 2
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 3
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 4
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 5
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
        case 6
            errorarea_as(X, Y, E , colour(iStim,:), alpha);
    end    
            
    hold on;
end

set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis tight
line(xlim, [0 0], 'linestyle','--','color','k')
xlabel('Time from stimulus onset (s)')
% legend('0% contrast', '10%','25%','50%', '75%', '100%')
ylabel('\Delta log power')
title('Narrow band gamma power (59-62Hz)')
line([0 0], ylim, 'linestyle','--','color','k')
line([1 1], ylim, 'linestyle','--','color','k')

contrasts = [0 10 25 50 75 100];


% subplot(121)
% plot(contrasts, broad_band_pow-pre_broad_band_pow, 'ko-', ...
%     contrasts, narrow_band_pow-pre_narrow_band_pow, 'ro-', ...
%     'linewidth',2);
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis tight
% xlabel('Contrast %')
% legend('Broad band (30-40 Hz)', 'Narrow-band (60-63 Hz)')
% line(xlim, narrow_band_pow(1).*[1 1],'color','r','linestyle','--')
% ylabel('mean log power')

%%
frequency_pow1 = sq(mean(frequency_pow(:,:,:,times>0.8 & times<1.1),4));
% frequency_diff1 = sq(mean(frequency_diff(:,:,:,times>0.8 & times<1.1),4));

% f = es_LFP.freq>51 & es_LFP.freq<62;
% freq_considered = es_LFP.freq(f);
% for n = 1:6
%     [~,fp_pos(n)] = min(sq(frequency_diff1(n,5,f)));
% end
% nbg_peak = freq_considered(fp_pos);

clear narrow_band_pow broad_band_pow

% narrow_band_pow = sq(mean(frequency_pow1(:,1:end,es_LFP.freq>=59 & es_LFP.freq<=62),3));
% broad_band_pow  = sq(mean(frequency_pow1(:,1:end,es_LFP.freq>=30 & es_LFP.freq<=40),3));


broad_band_pow   = sq(mean(broad_freq_pow_diff(:, 1:end,times>0.8 & times<1.1),3));
narrow_band_pow  = sq(mean(narrow_freq_pow_diff(:,1:end,times>0.8 & times<1.1),3));

subplot(233)
Xn = contrasts;
Yn = mean(narrow_band_pow);
En = sem(narrow_band_pow);
errorbar(Xn,Yn,En,'b')
hold on;
Xb = contrasts;
Yb = mean(broad_band_pow);
Eb = sem(broad_band_pow);
errorbar(Xb,Yb,Eb,'r')
axis([-1 101 -0.3 0.7])
hold off;
line(xlim,[0 0],'linestyle','--','color','k')
% plot(contrasts, mean(broad_band_pow),'b')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlabel('Contrast')
ylabel('\Delta Power (Stim - Pre stim)')
legend('Narrow-band \gamma (0.5-1s)', 'Broad-band \gamma (at 1s)')
% set(gca,'XScale','log')

% Single example 100% contrast
subplot(2,3,1)
imagesc(times, es_LFP.freq(f), squeeze(nanmean(frequency_pow_diff(3,6,f,:),1)))
title(['Example, 100% contrast']);

% Average n = x, 100% contrast
subplot(2,3,2)
imagesc(times, es_LFP.freq(f), squeeze(nanmean(frequency_pow_diff(:,6,f,:),1)))
title(['Average (n=7), 100% contrast']);

% Single example 50% contrast
subplot(2,3,4)
imagesc(times, es_LFP.freq(f), squeeze(nanmean(frequency_pow_diff(3,4,f,:),1)))
title(['Example, 50% contrast']);

% Average n = x, 50% contrast
subplot(2,3,5)
imagesc(times, es_LFP.freq(f), squeeze(nanmean(frequency_pow_diff(:,4,f,:),1)))
title(['Average (n=7), 50% contrast']);

for n = [1 2 4 5]
    subplot(2,3,n)
    set(gca, 'Clim', [-0.5 0.5]); RedWhiteBlue;
    axis xy;
    line([times(end) times(end)], [bLow bHigh], 'color', 'r', 'linewidth',4);
    line([times(end) times(end)], [nLow nHigh], 'color', 'b', 'linewidth',4);
    line([times(end) times(end)], [bLow2 bHigh2], 'color', 'r', 'linewidth',4);
    line([0 0], ylim, 'linestyle','--','color','k')
    line([1 1], ylim, 'linestyle','--','color','k')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none'); colorbar;
end
for n = 1:7
    [Cn(n), rho_n(n)] = corr(Xn', narrow_band_pow(n,:)');
    [Cb(n), rho_b(n)] = corr(Xb', broad_band_pow(n,:)');
end
subplot(2,3,6)
hold off
bar([1 2],[mean(Cb) mean(Cn)],'k')
hold on;
errorbar([2],[mean(Cn)], [sem(Cn)],'b')
errorbar([1],[mean(Cb)], [sem(Cb)],'r')

set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
