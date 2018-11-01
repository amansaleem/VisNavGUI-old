function [es, es_behav] = plot_pow_spd_rew(animal, iseries, exp_list, freq_range, chn)

es = VRLoadMultipleExpts(animal,   iseries,   exp_list,'BEHAV_LFP',[4 chn]);

es_behav = VRLoadMultipleExpts(animal,   iseries,   exp_list);

frange = es.freq>freq_range(1) & es.freq<freq_range(2);

offset = (freq_range(1)+3);
gain   = (freq_range(2)-freq_range(1) - 3)*0.8;

lick = zeros(size(es_behav.traj));
lick(es_behav.lick>0) = 3;

rew = zeros(size(es_behav.reward));
rew(es_behav.reward>1) = 1;

figure;
subplot(211)
imagesc(es.sampleTimes, es.freq(frange), log(es.powB(:,frange)')); axis xy
hold on;
plot(es_behav.sampleTimes, offset + gain*(es_behav.smthBallSpd./max(es_behav.smthBallSpd)),'k', 'linewidth',2);
set(gca,'color','none','box','off','TickDir','out','fontsize',14)
ylabel('Freq (Hz)');
title([animal ' ' num2str(iseries) ' ' num2str(exp_list) '  Run speed and spectrogram']);

subplot(212)
imagesc(es.sampleTimes, es.freq(frange), log(es.powB(:,frange)')); axis xy
hold on;
plot(es_behav.sampleTimes, offset + gain*rew,'k', 'linewidth',3);
plot(es_behav.sampleTimes, offset + lick,'b', 'linewidth',1.5);
plot(es_behav.sampleTimes, offset - lick,'b', 'linewidth',1.5);
set(gca,'color','none','box','off','TickDir','out','fontsize',14)
title('Reward, Licks and Spectrogram');
ylabel('Freq (Hz)');
xlabel('Time (s)');
