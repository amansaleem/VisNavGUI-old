function [comod_V1, pow_V1, freq] = plot_comodulation(animal, iseries, iexp, channels)

if nargin<4
    channels = [20 40];
end

es = VRLoadMultipleExpts(animal,iseries,iexp,'BEHAV_LFP',channels);

% frange = es.freq>54 & es.freq<80;
frange = es.freq>20 & es.freq<90 & ~(es.freq>47 & es.freq<53);

figure('Position',[14 549 4000 400]);
subplot(131)
imagesc(es.freq(frange), es.freq(frange), corr(es.powA(:,frange),es.powA(:,frange)))
axis xy, axis equal; axis tight;
set(gca,'color','none','box','off','TickDir','out','fontsize',14)
xlabel('Freq, CA1')
ylabel('Freq, CA1')
title('CA1 - CA1')

subplot(132)
imagesc(es.freq(frange), es.freq(frange), corr(es.powB(:,frange),es.powB(:,frange)))
axis xy, axis equal; axis tight;
set(gca,'color','none','box','off','TickDir','out','fontsize',14)
xlabel('Freq, V1')
ylabel('Freq, V1')
title('V1 - V1')

freq = es.freq(frange);
pow_V1 = es.powB(:,frange);
comod_V1 = corr(es.powB(:,frange),es.powB(:,frange));

subplot(133)
imagesc(es.freq(frange), es.freq(frange), corr(es.powA(:,frange),es.powB(:,frange)))
axis xy, axis equal; axis tight;
set(gca,'color','none','box','off','TickDir','out','fontsize',14)
xlabel('Freq, CA1')
ylabel('Freq, V1')
title('CA1 - V1')