function plotTrigAnalysis(es, type, frange_lims, plot_log, run_sep)

if nargin<4
    plot_log = 0;
end

if nargin<5
    run_sep = 0;
end

frange = (es.freq>frange_lims(1) & es.freq<frange_lims(2)) & ~(es.freq>48 & es.freq<52);
input = eval(['es.' type]);
if strcmp(type,'powAB')
    input = abs(input);
end
t = es.contrast~=0;

if plot_log
    pstring = '(ones(25,1)*es.freq(frange)).*';
else
    pstring = 'log';
end
if strcmp(type,'cohPhi')
    pstring = '';
end

% 
grays = es.traj;
grays(es.traj==0) = 0;
grays(es.traj~=0) = 1;
diffGrays = ([0 diff(grays)']');

running = zeros(size(es.traj));
running(es.smthBallSpd>2) = 1;
diffRun = ([0 diff(running)']');

VR_starts   = find(diffGrays==1 & t);
VR_ends     = find(diffGrays==-1 & circshift(t,2));
lickTimes   = find(es.lick & t);
rewTimes    = find(~isnan(es.reward) & t);
run_starts  = find(diffRun==1 & (es.lick==0) & isnan(es.reward));
run_ends    = find(diffRun==-1 & (es.lick==0) & isnan(es.reward));

[trig_pow_lick, alltrig_pow_lick]       = calculateTrigPower(input, lickTimes);
[trig_pow_reward, alltrig_pow_reward]   = calculateTrigPower(input, rewTimes);
[trig_pow_end, alltrig_pow_end]         = calculateTrigPower(input, VR_ends);
[trig_pow_start, alltrig_pow_start]     = calculateTrigPower(input, VR_starts);
[trig_pow_rstart, alltrig_pow_rstart]   = calculateTrigPower(input, run_starts);
[trig_pow_rend, alltrig_pow_rend]       = calculateTrigPower(input, run_ends);

[trig_spd_starts, all_trig_spd_starts]  = calculateTrigPower(es.smthBallSpd, VR_starts);
[trig_spd_rew, all_trig_spd_rew]        = calculateTrigPower(es.smthBallSpd, rewTimes);
[trig_spd_lick, all_trig_spd_lick]      = calculateTrigPower(es.smthBallSpd, lickTimes);
[trig_spd_ends, all_trig_spd_ends]      = calculateTrigPower(es.smthBallSpd, VR_ends);
[trig_spd_rstarts, all_trig_spd_rstarts]= calculateTrigPower(es.smthBallSpd, run_starts);
[trig_spd_rends, all_trig_spd_rends]    = calculateTrigPower(es.smthBallSpd, run_ends);

%% Get the separable predictions
if run_sep
temp = all_trig_spd_ends;
temp(isnan(temp)) = 0;
[sRow_ends,sCol_ends,sScl_ends,~,sRes_ends] = MakeSeparable(temp,0);

temp = all_trig_spd_starts;
temp(isnan(temp)) = 0;
[sRow_starts,sCol_starts,sScl_starts,~,sRes_starts] = MakeSeparable(temp,0);

temp = all_trig_spd_rends;
temp(isnan(temp)) = 0;
[sRow_rends,sCol_rends,sScl_rends,~,sRes_rends] = MakeSeparable(temp,0);

temp = all_trig_spd_rstarts;
temp(isnan(temp)) = 0;
[sRow_rstarts,sCol_rstarts,sScl_rstarts,~,sRes_rstarts] = MakeSeparable(temp,0);

temp = all_trig_spd_lick;
temp(isnan(temp)) = 0;
[sRow_lick,sCol_lick,sScl_lick,~,sRes_lick] = MakeSeparable(temp,0);

temp = all_trig_spd_rew;
temp(isnan(temp)) = 0;
[sRow_rew,sCol_rew,sScl_rew,~,sRes_rew] = MakeSeparable(temp,0);

for ifreq = 1:size(alltrig_pow_end,3)
    temp = (squeeze(alltrig_pow_end(:,:,ifreq)));
    [fRow_ends(ifreq,:),fCol_ends(ifreq,:),fScl_ends(ifreq,:),~,residual] = MakeSeparable(temp,0);
    fRes_ends(ifreq,:) = mean(residual);
    corr_ends(ifreq) = corr(fCol_ends(ifreq,:)', sCol_ends);
    
    temp = (squeeze(alltrig_pow_start(:,:,ifreq)));
    [fRow_starts(ifreq,:),fCol_starts(ifreq,:),fScl_starts(ifreq,:),~,residual] = MakeSeparable(temp,0);
    fRes_starts(ifreq,:) = mean(residual);
    corr_starts(ifreq) = corr(fCol_starts(ifreq,:)', sCol_starts);
    
    temp = (squeeze(alltrig_pow_rend(:,:,ifreq)));
    [fRow_rends(ifreq,:),fCol_rends(ifreq,:),fScl_rends(ifreq,:),~,residual] = MakeSeparable(temp,0);
    fRes_rends(ifreq,:) = mean(residual);
    corr_rends(ifreq) = corr(fCol_rends(ifreq,:)', sCol_rends);
    
    temp = (squeeze(alltrig_pow_rstart(:,:,ifreq)));
    [fRow_rstarts(ifreq,:),fCol_rstarts(ifreq,:),fScl_rstarts(ifreq,:),~,residual] = MakeSeparable(temp,0);
    fRes_rstarts(ifreq,:) = mean(residual);
    corr_rstarts(ifreq) = corr(fCol_rstarts(ifreq,:)', sCol_rstarts);
    
    temp = (squeeze(alltrig_pow_lick(:,:,ifreq)));
    [fRow_lick(ifreq,:),fCol_lick(ifreq,:),fScl_lick(ifreq,:),~,residual] = MakeSeparable(temp,0);
    fRes_lick(ifreq,:) = mean(residual);
    corr_lick(ifreq) = corr(fCol_lick(ifreq,:)', sCol_lick);
    
    temp = (squeeze(alltrig_pow_reward(:,:,ifreq)));
    [fRow_rew(ifreq,:),fCol_rew(ifreq,:),fScl_rew(ifreq,:),~,residual] = MakeSeparable(temp,0);
    fRes_rew(ifreq,:) = mean(residual);
    corr_rew(ifreq) = corr(fCol_rew(ifreq,:)', sCol_rew);
end
end
%% Figure with triggered means
figure('Position',[50 100 750 1000]);

subplot(6,3,[1 4]); 
imagesc(-10:14,es.freq(frange),eval(['transpose(' pstring '(trig_pow_start(:,frange)))'])); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('VR start triggered')
ylabel('Frequency (Hz)')
axis tight
subplot(6,3,7)
% plot(-10:14,min(es.freq(frange)) + 0.33*(max(es.freq(frange))-min(es.freq(frange)))*trig_spd_starts./max(trig_spd_starts),...
errorbar(-10:14,trig_spd_starts,nanstd(all_trig_spd_starts),...
    'k','linewidth',1.5);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis tight;
ylims = ylim;
set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
line([0 0],ylim,'linestyle','--', 'color','k')
ylabel('Speed (cm/s)')
hold off;

% (ones(25,1)*es.freq(frange))'.*
subplot(6,3,[2 5]); 
imagesc(-10:14,es.freq(frange),eval(['transpose(' pstring '(trig_pow_end(:,frange)))'])); axis xy
hold off
line([0 0],ylim,'linestyle','--', 'color','k')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('VR end triggered')
axis tight
subplot(6,3,8)
% plot(-10:14,min(es.freq(frange)) + 0.33*(max(es.freq(frange))-min(es.freq(frange)))*trig_spd_ends./max(trig_spd_ends)...
errorbar(-10:14,trig_spd_ends,nanstd(all_trig_spd_ends)...
    , 'k','linewidth',1.5);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis tight;
ylims = ylim;
set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
line([0 0],ylim,'linestyle','--', 'color','k')
% ylabel('Speed (cm/s)');
hold off;

subplot(6,3,[10 13]); 
imagesc(-10:14,es.freq(frange),eval(['transpose(' pstring '(trig_pow_lick(:,frange)))'])); axis xy
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
line([0 0],ylim,'linestyle','--', 'color','k')
title('LICK triggered')
ylabel('Frequency (Hz)')
axis tight
subplot(6,3,16)
% plot(-10:14,min(es.freq(frange)) + 0.33*(max(es.freq(frange))-min(es.freq(frange)))*trig_spd_lick./max(trig_spd_lick)...
errorbar(-10:14,trig_spd_lick,nanstd(all_trig_spd_lick)...
    , 'k','linewidth',1.5);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis tight;
ylims = ylim;
set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
line([0 0],ylim,'linestyle','--', 'color','k')
hold off;
xlabel('Time bin from trigger')

subplot(6,3,[11 14]); 
imagesc(-10:14,es.freq(frange),eval(['transpose(' pstring '(trig_pow_reward(:,frange)))'])); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('REWARD triggered')
axis tight
subplot(6,3,17)
% plot(-10:14,min(es.freq(frange)) + 0.33*(max(es.freq(frange))-min(es.freq(frange)))*trig_spd_rew./max(trig_spd_rew)...
    errorbar(-10:14,trig_spd_rew,nanstd(all_trig_spd_rew)...
    , 'k','linewidth',1.5);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis tight;
ylims = ylim;
set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
line([0 0],ylim,'linestyle','--', 'color','k')
hold off;

% run stuff
subplot(6,3,[3 6]); 
imagesc(-10:14,es.freq(frange),eval(['transpose(' pstring '(trig_pow_rstart(:,frange)))'])); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('Run start triggered')
% ylabel('Frequency (Hz)')
axis tight
subplot(6,3,9)
% plot(-10:14,min(es.freq(frange)) + 0.33*(max(es.freq(frange))-min(es.freq(frange)))*trig_spd_starts./max(trig_spd_starts),...
errorbar(-10:14,trig_spd_rstarts,nanstd(all_trig_spd_rstarts),...
    'k','linewidth',1.5);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis tight;
ylims = ylim;
set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
line([0 0],ylim,'linestyle','--', 'color','k')
% ylabel('Speed (cm/s)')
hold off;

subplot(6,3,[12 15]); 
imagesc(-10:14,es.freq(frange),eval(['transpose(' pstring '(trig_pow_rend(:,frange)))'])); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('Run end triggered')
axis tight
subplot(6,3,18)
% plot(-10:14,min(es.freq(frange)) + 0.33*(max(es.freq(frange))-min(es.freq(frange)))*trig_spd_ends./max(trig_spd_ends)...
errorbar(-10:14,trig_spd_rends,nanstd(all_trig_spd_rends)...
    , 'k','linewidth',1.5);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis tight;
ylims = ylim;
set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
line([0 0],ylim,'linestyle','--', 'color','k')
% ylabel('Speed (cm/s)');
hold off;

if run_sep

% %% Figure with SVD stuff, lick triggered
% figure('Position',[674   100   560   420]);
% 
% subplot(3,4,[1 5])
% imagesc(-10:14,es.freq(frange),log(trig_pow_lick(:,frange))'); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('VR lick triggered')
% ylabel('Frequency (Hz)')
% axis tight
% 
% subplot(3,4,[4 8])
% imagesc(-10:14,es.freq(frange),fRes_lick(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR start triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,9)
% errorbar(-10:14,trig_spd_lick,nanstd(all_trig_spd_lick),...
%     'k','linewidth',1.5);
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis tight;
% ylims = ylim;
% set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
% line([0 0],ylim,'linestyle','--', 'color','k')
% ylabel('Speed (cm/s)')
% hold off;
% 
% subplot(3,4,[2 6])
% imagesc(-10:14,es.freq(frange),-fRow_lick(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR lick triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,10)
% plot(-10:14, -sRow_lick);
% axis tight
% line([0 0],ylim,'linestyle','--', 'color','k')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis off
% 
% subplot(3,4,[3 7])
% plot(abs(corr_lick(frange)),es.freq(frange),'.-k','linewidth',1.5)
% axis tight
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% set(gca,'xlim',[0 1]);
% % axis off
% % subplot(3,3,9)
% 
% %% Figure with SVD stuff, reward triggered
% figure('Position',[ 1322         100         560         420]);;
% 
% subplot(3,4,[1 5])
% imagesc(-10:14,es.freq(frange),log(trig_pow_reward(:,frange))'); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('VR rew triggered')
% ylabel('Frequency (Hz)')
% axis tight
% 
% subplot(3,4,9)
% errorbar(-10:14,trig_spd_rew,nanstd(all_trig_spd_rew),...
%     'k','linewidth',1.5);
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis tight;
% ylims = ylim;
% set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
% line([0 0],ylim,'linestyle','--', 'color','k')
% ylabel('Speed (cm/s)')
% hold off;
% 
% subplot(3,4,[2 6])
% imagesc(-10:14,es.freq(frange),-fRow_rew(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR rew triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,[4 8])
% imagesc(-10:14,es.freq(frange),fRes_rew(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR start triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,10)
% plot(-10:14, -sRow_rew);
% axis tight
% line([0 0],ylim,'linestyle','--', 'color','k')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis off
% 
% subplot(3,4,[3 7])
% plot(abs(corr_rew(frange)),es.freq(frange),'.-k','linewidth',1.5)
% axis tight
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% set(gca,'xlim',[0 1]);
% % axis off
% % subplot(3,3,9)
% 
% %% Figure with SVD stuff, start triggered
% figure('Position',[678   643   560   420]);
% 
% subplot(3,4,[1 5])
% imagesc(-10:14,es.freq(frange),log(trig_pow_start(:,frange))'); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('VR start triggered')
% ylabel('Frequency (Hz)')
% axis tight
% 
% subplot(3,4,9)
% errorbar(-10:14,trig_spd_starts,nanstd(all_trig_spd_starts),...
%     'k','linewidth',1.5);
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis tight;
% ylims = ylim;
% set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
% line([0 0],ylim,'linestyle','--', 'color','k')
% ylabel('Speed (cm/s)')
% hold off;
% 
% subplot(3,4,[2 6])
% imagesc(-10:14,es.freq(frange),-fRow_starts(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR start triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,[4 8])
% imagesc(-10:14,es.freq(frange),fRes_starts(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR start triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,10)
% plot(-10:14, -sRow_starts);
% axis tight
% line([0 0],ylim,'linestyle','--', 'color','k')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis off
% 
% subplot(3,4,[3 7])
% plot(abs(corr_starts(frange)),es.freq(frange),'.-k','linewidth',1.5)
% axis tight
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% set(gca,'xlim',[0 1]);
% % axis off
% % subplot(3,3,9)
% 
% %% Figure with SVD stuff, end triggered
% figure('Position',[1321         643         560         420]);
% 
% subplot(3,4,[1 5])
% imagesc(-10:14,es.freq(frange),log(trig_pow_end(:,frange))'); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('VR end triggered')
% ylabel('Frequency (Hz)')
% axis tight
% 
% subplot(3,4,[4 8])
% imagesc(-10:14,es.freq(frange),fRes_ends(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR start triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,9)
% errorbar(-10:14,trig_spd_ends,nanstd(all_trig_spd_ends),...
%     'k','linewidth',1.5);
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis tight;
% ylims = ylim;
% set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
% line([0 0],ylim,'linestyle','--', 'color','k')
% ylabel('Speed (cm/s)')
% hold off;
% 
% subplot(3,4,[2 6])
% imagesc(-10:14,es.freq(frange),-fRow_ends(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR end triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,10)
% plot(-10:14, -sRow_ends);
% axis tight
% line([0 0],ylim,'linestyle','--', 'color','k')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis off
% 
% subplot(3,4,[3 7])
% plot(abs(corr_ends(frange)),es.freq(frange),'.-k','linewidth',1.5)
% axis tight
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% set(gca,'xlim',[0 1]);
% % axis off
% % subplot(3,3,9)
% 
% %% Figure with SVD stuff, run start triggered
% figure('Position',[100   643   560   420]);
% 
% subplot(3,4,[1 5])
% imagesc(-10:14,es.freq(frange),log(trig_pow_rstart(:,frange))'); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('Run start triggered')
% ylabel('Frequency (Hz)')
% axis tight
% 
% subplot(3,4,9)
% errorbar(-10:14,trig_spd_rstarts,nanstd(all_trig_spd_rstarts),...
%     'k','linewidth',1.5);
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis tight;
% ylims = ylim;
% set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
% line([0 0],ylim,'linestyle','--', 'color','k')
% ylabel('Speed (cm/s)')
% hold off;
% 
% subplot(3,4,[2 6])
% imagesc(-10:14,es.freq(frange),-fRow_rstarts(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR rstart triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,[4 8])
% imagesc(-10:14,es.freq(frange),fRes_rstarts(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR rstart triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,10)
% plot(-10:14, -sRow_rstarts);
% axis tight
% line([0 0],ylim,'linestyle','--', 'color','k')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis off
% 
% subplot(3,4,[3 7])
% plot(abs(corr_rstarts(frange)),es.freq(frange),'.-k','linewidth',1.5)
% axis tight
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% set(gca,'xlim',[0 1]);
% % axis off
% % subplot(3,3,9)
% 
% %% Figure with SVD stuff, end triggered
% figure('Position',[100         100         560         420]);
% 
% subplot(3,4,[1 5])
% imagesc(-10:14,es.freq(frange),log(trig_pow_rend(:,frange))'); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('Run end triggered')
% ylabel('Frequency (Hz)')
% axis tight
% 
% subplot(3,4,[4 8])
% imagesc(-10:14,es.freq(frange),fRes_rends(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR start triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,9)
% errorbar(-10:14,trig_spd_rends,nanstd(all_trig_spd_rends),...
%     'k','linewidth',1.5);
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis tight;
% ylims = ylim;
% set(gca,'YLim',[0 ylims(2)], 'XTick',[]);
% line([0 0],ylim,'linestyle','--', 'color','k')
% ylabel('Speed (cm/s)')
% hold off;
% 
% subplot(3,4,[2 6])
% imagesc(-10:14,es.freq(frange),-fRow_rends(frange,:)); axis xy
% line([0 0],ylim,'linestyle','--', 'color','k')
% hold off
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('VR rend triggered')
% % ylabel('Frequency (Hz)')
% axis tight
% axis off
% 
% subplot(3,4,10)
% plot(-10:14, -sRow_rends);
% axis tight
% line([0 0],ylim,'linestyle','--', 'color','k')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis off
% 
% subplot(3,4,[3 7])
% plot(abs(corr_rends(frange)),es.freq(frange),'.-k','linewidth',1.5)
% axis tight
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% set(gca,'xlim',[0 1]);
% % axis off
% % subplot(3,3,9)

%% Figure with just the BestRows
figure('Position',[ 899   305   717   705]);
% Run Starts
subplot(2,3,3)
imagesc(-10:14,es.freq(frange),-fRow_rstarts(frange,:)); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('Run Start triggered')
axis tight
% axis off
% Run Stops
subplot(2,3,6)
imagesc(-10:14,es.freq(frange),-fRow_rends(frange,:)); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('Run end triggered')
axis tight
% axis off
% VR starts
subplot(2,3,1)
imagesc(-10:14,es.freq(frange),-fRow_starts(frange,:)); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('VR Start triggered')
axis tight
% axis off

% VR ends
subplot(2,3,2)
imagesc(-10:14,es.freq(frange),-fRow_ends(frange,:)); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('VR end triggered')
axis tight
% axis off

% Licks
subplot(2,3,4)
imagesc(-10:14,es.freq(frange),-fRow_lick(frange,:)); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('Lick triggered')
axis tight
% axis off

% Reward
subplot(2,3,5)
imagesc(-10:14,es.freq(frange),-fRow_rew(frange,:)); axis xy
line([0 0],ylim,'linestyle','--', 'color','k')
hold off
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
title('Reward triggered')
axis tight
% axis off
end