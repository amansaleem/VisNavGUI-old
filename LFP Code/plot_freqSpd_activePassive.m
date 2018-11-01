pow_es_CL_active = getPowSpd(es, es.traj~=0 & es.contrast~=0 & es.outcome==2);
freqSpd_CL_active = freqVsSpd(pow_es_CL_active);

pow_es_CL_pass = getPowSpd(es, es.traj~=0 & es.contrast~=0 & (es.outcome==1));
freqSpd_CL_pass = freqVsSpd(pow_es_CL_pass);

pow_es_CL_miss = getPowSpd(es, es.traj~=0 & es.contrast~=0 & (es.outcome==0));
freqSpd_CL_miss = freqVsSpd(pow_es_CL_miss);

pow_es_CL_gray = getPowSpd(es, es.traj==0 | es.contrast==0);
freqSpd_CL_gray = freqVsSpd(pow_es_CL_gray);

figure;
plot(freqSpd_CL_pass.spd(2:end), freqSpd_CL_pass.freqA(2:end),'s-r')
hold on;
plot(freqSpd_CL_active.spd(2:end), freqSpd_CL_active.freqA(2:end), '*-g')
plot(freqSpd_CL_miss.spd(2:end), freqSpd_CL_miss.freqA(2:end), 'o-k')
plot(freqSpd_CL_gray.spd(2:end), freqSpd_CL_gray.freqA(2:end), 'o-', 'color',[.5 .5 .5])
set(gca, 'box','off','TickDir','out','fontsize',14)
legend('Passive','Active','Missed','gray','location','best')
xlabel('Speed (cm/s)')
ylabel('Frequency (Hz)')
% 
% figure
% plot(freqSpd_CL_pass.spd(2:end), freqSpd_CL_pass.freqA(2:end)-min(freqSpd_CL_pass.freqA(2:end)),'o-b')
% hold on;
% plot(freqSpd_CL_active.spd(2:end), freqSpd_CL_active.freqA(2:end)-min(freqSpd_CL_active.freqA(2:end)), 'o-r')
% plot(freqSpd_CL_gray.spd(2:end), freqSpd_CL_gray.freqA(2:end)-miss(freqSpd_CL_gray.freqA(2:end)), 'o-k')
% plot(freqSpd_CL_miss.spd(2:end), freqSpd_CL_miss.freqA(2:end)-min(freqSpd_CL_miss.freqA(2:end)), 'o-', 'color',[.5 .5 .5])
% set(gca, 'box','off','TickDir','out','fontsize',14)
% legend('Passive','Active','Missed','location','best')
% xlabel('Speed (cm/s)')
% ylabel('\Delta Frequency (Hz)')
