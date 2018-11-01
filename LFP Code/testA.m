global pepNEV;
global DIRS;

addpath('\\zserver\Code\Spikes')
addpath('\\zuser\aman\Work\Code\Behaviour analysis')
addpath('\\zuser\aman\Work\Code\LFP analysis')
addpath('\\zuser\aman\Work\Code\General')

SetDefaultDirs;
animal = 'M130207_BALL'
iseries = 313
CHANNELS_ORDER = MichiganGetLayout(animal,iseries)

iexp1 = 106;
expName1 = [animal '_' num2str(iseries) '_' num2str(iexp1) '.ns5']
ns5file1 = ['\\ZSERVER\Data\Cerebus\' animal filesep num2str(iseries) filesep num2str(iexp1) filesep expName1]
[~,SamplingRateInKHZ,nchan] = nsopen2(ns5file1);

data_run = zeros(1, size(pepNEV.ns.Data.data,2));
if size(pepNEV.ns.Data.data,2)>30000*240
    maxC = max(pepNEV.ns.Data.data(35,1:30000*240));
    minC = min(pepNEV.ns.Data.data(35,1:30000*240));
else
    maxC = max(pepNEV.ns.Data.data(35,:));
    minC = min(pepNEV.ns.Data.data(35,:));
end
data_run((abs(diff(pepNEV.ns.Data.data(35,:)-minC)./(maxC-minC)))>0.5) = 1;
data_run_down = sum(reshape(data_run(1:end-rem(size(pepNEV.ns.Data.data,2),500)),500,[]),1);

data_lfp = zeros(1, size(pepNEV.ns.Data.data,2));
% data_lfp(:) = int16(pepNEV.ns.Data.data(3,:));

k_lfp = decimate(double(pepNEV.ns.Data.data(3,:)), 30);
% k_lfp2 = decimate(double(pepNEV.ns.Data.data(18,:)), 30);
% k_run = decimate(data_run, 30);

[f_lfp,pow_lfp] = powerSpectrum(k_lfp,1000,0,'b');
[f_run,pow_run] = powerSpectrum(data_run_down,60,0,'b');

range_lfp = f_lfp<20;
range_run = f_run<20;

% figure;
% plot(f_lfp(range_lfp),pow_lfp(range_lfp).*f_lfp(range_lfp)./max(pow_lfp(range_lfp).*f_lfp(range_lfp)),...
%     f_run(range_run), pow_run(range_run)./max(pow_run(range_run)));
% xlabel('Hz')
% ylabel('Power')
% title(num2str(iexp1));

figure;
plot(f_lfp(range_lfp),smooth(pow_lfp(range_lfp))./max(pow_lfp(range_lfp)),...
    f_run(range_run), smooth(pow_run(range_run))./max(pow_run(range_run)));
xlabel('Hz')
ylabel('Power')
title(num2str(iexp1));