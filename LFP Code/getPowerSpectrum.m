function [f_exp1, f_exp2, pow_exp1, pow_exp2] = getPowerSpectrum(animal,iseries,iexp1,iexp2,chn)

global pepNEV;
global DIRS;

SetDefaultDirs;

if nargin<5
    chn = 16;
    flag_chn = 0;
else
    flag_chn = 1;
end

CHANNELS_ORDER = MichiganGetLayout(animal,iseries,16);

expName1 = [animal '_' num2str(iseries) '_' num2str(iexp1) '.ns5'];
ns5file1 = ['\\ZSERVER\Data\Cerebus\' animal filesep num2str(iseries) filesep num2str(iexp1) filesep expName1];

[~,SamplingRateInKHZ,nchan] = nsopen2(ns5file1);
for ichan = 1:nchan
    data_exp1(ichan,:) = int16(pepNEV.ns.Data.data(ichan,:));
end
data_exp1 = data_exp1(CHANNELS_ORDER,:);

expName2 = [animal '_' num2str(iseries) '_' num2str(iexp2) '.ns5'];
ns5file2 = ['\\ZSERVER\Data\Cerebus\' animal filesep num2str(iseries) filesep num2str(iexp2) filesep expName2];

[~,SamplingRateInKHZ,nchan] = nsopen2(ns5file2);
for ichan = 1:nchan
    data_exp2(ichan,:) = int16(pepNEV.ns.Data.data(ichan,:));
end
data_exp2 = data_exp2(CHANNELS_ORDER,:);

figure; 
subplot(1,2,1)
if ~flag_chn
    [f_exp1, pow_exp1] = powerSpectrum(data_exp1(chn,:),30000,0,'b');
else
    for chnIdx = 1:length(chn)
        chn(chnIdx)
        [f_exp1(chnIdx,:), pow_exp1(chnIdx,:)] = powerSpectrum(data_exp1(chn(chnIdx),:),30000,0,'b');
    end
end

clear data_cl;

if ~flag_chn
    [f_exp2, pow_exp2] = powerSpectrum(data_exp2(chn,:),30000,0,'r');
else
    for chnIdx = 1:length(chn)
        [f_exp2(chnIdx,:), pow_exp2(chnIdx,:)] = powerSpectrum(data_exp2(chn(chnIdx),:),30000,0,'r');
    end
end

if ~flag_chn
    legend('Closed-loop', 'Open-loop');
    title([expName1 '-' num2str(iexp1) '-' num2str(iexp2)]);
    xlabel('Frequency(Hz)')
    ylabel('Power');
    set(gca,'XLim',[1 101])
    
    subplot(1,2,2)
    for n = 1:size(f_exp1,1)
    loglog(f_exp1(n,:),  fastsmooth(pow_exp1(n,:).*f_exp1(n,:),50,3,1), 'b.', f_exp2(n,:),  fastsmooth(pow_exp2(n,:).*f_exp2(n,:),50,3,1), 'r.', 'linewidth', 1.5);
    legend('Closed-loop', 'Open-loop');
    xlabel('Frequency(Hz)')
    ylabel('Power * freq');
    set(gca,'XLim',[1 101])
    end
end