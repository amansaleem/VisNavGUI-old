%getPowerSpectrum(animal,iseries,iexp)
%close all; 
%clear all;

%clc


animal='M120307_BALL';
iseries=322;
iexp=107;

addpath \\zserver\Code\Spikes
global pepNEV;
global DIRS;

SetDefaultDirs;


rsamp=10;
CHANNELS_ORDER = MichiganGetLayout(animal,iseries,33);

expName1 = [animal '_' num2str(iseries) '_' num2str(iexp) '.ns5'];
ns5file1 = ['\\ZSERVER\Data\Cerebus\' animal filesep num2str(iseries) filesep num2str(iexp) filesep expName1];
expName2 = [animal '_' num2str(iseries) '_session_' num2str(iexp) '_trial001.mat'];
VRData = ['\\zserver\Data\ball\' animal filesep num2str(iseries) filesep expName2];

load(VRData)

[~,SamplingRateInKHZ,nchan] = nsopen2(ns5file1);


for ichan = 1:nchan
    data_exp(ichan,:) = resample(double(pepNEV.ns.Data.data(ichan,:)),1,rsamp);
end
data_exp = data_exp(CHANNELS_ORDER,:);


for chn=1:16
[f_exp(chn,:), pow_exp(chn,:)] = powerSpectrum(data_exp(chn,:),rsamp/30000,0,'b');
end

for chn=17:32
[f_exp2(chn-16,:), pow_exp2(chn-16,:)] = powerSpectrum(data_exp(chn,:),rsamp/30000,0,'r');
end


figure;
loglog(f_exp2(13,:), f_exp2(13,:)'.*smooth(pow_exp2(13,:),1),'r')
hold on
loglog(f_exp(12,:),f_exp(12,:)'.*smooth(pow_exp(12,:),1),'b')
axis tight
title(['powerspectrum in hippocampus and V1 logarithmic scale, session ' num2str(iexp) ',channel 12 and 29'])

figure;
for i=1:16
loglog(f_exp2(i,:), f_exp2(i,:)'.*smooth(pow_exp2(i,:),1),'r')
hold on
end
title('hippocampus all channels')
axis tight

figure;  
for i=1:16
loglog(f_exp(i,:),f_exp(10,:)'.*smooth(pow_exp(i,:),1),'b')
hold on
end
title('V1 all channels')
axis tight

figure;
for i=4:16
loglog(f_exp2(i,:), f_exp2(i,:)'.*smooth(pow_exp2(i,:),20),'r')
hold on
loglog(f_exp(i,:),f_exp(10,:)'.*smooth(pow_exp(i,:),20),'b')
end
title('hippocampus red, V1 blue, channels from 4 to 16')
axis tight



%%%%% to have same frequency resolution in power spectrum f_exp,pow_exp


% HighLim=find(f_exp<300, 1, 'last' );
% Freq=logspace(log10(min(f_exp)),log10(300),HighLim);
% pow_expF=smooth(f_exp.*pow_exp,10);
% Pow=interp1(f_exp(1:HighLim),pow_expF(1:HighLim),Freq);
% 
% figure;
% loglog(f_exp(1:HighLim),pow_expF(1:HighLim),'k')
% hold on
% loglog(Freq,Pow,'r')
% axis tight















