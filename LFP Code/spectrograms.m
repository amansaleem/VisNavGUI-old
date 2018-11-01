
close all 
clear all
clc

%info on which data to analyze
animal='M120307_BALL';
iseries=322;
iexp=107;

V1chn=10;
HPchn=27;

%load the data
addpath \\zserver\Code\Spikes
global pepNEV;
global DIRS;
SetDefaultDirs;

disp('Loading VR data')
VRdata=['VRdatasession' num2str(iexp)]; 
load(VRdata)


disp('Loading LFP data')
CHANNELS_ORDER = MichiganGetLayout(animal,iseries,33);
expName1 = [animal '_' num2str(iseries) '_' num2str(iexp) '.ns5'];
ns5file1 = ['\\ZSERVER\Data\Cerebus\' animal filesep num2str(iseries) filesep num2str(iexp) filesep expName1];
[~,SamplingRateInKHZ,nchan] = nsopen2(ns5file1);
V1Chn=find(CHANNELS_ORDER==V1chn);
HPChn=find(CHANNELS_ORDER==HPchn);
dataHP= double(pepNEV.ns.Data.data(HPChn,:)); 
dataV1= double(pepNEV.ns.Data.data(V1Chn,:));

disp('Loading photodiode data')
NameTimes = [animal '_' num2str(iseries) '_' num2str(iexp) '_screenTimes.mat'];
filenameScreentime=['\\ZSERVER\Data\multichanspikes\' animal filesep num2str(iseries) filesep NameTimes];
load(filenameScreentime)

%parameters for the spectrograms
params.Fs=30000;
params.fpass=[0 300];
params.tapers=[5 9];
movingwin=[8  0.8];

%spectra and coherence analysis
disp('Calculating spectra and coherence')
[C,phi,S12,S1,S2,t,f]=cohgramc(dataV1',dataHP',movingwin,params);

%whitening
disp('Whitening the spectra')
for i=1:size(S1,1)
SV1(i,:)=f.*S1(i,:);
SHP(i,:)=f.*S2(i,:);
end


%Get the speed of the animal in the VR
disp('Interpolating speed values')
speed=interp1(VRdata.TRIAL.time(2:length(VRdata.TRIAL.time))-VRdata.TRIAL.time(2)+dig_pho(1)/30000, VRdata.ballspeed(2:length(VRdata.ballspeed)),t);
Speed=smooth(speed,30);

%sort the timepoints
disp('Ordering spectra for ascending speed')
Speed_ord=sort(Speed);
nf=length(f);
count=1;
for S=1:length(t)
    if S==1 || Speed_ord(S)~=Speed_ord(S-1)
    times=find(Speed==Speed_ord(S));
    SV1_ord(count:count+length(times)-1,1:nf)=SV1(times,1:nf);
    SHP_ord(count:count+length(times)-1,1:nf)=SHP(times,1:nf);
    T(count:count+length(times)-1)=times;
    count=count+length(times);
    clear times
    end
end


% Plotting
%spectrogram of V1
figure;
plot_matrix(SV1,t,f)
title(['Spectrogram of V1, channel ' num2str(V1chn) ', session ' num2str(iexp)])
axis tight

%spectrogram of HP
figure;
plot_matrix(SHP,t,f) 
title(['Spectrogram of hippocampus, channel ' num2str(HPchn) ', session ' num2str(iexp)])
axis tight

%Coherency magnitude
figure;
plot_matrix(C,t,f,'n')  
title(['Coherency magnitude, session' num2str(iexp)])
axis tight

%Coherency phase
PHI=angle(phi);
figure;
plot_matrix(PHI,t,f)
title(['Coherency phase, session' num2str(iexp)])
axis tight

%spectrogram and animal speed
figure;
subplot(3,1,1)
plot_matrix(SV1,t,f)
title(['Spectrogram of V1, channel ' num2str(V1chn) ', session ' num2str(iexp)])
axis tight

subplot(3,1,2)
plot_matrix(SHP,t,f) 
title(['Spectrogram of hippocampus, channel ' num2str(HPchn)])
axis tight

subplot(3,1,3)
plot(t,Speed)
title('Speed of the animal')
axis tight

%Spectrogram with increasing speed
figure;
subplot(3,1,1)
plot_matrix(SV1_ord,t,f)
title(['Spectrogram of V1, channel ' num2str(V1chn) ', session ' num2str(iexp)])
axis tight 
set(gca,'XTick',t(1):50:t(end))
set(gca,'XTickLabel',T(1:50:end))

subplot(3,1,2)
plot_matrix(SHP_ord,t,f) 
title(['Spectrogram of hippocampus, channel ' num2str(HPchn)])
axis tight
set(gca,'XTick',t(1):50:t(end))
set(gca,'XTickLabel',T(1:50:end))

subplot(3,1,3)
plot(t,Speed_ord)
title('Speed of the animal')
axis tight
set(gca,'XTick',t(1):50:t(end))
set(gca,'XTickLabel',T(1:50:end))


%cross-spectrum
SCross=abs(sqrt(real(S12).^2+imag(S12).^2));
figure;
plot_matrix(SCross,t,f) 
title('Cross-spectrum of hippocampusand V1')
axis tight



%save graphs and data

cd \\zserver\Lab\Tmp\Chiara\channels10and27_session107\Spectrograms_finerFreqRes
 
workspaceName=['SpectrogramSession' num2str(iexp) '_channels' num2str(V1chn) '_' num2str(HPchn) '.mat'];
save workspaceName
saveas(1,'V1')
saveas(2,'HP')
saveas(3,'coherencyMag')
saveas(4,'coherencyPhase')
saveas(5,'V1HPSpeed')
saveas(6,'V1HPSpeed_Ascending')
saveas(7,'Cross-spectrum')

disp('SAVED!!')



