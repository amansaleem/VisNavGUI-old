
%if flag=1, calculates spectrograms, sorted by ascending speed, .. if
%flag=0 -> only coherence analysis as output

function [t,f,SV1,SHP,C,PHI,Speed,SV1_ord,SHP_ord,Speed_ord,T,SCross]=spectrogramFun(flag,dataV1,dataHP,VRdata,dig_pho)

%parameters for the spectrograms
params.Fs=30000;
params.fpass=[0 150];
params.tapers=[5 9];
movingwin=[8  0.8];

%spectra and coherence analysis
disp('Calculating spectra and coherence')
[C,phi,S12,S1,S2,t,f]=cohgramc(dataV1',dataHP',movingwin,params);


%Coherency phase
PHI=angle(phi);


if flag==1
%whitening
disp('Whitening the spectra')
for i=1:size(S1,1)
SV1(i,:)=f.*S1(i,:);
SHP(i,:)=f.*S2(i,:);
end


%Get the speed of the animal in the VR
disp('Interpolating speed values')
speed=interp1(VRdata.TRIAL.time(2:length(VRdata.TRIAL.time))-VRdata.TRIAL.time(2)+dig_pho(1)/30000, VRdata.ballspeed(2:length(VRdata.ballspeed)),t);
Speed=smooth(speed,20);

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

%cross-spectrum
SCross=abs(sqrt(real(S12).^2+imag(S12).^2));

else
   SV1=0;
   SHP=0;
   Speed=0;
   SV1_ord=0;
   SHP_ord=0;
   Speed_ord=0;
   SCross=0;
   T=0;
    
end



end
