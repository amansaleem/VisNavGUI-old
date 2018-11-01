function  screenTimes = VRcheckscreenTimes(VRdata, screenTimes)
screenTimes2 = sort(VRdata.TRIAL.time(:));
maxTime = min(find(isnan(screenTimes2)))-1;
screenTimes2 = screenTimes2(1:maxTime);


screenTimes = screenTimes;
display('WARNING!!!! temporary fix. Check the photodiode signalling');
screenTimes = screenTimes'-screenTimes(1);
%condition added to deal with cases where the recording has more
%screentimes than VR at the end. Please check that additional
%screentimes are indeed at the end of the recording and not anywhere
%else.
if numel(screenTimes) < numel(screenTimes2)
    maxTime  = min(find((screenTimes2-max(screenTimes))>0));
    screenTimes = screenTimes(1:maxTime);
else
    display('WARNING!!!! TOO MANY SCREENTIMES IN THE RECORDING: please check that they are at the end')
    maxTime  = min(find((screenTimes-max(screenTimes2))>0));
    screenTimes = screenTimes(1:maxTime);
end
end