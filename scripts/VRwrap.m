function es = VRwrap(es,wrapfactor,Frecenter)
roomlength = round(max(es.traj));
newroomlength = floor(roomlength/wrapfactor);
newroomPercent = floor(round(max(es.trajPercent))/wrapfactor);
if Frecenter
    rewPos = min(es.rewardPos(1:min(end,10)));
    es.rewardPos = mod(es.rewardPos - rewPos,newroomlength);%
else
    rewPos = 0;
end
es.halfID = single(mod(es.traj - rewPos,roomlength) < newroomlength) + 2*single(mod(es.traj - rewPos,roomlength) >= newroomlength);
es.traj = mod(es.traj - rewPos,newroomlength);%mod(obj.data.es.traj-40,100);
es.trajPercent = mod(es.trajPercent - rewPos,newroomPercent);%mod(obj.data.es.trajPercent-40,100);

es.rewardPos = es.rewardPos/newroomlength*100;%
es.traj = es.traj/newroomlength*100;
es.trajPercent = es.trajPercent/newroomPercent*100;
newroomlength = 100;

es.trialID = ones(size(es.traj));
triallim = find(abs(diff(medfilt1(es.traj,10))) > newroomlength/4 & abs(diff(es.traj)) > newroomlength/4);
falsetrial = find(diff(triallim) < 60) + 1;
triallim(falsetrial) = [];
for tt = 1:numel(triallim)
    es.trialID(triallim(tt):end) = es.trialID(triallim(tt):end) + 1;
end

es.trialgainchange = zeros(size(es.trialID));
idx = [1; find(diff(es.gain) ~= 0); numel(es.trialID)];
for i = 1:numel(idx)-1
    es.trialgainchange(idx(i):idx(i+1)) = es.trialID(idx(i):idx(i+1)) - es.trialID(idx(i)) + 1;
end
end