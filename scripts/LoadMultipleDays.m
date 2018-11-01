function es = LoadMultipleDays(animal,Fall)
if nargin < 2
    Fall = false;
end

infoAll = getDataInfo(animal);
date_list = [];
for n = 1:length(infoAll)
    date_list = [date_list, {infoAll(n).date}];
end

if ismac
    date_list{1} = 'bloddy mac';
end
if ~Fall
    [Selection,ok] = listdlg('PromptString','Select a series:',...
        'SelectionMode','multiple',...
        'ListString',date_list);
    if ok
        date = date_list(Selection);
        iseries_list = zeros(1,numel(date));
        for n = 1:numel(date)
            iseries_list(n) = str2num(date{n});
        end    
    end
else
    iseries_list = zeros(1,numel(date_list));
    for n = 1:numel(date_list)
        iseries_list(n) = str2num(date_list{n});
    end
end

wrapfactor = 2;
tolerance = 0.2;

[VRdata, ~, es] = VRWheelLoad(animal, iseries_list(1), 101);

roomlength = (floor(max(es.traj)/50)+1) * 50 ;
newroomlength = floor(roomlength/2);

% rewPos = es.rewardPos(1);
% es.rewardPos = mod(es.rewardPos - rewPos,newroomlength);
% es.traj = mod(es.traj - rewPos,newroomlength);
% es.trajPercent = mod(es.trajPercent - rewPos,newroomlength);

es.dayID = iseries_list(1) * ones(numel(es.trialID),1);
es.dayNum = 1 * ones(numel(es.trialID),1);
if isfield(VRdata.EXP, 'nTrialChange')
   es.nTrialChange = VRdata.EXP.nTrialChange * ones(numel(es.trialID),1);
end
VRdata.TRIAL.goodlick(VRdata.TRIAL.goodlick(:,1) == 0,2) = 1000;
es.badtrials = sum(sum((VRdata.TRIAL.goodlick(1:max(es.trialID),:) == 0))) * ones(numel(es.trialID),1);
es.maxbadlicks = VRdata.EXP.maxBadLicks * ones(numel(es.trialID),1);
for iseries = 1:length(iseries_list)
    if iseries == 1
        expstart = 102;
    else
        expstart = 101;
    end
    for iexp = expstart:115
        try
        [VRdata, ~, esX] = VRWheelLoad(animal, iseries_list(iseries), iexp);
                
%         rewPos = esX.rewardPos(1);
%         esX.rewardPos = mod(esX.rewardPos - rewPos,newroomlength);%
%         esX.traj = mod(esX.traj - rewPos,newroomlength);
%         esX.trajPercent = mod(esX.trajPercent - rewPos,newroomlength);

        if isfield(es, 'nTrialChange')
           es.nTrialChange = [es.nTrialChange ; VRdata.EXP.nTrialChange * ones(numel(esX.trialID),1)];
        end
        VRdata.TRIAL.goodlick(VRdata.TRIAL.goodlick(:,1) == 0,2) = 1000;
        es.badtrials = 0;%[es.badtrials ; sum(sum((esX.goodlick(1:max(esX.trialID),:) == 0))) * ones(numel(esX.trialID),1)];
        es.maxbadlicks = [es.maxbadlicks ; VRdata.EXP.maxBadLicks * ones(numel(esX.trialID),1)];
        
        es = combineTwoVRexpts(es, esX);
        
        es.dayID = [es.dayID ; iseries_list(iseries) * ones(numel(esX.trialID),1)];
        es.dayNum = [es.dayNum ; iseries * ones(numel(esX.trialID),1)];
        
        catch
            warning(['no exp ' num2str(iexp) ' for session ' num2str(iseries_list(iseries))]);
        end
    end
end


es.RecDay = false(size(es.dayID));

end