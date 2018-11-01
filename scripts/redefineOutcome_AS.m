function [es] = redefineOutcome_AS(es)
% Give in the input es and this is redefined as the following
% 0: Passive
% 1: Early
% 2: Correct
% 3: Late
% 4: Miss
% 5: TimeOut

%% Finding the trials as correct, early, miss and late
trialEnds = [find(diff(es.trialID)>=1)-15];
temp = zeros(size(es.traj));
temp(trialEnds) = 1;
trialEnds = temp;

outcome.complete        = (es.trajPercent(trialEnds>0)>80) | es.outcome(trialEnds>0)~=1;
outcome.correctTrials   = es.trialID(trialEnds & es.outcome==2);
outcome.earlyTrials     = es.trialID(trialEnds & (es.outcome==0 & es.trajPercent<80));
outcome.misslateTrials  = es.trialID(trialEnds  & es.outcome==0 & es.trajPercent>80);

es.outcome(es.outcome==1) = 0;

% Splitting miss and late trials
outcome.missTrials = [];
outcome.lateTrials = [];
for tIdx = outcome.misslateTrials'
    % Does the animal lick after the reward position?
    if sum(es.lick((es.trialID==tIdx) & (es.trajPercent>70)))>0
        outcome.lateTrials = [outcome.lateTrials tIdx];
    else
        outcome.missTrials = [outcome.missTrials tIdx];
    end
end

% Early
for tIdx = [outcome.earlyTrials']
    es.outcome(es.trialID==tIdx) = 1;
end

% Correct
for tIdx = [outcome.correctTrials']
    es.outcome(es.trialID==tIdx) = 2;
end
% Late
for tIdx = [outcome.lateTrials]
    es.outcome(es.trialID==tIdx) = 3;
end
% Miss
for tIdx = [outcome.missTrials]
    es.outcome(es.trialID==tIdx) = 4;
end

% %% Finding the trials as correct, early, miss and late
% trialEnds = [find(diff(es.trialID)>=1)-15];
% temp = zeros(size(es.traj));
% temp(trialEnds) = 1;
% trialEnds = temp;
% 
% t_low   = es.contrast>0 & es.contrast<0.6;
% t       = es.contrast==0.6 & es.gain==1 & es.roomLength==1 ;
% t_high  = es.contrast>0.6;
% 
% outcome.complete        = (es.trajPercent(trialEnds>0)>80) | es.outcome(trialEnds>0)~=1;
% outcome.correctTrials   = es.trialID(trialEnds & es.outcome==2);
% outcome.earlyTrials     = es.trialID(trialEnds & (es.outcome==0 & es.trajPercent<80));
% outcome.misslateTrials  = es.trialID(trialEnds & es.outcome==0 & es.trajPercent>80);
% outcome.lowContrast     = es.trialID(trialEnds & t_low);
% outcome.normContrast    = es.trialID(trialEnds & t);
% outcome.highContrast    = es.trialID(trialEnds & t_high);
% 
% es.outcome(es.outcome==1) = NaN;
% 
% % Splitting miss and late trials
% outcome.missTrials = [];
% outcome.lateTrials = [];
% for tIdx = outcome.misslateTrials'
%     % Does the animal lick after the reward position?
%     if sum(es.lick((es.trialID==tIdx) & (es.trajPercent>70)))>0
%         outcome.lateTrials = [outcome.lateTrials tIdx];
%     else
%         outcome.missTrials = [outcome.missTrials tIdx];
%     end
% end
% 
% % Early
% for tIdx = [outcome.earlyTrials']
%     es.outcome(es.trialID==tIdx) = 1;
% end
% 
% % Correct
% for tIdx = [outcome.correctTrials']
%     es.outcome(es.trialID==tIdx) = 2;
% end
% % Late
% for tIdx = [outcome.lateTrials]
%     es.outcome(es.trialID==tIdx) = 3;
% end
% % Miss
% for tIdx = [outcome.missTrials]
%     es.outcome(es.trialID==tIdx) = 4;
% end
% 
% if nargout>1
%     outcome.numbers.low.early = length(intersect(outcome.lowContrast, outcome.earlyTrials));
%     outcome.numbers.low.correct = length(intersect(outcome.lowContrast, outcome.correctTrials));
%     outcome.numbers.low.late = length(intersect(outcome.lowContrast, outcome.lateTrials));
%     outcome.numbers.low.miss = length(intersect(outcome.lowContrast, outcome.missTrials));
%     
%     outcome.numbers.norm.early = length(intersect(outcome.normContrast, outcome.earlyTrials));
%     outcome.numbers.norm.correct = length(intersect(outcome.normContrast, outcome.correctTrials));
%     outcome.numbers.norm.late = length(intersect(outcome.normContrast, outcome.lateTrials));
%     outcome.numbers.norm.miss = length(intersect(outcome.normContrast, outcome.missTrials));
%     
%     outcome.numbers.high.early = length(intersect(outcome.highContrast, outcome.earlyTrials));
%     outcome.numbers.high.correct = length(intersect(outcome.highContrast, outcome.correctTrials));
%     outcome.numbers.high.late = length(intersect(outcome.highContrast, outcome.lateTrials));
%     outcome.numbers.high.miss = length(intersect(outcome.highContrast, outcome.missTrials));
% end