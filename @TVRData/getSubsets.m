function tidx = getSubsets(obj,contrast,gain,roomlength,outcome,speed_th,FnoBlanks,FnoAfterBlanks)
if nargin < 2 || isempty(contrast)
    contrast = 1:numel(obj.SubsetVal.contrast);
    gain = 1:numel(obj.SubsetVal.gain);
    roomlength = 1:numel(obj.SubsetVal.roomlength);
    outcome = 1:numel(obj.SubsetVal.outcome);
end
if nargin < 3 || isempty(gain)
    gain = 1:numel(obj.SubsetVal.gain);
    roomlength = 1:numel(obj.SubsetVal.roomlength);
    outcome = 1:numel(obj.SubsetVal.outcome);
end
if nargin < 4 || isempty(roomlength)
    roomlength = 1:numel(obj.SubsetVal.roomlength);
    outcome = 1:numel(obj.SubsetVal.outcome);
end
if nargin < 5 || isempty(outcome)
    outcome = 1:numel(obj.SubsetVal.outcome);
end
if nargin < 6
    speed_th = obj.SpeedThresh;
end
if nargin < 7
    FnoBlanks = true;
end
if nargin < 8
    FnoAfterBlanks = true;%false;%
end

contrast = contrast(:);
gain = gain(:);
roomlength = roomlength(:);
outcome = outcome(:);

if max(contrast) <= size(obj.Subset.Xall,1) && max(gain) <= size(obj.Subset.Xall,2) && max(roomlength) <= size(obj.Subset.Xall,3)
    subset = cell2mat(obj.Subset.Xall(contrast, gain, roomlength, min(outcome,size(obj.Subset.Xall,4))));
    subset = reshape(subset,[size(obj.data.es.trialID,1) size(contrast,1)*size(gain,1)*size(roomlength,1)*size(outcome,1)]);
    tidx = logical(sum(subset,2));
else
    tidx = false(size(obj.data.es.smthBallSpd));
end

tidx = tidx & obj.data.es.smthBallSpd > speed_th & ~isnan(obj.data.es.smthBallSpd);% & ~isnan(obj.data.es.smthTrajSpd);

if FnoBlanks && isfield(obj.data.es,'blanks')
    tidx = tidx & obj.Subset.noBlanks;%~(obj.data.es.blanks) & ~(obj.data.es.afterblanks);
end
if FnoAfterBlanks && isfield(obj.data.es,'afterblanks')
    tidx = tidx & obj.Subset.noAfterBlanks;
end
end