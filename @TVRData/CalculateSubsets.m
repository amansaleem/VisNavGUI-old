function obj = CalculateSubsets(obj)
es = obj.data.es;

speed_th = obj.SpeedThresh;

obj.SubsetVal.contrast = sort(unique(es.contrast(~isnan(es.contrast))),'ascend');
obj.SubsetVal.gain = sort(unique(es.gain(~isnan(es.gain))),'ascend');
obj.SubsetVal.roomlength = sort(unique(es.roomLength(~isnan(es.roomLength))),'ascend');
obj.SubsetVal.outcome = sort(unique(es.outcome(~isnan(es.outcome))),'ascend');

obj.SubsetVal.gain = obj.SubsetVal.gain(:)';
obj.SubsetVal.contrast = obj.SubsetVal.contrast(:)';
obj.SubsetVal.roomlength = obj.SubsetVal.roomlength(:)';
obj.SubsetVal.outcome = obj.SubsetVal.outcome(:)';

nbcont = numel(obj.SubsetVal.contrast);
nbgain = numel(obj.SubsetVal.gain);
nbroomlength = numel(obj.SubsetVal.roomlength);
nboutcome = numel(obj.SubsetVal.outcome);

obj.Subset.Spd = es.smthBallSpd > speed_th & ~isnan(es.smthBallSpd) & ~isnan(es.smthTrajSpd);

obj.Subset.Spdprofile = true(size(obj.Subset.Spd));

qlimit = 0.01;
itraj = floor(es.trajPercent)+1;
ntrajbins = max(itraj);
spdquantilelim = zeros(ntrajbins,2);
obj.Bayes.Speeds = NaN(size(obj.data.es.trajspeed));
idxref = es.outcome==2;
for xx = 1:ntrajbins
    spdquantilelim(xx,1) = quantile(es.smthBallSpd(idxref & itraj == xx & es.smthBallSpd > speed_th),qlimit);
    spdquantilelim(xx,2) = quantile(es.smthBallSpd(idxref & itraj == xx & es.smthBallSpd > speed_th),1-qlimit);
end
if es.CircularMaze
    spdquanttemp = repmat(spdquantilelim,[3 1]);
    spdquanttemp(:,1) = smooth(spdquanttemp(:,1),5);
    spdquanttemp(:,2) = smooth(spdquanttemp(:,2),5);
    spdquantilelim(:,1) = spdquanttemp(ntrajbins+1:2*ntrajbins,1);
    spdquantilelim(:,2) = spdquanttemp(ntrajbins+1:2*ntrajbins,2);
else
    spdquanttemp = spdquantilelim;
    spdquantilelim(:,1) = smooth(spdquanttemp(:,1),5);
    spdquantilelim(:,2) = smooth(spdquanttemp(:,2),5);
end
spdquanttemp = [];

obj.Subset.Spdprofile(es.smthBallSpd < spdquantilelim(itraj,1) | es.smthBallSpd > spdquantilelim(itraj,2)) = false;

if islogical(es.blanks(1))
    obj.Subset.noBlanks = ~es.blanks;
    obj.Subset.noAfterBlanks = ~es.afterblanks;
else
    obj.Subset.noBlanks = true(size(es.blanks));%es.smthBallSpd > speed_th & ~isnan(es.smthBallSpd) & ~isnan(es.smthTrajSpd);
    obj.Subset.noAfterBlanks = true(size(es.blanks));
end
obj.Subset.contrast = cell(nbcont,1);
obj.Subset.gain = cell(nbgain,1);
obj.Subset.roomlength = cell(nbroomlength,1);
obj.Subset.outcome = cell(nboutcome,1);
obj.Subset.Xall = cell(nbcont, nbgain, nbroomlength, nboutcome);

for i = 1:nbcont
    obj.Subset.contrast{i} = es.contrast==obj.SubsetVal.contrast(i);
end
for i = 1:nbgain
    obj.Subset.gain{i} = es.gain==obj.SubsetVal.gain(i);
end
for i = 1:nbroomlength
    obj.Subset.roomlength{i} = es.roomLength==obj.SubsetVal.roomlength(i);
end
for i = 1:nboutcome
    obj.Subset.outcome{i} = es.outcome==obj.SubsetVal.outcome(i);
end

for c = 1:nbcont
    for g = 1:nbgain
        for r = 1:nbroomlength
            for o = 1:nboutcome
                obj.Subset.Xall{c, g, r, o} = es.contrast==obj.SubsetVal.contrast(c) & ...
                                              es.gain==obj.SubsetVal.gain(g) & ...
                                              es.roomLength==obj.SubsetVal.roomlength(r) & ...
                                              es.outcome==obj.SubsetVal.outcome(o);
            end
        end
    end
end
end