function obj = Calculate2Dmaps(obj, varXname, varYname, Tsmthwin, Xbinsize, Ybinsize, Xsmthwin, Ysmthwin, delayT, FcircularX, FcircularY)
size_th = 100;
if nargin < 7
    Xsmthwin = 2*Xbinsize;
    Ysmthwin = 2*Ybinsize;
end
if nargin < 9
    delayT = 0;
end
if nargin < 10 
    FcircularX = true;
    FcircularY = true;
end
disp(['time window = ' num2str(Tsmthwin) ' ms']);

nbcont = numel(obj.SubsetVal.contrast) + 1;
nbgain = numel(obj.SubsetVal.gain);
nbroomlength = numel(obj.SubsetVal.roomlength);
nboutcome = numel(obj.SubsetVal.outcome);

spiketrain = circshift(obj.data.es.spikeTrain,[-delayT 0]);
sampleRate = mean(1./obj.data.es.sampleSize);
for icell = 1:size(spiketrain,2)
    spiketrain(:,icell) = smthInTime(spiketrain(:,icell), sampleRate, Tsmthwin, 'same', [], 'boxcarsum_centered');
end

if FcircularX
    if strcmp(varXname,'trajPercent')
        minvarXname = 0;
        maxvarXname = 100;
    else
        minvarXname = min(obj.data.es.(varXname));
        maxvarXname = max(obj.data.es.(varXname));
    end
    for k = 1:size(obj.data.es.(varXname),2)
        Xtemp = unwrap(obj.data.es.(varXname)(:,k)/maxvarXname*2*pi)*maxvarXname/(2*pi);
        Xtemp = smthInTime(Xtemp, sampleRate, Tsmthwin, 'same', [], 'boxcar_centered');
        varX(:,k) = mod(Xtemp,maxvarXname);
    end
else
    for k = 1:size(obj.data.es.(varXname),2)
        varX(:,k) = smthInTime(obj.data.es.(varXname)(:,k), sampleRate, Tsmthwin, 'same', [], 'boxcar_centered');
    end
end
if FcircularY
    if strcmp(varYname,'LFPphase') || strcmp(varYname,'LFPphase2')
        minvarYname = 0;
        maxvarYname = 360;
    else
        minvarYname = min(obj.data.es.(varYname));
        maxvarYname = max(obj.data.es.(varYname));
    end
    for k = 1:size(obj.data.es.(varYname),2)
        Ytemp = unwrap(obj.data.es.(varYname)(:,k)/maxvarYname*2*pi)*maxvarYname/(2*pi);
        Ytemp = smthInTime(Ytemp, sampleRate, Tsmthwin, 'same', [], 'boxcar_centered');
        varY(:,k) = mod(Ytemp,maxvarYname);
    end
else
    for k = 1:size(obj.data.es.(varYname),2)
        varY(:,k) = obj.data.es.(varYname)(:,k);%smthInTime(obj.data.es.(varYname), sampleRate, Tsmthwin, 'same', [], 'boxcar_centered');
    end
end

numBinsX = floor(max(varX(:,1))/Xbinsize) + 1;
for k = 1:size(obj.data.es.(varXname),2)
    [varX(:,k), binsX] = normalise1var(varX(:,k), numBinsX, [], [minvarXname maxvarXname]);
end
if size(varX,2) == 2
    varX_probe1 = repmat(varX(:,1),[1 sum(obj.CellInfo.Probe == 1)]);
    varX_probe2 = repmat(varX(:,2),[1 sum(obj.CellInfo.Probe == 2)]);
    varX = cat(2,varX_probe1,varX_probe2);
end

if strcmp(varYname,'ballspeed') || strcmp(varYname,'smthBallSpd')
    cbase = find(obj.SubsetVal.contrast == mode(obj.data.es.contrast));
    gbase = find(obj.SubsetVal.gain == mode(obj.data.es.gain));
    rbase = find(obj.SubsetVal.roomlength == mode(obj.data.es.roomLength));
    obase = find(obj.SubsetVal.outcome == 2);
    idxref = obj.getSubsets(cbase, gbase, rbase, obase, obj.SpeedThresh);
    
    %check that the following works
    speedprofile = NaN(numBinsX,1);
    varYbinned = NaN(size(varY));
    for xx = 1:numBinsX
        speedprofile(xx) = median(varY(idxref & varX==xx));
    end
    speed = (varY-speedprofile(varX));%./(speedprofile(X));
    spdrange = 10;
    numBinsY = floor(spdrange/Ybinsize)+1 + 2;
    if numBinsY > 2
        binsize = 10/(numBinsY-2-1);
    else
        binsize = 0;
    end
    minbinval = -5 - binsize/2;
    maxbinval = 5 + binsize/2;
    speed(speed < minbinval) = minbinval - 1;
    speed(speed > maxbinval) = maxbinval + 1;
    [varYbinned(speed >= minbinval & speed <= maxbinval), ~] = normalise1var(speed(speed >= minbinval & speed <= maxbinval), numBinsY-2,[],[minbinval maxbinval]);
    varYbinned(speed >= minbinval & speed <= maxbinval) = varYbinned(speed >= minbinval & speed <= maxbinval) + 1;
    varYbinned(speed < minbinval) = 1;
    varYbinned(speed > maxbinval) = numBinsY;
    varY = varYbinned;
    
%     spdprofile = fast1Dmap(varX(idxref),varY(idxref),1,1,1/(floor(Xsmthwin/Xbinsize)/numBinsX),FcircularX);
%     varY = varY./(spdprofile(varX));
%     varY(varY < 0.5) = 0.4;
%     varY(varY > 1.5) = 1.6;
%     numBinsY = floor((max(varY(idxref))-min(varY(idxref)))/Ybinsize) + 1;
%     [varY, binsY] = normalise1var(varY, numBinsY);

    
%     spdquantilelim = zeros(numBinsX,2);
%     varYBinned = NaN(size(varY));
%     numBinsY = Ybinsize;
%     binsY = zeros(1,numBinsY);    
%     if numBinsY > 1
%         for spd = 1:numBinsY
%             for xx = 1:numBinsX
%                 spdquantilelim(xx,1) = quantile(varY(idxref & varX == xx),max(0,(spd-1)/numBinsY));
%                 spdquantilelim(xx,2) = quantile(varY(idxref & varX == xx),min(1,(spd)/numBinsY));
%             end
%             varYBinned(varY(:) >= spdquantilelim(varX(:),1) & varY(:) < spdquantilelim(varX(:),2)) = spd;
%             if spd == 1
%                 varYBinned(varY(:) <= spdquantilelim(varX(:),1)) = spd;
%             end
%             if spd == numBinsY
%                 varYBinned(varY(:) >= spdquantilelim(varX(:),2)) = spd;
%             end
%             binsY(spd) = mean(0.5*(spdquantilelim(:,1)+spdquantilelim(:,2)));
%         end
%     else
%         varYBinned = varY;
%     end
%     varY = varYBinned;
else
    numBinsY = floor(max(varY(:,1))/Ybinsize) + 1;
end
for k = 1:size(obj.data.es.(varYname),2)
    [varY(:,k), binsY] = normalise1var(varY(:,k), numBinsY, [], [minvarYname maxvarYname]);
end
if size(varY,2) == 2
    varY_probe1 = repmat(varY(:,1),[1 sum(obj.CellInfo.Probe == 1)]);
    varY_probe2 = repmat(varY(:,2),[1 sum(obj.CellInfo.Probe == 2)]);
    varY = cat(2,varY_probe1,varY_probe2);
end

win = max(1,round(Tsmthwin/(1000/sampleRate)));
T = obj.data.es.sampleSize.*win;
try
obj.maps2d.([varXname '_' varYname]) = [];
for c = 1:nbcont
    for g = 1:nbgain
        for r = 1:nbroomlength
            for o = 1:nboutcome
                obj.maps2d.([varXname '_' varYname]){c, g, r, o} = TtwoDimMap('Xsmth_win', floor(Xsmthwin/Xbinsize), 'Ysmth_win', floor(Ysmthwin/Ybinsize));
                obj.maps2d.([varXname '_' varYname]){c, g, r, o}.FcircularX = FcircularX;
                obj.maps2d.([varXname '_' varYname]){c, g, r, o}.FcircularY = FcircularY;
                obj.maps2d.([varXname '_' varYname]){c, g, r, o}.numBinsX = numBinsX;
                obj.maps2d.([varXname '_' varYname]){c, g, r, o}.numBinsY = numBinsY;
                obj.maps2d.([varXname '_' varYname]){c, g, r, o}.binsX = binsX;
                obj.maps2d.([varXname '_' varYname]){c, g, r, o}.binsY = binsY;
                obj.maps2d.([varXname '_' varYname]){c, g, r, o}.qthreshold = 1;
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.Fdiscarditer = true;
                if c > numel(obj.SubsetVal.contrast)
                    contidx = find(obj.SubsetVal.contrast>0);
                else
                    contidx = c;
                end
                idx = obj.getSubsets(contidx, g, r, o, obj.SpeedThresh, true, true);
                if o == 3
                    if sum(idx) > size_th
                        obj.maps2d.([varXname '_' varYname]){c, g, r, o}.trainSpikeMap(varX(idx,:), varY(idx,:), spiketrain(idx,:), T(idx));
                    end
                end
                
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o} = TtwoDimMap('Xsmth_win', floor(Xsmthwin/Xbinsize), 'Ysmth_win', floor(Ysmthwin/Ybinsize));
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.FcircularX = FcircularX;
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.FcircularY = FcircularY;
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.numBinsX = numBinsX;
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.numBinsY = numBinsY;
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.binsX = binsX;
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.binsY = binsY;
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.qthreshold = 1;
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.kfold = 2;
                obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.Fdiscarditer = false;
                if c > numel(obj.SubsetVal.contrast)
                    contidx = find(obj.SubsetVal.contrast>0);
                else
                    contidx = c;
                end
                idx = obj.getSubsets(contidx, g, r, o, obj.SpeedThresh, true, true);
                if o == 3
                    if sum(idx) > size_th
                        obj.maps2d.([varXname '_' varYname '_2fold']){c, g, r, o}.trainSpikeMap(varX(idx,:), varY(idx,:), spiketrain(idx,:), T(idx));
                    end
                end
            end
        end
    end
end
catch
    keyboard
end
end