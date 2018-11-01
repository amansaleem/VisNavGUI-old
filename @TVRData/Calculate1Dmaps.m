function obj = Calculate1Dmaps(obj, varXname, Tsmthwin, Xbinsize, Xsmthwin, delayT, Fcircular)
size_th = 100;
if nargin < 5
    Xsmthwin = 2*Xbinsize;
end
if nargin < 6
    delayT = 0;
end
if nargin < 7
    Fcircular = true;
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

if strcmp(varXname,'trajPercent')
    minvarXname = 0;
    maxvarXname = 100;
elseif strcmp(varXname,'LFPphase') || strcmp(varXname,'LFPphase2')
    minvarXname = 0;
    maxvarXname = 360;
else
    minvarXname = min(obj.data.es.(varXname));
    maxvarXname = max(obj.data.es.(varXname));
end
if Fcircular
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

for k = 1:size(obj.data.es.(varXname),2)
    numBins = floor(max(obj.data.es.(varXname)(:,k))/Xbinsize) + 1;
    [varX(:,k), bins] = normalise1var(varX(:,k), numBins, [], [minvarXname maxvarXname]);
end

if size(varX,2) == 2
    varX_probe1 = repmat(varX(:,1),[1 sum(obj.CellInfo.Probe == 1)]);
    varX_probe2 = repmat(varX(:,2),[1 sum(obj.CellInfo.Probe == 2)]);
    varX = cat(2,varX_probe1,varX_probe2);
end

win = max(1,round(Tsmthwin/(1000/sampleRate)));
T = obj.data.es.sampleSize.*win;

obj.maps1d.(varXname) = [];
for c = 1:nbcont
    for g = 1:nbgain
        for r = 1:nbroomlength
            for o = 1:nboutcome
                obj.maps1d.(varXname){c, g, r, o} = ToneDimMap('Xsmth_win', floor(Xsmthwin/Xbinsize));
                obj.maps1d.(varXname){c, g, r, o}.Fcircular = Fcircular;
                obj.maps1d.(varXname){c, g, r, o}.numBins = numBins;
                obj.maps1d.(varXname){c, g, r, o}.bins = bins;
                obj.maps1d.(varXname){c, g, r, o}.qthreshold = 1;
                obj.maps1d.(varXname){c, g, r, o}.Fdiscarditer = false;
                if c > numel(obj.SubsetVal.contrast)
                    contidx = find(obj.SubsetVal.contrast>0);
                else
                    contidx = c;
                end
                idx = obj.getSubsets(contidx, g, r, o);
                if o == 3
                    if sum(idx) > size_th
                        obj.maps1d.(varXname){c, g, r, o}.trainSpikeMap(varX(idx,:), spiketrain(idx,:), T(idx));
                    end
                end
                
                obj.maps1d.([varXname '_2fold']){c, g, r, o} = ToneDimMap('Xsmth_win', floor(Xsmthwin/Xbinsize));
                obj.maps1d.([varXname '_2fold']){c, g, r, o}.Fcircular = Fcircular;
                obj.maps1d.([varXname '_2fold']){c, g, r, o}.numBins = numBins;
                obj.maps1d.([varXname '_2fold']){c, g, r, o}.bins = bins;
                obj.maps1d.([varXname '_2fold']){c, g, r, o}.qthreshold = 1;
                obj.maps1d.([varXname '_2fold']){c, g, r, o}.kfold = 2;
                obj.maps1d.(varXname){c, g, r, o}.Fdiscarditer = false;
                if c > numel(obj.SubsetVal.contrast)
                    contidx = find(obj.SubsetVal.contrast>0);
                else
                    contidx = c;
                end
                idx = obj.getSubsets(contidx, g, r, o);
                if o == 3
                    if sum(idx) > size_th
                        obj.maps1d.([varXname '_2fold']){c, g, r, o}.trainSpikeMap(varX(idx,:), spiketrain(idx,:), T(idx));
                    end
                end
            end
        end
    end
end
end