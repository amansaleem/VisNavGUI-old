function [outputMatrixTime,outputMatrix,outputMatrixSpont, p, cellIDs, outputTimes] = getStimSpiketimes(animal, iseries, iexp, igroup, iaddinfo)

if nargin<4
    igroup = 8;
end
if nargout>3
   getTimes = true;
else
   getTimes = false;
end
if nargin<5
    iaddinfo = '1';
end

if numel(iaddinfo) == 1 && numel(igroup) > 1
    iaddinfotemp = iaddinfo{1};
    iaddinfo = cell(1,numel(igroup));
    for group_idx = 1:numel(igroup)
        iaddinfo{group_idx} = iaddinfotemp;
    end
end

% Load spike times and stim times
for igroup_idx = 1:length(igroup)
    try
    temp = getKwikSpikes(animal, iseries, iexp, igroup(igroup_idx),iaddinfo{igroup_idx});
    tempNumCells = length(temp);
    if igroup_idx == 1
        chans=temp;
        numCells = tempNumCells;
    else
        for icell = 1:tempNumCells
            chans(numCells+1) = temp(icell);
            numCells = numCells + 1;
        end
    end
    catch
    end
end

[stimTimes, p] = createStimMatrix(animal, iseries, iexp);

stimTimes.on = stimTimes.on./30000;
stimTimes.off = stimTimes.off./30000;

numStim = size(stimTimes.on,1);
numReps = size(stimTimes.on,2);
cellIDs = zeros(1,length(chans));

% stimTimes.on is istim x irep
% figure
wingauss = gausswin(5);
wingauss = wingauss/sum(wingauss);
for icell = 1:length(chans)
    st = chans(icell).spiketimes;
    cellIDs(icell) = chans(icell).icell;
    outputMatrix{icell} = zeros(numStim, numReps);
    outputMatrixTime{icell} = zeros(numStim, numReps,61);
    outputMatrixSpont{icell} = zeros(numStim, numReps);
    for istim = 1:numStim
        for irep = 1:numReps
            startTime = stimTimes.on(istim,irep);
            endTime   = stimTimes.off(istim,irep);
            if getTimes
                outputTimes(icell, istim, irep).times = st(st>(startTime-0.5) & st<=(endTime+0.5)) - startTime;
            end
            outputMatrix{icell}(istim, irep) = sum(st>startTime & st<=endTime);
            for tt = -20:40
                outputMatrixTime{icell}(istim, irep,tt+21) = sum(st>=startTime + tt*0.05 & st<startTime + (tt+1)*0.05);
            end
            outputMatrixTime{icell}(istim, irep, :) = conv(squeeze(outputMatrixTime{icell}(istim, irep, :)),wingauss,'same');
            outputMatrixSpont{icell}(istim, irep) = sum(st>startTime - (endTime - startTime) & st<startTime);
        end
    end
%     for pos = 1:size(outputMatrixTime{icell},1)
%         subplot(1,size(outputMatrixTime{icell},1),pos);
%         timeX = linspace(-20*0.05,40*0.05,size(outputMatrixTime{icell},3));
%         plot(timeX,squeeze(mean(outputMatrixTime{icell}(pos,:,:),2))/max(outputMatrixTime{icell}(:)));
%         set(gca,'Xlim',[timeX(1) timeX(end)],'Ylim',[0 1])
%     end
%     set(gcf,'Name',['cell #' num2str(icell)]);    
%     pause
end
    
end