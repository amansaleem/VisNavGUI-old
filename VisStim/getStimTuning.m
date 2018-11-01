function [tuning, stim, p, cellIDs] = getStimTuning(animal, iseries, iexplist, igroup, iaddinfo)

global DIRS

if nargin < 3
    iexplist = [];
end
if nargin<4
    igroup = 8;
end
if nargin<5
    iaddinfo = 'V1';
end
tuning = [];
stim = [];
p = [];
cellIDs = [];
if isempty(iexplist)
    iexplist = 1:10;
end
ivs = 0;
for ixp = 1:numel(iexplist)
    try
    iexp = iexplist(ixp);
    dirname = fullfile(DIRS.Cerebus,animal,num2str(iseries),num2str(iexp));
    if exist(dirname,'dir')

        [outputMatrixTime,outputMatrix, outputMatrixSpont, p, cellIDs] = getStimSpiketimes(animal, iseries, iexp, igroup, iaddinfo);
        ivs = ivs + 1;
%         stim(ixp) = p.pars(p.activepars{1},:);
        kstr = strfind(p.description,'(');
        tuning(ivs).name = p.description(1:kstr-2);
        tuning(ivs).time = (-20:40)*50;
        tuning(ivs).Xpos = p.pars(10,1:end-1);
        tuning(ivs).Ypos = p.pars(11,1:end-1);
        tuning(ivs).rep = zeros(length(outputMatrixTime),size(p.seqnums,1),size(outputMatrixTime{1},3));
        tuning(ivs).globalrep = zeros(length(outputMatrixTime),size(outputMatrixTime{1},3));
        tuning(ivs).repsem = zeros(length(outputMatrixTime),size(p.seqnums,1),size(outputMatrixTime{1},3));
        tuning(ivs).mean = zeros(length(outputMatrix),size(p.seqnums,1));
        tuning(ivs).meanNorm = zeros(length(outputMatrix),size(p.seqnums,1));
        tuning(ivs).median = zeros(length(outputMatrix),size(p.seqnums,1));
        tuning(ivs).zscore = zeros(length(outputMatrix),size(p.seqnums,1));
        tuning(ivs).sem  = zeros(length(outputMatrix),size(p.seqnums,1));
        tuning(ivs).active = zeros(1,length(outputMatrix));

        for icell = 1:length(outputMatrix)
            if ~isempty(outputMatrix{icell})
                tuning(ivs).rep(icell,:,:) = mean(outputMatrixTime{icell},2);%mean(outputMatrix{icell}') - mean(outputMatrixSpont{icell}(:));
                tuning(ivs).repsem(icell,:,:)  = std(outputMatrixTime{icell},0,2)/(size(outputMatrixTime{icell},2)^0.5);%std(outputMatrixSpont{icell}(:));
                tuning(ivs).globalrep(icell,:) = mean(mean(outputMatrixTime{icell},1),2);
                tuning(ivs).mean(icell,:,:) = mean(outputMatrix{icell}');
                tuning(ivs).median(icell,:,:) = median(outputMatrix{icell}');
                tuning(ivs).sem(icell,:,:)  = std(mean(outputMatrix{icell}'));
                tuning(ivs).zscore(icell,:,:) = (mean(outputMatrix{icell}') - mean(outputMatrixSpont{icell}(:)))./max(0.01,tuning(ivs).sem(icell,:,:));
                tuning(ivs).meanNorm(icell,:,:) = (mean(outputMatrix{icell}') - mean(outputMatrixSpont{icell}(:)))/max(0.01,mean(outputMatrixSpont{icell}(:)));
                tuning(ivs).active(icell) = 1;
            end
        end
    end
    catch
        warning(['something went wrong with VS #' num2str(iexp)]);
    end
end
end