function [output, es] = calculateCoherence(animal, iseries, iexp, chnA, chnB, groupID, addInfo)

if nargin<7
    addInfo = '1';
end

movingwin=[3 1];% [3 1];%[window_size window_shift] [3 1];
params.Fs = 1000;
params.tapers=[5 9];
params.fpass=[1 95];

[es, procSpec] = VR_LFP_power_safe(animal, iseries, iexp, 0, chnA, chnB);

for igroup = 1:length(groupID)
    display('Removing cell 30')
%     chans = getKwikSpikes_copy(animal, iseries, iexp, groupID(igroup), addInfo);
    chans = getKwikSpikes(animal, iseries, iexp, groupID(igroup), addInfo);
    allCells = [];
    for icell = 1:length(chans)
        numSpikes(icell) = length(chans(icell).spiketimes);
        allCells = [allCells' chans(icell).spiketimes']';
        output(igroup).cellIDs{icell} = chans(icell).id;
    end
    
    % coh = zeros(length(es.freq),sum(numSpikes>100));
    % coh_B = zeros(length(es.freq),sum(numSpikes>100));
    
    idx = 1;
    
    for icell = 1:length(chans)
        if numSpikes(icell)>100
            [output(igroup).A.C{icell},output(igroup).A.phi{icell},output(igroup).A.SLfpSpikes{icell},output(igroup).A.SLfp{icell},output(igroup).A.SSpikes{icell},t,f]=cohgramcpt(procSpec.ChnA,chans(icell).spiketimes,movingwin,params);
            [output(igroup).B.C{icell},output(igroup).B.phi{icell},output(igroup).B.SLfpSpikes{icell},output(igroup).B.SLfp{icell},output(igroup).B.SSpikes{icell}]=cohgramcpt(procSpec.ChnB,chans(icell).spiketimes,movingwin,params);
            coh(:,idx) = nanmean(output(igroup).A.C{icell},1);
            coh_B(:,idx) = nanmean(output(igroup).B.C{icell},1);
            idx = idx + 1; display(num2str(icell));
        end
    end
    
    if length(allCells)>100
        [output(igroup).A.C_allCells,output(igroup).A.phi_allCells,output(igroup).A.SLfpSpikes_allCells,output(igroup).A.SLfp_allCells,output(igroup).A.SSpikes_allCells]...
            =cohgramcpt(procSpec.ChnB,allCells,movingwin,params);
        [output(igroup).B.C_allCells,output(igroup).B.phi_allCells,output(igroup).B.SLfpSpikes_allCells,output(igroup).B.SLfp_allCells,output(igroup).B.SSpikes_allCells]...
            =cohgramcpt(procSpec.ChnB,allCells,movingwin,params);
    end
    
    output(igroup).t = t;
    output(igroup).freq = f;
    output(igroup).iexp = iexp;
    output(igroup).coh_A = coh;
    output(igroup).coh_B = coh_B;
    output(igroup).chnA = chnA;
    output(igroup).chnB = chnB;
    output(igroup).groupID = groupID(igroup);
    
end