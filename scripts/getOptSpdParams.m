function [nspdbins, smth_spd] = getOptSpdParams(EXP)
global DIRS
nProbe = numel(unique(EXP.CellInfo.Probe));
nspdbins = zeros(1,nProbe);
smth_spd = zeros(1,nProbe);
paramsfilename = 'decoderVisSpdParams';%'decoderRunSpdParams';%['decoderSpdParams_' model];
if ~iscell(EXP.series)
    paramsfilepath = [DIRS.multichanspikes filesep EXP.animal filesep EXP.series filesep paramsfilename '.mat'];
else
    paramsfilepath = [DIRS.multichanspikes filesep EXP.animal filesep EXP.series{1} filesep paramsfilename '.mat'];
end
if exist(paramsfilepath,'file')
    S = load(paramsfilepath);
    traingain = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
    maxPostparams = cell(1,nProbe);
    for iprobe = 1:nProbe
        maxPostparams{iprobe} = squeeze(max(S.ParamsMeanErr{1}(iprobe,:,1:2,:,traingain,:),[],6));
        [imax,jmax] = find(maxPostparams{iprobe} == max(maxPostparams{iprobe}(:)));
        nspdbins(iprobe) = S.nspdbinslist(2);%(imax);
        smth_spd(iprobe) = S.smth_spdlist(jmax);
    end    
else
    nspdbins(:) = 5;
    smth_spd(:) = 150;
end
% ParamsMeanErr
% ParamsCorr
% nspdbinslist
% smth_spdlist
% latencylist