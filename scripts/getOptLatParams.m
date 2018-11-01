function Latcorrection = getOptLatParams(EXP)
global DIRS
nProbe = numel(unique(EXP.CellInfo.Probe));
Latcorrection = zeros(1,nProbe);
paramsfilename = 'decoderLatParams';%['decoderSpdParams_' model];
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
        maxPostparams{iprobe} = squeeze(max(S.ParamsMeanErr{1}(iprobe,:,:,:,traingain,:),[],6));
        [~,jmax] = find(maxPostparams{iprobe} == max(maxPostparams{iprobe}(:)));
        Latcorrection(iprobe) = S.latcorrectionlist(jmax);
    end    
else
    Latcorrection(:) = 0;
end
% ParamsMeanErr
% ParamsCorr
% nspdbinslist
% smth_spdlist
% latencylist