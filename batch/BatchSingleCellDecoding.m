function resdec = BatchSingleCellDecoding

expt = getExperimentList;
area_str = 'CA1V1';
nanimals = numel(expt);%1;

Tsmthwin = 120;%300;%120;%
Xsmthwin = 4;
SpeedThreshold = 5;
nspeedbins = 5;
nthetaphsbins = 1;%0;%1;%
cellstr = 'All';
filesuffix_EXP = ['Twin' num2str(Tsmthwin) '_' 'Xwin' num2str(Xsmthwin) '_' 'spdth' num2str(SpeedThreshold) '_' num2str(nspeedbins) 'speedbins' '_' num2str(nthetaphsbins) 'thetabins' '_' cellstr];

nprobe = 2;
ngain = 4;
outcomelist = [1 2 3];
noutcome = numel(outcomelist);

resdec.SSI_avedec = cell(nprobe,ngain,noutcome);
resdec.SSI_maxdec = cell(nprobe,ngain,noutcome);
resdec.SSI_avedecAll = cell(nprobe,ngain,noutcome);
resdec.SSI_maxdecAll = cell(nprobe,ngain,noutcome);
resdec.SSI_avedecCross = cell(nprobe,ngain,noutcome);
resdec.SSI_maxdecCross = cell(nprobe,ngain,noutcome);
resdec.SSI_traj = cell(nprobe,ngain,noutcome);
resdec.SpatialInfoPerSpike_avedec = cell(nprobe,ngain,noutcome);
resdec.SpatialInfoPerSpike_maxdec = cell(nprobe,ngain,noutcome);
resdec.SpatialInfoPerSpike_avedecAll = cell(nprobe,ngain,noutcome);
resdec.SpatialInfoPerSpike_maxdecAll = cell(nprobe,ngain,noutcome);
resdec.SpatialInfoPerSpike_avedecCross = cell(nprobe,ngain,noutcome);
resdec.SpatialInfoPerSpike_maxdecCross = cell(nprobe,ngain,noutcome);
resdec.SpatialInfoPerSpike_traj = cell(nprobe,ngain,noutcome);
resdec.fieldpos = cell(1,nprobe);
resdec.bestchan = cell(1,nprobe);
resdec.animal = cell(1,nprobe);

for ianimal = 1:nanimals
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
            dDIRname = ['D:\DATA\batch'  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];%[DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];
            if exist([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' filesuffix_EXP '_DecodingPopCell.mat'],'file')
                S = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' filesuffix_EXP '_DecodingPopCell.mat']);
                allcont = size(S.CellInfo.SSI_avedec,1);
                for iprobe = 1:size(S.resAll.DecCells,2)
                    if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                        if size(S.CellInfo.SSI_avedec,2)<4
                            gainlist = [1 2 0 3];
                        else
                            gainlist = 1:4;
                        end
                        for g = 1:4
                            if gainlist(g) > 0
                                for o = 1:numel(outcomelist)
                                    resdec.SSI_avedec{iprobe,g,o} = [resdec.SSI_avedec{iprobe,g,o} S.CellInfo.SSI_avedec{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    resdec.SSI_maxdec{iprobe,g,o} = [resdec.SSI_maxdec{iprobe,g,o} S.CellInfo.SSI_maxdec{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    
                                    resdec.SSI_avedecAll{iprobe,g,o} = [resdec.SSI_avedecAll{iprobe,g,o} S.CellInfo.SSI_avedecAll{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    resdec.SSI_maxdecAll{iprobe,g,o} = [resdec.SSI_maxdecAll{iprobe,g,o} S.CellInfo.SSI_maxdecAll{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    
                                    resdec.SSI_avedecCross{iprobe,g,o} = [resdec.SSI_avedecCross{iprobe,g,o} S.CellInfo.SSI_avedecCross{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    resdec.SSI_maxdecCross{iprobe,g,o} = [resdec.SSI_maxdecCross{iprobe,g,o} S.CellInfo.SSI_maxdecCross{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    
                                    resdec.SSI_traj{iprobe,g,o} = [resdec.SSI_traj{iprobe,g,o} S.CellInfo.SSI_traj{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    
                                    resdec.SpatialInfoPerSpike_avedec{iprobe,g,o} = [resdec.SpatialInfoPerSpike_avedec{iprobe,g,o} S.CellInfo.SpatialInfoPerSpike_avedec{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    resdec.SpatialInfoPerSpike_maxdec{iprobe,g,o} = [resdec.SpatialInfoPerSpike_maxdec{iprobe,g,o} S.CellInfo.SpatialInfoPerSpike_maxdec{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    
                                    resdec.SpatialInfoPerSpike_avedecAll{iprobe,g,o} = [resdec.SpatialInfoPerSpike_avedecAll{iprobe,g,o} S.CellInfo.SpatialInfoPerSpike_avedecAll{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    resdec.SpatialInfoPerSpike_maxdecAll{iprobe,g,o} = [resdec.SpatialInfoPerSpike_maxdecAll{iprobe,g,o} S.CellInfo.SpatialInfoPerSpike_maxdecAll{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    
                                    resdec.SpatialInfoPerSpike_avedecCross{iprobe,g,o} = [resdec.SpatialInfoPerSpike_avedecCross{iprobe,g,o} S.CellInfo.SpatialInfoPerSpike_avedecCross{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    resdec.SpatialInfoPerSpike_maxdecCross{iprobe,g,o} = [resdec.SpatialInfoPerSpike_maxdecCross{iprobe,g,o} S.CellInfo.SpatialInfoPerSpike_maxdecCross{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                    
                                    resdec.SpatialInfoPerSpike_traj{iprobe,g,o} = [resdec.SpatialInfoPerSpike_traj{iprobe,g,o} S.CellInfo.SpatialInfoPerSpike_traj{allcont, gainlist(g), 1, outcomelist(o)}(S.resAll.DecCells{iprobe})];
                                end
                            else
                                for o = 1:numel(outcomelist)
                                    resdec.SSI_avedec{iprobe,g,o} = [resdec.SSI_avedec{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    resdec.SSI_maxdec{iprobe,g,o} = [resdec.SSI_maxdec{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    
                                    resdec.SSI_avedecAll{iprobe,g,o} = [resdec.SSI_avedecAll{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    resdec.SSI_maxdecAll{iprobe,g,o} = [resdec.SSI_maxdecAll{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    
                                    resdec.SSI_avedecCross{iprobe,g,o} = [resdec.SSI_avedecCross{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    resdec.SSI_maxdecCross{iprobe,g,o} = [resdec.SSI_maxdecCross{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    
                                    resdec.SSI_traj{iprobe,g,o} = [resdec.SSI_traj{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    
                                    resdec.SpatialInfoPerSpike_avedec{iprobe,g,o} = [resdec.SpatialInfoPerSpike_avedec{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    resdec.SpatialInfoPerSpike_maxdec{iprobe,g,o} = [resdec.SpatialInfoPerSpike_maxdec{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    
                                    resdec.SpatialInfoPerSpike_avedecAll{iprobe,g,o} = [resdec.SpatialInfoPerSpike_avedecAll{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    resdec.SpatialInfoPerSpike_maxdecAll{iprobe,g,o} = [resdec.SpatialInfoPerSpike_maxdecAll{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    
                                    resdec.SpatialInfoPerSpike_avedecCross{iprobe,g,o} = [resdec.SpatialInfoPerSpike_avedecCross{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    resdec.SpatialInfoPerSpike_maxdecCross{iprobe,g,o} = [resdec.SpatialInfoPerSpike_maxdecCross{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                    
                                    resdec.SpatialInfoPerSpike_traj{iprobe,g,o} = [resdec.SpatialInfoPerSpike_traj{iprobe,g,o} NaN(1,sum(S.resAll.DecCells{iprobe}))];
                                end
                            end
                        end
                        [~, pos] = max(S.CellInfo.field_traj{allcont, 2, 1, 3}(S.resAll.DecCells{iprobe},:),[],2);
                        resdec.fieldpos{iprobe} = [resdec.fieldpos{iprobe} pos(:)'];
                        resdec.bestchan{iprobe} = [resdec.bestchan{iprobe} S.CellInfo.bestchan(S.resAll.DecCells{iprobe})];
                        resdec.animal{iprobe} = [resdec.animal{iprobe} ianimal*ones(1,numel(S.CellInfo.SSI_avedec{allcont, 1, 1, 3}(S.resAll.DecCells{iprobe})))];
                    end
                end
            end
        end
    end
end
end