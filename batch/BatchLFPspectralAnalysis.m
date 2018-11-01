function BatchLFPspectralAnalysis
SetDirs;

expt = getExperimentList;
datadir = DIRS.data;
          
popresCorr.Tsmthwin = 50;%250;%250;%150;%300;%40;%120;%50
popresCorr.Xsmthwin = 4;%1;%
popresCorr.SpdSmthWin = popresCorr.Tsmthwin;
popresCorr.SpeedThreshold = 5;
popresCorr.nspeedbins = 3;
popresCorr.nthetaphsbins = 0;%1;%
popresCorr.cellstr = 'goodonly';%'All_50bins';%'goodonly';%'goodonly_unwrapped';%'goodonly';%'All';%
filesuffix_EXP = ['Twin' num2str(popresCorr.Tsmthwin) '_' 'Xwin' num2str(popresCorr.Xsmthwin) '_' 'spdth' num2str(popresCorr.SpeedThreshold) '_' num2str(popresCorr.nthetaphsbins) 'thetabins' '_' popresCorr.cellstr];
disp(filesuffix_EXP);

popresCorr.sampleRate = 60;
popresCorr.nSpdbins = 5;%1;%

lambdaSmooth = 2;
corrmaxlag = 180;
nanimal = numel(expt);

contval = [0.2:0.05:0.9];%[0.2 0.3 0.4];%[0.8 0.9];%
outvalcorr = 2;%[0 1 2 3 4 5];%[0 1 2 3 4];%

for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        dDIRname = ['D:\DATA\batch'  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];
        savedfile = [dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' 'goodonly' '_LFPSpikeCoherence.mat'];

        S = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' filesuffix_EXP '_cellProperties.mat']);
        Fgoodunits = S.CellInfo.Goodcluster;
        Finterneuron = S.CellInfo.Finterneuron;
        if expt(ianimal).CA1{iseries} == 1 && (expt(ianimal).V1{iseries} == 0)
            resCA1 = preprocessLFP(animalname, expt(ianimal).series{iseries}, expt(ianimal).exp{iseries}, 'CA1');
            resCA1.LFPfilt = resCA1.LFPfilt(:,resCA1.CA1chref);
%             resCA1_VS = preprocessLFP(animalname, expt(ianimal).series{iseries}, 1:10, 'CA1');
%             resCA1_VS.LFPfilt = resCA1_VS.LFPfilt(:,resCA1_VS.CA1chref);
            
            resCA1.Fgoodunit = Fgoodunits(S.CellInfo.Probe == 1);
%             resCA1_VS.Fgoodunit = Fgoodunits(S.CellInfo.Probe == 1);
            
            resCA1V1 = frequencyAnalysis(resCA1, []);
%             if ~isempty(resCA1_VS.LFPfilt)
%                 resCA1V1_VS = frequencyAnalysis(resCA1_VS, []);
%             else
%                 resCA1V1_VS = [];
%             end
            
            resCA1V1(1).Fgoodunit = Fgoodunits(S.CellInfo.Probe == 1) & ~Finterneuron(S.CellInfo.Probe == 1);
%             resCA1V1_VS(1).Fgoodunit = Fgoodunits(S.CellInfo.Probe == 1) & ~Finterneuron(S.CellInfo.Probe == 1);
            
%             popresLFP.f = resCA1V1(1).f;
%             for g = [2 1 3]
%                 popresLFP.CA1lfp_Spec{ianimal,iseries,g} = resCA1V1(1).lfp_Spec{g};
%                 popresLFP.meanBallSpd{ianimal,iseries,g} = resCA1V1(1).meanBallSpd{g};
%             end
%             popresLFP.repVS_CA1{ianimal,iseries} = resCA1V1(1).repVS;
%             popresLFP.StimType{ianimal,iseries} = resCA1V1(1).StimType;
        end
        if expt(ianimal).CA1{iseries} == 0 && (expt(ianimal).V1{iseries} == 1)
            resV1 = preprocessLFP(animalname, expt(ianimal).series{iseries}, expt(ianimal).exp{iseries}, 'V1');
            resV1.LFPfilt = resV1.LFPfilt(:,1);
%             resV1_VS = preprocessLFP(animalname, expt(ianimal).series{iseries}, 1:10, 'V1');
%             resV1_VS.LFPfilt = resV1_VS.LFPfilt(:,1);
            
            resV1.Fgoodunit = Fgoodunits(S.CellInfo.Probe == 2);
%             resV1_VS.Fgoodunit = Fgoodunits(S.CellInfo.Probe == 2);
            
            resCA1V1 = frequencyAnalysis([], resV1);
%             if ~isempty(resV1_VS.LFPfilt)
%                 resCA1V1_VS = frequencyAnalysis([], resV1_VS);
%             else
%                 resCA1V1_VS = [];
%             end
            
            resCA1V1(2).Fgoodunit = Fgoodunits(S.CellInfo.Probe == 2);
%             resCA1V1_VS(2).Fgoodunit = Fgoodunits(S.CellInfo.Probe == 2);
            
%             popresLFP.f = resCA1V1(2).f;
%             for g = [2 1 3]
%                 popresLFP.V1lfp_Spec{ianimal,iseries,g} = resCA1V1(2).lfp_Spec{g};
%                 popresLFP.meanBallSpd{ianimal,iseries,g} = resCA1V1(2).meanBallSpd{g};
%             end
%             popresLFP.repVS_V1{ianimal,iseries} = resCA1V1(2).repVS;
%             popresLFP.StimType{ianimal,iseries} = resCA1V1(2).StimType;
        end
        if expt(ianimal).CA1{iseries} == 1 && (expt(ianimal).V1{iseries} == 1)
            resCA1 = preprocessLFP(animalname, expt(ianimal).series{iseries}, expt(ianimal).exp{iseries}, 'CA1');
            resCA1.LFPfilt = resCA1.LFPfilt(:,resCA1.CA1chref);
%             resCA1_VS = preprocessLFP(animalname, expt(ianimal).series{iseries}, 1:10, 'CA1');
%             resCA1_VS.LFPfilt = resCA1_VS.LFPfilt(:,resCA1_VS.CA1chref);
            
            resCA1.Fgoodunit = Fgoodunits(S.CellInfo.Probe == 1);% & ~Finterneuron(S.CellInfo.Probe == 1);
%             resCA1_VS.Fgoodunit = Fgoodunits(S.CellInfo.Probe == 1);% & ~Finterneuron(S.CellInfo.Probe == 1);
            
            resV1 = preprocessLFP(animalname, expt(ianimal).series{iseries}, expt(ianimal).exp{iseries}, 'V1');
            resV1.LFPfilt = resV1.LFPfilt(:,1);
%             resV1_VS = preprocessLFP(animalname, expt(ianimal).series{iseries}, 1:10, 'V1');
%             resV1_VS.LFPfilt = resV1_VS.LFPfilt(:,1);
            
            resV1.Fgoodunit = Fgoodunits(S.CellInfo.Probe == 2);
%             resV1_VS.Fgoodunit = Fgoodunits(S.CellInfo.Probe == 2);
            
            resCA1V1 = frequencyAnalysis(resCA1, resV1);
%             if ~isempty(resV1_VS.LFPfilt)
%                 resCA1V1_VS = frequencyAnalysis(resCA1_VS, resV1_VS);
%             else
%                 resCA1V1_VS = [];
%             end
            
            resCA1V1(1).Fgoodunit = Fgoodunits(S.CellInfo.Probe == 1) & ~Finterneuron(S.CellInfo.Probe == 1);
            resCA1V1(2).Fgoodunit = Fgoodunits(S.CellInfo.Probe == 2);
%             resCA1V1_VS(1).Fgoodunit = Fgoodunits(S.CellInfo.Probe == 1) & ~Finterneuron(S.CellInfo.Probe == 1);
%             resCA1V1_VS(2).Fgoodunit = Fgoodunits(S.CellInfo.Probe == 2);
            
%             popresLFP.f = resCA1V1(1).f;
%             for g = [2 1 3]
%                 popresLFP.CA1lfp_Spec{ianimal,iseries,g} = resCA1V1(1).lfp_Spec{g};
%                 popresLFP.V1lfp_Spec{ianimal,iseries,g} = resCA1V1(2).lfp_Spec{g};
%                 
%                 popresLFP.CA1spk_SpecCh{ianimal,iseries,g} = resCA1V1(1).spk_SpecCh{g};
%                 popresLFP.CA1spk_Spec{ianimal,iseries,g} = resCA1V1(1).spk_Spec{g};
%                 popresLFP.CA1spk_CohSpec{ianimal,iseries,g} = resCA1V1(1).spk_CohSpec{g};
%                 popresLFP.CA1spk_PhsCohSpec{ianimal,iseries,g} = resCA1V1(1).spk_PhsCohSpec{g};
%                 popresLFP.CA1spk_CohSpecCh{ianimal,iseries,g} = resCA1V1(1).spk_CohSpecCh{g};
%                 popresLFP.CA1spk_PhsCohSpecCh{ianimal,iseries,g} = resCA1V1(1).spk_PhsCohSpecCh{g};
%                 
%                 popresLFP.V1spk_SpecCh{ianimal,iseries,g} = resCA1V1(2).spk_SpecCh{g};
%                 popresLFP.V1spk_Spec{ianimal,iseries,g} = resCA1V1(2).spk_Spec{g};
%                 popresLFP.V1spk_CohSpec{ianimal,iseries,g} = resCA1V1(2).spk_CohSpec{g};
%                 popresLFP.V1spk_PhsCohSpec{ianimal,iseries,g} = resCA1V1(2).spk_PhsCohSpec{g};
%                 popresLFP.V1spk_CohSpecCh{ianimal,iseries,g} = resCA1V1(2).spk_CohSpecCh{g};
%                 popresLFP.V1spk_PhsCohSpecCh{ianimal,iseries,g} = resCA1V1(2).spk_PhsCohSpecCh{g};
%                 
%                 popresLFP.POthetaCA1{ianimal,iseries,g} = resCA1V1(1).POthetaAve{g};
%                 popresLFP.POthetaV1{ianimal,iseries,g} = resCA1V1(2).POthetaAve{g};
%                 
%                 popresLFP.meanBallSpd{ianimal,iseries,g} = resCA1V1(2).meanBallSpd{g};
%             end
%             popresLFP.CA1lfp_SpecAll{ianimal,iseries} = resCA1V1(1).lfp_SpecAll;
%             popresLFP.V1lfp_SpecAll{ianimal,iseries} = resCA1V1(2).lfp_SpecAll;
%             
%             popresLFP.CA1spk_SpecChAll{ianimal,iseries} = resCA1V1(1).spk_SpecChAll;
%             popresLFP.CA1spk_SpecAll{ianimal,iseries} = resCA1V1(1).spk_SpecAll;
%             popresLFP.CA1spk_CohSpecAll{ianimal,iseries} = resCA1V1(1).spk_CohSpecAll;
%             popresLFP.CA1spk_PhsCohSpecAll{ianimal,iseries} = resCA1V1(1).spk_PhsCohSpecAll;
%             popresLFP.CA1spk_CohSpecChAll{ianimal,iseries} = resCA1V1(1).spk_CohSpecChAll;
%             popresLFP.CA1spk_PhsCohSpecChAll{ianimal,iseries} = resCA1V1(1).spk_PhsCohSpecChAll;
%             
%             popresLFP.V1spk_SpecChAll{ianimal,iseries} = resCA1V1(2).spk_SpecChAll;
%             popresLFP.V1spk_SpecAll{ianimal,iseries} = resCA1V1(2).spk_SpecAll;
%             popresLFP.V1spk_CohSpecAll{ianimal,iseries} = resCA1V1(2).spk_CohSpecAll;
%             popresLFP.V1spk_PhsCohSpecAll{ianimal,iseries} = resCA1V1(2).spk_PhsCohSpecAll;
%             popresLFP.V1spk_CohSpecChAll{ianimal,iseries} = resCA1V1(2).spk_CohSpecChAll;
%             popresLFP.V1spk_PhsCohSpecChAll{ianimal,iseries} = resCA1V1(2).spk_PhsCohSpecChAll;
%             
%             popresLFP.CA1lfp_SpecAll_VS{ianimal,iseries} = resCA1V1_VS(1).lfp_SpecAll;
%             popresLFP.V1lfp_SpecAll_VS{ianimal,iseries} = resCA1V1_VS(2).lfp_SpecAll;
%             
%             popresLFP.CA1spk_SpecChAll_VS{ianimal,iseries} = resCA1V1_VS(1).spk_SpecChAll;
%             popresLFP.CA1spk_SpecAll_VS{ianimal,iseries} = resCA1V1_VS(1).spk_SpecAll;
%             popresLFP.CA1spk_CohSpecAll_VS{ianimal,iseries} = resCA1V1_VS(1).spk_CohSpecAll;
%             popresLFP.CA1spk_PhsCohSpecAll_VS{ianimal,iseries} = resCA1V1_VS(1).spk_PhsCohSpecAll;
%             popresLFP.CA1spk_CohSpecChAll_VS{ianimal,iseries} = resCA1V1_VS(1).spk_CohSpecChAll;
%             popresLFP.CA1spk_PhsCohSpecChAll_VS{ianimal,iseries} = resCA1V1_VS(1).spk_PhsCohSpecChAll;
%             
%             popresLFP.V1spk_SpecChAll_VS{ianimal,iseries} = resCA1V1_VS(2).spk_SpecChAll;
%             popresLFP.V1spk_SpecAll_VS{ianimal,iseries} = resCA1V1_VS(2).spk_SpecAll;
%             popresLFP.V1spk_CohSpecAll_VS{ianimal,iseries} = resCA1V1_VS(2).spk_CohSpecAll;
%             popresLFP.V1spk_PhsCohSpecAll_VS{ianimal,iseries} = resCA1V1_VS(2).spk_PhsCohSpecAll;
%             popresLFP.V1spk_CohSpecChAll_VS{ianimal,iseries} = resCA1V1_VS(2).spk_CohSpecChAll;
%             popresLFP.V1spk_PhsCohSpecChAll_VS{ianimal,iseries} = resCA1V1_VS(2).spk_PhsCohSpecChAll;
%             
%             popresLFP.repVS_CA1{ianimal,iseries} = resCA1V1(1).repVS;
%             popresLFP.repVS_V1{ianimal,iseries} = resCA1V1(2).repVS;
%             popresLFP.StimType{ianimal,iseries} = resCA1V1(2).StimType;
        end
        
        save(savedfile,'resCA1V1','-v7.3');
%         save(savedfile,'resCA1V1','resCA1V1_VS','-v7.3');
        
        resCA1 = [];
        resV1 = [];
        resCA1V1 = [];
        resCA1_VS = [];
        resV1_VS = [];
        resCA1V1_VS = [];
    end
end
end

% to save popresLFP, run the following in the command window:
%save(['D:\DATA\batch\All' filesep 'popresLFP.mat'], 'popresLFP','-v7.3');

function mat_out = smooth2D(mat,lambdaSmooth)
mat(isnan(mat)) = 0;
G = smooth1D(repmat(mat,3,3),lambdaSmooth);
H = smooth1D(G',lambdaSmooth)';
mat_out = H(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
end

%to plot single sessions
% figure;
% for ianimal = 1:10
% for iseries = 1:6
% subplot(2,4,1);
% if ~isempty(popresLFP.SpecCA1{ianimal,iseries,2})
% plot(popresLFP.f,log(popresLFP.SpecCA1{ianimal,iseries,2}))
% set(gca,'Xscale','log');
% end
% subplot(2,4,5);
% if ~isempty(popresLFP.SpecV1{ianimal,iseries,2})
% plot(popresLFP.f,popresLFP.SpecV1{ianimal,iseries,2})
% set(gca,'Xscale','log');
% end
% subplot(4,4,[6,10]);
% if ~isempty(popresLFP.CohSpec{ianimal,iseries,2})
% plot(popresLFP.f,popresLFP.CohSpec{ianimal,iseries,2})
% set(gca,'Xscale','log');
% end
% subplot(4,4,3);
% if ~isempty(popresLFP.POthetaCA1{ianimal,iseries,2})
% imagesc(popresLFP.POthetaCA1{ianimal,iseries,2}')
% end
% subplot(4,4,[7,11]);
% if ~isempty(popresLFP.POthetaV1{ianimal,iseries,2})
% imagesc(LFP2iCSD(popresLFP.POthetaV1{ianimal,iseries,2}', 2e-6:2e-6:(32*2e-6)))
% end
% subplot(4,4,[8,12]);
% vsmap = 0;
% if ~isempty(popresLFP.repVS_V1{ianimal,iseries})
% for ivs = 1:numel(popresLFP.repVS_V1{ianimal,iseries})
% if ~isempty(popresLFP.repVS_V1{ianimal,iseries}{ivs})
% if size(popresLFP.repVS_V1{ianimal,iseries}{ivs},2) == 32
% vsmap = vsmap + popresLFP.repVS_V1{ianimal,iseries}{ivs};
% end
% end
% end
% end
% if numel(vsmap) > 1
% imagesc(LFP2iCSD(conv2(vsmap,ones(200,1),'same')', 2e-6:2e-6:(32*2e-6)))
% end
% pause
% end
% end