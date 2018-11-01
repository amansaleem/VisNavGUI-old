function popresLFP = BatchLFP_VS(popresLFP)
SetDirs;

expt = getExperimentList;
nanimal = numel(expt);

for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        if expt(ianimal).CA1{iseries} == 1 && (expt(ianimal).V1{iseries} == 0)
            res = preprocessLFP_VS(animalname, expt(ianimal).series{iseries}, [], 'CA1');
            popresLFP.repVS_CA1{ianimal,iseries} = res.LFPfilt_VS;
        end
        if expt(ianimal).CA1{iseries} == 0 && (expt(ianimal).V1{iseries} == 1)
            res = preprocessLFP_VS(animalname, expt(ianimal).series{iseries}, [], 'V1');
            popresLFP.repVS_V1{ianimal,iseries} = res.LFPfilt_VS;
        end
        if expt(ianimal).CA1{iseries} == 1 && (expt(ianimal).V1{iseries} == 1)
            resCA1 = preprocessLFP_VS(animalname, expt(ianimal).series{iseries}, [], 'CA1');
            resV1 = preprocessLFP_VS(animalname, expt(ianimal).series{iseries}, [], 'V1');
            popresLFP.repVS_CA1{ianimal,iseries} = resCA1.LFPfilt_VS;
            popresLFP.repVS_V1{ianimal,iseries} = resV1.LFPfilt_VS;
        end
        
        resCA1 = [];
        resV1 = [];
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