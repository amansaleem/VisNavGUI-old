function [ParamsMaxPost] = TestDecodingParams2(EXP, nspdbinslist, smth_spdlist, latcorrectionlist)

traincont = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
traingain = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
trainroomlength = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
trainoutcome = find(EXP.SubsetVal.outcome == 2);
type = 'mean';
nruns = 1;
goodidx_th = 30;
smth_win = 50;%20;%
numbins = 100;
thetachannel = 34;
nthetabins = 1;%1;%6;%
kfold = 20;
FoptiSmooth =0;
FGoodcluster = 1;
FUnsortedcluster = 1;
FMUAcluster = 1;
maxRate = 8;
zth = 1.96;
SSImin = 0;

ngain = numel(unique(EXP.data.es.gain(~isnan(EXP.data.es.gain))));
ParamsMaxPost = [];

for latcorr = 1:size(latcorrectionlist,2)
    for smthspd = 1:size(smth_spdlist,2)
        for spdbins = 1:size(nspdbinslist,2)
            smth_spd = smth_spdlist(:,smthspd);
            nspdbins = nspdbinslist(:,spdbins);
            latcorrection = latcorrectionlist(:,latcorr);
            
            if nspdbins ~= 1 || smthspd == 1
                EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,...
                    'smth_win', smth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'smth_spd', smth_spd, 'latcorrection', latcorrection,...
                    'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                    'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin);
                nProbe = numel(EXP.Bayes.PosError0);
                for iprobe = 1:nProbe
                    for g = 1:ngain
                        tidx = EXP.getSubsets(traincont,g,trainroomlength,trainoutcome,false);
                        Err =  EXP.Bayes.PosError0{iprobe};
                        maxtol = 1;%0.1;
                        [~, Xpred, ErrVec] = getcleandecidx(EXP.Bayes.Posterior0{iprobe}, EXP.Bayes.X,EXP.Bayes.goodidx_th,maxtol,EXP.data.es.trialID);
                        if isempty(ParamsMaxPost)
                            ParamsMaxPost = zeros(nProbe,size(latcorrectionlist,2),size(nspdbinslist,2),size(smth_spdlist,2),ngain,size(Err,2));
                        end
                        if ~isempty(Err)
                            for i = 1:100
                                MeanErr = max(nanmean(Err(tidx & EXP.Bayes.X == i,:),1),[],2);
                                ParamsMaxPost(iprobe,latcorr,spdbins,smthspd,g,i) = MeanErr;
                            end
                        end
                    end
                end                
            end
        end
    end
end
end