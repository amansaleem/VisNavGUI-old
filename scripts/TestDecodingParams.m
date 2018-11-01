function [xParams, ParamsMeanErr0, ParamsMeanErr1, ParamsMeanErr2, ParamsCorr, ParamsErr1, ParamsErr2, ParamsLickdistri, ParamsLickdistriXpred, ParamsPostbeforelick, ParamsMeanPostbeforelick] = TestDecodingParams(EXP, nspdbinslist, smth_spdlist, latcorrectionlist, alphalist, deltalist)

SpeedThreshold = 5;

traincont = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
traingain = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
trainroomlength = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
trainoutcome = find(EXP.SubsetVal.outcome == 2);
type = 'mean';
nruns = 1;
goodidx_th = 30;
Tsmth_win = 250;%35;%150;%20;%
Xsmth_win = 4;%
Flookuptable = 0;
numbins = 100;
thetachannel = 34;
nthetabins = 1;%1;%6;%
kfold = 20;
FoptiSmooth =0;
FGoodcluster = 1;
FUnsortedcluster = 0;
FMUAcluster = 0;
maxRate = inf;
zth = -inf;
SSImin = -inf;
speed_th = SpeedThreshold;

gainval = sort(unique(EXP.data.es.gain(~isnan(EXP.data.es.gain))),'ascend');
ngain = numel(unique(EXP.data.es.gain(~isnan(EXP.data.es.gain))))+1;
ParamsMeanErr0 = [];
ParamsMeanErr1 = [];
ParamsMeanErr2 = [];
ParamsCorr = [];
ParamsLickdistri = [];
ParamsLickdistriXpred = [];
ParamsMeanPostbeforelick = [];
ParamsPostbeforelick = [];
xParams = 0:0.1:100;

lambdaSmooth = 2;
for latcorr = 1:size(latcorrectionlist,2)
    smth_spd = Tsmth_win;%0;
    nspdbins = nspdbinslist(:,1);
    neyebins = nspdbinslist(:,1);
    latcorrection = latcorrectionlist(:,latcorr);
    alpha = 0;
    delta = 0;
    disp(['latency = ' num2str(latcorrection)]);
    EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'smth_spd', smth_spd, 'latcorrection', latcorrection,...
                'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin);
    
    for smthspd = 1:size(smth_spdlist,2)
        for spdbins = 1:size(nspdbinslist,2)
            for idelta = 1:size(deltalist,2)
                for ialpha = 1:size(alphalist,2)
                    smth_spd = smth_spdlist(:,smthspd);%0;%
                    nspdbins = nspdbinslist(:,spdbins);
                    latcorrection = latcorrectionlist(:,latcorr);
                    alpha = alphalist(:,ialpha);%0.5;
                    if alpha == 0
                        delta = 0;
                    else
                        delta = 1;
                    end
                    
%                     disp(['latency = ' num2str(delta) ' ; alpha = ' num2str(alpha)]);
                    
%                     EXP.Bayes.X = VisDistXmodel(EXP.data.es.(EXP.Bayes.varname),delta,alpha,latcorrection,EXP.data.es,EXP.Bayes.Tsmth_win,'boxcar_centered',EXP.Bayes.numBins);
                    EXP.Bayes.X = EXP.Bayes.X0;%EXP.Bayes.X = KalmanX2model(EXP.Bayes.X0,alpha,delta,1,0,EXP.data.es,EXP.Bayes.Tsmth_win,'boxcar_centered',EXP.Bayes.numBins);
                    
                    %             Tsmth_win = latcorrectionlist(:,latcorr);
                    %             smth_spd = Tsmth_win;
                    
                    %             if nspdbins ~= 1 || smthspd == 1%
                    %                 EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,...
                    %                     'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'smth_spd', smth_spd,...
                    %                     'latcorrection', latcorrection, 'alpha', alpha, 'delta', delta,...
                    %                     'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                    %                     'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin);
                    nProbe = numel(EXP.Bayes.PosError0);
                    for iprobe = 1:nProbe
                        for g = 1:ngain
                            allcontidx = find(EXP.SubsetVal.contrast > 0);
                            if g <= numel(gainval)
                                tidx = EXP.getSubsets(allcontidx,g,trainroomlength,trainoutcome,false);
                                %                             tidx = EXP.data.es.gain == gainval(g);%true(size(tidx));
                            else
                                tidx = EXP.getSubsets(allcontidx,1:numel(gainval),trainroomlength,trainoutcome,false);
                                tidx = true(size(tidx));
                            end
                            if sum(tidx) > 0
                                
                                Err =  EXP.Bayes.PosError0{iprobe};%sqrt(ErrVec.^2);%
                                Prange = size(EXP.Bayes.Posterior0{iprobe},2);
                                Xrange = max(EXP.Bayes.X);
                                for i = 1:Xrange
                                    Err(round(EXP.Bayes.X) == i,:) = circshift(EXP.Bayes.Posterior0{iprobe}(round(EXP.Bayes.X) == i,:),floor(Prange/2)-i,2);
                                end
                                
                                if isempty(ParamsMeanErr1)
                                    ParamsMeanErr0 = zeros(nProbe,size(latcorrectionlist,2),size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain,100,100);%size(Err,2));
                                    ParamsMeanErr1 = zeros(nProbe,size(latcorrectionlist,2),size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain,1001);%size(Err,2));
                                    ParamsMeanErr2 = zeros(nProbe,size(latcorrectionlist,2),size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain,1001);%size(Err,2));
                                    ParamsErr1 = zeros(nProbe,size(latcorrectionlist,2),size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain);
                                    ParamsErr2 = zeros(nProbe,size(latcorrectionlist,2),size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain);
                                end
                                
                                if ~isempty(Err)
                                    G = smooth1D(repmat(Err,1,3)',lambdaSmooth)';
                                    Err = G(:,size(Err,2)+1:2*size(Err,2));
                                    %                             for i = 1:nspdbins%Xrange
                                    Xidx = true(size(tidx));%EXP.Bayes.Spdbin == i;%EXP.Bayes.X == i;
                                    
                                    Xrange = max(EXP.Bayes.X);
                                    Prange = size(Err,2);
                                    MeanErr = zeros(Xrange,Prange);
                                    for xx = 1:Xrange
                                        MeanErr(xx,:) = nanmean(Err(tidx & Xidx & EXP.Bayes.X==xx,:),1);
                                    end
%                                     MeanErr = nanmean(MeanErr,1);
                                    
                                    %                                 MeanErr = nanmean(Err(tidx & Xidx,:),1);
%                                     vec1 = interp1(1:100,MeanErr,xParams,'spline');
                                    ParamsMeanErr0(iprobe,latcorr,spdbins,smthspd,ialpha,idelta,g,:,:) = MeanErr;%vec1;
                                    %                             end
                                end
                                
%                                 maxtol = 1;%0.1;%
%                                 Err =  EXP.Bayes.PosError0{iprobe};%sqrt(ErrVec.^2);%
%                                 [~, Xpred, ErrVec] = getcleandecidx(EXP.Bayes.Posterior0{iprobe}, EXP.Bayes.X,EXP.Bayes.goodidx_th,maxtol,EXP.data.es.trialID);
%                                 Prange = size(EXP.Bayes.Posterior0{iprobe},2);
%                                 Xrange = max(EXP.Bayes.X);
%                                 for i = 1:Xrange
%                                     Err(round(Xpred) == i,:) = circshift(EXP.Bayes.Posterior0{iprobe}(round(Xpred) == i,:),floor(Prange/2)-i,2);
%                                 end
%                                 
%                                 if ~isempty(Err)
%                                     G = smooth1D(repmat(Err,1,3)',lambdaSmooth)';
%                                     Err = G(:,size(Err,2)+1:2*size(Err,2));
%                                     %                             for i = 1:nspdbins%Xrange
%                                     Xidx = true(size(tidx));%EXP.Bayes.Spdbin == i;%EXP.Bayes.X == i;
%                                     
%                                     Xrange = max(EXP.Bayes.X);
%                                     Prange = size(Err,2);
%                                     MeanErr = zeros(Xrange,Prange);
%                                     MeanErrVec = zeros(Xrange,1);
%                                     for xx = 1:Xrange
%                                         MeanErr(xx,:) = nanmean(Err(tidx & Xidx & EXP.Bayes.X==xx,:),1);
%                                         MeanErrVec(xx) = nanmean(ErrVec(tidx & Xidx & EXP.Bayes.X==xx));
%                                     end
%                                     MeanErr = nanmean(MeanErr,1);
%                                     
%                                     %                                 MeanErr = nanmean(Err(tidx & Xidx,:),1);
%                                     vec1 = interp1(1:100,MeanErr,xParams,'spline');
%                                     ParamsMeanErr1(iprobe,latcorr,spdbins,smthspd,ialpha,idelta,g,:) = vec1;
%                                     ParamsErr1(iprobe,latcorr,spdbins,smthspd,ialpha,idelta,g) = sqrt(nanmean(MeanErrVec.^2));%(nanmean(MeanErrVec));%
%                                     %                             end
%                                 end
%                                 
%                                 
%                                 maxtol = 0.1;%0.1;%
%                                 Err =  EXP.Bayes.PosError0{iprobe};%sqrt(ErrVec.^2);%
%                                 [~, Xpred, ErrVec] = getcleandecidx(EXP.Bayes.Posterior0{iprobe}, EXP.Bayes.X,EXP.Bayes.goodidx_th,maxtol,EXP.data.es.trialID);
%                                 Prange = size(EXP.Bayes.Posterior0{iprobe},2);
%                                 Xrange = max(EXP.Bayes.X);
%                                 for i = 1:Xrange
%                                     Err(round(Xpred) == i,:) = circshift(EXP.Bayes.Posterior0{iprobe}(round(Xpred) == i,:),floor(Prange/2)-i,2);
%                                 end
%                                 
%                                 if ~isempty(Err)
%                                     G = smooth1D(repmat(Err,1,3)',lambdaSmooth)';
%                                     Err = G(:,size(Err,2)+1:2*size(Err,2));
%                                     %                             for i = 1:nspdbins%Xrange
%                                     Xidx = true(size(tidx));%EXP.Bayes.Spdbin == i;%EXP.Bayes.X == i;
%                                     
%                                     Xrange = max(EXP.Bayes.X);
%                                     Prange = size(Err,2);
%                                     MeanErr = zeros(Xrange,Prange);
%                                     MeanErrVec = zeros(Xrange,1);
%                                     for xx = 1:Xrange
%                                         MeanErr(xx,:) = nanmean(Err(tidx & Xidx & EXP.Bayes.X==xx,:),1);
%                                         MeanErrVec(xx) = nanmean(ErrVec(tidx & Xidx & EXP.Bayes.X==xx));
%                                     end
%                                     MeanErr = nanmean(MeanErr,1);
%                                     
%                                     %                                 MeanErr = nanmean(Err(tidx & Xidx,:),1);
%                                     vec1 = interp1(1:100,MeanErr,xParams,'spline');
%                                     ParamsMeanErr2(iprobe,latcorr,spdbins,smthspd,ialpha,idelta,g,:) = vec1;
%                                     ParamsErr2(iprobe,latcorr,spdbins,smthspd,ialpha,idelta,g) = sqrt(nanmean(MeanErrVec.^2));%(nanmean(MeanErrVec));%
%                                     %                             end
%                                 end
                                if latcorr == 1
                                    if isempty(ParamsLickdistri)
                                        ParamsLickdistri = zeros(nProbe,1,size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain,101);
                                        ParamsLickdistriXpred = zeros(nProbe,1,size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain,101);
                                        maxNTau = 60;
                                        ParamsMeanPostbeforelick = zeros(nProbe,1,size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain,maxNTau,100);
                                        ParamsPostbeforelick = zeros(nProbe,1,size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain,maxNTau,100);
                                    end

                                    maxtol = 0.1;%0.1;%
                                    [~, Xpred, ErrVec] = getcleandecidx(EXP.Bayes.Posterior0{iprobe}, EXP.Bayes.X,EXP.Bayes.goodidx_th,maxtol,EXP.data.es.trialID);
                                    licks = EXP.data.es.firstgoodlick;
                                    if sum(tidx) > 0
                                        [lickdistri,x] = fast1Dmap(EXP.Bayes.X(tidx),licks(tidx),1,60,[],EXP.data.es.CircularMaze);
                                        %                                     lickdistri = interp1([x x(end) + x(2) - x(1)],[lickdistri' lickdistri(1)],xParams,'spline');
                                        lickdistri = lickdistri./sum(lickdistri);
                                        [lickdistriXpred,x] = fast1Dmap(floor(Xpred(tidx & ~isnan(Xpred))),licks(tidx  & ~isnan(Xpred)),1,60,[],EXP.data.es.CircularMaze);
                                        %                                     lickdistriXpred = interp1([x x(end) + x(2) - x(1)],[lickdistriXpred' lickdistriXpred(1)],xParams,'spline');
                                        lickdistriXpred = lickdistriXpred./sum(lickdistriXpred);
                                    end
                                    ParamsLickdistri(iprobe,latcorr,spdbins,smthspd,ialpha,idelta,g,:) = lickdistri;
                                    ParamsLickdistriXpred(iprobe,latcorr,spdbins,smthspd,ialpha,idelta,g,:) = lickdistriXpred;


                                    licksidx = find(EXP.data.es.firstgoodlick);
                                    for Ntau = 1:maxNTau%i.e. we average up to 1 second before
                                        beforelick = zeros(size(EXP.data.es.firstrewlick));
                                        meanbeforelick = zeros(size(EXP.data.es.firstrewlick));
                                        for tt = 0:10
                                            beforelick(max(1,licksidx - Ntau - tt)) = 1;
                                        end
                                        for tt = 0:Ntau
                                            meanbeforelick(max(1,licksidx - tt)) = 1;
                                        end
                                        Post = nanmean(EXP.Bayes.Posterior0{iprobe}(tidx & beforelick,:),1);
                                        %                                     vec1 = interp1(1:100,Post,xParams,'spline');
                                        ParamsPostbeforelick(iprobe,latcorr,spdbins,smthspd,ialpha,idelta,g,Ntau,:) = Post;%vec1;

                                        meanPost = nanmean(EXP.Bayes.Posterior0{iprobe}(tidx & meanbeforelick,:),1);
                                        %                                     vec1 = interp1(1:100,meanPost,xParams,'spline');
                                        ParamsMeanPostbeforelick(iprobe,latcorr,spdbins,smthspd,ialpha,idelta,g,Ntau,:) = meanPost;%vec1;
                                    end
                                end
                            end
                        end
                    end
                    if isempty(ParamsCorr)
                        ParamsCorr = zeros(nProbe,size(latcorrectionlist,2),size(nspdbinslist,2),size(smth_spdlist,2),size(alphalist,2),size(deltalist,2),ngain-1);
                    end
                    %                 for iprobe = 1:nProbe
                    %                     if max(squeeze(ParamsMeanErr0(iprobe,latcorr,spdbins,smthspd,2,:))) > 0
                    %                         vec1 = interp1(1:100,squeeze(ParamsMeanErr0(iprobe,latcorr,spdbins,smthspd,1,:)),[0:0.1:100],'spline');
                    %                         vec2 = interp1(1:100,squeeze(ParamsMeanErr0(iprobe,latcorr,spdbins,smthspd,2,:)),[0:0.1:100],'spline');
                    %                         vec3 = interp1(1:100,squeeze(ParamsMeanErr0(iprobe,latcorr,spdbins,smthspd,3,:)),[0:0.1:100],'spline');
                    %                         [~,idx] = max(xcorr(vec1,vec2,200,'coeff'));
                    %                         ParamsCorr(iprobe,latcorr,spdbins,smthspd,1) = (idx - 201)*0.1;
                    %                         [~,idx] = max(xcorr(vec3,vec2,200,'coeff'));
                    %                         ParamsCorr(iprobe,latcorr,spdbins,smthspd,2) = (idx - 201)*0.1;
                    %                     end
                    %                 end
                    %             end
                    
                    
                end
            end
        end
    end
end
end