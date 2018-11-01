function res = BatchBayesGrandAverage(batch2p)
SetDirs;
if nargin < 1
    strlistvarname = {'2p data','electrophys data'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
    if ok && varnamesel == 1
        batch2p = true;
        
    elseif ok
        batch2p = false;
        
    end
end
if batch2p
    expt = getExperimentList2p;
    datadir = DIRS.data2p;
else
    expt = getExperimentList;
    datadir = DIRS.multichanspikes;
end
suffix = 'win150_speedThresh1_5speedbins_goodonly';%'win150_5speedbins_goodonly';%'win150_5speedbins_VisualSpd_goodonly';%
lambdaSmooth = 2;
corrmaxlag = 180;
nanimal = numel(expt);

for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        if exist([datadir filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'EXP_' suffix '.mat'],'file')
            S = load([datadir filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'EXP_' suffix '.mat']);
            EXP = TVRData;
            EXP.Copyobj(S.EXP);
%             Slat = load([datadir filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'decoderParams_latency.mat']);
            
            res.gainVal{ianimal,iseries} = EXP.SubsetVal.gain;
            res.roomlengthVal{ianimal,iseries} = EXP.SubsetVal.roomlength;
            res.outcomeVal{ianimal,iseries} = EXP.SubsetVal.outcome;
            res.gain{ianimal,iseries} = EXP.data.es.gain;
            res.contrast{ianimal,iseries} = EXP.data.es.contrast;
            res.trialID{ianimal,iseries} = EXP.data.es.trialID;
            res.blanks{ianimal,iseries} = EXP.data.es.blanks;
            res.afterblanks{ianimal,iseries} = EXP.data.es.afterblanks;
            res.trialgainchange{ianimal,iseries} = EXP.data.es.trialgainchange;
            
            res.traj{ianimal,iseries} = EXP.data.es.traj;
            res.trajspeed{ianimal,iseries} = EXP.data.es.trajspeed;
            res.ballspeed{ianimal,iseries} = EXP.data.es.ballspeed;
            
            res.lick{ianimal,iseries} = EXP.data.es.lick;
            res.firstgoodlick{ianimal,iseries} = EXP.data.es.firstgoodlick;
            res.firstrewlick{ianimal,iseries} = EXP.data.es.firstrewlick;
            res.goodlick{ianimal,iseries} = EXP.data.es.goodlick;
            res.preRewlick{ianimal,iseries} = EXP.data.es.preRewlick;
            res.postRewlick{ianimal,iseries} = EXP.data.es.postRewlick;
            res.passivelick{ianimal,iseries} = EXP.data.es.postRewlick;
            res.firstbadlick{ianimal,iseries} = EXP.data.es.postRewlick;
            res.badlick{ianimal,iseries} = EXP.data.es.postRewlick;
            
            res.pupilSize{ianimal,iseries} = EXP.data.es.pupilSize;
            res.eyeYpos{ianimal,iseries} = EXP.data.es.eyeYpos;
            res.eyeXpos{ianimal,iseries} = EXP.data.es.eyeXpos;
            
            runspeed = NaN(size(EXP.data.es.ballspeed));
            runspeed(~isnan(EXP.data.es.ballspeed)) = smthInTime(EXP.data.es.ballspeed(~isnan(EXP.data.es.ballspeed)), EXP.Bayes.sampleRate, EXP.Bayes.smth_spd(1), 'same', [], 'boxcar_centered');
            res.runSpeed{ianimal,iseries} = runspeed;
            visspeed = NaN(size(EXP.data.es.trajspeed));
            visspeed(~isnan(EXP.data.es.trajspeed)) = smthInTime(EXP.data.es.trajspeed(~isnan(EXP.data.es.trajspeed)), EXP.Bayes.sampleRate, EXP.Bayes.smth_spd(1), 'same', [], 'boxcar_centered');            
            res.visSpeed{ianimal,iseries} = visspeed;
            res.Spdbin{ianimal,iseries} = EXP.Bayes.Spdbin;
            
            res.X{ianimal,iseries} = EXP.Bayes.X;
            res.firstgoodlick{ianimal,iseries} = EXP.data.es.firstgoodlick;
            res.contrastVal{ianimal,iseries} = EXP.SubsetVal.contrast;            
            
            
            for iprobe = 1:numel(EXP.Bayes.Posterior0)
                res.DecCells{ianimal,iseries,iprobe} = EXP.Bayes.DecCellidx{iprobe};
                amp_th = 1;maxtol = 1;
                res.Xpred_ave{ianimal,iseries,iprobe} = getCircularAverage(EXP.Bayes.Posterior0{iprobe}',amp_th,maxtol);
                amp_th = 1;maxtol = 0.1;
                res.Xpred_max{ianimal,iseries,iprobe} = getCircularAverage(EXP.Bayes.Posterior0{iprobe}',amp_th,maxtol);
                
                for c = 1:numel(EXP.SubsetVal.contrast)
                    for g = 1:numel(EXP.SubsetVal.gain)
                        for r = 1:numel(EXP.SubsetVal.roomlength)
                            for o = 1:numel(EXP.SubsetVal.outcome)
                                tidx = EXP.getSubsets(c, g, r, o, EXP.Bayes.speed_th);
                                X = EXP.Bayes.X;
                                Xrange = max(floor(X));
                                Post = EXP.Bayes.Posterior0{iprobe};
                                
                                res.tidx{ianimal,iseries,iprobe,c,g,r,o} = tidx;
                                
                                res.PostXSum{ianimal,iseries,iprobe,c,g,r,o} = zeros(size(Post,2), Xrange);
                                res.XSum{ianimal,iseries,iprobe,c,g,r,o} = zeros(1, Xrange);
                                
                                for i = 1:Xrange
                                    res.PostXSum{ianimal,iseries,iprobe,c,g,r,o}(:,i) = nansum(Post(X == i & tidx,:),1)';
                                    res.XSum{ianimal,iseries,iprobe,c,g,r,o}(i) = sum(X == i & tidx & ~isnan(sum(Post,2)));
                                end
                            end
                        end
                    end
                end
            end
            
            
            
            
%             allcontidx = size(S.EXP.Bayes.Ave.PostX,2);
%             outcomeidx = find(S.EXP.Bayes.Ave.outcomeVal == 2);
%             roomlengthidx = 1;
%             for iprobe = 1:2
%                 if ~isempty(S.EXP.CellInfo.RFXpos)
%                     [~,xposmax] = max(mean(S.EXP.CellInfo.RFXpos(S.EXP.CellInfo.XposZmax>2 & S.EXP.CellInfo.Probe == 2,:),1));
%                     res.XposPop{ianimal,iprobe} = [res.XposPop{ianimal,iprobe} xposmax];
%                 else
%                     res.XposPop{ianimal,iprobe} = [res.XposPop{ianimal,iprobe} NaN];
%                 end
%                 if iprobe <= size(Slat.latEstim.ParamsMeanErr0,1)
%                     imax_med = getCircularAverage(squeeze(S.latEstim.ParamsMeanErr0(iprobe,:,1,1,1,1,2,:))',maxtol);
%                     imax_low = getCircularAverage(squeeze(S.latEstim.ParamsMeanErr0(iprobe,:,1,1,1,1,1,:))',maxtol);
%                     imax_high = getCircularAverage(squeeze(S.latEstim.ParamsMeanErr0(iprobe,:,1,1,1,1,3,:))',maxtol);
%                     
%                     latencylist = Slat.latparams.latcorrectionlist(Slat.latparams.latcorrectionlist > 0);
%                     shift_low = imax_low - imax_med;
%                     shift_high = imax_high - imax_med;
%                     optLat{1} = mean(latencylist(abs(shift_low) == min(abs(shift_low))));
%                     optLat{3} = mean(latencylist(abs(shift_high) == min(abs(shift_high))));
%                     optLat{2} = Slat.latEstim.OptLatency;
%                 else
%                     optLat{1} = NaN;
%                     optLat{3} = NaN;
%                     optLat{2} = NaN;
%                 end
%                 
%                 
%                 for g = 1:3
%                     if iprobe <= size(S.EXP.Bayes.Ave.PostX,1) && g <= size(S.EXP.Bayes.Ave.PostX,3) && ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
%                         n{ianimal,iprobe,g} = n{ianimal,iprobe,g}+1;
%                         res.PostX{ianimal,iprobe,g} = res.PostX{ianimal,iprobe,g} + S.EXP.Bayes.Ave.PostX{iprobe,allcontidx,g,roomlengthidx,outcomeidx};
%                         res.PostErrX{ianimal,iprobe,g} = res.PostErrX{ianimal,iprobe,g} + S.EXP.Bayes.Ave.PostErrX{iprobe,allcontidx,g,roomlengthidx,outcomeidx};
%                         res.DistriX{ianimal,iprobe,g} = res.DistriX{ianimal,iprobe,g} + S.EXP.Bayes.Ave.DistriX{iprobe,allcontidx,g,roomlengthidx,outcomeidx};
%                         res.DistriErrX{ianimal,iprobe,g} = res.DistriErrX{ianimal,iprobe,g} + S.EXP.Bayes.Ave.DistriErrX{iprobe,allcontidx,g,roomlengthidx,outcomeidx};
%                         
%                         [~, imaxcorr] = max(S.EXP.Bayes.Ave.MeanXerrCorr{iprobe,allcontidx,g,roomlengthidx,outcomeidx});
%                         res.MeanXerrCorr{ianimal,iprobe,g} = [res.MeanXerrCorr{ianimal,iprobe,g} imaxcorr];
%                         [~, imaxcorr] = max(S.EXP.Bayes.Ave.MeanDistriErrCorr{iprobe,allcontidx,g,roomlengthidx,outcomeidx});
%                         res.MeanDistriErrCorr{ianimal,iprobe,g} = [res.MeanDistriErrCorr{ianimal,iprobe,g} imaxcorr];
%                         
%                         res.MeanXErrMax{ianimal,iprobe,g} = [res.MeanXErrMax{ianimal,iprobe,g} S.EXP.Bayes.Ave.MeanXErrMax{iprobe,allcontidx,g,roomlengthidx,outcomeidx}];
%                         res.MeanXErrMaxPos{ianimal,iprobe,g} = [res.MeanXErrMaxPos{ianimal,iprobe,g} S.EXP.Bayes.Ave.MeanXErrMaxPos{iprobe,allcontidx,g,roomlengthidx,outcomeidx}];
%                         
%                         res.MeanDistriErrMax{ianimal,iprobe,g} = [res.MeanDistriErrMax{ianimal,iprobe,g} S.EXP.Bayes.Ave.MeanDistriErrMax{iprobe,allcontidx,g,roomlengthidx,outcomeidx}];
%                         res.MeanDistriErrMaxPos{ianimal,iprobe,g} = [res.MeanDistriErrMaxPos{ianimal,iprobe,g} S.EXP.Bayes.Ave.MeanDistriErrMaxPos{iprobe,allcontidx,g,roomlengthidx,outcomeidx}];
%                         
%                         res.Latopt{ianimal,iprobe,g} = [res.Latopt{ianimal,iprobe,g} optLat{g}];
%                     else
%                         res.MeanXerrCorr{ianimal,iprobe,g} = [res.MeanXerrCorr{ianimal,iprobe,g} NaN];
%                         res.MeanDistriErrCorr{ianimal,iprobe,g} = [res.MeanDistriErrCorr{ianimal,iprobe,g} NaN];
%                         
%                         res.MeanXErrMax{ianimal,iprobe,g} = [res.MeanXErrMax{ianimal,iprobe,g} NaN];
%                         res.MeanXErrMaxPos{ianimal,iprobe,g} = [res.MeanXErrMaxPos{ianimal,iprobe,g} NaN];
%                         
%                         res.MeanDistriErrMax{ianimal,iprobe,g} = [res.MeanDistriErrMax{ianimal,iprobe,g} NaN];
%                         res.MeanDistriErrMaxPos{ianimal,iprobe,g} = [res.MeanDistriErrMaxPos{ianimal,iprobe,g} NaN];
%                         
%                         res.Latopt{ianimal,iprobe,g} = [res.Latopt{ianimal,iprobe,g} NaN];
%                     end
%                 end
%             end
%             if (expt(ianimal).goodCA1{iseries} == 1) && (expt(ianimal).goodV1{iseries} == 1)
%                 Prange = size(S.EXP.Bayes.Posterior0{1},2);
%                 ErrCA1 = S.EXP.Bayes.Posterior0{1};
%                 G = smooth1D(repmat(ErrCA1,1,3)',lambdaSmooth)';
%                 ErrCA1 = G(:,Prange+1:2*Prange);
%                 [~, ~, XErrCA1] = getcleandecidx(ErrCA1, S.EXP.Bayes.X, 30, 1, S.EXP.data.es.trialID);
%                 ErrV1 = S.EXP.Bayes.Posterior0{2};
%                 G = smooth1D(repmat(ErrV1,1,3)',lambdaSmooth)';
%                 ErrV1 = G(:,Prange+1:2*Prange);
%                 [~, ~, XErrV1] = getcleandecidx(ErrV1, S.EXP.Bayes.X, 30, 1, S.EXP.data.es.trialID);
%                 speeds = EXP.data.es.smthBallSpd;
%                 for g = 1:3
%                     if g <= size(S.EXP.Bayes.Ave.PostX,3)
%                         ncorr{ianimal,g} = ncorr{ianimal,g}+1;
%                         tidx = EXP.getSubsets(find(EXP.SubsetVal.contrast>0), g, 1, find(EXP.SubsetVal.outcome == 2));
%                         tidx = tidx & ~isnan(XErrCA1) & ~isnan(XErrV1);
%                         idx = find(tidx);
%                         res.CA1V1corr{ianimal,g} = res.CA1V1corr{ianimal,g} + smooth(xcorr(XErrCA1(idx),XErrV1(idx),corrmaxlag,'coeff'),3);
%                         
%                         idx_rand = idx(randperm(numel(idx)));
%                         res.CA1V1corr_rand{ianimal,g} = res.CA1V1corr_rand{ianimal,g} + smooth(xcorr(XErrCA1(idx),XErrV1(idx_rand),corrmaxlag,'coeff'),3);
%                         
%                         idx_randX = zeros(size(idx));
%                         for ix = 1:max(S.EXP.Bayes.X)
%                             idxtemp = idx(S.EXP.Bayes.X(idx) == ix);
%                             idx_randX(S.EXP.Bayes.X(idx) == ix) = idxtemp(randperm(numel(idxtemp)));
%                         end
%                         res.CA1V1corr_randX{ianimal,g} = res.CA1V1corr_randX{ianimal,g} + smooth(xcorr(XErrCA1(idx),XErrV1(idx_randX),corrmaxlag,'coeff'),3);
%                         
%                         idx_randSpeed = zeros(size(idx));
%                         Speedrange = [min(speeds(idx)):(max(speeds(idx))-min(speeds(idx)))/10:max(speeds(idx))];
%                         for ispeed = 1:10
%                             idxtemp = idx(speeds(idx) >= Speedrange(ispeed) & speeds(idx) <= Speedrange(ispeed+1));
%                             idx_randSpeed(speeds(idx) >= Speedrange(ispeed) & speeds(idx) <= Speedrange(ispeed+1)) = idxtemp(randperm(numel(idxtemp)));
%                         end
%                         res.CA1V1corr_randS{ianimal,g} = res.CA1V1corr_randS{ianimal,g} + smooth(xcorr(XErrCA1(idx),XErrV1(idx_randSpeed),corrmaxlag,'coeff'),3);
%                         
%                     end
%                 end
%             end
        end
    end
%     for g = 1:3
%         res.CA1V1corr{ianimal,g} = res.CA1V1corr{ianimal,g}/ncorr{ianimal,g};
%         res.CA1V1corr_rand{ianimal,g} = res.CA1V1corr_rand{ianimal,g}/ncorr{ianimal,g};
%         res.CA1V1corr_randX{ianimal,g} = res.CA1V1corr_randX{ianimal,g}/ncorr{ianimal,g};
%         res.CA1V1corr_randS{ianimal,g} = res.CA1V1corr_randS{ianimal,g}/ncorr{ianimal,g};
%     end
%     for iprobe = 1:2
%         for g = 1:3
%             if n{ianimal,iprobe,g} > 0
%                 res.PostX{ianimal,iprobe,g} = res.PostX{ianimal,iprobe,g}/n{ianimal,iprobe,g};
%                 res.PostErrX{ianimal,iprobe,g} = res.PostErrX{ianimal,iprobe,g}/n{ianimal,iprobe,g};
%                 res.DistriX{ianimal,iprobe,g} = res.DistriX{ianimal,iprobe,g}/n{ianimal,iprobe,g};
%                 res.DistriErrX{ianimal,iprobe,g} = res.DistriErrX{ianimal,iprobe,g}/n{ianimal,iprobe,g};
%                 
%                 res.PostX{ianimal,iprobe,g} = smooth2D(res.PostX{ianimal,iprobe,g},lambdaSmooth);
%                 res.PostErrX{ianimal,iprobe,g} = smooth2D(res.PostErrX{ianimal,iprobe,g},lambdaSmooth);
% %                 res.DistriX{ianimal,iprobe,g} = smooth2D(res.DistriX{ianimal,iprobe,g},lambdaSmooth);
% %                 res.DistriErrX{ianimal,iprobe,g} = smooth2D(res.DistriErrX{ianimal,iprobe,g},lambdaSmooth);
%             end
%         end
%     end
end

% for iprobe = 1:2
%     for g = 1:3
%         for ianimal = 1:nanimal
%             if sum(res.PostX{ianimal,iprobe,g}(:))>0
%                 nAll{iprobe,g} = nAll{iprobe,g} + 1;
%                 res.PostXAll{iprobe,g} = res.PostXAll{iprobe,g} + res.PostX{ianimal,iprobe,g};
%                 res.PostErrXAll{iprobe,g} = res.PostErrXAll{iprobe,g} + res.PostErrX{ianimal,iprobe,g};
%                 res.DistriXAll{iprobe,g} = res.DistriXAll{iprobe,g} + res.DistriX{ianimal,iprobe,g};
%                 res.DistriErrXAll{iprobe,g} = res.DistriErrXAll{iprobe,g} + res.DistriErrX{ianimal,iprobe,g};
%             end
%         end
%         res.PostXAll{iprobe,g} = res.PostXAll{iprobe,g}/nAll{iprobe,g};
%         res.PostErrXAll{iprobe,g} = res.PostErrXAll{iprobe,g}/nAll{iprobe,g};
%         res.DistriXAll{iprobe,g} = res.DistriXAll{iprobe,g}/nAll{iprobe,g};
%         res.DistriErrXAll{iprobe,g} = res.DistriErrXAll{iprobe,g}/nAll{iprobe,g};
%     end
% end
% for g = 1:3
%     for ianimal = 1:nanimal
%         if sum(res.CA1V1corr{ianimal,g})>0
%             ncorrAll{g} = ncorrAll{g} + 1;
%             res.CA1V1corrAll{g} = res.CA1V1corrAll{g} + res.CA1V1corr{ianimal,g};
%             res.CA1V1corrAll_rand{g} = res.CA1V1corrAll_rand{g} + res.CA1V1corr_rand{ianimal,g};
%             res.CA1V1corrAll_randX{g} = res.CA1V1corrAll_randX{g} + res.CA1V1corr_randX{ianimal,g};
%             res.CA1V1corrAll_randS{g} = res.CA1V1corrAll_randS{g} + res.CA1V1corr_randS{ianimal,g};
%         end
%     end
%     res.CA1V1corrAll{g} = res.CA1V1corrAll{g}/ncorrAll{g};
%     res.CA1V1corrAll_rand{g} = res.CA1V1corrAll_rand{g}/ncorrAll{g};
%     res.CA1V1corrAll_randX{g} = res.CA1V1corrAll_randX{g}/ncorrAll{g};
%     res.CA1V1corrAll_randS{g} = res.CA1V1corrAll_randS{g}/ncorrAll{g};
% end
end

function mat_out = smooth2D(mat,lambdaSmooth)
mat(isnan(mat)) = 0;
G = smooth1D(repmat(mat,3,3),lambdaSmooth);
H = smooth1D(G',lambdaSmooth)';
mat_out = H(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
end