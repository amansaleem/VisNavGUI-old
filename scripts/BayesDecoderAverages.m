function EXP = BayesDecoderAverages(EXP,maxtol)
es = EXP.data.es;
EXP.Bayes.Ave = [];

cbase = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
gbase = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
rbase = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
obase = find(EXP.SubsetVal.outcome == 2);

clean_th = 0.3*size(EXP.Bayes.Posterior0{1},2);
lambdaSmooth = 2;

EXP.Bayes.X = EXP.Bayes.X0;

for iprobe = 1:numel(EXP.Bayes.Posterior0)
    Prange = size(EXP.Bayes.Posterior0{iprobe},2);
    Xrange = max(EXP.Bayes.X);
    for i = 1:Xrange
        EXP.Bayes.PosError0{iprobe}(round(EXP.Bayes.X) == i,:) = circshift(EXP.Bayes.Posterior0{iprobe}(round(EXP.Bayes.X) == i,:),floor(Prange/2)-i,2);
    end
end

EXP.Bayes.Ave.contrastVal = EXP.SubsetVal.contrast;
EXP.Bayes.Ave.gainVal = EXP.SubsetVal.gain;
EXP.Bayes.Ave.roomlengthVal = EXP.SubsetVal.roomlength;
EXP.Bayes.Ave.outcomeVal = EXP.SubsetVal.outcome;

for iprobe = 1:numel(EXP.Bayes.Posterior0)
    for c = 1:numel(EXP.SubsetVal.contrast) +1
        for g = 1:numel(EXP.SubsetVal.gain)
            for r = 1:numel(EXP.SubsetVal.roomlength)
                for o = 1:numel(EXP.SubsetVal.outcome)
                    if c > numel(EXP.SubsetVal.contrast)
                        contidx = find(EXP.SubsetVal.contrast>0);
                    else
                        contidx = c;
                    end
                    tidx = EXP.getSubsets(contidx, g, r, o, EXP.Bayes.speed_th);
                    X = EXP.Bayes.X;
                    Xrange = max(floor(X));
                    Prange = size(EXP.Bayes.Posterior0{iprobe},2);
                    Post = EXP.Bayes.Posterior0{iprobe};
                    Err = EXP.Bayes.PosError0{iprobe};
                    
                    EXP.Bayes.Ave.PostX{iprobe,c,g,r,o} = zeros(size(Post,2), Xrange);
                    EXP.Bayes.Ave.PostErrX{iprobe,c,g,r,o} = zeros(size(Post,2), Xrange);
                    EXP.Bayes.Ave.DistriX{iprobe,c,g,r,o} = zeros(size(Post,2), Xrange);
                    EXP.Bayes.Ave.DistriErrX{iprobe,c,g,r,o} = zeros(size(Post,2), Xrange);
                    
                    for i = 1:Xrange
                        EXP.Bayes.Ave.PostX{iprobe,c,g,r,o}(:,i) = nanmean(Post(X == i  & tidx,:),1)';
                        EXP.Bayes.Ave.PostErrX{iprobe,c,g,r,o}(:,i) = nanmean(Err(X == i  & tidx,:),1)';
                    end
                    
                    [~, Xpred,Err] = getcleandecidx(Post, X, clean_th, maxtol, es.trialID);
                    Err = Err + 50;
                    Xpred = floor(mod(Xpred,Xrange))+1;
                    [F, ctrs1, ctrs2, H] = smoothhist2D_corrected([X(tidx & ~isnan(Xpred)) Xpred(tidx & ~isnan(Xpred))], lambdaSmooth, [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                    EXP.Bayes.Ave.DistriX{iprobe,c,g,r,o} = F;
                    [F, ctrs1, ctrs2, H] = smoothhist2D_corrected([X(tidx & ~isnan(Err)) Err(tidx & ~isnan(Err))], lambdaSmooth, [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                    EXP.Bayes.Ave.DistriErrX{iprobe,c,g,r,o} = F;
                    
                    EXP.Bayes.Ave.PostXave{iprobe,c,g,r,o} = getCircularAverage(EXP.Bayes.Ave.PostX{iprobe,c,g,r,o},1,maxtol);
                    EXP.Bayes.Ave.ErrXave{iprobe,c,g,r,o} = getCircularAverage(EXP.Bayes.Ave.PostErrX{iprobe,c,g,r,o},1,maxtol);
                    EXP.Bayes.Ave.DistriXave{iprobe,c,g,r,o} = getCircularAverage(EXP.Bayes.Ave.DistriX{iprobe,c,g,r,o},0,maxtol);
                    EXP.Bayes.Ave.DistriErrXave{iprobe,c,g,r,o} = getCircularAverage(EXP.Bayes.Ave.DistriErrX{iprobe,c,g,r,o},0,maxtol);
                    
                    EXP.Bayes.Ave.MeanXerr{iprobe,c,g,r,o} = mean(EXP.Bayes.Ave.PostErrX{iprobe,c,g,r,o},2);
                    EXP.Bayes.Ave.MeanDistriErr{iprobe,c,g,r,o} = mean(EXP.Bayes.Ave.DistriX{iprobe,c,g,r,o},2);
                end
            end
        end
    end
    for c = 1:numel(EXP.SubsetVal.contrast) +1
        for g = 1:numel(EXP.SubsetVal.gain)
            for r = 1:numel(EXP.SubsetVal.roomlength)
                for o = 1:numel(EXP.SubsetVal.outcome)
                    EXP.Bayes.Ave.MeanXerrCorr{iprobe,c,g,r,o} = xcorr(EXP.Bayes.Ave.MeanXerr{iprobe,c,g,r,o}, EXP.Bayes.Ave.MeanXerr{iprobe,cbase,gbase,rbase,obase},50,'coeff');
                    EXP.Bayes.Ave.MeanDistriErrCorr{iprobe,c,g,r,o} = xcorr(EXP.Bayes.Ave.MeanDistriErr{iprobe,c,g,r,o}, EXP.Bayes.Ave.MeanDistriErr{iprobe,cbase,gbase,rbase,obase},50,'coeff');
                    if sum(isnan(EXP.Bayes.Ave.MeanXerr{iprobe,c,g,r,o})) < numel(EXP.Bayes.Ave.MeanXerr{iprobe,c,g,r,o}) - 2
                        x_interp = 0.1:0.1:100;
                        x_init = 1:numel(EXP.Bayes.Ave.MeanXerr{iprobe,c,g,r,o});
                        meanErr_interp = interp1(x_init,EXP.Bayes.Ave.MeanXerr{iprobe,c,g,r,o},x_interp,'spline');
                        [vmax, imax] = max(meanErr_interp);
                        EXP.Bayes.Ave.MeanXErrMax{iprobe,c,g,r,o} = vmax;
                        EXP.Bayes.Ave.MeanXErrMaxPos{iprobe,c,g,r,o} = x_interp(imax);

                        x_interp = 0:0.1:100;
                        x_init = 1:numel(EXP.Bayes.Ave.MeanDistriErr{iprobe,c,g,r,o});
                        meanErr_interp = interp1(x_init,EXP.Bayes.Ave.MeanDistriErr{iprobe,c,g,r,o},x_interp,'spline');
                        [vmax, imax] = max(meanErr_interp);
                        EXP.Bayes.Ave.MeanDistriErrMax{iprobe,c,g,r,o} = vmax;
                        EXP.Bayes.Ave.MeanDistriErrMaxPos{iprobe,c,g,r,o} = x_interp(imax);
                    end
                end
            end
        end
    end
end
end