function EXP = ComputeBayesAverage(EXP,Nperm)
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

maxtol = 0.1;%0.1;
EXP.Bayes.Ave.contrastVal = EXP.SubsetVal.contrast;
EXP.Bayes.Ave.gainVal = EXP.SubsetVal.gain;
EXP.Bayes.Ave.roomlengthVal = EXP.SubsetVal.roomlength;
EXP.Bayes.Ave.outcomeVal = EXP.SubsetVal.outcome;
EXP.Bayes.Ave.Nperm = Nperm;


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
                    [F, ctrs1, ctrs2, H] = smoothhist2D_corrected([X(tidx & ~isnan(Xpred)) Xpred(tidx & ~isnan(Xpred))], lambdaSmooth, [Xrange 100], 1:Xrange, 1:100, true, true);
                    EXP.Bayes.Ave.DistriX{iprobe,c,g,r,o} = F;
                    [F, ctrs1, ctrs2, H] = smoothhist2D_corrected([X(tidx & ~isnan(Err)) Err(tidx & ~isnan(Err))], lambdaSmooth, [Xrange 100], 1:Xrange, 1:100, true, true);
                    EXP.Bayes.Ave.DistriErrX{iprobe,c,g,r,o} = F;
                    
                    allidx = find(tidx);
                    postpredaveMean = 0;
                    postpredaveStd = 0;
                    distripredaveMean = 0;
                    distripredaveStd = 0;
                    errpredaveMean = 0;
                    errpredaveStd = 0;
                    distrierrpredaveMean = 0;
                    distrierrpredaveStd = 0;
                    errXmeanMean = 0;
                    errXmeanStd = 0;
                    errdistrimeanMean = 0;
                    errdistrimeanStd = 0;
                    for iperm = 0:Nperm
                        outidx = false(size(tidx));
                        if iperm > 0
                            idx = allidx(1+(iperm-1)*floor(numel(allidx)/Nperm):(iperm)*floor(numel(allidx)/Nperm));
                            outidx(idx) = true;
                            posttemp = zeros(size(EXP.Bayes.Ave.PostX{iprobe,c,g,r,o}));
                            errtemp = zeros(size(EXP.Bayes.Ave.PostX{iprobe,c,g,r,o}));
                            for i = 1:Xrange
                                posttemp(:,i) = nanmean(Post(X == i & tidx & ~outidx ,:),1)';
                                errtemp(:,i) = nanmean(Err(X == i & tidx & ~outidx ,:),1)';
                            end
                            [~, Xpred,Err] = getcleandecidx(Post, X, clean_th, maxtol, es.trialID);
                            Xpred = floor(mod(Xpred,Xrange))+1;
                            [F, ctrs1, ctrs2, H] = smoothhist2D_corrected([X(tidx & ~isnan(Xpred) & ~outidx) Xpred(tidx & ~isnan(Xpred) & ~outidx)], lambdaSmooth, [Xrange 100], 1:Xrange, 1:100, true, true);
                            distritemp = F;
                            
                            [F, ctrs1, ctrs2, H] = smoothhist2D_corrected([X(tidx & ~isnan(Xpred) & ~outidx) Xpred(tidx & ~isnan(Xpred) & ~outidx)], lambdaSmooth, [Xrange 100], 1:Xrange, 1:100, true, true);
                            distrierrtemp = F;
                        else
                            posttemp = EXP.Bayes.Ave.PostX{iprobe,c,g,r,o};
                            errtemp = EXP.Bayes.Ave.PostErrX{iprobe,c,g,r,o};
                            distritemp = EXP.Bayes.Ave.DistriX{iprobe,c,g,r,o};
                            distrierrtemp = EXP.Bayes.Ave.DistriErrX{iprobe,c,g,r,o};
                        end
                        
                        postpredave = getCircularAverage(posttemp,1,maxtol);
                        errpredave = getCircularAverage(errtemp,1,maxtol);
                        distripredave = getCircularAverage(distritemp,0,maxtol);
                        distrierrpredave = getCircularAverage(distrierrtemp,0,maxtol);
                        
                        errXmean = mean(errtemp,2);
                        errdistrimean = mean(distritemp,2);
                        
                        if iperm == 0
                            EXP.Bayes.Ave.PostXave{iprobe,c,g,r,o} = postpredave;
                            EXP.Bayes.Ave.ErrXave{iprobe,c,g,r,o} = errpredave;
                            EXP.Bayes.Ave.MeanXerr{iprobe,c,g,r,o} = errXmean;
                            EXP.Bayes.Ave.DistriXave{iprobe,c,g,r,o} = distripredave;
                            EXP.Bayes.Ave.DistriErrXave{iprobe,c,g,r,o} = distrierrpredave;
                            EXP.Bayes.Ave.MeanDistriErr{iprobe,c,g,r,o} = errdistrimean;
                        end
                                                
                        postpredaveMean = postpredaveMean + postpredave/Nperm;
                        postpredaveStd = postpredaveStd + (postpredave - EXP.Bayes.Ave.PostXave{iprobe,c,g,r,o}).^2*(Nperm-1)/Nperm;
                        errpredaveMean = errpredaveMean + errpredave/Nperm;
                        errpredaveStd = errpredaveStd + (errpredave - EXP.Bayes.Ave.ErrXave{iprobe,c,g,r,o}).^2*(Nperm-1)/Nperm;
                        distripredaveMean = distripredaveMean + distripredave/Nperm;
                        distripredaveStd = distripredaveStd + (distripredave - EXP.Bayes.Ave.DistriErrXave{iprobe,c,g,r,o}).^2*(Nperm-1)/Nperm;
                        distrierrpredaveMean = distrierrpredaveMean + distrierrpredave/Nperm;
                        distrierrpredaveStd = distrierrpredaveStd + (distrierrpredave - EXP.Bayes.Ave.DistriErrXave{iprobe,c,g,r,o}).^2*(Nperm-1)/Nperm;
                        
                        errXmeanMean = errXmeanMean + errXmean/Nperm;
                        errXmeanStd = errXmeanStd + (errXmean - EXP.Bayes.Ave.MeanXerr{iprobe,c,g,r,o}).^2*(Nperm-1)/Nperm;
                        errdistrimeanMean = errdistrimeanMean + errdistrimean/Nperm;
                        errdistrimeanStd = errdistrimeanStd + (errdistrimean - EXP.Bayes.Ave.MeanDistriErr{iprobe,c,g,r,o}).^2*(Nperm-1)/Nperm;
                    end
                    EXP.Bayes.Ave.PostXstd{iprobe,c,g,r,o} = (postpredaveStd/Nperm).^0.5;
                    EXP.Bayes.Ave.ErrXstd{iprobe,c,g,r,o} = (errXmeanStd/Nperm).^0.5;
                    EXP.Bayes.Ave.MeanErrXstd{iprobe,c,g,r,o} = (errpredaveStd/Nperm).^0.5;
                    EXP.Bayes.Ave.DistriXstd{iprobe,c,g,r,o} = (distripredaveStd/Nperm).^0.5;
                    EXP.Bayes.Ave.DistriErrXstd{iprobe,c,g,r,o} = (distrierrpredaveStd/Nperm).^0.5;
                    EXP.Bayes.Ave.DistriMeanErrXstd{iprobe,c,g,r,o} = (errdistrimeanStd/Nperm).^0.5;
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

function predave = getCircularAverageXdec(mat,maxtol)
    Prange = size(mat,1);
    Xrange = size(mat,2);
    postfilt = zeros(Prange,1);
    postfilt(1:Prange) = (0:(Prange-1))'/(Prange)*(2*pi)-pi/2;

    [maxval, ~] = max(mat,[],1);
    predmax = zeros(size(maxval));
    for i = 1:Xrange
        mat(mat(:,i) < (1-maxtol)*maxval(i),i) = 0;
    end

    a = sum((mat'*cos(postfilt)),2);b = sum((mat'*sin(postfilt)),2);
    predave = atan2(-a,b)+pi;%./sum(exp(log(2)*Post(tidx,:)),2);
    predave = predave/(2*pi)*Prange+1;
    predplus = predave + Prange;
    predminus = predave - Prange;
    xx = (1:Xrange)';
    predave = predave.* (abs(xx-predave)<abs(xx-predplus) & abs(xx-predave)<abs(xx-predminus)) + predplus.* (abs(xx-predplus)<abs(xx-predave) & abs(xx-predplus)<abs(xx-predminus)) + predminus.* (abs(xx-predminus)<abs(xx-predave) & abs(xx-predminus)<abs(xx-predplus));
    predave = smooth(predave,3);
end