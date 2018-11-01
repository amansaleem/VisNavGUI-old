function PlotVar = UpdateplotDecodSpks(PlotVar,EXP,Layout)
PlotVar.Plots = [];
es = EXP.data.es;
clean_th = 0.3*size(EXP.Bayes.Posterior0{1},2);
maxtol = PlotVar.maxtolerance;
interpol = 1;%0.05;
ampth = 0;%0;
lambdaSmooth = 2;
page = 2;
win = 1;

Posttype = 1;
Xtype = 1;

[goodidx, Xpred, Err] = getcleandecidx(EXP.Bayes.Posterior0{1}, EXP.Bayes.Xsmth0{Xtype},EXP.Bayes.goodidx_th,maxtol,EXP.data.es.trialID);
theta = LFPsignals((Err)', 60, 6, 9);
theta.phase = mod(unwrap(theta.phase)/(2*pi)*360,360);

% for iprobe = 1:size(EXP.Bayes.Posterior0,1)
%     [~,EXP.Bayes.MaxDecodedPosition0{iprobe,1}] = max(EXP.Bayes.Posterior0{iprobe,1},[],2);
% end

% Prange = size(EXP.Bayes.Posterior0{1,1},2);
% Postmin = 1.5;
% for iprobe = 1:size(EXP.Bayes.Posterior0,1)
%     EXP.Bayes.MaxDecodedPosition0{iprobe,1} = NaN(size(EXP.Bayes.Posterior0{iprobe,1},1),1);
%     Post = EXP.Bayes.Posterior0{iprobe,1};
%     Post(Post<=Postmin) = 0;
%     for tt = 1:size(EXP.Bayes.Posterior0{iprobe,1},1)
%         [pks,loc] = findpeaks(repmat(Post(tt,:),[1 3]));
%         loc = loc - Prange;
%         loc(pks<=Postmin | loc < 1 | loc > Prange) = [];
%         if ~isempty(loc)
%             [~,locminidx] = min(abs(Prange/(2*pi)*circ_dist(2*pi/Prange*loc,2*pi/Prange*EXP.Bayes.Xsmth0{Xtype}(tt))));
% %             EXP.Bayes.MaxDecodedPosition0{iprobe,1}(tt) = loc(locminidx);
%             [outfieldidx,fieldix,amp_th_i] = findfield(EXP.Bayes.Posterior0{iprobe,1}(tt,:),Postmin,loc(locminidx),false);
%             post_tt = EXP.Bayes.Posterior0{iprobe,1}(tt,:);
%             post_tt(outfieldidx) = 1;
%             EXP.Bayes.MaxDecodedPosition0{iprobe,1}(tt) = getCircularAverage(post_tt',0,1);%loc(locminidx);
%         else
%             EXP.Bayes.MaxDecodedPosition0{iprobe,1}(tt) = NaN;
%         end
%     end
% end

if ~PlotVar.Foverlap
    Layout.DivideWindow(page, win, numel(PlotVar.ChosenProbe)*size(PlotVar.ChosenObj,1)*size(PlotVar.ChosenContrast,2), size(PlotVar.ChosenGain,2)*size(PlotVar.ChosenRoomlength,2)*size(PlotVar.ChosenOutcome,2));
else
    Layout.DivideWindow(page, win, numel(PlotVar.ChosenProbe)*size(PlotVar.ChosenObj,1), 1);
end
colorcond = cell(numel(EXP.SubsetVal.gain),numel(EXP.SubsetVal.contrast));

for g = 1:numel(EXP.SubsetVal.gain)
    for c = 1:numel(EXP.SubsetVal.contrast)
        rgb = zeros(1,3);
        rgb(min(g,3)) = 1 * c/numel(EXP.SubsetVal.contrast);
        rgb(2) = rgb(2)+0.5;
        colorcond{g,c} = min(rgb,[1 1 1]);
    end
end
cidx = find(ismember(EXP.SubsetVal.contrast,[0.2 0.5 0.7 0.8]));
for c = 1:numel(cidx)
    if sum(EXP.SubsetVal.gain == 0.4)>0
        colorcond{EXP.SubsetVal.gain == 0.4,cidx(c)} = 'c';
    end
    if sum(EXP.SubsetVal.gain == 0.5)>0
        colorcond{EXP.SubsetVal.gain == 0.5,cidx(c)} = [0.5 0.5 0.5];
    end
    if sum(EXP.SubsetVal.gain == 0.6)>0
        colorcond{EXP.SubsetVal.gain == 0.6,cidx(c)} = 'm';
    end
end
lwidth = 1;
nplot = 0;
for iprobe = 1:numel(PlotVar.ChosenProbe)
    selprobe = PlotVar.ChosenProbe(iprobe);
    probestr = {'CA1','V1'};
    for c = 1:size(PlotVar.ChosenContrast,2)
        for g = 1:size(PlotVar.ChosenGain,2)
            for r = 1:size(PlotVar.ChosenRoomlength,2)
                for o = 1:size(PlotVar.ChosenOutcome,2)
                    for k = 1:size(PlotVar.ChosenObj,1)
                        tidx = EXP.getSubsets(PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o),EXP.Bayes.speed_th);
                        tidxref = EXP.getSubsets(EXP.Bayes.TrainContrast,EXP.Bayes.TrainGain,EXP.Bayes.TrainRoomlength,EXP.Bayes.TrainOutcome,EXP.Bayes.speed_th);
                        Ntrials = numel(unique((es.trialID(tidx))));
                        wintitle = [probestr{iprobe} ' ' 'Contrast: ' num2str(EXP.SubsetVal.contrast(PlotVar.ChosenContrast(:,c))) ' Gain: ' num2str(EXP.SubsetVal.gain(PlotVar.ChosenGain(:,g))) ' (' num2str(Ntrials) ' trials)'];
                        if ~PlotVar.Foverlap
                            Layout.subwindow{page, win, (iprobe-1)*size(PlotVar.ChosenContrast,2) + (c-1)*size(PlotVar.ChosenObj,2) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}.Title = wintitle;
                        else
                            Layout.subwindow{page, win, (iprobe-1) + k, 1}.Title = wintitle;
                        end
                        
                        exptidx = 0;
                        seriesid = unique(es.series);
                        for sid = 1:numel(seriesid)
                            exptidx = exptidx + (es.series == seriesid(sid) & ismember(es.iexp,PlotVar.explist(sid,:)));
                        end
                        
                        tidx = tidx & exptidx;
%                         tidx = tidx & ismember(es.series,PlotVar.Series) & ismember(es.trialgainchange,PlotVar.nthetaphsbins);
                        tidxref = tidxref & exptidx;
                        
                        phsrefidx = floor(PlotVar.thetaDecphase/(360/EXP.Bayes.nthetaphsbins)) + 1;                                                
                        if PlotVar.Fthetabins
                            nthetaphase = PlotVar.nthetaphsbins;%EXP.Bayes.nthetaphsbins;
                            phsval = mod(EXP.Bayes.LFPphase{1},360);%mod(es.LFPphase{selprobe}(:,PlotVar.thetaChannel)-180,360);
                            phsidx = ismember(phsval,mod(PlotVar.thetaphase:(PlotVar.thetaphase + 360/nthetaphase),360));
                        else
                            if ~PlotVar.FthetaPost
                                phsidx = true(size(tidx));%spdidx;%
                            else
                                phsval = mod(EXP.Bayes.LFPphase{1},360);%mod(es.LFPphase{selprobe}(:,PlotVar.thetaChannel)-180,360);
                                phs0 = (phsrefidx-1)*360/EXP.Bayes.nthetaphsbins;
                                phsidx = ismember(phsval,mod(phs0:(phs0 + 360/EXP.Bayes.nthetaphsbins),360));
                            end
                        end
                        
                        if PlotVar.Fspdbins
                            spdidx = ismember(EXP.Bayes.Spdbin,PlotVar.spdbin:PlotVar.spdbin);%ismember(EXP.Bayes.Eyebin,PlotVar.spdbin:PlotVar.spdbin);%ismember(es.trialgainchange,PlotVar.spdbin);%
                        else
                            spdidx = true(size(tidx));%spdidx;%
                        end
                        
                        if strcmp(PlotVar.ChosenObj{k}, 'Posterior x X') || strcmp(PlotVar.ChosenObj{k}, 'Posterior x D')
                            if ~PlotVar.Foverlap
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1)*size(PlotVar.ChosenContrast,2) + (c-1)*size(PlotVar.ChosenObj,2) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            else
                                nplot = (iprobe-1) + k;
                                if numel(PlotVar.Plots) < nplot
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1) + k, 1});
                                end
                            end
                            
                            if ~PlotVar.FthetaPost
                                Post = EXP.Bayes.Posterior0{selprobe,Posttype};%EXP.Bayes.Posterior0{selprobe};%evalin('base','Posterior0Ave');%EXP.Bayes.Posterior0{selprobe};
                            else
                                Post = EXP.Bayes.Posterior{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                            end
                            
                            Prange = size(Post,2);
                            if strcmp(PlotVar.ChosenObj{k}, 'Posterior x X')
                                X = EXP.Bayes.Xsmth0{Xtype};%EXP.Bayes.predicted{Xtype};%
                                gainratio = unique(es.gain(tidx)/mode(es.gain));
                                tex2 = 0.416 * Prange;
                                tex3 = 0.584 * Prange;
                            elseif strcmp(PlotVar.ChosenObj{k}, 'Posterior x D')
                               X = EXP.Bayes.D;
                                gainratio = unique(es.gain(tidx)/mode(es.gain));
                                tex2 = 0.416 / unique(es.gain(tidx)/mode(es.gain)) * Prange;
                                tex3 = 0.584 / unique(es.gain(tidx)/mode(es.gain)) * Prange;
                            end
                            Xrange = max(floor(X));
                            
                            
                            
                            if PlotVar.Fgoodtimebins
                                [goodidx, Xpred, Err] = getcleandecidx(Post, X, 2, maxtol, es.trialID);%getcleandecidx(Post, X, clean_th, maxtol, es.trialID);
                            else
                                goodidx = true(size(tidx));
                            end
                            goodidx = ~isnan(EXP.Bayes.MeanDecodedPosition0{selprobe,1}) & ~isnan(EXP.Bayes.MaxDecodedPosition0{selprobe,1});
%                             goodidx = goodidx & sum(EXP.Bayes.PosError0{selprobe,1}(:,floor(Prange/4):floor(3*Prange/4)),2)>floor(Prange/2);%abs(Prange/(2*pi)*circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{selprobe,1},2*pi/Prange*X))<=25;%
                            
                            
                            mat = zeros(size(Post,2), Xrange);
                            matdist = zeros(size(Post,2), Xrange);
                            matref = zeros(size(Post,2), Xrange);
                            matrefdist = zeros(size(Post,2), Xrange);
                            
                            if ~PlotVar.FdecXdistribution
                                for i = 1:Xrange
                                    mat(:,i) = nanmean(Post(X == i  & tidx & spdidx & phsidx & goodidx,:),1)';%nanmean(Post(X == i  & tidx & spdidx & phsidx & goodidx,:),1)';
%                                     maxdecErr(i) = Prange/(2*pi)*sqrt(circ_mean(circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{selprobe,1}(X == i  & tidx & spdidx & phsidx & goodidx),2*pi/Prange*X(X == i  & tidx & spdidx & phsidx & goodidx)).^2));
%                                     meandecErr(i) = Prange/(2*pi)*sqrt(circ_mean(circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{selprobe,1}(X == i  & tidx & spdidx & phsidx & goodidx),2*pi/Prange*X(X == i  & tidx & spdidx & phsidx & goodidx)).^2));
                                end
%                                 maxdecErr = Prange/(2*pi)*sqrt(circ_mean(circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{selprobe,1}(tidx & spdidx & phsidx & goodidx),2*pi/Prange*X(tidx & spdidx & phsidx & goodidx)).^2));
%                                 meandecErr = Prange/(2*pi)*sqrt(circ_mean(circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{selprobe,1}(tidx & spdidx & phsidx & goodidx),2*pi/Prange*X(tidx & spdidx & phsidx & goodidx)).^2));
%                                 disp(['max decoding error (electrode#' num2str(selprobe) '): ' num2str(nanmean(maxdecErr)) ' cm'])
%                                 disp(['mean decoding error(electrode#' num2str(selprobe) '): ' num2str(nanmean(meandecErr)) ' cm'])
                                
%                                 trialsidx = unique(es.trialID(tidx & spdidx & phsidx & goodidx));
%                                 matSum = NaN(size(Post,2), Xrange, numel(trialsidx));
%                                 XSum = zeros(size(Post,2), Xrange, numel(trialsidx));
%                                 for itrial = 1:numel(trialsidx)
%                                 for i = 1:Xrange
%                                     if sum(X == i  & tidx & spdidx & phsidx & goodidx & es.trialID == trialsidx(itrial))>0
% %                                     mat(:,i) = nanmean(Post(X == i  & tidx & spdidx & phsidx & goodidx,:),1)';
%                                     matSum(:,i,itrial) = nansum(Post(X == i  & tidx & spdidx & phsidx & goodidx & es.trialID == trialsidx(itrial),:),1)';
%                                     XSum(:,i,itrial) = nansum(~isnan(Post(X == i  & tidx & spdidx & phsidx & goodidx  & es.trialID == trialsidx(itrial),:)),1)';
% %                                     matref(:,i) = nanmean(Post(X == i  & tidxref & spdidx & phsidx & goodidx,:),1)';
%                                     end
%                                 end
%                                 end
%                                 mat = nanmean(matSum./XSum,3);
%                                 
%                                 for j = 1:Prange
%                                     mat(j,:) = mat(j,:) + matSum(j,:)./XSum(j,:);%matSum(j,:)./special_smooth_1d(XSum(j,:), lambdaSmooth/Xrange, [], Xrange);%special_smooth_1d(matSum(j,:), lambdaSmooth/Xrange, [], Xrange)./special_smooth_1d(XSum(j,:), lambdaSmooth/Xrange, [], Xrange);%
%                                 end
%                                 end
                            else
%                                 Posttemp = Post;
%                                 Xpred = EXP.Bayes.MaxDecodedPosition0{iprobe,1};%EXP.Bayes.DecodingError0{iprobe,1}+50;%EXP.Bayes.DecodingBias0{iprobe,1}+50;%getCircularAverage(Posttemp',0,0.01);
%                                 [F, ctrs1, ctrs2, H] = smoothhist2D_corrected([X(tidx & spdidx & phsidx & goodidx & ~isnan(Xpred)) Xpred(tidx & spdidx & phsidx & goodidx & ~isnan(Xpred))], [lambdaSmooth lambdaSmooth], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
                                
                                XSum = zeros(1,Xrange);
                                for i = 1:Xrange
                                    XSum(i) = nansum(X(tidx & spdidx & phsidx & goodidx & ~isnan(Xpred))==i);
                                end
%                                 Postcorr = NaN(size(Post));
%                                 Postref = NaN(size(EXP.Bayes.Posterior0{iprobe,1}));
%                                 for tt = 1:size(EXP.Bayes.Posterior0{iprobe,1},1)
%                                     if ~isnan(EXP.Bayes.MaxDecodedPosition0{iprobe,1}(tt))
%                                         Postref(tt,:) = circshift(EXP.Bayes.Posterior0{iprobe,1}(tt,:),50-round(EXP.Bayes.MaxDecodedPosition0{iprobe,1}(tt)),2);
%                                     end
%                                 end
%                                 Postref = nanmean(Postref,1);%nanmean(EXP.Bayes.PosError0{iprobe,1},1);
%                                 Postref = Postref-nanmean(Postref);
%                                 Post = Post - repmat(nanmean(Post,2),[1 Prange]);
%                                 ishift = 0;
%                                 for xshift = -49:50
%                                     ishift = ishift+1;
%                                     Postcorr(:,ishift) = sum(Post(:,:).*repmat(circshift(Postref,xshift),[size(Post,1) 1]),2);
%                                 end
%                                 Xpred = getCircularAverage(Postcorr',0,0.1,0.05);

%                                 Xpred = NaN(size(Post,1),1);
%                                 Xpred(~isnan(sum(Post,2))) = getCircularAverage(Post(~isnan(sum(Post,2)),:)',1,1,1);
%                                 Xpred = Prange/2+Prange/(2*pi)*circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{selprobe,1},2*pi/Prange*X);%
%                                 [~,Xpred] = max(EXP.Bayes.Posterior0{selprobe,1},[],2);
                                Xpred = EXP.Bayes.MaxDecodedPosition0{selprobe,1};
%                                 Xpred = Prange/2+Prange/(2*pi)*circ_dist(2*pi/Prange*Xpred,2*pi/Prange*X);%
%                                 Xpred = EXP.Bayes.MaxDecodedPosition0{selprobe,1};%(EXP.Bayes.DecodingError0{selprobe,1}+floor(Prange/2));%EXP.Bayes.DecodingBias0{iprobe,1}+floor(Prange/2);%
                                [F, ~, ~, ~] = smoothhist2D_AS([X(tidx & spdidx & phsidx & goodidx & ~isnan(Xpred)) Xpred(tidx & spdidx & phsidx & goodidx & ~isnan(Xpred))],...
                                        [NaN NaN], [Xrange Prange], 1:Xrange, 1:Prange, true, true);
%                                 F = F./repmat(XSum/(1/100),[Prange 1]);
                                
                                F = F./repmat(XSum*(1/100),[Prange 1]);%special_smooth2D(F, [lambdaSmooth/Prange lambdaSmooth/Xrange],[true true])./special_smooth2D(repmat(XSum*(1/100),[Prange 1]), [lambdaSmooth/Prange lambdaSmooth/Xrange],[true true]);%
                                
                                mat = F;
                            end
                            
                            Fcorrmat = false;
                            if Fcorrmat
                                matref = matref(1:Prange,1:Prange);
                                if PlotVar.Fspatialsmooth
                                    mat = ((mat'*matref)./(sqrt(sum(mat'.^2,2))*sqrt(sum(matref.^2,1))))';
                                    mattemp = conv2(repmat(mat,3,3), 1/9*ones(3), 'same');
                                end
                                mat = mattemp(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
                            end
                            
                            
                            if PlotVar.Fspatialsmooth
                                mat = special_smooth2D(mat, [lambdaSmooth/Prange lambdaSmooth/Xrange],[true true]);
%                                 mat(isnan(mat)) = 0;
%                                 G = smooth1D(repmat(mat,3,3),lambdaSmooth);
%                                 H = smooth1D(G',lambdaSmooth)';
%                                 mat = H(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
                                
                                %                             Z = filter2D(repmat(mat,3,3),lambdaSmooth);
                                %                             mat = Z(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
                                
                                %                             mattemp = conv2(repmat(mat,3,3), 1/9*ones(3), 'same');
                                %                             mat = mattemp(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
                            end
                            
                            if PlotVar.Fnormalize
                                mat = (mat - repmat(min(mat,[],1),size(mat,1),1))./(ones(size(mat,1),1)*(max(mat,[],1)-min(mat,[],1)));
                            end
                            
                            assignin('base',['mat' num2str(g)],mat);
                            
                            if PlotVar.FdisplayPredAve || PlotVar.FdisplayPredMax
                                mattemp = mat;
                                
                                predave = getCircularAverage(mattemp,ampth,maxtol,interpol);
%                                 mattemp([1:20 80:end],:) = 0;
%                                 mattemp = mattemp./repmat(sum(mattemp,1),[Prange 1]);
%                                 mattemp = cumsum(mattemp,1);%mat;%special_smooth2D(mat, [lambdaSmooth/Prange lambdaSmooth/Xrange],[true true]);%mat;
%                                 predave = 0;
%                                 for i = 1:Xrange
%                                     predave(i) = find(diff(sign(mattemp(:,i)-0.5)>0),1,'first');%
%                                 end

%                                 predave = repmat(predave,[3 1]);
%                                 predave = mod(special_smooth_1d(unwrap(predave/Prange*2*pi), lambdaSmooth/Xrange, [], Xrange)*Prange/(2*pi),Prange);
%                                 predave = predave(Xrange+1:2*Xrange);
                                assignin('base',['predave' num2str(g)],predave);
                            end
                            if PlotVar.FdisplayPredMax
                                [~, predmax] = max(mat,[],1);
                                predmax = predmax';
                                predplus = predmax + Prange;
                                predminus = predmax - Prange;
                                xx = (1:Xrange)';
                                predmax = predmax.* (abs(xx-predmax)<abs(xx-predplus) & abs(xx-predmax)<abs(xx-predminus)) + predplus.* (abs(xx-predplus)<abs(xx-predmax) & abs(xx-predplus)<abs(xx-predminus)) + predminus.* (abs(xx-predminus)<abs(xx-predmax) & abs(xx-predminus)<abs(xx-predplus));
                            end
                            
                            if PlotVar.Fdisplaylog
                                mat = log2(mat);
                            end
                            
                            numbins = Xrange;
                            x = 1:numbins;%predaveref;%
                            if ~isempty(mat)
                                if PlotVar.FdisplayMat
                                    PlotVar.Plots{nplot}.PlotMatrix(x,[], mat, [], true, 'CLim', PlotVar.Clim, 'PlotBoxAspectRatio', [1 1 1]);
                                    PlotVar.Plots{nplot}.palette = PlotVar.Palettename;
                                end
                                PlotVar.Plots{nplot}.PlotVector([0 Xrange], [0 Prange], [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k');
                                PlotVar.Plots{nplot}.PlotVector([0 Xrange], [Prange/2 Prange/2], [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'linestyle', '--');
                                try
                                    PlotVar.Plots{nplot}.PlotVector([0 Xrange*min(1,gainratio)], [0 Xrange/max(1,gainratio)], [],[],  true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'linestyle', '--');
                                catch
                                end
                                if PlotVar.FdisplayPredAve
                                    numbins = numel(predave);
                                    PlotVar.Plots{nplot}.PlotVector(x, predave, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'Ylim', [1 Prange], 'linewidth', lwidth);
                                    PlotVar.Plots{nplot}.PlotVector(x(predave < 1), predave(predave < 1) + numbins, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 Prange], 'linewidth', lwidth);
                                    PlotVar.Plots{nplot}.PlotVector(x(predave > numbins), predave(predave > numbins) - numbins, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 Prange], 'linewidth', lwidth);
                                end
                                if PlotVar.FdisplayPredMax
                                    numbins = numel(predmax);
                                    PlotVar.Plots{nplot}.PlotVector(x, predmax, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 Prange], 'linewidth', lwidth);
                                    PlotVar.Plots{nplot}.PlotVector(x(predmax < 1), predmax(predmax < 1) + numbins, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 Prange], 'linewidth', lwidth);
                                    PlotVar.Plots{nplot}.PlotVector(x(predmax > numbins), predmax(predmax > numbins) - numbins, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 Prange], 'linewidth', lwidth);
                                end
                                set(gca,'Xlim', [1 max(X)]);
                            end
                        end
                        if strcmp(PlotVar.ChosenObj{k}, 'Error x X') || strcmp(PlotVar.ChosenObj{k}, 'Error x D')
                            if ~PlotVar.Foverlap
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1)*size(PlotVar.ChosenContrast,2) + (c-1)*size(PlotVar.ChosenObj,2) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            else
                                nplot = (iprobe-1) + k;
                                if numel(PlotVar.Plots) < nplot
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1) + k, 1});
                                end
                            end
                            
                            if strcmp(PlotVar.ChosenObj{k}, 'Error x X')
                                X = EXP.Bayes.Xsmth0{Xtype};%EXP.Bayes.predicted{Xtype};%
                                if ~PlotVar.FthetaPost
                                    Err = EXP.Bayes.PosError0{selprobe,Posttype};%circshift(EXP.Bayes.PosError0{selprobe},[0 1]);%medfilt1(EXP.Bayes.Posterior,16);%
                                else
                                    Err = EXP.Bayes.PosError{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                                end
                                Prange = size(Err,2);
                                tex2 = 0.416 * Prange;
                                tex3 = 0.584 * Prange;
                            elseif strcmp(PlotVar.ChosenObj{k}, 'Error x D')
                                X = EXP.Bayes.D;
                                if ~PlotVar.FthetaPost
                                    Err = EXP.Bayes.DistError0{selprobe};%medfilt1(EXP.Bayes.Posterior,16);%
                                else
                                    Err = EXP.Bayes.DistError{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                                end
                                Prange = size(Err,2);
                                tex2 = 0.416 / unique(es.gain(tidx)/mode(es.gain)) * Prange;
                                tex3 = 0.584 / unique(es.gain(tidx)/mode(es.gain)) * Prange;
                            end
                            
                            Xrange = max(X(tidx));
                            if PlotVar.Fgoodtimebins
                                [goodidx, Xpred] = getcleandecidx(Err, 0*X + floor(Xrange/2), clean_th, maxtol, es.trialID);
                            else
                                goodidx = true(size(tidx));
                            end
                            
%                             [~,pred] = max(Err,[],2);
%                             for tt = 1:size(Err,1)
%                                 Err(tt,:) = circshift(Err(tt,:),floor(Prange/2)-pred(tt));
%                             end
                            
                            mat = zeros(Prange,Xrange);
                            aheaddist = zeros(1, Xrange);
                            behinddist = zeros(1, Xrange);
                            for i = 1:Xrange
                                mat(:,i) = nanmean(Err(X == i  & tidx & spdidx & phsidx & goodidx,:),1)';
                            end
                            
                            if PlotVar.Fspatialsmooth
                                G = smooth1D(repmat(mat,3,3),lambdaSmooth);
                                H = smooth1D(G',lambdaSmooth)';
                                mat = H(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
                            end
                            
                            assignin('base',['mat' num2str(g)],mat);
                            
                            for i = 1:Xrange
                                matPos = mat(:,i);
                                ahead = find(matPos > 1, 1, 'last');%find(matPos > min(matPos) + 0.5*(max(matPos)-min(matPos)), 1, 'last');%find((cumsum(matPos-min(matPos))/sum(matPos-min(matPos)))>=0.95,1,'first');
                                if ~isempty(ahead)
                                    aheaddist(i) = ahead;
                                end
                                behind = find(matPos > 1, 1, 'first');%find(matPos > min(matPos) + 0.5*(max(matPos)-min(matPos)), 1, 'first');%find((cumsum(matPos-min(matPos))/sum(matPos-min(matPos)))>=0.05,1,'first');
                                if ~isempty(behind)
                                    behinddist(i) = behind;
                                end
                            end
                            if PlotVar.FdisplayPredAve
                                mattemp = mat;
                                predave = getCircularAverage(mattemp,ampth,maxtol,interpol);
                            end
                            if PlotVar.FdisplayPredMax
                                [~, predmax] = max(mat,[],1);predmax = predmax';
                                predplus = predmax + Prange;
                                predminus = predmax - Prange;
                                X = (1:Xrange)';
                            end
                            aheaddist = aheaddist - predave'+60;
                            behinddist = predave - behinddist';
                            aheaddist = smooth(aheaddist,3);
                            behinddist = smooth(behinddist,3);
                            
                            if PlotVar.Fdisplaylog
                                mat = log2(mat);
                            end
                            if PlotVar.Fnormalize
                                mat = mat./(ones(size(mat,1),1)*max(mat,[],1));
                            end
                            if PlotVar.FdisplayMat
                                PlotVar.Plots{nplot}.PlotMatrix([],[], mat, [], true, 'CLim', PlotVar.Clim, 'PlotBoxAspectRatio', [1 1 1], 'Xlim', [1 Xrange]);
                                PlotVar.Plots{nplot}.palette = PlotVar.Palettename;
                            end
                            PlotVar.Plots{nplot}.PlotVector([tex2 tex2], [0 Prange], [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'linestyle', '--');
                            PlotVar.Plots{nplot}.PlotVector([tex3 tex3], [0 Prange], [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'linestyle', '--');
                            PlotVar.Plots{nplot}.PlotVector([0 Xrange], [floor(Prange/2) floor(Prange/2)], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'linestyle', '--');
                            
                            if PlotVar.FdisplayPredAve
                                PlotVar.Plots{nplot}.PlotVector([], predave, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 Prange]);
                            end
                            if PlotVar.FdisplayPredMax
                                PlotVar.Plots{nplot}.PlotVector([], predmax, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 Prange]);
                            end
                        end
                        if strcmp(PlotVar.ChosenObj{k}, 'Posterior x Time') || strcmp(PlotVar.ChosenObj{k}, 'Error x Time')
                            if ~PlotVar.Foverlap
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1)*size(PlotVar.ChosenContrast,2) + (c-1)*size(PlotVar.ChosenObj,2) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            else
                                nplot = (iprobe-1) + k;
                                if numel(PlotVar.Plots) < nplot
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1) + k, 1});
                                end
                            end
                            if strcmp(PlotVar.ChosenObj{k}, 'Posterior x Time')
                                if ~PlotVar.FthetaPost
                                    Post = EXP.Bayes.Posterior0{selprobe,Posttype};%medfilt1(EXP.Bayes.Posterior,16);%EXP.data.es.spikeTrain(:,EXP.Bayes.DecCellidx{selprobe});%
                                else
                                    Post = EXP.Bayes.Posterior{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                                end
                                tex2 = 0.416 * size(EXP.Bayes.PosError0{selprobe,Posttype},2);
                                tex3 = 0.584 * size(EXP.Bayes.PosError0{selprobe,Posttype},2);
                            elseif strcmp(PlotVar.ChosenObj{k}, 'Error x Time')
                                if ~PlotVar.FthetaPost
                                    Post = EXP.Bayes.PosError0{selprobe,Posttype};%medfilt1(EXP.Bayes.Posterior,16);%
                                else
                                    Post = EXP.Bayes.PosError{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                                end
                                tex2 = floor(size(EXP.Bayes.PosError0{selprobe},2)/2);
                                tex3 = floor(size(EXP.Bayes.PosError0{selprobe},2)/2);
                            end
                            Prange = size(Post,2);
                            
                            if PlotVar.Fspatialsmooth
                                G = smooth1D(repmat(Post,1,3)',lambdaSmooth)';
                                Post = G(:,Prange+1:2*Prange);
                            end
                            
                            if PlotVar.FdisplayPredAve
                                [~, Xpred] = getcleandecidx(Post, EXP.Bayes.Xsmth0{Xtype}, clean_th, maxtol, es.trialID);
                                predave = Xpred(tidx);
                            end
                            
                            if PlotVar.FdisplayPredMax
                                X = EXP.Bayes.Xsmth0{Xtype};
                                Xrange = min(max(X),size(EXP.Bayes.Posterior0{selprobe},2));
                                [~, predmax] = max(Post(tidx,:),[],2);%pred = pred';
                                if strcmp(PlotVar.ChosenObj{k}, 'Posterior x Time')
                                    predplus = predmax + Xrange;
                                    predminus = predmax - Xrange;
                                    Xtemp = X(tidx);
                                    predmax = predmax.* (abs(Xtemp-predmax)<abs(Xtemp-predplus) & abs(Xtemp-predmax)<abs(Xtemp-predminus)) + predplus.* (abs(Xtemp-predplus)<abs(Xtemp-predmax) & abs(Xtemp-predplus)<abs(Xtemp-predminus)) + predminus.* (abs(Xtemp-predminus)<abs(Xtemp-predmax) & abs(Xtemp-predminus)<abs(Xtemp-predplus));
                                end
                            end
                            
                            if PlotVar.Fdisplaylog
                                Post = log2(Post);
                            end
                            if PlotVar.Fnormalize
                                Post = (Post-repmat(min(Post,[],2),1,size(Post,2)))./((max(Post,[],2)-min(Post,[],2))*ones(1,size(Post,2)));
                            end
                            
                            if PlotVar.FdisplayMat
                                PlotVar.Plots{nplot}.PlotMatrix([],[], Post(tidx,:)', [], true, 'CLim', PlotVar.Clim, 'PlotBoxAspectRatio', [1 1 1]);
                                PlotVar.Plots{nplot}.palette = PlotVar.Palettename;
                            end
                            PlotVar.Plots{nplot}.PlotVector([], EXP.Bayes.Xsmth0{Xtype}(tidx), [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k');
                            PlotVar.Plots{nplot}.PlotVector([], EXP.Bayes.Xsmth0{Xtype}(tidx)./(EXP.data.es.gain(tidx)/mode(EXP.data.es.gain)), [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k');
                            if isfield(es,'lick')
                                PlotVar.Plots{nplot}.PlotVector((find(es.lick(tidx) == 1 & es.outcome(tidx) == 2)), EXP.Bayes.Xsmth0{Xtype}(tidx & es.lick == 1 & es.outcome == 2), [], [], true, 'Marker', '+', 'MarkerEdgeColor', 'g', 'LineStyle', 'none');
                                PlotVar.Plots{nplot}.PlotVector((find(es.lick(tidx) == 1 & es.outcome(tidx) == 0)), EXP.Bayes.Xsmth0{Xtype}(tidx & es.lick == 1 & es.outcome == 0), [], [], true, 'Marker', '+', 'MarkerEdgeColor', 'r', 'LineStyle', 'none');
                            end
                            PlotVar.Plots{nplot}.PlotVector([0 sum(tidx)], [tex2 tex2], [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'linestyle', '--');
                            PlotVar.Plots{nplot}.PlotVector([0 sum(tidx)], [tex3 tex3], [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'linestyle', '--');
                            
                            if PlotVar.FdisplayPredAve
                                PlotVar.Plots{nplot}.PlotVector([], predave, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});%, 'linewidth', 2);
                                
                            end
                            if PlotVar.FdisplayPredMax
                                PlotVar.Plots{nplot}.PlotVector([], predmax, [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                            end
                            if PlotVar.FdisplayLFP
                                PlotVar.Plots{nplot}.PlotVector([], mod(es.LFPphase(tidx,PlotVar.thetaChannel)-180,360)/360*Prange, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'c','linewidth',1);
                            end
                        end
                        if strcmp(PlotVar.ChosenObj{k}, 'Theta x Error x X') || strcmp(PlotVar.ChosenObj{k}, 'Theta x Post x X')
                            xbins = PlotVar.nthetaXbins;
                            roomlength = size(EXP.Bayes.Posterior0{selprobe},2);
                            if ~PlotVar.Foverlap
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1)*size(PlotVar.ChosenContrast,2) + (c-1)*size(PlotVar.ChosenObj,2) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            else
                                nplot = (iprobe-1) + k;
                                if numel(PlotVar.Plots) < nplot
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1) + k, 1});
                                end
                            end
                            
                            if strcmp(PlotVar.ChosenObj{k}, 'Theta x Post x X') || strcmp(PlotVar.ChosenObj{k}, 'Theta time x Post x X')
                                FErr = false;
                                if ~PlotVar.FthetaPost
                                    Post = EXP.Bayes.Posterior0{selprobe,Posttype};%medfilt1(EXP.Bayes.Posterior,16);%
                                else
                                    Post = EXP.Bayes.Posterior{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                                end
                                tex2 = 0.416 * size(EXP.Bayes.PosError0{selprobe},2);
                                tex3 = 0.584 * size(EXP.Bayes.PosError0{selprobe},2);
                                if strcmp(PlotVar.ChosenObj{k}, 'Theta x Post x X')
                                    Fusephase = true;
                                else
                                    Fusephase = false;
                                end
                            elseif strcmp(PlotVar.ChosenObj{k}, 'Theta x Error x X') || strcmp(PlotVar.ChosenObj{k}, 'Theta time x Error x X')
                                FErr = true;
                                if ~PlotVar.FthetaPost
                                    Post = EXP.Bayes.PosError0{selprobe,Posttype};%medfilt1(EXP.Bayes.Posterior,16);%
                                else
                                    Post = EXP.Bayes.PosError{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                                end
                                
                                tex2 = floor(size(EXP.Bayes.PosError0{selprobe},2)/2);
                                tex3 = floor(size(EXP.Bayes.PosError0{selprobe},2)/2);
                                if strcmp(PlotVar.ChosenObj{k}, 'Theta x Error x X')
                                    Fusephase = false;%true;%
                                else
                                    Fusephase = false;
                                end
                            end
                            
                            X = EXP.Bayes.Xsmth0{Xtype};
                            Prange = size(Post,2);
                            Xrange = max(X);
                            if PlotVar.Fgoodtimebins
                                if ~FErr
                                    [goodidx, Xpred, Err] = getcleandecidx(Post, X,clean_th, maxtol, es.trialID);
                                else
                                    [goodidx, Xpred, Err] = getcleandecidx(Post, 0*X +  + floor(max(X)/2),clean_th, maxtol, es.trialID);
                                end
                            else
                                goodidx = true(size(tidx));
                            end
                            if Fusephase
                                nbphsbins = PlotVar.nthetaphsbins;
                                phsval = round(mod(EXP.Bayes.LFPphase{1},360));%mod(EXP.Bayes.LFPphase{selprobe},360);%mod(es.LFPphase(:,PlotVar.thetaChannel),360);%round(mod(EXP.Bayes.LFPphase{selprobe},360));%round(mod(theta.phase',360));%mod(es.LFPphase{selprobe}(:,PlotVar.thetaChannel)-180,360);%mod(es.LFPphase{selprobe}(:,PlotVar.thetaChannel),360);%mod(EXP.Bayes.LFPphase{selprobe}-180,360);mod(EXP.Bayes.LFPphase{selprobe},360);%
                            else
                                nbphsbins = PlotVar.nthetaphsbins;
                                nbtimebins = 18;
                                if selprobe == 1
                                    LFPtheta = es.LFPtheta(:,min(end,EXP.Bayes.thetaChannel));
                                elseif selprobe == 2
                                    LFPtheta = es.LFPtheta(:,min(end,EXP.Bayes.thetaChannel));
                                end
                                minidx0 = find([0;abs(diff(mod(EXP.Bayes.LFPphase{1},360)))>180]);%find([0;abs(diff(mod(EXP.Bayes.LFPphase{selprobe},360)))>180]);
                                thetaperiod0 = diff(minidx0);
                                thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
                                minidx0_2 = find([0;abs(diff(mod(EXP.Bayes.LFPphase{1},360)))>180]);
                                thetaAmp = zeros(size(minidx0));
                                for tt = 1:numel(minidx0)
%                                     nexttrough = find(minidx0_2>minidx0(tt),1,'first');
                                    if tt ~= numel(minidx0)
                                        nexttrough = find(minidx0_2>minidx0(tt) & minidx0_2<minidx0(tt+1));
                                    else
                                        nexttrough = find(minidx0_2>minidx0(tt),1,'first');
                                    end
                                    nexttrough = nexttrough(abs(LFPtheta(nexttrough)) == max(abs(LFPtheta(nexttrough))));
                                    if ~isempty(nexttrough)
                                        thetaAmp(tt) = LFPtheta(minidx0(tt)) - LFPtheta(nexttrough);
                                    end
                                end
                                thetamean = mean(thetaAmp);
                                thetasd = std(thetaAmp);
                                if selprobe == 3%2
                                    zth = 1;
                                    minidx0(abs(thetaAmp)<thetamean + zth*thetasd) = [];
                                end
                                
                                phsval = inf(size(EXP.Bayes.LFPphase{selprobe},1),1);
                                for tphs = 1:numel(minidx0)-1
                                    phsval(minidx0(tphs):minidx0(tphs)+nbphsbins-1) = (0:nbphsbins-1)/nbphsbins*360;
                                end
                            end
                            Frepthetha = false;
                            if Frepthetha
                                ttbins = [0:2*nbphsbins]*360/nbphsbins;
                            else
                                ttbins = [0:nbphsbins]*360/nbphsbins;
                            end
                            POdecfull = [];
                            POdecmean = 0;
                            POthetafull = [];
                            Xfull = [];
                            ttbinsfull = [];
                            lastttbins = 0;
                            predavefull = [];
                            predavestdfull = [];
                            predmaxfull = [];
                            [~, Xpred, ~] = getcleandecidx(Post, 0*X +  + floor(max(X)/2),clean_th, maxtol, es.trialID);
                            xbins0 = xbins;
%                             xbins = 1;
                            for xx = 1:xbins%PlotVar.spdbin:PlotVar.spdbin%
                                x1 = (xx-1)*(roomlength/xbins0);x2 = xx*(roomlength/xbins0);
                                if Fusephase
                                    POdec = zeros(nbphsbins,size(Post,2));
                                    POtheta = zeros(nbphsbins,1);
                                    %                                 [POdec, ctrs1, ctrs2, H] = smoothhist2D_corrected([Xpred(tidx & goodidx & X > x1 & X <= x2) phsval(tidx & goodidx & X > x1 & X <= x2)], lambdaSmooth, [Xrange nbphsbins], 1:Xrange, linspace(0,360,nbphsbins));
                                    for phs = 1:nbphsbins
                                        try
                                            phs0 = (phs-1)*360/nbphsbins;
                                            POdec(phs,:) = nanmean(Post(ismember(phsval,mod(phs0:(phs0 + 360/nbphsbins),360)) & tidx & spdidx & goodidx & X > x1 & X <= x2,:));
                                            POtheta(phs) = nanmean(EXP.Bayes.LFPphase{1}(ismember(phsval,mod(phs0:(phs0 + 360/nbphsbins),360)) & tidx & spdidx & goodidx & X > x1 & X <= x2));
                                        catch
                                            keyboard
                                        end
                                    end
                                else
                                    POdec = [];
                                    POtheta = 0;
                                    phsidx = false(size(tidx));
                                    phsidx(minidx0) = true;
                                    idx = find(tidx & spdidx & X > x1 & X <= x2);
                                    minidx = minidx0(ismember(minidx0,idx));
                                    thetaperiod =  thetaperiod0(ismember(minidx0,idx));
                                    goodidx = single(~isnan(sum(Post,2)));
                                    for tt = 1:numel(minidx)
                                        if minidx(tt) > 1 && minidx(tt)+nbtimebins-1 < size(Post,1) %&& sum(es.badlick(max(1,minidx(tt)):min(minidx(tt)+20*nbphsbins,end))) > 0 %&& sum(es.goodlick(max(1,minidx(tt)-24*nbphsbins):min(minidx(tt),end))) > 0
                                            POmap = Post(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),:);%.*(goodidx(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end))*ones(1,roomlength));
                                            POdec(:,:,tt) = interp1((0:nbtimebins-1)/thetaperiod(tt),POmap,0:1/nbphsbins:(1-(1/nbphsbins)));%POmap;
                                            POthetavec = es.LFPtheta(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end),min(end,EXP.Bayes.thetaChannel));%.*goodidx(max(1,minidx(tt)):min(minidx(tt)+nbtimebins-1,end));%/sum(minidx > 1 & minidx+nbtimebins-1 < size(Post,1));
                                            POtheta = POtheta + (interp1((0:nbtimebins-1)/thetaperiod(tt),POthetavec,0:1/nbphsbins:(1-(1/nbphsbins))))';%POthetavec;
%                                           POtheta = POtheta + ((EXP.Bayes.LFPphase{selprobe}(max(1,minidx(tt)):min(minidx(tt)+nbphsbins-1,end)).*goodidx(max(1,minidx(tt)):min(minidx(tt)+nbphsbins-1,end)))/sum(minidx > 1 & minidx+nbphsbins-1 < size(Post,1)));
                                            %                                     POtheta = POtheta + (l(max(1,minidx(tt)-nbphsbins):min(minidx(tt)+nbphsbins,end))/sum(minidx-nbphsbins > 1 & minidx+nbphsbins < size(p,1)));
                                        end
                                    end
                                    POdecall = POdec;
                                    POdec = nanmean(POdecall,3);
                                    POdecSum = nansum(POdecall,3);
                                    XdecSum = nansum(~isnan(POdecall),3);
%                                     POdec = special_smooth_2d(POdecSum, [1/nbphsbins 1/Prange], [], [], [nbphsbins Prange])./special_smooth_2d(XdecSum, [1/nbphsbins 1/Prange], [], [], [nbphsbins Prange]);
                                    POdec = POdecSum./XdecSum;
                                end
                                %                             [~, maxidx] = max(POdec,[],2);
                                %                             [~,idxascend] = sort(maxidx,'ascend');
                                %                             POdec = POdec(idxascend,:);
                                
                                if Frepthetha
                                    POdec = [POdec(floor(nbphsbins/2):end,:); POdec; POdec(1:floor(nbphsbins/2),:)];
                                    POtheta = [POtheta(floor(nbphsbins/2):end); POtheta; POtheta(1:floor(nbphsbins/2))];
                                else
                                    POdec = [POdec; POdec(1,:)];
                                    POtheta = [POtheta; POtheta(1)];
                                end
                                
                                if PlotVar.Fdisplaylog
                                    POdec = log2(POdec);
                                end
                                if PlotVar.Fspatialsmooth
                                    G = smooth1D(repmat(POdec,3,3),lambdaSmooth);%/size(POdec,2)*size(POdec,1));
                                    H = smooth1D(G',lambdaSmooth)';
                                    POdec = H(size(POdec,1)+1:2*size(POdec,1),size(POdec,2)+1:2*size(POdec,2));
                                    %                                 POdec = conv2(repmat(POdec,3,3),1/9*ones(3),'same');
                                    %                                 POdec = POdec(floor(size(POdec,1)/3)+1:2*floor(size(POdec,1)/3),floor(size(POdec,2)/3)+1:2*floor(size(POdec,2)/3));
                                end
                                if PlotVar.Fnormalize
                                    for i = 1:size(POdec,1)
                                        POdec(i,:) = (POdec(i,:)-min(POdec(i,:)))/(max(POdec(i,:))-min(POdec(i,:)));
                                    end
                                end
                                
                                %                             for i = 1:size(POdec,1)
                                %                                 POdec(i,:) = circshift(POdec(i,:),50 - round(nanmean(find(POdec(i,:) >= 0.9*max(POdec(i,:))))),2);
                                %                             end
                                
%                                 POdecmean = POdecmean+POdec/xbins;
%                                 POcorr = zeros(size(POdec));
%                                 POdecref = nanmean(POdec,1);
%                                 POdecref = POdecref-nanmean(POdecref);
%                                 POdec = POdec - repmat(nanmean(POdec,2),[1 size(POdec,2)]);
%                                 ishift = 0;
%                                 for xshift = -49:50
%                                     ishift = ishift+1;
%                                     POcorr(:,ishift) = sum(POdec.*repmat(circshift(POdecref,xshift),[size(POdec,1) 1]),2);
%                                 end
%                                 POdec = POcorr;
                                
                                POdecfull = [POdecfull POdec'];%[POdecfull POdecmean'];%
                                Xfull = [Xfull (x1+(x2-x1)/2)*ones(1,size(POdec',2))];
                                ttbinsfull = [ttbinsfull lastttbins + ttbins];
                                lastttbins = max(ttbinsfull);
                                POthetafull = [POthetafull POtheta'];
                                
                                if PlotVar.FdisplayPredAve
                                    predave = getCircularAverage(POdec',ampth,maxtol,interpol);
%                                     predave = (predave-mean(predave))/std(predave);
                                    
%                                     nfolds = 5;
%                                     allidx = 1:size(POdecall,3);
                                    predavestd = 0;
%                                     for kiter = 1:nfolds
%                                         idxtest = ((kiter-1)*floor(numel(allidx)/nfolds)+1):(kiter)*floor(numel(allidx)/nfolds);
%                                         idxtrain = allidx(~ismember(allidx,idxtest));
%                                         POdeciter = nanmean(POdecall(:,:,idxtrain),3);
%                                         POdeciterSum = nansum(POdecall(:,:,idxtrain),3);
%                                         XdeciterSum = nansum(~isnan(POdecall(:,:,idxtrain)),3);
%                                         POdeciter = special_smooth_2d(POdeciterSum, [1/nbphsbins 1/Prange], [], [], [nbphsbins Prange])./special_smooth_2d(XdeciterSum, [1/nbphsbins 1/Prange], [], [], [nbphsbins Prange]);
%                                         if PlotVar.Fspatialsmooth
%                                             G = smooth1D(repmat(POdeciter,3,3),lambdaSmooth);%/size(POdec,2)*size(POdec,1));
%                                             H = smooth1D(G',lambdaSmooth)';
%                                             POdeciter = H(size(POdeciter,1)+1:2*size(POdeciter,1),size(POdeciter,2)+1:2*size(POdeciter,2));
%                                         end
%                                         POdeciter = [POdeciter; POdeciter(1,:)];
%                                         predaveiter = getCircularAverage(POdeciter',0,maxtol);
% %                                         predaveiter = (predaveiter-mean(predaveiter))/std(predaveiter);
%                                         predavestd = predavestd + (nfolds - 1)/nfolds*(predaveiter - predave).^2;
%                                     end
%                                     predavestd = sqrt(predavestd);

%                                     
%                                     Xrange = size(POdec,2);
%                                     Prange = Xrange;
%                                     postfilt = zeros(Prange,1);
%                                     POdectemp = POdec;%POdecmean;%
%                                     
%                                     [maxval, ~] = max(POdectemp,[],2);
%                                     for i = 1:size(POdectemp,1)
%                                         POdectemp(i,POdectemp(i,:) < (1-maxtol)*maxval(i)) = 0;
%                                     end
%                                     postfilt(1:Prange) = (0:(Prange-1))'/(Prange)*(2*pi)-pi/2;
%                                     a = sum((POdectemp*cos(postfilt)),2);b = sum((POdectemp*sin(postfilt)),2);
%                                     predave = atan2(-a,b)+pi;%./sum(exp(log(2)*Post(tidx,:)),2);
%                                     predave = predave/(2*pi)*Prange + 1;
%                                     predplus = predave + Prange;
%                                     predminus = predave - Prange;
%                                     if ~FErr
%                                         Xsmth0 = (x1+x2)/2;
%                                         predave = predave.* (abs(Xsmth0-predave)<abs(Xsmth0-predplus) & abs(Xsmth0-predave)<abs(Xsmth0-predminus)) + predplus.* (abs(Xsmth0-predplus)<abs(Xsmth0-predave) & abs(Xsmth0-predplus)<abs(Xsmth0-predminus)) + predminus.* (abs(Xsmth0-predminus)<abs(Xsmth0-predave) & abs(Xsmth0-predminus)<abs(Xsmth0-predplus));
%                                     end
                                    
%                                     predave = getCircularAverage(POdectemp',0,maxtol);
                                    
                                    %                                 predavetemp = smooth(repmat(predave,3,1),3);
                                    %                                 predave = predavetemp(size(predave,1)+1:2*size(predave,1));
                                    predavefull = [predavefull predave'];
                                    predavestdfull = [predavestdfull predavestd'];
                                    behindpos(xx) = min(predave);
                                    aheadpos(xx) = max(predave);
                                end
                                if PlotVar.FdisplayPredMax
                                    [~, predmax] = max(POdec,[],2);predmax = predmax';
                                    predplus = predmax + Prange;
                                    predminus = predmax - Prange;
                                    if ~FErr
                                        Xsmth0 = (x1+x2)/2;
                                        predmax = predmax.* (abs(Xsmth0-predmax)<abs(Xsmth0-predplus) & abs(Xsmth0-predmax)<abs(Xsmth0-predminus)) + predplus.* (abs(Xsmth0-predplus)<abs(Xsmth0-predmax) & abs(Xsmth0-predplus)<abs(Xsmth0-predminus)) + predminus.* (abs(Xsmth0-predminus)<abs(Xsmth0-predmax) & abs(Xsmth0-predminus)<abs(Xsmth0-predplus));
                                    end
                                    %                                 predmax = medfilt1(predmax,3);
                                    predmaxfull = [predmaxfull predmax];
                                end
                            end
%                             predavefull = nanmean(reshape(predavefull,[nbphsbins+1 xbins]),2);
                            ttbinsfull = ttbinsfull(1:nbphsbins+1);
                            Xfull = Xfull(1:nbphsbins+1);
                            POthetafull = nanmean(reshape(POthetafull,[nbphsbins+1 xbins]),2);
                            POdecfull = nanmean(reshape(POdecfull,[Xrange nbphsbins+1 xbins]),3);
                            predavefull = getCircularAverage(POdecfull,ampth,maxtol,interpol);
                            
                            assignin('base',['predavetheta' num2str(g)],predavefull);
                            assignin('base',['predavestdtheta' num2str(g)],predavestdfull);
                            %                         assignin('base',['behindpos' num2str(g)],behindpos);
                            %                         assignin('base',['aheadpos' num2str(g)],aheadpos);
                            
                            %                         predavefullref = evalin('base','predavefull');
                            %                         predavefull = predavefull - predavefullref +50;
                            %                         assignin('base','predavefull',predavefull);
                            %                         assignin('base',['predavetheta_corr' num2str(g)],predavefull);
                            %                         assignin('base',['predavetheta' num2str(g)],predavefull);
                            if numel(POdec) > 1
                                if PlotVar.FdisplayMat
                                    PlotVar.Plots{nplot}.PlotMatrix(ttbinsfull,[], POdecfull, [], true, 'CLim', PlotVar.Clim, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20]);
                                    PlotVar.Plots{nplot}.palette = PlotVar.Palettename;
                                end
                                if PlotVar.FdisplayLFP
                                    PlotVar.Plots{nplot}.PlotVector(ttbinsfull,(POthetafull - mean(POthetafull))/max(POthetafull)*10+ size(POdecfull,1) + 10, [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20]);
                                    %                                 [~,maxidx] = findpeaks(POthetafull);
                                    %                                 maxidx(POthetafull(maxidx)< 0) = [];
                                    %                                 for i = 1:numel(maxidx)
                                    %                                     PlotVar.Plots{nplot}.PlotVector([ttbinsfull(maxidx(i)) ttbinsfull(maxidx(i))],[0 size(POdecfull,1) + 20], [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20], 'color',[0.5 0.5 0.5], 'PlotBoxAspectRatio', [1 1 1]);
                                    %                                 end
                                end
                                if strcmp(PlotVar.ChosenObj{k}, 'Theta x Post x X')
                                    PlotVar.Plots{nplot}.PlotVector(ttbinsfull, Xfull, [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20], 'Color','k', 'PlotBoxAspectRatio', [1 1 1]);
                                elseif strcmp(PlotVar.ChosenObj{k}, 'Theta x Error x X')
                                    PlotVar.Plots{nplot}.PlotVector(ttbinsfull, Xfull*0 + size(POdecfull,1)/2, [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20], 'Color','k', 'PlotBoxAspectRatio', [1 1 1]);
                                end
                                if PlotVar.FdisplayPredAve
                                    %                                 Xrange = size(POdec,2);
                                    %                                 Prange = Xrange;
                                    %                                 postfilt = zeros(Prange,1);
                                    %                                 postfilt(1:Prange) = (0:(Prange-1))'/(Prange)*(2*pi)-pi/2;
                                    %                                 a = sum((POdec*cos(postfilt)),2);b = sum((POdec*sin(postfilt)),2);
                                    %                                 predave = atan2(-a,b)+pi;%./sum(exp(log(2)*Post(tidx,:)),2);
                                    %                                 predave = predave/(2*pi)*Prange;
                                    %                                 predplus = predave + Prange;
                                    %                                 predminus = predave - Prange;
                                    %                                 Xsmth0 = (x1+x2)/2;
                                    %                                 predave = predave.* (abs(Xsmth0-predave)<abs(Xsmth0-predplus) & abs(Xsmth0-predave)<abs(Xsmth0-predminus)) + predplus.* (abs(Xsmth0-predplus)<abs(Xsmth0-predave) & abs(Xsmth0-predplus)<abs(Xsmth0-predminus)) + predminus.* (abs(Xsmth0-predminus)<abs(Xsmth0-predave) & abs(Xsmth0-predminus)<abs(Xsmth0-predplus));
                                    %                                 predave = medfilt1(predave,3);
                                    PlotVar.Plots{nplot}.PlotVector(ttbinsfull,predavefull, [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                                    PlotVar.Plots{nplot}.PlotVector(ttbinsfull(predavefull <= 0),predavefull(predavefull <= 0) + Xrange, [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                                    PlotVar.Plots{nplot}.PlotVector(ttbinsfull(predavefull >= Xrange),predavefull(predavefull >= Xrange) - Xrange, [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                                end
                                
                                if PlotVar.FdisplayPredMax
                                    %                                 [~, predmax] = max(POdec,[],2);predmax = predmax';
                                    %                                 predplus = predmax + Prange;
                                    %                                 predminus = predmax - Prange;
                                    %                                 Xsmth0 = (x1+x2)/2;
                                    %                                 predmax = predmax.* (abs(Xsmth0-predmax)<abs(Xsmth0-predplus) & abs(Xsmth0-predmax)<abs(Xsmth0-predminus)) + predplus.* (abs(Xsmth0-predplus)<abs(Xsmth0-predmax) & abs(Xsmth0-predplus)<abs(Xsmth0-predminus)) + predminus.* (abs(Xsmth0-predminus)<abs(Xsmth0-predmax) & abs(Xsmth0-predminus)<abs(Xsmth0-predplus));
                                    %                                 predmax = medfilt1(predmax,3);
                                    PlotVar.Plots{nplot}.PlotVector(ttbinsfull,predmaxfull, [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                                    %                                 PlotVar.Plots{nplot}.PlotVector(ttbinsfull(predmaxfull <= 0),predmaxfull(predmaxfull <= 0) + Xrange, [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                                    %                                 PlotVar.Plots{nplot}.PlotVector(ttbinsfull(predmaxfull >= Xrange),predmaxfull(predmaxfull >= Xrange) - Xrange, [], [], true, 'Xlim',[ttbinsfull(1) ttbinsfull(end)], 'Ylim',[0 size(POdecfull,1) + 20], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                                end
                            end
                            %                         set(gca,'Ylim',[42 58]);
                        end
                        if strcmp(PlotVar.ChosenObj{k}, 'Projection 45')
                            if ~PlotVar.Foverlap
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1)*size(PlotVar.ChosenContrast,2) + (c-1)*size(PlotVar.ChosenObj,2) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            else
                               nplot = (iprobe-1) + k;
                                if numel(PlotVar.Plots) < nplot
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1) + k, 1});
                                end
                            end
                            predave1 = evalin('base','predave1');
                            predave2 = evalin('base','predave2');
                            predave3 = evalin('base','predave3');
                            
                            predave1_interp = interp1(linspace(0,numel(predave1),numel(predave1)+1), [predave1(1);predave1], 0:0.1:99.9);predave1_interp(isnan(predave1_interp)) = 0;
                            predave2_interp = interp1(linspace(0,numel(predave2),numel(predave2)+1), [predave2(1);predave2], 0:0.1:99.9);predave2_interp(isnan(predave2_interp)) = 0;
                            predave3_interp = interp1(linspace(0,numel(predave3),numel(predave3)+1), [predave3(1);predave3], 0:0.1:99.9);predave3_interp(isnan(predave3_interp)) = 0;
                            predave_1 = zeros(size(predave1_interp));
                            predave_2 = zeros(size(predave2_interp));
                            predave_3 = zeros(size(predave3_interp));
                            x = 0:0.1:99.9;
                            for i = 1:numel(x)
                                predave_1(i) = x(round(mean(find(abs(predave2_interp-predave1_interp(i)) <= min(abs(predave2_interp-predave1_interp(i))))))) - x(i);
                                predave_2(i) = x(round(mean(find(abs(predave2_interp-predave2_interp(i)) <= min(abs(predave2_interp-predave2_interp(i))))))) - x(i);
                                predave_3(i) = x(round(mean(find(abs(predave2_interp-predave3_interp(i)) <= min(abs(predave2_interp-predave3_interp(i))))))) - x(i);
                            end
                            %                         predave_1 = smooth(predave_1,20);
                            %                         predave_2 = smooth(predave_2,20);
                            %                         predave_3 = smooth(predave_3,20);
                            assignin('base','predave_1',predave_1);
                            assignin('base','predave_2',predave_2);
                            assignin('base','predave_3',predave_3);
                            
                            PlotVar.Plots{nplot}.PlotVector(x, predave_1, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,1),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 100], 'linewidth', lwidth);
                            PlotVar.Plots{nplot}.PlotVector(x, predave_2, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,2),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 100], 'linewidth', lwidth);
                            PlotVar.Plots{nplot}.PlotVector(x, predave_3, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,3),PlotVar.ChosenContrast(:,c)}, 'Ylim', [1 100], 'linewidth', lwidth);
                            %                         PlotVar.Plots{nplot}.PlotVector([0 100], [50 70], [],[],  true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'Xlim', [0 100], 'Ylim', [0 100]);
                            %                         PlotVar.Plots{nplot}.PlotVector([0 100], [50 30], [],[],  true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'Xlim', [0 100], 'Ylim', [0 100]);
                            PlotVar.Plots{nplot}.PlotVector([0 100], [0 20], [],[],  true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'Xlim', [0 100], 'Ylim', [-50 50]);
                            PlotVar.Plots{nplot}.PlotVector([0 100], [0 -20], [],[],  true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k', 'Xlim', [0 100], 'Ylim', [-50 50]);
                            
                            %                         Xrange = max(EXP.Bayes.Xsmth0{Xtype});
                            %                         nidx = size(EXP.Bayes.Xsmth0{Xtype},1);
                            %                         mat = zeros(Xrange);
                            %                         if ~PlotVar.FthetaPost
                            %                             Post = EXP.Bayes.Posterior0{selprobe};%medfilt1(EXP.Bayes.Posterior,16);%
                            %                         else
                            %                             Post = EXP.Bayes.Posterior{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                            %                         end
                            %                         for i = 1:Xrange
                            %                             mat(:,i) = mean(Post(EXP.Bayes.Xsmth0{Xtype}(1:nidx) == i & tidx,:),1)';
                            %                         end
                            %                         [mat, X, distances] = get45Marginal(mat, 100);
                            %                         PlotVar.Plots{nplot}.PlotVector(X, mat, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                        end
                        if strcmp(PlotVar.ChosenObj{k}, 'Mean X-Error') || strcmp(PlotVar.ChosenObj{k}, 'Mean D-Error')
                            if ~PlotVar.Foverlap
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1)*size(PlotVar.ChosenContrast,2) + (c-1)*size(PlotVar.ChosenObj,2) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            else
                                nplot = (iprobe-1) + k;
                                if numel(PlotVar.Plots) < nplot
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1) + k, 1});
                                end
                            end
                            Prange = size(EXP.Bayes.PosError0{selprobe},2);
                            if strcmp(PlotVar.ChosenObj{k}, 'Mean X-Error')
                                if ~PlotVar.FthetaPost
                                    Err =  EXP.Bayes.PosError0{selprobe,Posttype};%./repmat(max(EXP.Bayes.PosError0{selprobe},[],2),[1 size(EXP.Bayes.PosError{selprobe},2)]);%EXP.Bayes.PosError0{selprobe};%medfilt1(EXP.Bayes.Posterior,16);%evalin('base','PosErrorCentered');%
                                else
                                    Err = EXP.Bayes.PosError{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                                end
                                
                                %                             Prange = size(EXP.Bayes.PosError0{selprobe},2);
                                %                             for tt = 1:size(EXP.Bayes.PosError0{selprobe},1)
                                %                                 [~,maxpos] = max(EXP.Bayes.PosError0{selprobe}(tt,:));
                                %                                 EXP.Bayes.PosError0{selprobe}(tt,:) = circshift(EXP.Bayes.PosError0{selprobe}(tt,:),floor(Prange/2)-maxpos,2);
                                %                             end
                                %                             assignin('base','PosErrorCentered',EXP.Bayes.PosError0{selprobe});
                                
                                
                                %                             idx = find(es.firstgoodlick);
                                %                             beforelick = zeros(size(es.firstrewlick));
                                %                             beforelick(max(1,idx - PlotVar.nthetaXbins)) = 1;
                                %                             tidx = tidx & beforelick;%es.firstrewlick;
                                
                                %                               tidx = tidx & ismember(EXP.Bayes.Xsmth0{Xtype},PlotVar.nthetaXbins:PlotVar.nthetaXbins+5);
                                
                                %                             dx = 4;
                                %                             if unique(EXP.data.es.gain(tidx))/mode(EXP.data.es.gain) > 1
                                %                                 x0 = dx;
                                %                             elseif unique(EXP.data.es.gain(tidx))/mode(EXP.data.es.gain) < 1
                                %                                 x0 = -dx;
                                %                             else
                                %                                 x0 = 0;
                                %                             end
                                %                             Err = circshift(Err,floor(x0),2);
                                
                                if PlotVar.Fgoodtimebins
                                    [goodidx, Xpred] = getcleandecidx(Err, 0*EXP.Bayes.Xsmth0{Xtype} + floor(max(EXP.Bayes.Xsmth0{Xtype})/2),clean_th, maxtol, es.trialID);
                                else
                                    goodidx = true(size(tidx));
                                end
                                if PlotVar.Fspatialsmooth
                                    G = smooth1D(repmat(Err,1,3)',lambdaSmooth)';
                                    Err = G(:,size(Err,2)+1:2*size(Err,2));
                                    %                                 Err = conv2(repmat(Err,3,3),1/9*ones(3),'same');
                                    %                                 Err = Err(:,floor(size(Err,2)/3)+1:2*floor(size(Err,2)/3));
                                end
                                
%                                 Err = Err./repmat(max(Err,[],2),[1 size(Err,2)]);
%                                 for xx = 1:size(Err,2)
%                                     Xpred = nanmean(Err(tidx & spdidx & phsidx & goodidx & EXP.Bayes.Xsmth0{Xtype}==xx,:),1);
%                                     Xpred = getCircularAverage(Xpred',0,maxtol);
%                                     Err(EXP.Bayes.Xsmth0{Xtype}==xx,:) = circshift(Err(EXP.Bayes.Xsmth0{Xtype}==xx,:),round(floor(Prange/2)-Xpred));
%                                 end

%                                 for tt = 1:size(Err,1)
%                                     Xpred = getCircularAverage(Err(tt,:)',0,maxtol);
%                                     Err(tt,:) = circshift(Err(tt,:),round(floor(Prange/2)-Xpred));
%                                 end

                                
%                                 if PlotVar.Fspatialsmooth
%                                     Post = EXP.Bayes.Posterior0{selprobe};
%                                     G = smooth1D(repmat(Post,1,3)',lambdaSmooth)';
%                                     Post = G(:,size(Post,2)+1:2*size(Post,2));
%                                     %                                 Err = conv2(repmat(Err,3,3),1/9*ones(3),'same');
%                                     %                                 Err = Err(:,floor(size(Err,2)/3)+1:2*floor(size(Err,2)/3));
%                                 end
                                [~, Xpred, ErrVec] = getcleandecidx(EXP.Bayes.Posterior0{selprobe}, EXP.Bayes.Xsmth0{Xtype},clean_th, maxtol, es.trialID);
%                                 Prange = size(EXP.Bayes.Posterior0{selprobe},2);
%                                 Xrange = max(EXP.Bayes.Xsmth0{Xtype});
%                                 for i = 1:Xrange
%                                     Err(round(Xpred) == i,:) = circshift(EXP.Bayes.Posterior0{selprobe}(round(Xpred) == i,:),floor(Prange/2)-i,2);
%                                 end
%                                 
                                cc = 0;
                                
                                
                                %                             phsval = mod(round(EXP.Bayes.LFPphase{selprobe})-180,360);%mod(es.LFPphase{selprobe}(:,PlotVar.thetaChannel)-180,360);
                                %                             nthetaphase = PlotVar.nthetaphsbins;
                                %                             MeanErr = zeros(1,nthetaphase);
                                %                             for phs = 1:nthetaphase
                                %                                 pidx = ismember(phsval,mod((phs-1)*360/nthetaphase:(phs*360/nthetaphase),360));
                                %                                 MeanErrRef = nanmean(circshift(Err(tidxref & pidx & goodidx,:),cc,2),1);
                                %                                 MeanErrtemp = nanmean(circshift(Err(tidx & pidx & goodidx,:),cc,2),1);
                                %                                 MeanErr(phs) = find(MeanErrtemp == max(MeanErrtemp), 1, 'last') - find(MeanErrRef == max(MeanErrRef), 1, 'last');
                                % %                                 [~,MeanErr(phs)] = max(xcorr(MeanErrtemp,MeanErrRef,50,'coeff'));
                                % %                                 MeanErr(phs) = MeanErr(phs) - 51;
                                %                             end
                                
                                %                             MeanErrRef = nanmean(Err(tidxref & phsidx & goodidx,:),1);
                                Xseg = 0:100;%0:20;
                                %                             goodidx = goodidx & es.firstgoodlick;
                                %                             idx = find(goodidx);
                                %                             for i = 1:10
                                %                                goodidx(max(1,idx-i)) = true;
                                %                             end
                                Xrange = max(EXP.Bayes.Xsmth0{Xtype});
                                Prange = size(Err,2);
                                MeanErr = zeros(Xrange,Prange);%zeros(Xrange,1);%
                                for xx = 1:Xrange
%                                     MeanErr(xx) = sum(round(ErrVec(tidx & spdidx & phsidx & goodidx))+51 == xx)/sum(tidx & spdidx & phsidx & goodidx);
                                    MeanErr(xx,:) = nanmean(Err(tidx & spdidx & phsidx & goodidx & EXP.Bayes.Xsmth0{Xtype}==xx,:),1);
%                                     Xpred = getCircularAverage(MeanErr(xx,:)',0,maxtol);
%                                     MeanErr(xx,:) = circshift(MeanErr(xx,:),round(floor(Prange/2)-Xpred));
                                end
%                                 MeanErr = MeanErr./repmat(max(MeanErr,[],2),[1 size(MeanErr,2)]);
                                MeanErr = nanmean(MeanErr,1);
                                %                             mat = evalin('base',['mat' num2str(g)]);
                                %                             MeanErr = mat(:, PlotVar.thetaphase);%nanmean(Err(tidx & spdidx & phsidx & ismember(EXP.Bayes.Xsmth0{Xtype},Xseg),:),1);
                                
                                %                             Xshift  = find(MeanErr == max(MeanErr)) - find(MeanErrRef == max(MeanErrRef));
                                %                             [~,Xshift] = max(xcorr(MeanErr,MeanErrRef,10,'coeff'));Xshift = Xshift - 11;
                                %                             MeanErr = circshift(MeanErr,-Xshift,2);
                                if PlotVar.Fnormalize
                                    maxval = max(MeanErr);
                                    minval = min(MeanErr);
                                    for i = 1:numel(MeanErr)
                                        MeanErr(i) = (MeanErr(i)-minval)/(maxval-minval);
                                    end
                                end
                                PlotVar.Plots{nplot}.PlotVector([], MeanErr, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)},'Ylim', [-4 4]);
                                PlotVar.Plots{nplot}.PlotVector([floor(size(Err,2)/2) floor(size(Err,2)/2)], [0 1.5], [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'linestyle', '--', 'Color', 'k');
                                
                            elseif strcmp(PlotVar.ChosenObj{k}, 'Mean D-Error')
                                if ~PlotVar.FthetaPost
                                    Err = EXP.Bayes.DistError0{selprobe};%medfilt1(EXP.Bayes.Posterior,16);%
                                else
                                    Err = EXP.Bayes.DistError{selprobe,phsrefidx};%medfilt1(EXP.Bayes.Posterior,16);%
                                end
                                if PlotVar.Fgoodtimebins
                                    [goodidx, Xpred] = getcleandecidx(Err, EXP.Bayes.D + floor(max(EXP.Bayes.D)/2),clean_th, maxtol, es.trialID);
                                else
                                    goodidx = true(size(tidx));
                                end
                                if PlotVar.Fspatialsmooth
                                    G = smooth1D(repmat(Err,1,3)',lambdaSmooth)';
                                    Err = G(:,size(Err,2)+1:2*size(Err,2));
                                    %                                 Err = conv2(repmat(Err,1,3),1/9*ones(1,3),'same');
                                    %                                 Err = Err(:,floor(size(Err,2)/3)+1:2*floor(size(Err,2)/3));
                                end
                                MeanErr = nanmean(Err(tidx & spdidx & phsidx & goodidx,:),1);
                                if PlotVar.Fnormalize
                                    maxval = max(MeanErr);
                                    for i = 1:numel(MeanErr)
                                        MeanErr(i) = MeanErr(i)/maxval;
                                    end
                                end
                                PlotVar.Plots{nplot}.PlotVector([], MeanErr, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                                PlotVar.Plots{nplot}.PlotVector([floor(size(Err,2)/2) floor(size(Err,2)/2)], [0 2], [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'linestyle', '--', 'Color', 'k');
                            end
                        end
                        if strcmp(PlotVar.ChosenObj{k}, 'Theta x Error phase')
                            if ~PlotVar.Foverlap
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1)*size(PlotVar.ChosenContrast,2) + (c-1)*size(PlotVar.ChosenObj,2) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            else
                                nplot = (iprobe-1) + k;
                                if numel(PlotVar.Plots) < nplot
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (iprobe-1) + k, 1});
                                end
                            end
                            X = EXP.Bayes.Xsmth0{Xtype};
                            if PlotVar.Fgoodtimebins
                                [goodidx, Xpred, Err] = getcleandecidx(EXP.Bayes.Posterior0{selprobe}, X, clean_th, maxtol, es.trialID);
                            else
                                goodidx = true(size(tidx));
                            end
                            Xseg = 0:100;%mod(PlotVar.thetaphase:PlotVar.thetaphase+20,max(X));
                            [F,ctrs1, ctrs2, H] = smoothhist2D_corrected([EXP.Bayes.LFPphase{1}(tidx & goodidx & ismember(X, Xseg)) theta.phase(tidx & goodidx & ismember(X, Xseg))'],2, [PlotVar.nthetaphsbins PlotVar.nthetaphsbins], 0:360/PlotVar.nthetaphsbins:360, 0:360/PlotVar.nthetaphsbins:360, true, true);
                            mat = F;
                            if PlotVar.Fspatialsmooth
                                G = smooth1D(repmat(mat,3,3),lambdaSmooth);
                                H = smooth1D(G',lambdaSmooth)';
                                mat = H(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
                            end
                            if PlotVar.FdisplayMat
                                phsbins = 0:360/PlotVar.nthetaphsbins:360;
                                PlotVar.Plots{nplot}.PlotMatrix(phsbins,phsbins, mat, [], true, 'Xlim', [phsbins(1) phsbins(end)], 'Ylim', [phsbins(1) phsbins(end)], 'PlotBoxAspectRatio', [1 1 1]);
                                PlotVar.Plots{nplot}.PlotVector(phsbins, phsbins, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'k');
                                PlotVar.Plots{nplot}.palette = PlotVar.Palettename;
                                %                             [vec, ~, distances] = get45Marginal(repmat(mat,3,3), 3*PlotVar.nthetaphsbins);
                                %                             PlotVar.Plots{nplot}.PlotVector([phsbins(1:end-1)-180 360+phsbins(1:end-1)-180 360*2+phsbins(1:end-1)-180], vec, [], [], true, 'PlotBoxAspectRatio', [1 1 1], 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                            end
                        end
                    end
                end
            end
        end
    end
end
end

function mat_oo = special_smooth2D(mat_i,win,Fcircular)
% lambdaSmooth = 2;
% nonvalid = isnan(mat_i);
% mat_i(nonvalid) = 1;%1
% G = smooth1D(repmat(mat_i,3,3),lambdaSmooth);
% H = smooth1D(G',lambdaSmooth)';
% mat_o = H(size(mat_i,1)+1:2*size(mat_i,1),size(mat_i,2)+1:2*size(mat_i,2));
% mat_o(nonvalid) = NaN;

mat_oo = mat_i;
for k = 1:size(mat_i,3)
    mat_o = mat_i(:,:,k);
    mat_o(isnan(mat_i(:,:,k))) = 0;
    if Fcircular(1)
        mat_o = repmat(mat_o,[3 1]);
    end
    if Fcircular(2)
        mat_o = repmat(mat_o,[1 3]);
    end
    if sum(isnan(win)) > 0
        if isnan(win(1)) && ~isnan(win(2))
            for i = 1:size(mat_o,1)
                mat_o(i,:) = special_smooth_1d(mat_o(i,:), win(2), [], size(mat_i,2));
            end
        end
        if isnan(win(2)) && ~isnan(win(1))
            for j = 1:size(mat_o,2)
                mat_o(:,j) = special_smooth_1d(mat_o(:,j), win(1), [], size(mat_i,1));
            end
        end
    else
        mat_o = special_smooth_2d(mat_o, win, [], [], [size(mat_i,1) size(mat_i,2)]);
    end
    if Fcircular(1)
        mat_o = mat_o(floor(size(mat_o,1)/3)+1:2*floor(size(mat_o,1)/3),:);
    end
    if Fcircular(2)
        mat_o = mat_o(:,floor(size(mat_o,2)/3)+1:2*floor(size(mat_o,2)/3));
    end
    
    mat_o(isnan(mat_i(:,:,k))) = NaN;
    mat_oo(:,:,k) = mat_o;
end
end