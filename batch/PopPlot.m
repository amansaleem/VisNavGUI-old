function PopPlot(popres,cellprop,reslat,maptype)

titlestr{1} = 'low';
titlestr{2} = 'med';
titlestr{3} = 'high';

gainVal = [0.8 1 1.2];
cl{1} = 'c';
cl{2} = 'g';
cl{3} = 'm';

lambdaSmooth = 2;
probestr{1} = popres.probestr{1};
probestr{2} = popres.probestr{2};

switch maptype
    case 'Posterior'
        maps = popres.PostXAll;
        Xdec = popres.PostXpred;
        Errdec = popres.PostErrpred;
        Errdecstd = popres.PostErrpred_std;
        Errdec_rel = popres.PostErrpred;
        Errdec_relstd = popres.PostErrpred;
        Errdec_reldec = popres.PostErrpred;
        Errdec_reldecstd = popres.PostErrpred;
        Errdec_reg = popres.PostErrpred_reg;
        meanErr = popres.s_meanErrX_Post;
        meanErrSeg = popres.s_meanErrX_Post;%popres.s_meanErrX_aveSeg;
    case 'DistriXdecMax'
        maps = popres.DistriXmaxAll;
        Xdec = popres.DistriXmaxpred;
        Errdec = popres.DistrimaxErrpred;
        Errdecstd = popres.DistrimaxErrpred;
        Errdec_rel = popres.DistrimaxErrpred;
        Errdec_relstd = popres.DistrimaxErrpred;
        Errdec_reldec = popres.DistrimaxErrpred;
        Errdec_reldecstd = popres.DistrimaxErrpred;
        Errdec_reg = popres.DistrimaxErrpred_reg;
        meanErr = popres.s_meanErrX_max;
        meanErrSeg = popres.s_meanErrX_aveSeg;
    case 'DistriXdecAve'
        maps = popres.DistriXaveAll;
        Xdec = popres.DistriXavepred;
        Errdec = popres.DistriaveErrpred;
        Errdecstd = popres.DistriaveErrpred_std;
        Errdec_rel = popres.DistriaveErrpred_rel;
        Errdec_relstd = popres.DistriaveErrpred_relstd;
        Errdec_reldec = popres.DistriaveErrpred_reldec;
        Errdec_reldecstd = popres.DistriaveErrpred_reldecstd;
        Errdec_reg = popres.DistriaveErrpred_reg;
        meanErr = popres.s_meanErrX_ave;
        meanErrSeg = popres.s_meanErrX_aveSeg;
end

nanimal = size(meanErr,1);
nProbe = size(meanErr,3);
ngain = size(popres.PostXAll,2);

meanErrAll = cell(nProbe,ngain);
meanErrAnimal = cell(nProbe,ngain);
meanSpdAll = cell(nProbe,ngain);
meanLickAll = cell(nProbe,ngain);
meanErrAllSeg = cell(nProbe,ngain);
nSeg = 5;
meanLickAnimal = cell(nProbe,ngain);
for iprobe = 1:nProbe
    for g = 1:3
        for ianimal = 1:nanimal
            meanErrSession = 0;
            meanLickSession = 0;
            nSession = 0;
            for iseries = 1:size(meanErr,2)
                if ~isempty(max([popres.s_meanlickX{ianimal,iseries,1,g} popres.s_meanlickX{ianimal,iseries,2,g}]))
                    meanLickAll{iprobe,g} = [meanLickAll{iprobe,g}  max([popres.s_meanlickX{ianimal,iseries,1,g} popres.s_meanlickX{ianimal,iseries,2,g}])];
                else
                    meanLickAll{iprobe,g} = [meanLickAll{iprobe,g}  NaN];
                end
                if ~isempty(meanErr{ianimal,iseries,iprobe,g})
                    meanErrAll{iprobe,g} = [meanErrAll{iprobe,g}  meanErr{ianimal,iseries,iprobe,g}];
                    meanErrAllSeg{iprobe,g} = [meanErrAllSeg{iprobe,g}  meanErrSeg{ianimal,iseries,iprobe,g}'];
                    [speedpro,~,~] = fast1Dmap(popres.s_X{ianimal,iseries,iprobe,g},popres.s_speed{ianimal,iseries,iprobe,g},1,1,[], true);
                    meanSpdAll{iprobe,g} = [meanSpdAll{iprobe,g}  nanmean(speedpro)];
                    if ~isnan(meanErr{ianimal,iseries,iprobe,g})
                        nSession = nSession + 1;
                        meanErrSession = meanErrSession + meanErr{ianimal,iseries,iprobe,g};
                        meanLickSession = meanLickSession + popres.s_meanlickX{ianimal,iseries,iprobe,g};
                    end
                else
                    meanErrAll{iprobe,g} = [meanErrAll{iprobe,g}  NaN];
                    meanSpdAll{iprobe,g} = [meanSpdAll{iprobe,g} NaN];
%                     meanErrAllSeg{iprobe,g} = [meanErrAllSeg{iprobe,g} NaN(nSeg,1)];
                end
            end
            meanErrAnimal{iprobe,g} = [meanErrAnimal{iprobe,g} meanErrSession/nSession];
            meanLickAnimal{iprobe,g} = [meanLickAnimal{iprobe,g} meanLickSession/nSession];
        end
    end
end
X = popres.dx:popres.dx:(size(maps{1,2},2)*popres.dx);
Y = popres.dx:popres.dx:(size(maps{1,2},2)*popres.dx);
Xreg = (popres.dx*popres.dx):(popres.dx*popres.dx):(size(maps{1,2},2)*popres.dx);
for iprobe = 1:nProbe
    figure('name',['Decoding : ' probestr{iprobe}]);
    for g = 1:3
        subplot(3,3,g);
        imagesc(X,Y,(maps{iprobe,g}));
        hold on;plot(X,Xdec{iprobe,2}*popres.dx,'k');plot([0 100],[0 100],'k');
%         hold on;plot(X,Xdec{iprobe,g}*popres.dx,cl{g});
%         hold on;plot(Xdec{iprobe,g},'k');
        colormap(parula);%colormap(RedWhiteBlue);%
        set(gca,'Clim',[0 2],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
        xlabel('actual');
        ylabel('decoded');
        title([probestr{iprobe} ' ' titlestr{g}]);
        
        subplot(3,4,5);
        ax = gca;
        hold on;ciplot(ax,X,Errdec{iprobe,g}*popres.dx  - 50,popres.dx*(Errdecstd{iprobe,g}).^0.5,0.4,cl{g});
%         hold on;plot(X,Errdec{iprobe,g}*popres.dx-Errdec{iprobe,2}*popres.dx,cl{g});
        set(gca,'Xlim',[0 max(X)],'Ylim',[-20 20]);
%         hold on;plot(X/gainVal(g),Errdec{iprobe,g}*popres.dx + X'*(1-1/gainVal(g)),cl{g});
%         set(gca,'Xlim',[X(1)/gainVal(1) X(end)/gainVal(1)],'Ylim',[-20 20]);
        xlabel('position');
        ylabel([probestr{iprobe} ' decoding error']);
        
        subplot(3,4,6);
        ax = gca;
%         hold on;ciplot(ax,X,(Errdec{iprobe,g} - Errdec{iprobe,2})*popres.dx,0*(Errdecstd{iprobe,g}).^0.5,0.4,cl{g});
        hold on;ciplot(ax,X,Errdec_reldec{iprobe,g}*popres.dx,popres.dx*(Errdec_reldecstd{iprobe,g}).^0.5,0.4,cl{g});
        set(gca,'Xlim',[0 max(X)],'Ylim',[-20 20]);
%         [~,isort] = sort(Xdec{iprobe,g},'ascend');
%         err = interp1(Xdec{iprobe,g}(isort),Errdec{iprobe,g}(isort),0:max(Xdec{iprobe,2}),'spline');
%         [~,isort] = sort(Xdec{iprobe,2},'ascend');
%         err_ref = interp1(Xdec{iprobe,2}(isort),Errdec{iprobe,2}(isort),0:max(Xdec{iprobe,2}),'spline');
%         hold on;plot((0:max(Xdec{iprobe,2}))*popres.dx,(err-(err_ref))*popres.dx,cl{g});
%         hold on;plot([0 max(Xdec{iprobe,2})*popres.dx],[0 max(Xdec{iprobe,2})*(1 -gainVal(g))*popres.dx],cl{g},'LineStyle','--');
%         set(gca,'Xlim',[0 max(Xdec{iprobe,2})*popres.dx],'Ylim',[-20 20]);
        
%         hold on;plot(X/gainVal(g),Errdec{iprobe,g}*popres.dx + X'*(1-1/gainVal(g)),cl{g});
%         set(gca,'Xlim',[X(1)/gainVal(1) X(end)/gainVal(1)],'Ylim',[-20 20]);
        xlabel('decoded position');
        ylabel([probestr{iprobe} ' decoding error']);
        
        subplot(3,4,7);
        hold on;plot(Xreg,Errdec_reg{iprobe,g}*popres.dx-Errdec_reg{iprobe,2}*popres.dx,cl{g});
        hold on;plot([min(Xreg) max(Xreg)],[min(Xreg) max(Xreg)*(1 -gainVal(g))],cl{g},'LineStyle','--');
        set(gca,'Xlim',[min(Xreg) max(Xreg)],'Ylim',[-20 20]);
%         hold on;plot(Xreg/gainVal(g),Errdec_reg{iprobe,g}*popres.dx+ Xreg*(1-1/gainVal(g)),cl{g});
%         hold on;plot([Xreg(1)/gainVal(g) Xreg(end)/gainVal(g)],[0 -(1 -gainVal(g))*Xreg(end)/gainVal(g)],cl{g},'LineStyle','--');
%         set(gca,'Xlim',[Xreg(1)/gainVal(1) Xreg(end)/gainVal(1)],'Ylim',[-20 20]);
        xlabel('position');
        ylabel([probestr{iprobe} ' relative error']);
        
        if ismember(g,[1 3])
            subplot(3,4,8);
            hold on;
            hshift = histogram(meanErrAll{iprobe,g}-meanErrAll{iprobe,2},[-inf -10:1:10 inf],'EdgeColor','none','FaceColor',cl{g});
            set(gca,'Xlim',[-20 20]);
            hold on;plot([0 0],[0 15],'k');
            xlabel([probestr{iprobe} ' decoding error']);
            ylabel('# of sessions');
        end
        
        if iprobe == 1
            subplot(3,3,7);
            [speedpro,x,speedprostd] = fast1Dmap(popres.X{1,g},popres.speed{1,g},1,1,[], true);
            hold on;plot(x,speedpro,cl{g});
            hold on;ciplot(gca,x,speedpro,speedprostd,0.1,cl{g});
            xlabel('position');
            ylabel('running speed');
        else
            subplot(3,3,7);
            [eyepro,x,eyeprostd] = fast1Dmap(popres.X{1,g}(~isnan(popres.eyeXpos{1,g})),popres.eyeXpos{1,g}(~isnan(popres.eyeXpos{1,g})),1,1,[], true);
            hold on;plot(x,eyepro,cl{g});
            xlabel('position');
            ylabel('pupil X position');
%             hold on;ciplot(gca,x,eyepro,eyeprostd,0.1,cl{g});
        end
        
        subplot(3,3,9);
        hold on;
        barh(meanErrAnimal{iprobe,g} - meanErrAnimal{iprobe,2},'Facecolor',cl{g},'Linestyle','none');
        set(gca,'Xlim',[-20 20],'Ylim',[0 11]);
        xlabel([probestr{iprobe} ' decoding error']);
        ylabel('animal #');
        
        if ismember(g,[1 3])
            subplot(3,6,16);
            hold on;scatter(meanSpdAll{iprobe,g},meanErrAll{iprobe,g}-meanErrAll{iprobe,2},cl{g});
            hold on;plot([1 70],[0 0]);
            set(gca,'Xlim',[1 70],'Ylim',[-20 20],'PlotBoxAspectRatio', [1 1 1]);
            xlabel('nanmean run speed');
            ylabel([probestr{iprobe} ' error']);
        end
        
        if iprobe == 2 && ismember(g,[1 3])
            subplot(3,6,15);
            hold on;
            scatter(meanErrAll{1,g}-meanErrAll{1,2},meanErrAll{2,g}-meanErrAll{2,2},cl{g});
            
%             CA1gainshift{g} = meanErrAll{1,g}-meanErrAll{1,2};
%             V1gainshift{g} = meanErrAll{2,g}-meanErrAll{2,2};
%             CA1nanidx = isnan(meanErrAll{1,g}-meanErrAll{1,2});
%             V1nanidx = isnan(meanErrAll{2,g}-meanErrAll{2,2});
%             CA1gainshift{g} = CA1gainshift{g}(~CA1nanidx & ~V1nanidx);
%             V1gainshift{g} = V1gainshift{g}(~CA1nanidx & ~V1nanidx);
%             [r,p] = corr(CA1gainshift{g}',V1gainshift{g}');
%             disp(r)
%             disp(p)
            
            hold on;plot([-20 20],[-20 20],'k');
            hold on;plot([0 0],[-20 20],'k');
            hold on;plot([-20 20],[0 0],'k');
            set(gca,'Xlim',[-20 20],'Ylim',[-20 20],'PlotBoxAspectRatio', [1 1 1]);
            xlabel('CA1 error');
            ylabel('V1 error');
        end
    end
end 

if nProbe > 1
    figure('Name','V1CA1 noise correlation');
    for g = [2 1 3]
        CA1V1corr{g} = 0;
        CA1V1corr_randXS{g} = 0;
        CA1V12Dcorr{g} = 0;
        CA1V12Dcorr_randXS{g} = 0;
        CA1V12Dcorr_randXSsqr{g} = 0;
    end
    for g = [2 1 3]
        subplot(5,3,g);
        nX = size(popres.CA1V1corr,2);
        for xx = 1:nX
            CA1V1corr{g} = CA1V1corr{g} + popres.CA1V1corr{g,xx}/nX;
            CA1V1corr_randXS{g} = CA1V1corr_randXS{g} + popres.CA1V1corr_randXS{g,xx}/nX;
            CA1V12Dcorr{g} = CA1V12Dcorr{g} + popres.CA1V12Dcorr{g,xx}/nX;
            CA1V12Dcorr_randXS{g} = CA1V12Dcorr_randXS{g} + popres.CA1V12Dcorr_randXS{g,xx}/nX;
            CA1V12Dcorr_randXSsqr{g} = CA1V12Dcorr_randXSsqr{g} + popres.CA1V12Dcorr_randXSsqr{g,xx}/nX;
        end
        zval = (popres.CA1V1corr{g} - popres.CA1V1corr_randXS{g})./(popres.CA1V1corr_randXSsqr{g} - popres.CA1V1corr_randXS{g}.^2).^0.5;
%         scatter(CA1V1corr_randXS{g},CA1V1corr{g},'MarkerEdgeColor',cl{g},'SizeData',10);
        scatter(CA1V1corr_randXS{g}(zval>=2),CA1V1corr{g}(zval>=2),'MarkerEdgeColor',cl{g},'SizeData',10);
        hold on;scatter(CA1V1corr_randXS{g}(zval<2),CA1V1corr{g}(zval<2),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',10);
        hold on;plot([-0.1 0.35],[-0.1 0.35],'k');
        hold on;plot([0 0],[-0.1 0.35],'k');
        hold on;plot([-0.1 0.35],[0 0],'k');
        set(gca,'Xlim',[-0.1 0.35],'Ylim',[-0.1 0.35],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
        title([titlestr{g}]);
        xlabel('due to X and speed');
        ylabel('error correlation');
        
        subplot(5,3,3+g);
        imagesc(CA1V12Dcorr{g})
        %     imagesc(popres.CA1V12Dcorr{g}-popres.CA1V12Dcorr_randS{g})
        set(gca,'Xlim',[1 100],'Ylim',[1 100],'Clim',[0 2],'Ydir','normal','XAxisLocation','origin','YAxisLocation','origin','PlotBoxAspectRatio', [1 1 1])
        hold on;plot([50 50],[1 100],'k');
        hold on;plot([1 100],[50 50],'k');
        hold on;plot([1 100],[1 100],'k');
        title([titlestr{g} '(raw)']);
        %     title([titlestr{g} '(- spd)']);
        xlabel('CA1');
        ylabel('V1');
        
        subplot(5,3,6+g);
        imagesc(CA1V12Dcorr_randXS{g})
        %     imagesc(popres.CA1V12Dcorr{g}-popres.CA1V12Dcorr_randX{g})
        set(gca,'Xlim',[1 100],'Ylim',[1 100],'Clim',[0 2],'Ydir','normal','XAxisLocation','origin','YAxisLocation','origin','PlotBoxAspectRatio', [1 1 1])
        hold on;plot([50 50],[1 100],'k');
        hold on;plot([1 100],[50 50],'k');
        hold on;plot([1 100],[1 100],'k');
        title([titlestr{g} '(X and speed)']);
        %     title([titlestr{g} '(- X)']);
        xlabel('CA1');
        ylabel('V1');
        
        subplot(5,3,9+g);
        imagesc((CA1V12Dcorr{g}-CA1V12Dcorr_randXS{g}));%./(popres.CA1V12Dcorr_randXSsqr{g}-popres.CA1V12Dcorr_randXS{g}.^2).^0.5)
        set(gca,'Xlim',[1 100],'Ylim',[1 100],'Clim',[-0.4 0.4],'Ydir','normal','XAxisLocation','origin','YAxisLocation','origin','PlotBoxAspectRatio', [1 1 1])
        hold on;plot([50 50],[1 100],'k');
        hold on;plot([1 100],[50 50],'k');
        hold on;plot([1 100],[1 100],'k');
        title([titlestr{g} '(- X x Spd)']);
        xlabel('CA1');
        ylabel('V1');
        
        subplot(5,3,12+g);
        plot(-floor(numel(popres.CA1V11Dcorr{g})/2)*1/60:1/60:floor(numel(popres.CA1V11Dcorr{g})/2)*1/60,popres.CA1V11Dcorr{g});
        hold on;plot(-floor(numel(popres.CA1V11Dcorr{g})/2)*1/60:1/60:floor(numel(popres.CA1V11Dcorr{g})/2)*1/60,popres.CA1V11Dcorr_randXS{g});
        hold on;plot(-floor(numel(popres.CA1V11Dcorr{g})/2)*1/60:1/60:floor(numel(popres.CA1V11Dcorr{g})/2)*1/60,popres.CA1V11Dcorr{g} - popres.CA1V11Dcorr_randXS{g});
        hold on;plot([0 0],[0 0.1])
        axis tight;
        set(gca,'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
        title([titlestr{g}]);
        xlabel('time lag (s)');
        ylabel('error correlation');
    end
end

figure('name','Licks');
for g = [1 3]
    for iprobe = 1:nProbe
    for xx = 1:size(meanErrAllSeg{iprobe,g},1)
        subplot(nProbe+1,nSeg,iprobe*nSeg+xx);
        hold on;
        scatter(meanErrAllSeg{iprobe,g}(xx,:)-meanErrAllSeg{iprobe,2}(xx,:),meanLickAll{1,g}-meanLickAll{1,2},cl{g});
        vecErr = meanErrAllSeg{iprobe,g}(xx,:)-meanErrAllSeg{iprobe,2}(xx,:);
        vecLick = meanLickAll{1,g}-meanLickAll{1,2};
%         [rho,p] = corr(vecErr(~isnan(vecErr) & ~isnan(vecLick))',vecLick(~isnan(vecErr) & ~isnan(vecLick))');
%         text(-10,-7 + 3*(g-1),[num2str(rho) '               ' num2str(p)]);
        hold on;plot([-10 10],[10 -10],'k');
    end
    end
    subplot(3,5,1);
    hshift = histogram(meanLickAll{1,g}-meanLickAll{1,2},[-inf -10.5:1:10.5 inf],'EdgeColor','none','FaceColor',cl{g});
    set(gca,'Xlim',[-20 20],'PlotBoxAspectRatio', [1 1 1]);
    hold on;plot([0 0],[0 15],'k');
    xlabel('Spatial shift (cm)');
    ylabel('# of sessions');
end

figure('Name','Optimal delay')
if ~isempty(reslat)
    for iprobe = 1:size(meanErr,3)
        subplot(1,2,iprobe);
        title(probestr{iprobe});
        hold on;
        hlow = histogram([reslat.Latopt{iprobe,1} reslat.Latopt{iprobe,3}],-50:100:1550,'EdgeColor','none','FaceColor',cl{1});
        hhigh = histogram(reslat.Latopt{iprobe,3},-50:100:1550,'EdgeColor','none','FaceColor',cl{3});
        set(gca,'Xlim',[-100 1550],'PlotBoxAspectRatio', [1 1 1]);
    end
end
end


% maxtol = 1;
% nanimal = size(res.PostX,1);
% res = getresSEM(res,maxtol);
% for iprobe = 1:2
%     for g = 1:3
%         res.XpredAll{iprobe,g} = getCircularAverage(res.PostXAll{iprobe,g},maxtol);%getCircularAverage(res.DistriXAll{iprobe,g},maxtol);
%         res.MeanXErrAll{iprobe,g} = nanmean(res.PostErrXAll{iprobe,g},2);
%     end
% end
% 
% % for ianimal = 1:nanimal
% %     figure;
% %     for iprobe = 1:2
% %         for g = 1:3
% %             if sum(res.PostX{ianimal,iprobe,g}(:))>0
% %                 subplot(2,3,(iprobe-1)*3 + g);
% %                 imagesc(res.PostX{ianimal,iprobe,g});%imagesc(res.DistriXAll{iprobe,g});%
% %                 set(gca,'Ydir','normal','Clim',[0.5 1.5],'PlotBoxAspectRatio', [1 1 1])
% %                 colormap(RedWhiteBlue)
% %                 hold on;plot([0 100],[0 100],'w')
% %                 Xpred = getCircularAverage(res.PostX{ianimal,iprobe,g},maxtol);
% %                 hold on;plot(Xpred,'k');
% %             end
% %         end
% %     end
% % end
% c{1} = 'c';
% c{2} = 'k';
% c{3} = 'm';
% titlestr{1} = 'low';
% titlestr{2} = 'med';
% titlestr{3} = 'high';
% probestr{1} = 'CA1';
% probestr{2} = 'V1';
% for iprobe = 1:2
%     figure('name',['Decoding : ' probestr{iprobe}]);
%     for g = 1:3
%         subplot(3,3,g);
%         imagesc(res.PostXAll{iprobe,g});%imagesc(res.DistriXAll{iprobe,g});%
%         set(gca,'Ydir','normal','Clim',[0.5 1.5],'PlotBoxAspectRatio', [1 1 1])
%         colormap(RedWhiteBlue)
%         hold on;plot([0 100],[0 100],'w')
%         hold on;plot(res.XpredAll{iprobe,g},'k');
%         title([probestr{iprobe} ' ' titlestr{g}]);
%         
%         subplot(3,3,4);
%         XErr = res.XpredAll{iprobe,g} - (1:100)';
%         XErr(XErr>50) = XErr(XErr>50) - 100;
%         XErr(XErr<-50) = XErr(XErr<-50) + 100;
%         hold on;ciplot(gca,1:100,XErr,res.XpredAllsem{iprobe,g}(:)',0.2,c{g});%plot(Xpred(:) - (1:100)');
% %         hold on;plot(XErr,c{g});
%         set(gca,'Ylim',[-20 20]);
%         
%         subplot(3,3,5);
%         [~,imax] = max(res.MeanXErrAll{iprobe,2});
%         imax = imax-1;
% %         hold on;ciplot(gca,-50:49,circshift(res.MeanXErrAll{iprobe,g},50-imax),circshift(res.MeanXErrAllsem{iprobe,g},50-imax),0.2,c{g});
%         hold on;plot(-50:49,circshift(res.MeanXErrAll{iprobe,g},50-imax),c{g});
%         set(gca,'Xlim',[-20 20]);
%         [~,imax] = max(circshift(res.MeanXErrAll{iprobe,g},50-imax));
%         hold on;plot([imax-51 imax-51],[0 4],c{g});
%         set(gca,'Xlim',[-20 20]);
%     end
%     shift_low = [];
%     shift_high = [];
%     OptLat_low = [];
%     OptLat_high = [];
%     for ianimal = 1:size(res.MeanXErrMaxPos,1)
%         shift_low = [shift_low res.MeanXErrMaxPos{ianimal,iprobe,1}-res.MeanXErrMaxPos{ianimal,iprobe,2}];
%         shift_high = [shift_high res.MeanXErrMaxPos{ianimal,iprobe,3}-res.MeanXErrMaxPos{ianimal,iprobe,2}];
%         OptLat_low = [OptLat_low res.Latopt{ianimal,iprobe,1}];
%         OptLat_high = [OptLat_high res.Latopt{ianimal,iprobe,3}];
%     end
%     subplot(3,3,6);
%     hold on;
%     hlow = histogram(shift_low,-5:0.5:5,'EdgeColor','none','FaceColor',c{1});
%     hhigh = histogram(shift_high,-5:0.5:5,'EdgeColor','none','FaceColor',c{3});
%     set(gca,'Xlim',[-5 5],'Ylim',[0 18],'PlotBoxAspectRatio', [1 1 1]);
%     hold on;plot([0 0],[0 18],'k');
%     
%     subplot(3,3,8);
%     hold on;
%     hlow = histogram([OptLat_low OptLat_high],-150:100:750,'EdgeColor','none','FaceColor',c{1});
%     hhigh = histogram(OptLat_high,-150:100:7500,'EdgeColor','none','FaceColor',c{3});
%     set(gca,'Xlim',[-200 1000]);
%     
%     subplot(3,3,7);
%     hold on;
%     meanVSrep = nanmean(cellprop.globalXposrep(cellprop.Probe == iprobe & cellprop.Goodcluster ,16:40)+cellprop.globalYposrep(cellprop.Probe == iprobe & cellprop.Goodcluster ,16:40),1);
%     xtimes = (-4:20)*50;
%     stim = zeros(size(meanVSrep));
%     stim(5:end) = sin(xtimes(5:end)/1000*2*pi);
%     plot(xtimes,meanVSrep);
%     plot(xtimes,nanmean(meanVSrep)+(max(meanVSrep)-nanmean(meanVSrep))*stim);
%     set(gca,'Xlim',[xtimes(1) xtimes(end)]);
%     
%     if iprobe == 2
%         V1gainshift = cell(1,3);
%         CA1gainshift = cell(1,3);
%         subplot(3,3,9);
%         for ianimal = 1:size(res.MeanXErrMaxPos,1)
%             CA1V1series = ~isnan(res.MeanXErrMaxPos{ianimal,1,2}) & ~isnan(res.MeanXErrMaxPos{ianimal,2,2});
%             if sum(CA1V1series)>0
%                 for g = [1 3]
%                     hold on;scatter(res.MeanXErrMaxPos{ianimal,1,g}(CA1V1series) - res.MeanXErrMaxPos{ianimal,1,2}(CA1V1series),res.MeanXErrMaxPos{ianimal,2,g}(CA1V1series) - res.MeanXErrMaxPos{ianimal,2,2}(CA1V1series),...
%                         'MarkerEdgeColor',c{g},'MarkerFaceColor',c{g},'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.2);
%                     V1gainshift{g} = [V1gainshift{g} res.MeanXErrMaxPos{ianimal,2,g}(CA1V1series) - res.MeanXErrMaxPos{ianimal,2,2}(CA1V1series)];
%                     CA1gainshift{g} = [CA1gainshift{g} res.MeanXErrMaxPos{ianimal,1,g}(CA1V1series) - res.MeanXErrMaxPos{ianimal,1,2}(CA1V1series)];
%                 end
%             end
%         end
%         CA1nanidx = isnan(CA1gainshift{g});
%         V1nanidx = isnan(V1gainshift{g});
%         CA1gainshift{g} = CA1gainshift{g}(~CA1nanidx & ~V1nanidx);
%         V1gainshift{g} = V1gainshift{g}(~CA1nanidx & ~V1nanidx);
%         for g = [1 3]
%             [r,p] = corr(CA1gainshift{g}',V1gainshift{g}');
%             rho{g} = r; pval{g} = p;
%         end
%         hold on;plot([-5 5],[-5 5],'k');
%         title(['rlow = ' num2str(rho{1}) ' rhigh = ' num2str(rho{3})]);
%         set(gca,'Xlim',[-5 5],'Ylim',[-5 5],'XAxisLocation','origin','YAxisLocation','origin','PlotBoxAspectRatio', [1 1 1]);
%     end
% end
% 
% figure('Name','CA1V1 correlation')
% Xxcorr = ((1:numel(res.CA1V1corrAll{2})) - floor(numel(res.CA1V1corrAll{2})/2)+1)*1000/60;
% for g = 1:3
%     subplot(1,3,g);
%     hold on;plot(Xxcorr,res.CA1V1corrAll{g});
%     hold on;plot(Xxcorr,res.CA1V1corrAll_rand{g});
%     hold on;plot(Xxcorr,res.CA1V1corrAll_randX{g});
%     hold on;plot(Xxcorr,res.CA1V1corrAll_randS{g});
%     set(gca,'Xlim',[Xxcorr(1) Xxcorr(end)],'Ylim', [-0.02 0.07]);
%     title(titlestr{g});
% end


% end

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

function res = getresSEM(res,maxtol)
nanimal = size(res.PostX,1);
for iprobe = 1:2
    for g = 1:3
        for iperm = 0:nanimal
            PostXperm = 0;
            PostErrXperm = 0;
            DistriXperm = 0;
            DistriErrXperm = 0;
            nAll = 0;
            for ianimal = 1:nanimal%[1 2 3 4 5 8 9]%
                if sum(res.PostX{ianimal,iprobe,g}(:))>0 && ianimal ~= iperm
                    nAll = nAll + 1;
                    PostXperm = PostXperm + res.PostX{ianimal,iprobe,g};
                    PostErrXperm = PostErrXperm + res.PostErrX{ianimal,iprobe,g};
                    DistriXperm = DistriXperm + res.DistriX{ianimal,iprobe,g};
                    DistriErrXperm = DistriErrXperm + res.DistriErrX{ianimal,iprobe,g};
                end
            end
            PostXperm = PostXperm/nAll;
            Xpredperm = getCircularAverage(PostXperm,maxtol);
            XErr = Xpredperm - (1:100)';
            XErr(XErr>50) = XErr(XErr>50) - 100;
            XErr(XErr<-50) = XErr(XErr<-50) + 100;
            
            PostErrXperm = PostErrXperm/nAll;
            MeanXErrperm = nanmean(PostErrXperm,2);
            DistriXperm = DistriXperm/nAll;
            DistriErrXperm = DistriErrXperm/nAll;
            
            if iperm == 0
                res.PostXAll{iprobe,g} = PostXperm;%res.PostXAll{iprobe,g} = PostXperm;
                res.PostErrXAll{iprobe,g} = PostErrXperm;
                res.DistriXAll{iprobe,g} = DistriXperm;
                res.DistriErrXAll{iprobe,g} = DistriErrXperm;
                
                res.XpredAll{iprobe,g} = getCircularAverage(res.PostXAll{iprobe,g},maxtol);
                XErrAll{iprobe,g} = XErr;
                res.XpredAll{iprobe,g} = res.XpredAll{iprobe,g};
                res.MeanXErrAll{iprobe,g} = nanmean(res.PostErrXAll{iprobe,g},2);
                
                res.XpredAllsem{iprobe,g} = 0;
                res.MeanXErrAllsem{iprobe,g} = 0;
            else
                res.XpredAllsem{iprobe,g} = res.XpredAllsem{iprobe,g} + (XErr - XErrAll{iprobe,g}).^2;
                res.MeanXErrAllsem{iprobe,g} = res.MeanXErrAllsem{iprobe,g} + (MeanXErrperm - res.MeanXErrAll{iprobe,g}).^2;
            end
        end
        res.XpredAllsem{iprobe,g} = ((res.XpredAllsem{iprobe,g}*(nAll - 1)/nAll).^0.5)/nAll^0.5;
        res.MeanXErrAllsem{iprobe,g} = ((res.MeanXErrAllsem{iprobe,g}*(nAll - 1)/nAll).^0.5)/nAll^0.5;
    end
end


end

function mat_out = smooth2D(mat,lambdaSmooth)
mat(isnan(mat)) = 0;
G = smooth1D(repmat(mat,3,3),lambdaSmooth);
H = smooth1D(G',lambdaSmooth)';
mat_out = H(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
end

function [predave_out,x] = Xpredcorrection(predave,predave_ref)
predave_ref(predave_ref - (1:100)' > 50) = predave_ref(predave_ref - (1:100)' > 50) - 100;
predave_ref(predave_ref - (1:100)' < -50) = predave_ref(predave_ref - (1:100)' < -50) + 100;
predave(predave - (1:100)' > 50) = predave(predave - (1:100)' > 50) - 100;
predave(predave - (1:100)' < -50) = predave(predave - (1:100)' < -50) + 100;

predave_interp = interp1(linspace(0,numel(predave),numel(predave)+1), [predave(1);predave], 0:0.1:99.9);predave_interp(isnan(predave_interp)) = 0;
predaveref_interp = interp1(linspace(0,numel(predave_ref),numel(predave_ref)+1), [predave_ref(1);predave_ref], 0:0.1:99.9);predaveref_interp(isnan(predaveref_interp)) = 0;
predave_out = zeros(size(predave_interp));
x = 0:0.1:99.9;
predaveref_interp = [predaveref_interp(1:end-1)-100 predaveref_interp predaveref_interp(2:end)+100];
xrep = [x(1:end-1)-100 x x(2:end)+100];
for i = 1:numel(x)
    idxmatch = find(abs(predaveref_interp-predave_interp(i)) <= min(abs(predaveref_interp-predave_interp(i))));
    idxmatch = idxmatch(abs(idxmatch-(numel(x)-1) - i) == min(abs(idxmatch-(numel(x)-1) - i)));
    predave_out(i) = xrep(round(idxmatch)) - x(i);
end
% predave_out(predave_out>50) = predave_out(predave_out>50) - 100;
% predave_out(predave_out<-50) = predave_out(predave_out<-50) + 100;
end