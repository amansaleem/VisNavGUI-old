function PopPlot2(popres,popresCorr,reslat)

titlestr{1} = 'low';
titlestr{2} = 'med';
titlestr{3} = 'high';

gainVal = [0.8 1 1.2];
cl{1} = 'c';
cl{2} = 'g';
cl{3} = 'm';

probestr{1} = popres.probestr{1};
probestr{2} = popres.probestr{2};

errname = 'Err';%'Err';%'Bias';%
errnameTheta = 'Err';
errnamediff = 'Err';

maps = popres.PostXAll;
Xdec = popres.PostXpred;
Errdec = popres.([errname 'Xpred']);
Errdecstd = popres.([errname 'Xpred_SE']);
Errdec_rel = popres.([errnamediff 'Xpreddiff']);
Errdec_relstd = popres.([errnamediff 'Xpreddiff_SE']);
Errdec_reldec = popres.([errname 'Xpred']);
Errdec_reldecstd = popres.([errname 'Xpred']);
% Errdec_reg = popres.(['Post' errname 'pred_reg']);
meanErr = popres.(['s_mean' errname 'X_Post']);
meanErrSpd = popres.(['mean' errname 'X_PostSpd']);
meanErrEye = popres.(['mean' errname 'X_PostEye']);
meanErrVis = popres.(['mean' errname 'X_PostVis']);
if isfield(popres, (['mean' errname 'X_PostSpd_SE']))
    meanErrSpd_std = popres.(['mean' errname 'X_PostSpd_SE']);
    meanErrEye_std = popres.(['mean' errname 'X_PostEye_SE']);
    meanErrVis_std = popres.(['mean' errname 'X_PostVis_SE']);
else
    meanErrSpd_std = cell(size(popres.(['mean' errname 'X_PostSpd'])));
    meanErrEye_std = cell(size(popres.(['mean' errname 'X_PostEye'])));
    meanErrVis_std = cell(size(popres.(['mean' errname 'X_PostSpd'])));
end

meanErrSeg = popres.(['s_mean' errname 'X_Post']);%popres.s_meanErrX_aveSeg;
    
nanimal = size(meanErr,1);
nProbe = size(meanErr,3);
ngain = size(popres.PostXAll,2);

meanErrAll = cell(nProbe,ngain);
meanErrAllZscore = cell(nProbe,ngain);
meanErrAnimal = cell(nProbe,ngain);
meanErrAllSpd = cell(nProbe,ngain);
meanErrAllEye = cell(nProbe,ngain);
meanErrAllVis = cell(nProbe,ngain);
meanRunSpdAll = cell(nProbe,ngain);
meanVisSpdAll = cell(nProbe,ngain);
meanLickAll = cell(nProbe,ngain);
meanErrAllSeg = cell(nProbe,ngain);

thetaZscoreAll = cell(nProbe,ngain);

nSeg = 5;

nSpdbins = popres.nSpdbins;
nVisbins = popres.nVisbins;
nEyebins = popres.nEyebins;
nPhsbins = popres.nPhsbins;

meanLickAnimal = cell(nProbe,ngain);
for iprobe = 1:nProbe
    for g = 1:3
        for ianimal = 1:nanimal
            meanErrSession = 0;
            meanLickSession = 0;
            nSession = 0;
            for iseries = 1:size(meanErr,2)
                if ~isempty(meanErr{ianimal,iseries,iprobe,g}) && ~isnan(meanErr{ianimal,iseries,iprobe,g})
                    meanLickAll{iprobe,g} = [meanLickAll{iprobe,g} max([popres.s_meanlickX{ianimal,iseries,1,g} popres.s_meanlickX{ianimal,iseries,size(popres.s_meanlickX,3),g}])];%[meanLickAll{iprobe,g}  max([popres.s_meanlickX{ianimal,iseries,1,g} popres.s_meanlickX{ianimal,iseries,size(popres.s_meanlickX,3),g}])];
%                     if g~=2
                        meanErrAll{iprobe,g} = [meanErrAll{iprobe,g}  meanErr{ianimal,iseries,iprobe,g}];%[meanErrAll{iprobe,g}  getCircularAverage(nanmean(popres.s_ErrXAll{ianimal,iseries,iprobe,g},2),1,1,0.05)];%
%                     else
%                         meanErrAll{iprobe,g} = [meanErrAll{iprobe,g}  0];%
%                     end
%                     if g~=2
%                         meanErrAll{iprobe,g} = [meanErrAll{iprobe,g}  100/(2*pi)*circ_mean(2*pi/100*(popres.s_ErrXpred{ianimal,iseries,iprobe,g}-50))];%meanErr{ianimal,iseries,iprobe,g}];
%                     else
%                         meanErrAll{iprobe,g} = [meanErrAll{iprobe,g}  0];
%                     end

                    meanErrAllSpd{iprobe,g} = [meanErrAllSpd{iprobe,g}  popres.s_meanErrX_PostSpd{ianimal,iseries,iprobe,g}(:)-nanmean(popres.s_meanErrX_PostSpd{ianimal,iseries,iprobe,g}(:))];
                    meanErrAllEye{iprobe,g} = [meanErrAllEye{iprobe,g}  popres.s_meanErrX_PostEye{ianimal,iseries,iprobe,g}(:)-nanmean(popres.s_meanErrX_PostEye{ianimal,iseries,iprobe,g}(:))];
                    meanErrAllVis{iprobe,g} = [meanErrAllVis{iprobe,g}  popres.s_meanErrX_PostVis{ianimal,iseries,iprobe,g}(:)-nanmean(popres.s_meanErrX_PostVis{ianimal,iseries,iprobe,g}(:))];
                    meanErrAllSeg{iprobe,g} = [meanErrAllSeg{iprobe,g}  meanErrSeg{ianimal,iseries,iprobe,g}'];
                    [speedpro,~,~] = fast1Dmap(popres.s_X{ianimal,iseries,iprobe,g},popres.s_runspeed{ianimal,iseries,iprobe,g},1,1,[], true);
                    meanRunSpdAll{iprobe,g} = [meanRunSpdAll{iprobe,g}  nanmean(speedpro)];
                    [speedpro,~,~] = fast1Dmap(popres.s_X{ianimal,iseries,iprobe,g},popres.s_visspeed{ianimal,iseries,iprobe,g},1,1,[], true);
                    meanVisSpdAll{iprobe,g} = [meanVisSpdAll{iprobe,g}  nanmean(speedpro)];
                    if ~isnan(meanErr{ianimal,iseries,iprobe,g})
                        nSession = nSession + 1;
                        meanErrSession = meanErrSession + meanErr{ianimal,iseries,iprobe,g};
                        meanLickSession = meanLickSession + popres.s_meanlickX{ianimal,iseries,iprobe,g};
                    end
                    if isfield(popres,'ErrXpredPhsAllZscore')
                        thetaZscoreAll{iprobe,g} = [thetaZscoreAll{iprobe,g} popres.ErrXpredNormPhsAllZscore{ianimal,iseries,iprobe,g}];
                    end
                    meanErrAllZscore{iprobe,g} = [meanErrAllZscore{iprobe,g}  popres.meanErrXAllZscore{ianimal,iseries,iprobe,g}];
                else
                    meanLickAll{iprobe,g} = [meanLickAll{iprobe,g}  NaN];
                    
                    meanErrAll{iprobe,g} = [meanErrAll{iprobe,g}  NaN];
                    meanErrAllZscore{iprobe,g} = [meanErrAllZscore{iprobe,g}  NaN];
                    meanErrAllSpd{iprobe,g} = [meanErrAllSpd{iprobe,g}  NaN(nSpdbins,1)];
                    meanErrAllEye{iprobe,g} = [meanErrAllEye{iprobe,g}  NaN(nEyebins,1)];
                    meanErrAllVis{iprobe,g} = [meanErrAllVis{iprobe,g}  NaN(nVisbins,1)];
                    
                    meanRunSpdAll{iprobe,g} = [meanRunSpdAll{iprobe,g} NaN];
                    meanVisSpdAll{iprobe,g} = [meanVisSpdAll{iprobe,g} NaN];
                    meanErrAllSeg{iprobe,g} = [meanErrAllSeg{iprobe,g} NaN];
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
Clim = [0 2];%[0.5 1.5];%[-0.5 0.5];%[0.5 1.5];
for iprobe = 1:nProbe
    figure('name',['Decoding : ' probestr{iprobe}]);
    for g = [2 1 3]
        subplot(3,3,g);
        disp([num2str(nansum(abs(meanErrAllZscore{iprobe,g})>=2)./sum(~isnan(meanErrAllZscore{iprobe,g}))) '  ' num2str(nanmean(meanErrAll{iprobe,g}-meanErrAll{iprobe,2}))]);
%         errmap = popres.PostErrAll{iprobe,g};
%         errmapref = popres.PostErrAll{iprobe,2};
%         for i=1:100
%             errmap(:,i) = errmap(:,i)-mean(errmap(:,i));
%             errmap(:,i) = errmap(:,i)/(mean(errmap(:,i).^2))^0.5;
%             errmapref(:,i) = errmapref(:,i)-mean(errmapref(:,i));
%             errmapref(:,i) = errmapref(:,i)/(mean(errmapref(:,i).^2))^0.5;
%         end
%         errCorrmap = zeros(size(errmap));
%         ishift = 0;
%         for xshift = -49:50
%             ishift = ishift+1;
%             errCorrmap(ishift,:) = diag(errmap'*circshift(errmapref,xshift,1)/size(errmap,1));
%         end
%         Errdec{iprobe,g} = getCircularAverage(errCorrmap,0,1);
%         imagesc(errCorrmap);
%         
        imagesc(X,Y,(maps{iprobe,g}));
%         hold on;plot(X,Xdec{iprobe,2},'w');
        hold on;plot([0 100],[0 100],'w--');
%         hold on;plot(X,Xdec{iprobe,g}*popres.dx,cl{g});
%         hold on;plot(Xdec{iprobe,g},'k');
        colormap(jet);%colormap(parula);%colormap(RedWhiteBlue);%
        set(gca,'Clim',Clim,'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
        xlabel('actual');
        ylabel('decoded');
        title([probestr{iprobe} ' ' titlestr{g}]);
        
        subplot(3,4,5);
        ax = gca;
        if g~=2 || ~strcmp(errname,'Corr')
            hold on;ciplot(ax,X,Errdec{iprobe,g} - 50,Errdecstd{iprobe,g},0.4,cl{g});
            hold on;plot([X(1) X(end)],[0 0],'k--');
        end
%         hold on;plot(X,Errdec{iprobe,g}*popres.dx-Errdec{iprobe,2}*popres.dx,cl{g});
        set(gca,'Xlim',[0 max(X)],'Ylim',[-10 10]);
%         hold on;plot(X/gainVal(g),Errdec{iprobe,g}*popres.dx + X'*(1-1/gainVal(g)),cl{g});
%         set(gca,'Xlim',[X(1)/gainVal(1) X(end)/gainVal(1)],'Ylim',[-20 20]);
        xlabel('position');
        ylabel([probestr{iprobe} ' decoding error']);
        
        if ismember(g,[1 3])
            subplot(3,4,6);
            hold on;ciplot(gca,X,Errdec_rel{iprobe,g},Errdec_relstd{iprobe,g},0.4,cl{g});%plot(X,Errdec{iprobe,g} - Errdec{iprobe,2},cl{g});
            hold on;plot([X(1) X(end)],[0 0],'k--');
            hold on;plot([X(1) X(end)],[(X(1)-1) (X(end)-1)*(1 -gainVal(g))],cl{g},'LineStyle','--');
            %         hold on;plot(Xreg,Errdec_reg{iprobe,g}-Errdec_reg{iprobe,2},cl{g});
            %         hold on;plot([min(Xreg) max(Xreg)],[min(Xreg) max(Xreg)*(1 -gainVal(g))],cl{g},'LineStyle','--');
            set(gca,'Xlim',[0 max(X)],'Ylim',[-10 10]);
            %         hold on;plot(Xreg/gainVal(g),Errdec_reg{iprobe,g}*popres.dx+ Xreg*(1-1/gainVal(g)),cl{g});
            %         hold on;plot([Xreg(1)/gainVal(g) Xreg(end)/gainVal(g)],[0 -(1 -gainVal(g))*Xreg(end)/gainVal(g)],cl{g},'LineStyle','--');
            %         set(gca,'Xlim',[Xreg(1)/gainVal(1) Xreg(end)/gainVal(1)],'Ylim',[-20 20]);
            xlabel('position');
            ylabel([probestr{iprobe} ' relative error']);
        end
        
        if ismember(g,[1 3])
            subplot(3,4,7);
            hold on;
            hshift = histogram(meanErrAll{iprobe,g}-meanErrAll{iprobe,2},[-inf -10.5:1:10.25 inf],'EdgeColor','none','FaceColor',cl{g});
            set(gca,'Xlim',[-11 11]);
            hold on;plot([0 0],[0 15],'k');
            xlabel([probestr{iprobe} ' decoding error']);
            ylabel('# of sessions');
            
            subplot(3,4,8);
            hold on;
            scatter(meanErrAll{iprobe,1}-meanErrAll{iprobe,2},meanErrAll{iprobe,3}-meanErrAll{iprobe,2},'k.');
            scatter(meanErrAll{iprobe,1}(abs(meanErrAllZscore{iprobe,1})>=2)-meanErrAll{iprobe,2}(abs(meanErrAllZscore{iprobe,1})>=2),meanErrAll{iprobe,3}(abs(meanErrAllZscore{iprobe,1})>=2)-meanErrAll{iprobe,2}(abs(meanErrAllZscore{iprobe,1})>=2),'co');
            scatter(meanErrAll{iprobe,1}(abs(meanErrAllZscore{iprobe,3})>=2)-meanErrAll{iprobe,2}(abs(meanErrAllZscore{iprobe,3})>=2),meanErrAll{iprobe,3}(abs(meanErrAllZscore{iprobe,3})>=2)-meanErrAll{iprobe,2}(abs(meanErrAllZscore{iprobe,3})>=2),'ms');
            set(gca,'Xlim',[-10 10]);
            set(gca,'Ylim',[-10 10]);
            hold on;plot([-10 10],[10 -10],'k');
            hold on;plot([-10 10],[0 0],'k');
            hold on;plot([0 0],[-10 10],'k');
            ylabel([probestr{iprobe} ' relative error high']);
            xlabel([probestr{iprobe} ' relative error low']);
        end
        
        if iprobe == 1
            subplot(3,3,7);
            [speedpro,x,speedprostd] = fast1Dmap(popres.X{1,g},popres.runspeed{1,g},1,1,[], true);
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
            hold on;scatter(meanRunSpdAll{iprobe,g},meanErrAll{iprobe,g}-meanErrAll{iprobe,2},cl{g});
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

if ~isempty(popresCorr)
    figure('Name','V1CA1 noise correlation');
    for g = [2 1 3]
        subplot(5,3,g);
        imagesc(popresCorr.CA1V1Cov{g})
        set(gca,'Xlim',[1 100],'Clim',[-0.4 0.4],'Ylim',[1 100],'Ydir','normal','XAxisLocation','origin','YAxisLocation','origin','PlotBoxAspectRatio', [1 1 1])
        hold on;plot([50 50],[1 100],'k');
        hold on;plot([1 100],[50 50],'k');
        hold on;plot([1 100],[1 100],'k');
        title([titlestr{g} 'Cov']);
        %     title([titlestr{g} '(- spd)']);
        xlabel('CA1');
        ylabel('V1');
        
        subplot(5,3,3+g);
        imagesc(popresCorr.CA1V1Cov_shuffled{g})
        set(gca,'Xlim',[1 100],'Clim',[-0.4 0.4],'Ylim',[1 100],'Ydir','normal','XAxisLocation','origin','YAxisLocation','origin','PlotBoxAspectRatio', [1 1 1])
        hold on;plot([50 50],[1 100],'k');
        hold on;plot([1 100],[50 50],'k');
        hold on;plot([1 100],[1 100],'k');
        title([titlestr{g} 'Cov_shf']);
        %     title([titlestr{g} '(- spd)']);
        xlabel('CA1');
        ylabel('V1');
        
        subplot(5,3,6+g);
        imagesc(popresCorr.CA1V1Cross{g})
        set(gca,'Xlim',[1 100],'Clim',[-0.05 0.05],'Ylim',[1 100],'Ydir','normal','XAxisLocation','origin','YAxisLocation','origin','PlotBoxAspectRatio', [1 1 1])
        hold on;plot([50 50],[1 100],'k');
        hold on;plot([1 100],[50 50],'k');
        hold on;plot([1 100],[1 100],'k');
        title([titlestr{g} 'Cross']);
        %     title([titlestr{g} '(- spd)']);
        xlabel('CA1');
        ylabel('V1');
        
        subplot(5,3,9+g);
        imagesc(popresCorr.CA1V1Cross_shuffled{g})
        set(gca,'Xlim',[1 100],'Clim',[-0.05 0.05],'Ylim',[1 100],'Ydir','normal','XAxisLocation','origin','YAxisLocation','origin','PlotBoxAspectRatio', [1 1 1])
        hold on;plot([50 50],[1 100],'k');
        hold on;plot([1 100],[50 50],'k');
        hold on;plot([1 100],[1 100],'k');
        title([titlestr{g} 'Cross_shf']);
        %     title([titlestr{g} '(- spd)']);
        xlabel('CA1');
        ylabel('V1');
        
        if isfield(popresCorr,'CA1V1Coh')
            subplot(5,3,12+g);
            plot(nanmean(popresCorr.CA1V1Coh{g}));
            hold on;
            plot(popresCorr.CA1V1Coh_rand{g});
        end
    end
end

figure('name','Licks');
for g = [1 3]
    for iprobe = 1:nProbe
    for xx = 1:size(meanErrAllSeg{iprobe,g},1)
        subplot(nProbe+1,nSeg,iprobe*nSeg+xx);
        hold on;
        scatter(meanErrAll{iprobe,g}-meanErrAll{iprobe,2},meanLickAll{1,g}-meanLickAll{1,2},cl{g});
        vecErr = meanErrAll{iprobe,g}-meanErrAll{iprobe,2};
        vecLick = meanLickAll{1,g}-meanLickAll{1,2};
        [rho,p] = corr(vecErr(~isnan(vecErr) & ~isnan(vecLick))',vecLick(~isnan(vecErr) & ~isnan(vecLick))');
        text(-10,-7 + 3*(g-1),[num2str(rho) '               ' num2str(p)]);
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

if isfield(popres,['s_' errnameTheta 'XPhsAll'])
    figure('Name','Theta precession')
    for iprobe = 1:nProbe
        for g = [2 1 3]
            subplot(2,4,(iprobe-1)*4+1);
            if g == 2
                imagesc(popres.([errnameTheta 'XPhsAll']){iprobe,g});%./repmat(max(popres.PostErrPhsAll{iprobe,g},[],1),[size(popres.PostErrPhsAll{iprobe,g},1) 1]));
            end
            set(gca,'Clim',[0 1],'Ylim',[30 70],'Ydir','normal');
            subplot(2,4,(iprobe-1)*4+2);
            hold on;
            if g == 2
                ciplot(gca,linspace(0,360,size(popres.([errnameTheta 'XpredPhsAll']){iprobe,g},1)),popres.([errnameTheta 'XpredPhsAll']){iprobe,g}-mean(popres.([errnameTheta 'XpredPhsAll']){iprobe,g}),popres.([errnameTheta 'XpredPhsAll_SE']){iprobe,g},0.5,cl{g});
            else
                plot(linspace(0,360,size(popres.([errnameTheta 'XpredPhsAll']){iprobe,g},1)),popres.([errnameTheta 'XpredPhsAll']){iprobe,g}-mean(popres.([errnameTheta 'XpredPhsAll']){iprobe,g}),cl{g});
            end
            subplot(2,4,(iprobe-1)*4+3);
            hold on;
            if g == 2
                ciplot(gca,linspace(0,360,size(popres.([errnameTheta 'XpredNormPhsAll']){iprobe,g},1)),popres.([errnameTheta 'XpredNormPhsAll']){iprobe,g},popres.([errnameTheta 'XpredNormPhsAll_SE']){iprobe,g},0.5,cl{g});
            else
                plot(gca,linspace(0,360,size(popres.([errnameTheta 'XpredNormPhsAll']){iprobe,g},1)),popres.([errnameTheta 'XpredNormPhsAll']){iprobe,g},cl{g});
            end
            subplot(2,4,(iprobe-1)*4+4);
            histogram(thetaZscoreAll{iprobe,2},[0:0.5:10 +inf],'EdgeColor','none','FaceColor',cl{g});
%             set(gca,'Ylim',[45 55]);
        end
    end
end

figure('Name','Controls')
for iprobe = 1:nProbe
    for g = [1 3]
        subplot(4,4,(iprobe-1)*2+1);
        hold on;
        bar(meanErrSpd{iprobe,g}-meanErrSpd{iprobe,2},'EdgeColor','none','FaceColor',cl{g});
        errorbar(1:numel(meanErrSpd{iprobe,g}),meanErrSpd{iprobe,g}-meanErrSpd{iprobe,2},meanErrSpd_std{iprobe,g},'linestyle','none','color','k');
%         plot(nanmean(meanErrAllSpd{iprobe,g}-meanErrAllSpd{iprobe,2},2),cl{g});
        set(gca,'Xlim',[0 nSpdbins+1]);
        set(gca,'Ylim',[-10 10]);
        ylabel([probestr{iprobe} ' decoding error']);
        xlabel('Run speed bins');
        
        subplot(4,4,(iprobe-1)*2+2);
        hold on;scatter(meanRunSpdAll{iprobe,g},meanErrAll{iprobe,g}-meanErrAll{iprobe,2},cl{g});
        hold on;plot([1 70],[0 0]);
        set(gca,'Xlim',[1 70],'Ylim',[-10 10],'PlotBoxAspectRatio', [1 1 1]);
        xlabel('mean run speed');
        
        subplot(4,4,4+(iprobe-1)*2+1);
        hold on;
        bar(meanErrVis{iprobe,g} - meanErrVis{iprobe,2},'EdgeColor','none','FaceColor',cl{g});
        errorbar(1:numel(meanErrVis{iprobe,g}),meanErrVis{iprobe,g}-meanErrVis{iprobe,2},meanErrVis_std{iprobe,g},'linestyle','none','color','k');
%         plot(nanmean(meanErrAllVis{iprobe,g}-meanErrAllVis{iprobe,2},2),cl{g});
        set(gca,'Xlim',[0 nVisbins+1]);
        set(gca,'Ylim',[-10 10]);
        ylabel([probestr{iprobe} ' decoding error']);
        xlabel('Visual speed bins');
        
        subplot(4,4,4+(iprobe-1)*2+2);
        hold on;scatter(meanVisSpdAll{iprobe,g},meanErrAll{iprobe,g}-meanErrAll{iprobe,2},cl{g});
        hold on;plot([1 70],[0 0]);
        set(gca,'Xlim',[1 70],'Ylim',[-10 10],'PlotBoxAspectRatio', [1 1 1]);
        xlabel('mean visual speed');
        
        subplot(4,4,8+(iprobe-1)*2+1);
        hold on;
        bar(meanErrEye{iprobe,g}-meanErrEye{iprobe,2},'EdgeColor','none','FaceColor',cl{g});
        errorbar(1:numel(meanErrEye{iprobe,g}),meanErrEye{iprobe,g}-meanErrEye{iprobe,2},meanErrEye_std{iprobe,g},'linestyle','none','color','k');
%         plot(nanmean(meanErrAllEye{iprobe,g}-meanErrAllEye{iprobe,2},2),cl{g});
        set(gca,'Xlim',[0 nEyebins+1],'Ylim',[-10 10]);
        ylabel([probestr{iprobe} ' decoding error']);
        xlabel('Eye Xpos bins');
    end
end
if ~isempty(reslat)
    for iprobe = 1:size(meanErr,3)
        subplot(4,4,12 + (iprobe-1)*2 + 1);
        title(probestr{iprobe});
        hold on;
        hlow = histogram([reslat.Latopt{iprobe,1} reslat.Latopt{iprobe,3}],-50:100:1550,'EdgeColor','none','FaceColor',cl{1});
        hhigh = histogram(reslat.Latopt{iprobe,3},-50:100:1550,'EdgeColor','none','FaceColor',cl{3});
        set(gca,'Xlim',[-100 1550]);
    end
end
end