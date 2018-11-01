function PopPlotcellprop(cellprop)
batch2p = true;
strlistvarname = {'2p data','electrophys data'};
[varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
if ok && varnamesel == 1
    batch2p = true;
    strlistvarname = {'V1medial','V1lateral','PPC', 'AL'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'area','InitialValue', 1);
    Cellpos2p = varnamesel;
    probestr{1} = cell2mat(strlistvarname(varnamesel));
    probestr{2} = 'none';
    nProbe = 1;
elseif ok
    batch2p = false;
    Cellpos2p = 0;
    if ~isfield(cellprop, 'Cellpos2p')
        cellprop.Cellpos2p = zeros(size(cellprop.Probe));
    end
    probestr{1} = 'CA1';
    probestr{2} = 'V1';
    nProbe = 2;
end


c{1} = 'c';
c{2} = 'k';
c{3} = 'm';
titlestr{1} = 'low';
titlestr{2} = 'med';
titlestr{3} = 'high';

field_orig = cellprop.field;

SSImin = 0;%0.75;%
Xcorrmin = 0.7;%-inf;%0.75;%-inf;%
maxshift = +inf;%30;%30;%30;%+inf;%
Reliabilitymin = -inf;%0.7;%
thetaReliabilitymin = -inf;%0.5;%
zth = 2;%3;%
zth_lowhigh = 2;%-inf;%2;

for iprobe = 1:nProbe%1:nProbe
    if iprobe == 1
        chfocus = 1:32;
        xposfocus = 1:10;
        if batch2p
            maxminrate = +inf;%5;%
            minmaxrate = -inf;%3;%
            cellprop.bestchan = ones(size(cellprop.bestchan));
        else
            maxminrate = +inf;%8;%+inf;%10;%5;%
            minmaxrate = 0;%0;%3;%
        end
    else
        chfocus = 1:32;
        xposfocus = 1:10;
        if batch2p
            maxminrate = +inf;%5;%
            minmaxrate = -inf;%3;%cellprop.bestchan = ones(size(cellprop.bestchan));
        else
            maxminrate = inf;%+inf;%10;%5;%
            minmaxrate = 0;%3;%
        end
    end
    figure('name',['Single cells : ' probestr{iprobe}]);
    iplot = 0;
    
    cellprop.field = field_orig;
    
    shift_low0 = (cellprop.fieldXcorrMax{1}-cellprop.fieldXcorrMax{2})';
    shift_high0 = (cellprop.fieldXcorrMax{3}-cellprop.fieldXcorrMax{2})';
    for g = 1:3
    fieldZ{g} = zeros(1,size(cellprop.fieldSE{2},1));
    for icell = 1:size(cellprop.fieldSE{2},1)
        fieldZ{g}(icell) = (cellprop.field{g}(icell,min(floor(cellprop.fieldPos{g}(icell))+1,size(cellprop.field{g},2)))-mean(cellprop.field{g}(icell,:)))/cellprop.fieldSE{g}(icell,min(floor(cellprop.fieldPos{g}(icell))+1,size(cellprop.field{g},2)));
    end
    end
    numBinsX = size(cellprop.fieldXcorr{2},2);
    goodcells = (cellprop.Goodcluster & ~cellprop.Finterneuron & cellprop.Probe == iprobe & ismember(cellprop.bestchan,chfocus)... % & ismember(cellprop.XposPop,xposfocus)...
                     & cellprop.SSI{2} >=SSImin & cellprop.fieldMin{2} <= maxminrate & cellprop.fieldMax{2} >= minmaxrate...
                     & max(cellprop.fieldXcorr{1}(:,floor(numBinsX/4):floor(3*numBinsX/4)),[],2)'>Xcorrmin  & max(cellprop.fieldXcorr{3}(:,floor(numBinsX/4):floor(3*numBinsX/4)),[],2)'>Xcorrmin...
                     & cellprop.reliabilityCorr{2}' >= Reliabilitymin...
                     & abs(shift_low0')<maxshift & abs(shift_high0')<maxshift & ismember(cellprop.Cellpos2p,Cellpos2p) & fieldZ{2} >= zth & fieldZ{1} >= zth_lowhigh  & fieldZ{3} >= zth_lowhigh);
                 
    maxpos = cellprop.fieldPos{2};
%     maxpos = [];
%     for icell = 1:size(cellprop.field_set2{2},1)
%         xorig = 1:size(cellprop.field_set2{2},2);
%         xinterp = 0:0.1:size(cellprop.field_set2{2},2);
%         mapinterp = interp1(xorig,cellprop.field_set2{2}(icell,:),xinterp,'spline');
%         [~,imax] = max(mapinterp);
%         maxpos(icell) = xinterp(imax);
%     end
%     cellprop.field = cellprop.field_set1;
    
    [fPos,isort] = sort(maxpos,'ascend');
    
    for g = [2 1 3]
        iplot = iplot + 1;
        subplot(3,9,[(g-1)*3+1 (g-1)*3+2 (g-1)*3+1+9 (g-1)*3+2+9]);
        normfields = (cellprop.field{g}(isort(goodcells(isort)),:) - repmat(min(cellprop.field{g}(isort(goodcells(isort)),:),[],2),[1 size(cellprop.field{g},2)]))./repmat(max(cellprop.field{g}(isort(goodcells(isort)),:),[],2)-min(cellprop.field{g}(isort(goodcells(isort)),:),[],2),[1 size(cellprop.field{g},2)]);
%         normfields = (cellprop.field{g}(isort(goodcells(isort)),:))./repmat(max(cellprop.field{g}(isort(goodcells(isort)),:),[],2),[1 size(cellprop.field{g},2)]);
        centeredfields = zeros(size(normfields));
        for icell = 1:size(normfields,1)
            centeredfields(icell,:) = circshift(normfields(icell,:),-round(fPos(icell))+50);
        end        
        

        imagesc(-normfields);
        hold on;plot(maxpos(isort(goodcells(isort))),1:sum(goodcells),'w');
%         imagesc(-centeredfields);
%         hold on;plot([50 50],[1 numel(goodcells)],'w');
        
        colormap(bone);%colormap(parula);%(RedWhiteBlue);
        
%         hold on;plot(mod(fPos+50,100),1:numel(goodcells),'r');
        set(gca,'Ydir','normal','Clim', [-1 0],'PlotBoxAspectRatio', [2 3 1]);
        title([probestr{iprobe} ' ' titlestr{g}]);
        xlabel('position');
        ylabel('cell#');
        
%         normfields = (cellprop.field{g}(goodcells(isort),:))./repmat(max(cellprop.field{g}(goodcells(isort),:),[],2),[1 size(cellprop.field{g},2)]);
        
        meanfield = nanmean(centeredfields);
        
        subplot(3,9,[(g-1)*3+3 (g-1)*3+3+9]);
        nbin = 10;
        vec{iprobe,g} = zeros(nbin,100);
        cellave = cell(1,nbin);
        for xx = 1:nbin
            cellave{xx} = find(ismember(floor(maxpos(isort(goodcells(isort)))),mod(((xx-1)*(100/nbin)-(100/nbin)/2):((xx-1)*(100/nbin)+(100/nbin)/2),100)));
            vec{iprobe,g}(xx,:) = mean(normfields(ismember(floor(maxpos(isort(goodcells(isort)))),mod(((xx-1)*(100/nbin)-(100/nbin)/2):((xx-1)*(100/nbin)+(100/nbin)/2),100)),:));
            hold on;
            if g == 2 && iprobe == 1
                plot(vec{iprobe,g}(xx,:)+xx,'k')
            elseif g == 2 && iprobe == 2
                plot(vec{iprobe,g}(xx,:)+xx,'k')
                plot(vec{1,g}(xx,:)+xx,'--','Color',[0.5 0.5 0.5])
            else
                plot(vec{iprobe,g}(xx,:)+xx,c{g})
                plot(vec{iprobe,2}(xx,:)+xx,c{2})
            end
        end
        set(gca,'Xlim', [1 100],'Ylim', [1 nbin+1],'PlotBoxAspectRatio', [1 3.5 1]);
        
        subplot(3,9,[(g-1)*3+1 (g-1)*3+2 (g-1)*3+1+9 (g-1)*3+2+9]);
        for xx = 1:nbin
            hold on;plot(100*ones(1,numel(cellave{xx})),cellave{xx},'linewidth',2);
        end
        
%         disp(num2str(max(meanfield(1),meanfield(end))/meanfield(50)));
        
%         subplot(3,3,9);
%         hold on;plot(-50:49,mean(centeredfields,1),c{g});
%         [~,imax] = max(mean(centeredfields,1));
%         hold on;plot([imax-51 imax-51],[0 1],c{g});
%         set(gca,'Xlim',[-20 20],'Ylim',[0.3 1]);
    end    
    
%     [~,shift_low] = max(cellprop.fieldXcorr{1}(goodcells,:),[],2);
%     shift_low = shift_low - (floor(size(cellprop.fieldXcorr{1},2)/2) + 1);
    shift_low = shift_low0(goodcells);%cellprop.fieldXcorrMax{1}(goodcells)-cellprop.fieldXcorrMax{2}(goodcells);
    shift_low_SE = cellprop.fieldXcorrMaxSE{1}(goodcells);
%     [~,shift_high] = max(cellprop.fieldXcorr{3}(goodcells,:),[],2);
%     shift_high = shift_high - (floor(size(cellprop.fieldXcorr{3},2)/2) + 1);
    shift_high = shift_high0(goodcells);%cellprop.fieldXcorrMax{3}(goodcells)-cellprop.fieldXcorrMax{2}(goodcells);
    shift_high_SE = cellprop.fieldXcorrMaxSE{3}(goodcells);
    
    disp(num2str(sum(shift_low(:)./shift_low_SE(:)<-zth)/numel(shift_low)))
    disp(num2str(sum(shift_high(:)./shift_high_SE(:)>zth)/numel(shift_high)))
    disp(num2str(sum(shift_low(:)./shift_low_SE(:)<-zth | shift_high(:)./shift_high_SE(:)>zth)/numel(shift_low)))
%     for k = 1:floor(numel(goodcells)/20)+1
%         figure;
%         n = 0;
%         for icell = (1+(k-1)*20):k*20
%             if icell <= numel(goodcells)
%                 n = n+1;
%                 subplot(4,5,n);
%                 for g = 1:3
%                     hold on;plot(cellprop.field{g}(goodcells(icell),:),c{g})
%                     axis tight;
%                     text(1,0.9,[num2str(cellprop.fieldPos{2}(goodcells(icell))) '   ' num2str(shift_low(icell)) '   ' num2str(shift_high(icell))])
%                 end
%             end
%         end
%     end

%     shift_low = cellprop.fieldCOM{1}(goodcells) - cellprop.fieldCOM{2}(goodcells);%imax_low;%
%     shift_low(shift_low > 50) = shift_low(shift_low > 50) - 100;
%     shift_low(shift_low < -50) = shift_low(shift_low < -50) + 100;
%     shift_high = cellprop.fieldCOM{3}(goodcells) - cellprop.fieldCOM{2}(goodcells);%imax_high;%
%     shift_high(shift_high > 50) = shift_high(shift_high > 50) - 100;
%     shift_high(shift_high < -50) = shift_high(shift_high < -50) + 100;
    fieldpos = cellprop.fieldPos{2}(goodcells);
    subplot(3,6,13);
    scatter(fieldpos(shift_low>=0),shift_low(shift_low>=0),'MarkerEdgeColor',c{1},'MarkerFaceColor',c{1},'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',1,'SizeData',10,'Marker','.');
    hold on;scatter(fieldpos(shift_high<=0),shift_high(shift_high<=0),'MarkerEdgeColor',c{3},'MarkerFaceColor',c{3},'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',1,'SizeData',10,'Marker','.');
    
    hold on;scatter(fieldpos(shift_low<0),shift_low(shift_low<0),'MarkerEdgeColor',c{1},'MarkerFaceColor',c{1},'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',1,'SizeData',10,'Marker','.');
    hold on;scatter(fieldpos(shift_low(:)./shift_low_SE(:)<-zth),shift_low(shift_low(:)./shift_low_SE(:)<-zth),'MarkerEdgeColor','k','MarkerEdgeAlpha',0.3,'MarkerFaceColor','none','SizeData',20);
    hold on;plot([0 100],[0 0],'k');
    set(gca,'Ylim',[-20 20]);
    xlabel('position');
    ylabel('response shift');
%     subplot(3,6,15);
    scatter(fieldpos(shift_high>0),shift_high(shift_high>0),'MarkerEdgeColor',c{3},'MarkerFaceColor',c{3},'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',1,'SizeData',10,'Marker','.');
    hold on;scatter(fieldpos(shift_high(:)./shift_high_SE(:)>zth),shift_high(shift_high(:)./shift_high_SE(:)>zth),'MarkerEdgeColor','k','MarkerEdgeAlpha',0.3,'MarkerFaceColor','none','SizeData',20);
    hold on;plot([0 100],[0 0],'k');
    set(gca,'Ylim',[-20 20]);
    xlabel('position');
    ylabel('response shift');
    
    
    subplot(3,6,14);
    hlow = histogram(shift_low,-50.5:50.5);
    Vlow = hlow.Values;
    hhigh = histogram(shift_high,-50.5:50.5);
    Vhigh = hhigh.Values;
    xbins = hhigh.BinEdges(1:end-1)+(hhigh.BinEdges(2)-hhigh.BinEdges(1))/2;
    
    subplot(3,6,14);
    b = barh(xbins,Vlow);
    b(1).FaceColor = c{1};
    b(1).FaceAlpha = 0.5;
    b(1).LineStyle = 'none';
    hold on;plot([0 100],[0 0],'k');
    set(gca,'Ylim',[-20 20],'Xlim',[0 max(max(Vlow),max(Vhigh))]);
    xlabel('# of cells');
    ylabel('response shift');
%     subplot(3,6,16);
    b = barh(xbins,Vhigh);
    b(1).FaceColor = c{3};
    b(1).FaceAlpha = 0.5;
    b(1).LineStyle = 'none';
    hold on;plot([0 100],[0 0],'k');
    set(gca,'Ylim',[-20 20],'Xlim',[0 max(max(Vlow),max(Vhigh))]);
    xlabel('# of cells');
    ylabel('response shift');
    
    
    subplot(3,3,9);
    animal_list = unique(cellprop.animal);
    meanshift_low = zeros(1,numel(animal_list));
    meanshift_high = zeros(1,numel(animal_list));
    semshift_low = zeros(1,numel(animal_list));
    semshift_high = zeros(1,numel(animal_list));
    for ianimal = 1:numel(animal_list)
        meanshift_low(ianimal) = mean(shift_low(cellprop.animal(goodcells) == animal_list(ianimal)));
        meanshift_high(ianimal) = mean(shift_high(cellprop.animal(goodcells) == animal_list(ianimal)));
        semshift_low(ianimal) = sem(shift_low(cellprop.animal(goodcells) == animal_list(ianimal)));
        semshift_high(ianimal) = sem(shift_high(cellprop.animal(goodcells) == animal_list(ianimal)));
    end
    b = bar(1:numel(animal_list),meanshift_low);
    b(1).FaceColor = c{1};
    b(1).LineStyle = 'none';
    set(gca,'Ylim',[-20 20],'Xlim',[0 numel(animal_list)+1]);
    hold on;
    errorbar(1:numel(animal_list),meanshift_low,semshift_low,'Color',c{1},'LineStyle','none');
    b = bar(1:numel(animal_list),meanshift_high);
    b(1).FaceColor = c{3};
    b(1).LineStyle = 'none';
    errorbar(1:numel(animal_list),meanshift_high,semshift_high,'Color',c{3},'LineStyle','none');
    set(gca,'Ylim',[-20 20],'Xlim',[0 numel(animal_list)+1]);
    xlabel('animal #');
    ylabel('mean response shift');
    
%     %to plot mean peak firing across cells
%     subplot(3,3,9);
%     for g = 1:3
%         meanrate(g) = mean(cellprop.fieldMax{g}(goodcells));
%         semrate(g) = sem(cellprop.fieldMax{g}(goodcells));
%     end
%     b = bar(1:3,meanrate);
%     b(1).FaceColor = 'k';
%     b(1).LineStyle = 'none';
%     set(gca,'Ylim',[0 10],'Xlim',[0 4]);
%     hold on;
%     errorbar(1:3,meanrate,semrate,'Color','k','LineStyle','none');

%to plot preference for half1 or 2 in unit of SD
% max1 = max(cellprop.field{2}(:,1:50),[],2);
% max2 = max(cellprop.field{2}(:,51:100),[],2);
% for icell = 1:size(cellprop.field{2},1)
% SEfield(icell) = cellprop.fieldSE{2}(icell,min(size(cellprop.field{2},2),floor(cellprop.fieldPos{2}(icell))+1));
% end
% figure;
% histogram((max2(goodcells)-max1(goodcells))./SEfield(goodcells)',-5.1:0.2:5.1,'EdgeColor','k','FaceColor','k');
% hold on;histogram((max2(goodcells & (mod(maxpos(goodcells),50)<5 | mod(maxpos(goodcells),50)>45))-max1(goodcells & (mod(maxpos(goodcells),50)<5 | mod(maxpos(goodcells),50)>45)))./SEfield(goodcells & (mod(maxpos(goodcells),50)<5 | mod(maxpos(goodcells),50)>45))',-5.1:0.2:5.1,'EdgeColor','g','FaceColor','g');
% set(gca,'Ylim',[-5 120],'Xlim',[-5 5])

nphs = size(cellprop.field2dXthetapos{2},2);
goodthetacells = (cellprop.ZXmaxthetapos{2}'>=zth  & cellprop.ZXminthetapos{2}'>=zth) &...
                 (cellprop.ZXmaxthetapos{1}'>=zth_lowhigh  & cellprop.ZXminthetapos{1}'>=zth_lowhigh) &...
                 (cellprop.ZXmaxthetapos{3}'>=zth_lowhigh  & cellprop.ZXminthetapos{3}'>=zth_lowhigh);% & sum(cellprop.fieldXcorrtheta{1}>Xcorrmin,2)'==size(cellprop.fieldXcorrtheta{1},2)  & sum(cellprop.fieldXcorrtheta{3}>Xcorrmin,2)'==size(cellprop.fieldXcorrtheta{3},2));% & cellprop.thetareliabilityCorr{2}' >=thetaReliabilitymin);

goodthetacells_low = goodthetacells & (cellprop.ZXmaxthetapos{1}'>=zth  & cellprop.ZXminthetapos{1}'>=zth);%
goodcells_low = ((cellprop.fieldXcorrMax{1}-cellprop.fieldXcorrMax{2})./cellprop.fieldXcorrMaxSE{1})<-zth;
goodthetacells_high = goodthetacells & (cellprop.ZXmaxthetapos{3}'>=zth  & cellprop.ZXminthetapos{3}'>=zth);%
goodcells_high = ((cellprop.fieldXcorrMax{3}-cellprop.fieldXcorrMax{2})./cellprop.fieldXcorrMaxSE{3})>zth;

goodthetaphscells = (cellprop.ZPhsmaxtheta{2}'>=zth & cellprop.ZPhsmintheta{2}'>=zth);
% goodthetacells = goodthetacells & goodthetaphscells;

crossthetaprec = cellprop.PhsXcrossthetapos{2}';%mean(cellprop.field2dXthetaposNorm{2}(:,3:7),2)-mean(cellprop.field2dXthetaposNorm{2}(:,12:16),2);
[~,thetasort] = sort(crossthetaprec,'ascend');
thetasort = isort;

for g = [2 1 3]
    thetafield{g} = zeros(size(cellprop.field2dXPhstheta{g}));
    for icell = 1:size(cellprop.field2dXPhstheta{g},1)
        if ~isnan(floor(cellprop.fieldCOM{2}(icell)))
            fieldtheta = squeeze(cellprop.field2dXPhstheta{g}(icell,:,:));
            fieldthetaref = squeeze(cellprop.field2dXPhstheta{2}(icell,:,:));
            thetafield{g}(icell,:,:) = circshift((fieldtheta-min(fieldthetaref(:)))/(max(fieldthetaref(:))-min(fieldthetaref(:))),50 - floor(cellprop.fieldCOM{2}(icell)),2);
            %thetafield{g}(icell,:,:) = circshift((fieldtheta-min(fieldthetaref(:)))/(max(fieldthetaref(:))-min(fieldthetaref(:))),50 - floor(cellprop.fieldCOM{2}(icell)),2);
        end
    end
    thetafieldcorr{g} = zeros(size(cellprop.field2dXcorrtheta{g}));
    for icell = 1:size(cellprop.field2dXcorrtheta{g},1)
        if ~isnan(floor(cellprop.fieldCOM{2}(icell)))
            fieldtheta = squeeze(cellprop.field2dXcorrtheta{g}(icell,:,:));
            thetafieldcorr{g}(icell,:,:) = fieldtheta;
        end
    end
end

figure('Name','Theta precession1')
for g = [2 1 3]
    thetamod = (cellprop.field2dXcorrthetamax{g}(thetasort(goodcells(thetasort)  & goodthetacells(thetasort)),:)-repmat(mean(cellprop.field2dXcorrthetamax{g}(thetasort(goodcells(thetasort) & goodthetacells(thetasort)),:),2),[1 nphs]))';
    thetamodref = (cellprop.field2dXcorrthetamax{2}(thetasort(goodcells(thetasort)  & goodthetacells(thetasort)),:)-repmat(mean(cellprop.field2dXcorrthetamax{2}(thetasort(goodcells(thetasort) & goodthetacells(thetasort)),:),2),[1 nphs]))';
    subplot(5,9,(g-1)*3+1:(g-1)*3+2);
    imagesc(thetamod);
    set(gca,'Clim',[-2 2],'Ydir','normal');
    colormap(RedWhiteBlue);
    subplot(5,9,(g-1)*3+3);
    plot(nanmean(thetamod,2),1:nphs,'k');
    hold on;plot(nanmean(thetamodref,2),1:nphs,'k--');
    set(gca,'Ylim',[1 nphs]);
    subplot(5,9,9+((g-1)*3+1:(g-1)*3+2));
    h = histogram(cellprop.PhsXcrossthetapos{g}(goodcells  & goodthetacells),0:20:360);
    thetaprecPhase = h.Values./sum(goodcells);
    bar(thetaprecPhase,'k');
    set(gca,'Xlim',[-1 19]);
    subplot(5,9,36+((g-1)*3+1));
    h = histogram(360/18*cellprop.phsfieldPos{g}(goodcells  & goodthetaphscells),0:20:360);
    thetaPhase = h.Values./sum(goodcells);
    bar(thetaPhase,'k');
    subplot(5,9,36+((g-1)*3+2));
    normthetaPhsfield = (cellprop.field2dPhstheta{g}(goodcells & goodthetaphscells,:)./repmat(mean(cellprop.field2dPhstheta{g}(goodcells & goodthetaphscells,:),2),[1 size(cellprop.field2dPhstheta{g},2)]));
    ciplot(gca,[],nanmean(normthetaPhsfield,1),nanstd(normthetaPhsfield,[],1)/sqrt(size(normthetaPhsfield,1)),0.5,'k');
    
%     h = histogram(cellprop.PhsthetaMod{g}(goodcells  & goodthetaphscells)-1,0:0.05:2);
%     thetaPhase = h.Values./sum(goodcells);
%     bar(0:0.05:1.95,thetaPhase,'k');
%     set(gca,'Xlim',[-0.05 2]);
    
    subplot(5,9,18+((g-1)*3+1:(g-1)*3+2));
    h = histogram(cellprop.Xmaxthetapos{g}(goodcells  & goodthetacells)-cellprop.Xminthetapos{g}(goodcells  & goodthetacells),0:20);
    thetaPhase = h.Values./sum(goodcells);
    bar(thetaPhase,'k');
    set(gca,'Xlim',[-1 21]);
    subplot(5,9,27+((g-1)*3+1:(g-1)*3+2));
    scatter(360/18*cellprop.phsfieldPos{g}(goodcells  & goodthetacells & goodthetaphscells),cellprop.PhsXcrossthetapos{g}(goodcells  & goodthetacells & goodthetaphscells))
    set(gca,'Xlim',[0 360],'Ylim',[0 360]);
%     thetamax = max(cellprop.field2dXcorrthetamax{g}(goodcells  & goodthetacells,:),[],2);
%     thetamin = min(cellprop.field2dXcorrthetamax{g}(goodcells  & goodthetacells,:),[],2);
%     polarplot(cellprop.PhsXcrossthetapos{g}(goodcells  & goodthetacells),thetamax-thetamin,'.')
end

figure('Name','Theta precession2')
for g = [2 1 3]
    subplot(4,3,g);
    thetamap = squeeze(nanmean(thetafield{g}(goodcells  & goodthetacells & crossthetaprec >= 0 & crossthetaprec < 120,:,:),1));
    if g==2
        disp(sum(goodcells  & goodthetacells & crossthetaprec >= 0 & crossthetaprec < 120)/sum(goodcells));
    end
    imagesc(-thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.01,0.05),1:nphs,'w');
    set(gca,'Clim',[-0.8 -0.2],'Xlim',[30 70],'Ydir','normal');
    subplot(4,3,g+3);
    thetamap = squeeze(nanmean(thetafield{g}(goodcells  & goodthetacells & crossthetaprec >= 120 & crossthetaprec < 240,:,:),1));
    if g==2
        disp(sum(goodcells  & goodthetacells & crossthetaprec >= 120 & crossthetaprec < 240)/sum(goodcells));
    end
    imagesc(-thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.01,0.05),1:nphs,'w');
    set(gca,'Clim',[-0.8 -0.2],'Xlim',[30 70],'Ydir','normal');
    subplot(4,3,g+6);
    thetamap = squeeze(nanmean(thetafield{g}(goodcells  & goodthetacells & crossthetaprec >= 240 & crossthetaprec < 360,:,:),1));
    if g==2
        disp(sum(goodcells  & goodthetacells & crossthetaprec >= 240 & crossthetaprec < 360)/sum(goodcells));
    end
    imagesc(-thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.01,0.05),1:nphs,'w');
    set(gca,'Clim',[-0.8 -0.2],'Xlim',[30 70],'Ydir','normal');
    subplot(4,3,g+9);
    thetamap = squeeze(nanmean(thetafield{g}(goodcells  & goodthetacells,:,:),1));
    if g==2
        disp(sum(goodcells  & goodthetacells)/sum(goodcells));
    end
    imagesc(-thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.01,0.05),1:nphs,'w');
    set(gca,'Clim',[-0.8 -0.2],'Xlim',[30 70],'Ydir','normal');
    colormap('bone')
end

figure('Name','Theta precession3')
for g = [2 1 3]
    subplot(4,3,g);
    thetamap = squeeze(nanmean(thetafieldcorr{g}(goodcells  & goodthetacells & crossthetaprec >= 0 & crossthetaprec < 120,:,:),1));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.01,0.05),1:nphs,'w');
    set(gca,'Clim',[0 1],'Xlim',[30 70],'Ydir','normal');
    subplot(4,3,g+3);
    thetamap = squeeze(nanmean(thetafieldcorr{g}(goodcells  & goodthetacells & crossthetaprec >= 120 & crossthetaprec < 240,:,:),1));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.01,0.05),1:nphs,'w');
    set(gca,'Clim',[0 1],'Xlim',[30 70],'Ydir','normal');
    subplot(4,3,g+6);
    thetamap = squeeze(nanmean(thetafieldcorr{g}(goodcells  & goodthetacells & crossthetaprec >= 240 & crossthetaprec <= 360,:,:),1));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,1),1:nphs,'w');
    set(gca,'Clim',[0 1],'Xlim',[30 70],'Ydir','normal');
    subplot(4,3,g+9);
    thetamap = squeeze(nanmean(thetafieldcorr{g}(goodcells  & goodthetacells,:,:),1));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.01,0.05),1:nphs,'w');
    set(gca,'Clim',[0 1],'Xlim',[30 70],'Ydir','normal');
    colormap('jet')
end

end
end

%to look at single cell precession: run this in debug mode
% f1 = figure;f2 = figure;
% for icell = 1:17000
% if goodthetacells(icell) && goodcells(icell) && cellprop.animal(icell) == 7 && cellprop.iseries(icell) == 3
% map2d = squeeze(cellprop.field2dXPhstheta{2}(icell,:,:));
% figure(f1);
% colormap(bone)
% subplot(4,1,[2 3]);
% hold off;
% imagesc(-(map2d - min(map2d(:)))./(max(map2d(:)) - min(map2d(:))));
% hold on;plot(cellprop.field2dXthetapos{2}(icell,:),1:18,'w');
% set(gca,'Clim',[-1 0],'Ydir','normal');
% subplot(4,1,1);
% hold off;
% plot(mean(map2d,1));
% subplot(4,1,4);
% hold off;
% plot(mean(map2d,2));
% figure(f2);
% colormap(parula)
% hold off;
% imagesc(squeeze(cellprop.field2dXcorrtheta{2}(icell,:,:)));
% hold on;plot(cellprop.field2dXcorrthetamax{2}(icell,:),1:18,'w');
% hold on;plot(cellprop.field2dXcorrthetamax{2}(icell,:)+cellprop.field2dXcorrthetamaxSE{2}(icell,:),1:18,'y');
% hold on;plot(cellprop.field2dXcorrthetamax{2}(icell,:)-cellprop.field2dXcorrthetamaxSE{2}(icell,:),1:18,'y');
% set(gca,'Clim',[-1 1],'Ydir','normal');
% pause;
% end
% end