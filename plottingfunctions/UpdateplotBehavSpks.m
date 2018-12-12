function PlotVar = UpdateplotBehavSpks(PlotVar,EXP,Layout)
PlotVar.Plots = [];
es = EXP.data.es;

es.spikeTrain = circshift(es.spikeTrain,[-PlotVar.delayT 0]);

page = 1;
win = 1;
Layout.DivideWindow(page, win, size(PlotVar.ChosenObj,1)*size(PlotVar.ChosenCell,2)*size(PlotVar.ChosenContrast,2), size(PlotVar.ChosenGain,2)*size(PlotVar.ChosenRoomlength,2)*size(PlotVar.ChosenOutcome,2));

nplot =  0;
nplotfast1Dmap = 0;
maxfast1Dmap = 0;

trialID = es.trialID;
for c = 1:size(PlotVar.ChosenContrast,2)
    for g = 1:size(PlotVar.ChosenGain,2)
        for r = 1:size(PlotVar.ChosenRoomlength,2)
            for o = 1:size(PlotVar.ChosenOutcome,2)
                es.trialID = trialID;
                Fnoblanks = true;
                Fnoafterblanks = true;
                tidx = EXP.getSubsets(PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o),PlotVar.speed_th, Fnoblanks, Fnoafterblanks);
                exptidx = 0;
                seriesid = unique(es.series);
                for sid = 1:numel(seriesid)
                    exptidx = exptidx + (es.series == seriesid(sid) & ismember(es.iexp,PlotVar.explist(sid,:)));
                end
                halfidx = ones(size(es.traj));%ismember(es.halfID,[1]);%es.completetrial;%
                tidx = tidx & exptidx & halfidx;
%                 tidx = tidx & EXP.Bayes.Spdbin==PlotVar.thetaChannel;
%                 tidx = tidx & EXP.data.es.smthBallSpd>0;
                
                tidxAllphs = tidx;
                trialIDs = unique(es.trialID(tidx));
                es.trialID(~ismember(es.trialID,trialIDs)) = 0;
                for itr = 1:numel(trialIDs)
                    es.trialID(es.trialID==trialIDs(itr)) = itr;
                end
                
                for cell = 1:size(PlotVar.ChosenCell,2)                    
                    icell = PlotVar.ChosenCell(:,cell);
                    spktrain = es.spikeTrain(:,icell);
                    spktrain = sum(spktrain,2);
                    
                    Fdispmat = PlotVar.Fdispmat;
                    if Fdispmat
                        spkidx = true(size(spktrain));
                    else
                        spkidx = spktrain>0;
                    end
                    
                    try
                        varX = es.(PlotVar.ChosenVarX{1})(:);
                        varY = es.(PlotVar.ChosenVarY{1})(:);
                    catch
                        error('selected variable not a field of the EXP.data.es structure')
                    end
                    
                    for k = 1:size(PlotVar.ChosenObj,1)
                        wintitle = ['Cell: ' num2str(PlotVar.ChosenCell(:,cell)') ' Contrast: ' num2str(EXP.SubsetVal.contrast(PlotVar.ChosenContrast(:,c)')) ' Gain: ' num2str(EXP.SubsetVal.gain(PlotVar.ChosenGain(:,g)'))];
                        Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}.Title = wintitle;
                                                
                        switch PlotVar.ChosenObj{k}                            
                            case 'Raster'
                                nplot = nplot + 1;
%                                 if ~isempty(EXP.CellInfo.trajdecCell_Xave{icell})
%                                     varX = EXP.CellInfo.trajdecCell_Xmax{icell}-1;
%                                 else
%                                     varX = 0*ones(size(varX));
%                                 end
%                               
%                                 varX = EXP.Bayes.predicted{1};%
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                                idxsegments = find(abs(diff(varX(tidx & es.outcome == 2)))>10);
                                varxx = varX(tidx & es.outcome == 2);
                                varyy = varY(tidx & es.outcome == 2);
                                for itrial = 1:numel(idxsegments)-1
                                    hold on;plot(varxx(idxsegments(itrial)+1:idxsegments(itrial+1)-1),varyy(idxsegments(itrial)+1:idxsegments(itrial+1)-1),'color',[0.9 0.9 0.9]);
                                end
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & es.outcome == 2, varX, varY, spktrain, Fdispmat, [0 0 0]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & es.outcome == 3, varX, varY, spktrain, Fdispmat, [1 0 0]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & es.outcome == 4, varX, varY, spktrain, Fdispmat, [0 1 0]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & es.outcome == 0, varX, varY, spktrain, Fdispmat, [1 0 0]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & es.outcome == 1, varX, varY, spktrain, Fdispmat, [0 0 1]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & es.outcome == 5, varX, varY, spktrain, Fdispmat, [0.5 0.5 0.5]);
                            case 'fast1Dmap' 
                                varX = es.(PlotVar.ChosenVarX{1})(:);%smthInTime(es.(PlotVar.ChosenVarX{1})(:), 60, 150, 'same', [], 'boxcar_centered');
                                varX = normalise1var(varX,100,[],[0 100]);
                                dx = 1;
                                samplerate = 1./es.sampleSize;%60;
                                nbXbinsmth = round(1/(4/100));%round(1/(PlotVar.Xbinsize/100));
                                spktrain(~isnan(spktrain)) = smthInTime(spktrain(~isnan(spktrain)), mean(samplerate), 15, 'same', [], 'boxcar_centered');
                                
                                [map,x] = fast1Dmap(varX(tidx), spktrain(tidx), dx, samplerate,nbXbinsmth,es.CircularMaze);
                                maxfast1Dmap = max(max(map),maxfast1Dmap);
                                
%                                 mat = zeros(1001);
%                                 for i = 1:1001
%                                     for j = 1:1001
%                                         mat(i,j) = map(mod(floor(((180+atan2d((j-501),(i-501)))/360)*100)-1,100)+1);
%                                     end
%                                 end
%                                 figure;imagesc(fliplr(mat))
%                                 set(gca,'Clim',[0 max(mat(:))],'Ydir','normal')
                                
                                nplot = nplot + 1;
                                if nplotfast1Dmap == 0
                                    nplotfast1Dmap = nplot;
                                    PlotVar.Plots{nplotfast1Dmap} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                end                                
%                                 PlotVar.Plots{nplotfast1Dmap}.PlotVector(x, map, [], [], true, 'Color', 'k', 'Ylim', [0 1.1*max(map)+1*(max(map)==0)]);

                                PlotVar.Plots{nplotfast1Dmap}.PlotVector(x, map, [], [], true, 'Ylim', [0 1.1*maxfast1Dmap+1*(maxfast1Dmap==0)]);
                                
%                                 [SInfo, SInfoperSpk, ~] = SpatialInformation(varX(tidx), spktrain(tidx), dx, samplerate,nbXbinsmth,es.CircularMaze);
%                                 disp(SInfoperSpk)
%                                 
%                                 if ~isempty(EXP.CellInfo.trajdecCell_Xave{icell})
%                                     varX = EXP.CellInfo.trajdecCell_Xave{icell}-1;%EXP.CellInfo.trajdecAll_Xmax{1}-1;%
%                                 else
%                                     varX = 0*ones(size(varX));
%                                 end
%                                 [map,x] = fast1Dmap(varX(tidx), spktrain(tidx), dx, samplerate,nbXbinsmth,es.CircularMaze);
%                                 maxfast1Dmap = max(max(map),maxfast1Dmap);
%                                 PlotVar.Plots{nplotfast1Dmap}.PlotVector(x, map, [], [], true, 'Ylim', [0 1.1*maxfast1Dmap+1*(maxfast1Dmap==0)]);
%                                 
%                                 [SInfo, SInfoperSpk, ~] = SpatialInformation(varX(tidx), spktrain(tidx), dx, samplerate,nbXbinsmth,es.CircularMaze);
%                                 disp(SInfoperSpk)
                                
                                [~,imax] = max(EXP.CellInfo.fieldXcorr{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(icell,:),[],2);
                                text(PlotVar.Plots{nplotfast1Dmap}.ax{1},0.5,0.1+g*0.1,['Xshift:' num2str(imax(1)-51)]);
                                
%                                 [map,x] = fast1Dmap(EXP.Bayes.predicted{1}(tidx), spktrain(tidx), dx, samplerate,nbXbinsmth,es.CircularMaze);
%                                 hold on;plot(x,map);

%                                 [map,~,~] = fast2Dmap(varX(tidx),EXP.Bayes.predicted{1}(tidx),spktrain(tidx),dx,dx,samplerate,nbXbinsmth,nbXbinsmth, es.CircularMaze);

                                    
                            case '1Dmap'
                                nplot = nplot + 1;
                                contidx = PlotVar.ChosenContrast(:,c);
                                gainidx = PlotVar.ChosenGain(:,g);
                                roomlengthidx = PlotVar.ChosenRoomlength(:,r);
                                outcomeidx = PlotVar.ChosenOutcome(:,o);
                                distri_Th = 0.99;
%                                 SinfoperSpike = EXP.CellInfo.SpatialInfoPerSpike{contidx,gainidx,roomlengthidx,outcomeidx}(icell);
%                                 SinfoperSpikeZscore = (SinfoperSpike)/quantile(EXP.CellInfo.SpatialInfoPerSpikeRef{contidx,gainidx,roomlengthidx,outcomeidx}(icell,:),distri_Th);
%                                 Sinfo = EXP.CellInfo.SpatialInfo{contidx,gainidx,roomlengthidx,outcomeidx}(icell);
%                                 SinfoZscore = (EXP.CellInfo.SpatialInfo{contidx,gainidx,roomlengthidx,outcomeidx}(icell))/quantile(EXP.CellInfo.SpatialInfoRef{contidx,gainidx,roomlengthidx,outcomeidx}(icell,:),distri_Th);
%                                 SSIZscore = (EXP.CellInfo.SSI{contidx,gainidx,roomlengthidx,outcomeidx}(icell))/quantile(EXP.CellInfo.SSIRef{contidx,gainidx,roomlengthidx,outcomeidx}(icell,:),distri_Th);
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, 2});%((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                for cc = 1:size(PlotVar.ChosenContrast,1)
                                    for gg = 1:size(PlotVar.ChosenGain,1)
                                        for rr = 1:size(PlotVar.ChosenRoomlength,1)
                                            for oo = 1:size(PlotVar.ChosenOutcome,1)
                                                Plot1DmapSpks(PlotVar.Plots{nplot}, EXP.maps1d.trajPercent{PlotVar.ChosenContrast(cc,c),PlotVar.ChosenGain(gg,g),PlotVar.ChosenRoomlength(rr,r),PlotVar.ChosenOutcome(oo,o)}, icell);                                                
                                            end
                                        end
                                    end
                                end
                                try
%                                 zscore = (EXP.CellInfo.field{contidx,gainidx,roomlengthidx,outcomeidx}(icell,floor(EXP.CellInfo.fieldPos{contidx,gainidx,roomlengthidx,outcomeidx}(icell))+1)...
%                                          - EXP.CellInfo.rate{contidx,gainidx,roomlengthidx,outcomeidx}(icell))/max(EXP.CellInfo.fieldStdRef{contidx,gainidx,roomlengthidx,outcomeidx}(icell,:));
                                text(0.5,0.1,['zscore:' num2str(zscore) ...
                                              ' / mean rate:' num2str(EXP.CellInfo.rate{contidx,gainidx,roomlengthidx,outcomeidx}(icell)) ...
                                              ' / SSI:' num2str(EXP.CellInfo.SSI{contidx,gainidx,roomlengthidx,outcomeidx}(icell)) ...
                                              ' / SSI zscore:' num2str(SSIZscore) ...
                                              ' / Spatial info:' num2str(SinfoperSpike) ...
                                              ' / Spatial info zscore:' num2str(SinfoperSpikeZscore) ...
                                              ' / COM:' num2str(EXP.CellInfo.fieldCOM{contidx,gainidx,roomlengthidx,outcomeidx}(icell)) ...
                                              ' / Pos:' num2str(EXP.CellInfo.fieldPos{contidx,gainidx,roomlengthidx,outcomeidx}(icell)) ...
                                              ' / PosSE:' num2str(EXP.CellInfo.fieldPosSE{contidx,gainidx,roomlengthidx,outcomeidx}(icell))]);
                                catch
                                end
                            case '2DmapXSpd'
                                map = EXP.maps2d.trajPercent_smthBallSpd{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModel;
                                vecmean = EXP.maps2d.trajPercent_smthBallSpd{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelXpos;
                                vecSE = EXP.maps2d.trajPercent_smthBallSpd{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelXpos;
                                
                                x = 1:size(map,2);
                                y = 1:size(map,1);
                                
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                                PlotVar.Plots{nplot}.palette = 'parula';%'hot';%
                                if max(map(:)) > 0
                                    PlotVar.Plots{nplot}.PlotMatrix(x, (y), map, [], true,'Clim',[0 max(map(:))]);
                                    PlotVar.Plots{nplot}.PlotVector(vecmean, y, [], [], true);
                                    PlotVar.Plots{nplot}.PlotVector(mod(vecmean+vecSE,size(map,2)), y, [], [], true);
                                    PlotVar.Plots{nplot}.PlotVector(mod(vecmean-vecSE,size(map,2)), y, [], [], true);
                                end
                            case '2DmapXTheta'
                                if isfield(es, 'LFPphase')
                                    map = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModel;
                                    vecmean = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelXmax;                                    
                                    vecSE = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelXmax;                                    
%                                     vecmeanphs = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelYpos;                                    
%                                     vecmeanphsSE = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelYpos;                                    
                                    [~,imaxprec] = max(abs(vecmean-mean(vecmean)));
                                    zCOM = (vecmean(imaxprec) - mean(vecmean))/vecSE(imaxprec);
                                    x = 1:size(map,2);
                                    y = linspace(0,360,size(map,1));
                                    
                                    nplot = nplot + 1;
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                    PlotVar.Plots{nplot}.palette = 'parula';%'hot';%
                                    if max(map(:)) > 0
                                        PlotVar.Plots{nplot}.PlotMatrix(x, (y), map, [], true,'Clim',[0 max(map(:))]);
                                        PlotVar.Plots{nplot}.PlotVector(vecmean, y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(mod(vecmean+vecSE,size(map,2)), y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(mod(vecmean-vecSE,size(map,2)), y, [], [], true);
                                        
%                                         PlotVar.Plots{nplot}.PlotVector(x, vecmeanphs*360/size(map,1), [], [], true);
%                                         PlotVar.Plots{nplot}.PlotVector(x, mod((vecmeanphs+vecmeanphsSE)*360/size(map,1),360), [], [], true);
%                                         PlotVar.Plots{nplot}.PlotVector(x, mod((vecmeanphs-vecmeanphsSE)*360/size(map,1),360), [], [], true);
                                    end
                                    
                                    xmap = mean(map,1);
                                    [~,imax] = max(xmap);
                                    [~,fieldidx] = findfield(xmap,1);
                                    fieldimax = find(fieldidx == imax);
                                    
                                    fielddist_to_max = unwrap(fieldidx/numel(xmap)*2*pi)*numel(xmap)/(2*pi);
                                    fielddist_to_max = fielddist_to_max - fielddist_to_max(fieldimax);
                                    fielddist_to_max(1:fieldimax) = fielddist_to_max(1:fieldimax)/abs(min(fielddist_to_max));
                                    fielddist_to_max(fieldimax+1:end) = fielddist_to_max(fieldimax+1:end)/abs(max(fielddist_to_max));
                                    
                                    infielddist = NaN(numel(xmap),1);
                                    infielddist(fieldidx) = fielddist_to_max;
                                    infielddist = [infielddist;infielddist(1)];
                                    infielddist_diff  = diff(infielddist);
                                    varDist = infielddist(floor(varX)+1)+mod(varX,1).*infielddist_diff(floor(varX)+1);
                                    spktraintemp = spktrain;
                                    spktraintemp(isnan(varDist)) = 0;
                                    
%                                     map_xphs = map/sum(map(:))*sum(spktrain(tidxAllphs));
%                                     [rho,pval] = circ_corrcc2(y/(numel(y)*dy)*2*pi, x/(numel(x)*dx)*2*pi, map_xphs);
%                                     [rho,pval] = circ_corrcl(es.spikePhase(tidxAllphs & spktraintemp>0.5,icell)/360*2*pi, varXtemp(tidxAllphs & spktraintemp>0.5));

%                                     kfold = 20;
%                                     [slope,phi0,rho,slopeSE,phi0SE,rhoSE] = circularlinearfit(es.spikePhase(tidxAllphs & spktraintemp>0.5,icell)/360*2*pi,varDist(tidxAllphs & spktraintemp>0.5),[],kfold);
%                                     phi0 = mod(phi0+2*pi,2*pi)*360/(2*pi);
%                                     phi0SE = mod(phi0SE+2*pi,2*pi)*360/(2*pi);
%                                     
%                                     phs = repmat(2*pi/360*y',[1 numel(x)]);
%                                     xphs = repmat(-floor(numel(xmap)/2)+1:floor(numel(xmap)/2),[numel(y) 1])/numel(fieldidx);
%                                     phsmap = map;
%                                     phsmap(:,~ismember(1:size(map,2),fieldidx)) = 0;
%                                     phsmap = circshift(phsmap,-imax+floor(numel(xmap)/2),2);
%                                     [slope,phi0,rho] = circularlinearfit(phs(:),xphs(:),phsmap(:));
                                    
                                    slope = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelslopeXY;
                                    phi0 = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelphi0XY;
                                    rho = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelrhoXY;
                                    slopeSE = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelslopeXY;
                                    phi0SE = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelphi0XY;
                                    rhoSE = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelrhoXY;
                                    
                                    
%                                     [map1d,x] = fast1Dmap(theta(tidxAllphs),spktrain(tidxAllphs),dy,samplerate,nbYbinsmth,false);
                                    map1d = mean(map,2);
                                    map_phs = map1d/sum(map1d)*sum(spktrain(tidxAllphs));
                                    text(1,1,['rho = ' num2str(rho) ';z=' num2str(rho/rhoSE) '  slope = ' num2str(slope) ';z=' num2str(slope/slopeSE) '  offset = ' num2str(phi0) ';SE=' num2str(phi0SE) '     zCOM = ' num2str(zCOM)]);
                                end
                            case 'Theta'
                                if isfield(es, 'LFPphase')
                                    xcorrtheta = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelX;
                                    Xpos = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelXmax;
                                    XposSE = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SEcorrModelXmax;
%                                     xcorrthetatemp = xcorrtheta;
%                                     xcorrthetatemp(:,[1:25 75:100]) = 0;
%                                     Xpos = getCircularAverage(xcorrthetatemp',0,1);
                                    [~,imaxprec] = max(abs(Xpos-mean(Xpos)));
                                    zCOM = (Xpos(imaxprec)-mean(Xpos))/XposSE(imaxprec);
%                                     dx = 1;
%                                     dy = 20;%5;%
%                                     samplerate = 1./es.sampleSize(tidxAllphs);%60;
%                                     nbXbinsmth = round(1/(PlotVar.Xbinsize/100));
%                                     nbYbinsmth = round(1/(2/(360/dy)));
%                                     theta = mod(es.LFPphase(:,min(end,PlotVar.thetaChannel)),360);%mod(EXP.Bayes.LFPphase,360);%es.smthBallSpd; %mod(EXP.Bayes.LFPphase(tidxAllphs),360);%mod(es.LFPphase(tidxAllphs,PlotVar.thetaChannel)-180,360)
%                                     [map,x,y] = fast2Dmap(varX(tidxAllphs),theta(tidxAllphs),spktrain(tidxAllphs),dx,dy,samplerate,nbXbinsmth, nbYbinsmth,es.CircularMaze);
%                                     xmap = mean(map,1);
%                                     [~,imax] = max(xmap);
%                                     varXtemp = mod(varX - imax + max(varX)/2,max(varX))- max(varX)/2;
%                                     xmap = circshift(xmap,-imax+floor(numel(xmap)/2));
%                                     qth = quantile(xmap,0.5);
%                                     fieldstart = x(find(xmap(1:floor(numel(xmap)/2))<=mean(xmap),1,'last'))-floor(numel(xmap)/2);
%                                     fieldend = x(find(xmap(floor(numel(xmap)/2)+1:end)<=mean(xmap),1,'first'));
%                                     spktraintemp = spktrain;
%                                     spktraintemp(varXtemp<=fieldstart | varXtemp>=fieldend) = 0;
% %                                     map_xphs = map/sum(map(:))*sum(spktrain(tidxAllphs));
% %                                     [rho,pval] = circ_corrcc2(y/(numel(y)*dy)*2*pi, x/(numel(x)*dx)*2*pi, map_xphs);
%                                     [rho,pval] = circ_corrcc(es.spikePhase(tidxAllphs & spktraintemp>0.5,icell)/360*2*pi, varXtemp(tidxAllphs & spktraintemp>0.5)/max(varX)*2*pi);
% %                                     [map1d,x] = fast1Dmap(theta(tidxAllphs),spktrain(tidxAllphs),dy,samplerate,nbYbinsmth,false);
%                                     map1d = mean(map,2);
%                                     map_phs = map1d/sum(map1d)*sum(spktrain(tidxAllphs));
%                                     [rtest_pval,rtest_z] = circ_rtest(y/360*2*pi,map_phs,dy/360*2*pi);
%                                     
% %                                     map = repmat(map,[3 1]);
% %                                     y = [y 360+y 2*360+y];
%                                     
% %                                     map = EXP.maps2d.trajPercent_smthBallSpd{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModel;
% %                                     y = EXP.maps2d.trajPercent_smthBallSpd{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.binsY;
% %                                     
%                                     map = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModel;
%                                     vecmean = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelXpos;                                    
%                                     vecSE = EXP.maps2d.trajPercent_LFPphase2{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelXpos;                                    
                                    
                                    nplot = nplot + 1;
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                    PlotVar.Plots{nplot}.palette = 'parula';%'hot';%
                                    if max(map(:)) > 0
                                        PlotVar.Plots{nplot}.PlotMatrix(x, (y), xcorrtheta, [], true,'Clim',[0 1]);
%                                         PlotVar.Plots{nplot}.PlotMatrix(x, (y), map, [], true,'Clim',[0 max(map(:))]);
                                        PlotVar.Plots{nplot}.PlotVector(Xpos, y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(Xpos+XposSE, y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(Xpos-XposSE, y, [], [], true);
                                    end
                                    text(1,1,['zCOM = ' num2str(zCOM)]);
                                end
                            case 'POTH'
                                if isfield(es, 'LFPphase')
                                    roomlength = 100;
                                    dx = 1;
                                    dy = 36;
                                    samplerate = 1./es.sampleSize(tidxAllphs);%60;
                                    xbins = floor(roomlength/dx);
                                    x = linspace(0, roomlength, xbins);
                                    tau = 20;
                                    nbphsbins = 18;
                                    nbphsbinsmth = round(1/(1/nbphsbins));
                                    nbXbinsmth = round(1/(PlotVar.Xbinsize/100));
                                    nbYbinsmth = round(1/(1/(360/dy)));
                                    POTH = zeros(nbphsbins,xbins);
                                    
                                    minidx0 = find([0;abs(diff(mod(es.LFPphase(:,min(end,PlotVar.thetaChannel)),360)))>180]);%
                                    thetaperiod0 = diff(minidx0);
                                    thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
                                    
%                                     minidx0 = min(minidx0 + round(thetaperiod/2),numel(tidx));
                                    
                                    theta = zeros(size(es.LFPphase(:,min(end,PlotVar.thetaChannel))));
                                    for tt = 1:numel(minidx0)-1
                                        theta(minidx0(tt):minidx0(tt+1)) = 360*(0:(minidx0(tt+1)-minidx0(tt)))/thetaperiod0(tt);
                                    end
                                    Xrange = max(es.trajPercent);
                                    
                                    idx = find(tidx);
                                    minidx = minidx0(ismember(minidx0,idx));
                                    thetaperiod = thetaperiod0(ismember(minidx0,idx));
                                    phsidx = false(size(tidx));
                                    phsidx(minidx) = true;
                                    minidx = find(phsidx & tidx);
                                    spkPOTH_All = zeros(nbphsbins,numel(minidx));
                                    XPOTH_All = zeros(nbphsbins,numel(minidx));
                                    theta_All = zeros(nbphsbins,numel(minidx));
                                    for tt = 1:numel(minidx)
                                        if minidx(tt) > 1 && minidx(tt)+tau-1 < size(spktrain,1)
                                            spkthetavec = spktrain(max(1,minidx(tt)):min(minidx(tt)+tau-1,end),:);
                                            spkthetavec = (interp1((0:tau-1)/thetaperiod(tt),spkthetavec,0:1/nbphsbins:(1-(1/nbphsbins))))';
                                            Xthetavec = es.trajPercent(max(1,minidx(tt)):min(minidx(tt)+tau-1,end),:);
                                            Xthetavec = unwrap(Xthetavec/(Xrange)*2*pi)*Xrange/(2*pi);
                                            Xthetavec = mod((interp1((0:tau-1)/thetaperiod(tt),Xthetavec,0:1/nbphsbins:(1-(1/nbphsbins)))),Xrange)';
                                            spkPOTH_All(:,tt) = spkthetavec;
                                            XPOTH_All(:,tt) = Xthetavec;
                                            theta_All(:,tt) = 360*(0:1/nbphsbins:(1-(1/nbphsbins)));
                                        end
                                    end
                                    [map,~,~] = fast2Dmap(XPOTH_All(:),theta_All(:),spkPOTH_All(:),dx,360/nbphsbins,nanmean(samplerate),nbXbinsmth,nbphsbinsmth, es.CircularMaze);
                                    
                                    
%                                     [map,~,~] = fast2Dmap(varX(tidxAllphs),theta(tidxAllphs),spktrain(tidxAllphs),dx,dy,samplerate,nbXbinsmth,nbYbinsmth, es.CircularMaze);
                                    [~,imax] = max(sum(map,1));
                                    varXtemp = mod(varX - imax + max(varX)/2,max(varX))- max(varX)/2;
                                    spktraintemp = spktrain;
%                                     spktraintemp(varXtemp<=-10 | varXtemp>=10) = 0;
                                    [rho,pval] = circ_corrcc(theta(tidxAllphs & spktraintemp>0)/360*2*pi, varXtemp(tidxAllphs & spktraintemp>0)/max(varX)*2*pi);
%                                     disp([num2str(rho) ' ' num2str(pval)]);
                                    %                                         minidx = find([0;abs(diff(mod(es.LFPphase(:,PlotVar.thetaChannel)-180,360)))>180]);
                                    POTH = map;
%                                     for xx = 1:xbins
%                                         x1 = (xx-1)*dx;x2 = xx*dx;
%                                         minidx = minidx0;
% %                                         minidx0_2 = find([0;abs(diff(mod(EXP.Bayes.LFPphase{1}+180,360)))>180]);
% %                                         thetaAmp = zeros(size(minidx0));
% %                                         for tt = 1:numel(minidx0)
% %                                             %                                     nexttrough = find(minidx0_2>minidx0(tt),1,'first');
% %                                             if tt ~= numel(minidx0)
% %                                                 nexttrough = find(minidx0_2>minidx0(tt) & minidx0_2<minidx0(tt+1));
% %                                             else
% %                                                 nexttrough = find(minidx0_2>minidx0(tt),1,'first');
% %                                             end
% %                                             nexttrough = nexttrough(abs(LFPtheta(nexttrough)) == max(abs(LFPtheta(nexttrough))));
% %                                             if ~isempty(nexttrough)
% %                                                 thetaAmp(tt) = LFPtheta(minidx0(tt)) - LFPtheta(nexttrough);
% %                                             end
% %                                         end
% %                                         thetamean = mean(thetaAmp);
% %                                         thetasd = std(thetaAmp);
% %                                         if selprobe == 2
% %                                             zth = 1;
% %                                             minidx0(abs(thetaAmp)<thetamean + zth*thetasd) = [];
% %                                         end
% %                                         
%                                         phsidx = false(size(tidx));
%                                         phsidx(minidx) = true;
%                                         minidx = find(phsidx & tidx & es.traj > x1 & es.traj <= x2);   
%                                         POTH_All = zeros(size(POTH,1),size(POTH,2),numel(minidx));
%                                         
%                                         for tt = 1:numel(minidx)
%                                             if minidx(tt)-tau > 1 && minidx(tt)+tau < size(spktrain,1) %&& sum(es.badlick(max(1,minidx(tt)):min(minidx(tt)+24*tau,end))) > 0 %&& sum(es.goodlick(max(1,minidx(tt)-24*tau):min(minidx(tt),end))) > 0
%                                                 thetavec = spktrain(max(1,minidx(tt)-tau):min(minidx(tt)+tau,end),:);%.*(tidx(max(1,minidx(tt)-tau):min(minidx(tt)+tau,end)))));
%                                                 thetavec = (interp1((-tau:tau)/thetaperiod(tt),thetavec,-1/2:1/nbphsbins:(1-(1/nbphsbins))/2))';
%                                                 POTH_All(:,xx,tt) = thetavec;
%                                             end
%                                         end
%                                         POTH(:,xx) = nansum(POTH_All(:,xx,:),3);
%                                         Xsum(:,xx) = sum(~isnan(POTH_All(:,xx,:)),3);
%                                     end
%                                     POTHtemp = special_smooth_2d(repmat(POTH, 3, 3),[1/nbphsbinsmth 1/nbXbinsmth],0,0,[nbphsbins xbins])./special_smooth_2d(repmat(Xsum, 3, 3),[1/nbphsbinsmth 1/nbXbinsmth],0,0,[nbphsbins xbins]);
%                                     POTH = POTHtemp(size(POTH,1) + 1:2*size(POTH,1), size(POTH,2) + 1:2*size(POTH,2));
%                                     POTH = POTH*nanmean(samplerate);

%                                     POTH = repmat(POTH,[2 1]);
                                    y = 360*(0:1/nbphsbins:(1-(1/nbphsbins)));
%                                     y = [y y+max(y)];
%                                     POTHtemp = conv2(repmat(POTH, 3, 3), 1/9*ones(3), 'same');
%                                     POTH = POTHtemp(size(POTH,1) + 1:2*size(POTH,1), size(POTH,2) + 1:2*size(POTH,2));

                                    nplot = nplot + 1;
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                    PlotVar.Plots{nplot}.palette = 'parula';%'hot';
                                    PlotVar.Plots{nplot}.PlotMatrix(x, y, POTH, [], 'Ylabel',num2str(fliplr(POTH)));
                                end
                            case 'POTH2'%'X x Speed map'
                                dx = PlotVar.Xbinsize;
                                samplerate = 1./es.sampleSize;%60;
                                nbXbinsmth = round(1/(PlotVar.Xbinsize/100));
                                spktrain = es.spikeTrain(:,icell);
                                ndT = 100;
                                mapdT = [];
                                for dT = -ndT:ndT
                                    spktrain = circshift(es.spikeTrain(:,icell),[-dT 0]);
                                    [map,x] = fast1Dmap(varX(tidx), spktrain(tidx), dx, samplerate,nbXbinsmth,es.CircularMaze);
                                    if isempty(mapdT)
                                        mapdT = zeros(ndT*2+1,numel(map));
                                    end
                                    mapdT(dT+ndT+1,:) = map;
                                end
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                                PlotVar.Plots{nplot}.palette = 'parula';%'hot';%
                                y = -ndT:ndT;
                                if max(mapdT(:)) > 0
                                    PlotVar.Plots{nplot}.PlotMatrix(x, (y), mapdT/max(mapdT(:)), [], 'Ylabel',num2str(y),'Clim',[0 1]);
                                end
                                PlotVar.Plots{nplot}.PlotVector([x(1) x(end)], [0 0],[],[],true, 'color', 'k');
                                [ymax,xmax] = find(mapdT == max(mapdT(:)));
                                PlotVar.Plots{nplot}.PlotVector(x(xmax), y(ymax),[],[],true, 'Marker','+', 'linestyle', 'none', 'color', 'k');
                                
%                                 if isfield(es, 'LFPphase')
%                                     dx = 1;%PlotVar.Xbinsize;
%                                     dy = 36;
%                                     samplerate = 60;
%                                     nbXbinsmth = round(1/(PlotVar.Xbinsize/100));
%                                     visualspeed = NaN(size(es.trajspeed));
%                                     visualspeed(~isnan(es.trajspeed)) = smthInTime(es.trajspeed(~isnan(es.trajspeed)), samplerate, 0, 'same', [], 'boxcar');
%                                     [map,x] = fast1Ddelaymap(varX(tidxAllphs),spktrain(tidxAllphs),dx,samplerate,nbXbinsmth,visualspeed(tidxAllphs),20);
%                                     y = 0:20;
% %                                     map = map./repmat(max(map,[],1),size(map,1),1);
%                                     nplot = nplot + 1;
%                                     PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
%                                     PlotVar.Plots{nplot}.palette = 'hot';%'hot';%
%                                     if max(map(:)) > 0
%                                         PlotVar.Plots{nplot}.PlotMatrix(x, (y), map, [], 'Ylabel',num2str(fliplr(y)),'Clim',[0 max(map(:))]);
%                                     end
%                                 end
                            case 'Xpos'
                                iexp = 1;
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                tuningCurve = EXP.data.VStuning(iexp).zscore(icell,1:end-1);%EXP.data.VStuning(iexp).mean(icell,1:end-1);
                                errorset = EXP.data.VStuning(iexp).sem(icell,1:end-1);
                                w = gausswin(5);
                                PlotVar.Plots{nplot}.PlotVector(EXP.data.VStuning(iexp).Xpos, conv(tuningCurve,w,'same'), errorset, true, 'Ylim', [1.1*min(tuningCurve - errorset)+1*(max(tuningCurve  + errorset)==0) 1.1*max(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0)]);
%                                 set(gca,'Ylim', [-1.1*max(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0) 1.1*max(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0)]); 
                            case 'Ypos'
                                iexp = 2;
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                tuningCurve = EXP.data.VStuning(iexp).zscore(icell,1:end-1);%EXP.data.VStuning(iexp).mean(icell,1:end-1);
                                errorset = EXP.data.VStuning(iexp).sem(icell,1:end-1);
                                w = gausswin(5);
                                PlotVar.Plots{nplot}.PlotVector(EXP.data.VStuning(iexp).Ypos, conv(tuningCurve,w,'same'), errorset, true, 'Ylim', [1.1*min(tuningCurve - errorset)+1*(max(tuningCurve  + errorset)==0) 1.1*max(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0)]);
%                                 set(gca,'Ylim', [0 1.1*max(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0)]); 
                            case 'Ori'
                                iexp = 3;
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                tuningCurve = EXP.data.VStuning(iexp).mean(icell,1:end-1);
                                errorset = EXP.data.VStuning(iexp).sem(icell,1:end-1);
                                w = gausswin(5);
                                PlotVar.Plots{nplot}.PlotVector(0:numel(tuningCurve)-1, conv(tuningCurve,w,'same'), errorset, true, 'Ylim', [1.1*min(tuningCurve - errorset)+1*(max(tuningCurve  + errorset)==0) 1.1*max(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0)]);
%                                 set(gca,'Ylim', [0 1.1*max(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0)]); 
                            case 'Xpos x Ypos'
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                tuningCurveX = EXP.data.VStuning(1).meanNorm(icell,1:end-1);
                                tuningCurveY = EXP.data.VStuning(2).meanNorm(icell,1:end-1);
                                XYmap = tuningCurveY'*tuningCurveX;
                                w = gausswin(5);
                                PlotVar.Plots{nplot}.PlotMatrix(0:numel(tuningCurveX)-1, 0:numel(tuningCurveY)-1, interp2(conv2(XYmap,w*w','same'),2), [], 'YDir', 'normal');
                                PlotVar.Plots{nplot}.palette = 'hot';%'RedWhiteBlue';%'parula';                                
                        end
                    end
                end
            end
        end
    end
end
end