function PlotVar = UpdateplotPopSpks(PlotVar,EXP,Layout)
PlotVar.Plots = [];
roomlength = 100;
page = 3;
win = 1;
Layout.DivideWindow(page, win, size(PlotVar.ChosenObj,1)*size(PlotVar.ChosenContrast,2), size(PlotVar.ChosenGain,2)*size(PlotVar.ChosenRoomlength,2)*size(PlotVar.ChosenOutcome,2));

cbase = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
gbase = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
rbase = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
obase = find(EXP.SubsetVal.outcome == 2);

SSIrangeidx = find(EXP.CellInfo.SSI{cbase,gbase,rbase,obase} >= PlotVar.SSIrange(1) & EXP.CellInfo.SSI{cbase,gbase,rbase,obase} <= PlotVar.SSIrange(2));
raterangeidx = find(EXP.CellInfo.fieldMin{cbase,gbase,rbase,obase} >= PlotVar.Raterange(1) & EXP.CellInfo.fieldMin{cbase,gbase,rbase,obase} <= PlotVar.Raterange(2));
Probeidx = find(ismember(EXP.CellInfo.Probe, PlotVar.Probes));
goodunits = find(EXP.CellInfo.Goodcluster | EXP.CellInfo.Unsortedcluster);
interneurons = find(EXP.CellInfo.Finterneuron);
[maxZ,~] = max(EXP.CellInfo.fieldZ{cbase,gbase,rbase,obase},[],2);
goodplacefields = find(EXP.CellInfo.fieldPosSE{cbase,gbase,rbase,obase} < roomlength & abs(maxZ') > PlotVar.zth);
goodplacefields = 1:size(EXP.CellInfo.fieldPosSE{cbase,gbase,rbase,obase},2);

distri_Th = 0.99;
% SinfoperSpike = EXP.CellInfo.SpatialInfoPerSpike{cbase, gbase, rbase, obase};
% SinfoperSpikeZscore = (SinfoperSpike - quantile(EXP.CellInfo.SpatialInfoPerSpikeRef{cbase, gbase, rbase, obase},distri_Th,2)')./quantile(EXP.CellInfo.SpatialInfoPerSpikeRef{cbase, gbase, rbase, obase},distri_Th,2)';
% signiSinfoperSpike = find(SinfoperSpikeZscore>0);


[~,cellidx] = sort(EXP.CellInfo.fieldPos{cbase,gbase,rbase,obase},'ascend');
% cellidx = cellidx(ismember(cellidx,SSIrangeidx) & ismember(cellidx,raterangeidx) & ismember(cellidx,goodunits) & ismember(cellidx,goodplacefields) & ismember(cellidx,Probeidx));
cellidx = cellidx(ismember(cellidx,goodunits) & ~ismember(cellidx,interneurons));

colorcond = cell(numel(EXP.SubsetVal.gain),numel(EXP.SubsetVal.contrast));
for g = 1:numel(EXP.SubsetVal.gain)
    for c = 1:numel(EXP.SubsetVal.contrast)
        rgb = zeros(1,3);
        rgb(g) = 1 * c/numel(EXP.SubsetVal.contrast);
        colorcond{g,c} = rgb;
    end
end
cidx = find(ismember(EXP.SubsetVal.contrast,[0.5 0.7]));
try
for c = 1:numel(cidx)
    colorcond{EXP.SubsetVal.gain == 0.4,cidx(c)} = 'c';
    colorcond{EXP.SubsetVal.gain == 0.5,cidx(c)} = [0.5 0.5 0.5];
    colorcond{EXP.SubsetVal.gain == 0.6,cidx(c)} = 'm';
end
catch
end

nplot =  0;
for c = 1:size(PlotVar.ChosenContrast,2)
    for g = 1:size(PlotVar.ChosenGain,2)
        for r = 1:size(PlotVar.ChosenRoomlength,2)
            for o = 1:size(PlotVar.ChosenOutcome,2)
                tidx = EXP.getSubsets(PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o));
                
                exptidx = 0;
                seriesid = unique(EXP.data.es.series);
                for sid = 1:numel(seriesid)
                    exptidx = exptidx + (EXP.data.es.series == seriesid(sid) & ismember(EXP.data.es.iexp,PlotVar.explist(sid,:)));
                end
                
                tidx = tidx & exptidx;
                
                for k = 1:size(PlotVar.ChosenObj,1)
                    switch PlotVar.ChosenObj{k}
                        case 'fields'
                            nplot = nplot + 1;
                            PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (c-1)*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});                            
                            mat = EXP.CellInfo.field{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx,:)';
                            matCOMref = EXP.CellInfo.fieldCOM{cbase,gbase,rbase,obase}(cellidx)';
                            matPosref = EXP.CellInfo.fieldPos{cbase,gbase,rbase,obase}(cellidx)';
                            matPosSE = EXP.CellInfo.fieldPosSE{cbase,gbase,rbase,obase}(cellidx)';
%                             matZ = EXP.CellInfo.fieldZ{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx,:)';
%                             mat(abs(matZ)<1.96) = 0;                            
                            if PlotVar.FNormalize
                                mat = mat./repmat(max(mat),[size(mat,1) 1]);
                            end
                            
                            
                            matPos = matPosref;
                            
%                             matPos = zeros(size(matPosref));
%                             for icell = 1:numel(cellidx)
%                                 mat(:,icell) = circshift(mat(:,icell),-round(matPosref(icell)+size(mat,1)/2));
%                                 matPos(icell) = round(size(mat,1)/2);
%                             end
                            
%                             PlotVar.Plots{nplot}.PlotVector([],mean(mat,2), [], [], true, 'color','k','linewidth', 2,'Xlim',[1 size(mat,1)],'Ylim',[0 1]);
                            
                            PlotVar.Plots{nplot}.palette = 'parula';
                            PlotVar.Plots{nplot}.PlotMatrix(linspace(1,roomlength,size(mat,1)), [], mat', [], true, 'PlotBoxAspectRatio', [1 1 1], 'Clim', PlotVar.Clim);
                            
                                                        
                            PlotVar.Plots{nplot}.PlotVector(matPos, 1:numel(matPos), [], [], true, 'color',[1 1 1],'linewidth', 2,'Xlim',[1 size(mat,1)]);
                            PlotVar.Plots{nplot}.PlotVector(matCOMref, 1:numel(matCOMref), [], [], true, 'color',[0 0 0],'linewidth', 2,'Xlim',[1 size(mat,1)]);
                            
%                             PlotVar.Plots{nplot}.PlotVector(matPos + matPosSE, 1:numel(matPos), [], [], true, 'color',[0.5 0.5 0.5],'linewidth', 2,'Xlim',[1 size(mat,1)]);
%                             PlotVar.Plots{nplot}.PlotVector(matPos - matPosSE, 1:numel(matPos), [], [], true, 'color',[0.5 0.5 0.5],'linewidth', 2,'Xlim',[1 size(mat,1)]);
                            
                        case 'fieldZ'
                            nplot = nplot + 1;
                            PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (c-1)*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});                            
                            mat = EXP.CellInfo.field{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx,:)';
                            if PlotVar.FNormalize
                                mat = mat./repmat(max(mat),[size(mat,1) 1]);
                            end
                            PlotVar.Plots{nplot}.palette = 'parula';
                            PlotVar.Plots{nplot}.PlotMatrix(linspace(1,roomlength,size(mat,1)), [], mat', [], true, 'PlotBoxAspectRatio', [1 1 1], 'Clim', PlotVar.Clim);
                            
                            matPos = EXP.CellInfo.fieldPos{cbase,gbase,rbase,obase}(cellidx)';
                            matPosSE = EXP.CellInfo.fieldPosSE{cbase,gbase,rbase,obase}(cellidx)';
                            PlotVar.Plots{nplot}.PlotVector(matPos, 1:numel(matPos), [], [], true, 'color',[1 1 1],'linewidth', 2);
                            PlotVar.Plots{nplot}.PlotVector(matPos + matPosSE, 1:numel(matPos), [], [], true, 'color',[0.5 0.5 0.5],'linewidth', 2);
                            PlotVar.Plots{nplot}.PlotVector(matPos - matPosSE, 1:numel(matPos), [], [], true, 'color',[0.5 0.5 0.5],'linewidth', 2);
                            
                        case 'theta precession'
                            nplot = nplot + 1;
                            PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (c-1)*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});                            
                            mat = EXP.CellInfo.field{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx,:)';
                            matPosref = EXP.CellInfo.fieldPos{cbase,gbase,rbase,obase}(cellidx)';
                            matCOMref = EXP.CellInfo.fieldCOM{cbase,gbase,rbase,obase}(cellidx)';
                            matPosSE = EXP.CellInfo.fieldPosSE{cbase,gbase,rbase,obase}(cellidx)';
%                             matZ = EXP.CellInfo.fieldZ{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx,:)';
%                             mat(abs(matZ)<1.96) = 0;                            
                            if PlotVar.FNormalize
                                mat = mat./repmat(max(mat),[size(mat,1) 1]);
                            end
                            
                            for icell = 1:size(mat,2)
                                mat(:,icell) = circshift(mat(:,icell),50 - floor(matPosref(icell))-1);
                            end
                            
                            matPos = matCOMref;
                            
%                             matPos = zeros(size(matPosref));
%                             for icell = 1:numel(cellidx)
%                                 mat(:,icell) = circshift(mat(:,icell),-round(matPosref(icell)+size(mat,1)/2));
%                                 matPos(icell) = round(size(mat,1)/2);
%                             end
                            
%                             PlotVar.Plots{nplot}.PlotVector([],mean(mat,2), [], [], true, 'color','k','linewidth', 2,'Xlim',[1 size(mat,1)],'Ylim',[0 1]);
                            
                            PlotVar.Plots{nplot}.palette = 'parula';
                            PlotVar.Plots{nplot}.PlotMatrix(linspace(1,roomlength,size(mat,1)), [], mat', [], true, 'PlotBoxAspectRatio', [1 1 1], 'Clim', PlotVar.Clim);
%                             PlotVar.Plots{nplot}.PlotVector([],sum(mat,2), [], [], true, 'color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                            
                                                        
%                             PlotVar.Plots{nplot}.PlotVector(matPos, 1:numel(matPos), [], [], true, 'color',[1 1 1],'linewidth', 2,'Xlim',[1 size(mat,1)]);
                            PlotVar.Plots{nplot}.PlotVector(50*ones(1,numel(matPos)), 1:numel(matPos), [], [], true, 'color',[1 1 1],'linewidth', 2,'Xlim',[1 size(mat,1)]);
                            
%                             PlotVar.Plots{nplot}.PlotVector(matPos + matPosSE, 1:numel(matPos), [], [], true, 'color',[0.5 0.5 0.5],'linewidth', 2,'Xlim',[1 size(mat,1)]);
%                             PlotVar.Plots{nplot}.PlotVector(matPos - matPosSE, 1:numel(matPos), [], [], true, 'color',[0.5 0.5 0.5],'linewidth', 2,'Xlim',[1 size(mat,1)]);
                            
                            
                        case 'theta precession '                            
                            nplot = nplot + 1;
                            PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (c-1)*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            h = histogram(EXP.CellInfo.PhaseMax{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx), floor(roomlength/5));
                            
                            PlotVar.Plots{nplot}.PlotVector(h.BinEdges(1:end-1), h.Values, [], [], true);  
                        case 'distri fieldpos'
                            nplot = nplot + 1;
                            PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (c-1)*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            h = histogram(EXP.CellInfo.fieldPos{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx), floor(roomlength/5));
                            
                            PlotVar.Plots{nplot}.PlotVector(h.BinEdges(1:end-1), h.Values, [], [], true);  
                            
                        case 'distri SSI'
                            nplot = nplot + 1;
                            PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (c-1)*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
%                             [~,imax] = max(EXP.CellInfo.fieldXcorr{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx,:),[],2);
%                             h = histogram(imax-51,-50.5:50.5);
%                             PlotVar.Plots{nplot}.PlotVector([0 0], [0 max(h.Values)], [], [], true, 'Xlim', [-20 20],'linewidth',2); 
%                             PlotVar.Plots{nplot}.PlotVector(h.BinEdges(1:end-1), h.Values, [], [], true, 'Xlim', [-20 20]); 
                            h = histogram(EXP.CellInfo.SSI{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx), floor(roomlength/5));
                            PlotVar.Plots{nplot}.PlotVector(h.BinEdges(1:end-1), h.Values, [], [], true); 
                        case 'distri X shift'
                            nplot = nplot + 1;
                            PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (c-1)*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            [~,imax] = max(EXP.CellInfo.fieldXcorr{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx,:),[],2);
                            h = histogram(imax-51,-50.5:50.5);
                            PlotVar.Plots{nplot}.PlotVector([0 0], [0 max(h.Values)], [], [], true, 'Xlim', [-20 20],'linewidth',2); 
                            PlotVar.Plots{nplot}.PlotVector(h.BinEdges(1:end-1), h.Values, [], [], true, 'Xlim', [-20 20]);
                        case 'distri max rate'
                            nplot = nplot + 1;
                            PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (c-1)*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            h = histogram(EXP.CellInfo.rate{PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o)}(cellidx), floor(roomlength/5));
                            
                            PlotVar.Plots{nplot}.PlotVector(h.BinEdges(1:end-1), h.Values, [], [], true); 
                        case 'distri theta coupling'
                            nplot = nplot + 1;
                            PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, (c-1)*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                            XYmap = EXP.data.VStuning(2).meanNorm(cellidx,1:end-1)'*EXP.data.VStuning(1).meanNorm(cellidx,1:end-1);%EXP.data.VStuning(4).mean(cellidx,1:end-1)'*EXP.data.VStuning(3).mean(cellidx,1:end-1);
                            w = gausswin(5);
                            PlotVar.Plots{nplot}.PlotMatrix(0:size(XYmap,2)-1, 0:size(XYmap,1)-1, interp2(conv2(XYmap,w*w','same'),2), []);
                            PlotVar.Plots{nplot}.palette = 'hot';%'RedWhiteBlue';%'parula';
                    end
                end
            end
        end
    end
end

end
