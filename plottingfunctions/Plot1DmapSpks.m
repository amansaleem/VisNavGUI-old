function Plot1DmapSpks(dataplot,map1D,icell)
if ~isempty(map1D.model)
    sampleRate = map1D.sampleRate;
    meanresp = map1D.model.tuning(icell).meanrespModel;
    stdresp = map1D.model.tuning(icell).SErespModel;
%     kfold = size(map1D.model.tuning(icell).respModel,1);
%     for i = 1:kfold
%         stdresp = stdresp + (kfold - 1)/kfold*(map1D.model.tuning(icell).respModel(i,:) - meanresp).^2;
%     end
%     stdresp = stdresp.^0.5;
% %     stdresp = std(map1D.model.tuning(icell).respModel,1);
    tuningCurve = meanresp;%.*sampleRate;
    errorset    = stdresp;%.*sampleRate;
    
    dataplot.PlotVector([], tuningCurve, errorset, true, 'Ylim', [0 1.1*max(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0)]);
    set(gca,'Ylim', [min(0,1.1*min(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0)) 1.1*max(tuningCurve + errorset)+1*(max(tuningCurve  + errorset)==0)]);
    
    xCOM = map1D.model.tuning(icell).meanrespModelXpos;
    xmax = map1D.model.tuning(icell).meanrespModelXmax;
    hold on;plot([xCOM xCOM],[0 max(tuningCurve)]);
    hold on;plot([xmax xmax],[0 max(tuningCurve)]);
    
    m = mean(map1D.model.tuning(icell).meanrespModel);
    mm = mean(map1D.model.tuning(icell).meanrespModel(map1D.model.tuning(icell).meanrespModel<m));
    hold on;plot([1 numel(map1D.model.tuning(icell).meanrespModel)],[m m]);

%     dataplot.PlotVector(map1D.bins, tuningCurve./errorset, true);    
%     set(gca,'Ylim', [0 1.1*max(tuningCurve./errorset)+1*(max(tuningCurve./errorset)==0)]);

%     zth = 1.96;
%     if max(meanresp) >= min(meanresp) + zth*stdresp(meanresp == max(meanresp))
%         dataplot.PlotVector(map1D.bins, tuningCurve, errorset, true, [], 'color', 'g');
%     else
%         dataplot.PlotVector(map1D.bins, tuningCurve, errorset, true, [], 'color', 'r');
%     end
end
end