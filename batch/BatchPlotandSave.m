function BatchPlotandSave(PlotVarDecoder,Layout,EXP,objStr,foldername,figname)
    dirname = foldername;
    if ~isdir(dirname)
        mkdir(dirname);
    end
    
    if ~strcmp(objStr,'Mean X-Error')
        PlotVarDecoder.ChosenObj = {objStr};
        PlotVarDecoder.Fdisplaylog = false;
        PlotVarDecoder.Foverlap = false;
        PlotVarDecoder.FdisplayMat = true;
        PlotVarDecoder.FdisplayPredMax = false;
        PlotVarDecoder.FdisplayPredAve = true;
        PlotVarDecoder.Fnormalize = false;
        PlotVarDecoder.Clim = [0.5 1.5];
        PlotVarDecoder.Palettename = 'RedWhiteBlue';
        PlotVarDecoder.ChosenProbe = 1:numel(EXP.Bayes.Posterior0);
        
        PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
        Layout.fig2pdf([dirname filesep EXP.animal '_' num2str(EXP.iseries) '_' PlotVarDecoder.ChosenObj{1} '_' figname '_mat.pdf']);
        

%         PlotVarDecoder.ChosenObj = {objStr};
%         PlotVarDecoder.Fdisplaylog = false;
%         PlotVarDecoder.Foverlap = false;
%         PlotVarDecoder.FdisplayMat = true;
%         PlotVarDecoder.FdisplayPredMax = false;
%         PlotVarDecoder.FdisplayPredAve = true;
%         PlotVarDecoder.Fnormalize = true;
%         PlotVarDecoder.Clim = [0 1];
%         PlotVarDecoder.Palettename = 'parula';    
%         if ~Fsplit
%             cont = PlotVarDecoder.ChosenContrast;
%             for c = 1:numel(cont)
%                 PlotVarDecoder.ChosenContrast = cont(c);
%                 if EXP.SubsetVal.contrast(cont(c)) ~= 0
%                     PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);            
%                     Layout.fig2pdf([dirname filesep EXP.animal '_' num2str(EXP.iseries) '_' PlotVarDecoder.ChosenObj{1} '_cont' num2str(c) '_matNorm.pdf']);
%                 end
%             end
%             PlotVarDecoder.ChosenContrast = cont;
%         else
%             gain = PlotVarDecoder.ChosenGain;
%             cont = PlotVarDecoder.ChosenContrast;
%             for g = 1:numel(gain)
%                 for c = 1:numel(cont)
%                     PlotVarDecoder.ChosenGain = gain(g);
%                     PlotVarDecoder.ChosenContrast = cont(c);
%                     if EXP.SubsetVal.contrast(cont(c)) ~= 0
%                         PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
%                         Layout.fig2pdf([dirname filesep EXP.animal '_' num2str(EXP.iseries) '_' PlotVarDecoder.ChosenObj{1} '_' num2str(g) '_' num2str(c) '_matNorm.pdf']);
%                     end
%                 end
%             end
%             PlotVarDecoder.ChosenGain = gain;
%             PlotVarDecoder.ChosenContrast = cont;
%         end
    end

    PlotVarDecoder.ChosenObj = {objStr};
    PlotVarDecoder.Fdisplaylog = false;
    PlotVarDecoder.Foverlap = true;
    PlotVarDecoder.FdisplayMat = false;
    PlotVarDecoder.FdisplayPredMax = false;
    PlotVarDecoder.FdisplayPredAve = true;
    PlotVarDecoder.Fnormalize = true;
    PlotVarDecoder.Clim = [0.5 1.5];
    PlotVarDecoder.Palettename = 'RedWhiteBlue';
    PlotVarDecoder.Foverlap = true;
    PlotVarDecoder.ChosenProbe = 1:numel(EXP.Bayes.Posterior0);
    
    PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
    Layout.fig2pdf([dirname filesep EXP.animal '_' num2str(EXP.iseries) '_' PlotVarDecoder.ChosenObj{1} '_' figname '_pred.pdf']);
end