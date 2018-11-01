function VisNavGUI
%path of the directories where your data are
SetDirs;
%create the main data object
EXP = TVRData;

%create fields of the parameter structures used in the GUI
[PlotVar1DMaps, PlotVarDecoder, PlotVarPop, PlotVarBehavior] = DefinePlotVariable;

%create the GUI windows
Layout = TMultigraph('VisNav');
Layout.RenamePage(1,'Spike1DMaps');
Layout.addPage('Decoding');
Layout.addPage('Population');
Layout.addPage('Behavior');

Layout.Hdividepage(1, 2, [-9 -1]);
Layout.Hdividepage(2, 2, [-9 -1]);
Layout.Hdividepage(3, 2, [-9 -1]);
Layout.Hdividepage(4, 2, [-9 -1]);
Layout.addPermWindow('Vertical',[-1 -9]);

CellDialog0 = Tdialog(Layout.permwindow{1});
SpkMapsDialog = Tdialog(Layout.window{1,2});
DecodingDialog = Tdialog(Layout.window{2,2});
PopDialog = Tdialog(Layout.window{3,2});
BehaviorDialog = Tdialog(Layout.window{4,2});

%creating the GUI dialogs
pagenum = 0;%page 0 is always displayed, regardless of the panel that is selected
winnum = 1;
relsize = [];
CellDialog0.getpushbutton('Animal',@(hObject, eventdata)SelectAnimal_callback(hObject, eventdata, DIRS, CellDialog0, EXP)); relsize = [relsize 1];
CellDialog0.getpushbutton('Series',@(hObject, eventdata)SelectSeries_callback(hObject, eventdata, DIRS, CellDialog0, EXP)); relsize = [relsize 1];
CellDialog0.getvariable('Shank #',@(hObject, eventdata)GetShanknum_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.Shanknum); relsize = [relsize 0.5 0.5];
CellDialog0.getvariable('Shank suffix',@(hObject, eventdata)GetShanksuffix_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),strjoin(PlotVar1DMaps.Shanksuffix)); relsize = [relsize 0.5 0.5];
CellDialog0.getvariable('Speed thresh.',@(hObject, eventdata)GetSpeedThreshold_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.SpeedThreshold); relsize = [relsize 0.5 0.5];
CellDialog0.getvariable('# theta phs bins',@(hObject, eventdata)GetNthetaphsbins_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.Nthetaphsbins); relsize = [relsize 0.5 0.5];
CellDialog0.getvariable('Time smth',@(hObject, eventdata)GetSmthTimeWindow_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.SmthTimeWindow); relsize = [relsize 0.5 0.5];
CellDialog0.getvariable('Spatial smth',@(hObject, eventdata)GetSmthSpatialWindow_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.SmthSpatialWindow); relsize = [relsize 0.5 0.5];
CellDialog0.getpushbutton('Exp',@(hObject, eventdata)SelectExpt_callback(hObject, eventdata, DIRS, CellDialog0, SpkMapsDialog, BehaviorDialog, PopDialog, EXP, PlotVarDecoder, PlotVar1DMaps, PlotVarBehavior, PlotVarPop)); relsize = [relsize 1];
CellDialog0.getpushbutton('Load Saved file',@(hObject, eventdata)LoadSavedExpt_callback(hObject, eventdata, DIRS, CellDialog0, SpkMapsDialog, BehaviorDialog, PopDialog, DecodingDialog, EXP, PlotVarDecoder, PlotVar1DMaps, PlotVarBehavior, PlotVarPop)); relsize = [relsize 1];
CellDialog0.getcheckbox('Load Saved maps',@(hObject, eventdata)FloadfromfileCellinfo_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.FloadfromfileCellinfo); relsize = [relsize 1];
CellDialog0.getpushbutton('Load multiple Series',@(hObject, eventdata)SelectMultipleSeries_callback(hObject, eventdata, DIRS, CellDialog0, SpkMapsDialog, BehaviorDialog, PopDialog, EXP, PlotVar1DMaps, PlotVarBehavior, PlotVarPop, DecodingDialog, PlotVarDecoder,Layout)); relsize = [relsize 1];
CellDialog0.getpushbutton('Compute 1Dmap',@(hObject, eventdata)Compute1Dmap_callback(hObject, eventdata, EXP, PlotVar1DMaps)); relsize = [relsize 1];
CellDialog0.getpushbutton('Bayesian Decoder',@(hObject, eventdata)BayesDecoder_callback(hObject, eventdata, EXP, PlotVarDecoder, DecodingDialog, Layout)); relsize = [relsize 1];
CellDialog0.getpushbutton('Batch Bayes Decoder',@(hObject, eventdata)BatchBayesDecoder_callback(hObject, eventdata, DIRS, EXP, PlotVarDecoder, DecodingDialog, Layout)); relsize = [relsize 1];
CellDialog0.getpushbutton('Search decoder Opt. Params',@(hObject, eventdata)BatchBayesDecoderParameters_callback(hObject, eventdata, DIRS, EXP, PlotVarDecoder, DecodingDialog, Layout)); relsize = [relsize 1];
CellDialog0.getpushbutton('Batch decoder Single cells',@(hObject, eventdata)BatchBayesDecoderSingleCells_callback(hObject, eventdata, DIRS, EXP, PlotVarDecoder, DecodingDialog, Layout)); relsize = [relsize 1];
CellDialog0.getlistbox('CA1 - cell#',{'null'},1,1,@(hObject, eventdata)SelectCellCA1_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 0.5 12];
CellDialog0.getlistbox('V1 - cell#',{'null'},1,1,@(hObject, eventdata)SelectCellV1_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 0.5 8];
CellDialog0.getcheckbox('Pool cells',@(hObject, eventdata)FPoolCell_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.FpoolCell); relsize = [relsize 1];

Layout.SetWindowLayout(pagenum, winnum, -relsize, 'Vertical');


pagenum = 1;
winnum = 2;
relsize = [];
SpkMapsDialog.getpushbutton('Series/Exp selection',@(hObject, eventdata)GetSeriesrange_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 1];
SpkMapsDialog.getlistbox('X variable',{'null'},1,1,@(hObject, eventdata)SelectVarX_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 0.3 2];
SpkMapsDialog.getlistbox('Y variable',{'null'},1,1,@(hObject, eventdata)SelectVarY_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 0.3 2];
SpkMapsDialog.getlistbox('Plots',{'null'},1,1,@(hObject, eventdata)SelectPlotObj_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 0.3 2];
SpkMapsDialog.getcheckbox('Disp. Mat.',@(hObject, eventdata)Fdispmat_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.Fdispmat); relsize = [relsize 0.3];
SpkMapsDialog.getvariable('Speed Thresh.',@(hObject, eventdata)Getspeedth_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.speed_th); relsize = [relsize 0.3 0.3];
SpkMapsDialog.getvariable('delayT',@(hObject, eventdata)GetdelayT_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.delayT); relsize = [relsize 0.3 0.3];
SpkMapsDialog.addslider(-100,100,1,5,PlotVar1DMaps.delayT); relsize = [relsize 0.2];
SpkMapsDialog.getvariable('Xbin size',@(hObject, eventdata)GetXbinsize_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.Xbinsize); relsize = [relsize 0.3 0.3];
SpkMapsDialog.addslider(0,100,1,2.5,PlotVar1DMaps.Xbinsize); relsize = [relsize 0.2];
SpkMapsDialog.getvariable('Theta Channel',@(hObject, eventdata)GetthetaChannel_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.thetaChannel); relsize = [relsize 0.3 0.3];
SpkMapsDialog.addslider(1,37,1,4,PlotVar1DMaps.thetaChannel); relsize = [relsize 0.2];
SpkMapsDialog.getcheckbox('Combine conditions',@(hObject, eventdata)FCombinedCond_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.FXcond); relsize = [relsize 0.3];
SpkMapsDialog.getlistbox('Contrast',{'null'},1,1,@(hObject, eventdata)SelectContrast_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 0.3 2];
SpkMapsDialog.getcheckbox('Pool Contrasts',@(hObject, eventdata)FPoolContrast_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.FpoolContrast); relsize = [relsize 0.5];
SpkMapsDialog.getlistbox('Gain',{'null'},1,1,@(hObject, eventdata)SelectGain_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 0.3 2];
SpkMapsDialog.getcheckbox('Pool Gains',@(hObject, eventdata)FPoolGain_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.FpoolGain); relsize = [relsize 0.5];
SpkMapsDialog.getlistbox('Roomlength',{'null'},1,1,@(hObject, eventdata)SelectRoomlength_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 0.3 2];
SpkMapsDialog.getcheckbox('Pool Roomlengths',@(hObject, eventdata)FPoolRoomlength_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.FpoolRoomlength); relsize = [relsize 0.5];
SpkMapsDialog.getlistbox('Outcome',{'null'},1,1,@(hObject, eventdata)SelectOutcome_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)); relsize = [relsize 0.3 2];
SpkMapsDialog.getcheckbox('Pool Outcomes',@(hObject, eventdata)FPoolOutcome_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout),PlotVar1DMaps.FpoolOutcome); relsize = [relsize 0.5];

Layout.SetWindowLayout(pagenum, winnum, -relsize, 'Vertical');

pagenum = 2;
winnum = 2;
relsize = [];
DecodingDialog.getpushbutton('Series/Exp selection',@(hObject, eventdata)GetDecSeriesrange_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)); relsize = [relsize 1];
DecodingDialog.getlistbox('Probe',{'CA1','V1'},0,2,@(hObject, eventdata)GetChosenProbe_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)); relsize = [relsize 0.5 1.5];
DecodingDialog.getvariable('Tolerance for Max.',@(hObject, eventdata)Getmaxtolerance_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.maxtolerance); relsize = [relsize 0.5 0.5];
DecodingDialog.getcheckbox('use Xdec distri',@(hObject, eventdata)FdecXdistribution_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FdecXdistribution); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('only good time bins',@(hObject, eventdata)Fgoodtimebins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.Fgoodtimebins); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('Smooth spatiallly',@(hObject, eventdata)Fspatialsmooth_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.Fspatialsmooth); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('Use Theta Post.',@(hObject, eventdata)FthetaPost_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FthetaPost); relsize = [relsize 0.5];
DecodingDialog.getvariable('Theta Post Ref.',@(hObject, eventdata)GetthetaDecphase_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.thetaDecphase); relsize = [relsize 0.5 0.5];
DecodingDialog.addslider(0,360,1,30,PlotVarDecoder.thetaDecphase); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('Use Theta bins',@(hObject, eventdata)Fthetabins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.Fthetabins); relsize = [relsize 0.5];
DecodingDialog.getvariable('Theta phase',@(hObject, eventdata)Getdecthetaphase_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.thetaphase); relsize = [relsize 0.5 0.5];
DecodingDialog.addslider(0,360,1,30,PlotVarDecoder.thetaphase); relsize = [relsize 0.5];
DecodingDialog.getvariable('# Theta phase bins',@(hObject, eventdata)Getdecnthetaphsbins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.nthetaphsbins); relsize = [relsize 0.5 0.5];
DecodingDialog.addslider(1,24,1,2,PlotVarDecoder.nthetaphsbins); relsize = [relsize 0.5];
DecodingDialog.getvariable('Theta Channel',@(hObject, eventdata)GetdecthetaChannel_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.thetaChannel); relsize = [relsize 0.5 0.5];
DecodingDialog.addslider(1,37,1,4,PlotVarDecoder.thetaChannel); relsize = [relsize 0.5];
DecodingDialog.getvariable('# Theta X bins',@(hObject, eventdata)GetdecnthetaXbins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.nthetaXbins); relsize = [relsize 0.5 0.5];
DecodingDialog.addslider(1,100,1,2,PlotVarDecoder.nthetaXbins); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('Use Speed bins',@(hObject, eventdata)Fspdbins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.Fspdbins); relsize = [relsize 0.5];
DecodingDialog.getvariable('Speed bin',@(hObject, eventdata)Getspdbin_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.spdbin); relsize = [relsize 0.5 0.5];
DecodingDialog.addslider(1,10,1,2,PlotVarDecoder.spdbin); relsize = [relsize 0.5];
DecodingDialog.getlistbox('Plots',{'null'},1,1,@(hObject, eventdata)SelectDecPlotObj_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)); relsize = [relsize 0.5 3];
DecodingDialog.getlistbox('Contrast',{'null'},1,1,@(hObject, eventdata)SelectDecContrast_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)); relsize = [relsize 0.5 2];
DecodingDialog.getcheckbox('Pool Contrasts',@(hObject, eventdata)FPoolDecContrast_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FpoolContrast); relsize = [relsize 0.5];
DecodingDialog.getlistbox('Gain',{'null'},1,1,@(hObject, eventdata)SelectDecGain_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)); relsize = [relsize 0.5 2];
DecodingDialog.getcheckbox('Pool Gain',@(hObject, eventdata)FPoolDecGain_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FpoolGain); relsize = [relsize 0.5];
DecodingDialog.getlistbox('Roomlength',{'null'},1,1,@(hObject, eventdata)SelectDecRoomlength_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)); relsize = [relsize 0.5 1];
DecodingDialog.getcheckbox('Pool Roomlength',@(hObject, eventdata)FPoolDecRoomlength_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FpoolRoomlength); relsize = [relsize 0.5];
DecodingDialog.getlistbox('Outcome',{'null'},1,1,@(hObject, eventdata)SelectDecOutcome_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)); relsize = [relsize 0.5 2];
DecodingDialog.getcheckbox('Pool Outcome',@(hObject, eventdata)FPoolDecOutcome_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FpoolOutcome); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('Display Log',@(hObject, eventdata)Fdisplaylog_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.Fdisplaylog); relsize = [relsize 0.5];
DecodingDialog.getvariable('Climits',@(hObject, eventdata)GetClim_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.Clim); relsize = [relsize 0.5 0.5];
DecodingDialog.getcheckbox('Normalize',@(hObject, eventdata)Fnormalize_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.Fnormalize); relsize = [relsize 0.5];
DecodingDialog.getpopupmenu('Palette',{'RedWhiteBlue','parula','hot','hotcoldlin'},1,@(hObject, eventdata)GetPalettename_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)); relsize = [relsize 0.5 0.7];
DecodingDialog.getcheckbox('overlap',@(hObject, eventdata)Foverlap_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.Foverlap); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('display PredMax',@(hObject, eventdata)FdisplayPredMax_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FdisplayPredMax); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('display PredAve',@(hObject, eventdata)FdisplayPredAve_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FdisplayPredAve); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('display Mat',@(hObject, eventdata)FdisplayMat_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FdisplayMat); relsize = [relsize 0.5];
DecodingDialog.getcheckbox('display LFP',@(hObject, eventdata)FdisplayLFP_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout),PlotVarDecoder.FdisplayLFP); relsize = [relsize 0.5];

Layout.SetWindowLayout(pagenum, winnum, -relsize, 'Vertical');

pagenum = 3;
winnum = 2;
relsize = [];
PopDialog.getpushbutton('Series/Exp selection',@(hObject, eventdata)GetPopSeriesrange_callback(hObject, eventdata, EXP, PlotVarPop, Layout)); relsize = [relsize 1];
PopDialog.getlistbox('Plots',{'null'},1,1,@(hObject, eventdata)SelectPopPlotObj_callback(hObject, eventdata, EXP, PlotVarPop, Layout)); relsize = [relsize 0.3 2];
PopDialog.getvariable('SSI range',@(hObject, eventdata)GetPopSSIrange_callback(hObject, eventdata, EXP, PlotVarPop, Layout),PlotVarPop.SSIrange); relsize = [relsize 0.3 0.3];
PopDialog.getvariable('rate range',@(hObject, eventdata)GetPopraterange_callback(hObject, eventdata, EXP, PlotVarPop, Layout),PlotVarPop.Raterange); relsize = [relsize 0.3 0.3];
PopDialog.getvariable('zth max >',@(hObject, eventdata)GetPopzth_callback(hObject, eventdata, EXP, PlotVarPop, Layout),PlotVarPop.zth); relsize = [relsize 0.3 0.3];
PopDialog.getlistbox('Probe',{'CA1','V1'},0,2,@(hObject, eventdata)GetPopProbe_callback(hObject, eventdata, EXP, PlotVarPop, Layout)); relsize = [relsize 0.5 1.5];
PopDialog.getlistbox('Contrast',{'null'},1,1,@(hObject, eventdata)SelectPopContrast_callback(hObject, eventdata, EXP, PlotVarPop, Layout)); relsize = [relsize 0.3 2];
PopDialog.getcheckbox('Pool Contrasts',@(hObject, eventdata)FPoolPopContrast_callback(hObject, eventdata, EXP, PlotVarPop, Layout),PlotVarPop.FpoolContrast); relsize = [relsize 0.5];
PopDialog.getlistbox('Gain',{'null'},1,1,@(hObject, eventdata)SelectPopGain_callback(hObject, eventdata, EXP, PlotVarPop, Layout)); relsize = [relsize 0.3 2];
PopDialog.getcheckbox('Pool Gains',@(hObject, eventdata)FPoolPopGain_callback(hObject, eventdata, EXP, PlotVarPop, Layout),PlotVarPop.FpoolGain); relsize = [relsize 0.5];
PopDialog.getlistbox('Roomlength',{'null'},1,1,@(hObject, eventdata)SelectPopRoomlength_callback(hObject, eventdata, EXP, PlotVarPop, Layout)); relsize = [relsize 0.3 2];
PopDialog.getcheckbox('Pool Roomlengths',@(hObject, eventdata)FPoolPopRoomlength_callback(hObject, eventdata, EXP, PlotVarPop, Layout),PlotVarPop.FpoolRoomlength); relsize = [relsize 0.5];
PopDialog.getlistbox('Outcome',{'null'},1,1,@(hObject, eventdata)SelectPopOutcome_callback(hObject, eventdata, EXP, PlotVarPop, Layout)); relsize = [relsize 0.3 2];
PopDialog.getcheckbox('Pool Outcomes',@(hObject, eventdata)FPoolPopOutcome_callback(hObject, eventdata, EXP, PlotVarPop, Layout),PlotVarPop.FpoolOutcome); relsize = [relsize 0.5];
PopDialog.getvariable('Climits',@(hObject, eventdata)GetPopClim_callback(hObject, eventdata, EXP, PlotVarPop, Layout),PlotVarPop.Clim); relsize = [relsize 0.3 0.3];
PopDialog.getcheckbox('Normalize',@(hObject, eventdata)PopFnormalize_callback(hObject, eventdata, EXP, PlotVarPop, Layout),PlotVarPop.FNormalize); relsize = [relsize 0.5];

Layout.SetWindowLayout(pagenum, winnum, -relsize, 'Vertical');

pagenum = 4;
winnum = 2;
relsize = [];
BehaviorDialog.getpushbutton('Series/Exp selection',@(hObject, eventdata)GetBehSeriesrange_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)); relsize = [relsize 1];
BehaviorDialog.getlistbox('Plots',{'null'},1,1,@(hObject, eventdata)SelectBehPlotObj_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)); relsize = [relsize 0.3 5];
BehaviorDialog.getvariable('Speed Thresh.',@(hObject, eventdata)GetBehspeedth_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout),PlotVarBehavior.speed_th); relsize = [relsize 0.3 0.3];
BehaviorDialog.getlistbox('Contrast',{'null'},1,1,@(hObject, eventdata)SelectBehContrast_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)); relsize = [relsize 0.3 3];
BehaviorDialog.getcheckbox('Pool Contrasts',@(hObject, eventdata)FPoolBehContrast_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout),PlotVarBehavior.FpoolContrast); relsize = [relsize 0.5];
BehaviorDialog.getlistbox('Gain',{'null'},1,1,@(hObject, eventdata)SelectBehGain_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)); relsize = [relsize 0.3 2];
BehaviorDialog.getcheckbox('Pool Gain',@(hObject, eventdata)FPoolBehGain_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout),PlotVarBehavior.FpoolGain); relsize = [relsize 0.5];
BehaviorDialog.getlistbox('Roomlength',{'null'},1,1,@(hObject, eventdata)SelectBehRoomlength_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)); relsize = [relsize 0.3 2];
BehaviorDialog.getcheckbox('Pool Roomlength',@(hObject, eventdata)FPoolBehRoomlength_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout),PlotVarBehavior.FpoolRoomlength); relsize = [relsize 0.5];
BehaviorDialog.getlistbox('Outcome',{'null'},1,1,@(hObject, eventdata)SelectBehOutcome_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)); relsize = [relsize 0.3 2];
BehaviorDialog.getcheckbox('Pool Outcome',@(hObject, eventdata)FPoolBehOutcome_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout),PlotVarBehavior.FpoolOutcome); relsize = [relsize 0.5];

Layout.SetWindowLayout(pagenum, winnum, -relsize, 'Vertical');

assignin('base','Layout',Layout);

end


% *************************Callback functions***************************** %
function GetSpeedThreshold_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = str2num(get(hObject,'String'));
PlotVar1DMaps.SpeedThreshold = val;
end

function GetNthetaphsbins_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = str2num(get(hObject,'String'));
PlotVar1DMaps.Nthetaphsbins = val;
end

function GetShanknum_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = str2num(get(hObject,'String'));
PlotVar1DMaps.Shanknum = val;
end

function GetShanksuffix_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
str = get(hObject,'String');
PlotVar1DMaps.Shanksuffix = strsplit(str);
end

function GetSmthTimeWindow_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = str2num(get(hObject,'String'));
PlotVar1DMaps.SmthTimeWindow = val;
end

function GetSmthSpatialWindow_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = str2num(get(hObject,'String'));
PlotVar1DMaps.SmthSpatialWindow = val;
end

function FloadfromfileCellinfo_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.FloadfromfileCellinfo = get(hObject,'Value');
end

function SelectCellCA1_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
sel = get(hObject,'Value');
PlotVar1DMaps.ChosenCell = sel;
if PlotVar1DMaps.FpoolCell
    PlotVar1DMaps.ChosenCell = PlotVar1DMaps.ChosenCell(:);
else
    PlotVar1DMaps.ChosenCell = PlotVar1DMaps.ChosenCell(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
uicontrol(hObject);
end

function SelectCellV1_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
sel = get(hObject,'Value');
PlotVar1DMaps.ChosenCell = sel + sum(EXP.CellInfo.Probe == 1);
if PlotVar1DMaps.FpoolCell
    PlotVar1DMaps.ChosenCell = PlotVar1DMaps.ChosenCell(:);
else
    PlotVar1DMaps.ChosenCell = PlotVar1DMaps.ChosenCell(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
uicontrol(hObject);
end


function FPoolCell_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.FpoolCell = get(hObject,'Value');
if PlotVar1DMaps.FpoolCell
    PlotVar1DMaps.ChosenCell = PlotVar1DMaps.ChosenCell(:);
else
    PlotVar1DMaps.ChosenCell = PlotVar1DMaps.ChosenCell(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function SelectAnimal_callback(hObject, eventdata, DIRS, dialog, EXP)
EXP.SelectAnimal(DIRS,{'BALL','JF','MIK'});
dialog.UpdateUIcontrol('Animal','String', EXP.animal);
end


function GetSeriesrange_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.explist = zeros(size(EXP.data.es.explist));
for i = 1:size(EXP.data.es.explist,1)
    [expsel,ok] = listdlg('ListString',strsplit(num2str(EXP.data.es.explist(i,1:find(EXP.data.es.explist(i,:)>0,1,'last')))), 'Name', ['selected exp for series' num2str(EXP.data.es.serieslist(i))], 'InitialValue', 1:find(EXP.data.es.explist(i,:)>0,1,'last'));
    if ok
        PlotVar1DMaps.explist(i,1:numel(expsel)) = EXP.data.es.explist(i,expsel);
    end
end
end

function GetDecSeriesrange_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.explist = zeros(size(EXP.data.es.explist));
for i = 1:size(EXP.data.es.explist,1)
    [expsel,ok] = listdlg('ListString',strsplit(num2str(EXP.data.es.explist(i,1:find(EXP.data.es.explist(i,:)>0,1,'last')))), 'Name', ['selected exp for series' num2str(EXP.data.es.serieslist(i))], 'InitialValue', 1:find(EXP.data.es.explist(i,:)>0,1,'last'));
    if ok
        PlotVarDecoder.explist(i,1:numel(expsel)) = EXP.data.es.explist(i,expsel);
    end
end
end

function GetPopSeriesrange_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.explist = zeros(size(EXP.data.es.explist));
for i = 1:size(EXP.data.es.explist,1)
    [expsel,ok] = listdlg('ListString',strsplit(num2str(EXP.data.es.explist(i,1:find(EXP.data.es.explist(i,:)>0,1,'last')))), 'Name', ['selected exp for series' num2str(EXP.data.es.serieslist(i))], 'InitialValue', 1:find(EXP.data.es.explist(i,:)>0,1,'last'));
    if ok
        PlotVarPop.explist(i,1:numel(expsel)) = EXP.data.es.explist(i,expsel);
    end
end
end

function GetBehSeriesrange_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
PlotVarBehavior.explist = zeros(size(EXP.data.es.explist));
for i = 1:size(EXP.data.es.explist,1)
    [expsel,ok] = listdlg('ListString',strsplit(num2str(EXP.data.es.explist(i,1:find(EXP.data.es.explist(i,:)>0,1,'last')))), 'Name', ['selected exp for series' num2str(EXP.data.es.serieslist(i))], 'InitialValue', 1:find(EXP.data.es.explist(i,:)>0,1,'last'));
    if ok
        PlotVarBehavior.explist(i,1:numel(expsel)) = EXP.data.es.explist(i,expsel);
    end
end
end


function SelectSeries_callback(hObject, eventdata, DIRS, dialog, EXP)
EXP.SelectSeries();
dialog.UpdateUIcontrol('Series','String', EXP.series);
end

function SelectExpt_callback(hObject, eventdata, DIRS, dialog1, dialog2, dialogbeh, dialogpop, EXP, PlotVarDecoder, PlotVar1DMaps, PlotVarBehavior, PlotVarPop)
%Loading callback function
EXP.SelectExpt();
dialog1.UpdateUIcontrol('Exp','String', num2str(EXP.exptList));
PlotVar1DMaps.samplerate = 60;
EXP.LoadVRData(PlotVar1DMaps.Shanknum, PlotVar1DMaps.Shanksuffix,PlotVar1DMaps.SpeedThreshold,PlotVar1DMaps.Nthetaphsbins,PlotVar1DMaps.SmthTimeWindow,PlotVar1DMaps.samplerate);

EXP.data.es.series = 0*EXP.data.es.series + EXP.iseries;
EXP.data.es.explist = zeros(numel(unique(EXP.data.es.series)),numel(unique(EXP.data.es.iexp)));

seriesnum = unique(EXP.data.es.series);
EXP.data.es.serieslist = seriesnum(:);
for i = 1:numel(seriesnum)
    iexpnum = unique(EXP.data.es.iexp(EXP.data.es.series == seriesnum(i)));
    EXP.data.es.explist(i,1:numel(iexpnum)) = iexpnum;
end

PlotVarDecoder.explist = EXP.data.es.explist;
PlotVar1DMaps.explist = EXP.data.es.explist;
PlotVarBehavior.explist = EXP.data.es.explist;
PlotVarPop.explist = EXP.data.es.explist;

assignin('base',[EXP.animal '_es_' num2str(EXP.iseries)],EXP.data.es);

if EXP.data.ephys || EXP.data.twophoton
    if ~EXP.data.twophoton
        if ~isdir([DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'])
            mkdir([DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed']);
        end
        savedcellinfo = [DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed' filesep EXP.animal '_' num2str(EXP.iseries) '_cellProperties.mat'];
        savedmaps1d = [DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed' filesep EXP.animal '_' num2str(EXP.iseries) '_maps1d.mat'];
        savedVS = [DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed' filesep EXP.animal '_' num2str(EXP.iseries) '_VS.mat'];
        saveddecsinglecell = [DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed' filesep EXP.animal '_' num2str(EXP.iseries) '_DecodingPopCell.mat'];
    else
        if ~isdir([DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'])
            mkdir([DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed']);
        end
        savedcellinfo = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep EXP.animal '_' num2str(EXP.iseries) '_cellProperties.mat'];
        savedmaps1d = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep EXP.animal '_' num2str(EXP.iseries) '_maps1d.mat'];
        savedmaps1dphs = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep EXP.animal '_' num2str(EXP.iseries) '_maps1d_phs.mat'];
        savedmaps2dphs = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep EXP.animal '_' num2str(EXP.iseries) '_maps2d_phs.mat'];
        savedmaps1d_2fold = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep EXP.animal '_' num2str(EXP.iseries) '_maps1d_2fold.mat'];
        savedmaps1dphs_2fold = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep EXP.animal '_' num2str(EXP.iseries) '_maps1dphs_2fold.mat'];
        savedmaps2dphs_2fold = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep EXP.animal '_' num2str(EXP.iseries) '_maps2d_phs_2fold.mat'];
        savedVS = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep EXP.animal '_' num2str(EXP.iseries) '_VS.mat'];
        saveddecsinglecell = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep EXP.animal '_' num2str(EXP.iseries) '_DecodingPopCell.mat'];
    end
    
    %number of random permutations for computing significance of single cell responses
    Nperm_cellprop = 0;%100; %takes a long time so use >2 only if useful
    %whether or not to overwrite the previously saved analysis
    FloadfromfileCellinfo = PlotVar1DMaps.FloadfromfileCellinfo;%true;%false;

    if exist(savedcellinfo) && exist(savedmaps1d) && exist(savedVS) && FloadfromfileCellinfo
        disp('loading cellInfo from saved file');
        S = load(savedcellinfo);
        EXP.CellInfo = S.CellInfo;
        S = load(savedVS);
        EXP.data.VStuning = S.VStuning;
        EXP.data.VSstim = S.VSstim;
        S = load(savedmaps1d);
        EXP.maps1d.trajPercent = S.maps1d;
        S = load(savedmaps1dphs);
        EXP.maps1d.LFPphase2 = S.maps1d;
        S = load(savedmaps2dphs);
        EXP.maps2d.trajPercent_LFPphase2 = S.maps1d;
    else
        EXP.maps1d = [];
        EXP.maps2d = [];
        EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs,EXP.data.es.ProbeIDs);
        %we call definecellprop here to get the correction factor for the LFPphase
        if isfield(EXP.data.es,'LFPphase')
            EXP.defineCellProp(Nperm_cellprop);
            nProbe = numel(unique(EXP.CellInfo.Probe));
            EXP.data.es.LFPphase2 = NaN(size(EXP.data.es.LFPphase,1),nProbe);
            for iprobe = 1:nProbe
                EXP.data.es.LFPphase2(:,iprobe) = mod(EXP.data.es.LFPphase - round(EXP.CellInfo.LFP2Spike_phscorrMUAnorm{1}),360);% - EXP.CellInfo.LFP2Spike_phscorrMUAnorm{iprobe},360);%- EXP.CellInfo.LFP2Spike_phscorrMUA{iprobe},360);
            end
        end
        EXP.CalculateStimTuning([], PlotVar1DMaps.Shanknum, PlotVar1DMaps.Shanksuffix);
        delay = 0;
        Xbinsize = 1;
        tic
        EXP.Calculate1Dmaps('trajPercent', PlotVar1DMaps.SmthTimeWindow, Xbinsize, PlotVar1DMaps.SmthSpatialWindow,delay, EXP.data.es.CircularMaze);
        toc
%         Spdbinsize = 0.1;
%         SmthSpdWindow = 0.1;
%         EXP.Calculate2Dmaps('trajPercent', 'smthBallSpd', PlotVar1DMaps.SmthTimeWindow, Xbinsize, Spdbinsize, PlotVar1DMaps.SmthSpatialWindow, SmthSpdWindow,delay, EXP.data.es.CircularMaze, false);
        if isfield(EXP.data.es,'LFPphase2')
            Phsbinsize = 20;
            SmthPhsWindow = 40;
            tic
            EXP.Calculate1Dmaps('LFPphase2', PlotVar1DMaps.SmthTimeWindow, Phsbinsize, SmthPhsWindow,delay, true);
            toc
            tic
            EXP.Calculate2Dmaps('trajPercent', 'LFPphase2', PlotVar1DMaps.SmthTimeWindow, Xbinsize, Phsbinsize, PlotVar1DMaps.SmthSpatialWindow, SmthPhsWindow,delay, EXP.data.es.CircularMaze, true);
            toc
        end
        EXP.defineCellProp(Nperm_cellprop);
    end
%     if ~(exist(savedcellinfo) && exist(savedmaps1d) && exist(savedVS))
%         VStuning = EXP.data.VStuning;
%         VSstim = EXP.data.VSstim;
%         save(savedVS,'VStuning','VSstim');
%         maps1d = EXP.maps1d.trajPercent;
%         save(savedmaps1d,'maps1d');
%         maps1d = EXP.maps1d.trajPercent_2fold;
%         save(savedmaps1d_2fold,'maps1d');        
%         maps2d = EXP.maps2d.trajPercent_LFPphase;
%         save(savedmaps2dphs,'maps2d');
%         maps2d = EXP.maps2d.trajPercent_LFPphase_2fold;
%         save(savedmaps2d_2fold,'maps2d');        
%         CellInfo = EXP.CellInfo;
%         save(savedcellinfo,'CellInfo');
%     end
    
    if exist(saveddecsinglecell)
        disp('loading decoding trajectories from saved file');
        S = load(saveddecsinglecell);
        EXP.CellInfo.trajdecCell_Xave = S.res.Xpred_ave;
        EXP.CellInfo.trajdecCell_Xmax = S.res.Xpred_max;
        EXP.CellInfo.trajdecAll_Xave = S.resAll.Xpred_ave;
        EXP.CellInfo.trajdecAll_Xmax = S.resAll.Xpred_max;
    else
        EXP.CellInfo.trajdecCell_Xave = cell(1,EXP.CellInfo.NumCells);
        EXP.CellInfo.trajdecCell_Xmax = cell(1,EXP.CellInfo.NumCells);
        EXP.CellInfo.trajdecAll_Xave = cell(1,EXP.CellInfo.NumCells);
        EXP.CellInfo.trajdecAll_Xmax = cell(1,EXP.CellInfo.NumCells);
    end
end


%****updating the dialogs parameters****%
cbase = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
gbase = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
rbase = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
obase = find(EXP.SubsetVal.outcome == 2);
if isempty(obase)
    obase = 1;
end
if ~isempty(EXP.CellInfo)
    dialog1.UpdateUIcontrol('CA1 - cell#','String', EXP.CellInfo.CellListString(EXP.CellInfo.Probe == 1), 'max', sum(EXP.CellInfo.Probe == 1));
    dialog1.UpdateUIcontrol('V1 - cell#','String', EXP.CellInfo.CellListString(EXP.CellInfo.Probe == 2), 'max', sum(EXP.CellInfo.Probe == 2));
end
dialog2.UpdateUIcontrol('X variable','String', PlotVar1DMaps.PlotVarList, 'max', 1, 'min', 0, 'Val', 1); 
PlotVar1DMaps.ChosenVarX = PlotVar1DMaps.PlotVarList(1);
dialog2.UpdateUIcontrol('Y variable','String', PlotVar1DMaps.PlotVarList, 'max', 1, 'min', 0, 'Val', 3);
PlotVar1DMaps.ChosenVarY = PlotVar1DMaps.PlotVarList(3);
dialog2.UpdateUIcontrol('Plots','String', PlotVar1DMaps.PlotObjList, 'max', numel(PlotVar1DMaps.PlotObjList), 'min', 0, 'Val', 1); 
PlotVar1DMaps.ChosenObj = PlotVar1DMaps.PlotObjList(1);
dialog2.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', cbase); 
PlotVar1DMaps.ChosenContrast = cbase;
dialog2.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', gbase); 
PlotVar1DMaps.ChosenGain = gbase;
dialog2.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', rbase);
PlotVar1DMaps.ChosenRoomlength = rbase;
dialog2.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', obase);
PlotVar1DMaps.ChosenOutcome = obase;

dialogbeh.UpdateUIcontrol('Plots','String', PlotVarBehavior.PlotObjList, 'max', numel(PlotVarBehavior.PlotObjList), 'min', 0, 'Val', 1); 
PlotVarBehavior.ChosenObj = PlotVarBehavior.PlotObjList(1);
dialogbeh.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', cbase); 
PlotVarBehavior.ChosenContrast = cbase;
dialogbeh.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', gbase); 
PlotVarBehavior.ChosenGain = gbase;
dialogbeh.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', rbase);
PlotVarBehavior.ChosenRoomlength = rbase;
dialogbeh.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', obase);
PlotVarBehavior.ChosenOutcome = obase;

dialogpop.UpdateUIcontrol('Plots','String', PlotVarPop.PlotObjList, 'max', numel(PlotVarPop.PlotObjList), 'min', 0, 'Val', 1); 
PlotVarPop.ChosenObj = PlotVarPop.PlotObjList(1);
dialogpop.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', cbase); 
PlotVarPop.ChosenContrast = cbase;
dialogpop.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', gbase); 
PlotVarPop.ChosenGain = gbase;
dialogpop.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', rbase);
PlotVarPop.ChosenRoomlength = rbase;
dialogpop.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', obase);
PlotVarPop.ChosenOutcome = obase;
end

function LoadSavedExpt_callback(hObject, eventdata, DIRS, dialog1, dialog2, dialogbeh, dialogpop, dialogdec, EXP, PlotVarDecoder, PlotVar1DMaps, PlotVarBehavior, PlotVarPop)
%Loading callback function
if isdir(['D:\DATA\batch' filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'])
    dirname = ['D:\DATA\batch' filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];%[DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];
elseif isdir([DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries)])
    dirname = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries)];
else
    error('no processed file yet');
end
filename = uigetfile(dirname);
S = load([dirname filesep filename]);
EXP.Copyobj(S.EXP);

dialog1.UpdateUIcontrol('Exp','String', num2str(EXP.exptList));
EXP.data.es.series = 0*EXP.data.es.series + EXP.iseries;
EXP.data.es.explist = zeros(numel(unique(EXP.data.es.series)),numel(unique(EXP.data.es.iexp)));

seriesnum = unique(EXP.data.es.series);
EXP.data.es.serieslist = seriesnum(:);
for i = 1:numel(seriesnum)
    iexpnum = unique(EXP.data.es.iexp(EXP.data.es.series == seriesnum(i)));
    EXP.data.es.explist(i,1:numel(iexpnum)) = iexpnum;
end

PlotVarDecoder.explist = EXP.data.es.explist;
PlotVar1DMaps.explist = EXP.data.es.explist;
PlotVarBehavior.explist = EXP.data.es.explist;
PlotVarPop.explist = EXP.data.es.explist;

% assignin('base',[EXP.animal '_es_' num2str(EXP.iseries)],EXP.data.es);


%****updating the dialogs parameters****%
cbase = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
gbase = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
rbase = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
obase = find(EXP.SubsetVal.outcome == 2);
if isempty(obase)
    obase = 1;
end
if ~isempty(EXP.CellInfo)
    dialog1.UpdateUIcontrol('CA1 - cell#','String', EXP.CellInfo.CellListString(EXP.CellInfo.Probe == 1), 'max', sum(EXP.CellInfo.Probe == 1));
    dialog1.UpdateUIcontrol('V1 - cell#','String', EXP.CellInfo.CellListString(EXP.CellInfo.Probe == 2), 'max', sum(EXP.CellInfo.Probe == 2));
end
dialog2.UpdateUIcontrol('X variable','String', PlotVar1DMaps.PlotVarList, 'max', 1, 'min', 0, 'Val', 1); 
PlotVar1DMaps.ChosenVarX = PlotVar1DMaps.PlotVarList(1);
dialog2.UpdateUIcontrol('Y variable','String', PlotVar1DMaps.PlotVarList, 'max', 1, 'min', 0, 'Val', 3);
PlotVar1DMaps.ChosenVarY = PlotVar1DMaps.PlotVarList(3);
dialog2.UpdateUIcontrol('Plots','String', PlotVar1DMaps.PlotObjList, 'max', numel(PlotVar1DMaps.PlotObjList), 'min', 0, 'Val', 1); 
PlotVar1DMaps.ChosenObj = PlotVar1DMaps.PlotObjList(1);
dialog2.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', cbase); 
PlotVar1DMaps.ChosenContrast = cbase;
dialog2.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', gbase); 
PlotVar1DMaps.ChosenGain = gbase;
dialog2.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', rbase);
PlotVar1DMaps.ChosenRoomlength = rbase;
dialog2.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', obase);
PlotVar1DMaps.ChosenOutcome = obase;

dialogbeh.UpdateUIcontrol('Plots','String', PlotVarBehavior.PlotObjList, 'max', numel(PlotVarBehavior.PlotObjList), 'min', 0, 'Val', 1); 
PlotVarBehavior.ChosenObj = PlotVarBehavior.PlotObjList(1);
dialogbeh.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', cbase); 
PlotVarBehavior.ChosenContrast = cbase;
dialogbeh.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', gbase); 
PlotVarBehavior.ChosenGain = gbase;
dialogbeh.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', rbase);
PlotVarBehavior.ChosenRoomlength = rbase;
dialogbeh.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', obase);
PlotVarBehavior.ChosenOutcome = obase;

dialogpop.UpdateUIcontrol('Plots','String', PlotVarPop.PlotObjList, 'max', numel(PlotVarPop.PlotObjList), 'min', 0, 'Val', 1); 
PlotVarPop.ChosenObj = PlotVarPop.PlotObjList(1);
dialogpop.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', cbase); 
PlotVarPop.ChosenContrast = cbase;
dialogpop.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', gbase); 
PlotVarPop.ChosenGain = gbase;
dialogpop.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', rbase);
PlotVarPop.ChosenRoomlength = rbase;
dialogpop.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', obase);
PlotVarPop.ChosenOutcome = obase;

cbase = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
gbase = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
rbase = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
obase = find(EXP.SubsetVal.outcome == 2);

dialogdec.UpdateUIcontrol('Plots','String', PlotVarDecoder.PlotObjList, 'max', numel(PlotVarDecoder.PlotObjList), 'min', 0, 'Val', 1); 
PlotVarDecoder.ChosenObj = PlotVarDecoder.PlotObjList(1);
dialogdec.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', cbase); 
PlotVarDecoder.ChosenContrast = cbase;
dialogdec.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', gbase); 
PlotVarDecoder.ChosenGain = gbase;
dialogdec.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', rbase);
PlotVarDecoder.ChosenRoomlength = rbase;
dialogdec.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', obase);
PlotVarDecoder.ChosenOutcome = obase;
end

function SelectMultipleSeries_callback(hObject, eventdata, DIRS, dialog1, dialog2, dialogbeh, dialogpop, EXP, PlotVar1DMaps, PlotVarBehavior, PlotVarPop, DecodingDialog, PlotVarDecoder,Layout)
%To concatenate multiple exp which have already been processed and saved
prompt = {'# of animals','file suffix'};
dlg_title = 'how many animals?';
num_lines = 1;
def = {'1','win150'};
nani = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(nani)
    nanimals = str2num(nani{1});
    suffix = nani{2};
else
    error('no valid animal number');
end
ianimal = cell(nanimals,1);
iseries = cell(nanimals,1);
series = cell(nanimals,1);

for k = 1:nanimals
    EXP.SelectAnimal(DIRS,{'BALL','JF'});
    EXP.SelectSeries('multiple');
    ianimal{k} = EXP.animal;
    iseries{k} = EXP.iseries;
    series{k} = EXP.series;
end

n = 0;
MaxLickdistri = cell(2,1);
MaxMeanPostLick = cell(2,1);
MaxPostLick = cell(2,1);
for iprobe = 1:2
    MaxLickdistri{iprobe} = zeros(numel(ianimal),4);
    MaxMeanPostLick{iprobe} = zeros(numel(ianimal),4,60);
    MaxPostLick{iprobe} = zeros(numel(ianimal),4,60);
end

for k = 1:numel(ianimal)
    for i = 1:numel(iseries{k})
        n = n+1;
        disp(num2str(iseries{k}(i)));
        EXP.animal = ianimal{k};
        EXP.series = series{k}(i);
        EXP.iseries = iseries{k}(i);
        
        if k == 1 && i == 1
            S = load([DIRS.multichanspikes filesep EXP.animal filesep num2str(iseries{k}(i)) filesep 'EXP_' suffix '.mat']);
            if size(S.EXP.Bayes.Posterior,1) < 2
                for phsidx = 1:size(S.EXP.Bayes.Posterior,2)
                    S.EXP.Bayes.Posterior{2,phsidx} = NaN(size(S.EXP.Bayes.Posterior{1,phsidx}));
                    S.EXP.Bayes.PosError{2,phsidx} = NaN(size(S.EXP.Bayes.PosError{1,phsidx}));
                    S.EXP.Bayes.DistError{2,phsidx} = NaN(size(S.EXP.Bayes.DistError{1,phsidx}));
                    S.EXP.Bayes.prediction{2,phsidx} = NaN(size(S.EXP.Bayes.prediction{1,phsidx}));
                end
                S.EXP.Bayes.Posterior0{2} = NaN(size(S.EXP.Bayes.Posterior0{1}));
                S.EXP.Bayes.PosError0{2} = NaN(size(S.EXP.Bayes.PosError0{1}));
                S.EXP.Bayes.DistError0{2} = NaN(size(S.EXP.Bayes.DistError0{1}));
                S.EXP.Bayes.LFPphase{2} = NaN(size(S.EXP.Bayes.LFPphase{1}));
            end
            S.EXP.data.es.series = 0*S.EXP.data.es.series + n;
            EXP.Copyobj(S.EXP);
        else
            S = load([DIRS.multichanspikes filesep EXP.animal filesep num2str(iseries{k}(i)) filesep 'EXP_' suffix '.mat']);
            if size(S.EXP.Bayes.Posterior,1) < 2
                for phsidx = 1:size(S.EXP.Bayes.Posterior,2)
                    S.EXP.Bayes.Posterior{2,phsidx} = NaN(size(S.EXP.Bayes.Posterior{1,phsidx}));
                    S.EXP.Bayes.PosError{2,phsidx} = NaN(size(S.EXP.Bayes.PosError{1,phsidx}));
                    S.EXP.Bayes.DistError{2,phsidx} = NaN(size(S.EXP.Bayes.DistError{1,phsidx}));
                    S.EXP.Bayes.prediction{2,phsidx} = NaN(size(S.EXP.Bayes.prediction{1,phsidx}));
                end
                S.EXP.Bayes.Posterior0{2} = NaN(size(S.EXP.Bayes.Posterior0{1}));
                S.EXP.Bayes.PosError0{2} = NaN(size(S.EXP.Bayes.PosError0{1}));
                S.EXP.Bayes.DistError0{2} = NaN(size(S.EXP.Bayes.DistError0{1}));
                S.EXP.Bayes.LFPphase{2} = NaN(size(S.EXP.Bayes.LFPphase{1}));
            end
            S.EXP.data.es.series = 0*S.EXP.data.es.series + n;
            EXP.Appenobj(S.EXP);
        end
    end
end
EXP.data.es = getESDataSubset(EXP.data.es, 'smthBallSpd', 5, []);

%to simplify display selection
% EXP.data.es.contrast(EXP.data.es.contrast < 0.5 & EXP.data.es.contrast > 0) = 0.2;
% EXP.data.es.contrast(EXP.data.es.contrast == 0.5 | EXP.data.es.contrast == 0.7) = 0.5;
% EXP.data.es.contrast(EXP.data.es.contrast > 0.7) = 0.8;

EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs, EXP.data.es.ProbeIDs);
EXP.CalculateSubsets();


%****Updating dialogs parameters****%

EXP.Bayes.TrainContrast = find(EXP.SubsetVal.contrast == 0.5);
EXP.Bayes.TrainGain = find(EXP.SubsetVal.gain == 0.5);
EXP.Bayes.TrainRoomlength = find(EXP.SubsetVal.roomlength == 1);
EXP.Bayes.TrainOutcome = find(EXP.SubsetVal.outcome == 2);
% EXP.defineCellProp;%defineCellProp only works for single series for now

if ~isempty(EXP.CellInfo)
    dialog1.UpdateUIcontrol('CA1 - cell#','String', EXP.CellInfo.CellListString, 'max', EXP.CellInfo.NumCells);
end
dialog2.UpdateUIcontrol('X variable','String', PlotVar1DMaps.PlotVarList, 'max', 1, 'min', 0, 'Val', 1); 
PlotVar1DMaps.ChosenVarX = PlotVar1DMaps.PlotVarList(1);
dialog2.UpdateUIcontrol('Y variable','String', PlotVar1DMaps.PlotVarList, 'max', 1, 'min', 0, 'Val', 3);
PlotVar1DMaps.ChosenVarY = PlotVar1DMaps.PlotVarList(3);
dialog2.UpdateUIcontrol('Plots','String', PlotVar1DMaps.PlotObjList, 'max', numel(PlotVar1DMaps.PlotObjList), 'min', 0, 'Val', 1); 
PlotVar1DMaps.ChosenObj = PlotVar1DMaps.PlotObjList(1);
dialog2.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', 1); 
PlotVar1DMaps.ChosenContrast = 1:numel(EXP.SubsetVal.contrast);
dialog2.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', 1); 
PlotVar1DMaps.ChosenGain = 1:numel(EXP.SubsetVal.gain);
dialog2.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', 1);
PlotVar1DMaps.ChosenRoomlength = 1:numel(EXP.SubsetVal.roomlength);
dialog2.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', 1);
PlotVar1DMaps.ChosenOutcome = 1:numel(EXP.SubsetVal.outcome);

dialogbeh.UpdateUIcontrol('Plots','String', PlotVarBehavior.PlotObjList, 'max', numel(PlotVarBehavior.PlotObjList), 'min', 0, 'Val', 1); 
PlotVarBehavior.ChosenObj = PlotVarBehavior.PlotObjList(1);
dialogbeh.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', 1); 
PlotVarBehavior.ChosenContrast = 1:numel(EXP.SubsetVal.contrast);
dialogbeh.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', 1); 
PlotVarBehavior.ChosenGain = 1:numel(EXP.SubsetVal.gain);
dialogbeh.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', 1);
PlotVarBehavior.ChosenRoomlength = 1:numel(EXP.SubsetVal.roomlength);
dialogbeh.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', 1);
PlotVarBehavior.ChosenOutcome = 1:numel(EXP.SubsetVal.outcome);

dialogpop.UpdateUIcontrol('Plots','String', PlotVarPop.PlotObjList, 'max', numel(PlotVarPop.PlotObjList), 'min', 0, 'Val', 1); 
PlotVarPop.ChosenObj = PlotVarPop.PlotObjList(1);
dialogpop.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', 1); 
PlotVarPop.ChosenContrast = 1:numel(EXP.SubsetVal.contrast);
dialogpop.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', 1); 
PlotVarPop.ChosenGain = 1:numel(EXP.SubsetVal.gain);
dialogpop.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', 1);
PlotVarPop.ChosenRoomlength = 1:numel(EXP.SubsetVal.roomlength);
dialogpop.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', 1);
PlotVarPop.ChosenOutcome = 1:numel(EXP.SubsetVal.outcome);

DecodingDialog.UpdateUIcontrol('Plots','String', PlotVarDecoder.PlotObjList, 'max', numel(PlotVarDecoder.PlotObjList), 'min', 0, 'Val', 1); 
PlotVarDecoder.ChosenObj = PlotVarDecoder.PlotObjList(1);
DecodingDialog.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', 1:numel(EXP.SubsetVal.contrast)); 
PlotVarDecoder.ChosenContrast = (1:numel(EXP.SubsetVal.contrast));
DecodingDialog.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', 1:numel(EXP.SubsetVal.gain)); 
PlotVarDecoder.ChosenGain = (1:numel(EXP.SubsetVal.gain));
DecodingDialog.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', 1);
PlotVarDecoder.ChosenRoomlength = (1:EXP.SubsetVal.roomlength);
DecodingDialog.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', find(EXP.SubsetVal.outcome == 2));
PlotVarDecoder.ChosenOutcome = find(EXP.SubsetVal.outcome == 2);
end

function Compute1Dmap_callback(hObject, eventdata, EXP, PlotVar1DMaps)
%Computes cross-validated 1D maps
delay = 0;
EXP.Calculate1Dmaps('trajPercent',PlotVar1DMaps.SmthTimeWindow, PlotVar1DMaps.SmthSpatialWindow,delay);
end

function BayesDecoder_callback(hObject, eventdata, EXP, PlotVarDecoder, DecodingDialog, Layout)
%Runs the Bayesian decoder
EXP.RunBayesDecoder;%EXP.RunBayesDecoder2D;%('trajPercent','spikeTrain');%('distTrav','spikeTrain');%


%Updating the parameters of the dialogs%
cbase = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
gbase = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
rbase = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
obase = find(EXP.SubsetVal.outcome == 2);

DecodingDialog.UpdateUIcontrol('Plots','String', PlotVarDecoder.PlotObjList, 'max', numel(PlotVarDecoder.PlotObjList), 'min', 0, 'Val', 1); 
PlotVarDecoder.ChosenObj = PlotVarDecoder.PlotObjList(1);
DecodingDialog.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', cbase); 
PlotVarDecoder.ChosenContrast = cbase;
DecodingDialog.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', gbase); 
PlotVarDecoder.ChosenGain = gbase;
DecodingDialog.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', rbase);
PlotVarDecoder.ChosenRoomlength = rbase;
DecodingDialog.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', obase);
PlotVarDecoder.ChosenOutcome = obase;

PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function BatchBayesDecoder_callback(hObject, eventdata, DIRS, EXP, PlotVarDecoder, DecodingDialog, Layout)
%Batch Bayesian decoding over multiple experiments. Edit
%BatchBayesDecoder.m before using
[EXP, PlotVarDecoder, DecodingDialog] = BatchBayesDecoder(EXP, DIRS, PlotVarDecoder, DecodingDialog);
end

function BatchBayesDecoderParameters_callback(hObject, eventdata, DIRS, EXP, PlotVarDecoder, DecodingDialog, Layout)
%Batch the Bayesian decoding over values of multiple parameters 
%edit BatchBayesDecoderParameters.m before using
EXP = BatchBayesDecoderParameters(EXP, DIRS);
end

function BatchBayesDecoderSingleCells_callback(hObject, eventdata, DIRS, EXP, PlotVarDecoder, DecodingDialog, Layout)
%Batch the Bayesian decoding omitting one cell at every iteration 
%edit BatchBayesDecoderSingleCells.m before using
% EXP = BatchBayesSingleCells(EXP, DIRS);
end

function GetChosenProbe_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.ChosenProbe = get(hObject,'Val');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Getmaxtolerance_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = str2num(get(hObject,'String'));
PlotVarDecoder.maxtolerance = val;
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function FdecXdistribution_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FdecXdistribution = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end


function Fgoodtimebins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.Fgoodtimebins = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Fspatialsmooth_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.Fspatialsmooth = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end


function FthetaPost_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FthetaPost = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Fthetabins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.Fthetabins = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Fspdbins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.Fspdbins = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function GetdecthetaChannel_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = str2num(get(hObject,'String'));
PlotVarDecoder.thetaChannel = val;
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function GetthetaDecphase_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = str2num(get(hObject,'String'));
PlotVarDecoder.thetaDecphase = val;
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Getspdbin_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = str2num(get(hObject,'String'));
PlotVarDecoder.spdbin = val;
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Getdecthetaphase_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = str2num(get(hObject,'String'));
PlotVarDecoder.thetaphase = val;
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Getdecnthetaphsbins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = str2num(get(hObject,'String'));
PlotVarDecoder.nthetaphsbins = val;
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function GetdecnthetaXbins_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = str2num(get(hObject,'String'));
PlotVarDecoder.nthetaXbins = val;
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function SelectDecX_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = get(hObject,'Value');
str = get(hObject,'String');
PlotVarDecoder.ChosenVarX = str(val);
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function SelectDecY_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = get(hObject,'Value');
str = get(hObject,'String');
PlotVarDecoder.ChosenVarY = str(val);
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function SelectDecPlotObj_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = get(hObject,'Value');
str = get(hObject,'String');
PlotVarDecoder.ChosenObj = str(val);
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function SelectDecContrast_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.ChosenContrast = get(hObject,'Value');
if PlotVarDecoder.FpoolContrast
    PlotVarDecoder.ChosenContrast = PlotVarDecoder.ChosenContrast(:);
else
    PlotVarDecoder.ChosenContrast = PlotVarDecoder.ChosenContrast(:)';
end
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function SelectDecGain_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.ChosenGain = get(hObject,'Value');
if PlotVarDecoder.FpoolGain
    PlotVarDecoder.ChosenGain = PlotVarDecoder.ChosenGain(:);
else
    PlotVarDecoder.ChosenGain = PlotVarDecoder.ChosenGain(:)';
end
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function SelectDecRoomlength_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.ChosenRoomlength = get(hObject,'Value');
if PlotVarDecoder.FpoolRoomlength
    PlotVarDecoder.ChosenRoomlength = PlotVarDecoder.ChosenRoomlength(:);
else
    PlotVarDecoder.ChosenRoomlength = PlotVarDecoder.ChosenRoomlength(:)';
end
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function SelectDecOutcome_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.ChosenOutcome = get(hObject,'Value');
if PlotVarDecoder.FpoolOutcome
    PlotVarDecoder.ChosenOutcome = PlotVarDecoder.ChosenOutcome(:);
else
    PlotVarDecoder.ChosenOutcome = PlotVarDecoder.ChosenOutcome(:)';
end
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end


function FPoolDecContrast_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FpoolContrast = get(hObject,'Value');
if PlotVarDecoder.FpoolContrast
    PlotVarDecoder.ChosenContrast = PlotVarDecoder.ChosenContrast(:);
else
    PlotVarDecoder.ChosenContrast = PlotVarDecoder.ChosenContrast(:)';
end
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function FPoolDecGain_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FpoolGain = get(hObject,'Value');
if PlotVarDecoder.FpoolGain
    PlotVarDecoder.ChosenGain = PlotVarDecoder.ChosenGain(:);
else
    PlotVarDecoder.ChosenGain = PlotVarDecoder.ChosenGain(:)';
end
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function FPoolDecRoomlength_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FpoolRoomlength = get(hObject,'Value');
if PlotVarDecoder.FpoolRoomlength
    PlotVarDecoder.ChosenRoomlength = PlotVarDecoder.ChosenRoomlength(:);
else
    PlotVarDecoder.ChosenRoomlength = PlotVarDecoder.ChosenRoomlength(:)';
end
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function FPoolDecOutcome_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FpoolOutcome = get(hObject,'Value');
if PlotVarDecoder.FpoolOutcome
    PlotVarDecoder.ChosenOutcome = PlotVarDecoder.ChosenOutcome(:);
else
    PlotVarDecoder.ChosenOutcome = PlotVarDecoder.ChosenOutcome(:)';
end
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Fdisplaylog_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.Fdisplaylog = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function GetClim_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = str2num(get(hObject,'String'));
PlotVarDecoder.Clim = val;
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Fnormalize_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.Fnormalize = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function GetPalettename_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
val = get(hObject,'Value');
str = (get(hObject,'String'));
PlotVarDecoder.Palettename = str{val};
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function Foverlap_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.Foverlap = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function FdisplayPredMax_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FdisplayPredMax = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function FdisplayPredAve_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FdisplayPredAve = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function FdisplayMat_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FdisplayMat = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function FdisplayLFP_callback(hObject, eventdata, EXP, PlotVarDecoder, Layout)
PlotVarDecoder.FdisplayLFP = get(hObject,'Value');
PlotVarDecoder = UpdateplotDecodSpks(PlotVarDecoder,EXP,Layout);
end

function SelectPopPlotObj_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
val = get(hObject,'Value');
str = get(hObject,'String');
PlotVarPop.ChosenObj = str(val);
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function GetPopSSIrange_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
val = str2num(get(hObject,'String'));
PlotVarPop.SSIrange = val;
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function GetPopraterange_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
val = str2num(get(hObject,'String'));
PlotVarPop.Raterange = val;
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function GetPopzth_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
val = str2num(get(hObject,'String'));
PlotVarPop.zth = val;
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function GetPopProbe_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.Probes = get(hObject,'Val');
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function SelectPopContrast_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.ChosenContrast = get(hObject,'Value');
if PlotVarPop.FpoolContrast
    PlotVarPop.ChosenContrast = PlotVarPop.ChosenContrast(:);
else
    PlotVarPop.ChosenContrast = PlotVarPop.ChosenContrast(:)';
end
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function FPoolPopContrast_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.FpoolContrast = get(hObject,'Value');
if PlotVarPop.FpoolContrast
    PlotVarPop.ChosenContrast = PlotVarPop.ChosenContrast(:);
else
    PlotVarPop.ChosenContrast = PlotVarPop.ChosenContrast(:)';
end
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function SelectPopGain_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.ChosenGain = get(hObject,'Value');
if PlotVarPop.FpoolGain
    PlotVarPop.ChosenGain = PlotVarPop.ChosenGain(:);
else
    PlotVarPop.ChosenGain = PlotVarPop.ChosenGain(:)';
end
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function FPoolPopGain_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.FpoolGain = get(hObject,'Value');
if PlotVarPop.FpoolGain
    PlotVarPop.ChosenGain = PlotVarPop.ChosenGain(:);
else
    PlotVarPop.ChosenGain = PlotVarPop.ChosenGain(:)';
end
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function SelectPopRoomlength_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.ChosenRoomlength = get(hObject,'Value');
if PlotVarPop.FpoolRoomlength
    PlotVarPop.ChosenRoomlength = PlotVarPop.ChosenRoomlength(:);
else
    PlotVarPop.ChosenRoomlength = PlotVarPop.ChosenRoomlength(:)';
end
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function FPoolPopRoomlength_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.FpoolRoomlength = get(hObject,'Value');
if PlotVarPop.FpoolRoomlength
    PlotVarPop.ChosenRoomlength = PlotVarPop.ChosenRoomlength(:);
else
    PlotVarPop.ChosenRoomlength = PlotVarPop.ChosenRoomlength(:)';
end
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function SelectPopOutcome_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.ChosenOutcome = get(hObject,'Value');
if PlotVarPop.FpoolOutcome
    PlotVarPop.ChosenOutcome = PlotVarPop.ChosenOutcome(:);
else
    PlotVarPop.ChosenOutcome = PlotVarPop.ChosenOutcome(:)';
end
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function FPoolPopOutcome_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.FpoolOutcome = get(hObject,'Value');
if PlotVarPop.FpoolOutcome
    PlotVarPop.ChosenOutcome = PlotVarPop.ChosenOutcome(:);
else
    PlotVarPop.ChosenOutcome = PlotVarPop.ChosenOutcome(:)';
end
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function PopFnormalize_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
PlotVarPop.FNormalize = get(hObject,'Value');
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end

function GetPopClim_callback(hObject, eventdata, EXP, PlotVarPop, Layout)
val = str2num(get(hObject,'String'));
PlotVarPop.Clim = val;
PlotVarPop = UpdateplotPopSpks(PlotVarPop,EXP,Layout);
end


function SelectVarX_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = get(hObject,'Value');
str = get(hObject,'String');
PlotVar1DMaps.ChosenVarX = str(val);
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function SelectVarY_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = get(hObject,'Value');
str = get(hObject,'String');
PlotVar1DMaps.ChosenVarY = str(val);
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function SelectPlotObj_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = get(hObject,'Value');
str = get(hObject,'String');
PlotVar1DMaps.ChosenObj = str(val);
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function Getspeedth_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = str2num(get(hObject,'String'));
PlotVar1DMaps.speed_th = val;
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function GetdelayT_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = str2num(get(hObject,'String'));
PlotVar1DMaps.delayT = val;
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function GetXbinsize_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = str2num(get(hObject,'String'));
PlotVar1DMaps.Xbinsize = val;
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function GetthetaChannel_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
val = str2num(get(hObject,'String'));
PlotVar1DMaps.thetaChannel = val;
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end


function Fdispmat_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.Fdispmat = get(hObject,'Value');
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function FCombinedCond_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.FXcond = get(hObject,'Value');
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function SelectContrast_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.ChosenContrast = get(hObject,'Value');
if PlotVar1DMaps.FpoolContrast
    PlotVar1DMaps.ChosenContrast = PlotVar1DMaps.ChosenContrast(:);
else
    PlotVar1DMaps.ChosenContrast = PlotVar1DMaps.ChosenContrast(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function FPoolContrast_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.FpoolContrast = get(hObject,'Value');
if PlotVar1DMaps.FpoolContrast
    PlotVar1DMaps.ChosenContrast = PlotVar1DMaps.ChosenContrast(:);
else
    PlotVar1DMaps.ChosenContrast = PlotVar1DMaps.ChosenContrast(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function SelectGain_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.ChosenGain = get(hObject,'Value');
if PlotVar1DMaps.FpoolGain
    PlotVar1DMaps.ChosenGain = PlotVar1DMaps.ChosenGain(:);
else
    PlotVar1DMaps.ChosenGain = PlotVar1DMaps.ChosenGain(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function FPoolGain_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.FpoolGain = get(hObject,'Value');
if PlotVar1DMaps.FpoolGain
    PlotVar1DMaps.ChosenGain = PlotVar1DMaps.ChosenGain(:);
else
    PlotVar1DMaps.ChosenGain = PlotVar1DMaps.ChosenGain(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function SelectRoomlength_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.ChosenRoomlength = get(hObject,'Value');
if PlotVar1DMaps.FpoolRoomlength
    PlotVar1DMaps.ChosenRoomlength = PlotVar1DMaps.ChosenRoomlength(:);
else
    PlotVar1DMaps.ChosenRoomlength = PlotVar1DMaps.ChosenRoomlength(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function FPoolRoomlength_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.FpoolRoomlength = get(hObject,'Value');
if PlotVar1DMaps.FpoolRoomlength
    PlotVar1DMaps.ChosenRoomlength = PlotVar1DMaps.ChosenRoomlength(:);
else
    PlotVar1DMaps.ChosenRoomlength = PlotVar1DMaps.ChosenRoomlength(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function SelectOutcome_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.ChosenOutcome = get(hObject,'Value');
if PlotVar1DMaps.FpoolOutcome
    PlotVar1DMaps.ChosenOutcome = PlotVar1DMaps.ChosenOutcome(:);
else
    PlotVar1DMaps.ChosenOutcome = PlotVar1DMaps.ChosenOutcome(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end

function FPoolOutcome_callback(hObject, eventdata, EXP, PlotVar1DMaps, Layout)
PlotVar1DMaps.FpoolOutcome = get(hObject,'Value');
if PlotVar1DMaps.FpoolOutcome
    PlotVar1DMaps.ChosenOutcome = PlotVar1DMaps.ChosenOutcome(:);
else
    PlotVar1DMaps.ChosenOutcome = PlotVar1DMaps.ChosenOutcome(:)';
end
PlotVar1DMaps = UpdateplotBehavSpks(PlotVar1DMaps,EXP,Layout);
end



function SelectBehPlotObj_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
val = get(hObject,'Value');
str = get(hObject,'String');
PlotVarBehavior.ChosenObj = str(val);
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end

function GetBehspeedth_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
val = str2num(get(hObject,'String'));
PlotVarBehavior.speed_th = val;
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end

function SelectBehContrast_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
PlotVarBehavior.ChosenContrast = get(hObject,'Value');
if PlotVarBehavior.FpoolContrast
    PlotVarBehavior.ChosenContrast = PlotVarBehavior.ChosenContrast(:);
else
    PlotVarBehavior.ChosenContrast = PlotVarBehavior.ChosenContrast(:)';
end
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end

function SelectBehGain_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
PlotVarBehavior.ChosenGain = get(hObject,'Value');
if PlotVarBehavior.FpoolGain
    PlotVarBehavior.ChosenGain = PlotVarBehavior.ChosenGain(:);
else
    PlotVarBehavior.ChosenGain = PlotVarBehavior.ChosenGain(:)';
end
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end

function SelectBehRoomlength_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
PlotVarBehavior.ChosenRoomlength = get(hObject,'Value');
if PlotVarBehavior.FpoolRoomlength
    PlotVarBehavior.ChosenRoomlength = PlotVarBehavior.ChosenRoomlength(:);
else
    PlotVarBehavior.ChosenRoomlength = PlotVarBehavior.ChosenRoomlength(:)';
end
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end

function SelectBehOutcome_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
PlotVarBehavior.ChosenOutcome = get(hObject,'Value');
if PlotVarBehavior.FpoolOutcome
    PlotVarBehavior.ChosenOutcome = PlotVarBehavior.ChosenOutcome(:);
else
    PlotVarBehavior.ChosenOutcome = PlotVarBehavior.ChosenOutcome(:)';
end
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end

function FPoolBehContrast_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
PlotVarBehavior.FpoolContrast = get(hObject,'Value');
if PlotVarBehavior.FpoolContrast
    PlotVarBehavior.ChosenContrast = PlotVarBehavior.ChosenContrast(:);
else
    PlotVarBehavior.ChosenContrast = PlotVarBehavior.ChosenContrast(:)';
end
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end

function FPoolBehGain_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
PlotVarBehavior.FpoolGain = get(hObject,'Value');
if PlotVarBehavior.FpoolGain
    PlotVarBehavior.ChosenGain = PlotVarBehavior.ChosenGain(:);
else
    PlotVarBehavior.ChosenGain = PlotVarBehavior.ChosenGain(:)';
end
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end

function FPoolBehRoomlength_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
PlotVarBehavior.FpoolRoomlength = get(hObject,'Value');
if PlotVarBehavior.FpoolRoomlength
    PlotVarBehavior.ChosenRoomlength = PlotVarBehavior.ChosenRoomlength(:);
else
    PlotVarBehavior.ChosenRoomlength = PlotVarBehavior.ChosenRoomlength(:)';
end
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end

function FPoolBehOutcome_callback(hObject, eventdata, EXP, PlotVarBehavior, Layout)
PlotVarBehavior.FpoolOutcome = get(hObject,'Value');
if PlotVarBehavior.FpoolOutcome
    PlotVarBehavior.ChosenOutcome = PlotVarBehavior.ChosenOutcome(:);
else
    PlotVarBehavior.ChosenOutcome = PlotVarBehavior.ChosenOutcome(:)';
end
PlotVarBehavior = UpdateplotBehavior(PlotVarBehavior,EXP,Layout);
end
