function EXP = ProcessVRData(animalname, series, exp)
%Process VR data using the basic functions run by VisNavGUI
%(The point of this file is that it can be run independently from the GUI) 
%To batch multiple files and save the processed data, 
%run EXP = BatchVRdata in the command window.

SetDirs;

%Parameters
Shanknum = [0:7 0];
Shanksuffix = {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
SpeedThreshold = 1;%cm/sec
Nthetaphsbins = 0;
SmthTimeWindow = 150;%ms 
SmthSpatialWindow = 4;%percent of the total maze length
delay = 0;
Nperm_cellprop = 2;

%create the main data object
EXP = TVRData;

if nargin <1
    EXP.animal = animalname;
    EXP.series = num2str(series);
    EXP.iseries = series;
    EXP.exptList = exp;
else
    EXP.SelectAnimal(DIRS,{'BALL','JF','MIK'});
    EXP.SelectSeries();
    EXP.SelectExpt();
end

%Load data
EXP.LoadVRData(Shanknum, Shanksuffix, SpeedThreshold, Nthetaphsbins);

%Tuning for VS
EXP.CalculateStimTuning([], Shanknum, Shanksuffix);

%Compute cross-validated 1D response profiles
EXP.Calculate1Dmaps('trajPercent', SmthTimeWindow, SmthSpatialWindow,delay);
EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs,EXP.data.es.ProbeIDs);
EXP.defineCellProp(Nperm_cellprop);

%Run the Bayesian decoder
%Call directly for EXP.RunBayesDecoder() if you want the dialogs to be displayed
traincont = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
traingain = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
trainroomlength = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
trainoutcome = find(EXP.SubsetVal.outcome == 2);
type = 'mean';
Flookuptable = false;
nruns = 1;
goodidx_th = 30;
Tsmth_win = 150;%20;%
Xsmth_win = 4;
numbins = 100;
thetachannel = 34;
nthetabins = 1;%1;%6;%
nspdbins = 5;%1;%
smth_spd = Tsmth_win;
speed_th = 0;
kfold = 20;
FoptiSmooth =0;
FGoodcluster = 1;
FUnsortedcluster = 0;
FMUAcluster = 0;
maxRate = inf;
zth = -inf;
SSImin = -inf;
latcorrection = 0;
alpha = 0;
delta = 0;

EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
    'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'smth_spd', smth_spd, 'latcorrection', latcorrection,...
    'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
    'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin);

%or call directly for this function if you want to select parameters from
%the dialogs
% EXP.RunBayesDecoder()

maxtol = 0.1;%range of tolerance from the max to compute the circular averages
EXP = BayesDecoderAverages(EXP,maxtol);
end