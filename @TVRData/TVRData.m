classdef TVRData < TExpData
    %
    properties
        Subset
        SubsetVal
        SpeedThresh
        SmthTimeWindow
        maps1d
        maps2d
        Bayes
    end
    
    methods
        function obj = TVRData()
            obj@TExpData();
            
            obj.SpeedThresh = 5;
            obj.SmthTimeWindow = 150;
            obj.Subset.all = [];
            obj.Subset.gain = [];
            obj.Subset.contrast = [];
            obj.Subset.roomlength = [];
            obj.Subset.outcome = [];
            obj.Subset.Xall = [];
            
            obj.SubsetVal.gain = [];
            obj.SubsetVal.contrast = [];
            obj.SubsetVal.roomlength = [];
            obj.SubsetVal.outcome = [];
            
            obj.maps1d = [];
            obj.maps2d = [];
            obj.Bayes = [];
        end        
        
        obj = LoadVRData(obj, shank_list, suffix, speed_th, nthetaphsbins, SmthTimeWindow, samplerate);
        obj = CalculateSubsets(obj); 
        idx = getSubsets(obj, contrast, gain, roomlength, outcome, speed_th, FnoBlanks, FnoAfterBlanks);
        obj = defineCellProp(obj,Nperm);
        obj = Calculate1Dmaps(obj, varXname, Tsmthwin, Xbinsize, Xsmthwin, delayT, Fcircular);
        obj = Calculate2Dmaps(obj, varXname, varYname, Tsmthwin, Xbinsize, Ybinsize, Xsmthwin, Ysmthwin, delayT, FcircularX, FcircularY);
        obj = RunBayesDecoder(obj, varname, predictorname, varargin);
        obj = RunBayesDecoder2D(obj, varname, predictorname, varargin);
        obj = CalculateStimTuning(obj, exptlist, shank_list, suffix);
        obj = SimulPlaceFields(obj,delay,modeltype,VistoDistfactor);
    end
end