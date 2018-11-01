classdef TDecoder
    %     And abstract class that can decode neural data
    %       Methods include: train, predict (decode), ...
    %       Current decoders include: linear, treebagger and bayes
    % Aman Saleem
    % Oct 2013
    
    properties
        type;       % linear / treebagger / bayes
        variable;   % the description of the variable decoder
        variable_range; % min and max values of the variable used while training 
        Tsmth_win = 150;
        Xsmth_win = 1;
        Ysmth_win = 1;
        model;
        performance; % cross-validated (if kfold>1) performance at training
        meanPerformance; % mean performance
        sampleRate = 60;
    end
    
    methods (Abstract)
        [obj, prediction] = trainDecoder(obj, X, Y, Tsmth_win);
        [X] = predictVariable(obj, Y, Tsmth_win);
    end
    methods
        function obj = TDecoder(varargin)
            
            pnames = {'type' 'variable' 'variable_range'...
                'Tsmth_win' 'Xsmth_win' 'model' 'sampleRate'};
            dflts  = {'linear' 'P' []...
                150 1 [] 60};
            [obj.type, obj.variable, obj.variable_range, ...
                obj.Tsmth_win, obj.Xsmth_win, obj.model, obj.sampleRate] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
    end
end