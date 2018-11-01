classdef bayes2dDecoder < Decoder
    %     Create a Decoder class object of type '2D Bayes'. The decoder can be trained using the
    %     function 'trainDecoder' and used to decode with the function
    %     'predictVariable'
    %     Usage:
    % obj = bayes2dDecoder;
    % optional arguments
    %     'variable': variable being coded/decoded
    %     'variable_range': min and max values of variable during training
    %     'smth_win': smth_win used for training
    %     'kfold': number of cross-validation windows, 1 fors no
    %     cross-validation (dfl: 5)
    %     'CVO': cross-validation structure
    %     'sampleRate': the sampling rate of the data (dfl: 60Hz).
    %     'numBins':  number of bins by which the variable is discretized
    % [obj, prediction, X, Posterior, nPosterior] = obj.trainDecoder( X, Y, smth_win);
    % inputs: obj: object
    %         X: variable encoded/decoded
    %         Y: firing rate / signal
    %         smth_win(optional): the smoothing window to use on the signal (dfl: 150ms)
    % Ouputs: obj
    %         cross-validated prediction
    %         X (if there are NaN's in the data, these time points are taken out
    % [X] = obj.predictVariable(Y, smth_win, type);
    % inputs: obj: object
    %         Y: firing rate / signal
    %         smth_win(optional): the smoothing window to use on the signal (dfl: 150ms)
    %         type: options 'best', 'mean' (dfl: 'best')
    % output: X: variable encoded/decoded
    %     [pred, Posterior, nPosterior] = obj.predictBayesDecoder(Y, smth_win);
    % The trained object properties:
    %       train_mean: this is necessary to calculate the explained variance in general
    %       CVO: cross-validation object, the sets into which the data was seperated
    %       bins: the values of each bins, this is useful for plotting
    %       relPerformance: This is the (prob(correct postition) - prob(baseline)) ./ (prob(max at that time) - prob(baseline))
    %       confidence: this is the mean of the max at each time point
    %       numBins: the number of bins the data was split into to fit the models
    %       performance: This is explained variance on each cross-validation iteration
    %       meanPerformance: mean of performance
    %       model.EV: The explained variance of each cell on each iteration
    %       model.train_mean: this is necessary to calculate the explained variance in general
    %       model.meanModel: the mean fit model across iterations (units are normalised)
    %       model.bestModel: the model on the iteration that gave the best performance on the population decoding (units are normalised)
    %       model.trained : structure with each of the iterations
    %         respModel: normalised model of each cell
    %         respModel_orig: non-normalised version. Firing rate can be gotten by multiplying with sampleRate
    %
    % Aman Saleem
    % April 2014
    properties
        kfold;       % number of partitions made to train the decoder
        train_mean;  % mean response used for training
        CVO;         % cross-validation structure
        bins;        % the values of the bin limits
        binsB;
        relPerformance;
        confidence;
        numBins;     % number of bins by which the variable is discretized
        numBinsB;     % number of bins by which the variable is discretized
    end
    
    methods
        function obj = bayes2dDecoder(varargin)
            
            pnames = {'type' 'variable' 'variable_range'...
                'smth_win' 'model'...
                'kfold' 'performance'...
                'train_mean' 'CVO' 'sampleRate' 'numBins' 'numBinsB'};
            dflts  = {'2D Bayes' 'PS' []...
                150 []...
                5 []...
                [] [] 60 50 15};
            [obj.type, obj.variable, obj.variable_range, ...
                obj.smth_win, obj.model,...
                obj.kfold, obj.performance,...
                obj.train_mean, obj.CVO, obj.sampleRate, obj.numBins, obj.numBinsB] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        
        function [obj, prediction, X, Posterior, nPosterior] = trainDecoder(obj, X, Y, smth_win, subset);
            if nargin<4
                smth_win = obj.smth_win;
            else
                obj.smth_win = smth_win;
            end
            if nargin<5
                subset = true(size(X));
            end
            [obj, prediction, X, Posterior, nPosterior] = train2dBayesDecoder(obj, X, Y, smth_win, subset);
        end
        
        function [X] = predictVariable(obj, Y, smth_win, type);
            if nargin<3
                smth_win = obj.smth_win;
            end
            if nargin<4
                type = 'best';
            end
            
            t = ones(size(Y,1),1);
            t(isnan(sum(Y,2))) = 0;
            t = t>0;
            X = NaN*ones(size(Y,1),1);
            
            [X(t)] = predictBayesDecoder(obj, Y(t,:), smth_win, type);
        end
        
        function [pred, Posterior, nPosterior] = predictBayesDecoder(obj, Y, smth_win, type);
            for icell = 1:size(Y, 2);
                Y(:,icell) = smthInTime(Y(:,icell), obj.sampleRate, smth_win);
            end
            switch type
                case 'mean'
                    Posterior = calcPosterior(obj.model.meanModel,Y);
                case 'best'
                    Posterior = calcPosterior(obj.model.bestModel,Y);
            end
            
            nPosterior = Posterior./repmat(median(Posterior),size(Posterior,1),1);
            nPosterior = nPosterior./(sum(nPosterior,2)*ones(1,size(nPosterior,2)));
    
            % Converting to log likelihood
            baseline    = 1./(obj.numBins*obj.numBinsB);
            Posterior   = log2(Posterior/baseline);
            nPosterior   = log2(nPosterior/baseline);
            
            Posterior = permute(reshape(Posterior, obj.numBins, obj.numBinsB,[]),[3 1 2]);
        end
        
    end
end

