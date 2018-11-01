classdef TbayesDecoder2D < TDecoder
    %     Create a Decoder class object of type 'linear'. The decoder can be trained using the
    %     function 'trainDecoder' and used to decode with the function
    %     'predictVariable'
    %     Usage:
    % obj = bayesDecoder;
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
    % October 2013
    properties
        kfold;       % number of partitions made to train the decoder
        train_mean;  % mean response used for training
        CVO;         % cross-validation structure
        binsX;        % the values of the bin limits
        binsY;
        relPerformance;
        confidence;
        numBinsX;     % number of bins by which the variable is discretized
        numBinsY;
        FoptiSmooth;
        Fcircular;
        Flookuptable;
    end
    
    methods
        function obj = TbayesDecoder2D(varargin)
            
            pnames = {'type' 'variable' 'variable_range'...
                'Tsmth_win' 'Xsmth_win' 'Ysmth_win' 'model'...
                'kfold' 'performance'...
                'train_mean' 'CVO' 'sampleRate' 'numBinsX' 'numBinsY' 'FoptiSmooth' 'Fcircular' 'Flookuptable'};
            dflts  = {'bayes' 'P' []...
                150 1 1 []...
                5 []...
                [] [] 60 100 10 true true false};
            [obj.type, obj.variable, obj.variable_range, ...
                obj.Tsmth_win, obj.Xsmth_win, obj.Ysmth_win, obj.model,...
                obj.kfold, obj.performance,...
                obj.train_mean, obj.CVO, obj.sampleRate, obj.numBinsX, obj.numBinsY, obj.FoptiSmooth, obj.Fcircular, obj.Flookuptable] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        
       function [obj, prediction, X, Posterior, nPosterior, nonNormPosterior] = trainDecoder(obj, X, Y, Z, T, Xsmth_win, Ysmth_win, subset)
            if nargin<7
                Xsmth_win = obj.Xsmth_win;
                Ysmth_win = obj.Ysmth_win;
            else
                obj.Xsmth_win = Xsmth_win;
                obj.Ysmth_win = Ysmth_win;
            end
            if nargin<9
                subset = true(size(X));
            end
            obj.CVO = [];
            [obj, prediction, X, Posterior, ~] = trainBayesDecoder2D(obj, X, Y, Z, T, Xsmth_win, Ysmth_win, subset);
            
            nPosterior = [];
            nonNormPosterior = [];
        end
        
        function [prediction] = predictVariable(obj, X, Y, type)
            if nargin<4
                type = 'best';
            end
            
            t = ones(size(Y,1),1);
            t(isnan(sum(Y,2))) = 0;
            t = t>0;
            prediction = NaN*ones(size(Y,1),1);
            
            [prediction(t)] = predictBayesDecoder(obj, X, Y(t,:), type);
        end
        
        function [pred, X, Posterior, nPosterior, nonNormPosterior] = predictBayesDecoder(obj, X, Y, Z, T, type)
            baseline = 1./(obj.numBinsX*obj.numBinsY);
            if ~ischar(type)
                custommodel = type;
            end
            switch type
                case 'mean'
                    [Posterior, ~] = calcPosterior2D(obj.model.meanModel, Z, T, baseline);
%                 case 'best'
%                     [Posterior, nonNormPosterior] = calcPosterior(obj.model.bestModel, Y, T, baseline);
%                 case 'All'
%                     Posterior = [];
%                     nonNormPosterior = [];
%                     for iter = 1:numel(obj.model.trained)
%                         [Posteriortemp, nonNormPosteriortemp] = calcPosterior(obj.model.trained(iter).respModel, Y, T, baseline);
%                         if isempty(Posterior)
%                             Posterior = zeros(size(Posteriortemp,1),size(Posteriortemp,2),numel(obj.model.trained));
%                             nonNormPosterior = zeros(size(nonNormPosteriortemp,1),size(nonNormPosteriortemp,2),numel(obj.model.trained));
%                         end
%                         Posterior(:,:,iter) = Posteriortemp;
%                         nonNormPosterior(:,:,iter) = nonNormPosteriortemp;
%                     end
%                 case 'custom'
%                     [Posterior, nonNormPosterior] = calcPosterior(custommodel, Y, T, baseline);
            end
            
            X = [X(:) Y(:)];
            nPosterior = [];
            nonNormPosterior = [];
    
            % Converting to log likelihood                        
%             Posterior   = log2(Posterior);
%             nPosterior   = log2(nPosterior);
            pred = zeros(size(Posterior,1),2);
            [~, predtempX] = max(squeeze(sum(Posterior,2)),[],2);
            [~, predtempY] = max(squeeze(sum(Posterior,3)),[],2);
            pred(:,1) = predtempX(:);
            pred(:,2) = predtempY(:);
            
%             [X, ~] = normalise1var(X, obj.numBins);
        end        
    end
end

