classdef TbayesDecoder < TDecoder
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
        bins;        % the values of the bin limits
        relPerformance;
        confidence;
        numBins;     % number of bins by which the variable is discretized
        FoptiSmooth;
        Fcircular;
        Flookuptable;
    end
    
    methods
        function obj = TbayesDecoder(varargin)
            
            pnames = {'type' 'variable' 'variable_range'...
                'Tsmth_win' 'Xsmth_win' 'model'...
                'kfold' 'performance'...
                'train_mean' 'CVO' 'sampleRate' 'numBins' 'FoptiSmooth' 'Fcircular' 'Flookuptable'};
            dflts  = {'bayes' 'P' []...
                150 1 []...
                5 []...
                [] [] 60 100 true true false};
            [obj.type, obj.variable, obj.variable_range, ...
                obj.Tsmth_win, obj.Xsmth_win, obj.model,...
                obj.kfold, obj.performance,...
                obj.train_mean, obj.CVO, obj.sampleRate, obj.numBins, obj.FoptiSmooth, obj.Fcircular, obj.Flookuptable] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        
        function [obj, prediction, X, Posterior, nPosterior, nonNormPosterior] = trainDecoder(obj, X, Y, Yfilt, T, Xsmth_win, subset)
            if nargin<6
                Xsmth_win = obj.Xsmth_win;
            else
                obj.Xsmth_win = Xsmth_win;
            end
            if nargin<7
                subset = true(size(X));
            end
            obj.CVO = [];
            baseline = 1./obj.numBins;
            if ~obj.Flookuptable
                [obj, prediction, X, Posterior, nonNormPosterior] = trainBayesDecoder(obj, X, Y, Yfilt, T, Xsmth_win, subset);
            else
                [obj, prediction, X, Posterior, nonNormPosterior] = trainBayesDecoder2(obj, X, Y, T, Xsmth_win, subset);
            end
            
            nPosterior = Posterior./repmat(median(Posterior),size(Posterior,1),1);
            nPosterior = nPosterior./(sum(nPosterior,2)*ones(1,size(nPosterior,2)));
            nPosterior   = nPosterior/baseline;
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
        
        function [pred, X, Posterior, nPosterior, nonNormPosterior] = predictBayesDecoder(obj, X, Y, T, type)
            baseline = 1./obj.numBins;
            if ~ischar(type)
                custommodel = type;
            end
            switch type
                case 'mean'
                    if ~obj.Flookuptable
                        [Posterior, nonNormPosterior] = calcPosterior(obj.model.meanModel, Y, T, baseline);
                    else
                        [Posterior, nonNormPosterior] = calcPosterior2(obj.model.meanModel, Y, T, baseline);
                    end
                case 'best'
                    [Posterior, nonNormPosterior] = calcPosterior(obj.model.bestModel, Y, T, baseline);
                case 'All'
                    Posterior = [];
                    nonNormPosterior = [];
                    for iter = 1:numel(obj.model.trained)
                        [Posteriortemp, nonNormPosteriortemp] = calcPosterior(obj.model.trained(iter).respModel, Y, T, baseline);
                        if isempty(Posterior)
                            Posterior = zeros(size(Posteriortemp,1),size(Posteriortemp,2),numel(obj.model.trained));
                            nonNormPosterior = zeros(size(nonNormPosteriortemp,1),size(nonNormPosteriortemp,2),numel(obj.model.trained));
                        end
                        Posterior(:,:,iter) = Posteriortemp;
                        nonNormPosterior(:,:,iter) = nonNormPosteriortemp;
                    end
                case 'custom'
                    [Posterior, nonNormPosterior] = calcPosterior(custommodel, Y, T, baseline);
            end
            nPosterior = zeros(size(Posterior));
            for iter = 1:size(Posterior,3)
                nPosterior(:,:,iter) = Posterior(:,:,iter)./repmat(median(Posterior(:,:,iter)),size(Posterior(:,:,iter),1),1);
                nPosterior(:,:,iter) = nPosterior(:,:,iter)./(sum(nPosterior(:,:,iter),2)*ones(1,size(nPosterior(:,:,iter),2)));
                nPosterior(:,:,iter)   = nPosterior(:,:,iter)/baseline;
            end
    
            % Converting to log likelihood                        
%             Posterior   = log2(Posterior);
%             nPosterior   = log2(nPosterior);
            pred = zeros(size(Posterior,1),size(Posterior,3));
            for iter = 1:size(Posterior,3)
                [~, predtemp] = max(Posterior(:,:,iter),[],2);
                pred(:,iter) = predtemp';
            end
            
%             [X, ~] = normalise1var(X, obj.numBins);
        end        
    end
end

