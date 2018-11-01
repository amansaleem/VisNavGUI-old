function [EV, corrError, L, Q, train_mean]      = calCrossValExpVar(train, stest, spred, test, pred)

% Usage: [EV, corrError, L, Q, train_mean]      = calCrossValExpVar(train, stest, spred, test, pred)
% function to calculate the cross-validation correlation and the explained variance
% between the data and a prediction

avar = nansum((stest-nanmean(train)).^2)./sum(~isnan(stest));
spred = (spred - nanmean(spred) + nanmean(train));

train_mean = nanmean(train);
MSE  = nansum((stest-spred).^2)./sum(~isnan(stest) & ~isnan(spred));
EV = 1 - (MSE/avar);
if nargout>1
    corrError = corr(spred,stest);
    index = (~isnan(test) & ~isnan(pred));
    test = test(index);
    pred = pred(index);
    T = sum(index);
    Q = (2*(test'*pred) - sum(pred.^2))./T;
    L = (sum(log(eps+(test.*pred))) - sum(pred))./T;
end


%     function [EV, corrError, L, Q]                  = calculateModelError(train_mean, stest, spred, test, pred)
%         % function to calculate the cross-validation correlation and the explained variance
%         % between the data and a prediction
%         avar = nansum((stest-train_mean).^2)./sum(~isnan(stest));
%         spred = (spred - nanmean(spred) + train_mean);
%
%         MSE  = nansum((stest-spred).^2)./sum(~isnan(stest) & ~isnan(spred));
%         EV = 1 - (MSE/avar);
%         if nargout>1
%             corrError = corr(spred,stest);
%             index = (~isnan(test) & ~isnan(pred));
%             test = test(index);
%             pred = pred(index);
%             T = sum(index);
%             Q = (2*(test'*pred) - sum(pred.^2))./T;
%             L = (sum(log(eps+(test.*pred))) - sum(pred))./T;
%         end
%
%     end