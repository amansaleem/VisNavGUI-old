function [B,bins] = normalise1var(A, aGrid, c, minmax)

% Usage: [B bins] = normalise1var(A, aGrid, c);
% 
% [B] = normalise1var(A) or [B] = normalise1var(A,[],c);
% Normalize the variable so it goes from 0 to 1
% 
% [B bins] = normalise1var(A, aGrid, c);
% Normalize the variable so it is discretized into 'aGrid' bins with centres 'bins'
% 'c' is optional logical array, given this the discretizasion is limited to the specified subset 

if nargin<3
    c = true(size(A));
else
    c = logical(c);
end
if isempty(c)
    c = true(size(A));
end

if nargin < 4
    maxA = max(A(c & ~isnan(A) & ~isinf(A)));
    minA = min(A(c & ~isnan(A) & ~isinf(A)));
else
    maxA = minmax(2);
    minA = minmax(1);
end

if sum(~isnan(A)) == 0
    B = A;
    bins = [];
elseif nargin<2 || isempty(aGrid)
    B = (A-minA)/(maxA-minA);
    bins = [];
else
    bins = minA:((maxA-minA)/(aGrid+eps)):maxA;
    A = (A-minA)/(maxA-minA);
    B = 1+floor(aGrid*A/(1+eps));%1+floor(aGrid*A/(1+eps));
    B(B<1) = 1;
    B(B>aGrid) = aGrid;
end