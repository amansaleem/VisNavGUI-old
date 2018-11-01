function [out, X, distances] = get45Marginal(in, steps)

if nargin<2
    steps = 50;
end

[lenA, lenB] = size(in);
B = repmat( (1:lenA)',[1 lenB]);
A = repmat( (1:lenB) ,[lenA 1]);

distances = sqrt(A.^2+B.^2).*sin( atan(B./A) - pi/4); 
minD = min(distances(:));
maxD = max(distances(:));
stepSize = (maxD-minD)./(steps-1);

X = minD:stepSize:maxD;
for n = 1:length(X)-1
    tmp = (distances>=X(n) & distances<X(n+1));
    out(n) = nanmean(in(tmp));
end
out(n+1) = nanmean(in(distances>X(n)));