function [F, ctrs1, ctrs2, H] = smoothhist2D_AS(X,lambda,nbins, bins1, bins2,Fcircular1,Fcircular2,outliercutoff,plottype)
% SMOOTHHIST2D Plot a smoothed histogram of bivariate data.
%   SMOOTHHIST2D(X,LAMBDA,NBINS) plots a smoothed histogram of the bivariate
%   data in the N-by-2 matrix X.  Rows of X correspond to observations.  The
%   first column of X corresponds to the horizontal axis of the figure, the
%   second to the vertical. LAMBDA is a positive scalar smoothing parameter;
%   higher values lead to more smoothing, values close to zero lead to a plot
%   that is essentially just the raw data.  NBINS is a two-element vector
%   that determines the number of histogram bins in the horizontal and
%   vertical directions.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF) plots outliers in the data as points
%   overlaid on the smoothed histogram.  Outliers are defined as points in
%   regions where the smoothed density is less than (100*CUTOFF)% of the
%   maximum density.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,[],'surf') plots a smoothed histogram as a
%   surface plot.  SMOOTHHIST2D ignores the CUTOFF input in this case, and
%   the surface plot does not include outliers.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF,'image') plots the histogram as an
%   image plot, the default.
%
%   Example:
%       X = [mvnrnd([0 5], [3 0; 0 3], 2000);
%            mvnrnd([0 8], [1 0; 0 5], 2000);
%            mvnrnd([3 5], [5 0; 0 1], 2000)];
%       smoothhist2D(X,5,[100, 100],.05);
%       smoothhist2D(X,5,[100, 100],[],'surf');
%
%   Reference:
%      Eilers, P.H.C. and Goeman, J.J (2004) "Enhancing scaterplots with
%      smoothed densities", Bioinformatics 20(5):623-628.

%   Copyright 2009 The MathWorks, Inc.
%   Revision: 1.0  Date: 2006/12/12
%
%   Requires MATLAB® R14.

if nargin < 8 || isempty(outliercutoff), outliercutoff = .05; end
if nargin < 9, plottype = 'none'; end

edges1 = bins1;
ctrs1 = edges1(1:end) + .5*[diff(edges1) edges1(end)-edges1(end-1)];
edges1 = [-Inf bins1(2:end) Inf];
edges2 = bins2;
ctrs2 = edges2(1:end) + .5*[diff(edges2) edges2(end)-edges2(end-1)];
edges2 = [-Inf bins2(2:end) Inf];

[n,p] = size(X);
bin = zeros(n,2);
% Reverse the columns of H to put the first column of X along the
% horizontal axis, the second along the vertical.
[dum,~,bin(:,2)] = histcounts(X(:,1),edges1);%
[dum,~,bin(:,1)] = histcounts(X(:,2),edges2);%
% [dum,bin(:,2)] = histc(X(:,1),edges1);
% [dum,bin(:,1)] = histc(X(:,2),edges2);
H = accumarray(bin,1,nbins([2 1]));% ./ n;

% Eiler's 1D smooth, twice
if Fcircular1
    repX = 3;
else
    repX = 1;
end
if Fcircular2
    repY = 3;
else
    repY = 1;
end
K = repmat(H,repY,repX);
F = special_smooth2D(K,[lambda(2)/repY/size(H,1) lambda(1)/repX/size(H,2)]);
% if ~isnan(lambda(1))
%     G = smooth1D(K,lambda(1));
% else
%     G = K;
% end
% if ~isnan(lambda(2))
%     F = smooth1D(G',lambda(2))';
% else
%     F = G;
% end
if Fcircular1
    F = F(:,floor(size(F,2)/repX)+1:2*floor(size(F,2)/repX));
end
if Fcircular2
    F = F(floor(size(F,1)/repY)+1:2*floor(size(F,1)/repY),:);
end
% % An alternative, using filter2.  However, lambda means totally different
% % things in this case: for smooth1D, it is a smoothness penalty parameter,
% % while for filter2D, it is a window halfwidth
% F = filter2D(H,lambda);

relF = F./max(F(:));
if outliercutoff > 0
    outliers = (relF(nbins(2)*(bin(:,2)-1)+bin(:,1)) < outliercutoff);
end

nc = 256;
switch plottype
    case 'surf'
        colormap(hot(nc));
        surf(ctrs1,ctrs2,F,'edgealpha',0);
    case 'image'
        colormap(hot(nc));
        image(ctrs1,ctrs2,floor(nc.*relF) + 1);
        hold on
        % plot the outliers
        if outliercutoff > 0
            plot(X(outliers,1),X(outliers,2),'.','MarkerEdgeColor',[.8 .8 .8]);
        end
        %     % plot a subsample of the data
        %     Xsample = X(randsample(n,n/10),:);
        %     plot(Xsample(:,1),Xsample(:,2),'bo');
        hold off
    case 'none'
end
end

%-----------------------------------------------------------------------------
function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
% This is a better solution, but takes a bit longer for n and m large
% opts.RECT = true;
% D1 = [diff(E,1); zeros(1,n)];
% D2 = [diff(D1,1); zeros(1,n)];
% Z = linsolve([E; 2.*sqrt(lambda).*D1; lambda.*D2],[Y; zeros(2*m,n)],opts);
end

%-----------------------------------------------------------------------------
function Z = filter2D(Y,bw)
z = -1:(1/bw):1;
k = .75 * (1 - z.^2); % epanechnikov-like weights
k = k ./ sum(k);
Z = filter2(k'*k,Y);
end

function mat_o = special_smooth2D(mat_i,win)
mat_o = mat_i;
mat_o(isnan(mat_i)) = 0;
if sum(isnan(win)) > 0
    if isnan(win(1)) && ~isnan(win(2))
        for i = 1:size(mat_o,1)
            mat_o(i,:) = special_smooth_1d(mat_o(i,:), win(2), [], size(mat_i,2));
        end
    end
    if isnan(win(2)) && ~isnan(win(1))
        for j = 1:size(mat_o,2)
            mat_o(:,j) = special_smooth_1d(mat_o(:,j), win(1), [], size(mat_i,1));
        end
    end
else
    mat_o = special_smooth_2d(mat_o, win, [], [], [size(mat_i,1) size(mat_i,2)]);
end
mat_o(isnan(mat_i)) = NaN;
end

