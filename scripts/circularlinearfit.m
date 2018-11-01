function [s,phi0,rho,s_SE,phi0_SE,rho_SE] = circularlinearfit(Phi,X,map,kfold)
% a = 1/n * sum(cos(Phi - 2*pi*a*X))^2;
% b = 1/n * sum(sin(Phi - 2*pi*a*X))^2;
% R = sqrt(a + b);
if nargin < 3
    map = ones(size(Phi));
end
if isempty(map)
    map = ones(size(Phi));
end
if nargin < 4 || max(map(:)) ~= min(map(:))
    kfold = 1;
end

[s,phi0,rho,R] = circlincorr(Phi,X,map);

if kfold > 1
    s_std = 0;
    phi0_std = 0;
    rho_std = 0;
    s_iter = zeros(1,kfold);
    phi0_iter = zeros(1,kfold);
    rho_iter = zeros(1,kfold);
    n = numel(Phi);
    idxshf = randperm(n);
    if n <= 10*kfold
        kfold = n;
    end
    for kiter = 1:kfold
        idx_test = idxshf((kiter-1)*round(n/kfold)+1:min(kiter*round(n/kfold),end));
        idx_train = idxshf(~ismember(idxshf,idx_test));
        Phi_iter = Phi(idx_train);
        X_iter = X(idx_train);
        map_iter = ones(size(Phi_iter));
        [s_iter(kiter),phi0_iter(kiter),rho_iter(kiter)] = circlincorr(Phi_iter,X_iter,map_iter,R);
        %[s_iter(kiter),phi0_iter(kiter),rho_iter(kiter)] = circlincorr(Phi_iter,X_iter,map_iter,R);
    end
    for kiter = 1:kfold
        s_fold = sum(~isnan(s_iter));
        s_SE = (s_fold- 1)/s_fold*nansum((s_iter - s).^2);
        phi0_fold = sum(~isnan(phi0_iter));
        phi0_SE = (phi0_fold- 1)/phi0_fold*nansum(circ_dist(phi0_iter,phi0).^2);
        rho_fold = sum(~isnan(rho_iter));
        rho_SE = (rho_fold- 1)/rho_fold*nansum((rho_iter - rho).^2);
    end
    s_SE = sqrt(s_SE);
    phi0_SE = sqrt(phi0_SE);
    rho_SE = sqrt(rho_SE);
else
    s_SE = NaN;
    phi0_SE = NaN;
    rho_SE = NaN;
end
end

function [s,phi0,rho,R] = circlincorr(Phi,X,map,Rref)
if nargin < 4
    Rref = [];
end
n = sum(map);
R = NaN(361,1);
A = (-90:90)/360*2*pi;
for i = 1:numel(A)
    a = A(i);
    R(i) = sqrt((1/n * sum(map.*cos(Phi - 2*pi*a*X)))^2 + (1/n * sum(map.*sin(Phi - 2*pi*a*X)))^2);
end
if ~isempty(Rref)
%     Rcorr = xcorr(R,Rref,90,'coeff');
%     [~,ishift] = max(Rcorr);
%     ishift = ishift-(90+1);
    [~,imaxref] = max(Rref);
%     imax = imaxref + ishift;
%     s = A(imax);
    
    R(1:(floor(numel(A)/2) - (abs(floor(numel(A)/2)-imaxref)+45))) = 0;
    R((floor(numel(A)/2) + (abs(floor(numel(A)/2)-imaxref)+45)):end) = 0;
    [~,imax] = max(R);
    s = A(imax);
else
    [~,imax] = max(R);
    s = A(imax);
end


phi0_num = sum(map.*sin(Phi - 2*pi*s*X));
phi0_den = sum(map.*cos(Phi - 2*pi*s*X));
phi0 = atan2(phi0_num,phi0_den);

theta = 2*pi*abs(s)*X;
theta_num = sum(map.*sin(theta));
theta_den = sum(map.*cos(theta));
theta_bar = atan2(theta_num,theta_den);

Phi_num = sum(map.*sin(Phi));
Phi_den = sum(map.*cos(Phi));
Phi_bar = atan2(Phi_num,Phi_den);

rho_num = sum(map.*sin(Phi - Phi_bar).*sin(theta - theta_bar));
rho_den = sqrt(sum(map.*sin(Phi - Phi_bar).^2)*sum(map.*sin(theta - theta_bar).^2));
rho = rho_num/rho_den;
end