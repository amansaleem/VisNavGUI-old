function [F, c1, c2, con, H] = smoothhist2D_corrected(Z, smth_win, M, bins1, bins2, Fcircular1, Fcircular2,FnormXonly)
if isempty(Z)
    F = nan*ones(M);
    c1 = bins1;
    c2 = bins2;
    con = F;
    H = F;
    return;
end

if nargin < 8
    FnormXonly = true;
end

X = Z(:,1);
Y = Z(:,2);
[F, c1, c2, H] = smoothhist2D_AS(Z, smth_win, M, bins1, bins2, Fcircular1, Fcircular2);%[F, c1, c2] = smoothhist2D(Z, smth_win, M);%

[KX, ~, ~, ~] = smoothhist2D_AS([X ones(size(Y))], [smth_win(1) NaN], [M(1) 2], bins1, 1:2, Fcircular1, false);%smoothhist2D([X ones(size(Y))], [smth_win(2) NaN], [M(1) 2], Fcircular1, false);
[KY, ~, ~, ~] = smoothhist2D_AS([ones(size(Y)) Y], [NaN smth_win(2)], [2 M(2)], 1:2, bins2, false, Fcircular2);%smoothhist2D([ones(size(Y)) Y], [NaN smth_win(1)], [2 M(2)], false, Fcircular2);
KX = KX(1,:);
KY = KY(:,1);

if FnormXonly
    con = repmat(KX,size(F,1),1);
    F = (F ./ con);
    F = F ./ repmat(sum(F,1),size(F,1),1);
    baseline = 1/size(F,2);
    F = F ./ baseline;
else
    con = KY*KX;
    F = (F ./ sqrt(con));
    F = F/sum(F(:));
    baseline = 1/(size(F,1)*size(F,2));
    F = F ./ baseline;
end
% con = con ./ sum(con(:)) * sum(F(:));


% H = (H ./ con);
% imagesc(c1, c2, F)
end