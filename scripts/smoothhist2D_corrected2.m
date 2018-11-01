function [F, c1, c2, con, H] = smoothhist2D_corrected2(Z, smth_win, M, bins1, bins2, Fcircular1, Fcircular2)
if isempty(Z)
    F = nan*ones(M);
    c1 = bins1;
    c2 = bins2;
    con = F;
    H = F;
    return;
end

X = Z(:,1);
Y = Z(:,2);
[F, c1, c2, H] = smoothhist2D_AS(Z, smth_win, M, bins1, bins2, Fcircular1, Fcircular2);%[F, c1, c2] = smoothhist2D(Z, smth_win, M);%

KX = smoothhist2D([X ones(size(Y))], smth_win, [M(1) 2], Fcircular1, false);
KY = smoothhist2D([ones(size(Y)) Y], smth_win, [2 M(2)], false, Fcircular2);
KX = KX(2,:);
KY = KY(:,2);

con = KY*KX;
% con = con ./ sum(con(:)) * sum(F(:));

% F = (F ./ con);
% F = F ./ repmat(sum(F,1),size(F,1),1);

%commented out just for testing purpose
% F = F / sum(F(:));
baseline = 1/(numel(F(:)));
F = F ./ baseline;

% H = (H ./ con);
% imagesc(c1, c2, F)
end