function [F, c1, c2, con, H] = smoothhist2D_SVDcorrected(Z, smth_win, M, bins1, bins2,Fcircular1,Fcircular2)
if isempty(Z)
    F = nan*ones(M);
    c1 = bins;
    c2 = bins;
    con = F;
    H = F;
    return;
end

X = Z(:,1);
Y = Z(:,2);
[F, c1, c2, ~] = smoothhist2D_AS(Z, smth_win, M, bins1, bins2, Fcircular1, Fcircular2);

[~,~,~,BestModel,Residual] = MakeSeparable(F);

% imagesc(c1, c2, F)

H = F;
F = Residual;
con = BestModel;

end