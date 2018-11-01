function predmed = getCircularMedian(mat)
Prange = size(mat,1);
Xrange = size(mat,2);
postfilt = zeros(Prange,1);
postfilt(1:Prange+1) = (0:(Prange))'/(Prange)*(2*pi);
mat = mat*1000;
mat = [mat;mat(1,:)];
dd0 = circ_dist2(postfilt,postfilt);
n = size(mat,2);
predmed = NaN(1,n);
for i = 1:n
dd = dd0 .* (mat(:,i)*(mat(:,i))');
ddpos = dd;
ddpos(ddpos<0) = 0;
ddneg = dd;
ddneg(ddneg>0) = 0;
m1 = sum(ddpos(1:end,1:end),1);                                                    
m2 = sum(ddneg(1:end,1:end),1);

dm = abs(m1+m2);
predmean = mod(circ_mean(postfilt,mat(:,i)),2*pi);
dm(abs(circ_dist(postfilt,predmean))>pi/2) = inf;
m = min(dm);
idx = find(dm==m,2);
predmed(i) = postfilt(idx);
% predmed(i) = circ_mean(mod(postfilt(idx),pi));
% if abs(circ_dist(predmean,predmed(i))) > abs(circ_dist(predmean,predmed(i)+pi))
%     predmed(i) = mod(predmed(i)+pi,2*pi);
% end
% predmed(i) = predmean;
predmed(i) = mod(predmed(i) + 2*pi, 2*pi);
predmed(i) = (predmed(i)/(2*pi)*Prange)+1;
end
