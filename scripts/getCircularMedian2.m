function predmed = getCircularMedian2(mat)
Prange = size(mat,1);
Xrange = size(mat,2);
postfilt = zeros(Prange,1);
postfilt(1:Prange+1) = (0:(Prange))'/(Prange)*(2*pi);
% mat = [mat;mat(1,:)];
n1 = size(mat,1);
n2 = size(mat,2);
mat0 = mat;
m1 = inf(n1,n2);
m2 = inf(n1,n2);
for i = 1:n1
    mat = circshift(mat0,-i+floor(Prange/2),1);
    m1(i,:) = sum(mat(1:floor(Prange/2),:));
    m2(i,:) = sum(mat(floor(Prange/2):end,:));
end
predmean = getCircularAverage(mat0,0,1);
dm = abs(m1-m2);
[~, predmed] = min(dm,[],1);
predmed = predmed(:);
predmed((predmed-predmean) > floor(Prange/4)) = predmed((predmed-predmean) > floor(Prange/4)) - floor(Prange/2);
predmed((predmed-predmean)< - floor(Prange/4)) = predmed((predmed-predmean) < -floor(Prange/4)) + floor(Prange/2);

% for i = 1:n2
%     m1 = inf(1,n1);
%     m2 = inf(1,n1);
%     predmean = mod(circ_mean(postfilt,mat(:,i)),2*pi);
%     krange = find(abs(circ_dist(postfilt,predmean))<pi/2);
%     krange = krange(:)';
%     for k = krange
%         m2(k) = sum(mat(mod((k-500:k)-1,1000)+1,i));
%         m1(k) = sum(mat(mod((k:k+500)-1,1000)+1,i));
%     end
%     dm = abs(m1-m2);
%     % predmean = mod(circ_mean(postfilt,mat(:,i)),2*pi);
%     % dm(abs(circ_dist(postfilt,predmean))>pi/2) = inf;
%     m = min(dm);
%     idx = find(dm(1:end-1)==m,2);
%     predmed(i) = postfilt(idx);
%     predmed(i) = mod(predmed(i) + 2*pi, 2*pi);
%     predmed(i) = (predmed(i)/(2*pi)*Prange)+1;
%     % pause
% end