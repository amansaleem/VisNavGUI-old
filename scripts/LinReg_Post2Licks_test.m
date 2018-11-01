%% instantaneous linear model
kfold = 5;
lambda = 0.1;
idx = EXP.getSubsets(2,2,1,3,false);
licks = smthInTime(es.preRewlick + es.badlick, 60, 150, 'same');
Have = 0;
lickpred = zeros(size(EXP.Bayes.Posterior0,1),1);
allidx = find(idx);
Post = [ones(size(EXP.Bayes.Posterior0,1),1) EXP.Bayes.Posterior0(:,1:end)];%medfilt1(EXP.Bayes.Posterior,16);%
Post(Post<1) = 0;

% Post = 0*Post;
% Xdec = smthInTime(Xpred, 60, 150, 'same');
% for x=0:99
% Post(floor(es.traj)==x,x+1)=1;
% end

for iter = 1:kfold
    if iter < kfold
        testidx = allidx((iter-1)*floor(numel(allidx)/kfold)+1:iter*floor(numel(allidx)/kfold));
    else
        testidx = allidx((iter-1)*floor(numel(allidx)/kfold)+1:end);
    end
    trainidx = allidx(~ismember(allidx,testidx));        
    Posttrain = Post(trainidx,:);
    Posttrain(Posttrain<0) = 0;
    lickstrain = licks(trainidx);
    Hess = Posttrain'*Posttrain;
    Q = (Posttrain'*(lickstrain));
    H = (Hess + lambda*trace(Hess)*eye(size(Posttrain,2)))\Q;
    Have = Have + H/kfold;
    lickpred(testidx) = (Post(testidx,:)*H);
end
lickpred(~idx) = (Post(~idx,:)*Have);
figure;plot(licks)
hold on;plot(lickpred)
figure;plot(xcorr(licks(es.outcome ~= 2),lickpred(es.outcome ~= 2),50,'coeff'))
figure;plot(xcorr(licks(es.outcome == 2),lickpred(es.outcome == 2),50,'coeff'))
figure;plot(Have)

%%instantaneous linear model using position
licks = smthInTime(es.firstgoodlick(tidx) + es.firstbadlick(tidx), 60, 150, 'same');
Post = EXP.Bayes.Posterior0;%medfilt1(EXP.Bayes.Posterior,16);%
Post = Post(tidx,:);
Post = 0*Post;
for x=0:99
Post(floor(es.traj(tidx))==x,x+1)=1;
end
Post = [ones(size(Post,1),1) Post];
Post0 = Post;
licks0 = licks;
Post = Post(1:floor(0.8*size(Post,1)),:);licks = licks(1:floor(0.8*size(licks,1)));
H = (Post'*Post + 10*eye(size(Post,2)))\(Post'*(licks));
Post = Post0;
licks = licks0;
figure;plot(H(2:end))
figure;plot(licks)
hold on;plot(Post*H)
hold on;plot(floor(0.8*size(Post,1))+1:size(Post,1),Post(floor(0.8*size(Post,1))+1:end,:)*H)
hold on;plot(es.traj(tidx)/1000);
figure;scatter(licks(floor(0.8*size(Post,1))+1:end),Post(floor(0.8*size(Post,1))+1:end,:)*H)
set(gca,'Xlim',[0 0.2])

%%linear model with memory
licks = smthInTime(es.firstgoodlick(tidx) + es.firstbadlick(tidx), 60, 150, 'same');
Post = EXP.Bayes.Posterior0(tidx,:);
Posttemp = Post;
ntau = 5;
for tt=1:ntau
Post = [Post circshift(Posttemp,tt,1)];
end
Post(Post<0) = 0;
Post = [ones(size(Post,1),1) Post];
Post0 = Post;
licks0 = licks;
Post = Post(1:floor(0.8*size(Post,1)),:);licks = licks(1:floor(0.8*size(licks,1)));
H = (Post'*Post + 70000*eye(size(Post,2)))\(Post'*(licks));
Post = Post0;
licks = licks0;
figure;
for i = 1:ntau
plot(H(2+(i-1)*size(Posttemp,2):2+i*size(Posttemp,2)-1));
hold on;
end
figure;plot(H(2:end))
figure;plot(licks)
hold on;plot(Post*H)
hold on;plot(floor(0.8*size(Post,1))+1:size(Post,1),Post(floor(0.8*size(Post,1))+1:end,:)*H)
hold on;plot(es.traj(tidx)/1000);
figure;scatter(licks(floor(0.8*size(Post,1))+1:end),Post(floor(0.8*size(Post,1))+1:end,:)*H)
set(gca,'Xlim',[0 0.2],'Ylim',[0 0.2])

%%
