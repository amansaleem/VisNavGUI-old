function beh = BehaviorSummary(es,F1stlick,roomlength)
%lick distribution
gainval = [0.4 0.5 0.6];
outcomeval = [0 2 3];
if nargin < 3
    roomlength = 100;%50;%round(max(es.trajPercent));
end
% idx1 = (es.trajPercent < floor(roomlength/2));
% idx2 = (es.trajPercent >= floor(roomlength/2));
% es.trajPercent(idx1) = es.trajPercent(idx1) + floor(roomlength/2);
% es.trajPercent(idx2) = es.trajPercent(idx2) - floor(roomlength/2);

% es.trajPercent = es.trajPercent.*es.gain/mode(es.gain);
es.trajPercent(es.trajPercent < 0) = min(es.trajPercent(es.trajPercent > 0));

if F1stlick
    goodlicks = es.firstgoodlick;
    badlicks = es.firstbadlick;
else
    goodlicks = es.preRewlick;%es.firstgoodlick;%
    badlicks = es.badlick;
end


% if outcomeval == 2
%     lim1 = round(max(es.trajPercent) * 0.3);
%     lim2 = round(max(es.trajPercent) * 0.9);
%     idx = find(es.reward == 2);
%     for i = 1:numel(idx)-1
% %         licks(idx(i)-1:(idx(i)-1+find(es.trajPercent((idx(i)-1):idx(i+1))>lim1 & es.trajPercent((idx(i)-1):idx(i+1))<lim2,1,'first')-1)) = 0;
%     end
% end

sessionlist = unique(es.dayID);
% sessionlist(1:40) = [];
nbins = roomlength;%(floor(max(es.trajPercent)/50) + 1)*50;
speed_th = 5;

beh.lickMat = cell(numel(gainval),numel(outcomeval));
for g = 1:numel(gainval)
    for o = 1:numel(outcomeval)
        beh.lickMat{g,o} = zeros(numel(sessionlist), 2*nbins + 1);
    end
end
beh.delatgain = zeros(3,numel(sessionlist));
beh.nbtrials = zeros(3,numel(sessionlist));
beh.correct = zeros(3,numel(sessionlist));
beh.ntrialchange = zeros(1,numel(sessionlist));
beh.maxbadlicks = zeros(1,numel(sessionlist));
beh.recday = zeros(1,numel(sessionlist));

beh.perf = cell(1,numel(gainval));
beh.perfSE = cell(1,numel(gainval));
beh.lickshift = cell(1,numel(gainval));
for g = 1:numel(gainval)
    beh.perf{g} = zeros(1,numel(sessionlist));
    beh.perfSE{g} = zeros(1,numel(sessionlist));
    beh.lickshift{g} = zeros(1,numel(sessionlist));
end
beh.speedAve = cell(numel(gainval),numel(outcomeval));
beh.speedStd = cell(numel(gainval),numel(outcomeval));
for g = 1:numel(gainval)
    for o = 1:numel(outcomeval)
        beh.speedAve{g,o} = zeros(numel(sessionlist), 2*nbins + 1);
        beh.speedStd{g,o} = zeros(numel(sessionlist), 2*nbins + 1);
    end
end

for i = 1:numel(sessionlist)
    for g = 1:numel(gainval)
        for o = 1:numel(outcomeval)
            if outcomeval(o) == 2
                idx = (es.blanks == 0 & es.gain == gainval(g) & es.outcome == outcomeval(o) &  es.dayID == sessionlist(i));% & es.smthBallSpd > speed_th);
            elseif outcomeval(o) == 0
                idx = (es.blanks == 0 & es.gain == gainval(g) & es.outcome == outcomeval(o) &  es.dayID == sessionlist(i));% & es.smthBallSpd > speed_th);
            elseif outcomeval(o) == 3
                idx = (es.blanks == 0 & es.gain == gainval(g) &  es.dayID == sessionlist(i));% & es.smthBallSpd > speed_th);
            end
            licks = smthInTime(double(goodlicks | badlicks), 60, 0,[],[],'boxcarsum');
            if sum(idx) > 0 && sum(licks)>0
%                [f,x] = ecdf(es.trajPercent(idx));
%                [n,~] = ecdfhist(f,x,0:nbins-1);
%                 beh.lickMat{g,o}(i,:) = smooth(n,3)/sum(smooth(n,3));                
                [lickdistri,x] = fast1Dmap(es.trajPercent(idx),licks(idx),2,60);
                lickdistri = interp1([x x(end) + x(2) - x(1)],[lickdistri' lickdistri(1)],linspace(x(1),(x(end) + x(2) - x(1)),2*nbins + 1),'spline');
                beh.lickMat{g,o}(i,1:min(numel(beh.lickMat{g,o}),numel(lickdistri))) = lickdistri(1:min(numel(beh.lickMat{g,o}),numel(lickdistri)))./sum(lickdistri);
            end
            if outcomeval(o) == 2 || outcomeval(o) == 0
                idx = (es.outcome == outcomeval(o)  & es.blanks == 0 & es.gain == gainval(g) & es.smthBallSpd > speed_th & es.dayID == sessionlist(i));
            elseif outcomeval(o) == 3
                idx = (ismember(es.outcome,[0 2]) & es.blanks == 0 & es.gain == gainval(g) & es.smthBallSpd > speed_th & es.dayID == sessionlist(i));
            end
            if sum(idx) > 0
                [speedpro,x,~] = easyPlaceProfile(es.trajPercent(idx),es.ballspeed(idx),2);
                speedpro = interp1([x x(end) + x(2) - x(1)],[speedpro' speedpro(1)],x:0.5:(x(end) + x(2) - x(1)),'spline');
                beh.speedAve{g,o}(i,1:min(numel(beh.speedAve{g,o}),numel(speedpro))) = speedpro;  
            end
        end
    end
%     for xx = 0:nbins-1        
%         beh.speedAve(i, xx + 1) = mean(es.smthBallSpd(idx));
%         beh.speedStd(i, xx + 1) = std(es.smthBallSpd(idx));
%     end
    
    beh.delatgain(1,i) = sum(ismember(es.gain(es.dayID == sessionlist(i) & es.outcome == 2),0.4))/sum(ismember(es.gain(es.dayID == sessionlist(i) & es.outcome == 2),[0.4 0.5 0.6]));
    beh.delatgain(2,i) = sum(ismember(es.gain(es.dayID == sessionlist(i) & es.outcome == 2),0.5))/sum(ismember(es.gain(es.dayID == sessionlist(i) & es.outcome == 2),[0.4 0.5 0.6]));
    beh.delatgain(3,i) = sum(ismember(es.gain(es.dayID == sessionlist(i) & es.outcome == 2),0.6))/sum(ismember(es.gain(es.dayID == sessionlist(i) & es.outcome == 2),[0.4 0.5 0.6]));
    
    for g = 1:numel(gainval)
        beh.nbtrials(g,i) = numel(unique(es.trialID(es.dayID == sessionlist(i) & es.gain == gainval(g) & es.outcome < 5)));
        correcttrials = unique(es.trialID(es.outcome == 2 & es.dayID == sessionlist(i) & es.gain == gainval(g)));
        incorrecttrials = unique(es.trialID(es.outcome == 0 & es.dayID == sessionlist(i) & es.gain == gainval(g)));
        if numel(unique([correcttrials;incorrecttrials])) > 0
            beh.correct(g,i) = sum(~ismember(correcttrials,incorrecttrials))/numel(unique([correcttrials;incorrecttrials]));
        end
    end
    beh.ntrialchange(i) = mean(es.nTrialChange(es.dayID == sessionlist(i)));
    beh.maxbadlicks(i) = mean(es.maxbadlicks(es.dayID == sessionlist(i)));
    beh.badtrials(i) = numel(unique(es.trialID(es.dayID == sessionlist(i) & es.outcome == 0)));%divided by 2 cause it's in reward counts here
    beh.recday(i) = max(es.RecDay(es.dayID == sessionlist(i))); 
    
    kfold = 20;
    for g = 1:numel(gainval)
        tidx = find(es.gain == gainval(g) & ismember(es.outcome,1:4) & ~es.blanks & es.smthBallSpd >= 5 & es.dayID == sessionlist(i));
        beh.perf{g}(i) = sum(es.firstgoodlick(tidx))/(sum(es.firstgoodlick(tidx)) + sum(es.firstbadlick(tidx)));
        for k = 1:kfold
            outtidx = tidx((k-1)*floor(numel(tidx)/kfold)+1:k*floor(numel(tidx)/kfold));
            ktidx = tidx(~ismember(tidx,outtidx));
            perffold = sum(es.firstgoodlick(ktidx))/(sum(es.firstgoodlick(ktidx)) + sum(es.firstbadlick(ktidx)));
            beh.perfSE{g}(i) = beh.perfSE{g}(i) + (kfold - 1)/kfold*((perffold - beh.perf{g}(i))^2);
        end
        beh.perfSE{g}(i) = sqrt(beh.perfSE{g}(i));
        
        tidx = (es.gain == gainval(g) & ismember(es.outcome,2) & es.blanks == 0 & es.dayID == sessionlist(i));
        if sum(tidx) > 0 && beh.perf{g}(i) > 0.65
            tidxref = (es.gain == 0.5 & ismember(es.outcome,2) & es.blanks == 0 & es.dayID == sessionlist(i));
            licks = smthInTime(double(goodlicks | badlicks), 60, 0,[],[],'boxcarsum');
            [lickdistri,x] = fast1Dmap(es.trajPercent(tidx),licks(tidx),2,60);
            lickdistri = interp1([x x(end) + x(2) - x(1)],[lickdistri' lickdistri(1)],x:0.5:(x(end) + x(2) - x(1)),'spline');
            [lickdistriref,x] = fast1Dmap(es.trajPercent(tidxref),licks(tidxref),2,60);
            lickdistriref = interp1([x x(end) + x(2) - x(1)],[lickdistriref' lickdistriref(1)],x:0.5:(x(end) + x(2) - x(1)),'spline');
            x = x:0.5:(x(end) + x(2) - x(1));
            lickshift = x(find(lickdistri == max(lickdistri),1,'first')) - x(find(lickdistriref == max(lickdistriref),1,'first'));            
            if lickshift < -max(x)/2
                lickshift = lickshift + max(x);
            end
            if lickshift > max(x)/2
                lickshift = lickshift - max(x);
            end
            lickshift(abs(lickshift) > 15) = nan;
            if lickshift < 0 && g == 3
                keyboard
            end
            beh.lickshift{g}(i) = lickshift;            
        else
            beh.lickshift{g}(i) = nan;
        end        
    end
end

x = linspace(x(1),(x(end) + x(2) - x(1)),2*nbins + 1);

figure;
nsubplot = 3;
day = numel(beh.recday);
alphaidx = 0.2;
h = subplot(2,nsubplot,1);
ciplot(h, 1:numel(beh.perf{gainval == 0.5}),beh.perf{gainval == 0.5},beh.perfSE{gainval == 0.5},alphaidx);
hold on;plot(beh.perf{gainval == 0.5},'k');
set(gca,'PlotBoxaspectratio',[1 1 1])

subplot(2,nsubplot,2);
GMModel = fitgmdist([beh.lickshift{3}(~isnan(beh.lickshift{3}) & ~isnan(beh.lickshift{1}));beh.lickshift{1}(~isnan(beh.lickshift{3}) & ~isnan(beh.lickshift{1}))]',1);
hold on;
x1 = -20:0.2:20;
y1 = -20:0.2:20;
[X,Y] = meshgrid(x1,y1);
Z = mvnpdf([X(:) Y(:)],GMModel.mu,GMModel.Sigma);
Z = reshape(Z,length(X),length(Y));
v = [0.5*max(Z(:)),0.5*max(Z(:))];
[C,h] = contourf(X,Y,Z,v,'LineStyle','none','FaceColor',[0.8 0.8 0.8]);
scatter(beh.lickshift{3},beh.lickshift{1},'.k');
scatter(beh.lickshift{3}(day),beh.lickshift{1}(day),'.r');

set(gca,'Xlim',[-20 20],'Ylim',[-20 20]);
hold on;plot([-20 20],[0 0],'k');
hold on;plot([0 0],[-20 20],'k');
set(gca,'PlotBoxaspectratio',[1 1 1])

subplot(2,nsubplot,3);
bar(1:3,[beh.perf{2}(day) beh.perf{1}(day) beh.perf{3}(day)],'k');
hold on;errorbar(1:3,[beh.perf{2}(day) beh.perf{1}(day) beh.perf{3}(day)],[beh.perfSE{2}(day) beh.perfSE{1}(day) beh.perfSE{3}(day)],'.');
set(gca,'PlotBoxaspectratio',[1 1 1])

subplot(2,nsubplot,4);
color{1} = 'c';color{2} = [0.5 0.5 0.5];color{3} = 'm';
maxP = 0;
for g = 1:numel(gainval)
    hold on;plot(x,beh.lickMat{g,2}(day,:),'color',color{g});
    maxP = max(maxP,max(beh.lickMat{g,2}(day,:)));
end
set(gca,'Xlim',[min(x) max(x)],'Ylim',[0 1.1*maxP]);
set(gca,'PlotBoxaspectratio',[1 1 1]);

subplot(2,nsubplot,5);
color{1} = 'c';color{2} = [0.5 0.5 0.5];color{3} = 'm';
maxP = 0;
for g = 1:numel(gainval)
    hold on;plot(x,beh.lickMat{g,3}(day,:),'color',color{g});
    maxP = max(maxP,max(beh.lickMat{g,3}(day,:)));
end
set(gca,'Xlim',[min(x) max(x)],'Ylim',[0 1.1*maxP]);
set(gca,'PlotBoxaspectratio',[1 1 1])

subplot(2,nsubplot,6);
color{1} = 'c';color{2} = [0.5 0.5 0.5];color{3} = 'm';
maxP = 0;
for g = 1:numel(gainval)
    hold on;plot(x,beh.speedAve{g,2}(day,:),'color',color{g});
    maxP = max(maxP,max(beh.speedAve{g,2}(day,:)));
end
set(gca,'Xlim',[min(x) max(x)],'Ylim',[0 1.1*maxP]);
set(gca,'PlotBoxaspectratio',[1 1 1])
