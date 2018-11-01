function plotBehaviorSummary(beh,roomlength)
outcomeval = [1 2 3];
figure;
if nargin < 2
    roomlength = 100;%(2nd cell is good trials)
end
x = -floor(roomlength/2):floor(roomlength/2)-1;%-25:24;%
nlastexp = 4;

if nargin < 2
    outcomeval = 2;%(2nd cell is good trials)
end

for i = 1:4
    for j = 1:5
        h(i,j) = subplot(4,5,(i-1)*5+j);
    end
end
color{1} = 'b';
color{2} = 'm';
color{3} = 'r';

for g = 1:size(beh.lickMat,1)
    for o = 1:numel(outcomeval)
        recday1 = find(beh.recday == 1, 1, 'first'); 
        imagesc(x, [], beh.lickMat{g,outcomeval(3 - o+1)},'parent',h(3-(g-1),o));
        hold(h(3-(g-1),o),'on');
        plot(h(3-(g-1),o), [x(1) x(end)],[recday1-nlastexp recday1-nlastexp],'r');
        plot(h(3-(g-1),o), [x(1) x(end)],[recday1 recday1],'r');
        hold(h(4,o),'on');
        plot(h(4,o), x, mean(beh.lickMat{g,outcomeval(3 - o+1)}(recday1-nlastexp:recday1,:))/sum(mean(beh.lickMat{g,outcomeval(3 - o+1)}(recday1-nlastexp:recday1,:))),color{g});
        set(h(3-(g-1),o),'Clim',[0 0.1]);
        set(h(4,o),'Ylim',[0 0.2])
        set(h(4,o),'Xlim',[x(1) x(end)])
    end
    imagesc(x, [], beh.speedAve{g,2},'parent',h(3-(g-1),4));
    hold(h(3-(g-1),4),'on');
    plot(h(3-(g-1),4), [x(1) x(end)],[recday1-nlastexp recday1-nlastexp],'r');
    plot(h(3-(g-1),4), [x(1) x(end)],[recday1 recday1],'r');
    beh.speedAve{g,outcomeval(2)}(beh.speedAve{g,outcomeval(2)}==0) = NaN;
    hold(h(4,4),'on');
    plot(h(4,4), x, mean(beh.speedAve{g,outcomeval(2)}(recday1-nlastexp:recday1,:),'omitnan'),color{g});
    set(h(3-(g-1),4),'Clim',[0 60])
    set(h(4,4),'Ylim',[0 100])
    set(h(4,4),'Xlim',[x(1) x(end)])
end

hold(h(1,5),'on');b = bar(h(1,5),beh.correct',color{g});
for g = 1:3
    b(g).FaceColor = color{g};b(g).EdgeColor = color{g};
end
hold(h(1,5),'on');scatter(h(1,5),find(beh.recday), zeros(1,sum(beh.recday)) ,'g','linewidth',5);
hold(h(1,5),'on');plot(h(1,5),sum(beh.nbtrials.*beh.correct)./sum(beh.nbtrials),'k');
set(h(1,5),'Ylim',[0 1])

hold(h(2,5),'on');b = bar(h(2,5),beh.delatgain' * 40,color{g});
for g = 1:3
    b(g).FaceColor = color{g};b(g).EdgeColor = color{g};
end
hold(h(2,5),'on');scatter(h(2,5),find(beh.recday), zeros(1,sum(beh.recday)) ,'g','linewidth',5);
hold(h(2,5),'on');plot(h(2,5),beh.ntrialchange,'c');
hold(h(2,5),'on');plot(h(2,5),beh.maxbadlicks,'y');
set(h(2,5),'Ylim',[0 40])

hold(h(3,5),'on');scatter(h(3,5),find(beh.recday), zeros(1,sum(beh.recday)) ,'g','linewidth',5);
hold(h(3,5),'on');plot(h(3,5),sum(beh.nbtrials),'c');
hold(h(3,5),'on');plot(h(3,5),beh.badtrials,'r');
set(h(3,5),'Ylim',[0 300])

end