function [lickpos_actual,lickpos_dec] = PopplotLicks(res,dx,Tsmth_win,Tshift)
expt = getExperimentList;
figure;
sampleRate = 60;
delay = round(Tshift/1000*sampleRate);
iprobe = 1;
amp_th = 1;
maxtol = 1;
lickpos_actual = cell(1,3);
lickpos_dec = cell(1,3);
lickdistri_mean = cell(1,3);
lickdistridec_mean = cell(1,3);
for g = 1:3
    lickdistri_mean{g} = 0;
    lickdistridec_mean{g} = 0;
end
for ianimal = 1:size(res.lick,1)
    for iseries = 1:size(res.lick,2)
        if ~isempty(res.Xpred_ave{ianimal,iseries,iprobe})
            if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                contref = find(ismember(res.contrastVal{ianimal,iseries}, [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]));
                RLref = find(ismember(res.roomlengthVal{ianimal,iseries}, [1]));
                outcomeref = find(ismember(res.outcomeVal{ianimal,iseries}, [2]));
                gainlist = [0.4 0.5 0.6];
                tidx = cell(1,3);
                for g = 1:3
                    tidx{g} = false(size(res.tidx{ianimal,iseries,iprobe,1,2,1,3}));
                    for c = 1:numel(contref)
                        for rl = 1:numel(RLref)
                            for o = 1:numel(outcomeref)
                                if ~isempty(res.tidx{ianimal,iseries,iprobe,contref(c),g,RLref(rl),outcomeref(o)})
                                    tidx{g} = tidx{g} | res.tidx{ianimal,iseries,iprobe,contref(c),g,RLref(rl),outcomeref(o)};
                                end
                            end
                        end
                    end
                end
                lick = res.firstgoodlick{ianimal,iseries};%smthInTime(res.firstgoodlick{ianimal,iseries}, sampleRate, Tsmth_win, 'same', [], 'boxcar_centered');
%                 Xtemp = unwrap(res.X{ianimal,iseries}/max(res.X{ianimal,iseries})*2*pi)*max(res.X{ianimal,iseries})/(2*pi);
%                 Xtemp(~isnan(Xtemp)) = smthInTime(Xtemp(~isnan(Xtemp)), sampleRate, Tsmth_win, 'same', [], 'boxcar_centered');
%                 Xtemp(~isnan(Xtemp)) = mod(Xtemp(~isnan(Xtemp)),max(res.X{ianimal,iseries}));
%                 res.X{ianimal,iseries} = Xtemp;
%                 
                Xtemp = unwrap(res.Xpred_ave{ianimal,iseries,iprobe}/max(res.Xpred_ave{ianimal,iseries,iprobe})*2*pi)*max(res.Xpred_ave{ianimal,iseries,iprobe})/(2*pi);
                Xtemp(~isnan(Xtemp)) = smthInTime(Xtemp(~isnan(Xtemp)), sampleRate, Tsmth_win, 'same', [], 'boxcar_centered');
                Xtemp(~isnan(Xtemp)) = mod(Xtemp(~isnan(Xtemp)),max(res.Xpred_ave{ianimal,iseries,iprobe}));
                res.Xpred_ave{ianimal,iseries,iprobe} = Xtemp;
                
                Xshift = res.Xpred_ave{ianimal,iseries,iprobe}-res.X{ianimal,iseries};
                Xshift(Xshift>50) = Xshift(Xshift>50) - 100;
                Xshift(Xshift<-50) = Xshift(Xshift<-50) + 100;
                gooddecidx = abs(Xshift) < 50;
                
                if sum(tidx{1}) > 0 && sum(tidx{3}) > 0
                    subplot(1,2,1);
                    X = circshift(res.X{ianimal,iseries}(tidx{2} & gooddecidx),delay);
                    Y = lick(tidx{2} & gooddecidx);
                    [map,x,mapstd] = fast1Dmap(X,Y,dx,1,25,1);
                    map = circshift(map/sum(map)*100/dx,round(100/dx/2));
                    lickdistri_mean{2} = lickdistri_mean{2} + map;
                    [~,imax] = max(map);
                    istart = find(map(1:imax)<amp_th,1,'last');
                    iend =  imax -1 + find(map(imax:end)<amp_th,1,'first');
                    map([1:istart iend:end]) = 0;
                    hold off;plot(x,map);
                    lickpos = getCircularAverage(map,amp_th,maxtol);
                    lickpos = x(floor(lickpos)) + rem(lickpos,1)*(x(2) - x(1));
                    lickpos_actual{2} = [lickpos_actual{2} lickpos];
                    
                    X = circshift(res.X{ianimal,iseries}(tidx{1} & gooddecidx),delay);
                    Y = lick(tidx{1} & gooddecidx);
                    [map,x,mapstd] = fast1Dmap(X,Y,dx,1,25,1);
                    map = circshift(map/sum(map)*100/dx,round(100/dx/2));
                    lickdistri_mean{1} = lickdistri_mean{1} + map;
                    [~,imax] = max(map);
                    istart = find(map(1:imax)<amp_th,1,'last');
                    iend =  imax -1 + find(map(imax:end)<amp_th,1,'first');
                    map([1:istart iend:end]) = 0;
                    hold on;plot(x,map,'c');
                    lickpos = getCircularAverage(map,amp_th,maxtol);
                    lickpos = x(floor(lickpos)) + rem(lickpos,1)*(x(2) - x(1));
                    lickpos_actual{1} = [lickpos_actual{1} lickpos-lickpos_actual{2}(end)];
                    
                    
                    X = circshift(res.X{ianimal,iseries}(tidx{3} & gooddecidx),delay);
                    Y = lick(tidx{3} & gooddecidx);
                    [map,x,mapstd] = fast1Dmap(X,Y,dx,1,25,1);
                    map = circshift(map/sum(map)*100/dx,round(100/dx/2));
                    lickdistri_mean{3} = lickdistri_mean{3} + map;
                    [~,imax] = max(map);
                    istart = find(map(1:imax)<amp_th,1,'last');
                    iend =  imax -1 + find(map(imax:end)<amp_th,1,'first');
                    map([1:istart iend:end]) = 0;
                    hold on;plot(x,map,'m');
                    lickpos = getCircularAverage(map,amp_th,maxtol);
                    lickpos = x(floor(lickpos)) + rem(lickpos,1)*(x(2) - x(1));
                    lickpos_actual{3} = [lickpos_actual{3} lickpos-lickpos_actual{2}(end)];
                    
                    
                    set(gca,'Xlim',[1 100]);
                    
                    subplot(1,2,2);
                    try
                        X = circshift(res.Xpred_ave{ianimal,iseries,iprobe}(tidx{2} & gooddecidx),delay);
                        Y = lick(tidx{2} & gooddecidx);
                        [map,x,mapstd] = fast1Dmap(X,Y,dx,1,25,1);
                        map = circshift(map/sum(map)*100/dx,round(100/dx/2));
                        lickdistridec_mean{2} = lickdistridec_mean{2} + map;
                        [~,imax] = max(map);
                        istart = find(map(1:imax)<amp_th,1,'last');
                        iend =  imax -1 + find(map(imax:end)<amp_th,1,'first');
                        map([1:istart iend:end]) = 0;
                        hold off;plot(x,map);
                        lickpos = getCircularAverage(map,amp_th,maxtol);
                        lickpos = x(floor(lickpos)) + rem(lickpos,1)*(x(2) - x(1));
                        lickpos_dec{2} = [lickpos_dec{2} lickpos];
                    catch
                        lickpos_dec{2} = [lickpos_dec{2} NaN];
                    end
                    try
                        X = circshift(res.Xpred_ave{ianimal,iseries,iprobe}(tidx{1} & gooddecidx),delay);
                        Y = lick(tidx{1} & gooddecidx);
                        [map,x,mapstd] = fast1Dmap(X,Y,dx,1,25,1);
                        map = circshift(map/sum(map)*100/dx,round(100/dx/2));
                        lickdistridec_mean{1} = lickdistridec_mean{1} + map;
                        [~,imax] = max(map);
                        istart = find(map(1:imax)<amp_th,1,'last');
                        iend =  imax -1 + find(map(imax:end)<amp_th,1,'first');
                        map([1:istart iend:end]) = 0;
                        hold on;plot(x,map,'c');
                        lickpos = getCircularAverage(map,amp_th,maxtol);
                        lickpos = x(floor(lickpos)) + rem(lickpos,1)*(x(2) - x(1));
                        lickpos_dec{1} = [lickpos_dec{1} lickpos-lickpos_dec{2}(end)];
                    catch
                        lickpos_dec{1} = [lickpos_dec{1} NaN];
                    end
                    try
                        X = circshift(res.Xpred_ave{ianimal,iseries,iprobe}(tidx{3} & gooddecidx),delay);
                        Y = lick(tidx{3} & gooddecidx);
                        [map,x,mapstd] = fast1Dmap(X,Y,dx,1,25,1);
                        map = circshift(map/sum(map)*100/dx,round(100/dx/2));
                        lickdistridec_mean{3} = lickdistridec_mean{3} + map;
                        [~,imax] = max(map);
                        istart = find(map(1:imax)<amp_th,1,'last');
                        iend =  imax -1 + find(map(imax:end)<amp_th,1,'first');
                        map([1:istart iend:end]) = 0;
                        hold on;plot(x,map,'m');
                        lickpos = getCircularAverage(map,amp_th,maxtol);
                        lickpos = x(floor(lickpos)) + rem(lickpos,1)*(x(2) - x(1));
                        lickpos_dec{3} = [lickpos_dec{3} lickpos-lickpos_dec{2}(end)];
                    catch
                        lickpos_dec{3} = [lickpos_dec{3} NaN];
                    end
                    set(gca,'Xlim',[1 100]);
                end
                delete(gcf)
%                 pause;
            end
        end
    end
end
figure;
subplot(2,2,1);
hlow = histogram(lickpos_actual{1},-20.25:0.5:20.25,'EdgeColor','none','FaceColor','c');
hold on;hhigh = histogram(lickpos_actual{3},-20.25:0.5:20.25,'EdgeColor','none','FaceColor','m');
subplot(2,2,2);
hlow = histogram(lickpos_dec{1},-20.25:0.5:20.25,'EdgeColor','none','FaceColor','c');
hold on;hhigh = histogram(lickpos_dec{3},-20.25:0.5:20.25,'EdgeColor','none','FaceColor','m');
subplot(2,2,3);
plot(lickdistri_mean{1}/sum(lickdistri_mean{1}),'c');
hold on;plot(lickdistri_mean{2}/sum(lickdistri_mean{2}),'k');
hold on;plot(lickdistri_mean{3}/sum(lickdistri_mean{3}),'m');
subplot(2,2,4);
plot(lickdistridec_mean{1}/sum(lickdistridec_mean{1}),'c');
hold on;plot(lickdistridec_mean{2}/sum(lickdistridec_mean{2}),'k');
hold on;plot(lickdistridec_mean{3}/sum(lickdistridec_mean{3}),'m');
end