function PlotVar = UpdateplotBehavior(PlotVar,EXP,Layout)
PlotVar.Plots = [];
es = EXP.data.es;
page = 4;
win = 1;
Layout.DivideWindow(page, win,  2, 2);
nplot = 0;

roomlength = 100;
ctraj = es.traj;

ctraj_rewcentered = es.traj;
idx1 = (ctraj_rewcentered < floor(roomlength/2));
idx2 = (ctraj_rewcentered >= floor(roomlength/2));
ctraj_rewcentered(idx1) = ctraj_rewcentered(idx1) + floor(roomlength/2);
ctraj_rewcentered(idx2) = ctraj_rewcentered(idx2) - floor(roomlength/2);

nplotlick = 0;
nplotspeedX = 0;
nploteyeX = 0;
nplottraj = 0;
colorcond = cell(numel(EXP.SubsetVal.gain),numel(EXP.SubsetVal.contrast));
for g = 1:numel(EXP.SubsetVal.gain)
    for c = 1:numel(EXP.SubsetVal.contrast)
        rgb = zeros(1,3);
        rgb(g) = 1 * c/numel(EXP.SubsetVal.contrast);
        colorcond{g,c} = rgb;
    end
end
cidx = find(ismember(EXP.SubsetVal.contrast,[0.5 0.7]));
try
for c = 1:numel(cidx)
    colorcond{EXP.SubsetVal.gain == 0.4,cidx(c)} = 'c';
    colorcond{EXP.SubsetVal.gain == 0.5,cidx(c)} = [0.5 0.5 0.5];
    colorcond{EXP.SubsetVal.gain == 0.6,cidx(c)} = 'm';
end
catch
end

for c = 1:size(PlotVar.ChosenContrast,2)
    for g = 1:size(PlotVar.ChosenGain,2)
        for r = 1:size(PlotVar.ChosenRoomlength,2)
            for o = 1:size(PlotVar.ChosenOutcome,2)
                for k = 1:size(PlotVar.ChosenObj,1)
                    tidx = EXP.getSubsets(PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o),PlotVar.speed_th, true, true);
                    exptidx = 0;
                    seriesid = unique(es.series);
                    for sid = 1:numel(seriesid)
                        exptidx = exptidx + (es.series == seriesid(sid) & ismember(es.iexp,PlotVar.explist(sid,:)));
                    end
                    tidx = tidx & exptidx;
                    if strcmp(PlotVar.ChosenObj{k}, 'Traj x Time')
                        if nplottraj == 0
                            nplot = nplot + 1;
                            nplottraj = nplot;
                            PlotVar.Plots{nplottraj} = TDataPlot(Layout.subwindow{page, win, 1, 1});
                            PlotVar.Plots{nplottraj}.PlotVector((es.sampleTimes')./60, es.traj', [], [], true, 'Color', 'k');
                            if isfield(es,'lick')
                                PlotVar.Plots{nplottraj}.PlotVector((es.sampleTimes(es.lick == 1 & ismember(es.outcome,2)))./60, es.traj(es.lick == 1 & ismember(es.outcome,2)), [], [], true, 'Marker', 'o', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g','LineStyle', 'none');
                                PlotVar.Plots{nplottraj}.PlotVector((es.sampleTimes(es.lick == 1 & ismember(es.outcome,[0 1 3 4 5])))./60, es.traj(es.lick == 1 & ismember(es.outcome,[0 1 3 4 5])), [], [], true, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r','LineStyle', 'none');
                            end       
                        end
                    end
                    if strcmp(PlotVar.ChosenObj{k}, 'run Speed x X')
                        if nplotspeedX == 0
                            nplot = nplot + 1;
                            nplotspeedX = nplot;
                            PlotVar.Plots{nplotspeedX} = TDataPlot(Layout.subwindow{page, win, 1, 2});
                        end
                        runspeed = NaN(size(es.ballspeed));
                        sampleRate = mean(1./es.sampleSize);
%                         runspeed(~isnan(es.ballspeed)) = smthInTime(es.ballspeed(~isnan(es.ballspeed)), sampleRate, 150, 'same', [], 'boxcar_centered');
                        runspeed = es.smthBallSpd;
                        [speedpro,x,speedprostd] = fast1Dmap(ctraj(tidx),runspeed(tidx),1,sampleRate,[], es.CircularMaze);
                        speedpro = interp1([x x(end) + x(2) - x(1)],[speedpro' speedpro(1)],0:0.5:x(end)-0.5,'spline');
                        speedprostd = interp1([x x(end) + x(2) - x(1)],[speedprostd' speedprostd(1)],0:0.5:x(end)-0.5,'spline');
                        speedpro = speedpro/60;
                        speedprostd = speedprostd/60;
                        x = 0:0.5:x(end)-0.5;
                        X = x;%x - floor(roomlength/2);
                        PlotVar.Plots{nplotspeedX}.PlotVector(X, speedpro, speedprostd, [], true,'Xlim', [0 roomlength], 'Ylim', [0 1.5*max(speedpro)],  'linewidth', 2,...
                                                              'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)},...
                                                              'FaceColor', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                    end
                    if strcmp(PlotVar.ChosenObj{k}, 'visual Speed x X')
                        if nplotspeedX == 0
                            nplot = nplot + 1;
                            nplotspeedX = nplot;
                            PlotVar.Plots{nplotspeedX} = TDataPlot(Layout.subwindow{page, win, 1, 2});
                        end
                        sampleRate = mean(1./es.sampleSize);
                        visualspeed = NaN(size(es.trajspeed));
                        visualspeed(~isnan(es.trajspeed)) = smthInTime(es.trajspeed(~isnan(es.trajspeed))/0.5, sampleRate, 150, 'same', [], 'boxcar_centered');
                        [speedpro,x,speedprostd] = fast1Dmap(ctraj(tidx),visualspeed(tidx),1,sampleRate,[], es.CircularMaze);%[map,x,y] = fast2Dmap(X,Y,Z,dx,dy,scaling,nbXsmthbins,Fcircular)
                        speedpro = interp1([x x(end) + x(2) - x(1)],[speedpro' speedpro(1)],0:0.5:x(end)-0.5,'spline');
                        speedprostd = interp1([x x(end) + x(2) - x(1)],[speedprostd' speedprostd(1)],0:0.5:x(end)-0.5,'spline');
                        speedpro = speedpro/60;
                        speedprostd = speedprostd/60;
                        x = 0:0.5:x(end)-0.5;
                        X = x;%x - floor(roomlength/2);
                        PlotVar.Plots{nplotspeedX}.PlotVector(X, speedpro, speedprostd, [], true, 'Xlim', [0 roomlength], 'Ylim', [0 1.2*max(speedpro)],  'linewidth', 2,...
                                                              'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)},...
                                                              'FaceColor', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                    end
                    if strcmp(PlotVar.ChosenObj{k}, 'eye Xpos')
                        if nploteyeX == 0
                            nplot = nplot + 1;
                            nploteyeX = nplot;
                            PlotVar.Plots{nploteyeX} = TDataPlot(Layout.subwindow{page, win, 2, 2});
                        end
                        sampleRate = mean(1./es.sampleSize);
                        eyeX = NaN(size(es.eyeXpos));
                        eyeX(~isnan(es.eyeXpos)) = smthInTime(es.eyeXpos(~isnan(es.eyeXpos)), sampleRate, 150, 'same', [], 'boxcar_centered');
                        [eyepro,x,eyeprostd] = fast1Dmap(ctraj(tidx & ~isnan(eyeX)),eyeX(tidx & ~isnan(eyeX)),1,sampleRate,[],es.CircularMaze);
%                         eyepro = interp1([x x(end) + x(2) - x(1)],[eyepro' eyepro(1)],0:0.5:x(end)-0.5,'spline');
%                         eyeprostd = interp1([x x(end) + x(2) - x(1)],[eyeprostd' eyeprostd(1)],0:0.5:x(end)-0.5,'spline');
                        eyepro = eyepro/60;
                        eyeprostd = eyeprostd/60;
%                         x = 0:0.5:x(end)-0.5;
                        X = x;%x - floor(roomlength/2);
                        PlotVar.Plots{nploteyeX}.PlotVector(X, eyepro, eyeprostd, [], true, 'Xlim', [0 roomlength], 'Ylim', [1.2*min(eyepro-eyeprostd) 1.2*max(eyepro+eyeprostd)],  'linewidth', 2,...
                                                              'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)},...
                                                              'FaceColor', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                    end
                    if strcmp(PlotVar.ChosenObj{k}, 'eye Ypos')
                        if nploteyeX == 0
                            nplot = nplot + 1;
                            nploteyeX = nplot;
                            PlotVar.Plots{nploteyeX} = TDataPlot(Layout.subwindow{page, win, 2, 2});
                        end
                        sampleRate = mean(1./es.sampleSize);
                        [eyepro,x,eyeprostd] = fast1Dmap(ctraj(tidx),es.eyeYpos(tidx),1,sampleRate,[], es.CircularMaze);
                        eyepro = interp1([x x(end) + x(2) - x(1)],[eyepro' eyepro(1)],0:0.5:x(end)-0.5,'spline');
                        eyeprostd = interp1([x x(end) + x(2) - x(1)],[eyeprostd' eyeprostd(1)],0:0.5:x(end)-0.5,'spline');
                        eyepro = eyepro/60;
                        eyeprostd = eyeprostd/60;
                        x = 0:0.5:x(end)-0.5;
                        X = x;%x - floor(roomlength/2);
                        PlotVar.Plots{nploteyeX}.PlotVector(X, eyepro, eyeprostd, [], true, 'Xlim', [0 roomlength], 'Ylim', [1.2*min(eyepro-eyeprostd) 1.2*max(eyepro+eyeprostd)],  'linewidth', 2,...
                                                              'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)},...
                                                              'FaceColor', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                    end
                    if strcmp(PlotVar.ChosenObj{k}, 'pupil Size')
                        if nploteyeX == 0
                            nplot = nplot + 1;
                            nploteyeX = nplot;
                            PlotVar.Plots{nploteyeX} = TDataPlot(Layout.subwindow{page, win, 2, 2});
                        end
                        sampleRate = mean(1./es.sampleSize);
                        [eyepro,x,eyeprostd] = fast1Dmap(ctraj(tidx),es.pupilSize(tidx),1,sampleRate,[], es.CircularMaze);
                        eyepro = interp1([x x(end) + x(2) - x(1)],[eyepro' eyepro(1)],0:0.5:x(end)-0.5,'spline');
                        eyeprostd = interp1([x x(end) + x(2) - x(1)],[eyeprostd' eyeprostd(1)],0:0.5:x(end)-0.5,'spline');
                        eyepro = eyepro/60;
                        eyeprostd = eyeprostd/60;
                        x = 0:0.5:x(end)-0.5;
                        X = x;%x - floor(roomlength/2);
                        PlotVar.Plots{nploteyeX}.PlotVector(X, eyepro, eyeprostd, [], true, 'Xlim', [0 roomlength], 'Ylim', [min(eyepro-eyeprostd) max(eyepro+eyeprostd)],  'linewidth', 2,...
                                                              'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)},...
                                                              'FaceColor', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                    end
                    if strcmp(PlotVar.ChosenObj{k}, 'Licks x X')
                        if nplotlick == 0
                            nplot = nplot + 1;
                            nplotlick = nplot;
                            PlotVar.Plots{nplotlick} = TDataPlot(Layout.subwindow{page, win, 2, 1});
                        end
                        sampleRate = mean(1./es.sampleSize);
                        licks = smthInTime(double(es.firstgoodlick | es.firstbadlick), sampleRate, 0,[],[],'boxcarsum');
                        
%                         perf = sum(es.firstgoodlick(tidx))/(sum(es.firstgoodlick(tidx)) + sum(es.firstbadlick(tidx)))
                        
                        if sum(tidx) > 0
                           [f,x] = ecdf(ctraj(licks > 0 & tidx));
                           [n,~] = ecdfhist(f,x,0:roomlength-1);
                           lickdistri = n;
                           
                           [lickdistri,x] = fast1Dmap(ctraj_rewcentered(tidx),licks(tidx),2,60,[], es.CircularMaze);%
                           lickdistri = lickdistri./sum(lickdistri)/(1/numel(x));
                           lickdistri = interp1([x x(end) + x(2) - x(1)],[lickdistri' lickdistri(1)],0:0.5:(x(end) + x(2) - x(1))-0.5,'spline');
                           
%                            lickdistri = smooth(lickdistri,5)/sum(smooth(lickdistri,5)); 
%                            lickdistri = lickdistri./max(lickdistri);
%                            lickdistri = f;                           
                        end
%                         lickdistri = [0 lickdistri' 1];
%                         x = [0 x' 100];
                        x = 0:0.5:(x(end) + x(2) - x(1))-0.5;%0:roomlength-1;
                        X = x - floor(roomlength/2);
                        PlotVar.Plots{nplotlick}.PlotVector(X, lickdistri, [], [], true, 'Xlim', [X(1) X(end)], 'linewidth', 2, 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                    end
                    if strcmp(PlotVar.ChosenObj{k}, 'Cum. Licks x X')
                        if nplotlick == 0
                            nplot = nplot + 1;
                            nplotlick = nplot;
                            PlotVar.Plots{nplotlick} = TDataPlot(Layout.subwindow{page, win, 2, 1});
                        end
                        sampleRate = mean(1./es.sampleSize);
                        licks = smthInTime(double(es.firstgoodlick | es.firstbadlick), sampleRate, 0,[],[],'boxcarsum');
                        
%                         perf = sum(es.firstgoodlick(tidx))/(sum(es.firstgoodlick(tidx)) + sum(es.firstbadlick(tidx)))
                        
                        if sum(tidx) > 0
                           [f,x] = ecdf(ctraj(licks > 0 & tidx));
                           [n,~] = ecdfhist(f,x,0:roomlength-1);
                           lickdistri = n;
                           
                           [lickdistri,x] = fast1Dmap(ctraj_rewcentered(tidx),licks(tidx),2,60,[], es.CircularMaze);
                           lickdistri = interp1([x x(end) + x(2) - x(1)],[lickdistri' lickdistri(1)],0:0.5:x(end)-0.5,'spline');
                           lickdistri = lickdistri./sum(lickdistri);
                           lickdistri = cumsum(lickdistri);
%                            lickdistri = smooth(lickdistri,5)/sum(smooth(lickdistri,5)); 
%                            lickdistri = lickdistri./max(lickdistri);
%                            lickdistri = f;                           
                        end
%                         lickdistri = [0 lickdistri' 1];
%                         x = [0 x' 100];
                        x = 0:0.5:x(end)-0.5;%0:roomlength-1;
                        X = x - floor(roomlength/2);
                        PlotVar.Plots{nplotlick}.PlotVector(X, lickdistri, [], [], true, 'Xlim', [X(1) X(end)], 'linewidth', 2, 'Color', colorcond{PlotVar.ChosenGain(:,g),PlotVar.ChosenContrast(:,c)});
                    end
                end
            end
        end
    end
end
end