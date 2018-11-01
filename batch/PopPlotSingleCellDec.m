function PopPlotSingleCellDec(resdec,varname,Xvar,Yvar)
if nargin < 2
    varname = 'SpatialInfoPerSpike';%'SSI';
    Xvar = 'avedec';%'_traj';
    Yvar = 'avedecCross';
end
probestr = {'CA1','V1'};
gainstr = {'low gain','medium gain','high gain','All'};
outcomestr = {'error','All','Correct'};
for iprobe = 1:2
    goodcellidx = resdec.fieldpos{iprobe}>=0 & resdec.fieldpos{iprobe}<=100;% & resdec.SSI_traj{iprobe,2,3}>0.5;
    figure('Name',probestr{iprobe});
    for o = 1:3
        for g = 1:4
            subplot(3,4,(3-o)*4+g)
            scatter(resdec.SSI_traj{iprobe,g,o}(goodcellidx),resdec.SSI_avedec{iprobe,g,o}(goodcellidx));
            scatter(resdec.([varname '_' Xvar]){iprobe,g,o}(goodcellidx),resdec.([varname '_' Yvar]){iprobe,g,o}(goodcellidx));
            maxXY = max([resdec.SSI_traj{iprobe,g,o}(goodcellidx) resdec.SSI_avedec{iprobe,g,o}(goodcellidx)]);
            set(gca,'Xlim',[0 maxXY],'Ylim',[0 maxXY]);
            hold on;plot([0 maxXY],[0 maxXY]);
            if o == 3
                title(gainstr{g});
            end
            if o == 1
                xlabel([varname ' - actual']);
            end
            if g == 1
                ylabel([varname ' - decoded - ' outcomestr{o}]);
            end
        end
    end
end
end