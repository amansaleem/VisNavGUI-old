function PlotBehavSpks(dataplot,idx,varX,varY,varZ,Fdispmat, color)
if Fdispmat
    nbbins = 100;
    if sum(idx) > 0
        map = full(sparse(varY(idx), floor(varX(idx))+1, varZ(idx), max(varY(idx)), nbbins));
        dataplot.palette = 'RedWhiteBlue';
        Clim = [-max(map(:)) max(map(:))];
        dataplot.PlotMatrix(1:nbbins, 1:max(varY(idx)), flipud(map), [], true, 'Clim',Clim);
    end
else
    dataplot.PlotVector(varX(idx), varY(idx), [], [], true, 'LineStyle', 'none', 'Marker', '.', 'MarkerEdgeColor', color, 'Xlim', [min(varX) max(varX)]);
end
end