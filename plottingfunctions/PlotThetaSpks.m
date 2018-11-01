function PlotThetaSpks(dataplot, idx, varX, Phs)
colormap_val = colormap(hsv);
dataplot.PlotVector(varX(idx), Phs(idx), [], true, 'k.', 'MarkerSize', 5);
dataplot.SetAxesParams('XLim', [min(varX(idx)) max(varX(idx))], 'YLim', [0 360*2]);
end