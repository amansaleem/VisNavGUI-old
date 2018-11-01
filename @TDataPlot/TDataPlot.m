classdef TDataPlot < matlab.mixin.SetGet

properties
    %axes
    ax
    %Line or surface object
    Tplot
    %title
    title
    %
    palette
    %
    Colorbar
    %
    Fhold
end

methods
    function obj = TDataPlot(winhandle, ijsub, str)
        if nargin < 2
            ijsub = [1 1];
        end
        if nargin < 3
            str = '';
        end
        ksub = 0;
        for i = 1:ijsub(1)
            for j = 1:ijsub(2)
                ksub = ksub+ 1;
                obj.ax{ksub} = subplot(ijsub(1),ijsub(2),ksub,'Parent',winhandle);
            end
        end
        
        obj.Tplot = [];
        obj.title = str;
        obj.palette = 'RedWhiteBlue';%'default';%
        
        if ~isempty(str)
            title(str);
        end
        
        obj.addlistener('ObjectBeingDestroyed',@deleteaxs);
        
        function deleteaxs(~,~)
            if ishandle(obj.ax)==1
                delete(obj.ax);
                delete(obj);
            end
        end
    end
    
    function obj = set.palette(obj,palettename)
        obj.palette = palettename;
        if ~strcmp(obj.palette,'default')
            fcolor = str2func(obj.palette);
            for ksub = 1:numel(obj.ax)
                colormap(obj.ax{ksub},fcolor());
            end
        end
    end
    
    function obj = set.Colorbar(obj, onoffstr)
        obj.Colorbar = onoffstr;
        for ksub = 1:numel(obj.ax)
            colorbar(obj.ax{ksub},obj.Colorbar);
        end
    end
    
    function obj = PlotVector(obj, Xdata, Ydata, Errordata, ksub, Fhold, varargin)
        if nargin < 4
            Errordata = [];
        end
        if nargin < 5
            ksub = 1;
        end
        if isempty(ksub)
            ksub = 1;
        end
        if nargin < 6
            Fhold = false;
        end
        if isempty(Xdata)
            Xdata = 1:numel(Ydata);
        end
        Ydata = Ydata(:);
        Xdata = Xdata(:);
        if ~isempty(Xdata) && ~isempty(Ydata)
            if ~Fhold
                delete(obj.ax{ksub}.Children);
                obj.Tplot = [];
            else
                hold(obj.ax{ksub},'on');
            end
            if isempty(obj.Tplot)
                if isempty(Errordata)
                    obj.Tplot = {plot(obj.ax{ksub},Xdata,Ydata)};
                    newplot = 1;
                else
%                     obj.Tplot = {errorbar(obj.ax{ksub},Xdata,Ydata,Errordata)};
                    [hline, hpatch] = ciplot(obj.ax{ksub},Xdata,Ydata,Errordata,0.2);
                    obj.Tplot = [{hline, hpatch}];
                    newplot = 2;
                end
            else
                if isempty(Errordata)
                    obj.Tplot = [obj.Tplot {plot(obj.ax{ksub},Xdata,Ydata)}];
                    newplot = 1;
                else
%                     obj.Tplot = [obj.Tplot {errorbar(obj.ax{ksub},Xdata,Ydata,Errordata)}];
                    [hline, hpatch] = ciplot(obj.ax{ksub},Xdata,Ydata,Errordata,0.2);
                    obj.Tplot = [obj.Tplot {hline, hpatch}];
                    newplot = 2;
                end
            end
            
            axis tight;
            for p = numel(obj.Tplot) - newplot + 1:numel(obj.Tplot)
                [varplot, varax] = obj.SortProperties(p, varargin);
                obj.SetPlotParams(p,varplot{:});
            end
            obj.SetAxesParams(varax{:});
        end
    end
    
    function obj = PlotMatrix(obj,Xdata,Ydata,Zdata, ksub,Fhold,varargin)
        if nargin < 5
            ksub = 1;
        end
        if isempty(ksub)
            ksub = 1;
        end
        if nargin < 6
            Fhold = false;
        end
        if isempty(Xdata)
            Xdata = 1:size(Zdata,2);
        end
        if isempty(Ydata)
            Ydata = 1:size(Zdata,1);
        end
        if ~isempty(Xdata) && ~isempty(Ydata) && ~isempty(Zdata)
            if ~Fhold
                delete(obj.ax{ksub}.Children);
                obj.Tplot = [];
            else
                hold(obj.ax{ksub},'on');
            end
            
            if isempty(obj.Tplot)
                obj.Tplot = {imagesc(Xdata,Ydata,Zdata,'Parent',obj.ax{ksub})};
            else
                obj.Tplot = [obj.Tplot {imagesc(Xdata,Ydata,Zdata,'Parent',obj.ax{ksub})}];
            end
            
            az=0;el=90;
            view(obj.ax{ksub},az,el);
            axis tight;
            
            [varplot, varax] = obj.SortProperties(numel(obj.Tplot), varargin);
            obj.SetPlotParams(numel(obj.Tplot), varplot{:});
            obj.SetAxesParams(varax{:});                        
            
            if ~strcmp(obj.palette,'default')
                fcolor = str2func(obj.palette);
                colormap(obj.ax{ksub},fcolor());
            end
        end
    end
    
    function obj = SetPlotParams(obj, p, varargin)
        if ~isempty(varargin)
            set(obj.Tplot{p}, varargin{:});
        end
    end
    
    function obj = SetAxesParams(obj, varargin)
        if ~isempty(varargin)
            for ksub = 1:numel(obj.ax)
                set(obj.ax{ksub}, varargin{:});
            end
        end
    end
    
    function [varplotout, varasxout] = SortProperties(obj, p, varin)
        propidx = false(1,size(varin,2));
        for v = 1:size(varin,2)
            if ischar(varin{v})
                propidx(v) = isprop(obj.Tplot{p},varin{v});
            end
        end
        propidx(find(propidx) + 1) = true;
        varplotout = varin(propidx);
        
        varin = varin(~propidx);
        propidx = false(1,size(varin,2));
        for v = 1:size(varin,2)
            if ischar(varin{v})
                if ~strcmpi(varin{v},'color')
                    propidx(v) = isprop(obj.ax{1}, varin{v});
                end
            end
        end
        propidx(find(propidx) + 1) = true;
        varasxout = varin(propidx);
    end
end

end