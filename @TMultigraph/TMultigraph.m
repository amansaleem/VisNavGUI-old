classdef TMultigraph < handle

properties
    %parent figure
    parentfig
    %
    parentgrid
    %uix.TabPanel object
    tabpanel
    %uix.CardPanel object
    page
    %uix.HBoxFlex or VBoxFlex object
    pageLayout
    %uix.GridFlex object
    window
    %
    subwindow
    %
    permwindow
end

 methods
     function obj = TMultigraph(Name, title, figurehandle) 
         if nargin < 3
             figurehandle = figure('Units','normalized','Position',[0 0 1 1],'Name',Name,'NumberTitle','off');
         end
         if nargin < 2
             title = 'Page 1';
         end
         obj.parentfig = figurehandle;
         obj.parentgrid = uix.GridFlex('Parent',obj.parentfig, 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1] );
         obj.tabpanel = uix.TabPanel('Parent', obj.parentgrid, 'Padding', 5, 'TabWidth', 100, 'BackgroundColor', [0.8 0.8 0.8] );
         pagecount = 1;
         obj.page = uix.CardPanel('Parent',obj.tabpanel);
         obj.pageLayout{pagecount} = uix.HBoxFlex('Parent',obj.page(pagecount), 'Spacing', 5, 'Padding', 5 , 'BackgroundColor', [1 1 1] );
         obj.window{1,1} = uix.GridFlex('Parent',obj.pageLayout{pagecount}, 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1] );
         obj.tabpanel.TabTitles(1) = {title};                  
     end
     
     function obj = addPage(obj, title)
         pagecount = numel(obj.page) + 1;
         if nargin < 2
             title = ['Page ' num2str(obj.pagecount)];
         end
         obj.page(pagecount) = uix.CardPanel('Parent',obj.tabpanel);
         obj.pageLayout{pagecount} = uix.HBoxFlex('Parent',obj.page(pagecount), 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1] );
         obj.window{pagecount,1} = uix.GridFlex('Parent',obj.pageLayout{pagecount}, 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1] );
         obj.tabpanel.TabTitles(pagecount) = {title};
     end
     
     function obj = deletePage(obj, pagnum)
         delete(obj.page(pagnum));
     end
     
     function obj = addPermWindow(obj,str,dim)
         permwincount = numel(obj.permwindow) + 1;
         switch str
             case {'Horizontal','H'}
                 obj.permwindow{permwincount} = uix.HBoxFlex('Parent',obj.parentgrid, 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1] );
                 obj.parentgrid.Contents = [obj.parentgrid.Contents(end) obj.parentgrid.Contents(1:end-1)];
                 set(obj.parentgrid, 'Heights', dim, 'Spacing', 5 );
             case {'Vertical','V'}
                 obj.permwindow{permwincount} = uix.VBoxFlex('Parent',obj.parentgrid, 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1] );
                 obj.parentgrid.Contents = [obj.parentgrid.Contents(end) obj.parentgrid.Contents(1:end-1)];
                 set(obj.parentgrid, 'Widths', dim, 'Spacing', 5 );
         end
     end
     
     function obj = deletePermWindow(obj, winnum)
         delete(obj.permwindow(winnum));
     end
     
     function obj = RenamePage(obj, pagnum, title)
         obj.tabpanel.TabTitles(pagnum) = {title};
     end
     
     function obj = Hdividepage(obj, pagenum, nwin, dim)
         delete(obj.page(pagenum).Children);
         obj.pageLayout{pagenum} = uix.HBoxFlex('Parent',obj.page(pagenum), 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1]);
         for w = 1:nwin
             obj.window{pagenum,w} = uix.GridFlex('Parent',obj.pageLayout{pagenum}, 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1]);
         end
         set( obj.pageLayout{pagenum}, 'Widths', dim, 'Spacing', 5 );
     end 
     
     function obj = Vdividepage(obj, pagenum, nwin, dim)
         delete(obj.page(pagenum).Children);
         obj.pageLayout{pagenum} = uix.VBoxFlex('Parent',obj.page(pagenum), 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1]);
         for w = 1:nwin
             obj.window{pagenum,w} = uix.GridFlex('Parent',obj.pageLayout{pagenum}, 'Spacing', 5, 'Padding', 5, 'BackgroundColor', [1 1 1]);
         end
         set( obj.pageLayout{pagenum}, 'Heights', dim, 'Spacing', 5 , 'BackgroundColor', [1 1 1]);
     end 
     
     function obj = SetWindowType(obj,pagenum, winnum, str)
         if isprop(obj.pageLayout{pagenum},'Heights')
            oldHeights = get(obj.pageLayout{pagenum},'Heights');
         end
         if isprop(obj.pageLayout{pagenum},'Widths')
            oldWidths = get(obj.pageLayout{pagenum},'Widths');
         end
         delete(obj.window{pagenum,winnum});
         switch str
             case {'Horizontal','H'}
                 obj.window{pagenum,winnum} = uix.HBoxFlex('Parent',obj.pageLayout{pagenum}, 'Spacing', 3, 'Padding', 3, 'BackgroundColor', [1 1 1]);                 
             case {'Vertical','V'}
                 obj.window{pagenum,winnum} = uix.VBoxFlex('Parent',obj.pageLayout{pagenum}, 'Spacing', 3, 'Padding', 3, 'BackgroundColor', [1 1 1]);
             case {'Grid','G'}
                 obj.window{pagenum,winnum} = uix.GridFlex('Parent',obj.pageLayout{pagenum}, 'Spacing', 3, 'Padding', 3, 'BackgroundColor', [1 1 1]);
         end         
         ncontent = numel(get(obj.pageLayout{pagenum},'Contents'));
         if winnum < ncontent
             contentorder = 1:ncontent;
             contentorder(winnum) = ncontent;
             contentorder(winnum+1:end) = winnum:(ncontent - 1);
             obj.pageLayout{pagenum}.Contents = obj.pageLayout{pagenum}.Contents(contentorder);
         end
         if isprop(obj.pageLayout{pagenum},'Heights')
            set( obj.pageLayout{pagenum}, 'Heights', oldHeights, 'BackgroundColor', [1 1 1]);
         end
         if isprop(obj.pageLayout{pagenum},'Widths')
            set( obj.pageLayout{pagenum}, 'Widths', oldWidths, 'BackgroundColor', [1 1 1]);
         end
     end
     
     function obj = SetWindowLayout(obj,pagenum, winnum, dim, str)
         switch str
             case {'Vertical','V'}
                 if pagenum > 0
                    set( obj.window{pagenum,winnum}, 'Heights', dim, 'BackgroundColor', [1 1 1]);
                 else
                    set( obj.permwindow{winnum}, 'Heights', dim, 'BackgroundColor', [1 1 1]);
                 end
             case {'Horizontal','H'}
                 if pagenum > 0
                     set( obj.window{pagenum,winnum}, 'Widths', dim, 'BackgroundColor', [1 1 1]);
                 else
                     set( obj.permwindow{winnum}, 'Widths', dim, 'BackgroundColor', [1 1 1]);
                 end
         end
     end
     
     function obj = DivideWindow(obj,pagenum, winnum, nlines, ncolumns)
         delete(obj.window{pagenum,winnum}.Children);
         for j = 1:ncolumns
             for i = 1:nlines             
                 obj.subwindow{pagenum,winnum,i,j} = uipanel('Parent',obj.window{pagenum,winnum});
             end
         end
         obj.SetWindowLayout(pagenum, winnum, -1*ones(1, nlines), 'V');
         obj.SetWindowLayout(pagenum, winnum, -1*ones(1, ncolumns), 'H');
     end
     
     function obj = fig2pdf(obj,filename,Fappend,figname)
         if nargin<2
             [FileName,PathName] = uiputfile;
             if FileName~=0
                 filename = [PathName FileName];
             end
         end
         if nargin<3
             Fappend = false;
         end
         if nargin<4
             figname = '';
         end
         
         dot=strfind(filename,'.');
         if isempty(dot)
             filename = [filename '.pdf'];
         end
         pagenum = obj.tabpanel.Selection;
         
         Scrsize = get(0,'ScreenSize');
         if ~isempty(figname)
             fig2save = figure('Visible','on','Position',[Scrsize(1) Scrsize(2) Scrsize(3) Scrsize(4)],'name',figname,'renderer','painters');
         else
             fig2save = figure('Visible','on','Position',[Scrsize(1) Scrsize(2) Scrsize(3) Scrsize(4)],'renderer','painters');
         end
         child = get(obj.page(pagenum),'Children');
         child = get(get(child,'Children'),'Children');
         for c = 1:numel(child)
             if ~isempty(child{c})
                 if ~isequal(child{c}.Type,'uicontrol')
                     copyobj(child{c},fig2save);
                 end
             end
         end
         set(fig2save,'Units','Inches');
         pos = get(fig2save,'Position');
         set(fig2save,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
         print(fig2save,filename,'-dpdf','-r0')
         
%             if Fappend
%                 export_fig(fig2save,filename, '-painters','-bookmark','-append');
%             else
%                 export_fig(fig2save,filename, '-painters','-bookmark');
% %                 export_fig(fig2save,filename);
%             end
            delete(fig2save);
     end
 end

end
