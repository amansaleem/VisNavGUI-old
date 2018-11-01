classdef Tdialog < handle
properties
    %number of columns which must divide the parent panel of the dialog
    %Default is 1.
    nbcolumn
    %RGB MATLAB vector defining the background color of the dialog.
    background
    %
    backgroundTxt
end

properties (Hidden=true,Dependent=true)
    %Size of a column of the dialog parent panel
    %all units are expressed in normalized units of the parent panel
    %size
    columnsize
    %
    nbYuicontrol
end

properties (Hidden=true,SetAccess=protected)        
    %handle of the dialog parent panel
    parentpanel
    %width of the dialog parent panel
    parentwidth
    %height of the dialog parent panel
    parentheight
    %vector of handles containing all handles of the dialog uicontrol
    uicontrolhandle
    %cell array of names of the dialog uicontrol
    uicontrolname
    %vertical spacing interval (normalized units)
    yinterfactor
    %horizontal spacing interval (normalized units)
    xinterfactor
    %height of one given uicontrol (normalized units)
    linesize
    %figure handle created when using shownonmodal and showmodal methods
    modalfigure
end

methods
    function obj=Tdialog(parenthandle)        
        obj.parentpanel = parenthandle;
        psize = get(obj.parentpanel,'Position');
        obj.parentwidth = psize(3);
        obj.parentheight = psize(4);
        obj.xinterfactor = 1/60;
        obj.yinterfactor = 1/500;
        obj.linesize = 1/45;
        obj.nbcolumn = 1;
        obj.background = [0.9 0.9 0.9];
        obj.backgroundTxt = [1 1 1];
    end

    function columnsize=get.columnsize(obj)
        if obj.nbcolumn ~= 0
            columnsize = 1/obj.nbcolumn;
        else
            error('bad Tdialog.nbcolumn property value');
        end
    end

    function obj=set.background(obj,RGB)
        obj.background = RGB;
        set(obj.parentpanel,'BackgroundColor',obj.background);
    end

    function nbYuicontrol = get.nbYuicontrol(obj)
        nbYuicontrol = numel(obj.uicontrolhandle);
    end

    function obj=getvariable(obj,str,f_callback,varval)
        htxtgetvar = uicontrol('Parent',obj.parentpanel,'Style','text',...
                                 'String',str,'BackgroundColor',obj.backgroundTxt);

        if nargin<4
            varval=0;
        end
        hgetvar = uicontrol('Parent',obj.parentpanel,'Style','edit',...
                            'String',num2str(varval),'BackgroundColor',obj.background,...
                            'Callback',f_callback);

        if isempty(obj.uicontrolhandle)
            obj.uicontrolhandle = {hgetvar};
            obj.uicontrolname = {str};
        else
            obj.uicontrolhandle = [obj.uicontrolhandle {hgetvar}];
            obj.uicontrolname = [obj.uicontrolname {str}];
        end

    end

    function obj=addslider(obj,min1,max1,dx1,dx2,slideval)
        lasthandle=length(obj.uicontrolhandle);
        lastcallback=get(obj.uicontrolhandle{lasthandle},'Callback');
        if nargin<6
            slideval=min1;
        end
        haddslider=uicontrol('Parent',obj.parentpanel,'Style','slider','Min',min1,'Max',max1,'Value',slideval,'SliderStep',[dx1/(max1-min1) dx2/(max1-min1)],'BackgroundColor',obj.background,...
                             'Callback',{@slider_callback});

        function slider_callback(source,event)
            slidevalue=get(source,'Value');
            slidevalue=round(slidevalue/dx1)*dx1;
            set(source,'Value',slidevalue);
            set(obj.uicontrolhandle{lasthandle},'String',slidevalue);
            lastcallback(obj.uicontrolhandle{lasthandle},event);
        end

        if isempty(obj.uicontrolhandle)
            obj.uicontrolhandle = {haddslider};
        else
            obj.uicontrolhandle = [obj.uicontrolhandle {haddslider}];
            obj.uicontrolname = [obj.uicontrolname {[obj.uicontrolname{lasthandle} '_slider']}];
        end            
    end

    function obj=addgetpathbutton(obj,str)
        lasthandle=length(obj.uicontrolhandle);
        lastcallback=get(obj.uicontrolhandle{lasthandle},'Callback');
        haddpathbutton = uicontrol('Parent',obj.parentpanel,'Style','pushbutton',...
                                 'String',str,'BackgroundColor',obj.background-[0.1 0.1 0.1],...
                                 'Callback',{@path_callback});

        function path_callback(source,event)
            Dirpath=uigetdir(get(obj.uicontrolhandle{lasthandle},'String'));
            if ischar(Dirpath)
                set(obj.uicontrolhandle{lasthandle},'String',Dirpath);
                lastcallback(obj.uicontrolhandle{lasthandle},event);
            end
        end
        if isempty(obj.uicontrolhandle)
            obj.uicontrolhandle={haddpathbutton};
            obj.uicontrolname = {str};
        else
            obj.uicontrolhandle=[obj.uicontrolhandle {haddpathbutton}];
            obj.uicontrolname = [obj.uicontrolname {str}];
        end

    end

    function obj=addgetfilebutton(obj,str)
        lasthandle=length(obj.uicontrolhandle);
        lastcallback=get(obj.uicontrolhandle{lasthandle},'Callback');
        haddfilebutton = uicontrol('Parent',obj.parentpanel,'Style','pushbutton',...
                                 'String',str,'BackgroundColor',obj.background-[0.1 0.1 0.1],...
                                 'Callback',{@path_callback});

        function path_callback(source,event)
            [FileName,PathName,~]=uigetfile(get(obj.uicontrolhandle{lasthandle},'String'));
            if ischar(FileName)
                set(obj.uicontrolhandle{lasthandle},'String',[PathName FileName]);
                lastcallback(obj.uicontrolhandle{lasthandle},event);
            end
        end
        if isempty(obj.uicontrolhandle)
            obj.uicontrolhandle={haddfilebutton};
            obj.uicontrolname = {str};
        else
            obj.uicontrolhandle=[obj.uicontrolhandle {haddfilebutton}];
            obj.uicontrolname = [obj.uicontrolname {str}];
        end

    end


    function obj=getpushbutton(obj,str,f_callback)
        hgetpushbutton = uicontrol('Parent',obj.parentpanel,'Style','pushbutton',...
                                 'String',str,'BackgroundColor',obj.background-[0.1 0.1 0.1],...
                                 'Callback',f_callback);            

        if isempty(obj.uicontrolhandle)
            obj.uicontrolhandle={hgetpushbutton};
            obj.uicontrolname = {str};
        else
            obj.uicontrolhandle=[obj.uicontrolhandle {hgetpushbutton}];
            obj.uicontrolname = [obj.uicontrolname {str}];
        end
    end

    function obj=getpopupmenu(obj,title,strval,choice0,f_callback)
        choice1=0;
        for k=1:size(strval,2)
            if strcmp(choice0,strval{k})
                choice1=k;
            end
        end
        if choice1==0;
            choice1=1;
        end

        htxtpopupmenu = uicontrol('Parent',obj.parentpanel,'Style','text',...
                                 'String',title,'BackgroundColor',obj.backgroundTxt);

        hpopupmenu=uicontrol('Parent',obj.parentpanel,'Style','popupmenu',...
                            'String',strval,'Value',choice1,'BackgroundColor',obj.background,...
                            'Callback',f_callback);

        if isempty(obj.uicontrolhandle)
            obj.uicontrolhandle={hpopupmenu};
            obj.uicontrolname = {title};
        else
            obj.uicontrolhandle=[obj.uicontrolhandle {hpopupmenu}];
            obj.uicontrolname = [obj.uicontrolname {title}];
        end   

    end

    function obj=getlistbox(obj,title,strval,selmin,selmax,f_callback)
       htxtlistbox = uicontrol('Parent',obj.parentpanel,'Style','text',...
                                 'String',title,'BackgroundColor',obj.backgroundTxt);

       hlistbox=uicontrol('Parent',obj.parentpanel,'Style','listbox','min',selmin,'max',selmax,...
                            'String',strval,'BackgroundColor',obj.background,...
                            'Callback',f_callback);

        if isempty(obj.uicontrolhandle)
            obj.uicontrolhandle={hlistbox};
            obj.uicontrolname = {title};
        else
            obj.uicontrolhandle=[obj.uicontrolhandle {hlistbox}];
            obj.uicontrolname = [obj.uicontrolname {title}];
        end       
    end

    function obj=getcheckbox(obj,title,f_callback,checkval)        
        if nargin<4
            checkval=0;
        end

        hgetcheckbox = uicontrol('Parent',obj.parentpanel,'Style','checkbox',...
                                 'String',title,'Min',0,'Max',1,'Value',checkval,'BackgroundColor',obj.background,...
                                 'Callback',f_callback);            

        if isempty(obj.uicontrolhandle)
            obj.uicontrolhandle = {hgetcheckbox};
            obj.uicontrolname = {title};
        else
            obj.uicontrolhandle = [obj.uicontrolhandle {hgetcheckbox}];
            obj.uicontrolname = [obj.uicontrolname {title}];
        end

    end


    function obj=settext(obj,title,str)
        htxtsettext = uicontrol('Parent',obj.parentpanel,'Style','text',...
                            'String',title);
        
        hsettext=uicontrol('Parent',obj.parentpanel,'Style','listbox','Enable','off',...
                            'String',str,'BackgroundColor',obj.background);

        if isempty(obj.uicontrolhandle)
            obj.uicontrolhandle = {hsettext};
            obj.uicontrolname = {title};
        else
            obj.uicontrolhandle = [obj.uicontrolhandle {hsettext}];
            obj.uicontrolname = [obj.uicontrolname {title}];
        end
    end  
    
    function obj=UpdateUIcontrol(obj,line,varargin)
        if ischar(line)
            num = 0;
            Fmatch = false;
            while ~Fmatch
                num = num + 1;
                Fmatch = isequal(obj.uicontrolname{num},line);                
            end
            set(obj.uicontrolhandle{num},varargin{:});
        else
            set(obj.uicontrolhandle{line},varargin{:});
        end        
    end
end
    
methods (Hidden = true)
    function obj=updatepanelHeight(obj)
        last=length(obj.uicontrolhandle);
        set(obj.parentpanel,'Units','pixels');
        for k=1:last
            set(obj.uicontrolhandle(k),'Units','pixels');
        end
        poslast=get(obj.uicontrolhandle(last),'Position');
        bottom=poslast(2);
        for k=1:length(obj.uicontrolhandle)
            pos=get(obj.uicontrolhandle(k),'Position');
            pos(2)=pos(2)-bottom;
            set(obj.uicontrolhandle(k),'Position',pos);
        end

        pos=get(obj.parentpanel,'Position');
        pos(4)=pos(4)-bottom;
        set(obj.parentpanel,'Position',pos);

        set(obj.parentpanel,'Units','normalized');
        for k=1:length(obj.uicontrolhandle)
            set(obj.uicontrolhandle(k),'Units','normalized')
        end
    end
end    
end