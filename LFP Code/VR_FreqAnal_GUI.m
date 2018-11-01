function varargout = VR_FreqAnal_GUI(varargin)
% VR_FREQANAL_GUI MATLAB code for VR_FreqAnal_GUI.fig
%      VR_FREQANAL_GUI, by itself, creates a new VR_FREQANAL_GUI or raises the existing
%      singleton*.
%
%      H = VR_FREQANAL_GUI returns the handle to a new VR_FREQANAL_GUI or the handle to
%      the existing singleton*.
%
%      VR_FREQANAL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VR_FREQANAL_GUI.M with the given input arguments.
%
%      VR_FREQANAL_GUI('Property','Value',...) creates a new VR_FREQANAL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VR_FreqAnal_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VR_FreqAnal_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VR_FreqAnal_GUI

% Last Modified by GUIDE v2.5 26-Apr-2013 11:43:37

% Begin initialization code - DO NOT EDIT

addpath('\\zserver\Code\Spikes\')
addpath('\\zuser\aman\work\Code\General\')
addpath('\\zuser\aman\work\Code\General\ChronuxToolbox\spectral_analysis\continuous\')
addpath('\\zuser\aman\work\Code\Behaviour analysis\')
addpath('\\zuser\aman\work\Code\LFP analysis\')

SetDefaultDirs

global info

% info.range_low = 1;
% info.range_high= 90;
% info.downsample = 1000;
% info.twoChn = 0;

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @VR_FreqAnal_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @VR_FreqAnal_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before VR_FreqAnal_GUI is made visible.
function VR_FreqAnal_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VR_FreqAnal_GUI (see VARARGIN)

% Choose default command line output for VR_FreqAnal_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global info

info.range_low = 52;
info.range_high= 90;
info.downsample = 1000;
info.twoChn = 0;

% UIWAIT makes VR_FreqAnal_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VR_FreqAnal_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global info
info.animal = get(hObject,'String');


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
global info
info.iseries = str2num(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function handles = edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
global es;
global info
info.iexp = str2num(get(hObject,'String'));

info.specialB = 0;

% [~,~,es] = VRWheelLoad(info.animal, info.iseries, info.iexp);
set(handles.radiobutton5, 'Value', 0);
set(handles.radiobutton4, 'Value', 0);
set(handles.radiobutton2, 'Value', 0);
set(handles.radiobutton6, 'Value', 0);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
global info
global pepNEV
global ChnA
global timestamps

list_entries = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');
info.ChnA = list_entries(index_selected);
info.ChnA = str2num(info.ChnA{1});

oldpointer = get(handles.figure1, 'pointer');
set(handles.figure1, 'pointer', 'watch')
drawnow;
% set(handles.figure1, 'pointer', oldpointer)
timestamps.fulllength = length(pepNEV.ns.Data.data(info.ChnA,:));
timestamps.samplingRateOrig = info.SamplingRateInKHZ * 1000;

ChnA = decimate(double(pepNEV.ns.Data.data(info.ChnA,:)), timestamps.samplingRateOrig/info.downsample);
timestamps.ChnAsampInt = 1./info.downsample;
timestamps.ChnAorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
timestamps.ChnA = 0:timestamps.ChnAsampInt:...
    (timestamps.ChnAsampInt*(timestamps.fulllength*(info.downsample/timestamps.samplingRateOrig)));

set(handles.figure1, 'pointer', oldpointer)
if strcmp(get(handles.figure1,'pointer'),'watch')
    set(handles.figure1, 'pointer', 'arrow');
end
% disp('Done loading ChannelA');

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
global info
global pepNEV
global ChnB
global ChnB_down
global timestamps

list_entries = get(handles.listbox2,'String');
index_selected = get(handles.listbox2,'Value');
info.ChnB = list_entries(index_selected);
info.ChnB = str2num(info.ChnB{1});

oldpointer = get(handles.figure1, 'pointer');
set(handles.figure1, 'pointer', 'watch')
drawnow;
if ~info.specialB
    % set(handles.figure1, 'pointer', oldpointer)
    timestamps.fulllength = length(pepNEV.ns.Data.data(info.ChnB,:));
    timestamps.samplingRateOrig = info.SamplingRateInKHZ * 1000;

    ChnB = decimate(double(pepNEV.ns.Data.data(info.ChnB,:)), timestamps.samplingRateOrig/info.downsample);
    timestamps.ChnBsampInt = 1./info.downsample;
    timestamps.ChnBorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
    timestamps.ChnB = 0:timestamps.ChnBsampInt:...
    (timestamps.ChnBsampInt*(timestamps.fulllength*(info.downsample/timestamps.samplingRateOrig)));
    if length(timestamps.ChnB)>length(ChnB)
        timestamps.ChnB = timestamps.ChnB(1:length(ChnB));
    end
    ChnB_down = [];
    timestamps.ChnB_down = [];
else
    timestamps.fulllength = length(pepNEV.ns.Data.data(info.ChnB,:));
    timestamps.samplingRateOrig = info.SamplingRateInKHZ * 1000;

    ChnB = zeros(1, size(pepNEV.ns.Data.data,2));
    maxC = max(pepNEV.ns.Data.data(info.ChnB,:));
    minC = min(pepNEV.ns.Data.data(info.ChnB,:));
    ChnB((abs(diff(pepNEV.ns.Data.data(info.ChnB,:)-minC)./(maxC-minC)))>0.5) = 1;
    ChnB_down = sum(reshape(ChnB(1:end-rem(size(pepNEV.ns.Data.data,2),250)),250,[]),1);
    
    timestamps.ChnBsampInt = 1./info.downsample;
    timestamps.ChnBorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
    timestamps.ChnB = 0:timestamps.ChnBsampInt:...
    (timestamps.ChnBsampInt*(timestamps.fulllength*(info.downsample/timestamps.samplingRateOrig)));
    if length(timestamps.ChnB)>length(ChnB)
        timestamps.ChnB = timestamps.ChnB(1:length(ChnB));
    end
    timestamps.ChnBsampInt_down = 1./(timestamps.samplingRateOrig./250);
    timestamps.ChnB_down = 0:timestamps.ChnBsampInt_down:...
    (timestamps.ChnBsampInt_down*(timestamps.fulllength*(1/250)));
    if length(timestamps.ChnB_down)>length(ChnB_down)
        timestamps.ChnB_down = timestamps.ChnB_down(1:length(ChnB_down));
    end
    info.runSpd = ChnB_down;
    info.smthRunSpd = smthInTime(info.runSpd, 1./timestamps.ChnBsampInt_down, info.runSmthWin);
    info.time = timestamps.ChnB_down;
end
set(handles.figure1, 'pointer', oldpointer)
if strcmp(get(handles.figure1,'pointer'),'watch')
    set(handles.figure1, 'pointer', 'arrow');
end
% disp('Done loading ChannelB');

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
global info
info.range_low = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
global info
info.range_high = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global info
global ChnA
global ChnB
global ChnB_down
% clear global pepNEV

if ~isnan(info.ChnA) & ~isnan(info.ChnB)
    info.twoChn = 1;
else
    info.twoChn = 0;
end

if info.twoChn
    [f_A,pow_A] = powerSpectrum(ChnA,info.downsample,0,'b');
    if info.specialB
        [f_B,pow_B] = powerSpectrum(ChnB_down,120,0,'b');
    else
        [f_B,pow_B] = powerSpectrum(ChnB,info.downsample,0,'b');
    end
    
    range_A = f_A<info.range_high & f_A>info.range_low;
    range_B = f_B<info.range_high & f_B>info.range_low;
    
%     figure;
    plot(handles.axes1,f_A(range_A), smooth(f_A(range_A).*pow_A(range_A))./max(f_A(range_A).*pow_A(range_A)),...
        f_B(range_B), smooth(f_B(range_B).*pow_B(range_B))./max(f_B(range_B).*pow_B(range_B)));
else
    [f_B,pow_B] = powerSpectrum(ChnB,60,0,'b');
    range_B = f_B<info.range_high & f_B>info.range_low;
%     figure;
    plot(handles.axes1,f_B(range_B), smooth(f_B(range_B).*pow_B(range_B))./max(f_B(range_B).*pow_B(range_B)));
end
xlabel('Hz')
ylabel('Normalised Power')
title('Whitened power spectrum');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global info
global ChnA
global ChnB
global ChnB_down
global procSpec
global timestamps
% clear global pepNEV

oldpointer = get(handles.figure1, 'pointer');
set(handles.figure1, 'pointer', 'watch')
drawnow;

if ~isnan(info.ChnA) & ~isnan(info.ChnB)
    info.twoChn = 1;
    if ~info.specialB
        params.Fs=info.downsample;
        params.fpass=[info.range_low info.range_high];
        params.tapers=[3 1];
        movingwin=[5  0.5];
        procSpec.ChnA = ChnA';
        procSpec.ChnB = ChnB';
        
        [procSpec.C,procSpec.phi,procSpec.SAB,procSpec.SA,procSpec.SB,procSpec.t,procSpec.f]=cohgramc(procSpec.ChnA,procSpec.ChnB,movingwin,params);
        axes(handles.axes4);
        imagesc(procSpec.t,procSpec.f,abs(procSpec.SAB)'); axis xy
        axes(handles.axes5);
        imagesc(procSpec.t,procSpec.f,abs(procSpec.SB')); axis xy
        axes(handles.axes6);
        imagesc(procSpec.t,procSpec.f,abs(procSpec.SA')); axis xy
    else
        downRate = floor(length(timestamps.ChnA)/length(timestamps.ChnB_down));
        nChnA = decimate(ChnA,downRate);
        tChnA = decimate(timestamps.ChnA,downRate);
        nChnA = interp1q(tChnA', nChnA', timestamps.ChnB_down')';
        
        params.Fs=(1./timestamps.ChnBsampInt_down);
        params.fpass=[info.range_low min(info.range_high, params.Fs/2)];
        params.tapers=[5 9];
        movingwin=[8  0.8];
        procSpec.ChnA = nChnA';
        procSpec.ChnB = ChnB_down';
        [procSpec.C,procSpec.phi,procSpec.SAB,procSpec.SA,procSpec.SB,procSpec.t,procSpec.f]=cohgramc(procSpec.ChnA,procSpec.ChnB,movingwin,params);
        axes(handles.axes4);
        imagesc(procSpec.t,procSpec.f,abs(procSpec.SAB)'); axis xy
        axes(handles.axes5);
        imagesc(procSpec.t,procSpec.f,abs(procSpec.SB)'); axis xy
        axes(handles.axes6);
        imagesc(procSpec.t,procSpec.f,abs(procSpec.SA)'); axis xy
    end
else
    info.twoChn = 0;
end
set(handles.figure1, 'pointer', oldpointer)
if strcmp(get(handles.figure1,'pointer'),'watch')
    set(handles.figure1, 'pointer', 'arrow');
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global pepNEV

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global info
global ChnA
global ChnB
global ChnB_down
% clear global pepNEV

if ~isnan(info.ChnA) & ~isnan(info.ChnB)
    info.twoChn = 1;
else
    info.twoChn = 0;
end

if info.twoChn
    [f_A,pow_A] = powerSpectrum(ChnA,info.downsample,0,'b');
    if info.specialB
        [f_B,pow_B] = powerSpectrum(ChnB_down,120,0,'b');
    else
        [f_B,pow_B] = powerSpectrum(ChnB,info.downsample,0,'b');
    end
    
    range_A = f_A<info.range_high & f_A>info.range_low;
    range_B = f_B<info.range_high & f_B>info.range_low;
    
%     figure;
    plot(handles.axes2,f_A(range_A), smooth(pow_A(range_A))./max(pow_A(range_A)),...
        f_B(range_B), smooth(pow_B(range_B))./max(pow_B(range_B)));
else
    [f_B,pow_B] = powerSpectrum(ChnB,60,0,'b');
    range_B = f_B<info.range_high & f_B>info.range_low;
%     figure;
    plot(handles.axes2,f_B(range_B), smooth(pow_B(range_B))./max(pow_B(range_B)));
end
xlabel('Hz')
ylabel('Normalized Power')
title('Power spectrum');


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
global info
info.specialB = get(hObject,'Value');

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
global info
info.downsample = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
global pepNEV
global info
global ChnB

loadNS = ~get(hObject,'Value');
[~,~,es] = VRWheelLoad(info.animal,info.iseries,info.iexp);
info.runSpd = es.ballspeed(~isnan(es.ballspeed));
info.time   = es.sampleTimes(~isnan(es.ballspeed));
info.runSmthWin = 300;

if loadNS
    set(handles.radiobutton5,'Value',1);
    expName = [info.animal '_' num2str(info.iseries) '_' num2str(info.iexp) '.ns5'];
    ns5file = ['\\ZSERVER\Data\Cerebus\' info.animal filesep num2str(info.iseries) filesep num2str(info.iexp) filesep expName]
    [~,info.SamplingRateInKHZ,info.nchan] = nsopen2(ns5file);
    info.ChnA = NaN;
    info.ChnB = 0;
    chan_list = [];
    for ichan = 1:info.nchan
        chan_list = [chan_list {num2str(ichan)}];
    end
    set(handles.listbox1,'String',chan_list);
    set(handles.listbox2,'String',chan_list);
    info.specialB = 0;
else
    set(handles.radiobutton5,'Value',0);
    info.ChnA = NaN;
    info.ChnB = 0;
    info.specialB = 1;
    [~,~,es] = VRWheelLoad(info.animal,info.iseries,info.iexp);
    ChnB = es.ballspeed;
end
info.smthRunSpd = smthInTime(info.runSpd,60,info.runSmthWin);
plot(handles.axes3,info.time, info.smthRunSpd);
axes(handles.axes3); axis tight
if ~isnan(info.ChnA)
    info.twoChn = 1;
else
    info.twoChn = 0;
end


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
global pepNEV
global info
global ChnB

loadNS = get(hObject,'Value');
[~,~,es] = VRWheelLoad(info.animal,info.iseries,info.iexp);
info.runSpd = es.ballspeed(~isnan(es.ballspeed));
info.time   = es.sampleTimes(~isnan(es.ballspeed));
info.runSmthWin = 300;

if loadNS
    set(handles.radiobutton4,'Value',0);
    expName = [info.animal '_' num2str(info.iseries) '_' num2str(info.iexp) '.ns5'];
    ns5file = ['\\ZSERVER\Data\Cerebus\' info.animal filesep num2str(info.iseries) filesep num2str(info.iexp) filesep expName]
    [~,info.SamplingRateInKHZ,info.nchan] = nsopen2(ns5file);
    info.ChnA = NaN;
    info.ChnB = 0;
    chan_list = [];
    for ichan = 1:info.nchan
        chan_list = [chan_list {num2str(ichan)}];
    end
    chan_list;
    set(handles.listbox1,'String',chan_list);
    set(handles.listbox2,'String',chan_list);
    info.specialB = 0;
else
    set(handles.radiobutton4,'Value',1);
    info.ChnA = NaN;
    info.ChnB = 0;
    info.specialB = 1;
    ChnB = info.runSpd;
end
info.smthRunSpd = smthInTime(info.runSpd,60,info.runSmthWin);
plot(handles.axes3,info.time, info.smthRunSpd);
axes(handles.axes3); axis tight
if ~isnan(info.ChnA)
    info.twoChn = 1;
else
    info.twoChn = 0;
end


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6

% Sort in time. 
global info
global procSpec

downRate = floor(length(info.time)/length(procSpec.t));
nRun = decimate(info.smthRunSpd,downRate);
tRun = decimate(info.time,downRate);
procSpec.runSpd = interp1q(tRun', nRun', procSpec.t')';
[~,procSpec.order] = sort(procSpec.runSpd);

plot(handles.axes3, procSpec.runSpd(procSpec.order)); 
axes(handles.axes3); axis tight;
axes(handles.axes4);
imagesc(1:length(procSpec.t),procSpec.f,abs(procSpec.SAB(procSpec.order,:))'); axis xy
axes(handles.axes5);
imagesc(1:length(procSpec.t),procSpec.f,abs(procSpec.SB(procSpec.order,:)')); axis xy
axes(handles.axes6);
imagesc(1:length(procSpec.t),procSpec.f,abs(procSpec.SA(procSpec.order,:)')); axis xy

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
global info
info.runSmthWin = str2double(get(hObject,'String'));
info.smthRunSpd = smthInTime(info.runSpd,60,info.runSmthWin);
plot(handles.axes3,info.time, info.smthRunSpd);
axes(handles.axes3); axis tight;

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global info
global ChnA
global ChnB
global ChnB_down
global procSpec
global timestamps

filename = [info.animal '_' num2str(info.iseries) '_' num2str(info.iexp) '_'];
assignin('base',[filename 'info'],info)
assignin('base',[filename 'timestamps'],timestamps)
assignin('base',[filename 'Chn' num2str(info.ChnA)],ChnA);
assignin('base',[filename 'Chn' num2str(info.ChnB)],ChnB);
assignin('base',[filename 'Chn' num2str(info.ChnB) '_down'],ChnB_down);
assignin('base',[filename 'Chn' num2str(info.ChnA) '_' 'Chn' num2str(info.ChnB) '_coh'],procSpec);