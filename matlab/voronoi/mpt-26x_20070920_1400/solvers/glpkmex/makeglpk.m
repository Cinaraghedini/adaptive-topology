function varargout = makeglpk(varargin)
% MAKEGLPK M-file for makeglpk.fig
%      MAKEGLPK, by itself, creates a new MAKEGLPK or raises the existing
%      singleton*.
%
%      H = MAKEGLPK returns the handle to a new MAKEGLPK or the handle to
%      the existing singleton*.
%
%      MAKEGLPK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAKEGLPK.M with the given input arguments.
%
%      MAKEGLPK('Property','Value',...) creates a new MAKEGLPK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before makeglpk_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to makeglpk_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help makeglpk

% Last Modified by GUIDE v2.5 28-Nov-2004 18:51:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @makeglpk_OpeningFcn, ...
                   'gui_OutputFcn',  @makeglpk_OutputFcn, ...
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


% --- Executes just before makeglpk is made visible.
function makeglpk_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to makeglpk (see VARARGIN)

current=pwd;
handles.GLPKpath=current;
if isunix 
    handles.GLPKlib=[current '/src/libglpk.a'];
    handle.GLPKinclude=[current '/include'];
else
    handles.GLPKlib=[current '\src\libglpk.a'];
    handles.GLPKinclude=[current '\include'];
end

handles.filename='glpkmex.c';
handles.cygwin=1;
if isunix,
    handles.mexopts=[current '/mexopts.bat'];
else
    handles.mexopts=[current '\mexopts.bat'];
end

set(handles.glpkpathtext,'String',handles.GLPKpath);
set(handles.libglpkatext,'String',handles.GLPKlib);

% Choose default command line output for makeglpk
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes makeglpk wait for user response (see UIRESUME)
% uiwait(handles.makeglpkfig);


% --- Outputs from this function are returned to the command line.
function varargout = makeglpk_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function glpkpathtext_Callback(hObject, eventdata, handles)
% hObject    handle to glpkpathtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of glpkpathtext as text
%        str2double(get(hObject,'String')) returns contents of glpkpathtext as a double


% --- Executes during object creation, after setting all properties.
function glpkpathtext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to glpkpathtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% glpkpathtext browse callback
%pbrowsecb = ['f=get(gcbo,''UserData'');' ...
%	     'p=uigetpath(get(f,''UserData''));' ...
%	     'if ~isempty(p),set(f,''String'',p),end'];

%eval(pbrowsecb);
pbrowsecb=uigetdir('','Select GLPK Folder');
path=pbrowsecb;

if path ~= 0
    handles.GLPKpath=path;
    if isunix,
        handles.GLPKlib=[handles.GLPKpath '/src/libglpk.a'];
        handles.GLPKinclude=[handles.GLPKpath '/include'];
    else
        handles.GLPKlib=[handles.GLPKpath '\src\libglpk.a'];
        handles.GLPKinclude=[handles.GLPKpath '\include'];
    end

    set(handles.glpkpathtext,'String',handles.GLPKpath);
    set(handles.libglpkatext,'String',handles.GLPKlib);
end

% Update handles structure
guidata(hObject, handles);


function libglpkatext_Callback(hObject, eventdata, handles)
% hObject    handle to libglpkatext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of libglpkatext as text
%        str2double(get(hObject,'String')) returns contents of libglpkatext as a double


% --- Executes during object creation, after setting all properties.
function libglpkatext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to libglpkatext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% file browse callback
[filename, pathname, filt]=uigetfile('*.*','Select the library...');

if filename ~= 0
    handles.GLPKlib=[pathname filename];
    set(handles.libglpkatext,'String',handles.GLPKlib);
end
     
guidata(hObject, handles);





% --- Executes on button press in mexopts.
function mexopts_Callback(hObject, eventdata, handles)
% hObject    handle to mexopts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% file browse callback
[filename, pathname, filt]=uigetfile('*.bat','Select the mexopts.bat...');

handles.mexopts=[pathname filename];
set(handles.libglpkatext,'String',handles.GLPKlib);
     
guidata(hObject, handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=findobj('Tag','makeglpkfig');

close(h);

% --- Executes on button press in cygwin.
function cygwin_Callback(hObject, eventdata, handles)
% hObject    handle to cygwin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cygwin

handles.cygwin=1;

guidata(hObject, handles);



% --- Executes on button press in defcompiler.
function defcompiler_Callback(hObject, eventdata, handles)
% hObject    handle to defcompiler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of defcompiler

handles.cygwin=0;

guidata(hObject, handles);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cmd=['-I' handles.GLPKinclude ' ' handles.filename ' ' handles.GLPKlib];

cygwin=[];
if handles.cygwin
    cygwin=['-f ' handles.mexopts ' '];
end
eval(['mex ' cygwin cmd]);

msgbox('GLPKMEX has been compiled successfully');


