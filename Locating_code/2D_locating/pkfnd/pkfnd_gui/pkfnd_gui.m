function varargout = pkfnd_gui(varargin)
% PKFND_GUI MATLAB code for pkfnd_gui.fig
%      This is a gui that helps you optimise the pkfnd settings
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pkfnd_gui

% Last Modified by GUIDE v2.5 25-Nov-2016 12:22:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pkfnd_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pkfnd_gui_OutputFcn, ...
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


% --- Executes just before pkfnd_gui is made visible.
function pkfnd_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pkfnd_gui (see VARARGIN)

handles.im=imread(varargin{1});
%handles.im2=handles.im;
axes(handles.image_for_pkfnd)
imagesc(handles.im)

numSteps1 = 255;
set(handles.slider1, 'Min', 0);
set(handles.slider1, 'Max', numSteps1-1);
set(handles.slider1, 'Value', 0);
set(handles.slider1, 'SliderStep', [1/(numSteps1-1) , 1/(numSteps1-1) ]);
% save the current/last slider value
handles.lastSlider1Val = get(handles.slider1,'Value');

numSteps2 = 100;
set(handles.slider2, 'Min', 0);
set(handles.slider2, 'Max', numSteps2-1);
set(handles.slider2, 'Value', 0);
set(handles.slider2, 'SliderStep', [1/(numSteps2-1) , 1/(numSteps2-1) ]);
% save the current/last slider value
handles.lastSlider2Val = get(handles.slider2,'Value');

% Choose default command line output for pkfnd_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pkfnd_gui wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pkfnd_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = [get(handles.slider1,'Value'),get(handles.slider2,'Value')];
delete(hObject)

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles = guidata(hObject);
sliderValue = get(handles.slider1,'Value');
set(handles.slider1,'Value',round(sliderValue));
set(handles.text13,'String',num2str(round(sliderValue)));
handles.coords=pkfnd(handles.im,get(handles.slider1,'Value'),get(handles.slider2,'Value'));
axes(handles.image_for_pkfnd)
imagesc(handles.im), hold on
scatter(handles.coords(:,1),handles.coords(:,2),'xr')
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles = guidata(hObject);
sliderValue = get(handles.slider2,'Value');
set(handles.slider2,'Value',round(sliderValue));
set(handles.text14,'String',num2str(round(sliderValue)));
handles.coords=pkfnd(handles.im,get(handles.slider1,'Value'),get(handles.slider2,'Value'));
axes(handles.image_for_pkfnd)
imagesc(handles.im), hold on
scatter(handles.coords(:,1),handles.coords(:,2),'xr')
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%close(gcf)
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1);
%delete(hObject);
