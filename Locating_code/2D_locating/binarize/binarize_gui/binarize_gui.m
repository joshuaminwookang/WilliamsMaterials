function varargout = binarize_gui(varargin)
% 
% NAME:
%               binarize_gui
% PURPOSE:
%               Provides a visual way to optimise the parameters when
%               you're using binarize on an image. Creates a gui with the
%               image that you would like to bbinarize. Then has three slider
%               bars for the three bbinarize options (threshold,erosion,dilation)
%               that you can play around with. Finally, when you close the
%               figure, it returns the optimal parameters.
% 
% CATEGORY:
%               Image Processing
% CALLING SEQUENCE:
%               b=binarize_gui(image_name)
% INPUTS:
%               image_name:  The file name of the image to be filtered.
%
% OUTPUTS:
%               b:    A vector containing the three values of threshold,
%               erosion and dilation to be used with binarize.
% NOTES:
% As an example, there is a TFM image in this folder.
% binarize can be used to find the particles in this figure by appropriate choices
% of the parameters. To run the example, use
% b=binarize_gui('image1.tif')
% 
% POSSIBLE IMPROVEMENTS NEEDED:
% 1) Modify so that if image_name is an image array (rather than the string
% containing the image's filename), the code still works.
%
% MODIFICATION HISTORY:
%               Written by Hendrik Spanke, ETH, Dec 2016
%               Slider values edited by Katrina Smith-Mannschott, ETH, Apr
%               2018

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @binarize_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @binarize_gui_OutputFcn, ...
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


% --- Executes just before binarize_gui is made visible.
function binarize_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to binarize_gui (see VARARGIN)

if isa(varargin{1},'char')==1
    handles.im=imread(varargin{1});
else
    handles.im=varargin{1};
end
handles.im2=handles.im;
axes(handles.image_for_binarize)
imagesc(handles.im)

numSteps1 = max(handles.im(:));
set(handles.slider1, 'Min', 0);
set(handles.slider1, 'Max', numSteps1-1);
set(handles.slider1, 'Value', 0);
set(handles.slider1, 'SliderStep', [1/(numSteps1-1) , 1/(numSteps1-1) ]);
set(handles.text9,'String',num2str(numSteps1-1));
% save the current/last slider value
handles.lastSlider1Val = get(handles.slider1,'Value');

numSteps2 = 20;
set(handles.slider2, 'Min', 0);
set(handles.slider2, 'Max', numSteps2-1);
set(handles.slider2, 'Value', 0);
set(handles.slider2, 'SliderStep', [1/(numSteps2-1) , 1/(numSteps2-1) ]);
% save the current/last slider value
handles.lastSlider2Val = get(handles.slider2,'Value');

numSteps3 = 20;
set(handles.slider3, 'Min', 0);
set(handles.slider3, 'Max', numSteps3-1);
set(handles.slider3, 'Value', 0);
set(handles.slider3, 'SliderStep', [1/(numSteps3-1) , 1/(numSteps3-1) ]);
% save the current/last slider value
handles.lastSlider3Val = get(handles.slider3,'Value');

% Choose default command line output for binarize_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes binarize_gui wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = binarize_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = [get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value')];
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
handles.im2=binarize(handles.im,get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value'));
axes(handles.image_for_binarize)
imagesc(handles.im2)
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
handles.im2=binarize(handles.im,get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value'));
axes(handles.image_for_binarize)
imagesc(handles.im2)
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


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles = guidata(hObject);
sliderValue = get(handles.slider3,'Value');
set(handles.slider3,'Value',round(sliderValue));
set(handles.text15,'String',num2str(round(sliderValue)));
handles.im2=binarize(handles.im,get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value'));
axes(handles.image_for_binarize)
imagesc(handles.im2)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
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


% --- Executes during object creation, after setting all properties.
function text5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
