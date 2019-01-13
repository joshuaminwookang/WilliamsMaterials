function varargout = master_tracking_gui(varargin)
% 
% NAME:
%               master_tracking_gui
% PURPOSE:
%               Provides a visual way to quickly find the best way to track
%               spherical objects in an image. First always gives you an
%               option to bpass the image (e.g. to get rid of background of
%               of pixel scale noise). After this step, offers three
%               options:
%               1) pkfnd and cntrd. pkfnd finds a rough estimate of the
%               particle position. Then cntrd refines the fit based on the
%               positions from pkfnd. NB cntrd relies on there being a
%               maximum/minimum(?) in the centre of each particle, so may
%               not work well for particles with rings of alternating
%               dark/light.
%               2) pkfnd and radialcenter. pkfnd finds a rough estimate of
%               the particle position. Then radialcenter (from Raghu)
%               refines the position by finding the center of radial
%               symmetry in the vicinity of the rough point estimates. This
%               will be better for particles that have rings of colours.
%               3) regionprops/weightedcentroid. The image is binarized
%               with matlab's binarize code.. this finds bright/dark spots
%               in the image. Then regionprops is used to find the weighted
%               centroids of these spots.
% 
% CATEGORY:
%               Image Processing
% CALLING SEQUENCE:
%               a=master_tracking_gui(image_name)
% INPUTS:
%               image_name:  The file name of the image to be filtered. Or
%               a matrix containing the image.
%
% OUTPUTS:
%               a: A structure containing the details of the chosen
%               particle detecting algorithms, and the positions of the
%               particles that have been found (a.pks). e.g.:
%               a = 
%
%                 bpass_params: [0 11 0]
%               roughfind_method: 'pkfnd'
%               roughfind_params: [36 13]
%                  refine_method: 'cntrd'
%                  refine_params: 13
%                            pks: [2496x4 double]
%
% NOTES:
% As an example, there is a figure of a sessile droplet in this folder.
% bpass can be used to find the edges in this figure by appropriate choices
% of the parameters. To run the example, use
% a=bpass_gui('image1.tif')
%
% My intention is to make a master code that will take the structure from
% above and apply it to a sequence of images for TFM/particle tracking.
% 
% POSSIBLE IMPROVEMENTS NEEDED:
% 1) I think we can do a better job with the regionprops code. At the
% moment, it just finds light blobs and then uses weightedcentroid to find
% the centre of the particles. Should we instead convert the blobs into
% circular areas centred around the centroid of the blobs and then run
% weightedcentroid on these? Perhaps needs thinking about the best way to
% avoid bias
%
% MODIFICATION HISTORY:
%               Written by Rob Style, ETH, Dec 2016


% Edit the above text to modify the response to help master_tracking_gui

% Last Modified by GUIDE v2.5 28-Dec-2016 11:19:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @master_tracking_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @master_tracking_gui_OutputFcn, ...
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


% --- Executes just before master_tracking_gui is made visible.
function master_tracking_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to master_tracking_gui (see VARARGIN)

if isa(varargin{1},'char')==1 % Work out if the image input is a matrix or a filename
    handles.im=imread(varargin{1});
else
    handles.im=varargin{1};
end
handles.im2=handles.im; % Save an unadulterated copy of the image
axes(handles.axes1)
imagesc(handles.im) % Display the image
set( gcf, 'toolbar', 'none' ) % Make a toolbar so that you can zoom in on the figure if needed
set( gcf, 'toolbar', 'figure' )

% Make the sliders for the bpass
numSteps1 = 32; 
set(handles.slider1, 'Min', 0);
set(handles.slider1, 'Max', numSteps1-1);
set(handles.slider1, 'Value', 0);
set(handles.slider1, 'SliderStep', [1/(numSteps1-1) , 1/(numSteps1-1) ]);
% save the current/last slider value
handles.lastSlider1Val = get(handles.slider1,'Value');

numSteps2 = 102;
set(handles.slider2, 'Min', 0);
set(handles.slider2, 'Max', numSteps2-1);
set(handles.slider2, 'Value', 0);
set(handles.slider2, 'SliderStep', [1/(numSteps2-1) , 1/(numSteps2-1) ]);
% save the current/last slider value
handles.lastSlider2Val = get(handles.slider2,'Value');

numSteps3 = 16;
set(handles.slider3, 'Min', 0);
set(handles.slider3, 'Max', numSteps3-1);
set(handles.slider3, 'Value', 0);
set(handles.slider3, 'SliderStep', [1/(numSteps3-1) , 1/(numSteps3-1) ]);
% save the current/last slider value
handles.lastSlider3Val = get(handles.slider3,'Value');

% Make invisible the things that we don't need (e.g. buttons etc)
set(handles.text3,'Visible','Off')
set(handles.popupmenu1,'Visible','Off')
set(handles.text7,'Visible','Off')
set(handles.edit1,'Visible','Off')
set(handles.pushbutton2,'Visible','Off')
set(handles.pushbutton1,'Visible','Off')

handles.stepnumber=1; % This is the marker that we're in the first (bpass step)

% Choose default command line output for master_tracking_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes master_tracking_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

uiwait(handles.figure1); % When you do uiresume(handles.figure1) later, it will return the results and close the figure

% --- Outputs from this function are returned to the command line.
function varargout = master_tracking_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

% Make a structure with all the information that we need
tracking_struct=struct;
if handles.stepnumber==1 % if you close before finishing the bpass step you just get the bpass parameters
    tracking_struct.bpass_params=[get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value')];
else
    tracking_struct.bpass_params=handles.bpass_params;
    if get(handles.popupmenu1,'Value')==1
        tracking_struct.roughfind_method='pkfnd';
        tracking_struct.roughfind_params=[get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value')];
        tracking_struct.refine_method='cntrd';
        tracking_struct.refine_params=[1+2*get(handles.slider3,'Value')];
        tracking_struct.pks=handles.pks;
    elseif get(handles.popupmenu1,'Value')==2
        tracking_struct.roughfind_method='pkfnd';
        tracking_struct.roughfind_params=[get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value')];
        tracking_struct.refine_method='radialcenter';
        tracking_struct.refine_params=[1+2*get(handles.slider3,'Value')];
        tracking_struct.pks=handles.pks;
    else
        regions = bwlabel(handles.im2,8);
        props = regionprops(regions,handles.im,'WeightedCentroid','Area'); % find features in BW (Particles)

        X=[]; Y=[]; Area=[];
        for n = 1:1:size(props)    % go through features in BW
            % Particles:
            X = [X; props(n).WeightedCentroid(1)];
            Y = [Y; props(n).WeightedCentroid(2)];
            Area = [Area; props(n).Area];   
        end

        handles.pks=cat(2,X,Y,Area);
        tracking_struct.roughfind_method='binarize';
        tracking_struct.roughfind_params=[get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value')];
        tracking_struct.refine_method='regionprops_weightedcentroid';
        tracking_struct.pks=handles.pks;
    end
end

varargout{1} = tracking_struct;
delete(hObject)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

handles = guidata(hObject);
popupvalue = get(handles.popupmenu1,'Value');
if popupvalue ==1 % pkfnd and cntrd chosen
    % Stuff here for doing pkfnd and cntrd particle finding
    set(handles.text2,'String','Choose optimal parameters for pkfnd and cntrd');
    
    set(handles.slider1,'Visible','On')
    set(handles.slider2,'Visible','On')
    set(handles.slider3,'Visible','On')
    
    set(handles.text4,'Visible','On')
    set(handles.text4,'String','Threshold pixel intensity')
    set(handles.text5,'Visible','On')
    set(handles.text5,'String','Particle diameter (for pkfnd - slightly larger than the size of the particle)')
    set(handles.text6,'Visible','On')
    set(handles.text6,'String','Particle diameter (for cntrd - if in doubt, copy above value)')
    set(handles.text9,'Visible','On')
    set(handles.text10,'Visible','On')
    set(handles.text11,'Visible','On')
    set(handles.text12,'Visible','On')
    set(handles.text13,'Visible','On')
    set(handles.text14,'Visible','On')
    set(handles.text15,'Visible','On')
    set(handles.text16,'Visible','On')
    set(handles.text17,'Visible','On')
    set(handles.pushbutton1,'Visible','On')
    set(handles.pushbutton2,'Visible','On')
    
    numSteps1 = 512; % Upper limit for the thresholding number for pkfnd.
    set(handles.slider1, 'Min', 0);
    set(handles.slider1, 'Max', numSteps1-1);
    guess_best_thresh=sum(sum(handles.im2.^2))/sum(sum(handles.im2))/3;
    set(handles.slider1, 'Value', round(guess_best_thresh));
    set(handles.text15, 'String', num2str(round(guess_best_thresh)));
    set(handles.text12, 'String', numSteps1-1);
    set(handles.slider1, 'SliderStep', [1/(numSteps1-1) , 1/(numSteps1-1) ]);
    % save the current/last slider value
    handles.lastSlider1Val = get(handles.slider1,'Value');

    numSteps2 = 21; % Particle size for pkfnd. THIS MUST BE ODD
    set(handles.slider2, 'Min', 1);
    set(handles.slider2, 'Max', (numSteps2-1)/2);
    % I've set it up so that the slider goes 1...x...N, but then the value
    %that is displayed (and used) is 2x+1, to ensure that the number is
    %odd. This is important to understand if modifying the code.
    estimate_from_bpass=round(handles.bpass_params(2)/2)*2+1;
    set(handles.slider2, 'Value', round((estimate_from_bpass-1)/2));
    set(handles.text16, 'String', num2str(estimate_from_bpass));
    set(handles.text13, 'String', numSteps2);
    set(handles.slider2, 'SliderStep', [ 2/(numSteps2-1) , 2/(numSteps2-1) ]);
    % save the current/last slider value
    handles.lastSlider2Val = get(handles.slider2,'Value');

    numSteps3 = 21; % Particle size for cntrd. THIS MUST BE ODD
    set(handles.slider3, 'Min', 1);
    set(handles.slider3, 'Max', (numSteps3-1)/2);
    % I've set it up so that the slider goes 1...x...N, but then the value
    %that is displayed (and used) is 2x+1, to ensure that the number is
    %odd. This is important to understand if modifying the code.
    set(handles.slider3, 'Value', round((estimate_from_bpass-1)/2));
    set(handles.text17, 'String', num2str(estimate_from_bpass));
    set(handles.text14, 'String', numSteps3);
    set(handles.slider3, 'SliderStep', [ 2/(numSteps2-1) , 2/(numSteps2-1) ]);
    % save the current/last slider value
    handles.lastSlider3Val = get(handles.slider3,'Value');
    
     handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'));
     handles.pks=cntrd(double(handles.im),handles.coords,1+2*round(get(handles.slider3,'Value')));
     axes(handles.axes1)
     imagesc(handles.im2)
     hold on
     scatter(handles.pks(:,1),handles.pks(:,2),'xr')
     hold off
    
elseif popupvalue ==2 % pkfnd and radialcenter
    % Stuff here for doing pkfnd and radial symmetry finding
    set(handles.text2,'String','Choose optimal parameters for pkfnd and symmetry finding');
    
    set(handles.slider1,'Visible','On')
    set(handles.slider2,'Visible','On')
    set(handles.slider3,'Visible','On')
    
    set(handles.text4,'Visible','On')
    set(handles.text4,'String','Threshold pixel intensity')
    set(handles.text5,'Visible','On')
    set(handles.text5,'String','Particle diameter (for pkfnd - slightly larger than the size of the particle)')
    set(handles.text6,'Visible','On')
    set(handles.text6,'String','Particle diameter (for radialcenter - if in doubt, copy above value)')
    set(handles.text9,'Visible','On')
    set(handles.text10,'Visible','On')
    set(handles.text11,'Visible','On')
    set(handles.text12,'Visible','On')
    %set(handles.text12,'String',num2str(255))
    set(handles.text13,'Visible','On')
    %set(handles.text13,'String',num2str(21))
    set(handles.text14,'Visible','On')
    %set(handles.text14,'String',num2str(20))
    set(handles.text15,'Visible','On')
    set(handles.text16,'Visible','On')
    set(handles.text17,'Visible','On')
    set(handles.pushbutton1,'Visible','On')
    set(handles.pushbutton2,'Visible','On')
    
    numSteps1 = 512;
    set(handles.slider1, 'Min', 0);
    set(handles.slider1, 'Max', numSteps1-1);
    guess_best_thresh=sum(sum(handles.im2.^2))/sum(sum(handles.im2))/3;
    set(handles.slider1, 'Value', round(guess_best_thresh));
    set(handles.text15, 'String', num2str(round(guess_best_thresh)));
    set(handles.text12, 'String', numSteps1-1);
    set(handles.slider1, 'SliderStep', [1/(numSteps1-1) , 1/(numSteps1-1) ]);
    % save the current/last slider value
    handles.lastSlider1Val = get(handles.slider1,'Value');

    numSteps2 = 400; % particle size for pkfnd. MUST BE ODD
    set(handles.slider2, 'Min', 1);
    set(handles.slider2, 'Max', (numSteps2-1)/2);
    estimate_from_bpass=round(handles.bpass_params(2)/2)*2+1;
    % I've set it up so that the slider goes 1...x...N, but then the value
    %that is displayed (and used) is 2x+1, to ensure that the number is
    %odd. This is important to understand if modifying the code.
    set(handles.slider2, 'Value', round((estimate_from_bpass-1)/2));
    set(handles.text16, 'String', num2str(estimate_from_bpass));
    set(handles.text13, 'String', numSteps2);
    set(handles.slider2, 'SliderStep', [ 2/(numSteps2-1) , 2/(numSteps2-1) ]);
    % save the current/last slider value
    handles.lastSlider2Val = get(handles.slider2,'Value');

    numSteps3 = 400; % particle size fo radialcenter. MUST BE ODD.
    set(handles.slider3, 'Min', 1);
    set(handles.slider3, 'Max', (numSteps3-1)/2);
    set(handles.slider3, 'Value', round((estimate_from_bpass-1)/2));
    % I've set it up so that the slider goes 1...x...N, but then the value
    %that is displayed (and used) is 2x+1, to ensure that the number is
    %odd. This is important to understand if modifying the code.
    set(handles.text17, 'String', num2str(estimate_from_bpass));
    set(handles.text14, 'String', numSteps3);
    set(handles.slider3, 'SliderStep', [ 2/(numSteps2-1) , 2/(numSteps2-1) ]);
    % save the current/last slider value
    handles.lastSlider3Val = get(handles.slider3,'Value');
    
    handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'));
    p=handles.coords;
    particle_sz=1+2*get(handles.slider3,'Value');
    p_rs=[];
    %zeros(length(p),3);
    for i=1:length(p)
        xcentre=round(p(i,1)); % Centre of the particle to pixel level resolution
        ycentre=round(p(i,2));
        if ycentre-(particle_sz-1)/2>0 && ycentre+(particle_sz-1)/2<=length(handles.im(:,1)) && xcentre-(particle_sz-1)/2 >0 && xcentre+(particle_sz-1)/2 <= length(handles.im(1,:))
            % Here we take an area of the image that is particle_sz by
            % particle_sz in dimensions around the position of each
            % particle. Then we find the centre of radial symmetry in this
            % test area.
            testim=handles.im(ycentre-(particle_sz-1)/2:ycentre+(particle_sz-1)/2,xcentre-(particle_sz-1)/2:xcentre+(particle_sz-1)/2);
            [c,d,f]=radialcenter(double(testim));
            p_rs=[p_rs ; xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
            %p_rs(i,:)=[xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
        end
    end
    handles.pks=p_rs;
    axes(handles.axes1)
    imagesc(handles.im2)
    hold on
    scatter(handles.pks(:,1),handles.pks(:,2),'xr')
    hold off
    
else
    % Stuff here for doing regionprops
    set(handles.text2,'String','Choose optimal parameters for regionprops thresholding');
    set(handles.slider1,'Visible','On')
    set(handles.slider2,'Visible','On')
    set(handles.slider3,'Visible','On')
    set(handles.text4,'Visible','On')
    set(handles.text4,'String','Threshold (for binarizing)')
    set(handles.text5,'Visible','On')
    set(handles.text5,'String','Erosion (number of pixels that should be eroded from the particles)')
    set(handles.text6,'Visible','On')
    set(handles.text6,'String','Dilation (number of pixels that each particle should be dilated by)')
    set(handles.text9,'Visible','On')
    set(handles.text10,'Visible','On')
    set(handles.text11,'Visible','On')
    set(handles.text12,'Visible','On')
    set(handles.text12,'String',num2str(511))
    set(handles.text13,'Visible','On')
    set(handles.text13,'String',num2str(20))
    set(handles.text14,'Visible','On')
    set(handles.text14,'String',num2str(20))
    set(handles.text15,'Visible','On')
    set(handles.text16,'Visible','On')
    set(handles.text17,'Visible','On')
    set(handles.pushbutton1,'Visible','On')
    set(handles.pushbutton2,'Visible','On')
    
    numSteps1 = 256; % Can be either even or odd
    set(handles.slider1, 'Min', 0);
    set(handles.slider1, 'Max', numSteps1-1);
    guess_best_thresh=sum(sum(handles.im2.^2))/sum(sum(handles.im2))/3;
    set(handles.slider1, 'Value', round(guess_best_thresh));
    set(handles.text15, 'String', num2str(round(guess_best_thresh)));
    set(handles.slider1, 'SliderStep', [1/(numSteps1-1) , 1/(numSteps1-1) ]);
    % save the current/last slider value
    handles.lastSlider1Val = get(handles.slider1,'Value');

    numSteps2 = 21; % Can be either even or odd
    set(handles.slider2, 'Min', 0);
    set(handles.slider2, 'Max', numSteps2-1);
    set(handles.slider2, 'Value', 0);
    set(handles.text16, 'String', num2str(0));
    set(handles.slider2, 'SliderStep', [1/(numSteps2-1) , 1/(numSteps2-1) ]);
    % save the current/last slider value
    handles.lastSlider2Val = get(handles.slider2,'Value');

    numSteps3 = 21; % Can be either even or odd
    set(handles.slider3, 'Min', 0);
    set(handles.slider3, 'Max', numSteps3-1);
    set(handles.slider3, 'Value', 0);
    set(handles.text17, 'String', num2str(0));
    set(handles.slider3, 'SliderStep', [1/(numSteps3-1) , 1/(numSteps3-1) ]);
    % save the current/last slider value
    handles.lastSlider3Val = get(handles.slider3,'Value');
    
    handles.im2=binarize(handles.im_pass,get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value'));
    axes(handles.axes1)
    imagesc(handles.im2)
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
set(handles.text15,'String',num2str(round(sliderValue))); 

    
if handles.stepnumber==1 % Slider being used for bpass step
    handles.im2=bpass(handles.im,get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value'));
    axes(handles.axes1)
    imagesc(handles.im2)
elseif get(handles.popupmenu1,'Value')==1 % Slider being used for pkfnd threshold
    handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'));
    handles.pks=cntrd(double(handles.im),handles.coords,1+2*round(get(handles.slider3,'Value')));
    axes(handles.axes1)
    imagesc(handles.im2)
    hold on
    scatter(handles.pks(:,1),handles.pks(:,2),'xr')
    hold off
elseif get(handles.popupmenu1,'Value')==2 % Slider being used for pkfnd threshold
    handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'));
    p=handles.coords;
    particle_sz=1+2*get(handles.slider3,'Value');
    p_rs=[];
    %zeros(length(p),3);
    for i=1:length(p)
        xcentre=round(p(i,1));
        ycentre=round(p(i,2));
        if ycentre-(particle_sz-1)/2>0 && ycentre+(particle_sz-1)/2<=length(handles.im(:,1)) && xcentre-(particle_sz-1)/2 >0 && xcentre+(particle_sz-1)/2 <= length(handles.im(1,:))
            testim=handles.im(ycentre-(particle_sz-1)/2:ycentre+(particle_sz-1)/2,xcentre-(particle_sz-1)/2:xcentre+(particle_sz-1)/2);
            [c,d,f]=radialcenter(double(testim));
            p_rs=[p_rs ; xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
            %p_rs(i,:)=[xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
        end
    end
    handles.pks=p_rs;
    axes(handles.axes1)
    imagesc(handles.im2)
    hold on
    scatter(handles.pks(:,1),handles.pks(:,2),'xr')
    hold off
else % Slider being used for binarize threshold
    handles.im2=binarize(handles.im_pass,get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value'));
    axes(handles.axes1)
    imagesc(handles.im2)
end

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
if handles.stepnumber==1 % Slider being used for bpass step
    set(handles.text16,'String',num2str(round(sliderValue)));
elseif get(handles.popupmenu1,'Value')==1 % Slider being used for pkfnd particle size. Note the 2x+1 factor due to requirement to be odd.
    set(handles.text16,'String',num2str(1+2*round(sliderValue)));
elseif get(handles.popupmenu1,'Value')==2 % Slider being used for pkfnd particle size. Note the 2x+1 factor due to requirement to be odd.
    set(handles.text16,'String',num2str(1+2*round(sliderValue)));
else % Slider being used for binarize step
   set(handles.text16,'String',num2str(round(sliderValue))); 
end

if handles.stepnumber==1
    handles.im2=bpass(handles.im,get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value'));
    axes(handles.axes1)
    imagesc(handles.im2)
elseif get(handles.popupmenu1,'Value')==1
	handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'));
    handles.pks=cntrd(double(handles.im),handles.coords,1+2*round(get(handles.slider3,'Value')));
    axes(handles.axes1)
    imagesc(handles.im2)
    hold on
    scatter(handles.pks(:,1),handles.pks(:,2),'xr')
    hold off
elseif get(handles.popupmenu1,'Value')==2

    handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'));
    p=handles.coords;
    particle_sz=1+2*get(handles.slider3,'Value');
    p_rs=[];
    %zeros(length(p),3);
    for i=1:length(p)
        xcentre=round(p(i,1));
        ycentre=round(p(i,2));
        if ycentre-(particle_sz-1)/2>0 && ycentre+(particle_sz-1)/2<=length(handles.im(:,1)) && xcentre-(particle_sz-1)/2 >0 && xcentre+(particle_sz-1)/2 <= length(handles.im(1,:))
            testim=handles.im(ycentre-(particle_sz-1)/2:ycentre+(particle_sz-1)/2,xcentre-(particle_sz-1)/2:xcentre+(particle_sz-1)/2);
            [c,d,f]=radialcenter(double(testim));
            p_rs=[p_rs ; xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
            %p_rs(i,:)=[xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
        end
    end
    handles.pks=p_rs;
    axes(handles.axes1)
    imagesc(handles.im2)
    hold on
    scatter(handles.pks(:,1),handles.pks(:,2),'xr')
    hold off
    
else
    handles.im2=binarize(handles.im_pass,get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value'));
    axes(handles.axes1)
    imagesc(handles.im2)
end

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
if handles.stepnumber==1
    set(handles.text17,'String',num2str(round(sliderValue)));
elseif get(handles.popupmenu1,'Value')==1 % Slider being used for cntrd particle size. Note the 2x+1 factor due to requirement to be odd.
    set(handles.text17,'String',num2str(1+2*round(sliderValue)));
elseif get(handles.popupmenu1,'Value')==2 % Slider being used for radialcenter particle size. Note the 2x+1 factor due to requirement to be odd.
    set(handles.text17,'String',num2str(1+2*round(sliderValue)));
else
   set(handles.text17,'String',num2str(round(sliderValue))); 
end
    
if handles.stepnumber==1
    handles.im2=bpass(handles.im,get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'),get(handles.slider3,'Value'));
    axes(handles.axes1)
    imagesc(handles.im2)
elseif get(handles.popupmenu1,'Value')==1
    handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),get(handles.slider2,'Value'));
	handles.pks=cntrd(double(handles.im),handles.coords,1+2*round(get(handles.slider3,'Value')));
    axes(handles.axes1)
    imagesc(handles.im2)
    hold on
    scatter(handles.pks(:,1),handles.pks(:,2),'xr')
    hold off
elseif get(handles.popupmenu1,'Value')==2

    handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'));
    p=handles.coords;
    particle_sz=1+2*get(handles.slider3,'Value');
    p_rs=[];
    %zeros(length(p),3);
    for i=1:length(p)
        xcentre=round(p(i,1));
        ycentre=round(p(i,2));
        if ycentre-(particle_sz-1)/2>0 && ycentre+(particle_sz-1)/2<=length(handles.im(:,1)) && xcentre-(particle_sz-1)/2 >0 && xcentre+(particle_sz-1)/2 <= length(handles.im(1,:))
            testim=handles.im(ycentre-(particle_sz-1)/2:ycentre+(particle_sz-1)/2,xcentre-(particle_sz-1)/2:xcentre+(particle_sz-1)/2);
            [c,d,f]=radialcenter(double(testim));
            p_rs=[p_rs ; xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
            %p_rs(i,:)=[xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
        end
    end
    handles.pks=p_rs;
    axes(handles.axes1)
    imagesc(handles.im2)
    hold on
    scatter(handles.pks(:,1),handles.pks(:,2),'xr')
    hold off
    
else    
    handles.im2=binarize(handles.im_pass,get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value'));
    axes(handles.axes1)
    imagesc(handles.im2)
end

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



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Button to end
uiresume(handles.figure1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is the generate figure button
handles = guidata(hObject);
if get(handles.popupmenu1,'Value')==1
    %handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),get(handles.slider2,'Value'));
    coords=handles.coords;
    im=double(handles.im);
    sz=1+2*get(handles.slider3,'Value');
    handles.pks=cntrd(im,coords,sz);
    handles.pkfnd_params=[get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'),1+2*get(handles.slider3,'Value')];
elseif get(handles.popupmenu1,'Value')==2
    handles.coords=pkfnd(handles.im_pass,get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'));
    p=handles.coords;
    particle_sz=1+2*get(handles.slider3,'Value');
    p_rs=[];
    %zeros(length(p),3);
    for i=1:length(p)
        xcentre=round(p(i,1));
        ycentre=round(p(i,2));
        if ycentre-(particle_sz-1)/2>0 && ycentre+(particle_sz-1)/2<=length(handles.im(:,1)) && xcentre-(particle_sz-1)/2 >0 && xcentre+(particle_sz-1)/2 <= length(handles.im(1,:))
            testim=handles.im(ycentre-(particle_sz-1)/2:ycentre+(particle_sz-1)/2,xcentre-(particle_sz-1)/2:xcentre+(particle_sz-1)/2);
            [c,d,f]=radialcenter(double(testim));
            p_rs=[p_rs ; xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
            %p_rs(i,:)=[xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
        end
    end
    handles.pks=p_rs;
    handles.radcent_params=[get(handles.slider1,'Value'),1+2*get(handles.slider2,'Value'),1+2*get(handles.slider3,'Value')];

else
    
    regions = bwlabel(handles.im2,8);
    props = regionprops(regions,handles.im,'WeightedCentroid','Area'); % find features in BW (Particles)

    X=[]; Y=[]; Area=[];
    for n = 1:1:size(props)    % go through features in BW
        % Particles:
            X = [X; props(n).WeightedCentroid(1)];
            Y = [Y; props(n).WeightedCentroid(2)];
            Area = [Area; props(n).Area];   
    end

    handles.pks=cat(2,X,Y,Area);
    handles.regionprops_params=[get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value')];
end

figure
imagesc(handles.im);
hold on
plot(handles.pks(:,1),handles.pks(:,2),'rx')
if get(handles.popupmenu1,'Value')==1
    title(['pkfnd cntrd, params: ',num2str(handles.pkfnd_params(1)),' ',num2str(handles.pkfnd_params(2)),' ',num2str(handles.pkfnd_params(3))])
elseif get(handles.popupmenu1,'Value')==2
    title(['pkfnd radialcenter, params: ',num2str(handles.radcent_params(1)),' ',num2str(handles.radcent_params(2)),' ',num2str(handles.radcent_params(3))])
else
    title(['Regionprops with WeightedCentroid, params: ',num2str(handles.regionprops_params(1)),' ',num2str(handles.regionprops_params(2)),' ',num2str(handles.regionprops_params(3))])
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is the finished bpass button. Uncover the other options as necessary
% and move onto the pkfnd etc part

handles = guidata(hObject);
% Record the bpass parameters the final output
handles.im_pass=handles.im2;
handles.bpass_params=[get(handles.slider1,'Value'),get(handles.slider2,'Value'),get(handles.slider3,'Value')];
% Hide all the options except the option to choose the locating method.
set(handles.slider1,'Visible','Off')
set(handles.slider2,'Visible','Off')
set(handles.slider3,'Visible','Off')
set(handles.text4,'Visible','Off')
set(handles.text5,'Visible','Off')
set(handles.text6,'Visible','Off')
set(handles.text9,'Visible','Off')
set(handles.text10,'Visible','Off')
set(handles.text11,'Visible','Off')
set(handles.text12,'Visible','Off')
set(handles.text13,'Visible','Off')
set(handles.text14,'Visible','Off')
set(handles.text15,'Visible','Off')
set(handles.text16,'Visible','Off')
set(handles.text17,'Visible','Off')
set(handles.pushbutton3,'Visible','Off')

set(handles.text3,'Visible','On')
set(handles.popupmenu1,'Visible','On')

% Tell it that you've moved onto the second step
handles.stepnumber=2;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function text15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% This is the title


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1);
%delete(hObject);
