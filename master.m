% Master script to process double-dyed PDMS confocal data for 
% the Adhesion Induced Phase Separation project @ Williams Collage Physics

% Release: v.1.0 (July 2019)
% (c) Josh Kang 

% Original script described in K. E. Jensen & N. Nakamura (RSI 2016), and also described in
% Katharine E. Jensen's dissertation ("Structure and defects of hard-sphere
% colloidal crystals and glasses", Harvard University, 2013; section 2.3).

%% Step 0. Set Variables and Constants

disp('------------------------------------------------------------------')
disp('             Welcome to AIPS analysis v.1.0  ')
disp('------------------------------------------------------------------')
disp([newline 'Getting started ...' ])

% 0-1) INTERACTIVE settings update

% Ask user for module name 
module_name = input([newline 'Set current module name: ' ],'s');
file_type = '.tif'; %default setting

% Should we run interactively?
run_interactively = input([newline 'Run interactively? (y/n) :'],'s');
    if strcmp(run_interactively, 'y')
        run_interactively = 1;
    else 
        run_interactively = 0;
    end
    
% Ask for scale of confocal image data 
disp([newline 'What were the xyz scales (in �m) of your images? (Enter in the format [ __ , __ , __ ])']);
scale = input([': ']);

% Ask for image location
    if ~exist('image_folder','var')
        image_folder = input([newline 'No data folder defined. Read from current directory? (y/n) '],'s');
            if strcmp (image_folder, 'y')
                image_folder = pwd;
            else 
                image_folder = input('Enter image location directory:' , 's');
            end
    end
disp([newline 'Data will be read from the directory: ' image_folder newline])


% AUTOMATED settings update

% Remember where code is located at
code_folder = pwd; 

% Read all image files from directory and sort alphabetically
% TODO: set search keyword according to OME TIF RGB output names
cd(image_folder);
addpath(image_folder);

image_beads_list = strsplit(ls('*.tif'));
image_bckgd_list = strsplit(ls('*.tif'));
image_beads_list(end) = []; 
image_bckgd_list(end) = [];
image_beads_list = sort(image_beads_list);
image_bckgd_list = sort(image_bckgd_list);
cd(code_folder);

% How many samples did we find?
image_count = length(image_beads_list);
disp(['There were a total of ' num2str(image_count) ' sets of image data found in ' image_folder])

% User-defined settings for iterative locating (for info look at iterative_locating_parameters.m) 
false_particle_size = [9 9 15];
range = [1 1e5]; %to load in all available image slices, set range(2) to a large number (or Inf); software 
                % will automatically move on when it runs out of images to load

%the approximate size of the image noise, in pixels (x, y,and z directions):
lnoise = [2 2 4];
%the approximate particle diameter in pixels (x, y, and z directions):
lobject = [15 15 30];

%numbers somewhat larger than the largest particle diameter:
diameters = [16 16 32];
%numbers somewhat smaller than the smallest particle diameter:
mask_size = [7 7 15];

%a minimum separation cutoff to help prevent doubly-located particles 
%(set conservatively; this is to avoid a single particle being detected
%twice in any given particle locating iteration)
min_separation = [5 5 10];

masscut_initial = 2e5; % first iteration generally only finds real particles to start with
masscut_residuals = 1e5; % required nonzero to prevent residual noise from being detected as "particles" once all of the real particles have been found

% Arrays of data to be saved and will be used throughout this analysis
r_um_stack = [];
r_zeroed_stack = [];
beads_Rd_stack = [];
bckgd_rz_stack = [];

disp(['[Success] finished configuring initial parameters ...' newline])


%% Step 1. Locating fluorescent beads

disp('------------------------------------------------------------------')
disp('             Step 1: Locating fluorescent beads')
disp('------------------------------------------------------------------')

for counter = 1:image_count
    current_filename = image_beads_list{counter};
    disp(['Starting particle locating for ' current_filename newline])
% 1-1) Preliminaries - check that needed variables are defined, and save out the locating parameters
    if ~exist('invert_image','var')
        invert_image = 0; %default is NOT to invert the image; this is appropriate for dyed particles; for dye in the fluid phase, invert the image!
    end

% 1-2) Load in the raw data image stack
    disp(['Image processing and particle locating started at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '.'])
    disp([newline 'reading in files...'])
    clear raw %ensure doesn't already exist in the workspace

    try %if range is set to be bigger than the actual stack, will read in files until there are no more z-slices, then continue with locating
        for i=1:1+range(2)-range(1)
            I = imread([image_folder '/' current_filename],'Index',range(1)+i-1);
            raw(:,:,i) = I;
        end
    catch %if errors because it's out of images, just continue with what we have
        disp(['There were ' num2str(i-1) ' images in stack ' current_filename])
    end

    if invert_image
        raw = 255 - raw; %invert the 8-bit image; use when the fluorescent dye is in the fluid phase, not in the particles
    end
    disp(['Image loading finished at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '.'])

    if run_interactively
        beep; pause(0.2); beep
        % if running interactively, scan through the images by eye to check
        % that they were loaded properly:
        figure; for i=1:size(raw,3); imagesc(raw(:,:,i)); axis image; pause(0.1); end
        
        % alternate version: if you would prefer to go through the images one
        % by one (hitting enter to go to the next image), uncomment the following:
        %figure; for i=1:size(raw,3); imagesc(raw(:,:,i)); axis image; input(''); end
    end

    % 1-3) Processing step 1: bandpass filter remove high-frequency noise and flatten the image background

    disp([newline 'bandpass filtering the raw data...'])
    res = bpass3dMB(raw, lnoise, lobject, [0 0]); %the [0 0] just means that there are additional options in the function we are not using

    disp(['Image bandpass filtering finished at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '.'])

    if run_interactively
        beep; pause(0.2); beep
        % if running interactively, look through the filtered images to see how they look compared to the original:
        figure; for i=1:size(raw,3); subplot(2,1,1); imagesc(raw(:,:,i)); axis image; subplot(2,1,2); imagesc(res(:,:,i)); axis image; pause(0.1); end
        % also scan the SIDE view (y-z) (very important to judge the z-filtering):
        figure; for i=1:size(raw,1); subplot(2,1,1); imagesc(permute(raw(i,:,:),[3 2 1])); axis image; subplot(2,1,2); imagesc(permute(res(i,:,:),[3 2 1])); axis image; pause(0.1); end
    end

    %1-4) Processing step 2: locate the particles from the filtered images
    disp('feature finding...')
    r=feature3dMB(res, diameters, mask_size, size(res), [1 1 0], min_separation, masscut_initial, 0); %the [1,0,0] here means we are using only the first optional setting: a minimum separation

    % during manual parameter-refining, at the end of this step, use show_slice to overlay the particle coordinates on the original images to check the locating.
    %save r_firstpass.mat 
    disp(['Initial particle locating complete at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '. Moving on to find any missed particles by running particle locating iteratively on the image residuals.'])

    % 1-5) Save the particle locations

    disp('recording coordinates...')
    
    % saves as 2-column, tab-delimited text file
    disp(['Particle locations saved in file ' module_name '_beads_sphere' num2str(counter) '_xyz_coordinates.txt'])
    disp(['Particle locating for ' current_filename ' entirely finished at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '.' newline])

    % convert xyz coordinates accordingly and save them in a txt file
    r3 = r(:,1:3);
    r_um = scale .* r3;
    scatter3(r_um(:,1),r_um(:,2),r_um(:,3),20,r_um(:,3),'filled')
    axis equal
    
    % save r_um
    r_um_stack(:,:,counter) = r_um;
    dlmwrite([module_name '_bead_sphere' num2str(counter) '_xyz_coordinates.txt'], r_um, '\t');
end

% Setup stacks to save data that have dimensions dependent on 'image_count'
center_stack = zeros(1,2,image_count);
rz_stack = zeros(1,2,image_count);

%% Step 2. Level and Collapse Flourescent Beads & Flourescent Background PDMS Data 
% @dependencies: select_region.m, fit_poly22_to_points.m
% @params: r_um
% @returns: select_vector, undeformed_plane_fit,zero_surface_to_subtract, r_zeroed
% @writes: ..._r_zeroed.txt

disp(newline)
disp('------------------------------------------------------------------')
disp('         Step 2: level and rz-collapse data')
disp('------------------------------------------------------------------')
    
for counter = 1:image_count
    current_filename = image_beads_list{counter};
    r_um = r_um_stack(:,:,counter);
    
    disp(['Leveling and collapsing bead data for: ' current_filename newline])
    
    % 2-1) zeroing / leveling...
    % zero/level the data by fitting a plane to the OUTER points, far from the
    % indenting sphere
    %
    % In this case, we always put the sphere in the top-left-ish, so the
    % bottom (high row number, large x) and right (high column number, large y)
    % edges are where the data are flat. We'll use the points there to sample
    % the flat plane.

    % If that changes, you'll want to adjust the select_region inputs below.
    % See >>help select_region for instructions on that function.

    %only fit the edges of the image so that get a sampling of the
    %undeformed plane around the particle
    dist_from_edge = 10; %CONST

    %for a symmetric box:     
    %select_vector = ~select_region(r,'box',[dist_from_edge dist_from_edge -inf],[max(r(:,1))-dist_from_edge max(r(:,2))-dist_from_edge inf]);
    %for only the right and bottom edges:
    select_vector = select_region(r_um,'box',[0 0 -inf],[max(r_um(:,1))-dist_from_edge max(r_um(:,2))-dist_from_edge inf]);

    %undeformed_plane_fit = fit_plane_to_points(r(~select_vector,:));
    undeformed_plane_fit = fit_poly22_to_points(r_um(~select_vector,:)); %works better - takes out any curvature in surface

    %generate zero surface to subtract at every point:
    zero_surface_to_subtract = undeformed_plane_fit.p00 + undeformed_plane_fit.p10*r_um(:,1) + undeformed_plane_fit.p01*r_um(:,2)...
        + undeformed_plane_fit.p20*r_um(:,1).^2 + undeformed_plane_fit.p11*r_um(:,1).*r_um(:,2) + undeformed_plane_fit.p02*r_um(:,2).^2;

    %subtract to level profile:
    r_zeroed = [r_um(:,1:2) -zero_surface_to_subtract + r_um(:,3)];

    %and save, put away
    r_zeroed_stack(:,:,counter) = r_zeroed;
    dlmwrite([module_name '_beads_sphere' num2str(counter) '_r_zeroed.txt'],r_zeroed,'\t')


    % 2-2) Now do the semi-manual collapse
    % this requires manual inputs -- follow the instructions on the screen!
    [center, rz] = radial_collapse_manual(r_zeroed);

    %save the center and rz data for this sphere as text and to stack:
    center_stack(:,:,counter) = center; rz_stack(:,:,counter) = rz;
    dlmwrite([module_name '_beads_sphere' num2str(counter) '_rz_collapsed.txt'],rz,'\t')

    % 2-4) make a nice plot of the collapsed data for this sphere
    figure('Name', [module_name ': beads collapse for sphere' num2str(counter)])
    plot(rz(:,1),rz(:,2),'.','DisplayName',['profile ' module_name '_sphere' num2str(counter)])
    hold all
    %make the plot that looks nice:
    set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on
    axis equal
    xlabel('r (�m) ')
    ylabel('z (�m) ')
    savefig([module_name '_beads_sphere' num2str(counter) '_collapsed_figure.fig'])
    
    % 2-5) Level and collapse background (Nile Red) PDMS 
    [leveled,zero] = background_level_z(raw, scale, undeformed_plane_fit);
    bg = background_rz_collapse(leveled, scale, center, zero);
    bckgd_rz_stack(:,:,counter) = bg; % save the leveled & rz-collapsed bckgd data
    
    % create plot for comparison/overlay
    figure('Name', ['Complete rz collapse for ' module_name '_sphere' num2str(counter)])
    plot(bg(:,1),bg(:,2),'.','DisplayName',['profile ' module_name '_beads_sphere' num2str(counter)])
    hold all
    plot(rz(:,1),rz(:,2),'.','DisplayName',['profile ' module_name '_bckgd_sphere' num2str(counter)])
    % make the plot that looks nice:
    set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on
    % axis equal
    axis ([0, 60, -15,5])
    xlabel('r (\mu m) ')
    ylabel('z (\mu m) ')
    savefig([module_name '_sphere' num2str(counter) '_overlay_figure.fig'])

    % 2-6) visualize-side view for sanity check
    disp('Santiy check...')
    for i=290:1024;
        imagesc(permute(leveled(i,:,:)>220,[3 2 1])); 
        hold all; 
        plot(rz(:,1)/0.14 + 330, rz(:,2)/0.225 + 61,'r.');
        hold off; 
        axis image; 
        axis xy; 
        title(i); 
        input(''); 
    end
end

%% Step 3. Circle fit the collapsed bead coordinates to find sphere Radius and (indentation) depth
% @dependencies: cicle_fit.m, circle_fit_error.m
% @params: rz, 
% @returns: profile_fit, r0z0R_fit, Rvsd
% @writes: r vs d.mat

disp(newline)
disp('------------------------------------------------------------------')
disp('       Step 3: Circle fit bead coordinates ')
disp('------------------------------------------------------------------')

for counter = 1:image_count 
    current_filename = image_beads_list{counter};
    rz = rz_stack(:,:,counter);
    disp(['Circle-fitting bead data for: ' module_name ' sphere ' num2str(counter)  newline])
    
    % 3-1) determine (manually, looking at a plot) the inner region (which will be fit to a circle)
    figure('name','Confocal profiles')
    plot(rz(:,1),rz(:,2),'.')
    grid on; box on
    set(gca,'LineWidth',1,'FontSize',16,'FontWeight','bold')

    xlabel('Distance from indenter center (µm)','FontSize',18,'FontWeight','bold')
    ylabel('Surface profile (µm)','FontSize',18,'FontWeight','bold')
    xlim([0 40])
    % this asks you to type in a number - what radial distance is definitely
    % just the inner circle?
    inner_region_limit = input('Inner region r-cutoff: ');

    % 3-2) Now fit a circle to the inner region of these to get the sphere position
    % select only those points inside the desired R range, and mirror over the y-axis:
    profile_to_fit = double(rz(~(rz(:,1)>inner_region_limit),:));
    profile_to_fit = [[-profile_to_fit(:,1) profile_to_fit(:,2)]; profile_to_fit]; %now mirrored

    % add these data dto the plot, highlighted:
    hold all
    plot(profile_to_fit(:,1), profile_to_fit(:,2),'.')

    % Circle fit the profile
    %try
    r0z0R_fit = circle_fit(profile_to_fit);
    %catch
    %    warning('error: data could not be fit')
    %    r0z0R_fit = [0 0 0];
    %end
    disp('Radius is: ')
    r0z0R_fit(1,3)
    disp('indentation depth is: ')
    r0z0R_fit(:,3)-r0z0R_fit(:,2)

    % Save R vs d data to array
    R_d = [r0z0R_fit(:,3) (r0z0R_fit(:,3)-r0z0R_fit(:,2))];
    beads_Rd_stack = [beads_Rd_stack;R_d];

    % visually check if circle fit was good
    hold all
    viscircles(r0z0R_fit(1,1:2),r0z0R_fit(1,3),'Color','k'); 
    axis equal

    % (3-3) Save R vs d data for future use (compiled figure)
    % change variable name !!!
    save([module_name '_beads_Rd.mat'], beads_Rd_stack);
end

%% Step 4. Calculate fluid phase volume
disp(newline)
disp('------------------------------------------------------------------')
disp('           Step 4: Calculate fluid phase PDMS volume ')
disp('------------------------------------------------------------------')


%% Step 5. Plot Final R vs d

% (5-1) Produce a d vs R plot based on data you just analyzed
% set appropriate name for figure
figure('name', [module_name '_' image_filename ' R vs d'])
scatter(final(:,1),final(:,2),'filled')
grid on; box on
set(gca,'LineWidth',1,'FontSize',12,'FontWeight','bold')

ylabel('Distance from Indenter Center, d (�m)','FontSize',14,'FontWeight','bold')
xlabel('Radius of Spehre, R (�m)','FontSize',14,'FontWeight','bold')
xlim([0 40])
ylim([0 20])

