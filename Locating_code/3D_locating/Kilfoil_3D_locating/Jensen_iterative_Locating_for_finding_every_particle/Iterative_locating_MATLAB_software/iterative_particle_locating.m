% iterative_particle_locating.m
%
% This is a script that implements the iterative particle locating algorithm
% described in K. E. Jensen & N. Nakamura (RSI 2016), and also described in
% Katharine E. Jensen's dissertation ("Structure and defects of hard-sphere
% colloidal crystals and glasses", Harvard University, 2013; section 2.3).
% The code below goes through the following steps:
%
% (1) Preliminaries - check that needed variables are defined, and save out the locating parameters
% (2) Load in the raw data image stack
% (3) Processing step 1: bandpass filter remove high-frequency noise and flatten the image background
% (4) Processing step 2: locate the particles from the filtered images
% (5) Processing step 3: Find any missed particles by subtracting the found particles and iterating the previous steps
% (6) Save the particle locations
%
% This software specifically interfaces with the MATLAB particle locating
% software developed by Gao and Kilfoil (see Y. Gao and M.L. Kilfoil.
% Opt. Express vol. 17, 4685 (2009)). The Gao and Kilfoil software is
% publicly available under ``MATLAB 3D feature finding algorithms'' at
% http://people.umass.edu/kilfoil/downloads.html
%
% However, the iterative algorithm implemented here is compatible with any
% particle locating software, and it should be straightforward to adapt
% this script for use with your favorite locating routines. In this case,
% skip down to Step 5 below where the iteration procedure is implemented.
% You will need one additional set of parameters for the iteration, which
% is the "false_particle_size" x-, y-, and z-values. These are used to
% determine how much of the image is deleted around each found particle
% prior to iteration.
%
% Input parameters for locating and iteration should be in the MATLAB
% workspace before this script is called. The script
% interative_particle_locating_input_parameters.m can be used to set these
% parameters, or they can be set manually prior to calling this script.
%
% In choosing the parameters for particle locating, one should err on the
% side of missing particles on any given locating iteration rather than
% risk double-counting a particle. These missed particles will be found
% during later iterations.
%
% Please direct any questions or comments to Katharine Jensen at
% kjensen@post.harvard.edu
%
% Algorithm developed by Nobotomo Nakamura and Katharine Jensen in 2011 and early 2012.
% Cleaned up and prepared for distribution in June 2015.



%% (1) Preliminaries - check that needed variables are defined, and save out the locating parameters

display('Checking whether the required parameters are defined, and saving a record of what parameters were used...')

run_interactively = 1; %if set to 1, will display figures and stacks for the user to view; if set to 0, nothing will be displayed

if ~exist('run_test_stack','var')
    %default is that we're running the entire image stack; but can use
    %run_test_stack to crop to a smaller image volume for checking input
    %parameters quickly
    run_test_stack = 0;
end

if ~exist('invert_image','var')
    invert_image = 0; %default is NOT to invert the image; this is appropriate for dyed particles; for dye in the fluid phase, invert the image!
end

%generate the output filename based on the image filename, removing the image file extension
output_filename = [image_filename(1:find(image_filename == '.',1,'last')-1) '_xyz_coordinates.txt'];

%if needed, set the "masscut" (minimum integrated intensity to be
%considered a real particle) parameters to default values, and warn:
if ~exist('masscut_initial','var')
    warning('The parameter masscut_initial was not defined. It will be set to a default value, but it is strongly advised that the user define this parameter in the future.')
    masscut_initial = 3e4;
end

if ~exist('masscut_residuals','var')
    warning('The parameter masscut_residuals was not defined. It will be set to a default value, but it is strongly advised that the user define this parameter in the future.')
    masscut_residuals = 1e5;
end

%save as a .mat file (and check that everything is in place):
save(['locating_parameters_' datestr(now,'YYmmDD_hhMMss') '.mat'],'range','image_filename','folder','lnoise','lobject','diameters','mask_size','min_separation','masscut_initial','masscut_residuals','false_particle_size','output_filename','invert_image','run_test_stack')
%this save statement will error out if everything isn't properly defined
%these parameters can then be reloaded directly to re-run this script in the future 


%% (2) Load in the raw data image stack

display(['Image processing and particle locating started at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '.'])

display('reading in files...')
clear raw %ensure doesn't already exist in the workspace

try %if range is set to be bigger than the actual stack, will read in files until there are no more z-slices, then continue with locating
    for i=1:1+range(2)-range(1)
        I = imread([folder '/' image_filename],'Index',range(1)+i-1);
        if run_test_stack
            I = I(1:100,1:100); %optionally, crop images for testing a small x-y region only, as when refining the locating parameters
        end
        raw(:,:,i) = I;
    end
catch %if errors because it's out of images, just continue with what we have
    display(['There were ' num2str(i-1) ' images in stack ' image_filename])
end

if invert_image
    raw = 255 - raw; %invert the 8-bit image; use when the fluorescent dye is in the fluid phase, not in the particles
end

display(['Image loading finished at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '.'])


if run_interactively
    beep; pause(0.2); beep
    % if running interactively, scan through the images by eye to check
    % that they were loaded properly:
    figure; for i=1:size(raw,3); imagesc(raw(:,:,i)); axis image; pause(0.1); end
    
    % alternate version: if you would prefer to go through the images one
    % by one (hitting enter to go to the next image), uncomment the following:
    %figure; for i=1:size(raw,3); imagesc(raw(:,:,i)); axis image; input(''); end
end



%% (3) Processing step 1: bandpass filter remove high-frequency noise and flatten the image background

display('bandpass filtering the raw data...')
res=bpass3dMB(raw, lnoise, lobject, [0 0]); %the [0 0] just means that there are additional options in the function we are not using

display(['Image bandpass filtering finished at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '.'])


if run_interactively
    beep; pause(0.2); beep
    % if running interactively, look through the filtered images to see how they look compared to the original:
    figure; for i=1:size(raw,3); subplot(2,1,1); imagesc(raw(:,:,i)); axis image; subplot(2,1,2); imagesc(res(:,:,i)); axis image; pause(0.1); end
    % also scan the SIDE view (y-z) (very important to judge the z-filtering):
    figure; for i=1:size(raw,1); subplot(2,1,1); imagesc(permute(raw(i,:,:),[3 2 1])); axis image; subplot(2,1,2); imagesc(permute(res(i,:,:),[3 2 1])); axis image; pause(0.1); end
end



%% (4) Processing step 2: locate the particles from the filtered images

display('feature finding...')
r=feature3dMB(res, diameters, mask_size, size(res), [1 1 0], min_separation, masscut_initial, 0); %the [1,0,0] here means we are using only the first optional setting: a minimum separation

% during manual parameter-refining, at the end of this step, use show_slice to overlay the particle coordinates on the original images to check the locating.
save r_firstpass.mat r
display(['Initial particle locating complete at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '. Moving on to find any missed particles by running particle locating iteratively on the image residuals.'])



%% (5) Processing step 3: Find any missed particles by subtracting the found particles and iterating the previous steps

% Create new "raw" data by replacing the already-found particles with a
% virtual/false particle of all zeros, thus erasing those particles from
% the original image. The remaining image only contains the residuals:
% those particles that have not yet been found by the particle locating.
% With the easily-found particles removed, these [usually sparse] remaining
% particles are now straightforward to locate.
%
% The method iterates until no new particles are found. The number of
% iterations required depends on the quiality original images. For
% evenly-illuminated image stacks containing ~50,000 particles, we find
% that the algorithm usually converges in about 4-6 iterations. For image
% stacks in which the particles vary significantly in brightness (either
% due to uneven illumination or variation in the dye concentration in each
% particle), more iterations are required.


%%%%%%%%%%%%%%% 
% Step A: Generate the "false_particle" of all zeros:

%The false_particle_size parameters must be odd integers so that the
%ellipsoid of zeros has a well-defined center:
false_particle_size = round(false_particle_size); %ensure integers
false_particle_size = false_particle_size + ~mod(false_particle_size,2); %ensure odd numbers

%Break out the size parameters so the code below is more easily human-readable: 
a = false_particle_size(1)/2;
b = false_particle_size(2)/2;
c = false_particle_size(3)/2;

% default the matrix elements to 1 so that they will leave the image unchanged when applied:
false_particle = ones(false_particle_size(1),false_particle_size(2),false_particle_size(3));

%created a virtual particle of zeros in the center of this matrix:
for i=1:false_particle_size(1)
    for j=1:false_particle_size(2)
        for k=1:false_particle_size(3)
            v = [i j k]-(false_particle_size./2+1/2); %subtract the center coordinates so looking radially out from center
            % volume inside an ellipse:
            if ((v(1)/a)^2 + (v(2)/b)^2 + (v(3)/c)^2)  <= 1
                false_particle(i,j,k) = 0; %so will erase when multiplied in
            end
        end
    end
end

raw_residuals = double(raw); %make a new "raw" data set while retaining the original images (this may be questionable wrt memory usage, but useful for debugging)

%how many pixels in each direction is the virtual particle matrix?
x_range = -(a-1/2):a-1/2;
y_range = -(b-1/2):b-1/2;
z_range = -(c-1/2):c-1/2;



%%%%%%%%%%%%%%% 
% Step B: Delete the found particles, and iterate the particle locating:

residual_start_index = 1; %keeps track of which found particles we've already deleted from the residuals image stack
still_searching = 1; %if == 1, particles are still being found 
locating_iteration = 2;

while still_searching == 1 %there are still particles to find!
    display(['Beginning particle locating iteration ' num2str(locating_iteration) '...'])
    
    %Use the false_particle matrix created above to delete the raw date everywhere a particle was found:
    %(this requires careful bookkeeping of matrix indices, especially close to the edges of the image)
    for i=residual_start_index:size(r,1) 
        this_particle_center = round(r(i,1:3)); %recall that the coordinates are still in pixels
        
        %figure out the ranges; and remove any array-out-of-bounds problems:
        this_xrange = this_particle_center(1) + x_range;
        invalid_xrange = this_xrange<1 | this_xrange>size(raw,1);
        this_xrange(invalid_xrange) = [];
        
        this_yrange = this_particle_center(2) + y_range;
        invalid_yrange = this_yrange<1 | this_yrange>size(raw,2);
        this_yrange(invalid_yrange) = [];
        
        this_zrange = this_particle_center(3) + z_range;
        invalid_zrange = this_zrange<1 | this_zrange>size(raw,3);
        this_zrange(invalid_zrange) = [];
        
        %remove this already-found particle from the image by multiplying in the virtual particle matrix:
        raw_residuals(this_xrange, this_yrange, this_zrange) = ...
            raw_residuals(this_xrange, this_yrange, this_zrange) .* false_particle(~invalid_xrange,~invalid_yrange,~invalid_zrange);
    end 
    %found particles have now all been deleted
    
    if run_interactively
        beep; pause(0.2); beep
        % if running interactively, scan through the new raw data by eye to check to see how they look:
        figure; for i=1:size(raw_residuals,3); imagesc(raw_residuals(:,:,i)); axis image; pause(0.1); end
        % or go through the images one by one (hit ENTER to go to the next image):
        %figure; for i=1:size(raw_residuals,3); imagesc(raw_residuals(:,:,i)); axis image; input(''); end
    end
    
    
    % Next, locate any particles that remain, using the same settings to filter
    % and locate as before, but setting a strict "masscut" threshold for
    % the minimum integrated intensity required for a bright spot to be
    % considered a "real" particle. This parameter is essential so that the
    % iterations knows when all real particles have been found, and stops
    % iterating at that point.
    
    display('bandpass filtering the residual raw data...')
    res_residuals=bpass3dMB(raw_residuals, lnoise, lobject, [0 0]);
    
    if run_interactively
        beep; pause(0.2); beep
        % look through the new filtered images to see how they look, both
        % top-view and side-view (y-z):
        figure; for i=1:size(raw_residuals,3); imagesc(res_residuals(:,:,i)); axis image; pause(0.1); end
        for i=1:size(raw_residuals,1); imagesc(permute(raw_residuals(i,:,:),[3 2 1])); axis image; pause(0.1); end
        % or go through the images one by one (hit ENTER to go to the next image):
        %figure; for i=1:size(raw_residuals,3); imagesc(res_residuals(:,:,i)); axis image; input(''); end
    end
    
    
    display('feature finding additional particles...')
    r_residuals=feature3dMB(res_residuals, diameters, mask_size, size(res_residuals), [1 1 0], min_separation, masscut_residuals, 0);
    
   
    if isempty(r_residuals) %we're done; no additional particles were found
        still_searching = 0;
    else %there are still particles to be found
        %add the new particles to the master list, and increment the locating iteration:
        residual_start_index = size(r,1)+1; %keep track of which are the new particles to be removed
        r = [r; r_residuals]; %add the new ones to the list of all found particles
        locating_iteration = locating_iteration + 1;
        %newly-found particles will be erased from raw_residuals at the start of the next loop
    end
    
end %of iterative locating



%% (6) Save the particle locations

display('recording coordinates...')

% saves as 2-column, tab-delimited text file
dlmwrite(output_filename, r, '\t');

display(['Particle locations saved in file ' output_filename])
display(['Particle locating for ' image_filename ' entirely finished at ' datestr(now,'HH:MM:ss') ' on ' datestr(now, 'mm-DD-YYYY') '.'])



%% (7) Enjoy your particle locations!

% KEJ recommends doing some eyeball-checking at this point of the particle
% locations overlaid on the images to make sure everything looks good. From
% KEJ's dissertation software, show_slice.m will do this.
%
% It's also not a bad idea to run through the locations and check for
% any double-found particles (i.e. particles that would be
% overlapping; these should be very rare but still happen), and consolidate
% any such pairs into a single particle location. From KEJ's disseration
% software, remove_and_consolidate_double_hits.m will do this.

