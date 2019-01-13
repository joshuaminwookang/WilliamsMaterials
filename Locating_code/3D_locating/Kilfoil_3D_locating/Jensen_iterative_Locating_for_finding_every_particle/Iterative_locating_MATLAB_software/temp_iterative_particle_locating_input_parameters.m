% iterative_particle_locating_input_parameters.m
% 
% This is a helper script designed to get all of the input parameters
% needed to run iterative_particle_locating.m into the MATLAB workspace,
% and to ensure that the locating software is in MATLAB's current directory
% structure (that is, MATLAB knows where the software is). It is designed
% to be run just prior to running iterative_particle_locating.m, but is not
% required for iterative_particle_locating to work. However, it is helpful
% to make sure all needed variables are defined in the workspace.
%
% This software specifically interfaces with the MATLAB particle locating
% software developed by Gao and Kilfoil (see Y. Gao and M.L. Kilfoil.
% Opt. Express vol. 17, 4685 (2009)). The Gao and Kilfoil software is
% publicly available under ``MATLAB 3D feature finding algorithms'' at
% http://people.umass.edu/kilfoil/downloads.html
%
% However, the iterative algorithm implemented here is compatible with any
% particle locating software. The only parameter set here that is critical
% to the iterative locating algorithm is the "false_particle_size". This
% vector of x-, y-, and z- particle diameters is used to determine how much
% of the original image stack will be deleted around each found particle 
% prior to iteration.
%
% In choosing the parameters for particle locating, one should err on the
% side of missing particles on any given locating iteration rather than
% risk double-counting a particle. These missed particles will be found
% during later iterations. Example parameters defined in this script were
% appropriate for Example 1 (colloid glass) included in the Supplement to
% Jensen & Nakamura (2015).
%
% For a more detailed description of the method of iterative particle
% locating, see K. E. Jensen & N. Nakamura (RSI 2015). It is also described in
% Katharine E. Jensen's dissertation ("Structure and defects of hard-sphere
% colloidal crystals and glasses", Harvard University, 2013; section 2.3).
%
% Katharine Jensen - cleaned up and prepared for distribution June 2015


%% Ensure that the locating software is on the current Matlab path structure:
% these paths should be edited to reflect the correct paths on your computer 

addpath ~/Documents/Matlab/ 
addpath '~/Documents/MATLAB/git_groupcode/Locating_code/3D_locating/Kilfoil_3D_locating/Jensen_iterative_Locating_for_finding_every_particle/Iterative_locating_MATLAB_software'
addpath '/Users/skatrina/Documents/MATLAB/git_groupcode/Locating_code/3D_locating/Kilfoil_3D_locating'



%% User-defined setting for iterative locating

%for creating residual images, define the extent of the virtual/false particle
%that will be used to erase already-found particles from the image data.
%(these should be odd integers; if not, the software will automatically round up to the nearest odd integer)
false_particle_size = [8 8 80];



%% User-defined settings for image loading and data saving:

%whether to run iterative locating
single_pass_locating = 1; %if set to 1, will NOT do any iterative locating, but rather will save results of single-pass locating

% whether the image will be inverted before processing; set to 1 when the dye is in the fluid phase, not in the particles: 
invert_image = 0; 

% the range of z-levels (aka image slices) that you want to load in:
range = [1 19000000]; %to load in all available image slices, set range(2) to a large number (or Inf); software will automatically move on when it runs out of images to load

% file name of 3D image stack that will be processed:
% (this can be done a few ways; examples are below)

%manually define the file name:
image_filename = 'sampleA_Capture10_series.tif';
%image_filename = 'Example2_crystal.tif';

%or, uncomment the following to assemble the file name automatically:
%(this is useful, for example, when this script will be called
%automatically by higher-level script that defines the time step for each stack)
% % 
% % image_filename_start = 'Image_Stack';
% % 
% % if exist('timestep','var')
% %     timestep_string = ['_t' sprintf('%03d',timestep)]; %generate the string automatically; set to '%02d' to get 2-digit timesteps
% % end
% % 
% % image_filename_end = '.tif'; %whatever is applicable for your data
% % 
% % %assemble the filename; make timestep_string optional
% % if exist('timestep_string','var')
% %     image_filename = [image_filename_start timestep_string image_filename_end];
% % else
% %     image_filename = [image_filename_start image_filename_end];
% % end
% % 
% % display(['Image file ' image_filename ' will be used.'])


%OPTIONAL: define the folder where the images to be processed are saved
%if 'folder' is not defined, the software will look for the images in the
%present working directory

folder = '/Users/skatrina/Documents/Data/Fonocal/20180215/sampleA_Capture10/';

if ~exist('folder','var')
    display('No data folder defined. The images will be read from the current directory.')
    folder = './';
else
    display(['Data will be read from the directory: ' folder])
end



%% User-defined settings for image processing and particle locating:
% These parameters are specific to be input into the Gao & Kilfoil particle
% locating software. If you are using the iterative locating algorithm with
% different base locating routines, the following will be different.


%the approximate size of the image noise, in pixels (x, y,and z directions):
lnoise = [2 2 2];

%the approximate particle diameter in pixels (x, y, and z directions):
lobject = false_particle_size; %KSM edited to make same as false_particle_size; from originally [x y z] format

%numbers somewhat larger than the largest particle diameter:
diameters = [10 10 25];
%numbers somewhat smaller than the smallest particle diameter:
mask_size = [4 4 70];

%a minimum separation cutoff to help prevent doubly-located particles 
%(set conservatively; this is to avoid a single particle being detected
%twice in any given particle locating iteration)
min_separation = [4 4 4];

% masscut parameters, determined by looking at the histogram of "masses"
% (i.e. the integrated intensity of the found particles); in the Gao &
% Kilfoil code feature3dMB.m, the masses have just been calculated by line 210, 
% so this is the appropriate place to pause and generate the histogram.
%
% Features having less integrated intensity than this cutoff will be
% excluded from future analysis. 
% 
% Note that during iterations, the histogram of masses has a long tail du
% to particles cut off by the edges of the sample; however, these will be
% excluded from later analysis as being too close to the edge, so there is
% no actual ambiguity between the low-mass noise and the high-mass "real"
% particles.

masscut_initial = 2.5e5; % first iteration generally only finds real particles to start with
masscut_residuals = 60000; % required nonzero to prevent residual noise from being detected as "particles" once all of the real particles have been found

%if the masscuts are not set, they will be set to default values by the locating
%program (useful if the user wants to run a first-pass of locating in order
%to plot the mass histogram)





