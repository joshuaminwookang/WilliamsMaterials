% Example script to extract the 3d positions of fluorescent beads attached
% to the surface of PDMS

% First need to load the image file into a 3D array: (nx,ny,nz)
fname = 'fluorinert_droplet.tif';
info = imfinfo(fname);
num_images = numel(info);
for k = 1:num_images
    im(:,:,k) = imread(fname, k, 'Info', info);
end

%% Get input parameters
% The original code is rather hard to understand. Kate Jensen has gone
% through and worked out exactly what each of the parameters mean. This is
% a summary of what she has worked out.

% I recommend that you open the image up in imageJ so that you can work out
% some of the parameters. Particularly important is to work out the typical
% size of the particle blobs in x,y and z. The image also needs to be white
% blobs on a black background.

% lobject is a typical size of the particle in pixels
lobject=[7 7 11];

% lnoise is the approximate noise in the image. I suspect this should be
% [1 1 1] if you're not sure.
lnoise=[1 1 1];

% diameters is the dimensions of a blob that is somewhat bigger than the
% largest particle diameter. When there are big gaps between features, this
% can be something like:
diameters=round(1.75*lobject);
% If features are closer together, check out the example in
% Jensen_iterative_Locating_for_finding_every_particle/Iterative_locating_MATLAB_software

% mask_size is somewhat smaller than the smallest particle diameter. The
% code first finds the pixel-level positions of the centroid positions, and
% then uses a mask of this size around that position to find centroid of the particle to
% sub-pixel acuracy. It's important that this isn't too big so that it
% doesn't extend outside the particles.
mask_size = round(lobject/2);

%a minimum separation cutoff to help prevent doubly-located particles 
%(set conservatively; this is to avoid a single particle being detected
%twice in any given particle locating iteration)
min_separation = lobject;

% We also need to set a masscut parameter - but this is easier to do by
% having a look at the histogram outputted by the code. The masscut
% parameter gives the integrated intensity that must be inside a particle
% for it to count as a particle (to avoid detecting noise). Thus if a
% particle was a sphere of radius 7 pixels, and each pixel was at least 100
% in intensity, we would set the intensity to pi*7^2 *100.

% To work out the correct masscut, we first run the code to see the mass
% histogram.

% Firstly we run the bpass:

res=bpass3dMB(im, lnoise, lobject,[0 0]); % Have a look at the FFindinstruct.txt for further details about the other options for bpass3dMB

%% Next we run the code without a masscut
r=feature3dMB(res, diameters, mask_size, size(res), [1 0 0], min_separation, 0, 0); %the [1,0,0] here means we are using only the first optional setting: a minimum separation

% The code will output a histogram showing the integrated intensity of all
% the particles found. You can choose an appropriate value from this.

masscut=2.8*10^5;
%% Next we run the code with a masscut
r=feature3dMB(res, diameters, mask_size, size(res), [1 1 0], min_separation, masscut, 0);

%% Analysing the data
%  OUTPUTS:
% 		r(:,1):	this contains the x centroid positions, in pixels.
% 		r(:,2): this contains the y centroid positions, in pixels. 
% 		r(:,3): this contains the z centroid positions, in pixels.
% 		r(:,4): this contains the integrated brightness.
% 		r(:,5): this contains the squared radius of gyration. 
% 		r(:,6): this contains the peak height of the feature.
% 		r(:,7): this contains the fraction of voxels above threshold.
figure
plot3(r(:,1),r(:,2),r(:,3),'.')