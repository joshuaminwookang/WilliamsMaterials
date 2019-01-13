% Cleaned up code to demonstrate point tracking with pkfnd/cntrd

% NAME:
%               raghu_radial_symmetry_example_code
% PURPOSE:
%               Demonstrates how to use radialcenter to find the
%               positions of a bunch of particles
% 
% CATEGORY:
%               Image Processing
%
% OUTPUTS:
%               p_sr:    A list of particle positions. First two columns are
%               particle positions. Third column is an estimate of the size
%               of the particle (see the description of 'sigma' in the help
%               file for radialcenter)
% NOTES:
%               Consists of three steps. 1) Run bpass to clean everything
%               up. 2) Use pkfnd to find the rough positions of the
%               particle. 3) Use radialcenter to find the centres of radial
%               symmetry near each of the rough particle positions.
%               - Needs pkfnd gui to make more user friendly
%               - There should be an example image1.tif in this folder that
%               this example works with.
%               - For a sequence of images taken in the same lighting
%               conditions, you should do this on a typical image. Then use
%               the parameters on every image in the sequence in a loop...

% MODIFICATION HISTORY:
%               Written by Rob Style Dec 2016
%
% POTENTIAL IMPROVEMENTS
%               At the moment, there will be an error if the particle
%               positions are too close to the edge of the image. When you
%               select a small area of the image (testim) to find the
%               radially-symmetric centre of, then you will have problems
%               if the window that you try to get goes outside the range of
%               indices of the original image. This should be easy to fix
%               by throwing away points that are too close to the edge.
%
%               Potentially find some way of throwing out blobs that aren't
%               too radially symmetric (e.g. where two particles point
%               spread functions are overlapping)?

%%

% Firstly work out the best parameters for the bpass
a=bpass_gui('image1.tif');

% I found a =[0 10  0] to give quite nice looking results.
% NB try and keep the second parameter in bpass_gui as small as possible
% without shrinking the blobs down too much.
%%
im_pass=bpass(imread('image1.tif'),a(1),a(2),a(3));

% Ideally here, we would put in a gui to optimise the pkfnd parameters too.
% For the moment, pkfnd has two parameters that you need to input. The
% first is a threshold - a minimum intensity that might be a peak. We'll
% just estimate this by taking the average intensity of the picture

th=mean(mean(im_pass));

% The second parameter is a size that is slightly larger than then diameter
% of the average blobs in the image. This is about the same as 

sz=9;

p=pkfnd(im_pass,th,sz); % Gives a list with all the particle positions

% Display the results of pkfnd
imagesc(imread('image1.tif'))
hold on
plot(p(:,1),p(:,2),'rx')

%% Now for each particle, make a little window around the particle centre
% and find the centre of radial symmetry of the particle.

% First we need an estimate of the size of the particle. This should be an
% odd number that is just bigger than the diameter of a particle in pixels, but not so big
% that this box overlaps with light from other particles.

im=double(imread('image1.tif'));
particle_sz=5;
clear p_rs
p_rs=zeros(length(p),3);

for i=1:length(p)
    xcentre=round(p(i,1));
    ycentre=round(p(i,2));
    testim=im(ycentre-(particle_sz-1)/2:ycentre+(particle_sz-1)/2,xcentre-(particle_sz-1)/2:xcentre+(particle_sz-1)/2);
    [c,d,f]=radialcenter(testim);
    p_rs(i,:)=[xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
end

figure
imagesc(im)
hold on
plot(p_rs(:,1),p_rs(:,2),'rx')