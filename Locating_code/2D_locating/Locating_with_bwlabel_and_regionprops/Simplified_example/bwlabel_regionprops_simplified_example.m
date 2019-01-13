% Cleaned up code to demonstrate point tracking with bwlabel/regionprops

% NAME:
%               bwlabel_regionprops_simplified_example
% PURPOSE:
%               Demonstrates how to use bwlabel and regionprops to find the
%               positions of a bunch of particles
%
% CATEGORY:
%               Image Processing
%
% OUTPUTS:
%               c:    A list of particle positions. First two columns are
%               particle positions. Third column is area of the features in
%               the image after erosion and dilation (if applied). This can
%               be used to filter the feature positions for specific size
%               ranges to exclude artifacts.
% NOTES:
%               - Needs  bpass_gui to make more user friendly
%               - There should be an example image1.tif in this folder that
%               this example works with.
%               - For a sequence of images taken in the same lighting
%               conditions, you should do this on a typical image. Then use
%               the parameters on every image in the sequence in a loop...

% MODIFICATION HISTORY:
%               Written by Hendrik Spanke Dec 2016
%
% POTENTIAL IMPROVEMENTS:
%               We could add something to detect how spherical the blobs
%               are, and throw them away if they're non-spherical (e.g. two
%               particles with point spread functions that are overlapping)

%%

% Firstly work out the best parameters for the bpass
a=bpass_gui('image1.tif');

% I found a =[0 10  0] to give quite nice looking results.
% NB try and keep the second parameter in bpass_gui as small as possible
% without shrinking the blobs down too much.
%
im_pass=bpass(imread('image1.tif'),a(1),a(2),a(3));

%% Prepare Image for Particles: ----------------------------------------
% bwlabel needs a binary image as input. Start GUI for binarization and
% work out best parameters.
b=binarize_gui(im_pass);

% Apply the selected parameters to the image
BW=binarize(im_pass,b(1),b(2),b(3));

%% Identify Particles and find Center/Area
% information on the area of the particle allows for the filtering
% of agglomerates or other features larger than the expected size
% of a normal particle.

im=double(imread('image1.tif'));
regions = bwlabel(BW,8);
%props = regionprops(regions, 'Centroid', 'Area'); % find features in BW (Particles)
props = regionprops(regions,im,'WeightedCentroid','Area'); % find features in BW (Particles)

X=[]; Y=[]; Area=[];
for n = 1:1:size(props)    % go through features in BW
    % Particles:
        X = [X; props(n).WeightedCentroid(1)];
        Y = [Y; props(n).WeightedCentroid(2)];
        Area = [Area; props(n).Area];   
end

c=cat(2,X,Y,Area);

%% Display the results

imagesc(imread('image1.tif'))
hold on
plot(c(:,1),c(:,2),'rx')
